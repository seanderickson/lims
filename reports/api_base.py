# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import base64
from functools import wraps
import logging
import re

import django.core.urlresolvers
from django.conf import settings
from django.conf.urls import url, patterns, include
from django.contrib.auth import authenticate
from django.core.cache import cache
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied,\
    ImproperlyConfigured
import django.core.exceptions
from django.core.signals import got_request_exception
from django.http.response import HttpResponseBase, HttpResponse, \
    HttpResponseNotFound, Http404, HttpResponseForbidden, HttpResponseBadRequest, \
    HttpResponseServerError
from django.middleware.csrf import _sanitize_token, REASON_NO_REFERER, \
    REASON_BAD_REFERER, REASON_BAD_TOKEN
from django.utils import six
from django.utils.crypto import constant_time_compare
from django.utils.encoding import force_text
from django.utils.http import same_origin
from django.views.decorators.csrf import csrf_exempt

from reports.utils.django_requests import convert_request_method_to_put
from db.support.data_converter import default_converter
from reports import ValidationError, InformationError, BadRequestError, \
    BackgroundJobImmediateResponse, LoginFailedException, API_RESULT_ERROR
from reports.serialize import XLSX_MIMETYPE, SDF_MIMETYPE, XLS_MIMETYPE, \
    CSV_MIMETYPE, JSON_MIMETYPE
from reports.serializers import BaseSerializer, LimsSerializer


logger = logging.getLogger(__name__)

DEBUG = False
DEBUG_WRAPPER = True

TRAILING_SLASH = '/?'

def un_cache(_func):
    '''
    Wrapper function to disable caching for 
    SQLAlchemyResource.stream_response_from_statement and other caches
    ''' 
    @wraps(_func)
    def _inner(self, *args, **kwargs):
        logger.debug('decorator un_cache: %r, %r, %r', 
            self, _func,args )
        self.clear_cache(None, **kwargs)
        self.set_caching(False)
        result = _func(self, *args, **kwargs)
        self.set_caching(True)
        logger.debug('decorator un_cache done: %s, %s', self, _func )
        return result

    return _inner

# Authentication Classes taken from tastypie.authentication

class Authentication(object):
    """
    A simple base class to establish the protocol for auth.
    """
    def __init__(self):
        pass

    def is_authenticated(self, request, **kwargs):
        return True

    def get_identifier(self, request):
        """
        Provides a unique string identifier for the requestor.

        This implementation returns a combination of IP address and hostname.
        """
        return "%s_%s" % (request.META.get('REMOTE_ADDR', 'noaddr'), 
            request.META.get('REMOTE_HOST', 'nohost'))

#     def check_active(self, user):
#         """
#         Ensures the user has an active account.
# 
#         Optimized for the ``django.contrib.auth.models.User`` case.
#         """
#         return user.is_active


class MultiAuthentication(object):
    """
    An authentication backend that tries a number of backends in order.
    """
    def __init__(self, *backends, **kwargs):
        super(MultiAuthentication, self).__init__(**kwargs)
        self.backends = backends

    def is_authenticated(self, request, **kwargs):
        """
        Identifies if the user is authenticated to continue or not.
        """
        unauthorized = False

        for backend in self.backends:
            check = backend.is_authenticated(request, **kwargs)

            if check:
                if isinstance(check, HttpResponse):
                    unauthorized = unauthorized or check
                else:
                    request._authentication_backend = backend
                    return check

        return unauthorized

    def get_identifier(self, request):
        """
        Provides a unique string identifier for the requestor.

        This implementation returns a combination of IP address and hostname.
        """
        try:
            return request._authentication_backend.get_identifier(request)
        except AttributeError:
            return 'nouser'

      
class IccblBasicAuthentication(Authentication):
    """
    Handles HTTP Basic auth against a specific auth backend if provided,
    or against all configured authentication backends using the
    ``authenticate`` method from ``django.contrib.auth``.

    Optional keyword arguments:

    ``backend``
        If specified, use a specific ``django.contrib.auth`` backend instead
        of checking all backends specified in the ``AUTHENTICATION_BACKENDS``
        setting.
    ``realm``
        The realm to use in the ``HttpUnauthorized`` response.  Default:
        ``screensaver``.
    NOTE: modified from tastypie.authentication.BasicAuthentication:
    - modified to support the format <superuser_username>:<username>:<password>
    format used by the LIMS to support login "as user" by the superuser.
    """
    def __init__(self, backend=None, realm='screensaver', **kwargs):
        super(IccblBasicAuthentication, self).__init__(**kwargs)
        self.backend = backend
        self.realm = realm

    def _unauthorized(self):
        response = HttpResponse(status=401)
        # FIXME: Sanitize realm.
        response['WWW-Authenticate'] = 'Basic Realm="%s"' % self.realm
        return response

    def is_authenticated(self, request, **kwargs):
        """
        Checks a user's basic auth credentials against the current
        Django auth backend.

        Should return either ``True`` if allowed, ``False`` if not or an
        ``HttpResponse`` if you need something custom.
        """
        if not request.META.get('HTTP_AUTHORIZATION'):
            return False
            # return self._unauthorized()

        try:
            (auth_type, data) = request.META['HTTP_AUTHORIZATION'].split()
            if auth_type.lower() != 'basic':
                return self._unauthorized()
            user_pass = base64.b64decode(data).decode('utf-8')
        except Exception,e:
            logger.exception('auth fail: %r', e)
            return self._unauthorized()
        
        bits = user_pass.split(':')
        if len(bits) > 3:
            logger.error('user-pass split returns > 3 strings: %d', len(bits))
            return self._unauthorized()
        elif len(bits) == 3:
            # Special case, supports proxy login by superuser "as" user using format:
            # <superuser_username>:<username>:<password>
            # recombine the <superuser_username>:<username> and allow the 
            # reports.auth backend to process using the 
            # reports.auth.USER_PROXY_LOGIN_PATTERN 
            username = '%s:%s' % (bits[0],bits[1])
            password = bits[2]
        elif len(bits) == 2:
            username = bits[0]
            password = bits[1]
        else:
            logger.error('user-pass split returns only one element')
            return self._unauthorized()

        if self.backend:
            user = self.backend.authenticate(username=username, password=password)
        else:
            user = authenticate(username=username, password=password)

        if user is None:
            logger.info('no user found for: %s', username)
            return self._unauthorized()

        if user.is_active is not True:
            logger.info('user is not active: %s', username)
            return False

        request.user = user
        return True

    def get_identifier(self, request):
        """
        Provides a unique string identifier for the requestor.

        This implementation returns the user's basic auth username.
        """
        return request.META.get('REMOTE_USER', 'nouser')

class IccblSessionAuthentication(Authentication):
    """
    Replaces Tastypie authentication.SessionAuthentication:
    - update to support Django >1.4 "csrfmiddlewaretoken"
    @see django.middleware.csrf
    
    From tastypie comments:
    
    Cargo-culted from Django 1.3/1.4's ``django/middleware/csrf.py``.
    An authentication mechanism that piggy-backs on Django sessions.

    This is useful when the API is talking to Javascript on the same site.
    Relies on the user being logged in through the standard Django login
    setup.

    Requires a valid CSRF token.
    """
    def is_authenticated(self, request, **kwargs):
        """
        Checks to make sure the user is logged in & has a Django session.
        """
        try:
            csrf_token = _sanitize_token(
                    request.COOKIES[settings.CSRF_COOKIE_NAME])
        except KeyError:
            logger.error('Check that basic auth is working; '
                'requires user.is_active flag be set')
            logger.error('cookies: %r', request.COOKIES)
            logger.error('reject: NO CSRF cookie: %r', settings.CSRF_COOKIE_NAME)
            return False

        if request.method in ('GET', 'HEAD', 'OPTIONS', 'TRACE'):
            return request.user.is_authenticated()

        if getattr(request, '_dont_enforce_csrf_checks', False):
            return request.user.is_authenticated()
        
        if request.is_secure():
            # Suppose user visits http://example.com/
            # An active network attacker (man-in-the-middle, MITM) sends a
            # POST form that targets https://example.com/detonate-bomb/ and
            # submits it via JavaScript.
            #
            # The attacker will need to provide a CSRF cookie and token, but
            # that's no problem for a MITM and the session-independent
            # nonce we're using. So the MITM can circumvent the CSRF
            # protection. This is true for any HTTP connection, but anyone
            # using HTTPS expects better! For this reason, for
            # https://example.com/ we need additional protection that treats
            # http://example.com/ as completely untrusted. Under HTTPS,
            # Barth et al. found that the Referer header is missing for
            # same-domain requests in only about 0.2% of cases or less, so
            # we can use strict Referer checking.
            referer = force_text(
                request.META.get('HTTP_REFERER'),
                strings_only=True,
                errors='replace'
            )
            if referer is None:
                logger.error('reject: %s', REASON_NO_REFERER)
                return False

            # Note that request.get_host() includes the port.
            good_referer = 'https://%s/' % request.get_host()
            if not same_origin(referer, good_referer):
                reason = REASON_BAD_REFERER % (referer, good_referer)
                logger.error('reject: %s', reason)
                return False

        if csrf_token is None:
            # No CSRF cookie. For POST requests, we insist on a CSRF cookie,
            # and in this way we can avoid all CSRF attacks, including login
            # CSRF.
            logger.error('reject: NO CSRF cookie')
            return False

        # Check non-cookie token for match.
        request_csrf_token = ""
        if request.method == "POST":
            try:
                request_csrf_token = request.POST.get('csrfmiddlewaretoken', '')
            except IOError:
                # Handle a broken connection before we've completed reading
                # the POST data. process_view shouldn't raise any
                # exceptions, so we'll ignore and serve the user a 403
                # (assuming they're still listening, which they probably
                # aren't because of the error).
                pass

        if request_csrf_token == "":
            # Fall back to X-CSRFToken, to make things easier for AJAX,
            # and possible for PUT/DELETE.
            request_csrf_token = request.META.get('HTTP_X_CSRFTOKEN', '')

        if not constant_time_compare(request_csrf_token, csrf_token):
            logger.error('%r != %r', request_csrf_token, csrf_token)
            logger.error('reject: %s', REASON_BAD_TOKEN)
            return False
        
        return request.user.is_authenticated()

    def get_identifier(self, request):
        """
        Provides a unique string identifier for the requestor.

        This implementation returns the user's username.
        """

        return getattr(request.user, get_username_field())       


class Authorization(object):

    def _is_resource_authorized(
        self, user, permission_type, resource_name=None, **kwargs):

        if resource_name is None:
            resource_name = self.resource_name
          
        raise NotImplementedError(
            '_is_resource_authorized must be implemented for %s, %s, %s',
            resource_name, user, permission_type)


class ResourceOptions(object):
    """
    A configuration class for ``Resource``.

    Provides sane defaults and the logic needed to augment these settings with
    the internal ``class Meta`` used on ``Resource`` subclasses.
    """
    serializer = BaseSerializer()
    authentication = Authentication()
    authorization = Authorization()
    api_name = None
    resource_name = None
    alt_resource_name = None
    object_class = None
    queryset = None
    always_return_data = False
    collection_name = 'objects'
    detail_uri_name = 'pk'

    def __new__(cls, meta=None):
        overrides = {}

        # Handle overrides.
        if meta:
            for override_name in dir(meta):
                # No internals please.
                if not override_name.startswith('_'):
                    overrides[override_name] = getattr(meta, override_name)

        if six.PY3:
            return object.__new__(type('ResourceOptions', (cls,), overrides))
        else:
            return object.__new__(type(b'ResourceOptions', (cls,), overrides))

class DeclarativeMetaclass(type):

    def __new__(cls, name, bases, attrs):
        new_class = \
            super(DeclarativeMetaclass, cls).__new__(cls, name, bases, attrs)
        opts = getattr(new_class, 'Meta', None)
        new_class._meta = ResourceOptions(opts)

        return new_class

class IccblBaseResource(six.with_metaclass(DeclarativeMetaclass)):
    """
    see tastypie.resources.Resource:
    -- use StreamingHttpResponse or the HttpResponse
    -- control application specific caching
    """

    def __init__(self, api_name=None):

        if not api_name is None:
            self._meta.api_name = api_name

    def base_urls(self):
        """
        The standard URLs this ``Resource`` should respond to.
        """
        return [
            url(r"^(?P<resource_name>%s)%s$" % (
                self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/schema%s$" % (
                self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name="api_get_schema"),
        ]

    def dispatch_clear_cache(self, request, **kwargs):
        self.clear_cache(request, **kwargs)
        return self.build_response(request, 'ok', **kwargs)

    def prepend_urls(self):
        """
        A hook for adding your own URLs or matching before the default URLs.
        """
        return []

    @property
    def urls(self):
        """
        The endpoints this ``Resource`` responds to.

        Mostly a standard URLconf, this is suitable for either automatic use
        when registered with an ``Api`` class or for including directly in
        a URLconf should you choose to.
        """
        urls = [            
            url(r"^(?P<resource_name>%s)/clear_cache%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_clear_cache'), name="api_clear_cache"),
        ]
        urls += self.prepend_urls()
        urls += self.base_urls()
        urlpatterns = patterns('',
            *urls
        )
        return urlpatterns
    
    def is_authenticated(self, request):
        """
        Note: authentication will raise a PermissionDenied error if if fails
        """
        
        auth_result = self._meta.authentication.is_authenticated(request)
        if auth_result is not True:
            raise PermissionDenied
        return auth_result

    def dispatch_list(self, request, **kwargs):
        """
        A view for handling the various HTTP methods (GET/POST/PUT/DELETE) over
        the entire list of resources.

        Relies on ``Resource.dispatch`` for the heavy-lifting.
        """
        return self.dispatch('list', request, **kwargs)

    def dispatch_detail(self, request, **kwargs):
        """
        A view for handling the various HTTP methods (GET/POST/PUT/DELETE) on
        a single resource.

        Relies on ``Resource.dispatch`` for the heavy-lifting.
        """
        return self.dispatch('detail', request, **kwargs)

    def dispatch(self, request_type, request, **kwargs):
        """
        From original Tastypie structure:
        - lookup the API method
        # TODO: does very little at this point; refactor into wrap_view
        """
        
        request_method = request.method.lower()
        method_name = "%s_%s" % (request_method, request_type)
        method = getattr(self, method_name, None)

        if method is None:
            raise Http404('"%s" is not implemented for "%s"' 
                % (method_name, self._meta.resource_name))

        convert_request_method_to_put(request)
        response = method(request, **kwargs)
        
        # FIXME: remove this, require all types to return a response
        if not isinstance(response, HttpResponseBase):
            return HttpResponseNotFound()
        
        return response

    def exception_handler(self,_func):
        '''
        Exception handler wrapper for IccblBaseResource classes
        '''
        
        @wraps(_func)
        def _inner(self, *args, **kwargs):
            logger.debug('self: %r, func: %r, args: %r, kwargs: %r', 
                self, _func, args, kwargs)
            
            request = args[0]
            try:
                response = _func(*args, **kwargs)
            
            except BackgroundJobImmediateResponse as e:
                logger.info('BackgroundJobImmediateResponse returned: %r', 
                    e.httpresponse)
                response = e.httpresponse
            except BadRequestError as e:
                # The message is the first/only arg
                logger.exception('Bad request exception: %r', e)
                response = self.build_error_response(
                    request, e.errors, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % API_RESULT_ERROR
            except InformationError as e:
                logger.exception('Information error: %r', e)
                response = self.build_error_response(
                    request, e.errors, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % API_RESULT_ERROR
            except ValidationError as e:
                logger.exception('Validation error: %r', e)
                response = self.build_error_response(
                    request, { API_RESULT_ERROR: e.errors }, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % API_RESULT_ERROR
            except django.core.exceptions.ValidationError as e:
                logger.exception('Django validation error: %s', e)
                response = self.build_error_response(
                    request, { API_RESULT_ERROR: e.message_dict }, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % API_RESULT_ERROR
            except LoginFailedException as e:
                logger.info('LoginFailedException ex: %r', e)
                data = {
                    'Login Failed: ': '%s' % e }
                response = self.build_error_response(
                    request, data, response_class=HttpResponseForbidden, **kwargs)
            except PermissionDenied as e:
                logger.info('PermissionDenied ex: %r', e)
                data = {
                    'Permission Denied: ': '%s' % e }
                response = self.build_error_response(
                    request, data, response_class=HttpResponseForbidden, **kwargs)
            except ObjectDoesNotExist as e:
                logger.info('not found: %r', e)
                response = self.build_error_response(
                    request, { 'msg': '%r' % e }, 
                    response_class=HttpResponseNotFound, **kwargs)
            except Http404 as e:
                logger.info('not found: %r', e)
                response = self.build_error_response(
                    request, { 'msg': '%r' % e }, 
                    response_class=HttpResponseNotFound, **kwargs)
            except Exception as e:
                logger.exception('Unhandled exception: %r', e)
                if hasattr(e, 'response'):
                    # A specific response was specified
                    response = e.response
                else:
                    logger.exception('Unhandled exception: %r', e)
    
                    # A real, non-expected exception.
                    # Handle the case where the full traceback is more helpful
                    # than the serialized error.
                    if settings.DEBUG:
                        
                        logger.warn('raise full exception for %r', e)
                        raise
    
                    # Rather than re-raising, we're going to things similar to
                    # what Django does. The difference is returning a serialized
                    # error message.
                    logger.exception('handle 500 error %r...', str(e))
                    response = self._handle_500(request, e)
            
            return response
        return _inner

    
    
    def wrap_view(self, view):
        """
        Handle common operations for views:
        - authentication
        - exception handling
        - client side cache headers
        - download cookie
        
        """

        # NOTE: csrf is enforced by SessionAuthentication, all resources
        @csrf_exempt
        def wrapper(request, *args, **kwargs):
            response = None
            callback = getattr(self, view)
            if DEBUG_WRAPPER:
                msg = ()
                if kwargs:
                    msg = [ (key,str(kwargs[key])[:100]) 
                        for key in kwargs.keys() ]
                logger.info(
                    'wrap_view: %r, method: %r, user: %r, request: %r, '
                    'kwargs: %r', 
                    view, request.method, request.user, request, msg)
            else:
                logger.info('wrap_view: %r, %r', view, request)
            
            # re-wrap to include authentication
            def _inner(*args, **kwargs):
                self.is_authenticated(request)
                return callback(*args,**kwargs)
            
            response = self.exception_handler(_inner)(self,request, *args, **kwargs)

            # Custom ICCB parameter: set cookie to tell the browser javascript
            # UI that the download request is finished
            downloadID = request.GET.get('downloadID', None)
            if downloadID:
                logger.info('set cookie "downloadID" %r', downloadID )
                response.set_cookie('downloadID', downloadID)
            else:
                logger.debug('no downloadID: %s' % request.GET )
             
            return response
        return wrapper
    
    def set_caching(self,use_cache):
        logger.debug('set_caching: %r, %r', use_cache, self._meta.resource_name)
        self.use_cache = use_cache

    def get_cache(self):
        raise NotImplementedError('get_cache must be implemented: %r', self._meta.resource_name)

    def deserialize(self, request, format=None, schema=None):
        
        # TODO: refactor to use delegate deserialize to the serializer class:
        # - allow Resource classes to compose functionality vs. inheritance
        
        if format is not None:
            content_type = self.get_serializer() \
                               .get_content_type_for_format(format)
        else:
            content_type = self.get_serializer().get_content_type(request)
        logger.info('deserialize for content_type: %r', content_type)
        
        if schema is None:
            schema = self.build_schema(user=request.user)
        # NOTE: Injecting information about the list fields, so that they
        # can be properly parsed (this is a custom serialization format)
        list_keys = [x for x,y in schema['fields'].items() 
            if y.get('data_type') == 'list']
            
        
        if content_type.startswith('multipart'):
            if not request.FILES:
                logger.error('No "FILES" found with multipart content.')
                raise BadRequestError({
                    'ContentType':
                    'No "FILES" found with multipart content: %r' % content_type
                })
                
            logger.info('request.Files.keys: %r', request.FILES.keys())
            
            # NOTE: may process *only* one attached file per upload
            if len(request.FILES.keys()) != 1:
                raise BadRequestError({ 
                    'FILES': 'File upload supports only one file at a time',
                    'filenames': request.FILES.keys(),
                })
             
            # FIXME: rework to use the multipart Content-Type here
            if 'sdf' in request.FILES:  
                file = request.FILES['sdf']
                return (
                    self.get_serializer().deserialize(
                        file.read(), SDF_MIMETYPE), 
                    { 'filename': file.name })
            elif 'xls' in request.FILES:
                file = request.FILES['xls']
                return (
                    self.get_serializer().deserialize(
                        file.read(), XLS_MIMETYPE, list_keys=list_keys), 
                    { 'filename': file.name } )
            elif 'csv' in request.FILES:
                file = request.FILES['csv']
                return (
                    self.get_serializer().deserialize(
                        file.read(), CSV_MIMETYPE, list_keys=list_keys), 
                    { 'filename': file.name } )
            else:
                raise BadRequestError({
                    'Files':
                    'Unsupported multipart file keys: %r' % request.FILES.keys()})
        
        elif content_type in [XLS_MIMETYPE,XLSX_MIMETYPE,CSV_MIMETYPE]:
            return (
                self.get_serializer().deserialize(
                    request.body,content_type, list_keys=list_keys), None )
        else:
            return self.get_serializer().deserialize(request.body,content_type), None
            
    def serialize(self, data, content_type):
        logger.debug('serialize to: %r', content_type)
        return self.get_serializer().serialize(data, content_type)

    def _get_filename(self, readable_filter_hash, schema, filename=None, **extra):
        MAX_VAL_LENGTH = 20
        file_elements = [self._meta.resource_name]
        if filename is not None:
            file_elements.append(filename)
        if extra is not None:
            for key,val in extra.items():
                file_elements.append(str(key))
                if val is not None:
                    val = default_converter(str(val))
                    val = val[:MAX_VAL_LENGTH]
                    file_elements.append(val)
        for key,val in readable_filter_hash.items():
            if key not in schema['id_attribute']:
                file_elements.append(str(key))
            val = default_converter(str(val))
            val = val[:MAX_VAL_LENGTH]
            file_elements.append(val)
                
        if len(file_elements) > 1:
            # Add an extra separator for the resource name
            file_elements.insert(1,'_')
            
        filename = '_'.join(file_elements)
        logger.info('filename: %r', filename)
        MAX_FILENAME_LENGTH = 128
        filename = filename[:128]
        return filename
    
    def get_serializer(self):
        return self._meta.serializer
        
    def build_response(self, request, data, response_class=HttpResponse, 
                       format=None, **kwargs):
        if format is not None:
            content_type = \
                self.get_serializer().get_content_type_for_format(format)
        else:
            content_type = \
                self.get_serializer().get_accept_content_type(request)
        logger.debug(
            'build response for data: %r, content type: %r', data, content_type)
        serialized = self.serialize(data, content_type)
        response = response_class(
            content=serialized, 
            content_type=content_type)
        
        # FIXME: filename is not being set well here:
        # - used for downloads; reports.api resources use
        # this method to serialize; all others use streaming serializers.

        format = self.get_serializer().get_format_for_content_type(content_type)
        if format != 'json':
            filename = kwargs.get('filename', None)
            if filename is None:
                filename = self._meta.resource_name
            if format == 'csv': 
                filename += '.csv'
            elif format == 'xls': 
                filename += '.xls'
            response['Content-Disposition'] = \
                'attachment; filename=%s' % filename
        if 'status_code' in kwargs:
            response.status_code = kwargs['status_code']
            
        logger.debug('response: %r: %r', response, response.status_code)
        return response 
    
    def build_error_response(
            self, request, data, response_class=HttpResponseBadRequest, **kwargs):

        try:
            return self.build_response(
                request, data, response_class=response_class, **kwargs)
        except Exception, e:
            logger.exception('On trying to serialize the error response: %r, %r',
                data, e)
            return HttpResponseBadRequest(content=data, content_type='text/plain')

    def _print_exception_message(text):
        # from tastypie.resources.sanitize
        return escape(six.text_type(exception))\
            .replace('&#39;', "'").replace('&quot;', '"')

    def _handle_500(self, request, exception):
        ''' 
        '''
        import traceback
        import sys
        the_trace = '\n'.join(traceback.format_exception(*(sys.exc_info())))
        response_class = HttpResponseServerError
        response_code = 500

        if settings.DEBUG:
            data = {
                "error_message": _print_exception_message(exception),
                "traceback": the_trace,
            }
            return self.build_error_response(
                request, data, response_class=response_class)
        
        # TODO: configure logging email on errors

        # Send the signal so other apps are aware of the exception.
        got_request_exception.send(self.__class__, request=request)

        data = {
            "SERVER_ERROR": 
                "Sorry, this request could not be processed. Please try again later.",
        }
        return self.build_error_response(
            request, data, response_class=response_class)
        


class Api(object):

    def __init__(self, api_name="v1"):
        self.api_name = api_name
        self._registry = {}

    def register(self, resource):
        resource_name = getattr(resource._meta, 'resource_name', None)

        if resource_name is None:
            raise ImproperlyConfigured(
                'Resource %r must define a "resource_name".' % resource)

        self._registry[resource_name] = resource


    def unregister(self, resource_name):
        if resource_name in self._registry:
            del(self._registry[resource_name])

    def wrap_view(self, view):
        def wrapper(request, *args, **kwargs):
            try:
                return getattr(self, view)(request, *args, **kwargs)
            except BadRequestError:
                return HttpResponseBadRequest()
        return wrapper

    def prepend_urls(self):
        """
        A hook for adding your own URLs or matching before the default URLs.
        """
        return []

    @property
    def urls(self):
        """
        Provides URLconf details for the ``Api`` and all registered
        ``Resources`` beneath it.
        """
        pattern_list = [
            url(
                r"^(?P<api_name>%s)%s$" % (self.api_name, TRAILING_SLASH), 
                self.wrap_view('top_level'), 
                name="api_%s_top_level" % self.api_name),
        ]

        for name in sorted(self._registry.keys()):
            self._registry[name].api_name = self.api_name
            pattern_list.append(
                url(r"^(?P<api_name>%s)/" % self.api_name, 
                    include(self._registry[name].urls)))

        urlpatterns = self.prepend_urls()


        urlpatterns += pattern_list
        return urlpatterns

    def top_level(self, request, api_name=None):
        """
        A view that returns a serialized list of all resources registers
        to the ``Api``. Useful for discovery.
        """
        fullschema = parse_val(
            request.GET.get('fullschema', False),
            'fullschema', 'boolean')

        available_resources = {}

        if api_name is None:
            api_name = self.api_name

        for name, resource in self._registry.items():
            if not fullschema:
                schema = self._build_reverse_url("api_get_schema", kwargs={
                    'api_name': api_name,
                    'resource_name': name,
                })
            else:
                schema = resource.build_schema()

            available_resources[name] = {
                'list_endpoint': 
                    self._build_reverse_url("api_dispatch_list", kwargs={
                        'api_name': api_name,
                        'resource_name': name,
                    }),
                'schema': schema,
            }

        serializer = LimsSerializer()
        content_type = serializer.get_accept_content_type(request)
        serialized =  serializer.serialize(available_resources, content_type)
        return HttpResponse(
            content=serialized, 
            content_type=content_type)

    def _build_reverse_url(self, name, args=None, kwargs=None):
        return django.core.urlresolvers.reverse(name, args=args, kwargs=kwargs)
