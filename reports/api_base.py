# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import base64
from functools import wraps
import logging
import re

from django.conf import settings
from django.conf.urls import url, include
import django.contrib.auth
from django.core.cache import cache
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied, \
    ImproperlyConfigured
import django.core.exceptions
from django.core.signals import got_request_exception
from django.http.response import HttpResponseBase, HttpResponse, \
    HttpResponseNotFound, Http404, HttpResponseForbidden, HttpResponseBadRequest, \
    HttpResponseServerError
from django.middleware.csrf import _sanitize_token, REASON_NO_REFERER, \
    REASON_BAD_REFERER, REASON_BAD_TOKEN, REASON_MALFORMED_REFERER, \
    REASON_INSECURE_REFERER
import django.urls
from django.utils import six
from django.utils.crypto import constant_time_compare
from django.utils.encoding import force_text
from django.views.decorators.csrf import csrf_exempt

from reports.utils import default_converter
from reports import ValidationError, InformationError, BadRequestError, \
    ApiNotImplemented, BackgroundJobImmediateResponse, LoginFailedException, \
    API_RESULT_ERROR
from reports.auth import DEBUG_AUTHENTICATION
from reports.serialize import XLSX_MIMETYPE, SDF_MIMETYPE, XLS_MIMETYPE, \
    CSV_MIMETYPE, JSON_MIMETYPE
from reports.serializers import BaseSerializer, LimsSerializer
from reports.utils.django_requests import convert_request_method_to_put
from django.utils.http import is_same_domain


# from django.utils.http import same_origin
# NOTE: API design is based loosely on the tastypie project: 
# see attributions, e.g.:
# From tastypie.compat
# Compatability for salted vs unsalted CSRF tokens;
# Django 1.10's _sanitize_token also hashes it, so it can't be compared directly.
# Solution is to call _sanitize_token on both tokens, then unsalt or noop both
try:
    from django.middleware.csrf import _unsalt_cipher_token
    def unsalt_token(token):
        return _unsalt_cipher_token(token)
except ImportError:
    def unsalt_token(token):
        return token


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

class Authentication(object):
    def __init__(self):
        pass

    def is_authenticated(self, request, **kwargs):
        '''
        @return True if successful authentication is performed
        '''
        return True

class MultiAuthentication(object):
    '''
    Authenticate using the given authentication_clients, in order.
    '''
    def __init__(self, *authentication_clients, **kwargs):
        '''
        @param authentication_clients in order, to try
        '''
        super(MultiAuthentication, self).__init__(**kwargs)
        self.authentication_clients = authentication_clients

    def is_authenticated(self, request, **kwargs):

        for authentication_client in self.authentication_clients:
            
            if DEBUG_AUTHENTICATION:
                logger.info('try authentication: %r', authentication_client)
                
            auth_result = authentication_client.is_authenticated(request, **kwargs)
            
            if DEBUG_AUTHENTICATION:
                logger.info('authentication_client: %r, result: %r',
                    authentication_client, auth_result)
            
            if auth_result is True:
                return True
            
        return False
    
class IccblBasicAuthentication(Authentication):
    '''
    Handles HTTP Basic auth using django.contrib.auth.authenticate
    
    NOTE: modified from tastypie.authentication.BasicAuthentication:
    - modified to support the format <superuser_username>:<username>:<password>
    format used by the LIMS to support login 'as user' by the superuser.
    '''
    def __init__(self, **kwargs):
        super(IccblBasicAuthentication, self).__init__(**kwargs)

    def is_authenticated(self, request, **kwargs):
        '''
        Checks a user's basic auth credentials using
        django.contrib.auth.authenticate
        
        @return True if successful authentication is performed
        '''
        
        header_authorization = request.META.get('HTTP_AUTHORIZATION')
        if not header_authorization:
            if DEBUG_AUTHENTICATION:
                logger.info('HTTP_AUTHORIZATION header not found')
            return False

        try:
            (auth_type, data) = header_authorization.split()
            if auth_type.lower() != 'basic':
                if DEBUG_AUTHENTICATION:
                    logger.info('HTTP_AUTHORIZATION type is not "basic"')
                return False
            user_pass = base64.b64decode(data).decode('utf-8')
        except Exception,e:
            logger.exception('auth fail on decoding user/password: %r', e)
            raise LoginFailedException({
                'HTTP_AUTHORIZATION': 'failed to decode basic auth header'})
        
        bits = user_pass.split(':')
        if len(bits) > 3:
            logger.error('user-pass split returns > 3 strings: %d', len(bits))
            raise LoginFailedException({
                'HTTP_AUTHORIZATION': 'failed to decode basic auth header: '
                    'user:password format'})
        elif len(bits) == 3:
            # Special case, supports proxy login by superuser 'as' user using format:
            # <superuser_username>:<username>:<password>
            # recombine the <superuser_username>:<username> and allow the 
            # django.contrib.auth.authenticate backend to process using the 
            # reports.auth.USER_PROXY_LOGIN_PATTERN 
            username = '%s:%s' % (bits[0],bits[1])
            password = bits[2]
        elif len(bits) == 2:
            username = bits[0]
            password = bits[1]
        else:
            logger.error('user-pass split returns only one element')
            raise LoginFailedException({
                'HTTP_AUTHORIZATION': 'failed to decode basic auth header: '
                    'user:password format'})

        if DEBUG_AUTHENTICATION:
            logger.info(
                'calling django.contrib.auth.authenticate for user %r', username)
        user = django.contrib.auth.authenticate(
            username=username, password=password)

        if user is None:
            logger.info('no user found for: "%s"', username)
            raise LoginFailedException({
                'HTTP_AUTHORIZATION': 'unknown user: %s' % username })

        if user.is_active is not True:
            logger.info('user is not active: %s', username)
            raise LoginFailedException({
                'HTTP_AUTHORIZATION': 'user is not active: "%s"' % username })

        request.user = user
        return True

class IccblSessionAuthentication(Authentication):
    '''
    Use the Django session to validate that the current user is authenticated.
    
    Note: see tastypie.authentication.SessionAuthentication for original 
    implementation - updated to support Django >1.4 "csrfmiddlewaretoken".
    
    Note: Requires a valid CSRF token, "csrfmiddlewaretoken" may be passed as a 
    POST parameter
    
    @see django.middleware.csrf v1.8 (used as a guide)
    
    '''
    def is_authenticated(self, request, **kwargs):
        '''
        Checks to make sure the user is logged in & has a Django session.
        '''
        if request.method in ('GET', 'HEAD', 'OPTIONS', 'TRACE'):
            return bool(request.user.is_authenticated)

        if getattr(request, '_dont_enforce_csrf_checks', False):
            logger.info('_dont_enforce_csrf_checks is set...')
            return bool(request.user.is_authenticated)

        csrf_token = _sanitize_token(request.COOKIES.get(settings.CSRF_COOKIE_NAME))
        if csrf_token is None:
            # No CSRF cookie. For POST requests, required.
            logger.error('reject: NO CSRF cookie')
            return False
        elif DEBUG_AUTHENTICATION:
            logger.info('Found cookie: %r: %r',
                settings.CSRF_COOKIE_NAME, csrf_token)

        if request.is_secure():
            if DEBUG_AUTHENTICATION:
                logger.info('perform secure session check.')
            
            if self._django_csrf_check(request) is not True:
                return False
#             referer = request.META.get('HTTP_REFERER')
# 
#             if referer is None:
#                 return False
# 
#             good_referer = 'https://%s/' % request.get_host()
# 
#             if not same_origin(referer, good_referer):
#                 return False

        request_csrf_token = ''
        if request.method == 'POST':
            # Look for POSTED csrf token:
            # Use the >1.4 Django Forms token key: "csrfmiddlewaretoken"
            try:
                request_csrf_token = request.POST.get('csrfmiddlewaretoken', '')
                if DEBUG_AUTHENTICATION:
                    logger.info(
                        'SessionAuthentication: POST csrf token (%r): %r', 
                        'csrfmiddlewaretoken', request_csrf_token)
            except IOError:
                # Handle a broken connection before we've completed reading
                # the POST data. 
                # (assuming they're still listening, which they probably
                # aren't because of the error).
                pass

        if request_csrf_token == '':
            # Fall back to X-CSRFToken, used by Ajax clients
            request_csrf_token = request.META.get('HTTP_X_CSRFTOKEN', '')
            if DEBUG_AUTHENTICATION:
                logger.info(
                    'SessionAuthentication: POST csrf token (%r): %r', 
                    'HTTP_X_CSRFTOKEN', request_csrf_token)

        request_csrf_token = _sanitize_token(request_csrf_token)

        if not constant_time_compare(unsalt_token(request_csrf_token),
                                     unsalt_token(csrf_token)):
            logger.warn('CSRF tokens do not match: %r', request)
            return False

        return bool(request.user.is_authenticated)

    def _django_csrf_check(self, request):
        '''
        Taken from django.middleware.csrf:
        @see django.middleware.csrf
        '''
        
        referer = force_text(
            request.META.get('HTTP_REFERER'),
            strings_only=True,
            errors='replace'
        )
        if referer is None:
            logger.error('csrf reject: %s', REASON_NO_REFERER)
            return False
        
        from django.utils.six.moves.urllib.parse import urlparse

        referer = urlparse(referer)

        # Make sure we have a valid URL for Referer.
        if '' in (referer.scheme, referer.netloc):
            logger.error('csrf reject: %s', REASON_MALFORMED_REFERER)
            return False

        # Ensure that our Referer is also secure.
        if referer.scheme != 'https':
            logger.error('csrf reject: %s', REASON_INSECURE_REFERER)
            return False

        # If there isn't a CSRF_COOKIE_DOMAIN, assume we need an exact
        # match on host:port. If not, obey the cookie rules.
        if settings.CSRF_COOKIE_DOMAIN is None:
            # request.get_host() includes the port.
            good_referer = request.get_host()
        else:
            good_referer = settings.CSRF_COOKIE_DOMAIN
            server_port = request.get_port()
            if server_port not in ('443', '80'):
                good_referer = '%s:%s' % (good_referer, server_port)

        # Here we generate a list of all acceptable HTTP referers,
        # including the current host since that has been validated
        # upstream.
        good_hosts = list(settings.CSRF_TRUSTED_ORIGINS)
        good_hosts.append(good_referer)

        if not any(is_same_domain(referer.netloc, host) for host in good_hosts):
            reason = REASON_BAD_REFERER % referer.geturl()
            logger.info('csrf reject: %r', reason)
            return False
        
        return True
        
class Authorization(object):

    def _is_resource_authorized(
        self, user, permission_type, resource_name=None, **kwargs):

        if resource_name is None:
            resource_name = self.resource_name
          
        raise NotImplementedError(
            '_is_resource_authorized must be implemented for %s, %s, %s',
            resource_name, user, permission_type)


class ResourceOptions(object):
    '''
    A configuration class for ``Resource``.

    Provides sane defaults and the logic needed to augment these settings with
    the internal ``class Meta`` used on ``Resource`` subclasses.
    '''
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
    '''
    see tastypie.resources.Resource:
    -- use StreamingHttpResponse or the HttpResponse
    -- control application specific caching
    '''

    def __init__(self, api_name=None):

        if not api_name is None:
            self._meta.api_name = api_name

    def base_urls(self):
        '''
        The standard URLs this ``Resource`` should respond to.
        '''
        return [
            url(r'^(?P<resource_name>%s)%s$' % (
                self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_list'), name='api_dispatch_list'),
            url(r'^(?P<resource_name>%s)/schema%s$' % (
                self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name='api_get_schema'),
        ]

    def dispatch_clear_cache(self, request, **kwargs):
        self.clear_cache(request, **kwargs)
        return self.build_response(request, 'ok', **kwargs)

    def dispatch_clear_all_caches(self, request, **kwargs):
        self.clear_cache(request, all=True)
        return self.build_response(request, 'ok', **kwargs)

    def prepend_urls(self):
        '''
        A hook for adding your own URLs or matching before the default URLs.
        '''
        return []

    @property
    def urls(self):
        '''
        The endpoints this ``Resource`` responds to.

        Mostly a standard URLconf, this is suitable for either automatic use
        when registered with an ``Api`` class or for including directly in
        a URLconf should you choose to.
        '''
        urls = [            
            url(r'^(?P<resource_name>%s)/clear_cache%s$' 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_clear_cache'), name='api_clear_cache'),
            url(r'^(?P<resource_name>%s)/clear_all_caches%s$' 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_clear_all_caches'), name='api_clear_cache'),
        ]
        urls += self.prepend_urls()
        urls += self.base_urls()
        return urls
    
    def is_authenticated(self, request):
        '''
        Returns True or raises an Exception for authentication failures
        '''
        
        auth_result = self._meta.authentication.is_authenticated(request)
        if auth_result is not True:
            logger.info('auth_result: %r', auth_result)
            raise LoginFailedException('auth result: %r' % auth_result)
        return auth_result

    def dispatch_list(self, request, **kwargs):
        '''
        A view for handling the various HTTP methods (GET/POST/PUT/DELETE) over
        the entire list of resources.

        Relies on ``Resource.dispatch`` for the heavy-lifting.
        '''
        return self.dispatch('list', request, **kwargs)

    def dispatch_detail(self, request, **kwargs):
        '''
        A view for handling the various HTTP methods (GET/POST/PUT/DELETE) on
        a single resource.

        Relies on ``Resource.dispatch`` for the heavy-lifting.
        '''
        return self.dispatch('detail', request, **kwargs)

    def dispatch(self, request_type, request, **kwargs):
        '''
        From original Tastypie structure:
        - lookup the API method
        # TODO: does very little at this point; refactor into wrap_view
        '''
        
        request_method = request.method.lower()
        method_name = '%s_%s' % (request_method, request_type)
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
            except ValidationError as e:
                logger.exception('Validation error: %r', e.errors)
                response = self.build_error_response(
                    request, { API_RESULT_ERROR: e.errors }, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % API_RESULT_ERROR
            except LoginFailedException as e:
                logger.info('LoginFailedException ex: %r', e)
                if hasattr(e, 'error_dict'):
                    data = { API_RESULT_ERROR: e.message_dict }
                else:
                    data = { 'Login failed': str(e) }
                response = self.build_error_response(
                    request, data, status_code=401, **kwargs)
                if request.META.get('HTTP_AUTHORIZATION'):
                    response['WWW-Authenticate'] = \
                        'Basic realm="%s"' % settings.BASIC_AUTH_REALM
            except django.core.exceptions.ValidationError as e:
                logger.exception('Django validation error: %s', e)
                if hasattr(e, 'error_dict'):
                    data = { API_RESULT_ERROR: e.message_dict }
                else:
                    data = { API_RESULT_ERROR: str(e) }
                response = self.build_error_response(
                    request, data, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % API_RESULT_ERROR
            except PermissionDenied as e:
                logger.exception('PermissionDenied ex: %r', e)
                data = {
                    'Permission Denied: ': '%s' % e }
                response = self.build_error_response(
                    request, data, response_class=HttpResponseForbidden, **kwargs)
            except ObjectDoesNotExist as e:
                logger.exception('not found: %r', e)
                response = self.build_error_response(
                    request, { 'msg': '%r' % e }, 
                    response_class=HttpResponseNotFound, **kwargs)
            except Http404 as e:
                logger.exception('not found: %r', e)
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
        '''
        Handle common operations for views:
        - authentication
        - exception handling
        - client side cache headers
        - download cookie
        
        '''

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
                username = str(request.user)
                if hasattr(request.user, 'username'):
                    username = getattr(request.user,'username')
                logger.info(
                    'wrap_view: %r, method: %r, username: %r, request: %r, '
                    'kwargs: %r', 
                    view, request.method, username, request, msg)
                logger.info('REMOTE_ADDR: %r', request.META.get('REMOTE_ADDR'))
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
                logger.debug('set cookie "downloadID" %r', downloadID )
                response.set_cookie('downloadID', downloadID)
            else:
                logger.debug('no downloadID: %s' % request.GET )
             
            return response
        return wrapper
    
    def set_caching(self,use_cache):
        logger.debug('set_caching: %r, %r', use_cache, self._meta.resource_name)
        self.use_cache = use_cache

    def get_cache(self):
        raise ApiNotImplemented(self._meta.resource_name, 'get_cache')

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
            elif 'xlsx' in request.FILES:
                file = request.FILES['xlsx']
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
                'error_message': _print_exception_message(exception),
                'traceback': the_trace,
            }
            return self.build_error_response(
                request, data, response_class=response_class)
        
        # TODO: configure logging email on errors

        # Send the signal so other apps are aware of the exception.
        got_request_exception.send(self.__class__, request=request)

        data = {
            'SERVER_ERROR': 
                'Sorry, this request could not be processed. Please try again later.',
        }
        return self.build_error_response(
            request, data, response_class=response_class)
        


class Api(object):

    def __init__(self, api_name='v1'):
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
        '''
        A hook for adding your own URLs or matching before the default URLs.
        '''
        return []

    @property
    def urls(self):
        '''
        Provides URLconf details for the ``Api`` and all registered
        ``Resources`` beneath it.
        '''
        pattern_list = [
            url(
                r'^(?P<api_name>%s)%s$' % (self.api_name, TRAILING_SLASH), 
                self.wrap_view('top_level'), 
                name='api_%s_top_level' % self.api_name),
        ]

        for name in sorted(self._registry.keys()):
            self._registry[name].api_name = self.api_name
            pattern_list.append(
                url(r'^(?P<api_name>%s)/' % self.api_name, 
                    include(self._registry[name].urls)))

        urlpatterns = self.prepend_urls()


        urlpatterns += pattern_list
        return urlpatterns

    def top_level(self, request, api_name=None):
        '''
        A view that returns a serialized list of all resources registers
        to the ``Api``. Useful for discovery.
        '''
        fullschema = parse_val(
            request.GET.get('fullschema', False),
            'fullschema', 'boolean')

        available_resources = {}

        if api_name is None:
            api_name = self.api_name

        for name, resource in self._registry.items():
            if not fullschema:
                schema = self._build_reverse_url('api_get_schema', kwargs={
                    'api_name': api_name,
                    'resource_name': name,
                })
            else:
                schema = resource.build_schema()

            available_resources[name] = {
                'list_endpoint': 
                    self._build_reverse_url('api_dispatch_list', kwargs={
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
        return django.urls.reverse(name, args=args, kwargs=kwargs)
