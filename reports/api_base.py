from __future__ import unicode_literals

from functools import wraps
import logging
import re

from django.conf import settings
from django.conf.urls import url, patterns
from django.core.cache import cache
from django.core.exceptions import ObjectDoesNotExist
import django.core.exceptions
from django.core.signals import got_request_exception
from django.http.response import HttpResponseBase, HttpResponse, \
    HttpResponseNotFound, Http404
from django.utils import six
from django.utils.cache import patch_cache_control, patch_vary_headers
from django.views.decorators.csrf import csrf_exempt
from tastypie.authentication import Authentication
from tastypie.cache import NoCache
from tastypie.exceptions import ImmediateHttpResponse, BadRequest, NotFound
from tastypie.http import HttpBadRequest, HttpNotImplemented, HttpNoContent, \
    HttpApplicationError
from tastypie.resources import Resource, convert_post_to_put, sanitize
from tastypie.utils.urls import trailing_slash

from db.support.data_converter import default_converter
from reports import ValidationError, InformationError, API_RESULT_ERROR
from reports.serialize import XLSX_MIMETYPE, SDF_MIMETYPE, XLS_MIMETYPE, \
    CSV_MIMETYPE, JSON_MIMETYPE
from reports.serializers import BaseSerializer


logger = logging.getLogger(__name__)


def un_cache(_func):
    '''
    Wrapper function to disable caching for 
    SQLAlchemyResource.stream_response_from_statement and other caches
    ''' 
    @wraps(_func)
    def _inner(self, *args, **kwargs):
        logger.debug('decorator un_cache: %s, %s', self, _func )
        self.clear_cache()
        self.set_caching(False)
        result = _func(self, *args, **kwargs)
        self.set_caching(True)
        logger.debug('decorator un_cache done: %s, %s', self, _func )
        return result

    return _inner

class Authorization(object):

    def _is_resource_authorized(self, resource_name, user, permission_type):
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
    cache = NoCache()
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
        new_class = super(DeclarativeMetaclass, cls).__new__(cls, name, bases, attrs)
        opts = getattr(new_class, 'Meta', None)
        new_class._meta = ResourceOptions(opts)

        return new_class


class IccblBaseResource(six.with_metaclass(DeclarativeMetaclass)):
    """
    Override tastypie.resources.Resource:
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
                self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/schema%s$" % (
                self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<%s>.*?)%s$" % (
                self._meta.resource_name, self._meta.detail_uri_name, 
                trailing_slash()), self.wrap_view('dispatch_detail'), 
                name="api_dispatch_detail"),
        ]
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
        urls = self.prepend_urls()
        urls += self.base_urls()
        urlpatterns = patterns('',
            *urls
        )
        return urlpatterns
    
    def is_authenticated(self, request):
        """
        Handles checking if the user is authenticated and dealing with
        unauthenticated users.

        Mostly a hook, this uses class assigned to ``authentication`` from
        ``Resource._meta``.
        """
        # Authenticate the request as needed.
        auth_result = self._meta.authentication.is_authenticated(request)

        if isinstance(auth_result, HttpResponse):
            raise ImmediateHttpResponse(response=auth_result)

        if not auth_result is True:
            raise ImmediateHttpResponse(response=http.HttpUnauthorized())

    def error_response(self, request, errors, response_class=None):
        """
        Extracts the common "which-format/serialize/return-error-response"
        cycle.

        Should be used as much as possible to return errors.
        """
        if response_class is None:
            response_class = http.HttpBadRequest

        desired_format = None

        if request:
            if request.GET.get('callback', None) is None:
                try:
                    desired_format = self.determine_format(request)
                except BadRequest:
                    pass  # Fall through to default handler below
            else:
                # JSONP can cause extra breakage.
                desired_format = 'application/json'

        if not desired_format:
            desired_format = self._meta.default_format

        try:
            serialized = self.serialize(request, errors, desired_format)
        except BadRequest as e:
            error = "Additional errors occurred, but serialization of those errors failed."

            if settings.DEBUG:
                error += " %s" % e

            return response_class(content=error, content_type='text/plain')

        return response_class(
            content=serialized, content_type=build_content_type(desired_format))
    
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
        method = getattr(self, "%s_%s" 
            % (request_method, request_type), None)

        if method is None:
            raise ImmediateHttpResponse(response=HttpNotImplemented())

#         self.is_authenticated(request)
        convert_post_to_put(request)
        logger.info('calling method: %s.%s_%s', 
            self._meta.resource_name, request_method, request_type)
        response = method(request, **kwargs)

        # # If what comes back isn't a ``HttpResponse``, assume that the
        # # request was accepted and that some action occurred. This also
        # # prevents Django from freaking out.
        if not isinstance(response, HttpResponseBase):
            return HttpNoContent()
        
#         # Custom ICCB parameter: set cookie to tell the browser javascript
#         # UI that the download request is finished
#         downloadID = request.GET.get('downloadID', None)
#         if downloadID:
#             logger.info('set cookie "downloadID" %r', downloadID )
#             response.set_cookie('downloadID', downloadID)
#         else:
#             logger.debug('no downloadID: %s' % request.GET )

        return response

    def wrap_view(self, view):
        """
        Override the tastypie implementation to handle our own ValidationErrors.
        """

        @csrf_exempt
        def wrapper(request, *args, **kwargs):
            DEBUG_WRAPPER = True
            response = None
            try:
                callback = getattr(self, view)
                if DEBUG_WRAPPER:
                    msg = ()
                    if kwargs:
                        msg = [ (key,kwargs[key]) for key in kwargs.keys() 
                            if len(str(kwargs[key]))<100]
                    logger.info('callback: %r, %r', view, msg)
                    logger.info('request: %r', request)
                else:
                    logger.info('callback: %r, %r', callback, view)
                
                self.is_authenticated(request)
        
                response = callback(request, *args, **kwargs)
                # Our response can vary based on a number of factors, use
                # the cache class to determine what we should ``Vary`` on so
                # caches won't return the wrong (cached) version.
                varies = getattr(self._meta.cache, "varies", [])

                if varies:
                    patch_vary_headers(response, varies)

                if self._meta.cache.cacheable(request, response):
                    if self._meta.cache.cache_control():
                        # If the request is cacheable and we have a
                        # ``Cache-Control`` available then patch the header.
                        patch_cache_control(response, **self._meta.cache.cache_control())

                if request.is_ajax() and not response.has_header("Cache-Control"):
                    # IE excessively caches XMLHttpRequests, so we're disabling
                    # the browser cache here.
                    # See http://www.enhanceie.com/ie/bugs.asp for details.
                    patch_cache_control(response, no_cache=True)

            except BadRequest as e:
                # The message is the first/only arg
                logger.exception('Bad request exception: %r', e)
                data = {"error": sanitize(e.args[0]) if getattr(e, 'args') else ''}
                response = self.build_error_response(request, data, **kwargs)
            except InformationError as e:
                logger.exception('Information error: %r', e)
                response = self.build_error_response(
                    request, { 'Messages': e.errors }, **kwargs)
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
                    if settings.DEBUG and getattr(settings, 'TASTYPIE_FULL_DEBUG', False):
                        logger.warn('raise tastypie full exception for %r', e)
                        raise
    
                    # Rather than re-raising, we're going to things similar to
                    # what Django does. The difference is returning a serialized
                    # error message.
                    logger.exception('handle 500 error %r...', str(e))
                    response = self._handle_500(request, e)

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

    def clear_cache(self):
        logger.debug('clearing the cache from resource: %s (all caches cleared)' 
            % self._meta.resource_name)
        cache.clear()

    def set_caching(self,use_cache):
        self.use_cache = use_cache

    def deserialize(self, request, format=None):
        
        content_type = self._meta.serializer.get_content_type(request, format)
        logger.info('content_type: %r', content_type)
        if content_type.startswith('multipart'):
            logger.info('request.Files.keys: %r', request.FILES.keys())
            # process *only* one attached file
            if len(request.FILES.keys()) != 1:
                raise BadRequest({ 
                    'FILES': 'File upload supports only one file at a time',
                    'filenames': request.FILES.keys(),
                })
             
            if 'sdf' in request.FILES:  
                file = request.FILES['sdf']
                return self._meta.serializer.deserialize(
                    file.read(), SDF_MIMETYPE)
            elif 'xls' in request.FILES:
                file = request.FILES['xls']
                
                schema = self.build_schema()
                list_keys = [x for x,y in schema['fields'].items() 
                    if y.get('data_type') == 'list']
                return self._meta.serializer.deserialize(
                    file.read(), XLS_MIMETYPE, **{ 'list_keys': list_keys})
            elif 'csv' in request.FILES:
                file = request.FILES['csv']

                schema = self.build_schema()
                list_keys = [x for x,y in schema['fields'].items() 
                    if y.get('data_type') == 'list']
                
                return self._meta.serializer.deserialize(
                    file.read(), CSV_MIMETYPE, **{ 'list_keys': list_keys})
            else:
                raise BadRequest(
                    'Unsupported multipart file key: %r', request.FILES.keys())
        
        logger.info('use deserializer: %r', content_type)
        
        if content_type in [XLS_MIMETYPE,XLSX_MIMETYPE,CSV_MIMETYPE]:
            # NOTE: must inject information about the list fields, so that they
            # can be properly parsed (this is a custom serialization format)
            schema = self.build_schema()
            list_keys = [x for x,y in schema['fields'].items() 
                if y.get('data_type') == 'list']
            return self._meta.serializer.deserialize(
                request.body,content_type,
                **{ 'list_keys': list_keys})
        else:
            return self._meta.serializer.deserialize(request.body,content_type)
            
    def serialize(self, request, data, format=None):
        content_type = self._meta.serializer.get_accept_content_type(request, format)
        return self._meta.serializer.serialize(data, content_type)

    def _get_filename(self, readable_filter_hash, filename=None, **extra):
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
    
    def _get_filename_old(self, schema, kwargs,filename=''):
        logger.info('filename: %r, kwargs: %r', filename, kwargs )
        filekeys = [filename]
        if 'id_attribute' in schema:
            filekeys.extend([ str(kwargs[key]) for 
                key in schema['id_attribute'] if key in kwargs ])
        else:
            _dict = {key:val for key,val in kwargs.items() 
                if key not in [
                    'visibilities','exact_fields','api_name','resource_name',
                    'includes','order_by']}
            for i,(x,y) in enumerate(_dict.items()):
                filekeys.append(str(x))
                filekeys.append(str(y))
                if i == 10:
                    break
        filekeys = [x for x in filekeys if x.strip() ]
        filekeys.insert(0,self._meta.resource_name)
        filename = '_'.join(filekeys)
        filename = re.sub(r'[\W]+','_',filename)
        logger.debug('get_filename: %r, %r' % (filename, kwargs))
        return filename
    
    def get_content_type(self, request, format=None):
         
        return self._meta.serializer.get_content_type(request,format=format)
 
    def get_accept_content_type(self, request, format=None):
         
        return self._meta.serializer.get_accept_content_type(request,format=format)
        
    def build_response(self, request, data, response_class=HttpResponse, **kwargs):
        
        content_type = self._meta.serializer.get_content_type(
            request, format=kwargs.get('format', None))
        # FIXME: "data" must be a dict {'objects': [data] } for xls, csv serialization
        serialized = self.serialize(request, data, format=kwargs.get('format', None))
        response = response_class(
            content=serialized, 
            content_type=content_type)
        # FIXME: filename is not being set well here:
        # - used for downloads; reports.api resources use
        # this method to serialize; all others use streaming serializers.
        format = self._meta.serializer.get_format_for_content_type(content_type)
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
        return response 
    
    def build_error_response(self, request, data, response_class=HttpBadRequest, **kwargs):
        format = 'json'
        if kwargs and 'format' in kwargs:
            format = kwargs['format']
        try:
            return self.build_response(
                request, data, response_class=response_class, format=format)
        except Exception, e:
            logger.exception('On trying to serialize the error response: %r, %r',
                data, e)
            return HttpBadRequest(content=data, content_type='text/plain')

    def _handle_500(self, request, exception):
        ''' Override Tastypie for serialization'''
        import traceback
        import sys
        the_trace = '\n'.join(traceback.format_exception(*(sys.exc_info())))
        response_class = HttpApplicationError
        response_code = 500

        NOT_FOUND_EXCEPTIONS = (NotFound, ObjectDoesNotExist, Http404)

        if isinstance(exception, NOT_FOUND_EXCEPTIONS):
            response_class = HttpResponseNotFound
            response_code = 404

        if settings.DEBUG:
            data = {
                "error_message": sanitize(six.text_type(exception)),
                "traceback": the_trace,
            }
            return self.build_error_response(request, data, response_class=response_class)

        # When DEBUG is False, send an error message to the admins (unless it's
        # a 404, in which case we check the setting).
        send_broken_links = getattr(settings, 'SEND_BROKEN_LINK_EMAILS', False)

        if not response_code == 404 or send_broken_links:
            log = logging.getLogger('django.request.tastypie')
            log.error('Internal Server Error: %s' % request.path, exc_info=True,
                      extra={'status_code': response_code, 'request': request})

        # Send the signal so other apps are aware of the exception.
        got_request_exception.send(self.__class__, request=request)

        # Prep the data going out.
        data = {
            "error_message": getattr(
                settings, 
                'TASTYPIE_CANNED_ERROR', 
                "Sorry, this request could not be processed. Please try again later."),
        }
        return self.build_error_response(request, data, response_class=response_class)
        
       
