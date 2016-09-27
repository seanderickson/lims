from __future__ import unicode_literals

from functools import wraps
import logging
import re

from django.conf import settings
from django.core.cache import cache
import django.core.exceptions
from django.http.response import HttpResponseBase, HttpResponse,\
    HttpResponseNotFound, Http404
from django.utils.cache import patch_cache_control, patch_vary_headers
from django.views.decorators.csrf import csrf_exempt
from tastypie.exceptions import ImmediateHttpResponse, BadRequest, NotFound
from tastypie.http import HttpBadRequest, HttpNotImplemented, HttpNoContent,\
    HttpApplicationError
from tastypie.resources import Resource, convert_post_to_put, sanitize

from reports import HEADER_APILOG_COMMENT, ValidationError, InformationError, _now
from reports.models import ApiLog
from reports.serialize import XLSX_MIMETYPE, SDF_MIMETYPE, XLS_MIMETYPE, \
    CSV_MIMETYPE, JSON_MIMETYPE
from django.core.exceptions import ObjectDoesNotExist
from django.core.signals import got_request_exception
import six


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


class IccblBaseResource(Resource):
    """
    Override tastypie.resources.Resource:
    -- use StreamingHttpResponse or the HttpResponse
    -- control application specific caching
    """

#     content_types = {
#                      'xls': XLS_MIMETYPE,
#                      'xlsx': XLSX_MIMETYPE,
#                      'csv': CSV_MIMETYPE,
#                      'sdf': SDF_MIMETYPE,
#                      'json': JSON_MIMETYPE,
#                      }

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
                return self._meta.serializer.deserialize(
                    file.read(), CSV_MIMETYPE)
            else:
                raise BadRequest(
                    'Unsupported multipart file key: %r', request.FILES.keys())
        
        logger.info('use deserializer: %r', content_type)
        return self._meta.serializer.deserialize(request.body,content_type)


    def serialize(self, request, data, format=None):
        content_type = self._meta.serializer.get_accept_content_type(request, format)
        return self._meta.serializer.serialize(data, content_type)

    def _get_filename(self, schema, kwargs, filename=''):
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
        
#     def get_serialize_format(self, request, format=None):
#         
#         return self._meta.serializer.get_serialize_format(request, format=format)

    def build_response(self, request, data, response_class=HttpResponse, **kwargs):
        
        content_type = self._meta.serializer.get_content_type(
            request, format=kwargs.get('format', None))
        serialized = self.serialize(request, data, format=kwargs.get('format', None))
        return response_class(
            content=serialized, 
            content_type=content_type)

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
            "error_message": getattr(settings, 'TASTYPIE_CANNED_ERROR', "Sorry, this request could not be processed. Please try again later."),
        }
        return self.build_error_response(request, data, response_class=response_class)
        
    def wrap_view(self, view):
        """
        Override the tastypie implementation to handle our own ValidationErrors.
        """

        @csrf_exempt
        def wrapper(request, *args, **kwargs):
            DEBUG_WRAPPER = True
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

                return response
            except BadRequest as e:
                logger.exception('bad request...')
                # for BadRequest, the message is the first/only arg
                logger.exception('Bad request exception: %r', e)
                data = {"error": sanitize(e.args[0]) if getattr(e, 'args') else ''}
                return self.build_error_response(request, data, **kwargs)
            except InformationError as e:
                logger.exception('Information error: %r', e)
                response = self.build_error_response(
                    request, { 'Messages': e.errors }, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % 'errors'
                    downloadID = request.GET.get('downloadID', None)
                    if downloadID:
                        logger.info('set cookie "downloadID" %r', downloadID )
                        response.set_cookie('downloadID', downloadID)
                    else:
                        logger.debug('no downloadID: %s' % request.GET )
                return response
            
            except ValidationError as e:
                logger.exception('Validation error')
                response = self.build_error_response(
                    request, { 'errors': e.errors }, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % 'errors'
                    downloadID = request.GET.get('downloadID', None)
                    if downloadID:
                        logger.info('set cookie "downloadID" %r', downloadID )
                        response.set_cookie('downloadID', downloadID)
                    else:
                        logger.debug('no downloadID: %s' % request.GET )
                return response
            except django.core.exceptions.ValidationError as e:
                logger.exception('Django validation error: %s', e)
                response = self.build_error_response(
                    request, { 'errors': e.message_dict }, **kwargs)
                if 'xls' in response['Content-Type']:
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.xlsx' % 'errors'
                    downloadID = request.GET.get('downloadID', None)
                    if downloadID:
                        logger.info('set cookie "downloadID" %r', downloadID )
                        response.set_cookie('downloadID', downloadID)
                    else:
                        logger.debug('no downloadID: %s' % request.GET )
                return response
                
            except Exception as e:
                if hasattr(e, 'response'):
                    # A specific response was specified
                    return e.response

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
                return self._handle_500(request, e)

        return wrapper

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
        Override tastypie method to eliminate much unneeded functionality:
        - 'allowed' methods
        - HTTP_X_HTTP_METHOD_OVERRIDE
        - throttling
        - verification check of the response class 
        Other modifications:
        - use of the "downloadID" cookie
        """
        
        # allowed_methods = getattr(
        #     self._meta, "%s_allowed_methods" % request_type, None)
        # 
        # if 'HTTP_X_HTTP_METHOD_OVERRIDE' in request.META:
        #     request.method = request.META['HTTP_X_HTTP_METHOD_OVERRIDE']
        # 
        # request_method = self.method_check(request, allowed=allowed_methods)
        request_method = request.method.lower()
        method = getattr(self, "%s_%s" 
            % (request_method, request_type), None)

        if method is None:
            raise ImmediateHttpResponse(response=HttpNotImplemented())

        self.is_authenticated(request)
        # self.throttle_check(request)

        # All clear. Process the request.
        convert_post_to_put(request)
        logger.info('calling method: %s_%s, kwargs: %r',
            request_method, request_type,kwargs)
        response = method(request, **kwargs)

        # # Add the throttled request.
        # self.log_throttled_access(request)

        # # If what comes back isn't a ``HttpResponse``, assume that the
        # # request was accepted and that some action occurred. This also
        # # prevents Django from freaking out.
        if not isinstance(response, HttpResponseBase):
            return HttpNoContent()
        
        # Custom ICCB parameter: set cookie to tell the browser javascript
        # UI that the download request is finished
        downloadID = request.GET.get('downloadID', None)
        if downloadID:
            logger.info('set cookie "downloadID" %r', downloadID )
            response.set_cookie('downloadID', downloadID)
        else:
            logger.debug('no downloadID: %s' % request.GET )

        return response
       
