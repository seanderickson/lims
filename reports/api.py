import datetime
import json
import logging
import sys
import os
import re
import traceback
from collections import defaultdict, OrderedDict
from copy import deepcopy
from django.conf.urls import url
from django.conf import settings
from django.core.cache import cache
from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.models import User
from django.db.models.aggregates import Max
from django.db.models import Q
from django.db import transaction
from django.utils.encoding import smart_text
from django.utils import timezone, importlib
from django.forms.models import model_to_dict
from django.http.response import HttpResponseBase
from tastypie.exceptions import NotFound, ImmediateHttpResponse, Unauthorized,\
    BadRequest
from tastypie.bundle import Bundle
from tastypie.authorization import Authorization, ReadOnlyAuthorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication,\
    MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
# NOTE: tastypie.fields is required for dynamic field instances using eval
from tastypie import fields 
import tastypie.resources
from tastypie.resources import Resource, ModelResource
from tastypie.utils.urls import trailing_slash

from reports.serializers import LimsSerializer, CsvBooleanField, CSVSerializer
from reports.models import MetaHash, Vocabularies, ApiLog, ListLog, Permission, \
                           UserGroup, UserProfile, Record
# import lims.settings 
from tastypie.utils.timezone import make_naive
from tastypie.http import HttpForbidden
from tastypie.validation import Validation
from reports.dump_obj import dumpObj
from tastypie.utils.dict import dict_strip_unicode_keys
from functools import wraps
import six
from django.http.response import HttpResponse
import csv
from tastypie.utils.mime import build_content_type

logger = logging.getLogger(__name__)


class UserGroupAuthorization(Authorization):
    
    def _is_resource_authorized(self, resource_name, user, permission_type):
        logger.debug(str(("_is_resource_authorized", resource_name, user, permission_type)))
        scope = 'resource'
        prefix = 'permission'
        uri_separator = '/'
        permission_str =  uri_separator.join([prefix,scope,resource_name,permission_type])       

        logger.debug(str(( 'authorization query:', permission_str, 
            'user', user, user.is_superuser)))
        
        if user.is_superuser:
            return True
        
        userprofile = user.userprofile
        
        permissions = [x for x in 
            userprofile.permissions.all().filter(
                scope=scope, key=resource_name, type=permission_type)]
        
        if permissions:
            if(logger.isEnabledFor(logging.DEBUG)):
                logger.debug(str(('user',user ,'auth query', permission_str, 
                    'found matching user permissions', permissions)))
            return True
        
        logger.debug(str(('user',user ,'auth query', permission_str, 
            'not found in user permissions', permissions)))
        
        permissions_group = [ permission 
                for group in userprofile.usergroup_set.all() 
                for permission in group.get_all_permissions(
                    scope=scope, key=resource_name, type=permission_type)]
        if permissions_group:
            if(logger.isEnabledFor(logging.DEBUG)):
                logger.debug(str(('user',user ,'auth query', permission_str,
                    'found matching usergroup permissions', permissions_group)))
            return True
        
        logger.info(str(('user',user ,'auth query', permission_str,
             'not found in group permissions', permissions_group)))
        
        # Note: the TP framework raises the "Unauthorized" error: it then 
        # translates this into the (incorrect) HttpUnauthorized (401) response
        # Instead, raise an immediate exception with the correct 403 error code
        raise ImmediateHttpResponse(response=HttpForbidden(
            str(('user',user ,'permission not found', permission_str))))
    
    def read_list(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'read'):
            return object_list

    def read_detail(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'read'):
            return True

    def create_list(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'write'):
            return object_list

    def create_detail(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'write'):
            return True

    def update_list(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'write'):
            return object_list

    def update_detail(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'write'):
            return True

    def delete_list(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'write'):
            return object_list

    def delete_detail(self, object_list, bundle):
        if self._is_resource_authorized(
            self.resource_meta.resource_name, bundle.request.user, 'write'):
            return True


class SuperUserAuthorization(ReadOnlyAuthorization):
    
    def delete_list(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return object_list
        raise Unauthorized("Only superuser may delete lists.")

    def delete_detail(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return object_list
        raise Unauthorized("Only superuser may delete.")
 
    def create_list(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return object_list
        raise Unauthorized("Only superuser may create lists.")

    def create_detail(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return True
        raise Unauthorized("Only superuser may create.")

    def update_list(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return object_list
        raise Unauthorized("Only superuser may update lists.")

    def update_detail(self, object_list, bundle):
        logger.info(str(('create detail', bundle.request.user)))
        if bundle.request.user.is_superuser:
            return True
        raise Unauthorized("Only superuser may update.")


# TODO: this class should be constructed as a Mixin, not inheritor of ModelResource
class PostgresSortingResource(ModelResource):

    def __init__(self, **kwargs):
        super(PostgresSortingResource,self).__init__( **kwargs)

    def apply_sorting(self, obj_list, options):
        '''
        Override sorting so that we can make postgres sort nulls as less than
         everything else.

        Caveat: this will not work with joined fields unless they have an alias.  
        This is because it creates a field like:
        (screensaver_user_id is null) AS "screensaver_user_id_null"
        - if this field is duplicated in two sides of a join, then it must be 
        referenced by an alias, or as "table".screensaver_user_id, 
        and we are not supporting table specifications in this method, so if 
        joined fields are used, they must be referenced by alias only.
        '''
        
        obj_list = super(PostgresSortingResource, self).apply_sorting(
            obj_list, options)
        logger.debug(str(('order_by', obj_list.query.order_by)))
        extra_select = {}
        non_null_fields = options.get('non_null_fields', [])
        new_ordering = []

        logger.info(str(('==== non null fields', non_null_fields))) 
        for field in obj_list.query.order_by:
            original_field = field
            is_null_dir = '-'  # default nulls first for ascending
            if field.startswith('-'):
                is_null_dir = ''
                field = field[1:]
            if field in non_null_fields:
                continue
            extra_select[field+"_null"]=field + ' is null'
            new_ordering.append(is_null_dir + field+"_null")
            new_ordering.append(original_field)
        logger.debug(str(('extra_select', extra_select, 
                          'new_ordering', new_ordering)))
        obj_list = obj_list.extra(extra_select)

        obj_list.query.clear_ordering(force_empty=True)
        obj_list.query.add_ordering(*new_ordering)
        
        return obj_list


def log_obj_create(obj_create_func):
    @transaction.atomic()
    @wraps(obj_create_func)
    def _inner(self, bundle, **kwargs):
        logger.debug(str(('decorator start; log_obj_create', kwargs)))
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('----log obj_create', bundle)))

        kwargs = kwargs or {}    
        
        # Note: "full" logs would log all of the created data in the resource;
        # whereas "not full" only logs the creation; with this strategy, logs 
        # must be played backwards to recreate an entity state.
        full = False
        
        bundle = obj_create_func(self, bundle=bundle, **kwargs)
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('object created', bundle.obj )))
        log = ApiLog()
        log.username = bundle.request.user.username 
        log.user_id = bundle.request.user.id 
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
        if full:
            log.diffs = json.dumps(bundle.obj)
            
        # user can specify any valid, escaped json for this field
        # FIXME: untested
        if 'apilog_json_field' in bundle.data:
            log.json_field = json.dumps(bundle.data['apilog_json_field'])
            
        log.uri = self.get_resource_uri(bundle)
        log.key = '/'.join([str(x) for x in self.detail_uri_kwargs(bundle).values()])

        # TODO: how do we feel about passing form data in the headers?
        # TODO: abstract the form field name
        if 'HTTP_APILOG_COMMENT' in bundle.request.META:
            log.comment = bundle.request.META['HTTP_APILOG_COMMENT']
            
        if 'parent_log' in kwargs:
            log.parent_log = kwargs.get('parent_log', None)
        
        log.save()
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('create, api log', log)) )

        logger.debug(str(('decorator done; log_obj_create')))
        return bundle    
                            
    return _inner

# TODO: make class decorator as follows
# http://stackoverflow.com/questions/666216/decorator-classes-in-python
# class NullDecl (object):
#     def __init__ (self, func):
#         self.func = func
#         for n in set(dir(func)) - set(dir(self)):
#             setattr(self, n, getattr(func, n))
# 
#     def __call__ (self, * args):
#         return self.func (*args)
#     def __repr__(self):
#         return self.func

def is_empty_diff(difflog):
    if not difflog:
     return True
    
    empty = True;
    for key, value in difflog.items():
        if value:
            empty = False;
    return empty

def compare_dicts(dict1, dict2, excludes=['resource_uri'], full=False):
    '''
    @param full (default False) - a full compare shows added keys as well as diff keys
    
    Note: "full" logs would log all of the created data in the resource;
    whereas "not full" only logs the creation; with this strategy, logs 
    must be played backwards to recreate an entity state.
    '''
    original_keys = set(dict1.keys())-set(excludes)
    updated_keys = set(dict2.keys())-set(excludes)
    
    intersect_keys = original_keys.intersection(updated_keys)
    log = {'diffs': {}}
    
    added_keys = list(updated_keys - intersect_keys)
    if len(added_keys)>0: 
        log['added_keys'] = added_keys
        if full:
            log['diffs'].update( 
                dict(zip( added_keys,([None,dict2[key]] for key in added_keys if dict2[key]) )) )
    
    removed_keys = list(original_keys- intersect_keys)
    if len(removed_keys)>0: 
        log['removed_keys'] = removed_keys
        if full:
            log['diffs'].update(
                dict(zip(removed_keys,([dict1[key],None] for key in removed_keys if dict1[key]) )) )
    
    diff_keys = list()
    for key in intersect_keys:
        val1 = dict1[key]
        val2 = dict2[key]
        # NOTE: Tastypie converts to tz naive on serialization; then it 
        # forces it to the default tz upon deserialization (in the the 
        # DateTimeField convert method); for the purpose of this comparison,
        # then, make both naive.
        if isinstance(val2, datetime.datetime):
            val2 = make_naive(val2)
        if val1 != val2: 
            diff_keys.append(key)
    # Note, simple equality not used, since the serialization isn't 
    # symmetric, e.g. see datetimes, where tz naive dates look like UTC 
    # upon serialization to the ISO 8601 format.
    #         diff_keys = \
    #             list( key for key in intersect_keys 
    #                     if dict1[key] != dict2[key])

    if len(diff_keys)>0: 
        log['diff_keys'] = diff_keys
        log['diffs'].update(
            dict(zip(diff_keys,([dict1[key],dict2[key]] for key in diff_keys ) )))
#             dict(zip(diff_keys,([dict1[key],dict2[key]] for key in diff_keys if (full or dict1[key]) ) )))
    
    return log


def log_obj_update(obj_update_func):
    
    @transaction.atomic()
    @wraps(obj_update_func)
    def _inner(self, bundle, **kwargs):
        logger.debug(str(('decorator start; log_obj_update', 
            self,self._meta.resource_name)))
        kwargs = kwargs or {}    
        
        original_bundle = Bundle(
            data={},
            request=bundle.request)
        if hasattr(bundle,'obj'): 
            original_bundle.obj = bundle.obj
        else:
            bundle = self._locate_obj(bundle, **kwargs)
            original_bundle.obj = bundle.obj
        
        # store and compare dehydrated outputs: 
        # the api logger is concerned with what's sent out of the system, i.e.
        # the dehydrated output, not the internal representations.
        
        ## filter out the fields that aren't actually on this resource (well fields on reagent)
        schema = self.build_schema()
        fields = schema['fields']
        original_bundle.data = { key: original_bundle.data[key] 
            for key in original_bundle.data.keys() if key in fields.keys() }
        original_bundle = self.full_dehydrate(original_bundle)
        updated_bundle = obj_update_func(self,bundle, **kwargs)
        updated_bundle = self.full_dehydrate(updated_bundle)
        difflog = compare_dicts(original_bundle.data, updated_bundle.data)

        log = ApiLog()
        log.username = bundle.request.user.username 
        log.user_id = bundle.request.user.id 
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
        log.uri = self.get_resource_uri(bundle)
        log.key = '/'.join([str(x) for x in self.detail_uri_kwargs(bundle).values()])
        
        log.diff_dict_to_api_log(difflog)

        # user can specify any valid, escaped json for this field
        if 'apilog_json_field' in bundle.data:
            log.json_field = bundle.data['apilog_json_field']
        
        # TODO: abstract the form field name
        if 'HTTP_APILOG_COMMENT' in bundle.request.META:
            log.comment = bundle.request.META['HTTP_APILOG_COMMENT']
            logger.info(str(('log comment', log.comment)))

        if 'parent_log' in kwargs:
            log.parent_log = kwargs.get('parent_log', None)
        
        log.save()
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('update, api log', log)) )
        
        return updated_bundle
    
    return _inner
    
def log_obj_delete(obj_delete_func):    
    
    @transaction.atomic()
    @wraps(obj_delete_func)
    def _inner(self, bundle, **kwargs):
        logger.debug('---log obj_delete')
        kwargs = kwargs or {}    
        
        obj_delete_func(self,bundle=bundle, **kwargs)
        
        log = ApiLog()
        log.username = bundle.request.user.username 
        log.user_id = bundle.request.user.id 
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
                    
        # user can specify any valid, escaped json for this field
        if 'apilog_json_field' in bundle.data:
            log.json_field = json.dumps(bundle.data['apilog_json_field'])
        log.uri = self.get_resource_uri(bundle)
        log.key = '/'.join([str(x) for x in self.detail_uri_kwargs(bundle).values()])

        # TODO: how do we feel about passing form data in the headers?
        # TODO: abstract the form field name
        if 'HTTP_APILOG_COMMENT' in bundle.request.META:
            log.comment = bundle.request.META['HTTP_APILOG_COMMENT']

        if 'parent_log' in kwargs:
            log.parent_log = kwargs.get('parent_log', None)
            
        log.save()
        logger.debug(str(('delete, api log', log)) )
        
        return bundle


def log_patch_list(patch_list_func):    
      
    @transaction.atomic()
    def _inner(self, request, **kwargs):
        logger.info(str(('create an apilog for the patch list')))
        listlog = ApiLog()
        listlog.username = request.user.username 
        listlog.user_id = request.user.id 
        listlog.date_time = timezone.now()
        listlog.ref_resource_name = self._meta.resource_name
        listlog.api_action = 'PATCH_LIST'
        listlog.uri = self.get_resource_uri()
        # TODO: how do we feel about passing form data in the headers?
        # TODO: abstract the form field name
        if 'HTTP_APILOG_COMMENT' in request.META:
            listlog.comment = request.META['HTTP_APILOG_COMMENT']
         
        listlog.save();
        listlog.key = listlog.id
        listlog.save()
        
        kwargs = kwargs or {}    
        if not kwargs.get('parent_log', None):
            kwargs['parent_log'] = listlog
                  
        response = patch_list_func(self, request, **kwargs) 
        
#         self.listlog = None
         
        return response        
    return _inner


class LoggingMixin(Resource):
    '''
    Intercepts obj_create, obj_update and creates an ApiLog entry for the action
    
    Note: whatever is being extended with the LoggingMixin must also define a
    "detail_uri_kwargs" method that returns an _ordered_dict_, since we log the 
    kwargs as ordered args.
    ** note: "detail_uri_kwargs" returns the set of lookup keys for the resource 
    URI construction.
    '''

    def make_log(self, request, **kwargs):
        log = ApiLog()
        log.username = request.user.username 
        log.user_id = request.user.id 
        log.date_time = timezone.now()
        log.api_action = str((request.method)).upper()

        # TODO: how do we feel about passing form data in the headers?
        # TODO: abstract the form field name
        if 'HTTP_APILOG_COMMENT' in request.META:
            log.comment = request.META['HTTP_APILOG_COMMENT']
    
        if kwargs:
            for key, value in kwargs.items():
                if hasattr(log, key):
                    setattr(log, key, value)
        
        return log
    
    def _locate_obj(self,bundle, **kwargs):
        # lookup the object, the same way that it would be looked up in 
        # ModelResource.obj_update
        if not bundle.obj or not self.get_bundle_detail_data(bundle):
            try:
                lookup_kwargs = self.lookup_kwargs_with_identifiers(bundle, kwargs)
            except:
                # if there is trouble hydrating the data, fall back to just
                # using kwargs by itself (usually it only contains a "pk" key
                # and this will work fine.
                lookup_kwargs = kwargs

            try:
                bundle.obj = self.obj_get(bundle=bundle, **lookup_kwargs)
            except ObjectDoesNotExist:
                raise NotFound(("A model instance matching the provided "
                                " arguments could not be found: ", lookup_kwargs))
        return bundle
    
    @log_obj_create
    def obj_create(self, bundle, **kwargs): 
        return super(LoggingMixin, self).obj_create(bundle, **kwargs)
    
    @log_obj_update
    def obj_update(self, bundle, **kwargs):
        return super(LoggingMixin, self).obj_update(bundle, **kwargs)
      
    @log_obj_delete
    def obj_delete(self, bundle, **kwargs):
        return super(LoggingMixin, self).obj_delete(bundle, **kwargs) 
    
    @log_patch_list
    def patch_list(self, request, **kwargs):
        return Resource.patch_list(self,request, **kwargs) 
         
#     @transaction.atomic()
#     def patch_list(self, request, **kwargs):
#         ''' Override
#         '''
#         # create an apilog for the patch list
# #         listlog = self.listlog = ApiLog()
#         listlog = ApiLog()
#         listlog.username = request.user.username 
#         listlog.user_id = request.user.id 
#         listlog.date_time = timezone.now()
#         listlog.ref_resource_name = self._meta.resource_name
#         listlog.api_action = 'PATCH_LIST'
#         listlog.uri = self.get_resource_uri()
#         # TODO: how do we feel about passing form data in the headers?
#         # TODO: abstract the form field name
#         if 'HTTP_APILOG_COMMENT' in request.META:
#             listlog.comment = request.META['HTTP_APILOG_COMMENT']
#             
#         response =  Resource.patch_list(self, request, **kwargs) 
#          
#         returned_objects = self._meta.serializer.deserialize(response.content, response['Content-Type'])
#  
#         uris = [ obj['resource_uri'] for obj in returned_objects['objects']]
#         listlog.save();
#  
#         logger.debug(str(('adding created obj to list', uris)))
# #         # TODO: this should be a linked to many field?
# #         listlog.added_keys = json.dumps(uris)
#         for uri in uris:
#              
#             ## Convert uri to ref_resource_name, key
#             ref_resource_name = listlog.ref_resource_name
#             key = None
#             if ref_resource_name in uri:
#                 index = uri.find(ref_resource_name)+len(ref_resource_name)+1
#                 if len(uri) >= index:
#                     key = uri[index:]
#              
#             if not key:
#                 logger.warn(str((
#                     "returned patch_list uris don't contain the resource name/key",
#                     ref_resource_name, uri)))
#             sub = ListLog(uri=uri, ref_resource_name=ref_resource_name, key=key, apilog=listlog)
#             sub.save()
#          
#         listlog.key = listlog.id
#         listlog.save()
# #         self.listlog = None
#     
#         return response        

class StreamingResource(Resource):
    """
    Override tastypie.resources.Resource to replace check:
     if not isinstance(response, HttpResponse):
        return http.HttpNoContent()
    with:
     if not isinstance(response, HttpResponseBase):
        return http.HttpNoContent()
    -- this allows for use of the StreamingHttpResponse or the HttpResponse
    """

    content_types = {'json': 'application/json',
                     'xls': 'application/xls',
                     'csv': 'text/csv',
                     'sdf': 'chemical/x-mdl-sdfile',
                     }

    def dispatch(self, request_type, request, **kwargs):
        """
        Override tastypie.resources.Resource to replace check:
         if not isinstance(response, HttpResponse):
            return http.HttpNoContent()
        with:
         if not isinstance(response, HttpResponseBase):
            return http.HttpNoContent()
        -- this allows for use of the StreamingHttpResponse or the HttpResponse
        """
        allowed_methods = getattr(self._meta, "%s_allowed_methods" % request_type, None)

        if 'HTTP_X_HTTP_METHOD_OVERRIDE' in request.META:
            request.method = request.META['HTTP_X_HTTP_METHOD_OVERRIDE']

        request_method = self.method_check(request, allowed=allowed_methods)
        method = getattr(self, "%s_%s" % (request_method, request_type), None)

        if method is None:
            raise ImmediateHttpResponse(response=http.HttpNotImplemented())

        self.is_authenticated(request)
        self.throttle_check(request)

        # All clear. Process the request.
        tastypie.resources.convert_post_to_put(request)
        response = method(request, **kwargs)

        # Add the throttled request.
        self.log_throttled_access(request)

        # If what comes back isn't a ``HttpResponse``, assume that the
        # request was accepted and that some action occurred. This also
        # prevents Django from freaking out.
        if not isinstance(response, HttpResponseBase):
            return http.HttpNoContent()

        return response

    
from reports.utils.profile_decorator import profile

class UnlimitedDownloadResource(Resource):
    ''' 
    Resource that will stream the entire endpoint (with params):
    - instead of paginating it
    - instead of doing it in memory
    - will set a filename to the response header if specified
    '''
    #     @profile("unlimited_get_list.prof")
    def get_list(self, request, **kwargs):
        from reports.serializers import csv_convert
        from django.utils.encoding import smart_str
        
        try:
            limit = request.GET.get('limit', None)
            if limit == '0':
                logger.error(str(('using download specific get_list')))
                base_bundle = self.build_bundle(request=request)
                objects = self.obj_get_list(bundle=base_bundle, **self.remove_api_resource_names(kwargs))
                sorted_objects = self.apply_sorting(objects, options=request.GET)
            
                desired_format = self.determine_format(request)
                content_type=build_content_type(desired_format)
    
                response = HttpResponse(content_type=content_type)
                name = self._meta.resource_name
                response['Content-Disposition'] = \
                    'attachment; filename=%s.csv' % unicode(name).replace('.', '_')
                
                if desired_format == 'application/json':
                    pass;
                elif desired_format == 'application/xls':
                    pass;
                elif desired_format == 'text/csv':
                    csvwriter = csv.writer(response)
                    i = 0
                    keys = None
                    for obj in sorted_objects.iterator():
                        bundle = self.build_bundle(obj=obj, request=request)
                        bundle = self.full_dehydrate(bundle, for_list=True)
                        item = bundle.data
                        if i == 0:
                            keys = item.keys()
                            logger.info(str(('create the header row', [smart_str(key) for key in keys])))
                            csvwriter.writerow([smart_str(key) for key in keys])
                        i += 1
                        _list = []
                        for key in keys:
                            if key in item:
                                _list.append(csv_convert(item[key]))
                        csvwriter.writerow(_list)
                        if i % 1000 == 0:
                            logger.info(str(('logged', i)))
                    logger.error(str(('return response')))
                else:
                    raise Exception(str(('unknown format', desired_format)))
                
                
    #             serialized = self.serialize(request, data, desired_format)
    
                return response
            
            else:
                logger.info('using superclass get_list')
                return super(UnlimitedDownloadResource, self).get_list(request, **kwargs)

        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list()', self._meta.resource_name, 
                msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e


    def get_list2(self, request, **kwargs):
        """
        Returns the entire collection specified by the URI & params
        """
        from reports.serializers import csv_convert
        from django.utils.encoding import smart_str
        
        limit = request.GET.get('limit', None)
        if limit == '0':
            logger.error(str(('using download specific get_list')))
            base_bundle = self.build_bundle(request=request)
            objects = self.obj_get_list(bundle=base_bundle, **self.remove_api_resource_names(kwargs))
            sorted_objects = self.apply_sorting(objects, options=request.GET)
            
            response = HttpResponse(mimetype='text/csv')
            name = self._meta.resource_name
            response['Content-Disposition'] = \
                'attachment; filename=%s.csv' % unicode(name).replace('.', '_')
    
            csvwriter = csv.writer(response)
            
            i = 0
            keys = None
            for obj in sorted_objects.iterator():
                bundle = self.build_bundle(obj=obj, request=request)
                bundle = self.full_dehydrate(bundle, for_list=True)
                item = bundle.data
                if i == 0:
                    keys = item.keys()
                    logger.info(str(('create the header row', [smart_str(key) for key in keys])))
                    csvwriter.writerow([smart_str(key) for key in keys])
                i += 1
                _list = []
                for key in item:
                    _list.append(csv_convert(item[key]))
                csvwriter.writerow(_list)
                if i % 1000 == 0:
                    logger.info(str(('logged', i)))
            logger.error(str(('return response')))
            
            return response
        else:
            logger.info('using superclass get_list')
            return super(UnlimitedDownloadResource, self).get_list(request, **kwargs)
        
def download_tmp_file(path, filename):
    """                                                                         
    Send a file through Django without loading the whole file into              
    memory at once. The FileWrapper will turn the file object into an           
    iterator for chunks of 8KB.                                                 
    """
    try:
        _file = file(_path)
        logger.debug(str(('download_attached_file',_path,_file)))
        wrapper = FileWrapper(_file)
        # use the same type for all files
        response = HttpResponse(wrapper, content_type='text/plain') 
        response['Content-Disposition'] = \
            'attachment; filename=%s' % unicode(filename)
        response['Content-Length'] = os.path.getsize(_path)
        return response
    except Exception,e:
        logger.error(str(('could not find attached file object for id', id, e)))
        raise e



# NOTE if using this class, must implement the "not implemented error" methods
# on Resource (these are implemented with ModelResource):
# detail_uri_kwargs
# get_obj_list
# apply_filters
# obj_get_list
# obj_get
# obj_create
# obj_update
# obj_delete
class ManagedResource(LoggingMixin):
    '''
    Uses the field and resource definitions in the Metahash store to determine 
    the fields to expose for a Resource
    '''
    resource_registry = {}
    
    def __init__(self, field_definition_scope='fields.metahash', **kwargs):
        self.resource = self._meta.resource_name
        self.scope = 'fields.' + self.resource
        self.field_definition_scope = field_definition_scope
        self.meta_bootstrap_fields = ['resource_uri']
        
        logger.debug(str(('---init resource', 
                          self.resource, self.scope, field_definition_scope)))
        
        ManagedResource.resource_registry[self.scope] = self;

        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
        metahash = MetaHash.objects.get_and_parse(
            scope=self.scope, 
            field_definition_scope=field_definition_scope)
        for key,fieldhash in metahash.items():
            if 'filtering' in fieldhash and fieldhash['filtering']:
                self.Meta.filtering[key] = ALL_WITH_RELATIONS
        
        for key,fieldhash in metahash.items():
            if 'ordering' in fieldhash and fieldhash['ordering']:
                self.Meta.ordering.append(key)
        
        super(ManagedResource,self).__init__(**kwargs)
        self.original_fields = deepcopy(self.fields)
        self.create_fields()
        
    # local method  
    def reset_field_defs(self, scope=None):
        if not scope:
            scope = self.scope
        if scope not in ManagedResource.resource_registry:
            msg = str((
                'resource for scope not found: ', scope, 
                'in resource registry',ManagedResource.resource_registry.keys(),
                'possible cause: resource not entered in urls.py' ))
            logger.warn(msg)
            raise Exception(msg)
        resource = ManagedResource.resource_registry[scope]

        resource.clear_cache();
                
        logger.info(str((
            '----------reset_field_defs, resource_name' , 
            resource._meta.resource_name, 'scope', scope, 'resource', resource )))
        resource.create_fields();
        resource.reset_filtering_and_ordering();
    
    def clear_cache(self):
        
        # provisional 20140825
        logger.info('clear cache')
        cache.delete(self._meta.resource_name + ':schema')
        self.field_alias_map = {}
        
        
    # local method    
    # TODO: allow turn on/off of the reset methods for faster loading.
    def reset_filtering_and_ordering(self):
        self._meta.filtering = {}
        self._meta.ordering = []
        metahash = MetaHash.objects.get_and_parse(scope=self.scope, clear=True)
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self._meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self._meta.ordering.append(key)
        logger.debug(str(('meta filtering', self._meta.filtering)))
    
    # local method
    def create_fields(self):
        
        logger.debug(str(('--create_fields', self._meta.resource_name, 
            self.scope, 'original fields', self.original_fields.keys() )))
        if hasattr(self._meta, 'bootstrap_fields'):
            logger.debug(str(('bootstrap fields', self._meta.bootstrap_fields)))

#         try:        
#             _fields = self.build_schema()['fields']
#         except Exception, e:
#             logger.info(str(('in create_fields: resource information not available', e)))
        _fields = MetaHash.objects.get_and_parse(scope=self.scope, 
                    field_definition_scope='fields.metahash', clear=True)

#         self.local_field_defs = _fields
        logger.debug(str(('managed fields to create', _fields.keys())))

        try:
            resource_def = MetaHash.objects.get(
                scope='resource', key=self._meta.resource_name)
            resource_definition = resource_def.model_to_dict(scope='fields.resource')
            
            ## Get the supertype fields
            # TODO: -could- get the schema from the supertype resource
            supertype = resource_definition.get('supertype', '')
            if supertype:
                logger.info(str(('========= supertype fields', self._meta.resource_name)))
                supertype_fields = deepcopy(
                    MetaHash.objects.get_and_parse(
                        scope='fields.' + supertype, field_definition_scope='fields.metahash'))
                logger.info(str(('======= supertype fields', supertype, supertype_fields.keys())))
                supertype_fields.update(_fields)
                _fields = supertype_fields
        except Exception, e:
            logger.info(str(('in create_fields: resource information not available',self._meta.resource_name, e)))

        ## build field alias table
        self.field_alias_map = {}
        for field_name, item in _fields.items():
            alias = item.get('alias', None)
            if alias:
                self.field_alias_map[alias] = field_name

        new_fields = {}
        for field_name, field_obj in self.original_fields.items():
            if field_name in _fields:
                new_fields[field_name] = deepcopy(field_obj)
            elif ( hasattr(self._meta, 'bootstrap_fields') 
                    and field_name in self._meta.bootstrap_fields ):
                logger.debug('====== bootstrapping field: ' + field_name)
                new_fields[field_name] = deepcopy(field_obj)
            elif field_name in self.meta_bootstrap_fields:
                new_fields[field_name] = deepcopy(field_obj)

        unknown_keys = set(_fields.keys()) - set(new_fields.keys())
        logger.debug(str(('managed keys not yet defined', unknown_keys)))
        for field_name in unknown_keys:
            field_def = _fields[field_name]
            if 'json_field_type' in field_def and field_def['json_field_type']:
                # TODO: use type to create class instances
                # JSON fields are read only because they are hydrated in the 
                # hydrate_json_field method
                if field_def['json_field_type'] == 'fields.BooleanField':
                    new_fields[field_name] = eval(field_def['json_field_type'])(
                        attribute=field_name,
                        readonly=True, blank=True, null=True, default=False ) 
                else:
                    new_fields[field_name] = eval(field_def['json_field_type'])(
                        attribute=field_name,readonly=True, blank=True, null=True) 
            elif 'linked_field_type' in field_def and field_def['linked_field_type']:
                new_fields[field_name] = eval(field_def['linked_field_type'])(
                    attribute=field_name,readonly=True, blank=True, null=True) 
            else:
                logger.debug('creating unknown field as a char: ' + field_name)
                new_fields[field_name] = fields.CharField(
                    attribute=field_name, readonly=True, blank=True, null=True)

        logger.debug(str((
            'resource', self._meta.resource_name, self.scope, 
            'create_fields done: fields created', new_fields.keys() )))
        self.fields = new_fields
        return self.fields

    def alias_item(self, item):
        if self.field_alias_map:
            new_dict =  dict(zip(
                (self.field_alias_map.get(key, key),val) 
                for (key,val) in item.items()))
            return new_dict
        else:
            return item 
        
    def create_aliasmapping_iterator(self, data):
        logger.info(str(('====testing data', data)))
        def data_generator(data):
            for item in data:
                yield alias_item(data)
            
        if self.field_alias_map:
            data_generator(data)
        else:
            return data 
        
    def build_schema(self):
        '''
        Override
        '''
        # FIXME: consider using the cache decorator or a custom memoize decorator?
        schema = cache.get(self._meta.resource_name + ":schema")
        if schema:
            #if logger.isEnabledFor(logging.DEBUG):
            #    logger.debug(str(('====schema:', self._meta.resource_name, schema['fields'])))
            return schema
        
        logger.debug('------build_schema: ' + self.scope)
        
        try:
            schema = {}
            schema['fields'] = deepcopy(
                MetaHash.objects.get_and_parse(
                    scope=self.scope, field_definition_scope='fields.metahash'))
            
            if 'json_field' in schema['fields']: 
                # because we don't want this serialized directly (see dehydrate)
                schema['fields'].pop('json_field')  
            # all TP resources have an resource_uri field
            if not 'resource_uri' in schema['fields']:
                schema['fields']['resource_uri'] = { 
                    'key': 'resource_uri','table': 'None', 'visibility':[] }
            # all TP resources have an id field
            if not 'id' in schema['fields']:
                schema['fields']['id'] = { 'key': 'id', 'table':'None', 'visibility':[] }
            
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            if isinstance(e, ImmediateHttpResponse):
                msg = str(e.response)
            logger.warn(str(('on build_schema()', self._meta.resource_name, 
                msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e
            
        try:
            ## FIXED: client can get the Resource definition from either the 
            ## schema (here), or from the Resource endpoint; 
            ## SO the "resource_definition" here is copied to the endpoint bundle.data
            resource_def = MetaHash.objects.get(
                scope='resource', key=self._meta.resource_name)
            schema['resource_definition'] = resource_def.model_to_dict(scope='fields.resource')
            # TODO: -could- get the schema from the supertype resource
            supertype = schema['resource_definition'].get('supertype', '')
            if supertype:
                supertype_fields = deepcopy(
                    MetaHash.objects.get_and_parse(
                        scope=self.scope, field_definition_scope='fields.metahash'))
                supertype_fields.update(schema['fields'])
                schema['fields'] = supertype_fields
            
            # content_types: all resources serve JSON and CSV
            content_types = schema['resource_definition'].get('content_types', None)
            if not content_types:
                content_types = []
            _temp = set(content_types)
            _temp.add('json')
            _temp.add('csv')
            schema['resource_definition']['content_types'] = list(_temp)
            
            # find the default table for the resource
            default_table = schema['resource_definition']['table']
            logger.info(str(('default_table', default_table)))
            for key,field in schema['fields'].items():
                if not field.get('table', None):
                    field['table'] = default_table
                    
#             # pick apart the ORM to find DB particulars for read optimizations
#             sql = schema.get('sql', {})
#             
#             linked_table_module = schema['resource_definition'].get('linked_table_module', '')
#             if linked_table_module:
#                 if '.' in linked_field_module:
#                     # Try to import.
#                     module_bits = linked_field_module.split('.')
#                     module_path, class_name = '.'.join(module_bits[:-1]), module_bits[-1]
# #                     logger.info(str(('====linked: ' , linked_field_module,'gives: ', module_path, class_name)))
#                     module = importlib.import_module(module_path)
#                     sql['db_table'] = getattr(module.objects._meta, 'db_table')
#                 else:
#                     # We've got a bare class name here, which won't work (No AppCache
#                     # to rely on). Try to throw a useful error.
#                     raise ImportError(
#                         "linked_field_module requires a Python-style path "
#                         "(<module.module.Class>) to lazy load related resources. "
#                         "Only given '%s'." % linked_field_module )
        except Exception, e:
            logger.info(str(('in create_fields: resource information not available',self._meta.resource_name, e)))
        
        logger.debug('------build_schema,done: ' + self.scope ) 
       
        cache.set(self._meta.resource_name + ':schema', schema)
        return schema
    
    def is_valid(self, bundle, request=None):
        """
        obj_create{ full_hydrate, save{ is_valid, save_related, save, save_m2m }}
         
        NOTE: not extending tastypie.Validation variations, since they don't do 
        much, and we need access to the meta inf here anyway.
        NOTE: since "is_valid" is called at save(), so post-hydrate, all of the 
        defined fields have already "hydrated" the data; 
        this is a fail-fast validation, essentially, and preempts this validation.
        Here we will validate contextually, based on information in the metahash;
        overridden in each resource for more specific needs.
        
        Performs a check on the data within the bundle (and optionally the
        request) to ensure it is valid.

        Should return a dictionary of error messages. If the dictionary has
        zero items, the data is considered valid. If there are errors, keys
        in the dictionary should be field names and the values should be a list
        of errors, even if there is only one.
        """
        
        schema = self.build_schema()
        fields = schema['fields']
        
        # cribbed from tastypie.validation.py - mesh data and obj values, then validate
        data = {}
        if bundle.obj.pk:
            data = model_to_dict(bundle.obj)
        if data is None:
            data = {}
        data.update(bundle.data)
        
        # do validations
        errors = defaultdict(list)
        
        
        for name, field in fields.items():
            keyerrors = []
            value = data.get(name, None)
            
            if field.get('required', False):
                logger.debug(str(('check required: ', name, value)))
                
                if not value:
                     keyerrors.append('required')
                if isinstance(value, basestring):
                    if len(value.strip()) == 0:
                        keyerrors.append('required')
                        
            if not value:
                if keyerrors:
                    errors[name] = keyerrors            
                continue
            
            ##FIXME: some vocab fields are not choices fields
            if 'choices' in field and field['choices']:
                logger.debug(str(('check choices: ', name, value, field['choices'])))
                if field['ui_type'] != 'Checkboxes':
                    if value not in field['choices']:
                        keyerrors.append(
                            str((value, 'is not one of', field['choices'])))
                else:
                    for x in value:
                        if x not in field['choices']:
                            keyerrors.append(
                                str((value, 'are not members of', field['choices'])))

            if 'regex' in field and field['regex']:
                logger.debug(str(('check regex: ', name, value, field['regex'] )))
                if not re.match(field['regex'], value):
                    msg = field.get('validation_message', None)
                    if not msg:
                        msg = str((value, 'failed to match the pattern', field['regex']))
                    keyerrors.append(msg)

            if keyerrors:
                errors[name] = keyerrors
                
        
        if errors:
            bundle.errors[self._meta.resource_name] = errors
            logger.warn(str((
                'bundle errors', bundle.errors, len(bundle.errors.keys()),
                'bundle_data', data)))

            return False
        return True
        
    def deserialize(self, *args, **kwargs):
        '''
        Because field alias mapping is specific to each resource, override the 
        deserialize method here to patch in an iterator that will use the alias
        map to convert posted fields to mapped fields.
         
        ## TODO: implement tests for this
        '''
        deserialized = super(ManagedResource, self).deserialize(*args, **kwargs)
        
        if self._meta.collection_name in deserialized: 
            # this is a list of data
            deserialized[self._meta.collection_name] = \
                self.create_aliasmapping_iterator(deserialized[self._meta.collection_name])
        else:   
            # this is a single item of data
            deserialized = self.alias_item(deserialized)
        return deserialized
        
    def dehydrate(self, bundle):
        ''' 
        Implementation hook method, override to augment bundle, post dehydrate
        by the superclass used here to do the "hydrate_json_field"
        '''
        if len(bundle.data) == 0 : return bundle
        
#         schema = self.build_schema()
#         _fields = schema['fields']
        _fields = MetaHash.objects.get_and_parse(
            scope=self.scope, field_definition_scope='fields.metahash')
        for key in [ 
                x for x,y in _fields.items() if y.get('json_field_type') ]:
            bundle.data[key] = bundle.obj.get_field(key);
        
        bundle.data['json_field'] = ''
        # json_field will not be part of the public API, it is for internal use
        bundle.data.pop('json_field') 
        # override the resource_uri, since we want to export the permanent composite key
#        bundle.data['resource_uri'] = 
#             self.build_resource_uri(self.resource, bundle.data) or 
#             bundle.data['resource_uri']
        
        return bundle
    
    # implementation hook, to deserialize the embedded json fields
    def hydrate_json_field(self, bundle):
        '''
        hydrate bundle data values that will be stuffed into the json_field
        -Note: as mentioned elsewhere, for the initial load of the 
        Metahash:fields, fields that are JSON fields (to be stuffed into 
        json_field) must be first defined as a record with a 
        scope='metahash:field'; then they can be set as attribute values on 
        other fields in the second step.
        '''
        
        json_obj = {}
        
        # FIXME: why is clear=True here?
        local_field_defs = MetaHash.objects.get_and_parse(
            scope=self.scope, field_definition_scope='fields.metahash', clear=True)
        
        # Use the tastypie field type that has been designated to convert each
        # field in the json stuffed field just like it were a real db field
        for key in [ 
            str(x) for x,y in local_field_defs.items() 
                if 'json_field_type' in y and y['json_field_type'] ]:
            if key not in self.fields:
                raise RuntimeError(str((
                    'for the resource', self._meta.resource_name, 
                    'the key to deserialize', key, 
                    'was not defined as a resource field: fields.', 
                    self.fields.keys() )))
            val = bundle.data.get(key,None)
            if val:
                try:
                    if hasattr(val, "strip"): # test if it is a string
                        val = self.fields[key].convert(
                            smart_text(val,'utf-8', errors='ignore'))
                    # test if it is a sequence
                    elif hasattr(val, "__getitem__") or hasattr(val, "__iter__"): 
                        val = [smart_text(x,'utf-8',errors='ignore') for x in val]
                    json_obj.update({ key:val })
                except Exception, e:
                    logger.error('ex', e)
                    extype, ex, tb = sys.exc_info()
                    formatted = traceback.format_exception_only(extype, ex)[-1]
                    msg = str((
                        'failed to convert', key, 'with value', val, 'message', 
                        formatted)).replace("'","")
                    if key in self.fields:
                        msg += str(('with tastypie field type', 
                                    type(self.fields[key]) ))
                    e =  RuntimeError, msg
                    logger.warn(str((
                        'throw', e, tb.tb_frame.f_code.co_filename, 
                        'error line', tb.tb_lineno)))
                    raise e
        bundle.data['json_field'] = json.dumps(json_obj);
        logger.debug(str(('--- hydrated:', bundle.data['json_field'])))
        return bundle;


    # override
    def obj_get(self, bundle, **kwargs):
        try:
            bundle = super(ManagedResource, self).obj_get(bundle, **kwargs);
            return bundle
        except Exception, e:
            extype, ex, tb = sys.exc_info()
            logger.warn(str((
                'throw', e, tb.tb_frame.f_code.co_filename, 'error line', 
                tb.tb_lineno, extype, ex)))
            msg = str((e));
            if hasattr(e, 'response'): 
                msg = str(e.response)
            logger.warn(str(('==ex on get',bundle.request.path,e, msg )))
            raise e
#             extype, ex, tb = sys.exc_info()
#             logger.warn(str((
#                 'throw', e, tb.tb_frame.f_code.co_filename, 'error line', 
#                 tb.tb_lineno, extype, ex)))
#             logger.warn(str(('==ex on get, kwargs', kwargs, e)))
#             raise e
# #             raise type(e), str((type(e), e,
# #                                 'request.path', bundle.request.path, kwargs))

    def _get_attribute(self, obj, attribute):
        '''
        get an attribute that is possibly defined by dot notation:
        so reagent.substance_id will first get the reagent, then the substance_id
        '''
        parts = attribute.split('.')
        
        current_obj = obj
        for part in parts:
            if hasattr(current_obj,part):
                current_obj = getattr(current_obj, part)
        
        return current_obj 
    
    def _get_hashvalue(self, dictionary, attribute):
        '''
        get an attribute that is possibly defined by dot notation:
        so reagent.substance_id will first get the reagent, then the substance_id
        '''
        parts = attribute.split('.')
        current_val = dictionary
        for part in parts:
            if part in current_val:
                current_val = current_val[part]
        
        return current_val
    
    def _handle_500(self, request, exception):
        logger.error(str(('handle_500 error', self._meta.resource_name, str(exception))))
        return super(ManagedResource, self)._handle_500(request, exception)
        
    # override
    def detail_uri_kwargs(self, bundle_or_obj):
        """
        Override resources.ModelResource
        Given a ``Bundle`` or an object (typically a ``Model`` instance),
        it returns the extra kwargs needed to generate a detail URI.

        By default, it uses the model's ``pk`` in order to create the URI.
        """
        if bundle_or_obj is None:
            return {}

        resource_name = self._meta.resource_name
        id_attribute = None
        try:
            schema = self.build_schema()
            if 'resource_definition' not in schema:
                self.clear_cache()
                schema = self.build_schema()
            resource = schema['resource_definition']
            
            # TODO: memoize
            # note use an ordered dict here so that the args can be returned as
            # a positional array for 
            kwargs = OrderedDict() 
            id_attribute = resource['id_attribute']
            for x in id_attribute:
                val = ''
                if isinstance(bundle_or_obj, Bundle):
#                     val = getattr(bundle_or_obj.obj,x)
                    val = self._get_attribute(bundle_or_obj.obj, x)
                else:
                    if hasattr(bundle_or_obj, x):
#                         val = getattr(bundle_or_obj,x)  
                        val = self._get_attribute(bundle_or_obj,x)  
                    elif isinstance(bundle_or_obj, dict):
#                         val = bundle_or_obj[x] # allows simple dicts
                        val = self._get_hashvalue(bundle_or_obj, x) # allows simple dicts
                    else:
                        raise Exception(str(('obj', type(obj), obj, 'does not contain', x)))
                if isinstance(val, datetime.datetime):
                    val = val.isoformat()
                else:
                    val = str(val)
                
                kwargs[x] = val
            
            return kwargs
            
        except Exception, e:
            #             extype, ex, tb = sys.exc_info()
            #             logger.warn(str((
            #                 'throw', e, tb.tb_frame.f_code.co_filename, 'error line', 
            #                 tb.tb_lineno, extype, ex)))
            #             logger.warn(str(('==cannot grab id_attribute', bundle_or_obj, id_attribute, e)))
            logger.warn(str(('cannot grab id_attribute for', resource_name,bundle_or_obj, id_attribute,  e)))
            
            try:
                if logger.isEnabledFor(logging.INFO):
                    logger.info(str((
                        'unable to locate resource: ', resource_name,
                        ' has it been loaded yet for this resource?',
                        'also note that this may not work with south, since model methods',
                        'are not available: ', e, 
                        'type', type(bundle_or_obj),
                        'bundle', bundle_or_obj,id_attribute,
                         )))
            except Exception, e:
                logger.info(str(('reporting exception', e)))
        # Fall back to base class implementation 
        # (using the declared primary key only, for ModelResource)
        # This is useful in order to bootstrap the ResourceResource
        logger.info(str(( 'use base class method for ', bundle_or_obj)))
        return super(ManagedResource,self).detail_uri_kwargs(bundle_or_obj)

    def get_via_uri(self, uri, request=None):
        '''
        Override the stock method to allow lookup of relative uri's:
        - a 'relative uri' - or 'local uri' is one that doesn't include the 
        api name ("v1" for instance), but rather, is truncated on the left, 
        so that "api/vi/resource_name/key1" becomes "resource_name/key1".  
        This is useful because input file records can have a shorter 
        "resource_uri" field.
        '''
        if self._meta.resource_name not in uri:
            raise Exception(str((
                'invalid URI', uri, 
                'must contain at least the resource name', 
                self._meta.resource_name)))
        
        if request and request.path:
            path = request.path
            # remove the parts after the api_name ("v1") because that part is 
            # the resource name, calling context may not be for this resource
            path = path[: path.find(
                self._meta.api_name)+len(self._meta.api_name)+1] 
            local_uri = uri
            if path not in local_uri:
                uri = path + local_uri
        
        return super(ManagedResource, self).get_via_uri(uri, request);

    def get_local_resource_uri(
            self, bundle_or_obj=None, url_name='api_dispatch_list'):
        '''
        special 'local' version of the uri - when creating the uri for 
        containment lists (user.permissionss, for example), convert 
        "reports/api/v1/permission/resource/read" to "permission/resource/read"
        '''
        uri = super(ManagedResource, self).get_resource_uri(
            bundle_or_obj=bundle_or_obj, url_name=url_name)
        return uri[uri.find(self._meta.resource_name):]
    
    def get_resource_uri(self,bundle_or_obj=None, url_name='api_dispatch_list'):
        uri = super(ManagedResource, self).get_resource_uri(
            bundle_or_obj=bundle_or_obj, url_name=url_name)
        return uri
        
    # implementation hook - URLS to match _before_ the default URLS
    # used here to allow the natural keys [scope, key] to be used
    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" 
        # [ any word, except "schema" ]
        # also note the double underscore "__" is because we also don't want to
        # match in the first clause. Don't want "schema" since that reserved
        # word is used by tastypie for the schema definition for the resource
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.-:]+)/(?P<key>[^/]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            # TODO: is this needed here on metahash? we aren't using just "key" 
            # as a key, which is what causes the conflict with "schema", so probably not
            url(r"^(?P<resource_name>%s)/(?P<key>((?=(schema))__|(?!(schema))[^/]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

    def create_response(self, request, data, response_class=HttpResponse, **response_kwargs):
        '''
        Override to set a Content-Disposition attachment header in the response;
        - this can be used by the browser to give the download a name.
        - TODO: this is an expedient solution to the download button in the 
        browser, a better solution will set the Content-Disposition from the client
        '''
        response = super(ManagedResource, self).create_response(request, data, response_class, **response_kwargs)

        if 'Content-Disposition' not in response:
            format = request.GET.get('format', None)
            if format and format.lower() != 'json':
                if format in self._meta.serializer.content_types:
                    desired_format = self._meta.serializer.content_types[format]
                    response['Content-Type'] = desired_format
                    response['Content-Disposition'] = \
                        'attachment; filename=%s.%s' % (self._meta.resource_name, format)
                

        
# FIXME: Setting a filename to the response:
# - only set this header on the Response if it has been set in the Request;
# otherwise, browser may save response as a file rather than loading it with JS
#         filename = self._meta.resource_name
#         if hasattr(response, '_headers'):
#             if 'content-type' in response._headers:
#                 content_type = response._headers['content-type'][1]
#                 if 'text/csv' in content_type:
#                     filename += '.csv'
#                 elif 'application/json' in content_type:
#                     filename += '.json'
#         else:
#             logger.warn(str(('no "_header" found in response',
#                              'could not determine extension', response)))
#         response['Content-Disposition'] = 'attachment; filename=%s' % filename
        return response
 
 
class ExtensibleModelResourceMixin(ModelResource):
    '''
    Tastypie ModelResource mixin that passes the full request/url parsing kwargs
    on to the underlying "get_obj_list method:    
    Tastypie sequence:
    Resource.get_list()-> 
        ModelResource.obj_get_list()->  (kwargs not passed further)
            ModelResource.build_filters()->
            ModelResource.apply_filters()->
                ModelResource.get_obj_list()
    Ordinarily, "get_obj_list" does not receive the kwargs from the base
    Resource class - it just returns a stock query for the model.  With this 
    modification, extra args can be passed to modfify the stock query (and 
    return extra columns, for instance).
    
    Note: TP is built for extensible hooks, but this is pointing to a custom
    Resource class implementation.
    '''

    # Override Resoure to enable (optional) kwargs to modify the base query
    # Provisional, move to base class
    # get_list->obj_get_list->build_filters->apply_filters->get_obj_list
    def obj_get_list(self, bundle, **kwargs):
        """
        A ORM-specific implementation of ``obj_get_list``.

        Takes an optional ``request`` object, whose ``GET`` dictionary can be
        used to narrow the query.
        """
        filters = {}

        if hasattr(bundle.request, 'GET'):
            # Grab a mutable copy.
            filters = bundle.request.GET.copy()

        # Update with the provided kwargs.
        filters.update(kwargs)
        applicable_filters = self.build_filters(filters=filters)

        try:
            # MODIFICATION: adding kwargs to the apply_filters call - sde
            logger.debug(str(('kwargs', kwargs)))
            _kwargs = kwargs
            if 'request' in _kwargs: 
                _kwargs = {}
            objects = self.apply_filters(bundle.request, applicable_filters, **_kwargs)
            
            if self._meta.resource_name == 'apilog':
                logger.debug(str(('kwargs', kwargs)))
                return self._meta.authorization.read_list(objects, bundle, **_kwargs)
            else:
                return self.authorized_read_list(objects, bundle)
        except ValueError, e:
            logger.warn(str(('on obj_get_list', e)))
            raise BadRequest(str(("Invalid resource lookup data provided (mismatched type).", e)))
        
        
    def apply_filters(self, request, applicable_filters, **kwargs): 
        '''
        Delegates to the parent ModelResource.apply_filters method.
        '''       
        query = self.get_object_list(request, **kwargs)
        return query.apply_filters(request, applicable_filters);


    def get_object_list(self, request, **kwargs):
        '''
        Override this method if the kwargs will be used to modify the base query
        returned by get_obj_list()
        '''
        return super(ExtensibleModelResourceMixin, self).get_object_list(request);


class FilterModelResourceMixin(ExtensibleModelResourceMixin):
    '''
    Tastypie ModelResource mixin to enable "exclude" filters as well as "include"
    filters.
    
    How:  Modifies the dict returned from build_filters;
    - the new dict contains a top level: "include" and "exclude" sub dicts.
    - "include is the same, "exlude" is for an exclude filter.
    '''
    
    def build_filters(self, filters=None):
        ''' 
        Override of Resource - 
        Enable excludes as well as regular filters.
        see https://github.com/toastdriven/django-tastypie/issues/524#issuecomment-34730169
        
        Note: call sequence: 
        get_list->obj_get_list->build_filters->apply_filters->get_obj_list
        
        '''
        if not filters:
            return filters
    
        applicable_filters = {}
    
        # Normal filtering
        filter_params = dict([(x, filters[x]) 
            for x in filter(lambda x: not x.endswith('__ne'), filters)])
        applicable_filters['filter'] = \
            super(FilterModelResourceMixin, self).build_filters(filter_params)
    
        # Exclude filtering
        exclude_params = dict([(x[:-4], filters[x]) 
            for x in filter(lambda x: x.endswith('__ne'), filters)])
        applicable_filters['exclude'] = \
            super(FilterModelResourceMixin, self).build_filters(exclude_params)
    
        return applicable_filters 
       
    def apply_filters(self, request, applicable_filters, **kwargs):
        ''' 
        Override of ModelResource - 
        "applicable_filters" now contains two sub dictionaries:
        - "filter" - normal filters
        - "exclude" - exclude filters, to be applied serially (AND'ed)
        '''

        query = self.get_object_list(request, **kwargs)
        
        f = applicable_filters.get('filter')
        if f:
            query = query.filter(**f)
            
        e = applicable_filters.get('exclude')
        if e:
            for exclusion_filter, value in e.items():
                query = query.exclude(**{exclusion_filter: value})

        return query 
    
class ManagedModelResource(FilterModelResourceMixin, 
                           ManagedResource, PostgresSortingResource):
    pass

class MetaHashResource(ManagedModelResource):
    '''
    This table serves as a relatively low volume triple/quad store;
    - values can be defined at will for the text-based, software-implemented
    json_value storage field:
      - these values are defined by "virtual" fields in the quad store, keyed by
        a "field.tablename" scope.
      - must define a json_field_type for any "virtual" field stored in the 
        json_value field.
    - key field
    - scope field; secondary key, for quad store usage patterns (table based)
    
    The fields of this class are determined by the ManagedResource mixin.
    '''
    
    class Meta:
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field_type', 
                            'json_field','linked_field_type']
        queryset = MetaHash.objects.filter(
            scope__startswith="fields.").order_by('scope','ordinal','key')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= UserGroupAuthorization()        
        ordering = []
        filtering = {} #{'scope':ALL, 'key':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 
        resource_name = 'metahash'

    def __init__(self, **kwargs):
        super(MetaHashResource,self).__init__(**kwargs)
    
    def obj_create(self, bundle, **kwargs):
        '''
        Override - because the metahash resource is both a resource and the 
        definer of json fields, reset_field_defs after each create/update, 
        in case, new json fields are defined,or in case ordering,filtering 
        groups are updated
        '''
        bundle = super(MetaHashResource, self).obj_create(bundle, **kwargs);
        if getattr(bundle.obj,'scope').find('fields') == 0: #'fields.metahash':
            self.reset_field_defs(getattr(bundle.obj,'scope'))
        return bundle

    def obj_update(self, bundle, **kwargs):
        bundle = super(MetaHashResource, self).obj_update(bundle, **kwargs);
        self.reset_field_defs(getattr(bundle.obj,'scope'))
        return bundle
    
    def is_valid(self, bundle, request=None):
        '''
        We need to override this to bypass when initializing
        '''
        if getattr(bundle.obj,'scope').find('fields') == 0: #'fields.metahash':
            return True
        
        return super(MetaHashResource, self).is_valid(bundle, request=request)
            

    def hydrate(self, bundle):
        bundle = super(MetaHashResource, self).hydrate(bundle);
        return bundle
    
    def build_schema(self):
        schema = super(MetaHashResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 
            'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
        return schema
        
    def build_key(self, resource_name, data):
        '''
        Override, because the metahash resource is special, and will always use
        a /scope/key/ as key
        '''    
        return data['scope'] + '/' + data['key']

class VocabulariesResource(ManagedModelResource):
    '''
    This resource extends the ManagedModelResource using a new table 
    (vocabularies) but has fields defined in the Metahash table.
    '''
    def __init__(self, **kwargs):
        super(VocabulariesResource,self).__init__(**kwargs)

    class Meta:
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field']
        queryset = Vocabularies.objects.all().order_by('scope', 'ordinal', 'key')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= UserGroupAuthorization() #SuperUserAuthorization()        
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name = 'vocabularies'
        max_limit = 10000
    
    def build_schema(self):
        schema = super(VocabulariesResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 
            'label': 'Vocabulary', 'searchColumn': 'scope', 'options': temp }
        return schema

class ResourceResource(ManagedModelResource):
    '''
    This resource extends the ManagedModelResource, uses the metahash table
    internally, and has fields defined in the Metahash table.
    '''
    def __init__(self, **kwargs):
        super(ResourceResource,self).__init__(
            field_definition_scope='fields.resource', **kwargs)

    class Meta:
        '''
        Note, does not need the 'json_field_type' since MetahashResource is 
        managing the fields
        '''
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field'] 
        queryset = MetaHash.objects.filter(
            scope='resource').order_by('key', 'ordinal', 'scope')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= SuperUserAuthorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='resource' 
    
    def build_schema(self):
        schema = super(ResourceResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('key')]
        schema['extraSelectorOptions'] = { 
            'label': 'Resource', 'searchColumn': 'key', 'options': temp }
        return schema
    
    def is_valid(self, bundle, request=None):
        '''
        We need to override this to bypass when initializing
        TODO: re-examine dehydrate, same issue there.
        '''
        try:
            return super(ResourceResource, self).is_valid(bundle, request=request)
        except ObjectDoesNotExist, e:
            # notify and bypass
            logger.warn(str(('Resources not defined', e, self._meta.resource_name)))
            
            return True;
            
    def dehydrate(self, bundle):
        bundle = super(ResourceResource,self).dehydrate(bundle)
        # Get the schema
        # FIXME: why is the resource registry keyed off of "field."+key ?
        resource = ManagedResource.resource_registry['fields.'+bundle.obj.key]
        if resource:
            ## FIXED: use the schema version of the "resource definition";
            ## This is because there are modifications ("table" and "content_types")
            ## applied to the resource_definition in build_schema
            _temp_schema = deepcopy(resource.build_schema());
            bundle.data = deepcopy(_temp_schema['resource_definition'])
#             del _temp_schema['resource_definition']
            bundle.data['schema'] = _temp_schema
        else:
            logger.error('no API resource found in the registry for ' + 
                         bundle.data['key'] + 
                         '.  Cannot build the schema for this resource.' )
        return bundle

    def obj_create(self, bundle, **kwargs):
        '''
        Override - because the metahash resource is both a resource and the 
        definer of json fields, reset_field_defs after each create/update, 
        in case, new json fields are defined,or in case ordering,filtering 
        groups are updated
        '''
        bundle = super(ResourceResource, self).obj_create(bundle, **kwargs);
        if getattr(bundle.obj,'scope').find('fields') == 0: #'fields.metahash':
            self.reset_field_defs(getattr(bundle.obj,'scope'))
        return bundle

    def obj_update(self, bundle, **kwargs):
        bundle = super(ResourceResource, self).obj_update(bundle, **kwargs);
        self.reset_field_defs(getattr(bundle.obj,'scope'))
        return bundle

class ApiLogAuthorization(UserGroupAuthorization):
    '''
    Specialized authorization, allows users to read logs for resources they are 
    authorized for.
    '''
    
    def read_list(self, object_list, bundle, ref_resource_name=None, **kwargs):
        if not ref_resource_name:
            ref_resource_name = self.resource_meta.resource_name;
        if self._is_resource_authorized(
            ref_resource_name, bundle.request.user, 'read'):
            return object_list

    def read_detail(self, object_list, bundle, ref_resource_name=None, **kwargs):
        if not ref_resource_name:
            ref_resource_name = self.resource_meta.resource_name;
        if self._is_resource_authorized(
            ref_resource_name, bundle.request.user, 'read'):
            return True


class ApiLogResource(ManagedModelResource):
    
    class Meta:
        queryset = ApiLog.objects.all().order_by(
            'ref_resource_name', 'username','date_time')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= ApiLogAuthorization() #Authorization()        
        ordering = []
        filtering = {'username':ALL, 'uri': ALL, 'ref_resource_name':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='apilog' 
        max_limit = 100000
    
    def __init__(self, **kwargs):
        self.scope = 'fields.apilog'
        super(ApiLogResource,self).__init__(**kwargs)

    def get_resource_uri(self,bundle_or_obj=None):
        '''
        for efficiency, return a localized URI:
        /apilog/k1/k2/.../kn
        '''
        parts = [self._meta.resource_name]
        if bundle_or_obj is not None:
            id_kwarg_ordered = self.detail_uri_kwargs(bundle_or_obj)
            parts.extend(id_kwarg_ordered.values())
        return '/'.join(parts)
        
    def dehydrate_child_logs(self, bundle):
        
        if bundle.obj.child_logs.exists():
            return len(bundle.obj.child_logs.all())
#             uris = list()
#             for child_log in bundle.obj.child_logs.all():
#                 uris.append(self.get_resource_uri(child_log)) 
#             return uris
        return None
    
#     def dehydrate_diff_keys(self, bundle):
#         diff_keys = bundle.obj.diff_keys
#         logger.info(str(('=== diff_keys: ', diff_keys)))
#         if diff_keys:
#             changed = json.loads(diff_keys)
#             logger.info(str(('=== diff_keys2: ', changed)))
#             return changed     
#     
#     def dehydrate_diffs(self, bundle):
#         if bundle.obj.diffs:
#             diffs = json.loads(bundle.obj.diffs)
#             logger.info(str(('=== diffs2: ', diffs)))
#             return diffs 
    
    def dehydrate_parent_log_uri(self, bundle):
        parent_log = bundle.obj.parent_log
        if parent_log:
            id_kwarg_ordered = self.detail_uri_kwargs(parent_log)
            return '/'.join(id_kwarg_ordered.values())
        return None
    
#     def dehydrate_added_keys(self, bundle):
#         
#         if ( not getattr(bundle.obj, 'added_keys', None) and
#              bundle.obj.listlog_set.exists()):
#             keys = list()
#             for x in bundle.obj.listlog_set.all():
#                 if x.key and x.ref_resource_name:
#                     keys.append("%s/%s" % (x.ref_resource_name,x.key) )
#                 elif x.uri:
#                     keys.append(x.uri)
#             return keys
#         else:
#             return bundle.obj.added_keys
    
    def build_schema(self):
        schema = super(ApiLogResource,self).build_schema()
        temp = [ x.key for x in 
            MetaHash.objects.all().filter(scope='resource').distinct('key')]
        schema['extraSelectorOptions'] = { 
            'label': 'Resource', 
            'searchColumn': 'ref_resource_name', 'options': temp }
        return schema

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)/children%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_apilog_childview'), name="api_dispatch_apilog_childview"),
            url((r"^(?P<resource_name>%s)/children/(?P<ref_resource_name>[\w\d_.\-:]+)"
                 r"/(?P<key>[\w\d_.\-\+: ]+)"
                 r"/(?P<date_time>[\w\d_.\-\+:]+)%s$")
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_apilog_childview2'), name="api_dispatch_apilog_childview2"),
            url((r"^(?P<resource_name>%s)/(?P<ref_resource_name>[\w\d_.\-:]+)"
                 r"/(?P<key>[\w\d_.\-\+: ]+)"
                 r"/(?P<date_time>[\w\d_.\-\+:]+)%s$")
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

    def dispatch_apilog_childview(self, request, **kwargs):
        kwargs['parent_log_id'] = kwargs.pop('id')
        return ApiLogResource().dispatch('list', request, **kwargs)    

    def dispatch_apilog_childview2(self, request, **kwargs):
        logger.info(str(('kwargs', kwargs)))
        ref_resource_name = kwargs.pop('ref_resource_name')
        key = kwargs.pop('key')
        date_time = kwargs.pop('date_time')
        logger.info(str(('childview2', ref_resource_name, key, date_time)))
        parent_log = ApiLog.objects.get(ref_resource_name=ref_resource_name, key=key, date_time=date_time)
        kwargs['parent_log_id'] = parent_log.id
        return ApiLogResource().dispatch('list', request, **kwargs)    


class UserResource(ManagedModelResource):

    username = fields.CharField('user__username', null=False, readonly=True)
    first_name = fields.CharField('user__first_name', null=False, readonly=True)
    last_name = fields.CharField('user__last_name', null=False, readonly=True)
    email = fields.CharField('user__email', null=False, readonly=True)
    is_staff = CsvBooleanField('user__is_staff', null=True, readonly=True)
    is_superuser = CsvBooleanField('user__is_superuser', null=True, readonly=True)

    usergroups = fields.ToManyField(
        'reports.api.UserGroupResource', 'usergroup_set', related_name='users', 
        blank=True, null=True)
    permissions = fields.ToManyField(
        'reports.api.PermissionResource', 'permissions', null=True) #, related_name='users', blank=True, null=True)
    
    all_permissions = fields.ListField(attribute='all_permissions', blank=True, null=True, readonly=True)

    is_for_group = fields.BooleanField(attribute='is_for_group', blank=True, null=True)

    def __init__(self, **kwargs):
        super(UserResource,self).__init__(**kwargs)

    class Meta:
        queryset = UserProfile.objects.all().order_by('username') 
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        
        # FIXME: should override UserGroupAuthorization, and should allow user to view
        # (1) record by default: their own.
        authorization = SuperUserAuthorization()
#         authorization= UserGroupAuthorization() #SuperUserAuthorization()        
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name = 'user'

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<username>((?=(schema))__|(?!(schema))[^/]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>((?=(schema))__|(?!(schema))[^/]+))/groups%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_user_groupview'), name="api_dispatch_user_groupview"),
            url(r"^(?P<resource_name>%s)/(?P<username>((?=(schema))__|(?!(schema))[^/]+))/permissions%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_user_permissionview'), name="api_dispatch_user_permissionview"),
            ]    

    def dispatch_user_groupview(self, request, **kwargs):
        # signal to include extra column
        return UserGroupResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_permissionview(self, request, **kwargs):
        # signal to include extra column
        return PermissionResource().dispatch('list', request, **kwargs)    
        
    def obj_update(self, bundle, skip_errors=False, **kwargs):
        bundle = super(UserResource, self).obj_update(bundle, **kwargs);
        
        # Update the auth.user 
        django_user = bundle.obj.user
        
        # TODO validate these fields
        django_user.first_name = bundle.data.get('first_name')
        django_user.last_name = bundle.data.get('last_name')
        django_user.email = bundle.data.get('email')
        # Note cannot update username
        django_user.save()
        
    def obj_create(self, bundle, **kwargs):
        bundle = super(UserResource, self).obj_create(bundle, **kwargs);
        return bundle

    def hydrate(self, bundle):
        ''' 
        Called by full_hydrate 
        sequence is obj_create->full_hydrate(hydrate, then full)->save
        
        Our custom implementation will create an auth_user for the input; so 
        there will be a reports_userprofile.user -> auth_user.
        '''
        bundle = super(UserResource, self).hydrate(bundle);
        
        # fixup the username; stock hydrate will set either, but if it's not 
        # specified, then we will use the ecommons        
        ecommons = bundle.data.get('ecommons_id')
        username = bundle.data.get('username')
        email=bundle.data.get('email')
        first_name=bundle.data.get('first_name')
        last_name=bundle.data.get('last_name')
        is_staff = self.is_staff.convert(bundle.data.get('is_staff'))
        
        if not username:
            username = ecommons;
        bundle.obj.username = username
        
        django_user = None
        try:
            django_user = bundle.obj.user
        except ObjectDoesNotExist, e:
            from django.contrib.auth.models import User as DjangoUser
            try:
                django_user = DjangoUser.objects.get(username=username)
            except ObjectDoesNotExist, e:
                # ok, will create
                pass;

        if django_user:            
            django_user.first_name = first_name
            django_user.last_name = last_name
            django_user.email = email
            django_user.save();
        else:
            django_user = DjangoUser.objects.create_user(
                username, 
                email=email, 
                first_name=first_name, 
                last_name=last_name)
            # NOTE: we'll use user.is_password_usable() to verify if the 
            # user has a staff/manual django password account
            #             logger.debug(str(('save django user', django_user)))
            # Note: don't save yet, since the userprofile should be saved first
            # django_user.save()
            # this has to be done to set the FK on obj; since we're the only
            # side maintaining this rel' with auth_user

        django_user.is_staff = is_staff
        bundle.obj.user=django_user 
        return bundle
    
    def is_valid(self, bundle):
        """
        Should return a dictionary of error messages. If the dictionary has
        zero items, the data is considered valid. If there are errors, keys
        in the dictionary should be field names and the values should be a list
        of errors, even if there is only one.
        """
        
        # cribbed from tastypie.validation.py:
        # - mesh data and obj values, then validate
        data = {}
        if bundle.obj.pk:
            data = model_to_dict(bundle.obj)
        if data is None:
            data = {}
        data.update(bundle.data)
        
        # do validations
        errors = defaultdict(list)
        
        # TODO: rework this to be driven by the metahash
        
        if not data.get('first_name'):
            errors['first_name'] = ['first_name must be specified']
        
        if not data.get('last_name'):
            errors['last_name'] = ['last_name must be specified']
        
        if not data.get('email'):
            errors['email'] = ['email must be specified']
        
        ecommons = data.get('ecommons_id')
        username = data.get('username')
        
        if ecommons and username and (ecommons != username) :
            errors['specify either username or ecommons, not both']
        elif ecommons:
            bundle.obj.username = ecommons;
            
        if errors:
            bundle.errors[self._meta.resource_name] = errors
            logger.warn(str(('bundle errors', bundle.errors, len(bundle.errors.keys()))))
            return False
        return True
        
    def get_object_list(self, request, groupname=None):
        ''' 
        Called immediately before filtering, actually grabs the (ModelResource) 
        query - 
        
        Override this and apply_filters, so that we can control the 
        extra column "is_for_group".  This extra column is present when 
        navigating to users from a usergroup; see prepend_urls.
        '''
        query = super(UserResource, self).get_object_list(request);
        
        logger.debug(str(('get_obj_list', groupname)))
        if groupname:
            query = query.extra(select = {
                'is_for_group': ( 
                    '(select count(*)>0 '
                    ' from reports_usergroup ug '
                    ' join reports_usergroup_users ruu on(ug.id=ruu.usergroup_id) '
                    ' where ruu.userprofile_id=reports_userprofile.id '
                    ' and ug.name like %s )' ),
              },
              select_params = [groupname] )
            query = query.order_by('-is_for_group','user__last_name', 'user__first_name')
        return query

    def apply_filters(self, request, applicable_filters, **kwargs):

        query = self.get_object_list(request, **kwargs)
        logger.info(str(('applicable_filters', applicable_filters)))
        filters = applicable_filters.get('filter')
        if filters:
            
            # Grab the groups/users filter out of the dict
            groups_filter_val = None
            for f in filters.keys():
                if 'usergroup' in f:
                    groups_filter_val = filters.pop(f)

            query = query.filter(**filters)
            
            # then add the groups filter back in
            if groups_filter_val:
                ids = [x.id for x in UserProfile.objects.filter(
                        usergroup__name__iexact=groups_filter_val)]
                query = query.filter(id__in=ids)
            
        e = applicable_filters.get('exclude')
        if e:
            groups_filter_val = None
            for x in e.keys():
                if 'usergroup' in x:
                    groups_filter_val = e.pop(x)
            for exclusion_filter, value in e.items():
                query = query.exclude(**{exclusion_filter: value})

            # then add the user/groups filter back in
            if groups_filter_val:
                ids = [x.id for x in UserProfile.objects.filter(
                        usergroup__name__iexact=groups_filter_val)]
                query = query.exclude(id__in=ids)

        return query                 

    def apply_sorting(self, obj_list, options):
        '''
        Override to exclude certain fields from the PostgresSortingResource
        ''' 
        options = options.copy()
        options['non_null_fields'] = ['groups','is_for_group','user__first_name', 'user__last_name'] 
        obj_list = super(UserResource, self).apply_sorting(obj_list, options)
        return obj_list

    def get_resource_uri(self, bundle_or_obj=None, url_name='api_dispatch_list'):
        ''' 
        Override - Either have to generate localized resource uri, or, 
        in client, equate localized uri with non-localized uri's (* see user.js,
        where _.without() is used).
        This modification represents the first choice
        '''
        return self.get_local_resource_uri(
            bundle_or_obj=bundle_or_obj, url_name=url_name)

    def dehydrate_permissions(self, bundle):
        uri_list = []
        P = PermissionResource()
        # todo https://docs.djangoproject.com/en/dev/topics/db/queries/#lookups-that-span-relationships        
         
        #         userprofile = bundle.obj.userprofile_set.all()[0]
        for p in bundle.obj.permissions.all():
            uri_list.append(P.get_local_resource_uri(p))
        return uri_list;
    
    def dehydrate_all_permissions(self, bundle):
        uri_list = set()
        uri_list.update(self.dehydrate_permissions(bundle))
        P = PermissionResource()
        for permission in [ permission
            for usergroup in bundle.obj.usergroup_set.all()
            for permission in usergroup.get_all_permissions()]:
            uri_list.add(P.get_local_resource_uri(permission))
        
        return list(uri_list)
       
    def dehydrate_usergroups(self, bundle):
        uri_list = []
        UR = UserGroupResource()
        for g in bundle.obj.usergroup_set.all():
            uri_list.append(UR.get_local_resource_uri(g))
        return uri_list;

       
class UserGroupResource(ManagedModelResource):
    
    # relational fields must be defined   
    permissions = fields.ToManyField(
        'reports.api.PermissionResource', 'permissions',related_name='groups', 
        null=True) #, related_name='users', blank=True, null=True)
    users = fields.ToManyField('reports.api.UserResource', 'users', 
        related_name='usergroups', blank=True, null=True)
    
    super_groups = fields.ToManyField('reports.api.UserGroupResource', 
        'super_groups', blank=True, null=True)
    sub_groups = fields.ToManyField('reports.api.UserGroupResource', 
        'sub_groups', blank=True, null=True)

    all_permissions = fields.ListField(attribute='all_permissions', 
        blank=True, null=True, readonly=True)
    all_users = fields.ListField(attribute='all_users', 
        blank=True, null=True, readonly=True)
    
    is_for_user = fields.BooleanField(attribute='is_for_user', blank=True, null=True)
    is_for_group = fields.BooleanField(attribute='is_for_group', blank=True, null=True)

    class Meta:
        queryset = UserGroup.objects.all();        
        
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization = UserGroupAuthorization() #SuperUserAuthorization()        

        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='usergroup' 
    
    def __init__(self, **kwargs):
        super(UserGroupResource,self).__init__(**kwargs)

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<name>((?=(schema))__|(?!(schema))[^/]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<name>((?=(schema))__|(?!(schema))[^/]+))/users%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_group_userview'), name="api_dispatch_group_userview"),
            url(r"^(?P<resource_name>%s)/(?P<name>((?=(schema))__|(?!(schema))[^/]+))/permissions%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_group_permissionview'), 
                name="api_dispatch_group_permissionview"),
            url(r"^(?P<resource_name>%s)/(?P<name>((?=(schema))__|(?!(schema))[^/]+))/supergroups%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_group_supergroupview'), 
                name="api_dispatch_group_supergroupview"),
            ]

    def dispatch_group_userview(self, request, **kwargs):
        # signal to include extra column
        kwargs['groupname'] = kwargs.pop('name')  
        return UserResource().dispatch('list', request, **kwargs)    
    
    def dispatch_group_permissionview(self, request, **kwargs):
        # signal to include extra column
        kwargs['groupname'] = kwargs.pop('name')  
        return PermissionResource().dispatch('list', request, **kwargs)       
   
    def dispatch_group_supergroupview(self, request, **kwargs):
        # signal to include extra column
        kwargs['groupname'] = kwargs.pop('name')  
        return UserGroupResource().dispatch('list', request, **kwargs)    

    def get_resource_uri(self, bundle_or_obj=None, url_name='api_dispatch_list'):
        '''Override to shorten the URI'''
        return self.get_local_resource_uri(
            bundle_or_obj=bundle_or_obj, url_name=url_name)
    
    def obj_create(self, bundle, **kwargs):
        bundle = super(UserGroupResource, self).obj_create(bundle=bundle, **kwargs)
        return bundle

    def hydrate(self, bundle):
        bundle = super(UserGroupResource, self).hydrate(bundle);
        return bundle;
    
    def build_schema(self):
        schema = super(UserGroupResource,self).build_schema()
        return schema
    
    def dehydrate_permission_list(self, bundle):
        permissions = [ [x.scope, x.key, x.type] 
                            for x in bundle.obj.permissions.all()]
        return permissions
        
    def dehydrate_user_list(self,bundle):
        users = []
        for user in bundle.obj.users.all():
             users.append(
                 '[ %s - %s %s ]' 
                 % (user.username, user.first_name, user.last_name))
        return users
    
    #     def dehydrate_users(self, bundle):
    #         uri_list = []
    #         U = UserResource()
    #         for user in bundle.obj.users.all():
    #              uri_list.append(U.get_local_resource_uri(
    #                  { 'screensaver_user_id':user.screensaver_user_id }))
    #         return uri_list;
        
    def dehydrate_permissions(self, bundle):
        uri_list = []
        P = PermissionResource()
        for p in bundle.obj.permissions.all():
            uri_list.append(P.get_local_resource_uri(p))
        return uri_list;
    
    def dehydrate_super_groups(self, bundle):
        '''
        shallow report of groups contained directly in this group
        '''
        uri_list = []
        for g in bundle.obj.super_groups.all():
            uri_list.append(self.get_local_resource_uri(g))
        
        return uri_list
        
    def dehydrate_sub_groups(self, bundle):
        '''
        shallow report of groups contained directly in this group
        '''
        uri_list = []
        for g in bundle.obj.sub_groups.all():
            uri_list.append(self.get_local_resource_uri(g))
        
        return uri_list
        
    def dehydrate_all_permissions(self, bundle):
        P = PermissionResource()
        return [P.get_local_resource_uri(permission) for permission in 
                    bundle.obj.get_all_permissions()]
    
    def dehydrate_all_users(self, bundle):
        U = UserResource()
        return [U.get_local_resource_uri(user) for user in 
                    bundle.obj.get_all_users()]

    
    def dehydrate(self,bundle):
        bundle.data['id'] = bundle.obj.id
        return bundle

    # ModelResource override
    # get_list->obj_get_list->build_filters->apply_filters->get_obj_list
    def get_object_list(self, request, **kwargs): #is_for_user=None):
        ''' 
        Called immediately before filtering, this method actually grabs the 
        (ModelResource) base query - 
        
        Because we are using the ExtensibleModelResourceMixin here, we are also
        getting the kwargs from the Request.GET.  We use extra kwargs,
        ("username") in this case, to signal that we want to include the extra 
        column, "is_for_user".  
        
        Note: This special case is served from the url:
        /user/<username>/groups.
        All of this is a convenience feature for the client code; having
        "is_for_user" allows for easy filtering on the client.
        '''
        logger.debug(str(('get_obj_list', kwargs)))
        query = super(UserGroupResource, self).get_object_list(request);
        
        if 'username' in kwargs:
            is_for_user = kwargs.pop('username')
            query = query.extra(select = {
                'is_for_user': ( 
                    '(select count(*)>0 '
                    ' from reports_userprofile up '
                    ' join reports_usergroup_users ruu on(up.id=ruu.userprofile_id) '
                    ' where ruu.usergroup_id=reports_usergroup.id '
                    ' and up.username = %s )' ),
              },
              select_params = [is_for_user] )
            query = query.order_by('-is_for_user', 'name')        
        if 'groupname' in kwargs:
            is_for_group = kwargs.pop('groupname')
            query = query.extra(select = {
                'is_for_group': ( 
                    '(select count(*)>0 '
                    ' from reports_usergroup ug '
                    ' join reports_usergroup_super_groups ugsg on(ug.id=ugsg.from_usergroup_id) '
                    ' where ugsg.to_usergroup_id=reports_usergroup.id '
                    ' and ug.name = %s )' ),
              },
              select_params = [is_for_group] )
            query = query.exclude(name=is_for_group)
            query = query.order_by('-is_for_group', 'name')
        return query
    
    def apply_filters(self, request, applicable_filters, **kwargs):
        '''
        ModelResource override - 
        because the FilterModelResource mixin is being used here, the 
        applicable_filters includes an extra level: {filter,excludes}
        '''
        logger.info(str(('apply_filters', applicable_filters, kwargs)))
        query = self.get_object_list(request, **kwargs)
        
        filters = applicable_filters.get('filter')
        if filters:

            # Grab the users filter out of the dict
            users_filter_val = None
            for x in filters.keys():
                if 'users' in x:
                    users_filter_val = filters.pop(x)

            query = query.filter(**filters)
            
            if users_filter_val:
                ids = [x.id for x in UserGroup.objects.filter(
                        users__username__iexact=users_filter_val)]
                query = query.filter(id__in=ids)          
            
        e = applicable_filters.get('exclude')
        if e:
            users_filter_val = None
            for x in e.keys():
                if 'users' in x:
                    users_filter_val = e.pop(x)
            
            for exclusion_filter, value in e.items():
                query = query.exclude(**{exclusion_filter: value})

            if users_filter_val:
                ids = [x.id for x in UserGroup.objects.filter(
                        users__username__iexact=users_filter_val)]
                query = query.exclude(id__in=ids)          

        return query

    def apply_sorting(self, obj_list, options):
        options = options.copy()
        
        # Override to exclude these fields in the PostgresSortingResource 
        options['non_null_fields'] = ['is_for_user', 'is_for_group'] 

        obj_list = super(UserGroupResource, self).apply_sorting(obj_list, options)
        return obj_list

    

class PermissionResource(ManagedModelResource):
    
    usergroups = fields.ToManyField(
        'reports.api.UserGroupResource', 'usergroup_set', 
        related_name='permissions', blank=True, null=True)
    users = fields.ToManyField(
        'reports.api.UserResource', 'userprofile_set', 
        related_name='permissions', blank=True, null=True)
    
    groups = fields.CharField(attribute='groups', blank=True, null=True)

    is_for_group = fields.BooleanField(
        attribute='is_for_group', blank=True, null=True)
    is_for_user = fields.BooleanField(
        attribute='is_for_user', blank=True, null=True)

    class Meta:
        # note: the queryset for this resource is actually the permissions
        queryset = Permission.objects.all().order_by('scope', 'key')

        # FIXME: creating a "groups" field that can be used to sort
        #         key = 'groups'
        #         if 'postgres' in lims.settings.DATABASES['default']['ENGINE'].lower():
        #             queryset = queryset.extra( select = {
        #               key: ( "( select array_to_string(array_agg(ug.name), ', ') " 
        #                      "  from reports_usergroup ug "
        #                      "  join reports_usergroup_permissions ugp "
        #                         " on(ug.id=ugp.usergroup_id) "
        #                      "  where ugp.permission_id=reports_permission.id)" )
        #             } ) 
        #         else:
        #             logger.warn(str((
        #                 '=========using the special sqllite lims.settings.DATABASES', 
        #                 lims.settings.DATABASES)))
        #         queryset = queryset.extra( select = {
        #           key: ( "( select group_concat(ug.name, ', ') " 
        #                  "  from reports_usergroup ug "
        #                  "  join reports_usergroup_permissions ugp "
        #                     " on(ug.id=ugp.usergroup_id) "
        #                  "  where ugp.permission_id=reports_permission.id)" )
        #         } ) 

        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= UserGroupAuthorization() #SuperUserAuthorization()        
        object_class = object
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        # note, use this so that the queryset fields are not all added by default
        includes = [] 
        always_return_data = True # this makes Backbone happy
        resource_name='permission' 
    
    def __init__(self, **kwargs):
        super(PermissionResource,self).__init__(**kwargs)
        
        # create all of the permissions on startup
        resources = MetaHash.objects.filter(
            Q(scope='resource')|Q(scope__contains='fields.'))
        query = self._meta.queryset._clone()
        permissionTypes = Vocabularies.objects.all().filter(
            scope='permission.type')
        for r in resources:
            found = False
            for perm in query:
                if perm.scope==r.scope and perm.key==r.key:
                    found = True
            if not found:
                logger.info(str(('initialize permission not found: ', 
                                 r.scope, r.key)))
                for ptype in permissionTypes:
                    p = Permission.objects.create(
                        scope=r.scope, key=r.key, type=ptype.key)
                    logger.info(str(('bootstrap created permission', p)))

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.\-:]+)/"
                 r"(?P<key>[\w\d_.\-\+:]+)/(?P<type>[\w\d_.\-\+:]+)%s$" ) 
                        % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]
    
    def get_object_list(self, request, **kwargs): #username=None, groupname=None):
        ''' 
        Called immediately before filtering, actually grabs the (ModelResource) 
        query - 
        Override this and apply_filters, so that we can control the 
        extra column "is_for_group":
        This extra column is present when navigating to permissions from a 
        usergroup; see prepend_urls.
        TODO: we could programmatically create the "is_for_group" column by 
        grabbing the entire queryset, converting to an array of dicts, and 
        adding this field    
        '''
        query = super(PermissionResource, self).get_object_list(request);
        if 'groupname' in kwargs:
            groupname = kwargs.pop('groupname')
            logger.info(str(('get_obj_list', groupname)))
            query = query.extra(select = {
                'is_for_group': (
                    '( select count(*)>0 '
                    '  from reports_usergroup ug '
                    '  join reports_usergroup_permissions rup '
                       '  on(ug.id=rup.usergroup_id) '
                    ' where rup.permission_id=reports_permission.id '
                    ' and ug.name = %s )' ),
              },
              select_params = [groupname] )
            query = query.order_by('-is_for_group')
        if 'username' in kwargs:
            username = kwargs.pop('username')
            query = query.extra(select = {
                'is_for_user': (
                    '( select count(*)>0 '
                    '  from reports_userprofile up '
                    '  join reports_userprofile_permissions rup '
                       '  on(up.id=rup.userprofile_id) '
                    ' where rup.permission_id=reports_permission.id '
                    ' and up.username = %s )' ),
              },
              select_params = [username] )
            query = query.order_by('-is_for_user')
        return query
    
    def apply_filters(self, request, applicable_filters, **kwargs):
        
        query = self.get_object_list(request, **kwargs)
        logger.info(str(('applicable_filters', applicable_filters)))
        filters = applicable_filters.get('filter')
        if filters:
            
            # Grab the groups/users filter out of the dict
            groups_filter_val = None
            users_filter_val = None
            for f in filters.keys():
                if 'groups' in f:
                    groups_filter_val = filters.pop(f)
                if 'userprofile' in f:
                    users_filter_val = filters.pop(f)

            query = query.filter(**filters)
            
            # then add the groups filter back in
            if groups_filter_val:
                ids = [x.id for x in Permission.objects.filter(
                        usergroup__name__iexact=groups_filter_val)]
                query = query.filter(id__in=ids)
            if users_filter_val:
                ids = [x.id for x in Permission.objects.filter(
                        userprofile__username__iexact=users_filter_val)]
                query = query.filter(id__in=ids)
            
        e = applicable_filters.get('exclude')
        if e:
            groups_filter_val = None
            users_filter_val = None
            for x in e.keys():
                if 'userprofile' in x:
                    users_filter_val = e.pop(x)
                if 'groups' in x:
                    groups_filter_val = e.pop(x)
            for exclusion_filter, value in e.items():
                query = query.exclude(**{exclusion_filter: value})

            # then add the user/groups filter back in
            if groups_filter_val:
                ids = [x.id for x in Permission.objects.filter(
                        usergroup__name__iexact=groups_filter_val)]
                query = query.exclude(id__in=ids)
            if users_filter_val:
                ids = [x.id for x in Permission.objects.filter(
                        userprofile__username__iexact=users_filter_val)]
                query = query.exclude(id__in=ids)

        return query         

    def apply_sorting(self, obj_list, options):
        options = options.copy()
        # Override to exclude this field in the PostgresSortingResource 
        options['non_null_fields'] = ['groups','is_for_group','users','is_for_user'] 
        obj_list = super(PermissionResource, self).apply_sorting(
            obj_list, options)
        return obj_list

    def get_resource_uri(self, bundle_or_obj=None, url_name='api_dispatch_list'):
        return self.get_local_resource_uri(
            bundle_or_obj=bundle_or_obj, url_name=url_name)
    
    def obj_get(self, bundle, **kwargs):
        ''' 
        basically, if a permission is requested that does not exist, 
        it is created
        '''
        try:
            return super(PermissionResource, self).obj_get(bundle, **kwargs)
        except ObjectDoesNotExist:
            logger.info(str(('create permission on the fly', kwargs)))
            p = Permission(**kwargs)
            p.save()
            return p
    
    def build_schema(self):
        schema = super(PermissionResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 
            'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
        return schema
        

KEY_QUERY_ALIAS_PATTERN = '_{key}'

class ManagedLinkedResource(ManagedModelResource):
    ''' store resource virtual fields in a related table
    '''

    def __init__(self, **kwargs):
        super(ManagedLinkedResource,self).__init__(**kwargs)
        self.linked_field_defs = None
            
    def get_linked_fields(self, scope=None):
        '''
        Generate the resource fields that will be stored in the linked table
        '''
        if not self.linked_field_defs:
            
            schema = self.build_schema()
            _fields = schema['fields']
            resource = schema['resource_definition']
            
            self.linked_field_defs = { x: _fields[x] 
                for x,y in _fields.items() 
                    if y.get('linked_field_type',None) }

            logger.debug(str(('lookup the module.model for each linked field', 
                self.linked_field_defs.keys() )))
            for key,field_def in self.linked_field_defs.items():
                
                # Turn off dehydration for any of these fields that correspond 
                # to the automatic ModelResource fields
                if key in self.fields:
                    self.fields[key].use_in = None
                
                linked_field_module = field_def.get('linked_field_module', None)
                if not linked_field_module:
                    linked_field_module = resource.get('linked_table_module', None)
                if not linked_field_module:
                    raise Exception(str((
                        'no "linked_field_module" found in the field def', 
                        field_def, 
                        'no "linked_table_module" found in the resource def', 
                        resource)))
                    
                if '.' in linked_field_module:
                    # Try to import.
                    module_bits = linked_field_module.split('.')
                    module_path, class_name = '.'.join(module_bits[:-1]), module_bits[-1]
                    module = importlib.import_module(module_path)
                else:
                    # We've got a bare class name here, which won't work (No AppCache
                    # to rely on). Try to throw a useful error.
                    raise ImportError(
                        "linked_field_module requires a Python-style path "
                        "(<module.module.Class>) to lazy load related resources. "
                        "Only given '%s'." % linked_field_module )
        
                module_class = getattr(module, class_name, None)
        
                if module_class is None:
                    raise ImportError(
                        "Module '%s' does not appear to have a class called '%s'."
                             % (module_path, class_name))
                else:
                    field_def['linked_field_model'] = module_class
                    field_def['meta_field_instance'] = \
                            MetaHash.objects.get(key=field_def['key'])

        return self.linked_field_defs
    
    @log_obj_create
    @transaction.atomic()
    def obj_create(self, bundle, **kwargs):
        logger.info(str(('=== obj_create', self._meta.resource_name, bundle.data)))
        
        bundle.obj = self._meta.object_class()

        for key, value in kwargs.items():
            setattr(bundle.obj, key, value)

        bundle = self.full_hydrate(bundle)
        
        # TODO: == make sure follows is in a transaction block
        ## OK if called in "patch list"; not in "put list", TP's implementation of
        ## "put list" has a "rollback" function instead of a tx block; so they
        ## are implementing their own tx: see:
        ## "Attempt to be transactional, deleting any previously created
        ##  objects if validation fails."
        
        bundle = self.save(bundle)

        logger.debug(str(('==== save_linked_fields', self.get_linked_fields().keys() )))
        
        simple_linked_fields = {
            k:v for (k,v) in self.get_linked_fields().items() if v.get('linked_field_module',None)}
        for key,item in simple_linked_fields.items():
#             logger.debug(str(('populating simple linked field', item)))
            linkedModel = item.get('linked_field_model')
            val = bundle.data.get(key,None)
            field = self.fields[key]
            if val:
                val = self._safe_get_field_val(key,field, val)
                if item['linked_field_type'] != 'fields.ListField':
                    linkedObj = linkedModel()
                    self._set_value_field(linkedObj, bundle.obj, item, val)
                else:
                    self._set_multivalue_field(linkedModel, bundle.obj, item, val)

        # complex fields: 
        # TODO: using a blank in 'linked_field_module' to indicate, this is abstruse
        complex_linked_fields = {
            k:v for (k,v) in self.get_linked_fields().items() if not v.get('linked_field_module',None)}
        if len(complex_linked_fields):
            # setup the linked model instance: some magic here - grab the model
            # from the *first* field, since -all- the complex fields have the same one
            linkedModel = complex_linked_fields.values()[0]['linked_field_model']
            linkedObj = linkedModel()
            setattr( linkedObj, complex_linked_fields.values()[0]['linked_field_parent'], bundle.obj)
            
            for key,item in complex_linked_fields.items():
                val = bundle.data.get(key,None)
                field = self.fields[key]
                if val:
                    val = self._safe_get_field_val(key,field, val)
                    setattr( linkedObj, item['linked_field_value_field'], val)
            linkedObj.save()
                
        
        return bundle
    
    def _safe_get_field_val(self, key,field, val):
        try:
            if hasattr(val, "strip"): # test if it is a string
                val = smart_text(val,'utf-8', errors='ignore')
                if isinstance( field, fields.ListField): 
                    val = (val,)
                val = field.convert(val)
            # test if it is a sequence - only string lists are supported
            elif hasattr(val, "__getitem__") or hasattr(val, "__iter__"): 
                val = [smart_text(x,'utf-8',errors='ignore') for x in val]
            return val
        except Exception, e:
            logger.error('ex', e)
            extype, ex, tb = sys.exc_info()
            formatted = traceback.format_exception_only(extype, ex)[-1]
            msg = str((
                'failed to convert', key, 'with value', val, 'message', 
                formatted)).replace("'","")
            if key in self.fields:
                msg += str(('with tastypie field type', 
                            type(field) ))
            e =  RuntimeError, msg
            logger.warn(str((
                'throw', e, tb.tb_frame.f_code.co_filename, 
                'error line', tb.tb_lineno)))
            raise e

    def _set_value_field(self, linkedObj, parent, item, val):
        ## TODO: updates should be able to set fields to None
        
#         logger.debug(str(('_set_value_field', linkedObj, parent, item['key'], val)))
        setattr( linkedObj, item['linked_field_parent'], parent)
        
        if item.get('linked_field_meta_field', None):
            setattr( linkedObj,item['linked_field_meta_field'], item['meta_field_instance'])
        
        setattr( linkedObj, item['linked_field_value_field'], val)
        linkedObj.save()

    def _set_multivalue_field(self, linkedModel, parent, item, val):
        logger.info(str(('_set_multivalue_field', item['key'], linkedModel, parent, item, val)))
        if isinstance(val, six.string_types):
            val = (val) 
        for i,entry in enumerate(val):
            linkedObj = linkedModel()
            setattr( linkedObj, item['linked_field_parent'], parent)
            if item.get('linked_field_meta_field', None):
                setattr( linkedObj,item['linked_field_meta_field'], item['meta_field_instance'])
            setattr( linkedObj, item['linked_field_value_field'], entry)
            if hasattr(linkedObj, 'ordinal'):
                linkedObj.ordinal = i
            linkedObj.save()
    
    @log_obj_update
    def obj_update(self, bundle, skip_errors=False, **kwargs):
        """
        A linked_field specific version
        """
        bundle = self._locate_obj(bundle)
        
        bundle = self.full_hydrate(bundle)
        self.is_valid(bundle)
        if bundle.errors and not skip_errors:
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))

        bundle.obj.save()
        
        logger.info(str(('==== update_linked_fields', self.get_linked_fields().keys() )))

        simple_linked_fields = {
            k:v for (k,v) in self.get_linked_fields().items() if v.get('linked_field_module',None)}
        for key,item in simple_linked_fields.items():
#             logger.debug(str(('populating simple linked field', item)))
            val = bundle.data.get(key,None)
            field = self.fields[key]
            
            if val:
                val = self._safe_get_field_val(key,field, val)
                linkedModel = item.get('linked_field_model')

                params = { item['linked_field_parent']: bundle.obj }
                if item.get('linked_field_meta_field', None):
                    params[item['linked_field_meta_field']] = item['meta_field_instance']

                if item['linked_field_type'] != 'fields.ListField':
                    linkedObj = None
                    try: 
                        linkedObj = linkedModel.objects.get(**params)
                    except ObjectDoesNotExist:
                        logger.warn(str((
                            'update, could not find extant linked field for', 
                            bundle.obj, item['meta_field_instance'])))
                        linkedObj = linkedModel()
                    
                    self._set_value_field(linkedObj, bundle.obj, item, val)
                else:
                    query = linkedModel.objects.filter( **params ) #.order_by('ordinal')
                    values = query.values_list(
                            item['linked_field_value_field'], flat=True)
                    if values == val:
                        pass
                    else:
                        query.delete()
                    self._set_multivalue_field(linkedModel, bundle.obj, item, val)
        
        # complex fields: 
        # TODO: using a blank in 'linked_field_module' to indicate, this is abstruse
        complex_linked_fields = {
            k:v for (k,v) in self.get_linked_fields().items() if not v.get('linked_field_module',None)}
        if len(complex_linked_fields):
            # setup the linked model instance: some magic here - grab the model
            # from the *first* field, since -all- the complex fields have the same one
            linkedModel = complex_linked_fields.values()[0]['linked_field_model']
            linkedObj = None
            try: 
                linkedObj = linkedModel.objects.get(**{ 
                    item['linked_field_parent']: bundle.obj })
            except ObjectDoesNotExist:
                logger.warn(str((
                    'update, could not find extant linked complex module for', 
                    bundle.obj)))
                linkedObj = linkedModel()
                setattr( linkedObj, item['linked_field_parent'], bundle.obj)
            
            for key,item in complex_linked_fields.items():
#                 logger.debug(str(('populating complex linked field', item)))
                val = bundle.data.get(key,None)
                field = self.fields[key]
                if val:
                    val = self._safe_get_field_val(key,field, val)
                    setattr( linkedObj, item['linked_field_value_field'], val)
            linkedObj.save()
        
        return bundle
                
    @log_obj_delete
    def obj_delete(self, bundle, **kwargs):
        if not hasattr(bundle.obj, 'delete'):
            try:
                bundle.obj = self.obj_get(bundle=bundle, **kwargs)
            except ObjectDoesNotExist:
                raise NotFound("A model instance matching the provided arguments could not be found.")

        self.authorized_delete_detail(self.get_object_list(bundle.request), bundle)

        # TODO: TEST!
        logger.info(str(('==== delete_linked_fields', self.get_linked_fields().keys() )))
        linkedModel = item.get('linked_field_model')
        linkedModel.objects.filter(**{
            item['linked_field_parent']: bundle.obj }).delete()
        bundle.obj.delete()

    def full_dehydrate(self, bundle, for_list=False):
        # trigger get_linked_fields to turn off "use_in" for model fields
        self.get_linked_fields()
        bundle =  ManagedModelResource.full_dehydrate(self, bundle, for_list=for_list)
        return bundle
    
    def get_object_list(self, request):
        query = super(ManagedLinkedResource,self).get_object_list(request)
        
        # FIXME: SQL injection attack through metadata.  (at least to avoid inadvertant actions)
        # TODO: use SqlAlchemy http://docs.sqlalchemy.org/en/latest/core/expression_api.html
        extra_select = OrderedDict()
        # NOTE extra_tables cannot be used, because it creates an inner join
        #         extra_tables = set()
        extra_where = []
        extra_params = []
        for key,item in self.get_linked_fields().items():
            key_query_alias = KEY_QUERY_ALIAS_PATTERN.format(key=key)
            
            field_name = item.get('field_name', None)
            if not field_name:
                field_name = item.get('linked_field_value_field',key)
            
            linkedModel = item.get('linked_field_model')
            field_table = linkedModel._meta.db_table
            
            parent_table = query.model._meta.db_table
            parent_table_key = query.model._meta.pk.name

            format_dict = {
                'field_name': field_name, 
                'field_table': field_table,
                'linked_field_parent': item['linked_field_parent'],
                'parent_table': parent_table,
                'parent_table_key': parent_table_key }
            
            if item['linked_field_type'] != 'fields.ListField':
                sql = ( 'select {field_name} from {field_table} {where}')
                where = ('WHERE {field_table}.{linked_field_parent}_id'
                            '={parent_table}.{parent_table_key} ')
                
                if item.get('linked_field_meta_field', None):
                    format_dict['meta_field'] = item['linked_field_meta_field']
                    meta_field_id = getattr(item['meta_field_instance'], 'pk')
                    where += ' and {field_table}.{meta_field}_id=%s '
                    extra_params.append(meta_field_id)
#                 else:
#                     extra_select[key_query_alias] = \
#                         '%s.%s' % (field_table, item['linked_field_value_field'])
#                     extra_tables.add(field_table)
#                     extra_where.append(
#                         '{field_table}.{linked_field_parent}_id'
#                             '={parent_table}.{parent_table_key}'.format(**format_dict))
                format_dict['where'] = where.format(**format_dict)
                sql = sql.format(**format_dict)
                extra_select[key_query_alias] = sql
            if item['linked_field_type'] == 'fields.ListField':
                sql = \
'''    (select $$["$$ || array_to_string(array_agg({field_name}), $$","$$) || $$"]$$
        from (select {field_name} from {field_table} 
        {where} {order_by}) a) '''
                
                where = ('WHERE {field_table}.{linked_field_parent}_id'
                            '={parent_table}.{parent_table_key} ')

                if item.get('linked_field_meta_field', None):
                    format_dict['meta_field'] = item['linked_field_meta_field']
                    meta_field_id = getattr(item['meta_field_instance'], 'pk')
                    where += ' and {field_table}.{meta_field}_id=%s '
                    extra_params.append(meta_field_id)
                else:
                    pass
                
                format_dict['order_by'] = ''
                ordinal_field = item.get('ordinal_field', None)
                if ordinal_field:
                    format_dict['order_by'] = ' ORDER BY %s ' % ordinal_field
                format_dict['where'] = where.format(**format_dict)
                sql = sql.format(**format_dict)
                extra_select[key_query_alias] = sql

        query = query.extra(
            select=extra_select, where=extra_where,select_params=extra_params )
        logger.info(str(('==== query', query.query.sql_with_params())))
        return query
     
    def dehydrate(self, bundle):
        try:
            keys_not_available = []
            for key,item in self.get_linked_fields().items():
                key_query_alias = KEY_QUERY_ALIAS_PATTERN.format(key=key)
                bundle.data[key] = None
                if hasattr(bundle.obj, key_query_alias):
                    bundle.data[key] = getattr(bundle.obj, key_query_alias)
                    if bundle.data[key] and item['linked_field_type'] == 'fields.ListField':
                        bundle.data[key] = json.loads(bundle.data[key])
                else:
                    keys_not_available.append(key_query_alias)
            if keys_not_available:
                logger.error(str(('keys not available', keys_not_available)))
            return bundle
        except Exception, e:
            extype, ex, tb = sys.exc_info()
            msg = str(e)
            if isinstance(e, ImmediateHttpResponse):
                msg = str(e.response)
            logger.warn(str((
                'throw', e, msg, tb.tb_frame.f_code.co_filename, 'error line', 
                tb.tb_lineno, extype, ex)))
            raise e
        return bundle
            
    def dehydrate_inefficient(self, bundle):
        '''
        Note - dehydrate only to be used for small sets. 
        Looks each of the the fields up as separare query; should be able to modify
        the parent query and find these fields in it using ORM methods 
        (TODO: we will implement obj-get-list methods)
        '''
        try:
            for key,item in self.get_linked_fields().items():
                bundle.data[key] = None
                linkedModel = item.get('linked_field_model')
                queryparams = { item['linked_field_parent']: bundle.obj }
                if item.get('linked_field_meta_field', None):
                    queryparams[item['linked_field_meta_field']] = item['meta_field_instance']
                if item['linked_field_type'] != 'fields.ListField':
                    try:
                        linkedObj = linkedModel.objects.get(**queryparams)
                        bundle.data[key] = getattr( linkedObj, item['linked_field_value_field'])
                    except ObjectDoesNotExist:
                        pass
                else:
                    query = linkedModel.objects.filter(**queryparams)
                    if hasattr(linkedModel, 'ordinal'):
                        query = query.order_by('ordinal')
                    values = query.values_list(
                            item['linked_field_value_field'], flat=True)
    #                 logger.debug(str((key,'multifield values', values)))
                    if values and len(values)>0:
                        bundle.data[key] = list(values)
            return bundle
        except Exception, e:
            extype, ex, tb = sys.exc_info()
            msg = str(e)
            if isinstance(e, ImmediateHttpResponse):
                msg = str(e.response)
            logger.warn(str((
                'throw', e, msg, tb.tb_frame.f_code.co_filename, 'error line', 
                tb.tb_lineno, extype, ex)))
        return bundle
    
    
class RecordResource(ManagedLinkedResource):
    ''' poc: store resource virtual fields in a related table
    '''
    class Meta:
        queryset = Record.objects.all()
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= SuperUserAuthorization()        

        ordering = []
        filtering = {'scope':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='record' 

    def __init__(self, **kwargs):
        super(RecordResource,self).__init__(**kwargs)
            

