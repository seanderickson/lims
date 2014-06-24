import datetime
import json
import logging
import sys
import os
import traceback
from collections import defaultdict, OrderedDict
from copy import deepcopy
from django.conf.urls import url
from django.conf import settings
from django.utils.encoding import smart_text
from django.utils import timezone
from django.forms.models import model_to_dict
from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.models import User
from django.db.models.aggregates import Max
from django.db.models import Q
from django.db import transaction
from tastypie.exceptions import NotFound, ImmediateHttpResponse, Unauthorized,\
    BadRequest
from tastypie.bundle import Bundle
from tastypie.authorization import Authorization, ReadOnlyAuthorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication,\
    MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
# NOTE: tastypie.fields is required for dynamic field instances using eval
from tastypie import fields 
from tastypie.resources import Resource, ModelResource
from tastypie.utils.urls import trailing_slash

from reports.serializers import LimsSerializer, CsvBooleanField
from reports.models import MetaHash, Vocabularies, ApiLog, Permission, \
                           UserGroup, UserProfile
# import lims.settings 
from tastypie.utils.timezone import make_naive
from django.db.utils import IntegrityError
        
logger = logging.getLogger(__name__)


# TODO: this is only a placeholder
class UserGroupAuthorization(Authorization):
    def read_list(self, object_list, bundle):
        
        user = bundle.request.user
        if user.is_superuser:
            return object_list
         
        userprofile = user.userprofile
        resource_name = self.resource_meta.resource_name
        
        permissions = [x for x in 
            userprofile.permissions.all().filter(scope='resource', key=resource_name)]
        logger.info(str(('user permissions', permissions)))
        permissions_group = [ permission 
                for group in userprofile.usergroup_set.all() 
                for permission in group.permissions.all()]
        logger.info(str(('group permissions', permissions)))
        
        
        raise Unauthorized("User requires permission to read this resource.")

    def read_detail(self, object_list, bundle):
        return True

    def create_list(self, object_list, bundle):
        return []

    def create_detail(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return True
        raise Unauthorized("Only superuser may create.")

    def update_list(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return object_list
        raise Unauthorized("Only superuser may update lists.")

    def update_detail(self, object_list, bundle):
        if bundle.request.user.is_superuser:
            return True
        raise Unauthorized("Only superuser may update.")

    def delete_list(self, object_list, bundle):
        return []

    def delete_detail(self, object_list, bundle):
        raise Unauthorized("You are not allowed to access that resource.")



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

        logger.debug(str(('==== non null fields', non_null_fields))) 
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

    
#     def apply_sorting1(self, obj_list, options):
#         """
#         Create a none-too-pretty workaround for the postgresql null sorting
#         issue - nulls sort higher than values, which is not desired.  
#         We want nulls to sort lower than values.
#         
#         Caveat: this will not work with joined fields unless they have an alias.  
#         This is because it creates a field like:
#         (screensaver_user_id is null) AS "screensaver_user_id_null"
#         - if this field is duplicated in two sides of a join, then it must be 
#         referenced by an alias, or as "table".screensaver_user_id, 
#         and we are not supporting table speciciations in this method, so if 
#         joined fields are used, they must be referenced by alias only.
# 
#         @param non_null_fields list - fields to ignore
#         """ 
#         obj_list = super(PostgresSortingResource, self).apply_sorting(
#             obj_list, options)
#         logger.debug(str(('order_by', obj_list.query.order_by)))
#         extra_select = {}
#         extra_ordering = []
#         
#         non_null_fields = options.get('non_null_fields', [])
#         logger.debug(str(('==== non null fields', non_null_fields))) 
#         for field in obj_list.query.order_by:
#             is_null_dir = '-'  # default nulls first for ascending
#             if field.startswith('-'):
#                 is_null_dir = ''
#                 field = field[1:]
#             if field in non_null_fields:
#                 continue
#             extra_select[field+"_null"]=field + ' is null'
#             extra_ordering.append(is_null_dir + field+"_null")
#         logger.debug(str(('extra_select', extra_select, 
#                           'extra_ordering', extra_ordering)))
#         obj_list = obj_list.extra(extra_select)
# 
#         # Note: the following doesn't work, something in the framework 
#         # deletes the extra order_by clause when apply_sorting, or, if this is
#         # run last, it deletes the sorting applied in apply_sorting...
#         #        obj_list = obj_list.extra(order_by=['-comments_null'])
# 
#         # Note: this doesn't work because the "is null" field order by clauses
#         # must be prepended so that they occur before their intended fields
#         #        obj_list.query.add_ordering('comments_null')
#         
#         temp = obj_list.query.order_by;
#         obj_list.query.clear_ordering(force_empty=True)
#         temp1 = []
#         for xfield in extra_ordering:
#             temp1.append(xfield)
#         
#         temp1.extend(temp)
#         logger.debug(str(('ordering', temp1)))
#         obj_list.query.add_ordering(*temp1)
#         
#         return obj_list



class LoggingMixin(Resource):
    '''
    intercepts obj_create, obj_update and creates an ApiLog entry for the action
    Note: whatever is being extended with the LoggingMixin must also define a
    "detail_uri_kwargs" method that returns an _ordered_dict_, since we log the 
    kwargs as ordered args.
    ** note: "detail_uri_kwargs" returns the set of lookup keys for the resource 
    URI construction.
    '''
        
    @transaction.commit_on_success()
    def patch_list(self, request, **kwargs):
        ''' Override
        '''
        # create an apilog for the patch list
        listlog = self.listlog = ApiLog()
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
        
        
        response =  super(LoggingMixin, self).patch_list(request, **kwargs) 
        
        listlog.save();
        listlog.key = listlog.id
        listlog.save()
        self.listlog = None

        return response        
    
    
    @transaction.commit_on_success()
    def obj_create(self, bundle, **kwargs):
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('----log obj_create', bundle)))
        
        bundle = super(LoggingMixin, self).obj_create(bundle=bundle, **kwargs)
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('object created', bundle.obj )))
        log = ApiLog()
        log.username = bundle.request.user.username 
        log.user_id = bundle.request.user.id 
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
        #        log.diffs = json.dumps(bundle.obj)
            
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
# 
#         if 'apilog_comment' in bundle.data:
#             log.comment = bundle.data['apilog_comment']
            
        log.save()
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(('create, api log', log)) )


        # TODO: create an analog of this in delete, update
        # if there is a listlog, it means "patch_list", or "put_list" were called
        if hasattr(self, 'listlog') and self.listlog:
            if(logger.isEnabledFor(logging.DEBUG)):
                logger.debug(str(('update listlog', self.listlog)))
            added_keys = []
            if self.listlog.added_keys:
                added_keys = json.loads(self.listlog.added_keys)
            added_keys.append(log.key); # TODO: append the log.id too?
            self.listlog.added_keys = json.dumps(added_keys)
        
        return bundle    
    
    # TODO: not tested
    @transaction.commit_on_success()
    def obj_delete(self, bundle, **kwargs):
        logger.info('---log obj_delete')
        
        super(LoggingMixin, self).obj_delete(bundle=bundle, **kwargs)
        
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
#         if 'apilog_comment' in bundle.data:
#             log.comment = bundle.data['apilog_comment']
            
        log.save()
        logger.info(str(('delete, api log', log)) )
        
        return bundle
    
    
    def _locate_obj(self, bundle, **kwargs):
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

    def compare_dicts(self, dict1, dict2, excludes=[]):
        original_keys = set(dict1.keys())-set(excludes)
        updated_keys = set(dict2.keys())-set(excludes)
        
        intersect_keys = original_keys.intersection(updated_keys)
        log = {}
        
        added_keys = list(updated_keys - intersect_keys)
        if len(added_keys)>0: 
            log['added_keys'] = added_keys
        
        removed_keys = list(original_keys- intersect_keys)
        if len(removed_keys)>0: 
            log['removed_keys'] = removed_keys
        
        diff_keys = list()
        
        s = self._meta.serializer
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
            log['diffs'] = dict(
                zip(diff_keys, 
                    ([dict1[key],dict2[key]] 
                        for key in diff_keys )  ))
        
        return log
    
    @transaction.commit_on_success()
    def obj_update(self, bundle, skip_errors=False, **kwargs):
        logger.info('--- log obj_update')
        original_bundle = Bundle(data=deepcopy(bundle.data))
        i=0;

        if hasattr(bundle,'obj'): original_bundle.obj = bundle.obj
        original_bundle = self._locate_obj(original_bundle, **kwargs)
        
        # store and compare dehydrated outputs: 
        # the api logger is concerned with what's sent out of the system, i.e.
        # the dehydrated output, not the internal representations.
        original_bundle = super(LoggingMixin, self).full_dehydrate(original_bundle)
        updated_bundle = super(LoggingMixin, self).obj_update(bundle=bundle, **kwargs)
        updated_bundle = super(LoggingMixin, self).full_dehydrate(updated_bundle)

        original_keys = set(original_bundle.data.keys())
        updated_keys = set(updated_bundle.data.keys())
        
        intersect_keys = original_keys.intersection(updated_keys)
        
        log = ApiLog()
        log.username = bundle.request.user.username 
        log.user_id = bundle.request.user.id 
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
        log.uri = self.get_resource_uri(bundle)
        log.key = '/'.join([str(x) for x in self.detail_uri_kwargs(bundle).values()])
        
        added_keys = list(updated_keys - intersect_keys)
        if len(added_keys)>0: 
            log.added_keys = json.dumps(added_keys)
        
        removed_keys = list(original_keys- intersect_keys)
        if len(removed_keys)>0: 
            log.removed_keys = json.dumps(removed_keys)
        
        diff_keys = list()
        
        for key in intersect_keys:
            val1 = original_bundle.data[key]
            val2 = updated_bundle.data[key]
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
        #                     if original_bundle.data[key] != updated_bundle.data[key])

        if len(diff_keys)>0: 
            log.diff_keys = json.dumps(diff_keys)
            log.diffs = json.dumps(dict(
                zip(diff_keys, 
                    ([original_bundle.data[key],updated_bundle.data[key]] 
                        for key in diff_keys )  )))
            
        # user can specify any valid, escaped json for this field
        # FIXME: untested
        if 'apilog_json_field' in bundle.data:
#             log.json_field = json.dumps(bundle.data['apilog_json_field'])
            log.json_field = bundle.data['apilog_json_field']
        
        # TODO: how do we feel about passing form data in the headers?
        # TODO: abstract the form field name
        if 'HTTP_APILOG_COMMENT' in bundle.request.META:
            log.comment = bundle.request.META['HTTP_APILOG_COMMENT']
            logger.info(str(('log comment', log.comment)))

        i +=1
        log.save()
        logger.info(str(('update, api log', log)) )
        
        return updated_bundle
                  

# NOTE if using this class, must implement the "not implemented error" methods
# on Resource (Are implemented with ModelResource)
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
    def reset_field_defs(self, scope):
        #         logger.info(str((
        #             '----------reset_field_defs, ' , scope, 'registry', 
        #             ManagedResource.resource_registry )))
        if scope not in ManagedResource.resource_registry:
            msg = str((
                'resource for scope not found: ', scope, 
                'in resource registry',ManagedResource.resource_registry.keys(),
                'possible cause: resource not entered in urls.py' ))
            logger.warn(msg)
            raise Exception(msg)
        resource = ManagedResource.resource_registry[scope]
        logger.info(str((
            '----------reset_field_defs, resource_name' , 
            resource._meta.resource_name, 'scope', scope, 'resource', resource )))
        resource.create_fields();
        resource.reset_filtering_and_ordering();
    
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
    
    # locally defined
    def get_field_def(self, name):
        return self.local_field_defs[name]
    
    # local method
    def create_fields(self):
        
        logger.debug(str(('--create_fields', self._meta.resource_name, 
            self.scope, 'original fields', self.original_fields.keys() )))
        if hasattr(self._meta, 'bootstrap_fields'):
            logger.debug(str(('bootstrap fields', self._meta.bootstrap_fields)))
        
        
        self.local_field_defs = local_field_defs = \
            MetaHash.objects.get_and_parse(scope=self.scope, 
                field_definition_scope='fields.metahash', clear=True)
        logger.debug(str(('managed fields to create', local_field_defs.keys())))
        new_fields = {}
        for field_name, field_obj in self.original_fields.items():
            if field_name in local_field_defs:
                new_fields[field_name] = deepcopy(field_obj)
            elif ( hasattr(self._meta, 'bootstrap_fields') 
                    and field_name in self._meta.bootstrap_fields ):
                logger.debug('====== bootstrapping field: ' + field_name)
                new_fields[field_name] = deepcopy(field_obj)
            elif field_name in self.meta_bootstrap_fields:
                logger.debug('====== meta bootstrapping field: ' + field_name)
                new_fields[field_name] = deepcopy(field_obj)

        unknown_keys = set(local_field_defs.keys()) - set(new_fields.keys())
        logger.debug(str(('managed keys not yet defined', unknown_keys)))
        for field_name in unknown_keys:
            field_def = local_field_defs[field_name]
            logger.debug(str(('virtual managed field:', field_name, field_def)))
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
            else:
                logger.debug('creating unknown field as a char: ' + field_name)
                new_fields[field_name] = fields.CharField(
                    attribute=field_name, readonly=True, blank=True, null=True)
                
        logger.debug(str((
            'resource', self._meta.resource_name, self.scope, 
            'create_fields done: fields created', new_fields.keys() )))
        self.fields = new_fields
        return self.fields

    def build_schema(self):
        '''
        Override
        '''
        logger.debug('------build_schema: ' + self.scope)
        try:
            schema = {}
            schema['fields'] = deepcopy(
                MetaHash.objects.get_and_parse(
                    scope=self.scope, field_definition_scope='fields.metahash'))
            
            if 'json_field' in schema['fields']: 
                # because we don't want this serialized directly (see dehydrate)
                schema['fields'].pop('json_field')  
            if not 'resource_uri' in schema['fields']:
                schema['fields']['resource_uri'] = { 'visibility':[] }
            if not 'id' in schema['fields']:
                schema['fields']['id'] = { 'visibility':[] }
            
            # FIXME: schema.resource_definition <=> resource.schema, which one?
            logger.debug(str((
                'trying to locate resource information', 
                self._meta.resource_name, self.scope)))
            resource_def = MetaHash.objects.get(
                scope='resource', key=self._meta.resource_name)
            schema['resource_definition'] = resource_def.model_to_dict(scope='fields.resource')

        except Exception, e:
            logger.warn(str(('on building schema', e, self._meta.resource_name)))
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            logger.error(str((exc_type, fname, exc_tb.tb_lineno)))
            raise e
            
        logger.debug('------build_schema,done: ' + self.scope)
        return schema
    
    def dehydrate(self, bundle):
        ''' 
        Implementation hook method, override to augment bundle, post dehydrate
        by the superclass used here to do the "hydrate_json_field"
        '''
        if len(bundle.data) == 0 : return bundle
        
        local_field_defs = MetaHash.objects.get_and_parse(
            scope=self.scope, field_definition_scope='fields.metahash')
        for key in [ 
                x for x,y in local_field_defs.items() if y.get('json_field_type') ]:
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
        logger.debug(str(('hydrate_json_field', bundle)))
        
        json_obj = {}
        local_field_defs = MetaHash.objects.get_and_parse(
            scope=self.scope, field_definition_scope='fields.metahash', clear=True)
        logger.debug(str(('local_field_defs',local_field_defs)))
        
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
                    # TODO: my my this is complicated, couldn't we just rethrow?
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
    def obj_create(self, bundle, **kwargs):
        try:
            bundle = super(ManagedResource, self).obj_create(bundle, **kwargs);
            return bundle
        except Exception, e:
            logger.warn(str(('==ex on create, kwargs', kwargs,
                             'request.path', bundle.request.path,e)))
            raise e
#             extype, ex, tb = sys.exc_info()
#             logger.warn(str((
#                 'throw', e, tb.tb_frame.f_code.co_filename, 'error line', 
#                 tb.tb_lineno, extype, ex)))
#             raise type(e), str(( type(e), e, 
#                                  'request.path', bundle.request.path, kwargs))

    # override
    def obj_update(self, bundle, **kwargs):
        try:
            bundle = super(ManagedResource, self).obj_update(bundle, **kwargs);
            return bundle
        except Exception, e:
            response = None;
            if hasattr(e, 'response'): response = e.response
            logger.warn(str(('==ex on update',bundle.request.path,e, response )))
            raise e

    # override
    def obj_get(self, bundle, **kwargs):
        try:
            bundle = super(ManagedResource, self).obj_get(bundle, **kwargs);
            return bundle
        except Exception, e:
            logger.warn(str(('==ex on get, kwargs', kwargs,
                             'request.path', bundle.request.path,e)))
            raise e
#             extype, ex, tb = sys.exc_info()
#             logger.warn(str((
#                 'throw', e, tb.tb_frame.f_code.co_filename, 'error line', 
#                 tb.tb_lineno, extype, ex)))
#             logger.warn(str(('==ex on get, kwargs', kwargs, e)))
#             raise type(e), str((type(e), e,
#                                 'request.path', bundle.request.path, kwargs))

    # override
    def detail_uri_kwargs(self, bundle_or_obj):
        """
        Override resources.ModelResource
        Given a ``Bundle`` or an object (typically a ``Model`` instance),
        it returns the extra kwargs needed to generate a detail URI.

        By default, it uses the model's ``pk`` in order to create the URI.
        """
        
        resource_name = self._meta.resource_name
        try:
            resource_def = MetaHash.objects.get(
                scope='resource', key=resource_name)
            resource = resource_def.model_to_dict(scope='fields.resource')
            
            # TODO: memoize
            # note use an ordered dict here so that the args can be returned as
            # a positional array for 
            kwargs = OrderedDict() 

            for x in resource['id_attribute']:
                val = ''
                if isinstance(bundle_or_obj, Bundle):
                    val = getattr(bundle_or_obj.obj,x)
                else:
                    if hasattr(bundle_or_obj, x):
                        val = getattr(bundle_or_obj,x) # it may be an object- 
                    else:
                        val = bundle_or_obj[x] # allows simple dicts
                kwargs[x] = str(val)
            
            return kwargs
            
        except Exception, e:
            logger.warn(str(('cannot grab id_attribute for', resource_name, e)))
            
            try:
                if logger.isEnabledFor(logging.INFO):
                    logger.info(str((
                        'unable to locate resource information[id_attribute];'
                        ' has it been loaded yet for this resource?',
                        'also note that this may not work with south, since model methods',
                        'are not available: ', resource_name, e, 
                        'type', type(bundle_or_obj), 'attr', resource['id_attribute'],
                        'bundle', bundle_or_obj,
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

    def create_response(self, *args, **kwargs):
        '''
        Override to set a Content-Disposition attachment header in the response;
        - this can be used by the browser to give the download a name.
        - TODO: this is an expedient solution to the download button in the 
        browser, but a better solution will probably be a javascript utility - sde 201404
        '''
        response = super(ManagedResource, self).create_response(*args, **kwargs)
        
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
            logger.info(str(('kwargs', kwargs)))
            _kwargs = kwargs
            if 'request' in _kwargs: 
                _kwargs = {}
            objects = self.apply_filters(bundle.request, applicable_filters, **_kwargs)
            return self.authorized_read_list(objects, bundle)
        except ValueError:
            raise BadRequest("Invalid resource lookup data provided (mismatched type).")
        
        
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
    
    class Meta:
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field_type', 
                            'json_field']
        queryset = MetaHash.objects.filter(
            scope__startswith="fields.").order_by('scope','ordinal','key')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= SuperUserAuthorization()        
        ordering = []
        filtering = {} #{'scope':ALL, 'key':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
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
        authorization= SuperUserAuthorization()        
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name = 'vocabularies'
    
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

    def dehydrate(self, bundle):
        bundle = super(ResourceResource,self).dehydrate(bundle)
        # Get the schema
        # FIXME: why is the resource registry keyed off of "field."+key ?
        resource = ManagedResource.resource_registry['fields.'+bundle.obj.key]
        if resource:
            bundle.data['schema'] = resource.build_schema();
        else:
            logger.error('no API resource found in the registry for ' + 
                         bundle.data['key'] + 
                         '.  Cannot build the schema for this resource.' )
        return bundle


class ApiLogResource(ManagedModelResource):
    
    class Meta:
        queryset = ApiLog.objects.all().order_by(
            'ref_resource_name', 'username','date_time')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        ordering = []
        filtering = {'username':ALL, 'uri': ALL, 'ref_resource_name':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='apilog' 
    
    def __init__(self, **kwargs):
        self.scope = 'fields.apilog'
        super(ApiLogResource,self).__init__(**kwargs)

    def build_schema(self):
        schema = super(ApiLogResource,self).build_schema()
        temp = [ 
            x.key for x in 
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
            url((r"^(?P<resource_name>%s)/(?P<ref_resource_name>[\w\d_.\-:]+)"
                 r"/(?P<key>[\w\d_.\-\+:]+)"
                 r"/(?P<date_time>[\w\d_.\-\+:]+)%s$")
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

# class CustomAuthentication(BasicAuthentication):
#     '''
#     Work-around authentication for dev on orchestra:
#     orchestra apache strips Basic authentication headers, and more work is needed
#     to store the csrf token.
#     '''
#     def is_authenticated(self, request, **kwargs):
#         '''
#         Use simple session authentication
#         NOTE: this does not perform csrf checks
#         '''
#         
#         logger.info(str(('=== in custom authentication', request.user.is_authenticated())))
#         return request.user.is_authenticated()    



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
    is_for_group = fields.BooleanField(attribute='is_for_group', blank=True, null=True)

    def __init__(self, **kwargs):
        super(UserResource,self).__init__(**kwargs)

    class Meta:
        queryset = UserProfile.objects.all().order_by('username') 
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
#         authorization= UserGroupAuthorization() #SuperUserAuthorization()        
        authorization= SuperUserAuthorization()        
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
            logger.info('save django user')
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
        
        logger.info(str(('AUTH_USER_MODEL', settings.AUTH_USER_MODEL)))
        logger.info(str(('get_obj_list', groupname)))
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
    is_for_user = fields.BooleanField(attribute='is_for_user', blank=True, null=True)

    class Meta:
        queryset = UserGroup.objects.all();        
        
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= SuperUserAuthorization()        

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
            ]

    def dispatch_group_userview(self, request, **kwargs):
        # signal to include extra column
        kwargs['groupname'] = kwargs.pop('name')  
        return UserResource().dispatch('list', request, **kwargs)    
    
    def dispatch_group_permissionview(self, request, **kwargs):
        # signal to include extra column
        kwargs['groupname'] = kwargs.pop('name')  
        return PermissionResource().dispatch('list', request, **kwargs)    

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
        logger.info(str(('get_obj_list', kwargs)))
        query = super(UserGroupResource, self).get_object_list(request);
        
        if 'username' in kwargs:
            is_for_user = kwargs.pop('username')
            logger.info(str(('get_obj_list', is_for_user)))
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
        options['non_null_fields'] = ['is_for_user'] # Override to exclude this field in the PostgresSortingResource 
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
        authorization= SuperUserAuthorization()        
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
        


