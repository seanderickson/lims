import json
import re
import sys
import os

from collections import defaultdict

from copy import deepcopy
from django.conf.urls import url
from django.utils.encoding import smart_str, smart_text
from django.utils import timezone
from django.forms.models import model_to_dict
from django.core.exceptions import ObjectDoesNotExist
from tastypie.exceptions import NotFound, ImmediateHttpResponse, TastypieError
from tastypie.bundle import Bundle
import traceback
#from django.contrib.auth.models import Group, User, Permission

from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
from tastypie import fields # NOTE: required for dynamic resource field definitions using eval
from tastypie.resources import Resource

from lims.api import PostgresSortingResource, LimsSerializer
from reports.models import MetaHash, Vocabularies, ApiLog, Permission, UserGroup

import logging
from tastypie.utils.urls import trailing_slash
from django.contrib.auth.models import User
from collections import OrderedDict
from db.models import ScreensaverUser
from django.db.models.aggregates import Max
#from reports import dump_obj
from lims.api import CsvBooleanField
        
logger = logging.getLogger(__name__)


class LoggingMixin(Resource):
    '''
    intercepts obj_create, obj_update and creates an ApiLog entry for the action
    Note: whatever is being extended with the LoggingMixin must also define a
    "detail_uri_kwargs" method that returns an _ordered_dict_, since we log the kwargs as ordered args.
    ** note: "detail_uri_kwargs" returns the set of lookup keys for the resource URI construction.
    '''
    def obj_create(self, bundle, **kwargs):
        logger.info('----log obj_create')
        
        bundle = super(LoggingMixin, self).obj_create(bundle=bundle, **kwargs)
        logger.info(str(('object created', bundle.obj )))
        log = ApiLog()
        log.username = bundle.request.user.username #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.user_id = bundle.request.user.id #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
#        log.diffs = json.dumps(bundle.obj)
            
        # user can specify any valid, escaped json for this field
        if 'apilog_json_field' in bundle.data:
            log.json_field = json.dumps(bundle.data['apilog_json_field'])
#        log.uri,log.key = self.get_uri(bundle.request.get_full_path(), self._meta.resource_name, bundle.obj, bundle.data)
        log.uri = self.get_resource_uri(bundle)
        log.key = '/'.join([str(x) for x in self.detail_uri_kwargs(bundle).values()])

        log.save()
        logger.info(str(('create, api log', log)) )
        
        return bundle    
    
    def obj_delete(self, bundle, **kwargs):
        logger.info('---log obj_delete')
        
        bundle = super(LoggingMixin, self).obj_create(bundle=bundle, **kwargs)
        
        log = ApiLog()
        log.username = bundle.request.user.username #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.user_id = bundle.request.user.id #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
#        log.diffs = json.dumps(bundle.obj)
                    
        # user can specify any valid, escaped json for this field
        if 'apilog_json_field' in bundle.data:
            log.json_field = json.dumps(bundle.data['apilog_json_field'])
#        log.uri,log.key = self.get_uri(bundle.request.get_full_path(), self._meta.resource_name, bundle.obj, bundle.data)
        log.uri = self.get_resource_uri(bundle)
        log.key = '/'.join([str(x) for x in self.detail_uri_kwargs(bundle).values()])

        log.save()
        logger.info(str(('delete, api log', log)) )
        
        return bundle
    
    
    def _locate_obj(self, bundle, **kwargs):
        # lookup the object, the same way that it would be looked up in ModelResource.obj_update
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
                raise NotFound("A model instance matching the provided arguments could not be found.")
        return bundle
    
    
#    def get_uri(self, full_path, resource_name, obj, data):
#        # replace the id with the "id attribute" ids (public ID's)
#        # this is necessary because users of the API may be using internal IDs not specified as the official IDs
#        #
#        # TODO: use Django's "reverse" along with named URLs in urls.py to lookup resource URLS with filled in values
#        # I saw it somewhere on stackoverflow
#        
#        
##        logger.info('--- get_uri: '+ full_path)
#        key = obj.id 
#        
#        # look up the resource definition
#                        
#        try:
##            logger.info(str(('trying to locate resource information', resource_name, self.scope)))
#            resource_def = MetaHash.objects.get(scope='resource', key=resource_name)
##            logger.info(str(('create dict for obj', resource_def)))
#            resource = resource_def.model_to_dict(scope='fields:resource')
##            logger.info(str(('--- got', resource, resource['id_attribute'] )))
#        except Exception, e:
#            logger.warn(str(('unable to locate resource information, has it been loaded yet for this resource?', resource_name, e)))
#            return (full_path, key)
#
#        try:            
#            # parse the key as whatever is after the resource name in the path
#            matchObject = re.match(r'(.*\/'+resource_name+'(\/)?)(.*)', full_path)
#            
#            provided_key = ''
#            if(matchObject):
#                if matchObject.group(2):
#                    provided_key = matchObject.group(2)
#                    if provided_key[len(provided_key)-1] == '/': 
#                        provided_key = provided_key[:len(provided_key)-1]
#            else:
#                raise Exception(str(('resource name not in path', resource_name, full_path )) )
#
#            # If the obj is provided, use it to get the key attributes    
#            if obj and obj.id:
#                new_public_key = "/".join([getattr(obj,x) for x in resource['id_attribute']])
#                logging.info(str(('created new public key', new_public_key, 'for object to log:', obj, 'resource', resource_name)))
#                full_path = matchObject.group(1) + new_public_key #+ '/'
#                key = new_public_key
#
#                if len(provided_key)> 0 and str(obj.id) != provided_key:
#                    logger.info(str(('provided key in the path is not equal to the obj id', full_path, provided_key, obj.id)))     
#            # otherwise, use the data bundle (i.e. not yet persisted data)
#            else:
#                new_public_key = "/".join([data[x] for x in resource['id_attribute']])
#                logging.info(str(('created new public key', new_public_key, 'for data to log:', data, 'resource', resource_name)))
#                full_path = matchObject.group(1) + new_public_key #+ '/'
#                key = new_public_key
#                    
##                    if obj and len(key) > 0:  # if the obj has a key, and you gave me a keyits equal to the key we got in the passed in path, then replace
##                        if (str((obj.id))) == key :
##                            # replace the key given with the key we want to record
##                            new_public_key = "/".join([getattr(obj,x) for x in resource['id_attribute']])
##                            logging.info(str(('created new public key', new_public_key, 'for object to log:', obj, 'resource', resource_name)))
##                            full_path = matchObject.group(1) + new_public_key # + '/'
##                            key = new_public_key
##                        else:
##                            logger.info(str(('id',obj.id,'doesnt equal match obj', key, full_path, resource_name, matchObject.group())))
##                    else:
##                        new_public_key = "/".join([data[x] for x in resource['id_attribute']])
##                        logging.info(str(('created new public key', new_public_key, 'for data to log:', data, 'resource', resource_name)))
##                        full_path = matchObject.group(1) + new_public_key #+ '/'
##                        key = new_public_key
##                        
##                else: # found the resource, but id was not given
##                    if obj:
##                        new_public_key = "/".join([getattr(obj,x) for x in resource['id_attribute']])
##                        logging.info(str(('created new public key', new_public_key, 'for object to log:', obj, 'resource', resource_name)))
##                        full_path = matchObject.group(1) + new_public_key #+ '/'
##                        key = new_public_key
##                    else: # no object, probably means this is a create
##                        new_public_key = "/".join([data[x] for x in resource['id_attribute']])
##                        logging.info(str(('created new public key', new_public_key, 'for data to log:', data, 'resource', resource_name)))
##                        full_path = matchObject.group(1) + new_public_key #+ '/'
##                        key = new_public_key
##                        
##            else:
##                logger.warn(str(('non-standard resource, does not contain resource_name:', resource_name, full_path)))
#        except Exception, e:
#            exc_type, exc_obj, exc_tb = sys.exc_info()
#            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
#            logger.error(str((exc_type, fname, exc_tb.tb_lineno)))
#            logger.error(str(('on trying to extract public key information from the object and data', full_path, resource_name, obj, data, resource, e)))
#            logger.error(str(('using uri,key:', full_path, key)))
#        logger.debug(str(('got uri,key:', full_path, key)))
#        return (full_path, key)
              
    
    def obj_update(self, bundle, skip_errors=False, **kwargs):
        logger.info('--- log obj_update')
        original_bundle = Bundle(data=deepcopy(bundle.data))
        i=0;

        if hasattr(bundle,'obj'): original_bundle.obj = bundle.obj
        original_bundle = self._locate_obj(original_bundle, **kwargs)

        i +=1
        logger.info('--- log obj_update fd' + str(i))
        
        original_bundle = super(LoggingMixin, self).full_dehydrate(original_bundle)
        
        i +=1
        logger.info('--- log obj_update u' + str(i))
        
        updated_bundle = super(LoggingMixin, self).obj_update(bundle=bundle, **kwargs)
        
        i +=1
        logger.info('--- log obj_update fd' + str(i))
        
        updated_bundle = super(LoggingMixin, self).full_dehydrate(updated_bundle)
        
        i +=1
        logger.info('--- log obj_update ' + str(i))
        
        # TODO: diff updated_bundle.data
        original_keys = set(original_bundle.data.keys())
        updated_keys = set(updated_bundle.data.keys())
        
        intersect_keys = original_keys.intersection(updated_keys)
        
        log = ApiLog()
        log.username = bundle.request.user.username #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.user_id = bundle.request.user.id #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.date_time = timezone.now()
        log.ref_resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
#        log.uri,log.key = self.get_uri(bundle.request.get_full_path(), self._meta.resource_name, bundle.obj, bundle.data)
        log.uri = self.get_resource_uri(bundle)
        log.key = '/'.join([str(x) for x in self.detail_uri_kwargs(bundle).values()])
        
        added_keys = list(updated_keys - intersect_keys)
        if len(added_keys)>0: 
            log.added_keys = json.dumps(added_keys)
        
        removed_keys = list(original_keys- intersect_keys)
        if len(removed_keys)>0: 
            log.removed_keys = json.dumps(removed_keys)
        
        diff_keys = list(key for key in intersect_keys if original_bundle.data[key] != updated_bundle.data[key])
        if len(diff_keys)>0: 
            log.diff_keys = json.dumps(diff_keys)
            log.diffs = json.dumps(dict(zip(diff_keys, ([original_bundle.data[key],updated_bundle.data[key]] for key in diff_keys)  )) )
            
        # user can specify any valid, escaped json for this field
        if 'apilog_json_field' in bundle.data:
            log.json_field = json.dumps(bundle.data['apilog_json_field'])
            
        i +=1
        logger.info('--- log obj_update apilog.save ' + str(i))
        log.save()
        logger.info(str(('update, api log', log)) )
        
        return updated_bundle
                  

# NOTE if using this class, must implement the "not implemented error" methods on Resource 
# (Are implemented with ModelResource)
# for instance, "detail_uri_kwargs" which returns the kwargs needed to construct the uri
class ManagedResource(LoggingMixin):
    '''
    Uses the field and resource definitions in the Metahash store to determine the fields to expose for a Resource
    '''
    resource_registry = {}
    
    def __init__(self, field_definition_scope='fields:metahash', **kwargs):
        self.resource = self._meta.resource_name
        self.scope = 'fields:' + self.resource
        self.field_definition_scope = field_definition_scope
        self.meta_bootstrap_fields = ['resource_uri']
        
        logger.debug(str(('---init resource', self.resource, self.scope, field_definition_scope)))
        
        ManagedResource.resource_registry[self.scope] = self;

        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
        metahash = MetaHash.objects.get_and_parse(scope=self.scope, 
                                                  field_definition_scope=field_definition_scope)
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)
        
        super(ManagedResource,self).__init__(**kwargs)
        self.original_fields = deepcopy(self.fields)
        self.create_fields()
        
    # local method  
    def reset_field_defs(self, scope):
        logger.info(str(('----------reset_field_defs, ' , scope, 'registry', ManagedResource.resource_registry )))
        resource = ManagedResource.resource_registry[scope]
        logger.info(str(('----------reset_field_defs, resource_name' , resource._meta.resource_name, 'scope', scope, 'resource', resource )))
        resource.create_fields();
        resource.reset_filtering_and_ordering();

#        self.create_fields()
#        self.reset_filtering_and_ordering()
    
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
        
        logger.debug(str(('--create_fields', self._meta.resource_name, self.scope, 'original fields', self.original_fields.keys() )))
        if hasattr(self._meta, 'bootstrap_fields'):
            logger.debug(str(( 'bootstrap fields', self._meta.bootstrap_fields )))
        
        
        local_field_defs = MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope='fields:metahash', clear=True)
        logger.debug(str(('managed fields to create', local_field_defs.keys())))
        new_fields = {}
        for field_name, field_obj in self.original_fields.items():
            if field_name in local_field_defs:
                new_fields[field_name] = deepcopy(field_obj)
            elif hasattr(self._meta, 'bootstrap_fields') and field_name in self._meta.bootstrap_fields:
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
                # JSON fields are read only because they are hydrated in the hydrate_json_field method
                if field_def['json_field_type'] == 'fields.BooleanField':
                    new_fields[field_name] = eval(field_def['json_field_type'])(attribute=field_name,readonly=True, blank=True, null=True, default=False) 
                else:
                    new_fields[field_name] = eval(field_def['json_field_type'])(attribute=field_name,readonly=True, blank=True, null=True) 
            
            else:
            
                logger.debug('creating unknown field as a char: ' + field_name)
                # TODO: better
                new_fields[field_name] = fields.CharField(attribute=field_name, readonly=True, blank=True, null=True)
#                    msg = 'no field defined for ' + name + ', at this time fields must be either json_fields, or defined normally as class fields'
#                    raise Exception(msg)
        
        
        logger.debug(str(('resource', self._meta.resource_name, self.scope, 'create_fields done: fields created', new_fields.keys() )))
        self.fields = new_fields
        return self.fields

    # override ModelResource      
    def build_schema(self):
        logger.debug('------build_schema: ' + self.scope)
        try:
            schema = {}
            schema['fields'] = deepcopy(MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope='fields:metahash'))
            
            if 'json_field' in schema['fields']: 
                schema['fields'].pop('json_field')  # because we don't want this serialized directly (see dehydrate)
            
            if not 'resource_uri' in schema['fields']:
                schema['fields']['resource_uri'] = { 'visibility':[] }
            if not 'id' in schema['fields']:
                schema['fields']['id'] = { 'visibility':[] }
            
            logger.debug(str(('trying to locate resource information', self._meta.resource_name, self.scope)))
            resource_def = MetaHash.objects.get(scope='resource', key=self._meta.resource_name)
        
            schema['resource_definition'] = resource_def.model_to_dict(scope='fields:resource')
        except Exception, e:
            logger.warn(str(('on building schema', e, self._meta.resource_name)))
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            logger.error(str((exc_type, fname, exc_tb.tb_lineno)))
            raise e
            
        logger.debug('------build_schema,done: ' + self.scope)
        return schema
    
    # implementation hook method, override to augment bundle, post dehydrate by the superclass
    # used here to do the "hydrate_json_field"
    def dehydrate(self, bundle):
        logger.debug(str(('------dehydrate: ', self.scope)) )

        if len(bundle.data) == 0 : return bundle
        
        local_field_defs = MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope='fields:metahash')
        for key in [ x for x,y in local_field_defs.items() if y.get('json_field_type') ]:
            bundle.data[key] = bundle.obj.get_field(key);
#            logger.info(str((key, bundle.data[key])))
        
        bundle.data['json_field'] = ''
        bundle.data.pop('json_field') # json_field will not be part of the public API, it is for internal use only
#        logger.info(str(('bundle.data', bundle.data)))
        # override the resource_uri, since we want to export the permanent composite key
#        bundle.data['resource_uri'] = self.build_resource_uri(self.resource, bundle.data) or bundle.data['resource_uri']
        
        return bundle
    
    # implementation hook, to deserialize the embedded json fields
    def hydrate_json_field(self, bundle):
        '''
        hydrate bundle data values that will be stuffed into the json_field
        -Note: as mentioned elsewhere, for the initial load of the Metahash:fields, fields that
        are JSON fields (to be stuffed into json_field) must be first defined as a scope='metahash:field',
        then they can be set as attribute values on other fields in the second step.
        '''
#        bundle = super(JsonAndDatabaseResource, self).hydrate(bundle);
        
        logger.debug(str(('hydrate_json_field', bundle)))
        
        json_obj = {}
        local_field_defs = MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope='fields:metahash', clear=True)
        logger.debug(str(('local_field_defs',local_field_defs)))
        
        # Use the tastypie field type that has been designated to convert each field in the json stuffed field just like it were a real db field
        for key in [ str(x) for x,y in local_field_defs.items() if 'json_field_type' in y and y['json_field_type'] ]:
            if key not in self.fields:
                raise RuntimeError(str(('for the resource', self._meta.resource_name, 
                    'the key to deserialize', key, 'was not defined as a resource field: fields:', self.fields.keys() )))
            val = bundle.data.get(key,None)
            if val:
                try:
#                    logger.info(str(('get value for key: ', key, val )))
                    if hasattr(val, "strip"): # test if it is a string
                        val = self.fields[key].convert(smart_text(val,'utf-8', errors='ignore'))
#                        logger.info(str(('got value for key: ', key, val)))
                    elif hasattr(val, "__getitem__") or hasattr(val, "__iter__"): # test if it is a sequence
                        val = [smart_text(x,'utf-8', errors='ignore') for x in val]
#                    logger.info(str(('got',key, val)))
                    json_obj.update({ key:val })
                except Exception, e:
                    # TODO: my my this is complicated, couldn't we just rethrow?
                    logger.error('ex', e)
                    extype, ex, tb = sys.exc_info()
                    formatted = traceback.format_exception_only(extype, ex)[-1]
                    msg = str(('failed to convert', key, 'with value', val, 'message', formatted)).replace("'","")
                    if key in self.fields:
                        msg += str(('with tastypie field type', type(self.fields[key]) ))
                    e =  RuntimeError, msg
                    logger.warn(str(('throw', e, tb.tb_frame.f_code.co_filename, 'error line', tb.tb_lineno)))
                    raise e
        bundle.data['json_field'] = json.dumps(json_obj);
        logger.debug(str(('--- hydrated:', bundle.data['json_field'])))
        return bundle;

    
    # override
    def obj_create(self, bundle, **kwargs):
        bundle = super(ManagedResource, self).obj_create(bundle, **kwargs);
        return bundle

    # override
    def obj_update(self, bundle, **kwargs):
        bundle = super(ManagedResource, self).obj_update(bundle, **kwargs);
        return bundle

#    # local method
#    def build_key(self, resource_name, data):
#        try:
#            resource_def = MetaHash.objects.get(scope='resource', key=resource_name)
#            resource = resource_def.model_to_dict(scope='fields:resource')
##            logger.info(str(('found resource', resource)))
#            
#            key_bits = []
#            for x in resource['id_attribute']:
#                if hasattr(data,x ):
#                    key_bits.append(getattr(data,x))
#                else:
#                    key_bits.append(data[x])
#            return  "/".join([str(x) for x in key_bits])
#        except Exception, e:
#            logger.warn(str(('unable to locate resource information[id_attribute]; has it been loaded yet for this resource?', resource_name, e)))
#    
#    URL_BASE = '/reports/api'
#
#    # local method
#    # TODO: REDO this using resources.Resource.get_resource_uri and resource_uri_kwargs methods
#    def build_resource_uri(self, resource_name, data):
#        new_key = self.build_key(resource_name, data)
#        if new_key:
##            base_uri = self.get_resource_uri()
##            if base_uri[len(base_uri)-1] != '/':
##                base_uri += '/'
##            return base_uri + new_key
#            return self.URL_BASE + '/' + self._meta.api_name + '/'  + resource_name + '/' + new_key
#
#    # local method
#    def obj_resource_uri(self, resource_name, obj):
#        new_key = self.build_key(resource_name, obj)
#        if new_key:
#            return self.URL_BASE + '/' + self._meta.api_name + '/'  + resource_name + '/' + new_key
##            base_uri = self.get_resource_uri()
##            if base_uri[len(base_uri)-1] != '/':
##                base_uri += '/'
##            return base_uri + new_key

    def detail_uri_kwargs(self, bundle_or_obj):
        resource_name = self._meta.resource_name
        try:
            resource_def = MetaHash.objects.get(scope='resource', key=resource_name)
            resource = resource_def.model_to_dict(scope='fields:resource')
#            logger.info(str(('found resource', resource)))
            
            # TODO: memoize
            kwargs = OrderedDict() # note use an ordered dict here so that the args can be returned as a positional array for 
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
            logger.warn(str(('unable to locate resource information[id_attribute]; has it been loaded yet for this resource?', resource_name, e)))
        # Fall back to base class implementation (using the declared primary key only, for ModelResource)
        # This is useful in order to bootstrap the ResourceResource
        return super(ManagedResource,self).detail_uri_kwargs(bundle_or_obj)

    def get_via_uri(self, uri, request=None):
        '''
        Override the stock method to allow lookup of relative uri's:
        - a 'relative uri' - or 'local uri' is one that doesn't include the api name ("v1" for instance), 
        but rather, is truncated on the left, so that "api/vi/resource_name/key1" becomes
        "resource_name/key1".  This is useful because input file records can have a
        shorter "resource_uri" field.
        '''
        if self._meta.resource_name not in uri:
            raise Exception(str(('invalid URI', uri, 'must contain at least the resource name', self._meta.resource_name)))
        
        if request and request.path:
            path = request.path
            # remove the parts after the api_name ("v1") because that part is the resource name, 
            # calling context may not be for this resource
            path = path[:path.find(self._meta.api_name)+len(self._meta.api_name)+1] 
            local_uri = uri
            if path not in local_uri:
                uri = path + local_uri
                #            uri = '/reports/api/v1/'+ uri  # TODO: poc, not permanent
                logger.info(str(('converted local URI', local_uri, ' to tastypie URI', uri)))
        
        return super(ManagedResource, self).get_via_uri(uri, request);

    def get_local_resource_uri(self, bundle_or_obj=None, url_name='api_dispatch_list'):
        '''
        special 'local' version of the uri - 
        when creating the uri for containment lists (user.permissionss, for example),
        convert "reports/api/v1/permission/resource/read" to "permission/resource/read"
        '''
        uri = super(ManagedResource, self).get_resource_uri(bundle_or_obj=bundle_or_obj, url_name=url_name)
#        return uri
        return uri[uri.find(self._meta.resource_name):]
    
    # implementation hook - URLS to match _before_ the default URLS
    # used here to allow the natural keys [scope, key] to be used
    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" 
        # [ any word, except "schema" ]
        # also note the double underscore "__" is because we also don't want to match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.-:]+)/(?P<key>[\w\d_.-]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            # TODO: is this needed here on metahash? we aren't using just "key" as a key, which is what causes the conflict with "schema", so probably not
            url(r"^(?P<resource_name>%s)/(?P<key>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
 
 
class ManagedModelResource(ManagedResource, PostgresSortingResource):
    pass

class MetaHashResource(ManagedModelResource):
    
    class Meta:
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field_type', 'json_field']
        queryset = MetaHash.objects.filter(scope__startswith="fields:").order_by('scope','ordinal','key')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        ordering = []
        filtering = {} #{'scope':ALL, 'key':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name = 'metahash'

    def __init__(self, **kwargs):
        super(MetaHashResource,self).__init__(**kwargs)
    
    # because the metahash resource is both a resource and the definer of json fields, 
    # reset_field_defs after each create/update, in case, new json fields are defined,
    # or in case ordering,filtering groups are updated
    def obj_create(self, bundle, **kwargs):
        bundle = super(MetaHashResource, self).obj_create(bundle, **kwargs);
        if getattr(bundle.obj,'scope').find('fields') == 0: #'fields:metahash':
            self.reset_field_defs(getattr(bundle.obj,'scope'))
#        if getattr(bundle.obj,'scope') == self.scope: #'fields:metahash':
#            logger.info(str(('post-obj_create, bundle.obj', bundle.obj, 'reset the hash')))
#            self.reset_field_defs();
#        elif getattr(bundle.obj,'scope') == 'fields:resource': #'fields:metahash':
#            logger.info(str(('post-obj_create, bundle.obj', bundle.obj, 'reset the hash')))
#            ResourceResource().reset_field_defs();
        return bundle

    def obj_update(self, bundle, **kwargs):
        bundle = super(MetaHashResource, self).obj_update(bundle, **kwargs);
        self.reset_field_defs(getattr(bundle.obj,'scope'))
#        if getattr(bundle.obj,'scope') == self.scope: #'fields:metahash':
#            if getattr(bundle.obj,'json_field_type'):
#                logger.info(str(('post-obj_update, bundle.obj', bundle.obj, 'reset the hash')))
#                self.reset_field_defs();
        return bundle

    def hydrate(self, bundle):
        bundle = super(MetaHashResource, self).hydrate(bundle);
        return bundle
    
    def build_schema(self):
        schema = super(MetaHashResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
        return schema
        
    # override, because the metahash resource is special, and will always use a /scope/key/ key    
    def build_key(self, resource_name, data):
        return data['scope'] + '/' + data['key']

class VocabulariesResource(ManagedModelResource):
    '''
    This resource extends the ManagedModelResource using a new table (vocabularies) but has
    fields defined in the Metahash table.
    '''
    def __init__(self, **kwargs):
        super(VocabulariesResource,self).__init__(**kwargs)

    class Meta:
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field']
        queryset = Vocabularies.objects.all().order_by('scope', 'ordinal', 'key')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name = 'vocabularies'
    
#    # called by ModelResource in ModelResource.apply_filters 
#    def get_object_list(self, request):
#        obj_list = super(VocabulariesResource, self).get_object_list(request);
#        # could apply filtering here
#        return obj_list
#    
#    def obj_get(self, bundle, **kwargs):
#        obj = super(VocabulariesResource, self).obj_get(bundle, **kwargs)
#        return obj
    
    def build_schema(self):
        schema = super(VocabulariesResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 'label': 'Vocabulary', 'searchColumn': 'scope', 'options': temp }
        return schema

class ResourceResource(ManagedModelResource):
    '''
    This resource extends the ManagedModelResource, uses the metahash table internally, and has
    fields defined in the Metahash table.
    '''
    def __init__(self, **kwargs):
        super(ResourceResource,self).__init__(field_definition_scope='fields:resource', **kwargs)

    class Meta:
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field'] # note, does not need the 'json_field_type' since MetahashResource is managing the fields
        queryset = MetaHash.objects.filter(scope='resource').order_by('key', 'ordinal', 'scope')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='resource' # appears that tastypie needs this with a resource named "resource"!
    
    def build_schema(self):
        schema = super(ResourceResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('key')]
        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'key', 'options': temp }
        return schema

class ApiLogResource(ManagedModelResource):
    '''
    '''
#    date_time = fields.DateTimeField()
    
    class Meta:
        queryset = ApiLog.objects.all().order_by('ref_resource_name', 'username','date_time')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'username':ALL, 'uri': ALL, 'ref_resource_name':ALL}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='apilog' 
    
    def __init__(self, **kwargs):
        self.scope = 'fields:apilog'
        super(ApiLogResource,self).__init__(**kwargs)

    def build_schema(self):
        schema = super(ApiLogResource,self).build_schema()
        temp = [ x.ref_resource_name for x in self.Meta.queryset.distinct('ref_resource_name')]
        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'ref_resource_name', 'options': temp }
        return schema

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<ref_resource_name>[\w\d_.\-:]+)/(?P<date_time>[\w\d_.\-\+:]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

# Resource to access the django auth_user model.
# Utilizes the local user model in db.screensaveruser as well.
# Note: screensaverusers are in the db model for historical reasons (chose not to migrate at genesis)
class UserResource(ManagedModelResource):
    
    # force the pk to be read only so that it doesn't try to create
    # TODO: store readonly attribute in the hash
    screensaver_user_id = fields.IntegerField('screensaver_user_id', readonly=True) 
    usergroups = fields.ToManyField('reports.api.UserGroupResource', 'usergroup_set', related_name='users', blank=True, null=True)
#    django_user = fields.OneToOneField('db.api.ScreensaverUserResource', 'screensaveruser', blank=True, null=True )
    permissions = fields.ToManyField('reports.api.PermissionResource', 'permissions', related_name='users', null=True) #, related_name='users', blank=True, null=True)
    
    class Meta:
        scope = 'fields:user'
        queryset = ScreensaverUser.objects.all()
        
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        resource_name='user' 
        authorization= Authorization()        
        
        excludes = ['password']
        ordering = []
        filtering = { }
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
        super(UserResource,self).__init__(**kwargs)

    # TODO: either have to generate localized resource uri, or, in client, equate localized uri with non-localized uri's 
    # (* see user.js, where _.without() is used) 
    # this modification represents the first choice
    def get_resource_uri(self, bundle_or_obj=None, url_name='api_dispatch_list'):
        return self.get_local_resource_uri(bundle_or_obj=bundle_or_obj, url_name=url_name)

    def hydrate(self, bundle):
        bundle = super(UserResource, self).hydrate(bundle);
        return bundle

    def dehydrate_permissions(self, bundle):
        uri_list = []
        P = PermissionResource()
        for p in bundle.obj.permissions.all():
            uri_list.append(P.get_local_resource_uri(p))
        logger.info(str(('-- deydrate_permissions', uri_list)))
        return uri_list;
    
    def obj_create(self, bundle, **kwargs):
        """
        A iccbl user specific implementation of ``obj_create``.
        Override so that we can call our own save?
        """
        logger.info(str(('++obj_create', bundle.data)))
        bundle.obj = self._meta.object_class()

        for key, value in kwargs.items():
            setattr(bundle.obj, key, value)

        self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        bundle = self.full_hydrate(bundle)
        return self.save(bundle)

    def is_valid(self, bundle):
        """
        Should return a dictionary of error messages. If the dictionary has
        zero items, the data is considered valid. If there are errors, keys
        in the dictionary should be field names and the values should be a list
        of errors, even if there is only one.
        """
        
        # cribbed from tastypie.validation.py - mesh data and obj values, then validate
        data = {}
        if bundle.obj.pk:
            data = model_to_dict(bundle.obj)
        if data is None:
            data = {}
        data.update(bundle.data)
        
        # do validations
        errors = defaultdict(list)
        
        # TODO: rework this to be driven by the metahash
        
        # TODO: clean up model, use only login_id
        if data.get('login_id') and data.get('ecommons_id'):
            errors['login_id'] = ['specify either ecommons or login_id']
        if not ( data.get('login_id') or data.get('ecommons_id') ):
            errors['login_id'] = ['specify either ecommons or login_id']
        
        if not data.get('first_name'):
            errors['first_name'] = ['first_name must be specified']
        
        if not data.get('last_name'):
            errors['last_name'] = ['last_name must be specified']
        
        if not data.get('email'):
            errors['email'] = ['email must be specified']
        
        username = data.get('ecommons_id')
        if not username:
            username = data.get('login_id')

        # TODO: verify that we should not check for extant here
#        # Special validations, because there are two user objects, the ScreensaverUser, django auth.User
#        try:
#            extant_user = User.objects.get(username=username)
#            errors['login_id'].append('login_id is already in use: ' + username)
#            errors['ecommons_id'].append('ecommons_id is already in use: ' + username)
#            
#            # TODO: email is not unique, should it be?
#            #            extant_user = User.objects.get(email=data['email'])
#            #            errors['email'].append('email is already in use: ' + username)
#        except ObjectDoesNotExist:
#            pass
        
        if errors:
            bundle.errors[self._meta.resource_name] = errors
            logger.warn(str(('bundle errors', bundle.errors, len(bundle.errors.keys()))))
            return False
        return True
        
    
    def save(self, bundle, skip_errors=False):
        ''' 
        overriding base save - so that we can create the django auth_user if needed.
        - everything else should be the same (todo: update if not)
        TODO: document: login_id overrides ecommons_id for assignment to the auth_user.login_id
        when creating the auth_user object
        '''
        logger.info(str(('+save', bundle.obj, bundle.obj.user)))
        self.is_valid(bundle)

        if bundle.errors and not skip_errors:
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))

        # Check if they're authorized.
        if bundle.obj.pk:
            self.authorized_update_detail(self.get_object_list(bundle.request), bundle)
        else:
            self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
        
        # create a new screensaver_user_id
        if not bundle.obj.screensaver_user_id:
            max_result = ScreensaverUser.objects.all().aggregate(Max('screensaver_user_id'))
            new_id = max_result.get('screensaver_user_id__max', 0) or 0
            new_id +=1
            logger.info(str(('creating new screensaver_user_id',new_id)))
            bundle.obj.screensaver_user_id = new_id
            
            bundle.obj.date_created = timezone.now()
        
        # Save FKs just in case.
        self.save_related(bundle)

        logger.info(str(('saving', bundle.obj, bundle.obj.user)))
        bundle.obj.save();
        
        # create a Django user
        username = bundle.obj.ecommons_id
        if bundle.obj.login_id:  # TODO: document: login_id overrides ecommons_id
            username = bundle.obj.login_id
        
        if not bundle.obj.user:
            logger.info('==========create a django user for username: ' + username )
            django_user = User.objects.create_user(username, 
                email=bundle.obj.email, 
                first_name=bundle.obj.first_name, 
                last_name=bundle.obj.last_name)
            django_user.screensaveruser = bundle.obj
            logger.info('save django user')
            django_user.save()
            # this has to be done to set the FK on obj; since we're the only side maintaining this rel' with auth_user
            bundle.obj.user=django_user 
            bundle.obj.save()
            
        bundle.objects_saved.add(self.create_identifier(bundle.obj))

        # Now pick up the M2M bits.
        m2m_bundle = self.hydrate_m2m(bundle)
        self.save_m2m(m2m_bundle)

        #TODO: set password as a separate step
        logger.info('user save done')
        return bundle


    def obj_update(self, bundle, skip_errors=False, **kwargs):
        """
        A iccbl user-specific implementation of ``obj_update``. NOTE: this is the same as superclass right now
        """
        logger.info(str(('++obj_update', bundle.data)))
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
                raise NotFound("A model instance matching the provided arguments could not be found.")

        bundle = self.full_hydrate(bundle)
        return self.save(bundle, skip_errors=skip_errors)
    
    # not tested
    def obj_delete(self, bundle, **kwargs):
        
        django_user = bundle.obj.user
        
        super(UserResource, self).obj_delete(bundle,**kwargs)
        
        if django_user:
            django_user.delete()
        
    def dehydrate_group_list(self,bundle):
        return [x.name for x in bundle.obj.usergroup_set.all()]

# TODO: require a custom list-o-objects viewer 
#    def dehydrate_permission_list(self, bundle):
#        return  [ [x.scope, x.key, x.type] for x in bundle.obj.permissions.all()]

    def prepend_urls(self):
        return [
#            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]    

# see http://www.maykinmedia.nl/blog/2012/oct/2/nested-resources-tastypie/
#    def dispatch_child(self, request, **kwargs):
#        username = kwargs.pop('username')
#        kwargs['users__username'] = username;
#        logger.info(str(('kwargs', kwargs)))
#        return PermissionResource().dispatch('list', request, **kwargs)    
#    


class UserGroupResource(ManagedModelResource):
    
    permissions = fields.ToManyField('reports.api.PermissionResource', 'permissions',related_name='groups', null=True) #, related_name='users', blank=True, null=True)
 
    # relational fields must be defined   
    users = fields.ToManyField('reports.api.UserResource', 'users', related_name='groups', blank=True, null=True)
    
    class Meta:
        queryset = UserGroup.objects.all();
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        

        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='usergroup' 
    
    def __init__(self, **kwargs):
        super(UserGroupResource,self).__init__(**kwargs)

    def get_resource_uri(self, bundle_or_obj=None, url_name='api_dispatch_list'):
        return self.get_local_resource_uri(bundle_or_obj=bundle_or_obj, url_name=url_name)
    
    def obj_create(self, bundle, **kwargs):
        logger.info(str(('----log obj_create', bundle.data)))
        bundle = super(UserGroupResource, self).obj_create(bundle=bundle, **kwargs)
        logger.info(str(('----log obj_created', model_to_dict(bundle.obj))))

    def build_schema(self):
        schema = super(UserGroupResource,self).build_schema()
#        temp = [ x.ref_resource_name for x in self.Meta.queryset.distinct('ref_resource_name')]
#        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'ref_resource_name', 'options': temp }
        return schema
    
    def dehydrate_permission_list(self, bundle):
        logger.info('dehydrate permissions')
        try:
            screensaverUser = bundle.obj.screensaveruser
            if screensaverUser:
#                permissions = [ model_to_dict(x, exclude='id') for x in screensaverUser.permissions.all()]
                permissions = [ [x.scope, x.key, x.type] for x in screensaverUser.permissions.all()]
                return permissions
        except Exception, e:
            logger.warn(str(('error accessing the screensaver user element', e)))
        
    def dehydrate_user_list(self,bundle):
        users = []
        for user in bundle.obj.users.all():
             users.append( '[ %d - %s %s ]' % (user.screensaver_user_id, user.first_name, user.last_name))
        return users
    
    def dehydrate_users(self, bundle):
        uri_list = []
        U = UserResource()
        for user in bundle.obj.users.all():
             uri_list.append(U.get_local_resource_uri({ 'screensaver_user_id':user.screensaver_user_id }))
        return uri_list;
        
    def dehydrate_permissions(self, bundle):
        uri_list = []
        P = PermissionResource()
        for p in bundle.obj.permissions.all():
             uri_list.append(P.get_local_resource_uri(p))
        return uri_list;
        
    def dehydrate(self,bundle):
        bundle.data['id'] = bundle.obj.id
        return bundle

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<name>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]

from django.db.models import Q
class PermissionResource(ManagedModelResource):
    
    users = fields.CharField()
    groups = fields.CharField()
    
    class Meta:
        # note: the queryset for this resource is actually the permissions
        queryset = Permission.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        object_class = object
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        includes = [] # note, use this so that the queryset fields are not all added by default
        always_return_data = True # this makes Backbone happy
        resource_name='permission' 
    
    def __init__(self, **kwargs):
        super(PermissionResource,self).__init__(**kwargs)
        
        # create all of the permissions on startup
        
        resources = MetaHash.objects.filter(Q(scope='resource')|Q(scope__contains='fields:')).order_by('key', 'ordinal', 'scope')
        query = self._meta.queryset._clone()
        permissionTypes = Vocabularies.objects.all().filter(scope='permission:type')
#        if(len(permissionTypes)==0):
#            raise TastypieError(str(('permission types have not been loaded (are vocabularies loaded into the database?)')))
        for r in resources:
            found = False
            for perm in query:
                if perm.scope==r.scope and perm.key==r.key:
                    found = True
            if not found:
                logger.info(str(('initialize permission not found: ', r.scope, r.key)))
                for ptype in permissionTypes:
                    p = Permission.objects.create(scope=r.scope, key=r.key, type=ptype.key)
                    logger.info(str(('bootstrap created permission', p)))

    def get_resource_uri(self, bundle_or_obj=None, url_name='api_dispatch_list'):
        return self.get_local_resource_uri(bundle_or_obj=bundle_or_obj, url_name=url_name)
    
    def obj_get(self, bundle, **kwargs):
        ''' basically, if a permission is requested that does not exist, it is created
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
        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
        return schema
        
    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.\-:]+)/(?P<key>[\w\d_.\-\+:]+)/(?P<type>[\w\d_.\-\+:]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]
        


## graveyard

# create-on-the-fly concept for permissions.
# (discarded in favor of creating all permissions ahead of time        
#    # Unfortunately, we cannot custom manage our m2m reln' using the provided "hydrate_foo" pattern;
#    # this is because tastypie is calling hydrate_permissions before saving the main object - 
#    # see: tastypie.resources.ModelResource.obj_update and obj_create.
#    # With this methodology, the actual hydration of the m2m uri's occurs _after_
#    # the main save() located at the end of ModelResource.obj_create; 
#    # therefore, we'll have to do the same thing.
#    # So, not this: def hydrate_permissions(self, bundle):
#    # but rather, this:
#    def obj_create(self, bundle, **kwargs):
#        bundle = super(UserGroupResource, self).obj_create(bundle, **kwargs)
#        
#        logger.info(str(('hydrate_permissions')))
#        if 'permissions' in bundle.data:
#            try:
#                permission_array = bundle.data['permissions']
#                for p in permission_array:
#                    try:
#                        permission = Permission.objects.get(**p)
#                    except ObjectDoesNotExist, e:
#                        # create the permission on the fly, if the resource it refers to exists
#                        # todo: make sure that the Resource can also be found
#                        resource = MetaHash.objects.get(scope=p['scope'], key=p['key'])
#                        permission = Permission(**p)
#                        permission.save()
#                    bundle.obj.permissions.add(permission)
#                    logger.info(str(('permission added to group', bundle.obj, permission)))
#            except Exception, e:
#                logger.warn(str(('error accessing the permission attribute on hydrate', e)))
#        return bundle
#    def _fill_in_missing(self):
#        queryset = MetaHash.objects.all().filter(Q(scope='resource')|Q(scope__contains='fields:'))
#        

#    def get_object_list(self, request):
#        query = self._meta.queryset._clone()
#        objs = []
#        i=0
#        for r in query:
#            perms = Permission.objects.filter(scope=r.scope, key=r.key )
#            if len(perms) < 1:
#                objs.append({ 'ordinal': i, 'users':[], 'groups':[], 'scope':r.scope, 'key':r.key, 'type': 'read' })
#                objs.append({  'ordinal': i, 'users':[], 'groups':[], 'scope':r.scope, 'key':r.key, 'type': 'write' })
#            else:
#                for p in perms:
#                    objs.append(model_to_dict(p))
#        logger.info(str(('objs', objs)))
#        return objs

#    def obj_get_list(self, request=None, **kwargs):
#        # Filtering disabled for brevity...
#        return self.get_object_list(request)
#    
#    def obj_get(self, request=None, **kwargs):
#        return Permission.objects.get(**kwargs)
#    
#    def apply_sorting(self, obj_list, options):
#        return obj_list
#    
#    def dehydrate(self,bundle):
#        bundle.data =s bundle.obj
#        return bundle


        
#    # permissions is a nested property - not exposing it as a resource - 
#    # so hydrate/dehydrate will be done custom
#    def dehydrate_permissions(self, bundle):
#        logger.info('dehydrate permissions')
#        try:
#            screensaverUser = bundle.obj.screensaveruser
#            if screensaverUser:
#                permissions = [ model_to_dict(x, exclude='id') for x in screensaverUser.permissions.all()]
#                logger.info(str(('created permissions', permissions)))
#                return permissions
#        except Exception, e:
#            logger.warn(str(('error accessing the screensaver user element', e)))
#            
#    
#    def hydrate_permissions(self, bundle):
#        logger.info(str(('hydrate_permissions')))
#        try:
#            screensaverUser = bundle.obj.screensaveruser
#            
#            permission_array = bundle.data['permissions']
#            for p in permission_array:
#                try:
#                    permission = Permission.objects.get(**p)
#                except ObjectDoesNotExist, e:
#                    # create the permission on the fly, if the resource it refers to exists
#                    # todo: make sure that the Resource can also be found
#                    resource = MetaHash.objects.get(scope=p['scope'], key=p['key'])
#                    permission = Permission(**p)
#                    permission.save()
#                screensaverUser.permissions.add(permission)
#                screensaverUser.save();
#                logger.info(str(('permission saved', permission)))
#        except Exception, e:
#            logger.warn(str(('error accessing the screensaver user element on hydrate', e)))
#        return bundle

#from django.db.models import Q
#class PermissionResource1(ManagedResource):
#    
#    users = fields.CharField()
#    groups = fields.CharField()
#    type = fields.CharField()
#    
#    class Meta:
#        # note: the queryset for this resource is actually the permissions
#        queryset = MetaHash.objects.all().filter(Q(scope='resource')|Q(scope__contains='fields:'))
#        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
#        authorization= Authorization()        
#        object_class = object
#        
#        ordering = []
#        filtering = {}
#        serializer = LimsSerializer()
#        excludes = [] #['json_field']
#        includes = [] # note, use this so that the queryset fields are not all added by default
#        always_return_data = True # this makes Backbone happy
#        resource_name='permission' 
#    
#    def __init__(self, **kwargs):
#        super(PermissionResource1,self).__init__(**kwargs)
#
#
#    def get_object_list(self, request):
#        query = self._meta.queryset._clone()
#        objs = []
#        i=0
#        for r in query:
#            perms = Permission.objects.filter(scope=r.scope, key=r.key )
#            if len(perms) < 1:
#                objs.append({ 'ordinal': i, 'users':[], 'groups':[], 'scope':r.scope, 'key':r.key, 'type': 'read' })
#                objs.append({  'ordinal': i, 'users':[], 'groups':[], 'scope':r.scope, 'key':r.key, 'type': 'write' })
#            else:
#                for p in perms:
#                    objs.append(model_to_dict(p))
#        logger.info(str(('objs', objs)))
#        return objs
#
#    def obj_get_list(self, request=None, **kwargs):
#        # Filtering disabled for brevity...
#        return self.get_object_list(request)
#    
#    def obj_get(self, request=None, **kwargs):
#        
#        return Permission.objects.get(**kwargs)
#    
#    def apply_sorting(self, obj_list, options):
#        return obj_list
#    
#    def dehydrate(self,bundle):
#        bundle.data = bundle.obj
#        return bundle
#
#    def prepend_urls(self):
#        return [
#            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            url(r"^(?P<resource_name>%s)/(?P<scope>((?=(schema))__|(?!(schema))[\w\d_.-]+))/(?P<key>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            ]
#
#class PermissionResource(MetahashManagedResource, PostgresSortingResource):
#
#    groups = fields.ToManyField('reports.api.UserGroupResource', 'group', related_name='permissions', blank=True, null=True)
#    users = fields.ToManyField('reports.api.UserResource', 'user', related_name='permissions', blank=True, null=True)
##    users_set = fields.ToManyField('reports.api.UserResource', 'user_set', related_name='permissions', blank=True, null=True)
#    
#    class Meta:
#        queryset = Permission.objects.all();
#        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
#        authorization= Authorization()        
#        # TODO: drive this from data
#        ordering = []
#        filtering = {}
#        serializer = LimsSerializer()
#        excludes = [] #['json_field']
#        always_return_data = True # this makes Backbone happy
#        resource_name='permission' 
#    
#    def __init__(self, **kwargs):
#        self.scope = 'fields:permission'
#        super(PermissionResource,self).__init__(**kwargs)
#
#    def dehydrate(self, bundle):
#        
#        # final dehydrate hook
#        group_query = bundle.obj.group_set;
#        logger.info(str(('groups', [str(x) for x in group_query.all()])))
#        # note there is an inconsistency in how the reverse m2m reln is traversed
#        # when there is no explicit foreign key attribute on this side of the reln
#        bundle.data['groups'] = [str(x) for x in bundle.obj.group_set.all()]  
#        
#        user_query = bundle.obj.user_set;
#        logger.info(str(('users', [str(x.username) for x in user_query.all()])))
#        bundle.data['users'] = [str(x.username) for x in bundle.obj.user_set.all()]
#        
#        return PostgresSortingResource.dehydrate(self, bundle)
#
#    def build_schema(self):
#        schema = super(PermissionResource,self).build_schema()
##        temp = [ x.ref_resource_name for x in self.Meta.queryset.distinct('ref_resource_name')]
##        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'ref_resource_name', 'options': temp }
#        return schema
#
#    def prepend_urls(self):
#        return [
#            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            url(r"^(?P<resource_name>%s)/(?P<name>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            ]



# ======= Graveyard
  
## TODO: obsolete
## Mixin class - note this class must be mixed in first, since the tastypie Resource class does not call super().__init__
## TODO: merge/refactor this with reports.JsonAndDatabaseResource
#class MetahashManagedResource(object):
#    ''' obsolete '''
#    def __init__(self, **kwargs):
##        if not self.scope:
##            self.scope = kwargs.pop('scope')
#        logger.info(str(('---------init MetahashManagedResource', self.scope)))
#        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
#        metahash = MetaHash.objects.get_and_parse(scope=self.scope)
#        
#        for key,hash in metahash.items():
#            if 'filtering' in hash and hash['filtering']:
#                self.Meta.filtering[key] = ALL_WITH_RELATIONS
#        logger.info(str(('filtering for ', self.scope, self.Meta.filtering )))
#        for key,hash in metahash.items():
#            if 'ordering' in hash and hash['ordering'] and key not in self.Meta.ordering:
#                self.Meta.ordering.append(key)
#        super(MetahashManagedResource,self).__init__( **kwargs)
#
#    def reset_filtering_and_ordering(self):
#        logger.info(str(('---------reset filtering and ordering', self.scope)))
#        self.Meta.filtering = {}
#        self.Meta.ordering = []
#        metahash = MetaHash.objects.get_and_parse(scope=self.scope)
#        for key,hash in metahash.items():
#            if 'filtering' in hash and hash['filtering']:
#                self.Meta.filtering[key] = ALL
#        
#        for key,hash in metahash.items():
#            if 'ordering' in hash and hash['ordering']:
#                self.Meta.ordering.append(key)
#        logger.info(str(('---------reset filtering and ordering, done', self.scope)))
#
#    def build_schema(self):
#        logger.info('--- build_schema: ' + self.scope )
#        schema = super(MetahashManagedResource,self).build_schema()  # obligatory super call, this framework does not utilize
#        metahash = MetaHash.objects.get_and_parse(scope=self.scope)
#
#        for key, value in metahash.items():
#            if key not in schema['fields']:
##                logger.info('creating a virtual field: ' + key)
#                schema['fields'][key] = {}
#            schema['fields'][key].update(value)
#            
#                
#        try:
#            logger.info(str(('trying to locate resource information', self._meta.resource_name, self.scope)))
#            resource_def = MetaHash.objects.get(scope='resource', key=self._meta.resource_name)
#        
#            schema['resource_definition'] = resource_def.model_to_dict(scope='fields:resource')
#        except Exception, e:
#            logger.warn(str(('on trying to locate resource information', e, self._meta.resource_name)))
#                
#        logger.info('--- build_schema, done: ' + self.scope )
#        return schema        
        
#class UserResource1(JsonAndDatabaseResource):
#    permissions = fields.CharField()
#    
#    class Meta:
#        queryset = User.objects.all()
#        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
#        resource_name='user' 
#        authorization= Authorization()        
#        
#        excludes = ['password']
#        ordering = []
#        filtering = { }
#        serializer = LimsSerializer()
#        
#    def __init__(self, **kwargs):
#        super(UserResource,self).__init__(resource='user',**kwargs)
#
#    def build_schema(self):
#        schema = super(UserResource,self).build_schema()
##        temp = [ x.ref_resource_name for x in self.Meta.queryset.distinct('ref_resource_name')]
##        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'ref_resource_name', 'options': temp }
#        return schema
#
#    # permissions is a nested property - not exposing it as a resource - 
#    # so hydrate/dehydrate will be done custom
#    def dehydrate_permissions(self, bundle):
#        logger.info('dehydrate permissions')
#        try:
#            screensaverUser = bundle.obj.screensaveruser
#            if screensaverUser:
#                permissions = [ model_to_dict(x, exclude='id') for x in screensaverUser.permissions.all()]
#                logger.info(str(('created permissions', permissions)))
#                return permissions
#        except Exception, e:
#            logger.warn(str(('error accessing the screensaver user element', e)))
#    
#    def hydrate_permissions(self, bundle):
#        logger.info(str(('hydrate_permissions')))
#        try:
#            screensaverUser = bundle.obj.screensaveruser
#            
#            permission_array = bundle.data['permissions']
#            for p in permission_array:
#                try:
#                    permission = Permission.objects.get(**p)
#                except ObjectDoesNotExist, e:
#                    # create the permission on the fly, if the resource it refers to exists
#                    # todo: make sure that the Resource can also be found
#                    resource = MetaHash.objects.get(scope=p['scope'], key=p['key'])
#                    permission = Permission(**p)
#                    permission.save()
#                screensaverUser.permissions.add(permission)
#                screensaverUser.save();
#                logger.info(str(('permission saved', permission)))
#        except Exception, e:
#            logger.warn(str(('error accessing the screensaver user element on hydrate', e)))
#        return bundle
#        
#    def prepend_urls(self):
#        return [
#            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            url(r"^(?P<resource_name>%s)/(?P<username>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            ]
    
#        
#class MetaHashResource1(LoggingMixin, JsonAndDatabaseResource):
#
#    def __init__(self, **kwargs):
#        super(MetaHashResource,self).__init__(resource='metahash', **kwargs)
#
#    class Meta:
#        scope='fields:metahash'
#        queryset = MetaHash.objects.filter(scope__startswith="fields:").order_by('scope','ordinal','key')
#        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
#        authorization= Authorization()        
#        ordering = []
#        filtering = {} #{'scope':ALL, 'key':ALL}
#        serializer = LimsSerializer()
#        excludes = [] #['json_field']
#        always_return_data = True # this makes Backbone happy
#
#    def build_schema(self):
#        schema = super(MetaHashResource,self).build_schema()
#        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
#        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
#        return schema
#    
#
#    def obj_get_list(self, request=None, **kwargs):
#        objs = self.get_object_list(request)
#        return objs
#    
#    # because the metahash resource is both a resource and the definer of json fields, 
#    # reset_field_defs after each create/update, in case, new json fields are defined,
#    # or in case ordering,filtering groups are updated
#    def obj_create(self, bundle, **kwargs):
#        bundle = super(MetaHashResource, self).obj_create(bundle, **kwargs);
#        # TODO: there may be a better post-create hook?
#        if getattr(bundle.obj,'scope') == self.scope: #'fields:metahash':
#            self.reset_field_defs();
##        elif getattr(bundle.obj,'scope') == 'screensaveruser:fields':
##            ScreensaverUserResource().reset_filtering_and_ordering();
#        return bundle
#
#    def obj_update(self, bundle, **kwargs):
#        bundle = super(MetaHashResource, self).obj_update(bundle, **kwargs);
#        # TODO: there may be a better post-create hook?
#        if getattr(bundle.obj,'scope') == self.scope: #'fields:metahash':
#            if getattr(bundle.obj,'json_field_type'):
#                self.reset_field_defs();
#        return bundle
#
#    def hydrate(self, bundle):
#        bundle = super(MetaHashResource, self).hydrate(bundle);
#        logger.info(str(('=============hydrate:', bundle.data)))
#        return bundle
#
#    def rollback(self, bundles):
#        pass        
            
#class JsonAndDatabaseResource(PostgresSortingResource):
#    '''
#    This is a Resource wherein the fields are specified in the "Fields" store 
#        (the "Fields" store is the endpoint: 
#        /reports/api/v1/metahash/fields:metahash/[field_name], implemented with
#         the  "MetahashResource" api endpoint defined in this file )
#    -- tastypie.resources.ModelResource creates Resource.fields for all fields in the 
#        underlying Meta.queryset (TODO: unless using "excludes" "includes")
#    -- fields may be defined on the subclass also, as usual. 
#    -- to be usable and included in the "schema" endpoint, fields must be 
#        defined in the "Fields" store (TODO: exclude undefined fields)
#    -- fields records in the Resource store that define the "json_field_type"
#        value are fields that exist in the "json_field" of the table. 
#    -- to be usable in the (javascript) UI, fields must be defined in the 
#        Resource store 
#        (endpoint: /reports/api/v1/metahash/fields:metahash/[field_name]) 
#    '''
#    
#    
#    def __init__(self, resource=None, field_definition_scope='fields:metahash', **kwargs):
#        assert resource != None, 'resource kwarg must be defined' 
#        self.resource = resource
#        self.scope = 'fields:' + resource
#        self.field_definition_scope = field_definition_scope
#        logger.info(str(('---init resource', self.resource, self.scope, field_definition_scope)))
#
#        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
#        metahash = MetaHash.objects.get_and_parse(scope=self.scope, 
#                                                  field_definition_scope=field_definition_scope)
#        for key,hash in metahash.items():
#            if 'filtering' in hash and hash['filtering']:
#                self.Meta.filtering[key] = ALL
#        
#        for key,hash in metahash.items():
#            if 'ordering' in hash and hash['ordering']:
#                self.Meta.ordering.append(key)
#        
#        super(JsonAndDatabaseResource,self).__init__(**kwargs)
#        self.original_fields = deepcopy(self.fields)
#
#        self.field_defs = {}
#
##        self.reset_field_defs()
##        self.get_field_defs(scope)
#        
#        
#        logger.info(str(('---init resource, done', self.resource, field_definition_scope)))
#        
#    def prepend_urls(self):
#        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" 
#        # [ any word, except "schema" ]
#        # also note the double underscore "__" is because we also don't want to match in the first clause.
#        # We don't want "schema" since that reserved word is used by tastypie 
#        # for the schema definition for the resource (used by the UI)
#        return [
#            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.-:]+)/(?P<key>[\w\d_.-]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            url(r"^(?P<resource_name>%s)/(?P<key>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#        ]    
#    
#    # TODO: allow turn on/off of the reset methods for faster loading.
#    def reset_field_defs(self):
##        if(settings.DO_NOT_RESET_FIELDS):
##            pass # TODO: turn off when loading for speed
#        logger.info('----------reset_field_defs, ' + self.scope)
#        self.fields = deepcopy(self.original_fields)
#        self.field_defs = {}
##        self.reset_filtering_and_ordering()
#        
#    # TODO: allow turn on/off of the reset methods for faster loading.
#    def reset_filtering_and_ordering(self):
#        self.Meta.filtering = {}
#        self.Meta.ordering = []
#        metahash = MetaHash.objects.get_and_parse(scope=self.scope, clear=True)
#        for key,hash in metahash.items():
#            if 'filtering' in hash and hash['filtering']:
#                self.Meta.filtering[key] = ALL
#        
#        for key,hash in metahash.items():
#            if 'ordering' in hash and hash['ordering']:
#                self.Meta.ordering.append(key)
#    
#
#    def get_field_defs(self, scope):
#        ''''
#        Accomplishes two things:
#        1. creates Tastypie resource.field definitions for "json" fields defined in the metahash for this table's scope.
#            - json fields are fields that exist in the json_field, or, fields composed of other fields (i.e. full name)
#        2. returns a hash of all field definitions for the tastypie fields on this resource.  
#            - used to generate the public "schema"
#            - used to determine "filterable" and "orderable" fields
#        '''
#        # dynamically define fields according to what is defined in the metahash:
#        # - includes fields that are stored in the JSON field.  
#        # Note that: - for a new database, there will be no fields in the JSON field initially, so they have 
#        # to be populated before they can be filled with data. (cannot define and fill at the same time).
#        # TODO: allow turn on/off of the reset methods for faster loading.
#        logger.info(str(('get_field_defs', self.field_defs.keys())))
#       
#        if not self.field_defs:
#            logger.debug('------get_field_defs: ' + scope)
#            self.field_defs = {}
#            # 1. query the metahash for fields defined for this scope
#            for fieldinformation in MetaHash.objects.all().filter(scope=scope):
#                logger.debug('---- meta field for scope: ' + scope + ', ' + fieldinformation.key)
#                self.field_defs[fieldinformation.key] = {}
#                if fieldinformation.is_json():
#                    if fieldinformation.key in self.fields:
#                        raise DatabaseError('Illegal definition of a json_field with the same name as a database field on this scope: ' + scope + ',' + fieldinformation.key)
#                    self.field_defs[fieldinformation.key].update({ 'json_field_type': fieldinformation.json_field_type })
##                else if fieldinformation.is_virtual():
##                    pass
#                else:
#                    pass
#                    # assume, for now, that the field is "virtual", i.e. it is calculated, composite or otherwise created in the Resource.
#                    # Note: may create a flag to designate fields as virtual
##                    if fieldinformation.key not in self.fields:
##                        raise DatabaseError('Illegal definition of a non-json_field, non-virtual field that is not a resource field on this scope: ' + scope + ',' + fieldinformation.key)
#            
#            # 2. add in model/resource fields not specified by a metahash entry
#            # this is useful because for the first load, nothing is defined for the scope=metash:fields
#            for resource_field_key in self.fields.keys():
#                if resource_field_key not in self.field_defs: # and not key == 'json_field': # TODO: should exclude the json_field from this list, right?
#                    logger.debug('--------------- create default field for ' + resource_field_key)
#                    self.field_defs[resource_field_key] = {}
#            
#            # 3. create the virtual fields as Tastypie resource fields        
#            for schema_key, schema_value in self.field_defs.items():
#                # if the field is a json type, eval the json_field_type into a tastypie fields.Field definition so that tastypie knows about it
#                if schema_value.get('json_field_type'):
#                    # one more, redundant check
#                    if schema_key in self.fields:
#                        # TODO: move these validation errors to the client create process
#                        raise DatabaseError('Illegal attempt to define a json_field with the same name as a database field on this scope: ' + scope + ',' + + schema_key)
#                    self.fields[schema_key] = eval(schema_value['json_field_type'])(attribute=schema_key,readonly=True, blank=True, null=True) # trying to pass the fields.Field object args
#            logger.info(str(('----- tastypie fields created', self.fields.keys())))
#            # query metahash for schema defs (attributes) of all the fields 
#            
#            # 4. fill in all of the field definition information for the schema specification
##            TODO: replace this with a call to reports.models.MetaManager.get_and_parse(scope='fields:metahash', field_definition_scope='fields:metahash');
#            # The preceding logic could be simplified if we defined a "base_fields =['key', 'scope', 'ordinal', 'json_field']", 
#            # since then we don't have to wait for them to be defined 
#            # 
#            schema_definitions = MetaHash.objects.get_and_parse(scope=scope, field_definition_scope='fields:metahash')
#            for schema_key,schema_value in self.field_defs.items():
#                if schema_key in schema_definitions:
#                    schema_value.update( schema_definitions[schema_key] )
#                else:
#                    logger.info(str(('no field def found for: ', schema_key)))
#
#            logger.debug(str(('---- tastypie fields created: ', self.field_defs)))
#            self.reset_filtering_and_ordering()
#            logger.debug('------get_field_defs, done: ' + scope)
#            
#        return self.field_defs
#    
#    def build_schema(self):
#        logger.info('------build_schema: ' + self.scope)
#        try:
#            local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
#            
#            # TODO: ought to just create the schema from scratch; here relying on tastypie to give the schema
#            ## caveat: this has the nice side effect of verifying that any field defined is actually known by tastypie (for serialization hooks)
#            schema = super(JsonAndDatabaseResource,self).build_schema()
#
#            for key, value in local_field_defs.items():
#                if not key in schema['fields']:
#                    schema['fields'][key] = {}
#                schema['fields'][key].update(value)
#                # help for fields not yet defined
#                if not schema['fields'][key].get('ui_type'):
#                    pass
##                    schema['fields'][key]['ui_type'] = schema['fields'][key].get('type') # TODO: this is vestigial, correct?
#            
#            if 'json_field' in schema['fields']: 
#                schema['fields'].pop('json_field')  # because we don't want this serialized directly (see dehydrate)
#        
#            logger.info(str(('trying to locate resource information', self._meta.resource_name, self.scope)))
#            resource_def = MetaHash.objects.get(scope='resource', key=self._meta.resource_name)
#        
#            schema['resource_definition'] = resource_def.model_to_dict(scope='fields:resource')
#        except Exception, e:
#            logger.warn(str(('on building schema', e, self._meta.resource_name)))
#            exc_type, exc_obj, exc_tb = sys.exc_info()
#            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
#            logger.error(str((exc_type, fname, exc_tb.tb_lineno)))
#            raise e
#            
#        logger.info('------build_schema,done: ' + self.scope)
#        return schema
#
#    
#    def build_key(self, resource_name, data):
#        try:
#            resource_def = MetaHash.objects.get(scope='resource', key=resource_name)
#            resource = resource_def.model_to_dict(scope='fields:resource')
#            logger.info(str(('found resource', resource)))
#            
#            key_bits = []
#            for x in resource['id_attribute']:
#                if hasattr(data,x ):
#                    key_bits.append(getattr(data,x))
#                else:
#                    key_bits.append(data[x])
#            return  "/".join(key_bits)
#        except Exception, e:
#            logger.warn(str(('unable to locate resource information[id_attribute]; has it been loaded yet for this resource?', resource_name, e)))
#
#    URL_BASE = '/reports/api'
#
#    def build_resource_uri(self, resource, data):
#        new_key = self.build_key(resource, data)
#        if new_key:
##            base_uri = self.get_resource_uri()
##            if base_uri[len(base_uri)-1] != '/':
##                base_uri += '/'
##            return base_uri + new_key
#            return self.URL_BASE + '/' + self.api_name + '/'  + resource + '/' + new_key
#
#    def obj_resource_uri(self, resource, obj):
#        new_key = self.build_key(resource, obj)
#        if new_key:
#            return self.URL_BASE + '/' + self.api_name + '/'  + resource + '/' + new_key
##            base_uri = self.get_resource_uri()
##            if base_uri[len(base_uri)-1] != '/':
##                base_uri += '/'
##            return base_uri + new_key
#
#    def dehydrate(self, bundle):
#        logger.info('------dehydrate: ' + self.scope)
##        bundle = super(JsonAndDatabaseResource, self).dehydrate(bundle);
#
#        if len(bundle.data) == 0 : return bundle
#        
#        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
#        for key in [ x for x,y in local_field_defs.items() if y.get('json_field_type') ]:
#            bundle.data[key] = bundle.obj.get_field(key);
#        
#        bundle.data['json_field'] = ''
#        bundle.data.pop('json_field') # json_field will not be part of the public API, it is for internal use only
#        
#        # override the resource_uri, since we want to export the permanent composite key
#        bundle.data['resource_uri'] = self.build_resource_uri(self.resource, bundle.data) or bundle.data['resource_uri']
#
#        return bundle
#    
#    def hydrate_json_field(self, bundle):
#        '''
#        hydrate bundle data values that will be stuffed into the json_field
#        -Note: as mentioned elsewhere, for the initial load of the Metahash:fields, fields that
#        are JSON fields (to be stuffed into json_field) must be first defined as a metahash:field,
#        then they can be set as attribute values on other fields in the second step.
#        '''
##        bundle = super(JsonAndDatabaseResource, self).hydrate(bundle);
#        
#        logger.debug(str(('hydrate', bundle)))
#        
#        json_obj = {}
#        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
#        logger.debug(str(('local_field_defs',local_field_defs)))
#        
#        # Use the tastypie field type that has been designated to convert each field in the json stuffed field just like it were a real db field
#        for key in [ str(x) for x,y in local_field_defs.items() if 'json_field_type' in y and y['json_field_type'] ]:
#            val = bundle.data.get(key,None)
#            if val:
#                try:
#                    logger.info(str(('get value for key: ', key, val )))
#                    if hasattr(val, "strip"): # test if it is a string
#                        val = self.fields[key].convert(smart_text(val,'utf-8', errors='ignore'))
#                        logger.info(str(('got value for key: ', key, val)))
#                    elif hasattr(val, "__getitem__") or hasattr(val, "__iter__"): # test if it is a sequence
#                        val = [smart_text(x,'utf-8', errors='ignore') for x in val]
#                    logger.info(unicode(('got',key, val)))
#                    json_obj.update({ key:val })
#                except Exception:
#                    extype, ex, tb = sys.exc_info()
#                    formatted = traceback.format_exception_only(extype, ex)[-1]
#                    e =  RuntimeError, unicode(('failed to convert', key, 'with value', val, 
#                                             'with tastypie field type', type(self.fields[key]), 'message', formatted)).replace("'",""), 
#                    logger.warn(unicode(('throw', e, tb.tb_frame.f_code.co_filename, 'line', tb.tb_lineno)))
#                    raise e
#        bundle.data['json_field'] = json.dumps(json_obj);
#        logger.debug(str(('--- hydrated:', bundle)))
#        return bundle;

#
#class ManagedResourceDeclarativeMetaclass(ModelDeclarativeMetaclass):
#    def __new__(cls, name, bases, attrs):
#        new_class = super(ManagedResourceDeclarativeMetaclass, cls).__new__(cls, name, bases, attrs)
#        logger.info(str(('++++++++++++', cls, name, bases, attrs)))
#        
#        managed_field_defs = new_class.get_managed_fields() # so, scope is needed in meta 
#    
#        # pick managed fields from the generated fields
#        if managed_field_defs:
#            new_base_fields = {}
#            old_base_fields = attrs['base_fields']
#            logger.info(str(('old_base_fields', old_base_fields)))
#            attrs['base_fields'] = new_base_fields
#            for field_name, field_object in old_base_fields.items():
#                if field_name in managed_field_defs: 
#                    new_base_fields[field_name] = deepcopy(field_object)
#            
#            # create the fields not yet created
#            
#            logger.info(str(('fields defined', new_base_fields.keys())))
#            keys = set(managed_field_defs.keys()).difference(set(new_base_fields.keys()))
#            logger.info(str(('fields to be defined', keys)))
#            for name in keys:
#                schemaDict = managed_field_defs[name]
#                logger.info(str(('processing', schemaDict)))
#                if 'json_field_type' in schemaDict and schemaDict['json_field_type']:
#                    # TODO: use type to create class instances
#                    new_base_fields[name] = eval(managed_field_defs[name]['json_field_type'])(attribute=name,readonly=True, blank=True, null=True) # trying to pass the fields.Field object args
#                else:
#                    logger.info('creating unknown field as a char: ' + name)
#                    new_base_fields[name] = fields.CharField(attribute=name, readonly=True, blank=True, null=True)
##                    msg = 'no field defined for ' + name + ', at this time fields must be either json_fields, or defined normally as class fields'
##                    raise Exception(msg)
#    
#            meta = attrs.get('Meta')
#            filtering = getattr(meta, 'filtering')    
#            for field_name, field_object in new_base_fields.items():
#                filtering[field_name] = ALL
#
#                
#        return new_class
#    
#class ManagedResource1(LoggingMixin, ModelResource):
#
#    __metaclass__ = ManagedResourceDeclarativeMetaclass
#    
#    @classmethod
#    def get_managed_fields(cls, fields=None, excludes=None):
#        """
#        """
#        
#        logger.info(str(('calling get_managed_fields for ', cls)))
#        opts = getattr(cls, 'Meta', None)
#        print 'options===========', dir(opts)
#        if opts and  hasattr(opts, 'scope'):        
#            scope = getattr(opts, 'scope')
#            print '=== scope', scope
#            final_fields = {}
#            fields = fields or []
#            excludes = excludes or []
#    
#            # 4. fill in all of the field definition information for the schema specification
#            # The preceding logic could be simplified if we defined a "base_fields =['key', 'scope', 'ordinal', 'json_field']", 
#            # since then we don't have to wait for them to be defined 
#            # 
#            schema_definitions = MetaHash.objects.get_and_parse(scope=scope, field_definition_scope='fields:metahash')
#            return schema_definitions
#        else: # TODO: better
#            logger.warn('"scope" must be defined in attributes')
#
#        
#    
#
#    def __init__(self, **kwargs):
#        self.scope = self._meta.scope
#        self.resource = self._meta.resource_name
#        logger.info('========scope: ' + self.scope)
#        super(ManagedResource1,self).__init__(**kwargs)   
#        
#        logger.info(str(('======= managed fields', self.fields))) 
#    
#    def build_schema(self):
#        logger.info('------build_schema: ' + self.scope)
#        try:
#            schema = {}
#            
#            local_field_defs = MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope='fields:metahash')
#            schema['fields'] = local_field_defs
#            
#            logger.info(str(('trying to locate resource information', self._meta.resource_name, self.scope)))
#            resource_def = MetaHash.objects.get(scope='resource', key=self._meta.resource_name)
#        
#            schema['resource_definition'] = resource_def.model_to_dict(scope='fields:resource')
#            logger.info('------build_schema,done: ' + self.scope)
#            return schema    
#        except Exception, e:
#            logger.warn(str(('on building schema', e, self._meta.resource_name)))
#            exc_type, exc_obj, exc_tb = sys.exc_info()
#            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
#            logger.error(str((exc_type, fname, exc_tb.tb_lineno)))
#            raise e
#            
#    def build_key(self, resource_name, data):
#        try:
#            resource_def = MetaHash.objects.get(scope='resource', key=resource_name)
#            resource = resource_def.model_to_dict(scope='fields:resource')
#            logger.info(str(('found resource', resource)))
#            
#            key_bits = []
#            for x in resource['id_attribute']:
#                if hasattr(data,x ):
#                    key_bits.append(getattr(data,x))
#                else:
#                    key_bits.append(data[x])
#            return  "/".join(key_bits)
#        except Exception, e:
#            logger.warn(str(('unable to locate resource information[id_attribute]; has it been loaded yet for this resource?', resource_name, e)))
#
#    URL_BASE = '/reports/api'
#
#    def build_resource_uri(self, resource, data):
#        new_key = self.build_key(resource, data)
#        if new_key:
##            base_uri = self.get_resource_uri()
##            if base_uri[len(base_uri)-1] != '/':
##                base_uri += '/'
##            return base_uri + new_key
#            return self.URL_BASE + '/' + self.api_name + '/'  + resource + '/' + new_key
#
#    def obj_resource_uri(self, resource, obj):
#        new_key = self.build_key(resource, obj)
#        if new_key:
#            return self.URL_BASE + '/' + self.api_name + '/'  + resource + '/' + new_key
##            base_uri = self.get_resource_uri()
##            if base_uri[len(base_uri)-1] != '/':
##                base_uri += '/'
##            return base_uri + new_key
#    
#    def prepend_urls(self):
#        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" 
#        # [ any word, except "schema" ]
#        # also note the double underscore "__" is because we also don't want to match in the first clause.
#        # We don't want "schema" since that reserved word is used by tastypie 
#        # for the schema definition for the resource (used by the UI)
#        return [
#            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.-:]+)/(?P<key>[\w\d_.-]+)%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#            url(r"^(?P<resource_name>%s)/(?P<key>((?=(schema))__|(?!(schema))[\w\d_.-]+))%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#        ]  
#        
#    def dehydrate(self, bundle):
#        logger.info('------dehydrate: ' + self.scope)
##        bundle = super(JsonAndDatabaseResource, self).dehydrate(bundle);
#
#        if len(bundle.data) == 0 : return bundle
#
#        local_field_defs = MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope='fields:metahash')
#        logger.info('2a')
#        for key in [ x for x,y in local_field_defs.items() if y.get('json_field_type') ]:
#            bundle.data[key] = bundle.obj.get_field(key);
#        logger.info('2b')
#        
#        bundle.data['json_field'] = ''
#        bundle.data.pop('json_field') # json_field will not be part of the public API, it is for internal use only
#        logger.info('2c')
#
#        # override the resource_uri, since we want to export the permanent composite key
#        bundle.data['resource_uri'] = self.build_resource_uri(self.resource, bundle.data) or bundle.data['resource_uri']
#        logger.info('2d')
#
#        return bundle
#    
#    def hydrate_json_field(self, bundle):
#        logger.info('------hydrate: ' + self.scope)
#        '''
#        hydrate bundle data values that will be stuffed into the json_field
#        -Note: as mentioned elsewhere, for the initial load of the Metahash:fields, fields that
#        are JSON fields (to be stuffed into json_field) must be first defined as a metahash:field,
#        then they can be set as attribute values on other fields in the second step.
#        '''
##        bundle = super(JsonAndDatabaseResource, self).hydrate(bundle);
#        
#        logger.debug(str(('hydrate', bundle)))
#        
#        json_obj = {}
#        local_field_defs = MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope='fields:metahash', clear=True)
##        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
#        logger.debug(str(('local_field_defs',local_field_defs)))
#        
#        # Use the tastypie field type that has been designated to convert each field in the json stuffed field just like it were a real db field
#        for key in [ str(x) for x,y in local_field_defs.items() if 'json_field_type' in y and y['json_field_type'] ]:
#            val = bundle.data.get(key,None)
#            if val:
#                try:
#                    logger.info(str(('get value for key: ', key, val )))
#                    if hasattr(val, "strip"): # test if it is a string
#                        val = self.fields[key].convert(smart_text(val,'utf-8', errors='ignore'))
#                        logger.info(str(('got value for key: ', key, val)))
#                    elif hasattr(val, "__getitem__") or hasattr(val, "__iter__"): # test if it is a sequence
#                        val = [smart_text(x,'utf-8', errors='ignore') for x in val]
#                    logger.info(unicode(('got',key, val)))
#                    json_obj.update({ key:val })
#                except Exception:
#                    extype, ex, tb = sys.exc_info()
#                    formatted = traceback.format_exception_only(extype, ex)[-1]
#                    e =  RuntimeError, unicode(('failed to convert', key, 'with value', val, 
#                                             'with tastypie field type', type(self.fields[key]), 'message', formatted)).replace("'",""), 
#                    logger.warn(unicode(('throw', e, tb.tb_frame.f_code.co_filename, 'line', tb.tb_lineno)))
#                    raise e
#        bundle.data['json_field'] = json.dumps(json_obj);
#        logger.debug(str(('--- hydrated:', bundle)))
#        return bundle; 
# 
        