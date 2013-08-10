import json
import re
import sys
import os

from copy import deepcopy
from django.db import DatabaseError
from django.conf.urls import url
from django.utils.encoding import smart_str
from django.utils import timezone

from django.forms.models import model_to_dict
from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL
from tastypie import fields # NOTE: required for dynamic resource field definitions using eval
from tastypie.resources import Resource

import time
from lims.api import PostgresSortingResource, CSVSerializer
from reports.models import MetaHash, Vocabularies, ApiLog

import logging
from django.core.exceptions import ObjectDoesNotExist
from tastypie.exceptions import NotFound
from tastypie.bundle import Bundle
        
logger = logging.getLogger(__name__)

class LoggingMixin(Resource):
    
    def obj_create(self, bundle, **kwargs):
        logger.info('----log obj_create')
        
        bundle = super(LoggingMixin, self).obj_create(bundle=bundle, **kwargs)
        
        log = ApiLog()
        log.username = bundle.request.user.username #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.user_id = bundle.request.user.id #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.date_time = timezone.now()
        log.resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
        log.uri,log.key = self.get_uri(bundle.request.get_full_path(), self._meta.resource_name, bundle.obj, bundle.data)
        log.save()
        
        return bundle    
    
    def obj_delete(self, bundle, **kwargs):
        logger.info('---log obj_delete')
        
        bundle = super(LoggingMixin, self).obj_create(bundle=bundle, **kwargs)
        
        log = ApiLog()
        log.username = bundle.request.user.username #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.user_id = bundle.request.user.id #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.date_time = timezone.now()
        log.resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
        log.uri,log.key = self.get_uri(bundle.request.get_full_path(), self._meta.resource_name, bundle.obj, bundle.data)
        log.save()
        
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
    
    
    def get_uri(self, full_path, resource_name, obj, data):
        # replace the id with the "id attribute" ids (public ID's)
        logger.info('--- get_uri: '+ full_path)
        key = obj.id 
        
        # look up the resource definition
                        
        try:
            logger.info(str(('trying to locate resource information', resource_name, self.scope)))
            resource_def = MetaHash.objects.get(scope='resource', key=resource_name)
            logger.info(str(('create dict for obj', resource_def)))
            resource = resource_def.model_to_dict(scope='fields:resource')
            logger.info(str(('--- got', resource, str((resource['id_attribute'])) )))
        except Exception, e:
            logger.warn(str(('unable to locate resource information, has it been loaded yet for this resource?', resource_name, e)))
            return (full_path, key)

        try:            
            matchObject = re.match(r'(.*\/'+resource_name+'\/)(.*)', full_path)
            if(matchObject):
                if matchObject.group(2):
                    key = matchObject.group(2)
                    if key[len(key)-1] == '/': key = key[:len(key)-1]
                    if obj:
                        if (str((obj.id))) == key :
                            # replace the string
                            new_public_key = "/".join([getattr(obj,x) for x in resource['id_attribute']])
                            logging.info(str(('created new public key', new_public_key, 'for object to log:', obj, 'resource', resource_name)))
                            full_path = matchObject.group(1) + new_public_key + '/'
                            key = new_public_key
                        else:
                            logger.info(str(('id',obj.id,'doesnt equal match obj', matchObject.group())))
                    else:
                        new_public_key = "/".join([data[x] for x in resource['id_attribute']])
                        logging.info(str(('created new public key', new_public_key, 'for data to log:', data, 'resource', resource_name)))
                        full_path = matchObject.group(1) + new_public_key + '/'
                        key = new_public_key
                        
                else: # found the resource, but id was not given
                    if obj:
                        new_public_key = "/".join([getattr(obj,x) for x in resource['id_attribute']])
                        logging.info(str(('created new public key', new_public_key, 'for object to log:', obj, 'resource', resource_name)))
                        full_path = matchObject.group(1) + new_public_key + '/'
                        key = new_public_key
                    else: # no object, probably means this is a create
                        new_public_key = "/".join([data[x] for x in resource['id_attribute']])
                        logging.info(str(('created new public key', new_public_key, 'for data to log:', data, 'resource', resource_name)))
                        full_path = matchObject.group(1) + new_public_key + '/'
                        key = new_public_key
                        
            else:
                logger.warn(str(('non-standard resource, does not contain resource_name:', resource_name, full_path)))
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            logger.error(str((exc_type, fname, exc_tb.tb_lineno)))
            logger.error(str(('on trying to extract public key information from the object and data', resource_name, obj, data, e)))
        logger.info(str(('get uri,key:', full_path, key)))
        return (full_path, key)
              
    
    def obj_update(self, bundle, skip_errors=False, **kwargs):
        logger.info('--- log obj_update')
        original_bundle = Bundle(data=deepcopy(bundle.data))
        if hasattr(bundle,'obj'): original_bundle.obj = bundle.obj
        original_bundle = self._locate_obj(original_bundle, **kwargs)
        
        original_bundle = super(LoggingMixin, self).full_dehydrate(original_bundle)
        
        updated_bundle = super(LoggingMixin, self).obj_update(bundle=bundle, **kwargs)
        updated_bundle = super(LoggingMixin, self).full_dehydrate(updated_bundle)
        
        # TODO: diff updated_bundle.data
        original_keys = set(original_bundle.data.keys())
        updated_keys = set(updated_bundle.data.keys())
        
        intersect_keys = original_keys.intersection(updated_keys)
        
        log = ApiLog()
        log.username = bundle.request.user.username #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.user_id = bundle.request.user.id #self._meta.authentication.get_identifier(bundle.request) # see tastypie.authentication.Authentication
        log.date_time = timezone.now()
        log.resource_name = self._meta.resource_name
        log.api_action = str((bundle.request.method)).upper()
        log.uri,log.key = self.get_uri(bundle.request.get_full_path(), self._meta.resource_name, bundle.obj, bundle.data)
        
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
            
        log.save()

        logger.info('obj_update done')
        
        return updated_bundle
                
#    def dispatch(self, request_type, request, **kwargs):
#        logger.debug('%s %s %s' % (request.method, request.get_full_path(), request.raw_post_data))
# 
#        try:
#            response = super(LoggingMixin, self).dispatch(request_type, request, **kwargs)
#            
#            log = ApiLog()
#            log.uri = request.get_full_path()
#            log.save()
#            
#        except Exception, e:
#            if hasattr(e, 'response'):
#                logger.debug(
#                    'Response %s %s' %
#                    (e.response.status_code, e.response.content))
#            else:
#                logger.debug(str(('failed dispatch', e)))
#            raise
  
# Mixin class - note this class must be mixed in first, since the tastypie Resource class does not call super().__init__
# TODO: merge/refactor this with reports.JsonAndDatabaseResource
class MetahashManagedResource(object):
    def __init__(self, **kwargs):
#        if not self.scope:
#            self.scope = kwargs.pop('scope')
        logger.info(str(('---------init MetahashManagedResource', self.scope)))
        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
        metahash = MetaHash.objects.get_and_parse(scope=self.scope)
        
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)
        super(MetahashManagedResource,self).__init__( **kwargs)

    def reset_filtering_and_ordering(self):
        logger.info(str(('---------reset filtering and ordering', self.scope)))
        self.Meta.filtering = {}
        self.Meta.ordering = []
        metahash = MetaHash.objects.get_and_parse(scope=self.scope)
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)
        logger.info(str(('---------reset filtering and ordering, done', self.scope)))

    def build_schema(self):
        logger.info('--- build_schema: ' + self.scope )
        schema = super(MetahashManagedResource,self).build_schema()  # obligatory super call, this framework does not utilize
        metahash = MetaHash.objects.get_and_parse(scope=self.scope)

        for key, value in metahash.items():
            if key not in schema['fields']:
#                logger.info('creating a virtual field: ' + key)
                schema['fields'][key] = {}
            schema['fields'][key].update(value)
            
                
        try:
            logger.info(str(('trying to locate resource information', self._meta.resource_name, self.scope)))
            resource_def = MetaHash.objects.get(scope='resource', key=self._meta.resource_name)
        
            schema['resource_definition'] = resource_def.model_to_dict(scope='fields:resource')
        except Exception, e:
            logger.warn(str(('on trying to locate resource information', e, self._meta.resource_name)))
                
        logger.info('--- build_schema, done: ' + self.scope )
        return schema  
  
        
class JsonAndDatabaseResource(PostgresSortingResource):
    def __init__(self, scope=None, field_definition_scope='fields:metahash', **kwargs):
        self.scope = scope
        assert scope != None, 'scope kwarg must be defined' 
        self.field_definition_scope = field_definition_scope
        logger.info(str(('---init resource', scope, field_definition_scope)))

        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
        metahash = MetaHash.objects.get_and_parse(scope=self.scope, field_definition_scope=field_definition_scope)
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)

        logger.info(str(('+++filtering', self.Meta.filtering)))
        logger.info(str(('ordering', self.Meta.ordering)))
        #        self.reset_filtering_and_ordering()

        super(JsonAndDatabaseResource,self).__init__(**kwargs)
        self.original_fields = deepcopy(self.fields)
        self.field_defs = {}
        logger.info(str(('---init resource, done', scope, field_definition_scope)))
        
    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" allows us to match any word, except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie for the schema definition for the resource (used by the UI)
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.-:]+)/(?P<key>[\w\d_.-]+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<key>((?=(schema))__|(?!(schema))[\w\d_.-]+))/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
    
    # TODO: allow turn on/off of the reset methods for faster loading.
    def reset_field_defs(self):
#        if(settings.DO_NOT_RESET_FIELDS):
#            pass # TODO: turn off when loading for speed
        logger.info('----------reset_field_defs, ' + self.scope)
        self.fields = deepcopy(self.original_fields)
        self.field_defs = {}
#        self.reset_filtering_and_ordering()
        
    # TODO: allow turn on/off of the reset methods for faster loading.
    def reset_filtering_and_ordering(self):
        self.Meta.filtering = {}
        self.Meta.ordering = []
        metahash = MetaHash.objects.get_and_parse(scope=self.scope)
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)

        
    def get_field_defs(self, scope):
        # dynamically define fields that are stored in the JSON field.  
        # Note that: - for a new database, there will be no fields in the JSON field initially, so they have 
        # to be populated before they can be filled with data. (cannot define and fill at the same time).
        # TODO: allow turn on/off of the reset methods for faster loading.
        if not self.field_defs:
            logger.debug('------get_field_defs: ' + scope)
            self.field_defs = {}
            # first, query the metahash for fields defined for this scope
            for fieldinformation in MetaHash.objects.all().filter(scope=scope):
                logger.debug('---- meta field for scope: ' + scope + ', ' + fieldinformation.key)
                self.field_defs[fieldinformation.key] = {}
                if fieldinformation.is_json():
                    if fieldinformation.key in self.fields:
                        # TODO: move these validation errors to the client create process
                        raise DatabaseError('Illegal definition of a json_field with the same name as a database field on this scope: ' + scope + ',' + fieldinformation.key)
                    self.field_defs[fieldinformation.key].update({ 'json_field_type': fieldinformation.json_field_type })
                else:
                    if fieldinformation.key not in self.fields:
                        raise DatabaseError('Illegal definition of a non-json_field that is not a database field on this scope: ' + scope + ',' + fieldinformation.key)
            
            # add in model/resource fields not specified by a metahash entry
            for resource_field_key in self.fields.keys():
                if resource_field_key not in self.field_defs: # and not key == 'json_field': # TODO: should exclude the json_field from this list, right?
                    logger.debug('--------------- create default field for ' + resource_field_key)
                    self.field_defs[resource_field_key] = {}
            
            # add the virtual API fields        
            for schema_key, schema_value in self.field_defs.items():
                # if the field is a json type, eval the json_field_type into a tastypie fields.Field definition so that tastypie knows about it
                if schema_value.get('json_field_type'):
                    # one more, redundant check
                    if schema_key in self.fields:
                        # TODO: move these validation errors to the client create process
                        raise DatabaseError('Illegal attempt to define a json_field with the same name as a database field on this scope: ' + scope + ',' + + schema_key)
                    self.fields[schema_key] = eval(schema_value['json_field_type'])(attribute=schema_key,readonly=True, blank=True, null=True) # trying to pass the fields.Field object args
            logger.debug(str(('----- tastypie fields created', self.fields)))
            # query metahash for schema defs (attributes) of all the fields 
            for schema_key,schema_value in self.field_defs.items():
                logger.debug(str(('-----make schema:', self.scope, schema_key)))
                # now fill in meta information for the schema report to UI; using data from this table itself, either in real or json fields!
                for meta_record in MetaHash.objects.all().filter(scope='fields:metahash'):  # fields:metahash are defined for all reports
                    if meta_record.key == 'key':
                        continue # don't put these into the schema definitions, its too recursive confusing
                    schema_value.update({
                          meta_record.key: MetaHash.objects.get_or_none(scope=scope, key=schema_key, function=lambda x : (x.get_field(meta_record.key)) )
                          })
                    
                # now check if the field uses controlled vocabulary, look that up now.  TODO: "vocabulary_scope_ref" should be a constant
                # TODO: "vocabulary_scope_ref" needs to be created by default as a field:metahash; this argues for making it a "real" field
                if schema_value.get(u'vocabulary_scope_ref'):
                    logger.debug(str(('looking for a vocabulary', schema_value['vocabulary_scope_ref'] )))
                    schema_value['choices'] = [x.key for x in Vocabularies.objects.all().filter(scope=schema_value['vocabulary_scope_ref'])]
                    logger.debug(str(('got', schema_value['choices'] )))
                logger.debug(str(('----defined schema: ', schema_key, schema_value)))

            self.reset_filtering_and_ordering()
            logger.debug('------get_field_defs, done: ' + scope)
            
        return self.field_defs
    
    def build_schema(self):
        logger.info('------build_schema: ' + self.scope)
        try:
            local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
            
            # TODO: probably ought to just create the schema from scratch; here relying on tastypie to give the schema
            ## caveat: this has the nice side effect of verifying that any field defined is actually known by tastypie (for serialization hooks)
            schema = super(JsonAndDatabaseResource,self).build_schema()
            for key, value in local_field_defs.items():
                schema['fields'][key].update(value)
                # help for fields not yet defined
                if not schema['fields'][key].get('ui_type'):
                    schema['fields'][key]['ui_type'] = schema['fields'][key].get('type')
            
            schema['fields'].pop('json_field')
        
            logger.info(str(('trying to locate resource information', self._meta.resource_name, self.scope)))
            resource_def = MetaHash.objects.get(scope='resource', key=self._meta.resource_name)
        
            schema['resource_definition'] = resource_def.model_to_dict(scope='fields:resource')
        except Exception, e:
            logger.warn(str(('on building schema', e, self._meta.resource_name)))
            raise e
            
        logger.info('------build_schema,done: ' + self.scope)
        return schema
    
    def dehydrate(self, bundle):
#        logger.info('------dehydrate: ' + self.scope)
#        bundle = super(JsonAndDatabaseResource, self).dehydrate(bundle);
        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
        for key in [ x for x,y in local_field_defs.items() if y.get('json_field_type') ]:
            bundle.data[key] = bundle.obj.get_field(key);
        
        bundle.data['json_field'] = ''
        bundle.data.pop('json_field') # json_field will not be part of the public API, it is for internal use only
        
        # override the resource_uri, since we want to export the permanent composite key
        base_uri = self.get_resource_uri()
        if base_uri[len(base_uri)-1] != '/':
            base_uri += '/'
        bundle.data['resource_uri'] =  base_uri + bundle.data['scope'] + '/' + bundle.data['key'] +'/'
                
        # and don't send the internal id out for PATCH uses, at least
        # But: we _do_ have to send it out for Backbone, since we don't know how to use things like composite keys yet - sde4
        # bundle.data.pop('id')
        
#        logger.info('------dehydrate, done: ' + self.scope)
        return bundle
    
    def hydrate_json_field(self, bundle):
#        bundle = super(JsonAndDatabaseResource, self).hydrate(bundle);
        
        logger.debug(str(('hydrate', bundle)))
        
        # not sure if the bundl obj is id'd yet
        #        old_json_obj = bundle.obj.get_field_hash()
        
        # otherwise have to get it using query
        #obj = self.queryset.get(scope=bundle.data['scope'], key=bundle.data['key'])  # Note: what if trying to change the key?
        #old_json_obj = obj.get_field_hash() if obj else {}
        
        json_obj = {}
        field_reset_required = False
        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
        logger.debug(str(('local_field_defs',local_field_defs)))
        for key in [ str(x) for x,y in local_field_defs.items() if 'json_field_type' in y and y['json_field_type'] ]:
            json_obj.update({ key: bundle.data.get(key,None)})

        bundle.data['json_field'] = json.dumps(json_obj);
        logger.debug(str(('--- hydrated:', bundle)))
        return bundle;

    
class MetaHashResource(LoggingMixin, JsonAndDatabaseResource):

    def __init__(self, **kwargs):
        super(MetaHashResource,self).__init__(scope='fields:metahash', **kwargs)

    class Meta:
        queryset = MetaHash.objects.filter(scope__startswith="fields:").order_by('scope','ordinal','key')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        ordering = []
        filtering = {} #{'scope':ALL, 'key':ALL}
        serializer = CSVSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy

    def build_schema(self):
        schema = super(MetaHashResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
        return schema
    
#    # called by ModelResource in ModelResource.apply_filters 
#    def get_object_list(self, request):
#        obj_list = super(MetaHashResource, self).get_object_list(request);
#        # could apply filtering here
#        return obj_list
#    
#    def obj_get(self, bundle, **kwargs):
#        obj = super(MetaHashResource, self).obj_get(bundle, **kwargs)
#        return obj
#    
#    def dehydrate(self, bundle):
#        bundle = super(MetaHashResource, self).dehydrate(bundle);
#        return bundle

    # because the metahash resource is both a resource and the definer of json fields, 
    # reset_field_defs after each create/update, in case, new json fields are defined,
    # or in case ordering,filtering groups are updated
    def obj_create(self, bundle, **kwargs):
        bundle = super(MetaHashResource, self).obj_create(bundle, **kwargs);
        # TODO: there may be a better post-create hook?
        if getattr(bundle.obj,'scope') == self.scope: #'fields:metahash':
            self.reset_field_defs();
#        elif getattr(bundle.obj,'scope') == 'screensaveruser:fields':
#            ScreensaverUserResource().reset_filtering_and_ordering();
        return bundle

    def obj_update(self, bundle, **kwargs):
        bundle = super(MetaHashResource, self).obj_update(bundle, **kwargs);
        # TODO: there may be a better post-create hook?
        if getattr(bundle.obj,'scope') == self.scope: #'fields:metahash':
            if getattr(bundle.obj,'json_field_type'):
                self.reset_field_defs();
        return bundle

    def hydrate(self, bundle):
        bundle = super(MetaHashResource, self).hydrate(bundle);
        return bundle

    def rollback(self, bundles):
        pass        


class VocabulariesResource(LoggingMixin, JsonAndDatabaseResource):
    '''
    This resource extends the JsonAndDatabaseResource using a new table (vocabularies) but has
    fields defined in the Metahash table.
    '''
    def __init__(self, **kwargs):
        super(VocabulariesResource,self).__init__(scope='fields:vocabularies', **kwargs)

    class Meta:
        queryset = Vocabularies.objects.all().order_by('scope', 'ordinal', 'key')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = CSVSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
    
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

class ResourceResource(LoggingMixin, JsonAndDatabaseResource):
    '''
    This resource extends the JsonAndDatabaseResource, uses the metahash table internally, and has
    fields defined in the Metahash table.
    '''
    def __init__(self, **kwargs):
        super(ResourceResource,self).__init__(scope='fields:resource', field_definition_scope='fields:resource', **kwargs)

    class Meta:
        queryset = MetaHash.objects.filter(scope='resource').order_by('key', 'ordinal', 'scope')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = CSVSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='resource' # appears that tastypie needs this with a resource named "resource"!
    
    def build_schema(self):
        schema = super(ResourceResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('key')]
        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'key', 'options': temp }
        return schema

class ApiLogResource(MetahashManagedResource, PostgresSortingResource):
    '''
    '''
    class Meta:
        queryset = ApiLog.objects.all().order_by('resource_name', 'username','date_time')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'username':ALL, 'uri': ALL, 'resource_name':ALL}
        serializer = CSVSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='apilog' 
    
    def __init__(self, **kwargs):
        self.scope = 'fields:apilog'
        super(ApiLogResource,self).__init__(**kwargs)

    def build_schema(self):
        schema = super(ApiLogResource,self).build_schema()
        temp = [ x.resource_name for x in self.Meta.queryset.distinct('resource_name')]
        schema['extraSelectorOptions'] = { 'label': 'Resource', 'searchColumn': 'resource_name', 'options': temp }
        return schema
