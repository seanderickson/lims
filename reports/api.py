from django.db import DatabaseError
from django.conf.urls import url
from django.conf.urls import url
import json
from copy import deepcopy

from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL
from tastypie import fields

from lims.api import PostgresSortingResource, CSVSerializer
from reports.models import FieldInformation, MetaHash, Vocabularies

import logging
from db.api import ScreensaverUserResource
        
logger = logging.getLogger(__name__)
        
class JsonAndDatabaseResource(PostgresSortingResource):
    def __init__(self, scope=None, **kwargs):
        self.scope = scope

        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
        metahash = MetaHash.objects.get_metahash(scope=scope)
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
        
        
    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" allows us to match any word, except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie for the schema definition for the resource (used by the UI)
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.-:]+)/(?P<key>[\w\d_.-]+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<key>((?=(schema))__|(?!(schema))[\w\d_.-]+))/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
    
    def reset_field_defs(self):
        logger.info('----------reset_field_defs, ' + self.scope)
        self.fields = deepcopy(self.original_fields)
        self.field_defs = {}
#        self.reset_filtering_and_ordering()
        
    def reset_filtering_and_ordering(self):
        self.Meta.filtering = {}
        self.Meta.ordering = []
        metahash = MetaHash.objects.get_metahash(scope=self.scope)
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)

        logger.info(str(('+++filtering', self.Meta.filtering)))
        logger.info(str(('ordering', self.Meta.ordering)))
        
    def get_field_defs(self, scope):
        # dynamically define fields that are stored in the JSON field.  
        # Note that: - for a new database, there will be no fields in the JSON field initially (naturally)
        #     - also, all JSON field values will be null to start (naturally)
        # ALSO NOTE: creating new JSON fields requires a restart (to reload) (todo: fix that)
        if not self.field_defs:
            logger.info('------get_field_defs: ' + scope)
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

            # query metahash for schema defs for all the fields
            for schema_key,schema_value in self.field_defs.items():
                logger.debug(str(('-----make schema:', self.scope, schema_key)))
                # now fill in meta information for the schema report to UI; using data from this table itself, either in real or json fields!
                for meta_record in MetaHash.objects.all().filter(scope='metahash:fields'):  # metahash:fields are defined for all reports
                    if meta_record.key == 'key':
                        continue # don't put these into the schema definitions, its too recursive confusing
                    schema_value.update({
                          meta_record.key: MetaHash.objects.get_or_none(scope=scope, key=schema_key, function=lambda x : (x.get_field(meta_record.key)) )
                          })
                    
                # now check if the field uses controlled vocabulary, look that up now.  TODO: "vocabulary_scope_ref" should be a constant
                # TODO: "vocabulary_scope_ref" needs to be created by default as a metahash:field; this argues for making it a "real" field
                if schema_value.get(u'vocabulary_scope_ref'):
                    logger.debug(str(('looking for a vocabulary', schema_value['vocabulary_scope_ref'] )))
                    schema_value['choices'] = [x.key for x in Vocabularies.objects.all().filter(scope=schema_value['vocabulary_scope_ref'])]
                    logger.debug(str(('got', schema_value['choices'] )))
                logger.debug(str(('----defined schema: ', schema_key, schema_value)))

            self.reset_filtering_and_ordering()
            
        return self.field_defs
    
    def build_schema(self):
        logger.info('------build_schema: ' + self.scope)
        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
        schema = super(JsonAndDatabaseResource,self).build_schema()
        
        for key, value in local_field_defs.items():
            schema['fields'][key].update(value)
            # help for fields not yet defined
            if not schema['fields'][key].get('ui_type'):
                schema['fields'][key]['ui_type'] = schema['fields'][key].get('type')
        
        schema['fields'].pop('json_field')
        
        return schema
    
    def dehydrate(self, bundle):
#        bundle = super(JsonAndDatabaseResource, self).dehydrate(bundle);
        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
        for key in [ x for x,y in local_field_defs.items() if y.get('json_field_type') ]:
            bundle.data[key] = bundle.obj.get_field(key);
        
        bundle.data['json_field'] = ''
        bundle.data.pop('json_field') # json_field will not be part of the public API, it is for internal use only
        
        # override the resource_uri, since we want to export the permanent composite key
        bundle.data['resource_uri'] = self.get_resource_uri() + bundle.data['scope'] + '/' + bundle.data['key'] +'/'
                
        # and don't send the internal id out for PATCH uses, at least
        # But: we _do_ have to send it out for Backbone, since we don't know how to use things like composite keys yet - sde4
        # bundle.data.pop('id')
        
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
        for key in [ str(x) for x,y in local_field_defs.items() if 'json_field_type' in y ]:
            json_obj.update({ key: bundle.data.get(key,None)})

        bundle.data['json_field'] = json.dumps(json_obj);
        logger.debug(str(('--- hydrated:', bundle)))
        return bundle;

    
class MetaHashResource(JsonAndDatabaseResource):

    def __init__(self, **kwargs):
        super(MetaHashResource,self).__init__(scope='metahash:fields', **kwargs)

    class Meta:
        queryset = MetaHash.objects.all().order_by('scope','ordinal','key')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {} #{'scope':ALL, 'key':ALL}
        serializer = CSVSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy

    def build_schema(self):
        schema = super(MetaHashResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['searchTerms'] = temp
        return schema
    
    # called by ModelResource in ModelResource.apply_filters 
    def get_object_list(self, request):
        obj_list = super(MetaHashResource, self).get_object_list(request);
        # could apply filtering here
        return obj_list
    
    def obj_get(self, bundle, **kwargs):
        obj = super(MetaHashResource, self).obj_get(bundle, **kwargs)
        return obj
    
    def dehydrate(self, bundle):
        bundle = super(MetaHashResource, self).dehydrate(bundle);
        return bundle

    # because the metahash resource is both a resource and the definer of json fields, 
    # reset_field_defs after each create/update, in case, new json fields are defined,
    # or in case ordering,filtering groups are updated
    def obj_create(self, bundle, **kwargs):
        bundle = super(MetaHashResource, self).obj_create(bundle, **kwargs);
        # TODO: there may be a better post-create hook?
        if getattr(bundle.obj,'scope') == 'metahash:fields':
            self.reset_field_defs();
        elif getattr(bundle.obj,'scope') == 'screensaveruser:fields':
            ScreensaverUserResource().reset_filtering_and_ordering();
        return bundle

    def obj_update(self, bundle, **kwargs):
        bundle = super(MetaHashResource, self).obj_update(bundle, **kwargs);
        # TODO: there may be a better post-create hook?
        if getattr(bundle.obj,'scope') == 'metahash:fields':
            if getattr(bundle.obj,'json_field_type'):
                self.reset_field_defs();
        return bundle

    def hydrate(self, bundle):
        bundle = super(MetaHashResource, self).hydrate(bundle);
        return bundle

    def rollback(self, bundles):
        pass        


class VocabulariesResource(JsonAndDatabaseResource):

    def __init__(self, **kwargs):
        super(VocabulariesResource,self).__init__(scope='vocabularies:fields', **kwargs)

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
    
    # called by ModelResource in ModelResource.apply_filters 
    def get_object_list(self, request):
        obj_list = super(VocabulariesResource, self).get_object_list(request);
        # could apply filtering here
        return obj_list
    
    def obj_get(self, bundle, **kwargs):
        obj = super(VocabulariesResource, self).obj_get(bundle, **kwargs)
        return obj
    
    def build_schema(self):
        schema = super(VocabulariesResource,self).build_schema()
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['searchTerms'] = temp
        return schema
