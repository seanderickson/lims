from django.db import DatabaseError
from django.conf.urls import url

from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL
from tastypie import fields

from lims.api import PostgresSortingResource, CSVSerializer
from reports.models import FieldInformation, MetaHash, Vocabularies

import logging
#from django.core.serializers.json import DjangoJSONEncoder
#from django.core.serializers import json
import json
from django.conf.urls import url
        
logger = logging.getLogger(__name__)

#class FieldInformationResource(PostgresSortingResource):
#
#    class Meta:
#        queryset = FieldInformation.objects.all()
#        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
#        authorization= Authorization()        
#        # TODO: drive this from data
#        ordering = []
#        filtering = {}
#        serializer = BackboneSerializer()
#
#
#class MetaHashResource1(PostgresSortingResource):
#    json_data = fields.DictField()
#    
#    class Meta:
#        queryset = MetaHash.objects.all()
#        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
#        authorization= Authorization()        
#        # TODO: drive this from data
#        ordering = []
#        filtering = {}
#        serializer = BackboneSerializer()
#      
      
#UI_TYPE_CHOICE = 'choice'
#UI_TYPE_MULTISELECT = 'multiselect'  
        
from copy import deepcopy

class JsonAndDatabaseResource(PostgresSortingResource):
    def __init__(self, scope=None, **kwargs):
        self.scope = scope

        metahash = MetaHash.objects.get_metahash(scope=scope)
        
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)

        logger.info(str(('+++filtering', self.Meta.filtering)))
        logger.info(str(('ordering', self.Meta.ordering)))
        
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
        self.fields = deepcopy(self.original_fields)
        self.field_defs = {}
        
    def get_field_defs(self, scope):
        # dynamically define fields that are stored in the JSON field.  
        # Note that: - for a new database, there will be no fields in the JSON field initially (naturally)
        #     - also, all JSON field values will be null to start (naturally)
        # ALSO NOTE: creating new JSON fields requires a restart (to reload) (todo: fix that)
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


        return self.field_defs
    
    def build_schema(self):
        logger.info('------build_schema: ' + self.scope)
        self.reset_field_defs()
        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
        schema = super(JsonAndDatabaseResource,self).build_schema()
        
        for key, value in local_field_defs.items():
            schema['fields'][key].update(value)
            # help for fields not yet defined
            if not schema['fields'][key].get('ui_type'):
                schema['fields'][key]['ui_type'] = schema['fields'][key].get('type')
#                if schema['fields'][key].get('type') == 'string':
#                    schema['fields'][key]['ui_type'] = 'text'
#                elif schema['fields'][key].get('type') == 'integer':
#                    schema['fields'][key]['ui_type'] = 'numeric'
        
        schema['fields'].pop('json_field')
        

#        logger.info(str(('generated schema', schema)))
#        for key,value in self.field_defs.items():
#            if value.get('ui_type',None) in [UI_TYPE_CHOICE,UI_TYPE_MULTISELECT]:
#                schema['fields'][key].update({'type': value['ui_type'],
#                                              'choices': [x.key for x in Vocabularies.objects.all().filter(scope=value['vocabulary_term:scope'])]} )
#            schema['fields'][key].update({
#                  'visibility': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('visibility')) ), 
#                  'title': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('title')) ),
#                  'description': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('description')) ),
#                  'order': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('order')) )
#                  });
        
        # TODO: refactor the building of the schema here to a super class
#        schema['fields']['type'].update({'type': 'choice',
#                                         'choices': [x.key for x in Vocabularies.objects.all().filter(scope='field:type')],
#                                         'visibility': MetaHash.objects.get_or_none(scope='metahash:fields', key='type', function=lambda x : (x.get_field('visibility')) ), 
#                                         'title': MetaHash.objects.get_or_none(scope='metahash:fields', key='type', function=lambda x : (x.get_field('title')) ),
#                                         'description': MetaHash.objects.get_or_none(scope='metahash:fields', key='type', function=lambda x : (x.get_field('description')) ),
#                                        })
#        schema['fields']['visibility'].update({ 'type': 'multiselect',
#                                                'choices': [x.key for x in Vocabularies.objects.all().filter(scope='field:visibility')]})
#        schema['fields']['resource_uri'].update({ 'visibility': MetaHash.objects.get(scope='metahash:fields', key='resource_uri').get_field('visibility'), })
#        schema['fields']['json_field'].update({ 
#           'visibility': MetaHash.objects.get_or_none(scope='metahash:fields', key='json_field', function=lambda x : (x.get_field('visibility'))), })
#        logger.info(str(('schema', schema)))
        return schema

    def dehydrate(self, bundle):
#        logger.info(str(('dehydrate', bundle)))
#        bundle = super(JsonAndDatabaseResource, self).dehydrate(bundle);
        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
        for key in [ x for x,y in local_field_defs.items() if y.get('json_field_type') ]:
#            logger.info('------get field:' + key)
            bundle.data[key] = bundle.obj.get_field(key);
        
#        bundle.data['toString'] = '[' + bundle.obj.scope + ',' + bundle.obj.key +']'; # TODO: refactor this, and improve it
        bundle.data['json_field'] = ''
        bundle.data.pop('json_field') # json_field will not be part of the public API, it is for internal use only
#        logger.info(str(('deyhdrated', bundle)))
        
        # override the resource_uri, since we want to export the permanent composite key
        bundle.data['resource_uri'] = self.get_resource_uri() + bundle.data['scope'] + '/' + bundle.data['key'] +'/'
                
        # and don't send the internal id out for PATCH uses, at least
        # But: we _do_ have to send it out for Backbone, since we don't know how to use things like composite keys yet - sde4
        # bundle.data.pop('id')
        
        return bundle
    
    def hydrate_json_field(self, bundle):
#        bundle = super(JsonAndDatabaseResource, self).hydrate(bundle);
        
        logger.info(str(('hydrate', bundle)))
        
        # not sure if the bundl obj is id'd yet
        #        old_json_obj = bundle.obj.get_field_hash()
        
        # otherwise have to get it using query
        #obj = self.queryset.get(scope=bundle.data['scope'], key=bundle.data['key'])  # Note: what if trying to change the key?
        #old_json_obj = obj.get_field_hash() if obj else {}
        
        json_obj = {}
        field_reset_required = False
        local_field_defs = self.get_field_defs(self.scope) # trigger a get field defs before building the schema
        for key in [ str(x) for x,y in local_field_defs.items() if 'json_field_type' in y ]:
            json_obj.update({ key: bundle.data.get(key,None)})

#        json_obj = {    
#            'title': bundle.data['title'], 
#            'description': bundle.data['description'],
#            'comment': bundle.data['comment'],
#            'type': bundle.data['type'],
#            'visibility': bundle.data['visibility'],
#        };
        logger.info(str(('--- json obj hydrated:', json_obj)))
        bundle.data['json_field'] = json.dumps(json_obj);
        
        
        return bundle;
    
#    title = fields.CharField(attribute='title', readonly=True, blank=True, null=True)
#    description = fields.CharField(attribute='description', readonly=True, blank=True, null=True)
#    comment = fields.CharField(attribute='comment', readonly=True, blank=True, null=True)
#    type = fields.CharField(attribute='type', readonly=True, blank=True, null=True)
#    visibility = fields.ListField(attribute='visibility', readonly=True, blank=True, null=True)
    
#            # now fill in meta information for the schema report to UI; using data from this table itself, either in real or json fields!
#            for record in MetaHash.objects.all().filter(scope='metahash:fields'):
#                value.update({
#                      record.key: MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field(record.key)) )
#                      })
#
#            # now check if the field uses controlled vocabulary, look that up now.  TODO: "vocabulary_scope_ref" should be a constant
#            # TODO: "vocabulary_scope_ref" needs to be created by default as a metahash:field; this argues for making it a "real" field
#            if value.get(u'vocabulary_scope_ref'):
#                logger.info(str(('looking for a vocabulary', value['vocabulary_scope_ref'] )))
#                value['choices'] = [x.key for x in Vocabularies.objects.all().filter(scope=value['vocabulary_scope_ref'])]
#                logger.info(str(('got', value['choices'] )))
#            logger.info(str(('defined field def: ', key, value)))
#        logger.info(str(('generated fields for MetaHashResource', self.field_defs)))   
            
        # Now add in any "real" fields that aren't yet specified; this allows an empty 
#        self.field_defs = {
#                           'title': {'type': fields.CharField, 'dehydrate':'json', 'hydrate':'json'},
#                           'description': {'type': fields.CharField,'dehydrate':'json', 'hydrate':'json' },
#                           'comment': {'type': fields.CharField ,'dehydrate':'json', 'hydrate':'json'},
#                           'type': {'type': fields.CharField, 'ui_type': UI_TYPE_CHOICE, 'vocabulary_term:scope':'field:type','dehydrate':'json', 'hydrate':'json' }, 
#                           'visibility': {'type': fields.ListField, 'ui_type': UI_TYPE_MULTISELECT, 'vocabulary_term:scope':'field:visibility','dehydrate':'json', 'hydrate':'json' },
#                           'json_field': {}, # don't recreate existing fields, TODO: can we get rid of declaration of "real" fields?
#                           'resource_uri':{},
#                           'alias':{},
#                           'scope':{},
#                           'key':{},
#                           'order':{},
#                           }
#        for key,value in self.field_defs.items():
#            logger.info(str(('type', value)))
#            if(key not in self.fields): # don't redefine the existing fields
#                self.fields[key] = value['type'](attribute=key, readonly=True, blank=True, null=True)
    
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
    
#    def build_schema(self):
#        local_field_defs = self.get_field_defs() # trigger a get field defs before building the schema
#        schema = super(MetaHashResource,self).build_schema()
#        
#        for key, value in local_field_defs.items():
#            schema['fields'][key].update(value)
#        logger.info(str(('generated schema', schema)))
##        for key,value in self.field_defs.items():
##            if value.get('ui_type',None) in [UI_TYPE_CHOICE,UI_TYPE_MULTISELECT]:
##                schema['fields'][key].update({'type': value['ui_type'],
##                                              'choices': [x.key for x in Vocabularies.objects.all().filter(scope=value['vocabulary_term:scope'])]} )
##            schema['fields'][key].update({
##                  'visibility': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('visibility')) ), 
##                  'title': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('title')) ),
##                  'description': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('description')) ),
##                  'order': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('order')) )
##                  });
#        
#        # TODO: refactor the building of the schema here to a super class
##        schema['fields']['type'].update({'type': 'choice',
##                                         'choices': [x.key for x in Vocabularies.objects.all().filter(scope='field:type')],
##                                         'visibility': MetaHash.objects.get_or_none(scope='metahash:fields', key='type', function=lambda x : (x.get_field('visibility')) ), 
##                                         'title': MetaHash.objects.get_or_none(scope='metahash:fields', key='type', function=lambda x : (x.get_field('title')) ),
##                                         'description': MetaHash.objects.get_or_none(scope='metahash:fields', key='type', function=lambda x : (x.get_field('description')) ),
##                                        })
##        schema['fields']['visibility'].update({ 'type': 'multiselect',
##                                                'choices': [x.key for x in Vocabularies.objects.all().filter(scope='field:visibility')]})
##        schema['fields']['resource_uri'].update({ 'visibility': MetaHash.objects.get(scope='metahash:fields', key='resource_uri').get_field('visibility'), })
##        schema['fields']['json_field'].update({ 
##           'visibility': MetaHash.objects.get_or_none(scope='metahash:fields', key='json_field', function=lambda x : (x.get_field('visibility'))), })
##        logger.info(str(('schema', schema)))
#        return schema
    
    # called by Resource.get_list
    # def obj_get_list(self, bundle, **kwargs):
    #     calls apply_filters, which calls get_object_list
    
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
#        for key in [ x for x,y in self.field_defs.items() if 'json_field_type' in y ]:
#            bundle.data[key] = bundle.obj.get_field(key);
#        bundle.data['title'] = bundle.obj.get_field('title')
#        bundle.data['description'] = bundle.obj.get_field('description')
#        bundle.data['comment'] = bundle.obj.get_field('comment')
#        bundle.data['type'] = bundle.obj.get_field('type')
#        bundle.data['visibility'] = bundle.obj.get_field('visibility')
        
        # TODO: "toString" required by the UI dialogs dealing with this object (i.e. delete)       
#        bundle.data['toString'] = '[' + bundle.obj.scope + ',' + bundle.obj.key +']'; # TODO: refactor this, and improve it
        
        return bundle

    def hydrate(self, bundle):
        bundle = super(MetaHashResource, self).hydrate(bundle);
        return bundle
#    
#        logger.info(str(('hydrate', bundle)))
#        
#        # not sure if the bundl obj is id'd yet
#        old_json_obj = bundle.obj.get_field_hash()
#        # otherwise have to get it using query
#        #obj = self.queryset.get(scope=bundle.data['scope'], key=bundle.data['key'])  # Note: what if trying to change the key?
#        #old_json_obj = obj.get_field_hash() if obj else {}
#        
#        json_obj = {}
#        for key in [ x for x,y in self.field_defs.items() if 'json_field_type' in y ]:
#            json_obj.update({ key: bundle.data.get(key,None)})
#            if json_obj.get(key, None) != old_json_obj.get(key,None):
#                self.reset_field_defs() # force next query to rebuild the schema
#                
##        json_obj = {    
##            'title': bundle.data['title'], 
##            'description': bundle.data['description'],
##            'comment': bundle.data['comment'],
##            'type': bundle.data['type'],
##            'visibility': bundle.data['visibility'],
##        };
#        bundle.data['json_field'] = json.dumps(json_obj);
#        return bundle;
    
#    def obj_create(self, bundle, **kwargs):
#        """
#        A ORM-specific implementation of ``obj_create``.
#        """
#        bundle.obj = self._meta.object_class()
#
#        for key, value in kwargs.items():
#            setattr(bundle.obj, key, value)
#
#        self.authorized_create_detail(self.get_object_list(bundle.request), bundle)
#        bundle = self.full_hydrate(bundle)
#        return self.save(bundle)
    
#    def obj_update(self, bundle, request=None, **kwargs):
#        return self.obj_create(bundle, request, **kwargs)
#
#    def obj_delete_list(self, request=None, **kwargs):
#        bucket = self._bucket()
#
#        for key in bucket.get_keys():
#            obj = bucket.get(key)
#            obj.delete()
#
#    def obj_delete(self, request=None, **kwargs):
#        bucket = self._bucket()
#        obj = bucket.get(kwargs['pk'])
#        obj.delete()

    def rollback(self, bundles):
        pass        


class VocabulariesResource(JsonAndDatabaseResource):

#    title = fields.CharField(attribute='title', readonly=True, blank=True, null=True)
#    description = fields.CharField(attribute='description', readonly=True, blank=True, null=True)
#    comment = fields.CharField(attribute='comment', readonly=True, blank=True, null=True)
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
    
#    def build_schema(self, scope):
#        return super(VocabulariesResource,self).build_schema('metahash:fields')

    # called by Resource.get_list
    # def obj_get_list(self, bundle, **kwargs):
    #     calls apply_filters, which calls get_object_list
    
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
        
#    def hydrate(self, bundle):
##        bundle = super(MetaHashResource, self).hydrate(bundle);
#        
#        
#        logger.info(str(('hydrate', bundle)))
#        json_obj = {    
#            'title': bundle.data['title'], 
#            'description': bundle.data['description'],
#            'comment': bundle.data['comment'],
#        };
#        bundle.data['json_field'] = json.dumps(json_obj);
#        return bundle;
#
#    def dehydrate(self, bundle):
#        logger.info(str(('dehydrate', bundle)))
#        
#        bundle.data['title'] = bundle.obj.get_field('title')
#        bundle.data['description'] = bundle.obj.get_field('description')
#        bundle.data['comment'] = bundle.obj.get_field('comment')
#        
#        # TODO: "toString" required by the UI dialogs dealing with this object (i.e. delete)       
#        bundle.data['toString'] = '[' + bundle.obj.scope + ',' + bundle.obj.key +']'; # TODO: refactor this, and improve it
#        
#        return bundle
    
