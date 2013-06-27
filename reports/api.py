from django.db import connection, DatabaseError

from reports.models import FieldInformation, MetaHash, Vocabularies

from django.conf.urls import url

from tastypie.authorization import Authorization
from tastypie.bundle import Bundle
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
from tastypie import fields, utils

from db.api import PostgresSortingResource

import logging
from django.core.serializers.json import DjangoJSONEncoder
from django.core.serializers import json
import json,re
        
logger = logging.getLogger(__name__)


class BackboneSerializer(Serializer):
    
    def from_json(self, content):
        """
        Given some JSON data, returns a Python dictionary of the decoded data.
        Override to quote attributes - the backbone client doesn't want to do this.
        """
        logger.info(str(("loading content:", content)))
        content = content.replace(r'(\w+):', r'"\1" :')
        logger.info(str(("loading content:", content)))
        return json.loads(content)
    
class FieldInformationResource(PostgresSortingResource):

    class Meta:
        queryset = FieldInformation.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {}
        serializer = BackboneSerializer()


class MetaHashResource1(PostgresSortingResource):
    json_data = fields.DictField()
    
    class Meta:
        queryset = MetaHash.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {}
        serializer = BackboneSerializer()
      
      
UI_TYPE_CHOICE = 'choice'
UI_TYPE_MULTISELECT = 'multiselect'  
        
class MetaHashResource(PostgresSortingResource):

#    title = fields.CharField(attribute='title', readonly=True, blank=True, null=True)
#    description = fields.CharField(attribute='description', readonly=True, blank=True, null=True)
#    comment = fields.CharField(attribute='comment', readonly=True, blank=True, null=True)
#    type = fields.CharField(attribute='type', readonly=True, blank=True, null=True)
#    visibility = fields.ListField(attribute='visibility', readonly=True, blank=True, null=True)
    
    def __init__(self, **kwargs):
        super(MetaHashResource,self).__init__(**kwargs)
        
        self.field_defs = {
                           'title': {'type': fields.CharField, 'dehydrate':'json', 'hydrate':'json'},
                           'description': {'type': fields.CharField,'dehydrate':'json', 'hydrate':'json' },
                           'comment': {'type': fields.CharField ,'dehydrate':'json', 'hydrate':'json'},
                           'type': {'type': fields.CharField, 'ui_type': UI_TYPE_CHOICE, 'vocabulary_term:scope':'field:type','dehydrate':'json', 'hydrate':'json' }, 
                           'visibility': {'type': fields.ListField, 'ui_type': UI_TYPE_MULTISELECT, 'vocabulary_term:scope':'field:visibility','dehydrate':'json', 'hydrate':'json' },
                           'json_field': {}, # don't recreate existing fields, TODO: can we get rid of this field?
                           'resource_uri':{},
                           'alias':{},
                           'scope':{},
                           'key':{},
                           'order':{},
                           }
        for key,value in self.field_defs.items():
            logger.info(str(('type', value)))
            if(key not in self.fields): # don't redefine the existing fields
                self.fields[key] = value['type'](attribute=key, readonly=True, blank=True, null=True)
    
    class Meta:
        queryset = MetaHash.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'scope':ALL, 'key':ALL}
        serializer = BackboneSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
    
    
    def build_schema(self):
        schema = super(MetaHashResource,self).build_schema()
        
        for key,value in self.field_defs.items():
            if value.get('ui_type',None) in [UI_TYPE_CHOICE,UI_TYPE_MULTISELECT]:
                schema['fields'][key].update({'type': value['ui_type'],
                                              'choices': [x.key for x in Vocabularies.objects.all().filter(scope=value['vocabulary_term:scope'])]} )
            schema['fields'][key].update({
                  'visibility': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('visibility')) ), 
                  'title': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('title')) ),
                  'description': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('description')) ),
                  'order': MetaHash.objects.get_or_none(scope='metahash:fields', key=key, function=lambda x : (x.get_field('order')) )
                  });
        
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
    
    # called by Resource.get_list
    # def obj_get_list(self, bundle, **kwargs):
    #     calls apply_filters, which calls get_object_list
    
    # called by ModelResource in ModelResource.apply_filters 
    def get_object_list(self, request):
        obj_list = super(MetaHashResource, self).get_object_list(request);
        # could apply filtering here
        return obj_list
    
    def obj_get(self, bundle, **kwargs):
        logger.info(str(('kwargs',kwargs)))
        obj = super(MetaHashResource, self).obj_get(bundle, **kwargs)
        return obj
    
    def dehydrate(self, bundle):
        for key in [ x for x,y in self.field_defs.items() if y.get('dehydrate',None) == 'json']:
            bundle.data[key] = bundle.obj.get_field(key);
#        bundle.data['title'] = bundle.obj.get_field('title')
#        bundle.data['description'] = bundle.obj.get_field('description')
#        bundle.data['comment'] = bundle.obj.get_field('comment')
#        bundle.data['type'] = bundle.obj.get_field('type')
#        bundle.data['visibility'] = bundle.obj.get_field('visibility')
        
        # TODO: "toString" required by the UI dialogs dealing with this object (i.e. delete)       
        bundle.data['toString'] = '[' + bundle.obj.scope + ',' + bundle.obj.key +']'; # TODO: refactor this, and improve it
        
        return bundle

    def hydrate(self, bundle):
        # unnecessary to call the super
        #        bundle = super(MetaHashResource, self).hydrate(bundle);
        
        logger.info(str(('hydrate', bundle)))
        json_obj = {}
        for key in [ x for x,y in self.field_defs.items() if y.get('hydrate',None) == 'json']:
            json_obj.update({ key: bundle.data.get(key,None)})
        
#        json_obj = {    
#            'title': bundle.data['title'], 
#            'description': bundle.data['description'],
#            'comment': bundle.data['comment'],
#            'type': bundle.data['type'],
#            'visibility': bundle.data['visibility'],
#        };
        bundle.data['json_field'] = json.dumps(json_obj);
        return bundle;
    
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


class VocabulariesResource(PostgresSortingResource):

    title = fields.CharField(attribute='title', readonly=True, blank=True, null=True)
    description = fields.CharField(attribute='description', readonly=True, blank=True, null=True)
    comment = fields.CharField(attribute='comment', readonly=True, blank=True, null=True)

    class Meta:
        queryset = Vocabularies.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {'scope':ALL, 'key': ALL, 'alias':ALL}
        serializer = BackboneSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
    
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
    
    def hydrate(self, bundle):
#        bundle = super(MetaHashResource, self).hydrate(bundle);
        
        
        logger.info(str(('hydrate', bundle)))
        json_obj = {    
            'title': bundle.data['title'], 
            'description': bundle.data['description'],
            'comment': bundle.data['comment'],
        };
        bundle.data['json_field'] = json.dumps(json_obj);
        return bundle;

    def dehydrate(self, bundle):
        logger.info(str(('dehydrate', bundle)))
        
        bundle.data['title'] = bundle.obj.get_field('title')
        bundle.data['description'] = bundle.obj.get_field('description')
        bundle.data['comment'] = bundle.obj.get_field('comment')
        
        # TODO: "toString" required by the UI dialogs dealing with this object (i.e. delete)       
        bundle.data['toString'] = '[' + bundle.obj.scope + ',' + bundle.obj.key +']'; # TODO: refactor this, and improve it
        
        return bundle
    
