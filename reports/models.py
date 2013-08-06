from django.db import models
from lims.models import GetOrNoneManager

import json
import logging

logger = logging.getLogger(__name__)


from django.core.cache import cache

class MetaManager(GetOrNoneManager):

    def __init__(self, **kwargs):
        super(MetaManager,self).__init__(**kwargs)
#        self.metahash = {}

    # this is how you override a Manager's base QuerySet
    def get_query_set(self):
        return super(MetaManager, self).get_query_set()

    def get_and_parse(self, scope='', field_definition_scope='fields:metahash'):
        metahash = cache.get('metahash:'+scope)
        if not metahash:
            metahash = self.get_and_parse_int(scope=scope, field_definition_scope=field_definition_scope)
            cache.set('metahash:'+scope, metahash);
        return metahash


    def get_and_parse_int(self, scope='', field_definition_scope='fields:metahash'):
        '''
        Query the metahash table for field definitions for this table
        scope - metahash scope, i.e. "fields:screensaveruser", or "fields:screen"
        field_definition_scope - scope that defines the json fields for this hash, e.g. "fields:metahash", or "fields:resource, or fields:vocabularies"
        '''
        
        # try to make clear that the field definitions, though stored in the metahash as well, could be in a separate table
        field_definition_table = MetaHash.objects.all().filter(scope=field_definition_scope)
        
        # the objects themselves are stored in the metahash table as well
        unparsed_objects = MetaHash.objects.all().filter(scope=scope)
        
        # TODO: one cannot deny, a cache might do some good here    
        parsed_objects = {}
        # first, query the metahash for fields defined for this scope
        for unparsed_object in unparsed_objects:
            logger.debug('---- meta field for scope: ' + scope + ', ' + unparsed_object.key)
            
#            hash = schema['fields'][field_key]
            parsed_object = {}
            for field_key in [ x.key for x in field_definition_table]:  # only need the key from the field definition table
                parsed_object.update({ field_key : unparsed_object.get_field(field_key) })
#                        MetaHash.objects.get_or_none(scope=scope, 
#                                                     key=field_key, 
#                                                     function=lambda x : (x.get_field(field_key)) )
                    
            # now check if the field uses controlled vocabulary, look that up now.  TODO: "vocabulary_scope_ref" should be a constant
            # TODO: "vocabulary_scope_ref" needs to be created by default as a metahash:field; this argues for making it a "real" field
            if parsed_object.get(u'vocabulary_scope_ref'):
                vocab_ref = parsed_object['vocabulary_scope_ref']
                logger.debug(str(('looking for a vocabulary', vocab_ref )))
                parsed_object['choices'] = [x.key for x in Vocabularies.objects.all().filter(scope=vocab_ref)]
                logger.debug(str(('got', parsed_object['choices'] )))
            
            parsed_objects[unparsed_object.key] = parsed_object
            
        return parsed_objects

# store field information here
class MetaHash(models.Model):
    objects                 = MetaManager()
#    objects                 = models.Manager() # default manager
    
    scope                   = models.CharField(max_length=35, blank=True)
    key                     = models.CharField(max_length=35, blank=True)
    alias                   = models.CharField(max_length=35, blank=True)
    ordinal                 = models.IntegerField();
    json_field_type         = models.CharField(max_length=128, blank=True, null=True); # required if the field is a JSON field; choices are from the TastyPie field types
    
    json_field                   = models.TextField(blank=True) # This is the "meta" field, it contains "virtual" json fields
    
    loaded_field = None
    
    class Meta:
        unique_together = (('scope', 'key'))    
    
    def get_field_hash(self):
        if self.json_field:
            if not self.loaded_field: 
                self.loaded_field = json.loads(self.json_field)
            return self.loaded_field
        else:
            return {}
    
    def get_field(self, field):
        if field in self._meta.get_all_field_names():
            return getattr(self,field)
        temp = self.get_field_hash()
        if(field in temp):
            return temp[field]
        else:
            logger.debug('unknown field: ' + field + ' for ' + str(self))
            return None
    
    def set_field(self, field, value):
        temp = self.get_field_hash()
        temp[field] = value;
        self.json_field = json.dumps(temp)
    
    def is_json(self):
        """
        Determines if this Meta record references a JSON nested field or not
        """
        return True if self.json_field_type else False
            
#    def dehydrate(self, bundle):
#        return bundle
    
    def model_to_dict(self, scope=None):
        '''
        Specialized model_to_dict for JSON containing tables defined using the Metahash Manager.
        - scope = "fields:<model_name>" - the scope of the field definitions in the metahash table for this object.
        - the scope is used to query the "fields" definitions in the metahash - the construct we are using to define all 
        publicly available fields; json or 'real'
        '''
        logger.info('------1 model_to_dict')
        assert scope, 'Must define the scope (used to query the field definitions in the metahash)'
        
        logger.info('------2 model_to_dict')
        fields = MetaHash.objects.get_and_parse(scope=scope)
#        logger.info(str(('fields', fields)))
        dict = {}
        logger.info('------3 model_to_dict')
        for key in fields.keys():
            dict[key] = self.get_field(key)
        logger.info('------4 model_to_dict')
        return dict
    
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.id, self.alias)))
    
    
class Vocabularies(models.Model):
    objects                 = MetaManager()
#    objects                 = models.Manager() # default manager
    
    scope                   = models.CharField(max_length=35, blank=True)
    key                     = models.CharField(max_length=35, blank=True)
    alias                   = models.CharField(max_length=35, blank=True)
    ordinal                 = models.IntegerField();
#    json_field_type         = models.CharField(max_length=128, blank=True, null=True); # required if the field is a JSON field; choices are from the TastyPie field types
    
    # All other fields are "virtual" JSON stored fields, (unless we decide to migrate them out to real fields for rel db use cases)
    # NOTE: "json_type" for all virtual JSON fields in the entire database are defined in the MetaHash
    json_field                   = models.TextField(blank=True)
       
    class Meta:
        unique_together = (('scope', 'key'))    
    
    def get_field_hash(self):
        if self.json_field:
            return json.loads(self.json_field)
        else:
            return {}
    
    def get_field(self, field):
        temp = self.get_field_hash()
        if(field in temp):
            return temp[field]
        else:
            logger.info(str((self,'unknown field: ',field)))
            return None
            
    def set_field(self, field, value):
        temp = self.get_field_hash()
        temp[field] = value;
        self.json_field = json.dumps(temp)
    
    
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.id)))
    
    
