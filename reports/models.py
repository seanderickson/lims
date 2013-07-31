from django.db import models
from lims.models import GetOrNoneManager

import json
import logging

logger = logging.getLogger(__name__)

#class FieldInformation(models.Model):
#    manager                 = FieldsManager()
#    objects                 = models.Manager() # default manager
#    
#    table                   = models.CharField(max_length=35, blank=True)
#    field                   = models.CharField(max_length=35, blank=True)
#    alias                   = models.CharField(max_length=35, blank=True)
#    queryset                = models.CharField(max_length=35, blank=True)
#    field_name              = models.TextField(blank=True) # override the LINCS name for display
#    show_in_detail          = models.BooleanField(default=False, null=False) # Note: default=False are not set at the db level, only at the Db-api level
#    show_in_list            = models.BooleanField(default=False, null=False) # Note: default=False are not set at the db level, only at the Db-api level
#    order                   = models.IntegerField(null=False)
#    use_for_search_index    = models.BooleanField(default=False) # Note: default=False are not set at the db level, only at the Db-api level
#    related_to              = models.TextField(blank=True)
#    description             = models.TextField(blank=True)
#    importance              = models.TextField(blank=True)
#    comments                = models.TextField(blank=True)
#    is_unrestricted         = models.BooleanField(default=False, null=False)
#    class Meta:
#        unique_together = (('table', 'field','queryset'),('field','alias'))    
#    def __unicode__(self):
#        return unicode(str((self.table, self.field, self.id, self.field_name)))
#    
#    def get_camel_case_field_name(self):
#        logger.debug(str(('create a camelCase name for:', self)))
#        field_name = self.field_name.strip().title()
#        # TODO: convert a trailing "Id" to "ID"
#        field_name = ''.join(['ID' if x.lower()=='id' else x for x in re.split(r'[_\s]+',field_name)])
#        
#        #field_name = re.sub(r'[_\s]+','',field_name)
#        field_name = field_name[0].lower() + field_name[1:];
#        #        logger.info(str(('created camel case name', field_name, 'for', self)))
#        return field_name

class MetaManager(GetOrNoneManager):

    def __init__(self, **kwargs):
        super(MetaManager,self).__init__(**kwargs)
#        self.metahash = {}


    
    # this is how you override a Manager's base QuerySet
    def get_query_set(self):
        return super(MetaManager, self).get_query_set()

    def get_metahash(self, scope='', field_definition_scope='fields:metahash'):
        '''
        Query the metahash table for field definitions for this table
        scope - metahash scope, i.e. "fields:screensaveruser", or "fields:screen"
        field_definition_scope - scope that defines the json fields for this hash, e.g. "fields:metahash", or "fields:resource, or fields:vocabularies"
        '''
#        if not self.metahash:
        # TODO: one cannot deny, a cache might do some good here    
        metahash = {}
        # first, query the metahash for fields defined for this scope
        for fieldinformation in MetaHash.objects.all().filter(scope=scope):
            logger.debug('---- meta field for scope: ' + scope + ', ' + fieldinformation.key)
            
            field_key = fieldinformation.key
#            hash = schema['fields'][field_key]
            metahash[field_key] = {}
            hash = metahash[field_key]
            for meta_record in MetaHash.objects.all().filter(scope=field_definition_scope):  # metahash:fields are defined for all reports
                hash.update({
                      meta_record.key: MetaHash.objects.get_or_none(scope=scope, key=field_key, function=lambda x : (x.get_field(meta_record.key)) )
                      })
                    
            # now check if the field uses controlled vocabulary, look that up now.  TODO: "vocabulary_scope_ref" should be a constant
            # TODO: "vocabulary_scope_ref" needs to be created by default as a metahash:field; this argues for making it a "real" field
            if hash.get(u'vocabulary_scope_ref'):
                vocab_ref = hash['vocabulary_scope_ref']
                logger.debug(str(('looking for a vocabulary', vocab_ref )))
                hash['choices'] = [x.key for x in Vocabularies.objects.all().filter(scope=vocab_ref)]
                logger.debug(str(('got', hash['choices'] )))
            
        return metahash

#
#class MetaHash1(models.Model):
#    manager                 = MetaManager()
#    objects                 = models.Manager() # default manager
#    
#    scope                   = models.CharField(max_length=35, blank=True)
#    key                     = models.CharField(max_length=35, blank=True)
#    alias                   = models.CharField(max_length=35, blank=True)
#    
#    value                   = models.TextField(blank=True)
#    
#    description             = models.TextField(blank=True)
#    comment                 = models.TextField(blank=True)
#    labels                  = models.TextField(blank=True)
#    capabilities            = models.TextField(blank=True) # i.e. 'edit', 'textsearch'
#    visibilities            = models.TextField(blank=True) # i.e. 'detail', 'list', 'edit'
#    readPermissions         = models.TextField(blank=True) # i.e. roles: admin, etc.
#    writePermissions         = models.TextField(blank=True)
#    createPermissions         = models.TextField(blank=True)
#    type                    = models.CharField(max_length=35, blank=True)
#    
#    class Meta:
#        unique_together = (('scope', 'key'),('scope','alias'))    
#    def __unicode__(self):
#        return unicode(str((self.scope, self.key, self.id, self.alias)))
#
## Notes on MetaHash Model

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
       
    class Meta:
        unique_together = (('scope', 'key'))    
    
    def get_field_hash(self):
        if self.json_field:
            return json.loads(self.json_field)
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
        assert scope, 'Must define the scope (used to query the field definitions in the metahash)'
        fields = MetaHash.objects.get_metahash(scope=scope)
#        logger.info(str(('fields', fields)))
        dict = {}
        for key in fields.keys():
            dict[key] = self.get_field(key)
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
    
    
