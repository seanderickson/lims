from django.db import models
from lims.models import GetOrNoneManager

import json
import logging


logger = logging.getLogger(__name__)


from django.core.cache import cache

class MetaManager(GetOrNoneManager):

    def __init__(self, **kwargs):
        super(MetaManager,self).__init__(**kwargs)

    # this is how you override a Manager's base QuerySet
    def get_query_set(self):
        return super(MetaManager, self).get_query_set()

    def get_and_parse(self, scope='', field_definition_scope='fields.metahash', clear=False):
        '''
        Query the metahash table for data identified by "scope", and fields defined by "field_definition_scope"
            e.g. "fields.screensaveruser", or "fields.screen"
        field_definition_scope - also defines what is in the the json_field for this hash;
            e.g. "fields.metahash", or "fields.resource, or fields.vocabularies"
        '''
        metahash = {}
        if not clear:
            metahash = cache.get('metahash:'+scope)
        else:
            cache.delete('metahash:'+scope)
            
        if not metahash:
            metahash = self.get_and_parse_int(scope=scope, field_definition_scope=field_definition_scope)
            cache.set('metahash:'+scope, metahash);
        else:
            logger.debug(str(('retrieve the cached get_and_parse metahash table field definitions for ',scope)))
        logger.debug(str(('get_and_parse done, for ', scope, 'hash found', metahash)))
        return metahash


    def get_and_parse_int(self, scope='', field_definition_scope='fields.metahash'):
        '''
        non-cached Query the metahash table for data identified by "scope", and fields defined by "field_definition_scope"
            e.g. "fields.screensaveruser", or "fields.screen"
        field_definition_scope - also defines what is in the the json_field for this hash;
            e.g. "fields.metahash", or "fields.resource, or fields.vocabularies"
        '''
        logger.debug('get_and_parse table field definitions for ' + scope)
        # try to make clear that the field definitions, though stored in the metahash as well, could be in a separate table
        field_definition_table = MetaHash.objects.all().filter(scope=field_definition_scope)
        
        # the objects themselves are stored in the metahash table as well
        unparsed_objects = MetaHash.objects.all().filter(scope=scope)
        
        parsed_objects = {}
        for unparsed_object in unparsed_objects:
            logger.debug('---- meta field for scope: ' + scope + ', ' + unparsed_object.key)
            parsed_object = {}
            for field_key in [ x.key for x in field_definition_table]:  # only need the key from the field definition table
                parsed_object.update({ field_key : unparsed_object.get_field(field_key) })
            
            # NOTE: choices for the "vocabulary_scope_ref" are being stored here for convenience
            # now check if the field uses controlled vocabulary, look that up now.  TODO: "vocabulary_scope_ref" should be a constant
            # TODO: "vocabulary_scope_ref" needs to be created by default as a metahash:field; this argues for making it a "real" field
            if parsed_object.get(u'vocabulary_scope_ref'):
                vocab_ref = parsed_object['vocabulary_scope_ref']
                parsed_object['choices'] = [x.key for x in Vocabularies.objects.all().filter(scope=vocab_ref)]
            
            parsed_objects[unparsed_object.key] = parsed_object
        #logger.debug(str(('got metahash, for ', scope, 'hash found', parsed_objects)))
        return parsed_objects


API_ACTION_POST = 'POST'
API_ACTION_PUT = 'PUT'
API_ACTION_PATCH = 'PATCH'
API_ACTION_DELETE = 'DELETE'
API_ACTION_CHOICES = ((API_ACTION_POST,API_ACTION_POST),
                      (API_ACTION_PUT,API_ACTION_PUT),
                      (API_ACTION_PATCH,API_ACTION_PATCH),
                      (API_ACTION_DELETE,API_ACTION_DELETE))
class ApiLog(models.Model):
    objects = models.Manager()
    
    user_id = models.IntegerField(null=False, blank=False)
    username = models.CharField(null=False, max_length=35)
    # name of the resource, i.e. "apilog" or "screen", "user", etc.
    ref_resource_name = models.CharField(null=False, max_length=35)
    # full public key of the resource instance being logged (may be composite, separted by '/')
    key = models.CharField(null=False, max_length=128)
    # the full uri of the resource instance being logged, so a combination of [base api uri]/[resource_name]/[key]
    uri = models.TextField(null=False)
    # date and time of the update; this is the key for the apilog record
    date_time = models.DateTimeField(null=False)
    api_action = models.CharField(max_length=10, null=False, choices=API_ACTION_CHOICES)
    
    added_keys = models.TextField(blank=True)
    removed_keys = models.TextField(blank=True)
    diff_keys = models.TextField(blank=True)
    diffs = models.TextField(blank=True)
    
    
    json_field = models.TextField(blank=True) # This is the "meta" field, it contains "virtual" json fields, defined in the metahash

    
#    @property
#    def get_date_time(self):
#        return time(self.milliseconds)

    class Meta:
        unique_together = (('ref_resource_name', 'date_time'))    

    def __unicode__(self):
        return unicode(str((self.api_action, self.ref_resource_name, self.key, str(self.date_time), self.username, self.diffs)))
    


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
        logger.info(str(('===================setting value', temp, field, value, json.dumps(temp))))
        self.json_field = json.dumps(temp)
    
    def is_json(self):
        """
        Determines if this Meta record references a JSON nested field or not
        """
#        return True if ( self.json_field_type and self.json_field_type.upper() != 'VIRTUAL' ) else False
        return True if self.json_field_type else False
            
#    def dehydrate(self, bundle):
#        return bundle
    
    def model_to_dict(self, scope=None):
        '''
        Specialized model_to_dict for JSON containing tables defined using the Metahash Manager.
        - scope = "fields.<model_name>" - the scope of the field definitions in the metahash table for this object.
        - the scope is used to query the "fields" definitions in the metahash - the construct we are using to define all 
        publicly available fields; json or 'real'
        '''
        assert scope, 'Must define the scope (used to query the field definitions in the metahash)'
        fields = MetaHash.objects.get_and_parse(scope=scope)
#        logger.info(str(('create dict: scope for field lookup: ', scope, 'fields found:', fields)))
        dict = {}
        for key in fields.keys():
#            logger.info('key: ' + key)
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
            logger.debug(str((self,'field not found: ',field))) # Note, json_field is sparse, not padded with empty attributes
            return None
            
    def set_field(self, field, value):
        temp = self.get_field_hash()
        temp[field] = value;
        self.json_field = json.dumps(temp)
    
    
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.id)))

        
class Permission(models.Model):
    scope                   = models.CharField(max_length=35, blank=True)  # scope of the permission
    key                     = models.CharField(max_length=35, blank=True)  # key of the permission
    type = models.CharField(max_length=15)
    class Meta:
        unique_together = (('scope', 'key', 'type'))    
        
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.type)))

    
    
class UserGroup(models.Model):
    name = models.TextField(unique=True, blank=False)
#    users = models.ManyToManyField('auth.User')
    users = models.ManyToManyField('db.ScreensaverUser')
    permissions = models.ManyToManyField('reports.Permission')

    def __unicode__(self):
        return unicode(str((self.name, self.users)))
        
    
#    
#class User(models.Model):
#    username = models.TextField(unique=True, blank=False)
#
##    permissions = models.ManyToManyField('Permission')
#    # permissions are simply a reference to a meta concept
#    readPermissions = models.ManyToManyField('Metahash', related_name='read_user')
#    writePermissions = models.ManyToManyField('Metahash', related_name='write_user')
#
#    def __unicode__(self):
#        return unicode(str((self.username)))

#class Permission(models.Model):
#    name = models.TextField(unique=True, blank=False)
#    
#    groups = models.ManyToManyField('Group')
#    users = models.ManyToManyField('User')
#    
#    type = models.TextField();
#    
#    resource = models.ManyToManyField('Metahash')
#
#    def __unicode__(self):
#        return unicode(str((self.name)))

    
#class UserPermission(models.Model):
#    user = models.ForeignKey('User')
#    permission = models.ForeignKey('Permission')
#    
#class UserGroup(models.Model):
#    user = models.ForeignKey('User')
#    group = models.ForeignKey('Group')
    
    
    
