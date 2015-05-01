import json
import logging

from django.db import models
from django.conf import settings
from tastypie.utils.dict import dict_strip_unicode_keys
from django.utils import timezone

logger = logging.getLogger(__name__)

from django.core.cache import cache

class GetOrNoneManager(models.Manager):
    """Adds get_or_none method to objects
    """
    def get_or_none(self, function=None, **kwargs):
        try:
            x = self.get(**kwargs)
            if x and function:
                return function(x)
            else:
                return x
        except self.model.DoesNotExist: # todo: check for better err to catch
            return None

class MetaManager(GetOrNoneManager):

    def __init__(self, **kwargs):
        super(MetaManager,self).__init__(**kwargs)

    # this is how you override a Manager's base QuerySet
    def get_query_set(self):
        return super(MetaManager, self).get_query_set()

    def get_and_parse(self, scope='', field_definition_scope='fields.metahash', 
                      clear=False):
        '''
        @param scope - i.e. the table to get field definitions for
        @param field_definition_scope - i.e. where the field properties are defined
        @param clear to clear the cache

        Query the metahash table for data identified by "scope", and fields 
        defined by "field_definition_scope"
            e.g. "fields.screensaveruser", or "fields.screen"
        field_definition_scope - also defines what is in the the json_field for
        this hash;
            e.g. "fields.metahash", or "fields.resource, or fields.vocabularies"
        '''
        metahash = {}
        if not clear:
            metahash = cache.get('metahash:'+scope)
        else:
            cache.delete('metahash:'+scope)
            
        if not metahash:
            metahash = self.get_and_parse_int(
                scope=scope, field_definition_scope=field_definition_scope)
            cache.set('metahash:'+scope, metahash);
        else:
            logger.debug(str((
                'retrieve the cached field definitions for ',scope)))
        logger.debug(str((
            'get_and_parse done, for ', scope, 'hash found', metahash)))
        return metahash


    def get_and_parse_int(self, scope='', 
                          field_definition_scope='fields.metahash'):
        '''
        non-cached Query the metahash table for data identified by "scope", 
        and fields defined by "field_definition_scope"
            e.g. "fields.screensaveruser", or "fields.screen"
        field_definition_scope - also defines what is in the the json_field for 
        this hash;
            e.g. "fields.metahash", or "fields.resource, or fields.vocabularies"
        '''
        logger.debug('get_and_parse table field definitions for ' + scope)
        # try to make clear that the field definitions, though stored in the 
        # metahash as well, could be in a separate table
        field_definition_table = MetaHash.objects.all().filter(
            scope=field_definition_scope)
        
        # the objects themselves are stored in the metahash table as well
        unparsed_objects = MetaHash.objects.all().filter(scope=scope)
        
        parsed_objects = {}
        for unparsed_object in unparsed_objects:
            parsed_object = {}
            # only need the key from the field definition table
            for field_key in [ x.key for x in field_definition_table]:  
                parsed_object.update(
                    { field_key : unparsed_object.get_field(field_key) })
            
            # NOTE: choices for the "vocabulary_scope_ref" are being stored 
            # here for convenience
            # now check if the field uses controlled vocabulary, look that up now.  
            # TODO: "vocabulary_scope_ref" should be a constant
            # TODO: "vocabulary_scope_ref" needs to be created by default as a 
            # metahash:field; this argues for making it a "real" field
            
            if parsed_object.get(u'vocabulary_scope_ref'):
                vocab_ref = parsed_object['vocabulary_scope_ref']
                logger.debug(str(('vocab_ref', vocab_ref)))
                parsed_object['choices'] = [
                    x.key for x in Vocabularies.objects.all().filter(
                        scope=vocab_ref)]
            
            parsed_objects[unparsed_object.key] = dict_strip_unicode_keys(parsed_object)
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
    username = models.CharField(null=False, max_length=128)

    # name of the resource, i.e. "apilog" or "screen", "user", etc.
    ref_resource_name = models.CharField(null=False, max_length=128)

    # full public key of the resource instance being logged (may be composite, 
    # separted by '/')
    key = models.CharField(null=False, max_length=128)
    
#     #TODO How to link to parent/children?
#     ref_containing_resource = models.CharField(null=False, max_length=105)
#     containing_key = models.CharField(null=False, max_length=128)
#     contained_keys = models.ListField(...)
        
    # the full uri of the resource instance being logged, so a combination of 
    # [base api uri]/[resource_name]/[key]
    uri = models.TextField(null=False)
    
    # date and time of the update; this is the key for the apilog record
    date_time = models.DateTimeField(null=False)
    
    api_action = models.CharField(
        max_length=10, null=False, choices=API_ACTION_CHOICES)
    
    added_keys = models.TextField(blank=True, null=True)
    removed_keys = models.TextField(blank=True, null=True)
    diff_keys = models.TextField(blank=True, null=True)
    diffs = models.TextField(blank=True, null=True)
    comment = models.TextField(blank=True, null=True)
    
    parent_log = models.ForeignKey('self', related_name='child_logs', null=True)
    
    # This is the "meta" field, it contains "virtual" json fields, defined in 
    # the metahash
    json_field = models.TextField(blank=True, null=True)
    
    class Meta:
        unique_together = (('ref_resource_name', 'key', 'date_time'))    

    def __unicode__(self):
        return unicode(str((self.api_action, self.ref_resource_name, self.key, 
            str(self.date_time), self.username, self.added_keys, self.removed_keys, 
            self.diff_keys, self.diffs, self.comment, self.json_field)))
    
    def diff_dict_to_api_log(self, log):
        '''
        Set the diff fields from a dict
        @param log a dict version of the diff fields
        '''
        json_dumps = lambda y: json.dumps(
            y, skipkeys=False,check_circular=True, allow_nan=True, 
            default=lambda x: str(x))
        if 'added_keys' in log:
            self.added_keys = json_dumps(log['added_keys'])
        if 'removed_keys' in log:
            self.removed_keys = json_dumps(log['removed_keys'])
        if 'diff_keys' in log:
            self.diff_keys = json_dumps(log['diff_keys'])
            self.diffs = json_dumps(log['diffs'])


class ListLog(models.Model):
    '''
    A model that holds the keys for the items created in a "put_list"
    '''
    
    apilog = models.ForeignKey('ApiLog')

    # name of the resource, i.e. "apilog" or "screen", "user", etc.
    ref_resource_name = models.CharField(max_length=64)

    # full public key of the resource instance being logged (may be composite, 
    # separted by '/')
    key = models.CharField(max_length=128)

    uri = models.TextField()
    
    class Meta:
        # TODO: must be apilog and either of(('ref_resource_name', 'key'),'uri')
        # -- so probably better to handle in software
        unique_together = (('apilog', 'ref_resource_name', 'key','uri'))    
    
    def __unicode__(self):
        return unicode(str((self.ref_resource_name, self.key, self.uri )))
    
    
class MetaHash(models.Model):
    objects                 = MetaManager()
    #    objects                 = models.Manager() # default manager
    
    scope = models.CharField(max_length=64, blank=True)
    key = models.CharField(max_length=64, blank=True)
    alias = models.CharField(max_length=64, blank=True)
    ordinal = models.IntegerField();

    # required if the record represents a JSON field; choices are from the TastyPie 
    # field types
    json_field_type = models.CharField(
        max_length=128, blank=True, null=True); 
    
    # This is the "meta" field, it contains "virtual" json fields
    json_field = models.TextField(blank=True) 

    # required if the record represents a linked  field; choices are from the TastyPie 
    # field types
    linked_field_type = models.CharField(
        max_length=128, blank=True, null=True); 
    
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
        # return ( True if ( 
        # self.json_field_type and self.json_field_type.upper() != 'VIRTUAL' ) 
        # else False )
        return True if self.json_field_type else False
                
    def model_to_dict(self, scope=None):
        '''
        Specialized model_to_dict for JSON containing tables defined using the 
        Metahash Manager.
        - scope = "fields.<model_name>" - the scope of the field definitions in 
        the metahash table for this object.
        - the scope is used to query the "fields" definitions in the metahash - 
            the construct we are using to define all publicly available fields; 
            json or 'real'
        '''
        assert scope, (
            'Must define the scope (used to query the field definitions in the metahash)')
        fields = MetaHash.objects.get_and_parse(scope=scope)
        dict = {}
        for key in fields.keys():
            dict[key] = self.get_field(key)
        return dict
    
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.id, self.alias)))
    
    
class Vocabularies(models.Model):
    
    objects                 = MetaManager()
    #    objects                 = models.Manager() # default manager
    
    scope                   = models.CharField(max_length=128, blank=True)
    key                     = models.CharField(max_length=128, blank=True)
    alias                   = models.CharField(max_length=64, blank=True)
    ordinal                 = models.IntegerField();
    title                   = models.CharField(max_length=512, blank=True)
    
    # All other fields are "virtual" JSON stored fields, (unless we decide to 
    # migrate them out to real fields for rel db use cases)
    # NOTE: "json_type" for all virtual JSON fields in the entire database are
    # defined in the MetaHash
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
            # Note, json_field is sparse, not padded with empty attributes
            logger.debug(str((self,'field not found: ',field))) 
            return None
            
    def set_field(self, field, value):
        temp = self.get_field_hash()
        temp[field] = value;
        self.json_field = json.dumps(temp)
    
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.id)))

        
class Permission(models.Model):
    scope = models.CharField(max_length=64, blank=True) # scope of the permission
    key = models.CharField(max_length=64, blank=True)  # key of the permission
    type = models.CharField(max_length=35)

#     user_permissions = models.ManyToManyField('reports.UserProfile')
    
    
    class Meta:
        unique_together = (('scope', 'key', 'type'))    
        
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.type)))
   
    
class UserGroup(models.Model):
    name = models.TextField(unique=True, blank=False)
    users = models.ManyToManyField('reports.UserProfile')
    permissions = models.ManyToManyField('reports.Permission')
    
    # inherit permissions from super_groups, this group is a sub_group to them
    super_groups = models.ManyToManyField('self', symmetrical=False, 
        related_name='sub_groups')

    def get_all_permissions(self, sub_groups=None, **kwargs):
        '''
        get permissions of this group, and any inherited through super_groups
        @param kwargs to filter permissions by attributes
        '''
        if not sub_groups: 
            sub_groups = set()
        sub_groups.add(self) # TODO: test recursion check
        # start with this groups permissions
        permissions = set(self.permissions.filter(**kwargs));
        
        for group in self.super_groups.all():
            if group not in sub_groups:
                permissions.update(group.get_all_permissions(
                    sub_groups=sub_groups, **kwargs))
        
        return list(permissions)
    
    def get_all_users(self, super_groups=None, **kwargs):
        '''
        get users of this group, and any from sub_groups as well
        @param kwargs to filter users by attributes
        '''
        if not super_groups:
            super_groups = set()
        super_groups.add(self)
        # start with this groups users
        users = set(self.users.filter(**kwargs));
        
        for group in self.sub_groups.all():
            if group not in super_groups:
                users.update(group.get_all_users(
                    super_groups=super_groups, **kwargs))
        
        return list(users)
    
#     def get_group_members(self, **kwargs):
#         users = set(self.users.filter(**kwargs))
#         for group in self.groups.filter(sub_groups__super_group=self):
#             users.update(group.get_group_members(**kwargs));
#             
#         return users
    
    def __unicode__(self):
        return unicode(str((self.name)) )
#         return unicode(str((self.name, [x for x in self.users.all()])))
    
    
#     groups = models.ManyToManyField('self', through='reports.Membership', 
#                                            symmetrical=False, 
#                                            related_name='member_of')
# #     permission_super_groups = models.ManyToManyField('self', through='reports.Membership', 
# #                                            symmetrical=False, 
# #                                            related_name='sub_group')
# 
#     def get_group_permissions(self, **kwargs):
#         # start with this groups permissions
#         permissions = set(self.permissions.all().filter(**kwargs));
#         
#         for group in self.groups.filter(super_groups__sub_group=self):
#             permissions.update(group.get_group_permissions(**kwargs))
#         
#         return permissions
#     
#     def get_group_members(self, **kwargs):
#         # start with this groups permissions
#         members = set(self.users.all().filter(**kwargs));
#         
#         for group in self.groups.filter(sub_groups__super_group=self):
#             members.update(group.get_group_permissions(**kwargs))
#         
#         return members
#     
#     def __unicode__(self):
#         return unicode(str((self.name, self.users)))
# 
#         
# 
# class Membership(models.Model):
#     sub_group = models.ForeignKey(UserGroup, related_name='sub_groups')
#     super_group = models.ForeignKey(UserGroup, related_name='super_groups')
    
 
class UserProfile(models.Model):
    objects                 = MetaManager()
    # link to django.contrib.auth.models.User, note: allow null so that it
    # can be created at the same time, but it is not allowed to be null in practice
    user = models.OneToOneField(settings.AUTH_USER_MODEL) #, null=True, blank=True) 
    # will mirror the auth_user.username field
    username = models.TextField(null=False,blank=False) 
    
    # Harvard specific fields
    gender = models.CharField(null=True, max_length=15)
    phone = models.TextField(blank=True)
    mailing_address = models.TextField(blank=True)
    comments = models.TextField(blank=True)
    ecommons_id = models.TextField(blank=True)
    harvard_id = models.TextField(blank=True)
    harvard_id_expiration_date = models.DateField(null=True, blank=True)
    harvard_id_requested_expiration_date = models.DateField(null=True, blank=True)
    
    
    # permissions assigned directly to the user, as opposed to by group
    permissions = models.ManyToManyField('reports.Permission')

    # required if the field is a JSON field; choices are from the TastyPie 
    # field types
    json_field_type = models.CharField(
        max_length=128, blank=True, null=True); 
    
    # This is the "meta" field, it contains "virtual" json fields
    json_field = models.TextField(blank=True) 

    def __unicode__(self):
        return unicode(str((self.ecommons_id, self.username)))

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
            # Note, json_field is sparse, not padded with empty attributes
            logger.debug(str((self,'field not found: ',field))) 
            return None
            
    def set_field(self, field, value):
        temp = self.get_field_hash()
        temp[field] = value;
        self.json_field = json.dumps(temp)
        
    def _get_first_name(self):
        "Returns the person's full name."
        return self.user.first_name
    first_name = property(_get_first_name)    
    
    def _get_last_name(self):
        "Returns the person's full name."
        return self.user.last_name
    last_name = property(_get_last_name)    
    
#     def _get_full_name(self):
#         "Returns the person's full name."
#         return '%s %s' % (self.first_name, self.last_name)
#     full_name = property(_get_full_name)    
#  


# 
## proof-of-concept: Typed Record table with virtual field support:
# 
# This is a particular case of the Metahash:fields/resources tables.  
# Instead of having "virtual" fields be in the json_field, in this case they will
# stored in the child RecordValue table
## There will be one "Record" or Parent table for every node in the schema graph.
## Each RecordTable will have a RecordValue table
class Record(models.Model):
    # some fields will always be better to store on the Record table.  we'll want
    # to indicate this as well in the Metahash.  
    base_value1 = models.TextField()
    
    # the scope key points to the particular type of resource represented
    # when joining with the RecordValue table, we will get the field key we 
    # want by finding the "fields" for this scope in the Metahash:fields table
    scope = models.CharField(max_length=64, blank=True)
    
class RecordValue(models.Model):
    # name of the parent field will be stored in the meta hash
    parent = models.ForeignKey('Record')
    # this field links to the column definition
    field_meta = models.ForeignKey('Metahash')
    # name of the value field will be stored in the meta hash
    value = models.TextField(null=True)

    
class RecordMultiValue(models.Model):
    # name of the parent field will be stored in the meta hash
    parent = models.ForeignKey('Record')
    # this field links to the column definition
    field_meta = models.ForeignKey('Metahash')
    # name of the value field will be stored in the meta hash
    value = models.TextField()
    ordinal = models.IntegerField()

    class Meta:
        unique_together = (('field_meta', 'parent', 'ordinal'))    
    
class RecordValueComplex(models.Model):
    '''
    This class exists to model extant complex linked tables, i.e. SMR, RNAi
    '''
    
    # name of the parent field will be stored in the meta hash
    parent = models.ForeignKey('Record', unique=True)
    # name of the value field will be stored in the meta hash
    value1 = models.TextField(null=True)
    value2 = models.TextField(null=True)

class Job(models.Model):
    # POC: 20141128
    # Version 0.1: for processing large uploaded files
    
    # Client information
    remote_addr = models.TextField(null=True)
    request_method = models.TextField(null=True)
    
    # original path info, used by job process to submit full job to endpoint
    path_info = models.TextField(null=True)
    
    # user comment on post
    comment = models.TextField(null=True);
    
    # 0.1: original POST input is saved to this file
    input_filename = models.TextField(null=True)
    
    date_time_fullfilled = models.DateTimeField(null=True) 
    date_time_processing = models.DateTimeField(null=True) 
    date_time_requested = models.DateTimeField(null=False, default=timezone.now) 

    # if the response is large, save to a temp file
    response_filename = models.TextField(null=True)
    
    # store response to field
    response_content = models.TextField(null=True)
    
    # HTTP response code, saved from the endpoint job submission on completion
    response_code = models.IntegerField();
    
    def __unicode__(self):
        return unicode(str((self.id, self.path_info, self.date_time_requested)))
    
    
