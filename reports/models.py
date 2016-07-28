from __future__ import unicode_literals
import json
import logging

from django.conf import settings
from django.contrib.auth.models import User
from django.core.cache import cache
from django.db import models
from django.utils import timezone
from tastypie.utils.dict import dict_strip_unicode_keys
from django.utils.encoding import python_2_unicode_compatible
from collections import OrderedDict
from django.core.serializers.json import DjangoJSONEncoder


logger = logging.getLogger(__name__)


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

    def get_and_parse(self, scope='', field_definition_scope='fields.field', 
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
            e.g. "fields.field", or "fields.resource, or fields.vocabulary"
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
            logger.debug(
                'get_and_parse done, for %r, hash found: %r', scope, metahash.keys())
        else:
            logger.debug('retrieve the cached field definitions for %r',scope)
        return metahash


    def get_and_parse_int(self, scope='', 
                          field_definition_scope='fields.field'):
        '''
        non-cached Query the metahash table for data identified by "scope", 
        and fields defined by "field_definition_scope"
            e.g. "fields.screensaveruser", or "fields.screen"
        field_definition_scope - also defines what is in the the json_field for 
        this hash;
            e.g. "fields.field", or "fields.resource, or fields.vocabulary"
        '''
        logger.debug('get_and_parse table field definitions for scope: %r, fds: %r',
            scope, field_definition_scope)
        # try to make clear that the field definitions, though stored in the 
        # metahash as well, could be in a separate table
        field_definition_table = MetaHash.objects.all().filter(
            scope=field_definition_scope)
        if not field_definition_table:
            logger.warn('field definitions not found for: %r',
                field_definition_scope)
            return {}
        logger.debug('field_definition_table: %r', field_definition_table)
        # the objects themselves are stored in the metahash table as well
        unparsed_objects = MetaHash.objects.all().filter(scope=scope).order_by('ordinal')
        logger.debug('unparsed_objects: %r', unparsed_objects)
        parsed_objects = OrderedDict()
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
                    x.key for x in Vocabulary.objects.all().filter(
                        scope=vocab_ref)]
            
            parsed_objects[unparsed_object.key] = dict_strip_unicode_keys(parsed_object)

        return parsed_objects


API_ACTION_POST = 'POST'
API_ACTION_PUT = 'PUT'
# NOTE: "CREATE" is not a REST verb - use to distinguish PATCH/modify, PATCH/create
API_ACTION_CREATE = 'CREATE' 
API_ACTION_PATCH = 'PATCH'
API_ACTION_DELETE = 'DELETE'
API_ACTION_CHOICES = ((API_ACTION_POST,API_ACTION_POST),
                      (API_ACTION_PUT,API_ACTION_PUT),
                      (API_ACTION_CREATE,API_ACTION_CREATE),
                      (API_ACTION_PATCH,API_ACTION_PATCH),
                      (API_ACTION_DELETE,API_ACTION_DELETE))
class ApiLog(models.Model):
    
    objects = models.Manager()
    
    # FIXME: change to foreign key
    user_id = models.IntegerField(null=False)
    username = models.CharField(null=False, max_length=128)

    # name of the resource, i.e. "apilog" or "screen", "user", etc.
    ref_resource_name = models.CharField(null=False, max_length=128, db_index=True)

    # full public key of the resource instance being logged (may be composite, 
    # separted by '/')
    key = models.CharField(null=False, max_length=128, db_index=True)
    
    # the full uri of the resource instance being logged, so a combination of 
    # [base api uri]/[resource_name]/[key]
    uri = models.TextField(null=False)
    
    # date and time of the update; this is the key for the apilog record
    date_time = models.DateTimeField(null=False)
    
    api_action = models.CharField(
        max_length=10, null=False, choices=API_ACTION_CHOICES)
    
    added_keys = models.TextField(null=True)
    removed_keys = models.TextField(null=True)
    diff_keys = models.TextField(null=True)
    diffs = models.TextField(null=True)
    comment = models.TextField(null=True)
    
    parent_log = models.ForeignKey('self', related_name='child_logs', null=True)
    
    # Nested fields are defined in the json_field
    json_field = models.TextField(null=True)
    
    class Meta:
        unique_together = (('ref_resource_name', 'key', 'date_time'))    

    def __unicode__(self):
        return unicode(str((self.api_action, self.ref_resource_name, self.key, 
            str(self.date_time), self.username, self.added_keys, self.removed_keys, 
            self.diff_keys, self.diffs, self.comment, self.json_field)))

    @staticmethod   
    def json_dumps(obj):
        return json.dumps(
            obj, skipkeys=False,check_circular=True, allow_nan=True, 
            default=lambda x: str(x), cls=DjangoJSONEncoder)
    
    def save(self, **kwargs):
        ''' override to convert json fields '''
        if isinstance(self.added_keys, dict):
            self.added_keys = self.json_dumps(self.added_keys)
        
        if isinstance(self.removed_keys, dict):
            self.removed_keys = self.json_dumps(self.removed_keys)
        
        if isinstance(self.diff_keys, dict):
            self.diff_keys = self.json_dumps(self.diff_keys)
        
        if isinstance(self.diffs, dict):
            self.diffs = self.json_dumps(self.diffs)
        
        return models.Model.save(self, **kwargs)
    
    def diff_dict_to_api_log(self, log):
        '''
        Set the diff fields from a dict
        @param log a dict version of the diff fields
        '''
        if 'added_keys' in log:
            self.added_keys = self.json_dumps(log['added_keys'])
        if 'removed_keys' in log:
            self.removed_keys = self.json_dumps(log['removed_keys'])
        if 'diff_keys' in log:
            self.diff_keys = self.json_dumps(log['diff_keys'])
            self.diffs = self.json_dumps(log['diffs'])


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
    
    scope = models.CharField(max_length=64)
    key = models.CharField(max_length=64)
    alias = models.CharField(max_length=64)
    ordinal = models.IntegerField();

    # required if the record represents a JSON field; choices are from the TastyPie 
    # field types
    json_field_type = models.CharField(max_length=128, null=True); 
    
    # This is the "meta" field, it contains "virtual" json fields
    json_field = models.TextField(null=True) 

    # required if the record represents a linked  field; choices are from the TastyPie 
    # field types
    linked_field_type = models.CharField(
        max_length=128, null=True); 
    
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
    
    
class Vocabulary(models.Model):
    
    objects                 = MetaManager()
    #    objects                 = models.Manager() # default manager
    
    scope                   = models.CharField(max_length=128)
    key                     = models.CharField(max_length=128)
    alias                   = models.CharField(max_length=64)
    ordinal                 = models.IntegerField();
    title                   = models.CharField(max_length=512)
    
    # All other fields are "virtual" JSON stored fields, (unless we decide to 
    # migrate them out to real fields for rel db use cases)
    # NOTE: "json_type" for all virtual JSON fields in the entire database are
    # defined in the MetaHash
    json_field                   = models.TextField(null=True)
       
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
        return unicode(str((self.scope, self.key, self.ordinal)))

        
class Permission(models.Model):
    scope = models.CharField(max_length=64) # scope of the permission
    key = models.CharField(max_length=64)  # key of the permission
    type = models.CharField(max_length=35)

#     user_permissions = models.ManyToManyField('reports.UserProfile')
    
    
    class Meta:
        unique_together = (('scope', 'key', 'type'))    
        
    def __unicode__(self):
        return unicode(str((self.scope, self.key, self.type)))
   
    
class UserGroup(models.Model):
    name = models.TextField(unique=True)
    title = models.TextField(unique=True, null=True)
    users = models.ManyToManyField('reports.UserProfile')
    permissions = models.ManyToManyField('reports.Permission')
    
    # Super Groups: 
    # - inherit permissions from super_groups, this group is a sub_group to them
    # - inherit users from sub_groups, this group is a super_group to them
    # NOTE: we are creating an "adjacency-tree" here; this can be non-performant
    # for large trees - which is not expected here; and it also requires use
    # of postgres-specific "array" types and operators (see reports.api.UserGroupResource)
    # The trade-off here is in simplicity of maintenance.
    # see "SQL antipatterns" for more discussion.
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
        
    def __unicode__(self):
        return unicode(str((self.name)) )

# @python_2_unicode_compatible
class UserProfile(models.Model):
    objects = MetaManager()
    
    # link to django.contrib.auth.models.User, note: allow null so that it
    # can be created at the same time, but it is not allowed to be null in practice
    user = models.OneToOneField(settings.AUTH_USER_MODEL, null=True, on_delete=models.CASCADE) #, null=True) 
    # user = models.OneToOneField(User, null=True, on_delete=models.SET_NULL) #, null=True) 
    
    # will mirror the auth_user.username field
    username = models.TextField(null=False, unique=True) 
    
    # Harvard specific fields
    phone = models.TextField(null=True)
    mailing_address = models.TextField(null=True)
    comments = models.TextField(null=True)

    # TODO: make this unique
    ecommons_id = models.TextField(null=True)

    harvard_id = models.TextField(null=True)
    harvard_id_expiration_date = models.DateField(null=True)
    harvard_id_requested_expiration_date = models.DateField(null=True)
    
    created_by_username = models.TextField(null=True)

    gender = models.CharField(null=True, max_length=15)    

    # deprecated, move to auth.user
    #     first_name = models.TextField()
    #     last_name = models.TextField()
    # email = models.TextField(null=True)

    # permissions assigned directly to the user, as opposed to by group
    permissions = models.ManyToManyField('reports.Permission')

    # required if the field is a JSON field; choices are from the TastyPie 
    # field types
    json_field_type = models.CharField(
        max_length=128, null=True); 
    
    # This is the "meta" field, it contains "virtual" json fields
    json_field = models.TextField(null=True)     

    def __str__(self):
        return ('UserProfile: { ecommons_id: %r, username: %r, auth_user: %r }'
            % (self.ecommons_id, self.username, self.user )) 
    
#     def __unicode__(self):
#         return ('UserProfile: { ecommons_id: %r, username: %r, auth_user: %r }'
#             % (self.ecommons_id, self.username, self.user )) 
#     
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
    
    @property
    def email(self):
        return self.user.email  
        
    @property    
    def first_name(self):
        "Returns the person's full name."
        return self.user.first_name
    @property
    def last_name(self):
        "Returns the person's full name."
        return self.user.last_name
    


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
    scope = models.CharField(max_length=64)
    
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
    parent = models.OneToOneField('Record', unique=True)
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
    
    
