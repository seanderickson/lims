from __future__ import unicode_literals

from _collections import defaultdict
from collections import OrderedDict
import json
import logging

from aldjemy.core import get_engine
from django.conf import settings
from django.contrib.auth.models import User
from django.core.cache import cache
from django.core.serializers.json import DjangoJSONEncoder
from django.db import models
from django.db import transaction
from django.utils import timezone
from django.utils.encoding import smart_bytes
import six

from reports import strftime_log
import reports.schema as SCHEMA
from reports.serialize import LimsJSONEncoder


API_ACTION = SCHEMA.VOCAB.apilog.api_action

logger = logging.getLogger(__name__)

def dict_strip_unicode_keys(uni_dict):
    """
    Convert a dict of unicode keys into a dict of ascii keys.

    Useful for converting a dict to a kwarg-able format.
    
    # see tastypie.dict_strip_unicode_keys for original implementation.
    """
    if six.PY3:
        return uni_dict

    data = {}

    for key, value in uni_dict.items():
        data[smart_bytes(key)] = value

    return data


class LogDiff(models.Model):
    
    # reference to the parent ApiLog
    log = models.ForeignKey('ApiLog', on_delete=models.CASCADE)
    
    #field = models.ForeignKey('Metahash')
    field_key = models.TextField()
    field_scope = models.TextField()
    
    before = models.TextField(null=True)
    after = models.TextField(null=True)
    
    class Meta:
        unique_together = (('log','field_key','field_scope'))    

    def __repr__(self):
        return (
            "<LogDiff(ref_resource_name=%r, key=%r, field=%r)>" 
            % (self.log.ref_resource_name, self.log.key, self.field_key))

class ApiLog(models.Model):
    
    objects = models.Manager()
    
    # FIXME: change to foreign key
    user_id = models.IntegerField(null=False)
    username = models.CharField(null=False, max_length=128)

    # Name of the resource, i.e. "apilog" or "screen", "user", etc.
    ref_resource_name = models.CharField(
        null=False, max_length=128, db_index=True)

    # Full public key of the resource instance being logged (may be composite, 
    # separated by '/')
    key = models.CharField(null=False, max_length=128, db_index=True)
    
    # Full uri of the resource instance being logged, 
    # a combination of [base api uri]/[resource_name]/[key]
    uri = models.TextField(null=False)
    
    # Date and time of the update; this is the key for the apilog record
    date_time = models.DateTimeField(null=False)
    
    api_action = models.CharField(
        max_length=10, null=False, 
        choices=zip(
            API_ACTION.get_ordered_dict().keys(),
            API_ACTION.get_ordered_dict().values(),
            )
        )
    
    comment = models.TextField(null=True)
    
    parent_log = models.ForeignKey(
        'self', related_name='child_logs', null=True, on_delete=models.CASCADE)
    
    # json_field stores meta information
    json_field = models.TextField(null=True)
    
    is_preview = models.BooleanField(default=False)
    
    class Meta:
        unique_together = (('ref_resource_name', 'key', 'date_time'))    

    def __init__(self, *args, **kwargs):
        self.diffs = {}
        models.Model.__init__(self, *args, **kwargs)
    
    def __repr__(self):
        return (
            '<ApiLog(id=%r, api_action=%r, ref_resource_name=%r, '
            'key=%r, uri=%r, date_time=%s, su_id: %r, username=%s)>'
            % (self.id, self.api_action, self.ref_resource_name, self.key,
               self.uri, strftime_log(self.date_time), self.user_id, 
               self.username))

    @property
    def log_uri(self):
        ''' Return the URI of the ApiLog
        '''
        
        return '/'.join([
            self.ref_resource_name,self.key, strftime_log(self.date_time)])
    
    @staticmethod   
    def json_dumps(obj):
        
        obj_as_dict = { k:v for k,v in vars(obj).items() if k[0] != '_'}
        diffs = defaultdict(list)
        for dl in obj.logdiff_set.all():
            diffs[dl.field_key] = [dl.before,dl.after]
        obj_as_dict['diffs'] = dict(diffs)
        return json.dumps(obj_as_dict, cls=LimsJSONEncoder)
    
    @staticmethod
    def _encode_before_after(val):
        '''
        Encode LogDiff.before and after values:
        All diff values are stored as strings unless they represent list values:
        - in this case the JSON representation of the list is stored. 
        '''
        if val is None:
            return val
        if isinstance(val, (list,tuple)):
            val = json.dumps(val, cls=LimsJSONEncoder)
        elif not isinstance(val, six.string_types):
            val = str(val)
        return val
        
    def save(self, **kwargs):
        ''' 
        Override to store/encode the diffs:
        - before/after are stored as the string representation of the field 
        value.
        - if the field value is a list, store the JSON representation.
        '''
        
        is_new = self.pk is None
        
        logger.debug('encode json field: log: %r', self.json_field)
        if self.json_field:
            if isinstance(self.json_field, dict):
                try:
                    self.json_field = json.dumps(self.json_field, cls=LimsJSONEncoder)
                except:
                    logger.exception('error with json_field value encoding: %r - %r', 
                        e, json_field)
        models.Model.save(self, **kwargs)
        
        if is_new:
            logger.debug('logging new diffs: %r', self.diffs)
            bulk_create_diffs = []
            for key,diffs in self.diffs.items():
                assert isinstance(diffs, (list,tuple))
                assert len(diffs) == 2
                before = self._encode_before_after(diffs[0])
                after = self._encode_before_after(diffs[1])
                bulk_create_diffs.append(LogDiff(
                    log=self,
                    field_key = key,
                    field_scope = 'fields.%s' % self.ref_resource_name,
                    before=before,
                    after=after))
            LogDiff.objects.bulk_create(bulk_create_diffs)
        else:
            logger.debug('logging update diffs: %r', self.diffs)
            # Note: this option should not be used for bulk creation
            for key,diffs in self.diffs.items():
                assert isinstance(diffs, (list,tuple))
                assert len(diffs) == 2
                before = self._encode_before_after(diffs[0])
                after = self._encode_before_after(diffs[1])
                found = False
                for logdiff in self.logdiff_set.all():
                    if logdiff.field_key == key:
                        logdiff.before=before
                        logdiff.after=after
                        logdiff.save()
                        found = True
                if not found:
                    LogDiff.objects.create(
                        log=self,
                        field_key = key,
                        field_scope = 'fields.%s' % self.ref_resource_name,
                        before=before,
                        after=after)
                        
    @staticmethod
    def bulk_create(logs):
        '''
        Utility method - bulk create/save ApiLog instances
        '''

        logger.debug('bulk create logs: %r', logs)
        with transaction.atomic():
            with get_engine().connect() as conn:
                last_id = int(
                    conn.execute(
                        'select last_value from reports_apilog_id_seq;')
                        .scalar() or 0)
                
                for log in logs:
                    if log.json_field:
                        if isinstance(log.json_field, dict):
                            try:
                                log.json_field = json.dumps(log.json_field, cls=LimsJSONEncoder)
                            except:
                                logger.exception('error with json_field value encoding: %r - %r', 
                                    e, log.json_field)
                
                ApiLog.objects.bulk_create(logs)
                #NOTE: Before postgresql & django 1.10 only: 
                # ids must be manually created on bulk create
                for i,log in enumerate(logs):
                    log.id = last_id+i+1
            
                bulk_create_diffs = []
                for i,log in enumerate(logs):
                    for key, logdiffs in log.diffs.items():
                        bulk_create_diffs.append(
                            LogDiff(
                                log=log,
                                field_key = key,
                                field_scope = 'fields.%s' % log.ref_resource_name,
                                before=ApiLog._encode_before_after(logdiffs[0]),
                                after=ApiLog._encode_before_after(logdiffs[1]))
                        )
                LogDiff.objects.bulk_create(bulk_create_diffs)
            
            return logs


class MetaManager(models.Manager):
    ''' Manager for the MetaHash table '''

    def __init__(self, **kwargs):
        super(MetaManager,self).__init__(**kwargs)

    def get_or_none(self, function=None, **kwargs):
        try:
            x = self.get(**kwargs)
            if x and function:
                return function(x)
            else:
                return x
        except self.model.DoesNotExist: 
            return None

    def get_and_parse(self, scope='', field_definition_scope='fields.field', 
                      clear=False):
        '''
        Parse fields from the Metadata store using definitions specified by the
        "field_definition_scope" (cached version).

        @param scope defines the target resource to return fields for:
            e.g. "fields.screensaveruser", or "fields.screen"
        
        @param field_definition_scope defines the schema of the field objects 
            that will be parsed:
            e.g. "fields.field", or "fields.resource, or fields.vocabulary"
        @param clear to clear the cache

        '''
        metahash = {}
        if not clear:
            metahash = cache.get('metahash:%s'%scope)
        else:
            cache.delete('metahash:%s'%scope)
            
        if not metahash:
            metahash = self._get_and_parse(
                scope=scope, field_definition_scope=field_definition_scope)
            cache.set('metahash:'+scope, metahash);
            logger.debug(
                'get_and_parse done, for %r, hash found: %r', 
                scope, metahash.keys())
        else:
            logger.debug('retrieve the cached field definitions for %r',scope)
        return metahash


    def _get_and_parse(self, scope='', 
                          field_definition_scope='fields.field'):
        '''
        Parse fields from the Metadata store using definitions specified by the
        "field_definition_scope" (non-cached version).
        '''
        logger.debug(
            'get_and_parse table field definitions for scope: %r, fds: %r',
            scope, field_definition_scope)
        
        # Parse field schema definition first
        field_definition_table = MetaHash.objects.all().filter(
            scope=field_definition_scope)
        if not field_definition_table:
            # Bootstrap case
            logger.warn('field definitions not found for: %r',
                field_definition_scope)
            return {}
        logger.debug('field_definition_table: %r', 
            [field.key for field in field_definition_table])
        
        # Use the field schema definition to parse the row objects:
        # Row objects themselves are stored in the metahash table as well.
        unparsed_objects = \
            MetaHash.objects.all().filter(scope=scope).order_by('ordinal')
        logger.debug('unparsed_objects: %r', 
            [field.key for field in unparsed_objects])
        
        parsed_objects = OrderedDict()
        for unparsed_object in unparsed_objects:

            parsed_object = {}
            
            for field_key in [x.key for x in field_definition_table]:
                parsed_object[field_key] = unparsed_object.get_field(field_key)
                
            # NOTE: choices for the "vocabulary_scope_ref" are being stored 
            # here for convenience
            # TODO: restrict choices to retired != True
            if parsed_object.get(u'vocabulary_scope_ref'):
                vocab_ref = parsed_object['vocabulary_scope_ref']
                parsed_object['choices'] = [
                    x.key for x in Vocabulary.objects.all().filter(
                        scope=vocab_ref)]
            
            parsed_objects[unparsed_object.key] = \
                dict_strip_unicode_keys(parsed_object)

        return parsed_objects


class MetaHash(models.Model):
    '''
    Store API metadata for resources and fields.
    
    Uses a composite key consisting of (scope, key).

    Additional data are defined as JSON in the json_field and are parsed 
    using the schema defined via the "field" resource.
    '''
    
    objects = MetaManager()
    
    scope = models.CharField(max_length=64)
    key = models.CharField(max_length=64)
    alias = models.CharField(max_length=64)
    ordinal = models.IntegerField();

    # Required if the record represents a JSON field; 
    # choices are from the TastyPie field types
    json_field_type = models.CharField(max_length=128, null=True); 
    
    json_field = models.TextField(null=True) 

    # Required if the record represents a linked field; 
    # choices are from the TastyPie field types.
    # 20181214 - deprecate, used for Reagent resource subtypes.
    linked_field_type = models.CharField(
        max_length=128, null=True); 
    
    loaded_field = None
    
    class Meta:
        unique_together = (('scope', 'key'))    
    
    def get_json_field_hash(self):
        ''' Load the json field.'''
        if self.json_field:
            if not self.loaded_field: 
                self.loaded_field = json.loads(self.json_field)
            return self.loaded_field
        else:
            return {}
    
    def get_field(self, field):
        '''Retrieve a field value, parsing it from the json_field if needed.'''
        field_names = set([
            f.name for f in self._meta.get_fields()])
        if field in field_names:
            return getattr(self,field)
        temp = self.get_json_field_hash()
        if field in temp:
            return temp[field]
        else:
            logger.debug('unknown field: %s',field)
            return None
    
    def set_json_field(self, field, value):
        '''Set a value to a field, stored internally in the json_field.'''
        
        temp = self.get_json_field_hash()
        temp[field] = value;
        self.json_field = json.dumps(temp, cls=LimsJSONEncoder)
    
    def is_json(self):
        '''True if the current metadata represents a JSON nested field'''
        return True if self.json_field_type else False
                
    def model_to_dict(self, scope):
        '''
        Specialized model_to_dict that will parse data from the internal
        json_field.
        
        @param scope for finding field definitions in the metahash table.
        e.g. "fields.<model_name>" 
        '''
        fields = MetaHash.objects.get_and_parse(scope=scope)
        dict = {}
        for key in fields.keys():
            dict[key] = self.get_field(key)
        return dict
    
    def __repr__(self):
        return (
            '<MetaHash(id=%r, scope=%r, key=%r)>'
            % (self.id, self.scope, self.key))

class Vocabulary(models.Model):
    '''Store vocabularies for the API'''
    
    objects                 = MetaManager()
    
    scope                   = models.CharField(max_length=128)
    key                     = models.CharField(max_length=128)
    alias                   = models.CharField(max_length=64)
    ordinal                 = models.IntegerField();
    title                   = models.CharField(max_length=512)
    is_retired              = models.NullBooleanField()
    comment                 = models.TextField(null=True)
    description             = models.TextField(null=True)
    
    # Other data may be stored in the json_field; field metadata for 
    # the vocabulary resource may be used to define these fields.
    json_field                   = models.TextField(null=True)
       
    class Meta:
        unique_together = (('scope', 'key'))    
    
    def get_json_field_hash(self):
        if self.json_field:
            return json.loads(self.json_field)
        else:
            return {}
    
    def get_field(self, field):
        temp = self.get_json_field_hash()
        if(field in temp):
            return temp[field]
        else:
            # Note, json_field is sparse, not padded with empty attributes
            logger.debug('%r, field not found: %r',self, field) 
            return None
            
    def set_json_field(self, field, value):
        temp = self.get_json_field_hash()
        temp[field] = value;
        self.json_field = json.dumps(temp, cls=LimsJSONEncoder)
    
    def __repr__(self):
        return (
            '<Vocabulary(scope=%r, key=%r, ordinal=%r)>'
            % (self.scope, self.key, self.ordinal))

        
class Permission(models.Model):
    '''
    Permissions for API resources and fields:

    - scope of the resource or field,
    - key of the resource or field,
    - type, e.g. "read" or "write".
    '''

    scope = models.CharField(max_length=64) # scope of the permission
    key = models.CharField(max_length=64)  # key of the permission
    type = models.CharField(max_length=35)
    
    class Meta:
        unique_together = (('scope', 'key', 'type'))    
        
    def __repr__(self):
        return (
            '<Permission(scope=%r, key=%r, type=%r)>'
            % (self.scope, self.key, self.type))
   
    
class UserGroup(models.Model):
    '''
    Define UserGroups that will have permissions:
    - UserGroups contain users,
    - UserGroups may belong to other UserGroups through the super_groups,
        - A UserGroup inherits the super_group permissions,
        - A super_group inherits the UserGroup's users.
    '''
    
    name = models.TextField(unique=True)
    title = models.TextField(unique=True, null=True)
    users = models.ManyToManyField('reports.UserProfile')
    permissions = models.ManyToManyField('reports.Permission')
    description = models.TextField(null=True)
    
    # Super Groups: 
    # - inherit permissions from super_groups, this group is a sub_group to them
    # - inherit users from sub_groups, this group is a super_group to them
    # NOTE: Creates an "adjacency-tree" here; this can be non-performant
    # for large trees - which is not expected here; and it also requires use
    # of postgres-specific "array" types and operators 
    # (see reports.api.UserGroupResource).
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
        
    def __repr__(self):
        return (
            '<UserGroup(name=%r)>'
            % (self.name))

class UserProfile(models.Model):
    ''' Extend the Django auth.User for the reports API.'''

    objects = MetaManager()
    
    # Link to django.contrib.auth.models.User, note: allow null so that it
    # can be created at the same time, but not null in practice
    user = models.OneToOneField(
        settings.AUTH_USER_MODEL, null=True, on_delete=models.CASCADE) 
    
    # Note: username mirrors the auth_user.username field
    username = models.TextField(null=False, unique=True) 
    
    # Reports API specific fields
    phone = models.TextField(null=True)
    mailing_address = models.TextField(null=True)
    comments = models.TextField(null=True)

    ecommons_id = models.TextField(null=True)

    # Note: fields also found on ScreensaverUser
    harvard_id = models.TextField(null=True)
    harvard_id_expiration_date = models.DateField(null=True)
    harvard_id_requested_expiration_date = models.DateField(null=True)
    
    created_by_username = models.TextField(null=True)

    gender = models.CharField(null=True, max_length=15)    

    # These permissions assigned directly to the user, as opposed to by group
    permissions = models.ManyToManyField('reports.Permission')

    def __repr__(self):
        return (
            '<UserProfile(username=%r, id=%d, auth_user=%d)>'
            % (self.username, self.id, self.user.id))
    
    def get_all_groups(self):

        groups = set()
        for group in self.usergroup_set.all():
            groups.add(group)
            for supergroup in group.super_groups.all():
                if supergroup in groups:
                    continue
                groups.add(supergroup)
            
        return groups
    
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
    

class Job(models.Model):
    '''Manage API background job data'''

    user_requesting = models.ForeignKey(
        'UserProfile', on_delete=models.PROTECT)
    
    # Unique URI for the resource action being serviced
    uri = models.TextField()
    method = models.TextField()
    encoding = models.TextField()
    content_type = models.TextField()
    http_accept = models.TextField()
    # JSON encoded request params
    params = models.TextField()
    
    # user comment on post
    comment = models.TextField(null=True);
    # Extra posted context data (filenames, etc.); JSON encoded
    context_data = models.TextField(null=True)
      
    # Assigned when the job is running
    process_id = models.TextField(null=True)
    # Extra runtime information, json encoded
    process_env = models.TextField(null=True)
    process_messages = models.TextField(null=True)
    
    state = models.TextField(
        default=SCHEMA.VOCAB.job.state.PENDING, 
        choices=zip(
            SCHEMA.VOCAB.job.state.get_ordered_dict().keys(),
            SCHEMA.VOCAB.job.state.get_ordered_dict().values(),
            ))
    date_time_requested = models.DateTimeField(null=False, default=timezone.now) 
    date_time_submitted = models.DateTimeField(null=True) 
    date_time_processing = models.DateTimeField(null=True) 
    date_time_completed = models.DateTimeField(null=True) 
    
    response_status_code = models.IntegerField(null=True)
    #JSON encoded response content
    response_content = models.TextField(null=True)
      
    def __repr__(self):
        return (
            '<Job(id=%d, user_id=%d)>'
            % (self.id, self.user_requesting.id))
     
