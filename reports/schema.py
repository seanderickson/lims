from __future__ import unicode_literals 
''' Schema utils:

Define schema constants
- Field names used by the API and intended as constant over versions
- Vocabulary values used by the API and intended as constant over versions
- Provide a dependency free reference to API fields

Utilities for parsing the schema
 '''

import json
import logging
from django.conf import settings
from collections import OrderedDict

logger = logging.getLogger(__name__)


VERSION = 'v1'
REPORTS_API_URI = 'reports/api/%s' % VERSION

API_RESULT_OBJ = 'object'
API_RESULT_DATA = 'objects'
API_RESULT_META = 'meta'

API_MSG_SUBMIT_COUNT = 'Data submitted'
API_MSG_RESULT = 'Result'
API_MSG_WARNING = 'Warning'
API_MSG_NOT_ALLOWED = 'Action not allowed'
API_MSG_UPDATED = 'Updated'
API_MSG_CREATED = 'Created'
API_MSG_UNCHANGED = 'Unchanged'
API_MSG_COMMENTS = 'Comments'
API_MSG_ACTION = 'Action'
API_MSG_SUCCESS = 'Success'
API_MSG_RESTRICTED_DATA = '(restricted data)'

API_PARAM_SEARCH = 'search'
# Nested search data; not encoded, a hash of data being passed internally 
# If passed from the client, nested search data will be ANDed with other searches
API_PARAM_NESTED_SEARCH = 'nested_search_data'
# Complex search - a search data structure sent as a POST header 
API_PARAM_COMPLEX_SEARCH_ID = 'search_id'

URI_PATH_COMPLEX_SEARCH = 'csearch'

# Date format for API - time zone is not used for dates
DATE_FORMAT = "%Y-%m-%d"
# Date Time format to use for serialized date times
DATE_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"
DATE_TIME_FILE_FORMAT = "%Y%m%d-%H%M%S"

class schema_obj(object):
    @classmethod
    def get_dict(cls):
        return { 
            k:v for k,v in vars(cls).items() 
                if k[0].isupper() }
    @classmethod    
    def get_ordered_dict(cls):
        return OrderedDict((
            (k,getattr(cls,k)) for k in dir(cls) 
                if k[0].isupper() and callable(getattr(cls, k)) is not True ))
    
    @classmethod
    def get_lookup_dict(cls):
        return { v:k for k,v in cls.get_dict().items() }

class ERROR(schema_obj):
    # Note: there is no ErrorResource
    resource_name = 'errors'
    
    LINE = 'line'

class API_PARAM(schema_obj):
    ''' API resource request parameters.'''
    
    VISIBILITIES = 'visibilities'

   
##### Define API resources

class RESOURCE(schema_obj):
    resource_name = 'resource'
    
    KEY = 'key'
    SCOPE = 'scope'
    ORDINAL = 'ordinal'
    TITLE = 'title'
    DESCRIPTION = 'description'
    COMMENT = 'comment'
    RESOURCE_URI = 'resource_uri'
    ID_ATTRIBUTE = 'id_attribute'
    TITLE_ATTRIBUTE = 'title_attribute'
    
    FIELDS = 'fields'
    
class FIELD(schema_obj):
    resource_name = 'field'
    
    KEY = 'key'
    SCOPE = 'scope'
    ORDINAL = 'ordinal'
    TITLE = 'title'
    DESCRIPTION = 'description'
    COMMENT = 'comment'
    DATA_TYPE = 'data_type'
    DISPLAY_TYPE = 'display_type'
    DISPLAY_OPTIONS = 'display_options'
    EDIT_TYPE = 'edit_type'
    VISIBILITY = 'visibility'
    EDITABILITY = 'editability'
    ORDERING = 'ordering'
    FILTERING = 'filtering'
    VOCAB_SCOPE_REF = 'vocabulary_scope_ref'
    RESOURCE_URI = 'resource_uri'
    
    DEPENDENCIES = 'dependencies'
    VALUE_TEMPLATE = 'value_template'
    IS_ALPHANUMERIC = 'alphanumeric_sort'
    
    DATA_ACCESS_LEVEL = 'data_access_level'
    VIEW_GROUPS = 'view_groups' 

class VOCABULARY(schema_obj):
    resource_name = 'vocabulary'
    
    TITLE = 'title'
    KEY = 'key'
    SCOPE = 'scope'
    ORDINAL = 'ordinal'
    TITLE = 'title'
    COMMENT = 'comment'
    DESCRIPTION = 'description'
    IS_RETIRED = 'is_retired'
    #'expire_interval_days'
    #'expire_notifiy_days'  

class USER(schema_obj):
    resource_name = 'user'
    
    USERNAME = 'username'
    FIRST_NAME = 'first_name'
    LAST_NAME = 'last_name'
    EMAIL = 'email'
    PERMISSIONS = 'permissions'
    USERGROUPS = 'usergroups'
    IS_ACTIVE = 'is_active'
    IS_STAFF = 'is_staff'
    IS_SUPERUSER = 'is_superuser'
    
class USERGROUP(schema_obj):
    resource_name = 'usergroup'
    
    NAME = 'name'
    DESCRIPTION = 'description'
    USERS = 'users'
    SUPER_GROUPS = 'super_groups'
    SUB_GROUPS = 'sub_groups'
    PERMISSIONS = 'permissions'
    ALL_SUB_GROUPS = 'all_sub_groups'
    ALL_SUPER_GROUPS = 'all_super_groups'
    ALL_USERS = 'all_users'
    ALL_PERMISSIONS = 'all_permissions'
    
    
class PERMISSION(schema_obj):
    resource_name = 'permission'
    
    SCOPE = 'scope'
    KEY = 'key'
    TYPE = 'type'
    
class JOB(schema_obj):
    resource_name = 'job'
    
    ID = 'id'

    USERNAME = 'username'
    URI = 'uri'
    METHOD = 'method'
    ENCODING = 'encoding'
    CONTENT_TYPE = 'content_type'
    HTTP_ACCEPT = 'http_accept'
    PARAMS = 'params'
    COMMENT = 'comment'
    CONTEXT_DATA = 'context_data'
    
    PROCESS_ID = 'process_id'
    PROCESS_ENV = 'process_env'
    PROCESS_MESSAGES = 'process_messages'
    STATE = 'state'
    DATE_TIME_REQUESTED = 'date_time_requested'
    DATE_TIME_SUBMITTED = 'date_time_submitted'
    DATE_TIME_PROCESSING = 'date_time_processing'
    DATE_TIME_COMPLETED = 'date_time_completed'
    
    RESPONSE_CONTENT = 'response_content'
    RESPONSE_STATUS_CODE = 'response_status_code'
    
    JOB_PROCESSING_FLAG = 'job_id'
    
class APILOG(schema_obj):
    resource_name = 'apilog'
    
    USERNAME = 'username'
    USER_ID = 'user_id'
    URI = 'uri'
    KEY = 'key'
    REF_RESOURCE_NAME = 'ref_resource_name'
    API_ACTION = 'api_action'
    DATE_TIME = 'date_time'
    DIFF_KEYS = 'diff_keys'
    DIFFS = 'diffs'
    JSON_FIELD = 'json_field'
    COMMENT = 'comment'
    CHILD_LOGS = 'child_logs'
    ID = 'id'
    PARENT_LOG_URI = 'parent_log_uri'
    PARENT_LOG_ID = 'parent_log_id'
    LOG_URI = 'log_uri'
    IS_PREVIEW = 'is_preview'

    
class VOCAB(schema_obj):
    ''' Define selected vocabulary constants used by the API.'''

    class filter_type(schema_obj):
        ''' Define API list filters'''
        
        # EXACT and EQUAL are synonyms
        EXACT = 'exact'
        EQUAL = 'eq'
        
        IS_NULL = 'is_null'
        IS_BLANK = 'is_blank'
        
        CONTAINS = 'contains'
        # Case insentive contains
        ICONTAINS = 'icontains'
        
        # Filter for items in the given list
        IN = 'in'
        
        # Numerical filters
        LESS_THAN = 'lt'
        LESS_THAN_EQUAL = 'lte'
        GREATER_THAN = 'gt'
        GREATER_THAN_EQUAL = 'gte'
        NOT_EQUAL = 'ne'
        
        # RANGE uses two part [begin,end] values
        RANGE = 'range'
        # ABOUT compares the filter value to the target, rounded to the same precision
        ABOUT = 'about'
        
        # INVERTED may be prepended to any field name to invert the type
        INVERTED = '-'
        
        # Separator used to split filter strings apart.
        LOOKUP_SEP = '__'
        
    class resource(schema_obj):
        class visibility(schema_obj):
            LIST = 'l'
            DETAIL = 'd'
            EDIT = 'e'
            CREATE = 'c'
            UPDATE = 'u'
    class field(schema_obj):
        class data_type(schema_obj):
            STRING = 'string'
            BOOLEAN = 'boolean'
            DATE = 'date'
            DATETIME = 'datetime'
            LIST = 'list'
            FLOAT = 'float'
            DECIMAL = 'decimal'
            INTEGER = 'integer'
    
            NUMERIC_TYPES = [FLOAT, DECIMAL, INTEGER]
        
        class display_type(schema_obj):
            LINK = 'link'
            IMAGE = 'image'
            SIUNIT = 'siunit'
            COMMENT_ARRAY = 'comment_array'
            FULL_STRING = 'full_string'
        class edit_type(schema_obj):
            SELECT = 'select'
            MULTISELECT = 'multiselect'
            MULTISELECT2 = 'multiselect2'
            MULTISELECT3 = 'multiselect3'
            TEXTAREA = 'textarea'
            TYPEAHEAD = 'typeahead'
            CUSTOM = 'custom'
        class visibility(schema_obj):
            LIST = 'l'
            DETAIL = 'd'
            EDIT = 'e'
            SUMMARY = 'summary'
            BILLING = 'billing'
            PROTOCOL = 'protocol'
            API = 'api'
            NONE = 'none'
            
            hidden_fields = [API, NONE]
            
        class editability(schema_obj):
            CREATE = 'c'
            UPDATE = 'u'
            LIST_UPDATE = 'l'
            
    class permission(schema_obj):
        class type(schema_obj):
            READ = 'read'
            WRITE = 'write'
            
    class job(schema_obj):
        class state(schema_obj):
            PENDING = 'pending'
            SUBMITTED = 'submitted'
            PROCESSING = 'processing'
            COMPLETED = 'completed'
            FAILED = 'failed'
    
    class apilog(schema_obj):
        class api_action(schema_obj):
            POST = 'POST'
            PUT = 'PUT'
            # NOTE: "CREATE" - to distinguish PATCH/modify, PATCH/create
            CREATE = 'CREATE' 
            PATCH = 'PATCH'
            DELETE = 'DELETE'
            

def parse_display_options(field):
    display_options = field.get('display_options',None)
    if display_options:
        try:
            display_options = json.loads(display_options.replace("'",'"'))
            field['display_options'] = display_options
            
        except Exception,e:
            logger.exception('schema error for field: %r, display_options: %r',
                key, display_options)
            display_options = None
            del field['display_options']
    return display_options

def parse_schema(schema):
    '''
    This is a tool for parsing the display_options fields on the client side.
    '''

    fields = schema[RESOURCE.FIELDS]
    
    # Parse the display options
    for key,field in fields.items():
        parse_display_options(field)
    logger.info('retrieved schema with fields: %r', fields.keys())
    return schema

def get_title(key, schema):
    if key not in schema[RESOURCE.FIELDS]:
        logger.warn('no field: %s in schema: %r', key, schema)
        return key
    return str(schema[RESOURCE.FIELDS][key]['title'])
    

def get_href(value, key, schema, record):
    ''' Generate the link for the given field, it possible, or return the value'''
    if value is None:
        return None
    if isinstance(value, (list,tuple)):
        logger.warn('list values are not supported (yet): %s: %r', key, value)
        return value
    val = value
    if key not in schema[RESOURCE.FIELDS]:
        return value
    field = schema[RESOURCE.FIELDS][key]
    if 'display_options' in field:
        display_options = field['display_options']
        if 'hrefTemplate' in display_options:
            hrefTemplate = display_options['hrefTemplate']
            try:
                href = hrefTemplate.format(**record)
                href = '/'.join([
                    settings.APP_PUBLIC_DATA.site_url, href])
                val = '<a href="%s">%s</a>' % ( 
                    href,value)
            except:
                logger.exception('hrefTemplate error: %r (%s), %r, %r', 
                    value, field['key'], record, display_options)
    return val

def get_vocab_title(key, value, schema):
    ''' get the vocab "title" for the value if "vocabulary" is available
    for the field in the schema'''
    if value is None:
        return None
    if isinstance(value, (list,tuple)):
        logger.warn('list values are not supported (yet): %s: %r', key, value)
        return value
    title = value
    if key in schema[RESOURCE.FIELDS]:
        field = schema[RESOURCE.FIELDS][key]
        if 'vocabulary' in field:
            vocab = field['vocabulary']
            if value in vocab:
                title = vocab[value]['title']
            elif str(value) in vocab:
                title = vocab[str(value)]['title']
            else:
                logger.warn('vocab not found: %s: %r, in %r',
                    key, value, vocab)
        else:
            logger.debug('get_vocab_title: no "vocabulary" for %s '
                'in the schema field: %s',key, field)
    else:
        logger.warn('get_vocab_title: no %r field in the schema: %r',
            key, schema['key'])
    return title

def replace_vocabularies(record, schema):
    ''' replace vocabularies for the given record '''
    new_record = dict(record)
    for k,v in new_record.items():
        new_record[k] = get_vocab_title(k, v, schema)
    return new_record

def replace_html_values(record, schema):
    ''' replace vocabularies, and create hrefs for the given record '''
    html_dict = dict(record)
    for k,v in html_dict.items():
        
        html_dict[k] = get_href(
            get_vocab_title(k,v,schema),
            k, schema, record)
    return html_dict