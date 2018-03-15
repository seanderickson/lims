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

API_PARAM_SEARCH = 'search'
# URI_PARAM_RAW_SEARCH = 'raw_search_data'
# Complex search - a search data structure sent as a POST header 
API_PARAM_COMPLEX_SEARCH_ID = 'search_id'
URI_PATH_COMPLEX_SEARCH = 'csearch'

# Date format for API - time zone is not used for dates
DATE_FORMAT = "%Y-%m-%d"

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

# Define API resources
class RESOURCE(schema_obj):
    ''' Constants for the Resource resource: constant field names'''
    resource_name = 'resource'

    FIELDS = 'fields'
    
class VOCAB(schema_obj):
    
    class resource(schema_obj):
        class visibility(schema_obj):
            LIST = 'l'
            DETAIL = 'd'
            EDIT = 'e'
            CREATE = 'c'
            UPDATE = 'u'
    

def parse_schema(schema):

    fields = schema[RESOURCE.FIELDS]
    
    # Parse the display options
    for key,field in fields.items():
        display_options = field.get('display_options',None)
        if display_options:
            try:
                display_options = json.loads(display_options.replace("'",'"'))
                field['display_options'] = display_options
            except Exception,e:
                logger.exception('schema error for field: %r, display_options: %r',
                    key, display_options)
                del field['display_options']
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