from __future__ import unicode_literals

import django.core.exceptions
import django.core.serializers.json
import django.utils.timezone

from collections import OrderedDict
import logging

# Note the csv package does not allow multibyte delimiters 
CSV_DELIMITER = b','  
LIST_DELIMITER_SQL_ARRAY = '~^'
LIST_DELIMITER_URL_PARAM = ','
LIST_DELIMITER_SUB_ARRAY = '$'
MAX_ROWS_PER_XLS_FILE = 100000
MAX_IMAGE_ROWS_PER_XLS_FILE = 1000

HTTP_PARAM_USE_VOCAB = 'use_vocabularies'
HTTP_PARAM_USE_TITLES = 'use_titles'
HTTP_PARAM_RAW_LISTS = 'raw_lists'
HTTP_PARAM_DATA_INTERCHANGE = 'data_interchange'

API_RESULT_ERROR = 'errors'

# Header custom comment field
HEADER_APILOG_COMMENT = 'HTTP_X_APILOG_COMMENT'

LIST_BRACKETS = '[]' # default char to surround nested list in xls, csv

logger = logging.getLogger(__name__)

class ValidationError(Exception):
    def __init__(self,errors=None, key=None, msg=None):
        
        assert errors is not None or (key and msg),( 
            'ValidationError initialization requires "errors" parameter')
         
        self.errors = errors or {}
         
        if key:
            if not isinstance(msg, (list,tuple)):
                msg = [msg]
            self.errors[key] = msg
     
    def __repr__(self, *args, **kwargs):
        return 'validation errors: %r' % self.errors

class InformationError(ValidationError):
    pass

class ParseError(ValidationError):
    pass

def _now():
    d = django.utils.timezone.now()
    if d.tzinfo:
        d = django.utils.timezone.localtime(d)
    logger.debug('timezone: %r, %r', d.tzinfo, d )
    return d

def strftime_log_manual(d):
    ''' format a datetime d as an ApiLog time:
        - the same as isoformat with millisecond precision
        - logs support millisecond precision due to postgresql limitations
    '''
    date_time_part = d.strftime('%Y-%m-%dT%H:%M:%S')
    milliseconds = d.microsecond
    if milliseconds > 0:
        milliseconds = milliseconds/1000
    utc_offset = d.strftime('%z')
    utc_offset = '%s:%s' % (utc_offset[:-2],utc_offset[-2:])
    return '%s.%d%s' % (date_time_part, milliseconds, utc_offset)

djangoJsonEncoder = django.core.serializers.json.DjangoJSONEncoder()

def strftime_log(d):
    ''' format a datetime d as an ApiLog time:
        - the same as isoformat with millisecond precision
        - logs support millisecond precision due to postgresql limitations
        - use the DjangoJSONEncoder, to give the same result as the serializers
    '''
    return djangoJsonEncoder.default(d)
    
    

