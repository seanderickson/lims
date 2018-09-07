# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from collections import defaultdict
import logging

import django.core.exceptions
import django.core.serializers.json
import django.utils.timezone

import reports.schema as SCHEMA


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

API_RESULT_ERROR = SCHEMA.ERROR.resource_name

# NOTE: the simulated django client requests expect the HTTP Header "Accept" to 
# be stored in the variable "HTTP_ACCEPT"
DJANGO_ACCEPT_PARAM = 'HTTP_ACCEPT'

HTTP_PARAM_AUTH = 'HTTP_AUTHORIZATION'
HTTP_PARAM_CONTENT_TYPE = 'CONTENT_TYPE'

# Header custom comment field
HEADER_APILOG_COMMENT = 'HTTP_X_APILOG_COMMENT'
# Header for custom comment to be used from clients 
# (Django translates custom headers, adding "HTTP_" and converting to underscores)
HEADER_APILOG_COMMENT_CLIENT = 'X-APILOG-COMMENT'

LIST_BRACKETS = '[]' # default char to surround nested list in xls, csv

logger = logging.getLogger(__name__)

class ValidationError(Exception):
    def __init__(self,errors=None, key=None, msg=None, input_id=None):
        
        assert errors is not None or (key and msg),( 
            'ValidationError initialization requires "errors" parameter')
         
        if errors is not None:
            if not isinstance(errors, dict):
                self.errors = { 'errors': errors }
            else:
                self.errors = errors
        else:
            self.errors = {}
        if key:
            if not isinstance(msg, (list,tuple,dict)):
                msg = [msg]
            self.errors[key] = msg
        if input_id:
            self.errors['input_id'] = input_id
            
    def __repr__(self, *args, **kwargs):
        return 'validation errors: %r' % self.errors
    def __str__(self, *args, **kwargs):
        return self.__repr__(*args, **kwargs)

class CumulativeError(ValidationError):
    
    def __init__(self):
        
        ValidationError.__init__(self, errors=defaultdict(dict))
        
    def add_error(self, key, new_errors, line=None):
        if line:
            new_errors[SCHEMA.ERROR.LINE] = line
        self.errors[key].update(new_errors)

    def _update_from(self, new_errors):
        # update a two-level dict
        if self.errors:
            for key, error_dict in new_errors.items():
                self.errors[key].update(error_dict)
        else:
            self.errors.update(new_errors)
    
    def update_from(self, cumulative_error):
        self._update_from(cumulative_error.errors)
        

class InformationError(ValidationError):
    pass

class ParseError(ValidationError):
    pass

class BadRequestError(ValidationError):
    pass

class MissingParam(ValidationError):
    
    def __init__(self, param_name):
        ValidationError.__init__(key=param_name, msg='required')
        
class ApiNotImplemented(ValidationError):

    def __init__(self, resource_name, method_name, errors=None):
        
        errors = errors or {}
        errors.setdefault('message','API is not implemented')
        errors['resource_name'] = resource_name
        errors['method_name'] = method_name
        ValidationError.__init__(self, errors=errors)

class BackgroundJobImmediateResponse(Exception):
    
    def __init__(self, httpresponse):
        self.httpresponse = httpresponse

class LoginFailedException(django.core.exceptions.ValidationError):
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
    
    

