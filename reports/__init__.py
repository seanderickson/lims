from __future__ import unicode_literals

import django.core.exceptions
from collections import OrderedDict
import logging

# Note the csv package does not allow multibyte delimiters 
CSV_DELIMITER = b','  
LIST_DELIMITER_SQL_ARRAY = ';'
LIST_DELIMITER_URL_PARAM = ','
MAX_ROWS_PER_XLS_FILE = 100000
MAX_IMAGE_ROWS_PER_XLS_FILE = 2000

HTTP_PARAM_USE_VOCAB = 'use_vocabularies'
HTTP_PARAM_USE_TITLES = 'use_titles'
HTTP_PARAM_RAW_LISTS = 'raw_lists'

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
            self.errors[key] = [msg]
     
    def __repr__(self, *args, **kwargs):
        return 'validation errors: %r' % self.errors
    
class ParseError(ValidationError):
    pass

