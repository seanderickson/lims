from __future__ import unicode_literals

LIST_DELIMITER_CSV = ';'
LIST_DELIMITER_XLS = ';'
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

