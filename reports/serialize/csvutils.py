from __future__ import unicode_literals

import csv
import logging

from reports.serialize import to_simple
from django.utils.encoding import smart_text, force_text


logger = logging.getLogger(__name__)

LIST_DELIMITER_CSV = ';'

# Use the csv package example for reading data:
# - Note: data are already presented as unicode at this point, so the
# unicodecsv package (expects bytes) is not appropriate.
def unicode_csv_reader(unicode_csv_data, dialect=csv.excel, **kwargs):
    # csv.py doesn't do Unicode; encode temporarily as UTF-8:
    csv_reader = csv.reader(utf_8_encoder(unicode_csv_data),
                            dialect=dialect, **kwargs)
    for row in csv_reader:
        # decode UTF-8 back to Unicode, cell by cell:
        yield [unicode(cell, 'utf-8') for cell in row]

def utf_8_encoder(unicode_csv_data):
    for line in unicode_csv_data:
        yield line.encode('utf-8')

def from_csv(csvfile, list_delimiter=LIST_DELIMITER_CSV, list_keys=None):
    '''
    Returns an in memory matrix (array of arrays) for the input file
    
    @param list_keys overrides nested list eval for column keys; no brackets '[]' are 
        needed to denote these columns as list columns - however, to comply with 
        the csv standard, they still have to be quoted (if list_delimiter=csv_delimiter)
    NOTES: 
    - nested lists are denoted by brackets, i.e. '[]',
    - to escape use '\[...' (i.e. when embedding a regex expression)
    TODO: version 2 - read from a stream
    '''
    # reader = csv.reader(csvfile, encoding='utf-8')
    reader = unicode_csv_reader(csvfile)
    return from_csv_iterate(reader, list_delimiter=list_delimiter, list_keys=list_keys)
    
def csv_generator(iterable, list_delimiter=LIST_DELIMITER_CSV, list_keys=None):
    '''
    Perform LIMS specific nested list serialization:
    - otherwise yields the iterable items in a generator.
    - FIXME: the name "csv_generator" is misleading; no csv serialization 
    operations are being performed here. propose: "nested_list_generator"
    '''
    list_keys = list_keys or []
    list_keys = list(list_keys)
    i = 0 
    for row in iterable:
        if i == 0:
            keys = [x for x in row]
        else:
            item = dict(zip(keys,row))
            for key in item.keys():
                val = item[key]
                if val and len(val)> 1:
                    if val[0] == '\\' and val[1] == '[':
                        # this could denote an escaped bracket, i.e. for a regex
                        item[key] = val[1:]
                    elif key in list_keys or val[0]=='[':
                        # due to the simplicity of the serializer, above, any 
                        # quoted string is a nested list
                        list_keys.append(key)
                        item[key] = [
                            x.strip() 
                            for x in val.strip('"[]').split(list_delimiter)]
                        logger.debug('%r, val: %r converted to %r', key, val, item[key])
            yield item
        i += 1
    logger.debug('read in data, count: %d', i )   
    
def from_csv_iterate(iterable, list_delimiter=LIST_DELIMITER_CSV, list_keys=None):
    '''
    Returns an in memory array of dicts for the iterable, representing a 
    csv-like input matrix.
    - the first row is interpreted as the dict keys, unless a list_keys param is 
    specified 
    
    FIXME: the name "from_csv_iterate" is misleading, as no csv-specific
    deserialization operations are going on here. propose: "nested_list_reader"
    '''
    list_keys = list_keys or []
    data_result = []
    i = 0
    keys = []
    list_keys = list(list_keys) 
    logger.debug('list_keys: %r', list_keys)
    for row in iterable:
        if i == 0:
            keys = [x for x in row]
        else:
            item = dict(zip(keys,row))
            for key in item.keys():
                val = item[key]
                if val and len(val)> 1:
                    if val[0] == '\\' and val[1] == '[':
                        # this could denote an escaped bracket, i.e. for a regex
                        item[key] = val[1:]
                    elif key in list_keys or val[0]=='[':
                        # due to the simplicity of the serializer, above, any 
                        # quoted string is a nested list
                        list_keys.append(key)
                        item[key] = []
                        for x in val.strip('"[]').split(list_delimiter):
                            x = x.strip()
                            if x:
                                item[key].append(x)
            data_result.append(item)
        i += 1
    logger.debug('read in data, count: ' + str(len(data_result)) )   
    return data_result

def dict_to_rows(_dict):
    ''' Utility that converts a dict into a table for writing to a spreadsheet
    '''
    
    logger.debug('_dict: %r', _dict)
    values = []
    if isinstance(_dict, dict):
        for key,val in _dict.items():
            for row in dict_to_rows(val):
                if not row:
                    values.append([key,None])
                else:
                    keyrow = [key]
                    if isinstance(row, basestring):
                        keyrow.append(row)
                    else:
                        keyrow.extend(row)
                    values.append(keyrow)
    else:
        values = (convert_list_vals(_dict),)
    return values

def convert_list_vals(val, delimiter=LIST_DELIMITER_CSV, list_brackets='[]'):
    delimiter = delimiter + ' '
    if isinstance(val, (list,tuple)):
        if list_brackets:
            return ( list_brackets[0] 
                + delimiter.join([smart_text(to_simple(x)) for x in val]) 
                + list_brackets[1] )
        else: 
            return delimiter.join([smart_text(to_simple(x)) for x in val]) 
    elif val != None:
        if type(val) == bool:
            if val:
                return 'TRUE'
            else:
                return 'FALSE'
        else:
            return force_text(to_simple(val))
    else:
        return None

def csv_convert_list_vals(val, delimiter=LIST_DELIMITER_CSV, list_brackets='[]'):
    '''
    Convert values for writing (to csv)
    NOTE: python 2 csv.writer does not support Unicode: all values must be 
    converted to 8-bit bytestrings to write.
    '''
    delimiter = delimiter + ' '
    if isinstance(val, (list,tuple)):
        if list_brackets:
            return ( list_brackets[0] 
                + delimiter.join([smart_text(to_simple(x)).encode('utf-8') for x in val]) 
                + list_brackets[1] )
        else: 
            return delimiter.join([smart_text(to_simple(x)).encode('utf-8') for x in val]) 
    elif val != None:
        if type(val) == bool:
            if val:
                return 'TRUE'
            else:
                return 'FALSE'
        else:
            x = force_text(to_simple(val))
            # MUST encode for csv (understands bytes only)
            y = x.encode('utf-8')
            return y
    else:
        return None
