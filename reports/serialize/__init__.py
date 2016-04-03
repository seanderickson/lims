from __future__ import unicode_literals

import dateutil
import six
import logging
from reports import ValidationError
from reports.serialize import csvutils

logger = logging.getLogger(__name__)

JSON_MIMETYPE = 'application/json'
CSV_MIMETYPE = 'text/csv'
XLS_MIMETYPE = 'application/xls'
XLSX_MIMETYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
SDF_MIMETYPE = 'chemical/x-mdl-sdfile'

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
        values = (csvutils.csv_convert(_dict),)
    return values


def parse_val(value, key, data_type):
    """
    All values are read as strings from the input files, so this function 
    converts them as directed.
    TODO: validation
    """
    try:
        if ( value is None 
            or value == '' 
            or value == 'None' 
            or value == u'None' 
            or value == u'n/a'):
            if data_type == 'string': 
                return ''
            else:  
                return None
        if data_type == 'string':
            return value
        elif data_type == 'integer':
            # todo: this is a kludge, create an integer from values like "5.0"
            return int(float(value))
        elif data_type == 'date':
            return dateutil.parser.parse(value)
        elif data_type == 'datetime':
            return dateutil.parser.parse(value)
        elif data_type == 'boolean':
            if value is True or value is False:
                 return value
            value = str(value)
            if(value.lower() == 'true' 
                or value.lower() == 't' or value == '1'): return True
            return False
        elif data_type in ['float','decimal']:
            return float(value)
        elif data_type == 'list':
            if isinstance(value, six.string_types):
                if value.strip():
                    return (value,) # convert string to list
                else:
                    return None
            return value # otherwise, better be a list
        else:
            raise Exception('unknown data type: %s: "%s"' % (key,data_type))
    except Exception, e:
        logger.exception('value not parsed %r:%r',key, value)
        raise ValidationError(key=key,msg='parse error: %r' % str(e))

def parse_json_field(val, key, json_field_type):
    'Parse a field nested in the json_obj'
     
    # FIXME: now that tastypie is removed, 
    # json_field_type should equal data_type
    if json_field_type == 'fields.CharField':
        return parse_val(val, key, 'string')
    elif json_field_type == 'fields.ListField':
        return parse_val(val, key, 'list')
    elif json_field_type == 'CsvBooleanField':                    
        return parse_val(val, key, 'boolean')
    elif json_field_type == 'fields.BooleanField':
        return parse_val(val, key, 'boolean')
    elif json_field_type == 'fields.IntegerField':
        return parse_val(val, key, 'integer')
    elif json_field_type == 'fields.DateField':
        raise NotImplementedError(
            'json DateField field is not implemented')
    elif json_field_type == 'fields.DecimalField':
        raise NotImplementedError(
            'json decimal field is not implemented')
    elif json_field_type == 'fields.FloatField':
        raise NotImplementedError(
            'json float field is not implemented')
    else:
        raise NotImplementedError(
            'unknown json_field_type: %s' % json_field_type)

