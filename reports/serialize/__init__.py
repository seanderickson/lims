from __future__ import unicode_literals

import datetime
import io
import logging

from PIL import Image
import dateutil
from django.conf import settings
from django.core.serializers.json import DjangoJSONEncoder
from django.core.urlresolvers import resolve
from django.utils.encoding import force_text
import pytz
import six

from reports import ValidationError
from decimal import Decimal


logger = logging.getLogger(__name__)

JSON_MIMETYPE = 'application/json'
CSV_MIMETYPE = 'text/csv'
XLS_MIMETYPE = 'application/xls'
XLSX_MIMETYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
SDF_MIMETYPE = 'chemical/x-mdl-sdfile'

def make_local_time(datetime_obj):
    '''
    Motivation: convert UST datetimes to the local timezone
    '''
    if settings.USE_TZ and settings.TIME_ZONE:
        if isinstance(datetime_obj, datetime.datetime):
            return datetime_obj.astimezone(pytz.timezone(settings.TIME_ZONE))
    return datetime_obj


class LimsJSONEncoder(DjangoJSONEncoder):

    def default(self, o):
        if isinstance(o, datetime.datetime):
            o = make_local_time(o)
        return DjangoJSONEncoder.default(self, o)
    
json_encoder = LimsJSONEncoder()

def to_simple(data):
    """
    This brings complex Python data structures down to native types of the
    serialization format(s).
    """
    if isinstance(data, (list, tuple)):
        return [to_simple(item) for item in data]
    if isinstance(data, dict):
        return dict((key, to_simple(val)) for (key, val) in data.items())
    elif isinstance(data, datetime.datetime):
        return json_encoder.default(data)
    elif isinstance(data, datetime.date):
        return json_encoder.default(data)
    elif isinstance(data, datetime.time):
        return json_encoder.default(data)
    elif isinstance(data, bool):
        return data
    elif isinstance(data, (six.integer_types, float)):
        return data
    elif data is None:
        return None
    else:
        return force_text(data)

def parse_val(value, key, data_type, options=None):
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
            or value == 'null'
            or value == u'n/a'):
            if data_type == 'string': 
                return ''
            elif data_type == 'list':
                return []
            else:  
                return None
        if data_type == 'string':
            return value
        elif data_type == 'integer':
            # todo: this is a kludge, create an integer from values like "5.0"
            return int(float(value))
        elif data_type == 'date':
            return dateutil.parser.parse(value).date()
        elif data_type == 'datetime':
            return dateutil.parser.parse(value)
        elif data_type == 'boolean':
            if value is True or value is False:
                 return value
            value = str(value)
            if(value.lower() == 'true' 
                or value.lower() == 't' or value == '1'): return True
            return False
        elif data_type == 'float':
            return float(value)
        elif data_type == 'decimal':
            if isinstance(value, float):
                logger.warn('converting float: %r to decimal: %r',
                    value, Decimal(str(value)))
                value = str(value)
            return Decimal(value)
        elif data_type == 'list':
            if isinstance(value, six.string_types):
                if value.strip():
                    return (value,) # convert string to list
                else:
                    return []
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

def resolve_image(request, uri):
    logger.debug('find image at %r', uri)
    view, args, kwargs = resolve(uri)
    kwargs['request'] = request
    response = view(*args, **kwargs)
    return Image.open(io.BytesIO(response.content))
