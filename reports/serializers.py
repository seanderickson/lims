import csv
import copy
from collections import OrderedDict
import cStringIO
import StringIO
import json
import logging

from django.utils.encoding import smart_str
from django.core.serializers.json import DjangoJSONEncoder
from django.db.backends.util import CursorDebugWrapper
from psycopg2.psycopg1 import cursor
from tastypie import fields
from tastypie.bundle import Bundle
from tastypie.serializers import Serializer

import reports.utils.sdf2py as s2p
import reports.utils.serialize

logger = logging.getLogger(__name__)

class CsvBooleanField(fields.ApiField):
    """
    Because csv is not json, have to do some fudging with booleans;
    basically, for strings, case insensitive "true" is True, 
    all other values are False.
    Non-strings are interpreted as usual, using bool(val).
    """
    dehydrated_type = 'boolean'
    help_text = 'Boolean data. Ex: True'

    def convert(self, value):
        if value is None:
            return None
        if isinstance(value, basestring):
            if value.lower() == 'true':
                return True
            return False
        else:
            return bool(value)


class BackboneSerializer(Serializer):
    
    def from_json(self, content):
        """
        Given some JSON data, returns a Python dictionary of the decoded data.
        Override to quote attributes - the client is failing to.
        """
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(("loading content:", content)))
        content = content.replace(r'(\w+):', r'"\1" :')
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(("loading content:", content)))
        return json.loads(content)

# NOTE: removed this class as this is the stock behavior in newer tastypie
# class TimeZoneAwareDateSerializer(Serializer):
#     """
#     Tastypie converts all datetimes to timezone-naive (UTC); this serializer 
#     uses the full ISO 8601 format if the value isn't already naive.
#     Note: even when the date is stripped of timezone information (tz naive), 
#     when it is formatted for ISO 8601, it becomes a UTC time: as in:
#     2010-12-16T00:00:00
#     """
#     def format_datetime(self, data): 
#         if is_naive(data) or self.datetime_formatting == 'rfc-2822':
#             return super(TimeZoneAwareDateSerializer, self).format_datetime(data)
#  
#         return data.isoformat()

class PrettyJSONSerializer(Serializer):
    json_indent = 2

    def to_json(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        return json.dumps(data, cls=DjangoJSONEncoder,
                sort_keys=True, ensure_ascii=False, indent=self.json_indent)

class SDFSerializer(Serializer):
    
    def __init__(self, content_types=None, formats=None, **kwargs):

        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['sdf'] = 'chemical/x-mdl-sdfile'
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            _formats.append('sdf')
            formats = _formats
            
        super(SDFSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);
        
    def to_sdf(self, data):
        logger.warn(str(('to_sdf', data)))
        output = cStringIO.StringIO()
        MOLDATAKEY = s2p.MOLDATAKEY
        for d in data:
            output.write('$$$$\n')
            
            if MOLDATAKEY in d:
                output.write(str(d[MOLDATAKEY]))
                output.write('\n') 
                del d[MOLDATAKEY]
            logger.info(str(('d',d)))
            for k,v in d.items():
                output.write('> <%s>\n' % k)
                output.write(str(v))
                output.write('\n\n')
        return output.getvalue()
    
    def from_sdf(self, content, root=None):
        return { 'objects': s2p.parse_sdf(content) }

# class MultiPartDeserializser():
#     
#     def deserialize(self, content, format='application/json'):
#     """
#     Given some data and a format, calls the correct method to deserialize
#     the data and returns the result.
#     """
#     desired_format = None
# 
#     format = format.split(';')[0]
# 
#     for short_format, long_format in self.content_types.items():
#         if format == long_format:
#             if hasattr(self, "from_%s" % short_format):
#                 desired_format = short_format
#                 break
# 
#     if desired_format is None:
#         raise UnsupportedFormat("The format indicated '%s' had no available 
#             deserialization method. Please check your ``formats`` and
#             ``content_types`` on your Serializer." % format)
# 
#     if isinstance(content, six.binary_type):
#         content = force_text(content)
# 
#     deserialized = getattr(self, "from_%s" % desired_format)(content)
#     return deserialized



class CSVSerializer(Serializer):
    
    def __init__(self, content_types=None, formats=None, **kwargs):
        
        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['csv'] = 'text/csv'
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            _formats.append('csv')
            formats = _formats
            
        super(CSVSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);
        
        
    def to_csv(self, data, root='objects', options=None):
        '''
        @param root where the return object iterable is nested in the data
            object (None if no nesting, and data is iterable).
        
        '''
        
        options = options or {}
        data = self.to_simple(data, options)
        raw_data = StringIO.StringIO()
        # default: delimiter = ',' quotechar='"'
        writer = csv.writer(raw_data) 

        if 'error' in data:
            writer.writerow(['error'])
            writer.writerow([data['error']])
            logger.warn(str(('error', data)))
            return raw_data.getvalue()
            
        # TODO: stream this, don't do the whole dict at once 
        if 'objects' in data:
            data = data['objects']
        if len(data) == 0:
            return raw_data

        if isinstance(data, dict):
            # usually, this happens when the data is actually an error message;
            # but also, it could be just one item being returned
            keys = data.keys()
            writer.writerow([smart_str(key) for key in keys])
            writer.writerow(self.get_list(data))
        else:    
            # default 
            i = 0
            keys = None
            for item in data:
                if i == 0:
                    keys = item.keys()
                    writer.writerow([smart_str(key) for key in keys])
                i += 1
                writer.writerow(self.get_list(item))

        return raw_data.getvalue()
    
    def get_list(self,item):
        '''
        Convert a csv row into a list of values
        '''
        _list = []
        for key in item:
            logger.debug(str(('item', item)))
            if item[key] and isinstance(item[key], (list, tuple)):
                _list.append(
                    '[' + ','.join([smart_str(x) for x in item[key]]) + ']' )
            elif item[key] != None:
                val = item[key]
                if type(val) == bool:
                    if val:
                        _list.append('TRUE')
                    else:
                        _list.append('FALSE')
                else:
                    _list.append(smart_str(item[key]))
            else:
                _list.append(None)
        return _list
    
    def from_csv(self, content, root='objects'):
        '''
        @param root - property to nest the return object iterable in for the 
            response (None if no nesting, and return object will be an iterable)
        TODO: version 2 - read from a stream
        '''
        objects = reports.utils.serialize.from_csv(StringIO.StringIO(content))
        if root:
            return { root: objects }
        else:
            return objects
                
    
#     def from_csv(self, content, root='objects'):
#         '''
#         @param root - property to nest the return object iterable in for the 
#             response (None if no nesting, and return object will be an iterable)
#         
#         '''
#         raw_data = StringIO.StringIO(content)
#         
#         # TODO: also, stream this
#         # default: delimiter = ',' quotechar='"'
#         logger.debug('reading...')
#         reader = csv.reader(raw_data)
#         return self.from_csv_iterate(reader, root)
#         
#     def from_csv_iterate(self, iterable, root='objects'):
#         data_result = []
#         i = 0
#         keys = []
#         list_keys = [] # cache
#         for row in iterable:
#             if i == 0:
#                 keys = row
#             else:
#                 item = dict(zip(keys,row))
#                 for key in item.keys():
#                     val = item[key]
#                     if ( val and len(val)> 1 and 
#                             (key in list_keys or val[0]=='[') ):
#                         # due to the simplicity of the serializer, above, any 
#                         # quoted string is a nested list
#                         list_keys.append(key)
#                         item[key] = val.strip('"[]').split(',')
#                 data_result.append(item)
#             i += 1
#         logger.debug('read in data, count: ' + str(len(data_result)) )   
# 
#         if root:
#             return { root: data_result }
#         else:
#             return data_result



class CursorSerializer(Serializer):
    """
    A simple serializer that takes a cursor, queries it for its columns, and
    outputs this as either CSV or JSON.
    (The CSV output is used with SAF)
    """
    
    formats = ['json', 'jsonp', 'xml', 'yaml', 'csv']
    content_types = {
        'json': 'application/json',
        'jsonp': 'text/javascript',
        'xml': 'application/xml',
        'yaml': 'text/yaml',
        'csv': 'text/csv',
    }
  
    def to_csv(self, bundle_or_obj, options=None):
        '''
        NOTE: ignores all content except for "objects"
        '''
        
        raw_data = StringIO.StringIO()

        obj = bundle_or_obj            
        if isinstance(bundle_or_obj, Bundle):
            obj = bundle_or_obj.obj
            
        # this is an unexpected way to get this error, look into tastypie
        if(isinstance(obj,dict) and 'error_message' in obj):
            logger.warn(str(('report error', obj)))
            raw_data.writelines(('error_message\n',obj['error_message'],'\n'))
            return raw_data.getvalue() 
        elif isinstance(obj,dict) :
            
            writer = csv.writer(raw_data)
            wrote_to_csv = False
            for key, value in obj.items():
                if isinstance(value, cursor):
                    self._cursor_to_csv(value, writer)
                    wrote_to_csv = True
            
            for key, value in obj.items():
                if not isinstance(value, cursor):
                    if not wrote_to_csv:
                        writer.writewrow([key,json.dumps(
                            value,skipkeys=False,check_circular=True, 
                            allow_nan=True, default=lambda x: str(x))] )
                    else:
                        logger.warn(
                            'non-cursor data will not be written to csv: "' + 
                            key +'": ' + json.dumps(
                                value,skipkeys=False,check_circular=True, 
                                allow_nan=True, default=lambda x: str(x)))
            
        return raw_data.getvalue()
    
    def _cursor_to_csv(self, cursor, writer):

        i=0
        cols = [col[0] for col in cursor.description]
        # TODO: grab the column names here from the meta information store?
        writer.writerow(cols)

        for row in cursor.fetchall():
            writer.writerow(
                [smart_str(val, 'utf-8', errors='ignore') for val in row])
            i += 1
        logger.info('_cursor_to_csv done, wrote: %d' % i)

    def to_json(self,bundle_or_obj, options=None):
        
        logger.info(str(('typeof the object sent to_json',type(bundle_or_obj))))
        raw_data = StringIO.StringIO()
         
        if isinstance(bundle_or_obj, Bundle):
            obj = bundle_or_obj.obj
        else:
            obj = bundle_or_obj            
            
        if isinstance(obj,dict) and 'error_message' in obj :
            logger.warn(str(('report error', obj)))
            raw_data.writelines(('error_message\n',obj['error_message'],'\n'))
            return raw_data.getvalue() 
        elif isinstance(obj,dict) :
            raw_data.write('{');
            count = 0
            for key, value in obj.items():
                if count > 0: raw_data.write(', ')
                if isinstance(value, (cursor, CursorDebugWrapper) ):
                    raw_data.write('"' + key + '": [') # key should be 'objects'
                    self._cursor_to_json(value, raw_data)
                    raw_data.write(']')
                else:
                    raw_data.write('"' + key + '": ')
                    raw_data.write(json.dumps(
                        value,skipkeys=False,check_circular=True, 
                        allow_nan=True, default=lambda x: str(x)))
                count += 1
                
            raw_data.write('}')
                    
        return raw_data.getvalue() 
        
    def _cursor_to_json(self, _cursor, raw_data):
        if not isinstance(_cursor, (cursor, CursorDebugWrapper) ):
            raise Exception(
                str(('obj for serialization is not a "cursor": ', 
                     type(_cursor) )))
        
        i=0
        cols = [col[0] for col in _cursor.description]
        
        logger.info('begin serializing')
        for row in _cursor.fetchall():
            if i!=0:
                raw_data.write(',\n')
            raw_data.write(json.dumps(
                OrderedDict(zip(cols, row)),
                skipkeys=False,ensure_ascii=True,check_circular=True, 
                allow_nan=True, cls=DjangoJSONEncoder))
            i += 1

        logger.info('done, wrote: %d' % i)



class LimsSerializer(PrettyJSONSerializer, BackboneSerializer,CSVSerializer):
    ''' Combine all of the Serializers used by the API
    '''

class SmallMoleculeSerializer(LimsSerializer, SDFSerializer):
    ''' Combine all of the Serializers used by the API
    '''

