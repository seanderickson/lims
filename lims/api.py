
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer
from django.utils.encoding import smart_str
from copy import deepcopy

import csv
import StringIO
import json
import logging
        
from django.core.exceptions import ObjectDoesNotExist
from tastypie.exceptions import NotFound
from django.utils.timezone import is_naive
from collections import OrderedDict
from django.core.serializers.json import DjangoJSONEncoder
from tastypie.bundle import Bundle
from psycopg2.psycopg1 import cursor
from django.db.backends.util import CursorDebugWrapper


logger = logging.getLogger(__name__)


class BackboneSerializer(Serializer):
    
    def from_json(self, content):
        """
        Given some JSON data, returns a Python dictionary of the decoded data.
        Override to quote attributes - the backbone client doesn't want to do this.
        """
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(("loading content:", content)))
        content = content.replace(r'(\w+):', r'"\1" :')
        if(logger.isEnabledFor(logging.DEBUG)):
            logger.debug(str(("loading content:", content)))
        return json.loads(content)


class TimeZoneAwareDateSerializer(Serializer):
    """
    Our own serializer to format datetimes in ISO 8601 but with timezone
    offset.
    credit: http://www.tryolabs.com/Blog/2013/03/16/displaying-timezone-aware-dates-tastypie/
    """
    def format_datetime(self, data):
#        logger.info(str(('formatting date', data)))
        # If naive or rfc-2822, default behavior...
        if is_naive(data) or self.datetime_formatting == 'rfc-2822':
            return super(TimeZoneAwareDateSerializer, self).format_datetime(data)
 
        return data.isoformat()
 

class CSVSerializer(Serializer):
    formats = ['json', 'jsonp', 'xml', 'yaml', 'html', 'plist', 'csv']
    content_types = {
        'json': 'application/json',
        'jsonp': 'text/javascript',
        'xml': 'application/xml',
        'yaml': 'text/yaml',
        'html': 'text/html',
        'plist': 'application/x-plist',
        'csv': 'text/csv',
    }

    def to_csv(self, data, options=None):
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
        _list = []
        for key in item:
            logger.info(str(('item', item)))
            if item[key] and isinstance(item[key], (list, tuple)):
                _list.append( '[' + ','.join([smart_str(x) for x in item[key]]) + ']' )
            elif item[key] != None:
                _list.append(smart_str(item[key]))
            else:
                _list.append(None)
        return _list
    
    def from_csv(self, content):
        raw_data = StringIO.StringIO(content)
        data = { 'objects': [] }
        # TODO: also, stream this
        # default: delimiter = ',' quotechar='"'
        logger.info('reading...')
        reader = csv.reader(raw_data)
        
        i = 0
        keys = []
        list_keys = [] # cache
        for row in reader:
            if i == 0:
                keys = row
            else:
                item = dict(zip(keys,row))
                logger.debug(str(('read row', item)))
                for key in item.keys():
                    val = item[key]
                    if val and len(val)> 1 and (key in list_keys or val[0] == '['):
                        # due to the simplicity of the serializer, above, any quoted string is a nested list
                        list_keys.append(key)
                        item[key] = val.strip('"[]').split(',')
                        logger.debug(str(('converted',val,item[key])))
                data['objects'].append(item)
            i += 1
                
        return data



class CursorSerializer(Serializer):
    """
    A simple serializer that takes a cursor, queries it for its columns, and outputs 
    this as either CSV or JSON.
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

#    def get_saf_columns(self, query):
#        return ['one','two', 'three']
    
    def to_csv(self, bundle_or_obj, options=None):
        '''
        NOTE: ignores all content except for "objects"
        '''
        
        logger.info(str(('typeof the object sent to_csv',type(cursor))))
        raw_data = StringIO.StringIO()

        obj = bundle_or_obj            
        if isinstance(bundle_or_obj, Bundle):
            obj = bundle_or_obj.obj
            
        # this is an unexpected way to get this error, look into tastypie sequence calls
        if(isinstance(obj,dict) and 'error_message' in obj):
            logger.info(str(('report error', obj)))
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
                        writer.writewrow([key,json.dumps(value,skipkeys=False,check_circular=True, allow_nan=True, default=lambda x: str(x))] )
                    else:
                        logger.info('non-cursor data will not be written to csv: "' + key 
                                    + '": ' + json.dumps(value,skipkeys=False,check_circular=True, allow_nan=True, default=lambda x: str(x)))
            
        return raw_data.getvalue()
    
    def _cursor_to_csv(self, cursor, writer):
        # no response header needed?
        #        response = HttpResponse(mimetype='text/csv')
        #        response['Content-Disposition'] = 'attachment; filename=%s.csv' % unicode('test.output').replace('.', '_')
        #        raw_data.write(response)
        i=0
        cols = [col[0] for col in cursor.description]
        
        # TODO: grab the column names here
        writer.writerow(cols)

        for row in cursor.fetchall():
            writer.writerow([smart_str(val, 'utf-8', errors='ignore') for val in row])
            i += 1
        logger.info('_cursor_to_csv done, wrote: %d' % i)

    def to_json(self,bundle_or_obj, options=None):
        
        logger.info(str(('typeof the object sent to_json',type(bundle_or_obj))))
#        logger.info(str(('to_csv for SAF for cursor', cursor)))
        raw_data = StringIO.StringIO()
         
        if isinstance(bundle_or_obj, Bundle):
            obj = bundle_or_obj.obj
        else:
            obj = bundle_or_obj            
            
        # this is an unexpected way to get this error, look into tastypie sequence calls
        if isinstance(obj,dict) and 'error_message' in obj :
            logger.info(str(('report error', obj)))
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
                    raw_data.write(json.dumps(value,skipkeys=False,check_circular=True, allow_nan=True, default=lambda x: str(x)))
                count += 1
                
            raw_data.write('}')
            
            # then in this case, this is a non-error dict, probably for the schema, dump and return.
#            raw_data.writelines(json.dumps(obj))
#            return raw_data.getvalue()
#            return json.dumps(obj,skipkeys=False,check_circular=True, allow_nan=True, default=lambda x: str(x))
        
        return raw_data.getvalue() # TODO: how to stream entire set
        
    def _cursor_to_json(self, _cursor, raw_data):
        if not isinstance(_cursor, (cursor, CursorDebugWrapper) ):
            raise Exception(unicode(('obj for serialization is not a "cursor": ', type(_cursor) )))
        
        i=0
        cols = [col[0] for col in _cursor.description]
        
        logger.info('begin serializing')
        for row in _cursor.fetchall():
#            logger.info(unicode((row)))
            if i!=0:
                raw_data.write(',\n')
            raw_data.write(json.dumps(OrderedDict(zip(cols, row)),skipkeys=False,ensure_ascii=True,check_circular=True, allow_nan=True, cls=DjangoJSONEncoder))
            #            raw_data.write(json.dumps(dict(zip(cols, row))))
            i += 1

        logger.info('done, wrote: %d' % i)


class LimsSerializer(BackboneSerializer,TimeZoneAwareDateSerializer,CSVSerializer):
    ''' Combine all of the Serializers used by the API
    '''
    
    

# TODO: this class should be constructed as a Mixin, not inheritor of ModelResource
class PostgresSortingResource(ModelResource):

    def __init__(self, **kwargs):
        super(PostgresSortingResource,self).__init__( **kwargs)

    def apply_sorting(self, obj_list, options):
        """
        Create a none-too-pretty workaround for the postgresql null sorting issue - nulls sort higher than values, 
        which is not desired.  We want nulls to sort lower than values.
        """ 
        
        obj_list = super(PostgresSortingResource, self).apply_sorting(obj_list, options)
        logger.info(str(('order_by', obj_list.query.order_by)))
        extra_select = {}
        extra_ordering = []
        for field in obj_list.query.order_by:
            is_null_dir = '-'  # default nulls first for ascending
            if field.startswith('-'):
                is_null_dir = ''
                field = field[1:]
            extra_select[field+"_null"]=field + ' is null'
            extra_ordering.append(is_null_dir + field+"_null")
        logger.info(str(('extra_select', extra_select, 'extra_ordering', extra_ordering)))
        obj_list = obj_list.extra(extra_select)

        # Note that the following doesn't work, something in the framework deletes the extra 
        # order_by clause when apply_sorting, or, if this is run last, it deletes the sorting applied in apply_sorting...
        #        obj_list = obj_list.extra(order_by=['-comments_null'])

        # note: this doesn't work because the "is null" field order by clauses
        # must be prepended so that they occur before their intended fields
        #        obj_list.query.add_ordering('comments_null')
        
        temp = obj_list.query.order_by;
        obj_list.query.clear_ordering()
        for xfield in extra_ordering:
            temp.insert(0,xfield)
        logger.info(unicode(('ordering', temp)))
        obj_list.query.add_ordering(*temp)
        
#        logger.info(str(('obj_list.query', obj_list.query.as_sql())))
        return obj_list

