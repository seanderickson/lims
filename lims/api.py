
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


class CSVSerializer(BackboneSerializer):
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
        obj_list.query.add_ordering(*temp)
        
#        logger.info(str(('obj_list.query', obj_list.query.as_sql())))
        return obj_list

