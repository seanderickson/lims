import csv
import copy
from collections import OrderedDict
import cStringIO
import StringIO
import json
import logging
import re

from django.utils.encoding import smart_str
from django.core.serializers.json import DjangoJSONEncoder
from django.db.backends.util import CursorDebugWrapper
from psycopg2.psycopg1 import cursor
from tastypie import fields
from tastypie.bundle import Bundle
from tastypie.serializers import Serializer

import reports.utils.sdf2py as s2p
import reports.utils.serialize
import xlwt
import xlrd
from reports.utils.serialize import from_csv_iterate

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


class TextIntegerField(fields.IntegerField):
    """
    Special hydration for text-to-integer interpretation; 
    - convert empty string to None
    """
    dehydrated_type = 'integer'
    help_text = 'Integer data. Ex: 2673'

    def hydrate(self, bundle):
        val = super(TextIntegerField,self).hydrate(bundle)
        if not val or val == '': 
            return None
        return val


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
    
    # FIXME: init could call super.__init__ first, then just modify
    # "self.content_types" and "self.formats"
    def __init__(self, content_types=None, formats=None, **kwargs):

        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['sdf'] = 'chemical/x-mdl-sdfile'
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            formats = _formats
        formats.append('sdf')
            
        super(SDFSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);

        
    def to_sdf(self, data, options=None):
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(str(('to_sdf', data, options)))

        data = self.to_simple(data, options)
        output = cStringIO.StringIO()

        # TODO: smarter way to ignore 'objects'
        if 'objects' in data:
            data = data['objects']
        if len(data) == 0:
            return data
        
        if isinstance(data,dict):
            data = [data]
            
        MOLDATAKEY = s2p.MOLDATAKEY
        for d in data:
            
            if d.get(MOLDATAKEY, None):
                output.write(str(d[MOLDATAKEY]))
                output.write('\n') 
                # because we've not copied the data, don't delete it
                # future optimize: implement data as iterable
                #                 del d[MOLDATAKEY]
            for k,v in d.items():
                if k == MOLDATAKEY: 
                    continue
                output.write('> <%s>\n' % k)
                # according to 
                # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
                # "only one blank line should terminate a data item"
                if v:
                    # find lists, but not strings (or dicts)
                    # Note: a dict here will be non-standard; probably an error 
                    # report, so just stringify dicts as is.
                    if not hasattr(v, "strip") and isinstance(v, (list,tuple)): 
                        for x in v:
                            # DB should be UTF-8, so this should not be necessary,
                            # however, it appears we have legacy non-utf data in 
                            # some tables (i.e. small_molecule_compound_name 193090
                            output.write(unicode.encode(x,'utf-8'))
#                             output.write(str(x))
                            output.write('\n')
                    else:
                        output.write(str(v))
                        output.write('\n')

                output.write('\n')
            output.write('$$$$\n')
        return output.getvalue()
    
    def from_sdf(self, content, root='objects'):
        '''
        @param root - property to nest the return object iterable in for the 
            response (None if no nesting, and return object will be an iterable)

        NOTE: the use of "objects" only makes sense for lists - and cannot be used
        when there is only one item in the list; therefore this will not even work
        with a "POST"; tastypie isn't aware that the object is nested in "objects".
        This shouldn't be used at all, but we need to dig further into tastypie

        TODO: version 2 - read from a stream
        '''
        objects = s2p.parse_sdf(content,
            _delimre=re.compile(ur'(?<=\n)\$\$\$\$'))
        if root and not isinstance(objects, dict):
            return { root: objects }
        else:
            return objects

# class SmilesPNGSerializer(Serializer):
#     
#     def __init__(self, content_types=None, formats=None, **kwargs):
# 
#         if not content_types:
#             content_types = Serializer.content_types.copy();
#         content_types['png'] = 'image/png'
#         
#         if not formats:
#             _formats = Serializer.formats # or []
#             _formats = copy.copy(_formats)
#             formats = _formats
#         formats.append('png')
#             
#         super(SmilesPNGSerializer,self).__init__(
#             formats=formats, 
#             content_types=content_types,**kwargs);
# 
#         
#     def to_png(self, data, options=None):
#         import rdkit.Chem
#         import rdkit.Chem.AllChem
#         import rdkit.Chem.Draw
# #         import matplotlib
# 
# 
#         m = rdkit.Chem.MolFromSmiles('Cc1ccccc1')
#         rdkit.Chem.AllChem.Compute2DCoords(m)
#         im = rdkit.Chem.Draw.MolToImage(m)
#         return im
# #         matplotlib.pyplot.imshow(im)
#         
# #         response = HttpResponse(mimetype="image/png")
# #         img.save(response, "PNG")
# #         return response

# class MultiPartDeserializer(Serializer):
#      
#     def __init__(self, content_types=None, formats=None, **kwargs):
#         
#         if not content_types:
#             content_types = Serializer.content_types.copy();
#         content_types['multipart_form_data'] = 'multipart/form-data'
#         
#         if not formats:
#             _formats = Serializer.formats # or []
#             _formats = copy.copy(_formats)
#             formats = _formats
#         formats.append('multipart_form_data')
#             
#         super(MultiPartDeserializer,self).__init__(
#             formats=formats, 
#             content_types=content_types,**kwargs);
# 
#     def from_multipart_form_data(self, content, **kwargs):
# 
#         logger.info(str(('content', content)))

def csv_convert(val):
    if isinstance(val, (list,tuple)):
        return '[' + ','.join([smart_str(x) for x in val]) + ']' 
    elif val != None:
        if type(val) == bool:
            if val:
                return 'TRUE'
            else:
                return 'FALSE'
        else:
            return smart_str(val)
    else:
        return None

class XLSSerializer(Serializer):
    
    def __init__(self,content_types=None, formats=None, **kwargs):
        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['xls'] = 'application/xls'
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            formats = _formats
        formats.append('xls')
            
        super(XLSSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);


    def to_xls(self, data, options=None):

        options = options or {}
        data = self.to_simple(data, options)
        
        raw_data = StringIO.StringIO()
        book = xlwt.Workbook(encoding='utf8')
        
        if 'error' in data or 'error_messsage' in data:
            sheet = book.add_sheet('error')
            sheet.write(0, 0, 'error')
            sheet.write(1, 0, data.get('error', data.get('error_message', 'unknown error')))
            if data.get('traceback', None):
                sheet.write(0,1, 'traceback')
                sheet.write(1,1, data.get('traceback', ''))
            book.save(raw_data)
            return raw_data.getvalue()
        
        # TODO: smarter way to ignore 'objects'
        if 'objects' in data:
            data = data['objects']
        if len(data) == 0:
            return data

        sheet = book.add_sheet('objects')

        if isinstance(data, dict):
            # usually, this happens when the data is actually an error message;
            # but also, it could be just one item being returned
            for i,(key,item) in enumerate(data.items()):
                sheet.write(0,i,smart_str(key))
                sheet.write(1,i,csv_convert(item))
        else:    
            # default 
            keys = None
            for row,item in enumerate(data):
                if row == 0:
                    for i, key in enumerate(item.keys()):
                        sheet.write(0,i,smart_str(key))
                for i, val in enumerate(item.values()):
                    if val and len(csv_convert(val)) > 32767: 
                        logger.error(str(('warn, row too long', row,key, csv_convert(val))))
                    sheet.write(row+1,i,csv_convert(val))
        
        book.save(raw_data)
        
        return raw_data.getvalue()

    def from_xls(self, content, root='objects'):
        
        wb = xlrd.open_workbook(file_contents=content)
        
        if wb.nsheets > 1:
            logger.warn('only first page of workbooks supported')
        
        # TODO: if root is specified, then get the sheet by name
        sheet = wb.sheet_by_index(0)

        # convert sheet to a flat array
#         rows = []
#         for row in range(sheet.nrows):
#             values = []
#             for col in range(sheet.ncols):
#                 values.append(sheet.cell(row,col).value)
#             rows.append(values)

        def read_sheet(sheet):
            def read_row(row):
                for col in range(sheet.ncols):
                    cell = sheet.cell(row,col)
                    value = cell.value
                    if not value:
                        yield None
                    elif cell.ctype == xlrd.XL_CELL_NUMBER:
                        ival = int(value)
                        if value == ival:
                            value = ival
                        yield str(value)
                    else:
                        yield str(value)
            for row in range(sheet.nrows):
                yield read_row(row)

        # because workbooks are treated like sets of csv sheets, now convert
        # as if this were a csv sheet
        data = from_csv_iterate(read_sheet(sheet))

        if root:
            return { root: data }
        else:
            return data
                
class CSVSerializer(Serializer):
    
    def __init__(self, content_types=None, formats=None, **kwargs):
        
        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['csv'] = 'text/csv'
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            formats = _formats
        formats.append('csv')
            
        super(CSVSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);
        
        
    def to_csv(self, data, root='objects', options=None):
        '''
        @param root ignored for csv!
        
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
        
        # TODO: smarter way to ignore 'objects'
        # TODO: if error is thrown all the way to tp.wrapper, then root may not be a string...
        if 'objects' in data:
            data = data['objects']

        if len(data) == 0:
            return data

        if isinstance(data, dict):
            # usually, this happens when the data is actually an error message;
            # but also, it could be just one item being returned
            logger.error(str(('non-standard data', data)))
            raise Exception(str(('non-standard data', data)))
#             keys = data.keys()
#             writer.writerow([smart_str(key) for key in keys])
#             writer.writerow(self.get_list(data))
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

    def to_csv_stream(self, query, options=None, outputstream=None):
        '''
        @param root ignored for csv!
        
        '''
        
        options = options or {}
        data = self.to_simple(data, options)
        
        # default: delimiter = ',' quotechar='"'
        writer = csv.writer(raw_data) 

        if 'error' in data:
            writer.writerow(['error'])
            writer.writerow([data['error']])
            logger.warn(str(('error', data)))
            return raw_data.getvalue()
            
        i = 0
        keys = None
        for item in query:
            if i == 0:
                keys = item.keys()
                writer.writerow([smart_str(key) for key in keys])
            i += 1
            writer.writerow(self.get_list(item))

        
    def get_list(self,item):
        '''
        Convert a csv row into a list of values
        '''
        _list = []
        for key in item:
            _list.append(csv_convert(item[key]))
        return _list
    
    def from_csv(self, content, root='objects'):
        '''
        @param root - property to nest the return object iterable in for the 
            response (None if no nesting, and return object will be an iterable)

        NOTE: the use of "objects" only makes sense for lists - and cannot be used
        when there is only one item in the list; therefore this will not even work
        with a "POST"; tastypie isn't aware that the object is nested in "objects".
        This shouldn't be used at all, but we need to dig further into tastypie

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


class LimsSerializer(PrettyJSONSerializer, BackboneSerializer,CSVSerializer, 
                        SDFSerializer, XLSSerializer):
    ''' 
    Combine all of the Serializers used by the API
    '''
    

# class SmallMoleculeSerializer(LimsSerializer, SDFSerializer):
#     ''' Combine all of the Serializers used by the API
#     '''

