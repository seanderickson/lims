# Non-streaming implementations of Tastypie Serializer
from __future__ import unicode_literals

import StringIO
import cStringIO
from collections import OrderedDict
import copy
import csv
import json
import logging
import re

from django.core.serializers.json import DjangoJSONEncoder
from django.db.backends.utils import CursorDebugWrapper
from django.http.response import StreamingHttpResponse
from django.utils.encoding import smart_str
import mimeparse
from openpyxl.workbook.workbook import Workbook
from psycopg2.psycopg1 import cursor
import six
from tastypie.exceptions import UnsupportedFormat, BadRequest, \
    ImmediateHttpResponse
from tastypie.serializers import Serializer
import xlrd

from db.support import screen_result_importer
from reports.serialize import XLSX_MIMETYPE, XLS_MIMETYPE, SDF_MIMETYPE
from reports.serialize import dict_to_rows
from reports.serialize.csvutils import LIST_DELIMITER_CSV
import reports.serialize.csvutils as csvutils
import reports.serialize.sdfutils as sdfutils
from reports.serialize.streaming_serializers import generic_xlsx_response
from reports.serialize.xlsutils import LIST_DELIMITER_XLS
import reports.serialize.xlsutils as xlsutils


logger = logging.getLogger(__name__)

class BaseSerializer(Serializer):

    @staticmethod
    def get_content(resp):
        if isinstance(resp, StreamingHttpResponse):
            # Response is in an iterator, 
            # and can can possibly be iterated only once (resultsets)
            # Use this utility to cache the result (e.g. for testing)
            if not hasattr(resp,'cached_content') or not getattr(resp, 'cached_content'):
                buffer = cStringIO.StringIO()
                for line in resp.streaming_content:
                    buffer.write(line)
                resp.cached_content = buffer.getvalue()
            return resp.cached_content
        else:
            return resp.content

    def get_format(self, request, **kwargs):
        '''
        Return mime-type "format" set on the request, or
        use mimeparse.best_match:
        "Return the mime-type with the highest quality ('q') from list of 
        candidates."
        - uses Resource.content_types
        - replaces tastypie.utils.mime.determine_format
        '''
        
        format = None
        if kwargs and 'format' in kwargs:
            format = kwargs['format']
        
        if not format and request.GET.get('format',None):
            format = request.GET.get('format')
            
        if format:
            logger.debug('mapping format: %r', format)
            if format in self.content_types:
                format = self.content_types[format]
                logger.debug('format: %r', format)
            else:
                logger.error('unknown format: %r, options: %r', 
                    format, self.content_types)
                raise ImmediateHttpResponse("unknown format: %s" % format)
                
        else:
            logger.info('get_format: Try to fallback on the Accepts header')
            if request.META.get('HTTP_ACCEPT', '*/*') != '*/*':
                try:
                    format = mimeparse.best_match(
                        self.content_types.values(), 
                        request.META['HTTP_ACCEPT'])
                    if format == 'text/javascript':
                        # NOTE - 
                        # if the HTTP_ACCEPT header contains multiple entries 
                        # with equal weighting, mimeparse.best_match returns
                        # the last match. This results in the request header:
                        # "application/json, text/javascript, */*; q=0.01"
                        # (sent from jquery ajax call) returning 'text/javascript'
                        # because the tastypie wrapper interprets this as 
                        # a JSONP request, override here and set to 
                        # 'application/json' 
                        if 'application/json' in  request.META['HTTP_ACCEPT']:
                            format = 'application/json'

                    if not format:
                        raise ImmediateHttpResponse(
                            "no best match format for HTTP_ACCEPT: %r"
                            % request.META['HTTP_ACCEPT'])
                        
                    logger.info('format %s, HTTP_ACCEPT: %s', format, 
                        request.META['HTTP_ACCEPT'])
                except ValueError:
                    logger.exception('Invalid Accept header: %r',
                        request.META['HTTP_ACCEPT'])
                    raise ImmediateHttpResponse('Invalid Accept header')
            elif request.META.get('CONTENT_TYPE', '*/*') != '*/*':
                format = request.META.get('CONTENT_TYPE', '*/*')

        logger.debug('got format: %s', format)
        return format

    def serialize(self, data, format='application/json', options=None):
        """
        Override Tastypie - 
        """

        desired_format = None
        if options is None:
            options = {}

        for short_format, long_format in self.content_types.items():
            if format == long_format:
                if hasattr(self, "to_%s" % short_format):
                    desired_format = short_format
                    break

        if desired_format is None:
            raise UnsupportedFormat(
                "The format indicated '%s' had no available serialization "
                "method. Please check your ``formats`` and ``content_types`` "
                "on your Serializer." % format)

        serialized = getattr(self, "to_%s" % desired_format)(data, options)
        return serialized

    def deserialize(self, request=None, content=None, format='application/json'):
        """
        Override Tastypie - 
        - pass the request, to allow multipart form processing in subclasses
        """
        assert (request or content), ('must specify either request or content')

        desired_format = None

        format = format.split(';')[0]

        for short_format, long_format in self.content_types.items():
            if format == long_format:
                if hasattr(self, "from_%s" % short_format):
                    desired_format = short_format
                    break

        if desired_format is None:
            raise UnsupportedFormat(
                ("The format indicated '%s' had no available deserialization "
                 "method. Please check your ``formats`` and ``content_types`` "
                 "on your Serializer.") % format)

        if not content:
            content = request.body
        
        if isinstance(content, six.binary_type):
            content = force_text(content)

        deserialized = getattr(self, "from_%s" % desired_format)(content)
        return deserialized

class BackboneSerializer(BaseSerializer):
    
    def from_json(self, content):
        """
        Override to quote attributes from the client.
        """
        content = content.decode('utf-8').replace(r'(\w+):', r'"\1" :')
        if content:
            return json.loads(content)
        else:
            return None

class PrettyJSONSerializer(BaseSerializer):
    json_indent = 2

    def to_json(self, data, options=None):
        options = options or {}
        data = self.to_simple(data, options)
        # NOTE, using "ensure_ascii" = True to force encoding of all 
        # chars to the ascii charset; otherwise, cStringIO has problems
        return json.dumps(
            data, cls=DjangoJSONEncoder,
            sort_keys=True, ensure_ascii=True, indent=self.json_indent, 
            encoding="utf-8")


class SDFSerializer(BaseSerializer):
    
    def __init__(self, content_types=None, formats=None, **kwargs):

        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['sdf'] = SDF_MIMETYPE
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            formats = _formats
        formats.append('sdf')
            
        super(SDFSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);
        
    def to_sdf(self, data, options=None):
        
        data = self.to_simple(data, options)

        if 'objects' in data:
            data = data['objects']
        if len(data) == 0:
            return data
        
        output = cStringIO.StringIO()
        sdfutils.to_sdf(data,output)
        return output.getvalue()
    
    def from_sdf(self, content, root='objects'):
        '''
        @param root - property to nest the return object iterable in for the 
            response (None if no nesting, and return object will be an iterable)

        '''
        objects = sdfutils.parse_sdf(content)
        if root and not isinstance(objects, dict):
            return { root: objects }
        else:
            return objects


class XLSSerializer(BaseSerializer):
    
    def __init__(self,content_types=None, formats=None, **kwargs):
        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['xls'] = XLS_MIMETYPE
        content_types['xlsx'] = XLSX_MIMETYPE
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            formats = _formats
        formats.append('xls')
        formats.append('xlsx')
            
        super(XLSSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);


    def deserialize(self, request=None, content=None, format='application/json'):
        """
        Override - 20150304 - to remove "force_text"
        """
        assert (request or content), ('must specify either request or content')

        if not format:
            format = request.META.get('CONTENT_TYPE', 'application/json')

        desired_format = None

        format = format.split(';')[0]

        for short_format, long_format in self.content_types.items():
            if format == long_format:
                if hasattr(self, "from_%s" % short_format):
                    desired_format = short_format
                    break

        if desired_format is None:
            raise UnsupportedFormat(
                ("The format indicated '%s' had no available deserialization "
                "method. Please check your ``formats`` and ``content_types`` "
                "on your Serializer.") % format)
        
        # override for xls
        # if isinstance(content, six.binary_type):
        #     content = force_text(content)

        if not content:
            content = request.body

        deserialized = getattr(self, "from_%s" % desired_format)(content)
        return deserialized
    
    def to_xls(self,data, options=None):
        # ALL XLS serialization is to xlsx
        return self.to_xlsx(data, options)
    
    def to_xlsx(self, data, options=None):
        logger.info('Non-streamed xlsx data using generic serialization ')
        
        def sheet_rows(list_of_objects):
            ''' write a header row using the object keys '''
            for i,item in enumerate(list_of_objects):
                if i == 0:
                    yield item.keys()
                yield item.values()
        if 'objects' in data:
            data['objects'] = sheet_rows(data['objects'])
        else:
            data = { 'objects': sheet_rows(data) }
        
        response = generic_xlsx_response(data)
        return self.get_content(response)
        
    def from_xlsx(self, content, root='objects'):
        return self.from_xls(content, root=root)

    def from_xls(self, content, root='objects'):
        
        if isinstance(content, six.string_types):
            wb = xlrd.open_workbook(file_contents=content)
        else:
            wb = xlrd.open_workbook(cStringIO.StringIO(content))
            
        if wb.nsheets > 1:
            logger.warn('only first page of workbooks supported')
        
        # TODO: if root is specified, then get the sheet by name
        sheet = wb.sheet_by_index(0)
        
        if sheet.name.lower() in ['error', 'errors']:
            return xlsutils.sheet_rows(sheet)
            
        # because workbooks are treated like sets of csv sheets, now convert
        # as if this were a csv sheet
        data = csvutils.from_csv_iterate(
            xlsutils.sheet_rows(sheet),list_delimiter=LIST_DELIMITER_XLS)

        if root:
            return { root: data }
        else:
            return data

#     def to_xlsx_old(self, data, options=None):
#         '''
#         TODO: replace xlsxwriter.Workbook with openpyxl.Workbook
#         - 20160411
#         NOTE: not used by the API - streaming_serializers used instead
#         '''
#         
#         options = options or {}
#         data = self.to_simple(data, options)
#         
#         raw_data = StringIO.StringIO()
#         workbook = Workbook(raw_data)
#         
#         if 'error' in data or 'error_messsage' in data:
#             sheet = workbook.add_worksheet('error')
#             sheet.write(0, 0, 'error')
#             sheet.write(1, 0, data.get('error', data.get('error_message', 'unknown error')))
#             if data.get('traceback', None):
#                 sheet.write(0,1, 'traceback')
#                 sheet.write(1,1, data.get('traceback', ''))
#             book.close()
#             return raw_data.getvalue()
#         
#         if 'objects' in data:
#             data = data['objects']
#         if len(data) == 0:
#             return data
# 
#         sheet = workbook.add_worksheet('objects')
# 
#         if isinstance(data, dict):
#             # usually, this happens when the data is actually an error message;
#             # but also, it could be just one item being returned
#             for i,(key,item) in enumerate(data.items()):
#                 sheet.write(0,i,smart_str(key))
#                 sheet.write(1,i,csvutils.csv_convert(item))
#         else:    
#             # default 
#             keys = None
#             for row,item in enumerate(data):
#                 if row == 0:
#                     for i, key in enumerate(item.keys()):
#                         sheet.write(0,i,smart_str(key))
#                 for i, val in enumerate(item.values()):
#                     val = csvutils.csv_convert(val, delimiter=LIST_DELIMITER_XLS)
#                     if val and len(val) > 32767: 
#                         logger.error('warn, row too long, %d, key: %r, len: %d', 
#                             row,key,len(val) )
#                     sheet.write(row+1,i,val)
#         
#         workbook.close()
#         
#         return raw_data.getvalue()
# 
#     def to_xls_old(self, data, options=None):
# 
#         options = options or {}
#         data = self.to_simple(data, options)
#         
#         raw_data = StringIO.StringIO()
#         book = xlwt.Workbook(encoding='utf8')
#         
#         if 'error' in data or 'error_messsage' in data:
#             sheet = book.add_sheet('error')
#             sheet.write(0, 0, 'error')
#             sheet.write(1, 0, data.get(
#                 'error', data.get('error_message', 'unknown error')))
#             if data.get('traceback', None):
#                 sheet.write(0,1, 'traceback')
#                 sheet.write(1,1, data.get('traceback', ''))
#             book.save(raw_data)
#             return raw_data.getvalue()
#         
#         if 'objects' in data:
#             data = data['objects']
#         if len(data) == 0:
#             return data
# 
#         sheet = book.add_sheet('objects')
# 
#         if isinstance(data, dict):
#             # usually, this happens when the data is actually an error message;
#             # but also, it could be just one item being returned
#             for i,(key,item) in enumerate(data.items()):
#                 sheet.write(0,i,smart_str(key))
#                 sheet.write(1,i,csvutils.csv_convert(item))
#         else:    
#             # default 
#             keys = None
#             for row,item in enumerate(data):
#                 if row == 0:
#                     for i, key in enumerate(item.keys()):
#                         sheet.write(0,i,smart_str(key))
#                 for i, val in enumerate(item.values()):
#                     val = csvutils.csv_convert(val, delimiter=LIST_DELIMITER_XLS)
#                     if val and len(val) > 32767: 
#                         logger.error('warn, row too long, %d, key: %r, len: %d', 
#                             row,key,len(val) )
#                     sheet.write(row+1,i,val)
#         
#         book.save(raw_data)
#         
#         return raw_data.getvalue()



class ScreenResultSerializer(XLSSerializer):

    def __init__(self,content_types=None, formats=None, **kwargs):
        if not content_types:
            content_types = Serializer.content_types.copy();
        content_types['xls'] = 'application/xls'
        content_types['xlsx'] = XLSX_MIMETYPE
        
        if not formats:
            _formats = Serializer.formats # or []
            _formats = copy.copy(_formats)
            formats = _formats
        formats.append('xls')
        formats.append('xlsx')
            
        super(ScreenResultSerializer,self).__init__(
            formats=formats, 
            content_types=content_types,**kwargs);

    def get_format(self, request, **kwargs):
        format = super(ScreenResultSerializer,self).get_format(request,**kwargs)
        if format.startswith('multipart'):
            logger.info('request.Files.keys: %r', request.FILES.keys())
            if len(request.FILES.keys()) != 1:
                raise ImmediateHttpResponse(
                    response=self.error_response(
                        request,
                        { 'FILES': 'File upload supports only one file at a time',
                          'filenames': request.FILES.keys(),
                        }
                ))
            if 'xls' in request.FILES:
                format = 'application/xls'
            else:
                logger.error(
                    'Unsupported multipart file key: %r',request.FILES.keys())
                raise UnsupportedFormat(
                    'Unsupported multipart file key: %r' % request.FILES.keys())

        logger.info('screen result format determined: %r', format)
        return format 
    
    def deserialize(self, request=None, content=None, format='application/json'):
        '''
        Override deserialize:
        - pass in request
        - get content from the multipart form
        '''
        assert (request or content), ('must specify either request or content')
        logger.info('format: %r', format)
        if not format:
            format = request.META.get('CONTENT_TYPE', 'application/json')
        logger.info('deserialize: %r', format)
        if format.startswith('multipart'):
            logger.info('multipart, files: %r', dir(request))
            logger.info('request.Files.keys: %r', request.FILES.keys())
            if len(request.FILES.keys()) != 1:
                raise ImmediateHttpResponse(
                    response=self.error_response(
                        request,
                        { 'FILES': 'File upload supports only one file at a time',
                          'filenames': request.FILES.keys(),
                        }
                ))
            if 'xls' in request.FILES:
                file = request.FILES['xls']
                deserialized = screen_result_importer.read_file(file)
            else:
                logger.error(
                    'Unsupported multipart file key: %r',request.FILES.keys())
                raise UnsupportedFormat(
                    'Unsupported multipart file key: %r' % request.FILES.keys())
        else:
            desired_format = None
            format = format.split(';')[0]
            for short_format, long_format in self.content_types.items():
                if format == long_format:
                    if hasattr(self, "from_%s" % short_format):
                        desired_format = short_format
                        break
            if not content:
                content = request.body
            if desired_format is None:
                raise UnsupportedFormat(
                    ("The format indicated '%s' had no available deserialization "
                    "method. Please check your ``formats`` and ``content_types`` "
                    "on your Serializer.") % format)

            try:
                deserialized = getattr(self, "from_%s" % desired_format)(content)
            except Exception, e:
                logger.exception('on deserialize: %r, %r', desired_format, content)
                raise
        if not deserialized:
            raise BadRequest('unknown format')
        return deserialized

    def to_xlsx(self, data, options=None):
        logger.info('Non-streamed Screenresult using generic serialization')
        response = generic_xlsx_response(data)
        return self.get_content(response)
    
    def to_xls(self, data, options=None):
        return self.to_xlsx(data, options)
    
    def from_xlsx(self, content):
        return self.from_xls(content)

    def from_xls(self, content):
        
        if isinstance(content, six.string_types):
            wb = xlrd.open_workbook(file_contents=content)
        else:
            wb = xlrd.open_workbook(cStringIO.StringIO(content))
        return screen_result_importer.read_workbook(wb)

                
class CSVSerializer(BaseSerializer):
    
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
            for row in dict_to_rows(data['error']):
                writer.writerow(row)
            # writer.writerow(['error'])
            # writer.writerow([data['error']])
            # logger.warn(str(('error', data)))
            return raw_data.getvalue()
            
        if 'objects' in data:
            data = data['objects']

        if len(data) == 0:
            return data

        if isinstance(data, dict):
            # usually, this happens when the data is actually an error message;
            # but also, it could be just one item being returned
            raise Exception('non-standard data: embedded dict: %r', data)
        else:    
            i = 0
            keys = None
            for item in data:
                if i == 0:
                    keys = item.keys()
                    writer.writerow([smart_str(key) for key in keys])
                i += 1
                writer.writerow([csvutils.csv_convert(val) for val in item.values()])

        return raw_data.getvalue()

    def from_csv(self, content, root='objects'):
        '''
        @param root - property to nest the return object iterable in for the 
            response (None if no nesting, and return object will be an iterable)

        '''
        objects = csvutils.from_csv(
            StringIO.StringIO(content),list_delimiter=LIST_DELIMITER_CSV)
        if root:
            return { root: objects }
        else:
            return objects
                

class CursorSerializer(BaseSerializer):
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
  
    def to_csv(self, obj, options=None):
        raw_data = StringIO.StringIO()

        # this is an unexpected way to get this error, look into tastypie
        if(isinstance(obj,dict) and 'error_message' in obj):
        
            logger.warn('report error: %r', obj)
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
                        # NOTE, using "ensure_ascii" = True to force encoding of all 
                        # chars to the ascii charset; otherwise, cStringIO has problems
                        writer.writewrow([key,json.dumps( 
                            value,skipkeys=False,check_circular=True,ensure_ascii=True, 
                            allow_nan=True, default=lambda x: str(x), encoding="utf-8")] )
                    else:
                        logger.warn(
                            'non-cursor data will not be written to csv: "' + 
                            key +'": ' + json.dumps(
                                value,skipkeys=False,check_circular=True,ensure_ascii=True,
                                allow_nan=True, default=lambda x: str(x)), encoding="utf-8")
            
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

    def to_json(self,obj, options=None):
        
        raw_data = StringIO.StringIO()
         
        if isinstance(obj,dict) and 'error_message' in obj :
            logger.warn('report error: %r', obj)
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
                    # NOTE, using "ensure_ascii" = True to force encoding of all 
                    # chars to the ascii charset; otherwise, cStringIO has problems
                    raw_data.write(json.dumps(
                        value,skipkeys=False,check_circular=True,ensure_ascii=True,
                        allow_nan=True, default=lambda x: str(x), encoding="utf-8"))
                count += 1
                
            raw_data.write('}')
                    
        return raw_data.getvalue() 
        
    def _cursor_to_json(self, _cursor, raw_data):
        if not isinstance(_cursor, (cursor, CursorDebugWrapper) ):
            raise Exception(
                'obj for serialization is not a "cursor": %r'
                % type(_cursor) )
        
        i=0
        cols = [col[0] for col in _cursor.description]
        
        logger.info('begin serializing')
        for row in _cursor.fetchall():
            if i!=0:
                raw_data.write(',\n')
            # NOTE, using "ensure_ascii" = True to force encoding of all 
            # chars to the ascii charset; otherwise, cStringIO has problems
            raw_data.write(json.dumps(
                OrderedDict(zip(cols, row)),
                skipkeys=False,ensure_ascii=True,check_circular=True,
                allow_nan=True, cls=DjangoJSONEncoder,encoding="utf-8"))
            i += 1

        logger.info('done, wrote: %d' % i)

class LimsSerializer(PrettyJSONSerializer, BackboneSerializer,CSVSerializer, 
                        SDFSerializer, XLSSerializer):
    ''' 
    Combine all of the Serializers used by the API
    '''


