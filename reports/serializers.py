# Non-streaming implementations of Tastypie Serializer
from __future__ import unicode_literals

import cStringIO
from collections import OrderedDict
import csv
import json
import logging

from django.conf import settings
from django.db.backends.utils import CursorDebugWrapper
from django.http.response import StreamingHttpResponse
from django.utils.encoding import smart_text, force_text
import mimeparse
import six
from tastypie.exceptions import BadRequest
import xlrd

from db.support import screen_result_importer
from reports.serialize import XLSX_MIMETYPE, XLS_MIMETYPE, SDF_MIMETYPE, \
    JSON_MIMETYPE, CSV_MIMETYPE, to_simple, LimsJSONEncoder
from reports.serialize.csvutils import LIST_DELIMITER_CSV, dict_to_rows
import reports.serialize.csvutils as csvutils
import reports.serialize.sdfutils as sdfutils
from reports.serialize.streaming_serializers import generic_xlsx_response, \
    get_xls_response
from reports.serialize.xlsutils import LIST_DELIMITER_XLS
import reports.serialize.xlsutils as xlsutils


logger = logging.getLogger(__name__)

    
class BaseSerializer(object):
    
    content_types = {'json': JSON_MIMETYPE,
                     'html': 'text/html',
                     'csv': CSV_MIMETYPE,
                     'xls': XLS_MIMETYPE,
                     'xlsx': XLSX_MIMETYPE,
                     'sdf': SDF_MIMETYPE }

    def __init__(self, content_types=None):

        self.content_types = content_types or {}    
        content_types['json'] = JSON_MIMETYPE

    @staticmethod
    def get_content(resp):
        ''' 
        Returns the content of the response:
        - for StreamingResponse, the response may be an iterator,
        and can be iterated only once (resultsets), this method will 
        cache the result for subsequent accesses (e.g. for testing).
        '''
        if isinstance(resp, StreamingHttpResponse):
            if (not hasattr(resp,'cached_content') 
                    or not getattr(resp, 'cached_content')):
                buffer = cStringIO.StringIO()
                for line in resp.streaming_content:
                    buffer.write(line)
                resp.cached_content = buffer.getvalue()
            return resp.cached_content
        else:
            return resp.content

    def get_content_type_for_format(self, format):

        if format in self.content_types:
            return self.content_types[format]
        else:
            logger.warn(
                'unknown format: %r, options: %r',
                format, self.content_types.keys())
        return None
    
    def get_format_for_content_type(self,content_type):
        
        # pull the charset off
        # e.g. application/vnd.openxmlformats-officedocument.spreadsheetml.sheet; charset=utf-8
        content_type = content_type.split(';')[0]
        
        desired_format = None
        for short_format, long_format in self.content_types.items():
            if content_type == long_format:
                if hasattr(self, "to_%s" % short_format):
                    desired_format = short_format
                    break
        return desired_format
    
    def get_accept_content_type(self, request, format=None):
        
        logger.debug('get_accept_content_type: %r, %r', request, format)

        content_type = None
        
        if format is None:
            if request.GET and request.GET.get('format',None):
                format = request.GET.get('format')
        if format is not None:
            content_type = self.get_content_type_for_format(format)
        
        if content_type is None and request is not None:
            if request.META and request.META.get('HTTP_ACCEPT', '*/*') != '*/*':
                try:
                    content_type = mimeparse.best_match(
                        self.content_types.values(), 
                        request.META['HTTP_ACCEPT'])
                    if content_type == 'text/javascript':
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
                            content_type = 'application/json'
                    
                    logger.info('"HTTP_ACCEPT" - content_type: %r', content_type)
    
                    if not content_type:
                        raise BadRequest(
                            "no best match format for HTTP_ACCEPT: %r"
                            % request.META['HTTP_ACCEPT'])
                except ValueError:
                    raise BadRequest('Invalid Accept header: %r',
                        request.META['HTTP_ACCEPT'])
            elif request.META and request.META.get('CONTENT_TYPE', '*/*') != '*/*':
                content_type = request.META.get('CONTENT_TYPE', '*/*')
                logger.info('fallback to "CONTENT_TYPE": %r', content_type)
            else:
                raise BadRequest(
                    'no CONTENT_TYPE or HTTP_ACCEPT header found: %r, %r', 
                    request, format)
        if not content_type:
            logger.info('not content type found')
        return content_type
        
    def get_content_type(self, request, format=None):    

        logger.debug('get_content_type: %r, %r', request, format)

        content_type = None
        
        if format is None:
            if request.GET and request.GET.get('format',None):
                format = request.GET.get('format')
        if format is not None:
            content_type = self.get_content_type_for_format(format)
        
        if content_type is None and request is not None:
            if request.META and request.META.get('CONTENT_TYPE', '*/*') != '*/*':
                content_type = request.META.get('CONTENT_TYPE', '*/*')
                logger.debug('"CONTENT_TYPE": %r', content_type)
            elif request.META and request.META.get('HTTP_ACCEPT', '*/*') != '*/*':
                logger.info(
                    'no "CONTENT_TYPE" found, fallback "HTTP_ACCEPT" header %r',
                    request.META.get('HTTP_ACCEPT', '*/*'))
                content_type = mimeparse.best_match(
                    self.content_types.values(), 
                    request.META['HTTP_ACCEPT'])
                if content_type == 'text/javascript':
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
                        content_type = 'application/json'
                
                logger.debug('"HTTP_ACCEPT" - content_type: %r', content_type)
        if not content_type:
            msg = 'no best match format for CONTENT_TYPE: '
            if request:
                msg += request.META.get('CONTENT_TYPE','-no content type specified')
            else:
                msg += ', request not specified'
            if format:
                msg += ', format: %r' % format
            else:
                msg += ', format not specified. '
            msg += 'ser: %r' % self
            raise BadRequest(msg)
        return content_type

    def serialize(self, data, content_type):
        
        desired_format = self.get_format_for_content_type(content_type)
        if desired_format is None:
            msg = ( 'unknown serialize content_type: %r or format: %r, options: %r'
                    % (content_type, desired_format, self.content_types.values()))
            raise BadRequest(msg)
        serialized = getattr(self, "to_%s" % desired_format)(data)
        return serialized
    
    def deserialize(self, content, content_type, **kwargs):

        desired_format = self.get_format_for_content_type(content_type)
        if desired_format is None:
            msg = ( 'unknown deserialize content_type: %r or format: %r, options: %r'
                    % (content_type, desired_format, self.content_types.values()))
            raise BadRequest(msg)

        if not content:
            return {}
    
        logger.debug('deserializing for %r', desired_format)

        deserialized = getattr(self, "from_%s" % desired_format)(content,**kwargs)
        return deserialized

    def to_html(self, data, options=None):
        ''' For error reporting '''
        
        if isinstance(data, six.string_types):
            return data
        else:
            return '<br>'.join([str(row) for row in dict_to_rows(data)])

    def to_json(self, data, options=None):
        json_indent = 2
        # NOTE, using "ensure_ascii" = True to force encoding of all 
        # chars to the ascii charset; otherwise, cStringIO has problems
        return json.dumps(
            data, cls=LimsJSONEncoder,
            sort_keys=True, ensure_ascii=True, indent=json_indent, 
            encoding="utf-8")

    def from_json(self, content, **kwargs):
        """
        Override to quote attributes from the client.
        """
        if isinstance(content, six.binary_type):
            content = force_text(content)
         
        content = content.decode('utf-8').replace(r'(\w+):', r'"\1" :')
        if content:
            return json.loads(content)
        else:
            return None


class SDFSerializer(BaseSerializer):
    
    def __init__(self, content_types=None):

        content_types = content_types or {}    
        content_types['sdf'] = SDF_MIMETYPE
        
        super(SDFSerializer,self).__init__(content_types=content_types);
        
    def to_sdf(self, data, options=None):
        
        data = to_simple(data)

        if 'objects' in data:
            data = data['objects']
        if len(data) == 0:
            return data
        
        output = cStringIO.StringIO()
        sdfutils.to_sdf(data,output)
        return output.getvalue()
    
    def from_sdf(self, content, root='objects', **kwargs):
        '''
        @param root - property to nest the return object iterable in for the 
            response (None if no nesting, and return object will be an iterable)

        '''
        if isinstance(content, six.binary_type):
            content = force_text(content)
         
        objects = sdfutils.parse_sdf(content)
        if root and not isinstance(objects, dict):
            return { root: objects }
        else:
            return objects


class XLSSerializer(BaseSerializer):
    
    def __init__(self,content_types=None):

        content_types = content_types or {}    
        content_types['xls'] = XLS_MIMETYPE
        content_types['xlsx'] = XLSX_MIMETYPE
        
        super(XLSSerializer,self).__init__(content_types=content_types);

    def to_xls(self,data, options=None, **kwargs):

        # Note: all XLS serialization is to xlsx
        return self.to_xlsx(data, options)
    
    def to_xlsx(self, data, options=None):

        logger.info('Non-streamed xlsx data using generic serialization')
        def sheet_rows(list_of_objects):
            ''' write a header row using the object keys '''
            for i,item in enumerate(list_of_objects):
                if i == 0:
                    yield item.keys()
                yield item.values()
        if 'objects' in data:
            data['objects'] = sheet_rows(data['objects'])
        elif isinstance(data, (list,tuple)):
            data = { 'objects': sheet_rows(data) }
        
        response = generic_xlsx_response(data)
        return self.get_content(response)
        
    def from_xlsx(self, content, root='objects',**kwargs):

        return self.from_xls(content, root=root, **kwargs)

    def from_xls(self, content, root='objects',**kwargs):
        
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
            xlsutils.sheet_rows(sheet),
            list_delimiter=LIST_DELIMITER_XLS, 
            list_keys=kwargs.get('list_keys', None))
 
        if root:
            return { root: data }
        else:
            return data

                
class CSVSerializer(BaseSerializer):
    
    def __init__(self, content_types=None):
        
        content_types = content_types or {}    
        content_types['csv'] = CSV_MIMETYPE

        super(CSVSerializer,self).__init__(content_types=content_types)
        
    def to_csv(self, data, root='objects', options=None):

        data = to_simple(data)

        raw_data = cStringIO.StringIO()
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
            data = dict_to_rows(data)
            for item in data:
                writer.writerow(item)
        else:
            i = 0
            keys = None
            for item in data:
                if i == 0:
                    keys = item.keys()
                    writer.writerow([smart_text(key) for key in keys])
                i += 1
                writer.writerow([csvutils.csv_convert(val) for val in item.values()])

        return raw_data.getvalue()

    def from_csv(self, content, root='objects', **kwargs):
        '''
        @param root - property to nest the return object iterable in for the 
            response (None if no nesting, and return object will be an iterable)

        '''
        if isinstance(content, six.binary_type):
            content = force_text(content)
         
        objects = csvutils.from_csv(
            cStringIO.StringIO(content),
            list_delimiter=LIST_DELIMITER_CSV,
            list_keys=kwargs.get('list_keys', None))
        if root:
            return { root: objects }
        else:
            return objects

class ScreenResultSerializer(XLSSerializer,SDFSerializer,CSVSerializer):

    def __init__(self, content_types=None):
        
        content_types = content_types or {}    
        content_types['xls'] = XLS_MIMETYPE
        content_types['xlsx'] = XLSX_MIMETYPE
        content_types['json'] = JSON_MIMETYPE
        
        super(ScreenResultSerializer,self).__init__(content_types=content_types);

    def to_xlsx(self, data, options=None):
        logger.debug(
            'serialize Non-streamed Screenresult using generic serialization')
        response = get_xls_response(data, 'generic_file')
        return self.get_content(response)
    
    def to_xls(self, data, options=None):
        return self.to_xlsx(data, options)
    
    def from_xlsx(self, content, **kwargs):
        return self.from_xls(content, **kwargs)

    def from_xls(self, content, **kwargs):
        if isinstance(content, six.string_types):
            wb = xlrd.open_workbook(file_contents=content)
        else:
            wb = xlrd.open_workbook(cStringIO.StringIO(content))
        return screen_result_importer.read_workbook(wb)

    def to_json(self, data, options=None):
        return XLSSerializer.to_json(self, data, options=options)

    
class CursorSerializer(BaseSerializer):
    """
    A simple serializer that takes a cursor, queries it for its columns, and
    outputs this as either CSV or JSON.
    (The CSV output is used with SAF)
    """
    
    def __init__(self, content_types=None):

        content_types = content_types or {}    
        content_types['csv'] = CSV_MIMETYPE
        content_types['json'] = JSON_MIMETYPE
        
        super(CursorSerializer,self).__init__(content_types=content_types);

    def to_csv(self, obj, options=None):
        raw_data = cStringIO.StringIO()

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
                [smart_text(val, 'utf-8', errors='ignore') for val in row])
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
                allow_nan=True, cls=LimsJSONEncoder,encoding="utf-8"))
            i += 1

        logger.info('done, wrote: %d' % i)


class LimsSerializer(CSVSerializer,SDFSerializer, XLSSerializer):
    ''' 
    Combine all of the Serializers used by the API
    '''


