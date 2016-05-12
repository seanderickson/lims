'''
Utilities for streaming sql connection cursors in different serialization formats
'''

from __future__ import unicode_literals

import cStringIO
from collections import OrderedDict
import csv
import io
import json
import logging
import os.path
import re
import shutil
import sys
from tempfile import SpooledTemporaryFile, NamedTemporaryFile
import time
from wsgiref.util import FileWrapper
from zipfile import ZipFile

from PIL import Image
from django.conf import settings
from django.core.serializers.json import DjangoJSONEncoder
from django.core.urlresolvers import resolve
from django.http.response import StreamingHttpResponse
from tastypie.exceptions import ImmediateHttpResponse
# from django.utils.encoding import smart_str
# import openpyxl
# using XlsxWriter for constant memory usage
import xlsxwriter

from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    LIST_BRACKETS, MAX_IMAGE_ROWS_PER_XLS_FILE, MAX_ROWS_PER_XLS_FILE, \
    CSV_DELIMITER, HTTP_PARAM_RAW_LISTS


import reports.serialize.csvutils as csvutils
import reports.serialize.sdfutils as sdfutils
from reports.serialize.csvutils import LIST_DELIMITER_CSV
from reports.serialize.xlsutils import generic_xls_write_workbook, screenresult_xls_write_workbook, LIST_DELIMITER_XLS
from reports.serialize import XLSX_MIMETYPE
from reports.serialize import dict_to_rows
from db.support.data_converter import default_converter
from django.utils.timezone import localtime
logger = logging.getLogger(__name__)

MOLDATAKEY = sdfutils.MOLDATAKEY

DEBUG_STREAMING = False or logger.isEnabledFor(logging.DEBUG)

class ChunkIterWrapper(object):
    ''' 
    Iterate in "chunks" of chunk_size chars.
    '''
    def __init__(self, iter_stream, chunk_size = 1024**2):
        self.iter_stream = iter_stream
        self.chunk_size = chunk_size
        self.fragment = None
    
    def next(self):
        logger.debug('chunking')
        bytes = cStringIO.StringIO()
        
        bytecount = 0
        try:
            if self.fragment:
                #  FIXME: 
                # if fragment (row remainder) is still > chunk_size, will just serve it here.
                bytecount = len(self.fragment)
                bytes.write(self.fragment)
                self.fragment = None

            while bytecount < self.chunk_size:
                row = self.iter_stream.next()
                rowlen = len(row)
                if bytecount + rowlen < self.chunk_size:
                    bytes.write(row)
                    bytecount += rowlen
                else:
                    nextbytes = self.chunk_size-bytecount
                    bytes.write(row[:nextbytes])
                    self.fragment = row[nextbytes:]
                    bytecount += nextbytes
        except StopIteration:
            if self.fragment:
                bytes.write(self.fragment)
                self.fragment = None
            else:
                raise StopIteration
        except Exception, e:
            logger.exception('streaming exception')
            raise e   
        finally:
            if bytes.getvalue():
                logger.debug('chunk size %s chars' % len(bytes.getvalue()))
                return bytes.getvalue()
            else:
                logger.debug('normal stop iteration, fragment: %s, bytes: %s' 
                    % ( self.fragment, bytes.getvalue()))
                raise StopIteration

    def __iter__(self):
        return self


def interpolate_value_template(value_template, row):
    ''' 
    Utility function for transforming cell values:
    a "value_template" is of the form:
    "text... {field_name} ..text"
    wherein the {field_name} is replaced with the field-value
    '''                
    def get_value_from_template(matchobj):
        val = matchobj.group()
        val = re.sub(r'[{}]+', '', val)
        if DEBUG_STREAMING:
            logger.info('val from value template: %r, %r, %r',
                val,row.has_key(val), row[val])
        if row.has_key(val):
            return str(row[val])
        else:
            logger.error(
                'field needed for value template %r is not available: %r', 
                val, row)
            return ''
    return re.sub(r'{([^}]+)}', get_value_from_template, value_template)


def json_generator(cursor,meta,request, is_for_detail=False,field_hash=None):
    if DEBUG_STREAMING: logger.info('meta: %r', meta )
    
    # NOTE, using "ensure_ascii" = True to force encoding of all 
    # chars to be encoded using ascii or escaped unicode; 
    # because some chars in db might be non-UTF8
    # and downstream programs have trouble with mixed encoding (cStringIO)
    if not is_for_detail:
        yield ( '{ "meta": %s, "objects": [' 
            % json.dumps(meta, ensure_ascii=True, encoding="utf-8"))
    i=0
    try:
        for row in cursor:
            if DEBUG_STREAMING: logger.info('row: %r', row)
            if field_hash:
                _dict = OrderedDict()
                for key, field in field_hash.iteritems():
                    value = None
                    if row.has_key(key):
                        value = row[key]
                    if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                         or field.get('linked_field_type',None) == 'fields.ListField'
                         or field.get('data_type', None) == 'list' ):
                        # FIXME: need to do an escaped split
                        if hasattr(value, 'split'):
                            value = value.split(LIST_DELIMITER_SQL_ARRAY)

                    if DEBUG_STREAMING: 
                        logger.info('key %r val %r, %r', key, value, type(value))
                    
                    _dict[key] = value
                    
                    if field.get('value_template', None):
                        value_template = field['value_template']
                        if DEBUG_STREAMING: 
                            logger.info('field: %r, value_template: %r', key, value_template)
                        newval = interpolate_value_template(value_template, row)
                        if field['display_type'] == 'image':
                            # hack to speed things up:
                            if ( key == 'structure_image' and
                                    row.has_key('library_well_type') and
                                    row['library_well_type'] == 'empty' ):
                                continue
                            # see if the specified url is available
                            try:
                                view, args, kwargs = resolve(newval)
                                kwargs['request'] = request
                                view(*args, **kwargs)
                                _dict[key] = newval
                            except Exception:
                                logger.info('no image found at %s',newval)
                        else:
                            _dict[key]=newval

            else:
                if DEBUG_STREAMING: logger.info('raw: %r', row)
                _dict = dict((x,y) for x, y in row.items())

            i += 1
            try:
                if i == 1:
                    # NOTE, using "ensure_ascii" = True to force encoding of all 
                    # chars to the ascii charset; otherwise, cStringIO has problems
                    yield json.dumps(_dict, cls=DjangoJSONEncoder,
                        sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
                else:
                    # NOTE, using "ensure_ascii" = True to force encoding of all 
                    # chars to the ascii charset; otherwise, cStringIO has problems
                    # e.g. "tm" becomes \u2122
                    # Upon fp.write  the unicode is converted to the default charset?
                    # also, CStringIO doesn't support UTF-8 mixed with ascii, for instance
                    # NOTE2: control characters are not allowed: should convert \n to "\\n"
                    yield ', ' + json.dumps(_dict, cls=DjangoJSONEncoder,
                        sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
            except Exception, e:
                logger.exception('dict: %r', _dict)
                raise e
        logger.debug('streaming finished')
        
        if not is_for_detail:
            yield ' ] }'

    except Exception, e:
        logger.exception('json streaming')
        raise e                      


class Echo(object):
    """An object that implements just the write method of the file-like
    interface.
    """
    def write(self, value):
        return value

def csv_generator(
        cursor,request,field_hash=None,title_function=None,list_brackets=None):    
    
    pseudo_buffer = Echo()
    quotechar = b'"' # note that csv under python 2.7 doesn't allow multibyte quote char
    csvwriter = csv.writer(
        pseudo_buffer, delimiter=CSV_DELIMITER, quotechar=quotechar, 
        quoting=csv.QUOTE_ALL, lineterminator="\n")
    try:
        i=0
        for row in cursor:
            i += 1
            if i == 1:
                if field_hash:
                    titles = field_hash.keys()
                    logger.info('keys: %r, title_function: %r', titles, title_function)
                    if title_function:
                        titles = [title_function(key) for key in titles]
                    logger.info('titles: %r', titles)
                    yield csvwriter.writerow(titles)
                else:
                    yield csvwriter.writerow(row.keys())
    
            if field_hash:
                values = []
                for col, (key, field) in enumerate(field_hash.iteritems()):
                    value = ''
                    if row.has_key(key):
                        value = row[key]
                    if value and ( 
                        field.get('json_field_type',None) == 'fields.ListField' 
                         or field.get('linked_field_type',None) == 'fields.ListField'
                         or field.get('data_type', None) == 'list' ):
                        # check if it is a string (already may be a list)
                        if hasattr(value, 'split'):
                            value = value.split(LIST_DELIMITER_SQL_ARRAY)
                                            
                    if field.get('value_template', None):
                        value_template = field['value_template']
                        newval = interpolate_value_template(value_template, row)
                        if field['display_type'] == 'image': 
                            if ( key == 'structure_image' and
                                    row.has_key('library_well_type') and
                                    row['library_well_type'] == 'empty' ):
                                # hack to speed things up:
                                continue
                            # see if the specified url is available
                            try:
                                view, args, kwargs = resolve(newval)
                                kwargs['request'] = request
                                response = view(*args, **kwargs)
                                value = newval
                            except Exception, e:
                                logger.info('no image at: %r, %r', newval,e)
                        else:
                            value = newval

                    values.append(value)
                yield csvwriter.writerow([
                    csvutils.csv_convert(val, list_brackets=list_brackets) 
                        for val in values])
            
            else:
                yield csvwriter.writerow([
                    csvutils.csv_convert(val, list_brackets=list_brackets) 
                        for val in row.values()])

    except Exception, e:
        logger.exception('csv streaming error')
        raise e                      


def sdf_generator(cursor, field_hash=None):
    '''
    Yield each field of each row as a string in SDF molfile format
    - chemical/x-mdl-sdfile
    see: http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
    
    @param field_hash optional field definitions:
    key - sql recordset field key
    value - hash defining:
    - data_type
    - value_template
    '''
    try:
    
        i = 0
        for row in cursor:
            i += 1
            if row.has_key(MOLDATAKEY) and row[MOLDATAKEY]:
                yield str(row[MOLDATAKEY])
                yield '\n' 
    
            if field_hash:
                for col, (key, field) in enumerate(field_hash.iteritems()):
                    if key == MOLDATAKEY: 
                        continue
                    yield '> <%s>\n' % key
                    # according to 
                    # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
                    # "only one blank line should terminate a data item"
                    value = None
                    
                    if row.has_key(key):
                        value = row[key]
                    
                    if field.get('value_template', None):
                        value_template = field['value_template']
                        if DEBUG_STREAMING:     
                            logger.info('field: %r:%r, %r', key, value_template)
                        value = interpolate_value_template(value_template, row)
                    
                    if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                         or field.get('linked_field_type',None) == 'fields.ListField'
                         or field.get('data_type', None) == 'list' ):
                        # check if it is a string (already may be a list)
                        if hasattr(value, 'split'):
                            value = value.split(LIST_DELIMITER_SQL_ARRAY)
    
                    if value:
                        # find lists, but not strings (or dicts)
                        # Note: a dict here will be non-standard; probably an error 
                        # report, so just stringify dicts as is.
                        if not hasattr(value, "strip") and isinstance(value, (list,tuple)): 
                            for x in value:
                                # DB should be UTF-8, so this should not be necessary,
                                # however, it appears we have legacy non-utf data in 
                                # some tables (i.e. small_molecule_compound_name 193090
                                yield unicode.encode(x,'utf-8')
                                yield '\n'
                        else:
                            yield str(value)
                            yield '\n'
    
                    yield '\n'
                yield '$$$$\n'
                        
            else:
                logger.info('no field hash defined')
                for k,v in row.items():
                    if k == MOLDATAKEY: 
                        continue
                    yield '> <%s>\n' % k
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
                                yield unicode.encode(x,'utf-8')
                                yield '\n'
                        else:
                            yield str(v)
                            yield '\n'
    
                    yield '\n'
                yield '$$$$\n'
    
        logger.info('wrote %d', i)
    except Exception, e:
        logger.exception('streaming error')
        raise ImmediateHttpResponse('ex during sdf export: %r' % e)     


def cursorGenerator(cursor, field_hash, visible_fields):
    '''
    Generate dicts from cursor rows and visible fields 
    '''
    assert set(visible_fields) <= set(field_hash.keys()), (
        'programming error: field hash: %r must contain all visible fields: %r'
        % (field_hash.keys(), visible_fields))
    
    row_count = 0
    for row in cursor:
        output_row = []
        for key in visible_fields:
            field = field_hash[key]
            value = None
            if row.has_key(key):
                value = row[key]
                     
            if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                 or field.get('linked_field_type',None) == 'fields.ListField'
                 or field.get('data_type', None) == 'list' ):
                value = value.split(LIST_DELIMITER_SQL_ARRAY)
            if field.get('value_template', None):
                value_template = field['value_template']
                if DEBUG_STREAMING: 
                    logger.info('field: %r, value_template: %r', key, value_template)
                value = interpolate_value_template(value_template, row)
            output_row.append(value)
        row_count += 1
        yield dict(zip(visible_fields,output_row))

def CursorGenerator(object):
    '''
    Generate dicts from cursor rows and visible fields 
    '''
    
    def __init__(cursor, schema, visible_fields):
        self.cursor = cursor
        self.schema = schema
        self.visible_fields = visible_fields
        
    def __iter__(self):
        return self
    
    def __next__(self):
        return self.next()
    
    def next(self):
        fields = schema['fields']
        row_count = 0
        for row in cursor:
            if row_count == 0:
                yield self.visible_fields
            
            output_row = []
            for key in self.visible_fields:
                field = fields[key]
                value = None
                if row.has_key(key):
                    value = row[key]
                         
                if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                     or field.get('linked_field_type',None) == 'fields.ListField'
                     or field.get('data_type', None) == 'list' ):
                    value = value.split(LIST_DELIMITER_SQL_ARRAY)
                if field.get('value_template', None):
                    value_template = field['value_template']
                    if DEBUG_STREAMING: 
                        logger.info('field: %r, value_template: %r', key, value_template)
                    value = interpolate_value_template(value_template, row)
                output_row.append(value)
            row_count += 1
            yield output_row


def get_xls_response(
        cursor, output_filename,request,field_hash=None,
        title_function=None, list_brackets=None):
    '''
    Create an xlsx file that will be streamed through the StreamingHttpResponse.
    - if length exceeds MAX_ROWS_PER_XLS_FILE, create multiple files and zip them.
    - TODO: when using xlsx, can simply add extra sheets to the file.
    @param output_filename - for naming temp files

    FIXME: wrap cursor with cursorgenerator; pass in the image columns as arg
    FIXME: rework this using the generic_xlsx_response as a template:
    - this method is used for all xlsx serialization at this time, except 
    for in testing, and in ScreenResultSerializer - 20190419.
    '''

    def write_xls_image(worksheet, filerow, col, val, request):
        '''
        Retrieve and write an image to the worksheet row, col
        '''
        logger.info('write image %r to row: %d col %d', val, filerow, col)
        newval = val
        try:
            view, args, kwargs = resolve(newval)
            kwargs['request'] = request
            response = view(*args, **kwargs)
            image = Image.open(io.BytesIO(response.content))
            height = image.size[1]
            width = image.size[0]
            worksheet.set_row(filerow, height)
            scaling = 0.130 # trial and error width in default excel font
            worksheet.set_column(col,col, width*scaling)
            worksheet.insert_image(
                filerow, col, newval, 
                {'image_data': io.BytesIO(response.content)})
        except Exception, e:
            logger.info('no image at: %r, %r', newval,e)

    # create a temp dir
    # with TemporaryFile() as f:
    temp_dir = os.path.join(
        settings.TEMP_FILE_DIR, str(time.clock()).replace('.', '_'))
    os.mkdir(temp_dir)
    try:
        irow=0
        # Create an new Excel file and add a worksheet.
        filename = '%s.xlsx' % (output_filename)
        temp_file = os.path.join(temp_dir, filename)
        if DEBUG_STREAMING: logger.info('temp file: %r', temp_file)
        workbook = xlsxwriter.Workbook(temp_file)
        worksheet = workbook.add_worksheet()
         
        # FIXME: only need max rows if the file will be too big (structure images)
        # or too long (>65535, for some versions of xls; for that case
        # should implement a mult-sheet solution.
        max_rows_per_file = MAX_ROWS_PER_XLS_FILE
        file_names_to_zip = [temp_file]
        filerow = 0
                             
        for row in cursor:
            if filerow == 0:
                if field_hash:
                    for col, (key, field) in enumerate(field_hash.items()):
                        if title_function:
                            key = title_function(key) 
                        worksheet.write(filerow,col,key)
                else:
                    for col,name in enumerate(row.keys()):
                        worksheet.write(filerow,col,name)
                filerow += 1
                     
            if field_hash:
                for col, (key, field) in enumerate(field_hash.iteritems()):
                    value = None
                    if row.has_key(key):
                        value = row[key]
                             
                    if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                         or field.get('linked_field_type',None) == 'fields.ListField'
                         or field.get('data_type', None) == 'list' ):
                        value = value.split(LIST_DELIMITER_SQL_ARRAY)
                    worksheet.write(filerow, col, 
                        csvutils.csv_convert(value, delimiter=LIST_DELIMITER_XLS, 
                            list_brackets=list_brackets))

                    if field.get('value_template', None):
                        value_template = field['value_template']
                        if DEBUG_STREAMING: logger.info('field: %r,%r,%r', key, value_template)
                        newval = interpolate_value_template(value_template, row)
                        if field['display_type'] == 'image':
                            max_rows_per_file = MAX_IMAGE_ROWS_PER_XLS_FILE
                            # hack to speed things up:
                            if ( key == 'structure_image' and
                                    row.has_key('library_well_type') and
                                    row['library_well_type'] == 'empty' ):
                                continue
                            write_xls_image(worksheet, filerow, col, newval, request)
                        else:
                            worksheet.write(filerow, col, newval)
            else:
                for col,val in enumerate(row.values()):
                    worksheet.write(filerow, col, val)
 
            irow +=1
            filerow +=1
                 
            if irow % max_rows_per_file == 0:
                workbook.close()
                logger.info('wrote file: %r', temp_file)
 
                # Create an new Excel file and add a worksheet.
                filename = '%s_%s.xlsx' % (output_filename, irow)
                temp_file = os.path.join(temp_dir, filename)
                workbook = xlsxwriter.Workbook(temp_file)
                worksheet = workbook.add_worksheet()
                         
                file_names_to_zip.append(temp_file)
                filerow = 0
                     
        workbook.close()
        logger.info('wrote file: %r', temp_file)
 
        content_type = '%s; charset=utf-8' % XLSX_MIMETYPE
        if len(file_names_to_zip) >1:
            # create a temp zip file
            content_type='application/zip; charset=utf-8'
            temp_file = os.path.join('/tmp',str(time.clock()))
            logger.info('temp ZIP file: %r', temp_file)
 
            with ZipFile(temp_file, 'w') as zip_file:
                for _file in file_names_to_zip:
                    zip_file.write(_file, os.path.basename(_file))
            logger.info('wrote file %r', temp_file)
            filename = '%s.zip' % output_filename

        _file = file(temp_file)
        logger.info('download tmp file: %r, %r',temp_file,_file)
        wrapper = FileWrapper(_file)
        response = StreamingHttpResponse(
            wrapper, content_type=content_type) 
        response['Content-Length'] = os.path.getsize(temp_file)
        response['Content-Disposition'] = \
            'attachment; filename=%s' % filename
        return response
    except Exception, e:
        logger.exception('xls streaming error')
        raise e   
    finally:
        try:
            logger.info('rmdir: %r', temp_dir)
            shutil.rmtree(temp_dir)
            if os.path.exists(temp_file):
                logger.info('remove: %r', temp_file)
                os.remove(temp_file)     
        except Exception, e:
            logger.exception('on xlsx & zip file process file: %s' % output_filename)
            raise ImmediateHttpResponse('ex during rmdir: %r, %r' 
                %(output_filename, e))
    

class FileWrapper1:
    """
    Wrapper to convert file-like objects to iterables
        - modified to delete the file after iterating
    """

    def __init__(self, filelike, blksize=8192):
        self.filelike = filelike
        self.blksize = blksize

    def __getitem__(self,key):
        data = self.filelike.read(self.blksize)
        if data:
            return data
        
        # Modify: delete file after iterating
        logger.info('1done...')
        self.filelike.close()
        os.remove(self.filelike.name)

        raise IndexError

    def __iter__(self):
        return self

    def next(self):
        data = self.filelike.read(self.blksize)
        if data:
            return data
        
        # Modify: delete file after iterating
        logger.info('done writing to response...')
        self.filelike.close()
        os.remove(self.filelike.name)
        
        raise StopIteration

def generic_xlsx_response(data):
    '''
    Write out a data dictionary:
    dict keys: named sheets
    values:
    - if dict, convert to rows using dict_to_rows
    - if list, write directly as sheet rows
    - otherwise write as string
    '''
    # using XlsxWriter for constant memory usage
    with  NamedTemporaryFile(delete=False) as temp_file:
        generic_xls_write_workbook(temp_file, data)
        temp_file.seek(0, os.SEEK_END)
        size = temp_file.tell()
        temp_file.seek(0)   
    logger.info('stream to response')
    _file = file(temp_file.name)
    response = StreamingHttpResponse(FileWrapper1(_file)) 
    response['Content-Length'] = size
    response['Content-Type'] = XLSX_MIMETYPE
    return response
               
