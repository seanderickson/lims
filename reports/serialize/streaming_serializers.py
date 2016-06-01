# -*- coding: utf-8 -*-
'''
Utilities for streaming sql connection cursors in different serialization formats
'''

from __future__ import unicode_literals

import cStringIO
from collections import OrderedDict
import csv
import json
import logging
import os.path
import re
import shutil
import sys
from tempfile import NamedTemporaryFile
import time
from wsgiref.util import FileWrapper
from zipfile import ZipFile

from django.conf import settings
from django.core.serializers.json import DjangoJSONEncoder
from django.core.urlresolvers import resolve
from django.http.response import StreamingHttpResponse
import six
import xlsxwriter

from db.support.data_converter import default_converter
from reports import LIST_DELIMITER_SQL_ARRAY, \
    MAX_IMAGE_ROWS_PER_XLS_FILE, MAX_ROWS_PER_XLS_FILE, \
    CSV_DELIMITER
from reports.serialize import XLSX_MIMETYPE
from reports.serialize import dict_to_rows
import reports.serialize
import reports.serialize.csvutils as csvutils
import reports.serialize.sdfutils as sdfutils
from reports.serialize.xlsutils import generic_xls_write_workbook, \
    xls_write_workbook, write_xls_image, LIST_DELIMITER_XLS


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
                logger.debug('chunk size %s bytes' % len(bytes.getvalue()))
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


def image_generator(rows, image_keys, request):
    '''
    Check that any image values in the rows can be fetched:
    - replace the raw value given with the absolute URI
    @param rows an iterator that returns a dict for each row
    
    '''
    for row in rows:
        for key,val in row.items():
            if not val:
                continue
            if key in image_keys:
                # hack to speed things up:
                if ( key == 'structure_image' and
                        'library_well_type' in row and
                        row['library_well_type'].lower() == 'empty' ):
                    row[key] = None
                else:
                    try:
                        # Test whether image exists
                        image = reports.serialize.resolve_image(request, val)
                        # If it exists, write the fullpath to the file
                        fullpath = request.build_absolute_uri(val)
                        logger.debug('image exists: %r, abs_uri: %r', val, fullpath)
                        row[key] = fullpath
                    except Exception, e:
                        logger.info('no image at: %r, %r', val,e)
                        row[key] = None
        yield row


def json_generator(data, meta, is_for_detail=False):
    
    if DEBUG_STREAMING: logger.info('meta: %r', meta )
    
    # NOTE, using "ensure_ascii" = True to force encoding of all 
    # chars to be encoded using ascii or escaped unicode; 
    # because some chars in db might be non-UTF8
    # and downstream programs have trouble with mixed encoding (cStringIO)
    if not is_for_detail:
        yield ( '{ "meta": %s, "objects": [' 
            % json.dumps(meta, ensure_ascii=True, encoding="utf-8"))
    try:
        for rownum, row in enumerate(data):
            try:
                if rownum == 0:
                    # NOTE, using "ensure_ascii" = True to force encoding of all 
                    # chars to the ascii charset; otherwise, cStringIO has problems
                    yield json.dumps(row, cls=DjangoJSONEncoder,
                        sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
                else:
                    # NOTE, using "ensure_ascii" = True to force encoding of all 
                    # chars to the ascii charset; otherwise, cStringIO has problems
                    # e.g. "tm" becomes \u2122
                    # CStringIO doesn't support UTF-8 mixed with ascii, for instance
                    # NOTE2: control characters are not allowed: should convert \n to "\\n"
                    yield ', ' + json.dumps(row, cls=DjangoJSONEncoder,
                        sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
            except Exception, e:
                logger.exception('dict: %r', row)
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


def csv_generator(data, title_function=None, list_brackets=None):    
    
    pseudo_buffer = Echo()
    quotechar = b'"' # note that csv under python 2.7 doesn't allow multibyte quote char
    csvwriter = csv.writer(
        pseudo_buffer, delimiter=CSV_DELIMITER, quotechar=quotechar, 
        quoting=csv.QUOTE_ALL, lineterminator="\n")
    try:
        for rownum, row in enumerate(data):
            if rownum == 0:
                titles = row.keys()
                if title_function:
                    titles = [title_function(key) for key in titles]
                yield csvwriter.writerow(titles)

            yield csvwriter.writerow([
                csvutils.csv_convert(val, list_brackets=list_brackets) 
                    for val in row.values()])
        logger.debug('wrote %d rows to csv', rownum)
    except Exception, e:
        logger.exception('CSV streaming error')
        raise e                      


def sdf_generator(data, title_function=None):    
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
        
        for rownum,row in enumerate(data):

            if row.get(MOLDATAKEY, None):
                yield str(row[MOLDATAKEY])
                yield '\n' 

            for i, (key,val) in enumerate(row.items()):
        
                if key == MOLDATAKEY:
                    continue
                title = key
                if title_function:
                    title = title_function(key)
                yield '> <%s>\n' % title

                if val:
                    # Note: find lists, but not strings (or dicts)
                    # Note: a dict here will be non-standard; probably an error 
                    # report, so just stringify dicts as is.
                    if not hasattr(val, "strip") and isinstance(val, (list,tuple)): 
                        for x in val:
                            # DB should be UTF-8, so this should not be necessary,
                            # however, it appears we have legacy non-utf data in 
                            # some tables (i.e. small_molecule_compound_name 193090
                            yield unicode.encode(x,'utf-8')
                            yield '\n'
                    else:
                        yield str(val)
                        yield '\n'
                yield '\n'
            yield '$$$$\n'

        logger.info('wrote %d', rownum)
    except Exception, e:
        logger.exception('SDF streaming error')
        raise


def cursor_generator(cursor, visible_fields, list_fields=[], value_templates=[]):
    '''
    Generate dicts from cursor rows and visible fields 
    '''
    for row in cursor:
        output_row = []
        for key in visible_fields:
            value = None
            if row.has_key(key):
                value = row[key]
                     
            if value and key in list_fields:
                if isinstance(value, six.string_types):
                    value = value.split(LIST_DELIMITER_SQL_ARRAY)
            if key in value_templates:
                value_template = value_templates[key]
                if DEBUG_STREAMING: 
                    logger.info('field: %r, value_template: %r', key, value_template)
                value = interpolate_value_template(value_template, row)
            output_row.append(value)
        yield OrderedDict(zip(visible_fields,output_row))


def get_xls_response(
        data, output_filename,request=None,image_keys=None,
        title_function=None, list_brackets=None):
    '''
    Create an xlsx file that will be streamed through the StreamingHttpResponse.
    - if length exceeds MAX_ROWS_PER_XLS_FILE, create multiple files and zip them.
    - TODO: when using xlsx, can simply add extra sheets to the file.
    @param output_filename - for naming temp files
 
    FIXME: wrap cursor with cursorgenerator; pass in the image columns as arg
    FIXME: rework this using the generic_xlsx_response as a template:
    - this method is used for all xlsx serialization at this time, except 
    for in testing, and in ScreenResultSerializer - 20160419.
    '''
    if not isinstance(data, dict):
        raise BadRequest(
            'unknown data for xls serialization: %r, must be a dict of '
            'sheet_row entries' % type(data))
 
    # create a temp dir
    # with TemporaryFile() as f:
    temp_dir = os.path.join(
        settings.TEMP_FILE_DIR, str(time.clock()).replace('.', '_'))
    os.mkdir(temp_dir)
    try:
        # Create an new Excel file and add a worksheet.
        filename = '%s.xlsx' % (output_filename)
        temp_file = os.path.join(temp_dir, filename)
        file_names_to_zip = [temp_file]
        if DEBUG_STREAMING: logger.info('temp file: %r', temp_file)

        workbook = xlsxwriter.Workbook(temp_file, {'constant_memory': True})
        
        for key, sheet_rows in data.items():
            logger.info('type sheet_rows: %r', type(sheet_rows))
            if isinstance(sheet_rows, (dict, OrderedDict)):
                sheet_name = default_converter(key)
                logger.info('writing sheet %r...', sheet_name)
                sheet = workbook.add_worksheet(sheet_name)
                for i, row in enumerate(dict_to_rows(sheet_rows)):
                    sheet.write_row(i,0,row)
            elif isinstance(sheet_rows, basestring):
                sheet_name = default_converter(key)
                logger.info('writing single string sheet %r...', sheet_name)
                sheet = workbook.add_worksheet(sheet_name)
                sheet.write_string(0,0,sheet_rows)
            else:
                sheet_name = default_converter(key)
                logger.info('writing sheets for base name %r...', sheet_name)

                max_rows_per_sheet = 2**20
                sheet = workbook.add_worksheet(sheet_name)
                filerow = 0
                sheets = 1
                for row,values in enumerate(sheet_rows):
                    if filerow == 0:
                        for i,(key,val) in enumerate(values.items()):
                            title = key
                            if title_function:
                                title = title_function(key)
                            sheet.write_string(filerow,i,title)
                        filerow += 1
                    for i, (key,val) in enumerate(values.items()):
                        val = csvutils.csv_convert(
                            val, delimiter=LIST_DELIMITER_XLS,
                            list_brackets=list_brackets)
                        if val is not None:
                            if len(val) > 32767: 
                                logger.error('warn, row too long, %d, key: %r, len: %d', 
                                    row,key,len(val) )
                            if image_keys and key in image_keys:
                                max_rows_per_sheet = MAX_IMAGE_ROWS_PER_XLS_FILE
                                if not request:
                                    raise Exception(
                                        'must specify the request parameter for image export')
                                # hack to speed things up:
                                if ( key == 'structure_image' and
                                        'library_well_type' in values and
                                        values['library_well_type'].lower() == 'empty' ):
                                    continue
                                write_xls_image(sheet, filerow, i, val, request)
                            else:
                                sheet.write_string(filerow,i,val)
                    filerow += 1
                    if row % 10000 == 0:
                        logger.info('wrote %d rows to temp file', row)
                
                    if filerow > max_rows_per_sheet:
                        workbook.close()
                        logger.info('wrote file: %r', temp_file)
          
                        # Create an new Excel file and add a worksheet.
                        filename = '%s_%s.xlsx' % (output_filename, filerow)
                        temp_file = os.path.join(temp_dir, filename)
                        workbook = xlsxwriter.Workbook(temp_file, {'constant_memory': True})
                        sheet = workbook.add_worksheet(sheet_name)
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
            raise


def get_xls_response_all_images_to_one_file(
        data, output_filename, request=None, image_keys=None,
        title_function=None, list_brackets=None):
    '''
    '''
    # using XlsxWriter for constant memory usage
    max_rows_per_sheet = 2**20

    with  NamedTemporaryFile(delete=False) as temp_file:

        logger.info('save to file; %r...', output_filename)
        xls_write_workbook(temp_file, data, request=request, 
            image_keys=image_keys, title_function=title_function, 
            list_brackets=list_brackets)
        logger.info('saved temp file for; %r', output_filename)
    
        temp_file.seek(0, os.SEEK_END)
        size = temp_file.tell()
        temp_file.seek(0)   

    logger.info('stream to response: file: %r...', output_filename)
    _file = file(temp_file.name)
    response = StreamingHttpResponse(FileWrapper1(_file)) 
    response['Content-Length'] = size
    response['Content-Type'] = XLSX_MIMETYPE
    response['Content-Disposition'] = \
        'attachment; filename=%s.xlsx' % output_filename
    return response


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
        
        # Modified: delete file after iterating
        logger.info('Filewrapper iteration finished, delete temp file %r...',
            self.filelike.name)
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


# def json_generator(cursor,meta,request, is_for_detail=False,field_hash=None):
#     if DEBUG_STREAMING: logger.info('meta: %r', meta )
#     
#     # NOTE, using "ensure_ascii" = True to force encoding of all 
#     # chars to be encoded using ascii or escaped unicode; 
#     # because some chars in db might be non-UTF8
#     # and downstream programs have trouble with mixed encoding (cStringIO)
#     if not is_for_detail:
#         yield ( '{ "meta": %s, "objects": [' 
#             % json.dumps(meta, ensure_ascii=True, encoding="utf-8"))
#     i=0
#     try:
#         for row in cursor:
#             if DEBUG_STREAMING: logger.info('row: %r', row)
#             if field_hash:
#                 _dict = OrderedDict()
#                 for key, field in field_hash.iteritems():
#                     value = None
#                     if row.has_key(key):
#                         value = row[key]
#                     if value and ( field.get('json_field_type',None) == 'fields.ListField' 
#                          or field.get('linked_field_type',None) == 'fields.ListField'
#                          or field.get('data_type', None) == 'list' ):
#                         # FIXME: need to do an escaped split
#                         if hasattr(value, 'split'):
#                             value = value.split(LIST_DELIMITER_SQL_ARRAY)
# 
#                     if DEBUG_STREAMING: 
#                         logger.info('key %r val %r, %r', key, value, type(value))
#                     
#                     _dict[key] = value
#                     
#                     if field.get('value_template', None):
#                         value_template = field['value_template']
#                         if DEBUG_STREAMING: 
#                             logger.info('field: %r, value_template: %r', key, value_template)
#                         newval = interpolate_value_template(value_template, row)
#                         if field['display_type'] == 'image':
#                             # hack to speed things up:
#                             if ( key == 'structure_image' and
#                                     row.has_key('library_well_type') and
#                                     row['library_well_type'] == 'empty' ):
#                                 continue
#                             # see if the specified url is available
#                             try:
#                                 view, args, kwargs = resolve(newval)
#                                 kwargs['request'] = request
#                                 view(*args, **kwargs)
#                                 _dict[key] = newval
#                             except Exception:
#                                 logger.info('no image found at %s',newval)
#                         else:
#                             _dict[key]=newval
# 
#             else:
#                 if DEBUG_STREAMING: logger.info('raw: %r', row)
#                 _dict = dict((x,y) for x, y in row.items())
# 
#             i += 1
#             try:
#                 if i == 1:
#                     # NOTE, using "ensure_ascii" = True to force encoding of all 
#                     # chars to the ascii charset; otherwise, cStringIO has problems
#                     yield json.dumps(_dict, cls=DjangoJSONEncoder,
#                         sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
#                 else:
#                     # NOTE, using "ensure_ascii" = True to force encoding of all 
#                     # chars to the ascii charset; otherwise, cStringIO has problems
#                     # e.g. "tm" becomes \u2122
#                     # Upon fp.write  the unicode is converted to the default charset?
#                     # also, CStringIO doesn't support UTF-8 mixed with ascii, for instance
#                     # NOTE2: control characters are not allowed: should convert \n to "\\n"
#                     yield ', ' + json.dumps(_dict, cls=DjangoJSONEncoder,
#                         sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
#             except Exception, e:
#                 logger.exception('dict: %r', _dict)
#                 raise e
#         logger.debug('streaming finished')
#         
#         if not is_for_detail:
#             yield ' ] }'
# 
#     except Exception, e:
#         logger.exception('json streaming')
#         raise e                      

# def csv_generator(
#         cursor,request,field_hash=None,title_function=None,list_brackets=None):    
#     
#     pseudo_buffer = Echo()
#     quotechar = b'"' # note that csv under python 2.7 doesn't allow multibyte quote char
#     csvwriter = csv.writer(
#         pseudo_buffer, delimiter=CSV_DELIMITER, quotechar=quotechar, 
#         quoting=csv.QUOTE_ALL, lineterminator="\n")
#     try:
#         i=0
#         for row in cursor:
#             i += 1
#             if i == 1:
#                 if field_hash:
#                     titles = field_hash.keys()
#                     logger.info('keys: %r, title_function: %r', titles, title_function)
#                     if title_function:
#                         titles = [title_function(key) for key in titles]
#                     logger.info('titles: %r', titles)
#                     yield csvwriter.writerow(titles)
#                 else:
#                     yield csvwriter.writerow(row.keys())
#     
#             if field_hash:
#                 values = []
#                 for col, (key, field) in enumerate(field_hash.iteritems()):
#                     value = ''
#                     if row.has_key(key):
#                         value = row[key]
#                     if value and ( 
#                         field.get('json_field_type',None) == 'fields.ListField' 
#                          or field.get('linked_field_type',None) == 'fields.ListField'
#                          or field.get('data_type', None) == 'list' ):
#                         # check if it is a string (already may be a list)
#                         if hasattr(value, 'split'):
#                             value = value.split(LIST_DELIMITER_SQL_ARRAY)
#                                             
#                     if field.get('value_template', None):
#                         value_template = field['value_template']
#                         newval = interpolate_value_template(value_template, row)
#                         if field['display_type'] == 'image': 
#                             if ( key == 'structure_image' and
#                                     row.has_key('library_well_type') and
#                                     row['library_well_type'] == 'empty' ):
#                                 # hack to speed things up:
#                                 continue
#                             # see if the specified url is available
#                             try:
#                                 view, args, kwargs = resolve(newval)
#                                 kwargs['request'] = request
#                                 response = view(*args, **kwargs)
#                                 value = newval
#                             except Exception, e:
#                                 logger.info('no image at: %r, %r', newval,e)
#                         else:
#                             value = newval
# 
#                     values.append(value)
#                 yield csvwriter.writerow([
#                     csvutils.csv_convert(val, list_brackets=list_brackets) 
#                         for val in values])
#             
#             else:
#                 yield csvwriter.writerow([
#                     csvutils.csv_convert(val, list_brackets=list_brackets) 
#                         for val in row.values()])
# 
#     except Exception, e:
#         logger.exception('csv streaming error')
#         raise e                      


# def sdf_generator(
#         cursor,request,field_hash=None,title_function=None,image_keys=None):    
#     '''
#     Yield each field of each row as a string in SDF molfile format
#     - chemical/x-mdl-sdfile
#     see: http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
#     
#     @param field_hash optional field definitions:
#     key - sql recordset field key
#     value - hash defining:
#     - data_type
#     - value_template
#     '''
#     try:
#     
#         i = 0
#         for row in cursor:
#             i += 1
#             if row.has_key(MOLDATAKEY) and row[MOLDATAKEY]:
#                 yield str(row[MOLDATAKEY])
#                 yield '\n' 
#     
#             if field_hash:
#                 for col, (key, field) in enumerate(field_hash.iteritems()):
#                     if key == MOLDATAKEY: 
#                         continue
#                     yield '> <%s>\n' % key
#                     # according to 
#                     # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
#                     # "only one blank line should terminate a data item"
#                     value = None
#                     
#                     if row.has_key(key):
#                         value = row[key]
#                     
#                     if field.get('value_template', None):
#                         value_template = field['value_template']
#                         if DEBUG_STREAMING:     
#                             logger.info('field: %r:%r, %r', key, value_template)
#                         value = interpolate_value_template(value_template, row)
#                     
#                     if value and ( field.get('json_field_type',None) == 'fields.ListField' 
#                          or field.get('linked_field_type',None) == 'fields.ListField'
#                          or field.get('data_type', None) == 'list' ):
#                         # check if it is a string (already may be a list)
#                         if hasattr(value, 'split'):
#                             value = value.split(LIST_DELIMITER_SQL_ARRAY)
#     
#                     if value:
#                         # find lists, but not strings (or dicts)
#                         # Note: a dict here will be non-standard; probably an error 
#                         # report, so just stringify dicts as is.
#                         if not hasattr(value, "strip") and isinstance(value, (list,tuple)): 
#                             for x in value:
#                                 # DB should be UTF-8, so this should not be necessary,
#                                 # however, it appears we have legacy non-utf data in 
#                                 # some tables (i.e. small_molecule_compound_name 193090
#                                 yield unicode.encode(x,'utf-8')
#                                 yield '\n'
#                         else:
#                             yield str(value)
#                             yield '\n'
#     
#                     yield '\n'
#                 yield '$$$$\n'
#                         
#             else:
#                 logger.info('no field hash defined')
#                 for k,v in row.items():
#                     if k == MOLDATAKEY: 
#                         continue
#                     yield '> <%s>\n' % k
#                     # according to 
#                     # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
#                     # "only one blank line should terminate a data item"
#                     if v:
#                         # find lists, but not strings (or dicts)
#                         # Note: a dict here will be non-standard; probably an error 
#                         # report, so just stringify dicts as is.
#                         if not hasattr(v, "strip") and isinstance(v, (list,tuple)): 
#                             for x in v:
#                                 # DB should be UTF-8, so this should not be necessary,
#                                 # however, it appears we have legacy non-utf data in 
#                                 # some tables (i.e. small_molecule_compound_name 193090
#                                 yield unicode.encode(x,'utf-8')
#                                 yield '\n'
#                         else:
#                             yield str(v)
#                             yield '\n'
#     
#                     yield '\n'
#                 yield '$$$$\n'
#     
#         logger.info('wrote %d', i)
#     except Exception, e:
#         logger.exception('streaming error')
#         raise ImmediateHttpResponse('ex during sdf export: %r' % e)     

               
# def cursor_generator(cursor, field_hash, visible_fields):
#     '''
#     Generate dicts from cursor rows and visible fields 
#     '''
#     assert set(visible_fields) <= set(field_hash.keys()), (
#         'programming error: field hash: %r must contain all visible fields: %r'
#         % (field_hash.keys(), visible_fields))
#     
#     for row in cursor:
#         output_row = []
#         for key in visible_fields:
#             field = field_hash[key]
#             value = None
#             if row.has_key(key):
#                 value = row[key]
#                      
#             if value and ( field.get('json_field_type',None) == 'fields.ListField' 
#                  or field.get('linked_field_type',None) == 'fields.ListField'
#                  or field.get('data_type', None) == 'list' ):
#                 value = value.split(LIST_DELIMITER_SQL_ARRAY)
#             if field.get('value_template', None):
#                 value_template = field['value_template']
#                 if DEBUG_STREAMING: 
#                     logger.info('field: %r, value_template: %r', key, value_template)
#                 value = interpolate_value_template(value_template, row)
#             output_row.append(value)
#         yield OrderedDict(zip(visible_fields,output_row))
               
               
# def get_xls_response(
#         cursor, output_filename,request,field_hash=None,
#         title_function=None, list_brackets=None):
#     '''
#     Create an xlsx file that will be streamed through the StreamingHttpResponse.
#     - if length exceeds MAX_ROWS_PER_XLS_FILE, create multiple files and zip them.
#     - TODO: when using xlsx, can simply add extra sheets to the file.
#     @param output_filename - for naming temp files
# 
#     FIXME: wrap cursor with cursorgenerator; pass in the image columns as arg
#     FIXME: rework this using the generic_xlsx_response as a template:
#     - this method is used for all xlsx serialization at this time, except 
#     for in testing, and in ScreenResultSerializer - 20190419.
#     '''
# 
#     # create a temp dir
#     # with TemporaryFile() as f:
#     temp_dir = os.path.join(
#         settings.TEMP_FILE_DIR, str(time.clock()).replace('.', '_'))
#     os.mkdir(temp_dir)
#     try:
#         irow=0
#         # Create an new Excel file and add a worksheet.
#         filename = '%s.xlsx' % (output_filename)
#         temp_file = os.path.join(temp_dir, filename)
#         if DEBUG_STREAMING: logger.info('temp file: %r', temp_file)
#         workbook = xlsxwriter.Workbook(temp_file)
#         worksheet = workbook.add_worksheet()
#          
#         # FIXME: only need max rows if the file will be too big (structure images)
#         # or too long (>65535, for some versions of xls; for that case
#         # should implement a mult-sheet solution.
#         max_rows_per_file = MAX_ROWS_PER_XLS_FILE
#         file_names_to_zip = [temp_file]
#         filerow = 0
#                              
#         for row in cursor:
#             if filerow == 0:
#                 if field_hash:
#                     for col, (key, field) in enumerate(field_hash.items()):
#                         if title_function:
#                             key = title_function(key) 
#                         worksheet.write(filerow,col,key)
#                 else:
#                     for col,name in enumerate(row.keys()):
#                         worksheet.write(filerow,col,name)
#                 filerow += 1
#                      
#             if field_hash:
#                 for col, (key, field) in enumerate(field_hash.iteritems()):
#                     value = None
#                     if row.has_key(key):
#                         value = row[key]
#                              
#                     if value and ( field.get('json_field_type',None) == 'fields.ListField' 
#                          or field.get('linked_field_type',None) == 'fields.ListField'
#                          or field.get('data_type', None) == 'list' ):
#                         value = value.split(LIST_DELIMITER_SQL_ARRAY)
#                     worksheet.write(filerow, col, 
#                         csvutils.csv_convert(value, delimiter=LIST_DELIMITER_XLS, 
#                             list_brackets=list_brackets))
# 
#                     if field.get('value_template', None):
#                         value_template = field['value_template']
#                         if DEBUG_STREAMING: logger.info('field: %r,%r,%r', key, value_template)
#                         newval = interpolate_value_template(value_template, row)
#                         if field['display_type'] == 'image':
#                             max_rows_per_file = MAX_IMAGE_ROWS_PER_XLS_FILE
#                             # hack to speed things up:
#                             if ( key == 'structure_image' and
#                                     row.has_key('library_well_type') and
#                                     row['library_well_type'] == 'empty' ):
#                                 continue
#                             write_xls_image(worksheet, filerow, col, newval, request)
#                         else:
#                             worksheet.write(filerow, col, newval)
#             else:
#                 for col,val in enumerate(row.values()):
#                     worksheet.write(filerow, col, val)
#  
#             irow +=1
#             filerow +=1
#                  
#             if irow % max_rows_per_file == 0:
#                 workbook.close()
#                 logger.info('wrote file: %r', temp_file)
#  
#                 # Create an new Excel file and add a worksheet.
#                 filename = '%s_%s.xlsx' % (output_filename, irow)
#                 temp_file = os.path.join(temp_dir, filename)
#                 workbook = xlsxwriter.Workbook(temp_file)
#                 worksheet = workbook.add_worksheet()
#                          
#                 file_names_to_zip.append(temp_file)
#                 filerow = 0
#                      
#         workbook.close()
#         logger.info('wrote file: %r', temp_file)
#  
#         content_type = '%s; charset=utf-8' % XLSX_MIMETYPE
#         if len(file_names_to_zip) >1:
#             # create a temp zip file
#             content_type='application/zip; charset=utf-8'
#             temp_file = os.path.join('/tmp',str(time.clock()))
#             logger.info('temp ZIP file: %r', temp_file)
#  
#             with ZipFile(temp_file, 'w') as zip_file:
#                 for _file in file_names_to_zip:
#                     zip_file.write(_file, os.path.basename(_file))
#             logger.info('wrote file %r', temp_file)
#             filename = '%s.zip' % output_filename
# 
#         _file = file(temp_file)
#         logger.info('download tmp file: %r, %r',temp_file,_file)
#         wrapper = FileWrapper(_file)
#         response = StreamingHttpResponse(
#             wrapper, content_type=content_type) 
#         response['Content-Length'] = os.path.getsize(temp_file)
#         response['Content-Disposition'] = \
#             'attachment; filename=%s' % filename
#         return response
#     except Exception, e:
#         logger.exception('xls streaming error')
#         raise e   
#     finally:
#         try:
#             logger.info('rmdir: %r', temp_dir)
#             shutil.rmtree(temp_dir)
#             if os.path.exists(temp_file):
#                 logger.info('remove: %r', temp_file)
#                 os.remove(temp_file)     
#         except Exception, e:
#             logger.exception('on xlsx & zip file process file: %s' % output_filename)
#             raise ImmediateHttpResponse('ex during rmdir: %r, %r' 
#                 %(output_filename, e))
