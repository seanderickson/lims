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
from tempfile import SpooledTemporaryFile
import time
from wsgiref.util import FileWrapper
from zipfile import ZipFile

from PIL import Image
from django.conf import settings
from django.core.serializers.json import DjangoJSONEncoder
from django.core.urlresolvers import resolve
from django.http.response import StreamingHttpResponse
from tastypie.exceptions import ImmediateHttpResponse
from xlsxwriter.workbook import Workbook

from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    LIST_BRACKETS, MAX_IMAGE_ROWS_PER_XLS_FILE, MAX_ROWS_PER_XLS_FILE, \
    LIST_DELIMITER_XLS, CSV_DELIMITER, HTTP_PARAM_RAW_LISTS
from reports.serializers import csv_convert
import reports.utils.sdf2py


logger = logging.getLogger(__name__)

MOLDATAKEY = reports.utils.sdf2py.MOLDATAKEY

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
                logger.debug(str(('fragment', self.fragment)))
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
            logger.info(str(('val from value template',val, 
                row.has_key(val), row[val])))
        if row.has_key(val):
            return str(row[val])
        else:
            logger.error(str((
                'field needed for value template is not available', 
                val, row)))
            return ''
    return re.sub(r'{([^}]+)}', get_value_from_template, value_template)


def json_generator(cursor,meta,request, is_for_detail=False,field_hash=None):
    
    if DEBUG_STREAMING: logger.info(str(('meta', meta )))
    
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
            if DEBUG_STREAMING: logger.info(str(('row', row, row.keys())))
            if field_hash:
                _dict = OrderedDict()
                for key, field in field_hash.iteritems():
                    if DEBUG_STREAMING: 
                        logger.info(str(('key', key,  str(row), 
                            row.has_key(key))))
                    value = None
                    if row.has_key(key):
                        value = row[key]
                    if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                         or field.get('linked_field_type',None) == 'fields.ListField'
                         or field.get('data_type', None) == 'list' ):
                        # FIXME: need to do an escaped split
                        if DEBUG_STREAMING: logger.info(str(('split', key, value)))
                        if hasattr(value, 'split'):
                            value = value.split(LIST_DELIMITER_SQL_ARRAY)

                    if DEBUG_STREAMING: logger.info(str(('key val', key, value, field)))
                    _dict[key] = value
                    
                    if field.get('value_template', None):
                        value_template = field['value_template']
                        if DEBUG_STREAMING: 
                            logger.info(str(('field', key, value_template)))
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
                if DEBUG_STREAMING: logger.info(str(('raw', row)))
                _dict = dict((x,y) for x, y in row.items())

            i += 1
            try:
                if DEBUG_STREAMING: logger.info(str(('_dict',i, _dict)))
                if i == 1:
                    # NOTE, using "ensure_ascii" = True to force encoding of all 
                    # chars to the ascii charset; otherwise, cStringIO has problems
                    yield json.dumps(_dict, cls=DjangoJSONEncoder,
                        sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
                else:
                    # NOTE, using "ensure_ascii" = True to force encoding of all 
                    # chars to the ascii charset; otherwise, cStringIO has problems
                    # so "tm" becomes \u2122
                    # Upon fp.write  the unicode is converted to the default charset?
                    # also, CStringIO doesn't support UTF-8 mixed with ascii, for instance
                    # NOTE2: control characters are not allowed: should convert \n to "\\n"
                    yield ', ' + json.dumps(_dict, cls=DjangoJSONEncoder,
                        sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
            except Exception, e:
                print 'Exception'
                logger.info(str(('exception')))
                logger.error(str(('ex', _dict, e)))
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

def csv_generator(cursor,request,
        field_hash=None,title_function=None,list_brackets=None):    
    
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
                    logger.info(str(('keys', titles, title_function)))
                    if title_function:
                        titles = [title_function(key) for key in titles]
                    logger.info(str(('titles', titles)))
                    yield csvwriter.writerow(titles)
                else:
                    yield csvwriter.writerow(row.keys())
    
            if field_hash:
                values = []
                for col, (key, field) in enumerate(field_hash.iteritems()):
                    value = ''
                    if row.has_key(key):
                        value = row[key]
                    if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                         or field.get('linked_field_type',None) == 'fields.ListField'
                         or field.get('data_type', None) == 'list' ):
                        # FIXME: must quote special strings?
#                                 value = '[' + ",".join(value.split(LIST_DELIMITER_SQL_ARRAY)) + ']'
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
                                logger.info(str(('no image at', newval,e)))
                        else:
                            value = newval

                    values.append(value)
                yield csvwriter.writerow([
                    csv_convert(val, list_brackets=list_brackets) 
                        for val in values])
            
            else:
                yield csvwriter.writerow([
                    csv_convert(val, list_brackets=list_brackets) 
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
            if DEBUG_STREAMING:
                logger.info(str(('row', i, row)))
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
                            logger.info(str(('field', key, value_template)))
                        value = interpolate_value_template(value_template, row)
                    if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                         or field.get('linked_field_type',None) == 'fields.ListField'
                         or field.get('data_type', None) == 'list' ):
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
                logger.info(str(('no field hash defined')))
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
    
        logger.info(str(('wrote i', i)))    
    except Exception, e:
        logger.exception('streaming error')
        raise ImmediateHttpResponse('ex during sdf export: %r' % e)     


def get_xls_response(cursor, output_filename,request,field_hash=None,
    title_function=None, list_brackets=None):
    '''
    @param output_filename - for naming temp files
    '''
    
    structure_image_dir = os.path.abspath(settings.WELL_STRUCTURE_IMAGE_DIR)
    # create a temp dir
    # FIXME: replace with tempfile.SpooledTemporaryFile
    # with TemporaryFile() as f:
    temp_dir = os.path.join(
        settings.TEMP_FILE_DIR, str(time.clock()).replace('.', '_'))
    os.mkdir(temp_dir)
    try:
        irow=0
        # Create an new Excel file and add a worksheet.
        filename = '%s.xlsx' % (output_filename)
        temp_file = os.path.join(temp_dir, filename)
        if DEBUG_STREAMING: logger.info(str(('temp file', temp_file)))
        workbook = Workbook(temp_file)
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
                        csv_convert(value, delimiter=LIST_DELIMITER_XLS, 
                            list_brackets=list_brackets))

                    if field.get('value_template', None):
                        value_template = field['value_template']
                        if DEBUG_STREAMING: logger.info(str(('field', key, value_template)))
                        newval = interpolate_value_template(value_template, row)
                        if field['display_type'] == 'image':
                            max_rows_per_file = MAX_IMAGE_ROWS_PER_XLS_FILE
                            # hack to speed things up:
                            if ( key == 'structure_image' and
                                    row.has_key('library_well_type') and
                                    row['library_well_type'] == 'empty' ):
                                continue
                            # see if the specified url is available
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
                                logger.info(str(('no image at', newval,e)))
                        else:
                            worksheet.write(filerow, col, newval)
            else:
                for col,val in enumerate(row.values()):
                    worksheet.write(filerow, col, val)
 
            irow +=1
            filerow +=1
                 
            if irow % max_rows_per_file == 0:
                workbook.close()
                logger.info(str(('wrote file', temp_file)))
 
                # Create an new Excel file and add a worksheet.
                filename = '%s_%s.xlsx' % (output_filename, irow)
                temp_file = os.path.join(temp_dir, filename)
                logger.info(str(('temp file', temp_file)))
                workbook = Workbook(temp_file)
                worksheet = workbook.add_worksheet()
                         
                file_names_to_zip.append(temp_file)
                filerow = 0
                     
        workbook.close()
        logger.info(str(('wrote file', temp_file)))
 
        content_type = (
            'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet; '
            'charset=utf-8' )
        if len(file_names_to_zip) >1:
            # create a temp zip file
            content_type='application/zip; charset=utf-8'
            temp_file = os.path.join('/tmp',str(time.clock()))
            logger.info(str(('temp ZIP file', temp_file)))
 
            with ZipFile(temp_file, 'w') as zip_file:
                for _file in file_names_to_zip:
                    zip_file.write(_file, os.path.basename(_file))
            logger.info(str(('wrote file', temp_file)))
            filename = '%s.zip' % output_filename

        _file = file(temp_file)
        logger.info(str(('download tmp file',temp_file,_file)))
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
            logger.info(str(('rmdir', temp_dir)))
            shutil.rmtree(temp_dir)
            if os.path.exists(temp_file):
                logger.info(str(('remove', temp_file)))
                os.remove(temp_file)     
            logger.info(str(('removed', temp_dir)))
        except Exception, e:
            logger.exception('on xlsx & zip file process file: %s' % output_filename)
            raise ImmediateHttpResponse(str(('ex during rmdir', e)))
    
