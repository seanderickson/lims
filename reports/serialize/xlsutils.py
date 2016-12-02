# TODO: convert to openpyxl
from __future__ import unicode_literals

from collections import OrderedDict
import io
import logging

from tastypie.exceptions import BadRequest
import xlrd
import xlsxwriter

from db.support.data_converter import default_converter
from reports import MAX_IMAGE_ROWS_PER_XLS_FILE
from reports.serialize import csvutils
import reports.serialize


# from PIL import Image
logger = logging.getLogger(__name__)


LIST_DELIMITER_XLS = ';'

def xls_write_workbook(file, data, request=None, image_keys=None, 
    title_function=None, list_brackets=None):
    '''
    ***WARNING*** xlsx files load fully into memory on display - if there are 
    >~ 2000 images, this will cause performance issues on the client.***
    @param sheet_rows iterable of dicts, one per row
    '''

    if not isinstance(data, dict):
        raise BadRequest(
            'unknown data for generic xls serialization: %r' % type(data))

    wb = xlsxwriter.Workbook(file, {'constant_memory': True})
    logger.info('xls_write_workbook for data: %r', data.keys())
    for key, sheet_rows in data.items():
        logger.info('type sheet_rows: %r', type(sheet_rows))
        if isinstance(sheet_rows, (dict, OrderedDict)):
            sheet_name = default_converter(key)
            logger.info('writing sheet %r...', sheet_name)
            sheet = wb.add_worksheet(sheet_name)
            for i, row in enumerate(csvutils.dict_to_rows(sheet_rows)):
                sheet.write_row(i,0,row)
        elif isinstance(sheet_rows, basestring):
            sheet_name = default_converter(key)
            logger.info('writing single string sheet %r...', sheet_name)
            sheet = wb.add_worksheet(sheet_name)
            sheet.write_string(0,0,sheet_rows)
        else:
            sheet_name = default_converter(key)
            logger.info('writing sheets for base name %r...', sheet_name)
            write_rows_to_sheet(wb, sheet_rows, sheet_name, request=request, 
                image_keys=image_keys, title_function=title_function, 
                list_brackets=list_brackets)
    wb.close()
    

def write_rows_to_sheet(wb, sheet_rows, sheet_basename,
    request=None, image_keys=None, title_function=None, list_brackets=None):
    '''
    ***WARNING*** xlsx files load fully into memory on display - if there are 
    >~ 2000 images, this will cause performance issues on the client.***
    @param sheet_rows iterable of dicts, one per row
    '''
    max_rows_per_sheet = 2**20
    
    sheet_name = sheet_basename 
    sheet = wb.add_worksheet(sheet_name)
    filerow = 0
    sheets = 1
    for row,values in enumerate(sheet_rows):
        if filerow >= max_rows_per_sheet:
            sheet_name = '%s_%d' % (sheet_basename,sheets)
            logger.info('rows: %d, max_rows_per_sheet: %d, adding sheet: %s', 
                row,max_rows_per_sheet, sheet_name)
            sheet = wb.add_worksheet(sheet_name)
            filerow = 0
            sheets += 1
        if filerow == 0:
            for i,(key,val) in enumerate(values.items()):
                if title_function:
                    key = title_function(key)
                sheet.write_string(filerow,i,key)
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

    

def generic_xls_write_workbook(file, data):
    '''
    Writes a dict of iterables to a workbook, 
    where key=sheet name, value=iterable that is ready to write (e.g. a list)
    '''
    wb = xlsxwriter.Workbook(file, {'constant_memory': True})
    
    if isinstance(data, dict):
        logger.info('generic_xls_write_workbook for data: %r', data.keys())
        for key, sheet_rows in data.items():
            sheet_name = default_converter(key)
            logger.info('writing sheet %r...', sheet_name)
            sheet = wb.add_worksheet(sheet_name)
            if isinstance(sheet_rows, dict):
                for i, row in enumerate(csvutils.dict_to_rows(sheet_rows)):
                    sheet.write_row(i,0,row)
            elif isinstance(sheet_rows, basestring):
                sheet.write_string(0,0,sheet_rows)
            else:
                generic_write_rows_to_sheet(sheet_rows, sheet)
    else:
        raise BadRequest(
            'unknown data for generic xls serialization: %r' % type(data))
    logger.info('save to file; %r', file.name)
    wb.close()


def generic_write_rows_to_sheet(rows, sheet):
    for row,values in enumerate(rows):
        for i, val in enumerate(values):
            val = csvutils.csv_convert(val, delimiter=LIST_DELIMITER_XLS)
            if val is not None:
                if len(val) > 32767: 
                    logger.error('warn, row too long, %d, key: %r, len: %d', 
                        row,key,len(val) )
                sheet.write_string(row,i,val)

def write_xls_image(worksheet, filerow, col, val, request):
    '''
    Retrieve and write an image to the worksheet row, col
    '''
    logger.debug('write image %r to row: %d col %d', val, filerow, col)
    try:
        image = reports.serialize.resolve_image(request, val)
        fullpath = request.build_absolute_uri(val)
        # view, args, kwargs = resolve(val)
        # kwargs['request'] = request
        # response = view(*args, **kwargs)
        # image = Image.open(io.BytesIO(response.content))
        height = image.size[1]
        width = image.size[0]
        worksheet.set_row(filerow, height)
        scaling = 0.130 # trial and error width in default excel font
        worksheet.set_column(col,col, width*scaling)
        bytes = io.BytesIO()
        image.save(bytes, image.format)
        worksheet.insert_image(filerow, col, fullpath, {'image_data': bytes })
    except Exception, e:
        logger.info('no image at: %r, %r', val,e)

def read_string(cell):
    value = cell.value
    if value is None:
        return None
    elif cell.ctype == xlrd.XL_CELL_NUMBER:
        ival = int(value)
        if value == ival:
            value = ival
        return str(value)
    else:
        value = str(value).strip()
        if not value:
            return None
        return value

def sheet_rows(workbook_sheet):
    def read_row(row):
        for col in range(workbook_sheet.ncols):
            yield read_string(workbook_sheet.cell(row,col))
    for row in range(workbook_sheet.nrows):
        yield read_row(row)

def sheet_rows_dicts(sheet):
    colnames = [xlrd.book.colname(i) for i in range(sheet.ncols)]
    rows = sheet_rows(sheet)
    header = rows.next()
    for row in rows:
        yield dict(zip(header,row))

def sheet_cols(workbook_sheet):
    def read_col(col):
        for row in range(workbook_sheet.nrows):
            yield read_string(workbook_sheet.cell(row,col))
    for col in range(workbook_sheet.ncols):
        yield read_col(col)
    
def workbook_sheets(workbook):
    for sheet_num in range(workbook.nsheets):
        yield workbook.sheet_by_index(sheet_num)

def workbook_as_datastructure(workbook):
    '''
    Create an ordered dict of the sheets in the workbook:
    {
        sheet_name: iterable of sheet rows
    }
    '''
    workbook_datastructure = OrderedDict()
    for sheet_num in range(workbook.nsheets):
        sheet = workbook.sheet_by_index(sheet_num)
        workbook_datastructure[sheet.name] = sheet_rows_dicts(sheet)
    return workbook_datastructure

def workbook_rows(workbook):
    for sheet in workbook_sheets(workbook):
        for row in sheet_rows(sheet):
            yield row