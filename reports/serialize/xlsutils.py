# TODO: convert to openpyxl
import xlrd
import xlsxwriter
from reports.serialize import dict_to_rows, csvutils
from tastypie.exceptions import BadRequest
import logging
from db.support.data_converter import default_converter

logger = logging.getLogger(__name__)


LIST_DELIMITER_XLS = ';'

def generic_xls_write_workbook(file, data):
    '''Writes a dict of iterables to a workbook, 
    where each iterable is ready to write
    '''
    wb = xlsxwriter.Workbook(file, {'constant_memory': True})
    
    if isinstance(data, dict):
        logger.info('generic_xls_write_workbook for data: %r', data.keys())
        for key, sheet_rows in data.items():
            sheet_name = default_converter(key)
            logger.info('writing sheet %r...', sheet_name)
            sheet = wb.add_worksheet(sheet_name)
            if isinstance(sheet_rows, dict):
                for i, row in enumerate(dict_to_rows(sheet_rows)):
                    sheet.write_row(i,0,row)
            elif isinstance(sheet_rows, basestring):
                sheet.write_string(0,0,sheet_rows)
            else:
                write_rows_to_sheet(sheet_rows, sheet)
    else:
        raise BadRequest(
            'unknown data for generic xls serialization: %r' % type(data))
    logger.info('save to file; %r', file.name)
    wb.close()

def write_rows_to_sheet(rows, sheet):
    for row,values in enumerate(rows):
        for i, val in enumerate(values):
            val = csvutils.csv_convert(val, delimiter=LIST_DELIMITER_XLS)
            if val:
                if len(val) > 32767: 
                    logger.error('warn, row too long, %d, key: %r, len: %d', 
                        row,key,len(val) )
                sheet.write_string(row,i,val)


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

def sheet_rows_by_colnames(sheets):
    for sheet in sheets:
#         colnames = [chr(ord('A')+i) for i in range(sheet.ncols)]
        colnames = [xlrd.book.colname(i) for i in range(sheet.ncols)]
        rows = sheet_rows(sheet)
        for row in rows:
            yield dict(zip(colnames,row))

def sheet_cols(workbook_sheet):
    def read_col(col):
        for row in range(workbook_sheet.nrows):
            yield read_string(workbook_sheet.cell(row,col))
    for col in range(workbook_sheet.ncols):
        yield read_col(col)
    
def workbook_sheets(workbook):
    for sheet_num in range(workbook.nsheets):
        yield workbook.sheet_by_index(sheet_num)

def workbook_rows(workbook):
    for sheet in workbook_sheets(workbook):
        for row in sheet_rows(sheet):
            yield row