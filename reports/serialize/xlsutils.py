# TODO: convert to openpyxl
import xlrd

LIST_DELIMITER_XLS = ';'

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