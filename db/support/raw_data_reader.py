from __future__ import unicode_literals

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
import logging
import os
import re

import xlrd
import reports.serialize.xlsutils as xlsutils
import xlsxwriter
from db.support import lims_utils


logger = logging.getLogger(__name__)

ALLOWED_COLS = [12,24,48]
ALLOWED_ROWS = [8,16,32]

ERROR_PLATE_SIZE = 'Plate size error'
ERROR_COL_SIZE = 'Column count error'
ERROR_ROW_SIZE = 'Row count error'
ERROR_ROW_PATTERN_MATCH = 'Value error'
      
def parse_to_numeric(matrices):
    errors = []
    for matrix_index, matrix2d in enumerate(matrices):
        for j,row in enumerate(matrix2d):
            for i,x in enumerate(row):
                try:
                    row[i] = Decimal(str(x))
                except Exception, e:
                    logger.warn('matrix parse error: %r', e)
                    errors.append('matrix: %d, row: %d, col: %d, non numeric value: %r'
                        % (matrix_index, j, i, x)) 
    return errors

def is_empty_plate_matrix(matrix2d):
    is_empty = True
    for row in matrix2d:
        row_empty = True
        for x in row:
            if x:
                x = x.strip()
            if x and x != '-':
                row_empty = False
                break
        if row_empty is False:
            is_empty = False
            break
        
    return is_empty

def read_rows(rowgenerator):
    plate_matrices = []
    errors = defaultdict(list)
    
    plate_matrix = None
    in_matrix = False
#     expected_cols = 0
    cols_detected = None
    rows_detected = None
    
    for i,cellgenerator in enumerate(rowgenerator):
        row = [x for x in cellgenerator]
        logger.debug('read row: %d: %r', i, row)

        # Matrix Header row: 
        # - empty cell followed by all number cells; 
        # - only consider if length is in ALLOWED_COLS 
        cell0 = None
        if row:
            cell0 = row[0]
            if cell0:
                cell0 = cell0.strip()
        logger.debug('cell0: %r', cell0)
        if not cell0:
            if in_matrix is False:
                header_row = []
                for possible_header_cell in row[1:]:
                    if ( possible_header_cell 
                        and re.match(r'\d{1,2}',possible_header_cell)):
                        header_row.append(x)
                    else:
                        break
                if len(header_row)>0:
                    _cols = len(header_row)
                    
                    if cols_detected is None:
                        if _cols not in ALLOWED_COLS:
                            msg = (
                                'row: %d, columns found: %d, '
                                'must be one of: %r, %r' 
                                % (i, _cols, ALLOWED_COLS, row))
                            logger.error(ERROR_PLATE_SIZE + ': ' + msg)
                            errors[ERROR_PLATE_SIZE].append(msg)
                            continue
                        else:
                            cols_detected = _cols
                            rows_detected = int(_cols*2/3)
                            logger.info('row: %d, cols detected: %r, rows: %r', 
                                i, cols_detected, rows_detected)
                    elif _cols != cols_detected:
                        msg = (
                            'row: %d, header columns found: %d, '
                            'columns detected: %d: %r'
                            % (i, _cols, cols_detected, row))
                        logger.error(ERROR_PLATE_SIZE + ': ' + msg)
                        errors[ERROR_PLATE_SIZE].append(msg)
                        continue
                    in_matrix = True
                    plate_matrix = []
                    plate_matrices.append(plate_matrix)
            else: 
                # in matrix; blank cell found; signals end of matrix
                if is_empty_plate_matrix(plate_matrix) is False:
                    _rows = len(plate_matrix)
                    if _rows != rows_detected:
                        msg = (
                            'row: %d, plate_matrix: %d, '
                            'rows found: %d, rows detected %d' 
                            % (i, len(plate_matrices), _rows, rows_detected))
                        logger.error(ERROR_ROW_SIZE + ': ' + msg)
                        errors[ERROR_ROW_SIZE].append(msg)
                    if _rows < rows_detected:
                        msg = (
                            'row: %d, plate_matrix: %d, must be '
                            'a row letter followed by numeric values only: %r'
                            % (i, len(plate_matrices),line))
                        errors[ERROR_ROW_PATTERN_MATCH].append(msg)
                    
                else:
                    logger.info('empty matrix found: %r', plate_matrix)
                    del plate_matrices[-1]
                in_matrix = False
        # Matrix grid rows:
        # - wellname in the first cell
        # - only 2 digit numbers in all other cells
        if cell0 and in_matrix is True:
            if re.match(r'^[A-Z]{1,2}$', cell0, re.IGNORECASE):
                # Once the header row is found, make this a provisional
                # matrix; not validated until the correct number of grid rows 
                # have been read in, with no empty rows, or only valid grid 
                # rows found after the header.
                _cols = len(row[1:])
                if _cols < cols_detected:
                    msg = (
                        'line: %d, matrix: %d, row: %d, not enough cols: %r'
                        % (i, len(plate_matrices), len(plate_matrix)+1,row))
                    logger.error(ERROR_COL_SIZE + ': ' + msg)
                    errors[ERROR_COL_SIZE].append(msg)
                else:
                    plate_matrix.append(row[1:cols_detected+1])
            else:
                # in matrix; found a cell0 that does not match row letter pattern
                msg = ('line: %d, matrix: %d, row letter not found for row: %d, %r'
                    % (i, len(plate_matrices), len(plate_matrix)+1, row))
                logger.error(ERROR_COL_SIZE + ': ' + msg)
                errors[ERROR_COL_SIZE].append(msg)
                in_matrix = False
                
    logger.info('plate matrices: %r', len(plate_matrices))        
    return (plate_matrices, errors)

def read_text(input_file):
    
    DEBUG = False or logger.isEnabledFor(logging.DEBUG)
    
    header_pattern = re.compile(r'^\s?(((\d{1,2})\s?)+)$')
    row_pattern = re.compile(
        r'^\s?([A-Z]{1,2}(\s+([\d\.E+-]+)\s?)+)$',flags=re.IGNORECASE)
    plate_matrices = []
    errors = defaultdict(list)
    plate_matrix = None
    in_matrix = False
    cols_detected = None
    rows_detected = None
    for i,line in enumerate(input_file):
        if DEBUG:
            logger.info('read line: %d: %r', i, line)
        if i < 10:
            logger.debug('read line: %d: %r', i, line)
        line = line.strip()
        if not line:
            continue
        header_match = header_pattern.match(line)
        if header_match:
            logger.debug('recognized header %d: %r', len(plate_matrices), line)
            _cols = len(re.split(r'\s+', line))
            logger.debug('_cols: %r', _cols)
            if cols_detected is None:
                if _cols not in ALLOWED_COLS:
                    # matrices are 2/3 ratio, so 12,24, or 48 in width
                    msg = (
                        'row: %d, columns found: %d, '
                        'must be one of: %r, %r' 
                        % (i, _cols, ALLOWED_COLS, line))
                    logger.error(ERROR_PLATE_SIZE + ': ' + msg)
                    errors[ERROR_PLATE_SIZE].append(msg)
                else:
                    cols_detected = _cols
                    rows_detected = int(_cols*2/3)
                    logger.info('row: %d, cols detected: %r, rows: %r', 
                        i, cols_detected, rows_detected)
            else:
                if _cols != cols_detected:
                    msg = (
                        'row: %d, header columns found: %d, '
                        'columns detected: %d: %r'
                        % (i, _cols, cols_detected, line))
                    logger.error(ERROR_PLATE_SIZE + ': ' + msg)
                    errors[ERROR_PLATE_SIZE].append(msg)
                    
            in_matrix = True
            plate_matrix = []
            plate_matrices.append(plate_matrix)
            continue
        if in_matrix:
            if row_pattern.match(line):
                row = re.split(r'\s+',line)
                _cols = len(row[1:])
                if _cols < cols_detected:
                    msg = (
                        'row: %d, matrix: %d, line: %d, '
                        'columns found: %d, columns detected: %d, %r'
                        % (i, len(plate_matrices), len(plate_matrix)+1,
                            len(row[1:]), cols_detected,row))
                    logger.error(ERROR_COL_SIZE + ': ' + msg)
                    errors[ERROR_COL_SIZE].append(msg)
                else:
                    plate_matrix.append(row[1:cols_detected+1])
            else:
                if is_empty_plate_matrix(plate_matrix) is False:
                    _rows = len(plate_matrix)
                    
                    if _rows != rows_detected:
                        msg = (
                            'row: %d, plate_matrix: %d, '
                            'rows found: %d, rows detected %d' 
                            % (i, len(plate_matrices), _rows, rows_detected))
                        logger.error(ERROR_ROW_SIZE + ': ' + msg)
                        errors[ERROR_ROW_SIZE].append(msg)
                    if _rows < rows_detected:
                        msg = (
                            'row: %d, plate_matrix: %d, must be '
                            'a row letter followed by numeric values only: %r'
                            % (i, len(plate_matrices),line))
                        errors[ERROR_ROW_PATTERN_MATCH].append(msg)
                else:
                    del plate_matrices[-1]
                    
                in_matrix = False
    if not plate_matrices:
        errors['parse error'] = 'no data read'
    logger.info('plate matrices: %d', len(plate_matrices))
    
    return (plate_matrices, errors)
    
def read_xlsx(input_file):

    wb = xlrd.open_workbook(file_contents=input_file.read())
    sheets = xlsutils.workbook_sheets(wb)
    plate_matrices = []
    errors = {}
    for sheet in sheets:
        logger.info('read sheet: %r', sheet)
        (sheet_plate_matrices,sheet_errors) = read_rows(xlsutils.sheet_rows(sheet))
        logger.info('for sheet: %r, %d matrices read', 
            sheet.name, len(sheet_plate_matrices))
        if sheet_plate_matrices:
            plate_matrices += sheet_plate_matrices
            logger.info('plate_matrices: %r', len(plate_matrices))
        if sheet_errors:
            logger.warn('sheet errors: %r: %r', sheet.name, sheet_errors)
            errors[sheet.name] = sheet_errors
            
    if not plate_matrices:
        errors['No matrices were found']
    return (plate_matrices, errors)

def read_csv(input_file):

    delimiter = str(',')
    try:
        dialect = csv.Sniffer().sniff(
            input_file.read(4096),delimiters=[',', '\t', ';', ' ', ':','^'])
        logger.info('delimiter determined: %r', dialect.delimiter)
        delimiter = dialect.delimiter
    except:
        logger.warn(
            'could not determine CSV delimiter, choosing default: %r', delimiter)
    input_file.seek(0)
    reader = csv.reader(input_file, delimiter=delimiter)
    (plate_matrices, errors) = read_rows(reader)
    
    return (plate_matrices, errors)

def read(input_file, filename):

    extension = os.path.splitext(filename)[-1]
    logger.info('file: %r, extension: %r',filename, extension)
    plate_matrices = []
    errors = {}
    if extension in ['.xls', '.xlsx']:
        logger.info('opening XLS file: %r', filename)
        (plate_matrices,errors) = read_xlsx(input_file)
    elif extension in ['.csv']:
        logger.info('opening CSV file: %r', filename)
        (plate_matrices,errors) = read_csv(input_file)
    else:
        logger.info('Assume file: %r is text', filename)
        (plate_matrices,errors) = read_text(input_file)
                        
    return (plate_matrices,errors)

def write_xlsx(workbook, matrices):
    '''
    Write the plate matrices directly to a spreadsheet:
    - transpose row/col into a single column by wellname
    '''
    DEBUG = False or logger.isEnabledFor(logging.DEBUG)

    logger.info('write workbook...')
    
    row_size = len(matrices[0])
    col_size = len(matrices[0][0])

    sheet = workbook.add_worksheet('Data')
    
    sheet.write_string(0,0,'Matrix')
    sheet.write_string(0,1,'Well')
    sheet.write_string(0,2,'Value')
    
    output_row = 0
    for matrix_number, matrix in enumerate(matrices):
        logger.info('write matrix: %r', matrix_number+1)
        # NOTE: to support xlsxwriter 'constant_memory': True - optimized write, 
        # rows must be written sequentially        
        for colnum in range(0, col_size):
            for rownum in range(0,row_size):
                output_row += 1
                sheet.write_number(output_row,0,matrix_number+1)
                sheet.write_string(
                    output_row,1,lims_utils.get_well_name(rownum, colnum))
                val = matrix[rownum][colnum]
                if val:
                    sheet.write_string(output_row,2,val)

        
parser = argparse.ArgumentParser(
    description='Parse raw data files into Python matrices')

parser.add_argument(
    '-f', '--file', required=True,
    help='''Raw data file;
            examples can be found in db/static/test_data/rawdata''')

parser.add_argument(
    '-of', '--outputfile', 
    help='''Output file (xlsx) to write''')

parser.add_argument(
    '-v', '--verbose', dest='verbose', action='count',
    help="Increase verbosity (specify multiple times for more)")    

if __name__ == "__main__":
    args = parser.parse_args()
    log_level = logging.WARNING # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
        DEBUG=True
    logging.basicConfig(
        level=log_level, 
        format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')        

    with open(args.file) as input_file:
        (plate_matrices,errors) = read(input_file, args.file)
        
        if errors:
            raise Exception('Parse errors: %r' % errors)
        
        if args.outputfile:
            logger.info('write to output file: %r', args.outputfile)
            workbook = xlsxwriter.Workbook(filename=args.outputfile)
            write_xlsx(workbook, plate_matrices)
            workbook.close()
            logger.info('done: output file: %r', args.outputfile)
        else:
            for i,matrix in enumerate(plate_matrices):
                print '\nmatrix #: %d\n' % i
                for j,row in enumerate(matrix):
                    print lims_utils.row_to_letter(j),row

