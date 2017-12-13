from __future__ import unicode_literals

import argparse
from collections import defaultdict
import csv
import logging
import os
import re

import xlrd

import reports.serialize.xlsutils as xlsutils


logger = logging.getLogger(__name__)

ALLOWED_COLS = [12,24,48]
ALLOWED_ROWS = [8,16,32]

ERROR_PLATE_SIZE = 'Plate size'
ERROR_COL_SIZE = 'Columns detected'
ERROR_ROW_SIZE = 'Rows detected'
            
def is_empty_plate_matrix(matrix2d):
    is_empty = True
    for row in matrix2d:
        row_empty = True
        for x in row:
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
    expected_cols = 0
    
    for i,cellgenerator in enumerate(rowgenerator):
        row = [x for x in cellgenerator]
        logger.debug('read row: %d: %r', i, row)

        # Matrix Header row: 
        # - empty cell followed by all number cells; 
        # - only consider if length is in ALLOWED_COLS 
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
                    expected_cols = len(header_row)
                    logger.debug('expected_cols: %r', expected_cols)
                    if expected_cols not in ALLOWED_COLS:
                        # matrices are 2/3 ratio, so 12,24, or 36 in width
                        msg = (
                            'Not enough cols: ineligible: row: %d, cols recognized: %d, '
                            'must be one of: %r, %r'
                            % (i, expected_cols, ALLOWED_COLS, row))
                        logger.info(msg)
                        continue
                        # raise Exception(msg)
                    in_matrix = True
                    plate_matrix = []
                    plate_matrices.append(plate_matrix)
            elif in_matrix:
                if is_empty_plate_matrix(plate_matrix) is False:
                    if len(plate_matrix) not in ALLOWED_ROWS:
                        msg = (
                            'row: %d, plate_matrix: %d read fail, '
                            'rows found: %r, must be one of %r' 
                            % (i, len(plate_matrices), len(plate_matrix), 
                                ALLOWED_ROWS))
                        logger.error(ERROR_ROW_SIZE + ': ' + msg)
                        errors[ERROR_ROW_SIZE].append(msg)
                        for i,row in enumerate(plate_matrix):
                            logger.error('%r: %r', i, row)
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
                # TODO: could verify that all values are numerical
                if expected_cols > 0:
                    if len(row[1:]) < expected_cols:
                        msg = (
                            'row: %d, matrix: %d, line: %d, not enough cols: %r'
                            % (i, len(plate_matrices), len(plate_matrix)+1,line))
                        logger.error(ERROR_COL_SIZE + ': ' + msg)
                        errors[ERROR_COL_SIZE].append(msg)
                    else:
                        plate_matrix.append(row[1:expected_cols+1])
                else:
                    plate_matrix.append(row[1:])
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
    expected_cols = 0
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
            expected_cols = len(re.split(r'\s+', line))
            logger.debug('expected_cols: %r', expected_cols)
            if expected_cols not in ALLOWED_COLS:
                # matrices are 2/3 ratio, so 12,24, or 48 in width
                msg = (
                    'row: %d, cols recognized: %d, '
                    'must be one of: %r, %r',
                    i, expected_cols, ALLOWED_COLS, line)
                logger.error(ERROR_PLATE_SIZE + ': ' + msg)
                errors[ERROR_PLATE_SIZE].append(msg)
            in_matrix = True
            plate_matrix = []
            plate_matrices.append(plate_matrix)
            continue
        if in_matrix:
            if row_pattern.match(line):
                row = re.split(r'\s+',line)
                if len(plate_matrix) == 0:
                    logger.debug('recognized row0: %r', line)
                if expected_cols > 0:
                    if len(row[1:]) < expected_cols:
                        msg = (
                            'row: %d, matrix: %d, line: %d, '
                            'cols recognized: %d, expected: %d, %r'
                            % (i, len(plate_matrices), len(plate_matrix)+1,
                                len(row[1:]), expected_cols,row))
                        logger.error(ERROR_COL_SIZE + ': ' + msg)
                        errors[ERROR_COL_SIZE].append(msg)
                    plate_matrix.append(row[1:expected_cols+1])
                else:
                    plate_matrix.append(row[1:])
            else:
                if is_empty_plate_matrix(plate_matrix) is False:
                    if len(plate_matrix) not in ALLOWED_ROWS:
                        msg = (
                            'row: %d, plate_matrix: %d, '
                            'rows found: %r, must be one of %r' 
                            % (i, len(plate_matrices), len(plate_matrix), 
                                ALLOWED_ROWS))
                        logger.error(ERROR_ROW_SIZE + ': ' + msg)
                        errors[ERROR_ROW_SIZE].append(msg)
                        for i,row in enumerate(plate_matrix):
                            logger.error('%r: %r', i, row)
                else:
                    del plate_matrices[-1]
                in_matrix = False
    logger.info('plate matrices: %d', len(plate_matrices))
    
    return (plate_matrices, errors)
    
def read_xlsx(input_file):

    wb = xlrd.open_workbook(file_contents=input_file.read())
    sheets = xlsutils.workbook_sheets(wb)
    plate_matrices = []
    errors = {}
    for sheet in sheets:
        
        (sheet_plate_matrices,sheet_errors) = read_rows(xlsutils.sheet_rows(sheet))
        logger.info('for sheet: %r, %d matrices read', 
            sheet.name, len(sheet_plate_matrices))
        if sheet_plate_matrices:
            plate_matrices += sheet_plate_matrices
            logger.info('plate_matrices: %r', len(plate_matrices))
        if sheet_errors:
            logger.warn('sheet errors: %r: %r', sheet.name, sheet_errors)
            errors[sheet.name] = sheet_errors
    return (plate_matrices, errors)

def read_csv(input_file):

    delimiter = ','
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
        
parser = argparse.ArgumentParser(
    description='Parse raw data files into Python matrices')

parser.add_argument(
    '-f', '--file', required=True,
    help='''Raw data file;
            examples can be found in db/static/test_data/rawdata''')
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
        
        for i,matrix in enumerate(plate_matrices):
            logger.info('matrix #: %r', i)
            logger.info('%r', matrix)
        
  
        