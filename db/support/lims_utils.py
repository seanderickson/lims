# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from collections import defaultdict
from decimal import Decimal
import decimal
from itertools import chain, combinations
import logging
import math
import re

from django.db.utils import ProgrammingError

from db import WELL_NAME_PATTERN, WELL_ID_PATTERN
from reports import ValidationError


## PLATE_SEARCH_LINE_SPLITTING_PATTERN:
# Each plate search range is split
# - on newline
# - on semicolon (*201705, to support GET request urlencoding)
# - on pipe "|" - 20180312, to support GET request urlencoding for compound names
# TODO: consider deprecating semicolon
PLATE_SEARCH_LINE_SPLITTING_PATTERN = re.compile(r'[\n;\|]+')

## PLATE_RANGE_SPLITTING_PATTERN:
# Split a raw plate range input into elements:
# - separated by space or comma, except, 
# - numbers separated by (spaces) and dash interpreted as a plate range
# - quoted strings preserved; to be interpreted as copy names
# - quoted strings may contain spaces and special chars if quoted
QUOTED_WORD_SPLITTING_PATTERN = \
    re.compile(r'''\'.*?\'|\".*?\"|\w+\s+\-\s+\w+|[^\,\s]+''')

logger = logging.getLogger(__name__)

ALLOWED_PLATE_SIZES = [96,384]

ERROR_DUPLICATE_WELLS = 'Duplicate Wells Found'
ERROR_PARSE_EQUALSIGN = 'May only have one equal sign per line'
DISALLOWED_WELL_RANGE_CHARS = re.compile(r'[^\s\-"\',a-zA-Z0-9]+')
ERROR_WELL_RANGE_DISALLOWED_CHARS = \
    'Disallowed chars found: &lt;' + DISALLOWED_WELL_RANGE_CHARS.pattern + '&gt;'    
ERROR_WELL_RANGE_PARTS_UNEQUAL = 'Both values of the range must be the same type'
ERROR_WELL_ROW_OUT_OF_RANGE = 'Row is out of range, max: %r'
ERROR_WELL_COL_OUT_OF_RANGE = 'Col is out of range: max: %r'
ERROR_ENTRY_NOT_RECOGNIZED = 'Unrecognized entry'
def get_rows(platesize):
    return int(math.sqrt(2*platesize/3))

def get_cols(platesize):
    return int(math.sqrt(3*platesize/2))
    
def get_well_name(row, col):
    ''' Convert zero-based row/col indexes to a well name:
    - row zero based index
    - col uses a zero based index to create a 1-based well_name column
    '''
    name= row_to_letter(row) + '{:0>2d}'.format(col+1)
    return name

def well_id(plate, well_name):    
    #return '%05d:%s' % (plate, well_name)
    return '%s:%s' % (str(plate).zfill(5), well_name)

def well_name_from_index_rows_first(index, platesize):
    '''
    DEPRECATED: 
    Well Name from the index, using Row first sorting to count:
    - e.g. A01, A02, A03... A24,B01,B02...
    @param  index: zero based index, counting from r0c0, r0c1, ... rNcN-1, rNcN
    ''' 
    rows = get_rows(platesize)
    cols = get_cols(platesize)
    name = get_well_name(int(index/cols),(index%cols))
    return name

def well_name_from_index(index, platesize):
    '''
    Well Name from the index, using Column first sorting to count:
    - e.g. A01, B01, C01... P01, A02, B02, ...
    - Note: this is the counting method used for Scrensaver 1.
    @param  index: zero based index, counting from c0r1, c1r2, ... cNrN-1, cNrN
    ''' 
    rows = get_rows(platesize)
    cols = get_cols(platesize)
    # Fill by cols just like the Screening lab does: col then row
    name = get_well_name(int(index/cols),(index%cols))
    logger.debug('index: %d, platesize: %d, name: %r', index, platesize, name)
    # Fill in cols by row then col
    #     name = get_well_name((index%rows),int(index/rows))
    return name

def index_from_get_well_name(well_name, platesize):
    
    row_index = well_name_row_index(well_name)
    col_index = well_name_col_index(well_name)
    rows = get_rows(platesize)
    cols = get_cols(platesize)

    if row_index >= rows:
        raise ValidationError(
            key='well_name', 
            msg='%r, row index: %d, is > rows: %d, for platesize: %d' 
                % (well_name, row_index, rows, platesize))
    if col_index >= cols:
        raise ValidationError(
            key='well_name', 
            msg='%r, col index: %d, is > cols: %d, for platesize: %d' 
                % (well_name, col_index, cols, platesize))
    
    # Fill in cols by row then col
    # index = col_index * rows + row_index

    # Fill by cols just like the Screening lab does: col then row
    index = row_index * cols + col_index
    
    if index > platesize:
        raise ValidationError(
            key='well_name', 
            msg='%r, index: %d, is > platesize: %d' 
                % (well_name, index, platesize))
    return index

def row_to_letter(index):
    ''' Convert 0 based row index to letter '''
    
    if index<26: 
        return 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'[index]
    else:
        rem = index%26
        part = int(index/26)-1
        return row_to_letter(part) + row_to_letter(rem)
    
def letter_to_row_index(rowletter):
    ''' Convert letter to zero-based row index '''
    
    if len(rowletter) == 1:
        return 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.index(rowletter.upper())
    else:
        return letter_to_row_index(rowletter[-1]) + 26*(letter_to_row_index(rowletter[:-1])+1)
    
def well_row_col(well_name):
    '''
    @return zero based (row_index,col_index)
    '''
    match = WELL_NAME_PATTERN.match(well_name)
    if not match:
        raise ValidationError(
            key='well_name', 
            msg='%r does not match pattern: %s' % (well_name,WELL_NAME_PATTERN.pattern))
    return (letter_to_row_index(match.group(1)), int(match.group(2))-1)

def well_name_row_index(well_name):
    ''' Convert well_name row letter to zero-based row index '''
    return well_row_col(well_name)[0]

def well_name_col_index(well_name):
    ''' Convert 1-base well_name column to zero-based column index '''
    return well_row_col(well_name)[1]


def well_id_plate_number(well_id):
    ''' Get the plate_number from the well_id '''
    match = WELL_ID_PATTERN.match(well_id)
    if not match:
        raise ValidationError(
            key='well_id', 
            msg='%r Does not match pattern: %s' % (well_id,WELL_ID_PATTERN.pattern))
    return int(match.group(1))

def plate_size_from_plate_type(plate_type):
    '''
    Get the plate size from the current plate_type vocabulary:
    eppendorf_384
    costar_96
    abgene_384
    genetix_384
    marsh_384
    nunc_96
    eppendorf_96
    Note: plate_type must end with the plate size integer for this to work:
    FIXME: plate size determined by magic value embedded in plate_types
    '''
    parts = plate_type.split('_')
    if len(parts) != 2:
        raise ValidationError(
            key='plate_type', 
            msg='not a recognized type: %r' % plate_type)
    plate_size = int(parts[1])
    if plate_size not in ALLOWED_PLATE_SIZES:
        raise ValidationError(
            key='plate_type',
            msg='plate_size: %d for plate_type: %r, not in allowed: %r'
                % (plate_size,plate_type,ALLOWED_PLATE_SIZES))
    return plate_size

def parse_named_well_ranges(raw_data, plate_size):
    logger.info('parse_named_well_ranges(%r, %r) ', raw_data, plate_size)
    
    errors = defaultdict(list)
#     wells_shared_between_ranges_error_msg = 'duplicate wells found in ranges: [%s]'
    named_well_ranges = {}
    
    if not raw_data:
        return (named_well_ranges,errors)
    
    ordinal = 1
    re_strip_quotes = re.compile(r'["\']+')
    
    for well_range in re.split(r'\n', raw_data):
        logger.debug('well_range: %r', well_range)
        well_range = well_range.strip()
        if not well_range:
            continue
        
        unparsed = well_range
        range_to_label = re.split(r'[=]+', well_range)
        logger.debug('range_to_label %r', range_to_label)
        label = ''
        if len(range_to_label) == 2:
            label = re_strip_quotes.sub('',range_to_label[1])
            unparsed = re_strip_quotes.sub('', range_to_label[0])
        elif len(range_to_label) > 2:
            errors[ERROR_PARSE_EQUALSIGN].append(well_range)
        
        named_well_range = named_well_ranges.get(label,None)
        if named_well_range is None:
            named_well_range = {
                'label': label,
                'ordinal': ordinal,
                'text': unparsed,
                'wells': []
                }
            named_well_ranges[label] = named_well_range
            ordinal += 1

        parsed_wells = parse_well_ranges(unparsed, plate_size, errors)
        
        logger.debug('parsed wells: %r, %r', label, parsed_wells)
        named_well_range['wells'].extend(parsed_wells)

    # find duplicates
    duplicate_wells = set()
    for test_label, test_well_range in named_well_ranges.items():
        
        for label, well_range in named_well_ranges.items():
            if test_well_range == well_range:
                continue
            duplicate_wells.update(
                set(test_well_range['wells']) & set(well_range['wells']))
    if duplicate_wells:
        errors[ERROR_DUPLICATE_WELLS] = sorted(duplicate_wells)
    
    if errors:
        logger.info('parse: %r returns errors: %r', raw_data,errors)
    return (named_well_ranges, errors)

def parse_well_ranges(raw_data, plate_size, errors):
    
    logger.info('parse_well_ranges: %r, %r', raw_data, plate_size);
    WELL_PATTERN = re.compile(r'^([a-zA-Z]{1,2})(\d{1,2})$')
    COL_PATTERN = re.compile('^(col:)?\s*(\d{1,2})$')
    ROW_PATTERN = re.compile('^(row:)?\s*([a-zA-Z]{1,2})$')
    
    n_cols = get_cols(plate_size)
    n_rows = get_rows(plate_size)
    wells = []
    
    raw_data = raw_data.strip()
    if not raw_data:
        return wells

    if DISALLOWED_WELL_RANGE_CHARS.match(raw_data):
        errors[ERROR_WELL_RANGE_DISALLOWED_CHARS].append(raw_data)
        return wells
  
    for input in re.split(r'\s*,\s*', raw_data):
        if not input:
            continue
        logger.info('parse well pattern: %r', input)
        parts = re.split(r'\s*-\s*', input)
        
        if len(parts) == 2:
            logger.info('parts: %r', parts)
            if ROW_PATTERN.match(parts[0]):
                match1 = ROW_PATTERN.match(parts[0])
                match2 = ROW_PATTERN.match(parts[1])
                if not match2:
                    errors[ERROR_WELL_RANGE_PARTS_UNEQUAL].append(input)
                start_row = letter_to_row_index(match1.group(2))
                stop_row = letter_to_row_index(match2.group(2))
                
                row_range = sorted([start_row,stop_row])
                logger.info('input: %r, row_range: %r', parts, row_range)
                if row_range[1] >= n_rows:
                    errors[ERROR_WELL_ROW_OUT_OF_RANGE%row_to_letter(n_rows-1)].append(input)
#                     errors.add('row is out of range: %r, max: %r', 
#                         input, row_to_letter(n_rows-1))
                    continue
                for i in range(0,n_cols):
                    for j in range(row_range[0],row_range[1]+1):
                        wells.append(get_well_name(j,i))
            elif COL_PATTERN.match(parts[0]):
                match1 = COL_PATTERN.match(parts[0])
                match2 = COL_PATTERN.match(parts[1])
                if not match2:
                    errors[ERROR_WELL_RANGE_PARTS_UNEQUAL].append(input)
                start_col = int(match1.group(2))-1
                stop_col = int(match2.group(2))-1
                
                col_range = sorted([start_col, stop_col])
                logger.info('col_range: %r, %r', parts, col_range)
                if col_range[1] >= n_cols:
                    errors[ERROR_WELL_COL_OUT_OF_RANGE%n_cols].append(input)
#                     errors.add('col is out of range: %r, max: %r', 
#                         input, n_cols)
                    continue
                for i in range(col_range[0],col_range[1]+1):
                    for j in range(0,n_rows):
                        wells.append(get_well_name(j,i))
            elif WELL_PATTERN.match(parts[0]):
                match1 = WELL_PATTERN.match(parts[0])
                match2 = WELL_PATTERN.match(parts[1])
                if not match2:
                    errors[ERROR_WELL_RANGE_PARTS_UNEQUAL].append(input)
                start_row = letter_to_row_index(match1.group(1))
                stop_row = letter_to_row_index(match2.group(1))
                start_col = int(match1.group(2))-1
                stop_col = int(match2.group(2))-1
                col_range = sorted([start_col, stop_col])
                if col_range[1] >= n_cols:
                    errors[ERROR_WELL_COL_OUT_OF_RANGE%n_cols].append(input)
#                     errors.add('col is out of range: %r, max: %r', 
#                         input, n_cols)
                    continue
                row_range = sorted([start_row,stop_row])
                logger.info('well range: %r, %r, %r', parts,row_range,col_range)
                if row_range[1] >= n_rows:
                    errors[ERROR_WELL_ROW_OUT_OF_RANGE%row_to_letter(n_rows-1)].append(input)
#                     errors.add('row is out of range: %r, max: %r', 
#                         input, row_to_letter(n_rows-1))
                    continue
                for i in range(col_range[0],col_range[1]+1):
                    for j in range(row_range[0],row_range[1]+1):
                        wells.append(get_well_name(j,i))
            else:
                errors[ERROR_ENTRY_NOT_RECOGNIZED].append(input)
#                 errors.append('unrecognized block range entry: %r'% input)
        elif len(parts) == 1:
            input = parts[0]
            if ROW_PATTERN.match(input):
                row = letter_to_row_index(ROW_PATTERN.match(input).group(2))
                if row >= n_rows:
                    errors[ERROR_WELL_ROW_OUT_OF_RANGE%row_to_letter(n_rows-1)].append(input)
#                     errors.add('row is out of range: %r, max: %r', 
#                         input, row_to_letter(n_rows-1))
                    continue
                for i in range(0,n_cols):
                    wells.append(get_well_name(row,i))
            elif COL_PATTERN.match(input):
                col = int(COL_PATTERN.match(input).group(2))-1
                if col >= n_cols:
                    errors[ERROR_WELL_COL_OUT_OF_RANGE%n_cols].append(input)
#                     errors.add('col is out of range: %r, max: %r', 
#                         input, n_cols)
                    continue
                for j in range(0,n_rows):
                    wells.append(get_well_name(j,col))
            elif WELL_PATTERN.match(input):
                row = letter_to_row_index(WELL_PATTERN.match(input).group(1))
                col = int(WELL_PATTERN.match(input).group(2))-1
                if row >= n_rows:
                    errors[ERROR_WELL_ROW_OUT_OF_RANGE%row_to_letter(n_rows-1)].append(input)
#                     errors.add('row is out of range: %r, max: %r', 
#                         input, row_to_letter(n_rows-1))
                    continue
                if col >= n_cols:
                    errors[ERROR_WELL_COL_OUT_OF_RANGE%n_cols].append(input)
#                     errors.add('col is out of range: %r, max: %r', 
#                         input, n_cols)
                    continue
                wells.append(get_well_name(row,col))
            else:
                errors[ERROR_ENTRY_NOT_RECOGNIZED].append(input)
#                 errors.append('unrecognized single specifier entry: %r'% input)
        else:
            errors[ERROR_ENTRY_NOT_RECOGNIZED].append(input)
#             errors.append('unrecognized entry: %r' % input)
    return sorted(wells)
  
   

def parse_wells_to_leave_empty(wells_to_leave_empty, plate_size):
    '''
    TODO: replace with parse_well_ranges
    Parse the wells to leave empty field of the Cherry Pick Request.
    '''
    
    
    logger.info('raw wells_to_leave_empty: %r, plate_size: %r', 
        wells_to_leave_empty, plate_size)

    ncols = get_cols(plate_size)
    nrows = get_rows(plate_size)
    row_pattern = re.compile('row:\s*([a-zA-Z]{1,2})', flags=re.IGNORECASE)
    col_pattern = re.compile(r'col:\s*(\d{1,2})', flags=re.IGNORECASE)
    selections = re.split(r'\s*,\s*', wells_to_leave_empty)
    new_selections = []
    for selection in selections:
        colmatch = col_pattern.match(selection)
        if colmatch: 
            col = int(colmatch.group(1))
            if col > ncols:
                raise ValidationError(
                    key='wells_to_leave_empty',
                    msg='column out of range: %d, %r' % (col, selection))
            new_selections.append('Col:%d' % col)
            continue
        rowmatch = row_pattern.match(selection)
        if rowmatch:
            row = letter_to_row_index(rowmatch.group(1))         
            if row >= nrows:
                raise ValidationError(
                    key='wells_to_leave_empty',
                    msg='row out of range: %r, %r' 
                        % (rowmatch.group(1), selection))
            new_selections.append('Row:%s' % rowmatch.group(1).upper())
            continue
        wellmatch = WELL_NAME_PATTERN.match(selection)
        if wellmatch:
            new_selections.append(selection.upper())
            continue
        
        raise ValidationError(
            key='wells_to_leave_empty',
            msg='unrecognized pattern: %r' % selection)
    logger.debug('new wells_to_leave_empty selections: %r', new_selections)

    decorated = []
    for wellname in new_selections:
        if 'Col:' in wellname:
            decorated.append((1,int(wellname.split(':')[1]),wellname))
        elif 'Row:' in wellname:
            decorated.append((2,wellname.split(':')[1],wellname))
        else:
            match = WELL_NAME_PATTERN.match(wellname)
            decorated.append((match.group(1),match.group(1), wellname))
    new_wells_to_leave_empty = [x[2] for x in sorted(decorated)]
    logger.info('wells_to_leave_empty: %r', new_wells_to_leave_empty)
    
    return new_wells_to_leave_empty

def assay_plate_available_wells(wells_to_leave_empty, plate_size):
    
    if plate_size not in ALLOWED_PLATE_SIZES:
        raise ValidationError(
            key='plate_size',
            msg=('plate_size: %d for assay_plate_available_wells, '
                'not in allowed: %r'
                % (plate_size,ALLOWED_PLATE_SIZES)))
    available_wells = []
    
    # Parse and sanitize the wells_to_leave_empty (should be already done)
    wells_to_leave_empty_list = []
    if wells_to_leave_empty:
        wells_to_leave_empty_list = parse_wells_to_leave_empty(
            wells_to_leave_empty, plate_size)
    row_specifier = 'Row:%s'
    col_specifier = 'Col:%d'
    for i in range(0,plate_size):
        well_name = well_name_from_index(i, plate_size)
        wellmatch = WELL_NAME_PATTERN.match(well_name)
        row = wellmatch.group(1)
        col = int(wellmatch.group(2))
        if row_specifier % row in wells_to_leave_empty_list:
            continue
        if col_specifier % col in wells_to_leave_empty_list:
            continue
        if well_name in wells_to_leave_empty_list:
            continue
        available_wells.append(well_name)
    return available_wells
    
def find_minimal_satisfying_set(complete_set,instance_sets):
    '''
    Given the complete set, find the minimal subset that intersects to form a 
    non-null set (at least one overlap) with all non-empty instance_sets
    '''    
    def powerset(iterable):
        "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        s = list(iterable)
        return set(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))
    
    # satisfying_sets = []
    min_satisfying_set = None
    # Look at all the powersets; until the first viable is found
    # Sort by len,name to find the first
    for subset in sorted(powerset(sorted(complete_set)), key=lambda x: (len(x),str(x))):
        logger.debug('consider: %r', subset)
        found = True
        for instance_set in instance_sets:
            if instance_set:
                if not set(subset)&set(instance_set):
                    found = False
                    break
        if found is True:
            logger.debug('found satisfying set: %r', subset)
            min_satisfying_set = subset
            break
    return min_satisfying_set

def get_siunit(default_unit_value=1e-6):
    '''
    Return the best match SI Unit symbol for the default_unit_value, such that:
    default_unit_value can be represented a number between 1 and 1000;
    (best_match_symbol_val)<=default_unit_value<(next_higher_symbol_val)
    '''
    siunits = [
      ['T', 1e12],
      ['G', 1e9],
      ['M', 1e6],
      ['k', 1e3],
      ['', 1],
      ['m', 1e-3,],
      ['μ', 1e-6,],
      ['n', 1e-9 ],
      ['p', 1e-12 ]
      ]
    for symbol,val in siunits:
        if val <= default_unit_value:
            return symbol
    return None

def convert_decimal(
    raw_val, default_unit=1e-6, decimals=1, multiplier=None):
    '''
    Convert a decimal by scaling to the default unit and rounding to the 
    given decimals (decimal digits), optionally multiplying by a multiplier.
    
    @param default_unit adjust raw_val to the "default_unit" 
    (as defined in the "display_options")
    - e.g. if default_unit = 1e-6:
        adjust the raw_val = raw_val.scaleb(6)
    - e.g. if default_unit = 1e6:
        adjust the raw_val = raw_val.scaleb(-6)
    
    @param decimals digits of precision to apply
    @param multiplier (Note: only the exponent of the multiplier is used, so
    only powers of 10 may be used)
    '''
    assert decimals >= 0, 'decimals must be >= 0'
    
    # get the scale (exponent) of the default unit
    # negate the scale for use with Decimal.scaleb()
    scale = -Decimal(str(default_unit)).adjusted()
    if multiplier is not None:
        # get the scale (exponent) of the multiplier
        multiplier = Decimal(str(multiplier)).adjusted()
        if multiplier != 0:
            scale = scale+multiplier
    decimals = Decimal('1e-%d'%int(decimals))
    val = Decimal(raw_val)
    if scale != 0:
        val = val.scaleb(scale)
    val = val.quantize(decimals, decimal.ROUND_HALF_UP)
    
    return val

def find_ranges(list_of_numbers):
    list_of_numbers = sorted(list_of_numbers)
    ranges = []
    range = []
    for num in list_of_numbers:
        if not range:
            range = [num,num]
        else:
            if range[1] < num-1:
                ranges.append(range)
                range = [num,num]
            else:
                range[1] = num
    if range:
        ranges.append(range)
    ranges = map(lambda x: 
        str(x[0]) if x[0]==x[1]
        else '%s-%s' % (x[0],x[1]), ranges)
    logger.debug('found ranges: %r for %r', ranges, list_of_numbers)
    return ranges

# Transform plate size by "quadrants"

def deconvolute_quadrant(source_ps, dest_ps, row, col):
    '''Map (1536,384)-well row,col to (384,96)-well output quadrant '''

    factor = source_ps/dest_ps
    if factor != 4:
        raise ProgrammingError('Deconvolute may only be used for '
            'source_ps/dest_ps == 4: %d/%d'
            % (source_ps, dest_ps))
        
    return col%(factor/2) +  (row%(factor/2))*(factor/2);

def deconvolute_row(source_ps, dest_ps,row, col):
    '''Map (1536,384-)well input row to (384-)well output row '''
    
    dest_matrix_number = deconvolute_quadrant(
            source_ps, dest_ps, row, col)
    factor = source_ps/dest_ps  
    if factor != 4:
        raise ProgrammingError('Deconvolute may only be used for '
            'source_ps/dest_ps == 4: %d/%d'
            % (source_ps, dest_ps))
    return row/(factor/2)+ row%(factor/2)-dest_matrix_number/(factor/2);

def deconvolute_col(source_ps, dest_ps, row, col):
    '''Map (1536,384)-well input col to (384,96)-well output col'''
    
    dest_matrix_number = deconvolute_quadrant(
            source_ps, dest_ps, row, col);
    factor = source_ps/dest_ps  
    if factor != 4:
        raise ProgrammingError('Deconvolute may only be used for '
            'source_ps/dest_ps == 4: %d/%d'
            % (source_ps, dest_ps))
    return col/(factor/2)+ col%(factor/2)-dest_matrix_number%(factor/2);

def convolute_row(source_ps, dest_ps, source_matrix_quadrant,row):
    '''Map (96,384)-well input row to (384,1536)-well output row '''
    
    factor = dest_ps/source_ps;  
    if factor != 4:
        raise ProgrammingError('convolute may only be used for '
            'dest_ps/source_ps == 4: %d/%d' % (dest_ps,source_ps))
    if source_matrix_quadrant not in [0,1,2,3]:
        raise ProgrammingError('source_matrix_quadrant must be 0<=n<4: %r'
            % source_matrix_quadrant)
    return row * factor/2 + source_matrix_quadrant/(factor/2)

def convolute_col(source_ps, dest_ps, source_matrix_quadrant, col):
    ''''Map (96,384)-well input col to (384,1536)-well output col '''

    factor = dest_ps/source_ps;  
    if factor != 4:
        raise ProgrammingError('convolute may only be used for '
            'dest_ps/source_ps == 4: %d/%d'% (dest_ps,source_ps))
    return col * factor/2 + source_matrix_quadrant%(factor/2)

def convolute_well(source_ps, dest_ps, wellname):
    '''
    Map a list of wells from (96,384) plate format to (384,1536) plate format.
    '''
    factor = dest_ps/source_ps  
    if factor != 4:
        raise ProgrammingError('convolute may only be used for '
            'dest_ps/source_ps == 4: %d/%d'% (dest_ps,source_ps))
    convoluted_wells = []
    (row,col) = well_row_col(wellname)
    for quadrant in range(0,factor):
        new_row = convolute_row(source_ps, dest_ps, quadrant, row)
        new_col = convolute_col(source_ps, dest_ps, quadrant, col)
        convoluted_wells.append(get_well_name(new_row,new_col))
    return convoluted_wells

def convolute_wells(source_ps, dest_ps, wells):
    '''
    Map a list of wells from (96,384) plate format to (384,1536) plate format.
    '''
    convoluted_wells = []
    factor = dest_ps/source_ps  
    if factor != 4:
        raise ProgrammingError('convolute may only be used for '
            'dest_ps/source_ps == 4: %d/%d'
            % (source_ps, dest_ps))
    for wellname in  wells:
        convoluted_wells.extend(convolute_well(source_ps, dest_ps, wellname))
    logger.debug('wells convoluted: %r, %r, %r, %r', 
        wells, source_ps,dest_ps, convoluted_wells)
    return convoluted_wells

def deconvolute_well(source_ps, dest_ps, wellname):
    '''
    Map a well from (384,1536) plate format to (96,384) plate format.
    '''
    (row,col) = well_row_col(wellname)
    quadrant = deconvolute_quadrant(source_ps, dest_ps, row, col)
    new_row = deconvolute_row(source_ps, dest_ps, row, col)
    new_col = deconvolute_col(source_ps, dest_ps, row, col)
    return get_well_name(new_row,new_col)

def deconvolute_wells(source_ps, dest_ps, wells):
    '''
    Map a list of wells from (384,1536) plate format to (96,384) plate format.

    @return map of the plateQuadrant->deconvolutedWells where plateQuadrant
    is in [0,1,2,3]
    '''
    
    deconvoluted_plate_quadrant_wells = [[] for q in range(0,4)]
    for wellname in wells:
        (row,col) = well_row_col(wellname)
        
        quadrant = deconvolute_quadrant(source_ps, dest_ps, row, col)
        new_row = deconvolute_row(source_ps, dest_ps, row, col)
        new_col = deconvolute_col(source_ps, dest_ps, row, col)
        
        deconvoluted_plate_quadrant_wells[quadrant].append(get_well_name(new_row,new_col))

    logger.info('wells deconvoluted: %r, %d, %d, %r', 
        wells, source_ps,dest_ps, deconvoluted_plate_quadrant_wells)
    return deconvoluted_plate_quadrant_wells


def create_blank_matrix(ps):
    matrix = []
    for rownum in range(0,get_rows(ps)):
        matrix.append([None]*get_cols(ps))
    return matrix    

def convolute_matrices(input_matrices, source_ps, dest_ps):
    '''Map (96,384)-well input matrices to (384,1536)-well output matrices '''
    
    factor = dest_ps/source_ps
    if factor != 4:
        raise ProgrammingError('convolute may only be used for '
            'dest_ps/source_ps == 4: %d/%d' % (source_ps, dest_ps))
    
    assert len(input_matrices)%4 == 0, 'input_matrices count must be a factor of 4'
      
    output_size = len(input_matrices)/4
    
    logger.info('convolute_matrices: convert %d input matrices to %d output',
        len(input_matrices),output_size)
    
    output_matrices = []
    for i in range(0,output_size):
        output_matrix = create_blank_matrix(dest_ps)
        output_matrices.append(output_matrix)
        
        start = i*4
        end = start+4
        quadrant_matrices = input_matrices[start:end]
        # fill the matrix    
        for quadrant, matrix in enumerate(quadrant_matrices):
            for rownum, row in enumerate(matrix):
                for colnum, val in enumerate(row):
                    dest_row = convolute_row(
                        source_ps, dest_ps, quadrant, rownum)
                    dest_col = convolute_col(
                        source_ps, dest_ps, quadrant, colnum)
                    output_matrix[dest_row][dest_col] = val 

    return output_matrices

def deconvolute_matrices(input_matrices, source_ps, dest_ps):
    '''Map (1536,384)-well input matrices to (384,96)-well output matrices '''

    factor = source_ps/dest_ps
    if factor != 4:
        raise ProgrammingError('deconvolute may only be used for '
            'source_ps/dest_ps == 4: %d/%d' % (source_ps, dest_ps))

    logger.debug('deconvolute_matrices: convert %d input matrices to %d output',
        len(input_matrices),len(input_matrices)*factor)

    output_matrices = []
    
    for input_matrix in input_matrices:
        new_output_matrices = [ create_blank_matrix(dest_ps) for i in range(0,4)]
        output_matrices += new_output_matrices
        for rownum, row in enumerate(input_matrix):
            for colnum, val in enumerate(row):
                output_quadrant = deconvolute_quadrant(
                    source_ps, dest_ps, rownum, colnum)
                output_row = deconvolute_row(source_ps, dest_ps, rownum, colnum)
                output_col = deconvolute_col(source_ps, dest_ps, rownum, colnum)
                
                new_output_matrices[output_quadrant][output_row][output_col] = val

    return output_matrices

