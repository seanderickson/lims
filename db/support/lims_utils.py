from __future__ import unicode_literals

import decimal
from decimal import Decimal
from itertools import chain, combinations
import logging
import math
import re

from db import WELL_NAME_PATTERN, WELL_ID_PATTERN
from reports import ValidationError

QUOTED_WORD_PATTERN = re.compile(r'''(\'.*?\'|\".*?\"|[^\s,;]+)''')

logger = logging.getLogger(__name__)

ALLOWED_PLATE_SIZES = [96,384]

def get_rows(platesize):
    return int(math.sqrt(2*platesize/3))

def get_cols(platesize):
    return int(math.sqrt(3*platesize/2))
    
def well_name(row, col):
    ''' Convert zero-based row/col indexes to a well name:
    - row is [A-Z]
    - col uses a zero based index to create a 1-based well_name column
    '''
    name= chr(ord('A')+row) + '{:0>2d}'.format(col+1)
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
    name = well_name(int(index/cols),(index%cols))
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
    name = well_name((index%rows),int(index/rows))
    return name

def index_from_well_name(well_name, platesize):
    
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
    
    index = col_index * rows + row_index
    
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
        return letter_to_row(rowletter[-1]) + 26*letter_to_row(rowletter[:-1])
    
def well_name_row_index(well_name):
    ''' Convert well_name row letter to zero-based row index '''
    match = WELL_NAME_PATTERN.match(well_name)
    if not match:
        raise ValidationError(
            key='well_name', 
            msg='%r does not match pattern: %s' % (well_name,WELL_NAME_PATTERN.pattern))
    return letter_to_row_index(match.group(1));

def well_name_col_index(well_name):
    ''' Convert 1-base well_name column to zero-based column index '''
    
    match = WELL_NAME_PATTERN.match(well_name)
    if not match:
        raise ValidationError(
            key='well_name', 
            msg='%r does not match pattern: %s' % (well_name,WELL_NAME_PATTERN.pattern))
    return int(match.group(2))-1

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
            msg='not a recognized type: %r' % assay_plate_type)
    plate_size = int(parts[1])
    if plate_size not in ALLOWED_PLATE_SIZES:
        raise ValidationError(
            key='plate_type',
            msg='plate_size: %d for plate_type: %r, not in allowed: %r'
                % (plate_size,plate_type,ALLOWED_PLATE_SIZES))
    return plate_size

def parse_wells_to_leave_empty(wells_to_leave_empty, plate_size):
    
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