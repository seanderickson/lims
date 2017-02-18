from __future__ import unicode_literals
import math
import logging
from db import WELL_NAME_PATTERN, WELL_ID_PATTERN
from reports import ValidationError
logger = logging.getLogger(__name__)


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

def well_name_from_index(index, platesize):
    '''
    @param  index: zero based index, counting from r0c0, r0c1, ... rNcN-1, rNcN
    ''' 
    rows = get_rows(platesize)
    cols = get_cols(platesize)
    name = well_name(int(index/cols),(index%cols))
    return name

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