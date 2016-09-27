from __future__ import unicode_literals
import math
import logging
logger = logging.getLogger(__name__)

def get_rows(platesize):
    return int(math.sqrt(2*platesize/3))

def get_cols(platesize):
    return int(math.sqrt(3*platesize/2))
    
def well_name(col, row):
    name= chr(ord('A')+row) + '{:0>2d}'.format(col)
#     logger.info(str((col,row,name)))
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
    name = well_name((index%cols)+1, int(index/cols))
#     logger.info(str((name, index,rows, cols, platesize)))
    return name