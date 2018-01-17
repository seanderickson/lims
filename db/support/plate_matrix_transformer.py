from __future__ import unicode_literals

import logging
import argparse

import db.support.raw_data_reader as raw_data_reader

from django.db.utils import ProgrammingError
from collections import OrderedDict
import xlsxwriter
from reports import ValidationError
from db.support import lims_utils

logger = logging.getLogger(__name__)

DEBUG = logger.isEnabledFor(logging.DEBUG)

ALLOWED_MATRIX_SIZES = [96,384,1536]

class Collation():
    
    PQ_C_REP_READ = 0
    PQ_REP_C_READ = 1
    C_PQ_REP_READ = 2
    C_REP_PQ_READ = 3
    REP_PQ_C_READ = 4
    REP_C_PQ_READ = 5
    
    ordered_members = ['PQ_C_REP_READ','PQ_REP_C_READ','C_PQ_REP_READ',
            'C_REP_PQ_READ','REP_PQ_C_READ','REP_C_PQ_READ',]
    
    @staticmethod
    def get_members():
        return Collation.ordered_members

    @staticmethod
    def get_value(member_name):
        logger.info('get_value: %r', member_name)
        if member_name.upper() in Collation.ordered_members:
            return Collation.ordered_members.index(member_name.upper())
        else:
            msg = 'must be one of %r' % Collation.ordered_members
            logger.warn('collation ' + msg)
            raise ValidationError(
                key='collation', msg=msg )
    @staticmethod
    def get_member_name(ivalue):
        return Collation.ordered_members[ivalue]

class Counter(object):
    ''' Implement a positional notation Counter:
    - each position contains an ordered array of counting choices 
    - each position has a positional notation value
    '''
    
    def __init__(self, counter_hash):
        ''' counter_hash: ordereddict of counter digits, e.g.
        OrderedDict(
            ('plate', (1,2,3)),
            ('condition', (1,2),
            ('replicate', (1,2),
            ('readout', (1,2),
        )
        '''
        self.counter_hash = counter_hash
    
    def size(self):
        size = 1
        for v in self.counter_hash.values():
            size *= len(v)
        return size
        
    def get_index(self,reading_hash):
        ''' Find the index for a counter "reading" '''
        
        if len(self.counter_hash) != len(reading_hash):
            raise ProgrammingError(
                'counter hash: %r does not match size: index hash %r '
                % (self.counter_hash.keys(), reading_hash.keys()))
        
        indexValue = 0   
        for key, value in reading_hash.items():
            position = self.counter_hash.keys().index(key)
            if value not in self.counter_hash.values()[position]:
                raise ProgrammingError(
                    '%r:%r not in digits: %r' 
                    % (key,value, self.counter_hash.values()[position])
                    )
            positionValue = self.counter_hash.values()[position].index(value)
            for lesserPositionDigits in self.counter_hash.values()[position+1:]:
                positionValue *= len(lesserPositionDigits)
            logger.debug('key: %r:%r, positionValue: %d', key,value, positionValue)
            indexValue += positionValue
                
        return indexValue
    
    def get_readout(self, index):
        ''' Generate the counter "reading" for the index '''
        
        reading_hash = OrderedDict()
        positionIndex = 1
        remainder = index
        for i,(key,digits) in enumerate(self.counter_hash.items()):
            positionValue = 1
            for lesserPositionDigits in self.counter_hash.values()[i+1:]:
                positionValue *= len(lesserPositionDigits)
            positionIndex = remainder/positionValue
            remainder = remainder % positionValue
            reading_hash[key] = digits[positionIndex]

        return reading_hash
    
    def iterate_counter_columns(self):
        '''
        Perform a standard iteration of the counter columns, per plate:
        - as specified by the LIMS functionality: iterate in the order:
        [conditions, replicates, readouts] for each plate.
        FIXME: this depends on the specific LIMS instance members for the counter.
        '''
        counter_hash = self.counter_hash
        for condition in counter_hash['condition']:
            for replicate in counter_hash['replicate']:
                for readout in counter_hash['readout']:
                    yield dict(zip(
                        ('condition','replicate','readout'),
                        (condition,replicate,readout)))
    
    def __repr__(self, *args, **kwargs):
        return '%r' % self.counter_hash
    
def transform(input_matrices, counter, aps, lps):
    
    assert aps in ALLOWED_MATRIX_SIZES, \
        ('assay_plate_size must be one of %r' % ALLOWED_MATRIX_SIZES)
    assert lps in ALLOWED_MATRIX_SIZES, \
        ('library_plate_size must be one of %r' % ALLOWED_MATRIX_SIZES)

    if aps < lps:
        logger.info('convolute matrices')
        factor = lps/aps  
        if factor != 4:
            msg = (
                'Convolute: library_plate_size/assay_plate_size != 4: %d/%d'
                    % (aps, lps))
            raise ValidationError({
                'assay_plate_size': msg,
                'library_plate_size': msg })
        if len(input_matrices) %4 != 0:
            msg = 'Convolute: input matrix array must contain a multiple of 4 members'
            raise ValidationError({
                'assay_plate_size': msg,
                'library_plate_size': msg })
        # Create an adjusted counter to match the input:
        # - add quadrant counter to the right of plate counter
        new_counter_hash = OrderedDict()
        for key,value in counter.counter_hash.items():
            new_counter_hash[key] = value
            if key == 'plate':
                new_counter_hash['quadrant'] = [0,1,2,3]
        counter96 = Counter(new_counter_hash)
        logger.info('counter96: %r', counter96)
        if counter96.size() != len(input_matrices):
            raise ProgrammingError('input_matrices length (%d) must match '
                'the counter length with 4 quadrants: (%d)'
                    % (len(input_matrices), counter96.size()))
        
        
        # - Create blank output matrices
        convoluted_matrices = [
            lims_utils.create_blank_matrix(lps) 
            for x in range(0,len(input_matrices)/4)]

        # Iterate through output (384) matrices and find the 96 matrix values
        # NOTE: could also start by iterating through input matrices
        for index,matrix in enumerate(convoluted_matrices):
            readout = counter.get_readout(index)
            for rownum, row in enumerate(matrix):
                for colnum in range(0,len(row)):
                    input_quadrant = lims_utils.deconvolute_quadrant(
                        lps, aps, 
                        rownum, colnum)
                    readout96 = dict(readout, quadrant=input_quadrant)
                    logger.debug('index: %d, 384 readout: %r, quadrant: %d, 96: %r',
                        index,readout,input_quadrant,readout96)
                    logger.debug('counter96: %r' % counter96.counter_hash)
                    input_index = counter96.get_index(readout96)
                    input_row = lims_utils.deconvolute_row(lps, aps,rownum, colnum)
                    input_col = lims_utils.deconvolute_col(lps, aps,rownum, colnum)
                    logger.debug('find: index: %d, cell: [%d][%d]',
                        input_index,input_row,input_col)
                    row[colnum] = input_matrices[input_index][input_row][input_col]
        
        return convoluted_matrices
        
    elif lps < aps:
        logger.info('deconvolute matrices')
        factor = aps/lps  
        if factor != 4:
            msg = (
                'Deconvolute: assay_plate_size/library_plate_size != 4: %d/%d'
                    % (aps, lps))
            raise ValidationError({
                'assay_plate_size': msg,
                'library_plate_size': msg })
        # Create an adjusted counter to match the input        
        plates = counter.counter_hash.get('plate')
        logger.info('plates: %r', plates)
        if len(plates) % 4 != 0:
            msg = 'Deconvolute: plate count must be a multiple of 4: %d' % len(plates)
            raise ValidationError({
                'plate_ranges': msg })
            
        plates_1536 = OrderedDict()
        for i,plate in enumerate(plates):
            plate_number_1536 = i/4
            if plate_number_1536 not in plates_1536:
                plates_1536[plate_number_1536] = []
            plates_1536[plate_number_1536].append(plate)
        logger.info('plates_1536: %r', plates_1536)
        new_counter_hash = counter.counter_hash.copy()
        new_counter_hash['plate'] = plates_1536.keys()
        counter1536 = Counter(new_counter_hash)

        # Create blank output matrices
        deconvoluted_matrices = [ None for x in range(0,len(input_matrices)*4) ] 
        # Iterate through input (1536) matrices and find the output 384 matrix value
        for index, matrix in enumerate(input_matrices):
            readout1536 = counter1536.get_readout(index)
            plate1536 = readout1536['plate']
            
            # Convert each 1536 plate separately, and find the output matrix position
            output_384_matrices = lims_utils.deconvolute_matrices(
                [matrix], aps, lps)
            
            for quadrant, matrix384 in enumerate(output_384_matrices):
                plate384 = plates_1536[plate1536][quadrant]
                readout384 = dict(readout1536, plate=plate384)
                index384 = counter.get_index(readout384)
                
                deconvoluted_matrices[index384] = matrix384    
        
        return deconvoluted_matrices
    
    else:
        return input_matrices
        
def write_xlsx(workbook, matrices, counter):
    '''
    Write the plate matrices directly to a spreadsheet, in collation order:
    - For testing, does not merge in library well data.
    
    '''
    DEBUG = False or logger.isEnabledFor(logging.DEBUG)

    logger.info('write workbook...')
    counter_hash = counter.counter_hash
    
    counter_readouts = []
    for condition in counter_hash['condition']:
        for replicate in counter_hash['replicate']:
            for readout in counter_hash['readout']:
                counter_readouts.append(dict(zip(
                    ('condition','replicate','readout'),
                    (condition,replicate,readout))))
    
    plate_size = len(matrices[0])*len(matrices[0][0])
    row_size = len(matrices[0])
    col_size = len(matrices[0][0])
    logger.info('row size: %d, col size: %d', row_size, col_size)
    for plate in counter_hash['plate']:
        logger.info('write plate: %r', plate)
        sheet = workbook.add_worksheet(str(plate))
        
        sheet.write_string(0,0,'Plate')
        sheet.write_string(0,1,'Well')
        col_to_matrix_index = []
        for i,counter_readout in enumerate(counter_readouts): 
            
            current_col = 2 + i
            collation_string = '{condition}_{readout}_{replicate}'.format(**counter_readout)
            logger.info('write collation_string: col: %d, %s', current_col, collation_string)
            sheet.write_string(0,current_col,collation_string )

            col_to_matrix_index.append(
                counter.get_index(dict(counter_readout, plate=plate)))

        # NOTE: to support xlsxwriter 'constant_memory': True - optimized write, 
        # rows must be written sequentially        
        for colnum in range(0, col_size):
            for rownum in range(0,row_size):
                output_row = 1 + rownum + colnum * row_size
                
                wellname = lims_utils.get_well_name(rownum, colnum)
                if DEBUG:
                    logger.info('write row: %d: %r', output_row, wellname)
                sheet.write_string(output_row,0,str(plate))
                sheet.write_string(output_row,1,wellname)

                for i,matrix_index in enumerate(col_to_matrix_index):
                    matrix = matrices[matrix_index]
                    current_col = 2 + i
                    val = matrix[rownum][colnum]
                    if DEBUG:
                        logger.info('write output_row: %d, col: %d,  val: %r', 
                            output_row, current_col, str(val))
                    sheet.write_string(output_row, current_col, str(val))


        # NOTE: to support xlsqriter 'constant_memory': True - optimized write,
        # rows must be written sequentially
        # for i,counter_readout in enumerate(counter_readouts): 
        #     
        #     current_col = 2 + i
        #     collation_string = '{condition}_{readout}_{replicate}'.format(**counter_readout)
        #     logger.info('write collation_string: col: %d, %s', current_col, collation_string)
        #     sheet.write_string(0,current_col,collation_string )
        #     
        #     matrix_index = counter.get_index(dict(counter_readout, plate=plate))
        #     logger.info('write matrix: %d', matrix_index)
        #     matrix = matrices[matrix_index]
        #     for colnum in range(0, col_size):
        #         for rownum in range(0,row_size):
        #             output_row = 1 + rownum + colnum * row_size
        #             val = matrix[rownum][colnum]
        #             if DEBUG:
        #                 logger.info('write output_row: %d, col: %d,  val: %r', 
        #                     output_row, current_col, str(val))
        #             sheet.write_string(output_row, current_col, str(val))
    logger.info('write xls finished')
    
def create_matrix_counter(collation, plates, conditions, replicates, readouts):
    
    if collation == Collation.PQ_C_REP_READ:
        return Counter(OrderedDict((
                ('plate', plates),
                ('condition', conditions),
                ('replicate', replicates),
                ('readout', readouts)
            )))
    elif collation == Collation.PQ_REP_C_READ:
        return Counter(OrderedDict((
                ('plate', plates),
                ('replicate', replicates),
                ('condition', conditions),
                ('readout', readouts)
            )))
    elif collation == Collation.C_PQ_REP_READ:
        return Counter(OrderedDict((
                ('condition', conditions),
                ('plate', plates),
                ('replicate', replicates),
                ('readout', readouts)
            )))
    elif collation == Collation.C_REP_PQ_READ:
        return Counter(OrderedDict((
                ('condition', conditions),
                ('replicate', replicates),
                ('plate', plates),
                ('readout', readouts)
            )))
    elif collation == Collation.REP_C_PQ_READ:
        return Counter(OrderedDict((
                ('replicate', replicates),
                ('condition', conditions),
                ('plate', plates),
                ('readout', readouts)
            )))
    elif collation == Collation.REP_PQ_C_READ:
        return Counter(OrderedDict((
                ('replicate', replicates),
                ('plate', plates),
                ('condition', conditions),
                ('readout', readouts)
            )))
    else:
        raise Exception('Collation not implemented %r', collation)

parser = argparse.ArgumentParser(
    description='Parse raw data files into Python matrices and transform')

parser.add_argument(
    '-f', '--file', required=True,
    help='''Raw data file;
            examples can be found in db/static/test_data/rawdata''')

parser.add_argument(
    '-of', '--outputfile', required=True,
    help='''Output file''')

parser.add_argument(
    '-c', '--collation', required=True,
    help='Collation Order one of %r' % Collation.ordered_members  )

parser.add_argument('-p','--plates', nargs='+', 
    help='<Required> plates separated by spaces', required=True)

parser.add_argument('-cond','--conditions', nargs='+', 
    help='<Required> conditions separated by spaces', required=True)

parser.add_argument('-rep','--replicates', nargs='+', 
    help='<Required> replicates separated by spaces', required=True)

parser.add_argument('-read','--readouts', nargs='+', 
    help='<Required> readouts separated by spaces', required=True)

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
    
    collation = Collation.get_value(args.collation)
    counter = create_matrix_counter(
        collation, args.plates, args.conditions, args.replicates, args.readouts)
    
    with open(args.file) as input_file:
        matrices = raw_data_reader.read(input_file, args.file)
        
        if len(matrices) == 0:
            raise Exception('no plates read')
        assay_plate_size = len(matrices[0])*len(matrices[0][0])
        
        library_plate_size = 384
        
        if assay_plate_size != library_plate_size:
            logger.info('transform from %d to %d', assay_plate_size, library_plate_size)
            
            matrices = transform(matrices, counter, assay_plate_size, library_plate_size)

        workbook = xlsxwriter.Workbook(filename=args.outputfile)
        
        write_xlsx(workbook, matrices, counter)
        
        workbook.close()
        
        