from __future__ import unicode_literals

import argparse
from collections import OrderedDict
import json
import logging
import re

import xlrd
from xlsxwriter.utility import xl_col_to_name

from db import WELL_NAME_PATTERN
import db.schema as SCHEMA
from db.support import lims_utils
from db.support.data_converter import default_converter
from reports import ParseError, ValidationError, LIST_DELIMITER_SQL_ARRAY
from reports.serialize import parse_val
from reports.serialize.xlsutils import sheet_cols, sheet_rows, \
    workbook_sheets, generic_xls_write_workbook


logger = logging.getLogger(__name__)

DEBUG_IMPORTER = False or logger.isEnabledFor(logging.DEBUG)

PARTITION_POSITIVE_MAPPING = \
    SCHEMA.VOCAB.resultvalue.partitioned_positive.get_dict()

CONFIRMED_POSITIVE_MAPPING = \
    SCHEMA.VOCAB.resultvalue.confirmed_positive.get_dict()
ASSAY_WELL_CONTROL_TYPES = {
    'p': SCHEMA.VOCAB.assaywell.control_type.ASSAY_POSITIVE_CONTROL,
    'n': SCHEMA.VOCAB.assaywell.control_type.ASSAY_CONTROL,
    # 20180315 - removed - "shared" is not used in the current database
    # 's': SCHEMA.VOCAB.assaywell.control_type.SHARED,
    'o': SCHEMA.VOCAB.assaywell.control_type.OTHER,
}
META_MAP = OrderedDict((
    ('Screen Number', 'screen_facility_id'),
))
DC_SCHEMA = SCHEMA.DATA_COLUMN
DATA_COLUMN_FIELD_MAP = OrderedDict((
    ('"Data" Worksheet Column',  'data_worksheet_column'),
    ('name',  DC_SCHEMA.NAME),
    ('title',  DC_SCHEMA.TITLE),
    ('data type',  DC_SCHEMA.DATA_TYPE),
    ('decimal places',  DC_SCHEMA.DECIMAL_PLACES),
    ('description',  DC_SCHEMA.DESCRIPTION),
    ('replicate number',  DC_SCHEMA.REPLICATE_ORDINAL),
    ('time point',  DC_SCHEMA.TIME_POINT),
    ('assay readout type',  DC_SCHEMA.ASSAY_READOUT_TYPE),
    ('if derived, how?',  DC_SCHEMA.HOW_DERIVED),
    ('if derived, from which columns?',  DC_SCHEMA.DERIVED_FROM_COLUMNS),
    ('primary or follow up?',  DC_SCHEMA.IS_FOLLOW_UP_DATA),
    ('which assay phenotype does it belong to?',  DC_SCHEMA.ASSAY_PHENOTYPE),
    ('comments',  DC_SCHEMA.COMMENTS),
    ('channel', DC_SCHEMA.CHANNEL),
    ('time point ordinal', DC_SCHEMA.TIME_POINT_ORDINAL),
    ('zdepth ordinal', DC_SCHEMA.ZDEPTH_ORDINAL),
    ('clo id', 'clo_id'), # 20180315 - Not used
    ('screen_facility_id', DC_SCHEMA.SCREEN_FACILITY_ID),
))

DATA_TYPE = SCHEMA.VOCAB.datacolumn.data_type
DATA_TYPE_VALUES = DATA_TYPE.get_ordered_dict().values()

RESULT_VALUE_FIELD_MAP = OrderedDict((
    ('plate', 'plate_number'),
    ('well', 'well_name'),
    ('type', 'assay_well_control_type'),
    ('exclude', 'exclude'),
))

# TODO: 20180315 - use db.schema to map constants for fields

# Extra columns that may appear in the exported report from the API
REPORTING_NON_RV_COLUMNS = [
    'well_id', 'plate_number','well_name','screen_facility_id', 
    'assay_well_control_type'
]

def data_column_field_mapper(fields):
    mapped_keys = []
    for key in fields:
        key = key.lower()
        mapped_key = None
        if key in DATA_COLUMN_FIELD_MAP:
            mapped_key = DATA_COLUMN_FIELD_MAP[key]
        else:
            for k in DATA_COLUMN_FIELD_MAP.values():
                if default_converter(key) == k:
                    mapped_key = k
                    break
        if not mapped_key:
            raise ParseError(
                key=key, 
                msg=('Key %r is not in the recognized datacolumn fields: %r'
                    % (key,DATA_COLUMN_FIELD_MAP.keys())))
        mapped_keys.append(mapped_key)
    return mapped_keys

def parse_columns(columns_sheet):
    logger.info('parsing columns sheet: %r', columns_sheet.name)
    columns = data_column_generator(sheet_cols(columns_sheet))
    parsed_cols = OrderedDict()
    errors = {}
    for i,column in enumerate(columns):
        parsed_col = {
            'is_derived': False,
            'is_follow_up_data': False,
            'ordinal': i
        }
        logger.debug('parsing column: %r', column['data_worksheet_column'])
        parsed_cols[column['data_worksheet_column']] = parsed_col
        for key,val in column.items():
            if key == 'is_follow_up_data':
                parsed_col[key] = ( val and val.lower() == 'follow up')
            elif key == 'data_type':
                val = default_converter(val)
                # handle validation errors in the api
                if val not in DATA_TYPE_VALUES:
                    key = '%s:%s' % (column['data_worksheet_column'],'data_type')
                    errors[key] = 'val: %r must be one of %r' % (val,DATA_TYPE_VALUES)
                parsed_col[key] = val
            elif key == 'assay_readout_type':
                parsed_col[key] = default_converter(val)
            else:
                if key == 'how_derived':
                    parsed_col['is_derived'] = ( 
                        val is not None and val.strip() is not '' )
                parsed_col[key] = val
        
        if parsed_col.get('decimal_places') is not None:
            try:
                key = '%s:%s' % (column['data_worksheet_column'],'data_type')
                column['decimal_places'] = parse_val(
                    column['decimal_places'],key,'integer')
            except ValidationError, e:
                errors.update(e.errors)
        logger.debug('parsed_col: %r', parsed_col)
    if errors:
        raise ValidationError(errors={'Data Columns': errors})
    
    logger.debug('parsed cols: %r', parsed_cols)
    return parsed_cols
        
def result_value_field_mapper(header_row, parsed_columns):
    if DEBUG_IMPORTER:
        logger.info('map result value header row... %r', parsed_columns.keys())
    mapped_row = []
    header_row = [x for x in header_row]
    for i,value in enumerate(header_row):
        if value.lower() in RESULT_VALUE_FIELD_MAP:
            mapped_row.append(RESULT_VALUE_FIELD_MAP[value.lower()])
        else:
            colname = xlrd.book.colname(i)
            mapped_row.append(colname)
    if DEBUG_IMPORTER:
        logger.info('mapped header row: %r', dict(zip(header_row, mapped_row)))
    unmapped = [key for key,value in RESULT_VALUE_FIELD_MAP.items() 
        if value not in mapped_row]
    if unmapped:
        msg=('Missing fields: %s in result values header row: %r'
            % (', '.join(['"%s"'%f for f in unmapped]), 
               ', '.join(['"%s"'% f for f in header_row])))
        logger.info(msg)
        raise ParseError(key='Header row',msg=msg)
    if DEBUG_IMPORTER:
        logger.info('mapped result value header row: %r', mapped_row) 
    return mapped_row
        
def parse_result_values(parsed_columns, sheets):
    logger.info('parse_result_values...')
    well_ids = set()
    parse_error = None
    for sheet in sheets:
        logger.info('parse result values sheet: %r...', sheet.name)
    
        rows = sheet_rows(sheet)
        try:
            header_row = result_value_field_mapper(rows.next(), parsed_columns)
        except ValidationError, e:
            logger.exception('error: %r', e)
            if not parse_error:
                parse_error = ParseError(errors={})
            if not sheet.name in parse_error.errors:
                parse_error.errors[str(sheet.name)] = {}
            parse_error.errors[sheet.name] = e.errors
            continue
        logger.info('output result values...')
        for i,row in enumerate(rows):
            try:
                result = parse_result_row(
                    i,parsed_columns,dict(zip(header_row,row)))
                if DEBUG_IMPORTER:
                    logger.info('parsed row: %d: %r',i,  result)
                if result['well_id'] in well_ids:
                    raise ParseError(
                        key=result['well_id'],
                        msg='duplicate')
                well_ids.add(result['well_id'])
                yield result
            except ValidationError,e:
                logger.exception('parse error: %r', e)
                if not parse_error:
                    parse_error = ParseError(errors={})
                if not sheet.name in parse_error.errors:
                    parse_error.errors[sheet.name] = {}
                parse_error.errors[sheet.name].update(e.errors)
    if parse_error:
        raise parse_error
    
def parse_result_row(i,parsed_columns,result_row):    
    
    logger.debug(
        'parse result row: %d, %r:  %r', i, parsed_columns.keys(), result_row)
    
    meta_columns = RESULT_VALUE_FIELD_MAP.values()
    parsed_row = {}
    excluded_cols = []
    
    meta_key = 'plate_number'
    val = result_row[meta_key]
    logger.debug('plate value to parse: %r', val)
    plate_number = parse_val(val, meta_key, 'integer')
    meta_key = 'well_name'
    val = result_row[meta_key]
    if WELL_NAME_PATTERN.match(val):
        wellname = val
    else:
        raise ParseError(
            key=i, 
            msg=('Well_name val %r does not follow the pattern: %r'
            % (val, WELL_NAME_PATTERN.pattern))) 
    parsed_row['well_id'] = \
        '%s:%s' % (str(plate_number).zfill(5), wellname)
    
    meta_key = 'assay_well_control_type'
    val = result_row.get(meta_key)
    parsed_row[meta_key] = None
    if val is not None:
        if val.lower() in ASSAY_WELL_CONTROL_TYPES:
            parsed_row[meta_key] = \
                ASSAY_WELL_CONTROL_TYPES[val.lower()]
        else:
            msg = ('%s: val %r is not one of the choices: %r'
                % (meta_key, val, ASSAY_WELL_CONTROL_TYPES))
            logger.error(msg)
            raise ValidationError(key=parsed_row['well_id'], msg=msg)

    meta_key = 'exclude'
    val = result_row.get(meta_key)
    if val is not None:
        if val.lower() == 'all':
            excluded_cols = parsed_columns.keys()
        else:
            excluded_cols = [x.strip().upper() for x in val.split(',')]
            unknown_excluded_cols = (
                set(excluded_cols) - set(parsed_columns.keys()))
            if unknown_excluded_cols:
                raise ValidationError(
                    key = parsed_row['well_id'],
                    msg = 'unknown excluded cols: %r' % unknown_excluded_cols )
            parsed_row[meta_key] = excluded_cols
            
    for colname, raw_val in result_row.items():
        logger.debug('colname: %r, raw_val: %r', colname, raw_val)
        if colname in meta_columns:
            continue
        if colname not in parsed_columns:
            # NOTE: this is no longer an error, as the result value sheet may
            # contain extra columns (selected by user on output)
            logger.debug(
                'result value column %r is not in recognized columns: %r', 
                colname, parsed_columns.keys())
            parsed_row[colname] = raw_val
            continue
        column = parsed_columns[colname]
        if raw_val is None:
            # 20180315 - verified with DJW, default values for
            # positive indicator columns
            if column['data_type'] == DATA_TYPE.BOOLEAN_POSITIVE:
                raw_val = False
            elif column['data_type'] == DATA_TYPE.PARTITIONED_POSITIVE:
                raw_val = 'NP'
            elif column['data_type'] == DATA_TYPE.CONFIRMED_POSITIVE:
                raw_val = 'NT'
            else:
                continue
        
        key = '%s-%s' % (parsed_row['well_id'],colname)
        parsed_row[colname] = raw_val
        
#         if column['data_type']  in ['numeric','decimal','integer']:
        if column['data_type']  in DATA_TYPE.numeric_types:
            if  column['decimal_places'] > 0:
                # parse, to validate only; use decimal for final parsing
                parse_val(raw_val, key, 'float')
            else:
                parsed_row[colname] = parse_val(raw_val, key, 'integer')
                
        elif column['data_type'] == DATA_TYPE.PARTITIONED_POSITIVE:
            val = raw_val.upper()
            if val not in PARTITION_POSITIVE_MAPPING:
                raise ValidationError(
                    key=key, 
                    msg='val: %r must be one of %r'
                        % (raw_val, PARTITION_POSITIVE_MAPPING.keys()))
            parsed_row[colname] = val
        elif column['data_type'] == DATA_TYPE.CONFIRMED_POSITIVE:
            val = raw_val.upper()
            if val not in CONFIRMED_POSITIVE_MAPPING:
                raise ValidationError(
                    key=key, 
                    msg='val: %r must be one of %r'
                        % (raw_val, CONFIRMED_POSITIVE_MAPPING.keys()))
            parsed_row[colname] = val
        elif column['data_type'] == DATA_TYPE.BOOLEAN_POSITIVE:
            val = parse_val(raw_val,key,'boolean')
            parsed_row[colname] = val
        logger.debug('parsed_row: %r', parsed_row)
    return parsed_row

def data_column_generator(data_columns):
    fields = [x for x in data_column_field_mapper(data_columns.next())]
    for col_def in data_columns:
        yield dict(zip(fields,col_def))

def read_workbook(wb):
    try:
        logger.info('read screen result file sheets...')
        sheets = workbook_sheets(wb)
        sheet = sheets.next()
        logger.info('read Screen Info sheet: %r', sheet.name)
        meta_cols = sheet_cols(sheet)
        meta_raw = dict(zip(meta_cols.next(),meta_cols.next()))
        meta = { META_MAP.get(key, key):val for key,val in meta_raw.items() }
        
        logger.info('meta: %r', meta)
        
        sheet = sheets.next()
        logger.info('read Data Columns sheet: %r...', sheet.name)
        columns = parse_columns(sheet)
        
        result_values = parse_result_values(columns,sheets)
        return OrderedDict((
            ('meta', meta),
            ('fields', columns),
            ('objects', result_values),
            ))
    except ValidationError,e:
        logger.exception('parse error: %r', e)
        raise 
   
    
def write_workbook(file, screen_facility_id, fields, result_values):
    generic_xls_write_workbook(
        file, 
        create_output_data(screen_facility_id, fields, result_values ))
    
def create_output_data(screen_facility_id, fields, result_values ):
    '''
    Translate Screen Result data into a data structure ready for Serialization:
    {
       'Screen Info': [ [ row1 ], [ row2 ]...].
       'Data Columns': [ [ row1 ], [ row2 ]...].
       'Data': [ [ row1 ], [ row2 ]...].
    }
    @param fields an iterable containing result_value data_column dicts and 
        field information dicts for the non-result value columns
    @param result_values an iterable containing result_value dicts
    '''
    logger.info('create screen result output data %r', screen_facility_id)

    data = OrderedDict()
    
    # 1. Meta sheet
    data['Screen Info'] = { 'Screen Number': screen_facility_id }
    
    # 2. Data Columns
    
    data_column_structure = []
    data_column_keys = []
    non_data_column_keys = []
    keys = sorted(fields.keys())
    for i,key in enumerate(keys):
        field = fields[key]
        if 'ordinal' not in field:
            field['ordinal'] = i
        if ( field.get('is_datacolumn',False) 
            or field.get('data_worksheet_column') is not None):
            data_column_keys.append(key)
        elif ( key not in RESULT_VALUE_FIELD_MAP.keys()
            and key not in REPORTING_NON_RV_COLUMNS ):
            non_data_column_keys.append(key)
    data_column_keys = sorted(
        data_column_keys, key=lambda x: fields[x]['ordinal'])
    non_data_column_keys = sorted(
        non_data_column_keys, key=lambda x: fields[x]['ordinal'])
    data_column_names_to_col_letter = {}
    for i, key in enumerate(data_column_keys):
        data_column_names_to_col_letter[fields[key]['name']] = \
            xl_col_to_name(len(RESULT_VALUE_FIELD_MAP)+i) 
    logger.info('data columns: %r, non_data_column_keys: %r', 
        data_column_keys, non_data_column_keys)
    logger.info('data_column_names_to_col_letter: %r', 
        data_column_names_to_col_letter)

    # Transpose/Pivot the field definitions into the output data_column sheet:
    # Row 0 - "Data" Worksheet Column
    # Row 1 - name
    # Row 2 - data_type
    # Row N - non data column fields
    # Column 0 - data column field label
    # Column 1-N data column values
    datacolumn_labels = DATA_COLUMN_FIELD_MAP.keys()
    header_row = [datacolumn_labels[0]]
    header_row.extend([
        xl_col_to_name(len(RESULT_VALUE_FIELD_MAP)+i) 
            for i in range(len(data_column_keys))])
    logger.info('Data Columns - header_row: %r', header_row)
    for i,(output_label,field_key) in enumerate(
            DATA_COLUMN_FIELD_MAP.items()[1:]):
        row = [output_label]
        for j,key in enumerate(data_column_keys):
            field = fields[key]
            val = field.get(field_key)
            if field_key == 'data_type':
                # TODO: 20170731: migrate the screenresult datacolumn to use 
                # "vocabulary_scope_ref" for the "positive" column types
                # This is a hack to preserve symmetry for read/write for now
                newval = None
                if val == 'string':
                    vocab_scope_ref = field.get('vocabulary_scope_ref')
                    if vocab_scope_ref == 'resultvalue.partitioned_positive':
                        newval = DATA_TYPE.PARTITIONED_POSITIVE
                    elif vocab_scope_ref \
                        == 'resultvalue.confirmed_positive_indicator':
                        newval = DATA_TYPE.CONFIRMED_POSITIVE
                elif val == 'boolean':
                    newval = DATA_TYPE.BOOLEAN_POSITIVE
                if newval:
                    logger.info(
                        'converted: %r:%r to %r: %r', key, field_key, val, newval)
                    val = newval
            if val:
                if field_key == 'is_follow_up_data':
                    if val == True:
                        val = 'Follow up'
                    elif val == False:
                        val = 'Primary'
                elif field_key == 'derived_from_columns':
                    logger.info('derived_from_columns: %r', val)
                    if field.get('screen_facility_id') == screen_facility_id:
                        logger.info('Translate derived_from_columns: %r', val)
                        if not set(data_column_names_to_col_letter.keys())\
                                .issuperset(set(val)):
                            raise ValidationError(
                                key='derived_from_columns', 
                                msg=('col: %r, values: %r are not in %r'
                                    %(key,val,
                                        data_column_names_to_col_letter.keys())))
                        val = ', '.join([
                            data_column_names_to_col_letter[dc_name] 
                                for dc_name in val])
                        logger.info('Translated derived_from_columns: %r', val)
                    else:
                        # Derived column for another screen
                        val = ', '.join(val)
                row.append(val)
            else:
                row.append(None)
                logger.debug(
                    'Note: datacolumn schema key %r is null in result value: %r', 
                    field_key, field)
        logger.debug('data column row: %r', row)
        data_column_structure.append(OrderedDict(zip(header_row,row)))

    data['Data Columns'] = data_column_structure

    # 3. Result Values sheet 
            
    def result_value_generator(result_values):
        
        logger.info('Write the result values sheet')
        header_row = []
        header_row.extend(RESULT_VALUE_FIELD_MAP.keys())
        header_row.extend([
            fields[key].get('title', key) for key in data_column_keys])
        header_row.extend(non_data_column_keys)
        logger.info('Result Values Header row: %r', header_row)

        control_type_mapping = {v:k for k,v in ASSAY_WELL_CONTROL_TYPES.items()}
        row_count = 0
        for result_value in result_values:
            row_count += 1
            if DEBUG_IMPORTER:
                logger.info('result_value: %d: %r', row_count, result_value)
            row = []
            
            row.extend(result_value['well_id'].split(':'))
            if ( result_value.has_key('assay_well_control_type')
                 and result_value['assay_well_control_type'] ):
                control_type = default_converter(
                    result_value['assay_well_control_type'])
                # note: "empty", "experimental", "buffer" are values that can be
                # found in this column, due to legacy data entry, but they are 
                # not valid
                if control_type in control_type_mapping:
                    row.append(control_type_mapping[control_type])
                else:
                    row.append(None)
            else:
                row.append(None)
            excluded_cols = []
            if result_value.has_key('exclude') and result_value['exclude']:
                temp = result_value['exclude']
                if hasattr(temp, 'split'):
                    temp = temp.split(LIST_DELIMITER_SQL_ARRAY)
                logger.debug('excluded data_column_keys: find %r, in %r', 
                    temp, data_column_keys)    
                for data_column_name in temp:
                    # excluded_cols.append(get_column_letter(
                    #     len(RESULT_VALUE_FIELD_MAP)+1
                    #         +data_column_keys.index(data_column_name)))
                    excluded_cols.append(xl_col_to_name(
                        len(RESULT_VALUE_FIELD_MAP)
                            + data_column_keys.index(data_column_name)))
                    excluded_cols = sorted(excluded_cols)
            row.append(','.join(excluded_cols))
            
            if DEBUG_IMPORTER:
                logger.info('write rvs: data_column_keys: %r', data_column_keys)
            for j,key in enumerate(data_column_keys):
                if result_value.has_key(key):
                    row.append(result_value[key])
                else:
                    row.append(None)
            # append the non-result value columns to the end of the row
            for j,key in enumerate(non_data_column_keys):
                if result_value.has_key(key):
                    row.append(result_value[key])
            
            if row_count % 10000 == 0:
                logger.info('generated %d rows', row_count)
            if DEBUG_IMPORTER:
                logger.info('generate row %d: %r',row_count, row)
            yield OrderedDict(zip(header_row,row))
    
    data['Data'] = result_value_generator(result_values)

    return data
        
def json_printer(data):
    
    meta = data['meta']
    columns = meta['fields']
    result_values = data['objects']
    
    return ('{ "screen": %s, "columns": %s, "results": %s }'
        % ( json.dumps(meta['screen_facility_id']), 
            json.dumps(columns),
            json.dumps([result_value for result_value in result_values])
        ))

parser = argparse.ArgumentParser(
    description='reads a serialized (xlsx) screen result')

parser.add_argument(
    '-f', '--screen_result_file', required=True,
    help='''Screen Result file [meta, datacolumns, results].
            examples can be found in reports/static/test_data/screens''')
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

    with open(args.screen_result_file) as input_file:
        wb = xlrd.open_workbook(file_contents=input_file.read())
        print json_printer(read_workbook(wb))
        
        