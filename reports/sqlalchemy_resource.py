from __future__ import unicode_literals

from collections import OrderedDict
from functools import wraps
import hashlib
import logging
import re

from aldjemy.core import get_engine, get_tables
from django.conf import settings
from django.core.cache import cache
import django.db.models.constants
from django.http.request import HttpRequest
from django.http.response import StreamingHttpResponse, HttpResponse, Http404
from sqlalchemy import select, asc, text
import sqlalchemy
from sqlalchemy.dialects import postgresql
from sqlalchemy.dialects.postgresql import array
from sqlalchemy.sql import and_, or_, not_, func
from sqlalchemy.sql import asc, desc, alias, Alias
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join, cast
from sqlalchemy.sql.expression import nullsfirst, nullslast
from sqlalchemy.sql.functions import func
import sqlalchemy.sql.sqltypes
from tastypie.exceptions import BadRequest, ImmediateHttpResponse

from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    LIST_BRACKETS, MAX_IMAGE_ROWS_PER_XLS_FILE, MAX_ROWS_PER_XLS_FILE, \
    HTTP_PARAM_RAW_LISTS,HTTP_PARAM_DATA_INTERCHANGE, HTTP_PARAM_USE_TITLES,\
    HTTP_PARAM_USE_VOCAB
from reports.api_base import IccblBaseResource, un_cache
from reports.serialize import XLSX_MIMETYPE, SDF_MIMETYPE, XLS_MIMETYPE,\
    JSON_MIMETYPE, CSV_MIMETYPE, parse_val
from reports.serialize.streaming_serializers import sdf_generator, \
    json_generator, get_xls_response, csv_generator, ChunkIterWrapper, \
    cursor_generator, image_generator, closing_iterator_wrapper
from reports.serializers import LimsSerializer
import json
import urllib
import six
from django.test.client import RequestFactory


logger = logging.getLogger(__name__)

DEBUG_FILTERS = False or logger.isEnabledFor(logging.DEBUG)

def _concat(*args):
    '''
    Use as a replacement for sqlalchemy.sql.functions.concat
    - "concat" is not available in postgresql 8.4
    '''
    return func.array_to_string(array([x for x in args]),'')
        
def _concat_with_sep(args=None,sep=None):
    '''
    Use as a replacement for sqlalchemy.sql.functions.concat
    - "concat" is not available in postgresql 8.4
    '''
    new_args = []
    for arg in args:
        new_args.append(arg)
        new_args.append(sep)
    new_args = new_args[:-1]
    return _concat(*new_args)

class SqlAlchemyResource(IccblBaseResource):
    '''
    A resource that uses SqlAlchemy to facilitate the read queries:
    - get_list
    - get_detail
    Note: write operations not implemented
    '''
    
    
    def __init__(self, *args, **kwargs):
        # store the Aldjemy tables in a "bridge" object, for legacy reasons
        self.bridge = get_tables()
        
        self.use_cache = True
        self.request_factory = RequestFactory()
        
        super(SqlAlchemyResource, self).__init__(*args, **kwargs)
    
    
    @classmethod
    def wrap_statement(cls, stmt, order_clauses, filter_expression):
        '''
        @param stmt - a sqlalchemy.sql.expression.join instance
        @param order_clauses - list of sqlalchemy.sql.expression.column
        @param filter_expression - sqlalchemy.whereclause
        '''
        if order_clauses:
            _alias = Alias(stmt)
            stmt = select([text('*')]).select_from(_alias)
            stmt = stmt.order_by(*order_clauses)
        if filter_expression is not None:
            logger.debug('filter_expression: %r' % filter_expression)
            if not order_clauses:
                _alias = Alias(stmt)
                stmt = select([text('*')]).select_from(_alias)
            stmt = stmt.where(filter_expression)

        count_stmt = select([func.count()]).select_from(stmt.alias())
        return (stmt,count_stmt)
    
    def _convert_request_to_dict(self, request):
        '''
        Transfer all values from GET, then POST to a dict
        Note: uses 'getlist' to retrieve all values:
        - if a value is single valued, then unwrap from the list
        - downstream methods expecting a list value must deal with non-list single values
        '''
        DEBUG = False or logger.isEnabledFor(logging.DEBUG)
        
        _dict = {}
        for key in request.GET.keys():
            val = request.GET.getlist(key)
            if DEBUG:
                logger.info('get key: %r, val: %r', key, val)
            # Jquery Ajax will send array list params with a "[]" suffix - 20151015
            if '[]' in key and key[-2:] == '[]':
                key = key[:-2]
            if len(val) == 1:
                _dict[key] = val[0]
            else:
                _dict[key] = val
            
        for key in request.POST.keys():
            val = request.POST.getlist(key)
            if DEBUG:
                logger.info('post key: %r, val: %r', key, val)
            # Jquery Ajax will post array list params with a "[]" suffix - 20151015
            key = key.replace('[]','')
            if len(val) == 1:
                _dict[key] = val[0]
            else:
                _dict[key] = val
        
        # check for single-valued known list values
        # Note: Jquery Ajax will post array list params with a "[]" suffix - 20151015
        known_list_values = [
            'includes','exact_fields', 'order_by', 'visibilities','other_screens',  
            'includes[]', 'order_by[]','exact_fields[]', 'visibilities[]',
            'other_screens[]']
        for key in known_list_values:
            val = _dict.get(key,[])
            if isinstance(val, basestring):
                _dict[key] = [val]
        
        # Parse known boolean params for convenience
        http_boolean_params = [
            HTTP_PARAM_DATA_INTERCHANGE,HTTP_PARAM_RAW_LISTS,
            HTTP_PARAM_USE_VOCAB, HTTP_PARAM_USE_TITLES]
        for key in http_boolean_params:
            _dict[key] = parse_val(
                _dict.get(key, False),key, 'boolean')
        if DEBUG:
            logger.info('params: %r', _dict)
        return _dict    
    
    def get_visible_fields(self, 
        schema_fields, filter_fields, manual_field_includes,
        visibilities, exact_fields=[], order_params=[]):
        '''
        Construct an ordered dict of schema fields that are visible, based on
        - the field["visibility"] of each field on the resource,
        - if the field is in the manual_field_includes
        - if the field is in the filter_fields
        - if the field key in another fields schema field['dependencies'] 
        
        
        TODO: this method is not SqlAlchemy specific
        '''
        DEBUG_VISIBILITY = False or logger.isEnabledFor(logging.DEBUG)
        visibilities = set(visibilities)
        if DEBUG_VISIBILITY:
            logger.info('get_visible_fields: field_hash initial: %r, manual: %r, exact: %r', 
                schema_fields.keys(),manual_field_includes, exact_fields )
            logger.info('field visibility settings: %r', 
                [ str((key,field['visibility'])) for key,field in schema_fields.items()])
        try:
            if exact_fields:
                temp = { key:field for key,field in schema_fields.items()
                    if key in exact_fields or key in filter_fields }
            else:
                temp = { key:field for key,field in schema_fields.items() 
                    if ((field.get('visibility', None) 
                            and visibilities & set(field['visibility'])) 
                        or field['key'] in manual_field_includes
                        or '*' in manual_field_includes ) }
            
                # manual excludes
                temp = { key:field for key,field in temp.items() 
                    if '-%s' % key not in manual_field_includes }
            
            # dependency fields
            dependency_fields = set()
            for field in temp.values():
        
                if field.get('value_template', None):
                    dependency_fields.update(
                        re.findall(r'{([a-zA-Z0-9_-]+)}', field['value_template']))
                if field.get('display_options', None):
                    dependency_fields.update(
                        re.findall(r'{([a-zA-Z0-9_-]+)}', field['display_options']))
                if field.get('dependencies',None):
                    dependency_fields.update(field.get('dependencies'))
                logger.debug('field: %s, dependencies: %s', field['key'],field.get('dependencies',[]))
            if DEBUG_VISIBILITY:
                logger.info('dependency_fields %s', dependency_fields)
            if dependency_fields:
                temp.update({ key:field 
                    for key,field in schema_fields.items() if key in dependency_fields })
            
            # filter_fields
            if filter_fields:
                temp.update({ key:field 
                    for key,field in schema_fields.items() if key in filter_fields })
            # order params
            if order_params:
                temp.update({ key:field 
                    for key,field in schema_fields.items() 
                        if ( key in order_params or '-%s'%key in order_params) })
            
            field_hash = OrderedDict(sorted(temp.iteritems(), 
                key=lambda x: x[1].get('ordinal',999))) 
    
            if DEBUG_VISIBILITY:
                logger.info('field_hash final: %s', field_hash.keys())
        
            if not field_hash:
                logger.info('no fields found: %r, %r', field_hash.keys(), manual_field_includes)
                response = HttpResponse('no fields specified')
                response.status_code = 400
                raise ImmediateHttpResponse(
                    response=response)
            
            return field_hash
        
        except ImmediateHttpResponse:
            raise
        
        except Exception, e:
            logger.exception('on get_visible_fields')
            raise e 

    def build_sqlalchemy_columns(
            self, fields, base_query_tables=None, custom_columns=None):
        '''
        Returns an ordered dict of sqlalchemy.sql.schema.Column objects, associated 
        with the sqlalchemy.sql.schema.Table definitions, which are bound to 
        the sqlalchemy.engine.Engine: 
        "Connects a Pool and Dialect together to provide a source of database 
        connectivity and behavior."
        
        @param fields - field definitions, from the resource schema
        @param bridge - a reports.utils.sqlalchemy_bridge.Bridge
        @param base_query_tables - if specified, the fields for these tables 
        will be available as part of the base query, so the column definitions
        become simpler, and do not need to be joined in. 
        @param manual_includes - columns to include even if the field 
        visibility is not set
        '''
        DEBUG_BUILD_COLUMNS = False or logger.isEnabledFor(logging.DEBUG)
        base_query_tables = base_query_tables or []
        custom_columns = custom_columns or []
        
        try:
            columns = OrderedDict()
            for field in fields:
                key = field['key']
                if key in custom_columns:
                    if DEBUG_BUILD_COLUMNS: 
                        logger.info(
                            'custom field: %r, %r', key,custom_columns[key])
                    columns[key] = custom_columns[key].label(key)
                    continue
                
                if DEBUG_BUILD_COLUMNS: 
                    logger.info('build column: %r, %r', field['key'], field)
                field_name = field.get('field', None)
                if not field_name:
                    field_name = field['key']
                
                field_table = field.get('table', None)
                
                if not field_table and DEBUG_BUILD_COLUMNS:
                    logger.info(
                        'field: %r, val: %r, skipping field because there is no '
                        '"field_table" value set',key,field)
                    continue
                if DEBUG_BUILD_COLUMNS: 
                    logger.info(
                        'field: %r, field_table: %r', field['key'], field_table )
                
                if field_table in base_query_tables:
                    # simple case: table.fields already selected in the base query:
                    # just need to specify them
                    if field_name in get_tables()[field_table].c:
                        col = get_tables()[field_table].c[field_name]
                    else:
                        raise Exception(
                            'field: %r, not found in table: %r'
                            % (field_name, field_table))
                    col = col.label(key)
                    columns[key] = col
                    
                # TODO: remove this; favor custom linking to subtables
                # used in reagent subclasses
                elif field.get('linked_field_value_field', None):
                    link_table = field['table']
                    link_table_def = get_tables()[link_table]
                    linked_field_parent = field['linked_field_parent']
                    link_field = linked_field_parent + '_id'
                    join_args = { 
                        'link_table': link_table, 'link_field': link_field,
                        'parent_table': linked_field_parent
                        }
                    
                    if field['linked_field_type'] != 'fields.ListField':
                        join_stmt = select([link_table_def.c[field_name]]).\
                            where(text('{link_table}.{link_field}='
                                    '{parent_table}.{link_field}'.format(**join_args)))
                        if field.get('linked_field_meta_field', None):
                            # TODO: test - the linked meta field is the "datacolumn type"
                            linked_field_meta_field = field['linked_field_meta_field']
                            meta_field_obj = MetaHash.objects.get(
                                key=field['key'], scope=field['scope'])
                            meta_table_def = get_tables()['metahash']
                            join_stmt.join(meta_table_def, 
                                link_table_def.c[linked_field_meta_field]==
                                    getattr(meta_field_obj,'pk') )
                        join_stmt = join_stmt.label(key)
                        columns[key] = join_stmt
                    elif field['linked_field_type'] == 'fields.ListField':
                        join_stmt = select([link_table_def.c[field_name]]).\
                            where(text('{link_table}.{link_field}='
                                    '{parent_table}.{link_field}'.format(**join_args)))
    
                        if field.get('linked_field_meta_field', None):
                            # TODO: test - the linked meta field is the "datacolumn type"
                            linked_field_meta_field = field['linked_field_meta_field']
                            meta_field_obj = MetaHash.objects.get(
                                key=field['key'], scope=field['scope'])
                            meta_table_def = get_tables()['metahash']
                            join_stmt.join(meta_table_def, 
                                link_table_def.c[linked_field_meta_field]==
                                    getattr(meta_field_obj,'pk') )
                        
                        ordinal_field = field.get('ordinal_field', None)
                        if ordinal_field:
                            join_stmt = join_stmt.order_by(link_table_def.c[ordinal_field])
                        join_stmt = join_stmt.alias('a')
                        stmt2 = select([func.array_to_string(
                                        func.array_agg(column(field_name)),
                                                       LIST_DELIMITER_SQL_ARRAY)])
                        stmt2 = stmt2.select_from(join_stmt).label(key)
                        columns[key] = stmt2
                else:
                    if DEBUG_BUILD_COLUMNS:        
                        logger.info(
                            'field is not in the base tables %r, nor in a linked field, '
                            'and is not custom: %s', base_query_tables, key)
            if DEBUG_BUILD_COLUMNS: 
                logger.info('columns: %r', columns.keys())
            return columns
        except Exception, e:
            logger.exception('on build sqlalchemy columns')
            raise e   

    @staticmethod
    def build_sqlalchemy_ordering(order_params, visible_fields):
        '''
        returns a scalar or list of ClauseElement objects which will comprise 
        the ORDER BY clause of the resulting select.
        @param order_params passed as list in the request.GET hash
        '''
        DEBUG_ORDERING = False or logger.isEnabledFor(logging.DEBUG)
        
        if DEBUG_ORDERING:
            logger.info('build sqlalchemy ordering: %s, visible fields: %s',
                order_params,visible_fields.keys())
        if order_params and isinstance(order_params, basestring):
            # standard, convert single valued list params
            order_params = [order_params]
        order_clauses = []
        for order_by in order_params:
            field_name = order_by
            order_clause = None
            if order_by.startswith('-'):
                field_name = order_by[1:]
                order_clause = nullslast(desc(column(field_name)))
                if ( field_name in visible_fields 
                    and visible_fields[field_name]['data_type'] == 'string'):
                    # For string field ordering, double sort as numeric and text
                    order_clause = text(
                        "(substring({field_name}, '^[0-9]+'))::int desc " # cast to integer
                        ",substring({field_name}, '[^0-9_].*$')  desc"  # works as text
                        .format(field_name=field_name))
            else:
                order_clause = nullsfirst(asc(column(field_name)))
                if ( field_name in visible_fields 
                    and visible_fields[field_name]['data_type'] == 'string'):
                    order_clause = text(
                        "(substring({field_name}, '^[0-9]+'))::int "
                        ",substring({field_name}, '[^0-9_].*$') "
                        .format(field_name=field_name))
            if field_name in visible_fields:
                order_clauses.append(order_clause)
            else:
                logger.warn(
                    'order_by field %r not in visible fields, skipping: ', 
                    order_by)
        if DEBUG_ORDERING:
            logger.info('order_clauses %s',order_clauses)     
        return order_clauses
    
    @staticmethod
    def filter_value_to_python(value, filter_type):
        """
        Turn the string ``value`` into a python object.
        """
        original_value = value
        if isinstance(value, six.string_types):
            value = urllib.unquote(value).decode('utf-8')
        # Simple values
        if value is 1:
            value = 1
        elif value in ['true', 'True', True]:
            value = True
        elif value in ['false', 'False', False]:
            value = False
        elif value in ('nil', 'none', 'None', None):
            value = None

        if filter_type in ('in', 'range') and len(value):
            if value and hasattr(value, '__iter__'):
                # value is already a list
                pass
            # Removed 20161108
            # elif hasattr(param_hash, 'getlist'):
            #     value = []
            #     for part in param_hash.getlist(filter_expr):
            #         value.extend(part.split(LIST_DELIMITER_URL_PARAM))
            else:
                value = value.split(LIST_DELIMITER_URL_PARAM)
        if DEBUG_FILTERS:
            logger.info('filter_value_to_python: %r: %r, %r', 
                original_value, value, filter_type)
        return value

    @staticmethod
    def build_sqlalchemy_filters(schema, param_hash):
        '''
        Build the full SqlAlchemy filter expression for all filters and search:
        
        @param param_hash: a hash of filter data:
        - filters defined as filter_key, filter_value; combined using "AND"
        - nested "nested_search_data":
            - an array of hashes of filter_data, to be OR'd together
            - each hash consists of filter_key, filter_value
        @return (
            search_expression - full combined SqlAlchemy search expression
            combined_filter_hash - field_keyed hash of expressions
            )
        '''
        
        if DEBUG_FILTERS: 
            logger.info('build_sqlalchemy_filters: param_hash %s', param_hash)

        # ordinary filters
        (filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filter_hash(schema, param_hash)
        combined_filter_hash = filter_hash
        filter_expression = and_(*filter_hash.values())

        # 20170511 - nested search_data not used (for well, plate, or screening inquiry)
        # Treat the nested "nested_search_data" as sets of params to be OR'd together,
        # then AND'd with the regular filters (if any)
        nested_search_data = param_hash.get('nested_search_data', None)
        if nested_search_data:
            logger.info('nested_search_data: %r', nested_search_data)
            if isinstance(nested_search_data, basestring):
                nested_search_data = json.loads(nested_search_data)
            if isinstance(nested_search_data, dict):
                # a standard dict of filters, to be OR'd with the current filter expression
                nested_search_data = [nested_search_data]
            if not isinstance(nested_search_data, (list,tuple)):
                raise Exception('nested_search_data must be a list of dicts')
            # each item in the search data array is a search hash, to be or'd
            search_expressions = []
            filter_fields = set(filter_hash.keys())
            for search_hash in nested_search_data:
                logger.info('search_hash: %s' % search_hash)
                
                (search_filter_hash,readable_search_filter_hash) = \
                    SqlAlchemyResource.\
                        build_sqlalchemy_filter_hash(schema,search_hash)
                search_expressions.append(and_(*search_filter_hash.values()))
                
                # Append search expressions for each field to a combined hash
                for field,expression in search_filter_hash.items():
                    logger.info(
                        'search filter to combine: %r, %r', field, expression)
                    if field in combined_filter_hash:
                        combined_filter_hash[field] = \
                            or_(combined_filter_hash[field], expression)
                    else:
                        combined_filter_hash[field] = expression
            if len(search_expressions) > 1:
                search_expressions = or_(*search_expressions)
            else:
                search_expressions = search_expressions[0]
            if len(filter_hash) > 0:
                filter_expression = and_(
                    search_expressions,
                    filter_expression)
            else: 
                filter_expression = search_expressions
            readable_filter_hash['search'] = '_'.join(search_hash.keys())    
        if DEBUG_FILTERS: 
            logger.info('filter_expression: %s, filter_fields: %s',
                filter_expression, combined_filter_hash.keys())
            logger.info(
                'readable_filter_hash: %r', readable_filter_hash)
        
        return (filter_expression,combined_filter_hash, readable_filter_hash)
    
    @staticmethod
    def parse_filter(filter_expr, value):
        if DEBUG_FILTERS:
            logger.info('parse filter: %r, %r', filter_expr, value)
        
        lookup_sep = django.db.models.constants.LOOKUP_SEP
        if lookup_sep not in filter_expr:
            # treat as a field__eq
            field_name = filter_expr
            filter_type = 'exact'
            if DEBUG_FILTERS:
                logger.info('interpret: %r as %r for %r:%r', 
                    filter_expr, filter_type,field_name, value)
        else:
            filter_bits = filter_expr.split(lookup_sep)
            if len(filter_bits) != 2:
                logger.warn(
                    'filter expression %r must be of the form '
                    '"field_name__expression"' % filter_expr )
            field_name = filter_bits[0]
            filter_type = filter_bits[1]
        inverted = False
        if field_name and field_name[0] == '-':
            inverted = True
            field_name = field_name[1:]
        if DEBUG_FILTERS:
            logger.info('build filter expr: field_name, %r, '
                'filter_type: %r value: %r', 
                field_name, filter_type, value)

        # 20161108: "param_hash" is sent to filter_value_to_python for obsoleted reason:
        # - it was the GET.param hash, and was checking if value was a list...
        #         value = SqlAlchemyResource.filter_value_to_python(
        #             value, param_hash, filter_expr, filter_type)
        value = SqlAlchemyResource.filter_value_to_python(
            value, filter_type)
            
        return (field_name, value, filter_type, inverted)
    
    @staticmethod
    def build_filter( 
        field_name, data_type, filter_type, inverted, value ):

        if DEBUG_FILTERS:
            logger.info('build filter: %r, %r, %r, %r, %r', 
                field_name, data_type, filter_type, inverted, value)
        
        expression = None
        col = column(field_name)
        if data_type in ['integer', 'float', 'decimal']:
            col = cast(col, sqlalchemy.sql.sqltypes.Numeric)
        elif data_type == 'boolean':
            if filter_type != 'is_null':
                col = cast(func.coalesce(col,False), sqlalchemy.sql.sqltypes.Boolean)
        if data_type == 'string':
            col = cast(col, sqlalchemy.sql.sqltypes.Text)
        if filter_type in ['exact', 'eq']:
            if data_type == 'string':
                value = str(value)
            expression = col == value
            if data_type == 'list':
                expression = text(
                    "'%s'=any(string_to_array(%s,'%s'))" 
                        % (value, field_name, LIST_DELIMITER_SQL_ARRAY))
        elif filter_type == 'about':
            decimals = 0
            if '.' in value:
                decimals = len(value.split('.')[1])
            expression = func.round(col, decimals) == value
            if DEBUG_FILTERS:
                logger.info(
                    'create "about" expression for: %r, %r, decimals %r', 
                    field_name, value, decimals)
        elif filter_type == 'contains':
            if data_type == 'string':
                value = str(value)
            expression = col.contains(value)
        elif filter_type == 'icontains':
            if data_type == 'string':
                value = str(value)
            expression = col.ilike('%{value}%'.format(value=value))
        elif filter_type == 'lt':
            expression = col < value
        elif filter_type == 'lte':
            expression = col <= value
        elif filter_type == 'gt':
            expression = col > value
        elif filter_type == 'gte':
            expression = col >= value
        elif filter_type == 'is_blank':
            if data_type == 'string':
                col = func.trim(col)
            if value and str(value).lower() == 'true':
                expression = col == None
                if data_type == 'string':
                    expression = col == ''
            else:
                expression = col != None
                if data_type == 'string':
                    col = func.trim(col)
                if (data_type == 'string' or 
                    data_type == 'list'):
                    expression = col != ''
        elif filter_type == 'is_null':
            if value and str(value).lower() == 'true':
                expression = col == None
            else:
                expression = col != None
            # TODO: test that col <> '' expression is created
        elif filter_type == 'in':
            if data_type == 'list': # NOTE: for the list type, interpret "in" as any of the
                # given values are in the field
                temp_expressions = []
                for _val in value:
                    temp_expressions.append(col.ilike('%{value}%'.format(value=_val)))
                
                expression = or_(*temp_expressions)
            else:
                expression = col.in_(value)
        elif filter_type == 'ne':
            if data_type == 'string':
                value = str(value)
            expression = col != value
        elif filter_type == 'range':
            if len(value) != 2:
                logger.error('field: %r, val: %r, '
                    'range expression must be list of length 2', 
                    field_name, value)
            else:
                expression = col.between(value[0], value[1], symmetric=True)
        else:
            logger.error(
                'field: %r, unknown filter type: %r for value: %r', 
                field_name, filter_type, value)
        if inverted:
            expression = not_(expression)
        return expression

    @staticmethod
    def build_sqlalchemy_filter_hash(schema, param_hash):
        '''
        Create a SqlAlchemy whereclause out of django style filters:
        - field_name__filter_expression
        
        @return field_keyed hash of expressions for the AND clause
        '''
        logger.debug('build_sqlalchemy_filter_hash %r' % param_hash)

        if param_hash is None:
            return (None,None)
        
        if not isinstance(param_hash, dict):
            raise Exception('filter hash must be a dict: %r', param_hash)
        
        filter_hash = {}
        readable_filter_hash = {}
        try:
            for filter_expr, value in param_hash.items():
                
                (field_name, value, filter_type, inverted) = \
                    SqlAlchemyResource.parse_filter(filter_expr, value)
                
                if not field_name in schema['fields']:
                    logger.debug('unknown filter field: %r, %r', field_name, filter_expr)
                    continue
                
                field = schema['fields'][field_name]
                expression = SqlAlchemyResource.build_filter(
                    field_name, field['data_type'], filter_type, inverted, value)
                if expression is not None:
                    filter_hash[field_name] = expression
                    
                    readable_expression = []
                    if field_name in filter_expr:
                        if inverted is True:
                            readable_expression.append('not')
                        if filter_type not in ('eq','exact'):
                            readable_expression.append(filter_type)
                    readable_value = value
                    if isinstance(readable_value,(list,tuple)):
                        readable_value = '_'.join([str(x) for x in readable_value]) 
                    readable_expression.append(str(readable_value))
                    readable_filter_hash[field_name] = '_'.join(readable_expression)
            if DEBUG_FILTERS:
                logger.info('filtered_fields: %s', filter_hash.keys())
                logger.info('readable_filter_hash: %r', readable_filter_hash)
            return (filter_hash, readable_filter_hash)
        except Exception, e:
            logger.exception('on build_sqlalchemy_filter_hash')
            raise e   

    def _get_list_response(self,request,**kwargs):
        '''
        Return a deserialized list of dicts
        '''
        logger.debug('_get_list_response: %r, %r', 
            self._meta.resource_name, {k:v for k,v in kwargs.items() if k !='schema'})
        includes = kwargs.pop('includes', '*')
        try:
            kwargs.setdefault('limit', 0)
            response = self.get_list(
                request,
                format='json',
                includes=includes,
                **kwargs)
            _data = self.get_serializer().deserialize(
                LimsSerializer.get_content(response), JSON_MIMETYPE)
            if self._meta.collection_name in _data:
                _data = _data[self._meta.collection_name]
            logger.debug(' data: %r', _data)
            return _data
        except Http404:
            return []
        except Exception as e:
            logger.exception('on get list: %r', e)
            raise
        
    def _get_detail_response(self,request,**kwargs):
        '''
        Return the detail response as a dict
        '''
        logger.info('_get_detail_response: %r, %r', 
            self._meta.resource_name, {k:v for k,v in kwargs.items() if k !='schema'})
        includes = kwargs.pop('includes', '*')
        try:
            response = self.get_detail(
                request,
                format='json',
                includes=includes,
                **kwargs)
            _data = {}
            if response.status_code == 200:
                _data = self._meta.serializer.deserialize(
                    LimsSerializer.get_content(response), JSON_MIMETYPE)
            else:
                logger.info(
                    'no data found for %r, %r, %r', 
                    self._meta.resource_name, kwargs, response.status_code)
            return _data
        except Http404:
            return {}
        
    #@un_cache
    def _get_detail_response_internal(self, **kwargs):
        request = self.request_factory.generic('GET', '.', 
            data={ 'CONTENT_TYPE': JSON_MIMETYPE, 'HTTP_ACCEPT': JSON_MIMETYPE },
            content_type=JSON_MIMETYPE)
        class User:
            is_superuser = True
            # @staticmethod
            # def is_superuser():
            #     return true
        request.user = User
        result = self._get_detail_response(request, **kwargs)
        return result

    # @un_cache
    def _get_list_response_internal(self, **kwargs):
        
        request = self.request_factory.generic('GET', '.', 
            data={ 'CONTENT_TYPE': JSON_MIMETYPE, 'HTTP_ACCEPT': JSON_MIMETYPE },
            content_type=JSON_MIMETYPE)
        content_type = self.get_serializer().get_accept_content_type(request)
        logger.info('_get_list_response_internal: %r', content_type)
        class User:
            is_superuser = True
        request.user = User
        result = self._get_list_response(request, **kwargs)
        return result
    
    def _cached_resultproxy(self, conn, stmt, count_stmt, param_hash, limit, offset):
        ''' 
        Cache for resultsets:
        - Always returns the cache object with a resultset, either from the cache,
        or executed herein.
        
        NOTE: limit and offset are included because this version of sqlalchemy
        does not support printing of them with the select.compile() function.
        
        TODO: cache clearing on database writes
        '''
        DEBUG_CACHE = False or logger.isEnabledFor(logging.DEBUG)
        # Limit check removed with the use of "use_caching" flag
        # if limit == 0:
        #    raise Exception('limit for caching must be >0')
        
        prefetch_number = 5
        if limit <= 1:
            prefetch_number = 1
        
        try:
            # use a hexdigest because statements can contain problematic chars 
            # for the memcache
            m = hashlib.md5()
            compiled_stmt = str(stmt.compile(
                dialect=postgresql.dialect(), compile_kwargs={"literal_binds": True}))
            logger.debug('compiled_stmt %s',compiled_stmt)
            if 'limit' in compiled_stmt.lower():
                # remove limit and offset; will calculate
                compiled_stmt = compiled_stmt[:compiled_stmt.lower().rfind('limit')]
            logger.debug('compiled_stmt for hash key: %s', compiled_stmt)
            key_digest = '%s_%s_%s' %(compiled_stmt, str(limit), str(offset))
            logger.debug('key_digest: %r', key_digest)
            m.update(key_digest)
            key = m.hexdigest()
            logger.debug('hash key: digest: %s, key: %s, limit: %s, offset: %s', 
                key_digest, key, limit, offset)
            cache_hit = cache.get(key)
            if cache_hit is not None:
                if ('stmt' not in cache_hit or
                        cache_hit['stmt'] != compiled_stmt):
                    cache_hit = None
                    logger.warn('cache collision for key: %r, %r', key, stmt)
            
            if cache_hit is None:
                logger.info('no cache hit for key: %r, executing stmt', key)
                # Note: if no cache hit, then retrive limit*n results, and 
                # cache several iterations at once.
                new_limit = limit * prefetch_number 
                if DEBUG_CACHE:
                    logger.info(
                        'no cache hit, create cache, limit: %s, '
                        'new limit for caching: %s',
                        limit, new_limit)
                if new_limit > 0:
                    stmt = stmt.limit(new_limit)
                resultset = conn.execute(stmt)
                prefetched_result = [dict(row) for row in resultset] if resultset else []
                logger.info('executed stmt %d', len(prefetched_result))
                
                if len(prefetched_result) < limit & offset == 0:
                    # Optimize, skip count if first page and less than limit are found.
                    count = len(prefetched_result)
                else:
                    logger.info('no cache hit, execute count')
                    if limit == 1:
                        logger.info('set count to 1, detail view')
                        count = 1
                    else:
                        count = conn.execute(count_stmt).scalar()
                logger.info('count: %s', count)
                
#                 if count == 1 or count < limit:
#                     return {
#                         'stmt': compiled_stmt,
#                         'cached_result': prefetched_result,
#                         'count': count,
#                     }
                
                if limit==0 and count > settings.MAX_ROWS_FOR_CACHE_RESULTPROXY:
                    logger.warn('too many rows to cache: %r, limit: %r, '
                        'see setting"MAX_ROWS_FOR_CACHE_RESULTPROXY"',
                        count, settings.MAX_ROWS_FOR_CACHE_RESULTPROXY)
                    return None
                
                # now fill in the cache with the prefetched sets or rows
                for y in range(prefetch_number):
                    new_offset = offset + limit*y;
                    _start = limit*y
                    if DEBUG_CACHE:
                        logger.info('new_offset: %d, start: %d, len: %d', 
                            new_offset, _start, len(prefetched_result))
                    if _start < len(prefetched_result):
                        key_digest = '%s_%s_%s' %(compiled_stmt, str(limit), str(new_offset))
                        m = hashlib.md5()
                        m.update(key_digest)
                        key = m.hexdigest()
                        rows_to_fetch = limit
                        if limit==0:
                            rows_to_fetch = count+1
                        _result = prefetched_result[_start:_start+rows_to_fetch]
                        _cache = {
                            'stmt': compiled_stmt,
                            'cached_result': _result,
                            'count': count,
                            'key': key }
                        if DEBUG_CACHE:
                            logger.info('add to cache, key: %s, limit: %s, offset: %s',
                                key, limit, new_offset)
                        cache.set( key, _cache, None)
                        if y == 0:
                            cache_hit = _cache
                    else:
                        logger.info('not caching: prefetched length: %s, start: %s, end: %s',
                            len(prefetched_result), _start, 'end')
                        break
                logger.info('store cached iterations: %s', y+1)
            else:
                logger.info('cache hit for key: %r', key)   
                
            return cache_hit
 
        except Exception, e:
            logger.exception('on conn execute')
            raise e   

    def stream_response_from_statement(self, request, stmt, count_stmt, 
            output_filename, field_hash={}, param_hash={}, 
            rowproxy_generator=None, is_for_detail=False,
            downloadID=None, title_function=None, use_caching=None, meta=None,
            format=None ):
        '''
        Execute the SQL stmt provided and stream the results to the response:
        
        Caching (for json responses only): resources will be cached if:
        - self.use_caching is True and use_caching is not False and limit > 0
        - limit == 0 and use_caching is True
        
        '''
        DEBUG_STREAMING = False or logger.isEnabledFor(logging.DEBUG)
        
        logger.info('stream_response_from_statement: %r %r', self._meta.resource_name, format )
        temp_param_hash = param_hash.copy()
        if 'schema' in temp_param_hash:
            del temp_param_hash['schema']
        if DEBUG_STREAMING:
            logger.info('stream_response_from_statement: %r, %r', 
                self._meta.resource_name,temp_param_hash)
        limit = param_hash.get('limit', 25)        
        try:
            limit = int(limit)
        except Exception:
            raise BadRequest(
                "Invalid limit '%s' provided. Please provide a positive integer." 
                % limit)
        if limit > 0:    
            stmt = stmt.limit(limit)
        if is_for_detail:
            limit = 1

        offset = param_hash.get('offset', 0 )
        try:
            offset = int(offset)
        except Exception:
            raise BadRequest(
                "Invalid offset '%s' provided. Please provide a positive integer." 
                % offset)
        if offset < 0:    
            offset = -offset
        stmt = stmt.offset(offset)
        
        conn = get_engine().connect()
        
        try:
            logger.debug('offset: %s, limit: %s', offset, limit)
        
            if DEBUG_STREAMING:
                logger.info('stmt: %s, param_hash: %s ', 
                    str(stmt.compile(
                            dialect=postgresql.dialect(), 
                            compile_kwargs={"literal_binds": True})), 
                    temp_param_hash)
                logger.info(
                    'count stmt %s', 
                    str(count_stmt.compile(
                        dialect=postgresql.dialect(), 
                        compile_kwargs={"literal_binds": True})))
           
            logger.info('format: %r', format) 
            if format is not None:
                content_type = self.get_serializer().get_content_type_for_format(format)
            else:
                content_type = self.get_serializer().get_accept_content_type(request)

            result = None
            if content_type == JSON_MIMETYPE:
                logger.info(
                    'streaming json, use_caching: %r, self.use_cache: %r, limit: %d, %r', 
                    use_caching, self.use_cache, limit, is_for_detail)
                if ((self.use_cache is True and use_caching is not False)
                        and ( use_caching is True or limit > 0)):
                    cache_hit = self._cached_resultproxy(
                        conn, stmt, count_stmt, param_hash, limit, offset)
                    if cache_hit:
                        result = cache_hit['cached_result']
                        count = cache_hit['count']
                    else:
                        # cache routine should always return a cache object
                        logger.info('cache not set: execute stmt')
                        count = conn.execute(count_stmt).scalar()
                        result = conn.execute(stmt)
                    logger.info('====count: %d====', count)
                    
                else:
                    logger.info('not cached, execute count stmt...')
                    # compiled_stmt = str(count_stmt.compile(
                    #     dialect=postgresql.dialect(),
                    #     compile_kwargs={"literal_binds": True}))
                    # logger.info('compiled_stmt %s', compiled_stmt)
                    
                    count = conn.execute(count_stmt).scalar()
                    logger.info('excuted count stmt: %d', count)
                    result = conn.execute(stmt)
                    logger.info('excuted stmt')

                if not meta:
                    meta = {
                        'limit': limit,
                        'offset': offset,
                        'total_count': count
                        }
                else:
                    temp = {
                        'limit': limit,
                        'offset': offset,
                        'total_count': count
                        }
                    temp.update(meta)    
                    meta = temp
                
                if rowproxy_generator:
                    result = rowproxy_generator(result)
                    
                logger.info('is for detail: %r, count: %r', is_for_detail, count)
                if is_for_detail and count == 0:
                    logger.info('detail not found')
                    conn.close()
                    return HttpResponse(status=404)
                
                if DEBUG_STREAMING:
                    logger.info('json setup done, meta: %r', meta)
    
            else: # not json
            
                logger.info('excute stmt')
                result = conn.execute(stmt)
                logger.info('excuted stmt')
                
                if rowproxy_generator:
                    result = rowproxy_generator(result)
                    # FIXME: test this for generators other than json generator        
            
            result = closing_iterator_wrapper(result, conn.close)
            return self.stream_response_from_cursor(request, result, output_filename, 
                field_hash=field_hash, 
                param_hash=param_hash, 
                is_for_detail=is_for_detail, 
                downloadID=downloadID, 
                title_function=title_function, 
                meta=meta, format=format)
        except Exception, e:
            logger.exception('on stream response')
            raise e
    
    def stream_response_from_cursor(
            self,request,result,output_filename, field_hash={}, param_hash={}, 
            is_for_detail=False, downloadID=None, title_function=None, 
            meta=None, format=None):
          
        try:

            list_brackets = LIST_BRACKETS
            if ( param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
                or request.GET.get(HTTP_PARAM_RAW_LISTS, False)):
                list_brackets = None
    
            if format is not None:
                content_type = \
                    self.get_serializer().get_content_type_for_format(format)
            else:
                content_type = \
                    self.get_serializer().get_accept_content_type(request)
            logger.debug('content_type: %s',content_type)
            
            image_keys = [key for key,field in field_hash.items()
                if field.get('display_type', None) == 'image']
            ordered_keys = sorted(field_hash.keys(), 
                key=lambda x: field_hash[x].get('ordinal',key))
            list_fields = [ key for (key,field) in field_hash.items() 
                if( field.get('json_field_type',None) == 'fields.ListField' 
                    or field.get('linked_field_type',None) == 'fields.ListField'
                    or field.get('data_type', None) == 'list' ) ]
            value_templates = {key:field['value_template'] 
                for key,field in field_hash.items() if field.get('value_template', None)}
            data = cursor_generator(
                result,ordered_keys,list_fields=list_fields,
                value_templates=value_templates)
                
            response = None
            if content_type == JSON_MIMETYPE:
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        json_generator(
                            image_generator(data, image_keys, request), 
                            meta, is_for_detail=is_for_detail)))
                response['Content-Type'] = content_type
            
            elif( content_type == XLS_MIMETYPE or
                content_type == XLSX_MIMETYPE ): 

                data = {
                    'data': data 
                }
                response = get_xls_response(
                    data, output_filename, request=request, 
                    title_function=title_function, image_keys=image_keys,
                    list_brackets=list_brackets)

            elif content_type == SDF_MIMETYPE:
                
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        sdf_generator(
                            image_generator(data,image_keys, request), 
                            title_function=title_function)),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.sdf' % output_filename
            
            elif content_type == CSV_MIMETYPE:
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        csv_generator(
                            image_generator(data, image_keys, request), 
                            title_function=title_function, 
                            list_brackets=list_brackets)),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.csv' % output_filename
            else:
                msg = 'unknown content_type: %r' % content_type
                raise BadRequest(msg)
            return response

        except Exception, e:
            logger.exception('on stream response')
            raise e  
