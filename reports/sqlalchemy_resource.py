from __future__ import unicode_literals

from collections import OrderedDict
from functools import wraps
import hashlib
import logging
import os.path
import re
import sys
import urllib

from django.core.cache import cache
import django.db.models.constants
import django.db.models.sql.constants
from django.http.response import StreamingHttpResponse, HttpResponse, Http404
from sqlalchemy import select, asc, text
import sqlalchemy
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.dialects.postgresql import array
from sqlalchemy.sql import and_, or_, not_          
from sqlalchemy.sql import asc, desc, alias, Alias
from sqlalchemy.sql import func
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join, cast
from sqlalchemy.sql.expression import nullsfirst, nullslast
from sqlalchemy.sql.functions import func
import sqlalchemy.sql.sqltypes
from tastypie.exceptions import BadRequest, ImmediateHttpResponse
from tastypie.utils.mime import build_content_type

from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    LIST_BRACKETS, MAX_IMAGE_ROWS_PER_XLS_FILE, MAX_ROWS_PER_XLS_FILE, \
    HTTP_PARAM_RAW_LISTS
from reports.api_base import IccblBaseResource, un_cache
from reports.serialize import XLSX_MIMETYPE, SDF_MIMETYPE, XLS_MIMETYPE
from reports.serialize.csvutils import LIST_DELIMITER_CSV, csv_convert
from reports.serialize.streaming_serializers import sdf_generator, \
    json_generator, get_xls_response, csv_generator, ChunkIterWrapper, \
    cursor_generator, image_generator, closing_iterator_wrapper
from reports.serializers import LimsSerializer
from reports.utils.sqlalchemy_bridge import Bridge

import django.core.signals 
from django.http.request import HttpRequest

logger = logging.getLogger(__name__)

unclosed_connections = []
def connection_close_callback(sender, **kwargs):
    logger.debug("Request finished! %r, %r", sender, kwargs)
    
#     Bridge().get_engine().dispose()
    
#     for c in unclosed_connections:
#         try:
#             c.close()
#         except Exception as e:
#             logger.exception('on conn close...')
#     del unclosed_connections[:]

django.core.signals.request_finished.connect(connection_close_callback)

def _concat(*args):
    '''
    Use as a replacement for sqlalchemy.sql.functions.concat
    - "concat" is not available in postgresql 8.4
    '''
    return func.array_to_string(array([x for x in args]),'')
        
class SqlAlchemyResource(IccblBaseResource):
    '''
    A resource that uses SqlAlchemy to facilitate the read queries:
    - get_list
    - get_detail
    Note: write operations not implemented
    '''
    
    # get a handle to the SqlAlchemy "bridge" for its table registry and 
    # connection handle
    bridge = Bridge()
    
    def __init__(self, *args, **kwargs):
        # get a handle to the SqlAlchemy "bridge" for its table registry and 
        # connection handle
        self.bridge = SqlAlchemyResource.bridge
        
        self.use_cache = True
        
        super(SqlAlchemyResource, self).__init__(*args, **kwargs)
    
    def get_connection(self):
        conn = self.bridge.get_engine().connect()
        unclosed_connections.append(conn)
        return conn
    
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
        _dict = {}
        
        for key in request.GET.keys():
            val = request.GET.getlist(key)
            # Jquery Ajax will post array list params with a "[]" suffix - 20151015
            if '[]' in key and key[-2:] == '[]':
                key = key[:-2]
            if len(val) == 1:
                _dict[key] = val[0]
            else:
                _dict[key] = val
            
        for key in request.POST.keys():
            val = request.POST.getlist(key)
            # Jquery Ajax will post array list params with a "[]" suffix - 20151015
            key = key.replace('[]','')
            if len(val) == 1:
                _dict[key] = val[0]
            else:
                _dict[key] = val
        
        # check for single-valued known list values
        # Note: Jquery Ajax will post array list params with a "[]" suffix - 20151015
        known_list_values = ['includes','exact_fields', 'order_by', 'visibilities',
            'includes[]', 'order_by[]','exact_fields[]', 'visibilities[]']
        for key in known_list_values:
            val = _dict.get(key,[])
            if isinstance(val, basestring):
                _dict[key] = [val]
        
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

    def build_sqlalchemy_columns(self, fields, base_query_tables=None, 
            custom_columns=None):
        '''
        returns an ordered dict of sqlalchemy.sql.schema.Column objects, associated 
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
                if DEBUG_BUILD_COLUMNS:
                    logger.info(str(('field:', key)))
                if key in custom_columns:
                    if DEBUG_BUILD_COLUMNS: 
                        logger.info(str(('custom field', key,custom_columns[key])))
                    columns[key] = custom_columns[key].label(key)
                    continue
                
                if DEBUG_BUILD_COLUMNS: 
                    logger.info(str(('build column', field['key'], field)))
                field_name = field.get('field', None)
                if not field_name:
                    field_name = field['key']
                
                field_table = field.get('table', None)
                
                if not field_table and DEBUG_BUILD_COLUMNS:
                    logger.info(str(('skipping field because there is no '
                        '"field_table" value set',key,field)))
                    continue
                if DEBUG_BUILD_COLUMNS: 
                    logger.info(str(('field', field['key'], 'field_table', field_table )))
                
                if field_table in base_query_tables:
                    # simple case: table.fields already selected in the base query:
                    # just need to specify them
                    if field_name in self.bridge[field_table].c:
                        col = self.bridge[field_table].c[field_name]
                    else:
                        raise Exception(str(('field', field_name, 
                            'not found in table', field_table)))
                    col = col.label(key)
                    columns[key] = col
                    
                elif field.get('linked_field_value_field', None):
                    link_table = field['table']
                    link_table_def = self.bridge[link_table]
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
                            meta_table_def = self.bridge['metahash']
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
                            meta_table_def = self.bridge['metahash']
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
            if DEBUG_BUILD_COLUMNS: logger.info(str(('columns', columns.keys())))
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
            else:
                order_clause = nullsfirst(asc(column(field_name)))
            if field_name in visible_fields:
                order_clauses.append(order_clause)
            else:
                logger.warn(str(('order_by field not in visible fields, skipping: ', 
                    order_by )))
        if DEBUG_ORDERING:
            logger.info('order_clauses %s',order_clauses)     
        return order_clauses
    
    @staticmethod
    def filter_value_to_python(value, param_hash, filter_expr, filter_type):
        """
        Turn the string ``value`` into a python object.
        - copied from TP
        """
        # Simple values
        if value is 1:
            value = 1
        elif value in ['true', 'True', True]:
            value = True
        elif value in ['false', 'False', False]:
            value = False
        elif value in ('nil', 'none', 'None', None):
            value = None

        # Split on ',' if not empty string and either an in or range filter.
        if filter_type in ('in', 'range') and len(value):
            if value and hasattr(value, '__iter__'):
                # value is already a list
                pass
            elif hasattr(param_hash, 'getlist'):
                value = []
                for part in param_hash.getlist(filter_expr):
                    value.extend(part.split(LIST_DELIMITER_URL_PARAM))
            else:
                value = value.split(LIST_DELIMITER_URL_PARAM)
        logger.debug(str(('filter value', filter_expr, value)))
        return value

    @staticmethod
    def build_sqlalchemy_filters(schema, param_hash={}):
        DEBUG_FILTERS = False or logger.isEnabledFor(logging.DEBUG)
        
        if DEBUG_FILTERS: 
            logger.info('build_sqlalchemy_filters: param_hash %s', param_hash)

        # ordinary filters
        (filter_expression, filter_fields) = \
            SqlAlchemyResource.build_sqlalchemy_filters_from_hash(schema, param_hash)
        
        # Treat the nested "search_data" as sets of params to be OR'd together,
        # then AND'd with the regular filters (if any)
        search_data = param_hash.get('search_data', None)
        if search_data:
            logger.info('search_data: %r', search_data)
            if search_data and isinstance(search_data, basestring):
                # standard, convert single valued list params
                search_data = [search_data]
        
            # each item in the search data array is a search hash, to be or'd
            search_expressions = []
            filter_fields = set(filter_fields)
            for search_hash in search_data:
                logger.info('search_hash: %s' % search_hash)
                (search_expression, search_fields) = SqlAlchemyResource.\
                    build_sqlalchemy_filters_from_hash(schema,search_hash)
                search_expressions.append(search_expression)
                filter_fields.update(search_fields)

            if len(search_expressions) > 1:
                search_expressions = or_(*search_expressions)
            else:
                search_expressions = search_expressions[0]
            if filter_expression is not None:
                filter_expression = and_(search_expressions,filter_expression)
            else: 
                filter_expression = search_expressions
                
        if DEBUG_FILTERS: 
            logger.info('filter_expression: %s, filter_fields: %s',
                filter_expression, filter_fields)
        
        return (filter_expression,filter_fields)
    
    @staticmethod
    def build_sqlalchemy_filters_from_hash(schema, param_hash):
        '''
        Attempt to create a SqlAlchemy whereclause out of django style filters:
        - field_name__filter_expression
        '''
        DEBUG_FILTERS = False or logger.isEnabledFor(logging.DEBUG)
        logger.debug('build_sqlalchemy_filters_from_hash %r' % param_hash)
        lookup_sep = django.db.models.constants.LOOKUP_SEP

        if param_hash is None:
            return (None,None)
        
        try:
            expressions = []
            filtered_fields = []
            
            _values = [] # store values for debug
            for filter_expr, value in param_hash.items():
                if lookup_sep not in filter_expr:
                    continue;
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
                if not field_name in schema['fields']:
                    logger.warn(str(('unknown filter field', field_name, filter_expr)))
                    continue
                if DEBUG_FILTERS:
                    logger.info('build filter expr: field_name, %r, '
                        'filter_type: %r value: %r', 
                        field_name, filter_type, value)
    
                field = schema['fields'][field_name]
                
                value = SqlAlchemyResource.filter_value_to_python(
                    value, param_hash, filter_expr, filter_type)
                if DEBUG_FILTERS:
                    logger.info('value : %r', value)
                _values.append(value)
                # TODO: all of the Django query terms:
                # QUERY_TERMS = set([
                #     'exact', 'iexact', 'contains', 'icontains', 'gt', 'gte', 
                # 'lt', 'lte', 'in',
                #     'startswith', 'istartswith', 'endswith', 'iendswith', 
                # 'range', 'year',
                #     'month', 'day', 'week_day', 'hour', 'minute', 'second', 
                # 'isnull', 'search',
                #     'regex', 'iregex',
                # ])
                
                if DEBUG_FILTERS:
                    logger.info(str(('find filter', field_name, filter_type, value, type(value))))
                
                expression = None
                col = column(field_name)
                if field['data_type'] in ['integer','float','decimal']:
                    col = cast(col,sqlalchemy.sql.sqltypes.Numeric)
                elif field['data_type'] == 'boolean':
                    col = cast(col,sqlalchemy.sql.sqltypes.Boolean)
                if field['data_type'] == 'string':
                    col = cast(col,sqlalchemy.sql.sqltypes.Text)

                if filter_type in ['exact','eq']:
                    if field['data_type'] == 'string':
                        value = str(value)
                    expression = col == value
                    if field['data_type'] == 'list':
                        expression = text(
                            "'%s'=any(string_to_array(%s,'%s'))"
                            % (value,field_name,LIST_DELIMITER_SQL_ARRAY))
                elif filter_type == 'about':
                    decimals = 0
                    if '.' in value:
                        decimals = len(value.split('.')[1])
                    expression = func.round(col,decimals) == value
                    if DEBUG_FILTERS:
                        logger.info(str(('create "about" expression for term:',
                            filter_expr,
                            value,'decimals',decimals)))
                elif filter_type == 'contains':
                    if field['data_type'] == 'string':
                        value = str(value)
                    expression = col.contains(value)
                elif filter_type == 'icontains':
                    if field['data_type'] == 'string':
                        value = str(value)
                    expression = col.ilike('%{value}%'.format(
                        value=value))
                elif filter_type == 'lt':
                    expression = col < value
                elif filter_type == 'lte':
                    expression = col <= value
                elif filter_type == 'gt':
                    expression = col > value
                elif filter_type == 'gte':
                    expression = col >= value
                elif filter_type == 'is_blank':
                    if field['data_type'] == 'string':
                        col = func.trim(col)
                    if value and str(value).lower() == 'true':
                        expression = col == None 
                        if field['data_type'] == 'string':
                            expression = col == ''
                    else:
                        expression = col != None
                        if field['data_type'] == 'string':
                            col = func.trim(col)
                        if ( field['data_type'] == 'string' 
                             or field['data_type'] == 'list' ):
                            expression = col != ''
                elif filter_type == 'is_null':
                    if value and str(value).lower() == 'true':
                        expression = col == None 
                    else:
                        expression = col != None
                elif filter_type == 'in':
                    if field['data_type'] == 'list':
                        # NOTE: for the list type, interpret "in" as any of the 
                        # given values are in the field
                        temp_expressions = [] 
                        for _val in value:
                            temp_expressions.append(col.ilike('%{value}%'.format(
                                value=_val)))
                        expression = or_(*temp_expressions)
                    else:
                        expression = col.in_(value)
                elif filter_type == 'ne':
                    if field['data_type'] == 'string':
                        value = str(value)
                    expression = col != value
                elif filter_type == 'range':
                    if len(value) != 2:
                        logger.error(str((
                            'value for range expression must be list of length 2', 
                            field_name, filter_expr, value)))
                        continue
                    else:
                        expression = col.between(
                            value[0],value[1],symmetric=True)
                else:
                    logger.error(str(('--- unknown filter type: ', 
                        field_name, filter_type,
                        'filter_expr',filter_expr )))
                    continue
    
                if inverted:
                    expression = not_(expression)
                
                logger.debug('filter_expr: %r' % filter_expr)
                expressions.append(expression)
                filtered_fields.append(field_name)
                
            logger.debug('filtered_fields: %s', filtered_fields)
            if DEBUG_FILTERS:
                logger.info(str(('values', _values)))
                
            if len(expressions) > 1: 
                return (and_(*expressions), filtered_fields)
            elif len(expressions) == 1:
                return (expressions[0], filtered_fields) 
            else:
                return (None, filtered_fields)
        except Exception, e:
            logger.exception('on build_sqlalchemy_filters_from_hash')
            raise e   
        

    def search(self, request, **kwargs):
        '''
        Implement a special search view to get around Tastypie deserialization
        methods on POST-list.
        Note: could implement a special tastypie serializer and "post_list"; but
        choosing not to couple with the framework in this way here.
        '''
         
        DEBUG_SEARCH = False or logger.isEnabledFor(logging.DEBUG)
         
        search_ID = kwargs['search_ID']
         
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
 
        search_data = param_hash.get('search_data', None)
        
        # NOTE: unquote serves the purpose of an application/x-www-form-urlencoded 
        # deserializer
        search_data = urllib.unquote(search_data)
         
        if search_data:
            # cache the search data on the session, to support subsequent requests
            # to download or modify
            request.session[search_ID] = search_data  
         
        if not search_data:
            if search_ID in request.session:
                search_data = request.session[search_ID]
            else:
                raise ImmediateHttpResponse(
                    response=self.error_response(request, 
                        { 'search_data for id missing: ' + search_ID: 
                            self._meta.resource_name + 
                                '.search requires a "search_data" param'}))
         
        search_data = json.loads(search_data)   
        param_hash['search_data'] = search_data
        if DEBUG_SEARCH:
            logger.info(str(('search', param_hash, kwargs)))
 
        response = self.build_list_response(request,param_hash=param_hash, **kwargs)
        
        # Because this view bypasses the IccblBaseResource.dispatch method, we
        # are implementing the downloadID cookie here, for now.
        downloadID = param_hash.get('downloadID', None)
        if downloadID:
            logger.info(str(('set cookie','downloadID', downloadID )))
            response.set_cookie('downloadID', downloadID)
        else:
            logger.info(str(('no downloadID')))
        
        return response

    def _get_list_response(self,request,**kwargs):
        '''
        Return a deserialized list of dicts
        '''
        includes = kwargs.pop('includes', '*')
        try:
            response = self.get_list(
                request,
                format='json',
                includes=includes,
                **kwargs)
            _data = self._meta.serializer.deserialize(
                request,
                LimsSerializer.get_content(response), format='application/json')
            _data = _data[self._meta.collection_name]
            logger.debug(' data: %s'% _data)
            return _data
        except Http404:
            return []
        except Exception as e:
            # FIXME: temporary travis debug
            logger.exception('on get list')
            return []
        
    @un_cache
    def _get_detail_response_internal(self, **kwargs):
        request = HttpRequest()
        class User:
            @staticmethod
            def is_superuser():
                return true
        request.user = User
#         temp = self.use_cache
#         self.use_cache = False
        result = self._get_detail_response(request, **kwargs)
#         self.use_cache = temp
        return result
    
    def _get_detail_response(self,request,**kwargs):
        '''
        Return the detail response as a dict
        '''
        logger.info('_get_detail_response: %r', kwargs)
        try:
            response = self.get_detail(
                request,
                format='json',
                includes='*',
                **kwargs)
            _data = self._meta.serializer.deserialize(
                request,
                LimsSerializer.get_content(response), format='application/json')
            logger.debug(' data: %s'% _data)
            return _data
        except Http404:
            return []
        
    def get_list(self, request, **kwargs):
        '''
        Override the Tastypie/Django ORM get_list method - list/reporting 
        operations will be handled using SqlAlchemy
        '''
        raise NotImplemented(str((
            'get_list must be implemented for the SqlAlchemyResource', 
            self._meta.resource_name)) )
        
    def build_list_response(self,request, param_hash={}, **kwargs):
        raise NotImplemented(str((
            'get_list_response must be implemented for the SqlAlchemyResource', 
            self._meta.resource_name)) )
    
    def _cached_resultproxy(self, conn, stmt, count_stmt, param_hash, limit, offset):
        ''' 
        Cache for resultsets:
        - Always returns the cache object with a resultset, either from the cache,
        or executed herein.
        
        NOTE: limit and offset are included because this version of sqlalchemy
        does not support printing of them with the select.compile() function.
        
        TODO: cache clearing on database writes
        '''
        DEBUG_CACHE = True or logger.isEnabledFor(logging.DEBUG)
        if limit == 0:
            raise Exception('limit for caching must be >0')
        
        prefetch_number = 5
        
        if limit==1:
            prefetch_number = 1
        
        try:
            # use a hexdigest because statements can contain problematic chars 
            # for the memcache
            m = hashlib.md5()
            compiled_stmt = str(stmt.compile(compile_kwargs={"literal_binds": True}))
            logger.debug('compiled_stmt %s',compiled_stmt)
            if 'limit' in compiled_stmt.lower():
                # remove limit and offset; will calculate
                compiled_stmt = compiled_stmt[:compiled_stmt.lower().rfind('limit')]
            logger.debug('compiled_stmt for hash key: %s', compiled_stmt)
            key_digest = '%s_%s_%s' %(compiled_stmt, str(limit), str(offset))
            m.update(key_digest)
            key = m.hexdigest()
            logger.debug('hash key: digest: %s, key: %s, limit: %s, offset: %s', 
                key_digest, key, limit, offset)
            cache_hit = cache.get(key)
            if cache_hit:
                if ('stmt' not in cache_hit or
                        cache_hit['stmt'] != compiled_stmt):
                    cache_hit = None
                    logger.warn(str(('cache collision for key', key, stmt)))
            
            if not cache_hit:
                # Note: if no cache hit, then retrive limit*n results, and 
                # cache several iterations at once.
                new_limit = limit * prefetch_number 
                if DEBUG_CACHE:
                    logger.info('no cache hit, create cache, limit: %s, new limit for caching: %s',
                        limit, new_limit)
                if new_limit > 0:
                    stmt = stmt.limit(new_limit)
                logger.info('no cache hit, executing stmt')
                resultset = conn.execute(stmt)
                prefetched_result = [dict(row) for row in resultset] if resultset else []
                logger.info('executed stmt')
                
#                 if DEBUG_CACHE:
                logger.info('no cache hit, execute count')
                if limit == 1:
                    count = 1
                else:
                    count = conn.execute(count_stmt).scalar()
#                 if DEBUG_CACHE:
                logger.info('count: %s', count)
                
                # now fill in the cache with the prefetched sets or rows
                for y in range(prefetch_number):
                    new_offset = offset + limit*y;
                    _start = limit*y
                    logger.info('new_offset: %d, start: %d, len: %d', 
                        new_offset, _start, len(prefetched_result))
                    if _start < len(prefetched_result):
                        key_digest = '%s_%s_%s' %(compiled_stmt, str(limit), str(new_offset))
                        m = hashlib.md5()
                        m.update(key_digest)
                        key = m.hexdigest()
                        _result = prefetched_result[_start:_start+limit]
                        _cache = {
                            'stmt': compiled_stmt,
                            'cached_result': _result,
                            'count': count }
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
                logger.info(str(('cache hit')))   
                
            return cache_hit
 
        except Exception, e:
            logger.exception('on conn execute')
            raise e   

    def stream_response_from_statement(self, request, stmt, count_stmt, 
            output_filename, field_hash={}, param_hash={}, 
            rowproxy_generator=None, is_for_detail=False,
            downloadID=None, title_function=None, use_caching=True, meta=None ):
        DEBUG_STREAMING = False or logger.isEnabledFor(logging.DEBUG)

        if DEBUG_STREAMING:
            logger.info('stream_response_from_statement: %r' % param_hash)
        limit = param_hash.get('limit', 0)        
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
        
#         conn = self.get_connection()
        conn = self.bridge.get_engine().connect()
        
        try:
            logger.debug('offset: %s, limit: %s', offset, limit)
        
            if DEBUG_STREAMING:
                logger.info('stmt: %s, param_hash: %s ', 
                    str(stmt.compile(compile_kwargs={"literal_binds": True})), 
                    param_hash)
                logger.info(str(('count stmt', str(count_stmt))))
            
            desired_format = self.get_serialize_format(request, **param_hash)
            logger.debug('---- desired_format: %r, hash: %r', desired_format, param_hash)
            result = None
            if desired_format == 'application/json':
                logger.info('streaming json, use_caching: %r, limit: %d', use_caching, limit)
#                 if not is_for_detail and use_caching and self.use_cache and limit > 0:
                if use_caching and self.use_cache and limit > 0:
                    cache_hit = self._cached_resultproxy(
                        conn, stmt, count_stmt, param_hash, limit, offset)
                    if cache_hit:
                        logger.info('cache hit')
                        result = cache_hit['cached_result']
                        count = cache_hit['count']
                    else:
                        # cache routine should always return a cache object
                        logger.error('error, cache not set: execute stmt')
                        count = conn.execute(count_stmt).scalar()
                        result = conn.execute(stmt)
                    logger.info(str(('====count====', count)))
                    
                else:
                    logger.info('execute count stmt...')
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
    
                # TODO: create a short-circuit if count==0
                # if count == 0:
                #    raise ImmediateHttpResponse(
                #        response=self.error_response(
                #            request, {'empty result': 'no records found'},
                #            response_class=HttpNotFound))
                if DEBUG_STREAMING:
                    logger.info('json setup done, meta: %r', meta)
    
            else: # not json
            
                logger.info('excute stmt')
                result = conn.execute(stmt)
                logger.info('excuted stmt')
                
                logger.info(str(('rowproxy_generator', rowproxy_generator)))
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
                meta=meta)
        except Exception, e:
            logger.exception('on stream response')
            raise e
#         finally:
#             conn.close()          
        
    def stream_response_from_cursor(self,request,result,output_filename,
            field_hash={}, param_hash={}, 
            is_for_detail=False,
            downloadID=None, title_function=None, meta=None):
          
        try:
                    
            list_brackets = LIST_BRACKETS
            if request.GET.get(HTTP_PARAM_RAW_LISTS, False):
                list_brackets = None
    
            desired_format = self.get_serialize_format(request, **param_hash)
            content_type=build_content_type(desired_format)
            logger.debug('desired_format: %s, content_type: %s', 
                desired_format, content_type)
            
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
            if desired_format == 'application/json':
                
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        json_generator(
                            image_generator(data, image_keys, request), 
                            meta, is_for_detail=is_for_detail)))
                response['Content-Type'] = content_type
            
            elif( desired_format == XLS_MIMETYPE or
                desired_format == XLSX_MIMETYPE ): 

                data = {
                    'data': data 
                }
                response = get_xls_response(
                    data, output_filename, request=request, 
                    title_function=title_function, image_keys=image_keys,
                    list_brackets=list_brackets)

            elif desired_format == SDF_MIMETYPE:
                
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        sdf_generator(
                            image_generator(data,image_keys, request), 
                            title_function=title_function)),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.sdf' % output_filename
            
            elif desired_format == 'text/csv':
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
                msg = 'unknown format: %r' % desired_format
                raise BadRequest(msg)
            return response

        except Exception, e:
            logger.exception('on stream response')
            raise e  
