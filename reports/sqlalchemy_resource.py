
from collections import OrderedDict
import hashlib
import logging
import os.path
import sys
import urllib
from functools import wraps

from django.core.cache import cache
import django.db.models.constants
import django.db.models.sql.constants
from django.http.response import StreamingHttpResponse, HttpResponse
from sqlalchemy import select, asc, text
import sqlalchemy
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.sql import and_, or_, not_          
from sqlalchemy.sql import asc, desc, alias, Alias
from sqlalchemy.sql import func
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join
from sqlalchemy.sql.expression import nullsfirst, nullslast
from tastypie.exceptions import BadRequest, ImmediateHttpResponse
from tastypie.http import HttpNotFound
from tastypie.resources import Resource
from tastypie.utils.mime import build_content_type

from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    LIST_BRACKETS, MAX_IMAGE_ROWS_PER_XLS_FILE, MAX_ROWS_PER_XLS_FILE, \
    LIST_DELIMITER_XLS, LIST_DELIMITER_CSV, HTTP_PARAM_RAW_LISTS
from reports.utils.sqlalchemy_bridge import Bridge
from reports.utils.streaming_serializers import sdf_generator, json_generator, \
    get_xls_response, csv_generator, ChunkIterWrapper


logger = logging.getLogger(__name__)

def un_cache(_func):
    '''
    Wrapper function to disable caching for 
    SQLAlchemyResource.stream_response_from_statement
    ''' 
    @wraps(_func)
    def _inner(self, *args, **kwargs):
        logger.warn(str(('decorator un_cache', self, kwargs )))
        SqlAlchemyResource.clear_cache(self)
        SqlAlchemyResource.set_caching(self,False)
        result = _func(self, *args, **kwargs)
        SqlAlchemyResource.set_caching(self,True)
        logger.warn(str(('decorator un_cache done', kwargs )))
        return result

    return _inner


class SqlAlchemyResource(Resource):
    '''
    A resource that uses SqlAlchemy to facilitate the read queries:
    - get_list
    - get_detail
    Note: write operations not implemented
    '''
    
#     class Meta:
#         queryset = Reagent.objects.all()
#         serializer = LimsSerializer() # still have a serializer for error response

    # get a handle to the SqlAlchemy "bridge" for its table registry and 
    # connection handle
    bridge = Bridge()
    
    def __init__(self, *args, **kwargs):
        # get a handle to the SqlAlchemy "bridge" for its table registry and 
        # connection handle
        self.bridge = SqlAlchemyResource.bridge
        
        self.use_cache = True
        
        super(SqlAlchemyResource, self).__init__(*args, **kwargs)
    
    def set_caching(self,use_cache):
        self.use_cache = use_cache

    @classmethod
    def wrap_statement(cls, stmt, order_clauses, filter_expression):
        '''
        @param stmt - a sqlalchemy.sql.expression.join instance
        @param order_clauses - list of sqlalchemy.sql.expression.column
        @param filter_expression - sqlalchemy.whereclause
        '''
        if order_clauses:
            logger.info(str(('order_clauses', [str(c) for c in order_clauses])))
            _alias = Alias(stmt)
            stmt = select([text('*')]).select_from(_alias)
            stmt = stmt.order_by(*order_clauses)
        if filter_expression is not None:
            logger.info(str(('filter_expression', str(filter_expression))))
            if not order_clauses:
                _alias = Alias(stmt)
                logger.info(str(('filter_expression', str(filter_expression))))
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
            if len(val) == 1:
                _dict[key] = val[0]
            else:
                _dict[key] = val
            
        for key in request.POST.keys():
            val = request.POST.getlist(key)
            if len(val) == 1:
                _dict[key] = val[0]
            else:
                _dict[key] = val
        
        # check for single-valued known list values
        known_list_values = ['includes', 'order_by']
        for key in known_list_values:
            val = _dict.get(key,[])
            if isinstance(val, basestring):
                _dict[key] = [val]
        
        _dict['desired_format'] = self.get_format(request)
        
        return _dict    
    
    def get_visible_fields(self, schema_fields, filter_fields, manual_field_includes,
                           is_for_detail=False, visibilities=[]):
        '''
        Construct an ordered dict of schema fields that are visible, based on
        - "list" in schema field["visibility"]
        - field key in (filter_fields or manual_field_includes)
        - field key in a schema field['dependencies'] 
        
        TODO: this method can be static
        TODO: this method is not SqlAlchemy specific
        '''
        logger.info(str(('get_visible_fields: field_hash initial: ', 
            schema_fields.keys() )))
        try:
            if visibilities:
                visibilities = set(visibilities)
            else:
                visibilities = set(['list'])
            if is_for_detail:
                visibilities.add('detail')
                # also return the edit fields, let the UI filter them
                # expedient, so the model does not have to be reloaded to edit
                # TODO: review; security issue
                visibilities.add('edit')
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
                dependency_fields.update(field.get('dependencies',[]))
            logger.info(str(('dependency_fields', dependency_fields)))
            if dependency_fields:
                temp.update({ key:field 
                    for key,field in schema_fields.items() if key in dependency_fields })
            
            # filter_fields
            if filter_fields:
                temp.update({ key:field 
                    for key,field in schema_fields.items() if key in filter_fields })
             
            field_hash = OrderedDict(sorted(temp.iteritems(), 
                key=lambda x: x[1].get('ordinal',999))) 
    
            logger.info(str(('field_hash final: ', field_hash.keys() )))
            return field_hash

        except Exception, e:
            logger.exception('on get_visible_fields')
            raise e 

    def build_sqlalchemy_columns(self, fields, base_query_tables=[], 
            custom_columns={}):
        '''
        returns an array of sqlalchemy.sql.schema.Column objects, associated 
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
                
                if not field_table:
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
                    logger.warn(str((
                        'field is not in the base tables or in a linked field, '
                        'and is not custom', key)))
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
        This method borrows from tastypie.resources.ModelResource.apply_sorting
        @param order_params passed as list in the request.GET hash
        '''
        
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
        logger.debug(str(('order_clauses', order_clauses)))     
        return order_clauses
    
    @staticmethod
    def filter_value_to_python(value, param_hash, filter_expr, filter_type):
        """
        Turn the string ``value`` into a python object.
        - copied from TP
        """
        # Simple values
        if value in ['true', 'True', True]:
            value = True
        elif value in ['false', 'False', False]:
            value = False
        elif value in ('nil', 'none', 'None', None):
            value = None

        # Split on ',' if not empty string and either an in or range filter.
        if filter_type in ('in', 'range') and len(value):
            if hasattr(param_hash, 'getlist'):
                value = []

                for part in param_hash.getlist(filter_expr):
                    value.extend(part.split(LIST_DELIMITER_URL_PARAM))
            else:
                value = value.split(LIST_DELIMITER_URL_PARAM)
        logger.debug(str(('filter value', filter_expr, value)))
        return value

    @staticmethod
    def build_sqlalchemy_filters(schema, param_hash={}, **kwargs):
        logger.info(str(('param_hash', param_hash, 'kwargs', kwargs)))

        # ordinary filters
        (filter_expression, filter_fields) = \
            SqlAlchemyResource.build_sqlalchemy_filters_from_hash(schema, param_hash)
        
        # Treat the nested "search_data" as sets of params to be OR'd together,
        # then AND'd with the regular filters (if any)
        search_data = param_hash.get('search_data', None)
        if search_data:
            logger.info(str(('search_data', search_data)))
            if search_data and isinstance(search_data, basestring):
                # standard, convert single valued list params
                search_data = [search_data]
        
            # each item in the search data array is a search hash, to be or'd
            search_expressions = []
            filter_fields = set(filter_fields)
            for search_hash in search_data:
                logger.info(str(('search_hash', search_hash)))
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
                
        logger.info(str(('filter_expression', filter_expression, filter_fields)))
        
        return (filter_expression,filter_fields)
    
    @staticmethod
    def build_sqlalchemy_filters_from_hash(schema, param_hash):
        '''
        Attempt to create a SqlAlchemy whereclause out of django style filters:
        - field_name__filter_expression
        '''
        DEBUG_FILTERS = False or logger.isEnabledFor(logging.DEBUG)
        logger.info('build_sqlalchemy_filters_from_hash %r' % param_hash)
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
                    logger.warn(str((
                        'filter expression must be of the form '
                        '"field_name__expression"',
                        filter_expr, filter_bits)))
                field_name = filter_bits[0]
                
                inverted = False
                if field_name and field_name[0] == '-':
                    inverted = True
                    field_name = field_name[1:]
                            
                filter_type = filter_bits[1]
                if not field_name in schema['fields']:
                    logger.warn(str(('unknown filter field', field_name, filter_expr)))
                    continue
    
                field = schema['fields'][field_name]
                
                value = SqlAlchemyResource.filter_value_to_python(
                    value, param_hash, filter_expr, filter_type)
                
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
                    logger.info(str(('find filter', field_name, filter_type, value)))
                
                expression = None
                if filter_type in ['exact','eq']:
                    expression = column(field_name)==value
                    if field['data_type'] == 'list':
                        expression = text(
                            "'%s'=any(string_to_array(%s,'%s'))"
                            % (value,field_name,LIST_DELIMITER_SQL_ARRAY))
                elif filter_type == 'about':
                    decimals = 0
                    if '.' in value:
                        decimals = len(value.split('.')[1])
                    expression = func.round(
                        sqlalchemy.sql.expression.cast(column(field_name),
                            sqlalchemy.types.Numeric),decimals) == value
                    if DEBUG_FILTERS:
                        logger.info(str(('create "about" expression for term:',
                            filter_expr,
                            value,'decimals',decimals)))
                elif filter_type == 'contains':
                    expression = column(field_name).contains(value)
                elif filter_type == 'icontains':
                    expression = column(field_name).ilike('%{value}%'.format(
                        value=value))
                elif filter_type == 'lt':
                    expression = column(field_name) < value
                elif filter_type == 'lte':
                    expression = column(field_name) <= value
                elif filter_type == 'gt':
                    expression = column(field_name) > value
                elif filter_type == 'gte':
                    expression = column(field_name) >= value
                elif filter_type == 'is_blank':
                    col = column(field_name)
                    if field['data_type'] == 'string':
                        col = func.trim(col)
                    if value and str(value).lower() == 'true':
                        expression = col == None 
                        if field['data_type'] == 'string':
                            expression = col == ''
                    else:
                        expression = col != None
                        if ( field['data_type'] == 'string' 
                             or field['data_type'] == 'list' ):
                            expression = col != ''
                elif filter_type == 'is_null':
                    col = column(field_name)
                    if field['data_type'] == 'string':
                        col = func.trim(col)
                    if value and str(value).lower() == 'true':
                        expression = col == None 
                    else:
                        expression = col != None
                elif filter_type == 'in':
                    expression = column(field_name).in_(value)
                elif filter_type == 'ne':
                    expression = column(field_name) != value
                elif filter_type == 'range':
                    if len(value) != 2:
                        logger.error(str((
                            'value for range expression must be list of length 2', 
                            field_name, filter_expr, value)))
                        continue
                    else:
                        expression = column(field_name).between(
                            value[0],value[1],symmetric=True)
                else:
                    logger.error(str(('--- unknown filter type: ', 
                        field_name, filter_type,
                        'filter_expr',filter_expr )))
                    continue
    
                if inverted:
                    expression = not_(expression)
                
                logger.info(str(('filter_expr',filter_expr,
                    'expression',str(expression) )))
                expressions.append(expression)
                filtered_fields.append(field_name)
                
            logger.info(str(('filtered_fields', filtered_fields)))
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
    
    def clear_cache(self):
        logger.debug('clearing the cache: resource: %s' % self._meta.resource_name)
        cache.clear()

    def _cached_resultproxy(self, stmt, count_stmt, param_hash, limit, offset):
        ''' 
        ad-hoc cache for some resultsets:
        - Always returns the cache object with a resultset, either from the cache,
        or executed herein.
        
        NOTE: limit and offset are included because this version of sqlalchemy
        does not support printing of them with the select.compile() function.
        
        TODO: cache clearing on database writes
        '''
        if limit == 0:
            raise Exception('limit for caching must be >0')
        
        prefetch_number = 5
        
        conn = self.bridge.get_engine().connect()
        try:
            # use a hexdigest because statements can contain problematic chars 
            # for the memcache
            m = hashlib.md5()
            compiled_stmt = str(stmt.compile(compile_kwargs={"literal_binds": True}))
            logger.debug(str(('compiled_stmt',compiled_stmt)))
            if 'limit' in compiled_stmt.lower():
                # remove limit and offset; will calculate
                compiled_stmt = compiled_stmt[:compiled_stmt.lower().rfind('limit')]
            logger.debug(str(('compiled_stmt for hash key', compiled_stmt)))
            key_digest = '%s_%s_%s' %(compiled_stmt, str(limit), str(offset))
            m.update(key_digest)
            key = m.hexdigest()
            logger.debug(str(('hash key:', key_digest, key, limit, offset)))
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
                logger.info(str(('limit', limit, 'new limit for caching', new_limit)))
                if new_limit > 0:
                    stmt = stmt.limit(new_limit)
                logger.info('no cache hit, executing stmt')
                resultset = conn.execute(stmt)
                prefetched_result = [dict(row) for row in resultset] if resultset else []
                logger.info('executed stmt')
                
                logger.info(str(('no cache hit, execute count')))
                count = conn.execute(count_stmt).scalar()
                logger.info(str(('count', count)))
                
                # now fill in the cache with the prefetched sets or rows
                for y in range(prefetch_number):
                    new_offset = offset + limit*y;
                    _start = limit*y
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
                        logger.info(str(('add to cache, key', key, limit, new_offset)))
                        cache.set( key, _cache, None)
                        if y == 0:
                            cache_hit = _cache
                    else:
                        logger.info(str(('not caching: prefetched length: ',
                            len(prefetched_result), _start, 'end')))
                        break
                logger.info(str(('cached iterations:', y+1)))
#                 cached_result = [dict(row) for row in resultset] if resultset else []
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
        logger.info(str(('stream_response_from_statement',param_hash)))

        try:
            limit = param_hash.get('limit', 0)        
            try:
                limit = int(limit)
            except ValueError:
                raise BadRequest(
                    "Invalid limit '%s' provided. Please provide a positive integer." 
                    % limit)
            if limit > 0:    
                stmt = stmt.limit(limit)
    
            offset = param_hash.get('offset', 0 )
            try:
                offset = int(offset)
            except ValueError:
                raise BadRequest(
                    "Invalid offset '%s' provided. Please provide a positive integer." 
                    % offset)
            if offset < 0:    
                offset = -offset
            stmt = stmt.offset(offset)
    
            logger.info(str(('offset', offset, 'limit', limit)))
            conn = self.bridge.get_engine().connect()
            
            logger.info(str(('stmt', 
                str(stmt.compile(compile_kwargs={"literal_binds": True})), 
                'param_hash', param_hash)))
            if DEBUG_STREAMING:
                logger.info(str(('count stmt', str(count_stmt))))
            
            desired_format = param_hash.get('desired_format',self.get_format(request))
            result = None
            if desired_format == 'application/json':
                logger.info(str(('streaming json')))
                if not is_for_detail and use_caching and self.use_cache and limit > 0:
                    cache_hit = self._cached_resultproxy(
                        stmt, count_stmt, param_hash, limit, offset)
                    if cache_hit:
                        logger.info('cache hit')
                        result = cache_hit['cached_result']
                        count = cache_hit['count']
                    else:
                        # cache routine should always return a cache object
                        logger.error(str(('error, cache not set: execute stmt')))
                        count = conn.execute(count_stmt).scalar()
                        result = conn.execute(stmt)
                    logger.info(str(('====count====', count)))
                    
                else:
                    logger.info(str(('execute stmt')))
                    count = conn.execute(count_stmt).scalar()
                    logger.info('excuted count stmt')
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
                    logger.info(str(('meta', meta)))
    
            else: # not json
            
                logger.info('excute stmt')
                result = conn.execute(stmt)
                logger.info('excuted stmt')
                
                logger.info(str(('rowproxy_generator', rowproxy_generator)))
                if rowproxy_generator:
                    result = rowproxy_generator(result)
                    # FIXME: test this for generators other than json generator        
    
        except Exception, e:
            logger.exception('on stream response')
            raise e          
        return self.stream_response_from_cursor(request, result, output_filename, 
            field_hash=field_hash, 
            param_hash=param_hash, 
            is_for_detail=is_for_detail, 
            downloadID=downloadID, 
            title_function=title_function, 
            meta=meta)
        
    def stream_response_from_cursor(self,request,result,output_filename,
            field_hash={}, param_hash={}, 
            is_for_detail=False,
            downloadID=None, title_function=None, meta=None):
          
        logger.info(str(('meta', meta, 'session', request.session, 
            request.session.session_key)))

        try:
                    
            list_brackets = LIST_BRACKETS
            if request.GET.get(HTTP_PARAM_RAW_LISTS, False):
                list_brackets = None
    
            desired_format = param_hash.get(
                'desired_format',self.get_format(request))
            content_type=build_content_type(desired_format)
            logger.info(str(('desired_format', desired_format, 
                'content_type', content_type)))
            
            response = None
            
            if desired_format == 'application/json':
                
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        json_generator(
                            result,meta,request,
                            is_for_detail=is_for_detail,field_hash=field_hash)))
                response['Content-Type'] = content_type
            
            
            elif( desired_format == 'application/xls' or
                desired_format == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'): 
                # NOTE: will respond with xlsx for both desired formats
                
                response = get_xls_response(result,output_filename,request, 
                    field_hash=field_hash, title_function=title_function,
                    list_brackets=list_brackets)
                    
            elif desired_format == 'chemical/x-mdl-sdfile':
                
                response = StreamingHttpResponse(
                    ChunkIterWrapper(sdf_generator(result, field_hash)),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.sdf' % output_filename
            
            elif desired_format == 'text/csv':
    
                response = StreamingHttpResponse(
                    csv_generator(
                        result,request,field_hash=field_hash,
                        title_function=title_function,list_brackets=list_brackets),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.csv' % output_filename
            else:
                msg = str(('unknown format', desired_format, output_filename))
                logger.error(msg)
                raise ImmediateHttpResponse(msg)
    
            return response

        except Exception, e:
            logger.exception('on stream response')
            raise e  
