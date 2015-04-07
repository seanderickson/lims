# Tastypie Resource extensions that use SqlAlchemy to generate get_list, get_detail data
import cStringIO
from collections import defaultdict, OrderedDict
import csv
import hashlib
import io
import json
import math
import os.path
import re
import shutil
import sys
import time

from PIL import Image

from django.core.urlresolvers import resolve
from django.http.response import HttpResponseBase
from django.http.response import StreamingHttpResponse, HttpResponse
import django.db.models.constants
from django.core.serializers.json import DjangoJSONEncoder
import django.db.models.sql.constants
from sqlalchemy import select, asc, text
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.sql import and_, or_, not_          
from sqlalchemy.sql import asc, desc, alias, Alias
from sqlalchemy.sql import func
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join
from sqlalchemy.sql.expression import nullsfirst, nullslast
from tastypie.resources import Resource
from tastypie.utils.mime import build_content_type
from tastypie.exceptions import BadRequest, ImmediateHttpResponse, \
    UnsupportedFormat
import sqlalchemy
from reports.utils.sqlalchemy_bridge import Bridge
from reports.serializers import csv_convert

import logging

from reports import CSV_DELIMITER,LIST_DELIMITER_SQL_ARRAY,LIST_DELIMITER_URL_PARAM,\
    HTTP_PARAM_RAW_LISTS,HTTP_PARAM_USE_TITLES,HTTP_PARAM_USE_VOCAB,\
    LIST_BRACKETS, MAX_IMAGE_ROWS_PER_XLS_FILE, MAX_ROWS_PER_XLS_FILE
from tastypie.http import HttpNotFound
    


logger = logging.getLogger(__name__)


# FIXME: Test
class ChunkIterWrapper(object):
    ''' 
    Iterate in "chunks" of chunk_size chars.
    '''
    def __init__(self, iter_stream, chunk_size = 1024**2):
        self.iter_stream = iter_stream
        self.chunk_size = chunk_size
        self.fragment = None
    
    def next(self):
        logger.info(str(('chunking')))
        # Note: cStringIO will not accept unicode strings that can not be encoded as 7-bit ascii
        bytes = cStringIO.StringIO()
        # NOTE: StringIO will have difficulty decoding UTF-8 mixed with ascii
        #         bytes = StringIO.StringIO()
        
        bytecount = 0
        try:
            if self.fragment:
                logger.info(str(('fragment', self.fragment)))
                #  FIXME: if fragment is still > chunk_size, will just serve it here.
                bytecount = len(self.fragment)
                bytes.write(self.fragment)
                self.fragment = None

            while bytecount < self.chunk_size:
                row = self.iter_stream.next()
                rowlen = len(row)
                if bytecount + rowlen < self.chunk_size:
                    bytes.write(row)
                    bytecount += rowlen
                else:
                    nextbytes = self.chunk_size-bytecount
                    bytes.write(row[:nextbytes])
                    self.fragment = row[nextbytes:]
                    bytecount += nextbytes
        except StopIteration:
            if self.fragment:
                bytes.write(self.fragment)
                self.fragment = None
            else:
                raise StopIteration
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on ChunkIterWrapper next', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e   
        finally:
            if bytes.getvalue():
                logger.info(str(('chunk size', len(bytes.getvalue()), 'chars')))
                return bytes.getvalue()
            else:
                logger.info(str(('stop iteration', self.fragment, bytes.getvalue())))
                raise StopIteration

    def __iter__(self):
        return self

class ChunkIterWrapper1(object):
    ''' 
    Iterate in "chunks" of min_iteration size (final size is still dependent on
    length of original iteration chunk size.
    FIXME: break streamed output into exact size (buffer remainder for next iteration)
    FIXME: 'min_iteration' replaced with bytes
    '''
    def __init__(self, iter_stream, min_iteration = 1000):
        self.iter_stream = iter_stream
        self.min_iteration = min_iteration
    
    def next(self):
        logger.info(str(('chunking')))
        bytes = cStringIO.StringIO()
        try:
            for i in range(self.min_iteration):
                bytes.write(self.iter_stream.next())
        except StopIteration:
            raise StopIteration
        finally:
            if bytes.getvalue():
                logger.info(str(('chunk size', len(bytes.getvalue()), 'chars')))
                return bytes.getvalue()
            else:
                logger.info(str(('stop iteration', bytes.getvalue())))
                raise StopIteration

    def __iter__(self):
        return self

class SqlAlchemyResource(Resource):
    '''
    A resource that uses SqlAlchemy to facilitate the read queries:
    - get_list
    - get_detail
    Note: write operations not implemented
    
    FIXME: serialization of list fields: they are sent as concatenated strings with commas
    FIXME: Reagent/SMR specific
    '''
    
#     class Meta:
#         queryset = Reagent.objects.all()
#         serializer = LimsSerializer() # still have a serializer for error response

    
    def __init__(self, *args, **kwargs):
        # get a handle to the SqlAlchemy "bridge" for its table registry and 
        # connection handle
        self.bridge = Bridge()
        
        super(SqlAlchemyResource, self).__init__(*args, **kwargs)

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
        
        return _dict    
       
    
    def get_visible_fields(self, schema_fields, filter_fields, manual_field_includes,
                           is_for_detail=False):
        '''
        Construct an ordered dict of schema fields that are visible, based on
        - "list" in schema field["visibility"]
        - field key in (filter_fields or manual_field_includes)
        - field key in a schema field['dependencies'] 
        
        TODO: this method can be static
        TODO: this method is not SqlAlchemy specific
        '''
        logger.info(str(('get_visible_fields: field_hash initial: ', schema_fields.keys() )))
        try:
            visibility_setting = 'list'
            if is_for_detail:
                visibility_setting = 'detail'
            temp = { key:field for key,field in schema_fields.items() 
                if ((field.get('visibility', None) and visibility_setting in field['visibility']) 
                    or field['key'] in manual_field_includes
                    or '*' in manual_field_includes )}
            
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
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_visible_fields', 
                msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e 

    def build_sqlalchemy_columns(self, fields, base_query_tables=[], custom_columns={}):
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
        @param manual_includes - columns to include even if the field visibility is not set
        '''
        DEBUG_BUILD_COLUMNS = False or logger.isEnabledFor(logging.DEBUG)
        
        try:
            columns = OrderedDict()
            for field in fields:
                key = field['key']
                if key in custom_columns:
                    columns[key] = custom_columns[key]
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
                    logger.warn(str(('field is not in the base tables or in a linked field, and is not custom', key)))
            if DEBUG_BUILD_COLUMNS: logger.info(str(('columns', columns.keys())))
            return columns.values()
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on build sqlalchemy columns', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
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
                logger.warn(str(('order_by field not in visible fields, skipping: ', order_by )))
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
        logger.info(str(('filter value', filter_expr, value)))
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
                (search_expression, search_fields) = \
                    SqlAlchemyResource.build_sqlalchemy_filters_from_hash(schema,search_hash)
                search_expressions.append(search_expression)
                filter_fields.update(search_fields)

            if len(search_expressions) > 1:
                search_expressions = or_(*search_expressions)
            else:
                search_expressions = search_expressions[0]
            logger.info(str(('testing...', search_expressions)))
            if filter_expression is not None:
                filter_expression = and_(search_expressions,filter_expression)
            else: 
                filter_expression = search_expressions
            logger.info(str(('testing2...', filter_expression)))
                
        logger.info(str(('filter_expression', filter_expression, filter_fields)))
        
        return (filter_expression,filter_fields)
        
    
    @staticmethod
    def build_sqlalchemy_filters_from_hash(schema, param_hash):
        '''
        Attempt to create a SqlAlchemy whereclause out of django style filters:
        - field_name__filter_expression
        '''
        DEBUG_FILTERS = False or logger.isEnabledFor(logging.DEBUG)
        logger.info(str(('build_sqlalchemy_filters_from_hash', param_hash)))
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
                    logger.warn(str(('filter expression must be of the form "field_name__expression"',
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
                #     'exact', 'iexact', 'contains', 'icontains', 'gt', 'gte', 'lt', 'lte', 'in',
                #     'startswith', 'istartswith', 'endswith', 'iendswith', 'range', 'year',
                #     'month', 'day', 'week_day', 'hour', 'minute', 'second', 'isnull', 'search',
                #     'regex', 'iregex',
                # ])
                
                logger.info(str(('find filter', field_name, filter_type, value)))
                
    #             if not value:
    #                 continue
    #             
                expression = None
                if filter_type in ['exact','eq']:
                    expression = column(field_name)==value
                elif filter_type == 'about':
                    decimals = 0
                    if '.' in value:
                        decimals = len(value.split('.')[1])
                    expression = func.round(
                        sqlalchemy.sql.expression.cast(
                            column(field_name),sqlalchemy.types.Numeric),decimals) == value
                    logger.info(str(('create "about" expression for term:',filter_expr,
                        value,'decimals',decimals)))
                elif filter_type == 'contains':
                    expression = column(field_name).contains(value)
                elif filter_type == 'icontains':
                    expression = column(field_name).ilike('%{value}%'.format(value=value))
                elif filter_type == 'lt':
                    expression = column(field_name) < value
                elif filter_type == 'lte':
                    expression = column(field_name) <= value
                elif filter_type == 'gt':
                    expression = column(field_name) > value
                elif filter_type == 'gte':
                    expression = column(field_name) >= value
                elif filter_type == 'is_null':
                    logger.info(str(('is_null', field_name, value, str(value) )))
                    col = column(field_name)
                    if field['ui_type'] == 'string':
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
                        logger.error(str(('value for range expression must be list of length 2', 
                            field_name, filter_expr, value)))
                        continue
                    else:
                        expression = column(field_name).between(value[0],value[1],symmetric=True)
                else:
                    logger.error(str(('--- unknown filter type: ', field_name, filter_type,
                        'filter_expr',filter_expr )))
                    continue
    
                if inverted:
                    expression = not_(expression)
                
                logger.info(str(('filter_expr',filter_expr,'expression',str(expression) )))
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
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on build_sqlalchemy_filters_from_hash', 
                msg, exc_type, fname, exc_tb.tb_lineno)))
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
                            self._meta.resource_name + '.search requires a "search_data" param'}))
         
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
        Override the Tastypie/Django ORM get_list method - list/reporting operations will be 
        handled using SqlAlchemy
        '''
        raise NotImplemented(str(('get_list must be implemented for the SqlAlchemyResource', 
            self._meta.resource_name)) )

        
    def build_list_response(self,request, param_hash={}, **kwargs):

        raise NotImplemented(str(('get_list_response must be implemented for the SqlAlchemyResource', 
            self._meta.resource_name)) )
        
    def _cached_resultproxy(self, stmt, count_stmt, limit, offset):
        ''' ad-hoc cache for some resultsets
        NOTE: limit and offset are included because this version of sqlalchemy
        does not support printing of them with the select.compile() function.
        TODO: cache clearing
        '''
        
        from django.core.cache import cache
        
        conn = self.bridge.get_engine().connect()
        try:
            # use a hexdigest because statements can contain problematic chars 
            # for the memcache
            m = hashlib.md5()
            compiled_stmt = str(stmt.compile(compile_kwargs={"literal_binds": True}))
            key_digest = '%s_%s_%s' %(compiled_stmt, str(limit), str(offset))
            m.update(key_digest)
            key = m.hexdigest()
#             logger.info(str(('hash key:', key_digest, key)))
            cache_hit = cache.get(key)
            if cache_hit:
                if ('stmt' not in cache_hit or
                        cache_hit['stmt'] != compiled_stmt):
                    cache_hit = None
                    logger.warn(str(('cache collision for key', key, stmt)))
            
            if not cache_hit:
                # Note: if no cache hit, then retrive limit*n results, and 
                # cache several iterations at once.
                prefetch_number = 5
                new_limit = limit * prefetch_number 
                logger.info(str(('limit', limit, 'new limit for caching', new_limit)))
                stmt = stmt.limit(new_limit)
                logger.info('executing stmt')
                resultset = conn.execute(stmt)
                prefetched_result = [dict(row) for row in resultset] if resultset else []
                logger.info('executed stmt')
                
                logger.info(str(('no cache hit, execute count')))
                count = conn.execute(count_stmt).scalar()
                logger.info(str(('count', count)))
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
                        cache.set( key, _cache, None)
                        if y == 0:
                            cache_hit = _cache
                    else:
                        logger.info(str(('prefetched length: ', len(prefetched_result), _start, 'end')))
                        break
                logger.info(str(('cached iterations:', y)))
#                 cached_result = [dict(row) for row in resultset] if resultset else []
            else:
                logger.info(str(('cache hit')))   
                
            return cache_hit
 
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on conn execute', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e   
        
          
    def stream_response_from_cursor(self, request, stmt, count_stmt, 
            output_filename, field_hash={}, param_hash={}, 
            rowproxy_generator=None, is_for_detail=False,
            downloadID=None, title_function=None ):
        
        DEBUG_STREAMING = False or logger.isEnabledFor(logging.DEBUG)

        limit = param_hash.get('limit', request.GET.get('limit', self._meta.limit))        
        try:
            limit = int(limit)
        except ValueError:
            raise BadRequest(
                "Invalid limit '%s' provided. Please provide a positive integer." % limit)
        if limit > 0:    
            stmt = stmt.limit(limit)

        offset = param_hash.get('offset', request.GET.get('offset', 0) )
        try:
            offset = int(offset)
        except ValueError:
            raise BadRequest(
                "Invalid offset '%s' provided. Please provide a positive integer." % offset)
        if offset < 0:    
            offset = -offset
        stmt = stmt.offset(offset)


        list_brackets = LIST_BRACKETS
        if request.GET.get(HTTP_PARAM_RAW_LISTS, False):
            list_brackets = None

        conn = self.bridge.get_engine().connect()
        
        logger.info(str(('stmt', str(stmt))))
        if DEBUG_STREAMING:
            logger.info(str(('count stmt', str(count_stmt))))
        
        logger.info('excute stmt')
        result = conn.execute(stmt)
        logger.info('excuted stmt')
        
        logger.info(str(('rowproxy_generator', rowproxy_generator)))
        if rowproxy_generator:
            result = rowproxy_generator(result)
            # FIXME: test this for generators other than json generator

        
        desired_format = self.get_format(request)
        content_type=build_content_type(desired_format)
        logger.info(str(('desired_format', desired_format, 
            'content_type', content_type)))

        
        def interpolate_value_template(value_template, row):
            ''' 
            Utility class for transforming cell values:
            a "value_template" is of the form:
            "text... {field_name} ..text"
            wherein the {field_name} is replaced with the field-value
            '''                
            def get_value_from_template(matchobj):
                val = matchobj.group()
                val = re.sub(r'[{}]+', '', val)
                if row.has_key(val):
                    return row[val]
                else:
                    logger.error(str(('field needed for value template is not available', val, row)))
                    return ''
            return re.sub(r'{([^}]+)}', get_value_from_template, value_template)
        
        response = None
        
        if desired_format == 'application/json':
            
            if not is_for_detail:
                cache_hit = self._cached_resultproxy(stmt, count_stmt, limit, offset)
                if cache_hit:
                    result = cache_hit['cached_result']
                    count = cache_hit['count']
                    if rowproxy_generator:
                        result = rowproxy_generator(result)
                else:
                    # use result from before
                    logger.info(str(('execute count')))
                    count = conn.execute(count_stmt).scalar()
                logger.info(str(('====count====', count)))
                
                if count == 0:
                    raise ImmediateHttpResponse(
                        response=self.error_response(
                            request, {'empty result': 'no records found'},
                            response_class=HttpNotFound))
                
                meta = {
                    'limit': limit,
                    'offset': offset,
                    'total_count': count
                    }    
                        
            def json_generator(cursor):
                if DEBUG_STREAMING: logger.info(str(('meta', meta)))
                # NOTE, using "ensure_ascii" = True to force encoding of all 
                # chars to be encoded using ascii or escaped unicode; 
                # because some chars in db might be non-UTF8
                # and downstream programs have trouble with mixed encoding (cStringIO)
                if not is_for_detail:
                    yield '{ "meta": %s, "objects": [' % json.dumps(meta, ensure_ascii=True, encoding="utf-8")
                i=0
                for row in cursor:
                    if DEBUG_STREAMING: logger.info(str(('row', row, row.keys())))
                    if field_hash:
                        _dict = OrderedDict()
                        for key, field in field_hash.iteritems():
                            if DEBUG_STREAMING: 
                                logger.info(str(('key', key,  row.has_key(key))))
                            value = None
                            try:
                                if row.has_key(key):
                                    value = row[key]
                                if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                                     or field.get('linked_field_type',None) == 'fields.ListField'
                                     or field.get('ui_type', None) == 'list' ):
                                    # FIXME: need to do an escaped split
                                    if DEBUG_STREAMING: logger.info(str(('split', key, value)))
                                    value = value.split(LIST_DELIMITER_SQL_ARRAY)
    
                                if DEBUG_STREAMING: logger.info(str(('key val', key, value, field)))
                                _dict[key] = value
                                
                                if field.get('value_template', None):
                                    value_template = field['value_template']
                                    if DEBUG_STREAMING: logger.info(str(('field', key, value_template)))
                                    newval = interpolate_value_template(value_template, row)
                                    if field['ui_type'] == 'image':
                                        # hack to speed things up:
                                        if ( key == 'structure_image' and
                                                row.has_key('library_well_type') and
                                                row['library_well_type'] == 'empty' ):
                                            continue
                                        # see if the specified url is available
                                        try:
                                            view, args, kwargs = resolve(newval)
                                            kwargs['request'] = request
                                            view(*args, **kwargs)
                                            _dict[key] = newval
                                        except Exception, e:
                                            logger.info(str(('no image at', newval,e)))
                                    else:
                                        _dict[key]=newval

                            except Exception, e:
                                exc_type, exc_obj, exc_tb = sys.exc_info()
                                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
                                msg = str(e)
                                logger.warn(str(('on get_vocabularies_by_scope', 
                                    msg, exc_type, fname, exc_tb.tb_lineno)))
                                raise e                      
                    else:
                        if DEBUG_STREAMING: logger.info(str(('raw', row)))
                        _dict = dict((x,y) for x, y in row.items())

                    i += 1
                    try:
                        if DEBUG_STREAMING: logger.info(str(('_dict',i, _dict)))
                        if i == 1:
                            # NOTE, using "ensure_ascii" = True to force encoding of all 
                            # chars to the ascii charset; otherwise, cStringIO has problems
                            yield json.dumps(_dict, cls=DjangoJSONEncoder,
                                sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
                        else:
                            # NOTE, using "ensure_ascii" = True to force encoding of all 
                            # chars to the ascii charset; otherwise, cStringIO has problems
                            # so "tm" becomes \u2122
                            # Upon fp.write  the unicode is converted to the default charset?
                            # also, CStringIO doesn't support UTF-8 mixed with ascii, for instance
                            yield ', ' + json.dumps(_dict, cls=DjangoJSONEncoder,
                                sort_keys=True, ensure_ascii=True, indent=2, encoding="utf-8")
                    except Exception, e:
                        print 'Exception'
                        logger.info(str(('exception')))
                        logger.error(str(('ex', _dict, e)))
                logger.info('streaming finished')
                
                if not is_for_detail:
                    yield ' ] }'
            
            if DEBUG_STREAMING: logger.info(str(('ready to stream 1...')))
            response = StreamingHttpResponse(ChunkIterWrapper(json_generator(result)))
            response['Content-Type'] = content_type
        
        
        elif( desired_format == 'application/xls' or
            desired_format == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'): 
            # NOTE: will respond with xlsx for both desired formats
            
            # FIXME: how to abstract images?
            
            structure_image_dir = os.path.abspath(settings.WELL_STRUCTURE_IMAGE_DIR)
            # create a temp dir
            # FIXME: temp directory as a setting
            temp_dir = os.path.join('/tmp', str(time.clock()).replace('.', '_'))
            os.mkdir(temp_dir)
            try:
                from xlsxwriter.workbook import Workbook
                irow=0
                # Create an new Excel file and add a worksheet.
                filename = '%s.xlsx' % (output_filename)
                temp_file = os.path.join(temp_dir, filename)
                if DEBUG_STREAMING: logger.info(str(('temp file', temp_file)))
                workbook = Workbook(temp_file)
                worksheet = workbook.add_worksheet()
                
                # FIXME: only need max rows if the file will be too big (structure images)
                # or too long (>65535, for some versions of xls; for that case
                # should implement a mult-sheet solution.
                max_rows_per_file = MAX_ROWS_PER_XLS_FILE
                file_names_to_zip = [temp_file]
                filerow = 0
                for row in result:
                    if filerow == 0:
                        if field_hash:
                            for col, (key, field) in enumerate(field_hash.items()):
                                if title_function:
                                    key = title_function(key) 
                                worksheet.write(filerow,col,key)
                                ## TODO: option to write titles
                                # worksheet.write(filerow,col,field['title'])
                        else:
                            for col,name in enumerate(row.keys()):
                                worksheet.write(filerow,col,name)
                        filerow += 1
                    
                    if field_hash:
                        for col, (key, field) in enumerate(field_hash.iteritems()):
                            value = None
                            if row.has_key(key):
                                value = row[key]
                            
                            if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                                 or field.get('linked_field_type',None) == 'fields.ListField'
                                 or field.get('ui_type', None) == 'list' ):
                                value = value.split(LIST_DELIMITER_SQL_ARRAY)
                            worksheet.write(filerow, col, 
                                csv_convert(value, delimiter=LIST_DELIMITER_XLS, 
                                    list_brackets=list_brackets))

                            if field.get('value_template', None):
                                value_template = field['value_template']
                                if DEBUG_STREAMING: logger.info(str(('field', key, value_template)))
                                newval = interpolate_value_template(value_template, row)
                                if field['ui_type'] == 'image':
                                    max_rows_per_file = MAX_IMAGE_ROWS_PER_XLS_FILE
                                    # hack to speed things up:
                                    if ( key == 'structure_image' and
                                            row.has_key('library_well_type') and
                                            row['library_well_type'] == 'empty' ):
                                        continue
                                    # see if the specified url is available
                                    try:
                                        view, args, kwargs = resolve(newval)
                                        kwargs['request'] = request
                                        response = view(*args, **kwargs)
                                        image = Image.open(io.BytesIO(response.content))
                                        height = image.size[1]
                                        width = image.size[0]
                                        worksheet.set_row(filerow, height)
                                        scaling = 0.130 # trial and error width in default excel font
                                        worksheet.set_column(col,col, width*scaling)
                                        worksheet.insert_image(
                                            filerow, col, newval, 
                                            {'image_data': io.BytesIO(response.content)})
                                    except Exception, e:
                                        logger.info(str(('no image at', newval,e)))
                                else:
                                    worksheet.write(filerow, col, newval)

                    else:
                        
                        for col,val in enumerate(row.values()):
                            worksheet.write(filerow, col, val)

                    irow +=1
                    filerow +=1
                    
                    if irow % max_rows_per_file == 0:
                        workbook.close()
                        logger.info(str(('wrote file', temp_file)))

                        # Create an new Excel file and add a worksheet.
                        filename = '%s_%s.xlsx' % (output_filename, irow)
                        temp_file = os.path.join(temp_dir, filename)
                        logger.info(str(('temp file', temp_file)))
                        workbook = Workbook(temp_file)
                        worksheet = workbook.add_worksheet()
                        
                        file_names_to_zip.append(temp_file)
                        filerow = 0
                    
                workbook.close()
                logger.info(str(('wrote file', temp_file)))

                content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet; charset=utf-8'
                if len(file_names_to_zip) >1:
                    # create a temp zip file
                    content_type='application/zip; charset=utf-8'
                    temp_file = os.path.join('/tmp',str(time.clock()))
                    logger.info(str(('temp ZIP file', temp_file)))

                    with ZipFile(temp_file, 'w') as zip_file:
                        for _file in file_names_to_zip:
                            zip_file.write(_file, os.path.basename(_file))
                    logger.info(str(('wrote file', temp_file)))
                    filename = '%s.zip' % output_filename

                _file = file(temp_file)
                logger.info(str(('download tmp file',temp_file,_file)))
                wrapper = FileWrapper(_file)
                response = StreamingHttpResponse(
                    wrapper, content_type=content_type) 
                response['Content-Disposition'] = \
                    'attachment; filename=%s' % filename
                response['Content-Length'] = os.path.getsize(temp_file)
            except Exception, e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
                msg = str(e)
                logger.warn(str(('on xlsx & zip file process', 
                    self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
                raise e   
            finally:
                try:
                    logger.info(str(('rmdir', temp_dir)))
                    shutil.rmtree(temp_dir)
                    if os.path.exists(temp_file):
                        logger.info(str(('remove', temp_file)))
                        os.remove(temp_file)     
                    logger.info(str(('removed', temp_dir)))
                except Exception, e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
                    msg = str(e)
                    logger.warn(str(('on xlsx & zip file process', self._meta.resource_name, 
                        msg, exc_type, fname, exc_tb.tb_lineno)))
                    raise ImmediateHttpResponse(str(('ex during rmdir', e)))
        elif desired_format == 'chemical/x-mdl-sdfile':
            import reports.utils.sdf2py
            MOLDATAKEY = reports.utils.sdf2py.MOLDATAKEY
            try:
                def sdf_generator(cursor):
                    i = 0
                    for row in cursor:
                        i += 1

                        if row.has_key(MOLDATAKEY) and row[MOLDATAKEY]:
                            yield str(row[MOLDATAKEY])
                            yield '\n' 

                        if field_hash:
                            for col, (key, field) in enumerate(field_hash.iteritems()):
                                if key == MOLDATAKEY: 
                                    continue
                                yield '> <%s>\n' % key
                                # according to 
                                # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
                                # "only one blank line should terminate a data item"
                                value = None
                                
                                if row.has_key(key):
                                    value = row[key]
                                
                                if field.get('value_template', None):
                                    value_template = field['value_template']
                                    if DEBUG_STREAMING:     
                                        logger.info(str(('field', key, value_template)))
                                    value = interpolate_value_template(value_template, row)
                                if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                                     or field.get('linked_field_type',None) == 'fields.ListField'
                                     or field.get('ui_type', None) == 'list' ):
                                    value = value.split(LIST_DELIMITER_SQL_ARRAY)

                                if value:
                                    # find lists, but not strings (or dicts)
                                    # Note: a dict here will be non-standard; probably an error 
                                    # report, so just stringify dicts as is.
                                    if not hasattr(value, "strip") and isinstance(value, (list,tuple)): 
                                        for x in value:
                                            # DB should be UTF-8, so this should not be necessary,
                                            # however, it appears we have legacy non-utf data in 
                                            # some tables (i.e. small_molecule_compound_name 193090
                                            yield unicode.encode(x,'utf-8')
                                            yield '\n'
                                    else:
                                        yield str(value)
                                        yield '\n'
            
                                yield '\n'
                            yield '$$$$\n'
                                    
                        else:
                        
                            for k,v in row.items():
                                if k == MOLDATAKEY: 
                                    continue
                                yield '> <%s>\n' % k
                                # according to 
                                # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
                                # "only one blank line should terminate a data item"
                                if v:
                                    # find lists, but not strings (or dicts)
                                    # Note: a dict here will be non-standard; probably an error 
                                    # report, so just stringify dicts as is.
                                    if not hasattr(v, "strip") and isinstance(v, (list,tuple)): 
                                        for x in v:
                                            # DB should be UTF-8, so this should not be necessary,
                                            # however, it appears we have legacy non-utf data in 
                                            # some tables (i.e. small_molecule_compound_name 193090
                                            yield unicode.encode(x,'utf-8')
                                            yield '\n'
                                    else:
                                        yield str(v)
                                        yield '\n'
                
                                yield '\n'
                            yield '$$$$\n'
        
                    logger.info(str(('wrote i', i)))    
    
                response = StreamingHttpResponse(
                    ChunkIterWrapper(sdf_generator(result)),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.sdf' % output_filename
            except Exception, e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
                msg = str(e)
                logger.warn(str(('sdf export process', self._meta.resource_name, 
                    msg, exc_type, fname, exc_tb.tb_lineno)))
                raise ImmediateHttpResponse(str(('ex during sdf export', e)))
                
        
        elif desired_format == 'text/csv':
            class Echo(object):
                """An object that implements just the write method of the file-like
                interface.
                """
                def write(self, value):
                    return value

            pseudo_buffer = Echo()
            csvwriter = csv.writer(
                pseudo_buffer, delimiter=CSV_DELIMITER, quotechar='"', 
                quoting=csv.QUOTE_ALL, lineterminator="\n")
            def csv_generator(cursor):
                i=0
                for row in cursor:
                    i += 1
                    if i == 1:
                        if field_hash:
                            titles = field_hash.keys()
                            logger.info(str(('keys', titles, title_function)))
                            if title_function:
                                titles = [title_function(key) for key in titles]
                            logger.info(str(('titles', titles)))
                            yield csvwriter.writerow(titles)
                            # TODO: option to write titles:
                            # yield csvwriter.writerow([field['title'] for field in field_hash.values()])
                        else:
                            yield csvwriter.writerow(row.keys())

                    if field_hash:
                        values = []
                        for col, (key, field) in enumerate(field_hash.iteritems()):
                            value = ''
                            if row.has_key(key):
                                value = row[key]
                            if value and ( field.get('json_field_type',None) == 'fields.ListField' 
                                 or field.get('linked_field_type',None) == 'fields.ListField'
                                 or field.get('ui_type', None) == 'list' ):
                                # FIXME: must quote special strings?
#                                 value = '[' + ",".join(value.split(LIST_DELIMITER_SQL_ARRAY)) + ']'
                                value = value.split(LIST_DELIMITER_SQL_ARRAY)
                            
#                             if field.get('data_type', None):
#                                 data_type = field['data_type']
#                                 if data_type == "float":
#                                     value = 
                                
                            if field.get('value_template', None):
                                value_template = field['value_template']
                                newval = interpolate_value_template(value_template, row)
                                if field['ui_type'] == 'image':
                                    # hack to speed things up:
                                    if ( key == 'structure_image' and
                                            row.has_key('library_well_type') and
                                            row['library_well_type'] == 'empty' ):
                                        continue
                                    # see if the specified url is available
                                    try:
                                        view, args, kwargs = resolve(newval)
                                        kwargs['request'] = request
                                        response = view(*args, **kwargs)
                                        value = newval
                                    except Exception, e:
                                        logger.info(str(('no image at', newval,e)))
                                else:
                                    value = newval

                            values.append(value)
                        yield csvwriter.writerow([
                            csv_convert(val, list_brackets=list_brackets) 
                                for val in values])
                    
                    else:
                        yield csvwriter.writerow([
                            csv_convert(val, list_brackets=list_brackets) 
                                for val in row.values()])

            response = StreamingHttpResponse(csv_generator(result),
                content_type=content_type)
            name = self._meta.resource_name
            response['Content-Disposition'] = \
                'attachment; filename=%s.csv' % output_filename
        else:
            msg = str(('unknown format', desired_format, output_filename))
            logger.error(msg)
            raise ImmediateHttpResponse(msg)

#         downloadID = request.GET.get('downloadID', None)
#         if downloadID:
#             logger.info(str(('set cookie','downloadID', downloadID )))
#             response.set_cookie('downloadID', downloadID)
#         else:
#             logger.info(str(('no downloadID', request.GET )))

        return response

