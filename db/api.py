# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import cStringIO
from collections import defaultdict, OrderedDict
from copy import deepcopy
from decimal import Decimal
import hashlib
import io
import json
import logging
from operator import itemgetter
import os
import random
import re
from tempfile import SpooledTemporaryFile, NamedTemporaryFile
import time
import urllib
from zipfile import ZipFile, ZipInfo

from aldjemy.core import get_engine, get_tables
import aldjemy.core
from django.conf import settings
from django.conf.urls import url
from django.core.cache import cache, caches
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.core.serializers.json import DjangoJSONEncoder
from django.db import connection
from django.db import transaction
from django.db.models import F, Q
from django.db.utils import ProgrammingError
from django.forms.models import model_to_dict
from django.http import Http404
from django.http.request import HttpRequest
from django.http.response import StreamingHttpResponse, HttpResponse
from django.utils import timezone
import six
from sqlalchemy import select, text, case
import sqlalchemy
from sqlalchemy.dialects import postgresql
from sqlalchemy.dialects.postgresql import array
from sqlalchemy.sql import and_, or_, not_, asc, desc, alias, Alias, func
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join, insert, delete, distinct, \
    exists, cast, union_all, union
from sqlalchemy.sql.expression import nullslast
import sqlalchemy.sql.schema
import sqlalchemy.sql.sqltypes
from tastypie.authentication import BasicAuthentication, MultiAuthentication
from tastypie.exceptions import BadRequest
from tastypie.resources import convert_post_to_put
from tastypie.utils.urls import trailing_slash
import unicodecsv

from db import WELL_ID_PATTERN, WELL_NAME_PATTERN, PLATE_PATTERN, \
    PLATE_RANGE_PATTERN, COPY_NAME_PATTERN
from db.models import ScreensaverUser, Screen, \
    ScreenResult, DataColumn, Library, Plate, Copy, \
    CopyWell, UserFacilityUsageRole, \
    PlateLocation, Reagent, Well, Activity, \
    AdministrativeActivity, SmallMoleculeReagent, SilencingReagent, GeneSymbol, \
    NaturalProductReagent, Molfile, Gene, GeneGenbankAccessionNumber, \
    CherryPickRequest, CherryPickAssayPlate, CherryPickLiquidTransfer, \
    CachedQuery, UserChecklist, AttachedFile, \
    ServiceActivity, LabActivity, Screening, LibraryScreening, AssayPlate, \
    SmallMoleculeChembankId, SmallMoleculePubchemCid, SmallMoleculeChemblId, \
    SmallMoleculeCompoundName, ScreenCellLines, ScreenFundingSupports, \
    ScreenKeyword, ResultValue, AssayWell, Publication, ScreenerCherryPick, \
    LabCherryPick, CherryPickRequestEmptyWell, CherryPickScreening, \
    LabAffiliation, UserAgreement, RawDataTransform, RawDataInputFile
from db.support import lims_utils, screen_result_importer, bin_packer,\
    raw_data_reader, plate_matrix_transformer
from db.support.data_converter import default_converter
from db.support.screen_result_importer import PARTITION_POSITIVE_MAPPING, \
    CONFIRMED_POSITIVE_MAPPING
from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    HTTP_PARAM_USE_TITLES, HTTP_PARAM_USE_VOCAB, HTTP_PARAM_DATA_INTERCHANGE, \
    LIST_BRACKETS, HTTP_PARAM_RAW_LISTS, HEADER_APILOG_COMMENT, ValidationError, \
    LIST_DELIMITER_SUB_ARRAY
from reports import ValidationError, InformationError, _now
from reports.api import API_MSG_COMMENTS, API_MSG_CREATED, \
    API_MSG_SUBMIT_COUNT, API_MSG_UNCHANGED, API_MSG_UPDATED, \
    API_MSG_ACTION, API_MSG_RESULT, API_MSG_WARNING, API_RESULT_DATA, \
    API_RESULT_META, API_RESULT_OBJ, API_MSG_NOT_ALLOWED, API_PARAM_OVERRIDE, \
    DEBUG_AUTHORIZATION, FieldResource
from reports.api import ApiLogResource, UserGroupAuthorization, \
    VocabularyResource, UserResource, UserGroupResource, ApiLogResource, \
    write_authorization, read_authorization, UserResourceAuthorization
import reports.api
from reports.api_base import un_cache, IccblSessionAuthentication
from reports.models import Vocabulary, ApiLog, UserProfile, \
    API_ACTION_DELETE, API_ACTION_PUT, API_ACTION_PATCH, API_ACTION_CREATE
from reports.serialize import parse_val, XLSX_MIMETYPE, SDF_MIMETYPE, \
    XLS_MIMETYPE, JSON_MIMETYPE, CSV_MIMETYPE, ZIP_MIMETYPE
from reports.serialize.csvutils import convert_list_vals
from reports.serialize.streaming_serializers import ChunkIterWrapper, \
    json_generator, cursor_generator, sdf_generator, generic_xlsx_response, \
    csv_generator, get_xls_response, image_generator, closing_iterator_wrapper, \
    FileWrapper1
from reports.serializers import LimsSerializer, \
    XLSSerializer, ScreenResultSerializer
from reports.sqlalchemy_resource import SqlAlchemyResource
from reports.sqlalchemy_resource import _concat, _concat_with_sep
from db.support.plate_matrix_transformer import Collation
import xlsxwriter

PLATE_NUMBER_SQL_FORMAT = 'FM9900000'
PSYCOPG_NULL = '\\N'
MAX_SPOOLFILE_SIZE = 1024*1024

##### API CONSTANTS
# TODO: move to an API-accessible properties file

API_MSG_COPYWELLS_DEALLOCATED = 'Copy wells deallocated'
API_MSG_COPYWELLS_ALLOCATED = 'Copy wells allocated'
API_MSG_SCREENING_PLATES_UPDATED = 'Library Plates updated'
API_MSG_SCREENING_EXTANT_PLATE_COUNT = 'Extant plates'
API_MSG_SCREENING_ADDED_PLATE_COUNT = 'Added plates'
API_MSG_SCREENING_DELETED_PLATE_COUNT = 'Deleted plates'
API_MSG_SCREENING_TOTAL_PLATE_COUNT = 'Library Plate Count'
API_MSG_SCP_CREATED = 'Created'
API_MSG_SCP_UNSELECTED = 'Unselected'
API_MSG_SCP_RESELECTED = 'Reselected'
API_MSG_SCPS_DELETED = 'Screener Cherry Picks Removed'      
API_MSG_SCPS_CREATED = 'Screener Cherry Picks Created'
API_MSG_LCP_CHANGED = 'Changed'
API_MSG_LCP_DESELECTED = 'Deselected'
API_MSG_LCP_SELECTED = 'Selected'
API_MSG_LCP_MULTIPLE_SELECTIONS_SUBMITTED = 'Multiple lab cherry pick selections submitted'
API_MSG_LCPS_CREATED = 'Lab Cherry Picks Created'
API_MSG_LCPS_MUST_BE_DELETED = "Lab Cherry Picks must be deleted"
API_MSG_LCPS_ASSIGNED = 'Assigned to Copies'
API_MSG_LCPS_UNFULFILLED = 'Unfulfilled'
API_MSG_LCPS_INSUFFICIENT_VOLUME = 'Insufficient volume'
API_MSG_LCPS_VOLUME_OVERRIDDEN = 'Insufficient volume overridden'
API_MSG_LCPS_REMOVED = 'Lab Cherry Picks Removed'
API_MSG_LCP_PLATES_ASSIGNED = 'Copy Plate assigned'
API_MSG_LCP_SOURCE_PLATES_ALLOCATED = 'Source plates having wells allocated'
API_MSG_LCP_SOURCE_PLATES_DEALLOCATED = 'Source plates having wells deallocated'
API_MSG_LCP_ASSAY_PLATES_PLATED = 'Cherry pick assay plates (plated)'
API_MSG_LCP_ASSAY_PLATES_CREATED = 'Cherry pick assay plates (created)'
API_MSG_CPR_ASSAY_PLATES_REMOVED = 'Cherry pick assay plates (removed)'
API_MSG_PLATING_CANCELED = 'Plating canceled'
API_MSG_CPR_PLATES_PLATED = 'Plates plated'
API_MSG_CPR_PLATES_SCREENED = 'Plates screened'
API_MSG_CPR_PLATED_CANCEL_DISALLOWED = 'Plating reservation may not be canceled after plates are plated'
VOCAB_LCP_STATUS_SELECTED = 'selected'
VOCAB_LCP_STATUS_NOT_SELECTED = 'not_selected'
VOCAB_LCP_STATUS_UNFULFILLED = 'unfulfilled'
VOCAB_LCP_STATUS_PLATED = 'plated'
VOCAB_USER_CLASSIFICATION_PI = 'principal_investigator'
VOCAB_SCREEN_TYPE_SM = 'small_molecule'
VOCAB_SCREEN_TYPE_RNAI = 'rnai'
# API_PARAM_PLATE_MAPPING_OVERRIDE = 'plate_mapping_override'
API_PARAM_SHOW_OTHER_REAGENTS = 'show_other_reagents'
API_PARAM_SHOW_ALTERNATE_SELECTIONS = 'show_alternate_selections'
API_PARAM_SHOW_COPY_WELLS = 'show_copy_wells'
API_PARAM_SHOW_RETIRED_COPY_WELlS = 'show_available_and_retired_copy_wells'
API_PARAM_SHOW_UNFULFILLED = 'show_unfulfilled'
API_PARAM_SHOW_INSUFFICIENT = 'show_insufficient'
API_PARAM_SHOW_MANUAL = 'show_manual'
API_PARAM_VOLUME_OVERRIDE = 'volume_override'
API_PARAM_SET_DESELECTED_TO_ZERO = 'set_deselected_to_zero'
logger = logging.getLogger(__name__)

DEBUG_SCREEN_ACCESS = False or logger.isEnabledFor(logging.DEBUG)
DEBUG_DC_ACCESS = False or logger.isEnabledFor(logging.DEBUG)
DEBUG_RV_CREATE = False or logger.isEnabledFor(logging.DEBUG)
    
def _get_raw_time_string():
  return timezone.now().strftime("%Y%m%d%H%M%S")

class DbApiResource(reports.api.ApiResource):

    def __init__(self, **kwargs):
        super(DbApiResource,self).__init__(**kwargs)
        self.resource_resource = None
        self.apilog_resource = None
        self.field_resource = None

        # Create index tables (shared by resources)
        with get_engine().connect() as conn:
            self._create_well_query_index_table(conn)
            self._create_well_data_column_positive_index_table(conn)
            self._create_screen_overlap_table(conn)
    
    @classmethod
    def get_table_def(cls,table_name):
        '''
        Work around to create table definitions for dynamic tables
        '''
        if table_name in get_tables():
            return get_tables()[table_name]
        else:
            if table_name == 'well_data_column_positive_index':
                return sqlalchemy.sql.schema.Table(
                    'well_data_column_positive_index', 
                    aldjemy.core.get_meta(),
                    sqlalchemy.Column('well_id', sqlalchemy.String),
                    sqlalchemy.Column('data_column_id', sqlalchemy.Integer),
                    sqlalchemy.Column('screen_id', sqlalchemy.Integer),
                    )
            elif table_name == 'screen_overlap':
                return sqlalchemy.sql.schema.Table(
                    'screen_overlap', 
                    aldjemy.core.get_meta(),
                    sqlalchemy.Column('screen_id', sqlalchemy.Integer),
                    sqlalchemy.Column('overlap_screen_id', sqlalchemy.Integer))
                
            else:
                raise Exception('unknown table: %r', table_name)
        

    def _create_well_query_index_table(self, conn):
        try:
            conn.execute(text('select * from "well_query_index"  limit 1; '));
            logger.debug('The well_query_index table exists')
            return
        except Exception as e:
            logger.exception('creating the well_query_index table')
       
        try:
            conn.execute(text(
                'CREATE TABLE "well_query_index" ('
                '"id" serial NOT NULL,  '
                '"well_id" text NOT NULL REFERENCES "well" ("well_id") '
                '    DEFERRABLE INITIALLY DEFERRED,'
                '"query_id" integer NOT NULL REFERENCES "cached_query" ("id") '
                '    DEFERRABLE INITIALLY DEFERRED'
                ');'
            ));
            logger.info('the well_query_index table has been created')
        except Exception, e:
            logger.info((
                'Exception: %r on trying to create the well_query_index table,'
                'note that this is normal if the table already exists '
                '(PostgreSQL <9.1 has no "CREATE TABLE IF NOT EXISTS"'), e)

    def _create_well_data_column_positive_index_table(self, conn):
        try:
            conn.execute(text(
                'select * from "well_data_column_positive_index" limit 1; '));
            logger.debug('The well_data_column_positive_index table exists')
            return
        except Exception as e:
            logger.info('creating the well_data_column_positive_index table')
        
        try:
            conn.execute(text(
                'CREATE TABLE well_data_column_positive_index ('
                ' "well_id" text NOT NULL, '
                ' "data_column_id" integer NOT NULL, ' 
                ' "screen_id" integer NOT NULL ' 
                ');'
            ));
            conn.execute(text(
                'CREATE INDEX wdc_screen_id '
                'on well_data_column_positive_index(screen_id);'))
            # TODO: determine if indexes would help.
            # Note: foreign keys are not needed and complicate delete
            # 'CREATE TABLE well_data_column_positive_index ('
            # ' "well_id" text NOT NULL REFERENCES "well" ("well_id") '
            # '    DEFERRABLE INITIALLY DEFERRED,'
            # ' "data_column_id" integer NOT NULL ' 
            # ' REFERENCES "data_column" ("data_column_id") '
            # '    DEFERRABLE INITIALLY DEFERRED, '
            # ' "screen_id" integer NOT NULL ' 
            # ' REFERENCES "screen" ("screen_id") '
            # '    DEFERRABLE INITIALLY DEFERRED'
            
            logger.info('the well_data_column_positive_index table created')
        except Exception, e:
            logger.info((
                'Exception: %r on trying to create the '
                'well_data_column_positive_index table,'
                'note that this is normal if the table already exists '
                '(PostgreSQL <9.1 has no "CREATE TABLE IF NOT EXISTS"'), e)
    
    def _create_screen_overlap_table(self, conn):
        try:
            conn.execute(text(
                'select * from "screen_overlap" limit 1; '));
            logger.debug('The screen_overlap table exists')
            return
        except Exception as e:
            logger.info('creating the screen_overlap table')
        
        try:
            conn.execute(text(
                'CREATE TABLE screen_overlap ('
                ' "screen_id" integer NOT NULL, '
                ' "overlap_screen_id" integer NOT NULL '
                ');'
            ));
            # Note: foreign keys are not needed and complicate delete
            # 'CREATE TABLE screen_overlap ('
            # ' "screen_id" integer NOT NULL '
            # 'REFERENCES "screen" ("screen_id") '
            # '    DEFERRABLE INITIALLY DEFERRED,'
            # ' "overlap_screen_id" integer NOT NULL '
            # 'REFERENCES "screen" ("screen_id") '
            # '    DEFERRABLE INITIALLY DEFERRED '

            conn.execute(text(
                'CREATE INDEX screen_overlap_screen_id '
                'on screen_overlap(screen_id);'))
            conn.execute(text(
                'CREATE INDEX screen_overlap_overlap_screen_id '
                'on screen_overlap(overlap_screen_id);'))
            logger.info('the screen_overlap table created')
        except Exception, e:
            logger.info((
                'Exception: %r on trying to create the '
                'screen_overlap table,'
                'note that this is normal if the table already exists '
                '(PostgreSQL <9.1 has no "CREATE TABLE IF NOT EXISTS"'), e)

    @classmethod
    def get_create_well_data_column_positive_index(cls):
        '''
        If the well_data_column_positive_index has been cleared, recreate.
        - recreate the entire index each time
        '''
        bridge = get_tables()
        _aw = bridge['assay_well']
        _sr = bridge['screen_result']
        _s = bridge['screen']
        _dc = bridge['data_column']
        _wdc = cls.get_table_def('well_data_column_positive_index')

        # Recreate the well_data_column_positive_index:
        # Note: because this is lazy recreated, we are creating the whole 
        # index each time; could just recreate for the specific screen
        count_stmt = select([func.count()]).select_from(_wdc)
        count = 0
        with get_engine().begin() as conn:
            count = int(conn.execute(count_stmt).scalar())        
            logger.info('well_data_column_positive_index count: %r', count)
        
        if count == 0:
            logger.info('well_data_column_positive_index count: %r, recreating', count)
        
            base_stmt = join(
                _aw, _dc, _aw.c.screen_result_id == _dc.c.screen_result_id)
            base_stmt = base_stmt.join(_sr, _sr.c.screen_result_id==_aw.c.screen_result_id)
            base_stmt = base_stmt.join(_s, _sr.c.screen_id==_s.c.screen_id)
            base_stmt = select([
                _aw.c.well_id,
                _dc.c.data_column_id,
                _sr.c.screen_id
                ]).select_from(base_stmt)            
            base_stmt = base_stmt.where(_aw.c.is_positive)
            base_stmt = base_stmt.where(_s.c.study_type==None)
            base_stmt = base_stmt.where(
                _dc.c.data_type.in_([
                    'boolean_positive_indicator',
                    'partition_positive_indicator', 
                    'confirmed_positive_indicator']))
            base_stmt = base_stmt.order_by(
                _dc.c.data_column_id, _aw.c.well_id)
            insert_statement = (
                insert(_wdc)
                    .from_select(['well_id', 'data_column_id','screen_id'], base_stmt))
            logger.info('execute mutual pos insert statement...')
            logger.debug(
                'mutual pos insert statement: %r',
                str(insert_statement.compile(
                    dialect=postgresql.dialect(),
                    compile_kwargs={"literal_binds": True})))
            get_engine().execute(insert_statement)
            logger.info('mutual pos insert statement, executed.')
        return _wdc
    
    @classmethod
    def get_create_screen_overlap_indexes(cls):
        logger.info('get_create_screen_overlap_indexes...')
        _wdc = cls.get_create_well_data_column_positive_index()
        _screen_overlap = cls.get_table_def('screen_overlap')
        wdc1 = _wdc.alias('wdc1')
        wdc2 = _wdc.alias('wdc2')
        
        count_stmt = select([func.count()]).select_from(_screen_overlap)
        count = 0
        with get_engine().begin() as conn:
            count = int(conn.execute(count_stmt).scalar())        
            logger.info('screen_overlap count: %r', count)
        
        if count == 0:
            logger.info('screen_overlap count:: %r, recreating', count)
        
            base_stmt = (
                select([wdc1.c.screen_id, wdc2.c.screen_id.label('overlap_screen_id')])
                .select_from(wdc1)
                .select_from(wdc2)
                .where(wdc1.c.screen_id!=wdc2.c.screen_id)
                .where(wdc1.c.well_id==wdc2.c.well_id)
                .group_by(wdc1.c.screen_id,wdc2.c.screen_id))
            insert_statement = (
                insert(_screen_overlap)
                    .from_select(['screen_id','overlap_screen_id'], base_stmt))
            logger.debug(
                'screen_overlap insert statement: %r',
                str(insert_statement.compile(
                    dialect=postgresql.dialect(),
                    compile_kwargs={"literal_binds": True})))
            get_engine().execute(insert_statement)
            logger.info('screen_overlap insert statement, executed.')
        return _screen_overlap
    
    def get_apilog_resource(self):
        if self.apilog_resource is None:
            self.apilog_resource = ApiLogResource()
        return self.apilog_resource
    
    def get_resource_resource(self):
        if self.resource_resource is None:
            self.resource_resource = ResourceResource()
        return self.resource_resource
    
    def get_field_resource(self):
        if self.field_resource is None:
            self.field_resource = FieldResource()
        return self.field_resource
    
    def get_request_user(self,request):
        
        return ScreensaverUser.objects.get(username=request.user.username)
    
    def build_schema(self, user=None):
        logger.debug('build schema for: %r', self._meta.resource_name)
        return self.get_resource_resource()._get_resource_schema(
            self._meta.resource_name, user=user)
    
class PlateLocationResource(DbApiResource):        
    
    class Meta:
        queryset = PlateLocation.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'platelocation'
        authorization = UserGroupAuthorization(resource_name)
        serializer = LimsSerializer()
        always_return_data = True 
        
    def __init__(self, **kwargs):
        super(PlateLocationResource, self).__init__(**kwargs)
        
        self.lcp_resource = None
        
    def get_librarycopyplate_resource(self):
        if not self.lcp_resource:
            self.lcp_resource = LibraryCopyPlateResource()
        return self.lcp_resource
    
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<room>[\w \-<>]+)"
                r"/(?P<freezer>[\w \-<>]+)"
                r"/(?P<shelf>[\w \-<>]+)"
                r"/(?P<bin>[\w <>]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    @read_authorization
    def get_detail(self, request, **kwargs):
        
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True

        kwargs['room__eq'] = kwargs.pop('room', None)
        kwargs['freezer__eq'] = kwargs.pop('freezer', None)
        kwargs['shelf__eq'] = kwargs.pop('shelf', None)
        kwargs['bin__eq'] = kwargs.pop('bin', None)
        
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, schema=None, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
            
        if schema is None:
            raise Exception('schema not initialized')
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            manual_field_includes.add('plate_location_id')
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                 
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            # specific setup 
 
            _p = self.bridge['plate']
            _pl = self.bridge['plate_location']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            
            plate_location_counts = (
                select([_pl.c.plate_location_id, func.count().label('plate_count')])
                    .select_from(
                        _pl.join(_p,_pl.c.plate_location_id
                            ==_p.c.plate_location_id))
                    .group_by(_pl.c.plate_location_id )
                ).cte('plate_counts')
            
            custom_columns = {
                'plate_count': literal_column('plate_counts.plate_count'),
                'libraries': (
                    select([func.array_to_string(
                        func.array_agg(distinct(_l.c.short_name)), 
                        LIST_DELIMITER_SQL_ARRAY)])
                        .select_from(
                            _l.join(_c,_l.c.library_id==_c.c.library_id)
                                .join(_p,_p.c.copy_id==_c.c.copy_id))
                        .where(_p.c.plate_location_id
                            ==literal_column('plate_location.plate_location_id')))
            };

            base_query_tables = ['plate', 'copy', 'plate_location', 'library']

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            # build the query statement
            
            j = _pl.join(
                    plate_location_counts, 
                    _pl.c.plate_location_id
                            ==plate_location_counts.c.plate_location_id,
                    isouter=True)
            stmt = select(columns.values()).select_from(j)

            if 'library_copy' in param_hash:
                param = param_hash.pop('library_copy')
                param = urllib.unquote(param).decode('utf-8')
                logger.info('param library_copy: %r', param)
                if '/' not in param or len(param.split('/')) !=2:
                    raise BadRequest(
                        'library_copy parameter must be of the form '
                        '"{library_short_name}/{copy_name}"')
                param = param.split('/')
                try:
                    library_copy = Copy.objects.get(
                        library__short_name=param[0],
                        name=param[1])
                    stmt = stmt.where(_pl.c.plate_location_id.in_(
                        select([_pl.c.plate_location_id])
                            .select_from(_pl.join(
                                _p,_pl.c.plate_location_id==_p.c.plate_location_id)
                                .join(_c,_p.c.copy_id==_c.c.copy_id))
                            .where(_c.c.copy_id==library_copy.copy_id)))
                except:
                    raise Http404('Copy not found: %r' % (param))
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
 
            if not order_clauses:
                stmt = stmt.order_by("room","freezer","shelf","bin")

            # logger.info(
            #     'stmt: %s',
            #     str(stmt.compile(
            #         dialect=postgresql.dialect(),
            #         compile_kwargs={"literal_binds": True})))

            # create a generator to wrap the cursor and expand copy_plates ranges
            def create_copy_plate_ranges_gen(generator):
                bridge = self.bridge
                _library = self.bridge['library']
                _lcp = self.bridge['plate']
                _cp = self.bridge['copy']
                lcp_query = (
                    select([ 
                        _library.c.short_name,
                        _cp.c.name,
                        _lcp.c.plate_number
                     ])
                    .select_from(
                        _lcp.join(_cp, _cp.c.copy_id == _lcp.c.copy_id)
                            .join(
                                _library,
                                _library.c.library_id == _cp.c.library_id))
                    .where(_lcp.c.plate_location_id == text(':plate_location_id'))
                    .group_by(
                        _library.c.short_name, _cp.c.name, _lcp.c.plate_number)
                    .order_by(
                        _library.c.short_name, _cp.c.name, _lcp.c.plate_number))
                logger.debug('lcp_query: %r', str(lcp_query.compile()))
                

                def copy_plate_ranges_generator(cursor):
                    if generator:
                        cursor = generator(cursor)
                    class Row:
                        def __init__(self, row):
                            self.row = row
                            self.entries = []
                            plate_location_id = row['plate_location_id']
                            query = conn.execute(
                                lcp_query, plate_location_id=plate_location_id)
                            copy = None
                            start_plate = None
                            end_plate = None
                            for x in query:
                                if not copy:
                                    copy = x[1]
                                    library = x[0]
                                if not start_plate:
                                    start_plate = end_plate = x[2]
                                if (x[0] != library 
                                    or x[1] != copy 
                                    or x[2] > end_plate + 1):
                                    # start a new range, save old range
                                    if start_plate != end_plate:
                                        self.entries.append('%s:%s:%s-%s'
                                            % (library, copy, start_plate, end_plate))
                                    else:
                                        self.entries.append('%s:%s:%s'
                                            % (library, copy, start_plate))
                                    start_plate = end_plate = x[2]
                                    copy = x[1]
                                    library = x[0]
                                else:
                                    end_plate = x[2]
                            if copy: 
                                self.entries.append('%s:%s:%s-%s'
                                    % (library, copy, start_plate, end_plate))
                                    
                        def has_key(self, key):
                            if key == 'copy_plate_ranges': 
                                return True
                            return self.row.has_key(key)
                        def keys(self):
                            return self.row.keys();
                        def __getitem__(self, key):
                            if key == 'copy_plate_ranges':
                                return self.entries
                            else:
                                return self.row[key]
                    conn = get_engine().connect()
                    # FIXME: use with get_engine().connect() as conn:
                    try:
                        for row in cursor:
                            yield Row(row)
                    finally:
                        conn.close()
                        
                return copy_plate_ranges_generator
            
            if 'copy_plate_ranges' in field_hash:
                
                rowproxy_generator = create_copy_plate_ranges_gen(rowproxy_generator)

            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])

            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash, param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None))

        except Exception, e:
            logger.exception('on get list')
            raise e   

    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_detail must be implemented')

    @write_authorization
    @un_cache  
    @transaction.atomic      
    def patch_detail(self, request, **kwargs):
        '''
        Override to generate informational summary for callee
        '''
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        deserialized = kwargs.pop('data', None)
        # allow for internal data to be passed
        if deserialized is None:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        logger.debug('patch detail %s, %s', deserialized,kwargs)

        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(deserialized, schema=schema, validate=True,**kwargs)

        original_data = None
        if kwargs_for_log:
            try:
                original_data = self._get_detail_response_internal(**kwargs_for_log)
            except Exception, e: 
                logger.exception('exception when querying for existing obj: %s', 
                    kwargs_for_log)
        try:
            log = kwargs.get('parent_log', None)
            if not log:
                log = self.make_log(request)
                log.save()
                kwargs['parent_log'] = log
            patch_response = self.patch_obj(request, deserialized, **kwargs)
            patched_plate_logs = patch_response[API_RESULT_OBJ]
        except ValidationError as e:
            logger.exception('Validation error: %r', e)
            raise e

        # get new state, for logging
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        update_log = self.log_patch(request, original_data,new_data,**kwargs)
        if update_log:
            update_log.save()
        patch_count = len(patched_plate_logs)
        update_count = len([x for x in patched_plate_logs if x.diffs ])
        unchanged_count = patch_count - update_count
        action = update_log.api_action if update_log else 'Unchanged'
        if action == API_ACTION_CREATE:
            action += ': ' + update_log.key
        meta = { 
            API_MSG_RESULT: { 
                'Plate location': action,
                API_MSG_SUBMIT_COUNT : patch_count, 
                API_MSG_UPDATED: update_count, 
                API_MSG_UNCHANGED: unchanged_count, 
                API_MSG_COMMENTS: log.comment
            }
        }
        if not self._meta.always_return_data:
            return self.build_response(
                request, {API_RESULT_META: meta }, response_class=HttpResponse, 
                **kwargs)
        else:
            return self.build_response(
                request,  {API_RESULT_META: meta }, response_class=HttpResponse, 
                **kwargs)
            
    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):

        logger.info('patch platelocation')
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_kwargs = self.get_id(deserialized, schema=schema, validate=True, **kwargs)

        create = False
        try:
            plate_location = PlateLocation.objects.get(**id_kwargs)
        except ObjectDoesNotExist:
            logger.info('plate location does not exist, creating: %r',
                id_kwargs)
            create=True
            plate_location = PlateLocation.objects.create(**id_kwargs)

        initializer_dict = self.parse(deserialized, schema=schema, create=create)
        errors = self.validate(initializer_dict, schema=schema,patch=not create)
        if errors:
            raise ValidationError(errors)
        
        copy_plate_ranges = deserialized.get(
            'copy_plate_ranges', [])
        all_plates = self._find_plates(
            schema['fields']['copy_plate_ranges']['regex'], copy_plate_ranges)
        
        # get the old plate locations, for logging
        plate_log_hash = {}
        lookup = ['room','freezer','shelf','bin']
        # combine all_plates (new) with plate_set.all() (previous)
        for plate in set(all_plates) | set(plate_location.plate_set.all()):
            plate_dict = {
                'library_short_name': plate.copy.library.short_name,
                'copy_name': plate.copy.name,
                'plate_number': str(plate.plate_number)
#                     'plate_number': str(plate.plate_number).zfill(5)
            }
            plate_log_hash[plate.plate_id] = [plate_dict, plate_dict.copy()]
            
            if plate.plate_location:
                plate_dict.update({
                    k:getattr(plate.plate_location,k, None)
                        for k in lookup })
        
        plate_location.plate_set = all_plates
        plate_location.save()
        
        for plate in set(all_plates) | set(plate_location.plate_set.all()):
            plate_dict = {
                'library_short_name': plate.copy.library.short_name,
                'copy_name': plate.copy.name,
                'plate_number': str(plate.plate_number)
#                     'plate_number': str(plate.plate_number).zfill(5)
            }
            if plate.plate_location:
                plate_dict.update({
                    k:getattr(plate.plate_location,k, None)
                        for k in lookup })
            plate_log_hash[plate.plate_id] = \
                [plate_log_hash[plate.plate_id][0],plate_dict]

        logger.info(
            'log copyplate patches for platelocation, '
            'len: %d...', len(plate_log_hash.items()))
        
        plate_logs = []
        lcp_id_attribute = \
            self.get_librarycopyplate_resource()\
                .build_schema(user=request.user)['id_attribute']
        for prev_dict,new_dict in plate_log_hash.values():
            log = self.get_librarycopyplate_resource().log_patch( 
                request,prev_dict,new_dict,
                **{ 'parent_log': kwargs.get('parent_log', None),
                    'full': True,
                    'id_attribute': lcp_id_attribute } )
            if log: 
                plate_logs.append(log)
        logger.info('logs created, saving...')
#             ApiLog.objects.bulk_create(plate_logs)
        ApiLog.bulk_create(plate_logs)
        logger.info('logs saved %r', plate_logs)
        
        return { API_RESULT_OBJ: plate_logs }
            
    # FIXME: replace with LibraryCopyPlateResource.find_plates (21070505)
    def _find_plates(self, regex_string, copy_plate_ranges):
        logger.info('find copy_plate_ranges: %r', copy_plate_ranges)
        
        if not regex_string:
            return []
        
        # parse library_plate_ranges
        # E.G. Regex: /(([^:]+):)?(\w+):(\d+)-(\d+)/
        matcher = re.compile(regex_string)
        
        # Expand the copy_plate_ranges
        new_copy_plate_ranges = []
        for copy_plate_range in copy_plate_ranges:
            match = matcher.match(copy_plate_range)
            if not match:
                raise ValidationError(
                    key='copy_plate_ranges',
                    msg=('%r does not match pattern: %s' 
                        % (copy_plate_range, regex_string)))
                break
            else:
                start_plate = int(match.group(4))
                if match.group(6):
                    end_plate = int(match.group(6))
                else:
                    end_plate = start_plate
                new_copy_plate_ranges.append({
                    'library_short_name': match.group(2),
                    'copy_name': match.group(3),
                    'start_plate': start_plate,
                    'end_plate': end_plate,
                     })
        copy_plate_ranges = new_copy_plate_ranges
        
        # validate/find the plate ranges
        logger.info(
            'get the referenced plates for: %r', copy_plate_ranges)
        all_plates = set()
        plates_changed_location = {}
        for _data in copy_plate_ranges:
            try:
                copy_name = _data['copy_name']
                library_short_name = _data['library_short_name']
                copy = Copy.objects.get(
                    name=copy_name,
                    library__short_name=library_short_name)
            except ObjectDoesNotExist:
                raise ValidationError(
                    key='copy_plate_ranges',
                    msg='{copy} does not exist: {val}'.format(
                        copy=copy_name,
                        val=str(_data)))
            try:
                start_plate = _data['start_plate']
                end_plate = _data['end_plate']
                plate_range = ( 
                    Plate.objects.all().filter(
                        copy=copy,
                        plate_number__range=[start_plate, end_plate])
                        .order_by('plate_number')
                    )
                if not plate_range:
                    raise ValidationError(
                        key='copy_plate_ranges',
                        msg=('plate range not found: {copy_name}{start_plate}-{end_plate}'
                            ).format(**_data))
                min = max = plate_range[0].plate_number
                for plate in plate_range:
                    if plate.plate_number < min:
                        min = plate.plate_number
                    if plate.plate_number > max:
                        max = plate.plate_number
                    all_plates.add(plate)
                if min != start_plate:
                    raise ValidationError(
                        key='copy_plate_ranges',
                        msg=('plate range start not found: {copy_name}:{start_plate}-{end_plate}'
                            ).format(**_data))
                if max != end_plate:
                    raise ValidationError(
                        key='copy_plate_ranges',
                        msg=('plate range end not found: {copy_name}:{start_plate}-{end_plate}'
                            ).format(**_data))
            except ObjectDoesNotExist:
                raise ValidationError(
                    key='copy_plate_ranges',
                    msg=('plate range not found: {copy_name}:{start_plate}-{end_plate}'
                        ).format(**_data))
        
        logger.debug('all plate keys: %r', [
            x for x in sorted(all_plates,key=lambda x: x.plate_number)])
        return all_plates

class LibraryCopyPlateResource(DbApiResource):
    retired_statuses = [
        'retired','discarded','lost','given_away','discarded_volume_transferred']

    class Meta:
        queryset = Plate.objects.all() 
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'librarycopyplate'
        authorization = UserGroupAuthorization(resource_name)
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):
        
        super(LibraryCopyPlateResource, self).__init__(**kwargs)
        self.plate_location_resource = None
        self.screen_resource = None
        self.library_screening_resource = None
        self.cpr_resource = None
        
    def get_cpr_resource(self):
        if self.cpr_resource is None:
            self.cpr_resource = CherryPickRequestResource()
        return self.cpr_resource
        
    def get_screen_resource(self):
        if not self.screen_resource:
            self.screen_resource = ScreenResource()
        return self.screen_resource
    
    def get_library_screening_resource(self):
        if not self.library_screening_resource:
            self.library_screening_resource = LibraryScreeningResource()
        return self.library_screening_resource
        
    def get_platelocation_resource(self):
        
        if not self.plate_location_resource:
            self.plate_location_resource = PlateLocationResource()
        return self.plate_location_resource
        
        
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/batch_edit%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('batch_edit'), name="api_lcp_batch_edit"),
            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('search'), name="api_search"),
            url(r"^(?P<resource_name>%s)/(?P<library_short_name>[\w.\-\+: ]+)"
                r"/(?P<copy_name>[\w.\-\+: ]+)"
                r"/(?P<plate_number>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<copy_name>[\w.\-\+: ]+)"
                r"/(?P<plate_number>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    @staticmethod
    def parse_plate_copy_search(plate_search_data):
        '''
        Create a plate search data structure:
            {   
                plates, plate_ranges, copies, 
                plate_numbers_expected, 
                plate_numbers_listed_in_order,
                plate_copy_keys_expected
            }
        '''
        logger.info('raw plate_search_data: %r', plate_search_data)
        # Use unquote to decode form data from a post
        if not isinstance(plate_search_data, (list,tuple)):
            plate_search_data = urllib.unquote(plate_search_data)
        else:
            plate_search_data = [urllib.unquote(x) for x in plate_search_data]
        logger.info('plate_search_data: %r', plate_search_data)
        parsed_searches = []
        
        # Process the patterns by line
        parsed_lines = plate_search_data
        if isinstance(parsed_lines, basestring):
            parsed_lines = re.split(
                lims_utils.PLATE_SEARCH_LINE_SPLITTING_PATTERN,parsed_lines)
            logger.info('parsed_lines: %r', parsed_lines)
        
        for _line in parsed_lines:
            _line = _line.strip()
            if not _line:
                continue

            parts = lims_utils.PLATE_RANGE_SPLITTING_PATTERN.findall(_line)
            logger.info('parse plate copy search: line parts: %r', parts)
            
            parsed_search = defaultdict(list)
            
            for part in parts:
                # unquote
                part = re.sub(r'["\']+','',part)

                if PLATE_PATTERN.match(part):
                    parsed_search['plate'].append(int(part))
                elif PLATE_RANGE_PATTERN.match(part):
                    match = PLATE_RANGE_PATTERN.match(part)
                    parsed_search['plate_range'].append(sorted([
                        int(match.group(1)), int(match.group(2))]))
                else:
                    # Must be a copy
                    if not COPY_NAME_PATTERN.match(part):
                        raise ValidationError(
                            key='plate_copy_search',
                            msg='unrecognized pattern: %r' % part)
                    logger.info('recognized copy: %r', part)
                    parsed_search['copy'].append(part)
                
            for k,v in parsed_search.items():
                parsed_search[k] = sorted(v)
            
            if parsed_search['copy']:
                plate_copy_keys_expected = set()
                for plate in parsed_search['plate']:
                    for copy in parsed_search['copy']:
                        plate_copy_keys_expected.add('%s/%s' % (copy,plate))
                for plate_range in parsed_search['plate_range']:
                    for plate in range(plate_range[0],plate_range[1]+1):
                        for copy in parsed_search['copy']:
                            plate_copy_keys_expected.add('%s/%s' % (copy,plate))
                parsed_search['plate_copy_keys_expected'] = \
                    list(plate_copy_keys_expected)
            else:
                plate_numbers_listed_in_order = list(parsed_search['plate'])
                plate_numbers_expected = set(parsed_search['plate'])
                plate_numbers_expected.update(*[
                    range(plate_range[0],plate_range[1]+1) 
                    for plate_range in parsed_search['plate_range']])
                parsed_search['plate_numbers_expected'] = \
                    sorted(plate_numbers_expected)
                for plate_range in parsed_search['plate_range']:
                    plate_numbers_listed_in_order.extend(
                        range(plate_range[0],plate_range[1]+1))
                parsed_search['plate_numbers_listed_in_order'] = \
                    plate_numbers_listed_in_order
            logger.debug('parsed: %r', parsed_search)
            parsed_searches.append(parsed_search)
        
        return parsed_searches
        
    @classmethod
    def find_plates(cls, plate_search_data):
        ''' 
        @return 
            set() of LibraryCopyPlate objects matching the
                plate_search_data (raw user search text)
            parsed_plate_search data structure 
                @see parse_plate_copy_search
            array of plates and/or plate_copy_keys expected but not found
        '''
        logger.info('find plates for patterns: %r', plate_search_data)

        if not plate_search_data:
            return []
        
        errors = set()
        plates = set()
        
        parsed_searches = cls.parse_plate_copy_search(plate_search_data)
        
        if not parsed_searches:
            return []
        
        for parsed_search in parsed_searches:
            logger.info('find plates for %r', parsed_search)
            plate_query = Plate.objects.all()
            
            plate_qs = []
            if parsed_search['plate']:
                plate_qs.append(Q(plate_number__in=parsed_search['plate']))
            if parsed_search['plate_range']:
                for plate_range in parsed_search['plate_range']:
                    plate_qs.append(Q(plate_number__range=plate_range))
            if plate_qs:
                plate_query = plate_query.filter(
                    reduce(lambda x,y: x|y,plate_qs))
            
            if parsed_search['copy']:
                plate_query = plate_query.filter(
                    copy__name__in=parsed_search['copy'])

            # check for not found:
            if parsed_search['copy']:
                plate_copy_keys_found = set([
                    '%s/%s' % (plate.copy.name, plate.plate_number)
                    for plate in plate_query.all()])
                not_found = (set(parsed_search['plate_copy_keys_expected'])
                    - plate_copy_keys_found)
                logger.info('plate-copies found: %r, expected: %r', 
                    sorted(plate_copy_keys_found), 
                    sorted(parsed_search['plate_copy_keys_expected']) )
                for plate_copy_key in not_found:
                    errors.add(plate_copy_key)
            else:
                plate_numbers_found = set([plate.plate_number 
                    for plate in plate_query.all() ])
                
                not_found = (set(parsed_search['plate_numbers_expected'])
                    - plate_numbers_found)
                logger.info('plates found: %r, expected: %r', 
                    sorted(plate_numbers_found), 
                    sorted(parsed_search['plate_numbers_expected']) )
                for plate_number in not_found:
                    errors.add(str(plate_number))
            plates.update(plate_query.all())
        
        if not plates:
#             raise ValidationError(key='Search', msg='No results found')
            errors.add('No plates found')
        logger.info('plates found: %d, errors: %r', len(plates), errors)
        return (plates, parsed_searches, [x for x in sorted(errors)])

    @read_authorization
    def get_detail(self, request, **kwargs):
        
        library_short_name = kwargs.get('library_short_name', None)
        if not library_short_name:
            logger.info('no library_short_name provided')
        copy_name = kwargs.get('copy_name', None)
        if not copy_name:
            raise Http404('must provide a copy_name parameter')
        plate_number = kwargs.get('plate_number', None)
        if not copy_name:
            raise Http404('must provide a plate_number parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    @classmethod
    def get_librarycopyplate_cte(cls, 
        filter_hash=None, library_id=None, copy_id=None, plate_ids=None):
        
        bridge = get_tables()
        _l = bridge['library']
        _c = bridge['copy']
        _p = bridge['plate']
        
        j = _p
        j = j.join(_c, _p.c.copy_id == _c.c.copy_id)
        j = j.join(_l, _c.c.library_id == _l.c.library_id)
        plate_table = (
            select([
                _p.c.plate_id,
                _p.c.plate_number,
                _p.c.copy_id,
                _c.c.library_id,
                _c.c.name.label('copy_name'),
                _concat(
                    _l.c.short_name, '/', _c.c.name, '/', 
                    cast(_p.c.plate_number, sqlalchemy.sql.sqltypes.Text)
                ).label('key'),
                ])
            .select_from(j))
        if library_id is not None:
            plate_table = plate_table.where(_c.c.library_id==library_id)
        if copy_id is not None:
            plate_table = plate_table.where(_c.c.copy_id==copy_id)
        if plate_ids is not None:
            plate_table = plate_table.where(_p.c.plate_id.in_(plate_ids))
        if filter_hash:
            extra_filters = [v for k,v in filter_hash.items()
                if k in ['plate_number', 'copy_name']]
            if extra_filters:
                logger.info('extra filters: %r', extra_filters)
                plate_table = (
                    select([literal_column(x) for x in plate_table.columns.keys()])
                    .select_from(Alias(plate_table))
                )
                # always OR at the subquery level
                plate_table = plate_table.where(or_(*extra_filters))

        return plate_table
        
    @classmethod
    def get_plate_copywell_statistics_cte(cls, 
        filter_hash=None, library_id=None, copy_id=None, plate_ids=None):  
        bridge = get_tables()
        _p = bridge['plate']
        _cw = bridge['copy_well']
        _c = bridge['copy']
        
        j = _cw
        if library_id is not None or copy_id is not None:
            j = j.join(_c, _cw.c.copy_id==_c.c.copy_id)
        cw_vols = (
            select([
                _cw.c.plate_id,
                func.count().label('count'),
                func.sum(_cw.c.volume).label('cum_well_vol'),
                func.min(_cw.c.volume).label('min_well_vol'),
                func.max(_cw.c.volume).label('max_well_vol'),
                func.min(_cw.c.mg_ml_concentration).label('min_well_mg_ml'),
                func.max(_cw.c.mg_ml_concentration).label('max_well_mg_ml'),
                func.min(_cw.c.molar_concentration).label('min_well_molar'),
                func.max(_cw.c.molar_concentration).label('max_well_molar'),
                ])
            .select_from(j)
            .group_by(_cw.c.plate_id))
        if library_id is not None:
            cw_vols = cw_vols.where(_c.c.library_id==library_id)
        if copy_id is not None:
            cw_vols = cw_vols.where(_c.c.copy_id==copy_id)
        # FIXME: adding plate_ids is not performant here
        # if plate_ids is not None:
        #     cw_vols = cw_vols.where(_p.c.plate_id.in_(plate_ids))
        cw_vols = cw_vols.cte('copy_well_volumes')
        
        # NOTE: this method of calculating plate volumes requires that plates
        # are never screened as "library_screenings" after being screened as
        # cherry pick screenings: (copy well statistics are set on the first
        # cherry pick screening, and are not reset on subsequent library_screenings)
        j2 = _p.outerjoin(cw_vols,_p.c.plate_id==cw_vols.c.plate_id)
        j2 = j2.join(_c, _p.c.copy_id==_c.c.copy_id)
        query = (
            select([
                _p.c.plate_id,
                _p.c.plate_number,
                _p.c.copy_id,
                _c.c.name.label('copy_name'),
                case([
                    (_p.c.experimental_well_count > 0,
                    ( func.coalesce(_p.c.remaining_well_volume,_p.c.well_volume) 
                        * (_p.c.experimental_well_count-cw_vols.c.count)
                        + cw_vols.c.cum_well_vol ) / _p.c.experimental_well_count)],
                    else_=None).label('avg_well_remaining_volume'),
                func.coalesce(_p.c.remaining_well_volume,_p.c.well_volume)
                    .label('remaining_well_volume'),
                func.coalesce(cw_vols.c.min_well_vol,
                    _p.c.remaining_well_volume,).label('min_well_remaining_volume'),
                func.coalesce(cw_vols.c.max_well_vol,
                    _p.c.remaining_well_volume).label('max_well_remaining_volume'),
                func.coalesce(cw_vols.c.min_well_mg_ml,
                    _p.c.mg_ml_concentration).label('min_mg_ml_concentration'),
                func.coalesce(cw_vols.c.max_well_mg_ml,
                    _p.c.mg_ml_concentration).label('max_mg_ml_concentration'),
                func.coalesce(cw_vols.c.min_well_molar,
                    _p.c.molar_concentration).label('min_molar_concentration'),
                func.coalesce(cw_vols.c.max_well_molar,
                    _p.c.molar_concentration).label('max_molar_concentration'),
                # (select([func.min(_cw.c.volume)])
                #     .select_from(_cw).where(_cw.c.plate_id==_p.c.plate_id)
                #     .label('min_well_remaining_volume')),
                # (select([func.max(_cw.c.volume)])
                #     .select_from(_cw).where(_cw.c.plate_id==_p.c.plate_id)
                #     .label('max_well_remaining_volume')),
                # (select([func.coalesce(
                #     func.min(_cw.c.mg_ml_concentration),
                #     _p.c.mg_ml_concentration)])
                #     .select_from(_cw).where(_cw.c.plate_id==_p.c.plate_id)
                #     .label('min_mg_ml_concentration')),
                # (select([func.coalesce(
                #     func.max(_cw.c.mg_ml_concentration),
                #     _p.c.mg_ml_concentration)])
                #     .select_from(_cw).where(_cw.c.plate_id==_p.c.plate_id)
                #     .label('max_mg_ml_concentration')),
                # (select([func.coalesce(
                #     func.min(_cw.c.molar_concentration),
                #     _p.c.molar_concentration)])
                #     .select_from(_cw).where(_cw.c.plate_id==_p.c.plate_id)
                #     .label('min_molar_concentration')),
                # (select([func.coalesce(
                #     func.max(_cw.c.molar_concentration),
                #     _p.c.molar_concentration)])
                #     .select_from(_cw).where(_cw.c.plate_id==_p.c.plate_id)
                #     .label('max_molar_concentration')),
                ])
            .select_from(j2)
            )
        if library_id is not None:
            query = query.where(_c.c.library_id==library_id)
        if copy_id is not None:
            query = query.where(_p.c.copy_id==copy_id)
        if plate_ids is not None:
            query = query.where(_p.c.plate_id.in_(plate_ids))

        if filter_hash:
            # Create a subquery and filter on the expressions that match the fields
            extra_filters = [v for k,v in filter_hash.items()
                if k in ['plate_number', 'copy_name']]
            if extra_filters:
                query = (
                    select([literal_column(x) for x in query.columns.keys()])
                    .select_from(Alias(query))
                )
                # Always OR at the subquery level
                query = query.where(or_(*extra_filters))
        
        return query
    
    @classmethod
    def get_plate_screening_statistics_cte(cls, 
        filter_hash=None, library_id=None, copy_id=None, plate_ids=None):
        bridge = get_tables()
        _p = bridge['plate']
        _a = bridge['activity']
        _ls = bridge['library_screening']
        _ap = bridge['assay_plate']
        _c = bridge['copy']

        j = ( _ap.join(_ls,_ls.c.activity_id==_ap.c.library_screening_id)
                .join(_a, _a.c.activity_id==_ls.c.activity_id)
                .join(_p, _ap.c.plate_id==_p.c.plate_id)
                .join(_c, _p.c.copy_id==_c.c.copy_id))
        
        query = (
            select([
                _p.c.copy_id,
                _c.c.name.label('copy_name'),
                _ap.c.plate_id,
                _ap.c.plate_number,
                func.min(_a.c.date_of_activity).label('first_date_screened'),
                func.max(_a.c.date_of_activity).label('last_date_screened')
                ])
            .select_from(j)
            .group_by(_ap.c.plate_id, _ap.c.plate_number, _p.c.copy_id,
                _c.c.name ) )
        if library_id is not None:
            query = query.where(_c.c.library_id == library_id)
        if copy_id is not None:
            query = query.where(_p.c.copy_id==copy_id)
        if plate_ids is not None:
            query = query.where(_p.c.plate_id.in_(plate_ids))
        if filter_hash:
            extra_filters = [v for k,v in filter_hash.items()
                if k in ['plate_number', 'copy_name']]
            if extra_filters:
                query = (
                    select([literal_column(x) for x in query.columns.keys()])
                    .select_from(Alias(query))
                )
                # Always OR at the subquery level
                query = query.where(or_(*extra_filters))
        ## TODO: this does not include cherry pick screenings, redo when 
        # cherry pick screenings are reworked:
        # rework: so that source plates are linked to screening activity
        
        return query
        
    def build_screened_plate_response(
        self, request, screen_facility_id=None, library_screening_id=None,
         **kwargs):    
        
        if screen_facility_id is None and library_screening_id is None:
            raise NotImplementedError(
                'must provide a "screen_facility_id" or "library_screening_id"')
        if screen_facility_id is not None and library_screening_id is not None:
            raise NotImplementedError(
                'must provide either a "screen_facility_id" or "library_screening_id"')
        if screen_facility_id is not None:
            if self.get_screen_resource()._meta.authorization\
                .has_screen_read_authorization(
                    request.user, screen_facility_id) is False:
                raise PermissionDenied
        if library_screening_id is not None:
            if self.get_library_screening_resource()._meta.authorization\
                .has_activity_read_authorization(
                    request.user, library_screening_id) is False:
                raise PermissionDenied
            
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False

        manual_field_includes = set(param_hash.get('includes', []))
        logger.info('manual_field_includes: %r', manual_field_includes)
        
        if screen_facility_id is not None:
            filename = 'plates_for_screen_%s' % screen_facility_id
        else:
            filename = 'plates_for_screening_%s' % library_screening_id
            
        # construct a limited plate view here 
        # start from the librarycopyplate schema, but limit to plate only fields
        schema = self.build_schema(user=request.user);
        
        fields_to_show = [
            'library_short_name', 'library_screening_status', 
            'library_comment_array','plate_number','comment_array', 
            'screening_count','assay_plate_count',
            'last_date_screened','first_date_screened']
        
        new_fields = {}
        for k,v in schema['fields'].items():
            if k in fields_to_show:
                new_fields[k] = v
                v['visibility'] = ['l']
        # Add a "copies_screened" field
        new_fields['copies_screened'] = schema['fields']['copy_name']
        new_fields['copies_screened']['key'] = 'copies_screened'
        new_fields['copies_screened']['title'] = 'Copies Screened'
        new_fields['copies_screened']['data_type'] = 'list'
        schema['fields'] = new_fields
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filter_expression = \
            self._meta.authorization.filter(request.user,filter_expression)
        visibilities = ['l','d']
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            visibilities,
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
         
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            # use "use_vocab" as a proxy to also adjust siunits for viewing
            rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                field_hash, rowproxy_generator)
        rowproxy_generator = \
            self._meta.authorization.get_row_property_generator(
                request.user, field_hash, rowproxy_generator)

        # specific setup 
        _p = self.bridge['plate']
        _c = self.bridge['copy']
        _l = self.bridge['library']
        _ls = self.bridge['library_screening']
        _la = self.bridge['lab_activity']
        _a = self.bridge['activity']
        _ap = self.bridge['assay_plate']
        _screen = self.bridge['screen']
        _apilog = self.bridge['reports_apilog']
        _assay_plates_query = (select([
            _ls.c.activity_id,
            _a.c.date_of_activity,
            _ap.c.plate_id,
            _c.c.copy_id,
            _c.c.name.label('copy_name'),
            _l.c.short_name,
            _c.c.library_id,
            _ap.c.plate_number,
            _concat(
                _l.c.short_name, '/', _c.c.name, '/', 
                cast(_p.c.plate_number, sqlalchemy.sql.sqltypes.Text)
            ).label('plate_key'),
            ])
            .select_from(
                _ap.join(_ls,_ap.c.library_screening_id==_ls.c.activity_id)
                   .join(_a,_ls.c.activity_id==_a.c.activity_id)
                   .join(_p, _ap.c.plate_id==_p.c.plate_id)
                   .join(_c, _c.c.copy_id==_p.c.copy_id)
                   .join(_l, _c.c.library_id==_l.c.library_id)
                   .join(_la,_ls.c.activity_id==_la.c.activity_id)
                   .join(_screen,_la.c.screen_id==_screen.c.screen_id)
                )
            .where(_ap.c.replicate_ordinal==0))
        
        if screen_facility_id is not None:
            _assay_plates_query = _assay_plates_query.where(
                _screen.c.facility_id==screen_facility_id)
        if library_screening_id is not None:
            _assay_plates_query = _assay_plates_query.where(
                _ls.c.activity_id==library_screening_id)
        _assay_plates = _assay_plates_query.cte('assay_plates')
        _assay_plates_inner = _assay_plates_query.cte('assay_plates_inner')
        _plate_comment_apilogs = ApiLogResource.get_resource_comment_subquery(
            self._meta.resource_name)
        _library_comment_apilogs = \
            ApiLogResource.get_resource_comment_subquery('library')
        
        # Get the library and plate keys to prefilter the logs
        with get_engine().connect() as conn:
            library_names = [x[0] for x in 
                conn.execute(select([distinct(_assay_plates.c.short_name)])
                    .select_from(_assay_plates))]
            plate_keys = [x[0] for x in 
                conn.execute(select([
                        distinct(_assay_plates.c.plate_key)])
                    .select_from(_assay_plates)) ]
            _plate_comment_apilogs = _plate_comment_apilogs.where(
                _apilog.c.key.in_(plate_keys))
            _library_comment_apilogs = _library_comment_apilogs.where(
                _apilog.c.key.in_(library_names))
        _library_comment_apilogs = \
            _library_comment_apilogs.cte('_library_comment_apilogs')
        _plate_comment_apilogs = _plate_comment_apilogs.cte('_comment_apilogs')

        stmt = ( 
            select([
                _assay_plates.c.short_name.label('library_short_name'),
                _l.c.library_name.label('library_name'),
                _l.c.screening_status.label('library_screening_status'),
                ( select([
                    func.array_to_string(
                        func.array_agg(literal_column('name')), 
                        LIST_DELIMITER_SQL_ARRAY),
                    ])
                    .select_from(
                        select([distinct(_c.c.name)])
                        .select_from(_c.join(
                            _assay_plates_inner,_c.c.copy_id
                                ==_assay_plates_inner.c.copy_id))
                        .where(_assay_plates_inner.c.plate_number
                            ==text('assay_plates.plate_number'))
                        .order_by(_c.c.name).alias('inner_copies'))
                    ).label('copies_screened'),
                ( select([func.min(_assay_plates_inner.c.date_of_activity)])
                    .select_from(_assay_plates_inner)
                    .where(_assay_plates_inner.c.plate_number
                        ==text('assay_plates.plate_number'))
                    ).label('first_date_screened'),
                ( select([func.max(_assay_plates_inner.c.date_of_activity)])
                    .select_from(_assay_plates_inner)
                    .where(_assay_plates_inner.c.plate_number
                        ==text('assay_plates.plate_number'))
                    ).label('last_date_screened'),
                _assay_plates.c.plate_number,
                func.count(None).label('assay_plate_count'),
                func.count(distinct(_assay_plates.c.activity_id))
                    .label('screening_count'),
                (
                    select([func.array_to_string(
                        func.array_agg(
                            _concat(                            
                                cast(_library_comment_apilogs.c.name,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                cast(_library_comment_apilogs.c.date_time,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                _library_comment_apilogs.c.comment)
                        ), 
                        LIST_DELIMITER_SQL_ARRAY) ])
                    .select_from(_library_comment_apilogs)
                    .where(_library_comment_apilogs.c.key
                        ==_assay_plates.c.short_name)
                    ).label('library_comment_array'),
                (
                    select([func.array_to_string(
                        func.array_agg(
                            _concat(                            
                                cast(_plate_comment_apilogs.c.name,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                cast(_plate_comment_apilogs.c.date_time,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                _plate_comment_apilogs.c.comment)
                        ), 
                        LIST_DELIMITER_SQL_ARRAY) ])
                    .select_from(_plate_comment_apilogs)
                    .where(_plate_comment_apilogs.c.key.ilike(
                        _concat(
                            _assay_plates.c.short_name,
                            '/%/',
                            cast(_assay_plates.c.plate_number, 
                                sqlalchemy.sql.sqltypes.Text))))
                    ).label('comment_array'),
            ])
            .select_from(
                _assay_plates.join(
                    _l,_assay_plates.c.library_id==_l.c.library_id))
            .group_by(
                _assay_plates.c.short_name,
                _l.c.library_name,
                _l.c.screening_status,
                _assay_plates.c.plate_number) )

        # general setup
        
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        
        if not order_clauses:
            stmt = stmt.order_by("plate_number")
        
        # compiled_stmt = str(stmt.compile(
        # dialect=postgresql.dialect(),
        # compile_kwargs={"literal_binds": True}))
        # logger.info('compiled_stmt %s', compiled_stmt)
        
        title_function = None
        if use_titles is True:
            def title_function(key):
                return field_hash[key]['title']
        if is_data_interchange:
            title_function = DbApiResource.datainterchange_title_function(
                field_hash,schema['id_attribute'])
        
        return self.stream_response_from_statement(
            request, stmt, count_stmt, filename,
            field_hash=field_hash, param_hash=param_hash,
            is_for_detail=False,
            rowproxy_generator=rowproxy_generator,
            title_function=title_function, meta=kwargs.get('meta', None))
        
    def build_list_response(self, request, schema=None, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False

        if schema is None:
            raise Exception('schema not initialized')
        
        manual_field_includes = set(param_hash.get('includes', []))
        logger.info('manual_field_includes: %r', manual_field_includes)

        # FIXME: plates_screened has been replaced by 
        # librarycopyplate.build_screened_plate_response?
        for_screen_facility_id = param_hash.pop('for_screen_facility_id', None)
        if for_screen_facility_id is not None:
            authorized = self.get_screen_resource()._meta.authorization\
                .has_screen_read_authorization(
                    request.user, for_screen_facility_id)
            if authorized is False:
                raise PermissionDenied
                
        library_screening_id = param_hash.pop('library_screening_id', None)
        if library_screening_id is not None:
            authorized = self.get_library_screening_resource()._meta\
                .authorization.has_activity_read_authorization(
                    request.user, library_screening_id)
            if authorized is False:
                raise PermissionDenied

        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id is not None:
            authorized = self.get_cpr_resource()._meta.authorization\
                .has_cherry_pick_read_authorization(request.user, cherry_pick_request_id)
            if authorized is False:
                raise PermissionDenied
            param_hash['cpr'] = cherry_pick_request_id
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        library_id = None
        library_short_name = param_hash.pop('library_short_name',
            param_hash.get('library_short_name__eq', None))
        if not library_short_name:
            logger.info('no library_short_name provided')
        else:
            param_hash['library_short_name__eq'] = library_short_name
            library_id = Library.objects.get(short_name=library_short_name).library_id
        
        copy_id = None
        copy_name = param_hash.pop('copy_name',
            param_hash.get('copy_name', None))
        if copy_name and library_id:
            param_hash['copy_name__eq'] = copy_name
            copy_id = Copy.objects.get(library_id=library_id, name=copy_name).copy_id
        plate_number = param_hash.pop('plate_number',
            param_hash.get('plate_number', None))
        if plate_number:
            param_hash['plate_number__eq'] = plate_number
        plate_ids = param_hash.pop('plate_ids', None)
        if plate_ids is not None:
            if isinstance(plate_ids,basestring):
                plate_ids = [int(x) for x in plate_ids.split(',')]
        plate_search_data = param_hash.pop('raw_search_data', None)
        logger.info('plate raw_search_data: %r', plate_search_data)
        
        if len(filter(lambda x: x is not None, 
            [cherry_pick_request_id, for_screen_facility_id, 
                library_screening_id, plate_ids, plate_search_data]))>1:
            raise NotImplementedError('Mutually exclusive params: %r'
                % ['cherry_pick_request_id', 'for_screen_facility_id', 
                    'library_screening_id', 'plate_ids','plate_search_data'])
            
        log_key = '/'.join(str(x) if x is not None else '%' 
            for x in [library_short_name,copy_name,plate_number])
        if log_key == '%/%/%':
            log_key = None

        try:
            # Use cherry_pick_request_id, screen_id, library_screening_id, or 
            # search data to pre-filter for plate_ids
            # TODO: grab keys as well for log query performance
            
            if cherry_pick_request_id is not None:
                _lcp = self.bridge['lab_cherry_pick']
                _well = self.bridge['well']
                _p = self.bridge['plate']
                cpr_plates = (
                    select([distinct(_p.c.plate_id).label('plate_id')])
                    .select_from(
                        _lcp.join(_well,_lcp.c.source_well_id==_well.c.well_id)
                            .join(_p,and_(
                                _p.c.plate_number==_well.c.plate_number,
                                _p.c.copy_id==_lcp.c.copy_id)))
                    .where(_lcp.c.cherry_pick_request_id==cherry_pick_request_id)
                    #.where(_lcp.c.selected==True)
                )
                with get_engine().connect() as conn:
                    plate_ids = [x[0] for x in 
                        conn.execute(cpr_plates)]
            
            # FIXME: plates_screened has been replaced by 
            # librarycopyplate.build_screened_plate_response?
            if for_screen_facility_id:
                manual_field_includes.add('assay_plate_count')
                # plates screened
                with get_engine().connect() as conn:
                    plate_ids = [x[0] for x in 
                        conn.execute(
                            self.get_screen_librarycopyplate_subquery(for_screen_facility_id))]
            if library_screening_id is not None:
                manual_field_includes.add('assay_plate_count')
                with get_engine().connect() as conn:
                    plate_ids = [x[0] for x in 
                        conn.execute(
                            self.get_libraryscreening_plate_subquery(library_screening_id))]

            # Parse plate search OR-clauses:
            # POST data: line delimited plate-range text entered by the user
            # Embedded plate_search url params: array of plate-ranges
            plate_search_errors = set()
            if plate_search_data is not None:
                (plates,parsed_searches, errors) = self.find_plates(plate_search_data)
                plate_search_errors.update(errors)
                plate_ids = [p.plate_id for p in plates]
                
            # general setup
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                 
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            # specific setup 
 
            _apilog = self.bridge['reports_apilog']
            _diff = self.bridge['reports_logdiff']
            _p = self.bridge['plate']
            _well = self.bridge['well']
            _pl = self.bridge['plate_location']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _ls = self.bridge['library_screening']
#             _la = self.bridge['lab_activity']
            _a = self.bridge['activity']
            _ap = self.bridge['assay_plate']
            _screen = self.bridge['screen']
            _plate_cte = (
                self.get_librarycopyplate_cte(
                    filter_hash=filter_hash,
                    library_id=library_id, copy_id=copy_id, plate_ids=plate_ids)
                        .cte('plate_cte'))
            _plate_statistics = (
                self.get_plate_copywell_statistics_cte(
                    filter_hash=filter_hash,
                    library_id=library_id, copy_id=copy_id, plate_ids=plate_ids)
                        .cte('plate_statistics'))
            _plate_screening_statistics = (
                self.get_plate_screening_statistics_cte(
                    filter_hash=filter_hash,
                    library_id=library_id, copy_id=copy_id, plate_ids=plate_ids)
                        .cte('plate_screening_statistics'))
            
            # Status Date/performedBy: Use window function for performance:
            #   SELECT DISTINCT ON (reports_apilog.key)
            #   last_value(date_time) OVER wnd as date_time,
            #   last_value(reports_apilog.id) OVER wnd as id,
            #   reports_apilog.key
            #  FROM reports_apilog
            #  where ref_resource_name = 'librarycopyplate'
            #  WINDOW wnd AS (
            #    PARTITION BY reports_apilog.id ORDER BY date_time
            #    ROWS BETWEEN UNBOUNDED PRECEDING AND UNBOUNDED FOLLOWING )
            w_params = { 
                'order_by': _apilog.c.date_time,
                'partition_by': _apilog.c.id }
            status_apilogs = (
                select([
                    func.last_value(_apilog.c.date_time).over(**w_params).label('date_time'),
                    func.last_value(_apilog.c.id).over(**w_params).label('id'),
                    _apilog.c.key
                    ])
                .select_from(_apilog.join(_diff, _apilog.c.id==_diff.c.log_id))
                .where(_apilog.c.ref_resource_name==self._meta.resource_name)
                .where(_diff.c.field_key=='status')
                .distinct(_apilog.c.key)
                )
            if log_key:
                status_apilogs = status_apilogs.where(_apilog.c.key.ilike(log_key))
            status_updated_apilogs = status_apilogs.cte('updated_apilogs')
            _user_cte = ScreensaverUserResource.get_user_cte().cte('user_cte')
            _status_apilogs = (
                select([
                    _apilog.c.date_time,
                    _apilog.c.key,
                    _user_cte.c.username,
                    _user_cte.c.name,                    
                    ])
                .select_from(_apilog.join(
                    status_updated_apilogs,
                        _apilog.c.id==status_updated_apilogs.c.id)
                    .join(_user_cte, _apilog.c.username==_user_cte.c.username))
                )
            _comment_apilogs = ApiLogResource.get_resource_comment_subquery(
                self._meta.resource_name)
            if log_key:
                _status_apilogs = _status_apilogs.where(_apilog.c.key.ilike(log_key))
                _comment_apilogs = _comment_apilogs.where(_apilog.c.key.ilike(log_key))
            _status_apilogs = _status_apilogs.cte('_status_apilogs')
            _comment_apilogs = _comment_apilogs.cte('_comment_apilogs')

            _library_comment_apilogs = \
                ApiLogResource.get_resource_comment_subquery('library')
            if library_short_name is not None:
                _library_comment_apilogs = _library_comment_apilogs.where(
                    _apilog.c.key==library_short_name)
            _library_comment_apilogs = \
                _library_comment_apilogs.cte('_library_comment_apilogs')

            
            custom_columns = {
                'assay_plate_count': (
                    select([func.count(None)])
                        .select_from(_ap)
                        .where(_ap.c.plate_id==_plate_cte.c.plate_id)),
                'librarycopyplate_id': (
                    _concat(_l.c.short_name,'/',_c.c.name,'/', 
                        cast(_p.c.plate_number, sqlalchemy.sql.sqltypes.Text))),
                'copyplate_id': (
                    _concat(_c.c.name,'/',
                        cast(_p.c.plate_number, sqlalchemy.sql.sqltypes.Text))),
                'remaining_well_volume': _plate_statistics.c.remaining_well_volume,
                'avg_remaining_volume': _plate_statistics.c.avg_well_remaining_volume,
                'min_remaining_volume': _plate_statistics.c.min_well_remaining_volume,
                'max_remaining_volume': _plate_statistics.c.max_well_remaining_volume,
                'min_molar_concentration': _plate_statistics.c.min_molar_concentration,
                'max_molar_concentration': _plate_statistics.c.max_molar_concentration,
                'min_mg_ml_concentration': _plate_statistics.c.min_mg_ml_concentration,
                'max_mg_ml_concentration': _plate_statistics.c.max_mg_ml_concentration,
#                 'plate_number': (literal_column(
#                     "to_char(plate.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
#                     .label('plate_number')),
                'first_date_screened': _plate_screening_statistics.c.first_date_screened,
                'last_date_screened': _plate_screening_statistics.c.last_date_screened,
                'status_date': _status_apilogs.c.date_time,
                'status_performed_by': _status_apilogs.c.name,
                'status_performed_by_username': _status_apilogs.c.username,
                'experimental_copy_well_count': (
                    select([func.count(None)])
                    .select_from(_well)
                    .where(_well.c.plate_number==_p.c.plate_number)
                    .where(_well.c.library_well_type=='experimental')),

                'comment_array': (
                    select([func.array_to_string(
                        func.array_agg(
                            _concat(                            
                                cast(_comment_apilogs.c.name,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                cast(_comment_apilogs.c.date_time,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                _comment_apilogs.c.comment)
                        ), 
                        LIST_DELIMITER_SQL_ARRAY) ])
                    .select_from(_comment_apilogs)
                    .where(_comment_apilogs.c.key==_plate_cte.c.key)
                    ),
                'library_comment_array': (
                    select([func.array_to_string(
                        func.array_agg(
                            _concat(                            
                                cast(_library_comment_apilogs.c.name,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                cast(_library_comment_apilogs.c.date_time,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                _library_comment_apilogs.c.comment)
                        ), 
                        LIST_DELIMITER_SQL_ARRAY) ])
                    .select_from(_library_comment_apilogs)
                    .where(_library_comment_apilogs.c.key==_l.c.short_name)
                    ),
            };

            base_query_tables = ['plate', 'copy', 'plate_location', 'library']

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            # build the query statement

            j = join(_plate_cte, _p, _p.c.plate_id == _plate_cte.c.plate_id)
            j = j.join(_c, _p.c.copy_id == _c.c.copy_id)
            j = j.join(
                _pl, _p.c.plate_location_id == _pl.c.plate_location_id,
                isouter=True)
            j = j.join(_l, _c.c.library_id == _l.c.library_id)

            if set(['status_date','status_performed_by',
                'status_performed_by_username']) | set(field_hash.keys()):
                j = j.outerjoin(_status_apilogs,_plate_cte.c.key==_status_apilogs.c.key)
            if set(['avg_remaining_volume','min_remaining_volume',
                'max_remaining_volume']) | set(field_hash.keys()):
                j = j.outerjoin(_plate_statistics,_p.c.plate_id==_plate_statistics.c.plate_id)     
            if (set(['first_date_screened', 'last_date_screened'])
                | set(field_hash.keys()) ):
                j = j.outerjoin(_plate_screening_statistics,
                    _p.c.plate_id==_plate_screening_statistics.c.plate_id)

            # FIXME: replaced by build_screened_plate_response? 20170912
            if for_screen_facility_id is not None:
                custom_columns['screening_count'] = (
                    select([distinct(_ls.c.activity_id)])
                        .select_from(_ap.join(
                            _screen, _ap.c.screen_id==_screen.c.screen_id))
                        .where(_ap.c.plate_id==_plate_cte.c.plate_id)
                        .where(_screen.c.facility_id==for_screen_facility_id)
                )
                custom_columns['assay_plate_count'] = (
                    select([func.count(None)])
                        .select_from(_ap.join(
                            _screen, _ap.c.screen_id==_screen.c.screen_id))
                        .where(_ap.c.plate_id==_plate_cte.c.plate_id)
                        .where(_screen.c.facility_id==for_screen_facility_id)
                )
                custom_columns['first_date_screened'] = (
                    select([func.min(_a.c.date_of_activity)])
                        .select_from(
                            _ap.join(_screen, _ap.c.screen_id==_screen.c.screen_id)
                               .join(_a, _ap.c.library_screening_id==_a.c.activity_id)
                            )
                        .where(_ap.c.plate_id==_plate_cte.c.plate_id)
                        .where(_screen.c.facility_id==for_screen_facility_id)
                    )
                custom_columns['last_date_screened'] = (
                    select([func.max(_a.c.date_of_activity)])
                        .select_from(
                            _ap.join(_screen, _ap.c.screen_id==_screen.c.screen_id)
                               .join(_a, _ap.c.library_screening_id==_a.c.activity_id)
                            )
                        .where(_ap.c.plate_id==_plate_cte.c.plate_id)
                        .where(_screen.c.facility_id==for_screen_facility_id)
                    )
                
            if library_screening_id is not None:
                custom_columns['assay_plate_count'] = (
                    select([func.count(None)])
                        .select_from(_ap)
                        .where(_ap.c.plate_id==_plate_cte.c.plate_id)
                        .where(_ap.c.library_screening_id==library_screening_id)
                ),
                

            
            stmt = select(columns.values()).select_from(j)
            
            if plate_ids is not None:
                stmt = stmt.where(_p.c.plate_id.in_(plate_ids))
            if library_id is not None:
                stmt = stmt.where(_l.c.library_id==library_id)
            if copy_id is not None:
                stmt = stmt.where(_p.c.copy_id==copy_id)

            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
 
            if not order_clauses:
                stmt = stmt.order_by("plate_number", "copy_name")

            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
         
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            meta = kwargs.get('meta', {})
            if plate_search_errors:
                meta['API_MSG_WARNING'] = \
                    'Plates not found: %s' % ', '.join(sorted(plate_search_errors))
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash, param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=meta)
                        
        except Exception, e:
            logger.exception('on get list')
            raise e   

    @classmethod
    def get_screen_librarycopyplate_subquery(cls, for_screen_facility_id):
        
        bridge = get_tables()
        _screen = bridge['screen']
        _assay_plate = bridge['assay_plate']
        _plate = bridge['plate']
        
        j = _plate
        j = j.join(_assay_plate, _plate.c.plate_id == _assay_plate.c.plate_id)
        j = j.join(_screen, _assay_plate.c.screen_id == _screen.c.screen_id)
        
        screen_lcps = (
            select([
                distinct(_plate.c.plate_id).label('plate_id')])
            .select_from(j)
            .where(_screen.c.facility_id == for_screen_facility_id)
            .where(_assay_plate.c.replicate_ordinal == 0))
        return screen_lcps

    @classmethod
    def get_libraryscreening_plate_subquery(cls, library_screening_id):
        
        bridge = get_tables()
        _assay_plate = bridge['assay_plate']
        _plate = bridge['plate']
        
        j = _plate
        j = j.join(_assay_plate, _plate.c.plate_id == _assay_plate.c.plate_id)
        
        _query = (
            select([
                distinct(_plate.c.plate_id).label('plate_id')])
            .select_from(j)
            .where(_assay_plate.c.library_screening_id == library_screening_id)
            .where(_assay_plate.c.replicate_ordinal == 0))
        return _query

    @write_authorization
    @un_cache        
    @transaction.atomic    
    def batch_edit(self, request, **kwargs):
        '''
        Batch edit is a POST operation:
        
        librarycopyplate batch_edit uses a POST form to send both
        the search data (3 types of search filter: "nested_search_data",  
        "raw_search_data", and GET search params),
        as well as the update data (in the form of "plate_info" and "plate_location")
        (Instead of sending all plates to be PATCHED);
        '''
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        # NOTE: authentication is being performed by wrap_view
        # self.is_authenticated(request)
        # if not self._meta.authorization._is_resource_authorized(
        #         self._meta.resource_name,request.user,'write'):
        #     raise ImmediateHttpResponse(
        #         response=HttpForbidden(
        #             'user: %s, permission: %s/%s not found' 
        #             % (request.user,self._meta.resource_name,'write')))
        convert_post_to_put(request)
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        logger.info('batch_edit plates...')
        logger.info('param_hash: %r', param_hash)
        
        plate_info_data = param_hash.pop('plate_info', None)
        if plate_info_data: 
            plate_info_data = json.loads(plate_info_data)
        logger.info('plate_info: %r', plate_info_data)
        
        location_data = param_hash.pop('plate_location', None)
        if location_data:
            location_data = json.loads(location_data)

        if plate_info_data is not None and location_data is not None:
            raise BadRequest('batch info: edit only one of %r'
                % ['plate_info', 'plate_location'])
        elif plate_info_data is None and location_data is None:
            raise BadRequest('"data" must contain one of %r'
                % ['plate_info', 'plate_location'])
        
        # Use the rest of the POST data for the plate search
        nested_search_data = param_hash.pop('nested_search_data', None)
        plate_kwargs = param_hash.copy()
        if nested_search_data:
            plate_kwargs['nested_search_data'] = json.loads(nested_search_data)
        
        logger.info('plate_kwargs: %r', plate_kwargs)
        
        if plate_info_data is not None:
            # Find all of the plates
            plate_kwargs['includes'] = ['plate_id']
            original_data = self._get_list_response_internal(**plate_kwargs)
            if not original_data:
                raise Http404
            logger.info('got the plate objects for search data: %d', len(original_data))
            logger.debug('plates %r', ['{copy_name}/{plate_number}'.format(**lcp)
                for lcp in original_data])
            
            plate_ids = [x['plate_id'] for x in original_data]

            query = Plate.objects.all().filter(plate_id__in=plate_ids)
            for key,value in plate_info_data.items():
                fi = schema['fields'].get(key,None)
                if fi is None or 'u' not in fi.get('editability',[]):
                    raise ValidationError(
                        key=key, msg='is not editable')
                query.update(**{ key: value })
                if key=='status':
                    if value == 'available':
                        query.update(**{ 'date_plated': _now()})
                        # query.update(**{ 'date_retired': None})
                    elif value in self.retired_statuses:
                        query.update(**{ 'date_retired': _now()})
                    
            logger.info('batch update complete, logging...')
            
            new_data = self._get_list_response_internal(**plate_kwargs)
            
            parent_log = self.make_log(request)
            parent_log.key = self._meta.resource_name
            parent_log.uri = self._meta.resource_name
            parent_log.save()
            plate_kwargs['parent_log'] = parent_log
            logs = self.log_patches(request, original_data, new_data,**plate_kwargs)
            patch_count = len(logs)
            update_count = len([x for x in logs if x.diffs ])
            unchanged_count = patch_count - update_count
            meta = { 
                API_MSG_RESULT: { 
                    API_MSG_SUBMIT_COUNT : patch_count, 
                    API_MSG_UPDATED: update_count, 
                    API_MSG_UNCHANGED: unchanged_count, 
                    API_MSG_COMMENTS: parent_log.comment
                }
            }
            return self.build_response(
                request, {API_RESULT_META: meta }, response_class=HttpResponse, 
                **kwargs)
            
        if location_data is not None:
            plate_location_fields = ['room', 'freezer', 'shelf', 'bin']
            if (not location_data 
                    or not set(plate_location_fields) | set(location_data.keys())):
                raise ValidationError(
                    key='data', msg='must contain plate location fields')
    
            logger.info('find or create the location: %r', location_data)
            original_location_data = (
                self.get_platelocation_resource()
                    ._get_detail_response_internal(**location_data))
            if not original_location_data:
                logger.info('plate location not found: %r', location_data)
            else:
                location_data['copy_plate_ranges'] = \
                    original_location_data.get('copy_plate_ranges',None)
    
            # Find all of the plates
            original_data = self._get_list_response_internal(**plate_kwargs)
            if not original_data:
                raise Http404
            logger.info('got the plate objects for search data: %d', len(original_data))
            logger.debug('plates for batch edit %r', 
                ['{copy_name}/{plate_number}'.format(**lcp)
                    for lcp in original_data])
            
            logger.info('convert the plates into ranges...')
            library_copy_ranges = {}
            for plate_data in original_data:
                library_copy = '{library_short_name}:{copy_name}'.format(**plate_data)
                plate_range = library_copy_ranges.get(library_copy, [])
                plate_range.append(int(plate_data['plate_number']))
                library_copy_ranges[library_copy] = plate_range
            copy_plate_ranges = []
            for library_copy,plate_range in library_copy_ranges.items():
                if not plate_range:
                    raise ValidationError(
                        key='library_copy', 
                        msg='no plates found for %r' % library_copy)
                plate_range = sorted(plate_range)
                copy_plate_ranges.append(
                    '%s:%s-%s' 
                    % (library_copy, 
                       str(plate_range[0]),
                       str(plate_range[-1]) ))
    
            logger.info(
                'patch original location: %r, with copy_ranges: %r', 
                location_data, copy_plate_ranges) 
            location_copy_ranges = \
                location_data.get('copy_plate_ranges', [])
            location_copy_ranges.extend(copy_plate_ranges)
            location_data['copy_plate_ranges'] = location_copy_ranges
            
            parent_log = self.make_log(request)
            parent_log.key = self._meta.resource_name
            parent_log.uri = self._meta.resource_name
            parent_log.save()
            kwargs['parent_log'] = parent_log
            
            # Consider: create PATCH logs for each of the plates patched as well
            
            response = self.get_platelocation_resource().patch_detail(
                request, 
                data=location_data,
                **kwargs)
            return response
    
    def validate(self, _dict, patch=False, schema=None):
        errors = DbApiResource.validate(self, _dict, patch=patch, schema=schema)
        
        if set(['status','is_active']) | set(_dict.keys()):
            
            if _dict.get('is_active',False) is True:
                if _dict.get('status', None) != 'available':
                    errors['is_active'] = 'Requires status == "available"'
        
        return errors
    
    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        
        logger.debug('patch obj, deserialized: %r', deserialized)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        logger.debug('id_kwargs: %r', id_kwargs)
        required_kwargs_check = ['copy_name', 'plate_number']
        if ( not id_kwargs 
                or not set(required_kwargs_check) & set(id_kwargs.keys())):
            raise ValidationError({
                k:'required' for k in required_kwargs_check if k not in id_kwargs })
        try:
            plate = Plate.objects.get(
                copy__name=id_kwargs['copy_name'],
                plate_number=id_kwargs['plate_number'])
        except ObjectDoesNotExist:
            logger.info('plate does not exist: %r',
                id_kwargs)
            raise
        
        initializer_dict = self.parse(deserialized, schema=schema, create=False)
        errors = self.validate(initializer_dict, schema=schema,patch=True)
        if errors:
            raise ValidationError(errors)
        for key, val in initializer_dict.items():
            if hasattr(plate, key):
                setattr(plate, key, val)
            if key=='status':
                if val == 'available':
                    plate.date_plated = _now()
                    # plate.date_retired = None
                elif val in self.retired_statuses:
                    plate.date_retired = _now() 
        if initializer_dict:
            
            # New 20170407 - is_active check
            if plate.status != 'available':
                plate.is_active = False
            
            plate.save()
            
        # Process plate location fields separately
        plate_location_fields = ['room', 'freezer', 'shelf', 'bin']
        location_data = {}
        update_plate_location = False
        for k in plate_location_fields:
            location_data[k] = deserialized.get(k,None)
            if location_data[k]:
                update_plate_location = True
        if update_plate_location:
            original_location_data = (
                self.get_platelocation_resource()
                    ._get_detail_response_internal(**location_data))
            try:
                plate_location = PlateLocation.objects.get(**location_data)
                logger.info('plate location found: %r', plate_location)
                
            except ObjectDoesNotExist:
                logger.info('plate location not found: %r', location_data)
                logger.info(
                    'plate location does not exist, creating: %r',
                    location_data)
                plate_location = PlateLocation.objects.create(**location_data)
                plate_location.save()
                
            if plate.plate_location != plate_location:
                logger.info('update plate location: %r, to %r', plate,  plate_location)
                plate.plate_location = plate_location
                plate.save()
            
            # FIXME: when using patch/post list, this will create a log for
            # each plate location update in the list, need to move this out to
            # the patch-post_list methods.
            new_location_data = (
                self.get_platelocation_resource()
                    ._get_detail_response_internal(**location_data))
            log = self.get_platelocation_resource().log_patch(
                request, original_location_data, new_location_data,**kwargs)
            if log:
                log.save()
                logger.info('log created; %r', log)
                    
        return { API_RESULT_OBJ: plate }

class UserAgreementResource(DbApiResource):
    
    VOCAB_USER_AGREEMENT_SM = 'sm'
    VOCAB_USER_AGREEMENT_RNAI = 'rnai'
    
    VOCAB_FILE_TYPE_SMUA = 'iccb_l_small_molecule_user_agreement'
    VOCAB_FILE_TYPE_RNAI = 'iccb_l_rnai_user_agreement'
    
    class Meta:

        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'useragreement'
        authorization = UserGroupAuthorization(resource_name)
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
        
        self.screensaveruser_resource = None
        self.attached_file_resource = None
        super(UserAgreementResource, self).__init__(**kwargs)

    def get_su_resource(self):
        if self.screensaveruser_resource is None:
            self.screensaveruser_resource = ScreensaverUserResource()
        return self.screensaveruser_resource
    
    def get_attached_file_resource(self):
        if self.attached_file_resource is None:
            self.attached_file_resource = AttachedFileResource()
        return self.attached_file_resource
    
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
#             url((r"^(?P<resource_name>%s)/" 
#                  r"(?P<user_agreement_id>([\d]+))%s$")
#                     % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<screensaver_user_id>([\d]+))/"
                 r"(?P<type>(\w+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ] 

    @read_authorization
    def get_detail(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)

    @read_authorization
    def get_list(self, request, **kwargs):
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        kwargs['file_type__in'] = [self.VOCAB_FILE_TYPE_RNAI, self.VOCAB_FILE_TYPE_SMUA]
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        manual_field_includes = set(param_hash.get('includes', []))
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, schema)
        filter_expression = \
            self._meta.authorization.filter(request.user,filter_expression)
              
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
         
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
        rowproxy_generator = \
            self._meta.authorization.get_row_property_generator(
                request.user, field_hash, rowproxy_generator)
        
        # specific setup
        _su = self.bridge['screensaver_user']
        _user_cte = ScreensaverUserResource.get_user_cte().cte('ua_user')
        _lab_head = ScreensaverUserResource.get_lab_head_cte('ualh').cte('ulab_head')
        _user_agreement = self.bridge['user_agreement']
        _lab_head_ua = _user_agreement.alias('lab_head_ua')
        _af = self.bridge['attached_file']
        _lhaf = _af.alias('lhaf')
        _vocab = self.bridge['reports_vocabulary']
        
        username = param_hash.pop('username', None)
        screensaver_user_id = param_hash.pop('screensaver_user_id', None)
        
        # Create an "all_users" table that contains one entry for each 
        # agreement type for each user in the system.
        agreement_type_table = (
            select([
                _vocab.c.key.label('type'),
                _vocab.c.ordinal,
                _vocab.c.is_retired.label('is_retired'),
            ])
            .select_from(_vocab)
            .where(_vocab.c.scope == 'useragreement.type')
            .where(func.coalesce(_vocab.c.is_retired,False) != True)
            .order_by(_vocab.c.ordinal))
        types = []
        compiled_stmt = str(agreement_type_table.compile(
            dialect=postgresql.dialect(),
            compile_kwargs={"literal_binds": True}))
        logger.info('user agreement temp table: %s',compiled_stmt)
        with get_engine().connect() as conn:
            result = conn.execute(agreement_type_table)
            types =  [x[0] for x in result ]
        logger.info('types: %r', types)
        queries = []
        for type in types:
            query = (
                select([
                    _su.c.screensaver_user_id,
                    _su.c.username,
                    _su.c.lab_head_id,
                    literal_column("'%s'"%type).label('type')])
                .select_from(_su))
            if username:
                query = query.where(
                    _su.c.username == username)
            if screensaver_user_id:
                query = query.where(
                    _su.c.screensaver_user_id == screensaver_user_id)
            queries.append(query)
        all_users = union_all(*queries).cte('all_users')
      
        # Build a table of entered agreements
        j = _su
        j = j.join(
            _user_agreement, _user_agreement.c.screensaver_user_id 
                == _su.c.screensaver_user_id, isouter=True)
        j = j.join(
            _af, _user_agreement.c.file_id 
                == _af.c.attached_file_id, isouter=True)
        entered_agreements = select([
            _su.c.screensaver_user_id,
            _user_agreement.c.data_sharing_level,
            _user_agreement.c.type,
            _user_agreement.c.date_active,
            _user_agreement.c.date_expired,
            _user_agreement.c.date_notified,
            _af.c.attached_file_id.label('file_id'),
            _af.c.filename,
            ]).select_from(j).cte('entered_agreements')
        _user_entered_agreements = entered_agreements.alias('user_entered_agreements')        
        _lh_entered_agreements = entered_agreements.alias('lh_entered_agreements')        
        
        custom_columns = {
            'screensaver_user_id': all_users.c.screensaver_user_id,
            'username': _user_cte.c.username,
            'user_name': _user_cte.c.name,
            'user_first_name': _user_cte.c.first_name,
            'user_last_name': _user_cte.c.last_name,
            'user_email': _user_cte.c.email,
            'type' : all_users.c.type,
            'data_sharing_level': _user_entered_agreements.c.data_sharing_level,
            'status': case([
                ( _user_entered_agreements.c.date_active != None, 
                  case([
                      (_user_entered_agreements.c.date_expired != None, 'expired')],
                      else_='active'))],
                else_='inactive'),  
            'date_active': _user_entered_agreements.c.date_active,
            'date_notified': _user_entered_agreements.c.date_notified,
            'date_expired': _user_entered_agreements.c.date_expired,
            'file_id': _user_entered_agreements.c.file_id,
            'filename': _user_entered_agreements.c.filename,
            'lab_head_id': _lh_entered_agreements.c.screensaver_user_id,
            'lab_head_name': _lab_head.c.name,
            'lab_affiliation_category': _lab_head.c.lab_affiliation_category,
            'lab_affiliation_name': _lab_head.c.lab_affiliation_name,
            'lab_head_data_sharing_level': _lh_entered_agreements.c.data_sharing_level,
            'lab_head_user_agreement_status': 
                case([
                ( _lh_entered_agreements.c.screensaver_user_id == None, None),
                ( _lh_entered_agreements.c.date_active != None, 
                    case([
                    (_lh_entered_agreements.c.date_expired != None, 'expired')],
                    else_='active'))],
                else_='inactive'),  
            'lab_head_file_id': _lh_entered_agreements.c.file_id,
            'lab_head_filename': _lh_entered_agreements.c.filename,
            }
        
        base_query_tables = [] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
         
        # Build the final query by joining all_users to entered agreements
        j = all_users
        j = j.join(_user_cte, all_users.c.screensaver_user_id
            ==_user_cte.c.screensaver_user_id)
        j = j.join(_user_entered_agreements, 
            and_( all_users.c.screensaver_user_id
                ==_user_entered_agreements.c.screensaver_user_id,
                all_users.c.type==_user_entered_agreements.c.type),
            isouter=True)
        j = j.join(
            _lh_entered_agreements,
            and_(
                all_users.c.lab_head_id != all_users.c.screensaver_user_id,
                all_users.c.lab_head_id == _lh_entered_agreements.c.screensaver_user_id,
                _lh_entered_agreements.c.type== all_users.c.type),
            isouter=True
            )
        j = j.join(_lab_head,
            all_users.c.lab_head_id == _lab_head.c.screensaver_user_id, isouter=True)
        
        stmt = select(columns.values()).select_from(j)
        stmt = stmt.order_by(
            all_users.c.screensaver_user_id, 
            all_users.c.type)

        # general setup

        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        
        # compiled_stmt = str(stmt.compile(
        #     dialect=postgresql.dialect(),
        #     compile_kwargs={"literal_binds": True}))
        # logger.info('compiled_stmt %s', compiled_stmt)
        
        title_function = None
        if use_titles is True:
            def title_function(key):
                return field_hash[key]['title']
        if is_data_interchange:
            title_function = DbApiResource.datainterchange_title_function(
                field_hash,schema['id_attribute'])
        
        return self.stream_response_from_statement(
            request, stmt, count_stmt, filename,
            field_hash=field_hash,
            param_hash=param_hash,
            is_for_detail=is_for_detail,
            rowproxy_generator=rowproxy_generator,
            title_function=title_function, meta=kwargs.get('meta', None),
            use_caching=True)

    @write_authorization
    @un_cache
    @transaction.atomic
    def post_list(self, request, schema=None, **kwargs):
        raise NotImplementedError
    
    @write_authorization
    @un_cache        
    @transaction.atomic
    def post_detail(self, request, **kwargs):
        '''
        Modified POST because:
        - Attached File may be included in request
        '''
        schema = kwargs.get('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        id_attribute = schema['id_attribute']
        
        # Perform manual deserialization; MULTIPART content is in POST dict
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.debug('param_hash: %r', 
            {k:v for k,v in param_hash.items() if k != 'schema'})
        deserialized = {}
        for key in fields.keys():
            if param_hash.get(key, None) is not None:
                deserialized[key] = param_hash[key]
        
        # POST is only used to create or reset to active; file is required
        # Note: see PATCH: status=deactivate
        attached_file = request.FILES.get('attached_file', None)
        if attached_file is None:
            status_update = deserialized.get('status', None)
            if status_update is None or status_update not in ('inactive','expire'):
                raise ValidationError(
                    key='attached_file',
                    msg='required to create or reactivate a user agreement')
        # NOT supporting normal deserialization with POST
        # if attached_file is None:
        #     # If no attached file, allow ordinary deserialization
        #     deserialized = self.deserialize(
        #         request, format=kwargs.get('format', None))
        #     logger.info('deserialized: %r', deserialized)

        # Validate Identity
        username = deserialized.pop('username', None)
        screensaver_user_id = deserialized.pop('screensaver_user_id', None)
        if username is None and screensaver_user_id is None:
            raise NotImplementedError(
                'must provide a screensaver_user_id or username parameter')
        if username is not None and screensaver_user_id is not None:
            raise NotImplementedError(
                'must provide either a screensaver_user_id or username parameter')
        if username is not None:
            su = ScreensaverUser.objects.get(username=username)
            screensaver_user_id = su.screensaver_user_id
        deserialized['screensaver_user_id'] = screensaver_user_id

        id_kwargs = self.get_id(deserialized,validate=True, **kwargs)

        # Create a parent log for the user
        user_resource = self.get_su_resource()
        # Cache the user data for logging
        kwargs_for_user = {
            'exact_fields': ['screensaver_user_id', 'is_active',
                'sm_data_sharing_level', 'rnai_data_sharing_level',
                'lab_head_id'] 
        }
        kwargs_for_user['screensaver_user_id'] = screensaver_user_id
        original_user_data = user_resource._get_detail_response_internal(**kwargs_for_user)
        if not original_user_data:
            msg = 'User not found: %r'
            raise ValidationError({
                'screensaver_user_id': msg % screensaver_user_id,
                'username': msg % username
            })
        parent_log = user_resource.make_log(
            request, attributes=original_user_data, api_action='PATCH')
        parent_log.save()

        log = self.make_log(
            request, id_kwargs, id_attribute=id_attribute, schema=schema)
        log.save()
        original_data = self._get_detail_response_internal(**id_kwargs)

        # POST is only used to create or reset to active
        # Note: see PATCH: status=deactivate, expire
        try:
            user_agreement = UserAgreement.objects.get(**id_kwargs)
            logger.info('UA exists: resetting: %r', user_agreement)
            user_agreement.date_active=_now()
            user_agreement.date_expired = None
            user_agreement.date_notified = None
            user_agreement.save()
        except ObjectDoesNotExist:
            user_agreement = UserAgreement.objects.create(**id_kwargs)
            if 'date_active' not in deserialized:
                deserialized['date_active'] = _now().date().strftime("%Y-%m-%d")
#             user_agreement.date_active=_now()
            logger.info('UA DNE, creating: %r', user_agreement)
        
        # Perform the PATCH
        logger.info('patch_obj: %r', deserialized)
        patch_result = self.patch_obj(request, deserialized, log=log, **kwargs)
        user_agreement = patch_result[API_RESULT_OBJ]
        user_agreement.save()

        # === Attached File ===
        # TODO: delete existing Attached File if replacing
        attached_type = None
        if user_agreement.type == self.VOCAB_USER_AGREEMENT_SM:
            attached_type = self.VOCAB_FILE_TYPE_SMUA
        elif user_agreement.type == self.VOCAB_USER_AGREEMENT_RNAI:
            attached_type = self.VOCAB_FILE_TYPE_RNAI
        else:
            raise ValidationError(
                key='type', 
                msg='must be one of %r' % [
                    self.VOCAB_USER_AGREEMENT_SM, self.VOCAB_USER_AGREEMENT_RNAI])
        kwargs_for_attachedfile = { 'screensaver_user_id': screensaver_user_id }
        kwargs_for_attachedfile['type'] = attached_type
        kwargs_for_attachedfile['filename'] = deserialized.get('filename')
        
        attached_file_resource = self.get_attached_file_resource()
        attached_file_response = \
            attached_file_resource.post_detail(request, **kwargs_for_attachedfile)
        logger.info('UserAgreement attached file result: %r', attached_file_response)
        if attached_file_response.status_code == 201:
            af_data = attached_file_resource._meta.serializer.deserialize(
                LimsSerializer.get_content(attached_file_response), JSON_MIMETYPE)
            logger.info('attached_file: %r', af_data)
            user_agreement.file_id = af_data['attached_file_id']
            user_agreement.save()
        else:
            error_resp = attached_file_resource._meta.serializer.deserialize(
                LimsSerializer.get_content(attached_file_response), JSON_MIMETYPE)
            logger.error('attached file error: %r, %r', 
                attached_file_response.status_code, error_resp)
            raise Exception('attached file resource error: %r', error_resp)

        # === Attached File - Done ===

        if user_agreement.date_active is not None:
            if user_agreement.date_expired is None:

                logger.info('UserAgreement is active: %r, add login access for user')
                # NOTE: is_active is a reports.user property and will only be 
                # set if the user has a ecommonsId or username set
                user_schema = user_resource.build_schema(request.user)
                user_patch_data = {
                    'screensaver_user_id': screensaver_user_id,
                    'is_active': True }
                patch_result = user_resource.patch_obj(
                    request, user_patch_data, schema=user_schema)
                logger.info('user patch result: %r', patch_result)

        
        new_user_data = user_resource._get_detail_response_internal(**kwargs_for_user)
        user_resource.log_patch(request, original_user_data, new_user_data, parent_log)
        parent_log.save()
        
        new_data = self._get_detail_response_internal(**id_kwargs)
        
        self.log_patch(request, original_data, new_data, log, 
            id_attribute=id_attribute, parent_log=parent_log, full_create_log=True)
        log.save()
        
        data = { API_RESULT_DATA: [new_data]}
        return self.build_response(request, data)
    
    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):

        logger.info('patch user agreement item: %r', deserialized)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        id_attribute = schema['id_attribute']

        id_kwargs = self.get_id(deserialized,validate=True,schema=schema,**kwargs)
        
        screensaver_user_id = id_kwargs['screensaver_user_id']
        screensaver_user = ScreensaverUser.objects.get(
            screensaver_user_id=screensaver_user_id)

        # PATCH may not be used to create a user agreement; see post_detail
        user_agreement = UserAgreement.objects.get(**id_kwargs)

        initializer_dict = self.parse(deserialized, schema=schema, create=True)
        logger.info('initializer: %r', initializer_dict)
        
        if not initializer_dict:
            raise Exception('Empty patch')
        
        errors = self.validate(initializer_dict, schema=schema)
        if errors:
            raise ValidationError(errors)
        
        new_status = initializer_dict.pop('status', None)
        if new_status is not None:
            if initializer_dict:
                logger.error(
                    'Other actions not allowed if expiring or deactivating: %r', 
                    initializer_dict)
                raise ValidationError(
                    key='status', 
                    msg='other values may not be set if expiring or deactivating')
            if user_agreement.date_active is None:
                raise ValidationError(
                    key='status', msg='Agreement is not active')
        else:
            if user_agreement.date_expired is not None:
                msg = 'Expired agreements may not be patched'
                raise ValidationError({
                    k:msg for k in initializer_dict.keys() })

        for key, val in initializer_dict.items():
            if hasattr(user_agreement, key):
                setattr(user_agreement, key, val)

        if user_agreement.date_active is None:
            raise ValidationError(
                key='date_active', msg='required')
        
        if user_agreement.data_sharing_level is None:
            raise ValidationError(
                key='data_sharing_level', msg='required')
        
        meta = {}
        if new_status is not None:
            if new_status == 'expired':
                if user_agreement.date_expired is not None:
                    raise ValidationError(
                        key='status', 
                        msg='Agreement may not be expired more than once')
                user_agreement.date_expired = _now()
                user_agreement.save()
                logger.info('UserAgreement %r expired', user_agreement)
#                 return { API_RESULT_OBJ: user_agreement }
            elif new_status == 'inactive':
                if user_agreement.date_expired is not None:
                    raise ValidationError(
                        key='status', 
                        msg='Agreement is already expired and may not be deactivated')
                user_agreement.date_active = None
                user_agreement.date_notified = None
                user_agreement.data_sharing_level = None
                user_agreement.file = None
                user_agreement.save()
                logger.info('UserAgreement %r deactivated', user_agreement)
#                 return { API_RESULT_OBJ: user_agreement }
            else:
                raise ValidationError(
                    key='status', msg='may only be set to "expired" or "inactive"')
        else:
            # Disallow setting DSL to a level different than the Lab Head
            lab_head_id = screensaver_user.lab_head_id
            if lab_head_id != screensaver_user.screensaver_user_id:
                try:
                    lh_user_agreement = UserAgreement.objects.get(
                        screensaver_user_id=lab_head_id, type=user_agreement.type)
                    if user_agreement.data_sharing_level is not None:
                        if user_agreement.data_sharing_level != \
                            lh_user_agreement.data_sharing_level:
                            raise ValidationError(
                                key='data_sharing_level',
                                msg='Must match Lab Head value: %r' 
                                    % lh_user_agreement.data_sharing_level)    
                except ObjectDoesNotExist:
                    logger.info('no lab head user agreement found')
                    raise ValidationError(
                        key='status',
                        msg='User may not be active until their lab head has a user agreement')
            # TODO: if user is a LabHead:
            # Display warning if Lab Head DSL is being set to something other than
            # Lab_member's DSL
            # TODO: if DSL does not match user screens, display warning?
            
            user_agreement.save()
        
        active_user_agreements = UserAgreement.objects\
            .filter(screensaver_user=screensaver_user)\
            .filter(date_active__isnull=False)\
            .filter(date_expired__isnull=True)
        # If all user agreements are inactive, turn off login capability
        if not active_user_agreements.exists():
            logger.info(
                'no active user agreements exist for %r, removing is_active', 
                screensaver_user)
            user_resource = self.get_su_resource()
            user_schema = user_resource.build_schema(request.user)
            user_patch_data = {
                'screensaver_user_id': screensaver_user_id,
                'is_active': False }
            patch_result = user_resource.patch_obj(
                request, user_patch_data, schema=user_schema)
            logger.info('user patch result: %r', patch_result)
        
        logger.info('UserAgreement %r patched', user_agreement)

        return { 
            API_RESULT_OBJ: user_agreement,
            API_RESULT_META: meta }

class ScreenAuthorization(UserGroupAuthorization):
    
    VOCAB_USER_AGREEMENT_SM = UserAgreementResource.VOCAB_USER_AGREEMENT_SM
    VOCAB_USER_AGREEMENT_RNAI = UserAgreementResource.VOCAB_USER_AGREEMENT_RNAI
    
    '''
    Implements the data viewing restrictions defined in the ICCB-L User 
    Agreement documents (Small Molecule and RNAi).
    
    Restriction mechanisms:
    1. Schema level: search and sort capability are removed for fields with 
    a Field.data_sharing_level of 1 or 2 (for restricted users).
    2. Query level: restricted users may only see Screens that they are 
    authorized to view per User Agreement viewing rules.
    3. Cursor Generator: restricted users may only see screen 
    properties having a Field.data_sharing_level matching the users access 
    level for the specific Screen.
    '''

    def __init__(self, *args, **kwargs):
        UserGroupAuthorization.__init__(self, *args, **kwargs)
        self._screen_overlap_table = None
        
#     def is_restricted_view(self, user):
#         return super(ScreenAuthorization, self).is_restricted_view(user)

    def _is_resource_authorized(
        self, user, permission_type, **kwargs):
        '''
        Override UserGroupAuthorization to give users permission for restricted
        read access to screens they are authorized to view.
        '''
        
        authorized = super(ScreenAuthorization, self)._is_resource_authorized(
            user, permission_type, **kwargs)
        if authorized is True:
            return True
        else:
            screensaver_user = ScreensaverUser.objects.get(username=user.username)
            authorized_screens = self.get_read_authorized_screens(screensaver_user)
            
            if authorized_screens:
                logger.info('_is_resource_authorized: %r - True', user)
                logger.debug(
                    'user: %r, authorized_screens: %r', 
                    screensaver_user.username, 
                    [x.facility_id for x in authorized_screens])
                return True
        
        return False

   
    def get_screen_overlapping_table(self):
        '''
        Create a current screen-overlap table from the database:
        - Fields: ['facility_id', 'data_sharing_level','screen_type',
        'overlapping_positive_screens']
        '''
        logger.info('get_screen_overlapping_table...')
        bridge = get_tables()
        _screen = bridge['screen']
        _overlap_screen = _screen.alias('overlap_screen')
        _screen_overlap = DbApiResource.get_create_screen_overlap_indexes()
        
        stmt = (
            select([
                _screen.c.facility_id,
                _screen.c.data_sharing_level,
                _screen.c.screen_type,
                func.array_agg(_overlap_screen.c.facility_id)
            ])
            .select_from(
                _screen
                    .join(
                        _screen_overlap, 
                        _screen_overlap.c.screen_id==_screen.c.screen_id, isouter=True)
                    .join(
                        _overlap_screen, 
                        _overlap_screen.c.screen_id
                            ==_screen_overlap.c.overlap_screen_id, isouter=True))
            .group_by(_screen.c.facility_id, _screen.c.data_sharing_level, 
                _screen.c.screen_type )
        )
        
        with get_engine().begin() as conn:
            _dict = {}
            rows = conn.execute(stmt).fetchall()
            for row in rows:
                _dict[row[0]] = dict(zip(
                    ['facility_id', 'data_sharing_level','screen_type',
                        'overlapping_positive_screens'],
                    row[0:]))
        return _dict
    
    def get_user_screens(self, screensaver_user):
        '''
        Return the screens that the user is a member of.
        '''
        if DEBUG_SCREEN_ACCESS:
            logger.info('get_user_screens %r', screensaver_user)
        my_screens = Screen.objects.all().filter(
            Q(lead_screener=screensaver_user)
            | Q(lab_head=screensaver_user)
            | Q(collaborators=screensaver_user))
        if DEBUG_SCREEN_ACCESS:
            logger.info('user: %r, screens: %r', 
                screensaver_user.username, [s.facility_id for s in my_screens])
        return set(my_screens)
    
    def get_user_data_sharing_level(self, screensaver_user, user_agreement_type):
        logger.info('get user data sharing agreements for %r, %r', 
            screensaver_user, user_agreement_type)
        active_agreements = \
            screensaver_user.useragreement_set.all()\
            .filter(type=user_agreement_type)\
            .filter(date_active__isnull=False)\
            .filter(date_expired__isnull=True)
        current_dsl = 0
        if active_agreements.exists():
            current_dsl = active_agreements[0].data_sharing_level
        else:
            logger.info('no active %r user agreements for user: %r', type, screensaver_user)
        logger.debug('user dsl: %r, %r, %r', 
            screensaver_user, user_agreement_type, current_dsl)
        return current_dsl
    
    def has_sm_data_deposited(self, screensaver_user):
        '''
        True if the user has qualifying data deposited
        - as defined in the User Agreement: level 1 or 2 data, depending on the 
        data_sharing_level for the user.
        '''
        current_dsl = self.get_user_data_sharing_level(
            screensaver_user, self.VOCAB_USER_AGREEMENT_SM)
        my_screens = self.get_user_screens(screensaver_user)
        for screen in my_screens:
            if hasattr(screen, 'screenresult'):
                if ( screen.screen_type == VOCAB_SCREEN_TYPE_SM
                     and screen.data_sharing_level < 3 ):
                    if screen.data_sharing_level == current_dsl:
                        logger.info(
                            'has_sm_data_deposited %r: True, dsl: %r', 
                            screensaver_user.screensaver_user_id, current_dsl)
                        return True
        if DEBUG_SCREEN_ACCESS:
            logger.info('has_sm_data_deposited %r: %r', 
                screensaver_user.screensaver_user_id, current_dsl)
        return False
    
    def has_rna_data_deposited(self, screensaver_user):
        '''
        True if the user has qualifying data deposited
        - as defined in the User Agreement: level 1 or 2 data, depending on the 
        data_sharing_level for the user.
        '''
        current_dsl = self.get_user_data_sharing_level(
            screensaver_user, self.VOCAB_USER_AGREEMENT_RNAI)
        my_screens = self.get_user_screens(screensaver_user)
        for screen in my_screens:
            if hasattr(screen, 'screenresult'):
                if ( screen.screen_type == VOCAB_SCREEN_TYPE_RNAI
                     and screen.data_sharing_level < 3 ):
                    if screen.data_sharing_level == current_dsl:
                        logger.info(
                            'has_rnai_data_deposited %r: True, dsl: %r', 
                            screensaver_user.screensaver_user_id, current_dsl)
                        return True
        if DEBUG_SCREEN_ACCESS:
            logger.info('has_rnai_data_deposited %r: %r', 
                screensaver_user.screensaver_user_id, current_dsl)
        return False
    
    def get_read_authorized_screens(self, screensaver_user):
        '''
        Returns the screens that the user has permission to view general details
        for.
        - either unrestricted "read" access or restricted access as defined in
        the User Agreement rules for data visibility. 
        '''
        
        if DEBUG_SCREEN_ACCESS:
            logger.info('get_read_authorized_screens %r', screensaver_user)
        authorized_screens = set()
        authorized_screens.update(self.get_user_screens(screensaver_user))
        public_screens = Screen.objects.all().filter(data_sharing_level=0)
        authorized_screens.update(public_screens)

        has_sm_data_deposited = self.has_sm_data_deposited(screensaver_user)
        has_rna_data_deposited = self.has_rna_data_deposited(screensaver_user)
        
        if has_sm_data_deposited:
            visible_screens = (
                Screen.objects.all().filter(screen_type=VOCAB_SCREEN_TYPE_SM)
                    .filter(data_sharing_level__lt=3))
            authorized_screens.update(visible_screens)
        if has_rna_data_deposited:
            visible_screens = (
                Screen.objects.all().filter(screen_type=VOCAB_SCREEN_TYPE_RNAI)
                    .filter(data_sharing_level__lt=3))
            authorized_screens.update(visible_screens)
        if DEBUG_SCREEN_ACCESS:
            logger.info(
                'user: %r, visible screens (dsl <3): %r', 
                screensaver_user.screensaver_user_id, 
                [x.facility_id for x in authorized_screens])
        return authorized_screens
    
    def has_screen_read_authorization(self, user, screen_facility_id):
        '''
        True if the user has permission to view general details for the resource;
        - User-Screen Access Level = 0, 1, 2, 3
        - either unrestricted "read" access or restricted access as defined in
        the User Agreement rules for data visibility. 
        '''
        logger.info('has_screen_read_authorization: %r, %r', user, screen_facility_id)

        is_restricted = self.is_restricted_view(user)
        if is_restricted is not True:
            return True
        
        screensaver_user = ScreensaverUser.objects.get(username=user.username)
        authorized_screens = self.get_read_authorized_screens(screensaver_user)
        result = screen_facility_id in set([screen.facility_id for screen in authorized_screens])
        logger.info('has_screen_read_authorization: %r, %r: %r', 
            user, screen_facility_id, result)
        return result
    
    def get_user_effective_data_sharing_level(self, screensaver_user, screen_type):
        effective_dsl = 0
        if screen_type==VOCAB_SCREEN_TYPE_RNAI \
            and self.has_rna_data_deposited(screensaver_user) is True:
            effective_dsl = self.get_user_data_sharing_level(
                screensaver_user, self.VOCAB_USER_AGREEMENT_RNAI)
        if screen_type==VOCAB_SCREEN_TYPE_SM \
            and self.has_sm_data_deposited(screensaver_user) is True:
            effective_dsl =  self.get_user_data_sharing_level(
                screensaver_user, self.VOCAB_USER_AGREEMENT_SM)
        if DEBUG_SCREEN_ACCESS:
            logger.info('effective dsl: %r, %r, %r', screensaver_user, screen_type, effective_dsl)
        return effective_dsl
        
        
#         if screensaver_user.useragreement_set.filter(type=screen_type).exists():
#             user_agreement = screensaver_user.useragreement_set.filter(type=screen_type)[0]
#             if user_agreement.date_expired is None:
#                 effective_dsl = user_agreement.data_sharing_level
# 
#         return effective_dsl
    
    def get_screen_access_level_table(self, username):
        '''
        Create a current user-screen access level table:
        - Fields: 
            'facility_id', 'data_sharing_level','screen_type',
            overlapping_positive_screens: 
                Calculated in two passess, filtered to only show screens that the 
                user has user_access_level_granted >= 2 from the unfiltered set
                of overlapping_positive_screens
            user_access_level_granted:
            Effective User-Screen Access Level - fields visible:
             
            0 - Field level 0 only; no overlapping data
            1 - (overlapping pos screens) Field level 0,1; no
                overlapping data (unless viewing from own
                screen results, not visible here)
            2 - (mutual shared screens) Field level 0,1,2, Screen Results, 
                Positives Summary; 
                overlapping_positive_screens: filtered to show users 
                user_access_level_granted 2,3 screens
            3 - (own screens) Field level 0,1,2,3, 
                Screen Results, CPRs, Activities, Visits; 
                overlapping_positive_screens: filtered to show users 
                user_access_level_granted 2,3 screens
        '''
        
        screen_overlapping_table = self.get_screen_overlapping_table()
        screensaver_user = ScreensaverUser.objects.get(username=username)
        my_screen_facility_ids = set([
            screen.facility_id for screen 
                in self.get_user_screens(screensaver_user)])
        
        if DEBUG_SCREEN_ACCESS:
            logger.info('my screens: %r, %r', username, my_screen_facility_ids)

        user_effective_sm_dsl = self.get_user_effective_data_sharing_level(
            screensaver_user, VOCAB_SCREEN_TYPE_SM)
        user_effective_rna_dsl = self.get_user_effective_data_sharing_level(
            screensaver_user, VOCAB_SCREEN_TYPE_RNAI)
        if DEBUG_SCREEN_ACCESS:
            logger.info('my dsls: %r, %r, %r', username, 
                user_effective_sm_dsl,user_effective_rna_dsl)
        
        my_qualified_sm_facility_ids = set([
            facility_id for facility_id in my_screen_facility_ids if
                screen_overlapping_table[facility_id]['data_sharing_level']
                    == user_effective_sm_dsl and
                screen_overlapping_table[facility_id]['screen_type'] == VOCAB_SCREEN_TYPE_SM            
            ])
        my_qualified_rnai_facility_ids = set([
            facility_id for facility_id in my_screen_facility_ids if
                screen_overlapping_table[facility_id]['data_sharing_level']
                    == user_effective_rna_dsl and 
                screen_overlapping_table[facility_id]['screen_type'] == VOCAB_SCREEN_TYPE_RNAI            
            ])

        
        if DEBUG_SCREEN_ACCESS:
            logger.info(
                'user: %s, '
                'user_effective_rna_dsl: %r, ' 
                'user_effective_sm_dsl: %r',
                username, user_effective_rna_dsl, user_effective_sm_dsl)
            
        def effective_access_level_function(
                facility_id=None, data_sharing_level=None, 
                screen_type=None, overlapping_positive_screens=None ):

            if screen_type == VOCAB_SCREEN_TYPE_RNAI:
                user_effective_dsl = user_effective_rna_dsl
                my_qualified_facility_ids = my_qualified_rnai_facility_ids
            elif screen_type == VOCAB_SCREEN_TYPE_SM:
                user_effective_dsl = user_effective_sm_dsl
                my_qualified_facility_ids = my_qualified_sm_facility_ids
            else:
                raise ProgrammingError('unknown screen_type: %r', screen_type)

            effective_access_level = None
            if facility_id in my_screen_facility_ids:
                effective_access_level = 3
            elif data_sharing_level == 0:
                effective_access_level = 2
            elif data_sharing_level == 3:
                # Note: level 3 screens other than own should not appear
                effective_access_level = None
            elif user_effective_dsl == 1 and data_sharing_level == 1:
                effective_access_level = 2
            elif user_effective_dsl in [1,2] and data_sharing_level in [1,2]:
                if  my_qualified_facility_ids & set(overlapping_positive_screens):
                    effective_access_level = 1
                else:
                    effective_access_level = 0
                    
            return effective_access_level
        
        # pass one: calculate the user_access_level_granted
        for facility_id,_dict in screen_overlapping_table.items():
            _dict['user_access_level_granted'] = \
                effective_access_level_function(**_dict)
        
        screens_by_user_access_level = defaultdict(set)
        for facility_id, screen in screen_overlapping_table.items():
            ual = screen['user_access_level_granted']
            screens_by_user_access_level[ual].add(facility_id)
        
        if DEBUG_SCREEN_ACCESS: 
            logger.info('screens_by_user_access_level: %r ',
                screens_by_user_access_level)
            logger.info('screen access levels for "%s" %r',
                username, [(k,len(v)) for k,v 
                    in screens_by_user_access_level.items()])

        # pass two: calculate overlapping screens visible
        overlapping_screens_visible_2 = (
            screens_by_user_access_level[2]
            | screens_by_user_access_level[3])
        overlapping_screens_visible_3 = (
            screens_by_user_access_level[1]
            | screens_by_user_access_level[2]
            | screens_by_user_access_level[3])
        for facility_id,_dict in screen_overlapping_table.items():
            if _dict['user_access_level_granted'] == 2:
                overlapping_positive_screens = \
                    set(_dict['overlapping_positive_screens'])
                overlapping_positive_screens &= overlapping_screens_visible_2
                _dict['overlapping_positive_screens'] = \
                    sorted(overlapping_positive_screens)
            elif _dict['user_access_level_granted'] == 3:
                overlapping_positive_screens = \
                    set(_dict['overlapping_positive_screens'])
                overlapping_positive_screens &= overlapping_screens_visible_3
                _dict['overlapping_positive_screens'] = \
                    sorted(overlapping_positive_screens)
            else:
                _dict['overlapping_positive_screens'] = None
        
        return screen_overlapping_table
    
    def filter(self, user, filter_expression):
        if self.is_restricted_view(user):
            if DEBUG_AUTHORIZATION:
                logger.info('create_authorized_screen_filter for %r', user)
            screensaver_user = ScreensaverUser.objects.get(username=user.username)
            authorized_screens = \
                self.get_read_authorized_screens(screensaver_user)
            # FIXME: should use screen_id, but it is not visible at the top level?
            auth_filter = column('facility_id').in_(
                [screen.facility_id for screen in authorized_screens])
            if filter_expression is not None:
                filter_expression = and_(filter_expression, auth_filter)
            else:
                filter_expression = auth_filter
 
        return filter_expression

    def filter_in_sql(self, user, stmt, screen_table):
        return stmt
        # if self.is_restricted_view(user):
        #     logger.info('create_authorized_screen_filter')
        #     screensaver_user = ScreensaverUser.objects.get(username=user.username)
        #     authorized_screens = \
        #         self.get_read_authorized_screens(screensaver_user)
        #     stmt = stmt.where(
        #         screen_table.c.screen_id.in_(
        #             [screen.screen_id for screen in authorized_screens]))
        # return stmt
    
    def get_row_property_generator(self, user, fields, extant_generator):
        '''
        Filter result properties based on authorization rules
        '''
        return self.get_access_level_property_generator(
           user, fields, extant_generator)
    
    def get_access_level_property_generator(self, user, fields, extant_generator):
        '''
        Override UserGroupAuthorization to filter based on User Agreement data
        sharing rules.
        '''
        
        is_restricted = self.is_restricted_view(user)
        if is_restricted is not True:
            return extant_generator
        else:
            logger.info('get_access_level_property_generator: %r for user: %r', 
                self.resource_name, user)
            
            screen_access_dict = self.get_screen_access_level_table(user.username)
            fields_by_level = self.get_fields_by_level(fields)
            logger.debug('fields by level: %r', fields_by_level)
            class Row:
                def __init__(self, row):
                    logger.debug(
                        'filter screen row: %r', 
                        [(key, row[key]) for key in row.keys()])
                    self.row = row
                    self.facility_id = facility_id = row['facility_id']
                    
                    effective_access_level = None
                    screen_data = screen_access_dict.get(facility_id, None)
                    if screen_data:
                        effective_access_level = \
                            screen_data.get('user_access_level_granted',None)
                    if DEBUG_SCREEN_ACCESS:
                        logger.info('screen: %r, effective_access_level: %r', 
                            facility_id, effective_access_level)
                    self.effective_access_level = effective_access_level
                    self.allowed_fields = set()
                    if effective_access_level is not None:
                        for level in range(0,effective_access_level+1):
                            self.allowed_fields.update(fields_by_level[level])
                            if DEBUG_SCREEN_ACCESS:
                                logger.info(
                                    'allow level: %r: %r', 
                                    level, fields_by_level[level])
                    else: 
                        logger.warn('user: %r effective_access_level is None, screen: %r',
                            user.username, facility_id)
                    if DEBUG_SCREEN_ACCESS:
                        logger.info('allowed fields: %r', self.allowed_fields)
                def has_key(self, key):
                    if key == 'user_access_level_granted':
                        return True
                    return self.row.has_key(key)
                def keys(self):
                    return self.row.keys();
                def __getitem__(self, key):
                    logger.debug(
                        'key: %r, allowed: %r', key, key in self.allowed_fields)
                    if key == 'user_access_level_granted':
                        return self.effective_access_level
                    if self.row[key] is None:
                        return None
                    else:
                        if key in self.allowed_fields:
                            if key == 'has_screen_result':
                                original_val = self.row[key]
                                val = original_val
                                if self.effective_access_level >= 2:
                                    # Change not_shared to available for sharing users
                                    if val == 2: # "not shared"
                                        val = 1 # "available"
                                else:
                                    if val == 1: # available
                                        val = 2 # not shared
                                return val
                            elif key == 'overlapping_positive_screens':
                                if self.effective_access_level >= 2:
                                    reference_screen = \
                                        screen_access_dict[self.facility_id]
                                    return reference_screen['overlapping_positive_screens']
                                # restricted users may not view
                                return None
                            else:
                                logger.debug('allow %r: %r', key, self.row[key])
                                return self.row[key]
                        else:
                            logger.debug(
                                '%r filter field: %r for restricted user: %r',
                                self.facility_id, key, user.username)
                            return None

            def screen_property_generator(cursor):
                if extant_generator is not None:
                    cursor = extant_generator(cursor)
                for row in cursor:
                    yield Row(row)
            
            return screen_property_generator
        

class ScreenResultAuthorization(ScreenAuthorization):
    
    def _is_resource_authorized(
        self, user, permission_type, screen_facility_id=None, **kwargs):
        '''
        Override UserGroupAuthorization to determine if the user is allowed 
        to view the screen results (level 2 or 3 access)
        '''
        authorized = super(ScreenAuthorization, self)._is_resource_authorized(
            user, permission_type, **kwargs)
        if authorized is True:
            return True
        else:
            if not screen_facility_id:
                raise NotImplementedError(
                    'must provide a screen_facility_id parameter')
            try:
                screen = Screen.objects.get(facility_id=screen_facility_id)

                screensaver_user = ScreensaverUser.objects.get(username=user.username)
                users_own_screens = self.get_user_screens(screensaver_user)
                user_effective_dsl = self.get_user_effective_data_sharing_level(
                    screensaver_user, screen.screen_type)
                can_view = False
                if screen in users_own_screens:
                    can_view = True
                elif user_effective_dsl in [0,1,2,3]:
                    if screen.data_sharing_level == 0:
                        can_view = True
                    elif screen.data_sharing_level == 1:
                        if  user_effective_dsl == 1:
                            can_view = True
                if DEBUG_SCREEN_ACCESS:
                    logger.info(
                        'user: %s, effective dsl: %r, screen: %r, type: %r, dsl: %r',
                        user.username, user_effective_dsl, screen_facility_id, 
                        screen.screen_type, screen.data_sharing_level)
                if can_view is not True:
                    raise PermissionDenied
                if not hasattr(screen, 'screenresult'):
                    raise Http404('No screen result for: %r'%screen_facility_id)
                return True
            
            except ObjectDoesNotExist:
                raise Http404
        return False

    def get_row_property_generator(self, user, fields, extant_generator):
        '''
        Filter result properties based on authorization rules
        '''
        return self.get_access_level_property_generator(
           user, fields, extant_generator)
    
    def get_access_level_property_generator(self, user, fields, extant_generator):
        
        is_restricted = self.is_restricted_view(user)
        if is_restricted is not True:
            return extant_generator
        else:
            
            access_level_1_fields = [ field['key'] 
                for field in fields.values() 
                    if field.get('user_access_level_granted') == 1]
            datacolumn_fields = [field for field in fields.values() 
                if field.get('is_datacolumn',False) is True ]
            access_level_1_cp = [field['key']
                for field in datacolumn_fields
                    if field['vocabulary_scope_ref'] 
                        == 'resultvalue.confirmed_positive_indicator']
            access_level_1_pp = [field['key']
                for field in datacolumn_fields
                    if field['vocabulary_scope_ref'] 
                        == 'resultvalue.partitioned_positive']
            access_level_1_boolean_positive = [field['key']
                for field in datacolumn_fields
                    if field['data_type'] == 'boolean']
            if DEBUG_AUTHORIZATION:
                logger.info('access level 1 confirmed positive: %r',
                    access_level_1_cp )
                logger.info('access level 1 partitioned positive: %r',
                    access_level_1_pp )
                logger.info('access level 1 boolean positive: %r',
                    access_level_1_boolean_positive )
            
            class Row:
                def __init__(self, row):
                    self.row = row
                    self.well_id = row['well_id']
                    self.is_positive = False
                    if row.has_key('is_positive'):
                        self.is_positive = row['is_positive']

                def has_key(self, key):
                    return self.row.has_key(key)
                def keys(self):
                    return self.row.keys();
                def __getitem__(self, key):
                    if self.row[key] is None:
                        return None
                    
                    elif key in access_level_1_fields:
                        value = self.row[key]
                        logger.info('level 1: %r:%r:%r',
                            self.well_id, key, value)
                        if self.is_positive is not True:
                            return None
                        final_value = None
                        if key in access_level_1_boolean_positive:
                            if value is True:
                                final_value = value
                        elif key in access_level_1_cp:
                            if int(value) == 3:
                                final_value = value
                        elif key in access_level_1_pp:
                            if int(value) != 0:
                                final_value = value
                        else:
                            logger.error(
                                'unknown level 1 result value field: %r:%r: %r',
                                self.well_id, key, value)
                        return final_value
                    else:
                        # NOTE: all users viewing screen results already have 
                        # granted access level 2 or 3 for this screen result
                        return self.row[key]

            def result_value_generator(cursor):
                if extant_generator is not None:
                    cursor = extant_generator(cursor)
                for row in cursor:
                    yield Row(row)
            
            return result_value_generator


class ScreenResultResource(DbApiResource):

    class Meta:
    
        queryset = ScreenResult.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'screenresult'
        authorization = ScreenResultAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = ScreenResultSerializer()
        object_class = dict
        max_limit = 10000
        
    def __init__(self, **kwargs):

        self.scope = 'fields.screenresult'
        super(ScreenResultResource, self).__init__(**kwargs)
        
        self.reagent_resource = None
        self.screen_resource = None
        self.datacolumn_resource = None
        
    def get_datacolumn_resource(self):
        if self.datacolumn_resource is None:
            self.datacolumn_resource = DataColumnResource()
        return self.datacolumn_resource

    def get_reagent_resource(self):
        if not self.reagent_resource:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource
    
    def get_screen_resource(self):
        if not self.screen_resource:
            self.screen_resource = ScreenResource()
        return self.screen_resource
    
    @transaction.atomic()        
    def clear_cache(self, all=False, by_date=None, by_uri=None, by_size=False):
        
        logger.info(
            'clear_cache called: all: %r, by_date: %r, by_uri: %r, by_size: %r',
            all, by_date, by_uri, by_size)

        DbApiResource.clear_cache(self)

        max_indexes_to_cache = getattr(
            settings, 'MAX_WELL_INDEXES_TO_CACHE', 2e+08)
        logger.debug('max_indexes_to_cache %s' % max_indexes_to_cache)

        _wellQueryIndex = self.bridge['well_query_index']

        try:
            query = CachedQuery.objects.filter(uri__contains='/screenresult/')
            if query.exists():
                logger.debug(
                    'clear_cache: screenresult queries to consider %s',
                    [(x.id, x.uri) for x in query])
            ids = set()
            if by_size:
                query = query.order_by('-datetime')
                cumulative_count = 0

                for q in query:
                    if q.count:
                        cumulative_count += q.count
                    if cumulative_count >= max_indexes_to_cache:
                        last_id_to_save = q.id
                        logger.info(
                            'cumulative_count: %d, last_id_to_save: %d',
                            cumulative_count, last_id_to_save)
                        query = query.filter(id__lte=last_id_to_save)
                        ids.update([q.id for q in query])
                        break
                
            if by_date:  # TODO: test
                query = query.filter(datetime__lte=by_date)
                ids.update([q.id for q in query])
            if by_uri: 
                query = query.filter(uri__exact=by_uri)
                ids.update([q.id for q in query])
                
            if ids or all:
                logger.info('clear cachedQueries: ids: %r, all: %r', ids, all)
                if all:
                    stmt = delete(_wellQueryIndex)
                    get_engine().execute(stmt)
                    logger.info('cleared all cached wellQueryIndexes')
                    CachedQuery.objects.all().delete()
                    # # TODO: delete the well_data_column_positive_indexes
                    # # related to this uri only
                    # stmt = delete(_well_data_column_positive_index)
                    # conn.execute(stmt)
                    # logger.info(
                    #     'cleared all cached well_data_column_positive_indexes')
                    # TODO: delete the screen_overlap references to this URI only
                else:
                    stmt = delete(_wellQueryIndex).where(
                        _wellQueryIndex.c.query_id.in_(ids))
                    get_engine().execute(stmt)
                    logger.info('cleared cached wellQueryIndexes: %r', ids)
                    CachedQuery.objects.filter(id__in=ids).delete()

                    # stmt = delete(_well_data_column_positive_index)
                    # conn.execute(stmt)
                    # logger.info(
                    #     'cleared all cached well_data_column_positive_indexes')
            else:
                logger.info('no CachedQuery values to be cleared')

        except Exception, e:
            logger.exception('on screenresult clear cache')
            raise e  

        # NOTE: for now, clear the wdcpi when a new screen_result is loaded:
        # (when "all" or "by_uri" is set)
        # TODO: incrementally clear the wdcpi for each set of ids, or all
        if all is True or by_uri is not None:
            logger.info('all: %r, by_uri: %r', all, by_uri)
            get_engine().execute(delete(self.get_table_def('well_data_column_positive_index')))
            logger.info(
                'cleared all cached well_data_column_positive_indexes')
            get_engine().execute(delete(self.get_table_def('screen_overlap')))
            logger.info(
                'cleared all cached screen_overlap entries')
        else:
            logger.info('do not clear the well_data_column_positive_index table...')
        
        # Manually clear the screen caches
        DbApiResource.clear_cache(self)
        caches['screen'].clear()
        # self.get_screen_resource().clear_cache()
            
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/"
                r"(?P<screen_facility_id>\w+)/" 
                r"(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/"
                r"(?P<screen_facility_id>\w+)/"
                r"schema%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/"
                r"(?P<screen_facility_id>\w+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]

    @read_authorization
    def get_detail(self, request, **kwargs):

        facility_id = kwargs.get('screen_facility_id', None)
        if not facility_id:
            raise NotImplementedError(
                'must provide a screen_facility_id parameter')

        well_id = kwargs.get('well_id', None)
        if not well_id:
            raise NotImplementedError('must provide a well_id parameter')

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):
        logger.info('get_list: %r', kwargs)
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)
        
    def _build_result_value_column(self, field_information):
        '''
        Each result value will be added to the query as a subquery select:
        (SELECT <value_field> 
         FROM result_value
         WHERE result_value.data_column_id=<id> 
         AND rv.well_id=assay_well.well_id limit 1) as <data_column_value_alias>
        '''
        key = field_information['key']
        data_column_id = field_information['data_column_id']
        data_column_type = field_information.get('data_type') 
        # TODO: column to select: use controlled vocabulary
        column_to_select = None
        if(data_column_type in ['numeric', 'decimal', 'integer']):  
            column_to_select = 'numeric_value'
        else:
            column_to_select = 'value'
        logger.debug('_build_result_value_column: %r, %r, %r, %r', 
            key, column_to_select, data_column_id, column_to_select)
        
        _rv = self.bridge['result_value']
        rv_select = select([column(column_to_select)]).select_from(_rv)
        rv_select = rv_select.where(
            _rv.c.data_column_id == field_information['data_column_id'])
        rv_select = rv_select.where(_rv.c.well_id == text('assay_well.well_id'))
        # FIXME: limit to rv's to 1 - due to the duplicated result value issue
        rv_select = rv_select.limit(1)  
        rv_select = rv_select.label(key)
        return rv_select

    def _build_result_value_cte(self, field_information):
        '''
        Not used - an alternate method of constructing the result values as 
        subqueries. Analyze indicates that this requires more memory and 
        scales poorly.
        '''
        _rv = self.bridge['result_value']
        field_name = field_information['key']
        data_column_type = field_information.get('data_type') 
        column_to_select = None
        if(data_column_type in ['numeric', 'decimal', 'integer']):  
            column_to_select = 'numeric_value'
        else:
            column_to_select = 'value'
        
        rv_select = (
            select([_rv.c.well_id, column(column_to_select).label('value')])
            .select_from(_rv)
            .where(
            _rv.c.data_column_id == field_information['data_column_id']))
        rv_select = rv_select.cte(field_information['key'] + '_cte')
        return rv_select

    def create_exclusions_cte(self, screenresult):    
        _rv = self.bridge['result_value']
        _dc = self.bridge['data_column']
        # Create the exclusions subquery: well_id:colnames_excluded
        # Note: this is used to enable filtering and for the UI
        # However, for export types, is overriden by the is_excluded_gen
        # because the generator creates the full col name
        excl_join = _rv.join(_dc, _rv.c.data_column_id==_dc.c.data_column_id)
        excluded_cols_select = (
            select([
                _rv.c.well_id,
                func.array_to_string(
                    func.array_agg(func.lower(_dc.c.name)), 
                    LIST_DELIMITER_SQL_ARRAY).label('exclude')
                ])
            .select_from(excl_join)
            .where(_dc.c.screen_result_id == screenresult.screen_result_id)
            .where(_rv.c.is_exclude)
            .group_by(_rv.c.well_id)
            .order_by(_rv.c.well_id)
            )
        return excluded_cols_select.cte('exclusions')
    
    def get_query(
            self, username, screenresult, param_hash, schema, limit, offset,
            **extra_params):

        DEBUG_SCREENRESULT = False or logger.isEnabledFor(logging.DEBUG)
        logger.info('build screenresult query')
        
        manual_field_includes = set(param_hash.get('includes', []))
        if screenresult.screen.study_type is None:
            manual_field_includes.add('assay_well_control_type')
            
        # general setup
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, schema, **extra_params)
                              
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params )
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
    
        # specific setup 
        
        _aw = self.bridge['assay_well']
        _w = self.bridge['well']
        _sr = self.bridge['screen_result']
        _s = self.bridge['screen']
        _rv = self.bridge['result_value']
        _dc = self.bridge['data_column']
        _reagent = self.bridge['reagent']
        _library = self.bridge['library']
        excluded_cols_select = self.create_exclusions_cte(screenresult)
                    
        # Strategy: 
        # 1. create the base clause, which will build a stored index in 
        # the table: well_query_index;
        # the base query uses only columns: well_id, +filters, +orderings
        base_clause = _aw
        base_clause = base_clause.join(_w, _aw.c.well_id == _w.c.well_id)
        base_clause = base_clause.join(
            _sr, _aw.c.screen_result_id == _sr.c.screen_result_id)
        base_clause = base_clause.join(_s, _sr.c.screen_id == _s.c.screen_id)
        base_clause = base_clause.join(
            _reagent,_w.c.well_id==_reagent.c.well_id, isouter=True)
        base_clause = base_clause.join(
            _library,_w.c.library_id==_library.c.library_id)

        base_custom_columns = {
            'well_id': literal_column('assay_well.well_id'),
            }

        filter_excluded = True
        if ('exclude' not in filter_hash 
                and 'exclude' not in order_params
                and '-exclude' not in order_params):
            filter_excluded = False
        if filter_excluded is True:
            logger.info('filter excluded...')
            base_clause = base_clause.join(
                excluded_cols_select, 
                _aw.c.well_id==excluded_cols_select.c.well_id,isouter=True)
            base_custom_columns['exclude'] = literal_column('exclusions.exclude')
            
        # The base_fields are always included in the query: 
        # - all (screenresult) fields not part of result value lookup
        # - any result_value lookups that are used in sort or filtering
        base_fields = [ fi for fi in field_hash.values() 
            if ( # everything except datacolumns
                fi['scope'] in ['fields.screenresult'] 
                 or fi['key'] in order_params
                 or '-%s' % fi['key'] in order_params
                 or fi['key'] in filter_hash)]
        logger.debug('base fields: %r', base_fields)
        # Using nested selects 
        for fi in [fi for fi in base_fields 
                if fi.get('is_datacolumn', None)]:
            rv_select = self._build_result_value_column(fi)
            base_custom_columns[fi['key']] = rv_select
        # USING CTEs (shown to be slower)           
        # for fi in [fi for fi in base_fields 
        #         if fi.get('is_datacolumn', None)]:
        #     rv_select = self._build_result_value_cte(fi)
        #     base_clause = base_clause.join(
        #         rv_select, _aw.c.well_id==rv_select.c.well_id,
        #         isouter=True)
        #     rv_selector = (
        #         select([rv_select.c.value])
        #         .select_from(rv_select)
        #         .where(rv_select.c.well_id==text('assay_well.well_id'))
        #         )
        #     rv_selector.label(fi['key'])
        #     base_custom_columns[fi['key']] = rv_selector            

        base_query_tables = [
            'assay_well','well','screen_result', 'screen', 'exclusions']
        
        base_columns = self.build_sqlalchemy_columns(
            base_fields, base_query_tables=base_query_tables,
            custom_columns=base_custom_columns) 
        logger.debug('base columns: %r', base_columns)
        # remove is_positive if not needed, to help the query planner
        if ('is_positive' not in filter_hash 
                and 'is_positive' in base_columns
                and 'is_positive' not in order_params
                and '-is_positive' not in order_params):
            del base_columns['is_positive']
            
        # Get the reagent columns and join them - only if in the base_fields
        reagent_resource = self.get_reagent_resource()
        base_reagent_tables = ['reagent','library','well']
        reagent_columns = []
        sub_columns = reagent_resource.build_sqlalchemy_columns(
            base_fields, base_reagent_tables,
            custom_columns=reagent_columns)
        for key,col in sub_columns.items():
            if key not in base_columns:
                if DEBUG_SCREENRESULT: 
                    logger.info('adding reagent column: %r...', key)
                base_columns[key] = col
        if not base_columns:
            raise ProgrammingError(
                'no base columns found in sub_columns: %r', sub_columns)
        base_stmt = select(base_columns.values()).select_from(base_clause)
        base_stmt = base_stmt.where(
            _aw.c.screen_result_id == screenresult.screen_result_id)
        (base_stmt, count_stmt) = \
            self.wrap_statement(base_stmt, order_clauses, filter_expression)
        base_stmt = base_stmt.order_by(asc(column('well_id')))  
        
        # 1.a insert the base statement well ids into the indexing table
        
        m = hashlib.md5()
        compiled_stmt = str(base_stmt.compile(
            dialect=postgresql.dialect(),
            compile_kwargs={"literal_binds": True}))
        
        if DEBUG_SCREENRESULT: logger.info('base_stmt: %r', compiled_stmt)
        
        m.update(compiled_stmt)
        key = m.hexdigest()
        if DEBUG_SCREENRESULT: logger.info('cached query key: %r', key)
        _wellQueryIndex = self.bridge['well_query_index']
        (cachedQuery, create_new_well_index_cache) = \
            CachedQuery.objects.all().get_or_create(key=key)
        if DEBUG_SCREENRESULT: 
            logger.info('create_new_well_index_cache: %r', 
                create_new_well_index_cache)
        if create_new_well_index_cache:
            with get_engine().begin() as conn:
                # NOTE: mixing Django connection with SQA connection
                # - thrown exceptions will rollback the nested SQA transaction
                # see: http://docs.sqlalchemy.org/en/latest/core/connections.html
                try:
                    cachedQuery.sql = compiled_stmt
                    cachedQuery.username = username
                    cachedQuery.uri = \
                        '/screenresult/%s' % screenresult.screen.facility_id
                    cachedQuery.params = json.dumps(
                        { k:v for k,v in param_hash.items() if k!='schema'})
                    cachedQuery.save()
    
                    base_stmt = base_stmt.alias('base_stmt')
                    insert_statement = insert(_wellQueryIndex).\
                        from_select(['well_id', 'query_id'],
                            select([
                                literal_column("well_id"),
                                literal_column(
                                    str(cachedQuery.id)).label('query_id')
                                ]).select_from(base_stmt))
                    if DEBUG_SCREENRESULT:
                        logger.info(
                            'insert stmt: %r',
                            str(insert_statement.compile(
                                dialect=postgresql.dialect(),
                                compile_kwargs={"literal_binds": True})))
                    logger.info('execute screenresult well index insert...')
                    conn.execute(insert_statement)
                    logger.info('screenresult well index complete')
                    count_stmt = _wellQueryIndex
                    count_stmt = select([func.count()]).select_from(count_stmt)
                    count_stmt = count_stmt.where(
                        _wellQueryIndex.c.query_id == cachedQuery.id)
                    logger.info('execute screenresult cached count query...')
                    cachedQuery.count = int(conn.execute(count_stmt).scalar())
                    if cachedQuery.count == 0:
                        cachedQuery.delete()
                        logger.warn('Query generates no results: %r', compiled_stmt)
                    else:
                        cachedQuery.save()
                    # clear out older cached query wells
                    self.clear_cache(by_size=True)
                except Exception, e:
                    cachedQuery.delete()
                    logger.exception('ex on get/create wellquery for screenresult')
                    raise e
        else:
            logger.info('using cached well_query: %r', cachedQuery)
            
        ######
        # Now that query_index table is populated, use it to build the 
        # output query
        logger.info('build screenresult output query...')
        # specific setup 
        base_query_tables = ['assay_well', 'screen']
        # force query to use well_query_index.well_id
        custom_columns = { 
            'well_id': 
                literal_column('well_query_index.well_id'),
            'exclude': literal_column('exclusions.exclude')
        }
        # Use the well_query_index well_ids as the central subquery loop
        _wqx = select([column('id'), column('well_id')]).select_from(_wellQueryIndex)
        _wqx = _wqx.where(_wellQueryIndex.c.query_id == cachedQuery.id)
        _wqx = _wqx.order_by(_wellQueryIndex.c.id)
        if limit > 0:    
            _wqx = _wqx.limit(limit)
        _wqx = _wqx.offset(offset)
        _wqx = _wqx.cte('wqx')
        # Join to well table first to take advantage of well_id foreign key
        # between well and assay_well
        j = join(_wqx, _w, _wqx.c.well_id == _w.c.well_id)
        j = j.join(_aw, _w.c.well_id == _aw.c.well_id )
        j = j.join(_sr, _aw.c.screen_result_id == _sr.c.screen_result_id)
        j = j.join(_s, _sr.c.screen_id == _s.c.screen_id)
        
        # JOIN THE REAGENT COLUMNS
        j = j.join(_reagent,_w.c.well_id==_reagent.c.well_id, isouter=True)
        j = j.join(_library,_w.c.library_id==_library.c.library_id)

        if 'exclude' in field_hash:
            # FIXME: 20171218: create an index for the exclude col on result_value
            # - ultimately needed is a refactor of the result value table:
            # - 1. divide into sm and rnai tables,
            # - 2. refactor (using assay_well?) to index each row in results 
            # (or create an assay_row table)
            j = j.join(excluded_cols_select, 
                excluded_cols_select.c.well_id == _aw.c.well_id, isouter=True)
        # Using nested selects
        for fi in [
            fi for fi in field_hash.values() 
                if fi.get('is_datacolumn', None)]:
            logger.debug('building rv column: %r', fi['key'])
            custom_columns[fi['key']] = self._build_result_value_column(fi)
            
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)

        # GET THE REAGENT COLUMNS AND JOIN THEM
        reagent_resource = self.get_reagent_resource()
        base_reagent_tables = ['reagent','library','well']
        reagent_columns = []
        sub_columns = reagent_resource.build_sqlalchemy_columns(
            field_hash.values(), base_reagent_tables,
            custom_columns=reagent_columns)
        for key,col in sub_columns.items():
            if key not in columns:
                columns[key] = col

        # Force the query to use well_query_index.well_id
        columns['well_id'] = literal_column('wqx.well_id').label('well_id')
        if DEBUG_SCREENRESULT: 
            logger.info('columns: %r', columns.keys())
        ######
        
        stmt = select(columns.values()).select_from(j)
        stmt = stmt.where(
            _aw.c.screen_result_id == screenresult.screen_result_id)
        stmt = stmt.order_by(_wqx.c.id)

        if DEBUG_SCREENRESULT: 
            logger.info(
                'stmt: %s',
                str(stmt.compile(
                    dialect=postgresql.dialect(),
                    compile_kwargs={"literal_binds": True})))
        logger.info('screenresult query built')
        return (field_hash, columns, stmt, count_stmt, cachedQuery, filename)
    
    def build_list_response(self, request, **kwargs):

        logger.info('build_list_response: %r', kwargs)
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        use_raw_lists = request.GET.get(HTTP_PARAM_RAW_LISTS, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
            use_raw_lists = True
        
        manual_field_includes = set(param_hash.get('includes', []))
        manual_field_includes.add('screen_facility_id')
        limit = param_hash.get('limit', 0)        
        try:
            limit = int(limit)
        except ValueError:
            raise BadRequest(
                "Invalid limit '%s' - positive integer required." % limit)
        param_hash['limit'] = limit
        
        offset = param_hash.get('offset', 0)
        try:
            offset = int(offset)
        except ValueError:
            raise BadRequest(
                "Invalid offset '%s' - positive integer required." % offset)
        if offset < 0:    
            offset = -offset
        param_hash['offset'] = offset
        
        is_for_detail = kwargs.pop('is_for_detail', False)

        extra_params = {}
        screen_facility_id = param_hash.pop('screen_facility_id', None)
        extra_params['Screen'] = screen_facility_id
        if not screen_facility_id:
            raise NotImplementedError(
                'must provide a screen_facility_id parameter')
            
        well_id = param_hash.pop('well_id', None)
        if well_id:
            param_hash['well_id__eq'] = well_id

        screenresult = ScreenResult.objects.get(
            screen__facility_id=screen_facility_id)              

        # Note: 20170905 - show_mutual_positives is deprecated: all extra columns
        # are added explicitly
        show_mutual_positives = bool(param_hash.get('show_mutual_positives', False))
        if show_mutual_positives is True:
            extra_params['mutual_positive'] = None
            
        extra_dc_ids = param_hash.get('dc_ids', None)
            
        logger.info('build screen_result schema...')
        schema = self.build_schema(
            screenresult=screenresult,
            show_mutual_positives=show_mutual_positives,
            user=request.user,
            extra_dc_ids=extra_dc_ids)
        logger.info('build screen_result schema - done')

        content_type = self.get_serializer().get_accept_content_type(
            request, format=kwargs.get('format', None))
        
        if content_type == SDF_MIMETYPE:
            manual_field_includes.add('molfile')
        else:
            manual_field_includes.discard('molfile')
        param_hash['includes'] = list(manual_field_includes)

        try:
            logger.info('build screenresult query...')
            (field_hash, columns, stmt, count_stmt, cachedQuery,filename) = \
                self.get_query(
                    request.user.username,
                    screenresult,
                    param_hash,
                    schema,
                    limit,
                    offset, **extra_params)
            # Custom screen result serialization
            
            def is_excluded_gen(rows):
                '''
                Collate all of the excluded columns by well_id.
                - This generator creates the full column names; allowing the 
                "exclusions" query to be simpler.
                '''
                excluded_well_to_datacolumn_map = {}
                for key,field in schema['fields'].items():
                    if ( field.get('is_datacolumn', False) 
                        and field['scope'] 
                            == 'datacolumn.screenresult-%s' % screen_facility_id):
                        
                        for well_id in (
                            ResultValue.objects.filter(
                                    data_column_id=field['data_column_id'])
                                .filter(is_exclude=True)
                                .values_list('well_id', flat=True)):
                            excluded_cols = excluded_well_to_datacolumn_map.get(well_id,[])
                            excluded_cols.append(key)
                            excluded_well_to_datacolumn_map[well_id] = excluded_cols
                for row in rows:
                    well_id = row['well_id']
                    if well_id in excluded_well_to_datacolumn_map:
                        excluded_cols = excluded_well_to_datacolumn_map[well_id]
                        row['exclude'] = sorted(excluded_cols)
                    else:
                        row['exclude'] = None
                    yield row                  
            # Perform custom serialization because each output format will be 
            # different.
            list_brackets = LIST_BRACKETS
            if use_raw_lists:
                list_brackets = None
            image_keys = [key for key,field in field_hash.items()
                if field.get('display_type', None) == 'image']
            ordered_keys = sorted(field_hash.keys(), 
                key=lambda x: field_hash[x].get('ordinal',key))
            list_fields = [ key for (key,field) in field_hash.items() 
                if( field.get('json_field_type',None) == 'fields.ListField' 
                    or field.get('linked_field_type',None) == 'fields.ListField'
                    or field.get('data_type', None) == 'list' ) ]
            value_templates = {key:field['value_template'] 
                for key,field in field_hash.items() 
                    if field.get('value_template', None)}

            title_function = None
            if use_titles is True:
                extra_titles = {
                    'plate': 'Plate Number',
                    'well': 'Well Name',
                    'type': 'Assay Well Control Type',
                    'exclude': 'Excluded'}
                def title_function(key):
                    if key in field_hash:
                        return field_hash[key]['title']
                    elif key in extra_titles:
                        return extra_titles[key]
                    else:
                        return key
            rowproxy_generator = None

            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
            
            logger.info('use_vocab: %r, content_type: %r', use_vocab, content_type)
            if ( use_vocab or content_type != JSON_MIMETYPE):
                # NOTE: xls export uses vocab values
                logger.info('use vocab generator...')
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            else:
                logger.info('do not use vocabularies')
            # FIXME: use closing wrapper
            conn = get_engine().connect()
            if logger.isEnabledFor(logging.DEBUG):
                logger.info(
                    'excute stmt %r...',
                    str(stmt.compile(
                        dialect=postgresql.dialect(),
                        compile_kwargs={"literal_binds": True})))
            logger.info('execute screenresult query...')
            result = conn.execute(stmt)
            logger.info('serialize screenresult response...')
            if rowproxy_generator:
                result = rowproxy_generator(result)

            data = cursor_generator(
                result,ordered_keys,list_fields=list_fields,
                value_templates=value_templates)
            data = closing_iterator_wrapper(data, conn.close)
            logger.info('content_type: %r', content_type)
            if content_type == JSON_MIMETYPE:
                meta = {
                    'limit': limit,
                    'offset': offset,
                    'screen_facility_id': screen_facility_id,
                    'total_count': cachedQuery.count,
                }    
                # Note: is_excluded generator not needed for JSON - already
                # part of the query, and no need to convert to col letters
                data = image_generator(data, image_keys, request) 
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        json_generator(data, meta, is_for_detail=is_for_detail)))
                response['Content-Type'] = content_type
                return response
            elif(content_type == XLS_MIMETYPE or
                content_type == XLSX_MIMETYPE): 
                data = is_excluded_gen(data)
                response = get_xls_response(
                    screen_result_importer.create_output_data(
                        screen_facility_id, field_hash, data),
                    filename, request=request, 
                    title_function=title_function, image_keys=image_keys,
                    list_brackets=list_brackets)
                return response
            elif content_type == SDF_MIMETYPE:
                data = image_generator(
                    is_excluded_gen(data),image_keys, request)
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        sdf_generator(data, title_function=title_function)),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.sdf' % filename
                return response
            elif content_type == CSV_MIMETYPE:
                data = image_generator(
                    is_excluded_gen(data),image_keys, request)
                response = StreamingHttpResponse(
                    ChunkIterWrapper(
                        csv_generator(
                            data,title_function=title_function, 
                            list_brackets=list_brackets)),
                    content_type=content_type)
                response['Content-Disposition'] = \
                    'attachment; filename=%s.csv' % filename
                return response
            else: 
                raise BadRequest('content_type not implemented: %r' % content_type)
        except Exception, e:
            logger.exception('on get list: %r', e)
            raise e  
    
    def get_mutual_positives_columns(self, screen_result_id):
        logger.info('get_mutual_positives_columns...')
        cache_key = '%s_mutual_positive_columns' % screen_result_id
        cached_ids = cache.get(cache_key)
        
        if not cached_ids:
        
            # Note: cache is cleared when any screen_results referenced by 
            # a datacolumn are re-loaded
            # 
            # SS1 Methodology:
            # select distinct(dc.data_column_id) 
            # from assay_well aw0
            # cross join assay_well aw1
            # inner join screen_result sr on aw1.screen_result_id=sr.screen_result_id
            # inner join data_column dc on sr.screen_result_id=dc.screen_result_id
            # where aw0.well_id=aw1.well_id
            # and aw0.is_positive
            # and aw1.is_positive
            # and ( dc.data_type = 'boolean_positive_indicator' 
            #       or dc.data_type = 'partition_positive_indicator' )
            # and aw0.screen_result_id = 941
            # and aw1.screen_result_id <> 941;
                    
            _aw = self.bridge['assay_well']
            _sr = self.bridge['screen_result']
            _dc = self.bridge['data_column']
            
            # Query to find mutual positive data columns:
            
            # OLD METHOD: USES EXISTS
            # - THIS HANGS POSTGRES (8.4 & 9.2) 
            # SELECT DISTINCT wdc.data_column_id 
            # FROM well_data_column_positive_index AS wdc 
            # JOIN data_column ON wdc.data_column_id = data_column.data_column_id 
            # WHERE data_column.screen_result_id != 1090 
            # AND (EXISTS (
            #     SELECT null FROM well_data_column_positive_index AS wdc1 
            #     JOIN data_column AS dc1 ON wdc1.data_column_id = dc1.data_column_id 
            #     WHERE wdc.well_id = wdc1.well_id AND dc1.screen_result_id = 1090));
            # 
            # NEW METHOD: USES TEMPORARY TABLE
            # with wdc1 as (
            # SELECT distinct(well_id) 
            #     FROM well_data_column_positive_index AS wdc1 
            #     JOIN data_column AS dc1 ON wdc1.data_column_id = dc1.data_column_id 
            #     WHERE dc1.screen_result_id = 1090
            # order by well_id
            # ) 
            # SELECT DISTINCT wdc.data_column_id 
            # FROM well_data_column_positive_index AS wdc 
            # JOIN data_column ON wdc.data_column_id = data_column.data_column_id
            # join wdc1 on wdc1.well_id=wdc.well_id 
            # WHERE data_column.screen_result_id != 1090;
    
            
            _aw = self.bridge['assay_well']
            _sr = self.bridge['screen_result']
            _dc = self.bridge['data_column']
            
            _wdc = self.get_create_well_data_column_positive_index().alias('wdc')
            _wdc1 = self.get_create_well_data_column_positive_index().alias('wdc1')
            
            _dc1 = _dc.alias('dc1')
            
            j2 = _wdc1.join(_dc1, _wdc1.c.data_column_id == _dc1.c.data_column_id)
            stmt2 = select([distinct(_wdc1.c.well_id).label('well_id')]).select_from(j2)
            stmt2 = stmt2.where(_dc1.c.screen_result_id == screen_result_id)
            stmt2.order_by('well_id')
            stmt2 = stmt2.cte('wdc1')
    
            j = _wdc.join(_dc, _wdc.c.data_column_id == _dc.c.data_column_id)
            j = j.join(stmt2, stmt2.c.well_id == _wdc.c.well_id)
            stmt = select([distinct(_wdc.c.data_column_id)]).select_from(j)
            stmt = stmt.where(_dc.c.screen_result_id != screen_result_id)
    
            if logger.isEnabledFor(logging.DEBUG):        
                sql = str(stmt.compile(
                    dialect=postgresql.dialect(),
                    compile_kwargs={"literal_binds": True}))
                logger.info('mutual positives statement: %r',sql)
            logger.info('execute mutual positives column query...')
            with get_engine().connect() as conn:
                result = conn.execute(stmt)
                cached_ids =  [x[0] for x in result ]
                logger.info('done, cols %r', cached_ids)
                cache.set(cache_key, cached_ids)
        else:
            logger.info('using cached mutual positive columns')
        return cached_ids

    def get_schema(self, request, **kwargs):

        if not 'screen_facility_id' in kwargs:
            raise Http404(
                'The screenresult schema requires a screen facility ID'
                ' in the URI, as in /screenresult/[facility_id]/schema/')
        facility_id = kwargs.pop('screen_facility_id')
        logger.info('get schema: %r', kwargs)
        if not self._meta.authorization._is_resource_authorized(
            request.user, 'read', screen_facility_id=facility_id):
            logger.info('Permission Denied: ScreenResult.schema, %r, %r',
                request.user, facility_id)
            raise PermissionDenied

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        show_all_other_screens = param_hash.get('show_all_other_screens', False)
        show_mutual_positives = bool(param_hash.get('show_mutual_positives', False))
        extra_dc_ids = param_hash.get('dc_ids', None)
        try:
            screenresult = ScreenResult.objects.get(
                screen__facility_id=facility_id)
            return self.build_response(
                request, self.build_schema(
                    screenresult, 
                    show_mutual_positives=show_mutual_positives,
                    user=request.user, extra_dc_ids=extra_dc_ids), 
                **kwargs)
        except ObjectDoesNotExist, e:
            raise Http404(
                'no screen result found for facility id: %r' % facility_id)

    def build_schema(
            self, screenresult=None, show_mutual_positives=False,
            user=None, extra_dc_ids=None):
        '''
        Note: 20170905 - show_mutual_positives is deprecated: all extra columns
        are added explicitly
        '''

        if user is None:
            raise NotImplementedError(
                'User must be provided for ScreenResult schema view')
        if screenresult is None:
            return super(ScreenResultResource, self).build_schema(user=user)
            
        screen_facility_id = screenresult.screen.facility_id
        
        cache_key = 'screenresult_schema_%s_%s_mutual_pos_%s' \
            % (user.username, screen_facility_id, show_mutual_positives)
        if extra_dc_ids:
            cache_key += '_dcs_.'  + '_'.join(extra_dc_ids)
        logger.info('build screenresult schema: %s', cache_key)
        data = cache.get(cache_key)
        
        def add_well_fields(current_fields):
            well_schema = self.get_reagent_resource().build_schema(
                screenresult.screen.screen_type, user=user)
            newfields = {}
            newfields.update(well_schema['fields'])
            for key,field in newfields.items():
                if 'l' in field['visibility']:
                    field['visibility'].remove('l')
                if 'd' in field['visibility']:
                    field['visibility'].remove('d')
            newfields.update(current_fields)
            
            if screenresult.screen.screen_type == VOCAB_SCREEN_TYPE_SM:
                if 'compound_name' in newfields:
                    newfields['compound_name']['visibility'] = ['l','d']
                    newfields['compound_name']['ordinal'] = 12
            elif screenresult.screen.screen_type == VOCAB_SCREEN_TYPE_RNAI:
                if 'facility_entrezgene_id' in newfields:
                    newfields['facility_entrezgene_id']['visibility'] = ['l','d']
                    newfields['facility_entrezgene_id']['ordinal'] = 10
                if 'facility_entrezgene_symbols' in newfields:
                    newfields['facility_entrezgene_symbols']['visibility'] = ['l','d']
                    newfields['facility_entrezgene_symbols']['ordinal'] = 11
            return newfields
            
            
        if data:
            logger.info('cached: %s', cache_key)
        else:
            logger.info('not cached: %s', cache_key)
            data = super(ScreenResultResource, self).build_schema(user=user)
            
            newfields = add_well_fields(data['fields'])
            if screenresult.screen.study_type is not None:
                del newfields['assay_well_control_type']
                
            max_ordinal = 0
            for fi in newfields.values():
                if fi.get('ordinal', 0) > max_ordinal:
                    max_ordinal = fi['ordinal']
            
            logger.info('map datacolumn definitions into field information definitions...')
            
            datacolumns = self.get_datacolumn_resource()\
                ._get_list_response_internal(
                    user, screen_facility_id=screen_facility_id)
            
            datacolumn_fields = {}
            for dc in datacolumns:
                dc['visibility'] = ['l','d']
                dc['is_datacolumn'] = True
                # NOTE: if user may view screenresults, access level > 1
                # - filtering and ordering are allowed
                dc['filtering'] = True
                dc['ordering'] = True
                datacolumn_fields[dc['key']] = dc
            max_ordinal += len(datacolumn_fields)
            
            newfields.update(datacolumn_fields)
            
            # Note: 20170905 - show_mutual_positives is deprecated: all extra columns
            # are added explicitly
            if show_mutual_positives is True or extra_dc_ids is not None:
                
                visible_screens = self.get_screen_resource()._get_list_response_internal(
                    user=user,
                    includes=[
                        'user_access_level_granted','overlapping_positive_screens',
                        'data_sharing_level','screen_type'])
                visible_screens = { screen['facility_id']:screen 
                    for screen in visible_screens }
                reported_screen = visible_screens[screen_facility_id]
                
                reference_datacolumns = self.get_datacolumn_resource()\
                    ._get_list_response_internal(
                        screen_type=reported_screen['screen_type'])
                reference_datacolumns = { dc['data_column_id']:dc 
                    for dc in reference_datacolumns }

                other_datacolumns = []
                # Note: 20170905 - show_mutual_positives is deprecated: all extra columns
                # are added explicitly
                if show_mutual_positives:
                    temp_datacolumns = self.get_datacolumn_resource()\
                        ._get_list_response_internal(
                            user,
                            screen_facility_id__in=
                                reported_screen['overlapping_positive_screens'])
                    for dc in temp_datacolumns:
                        reference_dc = reference_datacolumns[dc['data_column_id']]
                        if reference_dc['positives_count'] > 0:
                            other_datacolumns.append(dc)
                if extra_dc_ids:
                    extra_datacolumns = self.get_datacolumn_resource()\
                        ._get_list_response_internal(
                            user, 
                            data_column_id__in=extra_dc_ids)                            

                    if self._meta.authorization.is_restricted_view(user):
                        overlapping = reported_screen['overlapping_positive_screens']
                        for dc in extra_datacolumns:
                            if dc['user_access_level_granted'] == 1:
                                if reported_screen['user_access_level_granted'] < 3:
                                    continue
                                elif (dc['screen_facility_id'] 
                                    not in overlapping):
                                    continue
                                else:
                                    logger.info('allowed level 1 col'
                                        '%r: %r: %r', 
                                        dc['key'],dc['title'],dc['data_column_id'])
                                    other_datacolumns.append(dc)    
                            else:
                                other_datacolumns.append(dc)
                    else:
                        other_datacolumns.extend(extra_datacolumns)
                
                other_datacolumns = { dc['data_column_id']:dc 
                    for dc in other_datacolumns }
                decorated = [
                    (dc['screen_facility_id'],dc['ordinal'], dc) 
                        for dc in other_datacolumns.values()]
                decorated.sort(key=itemgetter(0,1))
                other_datacolumns = [dc for fid,ordinal,dc in decorated]
                                            
                datacolumn_fields = {}
                for i, dc in enumerate(other_datacolumns):
                    dc['visibility'] = ['l','d']
                    dc['is_datacolumn'] = True
                    dc['ordinal'] = max_ordinal + i
                    if dc['user_access_level_granted'] > 1:
                        dc['ordering'] = True
                        dc['filtering'] = True
                    datacolumn_fields[dc['key']] = dc
                    
                max_ordinal += len(datacolumn_fields)
                
                newfields.update(datacolumn_fields)
                
                
            data['fields'] = newfields
                
            logger.info('build screenresult schema done')
            cache.set(cache_key, data)
            
        return data

    def create_otherscreen_field(self,screen):
        columnName = "screen_%s" % screen.facility_id
        _dict = {}
        _dict['title'] = '%s (%s)' % (screen.facility_id, screen.title) 
        _dict['description'] =  _dict['title']
        _dict['is_screen_column'] = True
        _dict['key'] = columnName
        _dict['screen_facility_id'] = screen.facility_id
        _dict['visibility'] = ['']
        return (columnName, _dict)
        
    @write_authorization
    def put_list(self, request, **kwargs):
        return self.patch_detail(request, **kwargs)
    
    @write_authorization
    def post_list(self, request, **kwargs):
        return self.patch_detail(request, **kwargs)
    
    @write_authorization
    def patch_list(self, request, **kwargs):
        return self.patch_detail(request, **kwargs)

    @write_authorization
    @un_cache        
    def delete_list(self, request, **kwargs):
        return self.delete_detail(request, **kwargs)
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def delete_detail(self, request, **kwargs):
        if 'screen_facility_id' not in kwargs:
            raise BadRequest('screen_facility_id is required')
        screen_facility_id = kwargs['screen_facility_id']
        # Clear cache: note "all" because mutual positive columns can change
        # self.clear_cache(by_uri='/screenresult/%s' % screen_facility_id)
        self.clear_cache(all=True)
        screen = Screen.objects.get(facility_id=screen_facility_id)

        logger.info('delete screen result for %r', screen_facility_id)
        if hasattr(screen, 'screenresult'):
            screen_result = screen.screenresult
            logger.info('screen result: %r exists, deleting extant data',
                screen_result)
            screen_result.datacolumn_set.all().delete()
            screen_result.assaywell_set.all().delete()
            screen_result.screen.assayplate_set\
                .filter(library_screening__isnull=True).delete()
            screen.screenresult.delete()
            logger.info('screen_result deleted')
            screen_log = self.make_log(request, **kwargs)
            screen_log.ref_resource_name = 'screen'
            screen_log.key = screen.facility_id
            screen_log.uri = '/'.join(['screen', screen_log.key])
            screen_log.save()
        else:
            raise BadRequest('no screen result to delete')

        return HttpResponse(status=204)
            
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def patch_detail(self, request, **kwargs):
        
        if 'screen_facility_id' not in kwargs:
            raise BadRequest('screen_facility_id is required')
        screen_facility_id = kwargs['screen_facility_id']
        screen = Screen.objects.get(facility_id=screen_facility_id)

        data = self.deserialize(request)
        meta = data[API_RESULT_META]
        columns = data['fields']
        result_values = data[API_RESULT_DATA]
        
        if ('screen_facility_id' in meta 
            and meta['screen_facility_id'] != screen_facility_id):
            logger.warn(
                'screen_facility_id in file %r does not match url: %r'
                % (meta['screen_facility_id'], screen_facility_id))
        
        self.create_screen_result(request, screen, columns, result_values)
             
        if not self._meta.always_return_data:
            response_message = {'success': {
                API_MSG_RESULT: 'screen result loaded'}}
            response = self.build_response(request, response_message, **kwargs)
            response['Content-Disposition'] = (
                'attachment; filename=screen_result_loading_success-%s.xlsx' 
                % screen_facility_id )
            return response
        else:
            response = self.get_list(request, **kwargs)             
            response.status_code = 200
            return response 
    
    def create_screen_result(self, request, screen, columns, result_values, **kwargs):
        '''
        Create a screen results for the given screen (internal callers):
        - if results exist, replaces
        '''
        
        schema = self.build_schema(user=request.user)
        
        # Clear cache: note "all" because mutual positive columns can change
        # self.clear_cache(by_uri='/screenresult/%s' % screen_facility_id)
        self.clear_cache(all=True)
        
        id_attribute = schema['id_attribute']

        try:
            adminuser = ScreensaverUser.objects.get(username=request.user.username)
        except ObjectDoesNotExist as e:
            logger.error('admin user: %r does not exist', request.user.username )
            raise
        screen_log = self.make_log(request, **{ 'action': API_ACTION_PATCH })
        screen_log.ref_resource_name = 'screen'
        screen_log.key = screen.facility_id
        screen_log.uri = '/'.join(['screen', screen_log.key])
        screenresult_log = self.make_log(request, **kwargs)
        screenresult_log.ref_resource_name = 'screenresult'
        screenresult_log.key = screen.facility_id
        screenresult_log.uri = '/'.join(['screenresult', screen_log.key])
        screenresult_log.diffs = {}
        
        with transaction.atomic():

            logger.info('Create screen result for %r', screen.facility_id)
            if hasattr(screen, 'screenresult'):
                screen_result = screen.screenresult
                logger.info('screen result: %r exists, deleting extant data',
                    screen_result)
                screen_result.datacolumn_set.all().delete()
                screen_result.assaywell_set.all().delete()
                screen_result.screen.assayplate_set\
                    .filter(library_screening__isnull=True).delete()
                screen_log.diff_keys = ['last_data_loading_date']
                screen_log.diffs = { 
                    'last_data_loading_date': 
                        [screen_result.date_loaded, screen_log.date_time]
                }
                screen_result.date_loaded = screen_log.date_time
                screen_result.created_by = adminuser
            else:
                screen_result = ScreenResult.objects.create(
                    screen=screen,
                    experimental_well_count=0,
                    replicate_count=0,
                    date_loaded=screen_log.date_time,
                    created_by=adminuser)
                screen_log.diffs = {
                    'last_data_loading_date': [None, screen_log.date_time] }
                logger.info('created screen result: %r', screen_result)
            
            screen_log.save()
            screenresult_log.parent_log = screen_log
            logger.info('created log: %r', screen_log)
            
            logger.info(
                'Create screen result data columns for %r', screen.facility_id)
            sheet_col_to_datacolumn = {}
            derived_from_columns_map = {}
            errors = {}
            for i, (sheet_column, column_info) in enumerate(columns.items()):
                if ( column_info.get('screen_facility_id',None)
                    and column_info['screen_facility_id'] != screen.facility_id ):
                    logger.info('skipping column for screen_facility_id: %r', 
                        column_info.get('screen_facility_id') )
                    continue;
                if column_info.get('ordinal',None) is None:
                    column_info['ordinal'] = i
                try:
                    dc = self.get_datacolumn_resource().patch_obj( 
                    # TODO: chgd 20180105 dc = DataColumnResource().patch_obj(
                        request, column_info, screen_result=screen_result, **kwargs)
                    sheet_col_to_datacolumn[sheet_column] = dc
                    derived_from_columns = column_info.get(
                        'derived_from_columns', None)
                    if derived_from_columns:
                        derived_from_columns = [
                            x.strip().upper() 
                            for x in derived_from_columns.split(',')] 
                        derived_from_columns_map[sheet_column] = \
                            derived_from_columns
                except ValidationError, e:
                    errors.update({ sheet_column: e.errors }) 
            logger.info(
                'create derived_from_columns: %r',derived_from_columns_map)
            logger.debug(
                'sheet_col_to_datacolumn: %r', sheet_col_to_datacolumn)
            for sheet_column,derived_cols in derived_from_columns_map.items():
                if not set(derived_cols).issubset(columns.keys()):
                    col_errors = errors.get(sheet_column, {})
                    col_errors['derived_from_columns'] = (
                        '%s are not in available columns: %s'
                        % (convert_list_vals(derived_cols, delimiter=','),
                           convert_list_vals(columns.keys(), delimiter=',')))
                    errors[sheet_column] = col_errors
                else:
                    parent_column = sheet_col_to_datacolumn[sheet_column]
                    for colname in derived_cols:
                        parent_column.derived_from_columns.add(
                            sheet_col_to_datacolumn[colname])
                    logger.info(
                        'parent: %r, derived from columns: %r', parent_column, 
                        [col for col 
                            in parent_column.derived_from_columns.all()])
                    parent_column.save()
                
            if errors:
                logger.warn('errors: %r', errors)
                raise ValidationError( errors={ 'data_columns': errors })

            logger.debug(
                'sheet_col_to_datacolumn: %r', sheet_col_to_datacolumn.keys())
            try:
                self.create_result_values(
                    screen_result, result_values, sheet_col_to_datacolumn,
                    screenresult_log)
            except ValidationError, e:
                logger.exception('Validation error: %r', e)
                raise e
            # END of Transaction
            
        self.create_data_loading_statistics(screen_result)
    
        # Create indexes
        self.get_create_screen_overlap_indexes()
        
        screenresult_log.diffs.update({ 
            'created_by': [None,adminuser.username],  
            'experimental_well_count': [
                None,screen_result.experimental_well_count],
            'replicate_count': [None,screen_result.replicate_count],
            'channel_count': [None,screen_result.channel_count],
        })
        screenresult_log.save()
        
        if screen.study_type is None:
            logger.info('pre-generate the mutual positives index...')
            self.get_mutual_positives_columns(screen_result.screen_result_id)
            logger.info('done - pre-generate the mutual positives index')
        else:
            logger.info('screen is a study, do not pre-generate mutual positives index')
    
    def create_result_value(
            self, colname, value, dc, well, initializer_dict, 
            assay_well_initializer):
        
        if DEBUG_RV_CREATE:
            logger.info('create result value: %r: %r', colname, value)
        positive_types = [
            'confirmed_positive_indicator',
            'partition_positive_indicator',
            'boolean_positive_indicator' ]
        well_id = well.well_id
        
        rv_initializer = {}
        rv_initializer.update(initializer_dict)
        if colname in rv_initializer.get('exclude', []):
            rv_initializer['is_exclude'] = True
        if 'exclude' in rv_initializer:
            del rv_initializer['exclude']

        key = '%s-%s' % (well_id, colname)
        logger.debug(
            'create result value: %r:%r, colname: %r, dc: %r %r',
            key, value, colname, dc.name, dc.data_type)
        rv_initializer['data_column_id'] = dc.data_column_id

        # TODO: 20170731: migrate the data_type of the datacolumn:
        # "numeric" has been replaced by "integer" and "decimal";
        # "text" has been replaced by "string"
        if dc.data_type in ['numeric','decimal','integer']:
            if  dc.decimal_places > 0:
                # parse, to validate
                parse_val(value, key, 'float')
                # record the raw val, to save all digits (precision)
                rv_initializer['numeric_value'] = value
            else:
                rv_initializer['numeric_value'] = \
                    parse_val(value, key, 'integer')
        else:
            # Text value
            rv_initializer['value'] = value
        # TODO: 20170731: migrate the positive types to a separate 
        # integer only column
        if dc.data_type in positive_types:
            if dc.data_type == 'partition_positive_indicator':
                if value not in PARTITION_POSITIVE_MAPPING:
                    raise ValidationError(
                        key=key,
                        msg='%r val: %r must be one of %r'
                            % (dc.data_type, value, 
                                PARTITION_POSITIVE_MAPPING.keys()))
                value = PARTITION_POSITIVE_MAPPING[value]
                rv_initializer['value'] = value
            elif dc.data_type == 'confirmed_positive_indicator':
                if value not in CONFIRMED_POSITIVE_MAPPING:
                    raise ValidationError(
                        key=key,
                        msg='%r val: %r must be one of %r'
                            % (dc.data_type, value, 
                                CONFIRMED_POSITIVE_MAPPING.keys()))
                value = CONFIRMED_POSITIVE_MAPPING[value]
                rv_initializer['value'] = value
            elif  dc.data_type == 'boolean_positive_indicator':
                value = parse_val(value,key,'boolean')
                rv_initializer['value'] = value

            # NOTE: the positive values are recorded in all cases, but 
            # will only count if the well is experimental and not excluded
            if well.library_well_type != 'experimental':
                raise ValidationError(
                    key=key,
                    msg = ('non experimental well, not considered for positives:'
                     'well: %r, %r, col: %r, type: %r, val: %r'
                    % ( well_id, well.library_well_type, colname, 
                        dc.data_type, value)))
            elif rv_initializer['is_exclude'] is True:
                logger.warn(
                    ('excluded col, not considered for positives:'
                     'well: %r, col: %r, type: %r, val: %r'),
                    well_id, colname, dc.data_type, value)
            else:
                if dc.data_type == 'partition_positive_indicator':
                    if value == PARTITION_POSITIVE_MAPPING['W']:
                        dc.weak_positives_count += 1
                        dc.positives_count += 1
                        assay_well_initializer['is_positive'] = True
                    elif value == PARTITION_POSITIVE_MAPPING['M']:
                        dc.medium_positives_count += 1
                        dc.positives_count += 1
                        assay_well_initializer['is_positive'] = True
                    elif value == PARTITION_POSITIVE_MAPPING['S']:
                        dc.positives_count += 1
                        dc.strong_positives_count += 1
                        assay_well_initializer['is_positive'] = True
                elif dc.data_type == 'confirmed_positive_indicator':
                    if value == CONFIRMED_POSITIVE_MAPPING['CP']:
                        dc.positives_count += 1
                        assay_well_initializer['is_positive'] = True
                elif dc.data_type == 'boolean_positive_indicator':
                    if value is True:
                        dc.positives_count += 1
                        assay_well_initializer['is_positive'] = True
         
            if dc.data_type == 'confirmed_positive_indicator':
                if assay_well_initializer.get(
                        'confirmed_positive_value',PSYCOPG_NULL) != PSYCOPG_NULL:
                    raise ValidationError(
                        key=key,
                        msg=('only one "confirmed_positive_indicator" is'
                            'allowed per row'))
                assay_well_initializer['confirmed_positive_value'] = \
                    rv_initializer['value']
        
        if DEBUG_RV_CREATE:
            logger.info('rv_initializer: %r', rv_initializer)
        return rv_initializer
    
    @transaction.atomic
    def create_result_values(
            self, screen_result, result_values, sheet_col_to_datacolumn,
            screenresult_log):
        logger.info(
            'create result values for %r ...', screen_result.screen.facility_id)

        fieldnames = [
            'well_id', 'data_column_id',  # 'result_value_id',
            'value', 'numeric_value', 'is_positive',
            'is_exclude', 'assay_well_control_type',
            
        ]
        assay_well_fieldnames = [
            'screen_result_id', 'well_id', 'plate_number', 'is_positive',
            'confirmed_positive_value','assay_well_control_type',
        ]
        meta_columns = ['well_id', 'assay_well_control_type', 'exclude']
        
        with SpooledTemporaryFile(max_size=MAX_SPOOLFILE_SIZE) as f,\
            SpooledTemporaryFile(max_size=MAX_SPOOLFILE_SIZE) as assay_well_file:
        
            writer = unicodecsv.DictWriter(
                f, fieldnames=fieldnames, delimiter=str(','),
                lineterminator="\n")
            assay_well_writer = unicodecsv.DictWriter(
                assay_well_file, fieldnames=assay_well_fieldnames, 
                delimiter=str(','), lineterminator="\n")
            
            rows_created = 0
            rvs_to_create = 0
            # plates_max_replicate_loaded = {}
            logger.info('write temp file result values for screen: %r ...',
                screen_result.screen.facility_id)
            errors = {}
            while True:
                try: 
                    # iterating will trigger parsing
                    result_row = result_values.next()
                    initializer_dict = { 
                        fieldname:PSYCOPG_NULL for fieldname in fieldnames}
                    assay_well_initializer = { 
                        fieldname:PSYCOPG_NULL 
                            for fieldname in assay_well_fieldnames}
                    for meta_field in meta_columns:
                        if meta_field in result_row:
                            initializer_dict[meta_field] = result_row[meta_field]
                    try:
                        well = Well.objects.get(well_id=result_row['well_id']) 
                    except ObjectDoesNotExist, e:
                        logger.info('well not found: %r', result_row['well_id'])
                        raise ObjectDoesNotExist('well not found: %r' % result_row['well_id'])
                    # FIXME: check for duplicate wells
                    assay_well_initializer.update({
                        'screen_result_id': screen_result.screen_result_id,
                        'well_id': well.well_id,
                        'plate_number': well.plate_number,
                        'is_positive': False })
                    
                    if result_row['assay_well_control_type']:
                        allowed_control_well_types = ['empty', 'dmso']
                        if well.library_well_type in allowed_control_well_types:
                            assay_well_initializer['assay_well_control_type'] = \
                                result_row['assay_well_control_type']
                        else:
                            raise ValidationError(
                                key=well.well_id,
                                msg='control wells must be one of %r, found: %r'
                                 % (allowed_control_well_types, well.library_well_type))

                    for colname, val in result_row.items():
                        if DEBUG_RV_CREATE:
                            logger.info('result value to create: %r, %r', colname, val)
                        if colname in meta_columns:
                            continue
                        if val is None:
                            continue
                        if colname not in sheet_col_to_datacolumn:
                            logger.debug('extra col found in the Data sheet: %r', 
                                colname)
                            continue      
                        dc = sheet_col_to_datacolumn[colname]
                        rv_initializer = self.create_result_value(
                            colname, val, dc, well, initializer_dict, 
                            assay_well_initializer)
                        
                        # 20170424 - remove data loading replicate tracking - per JAS
                        # if (dc.is_derived is False
                        #     and well.library_well_type == 'experimental'
                        #     and rv_initializer['is_exclude'] is not True):
                        #     max_replicate = \
                        #         plates_max_replicate_loaded.get(well.plate_number, 0)
                        #     if dc.replicate_ordinal > max_replicate:
                        #         plates_max_replicate_loaded[well.plate_number] = \
                        #             dc.replicate_ordinal
                        # else:
                        #     logger.debug(('not counted for replicate: well: %r, '
                        #         'type: %r, initializer: %r'), 
                        #         well.well_id, well.library_well_type, rv_initializer)   
                        writer.writerow(rv_initializer)
                        rvs_to_create += 1
                
                    assay_well_writer.writerow(assay_well_initializer)
                    rows_created += 1
                    if rows_created % 10000 == 0:
                        logger.info(
                            'wrote %d result rows to temp file', rows_created)

                except ValidationError, e:
                    logger.info('validation error: %r', e)
                    errors.update(e.errors) 
                except StopIteration, e:
                    break
            
            if errors:
                logger.warn('errors: %r', errors)
                raise ValidationError(errors=errors)
            
            if not rvs_to_create:
                raise ValidationError( errors={ 
                    'result_values': 'no result values were parsed' })

            logger.info('result_values: rows: %d, result_values to create: %d',
                rows_created, rvs_to_create)

            logger.info(
                'use copy_from to create %d assay_wells...', rows_created)
            assay_well_file.seek(0)
            
            with connection.cursor() as conn:
                conn.copy_from(
                    assay_well_file, 'assay_well', sep=str(','), 
                    columns=assay_well_fieldnames, null=PSYCOPG_NULL)
                logger.info('assay_wells created.')
                
                logger.info(
                    'use copy_from to create %d result_values...', rvs_to_create)
                f.seek(0)
                conn.copy_from(
                    f, 'result_value', sep=str(','), columns=fieldnames, 
                    null=PSYCOPG_NULL)
                logger.info('result_values created.')
            screenresult_log.diffs.update({
                'result_values_created': [None,rvs_to_create],
                'assay_wells_loaded': [None, rows_created]
                })
        for dc in sheet_col_to_datacolumn.values():
            dc.save()
            
        # 20170424 - remove replicate tracking for data load - per JAS
        # if plates_max_replicate_loaded:
        #     plates_max_replicate_loaded = sorted(
        #         plates_max_replicate_loaded.values())
        #     logger.info(
        #         'plates_max_replicate_loaded: %r', plates_max_replicate_loaded)
        #     screen_result.screen.min_data_loaded_replicate_count = \
        #         plates_max_replicate_loaded[0]
        #     screen_result.screen.max_data_loaded_replicate_count = \
        #         plates_max_replicate_loaded[-1]
        # else:
        #     screen_result.screen.min_data_loaded_replicate_count = 1
        #     screen_result.screen.max_data_loaded_replicate_count = 1
        # screen_result.screen.save()
        # logger.info('screen_result.screen.max_data_loaded_replicate_count: %r',
        #     screen_result.screen.max_data_loaded_replicate_count)
        
        logger.info('create_result_values - done.')

    # REMOVED - find or create assay plates for screenresult load wells
    # This is removed because the only reason for creating these "assay_plates"
    # is to find the min/max replicates loaded for assay plates with data loaded.
    # NOTE 1: This stat could be calculated during load time instead - see: 
    # "plates_max_replicate_loaded" tracking hash in screen_result load process
    # NOTE 2: Per discussion, this stat will be dropped in SS2 - 201612, per
    # discussion with JenS   
    # def find_or_create_assay_plates_data_loaded(self,screen_result):
    # 
    #     # FIXME: create stats needed without creating assay_plates
    #     
    #     # Create assay plates if the don't exist
    #     # SS1 strategy:
    #     # see ScreenResult.findOrCreateAssayPlatesDataLoaded(plate_number, replicates_loaded)
    #     # a. find the max replicate for each plate in the result_values
    #     # 1. find most recent library screening for that plate
    #     # 2. find the assay plates (all replicates) for only that library screening
    #     # 3. if assay plate count < max replicate screened for that plate, then create ap's
    #     # 3.a in this case create the assay plates sans copy information
    #     # 4. re-run the stat to determine how many assay plates have data loaded:
    #     # - ap.replicate_ordinal in distinct(data_column.replicate_ordinal) for dc join assay_wells  
    #     
    #     sql = (
    #         'with replicates_plate_loaded as ( '
    #         'select plate_number, max(replicate_ordinal) plate_max '
    #         '  from screen  '
    #         '  join screen_result sr using(screen_id) '
    #         '  join data_column dc using(screen_result_id)'
    #         '  join result_value using(data_column_id)'
    #         '  join well using(well_id)'
    #         '  where screen.facility_id = %(facility_id)s '
    #         '  group by well.plate_number order by well.plate_number),'
    #         'replicates_assay_plate_screened as ( '
    #         'select plate_number, count(*) '
    #         '  from assay_plate ap '
    #         '  join library_screening ls on(ap.library_screening_id=ls.activity_id) ' 
    #         '  join screen using(screen_id) '
    #         '  where screen.facility_id = %(facility_id)s '
    #         '  group by plate_number order by plate_number )'
    #         'select '
    #         'rp.plate_number, '
    #         'rp.plate_max, '
    #         'rps.count '
    #         'from replicates_plate_loaded rp '
    #         'left join replicates_assay_plate_screened rps using(plate_number) '
    #         'where rp.plate_max > rps.count or rps.count is null; '
    #         )
    #     conn = self.bridge.get_engine().connect()
    #     result_proxy = conn.execute(
    #         sql, { 'facility_id': screen_result.screen.facility_id })
    #     
    #     assay_plates_created = []
    #     for (plate_number, replicates_needed, replicates_extant ) in result_proxy:
    #         if not replicates_extant:
    #             replicates_extant = 0
    #         replicates_to_create = (replicates_needed-replicates_extant)
    #         logger.info('plate_number: %r, replicates_to_create: %d', 
    #             plate_number, replicates_to_create)
    #         for i in range(replicates_extant, replicates_extant+replicates_needed):
    #             assay_plates_created.append(
    #                 AssayPlate.objects.create(
    #                     screen=screen_result.screen,
    #                     plate_number=plate_number,
    #                     replicate_ordinal=(i)))
    #     
    #     logger.info('TODO: create user message: created assay plates: %d, %r', 
    #         len(assay_plates_created), 
    #         [ap.plate_number for ap in assay_plates_created])
    #     
    #     # step two: set the new 'is_loaded' flag for all assay plates loaded 
    #     # in this screen result
        

    @transaction.atomic
    def create_data_loading_statistics(self, screen_result):
        with get_engine().connect() as conn:
            # NOTE: mixing Django connection with SQA connection
            # - thrown exceptions will rollback the nested SQA transaction
            # see: http://docs.sqlalchemy.org/en/latest/core/connections.html
            screen_facility_id = screen_result.screen.facility_id
            sql_experimental_wells_loaded = (
                'select count(*) '
                'from screen s '
                'join screen_result sr using(screen_id) '
                'join assay_well aw using(screen_result_id) '
                'join well w using(well_id) '
                'where w.library_well_type = %s '
                'and s.facility_id = %s; ')
            screen_result.experimental_well_count = int(
                conn.execute(
                    sql_experimental_wells_loaded,
                    ('experimental', screen_facility_id))
                .scalar() or 0)
            
            sql_replicate_count = (
                 'select max(replicate_ordinal) '
                 'from data_column ' 
                 'join screen_result using(screen_result_id) ' 
                 'join screen using (screen_id) '
                 'where facility_id = %s;')
            screen_result.replicate_count = int(
                conn.execute(
                    sql_replicate_count, screen_facility_id)
                .scalar() or 0)
            if screen_result.replicate_count == 0:
                screen_result.replicate_count = 1
                
            sql_channel_count = (
                 'select max(channel) '
                 'from data_column ' 
                 'join screen_result using(screen_result_id) ' 
                 'join screen using (screen_id) '
                 'where facility_id = %s;')
            screen_result.channel_count = int(
                conn.execute(
                    sql_channel_count, screen_facility_id)
                .scalar() or 0)
            if screen_result.channel_count == 0:
                screen_result.channels_count = 1
                
        screen_result.save()
            

class DataColumnAuthorization(ScreenAuthorization):
    
    def __init__(self, *args, **kwargs):
        ScreenAuthorization.__init__(self, *args, **kwargs)
        self.screen_resource = None
        
    def get_screen_resource(self):
        if self.screen_resource is None:
            self.screen_resource = ScreenResource()
        return self.screen_resource
    
    def filter(self, user, filter_expression):
        # TODO: test: replace with filter_in_sql for performance
        
        if self.is_restricted_view(user) is not True:
            return filter_expression
         
        logger.info('create authorized data columns filter for user: %r', user)
        screensaver_user = ScreensaverUser.objects.get(
            username=user.username)
         
        screen_access_dict = self.get_screen_access_level_table(user.username)
         
        # All fields are visible on data access level 2-3 
        # (user & screen dsl 1, screen dsl 0)
        or_clause = [column('screen_data_sharing_level') == 0]
        level_2_3_screen_ids = [ screen['facility_id']
            for screen in screen_access_dict.values() 
                if screen['user_access_level_granted'] in [2,3]]
        if level_2_3_screen_ids:
            or_clause.append(column('screen_facility_id').in_(level_2_3_screen_ids))
        # Positive columns only for data access level 1
        level_1_screen_ids = [ screen['facility_id']
            for screen in screen_access_dict.values() 
                if screen['user_access_level_granted']==1]
        if level_1_screen_ids:
            or_clause.append(and_(
                column('positives_count')>0,
                column('screen_facility_id').in_(level_1_screen_ids)))
        
        auth_filter = or_(*or_clause)
        # add or clauses conditionally instead
        # auth_filter = or_(
        #     column('screen_data_sharing_level') == 0,
        #     column('screen_facility_id').in_(level_2_3_screen_ids),
        #     and_(column('positives_count')>0,
        #         column('screen_facility_id').in_(level_1_screen_ids))
        #     )
        logger.debug('created data column filter: %r', auth_filter)   
        if filter_expression is not None:
            filter_expression = and_(filter_expression, auth_filter)
        else:
            filter_expression = auth_filter
          
        return filter_expression

    def filter_in_sql(self, user, stmt, screen_table, dc_table):
        return stmt
    
        # if self.is_restricted_view(user) is False:
        #     return stmt
        # 
        # logger.info('create authorized data columns filter for user: %r', user)
        # 
        # screen_access_dict = self.get_screen_access_level_table(user.username)
        # 
        # # All fields are visible on data access level 2-3 
        # # (user & screen dsl 1, screen dsl 0)
        # level_2_3_screen_ids = [ screen['facility_id']
        #     for screen in screen_access_dict.values() 
        #         if screen['user_access_level_granted'] in [2,3]]
        # # Positive columns only for data access level 1
        # level_1_screen_ids = [ screen['facility_id']
        #     for screen in screen_access_dict.values() 
        #         if screen['user_access_level_granted']==1]
        # stmt = stmt.where(or_(
        #     screen_table.c.data_sharing_level == 0,
        #     screen_table.c.facility_id.in_(level_2_3_screen_ids),
        #     and_(dc_table.c.positives_count>0,
        #         screen_table.c.facility_id.in_(level_1_screen_ids))
        #     ))
        # return stmt
    
    def get_row_property_generator(self, user, fields, extant_generator):
        '''
        Filter result properties based on authorization rules
        '''
        return self.get_access_level_property_generator(
           user, fields, extant_generator)
    
    def get_access_level_property_generator(self, user, fields, extant_generator):
        # Effective User-Screen Access Level:
        # 
        # 0 - Field level 0 only
        # 1 - Field level 0,1
        # 2 - Field level 0,1,2, Screen Results, Positives Summary
        # 3 - Field level 0,1,2,3, Screen Results, CPRs, Activities, Visits...
        
        if self.is_restricted_view(user) is not True:
            return extant_generator
        else:
            logger.debug(
                'get_access_level_property_generator for user: %r', user)
            screensaver_user = ScreensaverUser.objects.get(username=user.username)

            fields_by_level = self.get_fields_by_level(fields)
            screen_access_dict = self.get_screen_access_level_table(user.username)

            class Row:
                def __init__(self, row):
                    self.row = row
                    self.facility_id = facility_id = row['screen_facility_id']
                    
                    effective_access_level = None
                    screen_data = screen_access_dict.get(facility_id, None)
                    if screen_data:
                        effective_access_level = \
                            screen_data.get('user_access_level_granted',None)
                    if DEBUG_SCREEN_ACCESS:
                        logger.info('screen: %r, effective_access_level: %r', 
                            facility_id, effective_access_level)
                    self.effective_access_level = effective_access_level
                    
                    if ( self.effective_access_level is None
                        or self.effective_access_level < 1 ):
                        logger.error(
                            'dc shown for not visible screen: '
                            '%r, %r, access level: %r',
                            row['data_column_id'], facility_id, 
                            self.effective_access_level)
                            
                    self.allowed_fields = set()
                    if self.effective_access_level is not None:
                        for level in range(0,self.effective_access_level+1):
                            self.allowed_fields.update(fields_by_level[level])
                            if DEBUG_DC_ACCESS:
                                logger.info(
                                    'allow level: %r: %r, %r', 
                                    level, fields_by_level[level], 
                                    self.allowed_fields)
                    else: 
                        logger.warn('user: %r has no access level for screen: %r',
                            user.username, facility_id)
                    if DEBUG_DC_ACCESS:
                        logger.info(
                            'dc: %r:%r:%r user: %r, effective access level: %r',
                            row['data_column_id'],row['name'], facility_id, 
                            screensaver_user.username, self.effective_access_level )
                    if DEBUG_DC_ACCESS:
                        logger.info('allowed fields: %r', self.allowed_fields)
                def has_key(self, key):
                    if key == 'user_access_level_granted':
                        return True
                    return self.row.has_key(key)
                def keys(self):
                    return self.row.keys();
                def __getitem__(self, key):
                    logger.debug(
                        'key: %r, allowed: %r', key, key in self.allowed_fields)
                    if key == 'user_access_level_granted':
                        return self.effective_access_level
                    if self.row[key] is None:
                        return None
                    else:
                        if key in self.allowed_fields:
                            logger.debug('allow %r: %r', key, self.row[key])
                            return self.row[key]
                        elif DEBUG_DC_ACCESS:
                            logger.info(
                                '%r filter field: %r for restricted user: %r',
                                self.facility_id, key, user.username)
                        return None

            def dc_property_generator(cursor):
                if extant_generator is not None:
                    cursor = extant_generator(cursor)
                for row in cursor:
                    yield Row(row)
            
            return dc_property_generator
        
class DataColumnResource(DbApiResource):

    class Meta:

        queryset = DataColumn.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(
            BasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'datacolumn'
        authorization = DataColumnAuthorization(resource_name)
        ordering = []
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
        super(DataColumnResource, self).__init__(**kwargs)
        self.screen_resource = None
        
    def get_screen_resource(self):
        if self.screen_resource is None:
            self.screen_resource = ScreenResource()
        return self.screen_resource
    
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<data_column_id>\d+)%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/screen/"
                 r"(?P<screen_facility_id>([\w]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url((r"^(?P<resource_name>%s)/for_screen/"
                 r"(?P<for_screen_facility_id>([\w]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_datacolumn_other_screens_view'), 
                name="api_dispatch_datacolumn_other_screens_view"),
        ]    

    def dispatch_datacolumn_other_screens_view(self, request, **kwargs):
        return self.dispatch('list', request, **kwargs)    

    @read_authorization
    def get_detail(self, request, **kwargs):

        data_column_id = kwargs.get('data_column_id', None)
        if not data_column_id:
            logger.info('no data_column_id provided')
            raise NotImplementedError('must provide a data_column_id parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):
        
        logger.info('build datacolumn response...')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = super(DataColumnResource, self).build_schema(user=request.user)
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        
        try:
            # general setup
            screen_facility_id = param_hash.pop('screen_facility_id', None)
            for_screen_facility_id = param_hash.pop('for_screen_facility_id', None)
            manual_field_includes = set(param_hash.get('includes', []))
            # Add fields required to build the system representation
            manual_field_includes.update([
                'screen_id','screen_data_sharing_level',
                'name','data_type','decimal_places','ordinal',
                'screen_facility_id','data_column_id','user_access_level_granted'])
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)

            filter_expression = self._meta.authorization.filter(
                request.user, filter_expression)


            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            max_ordinal = len(field_hash)
            def create_dc_generator(extant_generator=None):
                ''' 
                Transform the DataColumn records into Resource.Fields for the UI
                '''
                
                class DataColumnRow:
                    def __init__(self, row):
                        self._dict = \
                            DataColumnResource._create_datacolumn_from_row(
                                max_ordinal, row)
                    def has_key(self, key):
                        return key in self._dict
                    def keys(self):
                        return self._dict.keys();
                    def __getitem__(self, key):
                        if key in self._dict:
                            return self._dict[key]

                def datacolumn_fields_generator(cursor):
                    if extant_generator:
                        cursor = extant_generator(cursor)
                    for row in cursor:
                        yield DataColumnRow(row)
                return datacolumn_fields_generator
            
            rowproxy_generator = create_dc_generator()
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(
                        field_hash, rowproxy_generator)

            # specific setup 

            _dc = self.bridge['data_column']
            _sr = self.bridge['screen_result']
            _screen = self.bridge['screen']
            _dc_derived_link = self.bridge['data_column_derived_from_columns']
            _dc_derived = _dc.alias('dc_derived')
            
            base_query_tables = [
                'data_column', 'screen']
            
            custom_columns = {
                # default to admin level; screen_property_generator will update
                'key': (
                    _concat('dc_',_screen.c.facility_id, '_', 
                        cast(_dc.c.data_column_id,sqlalchemy.sql.sqltypes.Text))),
                'user_access_level_granted': literal_column('3'),
                'derived_from_columns': (
                    select([func.array_to_string(
                        func.array_agg(literal_column('name')), 
                        LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([_dc_derived.c.name])
                        .select_from(_dc_derived.join(
                            _dc_derived_link,
                            _dc_derived_link.c.to_datacolumn_id
                                ==_dc_derived.c.data_column_id))
                        .where(_dc_derived_link.c.from_datacolumn_id
                            ==literal_column('data_column.data_column_id'))
                        .order_by(_dc_derived.c.ordinal).alias('inner'))
                    )
                }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            # build the query statement
            
            j = _dc
            j = j.join(_sr, _dc.c.screen_result_id == _sr.c.screen_result_id)
            j = j.join(_screen, _sr.c.screen_id == _screen.c.screen_id)
            stmt = select(columns.values()).select_from(j)
            
            # TODO: test if more efficient filtering in sql
            # stmt = self._meta.authorization.filter_in_sql(
            #     request.user, stmt, _screen, _dc)
            
            if screen_facility_id:
                stmt = stmt.where(_screen.c.facility_id == screen_facility_id)
            if for_screen_facility_id:
                for_screen = Screen.objects.get(facility_id=for_screen_facility_id)
                stmt = stmt.where(_screen.c.screen_type == for_screen.screen_type)
                stmt = stmt.where(_screen.c.facility_id != for_screen_facility_id)

            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)

                    
            stmt = stmt.order_by('ordinal')

            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True )
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

    data_type_lookup = {
        'partition_positive_indicator': {
            'vocabulary_scope_ref': 'resultvalue.partitioned_positive',
            'data_type': 'string',
            'edit_type': 'select'
            },
        'boolean_positive_indicator': {
            'data_type': 'boolean' 
            },
        'confirmed_positive_indicator': {
            'vocabulary_scope_ref': 'resultvalue.confirmed_positive_indicator',
            'data_type': 'string',
            'edit_type': 'select' 
            },
    }
            
    @staticmethod
    def _create_datacolumn_from_row(max_ordinal, row_or_dict ):
        ''' Transform a row as defined by the DataColumnResource query into a
        schema field descriptor.
        '''
        _dict = {
            'data_type': 'string',
            'display_options': None,
            'ordinal': 0,
            'visibility': ['l', 'd'],
            'filtering': True,
            'ordering': True,
            'is_datacolumn': True,
        }

        # TODO: 20170731: migrate the data_type of the datacolumn:
        # "numeric" has been replaced by "integer" and "decimal";
        # "text" has been replaced by "string"
        # TODO: 20170731: migrate the screenresult datacolumn format to use 
        # "vocabulary_scope_ref" for the "positive" column types
        key = row_or_dict['key']
        dc_data_type = row_or_dict['data_type']
        # NOTE: (hack) cache the "assay_data_type" here so that the 
        # screen_result_importer.create_output_data can use it as the "output"
        # data_type; the importer does not use vocabulary_scope_ref, and this 
        # is neeed to preserve symmetry of read/write files
        # (see TODO above)
        _dict['assay_data_type'] = dc_data_type
        if dc_data_type  in ['numeric','decimal','integer']:
            if row_or_dict.has_key('decimal_places'):
                decimal_places = row_or_dict['decimal_places']
                if decimal_places > 0:
                    _dict['data_type'] = 'decimal'
                    _dict['display_options'] = (
                         '{ "decimals": %s }' % decimal_places)
                else:
                    _dict['data_type'] = 'integer'
        elif dc_data_type in DataColumnResource.data_type_lookup:
            logger.debug(
                'update data type: %r, %r', 
                dc_data_type, DataColumnResource.data_type_lookup[dc_data_type])
            _dict.update(DataColumnResource.data_type_lookup[dc_data_type])
        elif dc_data_type in ['string','text']:
            _dict['display_type'] = 'full_string'
        logger.debug('_create_datacolumn_from_row: %r, %r, %r',
            key, dc_data_type, _dict )                            
        _dict['ordinal'] = max_ordinal + row_or_dict['ordinal']
        for k in row_or_dict.keys():
            if k not in _dict:
                _dict[k] = row_or_dict[k]
        logger.debug('created datacolumn from row: %r, %r', key, _dict)
        return _dict

    # FIXME: deprecated
    @staticmethod
    def _create_datacolumn_from_orm(dc):
        ''' Transform an ORM DataColumn record into a
        schema field descriptor.
        '''

        screen_facility_id = dc.screen_result.screen.facility_id
        screen = Screen.objects.get(facility_id=screen_facility_id)
        columnName = "dc_%s_%s" % (screen_facility_id, default_converter(dc.name))
        _dict = {}
        _dict.update(model_to_dict(dc))
        _dict['title'] = '%s [%s]' % (dc.name, screen_facility_id) 
        _dict['description'] = _dict['description'] or _dict['title']
        _dict['mouseover'] = '%s: %s - %s' % (screen_facility_id, screen.title, dc.name)
        _dict['comment'] = dc.comments
        _dict['is_datacolumn'] = True
        _dict['key'] = columnName
        _dict['scope'] = 'datacolumn.screenresult-%s' % screen_facility_id
        _dict['screen_facility_id'] = screen_facility_id
        _dict['assay_data_type'] = dc.data_type
        _dict['derived_from_columns'] = [x.name for x in dc.derived_from_columns.all()]
        _dict['visibility'] = ['api']
        _dict['filtering'] = True

        # TODO: 20170731: migrate the data_type of the datacolumn:
        # "numeric" has been replaced by "integer" and "decimal";
        # "text" has been replaced by "string"
        if dc.data_type in ['numeric','decimal','integer']:
            if _dict.get('decimal_places', 0) > 0:
                _dict['data_type'] = 'decimal'
                _dict['display_options'] = \
                    '{ "decimals": %s, "orderSeparator": "" }' % _dict['decimal_places']
            else:
                _dict['data_type'] = 'integer'
        elif dc.data_type in DataColumnResource.data_type_lookup:
            _dict.update(DataColumnResource.data_type_lookup[dc.data_type])
        else:
            _dict['data_type'] = 'string'
        _dict['data_column_id'] = dc.data_column_id
        logger.info('create dc from orm: %r: %r', columnName, _dict)
        return (columnName, _dict)

    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_detail(self, request, **kwargs):
        # TODO: 20170731: allow data_column_id to be passed as an arg so 
        # that the DataColumn may be patched external from the Screen Result
        raise NotImplementedError
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_list(self, request, **kwargs):
        # TODO: 20170731: allow data_column_id to be passed as an arg so 
        # that the DataColumn may be patched external from the Screen Result
        raise NotImplementedError
    
    @write_authorization
    @transaction.atomic    
    def patch_obj(self, request, deserialized, screen_result=None, **kwargs):
        
        # TODO: 20170731: allow data_column_id to be passed as an arg so 
        # that the DataColumn may be patched external from the Screen Result
        if screen_result is None:
            raise BadRequest('screen_result is required')
        logger.debug('patch_obj %s', deserialized)
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        # FIXME: default values
        initializer_dict = {
            'is_derived': False,
            'is_follow_up_data': False,
            'positives_count': 0,
            'strong_positives_count': 0,
            'medium_positives_count': 0,
            'weak_positives_count': 0
        }
        
        for key in fields.keys():
            if deserialized.get(key, None) is not None:
                val = deserialized[key]
                initializer_dict[key] = parse_val(
                    val, key, fields[key]['data_type']) 
                logger.debug('key: %r, val: %r, parsed: %r',
                    key, val, initializer_dict[key])

        errors = self.validate(initializer_dict, schema=schema, patch=False)
        if errors:
            raise ValidationError(errors)

        # TODO: 20170731: migrate the data_type of the datacolumn:
        # "numeric" has been replaced by "integer" and "decimal";
        # "text" has been replaced by "string"
        if initializer_dict['data_type'] == 'numeric':
            if initializer_dict.get('decimal_places',0) > 0:
                initializer_dict['data_type'] = 'decimal'
            else:
                initializer_dict['data_type'] = 'integer'
        if initializer_dict['data_type'] == 'text':
            initializer_dict['data_type'] = 'string'
                
        try:
            
            logger.debug('initializer dict: %s', initializer_dict)
            data_column = DataColumn(
                screen_result=screen_result)
            for key, val in initializer_dict.items():
                if hasattr(data_column, key):
                    setattr(data_column, key, val)
            data_column.save()
            return data_column
        except Exception, e:
            logger.exception('on patch detail')
            raise e  
    

class CopyWellResource(DbApiResource):
    
    COPYWELL_KEY = '{library_short_name}/{copy_name}/{well_id}'
    
    class Meta:
        
        queryset = CopyWell.objects.all().order_by('well_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'copywell'
        authorization = UserGroupAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):

        super(CopyWellResource, self).__init__(**kwargs)
        self.plate_resource = None
        
    def get_plate_resource(self):
        if self.plate_resource is None: 
            self.plate_resource = LibraryCopyPlateResource()
        return self.plate_resource

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('search'), name="api_search"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<copy_name>[\w.\-\+ ]+)" 
                r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<library_short_name>[\w.\-\+: ]+)"
                r"/(?P<copy_name>[\w.\-\+ ]+)" 
                r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<library_short_name>[\w.\-\+: ]+)"
                r"/(?P<copy_name>[\w.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]

    @read_authorization
    def get_detail(self, request, **kwargs):

        library_short_name = kwargs.get('library_short_name', None)
        if not library_short_name:
            logger.info('no library_short_name provided')
         
        copy_name = kwargs.get('copy_name', None)
        if not copy_name:
            raise NotImplementedError('must provide a copy_name parameter')
        
        well_id = kwargs.get('well_id', None)
        if not well_id:
            raise NotImplementedError('must provide a well_id parameter')

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, schema=None, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        is_for_detail = kwargs.pop('is_for_detail', False)

        if schema is None:
            raise Exception('schema not initialized')
        well_id = param_hash.pop('well_id', None)
        if well_id:
            param_hash['well_id__eq'] = well_id
        copy_name = param_hash.pop('copy_name', None)
        if copy_name:
            param_hash['copy_name__eq'] = copy_name
        library_short_name = param_hash.pop('library_short_name', None)
        if library_short_name:
            param_hash['library_short_name__eq'] = library_short_name
        library_screen_type = param_hash.pop('library_screen_type',None);

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
  
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)

            # TODO: remove this restriction if the query can be optimized
            # if filter_expression is None:
            #     raise InformationError(
            #         key='Input filters ',
            #         msg='Please enter a filter expression to begin')
            # else:
            logger.debug('filters: %r', readable_filter_hash)                      
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)

            # specific setup 
            base_query_tables = [
                'copy_well', 'copy', 'plate', 'well', 'library']
            
            _cw = self.bridge['copy_well']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _p = self.bridge['plate']
            _w = self.bridge['well']

            # TODO: optimize query join order copy-plate, then all-copy-plate-well, 
            # then copy-plate-well to copy_well
            # copy_plate = (
            #     select([
            #         _p.c.plate_id, _c.c.copy_id,
            #         _p.c.plate_number, _c.c.name,
            #         _p.c.status, _c.c.usage_type,
            #         _p.c.well_volume,_p.c.remaining_well_volume,
            #         _p.c.mg_ml_concentration, _p.c.molar_concentration,
            #         _p.c.screening_count, _p.c.cplt_screening_count,
            #         _l.c.short_name ])
            #     .select_from(_p.join(_c, _p.c.copy_id==_c.c.copy_id)
            #         .join(_l, _c.c.library_id==_l.c.library_id))
            #         ).cte('copy_plate')
            # all_copy_wells = (
            #     select([_w.c.well_id]))
            
            custom_columns = {
                'copywell_id': (
                    _concat(_l.c.short_name,'/',_c.c.name,'/',_w.c.well_id)),
                'volume': case([
                    (_w.c.library_well_type=='experimental', 
                         func.coalesce(_cw.c.volume, 
                             _p.c.remaining_well_volume, _p.c.well_volume) )],
                    else_=None),
                'initial_volume': case([
                    (_w.c.library_well_type=='experimental', 
                         func.coalesce(_cw.c.initial_volume,_p.c.well_volume) )],
                     else_=None),
                'consumed_volume': case([
                    (_w.c.library_well_type=='experimental', 
                        func.coalesce(_cw.c.initial_volume,_p.c.well_volume)-
                            func.coalesce(_cw.c.volume, _p.c.remaining_well_volume) )],
                    else_=None),
                'mg_ml_concentration': case([
                    (_w.c.library_well_type=='experimental', 
                        func.coalesce(
                            _cw.c.mg_ml_concentration,_p.c.mg_ml_concentration) )],
                    else_=None),
                'molar_concentration': case([
                    (_w.c.library_well_type=='experimental', 
                        func.coalesce(
                            _cw.c.molar_concentration,_p.c.molar_concentration) )],
                    else_=None),
                'cumulative_freeze_thaw_count': (
                    (func.coalesce(_p.c.screening_count,0) 
                        + func.coalesce(_p.c.cplt_screening_count,0))),
                
                
                # 'adjustments': case([
                #     (_w.c.library_well_type=='experimental', 
                #         func.coalesce(_cw.c.adjustments, 0) )],
                #     else_=None),
                # Note: the query plan makes this faster than the hash join of 
                # copy-copy_well
                # 'copy_name': literal_column(
                #    '( select copy.name from copy where copy.copy_id=copy_well.copy_id )'
                #    ).label('copy_name')
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            # build the query statement

            j = join(_w, _l, _w.c.library_id == _l.c.library_id)
            j = j.join(_c, _l.c.library_id == _c.c.library_id)
            j = j.join(_p, and_(
                _p.c.copy_id == _c.c.copy_id,
                _w.c.plate_number == _p.c.plate_number))
            j = j.join(_cw, and_(
                _cw.c.well_id == _w.c.well_id,
                _cw.c.copy_id == _c.c.copy_id), isouter=True )
            
            stmt = select(columns.values()).select_from(j)
            
            if library_screen_type is not None:
                stmt = stmt.where(_l.c.screen_type==library_screen_type)
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            if not order_clauses:
                stmt = stmt.order_by('copy_name', 'plate_number', 'well_id')
            
            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None))
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_list must be implemented')
                
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        raise NotImplementedError('delete_obj is not implemented')
    
    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        # TODO: optimize for list inputs (see well.patch)
        logger.debug('patch_obj %s', deserialized)

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        initializer_dict = {}
        id_kwargs = self.get_id(deserialized, schema=schema, validate=True, **kwargs)
        well_id = id_kwargs['well_id']

        try:
            well = Well.objects.get(well_id=well_id)
            library = well.library
        except ObjectDoesNotExist:
            msg = 'well not found: %r' % well_id
            logger.info(msg);
            raise Http404(msg)

        if well.library_well_type != 'experimental':
            logger.info('CopyWell patch: ignore non experimental well: %r', well_id)
            return None
        
        copy_name = id_kwargs['copy_name']
        try:
            librarycopy = Copy.objects.get(
                name=copy_name, library=library)
        except ObjectDoesNotExist:
            msg = 'copy_name not found: %r' % copy_name
            logger.info(msg);
            raise Http404(msg)
        
        try:
            plate = Plate.objects.get(
                plate_number=well.plate_number, copy=librarycopy)
        except ObjectDoesNotExist:
            msg = 'plate not found: %r:%r' % (library.short_name,copy_name)
            logger.info(msg);
            raise Http404(msg)

        # TODO: wrapper for parsing
        logger.debug('fields: %r, deserialized: %r', fields.keys(), deserialized)
        for key in fields.keys():
            if deserialized.get(key, None) is not None:
                initializer_dict[key] = parse_val(
                    deserialized.get(key, None), key, fields[key]['data_type']) 
        
        volume = initializer_dict.get('volume', None)        
        mg_ml_concentration = initializer_dict.get('mg_ml_concentration', None)
        molar_concentration = initializer_dict.get('molar_concentration', None)
#         adjustments = initializer_dict.get('adjustments', None)
        
        if volume is None and molar_concentration is None and mg_ml_concentration is None:
            msg = (
                'Must submit one of [volume, mg_ml_concentration, '
                'molar_concentration]: well: %r' % well_id)
            raise ValidationError({
                'volume': msg,
                'mg_ml_concentration': msg,
                'molar_concentration': msg,
                })
        try:
            copywell = CopyWell.objects.get(
                well=well, copy=librarycopy)
            
            # FIXME: only "set" adjustments if sent from the user
            # if volume is not None:
            #    copywell.adjustments += 1
        except ObjectDoesNotExist:
            # If creating, check that something is updated from the plate values
            if plate.remaining_well_volume == volume:
                volume = None
            if plate.mg_ml_concentration == mg_ml_concentration:
                mg_ml_concentration = None
            if plate.molar_concentration == molar_concentration:
                molar_concentration = None
            
            if all(v is None for v in 
                [volume,mg_ml_concentration,molar_concentration]):
                logger.info('Nothing to edit for: %r', deserialized)
                return None
            
            copywell = CopyWell.objects.create(
                well=well, copy=librarycopy, plate=plate,
                initial_volume = plate.remaining_well_volume)
            # FIXME: only "set" adjustments if sent from the user
            # if volume is not None:
            #    copywell.adjustments += 1
            copywell.save()
            logger.info('created cw: %r', copywell)

        if volume is not None:
            copywell.volume = volume
        # if adjustments is not None:
        #     copywell.adjustments = adjustments
        if mg_ml_concentration is not None:
            copywell.mg_ml_concentration = mg_ml_concentration
        if molar_concentration is not None:
            copywell.molar_concentration = molar_concentration
        
        copywell.save()
        logger.info('patch_obj done')
        return { API_RESULT_OBJ: copywell }
    
    @un_cache
    @transaction.atomic
    def _deallocate_well_volumes(
        self, volume, copywells_to_deallocate, parent_log):
        
        copywells_deallocated = set()
        meta = {}
        for copywell in copywells_to_deallocate:
            copy = copywell.copy
            new_volume = copywell.volume + volume
            
            log = self.make_child_log(parent_log)
            log.key = '/'.join([copy.library.short_name, copy.name, copywell.well_id])
            log.uri = '/'.join([log.ref_resource_name,log.key])
            log.diffs = {
                'volume': [copywell.volume, new_volume]}
            log.parent_log = parent_log
            logger.info('copywell adjusted: %r, %r', log, log.diffs)
            log.save()
            # adjust volume
            copywell.volume = new_volume
            copywell.save()
            copywells_deallocated.add(copywell)
            
        meta[API_MSG_COPYWELLS_DEALLOCATED] = len(copywells_deallocated)
        return meta
    
    @un_cache
    @transaction.atomic
    def _allocate_well_volumes(
        self, volume, copywells_to_allocate, parent_log):
        
        meta = {}
        copywell_volume_warnings = []
        copywells_allocated = set()
        for copywell in copywells_to_allocate:
            copy = copywell.copy
            key = '/'.join([copy.library.short_name, copy.name, copywell.well_id])
            logger.info('copywell: %r: vol: %r, requested: %r', 
                copywell, copywell.volume, volume)
            if copywell.volume < volume:
                copywell_volume_warnings.append(
                    'CopyWell: %s, '
                    '(available: %s uL)' 
                        % lims_utils.convert_decimal(
                            copywell.volume, 1e-6, 1),
                    '(requested: %s uL)' 
                        % lims_utils.convert_decimal(
                            volume, 1e-6, 1))

            new_volume = copywell.volume - volume
            
            log = self.make_child_log(parent_log)
            log.key = key
            log.uri = '/'.join([log.ref_resource_name,log.key])
            log.diffs = {
                'volume': [copywell.volume, new_volume]}
            log.parent_log = parent_log
            logger.info('copywell adjusted: %r, %r', log, log.diffs)
            log.save()
            # adjust volume
            copywell.volume = new_volume
            copywell.save()
            copywells_allocated.add(copywell)
        if copywell_volume_warnings:
            logger.info('%r:%r', 
                API_MSG_LCPS_INSUFFICIENT_VOLUME, copywell_volume_warnings)
            meta[API_MSG_LCPS_INSUFFICIENT_VOLUME] = copywell_volume_warnings
        meta[API_MSG_COPYWELLS_ALLOCATED] = len(copywells_allocated)
        return meta
            
            
    @un_cache
    @transaction.atomic
    def deallocate_cherry_pick_volumes(
        self, cpr, lab_cherry_picks_to_deallocate, parent_log,
        set_deselected_to_zero=False,
        update_screening_count=True):
        '''
        @param update_screening_count (default True) if false this deallocation
        does not affect the copywell.cherry_pick_screening_count;
        NOTE: if we want to track this, then should create a new 
        "update reservation" method; which will only adjust counts if all wells 
        for a plate are deselected.
        '''
        logger.debug('deallocate_cherry_pick_volumes: %r, %r, %r, %r',
            cpr, lab_cherry_picks_to_deallocate, set_deselected_to_zero, 
            update_screening_count)
        copywells_deallocated = []
        plates_adjusted = set()
        # find copy-wells
        for lcp in lab_cherry_picks_to_deallocate:
            logger.info('copywell to deallocate: %r', lcp)
            copy = lcp.copy
            if copy is None:
                raise ProgrammingError(
                    'LabCherryPick is already deallocated: %r', lcp)
            plate = Plate.objects.get(
                plate_number=lcp.source_well.plate_number, 
                copy=lcp.copy)
            plates_adjusted.add(plate)
            copywell_id = self.COPYWELL_KEY.format(
                well_id=lcp.source_well.well_id,
                library_short_name=copy.library.short_name,
                copy_name=copy.name)
            logger.info('copywell to deallocate: %r', copywell_id)
            try:
                copywell = CopyWell.objects.get(
                    well=lcp.source_well, copy=lcp.copy)
                
            except ObjectDoesNotExist:
                logger.error(
                    'copywell to deallocate not located: %r', copywell_id)
                
            log = self.make_child_log(parent_log)
            log.key = '/'.join([
                copy.library.short_name, copy.name, copywell.well_id])
            log.uri = '/'.join([log.ref_resource_name,log.key])
            if set_deselected_to_zero is False:
                new_volume = copywell.volume + cpr.transfer_volume_per_well_approved
            else:
                logger.info('set deallocated copywell to zero: %r', copywell_id)
                new_volume = 0
            log.diffs = {
                'volume': [copywell.volume, new_volume]}
            if update_screening_count is True:
                if copywell.cherry_pick_screening_count > 0:
                    current_cp_screening_count = copywell.cherry_pick_screening_count
                    new_cp_screening_count =  current_cp_screening_count - 1
                    log.diffs['cherry_pick_screening_count'] = [
                        current_cp_screening_count, 
                        new_cp_screening_count]
                    copywell.cherry_pick_screening_count = new_cp_screening_count
            log.parent_log = parent_log
            logger.info('copywell adjusted: %r, %r', log, log.diffs)
            log.save()
            # adjust volume
            copywell.volume = new_volume
            copywell.save()
            copywells_deallocated.append(copywell)
        logger.debug('copywells_deallocated: %r', copywells_deallocated)    
        logger.info('copywells_deallocated: %d', len(copywells_deallocated))    

        if update_screening_count is True:
            logger.debug('plates_adjusted: %r', plates_adjusted)
            logger.info('plates_adjusted: %d', len(plates_adjusted))
            for plate in plates_adjusted:
                plate_log = self.get_plate_resource().make_child_log(parent_log)
                plate_log.key = '/'.join([
                    copy.library.short_name, plate.copy.name, 
                    str(plate.plate_number)])
                plate_log.uri = '/'.join([log.ref_resource_name,log.key])
                if plate.cplt_screening_count < 1:
                    logger.warn(
                        'deallocation: plate: %r, cplt_screening_count already 0',
                        plate_log.key)
                    new_cp_screening_count = 0
                else:
                    new_plate_cplt_screening_count = plate.cplt_screening_count-1
                plate_log.diffs = {
                    'cplt_screening_count': [
                        plate.cplt_screening_count,
                        new_plate_cplt_screening_count ],
                    }
                plate_log.save()
                plate.cplt_screening_count = new_plate_cplt_screening_count
                plate.save()
            
        return { 
            'CPR #': cpr.cherry_pick_request_id,
            API_MSG_COPYWELLS_DEALLOCATED: len(copywells_deallocated),
            API_MSG_LCP_SOURCE_PLATES_DEALLOCATED: len(plates_adjusted)
         }
        
        
    def reserve_cherry_pick_volumes(
        self, cpr, fulfillable_lcps, parent_log, plates_to_ignore=None):
        '''
        @param plates_to_ignore (set) plates to ignore when adjusting the 
        cplt_screening_count
        '''
        copywells_adjusted = []
        plates_adjusted = set()
        # find copy-wells
        for lcp in fulfillable_lcps:
            copy = lcp.copy
            plate = Plate.objects.get(
                plate_number=lcp.source_well.plate_number, 
                copy=lcp.copy)
            plates_adjusted.add(plate)
            try:
                copywell = CopyWell.objects.get(
                    well=lcp.source_well, copy=lcp.copy)
                if copywell.initial_volume is None:
                    # Copywell may have no volumes set, if the copywell was
                    # created to track well-specific concentrations
                    # (the API reports the plate.remaining_well_volume for these)
                    copywell.initial_volume = plate.remaining_well_volume
                    if copywell.volume is not None:
                        logger.warn(
                            'copywell has initial_volume, but not volume: %r/%r, %r', 
                            copy.name, lcp.source_well, copywell.initial_volume)
                    copywell.volume = plate.remaining_well_volume
            except ObjectDoesNotExist:
                # create copy-wells that dne
                try:
                    logger.info('copywell dne: %r/%r, creating', 
                        copy.name, lcp.source_well)
                    logger.info('plate: %r, remaining_well_volume: %r', 
                        plate, plate.remaining_well_volume)
                    copywell = CopyWell.objects.create(
                        well=lcp.source_well, 
                        copy=lcp.copy, 
                        plate=plate,
                        volume=plate.remaining_well_volume,
                        initial_volume=plate.remaining_well_volume)
                    
                except ObjectDoesNotExist:
                    msg = ('plate not found: %r:%r' 
                        % (lcp.source_well.plate_number, lcp.copy.name))
                    logger.warn(msg);
                    raise ValidationError(
                        key='library_plate',
                        msg=msg)
            log = self.make_child_log(parent_log)
            log.key = '/'.join([copy.library.short_name, copy.name, copywell.well_id])
            log.uri = '/'.join([log.ref_resource_name,log.key])
            logger.debug('copywell: %r, volume: %r, cpr.transfer volume: %r', 
                copywell, copywell.volume, cpr.transfer_volume_per_well_approved)
            new_volume = copywell.volume - cpr.transfer_volume_per_well_approved
            current_cp_screening_count = copywell.cherry_pick_screening_count or 0
            new_cp_screening_count =  current_cp_screening_count + 1
            log.diffs = {
                'volume': [copywell.volume, new_volume],
                'cherry_pick_screening_count': [
                    current_cp_screening_count, 
                    new_cp_screening_count]
                }
            log.parent_log = parent_log
            log.save()
            # adjust volume
            copywell.volume = new_volume
            copywell.cherry_pick_screening_count = new_cp_screening_count
            copywell.save()
            copywells_adjusted.append(copywell)
            
        logger.info('plates_adjusted: %r', 
            [(plate.copy.name,plate.plate_number) for plate in plates_adjusted])
        for plate in plates_adjusted:
            
            if plates_to_ignore and plate in plates_to_ignore:
                logger.info('ignoring plate: %r, will not adjust cplt_screening_count',
                    plate)
                continue
            
            plate_log = self.get_plate_resource().make_child_log(parent_log)
            plate_log.key = '/'.join([
                copy.library.short_name, plate.copy.name, str(plate.plate_number)])
            plate_log.uri = '/'.join([log.ref_resource_name,log.key])
            new_plate_cplt_screening_count = plate.cplt_screening_count +1
            plate_log.diffs = {
                'cplt_screening_count': [
                    plate.cplt_screening_count,
                    new_plate_cplt_screening_count ],
                }
            plate_log.save()
            plate.cplt_screening_count = new_plate_cplt_screening_count
            plate.save()
            
        return { 
            'CPR #': cpr.cherry_pick_request_id,
            API_MSG_COPYWELLS_ALLOCATED: len(copywells_adjusted),
            API_MSG_LCP_SOURCE_PLATES_ALLOCATED: len(plates_adjusted)
         }

class CherryPickRequestAuthorization(ScreenAuthorization):        

    def _is_resource_authorized(
        self, user, permission_type, **kwargs):
        authorized = \
            super(CherryPickRequestAuthorization, self)\
                ._is_resource_authorized(user, permission_type, **kwargs)
        if authorized is True:
            return True
        
        return user.is_active

    def filter(self, user, filter_expression):
        
        if self.is_restricted_view(user) is False:
            return filter_expression
        
        screensaver_user = ScreensaverUser.objects.get(username=user.username)
        my_screens = self.get_user_screens(screensaver_user)
        # Can only see own screens
        
        auth_filter = column('screen_facility_id').in_(
            [screen.facility_id for screen in my_screens])
        if filter_expression is not None:
            filter_expression = and_(filter_expression, auth_filter)
        else:
            filter_expression = auth_filter

        return filter_expression

    def get_row_property_generator(self, user, fields, extant_generator):
        # If the user may see the CPR, there are no property restrictions
        return extant_generator
    
    def has_cherry_pick_read_authorization(self, user, cherry_pick_request_id):
        if self.is_restricted_view(user) is False:
            return True
        
        screensaver_user = ScreensaverUser.objects.get(username=user.username)
        my_screens = self.get_user_screens(screensaver_user)
        
        cpr = CherryPickRequest.objects.get(cherry_pick_request_id=cherry_pick_request_id)
        
        return cpr.screen in my_screens
        
class CherryPickRequestResource(DbApiResource):        
    
    class Meta:
    
        queryset = CherryPickRequest.objects.all().order_by('well_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'cherrypickrequest'
        authorization = CherryPickRequestAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):

        super(CherryPickRequestResource, self).__init__(**kwargs)
        
        self.well_resource = None
        self.screen_resource = None
        self.copywell_resource = None
        self.labcherrypick_resource = None
        self.screenercherrypick_resource = None
        self.cpp_resource = None
        self.screensaver_user_resource = None
        self.librarycopyplate_resource = None
        
    def get_librarycopyplate_resource(self):
        if self.librarycopyplate_resource is None:
            self.librarycopyplate_resource = LibraryCopyPlateResource()
        return self.librarycopyplate_resource
    
    def get_screen_resource(self):
        if self.screen_resource is None:
            self.screen_resource = ScreenResource()
        return self.screen_resource
    
    def get_user_resource(self):
        if self.screensaver_user_resource is None:
            self.screensaver_user_resource = ScreensaverUserResource()
        return self.screensaver_user_resource
    
    def get_well_resource(self):
        if self.well_resource is None:
            self.well_resource = WellResource()
        return self.well_resource
    
    def get_copywell_resource(self):
        if self.copywell_resource is None:
            self.copywell_resource = CopyWellResource()
        return self.copywell_resource
    
    def get_labcherrypick_resource(self):
        if self.labcherrypick_resource is None:
            self.labcherrypick_resource = LabCherryPickResource()
        return self.labcherrypick_resource
    
    def get_screenercherrypick_resource(self):
        if self.screenercherrypick_resource is None:
            self.screenercherrypick_resource = ScreenerCherryPickResource()
        return self.screenercherrypick_resource
    
    def get_cherrypickplate_resource(self):
        if self.cpp_resource is None:
            self.cpp_resource = CherryPickPlateResource()
        return self.cpp_resource
    
    def clear_cache(self):
        logger.info('clear_cache: CherryPickRequestResource...')
        DbApiResource.clear_cache(self)
        # NOTE: don't clear dependent resources to avoid circular refererence 
        # recursion; DbApiResource.clear_cache will clear all caches, for now
        # self.get_labcherrypick_resource().clear_cache()
        # self.get_screenercherrypick_resource().clear_cache()
        self.get_screen_resource().clear_cache()
        
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            
            url((r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)" 
                 r"/screener_cherry_pick%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screener_cherry_pick_view'),
                name="api_dispatch_screener_cherry_pick_view"),
            url((r"^(?P<resource_name>%s)/(?P<cherry_pick_request_id>[\d]+)"
                 r"/screener_cherry_pick/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_scp_schema'),
                name="api_get_scp_schema"),
            
            url((r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)" 
                 r"/lab_cherry_pick%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_lab_cherry_pick_view'),
                name="api_dispatch_lab_cherry_pick_view"),
            url((r"^(?P<resource_name>%s)/(?P<cherry_pick_request_id>[\d]+)"
                 r"/lab_cherry_pick/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_lcp_schema'),
                name="api_get_lcp_schema"),
            url((r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)" 
                 r"/source_plate%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_source_plate_view'),
                name="api_dispatch_source_plate_view"),

            url((r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)" 
                 r"/cherry_pick_plate%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_cherry_pick_plate_view'),
                name="api_dispatch_cherry_pick_plate_view"),
            
            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)/set_lab_cherry_picks%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_set_lab_cherry_picks'), 
                name="api_dispatch_set_lab_cherry_picks"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)/set_duplex_lab_cherry_picks%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_set_duplex_lab_cherry_picks'), 
                name="api_dispatch_set_duplex_lab_cherry_picks"),

#             url(r"^(?P<resource_name>%s)"
#                 r"/(?P<cherry_pick_request_id>[\d]+)/plate_lab_cherrypicks%s$" 
#                     % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('dispatch_plate_lab_cherrypicks'), 
#                 name="api_dispatch_plate_lab_cherrypicks"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)/reserve_map_lab_cherry_picks%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_reserve_map_lab_cherry_picks'), 
                name="api_dispatch_reserve_map_lab_cherry_picks"),
                
            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)/cancel_reservation%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_cancel_reservation'), 
                name="api_dispatch_cancel_reservation"),
                
            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)/delete_lab_cherry_picks%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_delete_lab_cherry_picks'), 
                name="api_delete_lab_cherry_picks"),
                
            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/plate_mapping_file%s$"
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view(
                    'get_plate_mapping_file'), 
                    name="api_get_plate_mapping_file"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/lab_cherry_pick_plating%s$"
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view(
                    'dispatch_lab_cherry_pick_plating'), 
                    name="api_dispatch_lab_cherry_pick_plating"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/lab_cherry_pick_plating/schema%s$"
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view(
                    'get_lab_cherry_pick_plating_schema'), 
                    name="get_get_lab_cherry_pick_plating_schema"),
        ]
    
    def dispatch_source_plate_view(self, request, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        manual_field_includes = set(param_hash.get('includes', []))
        if ('-copy_usage_type' not in manual_field_includes):
            manual_field_includes.add('copy_usage_type')
        kwargs['includes']=manual_field_includes
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        kwargs['schema'] = \
            self.get_librarycopyplate_resource().build_schema(user=request.user)
        # NOTE: authorization is performed in LibraryCopyPlateResource
        return self.get_librarycopyplate_resource()\
            .build_list_response(request, **kwargs)
        
    @read_authorization
    def get_plate_mapping_file(self, request, **kwargs):
        
        request_method = request.method.lower()
        if request_method != 'get':
            raise BadRequest('Only GET is allowed')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        schema = super(CherryPickRequestResource, self)\
            .build_schema(user=request.user)

        cpr_id = param_hash['cherry_pick_request_id']
        # add custom authorization for screening users
        if self._meta.authorization.has_cherry_pick_read_authorization(
                request.user, cpr_id) is False:
            raise PermissionDenied
         
        logger.info('get the cpr: %r', cpr_id) 
        cpr_obj = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        cpr = self._get_detail_response_internal(**kwargs)
        
        logger.info('get the cpps: %r', cpr_id) 
        cpps = self.get_cherrypickplate_resource()._get_list_response_internal(
            **{
                'cherry_pick_request_id': cpr_id
            })
        cpps = {cpp['plate_ordinal']:cpp for cpp in cpps }
        
        lcp_schema = self.get_labcherrypick_resource().build_schema(
            user=request.user, library_classification=cpr_obj.screen.screen_type)
        columns_to_include = [
            'library_plate','source_copy_name','source_well_name',
            'source_plate_type','destination_well','destination_plate_type']

        logger.info('get the lcps: %r', cpr_id) 
        lcps = self.get_labcherrypick_resource()._get_list_response_internal(
            **{
                'cherry_pick_request_id': cpr_id,
                'destination_well__is_null': False,
                'includes': [
                    'source_plate_type','destination_plate_type',
                    'location','-structure_image','-molfile',
                    '-library_plate_comment_array'],
                'order_by': ['destination_well']
            })
        plate_name_to_types = defaultdict(set)
        source_plate_to_dest_plates = defaultdict(set)
        plate_name_template = (
            '{requested_by_name}'
            ' ({screen_facility_id})'
            ' CP{cherry_pick_request_id}'
            '_plate_{{cherry_pick_plate_number}}'
            '_of_{number_plates}'
            ).format(**cpr)
        for lcp in lcps:
            plate_name = plate_name_template.format(**lcp)
            plate_name_to_types[plate_name].add(lcp['source_plate_type'])
            plate_copy = '%s/%s' % (
                str(lcp['library_plate']),lcp['source_copy_name'])
            source_plate_to_dest_plates[plate_copy].add(plate_name)
        plate_map = defaultdict(list)
        location_map = {}
        for lcp in lcps:
            plate_name = plate_name_template.format(**lcp)
            if len(plate_name_to_types[plate_name])>1:
                plate_name += ' ' + lcp['source_plate_type']
            plate_map[plate_name].append(lcp)    
            plate_copy = '%s/%s' % (
                str(lcp['library_plate']),lcp['source_copy_name'])
            location_map[plate_copy] = [
                lcp['library_plate'],lcp['source_copy_name'],lcp['location']]
        # Sort internally
        for plate_name in plate_map.keys():
            lcps = plate_map[plate_name]
            lcps = sorted(lcps, key=lambda lcp: lcp['source_well_id'])
            plate_map[plate_name] = lcps
        
        locations = [location_map[x] for x in sorted(location_map.keys())]
        
        vocabularies = DbApiResource.get_vocabularies(lcp_schema['fields'])
        def vocab_function(key,val):
            if key in vocabularies:
                return vocabularies[key][val]['title']
            else:
                return val
        def title_function(key):
            return lcp_schema['fields'][key]['title']

        extra_values = OrderedDict((
            ('Person Visiting', cpr['requested_by_name']),
            ('Screen Number', cpr['screen_facility_id']),
            (u'Volume (uL)', 
                Decimal(cpr['transfer_volume_per_well_approved']) * Decimal(1e6)),
            # NOTE: UTF-8 encodings are supported, but avoid using them because
            # Excel has poor support for reading properly opening encoded csv files
            # (u'Volume (\u00B5L)', 
            #     Decimal(cpr['transfer_volume_per_well_approved']) * Decimal(1e6)),
        ))
        
        # Outputs:
        # Zip file:
        # - assay plate files as: Jen Smith (729) CP44451  Plate 01 of 1 (Run1)
        # - plate locations file
        # - README.txt
        
        readme_types_warning = '\n'.join([
            'WARNING: Some cherry pick plates will be created from multiple ',
            'source plates of non-uniform plate types!',
            'The following cherry pick plates are specified across multiple files:',])

        reload_plates_warning = '\n'.join([
            'WARNING: Some cherry pick plates will be created from the same source plate!',
            'You will need to reload one or more source plates for each of the ',
            'following cherry pick plates:',])
        reload_plate_msg = 'Cherry pick plate %s requires reload of source plate: %s'

        
        readme_text = [
            'This zip file contains plate mappings for Cherry Pick Request %r' 
            % str(cpr_id) ]
        readme_text.append('Cherry pick plates:')
        plating_text = '{{filename}} Plated ({plating_date} by {plated_by_name})'
        # Open with delete=False; file will be closed and deleted 
        # by the FileWrapper1 instance when streaming is finished.
        zip_dir_name = \
            'screen_{screen_facility_id}_cp_{cherry_pick_request_id}'.format(**cpr)
        with  NamedTemporaryFile(delete=False) as temp_file:
            with ZipFile(temp_file, 'w') as zipfile:
                    
                # Plate maps
                for plate_name,lcps in plate_map.items():
                    raw_data = cStringIO.StringIO()
                    writer = unicodecsv.writer(
                        raw_data, encoding='utf-8',lineterminator='\r\n')
                    for i,lcp in enumerate(lcps):
                        if i==0:
                            title_row = [title_function(x) for x in columns_to_include]
                            title_row.extend(extra_values.keys())
                            writer.writerow(title_row)
                        values = []
                        for key in columns_to_include:
                            values.append(vocab_function(key, lcp[key]))
                        values.extend(extra_values.values())
                        writer.writerow(values)
                    filename = '%s.csv' % plate_name
                    logger.info('write; %r', plate_name)
                    cherry_pick_plate_number = lcps[0]['cherry_pick_plate_number']
                    cpp = cpps[cherry_pick_plate_number]
                    if cpp['plating_date']:
                        _text = plating_text.format(**cpp)
                        _text = _text.format(filename=filename)
                        readme_text.append(_text)
                    else:
                        readme_text.append(filename)
                    
                    zipi= ZipInfo()
                    # give full access to included file:
                    # NOTE: second high bit is used in external_attr (4 bytes)
                    zipi.external_attr = 0777 << 16L 
                    zipi.filename= '%s/%s' % (zip_dir_name,filename)
                    zipfile.writestr(zipi, raw_data.getvalue())
                
                # Location map
                raw_data = cStringIO.StringIO()
                writer = unicodecsv.writer(raw_data)
                writer.writerow(['Source Plate','Source Copy','Location'])
                for location in locations:
                    writer.writerow(location)
                zipi= ZipInfo()
                zipi.filename= '%s/%s' % (zip_dir_name,'plate-copy-location.csv')
                zipi.external_attr = 0777 << 16L
                zipfile.writestr(zipi, raw_data.getvalue())
                
                # Readme
                extra_plate_messages = []
                for plate_name,types in plate_name_to_types.items():
                    if len(types)>1:
                        extra_plate_messages.append(plate_name)
                if extra_plate_messages:
                    logger.info('extra_plate_messages: %r', extra_plate_messages)
                    extra_plate_messages.insert(0,readme_types_warning)
                    readme_text.append('\n')
                    readme_text.extend(extra_plate_messages)
                
                split_plate_messages = []
                dest_plates_to_source_plates = defaultdict(set)
                for plate_copy,dest_plates in source_plate_to_dest_plates.items():
                    if len(dest_plates) > 1:
                        for dest_plate in dest_plates:
                            logger.info('dest_plate: %r, has source plate: %r',
                                dest_plate,plate_copy)
                            dest_plates_to_source_plates[dest_plate].add(plate_copy)
                for dest_plate in sorted(dest_plates_to_source_plates.keys()):
                    source_plates = dest_plates_to_source_plates[dest_plate]
                    for source_plate in source_plates:
                        split_plate_messages.append(
                            reload_plate_msg % (dest_plate,source_plate))
                if split_plate_messages:
                    logger.info('split_plate_messages: %r', split_plate_messages)
                    split_plate_messages.insert(0,reload_plates_warning)
                    readme_text.append('\n')
                    readme_text.extend(split_plate_messages)
                            
                zipi= ZipInfo()
                zipi.filename= '%s/%s' % (zip_dir_name,'readme.txt')
                zipi.external_attr = 0777 << 16L
                zipfile.writestr(zipi, '\n'.join(readme_text))
            logger.info('wrote file %r', temp_file)
        
            temp_file.seek(0, os.SEEK_END)
            size = temp_file.tell()
            temp_file.seek(0)   
        
        filename = 'testCPR%d_plating_file.zip' % cpr['cherry_pick_request_id']
        logger.info('download zip file: %r',filename)
        _file = file(temp_file.name)
        response = StreamingHttpResponse(FileWrapper1(_file)) 
        response['Content-Length'] = size
        response['Content-Type'] = '%s; charset=utf-8' % ZIP_MIMETYPE
        response['Content-Disposition'] = \
            'attachment; filename=%s' % filename
        return response
        
    # TODO: implement the plate mapping using standard resource endpoints for 
    # plate_mapping view (of lab cherry pick) and plate location view (of cpr: tbi)        
    #         kwargs = {}
    #         kwargs['limit'] = 0
    #         kwargs['cherry_pick_request_id'] = cpr_id
    #         kwargs['includes'] = columns_to_include
    #                 for cherry_pick_plate_number in range(1,cpr['number_plates']+1):
    #                     kwargs['cherry_pick_plate_number'] = cherry_pick_plate_number
    #                     response = self.get_labcherrypick_resource().get_list(
    #                         request,
    #                         format='csv',
    #                         **kwargs)
    #                     logger.debug('response: %r', response)
    #                     
    #                     name = ((
    #                         'CP{{cherry_pick_request_id}}'
    #                         '_Screen{{screen_facility_id}}'
    #                         '_plate_{cherry_pick_plate_number}'
    #                         '_of_{{number_plates}}_{{requested_by_username}}.csv'
    #                         ).format(cherry_pick_plate_number=cherry_pick_plate_number)
    #                         .format(**cpr))
    #                     zipfile.writestr(name, LimsSerializer.get_content(response))

    def get_scp_schema(self, request, **kwargs):
        return self.get_screenercherrypick_resource().get_schema(request, **kwargs)    

    def get_lcp_schema(self, request, **kwargs):
        return self.get_labcherrypick_resource().get_schema(request, **kwargs)    
        
    @read_authorization
    def dispatch_screener_cherry_pick_view(self, request, **kwargs):
        return self.get_screenercherrypick_resource().dispatch('list', request, **kwargs)
    
    @read_authorization
    def dispatch_lab_cherry_pick_view(self, request, **kwargs):
        return self.get_labcherrypick_resource().dispatch('list', request, **kwargs)    

    @read_authorization
    def dispatch_lab_cherry_pick_plating(self, request, **kwargs):
        ''' 
        Show the lab cherry pick view after plating
        (Method used by the UI)
        - modify schema and colums to show plate-mapping fields 
        { "cherry_pick_assay_plate", "destination_well", etc. }
        - modify field ordering
        - filter: status=='plated' (TODO: implement "show_all")
        '''

        schema = self.get_labcherrypick_resource()\
            .build_lab_cherry_pick_plating_schema(request.user, **kwargs)
        kwargs['plating_schema'] = schema
        
        param_hash = self._convert_request_to_dict(request)

        show_copy_wells = parse_val(
            param_hash.get(API_PARAM_SHOW_COPY_WELLS, False),
            API_PARAM_SHOW_COPY_WELLS, 'boolean')
        logger.info('%r: %r', API_PARAM_SHOW_COPY_WELLS,show_copy_wells)
        show_available_and_retired_copy_wells = parse_val(
            param_hash.get(API_PARAM_SHOW_RETIRED_COPY_WELlS, False),
            API_PARAM_SHOW_RETIRED_COPY_WELlS, 'boolean')
        logger.info('%r: %r', 
            API_PARAM_SHOW_RETIRED_COPY_WELlS,show_available_and_retired_copy_wells)
        show_unfulfilled = parse_val(
            param_hash.get(API_PARAM_SHOW_UNFULFILLED, False),
            API_PARAM_SHOW_UNFULFILLED, 'boolean')
        if not any([show_copy_wells, show_available_and_retired_copy_wells, 
            show_unfulfilled]):
            kwargs['status__eq'] = 'plated'
            kwargs['order_by'] = param_hash.get(
                'order_by', ['cherry_pick_plate_number','destination_well'])
            kwargs['includes'] = param_hash.get(
                'includes', ['-library_plate_comment_array'])
        return self.get_labcherrypick_resource()\
            .dispatch('list',request, **kwargs)    

    @read_authorization
    def get_lab_cherry_pick_plating_schema(self, request, **kwargs):
        return self.get_labcherrypick_resource()\
            .get_lab_cherry_pick_plating_schema(request, **kwargs)    

    def dispatch_cherry_pick_plate_view(self, request, **kwargs):
        return self.get_cherrypickplate_resource()\
            .dispatch('list', request, **kwargs)    

    @read_authorization
    def get_detail(self, request, **kwargs):

        cherry_pick_request_id = kwargs.get('cherry_pick_request_id', None)
        if not cherry_pick_request_id:
            raise Http404('must provide a cherry_pick_request_id parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, schema=None, **kwargs):
        
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        if schema is None:
            raise Exception('schema not initialized')
        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id:
            param_hash['cherry_pick_request_id__eq'] = cherry_pick_request_id

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            manual_field_includes.add('has_pool_screener_cherry_picks')
            manual_field_includes.add('has_alternate_screener_cherry_pick_selections')
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
            
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)

            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
                    
            # specific setup 
            base_query_tables = ['cherry_pick_request', 'screen']
            _cpr = self.bridge['cherry_pick_request']
            _screen = self.bridge['screen']
            _su = self.bridge['screensaver_user']
            _lhsu = _su.alias('lhsu')
            _scp = self.bridge['screener_cherry_pick']
            _sirna_pool_duplexes = self.bridge['silencing_reagent_duplex_wells']
            _well = self.bridge['well']
            _library = self.bridge['library']
            _lcp = self.bridge['lab_cherry_pick']
            _cpp = self.bridge['cherry_pick_assay_plate']
            _cpp2 = _cpp.alias('cpp2')
            
            _user_cte = ScreensaverUserResource.get_user_cte().cte('users_cpr')
            _requestors = _user_cte.alias('requestors')
            _approvers = _user_cte.alias('approvers')
            _lead_screeners = _user_cte.alias('lead_screeners')
            lab_head_table = ScreensaverUserResource.get_lab_head_cte().cte('lab_heads')
            
            _lcp_subquery = (
                select([
                    _lcp.c.cherry_pick_request_id,
                    func.count().label('lcp_count'),
                    func.count(_lcp.c.copy_id).label('lcp_fullfilled_count')])
                .select_from(_lcp)
                .group_by(_lcp.c.cherry_pick_request_id)).cte('lcp_subquery')
                
            custom_columns = {
                'requested_by_id': _requestors.c.screensaver_user_id,
                'requested_by_name': _requestors.c.name,
                'volume_approved_by_username': _approvers.c.username,
                'volume_approved_by_name': _approvers.c.name,
                'lab_name': lab_head_table.c.lab_name_full,
                'lab_head_id': lab_head_table.c.screensaver_user_id,
                'lab_head_username': lab_head_table.c.username,
                'lead_screener_name': _lead_screeners.c.name,
                'lead_screener_id': _lead_screeners.c.screensaver_user_id,
                'lead_screener_username': _lead_screeners.c.username,
                'number_plates': (
                    select([func.count(None)])
                        .select_from(_cpp)
                        .where(_cpp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id)
                    ),
                'number_plates_completed': (
                    select([func.count(None)])
                        .select_from(_cpp)
                        .where(_cpp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id)
                        .where(_cpp.c.plating_date!=None)
                    ),
                'number_plates_screened': (
                    select([func.count(None)])
                        .select_from(_cpp)
                        .where(_cpp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id)
                        .where(_cpp.c.screening_date!=None)
                    ),
                'is_completed': (
                    case(
                        [(
                            (select([func.count(None)])
                             .select_from(_cpp)
                             .where(_cpp.c.cherry_pick_request_id==
                                 _cpr.c.cherry_pick_request_id).as_scalar()>0), 
                            (select([func.count(None)])
                             .select_from(_cpp)
                             .where(_cpp.c.cherry_pick_request_id
                                 ==_cpr.c.cherry_pick_request_id)
                             .where(_cpp.c.plating_date==None)).as_scalar()==0),
                        ],
                        else_=text('false'))
                    ),
                'number_unfulfilled_lab_cherry_picks': (
                    func.coalesce(
                        cast(_lcp_subquery.c.lcp_count - 
                                _lcp_subquery.c.lcp_fullfilled_count,
                             sqlalchemy.sql.sqltypes.Integer),0)),                    
                'total_number_lcps': func.coalesce(
                    cast(_lcp_subquery.c.lcp_count,
                        sqlalchemy.sql.sqltypes.Integer),0),
                'last_plating_activity_date': (
                    select([func.max(_cpp.c.plating_date)])
                    .select_from(_cpp)
                    .where(_cpp.c.cherry_pick_request_id
                        ==_cpr.c.cherry_pick_request_id)),
                'last_screening_activity_date': (
                    select([func.max(_cpp.c.screening_date)])
                    .select_from(_cpp)
                    .where(_cpp.c.cherry_pick_request_id
                        ==_cpr.c.cherry_pick_request_id)),
                
                # TODO: new: when the lcp's were reserved and mapped
                # 'date_volume_reserved': literal_column("'2016-12-07'"),
                
                'screener_cherry_picks': (
                    select([func.array_to_string(
                        func.array_agg(literal_column('screened_well_id')), 
                        LIST_DELIMITER_SQL_ARRAY) ])
                    .select_from(
                        select([_scp.c.screened_well_id])
                            .select_from(_scp)
                            # NOTE: when doing an inner select, must use literal_column,
                            # otherwise SQalchemy thinks it is another table to add
                            .where(_scp.c.cherry_pick_request_id
                                ==literal_column(
                                    'cherry_pick_request.cherry_pick_request_id'))
                            .where(_scp.c.selected)
                            .order_by('screened_well_id').alias('inner_cpwells'))
                        ),
                'screener_cherry_pick_count': (
                    select([func.count(None)])
                        .select_from(_scp)
                        .where(_scp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id)
                        .where(_scp.c.selected)),
                'has_pool_screener_cherry_picks': (
                    select([func.count(None)>0])
                        .select_from(
                            _scp.join(_well, 
                                _scp.c.screened_well_id==_well.c.well_id)
                            .join(_library, 
                                _well.c.library_id==_library.c.library_id))
                        .where(_library.c.is_pool==True)
                        .where(_scp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id).limit(1)),
                'has_alternate_screener_cherry_pick_selections': (
                    select([func.count(None)>0])
                        .select_from(_scp)
                        .where(_scp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id)
                        .where(_scp.c.searched_well_id!=_scp.c.screened_well_id)
                    )
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            j = join(_cpr, _screen, _cpr.c.screen_id == _screen.c.screen_id)
            j = j.join(
                _requestors, 
                _cpr.c.requested_by_id==_requestors.c.screensaver_user_id,
                isouter=True)
            j = j.join(
                _approvers, 
                _cpr.c.volume_approved_by_id==_approvers.c.screensaver_user_id,
                isouter=True)
            j = j.join(
                lab_head_table, 
                lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id,
                isouter=True)
            j = j.join(
                _lead_screeners,
                _lead_screeners.c.screensaver_user_id==_screen.c.lead_screener_id,
                isouter=True)

            if ('number_unfulfilled_lab_cherry_picks' in field_hash
                 or 'total_number_lcps' in field_hash):
                j = j.join(
                    _lcp_subquery, _cpr.c.cherry_pick_request_id
                        ==_lcp_subquery.c.cherry_pick_request_id, isouter=True)
            
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
#             compiled_stmt = str(stmt.compile(
#                 dialect=postgresql.dialect(),
#                 compile_kwargs={"literal_binds": True}))
#             logger.info('compiled_stmt %s', compiled_stmt)
            
            if not order_clauses:
                stmt = stmt.order_by(
                    nullslast(desc(column('cherry_pick_request_id'))))
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
             
        except Exception, e:
            logger.exception('on get list')
            raise e  
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def post_detail(self, request, **kwargs):
        # FIXME: Set the log URI using the containing screen URI
        return DbApiResource.post_detail(
            self, request, full_create_log=True, **kwargs)

    @write_authorization
    @un_cache  
    @transaction.atomic      
    def patch_detail(self, request, **kwargs):
        '''
        Override to generate informational summary for callee
        '''
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        deserialized = kwargs.pop('data', None)
        # allow for internal data to be passed
        if deserialized is None:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))
        
        if API_RESULT_DATA in deserialized:
            deserialized = deserialized[API_RESULT_DATA][0]
            
        logger.debug('patch detail %s, %s', deserialized,kwargs)

        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(
            deserialized, schema=schema, validate=True,**kwargs)

        original_data = None
        if kwargs_for_log:
            try:
                original_data = self._get_detail_response_internal(**kwargs_for_log)
                kwargs['original_data'] = original_data
            except Exception, e: 
                logger.exception('exception when querying for existing obj: %s', 
                    kwargs_for_log)
        try:
            parent_log = kwargs.get('parent_log', None)
            log = self.make_log(request)
            log.parent_log = parent_log
            log.save()
            kwargs['parent_log'] = log
            patch_response = self.patch_obj(request, deserialized, **kwargs)
        except ValidationError as e:
            logger.exception('Validation error: %r', e)
            raise e

        # get new state, for logging
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        logger.debug('original: %r, new: %r', original_data, new_data)
        log = self.log_patch(
            request, original_data,new_data,log=log, full_create_log=True, 
            excludes=['screener_cherry_picks'],
            **kwargs_for_log)
        # Set the log URI using the containing screen URI
        if log:
            log.uri = '/'.join([
                'screen',new_data['screen_facility_id'],
                log.ref_resource_name,log.key])
            log.save()
            logger.debug('log info: %r, %r', log, log.diffs )
        
        meta = {}
        if API_RESULT_META in patch_response:
            meta = patch_response[API_RESULT_META]
        return self.build_response(
            request,  { API_RESULT_META: meta }, 
            response_class=HttpResponse, **kwargs)
            
    @write_authorization
    @un_cache  
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        '''
        Create a new Cherry Pick Request
        - set the screener cherry picks if included
        '''
        
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        initializer_dict = {}

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        
        logger.info('patch CPR: %r', deserialized)
        
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.debug('param_hash: %r', param_hash)
        
        patch = bool(id_kwargs)
        initializer_dict = self.parse(deserialized, schema=schema, create=not patch)
        errors = self.validate(initializer_dict, schema=schema, patch=patch)
        if errors:
            raise ValidationError(errors)

        cpr = None
        if patch is True:
            try:
                cpr = CherryPickRequest.objects.get(**id_kwargs)
            except ObjectDoesNotExist:
                raise Http404(
                    'Cherry Pick Request does not exist for: %r', id_kwargs)
        
        if patch is not True:
            _key = 'screen_facility_id'
            _val = deserialized.get(_key, None)
            if _val is None:
                _val = kwargs.get(_key,None)
            if _val is None:
                raise ValidationError(
                    key='screen_facility_id',msg='required')
            try:
                screen = Screen.objects.get(facility_id=_val)
                initializer_dict['screen'] = screen
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='does not exist: {val}'.format(val=_val))
        else:
            screen = cpr.screen
            
        _key = 'requested_by_id'
        _val = deserialized.get(_key, None)
        if _val:
            try:
                requested_by_user = ScreensaverUser.objects.get(screensaver_user_id=_val)
                if requested_by_user not in screen.get_screen_users():
                    raise ValidationError(
                        key='requested_by_id',
                        msg='"%s" must be one of the screen users: %r'
                            %(_val, [x.screensaver_user_id 
                                for x in screen.get_screen_users()]))
                initializer_dict['requested_by'] = requested_by_user
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='does not exist: {val}'.format(val=_val))

        _key = 'volume_approved_by_username'
        _val = deserialized.get(_key, None)
        if _val:
            try:
                volume_approved_by_user = ScreensaverUser.objects.get(username=_val)
                
                if not self._meta.authorization._is_resource_authorized(
                        volume_approved_by_user.user.user,'write'):
                    raise ValidationError(
                        key='volume_approved_by_username',
                        msg='user %r does not have %r %r authorization' 
                            % (_val, self._meta.resource_name, 'write'))
                self.validate_volume_approver(volume_approved_by_user)
                initializer_dict['volume_approved_by'] = volume_approved_by_user
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='does not exist: {val}'.format(val=_val))
        logger.info('initializer_dict: %r', initializer_dict)
        _meta = {}
        try:
            if patch is not True:
                cpr = CherryPickRequest()

            allocated_lcps_query = cpr.lab_cherry_picks.filter(
                cherry_pick_assay_plate__isnull=False)
            if allocated_lcps_query.exists():
                disallowed_fields = [
                    'assay_plate_type',
                    'keep_source_plate_cherry_picks_together'
                    'is_randomized_assay_plate_layout'
                    'wells_to_leave_empty'
                    'transfer_volume_per_well_approved']
                if set(disallowed_fields) & set(initializer_dict.keys()):
                    raise ValidationError(
                        key=API_MSG_LCPS_MUST_BE_DELETED,
                        msg='Disallowed fields: %r' % disallowed_fields)
            model_field_names = [
                x.name for x in cpr._meta.get_fields()]
            for key, val in initializer_dict.items():
                if key in model_field_names:
                    setattr(cpr, key, val)

            cpr.save()
            logger.info('patch cpr created: %r', cpr)
            
            final_warn_msg = []
            wells_to_leave_empty = deserialized.get('wells_to_leave_empty', None)
            if wells_to_leave_empty is not None:
                parsed_wells_to_leave_empty = \
                    lims_utils.parse_wells_to_leave_empty(
                        wells_to_leave_empty, cpr.assay_plate_size)
                cpr.wells_to_leave_empty = ', '.join(parsed_wells_to_leave_empty)
                cpr.save()
                
            if 'screener_cherry_picks' in deserialized:
                plated_assay_plates_query = \
                    cpr.cherry_pick_assay_plates.filter(
                        plating_date__isnull=False)                
                if plated_assay_plates_query.exists():
                    raise ValidationError({
                        API_MSG_NOT_ALLOWED: 
                            ('Screener cherry picks may not be reassigned after '
                             'plates have been plated'),
                        API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count(),
                        API_MSG_LCP_ASSAY_PLATES_PLATED: plated_assay_plates_query.count()
                })
                
                if cpr.lab_cherry_picks.exists():
                    raise ValidationError({
                        'total_number_lcps': 
                            ('Lab cherry picks already assigned: (%d); '
                            'delete lab cherry picks to change '
                            'screener cherry pick selections' 
                            % cpr.lab_cherry_picks.count()),
                        API_MSG_LCPS_MUST_BE_DELETED: cpr.lab_cherry_picks.count()
                    })
                    
                # TODO: Test override
                override_param = parse_val(
                    param_hash.get(API_PARAM_OVERRIDE, False),
                        API_PARAM_OVERRIDE, 'boolean')
                
                raw_screener_cps = deserialized['screener_cherry_picks']
                if raw_screener_cps is None or len(raw_screener_cps) == 0:
                    logger.info('removing screener_cherry_picks')
                    # FIXME: remove lab_cherry_picks, iif:
                    # - not plated
                    # - replace allocate well volumes
                    # - create a log for the action
                    _meta[API_MSG_SCPS_DELETED] = cpr.screener_cherry_picks.all().count()
                    cpr.screener_cherry_picks.all().delete()
                    cpr.save()
                else:
                    cherry_pick_wells = self.find_wells(raw_screener_cps)
                    logger.info(
                        'found screener_cherry_picks: %d', len(cherry_pick_wells))
                    not_allowed_libraries = set()
                    discarded_libraries = set()
                    wrong_screen_type = set()
                    for well in cherry_pick_wells:
                        screen_type = cpr.screen.screen_type
                        if well.library.screen_type != screen_type:
                            wrong_screen_type.add(well.well_id)
                        if wrong_screen_type:
                            continue
                        if well.library.screening_status == 'discarded':
                            discarded_libraries.add(well.library)
                        if well.library.screening_status != 'allowed':
                            not_allowed_libraries.add(well.library)
                        if len(not_allowed_libraries)>0 and override_param is not True:
                            continue
                        screener_cherry_pick = ScreenerCherryPick.objects.create(
                            cherry_pick_request=cpr,
                            screened_well=well,
                            searched_well=well,
                            selected=True)
                    if wrong_screen_type:
                        wrong_screen_type = sorted(wrong_screen_type)
                        raise ValidationError(
                            key='Wrong Well screen_type, must be %s' 
                                % cpr.screen.screen_type,
                            msg=wrong_screen_type)
                    if discarded_libraries:
                        discarded_libraries = sorted([
                            '%s - status: %s' % (l.short_name,l.screening_status)
                                for l in discarded_libraries])
                        raise ValidationError({
                            'screener_cherry_picks': 'Discarded libraries',
                            'Libraries': discarded_libraries
                            })
                        
                    if not_allowed_libraries:
                        not_allowed_libraries = sorted([
                            '%s - status: %s' % (l.short_name,l.screening_status)
                                for l in not_allowed_libraries])
                    if len(not_allowed_libraries)>0 and override_param is not True:
                        raise ValidationError({
                            API_PARAM_OVERRIDE: 'required',
                            'screener_cherry_picks': (
                                'Override required to screen libraries that are '
                                'not allowed'),
                            'Libraries': not_allowed_libraries
                            }
                        )
                    if len(not_allowed_libraries)>0:
                        final_warn_msg.append(
                            ('Override used for libraries', not_allowed_libraries))
                    _meta[API_MSG_SCPS_CREATED] = cpr.screener_cherry_picks.all().count()    

            if final_warn_msg:
                _meta[API_MSG_WARNING] = final_warn_msg
                            
            response = { API_RESULT_OBJ: cpr, API_RESULT_META: _meta }
            logger.info('response: %r', response)
            return response
        except Exception, e:
            logger.exception('on patch_obj')
            raise e
    
    @staticmethod
    def _get_plate_size(assay_plate_type):
        parts = assay_plate_type.split('_')
        if len(parts) != 2:
            raise ValidationError(
                key='assay_plate_type', 
                msg='not a recognized type: %r' % assay_plate_type)
        plate_size = int(parts[1])
        return plate_size

    
    @classmethod
    def find_wells(cls, cherry_pick_well_patterns ):
        logger.debug('find wells for patterns: %r', cherry_pick_well_patterns)
        if not isinstance(cherry_pick_well_patterns, (list,tuple)):
            cherry_pick_well_patterns = (cherry_pick_well_patterns,)

        wells = set()
        (wells,errors) = WellResource.find_wells(cherry_pick_well_patterns)
        if errors:
            raise ValidationError(
                key='screener_cherry_picks',
                msg='%r' % errors )    
        
        non_experimental_wells = []
        for well in wells:
            if well.library_well_type != 'experimental':
                non_experimental_wells.append(well)
        if non_experimental_wells:
            raise ValidationError(   
                key='Can not screen non-experimental wells',
                msg=', '.join([well.well_id for well in non_experimental_wells]))
            
        logger.debug('found wells: %r', wells)
        return wells

    def validate_cpr_for_plating(self, cpr):
        logger.info('validating: %r', cpr)      
        if not cpr.transfer_volume_per_well_approved:
            raise ValidationError(
                key='transfer_volume_per_well_approved',
                msg='required')
        self.validate_volume_approver(cpr.volume_approved_by)
        
    def validate_volume_approver(self, user):
        logger.info('validating: %r', user)      
        if not user:
            raise ValidationError(
                key='volume_approved_by_username',
                msg='required')
        else:
            # find the user, verify that they have cherrypickrequest admin permission
            user_data = \
                self.get_user_resource()._get_detail_response_internal(**{
                'username': user.username 
                })
            if not user_data:
                raise ValidationError(
                    key='volume_approved_by_username', msg='not a valid user')
            if user_data.get('is_superuser', False) is True:
                return
            required_permission = 'resource/cherrypickrequest/write'
            if required_permission not in user_data['all_permissions']:
                logger.warn('user: %r, permissions: %r, does not contain: %r',
                    user_data['username'], user_data['all_permissions'], 
                    required_permission)
                raise ValidationError(
                    key='volume_approved_by_username', msg='user is not authorized')
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def dispatch_reserve_map_lab_cherry_picks(self, request, **kwargs):
        '''
        Reserve (allocate) the Copy Well volume for the Lab Cherry Picks 
        on the Cherry Pick Request:
        - only consider fulfilled LCPs (where a source copy has been assigned) 
        '''
        
        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')
        convert_post_to_put(request)
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        cherry_pick_request_id = param_hash['cherry_pick_request_id']
        logger.info(
            'dispatch_reserve_map_lab_cherry_picks for: %r...', 
            cherry_pick_request_id)
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cherry_pick_request_id)
        if not cpr.lab_cherry_picks.filter(copy__isnull=False).exists():
            raise ValidationError(
                key='total_number_lcps', 
                msg='No (fulfilled) Lab Cherry Picks have been created')
        
        self.validate_cpr_for_plating(cpr) 
        
        # Create parent log: cherry pick request: date_volume_reserved
        parent_log = self.make_log(request)
        parent_log.key = str(cpr.cherry_pick_request_id)
        parent_log.uri = '/'.join([
            'screen',cpr.screen.facility_id, 
            parent_log.ref_resource_name,parent_log.key])        
        parent_log.save()
        
        previous_date_reserved = cpr.date_volume_reserved
        cpr.date_volume_reserved = _now().date() 
        cpr.save()
        parent_log.diffs = {
            'date_volume_reserved': [previous_date_reserved, cpr.date_volume_reserved ]}
        
        status_messages = []
        copywell_deallocation_meta = None
        previous_number_of_plates = None           
        
        logger.info('check for previous plating assignments...')
        allocated_lcps_query = cpr.lab_cherry_picks.filter(
            cherry_pick_assay_plate__isnull=False)
        if allocated_lcps_query.exists():
            plated_assay_plates_query = \
                cpr.cherry_pick_assay_plates.filter(
                    plating_date__isnull=False)                
            if plated_assay_plates_query.exists():
                raise ValidationError({
                    API_MSG_NOT_ALLOWED: API_MSG_CPR_PLATED_CANCEL_DISALLOWED, 
                    API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count(),
                    API_MSG_LCP_ASSAY_PLATES_PLATED: plated_assay_plates_query.count()
            })
            raise ValidationError({
                API_MSG_NOT_ALLOWED: 
                    ('Lab cherry pick plates already assigned; '
                    'reservation must be canceled before reassignment is allowed'),
                API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count()
            })
        
        logger.info('Find the fulfillable lab cherry picks...')
        available_assay_plate_wells = cpr.assay_plate_available_wells
        logger.debug('available_assay_plate_wells: %r',available_assay_plate_wells)
        logger.info('available_assay_plate_wells len: %d', 
            len(available_assay_plate_wells))

        fulfillable_lcps = (
            cpr.lab_cherry_picks.filter(copy__isnull=False)
                .order_by('source_well__plate_number','copy__name') )

        logger.info('re-Check and reserve copy volumes...')
        lab_cherry_pick_copywells = \
            { lcp['source_copywell_id']: lcp for lcp in
                self.get_labcherrypick_resource()._get_list_response_internal(
                    **{
                        'cherry_pick_request_id': cpr.cherry_pick_request_id,
                        'source_copy_name__is_null': False,
                        'includes': [
                            'source_plate_type','destination_plate_type',
                            'source_copywell_id','source_copy_well_volume',
                            'volume_approved',
                            '-structure_image','-molfile','-library_plate_comment_array'],
                    })}
        logger.info('fetch output readable format...')
        lab_cherry_pick_copywells_output = \
            { lcp['source_copywell_id']: lcp for lcp in
                self.get_labcherrypick_resource()._get_list_response_internal(
                    **{
                        'cherry_pick_request_id': cpr.cherry_pick_request_id,
                        'source_copy_name__is_null': False,
                        'includes': [
                            'source_plate_type','destination_plate_type',
                            'source_copywell_id','source_copy_well_volume',
                            'volume_approved',
                            '-structure_image','-molfile','-library_plate_comment_array'],
                        HTTP_PARAM_USE_VOCAB: True,
                    })}
        logger.info('Check for insufficient well volumes...')
        unfulfillable_wells = []
        override_well_volume = parse_val(
            param_hash.get(API_PARAM_VOLUME_OVERRIDE, False),
            API_PARAM_VOLUME_OVERRIDE, 'boolean')
        
        for copywell_id,lcp_data in lab_cherry_pick_copywells.items():
            lcp_scwv = Decimal(lcp_data['source_copy_well_volume'])
            logger.debug('consider lcp_cw: %r, %r to %r', 
                copywell_id, lcp_scwv, cpr.transfer_volume_per_well_approved)
            if ( lcp_scwv < cpr.transfer_volume_per_well_approved ):
                lcp_output = lab_cherry_pick_copywells_output[copywell_id]
                unfulfillable_wells.append(
                    (copywell_id, 
                        '(available: %s uL)' % (lcp_output['source_copy_well_volume'] or 0),
                        '(requested: %s uL)' % lcp_output['volume_approved']))
        unfulfillable_wells = sorted(unfulfillable_wells, key=lambda x: x[0])
        if unfulfillable_wells and override_well_volume is not True:
            raise ValidationError({
                'transfer_volume_per_well_approved':
                    '%s: %r ' % (API_MSG_LCPS_INSUFFICIENT_VOLUME, 
                                 unfulfillable_wells),
                API_MSG_LCPS_INSUFFICIENT_VOLUME: unfulfillable_wells,
                API_PARAM_VOLUME_OVERRIDE: 'required',
                
            })
        elif unfulfillable_wells:
            logger.info(
                'Override: unfulfillable wells: %r' % unfulfillable_wells)
        else:
            logger.info('all wells are fulfillable')

        logger.info('Reserve copy well volumes...')
        copywell_reservation_meta = \
            self.get_copywell_resource().reserve_cherry_pick_volumes(
                cpr, fulfillable_lcps, parent_log)
        
        logger.info('Create the assay_plates...')  
#         # - if true, randomize the plate layout wells
#         if cpr.is_randomized_assay_plate_layout:
#             logger.debug('randomize the available assay plate wells: %r',
#                 available_assay_plate_wells)
#             random.shuffle(available_assay_plate_wells)
#             logger.debug('randomized available assay plate wells: %r',
#                 available_assay_plate_wells)
                  
        next_plate_ordinal = [1]
        def create_next_assay_plate():
            assay_plate = CherryPickAssayPlate.objects.create(
                cherry_pick_request=cpr,
                plate_ordinal=next_plate_ordinal[0],
                # TODO: deprecate attempt_ordinal
                attempt_ordinal=0,
                assay_plate_type=cpr.assay_plate_type,
                # TODO: deprecate cpapt
                cherry_pick_assay_plate_type='CherryPickAssayPlate',
            )
            next_plate_ordinal[0] = next_plate_ordinal[0]+1
            return assay_plate
        
        # Create a sortable multimap
        lcp_by_plate_copy = defaultdict(list)
        for lcp in fulfillable_lcps.all():
            plate_copy = '%s:%s' % (
                str(lcp.source_well.plate_number),lcp.copy.name)
#                 str(lcp.source_well.plate_number).zfill(5),lcp.copy.name)
            lcp_by_plate_copy[plate_copy].append(lcp)
        # Sort internally
        for plate_copy in lcp_by_plate_copy.keys():
            lcps = lcp_by_plate_copy[plate_copy]
            lcps = sorted(lcps, key=lambda lcp: lcp.source_well_id)
            lcp_by_plate_copy[plate_copy] = lcps
            
        if cpr.keep_source_plate_cherry_picks_together is True:

            logger.info("keep_source_plate_cherry_picks_together...")
            logger.info("use the bin packer to fit the lcp's to bins...")
            capacity = len(available_assay_plate_wells)
            packages = [ { 'name': plate_copy, 'size': len(lcps) } 
                for plate_copy,lcps in lcp_by_plate_copy.items() ]
            packed_bins = bin_packer.pack_bins(capacity, packages)
            packed_bins = sorted(
                packed_bins,
                key=lambda bin: (capacity-bin_packer.sum_bin(bin),bin[0]['name']),
                reverse=False)
            logger.info('packed bins: %r', 
                [ (bin_packer.sum_bin(bin),bin[0]['name']) for bin in packed_bins])
            
            logger.info(
                "assign the packed_bins to assay_plates, "
                "order by sum_bin, librarycopyplate.key")
            ordered_bins = []
            for packed_bin in packed_bins:
                if packed_bin not in ordered_bins:
                    ordered_bins.append(packed_bin)
                    partially_packed_plate_copy = [
                        package['name'] for package in packed_bin
                            if (len(lcp_by_plate_copy[package['name']]) 
                                > package['size']) ]
                    if len(partially_packed_plate_copy) > 1:
                        raise ProgrammingError(
                            'Packed bins can not contain more than one '
                            'partially packed plate_copies: %r' 
                            % partially_packed_plate_copy)
                    # NOTE: assume that source plate will never need > 2 assay plates
                    if partially_packed_plate_copy:
                        plate_copy = partially_packed_plate_copy[0]
                        for second_bin in packed_bins:
                            if second_bin != packed_bin:
                                if plate_copy in [ p['name'] for p in second_bin]:
                                    ordered_bins.append(second_bin)
            
            logger.info('create assay plates, in order...')
            assay_plates_created = []

            for packed_bin in ordered_bins:
            
                logger.info('using packed bin: %r', packed_bin)
                assay_plate_wells_to_use = \
                    available_assay_plate_wells[:bin_packer.sum_bin(packed_bin)]                                
                if cpr.is_randomized_assay_plate_layout:                    
                    random.shuffle(assay_plate_wells_to_use)
                assay_plate = create_next_assay_plate()
                assay_plates_created.append(assay_plate)
                assay_plate_well_index = 0
                    
                for package in packed_bin:

                    logger.info('package: %r', package)
                    
                    plate_copy = package['name']
                    size = package['size']
                    plate_copy_lcps = lcp_by_plate_copy[plate_copy]
                    lcps_to_plate = plate_copy_lcps[:size]
                    if size < len(plate_copy_lcps):
                        lcp_by_plate_copy[plate_copy] = plate_copy_lcps[size:]
                    
                    logger.info('lcps: %r, assay_plate.plate_ordinal: %r', 
                        len(lcps_to_plate),assay_plate.plate_ordinal)
                    
                    
                    for lcp in lcps_to_plate:
                        lcp.cherry_pick_assay_plate = assay_plate
                        well_name = assay_plate_wells_to_use[assay_plate_well_index]
                        lcp.assay_plate_row = lims_utils.well_name_row_index(well_name)
                        lcp.assay_plate_column = lims_utils.well_name_col_index(well_name)
                        assay_plate_well_index += 1
                        lcp.save()
        else: 
            logger.info('create assay plates, in order...')
            logger.info("do not keep_source_plate_cherry_picks_together")
            # keep_source_plate_cherry_picks_together == False
            # no bin packing
            capacity = len(available_assay_plate_wells)

            def randomize_lcps(lcps):
                current_well_assignments = [
                    (lcp.assay_plate_row,lcp.assay_plate_column) 
                    for lcp in lcps]
                random.shuffle(current_well_assignments)
                for i,lcp in enumerate(lcps):
                    shuffled_rc = current_well_assignments[i]
                    lcp.assay_plate_row = shuffled_rc[0]
                    lcp.assay_plate_column = shuffled_rc[1]
                    lcp.save()
            
            assay_plates_created = []
            assay_plate = None
            current_plate_lcps = []
            
            plate_copies = sorted(lcp_by_plate_copy.keys())
            for plate_copy in plate_copies:

                logger.info('plating: %r', plate_copy)
                lcps_to_plate = lcp_by_plate_copy[plate_copy]
                for lcp in lcps_to_plate:

                    if assay_plate is None:
                        assay_plate = create_next_assay_plate()
                        assay_plates_created.append(assay_plate)
                        assay_plate_well_index = 0

                    lcp.cherry_pick_assay_plate = assay_plate
                    well_name = available_assay_plate_wells[assay_plate_well_index]
                    lcp.assay_plate_row = lims_utils.well_name_row_index(well_name)
                    lcp.assay_plate_column = lims_utils.well_name_col_index(well_name)
                    lcp.save()
                    current_plate_lcps.append(lcp)
                    
                    assay_plate_well_index += 1
                    if assay_plate_well_index >= capacity:
                        assay_plate_well_index = 0
                        assay_plate = None
                    
            if cpr.is_randomized_assay_plate_layout:                    
                randomize_lcps(current_plate_lcps)
            
#             for lcp in fulfillable_lcps.all():
#                 
#                 copy_plate = '%s:%s' % (lcp.copy.name, lcp.source_well.plate_number)
#                 
#                 if assay_plate is None:
#                     assay_plate = create_next_assay_plate()
#                     assay_plates_created.append(assay_plate)
#                     assay_plate_well_index = 0
# 
#                 lcp.cherry_pick_assay_plate = assay_plate
#                 well_name = available_assay_plate_wells[assay_plate_well_index]
#                 lcp.assay_plate_row = lims_utils.well_name_row_index(well_name)
#                 lcp.assay_plate_column = lims_utils.well_name_col_index(well_name)
#                 lcp.save()
# 
#                 assay_plate_well_index += 1
#                 if assay_plate_well_index >= capacity:
#                     assay_plate_well_index = 0
#                     assay_plate = None

        # verify that all lcps have been assigned
        for lcp in fulfillable_lcps.all():
            copy_plate = '%s:%s' % (lcp.copy.name, lcp.source_well.plate_number)
            if lcp.cherry_pick_assay_plate is None:
                raise ProgrammingError(
                    'lcp has not been assigned: %s, %r' 
                    % lcp, copy_plate)
        
        # return some stats
        cpap_assignments = [ 
            'Plate ordinal: %d, Picks: %d' 
                % (ap.plate_ordinal, ap.labcherrypick_set.all().count()) 
                for ap in assay_plates_created ]
        
        plate_copies = sorted(lcp_by_plate_copy.keys())
        copy_plate_assigned_msg = [
            (plate_copy, len(lcp_by_plate_copy[plate_copy]))
                for plate_copy in plate_copies] 
        _meta = {
            API_MSG_LCP_PLATES_ASSIGNED: copy_plate_assigned_msg,
            API_MSG_LCP_ASSAY_PLATES_CREATED: cpap_assignments
        }
        _meta.update(copywell_reservation_meta)
        if unfulfillable_wells:
            _meta[API_MSG_LCPS_VOLUME_OVERRIDDEN] = unfulfillable_wells

        # log
        parent_log.diffs.update({
            'number_plates': [previous_number_of_plates, len(assay_plates_created)]
            })
        parent_log.json_field = _meta
        parent_log.save()
        
        return self.build_response(
            request, { API_RESULT_META: _meta }, 
            response_class=HttpResponse, **kwargs)
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def dispatch_delete_lab_cherry_picks(self, request, **kwargs):
        
        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')
        convert_post_to_put(request)

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        cpr_id = param_hash['cherry_pick_request_id']
        logger.info('dispatch_delete_lab_cherry_picks: %r', cpr_id) 
        schema = super(CherryPickRequestResource, self).build_schema(
            user=request.user)
         
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        lcp_query = cpr.lab_cherry_picks.all() 
        if lcp_query.exists():
            allocated_lcp_query = lcp_query.filter(cherry_pick_assay_plate__isnull=False)
            if allocated_lcp_query.exists():
                plated_assay_plates_query = \
                    cpr.cherry_pick_assay_plates.filter(
                        plating_date__isnull=False)                
                if plated_assay_plates_query.exists():
                    raise ValidationError({
                        API_MSG_NOT_ALLOWED: 
                            ('Lab Cherry Picks may not be deleted after plates are plated'),
                        API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count(),
                        API_MSG_LCP_ASSAY_PLATES_PLATED: plated_assay_plates_query.count()
                })
                raise ValidationError({
                    API_MSG_NOT_ALLOWED: 
                        ('Lab cherry pick plates already assigned; '
                        'cancel reservation to deallocate plates'),
                    API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count()
                })

            original_cpr = self._get_detail_response_internal(**{
                'cherry_pick_request_id': cpr_id,
                'includes': '-screener_cherry_picks'
             })
            logger.debug('original_cpr: %r', original_cpr)
            parent_log = self.make_log(request)
            parent_log.key = str(cpr_id)
            parent_log.uri = '/'.join([
                'screen',cpr.screen.facility_id,
                parent_log.ref_resource_name,parent_log.key])        
            parent_log.comment = API_MSG_LCPS_REMOVED
            
            meta = {}
            meta[API_MSG_LCPS_REMOVED] = cpr.lab_cherry_picks.all().count()
            logger.info('delete lcps...')
            LabCherryPick.objects.filter(cherry_pick_request=cpr).delete()
            logger.info('lcp delete done')
            new_cpr = self._get_detail_response_internal(**{
                'cherry_pick_request_id': cpr_id,
                'includes': '-screener_cherry_picks'
            })
            logger.debug('new_cpr: %r', new_cpr)
            parent_log = self.log_patch(
                request, original_cpr, new_cpr, parent_log, 
                excludes=['screener_cherry_picks'])
            parent_log.save()
            
            return self.build_response(
                request,  {API_RESULT_META: meta }, 
                response_class=HttpResponse, **kwargs)
        else:
            logger.info('no LCPs found to delete')
        # Empty response if no action taken
        return self.build_response(
            request,  {API_RESULT_META: 'no LCPs found to delete' }, 
            response_class=HttpResponse, **kwargs)

    @write_authorization
    @un_cache
    @transaction.atomic
    def dispatch_cancel_reservation(self, request, **kwargs):
        
        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')
        convert_post_to_put(request)

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        cpr_id = param_hash['cherry_pick_request_id']
        logger.info('dispatch_cancel_reservation: %r', cpr_id) 
        schema = super(CherryPickRequestResource, self).build_schema(
            user=request.user)
         
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        lcp_query = cpr.lab_cherry_picks.all() 
        if lcp_query.exists():
            allocated_lcp_query = lcp_query.filter(cherry_pick_assay_plate__isnull=False)
            if allocated_lcp_query.exists():
           
                original_cpr = self._get_detail_response_internal(**{
                    'cherry_pick_request_id': cpr_id })
                parent_log = self.make_log(request)
                parent_log.key = str(cpr_id)
                parent_log.uri = '/'.join([
                    'screen',cpr.screen.facility_id,
                    parent_log.ref_resource_name,parent_log.key])        
        
                meta = self._cancel_reservation(cpr, parent_log)

                new_cpr = self._get_detail_response_internal(**{
                    'cherry_pick_request_id': cpr_id })
                parent_log = self.log_patch(
                    request, original_cpr, new_cpr, parent_log, 
                    excludes=['screener_cherry_picks'])
                parent_log.save()
                
                return self.build_response(
                    request,  {API_RESULT_META: meta }, 
                    response_class=HttpResponse, **kwargs)
            else: 
                logger.warn('no allocated LCPs found to delete for CPR: %r', cpr_id)
        else: 
            logger.warn('no LCPs found to delete for CPR: %r', cpr_id)
        # Empty response if no action taken
        return self.build_response(
            request,  {API_RESULT_META: 'no LCPs found to cancel' }, 
            response_class=HttpResponse, **kwargs)
    
    def _cancel_reservation(self, cpr, parent_log):
        logger.info(
            '_cancel_reservation for: %r...', cpr)
        
        cpr_id = cpr.cherry_pick_request_id 
        meta = {}    
        
        if not cpr.screener_cherry_picks.filter(selected=True).exists():
            logger.warn('No screener cherry picks found for %r', cpr)

        lcp_query = cpr.lab_cherry_picks.all() 
        if lcp_query.exists():
            allocated_lcp_query = lcp_query.filter(cherry_pick_assay_plate__isnull=False)
            if allocated_lcp_query.exists():
                plated_assay_plates_query = \
                    cpr.cherry_pick_assay_plates.filter(plating_date__isnull=False)
                if plated_assay_plates_query.exists():
                    raise ValidationError({
                        API_MSG_NOT_ALLOWED: API_MSG_CPR_PLATED_CANCEL_DISALLOWED, 
                        API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count(),
                        API_MSG_LCP_ASSAY_PLATES_PLATED: plated_assay_plates_query.count()
                    })
                lab_cherry_picks_to_deallocate = (
                    cpr.lab_cherry_picks
                        .filter(cherry_pick_assay_plate__isnull=False)
                        .order_by('copy__name','source_well__plate_number') )
                logger.info('lcps to deallocate: %r', 
                    [str(lcp.source_well) for lcp in lab_cherry_picks_to_deallocate])
                parent_log.comment = API_MSG_PLATING_CANCELED
                parent_log.save()
                result_meta = \
                    self.get_copywell_resource().deallocate_cherry_pick_volumes(
                        cpr, lab_cherry_picks_to_deallocate, parent_log)
                meta.update(result_meta)
                
                cpr.date_volume_reserved = None
                cpr.save()

                logger.info('removing previous cherry pick assay plates for: %r, %d', 
                    cpr, cpr.cherry_pick_assay_plates.count())
                meta[API_MSG_CPR_ASSAY_PLATES_REMOVED] = \
                    cpr.cherry_pick_assay_plates.count()
                cpr.cherry_pick_assay_plates.all().delete()
                cpr.lab_cherry_picks.all().update(
                    assay_plate_row=None, assay_plate_column=None)
            else:
                logger.warn('No assay plates to cancel for %r', cpr)
        else:
            logger.warn('No lab cherry picks to cancel for %r', cpr)
        
        return meta
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def dispatch_set_duplex_lab_cherry_picks(self, request, **kwargs):
        '''
        Create Lab Cherry Picks for the Screener Cherry Picks that have already
        been created for the CPR:
        - in this case, select the duplex wells corresponding to the pool wells
        that have been selected as Screener Cherry Picks.
        - see dispatch_set_lab_cherry_picks for the automatic copy selection 
        process.
        '''
        
        logger.info('dispatch_set_duplex_lab_cherry_picks')
        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')
        convert_post_to_put(request)

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
         
        cpr_id = param_hash['cherry_pick_request_id']
        logger.info(
            'dispatch_set_duplex_lab_cherry_picks for: %r...', cpr_id)
         
        schema = super(CherryPickRequestResource, self).build_schema(
            user=request.user)
         
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        if not cpr.screener_cherry_picks.filter(selected=True).exists():
            raise ValidationError(
                key='screener_cherry_picks', 
                msg='No Screener Cherry Picks have been submitted')
        
        original_cpr = self._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id,
            'includes': '-screener_cherry_picks'
         })
        logger.debug('original_cpr: %r', original_cpr)
        
        meta = {}    
        
        parent_log = self.make_log(request)
        parent_log.key = str(cpr_id)
        parent_log.uri = '/'.join([
            'screen',cpr.screen.facility_id,
            parent_log.ref_resource_name,parent_log.key])        
        
        lcp_query = cpr.lab_cherry_picks.all() 
        if lcp_query.exists():
            if cpr.lab_cherry_picks.exists():
                raise ValidationError({
                    API_MSG_NOT_ALLOWED: 
                        ('Lab cherry picks already assigned: (%d); '
                        'delete lab cherry picks to change '
                        'screener cherry pick selections' 
                        % cpr.lab_cherry_picks.count()),
                    API_MSG_LCPS_MUST_BE_DELETED: cpr.lab_cherry_picks.count()
                })
        
        logger.info('create lcps (duplexes)...')
        # TODO: use bulk create to speed this up    
        lcps_created = []
        for scp in cpr.screener_cherry_picks.filter(selected=True):
            
            # Find duplex wells
            if ( not scp.screened_well.reagents.exists() 
                 or not hasattr(scp.screened_well.reagents.all()[0],'silencingreagent')):
                raise ValidationError(
                    key='screened_well',
                    msg='no silencing reagents found for screened well: %s' 
                        % scp.screened_well.well_id )
            sr = scp.screened_well.reagents.all()[0].silencingreagent
            if not sr.duplex_wells.exists():
                raise ValidationError(
                    key='screened_well',
                    msg='no duplex wells found for screened well: %s' 
                        % scp.screened_well.well_id )
            for duplex_well in sr.duplex_wells.all():
                logger.debug('creating lcp for scp: %r, %r', 
                    scp.screened_well.well_id, duplex_well.well_id)
                lab_cherry_pick = LabCherryPick.objects.create(
                    cherry_pick_request=cpr,
                    source_well=duplex_well,
                    screener_cherry_pick=scp)
                lcps_created.append(lab_cherry_pick)
            if len(lcps_created) %100 == 0:
                logger.info('created %d duplex lcps', len(lcps_created))
        logger.info('created %d duplex lcps', len(lcps_created))
        self.get_labcherrypick_resource().clear_cache()
        meta[API_MSG_LCPS_CREATED] = len(lcps_created)

        logger.info('Find and assign copies that are fulfillable...')
        meta_action = self._find_copies(cpr)
        if meta_action:
            meta.update(meta_action)
        
        new_cpr = self._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id,
            'includes': '-screener_cherry_picks'
         })
        logger.debug('new_cpr: %r', new_cpr)

        self.log_patch(request, original_cpr, new_cpr, parent_log, 
            excludes=['screener_cherry_picks'])
        parent_log.json_field = meta
        parent_log.save()
        logger.info('set_duplex_lab_cherry_picks: log: %r, %r', parent_log, parent_log.diffs)
        
        return self.build_response(
            request,  {API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def dispatch_set_lab_cherry_picks(self, request, **kwargs):
        '''
        Create Lab Cherry Picks for the Screener Cherry Picks that have already
        been created for the CPR.
        - find copies for the LCPs:
            copy.usage_type = cherry_pick_source_plates
            plate.status = available
        - assign appropriate copies for the LCPs, if no appropriate copy can be
        found, leave the LCP unfulfilled (see _find_copies).
        '''
        logger.info('dispatch_set_lab_cherry_picks')
        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')
        convert_post_to_put(request)

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
         
        cpr_id = param_hash['cherry_pick_request_id']
        logger.info(
            'dispatch_set_lab_cherry_picks for: %r...', cpr_id)
         
        schema = super(CherryPickRequestResource, self).build_schema(
            user=request.user)
         
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        if not cpr.screener_cherry_picks.filter(selected=True).exists():
            raise ValidationError(
                key='screener_cherry_picks', 
                msg='No Screener Cherry Picks have been submitted')
        
        original_cpr = self._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        
        meta = {}    
        
        parent_log = self.make_log(request)
        parent_log.key = str(cpr_id)
        parent_log.uri = '/'.join([
            'screen',cpr.screen.facility_id,
            parent_log.ref_resource_name,parent_log.key])        
        
        lcp_query = cpr.lab_cherry_picks.all() 
        if lcp_query.exists():
            if cpr.lab_cherry_picks.exists():
                raise ValidationError({
                    API_MSG_NOT_ALLOWED: 
                        ('Lab cherry picks already assigned: (%d); '
                        'delete lab cherry picks to change '
                        'screener cherry pick selections' 
                        % cpr.lab_cherry_picks.count()),
                    API_MSG_LCPS_MUST_BE_DELETED: cpr.lab_cherry_picks.count()
                })
            
        logger.info('create lcps...')    
        lcps_created = []
        for scp in cpr.screener_cherry_picks.filter(selected=True):
           lab_cherry_pick = LabCherryPick.objects.create(
                cherry_pick_request=cpr,
                source_well=scp.screened_well,
                screener_cherry_pick=scp)
           lcps_created.append(lab_cherry_pick)
        self.get_labcherrypick_resource().clear_cache()
        meta[API_MSG_LCPS_CREATED] = len(lcps_created)

        logger.info('Find and assign copies that are fulfillable...')
        meta_action = self._find_copies(cpr)
        if meta_action:
            meta.update(meta_action)
        
        new_cpr = self._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        self.log_patch(request, original_cpr, new_cpr, parent_log, 
            excludes=['screener_cherry_picks'])
        parent_log.json_field = meta
        parent_log.save()
        logger.info('set_lab_cherry_picks: log: %r, %r', parent_log, parent_log.diffs)
        
        return self.build_response(
            request,  {API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)

    def _find_copies(self,cpr):
        '''
        Find and assign copies that are fulfillable for the Lab Cherry Picks
        on the Cherry Pick Request:
        - well copy
        @param cpr Cherry Pick Request with Lab Cherry Picks already created
        '''
        logger.info('find copies for cpr: %d lab cherry picks...', 
            cpr.cherry_pick_request_id)
        #select the "best" copy for each contiguous plate   
        
        if not cpr.lab_cherry_picks.all().exists():
            raise ValidationError(key='lab_cherry_picks', msg='must be set')
        
        self.clear_cache()
        logger.info('find eligible lab_cherry_pick copy wells...')
        eligible_lab_cherry_pick_copywells = \
            self.get_labcherrypick_resource()._get_list_response_internal(
                **{
                    'cherry_pick_request_id': cpr.cherry_pick_request_id,
                    API_PARAM_SHOW_COPY_WELLS: True,
                    'source_copy_well_volume__gte': cpr.transfer_volume_per_well_approved, 
                    'includes': [
                        'source_plate_type','destination_plate_type','source_well_id'
                        'source_copywell_id','-structure_image','-molfile',
                        '-library_plate_comment_array'],
                })
            
        # NOTE: do not consider the retired cherry_pick_source_plates for
        # the server-generated copy assignments
        eligible_lab_cherry_pick_copywells = [
            lcp for lcp in eligible_lab_cherry_pick_copywells
            if not (lcp['source_plate_type']=='cherry_pick_source_plates'
                    and lcp['source_plate_status']=='retired' )]   
        
        logger.info('found %d eligible copy-wells for %d lab cherry picks', 
            len(eligible_lab_cherry_pick_copywells), 
            cpr.lab_cherry_picks.all().count())
        copy_wells_well_set = set([lcp['source_well_id'] for lcp in eligible_lab_cherry_pick_copywells])
        logger.info('lcp candidates found for source wells: %d', len(copy_wells_well_set))
        logger.debug('found eligible: %r', eligible_lab_cherry_pick_copywells)
        
        logger.info('Pick the best copy...')
        copy_sets_by_library = {}
        copy_instance_cache = {}
        pick_candidates_by_library = {}

        for pick_copy in eligible_lab_cherry_pick_copywells:
            
            source_well_id = pick_copy['source_well_id']
            copy_name = pick_copy['source_copy_name']
            library_short_name = pick_copy['library_short_name']
            copy_id = pick_copy['source_copy_id']
            copy_full_name = '%s:%s' % (library_short_name,copy_name)
            copy_instance_cache[copy_full_name] = Copy.objects.get(copy_id=copy_id)
            
            library_copy_set = copy_sets_by_library.get(library_short_name, set())
            library_copy_set.add(copy_full_name)
            copy_sets_by_library[library_short_name] = library_copy_set
            
            library_picks = pick_candidates_by_library.get(library_short_name,{})
            well_picks = library_picks.get(source_well_id, set())
            well_picks.add(copy_full_name)
            library_picks[source_well_id] = well_picks
            pick_candidates_by_library[library_short_name] = library_picks

        lcp_assigned_count = 0
        for library_short_name in pick_candidates_by_library.keys():
            
            copy_set = copy_sets_by_library.get(library_short_name,set())
            well_picks_for_library = pick_candidates_by_library[library_short_name]
            
            minimal_copy_set = lims_utils.find_minimal_satisfying_set(
                copy_set, well_picks_for_library.values())
            logger.info('library: %r, chosen minimal copy set: %r', 
                library_short_name, minimal_copy_set)
            if minimal_copy_set:
                # FIXME: only iterate the lcp's for the library
                for lcp in cpr.lab_cherry_picks.all():
                    if lcp.source_well_id in well_picks_for_library:
                        pick_copy_set = well_picks_for_library.get(lcp.source_well_id, None)
                        if pick_copy_set is None:
                            logger.info('no pick copy set for well: %r', lcp.source_well_id)
                            continue
                        eligible_copies = set(pick_copy_set) & set(minimal_copy_set)
                        if eligible_copies:
                            copy_full_name = sorted(eligible_copies)[0]
                            lcp.copy = copy_instance_cache[copy_full_name]
                            lcp.save()
                            lcp_assigned_count += 1
                        else:
                            logger.info('no eligible copies for %r, in %r',
                                lcp.source_well_id, pick_copy_set)
                    else:
                        logger.debug('no pick copy set for well %r', lcp.source_well_id)
            else:
                logger.info('no minimal copy sets found for library: %r', library_short_name)                
            
        meta = {}
        if lcp_assigned_count == 0:
            meta[API_MSG_WARNING] = 'No eligible copies found'
        else:
            total_count = cpr.lab_cherry_picks.all().count()
            meta = {
                API_MSG_LCPS_CREATED: total_count,
                API_MSG_LCPS_ASSIGNED: lcp_assigned_count,
                API_MSG_LCPS_UNFULFILLED: (total_count-lcp_assigned_count)
            }
        return meta


class ScreenerCherryPickResource(DbApiResource):        

    class Meta:

        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'screenercherrypick'
        authorization = CherryPickRequestAuthorization(resource_name)
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):

        super(ScreenerCherryPickResource, self).__init__(**kwargs)
        self.sr_resource = None
        self.smr_resource = None
        self.reagent_resource = None
        self.cpr_resource = None
        
    def get_sr_resource(self):
        if self.sr_resource is None:
            self.sr_resource = SilencingReagentResource()
        return self.sr_resource
    
    def get_smr_resource(self):
        if not self.smr_resource:
            self.smr_resource = SmallMoleculeReagentResource()
        return self.smr_resource
        
    def get_reagent_resource(self):
        if self.reagent_resource is None:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource
    
    def get_cpr_resource(self):
        if self.cpr_resource is None:
            self.cpr_resource = CherryPickRequestResource()
        return self.cpr_resource
    
    def clear_cache(self):
        DbApiResource.clear_cache(self)
        self.get_cpr_resource().clear_cache()
    
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/(?P<source_well_id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    @write_authorization
    @un_cache
    @transaction.atomic
    def post_list(self, request, schema=None, **kwargs):
        raise NotImplementedError()
    

    @write_authorization
    @transaction.atomic
    @un_cache        
    def patch_list(self, request, **kwargs):
        '''
        Allow patch list for changing selections only 
         - no new picks can be created through patch, picks can only be
         created on the "patch_detail" with the 
         cherry_pick_request.screener_cherry_picks attribute
        '''
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        if 'cherry_pick_request_id' not in kwargs:
            raise BadRequest('cherry_pick_request_id')
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=kwargs['cherry_pick_request_id'])
        logger.info(
            'patch_list: cpr: %r, screen: %r...', cpr, cpr.screen.facility_id)
         
        if cpr.lab_cherry_picks.exists():
            raise ValidationError({
                API_MSG_NOT_ALLOWED: 
                    ('Lab cherry picks already assigned: (%d); '
                    'delete lab cherry picks to change '
                    'screener cherry pick selections' 
                    % cpr.lab_cherry_picks.count()),
                API_MSG_LCPS_MUST_BE_DELETED: cpr.lab_cherry_picks.count()
            })
            
        deserialized = self.deserialize(request)
        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
 
        id_attribute = schema['id_attribute']

        original_cpr_data = self.get_cpr_resource()._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr.cherry_pick_request_id })
        
        original_data = self._get_list_response_internal(**{
            'cherry_pick_request_id': kwargs['cherry_pick_request_id'],
            API_PARAM_SHOW_OTHER_REAGENTS: True })
        if not original_data: 
            raise BadRequest(
                'Can not set Screener Cherry Picks using patch; '
                'must use the cherry pick request selection patch instead')
        
        original_selections = { scp_data['screened_well_id']: scp_data
            for scp_data in original_data }
        current_selections = {}
        for scp_data in original_data:
            screened_well_id = scp_data['screened_well_id']
            searched_well_id = scp_data['searched_well_id']
            if scp_data['selected'] is True:
                if searched_well_id in current_selections:
                    raise ProgrammingError(
                        'Original ScreenerCherryPicks have a duplicate:'
                        'searched_well: %r, has >1 selection: %r'
                        % (searched_well_id, original_data))
                current_selections[searched_well_id] = screened_well_id
        # Validate selections are valid and allowed:
        required_for_patch = set([
            'screened_well_id', 'searched_well_id', 'selected'])
        selection_updates = {}
        for selection_update in deserialized:
            if not required_for_patch.issubset(set(selection_update.keys())):
                raise ValidationError(
                    key='required_fields',
                    msg=(
                        'record: %r, missing fields: %r'
                        % (selection_update,
                        required_for_patch-set(selection_update.keys())))
                    )
            candidate_screened_well_id = selection_update['screened_well_id']
            if not candidate_screened_well_id in original_selections:
                raise ValidationError(
                    key='screened_well_id',
                    msg='not found in the current/alternate screener_cherry_picks: %r' 
                        % candidate_screened_well_id )
            selection_update['selected'] = parse_val(
                selection_update['selected'],'selected', 'boolean')
            selection_updates[selection_update['screened_well_id']]=selection_update
        # validate only one selection per searched_well
        for screened_well_id,selection_update in selection_updates.items():
            if selection_update['selected'] is True:
                searched_well_id = selection_update['searched_well_id']
                
                for other_well_id, su2 in selection_updates.items():
                    if other_well_id != screened_well_id:
                        if su2['searched_well_id'] == searched_well_id:
                            if su2['selected'] is True:
                                raise ValidationError(
                                    key='selected',
                                    msg=(
                                        'only one well can be selected for '
                                        'the searched_well_id: %r' % searched_well_id))
                if ( searched_well_id in current_selections
                    and searched_well_id not in selection_updates):
                    selection_updates[searched_well_id] = {
                        'selected': False,
                        'searched_well_id': searched_well_id,
                        'screened_well_id': screened_well_id }
        # make sure others in group are deselected. 
        # Fixme: tidy up this logic
        for scp_data in original_data:
            if scp_data['selected'] is True:
                for selection_update in selection_updates.values():
                    if selection_update['searched_well_id']==scp_data['searched_well_id']:
                        if selection_update['selected'] is True:
                            scp_data['selected'] = False
                            selection_updates[scp_data['screened_well_id']] = \
                                scp_data
                            break
        
        if not selection_updates:
            return self.build_response(
                request,  {API_RESULT_META: 'no new Selections found' }, 
                response_class=HttpResponse, **kwargs)
        
        logger.info('selection updates: %r', 
            [(scp['screened_well_id'],scp['selected'],scp['searched_well_id']) 
                for scp in selection_updates.values()])
                    
        original_scps = { 
            scp.screened_well_id: scp 
                for scp in cpr.screener_cherry_picks.all() }
        
        logger.debug('original_scps: %r', original_scps.keys())
        
        messages = []
        scps_to_create = []
        scps_to_reselect = []
        scps_to_unselect = []
        for screened_well_id, scp_selection_update in selection_updates.items():
            selected = scp_selection_update['selected']
            if selected is True:
                if screened_well_id in original_scps:
                    if original_scps[screened_well_id].selected is True:
                        messages.append(
                            'screened_well_id: %r is already selected'  % screened_well_id)
                    else:
                        scps_to_reselect.append(screened_well_id)
                else:
                    scps_to_create.append(screened_well_id)
            else:
                if screened_well_id in original_scps:
                    if original_scps[screened_well_id].selected is True:
                        scps_to_unselect.append(screened_well_id)
                    else:
                        messages.append(
                            'not currently selected: %r' % screened_well_id)
                else:
                    messages.append(
                        'not unselecting, well not found in original_scps: %r'
                        % screened_well_id)
        
        if not scps_to_create and not scps_to_reselect and not scps_to_unselect:
            _data = { API_RESULT_META: messages }
            return self.build_response(
                request, _data, response_class=HttpResponse, **kwargs)
        
        if scps_to_create:
            cherry_pick_wells = \
                CherryPickRequestResource.find_wells(scps_to_create)
                
            for well in cherry_pick_wells:
                full_selection = original_selections[well.well_id]
                searched_well_id = full_selection['searched_well_id']
                screener_cherry_pick = ScreenerCherryPick.objects.create(
                    cherry_pick_request=cpr,
                    screened_well=well,
                    searched_well_id=searched_well_id,
                    selected=True)
                
                # TODO: if selected well is restricted library, 'requires permission'
                # add user warning message
        if scps_to_unselect:
            for screened_well_id in scps_to_unselect:
                original_scps[screened_well_id].selected = False
                original_scps[screened_well_id].save()
        if scps_to_reselect:
            for screened_well_id in scps_to_reselect:
                original_scps[screened_well_id].selected = True
                original_scps[screened_well_id].save()
        
        result_message = {
            API_MSG_SCP_CREATED: ', '.join(scps_to_create),
            API_MSG_SCP_UNSELECTED: ', '.join(scps_to_unselect),
            API_MSG_SCP_RESELECTED: ', '.join(scps_to_reselect)
        }
        if messages:
            result_message[API_MSG_COMMENTS] = messages
        
        self.get_cpr_resource().clear_cache()
        new_cpr_data = self.get_cpr_resource()._get_detail_response_internal(
            **{ 'cherry_pick_request_id': cpr.cherry_pick_request_id })
        logger.debug('old cpr: %r, new cpr: %r', original_cpr_data, new_cpr_data)
        log = self.get_cpr_resource().make_log(request)
        log.key = str(cpr.cherry_pick_request_id)
        log.uri = '/'.join([
            'screen',cpr.screen.facility_id,log.ref_resource_name,log.key])        
        log.json_field = result_message
        self.get_cpr_resource().log_patch(
            request, original_cpr_data,new_cpr_data, 
            excludes=['screener_cherry_picks'])
        logger.info('cpr_log: %r, %r',log, log.diffs)
        # Note: because the actual screener cherry picks are not being logged, 
        # log the count, even if it has not changed
        log.diffs['screener_cherry_pick_count'] = [
            original_cpr_data['screener_cherry_pick_count'],
            new_cpr_data['screener_cherry_pick_count']]
        log.save()
        
        _data = { API_RESULT_META: result_message }
        return self.build_response(
            request, _data, response_class=HttpResponse, **kwargs)


    @read_authorization
    def get_detail(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def get_schema(self, request, **kwargs):
        if not 'cherry_pick_request_id' in kwargs:
            return self.build_response(request, 
                self.build_schema(user=request.user),**kwargs)
        
        cherry_pick_request_id = kwargs.pop('cherry_pick_request_id')
        try:
            cpr = CherryPickRequest.objects.get(
                cherry_pick_request_id=cherry_pick_request_id)
            # NOTE: does not support disinction between Small Molecule and Natural Product
            return self.build_response(
                request, 
                self.build_schema(
                    user=request.user, 
                    library_classification=cpr.screen.screen_type),
                **kwargs)
            
        except ObjectDoesNotExist, e:
            raise Http404(
                'Can not build schema - CherryPickRequest ID needed'
                'no cpr found for cherry_pick_request_id: %r' % cherry_pick_request_id)

    def build_schema(self, library_classification=None, user=None):
        logger.info('build sreenercherrypick schema for library_classification: %r',
            library_classification)
        schema = deepcopy(
            super(ScreenerCherryPickResource, self).build_schema(user=user))
        original_fields = schema['fields']
        if library_classification:
            # Add in reagent fields
            sub_data = self.get_reagent_resource().build_schema(
                library_classification=library_classification, user=user)
            
            newfields = {}
            newfields.update(sub_data['fields'])
            schema['fields'] = newfields
        
        # Add in well fields    
        well_schema = WellResource().build_schema(user=user)
        schema['fields'].update(well_schema['fields'])
        
        # Turn off the visibility of all inherited fields
        for key,field in schema['fields'].items():
            if 'l' in field['visibility']:
                field['visibility'].remove('l')
            if 'd' in field['visibility']:
                field['visibility'].remove('d')
        
        # Overlay the original scp fields on the top
        schema['fields'].update(original_fields)
        logger.info('schema  fields: %r', schema['fields'].keys())
        return schema

    def build_sqlalchemy_columns(self, fields, base_query_tables=None, 
            custom_columns=None):
        sub_columns = self.get_sr_resource().build_sqlalchemy_columns(fields)
        sub_columns.update(self.get_smr_resource().build_sqlalchemy_columns(fields))
#         sub_columns['plate_number'] = (literal_column(
#             "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
#             .label('plate_number'))
        if custom_columns is not None:
            sub_columns.update(custom_columns)
        return DbApiResource.build_sqlalchemy_columns(self, 
            fields, base_query_tables=base_query_tables, 
            custom_columns=sub_columns)

    def build_list_response(self, request, schema=None, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False

        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
        extra_params = {}
        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id is None:
            raise BadRequest('cherry_pick_request_id is required')
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cherry_pick_request_id)
        extra_params['CPR']=cherry_pick_request_id
        source_well_id = param_hash.pop('source_well_id', None)
        if source_well_id:
            param_hash['source_well_id__eq'] = source_well_id

        show_other_reagents = parse_val(
            param_hash.get(API_PARAM_SHOW_OTHER_REAGENTS, None),
            API_PARAM_SHOW_OTHER_REAGENTS, 'boolean')
        if show_other_reagents is True:
            extra_params[API_PARAM_SHOW_OTHER_REAGENTS] = None

        show_alternates = parse_val(
            param_hash.get(API_PARAM_SHOW_ALTERNATE_SELECTIONS, None),
            API_PARAM_SHOW_ALTERNATE_SELECTIONS, 'boolean')
        if show_alternates is True:
            extra_params[API_PARAM_SHOW_ALTERNATE_SELECTIONS] = None
        try:
            
            # Note: build schema for each request to use the subtype
            schema = self.build_schema(library_classification=cpr.screen.screen_type)
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            manual_field_includes.add('searched_well_id')
            manual_field_includes.add('selected')
            manual_field_includes.add('cherry_pick_request_id')
            
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema, **extra_params)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                                  
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            # specific setup 
            base_query_tables = [
                'screen',
                'cherry_pick_request',
                'well',
                'reagent',
                'library'
                ]
            # build the query statement
            _cpr = self.bridge['cherry_pick_request']
            _scp = self.bridge['screener_cherry_pick']
            _lcp = self.bridge['lab_cherry_pick']
            _screen = self.bridge['screen']
            _well = self.bridge['well']
            _reagent = self.bridge['reagent']
            _library = self.bridge['library']
            
            if show_other_reagents is True or show_alternates is True:

                _original_scps = (
                    select([
                        _scp.c.screener_cherry_pick_id,
                        _scp.c.screened_well_id,
                        _scp.c.searched_well_id,
                        _scp.c.selected,
                        _scp.c.cherry_pick_request_id,
                        _reagent.c.vendor_identifier, 
                        _reagent.c.vendor_name])
                    .select_from(
                        _scp.join(_reagent, _scp.c.screened_well_id==_reagent.c.well_id)
                            .join(_well,_scp.c.screened_well_id==_well.c.well_id)
                            .join(_library, _well.c.library_id==_library.c.library_id))
                    .where(_scp.c.cherry_pick_request_id==cherry_pick_request_id)
                    .where(_scp.c.searched_well_id==_scp.c.screened_well_id)
                    )
                _original_scps = _original_scps.cte('original_scps')

                _alternates = (
                    select([_reagent.c.well_id,_original_scps.c.searched_well_id])
                    .select_from(
                        _reagent.join(
                            _original_scps, 
                            and_(_reagent.c.vendor_identifier==_original_scps.c.vendor_identifier,
                                 _reagent.c.vendor_name==_original_scps.c.vendor_name)))
                    .where(_reagent.c.well_id!=_original_scps.c.searched_well_id)
                    ).cte('alternate_scps')

                combined_scps = union(
                    select([
                        _scp.c.screener_cherry_pick_id,
                        _scp.c.searched_well_id,
                        _scp.c.screened_well_id,
                        _scp.c.selected,
                        _scp.c.cherry_pick_request_id,
                        _reagent.c.vendor_identifier,
                        _reagent.c.vendor_name])
                    .select_from(
                        _scp.join(
                            _reagent,_scp.c.screened_well_id==_reagent.c.well_id))
                    .where(_scp.c.cherry_pick_request_id==cherry_pick_request_id),
                    select([
                        literal_column("0"),
                        _alternates.c.searched_well_id,
                        _alternates.c.well_id,
                        literal_column('false').label('selected'),
                        literal_column(cherry_pick_request_id).label('cherry_pick_request_id'),
                        _reagent.c.vendor_identifier,
                        _reagent.c.vendor_name
                        ])
                    .select_from(
                        _alternates.join(
                            _reagent,_alternates.c.well_id==_reagent.c.well_id))
                    .where(not_(_alternates.c.well_id.in_(
                        select([_scp.c.screened_well_id])
                        .select_from(_scp)
                        .where(_scp.c.cherry_pick_request_id==cherry_pick_request_id))))
                    )
                combined_scps = combined_scps.cte('combined')
                
                working_scp = combined_scps
            else:
                working_scp = _scp
                
            custom_columns = {
                'screener_cherry_pick_id': working_scp.c.screener_cherry_pick_id,
                'screened_well_id': working_scp.c.screened_well_id,
                'selected': working_scp.c.selected,
                'cherry_pick_request_id': working_scp.c.cherry_pick_request_id,
            }
            custom_columns['searched_well_id'] = working_scp.c.searched_well_id
                
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            
            j = join(working_scp,_cpr, 
                working_scp.c.cherry_pick_request_id == _cpr.c.cherry_pick_request_id)
            j = j.join(_screen,
                _cpr.c.screen_id == _screen.c.screen_id)
            j = j.join(_well, _well.c.well_id==working_scp.c.screened_well_id)
            j = j.join(_library, _well.c.library_id==_library.c.library_id)
            j = j.join(_reagent, working_scp.c.screened_well_id==_reagent.c.well_id)
            
            stmt = select(columns.values()).select_from(j)
            stmt = stmt.where(_cpr.c.cherry_pick_request_id==cherry_pick_request_id)

            if show_alternates is True:
                with get_engine().connect() as conn:
                    _alternates = (
                        select([distinct(_scp.c.searched_well_id)])
                        .select_from(_scp)
                        .where(_scp.c.cherry_pick_request_id==cherry_pick_request_id)
                        .where(_scp.c.searched_well_id!=_scp.c.screened_well_id)
                        .where(_scp.c.selected == True)
                        )
                    alternate_searched_well_ids = set([x[0] for x in 
                        conn.execute(_alternates)])
                    stmt = stmt.where(working_scp.c.searched_well_id.in_(alternate_searched_well_ids))
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)

            if not order_clauses:
                # Ordering for well_id must be alphanumeric
                # For string field ordering, double sort as numeric and text
                order_clause = text(
                    "(substring({field_name}, '^[0-9]+'))::int asc " # cast to integer
                    ",substring({field_name}, ':(.*$)') asc  "  # works as text
                    .format(field_name='searched_well_id'))
                if show_other_reagents is True or show_alternates is True:
                    
                    stmt = stmt.order_by(
                        order_clause,
                        desc(column('searched_well_id')==column('screened_well_id')),
                        desc(column('library_plate')),
                        desc(column('screened_well_id')),
                    )
                else:
                    stmt = stmt.order_by(
                        order_clause,
                        asc(column('library_plate')),
                        asc(column('screened_well_name'))
                        )
            
            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None))
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

class LabCherryPickResource(DbApiResource):        

    LCP_COPYWELL_KEY = '{library_short_name}/{source_copy_name}/{source_well_id}'

    class Meta:

        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'labcherrypick'
        authorization = CherryPickRequestAuthorization(resource_name)
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):

        super(LabCherryPickResource, self).__init__(**kwargs)
        self.sr_resource = None
        self.smr_resource = None
        self.reagent_resource = None
        self.cpr_resource = None
        self.copywell_resource = None
    
    def get_copywell_resource(self):
        if self.copywell_resource is None:
            self.copywell_resource = CopyWellResource()
        return self.copywell_resource
        
    def get_sr_resource(self):
        if self.sr_resource is None:
            self.sr_resource = SilencingReagentResource()
        return self.sr_resource
    
    def get_smr_resource(self):
        if not self.smr_resource:
            self.smr_resource = SmallMoleculeReagentResource()
        return self.smr_resource
        
    def get_reagent_resource(self):
        if self.reagent_resource is None:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource

    def get_cpr_resource(self):
        if self.cpr_resource is None:
            self.cpr_resource = CherryPickRequestResource()
        return self.cpr_resource
    
    def clear_cache(self):
        DbApiResource.clear_cache(self)
        self.get_cpr_resource().clear_cache()
    
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/(?P<source_well_id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
        
    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_list(self, request, **kwargs):
        
        DEBUG_LCP = False or logger.isEnabledFor(logging.DEBUG)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        cw_formatter = LabCherryPickResource.LCP_COPYWELL_KEY
        
        convert_post_to_put(request)
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        param_override = parse_val(
            param_hash.get(API_PARAM_OVERRIDE, False),
                API_PARAM_OVERRIDE, 'boolean')
        set_deselected_to_zero = parse_val(
            param_hash.get(API_PARAM_SET_DESELECTED_TO_ZERO, False),
                API_PARAM_SET_DESELECTED_TO_ZERO, 'boolean')
        id_attribute = schema['id_attribute']

        if 'cherry_pick_request_id' not in kwargs:
            raise BadRequest('cherry_pick_request_id')
        cpr_id = kwargs['cherry_pick_request_id']
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        logger.info(
            'patch_list: cpr: %r, screen: %r...', cpr, cpr.screen.facility_id)
        original_cpr = \
            self.get_cpr_resource()._get_detail_response_internal(
                **{'cherry_pick_request_id': cpr_id })
        # NOTE: not logging all LCP updates at this time, copywell logs can be used
        # original_lab_cherry_pick_copywells = \
        #     self._get_list_response_internal(
        #         **{
        #             'cherry_pick_request_id': cpr_id,
        #             'source_copy_name__is_null': False,
        #             'includes': '*'
        #         })
         
        deserialized = self.deserialize(request)
        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
        
        self.get_cpr_resource().validate_cpr_for_plating(cpr) 

        if not cpr.lab_cherry_picks.all().exists():
            raise ValidationError(
                key='lcp_selection_updates', 
                msg='No Lab Cherry Picks have been created')
        current_lcps = { 
            lcp.source_well_id:lcp for lcp in cpr.lab_cherry_picks.all() }
        parent_log = self.get_cpr_resource().make_log(request)
        parent_log.key = str(cpr.cherry_pick_request_id)
        parent_log.uri = '/'.join([
            'screen',cpr.screen.facility_id,
            parent_log.ref_resource_name,parent_log.key])        
        parent_log.save()

        is_mapped = cpr.lab_cherry_picks.filter(
            cherry_pick_assay_plate__isnull=False).exists()
        if is_mapped is True:
            plated_assay_plates_query = \
                cpr.cherry_pick_assay_plates.filter(
                    plating_date__isnull=False)                
            if plated_assay_plates_query.exists():
                raise ValidationError({
                    API_MSG_NOT_ALLOWED: API_MSG_CPR_PLATED_CANCEL_DISALLOWED, 
                    API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count(),
                    API_MSG_LCP_ASSAY_PLATES_PLATED: plated_assay_plates_query.count()
                })
            if param_override is not True:
                raise ValidationError({
                    API_PARAM_OVERRIDE: 'required',
                    API_MSG_WARNING: 
                        ('Lab cherry pick plates already assigned; '
                        'cancel reservation to deallocate plates, '
                        'or choose override to keep plating and change assignments'),
                    API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count()
                })
                
        # Validate selections are valid and allowed:
        required_for_patch = set([
            'source_well_id', 'selected'])
        copies = {}
        plates = {}
        selection_updates = {}
        selections_per_well = defaultdict(list)
        for selection_update in deserialized:
            if DEBUG_LCP:
                logger.info('selection update: %r', selection_update)
            if not required_for_patch.issubset(set(selection_update.keys())):
                raise ValidationError(
                    key='required_fields',
                    msg=(
                        'record: %r, missing fields: %r'
                        % (selection_update,
                        required_for_patch-set(selection_update.keys())))
                    )
            copy_id_fields = ['source_copy_id','source_copy_name']
            if 'source_copy_id' not in selection_update.keys():
                if 'source_copy_name' not in selection_update.keys():
                    raise ValidationError(
                        key='required_fields',
                        msg=('record: %r, must contain either: %r'
                            % (selection_update,copy_id_fields)))
            if set(copy_id_fields).issubset(set(selection_updates.keys())):
                raise ValidationError(
                    key='copy_id_fields',
                    msg='May only include one of: %r' % copy_id_fields)

            source_well_id = selection_update.get('source_well_id')
            if source_well_id not in current_lcps:
                error_dict = {
                    'source_well_id': 
                        'not a current lab cherry pick: %r' % source_well_id }
                if is_mapped:
                    error_dict[API_MSG_WARNING] = (
                        'lab cherry pick reservation must be canceled in '
                        'order to add new selections')
                raise ValidationError(error_dict)
            # locate copies:
            copy = None
            plate = None
            plate_number = lims_utils.well_id_plate_number(source_well_id)
            source_copy_id = selection_update.get('source_copy_id', None)
            if source_copy_id:
                # TODO: Raise a nice Validation Error
                copy = Copy.objects.get(copy_id=source_copy_id)
                plate = Plate.objects.get(
                    plate_number=plate_number,copy=copy)
            else:
                source_copy_name = selection_update.get('source_copy_name')
                plate = Plate.objects.get(
                    plate_number=plate_number,copy__name=source_copy_name)
                copy = plate.copy
            
            selection_update['copy'] = copy
            selection_update['plate'] = plate # cache the plate for later use
            selection_update['source_copy_name'] = copy.name
            selection_update['library_short_name'] = copy.library.short_name
            selection_copy_name = cw_formatter.format(**selection_update)
            selection_updates[selection_copy_name] = selection_update
            
            if selection_update['selected'] is True:
                selections_per_well[source_well_id].append(copy.name)
                
        logger.info('selection_updates: %r', selection_updates.keys())
        # Check that only one copy is selected per well
        errors = []
        for source_well_id, copies in selections_per_well.items():
            if len(copies) > 1:
                errors.append('%s: %s' %(source_well_id,','.join(copies)))
        if errors:
            logger.info('errors: %r', errors)
            raise ValidationError(
                key = API_MSG_LCP_MULTIPLE_SELECTIONS_SUBMITTED,
                msg = '\n'.join(errors))
        logger.debug('selection_updates: %r', selection_updates.keys())
        
        lcps_to_deselect = set()
        # First, find all of the deselections
        for selection_copy_name, selection_update in selection_updates.items():
            source_well_id = selection_update['source_well_id']
            current_lcp = current_lcps[source_well_id]
            current_copy_name = None
            if current_lcp.copy is not None:
                current_copy_name = cw_formatter.format(
                    library_short_name=current_lcp.copy.library.short_name,
                    source_copy_name=current_lcp.copy.name,
                    source_well_id=source_well_id)
            if selection_update['selected'] is True: 
                if selection_copy_name != current_copy_name:
                    logger.info('lcp to change: %r change %r to %r', 
                        source_well_id, current_copy_name, selection_copy_name)
                    lcps_to_deselect.add(current_lcp)
                else:
                    logger.info('lcp is already selected: %r, %r', 
                        current_copy_name, selection_update)
            else:
                if selection_copy_name == current_copy_name:
                    logger.info('lcp to deselect: %r', current_lcp)
                    lcps_to_deselect.add(current_lcp)
                else:
                    logger.info('lcp is already unselected: %r, %r', 
                        current_copy_name, selection_update)
        
        # Second, if mapped and deselections exist, deallocate
        result_meta_allocate = {}
        if is_mapped is True and len(lcps_to_deselect) > 0:
            logger.info('Already mapped, lcps to deallocate: %r', 
                lcps_to_deselect)
            # Note: signal update_screening_count=False
            # do not adjust the plate.cplt_screening_count or the 
            # copywell.cherry_pick_screening_count (see deallocate_cherry_pick_volumes for
            # rational.
            result_meta_allocate = \
                self.get_copywell_resource().deallocate_cherry_pick_volumes(
                    cpr, lcps_to_deselect, parent_log,
                    set_deselected_to_zero=set_deselected_to_zero,
                    update_screening_count=False)

        # Third, update the LCPS
        changed = []
        deselected = []
        selected = []
        lcps_to_allocate = set()
        for selection_copy_name, selection_update in selection_updates.items():
            source_well_id = selection_update['source_well_id']
            current_lcp = current_lcps[source_well_id]
            current_copy_name = None
            if current_lcp.copy is not None:
                current_copy_name = cw_formatter.format(
                    library_short_name=current_lcp.copy.library.short_name,
                    source_copy_name=current_lcp.copy.name,
                    source_well_id=source_well_id)
            logger.info('current lcp: %r, %r', current_copy_name, current_lcp)
            logger.info('selection_update: %r to %r', 
                selection_copy_name, selection_update['selected'])
            if selection_update['selected'] is True: 
                if selection_copy_name != current_copy_name:
                    if current_copy_name is None:
                        selected.append(selection_copy_name)
                    else:
                        changed.append([
                            current_copy_name, selection_copy_name])
                    current_lcp.copy = selection_update['copy']
                    current_lcp.is_manually_selected = True
                    current_lcp.save()
                    lcps_to_allocate.add(current_lcp)
                else:
                    logger.info('lcp is already selected: %r, %r', 
                        current_copy_name, selection_update)
            else:
                if selection_copy_name == current_copy_name:
                    deselected.append(current_copy_name)
                    current_lcp.copy = None
                    current_lcp.is_manually_selected = False
                    current_lcp.save()
                else:
                    logger.info('lcp is already unselected: %r, %r', 
                        current_copy_name, selection_update)

        # Fourth, if selection updates exist, and already is_mapped, allocate
        # NOTE: volume will be taken without requiring overrides 
        # for insufficient volume (User has already sent an override for 
        # mapping changes)
        if is_mapped is True and len(lcps_to_allocate) > 0:
            logger.info('Already mapped, extra lcps to allocate: %r', 
                lcps_to_deselect)
            # Plate cplt_screening_count should only be updated if this is the 
            # first lcp for the plate. If lcps already exist for the plate,
            # then the cplt_screening_count has already been adjusted.
            new_plate_assignments = \
                set([lcp['plate'] for lcp in selection_updates.values()])
            current_plate_assignments = set()
            for lcp in [lcp for 
                lcp in current_lcps.values() if lcp.copy is not None]:
                copy = lcp.copy
                try:
                    plate = Plate.objects.get(
                        plate_number=lcp.source_well.plate_number, 
                        copy=lcp.copy)
                except ObjectDoesNotExist:
                    logger.exception('Can not find plate: %r, copy: %r', 
                        lcp.source_well.plate_number, lcp.copy)
                    raise
                new_plate_assignments.add(plate)
            plates_to_ignore = current_plate_assignments-new_plate_assignments
            result_meta = \
                self.get_copywell_resource().reserve_cherry_pick_volumes(
                    cpr, lcps_to_allocate, parent_log, 
                    plates_to_ignore=plates_to_ignore)
            result_meta_allocate.update(result_meta)
            
            cpr.date_volume_reserved = _now().date() 
            cpr.save()
        
        changed = sorted(changed)
        deselected = sorted(deselected)
        selected = sorted(selected)
        
        # Check for insufficient well volumes
        unfulfillable_wells = []
        new_lab_cherry_pick_copywells = \
            self._get_list_response_internal(
                **{
                    'cherry_pick_request_id': cpr.cherry_pick_request_id,
                    'source_copy_name__is_null': False,
                    'includes': [
                        'source_plate_type','destination_plate_type',
                        'source_copywell_id','source_copy_well_volume',
                        'volume_approved',
                        '-structure_image','-molfile', 
                        '-library_plate_comment_array'],
                })
        for lcp_cw in new_lab_cherry_pick_copywells:
            name = cw_formatter.format(**lcp_cw)
            if ( Decimal(lcp_cw['source_copy_well_volume'])
                    < cpr.transfer_volume_per_well_approved ):
                logger.info(
                    'vol requires override: lcp_cw: %r, approved: %r, available: %r', 
                    name, cpr.transfer_volume_per_well_approved, 
                    Decimal(lcp_cw['source_copy_well_volume']))
                unfulfillable_wells.append(cw_formatter.format(**lcp_cw))
            if DEBUG_LCP:
                logger.info('checking: %r < %r',
                    Decimal(lcp_cw['source_copy_well_volume']),
                    cpr.transfer_volume_per_well_approved)
        warning_messages = {}
        if unfulfillable_wells:
            warning_messages[API_MSG_LCPS_INSUFFICIENT_VOLUME] = unfulfillable_wells
        if DEBUG_LCP:
            logger.info('warning msg: %r', warning_messages)
        # final tally:
        _meta = {
            API_MSG_LCP_CHANGED: changed,
            API_MSG_LCP_DESELECTED: deselected,
            API_MSG_LCP_SELECTED: selected,
            API_MSG_WARNING: warning_messages,
        }
        if result_meta_allocate:
            _meta.update(result_meta_allocate)
        
        logger.info('result_meta: %r', _meta)
        
        self.get_cpr_resource().clear_cache()
        new_cpr = self.get_cpr_resource()._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        self.get_cpr_resource().log_patch(
            request, original_cpr, new_cpr, log=parent_log)
        # Store the cpr.date_volume_reserved as a marker even if not changed
        if 'date_volume_reserved' not in parent_log.diffs:
            parent_log.diffs['date_volume_reserved'] = [
                cpr.date_volume_reserved,cpr.date_volume_reserved]
        parent_log.json_field = _meta
        parent_log.save()

        # TODO: not logging all lcp changes at this time; copywell child logs 
        # can be used 
        # original_lab_cherry_pick_copywells = {
        #     lcp['source_well_id']:lcp for lcp in original_lab_cherry_pick_copywells}
        # new_lab_cherry_pick_copywells = {
        #     lcp['source_well_id']:lcp for lcp in new_lab_cherry_pick_copywells}

        return self.build_response(
            request, { API_RESULT_META: _meta }, 
            response_class=HttpResponse, **kwargs)

    def get_schema(self, request, **kwargs):
        ''' Generate an HttpResponse for the schema '''
        
        return self.build_response(
            request, self.build_schema(user=request.user,**kwargs))

    def build_schema(self, library_classification=None, user=None, **kwargs):

        if library_classification is None:
            cpr_id = kwargs.get('cherry_pick_request_id', None)
            if cpr_id is not None:
                try:
                    cpr = CherryPickRequest.objects.get(pk=cpr_id)
                    library_classification=cpr.screen.screen_type
                except:
                    logger.exception('Note: building generic lab cherry pick schema'
                        '(no cherry_pick_request_id provided)')
        try:
            schema = deepcopy(
                super(LabCherryPickResource, self).build_schema(user=user))
            original_fields = schema['fields']
            # 20170516 - keep the well_id field; required for the structure_image field
            # omit_redundant_fields = [
            #     'well_id','plate_number','well_name','library_well_type',
            #     'library_short_name','library_name']
            omit_redundant_fields = [
                'plate_number','well_name','library_well_type',
                'library_short_name','library_name']
            if library_classification:
                # Add in reagent fields
                sub_data = self.get_reagent_resource().build_schema(
                    library_classification=library_classification, user=user)
                newfields = {}
                sub_fields = {key:field for key,field 
                    in sub_data['fields'].items() 
                        if key not in omit_redundant_fields}
                newfields.update(sub_fields)
                schema['fields'] = newfields
            
            # Add in well fields    
            well_schema = WellResource().build_schema(user=user)
            sub_fields = {key:field for key,field 
                in well_schema['fields'].items() 
                    if key not in omit_redundant_fields}
            schema['fields'].update(sub_fields)
            
            # Turn off the visibility of all inherited fields
            for key,field in schema['fields'].items():
                if 'l' in field['visibility']:
                    field['visibility'].remove('l')
                if 'd' in field['visibility']:
                    field['visibility'].remove('d')
            
            # Overlay the original lcp fields on the top
            schema['fields'].update(original_fields)
            
            
            logger.debug('new lcp fields: %r',
                [(field['key'],field['scope']) 
                    for field in schema['fields'].values()])
            return schema
        except Exception, e:
            logger.exception('xxx: %r', e)
            raise

    @read_authorization
    def get_lab_cherry_pick_plating_schema(self, request, **kwargs):
        ''' Generate an HttpResponse for the plate_mapping specific schema '''
        
        if not 'cherry_pick_request_id' in kwargs:
            raise BadRequest('cherry_pick_request_id is required')
        
        return self.build_response(
            request, 
            self.build_lab_cherry_pick_plating_schema(request.user, **kwargs),
            **kwargs)

    def build_lab_cherry_pick_plating_schema(self, user, **kwargs):

        # Note: build schema for each request to use the subtype
        cherry_pick_request_id = kwargs.get('cherry_pick_request_id', None)
        if cherry_pick_request_id is None:
            raise BadRequest('cherry_pick_request_id is required')
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cherry_pick_request_id)
        schema = self.build_schema(
            user=user,
            library_classification=cpr.screen.screen_type)
        
        # modify the schema for the plate mapping view
        
        visible_fields = [
            'cherry_pick_plate_number',
            'destination_well',
            'destination_plate_type',
            'library_short_name',
            'source_copy_name',
            'library_plate',
            'source_well_name',
            'source_copy_well_volume',
            'volume_approved',            
            'source_plate_status',
            'source_copy_usage_type',
            'source_plate_type',
            'location',               
           ]
        fields = schema['fields']
        for key,field in fields.items():
            if key in visible_fields:
                field['visibility'] = ['l']
                # NOTE: don't alter ordinal if possible, 
                # this makes it problematic included fields later
                # field['ordinal'] = visible_fields.index(key)
            else:
                if 'l' in field['visibility']:
                    field['visibility'].remove('l')
        fields['cherry_pick_plate_number']['ordinal'] = -10
        fields['destination_well']['ordinal'] = -9
        fields['source_copy_well_volume']['title'] = \
            'Source CopyWell Volume (after transfer)'
        fields['source_copy_well_volume']['description'] = \
            'Source CopyWell Volume (after transfer of '\
            'cherry pick volume to the destination well)'
            
        logger.debug('plate mapping visible fields: %r', 
            {k:{'key': v['key'], 'scope':v['scope'], 'title':v['title']}
                for k,v in fields.items()  })
        return schema
                            
    @read_authorization
    def get_detail(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        # Remove the default, user specific schema; lab cherry pick schema will
        # be specific to the CPR #
        kwargs['schema'] = None
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        # Remove the default, user specific schema; lab cherry pick schema will
        # be specific to the CPR #
        kwargs['schema'] = None
        return self.build_list_response(request, **kwargs)

    def build_sqlalchemy_columns(self, fields, base_query_tables=None, 
            custom_columns=None):
        sub_columns = self.get_sr_resource().build_sqlalchemy_columns(fields)
        sub_columns.update(self.get_smr_resource().build_sqlalchemy_columns(fields))
#         sub_columns['plate_number'] = (literal_column(
#             "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
#             .label('plate_number'))
        if custom_columns is not None:
            sub_columns.update(custom_columns)
        return DbApiResource.build_sqlalchemy_columns(self, 
            fields, base_query_tables=base_query_tables, 
            custom_columns=sub_columns)
        
    def build_list_response(self, request, schema=None, **kwargs):
        '''
        Optimized to show the lcp's for one cherry pick request at a time
        
        Options:
        - Default (no options): show only the LCPs that have been created. Do 
            not show the copy assignments.
        - API_PARAM_SHOW_COPY_WELLS:
            show copy wells for the LCPs where: 
            copy type is cherry_pick_source_plates and, 
            plate status available
        - API_PARAM_SHOW_RETIRED_COPY_WELlS:
            show copy wells for the LCPs where:
            copy type is cherry_pick_source_plates (available or retired) or 
            copy type is library_screening_plates (retired)
        - API_PARAM_SHOW_UNFULFILLED: show only LCPs with no copy assigned
        - API_PARAM_SHOW_INSUFFICIENT: 
            for LCPs that already have a copy assigned: show those LCPs where
            the copy well volume is less than the approved request transfer volume 
        - API_PARAM_SHOW_MANUAL: show manually selected LCPs
        
        '''

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False

        is_for_detail = kwargs.pop('is_for_detail', False)
        extra_params = {}
        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id is None:
            raise BadRequest('cherry_pick_request_id is required')
        extra_params['CPR'] = cherry_pick_request_id
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cherry_pick_request_id)
        show_copy_wells = parse_val(
            param_hash.get(API_PARAM_SHOW_COPY_WELLS, False),
            API_PARAM_SHOW_COPY_WELLS, 'boolean')
        if show_copy_wells is True:
            # NOTE: add a marker to the file_name extra_params
            extra_params[API_PARAM_SHOW_COPY_WELLS] = None
        show_available_and_retired_copy_wells = parse_val(
            param_hash.get(API_PARAM_SHOW_RETIRED_COPY_WELlS, False),
            API_PARAM_SHOW_RETIRED_COPY_WELlS, 'boolean')
        if show_available_and_retired_copy_wells is True:
            extra_params[API_PARAM_SHOW_RETIRED_COPY_WELlS] = None
        show_unfulfilled = parse_val(
            param_hash.get(API_PARAM_SHOW_UNFULFILLED, False),
            API_PARAM_SHOW_UNFULFILLED, 'boolean')
        if show_unfulfilled is True:
            extra_params[API_PARAM_SHOW_UNFULFILLED] = None
        show_insufficient = parse_val(
            param_hash.get(API_PARAM_SHOW_INSUFFICIENT, False),
            API_PARAM_SHOW_INSUFFICIENT, 'boolean')
        if show_insufficient is True:
            extra_params[API_PARAM_SHOW_INSUFFICIENT] = None
        show_manual = parse_val(
            param_hash.get(API_PARAM_SHOW_MANUAL, False),
            API_PARAM_SHOW_MANUAL, 'boolean')
        if show_manual is True:
            extra_params[API_PARAM_SHOW_MANUAL] = None
        
        try:
            
            # Note: build schema for each request to use the subtype
            schema = kwargs.get('plating_schema', None)
            if schema is None: 
                schema = self.build_schema(
                    user=request.user, 
                    library_classification=cpr.screen.screen_type)
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            manual_field_includes.add('source_copy_id')
            manual_field_includes.add('source_copy_comments')
            manual_field_includes.add('selected_copy_name')
            manual_field_includes.add('selected')
            manual_field_includes.add('source_copy_usage_type')
            manual_field_includes.add('cherry_pick_request_id')
            
            # FIXME: only add the comment array if selecting alternate copies
            if '-library_plate_comment_array' not in manual_field_includes:
                manual_field_includes.add('library_plate_comment_array')
            manual_field_includes.add('library_comment_array')
            
            if show_insufficient is True:
                manual_field_includes.add('volume_approved')
                manual_field_includes.add('source_copy_well_volume')

            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema, **extra_params)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
            
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)

            # specific setup 
            base_query_tables = [
                'screen',
                'cherry_pick_request',
                'lab_cherry_pick',
                'well',
                'copy',
                'plate',
                'plate_location',
                'reagent',
                'library',
                'cherry_pick_assay_plate'
                ]
            # build the query statement
            _cpr = self.bridge['cherry_pick_request']
            _scp = self.bridge['screener_cherry_pick']
            _lcp = self.bridge['lab_cherry_pick']
            _p = self.bridge['plate']
            _pl = self.bridge['plate_location']
            _copy = self.bridge['copy']
            _screen = self.bridge['screen']
            _well = self.bridge['well']
            _reagent = self.bridge['reagent']
            _library = self.bridge['library']
            _cpap = self.bridge['cherry_pick_assay_plate']
            _cw = self.bridge['copy_well']
            _apilog = self.bridge['reports_apilog']
            _logdiff = self.bridge['reports_logdiff']
            
            _plate_comment_apilogs = \
                ApiLogResource.get_resource_comment_subquery('librarycopyplate')
            
            _library_comment_apilogs = \
                ApiLogResource.get_resource_comment_subquery('library')
            _library_comment_apilogs = _library_comment_apilogs.cte('library_comment_apilogs')

            j = join(_lcp,_cpr, 
                _lcp.c.cherry_pick_request_id == _cpr.c.cherry_pick_request_id)
            j = j.join(_scp, _lcp.c.screener_cherry_pick_id
                ==_scp.c.screener_cherry_pick_id)
            j = j.join(_screen,
                _cpr.c.screen_id == _screen.c.screen_id)
            j = j.join(_well, _well.c.well_id==_lcp.c.source_well_id)
            j = j.join(_library, _well.c.library_id==_library.c.library_id)
            j = j.join(_reagent, _lcp.c.source_well_id==_reagent.c.well_id)
            j = j.join(_cpap, _lcp.c.cherry_pick_assay_plate_id
                ==_cpap.c.cherry_pick_assay_plate_id, isouter=True)
            
            custom_columns = {
                'status': case([
                    (and_(_lcp.c.copy_id==_copy.c.copy_id,
                          _lcp.c.cherry_pick_assay_plate_id==None,), 
                        text("'%s'"%VOCAB_LCP_STATUS_SELECTED) ),
                    (and_(_lcp.c.copy_id==_copy.c.copy_id,
                          _lcp.c.cherry_pick_assay_plate_id!=None,), 
                        text("'%s'"%VOCAB_LCP_STATUS_PLATED) )],
                    else_=text("'%s'"%VOCAB_LCP_STATUS_UNFULFILLED)),
                'destination_well': (
                    case([
                        (_lcp.c.assay_plate_row!=None,
                            _concat(
                                func.chr(_lcp.c.assay_plate_row+65),
                                func.trim(func.to_char(
                                    _lcp.c.assay_plate_column+1, '00')) )
                            )],
                        else_=None)),
                'selected': case([
                    (_lcp.c.copy_id==_copy.c.copy_id, text('true') )],
                        else_=text('false')),
                'source_copywell_id': (
                    case([(_lcp.c.copy_id!=None,
                            _concat(_library.c.short_name,'/',_copy.c.name,'/',
                                _lcp.c.source_well_id)
                        )],
                        else_=None )),
                'source_copy_well_volume': case([
                    (_well.c.library_well_type=='experimental', 
                         func.coalesce(_cw.c.volume, 
                             _p.c.remaining_well_volume, _p.c.well_volume) )],
                    else_=None),
                'source_copy_well_initial_volume': case([
                    (_well.c.library_well_type=='experimental', 
                         func.coalesce(_cw.c.initial_volume,_p.c.well_volume) )],
                     else_=None),
                'source_copy_well_consumed_volume': case([
                    (_well.c.library_well_type=='experimental', 
                        func.coalesce(_cw.c.initial_volume,_p.c.well_volume)-
                            func.coalesce(_cw.c.volume, _p.c.remaining_well_volume) )],
                    else_=None),
                'library_comment_array': (
                    select([func.array_to_string(
                        func.array_agg(
                            _concat(                            
                                cast(_library_comment_apilogs.c.name,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                cast(_library_comment_apilogs.c.date_time,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                _library_comment_apilogs.c.comment)
                        ), 
                        LIST_DELIMITER_SQL_ARRAY) ])
                    .select_from(_library_comment_apilogs)
                    .where(_library_comment_apilogs.c.key==_library.c.short_name)
                ),
            }
            
            # show screener selections (may be different if from pool wells)
            if ( 'screener_well_id' in field_hash
                or 'screener_library_short_name' in field_hash):
                _screener_well = _well.alias('screener_well')
                _screener_library = _library.alias('screener_library')
                j = j.join(_screener_well, 
                    _screener_well.c.well_id==_scp.c.screened_well_id)
                j = j.join(_screener_library, 
                    _screener_well.c.library_id==_screener_library.c.library_id)
                custom_columns['screener_well_id'] = _scp.c.screened_well_id
                custom_columns['screener_library_short_name'] = _screener_library.c.short_name
            
            # pull in the copy well volume
            if ( show_copy_wells is not True 
                 and show_available_and_retired_copy_wells is not True):
                j = j.join(_copy, _lcp.c.copy_id==_copy.c.copy_id, isouter=True)
                j = j.join(_p, and_(
                    _copy.c.copy_id == _p.c.copy_id,
                    _well.c.plate_number == _p.c.plate_number), isouter=True)
                j = j.join(
                    _pl, _p.c.plate_location_id == _pl.c.plate_location_id,
                    isouter=True)
                j = j.join(_cw, and_(
                    _cw.c.well_id == _well.c.well_id,
                    _cw.c.copy_id == _copy.c.copy_id), isouter=True )

                if 'library_plate_comment_array' in field_hash:
                    ### Performance hack: limit the apilogs for the query
                    with get_engine().connect() as conn:
                        temp = (
                            select([distinct(_concat(
                                _library.c.short_name,'/',_copy.c.name,'/', 
                                cast(_p.c.plate_number,sqlalchemy.sql.sqltypes.Text))
                                )])
                            .select_from(j)
                            .where(_lcp.c.cherry_pick_request_id==cherry_pick_request_id))
                        plate_keys = set([x[0] for x in 
                            conn.execute(temp)])
                        _plate_comment_apilogs = _plate_comment_apilogs.\
                            where(_apilog.c.key.in_(plate_keys))
                        _plate_comment_apilogs = \
                            _plate_comment_apilogs.cte('plate_comment_apilogs')

                custom_columns.update({
                    'selected_copy_name': case([
                        (_lcp.c.copy_id==_copy.c.copy_id, _copy.c.name )],
                            else_=None),
                })
            else: # show all / retired copies
                _original_copy = _copy.alias('c1')
                _copyplates_available = (
                    select([
                        _lcp.c.lab_cherry_pick_id,
                        _copy.c.copy_id,
                        _p.c.plate_id,
                        _concat(
                            _library.c.short_name,'/',_copy.c.name,'/', 
                            cast(_p.c.plate_number,sqlalchemy.sql.sqltypes.Text)).label('key')])
                    .select_from(
                        _lcp.join(_well,_well.c.well_id==_lcp.c.source_well_id)
                            .join(_p, _well.c.plate_number==_p.c.plate_number)
                            .join(_copy, _p.c.copy_id==_copy.c.copy_id)
                            .join(_library, _copy.c.library_id==_library.c.library_id))
                    .where(_lcp.c.cherry_pick_request_id==cherry_pick_request_id))
                if show_available_and_retired_copy_wells is not True:
                    _copyplates_available = \
                        _copyplates_available.where(and_(
                            _copy.c.usage_type=='cherry_pick_source_plates',
                            _p.c.status=='available'))
                else:
                    # NOTE: show retired cherry pick plates: make sure not to
                    # consider these when automatically mapping
                    # _copyplates_available = \
                    #     _copyplates_available.where(and_(
                    #         _copy.c.usage_type=='cherry_pick_source_plates',
                    #         _p.c.status=='available'))
                    _copyplates_available = \
                        _copyplates_available.where(or_(
                            _copy.c.usage_type=='cherry_pick_source_plates',
                            and_(_copy.c.usage_type=='library_screening_plates',
                                _p.c.status=='retired')))
                    
                _copyplates_available = _copyplates_available.where(
                    _p.c.status.in_(['available','retired']))
                _copyplates_available = _copyplates_available.cte('copy_plates_available')
                
                j = j.join(_copyplates_available,
                    _lcp.c.lab_cherry_pick_id
                    ==_copyplates_available.c.lab_cherry_pick_id)
                j = j.join(_original_copy, _lcp.c.copy_id
                    ==_original_copy.c.copy_id, isouter=True)
                j = j.join(_copy, _copy.c.copy_id
                    ==_copyplates_available.c.copy_id)
                j = j.join(_p, _p.c.plate_id==_copyplates_available.c.plate_id)
                j = j.join(
                    _pl, _p.c.plate_location_id == _pl.c.plate_location_id,
                    isouter=True)
                j = j.join(_cw, and_(
                    _cw.c.well_id == _well.c.well_id,
                    _cw.c.copy_id == _copy.c.copy_id), isouter=True )

                if 'library_plate_comment_array' in field_hash:
                    ### Performance hack: limit the apilogs for the query
                    with get_engine().connect() as conn:
                        temp = (
                            select([distinct(_concat(
                                _library.c.short_name,'/',_copy.c.name,'/', 
                                cast(_p.c.plate_number,sqlalchemy.sql.sqltypes.Text))
                                )])
                            .select_from(j)
                            .where(_lcp.c.cherry_pick_request_id==cherry_pick_request_id))
                        plate_keys = set([x[0] for x in 
                            conn.execute(temp)])
                        _plate_comment_apilogs = _plate_comment_apilogs.\
                            where(_apilog.c.key.in_(plate_keys))
                        _plate_comment_apilogs = \
                            _plate_comment_apilogs.cte('plate_comment_apilogs')
                    
                custom_columns.update({
                    'selected_copy_name': case([
                        (_lcp.c.copy_id==_copy.c.copy_id, _copy.c.name )],
                            else_=_original_copy.c.name),
                    'source_copywell_id': (
                        _concat(_library.c.short_name,'/',_copy.c.name,'/',
                            _lcp.c.source_well_id)
                            ),
                    'status': case([
                        (and_(_lcp.c.copy_id==_copy.c.copy_id,
                              _lcp.c.cherry_pick_assay_plate_id==None,), 
                            text("'%s'"%VOCAB_LCP_STATUS_SELECTED) ),
                        (and_(_lcp.c.copy_id==_copy.c.copy_id,
                              _lcp.c.cherry_pick_assay_plate_id!=None,), 
                            text("'%s'"%VOCAB_LCP_STATUS_PLATED) )],
                        else_=text("'%s'"%VOCAB_LCP_STATUS_NOT_SELECTED)),
                })

            custom_columns['library_plate_comment_array'] = (
                select([func.array_to_string(
                    func.array_agg(
                        _concat(                            
                            cast(_plate_comment_apilogs.c.name,
                                sqlalchemy.sql.sqltypes.Text),
                            LIST_DELIMITER_SUB_ARRAY,
                            cast(_plate_comment_apilogs.c.date_time,
                                sqlalchemy.sql.sqltypes.Text),
                            LIST_DELIMITER_SUB_ARRAY,
                            _plate_comment_apilogs.c.comment)
                    ), 
                    LIST_DELIMITER_SQL_ARRAY) ])
                .select_from(_plate_comment_apilogs)
                .where(_plate_comment_apilogs.c.key== _concat(
                    _library.c.short_name,'/',_copy.c.name,'/', 
                    cast(_p.c.plate_number,sqlalchemy.sql.sqltypes.Text)))
                )
#             custom_columns['library_plate_comment_array'] = (
#                 select([func.array_to_string(
#                     func.array_agg(
#                         _concat(                            
#                             cast(_plate_comment_apilogs.c.username,
#                                 sqlalchemy.sql.sqltypes.Text),
#                             LIST_DELIMITER_SUB_ARRAY,
#                             cast(_plate_comment_apilogs.c.date_time,
#                                 sqlalchemy.sql.sqltypes.Text),
#                             LIST_DELIMITER_SUB_ARRAY,
#                             _plate_comment_apilogs.c.comment)
#                     ), 
#                     LIST_DELIMITER_SQL_ARRAY) ])
#                 .select_from(_plate_comment_apilogs)
#                 .where(_plate_comment_apilogs.c.key==_p.c.key)
#                 )
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            
            stmt = select(columns.values()).select_from(j)
            stmt = stmt.where(_cpr.c.cherry_pick_request_id==cherry_pick_request_id)

            if show_insufficient is True:
                with get_engine().connect() as conn:
                    _lcps_insufficient = (
                        select([_lcp.c.lab_cherry_pick_id,
                            func.coalesce(_cw.c.volume, 
                                    _p.c.remaining_well_volume, _p.c.well_volume)])
                        .select_from(
                            _lcp.join(_copy, _lcp.c.copy_id==_copy.c.copy_id)
                                .join(_well, _lcp.c.source_well_id==_well.c.well_id)
                                .join(_p, and_(
                                    _p.c.plate_number==_well.c.plate_number,
                                    _p.c.copy_id==_copy.c.copy_id))
                                .join(_cw, and_(
                                    _cw.c.well_id==_well.c.well_id,
                                    _cw.c.copy_id==_copy.c.copy_id),isouter=True)
                                .join(_cpr, _lcp.c.cherry_pick_request_id
                                    ==_cpr.c.cherry_pick_request_id))
                        .where(_cpr.c.cherry_pick_request_id==cherry_pick_request_id)
                        )
                    if cpr.cherry_pick_assay_plates.exists():
                        _lcps_insufficient = _lcps_insufficient.where(_cw.c.volume < 0)
                    else:
                        _lcps_insufficient = _lcps_insufficient.where(
                            _cpr.c.transfer_volume_per_well_approved 
                                > func.coalesce(_cw.c.volume, 
                                    _p.c.remaining_well_volume, _p.c.well_volume))
                    lcp_ids = set([(x[0],x[1]) for x in 
                            conn.execute(_lcps_insufficient)])
                    logger.info('insufficient lcps: %r', lcp_ids)
                    lcp_ids = set([x[0] for x in lcp_ids])
                    stmt = stmt.where(_lcp.c.lab_cherry_pick_id.in_(lcp_ids))

            if show_manual is True:
                stmt = stmt.where(_lcp.c.is_manually_selected==True)
            
            if show_unfulfilled is True:
                if ( show_copy_wells is not True 
                     and show_available_and_retired_copy_wells is not True ):
                    stmt = stmt.where(_copy.c.name==None)
                else:
                    stmt = stmt.where(_original_copy.c.name==None)
            # Ordering for well_id must be alphanumeric
            # For string field ordering, double sort as numeric and text
            order_clauses.append(text(
                "(substring({field_name}, '^[0-9]+'))::int asc " # cast to integer
                ",substring({field_name}, ':(.*$)') asc  "  # works as text
                .format(field_name='source_well_id')))
            if ( show_copy_wells is True 
                 or show_available_and_retired_copy_wells is True):
                order_clauses.extend((
                    desc('selected'),asc('source_copy_usage_type'),
                    asc('source_copy_name'),))
            
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
           
#             compiled_stmt = str(stmt.compile(
#                 dialect=postgresql.dialect(),
#                 compile_kwargs={"literal_binds": True}))
#             logger.info('compiled_stmt %s', compiled_stmt)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None))
             
        except Exception, e:
            logger.exception('on get list')
            raise e  


class CherryPickPlateResource(DbApiResource):        

    class Meta:

        queryset = CherryPickAssayPlate.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'cherrypickassayplate'
        authorization = CherryPickRequestAuthorization(resource_name)
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):

        super(CherryPickPlateResource, self).__init__(**kwargs)
        self.cherry_pick_resource = None
        
    def get_cherry_pick_resource(self):
        if self.cherry_pick_resource is None:
            self.cherry_pick_resource = CherryPickRequestResource()
        return self.cherry_pick_resource
    
    def clear_cache(self):
        DbApiResource.clear_cache(self)
        self.get_cherry_pick_resource().clear_cache()
        
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/(?P<plate_ordinal>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    @read_authorization
    def get_detail(self, request, **kwargs):

        cherry_pick_request_id = kwargs.get('cherry_pick_request_id', None)
        if not cherry_pick_request_id:
            raise NotImplementedError(
                'must provide a cherry_pick_request_id parameter')
        plate_ordinal = kwargs.get('plate_ordinal', None)
        if not plate_ordinal:
            raise NotImplementedError(
                'must provide a plate_ordinal parameter')

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, schema=None, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        if schema is None:
            raise Exception('schema not initialized')

        is_for_detail = kwargs.pop('is_for_detail', False)
        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id:
            param_hash['cherry_pick_request_id__eq'] = cherry_pick_request_id
        plate_ordinal = param_hash.pop('plate_ordinal', None)
        if plate_ordinal:
            param_hash['plate_ordinal__eq'] = plate_ordinal

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
  
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)                                  
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
            
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            _cpap = self.bridge['cherry_pick_assay_plate']
            _cpr = self.bridge['cherry_pick_request']
            _lcp = self.bridge['lab_cherry_pick']
            _screen = self.bridge['screen']
            _user_cte = ScreensaverUserResource.get_user_cte().cte('users_cpr')
            _plated_by = _user_cte.alias('plated_by')
            _screened_by = _user_cte.alias('screened_by')


            _apilog = self.bridge['reports_apilog']
            _diff = self.bridge['reports_logdiff']


            # specific setup 
            base_query_tables = [
                'cherry_pick_assay_plate',
                'cherry_pick_request',
                'screen'
            ]
        
            custom_columns = {
                'plated_by_username': _plated_by.c.username,
                'plated_by_id': _plated_by.c.screensaver_user_id,
                'plated_by_name': _plated_by.c.name,
                'screened_by_id': _screened_by.c.screensaver_user_id,
                'screened_by_username': _screened_by.c.username,
                'screened_by_name': _screened_by.c.name,
                'lab_cherry_pick_count': (
                    select([func.count(None)])
                    .select_from(_lcp)
                    .where(_lcp.c.cherry_pick_request_id
                        ==_cpr.c.cherry_pick_request_id)
                    .where(_lcp.c.cherry_pick_assay_plate_id
                        ==_cpap.c.cherry_pick_assay_plate_id)),
                'plating_comments':(
                    select([_apilog.c.comment])
                        .select_from(
                            _apilog.join(_diff, _apilog.c.id==_diff.c.log_id))
                        .where(_apilog.c.ref_resource_name=='cherrypickassayplate')
                        .where(_apilog.c.key==
                            _concat(
                                cast(_cpap.c.cherry_pick_request_id,
                                     sqlalchemy.sql.sqltypes.Text),
                                '/',
                                cast(_cpap.c.plate_ordinal,sqlalchemy.sql.sqltypes.Text)))
                        .where(_diff.c.field_key=='plating_date')
                        .order_by(desc(_apilog.c.date_time))
                        .limit(1)
                    ),
                'screening_comments':(
                    select([_apilog.c.comment])
                        .select_from(
                            _apilog.join(_diff, _apilog.c.id==_diff.c.log_id))
                        .where(_apilog.c.ref_resource_name
                            =='cherrypickassayplate')
                        .where(_apilog.c.key==
                            _concat(
                                cast(_cpap.c.cherry_pick_request_id,
                                    sqlalchemy.sql.sqltypes.Text),
                                '/',
                                cast(_cpap.c.plate_ordinal,
                                    sqlalchemy.sql.sqltypes.Text)))
                        .where(_diff.c.field_key=='screening_date')
                        .order_by(desc(_apilog.c.date_time))
                        .limit(1)
                    ),
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            # build the query statement
            j = join(_cpap, _cpr,
                _cpap.c.cherry_pick_request_id == _cpr.c.cherry_pick_request_id)
            j = j.join(_screen, _screen.c.screen_id==_cpr.c.screen_id)
            j = j.join(_plated_by, _plated_by.c.screensaver_user_id
                ==_cpap.c.plated_by_id, isouter=True)
            j = j.join(_screened_by, _screened_by.c.screensaver_user_id
                ==_cpap.c.screened_by_id, isouter=True)
            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            # compiled_stmt = str(stmt.compile(
            # dialect=postgresql.dialect(),
            # compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
#             if not order_clauses:
            stmt = stmt.order_by(
                nullslast(asc(column('cherry_pick_request_id'))),
                nullslast(asc(column('plate_ordinal'))),)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None))
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_list(self, request, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        if 'cherry_pick_request_id' not in kwargs:
            raise BadRequest('cherry_pick_request_id')
        cpr_id = kwargs['cherry_pick_request_id']
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        logger.info(
            'patch_list: cpr: %r, screen: %r...', cpr, cpr.screen.facility_id)
         
        deserialized = self.deserialize(request)
        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
        logger.debug('patch cpaps: %r', deserialized)
        
        id_attribute = schema['id_attribute']
        
        assay_plates = { cpap.plate_ordinal:cpap 
            for cpap in cpr.cherry_pick_assay_plates.all() }
        
        if not assay_plates:
            raise ValidationError(
                key='Cherry Pick Assay Plates',
                msg='plates have not been created')
        
        parent_log = self.make_log(request)
        parent_log.ref_resource_name = 'cherrypickrequest'
        parent_log.key = str(cpr.cherry_pick_request_id)
        parent_log.uri = '/'.join([
            'screen',cpr.screen.facility_id,
            parent_log.ref_resource_name,parent_log.key])        
        parent_log.save()
        
        original_cpr = self.get_cherry_pick_resource()\
            ._get_detail_response_internal(**{
                'cherry_pick_request_id': cpr_id })
        original_cpap_data = self._get_list_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        
        plated_cpaps = []
        screened_cpaps = []
        for update_cpap in deserialized:
            initializer_dict = self.parse(update_cpap, schema=schema)
            plate_ordinal = parse_val(
                update_cpap.get('plate_ordinal', None), 'plate_ordinal', 'integer')
            if plate_ordinal is None:
                raise ValidationError(key='plate_ordinal', msg='required')
            
            cpap = assay_plates.get(plate_ordinal, None)
            if cpap is None:
                raise ValidationError(
                    key='plate_ordinal', 
                    msg='%d not found'%plate_ordinal)
            
            plating_date = initializer_dict.get('plating_date', None)
            screening_date = initializer_dict.get('screening_date', None)
            
            if plating_date is None and screening_date is None:
                msg = 'must submit either "plating_date" or "screening_date"'
                raise ValidationError({
                    'plating_date': msg, 'screening_date': msg })
            
            if plating_date is not None:
                plated_by_username = initializer_dict.get('plated_by_username', None)
                plated_by_id = initializer_dict.get('plated_by_id', None)
                if plated_by_username is None and plated_by_id is None:
                    raise ValidationError({
                        'plated_by_username': 'required',
                        'plated_by_id': 'required',
                    })
                try:
                    if plated_by_id is not None:
                        user = ScreensaverUser.objects.get(screensaver_user_id=plated_by_id)
                    else:
                        user = ScreensaverUser.objects.get(username=plated_by_username)
                    
                    # FIXME: check that user has cherrypickrequest/write permission
                    
                    cpap.plated_by = user
                    cpap.plating_date = plating_date
                    cpap.save()
                    plated_cpaps.append(cpap)
                except ObjectDoesNotExist:
                    if plated_by_id is not None:
                        raise ValidationError(
                            key='plated_by_id', 
                            msg='User: "%s" not found' % plated_by_id)
                    else:
                        raise ValidationError(
                            key='plated_by_username', 
                            msg='User: "%s" not found' % plated_by_username)
                
            if screening_date is not None:
                if cpap.plating_date is None:
                    raise ValidationError(
                        key='screening_date',
                        msg='"plating_date" must be set before "screening_date" can be set')
                screened_by_username = initializer_dict.get('screened_by_username', None)
                screened_by_id = initializer_dict.get('screened_by_id', None)
                if screened_by_username is None and screened_by_id is None:
                    raise ValidationError({
                        'screened_by_username': 'required',
                        'screened_by_id': 'required',
                    })
                try:
                    if screened_by_id is not None:
                        user = ScreensaverUser.objects.get(screensaver_user_id=screened_by_id)
                    else:
                        user = ScreensaverUser.objects.get(username=screened_by_username)
                    
                    if user not in cpr.screen.get_screen_users():
                        raise ValidationError(
                            key='screened_by_username',
                            msg='Not a member of this screen')
                    
                    cpap.screened_by = user
                    cpap.screening_date = screening_date
                    cpap.save()
                    screened_cpaps.append(cpap)
                except ObjectDoesNotExist:
                    if screened_by_id is not None:
                        raise ValidationError(
                            key='screened_by_id', 
                            msg='User: "%s" not found' % screened_by_id)
                    else:
                        raise ValidationError(
                            key='screened_by_username', 
                            msg='User: "%s" not found' % screened_by_username)
                # NOTE: last_screening_activity_date is handled in log_patch below
                # parent_log.diffs['last_screening_activity_date'] = [
                #     original_cpr.get('last_screening_activity_date', None),
                #     plating_date] 

        if plated_cpaps:
            # TODO: create a cherry pick liquid transfer (plating) activity 
            # - for legacy activity tracking system (201703)
            # TODO: test

            cplt = CherryPickLiquidTransfer.objects.create(
                status="Successful",
                screen=cpr.screen,
                cherry_pick_request = cpr,
                volume_transferred_per_well_from_library_plates
                    =cpr.transfer_volume_per_well_approved,
                apilog_uri=parent_log.log_uri,
                date_of_activity=plated_cpaps[0].plating_date,
                performed_by = plated_cpaps[0].plated_by,
                created_by = self.get_request_user(request)
            )
            
            # Set the legacy cplt link; NOTE: this is not used in V2
            cpap.cherry_pick_liquid_transfer = cplt
            cpap.save()
            # TODO: set api log to activity
            # cplt.log_id = parent_log
        if screened_cpaps:
            # TODO: create a cherry pick liquid transfer (plating) activity 
            # - for legacy activity tracking system (201703)
            # TODO: test
            cp_screening = CherryPickScreening.objects.create(
                screen = cpr.screen,
                cherry_pick_request = cpr,
                assay_protocol = '',
                assay_protocol_type = '',
                volume_transferred_per_well_to_assay_plates
                    =cpr.transfer_volume_per_well_approved,
                apilog_uri=parent_log.log_uri,
                date_of_activity=screened_cpaps[0].screening_date,
                performed_by = screened_cpaps[0].screened_by,
                created_by = self.get_request_user(request)
                )
            # TODO: set api log to activity
            # cplt.log_id = parent_log
            
        # Log
        self.get_cherry_pick_resource().clear_cache()
        new_cpr = self.get_cherry_pick_resource()._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        new_cpap_data = self._get_list_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        self.get_cherry_pick_resource().log_patch(
            request, original_cpr, new_cpr, log=parent_log)
        patch_logs = self.log_patches(
            request, original_cpap_data, new_cpap_data, schema=schema, 
            parent_log=parent_log)
        meta = {}
        plated_changed_count = len([
            x for x in patch_logs if 'plating_date' in x.diffs])
        if plated_changed_count > 0:
            meta[API_MSG_CPR_PLATES_PLATED] = plated_changed_count
        screened_changed_count = len([
            x for x in patch_logs if 'screening_date' in x.diffs])
        if screened_changed_count > 0:
            meta[API_MSG_CPR_PLATES_SCREENED] = screened_changed_count

        parent_log.json_field = meta
        parent_log.save()
        
        return self.build_response(
            request, { API_RESULT_META: meta }, 
            response_class=HttpResponse, **kwargs)
                
class LibraryCopyResource(DbApiResource):

    class Meta:

        queryset = Copy.objects.all().order_by('name')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'librarycopy'
        authorization = UserGroupAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):

        super(LibraryCopyResource, self).__init__(**kwargs)
        self.librarycopyplate_resource = None
        
    def get_plate_resource(self):
        if self.librarycopyplate_resource is None:
            self.librarycopyplate_resource = LibraryCopyPlateResource()
        return self.librarycopyplate_resource

    def clear_cache(self):
        logger.info('clear copy caches')
        DbApiResource.clear_cache(self)

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
           url((r"^(?P<resource_name>%s)"
                 r"/(?P<library_short_name>[\w.\-\+: ]+)"
                 r"/(?P<copy_name>[^/]+)%s$")  
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('search'), name="api_search"),
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<library_short_name>[\w.\-\+: ]+)"
                 r"/(?P<copy_name>[^/]+)"
                 r"/plate%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_librarycopyplateview'),
                name="api_dispatch_librarycopy_plateview"),
        ]    

    def dispatch_librarycopyplateview(self, request, **kwargs):

        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)    
        
    @read_authorization
    def get_detail(self, request, **kwargs):

        library_short_name = kwargs.get(u'library_short_name', None)
        if not library_short_name:
            raise NotImplementedError(
                'must provide a library_short_name parameter')
        copy_name = kwargs.get('copy_name', None)
        if not copy_name:
            raise NotImplementedError(
                'must provide a copy "copy_name" parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)
    
    def build_list_response(self, request, schema=None, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        if schema is None:
            raise Exception('schema not initialized')
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        library_short_name = param_hash.pop('library_short_name',
            param_hash.get('library_short_name__eq', None))
        if not library_short_name:
            logger.info('no library_short_name provided')
        else:
            param_hash['library_short_name__eq'] = library_short_name
        name = param_hash.pop('name', param_hash.get('name', None))
        if name:
            param_hash['name__eq'] = name
            
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            if is_for_detail:
                manual_field_includes.add('has_copywell_concentrations')
                manual_field_includes.add('has_copywell_volumes')

            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)

            # if filter_expression is None:
            #     raise InformationError(
            #         key='Input filters ',
            #         msg='Please enter a filter expression to begin')

            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
            rowproxy_generator = \
               self._meta.authorization.get_row_property_generator(
                   request.user, field_hash, rowproxy_generator)

            # specific setup 

            _c = self.bridge['copy']
            _l = self.bridge['library']
            _ap = self.bridge['assay_plate']
            _p = self.bridge['plate']
            _cw = self.bridge['copy_well']
            _ls = self.bridge['library_screening']
            _plate_statistics = (
                LibraryCopyPlateResource.get_plate_copywell_statistics_cte(
                    filter_hash=filter_hash)
                    .cte('plate_statistics'))
            _plate_screening_statistics = (
                LibraryCopyPlateResource.get_plate_screening_statistics_cte(
                    filter_hash=filter_hash)
                    .cte('plate_screening_statistics'))
            _copy_statistics = (
                select([
                    _plate_statistics.c.copy_id,
                    case([
                        ( func.avg(_plate_statistics.c.avg_well_remaining_volume) != None,
                           func.avg(_plate_statistics.c.avg_well_remaining_volume))], 
                        else_=func.avg(_plate_statistics.c.remaining_well_volume))
                        .label('avg_plate_volume'),
                    case([
                        ( func.min(_plate_statistics.c.min_well_remaining_volume) != None,
                           func.min(_plate_statistics.c.min_well_remaining_volume))
                        ], else_=func.min(_plate_statistics.c.remaining_well_volume))
                        .label('min_plate_volume'),
                    case([
                        ( func.max(_plate_statistics.c.max_well_remaining_volume) != None,
                           func.max(_plate_statistics.c.max_well_remaining_volume))
                        ],else_=func.max(_plate_statistics.c.remaining_well_volume))
                        .label('max_plate_volume'),
                    func.min(_plate_statistics.c.min_mg_ml_concentration)
                        .label('min_mg_ml_concentration'),
                    func.max(_plate_statistics.c.max_mg_ml_concentration)
                        .label('max_mg_ml_concentration'),
                    func.min(_plate_statistics.c.min_molar_concentration)
                        .label('min_molar_concentration'),
                    func.max(_plate_statistics.c.max_molar_concentration)
                        .label('max_molar_concentration'),
                     ])
                 .select_from(_plate_statistics)
                 .group_by(_plate_statistics.c.copy_id)
             ).cte('copy_statistics')

            _copy_screening_statistics = (
                select([
                    _plate_screening_statistics.c.copy_id,
                    func.max(_plate_screening_statistics.c.last_date_screened)
                        .label('last_date_screened'),
                    func.min(_plate_screening_statistics.c.first_date_screened)
                        .label('first_date_screened')])
                .select_from(_plate_screening_statistics)
                .group_by(_plate_screening_statistics.c.copy_id)
                ).cte('copy_screening_statistics')
            
            custom_columns = {
                'copy_plate_count': (
                    select([func.count(None)])
                    .select_from(_p).where(_p.c.copy_id==text('copy.copy_id'))),
                    
                'avg_plate_volume': _copy_statistics.c.avg_plate_volume,
                'min_plate_volume': _copy_statistics.c.min_plate_volume,
                'max_plate_volume': _copy_statistics.c.max_plate_volume,
                'avg_plate_screening_count': (
                    select([func.avg(_p.c.screening_count)])
                    .select_from(_p)
                    .where(_p.c.copy_id == literal_column('copy.copy_id'))),
                'avg_plate_cp_screening_count': (
                    select([func.avg(_p.c.cplt_screening_count)])
                    .select_from(_p)
                    .where(_p.c.copy_id == literal_column('copy.copy_id'))),
                'min_mg_ml_concentration': _copy_statistics.c.min_mg_ml_concentration,
                'max_mg_ml_concentration': _copy_statistics.c.max_mg_ml_concentration,
                'min_molar_concentration': _copy_statistics.c.min_molar_concentration,
                'max_molar_concentration': _copy_statistics.c.min_molar_concentration,
                
                'first_date_screened': _copy_screening_statistics.c.first_date_screened,
                'last_date_screened': _copy_screening_statistics.c.last_date_screened,

                'primary_plate_location': (
                    literal_column('\n'.join([
                    "( select room || '-' || freezer || '-' || shelf || '-' || bin ",
                    '    from plate_location pl ' ,
                    '    where pl.plate_location_id',
                    '=copy.primary_plate_location_id) ']))
                    .label('primary_plate_location')),
                'plate_locations': literal_column('\n'.join([
                    '(select count(distinct(plate_location_id)) ',
                    '    from plate p',
                    '    where p.copy_id = copy.copy_id ) '])).\
                    label('plate_locations'),
                'plates_available': literal_column('\n'.join([
                    '(select count(p)',
                    '    from plate p ',
                    '    where p.copy_id=copy.copy_id',
                    "    and p.status in ('available') ) "])).\
                    label('plates_available'),
                'has_copywell_volumes': (
                    exists(select([None])
                        .select_from(_cw)
                        .where(_cw.c.copy_id==_c.c.copy_id)
                        .where(_cw.c.initial_volume!=None)
                        )),
                'has_copywell_concentrations':
                    exists(select([None])
                        .select_from(_cw)
                        .where(_cw.c.copy_id==_c.c.copy_id)
                        .where(or_(
                            _cw.c.mg_ml_concentration!=None,
                            _cw.c.molar_concentration!=None)
                        )),
                }
            
            base_query_tables = ['copy', 'library']
 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            # build the query statement

            
            j = join(_c, _l, _c.c.library_id == _l.c.library_id)
            if set(['avg_remaining_volume','min_remaining_volume',
                'max_remaining_volume']) | set(field_hash.keys()):
                j = j.outerjoin(
                    _copy_statistics,_c.c.copy_id==_copy_statistics.c.copy_id)     
            j = j.outerjoin(_copy_screening_statistics,
                _c.c.copy_id==_copy_screening_statistics.c.copy_id)
            stmt = select(columns.values()).select_from(j)
            
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            if not order_clauses:
                stmt = stmt.order_by("library_short_name", "copy_name")
 
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])

            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)

            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash, param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
            
        except Exception, e:
            logger.exception('on get list')
            raise e   

    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_list must be implemented')
                
    @write_authorization
    @un_cache
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        ScreensaverUser.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        initializer_dict = {}

        # TODO: wrapper for parsing
        logger.debug('fields: %r, deserialized: %r', fields.keys(), deserialized)
        for key in fields.keys():
            if deserialized.get(key, None) is not None:
                initializer_dict[key] = parse_val(
                    deserialized.get(key, None), key, fields[key]['data_type']) 

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        
        try:
            short_name = id_kwargs['library_short_name']
            try:
                library = Library.objects.get(short_name=short_name)
            except ObjectDoesNotExist:
                msg = 'library_short_name not found: %r' % short_name
                logger.info(msg);
                raise Http404(msg)
            
            patch = False
            librarycopy = None
            try:
                librarycopy = Copy.objects.get(
                    name=id_kwargs['copy_name'], library=library)
                patch = True
                errors = self.validate(deserialized, schema=schema, patch=True)
                if errors:
                    raise ValidationError(errors)
            except ObjectDoesNotExist:
                librarycopy = Copy.objects.create(
                    name=id_kwargs['copy_name'], library=library)
                errors = self.validate(deserialized, schema=schema, patch=False)
                if errors:
                    raise ValidationError(errors)
                librarycopy.save()
            initializer_dict = {}
            for key in fields.keys():
                if key in deserialized:
                    initializer_dict[key] = parse_val(
                        deserialized.get(key, None), key,
                        fields[key]['data_type']) 
            if initializer_dict:
                logger.debug('initializer dict: %s', initializer_dict)
                for key, val in initializer_dict.items():
                    if hasattr(librarycopy, key):
                        setattr(librarycopy, key, val)
            else:
                logger.info(
                    'no (basic) library copy fields to update %s', deserialized)
            
            librarycopy.save()
            logger.info('librarycopy saved: %r', librarycopy)

            if patch is False:
            
                logger.info('create librarycopyplates for range: %d-%d, copy: %s',
                    librarycopy.library.start_plate, 
                    librarycopy.library.end_plate, librarycopy.name )
    
                initial_plate_status = deserialized.get(
                    'initial_plate_status', None)
                if initial_plate_status:
                    self.get_plate_resource().validate({
                        'status': initial_plate_status}, schema=schema)
            
                initial_plate_well_volume = deserialized.get(
                    'initial_plate_well_volume', None)
                if initial_plate_well_volume:
                    initial_plate_well_volume = parse_val(
                        initial_plate_well_volume, 
                        'initial_plate_well_volume', 'decimal')
                elif patch is False:
                    raise ValidationError(
                        key='initial_plate_well_volume',
                        msg='required')
                initial_plate_mg_ml_concentration = \
                    deserialized.get('initial_plate_mg_ml_concentration', None)
                if initial_plate_mg_ml_concentration:
                    initial_plate_mg_ml_concentration = parse_val(
                        initial_plate_mg_ml_concentration, 
                        'initial_plate_mg_ml_concentration', 'decimal')
                initial_plate_molar_concentration = \
                    deserialized.get('initial_plate_molar_concentration', None)
                if initial_plate_molar_concentration:
                    initial_plate_molar_concentration = parse_val(
                        initial_plate_mg_ml_concentration, 
                        'initial_plate_molar_concentration', 'decimal')
                if ( initial_plate_molar_concentration is not None and
                     initial_plate_mg_ml_concentration is not None ):
                    msg = ('May only set one of '
                        '["initial_plate_mg_ml_concentration",'
                        '"initial_plate_molar_concentration"]')
                    raise ValidationError({
                        'initial_plate_mg_ml_concentration': msg,
                        'initial_plate_molar_concentration': msg })     
                
                library = librarycopy.library
                
                if ( initial_plate_molar_concentration is None and 
                     initial_plate_mg_ml_concentration is None ):
                    # use the well concentrations to set the values
                    # case 1: all wells use the same concentration value
                    
                    mg_mls = ( 
                        library.well_set.all()
                            .filter(library_well_type='experimental')
                            .distinct('mg_ml_concentration')
                            .values_list('mg_ml_concentration', flat=True))
                    if len(mg_mls) == 1:
                        initial_plate_mg_ml_concentration = mg_mls[0]
                    molars = ( 
                        library.well_set.all()
                            .filter(library_well_type='experimental')
                            .distinct('molar_concentration')
                            .values_list('molar_concentration', flat=True))
                    if len(molars) == 1:
                        initial_plate_molar_concentration = molars[0]
                
                    if ( initial_plate_molar_concentration is not None and
                         initial_plate_mg_ml_concentration is not None ):
                        msg = ('Invalid library well concentrations: both '
                            '["mg_ml_concentration",'
                            '"molar_concentration"] are set')
                        logger.warn(msg)
                        initial_plate_molar_concentration = None
                        initial_plate_mg_ml_concentration = None
                        
                logger.debug('create plates start: %d, end: %d',
                    library.start_plate, library.end_plate)
                plates_created = 0
                copywells_created = 0
                for x in range(library.start_plate, library.end_plate + 1):
                    p = Plate.objects.create(
                        copy=librarycopy, 
                        plate_number=x,
                        screening_count=0,
                        status='not_specified',
                        well_volume=initial_plate_well_volume,
                        remaining_well_volume=initial_plate_well_volume,
                        mg_ml_concentration=initial_plate_mg_ml_concentration,
                        molar_concentration=initial_plate_molar_concentration )
                    
                    if initial_plate_status:
                        p.status = initial_plate_status
                        if p.status == 'available':
                            p.date_plated = _now()
                            # TODO: set retired if retired status
                    p.save()
                    plates_created += 1
                    if len(mg_mls)>1 or len(molars)>1:
                        # Library contains specific well concentrations, must 
                        # create copy_wells to reflect these concentrations
                        logger.info('creating copywells for plate: %r', p)
                        for well in library.well_set.all().filter(
                                plate_number=p.plate_number):
                            CopyWell.objects.create(
                                well=well, copy=librarycopy, plate=p,
                                mg_ml_concentration=well.mg_ml_concentration,
                                molar_concentration=well.molar_concentration)
                            copywells_created += 1
                    logger.debug('saved plate: %r', p.plate_number)
                logger.info('created %d plates, %d copywells', 
                    plates_created, copywells_created)
            
            # TODO: return Result meta for # plates, copies created
            logger.info('patch_obj done for librarycopy: %r', librarycopy)
            return { API_RESULT_OBJ: librarycopy }
            
        except Exception, e:
            logger.exception('on patch detail')
            raise e  

class PublicationAuthorization(ScreenAuthorization):
    # If the user can see screens, they can then see publications
    
    def filter(self, user, filter_expression):
        return UserGroupAuthorization.filter(self, user, filter_expression)
    def get_row_property_generator(self, user, fields, extant_generator):
        return extant_generator
    
    def has_publication_read_authorization(self, user, publication_id):
        return True
    
    # def filter(self, user, filter_expression):
    #     
    #     if self.is_restricted_view(user):
    #     
    #         screen_access_dict = self.get_screen_access_level_table(user.username)
    # 
    #         allowed_screens = { facility_id:_dict for facility_id,_dict
    #             in screen_access_dict if _dict['user_access_level_granted'] in [2,3] }    
    # 
    #         auth_filter = column('screen_facility_id').in_(allowed_screens.keys())
    #         if filter_expression is not None:
    #             filter_expression = and_(filter_expression, auth_filter)
    #         else:
    #             filter_expression = auth_filter
    #         
    #     return filter_expression
    
class PublicationResource(DbApiResource):

    class Meta:

        queryset = Publication.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'publication'
        authorization = PublicationAuthorization(resource_name)
        serializer = LimsSerializer()
        always_return_data = True

    def __init__(self, **kwargs):
        
        self.attached_file_resource = None
        self.screen_resource = None
        super(PublicationResource, self).__init__(**kwargs)
        
    def get_attached_file_resource(self):
        if not self.attached_file_resource:
            self.attached_file_resource = AttachedFileResource()
        return self.attached_file_resource

    def get_screen_resource(self):
        if not self.screen_resource:
            self.screen_resource = ScreenResource()
        return self.screen_resource

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<publication_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ] 
        
    @read_authorization
    def get_detail(self, request, **kwargs):

        publication_id = kwargs.get('publication_id', None)
        if not publication_id:
            raise NotImplementedError(
                'must provide a publication_id parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        logger.info('params: %r', param_hash.keys())
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        screen_facility_id = param_hash.pop('screen_facility_id', None)
        if screen_facility_id:
            param_hash['screen_facility_id__eq'] = screen_facility_id
        publication_id = param_hash.pop('publication_id', None)
        if publication_id:
            param_hash['publication_id__eq'] = publication_id
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)

            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)

                  
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
               self._meta.authorization.get_row_property_generator(
                   request.user, field_hash, rowproxy_generator)
 
            # specific setup
            _publication = self.bridge['publication']
            _af = self.bridge['attached_file']
            _screen = self.bridge['screen']
            
            j = _publication
            j = j.join(
                _af, _publication.c.publication_id==_af.c.publication_id, 
                isouter=True)
            j = j.join(
                _screen, _publication.c.screen_id==_screen.c.screen_id,
                isouter=True)
            
            custom_columns = {
                'lookup_pmid': literal_column("'lookup_pmid'"),
                }

            base_query_tables = ['publication', 'attached_file', 'screen' ] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            if not order_clauses and filter_expression is None:
                _alias = Alias(stmt)
                stmt = select([text('*')]).select_from(_alias)
            stmt = stmt.order_by('-publication_id')
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
             
        except Exception, e:
            logger.exception('on get_list %s' % self._meta.resource_name)
            raise e  

    def put_list(self, request, **kwargs):
        raise NotImplementedError(
            "Put list is not implemented for Publications")
    
    def put_detail(self, request, **kwargs):
        raise NotImplementedError(
            "Post detail is not implemented for Publications")
    
    def patch_list(self, request, **kwargs):
        raise NotImplementedError(
            "Patch list is not implemented for Publications")
    
    def patch_detail(self, request, **kwargs):
        raise NotImplementedError(
            "Patch detail is not implemented for Publications")
    
    @write_authorization
    @un_cache        
    @transaction.atomic
    def post_detail(self, request, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
                
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        fields = schema['fields']
        id_attribute = resource = schema['id_attribute']

        initializer_dict = {}
        for key in fields.keys():
            if key in param_hash:
                initializer_dict[key] = parse_val(
                    param_hash.get(key, None), key, fields[key]['data_type']) 

        parent_fields = set(['screen_facility_id', 'reagent_id'])
        initializer_keys = set(initializer_dict.keys())
        if ( parent_fields.isdisjoint(initializer_keys) ):
            msg='must provide one of: %r' % parent_fields
            raise ValidationError({
                'screen_facility_id': msg,
                'reagent_id': msg 
            })
        if (len(parent_fields & initializer_keys)>1):
            logger.warn('too many parent_fields: %r', 
                (parent_fields & initializer_keys))
            msg='Only one of parent fields allowed: %r' % parent_fields
            raise ValidationError({
                'screen_facility_id': msg,
                'reagent_id': msg
            })
        
        screen_facility_id = initializer_dict.pop('screen_facility_id', None)
        if screen_facility_id:
            try:
                screen = Screen.objects.get(
                    facility_id=screen_facility_id)
                initializer_dict['screen'] = screen
            except ObjectDoesNotExist:
                raise Http404('screen_facility_id %r does not exist' 
                    % screen_facility_id)
        reagent_id = initializer_dict.pop('reagent_id', None)
        if reagent_id: 
            try: 
                reagent = Reagent.objects.get(
                    reagent_id=reagent_id)
                initializer_dict['reagent'] = reagent
            except ObjectDoesNotExist:
                raise Http404('reagent_id does not exist: %s' % reagent_id)
        
        if ( 'pubmed_id' not in initializer_dict 
                and 'pubmed_central_id' not in initializer_dict):
            msg = 'must specify either "Pubmed ID" or "Pubmed CID"'
            raise ValidationError({
                'pubmed_id': msg,
                'pubmed_central_id': msg })

        publication = Publication.objects.create()
        for key, val in initializer_dict.items():
            if hasattr(publication, key):
                setattr(publication, key, val)
            else:
                logger.warn(
                    'no such attribute on publication: %s:%r' % (key, val))

        parent_log = kwargs.get('parent_log',None)
        if parent_log is None and publication.screen is not None:
            parent_log = self.log_to_screen(
                publication.screen, publication, request)
            kwargs['parent_log'] = parent_log
            
        publication.save()

        attached_file = request.FILES.get('attached_file', None)
        if attached_file:
            logger.info('create attached file for publication: %r', publication)
            
            af_request = HttpRequest()
            af_request.META['HTTP_ACCEPT'] = JSON_MIMETYPE
            af_request.user = request.user
            self.get_attached_file_resource().post_detail(
                af_request, **{ 
                    'attached_file': attached_file,
                    'publication_id': publication.publication_id,
                    'type': 'publication',
                    'parent_log': parent_log })
            
        kwargs_for_log = { 'parent_log': parent_log }
        for id_field in id_attribute:
            val = getattr(publication, id_field,None)
            if val:
                kwargs_for_log['%s' % id_field] = val
        logger.info('get new data: %r', kwargs_for_log)
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        log = self.log_patch(
            request, None, new_data, full_create_log=True, **kwargs_for_log)
        if log:
            log.save()
        logger.info('created log: %r', log)
        return self.build_response(request, new_data, status_code=201)

    def log_to_screen(self, screen, publication, request, is_delete=False):
        screen_resource = self.get_screen_resource()
        _screen_data = \
            screen_resource._get_detail_response_internal(**{
                'limit': 10,
                'facility_id': screen.facility_id,
                'exact_fields': ['publications'],
                })
        parent_log = screen_resource.make_log(request)
        parent_log.api_action = API_ACTION_PATCH
        parent_log.key = screen.facility_id
        parent_log.uri = '/'.join([parent_log.ref_resource_name,parent_log.key])
        parent_log.diff_keys = ['publications'];
        current_pubs = _screen_data['publications'] or []
        if is_delete:
            new_pubs = [x for x in current_pubs 
                if x != str(publication.title)]
        else:
            new_pubs = list(current_pubs)
            new_pubs.append(publication.title)
        parent_log.diffs = { 'publications': [current_pubs,new_pubs]}
        parent_log.save()
        logger.info('created parent log: %r', parent_log)
        return parent_log
    
    @write_authorization
    @un_cache   
    @transaction.atomic     
    def delete_detail(self, request, **kwargs):

        logger.info('delete publication...')
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_attribute = schema['id_attribute']
    
        publication_id = kwargs.get('publication_id', None)
        if not publication_id:
            NotImplementedError('must provide a publication_id parameter')
        
        publication = Publication.objects.get(publication_id=publication_id)

        parent_log = kwargs.get('parent_log',None)
        if parent_log is None and publication.screen is not None:
            parent_log = self.log_to_screen(
                publication.screen, publication, request, is_delete=True)

        kwargs_for_log = {}
        for id_field in id_attribute:
            if kwargs.get(id_field,None):
                kwargs_for_log[id_field] = kwargs[id_field]
        logger.debug('delete detail: %s' %(kwargs_for_log))
        if not kwargs_for_log:
            raise Exception('required id keys %s' % id_attribute)
        else:
            try:
                original_data = self._get_detail_response_internal(**kwargs_for_log)
            except Exception as e:
                logger.exception('original state not obtained')
                raise
#                 original_data = {}

        publication.delete()

        logger.info('deleted: %s' %kwargs_for_log)
        log_comment = None
        if HEADER_APILOG_COMMENT in request.META:
            log_comment = request.META[HEADER_APILOG_COMMENT]
        
        log = self.make_log(request)
        log.ref_resource_name = self._meta.resource_name
        log.key = '/'.join([str(original_data[x]) for x in id_attribute])
        log.uri = '/'.join([self._meta.resource_name,log.key])
        log.parent_log = parent_log
        log.api_action = API_ACTION_DELETE
        log.diffs = { k:[v,None] for k,v in original_data.items()}
        log.save()
        logger.info('delete, api log: %r', log)

        return HttpResponse(status=204)

class AttachedFileAuthorization(ScreenAuthorization):
    '''
    Extend ScreenAuthorization for convenient access to sharing levels
    '''
    
    def filter(self, user, filter_expression):
        if self.is_restricted_view(user) is False:
            return filter_expression

        screensaver_user = ScreensaverUser.objects.get(username=user.username)
        my_screens = self.get_user_screens(screensaver_user)
        
        auth_filter = or_(
            column('screen_facility_id').in_(
                [screen.facility_id for screen in my_screens]),
            column('username') == user.username)
            
        if filter_expression is not None:
            filter_expression = and_(filter_expression, auth_filter)
        else:
            filter_expression = auth_filter

        return filter_expression
        
    def get_row_property_generator(self, user, fields, extant_generator):
        return extant_generator
    
    def has_file_read_authorization(self, user, attached_file_id):
        logger.info('has_file_read_authorization: %r, %r', user, attached_file_id)
        is_restricted = self.is_restricted_view(user)
        logger.info('1...%r', is_restricted)
        if is_restricted is not True:
            return True
        screensaver_user = ScreensaverUser.objects.get(username=user.username)
        attached_file = AttachedFile.objects.get(attached_file_id=attached_file_id)
        
        logger.info('2...')
        if attached_file.screensaver_user == screensaver_user:
            return True
        elif attached_file.screen:
            authorized_screens = self.get_read_authorized_screens(screensaver_user)
            return attached_file.screen in authorized_screens
        if attached_file.screensaver_user is None and attached_file.screen is None:
            return True
        logger.info('return False')
        return False
        
class AttachedFileResource(DbApiResource):

    class Meta:

        queryset = AttachedFile.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'attachedfile'
        authorization = AttachedFileAuthorization(resource_name)
        serializer = LimsSerializer()
        always_return_data = True

    def __init__(self, **kwargs):
        
        super(AttachedFileResource, self).__init__(**kwargs)
        self.screen_resource = None
        self.user_resource = None
        
    def get_screen_resource(self):
        if self.screen_resource is None:
            self.screen_resource = ScreenResource()
        return self.screen_resource

    def get_user_resource(self):
        if self.user_resource is None:
            self.user_resource = ScreensaverUserResource()
        return self.user_resource
    
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<attached_file_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/user/(?P<screensaver_user_id>([\d]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/user/(?P<username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url((r"^(?P<resource_name>%s)/user/(?P<screensaver_user_id>([\d]+))" 
                 r"/(?P<attached_file_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/user/(?P<username>([\w]+))" 
                 r"/(?P<attached_file_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ] 
        
    def put_list(self, request, **kwargs):
        raise NotImplementedError(
            "Put list is not implemented for AttachedFiles")
    
    def put_detail(self, request, **kwargs):
        raise NotImplementedError(
            "Post detail is not implemented for AttachedFiles")
    
    def patch_list(self, request, **kwargs):
        raise NotImplementedError(
            "Patch list is not implemented for AttachedFiles")
    
    def patch_detail(self, request, **kwargs):
        raise NotImplementedError(
            "Patch detail is not implemented for AttachedFiles")
    
    def log_to_screen(self, screen, af, request, is_delete=False):
        '''
        Creates a parent log for attached_files for screens
        '''
        screen_resource = self.get_screen_resource()
        _screen_data = \
            screen_resource._get_detail_response_internal(**{
                'limit': 10,
                'facility_id': screen.facility_id,
                'exact_fields': ['attached_files'],
                })
        parent_log = screen_resource.make_log(request)
        parent_log.api_action = API_ACTION_PATCH
        parent_log.key = str(screen.facility_id)
        parent_log.uri = '/'.join([parent_log.ref_resource_name,parent_log.key])
        parent_log.diff_keys = ['attached_files'];
        current_files = _screen_data['attached_files'] or []
        if is_delete:
            new_files = [x for x in current_files 
                if x != str(af.filename)]
        else:
            new_files = list(current_files)
            new_files.append(af.filename)
        parent_log.diffs = { 'attached_files': [current_files,new_files]}
        parent_log.save()
        return parent_log

    def log_to_user(self, user, af, request, is_delete=False):
        '''
        Creates a parent log for attached_files for screens
        '''
        user_resource = self.get_screen_resource()
        parent_log = user_resource.make_log(request)
        parent_log.api_action = API_ACTION_PATCH
        parent_log.key = str(user.screensaver_user_id)
        parent_log.uri = '/'.join([parent_log.ref_resource_name,parent_log.key])
        parent_log.diff_keys = ['attached_files'];
        if is_delete:
            parent_log.diffs =  { 'attached_fiels': [af.filename,None] }
        else:
            parent_log.diffs =  { 'attached_fiels': [None, af.filename] }
        parent_log.save()
        return parent_log

    @write_authorization
    @un_cache        
    @transaction.atomic
    def post_detail(self, request, **kwargs):
        ''' 
        Custom post_detail: 
            - attached file is not editable, so no logging of former state
            - custom deserialization of the attached file from form parameters
        '''
        logger.info('post attached file')
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        fields = schema['fields']
        id_attribute = schema['id_attribute']
        initializer_dict = {}
        for key in fields.keys():
            if param_hash.get(key, None) is not None:
                initializer_dict[key] = parse_val(
                    param_hash.get(key, None), key, fields[key]['data_type']) 

        attached_file = request.FILES.get('attached_file', None)
        if not attached_file:
            attached_file = param_hash.get('attached_file', None)
        if not attached_file:
            # FIXME: "contents" has been eliminated
            contents = param_hash.get('contents', None)
            filename = param_hash.get('filename', None)
            if not (contents and filename):
                raise ValidationError(
                    key = 'attached_file',
                    msg = ('must provide either "attached_file" '
                           'or "contents+filename" parameters'))
            contents = contents.encode('utf-8')
            initializer_dict['contents'] = contents
            initializer_dict['filename'] = filename
        else:
            contents = attached_file.read()
            filename = attached_file.name
            if param_hash.get('filename', None):
                filename = param_hash.get('filename', None)
            initializer_dict['contents'] = contents
            initializer_dict['filename'] = filename

        if not 'created_by_username' in initializer_dict:
            initializer_dict['created_by_username'] = request.user.username
        try:
            admin_user = ScreensaverUser.objects.get(
                username=initializer_dict['created_by_username'])
            logger.debug('using admin_user %s' % admin_user)
            initializer_dict['created_by'] = admin_user
        except ObjectDoesNotExist:
            logger.exception('created_by_username does not exist: %s',
                initializer_dict['created_by_username'])
            raise ValidationError(
                key='created_by_username',
                msg='user: %r does not exist' % initializer_dict['created_by_username'])
        
        user = None
        if 'username' in initializer_dict:
            try:
                user = ScreensaverUser.objects.get(
                    username=initializer_dict['username'])
                initializer_dict['screensaver_user'] = user
            except ObjectDoesNotExist:
                raise Http404('username %r does not exist' 
                    % initializer_dict['username'])
        if 'screensaver_user_id' in initializer_dict:
            try:
                user1 = ScreensaverUser.objects.get(
                    screensaver_user_id=initializer_dict['screensaver_user_id'])
                if user is not None:
                    if user.screensaver_user_id != user1.screensaver_user_id:
                        raise ValidationError({
                            'username': 'does not match the specified screensaver_user_id',
                            'screensaver_user_id': 'does not match the specified username'
                        })
                user = user1
                initializer_dict['screensaver_user'] = user
            except ObjectDoesNotExist:
                raise Http404('screensaver_user_id %r does not exist' 
                    % initializer_dict['screensaver_user_id'])
            
        parent_fields = set(['screensaver_user', 'publication_id', 
            'screen_facility_id', 'reagent_id'])
        initializer_keys = set(initializer_dict.keys())
        if ( parent_fields.isdisjoint(initializer_keys) ):
            msg='must provide one of: %r' % parent_fields
            raise ValidationError({
                'username': msg, 'screen_facility_id': msg,
                'publication_id': msg, 'reagent_id': msg 
            })
        if (len(parent_fields & initializer_keys)>1):
            logger.warn('too many parent_fields: %r', 
                (parent_fields & initializer_keys))
            msg='Only one of parent fields allowed: %r' % parent_fields
            raise ValidationError({
                'username': msg, 'screensaver_user_id': msg, 'screen_facility_id': msg,
                'publication_id': msg, 'reagent_id': msg
            })
        elif 'screen_facility_id' in initializer_dict:
            try:
                screen = Screen.objects.get(
                    facility_id=initializer_dict['screen_facility_id'])
                initializer_dict['screen'] = screen
            except ObjectDoesNotExist:
                raise Http404('screen_facility_id %r does not exist' 
                    % initializer_dict['screen_facility_id'])
        elif 'publication_id' in initializer_dict: 
            try: 
                publication = Publication.objects.get(
                    publication_id=initializer_dict['publication_id'])
                initializer_dict['publication'] = publication
            except ObjectDoesNotExist:
                raise Http404('publication_id does not exist: %s' 
                    % initializer_dict['publication_id'])
        elif 'reagent_id' in initializer_dict: 
            try: 
                reagent = Reagent.objects.get(
                    reagent_id=initializer_dict['reagent_id'])
                initializer_dict['reagent'] = reagent
            except ObjectDoesNotExist:
                raise Http404('reagent_id does not exist: %s' 
                    % initializer_dict['reagent_id'])
        
        af = AttachedFile()
        for key, val in initializer_dict.items():
            if hasattr(af, key):
                setattr(af, key, val)
            else:
                logger.debug('no such attribute on attached_file: %s:%r' 
                    % (key, val))

        parent_log = kwargs.get('parent_log', None)
        if parent_log is None and af.screen is not None:
            parent_log = \
                self.log_to_screen(af.screen, af, request, is_delete=False)
        if parent_log is None and af.screensaver_user is not None:
            parent_log = \
                self.log_to_user(af.screensaver_user, af, request, is_delete=False)
        # TODO: log_to_reagent
        
        af.save()

        # get new state, for logging
        kwargs_for_log = { 
            'attached_file_id': af.attached_file_id, 
            'parent_log': parent_log 
        }
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        log = self.log_patch(
            request, None,new_data, 
            full_create_log=True,**kwargs_for_log)
        if log:
            log.save()

        if not self._meta.always_return_data:
            return HttpResponse(status=200)
        else:
            return self.build_response(request, new_data, status_code=201)
        
    @write_authorization
    @un_cache        
    @transaction.atomic
    def delete_detail(self, request, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        attached_file_id = kwargs.pop('attached_file_id', None)
        if not attached_file_id:
            NotImplementedError('must provide a attached_file_id parameter')
        
        af = AttachedFile.objects.get(attached_file_id=attached_file_id)
        
        try:
            ua_query = UserAgreement.objects.all().filter(file=af)
            if ua_query.exists():
                raise ValidationError(
                    key='%s' % af.filename,
                    msg='File is attached to the current (%s) User Agreement' % ua_query[0].type)
        except ObjectDoesNotExist:
            logger.info('ok to delete attached file, no associated user agreement')
        
        
        _dict = model_to_dict(af)

        # Create parent log
        parent_log = kwargs.get('parent_log', None)
        if parent_log is None and af.screen is not None:
            parent_log = \
                self.log_to_screen(af.screen, af, request, is_delete=True)
        if parent_log is None and af.screensaver_user is not None:
            parent_log = \
                self.log_to_user(af.screensaver_user, af, request, is_delete=True)

        af.delete()
        
        id_attribute = resource = schema['id_attribute']
        log = self.make_log(request, **kwargs)
        log.ref_resource_name = self._meta.resource_name
        log.key = '/'.join([str(_dict[x]) for x in id_attribute])
        log.uri = '/'.join([self._meta.resource_name, log.key])
        log.parent_log = parent_log
        log.api_action = API_ACTION_DELETE
        log.diffs = { k: [v,None] for k,v in _dict.items() }
        log.save()
        logger.debug('delete, api log: %r', log)

        return HttpResponse(status=204)

    @read_authorization
    def get_detail(self, request, **kwargs):

        attached_file_id = kwargs.get('attached_file_id', None)
        if not attached_file_id:
            raise NotImplementedError(
                'must provide a attached_file_id parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        username = param_hash.pop('username', None)
        if username:
            param_hash['username__eq'] = username
        screensaver_user_id = param_hash.pop('screensaver_user_id', None)
        if screensaver_user_id:
            param_hash['screensaver_user_id__eq'] = screensaver_user_id
        screen_facility_id = param_hash.pop('screen_facility_id', None)
        if screen_facility_id:
            param_hash['screen_facility_id__eq'] = screen_facility_id
        attached_file_id = param_hash.pop('attached_file_id', None)
        if attached_file_id:
            param_hash['attached_file_id__eq'] = attached_file_id
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            # add fields for authorization filters
            manual_field_includes.add('screen_facility_id')
            manual_field_includes.add('username')
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)

                  
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
               self._meta.authorization.get_row_property_generator(
                   request.user, field_hash, rowproxy_generator)
 
            # specific setup
            _af = self.bridge['attached_file']
            _su = self.bridge['screensaver_user']
            _screen = self.bridge['screen']
            _publication = self.bridge['publication']
            _reagent = self.bridge['reagent']
            _up = self.bridge['reports_userprofile']
            
            j = _af
            isouter = False
            username = param_hash.pop('username', None)
            if username:
                isouter = True
            screensaver_user_id = param_hash.pop('screensaver_user_id', None)
            if screensaver_user_id:
                isouter = True
            j = j.join(
                _su, _af.c.screensaver_user_id == _su.c.screensaver_user_id,
                isouter=True)
            j = j.join(
                _screen, _af.c.screen_id == _screen.c.screen_id,
                isouter=True)
            j = j.join(
                _publication, _af.c.publication_id == _publication.c.publication_id,
                isouter=True)
            j = j.join(
                _reagent, _af.c.reagent_id == _reagent.c.reagent_id,
                isouter=True)
            
            # This entire query doesn't fit the pattern, construct it manually
            # bleah
            custom_columns = {
                'user_fullname': literal_column(
                    "screensaver_user.last_name || ', ' || screensaver_user.first_name"),
                'created_by_username': literal_column(
                    '(Select au.username '
                    ' from screensaver_user au '
                    ' where au.screensaver_user_id=attached_file.created_by_id )'),
                }

            base_query_tables = [
                'attached_file', 'screensaver_user', 'screen','publication','reagent'] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            
            stmt = select(columns.values()).select_from(j)
            if username:
                stmt = stmt.where(_su.c.username == username)
            if screensaver_user_id:
                stmt = stmt.where(_su.c.screensaver_user_id == screensaver_user_id)
                
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            stmt = stmt.order_by('-attached_file_id')
            
            # compiled_stmt = str(stmt.compile(compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
             
        except Exception, e:
            logger.exception('on get_list %s' % self._meta.resource_name)
            raise e  


class ActivityResourceAuthorization(ScreenAuthorization):        

    def _is_resource_authorized(
        self, user, permission_type, **kwargs):
        authorized = \
            super(ActivityResourceAuthorization, self)\
                ._is_resource_authorized(user, permission_type, **kwargs)
        if authorized is True:
            return True
        
        return user.is_active
    
    def filter_in_sql(self, user, query, screen_table=None, serviced_user_table=None):
        return query
        # # NOTE: replaces filter for performance
        # if self.is_restricted_view(user) is False:
        #     return query
        # 
        # screensaver_user = ScreensaverUser.objects.get(username=user.username)
        # my_screens = self.get_user_screens(screensaver_user)
        # # Can only see own activities
        # 
        # if serviced_user_table is not None and screen_table is not None:
        #     query = query.where(or_(
        #         screen_table.c.screen_id.in_(
        #             [screen.screen_id for screen in my_screens]),
        #         serviced_user_table.c.username == user.username ))
        # elif serviced_user_table is not None:
        #     query = query.where(
        #         serviced_user_table.c.username == user.username )
        # elif screen_table is not None:
        #     query = query.where(
        #         screen_table.c.screen_id.in_(
        #             [screen.screen_id for screen in my_screens]))
        # else:
        #     raise ProgrammingError('must specify screen_table or serviced_user_table')
        # return query
    
    def filter(self, user, filter_expression):
        # TODO testing: use filter_in_sql for performance
        if self.is_restricted_view(user) is False:
            return filter_expression
         
        screensaver_user = ScreensaverUser.objects.get(username=user.username)
        my_screens = self.get_user_screens(screensaver_user)
        # Can only see own activities
         
        if self.resource_name == 'serviceactivity':
            auth_filter = or_(
                column('screen_facility_id').in_(
                    [screen.facility_id for screen in my_screens]),
                column('serviced_user_id') == screensaver_user.screensaver_user_id )
        else:
            auth_filter = column('screen_facility_id').in_(
                    [screen.facility_id for screen in my_screens])
             
        if filter_expression is not None:
            filter_expression = and_(filter_expression, auth_filter)
        else:
            filter_expression = auth_filter
         
        return filter_expression
    
    def get_row_property_generator(self, user, fields, extant_generator):
        # If the user may see the Activity, there are no property restrictions
        return extant_generator
    
    def has_activity_read_authorization(self, user, activity_id):

        if self.is_restricted_view(user) is False:
            return True
        else:
            screensaver_user = ScreensaverUser.objects.get(username=user.username)
            my_screens = self.get_user_screens(screensaver_user)
            
            activity = Activity.objects.get(activity_id=activity_id)
            if hasattr(activity, 'labactivity'):
                return activity.labactivity.screen.facility_id in set([
                    screen.facility_id for screen in my_screens])
            if hasattr(activity, 'serviceactivity'):
                return activity.serviceactivity.serviced_user.username == user.username
            
            return False        
        
class ActivityResource(DbApiResource):
    '''
    Activity Resource is a combination of the LabActivity and the ServiceActivity
    
    NOTE: 20170523
    ActivityResource needs to be reworked; the current design is from the 
    legacy Screensaver 1; 
    Activity is both a logging facility and a join class between User, Screen, 
    and CherryPickRequest.
    Refactor to make Activity a Reporting resource that logs:
    - performed by
    - serviced resource: Screen, User
    - serviced context: Screen, User, CherryPickRequest
    - referenced resources: Equipment
    - 
    '''

    class Meta:

        queryset = Activity.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'activity'
        authorization = ActivityResourceAuthorization(resource_name)
        serializer = LimsSerializer()
        ordering = []
        filtering = {}
        always_return_data = True 
        
    def __init__(self, **kwargs):

        self.service_activity_resource = None
        self.screen_resource = None
        self.su_resource = None
        super(ActivityResource, self).__init__(**kwargs)

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<activity_id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
    
    def get_su_resource(self):
        if self.su_resource is None:
            self.su_resource = ScreensaverUserResource()
        return self.su_resource

    def get_service_activity_resource(self):
        if self.service_activity_resource is None:
            self.service_activity_resource = ServiceActivityResource()
        return self.service_activity_resource

    def get_screen_resource(self):
        if self.screen_resource is None:
            self.screen_resource = ScreenResource()
        return self.screen_resource

    def build_schema(self, user=None):
         
        schema = super(ActivityResource, self).build_schema(user=user)
        return schema

    @read_authorization
    def get_detail(self, request, **kwargs):
 
        activity_id = kwargs.pop('activity_id', None)
        if not activity_id:
            raise Http404('must provide an activity_id parameter')
        else:
            kwargs['activity_id__eq'] = activity_id
 
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def get_custom_columns(self, alias_qualifier):
        '''
        Convenience method for subclasses: reusable custom columns
        @param alias_qualifier a sql compatible string used to name subqueries
            so that this method may be called multiple times to compose a query
        '''
        _screen = self.bridge['screen']
        _activity = self.bridge['activity']
        _su = self.bridge['screensaver_user']
        # perform some hacks to speed up the query
        _performed_by_users_cte = \
            ScreensaverUserResource.get_user_cte()\
            .where(
                exists(select([None]).select_from(_activity)
                    .where(_su.c.screensaver_user_id==_activity.c.performed_by_id)))\
            .cte('performers_%s' % alias_qualifier)
        _performed_by = _performed_by_users_cte.alias('performed_by_%s' % alias_qualifier)
        _performed_by1 = _performed_by_users_cte.alias('performed_by1_%s' % alias_qualifier)
        _performed_by2 = _performed_by_users_cte.alias('performed_by2_%s' % alias_qualifier)
        _created_by_users_cte = \
            ScreensaverUserResource.get_user_cte()\
            .where(
                exists(select([None]).select_from(_activity)
                    .where(_su.c.screensaver_user_id==_activity.c.created_by_id)))\
            .cte('creators_%s' % alias_qualifier)
        _created_by = _created_by_users_cte.alias('created_by_%s' % alias_qualifier)
        _created_by1 = _created_by_users_cte.alias('created_by1_%s' % alias_qualifier)
        _lhsu = _su.alias('lhsu_%s' % alias_qualifier)
        _sfs = self.bridge['screen_funding_supports']
        
        lab_head_table = \
            ScreensaverUserResource.get_lab_head_cte(alias_qualifier)\
                .cte('lab_heads_%s' % alias_qualifier)
        return {
            'performed_by_name': (
                select([_performed_by1.c.name])
                    .select_from(_performed_by1)
                    .where(_performed_by1.c.screensaver_user_id 
                        == _activity.c.performed_by_id)
                ),
            'performed_by_username': (
                select([_performed_by.c.username])
                    .select_from(_performed_by)
                    .where(_performed_by.c.screensaver_user_id 
                        == _activity.c.performed_by_id)
                ),
            'performed_by_user_id': (
                select([_performed_by2.c.screensaver_user_id])
                    .select_from(_performed_by2)
                    .where(_performed_by2.c.screensaver_user_id 
                        == _activity.c.performed_by_id)
                ),
            'created_by_name': (
                select([_created_by1.c.name])
                    .select_from(_created_by1)
                    .where(_created_by1.c.screensaver_user_id 
                        == _activity.c.created_by_id)
                ),
            'created_by_username': (
                select([_created_by.c.username])
                    .select_from(_created_by)
                    .where(_created_by.c.screensaver_user_id 
                        == _activity.c.created_by_id)
                ),
            'screen_lab_affiliation': (
                select([lab_head_table.c.lab_affiliation])
                .select_from(lab_head_table)
                .where(lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id)),
            'screen_lab_name': (
                select([lab_head_table.c.lab_name_full])
                .select_from(lab_head_table)
                .where(lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id)),
            'screen_lab_head_id': (
                select([lab_head_table.c.screensaver_user_id])
                .select_from(lab_head_table)
                .where(lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id)),
            'screen_funding_supports':
                select([func.array_to_string(
                    func.array_agg(literal_column('funding_support')
                    ), LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([_sfs.c.funding_support])
                            .select_from(_sfs)
                            .order_by(_sfs.c.funding_support)
                            .where(_sfs.c.screen_id 
                                == literal_column('screen.screen_id'))
                            .alias('inner')),
            'screen_lead_screener_name': (
                select([_concat(_su.c.first_name, ' ', _su.c.last_name)])
                    .select_from(_su)
                    .where(_su.c.screensaver_user_id 
                        == _screen.c.lead_screener_id)),
            'screen_lead_screener_id': (
                select([_su.c.screensaver_user_id])
                    .select_from(_su)
                    .where(_su.c.screensaver_user_id 
                        == _screen.c.lead_screener_id)),
            'screen_date_of_last_activity': literal_column(
                '( (select date_of_activity '
                '  from activity '
                '  join lab_activity la using(activity_id) '
                '  where la.screen_id=screen.screen_id '
                '  and not exists(select null from cherry_pick_liquid_transfer cplt'
                '     where cplt.activity_id = activity.activity_id) '
                '  UNION ALL'
                '  select date_of_activity '
                '  from activity '
                '  join service_activity sa using(activity_id) '
                '  where sa.serviced_screen_id=screen.screen_id )'
                '  order by date_of_activity desc LIMIT 1 )'),
        }

    def get_query(self, param_hash, user):
        
        # general setup
        schema = self.build_schema(user=user)
        
        manual_field_includes = set(param_hash.get('includes', []))
        # for join to screen query (TODO: only include if screen fields rqst'd)
        manual_field_includes.add('screen_id')
        manual_field_includes.add('activity_class')
        # for join to cherrypickrequest (TODO: req'd to complete activity id link)
        manual_field_includes.add('cherry_pick_request_id')
        param_hash['includes'] = list(manual_field_includes)
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, schema)
        
        # NOTE: try "filter_in_sql" for performance
        # NOTE: filters are done in the subquery clauses
#         filter_expression = \
#             self._meta.authorization.filter(user,filter_expression)
              
        order_params = param_hash.get('order_by', [])
        order_params.append('-date_of_activity')
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)

        # specific setup
        
        _activity = self.bridge['activity']
        # Create a UNION query of each subclass query:
        (field_hash_sa, columns_sa, stmt_sa, count_stmt_sa, filenamesa) = (
            self.get_service_activity_resource().get_query(param_hash, user))
        # if a field is not present in the subquery, create an empty field            
        sa_columns = []
        for key, field in field_hash.items():
            if field['scope'] == 'fields.activity':
                if key not in columns_sa:
                    if field['data_type'] == 'string':
                        sa_columns.append(cast(literal_column("null"), 
                            sqlalchemy.sql.sqltypes.Text).label(key))
                    else:
                        sa_columns.append(literal_column("null").label(key))
                else:
                    sa_columns.append(literal_column(key))
#         logger.info('sa_columns: %r', sa_columns)
        stmt_sa = stmt_sa.cte('serviceactivities')
        stmt1 = select(sa_columns).select_from(stmt_sa)

        (field_hash_la, columns_la, stmt_la, count_stmt_la) = (
            self.get_lab_activity_query(param_hash, user))
        # if a field is not present in the subquery, create an empty field            
        la_columns = []
        for key, field in field_hash.items():
            if field['scope'] == 'fields.activity':
                if key not in columns_la:
                    logger.info('create dummy col for %r', key)
                    if field['data_type'] == 'string':
                        la_columns.append(cast(literal_column("null"), 
                            sqlalchemy.sql.sqltypes.Text).label(key))
                    else: 
                        la_columns.append(literal_column("null").label(key))
                else:
                    la_columns.append(literal_column(key))
                    
#         logger.info('la_columns: %r', la_columns)
        stmt_la = stmt_la.cte('labactivities')
        stmt2 = select(la_columns).select_from(stmt_la)
        
        stmt = stmt1.union_all(stmt2)
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)

        if logger.isEnabledFor(logging.DEBUG):
            compiled_stmt = str(stmt.compile(
                dialect=postgresql.dialect(),
                compile_kwargs={"literal_binds": True}))
            logger.info('compiled_stmt %s', compiled_stmt)


        columns = { key:literal_column(key) for key in field_hash.keys()}
        return (field_hash, columns, stmt, count_stmt, filename)
    
    def get_custom_lab_activity_columns(self, alias_qualifier):

        _library_screening = self.bridge['library_screening']
        _cps = self.bridge['cherry_pick_screening']
        # _cplt = self.bridge['cherry_pick_liquid_transfer']
        _screen = self.bridge['screen']
        
        activity_type_column = cast(case([
            (_library_screening.c.activity_id != None,
               case([(_library_screening.c.is_for_external_library_plates,
                        'externallibraryscreening')],
                    else_='libraryscreening')),
            (_cps.c.activity_id != None,
                'cherrypickscreening')
            ],
            else_='cplt'), sqlalchemy.sql.sqltypes.Text)

        return { 
            'serviced_user_id': cast(
                literal_column("null"),sqlalchemy.sql.sqltypes.Numeric),
            'serviced_user': cast(
                literal_column("null"),sqlalchemy.sql.sqltypes.Text),
            'serviced_username': cast(
                literal_column("null"),sqlalchemy.sql.sqltypes.Text),
            'funding_support': cast(
                literal_column("null"),sqlalchemy.sql.sqltypes.Text),
            'cherry_pick_request_id': _cps.c.cherry_pick_request_id,
            # 20171108 - Do not show CPLT activities in the general activity report
            # 'cherry_pick_request_id': 
            #     cast(func.coalesce(
            #             _cps.c.cherry_pick_request_id, 
            #             _cplt.c.cherry_pick_request_id, None), 
            #         sqlalchemy.sql.sqltypes.Text),
            'type': activity_type_column,
            'activity_class': activity_type_column,
            }
        
    def get_lab_activity_query(self, param_hash, user):

        # general setup
        # schema for labactivity part of the query should only be the activity fields
        # (exclude the screen.fields)
        schema = deepcopy(self.build_schema(user=user))
        field_hash = schema['fields']
        field_hash = { key:val for key, val in field_hash.items() 
            if val['scope'] == 'fields.activity'}  
        schema['fields'] = field_hash
        
        manual_field_includes = set(param_hash.get('includes', []))
        manual_field_includes.add('screen_id')
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filter_expression = \
            self._meta.authorization.filter(user,filter_expression)
              
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
         
        # specific setup
        _sfs = self.bridge['screen_funding_supports']
        _a = self.bridge['activity']
        _la = self.bridge['lab_activity']
        _screening = self.bridge['screening']
        _screen = self.bridge['screen']

        j = _a
        j = j.join(_la, _a.c.activity_id == _la.c.activity_id)
        j = j.join(_screen, _la.c.screen_id == _screen.c.screen_id)

        # TODO: delegate to sub_classes (when built)
        _library_screening = self.bridge['library_screening']
        _cps = self.bridge['cherry_pick_screening']
        # 20171108 - Do not show CPLT activities in general report
        # _cplt = self.bridge['cherry_pick_liquid_transfer']
        j = j.join(
            _library_screening,
            _la.c.activity_id == _library_screening.c.activity_id, isouter=True)
        j = j.join(_cps, _la.c.activity_id == _cps.c.activity_id, isouter=True)
        # j = j.join(_cplt, _la.c.activity_id == _cplt.c.activity_id, isouter=True)
                
        custom_columns = self.get_custom_columns('la')
        custom_columns.update(self.get_custom_lab_activity_columns('lab_activity'))
        
        base_query_tables = ['activity', 'lab_activity', 'screening', 'screen'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        stmt = select(columns.values()).select_from(j)
        
        # stmt = self._meta.authorization.filter_in_sql(
        #     user, stmt, screen_table=_screen)
        # general setup
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        
        if logger.isEnabledFor(logging.DEBUG):
            compiled_stmt = str(stmt.compile(
                dialect=postgresql.dialect(),
                compile_kwargs={"literal_binds": True}))
            logger.info('compiled_stmt %s', compiled_stmt)
        
        return (field_hash, columns, stmt, count_stmt)

    def build_list_response(self, request, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        try:
            (field_hash, columns, stmt, count_stmt, filename) = \
                self.get_query(param_hash, request.user)
            
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
            rowproxy_generator = \
               self._meta.authorization.get_row_property_generator(
                   request.user, field_hash, rowproxy_generator)

            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

class ServiceActivityResource(ActivityResource):    

    class Meta:

        queryset = ServiceActivity.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'serviceactivity'
        authorization = ActivityResourceAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        super(ServiceActivityResource, self).__init__(**kwargs)

    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<activity_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/for_user/(?P<serviced_user_id>([\d]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/for_user/(?P<serviced_user_id>([\d]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]    

    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs.pop('schema', None)
        if schema is None:
            raise Exception('schema not initialized')
        fields = schema['fields']

        # NOTE: parse params only if needed for client overrides
        # param_hash = self._convert_request_to_dict(request)
        # param_hash.update(kwargs)
        # logger.debug('param_hash: %r', param_hash)

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)

        logger.info('patch ServiceActivity: %r', deserialized)
        
        patch = bool(id_kwargs)
        initializer_dict = self.parse(deserialized, create=not patch, schema=schema)
        errors = self.validate(initializer_dict, schema=schema, patch=patch)
        if errors:
            raise ValidationError(errors)
        
        activity_type = deserialized.get('type', None)
        if activity_type:
            initializer_dict['service_activity_type'] = activity_type

        serviced_user_id = deserialized.get('serviced_user_id', None)
        performed_by_user_id = deserialized.get('performed_by_user_id', None)
        serviced_screen_facility_id = deserialized.get('screen_facility_id', None)
        
        if not patch:
            if serviced_user_id is None and serviced_screen_facility_id is None:
                msg = 'Either serviced user or serviced screen_facility_id is required'
                raise ValidationError({
                    'serviced_user_id': msg,
                    'screen_facility_id': msg
                    })
            if not activity_type:
                raise ValidationError(
                    key='type',
                    msg='required')
            if not performed_by_user_id:
                raise ValidationError(
                    key='performed_by_user_id',
                    msg='required')

        if serviced_user_id:
            try:
                serviced_user = ScreensaverUser.objects.get(
                    screensaver_user_id=serviced_user_id)
                initializer_dict['serviced_user'] = serviced_user
            except ObjectDoesNotExist:
                raise ValidationError(
                    key='serviced_user_id',
                    msg='id does not exist: %r' % serviced_user_id)
        if serviced_screen_facility_id:
            try:
                serviced_screen = Screen.objects.get(
                    facility_id=serviced_screen_facility_id)
                initializer_dict['serviced_screen'] = serviced_screen
            except ObjectDoesNotExist:
                raise ValidationError(
                    key='screen_facility_id',
                    msg='screen does not exist: %r' % serviced_screen_facility_id)
        if performed_by_user_id:
            performed_by = \
                self.get_su_resource()._get_detail_response_internal(
                exact_fields=['screensaver_user_id','is_staff'],
                screensaver_user_id=performed_by_user_id)
            if not performed_by: 
                raise ValidationError(
                    key='performed_by_user_id',
                    msg='No such screensaver_user_id: %r' % performed_by_user_id)
            if performed_by.get('is_staff',False) != True:
                raise ValidationError(
                    key='performed_by_user_id',
                    msg='Must be a staff user')

        service_activity = None
        if patch:
            try:
                service_activity = ServiceActivity.objects.get(
                    pk=id_kwargs['activity_id'])
            except ObjectDoesNotExist:
                raise Http404(
                    'ServiceActivity does not exist for: %r', id_kwargs)
        else:
            # Set the created_by field:
            # NOTE: deprecate for SS V2
            try:
                adminuser = ScreensaverUser.objects.get(username=request.user.username)
            except ObjectDoesNotExist as e:
                logger.error('admin user: %r does not exist', request.user.username )
                raise
            
            service_activity = ServiceActivity()
            service_activity.created_by = adminuser
        
        service_activity.performed_by_id = performed_by_user_id
        
        model_field_names = [
            x.name for x in service_activity._meta.get_fields()]
        for key, val in initializer_dict.items():
            if key in model_field_names:
                setattr(service_activity, key, val)

        service_activity.save()
        logger.info('saved service_activity: %r', service_activity)
        return { API_RESULT_OBJ: service_activity }

    def make_log_key(self, log, attributes, id_attribute=None, schema=None, **kwargs):

        logger.debug('make_log_key: %r, %r, %r', attributes, id_attribute, kwargs)
        
        if attributes:
            ActivityResource.make_log_key(
                self, log, attributes, id_attribute=id_attribute, 
                schema=schema, **kwargs)
            logger.info('log key: %r, %r', log.key, log)
    
            keys = []
            # Always create the service activity log for the user, if available 
            if attributes.get('serviced_user_id', None):
                keys.append('screensaveruser')
                keys.append(str(attributes['serviced_user_id']))
            elif attributes.get('screen_facility_id', None):
                keys.append('screen')
                keys.append(attributes['screen_facility_id'])
            keys.append(self._meta.resource_name)
            log.uri = '%s/%s' % ('/'.join(keys), log.key)
        
            logger.info('log uri: %r, %r', log.uri, log)

    @write_authorization
    @un_cache
    @transaction.atomic
    def delete_obj(self, request, deserialized, **kwargs):

        activity_id = kwargs.get('activity_id', None)
        if activity_id:
            try:
                sa = ServiceActivity.objects.get(
                    activity_id=activity_id)
                sa.delete()
            except ObjectDoesNotExist:
                logger.warn('no such ServiceActivity: %s' % activity_id)
                raise Exception(
                    'ServiceActivity for activity_id: %s not found' % activity_id)
        else:
            raise Exception(
                'ServiceActivity delete action requires an activity_id %s' 
                % kwargs)
        
    def get_query(self, param_hash, user):

        schema = self.build_schema(user=user)
        logger.debug('serviceactivity fields: %r', schema['fields'].keys())
        
        # general setup
        alias_qualifier = 'sa'
        manual_field_includes = set(param_hash.get('includes', []))
        # for join to screen query (TODO: only include if screen fields rqst'd)
        manual_field_includes.add('screen_id')
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, schema)
        filter_expression = \
            self._meta.authorization.filter(user,filter_expression)

              
        order_params = param_hash.get('order_by', [])
        order_params.append('-date_of_activity')
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        field_hash = { key:val for key, val in field_hash.items() 
            if val['scope'] == 'fields.serviceactivity'}  
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
         
        # specific setup
        _a = self.bridge['activity']
        _sa = self.bridge['service_activity']
        _screen = self.bridge['screen']
        _user_cte = ScreensaverUserResource.get_user_cte().cte('users_serviced_sa')
        _serviced = _user_cte.alias('serviced_user')
        
        j = _a
        j = j.join(_sa, _a.c.activity_id == _sa.c.activity_id)
        j = j.join(
            _serviced,
            _sa.c.serviced_user_id == _serviced.c.screensaver_user_id, isouter=True)
        j = j.join(
            _screen,
            _sa.c.serviced_screen_id == _screen.c.screen_id, isouter=True)
        
        # Get custom columns from the parent (ActivityResource); the
        # ServiceActivity query rows contain the needed activity and screen ids
        custom_columns = \
            super(ServiceActivityResource, self).get_custom_columns('sa')
        custom_columns.update({
            'activity_class': cast(
                literal_column("'serviceactivity'"), 
                sqlalchemy.sql.sqltypes.Text),
            'serviced_user': _serviced.c.name,
            'serviced_username': _serviced.c.username,
            'serviced_user_id': _serviced.c.screensaver_user_id,
            })

        base_query_tables = ['activity', 'service_activity', 'screen'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        stmt = select(columns.values()).select_from(j)
        # stmt = self._meta.authorization.filter_in_sql(
        #     user, stmt, serviced_user_table=_serviced, screen_table=_screen)
        # general setup
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)

        # compiled_stmt = str(stmt.compile(
        #     dialect=postgresql.dialect(),
        #     compile_kwargs={"literal_binds": True}))
        # logger.info('compiled_stmt %s', compiled_stmt)
        
        return (field_hash, columns, stmt, count_stmt, filename)


class LibraryScreeningResource(ActivityResource):
    
    # Constants: may be overridden in the settings.py
    MIN_WELL_VOL_SMALL_MOLECULE = Decimal('0.0000069') # 6.9 uL
    MIN_WELL_VOL_RNAI = Decimal(0)
    
    ALLOWED_LIBRARY_SCREENING_STATUS = ('allowed',)
    WARN_LIBRARY_SCREENING_STATUS = (
        'requires_permission','not_recommended','retired',)
    # NOTE: all other library screening status are error statuses
    # ERROR_LIBRARY_SCREENING_STATUS = (
    #     'not_allowed','discarded',)
    ALLOWED_PLATE_STATUS = ('available', 'retired',)
    WARN_PLATE_STATUS = ('retired',)
    # NOTE: all other status are error statuses:
    # ERROR_PLATE_STATUS = [
    #     'discarded', 'given_away','not_specified', 'not_available', 
    #     'not_created', 'discarded_volume_transferred', 'lost']
    ALLOWED_COPY_USAGE_TYPE = ('library_screening_plates',)
    # NOTE: all other copy usage types are errors
    SCREENING_COUNT_THRESHOLD = 12

    # Messages for Screening Inquiry errors
    MSG_ALREADY_SCREENED = 'Plates have already been screened'
    MSG_PLATE_STATUS_ERROR = 'Plate status'
    MSG_SCREENING_COUNT = 'Screening count > %d' % SCREENING_COUNT_THRESHOLD
    MSG_PLATES_WELLS_ADJUSTED = 'Plate well volumes have been adjusted'
    MSG_COPY_USAGE_TYPE = 'Copy usage type'
    MSG_INSUFFICIENT_VOL = 'Insufficient vol: %s uL'
    MSG_NO_PLATE_VOLUME = 'No plate volume recorded'
    MSG_LIBRARY_SCREENING_STATUS = 'Library screening status'
    MSG_LIBRARY_SCREENING_TYPE = 'Library screening type'

    
    class Meta:

        queryset = LibraryScreening.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'libraryscreening'
        alt_resource_name = 'externallibraryscreening'
        authorization = ActivityResourceAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        super(LibraryScreeningResource, self).__init__(**kwargs)
        
        self.plate_resource = None
        self.copywell_resource = None
        self.library_resource = None
        self.lcp_resource = None
        
    def get_librarycopyplate_resource(self):
        if self.lcp_resource is None:
            self.lcp_resource = LibraryCopyPlateResource()
        return self.lcp_resource
    
    def get_library_resource(self):
        if self.library_resource is None:
            self.library_resource = LibraryResource()
        return self.library_resource
        
    def get_plate_resource(self):
        if self.plate_resource is None:
            self.plate_resource = LibraryCopyPlateResource()
        return self.plate_resource

    def get_copywell_resource(self):
        if self.copywell_resource is None:
            self.copywell_resource = CopyWellResource()
        return self.copywell_resource
    
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<activity_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<activity_id>([\d]+))/plates%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_plates_screened_view'),
                name="api_dispatch_plates_screened_view"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<activity_id>([\d]+))/libraries%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_libraries_screened_view'),
                name="api_dispatch_libraries_screened_view"),
        ]  
        
    def dispatch_plates_screened_view(self, request, **kwargs):

        library_screening_id = kwargs.pop('activity_id')
        # NOTE: authorization is performed in LibraryCopyPlateResource
        return self.get_plate_resource()\
            .build_screened_plate_response(
                request, library_screening_id=library_screening_id, **kwargs)    
        
    def dispatch_libraries_screened_view(self, request, **kwargs):

        library_screening_id = kwargs.pop('activity_id')
        # add custom authorization
        if self._meta.authorization.has_activity_read_authorization(
                request.user, library_screening_id) is False:
            raise PermissionDenied

        with get_engine().connect() as conn:
            _l = self.bridge['library']
            _ap = self.bridge['assay_plate']
            _c = self.bridge['copy']
            _p = self.bridge['plate']
            
            library_names = [x[0] for x in 
                conn.execute(
                    select([distinct(_l.c.short_name)])
                    .select_from(
                        _l.join(_c, _l.c.library_id==_c.c.library_id)
                          .join(_p, _p.c.copy_id==_c.c.copy_id)
                          .join(_ap, _ap.c.plate_id==_p.c.plate_id))
                    .where(_ap.c.library_screening_id==library_screening_id))]
            
            return self.get_library_resource().get_list(request,
                short_name__in=library_names, **kwargs)

    def dispatch_plate_range_search_view(self, request, **kwargs):
        ''' 
        Find: 
        - plates already asssociated with the libraryscreenings for the screen
        - plates matched by the "raw_search_data"
        Note: bypasses the "dispatch" framework call
        -- must be authenticated and authorized
        '''
        logger.info('dispatch_plate_range_search_view')
        self.is_authenticated(request)
        resource_name = kwargs.pop('resource_name', self._meta.resource_name)
        if not self._meta.authorization._is_resource_authorized(
                request.user,'read'):
            raise PermissionDenied(
                'user: %s, permission: %s/%s not found' 
                    % (request.user,resource_name,'read'))
        if request.method.lower() not in ['get','post']:
            return self.dispatch('detail', request, **kwargs)

        # With POST, params may come in the request
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        logger.info('params: %r', param_hash)
        
        facility_id = param_hash.get('facility_id', None)
        activity_id = param_hash.get('activity_id', None)
        plate_search_data = param_hash.get('raw_search_data', None)
        volume_required = parse_val(param_hash.get('volume_required', None),
            'volume_required', 'decimal')
        show_retired_plates = parse_val(
            param_hash.get('show_retired_plates'),
            'show_retired_plates','boolean')
        plate_status_types = ['available']
        if show_retired_plates:
            plate_status_types.append('retired')
        show_first_copy_only = parse_val(
            param_hash.get('show_first_copy_only',None),
            'show_first_copy_only','boolean')    
        hide_existing = parse_val(
            param_hash.get('hide_existing',None),
            'hide_existing','boolean')    
        
        if not facility_id:
            raise NotImplementedError(
                'must provide facility_id')
        if self.get_screen_resource()._meta.authorization\
            .has_screen_read_authorization(request.user, facility_id) is False:
            raise PermissionDenied
        try:
            screen = Screen.objects.get(facility_id=facility_id)
        except ObjectDoesNotExist:
            raise Http404(
                'screen does not exist for: %r', facility_id)
        filename = 'plate_search_for_%s' % screen.facility_id

        library_screening = None
        if activity_id is not None:
            filename += '_screening_' + activity_id
            if self._meta.authorization.has_activity_read_authorization(
                request.user, activity_id) is False:
                raise PermissionDenied
            library_screening = LibraryScreening.objects.get(activity_id=activity_id)
        
        searched_plate_ids = []
        plate_search_errors = []
        if plate_search_data:
            (plates,parsed_searches, errors) = \
                self.get_librarycopyplate_resource().find_plates(plate_search_data)
            if errors:
                logger.info('errors: %r', errors)
                plate_search_errors.append(
                    'Plates not found: %s' % ', '.join(errors))
            # library screening plates only
            library_screening_plates_only = []
            wrong_type_errors = defaultdict(set)
            for plate in plates:
                if plate.copy.usage_type in self.ALLOWED_COPY_USAGE_TYPE:
                    library_screening_plates_only.append(plate)
                else:
                    # Only show error if explicitly searching for the copy-plate
                    copy_plate_key = '%s/%s' %(plate.copy.name, plate.plate_number)
                    for parsed_search in parsed_searches:
                        if copy_plate_key in parsed_search['plate_copy_keys_expected']:
                            wrong_type_errors[plate.copy.usage_type].add(copy_plate_key)
            plates = library_screening_plates_only
            if wrong_type_errors:
                plate_search_errors.append(
                    self.MSG_COPY_USAGE_TYPE + ': %s' % '; '.join(
                        ['%s: %s' % (k, ', '.join(v)) 
                            for k,v in wrong_type_errors.items()]))
            # status types
            allowed_status_plates_only = []
            wrong_status_errors = defaultdict(set)
            for plate in plates:
                if plate.status in plate_status_types:
                    allowed_status_plates_only.append(plate)
                else:
                    # Only show error if explicitly searching for the copy-plate
                    copy_plate_key = '%s/%s' %(plate.copy.name, plate.plate_number)
                    for parsed_search in parsed_searches:
                        if copy_plate_key in parsed_search['plate_copy_keys_expected']:
                            wrong_status_errors[plate.status].add(copy_plate_key)
            plates = allowed_status_plates_only
            if wrong_status_errors:
                plate_search_errors.append(
                    self.MSG_PLATE_STATUS_ERROR + ': %s' % '; '.join(
                        ['%s: %s' % (k, ', '.join(v)) 
                            for k,v in wrong_status_errors.items()]))
                        
                        
            logger.info('show_first_copy_only: %r',show_first_copy_only)
            if show_first_copy_only is not True:
                searched_plate_ids = [x.plate_id for x in plates]
            else:
                plate_copy_map = defaultdict(list)
                for plate in plates:
                    plate_copy_map[plate.plate_number].append(
                        (plate.plate_location.freezer, plate.copy.name,plate.plate_id))
                for plate_number,copy_list in plate_copy_map.items():
                    searched_plate_ids.append(sorted(copy_list)[0][2])
        _data = []
        _data = self.get_plate_range_search_table(
            screen, searched_plate_ids,
            volume_required=volume_required)
        
        # If the library_screening ID is provided, filter out the extra rows
        # for other screenings as a convenience
        if library_screening is not None:
            new_data = []
            for x in _data:
                logger.info('filter data: %r', x)
                if x['library_screening_id'] == 0:
                    new_data.append(x)
                elif library_screening is not None:
                    if x['library_screening_id'] == library_screening.activity_id:
                        new_data.append(x)
            _data = new_data
        if hide_existing is True:
            logger.info('hide_existing for plate range search...')
            new_data = []
            for x in _data:
                if x['library_screening_id'] == 0:
                    new_data.append(x)
            _data = new_data
        
        meta = { 'total_count': len(_data) }
        if plate_search_errors:
            meta[API_MSG_WARNING] = plate_search_errors
        response_data = {
            API_RESULT_META: meta,
            API_RESULT_DATA: _data 
        }
        return self.build_response(request, response_data, filename=filename)

    def get_plate_range_search_table(
        self, screen, searched_plate_ids, volume_required=None):
        '''
        Generate a (screening inquiry) plate range table:
        - for the searched_plate_ids
        - find extant plate_ranges for the screen 
        - check for errors and warnings 
        '''    

        LSR = self
        logger.info('get_plate_range_search_table for ids: %r ...', searched_plate_ids)
        class ErrorDict():
            ''' Track errors and warnings for each plate range
            '''
            def __init__(self):
                self.errors = defaultdict(set)
                self.warnings = defaultdict(set)
                self.plate_errors = defaultdict(set)
                self.plate_warnings = defaultdict(set)
                
            def addPlateError(self, key, plate_number):
                self.plate_errors[key].add(plate_number)
            def addPlateWarning(self, key, plate_number):
                self.plate_warnings[key].add(plate_number)
            def addError(self, key, value):
                self.errors[key].add(value)
            def addWarning(self, key, value):
                self.warnings[key].add(value)
            
            def check_row(self,_row):
                logger.debug('check row: %r', _row)
                if  _row['activity_id'] == 0:
                    if _row['plate_number'] in extant_plate_numbers:
                        self.addPlateWarning(
                            LSR.MSG_ALREADY_SCREENED, plate_number)
                    if _row['status'] not in LSR.ALLOWED_PLATE_STATUS:
                        self.addPlateError(
                            LSR.MSG_PLATE_STATUS_ERROR + ': %s' % _row['status'],
                            plate_number)
                    if _row['status'] in LSR.WARN_PLATE_STATUS:
                        self.addPlateWarning(
                            LSR.MSG_PLATE_STATUS_ERROR + ': %s' % _row['status'],
                            plate_number)
                    if  _row['usage_type'] != 'library_screening_plates':
                        self.addError(LSR.MSG_COPY_USAGE_TYPE, _row['usage_type'])
                    if _row['remaining_well_volume']:
                        vol_min = LSR.MIN_WELL_VOL_RNAI
                        error_key = LSR.MSG_INSUFFICIENT_VOL
                        if screen.screen_type == VOCAB_SCREEN_TYPE_SM:
                            vol_min = LSR.MIN_WELL_VOL_SMALL_MOLECULE
                            error_key += (' (req %s)' 
                                % lims_utils.convert_decimal(vol_min,1e-6, 1))
                        if volume_required is not None:
                            vol_after_transfer = (
                                Decimal(_row['remaining_well_volume']) 
                                    - volume_required )
                            if vol_after_transfer < vol_min:
                                self.addPlateError(
                                    error_key % lims_utils.convert_decimal(
                                        vol_after_transfer,1e-6, 1),
                                    plate_number)
                    else:
                        self.addPlateWarning(LSR.MSG_NO_PLATE_VOLUME,plate_number)
                    if _row.get('cplt_screening_count',0) > 0:
                        self.addPlateWarning(
                            LSR.MSG_PLATES_WELLS_ADJUSTED, plate_number)
                    if _row.get('screening_count',0) > LSR.SCREENING_COUNT_THRESHOLD:
                        self.addPlateWarning(
                            LSR.MSG_SCREENING_COUNT, plate_number)

                if (_row['library_screening_status'] 
                        in LSR.WARN_LIBRARY_SCREENING_STATUS):
                    self.addWarning(
                        LSR.MSG_LIBRARY_SCREENING_STATUS,
                        _row['library_screening_status'])
                elif (_row['library_screening_status'] 
                        not in LSR.ALLOWED_LIBRARY_SCREENING_STATUS):
                    self.addError(
                        LSR.MSG_LIBRARY_SCREENING_STATUS,
                        _row['library_screening_status'])
                if _row['library_screen_type'] != screen.screen_type:
                    self.addError(LSR.MSG_LIBRARY_SCREENING_TYPE,
                        _row['library_screen_type'])
            
            def showErrors(self):
                full_errors = [ '%s: %s' % (k,', '.join(v))
                    for k,v in self.errors.items()]
                for k,v in self.plate_errors.items():
                    full_errors.append(
                        '%s: %s' % (k, ', '.join(lims_utils.find_ranges(v))))
                return sorted(full_errors)

            def showWarnings(self):
                full_warnings = [ '%s: %s' % (k,', '.join(v))
                    for k,v in self.warnings.items()]
                for k,v in self.plate_warnings.items():
                    full_warnings.append(
                        '%s: %s' % (k, ', '.join(lims_utils.find_ranges(v))))
                return sorted(full_warnings)

        logger.info('plate range search screen: %r, vol: %r', 
            screen.facility_id, volume_required)
        with get_engine().connect() as conn:
            _a = self.bridge['activity']
            _c = self.bridge['copy']
            _p = self.bridge['plate']
            _pl = self.bridge['plate_location']
            _l = self.bridge['library']
            _ls = self.bridge['library_screening']
            _la = self.bridge['lab_activity']
            _ap = self.bridge['assay_plate']
            _screen = self.bridge['screen']
            
            fields = [
                'activity_id', 'plate_key', 'library_short_name', 
                'library_screening_status', 'library_screen_type', 'copy_name', 
                'copy_comments', 'plate_location', 'usage_type', 
                'plate_number','cplt_screening_count', 'screening_count', 
                'remaining_well_volume','status']
            # 1. query for current library screening plates
            _assay_plates_query = (select([
                _ap.c.plate_id,
                _ap.c.plate_number,
                _ls.c.activity_id,
                _a.c.date_of_activity,
                _c.c.copy_id,
                _c.c.name.label('copy_name'),
                _l.c.short_name,
                _c.c.library_id,
                ])
                .select_from(
                    _ap.join(_ls,_ap.c.library_screening_id==_ls.c.activity_id)
                       .join(_a,_ls.c.activity_id==_a.c.activity_id)
                       .join(_p, _ap.c.plate_id==_p.c.plate_id)
                       .join(_c, _c.c.copy_id==_p.c.copy_id)
                       .join(_l, _c.c.library_id==_l.c.library_id)
                       .join(_la,_ls.c.activity_id==_la.c.activity_id)
                    )
                .where(_ap.c.replicate_ordinal==0))
            extant_plate_numbers = []
            if screen is not None:
                _assay_plates_query = _assay_plates_query.where(
                    _la.c.screen_id==screen.screen_id)
                extant_plate_numbers = [
                    x[1] for x in conn.execute(_assay_plates_query)]
            _assay_plates_query = _assay_plates_query.cte('assay_plates')
            j = _p
            j = j.join(_c, _c.c.copy_id==_p.c.copy_id)
            j = j.join(_l, _c.c.library_id==_l.c.library_id)
            j = j.join(_pl, _p.c.plate_location_id==_pl.c.plate_location_id,
                isouter=True)
            j = j.join(_assay_plates_query, _p.c.plate_id
                ==_assay_plates_query.c.plate_id)
            _extant_query = ( 
                select([
                    _assay_plates_query.c.activity_id,
                    _concat(
                        _l.c.short_name, '/', _c.c.name, '/', 
                        cast(_p.c.plate_number, sqlalchemy.sql.sqltypes.Text)
                    ).label('plate_key'),
                    _l.c.short_name.label('library_short_name'),
                    _l.c.screening_status.label('library_screening_status'),
                    _l.c.screen_type.label('library_screen_type'),
                    _c.c.name.label('copy_name'),
                    _c.c.comments.label('copy_comments'),
                    _concat_with_sep(
                        args=[_pl.c.room,_pl.c.freezer,_pl.c.shelf,_pl.c.bin],
                        sep='-').label('plate_location'),
                    _c.c.usage_type,
                    _p.c.plate_number,
                    _p.c.cplt_screening_count,
                    # TODO: consider using the cumulative_freeze_thaw_count
                    _p.c.screening_count,
                    _p.c.remaining_well_volume,
                    _p.c.status,
                ])
                .select_from(j)
                .order_by(
                    _assay_plates_query.c.activity_id,
                    _l.c.short_name, _c.c.name, _p.c.plate_number ))
            _extant_query = _extant_query.cte('extant')
            
            # 2. query for searched plates
            logger.info('searched_plate_ids: %r', searched_plate_ids)
            j = _p
            j = j.join(_c, _c.c.copy_id==_p.c.copy_id)
            j = j.join(_l, _c.c.library_id==_l.c.library_id)
            j = j.join(_pl, _p.c.plate_location_id==_pl.c.plate_location_id,
                isouter=True)
            _search_query = ( 
                select([
                    literal_column('0').label('activity_id'),
                    _concat(
                        _l.c.short_name, '/', _c.c.name, '/', 
                        cast(_p.c.plate_number, sqlalchemy.sql.sqltypes.Text)
                    ).label('plate_key'),
                    _l.c.short_name.label('library_short_name'),
                    _l.c.screening_status.label('library_screening_status'),
                    _l.c.screen_type.label('library_screen_type'),
                    _c.c.name.label('copy_name'),
                    _c.c.comments.label('copy_comments'),
                    _concat_with_sep(
                        args=[_pl.c.room,_pl.c.freezer,_pl.c.shelf,_pl.c.bin],
                        sep='-').label('plate_location'),
                    _c.c.usage_type,
                    _p.c.plate_number,
                    _p.c.cplt_screening_count,
                    # TODO: consider using the cumulative_freeze_thaw_count
                    _p.c.screening_count,
                    _p.c.remaining_well_volume,
                    _p.c.status,
                ])
                .select_from(j)
                .where(_p.c.plate_id.in_(searched_plate_ids))
                .order_by(_l.c.short_name, _c.c.name, _p.c.plate_number ))
            
            if extant_plate_numbers:
                _search_query = _search_query.cte('search_query')
                _combined_query = union(
                    select([literal_column(x) for x in fields])
                        .select_from(_extant_query),
                    select([literal_column(x) for x in fields])
                        .select_from(_search_query)
                    )
                _combined_query = _combined_query.order_by(
                    nullslast(desc(column('activity_id'))),
                    'library_short_name', 'copy_name','plate_number')    
            else:
                _combined_query = _search_query

            compiled_stmt = str(_combined_query.compile(
                dialect=postgresql.dialect(),
                compile_kwargs={"literal_binds": True}))
            logger.info('compiled_stmt %s', compiled_stmt)
                    
            _result = conn.execute(_combined_query)
            
            logger.info('executed, build table...')
            # Convert the query into plate-ranges
            _data = []
            librarycopy = None
            start_plate = 0
            end_plate = 0
            current_lc = None
            current_ls = None

            new_row = OrderedDict()
            errorDict = ErrorDict()
            plate_locations = set()
            plate_keys = set()
            for _row in cursor_generator(_result,fields):
                logger.debug('row: %r', _row)
                librarycopy = '{library_short_name}/{copy_name}'.format(**_row)
                plate_number = _row['plate_number']
                if ( librarycopy != current_lc 
                        or current_ls != _row['activity_id']
                        or end_plate < _row['plate_number']-1):
                    if new_row:
                        new_row['plate_locations'] = [x for x in plate_locations]
                        new_row['plate_keys'] = [x for x in plate_keys]
                        new_row['end_plate'] = end_plate
                        new_row['errors'] = errorDict.showErrors()
                        new_row['warnings'] = errorDict.showWarnings()
                        _data.append(new_row)
                    errorDict = ErrorDict()
                    plate_locations = set()
                    plate_keys = set()
                    new_row = OrderedDict((
                        ('library_screening_id', _row['activity_id']),
                        ('library_short_name', _row['library_short_name']),
                        ('library_screening_status', _row['library_screening_status']),
                        ('copy_name', _row['copy_name']),
                        ('copy_comments', _row['copy_comments']),
                        ('start_plate', _row['plate_number']),
                        ))
                    current_lc = librarycopy
                    current_ls = _row['activity_id']
                end_plate = _row['plate_number']
                errorDict.check_row(_row)
                plate_locations.add(_row['plate_location'])
                plate_keys.add(_row['plate_key'])
            if librarycopy:
                new_row['plate_locations'] = [x for x in plate_locations]
                new_row['plate_keys'] = [x for x in plate_keys]
                new_row['end_plate'] = end_plate
                new_row['errors'] = errorDict.showErrors()
                new_row['warnings'] = errorDict.showWarnings()
                _data.append(new_row)
        
        library_comments = self.get_library_comments(
            { _dict['library_short_name'] for _dict in _data})
        
        for _dict in _data:
            _dict['library_comment_array'] = \
                library_comments.get(_dict['library_short_name'],None)

        self.get_plate_comments_for_plate_range_data(_data, _combined_query)
        
        return _data

    def get_plate_comments_for_plate_range_data(self, _data, join_query):
        
        logger.info('get plate comments...')
        for _dict in _data:
            _dict['plate_comment_array'] = set()
        _comment_apilogs = \
            ApiLogResource.get_resource_comment_subquery('librarycopyplate')
        _apilogs = self.bridge['reports_apilog']
        _join_query = join_query.cte('join_query')
        _comment_apilogs = \
            _comment_apilogs.where(_apilogs.c.key.in_(
                select([_join_query.c.plate_key])
                .select_from(_join_query)
            ))
        _comment_apilogs = _comment_apilogs.cte('logs')
        query = (
            select([
                _comment_apilogs.c.key,
                func.array_agg(
                    _concat(                            
                        cast(_comment_apilogs.c.name,
                            sqlalchemy.sql.sqltypes.Text),
                        LIST_DELIMITER_SUB_ARRAY,
                        cast(_comment_apilogs.c.date_time,
                            sqlalchemy.sql.sqltypes.Text),
                        LIST_DELIMITER_SUB_ARRAY,
                        '(',_comment_apilogs.c.key, ') ',
                        _comment_apilogs.c.comment)
                    )
            ])
            .select_from(_comment_apilogs)
            .group_by(_comment_apilogs.c.key)
        )
        
        with get_engine().connect() as conn:
            for x in conn.execute(query):
                key = x[0]
                comment_array = x[1]
                for _dict in _data:
                    if key in _dict['plate_keys']:
                        _dict['plate_comment_array'].update(comment_array)
        for _dict in _data:
            _dict['plate_comment_array'] = list(_dict['plate_comment_array'])            

        logger.info('plate comments generated')

    def get_plate_comments_for_plate_range_data_bak(self, _data):
        
        logger.info('get plate comments...')
        cumulative_plate_keys = set()
        for _dict in _data:
            _dict['plate_comment_array'] = set()
            cumulative_plate_keys.update(_dict['plate_keys'])
        logger.info('cumulative_plate_keys: %r', cumulative_plate_keys)
        _comment_apilogs = \
            ApiLogResource.get_resource_comment_subquery('librarycopyplate').cte('logs')
        query = (
            select([
                _comment_apilogs.c.key,
                func.array_agg(
                    _concat(                            
                        cast(_comment_apilogs.c.name,
                            sqlalchemy.sql.sqltypes.Text),
                        LIST_DELIMITER_SUB_ARRAY,
                        cast(_comment_apilogs.c.date_time,
                            sqlalchemy.sql.sqltypes.Text),
                        LIST_DELIMITER_SUB_ARRAY,
                        '(',_comment_apilogs.c.key, ') ',
                        _comment_apilogs.c.comment)
                    )
            ])
            .select_from(_comment_apilogs)
            .group_by(_comment_apilogs.c.key)
            .where(_comment_apilogs.c.key.in_(cumulative_plate_keys))
        )
        
        with get_engine().connect() as conn:
            for x in conn.execute(query):
                key = x[0]
                comment_array = x[1]
                for _dict in _data:
                    if key in _dict['plate_keys']:
                        _dict['plate_comment_array'].update(comment_array)
        for _dict in _data:
            _dict['plate_comment_array'] = list(_dict['plate_comment_array'])            

        logger.info('plate comments generated')

    def get_library_comments(self, library_keys):

        logger.info('get library comments...')

        _comment_apilogs = \
            ApiLogResource.get_resource_comment_subquery('library').cte('logs')
        query = (
            select([
                _comment_apilogs.c.key,
                func.array_agg(
                    _concat(                            
                        cast(_comment_apilogs.c.name,
                            sqlalchemy.sql.sqltypes.Text),
                        LIST_DELIMITER_SUB_ARRAY,
                        cast(_comment_apilogs.c.date_time,
                            sqlalchemy.sql.sqltypes.Text),
                        LIST_DELIMITER_SUB_ARRAY,
                        _comment_apilogs.c.comment)
                    )
            ])
            .select_from(_comment_apilogs)
            .group_by(_comment_apilogs.c.key)
            .where(_comment_apilogs.c.key.in_(library_keys))
        )
        
        with get_engine().connect() as conn:
            comments = defaultdict(list)
            for x in conn.execute(query):
                comments[x[0]] = x[1]
            logger.info('library comments generated')
            return comments

    def get_query(self, schema, param_hash, user):
        '''  
        LibraryScreeningResource
        - NOTE: overrides the ActivityResource get_query:
        - TODO: LibraryScreeningResource does not need to extend ActivityResource
        
        special use case: 
        - use "library_plates_screened__contains" param to query for 
        LibraryScreenings that have screened a particular set of plate-copies
        '''

        manual_field_includes = set(param_hash.get('includes', []))
        
        library_plates_screened_search = param_hash.pop(
            'library_plates_screened__contains', None)
        if library_plates_screened_search:
            manual_field_includes.add('library_plates_screened')
            
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filter_expression = \
            self._meta.authorization.filter(user,filter_expression)

              
        order_params = param_hash.get('order_by', [])
        order_params.append('-date_of_activity')
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
         
        # specific setup
        _a = self.bridge['activity']
        _la = self.bridge['lab_activity']
        _screening = self.bridge['screening']
        _library_screening = self.bridge['library_screening']
        _screen = self.bridge['screen']
        _ap = self.bridge['assay_plate']
        _c = self.bridge['copy']
        _p = self.bridge['plate']
        j = _a
        j = j.join(_la, _a.c.activity_id == _la.c.activity_id)
        j = j.join(_screening, _screening.c.activity_id == _la.c.activity_id)
        j = j.join(
            _library_screening,
            _library_screening.c.activity_id == _la.c.activity_id)
        j = j.join(_screen, _la.c.screen_id == _screen.c.screen_id)

        # Get custom columns from the parent (ActivityResource); the
        # LibraryScreening query rows contain the needed activity and screen ids
        custom_columns = \
            super(LibraryScreeningResource, self).get_custom_columns('ls')
        custom_columns.update({
            'type': cast(case([
                (_library_screening.c.is_for_external_library_plates,
                        'externallibraryscreening')],
                    else_='libraryscreening'), sqlalchemy.sql.sqltypes.Text),
            'activity_class': cast(case([
                (_library_screening.c.is_for_external_library_plates,
                        'externallibraryscreening')],
                    else_='libraryscreening'), sqlalchemy.sql.sqltypes.Text),
            'libraries_screened_count': literal_column(
                '(select count(distinct(l.*)) from library l '
                'join copy using(library_id) join plate using(copy_id) '
                'join assay_plate ap using(plate_id) '
                'where library_screening_id='
                'library_screening.activity_id)' ),
             'library_plates_screened_count': literal_column(
                 '(select count(distinct(ap.plate_id)) '
                 'from assay_plate ap where library_screening_id='
                 'library_screening.activity_id)' ),
            })

        base_query_tables = [
            'activity', 'lab_activity', 'screening', 'library_screening',
            'screen'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        stmt = select(columns.values()).select_from(j)

        extra_params = {}
        if library_plates_screened_search:
            extra_params['plates_screened'] = None
            copy_plate_query = (
                select([
                    _p.c.plate_id,_p.c.plate_number, _c.c.name, 
                    _ap.c.library_screening_id])
                .select_from(_p.join(_c,_p.c.copy_id==_c.c.copy_id)
                    .join(_ap,_p.c.plate_id==_ap.c.plate_id)))
            
            logger.info('try to grep library_plates_screened_search: %r', 
                library_plates_screened_search)
            if isinstance(library_plates_screened_search, six.types.StringTypes):
                library_plates_screened_search = \
                    re.split(r'[,;]+', library_plates_screened_search)
            if not isinstance(library_plates_screened_search, (list,tuple)):
                library_plates_screened_search = (library_plates_screened_search,)
            or_clause = []
            for lps_search in library_plates_screened_search:
                logger.info('try: %r', lps_search)
                copy_plate_pattern = re.compile(r'([^\/]+)\/(\d+)')
                match = copy_plate_pattern.match(lps_search)
                if match: 
                    logger.info('found match %r for %r', match, lps_search)
                    or_clause.append(and_(
                        _c.c.name==match.group(1),_p.c.plate_number==match.group(2)))
            if or_clause:
                logger.info('library_copy_plate pattern subquery')
                copy_plate_query = copy_plate_query.where(or_(*or_clause))
                copy_plate_query = copy_plate_query.alias('cp_query')
                stmt = stmt.where(exists(
                    select([None]).select_from(copy_plate_query)
                        .where(copy_plate_query.c.library_screening_id
                            ==_library_screening.c.activity_id) ))
        
        if logger.isEnabledFor(logging.DEBUG):
            compiled_stmt = str(stmt.compile(
                dialect=postgresql.dialect(),
                compile_kwargs={"literal_binds": True}))
            logger.info('compiled_stmt %s', compiled_stmt)
        # general setup
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        logger.info('order_clauses: %r', order_clauses)
        stmt = stmt.order_by('activity_id')

        filename = self._get_filename(readable_filter_hash, schema, **extra_params)
        
        return (field_hash, columns, stmt, count_stmt,filename)
        
    @read_authorization
    def build_list_response(self, request, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
         
        is_for_detail = kwargs.pop('is_for_detail', False)
 
        try:
             
            (field_hash, columns, stmt, count_stmt,filename) = \
                self.get_query(schema, param_hash, request.user)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
            
            # create a generator to wrap the cursor and expand the library_plates_screened
            def create_lcp_gen(generator):
                bridge = self.bridge
                _library = self.bridge['library']
                _lcp = self.bridge['plate']
                _cp = self.bridge['copy']
                _ap = self.bridge['assay_plate']
                lcp_query = (
                    select([ 
                        _library.c.short_name,
                        _cp.c.name,
                        _ap.c.plate_number,
                     ])
                    .select_from(
                        _ap.join(_lcp, _ap.c.plate_id == _lcp.c.plate_id)
                            .join(_cp, _cp.c.copy_id == _lcp.c.copy_id)
                            .join(
                                _library,
                                _library.c.library_id == _cp.c.library_id))
                    .where(_ap.c.library_screening_id == text(':activity_id'))
                    .group_by(
                        _library.c.short_name, _cp.c.name, _ap.c.plate_number)
                    .order_by(
                        _library.c.short_name, _cp.c.name, _ap.c.plate_number))
                logger.debug('lcp_query: %r', str(lcp_query.compile()))
                
                def library_copy_plates_screened_generator(cursor):
                    if generator:
                        cursor = generator(cursor)
                    class Row:
                        def __init__(self, row):
                            self.row = row
                            self.entries = []
                            activity_id = row['activity_id']
                            query = conn.execute(
                                lcp_query, activity_id=activity_id)
                            copy = None
                            start_plate = None
                            end_plate = None
                            for x in query:
                                if not copy:
                                    copy = x[1]
                                    library = x[0]
                                if not start_plate:
                                    start_plate = end_plate = x[2]
                                if (x[0] != library 
                                    or x[1] != copy 
                                    or x[2] > end_plate + 1):
                                    # start a new range, save old range
                                    self.entries.append('%s:%s:%s-%s'
                                        % (library, copy, start_plate, end_plate))
                                    start_plate = end_plate = x[2]
                                    copy = x[1]
                                    library = x[0]
                                else:
                                    end_plate = x[2]
                            if copy: 
                                self.entries.append('%s:%s:%s-%s'
                                    % (library, copy, start_plate, end_plate))
                                    
                        def has_key(self, key):
                            if key == 'library_plates_screened': 
                                return True
                            return self.row.has_key(key)
                        def keys(self):
                            return self.row.keys();
                        def __getitem__(self, key):
                            if key == 'library_plates_screened':
                                return self.entries
                            else:
                                return self.row[key]
                    conn = get_engine().connect()
                    try:
                        for row in cursor:
                            yield Row(row)
                    finally:
                        conn.close()
                        
                return library_copy_plates_screened_generator
            
            if 'library_plates_screened' in field_hash:
                rowproxy_generator = create_lcp_gen(rowproxy_generator)

            rowproxy_generator = \
               self._meta.authorization.get_row_property_generator(
                   request.user, field_hash, rowproxy_generator)
                    
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
             
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
              
        except Exception, e:
            logger.exception('on get list')
            raise e  

    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_detail must be implemented')
    
    def validate(self, _dict, patch=False, schema=None):

        errors = ActivityResource.validate(self, _dict, schema=schema, patch=patch)
        if _dict.get('library_plates_screened', None):
            if bool(_dict.get('is_for_external_library_plates', False)):
                errors['library_plates_screened'] = (
                    'Can not specifiy library plates if '
                    '"is_for_external_library_plates"')

        return errors


    def build_patch_detail(self, request, deserialized, log=None, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(deserialized,schema=schema, validate=True,**kwargs)

        # NOTE: 20170321 not creating a screen "parent" log;
        # libraryScreening activities will stand on their own
        # original_screen_data = self.get_screen_resource()._get_detail_response_internal(**{
        #     'facility_id': original_data['screen_facility_id']})


        original_data = self._get_detail_response_internal(**kwargs_for_log)
        logger.debug('original libraryscreening data: %r', original_data)
        
        # NOTE: creating a log, even if no data have changed (may be comment only)
        log = self.make_log(request)
        # FIXME: Set the log URI using the containing screen URI
        log.key = '/'.join([str(kwargs_for_log[x]) for x in id_attribute])
        log.uri = '/'.join([log.ref_resource_name,log.key])
        log.save()

        obj = None
        plate_meta = None
        if deserialized:
            patch_result = self.patch_obj(request, deserialized, log=log, **kwargs)
            obj = patch_result[API_RESULT_OBJ]
            plate_meta = patch_result[API_RESULT_META]
        for id_field in id_attribute:
            val = getattr(obj, id_field,None)
            if val:
                kwargs_for_log['%s' % id_field] = val
        log.key = '/'.join([str(kwargs_for_log[x]) for x in id_attribute])
        screen_facility_id = obj.screen.facility_id
        log.uri = '/'.join([
            'screen', screen_facility_id,log.ref_resource_name,log.key])
        log.save()
        logger.info('log saved: %r', log)
        
        # TODO: create a log for the parent screen:
        # For now, the log.uri can be used to query for logs for screen

        new_data = self._get_detail_response_internal(**kwargs_for_log)
        kwargs_for_log[HTTP_PARAM_USE_VOCAB] = True
        new_data_display = self._get_detail_response_internal(**kwargs_for_log)
        self.log_patch(request, original_data,new_data,log=log, **kwargs)
        log.save()

        meta = { 
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 
                new_data_display['library_plates_screened_count'],
            'Volume per well transferred from Plates (nL)': 
                new_data_display['volume_transferred_per_well_from_library_plates'],
        }
        meta.update(plate_meta)
        _data = { API_RESULT_DATA: [new_data] }
        _data[API_RESULT_META] = { API_MSG_RESULT: meta }
        
        return _data

    def build_post_detail(self, request, deserialized, log=None, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        kwargs_for_log = self.get_id(
            deserialized,schema=schema, validate=False,**kwargs)
        
        id_attribute = schema['id_attribute']
        
        logger.info('post detail: %r, %r', kwargs_for_log, id_attribute)

        # NOTE: 20170321 not creating a screen "parent" log;
        # libraryScreening activities will stand on their own
        # original_screen_data = self.get_screen_resource()._get_detail_response_internal(**{
        #     'facility_id': original_data['screen_facility_id']})

        original_data = None
        log = self.make_log(request)
        log.save()
        if kwargs_for_log and len(kwargs_for_log.items())==len(id_attribute):
            # A full id exists, query for the existing state
            try:
                original_data = self._get_detail_response_internal(**kwargs_for_log)
            except Exception, e: 
                logger.exception('exception when querying for existing obj: %s', 
                    kwargs_for_log)
            
            if original_data is not None and len(original_data) != 0:
                raise ValidationError({ 
                    k: '%r Already exists' % v for k,v in kwargs_for_log.items() })
            original_data = None
        
            # NOTE: create a log if possible, with id_attribute, for downstream
            log.key = '/'.join([str(kwargs_for_log[x]) for x in id_attribute])
            log.uri = '/'.join([
                'screen', deserialized.get('screen_facility_id', None),
                log.ref_resource_name,log.key])
            log.save()
            logger.info('log saved: %r', log)

        patch_result = self.patch_obj(request, deserialized, log=log, **kwargs)
        obj = patch_result[API_RESULT_OBJ]
        plate_meta = patch_result[API_RESULT_META]
        logger.info('patch meta: %r', plate_meta)
        for id_field in id_attribute:
            val = getattr(obj, id_field,None)
            if val:
                kwargs_for_log['%s' % id_field] = val
        # Note: update the log with the id_attribute after object is created
        log.key = '/'.join([str(kwargs_for_log[x]) for x in id_attribute])
        log.uri = '/'.join([
            'screen', deserialized.get('screen_facility_id', None),
            log.ref_resource_name,log.key])
        log.save()
        
        obj.apilog_uri = log.log_uri
        obj.save()

        # get new state, for logging
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        if not new_data:
            raise BadRequest('no data found for the new obj created by post: %r', obj)
        self.log_patch(
            request, original_data,new_data,log=log, 
            full_create_log=True, **kwargs)
        log.save()

        # FIXME: create a log for the parent screen

        kwargs_for_log[HTTP_PARAM_USE_VOCAB] = True
        new_data_display = self._get_detail_response_internal(**kwargs_for_log)
        
        meta = { 
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 
                new_data_display['library_plates_screened_count'],
            'Volume per well transferred from Plates (nL)': 
                new_data_display['volume_transferred_per_well_from_library_plates'],
        }
        meta.update(plate_meta)
        _data = { API_RESULT_DATA: [new_data] }
        _data[API_RESULT_META] = { API_MSG_RESULT: meta }
        
        return _data

    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs.pop('schema', None)
        if schema is None:
            raise Exception('schema not initialized')
        fields = schema['fields']
        initializer_dict = {}
        ls_log = kwargs.get('log', None)
        if ls_log is None:
            raise BadRequest(
                'library screening log should be created by the callee of patch_obj')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.debug('param_hash: %r', param_hash)

        # FIXME: parse and validate only editable fields

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        patch = bool(id_kwargs)
        initializer_dict = self.parse(deserialized, create=not patch, schema=schema)
        errors = self.validate(initializer_dict, schema=schema, patch=patch)
        if errors:
            raise ValidationError(errors)
        
        if not patch:
            _key = 'screen_facility_id'
            _val = deserialized[_key]
            try:
                screen = Screen.objects.get(facility_id=_val)
                initializer_dict['screen'] = screen
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='does not exist: {val}'.format(val=_val))

        _key = 'performed_by_user_id'
        _val = deserialized.get(_key, None)
        if _val:
            try:
                performed_by_user = ScreensaverUser.objects.get(screensaver_user_id=_val)
                initializer_dict['performed_by'] = performed_by_user
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='does not exist: {val}'.format(val=_val))
        library_screening = None
        if patch:
            try:
                logger.info('%r', id_kwargs)
                library_screening = LibraryScreening.objects.get(
                    pk=id_kwargs['activity_id'])
                current_volume_transferred_per_well = \
                    library_screening.volume_transferred_per_well_from_library_plates
                new_volume_xfer = initializer_dict.get(
                    'volume_transferred_per_well_from_library_plates',None) 
                if ( new_volume_xfer is not None and 
                    current_volume_transferred_per_well
                        != new_volume_xfer ):
                    if library_screening.assayplate_set.exists():
                        raise ValidationError({
                            'volume_transferred_per_well_from_library_plates':
                                ( 
                                'Can not be changed if plates have been assigned; '
                                '(%r to %r' 
                                % (current_volume_transferred_per_well,new_volume_xfer)),
                            'library_plates_screened': (
                                'Remove all plates before changing volume assigned')
                            })
            except ObjectDoesNotExist:
                raise Http404(
                    'library_screening does not exist for: %r', id_kwargs)
        else:
            # Set the created_by field:
            # NOTE: after migration LibraryScreening will no longer have 
            # this field
            try:
                adminuser = ScreensaverUser.objects.get(username=request.user.username)
            except ObjectDoesNotExist as e:
                logger.error('admin user: %r does not exist', request.user.username )
                raise
            
            library_screening = LibraryScreening()
            library_screening.created_by = adminuser
        
        model_field_names = [
            x.name for x in library_screening._meta.get_fields()]
        for key, val in initializer_dict.items():
            if key in model_field_names:
                setattr(library_screening, key, val)

        logger.info('save library screening, fields: %r', library_screening)
        library_screening.save()
        library_plates_screened = deserialized.get(
            'library_plates_screened', None)
        if not library_plates_screened and not patch:
            raise ValidationError(
                key='library_plates_screened', 
                msg='required')
        # override_param = parse_val(
        #     param_hash.get(API_PARAM_OVERRIDE, False),
        #         API_PARAM_OVERRIDE, 'boolean')
        # override_vol_param = parse_val(
        #     param_hash.get(API_PARAM_VOLUME_OVERRIDE, False),
        #         API_PARAM_VOLUME_OVERRIDE, 'boolean')
        
        plate_meta = {}
        if library_plates_screened is not None:
            logger.debug('save library screening, set assay plates: %r', library_screening)
            plate_meta = self._set_assay_plates(
                request, schema, 
                library_screening, library_plates_screened, ls_log) 
            logger.info('save library screening, assay plates set: %r', library_screening)
        
            ls_log.json_field = json.dumps(plate_meta)
            logger.info('parent_log: %r', ls_log)
        self.create_screen_screening_statistics(library_screening.screen)
        self.create_screened_experimental_well_count(library_screening)
        
        return { API_RESULT_OBJ: library_screening, API_RESULT_META: plate_meta}
        
    def create_screened_experimental_well_count(self, library_screening):
        # TODO: this should be dynamic
        with get_engine().connect() as conn:
            sql = (
                'select count(*) '
                'from well '
                'where plate_number in ( '
                '    select distinct(plate_number) from assay_plate ' 
                '    where library_screening_id = %s ) '
                "and library_well_type='experimental'")
            library_screening.screened_experimental_well_count = int(
                conn.execute(
                    sql, (library_screening.activity_id))
                .scalar() or 0)
            library_screening.save()

    def create_screen_screening_statistics(self, screen):
        
        # TODO: these statistics must be updated when the library definitions are
        # reloaded, see LibrariesDAO in SS1 and LibraryContentsLoader
        
        with get_engine().connect() as conn:
            # NOTE: mixing Django connection with SQA connection
            # - thrown exceptions will rollback the nested SQA transaction
            # see: http://docs.sqlalchemy.org/en/latest/core/connections.html
            screen_facility_id = screen.facility_id
            sql = (
                'select count(w.well_id) '
                'from Well w, Screen s  '
                'join assay_plate ap using(screen_id) '
                'join plate p using(plate_id) '
                'join library_screening ls on(ap.library_screening_id=ls.activity_id) '
                'where   ap.replicate_ordinal = 0  '
                'and w.plate_number = p.plate_number ' 
                'and w.library_well_type = %s'
                'and s.facility_id = %s;')
            screen.screened_experimental_well_count = int(
                conn.execute(
                    sql, ('experimental', screen_facility_id))
                .scalar() or 0)
            sql = (
                'select count(distinct(w.well_id)) '
                'from Well w, Screen s  '
                'join assay_plate ap using(screen_id) '
                'join plate p using(plate_id) '
                'join library_screening ls on(ap.library_screening_id=ls.activity_id) '
                'where   ap.replicate_ordinal = 0  '
                'and w.plate_number = p.plate_number ' 
                'and w.library_well_type = %s'
                'and s.facility_id = %s;')
            screen.unique_screened_experimental_well_count = int(
                conn.execute(
                    sql, ('experimental', screen_facility_id))
                .scalar() or 0)
            logger.info('screen_screening_statistics: %d, %d',
                screen.screened_experimental_well_count, 
                screen.unique_screened_experimental_well_count)
            screen.save()
    
    @transaction.atomic    
    def _set_assay_plates(
            self, request, schema, library_screening, library_plates_screened,
            ls_log):
        '''
        - Create new assay plates
        - Adjust librarycopyplate volume
        - Adjust librarycopyplate screening count
        - Create librarycopyplate logs
        - TODO: create copy logs
        '''
        logger.info('set assay plates screened for: %r, %r', 
            library_screening.activity_id, library_plates_screened)
        
        # Parse library_plate_ranges
        # E.G. Regex: /(([^:]+):)?(\w+):(\d+)-(\d+)/
        regex_string = schema['fields']['library_plates_screened']['regex']
        matcher = re.compile(regex_string)
        def show_plates(plates):
            return ['%s/%d' % (plate.copy.name, plate.plate_number) 
                for plate in plates]
        
        min_plate_volume_after_transfer = 0
        if library_screening.screen.screen_type == VOCAB_SCREEN_TYPE_RNAI:
            min_plate_volume_after_transfer = \
                self.MIN_WELL_VOL_RNAI
        elif library_screening.screen.screen_type == VOCAB_SCREEN_TYPE_SM:
            min_plate_volume_after_transfer = \
                self.MIN_WELL_VOL_SMALL_MOLECULE
            
        # 1. Validate and parse the library_plates_screened input patterns
        new_library_plates_screened = []
        for lps in library_plates_screened:
            match = matcher.match(lps)
            if not match:
                raise ValidationError(
                    key='library_plates_screened',
                    msg=('%r does not match pattern: %s' 
                        % (lps, regex_string)))
                break
            else:
                start_plate = int(match.group(4))
                if match.group(6):
                    end_plate = int(match.group(6))
                else:
                    end_plate = start_plate
                new_library_plates_screened.append({
                    'library_short_name': match.group(2),
                    'copy_name': match.group(3),
                    'start_plate': start_plate,
                    'end_plate': end_plate,
                     })
        library_plates_screened = new_library_plates_screened
        
        # 2. Validate and find all the plates referenced in the ranges
        logger.info(
            'get the referenced plates for: %r', library_plates_screened)
        plate_ranges = []
        plate_keys = set()
        all_plate_ids = set() # extant ap's and new
        for _data in library_plates_screened:
            logger.info('lps data: %r', _data)
            try:
                copy_name = _data['copy_name']
                library_short_name = _data['library_short_name']
                copy = Copy.objects.get(
                    name=copy_name,
                    library__short_name=library_short_name)
                logger.info('found copy: %r', copy)
            except ObjectDoesNotExist:
                raise ValidationError(
                    key='library_plates_screened',
                    msg='{library_short_name}/{copy} does not exist: {val}'.format(
                        library_short_name=library_short_name,
                        copy=copy_name,
                        val=str(_data)))
            try:
                start_plate = Plate.objects.get(
                    copy=copy,
                    plate_number=_data['start_plate'])
                end_plate = Plate.objects.get(
                    copy=copy,
                    plate_number=_data['end_plate'])
                logger.info(
                    'found start: %r, end: %r plates', start_plate, end_plate)
                if start_plate.copy.library != end_plate.copy.library:
                    raise ValidationError(
                        key='library_plates_screened',
                        msg=('plate range must be for a single library: '
                             '{start_plate}-{end_plate}').format(**_data))
                if (start_plate.copy.library.screen_type 
                        != library_screening.screen.screen_type):
                    raise ValidationError(
                        key='library_plates_screened',
                        msg=('library.screen_type!=screen.screen_type: '
                             '{start_plate}-{end_plate} (%s)' 
                                % start_plate.copy.library.screen_type
                             ).format(**_data))
                plate_range = range(
                    start_plate.plate_number, end_plate.plate_number + 1)

                # REMOVED: 20170412 - per JAS/KR; allow multiple copies, same plate
                # to be screened in one screening
                # if plate_numbers & set(plate_range):
                #     raise ValidationError(
                #         key='library_plates_screened',
                #         msg=('A plate number can only be screened once per '
                #             'Library Screening: {start_plate}-{end_plate}'
                #             ).format(**_data))
#                 plate_numbers.update(plate_range)
                plate_keys.update([ '%s/%d' % (copy.name, plate_number) 
                    for plate_number in plate_range])
                logger.info('find the plate range: %s: %s-%s', 
                    copy.name, start_plate.plate_number, end_plate.plate_number)
                plate_range = Plate.objects.all().filter(
                    copy=copy,
                    plate_number__range=(start_plate.plate_number, end_plate.plate_number))
                plate_ranges.append(plate_range)
                all_plate_ids.update([x.plate_id for x in plate_range.all()])

            except ObjectDoesNotExist:
                logger.exception('plate range error')
                raise ValidationError({
                    'library_plates_screened':
                    ('plate range not found: {start_plate}-{end_plate}'
                        ).format(**_data),
                    'copy_name': _data['copy_name'] })
        logger.debug('plate keys: %r',plate_keys)
        
        logger.info('3. Cache current state for logging...')
        # Create a search criteria to poll the current plate state
        # TODO: cache and log the copy state as well
        _original_plate_data = []
        if library_screening.assayplate_set.exists():
            for ap in library_screening.assayplate_set.all():
                all_plate_ids.add(ap.plate.plate_id)

        logger.info('Cache plate data...')
        if all_plate_ids:
            _original_plate_data = \
                self.get_plate_resource()._get_list_response_internal(
                plate_ids=all_plate_ids,
                includes =['-library_comment_array', '-comment_array'])
        
        if not _original_plate_data:
            raise Exception('no original plate data found')        
        logger.info('Find extant assay plates that are kept, '
            'find and remove deleted assay plates...')
        extant_plates = set()
        deleted_plates = set()
        if library_screening.assayplate_set.exists():
            for ap in library_screening.assayplate_set.all():
                plate_key = '%s/%d' % (ap.plate.copy.name, ap.plate.plate_number)
                if plate_key in plate_keys:
                    extant_plates.add(ap.plate)
                else:
                    # 20161020: no longer tracking data_load actions for an 
                    # assay plate so this is removed:
                    # if ap.screen_result_data_loading:
                    #     raise ValidationError(
                    #         key='library_plates_screened',
                    #         msg=(
                    #             'Assay plate has data and cannot be removed: %d'
                    #             ) % ap.plate_number)
                    deleted_plates.add(ap.plate)
                    ap.delete()       
        logger.info('deleted plates: %r', show_plates(deleted_plates)) 

        logger.info('5. Create new assay plates...')
        created_plates = set()
        not_allowed_libraries = set()
        plate_errors = []
        plate_warnings = []
        for plate_range in plate_ranges:
            for plate in plate_range:
                if plate in created_plates:
                    plate_errors.append(
                        'plate: "%s/%s" overlaps an existing plate range'
                            % (plate.copy.name, plate.plate_number))    
                    continue
            
            for replicate in range(library_screening.number_of_replicates):
                for plate in plate_range:
                    if plate not in extant_plates:
                        # FIXME: 20170406 - allow other status types to be screened,
                        # show a warning if not "active"
                        if (plate.copy.library.screening_status 
                            in self.WARN_LIBRARY_SCREENING_STATUS):
                            not_allowed_libraries.add(plate.copy.library)
                        elif (plate.copy.library.screening_status 
                            not in self.ALLOWED_LIBRARY_SCREENING_STATUS):
                            plate_errors.append(
                                'plate: "%s/%s": Library status: "%s"' 
                                % (plate.copy.name, plate.plate_number,
                                    plate.library.screening_status))
                        if plate.status in self.WARN_PLATE_STATUS:
                            plate_warnings.append(
                                'plate: "%s/%s" status: "%s"'
                                    % (plate.copy.name,plate.plate_number,plate.status))
                        elif plate.status not in self.ALLOWED_PLATE_STATUS:
                            plate_errors.append(
                                'plate: "%s/%s" status: "%s"'
                                    % (plate.copy.name,plate.plate_number,plate.status))
                        if plate.copy.usage_type not in self.ALLOWED_COPY_USAGE_TYPE:
                            plate_errors.append(
                                'plate: "%s/%s": "%s"' 
                                % (plate.copy.name, plate.plate_number,
                                   plate.copy.usage_type))
                        # 20170407 allow libraryscreening even after 
                        #     # copywells have been created
                        # if plate.cplt_screening_count > 0:
                        #     # FIXME: 20170407 allow libraryscreening even after 
                        #     # copywells have been created
                        #     # implement
                        #     # CopyWellResource.allocateScreeningVolumes
                        #     
                        #     # FIXME:TODO: 20170412 - allow screening of 
                        #     # library_screning_plates after cherry pick volumes
                        #     # have been taken
                        #     
                        #     # Volume tracking will not work on the plate level
                        #     # after cherry pick volumes have been taken from the
                        #     # copy wells
                        #     plate_errors.append(
                        #         'plate: "%s/%s"; may not be screened after '
                        #         'cherry pick screenings (%d) have been performed'
                        #             % (plate.copy.name,plate.plate_number, 
                        #                plate.cplt_screening_count))
                        if plate_errors:
                            continue
                        # 2017041 - allow screening for "not_allowed" libraries
                        # if not_allowed_libraries and override_param is not True:
                        #    continue

                        ap = AssayPlate.objects.create(**{  
                            'plate': plate,
                            'plate_number': plate.plate_number,
                            'screen': library_screening.screen,
                            'library_screening': library_screening,
                            'replicate_ordinal': replicate
                            })
                        ap.save()
                        created_plates.add(plate)
                    else:
                        logger.debug('extant plate: %r', plate)
        if plate_errors:
            plate_errors = sorted(plate_errors)
            raise ValidationError(
                    key='library_plates_screened',
                    msg=sorted(plate_errors))

        if not_allowed_libraries:
            not_allowed_libraries = sorted([
                '%s - status: %s' % (l.short_name,l.screening_status)
                for l in not_allowed_libraries])
            # if override_param is not True:
            #     raise ValidationError({
            #         API_PARAM_OVERRIDE: 'required',
            #         'library_plates_screened': (
            #             'Override required to screen libraries that are '
            #             'not allowed'),
            #         'Libraries': not_allowed_libraries
            #         }
            #     )
        
        logger.info('Update the plate screening related plates (and copywells)...')
        # TODO: deprecate volume change after plates are created.
        # if current_volume_tranferred_per_well != new_volume_transferred_per_well:
        #     
        #     for plate in extant_plates:
        #         plate_key = '%s/%d' % (plate.copy.name, plate.plate_number)
        #         logger.debug('plate: %r, remaining vol: %r, %r', 
        #             plate_key, plate.remaining_well_volume, 
        #             new_volume_transferred_per_well)
        #         current_remaining_well_volume = \
        #             plate.remaining_well_volume or Decimal(0)
        #         current_remaining_well_volume += current_volume_transferred_per_well
        #         new_remaining_well_volume = \
        #             current_remaining_well_volume - new_volume_transferred_per_well
        #         
        #         if new_remaining_well_volume < 0:
        #             # 20170407 - per JAS/KR,
        #             # raise an Error instead for insufficient vol
        #             logger.info('plate: %r, insufficient vol: %r', 
        #                 plate_key, new_remaining_well_volume)
        #             plates_insufficient_volume.append(
        #                 (plate_key, 
        #                     '(available: %s uL)' 
        #                         % lims_utils.convert_decimal(
        #                             current_remaining_well_volume, 1e-6, 1),
        #                     '(requested: %s uL)' 
        #                         % lims_utils.convert_decimal(
        #                             new_volume_transferred_per_well, 1e-6, 1)))
        # 
        #         plate.remaining_well_volume = new_remaining_well_volume
        #         plate.save()
        #         
        #         # FIXME/TODO: implement copywell vol adj
        #         # NOTE: usually, screening copies should not have copywells
        #         # TODO: implement this if screening policy is changed to allow
        #         cw_check_query = (
        #             plate.copywell_set.exclude(volume=F('initial_volume'))
        #                 .exclude(volume__isnull=True))
        #         if cw_check_query.exists():
        #             logger.info('copywells found: %r', [x for x in cw_check_query.all()])
        # 
        #             self.get_copywell_resource().allocate_library_screening_volumes(
        #                 )
        #             
        #             # raise ValidationError(
        #             #     key='library_plates_screened',
        #             #     msg=('Can not create a library screening if copy wells have '
        #             #     'been adjusted, plate: %r' % plate_key))
        
        logger.info('Update and check created plates: %r' ,
            show_plates(created_plates))
        volume_to_transfer = library_screening.volume_transferred_per_well_from_library_plates
        if not volume_to_transfer:
            raise ValidationError(
                key='volume_transferred_per_well_from_library_plates',
                msg='required')
        plates_insufficient_volume = []
        plate_copywell_stats = {}
        plate_copywell_warnings = {}
        copywell_allocation_meta = {
            API_MSG_LCPS_INSUFFICIENT_VOLUME: plate_copywell_warnings,
            API_MSG_COPYWELLS_ALLOCATED: plate_copywell_stats
        }
        for plate in created_plates:
            plate_key = '%s/%d' % (plate.copy.name, plate.plate_number)
            logger.info('plate: %r, remaining vol: %r, %r', 
                plate_key, plate.remaining_well_volume, 
                volume_to_transfer)
            plate.screening_count = (plate.screening_count or 0) + 1
            remaining_well_volume = plate.remaining_well_volume or Decimal(0)
            remaining_well_volume -= volume_to_transfer
            if remaining_well_volume < min_plate_volume_after_transfer:
                # 20170605 - JAS - allowed, but show warning
                logger.info('plate: %r, insufficient vol: %r', 
                    plate_key, remaining_well_volume)
                plates_insufficient_volume.append(
                    (plate_key, 
                        '(available: %s uL)' 
                            % lims_utils.convert_decimal(
                                plate.remaining_well_volume, 1e-6, 1),
                        '(requested: %s nL)' 
                            % lims_utils.convert_decimal(
                                volume_to_transfer, 1e-9, 1)))
            cw_check_query = (
                plate.copywell_set.exclude(volume=F('initial_volume'))
                    .exclude(volume__isnull=True))
            if cw_check_query.exists():
                logger.info('libraryscreening copywells to allocate found: %r', 
                    [x for x in cw_check_query.all()])
                # NOTE: parent_log is set to the library_screening_log:
                # (the log tree is flattened to one level)
                meta = self.get_copywell_resource()._allocate_well_volumes(
                    volume_to_transfer, cw_check_query.all(), ls_log)
                if API_MSG_LCPS_INSUFFICIENT_VOLUME in meta:
                    plate_copywell_warnings[plate_key] = meta[API_MSG_LCPS_INSUFFICIENT_VOLUME]
                plate_copywell_stats[plate_key] = meta[API_MSG_COPYWELLS_ALLOCATED] 
            # if cw_check_query.exists():
            #     logger.info('copywells found: %r', [x for x in cw_check_query.all()])
            #     raise ValidationError(
            #         key='library_plates_screened',
            #         msg=('Can not create a library screening if copy wells have '
            #         'been adjusted, plate: %r' % plate_key))
            plate.remaining_well_volume = remaining_well_volume
            plate.save()
                        
        plate_copywell_deallocate_stats = {}
        for plate in deleted_plates:
            plate_key = '%s/%d' % (plate.copy.name, plate.plate_number)
            plate.screening_count -= 1
            remaining_well_volume = plate.remaining_well_volume or Decimal(0)
            remaining_well_volume += volume_to_transfer
            plate.remaining_well_volume = remaining_well_volume
            cw_check_query = (
                plate.copywell_set.exclude(volume=F('initial_volume'))
                    .exclude(volume__isnull=True))
            if cw_check_query.exists():
                logger.info('libraryscreening copywells to deallocate found: %r', 
                    [x for x in cw_check_query.all()])
                # NOTE: parent_log is set to the library_screening_log:
                # (the log tree is flattened to one level)
                meta = self.get_copywell_resource()\
                    ._deallocate_well_volumes(
                        volume_to_transfer, cw_check_query.all(), ls_log)
                plate_copywell_deallocate_stats[plate_key] = \
                    meta[API_MSG_COPYWELLS_DEALLOCATED]

            plate.save()
        # Volume check
        if plates_insufficient_volume:
            # Modified: 20170407 - per JAS/KR,
            # raise an Error instead for insufficient vol
            plates_insufficient_volume = sorted(plates_insufficient_volume)
            # msg = '%d plates' % len(plates_insufficient_volume)
            # if library_screening.screen.screen_type == VOCAB_SCREEN_TYPE_SM:
            #     extra_msg = (' (%s uL is required for Small Molecule)'
            #         %  lims_utils.convert_decimal(
            #             self.MIN_WELL_VOL_SMALL_MOLECULE, 1e-6, 1))
            #     msg += extra_msg
            # raise ValidationError({
            #     # API_PARAM_VOLUME_OVERRIDE: 'required',
            #     API_MSG_LCPS_INSUFFICIENT_VOLUME: msg,
            #     'library_plates_screened': plates_insufficient_volume
            #     })

        if all_plate_ids:
            logger.info('Fetch the new Plate state: '
                'log plate volume changes, screening count...')
            _new_plate_data = self.get_plate_resource()._get_list_response_internal(
                plate_ids=all_plate_ids,
                includes =['-library_comment_array', '-comment_array'])
            plate_logs = self.get_plate_resource().log_patches(
                request, _original_plate_data, _new_plate_data,
                parent_log=ls_log, api_action=API_ACTION_PATCH)
            logger.info('plate_logs created: %d', len(plate_logs))
        
        # TODO: log the copy state as well

        logger.info('plates: updated: %r, created: %r, deleted: %r', 
            show_plates(extant_plates), show_plates(created_plates), 
            show_plates(deleted_plates))
        
        meta = {
            API_MSG_SCREENING_ADDED_PLATE_COUNT: len(created_plates),
            API_MSG_SCREENING_EXTANT_PLATE_COUNT : len(extant_plates),
            API_MSG_SCREENING_DELETED_PLATE_COUNT: len(deleted_plates),
            API_MSG_SCREENING_PLATES_UPDATED: 
                (len(created_plates) + len(deleted_plates))
            }
        if created_plates:
            meta.update(copywell_allocation_meta)
        if plate_copywell_deallocate_stats:
            meta[API_MSG_COPYWELLS_DEALLOCATED] = plate_copywell_deallocate_stats
        warnings = {}
        if not_allowed_libraries:
            warnings['library_screening_status'] = \
                "library status: " \
                + ', '.join(not_allowed_libraries)
        if plate_warnings:
            warnings['plate_status'] = sorted(plate_warnings)
        if warnings:
            meta[API_MSG_WARNING] = warnings
        
        logger.info('return meta information: %r', meta)
        # MODIFIED: 20170605 - per JAS/KR, overdraw volume but warn
        if plates_insufficient_volume:
           meta[API_MSG_LCPS_INSUFFICIENT_VOLUME] = plates_insufficient_volume
        return meta
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, log=None, **kwargs):
        
        raise NotImplementedError(
            'Library Screening may not be deleted - remove plates to negate')
        #         
        #         if log is None:
        #             raise BadRequest('log must be set: %r', kwargs)
        #         
        #         activity_id = kwargs.get('activity_id', None)
        #         if activity_id:
        #             try:
        #                 ls = LibraryScreening.objects.get(activity_id=activity_id)
        #                 plates = set()
        #                 for ap in ls.assayplate_set.all():
        #                     plates.add(ap.plate)
        #                 
        #                 for plate in plates:
        #                     plate.screening_count -= 1
        #                     remaining_well_volume = plate.remaining_well_volume or 0
        #                     remaining_well_volume += ls.volume_transferred_per_well_from_library_plates
        #                     plate.remaining_well_volume = remaining_well_volume
        #                     plate.save()
        # 
        #                 # return meta information for logging
        #                 return {
        #                     API_MSG_SCREENING_ADDED_PLATE_COUNT: 0,
        #                     API_MSG_SCREENING_EXTANT_PLATE_COUNT : 0,
        #                     API_MSG_SCREENING_DELETED_PLATE_COUNT: len(plates)
        #                     }
        #             except ObjectDoesNotExist:
        #                 logger.warn('no such library_screening: %s' % activity_id)
        #                 raise Exception(
        #                     'library_screeningfor activity_id: %s not found' 
        #                     % activity_id)
        #         else:
        #             raise Exception(
        #                 'library_screening delete action requires an activity_id %s' 
        #                 % kwargs)

class RawDataTransformerResource(DbApiResource):

    ERROR_CONTROL_PARSE = 'Parse errors'
    ERROR_CONTROL_DUPLICATES = lims_utils.ERROR_DUPLICATE_WELLS
    ERROR_CONTROL_WELL_TYPE = (
        'Control wells must be \'empty\', \'DMSO\', or \'Library Control\'')
    ERROR_MATRIX_SIZE_DETECTED = \
        'Matrix size detected: %d, does not match Assay Plate Size: %d'
    ERROR_COLLATION_COUNT = 'Number of collations: %d, '\
        'must be a divisor of the number matrices read: %d'
    ERROR_PLATE_COUNT = (
        'Plates required (%d) does not match # of plates entered (%d): '
        'Matrices read: %d (transformed: %d), '
        'Collations: %d')
    
    control_type_abbreviations = {
        'assay_positive_controls': 'P',
        'assay_negative_controls': 'N',
        'assay_other_controls': 'O',
        'library_controls': 'C'
    }
    library_well_type_abbreviations = {
        'experimental': 'X',
        'empty': 'E',
        'dmso': 'D',
        'library_control': 'C',
        'rnai_buffer': 'B' 
    }

    class Meta:

        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'rawdatatransform'
        authorization = UserGroupAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 
        
    def __init__(self, **kwargs):
        super(RawDataTransformerResource, self).__init__(**kwargs)
        self.lcp_resource = None
        self.reagent_resource = None
        self.labcherrypick_resource = None

    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<screen_facility_id>([\w]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<screen_facility_id>([\w]+))/"
                 r"(?P<cherry_pick_request_id>(\d+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
    
    def get_lcp_resource(self):
        if self.lcp_resource is None:
            self.lcp_resource = LibraryCopyPlateResource()
        return self.lcp_resource
    
    def get_reagent_resource(self):
        if self.reagent_resource is None:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource
    
    def get_labcherrypick_resource(self):
        if self.labcherrypick_resource is None:
            self.labcherrypick_resource = LabCherryPickResource()
        return self.labcherrypick_resource
    
    @read_authorization
    def get_detail(self, request, **kwargs):

        cherry_pick_request_id = kwargs.get('cherry_pick_request_id', None)
        screen_facility_id = kwargs.get('screen_facility_id', None)
        
        if screen_facility_id is None and cherry_pick_request_id is None:
            raise Http404('No screen or Cherry Pick Request Specified')
        
        screen = None
        if screen_facility_id:
            screen = Screen.objects.get(facility_id=screen_facility_id)
        cpr = None
        if cherry_pick_request_id:
            cpr = CherryPickRequest.objects.get(cherry_pick_request_id=cherry_pick_request_id)
            
        schema = kwargs.pop('schema')
        logger.info('get detail: %r', kwargs)
        schema = self.build_schema()

        _data = {}
        query = RawDataTransform.objects.all()
        if cpr:
            query = query.filter(cherry_pick_request=cpr)
            _data['cherry_pick_request_id'] = cpr.cherry_pick_request_id
        else:
            query = query.filter(screen=screen)
            _data['screen_facility_id'] = screen.facility_id
#         else:
#             query = query.filter(screen__isnull=True)
#         if cpr:
#             query = query.filter(cherry_pick_request=cpr)
#             _data['cherry_pick_request_id'] = cpr.cherry_pick_request_id
#         else:
#             query = query.filter(cherry_pick_request__isnull=True)    
            
        if not query.exists():
            raise Http404
        
        if query.count() != 1:
            raise Http404('Wrong number of objects returned: %d', query.count())
        rdt = query.all()[0]
        
        for key in schema['fields'].keys():
            if hasattr(rdt, key):
                _data[key] = getattr(rdt,key)
        
        if rdt.rawdatainputfile_set.all().exists():
            input_files = []
            _data['input_files'] = input_files
            for rdif in rdt.rawdatainputfile_set.all():
                logger.info('found input file: %r', rdif)
                input_file = {}
                input_files.append(input_file)
                for key in schema['input_file_fields'].keys():
                    if hasattr(rdif, key):
                        input_file[key] = getattr(rdif,key)
        
        final_data = { API_RESULT_DATA: _data }
        
        return self.build_response(request, _data)
        
    def build_schema(self, user=None):
        schema = DbApiResource.build_schema(self, user=user)    
    
        extra_fields = self.get_field_resource()._build_fields(
            scopes=['fields.rawdatainputfile','fields.rawdataoutput'])
        
        schema['input_file_fields'] = { field['key']:field 
            for field in extra_fields if field['scope'] == 'fields.rawdatainputfile' }
        schema['output_fields'] = { field['key']:field 
            for field in extra_fields if field['scope'] == 'fields.rawdataoutput' }
    
        return schema
    
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        logger.info('params: %r', param_hash.keys())
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)

            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)

                  
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
               self._meta.authorization.get_row_property_generator(
                   request.user, field_hash, rowproxy_generator)
 
            # specific setup
            _rawdatatransform = self.bridge['raw_data_transform']
            _rawdatainput = self.bridge['raw_data_input_file']
            _screen = self.bridge['screen']
            _cpr = self.bridge['cherry_pick_request']
            
            j = _rawdatatransform
            j = j.join(
                _screen, _rawdatatransform.c.screen_id==_screen.c.screen_id,
                isouter=True)
            j = j.join(
                _cpr, _rawdatatransform.c.cherry_pick_request_id==_cpr.c.cherry_pick_request_id,
                isouter=True)
            
            custom_columns = {
#                 'lookup_pmid': literal_column("'lookup_pmid'"),
                }

            base_query_tables = ['raw_data_transform', 'screen','cherry_pick_request' ] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            if not order_clauses and filter_expression is None:
                _alias = Alias(stmt)
                stmt = select([text('*')]).select_from(_alias)
            stmt = stmt.order_by('-screen_facility_id')
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
             
        except Exception, e:
            logger.exception('on get_list %s' % self._meta.resource_name)
            raise e  

    def put_list(self, request, **kwargs):
        raise NotImplementedError(
            "Put list is not implemented for Raw Data Transform")
    
    def put_detail(self, request, **kwargs):
        raise NotImplementedError(
            "Post detail is not implemented for Raw Data Transform")
    
    def patch_list(self, request, **kwargs):
        raise NotImplementedError(
            "Patch list is not implemented for Raw Data Transform")
    
    def patch_detail(self, request, **kwargs):
        raise NotImplementedError(
            "Patch detail is not implemented for Raw Data Transform")
    
    @write_authorization
    @transaction.atomic
    def post_detail(self, request, **kwargs):
        '''
        Parse raw data input files containing plate read data:
        - Input matrices are collated based on the specified ordering of
        plates, conditions, readouts, and replicates.
        - Input matrix values are associated with library or cherry pick well
        data.
        - Assay and library control well data are parsed and associated with 
        raw data values
        - All input settings are stored on RawDataTransform and RawDataInputFile
        objects in the database.
        - Parsed data is returned in Excel spreadsheet format but not stored in 
        the database.
        '''
        
        logger.info('post_detail...')
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
                
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        if len(request.FILES) == 0:
            raise ValidationError(key='input_file', msg='Required')
        logger.info('request.FILES: %r', request.FILES.keys())
        
        fields = schema['fields']
        initializer_dict = self.parse(param_hash, fields=fields)
        logger.info('initializer: %r', initializer_dict)
        
        errors = self.validate(initializer_dict, patch=False, schema=schema)
        if errors:
            raise ValidationError(errors)
        # Expand plate ranges
        plate_ranges = initializer_dict['plate_ranges']
        lcp_resource = self.get_lcp_resource()
        (plates,parsed_searches, errors) = lcp_resource.find_plates(plate_ranges)
        if errors:
            raise ValidationError(
                key='plate_ranges', 
                msg = 'Plates not found: %s' % ', '.join(sorted(errors)))
        logger.info('parsed_searches: %r', parsed_searches)
        plate_numbers = parsed_searches[0]['plate_numbers_listed_in_order']
        # check for duplicates
        plate_check_set = set()
        duplicate_plates = set()
        for plate_number in plate_numbers:
            if plate_number not in plate_check_set:
                plate_check_set.add(plate_number)
            else:
                duplicate_plates.add(plate_number)
        if duplicate_plates:
            raise ValidationError(
                key='plate_ranges', 
                msg='plates specified more than once: %s' 
                    % ', '.join([str(x) for x in duplicate_plates]))
        logger.info('plate_numbers: %r', plate_numbers)
        
        screen_facility_id = initializer_dict.pop('screen_facility_id', None)
        cherry_pick_request_id = \
            initializer_dict.pop('cherry_pick_request_id', None)
        if screen_facility_id is None and cherry_pick_request_id is None:
            msg='must provide one of: %r' % [
                'screen_facility_id','cherry_pick_request_id']
            raise ValidationError({
                'screen_facility_id': msg,
                'cherry_pick_request_id': msg 
            })
        screen = None
        if screen_facility_id is not None and cherry_pick_request_id is None:
            try:
                screen = Screen.objects.get(
                    facility_id=screen_facility_id)
                initializer_dict['screen'] = screen
            except ObjectDoesNotExist:
                raise Http404('screen_facility_id %r does not exist' 
                    % screen_facility_id)
        cpr = None
        if cherry_pick_request_id is not None: 
            try: 
                cpr = CherryPickRequest.objects.get(
                    cherry_pick_request_id=cherry_pick_request_id)
                initializer_dict['cherry_pick_request'] = cpr
            except ObjectDoesNotExist:
                raise Http404('cherry_pick_request_id does not exist: %s' 
                    % cherry_pick_request_id)
        
        try:
            rdt = RawDataTransform.objects.get(
                screen=screen, cherry_pick_request=cpr)
            rdt.rawdatainputfile_set.all().delete()
        except ObjectDoesNotExist:
            logger.info('create a new raw data transform')
            rdt = RawDataTransform.objects.create(
                screen=screen, cherry_pick_request=cpr)
        
        for key, val in initializer_dict.items():
            if hasattr(rdt, key):
                setattr(rdt, key, val)

        _meta = {}

        aps = rdt.assay_plate_size
        lps = rdt.library_plate_size
        
        if aps > lps:
            if len(plate_numbers) % 4 != 0:
                raise ValidationError(key='plate_ranges', 
                    msg='Must be a multiple of 4 if assay_plate_size==1536')
        
        # Retrieve the library well / Lab cherry pick well data:
        # - also record the user entered control well information
        rnai_columns = [
            'vendor_entrezgene_symbols','vendor_entrezgene_id',
            'vendor_genbank_accession_numbers','vendor_gene_name', 
            'facility_entrezgene_symbols','facility_entrezgene_id',
            'facility_genbank_accession_numbers','facility_gene_name', 
            'vendor_name','vendor_identifier','is_deprecated']
        rnai_columns_to_write = []
        wells = {}
        lcp_copywells = {}
        control_well_named_ranges = None
        if screen:
            if rdt.screen.screen_type == 'rnai':
                rnai_columns_to_write = rnai_columns
                # TODO: pool
            exact_well_fields = [
                'library_well_type','well_id','plate_number','well_name']
            exact_well_fields.extend(rnai_columns_to_write)
            logger.info('get well information: %r, fields: %r ...',
                plate_numbers, exact_well_fields)
            wells = self.get_reagent_resource()._get_list_response_internal(**{
                    'plate_number__in': plate_numbers,
                    'exact_fields': exact_well_fields
                })
            if not wells:
                raise ProgrammingError(
                    'no wells found for the plates: %r', plate_numbers)
            logger.info('retrieved %d wells for transformation', len(wells))
            logger.info('sample well keys: %r', wells[0].keys())
            wells = { well['well_id']: well for well in wells }
            control_well_named_ranges = self.get_control_wells(rdt)
            
            self.set_control_wells_to_library_wells(
                rdt, control_well_named_ranges, plate_numbers, wells)
            
        elif cpr:
            includes = [
                'vendor_batch_id','vendor_name','vendor_identifier',
                '-structure_image','-molfile','-library_plate_comment_array']
            if cpr.screen.screen_type == 'rnai':
                rnai_columns_to_write = rnai_columns
                includes.extend(rnai_columns_to_write)
            lcp_copywells = \
                self.get_labcherrypick_resource()\
                    ._get_list_response_internal(**{
                        'cherry_pick_request_id': cpr.cherry_pick_request_id,
                        'source_copy_name__is_null': False,
                        'includes': includes,
                    })
            if not lcp_copywells:
                raise ProgrammingError('No lab cherry picks found for cpr: %r', 
                    cpr.cherry_pick_request_id)
            logger.info('retrieved %d lcp_copywells for transformation', 
                len(lcp_copywells))
            logger.info('sample lcp well keys: %r', lcp_copywells[0].keys())
            lcp_copywells = { 
                lims_utils.well_id(
                    lcp['cherry_pick_plate_number'],lcp['destination_well']):lcp
                        for lcp in lcp_copywells }
            control_well_named_ranges = self.get_control_wells(rdt)
            
            
        # Read in the input matrices
        input_file_fields = schema['input_file_fields']
        logger.info('input_file_fields: %r', input_file_fields.keys())
        
        vocab_scope = input_file_fields['readout_type']['vocabulary_scope_ref']
        readout_vocab = self.get_vocab_resource()._get_vocabularies_by_scope(vocab_scope)
        if not readout_vocab:
            logger.warn('no vocabulary found for scope: %r, field: %r', 
                vocab_scope, 'readout_type')
        def read_input_matrices(rdt, ordinal, filekey, input_file):
            logger.info('read matrices for file: %d, %s', ordinal, filekey)
            logger.info('data file: %r', input_file)
            input_file_initializer = {
                'ordinal': ordinal
                }
            for key,val in param_hash.items():
                if key.find(filekey) == 0:
                    field_key = key[len(filekey)+1:]
                    if field_key in input_file_fields:
                        field = input_file_fields[field_key]
                        input_file_initializer[field_key] = \
                            parse_val(val, field_key, field['data_type'])
            logger.info('input_file_initializer: %r', input_file_initializer)

            errors = self.validate(input_file_initializer, 
                patch=False, 
                schema={ 
                    'id_attribute': ['ordinal'],
                    'fields': input_file_fields })
            if errors:
                raise ValidationError(key=filekey, msg=errors)

            filename = initializer_dict.get('filename', None)
            if filename is None:
                filename = input_file.name
                input_file_initializer['filename'] = filename
            rdif = RawDataInputFile(**input_file_initializer)
            
            (matrices, errors) = raw_data_reader.read(input_file, filename)
            if errors: 
                raise ValidationError(key=filekey, msg=errors)
            if not matrices:
                raise ValidationError(key=filekey, msg='no matrices read')
            
            parse_errors = raw_data_reader.parse_to_numeric(matrices)
            if parse_errors:
                raise ValidationError(key=filekey, msg=parse_errors)
            
            assay_plate_size_read = len(matrices[0])*len(matrices[0][0])
            if assay_plate_size_read != aps:
                msg = self.ERROR_MATRIX_SIZE_DETECTED % (
                    assay_plate_size_read, aps)
                raise ValidationError({
                    filekey: msg, 'assay_plate_size': msg })
            
            collation = Collation.get_value(rdif.collation_order)
            logger.info('read collation: %r', collation)
            conditions = re.split(r'[\s,]+', rdif.conditions) \
                if rdif.conditions else ['C1',]
            replicates = [chr(ord('A')+x) for x in range(0,rdif.replicates)]
            readouts = re.split(r'[\s,]+', rdif.readouts) \
                if rdif.readouts else ['read1',]
            logger.info('conditions: %r, readouts: %r, replicates: %r',
                conditions, readouts, replicates)

            # Determine the plates to read and validate relative sizes
            
            collation_count = len(conditions)*len(replicates)*len(readouts)
            transformed_matrix_count = len(matrices)
            if aps > lps:
                transformed_matrix_count *=4
            elif lps > aps:
                collation_count *= 4

            if collation_count > len(matrices):
                msg = self.ERROR_COLLATION_COUNT % (
                    collation_count, len(matrices))
                raise ValidationError(key=filekey, msg=msg)
            if len(matrices)%collation_count != 0:
                msg = self.ERROR_COLLATION_COUNT % (
                    collation_count, len(matrices))
                raise ValidationError(key=filekey, msg=msg)
            
            plates_required = transformed_matrix_count/collation_count
            logger.info('collation count: %d, plates required: %d',
                collation_count, plates_required)
            if plates_required != len(plate_numbers):
                logger.info(str((len(matrices), transformed_matrix_count,
                    collation_count, plates_required, len(plate_numbers))))
                msg = self.ERROR_PLATE_COUNT % (
                    plates_required, len(plate_numbers),
                    len(matrices), transformed_matrix_count,
                    collation_count)
                raise ValidationError({
                    filekey: msg, 'plate_ranges': msg })
            
            counter = plate_matrix_transformer.create_matrix_counter(
                collation, plate_numbers, conditions, replicates, readouts)
            logger.info('counter: %r', counter)

            if aps != lps:
                matrices = plate_matrix_transformer.transform(
                    matrices, counter, aps,lps)
                logger.info('transformed matrices: %d', len(matrices))

            rdif.raw_data_transform = rdt
            rdif.save()
            
            # cache readout_type title
            if rdif.readout_type in readout_vocab:
                rdif.readout_title = readout_vocab[rdif.readout_type]['title']
            else:
                logger.warn('vocab not found for readout_type: %r, vocabs: %r',
                    rdif.readout_type, readout_vocab)
                rdif.readout_title = rdif.readout_type
            return (matrices, rdif, counter)
        
        
        with  NamedTemporaryFile(
            delete=False, suffix='%s' % request.user.username) as temp_file:
            
            workbook = xlsxwriter.Workbook(temp_file) #, {'constant_memory': True})
            
            input_file_data = []
            for ordinal, filekey in enumerate(sorted(request.FILES.keys())):
                logger.info('read file: %d, %r', ordinal, filekey)
                input_file = request.FILES[filekey]
        
                (matrices, rdif, counter) = \
                    read_input_matrices(rdt, ordinal, filekey, input_file)
                input_file_data.append((matrices, rdif, counter))

                _matrix_read_meta = OrderedDict((
                    ('Ordinal', rdif.ordinal), 
                    ('Filename', rdif.filename), 
                    ('Collation', rdif.collation_order),
                    ('Readout Type', rdif.readout_title)
                ))
                _matrix_read_meta['Matrices'] = len(matrices)
                if aps > lps:
                    _matrix_read_meta['Matrices read (1536 well)'] = len(matrices)/4
                elif aps < lps:
                    _matrix_read_meta['Matrices read (96 well)'] = len(matrices)*4
                logger.info('read matrices: %d', len(matrices))
                
                for k,v in counter.counter_hash.items():
                    _matrix_read_meta[k.title() + 's'] = ', '.join([str(x) for x in v])
                _meta['File %d' % (rdif.ordinal+1)] = _matrix_read_meta
                logger.info('Raw data transform file read meta: %r', _matrix_read_meta)
                
            if screen:
                self.write_screen_xlsx(
                    rdt, plate_numbers, input_file_data, workbook, wells, rnai_columns_to_write)
            elif cpr:
                self.write_cpr_xlsx(
                    rdt, plate_numbers, input_file_data, workbook, lcp_copywells, 
                    control_well_named_ranges,rnai_columns_to_write)
            else:
                raise ProgrammingError('no screen or cpr')

                
            workbook.close()
            temp_file.close()
            
            rdt.temp_output_filename = temp_file.name
            rdt.save()
            logger.info('wrote temp file: %r', temp_file.name)
        
        _meta[API_MSG_RESULT] = 'success'
        
        return self.build_response(
            request, {API_RESULT_META: _meta }, response_class=HttpResponse, **kwargs)

    def get_control_wells(self, rdt):
        '''
        Read in the control well fields:
        - parse well selection input
        - transform the assay plate control wells to library plate format
        - set the control type abbreviation to well['type_abbreviation']
        - set the control label to well['control_label']
        collate errors: 
        - parsing
        - check for duplicated wells between types
        - check for control labels assigned to experimental wells
        '''
        
        logger.info('get_control_wells: %r', rdt)
        aps = rdt.assay_plate_size
        lps = rdt.library_plate_size
        combined_errors = defaultdict(dict)
        assay_control_named_ranges = {}
        control_well_fields = [
            'assay_positive_controls', 'assay_negative_controls',
            'assay_other_controls']
        for cfield in control_well_fields:
            (named_ranges, errors) = \
                lims_utils.parse_named_well_ranges(
                    getattr(rdt, cfield), aps)
            if errors:
                combined_errors[cfield][self.ERROR_CONTROL_PARSE] = errors
            assay_control_named_ranges[cfield] = named_ranges
            logger.info('control_well_field: %r, %r, named_well_ranges: %r',
                cfield, getattr(rdt, cfield), named_ranges)

        logger.info('assay control named ranges: %r', assay_control_named_ranges)

        (library_control_named_ranges, errors) = \
            lims_utils.parse_named_well_ranges(
                getattr(rdt, 'library_controls'), lps)
        if errors:
            error_hash = combined_errors['library_controls']
            error_hash[self.ERROR_CONTROL_PARSE] = errors
        logger.info('library_control_named_ranges: %r', 
            library_control_named_ranges)
        
        # Check for duplicated control wells between types:
        DUP_ERROR = self.ERROR_CONTROL_DUPLICATES
        # Assay well duplicates:
        assay_controls_flattened = defaultdict(set)
        for ctype, named_ranges in assay_control_named_ranges.items():
            for label, named_range in named_ranges.items():
                assay_controls_flattened[ctype].update(named_range['wells'])
        
        for test_ctype, test_wells in assay_controls_flattened.items():
            duplicate_wells = set()
            for ctype, assay_wells in assay_controls_flattened.items():
                if test_ctype == ctype:
                    continue
                duplicate_wells.update(
                    set(test_wells) & set(assay_wells))
            if duplicate_wells:
                combined_errors[test_ctype][DUP_ERROR] = duplicate_wells

        # check for library control well duplicates
        library_control_wells = set()
        for named_range in library_control_named_ranges.values():
            library_control_wells.update(named_range['wells'])
        
        if aps < lps:
            duplicate_wells = set()
            for assay_ctype, assay_wells in assay_controls_flattened.items():
                convoluted_wells = lims_utils.convolute_wells(
                    aps, lps, assay_wells)
                duplicates = library_control_wells & set(convoluted_wells)
                if duplicates:
                    duplicate_wells |= duplicates
                    combined_errors[assay_ctype][DUP_ERROR] = duplicates
            if duplicate_wells:
                combined_errors['library_controls'][DUP_ERROR] = duplicate_wells 
        elif lps < aps:
            duplicate_wells = set()
            for assay_ctype, assay_wells in assay_controls_flattened.items():
                deconvoluted_wells_in_quadrants = lims_utils.deconvolute_wells(
                    aps, lps, assay_wells)
                duplicates = library_control_wells & set(
                    [well for sublist in deconvoluted_wells_in_quadrants 
                        for well in sublist])
                if duplicates:
                    duplicate_wells |= duplicates
                    combined_errors[assay_ctype][DUP_ERROR] = duplicates
            if duplicate_wells:
                combined_errors['library_controls'][DUP_ERROR] = duplicate_wells 
        else:
            duplicate_wells = set()
            for assay_ctype, assay_wells in assay_controls_flattened.items():
                duplicates = library_control_wells & set(assay_wells)
                if duplicates:
                    duplicate_wells |= duplicates
                    combined_errors[assay_ctype][DUP_ERROR] = duplicates
            if duplicate_wells:
                combined_errors['library_controls'][DUP_ERROR] = duplicate_wells 
        
        if combined_errors:
            raise ValidationError(combined_errors)
        
        assay_control_named_ranges['library_controls'] = library_control_named_ranges
        
        return assay_control_named_ranges 
        
    def set_control_wells_to_library_wells(
        self, rdt, control_well_named_ranges, plates_in_order,wells):

        aps = rdt.assay_plate_size
        lps = rdt.library_plate_size
        
        # Associate the control well label and type with the wells
        combined_errors = defaultdict(dict)
        control_well_exceptions = {}
        
        library_control_named_ranges = control_well_named_ranges['library_controls']
        assay_control_named_ranges = { ctype: nr for ctype,nr in 
            control_well_named_ranges.items() if ctype != 'library_controls' }
        for well in wells.values():
            library_well_type = well.get('library_well_type')
            plate_number = well['plate_number']
            well_name = well['well_name']
            well['type_abbreviation'] = \
                self.library_well_type_abbreviations[library_well_type]
            
            if plate_number not in plates_in_order:
                # NOT a validation error, plate numbers were used to find wells
                raise ProgrammingError('plate not found: %r', plate_number)
            plate_index = plates_in_order.index(plate_number)
            
            assay_control_wellname = well_name
            if aps > lps:
                # convolute the well to (1536)
                quadrant = plate_index % 4
                assay_control_wellname = lims_utils.convolute_well(
                    lps,aps,well_name)[quadrant]
            elif aps < lps:
                # deconvolute the well to (96)
                assay_control_wellname = lims_utils.deconvolute_well(
                    lps,aps,well_name)
            
            library_control_exceptions = control_well_exceptions.setdefault(
                'library_controls',defaultdict(set))
            found = False
            for label, named_range in library_control_named_ranges.items():
                if well_name in named_range['wells']:
                    well['control_label'] = label
                    found = True
                    if library_well_type == 'experimental':
                        library_control_exceptions[well_name].add(
                            well['well_id'])
            
            if found is False:
                for ctype, named_ranges in assay_control_named_ranges.items():
                    ctype_exceptions = \
                        control_well_exceptions.setdefault(
                            ctype,defaultdict(set))
                    for label, named_range in named_ranges.items():
                        if assay_control_wellname in named_range['wells']:
                            well['type_abbreviation'] = \
                                self.control_type_abbreviations[ctype]
                            well['control_label'] = label
                            if library_well_type == 'experimental':
                                ctype_exceptions[assay_control_wellname].add(
                                    well['well_id'])
                            break
        
        logger.info('control_well_exceptions: %r', control_well_exceptions)    
        for ctype, well_exception_dict in control_well_exceptions.items():
            if not well_exception_dict:
                continue
            well_exception_dict = { wellname:','.join(well_ids) 
                for wellname, well_ids in well_exception_dict.items() }
            combined_errors[ctype][self.ERROR_CONTROL_WELL_TYPE] = \
                well_exception_dict
        
        if combined_errors:
            raise ValidationError(combined_errors)
        
    def write_cpr_xlsx(
        self, rdt, plate_numbers, input_file_data, workbook, lcp_well_hash,
        control_well_named_ranges,rnai_columns_to_write ):
        '''
        Write the plate matrices directly to a spreadsheet, in collation order:
        - merge in lab_cherry_pick data.
        '''
        DEBUG = False or logger.isEnabledFor(logging.DEBUG)
        logger.info('write cpr worksheet for %r...', rdt)

        cpr_id = rdt.cherry_pick_request.cherry_pick_request_id
        aps = rdt.assay_plate_size
        lps = rdt.library_plate_size
        
        control_well_hash = {}
        for ctype, named_ranges in control_well_named_ranges.items():
            for label,named_range in named_ranges.items():
                for wellname in named_range['wells']:
                    control_well_hash[wellname] = (ctype, label)
        
        headers = ['Plate', 'Well', 'Type']
        if aps > lps:
            headers = ['Plate','Well','Source Plate','Quadrant','Source Well','Type']
        elif lps > aps:
            headers = ['Plate','Well','Quadrant','Source Well','Type']

        row_size = lims_utils.get_rows(lps)
        col_size = lims_utils.get_cols(lps)
        logger.info('row size: %d, col size: %d', row_size, col_size)
        
        # NOTE: Even with "all_plates_in_single_worksheet" option, each input 
        # file must be in its own sheet, because collations may differ
        single_sheet = None
        sheet_name = None
        quadrant = None
        source_plates = None
        cumulative_output_row = 0
        if rdt.output_sheet_option \
                == 'all_plates_in_single_worksheet':
            sheet_name = 'data_%d' % (rdif.ordinal+1)
            single_sheet = workbook.add_worksheet(sheet_name)
        
        for plate_index, plate in enumerate(plate_numbers):
            cpr_plate = 'CP%d_%d' % (cpr_id, plate)
            new_sheet_name = 'CP%d_%d' % (cpr_id, plate)
            if aps > lps:
                quadrant = plate_index % 4
                source_plates = plate_numbers[(plate_index/4):(plate_index/4+4)]
                new_sheet_name = 'CP%d_%s' % (
                    cpr_id, ','.join([str(p) for p in source_plates]))
            elif lps > aps:
                pass
            if single_sheet:
                sheet = single_sheet
            elif sheet_name != new_sheet_name:
                sheet_name = new_sheet_name
                sheet = workbook.add_worksheet(sheet_name)
            logger.info('write sheet: %r', sheet_name)

            # Write the header row for the plate:
            for i, header in enumerate(headers):
                sheet.write_string(0,i,header)
            
            current_col = len(headers)
            for (matrices, rdif, counter) in input_file_data:
                for i,counter_readout in enumerate(counter.iterate_counter_columns()): 
                    collation_string = \
                        '{readout}_{condition}_{replicate}'.format(
                            **counter_readout).title()
                    if rdif.readout_title:
                        collation_string = '%s_%s' % (rdif.readout_title, collation_string)
                    else:
                        collation_string = '%s_%s' % (rdif.readout_type, collation_string)
                        
                    if DEBUG:
                        logger.info('write collation_string: col: %d, %s', 
                            current_col, collation_string)
                    sheet.write_string(0,current_col+i,collation_string )
                current_col += i + 1
    
            sheet.write_string(0,current_col,'Pre-Loaded Controls')
            current_col += 1
            sheet.write_string(0,current_col,'Library Plate')
            current_col += 1
            sheet.write_string(0,current_col,'Source Well')
            current_col += 1
            sheet.write_string(0,current_col,'Library Name')
            current_col += 1
            sheet.write_string(0,current_col,'Vendor ID')
            current_col += 1
            sheet.write_string(0,current_col,'Vendor Batch ID')
                
            if rnai_columns_to_write:
                current_col += 1
                sheet.write_string(0,current_col,'Gene Symbol')
                current_col += 1
                sheet.write_string(0,current_col,'Gene IDs')
                current_col += 1
                sheet.write_string(0,current_col,'Genbank Accession Nos')
                current_col += 1
                sheet.write_string(0,current_col,'Gene Names')
                    
            # Write the values for the plate:
            
            # NOTE: to support xlsxwriter 'constant_memory': True - optimized write, 
            # rows must be written sequentially        
            for rownum in range(0,row_size):
                for colnum in range(0, col_size):
                    output_row = 1 + cumulative_output_row + colnum + rownum * col_size
                    current_col = 0      
                    
                    wellname = lims_utils.get_well_name(rownum, colnum)
                    well_id = lims_utils.well_id(plate, wellname)
                    lcp_well = lcp_well_hash.get(well_id, None)
                    if DEBUG:
                        logger.info('write row: %d: %r', output_row, well_id)
                        logger.info('found lcp well: %r', lcp_well)
                    
                    sheet.write_string(output_row,current_col,cpr_plate)
                    current_col += 1
                    sheet.write_string(output_row,current_col,wellname)
                    current_col += 1
                    
                    assay_plate_wellname = wellname
                    if aps > lps:
                        assay_plate_wellname = \
                            lims_utils.convolute_well(
                                lps, aps, wellname)[quadrant]
                        sheet.write_string(output_row,current_col,
                            ','.join([str(p) for p in source_plates]))
                        current_col += 1
                        sheet.write_string(output_row,current_col,str(quadrant+1))
                        current_col += 1
                        sheet.write_string(output_row,current_col,assay_plate_wellname);
                        current_col += 1
                    elif lps > aps:
                        assay_plate_wellname = \
                            lims_utils.deconvolute_well(lps, aps, wellname)
                        (row,col) = lims_utils.well_row_col(wellname)
                        source_quadrant = lims_utils.deconvolute_quadrant(lps, aps, row, col)
                        sheet.write_string(output_row,current_col,str(source_quadrant+1))
                        current_col += 1
                        sheet.write_string(output_row,current_col,assay_plate_wellname);
                        current_col += 1
                    
                    control = control_well_hash.get(assay_plate_wellname, None)
                    
                    if control:
                        abbrev = self.control_type_abbreviations[control[0]]
                        sheet.write_string(output_row,current_col,abbrev)
                    elif lcp_well:
                        sheet.write_string(output_row,current_col,'X')
                    else:
                        sheet.write_string(output_row,current_col,'E')
                    current_col += 1

                    for j,(matrices, rdif, counter) in enumerate(input_file_data):
                        for i,counter_readout in enumerate(counter.iterate_counter_columns()): 
                            matrix_index = counter.get_index(dict(counter_readout, plate=plate))
                            matrix = matrices[matrix_index]
                            val = matrix[rownum][colnum]
                            if DEBUG:
                                logger.info('write output_row: %d, col: %d,  val: %r', 
                                    output_row, current_col+i, str(val))
                            sheet.write_number(output_row, current_col+i, val)
                        current_col += i +1

                    if control:
                        sheet.write_string(output_row,current_col,control[1])
                    elif lcp_well:    
                        current_col += 1
                        sheet.write_string(output_row,current_col,
                            str(lcp_well.get('library_plate')))
                        current_col += 1
                        sheet.write_string(output_row,current_col,
                            lcp_well.get('source_well_name'))
                        current_col += 1
                        sheet.write_string(output_row,current_col,
                            lcp_well.get('library_short_name'))
                        current_col += 1
                        vals = [lcp_well.get('vendor_name'),
                            lcp_well.get('vendor_identifier')]
                        vals = [ str(v) for v in vals if v ]
                        sheet.write_string(output_row,current_col,' '.join(vals))
                        current_col += 1
                        val = lcp_well.get('vendor_batch_id',None)
                        if val:
                            sheet.write_string(output_row,current_col,)

                        if rnai_columns_to_write:
                            current_col += 1
                            vals = [lcp_well.get('vendor_entrezgene_symbols'),
                                    lcp_well.get('facility_entrezgene_symbols')]
                            vals = sorted(set(
                                [ v for vs in vals if vs for v in vs if v ]))
                            sheet.write_string(output_row,current_col,', '.join(vals))
                            
                            current_col += 1
                            vals = [lcp_well.get('vendor_entrezgene_id'),
                                    lcp_well.get('facility_entrezgene_id')]
                            vals = sorted(set([ str(v) for v in vals if v ]))
                            sheet.write_string(output_row,current_col,', '.join(vals))
                            
                            current_col += 1
                            vals = [lcp_well.get('vendor_genbank_accession_numbers'),
                                    lcp_well.get('facility_genbank_accession_numbers')]
                            vals = sorted(set(
                                [ str(v) for vs in vals if vs for v in vs if v ]))
                            sheet.write_string(output_row,current_col,', '.join(vals))
                            
                            current_col += 1
                            vals = [lcp_well.get('vendor_gene_name'),
                                    lcp_well.get('facility_gene_name')]
                            vals = sorted(set([ str(v) for v in vals if v ]))
                            sheet.write_string(output_row,current_col,', '.join(vals))
            if single_sheet:
                cumulative_output_row = output_row                
            if aps > lps and quadrant < 4:
                cumulative_output_row = output_row                


    def write_screen_xlsx(
        self, rdt, plate_numbers, input_file_data, workbook, library_well_hash,
        rnai_columns_to_write):
        '''
        Write the plate matrices directly to a spreadsheet, in collation order:
        - merge in library well data.
        
        '''
        DEBUG = False or logger.isEnabledFor(logging.DEBUG)

        aps = rdt.assay_plate_size
        lps = rdt.library_plate_size

        logger.info('write workbook...')
        
        headers = ['Plate', 'Well', 'Type','Exclude']
        if aps > lps:
            headers = ['Plate','Well','Source Plate','Quadrant','Source Well','Type','Exclude']
        elif lps > aps:
            headers = ['Plate','Well','Quadrant','Source Well','Type','Exclude']

        row_size = lims_utils.get_rows(lps)
        col_size = lims_utils.get_cols(lps)
        logger.info('row size: %d, col size: %d', row_size, col_size)

        single_sheet = None
        sheet_name = None
        quadrant = None
        source_plates = None
        cumulative_output_row = 0
        if rdt.output_sheet_option == 'all_plates_in_single_worksheet':
            sheet_name = 'data'
            single_sheet = workbook.add_worksheet(sheet_name)
        
        for plate_index, plate in enumerate(plate_numbers):
            new_sheet_name = str(plate)
            if aps > lps:
                quadrant = plate_index % 4
                source_plates = plate_numbers[(plate_index/4):(plate_index/4+4)]
                new_sheet_name = ','.join([str(p) for p in source_plates])
            elif lps > aps:
                pass
            if single_sheet:
                sheet = single_sheet
            elif sheet_name != new_sheet_name:
                sheet_name = new_sheet_name
                sheet = workbook.add_worksheet(sheet_name)
            logger.info('write sheet: %r', sheet_name)

            # Write the header row for the plate:
            for i, header in enumerate(headers):
                sheet.write_string(0,i,header)
                        
            current_col = len(headers)
            for (matrices, rdif, counter) in input_file_data:
                for i,counter_readout in enumerate(counter.iterate_counter_columns()): 
                    collation_string = \
                        '{readout}_{condition}_{replicate}'.format(
                            **counter_readout).title()
                    if rdif.readout_title:
                        collation_string = '%s_%s' % (rdif.readout_title, collation_string)
                    else:
                        collation_string = '%s_%s' % (rdif.readout_type, collation_string)
                    if DEBUG:
                        logger.info('write collation_string: col: %d, %s', 
                            current_col, collation_string)
                    sheet.write_string(0,current_col+i,collation_string )
                current_col += i + 1
    
            sheet.write_string(0,current_col,'Pre-Loaded Controls')
            if rnai_columns_to_write:
                current_col += 1
                sheet.write_string(0,current_col,'Entrezgene Symbol')
                current_col += 1
                sheet.write_string(0,current_col,'Entrezgene ID')
                current_col += 1
                sheet.write_string(0,current_col,'Genbank Accession No.')
                current_col += 1
                sheet.write_string(0,current_col,'Catalog No.')
                current_col += 1
                sheet.write_string(0,current_col,'Gene Name')
                current_col += 1
                sheet.write_string(0,current_col,'Deprecated Pool')
                
            # NOTE: to support xlsxwriter 'constant_memory': True - optimized write, 
            # rows must be written sequentially  
            for rownum in range(0,row_size):
                for colnum in range(0, col_size):
                    output_row = 1 + cumulative_output_row + colnum + rownum * col_size
                    current_col = 0      
                    
                    wellname = lims_utils.get_well_name(rownum, colnum)
                    well_id = lims_utils.well_id(plate, wellname)
                    
                    well = library_well_hash[well_id]
                    if DEBUG:
                        logger.info('write row: %d: %r', output_row, wellname)
                    sheet.write_string(output_row,current_col,str(plate))
                    current_col += 1
                    sheet.write_string(output_row,current_col,wellname)
                    current_col += 1
                    
                    if aps > lps:
                        sheet.write_string(output_row,current_col,
                            ','.join([str(p) for p in source_plates]))
                        current_col += 1
                        source_wells = lims_utils.convolute_well(lps, aps, wellname)
                        sheet.write_string(output_row,current_col,str(quadrant+1))
                        current_col += 1
                        sheet.write_string(output_row,current_col,source_wells[quadrant]);
                        current_col += 1
                    elif lps > aps:
                        (row,col) = lims_utils.well_row_col(wellname)
                        source_quadrant = lims_utils.deconvolute_quadrant(lps, aps, row, col)
                        source_wellname = lims_utils.deconvolute_well(lps, aps, wellname)
                        sheet.write_string(output_row,current_col,str(source_quadrant+1))
                        current_col += 1
                        sheet.write_string(output_row,current_col,source_wellname);
                        current_col += 1
                    
                    sheet.write_string(
                        output_row,current_col,well.get('type_abbreviation',None))
                    current_col += 1
                    # NOP sheet.write_string(output_row,current_col, 'exclude')
                    current_col += 1
                    
                    for (matrices, rdif, counter) in input_file_data:
                        for i,counter_readout in enumerate(counter.iterate_counter_columns()): 
                            matrix_index = counter.get_index(dict(counter_readout, plate=plate))
                            matrix = matrices[matrix_index]
                            val = matrix[rownum][colnum]
                            if DEBUG:
                                logger.info('write output_row: %d, col: %d,  val: %r', 
                                    output_row, current_col+i, str(val))
                            sheet.write_number(output_row, current_col+i, val)
                        current_col += i +1

                    control_label = well.get('control_label', None)
                    if control_label:
                        sheet.write_string(output_row,current_col,control_label)
                    elif rnai_columns_to_write:
                        current_col += 1
                        vals = [well.get('vendor_entrezgene_symbols'),
                                well.get('facility_entrezgene_symbols')]
                        vals = sorted(set(
                            [ v for vs in vals if vs for v in vs if v ]))
                        sheet.write_string(output_row,current_col,', '.join(vals))
                        
                        current_col += 1
                        vals = [well.get('vendor_entrezgene_id'),
                                well.get('facility_entrezgene_id')]
                        vals = sorted(set([ str(v) for v in vals if v ]))
                        sheet.write_string(output_row,current_col,', '.join(vals))
                        
                        current_col += 1
                        vals = [well.get('vendor_genbank_accession_numbers'),
                                well.get('facility_genbank_accession_numbers')]
                        vals = sorted(set(
                            [str(v) for vs in vals if vs for v in vs if v ]))
                        sheet.write_string(output_row,current_col,', '.join(vals))
                        
                        current_col += 1
                        vals = [well.get('vendor_name'),well.get('vendor_identifier')]
                        vals = [ str(v) for v in vals if v ]
                        sheet.write_string(output_row,current_col,' '.join(vals))

                        current_col += 1
                        vals = [well.get('vendor_gene_name'),
                                well.get('facility_gene_name')]
                        vals = sorted(set([ str(v) for v in vals if v ]))
                        sheet.write_string(output_row,current_col,', '.join(vals))

                        current_col += 1
                        sheet.write_string(output_row,current_col,
                            'Y' if well.get('is_deprecated') is True else '')

            if single_sheet:
                cumulative_output_row = output_row                
            if aps > lps and quadrant < 4:
                cumulative_output_row = output_row                
        

class ScreenResource(DbApiResource):
    
    class Meta:

        queryset = Screen.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'screen'
        authorization = ScreenAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 
        
    def __init__(self, **kwargs):
        super(ScreenResource, self).__init__(**kwargs)
        self.apilog_resource = None
        self.cpr_resource = None
        self.activity_resource = None
        self.lcp_resource = None
        self.libraryscreening_resource = None
        self.screenresult_resource = None
        self.su_resource = None
         
    def clear_cache(self):
        logger.info('clear screen caches')
        DbApiResource.clear_cache(self)
        caches['screen'].clear()
        self.get_screenresult_resource().clear_cache()
        
    def get_su_resource(self):
        if self.su_resource is None:
            self.su_resource = ScreensaverUserResource()
        return self.su_resource
        
    def get_screenresult_resource(self):
        if self.screenresult_resource is None:
            self.screenresult_resource = ScreenResultResource()
        return self.screenresult_resource
    
    def get_librarycopyplate_resource(self):
        if self.lcp_resource is None:
            self.lcp_resource = LibraryCopyPlateResource()
        return self.lcp_resource
    
    def get_apilog_resource(self):
        if self.apilog_resource is None:
            self.apilog_resource = ApiLogResource()
        return self.apilog_resource
    
    def get_cpr_resource(self):
        if self.cpr_resource is None:
            self.cpr_resource = CherryPickRequestResource()
        return self.cpr_resource
    
    def get_activity_resource(self):
        if self.activity_resource is None:
            self.activity_resource = ActivityResource()
        return self.activity_resource
    
    def get_library_screening_resource(self):
        if self.libraryscreening_resource is None:
            self.libraryscreening_resource = LibraryScreeningResource()
        return self.libraryscreening_resource
        
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/ui%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_detail_uiview'), 
                name="api_dispatch_screen_detail_uiview"),

            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/libraries%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_libraryview'),
                name="api_dispatch_screen_libraryview"),

            # FIXME: screen/{facility_id}/cherrypickrequest is the canonical form
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/cherrypicks%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_cherrypickview'),
                name="api_dispatch_screen_cherrypickview"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/cherrypickrequest%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_cherrypickview'),
                name="api_dispatch_screen_cherrypickrequest"),
                
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/plates_screened%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_plates_screened_view'),
                name="api_dispatch_plates_screened_view"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/copyplatesloaded%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_lcp_loadedview'),
                name="api_dispatch_screen_lcp_loadedview"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/billing%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_billingview'),
                name="api_dispatch_screen_billingview"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/datacolumns%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_datacolumnview'),
                name="api_dispatch_screen_datacolumnview"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/activities%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_activityview'),
                name="api_dispatch_screen_activityview"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/libraryscreening%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_screeningview'),
                name="api_dispatch_screen_screeningview"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/screenresult%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_screenresultview'),
                name="api_dispatch_screen_screenresultview"),
            url(r"^(?P<resource_name>%s)/(?P<facility_id>([\w]+))/attachedfiles%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_attachedfileview'),
                name="api_dispatch_screen_attachedfileview"),
            url((r"^(?P<resource_name>%s)/(?P<facility_id>([\w]+))"
                 r"/attachedfiles/(?P<attached_file_id>([\d]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_attachedfiledetailview'),
                name="api_dispatch_screen_attachedfiledetailview"),
            url(r"^(?P<resource_name>%s)/(?P<facility_id>([\w]+))/publications%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_publicationview'),
                name="api_dispatch_screen_publicationview"),
            url((r"^(?P<resource_name>%s)/(?P<facility_id>([\w]+))"
                 r"/publications/(?P<publication_id>([\d]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_publicationdetailview'),
                name="api_dispatch_screen_publicationdetailview"),
            url((r"^(?P<resource_name>%s)/(?P<facility_id>([\w]+))"
                 r"/plate_range_search%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_plate_range_search_view'), 
                name="api_dispatch_plate_range_search_view"),
            url((r"^(?P<resource_name>%s)/(?P<facility_id>([\w]+))"
                 r"/plate_range_search/(?P<activity_id>(\d+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_plate_range_search_view'), 
                name="api_dispatch_plate_range_search_view"),
        ]    

    def dispatch_plate_range_search_view(self, request, **kwargs):
        ''' 
        Find: 
        - plates already asssociated with the libraryscreenings for the screen
        - plates matched by the "raw_search_data"
        Note: bypasses the "dispatch" framework call
        -- must be authenticated and authorized
        '''
        logger.info('screen plate range search view')
        return self.get_library_screening_resource().\
            dispatch_plate_range_search_view(request, **kwargs)
    
    def build_schema(self, user=None):
        
        logger.info('build screen schema for user: %r', user)
        schema = DbApiResource.build_schema(self, user=user)
        
        if self._meta.authorization.is_restricted_view(user):
            if DEBUG_SCREEN_ACCESS:
                logger.info(
                    'Screen schema: Remove sort and search capability for '
                    'restricted access fields')
            restricted_fields = set()
            for key,field in schema['fields'].items():
                if 'data_access_level' in field:
                    if field['data_access_level'] != 0:
                        field['filtering'] = False
                        field['ordering'] = False
                        restricted_fields.add(key)
            if DEBUG_SCREEN_ACCESS:
                logger.info('schema restricted_fields: %r', restricted_fields)
        return schema
        
    def dispatch_screen_detail_uiview(self, request, **kwargs):
        ''' 
        Special method to populate nested entities for the UI 
        - bypasses the "dispatch" framework call
        -- must be authenticated and authorized
        '''
        logger.info('dispatch_screen_detail_uiview')
        
        facility_id = kwargs.get('facility_id', None)
        if not facility_id:
            raise NotImplementedError('must provide a facility_id parameter')
        
        self.is_authenticated(request)

        resource_name = kwargs.pop('resource_name', self._meta.resource_name)
        authorized = self._meta.authorization._is_resource_authorized(
            request.user, 'read', **kwargs)
        if authorized:
            authorized = self._meta.authorization.has_screen_read_authorization(
                request.user, facility_id)
        if authorized is not True:
            raise PermissionDenied
        
        is_restricted_view = self._meta.authorization.is_restricted_view(
            request.user)
        logger.info('is_restricted: %r: %r: %r', resource_name, request.user, 
            is_restricted_view)
                
        if request.method.lower() != 'get':
            return self.dispatch('detail', request, **kwargs)
        
        cache_key = 'detail_ui_%s_%s' % (facility_id, request.user.username) 
        screen_cache = caches['screen']
        _data = screen_cache.get(cache_key)
        
        if not _data:
            logger.info('cache key not set: %s', cache_key)
            _data = self._get_detail_response_internal(
                user=request.user, **kwargs)
            if not _data:
                raise Http404('no screen found for facility_id: %r' % facility_id)
            else:
                # response = self.dispatch('detail', request, format='json', **kwargs )
                # if response.status_code == 200:
                #     _data = self._meta.serializer.deserialize(
                #         JSON_MIMETYPE, 
                #         response['Content-Type'])
                if _data.get('user_access_level_granted') == 3:
                    logger.info('retrieve status data...')
                    _status_data = \
                        self.get_apilog_resource()._get_list_response_internal(**{
                            'key': _data['facility_id'],
                            'ref_resource_name': 'screen',
                            'diff_keys': 'status',
                            'order_by': ['-date_time'],
                            })
                    _data['status_data'] = _status_data
                
                    _cpr_data = \
                        self.get_cpr_resource()._get_list_response_internal(**{
                            'limit': 0,
                            'screen_facility_id__eq': _data['facility_id'],
                            'order_by': ['-date_requested'],
                            'exact_fields': [
                                'cherry_pick_request_id','date_requested', 
                                'requested_by_name'],
                            })
                    _data['cherry_pick_request_data'] = _cpr_data
                    
                    _latest_activities_data = \
                        self.get_activity_resource()._get_list_response_internal(**{
                            'screen_facility_id__eq': _data['facility_id'],
                            'limit': 1,
                            'order_by': ['-date_of_activity'],
                            'exact_fields': ['activity_id','type', 
                                'performed_by_name'],
                            })
                    _data['latest_activities_data'] = _latest_activities_data
                    
                    # TODO: attached files
                    # TODO: publications
                    screen_cache.set(cache_key, _data)
                else:
                    # do not cache if restricted
                    _data['is_restricted_view'] = True
        else:
            logger.info('cache key set: %s', cache_key)
              
        # Serialize
        # FIXME: refactor to generalize serialization:
        # see build_response method
        content_type = self.get_serializer().get_accept_content_type(
            request,format=kwargs.get('format', None))
        if content_type in [XLS_MIMETYPE,CSV_MIMETYPE]:
            _data = {'objects': [_data]}
            filename = 'screen_detail_%s' % facility_id
        response = HttpResponse(
            content=self.get_serializer().serialize(
                _data, content_type),
            content_type=content_type)
        if content_type == XLS_MIMETYPE:
            response['Content-Disposition'] = \
                'attachment; filename=%s.xlsx' % filename
        if content_type == CSV_MIMETYPE:
            response['Content-Disposition'] = \
                'attachment; filename=%s.csv' % filename
        downloadID = request.GET.get('downloadID', None)
        if downloadID:
            logger.info('set cookie "downloadID" %r', downloadID )
            response.set_cookie('downloadID', downloadID)
        else:
            logger.debug('no downloadID: %s' % request.GET )
    
        return response;
        
    def dispatch_screen_attachedfileview(self, request, **kwargs):
        kwargs['screen_facility_id'] = kwargs.pop('facility_id')
        method = 'list'
        if request.method.lower() == 'post':
            # if post is used, force to "post_detail"
            method = 'detail'
        return AttachedFileResource().dispatch(method, request, **kwargs)    

    def dispatch_screen_attachedfiledetailview(self, request, **kwargs):
        kwargs['screen_facility_id'] = kwargs.pop('facility_id')
        return AttachedFileResource().dispatch('detail', request, **kwargs)    
                
    def dispatch_screen_publicationview(self, request, **kwargs):
        kwargs['screen_facility_id'] = kwargs.pop('facility_id')
        method = 'list'
        if request.method.lower() == 'post':
            # if post is used, force to "post_detail"
            method = 'detail'
        return PublicationResource().dispatch(method, request, **kwargs)    

    def dispatch_screen_publicationdetailview(self, request, **kwargs):
        kwargs['screen_facility_id'] = kwargs.pop('facility_id')
        return PublicationResource().dispatch('detail', request, **kwargs)    
                
    def dispatch_screen_activityview(self, request, **kwargs):
        kwargs['screen_facility_id__eq'] = kwargs.pop('facility_id')
        return self.get_activity_resource().dispatch('list', request, **kwargs)    
    
    def dispatch_screen_screeningview(self, request, **kwargs):
        kwargs['screen_facility_id__eq'] = kwargs.pop('facility_id')
        return self.get_library_screening_resource().dispatch('list', request, **kwargs)    
    
    def dispatch_screen_screenresultview(self, request, **kwargs):
        logger.info('dispatch screenresultview...')
        kwargs['screen_facility_id'] = kwargs.pop('facility_id')
        return ScreenResultResource().dispatch('list', request, **kwargs)    
    
    def dispatch_screen_datacolumnview(self, request, **kwargs):
        kwargs['screen_facility_id__eq'] = kwargs.pop('facility_id')
        return DataColumnResource().dispatch('list', request, **kwargs)    
        
    def dispatch_screen_cherrypickview(self, request, **kwargs):
        kwargs['screen_facility_id'] = kwargs.pop('facility_id')
        return self.get_cpr_resource().dispatch('list', request, **kwargs)    
        
    def dispatch_screen_libraryview(self, request, **kwargs):
        kwargs['for_screen_facility_id'] = kwargs.pop('facility_id')
        # NOTE: authorization is checked on LibraryResource
        return LibraryResource().dispatch('list', request, **kwargs)    

    def dispatch_plates_screened_view(self, request, **kwargs):

        screen_facility_id = kwargs.pop('facility_id')
        # NOTE: authorization is performed in LibraryCopyPlateResource
        return self.get_librarycopyplate_resource()\
            .build_screened_plate_response(request, screen_facility_id, **kwargs)    

    def dispatch_screen_billingview(self, request, **kwargs):
        kwargs['visibilities'] = 'billing'
        return self.dispatch('detail', request, **kwargs)    

    @read_authorization
    def get_detail(self, request, **kwargs):
        facility_id = kwargs.get('facility_id', None)
        if not facility_id:
            raise NotImplementedError('must provide a facility_id parameter')
        
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])

        return self.build_list_response(request, **kwargs)

    def get_query(self, request, schema, param_hash):
        
        DEBUG_SCREEN = False or logger.isEnabledFor(logging.DEBUG)
        # screens_for_username = param_hash.get('screens_for_username', None)
        screens_for_userid = param_hash.get('screens_for_userid', None)
        # general setup

        facility_id = param_hash.pop('facility_id', None)
        if facility_id:
            param_hash['facility_id__eq'] = facility_id
        
        manual_field_includes = set(param_hash.get('includes', []))
        # for joins
        manual_field_includes.add('screen_id')
        manual_field_includes.add('has_screen_result')
        manual_field_includes.add('data_sharing_level')
        
        extra_params = {}
        if screens_for_userid:
            screener_role_cte = ScreenResource.get_screener_role_cte().cte(
                'screener_roles1')
            manual_field_includes.add('screensaver_user_role')
            extra_params['user'] = screens_for_userid
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, schema, **extra_params)
        
        filter_expression = \
            self._meta.authorization.filter(request.user,filter_expression)

        order_params = param_hash.get('order_by', [])
        order_params.append('-facility_id')
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
        logger.debug('visible fields: %r', field_hash.keys())
        
        # specific setup
        base_query_tables = ['screen', 'screen_result']
        _screen = self.bridge['screen']
        _child_screen = _screen.alias('child_screen')
        _parent_screen = _screen.alias('ps')
        _screen_result = self.bridge['screen_result']
        _ap = self.bridge['assay_plate']
        _aw = self.bridge['assay_well']
        _library = self.bridge['library']
        _copy = self.bridge['copy']
        _plate = self.bridge['plate']
        _cpr = self.bridge['cherry_pick_request']
        _lcp = self.bridge['lab_cherry_pick']
        _cpap = self.bridge['cherry_pick_assay_plate']
        _cplt = self.bridge['cherry_pick_liquid_transfer']
        _sfs = self.bridge['screen_funding_supports']
        _screen_collaborators = self.bridge['screen_collaborators']
        _su = self.bridge['screensaver_user']
        _lhsu = _su.alias('lhsu')
        _user_cte = ScreensaverUserResource.get_user_cte().cte('users_s1')
        _collaborator = _user_cte.alias('collaborator')
        _activity = self.bridge['activity']
        _screen_keyword = self.bridge['screen_keyword']
        _screen_cell_lines = self.bridge['screen_cell_lines']
        _library_screening = self.bridge['library_screening']
        _cp_screening = self.bridge['cherry_pick_screening']
        _lab_activity = self.bridge['lab_activity']
        _service_activity = self.bridge['service_activity']
        _publication = self.bridge['publication']
        _attached_file = self.bridge['attached_file']
        _dc = self.bridge['data_column']
        _overlap_screens = self.get_create_screen_overlap_indexes()
        _overlap_screen = _screen.alias('overlap_screen')
        # create CTEs -  Common Table Expressions for the intensive queries:
        
        collaborators = (
            select([
                _screen_collaborators.c.screen_id,
                _collaborator.c.name,
                _collaborator.c.screensaver_user_id,
                _collaborator.c.username,
                _collaborator.c.email,
                _concat(
                    _collaborator.c.name, '<', _collaborator.c.email, '>'
                    ).label('fullname')])
            .select_from(_collaborator.join(
                _screen_collaborators,
                _collaborator.c.screensaver_user_id 
                    == _screen_collaborators.c.screensaveruser_id))
            .order_by(_collaborator.c.username))
        collaborators = collaborators.cte('collaborators')

        pin_transfer_approval = (
            select([
                _screen.c.screen_id,
                _user_cte.c.username,_user_cte.c.name,
                _activity.c.date_of_activity, _activity.c.comments])
            .select_from(
                _activity.join(
                    _user_cte, 
                    _activity.c.performed_by_id==_user_cte.c.screensaver_user_id)
                .join(
                    _screen, 
                    _activity.c.activity_id==_screen.c.pin_transfer_admin_activity_id))
            )
        pin_transfer_approval = pin_transfer_approval.cte('pta')

        # create a cte for the max screened replicates_per_assay_plate
        # - cross join version:
        # select 
        # ap.screen_id, ap.plate_number, 
        # ap.replicate_ordinal,lesser.replicate_ordinal  
        # from assay_plate ap 
        # left outer join assay_plate lesser 
        # on ap.plate_number=lesser.plate_number 
        # and ap.screen_id=lesser.screen_id 
        # and lesser.replicate_ordinal > ap.replicate_ordinal  
        # where lesser.replicate_ordinal is null;
        # - aggregate version:
        # select 
        # ap.screen_id, 
        # ap.plate_number, 
        # max(replicate_ordinal) as max_ordinal 
        # from assay_plate ap 
        # group by screen_id, plate_number
        # order by screen_id, plate_number            
        
        aps = select([_ap.c.screen_id, func.count(1).label('count')]).\
            select_from(_ap).\
            where(_ap.c.library_screening_id != None).\
            group_by(_ap.c.screen_id).\
            order_by(_ap.c.screen_id)
        aps = aps.cte('aps')

        # Removed - 20160408 
        # per discussion with JS, no need to try to link presumed assay_plates
        # from the screen_result to (screening visit) assay_plates;
        # Instead, we can infer replicate loading from the screen_result import
        # where the data_column shows the replicate
        # apdl = select([_ap.c.screen_id, func.count(1).label('count')]).\
        #     select_from(_ap).\
        #     where(_ap.c.screen_result_data_loading_id != None).\
        #     group_by(_ap.c.screen_id).\
        #     order_by(_ap.c.screen_id)
        # apdl = apdl.cte('apdl')
        
        # create a cte for the max screened replicates_per_assay_plate
        apsrc = select([
            _ap.c.screen_id,
            _ap.c.plate_id,
            func.max(_ap.c.replicate_ordinal).label('max_per_plate') ]).\
                select_from(_ap).\
                group_by(_ap.c.screen_id, _ap.c.plate_id).\
                order_by(_ap.c.screen_id, _ap.c.plate_id)
        apsrc = apsrc.cte('apsrc')
        
        # Altered - 20160408 
        # per discussion with JS, no need to try to link presumed assay_plates
        # from the screen_result to (screening visit) assay_plates;
        # Instead, we can infer replicate loading from the screen_result import
        # where the data_column shows the replicate
        # create a cte for the max data loaded replicates per assay_plate
        # apdlrc = select([
        #     _ap.c.screen_id,
        #     _ap.c.plate_number,
        #     func.max(_ap.c.replicate_ordinal).label('max_per_plate') ]).\
        #         select_from(_ap).\
        #         where(_ap.c.screen_result_data_loading_id != None).\
        #         group_by(_ap.c.screen_id, _ap.c.plate_number).\
        #         order_by(_ap.c.screen_id, _ap.c.plate_number)
        # apdlrc = apdlrc.cte('apdlrc')
        
        lps = (
            select([
                _ap.c.screen_id,
                func.count(distinct(_ap.c.plate_id)).label('count')
            ])
            .select_from(_ap.join(_screen,_ap.c.screen_id==_screen.c.screen_id))
            .where(_ap.c.library_screening_id != None)
            .group_by(_ap.c.screen_id))
        if facility_id:
            lps = lps.where(_screen.c.facility_id == facility_id )
        lps = lps.cte('lps')    
        # Altered - 20160408 
        # per discussion with JS, no need to try to link presumed assay_plates
        # from the screen_result to (screening visit) assay_plates;
        # Instead, we can infer replicate loading from the screen_result import
        # where the data_column shows the replicate
        # create a cte for the max data loaded replicates per assay_plate
        # lpdl = select([
        #     _ap.c.screen_id,
        #     func.count(distinct(_ap.c.plate_number)).label('count')]).\
        #         select_from(_ap).\
        #         where(_ap.c.screen_result_data_loading_id != None).\
        #         group_by(_ap.c.screen_id).cte('lpdl')
        
        libraries_screened = (
            select([
                func.count(distinct(_library.c.library_id)).label('count'),
                _ap.c.screen_id])
            .select_from(
                _ap.join(_plate, _ap.c.plate_id == _plate.c.plate_id)
                    .join(_copy, _plate.c.copy_id == _copy.c.copy_id)
                    .join(_library, _copy.c.library_id == _library.c.library_id))
            .group_by(_ap.c.screen_id).cte('libraries_screened'))
                
        # FIXME: use cherry_pick_liquid_transfer status vocabulary 
        tplcps = (
            select([
                _cpr.c.screen_id,
                func.count(1).label('count')])
            .select_from(
                _cpr.join(
                    _lcp,
                    _cpr.c.cherry_pick_request_id 
                        == _lcp.c.cherry_pick_request_id)
                .join(
                    _cpap,
                    _lcp.c.cherry_pick_assay_plate_id 
                        == _cpap.c.cherry_pick_assay_plate_id)
                .join(
                    _cplt,
                    _cplt.c.activity_id 
                        == _cpap.c.cherry_pick_liquid_transfer_id))
            .group_by(_cpr.c.screen_id)
            .where(_cplt.c.status == 'Successful').cte('tplcps'))
        
        # Create an inner screen-screen_result query to prompt the  
        # query planner to index join not hash join
        new_screen_result = (select([
                _screen.c.screen_id,
                _screen_result.c.screen_result_id,
                _screen_result.c.date_loaded,
                _screen_result.c.experimental_well_count])
            .select_from(
                _screen.join(
                    _screen_result,
                    _screen.c.screen_id == _screen_result.c.screen_id,
                    isouter=True)))
        if facility_id:
            new_screen_result = new_screen_result.where(
                _screen.c.facility_id == facility_id)
        new_screen_result = new_screen_result.cte('screen_result')
            
        lab_head_table = ScreensaverUserResource.get_lab_head_cte().cte('lab_heads')

        try:
            custom_columns = {
                # default to admin level; screen_property_generator will update
                'user_access_level_granted': literal_column('3'),
                'overlapping_positive_screens': (
                    select([
                        func.array_to_string(
                            func.array_agg(_overlap_screen.c.facility_id),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        _overlap_screens.join(
                            _overlap_screen, 
                            _overlap_screens.c.overlap_screen_id
                                ==_overlap_screen.c.screen_id))
                    .where(_overlap_screens.c.screen_id
                        == literal_column('screen.screen_id'))
                    ),
                'reconfirmation_screens': (
                    select([
                        func.array_to_string(
                            func.array_agg(_child_screen.c.facility_id),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(_child_screen)
                    .where(_child_screen.c.parent_screen_id
                        == literal_column('screen.screen_id'))),
                'primary_screen': (
                    select([_parent_screen.c.facility_id])
                    .select_from(_parent_screen)
                    .where(_parent_screen.c.screen_id
                        == literal_column('screen.parent_screen_id'))),
                'collaborator_ids': (
                    select([
                        func.array_to_string(
                            func.array_agg(collaborators.c.screensaver_user_id),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(collaborators)
                    .where(collaborators.c.screen_id 
                        == literal_column('screen.screen_id'))),
                'collaborator_names': (
                    select([
                        func.array_to_string(
                            func.array_agg(collaborators.c.name),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(collaborators)
                    .where(collaborators.c.screen_id 
                        == literal_column('screen.screen_id'))),
                'lab_affiliation': lab_head_table.c.lab_affiliation,
                'lab_name': lab_head_table.c.lab_name_full,
                'lab_head_name': lab_head_table.c.name,
                'lab_head_id': lab_head_table.c.screensaver_user_id,
                'lead_screener_name': (
                    select([_concat(_su.c.first_name, ' ', _su.c.last_name)])
                    .select_from(_su)
                    .where(_su.c.screensaver_user_id == _screen.c.lead_screener_id)),
                'lead_screener_id': (
                    select([_su.c.screensaver_user_id])
                    .select_from(_su)
                    .where(_su.c.screensaver_user_id == _screen.c.lead_screener_id)),
                'has_screen_result': (
                    case([(new_screen_result.c.screen_result_id != None, 1)],
                        else_= 3 )                    
                    ),
                'cell_lines': (
                    select([
                        func.array_to_string(
                            func.array_agg(_screen_cell_lines.c.cell_line),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(_screen_cell_lines)
                    .where(_screen_cell_lines.c.screen_id 
                        == literal_column('screen.screen_id'))),
                'date_of_first_activity': literal_column(
                    '( select date_of_activity '
                    '  from activity '
                    '  join lab_activity la using(activity_id) '
                    '  where la.screen_id=screen.screen_id '
                    '  order by date_of_activity asc LIMIT 1 )'
                    ),
                'date_of_last_activity': literal_column(
                    '( (select date_of_activity '
                    '  from activity '
                    '  join lab_activity la using(activity_id) '
                    '  where la.screen_id=screen.screen_id '
                    '  UNION ALL'
                    '  select date_of_activity '
                    '  from activity '
                    '  join service_activity sa using(activity_id) '
                    '  where sa.serviced_screen_id=screen.screen_id )'
                    '  order by date_of_activity desc LIMIT 1 )'
                    ),
                # FIXME: 20170713 - rework to be performant - 
                'activity_count': (
                    select([func.count(_activity.c.activity_id)])
                    .select_from(
                        _activity
                            .join(_service_activity, _activity.c.activity_id
                                ==_service_activity.c.activity_id,
                                isouter=True)
                            .join(_lab_activity, _activity.c.activity_id
                                == _lab_activity.c.activity_id,
                                isouter=True))
                    .where(or_(
                        _service_activity.c.serviced_screen_id
                            == literal_column('screen.screen_id'),
                        _lab_activity.c.screen_id
                            == literal_column('screen.screen_id')))
                    ),
                # # TODO: rework the update activity
                # 'screenresult_last_imported': (
                #     select([screen_result_update_activity.c.date_of_activity])
                #     .select_from(screen_result_update_activity)
                #     .where(screen_result_update_activity.c.screen_id 
                #         == literal_column('screen.screen_id'))),
                'funding_supports': (
                    select([
                        func.array_to_string(
                            func.array_agg(literal_column('funding_support')),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([_sfs.c.funding_support])
                        .select_from(_sfs)
                        .order_by(_sfs.c.funding_support)
                        .where(_sfs.c.screen_id == literal_column('screen.screen_id'))
                    .alias('inner'))),
                'total_plated_lab_cherry_picks': (
                    select([tplcps.c.count])
                    .select_from(tplcps)
                    .where(tplcps.c.screen_id == _screen.c.screen_id)),
                
                # TODO: convert to vocabulary
                'assay_readout_types': literal_column(
                    "(select array_to_string(array_agg(f1.assay_readout_type),'%s') "
                    '    from ( select distinct(assay_readout_type) '
                    '        from data_column ds '
                    '        join screen_result using(screen_result_id) '
                    '        where screen_id=screen.screen_id'
                    "        and (assay_readout_type = '') is false  ) as f1 ) " 
                    % LIST_DELIMITER_SQL_ARRAY),
                'library_plates_screened': (
                    select([lps.c.count])
                    .select_from(lps).where(lps.c.screen_id == _screen.c.screen_id)),
                'library_screenings': (
                    select([func.count(_library_screening.c.activity_id)])
                    .select_from(
                        _lab_activity.join(
                            _library_screening,
                            _library_screening.c.activity_id 
                                == _lab_activity.c.activity_id))
                    .where(_lab_activity.c.screen_id == _screen.c.screen_id)),
                'cherry_pick_screenings': (
                    select([func.count(_cp_screening.c.activity_id)])
                    .select_from(
                        _lab_activity.join(
                            _cp_screening,
                            _cp_screening.c.activity_id 
                                == _lab_activity.c.activity_id))
                    .where(_lab_activity.c.screen_id == _screen.c.screen_id)),

                # Altered - 20160408 
                # per discussion with JS, no need to try to link presumed assay_plates
                # from the screen_result to (screening visit) assay_plates;
                # Instead, we can infer replicate loading from the screen_result import
                # where the data_column shows the replicate
                # 'library_plates_data_loaded': (
                #     select([lpdl.c.count])
                #     .select_from(lpdl)
                #     .where(lpdl.c.screen_id == _screen.c.screen_id)),
                'library_plates_data_loaded': (
                    select([func.count(distinct(_aw.c.plate_number))])
                    .select_from(
                        _aw.join(_screen_result,
                            _aw.c.screen_result_id
                                == _screen_result.c.screen_result_id))
                    .where(_screen_result.c.screen_id == _screen.c.screen_id)
                    ),
                'assay_plates_screened': (
                    select([aps.c.count]).select_from(aps)
                    .where(aps.c.screen_id == _screen.c.screen_id)),
                # Altered - 20160408 
                # per discussion with JS, no need to try to link presumed assay_plates
                # from the screen_result to (screening visit) assay_plates;
                # Instead, we can infer replicate loading from the screen_result import
                # where the data_column shows the replicate
                # 'assay_plates_data_loaded': ( 
                #     select([apdl.c.count])
                #     .select_from(apdl)
                #     .where(apdl.c.screen_id == _screen.c.screen_id)),
                'libraries_screened_count': (
                    select([libraries_screened.c.count])
                    .select_from(libraries_screened)
                    .where(libraries_screened.c.screen_id == _screen.c.screen_id)),
                
                # FIXME: update extant records
                # # FIXME: use administrative activity vocabulary
                # 'last_data_loading_date': literal_column(
                #     '( select activity.date_created '
                #     '  from activity '
                #     '  join administrative_activity aa using(activity_id) '
                #     '  join screen_update_activity on update_activity_id=activity_id  '
                #     "  where administrative_activity_type = 'Screen Result Data Loading' " 
                #     '  and screen_id=screen.screen_id '
                #     '  order by date_created desc limit 1 )'
                #     ),
                'min_screened_replicate_count': (
                    select([func.min(apsrc.c.max_per_plate) + 1])
                    .select_from(apsrc)
                    .where(apsrc.c.screen_id == _screen.c.screen_id)),
                'max_screened_replicate_count': (
                    select([func.max(apsrc.c.max_per_plate) + 1])
                    .select_from(apsrc)
                    .where(apsrc.c.screen_id == _screen.c.screen_id)),
                # Altered - 20160408 
                # per discussion with JS, no need to try to link presumed assay_plates
                # from the screen_result to (screening visit) assay_plates;
                # Instead, we can infer replicate loading from the screen_result import
                # where the data_column shows the replicate
                # 'min_data_loaded_replicate_count': (
                #     select([func.min(apdlrc.c.max_per_plate) + 1])
                #     .select_from(apdlrc)
                #     .where(apdlrc.c.screen_id == _screen.c.screen_id)),
                # 'max_data_loaded_replicate_count': (
                #     select([func.max(apdlrc.c.max_per_plate) + 1])
                #     .select_from(apdlrc)
                #     .where(apdlrc.c.screen_id == _screen.c.screen_id)),
                'last_data_loading_date':
                    new_screen_result.c.date_loaded,
                'experimental_well_count': 
                    literal_column('screen_result.experimental_well_count'),
                'pin_transfer_approved_by_name': 
                    pin_transfer_approval.c.name,
                'pin_transfer_approved_by_username': 
                    pin_transfer_approval.c.username,
                'pin_transfer_date_approved': 
                    pin_transfer_approval.c.date_of_activity,
                'pin_transfer_comments': 
                    pin_transfer_approval.c.comments,
                'keywords': (
                    select([
                        func.array_to_string(
                            func.array_agg(_screen_keyword.c.keyword),
                            LIST_DELIMITER_SQL_ARRAY)])
                   .select_from(_screen_keyword)
                   .where(_screen_keyword.c.screen_id == _screen.c.screen_id)),
                'screen_id': _screen.c.screen_id,
                'publications': (
                    select([
                        func.array_to_string(
                            func.array_agg(_publication.c.title),
                            LIST_DELIMITER_SQL_ARRAY)])
                   .select_from(_publication)
                   .where(_publication.c.screen_id == _screen.c.screen_id)),
                'attached_files': (
                    select([
                        func.array_to_string(
                            func.array_agg(_attached_file.c.filename),
                            LIST_DELIMITER_SQL_ARRAY)])
                   .select_from(_attached_file)
                   .where(_attached_file.c.screen_id == _screen.c.screen_id)),
            }
        except Exception, e:
            logger.exception('on custom columns creation')
            raise 
        
        if screens_for_userid:
            custom_columns['screensaver_user_role'] = \
                screener_role_cte.c.screensaver_user_role
            
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)

        # build the query statement

        j = _screen
        if screens_for_userid:
            j = j.join(
                screener_role_cte,
                screener_role_cte.c.screen_id == _screen.c.screen_id)
        
        j = j.join(new_screen_result,
            _screen.c.screen_id == new_screen_result.c.screen_id)
        
        j = j.join(pin_transfer_approval,
            _screen.c.screen_id == pin_transfer_approval.c.screen_id, isouter=True)
        
        j = j.join(
            lab_head_table, 
            lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id,
            isouter=True)
        
        stmt = select(columns.values()).select_from(j)

        if screens_for_userid:
            stmt = stmt.where(
                screener_role_cte.c.screensaver_user_id == screens_for_userid)

        # TODO: test if more efficient filtering in sql
        # stmt = self._meta.authorization.filter_in_sql(
        #     request.user, stmt, _screen)
        
        # general setup
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        if DEBUG_SCREEN: 
            logger.info(
                'stmt: %s',
                str(stmt.compile(
                    dialect=postgresql.dialect(),
                    compile_kwargs={"literal_binds": True})))
          
        return (field_hash, columns, stmt, count_stmt, filename)

    def build_list_response(self, request, **kwargs):
        ''' 
        ScreenResource
        '''
        logger.info('ScreenResource - build_list_response')
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
         
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        manual_field_includes = set(param_hash.get('includes', []))
        
        # add fields for authorization filters
        manual_field_includes.add('user_access_level_granted')
        if self._meta.authorization.is_restricted_view(request.user):
            manual_field_includes.add('overlapping_positive_screens')
            manual_field_includes.add('screen_id')
            
        param_hash['includes'] = manual_field_includes
             
        (field_hash, columns, stmt, count_stmt, filename) = \
            self.get_query(request, schema, param_hash)
         
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            # use "use_vocab" as a proxy to also adjust siunits for viewing
            rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                field_hash, rowproxy_generator)
        
        rowproxy_generator = \
            self._meta.authorization.get_row_property_generator(
                request.user, field_hash, rowproxy_generator)
                
        title_function = None
        if use_titles is True:
            def title_function(key):
                return field_hash[key]['title']
        if is_data_interchange:
            title_function = DbApiResource.datainterchange_title_function(
                field_hash,schema['id_attribute'])
         
        return self.stream_response_from_statement(
            request, stmt, count_stmt, filename,
            field_hash=field_hash,
            param_hash=param_hash,
            is_for_detail=is_for_detail,
            rowproxy_generator=rowproxy_generator,
            title_function=title_function, meta=kwargs.get('meta', None),
            use_caching=True)
          

    @classmethod
    def get_screener_role_cte(cls):
        
        bridge = get_tables()
        _su = bridge['screensaver_user']
        _screen = bridge['screen']
        _screen_collab = bridge['screen_collaborators']
        collab = _su.alias('collab')
        ls = _su.alias('ls')
        pi = _su.alias('pi')
        
        j = _screen
        j = j.join(
            _screen_collab,
            _screen_collab.c.screen_id == _screen.c.screen_id, isouter=True)
        j = j.join(
            collab,
            collab.c.screensaver_user_id == _screen_collab.c.screensaveruser_id,
            isouter=True)
        j = j.join(
            ls, ls.c.screensaver_user_id == _screen.c.lead_screener_id,
            isouter=True)
        j = j.join(
            pi, pi.c.screensaver_user_id == _screen.c.lab_head_id,
            isouter=True)
        
        sa = (
            select([
                _screen.c.facility_id,
                _screen.c.screen_id,
                ls.c.screensaver_user_id.label('lead_screener_id'),
                pi.c.screensaver_user_id.label('pi_id'),
                func.array_agg(collab.c.screensaver_user_id).label('collab_ids')])
            .select_from(j)
            .group_by(_screen.c.facility_id, _screen.c.screen_id,
                ls.c.screensaver_user_id, pi.c.screensaver_user_id)).cte('screen_associates')
        screener_roles = (
            select([
                _su.c.screensaver_user_id,
                sa.c.facility_id,
                sa.c.screen_id,
                case([
                    (_su.c.screensaver_user_id == sa.c.lead_screener_id,
                        'lead_screener'),
                    (_su.c.screensaver_user_id == sa.c.pi_id,
                        'principal_investigator')
                    ],
                    else_='collaborator').label('screensaver_user_role')
                ])
            .select_from(sa)
            .where(or_(
                # TODO: replace with "any_()" from sqlalchemy 1.1 when avail
                _su.c.screensaver_user_id == text(' any(collab_ids) '),
                _su.c.screensaver_user_id == sa.c.lead_screener_id,
                _su.c.screensaver_user_id == sa.c.pi_id))
        )
        return screener_roles
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        Screen.objects.get(**id_kwargs).delete()
    
#     def validate(self, _dict, patch=False, current_object=None):
#         errors = DbApiResource.validate(self, _dict, patch=patch)
#         # if not errors:
#         #     errors = {}
#         #     dped = _dict.get('data_privacy_expiration_date', None)
#         #     dped_notified = _dict.get('data_privacy_expiration_notified_date', None)
#         #     min_dped = _dict.get('min_allowed_data_privacy_expiration_date', None)
#         #     max_dped = _dict.get('max_allowed_data_privacy_expiration_date', None)
#         #     if not dped:
#         #         if min_dped or max_dped:
#         #             errs['data_privacy_expiration_date'] = \
#         #                 'can not be null if min/max dates are set'
#         #     if dped_notified:
#         #         if min_dped:
#         #             errs['min_allowed_data_privacy_expiration_date'] = \
#         #                 'can not be set if the expiration notified date is set'
#         #         if self.max_allowed_data_privacy_expiration_date:
#         #             errs['min_allowed_data_privacy_expiration_date'] = \
#         #                 'can not be set if the expiration notified date is set'
#         return errors

    @write_authorization
    @un_cache
    @transaction.atomic        
    def post_detail(self, request, **kwargs):
        '''
        POST is used to create or update a resource; not idempotent;
        - The LIMS client will use POST to create exclusively
        '''
        logger.info('post_detail, screen')
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        logger.debug('patch detail %s, %s', deserialized,kwargs)
        
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.debug('param_hash: %r', param_hash)
        override_param = parse_val(
            param_hash.get(API_PARAM_OVERRIDE, False),
                API_PARAM_OVERRIDE, 'boolean')
        patch_facility_id = deserialized.get('facility_id', None)
        if patch_facility_id and override_param is not True:
            raise ValidationError({
                'facility_id':
                    'May not be specified in screen creation: %s' % patch_facility_id,
                API_PARAM_OVERRIDE: 'required' })
        
        if not patch_facility_id:    
            # find a new facility id
            max_facility_id_sql = '''
                select facility_id::text from screen 
                where parent_screen_id is null
                and study_type is null 
                and facility_id ~ '^[[:digit:]]+$' 
                order by facility_id::integer desc
                limit 1;
            '''
            # NOTE: not using facility_id ~ E'^\\d+$' syntax, which is valid when
            # "standard_conforming_strings = "on";
            # the test harness sets the standard_conforming_strings = off for the 
            # test connections; and does not recognize the E'\\d syntax at all. 
            # (even if standard_conforming_strings is set to "on").
            with get_engine().connect() as conn:
                max_facility_id = int(
                    conn.execute(max_facility_id_sql).scalar() or 0)
                patch_facility_id = str(max_facility_id+1)

        logger.info('new screen facility id to be created: %r', patch_facility_id)
        kwargs['facility_id'] = patch_facility_id
        
        return self.patch_detail(request,**kwargs)
        
    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_kwargs = self.get_id(deserialized, validate=True, schema=schema, **kwargs)
        logger.info('patch screen: %r', id_kwargs)
        create = False
        try:
            screen = Screen.objects.get(**id_kwargs)
        except ObjectDoesNotExist, e:
            create = True
            screen = Screen(**id_kwargs)

        initializer_dict = self.parse(deserialized,  schema=schema, create=create)
        logger.info('initializer_dict: %r', initializer_dict)
            
        errors = self.validate(deserialized,  schema=schema, patch=not create)
        if errors:
            raise ValidationError(errors)
        
        _key = 'lab_head_id'
        if _key in initializer_dict:
            try: 
                # may not be null
                lab_head = ScreensaverUser.objects.get(
                    screensaver_user_id=initializer_dict[_key])
                initializer_dict['lab_head'] = lab_head
                
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='No such user: %r' % initializer_dict[_key])
        _key = 'lead_screener_id'
        if _key in initializer_dict:
            try:
                # may not be null
                initializer_dict['lead_screener'] = (
                    ScreensaverUser.objects.get(
                        screensaver_user_id=initializer_dict[_key]))
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='No such user: %r' % initializer_dict[_key])
        _key = 'collaborator_ids'
        if initializer_dict.get(_key, None) is not None:
            collaborator_ids = initializer_dict[_key]
            collaborators = []
            # may empty
            for collaborator_id in collaborator_ids:
                try:
                    collaborators.append(ScreensaverUser.objects.get(
                        screensaver_user_id=collaborator_id))
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key=_key,
                        msg='No such user: %r' % collaborator_username)
            initializer_dict['collaborators'] = collaborators
        _key = 'pin_transfer_approved_by_username'
        pin_transfer_approved_by = None
        if _key in initializer_dict:
            val = initializer_dict[_key]
            if val:
                pin_transfer_approved_by = \
                    self.get_su_resource()._get_detail_response_internal(
                    exact_fields=['screensaver_user_id','is_staff'],
                    username=val)
                if not pin_transfer_approved_by: 
                    raise ValidationError(
                        key=_key,
                        msg='No such username: %r' % val)
                if pin_transfer_approved_by.get('is_staff',False) != True:
                    raise ValidationError(
                        key='pin_transfer_approved_by_username',
                        msg='Must be a staff user')
        pin_transfer_date_approved = \
            initializer_dict.get('pin_transfer_date_approved',None)
        pin_transfer_comments = \
            initializer_dict.get('pin_transfer_comments', None)
        
        related_initializer = {}
        related_initializer['cell_lines'] = \
            initializer_dict.pop('cell_lines', None)
        related_initializer['funding_supports'] = \
            initializer_dict.pop('funding_supports', None)
        related_initializer['keywords'] = \
            initializer_dict.pop('keywords', None)
        related_initializer['publications'] = \
            initializer_dict.pop('publications', None)
            
        for key, val in initializer_dict.items():
            if hasattr(screen, key):
                setattr(screen, key, val)
        screen.clean()
        screen.save()
        logger.info('screen.study_type: %r', screen.study_type)
        # NOTE: collaborators cannot be set until after the object is saved:
        # the many-to-many related manager is not functional until then.
        if 'collaborators' in initializer_dict:
            screen.collaborators = initializer_dict.get('collaborators', None)
        
        logger.info('save/created screen: %r', screen)
        
        # related objects
        
        _key = 'cell_lines'
        _val = related_initializer.get(_key, None)
        if _val is not None:
            (ScreenCellLines.objects
                .filter(screen=screen)
                .exclude(cell_line__in=_val)
                .delete())
            current_cell_lines = (
                ScreenCellLines.objects.filter(screen=screen)
                    .values_list('cell_line', flat=True))
            for cell_line in _val:
                if cell_line not in current_cell_lines:
                    ScreenCellLines.objects.create(
                        screen=screen,
                        cell_line=cell_line)
        _key = 'funding_supports'
        _val = related_initializer.get(_key, None)
        if _val is not None:
            (ScreenFundingSupports.objects
                .filter(screen=screen)
                .exclude(funding_support__in=_val)
                .delete())
            current_funding_supports = (
                ScreenFundingSupports.objects.filter(screen=screen)
                    .values_list('funding_support', flat=True))
            for funding_support in _val:
                if funding_support not in current_funding_supports:
                    ScreenFundingSupports.objects.create(
                        screen=screen,
                        funding_support=funding_support)
        
        # Set the pin transfer approval data
        if pin_transfer_approved_by is not None:
            if screen.pin_transfer_admin_activity is None:
                activity = \
                    Activity(performed_by_id=
                        pin_transfer_approved_by['screensaver_user_id'])
                activity.date_of_activity = \
                    activity.date_created
                activity.save()
                screen.pin_transfer_admin_activity = activity
                screen.save()
                logger.info('created pta: %r', 
                    screen.pin_transfer_admin_activity)
            else:
                screen.pin_transfer_admin_activity.performed_by_id = \
                    pin_transfer_approved_by['screensaver_user_id']
                screen.pin_transfer_admin_activity.save()
        if screen.pin_transfer_admin_activity is None:
            # secondary pin transfer validation
            if pin_transfer_date_approved:
                raise ValidationError(
                    key='pin_transfer_date_approved',
                    msg='requires pin_transfer_approved_by_username')    
            if pin_transfer_comments:
                raise ValidationError(
                    key='pin_transfer_comments',
                    msg='requires pin_transfer_approved_by_username')    
        else:
            if pin_transfer_date_approved:
                screen.pin_transfer_admin_activity.date_of_activity = \
                    pin_transfer_date_approved
            if pin_transfer_comments is not None:
                screen.pin_transfer_admin_activity.comments = \
                    pin_transfer_comments
            screen.pin_transfer_admin_activity.save()
            
        # TODO: determine if this is still used
        _key = 'keywords'
        _val = related_initializer.get(_key, None)
        if _val is not None:
            (ScreenKeyword.objects
                .filter(screen=screen)
                .exclude(keyword__in=_val)
                .delete())
            current_keywords = (
                ScreenKeyword.objects.filter(screen=screen)
                    .values_list('keyword', flat=True))
            for keyword in _val:
                if keyword not in current_keywords:
                    ScreenKeyword.objects.create(
                        screen=screen,
                        keyword=keyword)
        
        _key = 'publications'
        _val = related_initializer.get(_key, None)
        if _val is not None:
            current_publications = set([ 
                str(x) for x in
                    screen.publication_set.all()
                        .values_list('publication_id', flat=True)])
            logger.info('current: %r, val: %r', current_publications, _val)
            publications_delete = current_publications - set(_val)
            logger.info('delete publications: %r', publications_delete)
            if publications_delete:
                query = Publication.objects.all().filter(
                    publication_id__in=publications_delete)
                logger.info('delete: %r', [x for x in query])
                query.delete()
            # Note: prefer create publications with screen reference set
            publications_add = set(_val) - current_publications
            logger.info('add publications: %r', publications_add)
            if publications_add:
                query = Publication.objects.all().filter(
                    publication_id__in=publications_add)
                screen.publication_set.add(query)
        
        _key = 'primary_screen'
        _val = deserialized.get(_key, None)
        if _val:
            if ( screen.parent_screen 
                and screen.parent_screen.facility_id != _val):
                raise ValidationError(
                    key=_key,
                    msg='May not be reassigned')
            elif screen.parent_screen is None:
                try:
                    parent_screen = Screen.objects.get(facility_id=_val)
                    
                    if screen.screen_type != parent_screen.screen_type:
                        raise ValidationError(
                            key=_key,
                            msg='Primary screen must have the same screen_type')
                    if screen.lab_head != parent_screen.lab_head:
                        raise ValidationError(
                            key=_key,
                            msg='Primary screen must have the same lab_head')
                    
                    screen.parent_screen = parent_screen
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key=_key,
                        msg='Does not exist: %r' % _val)

        screen.save()
        logger.info('patch_obj done')
        return { API_RESULT_OBJ: screen }

class StudyResource(ScreenResource):
    
    class Meta:
        resource_name = 'study'
        max_limit = 10000
        always_return_data = True
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        authorization = ScreenAuthorization(resource_name)
        serializer = LimsSerializer()
        queryset = Screen.objects.all()  # .order_by('facility_id')
        
    def __init__(self, **kwargs):
        super(StudyResource, self).__init__(**kwargs)
    
    def prepend_urls(self):

        urls = super(StudyResource, self).prepend_urls()
        
        urls = [
            url(r"^(?P<resource_name>%s)"
                r"/create_confirmed_positive_study%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_create_confirmed_positive_study'), 
                name="api_dispatch_create_confirmed_positive_study"),
            url(r"^(?P<resource_name>%s)"
                r"/create_screened_count_study%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_create_screened_count_study'), 
                name="api_dispatch_create_screened_count_study"),
        ] + urls
        return urls
            
    def build_schema(self, user=None):
        # Bypass Screen schema
        schema = DbApiResource.build_schema(self, user=user)
        return schema

    
    @read_authorization
    def build_list_response(self, request, **kwargs):
        
        kwargs['study_type__is_null'] = False
        
        return super(StudyResource,self).build_list_response(request, **kwargs)

    @write_authorization
    @un_cache
    @transaction.atomic        
    def dispatch_create_screened_count_study(self, request, **kwargs):
        '''
        Create/Update the Screened Count Study:
        - Count of positives for each well
        - Screening count for each well
        '''

        logger.info('create the create_screened_count_study...')

        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')

        convert_post_to_put(request)
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs.pop('data')
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))
        if 'facility_id' not in deserialized:
            raise ValidationError(key='facility_id', msg='required')
        
        COL_NAME_POSITIVES_COUNT = 'Screen Positives Count'
        COL_DESC_POSITIVES_COUNT = \
            'Number of times scored as positive across all %s Screens'
        COL_NAME_SCREENED_COUNT = 'Screened Count'
        COL_DESC_SCREENED_COUNT = 'Number of times screened for all %s Screens'
        
        facility_id = deserialized['facility_id']
        data_sharing_level = deserialized.get('data_sharing_level', None)
        if data_sharing_level is None:
            deserialized['data_sharing_level'] = 0
        patch = True
        try: 
            current_study = Screen.objects.get(facility_id=facility_id)
        except ObjectDoesNotExist:
            patch = False    
        if patch is True:
            if 'study_type' in deserialized:
                del deserialized['study_type']
            if 'screen_type' in deserialized:
                del deserialized['screen_type']    
            
        _data = self.build_patch_detail(request, deserialized, **kwargs)
        _data = _data[API_RESULT_DATA][0]
        logger.info('study patched: %r', _data)
        study_obj = Screen.objects.get(facility_id=_data['facility_id'])
        logger.info('study created/updated: %r', study_obj)
        study_obj.date_created = _now()
        study_obj.study_type = 'in_silico'
        study_obj.save()
        
        col_desc_screened_count = COL_DESC_SCREENED_COUNT
        col_desc_positives_count = COL_DESC_POSITIVES_COUNT
        if study_obj.screen_type == 'small_molecule':
            col_desc_screened_count = col_desc_screened_count % 'Small Molecule'
            col_desc_positives_count = col_desc_positives_count % 'Small Molecule'
        elif study_obj.screen_type == 'rnai':
            col_desc_screened_count = col_desc_screened_count % 'RNAi'
            col_desc_positives_count = col_desc_positives_count % 'RNAi'
        else:
            raise ValidationError(key='screen_type', msg='Unknown type')
        
        # Create the Data Columns
        result_columns = OrderedDict((
            ('E', {
                'name': COL_NAME_POSITIVES_COUNT,
                'data_worksheet_column': 'E',
                'data_type': 'integer',
                'description': col_desc_positives_count,
            }),
            ('F', {
                'name': COL_NAME_SCREENED_COUNT,
                'data_worksheet_column': 'F',
                'data_type': 'integer',
                'description': col_desc_screened_count,
            }),
        ))
        
        sql = '''
        with aws as (
            select 
            well_id,
            is_positive,
            screen_id
            from assay_well
            join well w using(well_id)
            join screen_result using(screen_result_id)
            join screen using(screen_id)
            where screen.screen_type = '{screen_type}'
            and screen.study_type is null
            and w.library_well_type = 'experimental'
        )
        select
            well_id,
            count(*) as screened_count,
            count(case when is_positive then 1 end ) as positives_count
        from aws
        group by well_id
        order by well_id
        '''
        sql = sql.format(screen_type=study_obj.screen_type)
        columns = ['well_id', 'screened_count', 'positives_count']
        
        # Create the study by iterating through the report:
        # Collate the number of confirmed positives for each Screen
        with get_engine().connect() as conn:
            
            logger.info('execute the screened_count_study...')
            logger.debug('execute sql: %s', sql)
            result = conn.execute(text(sql))
            logger.info('executed the confirmed_positive_study')
            
            def result_value_generator(result):
                for row in result:
                    _dict = dict(zip(columns,row))
                    yield {
                        'well_id': _dict['well_id'],
                        'assay_well_control_type': None,
                        'E': _dict['positives_count'],
                        'F': _dict['screened_count'],
                    }
            logger.info('creating study values for: %r', study_obj)
            self.get_screenresult_resource().create_screen_result(
                request, study_obj, result_columns, 
                result_value_generator(result))
            logger.info('study values created for: %r', study_obj)
        
        if not self._meta.always_return_data:
            response_message = {'success': {
                API_MSG_RESULT: 'screen result loaded'}}
            response = self.build_response(request, response_message, **kwargs)
            return response
        else:
            new_data = self._get_detail_response_internal(facility_id=facility_id)
            data = { API_RESULT_DATA: [new_data]}
            response = self.build_response(request, data, status_code=201)
            return response 

    @write_authorization
    @un_cache
    @transaction.atomic        
    def dispatch_create_confirmed_positive_study(self, request, **kwargs):
        '''
        Create/Update the Confirmed Positives Study:
        - Count of follow-up screens for well
        - M+1 columns named "N duplexes confirming positive", 
        where 0 <= N <= M, and M is the max number of duplexes per pool in any 
        library, currently = 4). The value in each column is the number of 
        follow-up screens that confirmed the well as a positive with N duplexes
        - "Weighted Average" is the average number of confirmed positives per
        screen: 
        sum((duplexes_confirming_positive)*(count of screens))/number of screens
        '''

        logger.info('create the create_confirmed_positive_study...')

        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')

        convert_post_to_put(request)
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs.pop('data')
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))
        if 'facility_id' not in deserialized:
            raise ValidationError(key='facility_id', msg='required')
        
        VOCAB_CONFIRMED_POSITIVE = '3'
        VOCAB_FALSE_POSITIVE = '2'
        COL_NAME_WEIGHTED_AVG = 'Weighted Average'
        COL_DESC_WEIGHTED_AVG = 'Average number of confirmed positives per screen'
        COL_NAME_NUMBER_SCREENS = 'Number of screens'
        COL_DESC_NUMBER_SCREENS = 'Number of follow up screens testing duplexes for the pool well'
        COL_NAME_DUPLEX_COUNT = 'Number of screens confirming with %d duplexes'
        COL_DESC_DUPLEX_COUNT = 'Number of screens confirming with %d duplexes'

        facility_id = deserialized['facility_id']
        data_sharing_level = deserialized.get('data_sharing_level', None)
        if data_sharing_level is None:
            deserialized['data_sharing_level'] = 0
        patch = True
        try: 
            current_study = Screen.objects.get(facility_id=facility_id)
        except ObjectDoesNotExist:
            patch = False    
        if patch is True:
            if 'study_type' in deserialized:
                del deserialized['study_type']
            if 'screen_type' in deserialized:
                del deserialized['screen_type']    
            
        _data = self.build_patch_detail(request, deserialized, **kwargs)
        _data = _data[API_RESULT_DATA][0]
        logger.info('study patched: %r', _data)
        study_obj = Screen.objects.get(facility_id=_data['facility_id'])
        logger.info('study created/updated: %r', study_obj)
        study_obj.study_type = 'in_silico'
        study_obj.date_created = _now()
        study_obj.save()
        
        # Create the Data Columns
        result_columns = OrderedDict((
            ('E', {
                'name': COL_NAME_WEIGHTED_AVG,
                'data_worksheet_column': 'E',
                'data_type': 'numeric',
                'decimal_places': 2, 
                'description': COL_DESC_WEIGHTED_AVG,
            }),
            ('F', {
                'name': COL_NAME_NUMBER_SCREENS,
                'data_worksheet_column': 'F',
                'data_type': 'integer',
                'description': COL_DESC_NUMBER_SCREENS,
            }),
        ))

        for i in range(0,5):
            data_worksheet_column_letter = chr(ord('F')+i+1)
            result_columns[data_worksheet_column_letter] = {
                'name': COL_NAME_DUPLEX_COUNT % i,
                'description': COL_DESC_DUPLEX_COUNT %i,
                'data_type': 'integer',
                'data_worksheet_column': data_worksheet_column_letter
            }
        
        # SQL report: all pool reagents -> duplex reagents -> screen -> confirmed_positive_value
        sql = '''
            with pool_reagents as (
            select r.well_id, pr.* from silencing_reagent pr
              join reagent r using(reagent_id)
              join well w using(well_id)
              join library l using(library_id)
              where l.is_pool is true
              and l.screen_type = 'rnai' )
            select 
                pr.well_id as pr_well_id, 
                pr.reagent_id as pr_id,
                dr.well_id as dr_well_id, 
                dr.reagent_id as dr_id, 
                aw.confirmed_positive_value, 
                sr.screen_id,
                s.facility_id
              from pool_reagents pr
              join reagent prr using(reagent_id)
              join well prw on(prw.well_id=prr.well_id) 
              join silencing_reagent_duplex_wells srdw on(pr.reagent_id=srdw.silencingreagent_id)
              join well dw on(dw.well_id=srdw.well_id)
              join reagent dr on(dr.well_id=dw.well_id)
              join assay_well aw on(dw.well_id=aw.well_id)
              join screen_result sr using(screen_result_id)
              join screen s using(screen_id)
              where aw.confirmed_positive_value 
                  in ('{vocab_confirmed_positive}','{vocab_false_positive}')
              order by pr_id, dr_id, sr.screen_id
        '''
            
        sql = sql.format(vocab_confirmed_positive=VOCAB_CONFIRMED_POSITIVE,
            vocab_false_positive=VOCAB_FALSE_POSITIVE)
        columns = ['pr_well_id', 'pr_id', 'dr_well_id', 'dr_id', 
            'confirmed_positive_value', 'screen_id','facility_id']
        
        # Create the study by iterating through the report:
        # Collate the number of confirmed positives for each Screen
        with get_engine().connect() as conn:
            
            logger.info('execute the confirmed_positive_study...')
            logger.debug('execute sql: %s', sql)
            result = conn.execute(text(sql))
            logger.info('executed the confirmed_positive_study')
            
            def result_value_generator(result):
                
                def create_result(well_id, screens, screen_confirmed_positives):
                    screen_count = len(screens)
                    screens_confirming_zero_duplexes = len(screens)-len(screen_confirmed_positives)
                    result_row = {
                        'well_id': well_id  ,
                        'assay_well_control_type': None,
                        'F': screen_count ,
                        'G': screens_confirming_zero_duplexes
                    }
                    
                    weighted_value_sum = 0
                    weighted_average = 0
                    pool_reagent_counts = defaultdict(int)
                    
                    for i in range(1,5):
                        data_worksheet_column_letter = chr(ord('G')+i)
                        pool_reagent_counts[i] = 0
                        screens_for_weight = len([
                            num for num in screen_confirmed_positives.values() 
                            if num==i ])
                        weighted_value_sum += i * screens_for_weight
                        result_row[data_worksheet_column_letter] = screens_for_weight
                    if weighted_value_sum > 0:
                        weighted_average = round(weighted_value_sum/(1.0*screen_count), 2)
                        result_row['E'] = weighted_average
                    logger.debug('yield row: %r', result_row)
                    return result_row
                
                screens = None
                screen_confirmed_positives = None
                current_pool_well = None
                for row in result:
                    _dict = dict(zip(columns,row))
                    logger.debug('dict: %r', _dict)
                    if current_pool_well != _dict['pr_well_id']:
                        if current_pool_well is not None:
                            yield create_result(
                                current_pool_well, screens, 
                                screen_confirmed_positives)
                        current_pool_well = _dict['pr_well_id']
                        screens = set()
                        screen_confirmed_positives = defaultdict(int)

                    screens.add(_dict['facility_id'])
                    if _dict['confirmed_positive_value'] == VOCAB_CONFIRMED_POSITIVE:
                        screen_confirmed_positives[_dict['facility_id']] += 1
                if screens:
                    yield create_result(
                        current_pool_well, screens, screen_confirmed_positives)

            logger.info('creating study values for: %r', study_obj)
            self.get_screenresult_resource().create_screen_result(
                request, study_obj, result_columns, 
                result_value_generator(result))
            logger.info('study values created for: %r', study_obj)
        
        if not self._meta.always_return_data:
            response_message = {'success': {
                API_MSG_RESULT: 'screen result loaded'}}
            response = self.build_response(request, response_message, **kwargs)
            return response
        else:
            new_data = self._get_detail_response_internal(facility_id=facility_id)
            data = { API_RESULT_DATA: [new_data]}
            response = self.build_response(request, data, status_code=201)
            return response 

    @write_authorization
    @un_cache
    @transaction.atomic        
    def post_detail(self, request, **kwargs):
        '''
        POST is used to create or update a resource; not idempotent;
        - The LIMS client will use POST to create exclusively
        '''
        logger.info('post_detail, study')
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        logger.debug('patch detail %s, %s', deserialized,kwargs)

        patch_facility_id = deserialized.get('facility_id', None)
        if not patch_facility_id:    
            raise ValidationError(
                key='facility_id', msg='required')
        logger.info('new study facility id to be created: %r', patch_facility_id)
        kwargs['facility_id'] = patch_facility_id
        
        return self.patch_detail(request,**kwargs)

    @write_authorization
    @un_cache  
    @transaction.atomic      
    def patch_detail(self, request, **kwargs):
        '''
        PATCH is used to create or update a resource; not idempotent
        '''
        
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs.pop('data')
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))
        data_sharing_level = deserialized.get('data_sharing_level', None)
        if data_sharing_level is None:
            deserialized['data_sharing_level'] = 0
        _data = self.build_patch_detail(request, deserialized, **kwargs)

        return self.build_response(
            request, _data, response_class=HttpResponse, **kwargs)

class ScreensaverUserResourceAuthorization(UserResourceAuthorization):

    def get_user_screens(self, screensaver_user):
        my_screens = Screen.objects.all().filter(
            Q(lead_screener=screensaver_user)
            | Q(lab_head=screensaver_user)
            | Q(collaborators=screensaver_user))
        logger.info('user: %r, screens: %r', 
            screensaver_user.username, [s.facility_id for s in my_screens])
        return [screen for screen in my_screens]
    
    def get_associated_users(self, user):
        ''' collaborators '''
        
        screensaver_user = ScreensaverUser.objects.get(username=user.username)
        associated_users = set([screensaver_user])
        
        for screen in self.get_user_screens(screensaver_user):
            associated_users.add(screen.lead_screener)
            associated_users.add(screen.lab_head)
            associated_users.update([su for su in screen.collaborators.all()])

        logger.info('user: %r, associated users: %r', user.username, associated_users)
        
        return associated_users
        
class UserChecklistAuthorization(ScreensaverUserResourceAuthorization):        

    def _is_resource_authorized(
        self, user, permission_type, **kwargs):
        authorized = \
            super(UserChecklistAuthorization, self)\
                ._is_resource_authorized(user, permission_type, **kwargs)
        if authorized is True:
            return True
        
        return user.is_active

    def filter(self, user, filter_expression):
        
        if self.is_restricted_view(user) is False:
            return filter_expression
        
        auth_filter = column('username') == user.username
        if filter_expression is not None:
            filter_expression = and_(filter_expression, auth_filter)
        else:
            filter_expression = auth_filter

        return filter_expression

    def get_row_property_generator(self, user, fields, extant_generator):
        # If the user may see the CPR, there are no property restrictions
        return extant_generator

    
class UserChecklistResource(DbApiResource):

    class Meta:
        
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'userchecklist'
        authorization = UserChecklistAuthorization(resource_name)
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
        
        self.screensaveruser_resource = None
        super(UserChecklistResource, self).__init__(**kwargs)

    def get_su_resource(self):
        if self.screensaveruser_resource is None:
            self.screensaveruser_resource = ScreensaverUserResource()
        return self.screensaveruser_resource
    
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/" 
                 r"(?P<name>([\w_]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]    

    @read_authorization
    def get_detail(self, request, **kwargs):

        name = kwargs.get('name', None)
        if not name:
            raise NotImplementedError('must provide a name parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        checklist_name = param_hash.pop('name', None)
        if checklist_name:
            param_hash['name__eq'] = checklist_name
        
        # general setup
        
        manual_field_includes = set(param_hash.get('includes', []))
        manual_field_includes.add('is_checked');
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, schema)
        filter_expression = \
            self._meta.authorization.filter(request.user,filter_expression)
              
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
         
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
        rowproxy_generator = \
            self._meta.authorization.get_row_property_generator(
                request.user, field_hash, rowproxy_generator)
        
        # specific setup
        #             _su = self.bridge['screensaver_user']
        _user_cte = ScreensaverUserResource.get_user_cte().cte('cl_user')
        _admin_cte = ScreensaverUserResource.get_user_cte().cte('cl_admin')
        #             _admin = _su.alias('admin')
        _up = self.bridge['reports_userprofile']
        _user_checklist = self.bridge['user_checklist']
        _vocab = self.bridge['reports_vocabulary']
        
        # get the checklist vocabulary
        checklist_table = (
            select([
                _vocab.c.ordinal,
                _vocab.c.key.label('name'),
                _vocab.c.is_retired.label('is_retired'),
            ])
            .select_from(_vocab)
            .where(_vocab.c.scope == 'userchecklist.name'))
        checklist_table = Alias(checklist_table)
        
        # build the entered checklists
        
        j = _user_checklist
        j = j.join(
            _user_cte, _user_checklist.c.screensaver_user_id 
                == _user_cte.c.screensaver_user_id)
        j = j.join(
            _admin_cte, _user_checklist.c.admin_user_id 
                == _admin_cte.c.screensaver_user_id)
        entered_checklists = select([
            _user_cte.c.screensaver_user_id,
            _user_cte.c.username,
            _user_cte.c.name.label('user_name'),
            _user_cte.c.first_name.label('user_first_name'),
            _user_cte.c.last_name.label('user_last_name'),
            _user_checklist.c.name,
            _user_checklist.c.is_checked,
            _user_checklist.c.date_effective,
            _user_checklist.c.date_notified,
            _admin_cte.c.username.label('admin_username'),
            _admin_cte.c.name.label('admin_name'),
            
            ]).select_from(j)
        username = param_hash.pop('username', None)
        if username:
            entered_checklists = entered_checklists.where(
                _user_cte.c.username == username)
        screensaver_user_id = param_hash.pop('screensaver_user_id', None)
        if screensaver_user_id:
            entered_checklists = entered_checklists.where(
                _user_cte.c.screensaver_user_id == screensaver_user_id)
        entered_checklists = entered_checklists.cte('entered_checklists')
        
        # This entire query doesn't fit the pattern, construct it manually
        custom_columns = {
            'screensaver_user_id': func.coalesce(
                entered_checklists.c.screensaver_user_id, screensaver_user_id),
            'username': func.coalesce(
                entered_checklists.c.username, username),
            'user_name': entered_checklists.c.user_name,
            'user_first_name': entered_checklists.c.user_first_name,
            'user_last_name': entered_checklists.c.user_last_name,
            'admin_username': entered_checklists.c.admin_username,
            'admin_name': entered_checklists.c.admin_name,
            'name' : checklist_table.c.name,
            'is_checked': func.coalesce(
                entered_checklists.c.is_checked, False),
            'status': case([
                ( entered_checklists.c.is_checked == True, 'activated')],
                else_=case([
                    (entered_checklists.c.date_effective != None, 'deactivated')],
                    else_='not_completed')) ,
            'date_effective': entered_checklists.c.date_effective,
            'date_notified': entered_checklists.c.date_notified,
            'is_retired': checklist_table.c.is_retired,
            }
        
        base_query_tables = ['user_checklist', 'screensaver_user'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        isouter = False
        if username is not None or screensaver_user_id is not None:
            # if username/suid, then this is a user specific view:
            # - outer join in the two so that a full list is generated
            isouter = True
            
        j = checklist_table
        j = j.join(
            entered_checklists,
            checklist_table.c.name == entered_checklists.c.name,
            isouter=isouter)
        
        stmt = select(columns.values()).select_from(j)
        if username is None and screensaver_user_id is None:
            stmt = stmt.order_by(
                entered_checklists.c.user_last_name, 
                entered_checklists.c.user_first_name)
        # general setup
        if 'is_retired' not in readable_filter_hash:
            stmt = stmt.where(
                func.coalesce(checklist_table.c.is_retired,False) != True)
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        
        # compiled_stmt = str(stmt.compile(
        #     dialect=postgresql.dialect(),
        #     compile_kwargs={"literal_binds": True}))
        # logger.info('compiled_stmt %s', compiled_stmt)
        
        title_function = None
        if use_titles is True:
            def title_function(key):
                return field_hash[key]['title']
        if is_data_interchange:
            title_function = DbApiResource.datainterchange_title_function(
                field_hash,schema['id_attribute'])
        
        return self.stream_response_from_statement(
            request, stmt, count_stmt, filename,
            field_hash=field_hash,
            param_hash=param_hash,
            is_for_detail=is_for_detail,
            rowproxy_generator=rowproxy_generator,
            title_function=title_function, meta=kwargs.get('meta', None),
            use_caching=True)

    @write_authorization
    @un_cache
    @transaction.atomic
    def delete_obj(self, request, deserialized, **kwargs):
        raise NotImplementedError(
            'delete obj is not implemented for UserChecklist')

    @write_authorization
    @transaction.atomic()
    def patch_obj(self, request, deserialized, **kwargs):
        
        logger.info('patch checklist item: %r', deserialized)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        
        screensaver_user_id = deserialized.get('screensaver_user_id', None)
        if not screensaver_user_id:
            screensaver_user_id = kwargs.get('screensaver_user_id', None)
        if not screensaver_user_id:
            raise ValidationError(key='screensaver_user_id', msg='required')
        item_name = deserialized.get('name', None)
        if not item_name:
            raise ValidationError(key='name', msg='required')

        try:
            user = ScreensaverUser.objects.get(screensaver_user_id=screensaver_user_id)
        except ObjectDoesNotExist:
            raise ValidationError(
                key='screensaver_user_id',
                msg='screensaver_user_id does not exist: %r' % screensaver_user_id)
        create = False
        try:
            uci = UserChecklist.objects.get(
                screensaver_user=user,
                name=item_name)
            patch = True
            logger.info('UserChecklist to patch: %r', uci)
        except ObjectDoesNotExist:
            logger.info(
                'UserChecklist does not exist: %s/%s, creating' 
                % (screensaver_user_id, item_name))
            uci = UserChecklist(
                screensaver_user=user,
                name=item_name)
            create = True

        initializer_dict = self.parse(deserialized, create=create, schema=schema)
        errors = self.validate(initializer_dict, patch=not create, schema=schema)
        if errors:
            raise ValidationError(errors)

        admin_username = deserialized.get('admin_username', None)
        if not admin_username:
            raise ValidationError(
                key='admin_username',
                msg='required')
        try:
            admin_user = ScreensaverUser.objects.get(username=admin_username)
            initializer_dict['admin_user_id'] = admin_user.pk
            # Note, hasattr does not work for foreign keys if not initialized, 
            # must use the key
        except ObjectDoesNotExist:
            raise ValidationError(
                key='admin_username',
                msg='username does not exist: %r' % admin_username)

        for key, val in initializer_dict.items():
            if hasattr(uci, key):
                setattr(uci, key, val)
        
        uci.save()
        logger.info('UserChecklist %r created: %r', uci, create)
        return { API_RESULT_OBJ: uci }


class LabAffiliationResource(DbApiResource):

    class Meta:

        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'labaffiliation'
        authorization = UserGroupAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        
        super(LabAffiliationResource, self).__init__(**kwargs)
        self.su_resource = None
        
    def get_screensaver_resource(self):
        if self.su_resource is None:
            self.su_resource = ScreensaverUserResource()
        return self.su_resource

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<lab_affiliation_id>([\d]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]

    @read_authorization
    def get_detail(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)

        try:
            
            # general setup
            
            manual_field_includes = set(param_hash.get('includes', []))
            exact_fields = set(param_hash.get('exact_fields', []))
        
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                  
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            # specific setup
            _la = self.bridge['lab_affiliation']
            _lab_su = self.bridge['screensaver_user']
            
            _lab_head_cte = self.get_screensaver_resource().get_user_cte().cte('la_lab_heads')
            
            custom_columns = {
                'lab_head_ids': (
                    select([
                        func.array_to_string(
                            func.array_agg(_lab_head_cte.c.screensaver_user_id),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(_lab_su.join(
                        _lab_head_cte, 
                        _lab_su.c.screensaver_user_id
                            ==_lab_head_cte.c.screensaver_user_id))
                    .where(_lab_su.c.lab_affiliation_id 
                        == literal_column('lab_affiliation.lab_affiliation_id'))),
                'lab_head_names': (
                    select([
                        func.array_to_string(
                            func.array_agg(_lab_head_cte.c.name),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(_lab_su.join(
                        _lab_head_cte, 
                        _lab_su.c.screensaver_user_id
                            ==_lab_head_cte.c.screensaver_user_id))
                    .where(_lab_su.c.lab_affiliation_id
                        == literal_column('lab_affiliation.lab_affiliation_id'))),
                }

            # delegate to the user resource
            base_query_tables = [
                'lab_affiliation'] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            # build the query statement
            
            j = _la
            
            stmt = select(columns.values()).select_from(j)
            # natural ordering
            stmt = stmt.order_by(_la.c.category, _la.c.name)
            
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            # logger.info(
            #     'stmt: %s',
            #     str(stmt.compile(
            #         dialect=postgresql.dialect(),
            #         compile_kwargs={"literal_binds": True})))
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
             
        except Exception, e:
            logger.exception('on get_list')
            raise e  
        
    @write_authorization
    @transaction.atomic()
    def patch_obj(self, request, deserialized, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        logger.info('patch lab affiliation: %r, %r', deserialized, kwargs)
        
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)

        patch = bool(id_kwargs)
        initializer_dict = self.parse(deserialized, schema=schema, create=not patch)
        errors = self.validate(initializer_dict, schema=schema, patch=patch)
        if errors:
            raise ValidationError(errors)

        lab_affiliation = None
        if patch is True:
            try:
                lab_affiliation = LabAffiliation.objects.get(**id_kwargs)
            except ObjectDoesNotExist:
                raise Http404(
                    'Lab Affiliation does not exist for: %r', id_kwargs)
        else:
            lab_affiliation = LabAffiliation.objects.create(**id_kwargs)

        for key, val in initializer_dict.items():
            if hasattr(lab_affiliation, key):
                setattr(lab_affiliation, key, val)
        
        lab_affiliation.save()
        logger.info('LabAffiliation %r patch: %r', lab_affiliation, patch)
        return { API_RESULT_OBJ: lab_affiliation }


class ScreensaverUserResource(DbApiResource):    
    
    VOCAB_USER_AGREEMENT_RNAI = UserAgreementResource.VOCAB_USER_AGREEMENT_RNAI
    VOCAB_USER_AGREEMENT_SM = UserAgreementResource.VOCAB_USER_AGREEMENT_SM

    class Meta:

        # queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'screensaveruser'
        authorization = ScreensaverUserResourceAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['digested_password']
        max_limit = 10000
        always_return_data = True
        # TODO: utilize the cache_control mechanism to signal cache status
        # to the client (max-age, etag)
        # cache = SimpleCache(timeout=10)

    def __init__(self, **kwargs):
        
        self.user_resource = None
        super(ScreensaverUserResource, self).__init__(**kwargs)

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/groups%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_groupview'),
                name="api_dispatch_user_groupview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/groups%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_groupview'),
                name="api_dispatch_user_groupview"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/checklist%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_checklistview'),
                name="api_dispatch_user_checklistview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/checklist%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_checklistview'),
                name="api_dispatch_user_checklistview"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/attachedfiles%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_attachedfileview'),
                name="api_dispatch_user_attachedfileview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/attachedfiles%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_attachedfileview'),
                name="api_dispatch_user_attachedfileview"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/useragreement%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_useragreement_view'),
                name="api_dispatch_useragreement_view"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/useragreement%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_useragreement_view'),
                name="api_dispatch_useragreement_view"),
            url((r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))"
                 r"/attachedfiles/(?P<attached_file_id>([\d]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_attachedfiledetailview'),
                name="api_dispatch_user_attachedfiledetailview"),
            url((r"^(?P<resource_name>%s)/(?P<username>([\w]+))"
                 r"/attachedfiles/(?P<attached_file_id>([\d]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_attachedfiledetailview'),
                name="api_dispatch_user_attachedfiledetailview"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/serviceactivities%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_serviceactivityview'),
                name="api_dispatch_user_serviceactivityview"),
            url((r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))"
                 r"/serviceactivities/(?P<activity_id>([\d]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_serviceactivitydetailview'),
                name="api_dispatch_user_serviceactivitydetailview"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/activities%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_activityview'),
                name="api_dispatch_user_activityview"),
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>([\d]+))/screens%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_screenview'),
                name="api_dispatch_user_screenview"),
        ]    

    def dispatch_user_groupview(self, request, **kwargs):
        username = kwargs.pop('username', None)
        if username:
            su = ScreensaverUser.objects.get(username=username)
            kwargs['screensaver_user_id'] = su.screensaver_user_id
            
        return UserGroupResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_checklistview(self, request, **kwargs):
        return UserChecklistResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_attachedfileview(self, request, **kwargs):
        method = 'list'
        if request.method.lower() == 'post':
            # if put is used, force to "post_detail"
            method = 'detail'
        return AttachedFileResource().dispatch(method, request, **kwargs)    

    def dispatch_useragreement_view(self, request, **kwargs):
        method = 'list'
        if request.method.lower() == 'post':
            # if put is used, force to "post_detail"
            method = 'detail'
        return UserAgreementResource().dispatch(method, request, **kwargs)    

    def dispatch_user_attachedfiledetailview(self, request, **kwargs):
        return AttachedFileResource().dispatch('detail', request, **kwargs)    
    
    def dispatch_user_serviceactivityview(self, request, **kwargs):
        kwargs['serviced_user_id__eq'] = kwargs.pop('screensaver_user_id')
        return ServiceActivityResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_serviceactivitydetailview(self, request, **kwargs):
        return ServiceActivityResource().dispatch('detail', request, **kwargs)    
    
    def dispatch_user_activityview(self, request, **kwargs):
        screensaver_user_id = kwargs.pop('screensaver_user_id')
        
        # TODO: 20171111 - performance
        # Can be made more performant by directly joining user to activity table
        nested_search_data = [
            { 'serviced_user_id__eq' : screensaver_user_id},
            {'performed_by_user_id__eq': screensaver_user_id},
        ]
        
        kwargs['nested_search_data'] = nested_search_data
        return ActivityResource().dispatch('list', request, **kwargs)    

    def dispatch_user_screenview(self, request, **kwargs):
        kwargs['screens_for_userid'] = kwargs.pop('screensaver_user_id')
        return ScreenResource().dispatch('list', request, **kwargs)    

    def build_schema(self, user=None):
        
        schema = super(ScreensaverUserResource, self).build_schema(user=user)
        sub_schema = self.get_user_resource().build_schema(user=user);
        fields = {}
        fields.update(sub_schema['fields'])
        for key, val in schema['fields'].items():
            if key in fields:
                fields[key].update(val)
            else:
                fields[key] = val
        schema['fields'] = fields
        logger.debug('=== final screensaver_user fields: %r', fields)
        return schema
    
    @classmethod
    def get_user_cte(cls):
        
        bridge = get_tables()
        _su = bridge['screensaver_user']
#         _up = bridge['reports_userprofile']
#         _au = bridge['auth_user']
        
        j = _su
        user_table = (
            select([
                _su.c.screensaver_user_id,
                _su.c.username,
                _concat(_su.c.first_name, ' ', _su.c.last_name).label('name'),
                _concat(_su.c.last_name, ', ', _su.c.first_name).label('last_first'),
                _su.c.first_name,
                _su.c.last_name,
                _su.c.email,
                _su.c.lab_head_id
                ])
            .select_from(j))
        return user_table
        
    @classmethod
    def get_lab_head_cte(cls, alias_qualifier=''):
        bridge = get_tables()
        _su = bridge['screensaver_user']
        _user = ScreensaverUserResource.get_user_cte().cte(
            'lab_head_users_%s'%alias_qualifier)
        affiliation_table = bridge['lab_affiliation']
        _vocab = bridge['reports_vocabulary']
        la_categories = (
            select([_vocab.c.key, _vocab.c.scope, _vocab.c.title ])
                .select_from(_vocab)
                .where(_vocab.c.scope=='labaffiliation.category')).cte(
                    'labaffiliation_category')
        
        
        lab_head_table = (
            select([
                _user.c.screensaver_user_id,
                _user.c.name,
                _user.c.username,
                affiliation_table.c.lab_affiliation_id,
                affiliation_table.c.name.label('lab_affiliation_name'),
                affiliation_table.c.category.label('lab_affiliation_category'),
                _concat(
                    affiliation_table.c.name, ' (',
                    la_categories.c.title, ')').label('lab_affiliation'),
                func.array_to_string(
                    array([
                        _user.c.last_first,
                        ' - ', affiliation_table.c.name,
                        ' (', la_categories.c.title, ')']), '').label('lab_name_full'),
                ])
            .select_from(
                _user.join(
                    _su, _user.c.screensaver_user_id==_su.c.screensaver_user_id)
                .join(affiliation_table, 
                    _su.c.lab_affiliation_id==affiliation_table.c.lab_affiliation_id,
                    isouter=True )
                .join(la_categories, affiliation_table.c.category==la_categories.c.key))
            .where(_su.c.classification=='principal_investigator')
            )
        return lab_head_table
        
    def get_user_resource(self):

        if not self.user_resource:
            self.user_resource = UserResource()
        return self.user_resource
        
    @read_authorization
    def get_detail(self, request, **kwargs):

        screensaver_user_id = kwargs.get('screensaver_user_id', None)
        username = kwargs.get('username', None)
        ecommons_id = kwargs.get('ecommons_id', None)
        if screensaver_user_id is None and username is None and ecommons_id is None:
            logger.info(
                'no screensaver_user_id, username or ecommons_id provided: %r', 
                kwargs.keys())
            raise NotImplementedError(
                'must provide a screensaver_user_id or username parameter')

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)

        # general setup
            
        manual_field_includes = set(param_hash.get('includes', []))
        exact_fields = set(param_hash.get('exact_fields', []))
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, schema)

        filter_expression = self._meta.authorization.filter(
            request.user, filter_expression)
              
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
             
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
        rowproxy_generator = \
            self._meta.authorization.get_row_property_generator(
                request.user, field_hash, rowproxy_generator)
 
        # specific setup
        _su = self.bridge['screensaver_user']
        _au = self.bridge['auth_user']
        _up = self.bridge['reports_userprofile']
        _s = self.bridge['screen']
        _screen_collab = self.bridge['screen_collaborators']
        _fur = self.bridge['user_facility_usage_role']
        _lhsu = _su.alias('lhsu')
        _user_agreement = self.bridge['user_agreement']

        lab_head_table = ScreensaverUserResource.get_lab_head_cte().cte('lab_heads')
        lab_member = ScreensaverUserResource.get_user_cte().cte('users_labmember')
            
        _all_screens = union(
            select([
                _screen_collab.c.screensaveruser_id.label('screensaver_user_id'),
                _screen_collab.c.screen_id,
                _s.c.facility_id, 
                literal_column("'collaborator'").label('role')])
            .select_from(_screen_collab.join(
                _s, _s.c.screen_id==_screen_collab.c.screen_id)),
            select([
                _s.c.lead_screener_id.label('screensaver_user_id'),
                _s.c.screen_id, 
                _s.c.facility_id, 
                literal_column("'lead_screener'").label('role')]),
            select([
                _s.c.lab_head_id.label('screensaver_user_id'),
                _s.c.screen_id,
                _s.c.facility_id, 
                literal_column("'lab_head'").label('role')])
        ).order_by('facility_id').cte('all_screens')
                
            
        custom_columns = {
            'is_staff': func.coalesce(_au.c.is_staff, False),
            'is_active': func.coalesce(_au.c.is_active, False),
            'is_superuser': func.coalesce(_au.c.is_superuser, False),
            'username': _au.c.username,
            'name': _concat_with_sep(
                args=[_su.c.last_name,_su.c.first_name,], sep=', '),
            'screens': (
                select([func.array_to_string(
                    func.array_agg(_all_screens.c.facility_id),
                    LIST_DELIMITER_SQL_ARRAY)])
                .select_from(_all_screens)
                .where(_all_screens.c.screensaver_user_id
                    ==_su.c.screensaver_user_id)
                ),
            'screens_lab_head': (
                select([func.array_to_string(
                    func.array_agg(_all_screens.c.facility_id),
                    LIST_DELIMITER_SQL_ARRAY)])
                .select_from(_all_screens)
                .where(_all_screens.c.screensaver_user_id
                    ==_su.c.screensaver_user_id)
                .where(_all_screens.c.role=='lab_head')
                ),
            'screens_lead': (
                select([func.array_to_string(
                    func.array_agg(_all_screens.c.facility_id),
                    LIST_DELIMITER_SQL_ARRAY)])
                .select_from(_all_screens)
                .where(_all_screens.c.screensaver_user_id
                    ==_su.c.screensaver_user_id)
                .where(_all_screens.c.role=='lead_screener')
                ),
            'screens_collaborator': (
                select([func.array_to_string(
                    func.array_agg(_all_screens.c.facility_id),
                    LIST_DELIMITER_SQL_ARRAY)])
                .select_from(_all_screens)
                .where(_all_screens.c.screensaver_user_id
                    ==_su.c.screensaver_user_id)
                .where(_all_screens.c.role=='collaborator')
                ),
            'lab_name': lab_head_table.c.lab_name_full,
            'lab_head_id': lab_head_table.c.screensaver_user_id,
            'lab_member_ids': (
                select([
                    func.array_to_string(
                        func.array_agg(literal_column('screensaver_user_id')),
                        LIST_DELIMITER_SQL_ARRAY)])
                .select_from(
                    select([lab_member.c.screensaver_user_id])
                    .select_from(lab_member)
                    .order_by(lab_member.c.last_first)
                    .where(lab_member.c.lab_head_id
                        ==literal_column('screensaver_user.screensaver_user_id'))
                    .where(lab_member.c.lab_head_id!=lab_member.c.screensaver_user_id)
                    .alias('inner'))        
                ),
            'lab_member_names': (
                select([
                    func.array_to_string(
                        func.array_agg(literal_column('name')),
                        LIST_DELIMITER_SQL_ARRAY)])
                .select_from(
                    select([lab_member.c.name, lab_member.c.email])
                    .select_from(lab_member)
                    .order_by(lab_member.c.last_first)
                    .where(lab_member.c.lab_head_id
                        ==literal_column('screensaver_user.screensaver_user_id'))
                    .where(lab_member.c.lab_head_id!=lab_member.c.screensaver_user_id)
                    .alias('inner'))        
                ),
            'lab_member_emails': (
                select([
                    func.array_to_string(
                        func.array_agg(func.coalesce(literal_column('email'),'null')),
                        LIST_DELIMITER_SQL_ARRAY)])
                .select_from(
                    select([lab_member.c.name,lab_member.c.email])
                    .select_from(lab_member)
                    .order_by(lab_member.c.last_first)
                    .where(lab_member.c.lab_head_id
                        ==literal_column('screensaver_user.screensaver_user_id'))
                    .where(lab_member.c.lab_head_id!=lab_member.c.screensaver_user_id)
                    .alias('inner'))        
                ),
            'lab_affiliation_name': lab_head_table.c.lab_affiliation_name,
            'lab_affiliation_category': lab_head_table.c.lab_affiliation_category,
            'facility_usage_roles': (
                select([
                    func.array_to_string(
                        func.array_agg(_fur.c.facility_usage_role),
                        LIST_DELIMITER_SQL_ARRAY)])
                .select_from(_fur)
                .where(_fur.c.screensaver_user_id 
                    == _su.c.screensaver_user_id)),
            'sm_data_sharing_level': (
                select([_user_agreement.c.data_sharing_level])
                .select_from(_user_agreement)
                .where(_user_agreement.c.type==self.VOCAB_USER_AGREEMENT_SM)
                .where(_user_agreement.c.date_expired==None)
                .where(_user_agreement.c.screensaver_user_id
                    ==_su.c.screensaver_user_id)),
            'rnai_data_sharing_level': (
                select([_user_agreement.c.data_sharing_level])
                .select_from(_user_agreement)
                .where(_user_agreement.c.type==self.VOCAB_USER_AGREEMENT_RNAI)
                .where(_user_agreement.c.date_expired==None)
                .where(_user_agreement.c.screensaver_user_id
                    ==_su.c.screensaver_user_id)),
        }

        # delegate to the user resource
        default_fields = ['fields.screensaveruser', 'fields.user']
        _temp = { key:field for key, field in field_hash.items() 
            if field.get('scope', None) in default_fields }
        field_hash = _temp
        logger.debug('final field hash: %s', field_hash.keys())
        sub_columns = self.get_user_resource().build_sqlalchemy_columns(
            field_hash.values(),
            custom_columns=custom_columns)
        base_query_tables = [
            'screensaver_user', 'reports_user', 'auth_user'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=sub_columns)

        # build the query statement
        
        j = _su
        j = j.join(_up, _su.c.user_id == _up.c.id, isouter=True)
        j = j.join(_au, _up.c.user_id == _au.c.id, isouter=True)
        
        j = j.join(
            lab_head_table, 
            _su.c.lab_head_id==lab_head_table.c.screensaver_user_id,
            isouter=True)
        
        stmt = select(columns.values()).select_from(j)
        # natural ordering
        stmt = stmt.order_by(_su.c.last_name, _su.c.first_name)
            
        if self._meta.authorization.is_restricted_view(request.user):
            logger.info('create_authorized_user_filter')
            associated_users = \
                self._meta.authorization.get_associated_users(request.user)
            stmt = stmt.where(
                _su.c.screensaver_user_id.in_([
                    user.screensaver_user_id for user in associated_users]))
        
        # general setup
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        # logger.info(
        #     'stmt: %s',
        #     str(stmt.compile(
        #         dialect=postgresql.dialect(),
        #         compile_kwargs={"literal_binds": True})))

        title_function = None
        if use_titles is True:
            def title_function(key):
                return field_hash[key]['title']
        if is_data_interchange:
            title_function = DbApiResource.datainterchange_title_function(
                field_hash,schema['id_attribute'])
        
        return self.stream_response_from_statement(
            request, stmt, count_stmt, filename,
            field_hash=field_hash,
            param_hash=param_hash,
            is_for_detail=is_for_detail,
            rowproxy_generator=rowproxy_generator,
            title_function=title_function, meta=kwargs.get('meta', None),
            use_caching=True)
        
    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_list must be implemented')
                
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        ScreensaverUser.objects.get(**id_kwargs).delete()

    def get_id(self, deserialized, validate=False, schema=None, **kwargs):
        
        logger.debug('su get_id: %r, %r', deserialized, kwargs)
        id_kwargs = {}
        
        screensaver_user_id = kwargs.get('screensaver_user_id', None)
        if screensaver_user_id is None:
            screensaver_user_id = deserialized.get('screensaver_user_id')
        if screensaver_user_id is not None:
            screensaver_user_id = parse_val(
                screensaver_user_id, 'screensaver_user_id', 'integer')
            if screensaver_user_id is not None:
                id_kwargs['screensaver_user_id'] = screensaver_user_id
        if not id_kwargs:
            username = kwargs.get('username', None)
            if username is None:
                username = deserialized.get('username', None)
            if username:
                username = username.strip()
                if len(username) > 0:
                    id_kwargs['username'] = username
        if not id_kwargs:
            ecommons_id = kwargs.get('ecommons_id', None)
            if ecommons_id is None:
                ecommons_id = deserialized.get('ecommons_id', None)
            if ecommons_id:
                ecommons_id = ecommons_id.strip()
                if len(ecommons_id) > 0:
                    id_kwargs['ecommons_id'] = ecommons_id
            
        if not id_kwargs:
            if validate is True:
                raise ValueError, (
                    'Neither screensaver_user_id, username nor ecommons_id'
                    ' were specified')
            
        logger.info('su get_id: result: %r', id_kwargs)
        return id_kwargs

    @write_authorization
    @un_cache
    @transaction.atomic        
    def post_detail(self, request, **kwargs):
        '''
        Special POST to allow for screensaver_user_id generation
        '''
        logger.info('post_detail, screensaveruser')
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs.pop('data')
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        logger.debug('patch detail %s, %s', deserialized,kwargs)
        
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.debug('param_hash: %r', param_hash)

        id_kwargs = self.get_id(deserialized, validate=False, schema=schema, **kwargs)
        
        if not id_kwargs:
           
            if 'first_name' not in deserialized:
                raise ValidationError(key='first_name', msg='required')
            if 'last_name' not in deserialized:
                raise ValidationError(key='last_name', msg='required')
            
            try:
                first_name=parse_val(
                    deserialized['first_name'],'first_name', 'string')
                last_name = parse_val(
                    deserialized['last_name'],'last_name', 'string')
                extant_users = ScreensaverUser.objects.filter(
                    first_name__iexact=first_name,
                    last_name__iexact=last_name)
                if extant_users.exists():
                    ids = [su.screensaver_user_id for su in extant_users.all()]
                    msg = ('Screensaver User: %r has already been created using '
                        'the given first and last names: %s %s'
                        % (ids, first_name,last_name ))
                    raise ValidationError({
                        'first_name': msg, 'last_name': msg })
            except ObjectDoesNotExist:
                logger.info('ok to create new user: %r', deserialized)
        else:
            try:
                extant_user = ScreensaverUser.objects.get(**id_kwargs)
                msg = 'User %r, identified by %r already exists'
                raise ValidationError({
                    key:msg%(extant_user.screensaver_user_id, value) 
                        for key, value in id_kwargs.items() 
                })
                
            except ObjectDoesNotExist:
                logger.info('ok to create new user %r', id_kwargs)
                
        return super(ScreensaverUserResource,self).post_detail(
            request, data=deserialized, **kwargs)
    
    @write_authorization
    @un_cache
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):

        DEBUG_USER_PATCH = False or logger.isEnabledFor(logging.DEBUG)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        screensaver_user = None
        is_patch = False
        if 'screensaver_user_id' in id_kwargs:
            screensaver_user = ScreensaverUser.objects.get(
                screensaver_user_id=id_kwargs['screensaver_user_id'])
        else:
            try:
                if 'username' in id_kwargs:
                    screensaver_user = ScreensaverUser.objects.get(
                        username=id_kwargs['username'])
                elif 'ecommons_id' in id_kwargs:
                    screensaver_user = ScreensaverUser.objects.get(
                        ecommons_id=id_kwargs['ecommons_id'])
                else:
                    logger.info('creating non-login user: %r', deserialized)
            except ObjectDoesNotExist:
                is_patch=False
        if screensaver_user:
            is_patch = True
            logger.info('patching screensaveruser: %r: %r',
                deserialized, screensaver_user)
        else:
            logger.info('creating Screensaver User: %r', deserialized)
            
        errors = self.validate(deserialized, schema=schema,patch=is_patch)
        if errors:
            raise ValidationError(errors)
        
        screensaveruser_fields = { name:val for name, val in fields.items() 
            if val['scope'] == 'fields.screensaveruser'}
        logger.debug(
            'fields.screensaveruser fields: %s', screensaveruser_fields.keys())
        initializer_dict = self.parse(
            deserialized, create=not is_patch, fields=screensaveruser_fields)
        
        logger.debug('initializer_dict: %r', initializer_dict)
        
        if screensaver_user is None:
            if id_kwargs:
                logger.info('Create log in Screensaver User for: %r', id_kwargs)
                screensaver_user = ScreensaverUser(**id_kwargs)
            else:
                logger.info(
                    'Create non-log in Screensaver User for: %r %r',
                    initializer_dict['first_name'],initializer_dict['last_name'])
                screensaver_user = ScreensaverUser(
                    first_name=initializer_dict['first_name'],
                    last_name=initializer_dict['last_name'])
            logger.info('saving new user...')
            screensaver_user.save()
            logger.info('new user saved: %r: %r', 
                screensaver_user, screensaver_user.date_created)
        messages = []
        if initializer_dict:
            # Validation rules for username, ecommons_id
            new_username = initializer_dict.pop('username',None)
            if new_username is not None:
                new_username = new_username.strip()
                if len(new_username)==0:
                    new_username = None
            new_ecommons = initializer_dict.pop('ecommons_id',None)
            if new_ecommons is not None:
                new_ecommons = new_ecommons.strip()
                if len(new_ecommons)==0:
                    new_ecommons = None
            logger.info('new_ecommons: %r', new_ecommons)
            if new_username is not None and new_ecommons is not None:
                if new_username != new_ecommons:
                    raise ValidationError({
                        'username': 'does not match ecommons_id',
                        'ecommons_id': 'does not match username' })
            if new_username is not None:
                if screensaver_user.username is not None:
                    if screensaver_user.username != new_username:
                        raise ValidationError(key='username', msg='immutable')
                else:
                    screensaver_user.username = new_username
            if new_ecommons is not None:
                if screensaver_user.ecommons_id is not None:
                    if screensaver_user.ecommons_id != new_ecommons:
                        raise ValidationError(key='username', msg='immutable')
                    if screensaver_user.username is not None:
                        if screensaver_user.username != new_ecommons:
                            raise ValidationError(
                                key='username', 
                                msg='immutable (via ecommons_id)')
                    else:
                        screensaver_user.username = new_ecommons
                else:
                    logger.info('screensaver_user: %r, setting ecommons: %r', 
                        screensaver_user.screensaver_user_id, new_ecommons)
                    screensaver_user.ecommons_id = new_ecommons
                    if screensaver_user.username is None:
                        screensaver_user.username=screensaver_user.ecommons_id
                    elif screensaver_user.username != new_ecommons:
                        raise ValidationError(
                            key='username', 
                            msg='immutable (via ecommons_id)')
            
            # Validation rules for first_name, last_name
            first_name = initializer_dict.pop('first_name', None)
            last_name = initializer_dict.pop('last_name', None)
            
            if first_name is not None or last_name is not None:
                if last_name is None:
                    last_name = screensaver_user.last_name
                if first_name is None:
                    first_name = screensaver_user.first_name
                if first_name != screensaver_user.first_name \
                    or last_name != screensaver_user.last_name:
                    extant_users = ScreensaverUser.objects.filter(
                        first_name__iexact=first_name,
                        last_name__iexact=last_name)
                    if extant_users.exists():
                        ids = [su.screensaver_user_id for su in extant_users.all()]
                        msg = ('Screensaver User: %r has already been created using '
                            'the given first and last names: %s %s'
                            % (ids, first_name,last_name ))
                        raise ValidationError({
                            'first_name': msg, 'last_name': msg })
                    screensaver_user.first_name = first_name
                    screensaver_user.last_name = last_name

            key = 'lab_head_id'
            lab_head_id = initializer_dict.get(key,None)
            if lab_head_id:
                try:
                    lab_head = ScreensaverUser.objects.get(
                        screensaver_user_id=lab_head_id)
                    if lab_head.classification != VOCAB_USER_CLASSIFICATION_PI:
                        raise ValidationError(
                            'Chosen lab head "user.classification must be %r '
                            % VOCAB_USER_CLASSIFICATION_PI)
                except ObjectDoesNotExist, e:
                    raise ValidationError(
                        key=_key,
                        msg='Screensaver user for %r not found' % lh_id)
            
            key = 'classification'
            classification = initializer_dict.get(key,None)
            if classification:
                if classification == VOCAB_USER_CLASSIFICATION_PI:
                    if lab_head_id is not None:
                        raise ValidationError(
                            key='lab_head_id',
                            msg='Classification may not be % for lab member'
                                % VOCAB_USER_CLASSIFICATION_PI)
                elif screensaver_user.classification == VOCAB_USER_CLASSIFICATION_PI:
                    raise ValidationError(
                        key=key, msg='May not be changed from %r'
                            % VOCAB_USER_CLASSIFICATION_PI)
                    
            for key, val in initializer_dict.items():
                if hasattr(screensaver_user, key):
                    setattr(screensaver_user, key, val)
        else:
            logger.info(
                'no (basic) screensaver_user fields to update %s',
                deserialized)
        logger.info('User saved: %r: %r', 
            screensaver_user, screensaver_user.date_created)
        
        # Tie ScreensaverUser entry to Reports.User:
        # If reports_userprofile exists, or username is set, patch the 
        # reports_userprofile
        reports_kwargs = {}
        if screensaver_user.username is not None:
            reports_kwargs['username'] = screensaver_user.username
        user = None
        if reports_kwargs:
            # create/get userprofile
            # NOTE: reports.user patch data must include first/last names
            if 'first_name' not in deserialized or 'last_name' not in deserialized:
                deserialized.update({
                    'first_name': screensaver_user.first_name,
                    'last_name': screensaver_user.last_name})
            
            deserialized['username'] = screensaver_user.username
            if DEBUG_USER_PATCH:
                logger.info('patch the reports user: %r, %r',
                    reports_kwargs, deserialized)
            patch_response = self.get_user_resource().patch_obj(
                request, deserialized, **reports_kwargs)
            logger.info('patched userprofile %s', patch_response)
            user = patch_response[API_RESULT_OBJ]

        is_staff = False
        if user:
            if screensaver_user.user is None:
                if DEBUG_USER_PATCH:
                    logger.info('set the reports userprofile: %r to the su: %r', 
                        user, screensaver_user)
                screensaver_user.user = user
            else:
                if screensaver_user.user != user:
                    raise ValidationError(
                        key='username',
                        msg='ss user found: %r, not equal to current: %r'\
                            % (user, screensaver_user.user))    
            if DEBUG_USER_PATCH:
                logger.info('set the username: %r: %r, %r', 
                    screensaver_user.screensaver_user_id, 
                    screensaver_user.username, user.username)
            screensaver_user.username = user.username

            auth_user = user.user
            is_staff = auth_user.is_staff

        logger.info('1-b user saved: %r: %r', 
            screensaver_user, screensaver_user.date_created)
        
        screensaver_user.save()
        
        # Validate business rules for user types
        if is_staff != True:
            usergroups = deserialized.get('usergroups', None)
            if usergroups:
                raise ValidationError(
                    key='usergroups',
                    msg='May only be set for staff users')
            permissions = deserialized.get('permissions', None)
            if usergroups:
                raise ValidationError(
                    key='permissions',
                    msg='May only be set for staff users')

            if screensaver_user.classification == VOCAB_USER_CLASSIFICATION_PI:
                if screensaver_user.lab_affiliation is None:
                    # should already be set
                    raise ValidationError(
                        key='lab_affiliation_id',
                        msg='Must be set if classification is %r'
                            % VOCAB_USER_CLASSIFICATION_PI)
                if screensaver_user.lab_head != screensaver_user:
                    screensaver_user.lab_head = screensaver_user
            
            if screensaver_user.classification != VOCAB_USER_CLASSIFICATION_PI:
                
                lab_head = screensaver_user.lab_head
                if lab_head is None:
                    raise ValidationError({
                        'lab_head_id': 'Must be set for non staff',
                        'classification': 
                            'Must be %r if non staff and not a lab member'
                                % VOCAB_USER_CLASSIFICATION_PI })
                
                lab_head_user_agreements = \
                    UserAgreement.objects.all().filter(screensaver_user=lab_head)
                
                for ua in screensaver_user.useragreement_set.all():
                    if lab_head_user_agreements.exists():
                        for lhua in lab_head_user_agreements.all():
                            if lhua.type == ua.type:
                                if lhua.data_sharing_level != ua.data_sharing_level:
                                    messages.append
                
                # Reset user dsls to lab_head values
                # TODO: warn and override                
#                 if screensaver_user.sm_data_sharing_level \
#                         != lab_head.sm_data_sharing_level:
#                     if lab_head.sm_data_sharing_level is not None:
#                         messages.append('%r updated from %r to %r'
#                             % ( 'sm_data_sharing_level',
#                                 screensaver_user.sm_data_sharing_level,
#                                 lab_head.sm_data_sharing_level ))
#                         screensaver_user.sm_data_sharing_level = \
#                             lab_head.sm_data_sharing_level
#                 if screensaver_user.rnai_data_sharing_level \
#                         != lab_head.rnai_data_sharing_level:
#                     if lab_head.rnai_data_sharing_level is not None:
#                         messages.append('%r updated from %r to %r'
#                             % ( 'rnai_data_sharing_level',
#                                 screensaver_user.rnai_data_sharing_level,
#                                 lab_head.rnai_data_sharing_level ))
#                         screensaver_user.rnai_data_sharing_level = \
#                             lab_head.rnai_data_sharing_level
            
        screensaver_user.save()
        
        _key = 'lab_member_ids'
        if _key in initializer_dict:
            lab_member_ids = initializer_dict[_key]
            lab_members = []
            # may empty
            if lab_member_ids:
                if screensaver_user.classification != VOCAB_USER_CLASSIFICATION_PI:
                    raise ValidationError(
                        key='lab_member_ids',
                        msg='User classification must be %r' % VOCAB_USER_CLASSIFICATION_PI)
            for lab_member_id in lab_member_ids:
                try:
                    lab_member = ScreensaverUser.objects.get(
                        screensaver_user_id=lab_member_id)
                    if lab_member.classification == VOCAB_USER_CLASSIFICATION_PI:
                        raise ValidationError(
                            key='lab_member_ids',
                            msg='User is a %r and cannot be added as a lab member'
                                % VOCAB_USER_CLASSIFICATION_PI)
                    
                    lab_members.append(lab_member)
                    
                    # FIXME: need override from UI to set
                    msg='lab member data sharing level (%r) must match lab head'
                    if screensaver_user.useragreement_set.exists():
                        if lab_member.useragreement_set.exists():
                            for ua in screensaver_user.useragreement_set.all():
                                for member_ua in lab_member.useragreement_set.all():
                                    if ua.type == member_ua.type:
                                        if ua.data_sharing_level is None:
                                            continue
                                        if member_ua.data_sharing_level is None:
                                            continue
                                        if ua.date_expired is not None:
                                            continue
                                        if member_ua.date_expired is not None:
                                            continue
                                        if ua.data_sharing_level != member_ua.data_sharing_level:
                                            raise ValidationError(
                                                key='lab_member_ids',
                                                msg=msg%member_ua.data_sharing_level)
                    
#                     lab_member.sm_data_sharing_level = \
#                         screensaver_user.sm_data_sharing_level
#                     lab_member.rnai_data_sharing_level = \
#                         screensaver_user.rnai_data_sharing_level
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key=_key,
                        msg='No such user: %r' % lab_member_id)
            screensaver_user.lab_members = lab_members
        screensaver_user.save()
                
        if 'facility_usage_roles' in initializer_dict:
            current_roles = set([r.facility_usage_role 
                for r in screensaver_user.userfacilityusagerole_set.all()])
            new_roles = initializer_dict['facility_usage_roles']
            if new_roles is None:
                new_roles = []
            new_roles = set(new_roles)
            if DEBUG_USER_PATCH:
                logger.info('roles to delete: %s', current_roles - new_roles)
            (screensaver_user.userfacilityusagerole_set
                .filter(facility_usage_role__in=current_roles - new_roles)
                .delete())
            for role in new_roles - current_roles:
                if DEBUG_USER_PATCH:
                    logger.info(
                        'create facility_usage_role: %s, %s', 
                        screensaver_user, role)
                ur = UserFacilityUsageRole.objects.create(
                    screensaver_user=screensaver_user,
                    facility_usage_role=role)            
                ur.save()

        logger.info('patch_obj done: %r', screensaver_user)
        
        
        meta = {
            API_MSG_RESULT: {
                'messages': messages
            },
        }
        
        return { 
            API_RESULT_OBJ: screensaver_user,
            API_RESULT_META: meta 
        }
            

class LibraryResourceAuthorization(UserGroupAuthorization):        

    def _is_resource_authorized(
        self, user, permission_type, **kwargs):
        authorized = super(LibraryResourceAuthorization, self)\
            ._is_resource_authorized(user, permission_type, **kwargs)
        if authorized is True:
            return True
        # Allow all libraries for all users at this time.
        if permission_type == 'read':
            return user.is_active
        logger.info('resource not authorized: %r, %r, %r', 
            self.resource_name, user, permission_type)
        return False
    
    def is_restricted_view(self, user):
        return not self._is_resource_authorized(user, 'read')
        
    # NOTE: 20170908 - Not used, all authenticated users may access - JAS/CS
    def get_authorized_library_screen_types(self, user):
        pass
        # screensaver_user = ScreensaverUser.objects.get(username=user.username)
        # 
        # # TODO: review requirements
        # screen_types = set()
        # SM_TYPE = VOCAB_SCREEN_TYPE_SM
        # RNAI_TYPE = 'rnai'
        # if screensaver_user.sm_data_sharing_level is not None:
        #     screen_types.add(SM_TYPE)
        # if screensaver_user.rnai_data_sharing_level is not None:
        #     screen_types.add(SM_TYPE)
        # 
        # return screen_types

class NaturalProductReagentResource(DbApiResource):
    # Consider folding the NaturalProductReagentResource into ReagentResource
    
    class Meta:

        queryset = Reagent.objects.all()
        
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'naturalproductreagent'
        authorization = LibraryResourceAuthorization(resource_name)
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(NaturalProductReagentResource, self).__init__(**kwargs)

    def delete_reagents(self, library):
        
        NaturalProductReagent.objects.all().filter(well__library=library).delete()

    def _patch_wells(self, request, deserialized):
        ''' For bulk update: 
        - deserialized has been loaded with the well
        '''
        schema = self.build_schema(user=request.user)
        for i,well_data in enumerate(deserialized):
            well = well_data['well']
            is_patch = False
            if not well.reagents.exists():
                reagent = NaturalProductReagent(well=well)
            else:
                is_patch = True
                # TODO: only works for a single reagent
                # can search for the reagent using id_kwargs
                # reagent = well.reagents.all().filter(**id_kwargs)
                # TODO: update reagent
                reagent = well.reagents.all()[0]
                reagent = reagent.naturalproductreagent
                logger.debug('found reagent: %r, %r', reagent.well_id, reagent)

            field_type = 'create_fields'
            if is_patch == True:
                field_type = 'update_fields'
            fields = { key:field for key,field in schema['fields'].items()
                if key in schema[field_type] }
                
            initializer_dict = self.parse(well_data, fields=fields)
            errors = self.validate(initializer_dict, patch=is_patch, schema=schema)
            if errors:
                raise ValidationError(errors)
            
            for key, val in initializer_dict.items():
                if hasattr(reagent, key):
                    setattr(reagent, key, val)
    
            reagent.save()
            if i % 999 == 0:
                logger.info('patched %d reagents', i+1)
        logger.info('patched %d reagents', i+1)


class SilencingReagentResource(DbApiResource):
    # NOTE: OO model for inheritance is broken in multiple places:
    # Consider folding the SilencingReagentResource into ReagentResource
    
    class Meta:
    
        queryset = Reagent.objects.all()
        authentication = MultiAuthentication(
            BasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'silencingreagent'
        authorization = LibraryResourceAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):

        super(SilencingReagentResource, self).__init__(**kwargs)
        
        # for debugging
        self.patch_elapsedtime1 = 0
        self.patch_elapsedtime2 = 0
        self.patch_elapsedtime3 = 0
        
    def get_debug_times(self):
        
        logger.info('silencing reagent times: %r, %r, %r', 
            self.patch_elapsedtime1, self.patch_elapsedtime2, self.patch_elapsedtime3)
        self.patch_elapsedtime1 = 0
        self.patch_elapsedtime2 = 0
        self.patch_elapsedtime3 = 0

    def delete_reagents(self, library):
        
        SilencingReagent.objects.all().filter(well__library=library).delete()
        
    def build_sqlalchemy_columns(self, fields):
        '''
        @return an array of sqlalchemy.sql.schema.Column objects
        @param fields - field definitions, from the resource schema
        
        '''
        logger.debug(
            'build silencing reagent columns for: %r', 
            [field['key'] for field in fields])
        DEBUG_BUILD_COLS = False or logger.isEnabledFor(logging.DEBUG)
        
        bridge = self.bridge
        
        columns = {}
        vendor_gene_columns = [
            'vendor_entrezgene_id', 'vendor_gene_name', 'vendor_gene_species']
        vendor_gene_symbols = 'vendor_entrezgene_symbols'
        vendor_genebank_accession_numbers = 'vendor_genbank_accession_numbers'
        facility_gene_columns = [
            'facility_entrezgene_id', 'facility_gene_name',
            'facility_gene_species']
        facility_gene_symbols = 'facility_entrezgene_symbols'
        facility_genebank_accession_numbers = 'facility_genbank_accession_numbers'
        
        vendor_columns = set(vendor_gene_columns)
        vendor_columns.add(vendor_gene_symbols)
        vendor_columns.add(vendor_genebank_accession_numbers)
        
        facility_columns = set(facility_gene_columns)
        facility_columns.add(facility_gene_symbols)
        facility_columns.add(facility_genebank_accession_numbers)
        
        gene_table = bridge['gene']
        sirna_table = bridge['silencing_reagent']
        gene_symbol = bridge['gene_symbol']
        genbank_acc = bridge['gene_genbank_accession_number']
        well_table = bridge['well']
        
        # Example:
        # (select gene_name from gene 
        #     join silencing_reagent on(gene_id=vendor_gene_id) 
        #     where silencing_reagent.reagent_id = reagent.reagent_id ) as vendor_gene_name
        
        for field in fields:
            field_name = field.get('field', None)
            if not field_name:
                field_name = field['key']
            key = field['key']
            if DEBUG_BUILD_COLS: 
                logger.info('field[key]: %r, %r, %r', field['key'])
            join_stmt = None
            join_column = None
            if field['key'] in vendor_columns:
                join_column = 'vendor_gene_id'
            if field['key'] in facility_columns:
                join_column = 'facility_gene_id'

            if key in vendor_gene_columns or \
                    key in facility_gene_columns:
                join_stmt = gene_table.join(sirna_table,
                    gene_table.c['gene_id'] == sirna_table.c[join_column])
                select_stmt = select([gene_table.c[field_name]]).\
                    select_from(join_stmt)
                select_stmt = select_stmt.where(
                    text('silencing_reagent.reagent_id=reagent.reagent_id'))
                select_stmt = select_stmt.label(key)
                columns[key] = select_stmt

            if key == vendor_gene_symbols or \
                    key == facility_gene_symbols:
                join_stmt = gene_symbol.join(gene_table,
                    gene_symbol.c['gene_id'] == gene_table.c['gene_id'])
                join_stmt = join_stmt.join(sirna_table,
                    gene_table.c['gene_id'] == sirna_table.c[join_column])
                
                select_inner = select([gene_symbol.c[field_name]]).\
                    select_from(join_stmt)
                ordinal_field = field.get('ordinal_field', None)
                if ordinal_field:
                    select_inner = select_inner.order_by(
                        gene_symbol.c[ordinal_field])
                select_inner = select_inner.where(
                    text('silencing_reagent.reagent_id=reagent.reagent_id'))
                select_inner = select_inner.alias('a')
                select_stmt = select([
                    func.array_to_string(func.array_agg(
                        column(field_name)),LIST_DELIMITER_SQL_ARRAY)])
                select_stmt = select_stmt.select_from(select_inner)
                select_stmt = select_stmt.label(key)
                columns[key] = select_stmt

            if key == vendor_genebank_accession_numbers or \
                    key == facility_genebank_accession_numbers:
                join_stmt = genbank_acc.join(gene_table,
                    genbank_acc.c['gene_id'] == gene_table.c['gene_id'])
                join_stmt = join_stmt.join(sirna_table,
                    gene_table.c['gene_id'] == sirna_table.c[join_column])
                
                select_inner = select([genbank_acc.c[field_name]]).\
                    select_from(join_stmt)
                select_inner = select_inner.where(
                    text('silencing_reagent.reagent_id=reagent.reagent_id'))
                select_inner = select_inner.alias('a')
                select_stmt = select([func.array_to_string(
                                func.array_agg(column(field_name)),
                                               LIST_DELIMITER_SQL_ARRAY)])
                select_stmt = select_stmt.select_from(select_inner)
                select_stmt = select_stmt.label(key)
                columns[key] = select_stmt
            
            if key == 'duplex_wells':
                _duplex_wells = bridge['silencing_reagent_duplex_wells']
                select_inner = select([_duplex_wells.c['well_id']]).\
                    select_from(_duplex_wells)
                select_inner = select_inner.where(
                    text('silencingreagent_id=reagent.reagent_id'))
                select_inner = select_inner.order_by(_duplex_wells.c['well_id'])
                select_inner = select_inner.alias('a')
                select_stmt = select([func.array_to_string(
                                func.array_agg(column(field_name)),
                                               LIST_DELIMITER_SQL_ARRAY)])
                select_stmt = select_stmt.select_from(select_inner)
                select_stmt = select_stmt.label(key)
                columns[key] = select_stmt
                
                if DEBUG_BUILD_COLS:
                    logger.info(str((select_stmt)))
            
            # NOT Performant - not used, see ReagentResource pool_well definition
            if key == 'pool_well':
                _duplex_wells = bridge['silencing_reagent_duplex_wells']
                pool_reagent = bridge['reagent'].alias('pool_reagent')
                 
                columns[key] = (
                    select([func.array_to_string(
                        func.array_agg(literal_column('well_id')),
                        LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([pool_reagent.c.well_id])
                        .select_from(_duplex_wells.join(
                            pool_reagent,pool_reagent.c.reagent_id==
                                _duplex_wells.c.silencingreagent_id))
                        .where(_duplex_wells.c.well_id==text('reagent.well_id')).alias('inner_pool')
                    )
                # Some duplex wells are pointed to by multiple pool wells:
                # CPR44311, lcp report shows multple rows 
                # returned from subquery
                # e.g.: 50471:G08 - has 50599:I17,50599:I17...
                # (Mouse4 Pools: Remaining Genome)
                # 50448:C18 - has 50599:I17,50599:I17... 
                # (Mouse4 Pools: Remaining Genome)
                )
                
                  
        if DEBUG_BUILD_COLS: 
            logger.info('sirna columns: %r', columns.keys())
        
        return columns 

    # NOTE: removed; patching will only be done in batch, from library/well
    # @transaction.atomic()    
    # def patch_obj(self, request, deserialized, **kwargs):
    # 
    #     start_time = time.time()
    #      
    #     well = kwargs.get('well', None)
    #     if not well:
    #         raise ValidationError(key='well', msg='required')
    # 
    #     schema = self.build_schema()
    #     fields = schema['fields']
    #     mutable_fields = []
    #     for field in fields.values():
    #         editability = field.get('editability', None)
    #         if editability and (
    #             'u' in editability or (create and 'c' in editability )):
    #             mutable_fields.append(field)
    #     id_kwargs = self.get_id(deserialized, **kwargs)
    # 
    #     patch = False
    #     if not well.reagents.exists():
    #         reagent = SilencingReagent(well=well)
    #          
    #         # only allow duplex_wells to be set on create
    #         if 'duplex_wells' in kwargs:
    #             reagent.save()
    #             reagent.duplex_wells = kwargs['duplex_wells']
    #      
    #     else:
    #         patch = True
    #         # TODO: only works for a single reagent
    #         # can search for the reagent using id_kwargs
    #         # reagent = well.reagents.all().filter(**id_kwargs)
    #         # TODO: update reagent
    #         reagent = well.reagents.all()[0]
    #         reagent = reagent.silencingreagent
    #         logger.info('found reagent: %r, %r', reagent.well_id, reagent)
    # 
    #     self._set_reagent_values(reagent, deserialized, patch)                
    #      
    #     return reagent        

    def _patch_wells(self, request, deserialized):
        ''' For bulk update: 
        - deserialized has been loaded with the well & duplex wells
        '''
        schema = self.build_schema(user=request.user)
        for i,well_data in enumerate(deserialized):
            well = well_data['well']
            is_patch = False
            if not well.reagents.exists():
                reagent = SilencingReagent(well=well)
                
                # only allow duplex_wells to be set on create
                if well_data.get('duplex_wells', None):
                    reagent.save()
                    reagent.duplex_wells = well_data['duplex_wells']
            else:
                is_patch = True
                # TODO: only works for a single reagent
                # can search for the reagent using id_kwargs
                # reagent = well.reagents.all().filter(**id_kwargs)
                # TODO: update reagent
                reagent = well.reagents.all()[0]
                reagent = reagent.silencingreagent
                logger.debug('found reagent: %r, %r', reagent.well_id, reagent)

            field_type = 'create_fields'
            if is_patch == True:
                field_type = 'update_fields'
            fields = { key:field for key,field in schema['fields'].items()
                if key in schema[field_type] }
                
            self._set_reagent_values(reagent, well_data, is_patch, schema, fields)

            if i % 999 == 0:
                logger.info('patched %d reagents', i+1)
        logger.info('patched %d reagents', i+1)

    def _set_reagent_values(self, reagent, deserialized, is_patch, schema, fields):

        start_time = time.time()

        initializer_dict = self.parse(deserialized, fields=fields)
        errors = self.validate(initializer_dict, patch=is_patch, schema=schema)
        if errors:
            raise ValidationError(errors)

        self.patch_elapsedtime1 += (time.time() - start_time)
        start_time = time.time()
            
        for key, val in initializer_dict.items():
            if hasattr(reagent, key):
                setattr(reagent, key, val)
        reagent.save()
        logger.debug('patch silencing reagent: %r', reagent)

        self.patch_elapsedtime2 += (time.time() - start_time)
        start_time = time.time()

        # Now do the gene tables
        
        gene_key = 'entrezgene_id'
        if deserialized.get('vendor_%s' % gene_key, None):
            reagent.vendor_gene = \
                self._create_gene(deserialized, 'vendor')
        if deserialized.get('facility_%s' % gene_key, None):
            reagent.facility_gene = \
                self._create_gene(deserialized, 'facility')
        reagent.save()

        self.patch_elapsedtime3 += (time.time() - start_time)

    def _create_gene(self, data, source_type):
        
        gene_keys = ['entrezgene_id', 'gene_name', 'species_name']
        gene = Gene()
        for key in gene_keys:
            api_key = '%s_%s' % (source_type, key)
            val = data.get(api_key, None)
            if val:
                setattr(gene, key, val)
        gene.save()
        
        _key = 'entrezgene_symbols'
        if data.get('%s_%s' % (source_type, _key), None):
            symbol_list = data['%s_%s' % (source_type, _key)]
            for i, symbol in enumerate(symbol_list):
                gene_symbol = GeneSymbol()
                setattr(gene_symbol, 'entrezgene_symbol', symbol)
                setattr(gene_symbol, 'ordinal', i)
                setattr(gene_symbol, 'gene', gene)
                gene_symbol.save()
    
        _key = 'genbank_accession_numbers'
        if data.get('%s_%s' % (source_type, _key), None):
            _list = data['%s_%s' % (source_type, _key)]
            for i, num in enumerate(_list):
                accession_number = GeneGenbankAccessionNumber()
                setattr(accession_number, 'genbank_accession_number', num)
                setattr(accession_number, 'gene', gene)
                accession_number.save()
        
        return gene
        

class SmallMoleculeReagentResource(DbApiResource):
    # NOTE: OO model for inheritance is broken in multiple places:
    # Consider folding the SmallMoleculeReagentResource into ReagentResource
        
    class Meta:

        queryset = Reagent.objects.all() 
        authentication = MultiAuthentication(
            BasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'smallmoleculereagent' 
        authorization = LibraryResourceAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = []  # ['json_field']
        always_return_data = True  # this makes Backbone happy

    def __init__(self, **kwargs):
        super(SmallMoleculeReagentResource, self).__init__(**kwargs)
        # for debugging
        self.patch_elapsedtime1 = 0
        self.patch_elapsedtime2 = 0
        self.patch_elapsedtime3 = 0
        
    def get_debug_times(self):
        
        logger.info('sm reagent times: %r, %r, %r', 
            self.patch_elapsedtime1, self.patch_elapsedtime2, self.patch_elapsedtime3)
        self.patch_elapsedtime1 = 0
        self.patch_elapsedtime2 = 0
        self.patch_elapsedtime3 = 0

    def build_sqlalchemy_columns(self, fields):
        DEBUG_BUILD_COLS = False or logger.isEnabledFor(logging.DEBUG)

        bridge = self.bridge
        
        columns = {}
        return columns
#         _smcn = self.bridge['small_molecule_compound_name']
#         columns['compound_name'] = (
#             select([func.array_to_string(
#                 func.array_agg(literal_column('compound_name')
#                 ), LIST_DELIMITER_SQL_ARRAY)])
#                 .select_from(
#                     select([_smcn.c.compound_name])
#                         .select_from(_smcn)
#                         .order_by(_smcn.c.ordinal)
#                         .where(_smcn.c.reagent_id 
#                             == literal_column('reagent.reagent_id'))
#                         .alias('inner'))
#             )
#                   
#         if DEBUG_BUILD_COLS: 
#             logger.info('sm columns: %r', columns.keys())
#         
#         return columns 

        
    def delete_reagents(self, library):
        
        SmallMoleculeReagent.objects.all().filter(well__library=library).delete()

    def _patch_wells(self, request, deserialized):
        ''' For bulk update: 
        - deserialized has been loaded with the wells
        '''
        logger.info('patch reagents...')
        
        schema = self.build_schema(user=request.user)
        for i,well_data in enumerate(deserialized):
            well = well_data['well']
            is_patch = False
            if not well.reagents.exists():
                reagent = SmallMoleculeReagent(well=well)
                reagent.save()
            else:
                is_patch = True
                # TODO: only works for a single reagent
                # can search for the reagent using id_kwargs
                # reagent = well.reagents.all().filter(**id_kwargs)
                # TODO: update reagent
                reagent = well.reagents.all()[0]
                reagent = reagent.smallmoleculereagent
                logger.debug('found reagent: %r, %r', reagent.well_id, reagent)
            
            field_type = 'create_fields'
            if is_patch == True:
                field_type = 'update_fields'
            fields = { key:field for key,field in schema['fields'].items()
                if key in schema[field_type] }
                
            self._set_reagent_values(reagent, well_data, is_patch, schema, fields)
            
            if i and i % 1000 == 0:
                logger.info('patched %d reagents', i+1)
        logger.info('patched %d reagents', i+1)

    def _set_reagent_values(self, reagent, deserialized, is_patch, schema, fields):

        start_time = time.time()

        initializer_dict = self.parse(deserialized, fields=fields)
        errors = self.validate(initializer_dict, patch=is_patch, schema=schema)
        if errors:
            raise ValidationError(errors)

        self.patch_elapsedtime1 += (time.time() - start_time)
        start_time = time.time()
        
        for key, val in initializer_dict.items():
            if hasattr(reagent, key):
                setattr(reagent, key, val)
        logger.debug('patch small molecule reagent: %r', reagent)
        reagent.save()

        self.patch_elapsedtime2 += (time.time() - start_time)
        start_time = time.time()
        
        if 'compound_name' in initializer_dict:
            reagent.smallmoleculecompoundname_set.all().delete()
            # TODO: does this delete the old name entries?
            values = initializer_dict['compound_name'] or []
            for ordinal, val in enumerate(values):
                cn = SmallMoleculeCompoundName.objects.create(
                    reagent=reagent, compound_name=val, ordinal=ordinal)
                cn.save()
        if 'chembank_id' in initializer_dict:
            reagent.smallmoleculechembankid_set.all().delete()
            values = initializer_dict['chembank_id'] or []
            for id in values:
                cid = SmallMoleculeChembankId.objects.create(
                    reagent=reagent, chembank_id=id)
                cid.save()
        if 'pubchem_cid' in initializer_dict:
            reagent.smallmoleculepubchemcid_set.all().delete()
            values = initializer_dict['pubchem_cid'] or []
            for id in values:
                cid = SmallMoleculePubchemCid.objects.create(
                    reagent=reagent, pubchem_cid=id)
                cid.save()
        if 'chembl_id' in initializer_dict:
            reagent.smallmoleculechemblid_set.all().delete()
            values = initializer_dict['chembl_id'] or []
            for id in values:
                cid = SmallMoleculeChemblId.objects.create(
                    reagent=reagent, chembl_id=id)
                cid.save()
        if deserialized.get('molfile', None):
            Molfile.objects.all().filter(reagent=reagent).delete()
             
            molfile = Molfile.objects.create(
                reagent=reagent, molfile=deserialized['molfile'])
            molfile.save()

        self.patch_elapsedtime3 += (time.time() - start_time)
                
        logger.debug('sm patch_obj done')
        return reagent
    

class ReagentResource(DbApiResource):
    
    class Meta:

        queryset = Reagent.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'reagent'
        authorization = LibraryResourceAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):

        self.library_resource = None
        self.sr_resource = None
        self.smr_resource = None
        self.npr_resource = None
        self.well_resource = None
        super(ReagentResource, self).__init__(**kwargs)
    
    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('search'), name="api_search"),
            url(r"^(?P<resource_name>%s)/(?P<substance_id>[^:]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]
    
    def get_debug_times(self):
        
        self.get_sr_resource().get_debug_times()
        self.get_smr_resource().get_debug_times()
        
    def get_sr_resource(self):
        if not self.sr_resource:
            self.sr_resource = SilencingReagentResource()
        return self.sr_resource
    
    def get_smr_resource(self):
        if not self.smr_resource:
            self.smr_resource = SmallMoleculeReagentResource()
        return self.smr_resource
    
    def get_npr_resource(self):
        if not self.npr_resource:
            self.npr_resource = NaturalProductReagentResource()
        return self.npr_resource
    
    def get_well_resource(self):
        if not self.well_resource:
            self.well_resource = WellResource()
        return self.well_resource
    
    def get_reagent_resource(self, library_classification):
        if library_classification == 'rnai':
            return self.get_sr_resource()
        else:
            if library_classification == 'natural_products':
                return self.get_npr_resource()
            else:
                return self.get_smr_resource()
    
    def get_library_resource(self):
        if not self.library_resource:
            self.library_resource = LibraryResource()
        return self.library_resource
                
    def get_schema(self, request, **kwargs):
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.info('param hash: %r', param_hash)
        
        if not 'library_short_name' in kwargs:
            return self.build_response(
                request, 
                self.build_schema(user=request.user,**param_hash), **kwargs)
        
        library_short_name = kwargs.pop('library_short_name')
        try:
            library = Library.objects.get(short_name=library_short_name)
            return self.build_response(
                request, 
                self.build_schema(
                    user=request.user,library_classification=library.classification), 
                **kwargs)
            
        except Library.DoesNotExist, e:
            raise Http404(
                'Can not build schema - library def needed'
                'no library found for short_name: %r' % library_short_name)
                
    def build_schema(self, library_classification=None, user=None, **kwargs):
        logger.info('build reagent schema for library_classification: %r',
            library_classification)
        schema = deepcopy(super(ReagentResource, self).build_schema(user=user))
        
#         for library_classification in ('small_molecule', 'rnai'):
#             sub_data = self.get_reagent_resource(
#                 library_classification).build_schema(user=user)
#             logger.debug('sub_schema: %r', sub_data['fields'].keys())
#             newfields = {}
#             newfields.update(sub_data['fields'])
#             newfields.update(schema['fields'])
#             schema['fields'] = newfields
#             
#             for k, v in schema.items():
#                 if k != 'fields' and k in sub_data:
#                     schema[k] = sub_data[k]
            
#         if library_classification is not None:
#             sub_data = self.get_reagent_resource(
#                 library_classification).build_schema(user=user)
#             logger.debug('sub_schema: %r', sub_data['fields'].keys())
#             newfields = {}
#             newfields.update(sub_data['fields'])
#             newfields.update(schema['fields'])
#             schema['fields'] = newfields
#             
#             for k, v in schema.items():
#                 if k != 'fields' and k in sub_data:
#                     schema[k] = sub_data[k]
#         if 'search' in kwargs:
        if True:
            # Build the full schema for search
            # FIXME could determine the schema as in ReagentResource.build_list_response
            sub_data = \
                self.get_reagent_resource(library_classification='small_molecule')\
                .build_schema(user=user)

            newfields = {}
            newfields.update(sub_data['fields'])
            
            sub_data = \
                self.get_reagent_resource(library_classification='rnai')\
                .build_schema(user=user)
            newfields.update(sub_data['fields'])
            # all sub-fields are set not visible
            for k,field in newfields.items():
                if 'l' in field['visibility']:
                    field['visibility'].remove('l')
                if 'd' in field['visibility']:
                    field['visibility'].remove('d')
            newfields.update(schema['fields'])
            schema['fields'] = newfields
            
        well_schema = WellResource().build_schema(user=user)
        schema['fields'].update(well_schema['fields'])
        logger.info('reagent schema built')
        return schema

    def get_list(self, request, param_hash={}, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def get_query(self, 
        param_hash, user, library_classification=None, library=None,
        cherry_pick_request_id_screener=None,
        cherry_pick_request_id_lab=None, wells=None
        ):
        logger.info('get query for library_classification %r', library_classification )
        schema = self.build_schema(
            user=user, library_classification=library_classification)

        try:
            
            # general setup
            manual_field_includes = set(param_hash.get('includes', []))
            
            plate_numbers = param_hash.pop('plate_number__in',None)

            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            
            if filter_expression is None and wells is None:
                raise InformationError(
                    key='Input filters ',
                    msg='Please enter a filter expression to begin')
            filter_expression = \
                self._meta.authorization.filter(user,filter_expression)
                 
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            logger.debug('field hash scopes: %r',
                set([field.get('scope', None) 
                    for field in field_hash.values()]))
            if library_classification:
                # NOTE: breaks OO inheritance: include fields for subclasses
                default_fields = ['fields.well', 'fields.reagent']
                if library_classification == 'rnai':
                    default_fields.append('fields.silencingreagent')
                elif library_classification == 'natural_products':
                    default_fields.append('fields.naturalproductreagent')
                else:
                    default_fields.append('fields.smallmoleculereagent')
                
                _temp = { key:field for key, field in field_hash.items() 
                    if field.get('scope', None) in default_fields }
                field_hash = _temp
                logger.debug('final field hash: %r', field_hash.keys())
            else:
                # consider limiting fields available
                pass
            
            logger.debug('visible fields: %r', field_hash.keys())
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            # specific setup 
        
            base_query_tables = ['well', 'reagent', 'library']

            _well = self.bridge['well']
            _reagent = self.bridge['reagent']
            _library = self.bridge['library']
            _scp = self.bridge['screener_cherry_pick']
            _lcp = self.bridge['lab_cherry_pick']

            j = _well.join(
                _reagent, _well.c.well_id == _reagent.c.well_id, isouter=True)
            j = j.join(_library, _well.c.library_id == _library.c.library_id)
            
            custom_columns = {}
                        
            if 'pool_well' in field_hash:
                # NOTE: breaks OO inhertance composition
                # Required for performance
                _duplex_wells = self.bridge['silencing_reagent_duplex_wells']
                _pool_reagent = self.bridge['reagent']
                pool_wells = (
                    select([
                        _pool_reagent.c.well_id.label('pool_well_id'),
                        _duplex_wells.c.well_id
                        ])
                    .select_from(_duplex_wells.join(
                        _pool_reagent, _pool_reagent.c.reagent_id==
                            _duplex_wells.c.silencingreagent_id))
                    ).cte('pool_wells')
                    
                # NOTE: array_agg must be used
                # Some duplex wells are pointed to by multiple pool wells:
                # CPR44311, lcp report shows multple rows 
                # returned from subquery
                # e.g.: 50471:G08 - has 50599:I17,50599:I17...
                # (Mouse4 Pools: Remaining Genome)
                # 50448:C18 - has 50599:I17,50599:I17... 
                # (Mouse4 Pools: Remaining Genome)

                # j = j.join(
                #     pool_wells, pool_wells.c.well_id==_well.c.well_id,
                #     isouter=True )                                
                custom_columns['pool_well'] = (
                    select([func.array_to_string(
                        func.array_agg(pool_wells.c.pool_well_id),
                        LIST_DELIMITER_SQL_ARRAY)])
                        .select_from(pool_wells)
                        .where(pool_wells.c.well_id==_well.c.well_id))
            
            if cherry_pick_request_id_screener is not None:
                j = j.join(_scp,
                    _scp.c.screened_well_id==_well.c.well_id)
            if cherry_pick_request_id_lab is not None:
                j = j.join(_lcp,
                    _lcp.c.source_well_id==_well.c.well_id)
    
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            
            stmt = select(columns.values()).select_from(j)
            if wells is not None:
                stmt = stmt.where(_well.c.well_id.in_([w.well_id for w in wells]))
            if plate_numbers:
                stmt = stmt.where(_well.c.plate_number.in_(plate_numbers))
            if library:
                stmt = stmt.where(_well.c.library_id == library.library_id) 
            if cherry_pick_request_id_screener is not None:
                stmt = stmt.where(
                    _scp.c.cherry_pick_request_id
                        ==cherry_pick_request_id_screener)
            if cherry_pick_request_id_lab is not None:
                stmt = stmt.where(
                    _lcp.c.cherry_pick_request_id
                        ==cherry_pick_request_id_lab)
                
            # Could use the library param to limit the column building exercise
            # to the sub-resource, but since all columns can be joined, just include
            # the SR columns, as above.            
            # sub_resource = None
            # if library:
            #     sub_resource = self.get_reagent_resource(library_screen_type=library.screen_type)
            # if sub_resource and hasattr(sub_resource, 'build_sqlalchemy_columns'):
            #     sub_columns = sub_resource.build_sqlalchemy_columns(
            #         field_hash.values(), self.bridge)
            #     if DEBUG_GET_LIST: 
            #         logger.info(str(('sub_columns', sub_columns.keys())))
            #     columns = self.build_sqlalchemy_columns(
            #         field_hash.values(), base_query_tables=base_query_tables,
            #         custom_columns=sub_columns)
            # else:
            # 
            #     sub_columns = sub_resource.build_sqlalchemy_columns(
            #         field_hash.values(), self.bridge)
            #     if DEBUG_GET_LIST: 
            #         logger.info(str(('sub_columns', sub_columns.keys())))
            #     columns = self.build_sqlalchemy_columns(
            #         field_hash.values(), base_query_tables=base_query_tables,
            #         custom_columns=sub_columns)
            #      
                

            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            if not order_clauses:
                stmt = stmt.order_by("plate_number", "well_name")

#             compiled_stmt = str(stmt.compile(
#                 dialect=postgresql.dialect(),
#                 compile_kwargs={"literal_binds": True}))
#             logger.info('compiled_stmt %s', compiled_stmt)
    
            return (field_hash, columns, stmt, count_stmt,filename)
                        
        except Exception, e:
            logger.exception('on get list')
            raise e  
    
    def build_sqlalchemy_columns(self, fields, base_query_tables=None, 
            custom_columns=None):
        sub_columns = self.get_sr_resource().build_sqlalchemy_columns(fields)
        sub_columns.update(self.get_smr_resource().build_sqlalchemy_columns(fields))
#         sub_columns['plate_number'] = (literal_column(
#             "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
#             .label('plate_number'))
        
        if custom_columns is not None:
            sub_columns.update(custom_columns)
        logger.debug('final columns: %r', sub_columns.keys())
        return DbApiResource.build_sqlalchemy_columns(self, 
            fields, base_query_tables=base_query_tables, 
            custom_columns=sub_columns)
        
    @read_authorization
    def build_list_response(self, request, **kwargs):
        
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        param_hash.pop('schema', None)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        library_classification = None
        
        # well search data is raw line based text entered by the user
        well_search_data = param_hash.pop('raw_search_data', None)
        wells = None
        if well_search_data is not None:
            (wells,errors) = WellResource.find_wells(well_search_data)
            logger.info('well search errors: %r', errors)
            if errors:
                raise ValidationError(
                    key='well_search',
                    msg='%r' % errors )    
            classifications = set([w.library.classification for w in wells])
            logger.info('classifications: %r', classifications)
            if wells and len(classifications) == 1:
                for well in wells:
                    library_classification = well.library.classification
                    logger.info('c: %r', library_classification)
                    break
            logger.info('found wells: %d', len(wells))
        # TODO: eliminate dependency on library (for schema determination)
        library = None
        library_short_name = param_hash.pop('library_short_name', None)
        if not library_short_name:
            logger.info('no library_short_name provided')
        else:
            param_hash['library_short_name__eq'] = library_short_name
            library = Library.objects.get(short_name=library_short_name)
        
        well_id = param_hash.pop('well_id', None)
        if well_id:
            param_hash['well_id__eq'] = well_id
            if not library:
                library = Well.objects.get(well_id=well_id).library
        if library:
            library_classification = library.classification

        if library_classification is None:
            logger.warn('building response with no library classification: %r',
                {k:v for k,v in param_hash.items() if k != 'schema' })
        cherry_pick_request_id_screener = \
            param_hash.pop('cherry_pick_request_id_screener', None)
        if cherry_pick_request_id_screener is not None:
            cpr = CherryPickRequest.objects.get(
                cherry_pick_request_id=cherry_pick_request_id_screener)
            # NOTE: breaks for viewing natural product reagents
            library_classification = cpr.screen.screen_type
        cherry_pick_request_id_lab = \
            param_hash.pop('cherry_pick_request_id_lab', None)
        if cherry_pick_request_id_lab is not None:
            cpr = CherryPickRequest.objects.get(
                cherry_pick_request_id=cherry_pick_request_id_lab)
            # NOTE: breaks for viewing natural product reagents
            library_classification = cpr.screen.screen_type
        
        substance_id = param_hash.pop('substance_id', None)
        if substance_id:
            param_hash['substance_id__eq'] = substance_id
            if not library:
                library = Reagent.objects.get(
                    substance_id=substance_id).well.library

        # Note: build schema for each request to use the subtype
        schema = self.build_schema(library_classification=library_classification)

        manual_field_includes = set(param_hash.get('includes', []))
        content_type = self.get_serializer().get_accept_content_type(
            request, format=kwargs.get('format', None))
        if content_type == SDF_MIMETYPE:
            manual_field_includes.add('molfile')
            param_hash['includes'] = manual_field_includes
            
        (field_hash, columns, stmt, count_stmt, filename) = \
            self.get_query(
                param_hash, request.user,
                library_classification=library_classification,
                library=library,
                cherry_pick_request_id_screener=cherry_pick_request_id_screener,
                cherry_pick_request_id_lab=cherry_pick_request_id_lab,
                wells = wells)
        
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            # use "use_vocab" as a proxy to also adjust siunits for viewing
            rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                field_hash, rowproxy_generator)
        rowproxy_generator = \
            self._meta.authorization.get_row_property_generator(
                request.user, field_hash, rowproxy_generator)

        title_function = None
        if use_titles is True:
            def title_function(key):
                return field_hash[key]['title']
        if is_data_interchange:
            title_function = DbApiResource.datainterchange_title_function(
                field_hash,schema['id_attribute'])

        return self.stream_response_from_statement(
            request, stmt, count_stmt, filename,
            field_hash=field_hash, param_hash=param_hash,
            is_for_detail=is_for_detail,
            rowproxy_generator=rowproxy_generator,
            title_function=title_function, meta=kwargs.get('meta', None))
    
    def delete_reagents_for_library(self, library):
        self.get_reagent_resource(library.classification).delete_reagents(library)

    # NOTE: removed; patching will only be done in batch, from library/well
    # @transaction.atomic()    
    # def patch_obj(self, request, deserialized, **kwargs):
    #     well = kwargs.get('well', None)
    #     if not well:
    #         raise ValidationError(key='well', msg='required')
    #     
    #     reagent_resource = self.get_reagent_resource(
    #         well.library.classification)
    #     reagent = reagent_resource.patch_obj(request, deserialized, **kwargs)
    #     
    #     reagent.well = well
    #     reagent.save()
    #     return well

class WellResource(DbApiResource):

    class Meta:

        queryset = Well.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             IccblSessionAuthentication())
        resource_name = 'well'
        authorization = LibraryResourceAuthorization(resource_name)
        ordering = []
        filtering = {}
        serializer = LimsSerializer()   
        always_return_data = True 
        max_limit = 10000

    def __init__(self, **kwargs):
        self.library_resource = None
        self.reagent_resource = None
        super(WellResource, self).__init__(**kwargs)

    def get_reagent_resource(self):
        if not self.reagent_resource:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource

    def get_library_resource(self):
        if not self.library_resource:
            self.library_resource = LibraryResource()
        return self.library_resource

    def set_caching(self,use_cache):
        super(WellResource, self).set_caching(use_cache)
        self.get_reagent_resource().set_caching(use_cache)


    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})"
                r"/annotations%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_well_annotations_view'), 
                    name="api_dispatch_well_annotations_view"),
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})"
                 r"/duplex_wells%s$" )
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_well_duplex_view'), 
                    name="api_dispatch_well_duplex_view"),
        ]

    @read_authorization
    def dispatch_well_annotations_view(self, request, **kwargs):
        '''
        All Study annotations for a well
        Format is:
        [ { study 1 information, dc_1, dc_2, etc. }, { study 2 information ... }...]
        - where each study conforms to the study schema
        - each dc has:
            { dc_schema, value }
        '''
        well_id = kwargs.get('well_id', None)
        if not well_id:
            raise NotImplementedError('must provide a well_id parameter')
        
        screens = {}        
        for rv in (ResultValue.objects
                .filter(well_id=well_id)
                .filter(data_column__screen_result__screen__study_type__isnull=False)):
            screen = rv.data_column.screen_result.screen
            if screen.facility_id not in screens:
                _screen_data = screens.setdefault(
                    screen.facility_id,
                    { 
                        '1-screen_facility_id': screen.facility_id,
                        'facility_id': screen.facility_id,
                        'title': screen.title,
                        'summary': screen.summary,
                        'lead_screener_id': screen.lead_screener.screensaver_user_id,
                        'lead_screener_name': '%s %s' % (
                            screen.lead_screener.first_name, screen.lead_screener.last_name),
                        'lab_head_id': screen.lab_head.screensaver_user_id,
                        'lab_name': '%s %s' % (
                            screen.lab_head.first_name, screen.lab_head.last_name),
                        'date_created': screen.date_created,
                        'study_type': screen.study_type,
                        'values': {},
                        'fields': {}
                    })
            _screen_data = screens[screen.facility_id]
            (column_name, _dict) = \
                DataColumnResource._create_datacolumn_from_orm(rv.data_column)
            _screen_data['values'][column_name] = rv.value
            if _dict['data_type'] in ['decimal','integer','numeric']:
                _screen_data['values'][column_name] = rv.numeric_value
            _dict['visibility'] = ['l','d']
            _dict['title'] = _dict['name']
            _screen_data['fields'][column_name] = _dict
        
        final_data = sorted(screens.values(), key=lambda x: x['facility_id'])
        content_type = self.get_serializer().get_accept_content_type(
            request,format=kwargs.get('format', None))
        return HttpResponse(
            content=self.get_serializer().serialize(
                final_data, content_type),
            content_type=content_type)
        
    @read_authorization
    def dispatch_well_duplex_view(self, request, **kwargs):
        
        well_id = kwargs.get('well_id', None)
        if not well_id:
            raise NotImplementedError('must provide a well_id parameter')
        
        sr = SilencingReagent.objects.get(well__well_id=well_id)
        
        _aw = self.bridge['assay_well']
        _sr = self.bridge['screen_result']
        _s = self.bridge['screen']
        srr = ScreenResultResource()
        
        stmt = (
            select([_s.c.facility_id, _aw.c.confirmed_positive_value, 
                    _aw.c.is_positive])
                .select_from(_aw
                    .join(_sr, _aw.c.screen_result_id==_sr.c.screen_result_id)
                    .join(_s, _sr.c.screen_id==_s.c.screen_id))
                .where(_s.c.study_type != None))
        fields = ['facility_id', 'confirmed_positive_value', 'is_positive']

        with get_engine().connect()as conn:
            
            data_per_screen = {}
            well_info = {}
            for dw in sr.duplex_wells.all():
                
                item = { 'well_id': dw.well_id }
                sr = SilencingReagent.objects.get(well__well_id=dw.well_id)
                
                item['sequence'] = sr.sequence
                item['vendor_id'] = '%s %s' %(sr.vendor_name,sr.vendor_identifier)
                well_info[dw.well_id] = item
                
                result = conn.execute(stmt.where(_aw.c.well_id==dw.well_id))
                if logger.isEnabledFor(logging.DEBUG):
                    compiled_stmt = str(stmt.where(_aw.c.well_id==dw.well_id).compile(
                        dialect=postgresql.dialect(),
                        compile_kwargs={"literal_binds": True}))
                    logger.info('stmt: %r', compiled_stmt)
                for x in cursor_generator(result, fields):
                    screen_row = data_per_screen.setdefault(
                        x['facility_id'],
                        { 'screen_facility_id': x['facility_id']})
                    screen_row[dw.well_id] = x['confirmed_positive_value']
                            
            data = { 'duplex_wells': well_info,
                     'confirmed_positive_values': data_per_screen.values() }
                
        content_type = self.get_serializer().get_accept_content_type(
            request,format=kwargs.get('format', None))
        return HttpResponse(
            content=self.get_serializer().serialize(
                data, content_type),
            content_type=content_type)
        
    def get_schema(self, request, **kwargs):
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.debug('param hash: %r', param_hash)
        
        if not 'library_short_name' in param_hash:
            return self.build_response(request, 
                self.build_schema(user=request.user,**kwargs),
                **param_hash)
        
        library_short_name = kwargs.pop('library_short_name')
        try:
            library = Library.objects.get(short_name=library_short_name)
            return self.build_response(
                request, 
                self.build_schema(
                    user=request.user,
                    library_classification=library.classification),
                **kwargs)
            
        except Library.DoesNotExist, e:
            raise Http404(
                'Can not build schema - library def needed'
                'no library found for short_name: %r' % library_short_name)
                
    def build_schema(self, library_classification=None, user=None, **kwargs):

        data = super(WellResource, self).build_schema(user=user)
        if library_classification:
            sub_data = self.get_reagent_resource().build_schema(
                library_classification=library_classification, user=user)
            newfields = {}
            newfields.update(sub_data['fields'])
            newfields.update(data['fields'])
            data['fields'] = newfields
            
            data['content_types'] = sub_data['content_types']
        elif 'search' in kwargs:
            # Build the full schema for search
            # FIXME could determine the schema as in ReagentResource.build_list_response
            sub_data = self.get_reagent_resource().build_schema(
                library_classification='small_molecule', user=user)
            newfields = {}
            newfields.update(sub_data['fields'])
            sub_data = self.get_reagent_resource().build_schema(
                library_classification='rnai', user=user)
            newfields.update(sub_data['fields'])
            # all sub-fields are set not visible
            for field in newfields:
                if 'l' in field['visibility']:
                    field['visibility'].remove('l')
                if 'd' in field['visibility']:
                    field['visibility'].remove('d')
            newfields.update(data['fields'])
            data['fields'] = newfields
            # data['content_types'] = sub_data['content_types']
        else:
            pass
        return data
    
    @read_authorization
    def get_detail(self, request, **kwargs):

        well_id = kwargs.get('well_id', None)
        if not well_id:
            raise NotImplementedError('must provide a well_id parameter')

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.get_list(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):
        return self.get_reagent_resource().get_list(request, **kwargs)

    @write_authorization
    @un_cache        
    @transaction.atomic
    def put_list(self, request, **kwargs):
        '''
        Put list will replace all the library reagents
        '''
        if 'library_short_name' not in kwargs:
            raise BadRequest('library_short_name is required')
        library = Library.objects.get(
            short_name=kwargs['library_short_name'])
        logger.info('put wells for library: %r, deleting reagents...', library)
        self.get_reagent_resource().delete_reagents_for_library(library=library)

        logger.info('resetting well data...')
        well_query = Well.objects.all().filter(library=library)
        well_query.update(library_well_type='undefined')
        well_query.update(facility_id=None)
        well_query.update(is_deprecated=False)
        well_query.update(molar_concentration=None)
        well_query.update(mg_ml_concentration=None)
        
        logger.info('library %s reagents removed, patching...', library)
        return self.patch_list(request, **kwargs)
         
    def post_list(self, request, **kwargs):
        return self.patch_list(request, **kwargs)
    
    @write_authorization
    @un_cache        
    @transaction.atomic
    def patch_list(self, request, **kwargs):
        # TODO: 20160219 workflow for well/reagent patching: 
        # SS1:
        # if a reagent changes, then create a new reagent entry and update the 
        # "latest_released_reagent"
        # SS2:
        # version 1: Ok to update reagent in place and use well id to find it, 
        # because only one reagent is associated with a well
        # FUTURE: use update of reagent by well id as special case, only allow if 
        # only one well is attached to the reagent
        # Future: reagents are created first, then well is assoc with reagent
        # so put is allowed to create, but patch must reference extant reagent 
        # (changes reagent assoc with well)
         
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        if 'library_short_name' not in kwargs:
            raise BadRequest('library_short_name is required')
        library = Library.objects.get(
            short_name=kwargs['library_short_name'])
        logger.info(
            'patch_list: WellResource: library: %r...', library.short_name)
         
        deserialized = self.deserialize(request)
        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
 
        id_attribute = schema['id_attribute']
        kwargs_for_log = kwargs.copy()
        for id_field in id_attribute:
            ids = set()
            # Test for each id key; it's ok on create for ids to be None
            # TODO: this could be optimized to make the query smaller
            for _dict in [x for x in deserialized if x.get(id_field, None)]:
                ids.add(_dict.get(id_field))
            if ids:
                kwargs_for_log['%s__in' % id_field] = \
                    LIST_DELIMITER_URL_PARAM.join(ids)

        logger.info('get original state, for logging')
        kwargs_for_log['includes'] = ['*', '-molfile','-structure_image']
        # NOTE: do not consider "undefined" wells for diff logs (create actions 
        # will not be logged.
        kwargs_for_log['library_well_type__ne'] = 'undefined'
        original_data = self._get_list_response_internal(**kwargs_for_log)

        prev_version = library.version_number
        if library.version_number:
            library.version_number += 1
        else:
            library.version_number = 1
        library.save()
         
        library_log = self.make_log(request, **kwargs)
        library_log.diff_keys = ['version_number']
        library_log.diffs = { 'version_number': [prev_version, library.version_number]}
        library_log.ref_resource_name = 'library'
        # library_log.uri = self.get_library_resource().get_resource_uri(
        #     model_to_dict(library))
        library_log.key = library.short_name
        library_log.uri = '/'.join([library_log.ref_resource_name, library_log.key])
        # library_log.key = (
        #     '/'.join(
        #         self.get_library_resource()
        #             .get_id(model_to_dict(library)).values()))
        kwargs.update({ 'parent_log': library_log })
 
        logger.info('Cache library wells for patch...') 
        well_map = dict((well.well_id, well) 
            for well in library.well_set.all())
        if len(well_map) == 0:
            raise BadRequest('Library wells have not been created')
         
        logger.info('patch wells, count: %d', len(deserialized))
        # Note: wells can only be created on library creation
        fields = { key:field for key,field in schema['fields'].items()
            if key in schema['update_fields'] }
        for i,well_data in enumerate(deserialized):
            well_data['library_short_name'] = kwargs['library_short_name']
             
            well_id = well_data.get('well_id', None)
            if not well_id:
                well_name = well_data.get('well_name', None)
                plate_number = well_data.get('plate_number', None)
                if well_name and plate_number:                
                    well_id = '%s:%s' % (
                        str(plate_number).zfill(5), well_name)
                    well_data['well_id'] = well_id
 
            if not well_id:
                raise ValidationError(
                    key='well_id',
                    msg='well_id is required')
            
            well = well_map.get(well_id, None)
            if not well:
                raise ValidationError(
                    key='well_id',
                    msg=('well %r not found for this library %r'
                        % (well_id, well_data['library_short_name']))
                )
            well_data['well_id'] = well_id
            well_data['well'] = well

            initializer_dict = self.parse(well_data, fields=fields)
            logger.debug('well: %r: %r', well_id, initializer_dict)
            for key, val in initializer_dict.items():
                if hasattr(well, key):
                    setattr(well, key, val)
            
            well.save()
        
            duplex_wells = []
            if well_data.get('duplex_wells', None):
                if not library.is_pool:
                    raise ValidationError(
                        key='duplex_wells',
                        msg='library is not a pool libary: %r' % library.short_name)
                
                # TODO: batch get/set all of the duplex wells as well
                well_ids = well_data['duplex_wells']
                for well_id in well_ids:
                    try:
                        duplex_wells.append(Well.objects.get(well_id=well_id))
                    except:
                        raise ValidationError(
                            key='duplex_well not found',
                            msg='well: %r, pool well: %r' % (well.well_id, well_id))
                well_data['duplex_wells'] = duplex_wells
            if i and i % 1000 == 0:
                logger.info('patched: %d wells', i+1)
        logger.info('patched %d wells', i+1)

        self.get_reagent_resource().\
            get_reagent_resource(library.classification)._patch_wells(request, deserialized)        
        
        library.save()
        
        logger.info(
            'put_list: WellResource: library: %r; patch completed: %d', 
            library.short_name, i+1)
             
        self.get_reagent_resource().get_debug_times()
         
        experimental_well_count = library.well_set.filter(
            library_well_type__iexact='experimental').count()
        if library.experimental_well_count != experimental_well_count:
            library_log.diff_keys.append('experimental_well_count')
            library_log.diffs['experimental_well_count'] = \
                [library.experimental_well_count, experimental_well_count]
            library.experimental_well_count = experimental_well_count
            library.save()
 
        library_log.save()
                 
        logger.info('get new wells state, for logging...')
        new_data = self._get_list_response_internal(**kwargs_for_log)
         
        original_data_patches_only = []
        new_data_patches_only = []
        for item in original_data:
            for new_item in new_data:
                if item['well_id'] == new_item['well_id']:
                    original_data_patches_only.append(item)
                    new_data_patches_only.append(new_item)
         
        logger.debug('new data: %s', new_data_patches_only)
        logger.info('patch list done, original_data: %d, new data: %d' 
            % (len(original_data_patches_only), len(new_data_patches_only)))
        logs = self.log_patches(
            request, original_data_patches_only, new_data_patches_only,
            **kwargs)

        patch_count = len(deserialized)
        # Update: for wells, only measure what has diffed
        update_count = len([x for x in logs if x.diffs ])
        # Create: measure what is reported created (new_data), subtract updates,
        # because create actions are not included in updates for well patching.
        create_count = 0
        if not original_data:
            create_count = len(new_data) - update_count
        
        # FIXME:
        # update screen_experimental_well_count on library_screenings and on screens
        # see LibrariesDAO in SS1
        if update_count:
            self.update_screening_stats(library)
        
        meta = { 
            API_MSG_RESULT: { 
                API_MSG_SUBMIT_COUNT: patch_count, 
                API_MSG_UPDATED: update_count, 
                API_MSG_CREATED: create_count, 
                API_MSG_UNCHANGED: patch_count-update_count-create_count,
                API_MSG_ACTION: library_log.api_action, 
                API_MSG_COMMENTS: library_log.comment
            }
        }
        if not self._meta.always_return_data:
            return self.build_response(
                request, {API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)
        else:
            return self.build_response(
                request,  {API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)

    def update_screening_stats(self, library):

        with get_engine().connect() as conn:
            sql = (
                'with count_comparison as ( '
                '    with screenings as ( '
                '        select distinct(library_screening_id)  '
                '        from assay_plate ap, library  '
                '        where ap.plate_number between start_plate and end_plate   '
                '        and library_id = %s ) '
                '    select '
                '    activity_id, '
                '    screened_experimental_well_count, '
                '    count(distinct(well_id)) current_count '
                '    from assay_plate ap '
                '    join library_screening ls on(activity_id=library_screening_id) '
                '    join screenings on(ls.activity_id=screenings.library_screening_id) '
                '    join well on(well.plate_number=ap.plate_number) '
                "    where well.library_well_type='experimental' "
                '    group by ls.activity_id, screened_experimental_well_count  '
                ') '
                'update library_screening '
                'set screened_experimental_well_count = current_count '
                'from count_comparison  '
                'where count_comparison.activity_id=library_screening.activity_id '
                'and library_screening.screened_experimental_well_count != current_count; '
                )
            conn.execute(sql, (library.library_id))
            sql = (
                'with screen_update as ( '
                '    with screens as ( '
                '            select distinct(screen_id)  '
                '            from library_screening ls  '
                '            join assay_plate on(activity_id=library_screening_id) '
                '            join plate using(plate_id) '
                '            join copy using(copy_id) '
                '            where copy.library_id = %s '
                '    ) '
                '    select '
                '    screen.screen_id, '
                '    screen.facility_id, '
                '    screen.screened_experimental_well_count, '
                '    screen.unique_screened_experimental_well_count, '
                '    count(well_id) as current_count, '
                '    count(distinct(well_id)) as current_unique_count '
                '    from well w,  '
                '    screen join assay_plate ap using(screen_id) '
                '    join plate p using(plate_id)  '
                '    join screens on(screens.screen_id=screen.screen_id) '
                "    where w.library_well_type = 'experimental' "
                '    and ap.replicate_ordinal = 0 '
                '    and w.plate_number = p.plate_number '
                '    group by '
                '    screen.screen_id, '
                '    screen.facility_id, '
                '    screen.screened_experimental_well_count, '
                '    screen.unique_screened_experimental_well_count ) '
                'update screen '
                'set screened_experimental_well_count = current_count, '
                '  unique_screened_experimental_well_count = current_unique_count '
                'from screen_update '
                'where screen_update.screen_id=screen.screen_id '
                'and ( screen_update.screened_experimental_well_count!= current_count or  '
                '      screen_update.unique_screened_experimental_well_count != current_unique_count ); '
                )
            conn.execute(sql, (library.library_id))


    def patch_obj(self, request, deserialized, **kwargs):
        
        raise NotImplementedError('patch obj must be implemented')
        
    @classmethod
    def find_wells(cls, well_search_data ):
        ''' return set() of Well objects matching the line based 
        well_search_data entered by the user
        '''
        
        logger.info('find wells for patterns: %r', well_search_data)
        if not isinstance(well_search_data, (list,tuple)):
            well_search_data = (well_search_data,)

        wells = set()
        errors = []
        
        # Process the patterns by line
        parsed_lines = well_search_data
        if isinstance(parsed_lines, basestring):
            parsed_lines = re.split(r'\n+',parsed_lines)
            logger.info('parsed_lines: %r', parsed_lines)
        for _line in parsed_lines:
            _line = _line.strip()
            if not _line:
                continue
            logger.info('find_wells: well pattern line: %r', _line)
            
            line_wells = set()
            plate_numbers = set()
            patterns = re.split(r'[\s,]+', _line)
            search_kwargs = defaultdict(set)
            for _pattern in patterns:
                if not _pattern:
                    continue
                # get rid of quotes
                _pattern = re.sub(r'["\']+','',_pattern)
                match = WELL_ID_PATTERN.match(_pattern)
                if match is not None:
                    plate_number = int(match.group(1))
                    search_kwargs['plate_number'].add(plate_number)
                    
                    plate_numbers.add(plate_number)
                    well_name = match.group(2).upper()
                    search_kwargs['well_name'].add(well_name)
                elif PLATE_PATTERN.match(_pattern):
                    logger.info('plate pattern match: verify plates exist: %r', _pattern)
                    plate_query = Plate.objects.all().filter(plate_number=_pattern)
                    if plate_query.exists():
                        logger.info('found %r for %r', plate_numbers, _pattern)
                        search_kwargs['plate_number'].add(int(_pattern))
                        continue
                    else:
                        errors.append(
                            'plate: %r is does not exist' % _pattern) 
                elif PLATE_RANGE_PATTERN.match(_pattern):
#                     errors.append('pattern: %r, plate range pattern not supported' % _pattern)
                    match = PLATE_RANGE_PATTERN.match(_pattern)
                    plate_query = Plate.objects.all().filter(plate_number__range=sorted([
                        int(match.group(1)), int(match.group(2))]))
                    search_kwargs['plate_number'].update([p.plate_number for p in plate_query ])
                elif WELL_NAME_PATTERN.match(_pattern):
                    search_kwargs['well_name'].add(_pattern)
                else:
                    errors.append(
                        'pattern: %r is not a recognized as a well id, or '
                        'a plate number followed by a well name' % _pattern)
            if search_kwargs:
                query = Well.objects.all()
                plate_numbers = search_kwargs.get('plate_number', None)
                if not plate_numbers:
                    raise Exception('no plate numbers found: %r for %r', search_kwargs, _line)
                plate_numbers = [p for p in plate_numbers]
                if len(plate_numbers) > 1:
                    query = query.filter(plate_number__in=plate_numbers)
                else:
                    query = query.filter(plate_number=plate_numbers[0])
                well_names = search_kwargs.get('well_name', None)
                if well_names:
                    well_names = [w for w in well_names]
                    if len(well_names) > 1:
                        query = query.filter(well_name__in=well_names)
                    else: 
                        query = query.filter(well_name=well_names[0])
                line_wells = [w for w in query]
            logger.info('line: %r, wells: %d', _line, len(line_wells))
            if not line_wells:
                errors.append('no matches found for line: %r', _line)
            
            wells.update(line_wells)
        return (wells, errors)
        
class LibraryResource(DbApiResource):
    
    class Meta:

        queryset = Library.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(
            BasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'library'
        authorization = LibraryResourceAuthorization(resource_name)
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 
        
    def __init__(self, **kwargs):
        
        self.well_resource = None
        self.apilog_resource = None
        self.reagent_resource = None
        self.screen_resource = None
        super(LibraryResource, self).__init__(**kwargs)
        
    def get_apilog_resource(self):
        if not self.apilog_resource:
            self.apilog_resource = ApiLogResource()
        return self.apilog_resource
    
    def get_well_resource(self):
        if not self.well_resource:
            self.well_resource = WellResource()
        return self.well_resource

    def get_reagent_resource(self):
        if not self.reagent_resource:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource
    
    def get_screen_resource(self):
        if not self.screen_resource:
            self.screen_resource = ScreenResource()
        return self.screen_resource

    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/plate/(?P<plate_number>[^/]+)%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copyplateview'),
                name="api_dispatch_library_copyplateview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/plate%s$") % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copyplateview'),
                name="api_dispatch_library_copyplateview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/copywell/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copywellview'),
                name="api_dispatch_library_copywellview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/copywellhistory/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copywellhistoryview'),
                name="api_dispatch_library_copywellhistoryview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/copywell%s$") % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copywellview'),
                name="api_dispatch_library_copywellview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)%s$") 
                 % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copyview'),
                name="api_dispatch_library_copyview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/copy%s$") % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copyview'),
                name="api_dispatch_library_copyview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/plate%s$") % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copyplateview'),
                name="api_dispatch_library_copyplateview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/plate/(?P<plate_number>[^/]+)%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_copyplateview'),
                name="api_dispatch_library_copyplateview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/well%s$") % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_wellview'),
                name="api_dispatch_library_wellview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/reagent%s$") % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_library_reagentview'),
                name="api_dispatch_library_reagentview"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/reagent/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_reagent_schema'),
                name="api_get_reagent_schema"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/well/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_well_schema'),
                name="api_get_well_schema"),
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w.\-\+: ]+)"
                 r"/version%s$") % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_libraryversionview'),
                name="api_dispatch_libraryversionview"),
        ]    
    
    def get_well_schema(self, request, **kwargs):
        if not 'short_name' in kwargs:
            raise Http404(
                'The well schema requires a library short name'
                ' in the URI, as in /library/[short_name]/well/schema/')
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_well_resource().get_schema(request, **kwargs)    
  
    def get_reagent_schema(self, request, **kwargs):
        if not 'short_name' in kwargs:
            raise Http404(
                'The reagent schema requires a library short name'
                ' in the URI, as in /library/[short_name]/well/schema/')
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_reagent_resource().get_schema(request, **kwargs)    
  
    def dispatch_library_copyview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyResource().dispatch('list', request, **kwargs)    
 
    def dispatch_library_copyplateview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)   

    def dispatch_library_copywellview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return CopyWellResource().dispatch('list', request, **kwargs)   

    def dispatch_library_copywellhistoryview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return CopyWellHistoryResource().dispatch('list', request, **kwargs)   

    def dispatch_library_wellview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_well_resource().dispatch('list', request, **kwargs)    
                    
    def dispatch_library_reagentview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_reagent_resource().dispatch('list', request, **kwargs)    

    @read_authorization
    def get_detail(self, request, **kwargs):

        library_short_name = kwargs.pop('short_name', None)
        if not library_short_name:
            raise NotImplementedError('must provide a short_name parameter')
        else:
            kwargs['short_name__eq'] = library_short_name

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, schema=None, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        if schema is None:
            raise Exception('schema not initialized')
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        for_screen_facility_id = param_hash.pop('for_screen_facility_id', None)
        if for_screen_facility_id is not None:
            if self.get_screen_resource()._meta.authorization\
                .has_screen_read_authorization(
                    request.user,for_screen_facility_id) is False:
                raise PermissionDenied
        
        try:
            # general setup
            
            manual_field_includes = set(param_hash.get('includes', []))
            manual_field_includes.add('comment_array')
            if is_for_detail:
                manual_field_includes.add('concentration_types')
 
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                
            order_params = param_hash.get('order_by', [])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes,
                param_hash.get('visibilities'),
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
            
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)

            # specific setup
            _l = self.bridge['library']
            
#             _apilog = self.bridge['reports_apilog']
#             _logdiff = self.bridge['reports_logdiff']
#             _comment_apilogs = (
#                 select([
#                     _apilog.c.date_time, 
#                     _apilog.c.key,
#                     _apilog.c.username, 
#                     _apilog.c.comment])
#                 .select_from(_apilog)
#                 .where(_apilog.c.ref_resource_name==self._meta.resource_name)
#                 .where(_apilog.c.comment!=None)
#                 .order_by(desc(_apilog.c.date_time)))
            _comment_apilogs = ApiLogResource.get_resource_comment_subquery(
                self._meta.resource_name)
            _comment_apilogs = _comment_apilogs.cte('_comment_apilogs')

            custom_columns = {
                'comment_array': (
                    select([func.array_to_string(
                        func.array_agg(
                            _concat(                            
                                cast(_comment_apilogs.c.name,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                cast(_comment_apilogs.c.date_time,
                                    sqlalchemy.sql.sqltypes.Text),
                                LIST_DELIMITER_SUB_ARRAY,
                                _comment_apilogs.c.comment)
                        ), 
                        LIST_DELIMITER_SQL_ARRAY) ])
                    .select_from(_comment_apilogs)
                    .where(_comment_apilogs.c.key==_l.c.short_name)),
#                 'comments': (
#                     select([
#                         func.array_to_string(func.array_agg(
#                             literal_column('comment')), LIST_DELIMITER_SQL_ARRAY)])
#                         .select_from(_apilog)
#                         .where(_apilog.c.ref_resource_name=='library')
#                         .where(_apilog.c.key==literal_column('library.short_name'))
#                         # .where(not_(exists(
#                         #     select([None]).select_from(_logdiff)
#                         #     .where(_logdiff.c.log_id==_apilog.c.id))))
#                         .where(_apilog.c.comment!=None)
#                         ),
                'copy_plate_count': literal_column(
                    '(select count(distinct(p.plate_id))'
                    '    from plate p join copy c using(copy_id)'
                    '    where c.library_id=library.library_id)'
                    ).label('plate_count'),
                'plate_count': literal_column(
                    '(select count(distinct(w.plate_number))'
                    '    from well w'
                    '    where w.library_id=library.library_id)'
                    ).label('plate_count'),
                'copies': literal_column(
                    "(select array_to_string(array_agg(c1.name),'%s') "
                    '    from ( select c.name from copy c '
                    '    where c.library_id=library.library_id '
                    '    order by c.name) as c1 )' % LIST_DELIMITER_SQL_ARRAY
                    ).label('copies'),
                'screening_copies': literal_column(
                    "(select array_to_string(array_agg(c1.name),'%s') "
                    '    from ( select distinct(c.name) from copy c '
                    '    join plate p using(copy_id) '
                    '    where c.library_id=library.library_id '
                    "    and c.usage_type='library_screening_plates' "
                    "    and p.status in ('available') "
                    '    order by c.name) as c1 )' % LIST_DELIMITER_SQL_ARRAY
                    ).label('copies'),
                # TODO: copies2 is the same in all respects, except that it is 
                # used differently in the UI - not displayed as a list of links
                'copies2': literal_column(
                    "(select array_to_string(array_agg(c1.name),'%s') "
                    '    from ( select c.name from copy c '
                    '    where c.library_id=library.library_id '
                    '    order by c.name) as c1 )' % LIST_DELIMITER_SQL_ARRAY
                    ).label('copies2'),
                'owner': literal_column(
                    "(select u.first_name || ' ' || u.last_name "
                    '    from screensaver_user u '
                    '    where u.screensaver_user_id=library.owner_screener_id)'
                    ).label('owner'),
                'concentration_types': literal_column(
                    '(select array_to_string(ARRAY['
                    ' case when exists(select null from well '
                    '  where well.library_id=library.library_id '
                    "  and well.mg_ml_concentration is not null limit 1) then 'mg_ml' end, "
                    ' case when exists(select null from well '
                    '  where well.library_id=library.library_id '
                    "  and well.molar_concentration is not null limit 1) then 'molar' end "
                    "], '%s' ) )" % LIST_DELIMITER_SQL_ARRAY
                    )
                }
            ##### FIXME: 20171106 unused
#             concentration_sql = (
#                 'select short_name, '
#                 'exists(select null from well '
#                 '  where well.library_id=l.library_id '
#                 '  and well.mg_ml_concentration is not null limit 1) mg_ml_measure, '
#                 'exists(select null from well '
#                 '  where well.library_id=l.library_id '
#                 '  and well.molar_concentration is not null limit 1) molar_measure '
#                 'from library l;')
#             with get_engine().connect() as conn:
#                 result = conn.execute(text(concentration_sql))
#                 cached_ids =  [x[0] for x in result ]
            
                     
            base_query_tables = ['library']

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            # build the query statement

            j = _l
            stmt = select(columns.values()).select_from(j)

            if for_screen_facility_id:
                stmt = stmt.where(_l.c.library_id.in_(
                    self.get_screen_library_ids(for_screen_facility_id)))

            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            if not order_clauses:
                stmt = stmt.order_by("short_name")

            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            if is_data_interchange:
                title_function = DbApiResource.datainterchange_title_function(
                    field_hash,schema['id_attribute'])
                
            if False:
                logger.info(
                    'stmt: %s',
                    str(stmt.compile(
                        dialect=postgresql.dialect(),
                        compile_kwargs={"literal_binds": True})))
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash,
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True)
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

    @classmethod
    def get_screen_library_ids(cls, for_screen_facility_id):
        
        bridge = get_tables()
        _screen = bridge['screen']
        _lab_activity = bridge['lab_activity']
        _library_screening = bridge['library_screening']
        _assay_plate = bridge['assay_plate']
        _plate = bridge['plate']
        _copy = bridge['copy']
        
        j = _screen
        j = j.join(
            _lab_activity,
            _lab_activity.c.screen_id == _screen.c.screen_id)
        j = j.join(
            _library_screening,
            _library_screening.c.activity_id
                == _lab_activity.c.activity_id)
        j = j.join(
            _assay_plate,
            _library_screening.c.activity_id
                == _assay_plate.c.library_screening_id)
        j = j.join(_plate, _assay_plate.c.plate_id == _plate.c.plate_id)
        j = j.join(_copy, _copy.c.copy_id == _plate.c.copy_id)
        with get_engine().connect() as conn:
            query = (
                select([
                    distinct(_copy.c.library_id).label('library_id')])
                .select_from(j)
                .where(_screen.c.facility_id == for_screen_facility_id))
            library_ids = [x[0] for x in 
                conn.execute(query)]
            return library_ids
       
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        Library.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        logger.info('patch library: %r', deserialized)
        id_kwargs = self.get_id(deserialized, validate=True, schema=schema, **kwargs)
        # create/update the library
        create = False
        try:
            library = None
            try:
                library = Library.objects.get(**id_kwargs)
            except ObjectDoesNotExist, e:
                create = True
                logger.info('Library %s does not exist, creating', id_kwargs)
                library = Library(**id_kwargs)

            initializer_dict = self.parse(deserialized, create=create, schema=schema)
            errors = self.validate(initializer_dict, schema=schema,patch=not create)
            if errors:
                raise ValidationError(errors)
            
            for key, val in initializer_dict.items():
                if hasattr(library, key):
                    setattr(library, key, val)
            
            library.save()

            # now create the wells
            if create is True:
                plate_size = int(library.plate_size)
        
                try:
                    i = 0
                    logger.info('bulk create wells: plates: %s-%s, plate_size: %d', 
                        library.start_plate, library.end_plate, plate_size)
                    for plate in range(
                            int(library.start_plate), int(library.end_plate) + 1):
                        bulk_create_wells = []
                        for index in range(0, plate_size):
                            well = Well()
                            well.well_name = lims_utils.well_name_from_index(
                                index, plate_size)
                            well.well_id = lims_utils.well_id(plate, well.well_name)
                            well.library = library
                            well.plate_number = plate
                            # FIXME: use vocabularies for well type
                            well.library_well_type = 'undefined'
                            bulk_create_wells.append(well)
                            i += 1
                            if i % 1000 == 0:
                                logger.info('created %d wells', i)
                        Well.objects.bulk_create(bulk_create_wells)
                    logger.info(
                        'created %d wells for library %r, %r',
                        i, library.short_name, library.library_id)
                except Exception, e:
                    logger.exception('on library wells create')
                    raise e

            logger.info('patch_obj done')

            # clear the cached schema because plate range have updated
            cache.delete(self._meta.resource_name + ':schema')
            
            return { API_RESULT_OBJ: library }
            
        except Exception, e:
            logger.exception('on patch detail')
            raise e  
    

class ResourceResource(reports.api.ResourceResource):
    '''
    Motivation: to extend the reports.ResourceResource with "db" specific
    extensions.
    '''
    
    def _build_resources(self, user=None):
        
        resources = None
        if self.use_cache:
            resources = cache.get('dbresources')
        if not resources:

            resources = super(ResourceResource, self)._build_resources(user=user, use_cache=False)
        
            for key,resource in resources.items():
                self.extend_resource_specific_data(resource)
                
            cache.set('dbresources', resources)
    
        return resources
    
    def extend_resource_specific_data(self, resource_data):
        try:
            key = resource_data['key']
            if key == 'library':
                ranges = (Library.objects.all()
                    .order_by('start_plate')
                    .values_list('start_plate', 'end_plate'))
                plate_ranges = []
                temp = 0
                for s, e in ranges:
                    if temp == 0:
                        plate_ranges.append(s)
                    if s > temp+1:
                        plate_ranges.append(temp)
                        plate_ranges.append(s)
                    temp = e
                plate_ranges.append(temp)
                resource_data['library_plate_ranges'] = plate_ranges
     
                temp = [ 
                    x.library_type 
                    for x in Library.objects.all().distinct('library_type')]
                resource_data['extraSelectorOptions'] = { 
                    'label': 'Type',
                    'searchColumn': 'library_type',
                    'options': temp }
                
            elif key == 'librarycopyplate':
                temp = [ x for x in 
                    Plate.objects.all().distinct('status')
                        .values_list('status', flat=True)]
                resource_data['extraSelectorOptions'] = { 
                    'label': 'Type', 'searchColumn': 'status', 'options': temp }
            
            elif key == 'screen':
                temp = [ x.screen_type 
                    for x in Screen.objects.all().distinct('screen_type')]
                resource_data['extraSelectorOptions'] = { 
                    'label': 'Type', 'searchColumn': 'screen_type', 'options': temp }
        except Exception, e:
            # TODO: remove after migration
            logger.exception('catch error on extend_resource_specific_data (migration)')
        
#         elif key == 'labcherrypick':
#             resource_data['extraSelectorOptions'] = {
#                 'label': 'Status',
#                 'searchColumn': 'status',
#                 'options': ['unfulfilled','selected','plated','not_selected']}
            
