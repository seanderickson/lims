from __future__ import unicode_literals

from copy import deepcopy
# import csv
import unicodecsv
import hashlib
import io
import json
import logging
import re
from tempfile import SpooledTemporaryFile, NamedTemporaryFile
import time
import urllib

from aldjemy.core import get_engine, get_tables
import aldjemy.core
from django.conf import settings
from django.conf.urls import url
from django.core.cache import cache
from django.core.cache import caches
from django.core.exceptions import ObjectDoesNotExist
from django.core.serializers.json import DjangoJSONEncoder
from django.db import connection
from django.db import transaction
from django.db.models import F
from django.forms.models import model_to_dict
from django.http import Http404
from django.http.request import HttpRequest
from django.http.response import StreamingHttpResponse, HttpResponse
from django.utils import timezone
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
from tastypie.authentication import BasicAuthentication, SessionAuthentication, \
    MultiAuthentication
from tastypie.exceptions import BadRequest
from tastypie.http import HttpNotFound
import tastypie.http
from tastypie.utils.urls import trailing_slash

from db.models import ScreensaverUser, Screen, \
    ScreenResult, DataColumn, Library, Plate, Copy, \
    CopyWell, UserFacilityUsageRole, \
    PlateLocation, Reagent, Well, Activity, \
    AdministrativeActivity, SmallMoleculeReagent, SilencingReagent, GeneSymbol, \
    NaturalProductReagent, Molfile, Gene, GeneGenbankAccessionNumber, \
    CherryPickRequest, CherryPickAssayPlate, CherryPickLiquidTransfer, \
    CachedQuery, UserChecklistItem, AttachedFile, \
    ServiceActivity, LabActivity, Screening, LibraryScreening, AssayPlate, \
    SmallMoleculeChembankId, SmallMoleculePubchemCid, SmallMoleculeChemblId, \
    SmallMoleculeCompoundName, ScreenCellLines, ScreenFundingSupports, \
    ScreenKeyword, ResultValue, AssayWell, Publication, ScreenerCherryPick,\
    LabCherryPick, CherryPickRequestEmptyWell
from db.support import lims_utils, screen_result_importer, bin_packer
from db.support.data_converter import default_converter
from db.support.screen_result_importer import PARTITION_POSITIVE_MAPPING, \
    CONFIRMED_POSITIVE_MAPPING
from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    HTTP_PARAM_USE_TITLES, HTTP_PARAM_USE_VOCAB, HTTP_PARAM_DATA_INTERCHANGE, \
    LIST_BRACKETS, HTTP_PARAM_RAW_LISTS, HEADER_APILOG_COMMENT, ValidationError
from reports import ValidationError, InformationError, _now
from reports.api import API_MSG_COMMENTS, API_MSG_CREATED, \
    API_MSG_SUBMIT_COUNT, API_MSG_UNCHANGED, API_MSG_UPDATED, \
    API_MSG_ACTION, API_MSG_RESULT, API_MSG_WARNING, API_RESULT_DATA, \
    API_RESULT_META, API_RESULT_OBJ, API_MSG_NOT_ALLOWED, API_PARAM_OVERRIDE
from reports.api import ApiLogResource, UserGroupAuthorization, \
    VocabularyResource, UserResource, UserGroupResource, ApiLogResource, \
    write_authorization, read_authorization
import reports.api
from reports.api_base import un_cache
from reports.models import Vocabulary, ApiLog, UserProfile, \
    API_ACTION_DELETE, API_ACTION_PUT, API_ACTION_PATCH, API_ACTION_CREATE
from reports.serialize import parse_val, XLSX_MIMETYPE, SDF_MIMETYPE, \
    XLS_MIMETYPE, JSON_MIMETYPE, CSV_MIMETYPE, ZIP_MIMETYPE
from reports.serialize.csvutils import convert_list_vals
# from reports.serialize import encode_utf8
from reports.serialize.streaming_serializers import ChunkIterWrapper, \
    json_generator, cursor_generator, sdf_generator, generic_xlsx_response, \
    csv_generator, get_xls_response, image_generator, closing_iterator_wrapper,\
    FileWrapper1
from reports.serializers import LimsSerializer, \
    XLSSerializer, ScreenResultSerializer
from reports.sqlalchemy_resource import SqlAlchemyResource
from reports.sqlalchemy_resource import _concat
from decimal import Decimal
import six
from db import WELL_ID_PATTERN, WELL_NAME_PATTERN
from tastypie.resources import convert_post_to_put
from itertools import chain, combinations
from docutils.parsers.rst.directives.html import Meta
from math import ceil
import random
from django.db.utils import ProgrammingError
from zipfile import ZipFile
from collections import defaultdict, OrderedDict
import cStringIO
import os

PLATE_NUMBER_SQL_FORMAT = 'FM9900000'
PSYCOPG_NULL = '\\N'
MAX_SPOOLFILE_SIZE = 100*1024

# TODO: API constants: move to an API-accessible properties file
API_MSG_COPYWELLS_DEALLOCATED = 'copywells deallocated'
API_MSG_SCREENING_PLATES_UPDATED = 'Library Plates updated'
API_MSG_SCREENING_EXTANT_PLATE_COUNT = 'Extant plates'
API_MSG_SCREENING_ADDED_PLATE_COUNT = 'Added plates'
API_MSG_SCREENING_DELETED_PLATE_COUNT = 'Deleted plates'
API_MSG_SCREENING_TOTAL_PLATE_COUNT = 'Library Plate Count'
API_MSG_SCP_CREATED = 'created'
API_MSG_SCP_UNSELECTED = 'unselected'
API_MSG_SCP_RESELECTED = 'reselected'
API_MSG_SCPS_DELETED = 'Screener Cherry Picks Removed'      
API_MSG_SCPS_CREATED = 'Screener Cherry Picks Created'
API_MSG_LCP_CHANGED = 'changed'
API_MSG_LCP_DESELECTED = 'deselected'
API_MSG_LCP_SELECTED = 'selected'
API_MSG_LCP_MULTIPLE_SELECTIONS_SUBMITTED = 'Multiple lab cherry pick selections submitted'
API_MSG_LCPS_CREATED = 'Lab Cherry Picks Created'
API_MSG_LCPS_MUST_BE_DELETED = "Lab Cherry Picks must be deleted"
API_MSG_LCPS_ASSIGNED = 'Assigned to Copies'
API_MSG_LCPS_UNFULFILLED = 'Unfulfilled'
API_MSG_LCPS_INSUFFICIENT_VOLUME = 'Insufficient volume'
API_MSG_LCPS_VOLUME_OVERRIDDEN = 'Insufficient volume overridden'
API_MSG_LCPS_REMOVED = 'Lab Cherry Picks Removed'
API_MSG_LCP_PLATES_ASSIGNED = 'Copy Plate assigned'
API_MSG_LCP_ASSAY_PLATES_PLATED = 'cherry_pick_assay_plates plated'
API_MSG_LCP_ASSAY_PLATES_CREATED = 'cherry_pick_assay_plates created'
API_MSG_CPR_ASSAY_PLATES_REMOVED = 'cherry_pick_assay_plates removed'
API_MSG_PLATING_CANCELED = 'plating canceled'
API_MSG_CPR_PLATES_PLATED = 'cherry pick assay plates plated'
API_MSG_CPR_PLATES_SCREENED = 'cherry pick assay plates screened'
API_MSG_CPR_PLATED_CANCEL_DISALLOWED = 'Plating reservation may not be canceled after plates are plated'
VOCAB_LCP_STATUS_SELECTED = 'selected'
VOCAB_LCP_STATUS_NOT_SELECTED = 'not_selected'
VOCAB_LCP_STATUS_UNFULFILLED = 'unfulfilled'
VOCAB_LCP_STATUS_PLATED = 'plated'

# API_PARAM_PLATE_MAPPING_OVERRIDE = 'plate_mapping_override'
API_PARAM_SHOW_OTHER_REAGENTS = 'show_other_reagents'
API_PARAM_SHOW_COPY_WELLS = 'show_copy_wells'
API_PARAM_SHOW_RETIRED_COPY_WELlS = 'show_available_and_retired_copy_wells'
API_PARAM_SHOW_UNFULFILLED = 'show_unfulfilled'
API_PARAM_VOLUME_OVERRIDE = 'volume_override'

logger = logging.getLogger(__name__)
    
def _get_raw_time_string():
  return timezone.now().strftime("%Y%m%d%H%M%S")

class DbApiResource(reports.api.ApiResource):
    '''
   '''
    def __init__(self, **kwargs):
        super(reports.api.ApiResource,self).__init__(**kwargs)
        self.resource_resource = None

    def get_resource_resource(self):
        if not self.resource_resource:
            self.resource_resource = ResourceResource()
        return self.resource_resource

    
    def build_schema(self, user=None):
        logger.debug('build schema for: %r', self._meta.resource_name)
        return self.get_resource_resource()._get_resource_schema(
            self._meta.resource_name, user=user)
    

class PlateLocationResource(DbApiResource):        
    
    class Meta:
        queryset = PlateLocation.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'platelocation'
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
            # schema = super(PlateLocationResource, self).build_schema()
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            manual_field_includes.add('plate_location_id')
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)
                 
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
                    raise Http404
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
 
            if not order_clauses:
                stmt = stmt.order_by("room","freezer","shelf","bin")

            if True: 
                logger.info(
                    'stmt: %s',
                    str(stmt.compile(
                        dialect=postgresql.dialect(),
                        compile_kwargs={"literal_binds": True})))


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

    @un_cache        
    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_detail must be implemented')

    @write_authorization
    @un_cache  
    @transaction.atomic      
    def patch_detail(self, request, **kwargs):
        '''
        Override to generate informational summary for callee
        '''
        
        deserialized = kwargs.pop('data', None)
        # allow for internal data to be passed
        if deserialized is None:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        logger.debug('patch detail %s, %s', deserialized,kwargs)

        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        schema = kwargs['schema']
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(deserialized,validate=True,**kwargs)

        original_data = None
        if kwargs_for_log:
            try:
                original_data = self._get_detail_response(request,**kwargs_for_log)
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
        new_data = self._get_detail_response(request,**kwargs_for_log)
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
    @un_cache
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):

        logger.info('patch platelocation')
        schema = kwargs['schema']
        id_kwargs = self.get_id(deserialized, validate=True, **kwargs)

        with transaction.atomic():
            create = False
            try:
                plate_location = PlateLocation.objects.get(**id_kwargs)
            except ObjectDoesNotExist:
                logger.info('plate location does not exist, creating: %r',
                    id_kwargs)
                create=True
                plate_location = PlateLocation.objects.create(**id_kwargs)

            initializer_dict = self.parse(deserialized, create=create)
            errors = self.validate(initializer_dict, patch=not create)
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
                    'plate_number': str(plate.plate_number).zfill(5)
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
                    'plate_number': str(plate.plate_number).zfill(5)
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
                self.get_librarycopyplate_resource().build_schema()['id_attribute']
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
            
    def _find_plates(self, regex_string, copy_plate_ranges):
        logger.info('find copy_plate_ranges: %r', copy_plate_ranges)
        # parse library_plate_ranges
#         schema = self.build_schema()
        
        # E.G. Regex: /(([^:]+):)?(\w+):(\d+)-(\d+)/
#         regex_string = schema['fields']['copy_plate_ranges']['regex']
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
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'librarycopyplate'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 

    def __init__(self, **kwargs):
        
        super(LibraryCopyPlateResource, self).__init__(**kwargs)
        self.plate_location_resource = None
        
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
                func.to_char(_p.c.plate_number,PLATE_NUMBER_SQL_FORMAT).label('plate_number'),
                _p.c.copy_id,
                _c.c.library_id,
                _c.c.name.label('copy_name'),
                _concat(_l.c.short_name, '/', _c.c.name, '/', 
                    func.to_char(_p.c.plate_number,PLATE_NUMBER_SQL_FORMAT))
                    .label('key'),
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
        # TODO: adding plate_ids is not performant here
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
                func.to_char(_ap.c.plate_number,PLATE_NUMBER_SQL_FORMAT).label('plate_number'),
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
        
        for_screen_id = param_hash.pop('for_screen_id', None)
        # NOTE: loaded plates not being shown anymore
        # # FIXME: copyplatesloaded no longer works - 20160607
        # # because we are not creating "assay_plates" for screen results anymore
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        library_id = None
        library_short_name = param_hash.pop('library_short_name',
            param_hash.get('library_short_name__eq', None))
        if not library_short_name:
            logger.info('no library_short_name provided')
        else:
            log_key = library_short_name
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
            plate_number = str(int(plate_number)).zfill(5)
            param_hash['plate_number__eq'] = plate_number
        
        log_key = '/'.join(x if x else '%' 
            for x in [library_short_name,copy_name,plate_number])
        if log_key == '%/%/%':
            log_key = None
        try:
            plate_ids = None
            if for_screen_id:
                with get_engine().connect() as conn:
                    plate_ids = [x[0] for x in 
                        conn.execute(
                            self.get_screen_librarycopyplate_subquery(for_screen_id))]
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))

            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)
            # if (filter_expression is None
            #         and for_screen_id is None and loaded_for_screen_id is None):
            #     raise InformationError(
            #         key='Filter Data',
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
 
            # specific setup 
 
            _apilog = self.bridge['reports_apilog']
            _diff = self.bridge['reports_logdiff']
            _p = self.bridge['plate']
            _pl = self.bridge['plate_location']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _ls = self.bridge['library_screening']
            _a = self.bridge['activity']
            _ap = self.bridge['assay_plate']
            _plate_cte = (
                self.get_librarycopyplate_cte(
                    filter_hash=filter_hash,
                    library_id=library_id, copy_id=copy_id, plate_ids=plate_ids)
                        .cte('plate_cte'))
            _user_cte = ScreensaverUserResource.get_user_cte().cte('user_cte')
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
            _status_apilogs = (
                select([
                    _apilog.c.date_time,
                    _apilog.c.key,
                    _user_cte.c.username,
                    _user_cte.c.name,                    
                    ])
                .select_from(
                    _apilog.join(status_updated_apilogs, _apilog.c.id==status_updated_apilogs.c.id)
                        .join(_user_cte, _apilog.c.username==_user_cte.c.username))
                )
            if log_key:
                _status_apilogs = _status_apilogs.where(_apilog.c.key.ilike(log_key))
            _status_apilogs = _status_apilogs.cte('_status_apilogs')
            
            
            # j = _apilog.join(_diff,_apilog.c.id==_diff.c.log_id)
            # j = j.join(_user_cte, _apilog.c.username==_user_cte.c.username)
            # status_apilogs = (
            #     select([
            #         _apilog.c.date_time,
            #         _apilog.c.key,
            #         _diff.c.before,
            #         _diff.c.after,
            #         _user_cte.c.username,
            #         _user_cte.c.name,
            #         ])
            #     .select_from(j)
            #     .where(_apilog.c.ref_resource_name==self._meta.resource_name)
            #     .where(_diff.c.field_key=='status')
            #     .cte('status_apilogs'))
            
            custom_columns = {
                'librarycopyplate_id': (
                    _concat(_l.c.short_name,'/',_c.c.name,'/',
                        literal_column("to_char(plate.plate_number,'%s')" 
                            % PLATE_NUMBER_SQL_FORMAT))),
                'copyplate_id': (
                    _concat(_c.c.name,'/',
                        literal_column("to_char(plate.plate_number,'%s')" 
                            % PLATE_NUMBER_SQL_FORMAT))),
                'remaining_well_volume': _plate_statistics.c.remaining_well_volume,
                'avg_remaining_volume': _plate_statistics.c.avg_well_remaining_volume,
                'min_remaining_volume': _plate_statistics.c.min_well_remaining_volume,
                'max_remaining_volume': _plate_statistics.c.max_well_remaining_volume,
                'min_molar_concentration': _plate_statistics.c.min_molar_concentration,
                'max_molar_concentration': _plate_statistics.c.max_molar_concentration,
                'min_mg_ml_concentration': _plate_statistics.c.min_mg_ml_concentration,
                'max_mg_ml_concentration': _plate_statistics.c.max_mg_ml_concentration,
                'plate_number': (literal_column(
                    "to_char(plate.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
                    .label('plate_number')),
                'first_date_screened': _plate_screening_statistics.c.first_date_screened,
                'last_date_screened': _plate_screening_statistics.c.last_date_screened,
                'status_date': _status_apilogs.c.date_time,
                'status_performed_by': _status_apilogs.c.name,
                'status_performed_by_username': _status_apilogs.c.username,
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
            # if for_screen_id:
            #     _subquery = \
            #         self.get_screen_librarycopyplate_subquery(for_screen_id)
            #     _subquery = _subquery.cte('screen_lcps')
            #     j = j.join(_subquery, _subquery.c.plate_id == _p.c.plate_id)
            # NOTE: not tracking loaded plates anymore
            # if loaded_for_screen_id:
            #     _subquery = \
            #         self.get_screen_loaded_librarycopyplate_subquery(
            #             loaded_for_screen_id)
            #     _subquery = _subquery.cte('screen_loaded_lcps')
            #     j = j.join(_subquery, _subquery.c.plate_id == _p.c.plate_id)
            if set(['avg_remaining_volume','min_remaining_volume',
                'max_remaining_volume']) | set(field_hash.keys()):
                j = j.outerjoin(_plate_statistics,_p.c.plate_id==_plate_statistics.c.plate_id)     
            if (set(['first_date_screened', 'last_date_screened'])
                | set(field_hash.keys()) ):
                j = j.outerjoin(_plate_screening_statistics,
                    _p.c.plate_id==_plate_screening_statistics.c.plate_id)
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

            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename,
                field_hash=field_hash, param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None))
                        
        except Exception, e:
            logger.exception('on get list')
            raise e   

    @classmethod
    def get_screen_librarycopyplate_subquery(cls, for_screen_id):
        
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
            .where(_screen.c.facility_id == for_screen_id))
        return screen_lcps

    # NOTE: not tracking loaded plates anymore
    # @classmethod
    # def get_screen_loaded_librarycopyplate_subquery(cls, for_screen_id):
    #     
    #     subquery = cls.get_screen_librarycopyplate_subquery(for_screen_id)
    #     _assay_plate = get_tables()['assay_plate']
    #     subquery = subquery.where(
    #         _assay_plate.c.screen_result_data_loading_id != None)
    #     return subquery

    def get_id(self, deserialized, validate=False, schema=None, **kwargs):
        id_kwargs = super(LibraryCopyPlateResource, self).get_id(
            deserialized, validate=validate, schema=schema, **kwargs)
        
        if 'plate_number' in id_kwargs:
            id_kwargs['plate_number'] = str(int(id_kwargs['plate_number'])).zfill(5)
        return id_kwargs
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def batch_edit(self, request, **kwargs):
        '''
        Batch edit is like patch, except:
        - bypasses the "dispatch" framework call
        - can use the "search_data" to find the plates
        -- must be authenticated and authorized
        '''
        self.is_authenticated(request)
        if not self._meta.authorization._is_resource_authorized(
                self._meta.resource_name,request.user,'write'):
            raise ImmediateHttpResponse(
                response=HttpForbidden(
                    'user: %s, permission: %s/%s not found' 
                    % (request.user,self._meta.resource_name,'write')))

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        logger.info('batch_edit plates...')
        
        schema = super(LibraryCopyPlateResource, self).build_schema()
        deserialized = self.deserialize(request)
        
        data = deserialized.get('data', None)
        if data is None:
            raise BadRequest('must submit "data" on the request')
        
        plate_info_data = data.get('plate_info', None)
        location_data = data.get('plate_location', None)
        if plate_info_data is not None and location_data is not None:
            raise BadRequest('batch info: edit only one of %r'
                % ['plate_info', 'plate_location'])
        elif plate_info_data is None and location_data is None:
            raise BadRequest('"data" must contain one of %r'
                % ['plate_info', 'plate_location'])
        
        if plate_info_data is not None:
            # Find all of the plates
            kwargs['search_data'] = deserialized.pop('search_data', {})
            kwargs['includes'] = ['plate_id']
            original_data = self._get_list_response(request, **kwargs)
            
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
            
            new_data = self._get_list_response(request, **kwargs)
            
            parent_log = self.make_log(request)
            parent_log.key = self._meta.resource_name
            parent_log.uri = self._meta.resource_name
            parent_log.save()
            kwargs['parent_log'] = parent_log
            logs = self.log_patches(request, original_data, new_data,**kwargs)
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
            
        location_data = data.get('plate_location', None)
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
            kwargs['search_data'] = deserialized.pop('search_data', {})
            original_data = self._get_list_response(request, **kwargs)
            
            if not original_data:
                raise HttpNotFound
            
            logger.info('get the plate objects for search data: %d', len(original_data))
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
                       str(plate_range[0]).zfill(5),
                       str(plate_range[-1]).zfill(5) ))
    
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
            response = self.get_platelocation_resource().patch_detail(
                request, 
                data=location_data,
                **kwargs)
            return response
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        
        logger.debug('patch obj, deserialized: %r', deserialized)
        schema = kwargs['schema']
        fields = schema['fields']

        id_kwargs = self.get_id(deserialized, **kwargs)
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
        
        initializer_dict = self.parse(deserialized,create=False)
        errors = self.validate(initializer_dict, patch=True)
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
    

class ScreenResultResource(DbApiResource):

    class Meta:
    
        queryset = ScreenResult.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'screenresult'
        ordering = []
        filtering = {}
        serializer = ScreenResultSerializer()
        object_class = dict
        max_limit = 10000
        
    def __init__(self, **kwargs):

        self.scope = 'fields.screenresult'
        super(ScreenResultResource, self).__init__(**kwargs)
        
        with get_engine().connect() as conn:
            try:
                conn.execute(text('select * from "well_query_index"; '));
                logger.debug('The well_query_index table exists')
            except Exception as e:
                logger.info('creating the well_query_index table')
                self._create_well_query_index_table(conn)
            try:
                conn.execute(text(
                    'select * from "well_data_column_positive_index"; '));
                logger.debug('The well_data_column_positive_index table exists')
            except Exception as e:
                logger.info('creating the well_data_column_positive_index table')
                self._create_well_data_column_positive_index_table(conn)
        self.reagent_resource = None
        self.screen_resource = None
        self._well_data_column_positive_index = None

    def get_well_data_column_positive_index_table(self):
        ''' Work around method; create the well_data_column_positive_index 
        table for the Aldjemy bridge '''
        if not 'well_data_column_positive_index' in get_tables():
            if not self._well_data_column_positive_index:
                self._well_data_column_positive_index = sqlalchemy.sql.schema.Table(
                    'well_data_column_positive_index', 
                    aldjemy.core.get_meta(),
                    sqlalchemy.Column('well_id', sqlalchemy.String),
                    sqlalchemy.Column('data_column_id', sqlalchemy.Integer))
            return self._well_data_column_positive_index
        return get_tables()['well_data_column_positive_index']

    def get_reagent_resource(self):
        if not self.reagent_resource:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource
    
    def get_screen_resource(self):
        if not self.screen_resource:
            self.screen_resource = ScreenResource()
        return self.screen_resource
    
    def _create_well_query_index_table(self, conn):
       
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
                'CREATE TABLE well_data_column_positive_index ('
                ' "well_id" text NOT NULL REFERENCES "well" ("well_id") '
                '    DEFERRABLE INITIALLY DEFERRED,'
                ' "data_column_id" integer NOT NULL ' 
                ' REFERENCES "data_column" ("data_column_id") '
                '    DEFERRABLE INITIALLY DEFERRED'
                ');'
            ));
            logger.info('the well_data_column_positive_index table created')
        except Exception, e:
            logger.info((
                'Exception: %r on trying to create the '
                'well_data_column_positive_index table,'
                'note that this is normal if the table already exists '
                '(PostgreSQL <9.1 has no "CREATE TABLE IF NOT EXISTS"'), e)
    
    @transaction.atomic()        
    def clear_cache(self, all=False, by_date=None, by_uri=None, by_size=False):
        
        logger.info(
            'clear_cache called: all: %r, by_date: %r, by_uri: %r, by_size: %r',
            all, by_date, by_uri, by_size)

        # FIXME: need test cases for this
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
            get_engine().execute(delete(self.get_well_data_column_positive_index_table()))
            logger.info(
                'cleared all cached well_data_column_positive_indexes')
        
        self.get_screen_resource().clear_cache()
            
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
        _rv = self.bridge['result_value']
        field_name = field_information['key']
        data_column_type = field_information.get('data_type') 
        # TODO: column to select: use controlled vocabulary
        column_to_select = None
        if(data_column_type in ['numeric', 'decimal', 'integer']):  
            column_to_select = 'numeric_value'
        else:
            column_to_select = 'value'
        
        rv_select = select([column(column_to_select)]).select_from(_rv)
        rv_select = rv_select.where(
            _rv.c.data_column_id == field_information['data_column_id'])
        rv_select = rv_select.where(_rv.c.well_id == text('assay_well.well_id'))
        # FIXME: limit to rv's to 1 - due to the duplicated result value issue
        rv_select = rv_select.limit(1)  
        rv_select = rv_select.label(field_name)
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
        logger.info('build screenresult query')
        DEBUG_SCREENRESULT = False or logger.isEnabledFor(logging.DEBUG)
        manual_field_includes = set(param_hash.get('includes', []))
        manual_field_includes.add('assay_well_control_type')
            
        # general setup
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, **extra_params)
                              
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
        if filter_excluded:
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
        if 'plate_number' in base_fields:
            sub_columns['plate_number'] = (literal_column(
                "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
                .label('plate_number'))
        for key,col in sub_columns.items():
            if key not in base_columns:
                if DEBUG_SCREENRESULT: 
                    logger.info('adding reagent column: %r...', key)
                base_columns[key] = col
        
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
        _wqx = select(['id', 'well_id']).select_from(_wellQueryIndex)
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

        j = j.join(excluded_cols_select, 
            excluded_cols_select.c.well_id == _aw.c.well_id, isouter=True)
        # Using nested selects
        for fi in [
            fi for fi in field_hash.values() 
                if fi.get('is_datacolumn', None)]:
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
        sub_columns['plate_number'] = (literal_column(
            "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
            .label('plate_number'))
        for key,col in sub_columns.items():
            if key not in columns:
                columns[key] = col

        # Force the query to use well_query_index.well_id
        columns['well_id'] = literal_column('wqx.well_id').label('well_id')
        if DEBUG_SCREENRESULT: logger.info('columns: %r', columns.keys())
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

        show_mutual_positives = bool(param_hash.get('show_mutual_positives', False))
        if show_mutual_positives is True:
            extra_params['mutual_positive'] = None
        other_screens = param_hash.get('other_screens', [])
        if other_screens:
            extra_params['other_screens'] = None
        schema = self.build_schema(
            screenresult=screenresult,
            show_mutual_positives=show_mutual_positives,
            other_screens=other_screens)

#         filename = self._get_filename(schema, kwargs)
        content_type = self.get_accept_content_type(
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
                def title_function(key):
                    return field_hash[key]['title']
            rowproxy_generator = None
            if ( use_vocab or content_type != JSON_MIMETYPE):
                # NOTE: xls export uses vocab values
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

    def create_dc_positive_index(self, screen_result_id):
        # the well_data_column_positive_index has been cleared, recreate
        _aw = self.bridge['assay_well']
        _sr = self.bridge['screen_result']
        _dc = self.bridge['data_column']
        base_stmt = join(
            _aw, _dc, _aw.c.screen_result_id == _dc.c.screen_result_id)
        base_stmt = select([
                literal_column("well_id"),
                literal_column('data_column_id')
            ]).select_from(base_stmt)            
        base_stmt = base_stmt.where(_aw.c.is_positive)
        base_stmt = base_stmt.where(
            _dc.c.data_type.in_([
                'boolean_positive_indicator',
                'partition_positive_indicator', 
                'confirmed_positive_indicator']))
        base_stmt = base_stmt.order_by(
            _dc.c.data_column_id, _aw.c.well_id)
        insert_statement = (
            insert(self.get_well_data_column_positive_index_table())
                .from_select(['well_id', 'data_column_id'], base_stmt))
        logger.info(
            'mutual pos insert statement: %r',
            str(insert_statement.compile(
                dialect=postgresql.dialect(),
                compile_kwargs={"literal_binds": True})))
        get_engine().execute(insert_statement)
        logger.info('mutual pos insert statement, executed.')

    
    def get_mutual_positives_columns(self, screen_result_id):
        
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
            
            # Recreate the well_data_column_positive_index:
            # Note: because this is lazy recreated, we are creating the whole 
            # index each time; could just recreate for the specific screen
            count_stmt = self.get_well_data_column_positive_index_table()
            count_stmt = select([func.count()]).select_from(count_stmt)
            count = 0
            with get_engine().begin() as conn:
                count = int(conn.execute(count_stmt).scalar())        
                logger.info('well_data_column_positive_index count: %r', count)
            if count == 0:
                self.create_dc_positive_index(screen_result_id)
            
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
            _wdc = self.get_well_data_column_positive_index_table().alias('wdc')
            _wdc1 = self.get_well_data_column_positive_index_table().alias('wdc1')
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

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        other_screens = param_hash.get('other_screens', [])
        show_mutual_positives = bool(param_hash.get('show_mutual_positives', False))

        try:
            screenresult = ScreenResult.objects.get(
                screen__facility_id=facility_id)
            return self.build_response(
                request, self.build_schema(
                    screenresult, 
                    show_mutual_positives=show_mutual_positives,
                    other_screens=other_screens), 
                **kwargs)
        except ObjectDoesNotExist, e:
            raise Http404(
                'no screen result found for facility id: %r' % facility_id)
    
    def build_schema(
            self, screenresult=None, show_mutual_positives=False,
            other_screens=None, user=None):
        
        screen_facility_id = None
        if screenresult:
            screen_facility_id = screenresult.screen.facility_id

        cache_key = 'screenresult_schema_%s_mutual_pos_%s' \
            % (screen_facility_id, show_mutual_positives)
        logger.info('build screenresult schema: %s', cache_key)
        data = cache.get(cache_key)
        
        if data is None:
            logger.info('not cached: %s', cache_key)
            data = super(ScreenResultResource, self).build_schema(user=user)
            
            if screenresult:
                try:
                    well_schema = self.get_reagent_resource().build_schema(
                        screenresult.screen.screen_type, user=user)
                    newfields = {}
                    newfields.update(well_schema['fields'])
                    for key,field in newfields.items():
                        field['visibility'] = []
                    newfields.update(data['fields'])
                    
                    if screenresult.screen.screen_type == 'small_molecule':
                        if 'compound_name' in newfields:
                            newfields['compound_name']['visibility'] = ['l','d']
                            newfields['compound_name']['ordinal'] = 12
                    elif screenresult.screen.screen_type == 'rnai':
                        if 'facility_entrezgene_id' in newfields:
                            newfields['facility_entrezgene_id']['visibility'] = ['l','d']
                            newfields['facility_entrezgene_id']['ordinal'] = 10
                        if 'facility_entrezgene_symbols' in newfields:
                            newfields['facility_entrezgene_symbols']['visibility'] = ['l','d']
                            newfields['facility_entrezgene_symbols']['ordinal'] = 11
                    max_ordinal = 0
                    for fi in newfields.values():
                        if fi.get('ordinal', 0) > max_ordinal:
                            max_ordinal = fi['ordinal']
                    logger.info('map datacolumn definitions into field information definitions...')
                    for i, dc in enumerate(
                        DataColumn.objects
                            .filter(screen_result=screenresult)
                            .order_by('ordinal')):
                        (columnName, _dict) = \
                            DataColumnResource._create_datacolumn_from_orm(dc)
                        _dict['ordinal'] = max_ordinal + dc.ordinal + 1
                        _dict['visibility'] = ['l', 'd']
                        newfields[columnName] = _dict
                    
                    max_ordinal += dc.ordinal + 1
                    _current_sr = None
                    _ordinal = max_ordinal
                    logger.info('get mutual positives columns...')
                    mutual_positives_columns = \
                        self.get_mutual_positives_columns(
                            screenresult.screen_result_id)
                    logger.info('iterate mutual positives columns %r...',
                        show_mutual_positives)
                    if show_mutual_positives is True:
                        for i, dc in enumerate(DataColumn.objects
                            .filter(data_column_id__in=mutual_positives_columns)
                            .order_by('screen_result__screen__facility_id', 'ordinal')):
        
                            if _current_sr != dc.screen_result.screen_result_id:
                                _current_sr = dc.screen_result.screen_result_id
                                max_ordinal = _ordinal
                            _ordinal = max_ordinal + dc.ordinal + 1
                            (columnName, _dict) = \
                                DataColumnResource._create_datacolumn_from_orm(dc)
                            _dict['ordinal'] = _ordinal
                            if show_mutual_positives:
                                _dict['visibility'] = ['l', 'd']
                            
                            newfields[columnName] = _dict
                            
                            
                            logger.info('created mutual positive column:%s: %r', dc.name, _dict)    
                    max_ordinal = _ordinal
                    # Column selectors for "other screens"
                    for i, screen in enumerate(Screen.objects
                        .exclude(screen_id=screenresult.screen_id)
                        .filter(screen_type=screenresult.screen.screen_type)
                        .exclude(screenresult__isnull=True)
                        .order_by('facility_id')):
                    
                        (columnName, _dict) = \
                            self.create_otherscreen_field(screen)
                        _dict['scope'] = 'otherscreen.datacolumns'
                        _dict['ordinal'] = max_ordinal + i
                        newfields[columnName] = _dict
                        logger.debug('created other screenresult column: %r', _dict)    
                        
                    data['fields'] = newfields
                except Exception, e:
                    logger.exception('on build schema')
                    raise e  
                
            logger.info('build screenresult schema done')
            cache.set(cache_key, data)
        else:
            logger.info('using cached: %s', cache_key)
            
        if other_screens is not None:
            logger.info('include other screen columns: %r', other_screens)
            newfields = {}
            newfields.update(data['fields'])
            max_field = max(newfields.itervalues(), key=lambda x: int(x['ordinal']))
            max_ordinal = int(max_field['ordinal'])

            for other_screen_facility_id in other_screens:
                colname = 'screen_%s' % other_screen_facility_id
                if colname in newfields and newfields[colname].get('is_screen_column',False):
                    for i,dc in enumerate(DataColumn.objects
                        .filter(screen_result__screen__facility_id=other_screen_facility_id)):
                        (columnName, _dict) = \
                            DataColumnResource._create_datacolumn_from_orm(dc)
                        _ordinal = max_ordinal + dc.ordinal + 1
                        _dict['ordinal'] = _ordinal
                        _dict['visibility'] = ['l','d']
                        
                        newfields[columnName] = _dict
            data['fields'] = newfields
            logger.info('newfields: %r', newfields.keys())
        return data
    
    def create_otherscreen_field(self,screen):
        columnName = "screen_%s" % screen.facility_id
        _dict = {}
        _dict['title'] = '%s (%s)' % (screen.facility_id, screen.title) 
        _dict['description'] =  _dict['title']
        _dict['is_screen_column'] = True
        _dict['key'] = columnName
        _dict['screen_facility_id'] = screen.facility_id
        _dict['visibility'] = ['api']
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
            
            screen_log = self.make_log(request, **kwargs)
            screen_log.ref_resource_name = 'screen'
            screen_log.key = screen.facility_id
            screen_log.uri = '/'.join(['screen', screen_log.key])
            screen_log.save()
        else:
            raise BadRequest('no screen result to delete')

        return tastypie.http.HttpNoContent()
            
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def patch_detail(self, request, **kwargs):
        if 'screen_facility_id' not in kwargs:
            raise BadRequest('screen_facility_id is required')
        screen_facility_id = kwargs['screen_facility_id']
        
        data = self.deserialize(request)
        meta = data[API_RESULT_META]
        columns = data['fields']
        result_values = data[API_RESULT_DATA]
        
        if ('screen_facility_id' in meta 
            and meta['screen_facility_id'] != screen_facility_id):
            logger.warn(
                'screen_facility_id in file %r does not match url: %r'
                % (meta['screen_facility_id'], screen_facility_id))
        
        # Clear cache: note "all" because mutual positive columns can change
        # self.clear_cache(by_uri='/screenresult/%s' % screen_facility_id)
        self.clear_cache(all=True)
        
        schema = kwargs['schema']
        id_attribute = schema['id_attribute']

        screen = Screen.objects.get(facility_id=screen_facility_id)
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

            logger.info('Create screen result for %r', screen_facility_id)
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
            
            logger.info('Create screen result data columns for %r', screen_facility_id)
            sheet_col_to_datacolumn = {}
            derived_from_columns_map = {}
            errors = {}
            for sheet_column, column_info in columns.items():
                if ( column_info.get('screen_facility_id',None)
                    and column_info['screen_facility_id'] != screen_facility_id ):
                    logger.info('skipping column for screen_facility_id: %r', 
                        column_info.get('screen_facility_id') )
                    continue;
                    
                column_info['screen_result'] = screen_result
                try:
                    dc = DataColumnResource().patch_obj(request, column_info, **kwargs)
                    sheet_col_to_datacolumn[sheet_column] = dc
                    derived_from_columns = column_info.get(
                        'derived_from_columns', None)
                    if derived_from_columns:
                        derived_from_columns = [
                            x.strip().upper() for x in derived_from_columns.split(',')] 
                        derived_from_columns_map[sheet_column] = \
                            derived_from_columns
                except ValidationError, e:
                    errors.update({ sheet_column: e.errors }) 
            logger.info(
                'create derived_from_columns: %r',derived_from_columns_map)
            logger.info(
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
        screenresult_log.diffs.update({ 
            'created_by': [None,adminuser.username],  
            'experimental_well_count': [None,screen_result.experimental_well_count],
            'replicate_count': [None,screen_result.replicate_count],
            'channel_count': [None,screen_result.channel_count],
        })
        screenresult_log.save()
        
        if screen.project_phase != 'annotation':
            logger.info('pre-generate the mutual positives index...')
            self.get_mutual_positives_columns(screen_result.screen_result_id)
            logger.info('done - pre-generate the mutual positives index')
             
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
    
    def create_result_value(
            self, colname, value, dc, well, initializer_dict, 
            assay_well_initializer):
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

        logger.debug(
            'create result value: %r, colname: %r, dc: %r %r',
            value, colname, dc.name, dc.data_type)
        key = '%s-%s' % (well_id, colname)
        rv_initializer['data_column_id'] = dc.data_column_id

        if dc.data_type == 'numeric':
            if  dc.decimal_places > 0:
                # parse, to validate
                parse_val(value, key, 'float')
                # record the raw val, to save all digits (precision)
                rv_initializer['numeric_value'] = value
            else:
                rv_initializer['numeric_value'] = \
                    parse_val(value, key, 'integer')
        else:
            rv_initializer['value'] = value

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
                logger.debug(
                    ('non experimental well, not considered for positives:'
                     'well: %r, col: %r, type: %r, val: %r'),
                    well_id, colname, dc.data_type, value)
            elif rv_initializer['is_exclude'] is True:
                logger.info(
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
        
        logger.debug('rv_initializer: %r', rv_initializer)
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
                assay_well_file, fieldnames=assay_well_fieldnames, delimiter=str(','),
                lineterminator="\n")
            
            rows_created = 0
            rvs_to_create = 0
            plates_max_replicate_loaded = {}
            logger.info('write temp file result values for screen: %r ...',
                screen_result.screen.facility_id)
            errors = {}
            while True:
                try: 
                    # iterating will trigger parsing
                    result_row = result_values.next()
                    logger.debug('result_row: %r', result_row)
                    initializer_dict = { 
                        fieldname:PSYCOPG_NULL for fieldname in fieldnames}
                    assay_well_initializer = { 
                        fieldname:PSYCOPG_NULL for fieldname in assay_well_fieldnames}
                    for meta_field in meta_columns:
                        if meta_field in result_row:
                            initializer_dict[meta_field] = result_row[meta_field]

                    well = Well.objects.get(well_id=result_row['well_id']) 
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
                        if colname in meta_columns:
                            continue
                        if val is None:
                            continue
                        if colname not in sheet_col_to_datacolumn:
                            logger.debug('extra column found in the Data sheet: %r', colname)
                            continue      
                        dc = sheet_col_to_datacolumn[colname]
                        rv_initializer = self.create_result_value(
                            colname, val, dc, well, initializer_dict, 
                            assay_well_initializer)
                        
                        
                        if (dc.is_derived is False
                            and well.library_well_type == 'experimental'
                            and rv_initializer['is_exclude'] is not True):
                            max_replicate = \
                                plates_max_replicate_loaded.get(well.plate_number, 0)
                            if dc.replicate_ordinal > max_replicate:
                                plates_max_replicate_loaded[well.plate_number] = \
                                    dc.replicate_ordinal
                        else:
                            logger.debug(('not counted for replicate: well: %r, '
                                'type: %r, initializer: %r'), 
                                well.well_id, well.library_well_type, rv_initializer)        
                        writer.writerow(rv_initializer)
                        rvs_to_create += 1
                
                    assay_well_writer.writerow(assay_well_initializer)
                    rows_created += 1
                    if rows_created % 10000 == 0:
                        logger.info('wrote %d result rows to temp file', rows_created)

                except ValidationError, e:
                    logger.info('validation error: %r', e)
                    errors.update(e.errors) 
                except StopIteration, e:
                    break
            
            if errors:
                logger.warn('errors: %r', errors)
#                 raise ValidationError( errors={ 'result_values': errors })
                raise ValidationError(errors=errors)
            
            if not rvs_to_create:
                raise ValidationError( errors={ 'result_values': 'no result values were parsed' })

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
            
        if plates_max_replicate_loaded:
            plates_max_replicate_loaded = sorted(
                plates_max_replicate_loaded.values())
            logger.info(
                'plates_max_replicate_loaded: %r', plates_max_replicate_loaded)
            screen_result.screen.min_data_loaded_replicate_count = \
                plates_max_replicate_loaded[0]
            screen_result.screen.max_data_loaded_replicate_count = \
                plates_max_replicate_loaded[-1]
        else:
            screen_result.screen.min_data_loaded_replicate_count = 1
            screen_result.screen.max_data_loaded_replicate_count = 1
        screen_result.screen.save()
        logger.info('screen_result.screen.max_data_loaded_replicate_count: %r',
            screen_result.screen.max_data_loaded_replicate_count)
        
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
            
        

class DataColumnResource(DbApiResource):

    class Meta:

        queryset = DataColumn.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'datacolumn'
        ordering = []
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
        super(DataColumnResource, self).__init__(**kwargs)

    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<data_column_id>\d+)%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

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
        schema = super(DataColumnResource, self).build_schema()
        
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)

        try:
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            # Add fields required to build the system representation
            manual_field_includes.update([
                'name','data_type','decimal_places','ordinal',
                'screen_facility_id'])
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)

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
            
            max_ordinal = len(field_hash)
            def create_dc_generator(generator):
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
                    if generator:
                        cursor = generator(cursor)
                    for row in cursor:
                        yield DataColumnRow(row)
                return datacolumn_fields_generator
            
            rowproxy_generator = create_dc_generator(rowproxy_generator)
            # specific setup 

            _dc = self.bridge['data_column']
            _sr = self.bridge['screen_result']
            _screen = self.bridge['screen']
            _dc_derived_link = self.bridge['data_column_derived_from_columns']
            _dc_derived = _dc.alias('dc_derived')
            
            base_query_tables = [
                'data_column', 'screen']
            
            custom_columns = {
                'derived_from_columns': (
                    select([func.array_to_string(
                        func.array_agg(literal_column('name')), LIST_DELIMITER_SQL_ARRAY)])
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

            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            stmt = stmt.order_by('ordinal')
            
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

        # FIXME: migrate the data_type of the datacolumn so that it is not 
        # determined at runtime
        dc_data_type = row_or_dict['data_type']
        _dict['assay_data_type'] = dc_data_type
        if dc_data_type == 'numeric':
            if row_or_dict.has_key('decimal_places'):
                decimal_places = row_or_dict['decimal_places']
                if decimal_places > 0:
                    _dict['data_type'] = 'decimal'
                    _dict['display_options'] = (
                         '{ "decimals": %s }' % decimal_places)
                else:
                    _dict['data_type'] = 'integer'
        elif dc_data_type in DataColumnResource.data_type_lookup:
            _dict['data_type'] = dc_data_type                            
        _dict['ordinal'] = max_ordinal + row_or_dict['ordinal']
        for key in row_or_dict.keys():
            if key not in _dict:
                _dict[key] = row_or_dict[key]
        return _dict

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
        # FIXME: migrate the data_type of the datacolumn so that it is not 
        # determined at runtime
        if dc.data_type == 'numeric':
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
        return (columnName, _dict)

    @write_authorization
    @un_cache
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):

        logger.debug('patch_obj %s', deserialized)
        
        schema = kwargs['schema']
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

        errors = self.validate(initializer_dict, patch=False)
        if errors:
            raise ValidationError(errors)
                
        try:
            
            logger.debug('initializer dict: %s', initializer_dict)
            if 'screen_result' not in deserialized:
                raise BadRequest('screen_result is required')
            data_column = DataColumn(
                screen_result=deserialized['screen_result'])
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
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'copywell'
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
            # schema = super(CopyWellResource, self).build_schema()
#         filename = self._get_filename(schema, kwargs)
        well_id = param_hash.pop('well_id', None)
        if well_id:
            param_hash['well_id__eq'] = well_id
        copy_name = param_hash.pop('copy_name', None)
        if copy_name:
            param_hash['copy_name__eq'] = copy_name
        library_short_name = param_hash.pop('library_short_name', None)
        if library_short_name:
            param_hash['library_short_name__eq'] = library_short_name

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
  
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)

            # TODO: remove this restriction if the query can be optimized
            if filter_expression is None:
                raise InformationError(
                    key='Input filters ',
                    msg='Please enter a filter expression to begin')
                                  
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
                'screening_count': (
                    (func.coalesce(_p.c.screening_count,0) 
                        + func.coalesce(_p.c.cplt_screening_count,0))),
                # TODO: 20170104 - track cherry_pick_screening_count as well
                #  
                
                
                
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
    @un_cache
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        # TODO: optimize for list inputs (see well.patch)
        logger.debug('patch_obj %s', deserialized)

        schema = kwargs['schema']
        fields = schema['fields']
        initializer_dict = {}
        id_kwargs = self.get_id(deserialized, validate=True, **kwargs)
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
    def deallocate_cherry_pick_volumes(
        self, cpr, lab_cherry_picks_to_deallocate, parent_log):

        copywells_deallocated = []
        plates_adjusted = set()
        # find copy-wells
        for lcp in lab_cherry_picks_to_deallocate:
            copy = lcp.copy
            if copy is None:
                raise ProgrammingError('LabCherryPick is already deallocated: %r', lcp)
            plate = Plate.objects.get(
                plate_number=lcp.source_well.plate_number, 
                copy=lcp.copy)
            plates_adjusted.add(plate)
            copy_name = self.COPYWELL_KEY.format(
                well_id=lcp.source_well.well_id,
                library_short_name=copy.library.short_name,
                copy_name=copy.name)
            try:
                copywell = CopyWell.objects.get(
                    well=lcp.source_well, copy=lcp.copy)
                
            except ObjectDoesNotExist:
                logger.error('copywell to deallocate not located: %r', copy_name)
                
            log = self.make_child_log(parent_log)
            log.key = '/'.join([copy.library.short_name, copy.name, copywell.well_id])
            log.uri = '/'.join([log.ref_resource_name,log.key])
            new_volume = copywell.volume + cpr.transfer_volume_per_well_approved
            if copywell.cherry_pick_screening_count > 0:
                current_cp_screening_count = copywell.cherry_pick_screening_count
                new_cp_screening_count =  current_cp_screening_count - 1
                log.diffs = {
                    'volume': [copywell.volume, new_volume],
                    'cherry_pick_screening_count': [
                        current_cp_screening_count, 
                        new_cp_screening_count]
                    }
                log.parent_log = parent_log
                logger.info('copywell adjusted: %r, %r', log, log.diffs)
                log.save()
                copywell.cherry_pick_screening_count = new_cp_screening_count
            # adjust volume
            copywell.volume = new_volume
            copywell.save()
            copywells_deallocated.append(copywell)
        logger.info('copywells_deallocated: %r', copywells_deallocated)    
        logger.info('plates_adjusted: %r', plates_adjusted)
        for plate in plates_adjusted:
            plate_log = self.get_plate_resource().make_child_log(parent_log)
            plate_log.key = '/'.join([
                copy.library.short_name, plate.copy.name, str(plate.plate_number)])
            plate_log.uri = '/'.join([log.ref_resource_name,log.key])
            if plate.cplt_screening_count < 1:
                logger.warn('deallocation: plate: %r, cplt_screening_count is already 0',
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
            'cpr': cpr.cherry_pick_request_id,
            API_MSG_COPYWELLS_DEALLOCATED: len(copywells_deallocated),
            'plates_adjusted': len(plates_adjusted)
         }
        
        
    def reserve_cherry_pick_volumes(self, cpr, fulfillable_lcps, parent_log):
        
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
            
        logger.info('plates_adjusted: %r', plates_adjusted)
        for plate in plates_adjusted:
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
            API_RESULT_META: { 
                'cpr': cpr.cherry_pick_request_id,
                'copywells allocated': len(copywells_adjusted),
                'plates_adjusted': len(plates_adjusted)
             }
         }

# class CherryPickWellResource(DbApiResource):        
# 
#     class Meta:
#     
#         queryset = CherryPickWell.objects.all().order_by('well_id')
#         authentication = MultiAuthentication(BasicAuthentication(),
#                                              SessionAuthentication())
#         authorization = UserGroupAuthorization()
#         resource_name = 'cherrypickwell'
#         ordering = []
#         filtering = {}
#         serializer = LimsSerializer()
#         always_return_data = True 
# 
#     def __init__(self, **kwargs):
# 
#         super(CherryPickWellResource, self).__init__(**kwargs)
# 
#     def prepend_urls(self):
#         
#         return [
#             url(r"^(?P<resource_name>%s)/schema%s$" 
#                 % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('get_schema'), name="api_get_schema"),
#             url(r"^(?P<resource_name>%s)"
#                 r"/(?P<cherry_pick_request_id>\d)%s$" 
#                     % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('dispatch_list'), name="api_dispatch_list"),
#             url(r"^(?P<resource_name>%s)"
#                 r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
#                     % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#         ]
# 
#     @read_authorization
#     def get_detail(self, request, **kwargs):
# 
#         kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
#         kwargs['is_for_detail'] = True
#         return self.build_list_response(request, **kwargs)
#         
#     @read_authorization
#     def get_list(self, request, **kwargs):
# 
#         kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
#         return self.build_list_response(request, **kwargs)
# 
#     def build_list_response(self, request, schema=None, **kwargs):
# 
#         param_hash = self._convert_request_to_dict(request)
#         param_hash.update(kwargs)
#         is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
#         use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
#         use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
#         if is_data_interchange:
#             use_vocab = False
#             use_titles = False
#         
#         is_for_detail = kwargs.pop('is_for_detail', False)
#         if schema is None:
#             raise Exception('schema not initialized')
#             # schema = super(CherryPickRequestResource, self).build_schema()
#         filename = self._get_filename(schema, kwargs)
#         cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
#         if cherry_pick_request_id:
#             param_hash['cherry_pick_request_id__eq'] = cherry_pick_request_id
# 
#         try:
#             
#             # general setup
#           
#             manual_field_includes = set(param_hash.get('includes', []))
# 
#             (filter_expression, filter_hash, readable_filter_hash) = \
#                 SqlAlchemyResource.build_sqlalchemy_filters(
#                     schema, param_hash=param_hash)
#                                   
#             order_params = param_hash.get('order_by', [])
# 
#             logger.info('fields: %r', schema['fields'].keys())
#             field_hash = self.get_visible_fields(
#                 schema['fields'], filter_hash.keys(), manual_field_includes,
#                 param_hash.get('visibilities'),
#                 exact_fields=set(param_hash.get('exact_fields', [])),
#                 order_params=order_params)
#             order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
#                 order_params, field_hash)
#              
#             rowproxy_generator = None
#             if use_vocab is True:
#                 rowproxy_generator = \
#                     DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
#  
#             # specific setup 
# 
#             base_query_tables = ['cherry_pick_well']
#             _cp_well = self.bridge['cherry_pick_well']
#             _well = self.bridge['well']
#             _copy = self.bridge['copy']
#             _cpr = self.bridge['cherry_pick_request2']
#             _screen = self.bridge['screen']
#             _library = self.bridge['library']
#             
#             custom_columns = {
#                 'cherry_pick_request_id': _cp_well.c.cherry_pick_request_id,
#                 'library_short_name': _library.c.short_name,
#                 'plate_number': _well.c.plate_number,
#                 'well_id': _cp_well.c.source_well_id,
#                 'well_name': _well.c.well_name,
#                 'library_well_type': _well.c.library_well_type, 
#                 'copy_usage_type': _copy.c.usage_type,
#                 'plate_status': literal_column("'tbd'"),
#                 'copy_name': _copy.c.name,
#                 'volume': _cpr.c.transfer_volume_per_well_approved, 
#             }
#             
#             columns = self.build_sqlalchemy_columns(
#                 field_hash.values(), base_query_tables=base_query_tables,
#                 custom_columns=custom_columns)
# 
#             j = _cp_well
#             j = j.join(_cpr, _cp_well.c.cherry_pick_request_id==_cpr.c.id)
#             j = j.join(_screen, _cpr.c.screen_id == _screen.c.screen_id)
#             j = j.join(_well, _cp_well.c.source_well_id==_well.c.well_id)
#             j = j.join(_copy, _cp_well.c.copy_id==_copy.c.copy_id, isouter=True)
#             j = j.join(_library, _well.c.library_id==_library.c.library_id)
#             
#                         
#             stmt = select(columns.values()).select_from(j)
#             # general setup
#              
#             (stmt, count_stmt) = self.wrap_statement(
#                 stmt, order_clauses, filter_expression)
# 
#             compiled_stmt = str(stmt.compile(
#                 dialect=postgresql.dialect(),
#                 compile_kwargs={"literal_binds": True}))
#             logger.info('compiled_stmt %s', compiled_stmt)
#             
#             if not order_clauses:
#                 stmt = stmt.order_by(
#                     nullslast(desc(column('cherry_pick_request_id'))))
#             
#             title_function = None
#             if use_titles is True:
#                 def title_function(key):
#                     return field_hash[key]['title']
#             if is_data_interchange:
#                 title_function = DbApiResource.datainterchange_title_function(
#                     field_hash,schema['id_attribute'])
#             
#             return self.stream_response_from_statement(
#                 request, stmt, count_stmt, filename,
#                 field_hash=field_hash,
#                 param_hash=param_hash,
#                 is_for_detail=is_for_detail,
#                 rowproxy_generator=rowproxy_generator,
#                 title_function=title_function, meta=kwargs.get('meta', None),
#                 use_caching=True)
# 
#         except Exception, e:
#             logger.exception('on get list')
#             raise e  

# class CherryPickRequest2Resource(DbApiResource):        
# 
#     class Meta:
#     
#         queryset = CherryPickRequest2.objects.all().order_by('well_id')
#         authentication = MultiAuthentication(BasicAuthentication(),
#                                              SessionAuthentication())
#         authorization = UserGroupAuthorization()
#         resource_name = 'cherrypickrequest2'
#         ordering = []
#         filtering = {}
#         serializer = LimsSerializer()
#         always_return_data = True 
# 
#     def __init__(self, **kwargs):
# 
#         super(CherryPickRequest2Resource, self).__init__(**kwargs)
#         
#         self.copy_well_resource = None
#         
#     def get_copywell_resource(self):
#         if self.copy_well_resource is None:
#             self.copy_well_resource = CopyWellResource()
#         return self.copy_well_resource
#         
# 
#     def prepend_urls(self):
#         
#         return [
#             url(r"^(?P<resource_name>%s)/schema%s$" 
#                 % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('get_schema'), name="api_get_schema"),
#             url(r"^(?P<resource_name>%s)"
#                 r"/(?P<id>[\d]+)%s$" 
#                     % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#             url(r"^(?P<resource_name>%s)"
#                 r"/(?P<id>[\d]+)/wells%s$" 
#                     % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('dispatch_cherry_pick_wells'), 
#                 name="api_dispatch_cherry_pick_wells"),
#             url(r"^(?P<resource_name>%s)"
#                 r"/(?P<id>[\d]+)/copy_search%s$" 
#                     % (self._meta.resource_name, trailing_slash()),
#                 self.wrap_view('dispatch_cherry_pick_copy_search'), 
#                 name="api_dispatch_cherry_pick_copy_search"),
# 
#         ]
# 
#     def dispatch_cherry_pick_wells(self, request, **kwargs):
#         kwargs['cherry_pick_request_id'] = kwargs.pop('id')
#         return CherryPickWellResource().dispatch('list', request, **kwargs)   
# 
#     def dispatch_cherry_pick_copy_search(self, request, **kwargs):
#       
#         # FIXME: only POST allowed here
#         
#           
#         param_hash = self._convert_request_to_dict(request)
#         param_hash.update(kwargs)
#         
#         cherry_pick_request_id = param_hash['id']
#         logger.info('copy_search for cherry_pick: %r...', cherry_pick_request_id)
#         
#         schema = super(CherryPickRequest2Resource, self).build_schema()
#         
#         cpr = CherryPickRequest2.objects.get(id=cherry_pick_request_id)
#         
#         if not cpr.wells.exists():
#             raise ValidationError(
#                 key='cherry_pick_wells', 
#                 msg='No wells have been submitted')
#             
#         copy_wells = self.get_copywell_resource()._get_list_response_internal(
#             **{
#                 'well_id__in': [w.source_well_id for w in cpr.wells.all()],
#                 'copy_usage_type': 'cherry_pick_source_plates',
#                 'plate_status': 'available',
#             })
#         
#         # TODO: perform the search algorithm
#         
#         return self.build_response(
#             request,  {'objects': copy_wells }, response_class=HttpResponse, **kwargs)
#         
#     @read_authorization
#     def get_detail(self, request, **kwargs):
# 
#         _id = kwargs.get('id', None)
#         if not _id:
#             raise Http404('must provide a cherry_pick_request _id parameter')
#         kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
#         kwargs['is_for_detail'] = True
#         return self.build_list_response(request, **kwargs)
#         
#     @read_authorization
#     def get_list(self, request, **kwargs):
# 
#         kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
#         return self.build_list_response(request, **kwargs)
# 
#     def build_list_response(self, request, schema=None, **kwargs):
# 
#         param_hash = self._convert_request_to_dict(request)
#         param_hash.update(kwargs)
#         is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
#         use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
#         use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
#         if is_data_interchange:
#             use_vocab = False
#             use_titles = False
#         
#         is_for_detail = kwargs.pop('is_for_detail', False)
#         if schema is None:
#             raise Exception('schema not initialized')
#         filename = self._get_filename(schema, kwargs)
#         _id = param_hash.pop('id', None)
#         if _id:
#             param_hash['id__eq'] = _id
# 
#         try:
#             
#             # general setup
#           
#             manual_field_includes = set(param_hash.get('includes', []))
#   
#             (filter_expression, filter_hash, readable_filter_hash) = \
#                 SqlAlchemyResource.build_sqlalchemy_filters(
#                     schema, param_hash=param_hash)
#                                   
#             order_params = param_hash.get('order_by', [])
#             field_hash = self.get_visible_fields(
#                 schema['fields'], filter_hash.keys(), manual_field_includes,
#                 param_hash.get('visibilities'),
#                 exact_fields=set(param_hash.get('exact_fields', [])),
#                 order_params=order_params)
#             order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
#                 order_params, field_hash)
#              
#             rowproxy_generator = None
#             if use_vocab is True:
#                 rowproxy_generator = \
#                     DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
#  
#             # specific setup 
#             base_query_tables = ['cherry_pick_request2', 'screen']
#             _cpr = self.bridge['cherry_pick_request2']
#             _cp_well = self.bridge['cherry_pick_well']
#             _screen = self.bridge['screen']
#             _user_cte = ScreensaverUserResource.get_user_cte().cte('users_cpr')
#             _lab_cte = ScreensaverUserResource.get_lab_affiliation_cte().cte('lab')
#             
#             _requestors = _user_cte.alias('requestors')
#             _lead_screeners = _user_cte.alias('lead_screeners')
#             lab_head_table = ScreensaverUserResource.get_lab_head_cte().cte('lab_heads')
#             
#             custom_columns = {
#                 'requested_by_username': _requestors.c.username,
#                 'requested_by_name': _requestors.c.name,
#                 'lab_name': lab_head_table.c.lab_name_full,
#                 'lab_head_username': lab_head_table.c.username,
#                 'lead_screener_name': _lead_screeners.c.name,
#                 'lead_screener_username': _lead_screeners.c.username,
# #                 'screen_type': literal_column(
#                 'is_completed': literal_column('0'),
#                 'number_plates': literal_column('0'),
#                 'number_plates_completed': literal_column('0'),
#                 'number_unfulfilled_lab_cherry_picks': literal_column('0'),
#                 'total_number_lcps': literal_column('0'),
#                 'plating_activity_date': literal_column('0'),
#                 'cherry_pick_wells': (
#                     select([func.array_to_string(
#                         func.array_agg(literal_column('source_well_id')), 
#                         LIST_DELIMITER_SQL_ARRAY) ])
#                     .select_from(
#                         select([_cp_well.c.source_well_id])
#                             .select_from(_cp_well)
#                             .where(_cp_well.c.cherry_pick_request_id==_cpr.c.id)
#                             .order_by('source_well_id').alias('inner_cpwells'))
#                         ),
#             }
#             
#             columns = self.build_sqlalchemy_columns(
#                 field_hash.values(), base_query_tables=base_query_tables,
#                 custom_columns=custom_columns)
# 
#             j = join(_cpr, _screen, _cpr.c.screen_id == _screen.c.screen_id)
#             j = j.join(
#                 _requestors, 
#                 _cpr.c.requested_by_id==_requestors.c.screensaver_user_id,
#                 isouter=True)
#             j = j.join(
#                 lab_head_table, 
#                 lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id,
#                 isouter=True)
#             j = j.join(
#                 _lead_screeners,
#                 _lead_screeners.c.screensaver_user_id==_screen.c.lead_screener_id,
#                 isouter=True)
#             
#             stmt = select(columns.values()).select_from(j)
#             # general setup
#              
#             (stmt, count_stmt) = self.wrap_statement(
#                 stmt, order_clauses, filter_expression)
#             
#             if not order_clauses:
#                 stmt = stmt.order_by(
#                     nullslast(desc(column('id'))))
#             
#             title_function = None
#             if use_titles is True:
#                 def title_function(key):
#                     return field_hash[key]['title']
#             if is_data_interchange:
#                 title_function = DbApiResource.datainterchange_title_function(
#                     field_hash,schema['id_attribute'])
#             
#             return self.stream_response_from_statement(
#                 request, stmt, count_stmt, filename,
#                 field_hash=field_hash,
#                 param_hash=param_hash,
#                 is_for_detail=is_for_detail,
#                 rowproxy_generator=rowproxy_generator,
#                 title_function=title_function, meta=kwargs.get('meta', None),
#                 use_caching=True)
#              
#         except Exception, e:
#             logger.exception('on get list')
#             raise e  
# 
#     @write_authorization
#     @un_cache
#     @transaction.atomic
#     def patch_obj(self, request, deserialized, **kwargs):
# 
#         schema = kwargs['schema']
#         fields = schema['fields']
#         initializer_dict = {}
# 
#         id_kwargs = self.get_id(deserialized, **kwargs)
#         patch = bool(id_kwargs)
#         initializer_dict = self.parse(deserialized, create=not patch)
#         errors = self.validate(initializer_dict, patch=patch)
#         if errors:
#             raise ValidationError(errors)
#         
#         if not patch:
#             _key = 'screen_facility_id'
#             _val = deserialized[_key]
#             try:
#                 screen = Screen.objects.get(facility_id=_val)
#                 initializer_dict['screen'] = screen
#             except ObjectDoesNotExist:
#                 raise ValidationError(
#                     key=_key,
#                     msg='does not exist: {val}'.format(val=_val))
# 
#         _key = 'requested_by_username'
#         _val = deserialized.get(_key, None)
#         if _val:
#             try:
#                 requested_by_user = ScreensaverUser.objects.get(username=_val)
#                 initializer_dict['requested_by'] = requested_by_user
#             except ObjectDoesNotExist:
#                 raise ValidationError(
#                     key=_key,
#                     msg='does not exist: {val}'.format(val=_val))
#         logger.debug('initializer_dict: %r', initializer_dict)
#         try:
#             cpr = None
#             if patch:
#                 try:
#                     cpr = CherryPickRequest2.objects.get(**id_kwargs)
#                 except ObjectDoesNotExist:
#                     raise Http404(
#                         'Cherry Pick Request does not exist for: %r', id_kwargs)
#             else:
#                 cpr = CherryPickRequest2()
#             
#             model_field_names = [
#                 x.name for x in cpr._meta.get_fields()]
#             for key, val in initializer_dict.items():
#                 if key in model_field_names:
#                     setattr(cpr, key, val)
# 
#             cpr.save()
#             
#             
#             if 'cherry_pick_wells' in deserialized:
#                 
#                 cherry_pick_wells = self.find_wells(deserialized['cherry_pick_wells'])
#                 logger.info('found cp wells: %r', cherry_pick_wells)
#                 for well in cherry_pick_wells:
#                     cp_well = CherryPickWell.objects.create(
#                         cherry_pick_request=cpr,
#                         source_well=well)
#             
#             
#             return cpr
#         except Exception, e:
#             logger.exception('on patch_obj')
#             raise e
#     
#     def find_wells(self, cherry_pick_well_patterns ):
#         
#         if not isinstance(cherry_pick_well_patterns, (list,tuple)):
#             raise ValidationError(
#                 key='cherry_pick_wells',
#                 msg='must be a list of well_ids')
#         
#         wells = []
#         
#         errors = []
#         for _pattern in cherry_pick_well_patterns:
#             
#             if not WELL_ID_PATTERN.match(_pattern):
#                 errors.append(
#                     'pattern: %r is not a recognized well pattern (%r)' 
#                     % (_pattern, WELL_ID_PATTERN))
#             # TODO: search for PLATE well_id1, well_id2, ... well_idN
#             
#             wells.append(Well.objects.get(well_id=_pattern))
#         if errors:
#             raise ValidationError(
#                 key='cherry_pick_wells',
#                 msg='multiple: %r' % errors )    
#         
#         return wells

class CherryPickRequestResource(DbApiResource):        
    LCP_COPYWELL_KEY = '{library_short_name}/{source_copy_name}/{source_well_id}'

    class Meta:
    
        queryset = CherryPickRequest.objects.all().order_by('well_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'cherrypickrequest'
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
                    'get_lab_cherry_pick_plating'), 
                    name="api_get_lab_cherry_pick_plating"),
            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/lab_cherry_pick_plating/schema%s$"
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view(
                    'get_lab_cherry_pick_plating_schema'), 
                    name="get_get_lab_cherry_pick_plating_schema"),
        ]
    
#     def make_log(self, request, **kwargs):
#         log = DbApiResource.make_log(self, request, **kwargs)
#         log.ref_resource_name = self._meta.resource_name
#         log.key = kwargs.get('cherry_pick_request_id', None)
#         cpr = CherryPickRequest.objects.get(pk=log.key)
#         log.uri = '/'.join([
#             'screen',cpr.screen.facility_id,log.ref_resource_name,log.key])
#         logger.info('make cpr log: %r', log)
#         return log;

    @read_authorization
    def get_plate_mapping_file(self, request, **kwargs):
        
        request_method = request.method.lower()
        if request_method != 'get':
            raise BadRequest('Only GET is allowed')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        schema = super(CherryPickRequestResource, self).build_schema()

        cpr_id = param_hash['cherry_pick_request_id']
        logger.info('dispatch_plate_mapping_file: %r', cpr_id) 
         
        cpr_obj = CherryPickRequest.objects.get(
            cherry_pick_request_id=cpr_id)
        cpr = self._get_detail_response_internal(**kwargs)
        
        cpps = self.get_cherrypickplate_resource()._get_list_response_internal(
            **{
                'cherry_pick_request_id': cpr_id
            })
        cpps = {cpp['plate_ordinal']:cpp for cpp in cpps }
        lcp_schema = self.get_labcherrypick_resource().build_schema(
            cpr_obj.screen.screen_type)
        columns_to_include = [
            'library_plate','source_copy_name','source_well_name',
            'source_plate_type','destination_well','destination_plate_type']
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
        lcps = self.get_labcherrypick_resource()._get_list_response_internal(
            **{
                'cherry_pick_request_id': cpr_id,
                'destination_well__is_null': False,
                'includes': ['*', '-structure_image','-molfile'],
                'order_by': ['destination_well']
            })
        plate_map = defaultdict(list)
        location_map = {}
        for lcp in lcps:
            plate_map[lcp['cherry_pick_plate_number']].append(lcp)
            plate_copy = '{library_plate}/{source_copy_name}'.format(**lcp)
            location_map[plate_copy] = [
                lcp['library_plate'],lcp['source_copy_name'],lcp['location']]
        locations = [location_map[x] for x in sorted(location_map.keys())]
        vocabularies = \
            DbApiResource.get_vocabularies(lcp_schema['fields'])
        def vocab_function(key,val):
            if key in vocabularies:
                return vocabularies[key][val]['title']
            else:
                return val
        def title_function(key):
            return lcp_schema['fields'][key]['title']
        # Outputs:
        # Zip file:
        # - assay plate files as: Jen Smith (729) CP44451  Plate 01 of 1 (Run1)
        # - plate locations file
        # - README.txt
        
        readme_text = [
            'This zip file contains plate mappings for Cherry Pick Request %r' 
            % str(cpr_id) ]
        readme_text.append('Cherry pick plates:')
        plating_text = '{{filename}} Plated ({plating_date} by {plated_by_name})'
        # Open with delete=False; file will be closed and deleted 
        # by the FileWrapper1 instance when streaming is finished.
        with  NamedTemporaryFile(delete=False) as temp_file:
            with ZipFile(temp_file, 'w') as zipfile:
                    
                # Plate maps
                for cherry_pick_plate_number,lcps in plate_map.items():
                    raw_data = cStringIO.StringIO()
                    writer = unicodecsv.writer(raw_data, encoding='utf-8')
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
                    logger.info('lcp: %r',lcps[0])
                    logger.info('cpr: %r', cpr)
                    name = ((
                        'CP{{cherry_pick_request_id}}'
                        '_Screen{{screen_facility_id}}'
                        '_plate_{cherry_pick_plate_number}'
                        '_of_{{number_plates}}_{{requested_by_username}}.csv'
                        ).format(**lcps[0])
                        .format(**cpr)
                    )
                    logger.info('write; %r', name)
                    _text = plating_text.format(**cpps[cherry_pick_plate_number])
                    _text = _text.format(filename=name)
                    readme_text.append(_text)
                    zipfile.writestr(name, raw_data.getvalue())
                
                # Location map
                raw_data = cStringIO.StringIO()
                writer = unicodecsv.writer(raw_data)
                writer.writerow(['Source Plate','Source Copy','Location'])
                for location in locations:
                    writer.writerow(location)
                zipfile.writestr('plate-copy-location.csv', raw_data.getvalue())
                
                zipfile.writestr('readme.txt', '\n'.join(readme_text))
            logger.info('wrote file %r', temp_file)
        
            temp_file.seek(0, os.SEEK_END)
            size = temp_file.tell()
            temp_file.seek(0)   
        
        filename = 'CPR%d_plating_file.zip' % cpr['cherry_pick_request_id']
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
    def get_lab_cherry_pick_plating(self, request, **kwargs):
        return self.get_labcherrypick_resource().get_lab_cherry_pick_plating(request, **kwargs)    

    @read_authorization
    def get_lab_cherry_pick_plating_schema(self, request, **kwargs):
        return self.get_labcherrypick_resource().get_lab_cherry_pick_plating_schema(request, **kwargs)    

    def dispatch_cherry_pick_plate_view(self, request, **kwargs):
        return self.get_cherrypickplate_resource().dispatch('list', request, **kwargs)    

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
            # schema = super(CherryPickRequestResource, self).build_schema()
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
            filename = self._get_filename(readable_filter_hash)
                                  
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
            
            custom_columns = {
                'requested_by_username': _requestors.c.username,
                'requested_by_name': _requestors.c.name,
                'volume_approved_by_username': _approvers.c.username,
                'volume_approved_by_name': _approvers.c.name,
                'lab_name': lab_head_table.c.lab_name_full,
                'lab_head_username': lab_head_table.c.username,
                'lead_screener_name': _lead_screeners.c.name,
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
                    select([func.count(None)])
                        .select_from(_lcp)
                        .where(_lcp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id)
                        .where(_lcp.c.copy_id==None)),
                'total_number_lcps': (
                    select([func.count(None)])
                        .select_from(_lcp)
                        .where(_lcp.c.cherry_pick_request_id
                            ==_cpr.c.cherry_pick_request_id)),
                'last_plating_activity_date': (
                    select([func.max(_cpp.c.plating_date)])
                    .select_from(_cpp)
                    .where(_cpp.c.cherry_pick_request_id==_cpr.c.cherry_pick_request_id)),
                'last_screening_activity_date': (
                    select([func.max(_cpp.c.screening_date)])
                    .select_from(_cpp)
                    .where(_cpp.c.cherry_pick_request_id==_cpr.c.cherry_pick_request_id)),
                
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
                                ==literal_column('cherry_pick_request.cherry_pick_request_id'))
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
                            .join(_library, _well.c.library_id==_library.c.library_id))
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
            
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
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
        return DbApiResource.post_detail(self, request, **kwargs)

    @write_authorization
    @transaction.atomic      
    @un_cache  
    def patch_detail(self, request, **kwargs):
        '''
        Override to generate informational summary for callee
        '''
        deserialized = kwargs.pop('data', None)
        # allow for internal data to be passed
        if deserialized is None:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        logger.debug('patch detail %s, %s', deserialized,kwargs)

        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        schema = kwargs['schema']
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(deserialized,validate=True,**kwargs)

        original_data = None
        if kwargs_for_log:
            try:
                original_data = self._get_detail_response(request,**kwargs_for_log)
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
        new_data = self._get_detail_response(request,**kwargs_for_log)
        logger.debug('original: %r, new: %r', original_data, new_data)
        log = self.log_patch(
            request, original_data,new_data,log=log,
            excludes=['screener_cherry_picks'],
            **kwargs_for_log)
        # Set the log URI using the containing screen URI
        if log:
            log.uri = '/'.join([
                'screen',new_data['screen_facility_id'],log.ref_resource_name,log.key])
            log.save()
            logger.debug('log info: %r, %r', log, log.diffs )
        
        meta = {}
        if API_RESULT_META in patch_response:
            meta = patch_response[API_RESULT_META]
        return self.build_response(
            request,  { API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)
            
    
    @write_authorization
    @transaction.atomic
    @un_cache
    def patch_obj(self, request, deserialized, **kwargs):
        

        schema = kwargs['schema']
        fields = schema['fields']
        initializer_dict = {}

        id_kwargs = self.get_id(deserialized, **kwargs)
        
        logger.info('patch CPR: %r', deserialized)
        
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        logger.debug('param_hash: %r', param_hash)
        
        patch = bool(id_kwargs)
        initializer_dict = self.parse(deserialized, create=not patch)
        errors = self.validate(initializer_dict, patch=patch)
        if errors:
            raise ValidationError(errors)

        cpr = None
        if patch is True:
            try:
                cpr = CherryPickRequest.objects.get(**id_kwargs)
            except ObjectDoesNotExist:
                raise Http404(
                    'Cherry Pick Request does not exist for: %r', id_kwargs)
        
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
        else:
            screen = cpr.screen
            
        _key = 'requested_by_username'
        _val = deserialized.get(_key, None)
        if _val:
            try:
                requested_by_user = ScreensaverUser.objects.get(username=_val)
                if requested_by_user not in screen.get_screen_users():
                    raise ValidationError(
                        key='requested_by_username',
                        msg='"%s" must be one of the screen users: %r'
                            %(_val, [x.username for x in screen.get_screen_users()]))
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
                        self._meta.resource_name,volume_approved_by_user.user.user,'write'):
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
                if len(raw_screener_cps) == 0:
                    logger.info('removing screener_cherry_picks')
                    # FIXME: remove lab_cherry_picks, iif:
                    # - not plated
                    # - replace allocate well volumes
                    # - create a log for the action
                    _meta[API_MSG_SCPS_DELETED] = cpr.screener_cherry_picks.all().count()
                    cpr.screener_cherry_picks.all().delete()
                    cpr.save()
                else:
                    (cherry_pick_wells, warn_msg) = self.find_wells(raw_screener_cps)
                    final_warn_msg.extend(warn_msg)
                    logger.info(
                        'found screener_cherry_picks: %d', len(cherry_pick_wells))
                    not_allowed_libraries = set()
                    for well in cherry_pick_wells:
                        if well.library.screening_status != 'allowed':
                            not_allowed_libraries.add(
                                (well.library.short_name, well.library.screening_status))
                        if len(not_allowed_libraries)>0 and override_param is not True:
                            continue
                        screener_cherry_pick = ScreenerCherryPick.objects.create(
                            cherry_pick_request=cpr,
                            screened_well=well,
                            searched_well=well,
                            selected=True)
                    if len(not_allowed_libraries)>0 and override_param is not True:
                        raise ValidationError({
                            API_PARAM_OVERRIDE: 'required',
                            'screener_cherry_picks': (
                                'Override required to screen libraries that are '
                                'not allowed: %r' % not_allowed_libraries)
                            }
                        )
                    else:
                        final_warn_msg.append(
                            'Override used for libraries: %r' % not_allowed_libraries)
                    _meta[API_MSG_SCPS_CREATED] = cpr.screener_cherry_picks.all().count()    

            if final_warn_msg:
                _meta[API_MSG_WARNING] = final_warn_msg
                            
            response = { API_RESULT_OBJ: cpr, API_RESULT_META: _meta }
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
    def find_wells(self, cherry_pick_well_patterns ):
        logger.info('find wells for patterns: %r', cherry_pick_well_patterns)
        if not isinstance(cherry_pick_well_patterns, (list,tuple)):
            cherry_pick_well_patterns = (cherry_pick_well_patterns,)

        PLATE_PATTERN = re.compile(r'^(\d{1,5})$')
        wells = []
        errors = []
        warnings = []
        
        # Process the patterns by line
        for _line in cherry_pick_well_patterns:
            _line = _line.strip()
            if not _line:
                continue
            logger.info('pattern: %r', _line)
            plate_number = None
            
            for _pattern in re.split(r'[\s,]+', _line):
                
                match = WELL_ID_PATTERN.match(_pattern)
                if match is not None:
                    plate_number = int(match.group(1))
                    well_name = match.group(2).upper()
                    try:
                        well = Well.objects.get(
                            plate_number=plate_number, well_name=well_name)
                        wells.append(well)
                        logger.debug('found %r for %r', well, _pattern)
                    except ObjectDoesNotExist:
                        plate_number = None
                        errors.append(
                            'well: %r is does not exist' % _pattern) 
                elif PLATE_PATTERN.match(_pattern):
                    try:
                        plate = Plate.objects.get(plate_number=_pattern)
                        plate_number = int(_pattern)    
                        logger.debug('found %r for %r', plate, _pattern)
                    except ObjectDoesNotExist:
                        errors.append(
                            'plate: %r is does not exist' % _pattern) 
                elif plate_number is not None and WELL_NAME_PATTERN.match(_pattern):
                    try:
                        well = Well.objects.get(
                            plate_number=plate_number, well_name=_pattern.upper())
                        wells.append(well)
                        logger.debug('found %r for %r', well, _pattern)
                    except ObjectDoesNotExist:
                        errors.append(
                            'well: "%d:%s" is does not exist' 
                            % (plate_number,_pattern))
                else:
                    errors.append(
                        'pattern: %r is not a recognized as a well id, or '
                        'a plate number followed by a well name' % _pattern)
            logger.info('wells: %d', len(wells))
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
                key='screener_cherry_picks',
                msg='Cannot screen empty wells: '
                'source wells: %r' 
                    % [well.well_id for well in non_experimental_wells])
        return (wells, warnings)

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
    @transaction.atomic
    @un_cache
    def dispatch_reserve_map_lab_cherry_picks(self, request, **kwargs):
        
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
        
        # Find the fulfillable lab cherry picks
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
                        'includes': ['*', '-structure_image','-molfile'],
                    })}
        # Check for insufficient well volumes
        unfulfillable_wells = []
        override_well_volume = parse_val(
            param_hash.get(API_PARAM_VOLUME_OVERRIDE, False),
            API_PARAM_VOLUME_OVERRIDE, 'boolean')
        
#         for lcp in fulfillable_lcps:
        for copywell_id,lcp_data in lab_cherry_pick_copywells.items():
#             copywell_name = self.LCP_COPYWELL_KEY.format(**lcp_cw)
#             lcp_cw = lab_cherry_pick_copywells[lcp.source_well_id]
            lcp_scwv = Decimal(lcp_data['source_copy_well_volume'])
            logger.debug('consider lcp_cw: %r, %r to %r', 
                copywell_id, lcp_scwv, cpr.transfer_volume_per_well_approved)
            if ( lcp_scwv < cpr.transfer_volume_per_well_approved ):
                unfulfillable_wells.append(
                    (copywell_id, lcp_data['source_copy_well_volume'],
                        cpr.transfer_volume_per_well_approved) )
        if unfulfillable_wells and override_well_volume is not True:
            raise ValidationError({
                'transfer_volume_per_well_approved':
                    '%s: %r, "%s" required' 
                        % (API_MSG_LCPS_INSUFFICIENT_VOLUME, 
                           unfulfillable_wells, API_PARAM_VOLUME_OVERRIDE),
                API_MSG_LCPS_INSUFFICIENT_VOLUME: unfulfillable_wells,
                API_PARAM_VOLUME_OVERRIDE: 'required',
                
            })
        elif unfulfillable_wells:
            logger.info(
                'Override: unfulfillable wells: %r' % unfulfillable_wells)
        else:
            logger.info('all wells are fulfillable')

        # Reserve copy well volumes
        logger.info('Reserve copy well volumes...')
        copywell_reservation_meta = \
            self.get_copywell_resource().reserve_cherry_pick_volumes(
                cpr, fulfillable_lcps, parent_log)
        
        # Create assay plate mapping
        logger.info('Create the assay_plates...')  
        
        # cpr.is_randomized_assay_plate_layout
        # - if true, randomize the plate layout wells
        if cpr.is_randomized_assay_plate_layout:
            logger.info('randomize the available assay plate wells: %r',
                available_assay_plate_wells)
            random.shuffle(available_assay_plate_wells)
            logger.info('randomized available assay plate wells: %r',
                available_assay_plate_wells)
                  
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
        
        # Assignment ordering: plate, copy, source_well
#         lcps_by_plate_copy_well = [
#             (lcp['library_plate'],lcp['source_copy_name'],
#                 lcp['source_well_id'], lcp) 
#                 for lcp in fulfillable_lcps.all()]
#         lcps_by_plate_copy_well = sorted(
#             lcps_by_plate_copy_well_dest_well)
        #######
        
#         lcp_by_copy_plate = {}
#         for lcp in fulfillable_lcps.all():
#             copy_plate = '%s:%s' % (lcp.copy.name, lcp.source_well.plate_number)
#             copy_plate_list = lcp_by_copy_plate.get(copy_plate, [])
#             copy_plate_list.append(lcp)
#             lcp_by_copy_plate[copy_plate] = copy_plate_list
        lcp_by_plate_copy = defaultdict(list)
        for lcp in fulfillable_lcps.all():
            plate_copy = '%s:%s' % (lcp.source_well.plate_number,lcp.copy.name)
            lcp_by_plate_copy[plate_copy].append(lcp)
        # Sort internally
        for plate_copy in lcp_by_plate_copy.keys():
            lcps = lcp_by_plate_copy[plate_copy]
            lcps = sorted(lcps, key=lambda lcp: lcp.source_well_id)
            lcp_by_plate_copy[plate_copy] = lcps
            
        if cpr.keep_source_plate_cherry_picks_together is True:

            logger.info("use the bin packer to fit the lcp's to bins...")
            capacity = len(available_assay_plate_wells)
            packages = [ { 'name': plate_copy, 'size': len(lcps) } 
                for plate_copy,lcps in lcp_by_plate_copy.items() ]
            packed_bins = bin_packer.pack_bins(capacity, packages)
            packed_bins = sorted(packed_bins,key=lambda bin: bin[0]['name'])
            logger.info('packed bins: %r', packed_bins)
            
            logger.info(
                "assign the packed_bins to assay_plates, keep split bins adjacent...")
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
                    
                    logger.info('lcps: %r, assay_plate_well_index: %r', 
                        len(lcps_to_plate),assay_plate_well_index)
                
                    for lcp in lcps_to_plate:
                        lcp.cherry_pick_assay_plate = assay_plate
                        well_name = available_assay_plate_wells[assay_plate_well_index]
                        lcp.assay_plate_row = lims_utils.well_name_row_index(well_name)
                        lcp.assay_plate_column = lims_utils.well_name_col_index(well_name)
                        assay_plate_well_index += 1
                        lcp.save()
        else: 
            # keep_source_plate_cherry_picks_together == False
            # no bin packing
            capacity = len(available_assay_plate_wells)
            assay_plates_created = []
            assay_plate = None
            
            plate_copies = sorted(lcp_by_plate_copy.keys())
            for plate_copy in plate_copies:
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
    
                    assay_plate_well_index += 1
                    if assay_plate_well_index >= capacity:
                        assay_plate_well_index = 0
                        assay_plate = None
                    
            
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
#         for copy_plate, lcps in lcp_by_copy_plate.items():
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
        
        _meta = {
            API_MSG_LCP_PLATES_ASSIGNED: [(plate_copy_name, len(lcps)) 
                for plate_copy_name,lcps in lcp_by_plate_copy.items() ],
            API_MSG_LCP_ASSAY_PLATES_CREATED: cpap_assignments
        }
        _meta.update(copywell_reservation_meta[API_RESULT_META])
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

#     @write_authorization
#     @transaction.atomic
#     @un_cache
#     def bak_dispatch_reserve_map_lab_cherry_picks(self, request, **kwargs):
#         
#         request_method = request.method.lower()
#         if request_method != 'post':
#             raise BadRequest('Only POST is allowed')
#         convert_post_to_put(request)
# 
#         param_hash = self._convert_request_to_dict(request)
#         param_hash.update(kwargs)
#         cherry_pick_request_id = param_hash['cherry_pick_request_id']
#         logger.info(
#             'dispatch_reserve_map_lab_cherry_picks for: %r...', 
#             cherry_pick_request_id)
#         cpr = CherryPickRequest.objects.get(
#             cherry_pick_request_id=cherry_pick_request_id)
#         if not cpr.lab_cherry_picks.filter(copy__isnull=False).exists():
#             raise ValidationError(
#                 key='total_number_lcps', 
#                 msg='No (fulfilled) Lab Cherry Picks have been created')
#         
#         self.validate_cpr_for_plating(cpr) 
#         
#         parent_log = self.make_log(request)
#         parent_log.key = str(cpr.cherry_pick_request_id)
#         parent_log.uri = '/'.join([
#             'screen',cpr.screen.facility_id, 
#             parent_log.ref_resource_name,parent_log.key])        
#         parent_log.save()
#         
#         previous_date_reserved = cpr.date_volume_reserved
#         cpr.date_volume_reserved = _now().date() 
#         cpr.save()
#         parent_log.diffs = {
#             'date_volume_reserved': [previous_date_reserved, cpr.date_volume_reserved ]}
#         
#         status_messages = []
#         copywell_deallocation_meta = None
#         previous_number_of_plates = None           
#         
#         logger.info('check for previous plating assignments...')
#         allocated_lcps_query = cpr.lab_cherry_picks.filter(
#             cherry_pick_assay_plate__isnull=False)
#         if allocated_lcps_query.exists():
#             plated_assay_plates_query = \
#                 cpr.cherry_pick_assay_plates.filter(
#                     plating_date__isnull=False)                
#             if plated_assay_plates_query.exists():
#                 raise ValidationError({
#                     API_MSG_NOT_ALLOWED: API_MSG_CPR_PLATED_CANCEL_DISALLOWED, 
#                     API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count(),
#                     API_MSG_LCP_ASSAY_PLATES_PLATED: plated_assay_plates_query.count()
#             })
#             raise ValidationError({
#                 API_MSG_NOT_ALLOWED: 
#                     ('Lab cherry pick plates already assigned; '
#                     'reservation must be canceled before reassignment is allowed'),
#                 API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count()
#             })
#         
#         # Find the fulfillable lab cherry picks
#         available_assay_plate_wells = cpr.assay_plate_available_wells
#         logger.debug('available_assay_plate_wells: %r',available_assay_plate_wells)
#         logger.info('available_assay_plate_wells len: %d', 
#             len(available_assay_plate_wells))
# 
#         fulfillable_lcps = (
#             cpr.lab_cherry_picks.filter(copy__isnull=False)
#                 .order_by('source_well__plate_number','copy__name') )
# 
#         logger.info('re-Check and reserve copy volumes...')
#         lab_cherry_pick_copywells = \
#             { lcp['source_well_id']: lcp for lcp in
#                 self.get_labcherrypick_resource()._get_list_response_internal(
#                     **{
#                         'cherry_pick_request_id': cpr.cherry_pick_request_id,
#                         'source_copy_name__is_null': False,
#                         'includes': ['*', '-structure_image','-molfile'],
#                     })}
#         # Check for insufficient well volumes
#         unfulfillable_wells = []
#         override_well_volume = parse_val(
#             param_hash.get(API_PARAM_VOLUME_OVERRIDE, False),
#             API_PARAM_VOLUME_OVERRIDE, 'boolean')
#         
#         for lcp in fulfillable_lcps:
#             lcp_cw = lab_cherry_pick_copywells[lcp.source_well_id]
#             copywell_name = self.LCP_COPYWELL_KEY.format(**lcp_cw)
#             lcp_scwv = Decimal(lcp_cw['source_copy_well_volume'])
#             logger.debug('consider lcp_cw: %r, %r to %r', 
#                 copywell_name, lcp_scwv, cpr.transfer_volume_per_well_approved)
#             if ( lcp_scwv < cpr.transfer_volume_per_well_approved ):
#                 unfulfillable_wells.append(
#                     (copywell_name, lcp_cw['source_copy_well_volume'],
#                         cpr.transfer_volume_per_well_approved) )
#         if unfulfillable_wells and override_well_volume is not True:
#             raise ValidationError({
#                 'transfer_volume_per_well_approved':
#                     '%s: %r, "%s" required' 
#                         % (API_MSG_LCPS_INSUFFICIENT_VOLUME, 
#                            unfulfillable_wells, API_PARAM_VOLUME_OVERRIDE),
#                 API_MSG_LCPS_INSUFFICIENT_VOLUME: unfulfillable_wells,
#                 API_PARAM_VOLUME_OVERRIDE: 'required',
#                 
#             })
#         elif unfulfillable_wells:
#             logger.info(
#                 'Override: unfulfillable wells: %r' % unfulfillable_wells)
#         else:
#             logger.info('all wells are fulfillable')
# 
#         # Reserve copy well volumes
#         logger.info('Reserve copy well volumes...')
#         copywell_reservation_meta = \
#             self.get_copywell_resource().reserve_cherry_pick_volumes(
#                 cpr, fulfillable_lcps, parent_log)
#         
#         # Create assay plate mapping
#         logger.info('Create the assay_plates...')  
#         
#         # cpr.is_randomized_assay_plate_layout
#         # - if true, randomize the plate layout wells
#         if cpr.is_randomized_assay_plate_layout:
#             logger.info('randomize the available assay plate wells: %r',
#                 available_assay_plate_wells)
#             random.shuffle(available_assay_plate_wells)
#             logger.info('randomized available assay plate wells: %r',
#                 available_assay_plate_wells)
#                   
#         next_plate_ordinal = [1]
#         def create_next_assay_plate():
#             assay_plate = CherryPickAssayPlate.objects.create(
#                 cherry_pick_request=cpr,
#                 plate_ordinal=next_plate_ordinal[0],
#                 # TODO: deprecate attempt_ordinal
#                 attempt_ordinal=0,
#                 assay_plate_type=cpr.assay_plate_type,
#                 # TODO: deprecate cpapt
#                 cherry_pick_assay_plate_type='CherryPickAssayPlate',
#             )
#             next_plate_ordinal[0] = next_plate_ordinal[0]+1
#             return assay_plate
#         
#         # Assignment ordering: plate, copy, source_well
#         lcps_by_plate_copy_well = [
#             (lcp['library_plate'],lcp['source_copy_name'],
#                 lcp['source_well_id'], lcp) 
#                 for lcp in fulfillable_lcps.all()]
#         lcps_by_plate_copy_well = sorted(
#             lcps_by_plate_copy_well_dest_well)
#         
#         lcp_by_copy_plate = {}
#         for lcp in fulfillable_lcps.all():
#             copy_plate = '%s:%s' % (lcp.copy.name, lcp.source_well.plate_number)
#             copy_plate_list = lcp_by_copy_plate.get(copy_plate, [])
#             copy_plate_list.append(lcp)
#             lcp_by_copy_plate[copy_plate] = copy_plate_list
#         
#         if cpr.keep_source_plate_cherry_picks_together is True:
# 
#             logger.info("use the bin packer to fit the lcp's to bins...")
#             capacity = len(available_assay_plate_wells)
#             packages = [ { 'name': copy_plate, 'size': len(lcps) } 
#                 for copy_plate,lcps in lcp_by_copy_plate.items() ]
#             
#             packed_bins = bin_packer.pack_bins(capacity, packages)
#             logger.info('packed bins: %r', packed_bins)
#             
#             logger.info(
#                 "assign the packed_bins to assay_plates, keeping split bins adjacent...")
#             ordered_bins = []
#             for packed_bin in packed_bins:
#                 if packed_bin not in ordered_bins:
#                     ordered_bins.append(packed_bin)
#                     partially_packed_copy_plate = [
#                         copy_plate for package in packed_bin
#                             if len(lcp_by_copy_plate[package['name']]) > package['size'] ]
#                     if len(partially_packed_copy_plate) > 1:
#                         raise ProgrammingError(
#                             'Packed bins should not contain more than one '
#                             'partially packed copy_plates: %r' 
#                             % partially_packed_copy_plate)
#                     # NOTE: assume that source plate will never need > 2 assay plates
#                     if partially_packed_copy_plate:
#                         copy_plate = partially_packed_copy_plate[0]
#                         for second_bin in packed_bins:
#                             if second_bin != packed_bin:
#                                 if copy_plate in [ p['name'] for p in second_bin]:
#                                     ordered_bins.append(second_bin)
#             
#             logger.info('create assay plates, in order...')
#             assay_plates_created = []
#             for packed_bin in ordered_bins:
#                 logger.info('using packed bin: %r', packed_bin)
#                 assay_plate = create_next_assay_plate()
#                 assay_plates_created.append(assay_plate)
#                 assay_plate_well_index = 0
#                 for package in packed_bin:
#                     
#                     logger.info('package: %r', package)
#                     
#                     copy_plate = package['name']
#                     size = package['size']
#                     copy_plate_lcps = lcp_by_copy_plate[copy_plate]
#                     lcps = copy_plate_lcps[:size]
#                     lcp_by_copy_plate[copy_plate] = lcps
#                     
#                     logger.info('lcps: %r, assay_plate_well_index: %r', 
#                         len(lcps),assay_plate_well_index)
#                 
#                     for lcp in lcps:
#                         lcp.cherry_pick_assay_plate = assay_plate
#                         well_name = available_assay_plate_wells[assay_plate_well_index]
#                         lcp.assay_plate_row = lims_utils.well_name_row_index(well_name)
#                         lcp.assay_plate_column = lims_utils.well_name_col_index(well_name)
#                         assay_plate_well_index += 1
#                         lcp.save()
#         else: # keep_source_plate_cherry_picks_together == False
#             capacity = len(available_assay_plate_wells)
#             assay_plates_created = []
#             assay_plate = None
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
# 
#         # verify that all lcps have been assigned
# #         for copy_plate, lcps in lcp_by_copy_plate.items():
#         for lcp in fulfillable_lcps.all():
#             copy_plate = '%s:%s' % (lcp.copy.name, lcp.source_well.plate_number)
#             if lcp.cherry_pick_assay_plate is None:
#                 raise ProgrammingError(
#                     'lcp has not been assigned: %s, %r' 
#                     % lcp, copy_plate)
#         
#         # return some stats
#         cpap_assignments = [ 
#             'Plate ordinal: %d, Picks: %d' 
#                 % (ap.plate_ordinal, ap.labcherrypick_set.all().count()) 
#                 for ap in assay_plates_created ]
#         
#         _meta = {
#             API_MSG_LCP_PLATES_ASSIGNED: [(copy_name, len(lcps)) 
#                 for copy_name,lcps in lcp_by_copy_plate.items() ],
#             API_MSG_LCP_ASSAY_PLATES_CREATED: cpap_assignments
#         }
#         _meta.update(copywell_reservation_meta[API_RESULT_META])
#         if unfulfillable_wells:
#             _meta[API_MSG_LCPS_VOLUME_OVERRIDDEN] = unfulfillable_wells
# 
#         # log
#         parent_log.diffs.update({
#             'number_plates': [previous_number_of_plates, len(assay_plates_created)]
#             })
#         parent_log.json_field = _meta
#         parent_log.save()
#         
#         return self.build_response(
#             request, { API_RESULT_META: _meta }, 
#             response_class=HttpResponse, **kwargs)

    
    @write_authorization
    @transaction.atomic
    @un_cache
    def dispatch_delete_lab_cherry_picks(self, request, **kwargs):
        
        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')
        convert_post_to_put(request)

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        cpr_id = param_hash['cherry_pick_request_id']
        logger.info('dispatch_cancel_reservation: %r', cpr_id) 
        schema = super(CherryPickRequestResource, self).build_schema()
         
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
                'cherry_pick_request_id': cpr_id })
            parent_log = self.make_log(request)
            parent_log.key = str(cpr_id)
            parent_log.uri = '/'.join([
                'screen',cpr.screen.facility_id,
                parent_log.ref_resource_name,parent_log.key])        
            parent_log.comment = API_MSG_LCPS_REMOVED
            
            meta = {}
            meta[API_MSG_LCPS_REMOVED] = cpr.lab_cherry_picks.all().count()
            cpr.lab_cherry_picks.all().delete()
    
            new_cpr = self._get_detail_response_internal(**{
                'cherry_pick_request_id': cpr_id })
            parent_log = self.log_patch(
                request, original_cpr, new_cpr, parent_log, 
                excludes=['screener_cherry_picks'])
            parent_log.save()
            
            return self.build_response(
                request,  {API_RESULT_META: meta }, 
                response_class=HttpResponse, **kwargs)
    
        # Empty response if no action taken
        return HttpRequest

    @write_authorization
    @transaction.atomic
    @un_cache
    def dispatch_cancel_reservation(self, request, **kwargs):
        
        request_method = request.method.lower()
        if request_method != 'post':
            raise BadRequest('Only POST is allowed')
        convert_post_to_put(request)

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        cpr_id = param_hash['cherry_pick_request_id']
        logger.info('dispatch_cancel_reservation: %r', cpr_id) 
        schema = super(CherryPickRequestResource, self).build_schema()
         
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
    
        # Empty response if no action taken
        return HttpRequest
    
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
    @transaction.atomic
    @un_cache
    def dispatch_set_duplex_lab_cherry_picks(self, request, **kwargs):
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
         
        schema = super(CherryPickRequestResource, self).build_schema()
         
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
        self.get_labcherrypick_resource().clear_cache()
        meta[API_MSG_LCPS_CREATED] = len(lcps_created)

        # Find and assign copies that are fulfillable
        meta_action = self._find_copies(cpr)
        if meta_action:
            meta.update(meta_action)
        
        new_cpr = self._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        self.log_patch(request, original_cpr, new_cpr, parent_log, 
            excludes=['screener_cherry_picks'])
        parent_log.json_field = meta
        parent_log.save()
        logger.info('set_duplex_lab_cherry_picks: log: %r, %r', parent_log, parent_log.diffs)
        
        return self.build_response(
            request,  {API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)
    
    @write_authorization
    @transaction.atomic
    @un_cache
    def dispatch_set_lab_cherry_picks(self, request, **kwargs):
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
         
        schema = super(CherryPickRequestResource, self).build_schema()
         
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
            
        lcps_created = []
        for scp in cpr.screener_cherry_picks.filter(selected=True):
           lab_cherry_pick = LabCherryPick.objects.create(
                cherry_pick_request=cpr,
                source_well=scp.screened_well,
                screener_cherry_pick=scp)
           lcps_created.append(lab_cherry_pick)
        self.get_labcherrypick_resource().clear_cache()
        meta[API_MSG_LCPS_CREATED] = len(lcps_created)

        # Find and assign copies that are fulfillable
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
        on the Cherry Pick Request
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
                    'includes': ['*', '-structure_image','-molfile'],
                })
        logger.info('found %d eligible copy-wells for %d lab cherry picks', 
            len(eligible_lab_cherry_pick_copywells), 
            cpr.lab_cherry_picks.all().count())
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
        
        logger.info('pick candidates by library: %r', pick_candidates_by_library)    
        lcp_assigned_count = 0
        for library_short_name in pick_candidates_by_library.keys():
            
            copy_set = copy_sets_by_library.get(library_short_name,set())
            well_picks_for_library = pick_candidates_by_library[library_short_name]
            minimal_copy_sets = self._get_minimal_copy_set(copy_set,well_picks_for_library)
        
            if minimal_copy_sets:
                chosen_minimal_copy_set = sorted(minimal_copy_sets)[0]
                logger.info('library: %r, chosen minimal copy set: %r', 
                    library_short_name, chosen_minimal_copy_set)
#                 pick_copy_sets = well_picks_for_library[library_short_name]
                
                for lcp in cpr.lab_cherry_picks.all():
                    if lcp.source_well_id in well_picks_for_library:
                        pick_copy_set = well_picks_for_library.get(lcp.source_well_id, None)
                        if pick_copy_set is None:
                            logger.info('no pick copy set for well: %r', lcp.source_well_id)
                            continue
                        
                        eligible_copies = pick_copy_set & chosen_minimal_copy_set
                        
                        if eligible_copies:
                            copy_full_name = sorted(eligible_copies)[0]
                            lcp.copy = copy_instance_cache[copy_full_name]
                            lcp.save()
                            lcp_assigned_count += 1
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
#             meta[API_MSG_RESULT] = meta
        return meta

    def _get_minimal_copy_set(self,copy_set, pick_copy_sets ):
        minimal_copy_sets = []
        logger.info('complete copy set: %d = %r', len(copy_set), copy_set)
        logger.info('Determine minimal copy set(s)...')
        def powerset(iterable):
            "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
            s = list(iterable)
            return set(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))
        minimal_matching_copyset_len = None
        for copy_subset in powerset(
                sorted(copy_set, 
                       key=lambda copy_set: '%d:%s' % (len(copy_set),copy_set))):
            if minimal_matching_copyset_len is not None:
                if len(pick_copy_set) > minimal_matching_copyset_len:
                    break
            logger.debug('considering set: %r', copy_subset)
            for pick_copy_set in pick_copy_sets.values():
                logger.debug('considering pick copy set: %r', pick_copy_set)
                if set(copy_subset) & pick_copy_set:
                    logger.debug('eligible minimum copy set: %r', pick_copy_set)
                    minimal_matching_copyset_len = len(pick_copy_set)
                    minimal_copy_sets.append(pick_copy_set) 
        
        logger.info('found %d minimal copy sets', len(minimal_copy_sets))
        
        return minimal_copy_sets
        

class ScreenerCherryPickResource(DbApiResource):        

    class Meta:

        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'screenercherrypick'
        
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
 
        schema = kwargs['schema']
        id_attribute = schema['id_attribute']

        original_cpr_data = self.get_cpr_resource()._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr.cherry_pick_request_id })
        
        original_data = self._get_list_response(request, **{
            'cherry_pick_request_id': kwargs['cherry_pick_request_id'],
            API_PARAM_SHOW_OTHER_REAGENTS: True })
        if not original_data: 
            raise BadRequest(
                'Cannot set Screener Cherry Picks using patch; '
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
                        'searched_well_id': searched_well_id }
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
            return HttpResponse
        
        logger.info('selection updates: %r', selection_updates.keys())
                    
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
                            'not currently selected: %r', screened_well_id)
                else:
                    messages.append(
                        'not unselecting, well not found in original_scps: %r', screened_well_id)
        
        if not scps_to_create and not scps_to_reselect and not scps_to_unselect:
            _data = { API_RESULT_META: messages }
            return self.build_response(
                request, _data, response_class=HttpResponse, **kwargs)
        
        if scps_to_create:
            (cherry_pick_wells, warn_msg) = \
                CherryPickRequestResource.find_wells(scps_to_create)
            messages.extend(warn_msg)
                
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
            API_MSG_SCP_CREATED: '%r' % scps_to_create,
            API_MSG_SCP_UNSELECTED: '%r' % scps_to_unselect,
            API_MSG_SCP_RESELECTED: '%r' % scps_to_reselect
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
            return self.build_response(request, self.build_schema(),**kwargs)
        
        cherry_pick_request_id = kwargs.pop('cherry_pick_request_id')
        try:
            cpr = CherryPickRequest.objects.get(
                cherry_pick_request_id=cherry_pick_request_id)
            # NOTE: does not support disinction between Small Molecule and Natural Product
            return self.build_response(
                request, self.build_schema(cpr.screen.screen_type),**kwargs)
            
        except ObjectDoesNotExist, e:
            raise Http404(
                'cannot build schema - CherryPickRequest ID needed'
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
            field['visibility'] = []
        
        # Overlay the original scp fields on the top
        schema['fields'].update(original_fields)
        logger.info('schema  fields: %r', schema['fields'].keys())
        return schema

    def build_sqlalchemy_columns(self, fields, base_query_tables=None, 
            custom_columns=None):
        sub_columns = self.get_sr_resource().build_sqlalchemy_columns(fields)
        sub_columns.update(self.get_smr_resource().build_sqlalchemy_columns(fields))
        sub_columns['plate_number'] = (literal_column(
            "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
            .label('plate_number'))
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
        try:
            
            # Note: build schema for each request to use the subtype
            schema = self.build_schema(library_classification=cpr.screen.screen_type)
#             filename = self._get_filename(schema, kwargs, filename=filename)
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
#             manual_field_includes.add('searched_library_plate')
            manual_field_includes.add('searched_well_id')
            manual_field_includes.add('selected')
            manual_field_includes.add('cherry_pick_request_id')
            
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, **extra_params)
                                  
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
            
            if show_other_reagents is True:
                _original_scps = (
                    select([
                        _scp.c.screener_cherry_pick_id,
                        _scp.c.screened_well_id,
                        _scp.c.searched_well_id,
#                         _well.c.plate_number.label('searched_library_plate'),
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
                    ).cte('original_scps')
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
                    ).cte('combined')

                working_scp = combined_scps
            else:
                working_scp = _scp
                
            custom_columns = {
                'screener_cherry_pick_id': working_scp.c.screener_cherry_pick_id,
#                 'selected': (
#                     select([_scp.c.selected]).select_from(_scp)
#                         .where(_scp.c.screener_cherry_pick_id
#                             ==working_scp.c.screener_cherry_pick_id)),
                'screened_well_id': working_scp.c.screened_well_id,
                'selected': working_scp.c.selected,
#                 'selected_on_server': working_scp.c.selected,
                'cherry_pick_request_id': working_scp.c.cherry_pick_request_id,
# Does not work for pool-duplex LCP's
#                 'status': case([
#                     (_lcp.c.cherry_pick_assay_plate_id!=None, 
#                         text("'plated'") ),
#                     (_lcp.c.lab_cherry_pick_id!=None, 
#                         text("'selected'") )],
#                     else_=text("''")),
                
            }
#             if show_other_reagents is True:
            custom_columns['searched_well_id'] = working_scp.c.searched_well_id
#             else:
#                 custom_columns['searched_well_id'] = working_scp.c.screened_well_id
                
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            
            j = join(working_scp,_cpr, 
                working_scp.c.cherry_pick_request_id == _cpr.c.cherry_pick_request_id)
            
            # Fixme: does not work for duplex lcp's
            # j = j.join(_lcp, working_scp.c.screener_cherry_pick_id 
            #         ==_lcp.c.screener_cherry_pick_id, isouter=True)
            
            
            j = j.join(_screen,
                _cpr.c.screen_id == _screen.c.screen_id)
            j = j.join(_well, _well.c.well_id==working_scp.c.screened_well_id)
            j = j.join(_library, _well.c.library_id==_library.c.library_id)
            j = j.join(_reagent, working_scp.c.screened_well_id==_reagent.c.well_id,
                isouter=True)
            
            stmt = select(columns.values()).select_from(j)
            stmt = stmt.where(_cpr.c.cherry_pick_request_id==cherry_pick_request_id)
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
                if show_other_reagents is True:
                    
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

    class Meta:

        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'labcherrypick'
        
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

        if 'cherry_pick_request_id' not in kwargs:
            raise BadRequest('cherry_pick_request_id')
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=kwargs['cherry_pick_request_id'])
        logger.info(
            'patch_list: cpr: %r, screen: %r...', cpr, cpr.screen.facility_id)
         
        deserialized = self.deserialize(request)
        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
        
        schema = kwargs['schema']
        id_attribute = schema['id_attribute']
        COPYWELL_KEY = CherryPickRequestResource.LCP_COPYWELL_KEY
        
        self.get_cpr_resource().validate_cpr_for_plating(cpr) 

        if not cpr.lab_cherry_picks.all().exists():
            raise ValidationError(
                key='lcp_selection_updates', 
                msg='No Lab Cherry Picks have been created')
        current_lcps = { 
            lcp.source_well_id:lcp for lcp in cpr.lab_cherry_picks.all() }

        parent_log = self.make_log(request)
        parent_log.key = str(cpr.cherry_pick_request_id)
        parent_log.uri = '/'.join([
            'screen',cpr.screen.facility_id,
            parent_log.ref_resource_name,parent_log.key])        
        parent_log.save()

        if cpr.lab_cherry_picks.filter(cherry_pick_assay_plate__isnull=False).exists():
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
                    'cancel reservation to deallocate plates'),
                API_MSG_LCP_PLATES_ASSIGNED: cpr.cherry_pick_assay_plates.count()
            })
                
        # Validate selections are valid and allowed:
        required_for_patch = set([
            'source_well_id', 'selected'])
        copies = {}
        selection_updates = {}
        selections_per_well = defaultdict(list)
        for selection_update in deserialized:
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
                raise ValidationError(
                    key='source_well_id', 
                    msg='not a current lab cherry pick: %r' % source_well_id)
            
            # locate copies:
            copy = None
            source_copy_id = selection_update.get('source_copy_id', None)
            if source_copy_id:
                copy = Copy.objects.get(copy_id=source_copy_id)
            else:
                source_copy_name = selection_update.get('source_copy_name')
                plate_number = lims_utils.well_id_plate_number(source_well_id)
                plate = Plate.objects.get(
                    plate_number=plate_number,copy__name=source_copy_name)
                copy = plate.copy
            selection_update['copy'] = copy
            selection_update['source_copy_name'] = copy.name
            selection_update['library_short_name'] = copy.library.short_name
            selection_copy_name = COPYWELL_KEY.format(**selection_update)
            selection_updates[selection_copy_name] = selection_update
            
            if selection_update['selected'] is True:
                selections_per_well[source_well_id].append(copy.name)
        
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
        logger.debug('selection_updates: %r', selection_updates)
        
        changed = []
        deselected = []
        selected = []
        for selection_copy_name, selection_update in selection_updates.items():
            source_well_id = selection_update['source_well_id']
            current_lcp = current_lcps[source_well_id]
            current_copy_name = None
            if current_lcp.copy is not None:
                current_copy_name = COPYWELL_KEY.format(
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
                    current_lcp.save()
                else:
                    logger.info('lcp is already selected: %r, %r', 
                        current_copy_name, selection_update)
            else:
                if selection_copy_name == current_copy_name:
                    deselected.append(current_copy_name)
                    current_lcp.copy = None
                    current_lcp.save()
                else:
                    logger.info('lcp is already unselected: %r, %r', 
                        current_copy_name, selection_update)
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
                    'includes': '*'
                })
        for lcp_cw in new_lab_cherry_pick_copywells:
            name = COPYWELL_KEY.format(**lcp_cw)
            if ( Decimal(lcp_cw['source_copy_well_volume'])
                < cpr.transfer_volume_per_well_approved ):
                logger.info(
                    'vol requires override: lcp_cw: %r, approved: %r, available: %r', 
                    name, cpr.transfer_volume_per_well_approved, 
                    Decimal(lcp_cw['source_copy_well_volume']))
                unfulfillable_wells.append(COPYWELL_KEY.format(**lcp_cw))
        warning_messages = {}
        if unfulfillable_wells:
            warning_messages[API_MSG_LCPS_INSUFFICIENT_VOLUME] = unfulfillable_wells

        _data = {
            API_RESULT_META: {
                API_MSG_LCP_CHANGED: changed,
                API_MSG_LCP_DESELECTED: deselected,
                API_MSG_LCP_SELECTED: selected,
                API_MSG_WARNING: warning_messages,
                }}
        parent_log.json_field = _data
        parent_log.save()
        
        return self.build_response(
            request, _data, response_class=HttpResponse, **kwargs)

    def get_schema(self, request, **kwargs):
        ''' Generate an HttpResponse for the schema '''
        
        return self.build_response(request, self.build_schema(**kwargs))

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
                field['visibility'] = []
            
            # Overlay the original scp fields on the top
            schema['fields'].update(original_fields)
            logger.debug('schema  fields: %r', schema['fields'].keys())
            return schema
        except Exception, e:
            logger.exception('xxx: %r', e)
            raise

    @read_authorization
    def get_lab_cherry_pick_plating(self, request, **kwargs):
        
        if not 'cherry_pick_request_id' in kwargs:
            raise BadRequest('cherry_pick_request_id is required')
        schema = self.build_lab_cherry_pick_plating_schema(**kwargs)
        kwargs['schema'] = schema
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        kwargs['status__eq'] = 'plated'
        kwargs['order_by'] = ['cherry_pick_plate_number','destination_well']
        return self.build_list_response(request, **kwargs)

    @read_authorization
    def get_lab_cherry_pick_plating_schema(self, request, **kwargs):
        ''' Generate an HttpResponse for the plate_mapping specific schema '''
        
        if not 'cherry_pick_request_id' in kwargs:
            raise BadRequest('cherry_pick_request_id is required')
        
        return self.build_response(
            request, 
            self.build_lab_cherry_pick_plating_schema(**kwargs),
            **kwargs)

    def build_lab_cherry_pick_plating_schema(self, **kwargs):
        schema = self.build_schema(**kwargs)
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
#            'volume_approved_ul',
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
                field['ordinal'] = visible_fields.index(key)
            else:
                field['visibility'] = []
        fields['source_copy_well_volume']['title'] = \
            'Source CopyWell Volume (after transfer)'
        fields['source_copy_well_volume']['description'] = \
            'Source CopyWell Volume (after transfer of '\
            'cherry pick volume to the destination well)'
            
        logger.info('plate mapping visible fields: %r', 
            {k:{'key': v['key'], 'scope':v['scope'], 'title':v['title']}
                for k,v in fields.items() if 'l' in v['visibility'] })
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
        sub_columns['plate_number'] = (literal_column(
            "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
            .label('plate_number'))
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
        extra_params['CPR'] = cherry_pick_request_id
        cpr = CherryPickRequest.objects.get(
            cherry_pick_request_id=cherry_pick_request_id)
        show_copy_wells = parse_val(
            param_hash.get(API_PARAM_SHOW_COPY_WELLS, False),
            API_PARAM_SHOW_COPY_WELLS, 'boolean')
        logger.info('%r: %r', API_PARAM_SHOW_COPY_WELLS,show_copy_wells)
        if show_copy_wells is True:
            extra_params[API_PARAM_SHOW_COPY_WELLS] = None
        show_available_and_retired_copy_wells = parse_val(
            param_hash.get(API_PARAM_SHOW_RETIRED_COPY_WELlS, False),
            API_PARAM_SHOW_RETIRED_COPY_WELlS, 'boolean')
        if show_available_and_retired_copy_wells is True:
            extra_params[API_PARAM_SHOW_RETIRED_COPY_WELlS] = None
        logger.info('%r: %r', 
            API_PARAM_SHOW_RETIRED_COPY_WELlS,show_available_and_retired_copy_wells)
        show_unfulfilled = parse_val(
            param_hash.get(API_PARAM_SHOW_UNFULFILLED, False),
            API_PARAM_SHOW_UNFULFILLED, 'boolean')
        if show_unfulfilled is True:
            extra_params[API_PARAM_SHOW_UNFULFILLED] = None
        try:
            
            # Note: build schema for each request to use the subtype
            if schema is None: 
                schema = self.build_schema(library_classification=cpr.screen.screen_type)
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            manual_field_includes.add('source_copy_id')
            manual_field_includes.add('selected_copy_name')

            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash, **extra_params)
            
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
#                 'volume_approved_ul': _cpr.c.transfer_volume_per_well_approved * 1e6
            }
            # pull in the copy well volume
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
            
            if show_copy_wells is not True and show_available_and_retired_copy_wells is not True:
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

                custom_columns.update({
                    'selected_copy_name': case([
                        (_lcp.c.copy_id==_copy.c.copy_id, _copy.c.name )],
                            else_=None),
                        })
            else:
                _original_copy = _copy.alias('c1')
                _copyplates_available = (
                    select([
                        _lcp.c.lab_cherry_pick_id,
                        _copy.c.copy_id,
                        _p.c.plate_id ])
                    .select_from(
                        _lcp.join(_well,_well.c.well_id==_lcp.c.source_well_id)
                            .join(_p, _well.c.plate_number==_p.c.plate_number)
                            .join(_copy, _p.c.copy_id==_copy.c.copy_id))
                    .where(_lcp.c.cherry_pick_request_id==cherry_pick_request_id))
                if show_available_and_retired_copy_wells is not True:    
                    _copyplates_available = \
                        _copyplates_available.where(
                            _copy.c.usage_type=='cherry_pick_source_plates')
                else:
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
                
                custom_columns.update({
                    'selected_copy_name': case([
                        (_lcp.c.copy_id==_copy.c.copy_id, _copy.c.name )],
                            else_=_original_copy.c.name),
                    'status': case([
                        (and_(_lcp.c.copy_id==_copy.c.copy_id,
                              _lcp.c.cherry_pick_assay_plate_id==None,), 
                            text("'%s'"%VOCAB_LCP_STATUS_SELECTED) ),
                        (and_(_lcp.c.copy_id==_copy.c.copy_id,
                              _lcp.c.cherry_pick_assay_plate_id!=None,), 
                            text("'%s'"%VOCAB_LCP_STATUS_PLATED) )],
                        else_=text("'%s'"%VOCAB_LCP_STATUS_NOT_SELECTED)),
#                     'source_copy_well_volume': case([
#                         (_well.c.library_well_type=='experimental', 
#                              func.coalesce(_cw.c.volume, 
#                                  _p.c.remaining_well_volume, _p.c.well_volume) )],
#                         else_=None),
#                     'source_copy_well_initial_volume': case([
#                         (_well.c.library_well_type=='experimental', 
#                              func.coalesce(_cw.c.initial_volume,_p.c.well_volume) )],
#                          else_=None),
#                     'source_copy_well_consumed_volume': case([
#                         (_well.c.library_well_type=='experimental', 
#                             func.coalesce(_cw.c.initial_volume,_p.c.well_volume)-
#                                 func.coalesce(_cw.c.volume, _p.c.remaining_well_volume) )],
#                         else_=None),
                })

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)
            
            stmt = select(columns.values()).select_from(j)
            stmt = stmt.where(_cpr.c.cherry_pick_request_id==cherry_pick_request_id)

            if show_unfulfilled:
                if show_copy_wells is not True and show_available_and_retired_copy_wells is not True:
                    stmt = stmt.where(_copy.c.name==None)
                else:
                    stmt = stmt.where(_original_copy.c.name==None)
            
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            if not order_clauses:
                # Ordering for well_id must be alphanumeric
                order_clause = text(
                    "(substring({field_name}, '^[0-9]+'))::int asc " # cast to integer
                    ",substring({field_name}, ':(.*$)') asc  "  # works as text
                    .format(field_name='source_well_id'))
                # For string field ordering, double sort as numeric and text
                stmt = stmt.order_by(order_clause)
                if show_copy_wells is True or show_available_and_retired_copy_wells is True:
                    # For string field ordering, double sort as numeric and text
                    stmt = stmt.order_by(
                        desc('selected'),desc('source_copy_name'),)
           
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


class CherryPickPlateResource(DbApiResource):        

    class Meta:

        queryset = CherryPickAssayPlate.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'cherrypickassayplate'
        
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
            # schema = super(CherryPickPlateResource, self).build_schema()

        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id:
            param_hash['cherry_pick_request_id__eq'] = cherry_pick_request_id
        plate_ordinal = param_hash.pop('plate_ordinal', None)
        if plate_ordinal:
            param_hash['plate_ordinal__eq'] = plate_ordinal
#         attempt_ordinal = param_hash.pop('attempt_ordinal', None)
#         if attempt_ordinal:
#             param_hash['attempt_ordinal__eq'] = attempt_ordinal

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
  
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)                                  
            
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
                'plated_by_name': _plated_by.c.name,
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
                        .where(_apilog.c.ref_resource_name=='cherrypickassayplate')
                        .where(_apilog.c.key==
                            _concat(
                                cast(_cpap.c.cherry_pick_request_id,sqlalchemy.sql.sqltypes.Text),
                                '/',
                                cast(_cpap.c.plate_ordinal,sqlalchemy.sql.sqltypes.Text)))
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
        logger.info('patch lcps: %r', deserialized)
        
        schema = kwargs['schema']
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
        
        original_cpr = self.get_cherry_pick_resource()._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        original_cpap_data = self._get_list_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        
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
                if plated_by_username is None:
                    raise ValidationError(key='plated_by_username', msg='required')
                try:
                    user = ScreensaverUser.objects.get(username=plated_by_username)
                    
                    if user not in cpr.screen.get_screen_users():
                        raise ValidationError(
                            key='plated_by_username',
                            msg='Not a member of this screen')
                    
                    cpap.plated_by = user
                    cpap.plating_date = plating_date
                    cpap.save()
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key='plated_by_username', 
                        msg='User: "%s" not found' % plated_by_username)

            if screening_date is not None:
                if cpap.plating_date is None:
                    raise ValidationError(
                        key='screening_date',
                        msg='"plating_date" must be set before "screening_date" can be set')
                screened_by_username = initializer_dict.get('screened_by_username', None)
                if screened_by_username is None:
                    raise ValidationError(key='screened_by_username', msg='required')
                try:
                    user = ScreensaverUser.objects.get(username=screened_by_username)
                    
                    if user not in cpr.screen.get_screen_users():
                        raise ValidationError(
                            key='screened_by_username',
                            msg='Not a member of this screen')
                    
                    cpap.screened_by = user
                    cpap.screening_date = screening_date
                    cpap.save()
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key='screened_by_username', 
                        msg='User: "%s" not found' % screened_by_username)

        # Log
        self.get_cherry_pick_resource().clear_cache()
        new_cpr = self.get_cherry_pick_resource()._get_detail_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        new_cpap_data = self._get_list_response_internal(**{
            'cherry_pick_request_id': cpr_id })
        self.get_cherry_pick_resource().log_patch(
            request, original_cpr, new_cpr, log=parent_log)
        patch_logs = self.log_patches(
            request, original_cpap_data, new_cpap_data, schema=schema, parent_log=parent_log)
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
            request, { API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)
                
class LibraryCopyResource(DbApiResource):

    class Meta:

        queryset = Copy.objects.all().order_by('name')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'librarycopy'
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
            # schema = super(LibraryCopyResource, self).build_schema()
        
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
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
            filename = self._get_filename(readable_filter_hash)

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
                    select([func.count(None)]).select_from(_p).where(_p.c.copy_id==text('copy.copy_id'))),
                    
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

    @un_cache        
    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_list must be implemented')
                
    @write_authorization
    @un_cache
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):

        id_kwargs = self.get_id(deserialized, **kwargs)
        ScreensaverUser.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @un_cache
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs['schema']
        fields = schema['fields']
        initializer_dict = {}

        # TODO: wrapper for parsing
        logger.debug('fields: %r, deserialized: %r', fields.keys(), deserialized)
        for key in fields.keys():
            if deserialized.get(key, None) is not None:
                initializer_dict[key] = parse_val(
                    deserialized.get(key, None), key, fields[key]['data_type']) 

        id_kwargs = self.get_id(deserialized, **kwargs)
        
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
                errors = self.validate(deserialized, patch=True)
                if errors:
                    raise ValidationError(errors)
            except ObjectDoesNotExist:
                librarycopy = Copy.objects.create(
                    name=id_kwargs['copy_name'], library=library)
                errors = self.validate(deserialized, patch=False)
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
                    librarycopy.library.start_plate, librarycopy.library.end_plate, librarycopy.name )
    
                initial_plate_status = deserialized.get('initial_plate_status', None)
                if initial_plate_status:
                    self.get_plate_resource().validate({'status': initial_plate_status})
            
                initial_plate_well_volume = deserialized.get('initial_plate_well_volume', None)
                if initial_plate_well_volume:
                    initial_plate_well_volume = parse_val(
                        initial_plate_well_volume, 'initial_plate_well_volume', 'decimal')
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
                        for well in library.well_set.all().filter(plate_number=p.plate_number):
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


class PublicationResource(DbApiResource):

    class Meta:

        queryset = Publication.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        serializer = LimsSerializer()
        resource_name = 'publication'
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

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
        
        logger.info('params: %r', param_hash.keys())
        
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
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
            filename = self._get_filename(readable_filter_hash)
                  
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
                
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        schema = kwargs['schema']
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
        logger.info('parent_log: %r', parent_log)
        if parent_log is None and publication.screen is not None:
            parent_log = self.log_to_screen(
                publication.screen, publication, request)
            kwargs['parent_log'] = parent_log
            logger.info('1parent_log: %r', parent_log)
        logger.info('parent_log: %r', parent_log)
            
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
        new_data = self._get_detail_response(request,**kwargs_for_log)
        log = self.log_patch(request, None, new_data, **kwargs_for_log)
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
        schema = kwargs['schema']
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
                original_data = self._get_detail_response(request,**kwargs_for_log)
            except Exception as e:
                logger.exception('original state not obtained')
                original_data = {}

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
#         log.added_keys = json.dumps(original_data.keys(),cls=DjangoJSONEncoder)
#         log.diffs = json.dumps(original_data,cls=DjangoJSONEncoder)
        log.diffs = { k:[v,None] for k,v in original_data.items()}
        log.save()
        logger.info('delete, api log: %r', log)

        return tastypie.http.HttpNoContent()
        

class AttachedFileResource(DbApiResource):

    class Meta:

        queryset = AttachedFile.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        serializer = LimsSerializer()
        resource_name = 'attachedfile'
        always_return_data = True

    def __init__(self, **kwargs):
        
        super(AttachedFileResource, self).__init__(**kwargs)
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
                 r"(?P<attached_file_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/user/(?P<username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
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
        screen_resource = self.get_screen_resource()
        _screen_data = \
            screen_resource._get_detail_response_internal(**{
                'limit': 10,
                'facility_id': screen.facility_id,
                'exact_fields': ['attached_files'],
                })
        parent_log = screen_resource.make_log(request)
        parent_log.api_action = API_ACTION_PATCH
        parent_log.key = screen.facility_id
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
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        schema = kwargs['schema']
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

        parent_fields = set(['username', 'publication_id', 'screen_facility_id', 'reagent_id'])
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
                'username': msg, 'screen_facility_id': msg,
                'publication_id': msg, 'reagent_id': msg
            })
        if 'username' in initializer_dict:
            try:
                user = ScreensaverUser.objects.get(
                    username=initializer_dict['username'])
                initializer_dict['screensaver_user'] = user
            except ObjectDoesNotExist:
                raise Http404('username %r does not exist' 
                    % initializer_dict['username'])
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
                logger.warn('no such attribute on attached_file: %s:%r' 
                    % (key, val))

        parent_log = kwargs.get('parent_log', None)
        if parent_log is None and af.screen is not None:
            parent_log = \
                self.log_to_screen(af.screen, af, request, is_delete=False)
        # TODO: log_to_reagent
        # TODO: log_to_screensaver_user
        
        af.save()

        # get new state, for logging
        kwargs_for_log = { 
            'attached_file_id': af.attached_file_id, 
            'parent_log': parent_log 
        }
        new_data = self._get_detail_response(request,**kwargs_for_log)
        log = self.log_patch(request, None,new_data,**kwargs_for_log)
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
        attached_file_id = kwargs.pop('attached_file_id', None)
        if not attached_file_id:
            NotImplementedError('must provide a attached_file_id parameter')
        
        af = AttachedFile.objects.get(attached_file_id=attached_file_id)
        _dict = model_to_dict(af)

        # Create parent log
        parent_log = kwargs.get('parent_log', None)
        if parent_log is None and af.screen is not None:
            parent_log = \
                self.log_to_screen(af.screen, af, request, is_delete=True)

        af.delete()
        
        schema = kwargs['schema']
        id_attribute = resource = schema['id_attribute']

        log = self.make_log(request, **kwargs)
        log.ref_resource_name = self._meta.resource_name
        log.key = '/'.join([str(_dict[x]) for x in id_attribute])
        log.uri = '/'.join([self._meta.resource_name, log.key])
        log.parent_log = parent_log
        log.api_action = API_ACTION_DELETE
        log.added_keys = json.dumps(_dict.keys(), cls=DjangoJSONEncoder)
        log.diffs = json.dumps(_dict, cls=DjangoJSONEncoder)
        log.save()
        logger.debug('delete, api log: %r', log)

        return tastypie.http.HttpNoContent()

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

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
        
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
        username = param_hash.pop('username', None)
        if username:
            param_hash['username__eq'] = username
        screen_facility_id = param_hash.pop('screen_facility_id', None)
        if screen_facility_id:
            param_hash['screen_facility_id__eq'] = screen_facility_id
        attached_file_id = param_hash.pop('attached_file_id', None)
        if attached_file_id:
            param_hash['attached_file_id__eq'] = attached_file_id
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)
                  
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
                stmt = stmt.where(
                    _su.c.username == username)
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


class UserAgreementResource(AttachedFileResource):

    class Meta:

        queryset = AttachedFile.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        serializer = LimsSerializer()
        resource_name = 'useragreement'

    def __init__(self, **kwargs):
        
        self.screensaveruser_resource = None
        self.userchecklistitem_resource = None
        super(UserAgreementResource, self).__init__(**kwargs)

    def get_screensaveruser_resource(self):
        if not self.screensaveruser_resource:
            self.screensaveruser_resource = ScreensaverUserResource()
        return self.screensaveruser_resource
    
    def get_userchecklistitem_resource(self):
        if not self.userchecklistitem_resource:
            self.userchecklistitem_resource = UserChecklistItemResource()
        return self.userchecklistitem_resource

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<attached_file_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/user/(?P<username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
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
        kwargs['file_type__in'] = ['iccb_l_nsrb_rnai_user_agreement',
            '2010_iccb_l_nsrb_small_molecule_user_agreement']
        return self.build_list_response(request, **kwargs)
    
    @write_authorization
    @un_cache        
    @transaction.atomic
    def post_detail(self, request, **kwargs):
                
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        
        type = param_hash.pop('type', None)
        if not type:
            raise ValidationError(key='type', msg='required')
        
        if type == 'sm':
            type = '2010_iccb_l_nsrb_small_molecule_user_agreement'
        elif type == 'rna':
            type = 'iccb_l_nsrb_rnai_user_agreement'
        kwargs['type'] = type
        
        username = param_hash.pop('username', None)
        if not username:
            raise NotImplementedError('must provide a username parameter')
        
        admin_username = param_hash.get('admin_user', None)
        if not admin_username:
            admin_username = request.user.username
        admin_user = ScreensaverUser.objects.get(username=admin_username)
        
        dsl_usergroup = param_hash.pop('usergroup', None)
        if not dsl_usergroup:
            raise ValidationError(key='usergroup', msg='must be set')

        user_resource = self.get_screensaveruser_resource()
        kwargs_for_user = {
            'username': username,
            'exact_fields': ['username', 'usergroups'] 
        }

        user_data = user_resource._get_detail_response_internal(**kwargs_for_user)
        current_groups = user_data.get('usergroups',[]) or []
        current_groups = set(current_groups)
        small_molecule_usergroups = set([
            'smDsl1MutualScreens','smDsl2MutualPositives','smDsl3SharedScreens'])
        rna_usergroups = set([
            'rnaiDsl1MutualScreens','rnaiDsl2MutualPositives','rnaiDsl3SharedScreens'])
        if (dsl_usergroup not in small_molecule_usergroups and 
                dsl_usergroup not in rna_usergroups):
            raise ValidationError(key='usergroup',
                msg='must be in %r or %r' % (small_molecule_usergroups,rna_usergroups))
        apilog = self.make_log(request)
        apilog.ref_resource_name = 'screensaveruser'
        apilog.uri = [apilog.ref_resource_name, username]
        apilog.key = username
#         apilog.diff_keys = ['data_sharing_level']
        
        # - the user agreement is an attached file to the user
        super(UserAgreementResource,self).post_detail(request, **kwargs)

        # - the data sharing level group is assigned to the user
        
        current_val = current_groups & (small_molecule_usergroups | rna_usergroups)
        if dsl_usergroup in small_molecule_usergroups:
            current_groups = current_groups - small_molecule_usergroups
            checklist_item_name = 'current_small_molecule_user_agreement_active'
        if dsl_usergroup in rna_usergroups:
            current_groups = current_groups - rna_usergroups
            checklist_item_name = 'current_rnai_user_agreement_active'
        current_groups.add(dsl_usergroup)
        new_val = current_groups & (small_molecule_usergroups | rna_usergroups)
        apilog.diffs = { 'data_sharing_level': [','.join(current_val), ','.join(new_val)] }
        apilog.save()
        
        user_data['usergroups'] = list(current_groups)
        logger.info('patch user with data: %r', user_data)
        user_resource.patch_obj(request, user_data)

        # create a checklist item for the user agreement
        user_checklist_item_resource = self.get_userchecklistitem_resource()
        checklist_data = {
            'username': username, 
            'item_name': checklist_item_name,
            'item_group': 'forms',
            'admin_username': admin_user.username,
            'status': 'completed',
            'status_date': timezone.now().strftime("%Y%m%d") 
            }
        
        user_checklist_item_resource.patch_obj(request, checklist_data)

        return tastypie.http.HttpCreated()

    # TODO: create a "UserAgreementResource.delete_detail"


# 
# select
# a.activity_id,
# (
# CASE WHEN (select activity_id from library_screening ls where ls.activity_id = a.activity_id) is not null THEN (
#    CASE WHEN (select is_for_external_library_plates from library_screening ls where ls.activity_id=a.activity_id) is true THEN 'External Library Screening' ELSE 'Library Screening' END )
#      WHEN (select activity_id from cherry_pick_screening cs where cs.activity_id = a.activity_id ) is not null THEN 'Cherry Pick Screening'
#      WHEN (select activity_id from cherry_pick_liquid_transfer cplt where cplt.activity_id = a.activity_id ) is not null THEN 'Cherry Pick Liquid Transfer'
# else 'unknown_activity_type_tell_sean'
# END 
# ) activity_type,
# a.date_of_activity date_performed,
# a.date_created date_recorded,
# (select first_name || ' ' || last_name from screensaver_user su where su.screensaver_user_id = a.performed_by_id) performed_by,
# s.facility_id,
# screen_type,
# lab.last_name || ', ' || lab.first_name as lab_head,
# '' as serviced_user,
# la.affiliation_name as lab_affiliation,
# min(fs.value) as funding_support,
# si.status as status,
# (select min(date_of_activity) from activity a join lab_activity la using(activity_id) where la.screen_id = s.screen_id) date_of_first_activity,
# (select max(date_of_activity) from activity a join lab_activity la using(activity_id) where la.screen_id = s.screen_id) date_of_last_activity,
# regexp_replace(a.comments, E'[\\s]+', ' ','g') as comments
# from activity a
# join lab_activity lac using(activity_id)
# join screen s using(screen_id)
# left join screensaver_user lab on(s.lab_head_id = lab.screensaver_user_id)
# left join lab_head lab2 on(lab.screensaver_user_id = lab2.screensaver_user_id)
# left join lab_affiliation la on(lab2.lab_affiliation_id = la.lab_affiliation_id)
# left join screen_funding_support_link fsl on(s.screen_id = fsl.screen_id)
# left join funding_support fs on(fs.funding_support_id = fsl.funding_support_id)
# left join screen_status_item si on(si.screen_id = s.screen_id)
# left join screen_result sr on(sr.screen_id = s.screen_id)
# where (si.status is null or si.status_date = (select max(status_date) from screen_status_item si1 where si1.screen_id = s.screen_id)) /*can result in > 1 row per screen, if 2 statuses have same max date*/
# group by s.facility_id,s.screen_id, screen_type, lab_head, s.date_created, si.status, si.status_date, sr.date_created, 
# a.activity_id, a.date_of_activity, a.date_created, a.performed_by_id, la.affiliation_name, a.comments
#  UNION ALL
# select
# a.activity_id,
# sa.service_activity_type as activity_type,
# a.date_of_activity date_performed,
# a.date_created date_recorded,
# (select first_name || ' ' || last_name from screensaver_user su where su.screensaver_user_id = a.performed_by_id) performed_by,
# s.facility_id,
# screen_type,
# ( CASE WHEN (select screensaver_user_id from lab_head lh where lh.screensaver_user_id = sa.serviced_user_id) is not null  THEN (
#     ( select lh.last_name || ', ' || lh.first_name
#         from lab_head lab join screensaver_user lh on lab.screensaver_user_id=lh.screensaver_user_id where lab.screensaver_user_id = sa.serviced_user_id ))
#     ELSE ( select lh.last_name || ', ' || lh.first_name 
#         from screening_room_user sru 
#         join screensaver_user lh on sru.lab_head_id=lh.screensaver_user_id 
#         where sru.screensaver_user_id = sa.serviced_user_id )
# END ) as lab_head,
# serviced.last_name || ', ' || serviced.first_name as serviced_user,
# ( CASE WHEN (select screensaver_user_id from lab_head lh where lh.screensaver_user_id = sa.serviced_user_id) is not null  THEN (
#      ( select la.affiliation_name
#         from lab_head lab join lab_affiliation la using(lab_affiliation_id) where lab.screensaver_user_id = sa.serviced_user_id ))
#        ELSE ( select la.affiliation_name from screening_room_user sru 
#         join lab_head lh on sru.lab_head_id=lh.screensaver_user_id
#         join lab_affiliation la using(lab_affiliation_id) 
#         where sru.screensaver_user_id = sa.serviced_user_id )
# END ) as lab_affiliation,
# min(fs.value) as funding_support,
# si.status as status,
# (select min(date_of_activity) from activity a join lab_activity la using(activity_id) where la.screen_id = s.screen_id) date_of_first_activity,
# (select max(date_of_activity) from activity a join lab_activity la using(activity_id) where la.screen_id = s.screen_id) date_of_last_activity,
# regexp_replace(a.comments, E'[\\s]+', ' ','g') as comments
# from activity a
# join service_activity sa using(activity_id)
# left join screen s on(serviced_screen_id=screen_id)
# left join screensaver_user serviced on (serviced.screensaver_user_id=sa.serviced_user_id)
# left join funding_support fs on(fs.funding_support_id = sa.funding_support_id)
# left join screen_status_item si on(si.screen_id = s.screen_id)
# left join screen_result sr on(sr.screen_id = s.screen_id)
# where (si.status is null or si.status_date = (select max(status_date) from screen_status_item si1 where si1.screen_id = s.screen_id)) /*can result in > 1 row per screen, if 2 statuses have same max date*/
# group by s.facility_id,s.screen_id, screen_type, lab_head, s.date_created, si.status, si.status_date, sr.date_created, 
# a.activity_id, a.date_of_activity, a.date_created, a.performed_by_id,serviced_user, sa.serviced_user_id, sa.service_activity_type, a.comments
# order by facility_id, activity_id

class ActivityResource(DbApiResource):
    '''
    Activity Resource is a combination of the LabActivity and the ServiceActivity
    '''

    class Meta:

        queryset = Activity.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'activity'
        serializer = LimsSerializer()
        ordering = []
        filtering = {}
        always_return_data = True 
        
    def __init__(self, **kwargs):

        self.service_activity_resource = None
        self.screen_resource = None
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

    def get_service_activity_resource(self):
        if not self.service_activity_resource:
            self.service_activity_resource = ServiceActivityResource()
        return self.service_activity_resource

    def get_screen_resource(self):
        if not self.screen_resource:
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
        _user_cte = ScreensaverUserResource.get_user_cte().cte(
            'users_%s' % alias_qualifier)
        _performed_by = _user_cte.alias('performed_by_%s' % alias_qualifier)
        _performed_by1 = _user_cte.alias('performed_by1_%s' % alias_qualifier)
        _activity = self.bridge['activity']
        _su = self.bridge['screensaver_user']
        _lhsu = _su.alias('lhsu_la')
        _sfs = self.bridge['screen_funding_supports']

#         affiliation_table = ScreensaverUserResource.get_lab_affiliation_cte()
#         affiliation_table = affiliation_table.cte('affiliation_%s' % alias_qualifier)
        lab_head_table = ScreensaverUserResource.get_lab_head_cte().cte('lab_heads')

        
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
            'screen_lab_affiliation': (
                select([lab_head_table.c.lab_affiliation])
                .select_from(lab_head_table)
                .where(lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id)),
            'screen_lab_name': (
                select([lab_head_table.c.lab_name_full])
                .select_from(lab_head_table)
                .where(lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id)),
            'screen_lab_head_username': (
                select([lab_head_table.c.username])
                .select_from(lab_head_table)
                .where(lab_head_table.c.screensaver_user_id==_screen.c.lab_head_id)),
#             'screen_lab_name':
#                 (select([func.array_to_string(array(
#                         [_lhsu.c.last_name, ', ', _lhsu.c.first_name, ' - ',
#                          affiliation_table.c.title,
#                          ' (', affiliation_table.c.category, ')']), '')])
#                     .select_from(
#                         _lhsu.join(
#                             affiliation_table,
#                             affiliation_table.c.affiliation_name 
#                                 == _lhsu.c.lab_head_affiliation))
#                     .where(_lhsu.c.screensaver_user_id == _screen.c.lab_head_id)),
#             'screen_lab_affiliation': 
#                 (select([_concat(
#                         affiliation_table.c.title, ' (',
#                         affiliation_table.c.category, ')')])
#                     .select_from(
#                         _lhsu.join(
#                             affiliation_table,
#                             affiliation_table.c.affiliation_name 
#                                 == _lhsu.c.lab_head_affiliation))
#                     .where(
#                         _lhsu.c.screensaver_user_id == _screen.c.lab_head_id)),
#             'screen_lab_head_username':
#                 (select([_lhsu.c.username])
#                     .select_from(_lhsu)
#                     .where(_lhsu.c.screensaver_user_id == _screen.c.lab_head_id)),
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
            'screen_lead_screener_username': (
                select([_su.c.username])
                    .select_from(_su)
                    .where(_su.c.screensaver_user_id 
                        == _screen.c.lead_screener_id)),
            'screen_date_of_last_activity': literal_column(
                '( (select date_of_activity '
                '  from activity '
                '  join lab_activity la using(activity_id) '
                '  where la.screen_id=screen.screen_id '
                '  UNION ALL'
                '  select date_of_activity '
                '  from activity '
                '  join service_activity sa using(activity_id) '
                '  where sa.serviced_screen_id=screen.screen_id )'
                '  order by date_of_activity desc LIMIT 1 )'),
        }

    def get_query(self, param_hash):
        
        # general setup
        schema = self.build_schema()
        
        manual_field_includes = set(param_hash.get('includes', []))
        # for join to screen query (TODO: only include if screen fields rqst'd)
        manual_field_includes.add('screen_id')
        manual_field_includes.add('activity_class')
        param_hash['includes'] = list(manual_field_includes)
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash)
              
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_params.append('date_of_activity')
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)

        # specific setup
        
        _activity = self.bridge['activity']
        # Create a UNION query of each subclass query:
        (field_hash_sa, columns_sa, stmt_sa, count_stmt_sa, filenamesa) = (
            self.get_service_activity_resource().get_query(param_hash))
        # if a field is not present in the subquery, create an empty field            
        sa_columns = []
        for key in [key for key, field in field_hash.items() 
            if field['scope'] == 'fields.activity'] :
            if key not in columns_sa:
                sa_columns.append(cast(literal_column("null"), 
                    sqlalchemy.sql.sqltypes.Text).label(key))
            else:
                sa_columns.append(literal_column(key))
        stmt_sa = stmt_sa.cte('serviceactivities')
        stmt1 = select(sa_columns).select_from(stmt_sa)

        (field_hash_la, columns_la, stmt_la, count_stmt_la) = (
            self.get_lab_activity_query(param_hash))
        # if a field is not present in the subquery, create an empty field            
        la_columns = []
        for key in [key for key, field in field_hash.items() 
            if field['scope'] == 'fields.activity'] :
            if key not in columns_la:
                logger.info(
                    ('programming error: get_lab_activity_query '
                     'is missing the col: %r'), key)
                la_columns.append(cast(
                    literal_column("null"), 
                    sqlalchemy.sql.sqltypes.Text).label(key))
            else:
                la_columns.append(literal_column(key))
        stmt_la = stmt_la.cte('labactivities')
        stmt2 = select(la_columns).select_from(stmt_la)
        
        if logger.isEnabledFor(logging.DEBUG):
            compiled_stmt = str(stmt2.compile(
                dialect=postgresql.dialect(),
                compile_kwargs={"literal_binds": True}))
            logger.info('compiled_stmt %s', compiled_stmt)

        stmt = stmt1.union_all(stmt2)
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)

        columns = { key:literal_column(key) for key in field_hash.keys()}
        return (field_hash, columns, stmt, count_stmt, filename)
    
    def get_custom_lab_activity_columns(self, alias_qualifier):

        _library_screening = self.bridge['library_screening']
        _cps = self.bridge['cherry_pick_screening']
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
            'type': activity_type_column,
            'activity_class': activity_type_column,
            }
        
    def get_lab_activity_query(self, param_hash):

        # general setup
        # schema for labactivity part of the query should only be the activity fields
        # (exclude the screen.fields)
        schema = deepcopy(self.build_schema())
        field_hash = schema['fields']
        field_hash = { key:val for key, val in field_hash.items() 
            if val['scope'] == 'fields.activity'}  
        schema['fields'] = field_hash
        
        manual_field_includes = set(param_hash.get('includes', []))
        manual_field_includes.add('screen_id')
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
              
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
        j = j.join(
            _library_screening,
            _la.c.activity_id == _library_screening.c.activity_id, isouter=True)
        j = j.join(_cps, _la.c.activity_id == _cps.c.activity_id, isouter=True)
                
        custom_columns = self.get_custom_lab_activity_columns('lab_activity')
        custom_columns.update(self.get_custom_columns('la'))
        
        base_query_tables = ['activity', 'lab_activity', 'screening', 'screen'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        stmt = select(columns.values()).select_from(j)
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

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
        
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
        
        try:
            (field_hash, columns, stmt, count_stmt, filename) = self.get_query(param_hash)
            
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)

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


class CherryPickLiquidTransferResource(ActivityResource):
    ''' DEPRECATED: new CherryPickRequest does not use cplt '''

    class Meta:

        queryset = Screening.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        resource_name = 'cplt'
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        super(CherryPickLiquidTransferResource, self).__init__(**kwargs)

    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<activity_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
    
    def get_custom_columns(self, alias_qualifier):
        '''
        Convenience method for subclasses: reusable custom columns
        @param alias_qualifier a sql compatible string used to name subqueries
            so that this method may be called multiple times to compose a query
        '''
        ccs = super(CherryPickLiquidTransferResource, self).get_custom_columns(
            alias_qualifier)
        return ccs

    def get_query(self, param_hash):

        # general setup
        schema = self.build_schema()
        
        manual_field_includes = set(param_hash.get('includes', []))
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash)
              
        order_params = param_hash.get('order_by', [])
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
        _cplt = self.bridge['cherry_pick_liquid_transfer']
        _screen = self.bridge['screen']
        _cpap = self.bridge['cherry_pick_assay_plate']
        _cherry_pick = self.bridge['cherry_pick_request']
        j = _a
        j = j.join(_la, _a.c.activity_id == _la.c.activity_id)
        j = j.join(_cplt, _cplt.c.activity_id == _la.c.activity_id)        
        j = j.join(_screen, _la.c.screen_id == _screen.c.screen_id)

        custom_columns = {
            'type': literal_column("'cplt'"),
            'activity_class': literal_column("'cplt'"),
            'cherry_pick_request_id': (
                select([_cpap.c.cherry_pick_request_id])
                    .select_from(_cpap)
                    .where(_cpap.c.cherry_pick_liquid_transfer_id 
                        == _a.c.activity_id)
                    .limit(1))
            }
        custom_columns.update(self.get_custom_columns('cplt'))
        
        base_query_tables = ['activity', 'lab_activity',
            'cherry_pick_liquid_transfer', 'screen', 'cherry_pick'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        stmt = select(columns.values()).select_from(j)
        compiled_stmt = str(stmt.compile(
            dialect=postgresql.dialect(),
            compile_kwargs={"literal_binds": True}))
        logger.info('compiled_stmt %s', compiled_stmt)
        # general setup
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        
        return (field_hash, columns, stmt, count_stmt,filename)


class CherryPickScreeningResource(ActivityResource):    
    ''' DEPRECATED: new CherryPickRequest does not use cherry pick screening '''

    class Meta:

        queryset = Screening.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        resource_name = 'cherrypickscreening'
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        super(CherryPickScreeningResource, self).__init__(**kwargs)

    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<activity_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

    def get_query(self, param_hash):
        
        # general setup
        schema = self.build_schema()
        
        manual_field_includes = set(param_hash.get('includes', []))
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash)
              
        order_params = param_hash.get('order_by', [])
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
        _cps = self.bridge['cherry_pick_screening']
        _cpr = self.bridge['cherry_pick_request']
        _screen = self.bridge['screen']
        j = _a
        j = j.join(_la, _a.c.activity_id == _la.c.activity_id)
        j = j.join(_screening, _screening.c.activity_id == _la.c.activity_id)
        j = j.join(_cps, _cps.c.activity_id == _la.c.activity_id)
        j = j.join(
            _cpr,
            _cpr.c.cherry_pick_request_id == _cps.c.cherry_pick_request_id)
        j = j.join(_screen, _la.c.screen_id == _screen.c.screen_id)
                
        custom_columns = \
            super(CherryPickScreeningResource, self).get_custom_columns('ls')
        custom_columns.update({
            'type': cast(
                literal_column("'cherrypickscreening'"), 
                sqlalchemy.sql.sqltypes.Text),
            'activity_class': 
                cast(
                    literal_column("'cherrypickscreening'"), 
                    sqlalchemy.sql.sqltypes.Text),
            })

        base_query_tables = ['activity', 'lab_activity', 'screening',
            'cherry_pick_screening', 'cherry_pick_request', 'screen'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        stmt = select(columns.values()).select_from(j)
        compiled_stmt = str(stmt.compile(
            dialect=postgresql.dialect(),
            compile_kwargs={"literal_binds": True}))
        logger.info('compiled_stmt %s', compiled_stmt)
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        
        return (field_hash, columns, stmt, count_stmt, filename)
        
    @read_authorization
    def build_list_response(self, request, **kwargs):
 
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
         
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
 
        try:
             
            (field_hash, columns, stmt, count_stmt,filename) = self.get_query(param_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
  
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

class LibraryScreeningResource(ActivityResource):    

    class Meta:

        queryset = LibraryScreening.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        resource_name = 'libraryscreening'
        alt_resource_name = 'externallibraryscreening'
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        super(LibraryScreeningResource, self).__init__(**kwargs)
        
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
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.alt_resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<activity_id>([\d]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/" 
                 r"(?P<activity_id>([\d]+))%s$")
                    % (self._meta.alt_resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

    def get_query(self, schema, param_hash):
        '''  LibraryScreeningResource
        '''

        # general setup
#         schema = self.build_schema()
        
        manual_field_includes = set(param_hash.get('includes', []))
        
        library_plates_screened_search = param_hash.pop('library_plates_screened__contains', None)
        if library_plates_screened_search:
            manual_field_includes.add('library_plates_screened')
            
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
              
        order_params = param_hash.get('order_by', [])
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

                 
        custom_columns = \
            super(LibraryScreeningResource, self).get_custom_columns('ls')
        custom_columns.update({
            'type': cast(
                literal_column("'libraryscreening'"), 
                sqlalchemy.sql.sqltypes.Text),
            'activity_class': cast(
                literal_column("'libraryscreening'"), 
                sqlalchemy.sql.sqltypes.Text),
            'libraries_screened_count': literal_column(
                '(select count(distinct(l.*)) from library l '
                'join copy using(library_id) join plate using(copy_id) '
                'join assay_plate ap using(plate_number) '
                'where library_screening_id='
                'library_screening.activity_id)' ),
             'library_plates_screened_count': literal_column(
                 '(select count(distinct(ap.plate_number)) '
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

        filename = self._get_filename(readable_filter_hash, **extra_params)
        
        return (field_hash, columns, stmt, count_stmt,filename)
        
    @read_authorization
    def build_list_response(self, request, **kwargs):
 
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
         
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
 
        try:
             
            (field_hash, columns, stmt, count_stmt,filename) = self.get_query(schema, param_hash)
             
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
                        _ap.c.plate_number
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

    @un_cache        
    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_detail must be implemented')
    
    def validate(self, _dict, patch=False):

        errors = ActivityResource.validate(self, _dict, patch=patch)
        if _dict.get('library_plates_screened', None):
            if bool(_dict.get('is_for_external_library_plates', False)):
                errors['library_plates_screened'] = (
                    'cannot specifiy library plates if '
                    '"is_for_external_library_plates"')
        
        return errors

#     @write_authorization
#     @un_cache
#     @transaction.atomic
#     def post_detail(self, request, **kwargs):
#         response = super(LibraryScreeningResource, self).post_detail(request, **kwargs)
#         if response.status_code == 200:
#             _data = self._meta.serializer.deserialize(
#                 LimsSerializer.get_content(response), response['Content-Type'])
#             meta = { 
#                 'Result': { 
#                     'Library Plates updated': log.child_logs.all().count(),
#                     'Volume per well transferred from Plates': 
#                         _data['volume_transferred_per_well_from_library_plates'],
#                 }
#             }
#             _data = { 'objects': _data }
#             _data['meta'] = meta
#             if not self._meta.always_return_data:
#                 return self.build_response(
#                     request, {'meta': meta }, response_class=HttpResponse, **kwargs)
#             else:
#                 return self.build_response(
#                     request, _data, response_class=HttpResponse, **kwargs)
#         else:
#             return response

    def build_patch_detail(self, request, deserialized, log=None, **kwargs):
        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        schema = kwargs['schema']
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(deserialized,validate=True,**kwargs)

        original_data = self._get_detail_response(request,**kwargs_for_log)
        logger.info('original data: %r', original_data)
        
        # How to log parent entity data? rather than do this, users will have
        # to consult the "Screening Summary" history table
        # original_screen_data = self.get_screen_resource()._get_detail_response_internal(**{
        #     'facility_id': original_data['screen_facility_id']})

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
        log.uri = '/'.join([log.ref_resource_name,log.key])
        log.save()
        logger.info('log saved: %r', log)
        
        # TODO: create a log for the parent screen

        new_data = self._get_detail_response(request,**kwargs_for_log)
        logger.info('new data: %r', new_data)
        self.log_patch(request, original_data,new_data,log=log, **kwargs)
        # FIXME: Set the log URI using the containing screen URI
        logger.info('log created: %r, %r', log, log.diffs)
        log.save()

        meta = { 
            API_MSG_SCREENING_PLATES_UPDATED: log.child_logs.all().count(),
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 
                new_data['library_plates_screened_count'],
            'Volume per well transferred from Plates': 
                new_data['volume_transferred_per_well_from_library_plates'],
        }
        meta.update(plate_meta)
        _data = { API_RESULT_DATA: [new_data] }
        _data[API_RESULT_META] = { API_MSG_RESULT: meta }
        
#         new_data = { API_RESULT_DATA: new_data, }
#         if API_RESULT_META in patch_result:
#             new_data[API_RESULT_META] = patch_result[API_RESULT_META]
        return _data

    def build_post_detail(self, request, deserialized, log=None, **kwargs):
        kwargs_for_log = self.get_id(deserialized,validate=False,**kwargs)
        
        schema = kwargs['schema']
        id_attribute = schema['id_attribute']
        
        logger.info('post detail: %r, %r', kwargs_for_log, id_attribute)

        original_data = None
        log = self.make_log(request)
        log.save()
        if kwargs_for_log and len(kwargs_for_log.items())==len(id_attribute):
            # A full id exists, query for the existing state
            try:
                original_data = self._get_detail_response(request,**kwargs_for_log)
            except Exception, e: 
                logger.exception('exception when querying for existing obj: %s', 
                    kwargs_for_log)
            
            if original_data is not None and len(original_data) != 0:
                raise ValidationError({ 
                    k: '%r Already exists' % v for k,v in kwargs_for_log.items() })
            original_data = None
        
            # NOTE: create a log if possible, with id_attribute, for downstream
            log.key = '/'.join([str(kwargs_for_log[x]) for x in id_attribute])
            log.uri = '/'.join([log.ref_resource_name,log.key])
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
        log.uri = '/'.join([log.ref_resource_name,log.key])
        log.save()

        # get new state, for logging
        new_data = self._get_detail_response(request,**kwargs_for_log)
        logger.info('new data: %r', new_data)
        if not new_data:
            raise BadRequest('no data found for the new obj created by post: %r', obj)
        self.log_patch(request, original_data,new_data,log=log, **kwargs)

        # FIXME: Set the log URI using the containing screen URI
        log.save()
        meta = { 
            API_MSG_SCREENING_PLATES_UPDATED: log.child_logs.all().count(),
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 
                new_data['library_plates_screened_count'],
            'Volume per well transferred from Plates': 
                obj.volume_transferred_per_well_from_library_plates,
        }
        meta.update(plate_meta)
        _data = { API_RESULT_DATA: [new_data] }
        _data[API_RESULT_META] = { API_MSG_RESULT: meta }
        
        return _data

    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs['schema']
        fields = schema['fields']
        initializer_dict = {}
        ls_log = kwargs.get('log', None)
        if ls_log is None:
            raise BadRequest(
                'library screening log should be created by the callee of patch_obj')

        # TODO: wrapper for parsing
        # FIXME: move parsing until after validation
        # FIXME: parse and validate only editable fields

        id_kwargs = self.get_id(deserialized, **kwargs)
        logger.info('id_kwargs: %r', id_kwargs)
        patch = bool(id_kwargs)
        initializer_dict = self.parse(deserialized, create=not patch)
        errors = self.validate(initializer_dict, patch=patch)
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

        _key = 'performed_by_username'
        _val = deserialized.get(_key, None)
        if _val:
            try:
                performed_by_user = ScreensaverUser.objects.get(username=_val)
                initializer_dict['performed_by'] = performed_by_user
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=_key,
                    msg='does not exist: {val}'.format(val=_val))
        logger.info('initializer_dict: %r', initializer_dict)
        try:
            library_screening = None
            current_volume_tranferred_per_well = 0
            if patch:
                try:
                    logger.info('%r', id_kwargs)
                    library_screening = LibraryScreening.objects.get(
                        pk=id_kwargs['activity_id'])
                    current_volume_tranferred_per_well = \
                        library_screening.volume_transferred_per_well_from_library_plates
                except ObjectDoesNotExist:
                    raise Http404(
                        'library_screening does not exist for: %r', id_kwargs)
            else:
                library_screening = LibraryScreening()
            
            model_field_names = [
                x.name for x in library_screening._meta.get_fields()]
            for key, val in initializer_dict.items():
                if key in model_field_names:
                    setattr(library_screening, key, val)

            library_screening.save()
            new_volume_transferred_per_well = \
                library_screening.volume_transferred_per_well_from_library_plates
            library_plates_screened = deserialized.get(
                'library_plates_screened', [])
            if not library_plates_screened:
                raise ValidationError(
                    key='library_plates_screened', 
                    msg='required')
            # TODO: Test override
            override_param = parse_val(
                kwargs.get(API_PARAM_OVERRIDE, False),
                    API_PARAM_OVERRIDE, 'boolean')
            plate_meta = self._set_assay_plates(
                request, schema, 
                library_screening, library_plates_screened,
                current_volume_tranferred_per_well, 
                new_volume_transferred_per_well,
                ls_log, override_param)
            
            ls_log.json_field = json.dumps(plate_meta)
            logger.info('parent_log: %r', ls_log)
            self.create_screen_screening_statistics(library_screening.screen)
            self.create_screened_experimental_well_count(library_screening)
            
            return { API_RESULT_OBJ: library_screening, API_RESULT_META: plate_meta}
        except Exception, e:
            logger.exception('on patch_obj')
            raise e
        
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
            # Not used:
            # library_plates_data_analyzed_count = models.IntegerField(null=False, default=0)
            screen.save()
    
    @transaction.atomic    
    def _set_assay_plates(
            self, request, schema, library_screening, library_plates_screened,
            current_volume_tranferred_per_well, new_volume_transferred_per_well,
            ls_log, override_param):
        '''
        - Create new assay plates
        - Adjust librarycopyplate volume
        - Adjust librarycopyplate screening count
        - Create librarycopyplate logs
        - TODO: create copy logs
        '''
        logger.info('set assay plates screened for: %r, %r', 
            library_screening, library_plates_screened)
        # parse library_plate_ranges
        # E.G. Regex: /(([^:]+):)?(\w+):(\d+)-(\d+)/
        regex_string = schema['fields']['library_plates_screened']['regex']
        matcher = re.compile(regex_string)
        
        # Validate the library_plates_screened
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
        
        # validate the plate ranges
        logger.info(
            'get the referenced plates for: %r', library_plates_screened)
        plate_ranges = []
        plate_keys = set()
        plate_numbers = set()
        plate_search = []
        for _data in library_plates_screened:
            logger.debug('lps data: %r', _data)
            try:
                copy_name = _data['copy_name']
                library_short_name = _data['library_short_name']
                copy = Copy.objects.get(
                    name=copy_name,
                    library__short_name=library_short_name)
            except ObjectDoesNotExist:
                raise ValidationError(
                    key='library_plates_screened',
                    msg='{copy} does not exist: {val}'.format(
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
                             '{library_short_name},{screen_facility_id}'
                             ).format(**range))
                plate_range = range(
                    start_plate.plate_number, end_plate.plate_number + 1)
                if plate_numbers & set(plate_range):
                    raise ValidationError(
                        key='library_plates_screened',
                        msg=('A plate number can only be screened once per '
                            'Library Screening: {start_plate}-{end_plate}'
                            ).format(**_data))
                plate_numbers.update(plate_range)
                plate_keys.update([ '%s/%d' % (copy.name, plate_number) 
                    for plate_number in plate_range])
                logger.info('find the plate range: %s-%s', 
                    start_plate.plate_number, end_plate.plate_number)
                plate_ranges.append(Plate.objects.all().filter(
                    copy=copy,
                    plate_number__range=(
                        start_plate.plate_number, end_plate.plate_number)))

                # Also, create a search criteria to poll current plate state
                plate_search.append({
                    'copy_id': copy.copy_id,
                    'plate_number__range': 
                        [start_plate.plate_number, end_plate.plate_number] })
            
            except ObjectDoesNotExist:
                logger.exception('plate range error')
                raise ValidationError({
                    'library_plates_screened':
                    ('plate range not found: {start_plate}-{end_plate}'
                        ).format(**_data),
                    'copy_name': _data['copy_name'] })
        logger.debug('plate keys: %r, plate_numbers: %r', plate_keys, plate_numbers)
        logger.debug('plate search 1: %r', plate_search)
        
        # Create a search criteria to poll the current plate state
        # TODO: cache and log the copy state as well
        existing_ranges = {}
        if library_screening.assayplate_set.exists():
            for ap in library_screening.assayplate_set.all():
                # also, add to the criteria for search
                plate_range = existing_ranges.get(ap.plate.copy_id, [])
                if not plate_range:
                    plate_range = [ap.plate.plate_number, ap.plate.plate_number]
                elif plate_range[0] > ap.plate.plate_number:
                    plate_range[0]=ap.plate.plate_number
                elif plate_range[1] < ap.plate.plate_number:
                    plate_range[1] = ap.plate.plate_number
                existing_ranges[ap.plate.copy_id] = plate_range
        # Cache plate data
        logger.debug('plate search 2: %r', plate_search)
        plate_search.extend([{'copy_id': k, 'plate_number__range': v} 
            for k,v in existing_ranges.items()])
        logger.info('plate_search: %r', plate_search)    
        _original_plate_data = self.get_plate_resource()._get_list_response(
            request, search_data=plate_search)

        # Find extant plates, find and remove deleted assay plates        
        extant_plates = set()
        deleted_plates = set()
        if library_screening.assayplate_set.exists():
            for ap in library_screening.assayplate_set.all():
                found = False
                plate_key = '%s/%d' % (ap.plate.copy.name, ap.plate.plate_number)
                if ap.plate_number in plate_numbers:
                    if plate_key in plate_keys:
                        found = True
                        extant_plates.add(ap.plate)
                    else:
                        # if not found, then it is a different copy, same number,
                        # should be caught by validation
                        raise Exception(
                            'programming error: overlapping plate range')
                if not found:
                    # 20161020: no longer tracking data_load actions for an assay plate
                    # so this is removed
                    # if ap.screen_result_data_loading:
                    #     raise ValidationError(
                    #         key='library_plates_screened',
                    #         msg=(
                    #             'Assay plate has data and cannot be removed: %d'
                    #             ) % ap.plate_number)
                    deleted_plates.add(ap.plate)
                    ap.delete()       
        logger.info('deleted plates: %r', deleted_plates)

        library_warnings = []
        # Create assay plates
        
        # TODO: review SS1 policy on screening plates
        created_plates = set()
        for plate_range in plate_ranges:
            for replicate in range(library_screening.number_of_replicates):
                for plate in plate_range:
                    if plate not in extant_plates:
                        if plate.copy.library.screening_status != 'allowed':
                            warn_msg = (
                                'plate: "%s/%s"; library: %r, status is not "allowed"'
                                    % (plate.copy.name, plate.plate_number, 
                                           plate.copy.library.short_name))

                            if override_param is not True:
                                raise ValidationError({
                                    'library_plates_screened':warn_msg,
                                    API_PARAM_OVERRIDE: 'required'
                                })
                            else:
                                library_warnings.append(warn_msg)
                        if plate.status != 'available':
                            raise ValidationError(
                                key='library_plates_screened',
                                msg=  'plate: "%s/%s"; status is not "available"'\
                                    % (plate.copy.name,plate.plate_number))
                        if plate.cplt_screening_count > 0:
                            # Volume tracking will not work on the plate level
                            # after cherry pick volumes have been taken from the
                            # copy wells
                            raise ValidationError(
                                key='library_plates_screened',
                                msg=  'plate: "%s/%s"; may not be screened after '
                                'cherry pick screenings (%d) have been performed'\
                                    % (plate.copy.name,plate.plate_number, 
                                        plate.cplt_screening_count))
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
        
        # Update the plate screening related fields:
        # TODO/Review policy: update the copy_wells if exist
        for plate in extant_plates:
            remaining_well_volume = plate.remaining_well_volume or Decimal(0)
            remaining_well_volume += current_volume_tranferred_per_well
            remaining_well_volume -= new_volume_transferred_per_well
            plate.remaining_well_volume = remaining_well_volume
            plate.save()
            
            # NOTE: usually, screening copies should not have copywells
            # TODO: implement this if screening policy is changed to allow
            if plate.copywell_set.exclude(volume=F('initial_volume')).exists():
                raise NotImplementedError(
                    'Cannot create a library screening if copy wells have been adjusted'
                    ', plate: %r' % plate)
                # for cw in p.copywell_set.all():
                #     remaining_well_volume = cw.volume or Decimal(0)
                #     remaining_well_volume += current_volume_tranferred_per_well
                #     remaining_well_volume -= new_volume_transferred_per_well
                #     cw.volume = remaining_well_volume
        
        # Update the Plates affected:
        # - plate.screening_count and plate.remaining_well_volume 
        for plate in created_plates:
            
            plate.screening_count = (plate.screening_count or 0) + 1
            remaining_well_volume = plate.remaining_well_volume or Decimal(0)
            remaining_well_volume -= new_volume_transferred_per_well
            plate.remaining_well_volume = remaining_well_volume
            plate.save()
            if plate.copywell_set.exclude(volume=F('initial_volume')).exists():
                raise NotImplementedError(
                    'Cannot create a library screening if copy wells have been adjusted'
                    ', plate: %r' % plate)
                        
        for plate in deleted_plates:
            plate.screening_count -= 1
            remaining_well_volume = plate.remaining_well_volume or Decimal(0)
            remaining_well_volume += current_volume_tranferred_per_well
            plate.remaining_well_volume = remaining_well_volume
            plate.save()
            if plate.copywell_set.exclude(volume=F('initial_volume')).exists():
                raise NotImplementedError(
                    'Cannot create a library screening if copy wells have been adjusted'
                    ', plate: %r' % plate)

        # Fetch the new Plate state: log plate volume changes, screening count
        _new_plate_data = self.get_plate_resource()._get_list_response(
            request, search_data=plate_search)
        plate_logs = self.get_plate_resource().log_patches(
            request, _original_plate_data, _new_plate_data,
            parent_log=ls_log, api_action=API_ACTION_PATCH)
        logger.info('plate_logs created: %r', plate_logs)
        
        # TODO: log the copy state as well
        
        logger.debug('plates: updated: %r, created: %r, deleted: %r', 
            extant_plates, created_plates, deleted_plates)
        
        # return meta information
        meta = {
            API_MSG_SCREENING_ADDED_PLATE_COUNT: len(created_plates),
            API_MSG_SCREENING_EXTANT_PLATE_COUNT : len(extant_plates),
            API_MSG_SCREENING_DELETED_PLATE_COUNT: len(deleted_plates)
            }
        if library_warnings:
            meta[API_MSG_WARNING] = library_warnings
        return meta
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, log=None, **kwargs):
        
        raise NotImplementedError(
            'Library Screening may not be deleted - remove all plate ranges to negate')
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
            

class ServiceActivityResource(ActivityResource):    

    class Meta:

        queryset = ServiceActivity.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['digested_password']
        resource_name = 'serviceactivity'
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
            url(r"^(?P<resource_name>%s)/for_user/(?P<serviced_username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            url(r"^(?P<resource_name>%s)/for_user/(?P<serviced_username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]    

    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs['schema']
        fields = schema['fields']
        initializer_dict = {}
        # TODO: wrapper for parsing
        for key in fields.keys():
            if deserialized.get(key, None):
                initializer_dict[key] = parse_val(
                    deserialized.get(key, None), key, fields[key]['data_type']) 
        
        serviced_username = deserialized.get('serviced_username', None)
        if not serviced_username:
            raise Exception('serviced_username not specified %s' % deserialized)
        activity_type = deserialized.get('type', None)
        if not activity_type:
            raise Exception('activity type not specified %s' % deserialized)
        # TODO: drive this from the metadata "field" property
        initializer_dict['service_activity_type'] = activity_type
        performed_by_username = deserialized.get('performed_by_username', None)
        if not performed_by_username:
            raise ValidationError(
                key='performed_by_username',
                msg='required')
        try:
            serviced_user = ScreensaverUser.objects.get(
                username=serviced_username)
            initializer_dict['serviced_user'] = serviced_user
        except ObjectDoesNotExist:
            logger.exception(
                'serviced_user/username does not exist: %s' % serviced_username)
            raise
        try:
            performed_by_user = ScreensaverUser.objects.get(
                username=performed_by_username)
            initializer_dict['performed_by_user_id'] = performed_by_user.pk
        except ObjectDoesNotExist:
            logger.exception(
                'admin_user/username does not exist: %s' % admin_username)
            raise

        try:
            activity = None
            service_activity = None
            if 'activity_id' in initializer_dict:
                try:
                    activity_id = initializer_dict['activity_id']
                    activity = Activity.objects.get(pk=activity_id)
                    service_activity = ServiceActivity.objects.get(
                        activity=activity)
                except ObjectDoesNotExist:
                    logger.error('Activity does not exist: %s' % activity_id)
                    raise Exception('Activity does not exist: %s' % activity_id)
            else:
                activity = Activity.objects.create(
                    performed_by=performed_by_user,
                    date_of_activity=initializer_dict['date_of_activity'])
            for key, val in initializer_dict.items():
                if hasattr(activity, key):
                    # note: setattr only works for simple attributes, not foreign keys
                    setattr(activity, key, val)
                else:
                    logger.warn(
                        'no such attribute on activity: %s:%r' % (key, val))
            activity.save()
            if not service_activity:
                service_activity = ServiceActivity.objects.create(
                    activity=activity,
                    serviced_user=serviced_user)
                logger.info('created service_activity: %s' % service_activity)
            logger.info('initializer dict %s' % initializer_dict)
            for key, val in initializer_dict.items():
                if hasattr(service_activity, key):
                    # note: setattr only works for simple attributes, not foreign keys
                    setattr(service_activity, key, val)
                else:
                    logger.warn('no such attribute on service_activity: %s:%r' 
                        % (key, val))
            service_activity.save()
            return { API_RESULT_OBJ: service_activity }
        except Exception, e:
            logger.exception('on patch_obj')
            raise e
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def delete_obj(self, request, deserialized, **kwargs):

        activity_id = kwargs.get('activity_id', None)
        if activity_id:
            try:
                ServiceActivity.objects.get(
                    activity__activity_id=activity_id).delete()
            except ObjectDoesNotExist:
                logger.warn('no such ServiceActivity: %s' % activity_id)
                raise Exception(
                    'ServiceActivity for activity_id: %s not found' % activity_id)
        else:
            raise Exception(
                'ServiceActivity delete action requires an activity_id %s' 
                % kwargs)
        
    def get_query(self, param_hash):

        schema = self.build_schema()
        logger.debug('serviceactivity fields: %r', schema['fields'].keys())
        
        # general setup
        alias_qualifier = 'sa'
        manual_field_includes = set(param_hash.get('includes', []))
        # for join to screen query (TODO: only include if screen fields rqst'd)
        manual_field_includes.add('screen_id')
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash)
              
        order_params = param_hash.get('order_by', [])
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
        _user_cte = ScreensaverUserResource.get_user_cte().cte('users_serviced')
        _serviced = _user_cte.alias('serviced_user')
        
        j = _a
        j = j.join(_sa, _a.c.activity_id == _sa.c.activity_id)
        j = j.join(
            _serviced,
            _sa.c.serviced_user_id == _serviced.c.screensaver_user_id)
        j = j.join(
            _screen,
            _sa.c.serviced_screen_id == _screen.c.screen_id, isouter=True)
        
        custom_columns = \
            super(ServiceActivityResource, self).get_custom_columns('sa')
        custom_columns.update({
            'activity_class': cast(
                literal_column("'serviceactivity'"), 
                sqlalchemy.sql.sqltypes.Text),
            'serviced_user': _serviced.c.name,
            'serviced_username': _serviced.c.username,
            })

        base_query_tables = ['activity', 'service_activity', 'screen'] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
        
        stmt = select(columns.values()).select_from(j)

        # general setup
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        
        return (field_hash, columns, stmt, count_stmt, filename)


class ScreenAuthorization(UserGroupAuthorization):
    
    def _is_screen_authorized(self, screen, screensaver_user, permission_type):
        
        authorized = super(ScreenAuthentication,self)._is_resource_authorized(
            'screen', screensaver_user.user.user, permission_type)
        if not authorized and permission_type == 'read':
            if screensaver_user in screen.collaborators:
                authorized = True
            elif screensaver_user == screen.lead_screener:
                authorized = True
            elif screensaver_user == screen.lab_head:
                authorized = True
        return authorized
    
class ScreenResource(DbApiResource):
    
    class Meta:

        queryset = Screen.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = ScreenAuthorization()
        resource_name = 'screen'
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 
        
    def __init__(self, **kwargs):
        super(ScreenResource, self).__init__(**kwargs)
        self.apilog_resource = None
        self.cpr_resource = None
        self.activity_resource = None
        self.libraryscreening_resource = None
         
    def clear_cache(self):
        logger.info('clear screen caches')
        DbApiResource.clear_cache(self)
        caches['screen'].clear()
   
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
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/cherrypicks%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_cherrypickview'),
                name="api_dispatch_screen_cherrypickview"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w]+))/copyplates%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_screen_librarycopyplateview'),
                name="api_dispatch_screen_librarycopyplateview"),
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
        ]    
        
    def dispatch_screen_detail_uiview(self, request, **kwargs):
        ''' 
        Special method to populate nested entities for the UI 
        - bypasses the "dispatch" framework call
        -- must be authenticated and authorized
        '''
        self.is_authenticated(request)
        if not self._meta.authorization._is_resource_authorized(
                self._meta.resource_name,request.user,'read'):
            raise ImmediateHttpResponse(
                response=HttpForbidden(
                    'user: %s, permission: %s/%s not found' 
                    % (request.user,self._meta.resource_name,'read')))
        
        logger.info('dispatch_screen_detail_uiview')
        
        if request.method.lower() != 'get':
            return self.dispatch('detail', request, **kwargs)
        
        facility_id = kwargs.get('facility_id', None)
        if not facility_id:
            raise NotImplementedError('must provide a facility_id parameter')
        
        cache_key = 'detail_ui_%s' % facility_id
        screen_cache = caches['screen']
        _data = screen_cache.get(cache_key)
        
        if not _data:
            logger.info('cache key not set: %s', cache_key)
            _data = self._get_detail_response(request, **kwargs)
            if not _data:
                return Http404
            else:
                # response = self.dispatch('detail', request, format='json', **kwargs )
                # if response.status_code == 200:
                #     _data = self._meta.serializer.deserialize(
                #         JSON_MIMETYPE, 
                #         response['Content-Type'])
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
                        'limit': 10,
                        'screen_facility_id__eq': _data['facility_id'],
                        'order_by': ['-date_requested'],
                        'exact_fields': [
                            'cherry_pick_request_id','date_requested', 'requested_by_name'],
                        })
                _data['cherry_pick_request_data'] = _cpr_data
                
                _latest_activities_data = \
                    self.get_activity_resource()._get_list_response_internal(**{
                        'screen_facility_id__eq': _data['facility_id'],
                        'limit': 1,
                        'order_by': ['-date_of_activity'],
                        'exact_fields': ['activity_id','type', 'performed_by_name'],
                        })
                _data['latest_activities_data'] = _latest_activities_data
                
                # TODO: attached files
                
                # TODO: publications
                
                screen_cache.set(cache_key, _data)
        else:
            logger.info('cache key set: %s', cache_key)
                
        # Serialize
        # FIXME: refactor to generalize serialization:
        # see build_response method (needs rework)
        content_type = self.get_accept_content_type(
            request,format=kwargs.get('format', None))
        if content_type in [XLS_MIMETYPE,CSV_MIMETYPE]:
            _data = {'objects': [_data]}
        response = HttpResponse(
            content=self.get_serializer().serialize(
                _data, content_type),
            content_type=content_type)
        if content_type == XLS_MIMETYPE:
            response['Content-Disposition'] = \
                'attachment; filename=%s.xlsx' % self._get_filename(
                    self.build_schema(), kwargs)
        if content_type == CSV_MIMETYPE:
            response['Content-Disposition'] = \
                'attachment; filename=%s.csv' % self._get_filename(
                    self.build_schema(), kwargs)
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
        kwargs['screen_facility_id__eq'] = kwargs.pop('facility_id')
        return self.get_cpr_resource().dispatch('list', request, **kwargs)    
        
    def dispatch_screen_libraryview(self, request, **kwargs):
        kwargs['for_screen_id'] = kwargs.pop('facility_id')
        return LibraryResource().dispatch('list', request, **kwargs)    

    def dispatch_screen_librarycopyplateview(self, request, **kwargs):
        kwargs['for_screen_id'] = kwargs.pop('facility_id')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)    

    # NOTE: no longer supporting plates loaded stats
    # def dispatch_screen_lcp_loadedview(self, request, **kwargs):
    #     # FIXME: copyplatesloaded no longer works - 20160607
    #     # because we are not creating "assay_plates" for screen results anymore
    #     # can this be modified to show virtual "plates" loaded?
    #     kwargs['loaded_for_screen_id'] = kwargs.pop('facility_id')
    #     return LibraryCopyPlateResource().dispatch('list', request, **kwargs)    

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

    def get_query(self, schema, param_hash):
        
        DEBUG_SCREEN = False or logger.isEnabledFor(logging.DEBUG)
        screens_for_username = param_hash.get('screens_for_username', None)
        # general setup

        facility_id = param_hash.pop('facility_id', None)
        if facility_id:
            param_hash['facility_id__eq'] = facility_id
        
        manual_field_includes = set(param_hash.get('includes', []))
        # for joins
        manual_field_includes.add('screen_id')
        manual_field_includes.add('has_screen_result')
        
        extra_params = {}
        if screens_for_username:
            screener_role_cte = ScreenResource.get_screener_role_cte().cte(
                'screener_roles1')
            manual_field_includes.add('screensaver_user_role')
            extra_params['user'] = screens_for_username
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(readable_filter_hash, **extra_params)
              
        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
         
        # specific setup
        base_query_tables = ['screen', 'screen_result']
        _screen = self.bridge['screen']
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
        # _srua = self.bridge['screen_result_update_activity']
        _screen_keyword = self.bridge['screen_keyword']
        _screen_cell_lines = self.bridge['screen_cell_lines']
        _library_screening = self.bridge['library_screening']
        _cp_screening = self.bridge['cherry_pick_screening']
        _lab_activity = self.bridge['lab_activity']
        _service_activity = self.bridge['service_activity']
        _publication = self.bridge['publication']
        _attached_file = self.bridge['attached_file']
        # create CTEs -  Common Table Expressions for the intensive queries:
        
        collaborators = (
            select([
                _screen_collaborators.c.screen_id,
                _collaborator.c.name,
                _collaborator.c.username,
                _collaborator.c.email,
                _concat(
                    _collaborator.c.name, ' <', _collaborator.c.email, '>'
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
            _ap.c.plate_number,
            func.max(_ap.c.replicate_ordinal).label('max_per_plate') ]).\
                select_from(_ap).\
                group_by(_ap.c.screen_id, _ap.c.plate_number).\
                order_by(_ap.c.screen_id, _ap.c.plate_number)
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
                func.count(distinct(_ap.c.plate_number)).label('count')
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
            
#         affiliation_table = ScreensaverUserResource.get_lab_affiliation_cte()
#         affiliation_table = affiliation_table.cte('la')

        lab_head_table = ScreensaverUserResource.get_lab_head_cte().cte('lab_heads')


        try:
            custom_columns = {
                'collaborator_usernames': (
                    select([
                        func.array_to_string(
                            func.array_agg(collaborators.c.username),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(collaborators)
                    .where(collaborators.c.screen_id 
                        == literal_column('screen.screen_id'))),
                'collaborator_names': (
                    select([
                        func.array_to_string(
                            func.array_agg(collaborators.c.fullname),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(collaborators)
                    .where(collaborators.c.screen_id 
                        == literal_column('screen.screen_id'))),
                'lab_affiliation': lab_head_table.c.lab_affiliation,
                'lab_name': lab_head_table.c.lab_name_full,
                'lab_head_name': lab_head_table.c.name,
                'lab_head_username': lab_head_table.c.username,
                'lead_screener_name': (
                    select([_concat(_su.c.first_name, ' ', _su.c.last_name)])
                    .select_from(_su)
                    .where(_su.c.screensaver_user_id == _screen.c.lead_screener_id)),
                'lead_screener_username': (
                    select([_su.c.username])
                    .select_from(_su)
                    .where(_su.c.screensaver_user_id == _screen.c.lead_screener_id)),
                'has_screen_result': literal_column(
                    '(select dc.data_column_id is not null '
                    '     from data_column dc join screen_result using(screen_result_id) '
                    '     where screen_id=screen.screen_id limit 1 ) '
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
        
        if screens_for_username:
            custom_columns['screensaver_user_role'] = \
                screener_role_cte.c.screensaver_user_role
            
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)

        # build the query statement

        j = _screen
        if screens_for_username:
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

        if screens_for_username:
            stmt = stmt.where(
                screener_role_cte.c.username == screens_for_username)

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
        logger.info('2 - build_list_response')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
         
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
 
        try:
             
            (field_hash, columns, stmt, count_stmt, filename) = self.get_query(schema, param_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
                # use "use_vocab" as a proxy to also adjust siunits for viewing
                rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                    field_hash, rowproxy_generator)
                    
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
                ls.c.username.label('lead_screener_username'),
                pi.c.username.label('pi_username'),
                func.array_agg(collab.c.username).label('collab_usernames')])
            .select_from(j)
            .group_by(_screen.c.facility_id, _screen.c.screen_id,
                ls.c.username, pi.c.username)).cte('screen_associates')
        screener_roles = (
            select([
                _su.c.username,
                sa.c.facility_id,
                sa.c.screen_id,
                case([
                    (_su.c.username == sa.c.lead_screener_username,
                        'lead_screener'),
                    (_su.c.username == sa.c.pi_username,
                        'principal_investigator')
                    ],
                    else_='collaborator').label('screensaver_user_role')
                ])
            .select_from(sa)
            .where(or_(
                # TODO: replace with "any_()" from sqlalchemy 1.1 when avail
                _su.c.username == text(' any(collab_usernames) '),
                _su.c.username == sa.c.lead_screener_username,
                _su.c.username == sa.c.pi_username))
        )
        return screener_roles
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        
        id_kwargs = self.get_id(deserialized, **kwargs)
        Screen.objects.get(**id_kwargs).delete()
    
    def validate(self, _dict, patch=False, current_object=None):
        errors = DbApiResource.validate(self, _dict, patch=patch)
        # if not errors:
        #     errors = {}
        #     dped = _dict.get('data_privacy_expiration_date', None)
        #     dped_notified = _dict.get('data_privacy_expiration_notified_date', None)
        #     min_dped = _dict.get('min_allowed_data_privacy_expiration_date', None)
        #     max_dped = _dict.get('max_allowed_data_privacy_expiration_date', None)
        #     if not dped:
        #         if min_dped or max_dped:
        #             errs['data_privacy_expiration_date'] = \
        #                 'can not be null if min/max dates are set'
        #     if dped_notified:
        #         if min_dped:
        #             errs['min_allowed_data_privacy_expiration_date'] = \
        #                 'can not be set if the expiration notified date is set'
        #         if self.max_allowed_data_privacy_expiration_date:
        #             errs['min_allowed_data_privacy_expiration_date'] = \
        #                 'can not be set if the expiration notified date is set'
        return errors

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
        id_kwargs = self.get_id(deserialized, **kwargs)
        if id_kwargs:
            raise BadRequest(
                'POST may not be used to modify a screen: %r' % id_kwargs )
            
        # find a new facility id
        max_facility_id_sql = '''
            select facility_id::text, project_phase from screen 
            where project_phase='primary_screen' 
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
            logger.info('new screen facility id to be created: %d', max_facility_id+1)
        kwargs['facility_id'] = str(max_facility_id+1)
        
        return self.patch_detail(request,**kwargs)
        
    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        
        id_kwargs = self.get_id(deserialized, validate=True, **kwargs)
        try:
            create = False
            try:
                screen = Screen.objects.get(**id_kwargs)
            except ObjectDoesNotExist, e:
                create = True
                screen = Screen(**id_kwargs)

            initializer_dict = self.parse(deserialized, create=create)
            logger.info('initializer_dict: %r', initializer_dict)
                
            errors = self.validate(deserialized, patch=not create)
            if errors:
                raise ValidationError(errors)
            
            _key = 'lab_head_username'
            if _key in initializer_dict:
                try:
                    initializer_dict['lab_head'] = ScreensaverUser.objects.get(
                        username=initializer_dict[_key])
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key=_key,
                        msg='No such username: %r' % initializer_dict[_key])
            _key = 'lead_screener_username'
            if _key in initializer_dict:
                try:
                    initializer_dict['lead_screener'] = (
                        ScreensaverUser.objects.get(
                            username=initializer_dict[_key]))
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key=_key,
                        msg='No such username: %r' % initializer_dict[_key])
            _key = 'collaborator_usernames'
            if initializer_dict.get(_key, None):
                collaborators = []
                for collaborator_username in initializer_dict[_key]:
                    try:
                        collaborators.append(ScreensaverUser.objects.get(
                            username=collaborator_username))
                    except ObjectDoesNotExist:
                        raise ValidationError(
                            key=_key,
                            msg='No such username: %r' % collaborator_username)
                initializer_dict['collaborators'] = collaborators
            _key = 'pin_transfer_approved_by_username'
            pin_transfer_approved_by = None
            if _key in initializer_dict:
                try:
                    pin_transfer_approved_by = (
                        ScreensaverUser.objects.get(
                            username=initializer_dict[_key]))
                except ObjectDoesNotExist:
                    raise ValidationError(
                        key=_key,
                        msg='No such username: %r' % initializer_dict[_key])
            
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
            
            # NOTE: collaborators cannot be set until after the object is saved:
            # the many-to-many related manager is not functional until then.
            if 'collaborators' in initializer_dict:
                logger.info('collaborators: %r', initializer_dict['collaborators'])
                screen.collaborators = initializer_dict.get('collaborators', None)
            
            logger.info('save/created screen: %r', screen)
            
            # related objects
            
            _key = 'cell_lines'
            _val = related_initializer.get(_key, None)
            if _val:
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
            if _val:
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
            pin_transfer_date_approved = \
                initializer_dict.get('pin_transfer_date_approved',None)
            pin_transfer_comments = \
                initializer_dict.get('pin_transfer_comments', None)

            if pin_transfer_approved_by is not None:
                if screen.pin_transfer_admin_activity is None:
                    activity = \
                        Activity(performed_by=pin_transfer_approved_by)
                    activity.date_of_activity = \
                        activity.date_created
                    activity.save()
                    screen.pin_transfer_admin_activity = activity
                    screen.save()
                    logger.info('created pta: %r', 
                        screen.pin_transfer_admin_activity)
                else:
                    screen.pin_transfer_admin_activity.performed_by = \
                        pin_transfer_approved_by
                    screen.pin_transfer_admin_activity.save()
            if screen.pin_transfer_admin_activity is None:
                # secondary pin transfer validation
                if pin_transfer_date_approved is not None:
                    raise ValidationError(
                        key='pin_transfer_date_approved',
                        msg='requires pin_transfer_approved_by_username')    
                if pin_transfer_comments is not None:
                    raise ValidationError(
                        key='pin_transfer_comments',
                        msg='requires pin_transfer_approved_by_username')    
            else:
                if pin_transfer_date_approved is not None:
                    screen.pin_transfer_admin_activity.date_of_activity = \
                        pin_transfer_date_approved
                if pin_transfer_comments is not None:
                    screen.pin_transfer_admin_activity.comments = \
                        pin_transfer_comments
                screen.pin_transfer_admin_activity.save()
                
            # TODO: determine if this is still used
            _key = 'keywords'
            _val = related_initializer.get(_key, None)
            if _val:
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
            if _val:
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
            
            screen.save()
            logger.info('patch_obj done')
            return { API_RESULT_OBJ: screen }
            
        except Exception, e:
            logger.exception('on patch detail')
            raise e  

class StudyResource(ScreenResource):
    
    class Meta:
        resource_name = 'study'
        max_limit = 10000
        always_return_data = True
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = ScreenAuthorization()
        serializer = ScreenResultSerializer()
        queryset = Screen.objects.all()  # .order_by('facility_id')
        
    def __init__(self, **kwargs):
        super(StudyResource, self).__init__(**kwargs)
            
    def build_schema(self, user=None):
        # Bypass Screen schema
        schema = DbApiResource.build_schema(self, user=user)
        return schema

    
    @read_authorization
    def build_list_response(self, request, **kwargs):
        
        kwargs['project_phase__exact'] = 'annotation'
        
        return super(StudyResource,self).build_list_response(request, **kwargs)


class UserChecklistItemResource(DbApiResource):    

    class Meta:

        queryset = UserChecklistItem.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['']
        resource_name = 'userchecklistitem'
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        
        self.user_resource = None
        super(UserChecklistItemResource, self).__init__(**kwargs)

    def prepend_urls(self):
        
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/(?P<username>([\d\w]+))/" 
                 r"(?P<item_group>([\d\w_]+))/"
                 r"(?P<item_name>([\d\w_]+))%s$")
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]    

    @read_authorization
    def get_detail(self, request, **kwargs):

        username = kwargs.get('username', None)
        if not username:
            raise NotImplementedError('must provide a username parameter')
        item_group = kwargs.get('item_group', None)
        if not item_group:
            raise NotImplementedError('must provide a item_group parameter')
        item_name = kwargs.get('item_name', None)
        if not item_name:
            raise NotImplementedError('must provide a item_name parameter')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
        
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)
        item_group = param_hash.pop('item_group', None)
        if item_group:
            param_hash['item_group__eq'] = item_group
        item_name = param_hash.pop('item_name', None)
        if item_name:
            param_hash['item_name__eq'] = item_name
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)
                  
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
 
            # specific setup
            _su = self.bridge['screensaver_user']
            _admin = _su.alias('admin')
            _up = self.bridge['reports_userprofile']
            _uci = self.bridge['user_checklist_item']
            _vocab = self.bridge['reports_vocabulary']
            
            # get the checklist items & groups
            cig_table = (
                select([
                    _vocab.c.ordinal,
                    _vocab.c.key.label('item_group'),
                    func.array_to_string(
                        array(['checklistitem', _vocab.c.key, 'name']), '.')
                        .label('checklistitemgroup')])
                .select_from(_vocab)
                .where(_vocab.c.scope == 'checklistitem.group'))
            cig_table = Alias(cig_table)
            ci_table = (
                select([
                    _vocab.c.ordinal,
                    cig_table.c.item_group,
                    _vocab.c.key.label('item_name')])
                .select_from(
                    _vocab.join(cig_table,
                        _vocab.c.scope == cig_table.c.checklistitemgroup))
                ).order_by(cig_table.c.ordinal, _vocab.c.ordinal)
            ci_table = ci_table.cte('ci')

            # build the entered checklists
            
            j = _uci
            j = j.join(
                _su, _uci.c.screensaver_user_id == _su.c.screensaver_user_id)
            j = j.join(
                _admin, _uci.c.admin_user_id == _admin.c.screensaver_user_id)
            entered_checklists = select([
                _su.c.username,
                func.array_to_string(
                    array([_su.c.last_name, _su.c.first_name]), ', ')
                    .label('user_fullname'),
                _uci.c.item_group,
                _uci.c.item_name,
                _uci.c.status,
                _uci.c.status_date,
                _admin.c.username.label('admin_username')
                ]).select_from(j)
            username = param_hash.pop('username', None)
            if username:
                entered_checklists = entered_checklists.where(
                    _su.c.username == username)
            entered_checklists = entered_checklists.cte('entered_checklists')
            
            # This entire query doesn't fit the pattern, construct it manually
            # bleah
            custom_columns = {
                'username': func.coalesce(
                    entered_checklists.c.username, username),
                'user_fullname': entered_checklists.c.user_fullname,
                'admin_username': entered_checklists.c.admin_username,
                'item_group': ci_table.c.item_group,
                'item_name' : ci_table.c.item_name,
                'status': func.coalesce(
                    entered_checklists.c.status, 'not_completed'),
                'status_date': entered_checklists.c.status_date
                }

            base_query_tables = ['user_checklist_item', 'screensaver_user'] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            isouter = False
            if username:
                # if username, then this is a user specific view:
                # - outer join in the two so that a full list is generated
                isouter = True
                
            j = ci_table
            j = j.join(
                entered_checklists,
                ci_table.c.item_name == entered_checklists.c.item_name,
                isouter=isouter)
            
            stmt = select(columns.values()).select_from(j)
            if not username:
                stmt = stmt.order_by(entered_checklists.c.username)
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
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
    def delete_obj(self, request, deserialized, **kwargs):
        raise NotImplementedError(
            'delete obj is not implemented for UserChecklistItem')
    
    @write_authorization
    @un_cache
    @transaction.atomic()
    def patch_obj(self, request, deserialized, **kwargs):
        logger.info('patch checklist item: %r', deserialized)
        schema = kwargs['schema']
        fields = schema['fields']
        initializer_dict = {}
        # TODO: wrapper for parsing
        for key in fields.keys():
            if deserialized.get(key, None):
                initializer_dict[key] = parse_val(
                    deserialized.get(key, None), key, fields[key]['data_type']) 
        
        username = deserialized.get('username', None)
        if not username:
            raise Exception('username not specified %s' % deserialized)
        item_group = deserialized.get('item_group', None)
        if not item_group:
            raise Exception('item_group not specified %s' % deserialized)
        item_name = deserialized.get('item_name', None)
        if not item_name:
            raise Exception('item_name not specified %s' % deserialized)
        admin_username = deserialized.get('admin_username', None)
        if not admin_username:
            raise Exception('admin_username not specified %s' % deserialized)

        try:
            user = ScreensaverUser.objects.get(username=username)
            initializer_dict['screensaver_user_id'] = user.pk
        except ObjectDoesNotExist:
            logger.exception('user/username does not exist: %s' % username)
            raise
        try:
            admin_user = ScreensaverUser.objects.get(username=admin_username)
            initializer_dict['admin_user_id'] = admin_user.pk
        except ObjectDoesNotExist:
            logger.exception(
                'admin_user/username does not exist: %s' % admin_username)
            raise Exception(
                'admin_user/username does not exist: %s' % admin_username)

        try:
            try:
                uci = UserChecklistItem.objects.get(
                    screensaver_user=user,
                    item_group=item_group, item_name=item_name)
            except ObjectDoesNotExist:
                logger.info(
                    'UserChecklistItem does not exist: %s/%s/%s, creating' 
                    % (username, item_group, item_name))
                uci = UserChecklistItem()
            for key, val in initializer_dict.items():
                if hasattr(uci, key):
                    # note: setattr does not work for foreign keys
                    setattr(uci, key, val)
                else:
                    logger.warn(
                        'no such attribute on user_checklist_item: %s:%r' 
                        % (key, val))
            
            uci.save()
            return { API_RESULT_OBJ: uci }
        except Exception, e:
            logger.error('on patch_obj')
            raise e


class ScreensaverUserResource(DbApiResource):    

    class Meta:

        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['digested_password']
        resource_name = 'screensaveruser'
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
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/groups%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_groupview'),
                name="api_dispatch_user_groupview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/checklistitems%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_checklistitemview'),
                name="api_dispatch_user_checklistitemview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/attachedfiles%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_attachedfileview'),
                name="api_dispatch_user_attachedfileview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/useragreement%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_useragreement_view'),
                name="api_dispatch_useragreement_view"),
            url((r"^(?P<resource_name>%s)/(?P<username>([\w]+))"
                 r"/attachedfiles/(?P<attached_file_id>([\d]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_attachedfiledetailview'),
                name="api_dispatch_user_attachedfiledetailview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/serviceactivities%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_serviceactivityview'),
                name="api_dispatch_user_serviceactivityview"),
            url((r"^(?P<resource_name>%s)/(?P<username>([\w]+))"
                 r"/serviceactivities/(?P<activity_id>([\d]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_serviceactivitydetailview'),
                name="api_dispatch_user_serviceactivitydetailview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w]+))/screens%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_user_screenview'),
                name="api_dispatch_user_screenview"),
        ]    

    def dispatch_user_groupview(self, request, **kwargs):
        return UserGroupResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_checklistitemview(self, request, **kwargs):
        return UserChecklistItemResource().dispatch('list', request, **kwargs)    
    
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
        kwargs['serviced_username__eq'] = kwargs.pop('username')
        return ServiceActivityResource().dispatch('list', request, **kwargs)    

    def dispatch_user_screenview(self, request, **kwargs):
        kwargs['screens_for_username'] = kwargs.pop('username')
        return ScreenResource().dispatch('list', request, **kwargs)    

    def dispatch_user_serviceactivitydetailview(self, request, **kwargs):
        return ServiceActivityResource().dispatch('detail', request, **kwargs)    
    
    def build_schema(self, user=None):
        
        schema = super(ScreensaverUserResource, self).build_schema(user=user)
        sub_schema = self.get_user_resource().build_schema();
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
    def get_lab_affiliation_cte(cls):
        
        _vocab = get_tables()['reports_vocabulary']
        affiliation_category_table = (
            select([
                _vocab.c.ordinal,
                _vocab.c.key.label('category_key'),
                _vocab.c.title.label('category'),
                func.array_to_string(array(['labaffiliation.category',
                    _vocab.c.key]), '.').label('scope')])
            .select_from(_vocab)
            .where(_vocab.c.scope == 'labaffiliation.category'))
        affiliation_category_table = Alias(affiliation_category_table)
        affiliation_table = (
            select([
                _vocab.c.ordinal,
                affiliation_category_table.c.scope,
                affiliation_category_table.c.category,
                _vocab.c.key.label('affiliation_name'),
                _vocab.c.title])
            .select_from(
                _vocab.join(affiliation_category_table,
                    _vocab.c.scope == affiliation_category_table.c.scope))
            ).order_by(affiliation_category_table.c.ordinal, _vocab.c.ordinal)
        return affiliation_table
    
    @classmethod
    def get_user_cte(cls):
        
        bridge = get_tables()
        _su = bridge['screensaver_user']
        _up = bridge['reports_userprofile']
        _au = bridge['auth_user']
        
        j = _su
        j = j.join(_up, _up.c.id == _su.c.user_id)
        j = j.join(_au, _au.c.id == _up.c.user_id)
        user_table = (
            select([
                _su.c.screensaver_user_id,
                _au.c.username,
                _concat(_au.c.first_name, ' ', _au.c.last_name).label('name'),
                _concat(_au.c.last_name, ', ', _au.c.first_name).label('last_first'),
                _au.c.email
                ])
            .select_from(j))
        return user_table
        
    @classmethod
    def get_lab_head_cte(cls):
        
        affiliation_table = ScreensaverUserResource.get_lab_affiliation_cte().cte('lab_affil')
        _user = ScreensaverUserResource.get_user_cte().cte('lab_head_users')
        bridge = get_tables()
        _su = bridge['screensaver_user']
        
        lab_head_table = (
            select([
                _user.c.screensaver_user_id,
                _user.c.name,
                _user.c.username,
                _concat(
                    affiliation_table.c.title, ' (',
                    affiliation_table.c.category, ')').label('lab_affiliation'),
                func.array_to_string(
                    array([
                        _user.c.last_first,
                        ' - ', affiliation_table.c.title,
                        ' (', affiliation_table.c.category, ')']), '').label('lab_name_full'),
                ])
            .select_from(
                _user.join(
                    _su, _user.c.screensaver_user_id==_su.c.screensaver_user_id)
                .join(affiliation_table, 
                    _su.c.lab_head_affiliation==affiliation_table.c.affiliation_name ))
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
                'no screensaver_user_id, username or ecommons_id provided: %r', kwargs.keys())
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

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        schema = kwargs['schema']
        
        is_for_detail = kwargs.pop('is_for_detail', False)
#         filename = self._get_filename(schema, kwargs)

        try:
            
            # general setup
            
            manual_field_includes = set(param_hash.get('includes', []))
            exact_fields = set(param_hash.get('exact_fields', []))
        
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)
                  
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
 
            # specific setup
            _su = self.bridge['screensaver_user']
            _au = self.bridge['auth_user']
            _up = self.bridge['reports_userprofile']
            _s = self.bridge['screen']
            _screen_collab = self.bridge['screen_collaborators']
            _fur = self.bridge['user_facility_usage_role']
            _lhsu = _su.alias('lhsu')
            
#             affiliation_table = ScreensaverUserResource.get_lab_affiliation_cte().cte('la')
            
            lab_head_table = ScreensaverUserResource.get_lab_head_cte().cte('lab_heads')
            
            custom_columns = {
                'name': literal_column(
                    "auth_user.last_name || ', ' || auth_user.first_name"),
                'screens_lab_head': (
                    select([
                        func.array_to_string(
                            func.array_agg(_s.c.facility_id),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(_s)
                    .where(_s.c.lab_head_id == _su.c.screensaver_user_id)),
                'screens_lead': (
                    select([
                        func.array_to_string(
                            func.array_agg(_s.c.facility_id),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(_s)
                    .where(_s.c.lead_screener_id == _su.c.screensaver_user_id)),
                'screens_collaborator': (
                    select([
                        func.array_to_string(
                            func.array_agg(_s.c.facility_id),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        _s.join(
                            _screen_collab,
                            _s.c.screen_id == _screen_collab.c.screen_id))
                    .where(_screen_collab.c.screensaveruser_id 
                        == _su.c.screensaver_user_id)),
                'lab_name': lab_head_table.c.lab_name_full,
                'lab_head_username': lab_head_table.c.username,
#                 'lab_name': (
#                     select([
#                         func.array_to_string(
#                             array([
#                                 _lhsu.c.last_name, ', ', _lhsu.c.first_name,
#                                 ' - ', affiliation_table.c.title,
#                                 ' (', affiliation_table.c.category, ')']), '')])
#                     .select_from(
#                         _lhsu.join(
#                             affiliation_table,
#                             affiliation_table.c.affiliation_name 
#                                 == _lhsu.c.lab_head_affiliation))
#                     .where(_lhsu.c.screensaver_user_id == _su.c.lab_head_id)),
#                 'lab_head_username': (
#                     select([_lhsu.c.username])
#                     .select_from(_lhsu)
#                     .where(_lhsu.c.screensaver_user_id == _su.c.lab_head_id)),
                'facility_usage_roles': (
                    select([
                        func.array_to_string(
                            func.array_agg(_fur.c.facility_usage_role),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(_fur)
                    .where(_fur.c.screensaver_user_id 
                        == _su.c.screensaver_user_id)),
                }

            # delegate to the user resource
            default_fields = ['fields.screensaveruser', 'fields.user']
            _temp = { key:field for key, field in field_hash.items() 
                if field.get('scope', None) in default_fields }
            field_hash = _temp
            logger.debug('final field hash: %s', field_hash.keys())
            logger.info(
                'TODO: passing screensaver_user fields to reports_userprofile '
                'causes warnings')
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
            j = j.join(_up, _su.c.user_id == _up.c.id)
            j = j.join(_au, _up.c.user_id == _au.c.id)
            
            j = j.join(
                lab_head_table, 
                _su.c.lab_head_id==lab_head_table.c.screensaver_user_id,
                isouter=True)
            
            stmt = select(columns.values()).select_from(j)
            # natural ordering
            stmt = stmt.order_by(_au.c.last_name, _au.c.first_name)
            
            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            logger.debug(
                'stmt: %s',
                str(stmt.compile(
                    dialect=postgresql.dialect(),
                    compile_kwargs={"literal_binds": True})))
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

    def put_detail(self, request, **kwargs):
        raise NotImplementedError('put_list must be implemented')
                
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        
        id_kwargs = self.get_id(deserialized, **kwargs)
        ScreensaverUser.objects.get(**id_kwargs).delete()

    def get_id(self, deserialized, **kwargs):

        # FIXME: this mirrors UserResource.get_id
        # - update the inheritance so that 
        # ScreensaveruserResource extends UserResource
        id_kwargs = DbApiResource.get_id(self, deserialized, **kwargs)
        if not id_kwargs:
            if deserialized and deserialized.get('ecommons_id', None):
                id_kwargs = { 'ecommons_id': deserialized['ecommons_id']}
            elif kwargs and kwargs.get('ecommons_id', None):
                id_kwargs = { 'ecommons_id': kwargs['ecommons_id']}
            else:
                raise ValueError, '%s, nor was an ecommons_id specified' % e
        return id_kwargs
                
    @write_authorization
    @un_cache
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs['schema']
        fields = schema['fields']
        initializer_dict = {}
        # TODO: wrapper for parsing
        logger.debug('fields: %r, deserialized: %r', fields.keys(), deserialized)
        for key in fields.keys():
            if deserialized.get(key, None) is not None:
                initializer_dict[key] = parse_val(
                    deserialized.get(key, None), key, fields[key]['data_type']) 

        id_kwargs = self.get_id(deserialized, **kwargs)

        username = id_kwargs.get('username', None)
        ecommons_id = id_kwargs.get('ecommons_id', None)
        fields = { name:val for name, val in fields.items() 
            if val['scope'] == 'fields.screensaveruser'}
        logger.debug('fields.screensaveruser fields: %s', fields.keys())
        try:
            # create/get userprofile
            patch_response = self.get_user_resource().patch_obj(request, deserialized, **kwargs)
            logger.info('patched userprofile %s', patch_response)
            user = patch_response[API_RESULT_OBJ]

            # create the screensaver_user
            screensaver_user = None
            try:
                screensaver_user = ScreensaverUser.objects.get(user=user)
                logger.info('found user to patch: %r', screensaver_user)
                errors = self.validate(deserialized, patch=True)
                if errors:
                    raise ValidationError(errors)
            except ObjectDoesNotExist:
                logger.info('create user: %r', deserialized)
                errors = self.validate(deserialized, patch=False)
                if errors:
                    raise ValidationError(errors)
                try:
                    if username:
                        screensaver_user = \
                            ScreensaverUser.objects.get(
                                username=username)
                    elif ecommons_id:
                        screensaver_user = \
                            ScreensaverUser.objects.get(
                                user__ecommons_id=ecommons_id)
                    else:
                        raise NotImplementedError(
                            'username or ecommons_id must be specified')
                except ObjectDoesNotExist, e:
                    if not username:
                        logger.info(
                            ('username not specified, '
                            'setting username to ecommons_id: %s'), ecommons_id)
                        username = ecommons_id
                    
                    if hasattr(user, 'screensaveruser'):
                        raise ValueError(
                            'user already exists: %s: %s' 
                            % (user, user.screensaveruser))
                        
                    logger.info(
                        'Screensaver User %s does not exist, creating',
                        username)
                    screensaver_user = \
                        ScreensaverUser.objects.create(username=username)
                    screensaver_user.save()
            initializer_dict = {}
            for key in fields.keys():
                if key in deserialized:
                    initializer_dict[key] = parse_val(
                        deserialized.get(key, None), key,
                        fields[key]['data_type']) 
            if initializer_dict:
                for key, val in initializer_dict.items():
                    if hasattr(screensaver_user, key):
                        setattr(screensaver_user, key, val)
            else:
                logger.info(
                    'no (basic) screensaver_user fields to update %s',
                    deserialized)
            
            screensaver_user.user = user

            # also set legacy screensaveruser fields, temporary convenience
            screensaver_user.username = user.username
            screensaver_user.first_name = user.first_name
            screensaver_user.last_name = user.last_name
            screensaver_user.email = user.email
            screensaver_user.ecommons_id = user.ecommons_id
            # already moved:
            # harvard_id
            # harvard_id_expiration_date
            # harvard_id_requested_expiration_date
            # mailing_address
            
            screensaver_user.save()
            
            if initializer_dict.get('lab_head_username', None):
                lh_username = initializer_dict['lab_head_username']
                if lh_username:
                    try:
                        lab_head = ScreensaverUser.objects.get(
                            username=lh_username)
                        screensaver_user.lab_head = lab_head
                        screensaver_user.save()
                    except ObjectDoesNotExist, e:
                        logger.info(
                            'Lab Head Screensaver User %s does not exist',
                            lh_username)
                        raise BadRequest(
                            'lab_head_username not found %s' % lh_username)
                else:
                    screensaver_user.lab_head = None
                    screensaver_user.save();
                    
            if initializer_dict.get('facility_usage_roles', None):
                current_roles = set([r.facility_usage_role 
                    for r in screensaver_user.userfacilityusagerole_set.all()])
                new_roles = set(initializer_dict['facility_usage_roles'])
                logger.info('roles to delete: %s', current_roles - new_roles)
                (screensaver_user.userfacilityusagerole_set
                    .filter(facility_usage_role__in=current_roles - new_roles)
                    .delete())
                for role in new_roles - current_roles:
                    logger.info(
                        'create f-usage role: %s, %s', screensaver_user, role)
                    ur = UserFacilityUsageRole.objects.create(
                        screensaver_user=screensaver_user,
                        facility_usage_role=role)            
                    ur.save()

            # if 'classification' in initializer_dict:
            #     try:
            #         sru = ScreensaverUser.objects.get(screensaver_user=screensaver_user)
            #     except ObjectDoesNotExist, e:
            #         logger.info('Screening Room User %s does not exist, creating' % username)
            #         sru = ScreensaverUser.objects.create(screensaver_user=screensaver_user)
            #     sru.user_classification = initializer_dict['classification']
            #     sru.save()
            logger.info('patch_obj done')
            return { API_RESULT_OBJ: screensaver_user }
            
        except Exception, e:
            logger.exception('on patch detail')
            raise e  


class NaturalProductReagentResource(DbApiResource):
    # NOTE: OO model for inheritance is broken in multiple places:
    # Consider folding the NaturalProductReagentResource into ReagentResource
    
    class Meta:

        queryset = Reagent.objects.all()
        
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'naturalproductreagent'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(NaturalProductReagentResource, self).__init__(**kwargs)

    def delete_reagents(self, library):
        
        NaturalProductReagent.objects.all().filter(well__library=library).delete()

    def _patch_wells(self, deserialized):
        ''' For bulk update: 
        - deserialized has been loaded with the well
        '''
        schema = self.build_schema()
        for i,well_data in enumerate(deserialized):
            well = well_data['well']
            is_patch = False
            fields = schema['create_fields']
            if not well.reagents.exists():
                reagent = NaturalProductReagent(well=well)
            else:
                is_patch = True
                fields = schema['update_fields']
                # TODO: only works for a single reagent
                # can search for the reagent using id_kwargs
                # reagent = well.reagents.all().filter(**id_kwargs)
                # TODO: update reagent
                reagent = well.reagents.all()[0]
                reagent = reagent.naturalproductreagent
                logger.debug('found reagent: %r, %r', reagent.well_id, reagent)
                
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
            BasicAuthentication(), SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'silencingreagent'
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
                select_stmt = select([func.array_to_string(
                                func.array_agg(column(field_name)),
                                               LIST_DELIMITER_SQL_ARRAY)])
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
            
            # NOT Performant
            if key == 'pool_well':
                _duplex_wells = bridge['silencing_reagent_duplex_wells']
                pool_reagent = bridge['reagent'].alias('pool_reagent')
                 
                columns[key] = (
                    select([pool_reagent.c.well_id])
                        .select_from(_duplex_wells.join(
                            pool_reagent,pool_reagent.c.reagent_id==
                                _duplex_wells.c.silencingreagent_id))
                        .where(_duplex_wells.c.well_id==text('reagent.well_id'))
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

    def _patch_wells(self, deserialized):
        ''' For bulk update: 
        - deserialized has been loaded with the well & duplex wells
        '''
        schema = self.build_schema()
        for i,well_data in enumerate(deserialized):
            well = well_data['well']
            is_patch = False
            fields = schema['create_fields']
            if not well.reagents.exists():
                reagent = SilencingReagent(well=well)
                
                # only allow duplex_wells to be set on create
                if well_data.get('duplex_wells', None):
                    reagent.save()
                    reagent.duplex_wells = well_data['duplex_wells']
            else:
                is_patch = True
                fields = schema['update_fields']
                # TODO: only works for a single reagent
                # can search for the reagent using id_kwargs
                # reagent = well.reagents.all().filter(**id_kwargs)
                # TODO: update reagent
                reagent = well.reagents.all()[0]
                reagent = reagent.silencingreagent
                logger.debug('found reagent: %r, %r', reagent.well_id, reagent)
                
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
            BasicAuthentication(), SessionAuthentication())
        authorization = UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = []  # ['json_field']
        always_return_data = True  # this makes Backbone happy
        resource_name = 'smallmoleculereagent' 

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

    def _patch_wells(self, deserialized):
        ''' For bulk update: 
        - deserialized has been loaded with the wells
        '''
        logger.info('patch reagents...')
        
        schema = self.build_schema()
        for i,well_data in enumerate(deserialized):
            well = well_data['well']
            is_patch = False
            fields = schema['create_fields']
            if not well.reagents.exists():
                reagent = SmallMoleculeReagent(well=well)
                reagent.save()
            else:
                is_patch = True
                fields = schema['update_fields']
                # TODO: only works for a single reagent
                # can search for the reagent using id_kwargs
                # reagent = well.reagents.all().filter(**id_kwargs)
                # TODO: update reagent
                reagent = well.reagents.all()[0]
                reagent = reagent.smallmoleculereagent
                logger.debug('found reagent: %r, %r', reagent.well_id, reagent)
                
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
#         reagent.save()

        self.patch_elapsedtime3 += (time.time() - start_time)
                
        logger.debug('sm patch_obj done')
        return reagent
    

class ReagentResource(DbApiResource):
    
    class Meta:

        queryset = Reagent.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(),
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'reagent'
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
    
    def get_list(self, request, param_hash={}, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def get_query(self, 
        param_hash, library_classification=None, library=None,
        cherry_pick_request_id_screener=None,
        cherry_pick_request_id_lab=None,
        ):
        logger.info('get query for library_classification %r', library_classification )
        schema = self.build_schema(library_classification=library_classification)

        try:
            
            # general setup
            manual_field_includes = set(param_hash.get('includes', []))
   
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)
            
            if filter_expression is None:
                raise InformationError(
                    key='Input filters ',
                    msg='Please enter a filter expression to begin')
                 
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
                    
                j = j.join(
                    pool_wells, pool_wells.c.well_id==_well.c.well_id,
                    isouter=True )
                                
                custom_columns['pool_well'] = pool_wells.c.pool_well_id
            
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
        sub_columns['plate_number'] = (literal_column(
            "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
            .label('plate_number'))
        
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
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)
       
        # TODO: eliminate dependency on library (for schema determination)
        library = None
        library_classification = None
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
        content_type = self.get_accept_content_type(
            request, format=kwargs.get('format', None))
        if content_type == SDF_MIMETYPE:
            manual_field_includes.add('molfile')
            param_hash['includes'] = manual_field_includes
            
        (field_hash, columns, stmt, count_stmt, filename) = \
            self.get_query(
                param_hash, 
                library_classification=library_classification,
                library=library,
                cherry_pick_request_id_screener=cherry_pick_request_id_screener,
                cherry_pick_request_id_lab=cherry_pick_request_id_lab)
        
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
            # use "use_vocab" as a proxy to also adjust siunits for viewing
            rowproxy_generator = DbApiResource.create_siunit_rowproxy_generator(
                field_hash, rowproxy_generator)

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
        if not 'library_short_name' in kwargs:
            return self.build_response(request, self.build_schema(), **kwargs)
        
        library_short_name = kwargs.pop('library_short_name')
        try:
            library = Library.objects.get(short_name=library_short_name)
            return self.build_response(
                request, self.build_schema(library.classification), **kwargs)
            
        except Library.DoesNotExist, e:
            raise Http404(
                'cannot build schema - library def needed'
                'no library found for short_name: %r' % library_short_name)
                
    def build_schema(self, library_classification=None, user=None):
        logger.info('build reagent schema for library_classification: %r',
            library_classification)
        schema = deepcopy(super(ReagentResource, self).build_schema(user=user))
        if library_classification is not None:
            sub_data = self.get_reagent_resource(
                library_classification).build_schema(user=user)
            logger.debug('sub_schema: %r', sub_data['fields'].keys())
            newfields = {}
            newfields.update(sub_data['fields'])
            newfields.update(schema['fields'])
            schema['fields'] = newfields
            
            for k, v in schema.items():
                if k != 'fields' and k in sub_data:
                    schema[k] = sub_data[k]
            
        well_schema = WellResource().build_schema(user=user)
        schema['fields'].update(well_schema['fields'])

        logger.debug('schema: %r', schema)
        return schema

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
                                             SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'well'
        ordering = []
        filtering = {}
        serializer = LimsSerializer()   
#         xls_serializer = XLSSerializer()
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
                .filter(data_column__screen_result__screen__project_phase='annotation')):
            screen = rv.data_column.screen_result.screen
            if screen.facility_id not in screens:
                _screen_data = screens.setdefault(
                    screen.facility_id,
                    { 
                        '1-screen_facility_id': screen.facility_id,
                        'facility_id': screen.facility_id,
                        'title': screen.title,
                        'summary': screen.summary,
                        'lead_screener_username': screen.lead_screener.username,
                        'lead_screener_name': '%s %s' % (
                            screen.lead_screener.first_name, screen.lead_screener.last_name),
                        'lab_head_username': screen.lab_head.username,
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
        content_type = self.get_accept_content_type(
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
                .where(_s.c.project_phase != 'annotation'))
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
                
        content_type = self.get_accept_content_type(
            request,format=kwargs.get('format', None))
        return HttpResponse(
            content=self.get_serializer().serialize(
                data, content_type),
            content_type=content_type)
        
    def get_schema(self, request, **kwargs):
        if not 'library_short_name' in kwargs:
            return self.build_response(request, self.build_schema(),**kwargs)
        
        library_short_name = kwargs.pop('library_short_name')
        try:
            library = Library.objects.get(short_name=library_short_name)
            return self.build_response(
                request, self.build_schema(library.classification),**kwargs)
            
        except Library.DoesNotExist, e:
            raise Http404(
                'cannot build schema - library def needed'
                'no library found for short_name: %r' % library_short_name)
                
    def build_schema(self, library_classification=None, user=None):

        data = super(WellResource, self).build_schema(user=user)
        
        if library_classification:
            sub_data = self.get_reagent_resource().build_schema(
                library_classification=library_classification, user=user)
            newfields = {}
            newfields.update(sub_data['fields'])
            newfields.update(data['fields'])
            data['fields'] = newfields
            
            data['content_types'] = sub_data['content_types']
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
        logger.info('put wells for library: %r', library)
        logger.info('deleting reagents...')
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
         
        if 'library_short_name' not in kwargs:
            raise BadRequest('library_short_name is required')
        library = Library.objects.get(
            short_name=kwargs['library_short_name'])
        logger.info(
            'patch_list: WellResource: library: %r...', library.short_name)
         
        deserialized = self.deserialize(request)
        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
 
        schema = kwargs['schema']
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
        original_data = self._get_list_response(request, **kwargs_for_log)

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
        library_log.uri = self.get_library_resource().get_resource_uri(
            model_to_dict(library))
        library_log.key = (
            '/'.join(
                self.get_library_resource()
                    .get_id(model_to_dict(library)).values()))
        kwargs.update({ 'parent_log': library_log })
 
        logger.info('Cache library wells for patch...') 
        well_map = dict((well.well_id, well) 
            for well in library.well_set.all())
        if len(well_map) == 0:
            raise BadRequest('Library wells have not been created')
         
        logger.info('patch wells, count: %d', len(deserialized))
        # Note: wells can only be created on library creation
        fields = schema['update_fields']
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
            if i % 999 == 0:
                logger.info('patched: %d wells', i+1)
        logger.info('patched %d wells', i+1)

        self.get_reagent_resource().\
            get_reagent_resource(library.classification)._patch_wells(deserialized)        
        
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
        new_data = self._get_list_response(request, **kwargs_for_log)
         
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
                request, {'meta': meta }, response_class=HttpResponse, **kwargs)
        else:
            return self.build_response(
                request,  {'meta': meta }, response_class=HttpResponse, **kwargs)

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
        
        # logger.debug('patch: %r', deserialized)
        # library = kwargs.get('library', None)
        # if not library:
        #     library_short_name = kwargs.get('library_short_name', None)
        #     if not library_short_name:
        #         raise ValidationError(key='library_short_name', msg='required')
        #     try:
        #         library = Library.objects.get(short_name=library_short_name)
        #         kwargs['library'] = library
        #     except ObjectDoesNotExist:
        #         raise Http404('library not found: %r' % library_short_name)
        # 
        # initializer_dict = self.parse(deserialized)
        # 
        # id_kwargs = self.get_id(deserialized, **kwargs)
        # 
        # well = kwargs.get('well', None)
        # if not well:
        #     # find the well, to allow for patching
        #     try:
        #         well = Well.obj.get(**id_kwargs)
        #         kwargs['well'] = well
        #     except ObjectDoesNotExist:
        #         raise Http404('well not found: %r' % id_kwargs)
        # 
        # errors = self.validate(initializer_dict, patch=True)
        # if errors:
        #     raise ValidationError(errors)
        # 
        # for key, val in initializer_dict.items():
        #     if hasattr(well, key):
        #         setattr(well, key, val)
        # 
        # well.save()
        # 
        # duplex_wells = []
        # if deserialized.get('duplex_wells', None):
        #     if not library.is_pool:
        #         raise ValidationError(
        #             key='duplex_wells',
        #             msg='library is not a pool libary: %r' % library.short_name)
        #     well_ids = deserialized['duplex_wells']  # .split(';')
        #     for well_id in well_ids:
        #         try:
        #             duplex_wells.append(Well.objects.get(well_id=well_id))
        #         except:
        #             raise ValidationError(
        #                 key='duplex_well not found',
        #                 msg='well: %r, pool well: %r' % (well.well_id, well_id))
        #     kwargs['duplex_wells'] = duplex_wells
        # # lookup/create the reagent
        # # TODO: delegate this to the ReagentResource
        # self.get_reagent_resource().patch_obj(request, deserialized, **kwargs)
        # 
        # return well

                
class LibraryResource(DbApiResource):
    
    class Meta:

        queryset = Library.objects.all()  # .order_by('facility_id')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization = UserGroupAuthorization()
        resource_name = 'library'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        always_return_data = True 
        
    def __init__(self, **kwargs):
        
        self.well_resource = None
        self.apilog_resource = None
        self.reagent_resource = None
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
#             schema = super(LibraryResource, self).build_schema()
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        for_screen_id = param_hash.pop('for_screen_id', None)
#         filename = self._get_filename(schema, kwargs)

        try:
            # general setup
            
            manual_field_includes = set(param_hash.get('includes', []))
            
            if is_for_detail:
                manual_field_includes.add('concentration_types')
 
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(readable_filter_hash)
                
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

            # specific setup
            
            _apilog = self.bridge['reports_apilog']
            _logdiff = self.bridge['reports_logdiff']
            _l = self.bridge['library']
                                     
            custom_columns = {
#                 'comments': (
#                     select([
#                         func.array_to_string(func.array_agg(
#                             literal_column('comment')), LIST_DELIMITER_SQL_ARRAY)])
#                         .select_from(_apilog)
#                         .where(_apilog.c.ref_resource_name=='library')
#                         .where(_apilog.c.key==literal_column('library.short_name'))
#                         .where(_apilog.c.diffs==None)
#                         .where(_apilog.c.comment!=None)
#                         ),
                'comments': (
                    select([
                        func.array_to_string(func.array_agg(
                            literal_column('comment')), LIST_DELIMITER_SQL_ARRAY)])
                        .select_from(_apilog)
                        .where(_apilog.c.ref_resource_name=='library')
                        .where(_apilog.c.key==literal_column('library.short_name'))
                        # .where(not_(exists(
                        #     select([None]).select_from(_logdiff)
                        #     .where(_logdiff.c.log_id==_apilog.c.id))))
                        .where(_apilog.c.comment!=None)
                        ),
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
            concentration_sql = (
                'select short_name, '
                'exists(select null from well '
                '  where well.library_id=l.library_id '
                '  and well.mg_ml_concentration is not null limit 1) mg_ml_measure, '
                'exists(select null from well '
                '  where well.library_id=l.library_id '
                '  and well.molar_concentration is not null limit 1) molar_measure '
                'from library l;')
            with get_engine().connect() as conn:
                result = conn.execute(text(concentration_sql))
                cached_ids =  [x[0] for x in result ]
            
                     
            base_query_tables = ['library']

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns)

            # build the query statement

            j = _l
            stmt = select(columns.values()).select_from(j)

            if for_screen_id:
                stmt = stmt.where(_l.c.library_id.in_(
                    self.get_screen_library_ids(for_screen_id)))

            # general setup
             
            (stmt, count_stmt) = self.wrap_statement(
                stmt, order_clauses, filter_expression)
            
            if not order_clauses:
                stmt = stmt.order_by("short_name")
            
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
    def get_screen_library_ids(cls, for_screen_id):
        
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
                .where(_screen.c.facility_id == for_screen_id))
            library_ids = [x[0] for x in 
                conn.execute(query)]
            return library_ids
       
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        
        id_kwargs = self.get_id(deserialized, **kwargs)
        Library.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):
        logger.info('patch library: %r', deserialized)
        id_kwargs = self.get_id(deserialized, validate=True, **kwargs)
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

            initializer_dict = self.parse(deserialized, create=create)
            errors = self.validate(initializer_dict, patch=not create)
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
                    logger.info('bulk create wells: %s-%s', 
                        library.start_plate, library.end_plate)
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
    
    def _build_resources(self):
        
        resources = None
        if self.use_cache:
            resources = cache.get('dbresources')
        if not resources:

            resources = super(ResourceResource, self)._build_resources(use_cache=False)
        
            for key,resource in resources.items():
                self.extend_resource_specific_data(resource)
                
            cache.set('dbresources', resources)
    
        return resources
    
    def extend_resource_specific_data(self, resource_data):
        
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
        elif key == 'labcherrypick':
            resource_data['extraSelectorOptions'] = {
                'label': 'Status',
                'searchColumn': 'status',
                'options': ['unfulfilled','selected','plated','not_selected']}
            
