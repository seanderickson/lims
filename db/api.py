
import StringIO
from __builtin__ import StopIteration
import cStringIO
from collections import defaultdict, OrderedDict
from copy import deepcopy
import csv
import hashlib
import io
import json
import logging
import math
import os.path
import re
import shutil
import sys
import time
from zipfile import ZipFile

from PIL import Image
from django.conf import settings
from django.conf.urls import url
from django.contrib.auth.models import User
from django.core.cache import cache
from django.core.exceptions import ObjectDoesNotExist, MultipleObjectsReturned
from django.core.serializers.json import DjangoJSONEncoder
from django.db import connection
from django.db import transaction
from django.db.models.aggregates import Max, Min
import django.db.models.constants
import django.db.models.sql.constants
from django.forms.models import model_to_dict
from django.http import Http404
from django.http.response import StreamingHttpResponse, HttpResponse
from sqlalchemy import select, asc, text
import sqlalchemy
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.sql import and_, or_, not_          
from sqlalchemy.sql import asc, desc, alias, Alias
from sqlalchemy.sql import func
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join, insert, delete, distinct, \
    exists
from sqlalchemy.sql.expression import nullsfirst, nullslast
from tastypie import fields
from tastypie import http
from tastypie.authentication import BasicAuthentication, SessionAuthentication, \
    MultiAuthentication
from tastypie.authorization import Authorization
from tastypie.bundle import Bundle
from tastypie.constants import ALL_WITH_RELATIONS
from tastypie.exceptions import BadRequest, ImmediateHttpResponse, \
    UnsupportedFormat, NotFound
from tastypie.http import HttpNotFound
from tastypie.resources import Resource
from tastypie.utils import timezone
from tastypie.utils.dict import dict_strip_unicode_keys
from tastypie.utils.urls import trailing_slash
from tastypie.validation import Validation

from db.models import ScreensaverUser, Screen, LabHead, LabAffiliation, \
    ScreeningRoomUser, ScreenResult, DataColumn, Library, Plate, Copy, \
    CopyWell, \
    PlateLocation, Reagent, Well, LibraryContentsVersion, Activity, \
    AdministrativeActivity, SmallMoleculeReagent, SilencingReagent, GeneSymbol, \
    NaturalProductReagent, Molfile, Gene, GeneGenbankAccessionNumber, \
    CherryPickRequest, CherryPickAssayPlate, CherryPickLiquidTransfer, \
    CachedQuery, ChecklistItemEvent, UserChecklistItem
from db.support import lims_utils
from db.support.data_converter import default_converter
from reports import LIST_DELIMITER_SQL_ARRAY, \
    HTTP_PARAM_USE_TITLES, HTTP_PARAM_USE_VOCAB, HEADER_APILOG_COMMENT
from reports.api import ManagedModelResource, ManagedResource, ApiLogResource, \
    UserGroupAuthorization, ManagedLinkedResource, log_obj_update, \
    UnlimitedDownloadResource, IccblBaseResource, VocabulariesResource, \
    MetaHashResource, UserResource, compare_dicts, parse_val, ManagedSqlAlchemyResourceMixin, \
    UserGroupResource
from reports.models import MetaHash, Vocabularies, ApiLog, UserProfile
from reports.serializers import CursorSerializer, LimsSerializer, XLSSerializer
from reports.sqlalchemy_resource import SqlAlchemyResource, un_cache
from reports.utils.sqlalchemy_bridge import Bridge

PLATE_NUMBER_SQL_FORMAT = 'FM9900000'

logger = logging.getLogger(__name__)
    
def _get_raw_time_string():
  return timezone.now().strftime("%Y%m%d%H%M%S")
    

class ScreeningRoomUserResource(ManagedModelResource):
    screensaver_user = fields.ToOneField(
        'db.api.ScreensaverUserResource', attribute='screensaver_user', 
        full=True, full_detail=True, full_list=False)
    class Meta:
        queryset = ScreeningRoomUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
            SessionAuthentication())
        authorization= UserGroupAuthorization()
    
class LabAffiliationResource(ManagedModelResource):   
    class Meta:
        queryset = LabAffiliation.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
    
class LabHeadResource(ManagedModelResource):

    screens = fields.ToManyField('db.api.ScreenResource', 'screens', 
        related_name='lab_head', blank=True, null=True)

    lab_affiliation = fields.ToOneField('db.api.LabAffiliationResource', 
        attribute='lab_affiliation',  full=True, null=True)
    
    # rather than walk the inheritance hierarchy, will flatten this hierarchy 
    # in the dehydrate method
    #    screening_room_user = fields.ToOneField('db.api.ScreeningRoomUserResource', 
    #        attribute='screensaver_user',  full=True)
    
    id = fields.IntegerField(attribute='screensaver_user_id')
    
    class Meta:
        queryset = LabHead.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        
    def dehydrate(self, bundle):
        # flatten the inheritance hierarchy, rather than show nested
        # "lab_head->screening_room_user->screensaver_user"
        bundle.data.update(model_to_dict(bundle.obj.screensaver_user))
        bundle.data.update(model_to_dict(
            bundle.obj.screensaver_user.screensaver_user))
        bundle.data['screens'] = [
            model_to_dict(x) 
            for x in Screen.objects.filter(
                lab_head_id=bundle.obj.screensaver_user.screensaver_user_id)]
        
        return bundle        
    

class PlateLocationResource(ManagedModelResource):

    class Meta:
        queryset = PlateLocation.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'platelocation'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

        
    def __init__(self, **kwargs):
        super(PlateLocationResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<plate_id>((?=(schema))__|(?!(schema))[^/]+))%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]    
        
class LibraryCopyPlateResource(SqlAlchemyResource,ManagedModelResource):

    class Meta:
        queryset = Plate.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'librarycopyplate'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

    def __init__(self, **kwargs):
        super(LibraryCopyPlateResource,self).__init__(**kwargs)

    def prepend_urls(self):

        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('search'), name="api_search"),

            url(r"^(?P<resource_name>%s)/(?P<library_short_name>[\w\d_.\-\+: ]+)"
                r"/(?P<copy_name>[\w\d_.\-\+: ]+)"
                r"/(?P<plate_number>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<copy_name>[\w\d_.\-\+: ]+)"
                r"/(?P<plate_number>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    def get_detail(self, request, **kwargs):
        logger.info(str(('get_detail')))

        library_short_name = kwargs.get('library_short_name', None)
        if not library_short_name:
            logger.info(str(('no library_short_name provided')))
        
        copy_name = kwargs.get('copy_name', None)
        if not copy_name:
            logger.info(str(('no copy_name provided')))
            raise NotImplementedError('must provide a copy_name parameter')
        
        plate_number = kwargs.get('plate_number', None)
        if not copy_name:
            logger.info(str(('no plate_number provided')))
            raise NotImplementedError('must provide a plate_number parameter')
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)

    def get_list(self, request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)
        
    def build_list_response(self,request, param_hash={}, **kwargs):
            
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)
        
        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = super(LibraryCopyPlateResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)

        library_short_name = param_hash.pop('library_short_name', 
            param_hash.get('library_short_name__eq',None))
        if not library_short_name:
            logger.info(str(('no library_short_name provided')))
        else:
            param_hash['library_short_name__eq'] = library_short_name

        copy_name = param_hash.pop('copy_name', 
            param_hash.get('copy_name', None))
        if copy_name:
            param_hash['copy_name__eq'] = copy_name
            
        plate_number = param_hash.pop('plate_number', 
            param_hash.get('plate_number', None))
        if plate_number:
            param_hash['plate_number__eq'] = plate_number
            
        logger.info(str(('get_list', filename, param_hash)))
 
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)

            if filter_expression is None:
                msgs = { 'Library copy plates resource': 'can only service requests with filter expressions' }
                logger.info(str((msgs)))
                raise ImmediateHttpResponse(response=self.error_response(request,msgs))
                 
                 
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
 
            custom_columns={
                'screening_count': literal_column('p1.screening_count'),
                'plate_number': ( literal_column(
                    "to_char(plate.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
                    .label('plate_number') ), 
                'ap_count': literal_column('p1.ap_count'), 
                'dl_count':literal_column('p1.dl_count'),
                'first_date_data_loaded':literal_column('p1.first_date_data_loaded'), 
                'last_date_data_loaded':literal_column('p1.last_date_data_loaded'), 
                'first_date_screened':literal_column('p1.first_date_screened'), 
                'last_date_screened':literal_column('p1.last_date_screened'), 
                'status_date': literal_column(
                    '(select date_of_activity'
                    ' from activity a'
                    ' join administrative_activity aa on(aa.activity_id=a.activity_id) '
                    ' join plate_update_activity pu on(a.activity_id=pu.update_activity_id)'
                    ' where pu.plate_id = plate.plate_id'
                    " and aa.administrative_activity_type='Plate Status Update' "
                    ' order by date_created desc limit 1 )').label('status_date'),
                'status_performed_by': literal_column(
                    "(select su.first_name || ' ' || su.last_name "
                    ' from activity a'
                    ' join screensaver_user su on(a.performed_by_id=su.screensaver_user_id) '
                    ' join administrative_activity aa on(aa.activity_id=a.activity_id) '
                    ' join plate_update_activity pu on(a.activity_id=pu.update_activity_id)'
                    ' where pu.plate_id = plate.plate_id'
                    " and aa.administrative_activity_type='Plate Status Update' "
                    ' order by a.date_created desc limit 1 )').label('status_performed_by'),
                'status_performed_by_id': literal_column(
                    "(select su.screensaver_user_id "     # TODO: replace with the final login id
                    ' from activity a'
                    ' join screensaver_user su on(a.performed_by_id=su.screensaver_user_id) '
                    ' join administrative_activity aa on(aa.activity_id=a.activity_id) '
                    ' join plate_update_activity pu on(a.activity_id=pu.update_activity_id)'
                    ' where pu.plate_id = plate.plate_id'
                    " and aa.administrative_activity_type='Plate Status Update' "
                    ' order by a.date_created desc limit 1 )').label('status_performed_by_id'),
                'date_plated': literal_column(
                    '(select date_of_activity '
                    ' from activity a'
                    ' where a.activity_id=plate.plated_activity_id )').label('date_plated'),
                'date_retired': literal_column(
                    '(select date_of_activity '
                    ' from activity a'
                    ' where a.activity_id=plate.retired_activity_id )').label('date_retired'),
                    };

            base_query_tables = ['plate', 'copy','plate_location', 'library']

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )
            logger.info(str(('columns', columns)))
            # build the query statement

            _p = self.bridge['plate']
            _pl = self.bridge['plate_location']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _ap = self.bridge['assay_plate']

#             # NOTE: precalculated version
#             plate_screening_statistics = \
#                 select([text('*')]).\
#                     select_from(text('plate_screening_statistics')).cte('plate_screening_statistics')
#         
#             p1 = plate_screening_statistics.alias('p1')
            p1 = self.bridge['plate_screening_statistics']
            p1 = p1.alias('p1')
            j = join(_p, _c, _p.c.copy_id == _c.c.copy_id )
            j = j.join(p1, _p.c.plate_id == text('p1.plate_id'), isouter=True)
            j = j.join(_pl, _p.c.plate_location_id == _pl.c.plate_location_id )
            j = j.join(_l, _c.c.library_id == _l.c.library_id )

            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
 
            if not order_clauses:
                stmt = stmt.order_by("plate_number","copy_name")

            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']

            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, param_hash=param_hash, is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
            
                        
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e   
    
    def build_schema(self):
        schema = cache.get(self._meta.resource_name + ":schema")
        if not schema:
            # FIXME: these options should be defined automatically from a vocabulary in build_schema
            schema = super(LibraryCopyPlateResource,self).build_schema()
            temp = [ x.status for x in self.Meta.queryset.distinct('status')]
            schema['extraSelectorOptions'] = { 
                'label': 'Type', 'searchColumn': 'status', 'options': temp }
        return schema
    
    def obj_create(self, bundle, **kwargs):
        bundle.data['date_created'] = timezone.now()
        
        bundle.data['version'] = 1
        logger.info(str(('===creating library copy plate', bundle.data)))

        return super(LibraryCopyPlateResource, self).obj_create(bundle, **kwargs)

 
class NaturalProductReagentResource(ManagedLinkedResource):
    
    class Meta:

        queryset = Reagent.objects.all()
        
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'naturalproductreagent'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(NaturalProductReagentResource,self).__init__(**kwargs)


class SilencingReagentResource(ManagedLinkedResource):
    reagent_id = fields.IntegerField(default=None)
    class Meta:
        queryset = Reagent.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'silencingreagent'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(SilencingReagentResource,self).__init__(**kwargs)

    def build_sqlalchemy_columns(self, fields, bridge):
        '''
        returns an array of sqlalchemy.sql.schema.Column objects, associated 
        with the sqlalchemy.sql.schema.Table definitions, which are bound to 
        the sqlalchemy.engine.Engine which: 
        "Connects a Pool and Dialect together to provide a source of database 
        connectivity and behavior."
        
        @param fields - field definitions, from the resource schema
        
        '''
        DEBUG_BUILD_COLS = False or logger.isEnabledFor(logging.DEBUG)
        
        columns = {}
        vendor_gene_columns=['vendor_entrezgene_id',
            'vendor_gene_name','vendor_gene_species']
        vendor_gene_symbols = 'vendor_entrezgene_symbols'
        vendor_genebank_accession_numbers = 'vendor_genbank_accession_numbers'
        facility_gene_columns=['facility_entrezgene_id',
            'facility_gene_name','facility_gene_species']
        facility_gene_symbols = 'facility_entrezgene_symbols'
        facility_genebank_accession_numbers = 'facility_genbank_accession_numbers'
        
        duplex_wells = 'duplex_wells'
        
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
            label = field['key']
            if DEBUG_BUILD_COLS: 
                logger.info(str(('field[key]', field['key'])))
            join_stmt = None
            join_column = None
            if field['key'] in vendor_columns:
                join_column = 'vendor_gene_id'
            if field['key'] in facility_columns:
                join_column = 'facility_gene_id'

            if field['key'] in vendor_gene_columns or \
                    field['key'] in facility_gene_columns:
                join_stmt = gene_table.join(sirna_table, 
                    gene_table.c['gene_id'] == sirna_table.c[join_column])
                select_stmt = select([gene_table.c[field_name]]).\
                    select_from(join_stmt)
                select_stmt = select_stmt.where(
                    text('silencing_reagent.reagent_id=reagent.reagent_id'))
                select_stmt = select_stmt.label(label)
                columns[label] = select_stmt

            if field['key'] == vendor_gene_symbols or \
                    field['key'] == facility_gene_symbols:
                join_stmt = gene_symbol.join(gene_table, 
                    gene_symbol.c['gene_id'] == gene_table.c['gene_id'])
                join_stmt = join_stmt.join(sirna_table, 
                    gene_table.c['gene_id'] == sirna_table.c[join_column])
                
                select_inner = select([gene_symbol.c[field_name]]).\
                    select_from(join_stmt)
                ordinal_field = field.get('ordinal_field', None)
                if ordinal_field:
                    select_inner = select_inner.order_by(gene_symbol.c[ordinal_field])
                select_inner = select_inner.where(
                    text('silencing_reagent.reagent_id=reagent.reagent_id'))
                select_inner = select_inner.alias('a')
                select_stmt = select([func.array_to_string(
                                func.array_agg(column(field_name)),
                                               LIST_DELIMITER_SQL_ARRAY)])
                select_stmt = select_stmt.select_from(select_inner)
                select_stmt = select_stmt.label(label)
                columns[label] = select_stmt

            if field['key'] == vendor_genebank_accession_numbers or \
                    field['key'] == facility_genebank_accession_numbers:
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
                select_stmt = select_stmt.label(label)
                columns[label] = select_stmt
            
            if field['key'] == duplex_wells:
                duplex_wells = bridge['silencing_reagent_duplex_wells']
                
#                 join_stmt = duplex_wells.join(sirna_table, 
#                     duplex_wells.c['silencingreagent_id'] == sirna_table.c['reagent_id'])
                select_inner = select([duplex_wells.c['well_id']]).\
                    select_from(duplex_wells)
                select_inner = select_inner.where(
                    text('silencingreagent_id=reagent.reagent_id'))
                select_inner = select_inner.order_by(duplex_wells.c['well_id'])
                select_inner = select_inner.alias('a')

                select_stmt = select([func.array_to_string(
                                func.array_agg(column(field_name)),
                                               LIST_DELIMITER_SQL_ARRAY)])
                select_stmt = select_stmt.select_from(select_inner)
                select_stmt = select_stmt.label(label)
                columns[label] = select_stmt
                
                if DEBUG_BUILD_COLS:
                    logger.info(str((select_stmt)))
                
                
        if DEBUG_BUILD_COLS: 
            logger.info(str(('sirna columns', columns.keys())))
        
        return columns 

    def obj_create(self, bundle, **kwargs):
        
        bundle = super(SilencingReagentResource, self).obj_create(bundle, **kwargs)
        
        if 'duplex_wells' in kwargs:
            bundle.obj.silencingreagent.duplex_wells = kwargs['duplex_wells']
        
        # Now do the gene tables
        ## nastiness ensues!
        
        gene_key = 'entrezgene_id'
        if bundle.data.get('vendor_%s'%gene_key, None):
            bundle.obj.silencingreagent.vendor_gene = \
                self._create_gene(bundle.data, 'vendor')
        if bundle.data.get('facility_%s'%gene_key, None):
            bundle.obj.silencingreagent.facility_gene = \
                self._create_gene(bundle.data, 'facility')
        bundle.obj.silencingreagent.save()
        
        return bundle
    
    def _create_gene(self, data, source_type):
        
        gene_keys = ['entrezgene_id', 'gene_name', 'species_name']
        gene = Gene()
        for key in gene_keys:
            api_key = '%s_%s' % (source_type,key)
            val = data.get(api_key, None)
            if val:
                setattr(gene,key,val)
        gene.save()
        
        _key = 'entrezgene_symbols'
        if data.get('%s_%s' % (source_type,_key), None):
            symbol_list = data['%s_%s' % (source_type,_key)] #.split(';')
            for i,symbol in enumerate(symbol_list):
                gene_symbol = GeneSymbol()
                setattr(gene_symbol, 'entrezgene_symbol', symbol)
                setattr(gene_symbol, 'ordinal', i)
                setattr(gene_symbol, 'gene', gene)
                gene_symbol.save()
    
        _key = 'genbank_accession_numbers'
        if data.get('%s_%s' % (source_type,_key), None):
            _list = data['%s_%s' % (source_type,_key)] #.split(';')
            for i,num in enumerate(_list):
                accession_number = GeneGenbankAccessionNumber()
                setattr(accession_number, 'genbank_accession_number', num)
                setattr(accession_number, 'gene', gene)
                accession_number.save()
        
        return gene
    
    def dehydrate(self, bundle):
        
        bundle = super(SilencingReagentResource, self).dehydrate(bundle)
        
        if bundle.obj and hasattr(bundle.obj,'silencingreagent'):
            if bundle.obj.silencingreagent.vendor_gene:
                gene = bundle.obj.silencingreagent.vendor_gene
                type = 'vendor'
                self._dehydrate_gene(gene, type, bundle)
            
            if bundle.obj.silencingreagent.facility_gene:
                gene = bundle.obj.silencingreagent.facility_gene
                type = 'facility'
                self._dehydrate_gene(gene, type, bundle)
            
            if bundle.obj.silencingreagent.duplex_wells.exists():
                bundle.data['duplex_wells'] = ';'.join(
                    [x.well_id for x in bundle.obj.silencingreagent.duplex_wells.all().order_by('well_id') ])
        return bundle
        
    def _dehydrate_gene(self, gene, type, bundle):
        
        gene_keys = ['entrezgene_id', 'gene_name', 'species_name']

        for key in gene_keys:
            bundle.data['%s_%s' %(type,key)] = getattr(gene, key)
        
        _key = 'entrezgene_symbols'
        if gene.genesymbol_set.exists():
            bundle.data['%s_%s'%(type,_key)] = ';'.join(
                [x.entrezgene_symbol for x in gene.genesymbol_set.all().order_by('ordinal')])
        _key = 'genbank_accession_numbers'
        if gene.genegenbankaccessionnumber_set.exists():
            bundle.data['%s_%s'%(type,_key)] = ';'.join(
                [x.genbank_accession_number for x in gene.genegenbankaccessionnumber_set.all()])
        

class SmallMoleculeReagentResource(ManagedLinkedResource):
        
    class Meta:
        queryset = Reagent.objects.all() 
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= UserGroupAuthorization()

        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='smallmoleculereagent' 

    def __init__(self, **kwargs):
        super(SmallMoleculeReagentResource,self).__init__(**kwargs)


class ScreenResource(SqlAlchemyResource,ManagedModelResource):
    
    class Meta:
        queryset = Screen.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'screen'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

        
    def __init__(self, **kwargs):
#        self.
        super(ScreenResource,self).__init__(**kwargs)

    def prepend_urls(self):
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>([\w\d_]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

    def get_detail(self, request, **kwargs):
        logger.info(str(('get_detail')))

        facility_id = kwargs.get('facility_id', None)
        if not facility_id:
            logger.info(str(('no facility_id provided')))
            raise NotImplementedError('must provide a facility_id parameter')
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)
        
    
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

        
    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)
             
        schema = super(ScreenResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)
        
        facility_id = param_hash.pop('facility_id', None)
        if facility_id:
            param_hash['facility_id__eq'] = facility_id

        screen_type = param_hash.get('screen_type__eq', None)

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
#             if screen_type:
#                 if screen_type == 'rnai':
#                     manual_field_includes.append('fields.silencingreagent')
#                 elif screen_type == 'natural_products':
#                     manual_field_includes.append('fields.naturalproductreagent')
#                 else:
#                     manual_field_includes.append('fields.smallmoleculereagent')
            
            
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)
                  
            visibilities = set()                
            if is_for_detail:
                visibilities.update(['detail','summary'])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail, visibilities=visibilities)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
            base_query_tables = ['screen','screen_result']
            _screen = self.bridge['screen']
            _screen_result = self.bridge['screen_result']
            _ap = self.bridge['assay_plate']
            _library = self.bridge['library']
            _copy = self.bridge['copy']
            _plate = self.bridge['plate']
            _cpr = self.bridge['cherry_pick_request']
            _lcp = self.bridge['lab_cherry_pick']
            _cpap = self.bridge['cherry_pick_assay_plate']
            _cplt = self.bridge['cherry_pick_liquid_transfer']

            # create CTEs -  Common Table Expressions for the intensive queries:
            
            # create a cte for the max screened replicates_per_assay_plate
            # - cross join version:
            # select ap.screen_id, ap.plate_number, ap.replicate_ordinal,lesser.replicate_ordinal  
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

            apdl = select([_ap.c.screen_id, func.count(1).label('count')]).\
                select_from(_ap).\
                where(_ap.c.screen_result_data_loading_id != None ).\
                group_by(_ap.c.screen_id).\
                order_by(_ap.c.screen_id)
            apdl = apdl.cte('apdl')
            
            # create a cte for the max screened replicates_per_assay_plate
            apsrc = select([
                _ap.c.screen_id,
                _ap.c.plate_number,
                func.max(_ap.c.replicate_ordinal).label('max_per_plate') ]).\
                    select_from(_ap).\
                    group_by(_ap.c.screen_id, _ap.c.plate_number).\
                    order_by(_ap.c.screen_id, _ap.c.plate_number)
            apsrc = apsrc.cte('apsrc')
            
            # similarly, create a cte for the max data loaded replicates per assay_plate
            apdlrc = select([
                _ap.c.screen_id,
                _ap.c.plate_number,
                func.max(_ap.c.replicate_ordinal).label('max_per_plate') ]).\
                    select_from(_ap).\
                    where(_ap.c.screen_result_data_loading_id != None ).\
                    group_by(_ap.c.screen_id, _ap.c.plate_number).\
                    order_by(_ap.c.screen_id, _ap.c.plate_number)
            apdlrc = apdlrc.cte('apdlrc')
            
            lps = select([
                _ap.c.screen_id,
                func.count(distinct(_ap.c.plate_number)).label('count')]).\
                    select_from(_ap).\
                    where(_ap.c.library_screening_id != None).\
                    group_by(_ap.c.screen_id).cte('lps')
            lpdl = select([
                _ap.c.screen_id,
                func.count(distinct(_ap.c.plate_number)).label('count')]).\
                    select_from(_ap).\
                    where(_ap.c.screen_result_data_loading_id != None).\
                    group_by(_ap.c.screen_id).cte('lpdl')

            libraries_screened = select([
                func.count(distinct(_library.c.library_id)).label('count'),
                _ap.c.screen_id]).\
                    select_from(
                        _ap.join(_plate,_ap.c.plate_id==_plate.c.plate_id).\
                        join(_copy, _plate.c.copy_id==_copy.c.copy_id).\
                        join(_library,_copy.c.library_id==_library.c.library_id)).\
                    group_by(_ap.c.screen_id).cte('libraries_screened')
                    
            tplcps = select([
                _cpr.c.screen_id,
                func.count(1).label('count')]).\
                    select_from(_cpr.join(
                        _lcp,_cpr.c.cherry_pick_request_id==_lcp.c.cherry_pick_request_id).\
                        join(_cpap, _lcp.c.cherry_pick_assay_plate_id==_cpap.c.cherry_pick_assay_plate_id).\
                        join(_cplt,_cplt.c.activity_id==_cpap.c.cherry_pick_liquid_transfer_id)).\
                    group_by(_cpr.c.screen_id).\
                    where(_cplt.c.status == 'Successful' ).cte('tplcps')
                        
            custom_columns = {
                'lab_head_name': literal_column(
                    '( select su.first_name || $$ $$ || su.last_name'
                    '  from screensaver_user su '
                    '  where su.screensaver_user_id=screen.lab_head_id )'
                    ),
                'lead_screener_name': literal_column(
                    '( select su.first_name || $$ $$ || su.last_name'
                    '  from screensaver_user su '
                    '  where su.screensaver_user_id=screen.lead_screener_id )'
                    ),
                'lab_affiliation': literal_column(
                    '( select lf.affiliation_name '
                    '  from lab_affiliation lf '
                    '  join lab_head using(lab_affiliation_id) '
                    '  where lab_head.screensaver_user_id=screen.lab_head_id )'
                    ).label('lab_affiliation'),
                'has_screen_result': literal_column(
                    '(select exists(select null from screen_result '
                    '     where screen_id=screen.screen_id ) ) '
                    ),
                # TODO: convert to vocabulary
                'cell_lines': literal_column("'tbd'"),
#                 'cell_lines': literal_column(
#                     "(select array_to_string(array_agg(c1.value),'%s') "
#                     '    from ( select c.value from cell_line c '
#                     '        join screen_cell_line using(cell_line_id) '
#                     '        where screen_id=screen.screen_id '
#                     '        order by c.value) as c1 )' % LIST_DELIMITER_SQL_ARRAY )
#                     .label('cell_lines'), 
                # TODO: convert to vocabulary
                'transfection_agent': literal_column(
                    '( select value from transfection_agent '
                    '  where screen.transfection_agent_id=transfection_agent_id) '
                    ),
                'date_of_first_activity': literal_column(
                    '( select date_of_activity '
                    '  from activity '
                    '  join lab_activity la using(activity_id) '
                    '  where la.screen_id=screen.screen_id '
                    '  order by date_of_activity asc LIMIT 1 )'
                    ),
                'date_of_last_activity': literal_column(
                    '( select date_of_activity '
                    '  from activity '
                    '  join lab_activity la using(activity_id) '
                    '  where la.screen_id=screen.screen_id '
                    '  order by date_of_activity desc LIMIT 1 )'
                    ),
                'screenresult_last_imported': literal_column(
                    '( select date_of_activity '
                    '  from activity '
                    '  join screen_result_update_activity srua on(update_activity_id=activity_id) '
                    '  join screen_result sr using(screen_result_id) '
                    '  where sr.screen_id=screen.screen_id '
                    '  order by date_of_activity desc LIMIT 1 )'
                    ),
                # TODO: convert to vocabulary
                'funding_supports': literal_column(
                    "(select array_to_string(array_agg(f1.value),'%s') "
                    '    from ( select fs.value from funding_support fs '
                    '        join screen_funding_support_link using(funding_support_id) '
                    '        where screen_id=screen.screen_id '
                    '        order by fs.value) as f1 )' % LIST_DELIMITER_SQL_ARRAY ), 
                # FIXME: use cherry_pick_liquid_transfer status vocabulary 
                'total_plated_lab_cherry_picks': 
                    select([tplcps.c.count]).\
                        select_from(tplcps).\
                        where(tplcps.c.screen_id==_screen.c.screen_id),
                # 'total_plated_lab_cherry_picks': literal_column(
                # '( select count(*) '                                             
                # '  from cherry_pick_request cpr '
                # '  join lab_cherry_pick lcp using(cherry_pick_request_id) '
                # '  join cherry_pick_assay_plate cpap using(cherry_pick_assay_plate_id) '
                # '  join cherry_pick_liquid_transfer cplt '
                # '    on(cherry_pick_liquid_transfer_id=activity_id)'
                # "  where cplt.status = 'Successful' "  
                # '  and cpr.screen_id=screen.screen_id )'),
                
                # TODO: convert to vocabulary
                'assay_readout_types': literal_column(
                    "(select array_to_string(array_agg(f1.assay_readout_type),'%s') "
                    '    from ( select distinct(assay_readout_type) '
                    '        from data_column ds '
                    '        join screen_result using(screen_result_id) '
                    '        where screen_id=screen.screen_id ) as f1 )' 
                    % LIST_DELIMITER_SQL_ARRAY ), 
                'library_plates_screened': 
                    select([lps.c.count]).\
                        select_from(lps).where(lps.c.screen_id==_screen.c.screen_id),
                # 'library_plates_screened': literal_column(
                #     '( select count(distinct(plate_number)) '
                #     '  from assay_plate where screen_id=screen.screen_id)'),
                'library_plates_data_loaded': 
                    select([lpdl.c.count]).\
                        select_from(lpdl).where(lpdl.c.screen_id==_screen.c.screen_id),
                # 'library_plates_data_loaded': literal_column(
                #     '( select count(distinct(plate_number)) '
                #     '  from assay_plate ' 
                #     '  where screen_result_data_loading_id is not null '
                #     '  and screen_id=screen.screen_id )'),
                'assay_plates_screened': 
                    select([aps.c.count]).\
                        select_from(aps).where(aps.c.screen_id==_screen.c.screen_id),
                # 'assay_plates_screened': literal_column(
                #     '( select count(*) '
                #     '  from assay_plate where screen_id=screen.screen_id )') ,
                'assay_plates_data_loaded': 
                    select([apdl.c.count]).\
                        select_from(apdl).where(apdl.c.screen_id==_screen.c.screen_id),
                # 'assay_plates_data_loaded': literal_column(
                #     '( select count(1) '
                #     '  from assay_plate '
                #     '  where screen_result_data_loading_id is not null '
                #     '  and screen_id=screen.screen_id )'),
                'libraries_screened_count': 
                    select([libraries_screened.c.count]).\
                        select_from(libraries_screened).\
                        where(libraries_screened.c.screen_id==_screen.c.screen_id ),
                # 'libraries_screened_count': literal_column(
                #     '( select count(distinct(l.library_id)) '
                #     '  from assay_plate ap '
                #     '  join plate p using(plate_id)' 
                #     '  join copy using(copy_id) '
                #     '  join library l using(library_id) '
                #     '  where screen_id=screen.screen_id )'),
                    
                # FIXME: use administrative activity vocabulary
                'last_data_loading_date': literal_column(
                    '( select activity.date_created '
                    '  from activity '
                    '  join administrative_activity aa using(activity_id) '
                    '  join screen_update_activity on update_activity_id=activity_id  '
                    "  where administrative_activity_type = 'Screen Result Data Loading' " 
                    '  and screen_id=screen.screen_id '
                    '  order by date_created desc limit 1 )'
                    ),
                'min_screened_replicate_count': 
                    select([func.min(apsrc.c.max_per_plate)+1]).\
                        select_from(apsrc).where(apsrc.c.screen_id==_screen.c.screen_id),
                'max_screened_replicate_count': 
                    select([func.max(apsrc.c.max_per_plate)+1]).\
                        select_from(apsrc).where(apsrc.c.screen_id==_screen.c.screen_id),
                'min_data_loaded_replicate_count': 
                    select([func.min(apdlrc.c.max_per_plate)+1]).\
                        select_from(apdlrc).where(apdlrc.c.screen_id==_screen.c.screen_id),
                'max_data_loaded_replicate_count': 
                    select([func.max(apdlrc.c.max_per_plate)+1]).\
                        select_from(apdlrc).where(apdlrc.c.screen_id==_screen.c.screen_id),
                
#                 'min_screened_replicate_count': literal_column(
#                     '( select max +1 from '
#                     '  ( select plate_number, max(replicate_ordinal) from '
#                     '    ( select plate_number, replicate_ordinal '
#                     '      from assay_plate where screen_id=screen.screen_id  '
#                     '      group by plate_number, replicate_ordinal '
#                     '      order by plate_number, replicate_ordinal asc nulls last) a '
#                     '  group by plate_number order by max asc nulls last) b limit 1 )'),
#                 'max_screened_replicate_count': literal_column(
#                     '( select max +1 from '
#                     '  ( select plate_number, max(replicate_ordinal) from '
#                     '    ( select plate_number, replicate_ordinal '
#                     '      from assay_plate where screen_id=screen.screen_id  '
#                     '      group by plate_number, replicate_ordinal '
#                     '      order by plate_number, replicate_ordinal desc nulls last) a '
#                     '  group by plate_number order by max desc nulls last) b limit 1 )'),
#                 'min_data_loaded_replicate_count': literal_column('\n'.join([
#                     '( select max +1 from ',
#                     '  ( select plate_number, max(replicate_ordinal) from ',
#                     '    ( select plate_number, replicate_ordinal ',
#                     '      from assay_plate ',
#                     '      where screen_id=screen.screen_id  ',
#                     '      and screen_result_data_loading_id is not null ',
#                     '      group by plate_number, replicate_ordinal ',
#                     '      order by plate_number, replicate_ordinal asc nulls last) a ',
#                     '  group by plate_number order by max asc nulls last) b limit 1 )'])),
#                 'max_data_loaded_replicate_count': literal_column('\n'.join([
#                     '( select max +1 from ',
#                     '  ( select plate_number, max(replicate_ordinal) from ',
#                     '    ( select plate_number, replicate_ordinal ',
#                     '      from assay_plate ',
#                     '      where screen_id=screen.screen_id  ',
#                     '      and screen_result_data_loading_id is not null ',
#                     '      group by plate_number, replicate_ordinal ',
#                     '      order by plate_number, replicate_ordinal desc nulls last) a ',
#                     '  group by plate_number order by max desc nulls last) b limit 1 )'])),
            }
            
            
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement

            j = _screen
            j = j.join(_screen_result, 
                _screen.c.screen_id==_screen_result.c.screen_id, isouter=True)
            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  

    def build_schema(self):
        schema = super(ScreenResource,self).build_schema()
        temp = [ x.screen_type for x in self.Meta.queryset.distinct('screen_type')]
        schema['extraSelectorOptions'] = { 
            'label': 'Type', 'searchColumn': 'screen_type', 'options': temp }
        return schema

    @transaction.atomic()
    def obj_create(self, bundle, **kwargs):
        bundle.data['date_created'] = timezone.now()
        bundle.data['version'] = 1
        bundle = super(ScreenResource, self).obj_create(bundle, **kwargs)
        logger.info(str(('created', bundle)))
        return bundle
    

class ScreenResultResource(SqlAlchemyResource,ManagedResource):

    username = fields.CharField('user__username', null=False, readonly=True)
    class Meta:
        queryset = ScreenResult.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'screenresult'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        allowed_methods = ['get']

        object_class = dict
        max_limit = 10000
        
    def __init__(self, **kwargs):
        self.scope = 'fields.screenresult'
        super(ScreenResultResource,self).__init__(**kwargs)
        
        conn = self.bridge.get_engine().connect()
        
        try:
            # create the well_query_index table if it does not exist
            conn.execute(text(
                'CREATE TABLE "well_query_index" ('
                '"well_id" text NOT NULL REFERENCES "well" ("well_id") DEFERRABLE INITIALLY DEFERRED,'
                '"query_id" integer NOT NULL REFERENCES "cached_query" ("id") DEFERRABLE INITIALLY DEFERRED'
                ');'
            ));
            logger.info(str(('the well_query_index table has been created')))
        except Exception, e:
            logger.warn(str(('on trying to create the well_query_index table', 
                str(e), e, 'note that this exception is normal if the table already exists',
                '(PostgreSQL <9.1 has no "CREATE TABLE IF NOT EXISTS" clause')))

        try:
            # create the well_data_column_positive_index table if it does not exist
            conn.execute(text(
                'CREATE TABLE well_data_column_positive_index ('
                ' "well_id" text NOT NULL REFERENCES "well" ("well_id") DEFERRABLE INITIALLY DEFERRED,'
                ' "data_column_id" integer NOT NULL ' 
                ' REFERENCES "data_column" ("data_column_id") DEFERRABLE INITIALLY DEFERRED'
                ');'
            ));
            logger.info(str(('the well_data_column_positive_index table has been created')))
        except Exception, e:
            logger.warn(str(('on trying to create the well_data_column_positive_index table', 
                str(e), e, 'note that this exception is normal if the table already exists',
                '(PostgreSQL <9.1 has no "CREATE TABLE IF NOT EXISTS" clause')))
            
    # FIXME: need test cases for this
    def clear_cache(self, all=False, by_date=None, by_uri=None, by_size=False):
        ManagedResource.clear_cache(self)

        max_indexes_to_cache = getattr(settings, 'MAX_WELL_INDEXES_TO_CACHE',2e+07)
        
        logger.info(str(('max_indexes_to_cache', max_indexes_to_cache)))
        conn = self.bridge.get_engine().connect()
        _wellQueryIndex = self.bridge['well_query_index']
        _well_data_column_positive_index = self.bridge['well_data_column_positive_index']
        
        try:
            query = CachedQuery.objects.filter(uri__contains='/screenresult/')
            logger.info(str(('queries to consider', [(x.id,x.uri) for x in query])))
            ids = set()
            if by_size:
                query = query.order_by('-datetime')
                cumulative_count = 0

                for q in query:
                    if q.count:
                        cumulative_count += q.count
                    if cumulative_count >= max_indexes_to_cache:
                        last_id_to_save = q.id
                        logger.info(str(('cumulative_count', cumulative_count, last_id_to_save)))
                        query = query.filter(id__lte=last_id_to_save)
                        ids.update([q.id for q in query])

                        logger.info(str(('ids', ids)))
                        break
                
            if by_date: # TODO: test
                logger.info(str(('by_date', by_date)))
                query = query.filter(datetime__lte=by_date)
                ids.update([q.id for q in query])
            if by_uri: # TODO: test
                logger.info(str(('by_uri', by_uri)))
                query = query.filter(uri__eq=by_uri)
                ids.update([q.id for q in query])
                
            logger.info(str(('clear CachedQuery by ids:', ids, 'all', all, 
                'by_date', by_date, 'by_uri', by_uri, 'by_size', by_size)))
            if ids or all:
                logger.info(str(('clear cachedQueries',ids, 'all', all)))
                if all:
                    stmt = delete(_wellQueryIndex)
                    conn.execute(stmt)
                    logger.info(str(('cleared all cached wellQueryIndexes')))
                    CachedQuery.objects.all().delete()
                    
                    stmt = delete(_well_data_column_positive_index)
                    conn.execute(stmt)
                    logger.info(str(('cleared all cached well_data_column_positive_indexes')))
                else:
                    
                    stmt = delete(_wellQueryIndex).where(_wellQueryIndex.c.query_id.in_(ids))
                    conn.execute(stmt)
                    logger.info(str(('cleared cached wellQueryIndexes', ids)))
                    CachedQuery.objects.filter(id__in=ids).delete()
                    
                    # delete the well_data_column_positive_indexes related to this uri
                    
            else:
                logger.info(str(('no CachedQueries or well_query_index values need to be cleared')))

        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on clear_cache', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  
            

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
        

    def get_detail(self, request, **kwargs):
        logger.info(str(('get_detail')))

        facility_id = kwargs.get('screen_facility_id', None)
        if not facility_id:
            logger.info(str(('no screen_facility_id provided')))
            raise NotImplementedError('must provide a screen_facility_id parameter')

        well_id = kwargs.get('well_id', None)
        if not well_id:
            logger.info(str(('no well_id provided')))
            raise NotImplementedError('must provide a well_id parameter')
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)
        
    
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

    def _build_result_value_column(self,field_information):
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
        column_to_select = None
        if(data_column_type in ['numeric','decimal','integer']): #TODO: use controlled vocabulary
            column_to_select = 'numeric_value'
        else:
            column_to_select = 'value'
        
        rv_select = select([column_to_select]).select_from(_rv)
        rv_select = rv_select.where(_rv.c.data_column_id == field_information['data_column_id'])
        rv_select = rv_select.where(_rv.c.well_id == text('assay_well.well_id'))
        rv_select = rv_select.limit(1) # FIXME: this is due to the duplicated result value issue
        rv_select = rv_select.label(field_name)
        return rv_select
    
    # store id's in a temp table version
    def build_list_response(self,request, param_hash={}, **kwargs):

        DEBUG_GET_LIST = True or logger.isEnabledFor(logging.DEBUG)
        
        if DEBUG_GET_LIST:
            logger.info(str(('kwargs',kwargs,'param_hash',param_hash)))

        is_for_detail = kwargs.pop('is_for_detail', False)

        screen_facility_id = param_hash.pop('screen_facility_id', None)
        if not screen_facility_id:
            logger.info(str(('no screen_facility_id provided')))
            raise NotImplementedError('must provide a screen_facility_id parameter')
            
        well_id = param_hash.pop('well_id', None)
        if well_id:
            param_hash['well_id__eq'] = well_id

        limit = param_hash.get('limit', 0)        
        try:
            limit = int(limit)
        except ValueError:
            raise BadRequest(
                "Invalid limit '%s' provided. Please provide a positive integer." % limit)
        param_hash['limit'] = 0
        
        offset = param_hash.get('offset', 0 )
        try:
            offset = int(offset)
        except ValueError:
            raise BadRequest(
                "Invalid offset '%s' provided. Please provide a positive integer." % offset)
        if offset < 0:    
            offset = -offset
        param_hash['offset'] = 0
        
        logger.info(str(('offset', offset, 'limit', limit )))
            
        show_mutual_positives = param_hash.get('show_mutual_positives', False)
        
        manual_field_includes = set(param_hash.get('includes', []))
                        
        try:
            screenresult = ScreenResult.objects.get(screen__facility_id=screen_facility_id)              
            logger.info(str(('screen_facility_id', screen_facility_id)))
            # general setup
             
            schema = self.build_schema(
                screenresult=screenresult,
                show_mutual_positives=show_mutual_positives)
            
            filename = self._get_filename(schema, kwargs)
        
            logger.info(str(('fields',schema['fields'].keys())))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
                                  
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = \
                    IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
            
            conn = self.bridge.get_engine().connect()
            _aw = self.bridge['assay_well']
            _sr = self.bridge['screen_result']
            _s = self.bridge['screen']
            _rv = self.bridge['result_value']
            
            # strategy: 
            # 1. create the base clause, which will build a stored index in 
            # the table: well_query_index;
            # the base query uses only columns: well_id, +filters, +orderings
            base_clause = _aw
            base_custom_columns = {
                'well_id': literal_column('assay_well.well_id'),
                }

            # The base all fields not part of result value lookups, and
            # any result_value lookups that are used in sort or where:
            # TODO: if doing mutual positives query, then need to grab those datacolumns.
            base_fields = [ fi for fi in field_hash.values() 
                if ( not fi.get('is_datacolumn',None) # don't include datacolumns, unless:
                     or fi['key'] in order_params
                     or '-%s'%fi['key'] in order_params
                     or fi['key'] in filter_fields )]
            
            for fi in [fi for fi in base_fields if fi.get('is_datacolumn',None)]:
                rv_select = self._build_result_value_column(fi)
                base_custom_columns[fi['key']] = rv_select
            
            base_query_tables = ['assay_well']
            
            logger.info(str(('=== build sqlalchemy columns', base_fields, base_custom_columns)))
            base_columns = self.build_sqlalchemy_columns(
                base_fields, base_query_tables=base_query_tables,
                custom_columns=base_custom_columns ) 
            logger.info(str(('base columns', base_columns)))
            # remove is_positive if not needed, to help the query planner
            if 'is_positive' not in filter_fields and 'is_positive' in base_columns:
                del base_columns['is_positive']
            
            base_stmt = select(base_columns.values()).select_from(base_clause)
            base_stmt = base_stmt.where(
                _aw.c.screen_result_id == screenresult.screen_result_id )
            # always add well_id order
            base_stmt = base_stmt.order_by(asc(column('well_id'))) # no need for nulls first
            
            (base_stmt,count_stmt) = \
                self.wrap_statement(base_stmt,order_clauses,filter_expression )

            meta = {
                'limit': limit,
                'offset': offset,
            }    

            # 1.a insert the base statement well ids into the indexing table
            
            m = hashlib.md5()
            compiled_stmt = str(base_stmt.compile(compile_kwargs={"literal_binds": True}))
            logger.info(str(('compiled_stmt', compiled_stmt)))
            m.update(compiled_stmt)
            key = m.hexdigest()
            
            cachedQuery = None
            _wellQueryIndex = self.bridge['well_query_index']
            index_was_created = False
            try:
                (cachedQuery,index_was_created) = CachedQuery.objects.all().get_or_create(key=key)
                if index_was_created:
                    cachedQuery.sql=compiled_stmt
                    cachedQuery.username = request.user.username
                    cachedQuery.uri = '/screenresult/%s' % screenresult.screen.facility_id
                    cachedQuery.params = json.dumps(param_hash)
                    cachedQuery.save()

                    base_stmt = base_stmt.alias('base_stmt')
                    insert_statement = insert(_wellQueryIndex).\
                        from_select(['well_id','query_id'],
                            select([
                                literal_column("well_id"),
                                literal_column(str(cachedQuery.id)).label('query_id')
                                ]).select_from(base_stmt))
                    
                    logger.info(str(('insert stmt', 
                        str(insert_statement.compile(compile_kwargs={"literal_binds": True})))))
                    conn.execute(insert_statement)
                    logger.info(str(('executed', compiled_stmt)))
                    
                    count_stmt = _wellQueryIndex
                    count_stmt = select([func.count()]).select_from(count_stmt)
                    count_stmt = count_stmt.where(_wellQueryIndex.c.query_id==cachedQuery.id)
                    
                    cachedQuery.count = int(conn.execute(count_stmt).scalar())
                    cachedQuery.save()
                else:
                    logger.info(str(('using cached well_query: ', cachedQuery)))
                
                meta['total_count'] = cachedQuery.count
                
            except Exception, e:
                logger.info(str(('ex on get/create wellquery for screenresult', e)))
                raise e
            
            #############            
            # now construct the "real" query - real columns, using the query_index table
            
            # specific setup 
            base_query_tables = ['assay_well','screen']
            
            # force query to use well_query_index.well_id
            custom_columns = { 'well_id': literal_column(
                'well_query_index.well_id').label('well_id')}
            
            for fi in [fi for fi in field_hash.values() if fi.get('is_datacolumn',None)]:
                custom_columns[fi['key']] = self._build_result_value_column(fi)
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            #### use the well_query_index well_ids as the central subquery loop
             
            _wqx = select(['well_id']).select_from(_wellQueryIndex)
            _wqx = _wqx.where(_wellQueryIndex.c.query_id==cachedQuery.id)
            if limit > 0:    
                _wqx = _wqx.limit(limit)
            _wqx = _wqx.offset(offset)

            logger.info(str(('base stmt', str(base_stmt.compile()))))
            _wqx = _wqx.alias('wqx')
            # force query to use well_query_index.well_id
            columns['well_id'] =  literal_column(
                'wqx.well_id').label('well_id')
            
            ######
            
            j = join(_wqx, _aw, _wqx.c.well_id==_aw.c.well_id)
            j = j.join(_sr, _aw.c.screen_result_id==_sr.c.screen_result_id)
            j = j.join(_s, _sr.c.screen_id==_s.c.screen_id)
            
            stmt = select(columns.values()).select_from(j)
            stmt = stmt.where(_aw.c.screen_result_id == screenresult.screen_result_id )
            logger.info(str(('excute stmt', str(stmt.compile(compile_kwargs={"literal_binds": True}) ))))
            result = conn.execute(stmt)
            logger.info('excuted stmt')
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            if index_was_created:
                self.clear_cache(by_size=True)
            
            if rowproxy_generator:
                result = rowproxy_generator(result)
            return self.stream_response_from_cursor(
                request, result, filename, 
                field_hash=field_hash, param_hash=param_hash,
                is_for_detail=is_for_detail,
                title_function=title_function,
                meta=meta  )
    
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  

    def get_mutual_positives_columns(self, screen_result_id):
        
        # TODO: cache this / clear cache when any screen_results referenced by 
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
        # and ( dc.data_type = 'boolean_positive_indicator' or dc.data_type = 'partition_positive_indicator' )
        # and aw0.screen_result_id = 941
        # and aw1.screen_result_id <> 941;
        try:
            conn = self.bridge.get_engine().connect()
            
            _well_data_column_positive_index = self.bridge['well_dc_positive_index']
            _aw = self.bridge['assay_well']
            _sr = self.bridge['screen_result']
            _dc = self.bridge['data_column']
    
            count_stmt = _well_data_column_positive_index
            count_stmt = select([func.count()]).select_from(count_stmt)
            count = int(conn.execute(count_stmt).scalar())        
            
            if count == 0:
                # cleared, recreate
                logger.info(str(('recreate well_data_column_positive_index')))
                base_stmt = join(_aw,_dc,_aw.c.screen_result_id==_dc.c.screen_result_id)
                base_stmt = select([
                    literal_column("well_id"),
                    literal_column('data_column_id')
                    ]).select_from(base_stmt)            
                base_stmt = base_stmt.where(_aw.c.is_positive)
                base_stmt = base_stmt.where(
                    _dc.c.data_type.in_(['boolean_positive_indicator','partition_positive_indicator']))
                base_stmt = base_stmt.order_by(_dc.c.data_column_id, _aw.c.well_id)
                insert_statement = insert(_well_data_column_positive_index).\
                    from_select(['well_id','data_column_id'], base_stmt)
                    
                logger.info(str(('mutual pos insert statement', str(insert_statement.compile()))))
                conn.execute(insert_statement)
    
            # now run the query
            # select distinct(wdc.data_column_id)
            # from
            # well_dc_positive_index wdc
            # join data_column dc using(data_column_id)
            # where 
            # dc.screen_result_id <> 941
            # and exists(
            # select null 
            # from well_dc_positive_index wdc1 
            # join data_column dc1 using(data_column_id) 
            # where wdc1.well_id=wdc.well_id 
            # and dc1.screen_result_id = 941 );        
            _wdc = _well_data_column_positive_index.alias('wdc')
            j = _wdc.join(_dc,_wdc.c.data_column_id==_dc.c.data_column_id)
            stmt = select([distinct(_wdc.c.data_column_id)]).select_from(j)
            stmt = stmt.where(_dc.c.screen_result_id != screen_result_id )
            
            _wdc1 = _well_data_column_positive_index.alias('wdc1')
            _dc1 = _dc.alias('dc1')
            j2 = _wdc1.join(_dc1, _wdc1.c.data_column_id==_dc1.c.data_column_id)
            stmt2 = select([text('null')]).select_from(j2)
            stmt2 = stmt2.where(_wdc.c.well_id==_wdc1.c.well_id)
            stmt2 = stmt2.where(_dc1.c.screen_result_id == screen_result_id )
            
            stmt = stmt.where(exists(stmt2))
            
            logger.info(str(('mutual positives statement', 
                str(stmt.compile(compile_kwargs={"literal_binds": True})))))
            
            return [x.data_column_id for x in conn.execute(stmt)]
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_mutual_positives_columns', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  

    def get_schema(self, request, **kwargs):
        if not 'screen_facility_id' in kwargs:
            raise Http404(unicode((
                'The screenresult schema requires a screen facility ID'
                ' in the URI, as in /screenresult/[facility_id]/schema/')))
        facility_id = kwargs.pop('screen_facility_id')
        try:
            logger.info(str(('find: ' , facility_id)))
            screenresult = ScreenResult.objects.get(
                screen__facility_id=facility_id)
            logger.info(str(('screenresult resource', 
                             facility_id,screenresult.screen)))
            
            if screenresult:
                return self.create_response(
                    request, self.build_schema(screenresult))
            else:
                raise Http404(unicode((
                    'no results for the screen: ', facility_id)))
        except ObjectDoesNotExist, e:
            raise Http404(unicode((
                'no screen found for facility id', facility_id)))
    
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
            
    def build_schema(self, screenresult=None,show_mutual_positives=False):
        logger.debug(str(('==========build schema for screen result', screenresult)))
        try:
            data = super(ScreenResultResource,self).build_schema()
            
            if screenresult:
                # find the highest ordinal in the non-data_column fields
                max_ordinal = 0
                for fi in data['fields'].values():
                    if fi.get('ordinal',0) > max_ordinal:
                        max_ordinal = fi['ordinal']
                # translate the datacolumn definitions into field information definitions
                field_defaults = {
                    'visibility': ['list','detail'],
                    'data_type': 'string',
                    'filtering': True,
                    'is_datacolumn': True,
                    'scope': 'datacolumn.screenresult-%s' % screenresult.screen.facility_id
                    }
                for i,dc in enumerate(DataColumn.objects.filter(
                    screen_result=screenresult).order_by('ordinal')):
                    
                    (columnName,_dict) = self.create_datacolumn(dc, 
                        field_defaults=field_defaults)
                    _dict['ordinal'] = max_ordinal + dc.ordinal + 1
                    data['fields'][columnName] = _dict
                    
                max_ordinal += dc.ordinal + 1
                field_defaults = {
                    'visibility': ['api'],
                    'data_type': 'string',
                    'filtering': True,
                    }
                
                _current_sr = None
                _ordinal = max_ordinal
                for i,dc in enumerate(DataColumn.objects.
                        filter(data_column_id__in=
                        self.get_mutual_positives_columns(screenresult.screen_result_id)).
                        order_by('screen_result__screen__facility_id','ordinal')):

                    if _current_sr != dc.screen_result.screen_result_id:
                        _current_sr = dc.screen_result.screen_result_id
                        max_ordinal = _ordinal
                    _ordinal = max_ordinal + dc.ordinal + 1
                    
                    field_defaults['scope'] = 'mutual_positive.%s' \
                        % dc.screen_result.screen.facility_id
                    (columnName,_dict) = self.create_datacolumn(dc, 
                        field_defaults=field_defaults)
                    _dict['ordinal'] = _ordinal
                    data['fields'][columnName] = _dict
                    if show_mutual_positives:
                        _visibility = set(_dict['visibility'])
                        _visibility.update(['list','detail'])
                        _dict['visibility'] = list(_visibility)
            return data
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on build_schema', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  
    
    def create_datacolumn(self,dc,field_defaults={}):
        screen_facility_id = dc.screen_result.screen.facility_id
        columnName = "dc_%s_%s" % (screen_facility_id, default_converter(dc.name))
        _dict = field_defaults.copy()
        _dict.update(model_to_dict(dc))
        _dict['title'] = '%s [%s]' % (dc.name, screen_facility_id) 
        _dict['description'] = _dict['description'] or _dict['title']
        _dict['comment'] = dc.comments
        _dict['key'] = columnName
        _dict['scope'] = 'datacolumn.screenresult-%s' % screen_facility_id
        if dc.data_type == 'numeric':
            if _dict.get('decimal_places', 0) > 0:
                _dict['data_type'] = 'decimal'
                _dict['display_options'] = \
                    '{ "decimals": %s }' % _dict['decimal_places']
            else:
                _dict['data_type'] = 'integer'
        if dc.data_type in self.data_type_lookup:
            _dict.update(self.data_type_lookup[dc.data_type])
        return (columnName,_dict)

class DataColumnResource(ManagedModelResource):
    '''
    A DataColumn is the metadata for a single column in a Screen result set.
    - extends the "Metahash Field" definition
    - this metadata is stored in the datacolumn table, for historic reasons 
    (DataColumns were defined for the ScreenResult design.)
    
    '''
    
    # included to allow querying like ?screen__facility_id=##
    screen = fields.ToOneField('db.api.ScreenResource', 'screen_result__screen')  
    screen_facility_id = fields.CharField('screen_result__screen__facility_id')
    
    class Meta:
        queryset = DataColumn.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'datacolumn'
        
        ordering = []
        filtering = { 'screen': ALL_WITH_RELATIONS}
        serializer = LimsSerializer()
        metahashResource = MetaHashResource()

    def __init__(self, **kwargs):
#        self.
        super(DataColumnResource,self).__init__(**kwargs)

    def build_schema(self):
        schema = ManagedModelResource.build_schema(self)
        
        meta_field_schema = self._meta.metahashResource.build_schema()
        
        # mix in metahash field definitions: datacolum extends the metahash field
        new_fields_schema = meta_field_schema['fields'].copy()
        for fi in new_fields_schema.values():
            # set all supertype visibility to detail - users will only want 
            # datacolumn definitions
            fi['visibility'] = ['detail'] 
        new_fields_schema.update(schema['fields'])
        schema['fields'] = new_fields_schema
        
        return schema
    
    def dehydrate(self, bundle):
        # FIXME: use the same method here and in the ScreenResultResource
        
        # Custom dehydrate: each datacolumn extends a Metahash Field
        # - augment each datacolumn with the Metahash field (supertype) fields
        # - these fields are critical to UI interpretation

        bundle =  ManagedModelResource.dehydrate(self, bundle)
        
        field_defaults = {
            'visibility': ['list','detail'],
            'data_type': 'string',
            'filtering': True,
            }
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
        
        dc = bundle.obj
        _dict = bundle.data
        _dict.update(field_defaults)
        columnName = "dc_%s" % (default_converter(dc.name))
        #         _dict = field_defaults.copy()
        #         _dict.update(model_to_dict(dc))
        
        _dict['title'] = dc.name
        _dict['comment'] = dc.comments
        _dict['key'] = columnName
        # so that the value columns come last
        _dict['ordinal'] += len(self.fields) + dc.ordinal 
        
        if dc.data_type == 'numeric':
            if _dict.get('decimal_places', 0) > 0:
                _dict['data_type'] = 'decimal'
                _dict['backgrid_cell_type'] = 'Iccbl.DecimalCell'
                _dict['display_options'] = \
                    '{ "decimals": %s }' % _dict['decimal_places']
            else:
                _dict['data_type'] = 'integer'
        if dc.data_type in data_type_lookup:
            _dict.update(data_type_lookup[dc.data_type])
            
        return bundle
       
    def prepend_urls(self):

        return [
            # override the parent "base_urls" so that schema is matched first
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            
            
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<data_column_id>\d+)%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            
            # TODO: implement a natural key for datacolumns
            #url((r"^(?P<resource_name>%s)/"
            #     r"(?P<screen_facility_id>\w+)/"
            #     r"(?P<key>[\w_]+)%s$") 
            #        % (self._meta.resource_name, trailing_slash()), 
            #    self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
    

class ScreenSummaryResource(ManagedModelResource):
        
    class Meta:
        queryset = Screen.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'screensummary'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
#        self.
        super(ScreenSummaryResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" allows us 
        # to match any word, except "schema", and use it as the key value to 
        # search for.
        # also note the double underscore "__" is because we also don't want to
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url(r"^(?P<resource_name>%s)/(?P<facility_id>((?=(schema))__|(?!(schema))[^/]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

    def dehydrate(self, bundle):
        screen = bundle.obj
        try:
            # TODO: this is an example of the old activity system; we'll want 
            # to refactor this to a generic entry in apilog and then the actual 
            # values
            activities = bundle.obj.screenupdateactivity_set.all().filter(
                update_activity__administrative_activity_type='Screen Result Data Loading')
            if len(activities) > 0: 
                bundle.data['screenresult_last_imported'] =  \
                    activities[:1][0].update_activity.activity.date_created;
        except ScreenResult.DoesNotExist, e:
            logger.info(unicode(('no screenresult for ', bundle.obj)))
        return bundle

# Deprecate - use apilog viewer
class CopyWellHistoryResource(SqlAlchemyResource, ManagedModelResource):
    class Meta:
        queryset = CopyWell.objects.all().order_by('well_id')
        
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'copywellhistory'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

        
    def __init__(self, **kwargs):
        super(CopyWellHistoryResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # Note: because this prepends the other list, we have to make sure 
        # "schema" is matched
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('search'), name="api_search"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<copy_name>[\w\d_.\-\+ ]+)" 
                r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<copy_name>[\w\d_.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]

    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with wells/
        logger.info(str(('get_detail')))

        copy_name = kwargs.get('copy_name', None)
        if not copy_name:
            logger.info(str(('no copy_name provided')))
            raise NotImplementedError('must provide a copy_name parameter')
        
        well_id = kwargs.get('well_id', None)
        if not well_id:
            logger.info(str(('no well_id provided')))
            raise NotImplementedError('must provide a well_id parameter')
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)
        
    
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

        
    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = super(CopyWellHistoryResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)
        
        well_id = param_hash.pop('well_id', None)
        if well_id:
            param_hash['well_id__eq'] = well_id

        copy_name = param_hash.pop('copy_name', None)
        if copy_name:
            param_hash['copy_name__eq'] = copy_name

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)

            if filter_expression is None:
                msgs = { 'Copy well resource': 'can only service requests with filter expressions' }
                logger.info(str((msgs)))
                raise ImmediateHttpResponse(response=self.error_response(request,msgs))
                                  
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
            base_query_tables = [
                'copy_well', 'copy', 'plate', 'well','library',
                'well_volume_adjustment','activity']
            
            # NOTE: date_time is included here as an exercise:
            # why db table structure needs to be redone
            custom_columns = {
                'consumed_volume': literal_column(
                    'initial_volume-copy_well.volume').label('consumed_volume'),
                'date_time': literal_column('\n'.join([
                    'case when wva.well_volume_correction_activity_id is not null then (', 
                    'select a1.date_created from activity a1', 
                    'where a1.activity_id = wva.well_volume_correction_activity_id )',  
                    'else ( select a2.date_created from activity a2', 
                    'join cherry_pick_assay_plate cpap on(cpap.cherry_pick_liquid_transfer_id=a2.activity_id)',
                    'join lab_cherry_pick lcp on(lcp.cherry_pick_assay_plate_id=cpap.cherry_pick_assay_plate_id)',
                    'where lcp.lab_cherry_pick_id = wva.lab_cherry_pick_id ) ',
                    'end',
                    ])).label('date_time'),
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement

            _cw = self.bridge['copy_well']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _p = self.bridge['plate']
            _w = self.bridge['well']
            _wva = self.bridge['well_volume_adjustment']
            _a = self.bridge['activity']
            
            _wva = _wva.alias('wva')
            j = join(_cw, _c, _c.c.copy_id == _cw.c.copy_id )
            j = j.join(_p, _cw.c.plate_id == _p.c.plate_id )
            j = j.join(_w, _cw.c.well_id == _w.c.well_id )
            j = j.join(_l, _w.c.library_id == _l.c.library_id )
            j = j.join(_wva,onclause=(and_(
                _cw.c.copy_id == _wva.c.copy_id,_cw.c.well_id == _wva.c.well_id)),
                isouter=True)
#             j = j.join(_wva,_cw.c.well_id == _wva.c.well_id)
            j = j.join(_a, _wva.c.well_volume_correction_activity_id == _a.c.activity_id, isouter=True )
            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  
   

class CopyWellResource(SqlAlchemyResource, ManagedModelResource):
    
    class Meta:
        queryset = CopyWell.objects.all().order_by('well_id')
        
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'copywell'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

    def __init__(self, **kwargs):
        super(CopyWellResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # Note: because this prepends the other list, we have to make sure 
        # "schema" is matched
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('search'), name="api_search"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<copy_name>[\w\d_.\-\+ ]+)" 
                r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),

            url(r"^(?P<resource_name>%s)/(?P<library_short_name>[\w\d_.\-\+: ]+)"
                r"/(?P<copy_name>[\w\d_.\-\+ ]+)" 
                r"/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),

            url(r"^(?P<resource_name>%s)/(?P<library_short_name>[\w\d_.\-\+: ]+)"
                r"/(?P<copy_name>[\w\d_.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]

    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with wells/
        logger.info(str(('get_detail', kwargs)))

        library_short_name = kwargs.get('library_short_name', None)
        if not library_short_name:
            logger.info(str(('no library_short_name provided')))
         
        copy_name = kwargs.get('copy_name', None)
        if not copy_name:
            logger.info(str(('no copy_name provided')))
            raise NotImplementedError('must provide a copy_name parameter')
        
        well_id = kwargs.get('well_id', None)
        if not well_id:
            logger.info(str(('no well_id provided')))
            raise NotImplementedError('must provide a well_id parameter')
        
        kwargs['is_for_detail']=True
        
        return self.get_list(request, **kwargs)
        
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = super(CopyWellResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)
        
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
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)

            if filter_expression is None:
                msgs = { 'Copy well resource': 'can only service requests with filter expressions' }
                logger.info(str((msgs)))
                raise ImmediateHttpResponse(response=self.error_response(request,msgs))
                                  
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
            base_query_tables = ['copy_well', 'copy', 'plate', 'well','library']
            
            custom_columns = {
                'consumed_volume': literal_column('initial_volume-volume').label('consumed_volume'),
                # query plan makes this faster than the hash join of copy-copy_well
#                 'copy_name': literal_column(
#                     '( select copy.name from copy where copy.copy_id=copy_well.copy_id )'
#                     ).label('copy_name')
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement

            _cw = self.bridge['copy_well']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _p = self.bridge['plate']
            _w = self.bridge['well']
            
#             j = join(_cw, _c, _c.c.copy_id == _cw.c.copy_id )
            j = join(_cw, _w, _cw.c.well_id == _w.c.well_id )
            j = j.join(_p, _cw.c.plate_id == _p.c.plate_id )
            j = j.join(_c, _cw.c.copy_id == _c.c.copy_id )
            j = j.join(_l, _w.c.library_id == _l.c.library_id )
            
            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            if not order_clauses:
                stmt = stmt.order_by('copy_name','plate_number', 'well_id')
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  
  

class CherryPickRequestResource(SqlAlchemyResource,ManagedModelResource):        
    class Meta:
        queryset = CherryPickRequest.objects.all().order_by('well_id')
        
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'cherrypickrequest'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

    def __init__(self, **kwargs):
        super(CherryPickRequestResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # Note: because this prepends the other list, we have to make sure 
        # "schema" is matched
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with wells/
        logger.info(str(('get_detail', kwargs)))

        cherry_pick_request_id = kwargs.get('cherry_pick_request_id', None)
        if not cherry_pick_request_id:
            logger.info(str(('no cherry_pick_request_id provided')))
            raise NotImplementedError('must provide a cherry_pick_request_id parameter')
        
        kwargs['is_for_detail']=True
        
        return self.get_list(request, **kwargs)
        
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = super(CherryPickRequestResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)
        
        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id:
            param_hash['cherry_pick_request_id__eq'] = cherry_pick_request_id

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)
                                  
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
            base_query_tables = ['cherry_pick_request']
            
            custom_columns = {
                'screen_facility_id': literal_column(
                    '( select facility_id '
                    '  from screen where screen.screen_id=cherry_pick_request.screen_id )'
                    ).label('screen_facility_id'),
                'requested_by_name': literal_column(
                    '( select su.first_name || $$ $$ || su.last_name'
                    '  from screensaver_user su '
                    '  where su.screensaver_user_id=cherry_pick_request.requested_by_id )'
                    ).label('requested_by_name'),
                'lab_head_name': literal_column(
                    '( select su.first_name || $$ $$ || su.last_name'
                    '  from screensaver_user su '
                    '  join screen s on(lab_head_id=su.screensaver_user_id) '
                    '  where s.screen_id=cherry_pick_request.screen_id )'
                    ).label('lab_head_name'),
                'lab_head_id': literal_column(
                    '( select s.lab_head_id'
                    '  from screen s  '
                    '  where s.screen_id=cherry_pick_request.screen_id )'
                    ).label('lab_head_id'),
                'lead_screener_name': literal_column(
                    '( select su.first_name || $$ $$ || su.last_name'
                    '  from screensaver_user su '
                    '  join screen s on(lead_screener_id=su.screensaver_user_id) '
                    '  where s.screen_id=cherry_pick_request.screen_id )'
                    ).label('lead_screener_name'),
                'lead_screener_id': literal_column(
                    '( select s.lead_screener_id'
                    '  from screen s '
                    '  where s.screen_id=cherry_pick_request.screen_id )'
                    ).label('lead_screener_id'),
                'screen_type': literal_column(
                    '( select s.screen_type'
                    '  from screen s  '
                    '  where s.screen_id=cherry_pick_request.screen_id )'
                    ).label('screen_type'),
                # REDO as cherry pick status - not correct yet
                # query is wrong because iscompleted means has cplt's but no cplt's that have no status
                'is_completed': literal_column('\n'.join([
                    '( select count(*) > 0 ',
                    '  from cherry_pick_assay_plate cpap ',
                    '  join cherry_pick_liquid_transfer cplt ',
                    '    on( cpap.cherry_pick_liquid_transfer_id=cplt.activity_id ) ',
                    '  where cpap.cherry_pick_request_id=cherry_pick_request.cherry_pick_request_id ',
                    '  AND cplt.status is not null )'])).label('is_completed'), 
                # 
                'number_plates': literal_column('\n'.join([
                    '( select count(distinct(plate_ordinal)) ',
                    '  from cherry_pick_assay_plate cpap ',
                    '  where cpap.cherry_pick_request_id=cherry_pick_request.cherry_pick_request_id )'])
                    ).label('number_plates'), 
                'number_plates_completed': literal_column('\n'.join([
                    '( select count(*) ',
                    '  from cherry_pick_assay_plate cpap ',
                    '  join cherry_pick_liquid_transfer cplt ',
                    '    on( cpap.cherry_pick_liquid_transfer_id=cplt.activity_id ) ',
                    '  where cpap.cherry_pick_request_id=cherry_pick_request.cherry_pick_request_id ',
                    '  AND cplt.status in ($$Successful$$,$$Canceled$$) )'])
                    ).label('number_plates_completed'), 
                'total_number_lcps': literal_column(
                    'lcp.count ' 
                    ).label('total_number_lcps'),
                # following not performant
#                 'total_number_lcps': literal_column(
#                     '( select count(*) from lab_cherry_pick lcp ' 
#                     '  where lcp.cherry_pick_request_id=cherry_pick_request.cherry_pick_request_id )'
#                     ).label('total_number_lcps'),
                'plating_activity_date': literal_column('\n'.join([
                    '( select date_of_activity ',
                    '  from activity ',
                    '  join cherry_pick_liquid_transfer cplt using(activity_id) ',
                    ('  join cherry_pick_assay_plate cpap on(cpap.cherry_pick_liquid_transfer_id='
                        'cplt.activity_id) '),
                    '  where cpap.cherry_pick_request_id=cherry_pick_request.cherry_pick_request_id ',
                    '  order by date_of_activity desc LIMIT 1 ) '])).label('plating_activity_date')
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement
            _cpr = self.bridge['cherry_pick_request']
            _count_lcp_stmt = text(
                '( select cherry_pick_request_id, count(*) '
                ' from lab_cherry_pick '
                ' group by cherry_pick_request_id '
                ' order by cherry_pick_request_id ) as lcp ' )
#             _count_lcp_stmt = _count_lcp_stmt.alias('lcp') 
            j = join(_cpr,_count_lcp_stmt, 
                _cpr.c.cherry_pick_request_id == literal_column('lcp.cherry_pick_request_id'), 
                isouter=True)
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            if not order_clauses:
                stmt = stmt.order_by(nullslast(desc(column('cherry_pick_request_id'))))
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  


class CherryPickPlateResource(SqlAlchemyResource,ManagedModelResource):        

    class Meta:
        queryset = CherryPickAssayPlate.objects.all().order_by('cherry_pick_request_id')
        
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'cherrypickassayplate'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

    def __init__(self, **kwargs):
        super(CherryPickPlateResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # Note: because this prepends the other list, we have to make sure 
        # "schema" is matched
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/(?P<plate_ordinal>[\d]+)"
                r"/(?P<attempt_ordinal>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with wells/
        logger.info(str(('get_detail', kwargs)))

        cherry_pick_request_id = kwargs.get('cherry_pick_request_id', None)
        if not cherry_pick_request_id:
            logger.info(str(('no cherry_pick_request_id provided')))
            raise NotImplementedError('must provide a cherry_pick_request_id parameter')

        plate_ordinal = kwargs.get('plate_ordinal', None)
        if not plate_ordinal:
            logger.info(str(('no plate_ordinal provided')))
            raise NotImplementedError('must provide a plate_ordinal parameter')

        attempt_ordinal = kwargs.get('attempt_ordinal', None)
        if not attempt_ordinal:
            logger.info(str(('no attempt_ordinal provided')))
            raise NotImplementedError('must provide a attempt_ordinal parameter')

        kwargs['is_for_detail']=True
        
        return self.get_list(request, **kwargs)
        
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = super(CherryPickPlateResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)
        
        cherry_pick_request_id = param_hash.pop('cherry_pick_request_id', None)
        if cherry_pick_request_id:
            param_hash['cherry_pick_request_id__eq'] = cherry_pick_request_id

        plate_ordinal = param_hash.pop('plate_ordinal', None)
        if plate_ordinal:
            param_hash['plate_ordinal__eq'] = plate_ordinal

        attempt_ordinal = param_hash.pop('attempt_ordinal', None)
        if attempt_ordinal:
            param_hash['attempt_ordinal__eq'] = attempt_ordinal

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)
                                  
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
            base_query_tables = [
                'cherry_pick_assay_plate',
                'cherry_pick_request',
                'cherry_pick_liquid_transfer',
                'activity']
        
            custom_columns = {
                'screen_facility_id': literal_column(
                    '( select facility_id '
                    '  from screen where screen.screen_id=cherry_pick_request.screen_id )'
                    ).label('screen_facilty_id'),
                'plated_by_name': literal_column(
                    '( select su.first_name || $$ $$ || su.last_name'
                    '  from screensaver_user su '
                    '  where su.screensaver_user_id=activity.performed_by_id )'
                    ).label('plated_by_name'),
                'plated_by_id': literal_column(
                    'activity.performed_by_id').label('plated_by_id'),
                 # plate name will be constructed further from other parts
                'plate_name': literal_column(
                    ('cherry_pick_assay_plate.cherry_pick_request_id '
                        '|| $$:$$ || plate_ordinal || $$:$$ || attempt_ordinal ')
                    ).label('plate_name'),
                'number_plates': literal_column('\n'.join([
                    '( select count(distinct(plate_ordinal)) ',
                    '  from cherry_pick_assay_plate cpap ',
                    '  where cpap.cherry_pick_request_id=cherry_pick_assay_plate.cherry_pick_request_id )'])
                    ).label('number_plates'), 
                # TODO: is_plated, is_screened replacing status label for now
                'is_screened': literal_column(
                    '(exists ( '
                    '     select null from cherry_pick_assay_plate_screening_link '
                    '     where cherry_pick_assay_plate_id'
                        '=cherry_pick_assay_plate.cherry_pick_assay_plate_id )) '
                        ).label('is_screened'),
                'last_screening_date': literal_column('\n'.join([
                    '( select date_of_activity ',
                    '  from activity ',
                    '  join cherry_pick_assay_plate_screening_link cpapsl ',
                    '     on(cherry_pick_screening_id=activity_id) ',
                    ('  where cpapsl.cherry_pick_assay_plate_id'
                        '=cherry_pick_assay_plate.cherry_pick_assay_plate_id '),
                    '  order by date_of_activity desc LIMIT 1 ) '])
                    ).label('plating_activity_date'),
                'last_screened_by_id': literal_column('\n'.join([
                    '( select performed_by_id ',
                    '  from activity ',
                    '  join cherry_pick_assay_plate_screening_link cpapsl ',
                    '     on(cherry_pick_screening_id=activity_id) ',
                    ('  where cpapsl.cherry_pick_assay_plate_id'
                        '=cherry_pick_assay_plate.cherry_pick_assay_plate_id '),
                    '  order by date_of_activity desc LIMIT 1 ) '])
                    ).label('last_screened_by_id'),
                'last_screened_by_name': literal_column('\n'.join([
                    '( select su.first_name || $$ $$ || su.last_name'
                    '  from screensaver_user su ',
                    '  join activity a on(a.performed_by_id=su.screensaver_user_id) ',
                    '  join cherry_pick_assay_plate_screening_link cpapsl ',
                    '     on(cherry_pick_screening_id=activity_id) ',
                    ('  where cpapsl.cherry_pick_assay_plate_id'
                        '=cherry_pick_assay_plate.cherry_pick_assay_plate_id '),
                    '  order by date_of_activity desc LIMIT 1 ) '])
                    ).label('last_screened_by_name'),
                'screening_activities': literal_column('\n'.join([
                    '(select array_to_string(array_agg(activity_id),$$,$$) ',
                    '  from (select activity_id from ', 
                    '    activity ',
                    '    join cherry_pick_assay_plate_screening_link cpapsl ',
                    '      on(cherry_pick_screening_id=activity_id) ',
                    ('  where cpapsl.cherry_pick_assay_plate_id'
                        '=cherry_pick_assay_plate.cherry_pick_assay_plate_id '),
                    '  order by date_of_activity desc ) a ) '])
                    ).label('screening_activities')
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement
            _cpap = self.bridge['cherry_pick_assay_plate']
            _cpr = self.bridge['cherry_pick_request']
            _cplt = self.bridge['cherry_pick_liquid_transfer']
            _cplta = self.bridge['activity']
            j = join(_cpap,_cpr,
                _cpap.c.cherry_pick_request_id==_cpr.c.cherry_pick_request_id)
            j = j.join(_cplt,
                _cpap.c.cherry_pick_liquid_transfer_id==_cplt.c.activity_id,
                isouter=True)
            j = j.join(_cplta,
                _cplt.c.activity_id==_cplta.c.activity_id)
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            if not order_clauses:
                stmt = stmt.order_by(nullslast(desc(column('cherry_pick_request_id'))))
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  


class LibraryCopyResource(SqlAlchemyResource, ManagedModelResource):

    library_short_name = fields.CharField('library__short_name',  null=True)
    
    class Meta:
        queryset = Copy.objects.all().order_by('name')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'librarycopy'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

    def __init__(self, **kwargs):
        super(LibraryCopyResource,self).__init__(**kwargs)

    def prepend_urls(self):
        return [
            # override the parent "base_urls" so that we don't need to worry 
            # about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url((r"^(?P<resource_name>%s)"
                 r"/(?P<library_short_name>[\w\d_.\-\+: ]+)"
                 r"/(?P<name>[^/]+)%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/search/(?P<search_ID>[\d]+)%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('search'), name="api_search"),

            url((r"^(?P<resource_name>%s)"
                 r"/(?P<library_short_name>[\w\d_.\-\+: ]+)"
                 r"/(?P<name>[^/]+)"
                 r"/plate%s$" ) 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_librarycopyplateview'), 
                name="api_dispatch_librarycopy_plateview"),
        ]    

    def dispatch_librarycopyplateview(self, request, **kwargs):
        logger.info(str(('dispatch_librarycopyplateview', kwargs)))
#         kwargs['library_short_name'] = kwargs.pop('library__short_name')  
        kwargs['copy_name'] = kwargs.pop('name')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)    
        
    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with wells/
        logger.info(str(('get_detail')))

        library_short_name = kwargs.get('library_short_name', None)
        if not library_short_name:
            logger.info(str(('no library_short_name provided')))
            raise NotImplementedError('must provide a library_short_name parameter')
        
        copy_name = kwargs.get('name', None)
        if not copy_name:
            logger.info(str(('no copy "name" provided')))
            raise NotImplementedError('must provide a copy "name" parameter')
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)
    
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)
        
    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)
             
        schema = super(LibraryCopyResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)
        
        library_short_name = param_hash.pop('library_short_name',
            param_hash.get('library_short_name__eq',None))
        
        if not library_short_name:
            logger.info(str(('no library_short_name provided')))
        else:
            param_hash['library_short_name__eq'] = library_short_name
        
        name = param_hash.pop('name', param_hash.get('name',None))
        if name:
            param_hash['name__eq'] = name
            
        logger.info(str(('get_list', filename, kwargs)))
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)

            if filter_expression is None:
                msgs = { 'Library copies resource': 
                    'can only service requests with filter expressions' }
                logger.info(str((msgs)))
                raise ImmediateHttpResponse(response=self.error_response(request,msgs))
                                  
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = \
                    IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 

            custom_columns = {
                'copy_id': literal_column('c1.copy_id'),
                'library_short_name': literal_column(
                    'c1.short_name').label('library_short_name'),
                'plate_screening_count': literal_column('c1.plate_screening_count'),
                'copy_plate_count': literal_column('c1.copy_plate_count'),
                'plate_screening_count_average': literal_column(
                    'c1.plate_screening_count::float/c1.copy_plate_count::float ').\
                    label('plate_screening_count_average'),
                'avg_plate_volume': literal_column('c2.avg_plate_volume'),
                'min_plate_volume': literal_column('c2.min_plate_volume'), 
                'max_plate_volume': literal_column('c2.max_plate_volume'), 
                'screening_count': literal_column('c3.screening_count'),
                'ap_count': literal_column('c3.ap_count'),
                'dl_count': literal_column('c3.dl_count'),
                'first_date_data_loaded': literal_column('c3.first_date_data_loaded'), 
                'last_date_data_loaded': literal_column('c3.last_date_data_loaded'), 
                'first_date_screened': literal_column('c3.first_date_screened'), 
                'last_date_screened': literal_column('c3.last_date_screened'),
                'primary_plate_location': literal_column('\n'.join([
                    "( select room || '-' || freezer || '-' || shelf || '-' || bin ", 
                    '    from plate_location pl ' ,
                    '    where pl.plate_location_id=copy.primary_plate_location_id) '])).\
                    label('primary_plate_location'),
                'plate_locations': literal_column('\n'.join([
                    '(select count(distinct(plate_location_id)) ',
                    '    from plate p',
                    '    where p.copy_id = copy.copy_id ) '])).\
                    label('plate_locations'),
                'plates_available': literal_column('\n'.join([
                    '(select count(p)', 
                    '    from plate p ',
                    '    where p.copy_id=copy.copy_id', 
                    "    and p.status = 'Available' ) "])).\
                    label('plates_available'),
                }
            
            base_query_tables = ['copy','library']
 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement

            _c = self.bridge['copy']
            _l = self.bridge['library']
            _ap = self.bridge['assay_plate']
            
            # = copy volume statistics = 

            text_cols = ','.join([
                'c.copy_id',
                'c.name',
                'l.short_name',
                'avg(p.avg_remaining_volume) avg_plate_volume', 
                'min(p.min_remaining_volume) min_plate_volume', 
                'max(p.max_remaining_volume) max_plate_volume'])
            text_select = '\n'.join([
                'plate p ' 
                'join copy c using(copy_id) '
                'join library l using(library_id) '
                ])
            copy_volume_statistics = select([text(text_cols)]).\
                select_from( text(text_select))
                
            if library_short_name:
                copy_volume_statistics = \
                    copy_volume_statistics.where(literal_column('l.short_name') == library_short_name )
            copy_volume_statistics = \
                copy_volume_statistics.group_by(text('c.copy_id, c.name, l.short_name '))
            copy_volume_statistics = copy_volume_statistics.order_by(text('c.name '))
            copy_volume_statistics = copy_volume_statistics.cte('copy_volume_statistics')
    
            # = copy plate screening statistics = 
            
            text_cols = ','.join([
                'c.copy_id',
                'c.name',
                'l.short_name',
                ('(select count(distinct(p)) from plate p where p.copy_id=c.copy_id) '
                    'as copy_plate_count'),
                'count(distinct(ls1)) as plate_screening_count' ])
            text_select = '\n'.join([
                'copy c', 
                'join plate p using(copy_id) ',
                'join library l using(library_id) ',
                'left join ( select ap.plate_id, ls.activity_id ',
                '    from assay_plate ap', 
                '    join library_screening ls on(ap.library_screening_id=ls.activity_id) ',
                '    where ap.replicate_ordinal=0 ) as ls1 on(ls1.plate_id=p.plate_id) ',
                ])
            copy_plate_screening_statistics = select([text(text_cols)]).\
                select_from( text(text_select))
            if library_short_name:
                copy_plate_screening_statistics = \
                    copy_plate_screening_statistics.where(
                        literal_column('l.short_name') == library_short_name )
            
            copy_plate_screening_statistics = \
                copy_plate_screening_statistics.group_by(
                    'c.copy_id, c.name, l.short_name')
            copy_plate_screening_statistics = \
                copy_plate_screening_statistics.cte('copy_plate_screening_statistics')
    
            # = copy screening statistics = 
            
            # NOTE: precalculated version
            copy_screening_statistics = select([text('*')])\
                    .select_from(text('copy_screening_statistics'))\

            copy_screening_statistics = copy_screening_statistics.where(
                literal_column('short_name') == library_short_name)

            copy_screening_statistics = copy_screening_statistics.cte('copy_screening_statistics')
            
            c1 = copy_plate_screening_statistics.alias('c1')
            c2 = copy_volume_statistics.alias('c2')
            c3 = copy_screening_statistics.alias('c3')
            
            # TODO: join only if columns are included!!
            
            j = join(_c, _l, _c.c.library_id == _l.c.library_id)
            j = j.outerjoin(c1, _c.c.copy_id == text('c1.copy_id'))
            j = j.outerjoin(c2,text('c1.copy_id = c2.copy_id') )
            j = j.outerjoin(c3, text('c1.copy_id = c3.copy_id') )
            
            logger.info(str(('====j', str(j))))
            
            stmt = select(columns.values()).select_from(j)
            
            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            if not order_clauses:
                stmt = stmt.order_by("library_short_name","name")
 
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']

            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, param_hash=param_hash, is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
            
                        
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e   

    def obj_create(self, bundle, **kwargs):
        bundle.data['date_created'] = timezone.now()
        
        bundle.data['version'] = 1
        logger.info(str(('===creating library copy', bundle.data)))

        return super(LibraryCopyResource, self).obj_create(bundle, **kwargs)

    def is_valid(self, bundle):
        """
        Should set a dictionary of error messages (in the bundle). 
        If the dictionary has
        zero items, the data is considered valid. If there are errors, keys
        in the dictionary should be field names and the values should be a list
        of errors, even if there is only one.
        """
        
        fields = MetaHash.objects.get_and_parse(
            scope='fields.librarycopy', field_definition_scope='fields.metahash')
        
        # cribbed from tastypie.validation.py - mesh data and obj values, then validate
        data = {}
        if bundle.obj.pk:
            data = model_to_dict(bundle.obj)
        if data is None:
            data = {}
        data.update(bundle.data)
        
        # do validations
        errors = defaultdict(list)
        
        usage_type = data.get('usage_type')
        if usage_type:
            field_def = fields['usage_type']
            if usage_type not in field_def['choices']:
                errors['usage_type'] = str(('value is not one of the choices', 
                    usage_type, field_def['choices']))
            
        
        if errors:
            bundle.errors[self._meta.resource_name] = errors
            logger.warn(str((
                'bundle errors', bundle.errors, len(bundle.errors.keys()))))
            return False
        return True

class UserChecklistItemResource(ManagedSqlAlchemyResourceMixin, ManagedModelResource):    

    class Meta:
        queryset = UserChecklistItem.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['digested_password']
        resource_name = 'userchecklistitem'
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        self.user_resource = None
        super(UserChecklistItemResource,self).__init__(**kwargs)

    def prepend_urls(self):
        return [
            # override the parent "base_urls" so that we don't need to worry 
            # about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<checklist_item_event_id>([\d_]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]    

    def get_detail(self, request, **kwargs):
        logger.info(str(('get_detail')))

#         screensaver_user_id = kwargs.get('screensaver_user_id', None)
#         username = kwargs.get('username', None)
#         if not (screensaver_user_id or username):
#             logger.info(str(('no screensaver_user_id or username provided',kwargs)))
#             raise NotImplementedError('must provide a screensaver_user_id or username parameter')
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)
       
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = self.build_schema()

        filename = self._get_filename(schema, kwargs)
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)
                  
            visibilities = set()                
            if is_for_detail:
                visibilities.update(['detail','summary'])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail, visibilities=visibilities)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.\
                    create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup
            _su = self.bridge['screensaver_user']
            _admin = _su.alias('admin')
            _sru = self.bridge['screening_room_user']
            _up = self.bridge['reports_userprofile']
            _uci = self.bridge['user_checklist_item']
            
            _vocab = self.bridge['reports_vocabularies']
            
            # get the checklist items & groups
            cig_table = ( 
                select([
                    _vocab.c.key.label('item_group'),
                    func.concat('checklistitemgroup','.',
                        _vocab.c.key,'.','name').label('checklistitemgroup')])
                .select_from(_vocab)
                .where(_vocab.c.scope=='checklistitem.group') )
            cig_table = Alias(cig_table)
            ci_table = (
                select([
                    cig_table.c.item_group,
                    _vocab.c.key.label('item_name')])
                .select_from(
                    _vocab.join(cig_table,
                        _vocab.c.scope==cig_table.c.checklistitemgroup))
                )
            ci_table = ci_table.cte('ci')

            # build the entered checklists
            
            j = _uci
            j = j.join(_su, _uci.c.screensaver_user_id==_su.c.screensaver_user_id)
            j = j.join(_admin, _uci.c.admin_user_id==_admin.c.screensaver_user_id)
            entered_checklists = select([
                _su.c.username,
                func.concat( _su.c.last_name,', ', _su.c.first_name ).label('user_fullname'),
                _uci.c.item_group,
                _uci.c.item_name,
                _uci.c.status,
                _uci.c.status_date,
                _admin.c.username.label('admin_username')
                ]).select_from(j)
            username = param_hash.pop('username', None)
            if username:
                entered_checklists = entered_checklists.where(
                    _su.c.username==username)
            entered_checklists = entered_checklists.cte('entered_checklists')
            
            # This entire query doesn't fit the pattern, so have to construct it manually
            # bleah
            custom_columns = {
                'username': func.coalesce(entered_checklists.c.username,username),
                'user_fullname': entered_checklists.c.user_fullname,
                'admin_username': entered_checklists.c.admin_username,
                'item_group': ci_table.c.item_group,
                'item_name' : ci_table.c.item_name,
                'status': func.coalesce(entered_checklists.c.status,'not_completed'),
                'status_date': entered_checklists.c.status_date
                }

            base_query_tables = ['user_checklist_item','screensaver_user'] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            isouter=False
            if username:
                # if username, then this is a user specific view:
                # - outer join in the two so that a full list is generated
                isouter=True
                
            j = ci_table
            j = j.join(entered_checklists,
                ci_table.c.item_name==entered_checklists.c.item_name,isouter=isouter)
            
            stmt = select(columns.values()).select_from(j)
            if username:
                # if username, then this is a user specific view:
                stmt = stmt.order_by(ci_table.c.item_group,ci_table.c.item_name)
            else:
                stmt = stmt.order_by(entered_checklists.c.username)
            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  

     
class ScreensaverUserResource(ManagedSqlAlchemyResourceMixin, ManagedModelResource):    

    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['digested_password']
        resource_name = 'screensaveruser'
        max_limit = 10000
        always_return_data = True

    def __init__(self, **kwargs):
        self.user_resource = None
        super(ScreensaverUserResource,self).__init__(**kwargs)

    def prepend_urls(self):
        return [
            # override the parent "base_urls" so that we don't need to worry 
            # about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))/groups%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_user_groupview'), name="api_dispatch_user_groupview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))/checklistitems%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_user_checklistitemview'), name="api_dispatch_user_checklistitemview"),
        ]    

    def dispatch_user_groupview(self, request, **kwargs):
        return UserGroupResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_checklistitemview(self, request, **kwargs):
        return UserChecklistItemResource().dispatch('list', request, **kwargs)    
    
    def build_schema(self):
        
        schema = super(ScreensaverUserResource,self).build_schema()
        
        sub_schema = self.get_user_resource().build_schema();
        fields = {}
        fields.update(sub_schema['fields'])
        for key,val in schema['fields'].items():
            if key in fields:
                fields[key].update(val)
            else:
                fields[key] = val
        
        schema['fields'] = fields;
        
        return schema

    def get_user_resource(self):
        if not self.user_resource:
            self.user_resource = UserResource()
        return self.user_resource
        
    def get_detail(self, request, **kwargs):
        logger.info(str(('get_detail')))

        screensaver_user_id = kwargs.get('screensaver_user_id', None)
        username = kwargs.get('username', None)
        if not (screensaver_user_id or username):
            logger.info(str(('no screensaver_user_id or username provided',kwargs)))
            raise NotImplementedError('must provide a screensaver_user_id or username parameter')
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)
       
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = self.build_schema()

        filename = self._get_filename(schema, kwargs)
        
        screensaver_user_id = param_hash.pop('screensaver_user_id', None)
        if screensaver_user_id:
            param_hash['screensaver_user_id__eq'] = screensaver_user_id

        username = param_hash.pop('username', None)
        if username:
            param_hash['username__eq'] = username

        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)
                  
            visibilities = set()                
            if is_for_detail:
                visibilities.update(['detail','summary'])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail, visibilities=visibilities)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup
            _s = self.bridge['screen']
            _cl = self.bridge['collaborator_link']
            _su = self.bridge['screensaver_user']
            _admin = self.bridge['administrator_user']
            _sru = self.bridge['screening_room_user']
            _lh = self.bridge['lab_head']
            _la = self.bridge['lab_affiliation']
            _sru_fr = self.bridge['screening_room_user_facility_usage_role']
            _su_r = self.bridge['screensaver_user_role']
            
            _lhsu = _su.alias('lhsu')
            
            custom_columns = {
                'name': literal_column(
                    "screensaver_user.last_name || ', ' || screensaver_user.first_name"),
                # TODO: redo classification
                'classification': literal_column(
                    "( select coalesce( user_classification, 'admin') "
                    ' from  screensaver_user s1'
                    ' left join screening_room_user using(screensaver_user_id) '
                    ' where s1.screensaver_user_id = screensaver_user.screensaver_user_id)'),
                'screens_lab_head':
                    select([func.array_to_string(
                            func.array_agg(_s.c.facility_id),LIST_DELIMITER_SQL_ARRAY)]).\
                        select_from(_s).\
                        where(_s.c.lab_head_id==_su.c.screensaver_user_id),
                'screens_lead':
                    select([func.array_to_string(
                            func.array_agg(_s.c.facility_id),LIST_DELIMITER_SQL_ARRAY)]).\
                        select_from(_s).\
                        where(_s.c.lead_screener_id==_su.c.screensaver_user_id),
                'screens_collaborator':
                    select([func.array_to_string(
                            func.array_agg(_s.c.facility_id),LIST_DELIMITER_SQL_ARRAY)]).\
                        select_from(_s.join(_cl,_s.c.screen_id==_cl.c.screen_id)).\
                        where(_cl.c.collaborator_id==_su.c.screensaver_user_id),
                'lab_name':
                    select([func.concat(_lhsu.c.last_name,', ',
                        _lhsu.c.first_name,' - ',_la.c.affiliation_name)]).\
                    select_from(
                        _la.join(_lh,_la.c.lab_affiliation_id==_lh.c.lab_affiliation_id).\
                        join(_lhsu, _lh.c.screensaver_user_id==_lhsu.c.screensaver_user_id).\
                        join(_sru,_lh.c.screensaver_user_id==_sru.c.lab_head_id )).\
                    where(_sru.c.screensaver_user_id==_su.c.screensaver_user_id),
                'facility_usage_roles': 
                    select([func.array_to_string(
                            func.array_agg(_sru_fr.c.facility_usage_role),
                                LIST_DELIMITER_SQL_ARRAY)]).\
                        select_from(_sru_fr).\
                        where(_sru_fr.c.screening_room_user_id==_su.c.screensaver_user_id),
                'data_access_roles': 
                    select([func.array_to_string(
                            func.array_agg(_su_r.c.screensaver_user_role),
                                LIST_DELIMITER_SQL_ARRAY)]).\
                        select_from(_su_r).\
                        where(_su_r.c.screensaver_user_id==_su.c.screensaver_user_id),
                }

            # delegate to the user resource
            default_fields = ['fields.screensaveruser','fields.user']
            _temp = { key:field for key,field in field_hash.items() 
                if field.get('scope', None) in default_fields }
            field_hash = _temp
            logger.info(str(('final field hash: ', field_hash.keys())))
            logger.info(str(('key','last_name',field_hash['last_name'])))
            sub_columns = self.get_user_resource().build_sqlalchemy_columns(
                field_hash.values(),
                custom_columns=custom_columns)

            base_query_tables = ['screensaver_user','reports_user','auth_user'] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=sub_columns )

            # build the query statement
            _au = self.bridge['auth_user']
            _up = self.bridge['reports_userprofile']
            
            j = _su
            j = _su.join(_up,_su.c.user_id==_up.c.id).\
                    join(_au,_up.c.user_id==_au.c.id)
            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  


    # reworked 20150706   
    @un_cache        
    @transaction.atomic()
    def put_list(self,request, **kwargs):

        deserialized = self._meta.serializer.deserialize(
            request.body, 
            format=request.META.get('CONTENT_TYPE', 'application/json'))
        if not self._meta.collection_name in deserialized:
            raise BadRequest("Invalid data sent, must be nested in '%s'" 
                % self._meta.collection_name)
        deserialized = deserialized[self._meta.collection_name]

        logger.info(str(('put list', deserialized,kwargs)))
        
        with transaction.atomic():
            
            # TODO: review REST actions:
            # PUT deletes the endpoint
            
            ScreensaverUser.objects.all().delete()
            
            for _dict in deserialized:
                self.put_obj(_dict)

        # get new state, for logging
        response = self.get_list(
            request,
            desired_format='application/json',
            includes='*',
            **kwargs)
        new_data = self._meta.serializer.deserialize(
            LimsSerializer.get_content(response), format='application/json')
        new_data = new_data[self._meta.collection_name]
        
        logger.info(str(('new data', new_data)))
        self.log_patches(request, [],new_data,**kwargs)

    
        if not self._meta.always_return_data:
            return http.HttpAccepted()
        else:
            response = self.get_list(request, **kwargs)             
            response.status_code = 202
            return response
        
    @un_cache        
    @transaction.atomic()
    def patch_list(self, request, **kwargs):

        deserialized = self._meta.serializer.deserialize(
            request.body, 
            format=request.META.get('CONTENT_TYPE', 'application/json'))
        if not self._meta.collection_name in deserialized:
            raise BadRequest("Invalid data sent, must be nested in '%s'" 
                % self._meta.collection_name)
        deserialized = deserialized[self._meta.collection_name]

        logger.info(str(('patch list', deserialized,kwargs)))

        # cache state, for logging
        response = self.get_list(
            request,
            desired_format='application/json',
            includes='*',
            **kwargs)
        original_data = self._meta.serializer.deserialize(
            LimsSerializer.get_content(response), format='application/json')
        original_data = original_data[self._meta.collection_name]
        logger.info(str(('original data', original_data)))

        with transaction.atomic():
            
            for _dict in deserialized:
                self.patch_obj(_dict)
                
        # get new state, for logging
        response = self.get_list(
            request,
            desired_format='application/json',
            includes='*',
            **kwargs)
        new_data = self._meta.serializer.deserialize(
            LimsSerializer.get_content(response), format='application/json')
        new_data = new_data[self._meta.collection_name]
        
        logger.info(str(('new data', new_data)))
        self.log_patches(request, original_data,new_data,**kwargs)
        
        if not self._meta.always_return_data:
            return http.HttpAccepted()
        else:
            response = self.get_list(request, **kwargs)             
            response.status_code = 202
            return response
                
    @un_cache        
    @transaction.atomic()
    def patch_detail(self, request, **kwargs):

        deserialized = self._meta.serializer.deserialize(
            request.body, 
            format=request.META.get('CONTENT_TYPE', 'application/json'))
        logger.info(str(('patch detail', deserialized,kwargs)))

        # cache state, for logging
        username = self.get_user_resource().find_username(deserialized, **kwargs)
        response = self.get_list(
            request,
            desired_format='application/json',
            includes='*',
            **kwargs)
        original_data = self._meta.serializer.deserialize(
            LimsSerializer.get_content(response), format='application/json')
        original_data = original_data[self._meta.collection_name]
        logger.info(str(('original data', original_data)))

        with transaction.atomic():
            
            self.patch_obj(deserialized, **kwargs)

        # get new state, for logging
        response = self.get_list(
            request,
            desired_format='application/json',
            includes='*',
            **kwargs)
        new_data = self._meta.serializer.deserialize(
            LimsSerializer.get_content(response), format='application/json')
        new_data = new_data[self._meta.collection_name]
        
        logger.info(str(('new data', new_data)))
        self.log_patches(request, original_data,new_data,**kwargs)

        
        if not self._meta.always_return_data:
            return http.HttpAccepted()
        else:
            response = self.get_detail(request, **kwargs) 
            response.status_code = 202
            return response

    @un_cache        
    @transaction.atomic()
    def put_detail(self, request, **kwargs):
                
        deserialized = self._meta.serializer.deserialize(
            request.body, 
            format=request.META.get('CONTENT_TYPE', 'application/json'))

        logger.info(str(('put detail', deserialized,kwargs)))
        
        with transaction.atomic():
            logger.info(str(('put_detail:', kwargs)))
            
            self.put_obj(deserialized, **kwargs)
        
        # FIXME: logging
        
        if not self._meta.always_return_data:
            return http.HttpAccepted()
        else:
            response = self.get_detail(request, **kwargs) 
            response.status_code = 202
            return response
        
    @transaction.atomic()    
    def put_obj(self,deserialized, **kwargs):
        self.clear_cache()
        
        try:
            self.delete_obj(deserialized, **kwargs)
        except ObjectDoesNotExist,e:
            pass 
        
        return self.patch_obj(deserialized, **kwargs)
    
    @un_cache        
    @transaction.atomic()
    def delete_detail(self,deserialized, **kwargs):
        deserialized = self._meta.serializer.deserialize(
            request.body, 
            format=request.META.get('CONTENT_TYPE', 'application/json'))
        try:
            self.delete_obj(deserialized, **kwargs)
            # FIXME: logging
            return HttpResponse(status=202)
        except ObjectDoesNotExist,e:
            return HttpResponse(status=404)
    
    @transaction.atomic()    
    def delete_obj(self, deserialized, **kwargs):
        self.clear_cache()
        username = self.get_user_resource().find_username(deserialized,**kwargs)
        ScreensaverUser.objects.get(username=username).delete()
    
    @transaction.atomic()    
    def patch_obj(self,deserialized, **kwargs):
        self.clear_cache()

        username = self.get_user_resource().find_username(deserialized,**kwargs)
        
        schema = self.build_schema()
        fields = schema['fields']

        fields = { name:val for name,val in fields.items() 
            if val['table'] and val['table']=='screensaver_user'}
        
        try:
            # create/get userprofile

            user = self.get_user_resource().patch_obj(deserialized,**kwargs)

            # create the screensaver_user
            
            initializer_dict = {}
            for key in fields.keys():
                if deserialized.get(key,None):
                    initializer_dict[key] = parse_val(
                        deserialized.get(key,None), key,fields[key]['data_type']) 

            screensaver_user = None
            try:
                # FIXME: add username field to the screensaver_user table
                screensaver_user = ScreensaverUser.objects.get(user__username=username)
            except ObjectDoesNotExist, e:
                logger.info('Screensaver User %s does not exist, creating' % username)
                # FIXME: add username field to the screensaver_user table
                screensaver_user = ScreensaverUser.objects.create(username=username)
                screensaver_user.save()
            
            logger.info(str(('initializer dict', initializer_dict)))
            for key,val in initializer_dict.items():
                logger.info(str(('set',key,val,hasattr(screensaver_user, key))))
                
                if hasattr(screensaver_user,key):
                    setattr(screensaver_user,key,val)
            
            # also set
            screensaver_user.username = user.username
            screensaver_user.email = user.email
            screensaver_user.user = user
            screensaver_user.save()
            return screensaver_user
            
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on put detail', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  

        
class ReagentResource(SqlAlchemyResource, ManagedModelResource):
    
    class Meta:

        queryset = Reagent.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
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
        super(ReagentResource,self).__init__(**kwargs)
    
    def prepend_urls(self):
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
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
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)
        
    def build_list_response(self,request, param_hash={}, **kwargs):
        
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        logger.info(str(('kwargs', kwargs)))

        # TODO: eliminate dependency on library (for schema determination)
        library = None
        
        library_short_name = param_hash.pop('library_short_name', None)
        if not library_short_name:
            logger.info(str(('no library_short_name provided')))
        else:
            param_hash['library_short_name__eq'] = library_short_name
            library = Library.objects.get(short_name=library_short_name)

        well_id = param_hash.pop('well_id', None)
        if well_id:
            param_hash['well_id__eq'] = well_id
            if not library:
                library = Well.objects.get(well_id=well_id).library

        #         plate_number = param_hash.pop('plate_number', 
        #             param_hash.get('plate_number__eq', None))
        #         if plate_number:
        #             param_hash['plate_number__eq'] = int(plate_number)

        substance_id = param_hash.pop('substance_id', None)
        if substance_id:
            param_hash['substance_id__eq'] = substance_id
            if not library:
                library = Reagent.objects.get(substance_id=substance_id).well.library

        schema = self.build_schema(library=library)
          
        filename = self._get_filename(schema, kwargs)

        logger.info(str(('get_list', filename, param_hash)))

        try:
            
            # general setup
             
            manual_field_includes = set(param_hash.get('includes', []))
            desired_format = self.get_format(request)
            if desired_format == 'chemical/x-mdl-sdfile':
                manual_field_includes.add('molfile')
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash, **kwargs)
            
            if filter_expression is None:
                msgs = { 'reagent resource': 'can only service requests with filter expressions' }
                logger.info(str((msgs)))
                raise ImmediateHttpResponse(response=self.error_response(request,msgs))
                 
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes,
                is_for_detail=is_for_detail)
            
            logger.info(str(('field hash scopes', 
                set([field.get('scope', None) 
                    for field in field_hash.values()]) )) )
            if library:
                default_fields = ['fields.well','fields.reagent']
                if library.screen_type == 'rnai':
                    default_fields.append('fields.silencingreagent')
                elif library.screen_type == 'natural_products':
                    default_fields.append('fields.naturalproductreagent')
                else:
                    default_fields.append('fields.smallmoleculereagent')
                    
                _temp = { key:field for key,field in field_hash.items() 
                    if field.get('scope', None) in default_fields }
                field_hash = _temp
                logger.info(str(('final field hash: ', field_hash.keys())))
            else:
                # consider limiting fields available
                pass
            
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
        
            base_query_tables = ['well', 'reagent', 'library']
            
            columns = []
            sub_columns = self.get_sr_resource().build_sqlalchemy_columns(
                field_hash.values(), self.bridge)
            sub_columns['plate_number'] = ( literal_column(
                "to_char(well.plate_number,'%s')" % PLATE_NUMBER_SQL_FORMAT)
                .label('plate_number') )
            
            if DEBUG_GET_LIST: 
                logger.info(str(('sub_columns', sub_columns.keys())))
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=sub_columns)
            
# Could use the library param to limit the column building exercise
# to the sub-resource, but since all columns can be joined, just include
# the SR columns, as above.            
#             sub_resource = None
#             if library:
#                 sub_resource = self.get_reagent_resource(library_screen_type=library.screen_type)
#             if sub_resource and hasattr(sub_resource, 'build_sqlalchemy_columns'):
#                 sub_columns = sub_resource.build_sqlalchemy_columns(
#                     field_hash.values(), self.bridge)
#                 if DEBUG_GET_LIST: 
#                     logger.info(str(('sub_columns', sub_columns.keys())))
#                 columns = self.build_sqlalchemy_columns(
#                     field_hash.values(), base_query_tables=base_query_tables,
#                     custom_columns=sub_columns)
#             else:
# 
#                 sub_columns = sub_resource.build_sqlalchemy_columns(
#                     field_hash.values(), self.bridge)
#                 if DEBUG_GET_LIST: 
#                     logger.info(str(('sub_columns', sub_columns.keys())))
#                 columns = self.build_sqlalchemy_columns(
#                     field_hash.values(), base_query_tables=base_query_tables,
#                     custom_columns=sub_columns)
#                 
                
#                 # Note: excludes smr,rnai,np,... tables if library not specified
#                 logger.info(str(('build generic resource columns')))
#                 columns = self.build_sqlalchemy_columns(
#                     field_hash.values(), base_query_tables=base_query_tables)
            
            if DEBUG_GET_LIST: 
                logger.info(str(('columns', [str(col) for col in columns])))
            
            # Start building a query; use the sqlalchemy.sql.schema.Table API:
            _well = self.bridge['well']
            _reagent = self.bridge['reagent']
            _library = self.bridge['library']
            j = _well.join(_reagent, _well.c.well_id==_reagent.c.well_id, isouter=True)
            j = j.join(_library, _well.c.library_id == _library.c.library_id )
            stmt = select(columns.values()).select_from(j)
            
            if library:
                stmt = stmt.where(_well.c.library_id == library.library_id) 

            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            if not order_clauses:
                stmt = stmt.order_by("plate_number","well_name")

            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']

            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, param_hash=param_hash, is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
                        
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  
        
 
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
    
    def get_reagent_resource(self, library_screen_type=None):
        # FIXME: we should store the "type" on the entity
        
        if library_screen_type == 'rnai':
            return self.get_sr_resource()
        else:
            if library_screen_type == 'natural_products':
                return self.get_npr_resource()
            else:
                return self.get_smr_resource()
    
    def get_library_resource(self):
        if not self.library_resource:
            self.library_resource = LibraryResource()
        return self.library_resource

    def get_object_list(self, request, library_short_name=None):
        ''' 
        Note: any extra kwargs are there because we are injecting them into the 
        global tastypie kwargs in one of the various "dispatch_" handlers assigned 
        through prepend_urls.  Here we can explicitly add them to the query. 
        
        '''
        library = Library.objects.get(short_name=library_short_name)
        sub_resource = self.get_reagent_resource(library_screen_type=library.screen_type)
        query = sub_resource.get_object_list(request)
        logger.info(str(('==== query', query.query.sql_with_params())))
        
        ## also add in the "supertype" fields:
        query.select_related('well')
    
        if library_short_name:
            query = query.filter(well__library=library)
        return query

    def full_dehydrate(self, bundle, for_list=False):
        
        well_bundle = self.build_bundle(bundle.obj.well, request=bundle.request)
        well_bundle = self.get_well_resource().full_dehydrate(well_bundle)
        bundle.data.update(well_bundle.data)
        
        library = bundle.obj.well.library
        sub_resource = self.get_reagent_resource(library_screen_type=library.screen_type)
        bundle = sub_resource.full_dehydrate(bundle, for_list=for_list)
        
        return bundle
                
    def get_schema(self, request, **kwargs):
        if not 'library_short_name' in kwargs:
            return self.create_response(request, self.build_schema())
        
        library_short_name = kwargs.pop('library_short_name')
        try:
            library = Library.objects.get(short_name=library_short_name)
            return self.create_response(request, self.build_schema(library))
            
        except Library.DoesNotExist, e:
            raise Http404(unicode(( 'cannot build schema - library def needed'
                'no library found for short_name', library_short_name)))
                
    def build_schema(self, library=None):
        
        schema = super(ReagentResource,self).build_schema()

        # grab all of the subtypes
        
        sub_schema = self.get_npr_resource().build_schema();
        schema['fields'].update(sub_schema['fields']);
        
        sub_schema = self.get_sr_resource().build_schema();
        schema['fields'].update(sub_schema['fields']);
        
        sub_schema = self.get_smr_resource().build_schema();
        schema['fields'].update(sub_schema['fields']);
        
        well_schema = WellResource().build_schema()
        schema['fields'].update(well_schema['fields'])

        return schema

    def build_schema_old(self, library=None):
        
        schema = deepcopy(super(ReagentResource,self).build_schema())
        
        if library:
            sub_data = self.get_reagent_resource(library_screen_type=library.screen_type).build_schema()
            
            newfields = {}
            newfields.update(sub_data['fields'])
            newfields.update(schema['fields'])
            schema['fields'] = newfields
            
            for k,v in schema.items():
                if k != 'fields' and k in sub_data:
                    schema[k] = sub_data[k]
            
        well_schema = WellResource().build_schema()
        schema['fields'].update(well_schema['fields'])

        return schema


class WellResource(SqlAlchemyResource, ManagedModelResource):

    library_short_name = fields.CharField('library__short_name',  null=True)
    library = fields.CharField(null=True)
    
    class Meta:

        queryset = Well.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'well'
        ordering = []
        filtering = {}
        serializer = LimsSerializer()   

        xls_serializer = XLSSerializer()
        # Backbone/JQuery likes to JSON.parse the returned data
        always_return_data = True 
        max_limit = 10000

    def __init__(self, **kwargs):
        self.library_resource = None
        self.sr_resource = None
        self.smr_resource = None
        self.npr_resource = None
        self.reagent_resource = None
        super(WellResource,self).__init__(**kwargs)

    def deserialize(self, request, data=None, format=None):
        '''
        Override deserialize so we can pull apart the multipart form and get the 
        uploaded content.
        Note: native TP doesn't support multipart uploads, this will support
        standard multipart form uploads in modern browsers
        '''
        logger.info(str(('deserialize', format)))
        if not format:
            format = request.META.get('CONTENT_TYPE', 'application/json')

        if format.startswith('multipart'):
            if len(request.FILES.keys()) != 1:
                raise ImmediateHttpResponse(
                    response=self.error_response(request, 
                        { 'FILES', 'File upload supports only one file at a time'}))
            
            if 'sdf' in request.FILES:  
                # process *only* the first file
                file = request.FILES['sdf']
                format = 'chemical/x-mdl-sdfile'
                
                # NOTE: have to override super, because it ignores the format and 
                # grabs it again from the Request headers (which is "multipart...")
                #  return super(ReagentResource, self).deserialize(request, file, format) 
                deserialized = self._meta.serializer.deserialize(file.read(), format=format)

            elif 'xls' in request.FILES:
                # TP cannot handle binary file formats - it is calling 
                # django.utils.encoding.force_text on all input
                file = request.FILES['sdf']
                deserialized = self._meta.xls_serializer.from_xls(file.read())
            else:
                logger.error(str(('UnsupportedFormat', request.FILES.keys() )))
                raise UnsupportedFormat(str(('Unknown file type: ', request.FILES.keys()) ) )
        
        elif format == 'application/xls':
            # TP cannot handle binary file formats - it is calling 
            # django.utils.encoding.force_text on all input
            deserialized = self._meta.xls_serializer.from_xls(request.body)
            
        else:
            deserialized = super(WellResource, self).deserialize(request, request.body, format)    
        
        if self._meta.collection_name in deserialized: 
            # this is a list of data
            deserialized[self._meta.collection_name] = \
                self.create_aliasmapping_iterator(deserialized[self._meta.collection_name])
        else:   
            # this is a single item of data
            deserialized = self.alias_item(deserialized)
            
        return deserialized
    
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
    
    def get_reagent_resource(self, library_screen_type=None):
        # FIXME: we should store the "type" on the entity
        
        if library_screen_type == 'rnai':
            return self.get_sr_resource()
        else:
            if library_screen_type == 'natural_products':
                return self.get_npr_resource()
            else:
                return self.get_smr_resource()
    
    def get_full_reagent_resource(self):
        if not self.reagent_resource:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource

    def get_library_resource(self):
        if not self.library_resource:
            self.library_resource = LibraryResource()
        return self.library_resource

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
                
    def get_schema(self, request, **kwargs):
        if not 'library_short_name' in kwargs:
            return self.create_response(request, self.build_schema())
        
        library_short_name = kwargs.pop('library_short_name')
        try:
            library = Library.objects.get(short_name=library_short_name)
            return self.create_response(request, self.build_schema(library))
            
        except Library.DoesNotExist, e:
            raise Http404(unicode(( 'cannot build schema - library def needed'
                'no library found for short_name', library_short_name)))
                
    def build_schema(self, library=None):
        data = super(WellResource,self).build_schema()
        
        if library:
            sub_data = self.get_reagent_resource(library_screen_type=library.screen_type).build_schema()
            data = deepcopy(data)
            
            newfields = {}
            newfields.update(sub_data['fields'])
            newfields.update(data['fields'])
            data.update(sub_data)
            data['fields'] = newfields

        temp = [ x.title.lower() 
            for x in Vocabularies.objects.all().filter(scope='library.well_type')]
        data['extraSelectorOptions'] = { 
            'label': 'Type', 'searchColumn': 'library_well_type', 'options': temp }

        return data

    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with all resources
        logger.info(str(('get_detail')))
 
        well_id = kwargs.get('well_id', None)
        if not well_id:
            logger.info(str(('no well_id provided')))
            raise NotImplementedError('must provide a well_id parameter')
         
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)

    def get_list(self, request, **kwargs):
        return self.get_full_reagent_resource().get_list(request, **kwargs)

    def post_list(self, request, **kwargs):
        raise NotImplementedError("Post is not implemented for ReagentResource, use patch instead")
    
    def patch_list(self, request, **kwargs):
        # TODO: NOT TESTED
        return self.put_list(request, **kwargs)
    
    @transaction.atomic()
    def put_list(self, request, **kwargs):

        if 'library_short_name' not in kwargs:
            raise BadRequest('library_short_name is required')
        
#         deserialized = self.deserialize(
#             request, 
#             format=request.META.get('CONTENT_TYPE', 'application/json'))
        deserialized = self._meta.serializer.deserialize(
            request.body, 
            format=request.META.get('CONTENT_TYPE', 'application/json'))
        if not self._meta.collection_name in deserialized:
            raise BadRequest(str(("Invalid data sent. missing: " , self._meta.collection_name)))
        
        basic_bundle = self.build_bundle(request=request)
 
        library = Library.objects.get(short_name=kwargs['library_short_name'])
        prev_version = library.version_number
        if library.version_number:
            library.version_number += 1
        else:
            library.version_number = 1
        library.save()
        
        library_log = self.make_log(request)
        library_log.diff_keys = ['version_number']
        library_log.diffs = {
            'version_number': [prev_version, library.version_number]}
        library_log.ref_resource_name = 'library'
        library_log.uri = self.get_library_resource().get_resource_uri(library)
        library_log.key = '/'.join(
            [str(x) for x in self.get_library_resource().detail_uri_kwargs(library).values()])
        library_log.save()
        
                               
        # Cache all the wells on the library for use with this process 
        wellMap = dict( (well.well_id, well) for well in library.well_set.all())
        if len(wellMap)==0:
            raise BadRequest(str(('library has not been created, no wells', library)))
        
        i=0
        bundles_seen = []
        skip_errors=False
        for well_data in deserialized[self._meta.collection_name]:
            
            well_data['library_short_name']=kwargs['library_short_name']
            
            logger.debug(str(('well_data', well_data)))
            well_id = well_data.get('well_id', None)
            if not well_id:
                well_name = well_data.get('well_name', None)
                plate_number = well_data.get('plate_number',None)
                if well_name and plate_number:                
                    well_id = '%s:%s' %(str(plate_number).zfill(5), well_name)

            if not well_id:
                raise ImmediateHttpResponse(
                    response=self.error_response(request, {'well_id': 'required'}))
            
            well = wellMap.get(well_id, None)
            if not well:
                raise ImmediateHttpResponse(
                    response=self.error_response(request, {
                        'well_id': str(('well not found for this library', well_id))}))
                
            well_bundle = self.build_bundle(
                obj=well, data=well_data, request=request);
            
            kwargs.update({ 'library': library })
            kwargs.update({ 'parent_log': library_log })
            well_bundle = self.obj_update(well_bundle, **kwargs)
                
            i = i+1
            bundles_seen.append(well_bundle)
        
        logger.debug(str(('put reagents', i, library_log)))
        
        if not self._meta.always_return_data:
            return http.HttpNoContent()
        else:
            to_be_serialized = {}
            to_be_serialized[self._meta.collection_name] = [
                self.full_dehydrate(bundle, for_list=True) for bundle in bundles_seen]
            to_be_serialized = self.alter_list_data_to_serialize(request, to_be_serialized)
            return self.create_response(request, to_be_serialized)

    @log_obj_update      
    @transaction.atomic()
    def obj_update(self, well_bundle, **kwargs):
        # called only from local put_list  
            
        library = kwargs.pop('library')
        well_data = well_bundle.data
        
        well_bundle = self.full_hydrate(well_bundle)
        self.is_valid(well_bundle)
        if well_bundle.errors and not skip_errors:
            raise ImmediateHttpResponse(response=self.error_response(
                well_bundle.request, well_bundle.errors))
        well_bundle.obj.save()
        
        duplex_wells = []
        if well_data.get('duplex_wells', None):
            if not library.is_pool:
                raise ImmediateHttpResponse(
                    response=self.error_response(request, {
                        'duplex_wells': str(('library is not a pool libary', library))}))
            well_ids = well_data['duplex_wells'] #.split(';')
#             logger.info(str(('well_ids', well_ids, well_data['duplex_wells'])))
            for well_id in well_ids:
                try:
                    duplex_wells.append(Well.objects.get(well_id=well_id))
                except:
                    raise ImmediateHttpResponse(
                        response=self.error_response(well_bundle.request, {
                            'duplex_well not found': str(('pool well', well_bundle.obj, well_id))}))
                    
        logger.debug(str(('updated well', well_bundle.obj)))

        # lookup/create the reagent
        sub_resource = self.get_reagent_resource(library_screen_type=library.screen_type)
        
        reagent_bundle = sub_resource.build_bundle(
            data=well_data, request=well_bundle.request)
        if not well_bundle.obj.reagent_set.exists():
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(str(('==== creating reagent for ', well_bundle.obj)) )
            sub_resource.obj_create(
                reagent_bundle, 
                **{ 'well': well_bundle.obj, 'library': library, 'duplex_wells': duplex_wells })
        else:
            # NOTE: this only works if there is only one reagent in the well:
            # TODO: implement update for specific reagent through ReagentResource
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(str(('==== updating *first* reagent for ', well_bundle.obj)) )
            # lookup and update the reagent
            reagent_bundle.obj = well_bundle.obj.reagent_set.all()[0]
            sub_resource.obj_update(reagent_bundle)
            logger.debug(str(('updated reagent', reagent_bundle.obj)))
        
        return well_bundle


# TODO: Eventually, replace much of this with the ApiLog resource; 
# after determining best way to handle m2m reln's
class ActivityResource(SqlAlchemyResource,ManagedModelResource):

    performed_by = fields.ToOneField(
        'db.api.ScreensaverUserResource', 
        attribute='performed_by', 
        full=True, full_detail=True, full_list=True,
        null=True)
    performed_by_id = fields.IntegerField(attribute='performed_by_id');

    class Meta:
        queryset = AdministrativeActivity.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'activity'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

        
    def __init__(self, **kwargs):
        super(ActivityResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # Note: because this prepends the other list, we have to make sure 
        # "schema" is matched
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)"
                r"/(?P<cherry_pick_request_id>[\d]+)"
                r"/(?P<plate_ordinal>[\d]+)"
                r"/(?P<attempt_ordinal>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]

    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with wells/
        logger.info(str(('get_detail', kwargs)))

        kwargs['is_for_detail']=True
        
        return self.get_list(request, **kwargs)
            
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)

    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns django.http.response.StreamingHttpResponse 
        '''
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)

        schema = super(ActivityResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
  
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)
                                  
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
              
            order_params = param_hash.get('order_by',[])
            order_clauses = \
                SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = \
                    IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)
 
            # specific setup 
            base_query_tables = ['activity']
        
            custom_columns = {
                'screen_id': literal_column(
                    '( select facility_id '
                    '  from screen where screen.screen_id=cherry_pick_request.screen_id )'
                    ).label('screen_id'),
                'performed_by_name': literal_column(
                    "(select su.first_name || ' ' || su.last_name "
                    ' from activity a'
                    ' join screensaver_user su on(a.performed_by_id=su.screensaver_user_id) '
                    ' where a.activity_id=activity.id )').label('performed_by_name'),
                'activity_type': literal_column(
                    '()'
                    ).label('activity_type')
                    
            }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement
            _a = self.bridge['activity']
            _u = self.bridge['screensaver_user']
            
            
            j = join(_a,_u,
                _a.c.performed_by_id==_u.c.screensaver_user_id)
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            if not order_clauses:
                stmt = stmt.order_by(nullslast(desc(column('date_performed'))))
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  


class LibraryResource(SqlAlchemyResource, ManagedModelResource):
    
    class Meta:
        queryset = Library.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'library'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 
        
    def __init__(self, **kwargs):
        
        self.well_resource = None
        self.apilog_resource = None
        self.reagent_resource = None
        
        super(LibraryResource,self).__init__(**kwargs)
        
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

    def get_detail(self, request, **kwargs):
        # TODO: this is a strategy for refactoring get_detail to use get_list:
        # follow this with wells/
        logger.info(str(('get_detail', kwargs)))

        library_short_name = kwargs.pop('short_name', None)
        if not library_short_name:
            logger.info(str(('no library_short_name provided')))
            raise NotImplementedError('must provide a short_name parameter')
        else:
            kwargs['short_name__eq'] = library_short_name
        
        kwargs['is_for_detail']=True
        return self.get_list(request, **kwargs)
        
    
    def get_list(self,request,**kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)

        return self.build_list_response(request,param_hash=param_hash, **kwargs)
        
    def build_list_response(self,request, param_hash={}, **kwargs):
        ''' 
        Overrides tastypie.resource.Resource.get_list for an SqlAlchemy implementation
        @returns djanog.http.response.StreamingHttpResponse 
        '''
        
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        is_for_detail = kwargs.pop('is_for_detail', False)
        
        schema = super(LibraryResource,self).build_schema()

        filename = self._get_filename(schema, kwargs)

        try:
            # general setup
            
            manual_field_includes = set(param_hash.get('includes', []))
            if DEBUG_GET_LIST: 
                logger.info(str(('manual_field_includes', manual_field_includes)))
 
            (filter_expression, filter_fields) = \
                SqlAlchemyResource.build_sqlalchemy_filters(schema, param_hash=param_hash)
                
            field_hash = self.get_visible_fields(
                schema['fields'], filter_fields, manual_field_includes, 
                is_for_detail=is_for_detail)
             
            order_params = param_hash.get('order_by',[])
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(order_params, field_hash)
            
            rowproxy_generator = None
            if param_hash.get(HTTP_PARAM_USE_VOCAB,False):
                rowproxy_generator = IccblBaseResource.create_vocabulary_rowproxy_generator(field_hash)

            # specific setup
                                     
            custom_columns={
                'plate_count': literal_column(
                    '(select count(distinct(p.plate_id))'
                    '    from plate p join copy c using(copy_id)'
                    '    where c.library_id=library.library_id)').label('plate_count'), 
                'copies': literal_column(
                    "(select array_to_string(array_agg(c1.name),'%s') "
                    '    from ( select c.name from copy c '
                    '    where c.library_id=library.library_id '
                    '    order by c.name) as c1 )' % LIST_DELIMITER_SQL_ARRAY ).label('copies'), 
                # TODO: copies2 is the same in all respects, except that it is 
                # used differently in the UI
                'copies2': literal_column(
                    "(select array_to_string(array_agg(c1.name),'%s') "
                    '    from ( select c.name from copy c '
                    '    where c.library_id=library.library_id '
                    '    order by c.name) as c1 )' % LIST_DELIMITER_SQL_ARRAY ).label('copies2'), 
                'owner': literal_column(
                    "(select u.first_name || ' ' || u.last_name "
                    '    from screensaver_user u '
                    '    where u.screensaver_user_id=library.owner_screener_id)').label('owner')
                    };
                     
            base_query_tables = ['library']

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement
            _l = self.bridge['library']
            j = _l
            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(stmt,order_clauses,filter_expression )
            
            title_function = None
            if param_hash.get(HTTP_PARAM_USE_TITLES, False):
                title_function = lambda key: field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function  )
             
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  

    def prepend_urls(self):

        return [
            # override the parent "base_urls" so that we don't need to worry 
            # about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
                
#             url((r"^(?P<resource_name>%s)"
#                  r"/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))/schema%s$") 
#                     % (self._meta.resource_name, trailing_slash()), 
#                 self.wrap_view('get_schema'), name="api_get_schema"),

            # TODO: rework the "((?=(schema))__|(?!(schema))[^/]+)" to "[\w\d_.\-\+: ]+" used below
            # or even "[\/]+"
            url(r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/plate/(?P<plate_number>[^/]+)%s$") % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copyplateview'), 
                name="api_dispatch_library_copyplateview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/plate%s$") % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copyplateview'), 
                name="api_dispatch_library_copyplateview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/copywell/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$") % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copywellview'), 
                name="api_dispatch_library_copywellview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/copywellhistory/(?P<well_id>\d{1,5}\:[a-zA-Z]{1,2}\d{1,2})%s$") % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copywellhistoryview'), 
                name="api_dispatch_library_copywellhistoryview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/copywell%s$") % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copywellview'), 
                name="api_dispatch_library_copywellview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/copy/(?P<name>[^/]+)%s$") % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copyview'), 
                name="api_dispatch_library_copyview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/copy%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copyview'), 
#                 self.wrap_view('dispatch_librarycopyview'), 
                name="api_dispatch_library_copyview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/plate%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copyplateview'), 
                name="api_dispatch_library_copyplateview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/plate/(?P<plate_number>[^/]+)%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_copyplateview'), 
                name="api_dispatch_library_copyplateview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/well%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_wellview'), 
                name="api_dispatch_library_wellview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/reagent%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_reagentview'), 
                name="api_dispatch_library_reagentview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/reagent/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_reagent_schema'), name="api_get_reagent_schema"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/well/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_well_schema'), name="api_get_well_schema"),
                
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/version%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_libraryversionview'), 
                name="api_dispatch_libraryversionview"),
        ]    
    
    def get_well_schema(self, request, **kwargs):
        if not 'short_name' in kwargs:
            raise Http404(unicode((
                'The well schema requires a library short name'
                ' in the URI, as in /library/[short_name]/well/schema/')))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_well_resource().get_schema(request, **kwargs)    
  
    def get_reagent_schema(self, request, **kwargs):
        if not 'short_name' in kwargs:
            raise Http404(unicode((
                'The reagent schema requires a library short name'
                ' in the URI, as in /library/[short_name]/well/schema/')))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_reagent_resource().get_schema(request, **kwargs)    
  
    def dispatch_library_copyview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyResource().dispatch('list', request, **kwargs)    
 
    def dispatch_library_copyplateview(self, request, **kwargs):
        logger.info(str(('dispatch_library_copyplateview', kwargs)))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)   

    def dispatch_library_copywellview(self, request, **kwargs):
        logger.info(str(('dispatch_library_copywellview', kwargs)))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return CopyWellResource().dispatch('list', request, **kwargs)   

    def dispatch_library_copywellhistoryview(self, request, **kwargs):
        logger.info(str(('dispatch_library_copywellhistoryview', kwargs)))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return CopyWellHistoryResource().dispatch('list', request, **kwargs)   

    def dispatch_library_wellview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_well_resource().dispatch('list', request, **kwargs)    
                    
    def dispatch_library_reagentview(self, request, **kwargs):
        logger.info(str(('dispatch_library_reagentview ', kwargs)))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_reagent_resource().dispatch('list', request, **kwargs)    
                    
    def dispatch_libraryversionview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryContentsVersionResource().dispatch('list', request, **kwargs)    
        
    def build_schema(self, librarytype=None):
        schema = cache.get(self._meta.resource_name + ":schema")
        if not schema:
            # FIXME: these options should be defined automatically from a vocabulary in build_schema
            schema = super(LibraryResource,self).build_schema()
            
            if 'start_plate' in schema['fields'] and 'end_plate' in schema['fields']:
                # (only run if library fields already initialized)
                # Exemplary section - set start/end plate ranges
                maxsmr = ( 
                    Library.objects
                        .filter(screen_type='small_molecule')
                        .exclude(library_type='natural_products')
                        .aggregate(Max('end_plate')) )
                maxsmr = maxsmr['end_plate__max']      
                minrnai = ( 
                    Library.objects
                        .filter(screen_type='rnai')
                        .aggregate(Min('end_plate')) )
                minrnai = minrnai['end_plate__min']
                maxrnai = ( 
                    Library.objects
                        .filter(screen_type='rnai')
                        .aggregate(Max('end_plate')) )
                maxrnai = maxrnai['end_plate__max']
                schema['library_plate_range'] = [maxsmr,minrnai,maxrnai]
                schema['fields']['start_plate']['range'] = [maxsmr,minrnai,maxrnai]
                schema['fields']['end_plate']['range'] = [maxsmr,minrnai,maxrnai]
            
            temp = [ x.library_type for x in self.Meta.queryset.distinct('library_type')]
            schema['extraSelectorOptions'] = { 
                'label': 'Type', 'searchColumn': 'library_type', 'options': temp }
        return schema
    
    ##
    ## Note: @transaction.atomic() cannot be nested in commit_on_success, because
    ## of version compatability issues in django:
    ## "Starting with Django 1.6, atomic() is the only supported API for 
    ##  defining a transaction. Unlike the deprecated APIs, it'snestable and
    ##  always guarantees atomicity.
    ##
    @transaction.atomic()
    def obj_create(self, bundle, **kwargs):
        bundle.data['date_created'] = timezone.now()
        
        bundle.data['version'] = 1
        
        logger.debug(str(('===creating library', bundle.data)))

        bundle = super(LibraryResource, self).obj_create(bundle, **kwargs)

        # clear the cached schema because plate range have updated
        cache.delete(self._meta.resource_name + ':schema')
        
        # now create the wells
        
        library = bundle.obj
        logger.debug(str((
            'created library', library, library.start_plate, type(library.start_plate))))
        plate_size = int(library.plate_size)

        try:
            i =0
            for plate in range(int(library.start_plate), int(library.end_plate)+1):
                for index in range(0,plate_size):
                    well = Well()
                    # FIXME: get rid of version
                    well.version = 1
                    well.well_name = lims_utils.well_name_from_index(index,plate_size)
                    well.well_id = lims_utils.well_id(plate,well.well_name)
                    well.library = library
                    well.plate_number = plate
                    # FIXME: use vocabularies for well type
                    well.library_well_type = 'undefined'
                    well.save()
                    i += 1
            logger.info(str(('created', i, 'wells for library', library.short_name, library.library_id )))
            return bundle
        except Exception, e:

            extype, ex, tb = sys.exc_info()
            msg = str(e)
            if isinstance(e, ImmediateHttpResponse):
                msg = str(e.response)
            logger.warn(str((
                'throw', e, msg, tb.tb_frame.f_code.co_filename, 'error line', 
                tb.tb_lineno, extype, ex)))
            
            
            msg = str(e)
            if isinstance(e, ImmediateHttpResponse):
                msg = str(e.response)
            errMsg = str(('on creating wells for Library', library.short_name, msg))
            logger.warn(errMsg)
            raise ImmediateHttpResponse(response=self.error_response(
                bundle.request, { 'errMsg': errMsg }))


class LibraryContentsVersionResource(ManagedModelResource):

    library_short_name = fields.CharField('library__short_name',  null=True)
    loading_activity = fields.ToOneField(
        'db.api.ActivityResource', 
        attribute='library_contents_loading_activity__activity', 
        full=True, full_detail=True, full_list=True,
        null=True)
    release_activity = fields.ToOneField(
        'db.api.ActivityResource', 
        attribute='library_contents_release_activity__activity', 
        full=True, full_detail=True, full_list=True,
        null=True)
     
    date_loaded = fields.DateField(
        'library_contents_loading_activity__activity__date_of_activity', null=True)
    date_released = fields.DateField(
        'library_contents_release_activity__activity__date_of_activity', null=True)
    load_commments = fields.CharField(
        'library_contents_loading_activity__activity__comments', null=True)
    loaded_by_id = fields.IntegerField(
        'library_contents_loading_activity__activity__performed_by__screensaver_user_id',
        null=True)
        
    class Meta:
        queryset = LibraryContentsVersion.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'librarycontentsversion'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(LibraryContentsVersionResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<library__short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/(?P<version_number>[^/]+)%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]    
    
    def get_object_list(self, request, library_short_name=None):
        ''' 
        Note: any extra kwargs are there because we are injecting them into the 
        global tastypie kwargs in one of the various "dispatch_" handlers assigned 
        through prepend_urls.  Here we can explicitly add them to the query. 
        
        '''
        query = super(LibraryContentsVersionResource, self).get_object_list(request);
        if library_short_name:
            query = query.filter(library__short_name=library_short_name)
        return query

    def dehydrate(self, bundle):
        if bundle.obj.library_contents_loading_activity:
            sru = bundle.obj.library_contents_loading_activity.activity.performed_by
            bundle.data['loaded_by'] =  sru.first_name + ' ' + sru.last_name
        if bundle.obj.library_contents_loading_activity:
            sru = bundle.obj.library_contents_release_activity.activity.performed_by
            bundle.data['released_by'] =  sru.first_name + ' ' + sru.last_name
        return bundle
        
    def build_schema(self):
        schema = super(LibraryContentsVersionResource,self).build_schema()
        return schema
    
    def obj_create(self, bundle, **kwargs):
        bundle.data['date_created'] = timezone.now()
        
        bundle.data['version'] = 1
        super(LibraryContentsVersionResource, self).obj_create(bundle, **kwargs)
    

# class BasicAuthenticationAjaxBrowsers(BasicAuthentication):
#     '''
#     Solves the issue:
#     The session key may not be timed out, but the browser has cleared the 
#     basic-auth credentials: the Django templates use session auth and report 
#     that the user is logged in, but when an ajax request is made, the server
#     asks for basic-auth credentials, and the browser has already cleared them.
#     
#     see: 
#     http://sysadminpy.com/programming/2011/11/14/ajax-and-tastypie---check-if-a-user-has-authenticated/
#     '''
#     
#     
#     def __init__(self, *args, **kwargs):
#         super(BasicAuthenticationAjaxBrowsers, self).__init__(*args, **kwargs)
#  
#     def is_authenticated(self, request, **kwargs):
#         from django.contrib.sessions.models import Session
#         if 'sessionid' in request.COOKIES:
#             s = Session.objects.get(pk=request.COOKIES['sessionid'])
#             if '_auth_user_id' in s.get_decoded():
#                 u = User.objects.get(id=s.get_decoded()['_auth_user_id'])
#                 request.user = u
#                 return True
#         return super(BasicAuthenticationAjaxBrowsers, self).is_authenticated(request, **kwargs)

