import math
import os.path
import sys
import logging
from collections import defaultdict, OrderedDict
import cStringIO
import json
from __builtin__ import StopIteration
import csv
from wsgiref.util import FileWrapper
from zipfile import ZipFile
import time
import shutil
from copy import deepcopy
import re
import urllib2
import io

from django.conf.urls import url
from django.conf import settings
from django.forms.models import model_to_dict
from django.http import Http404
from django.db import connection
from django.db import transaction
from django.contrib.auth.models import User
from django.core.cache import cache
from django.db.models.aggregates import Max, Min
from django.http.response import StreamingHttpResponse, HttpResponse
import django.db.models.sql.constants
import django.db.models.constants
from django.core.urlresolvers import resolve

from PIL import Image

from tastypie.validation import Validation
from tastypie.utils import timezone
from tastypie.exceptions import BadRequest, ImmediateHttpResponse
from tastypie.utils.urls import trailing_slash
from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication,\
    MultiAuthentication
from tastypie.constants import ALL_WITH_RELATIONS
from tastypie import fields
from tastypie.resources import Resource
from tastypie.utils.mime import build_content_type

from sqlalchemy import select, asc
from sqlalchemy import text
from sqlalchemy.sql.expression import column, join
from sqlalchemy.sql import func
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.sql import asc,desc, alias, Alias
from sqlalchemy.sql.expression import nullsfirst,nullslast

from db.models import ScreensaverUser,Screen, LabHead, LabAffiliation, \
    ScreeningRoomUser, ScreenResult, DataColumn, Library, Plate, Copy,\
    PlateLocation, Reagent, Well, LibraryContentsVersion, Activity,\
    AdministrativeActivity, SmallMoleculeReagent, SilencingReagent, GeneSymbol,\
    NaturalProductReagent, Molfile, Gene, GeneGenbankAccessionNumber
from db.support import lims_utils
from reports.serializers import CursorSerializer, LimsSerializer, XLSSerializer
from reports.models import MetaHash, Vocabularies, ApiLog
from reports.api import ManagedModelResource, ManagedResource, ApiLogResource,\
    UserGroupAuthorization, ManagedLinkedResource, log_obj_update,\
    UnlimitedDownloadResource, StreamingResource
from reports.utils.sqlalchemy_bridge import Bridge
from reports.serializers import LIST_DELIMITER_XLS
from reports.serializers import csv_convert
from sqlalchemy.sql.elements import literal_column
from django.core.serializers.json import DjangoJSONEncoder
import hashlib

        
logger = logging.getLogger(__name__)

LIST_DELIMITER_SQL_ARRAY = ','
LIST_DELIMITER_URL_PARAM = ','
MAX_ROWS_PER_XLS_FILE = 100000
MAX_IMAGE_ROWS_PER_XLS_FILE = 2000

def _get_raw_time_string():
  return timezone.now().strftime("%Y%m%d%H%M%S")
    
class ScreensaverUserResource(ManagedModelResource):
#    screens = fields.ToManyField('db.api.ScreenResource', 'screens', 
# related_name='lab_head_id', blank=True, null=True)

    version = fields.IntegerField(attribute='version', null=True)
    administratoruser = fields.ToOneField(
        'db.api.ScreensaverUserResource', 
        attribute='administratoruser', null=True, blank=True)
    screeningroomuser = fields.ToOneField(
        'db.api.ScreensaverUserResource', 
        'screeningroomuser', null=True, blank=True)
    permissions = fields.ToManyField(
        'reports.api.PermissionResource', 'permissions', null=True)
    
    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['digested_password']
        detail_uri_name = 'screensaver_user_id'
        resource_name = 'screensaveruser'
        max_limit = 10000
        
    def __init__(self, **kwargs):
        super(ScreensaverUserResource,self).__init__( **kwargs)
  
    def dehydrate(self, bundle):
        bundle = super(ScreensaverUserResource, self).dehydrate(bundle);
        bundle.data['screens'] = [ x.facility_id 
            for x in Screen.objects.filter(
                lab_head_id=bundle.obj.screensaver_user_id)]
        return bundle        
      
    def apply_sorting(self, obj_list, options):
        options = options.copy()
        options['non_null_fields'] = ['screensaver_user_id']
        obj_list = super(ScreensaverUserResource, self).apply_sorting(
            obj_list, options)
        return obj_list
    
    def apply_filters(self, request, applicable_filters):
        logger.info(str(('apply_filters', applicable_filters)))
        
        return super(ScreensaverUserResource, self).apply_filters(request, 
            applicable_filters)

    def build_filters(self, filters=None):
        logger.info(str(('build_filters', filters)))
        
        return super(ScreensaverUserResource, self).build_filters(filters)
              
    def build_schema(self):
        schema = super(ScreensaverUserResource,self).build_schema()
        schema['idAttribute'] = ['screensaver_user_id']
        return schema
    
    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

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
    
class ScreenResultResource(ManagedResource):

    class Meta:
        queryset = ScreenResult.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'screenresult'
        
        ordering = []
        filtering = {}
        serializer = CursorSerializer()
        allowed_methods = ['get']

        object_class = dict
        max_limit = 10000
        
    def __init__(self, **kwargs):
        self.scope = 'fields.screenresult'
        super(ScreenResultResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" 
        # [ any word, except "schema" ]
        # also note the double underscore "__" is because we also don't want to
        # match in the first clause. Don't want "schema" since that reserved
        # word is used by tastypie for the schema definition for the resource
        return [
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<facility_id>((?=(schema))__|(?!(schema))[^/]+))%s$" )
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<facility_id>((?=(schema))__|(?!(schema))[^/]+))/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
        ]
        
    def get_object_list(self, request):
        logger.warn('Screen result listing not implemented')
        raise Http404(str(('Screen result listing not implemented',
                           request.path)))
        
    def obj_get_list(self, request=None, **kwargs):
        logger.info(unicode(('============= obj_get_list: kwargs', kwargs)))
        # Filtering disabled for brevity...
        return self.get_object_list(request)
    

    def apply_sorting(self, obj_list, options):
        options = options.copy()
        logger.info(str(('=======apply_sorting', options)))
        
#    @staticmethod
#    def get_request_param(param_name, querydict):
    
    def obj_get(self, request=None, **kwargs):
        logger.info(unicode(('============= obj_get: kwargs', kwargs)))
        
        if('bundle' in kwargs and hasattr(kwargs['bundle'].request, 'GET')):
            kwargs.update(kwargs['bundle'].request.GET)
            
            if 'limit' in kwargs:
                limit = kwargs['limit']
                # TODO: why are request parameters being wrapped as lists?
                if not isinstance(limit, (str, unicode)):  
                    # try it as a seq
                    limit = limit[0]
                try:
                    kwargs['limit'] = int(limit)
                except ValueError:
                    raise BadRequest(
                        ("Invalid limit '%s' provided. "
                         "Please provide a positive integer.") 
                         % kwargs['limit'])

            if 'offset' in kwargs:
                offset = kwargs['offset']
                # TODO: why are request parameters being wrapped as lists?
                if not isinstance(offset, (str, unicode)):  
                    # try it as a seq
                    offset = offset[0]
                try:
                    kwargs['offset'] = int(offset)
                except ValueError:
                    raise BadRequest(
                        ("Invalid offset '%s' provided. "
                         "Please provide a positive integer." )
                         % kwargs['offset'])

        if 'facility_id' not in kwargs:
            raise Http404(unicode(('no facility id given',kwargs, request)))
            

        facility_id = kwargs['facility_id']

        try:
            screenresult = ScreenResult.objects.get(
                screen__facility_id=facility_id)            

            # TODO: Not sure what to return for get_obj_list ? since in the 
            # CursorSerializer, you can see that we have to either look 
            # at the passed arg as a cursor, or as a bundle with the cursor as 
            # the obj.  how does TP want this done?
            if screenresult:
                result = {}
                result['meta'] = kwargs.copy();
                result['meta']['total_count'] = self.get_total_count(
                    screenresult)
                
                result['objects'] = self.get_screenresult_cursor(
                    screenresult, **kwargs)
                return result;
            else:
                raise Http404(unicode((
                    'no results for the screen: ', facility_id)))
        
        except Screen.DoesNotExist, e:
            logger.error(str(('no screen found for facility id', facility_id)))
            raise e

        return self.get_object_list(request)

    @staticmethod
    def create_query(screenresult):
        '''    
        select well_id, 
            ...other_entity_columns,  
            (select numeric_value as col1 
             from result_value rv1 
             where data_column_id=4754 and rv1.well_id = w.well_id ) as col1 
        from assay_well w 
        join screen_result using (screen_result_id) 
        join screen using (screen_id) 
        where facility_id = '1003';
        '''
                
        sql = 'select well_id '
        for i,dc in enumerate(
                DataColumn.objects.filter(screen_result=screenresult)):
            column_to_select = None
            if(dc.data_type == 'Numeric'): #TODO: use controlled vocabulary
                column_to_select = 'numeric_value'
            else:
                column_to_select = 'value'

            sql +=  (
                ",(SELECT {col} FROM result_value {alias} "
                "  where {alias}.data_column_id={dc_id} "
                "  and {alias}.well_id=w.well_id) as {column_name} " ).format(
                    col = column_to_select, 
                    alias = "dp_"+str(dc.data_column_id), 
                    dc_id = str(dc.data_column_id), 
                    column_name = "col_"+str(dc.data_column_id) )
        sql += ' FROM assay_well w where w.screen_result_id=%s '
        return sql
    
    def get_total_count(self, screenresult, **kwargs):
        cursor = connection.cursor()
        
        sql = self.create_query(screenresult, **kwargs);
        sql = 'select count(*) from (' + sql + ') a'
        cursor.execute(sql, [screenresult.screen_result_id])
        return cursor.fetchone()[0]
        
    def get_screenresult_cursor(
            self, screenresult, limit=25, offset=0, order_by=[], **kwargs):

        logger.info(unicode((
            '---get_screenresult_cursor', 'limit, offset, order_by', limit, 
            offset, order_by, 'kwargs', kwargs)))
        
        sql = self.create_query(screenresult)
        
        if len(order_by) > 0:
            # TODO: postgres only 
            orderings = map(lambda x:(
                x[1:] + ' DESC NULLS LAST' if x[0]=='-' 
                else x + ' ASC NULLS FIRST' ), 
                order_by )
            sql = ( 'SELECT * FROM ( ' + sql + ') as order_inner ORDER BY ' + ', '.join(orderings) )
                     
        sql += ' OFFSET ' + str(offset)
        sql += ' LIMIT ' + str(limit)
        
        logger.info(str(('sql',sql)))
        cursor = connection.cursor()
        cursor.execute(sql, [screenresult.screen_result_id])
        return cursor

    def get_schema(self, request, **kwargs):
        if not 'facility_id' in kwargs:
            raise Http404(unicode((
                'The screenresult schema requires a screen facility ID'
                ' in the URI, as in /screenresult/[facility_id]/schema/')))
        facility_id = kwargs.pop('facility_id')
        try:
            screenresult = ScreenResult.objects.get(
                screen__facility_id=facility_id)
            logger.info(str(('screenresult resource', 
                             facility_id,screenresult.screen)))
            
            # TODO: Not sure what to return for get_obj_list ? since in the
            # CursorSerializer, you can see that we have to either look 
            # at the passed arg as a cursor, or as a bundle with the cursor as 
            # the obj.  how does TP want this done?
            
            if screenresult:
                return self.create_response(request, 
                                            self.build_schema(screenresult))
            else:
                raise Http404(unicode((
                    'no results for the screen: ', facility_id)))
        except Screen.DoesNotExist, e:
            raise Http404(unicode((
                'no screen found for facility id', facility_id)))
            
    def build_schema(self, screenresult=None):
        logger.debug(str(('==========build schema for screen result', screenresult)))
        data = super(ScreenResultResource,self).build_schema()
        
        if screenresult:
            # now the datacolumn fields
            field_defaults = {
                'visibility': ['list','detail'],
                'ui_type': 'string',
                'type': 'string',
                'filtering': True,
                }
            for i,dc in enumerate(
                    DataColumn.objects.filter(screen_result=screenresult)):
                alias = "dp_"+str(dc.data_column_id)
                columnName = "col_"+str(dc.data_column_id)
                _dict = field_defaults.copy()
                _dict.update(model_to_dict(dc))
                
                _dict['title'] = dc.name
                _dict['comment'] = dc.comments
                _dict['key'] = columnName
                # so that the value columns come last
                _dict['ordinal'] += len(self.fields) + dc.ordinal 
    #            if dc.data_type == 'Numeric':
    #                _dict['ui_type'] = 'numeric'
                
                data['fields'][columnName] = _dict
            # TODO: get the data columns; convert column aliases to real names
        return data
    

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


class DataColumnResource(ManagedModelResource):
    # included to allow querying like ?screen__facility_id=##
    screen = fields.ToOneField('db.api.ScreenResource', 'screen_result__screen')  
    facility_id = fields.CharField('screen_result__screen__facility_id')
    
    class Meta:
        queryset = DataColumn.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(
            BasicAuthentication(), SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'datacolumn'
        
        ordering = []
        filtering = { 'screen': ALL_WITH_RELATIONS}
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
#        self.
        super(DataColumnResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" 
        # [ any word, except "schema" ]
        # also note the double underscore "__" is because we also don't want to
        # match in the first clause. Don't want "schema" since that reserved
        # word is used by tastypie for the schema definition for the resource
        return [
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<data_column_id>((?=(schema))__|(?!(schema))[^/]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

class ScreenResource(ManagedModelResource):

#    lab_head_full = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  
#             full=True) #, full_list=False) #, blank=True, null=True)
    lab_head_link = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  
        full=True)
    lead_screener_link = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  
        full=True)
    
    lab_head_id = fields.IntegerField(attribute='lab_head_id');
    lead_screener_id = fields.IntegerField(attribute='lead_screener_id');
    
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
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" 
        # [ any word, except "schema" ]
        # also note the double underscore "__" is because we also don't want to
        # match in the first clause. Don't want "schema" since that reserved
        # word is used by tastypie for the schema definition for the resource
        return [
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<facility_id>((?=(schema))__|(?!(schema))[^/]+))%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
                    
    def dehydrate(self, bundle):
        if bundle.obj.lead_screener:
            sru = bundle.obj.lead_screener.screensaver_user
            bundle.data['lead_screener'] =  sru.first_name + ' ' + sru.last_name
        if bundle.obj.lab_head:
            lh = bundle.obj.lab_head.screensaver_user.screensaver_user
            bundle.data['lab_head'] =  lh.first_name + ' ' + lh.last_name
        # TODO: the status table does not utilize a primary key, thus it is 
        # incompatible with the standard manager
        #        status_item = ScreenStatusItem.objects.filter(
        #               screen=bundle.obj).order_by('status_date')[0]
        #        bundle.data['status'] = status_item.status
        #        bundle.data['status_date'] = status_item.status_date

        bundle.data['has_screen_result'] = False
        try:
            bundle.data['has_screen_result'] = bundle.obj.screenresult != None
        except ScreenResult.DoesNotExist, e:
            logger.debug(str(('no screenresult for ', bundle.obj)))
        return bundle
    
    def build_schema(self):
        schema = super(ScreenResource,self).build_schema()
        temp = [ x.screen_type for x in self.Meta.queryset.distinct('screen_type')]
        schema['extraSelectorOptions'] = { 
            'label': 'Type', 'searchColumn': 'screen_type', 'options': temp }
        return schema


    def apply_sorting(self, obj_list, options):
        options = options.copy()
        logger.info(str(('options', options)))
        
        extra_order_by = []
        order_by = options.getlist('order_by',None)
        if order_by:
            logger.info(str(('order_by',order_by)))
            for field in order_by:
                temp = field
                dir=''
                if field.startswith('-'):
                    dir = '-'
                    field = field[1:]
                if field == 'lead_screener':
                    order_by.remove(temp)
                    extra_order_by.append(
                        dir+'lead_screener__screensaver_user__last_name')
                    extra_order_by.append(
                        dir+'lead_screener__screensaver_user__first_name')
                if field == 'lab_head':
                    order_by.remove(temp)
                    extra_order_by.append(
                        dir+'lab_head__screensaver_user__screensaver_user__last_name')
                    extra_order_by.append(
                        dir+'lab_head__screensaver_user__screensaver_user__first_name')
                if field == 'has_screen_result':
                    order_by.remove(temp)
                    obj_list = obj_list.extra({
                        'screenresult_isnull': (
                            '(select sr.screen_id is null '
                            'from screen_result sr where sr.screen_id = screen.screen_id) ')})
                    is_null_dir = '-'
                    if dir == '-': is_null_dir = ''
                    extra_order_by.append(is_null_dir+'screenresult_isnull')
            if len(order_by) > 0:
                options.setlist('order_by', order_by)
            else:
                del options['order_by'] 
        logger.info(str(('options',options)))
        obj_list = super(ScreenResource, self).apply_sorting(obj_list, options)
        
        if len(extra_order_by)>0:
            logger.info(str(('extra_order_by', extra_order_by)))
            obj_list = obj_list.order_by(*extra_order_by)
        return obj_list
    
    def hydrate(self, bundle):
        bundle = super(ScreenResource, self).hydrate(bundle);
        return bundle

    def obj_create(self, bundle, **kwargs):
        bundle.data['date_created'] = timezone.now()
        
#         key = 'total_plated_lab_cherry_picks'
#         if key not in bundle.data:
#             field_def = self.get_field_def(key)
#             bundle.data['total_plated_lab_cherry_picks'] = int(field_def['default'])
        bundle.data['version'] = 1
            
        return super(ScreenResource, self).obj_create(bundle, **kwargs)
    
    def save(self, bundle, skip_errors=False):
        ''' returns bundle
        '''
        return super(ScreenResource, self).save(bundle, skip_errors=skip_errors)

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


class LibraryCopyResource(ManagedModelResource):

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
                 r"/(?P<name>[^/]+)%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<library__short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/(?P<name>[^/]+)"
                 r"/plate%s$" ) 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_librarycopyplateview'), 
                name="api_dispatch_librarycopy_plateview"),
        ]    

    def dispatch_librarycopyplateview(self, request, **kwargs):
      
        kwargs['library_short_name'] = kwargs.pop('library__short_name')  
        kwargs['copy_name'] = kwargs.pop('name')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)    
        
    def get_object_list(self, request, library_short_name=None):
        ''' 
        Called immediately before filtering, actually grabs the (ModelResource) 
        query - 
        
        Override this and apply_filters, so that we can control the extra 
        column "is_for_group":
        This extra column is present when navigating to permissions from a 
        usergroup; see prepend_urls.
        '''
        query = super(LibraryCopyResource, self).get_object_list(request);
        logger.info(str(('get_obj_list', len(query))))
        if library_short_name:
            query = query.filter(library__short_name=library_short_name)
        return query
    
        
                    
    def apply_sorting(self, obj_list, options):
        options = options.copy()
        logger.info(str(('options', options)))
        
        extra_order_by = []
        order_by = options.getlist('order_by',None)
        if order_by:
            logger.info(str(('order_by',order_by)))
            for field in order_by:
                temp = field
                dir=''
                if field.startswith('-'):
                    dir = '-'
                    field = field[1:]
                if field == 'created_by':
                    order_by.remove(temp)
                    extra_order_by.append(dir+'created_by__last_name')
                    extra_order_by.append(dir+'created_by__first_name')
            if len(order_by) > 0:
                options.setlist('order_by', order_by)
            else:
                del options['order_by'] 

        obj_list = super(LibraryCopyResource, self).apply_sorting(obj_list, options)
        
        if len(extra_order_by)>0:
            logger.info(str(('extra_order_by', extra_order_by)))
            obj_list = obj_list.order_by(*extra_order_by)
        return obj_list

    def dehydrate(self, bundle):
        if bundle.obj.created_by:
            user = bundle.obj.created_by
            bundle.data['created_by'] =  user.first_name + ' ' + user.last_name
        return bundle
    
    def build_schema(self):
        # FIXME: these options should be defined automatically from a vocabulary in build_schema
        schema = super(LibraryCopyResource,self).build_schema()
#         temp = [ x.usage_type for x in self.Meta.queryset.distinct('usage_type')]
#         schema['extraSelectorOptions'] = { 
#             'label': 'Type', 'searchColumn': 'usage_type', 'options': temp }
        return schema
    
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
        
class LibraryCopyPlateResource(ManagedModelResource):

    library_short_name = fields.CharField('copy__library__short_name',  null=True)
    copy_name = fields.CharField('copy__name',  null=True)
    plate_location = fields.ToOneField('db.api.PlateLocationResource', 
                                        attribute='plate_location', 
                                        full=True, full_detail=True, full_list=True,
                                        null=True)
    
    # TODO:
    # status_date
    
    # plate_location = 
    
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
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<copy__library__short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/(?P<copy__name>[^/]+)"
                 r"/(?P<plate_number>[^/]+)%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]    

    def get_object_list(self, request, library_short_name=None, copy_name=None):
        ''' 
        Called immediately before filtering, actually grabs the (ModelResource) 
        query - 

        Note: any extra kwargs are there because we are injecting them into the 
        global TP kwargs in one of the various "dispatch_" handlers assigned 
        through prepend_urls.  Here we can explicitly add them to the query. 
        
        '''
        query = super(LibraryCopyPlateResource, self).get_object_list(request);
        if library_short_name:
            query = query.filter(copy__library__short_name=library_short_name)
        if copy_name:
            query = query.filter(copy_name=copy_name)
        return query
                    
    def apply_sorting(self, obj_list, options):
        options = options.copy()
        logger.info(str(('options', options)))
        
        extra_order_by = []
        
        # handle joined table sorts
        order_by = options.getlist('order_by',None)
        if order_by:
            for field in order_by:
                if field == 'copy_name':
                    temp = field
                    _dir=''
                    if field.startswith('-'):
                        _dir = '-'
                        field = field[1:]
                    order_by.remove(temp)
                    extra_order_by.append(_dir+'copy__name')
            if len(order_by) > 0:
                options.setlist('order_by', order_by)
            else:
                del options['order_by'] 
        obj_list = super(LibraryCopyPlateResource, self).apply_sorting(
            obj_list, options)
        
        if len(extra_order_by)>0:
            logger.info(str(('extra_order_by', extra_order_by)))
            obj_list = obj_list.order_by(*extra_order_by)
        return obj_list

    def dehydrate(self, bundle):
        if bundle.obj.created_by:
            user = bundle.obj.created_by
            bundle.data['created_by'] =  user.first_name + ' ' + user.last_name
        return bundle
    
    
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

 
# class GeneResource(ManagedLinkedResource):
# 
#     class Meta:
#         queryset = Gene.objects.all() 
#         authentication = MultiAuthentication(
#             BasicAuthentication(), SessionAuthentication())
#         authorization= UserGroupAuthorization()
# 
#         ordering = []
#         filtering = {}
#         serializer = LimsSerializer()
#         excludes = [] #['json_field']
#         always_return_data = True # this makes Backbone happy
#         resource_name='gene' 
# 
#     def __init__(self, **kwargs):
#         super(SmallMoleculeReagentResource,self).__init__(**kwargs)


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
        
        #FIXME: this is a hack for SIRNA, to handle the gene table
        - sirna->(vendor,facility)gene->gene_genbank_accession_number
        '''
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
                
                logger.info(str((select_stmt)))
                
                
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
        bytes = cStringIO.StringIO()
        
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

class SqlAlchemyResource(StreamingResource):
    '''
    FIXME: serialization of list fields: they are sent as concatenated strings with commas
    FIXME: Reagent/SMR specific
    '''
    
    class Meta:
        queryset = Reagent.objects.all()
#         serializer = LimsSerializer() # still have a serializer for error response

    
    def __init__(self, *args, **kwargs):
        # get a handle to the SqlAlchemy "bridge" for its table registry and 
        # connection handle
        self.bridge = Bridge()
        
        super(SqlAlchemyResource, self).__init__(*args, **kwargs)


    def build_sqlalchemy_columns(self, fields, base_query_tables=[], already_defined_columns={}):
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

        try:
            columns = OrderedDict()
            if already_defined_columns:
                columns = OrderedDict(already_defined_columns)
            for field in fields:
                if field['key'] in columns:
                    continue
                
                logger.info(str(('build column', field['key'])))
                field_name = field.get('field', None)
                if not field_name:
                    field_name = field['key']
                
                field_table = field.get('table', None)
                
                label = field['key']
                
                if not field_table:
                    continue
                logger.info(str(('field', field['key'], 'field_table', field_table )))
                if field_table in base_query_tables:
                    # simple case: table.fields already selected in the base query:
                    # just need to specify them
                    if field_name in self.bridge[field_table].c:
                        col = self.bridge[field_table].c[field_name]
                    else:
                        raise Exception(str(('field', field_name, 'not found in table', field_table)))
                    col = col.label(label)
                    columns[label] = col
                    
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
                        join_stmt = join_stmt.label(label)
                        columns[label] = join_stmt
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
                        stmt2 = stmt2.select_from(join_stmt).label(label)
                        columns[label] = stmt2
                
            logger.info(str(('columns', columns.keys())))
            return columns
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on build sqlalchemy columns', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e   

    def build_sqlalchemy_ordering(self, request):
        '''
        returns a scalar or list of ClauseElement objects which will comprise 
        the ORDER BY clause of the resulting select.
        This method borrows from tastypie.resources.ModelResource.apply_sorting
        '''
        options = request.GET
        parameter_name = 'order_by'
        if not parameter_name in options:
            return []
        
        order_by_args = []

        if hasattr(options, 'getlist'):
            order_bits = options.getlist(parameter_name)
        else:
            order_bits = options.get(parameter_name)

            if not isinstance(order_bits, (list, tuple)):
                order_bits = [order_bits]
        
        order_clauses = []
        for order_by in order_bits:
            field_name = order_by
            if order_by.startswith('-'):
                field_name = order_by[1:]
                order_clauses.append(nullslast(desc(column(field_name))))
            else:
                order_clauses.append(nullsfirst(asc(column(field_name))))
                
        return order_clauses
    
    def filter_value_to_python(self, value, filters, filter_expr, filter_type):
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
            if hasattr(filters, 'getlist'):
                value = []

                for part in filters.getlist(filter_expr):
                    value.extend(part.split(LIST_DELIMITER_URL_PARAM))
            else:
                value = value.split(LIST_DELIMITER_URL_PARAM)

        return value

    def build_sqlalchemy_filters(self, schema, request, **kwargs):
        '''
        Attempt to create a SqlAlchemy whereclause out of django style filters.
        
        This method borrows from tastypie.resources.ModelResource.build_filters
        
        Valid values are either a list of Django filter types (i.e.
        ``['startswith', 'exact', 'lte']``), the ``ALL`` constant or the
        ``ALL_WITH_RELATIONS`` constant.
        
        @return - (sqlalchemy.whereclause, [field_name_list])
        '''
        
        #     filtering = {
        #         'resource_field_name': ['exact', 'startswith', 'endswith', 'contains'],
        #         'resource_field_name_2': ['exact', 'gt', 'gte', 'lt', 'lte', 'range'],
        #         'resource_field_name_3': ALL,
        #         'resource_field_name_4': ALL_WITH_RELATIONS,
        #         ...
        #     }
        query_terms = django.db.models.sql.constants.QUERY_TERMS
        lookup_sep = django.db.models.constants.LOOKUP_SEP

        filters = request.GET.copy()
        # Update with the provided kwargs.
        filters.update(kwargs)

        if filters is None:
            return (None,None)
        
        expressions = []
        filtered_fields = []
        for filter_expr, value in filters.items():
            filter_bits = filter_expr.split(lookup_sep)
            field_name = filter_bits.pop(0)
            filter_type = 'exact'

            if not field_name in schema['fields']:
                continue
            filtered_fields.append(field_name)
            field = schema['fields'][field_name]
            
            if len(filter_bits) and filter_bits[-1] in query_terms:
                filter_type = filter_bits.pop()

            value = self.filter_value_to_python(value, filters, filter_expr, filter_type)

                        # TODO all the types
            #             QUERY_TERMS = set([
            #                 'exact', 'iexact', 'contains', 'icontains', 'gt', 'gte', 'lt', 'lte', 'in',
            #                 'startswith', 'istartswith', 'endswith', 'iendswith', 'range', 'year',
            #                 'month', 'day', 'week_day', 'hour', 'minute', 'second', 'isnull', 'search',
            #                 'regex', 'iregex',
            #             ])

            if filter_type == 'exact':
                expressions.append(column(field_name)==value)
            elif filter_type == 'contains':
                expressions.append(column(field_name).contains(value))
            elif filter_type == 'icontains':
                expressions.append(column(field_name).ilike('%{value}%'.format(value=value)))
        if len(expressions) > 1: 
            from sqlalchemy.sql import and_, or_, not_          
            return (and_(*expressions), filtered_fields)
        elif len(expressions) == 1:
            return (expressions[0], filtered_fields) 
        else:
            return (None, filtered_fields)

    def get_format(self, request):
        format = request.GET.get('format', None)
        logger.info(str(('format', format)))
        if format:
            if format in self.content_types:
                format = self.content_types[format]
                logger.info(str(('format', format)))
            else:
                logger.error(str(('unknown format', desired_format)))
                raise ImmediateHttpResponse("unknown format: %s" % desired_format)
        else:
            # Try to fallback on the Accepts header.
            if request.META.get('HTTP_ACCEPT', '*/*') != '*/*':
                try:
                    import mimeparse
                    format = mimeparse.best_match(
                        self.content_types.values(), request.META['HTTP_ACCEPT'])
                    logger.info(str(('format', format, request.META['HTTP_ACCEPT'])))
                except ValueError:
                    logger.error(str(('Invalid Accept header')))
                    raise ImmediateHttpResponse('Invalid Accept header')
        return format
        
    def get_list(self, request, **kwargs):
        raise NotImplemented(str(('get_list must be implemented for the SqlAlchemyResource', 
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
#                         logger.info(str((y, 'hash key:', key_digest, m.hexdigest())))
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
            output_filename, field_hash={}):
        limit = request.GET.get('limit', self._meta.limit) 
        try:
            limit = int(limit)
        except ValueError:
            raise BadRequest(
                "Invalid limit '%s' provided. Please provide a positive integer." % limit)
        if limit > 0:    
            stmt = stmt.limit(limit)

        offset = request.GET.get('offset', 0) 
        try:
            offset = int(offset)
        except ValueError:
            raise BadRequest(
                "Invalid offset '%s' provided. Please provide a positive integer." % offset)
        if offset < 0:    
            offset = -offset
        
        stmt = stmt.offset(offset)
        conn = self.bridge.get_engine().connect()
        
        logger.info(str(('stmt', str(stmt))))
        logger.info(str(('count stmt', str(count_stmt))))
        
        logger.info('excute stmt')
        result = conn.execute(stmt)
        logger.info('excuted stmt')
        
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
                    return ''
            return re.sub(r'{([^}]+)}', get_value_from_template, value_template)

        if desired_format == 'application/json':
            
            cache_hit = self._cached_resultproxy(stmt, count_stmt, limit, offset)
            if cache_hit:
                result = cache_hit['cached_result']
                count = cache_hit['count']
            else:
                # use result from before
                logger.info(str(('execute count')))
                count = conn.execute(count_stmt).scalar()
            logger.info(str(('====count====', count)))
        
            meta = {
                'limit': limit,
                'offset': offset,
                'total_count': count
                }    
                        
            def json_generator(cursor):
#                 logger.info(str(('meta', meta)))
                yield '{ "meta": %s, "objects": [' % json.dumps(meta)
                i=0
                for row in cursor:
#                     logger.info(str(('row', row, row.keys())))
                    if field_hash:
                        _dict = OrderedDict()
                        for key, field in field_hash.iteritems():
#                             logger.info(str(('key', key,  row.has_key(key))))
                            value = None
                            if row.has_key(key):
                                value = row[key]
                            if ( field['json_field_type'] == 'fields.ListField' 
                                 or field['linked_field_type'] == 'fields.ListField' ) and value:
                                # FIXME: need to do an escaped split
#                                 logger.info(str(('split', key, value)))
                                value = value.split(LIST_DELIMITER_SQL_ARRAY)

#                             logger.info(str(('key val', key, value, field)))
                            _dict[key] = value
                            
                            if field.get('value_template', None):
                                value_template = field['value_template']
#                                 logger.info(str(('field', key, value_template)))
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
                                        logger.debug(str(('no image at', newval,e)))
                                else:
                                    _dict[key]=newval
                    else:
#                         logger.info(str(('raw', row)))
                        _dict = dict((x,y) for x, y in row.items())

                    i += 1
                    if i == 1:
                        try:
#                             logger.info(str(('_dict',i, _dict)))
#                             logger.info(str(('_dict', json.dumps(_dict))))
#                             yield json.dumps(_dict)
                            yield json.dumps(_dict, cls=DjangoJSONEncoder,
                                sort_keys=True, ensure_ascii=False, indent=2)

                        except Exception, e:
                            logger.error(str(('ex', _dict, e)))
                    else:
                        try:
#                             logger.info(str(('_dict',i, _dict)))
#                             yield ', ' + json.dumps(_dict)
                            yield ', ' + json.dumps(_dict, cls=DjangoJSONEncoder,
                                sort_keys=True, ensure_ascii=False, indent=2)
                        except Exception, e:
                            logger.error(str(('ex============', _dict, e)))
                
#                 logger.info(str(('i', i)))    
                yield ' ] }'
            
            logger.info(str(('ready to stream 1...')))
            # NOTE: sqlalchemy ResultProxy will automatically close the result after last row fetch
            response = StreamingHttpResponse(ChunkIterWrapper(json_generator(result)))
            response['Content-Type'] = content_type
            return response
        
        elif desired_format == 'application/xls':
            
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
                logger.info(str(('temp file', temp_file)))
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
                                worksheet.write(filerow,col,key)
                                ## TODO: option to write titles
                                # worksheet.write(filerow,col,field['title'])
                        else:
                            for col,name in enumerate(record.keys()):
                                worksheet.write(filerow,col,name)
                        filerow += 1
                    
                    if field_hash:
                        for col, (key, field) in enumerate(field_hash.iteritems()):
                            if row.has_key(key):
                                value = row[key]
                            if ( field['json_field_type'] == 'fields.ListField' 
                                 or field['linked_field_type'] == 'fields.ListField' ) and value:
                                value = value.split(LIST_DELIMITER_SQL_ARRAY)
                            worksheet.write(filerow, col, 
                                csv_convert(value, delimiter=LIST_DELIMITER_XLS))

                            if field.get('value_template', None):
                                value_template = field['value_template']
                                logger.info(str(('field', key, value_template)))
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
                                        worksheet.insert_image(filerow, col, newval, {'image_data': io.BytesIO(response.content)})
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

                if len(file_names_to_zip) >1:
                    # create a temp zip file
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
                    wrapper, content_type='application/zip; charset=utf-8') 
                response['Content-Disposition'] = \
                    'attachment; filename=%s' % filename
                response['Content-Length'] = os.path.getsize(temp_file)
                return response
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

                        if MOLDATAKEY in row and row[MOLDATAKEY]:
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
#                                     logger.info(str(('field', key, value_template)))
                                    value = interpolate_value_template(value_template, row)
                                if ( field['json_field_type'] == 'fields.ListField' 
                                     or field['linked_field_type'] == 'fields.ListField' ) and value:
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
                return response
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
                pseudo_buffer, delimiter=',', quotechar='"', 
                quoting=csv.QUOTE_ALL, lineterminator="\n")
            def csv_generator(cursor):
                i=0
                for row in cursor:
                    i += 1
                    if i == 1:
                        if field_hash:
                            yield csvwriter.writerow(field_hash.keys())
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
                            if ( field['json_field_type'] == 'fields.ListField' 
                                 or field['linked_field_type'] == 'fields.ListField' ) and value:
                                # FIXME: must quote special strings?
#                                 value = '[' + ",".join(value.split(LIST_DELIMITER_SQL_ARRAY)) + ']'
                                value = value.split(LIST_DELIMITER_SQL_ARRAY)
                            
#                             if field.get('data_type', None):
#                                 data_type = field['data_type']
#                                 if data_type == "float":
#                                     value = 
                                
                            if field.get('value_template', None):
                                value_template = field['value_template']
                                logger.info(str(('field', key, value_template)))
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
                        yield csvwriter.writerow([csv_convert(val) for val in values])
                    
                    else:
                        yield csvwriter.writerow([csv_convert(val) for val in row.values()])

            response = StreamingHttpResponse(csv_generator(result),
                content_type=content_type)
            name = self._meta.resource_name
            response['Content-Disposition'] = \
                'attachment; filename=%s.csv' % output_filename
            return response

        else:
            msg = str(('unknown format', desired_format, output_filename))
            logger.error(msg)
            raise ImmediateHttpResponse(msg)

class LibraryCopyPlatesResource(SqlAlchemyResource, ManagedModelResource):
    class Meta:
        queryset = Plate.objects.all().order_by('name')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'librarycopyplates'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

    def __init__(self, **kwargs):
        super(LibraryCopyPlatesResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # Note: because this prepends the other list, we have to make sure 
        # "schema" is matched
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)/(?P<library_short_name>[\w\d_.\-\+: ]+)/(?P<name>[\w\d_.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]

    def get_list(self,request,**kwargs):
        try:
            _p = self.bridge['plate']
            _pl = self.bridge['plate_location']
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _ap = self.bridge['assay_plate']
            
            library = None
            copy = None
            plate = None
            if 'library_short_name' in kwargs:
                library = Library.objects.get(short_name=kwargs['library_short_name'])
                filename = '%s_%s' % (self._meta.resource_name, library.short_name )
                if 'name' in kwargs:
                    copy = Copy.objects.get(name=kwargs['name'])
                    filename = '%s_%s_%s' % (self._meta.resource_name, library.short_name, copy.name )
            else:
                filename = '%s' % (self._meta.resource_name )
                logger.info(str(('no library_short_name provided')))
    

#         # define plate_screening_statistics subquery        
#         text_cols = ','.join([
#             'p.plate_id', 
#             'c.copy_id',
#             'c.name',
#             'l.short_name',
#             'count(distinct(ls)) screening_count', 
#             'count(distinct(ap)) ap_count', 
#             'count(distinct(srdl)) dl_count',
#             'min(srdl.date_of_activity) first_date_data_loaded', 
#             'max(srdl.date_of_activity) last_date_data_loaded', 
#             'min(lsa.date_of_activity) first_date_screened', 
#             'max(lsa.date_of_activity) last_date_screened', 
#             ])
#         text_select = '\n'.join([
#             'copy c', 
#             'join plate p using(copy_id) ',
#             'join library l using(library_id) ',
#             'join assay_plate ap on(ap.plate_id=p.plate_id)',  
#             'left outer join library_screening ls on(library_screening_id=ls.activity_id)',
#             'left outer join activity lsa on(ls.activity_id=lsa.activity_id)', 
#             'left outer join (',
#             '    select a.activity_id, a.date_of_activity', 
#             '    from activity a', 
#             '    join administrative_activity using(activity_id) ) srdl', 
#             '    on(screen_result_data_loading_id=srdl.activity_id)'  
#             ])
# 
#         plate_screening_statistics = select([text(text_cols)]).\
#             select_from( text(text_select))
#         if library:
#             plate_screening_statistics = \
#                 plate_screening_statistics.where(
#                     _l.c.library_id == library.library_id )
#         plate_screening_statistics = plate_screening_statistics.group_by(
#             text('p.plate_id, c.copy_id, c.name, l.short_name '))
#         plate_screening_statistics = \
#             plate_screening_statistics.order_by(text('p.plate_id '))
#         plate_screening_statistics = plate_screening_statistics.cte('plate_screening_statistics')
        
        
            # NOTE: precalculated version
            plate_screening_statistics = \
                select([text('*')]).\
                    select_from(text('plate_screening_statistics')).cte('plate_screening_statistics')
        


            p1 = plate_screening_statistics.alias('p1')
            
            j = join(_p, _c, _p.c.copy_id == _c.c.copy_id )
            j = j.join(p1, _p.c.plate_id == text('p1.plate_id'), isouter=True)
            j = j.join(_pl, _p.c.plate_location_id == _pl.c.plate_location_id )
            j = j.join(_l, _c.c.library_id == _l.c.library_id )
        
            logger.info(str(('get_list', kwargs)))
            
            schema = super(LibraryCopyPlatesResource,self).build_schema()
        
            # FIXME: 20150114 - includes not being sent by UI
            includes = request.GET.getlist('includes', None)
            logger.info(str(('includes', includes)))
            if includes:
                manual_field_includes = set(includes)
            else:    
                manual_field_includes = set()
            logger.info(str(('manual_field_includes', manual_field_includes)))

            (filter_expression, filter_fields) = \
                self.build_sqlalchemy_filters(schema, request, **kwargs)
            
            temp = { key:field for key,field in schema['fields'].items() 
                if ((field.get('visibility', None) and 'list' in field['visibility']) 
                    or field['key'] in filter_fields 
                    or field['key'] in manual_field_includes )}
            # manual excludes
            temp = { key:field for key,field in temp.items() 
                if '-%s' % key not in manual_field_includes }
            
            field_hash = OrderedDict(sorted(temp.iteritems(), key=lambda x: x[1]['ordinal'])) 
            logger.info(str(('field_hash final: ', field_hash.keys() )))
            
            base_query_tables = ['plate', 'copy','plate_location', 'library']
            
            already_defined_columns={
                'screening_count': literal_column('p1.screening_count'), 
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
                'date_plated': literal_column(
                    '(select date_of_activity '
                    ' from activity a'
                    ' where a.activity_id=plate.plated_activity_id )').label('date_plated'),
                'date_retired': literal_column(
                    '(select date_of_activity '
                    ' from activity a'
                    ' where a.activity_id=plate.retired_activity_id )').label('date_retired'),
                    };
                    
            # TODO: test "includes" and this, when working (ui is stripping args from post)
            already_defined_columns = { key:field for key,field in already_defined_columns.iteritems()
                if key in field_hash }
            
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                already_defined_columns=already_defined_columns
                )
            cols = list(already_defined_columns.values())
            cols.extend( [v for k,v in columns.iteritems() if k not in already_defined_columns] )
#             columns['status_date']= already_defined_columns['status_date']
            logger.info(str(('columns', cols)))
            
            stmt = select(cols).select_from(j)
            stmt = stmt.order_by(_l.c.short_name, _c.c.name, _p.c.plate_number )

            order_clauses = self.build_sqlalchemy_ordering(request)
            if order_clauses:
                logger.info(str(('order_clauses', [str(c) for c in order_clauses])))
                _alias = Alias(stmt)
                stmt = select([text('*')]).select_from(_alias)
                stmt = stmt.order_by(*order_clauses)
            
            logger.info(str(('filter_expression', str(filter_expression))))
            if filter_expression is not None:
                if not order_clauses:
                    _alias = Alias(stmt)
                    logger.info(str(('filter_expression', str(filter_expression))))
                    stmt = select([text('*')]).select_from(_alias)
                stmt = stmt.where(filter_expression)

            count_stmt = select([func.count()]).select_from(stmt.alias())

            return self.stream_response_from_cursor(
                request, stmt, count_stmt, filename, field_hash=field_hash  )
        
            
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
                self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e  
        
        
class LibraryCopiesResource(SqlAlchemyResource, ManagedModelResource):
    ''' 
    Testing class for the "freeze copy thaw" reports
    '''

    class Meta:
        queryset = Copy.objects.all().order_by('name')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'librarycopies'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

        
    def __init__(self, **kwargs):
        super(LibraryCopiesResource,self).__init__(**kwargs)

    def prepend_urls(self):
        # Note: because this prepends the other list, we have to make sure 
        # "schema" is matched
        
        return [
            # override the parent "base_urls" so that we don't need to worry about schema again
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),

            url(r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
        ]


    def get_list(self,request,**kwargs):


        try:
            _c = self.bridge['copy']
            _l = self.bridge['library']
            _ap = self.bridge['assay_plate']
            
            library = None
            if 'library_short_name' in kwargs:
                library = Library.objects.get(short_name=kwargs['library_short_name'])
                filename = '%s_%s' % (self._meta.resource_name, library.short_name )
            else:
                filename = '%s' % (self._meta.resource_name )
                logger.info(str(('no library_short_name provided')))
            
            
            # = copy volume statistics =
            
            # REMOVED 20150112 - using precalculated table see 0016 migration
            #         text_cols=[ literal_column('c.copy_id').label('copy_id'),
            # text('c.name'),
            # text('c.short_name'),
            # 'avg(c.plate_remaining_volume) avg_plate_volume', 
            # 'min(c.plate_remaining_volume) min_plate_volume', 
            # 'max(c.plate_remaining_volume) max_plate_volume']
            #         text_select =''' 
            # ( select 
            #     p.copy_id, 
            #     p.well_volume - sum(la.volume_transferred_per_well_from_library_plates) as plate_remaining_volume,
            #     co.name,
            #     l.short_name
            #     from plate p
            #     join copy co using(copy_id)
            #     join library l using(library_id) 
            #     join assay_plate ap using(plate_id) 
            #     join screening ls on(ls.activity_id = ap.library_screening_id) 
            #     join lab_activity la using(activity_id) 
            #     where co.library_id = 108
            #     and ap.replicate_ordinal = 0 
            #     group by p.copy_id, p.plate_id, p.well_volume, co.name, l.short_name ) as c 
            # group by c.copy_id, c.name, c.short_name'''
            
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
            if library:
                copy_volume_statistics = copy_volume_statistics.where(_l.c.library_id == library.library_id )
            copy_volume_statistics = copy_volume_statistics.group_by(text('c.copy_id, c.name, l.short_name '))
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
            if library:
                copy_plate_screening_statistics = copy_plate_screening_statistics.where(_l.c.library_id == library.library_id )
            copy_plate_screening_statistics = \
                copy_plate_screening_statistics.group_by('c.copy_id, c.name, l.short_name')
            copy_plate_screening_statistics = copy_plate_screening_statistics.cte('copy_plate_screening_statistics')
    
            # = copy screening statistics = 
            
            # REMOVED 20150112 - using precalculated stats, see 0016 migrations
            #         text_cols = '''c.copy_id,
            # c.name,
            # l.short_name
            # ,count(distinct(ls)) screening_count 
            # ,count(distinct(ap)) ap_count 
            # ,count(distinct(srdl)) dl_count
            # ,min(srdl.date_of_activity) first_date_data_loaded 
            # ,max(srdl.date_of_activity) last_date_data_loaded 
            # ,min(lsa.date_of_activity) first_date_screened 
            # ,max(lsa.date_of_activity) last_date_screened'''
            #         text_select = \
            #             '''copy c 
            # join plate p using(copy_id) 
            # join library l using(library_id) 
            # join assay_plate ap on(ap.plate_id=p.plate_id)  
            #     left outer join library_screening ls on(library_screening_id=ls.activity_id)
            #     left outer join activity lsa on(ls.activity_id=lsa.activity_id) 
            # left outer join (select a.activity_id, a.date_of_activity from activity a join administrative_activity using(activity_id) ) srdl on(screen_result_data_loading_id=srdl.activity_id)  
            # group by c.copy_id, c.name, l.short_name'''
            # 
            # TODO: dynamic version
            #         copy_screening_statistics = \
            #             select([text(text_cols)]).select_from( text(text_select)).cte('copy_screening_statistics')
            
            # NOTE: precalculated version
            copy_screening_statistics = \
                select([text('*')]).\
                    select_from(text('copy_screening_statistics')).cte('copy_screening_statistics')
    
            # FIXME: 20150114 - itemize colums like librarycopyplates
            
            already_defined_columns = {
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
            
            #         text_cols = ','.join([
            #             'c1.copy_id', 
            #             'c1.short_name library_short_name',
            #             'c1.plate_screening_count',
            #             'c1.copy_plate_count',
            #             ('c1.plate_screening_count::float/c1.copy_plate_count::float '
            #                 'as plate_screening_count_average'),
            #             'c2.avg_plate_volume',
            #             'c2.min_plate_volume', 
            #             'c2.max_plate_volume', 
            #             'c3.screening_count',
            #             'c3.ap_count',
            #             'c3.dl_count',
            #             'c3.first_date_data_loaded', 
            #             'c3.last_date_data_loaded', 
            #             'c3.first_date_screened', 
            #             'c3.last_date_screened',
            #             '\n'.join([
            #                 "( select room || '-' || freezer || '-' || shelf || '-' || bin ", 
            #                 '    from plate_location pl ' ,
            #                 '    where pl.plate_location_id=copy.primary_plate_location_id) ',
            #                 '        as primary_plate_location']),
            #             '\n'.join([
            #                 '(select count(distinct(plate_location_id)) ',
            #                 '    from plate p',
            #                 '    where p.copy_id = copy.copy_id ) ',
            #                 'as plate_locations']),
            #             '\n'.join([
            #                 '(select count(p)', 
            #                 '    from plate p ',
            #                 '    where p.copy_id=copy.copy_id', 
            #                 "    and p.status = 'Available' ) ",
            #                 '    as plates_available']),
            #                 ])

            c1 = copy_plate_screening_statistics.alias('c1')
            c2 = copy_volume_statistics.alias('c2')
            c3 = copy_screening_statistics.alias('c3')
    
            j = join(_c, c1, _c.c.copy_id == text('c1.copy_id'), isouter=True)
            j = j.outerjoin(c2,text('c1.copy_id = c2.copy_id') )
            j = j.outerjoin(c3, text('c1.copy_id = c3.copy_id') )
    
            logger.info(str(('====j', str(j))))

            logger.info(str(('get_list', kwargs)))
            
            schema = super(LibraryCopiesResource,self).build_schema()
    
            # FIXME: 20150114 - includes not being sent by UI
            includes = request.GET.getlist('includes', None)
            logger.info(str(('includes', includes)))
            if includes:
                manual_field_includes = set(includes)
#                 for include in includes:
#                     manual_field_includes.add(include) #.split(LIST_DELIMITER_URL_PARAM))
            else:    
                manual_field_includes = set()
            logger.info(str(('manual_field_includes', manual_field_includes)))

            (filter_expression, filter_fields) = \
                self.build_sqlalchemy_filters(schema, request, **kwargs)
            
            temp = { key:field for key,field in schema['fields'].items() 
                if ((field.get('visibility', None) and 'list' in field['visibility']) 
                    or field['key'] in filter_fields 
                    or field['key'] in manual_field_includes )}
            
            # manual excludes
            temp = { key:field for key,field in temp.items() 
                if '-%s' % key not in manual_field_includes }
            
            field_hash = OrderedDict(sorted(temp.iteritems(), key=lambda x: x[1]['ordinal'])) 
    
            base_query_tables = ['copy','library']
            
#             already_defined_columns={};
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                already_defined_columns=already_defined_columns
                )
            logger.info(str(('columns', columns)))
            cols = list(already_defined_columns.values())
            cols.extend( [v for k,v in columns.iteritems() if k not in already_defined_columns] )
            
#             cols = [text(text_cols)]
#             cols.extend( [v for k,v in columns.iteritems() if k not in already_defined_columns] )
            
            stmt = select(cols).select_from(j)
            stmt = stmt.order_by(text('c1.short_name,c1.name'))

            order_clauses = self.build_sqlalchemy_ordering(request)
            if order_clauses:
                logger.info(str(('order_clauses', [str(c) for c in order_clauses])))
                _alias = Alias(stmt)
                stmt = select([text('*')]).select_from(_alias)
                stmt = stmt.order_by(*order_clauses)
            
            logger.info(str(('filter_expression', str(filter_expression))))
            if filter_expression is not None:
                if not order_clauses:
                    _alias = Alias(stmt)
                    logger.info(str(('filter_expression', str(filter_expression))))
                    stmt = select([text('*')]).select_from(_alias)
                stmt = stmt.where(filter_expression)

            count_stmt = select([func.count()]).select_from(stmt.alias())

            return self.stream_response_from_cursor(
                request, stmt, count_stmt, filename, field_hash=field_hash  )
        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('on get_list', 
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
 
    def get_list(self, request, **kwargs):
    
        if 'library_short_name' in kwargs:
            library = Library.objects.get(short_name=kwargs['library_short_name'])
        else:
            raise NotImplementedError('must provide a library_short_name parameter')
        logger.info(str(('kwargs', kwargs)))
        
        filename = '%s_%s' % (self._meta.resource_name, library.short_name )
        
        # Build the query columns using directions from our schema
        # Specify the tables in the base query (*will not need to re-join them)
        base_query_tables = ['well', 'reagent', 'library']
        schema = self.build_schema(library=library)

        (filter_expression, filter_fields) = \
            self.build_sqlalchemy_filters(schema, request, **kwargs)

        # get manual field includes from kwargs
        includes = request.GET.get('includes', None)
        logger.info(str(('includes', includes)))
        if includes:
            manual_field_includes = set(includes.split(LIST_DELIMITER_URL_PARAM))
        else:    
            manual_field_includes = set()
        logger.info(str(('manual_field_includes', manual_field_includes)))

        desired_format = self.get_format(request)
        #         if desired_format == 'application/xls':
        manual_field_includes.add('structure_image')
        if desired_format == 'chemical/x-mdl-sdfile':
            manual_field_includes.add('molfile')

        # Filter fields and put in an ordered dict
        # TODO: if this is a detail search, include those columns then.
        temp = { key:field for key,field in schema['fields'].items() 
            if ((field.get('visibility', None) and 'list' in field['visibility']) 
                or field['key'] in filter_fields 
                or field['key'] in manual_field_includes )}
        field_hash = OrderedDict(sorted(temp.iteritems(), key=lambda x: x[1]['ordinal'])) 
        
        sub_resource = self.get_reagent_resource(library)
        if hasattr(sub_resource, 'build_sqlalchemy_columns'):
            sub_columns = sub_resource.build_sqlalchemy_columns(
                field_hash.values(), self.bridge)
            logger.info(str(('sub_columns', sub_columns.keys())))
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                already_defined_columns=sub_columns)
        else:
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables)
        
        
        # Start building a query; use the sqlalchemy.sql.schema.Table API:
        _well = self.bridge['well']
        _reagent = self.bridge['reagent']
        _library = self.bridge['library']
        j = _well.join(_reagent, _well.c.well_id==_reagent.c.well_id, isouter=True)
        j = j.join(_library, _well.c.library_id == _library.c.library_id )
        stmt = select(columns.values()).\
            select_from(j).\
            where(_well.c.library_id == library.library_id) 

        # perform ordering and filters     
        
        order_clauses = self.build_sqlalchemy_ordering(request)
        logger.info(str(('order_clauses', [str(c) for c in order_clauses])))
        if order_clauses:
            _alias = Alias(stmt)
            stmt = select([text('*')]).select_from(_alias)
            stmt = stmt.order_by(*order_clauses)
        
        logger.info(str(('filter_expression', str(filter_expression))))
        if filter_expression is not None:
            if not order_clauses:
                _alias = Alias(stmt)
                stmt = select([text('*')]).select_from(_alias)
            stmt = stmt.where(filter_expression)

        # need the count
        # TODO: select (and join) only the columns that are being shown
        if filter_fields is not None:
            count_fields = [field for field in schema['fields'].values() 
                if field['key'] in filter_fields ]
            logger.info(str(('filter_fields', filter_fields, 'count_fields', count_fields)))
            count_columns = self.build_sqlalchemy_columns(count_fields, base_query_tables)
            logger.info(str(('count_columns', count_columns)))
            if count_columns:
                count_stmt = select(count_columns.values()).\
                    select_from(j).\
                    where(_well.c.library_id == library.library_id) 
                _alias = Alias(count_stmt)
                count_stmt = select([text('*')]).select_from(_alias)
                count_stmt = count_stmt.where(filter_expression)
                logger.info(str(('count_stmt',str(count_stmt))))
            else:
                logger.error('no count columns')
        else:
            count_stmt = select(columns).\
                select_from(j).\
                where(_well.c.library_id == library.library_id) 
        count_stmt = select([func.count()]).select_from(count_stmt.alias())
        
        return self.stream_response_from_cursor(
            request, stmt, count_stmt, filename, field_hash=field_hash  )
 
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
    
    def get_reagent_resource(self, library):
        # FIXME: we should store the "type" on the entity
        
        if library.screen_type == 'rnai':
            return self.get_sr_resource()
        else:
            if library.library_type == 'natural_products':
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
        sub_resource = self.get_reagent_resource(library)
        query = sub_resource.get_object_list(request)
        logger.info(str(('==== query', query.query.sql_with_params())))
        
        ## also add in the "supertype" fields:
        query.select_related('well')
    
        if library_short_name:
            query = query.filter(well__library=library)
#             logger.debug(str(('get reagent/well list', library_short_name, len(query))))
        return query

    def full_dehydrate(self, bundle, for_list=False):
#         bundle = super(ReagentResource, self).full_dehydrate(bundle)
        
        well_bundle = self.build_bundle(bundle.obj.well, request=bundle.request)
        well_bundle = self.get_well_resource().full_dehydrate(well_bundle)
        bundle.data.update(well_bundle.data)
        
        library = bundle.obj.well.library
        sub_resource = self.get_reagent_resource(library)
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
        schema = deepcopy(super(ReagentResource,self).build_schema())
        
        if library:
            sub_data = self.get_reagent_resource(library).build_schema()
            
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
#         serializer = CursorSerializer()

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
    
    def get_reagent_resource(self, library):
        # FIXME: we should store the "type" on the entity
        
        if library.screen_type == 'rnai':
            return self.get_sr_resource()
        else:
            if library.library_type == 'natural_products':
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
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<well_id>((?=(schema))__|(?!(schema))[^/]+))%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]
                
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
            sub_data = self.get_reagent_resource(library).build_schema()
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

    def get_list(self, request, **kwargs):
    
        if 'library_short_name' in kwargs:
            library = Library.objects.get(short_name=kwargs['library_short_name'])
            return self.get_full_reagent_resource().get_list(request, **kwargs)
        else:
            raise NotImplementedError('must provide a library_short_name parameter')
    
#     def obj_get_list(self, bundle, **kwargs):    
# 
#         if 'library_short_name' in kwargs:
#             library_short_name = kwargs.pop('library_short_name')
#             try:
#                 library = Library.objects.get(short_name=library_short_name)
#                 reagent_resource = self.get_reagent_resource(library)
#                 
#                 reagent_query = reagent_resource.obj_get_list(bundle, **kwargs)
#                 
#                 
#                 sql = ('select w.*,r.* '
#                     'from well w left outer join ({reagent_query}) r on(r.well_id=w.well_id) '
#                     'where w.library=%s')
#                 sql = sql.format(reagent_query=reagent_query.query.sql_with_params())
#                 logger.info(str(('===sql', sql)))
#                 
#                 cursor = connection.cursor()
#                 cursor.execute(sql, library.library_id)
#                 logger.info(str(('===sql2', sql)))
#                 
#                 
#                 class CursorIterator:
#                     def __init__(self, cursor):
#                         self.cursor = cursor
#                 
#                     def __iter__(self):
#                         return self
#                 
#                     def next(self): # Python 3: def __next__(self)
#                         obj = cursor.fetchone()
#                         if obj:
#                             bundle = self.build_bundle(obj=obj,request=bundle.request)
#                             bundle = reagent_resource.full_dehydrate(bundle)
#                             bundle = self.full_dehydrate(bundle)
#                             return bundle
#                         else:
#                             raise StopIteration
#                 
#                 return CursorIterator(cursor)
#                 
#             except Exception, e:
#                 exc_type, exc_obj, exc_tb = sys.exc_info()
#                 fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
#                 msg = str(e)
#                 logger.warn(str(('on obj_get_list',msg,
#                     self._meta.resource_name, msg, exc_type, fname, exc_tb.tb_lineno)))
#                 raise e 
# #                 raise Http404(str(('err', e)))
#                 
#     def apply_sorting(self, obj_list, options=None):
#         # disabled, for now
#         return obj_list
    
        
    def get_object_list(self, request, library_short_name=None):
        ''' 
        Note: any extra kwargs are there because we are injecting them into the 
        global tastypie kwargs in one of the various "dispatch_" handlers assigned 
        through prepend_urls.  Here we can explicitly add them to the query. 
        
        '''
        query = super(WellResource, self).get_object_list(request);
        if library_short_name:
            query = query.filter(library__short_name=library_short_name)
            logger.debug(str(('get well list', library_short_name, len(query))))
            
        return query
                    

    def dehydrate(self, bundle):
        
#         library = bundle.obj.library
#         
#         # TODO: migrate to using "well.reagents"
#         # FIXME: need to create a migration script that will invalidate all of the
#         # reagent.well_id's for reagents other than the "latest released reagent"
#         reagent_resource = self.get_reagent_resource(library)
#         if bundle.obj.reagent_set.exists():
#             reagent = bundle.obj.reagent_set.all()[0]            
#             sub_bundle = reagent_resource.build_bundle(
#                 obj=reagent, request=bundle.request)
#             sub_bundle = reagent_resource.full_dehydrate(sub_bundle)
#             if 'resource_uri' in sub_bundle.data:
#                 del sub_bundle.data['resource_uri'] 
#             bundle.data.update(sub_bundle.data)
#         else:
#             sub_bundle = reagent_resource.build_bundle(
#                 obj=Reagent(), request=bundle.request)
#             sub_bundle = reagent_resource.full_dehydrate(sub_bundle)
#             if 'resource_uri' in sub_bundle.data:
#                 del sub_bundle.data['resource_uri'] 
#             bundle.data.update(sub_bundle.data)
        return bundle
    
    def dehydrate_library(self, bundle):
        
        return self.get_library_resource().get_resource_uri(
            bundle_or_obj=bundle, url_name='api_dispatch_list')


    def post_list(self, request, **kwargs):
        raise NotImplementedError("Post is not implemented for ReagentResource, use patch instead")
    
    def patch_list(self, request, **kwargs):
        # TODO: NOT TESTED
        return self.put_list(request, **kwargs)
    
    @transaction.atomic()
    def put_list(self, request, **kwargs):

        if 'library_short_name' not in kwargs:
            raise BadRequest('library_short_name is required')
        
        deserialized = self.deserialize(
            request, 
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
                    well_id = '%s:%s' %(plate_number, well_name)

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
        sub_resource = self.get_reagent_resource(library)
        
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

# FIXME: will replace this with the ApiLog resource?
class ActivityResource(ManagedModelResource):

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

        
    def __init__(self, **kwargs):
        super(ActivityResource,self).__init__(**kwargs)

    def dehydrate(self, bundle):
        # FIXME: hack for demo
        bundle.data['activity_type'] = 'library_contents_update' 
        return bundle
    
#     def prepend_urls(self):
#         # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
#         # allows us to match any word (any char except forward slash), 
#         # except "schema", and use it as the key value to search for.
#         # also note the double underscore "__" is because we also don't want to 
#         # match in the first clause.
#         # We don't want "schema" since that reserved word is used by tastypie 
#         # for the schema definition for the resource (used by the UI)
#         return [
#             url((r"^(?P<resource_name>%s)/(?P<plate_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
#                 )  % (self._meta.resource_name, trailing_slash()), 
#                 self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]

class LibraryResource(ManagedModelResource):
    
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
    #         self._meta.validation = ManagedValidation(scope=self.scope)

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
            # override the parent "base_urls" so that we don't need to worry about schema again
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
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/copy%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_librarycopyview'), 
                name="api_dispatch_librarycopyview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/librarycopies%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_librarycopiesview'), 
                name="api_dispatch_librarycopiesview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/plate%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_libraryplateview'), 
                name="api_dispatch_libraryplateview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>[\w\d_.\-\+: ]+)"
                 r"/librarycopyplates%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_librarycopyplatesview'), 
                name="api_dispatch_librarycopyplatesview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/well%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_wellview'), 
                name="api_dispatch_library_wellview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/reagent%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_reagentview'), 
                name="api_dispatch_library_reagentview"),
            
#             url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
#                  r"/reagent2%s$" ) % (self._meta.resource_name, trailing_slash()), 
#                 self.wrap_view('dispatch_library_reagentview2'), 
#                 name="api_dispatch_library_reagentview2"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/reagent/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_reagent_schema'), name="api_get_reagent_schema"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/well/schema%s$") 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_well_schema'), name="api_get_well_schema"),
                
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/copy/(?P<copy_name>[^/]+)"
                 r"/plate%s$") % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_librarycopyplateview'), 
                name="api_dispatch_library_copy_plateview"),

            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/version%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_libraryversionview'), 
                name="api_dispatch_libraryversionview"),
        ]    

#     def get_schema(self, request, **kwargs):
#         if not 'short_name' in kwargs:
#             raise Http404(unicode((
#                 'The well schema requires a library short name'
#                 ' in the URI, as in /library/[short_name]/well/schema/')))
#         short_name = kwargs['short_name']
#         try:
#             library = Library.objects.get(short_name=short_name)
#             return self.create_response(request, self.build_schema(library))
#         except Library.DoesNotExist, e:
#             raise Http404(unicode((
#                 'no library found for short_name', short_name)))

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
  
    def dispatch_librarycopyview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyResource().dispatch('list', request, **kwargs)    
 
    def dispatch_librarycopiesview(self, request, **kwargs):
        logger.info(str(('short_name',kwargs['short_name'])))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopiesResource().dispatch('list', request, **kwargs)    
 
    def dispatch_libraryplateview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)   
    
    def dispatch_librarycopyplatesview(self, request, **kwargs): 
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyPlatesResource().dispatch('list', request, **kwargs)   

    def dispatch_library_wellview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_well_resource().dispatch('list', request, **kwargs)    
                    
    def dispatch_library_reagentview(self, request, **kwargs):
        logger.info(str(('dispatch_library_reagentview ', kwargs)))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_reagent_resource().dispatch('list', request, **kwargs)    
                    
#     def dispatch_library_reagentview2(self, request, **kwargs):
#         logger.info(str(('dispatch_library_reagentview ', kwargs)))
#         kwargs['library_short_name'] = kwargs.pop('short_name')
#         return ReagentResource2().dispatch('list', request, **kwargs)    
                    
    def dispatch_library_copyplateview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)    
        
    def dispatch_libraryversionview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryContentsVersionResource().dispatch('list', request, **kwargs)    
        
    def dehydrate(self, bundle):
        # get the api comments
        
        
        # FIXME: just poc: gets_all_ apilog comments, at this time
        # TODO: how to limit the number of comments?
        # FIXME: how to bypass hydrating comments when in the LoggingMixin on update?
        comments = self.get_apilog_resource().obj_get_list(
            bundle, ref_resource_name='library', key=bundle.obj.short_name)
        comment_list = []
        if len(comments) > 0:
            for comment in comments[:10]:
                # manually build the comment bundle, 
                # because the apilog.dehydrate_child_logs is non-performant
                
                comment_bundle = {
                    'username': comment.username,
                    'date_time': comment.date_time,
                    'comment': comment.comment,
                    'ref_resource_name': comment.ref_resource_name,
                    'key': comment.key,
                    }
#                 comment_bundle = self.get_apilog_resource().build_bundle(obj=comment)
#                 comment_bundle = self.get_apilog_resource().full_dehydrate(comment_bundle);
#                 comment_list.append(comment_bundle.data);
                comment_list.append(comment_bundle)
        bundle.data['comments'] = comment_list;
        return bundle
    
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

#     def obj_update(self, bundle, **kwargs):
#         bundle = super(LibraryResource, self).object_update(bundle, **kwargs)
#         # clear the cached schema because plate range have updated
#         cache.delete(self._meta.resource_name + ':schema')
# 
#         return bundle;

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
    
#     def hydrate(self,bundle):        
#         library = Library.objects.get(short_name=bundle.data['library_short_name'])
#         
#         from django.db.models import Max
#         result = LibraryContentsVersion.objects.all()\
#             .filter(library=library).aggregate(Max('version_number'))
#         version_number = result['version_number__max'] or 0
#         bundle.obj.library = library;
#         bundle.obj.version_number = version_number + 1
#         
#         return bundle


