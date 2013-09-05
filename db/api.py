from django.conf.urls import url
from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
from tastypie import fields, utils

from db.models import ScreensaverUser,Screen, LabHead, LabAffiliation, ScreeningRoomUser,\
    ScreenStatusItem, ScreenResult, DataColumn
from lims.api import CursorSerializer, LimsSerializer, PostgresSortingResource
from reports.models import MetaHash, Vocabularies
from django.http import Http404, HttpResponse

import time
import re

import logging
from reports.api import MetahashManagedResource, JsonAndDatabaseResource
from django.db import connection
from tastypie.resources import Resource
from tastypie.exceptions import BadRequest
        
logger = logging.getLogger(__name__)

    
class ScreensaverUserResource(MetahashManagedResource, PostgresSortingResource):
#    screens = fields.ToManyField('db.api.ScreenResource', 'screens', related_name='lab_head_id', blank=True, null=True)

    version = fields.IntegerField(attribute='version', null=True)
        
    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = ['digested_password']
        
    def __init__(self, **kwargs):
        self.scope = 'fields:screensaveruser'
        super(ScreensaverUserResource,self).__init__( **kwargs)
  
    def dehydrate(self, bundle):
        bundle = super(ScreensaverUserResource, self).dehydrate(bundle);
#        _time = time.time();
        bundle.data['screens'] = [x.facility_id for x in Screen.objects.filter(lab_head_id=bundle.obj.screensaver_user_id)]
#        logger.info(str(('dehydrate time', time.time()-_time )))
        return bundle        
      
    #    def full_dehydrate(self, bundle, for_list=False):
    #        _time = time.time();
    #        bundle = super(ScreensaverUserResource, self).full_dehydrate(bundle, for_list=for_list);
    #        logger.info(str(('full dehydrate time', time.time()-_time )))
    #        return bundle        
    
    def build_schema(self):
        schema = super(ScreensaverUserResource,self).build_schema()
        schema['idAttribute'] = ['screensaver_user_id']
        return schema
    
    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<screensaver_user_id>[\d]+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
    

class ScreeningRoomUserResource(PostgresSortingResource):
    screensaver_user = fields.ToOneField('db.api.ScreensaverUserResource', attribute='screensaver_user', full=True, full_detail=True, full_list=False)
    class Meta:
        queryset = ScreeningRoomUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
    
class LabAffiliationResource(PostgresSortingResource):   
    class Meta:
        queryset = LabAffiliation.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
    
from django.forms.models import model_to_dict

class LabHeadResource(PostgresSortingResource):

    screens = fields.ToManyField('db.api.ScreenResource', 'screens', related_name='lab_head', blank=True, null=True)

    lab_affiliation = fields.ToOneField('db.api.LabAffiliationResource', attribute='lab_affiliation',  full=True, null=True)
    
    # rather than walk the inheritance hierarchy, will flatten this hierarchy in the dehydrate method
    #    screening_room_user = fields.ToOneField('db.api.ScreeningRoomUserResource', attribute='screensaver_user',  full=True)
    
    id = fields.IntegerField(attribute='screensaver_user_id')
    
    class Meta:
        queryset = LabHead.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
#        resource_name = 'lab_head'
        
    def dehydrate(self, bundle):
        # flatten the inheritance hierarchy, rather than show nested "lab_head->screening_room_user->screensaver_user"
        bundle.data.update(model_to_dict(bundle.obj.screensaver_user))
        bundle.data.update(model_to_dict(bundle.obj.screensaver_user.screensaver_user))
        bundle.data['screens'] = [model_to_dict(x) for x in Screen.objects.filter(lab_head_id=bundle.obj.screensaver_user.screensaver_user_id)]
        
        return bundle        
    
class ScreenResultResource(MetahashManagedResource,Resource):

    class Meta:
        queryset = ScreenResult.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        resource_name = 'screenresult'
        
        ordering = []
        filtering = {}
        serializer = CursorSerializer()
        allowed_methods = ['get']

        object_class = dict
        
    def __init__(self, **kwargs):
        self.scope = 'fields:screenresult'
        super(ScreenResultResource,self).__init__(**kwargs)

#    def detail_uri_kwargs(self, bundle_or_obj):
#        kwargs = {}
#
#        if isinstance(bundle_or_obj, Bundle):
#            kwargs['pk'] = bundle_or_obj.obj.uuid
#        else:
#            kwargs['pk'] = bundle_or_obj.uuid
#
#        return kwargs


    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" allows us to match any word, except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie for the schema definition for the resource (used by the UI)
        return [
            url(r"^(?P<resource_name>%s)/(?P<facility_id>((?=(schema))__|(?!(schema))[\w\d_.-]+))/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<facility_id>((?=(schema))__|(?!(schema))[\w\d_.-]+))/schema/$" % self._meta.resource_name, self.wrap_view('get_schema'), name="api_get_schema"),
        ]    
        
    def get_object_list(self, request):
        logger.warn('Screen result listing not implemented')
        raise Http404(str(('Screen result listing not implemented',request.path)))
        
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
                if not isinstance(limit, (str, unicode)):  # TODO: why are request parameters being wrapped as lists?
                    # try it as a seq
                    limit = limit[0]
                try:
                    kwargs['limit'] = int(limit)
                except ValueError:
                    raise BadRequest("Invalid limit '%s' provided. Please provide a positive integer." % kwargs['limit'])

            if 'offset' in kwargs:
                offset = kwargs['offset']
                if not isinstance(offset, (str, unicode)):  # TODO: why are request parameters being wrapped as lists?
                    # try it as a seq
                    offset = offset[0]
                try:
                    kwargs['offset'] = int(offset)
                except ValueError:
                    raise BadRequest("Invalid offset '%s' provided. Please provide a positive integer." % kwargs['offset'])

        logger.info(unicode(('============= obj_get: kwargs', kwargs)))
        if 'facility_id' not in kwargs:
            raise Http404(unicode(('no facility id given',kwargs, request)))
            

        facility_id = kwargs['facility_id']

        try:
            screenresult = ScreenResult.objects.get(screen__facility_id=facility_id)
            logger.info(str(('screenresult resource', facility_id,screenresult.screen)))
            

            # TODO: Not sure what to return for get_obj_list ? since in the CursorSerializer, you can see that we have to either look 
            # at the passed arg as a cursor, or as a bundle with the cursor as the obj.  how does TP want this done?
            # see also: http://django-tastypie.readthedocs.org/en/latest/non_orm_data_sources.html        
            if screenresult:
                result = {}
                result['meta'] = kwargs.copy();
                result['meta']['total_count'] = self.get_total_count(screenresult)
                
                result['objects'] = self.get_screenresult_cursor(screenresult, **kwargs)
                return result;
            else:
                raise Http404(unicode(('no results for the screen: ', facility_id)))
        
        except Screen.DoesNotExist, e:
            logger.error(str(('no screen found for facility id', facility_id)))
            raise e
        
        
        return self.get_object_list(request)

    @staticmethod
    def create_query(screenresult):
        '''    
        select well_id, ...other_entity_columns,  
            (select numeric_value as col1 from result_value rv1 where data_column_id=4754 and rv1.well_id = w.well_id) as col1 
            from assay_well w join screen_result using (screen_result_id) join screen using (screen_id) where facility_id = '1003';
        '''
                
        sql = 'select well_id '
        for i,dc in enumerate(DataColumn.objects.filter(screen_result=screenresult)):
            alias = "dp_"+str(dc.data_column_id)
            columnName = "col_"+str(dc.data_column_id)
            column_to_select = None
            if(dc.data_type == 'Numeric'): #TODO: use controlled vocabulary
                column_to_select = 'numeric_value'
            else:
                column_to_select = 'value'

            sql +=  (",(SELECT " + column_to_select + " FROM result_value " + alias + 
                                " where " + alias + ".data_column_id="+str(dc.data_column_id) + " and " + alias + ".well_id=w.well_id) as " + columnName )
        sql += ' FROM assay_well w where w.screen_result_id=%s '
        return sql
    
    def get_total_count(self, screenresult, **kwargs):
        cursor = connection.cursor()
        
        sql = self.create_query(screenresult, **kwargs);
        sql = 'select count(*) from (' + sql + ') a'
        cursor.execute(sql, [screenresult.screen_result_id])
        return cursor.fetchone()[0]
        
    def get_screenresult_cursor(self, screenresult, limit=25, offset=0, order_by=[], **kwargs):
        logger.info(unicode(('---get_screenresult_cursor', 'limit, offset, order_by', limit, offset, order_by, 'kwargs', kwargs)))
        
        sql = self.create_query(screenresult)
        
        if len(order_by) > 0:
            orderings = map(lambda x: x[1:] + ' DESC NULLS LAST' if x[0]=='-' else x + ' ASC NULLS FIRST',order_by) # TODO: postgres only 
            logger.info(str(('==== orderings', orderings)))
            sql = 'SELECT * FROM ( ' + sql + ') as order_inner ORDER BY ' + ', '.join(orderings) 
                     
        sql += ' OFFSET ' + str(offset)
        sql += ' LIMIT ' + str(limit)
        
        logger.info(str(('sql',sql)))
        cursor = connection.cursor()
        cursor.execute(sql, [screenresult.screen_result_id])
        return cursor

    # override
    def get_schema(self, request, **kwargs):
#        bundle = self.build_bundle(request=request)
#        self.authorized_read_detail(self.get_object_list(bundle.request), bundle)
        facility_id = kwargs.pop('facility_id')
        if not facility_id:
            raise Http404(unicode(('The screenresult schema requires a screen facility ID in the URI, as in /screenresult/[facility_id]/schema/')))
        try:
            screenresult = ScreenResult.objects.get(screen__facility_id=facility_id)
            logger.info(str(('screenresult resource', facility_id,screenresult.screen)))
            

            # TODO: Not sure what to return for get_obj_list ? since in the CursorSerializer, you can see that we have to either look 
            # at the passed arg as a cursor, or as a bundle with the cursor as the obj.  how does TP want this done?
            # see also: http://django-tastypie.readthedocs.org/en/latest/non_orm_data_sources.html        
            if screenresult:
                return self.create_response(request, self.build_schema(screenresult))
            else:
                raise Http404(unicode(('no results for the screen: ', facility_id)))
        except Screen.DoesNotExist, e:
            raise Http404(unicode(('no screen found for facility id', facility_id)))
            

    def build_schema(self, screenresult):
        logger.info(unicode(('==========build schema for ', screenresult)))
        data = super(ScreenResultResource,self).build_schema()
        
        # now the datacolumn fields
        field_defaults = {
            'visibility': ['list','detail'],
            'ui_type': 'string',
            'type': 'string',
            'filtering': True,
            }
        for i,dc in enumerate(DataColumn.objects.filter(screen_result=screenresult)):
            alias = "dp_"+str(dc.data_column_id)
            columnName = "col_"+str(dc.data_column_id)
            _dict = field_defaults.copy()
            _dict.update(model_to_dict(dc))
            
            _dict['title'] = dc.name
            _dict['comment'] = dc.comments
            _dict['key'] = columnName
            _dict['ordinal'] += len(self.fields) + dc.ordinal # so that the value columns come last
#            if dc.data_type == 'Numeric':
#                _dict['ui_type'] = 'numeric'
            
            data['fields'][columnName] = _dict
        # TODO: get the data columns; convert column aliases to real names
#        logger.info(unicode(('schema ', data)))
        return data
    

class ScreenSummaryResource(JsonAndDatabaseResource):

    experimental_wells_loaded = fields.IntegerField('screenresult__experimental_well_count')    
    
    class Meta:
        queryset = Screen.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        resource_name = 'screensummary'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

    def __init__(self, **kwargs):
#        self.
        super(ScreenSummaryResource,self).__init__(scope = 'fields:screensummary', **kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" allows us to match any word, except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie for the schema definition for the resource (used by the UI)
        return [
            url(r"^(?P<resource_name>%s)/(?P<facility_id>((?=(schema))__|(?!(schema))[\w\d_.-]+))/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    

class ScreenResource(JsonAndDatabaseResource):

#    lab_head_full = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=True) #, full_list=False) #, blank=True, null=True)
    lab_head_link = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=False)
    lead_screener_link = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=False)
    
    lab_head_id = fields.IntegerField(attribute='lab_head_id');
    lead_screener_id = fields.IntegerField(attribute='lead_screener_id');
    
    class Meta:
        queryset = Screen.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        resource_name = 'screen'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()

        
    def __init__(self, **kwargs):
#        self.
        super(ScreenResource,self).__init__(scope = 'fields:screen', **kwargs)

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[\w\d_.-]+)" allows us to match any word, except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie for the schema definition for the resource (used by the UI)
        return [
            url(r"^(?P<resource_name>%s)/(?P<facility_id>((?=(schema))__|(?!(schema))[\w\d_.-]+))/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]    
                    
    def dehydrate(self, bundle):
#        _time = time.time();
        sru = bundle.obj.lead_screener.screensaver_user
        bundle.data['lead_screener'] =  sru.first_name + ' ' + sru.last_name
        lh = bundle.obj.lab_head.screensaver_user.screensaver_user
        bundle.data['lab_head'] =  lh.first_name + ' ' + lh.last_name
        # TODO: the status table does not utilize a primary key, thus it is incompatible with the standard manager
        #        status_item = ScreenStatusItem.objects.filter(screen=bundle.obj).order_by('status_date')[0]
        #        bundle.data['status'] = status_item.status
        #        bundle.data['status_date'] = status_item.status_date

#        logger.info(str(('dehydrate time', time.time()-_time )))
        bundle.data['has_screen_result'] = False
        try:
            bundle.data['has_screen_result'] = bundle.obj.screenresult != None
        except ScreenResult.DoesNotExist, e:
            logger.info(unicode(('no screenresult for ', bundle.obj)))
        return bundle

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
#                    options['order_by']= dir + 'lead_screener.screensaver_user.last_name';
                    logger.info(str(('remove ', temp)))
                    order_by.remove(temp)
                    extra_order_by.append(dir+'lead_screener__screensaver_user__last_name')
                    extra_order_by.append(dir+'lead_screener__screensaver_user__first_name')
                if field == 'lab_head':
                    logger.info(str(('remove ', temp)))
                    order_by.remove(temp)
                    extra_order_by.append(dir+'lab_head__screensaver_user__screensaver_user__last_name')
                    extra_order_by.append(dir+'lab_head__screensaver_user__screensaver_user__first_name')
                if field == 'has_screen_result':
                    logger.info(str(('remove ', temp)))
                    order_by.remove(temp)
                    obj_list = obj_list.extra({'screenresult_isnull': '(select sr.screen_id is null from screen_result sr where sr.screen_id = screen.screen_id) '})
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
    