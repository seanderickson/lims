import logging
from collections import defaultdict

from django.conf.urls import url
from django.forms.models import model_to_dict
from django.http import Http404
from django.db import connection

from tastypie.utils import timezone
from tastypie.exceptions import BadRequest, ImmediateHttpResponse
from tastypie.utils.urls import trailing_slash
from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication,\
    MultiAuthentication
from tastypie.constants import ALL_WITH_RELATIONS
from tastypie import fields

from db.models import ScreensaverUser,Screen, LabHead, LabAffiliation, \
    ScreeningRoomUser, ScreenResult, DataColumn, Library, Plate, Copy,\
    PlateLocation, Reagent, Well, LibraryContentsVersion, Activity,\
    AdministrativeActivity, SmallMoleculeReagent, SilencingReagent, GeneSymbol,\
    NaturalProductReagent, Molfile
from db.support import lims_utils
from reports.serializers import CursorSerializer, LimsSerializer
from reports.models import MetaHash, Vocabularies, ApiLog
from reports.api import ManagedModelResource, ManagedResource, ApiLogResource,\
    UserGroupAuthorization
from django.contrib.auth.models import User
from tastypie.validation import Validation
import math

        
logger = logging.getLogger(__name__)

    
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
            orderings = map(
                lambda x: ( 
                    x[1:] + ' DESC NULLS LAST' 
                        if x[0]=='-' else x + ' ASC NULLS FIRST' ), order_by ) 
            sql = ( 'SELECT * FROM ( ' + sql + ') as order_inner ORDER BY ' +
                    ', '.join(orderings) )
                     
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


class NaturalProductReagentResource(ManagedModelResource):
    
    library_short_name = fields.CharField('library__short_name',  null=True)

    class Meta:
        queryset = NaturalProductReagent.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'naturalproductreagent'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(NaturalProductReagentResource,self).__init__(**kwargs)


    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            # Note this entry will direct to put_list
            url((r"^(?P<resource_name>%s)/(?P<library_short_name>[^/]+)%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_list'), name="api_dispatch_list"),
            # Use this entry for read only  
            # TODO: does this work?  
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<reagent__library__short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/(?P<reagent__library_contents_version__version_number>[^/]+)"
                 r"/(?P<reagent__well_id>[^/]+)%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]    


class SilencingReagentResource(ManagedModelResource):
    
    library_short_name = fields.CharField('library__short_name',  null=True)

    class Meta:
        queryset = SilencingReagent.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'silencingreagent'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(SilencingReagentResource,self).__init__(**kwargs)

    def dehydrate(self, bundle):
        gene = bundle.obj.reagentfacilitygenes_set.all()[0]
        if gene:
            gene = gene.gene
            bundle.data['facility_gene_id'] = gene.entrezgene_id
            bundle.data['facility_gene_name'] = gene.gene_name
            # NOTE: Requires db/migrations/manual/0013 to convert gene_symbol
            bundle.data['facility_gene_symbols'] = [
                x.entrezgene_symbol for x in 
                    GeneSymbol.objects.all()
                              .filter(gene=gene).order_by('ordinal')]
        gene = bundle.obj.reagentvendorgenes_set.all()[0]
        if gene:
            gene = gene.gene
            bundle.data['vendor_gene_id'] = gene.entrezgene_id
            bundle.data['vendor_gene_name'] = gene.gene_name
            # NOTE: Requires db/migrations/manual/0013 to convert gene_symbol
            bundle.data['vendor_gene_symbols'] = [
                x.entrezgene_symbol for x in 
                    GeneSymbol.objects.all()
                              .filter(gene=gene).order_by('ordinal')]
            
        return bundle

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url((r"^(?P<resource_name>%s)/(?P<reagent__substance_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]


class SmallMoleculeReagentResource(ManagedModelResource):
    
    reagent_substance_id = fields.CharField('reagent__substance_id')
#     reagent = fields.ToOneField('db.api.ReagentResource', 'reagent')
    moldata = fields.CharField(null=True, blank=True)
        
    class Meta:
        queryset = SmallMoleculeReagent.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'smallmoleculereagent'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
    def __init__(self, **kwargs):
        super(SmallMoleculeReagentResource,self).__init__(**kwargs)

    def dehydrate_moldata(self, bundle):
        if hasattr(bundle.obj, 'molfile'):
            logger.info(str(('===dehydrating moldata', bundle.obj.molfile.molfile )))
            return bundle.obj.molfile.molfile
        else:
            logger.info(str(('no moldata for ', bundle.obj)))
    
    
    def hydrate(self, bundle):
        '''
        obj_create{ full_hydrate{ hydrate, hydrate_foo }, save{ is_valid, save_related, save, save_m2m }}
        
        '''
        bundle.data['date_created'] = timezone.now()
        
        logger.info(str(('===hydrate smr', bundle.data)))

        # NOTE: see "put_list" for performant method
        # TODO: drive the key value from the schema
        fieldname = 'reagent_substance_id'
        if fieldname not in bundle.data:
            bundle.errors[self._meta.resource_name] = {fieldname: 'required'}
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))
        reagent = Reagent.objects.get(substance_id = bundle.data['reagent_substance_id'])
        bundle.obj.reagent = reagent
        
        return bundle

    def obj_create(self, bundle, **kwargs):
        '''
        obj_create{ full_hydrate{ hydrate, hydrate_foo }, save{ is_valid, save_related, save, save_m2m }}
        
        '''
        bundle = super(SmallMoleculeReagentResource, self).obj_create(bundle, **kwargs)
        logger.info(str(('=== created smr', bundle.obj)))

        # all related objs must be created after creating the reagent
        if 'moldata' in bundle.data:
            molfile = Molfile()
            molfile.molfile = bundle.data['moldata']
            molfile.reagent = bundle.obj
            molfile.ordinal = 0
            molfile.save()
            logger.info(str(('=== created moldata for reagent', bundle.obj)))
        else:
            logger.warn(str(('=======no moldata found============',bundle.data)))
            
    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            # Note this entry will direct to put_list
            url((r"^(?P<resource_name>%s)/(?P<library_short_name>[^/]+)%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_smr_libraryview'), name="api_dispatch_smr_libraryview"),
            # Use this entry for read only    
            url((r"^(?P<resource_name>%s)"
                 r"/(?P<reagent__library__short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/(?P<library_contents_version__version_number>[^/]+)"
                 r"/(?P<reagent__well_id>[^/]+)%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]    

    def dispatch_smr_libraryview(self, request, **kwargs):
#         kwargs['library_short_name'] = kwargs.pop('short_name')
        return super(SmallMoleculeReagentResource, self).dispatch('list', request, **kwargs)    

    def get_object_list(self, request, library_short_name=None):
        ''' 
        Called immediately before filtering, actually grabs the (ModelResource) 
        query - 
        
        Override this and apply_filters, so that we can control the extra 
        column "is_for_group":
        This extra column is present when navigating to permissions from a 
        usergroup; see prepend_urls.
        '''
        query = super(SmallMoleculeReagentResource, self).get_object_list(request);
        logger.info(str(('get_obj_list', len(query))))
        if library_short_name:
            query = query.filter(reagent__library__short_name=library_short_name)
        return query
    

class ReagentResource(ManagedModelResource):

#     well = fields.ToOneField('db.api.WellResource', 
#                                         attribute='well', 
#                                         full=True, full_detail=True, full_list=True,
#                                         null=True)
    well = fields.CharField(null=True)
    
    class Meta:
        queryset = Reagent.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'reagent'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

        
    def __init__(self, **kwargs):
        self.sr_resource = None
        self.smr_resource = None
        self.npr_resource = None
        self.well_resource = None
        super(ReagentResource,self).__init__(**kwargs)


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
            npr_resource = NaturalProductReagentResource()
        return self.npr_resource
    
    def get_well_resource(self):
        if not self.well_resource:
            self.well_resource = WellResource()
        return self.well_resource
    
    def get_sub_resource(self, library):
        # FIXME: we should store the "type" on the entity
        if library.screen_type == 'rnai':
            return self.get_sr_resource()
        else:
            if library.library_type == 'natural_products':
                return self.get_npr_resource()
            else:
                return self.get_smr_resource()

    def get_sub_reagent(self, reagent, library):
        # FIXME: we should store the "type" on the entity
        if library.screen_type == 'rnai':
            if hasattr(reagent, 'silencingreagent'):
                return reagent.silencingreagent
        else:
            if library.library_type == 'natural_products':
                if hasattr(reagent, 'naturalproductreagent'):
                    return  reagent.naturalproductreagent
            else:
                if hasattr(reagent, 'smallmoleculereagent'):
                    return reagent.smallmoleculereagent
        return None
                    
    def get_schema(self, request, **kwargs):
        # FIXME: we should store the "type" on the entity
        if not 'library_short_name' in kwargs:
            return self.create_response(request, self.build_schema())
        
        library_short_name = kwargs.pop('library_short_name')
        try:
            library = Library.objects.get(short_name=library_short_name)
            return self.create_response(request, self.build_schema(library))
            
        except Library.DoesNotExist, e:
            raise Http404(unicode((
                'no library found for short_name', library_short_name)))

    def build_schema(self, library=None):
        logger.debug(str(('==========build schema for library reagent', library)))
        data = super(ReagentResource,self).build_schema()
        
        if library:
            sub_data = self.get_sub_resource(library).build_schema()
            data['fields'].update(sub_data['fields'])

        return data
    
    def dehydrate(self, bundle):

        reagent = bundle.obj
        if not reagent.well:
            return bundle
        library = reagent.well.library
        if not library:
            return bundle

        bundle.data['well_id'] = reagent.well.well_id
        logger.info(str(('===well_id', reagent.well.well_id)))
        
        sub_reagent = self.get_sub_reagent(reagent, library)
        if sub_reagent:
            sub_resource = self.get_sub_resource(library)
            sub_bundle = sub_resource.build_bundle(
                obj=sub_reagent, request=bundle.request)
            sub_bundle = sub_resource.full_dehydrate(sub_bundle)
            bundle.data.update(sub_bundle.data)
        return bundle

    def prepend_urls(self):
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            url((r"^(?P<resource_name>%s)/(?P<substance_id>((?=(schema))__|(?!(schema))[^/]+))%s$"
                )  % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]
    
    
    def put_list(self, request, **kwargs):

        if 'library_short_name' not in kwargs:
            raise BadRequest('library_short_name is required')
        
        basic_bundle = self.build_bundle(request=request)
        logger.info(str(('b', basic_bundle)))
 
        # Look for library in kwargs, create a version        
        library = Library.objects.get(short_name=kwargs['library_short_name'])
                          
        well_resource = self.get_well_resource()  
        deserialized = self.deserialize(
            request, request.body, 
            format=request.META.get('CONTENT_TYPE', 'application/json'))
        
        # Cache all the wells on the library for use with this process 
        wellMap = dict( (well.well_id, well) for well in library.well_set.all())
        if len(wellMap)==0:
            raise BadRequest(str(('library has not been created, no wells', library)))
        
        i=0
        for well_data in deserialized[self._meta.collection_name]:
            
            well_data['library_short_name']=kwargs['library_short_name']
            
            logger.info(str(('well_data', well_data)))
            well_id = well_data.get('well_id', None)
            if not well_id:
                well_name = well_data.get('well_name', None)
                plate_number = well_data.get('plate_number',None)
                if well_name and plate_number:                
                    well_id = '%s:%s' %(plate_number, well_name)

            if not well_id:
                bundle.errors[self._meta.resource_name] = {'well_id': 'required'}
                raise ImmediateHttpResponse(
                    response=self.error_response(bundle.request, bundle.errors))
            
            well = wellMap.get(well_id, None)
            if not well:
                bundle.errors[self._meta.resource_name] = {
                    'well_id': str(('well not found', well_id))}
                raise ImmediateHttpResponse(
                    response=self.error_response(bundle.request, bundle.errors))
                
            well_bundle = well_resource.build_bundle(data=well_data, request=request);
            well_bundle.obj = well
            well_bundle = self.full_hydrate(well_bundle)
            well_resource.save(well_bundle)
            
            logger.info(str(('updated well', well_bundle.obj)))

            # lookup/create the reagent
            
            if not 'resource_uri' in well_data:
                
                reagent_bundle = self.build_bundle(data=well_data, request=request)
                reagent_bundle = self.full_hydrate(reagent_bundle)
                reagent_bundle.obj.well = well
                self.save(reagent_bundle)
                
                # todo: figure out if there is a SMR/etc. to hydrate/create
                logger.info(str(('created reagent', reagent_bundle.obj)))

                sub_resource = self.get_sub_resource(library)
                if library.screen_type == 'rnai':
                    sub_bundle = sub_resource.build_bundle(data=well_data, request=request)
                    sub_bundle.data['reagent_substance_id'] = reagent_bundle.obj.substance_id
                    sub_bundle = sub_resource.full_hydrate(sub_bundle)
                    sub_resource.save(sub_bundle)
                else:
                    if library.library_type == 'natural_products':
                        sub_bundle = sub_resource.build_bundle(data=well_data, request=request)
                        sub_bundle.data['reagent_substance_id'] = reagent_bundle.obj.substance_id
                        sub_bundle = sub_resource.full_hydrate(sub_bundle)
                        sub_resource.save(sub_bundle)
                    else:
                        sub_bundle = sub_resource.build_bundle(data=well_data, request=request)
                        sub_bundle.data['reagent_substance_id'] = reagent_bundle.obj.substance_id
                        sub_bundle = sub_resource.full_hydrate(sub_bundle)
                        sub_resource.save(sub_bundle)
                        if 'moldata' in sub_bundle.data:
                            if hasattr(sub_bundle.obj, 'molfile'):
                                logger.info(str(('===hydrating extant moldata', bundle.obj.molfile.molfile )))
                                bundle.obj.molfile.molfile = sub_bundle.data['moldata']
                            else:
                                molfile = Molfile()
                            molfile.molfile = sub_bundle.data['moldata']
                            molfile.reagent = sub_bundle.obj
                            molfile.ordinal = 0
                            molfile.save()
                            logger.info(str(('=== created moldata for reagent', sub_bundle.obj)))
                
            else:
                # FIXME: untested, not complete
                # lookup and update the reagent
                self.obj_update(reagent_bundle)
                logger.info(str(('updated reagent', reagent_bundle.obj)))
            i = i+1
        
        logger.info(str(('put reagent', i)))

    def obj_create(self, bundle, **kwargs):
        '''
        obj_create{ full_hydrate{ hydrate, hydrate_foo }, save{ is_valid, save_related, save, save_m2m }}
        
        '''
        # get the library: this is available because the reagent resource is only
        # available through a library resource
        if 'library_short_name' not in kwargs:
            raise BadRequest('library_short_name is required')
        # Look for library in kwargs, create a version        
        library = Library.objects.get(short_name=kwargs['library_short_name'])

        logger.info(str(('===creating reagent', bundle.data)))
        
        bundle = super(ReagentResource, self).obj_create(bundle, **kwargs)
        
        sub_resource = self.get_sub_resource(library)
        sub_bundle = sub_resource.build_bundle(
            data=bundle.data, request=bundle.request)
        sub_bundle.data['reagent_substance_id'] = bundle.obj.substance_id
        
        sub_bundle = sub_resource.obj_create(sub_bundle)

#         bundle.data.update(sub_bundle.data)        
        return bundle

    def hydrate(self, bundle):
        logger.info(str(('==== hydrating the reagent for bundle ', bundle.data)))
        
        if ( not hasattr(bundle.obj,'well') or
                not getattr(bundle.obj,'well')):
            well_id = bundle.data.get('well_id', None)
            if not well_id:
                well_name = bundle.data.get('well_name', None)
                plate_number = bundle.data.get('plate_number',None)
                if well_name and plate_number:                
                    well_id = '%s:%s' %(plate_number, well_name)
            if not well_id:
                bundle.errors[self._meta.resource_name] = {'well_id': 'required'}
                raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))

            well = Well.objects.get(well_id=well_id)
        else:
            well = bundle.obj.well

        bundle.obj.well = well
        
        return bundle



class WellResource(ManagedModelResource):
#     # TODO: no-no to cache this way; need to lazy grab it like get_apilog_resource()
#     reagent_resource = ReagentResource()

    library_short_name = fields.CharField('library__short_name',  null=True)
#     reagent = fields.ToOneField('db.api.ReagentResource', 
#                                         attribute='latest_released_reagent', 
#                                         full=True, full_detail=True, full_list=True,
#                                         null=True)
    ###
    ##FIXME 20140905
    ## unfortunately, tastypie instantiates the related resource with _every_ 
    ## field dehydrate!
    ## so we are going to forgoe related fields for read purposes
    #library = fields.ToOneField('db.api.LibraryResource', attribute='library')    
    library = fields.CharField(null=True)
    # TODO:
    # status_date
    
    # plate_location = 
    
    class Meta:
        queryset = Well.objects.all() #.order_by('facility_id')
        authentication = MultiAuthentication(BasicAuthentication(), 
                                             SessionAuthentication())
        authorization= UserGroupAuthorization()
        resource_name = 'well'
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        
        # this makes Backbone/JQuery happy because it likes to JSON.parse the returned data
        always_return_data = True 

        
    def __init__(self, **kwargs):
        self.reagent_resource = None
        self.library_resource = None
        super(WellResource,self).__init__(**kwargs)

    def get_reagent_resource(self):
        if not self.reagent_resource:
            self.reagent_resource = ReagentResource()
        return self.reagent_resource

    def get_library_resource(self):
        if not self.library_resource:
            self.library_resource = LibraryResource()
        return self.library_resource


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
                 r"/(?P<well_id>((?=(schema))__|(?!(schema))[^/]+))%s$")  
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),]
                
    def get_object_list(self, request, library_short_name=None):
        ''' 
        Note: any extra kwargs are there because we are injecting them into the 
        global tastypie kwargs in one of the various "dispatch_" handlers assigned 
        through prepend_urls.  Here we can explicitly add them to the query. 
        
        '''
        query = super(WellResource, self).get_object_list(request);
        
        if library_short_name:
            query = query.filter(library__short_name=library_short_name)
        logger.info(str(('get well list', library_short_name, len(query))))
        return query
                    

    def dehydrate(self, bundle):
        
        library = bundle.obj.library
        
        # FIXME: migrate to using "well.reagent", and later, "well.reagents"
        if bundle.obj.reagent_set.exists():
            reagent = bundle.obj.reagent_set.all()[0]            
            logger.info(str(('=== found well reagent', bundle.obj, reagent)))
            sub_bundle = self.get_reagent_resource().build_bundle(
                obj=reagent, request=bundle.request)
            sub_bundle = self.get_reagent_resource().full_dehydrate(sub_bundle)
            logger.info(str(('=== dehydrated reagent', sub_bundle.data)))
            bundle.data.update(sub_bundle.data)
        else:
            #             logger.debug(str(('==== no well reagent', bundle.obj)))
            pass
        return bundle
    
    def dehydrate_library(self, bundle):
        
        return self.get_library_resource().get_resource_uri(
            bundle_or_obj=bundle, url_name='api_dispatch_list')


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
        logger.debug(str(('==========build schema for library wells', library)))
        data = super(WellResource,self).build_schema()
        
        if library:
            sub_data = self.get_reagent_resource().build_schema(library=library)
            data['fields'].update(sub_data['fields'])
#             if library.screen_type == 'rnai':
#                 sub_data = SilencingReagentResource().build_schema()
#                 data['fields'].update(sub_data['fields'])
#             else:
#                 if library.library_type == 'natural_products':
#                     sub_data = NaturalProductReagentResource().build_schema()
#                     data['fields'].update(sub_data['fields'])
#                 else:
#                     sub_data = SmallMoleculeReagentResource().build_schema()
#                     data['fields'].update(sub_data['fields'])

        temp = [ x.title.lower() 
            for x in Vocabularies.objects.all().filter(scope='library.well_type')]
        data['extraSelectorOptions'] = { 
            'label': 'Type', 'searchColumn': 'library_well_type', 'options': temp }

        return data
    
    def obj_update(self,bundle,**kwargs):
        '''
        patch_list (if resource_uri present):
        put_detail
        patch_detail
        obj_update{
            if not bundle.obj:
                lookup_kwargs_with_identifiers
                obj.get
            obj <- full_hydrate { hydrate, hydrate_foo }
            save { is_valid, save_related, save, save_m2m }
        '''
        bundle = super(WellResource, self).obj_update(bundle, **kwargs)
        
        logger.info(str(('created well', bundle.obj)))
        return bundle
    
    def obj_create(self, bundle, **kwargs):
        '''
        patch_list:
        put_list:
        post_list:
        put_detail:
        
        obj_create{ 
            full_hydrate{ hydrate, hydrate_foo }, 
            save{ is_valid, save_related, save, save_m2m }}
        
        '''
        
        logger.info(str(('kwargs', kwargs)))
        if 'library_short_name' in kwargs:
            bundle.data['library_short_name'] = kwargs['library_short_name']
                    
        bundle = super(WellResource, self).obj_create(bundle, **kwargs)
        
        logger.info(str(('created well', bundle.obj)))
        return bundle
        
    def hydrate(self, bundle):
        """
        obj_create{ full_hydrate{ hydrate, hydrate_foo, field.hydrate }, 
                    save{ is_valid, save_related, save, save_m2m }}
         
        A hook to allow an initial manipulation of data before all methods/fields
        have built out the hydrated data.

        Useful if you need to access more than one hydrated field or want
        to annotate on additional data.

        Must return the modified bundle.
        """
        bundle.data['date_created'] = timezone.now()
        
#         bundle.data['version'] = 1
        logger.info(str(('===creating well', bundle.data)))

        # MIGRATE "OLD" SS1 fields to new fields
        from db.support.data_converter import convert_well_data
        bundle.data = convert_well_data(bundle.data)
        logger.info(str(('===creating well', bundle.data)))

        if 'library_short_name' not in bundle.data:
            bundle.errors[self._meta.resource_name] = {'library_short_name': 'required'}
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))

        # TODO: drive the key value from the schema
        # TODO: we may do this manually, for performance reasons, particularly
        # if we can signal that we are loading *many* wells for a particular library
        bundle.data['library'] = 'library/' + bundle.data['library_short_name']

        # Nest the whole dict to reprocess on the reagent
        # NOTE: this method uses TP to create/update the related reagent
        # TODO: need to test if contains a reagent:
        #         if 'substance_id' in bundle[data].keys():
        bundle.data['reagent'] = bundle.data 
                
        # other manual construction, grab the plate_number, well_name to make the well_id
        plate_number = bundle.data.get('plate_number', None)
        if not plate_number:
            bundle.errors[self._meta.resource_name] = {'plate_number': 'required'}
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))
        well_name = bundle.data.get('well_name', None)
        if not well_name:
            bundle.errors[self._meta.resource_name] = {'well_name': 'required'}
            raise ImmediateHttpResponse(response=self.error_response(bundle.request, bundle.errors))
        bundle.data['well_id'] = '%s:%s' % (plate_number, well_name)
        
        return bundle


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
        
        #FIXME apilog not inited?
#         self.apiLog = ApiLogResource()
        # TODO: 20140905 - initializing here causes the WellResource to not get all fields, 
        # not initiated the right way
#         ------self.well_resource = WellResource()
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
        ## TODO: try to fix the "schema" issues here by putting the schema match first,
#         i.e.
#             url(r"^(?P<resource_name>%s)/schema%s$" % (self._meta.resource_name, trailing_slash()), self.wrap_view('get_schema'), name="api_get_schema"),
        
        # NOTE: this match "((?=(schema))__|(?!(schema))[^/]+)" 
        # allows us to match any word (any char except forward slash), 
        # except "schema", and use it as the key value to search for.
        # also note the double underscore "__" is because we also don't want to 
        # match in the first clause.
        # We don't want "schema" since that reserved word is used by tastypie 
        # for the schema definition for the resource (used by the UI)
        return [
            
            url(r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/copy%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_librarycopyview'), 
                name="api_dispatch_librarycopyview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/plate%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_libraryplateview'), 
                name="api_dispatch_libraryplateview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/well%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_wellview'), 
                name="api_dispatch_library_wellview"),
            
            url((r"^(?P<resource_name>%s)/(?P<short_name>((?=(schema))__|(?!(schema))[^/]+))"
                 r"/reagent%s$" ) % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_library_reagentview'), 
                name="api_dispatch_library_reagentview"),
            
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

    def get_well_schema(self, request, **kwargs):
        if not 'short_name' in kwargs:
            raise Http404(unicode((
                'The well schema requires a library short name'
                ' in the URI, as in /library/[short_name]/well/schema/')))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return WellResource().get_schema(request, **kwargs)    
 
    def get_reagent_schema(self, request, **kwargs):
        if not 'short_name' in kwargs:
            raise Http404(unicode((
                'The well schema requires a library short name'
                ' in the URI, as in /library/[short_name]/well/schema/')))
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return ReagentResource().get_schema(request, **kwargs)    
 
    def dispatch_librarycopyview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyResource().dispatch('list', request, **kwargs)    
 
    def dispatch_libraryplateview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return LibraryCopyPlateResource().dispatch('list', request, **kwargs)    

    def dispatch_library_wellview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_well_resource().dispatch('list', request, **kwargs)    
                    
    def dispatch_library_reagentview(self, request, **kwargs):
        kwargs['library_short_name'] = kwargs.pop('short_name')
        return self.get_reagent_resource().dispatch('list', request, **kwargs)    
                    
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
            bundle, ref_resource_name='library', key=bundle.obj.short_name,
            comment__isnull=False)
        comment_list = []
        if len(comments) > 0:
            for comment in comments[:10]:
                comment_bundle = self.get_apilog_resource().build_bundle(obj=comment)
                comment_bundle = self.get_apilog_resource().full_dehydrate(comment_bundle);
                comment_list.append(comment_bundle.data);
        bundle.data['comments'] = comment_list;
        return bundle
    
    def build_schema(self):
        schema = super(LibraryResource,self).build_schema()
        temp = [ x.library_type for x in self.Meta.queryset.distinct('library_type')]
        schema['extraSelectorOptions'] = { 
            'label': 'Type', 'searchColumn': 'library_type', 'options': temp }
        return schema
    
    def obj_create(self, bundle, **kwargs):
        bundle.data['date_created'] = timezone.now()
        
        bundle.data['version'] = 1
        logger.info(str(('===creating library', bundle.data)))

        bundle = super(LibraryResource, self).obj_create(bundle, **kwargs)
        
        # now create the wells
        
        library = bundle.obj
        logger.info(str(('created library', library, library.start_plate, type(library.start_plate))))
        plate_size = int(library.plate_size)
#         rows = lims_utils.get_rows(plate_size)
#         cols = lims_utils.get_cols(plate_size)
        for plate in range(int(library.start_plate), int(library.end_plate)+1):
            for index in range(0,plate_size):
                well = Well()
                well.well_name = lims_utils.well_name_from_index(index,plate_size)
                well.well_id = lims_utils.well_id(plate,well.well_name)
                well.library = library
                well.plate_number = plate
                well.library_well_type = 'empty'
                well.save()
        
        return bundle
        

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
