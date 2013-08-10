from django.conf.urls import url
from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
from tastypie import fields, utils

from db.models import ScreensaverUser,Screen, LabHead, LabAffiliation, ScreeningRoomUser,\
    ScreenStatusItem
from lims.api import CSVSerializer, PostgresSortingResource
from reports.models import MetaHash, Vocabularies

import time

import logging
from reports.api import MetahashManagedResource
        
logger = logging.getLogger(__name__)

    
class ScreensaverUserResource(MetahashManagedResource, PostgresSortingResource):
#    screens = fields.ToManyField('db.api.ScreenResource', 'screens', related_name='lab_head_id', blank=True, null=True)

    version = fields.IntegerField(attribute='version', null=True)
        
    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        ordering = []
        filtering = {}
        serializer = CSVSerializer()
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
    

class ScreenResource(MetahashManagedResource,PostgresSortingResource):

#    lab_head_full = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=True) #, full_list=False) #, blank=True, null=True)
    lab_head_link = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=False)
    lead_screener_link = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=False)
    
    class Meta:
        queryset = Screen.objects.all().filter()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        resource_name = 'screen'
        
        ordering = []
        filtering = {}
        serializer = CSVSerializer()

        
    def __init__(self, **kwargs):
        self.scope = 'fields:screen'
        super(ScreenResource,self).__init__( **kwargs)

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
            if len(order_by) > 0:
                options.setlist('order_by', order_by)
            else:
                del options['order_by'] 
        logger.info(str(('options',options)))
        obj_list = super(ScreenResource, self).apply_sorting(obj_list, options)
        
        logger.info(str(('extra_order_by', extra_order_by)))
        obj_list = obj_list.order_by(*extra_order_by)
        return obj_list
