from django.conf.urls import url
from tastypie.authorization import Authorization
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
from tastypie import fields, utils

from db.models import ScreensaverUser,Screen, LabHead, LabAffiliation, ScreeningRoomUser
from lims.api import CSVSerializer, PostgresSortingResource
from reports.models import MetaHash, Vocabularies

import logging
        
logger = logging.getLogger(__name__)

# Mixin class - note this class must be mixed in first, since the tastypie Resource class does not call super().__init__
class MetahashManagedResource(object):
    def __init__(self, **kwargs):
#        if not self.scope:
#            self.scope = kwargs.pop('scope')
        logger.info(str(('---------init MetahashManagedResource', self.scope)))
        # TODO: research why calling reset_filtering_and_ordering, as below, fails        
        metahash = MetaHash.objects.get_metahash(scope=self.scope)
        
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)
        super(MetahashManagedResource,self).__init__( **kwargs)

    def reset_filtering_and_ordering(self):
        self.Meta.filtering = {}
        self.Meta.ordering = []
        metahash = MetaHash.objects.get_metahash(scope=self.scope)
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)

        logger.info(str(('+++filtering', self.Meta.filtering)))
        logger.info(str(('ordering', self.Meta.ordering)))

    def build_schema(self):
        logger.info('--- build_schema: ' + self.scope )
        schema = super(MetahashManagedResource,self).build_schema()  # obligatory super call, this framework does not utilize
        metahash = MetaHash.objects.get_metahash(scope=self.scope)

        for key, value in metahash.items():
            if key not in schema['fields']:
                logger.info('creating a virtual field: ' + key)
                schema['fields'][key] = {}
            schema['fields'][key].update(value)
            
                
        try:
            logger.info(str(('trying to locate resource information', self._meta.resource_name, self.scope)))
            resource_def = MetaHash.objects.get(scope='resource', key=self._meta.resource_name)
        
            schema['resource_definition'] = resource_def.model_to_dict(scope='fields:resource')
        except Exception, e:
            logger.warn(str(('on trying to locate resource information', e, self._meta.resource_name)))
                
        return schema
    
class ScreensaverUserResource(MetahashManagedResource, PostgresSortingResource):
    screens = fields.ToManyField('db.api.ScreenResource', 'screens', related_name='lab_head_id', blank=True, null=True)

    version = fields.IntegerField(attribute='version', null=True)
        
    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        ordering = []
        filtering = {}
        serializer = CSVSerializer()
        
    def __init__(self, **kwargs):
        self.scope = 'fields:screensaveruser'
        super(ScreensaverUserResource,self).__init__( **kwargs)
  
    def dehydrate(self, bundle):
        bundle = super(ScreensaverUserResource, self).dehydrate(bundle);
        bundle.data['screens'] = [x.facility_id for x in Screen.objects.filter(lab_head_id=bundle.obj.screensaver_user_id)]
        
        return bundle        
    
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

    lab_head_full = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=True) #, full_list=False) #, blank=True, null=True)
    lab_head = fields.ToOneField('db.api.LabHeadResource', 'lab_head',  full=True)
    
    class Meta:
        queryset = Screen.objects.all().filter(screen_type='Small Molecule')
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
#        bundle = super(ScreenResource, self).dehydrate(bundle);
        sru = bundle.obj.lead_screener.screensaver_user
        bundle.data['lead_screener'] =  sru.first_name + ' ' + sru.last_name
#        logger.info('lead_screener: ' + bundle.data['lead_screener'])
        lh = bundle.obj.lab_head.screensaver_user.screensaver_user
        bundle.data['lab_head'] =  lh.first_name + ' ' + lh.last_name
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
