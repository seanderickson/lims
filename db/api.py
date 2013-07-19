#from django.db import connection, DatabaseError

from db.models import ScreensaverUser,Screen, LabHead, LabAffiliation, ScreeningRoomUser
from reports.models import MetaHash, Vocabularies
from lims.api import BackboneSerializer, PostgresSortingResource

#from django.conf.urls import url

from tastypie.authorization import Authorization
#from tastypie.bundle import Bundle
#from tastypie.resources import ModelResource, Resource
#from tastypie.serializers import Serializer
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
from tastypie import fields, utils

import logging
        
logger = logging.getLogger(__name__)

class ScreensaverUserResource(PostgresSortingResource):
    screens = fields.ToManyField('db.api.ScreenResource', attribute='screens', related_name='lab_head', blank=True, null=True)

    version = fields.IntegerField(attribute='version', null=True)
#    date_created = fields.DateTimeField(readonly=True, help_text='When the person was created', null=True)
        
    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        ordering = []
        filtering = {}
        serializer = BackboneSerializer()
        
    def __init__(self, **kwargs):
        
        metahash = MetaHash.objects.get_metahash(scope='screensaveruser:fields')
        
        for key,hash in metahash.items():
            if 'filtering' in hash and hash['filtering']:
                self.Meta.filtering[key] = ALL
        
        for key,hash in metahash.items():
            if 'ordering' in hash and hash['ordering']:
                self.Meta.ordering.append(key)
        super(ScreensaverUserResource,self).__init__(**kwargs)
        
    
    def build_schema(self):
        
        schema = super(ScreensaverUserResource,self).build_schema()
        
        scope = 'screensaveruser:fields'

        metahash = MetaHash.objects.get_metahash(scope='screensaveruser:fields')
                
        for key, value in schema['fields'].items():
            if key in metahash:
                value.update(metahash[key])
        
            
#        # first, query the metahash for fields defined for this scope
#        for fieldinformation in MetaHash.objects.all().filter(scope=scope):
#            logger.info('---- meta field for scope: ' + scope + ', ' + fieldinformation.key)
#            
#            field_key = fieldinformation.key
#            hash = schema['fields'][field_key]
#            for meta_record in MetaHash.objects.all().filter(scope='metahash:fields'):  # metahash:fields are defined for all reports
#                hash.update({
#                      meta_record.key: MetaHash.objects.get_or_none(scope=scope, key=field_key, function=lambda x : (x.get_field(meta_record.key)) )
#                      })
#                    
#            # now check if the field uses controlled vocabulary, look that up now.  TODO: "vocabulary_scope_ref" should be a constant
#            # TODO: "vocabulary_scope_ref" needs to be created by default as a metahash:field; this argues for making it a "real" field
#            if hash.get(u'vocabulary_scope_ref'):
#                vocab_ref = hash['vocabulary_scope_ref']
#                logger.info(str(('looking for a vocabulary', vocab_ref )))
#                hash['choices'] = [x.key for x in Vocabularies.objects.all().filter(scope=vocab_ref)]
#                logger.info(str(('got', hash['choices'] )))
        

# OLDER
#        # TODO: drive this from data
#        field_meta = {
#            'first_name': { 'name':'First', 'detail_view': True, 'list_view': True, 'order': 1 },
#            'last_name': { 'name':'Last', 'detail_view': True, 'list_view': True, 'order': 2 },
#            'ecommons_id': { 'name':'eCommons', 'detail_view': True, 'list_view': True, 'order': 3 },
#            'login_id': { 'name':'Login', 'detail_view': True, 'list_view': True, 'order': 4 },
#            'harvard_id': { 'name':'HUID', 'detail_view': True, 'list_view': True, 'order': 5 },
#            'mailing_address': { 'name':'Address', 'detail_view': True, 'list_view': False, 'order': 6 },
#            'phone': { 'name':'Phone', 'detail_view': True, 'list_view': True, 'order': 7 },
#            'email': { 'name':'Email', 'detail_view': True, 'list_view': True, 'order': 8 },
#            'comments': { 'name':'Comments', 'detail_view': True, 'list_view': False, 'order': 9 },
#            'screensaver_user_id': { 'name':'ID', 'detail_view': True, 'list_view': True, 'order': 10 },
#            'date_created': { 'name':'Date Created', 'detail_view': True, 'list_view': False, 'order': 100 },
#            'date_loaded': { 'name':'Date Loaded', 'detail_view': True, 'list_view': False, 'order': 100 },
#            'date_publicly_available': { 'name':'Publicly Available', 'detail_view': True, 'list_view': False, 'order': 100 },
#            'digested_password': { 'name':'digested_password', 'detail_view': False, 'list_view': False, 'order': 100 },
#            'harvard_id_expiration_date': { 'name':'HUID exp.', 'detail_view': True, 'list_view': False, 'order': 100 },
#            'harvard_id_requested_expiration_date': { 'name':'HUID requested exp.', 'detail_view': True, 'list_view': False, 'order': 100 },
#            'resource_uri': { 'name':'resource_uri', 'detail_view': False, 'list_view': False, 'order': 100 },
#            'version': { 'name':'version', 'detail_view': False, 'list_view': False, 'order': 100 },
#        }
#        
#        for field in schema['fields'].keys():
#            if field in field_meta:
#                schema['fields'][field].update(field_meta[field])
        return schema
    
    
#    def alter_list_data_to_serialize(self, request, data):
#        """
#        A hook to alter list data just before it gets serialized & sent to the user.
#
#        Useful for restructuring/renaming aspects of the what's going to be
#        sent.
#
#        Should accommodate for a list of objects, generally also including
#        meta data.
#        """
#        
#        # TODO: this is for the datatables jquery component
#        sEcho = request.GET.get('sEcho','')
#        if sEcho:
#            data['meta']['sEcho'] = int(sEcho)
#        return data

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
    
    lab_affiliation = fields.ToOneField('db.api.LabAffiliationResource', attribute='lab_affiliation',  full=True)
    
    # rather than walk the inheritance hierarchy, will flatten this hierarchy in the dehydrate method
    #    screening_room_user = fields.ToOneField('db.api.ScreeningRoomUserResource', attribute='screensaver_user',  full=True)
    
    id = fields.IntegerField(attribute='screensaver_user_id')
    
    class Meta:
        queryset = LabHead.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
    
    def dehydrate(self, bundle):
        # flatten the inheritance hierarchy, rather than show nested "lab_head->screening_room_user->screensaver_user"
        bundle.data.update(model_to_dict(bundle.obj.screensaver_user))
        bundle.data.update(model_to_dict(bundle.obj.screensaver_user.screensaver_user))
        
        return bundle        
    

class ScreenResource(PostgresSortingResource):
    lab_head_full = fields.ToOneField('db.api.LabHeadResource', attribute = 'lab_head',  full=True) #, full_list=False) #, blank=True, null=True)
    lab_head_id = fields.IntegerField(attribute='lab_head_id')
    
    class Meta:
        queryset = Screen.objects.all().filter(screen_type='Small Molecule')
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        
        # TODO: drive this from data
        ordering = ['facility_id', 'title', 'lead_screener.screensaver_user.last_name', 'date_loaded', 'date_publicly_available']
        # TODO: drive this from data
        filtering = {'first_name':ALL, 'last_name':ALL, 'comments':ALL }
        
    def dehydrate(self, bundle):
        sru = bundle.obj.lead_screener.screensaver_user
        bundle.data['lead_screener'] =  sru.first_name + ' ' + sru.last_name
        lh = bundle.obj.lab_head.screensaver_user.screensaver_user
        bundle.data['lab_head'] =  lh.first_name + ' ' + lh.last_name
        
#        logger.info(str(('bundle: ',bundle)))
        return bundle

    def build_schema(self):
        schema = super(ScreenResource,self).build_schema()
        
        # TODO: drive this from data
        field_meta = {
            'screen_id': { 'name':'ID', 'detail_view': False, 'list_view': False, 'order': 100 },
            'facility_id': { 'name':'ID', 'detail_view': True, 'list_view': True, 'order': -1 },
            'screen_type': { 'name':'Screen Type', 'detail_view': True, 'list_view': True, 'order': 1 },
            'project_phase': { 'name':'Project Phase', 'detail_view': True, 'list_view': True, 'order': 2 },
            'project_id': { 'name':'Project ID', 'detail_view': True, 'list_view': True, 'order': 3 },
            'title': { 'name':'Title', 'detail_view': True, 'list_view': True, 'order': 4 },
            'lab_head': { 'name':'Lab Head', 'detail_view': True, 'list_view': True, 'order': 5 },
            'lead_screener': { 'name':'Lead Screener', 'detail_view': True, 'list_view': True, 'order': 6 },
            'date_recorded': { 'name':'Date Recorded', 'detail_view': True, 'list_view': True, 'order': 7 },
            'status': { 'name':'Status', 'detail_view': True, 'list_view': True, 'order': 8 },
            'status_date': { 'name':'Status', 'detail_view': True, 'list_view': True, 'order': 9 },
            'screen_result': { 'name':'Screen Result', 'detail_view': True, 'list_view': True, 'order': 10 },
            'total_plated_lab_cherry_picks': { 'name':'Total Cherry Picks', 'detail_view': True, 'list_view': True, 'order': 11 },
        }

        for field in field_meta.keys():
            if field not in schema['fields']:
                schema['fields'][field] = field_meta[field]
                
        for field in schema['fields'].keys():
            if field in field_meta:
                schema['fields'][field].update(field_meta[field])
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
