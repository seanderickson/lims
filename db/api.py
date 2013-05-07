from django.db import connection, DatabaseError

from db.models import ScreensaverUser
from django.conf.urls import url

from tastypie.authorization import Authorization
from tastypie.bundle import Bundle
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
import logging
        
logger = logging.getLogger(__name__)

class ScreensaverUserResource(ModelResource):
    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        ordering = ['first_name', 'last_name', 'comments', 'date_loaded', 'date_publicly_available']
        filtering = {'first_name':ALL, 'last_name':ALL, 'comments':ALL }
        
    def dehydrate(self, bundle):
        return bundle

    def build_schema(self):
        schema = super(ScreensaverUserResource,self).build_schema()
        
        field_meta = {
            'first_name': { 'name':'First', 'detail_view': True, 'list_view': True, 'order': 1 },
            'last_name': { 'name':'Last', 'detail_view': True, 'list_view': True, 'order': 2 },
            'ecommons_id': { 'name':'eCommons', 'detail_view': True, 'list_view': True, 'order': 3 },
            'login_id': { 'name':'Login', 'detail_view': True, 'list_view': True, 'order': 4 },
            'harvard_id': { 'name':'HUID', 'detail_view': True, 'list_view': True, 'order': 5 },
            'mailing_address': { 'name':'Address', 'detail_view': True, 'list_view': False, 'order': 6 },
            'phone': { 'name':'Phone', 'detail_view': True, 'list_view': True, 'order': 7 },
            'email': { 'name':'Email', 'detail_view': True, 'list_view': True, 'order': 8 },
            'comments': { 'name':'Comments', 'detail_view': True, 'list_view': False, 'order': 9 },
            'screensaver_user_id': { 'name':'ID', 'detail_view': True, 'list_view': True, 'order': 10 },
            'date_created': { 'name':'Date Created', 'detail_view': True, 'list_view': False, 'order': 100 },
            'date_loaded': { 'name':'Date Loaded', 'detail_view': True, 'list_view': False, 'order': 100 },
            'date_publicly_available': { 'name':'Publicly Available', 'detail_view': True, 'list_view': False, 'order': 100 },
            'digested_password': { 'name':'digested_password', 'detail_view': False, 'list_view': False, 'order': 100 },
            'harvard_id_expiration_date': { 'name':'HUID exp.', 'detail_view': True, 'list_view': False, 'order': 100 },
            'harvard_id_requested_expiration_date': { 'name':'HUID requested exp.', 'detail_view': True, 'list_view': False, 'order': 100 },
            'resource_uri': { 'name':'resource_uri', 'detail_view': False, 'list_view': False, 'order': 100 },
            'version': { 'name':'version', 'detail_view': False, 'list_view': False, 'order': 100 },
        }
        
        for field in schema['fields'].keys():
            if field in field_meta:
                schema['fields'][field].update(field_meta[field])
        return schema
    
#    def build_filters(self, filters=None):
#        logger.info('build_filters')
#        if filters is None:
#            filters = {}
#
#        orm_filters = super(ScreensaverUserResource, self).build_filters(filters)
#
#        return orm_filters
    
#    def get_object_list(self, request):
#        logger.info('obj_get_list')
#        q = super(ScreensaverUserResource, self).get_object_list(request)
#        return q
    
    def apply_sorting(self, obj_list, options):
        """
        Create a non-too-pretty workaround for the postgresql null sorting issue - nulls sort higher than values, 
        which is not desired.  We want nulls to sort lower than values.
        """ 
        
        obj_list = super(ScreensaverUserResource, self).apply_sorting(obj_list, options)
        
        extra_select = {}
        extra_ordering = []
        for field in obj_list.query.order_by:
            if field.startswith('-'):
                field = field[1:]
            extra_select[field+"_null"]=field + ' is null'
            extra_ordering.append(field+"_null")
        logger.info(str(('extra_select', extra_select)))
        obj_list = obj_list.extra(extra_select)

        # Note that this doesn't work, something in the framework deletes the extra 
        # order_by clause when apply_sorting, or, if this is run last, it deletes the sorting applied in apply_sorting...
        #        obj_list = obj_list.extra(order_by=['-comments_null'])

        # note: this doesn't work because the "is null" field order by clauses
        # must be prepended so that they occur before their intended fields
        #        obj_list.query.add_ordering('comments_null')
        
        temp = obj_list.query.order_by;
        obj_list.query.clear_ordering()
        for xfield in extra_ordering:
            temp.insert(0,xfield)
        obj_list.query.add_ordering(*temp)
        
        return obj_list
    
    
#    def dehydrate(self, bundle):
#        _filter = lambda field_information: not bundle.obj.is_restricted or field_information.is_unrestricted # or is_authorized
#        bundle.data = get_detail_bundle(bundle.obj, ['smallmolecule',''], _filter=_filter)
#        return bundle
#
#    def build_schema(self):
#        schema = super(SmallMoleculeResource,self).build_schema()
#        schema['fields'] = get_detail_schema(SmallMolecule(),['smallmolecule'])
#        return schema 
#    
#    def override_urls(self):
#        """ Note, will be deprecated in >0.9.12; delegate to new method, prepend_urls
#        """
#        return self.prepend_urls();
#    
#    def prepend_urls(self):
#
#        return [
#            url(r"^(?P<resource_name>%s)/(?P<facility_id>\d+)\-(?P<salt_id>\d+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
#        ]

    def alter_list_data_to_serialize(self, request, data):
        """
        A hook to alter list data just before it gets serialized & sent to the user.

        Useful for restructuring/renaming aspects of the what's going to be
        sent.

        Should accommodate for a list of objects, generally also including
        meta data.
        """
        sEcho = request.GET.get('sEcho','')
        if sEcho:
            data['meta']['sEcho'] = int(sEcho)
        return data