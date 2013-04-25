from django.db import connection, DatabaseError

from db.models import ScreensaverUser
from django.conf.urls import url

from tastypie.authorization import Authorization
from tastypie.bundle import Bundle
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS

class ScreensaverUserResource(ModelResource):
    class Meta:
        queryset = ScreensaverUser.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        ordering = ['first_name', 'last_name', 'comments']
        filtering = {'first_name':ALL, 'last_name':ALL, 'comments':ALL }
        
    def dehydrate(self, bundle):
        return bundle
        
        # to override: resource_name = 'sm'
    # def prepend_urls(self):
    #     return [
    #         url(r"^(?P<resource_name>%s)/(?P<facility_id>\d+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
    #     ]

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