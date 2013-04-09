from django.db import connection, DatabaseError

from db.models import ScreensaverUser
from django.conf.urls.defaults import url

from tastypie.authorization import Authorization
from tastypie.bundle import Bundle
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer

class ScreensaverUserResource(ModelResource):
    class Meta:
        queryset = ScreensaverUser.objects.all()
        # to override: resource_name = 'sm'
    # def prepend_urls(self):
    #     return [
    #         url(r"^(?P<resource_name>%s)/(?P<facility_id>\d+)/$" % self._meta.resource_name, self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
    #     ]
