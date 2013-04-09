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
