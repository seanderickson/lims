from django.db import connection, DatabaseError

from reports.models import FieldInformation

from django.conf.urls import url

from tastypie.authorization import Authorization
from tastypie.bundle import Bundle
from tastypie.resources import ModelResource, Resource
from tastypie.serializers import Serializer
from tastypie.authentication import BasicAuthentication, SessionAuthentication, MultiAuthentication
from tastypie.constants import ALL, ALL_WITH_RELATIONS
from tastypie import fields, utils

from db.api import PostgresSortingResource

import logging
        
logger = logging.getLogger(__name__)

class BackboneCompatibleResource(PostgresSortingResource):

    class Meta:
        always_return_data = True

    def alter_list_data_to_serialize(self, request, data):
        return data["objects"]

class FieldInformationResource(BackboneCompatibleResource):

    class Meta:
        queryset = FieldInformation.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {}
