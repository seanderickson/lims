from django.db import connection, DatabaseError

from reports.models import FieldInformation, MetaHash

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
from django.core.serializers.json import DjangoJSONEncoder
from django.core.serializers import json
import json,re
        
logger = logging.getLogger(__name__)


class BackboneSerializer(Serializer):
    
    def from_json(self, content):
        """
        Given some JSON data, returns a Python dictionary of the decoded data.
        """
        logger.info(str(("loading content:", content)))
        # quote attributes - the backbone client doesn't want to do this
        content = content.replace(r'(\w+):', r'"\1" :')
        logger.info(str(("loading content:", content)))
        return json.loads(content)
    
class FieldInformationResource(PostgresSortingResource):

    class Meta:
        queryset = FieldInformation.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {}
        serializer = BackboneSerializer()

class MetaHashResource(PostgresSortingResource):

    class Meta:
        queryset = MetaHash.objects.all()
        authentication = MultiAuthentication(BasicAuthentication(), SessionAuthentication())
        authorization= Authorization()        
        # TODO: drive this from data
        ordering = []
        filtering = {}
        serializer = BackboneSerializer()
