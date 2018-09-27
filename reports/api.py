# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from collections import defaultdict
from copy import deepcopy
from decimal import Decimal
import decimal
from functools import wraps
import importlib
import json
import logging
from operator import itemgetter
import os
import re
import subprocess
import sys
import time
import urllib

from aldjemy.core import get_tables, get_engine
from django.conf import settings
from django.conf.global_settings import DEFAULT_CHARSET
from django.conf.urls import url
from django.contrib.auth.models import User as DjangoUser
from django.contrib.sessions.models import Session
from django.core.cache import caches
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.core.serializers.json import DjangoJSONEncoder
from django.db import transaction
from django.db.models import Q
from django.db.utils import ProgrammingError
from django.forms.models import model_to_dict
from django.http.request import HttpRequest
from django.http.response import HttpResponse, Http404
from django.test.client import RequestFactory
import six
from sqlalchemy import select, asc, text
from sqlalchemy.dialects import postgresql
from sqlalchemy.dialects.postgresql import array
from sqlalchemy.sql import and_, or_, not_, asc, desc, func
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join, distinct, exists
from sqlalchemy.sql.selectable import Alias

from db.support import lims_utils
from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    HTTP_PARAM_USE_TITLES, HTTP_PARAM_USE_VOCAB, HEADER_APILOG_COMMENT, \
    HTTP_PARAM_DATA_INTERCHANGE, InformationError, API_RESULT_ERROR
from reports import ValidationError, BadRequestError, ApiNotImplemented, \
    MissingParam, BackgroundJobImmediateResponse, _now
from reports.api_base import IccblBaseResource, un_cache, Authorization, \
    MultiAuthentication, IccblSessionAuthentication, IccblBasicAuthentication, \
    TRAILING_SLASH
from reports.models import MetaHash, Vocabulary, ApiLog, LogDiff, Permission, \
                           UserGroup, UserProfile, Job
from reports.utils import default_converter
import reports.schema as SCHEMA
from reports.serialize import parse_val, parse_json_field, XLSX_MIMETYPE, \
    SDF_MIMETYPE, XLS_MIMETYPE, JSON_MIMETYPE, MULTIPART_MIMETYPE, \
    LimsJSONEncoder
from reports.serializers import LimsSerializer, DJANGO_ACCEPT_PARAM
from reports.sqlalchemy_resource import SqlAlchemyResource, _concat
import reports.utils.background_client_util as background_client_util
import reports.utils.background_processor as background_processor
import reports.utils.si_unit as si_unit

logger = logging.getLogger(__name__)

URI_VERSION = 'v1'
BASE_URI = '/reports/api/' + URI_VERSION

# from reports.schema import API_RESULT_META
RESOURCE = SCHEMA.RESOURCE
API_RESULT_OBJ = SCHEMA.API_RESULT_OBJ
API_RESULT_DATA = SCHEMA.API_RESULT_DATA
API_RESULT_META = SCHEMA.API_RESULT_META
API_ACTION = SCHEMA.VOCAB.apilog.api_action


API_PARAM_OVERRIDE = 'override'
API_PARAM_PATCH_PREVIEW_MODE = 'patch_with_preview'
API_PARAM_SHOW_PREVIEW = 'show_preview'
API_PARAM_PREVIEW_LOGS = 'preview_logs'
API_PARAM_NO_BACKGROUND = 'no_background'

DEBUG_RESOURCES = False or logger.isEnabledFor(logging.DEBUG)
DEBUG_FIELDS = False or logger.isEnabledFor(logging.DEBUG)
DEBUG_AUTHORIZATION = False or logger.isEnabledFor(logging.DEBUG)
DEBUG_PATCH_LOG = False or logger.isEnabledFor(logging.DEBUG)


class UserGroupAuthorization(Authorization):
    
    def __init__(self, resource_name, *args, **kwargs):
        Authorization.__init__(self, *args, **kwargs)
        self.resource_name = resource_name
        
    @staticmethod
    def get_authorized_resources(user, permission_type):
        userprofile = user.userprofile
        permission_types = [permission_type]
        if permission_type == 'read':
            permission_types.append('write')
        resources_user = ( userprofile.permissions.all()
            .filter(scope='resource', type__in=permission_types)
            .values_list('key', flat=True))
        resources_group = [ permission.key 
                for group in userprofile.usergroup_set.all() 
                for permission in group.get_all_permissions(
                    scope='resource', type__in=permission_types)]
        return set(resources_user) | set(resources_group)
    
    def _is_resource_authorized(
        self, user, permission_type, resource_name=None, **kwargs):
        
        if not user:
            return False
        if DEBUG_AUTHORIZATION:
            logger.info("_is_resource_authorized: %s, user: %s, type: %s",
                self.resource_name, user, permission_type)
        scope = 'resource'
        prefix = 'permission'
        uri_separator = '/'
        
        if resource_name is None:
            resource_name = self.resource_name
        
        permission_str =  uri_separator.join([
            prefix,scope,resource_name,permission_type])       

        if DEBUG_AUTHORIZATION:
            logger.info('authorization query: %s, user %s, %s' 
                % (permission_str, user, user.is_superuser))
        
        if user.is_superuser:
            if DEBUG_AUTHORIZATION:
                logger.info('%s:%s access allowed for super user: %s' 
                    % (resource_name,permission_type,user))
            return True
        
        userprofile = user.userprofile
        permission_types = [permission_type]
        if permission_type == 'read':
            permission_types.append('write')
        query = userprofile.permissions.all().filter(
            scope=scope, key=resource_name, type__in=permission_types)
        if query.exists():
            if DEBUG_AUTHORIZATION:
                logger.info(
                    'user %s, auth query: %s, found matching user permissions %s'
                    % (user,permission_str,[str(x) for x in query]))
            logger.info('%s:%s user explicit permission for: %s' 
                % (resource_name,permission_type,user))
            return True
        
        if DEBUG_AUTHORIZATION:
            logger.info(
                'user %s, auth query: %s, not found in user permissions %s'
                % (user,permission_str,[str(x) for x in query]))
        
        permissions_group = [ permission 
                for group in userprofile.usergroup_set.all() 
                for permission in group.get_all_permissions(
                    scope=scope, key=resource_name, type__in=permission_types)]
        if permissions_group:
            if DEBUG_AUTHORIZATION:
                logger.info(
                    'user: %r, auth query: %r, found usergroup permissions: %r'
                    ,user, permission_str, permissions_group)
            logger.debug('%s:%s usergroup permission for: %s' 
                % (resource_name,permission_type,user))
            return True
        
        logger.info(
            'user: %r, auth query: %r, not authorized'
            ,user.username, permission_str)
        return False

    def is_restricted_view(self, user):
        '''
        True if the user does not have full "read" privilege on the resource
        '''
        authorized = UserGroupAuthorization._is_resource_authorized(self, user, 'read')
        if DEBUG_AUTHORIZATION:
            logger.info('is restricted: %r:%r, %r', 
                self.resource_name, user, not authorized)
        return not authorized 
    
    def get_fields_by_level(self, fields):
        '''
        Group the field keys by the data_access_level assigned in the schema:
        - if a data_access_level has not been assigned, the default value is
        used ( level 3 - visible only for own screens)
        '''
        
        # If the data_access_level is not set, use the default value
        default_field_data_access_level = 3
        
        restricted_fields = set()
        fields_by_level = defaultdict(set)
        for key,field in fields.items():
            if field['view_groups']: # these fields should already be elided
                logger.warn('field %r is only visible to groups: %r',
                    key, field['view_groups'])
                continue
            field_data_access_level = field.get('data_access_level', None)
            if field_data_access_level is None:
                field_data_access_level = default_field_data_access_level
            fields_by_level[field_data_access_level].add(key)
            if field_data_access_level > 0:
                restricted_fields.add(key)
        if DEBUG_AUTHORIZATION:
            logger.info('fields by level: %r', fields_by_level)
            logger.info('property restricted_fields: %r', restricted_fields)
            
        return fields_by_level

    def get_effective_access_level(self, user, row):
        logger.warn('%r - get_effective_access_level must be implemented'
            'if implementing the access level property generator',
            self.resource_name)
        return None

    def filter(self, user, filter_expression):
        '''
        Filter the result rows returned
        @param filter_expression: a SqlAlchemy query filter
        TODO: override per resource if restricted access is granted for the 
        resource.
        '''
        return filter_expression

    def get_row_property_generator(self, user, fields, extant_generator):
        '''
        Filter result properties based on authorization rules
        '''
        return extant_generator
        # Must be implemented
        # e.g.
        # return self.get_access_level_property_generator(
        #    user, fields, extant_generator)
        
    def get_access_level_property_generator(self, user, fields, extant_generator):
        '''
        Filter results based on the access level granted for the user for each
        row
        '''

        if self.is_restricted_view(user) is False:
            return extant_generator
        else:
            logger.info('get_access_level_property_generator: %r user: %r', 
                self.resource_name, user)
    
            fields_by_level = self.get_fields_by_level(fields)
            if DEBUG_AUTHORIZATION:
                logger.info('fields by level: %r', fields_by_level)
            outer_self = self
            
            class Row:
                def __init__(self, row):
                    logger.debug(
                        'filter row: %r', 
                        [(key, row[key]) for key in row.keys()])
                    self.row = row
                    self.effective_access_level = outer_self.get_effective_access_level(user, row)

                    self.allowed_fields = set()
                    if self.effective_access_level is not None:
                        for level in range(0,self.effective_access_level+1):
                            self.allowed_fields.update(fields_by_level[level])
                            if DEBUG_AUTHORIZATION:
                                logger.info(
                                    'allow level: %r: %r, %r', 
                                    level, fields_by_level[level], 
                                    self.allowed_fields)
                    if DEBUG_AUTHORIZATION:
                        logger.info('allowed fields: %r', self.allowed_fields)
                def has_key(self, key):
                    if key == 'user_access_level_granted':
                        return True
                    return self.row.has_key(key)
                def keys(self):
                    return self.row.keys();
                def __getitem__(self, key):
                    logger.debug(
                        'key: %r, allowed: %r', key, key in self.allowed_fields)
                    if key == 'user_access_level_granted':
                        return self.effective_access_level
                    if self.row[key] is None:
                        return None
                    else:
                        if key in self.allowed_fields:
                            logger.debug('allow %r: %r', key, self.row[key])
                            return self.row[key]
                        else:
                            logger.debug(
                                'row: %r filter field: %r for restricted user: %r',
                                self.row, key, user.username)
                            return None

            def access_level_property_generator(cursor):
                if extant_generator is not None:
                    cursor = extant_generator(cursor)
                for row in cursor:
                    yield Row(row)
            
            return access_level_property_generator
        
    
class AllAuthenticatedReadAuthorization(UserGroupAuthorization):
    ''' Allow read for all authenticated users'''

    def _is_resource_authorized(
        self, user, permission_type, **kwargs):
        logger.debug('user: %r, %r, p: %r', user, self.resource_name, permission_type)
        if permission_type == 'read':
            if user is None:
                return False
            if DEBUG_AUTHORIZATION:
                logger.info('grant read access: %r: %r', user, self.resource_name)
            return True
        return super(AllAuthenticatedReadAuthorization,self)\
            ._is_resource_authorized(user, permission_type)

    
def write_authorization(_func):
    '''
    Wrapper function to verify write authorization
    ''' 
    @wraps(_func)
    def _inner(self, *args, **kwargs):
        request = args[0]
        resource_name = kwargs.pop('resource_name', self._meta.resource_name)
        if not self._meta.authorization._is_resource_authorized(
            request.user,'write', **kwargs):
            msg = 'write auth failed for user: %r, resource: %r' \
                % (request.user, resource_name)
            logger.warn(msg)
            raise PermissionDenied(msg)
        if resource_name != 'screenresult':
            kwargs['schema'] = self.build_schema(user=request.user)

        return _func(self, *args, **kwargs)

    return _inner


def read_authorization(_func):
    '''
    Wrapper function to verify read authorization
    ''' 
    @wraps(_func)
    def _inner(self, *args, **kwargs):
        request = args[0]
        resource_name = kwargs.pop('resource_name', self._meta.resource_name)
        if DEBUG_AUTHORIZATION:
            logger.info('read auth for: %r using: %r', 
                resource_name, self._meta.authorization)
        if not self._meta.authorization._is_resource_authorized(
                request.user,'read', **kwargs):
            msg = 'read auth failed for: %r, %r' % (request.user, resource_name)
            logger.warn(msg)
            raise PermissionDenied(msg)

        # Use the user specific settings to build the schema
        # NOTE: screenresult schema is generated by the resource class
        if resource_name not in ['screenresult']:
            logger.info('build schema for: %r', resource_name)
            kwargs['schema'] = self.build_schema(user=request.user)
        
        return _func(self, *args, **kwargs)

    return _inner


def background_job(_func):
    '''
    Wrapper function to defer request to offline background job processor:
    
    - if SCHEMA.JOB.JOB_PROCESSING_FLAG is set as a request parameter, 
    service the job_id indicated,
    - else create a new [pending, submitted] Job
    
    @see reports.api.JobResource
    
    raises BackgroundJobImmediateResponse with result job_data
    ''' 
    @wraps(_func)
    def _inner(self, *args, **kwargs):
        request = args[0]
        
        
        no_background = kwargs.get(API_PARAM_NO_BACKGROUND, False)
        no_background = parse_val(no_background, API_PARAM_NO_BACKGROUND, 'boolean')
        
        if settings.BACKGROUND_PROCESSING is not True:
            logger.info(
                'background processing is turned off because of '
                'settings.BACKGROUND_PROCESSING')
            return _func(self, *args, **kwargs)
        elif no_background is True:
            logger.info('background processing is turned off because '
                        'API_PARAM_NO_BACKGROUND: %r is set', API_PARAM_NO_BACKGROUND)
            return _func(self, *args, **kwargs)
        else:
            logger.info('background processing request...')
           
            job_resource = JobResource()
            JOB = SCHEMA.JOB

            # JOB_PROCESSING_FLAG: indicates the job_id to service
            job_id = request.GET.get(JOB.JOB_PROCESSING_FLAG, None)
            if job_id is None:
                job_id = request.POST.get(JOB.JOB_PROCESSING_FLAG, None)
            
            if job_id is not None:
                
                logger.info('JOB_PROCESSING_FLAG is set - Service the job: %r', job_id)
                job_data = \
                    job_resource.service_job(job_id, self, _func, **kwargs)

                meta = { JOB.resource_name: job_data }
                data = { API_RESULT_META: meta }
                job_response = self.build_response(
                    request, data, status_code=200)
                                
                raise BackgroundJobImmediateResponse(job_response)
            else:
                
                logger.info('JOB_PROCESSING_FLAG is not set: create a new job...')
                new_job_data = \
                    job_resource.create_new_job_from_request(request, **kwargs)
                job_id = new_job_data[JOB.ID]
                logger.info(
                    'job created: %r, submit to background processor...', job_id)
                
                new_job_data = job_resource.submit_job(new_job_data)
                meta = { 'job': new_job_data }
                data = { API_RESULT_META: meta }
                response = self.build_response(request, data, status_code=201)
                                 
                raise BackgroundJobImmediateResponse(response)
            
    return _inner


class SuperUserAuthorization(Authorization):

    def _is_resource_authorized(
        self, user, permission_type, **kwargs):

        if DEBUG_AUTHORIZATION:
            logger.info('Super User Authorization for: %r, %r', 
                user, user.is_superuser)
        if user.is_superuser:
            return True
        
        return False
    

def compare_dicts(dict1, dict2, excludes=None, exclude_patterns=None, 
    log_empty_strings=False):
    '''
    @param log_empty_strings if set, then log when a field is changed from
    None to an empty string, or vice versa:
    NOTE: results from 1) default fields in the DB (should be cleaned up),
    2) API string None fields returning string (should be cleaned up)
    '''
    if DEBUG_PATCH_LOG:
        logger.info('compare dicts: %r, %r, excludes: %r, %r, log_empty_strings: %r', 
            dict1, dict2, excludes, exclude_patterns, log_empty_strings)
    
    _excludes = set(['resource_uri'])
    if excludes:
        _excludes = _excludes.union(set(excludes))
    if exclude_patterns:
        for exclude_pattern in exclude_patterns:
            _excludes.update([key for key in dict1.keys() if exclude_pattern in key ])
            _excludes.update([key for key in dict2.keys() if exclude_pattern in key ])
    original_keys = set(dict1.keys())-_excludes
    updated_keys = set(dict2.keys())-_excludes
    logger.debug('compare dicts, updated keys: %r', updated_keys)
    
    union_keys = original_keys.union(updated_keys)
    if DEBUG_PATCH_LOG:
        logger.info('union keys to check: %r', union_keys)
    diffs = {}
    for key in union_keys:
        v1 = dict1.get(key, None)
        v2 = dict2.get(key, None)

        if v1 != v2:
            if log_empty_strings is True:
                diffs[key] = [v1,v2]
            else:
                if v1 is None and str(v2) == '':
                    pass
                elif v2 is None and str(v1) == '':
                    pass
                else:
                    diffs[key] = [v1, v2]
    if DEBUG_PATCH_LOG:
        logger.info('diffs: %r', diffs)
    return diffs

    
def is_empty_diff(difflog):
    if not difflog:
     return True
    
    empty = True;
    for key, value in difflog.items():
        if value:
            empty = False;
    return empty

        
def download_tmp_file(path, filename):
    """                                                                         
    Send a file through Django without loading the whole file into              
    memory at once. The FileWrapper will turn the file object into an           
    iterator for chunks of 8KB.                                                 
    """
    try:
        _file = file(_path)
        wrapper = FileWrapper(_file)

        # use the same type for all files
        response = HttpResponse(wrapper, content_type='text/plain') 
        response['Content-Disposition'] = \
            'attachment; filename=%s' % unicode(filename)
        response['Content-Length'] = os.path.getsize(_path)
        return response
    except Exception,e:
        logger.exception('could not find attached file object for id: %r', id)
        raise e


def get_supertype_fields(resource_definition):
    supertype = resource_definition.get('supertype', None)
    if supertype:
        temp = MetaHash.objects.get(
            scope='resource', key=supertype)
        super_resource_def = temp.model_to_dict(scope='fields.resource')
        fields = get_supertype_fields(super_resource_def)
        
        fields.update(deepcopy(
            MetaHash.objects.get_and_parse(
                scope='fields.%s' % supertype, 
                field_definition_scope='fields.field')))
        for field in fields.values():
            if not field['table']:
                field['table'] = super_resource_def['table']
        return fields
    else:
        return {}    


class ApiLogAuthorization(UserGroupAuthorization):
    '''
    FIXME: not used - rework ApiLog auth - 20160324
    Specialized authorization, allows users to read logs for resources they are 
    authorized for.
    '''
    def filter(self, user, filter_expression):
        # authorization filter
        if not user.is_superuser:
            # FIXME: "read" is too open
            # - grant read level access on a case-by-case basis
            # i.e. for screen.status updates
            resources = UserGroupAuthorization.get_authorized_resources(
                user, 'read')
            auth_filter = column('ref_resource_name').in_(resources)
            if filter_expression is not None:
                filter_expression = and_(filter_expression, auth_filter)
            else:
                filter_expression = auth_filter
        return filter_expression


#     def read_list(self, object_list, bundle, ref_resource_name=None, **kwargs):
#         if not ref_resource_name:
#             ref_resource_name = self.resource_meta.resource_name;
#         if self._is_resource_authorized(
#             ref_resource_name, bundle.request.user, 'read'):
#             return object_list
# 
#     def read_detail(
#             self, object_list, bundle, ref_resource_name=None, **kwargs):
#         if not ref_resource_name:
#             ref_resource_name = self.resource_meta.resource_name;
#         if self._is_resource_authorized(
#             ref_resource_name, bundle.request.user, 'read'):
#             return True


class ApiResource(SqlAlchemyResource):
    '''
    Provides framework for PATCH and PUT request methods:
    - wrap logging 
    '''
    
    def __init__(self, **kwargs):
        super(ApiResource,self).__init__(**kwargs)
        logger.debug('initialize ApiResource for %r', self._meta.resource_name)
        self.resource_resource = None
        self.vocab_resource = None
        
        self.request_factory = RequestFactory()
        
    def clear_cache(self, request, **kwargs):
        logger.debug('clearing the cache from resource: %s' 
            % self._meta.resource_name)
        ApiResource.get_cache(self).clear()
    
    def get_cache(self):
        return caches['reports_cache']

    def get_resource_resource(self):
        if self.resource_resource is None:
            self.resource_resource = ResourceResource()
        return self.resource_resource
    
    def get_vocab_resource(self):
        if self.vocab_resource is None:
            self.vocab_resource = VocabularyResource()
        return self.vocab_resource
    
    def get_schema(self, request, **kwargs):
        return self.build_response(
            request, 
            self.build_schema(user=request.user, full_schema=True, **kwargs ), 
            **kwargs)

    @read_authorization
    def get_detail(self, request, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'get_detail')
    
    @read_authorization
    def get_list(self, request, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'get_list')
        
    def build_list_response(self,request, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'build_list_response')

    def _get_list_response(self,request,**kwargs):
        '''
        Return a deserialized list of dicts
        '''
        logger.debug('_get_list_response: %r, %r', 
            self._meta.resource_name, 
            {k:v for k,v in kwargs.items() if k !='schema'})
        includes = kwargs.pop('includes', '')
        try:
            kwargs.setdefault('limit', 0)
            response = self.get_list(
                request,
                includes=includes,
                **kwargs)
            _data = self.get_serializer().deserialize(
                LimsSerializer.get_content(response), JSON_MIMETYPE)
            if self._meta.collection_name in _data:
                _data = _data[self._meta.collection_name]
            logger.debug(' data: %r', _data)
            return _data
        except Http404:
            return []
        except Exception as e:
            logger.exception('on get list: %r', e)
            raise
        
    def _get_detail_response(self,request,**kwargs):
        '''
        Return the detail response as a dict
        '''
        logger.info('_get_detail_response: %r, %r', 
            self._meta.resource_name, 
            {k:v for k,v in kwargs.items() if k !='schema'})
        includes = kwargs.pop('includes', '*')
        try:
            response = self.get_detail(
                request,
                includes=includes,
                **kwargs)
            _data = {}
            if response.status_code == 200:
                _data = self._meta.serializer.deserialize(
                    LimsSerializer.get_content(response), JSON_MIMETYPE)
            else:
                logger.info(
                    'no data found for %r, %r, %r', 
                    self._meta.resource_name, kwargs, response.status_code)
            return _data
        except Http404:
            return {}
        
    #@un_cache
    def _get_detail_response_internal(self, user=None, **kwargs):
        request = self.request_factory.generic(
            'GET', '.', HTTP_ACCEPT=JSON_MIMETYPE,
            **kwargs )
        if user is None:
            logger.debug('_get_detail_response_internal, no user')
            class User:
                is_superuser = True
                username = 'internal_request'
                def is_authenticated(self):
                    return True
            request.user = User()
        else:
            request.user = user
        result = self._get_detail_response(request, **kwargs)
        return result

    def _get_list_response_internal(self, user=None, **kwargs):
        request = self.request_factory.generic(
            'GET', '.', HTTP_ACCEPT=JSON_MIMETYPE )
        if user is None:
            logger.debug('_get_list_response_internal, no user')
            class User:
                is_superuser = True
                username = 'internal_request'
                def is_authenticated(self):
                    return True
            request.user = User()
        else:
            request.user = user
        result = self._get_list_response(request, **kwargs)
        return result

    def _patch_detail_internal(self, data, user=None, **kwargs):
        '''
        Note: requires the target resource to support the "data" kwarg; this
        bypasses deserialization for this request.
        '''

        request = self.request_factory.generic(
            'PATCH', '.', HTTP_ACCEPT=JSON_MIMETYPE )
        if user is None:
            logger.debug('_patch_detail_internal, no user')
            class User:
                id = 0
                is_superuser = True
                username = 'internal_request'
                def is_authenticated(self):
                    return True
            request.user = User()
        else:
            request.user = user
        
        logger.info('execute internal patch: %r', kwargs)
        
        if 'schema' not in kwargs:
            kwargs['schema'] = self.build_schema(user=user)
        
        kwargs['data'] = data
        
        response = self.patch_detail(
            request,
            format='json',
            **kwargs)
        _data = self.get_serializer().deserialize(
            LimsSerializer.get_content(response), JSON_MIMETYPE)
        return _data
    
    def _patch_list_internal(self, data, user=None, **kwargs):
        '''
        Note: requires the target resource to support the "data" kwarg; this
        bypasses deserialization for this request.
        '''

        request = self.request_factory.generic(
            'PATCH', '.', HTTP_ACCEPT=JSON_MIMETYPE )
        if user is None:
            logger.debug('_patch_detail_internal, no user')
            class User:
                id = 0
                is_superuser = True
                username = 'internal_request'
                def is_authenticated(self):
                    return True
            request.user = User()
        else:
            request.user = user
        
        logger.info('execute internal patch: %r', kwargs)
        
        if 'schema' not in kwargs:
            kwargs['schema'] = self.build_schema(user=user)
        
        kwargs['data'] = data
        
        response = self.patch_list(
            request,
            format='json',
            **kwargs)
        _data = self.get_serializer().deserialize(
            LimsSerializer.get_content(response), JSON_MIMETYPE)
        return _data
    
    def build_schema(self, user=None, **kwargs):
        
        logger.debug('build schema for: %r: %r', self._meta.resource_name, user)
        schema = self.get_resource_resource()._get_resource_schema(
            self._meta.resource_name, user, **kwargs)
        if DEBUG_RESOURCES:
            logger.info('schema fields: %r', schema[RESOURCE.FIELDS].keys())
        return schema

    
    def _get_filename(self, data_hash, schema, is_for_detail=False, **kwargs):
        
        id_parts = self.get_file_id(data_hash, schema, is_for_detail)
        
        if kwargs is not None:
            MAX_VAL_LENGTH = 20
            for key,val in kwargs.items():
                id_parts.append(str(key))
                val = default_converter(str(val))
                if val:
                    val = val[:MAX_VAL_LENGTH]
                    id_parts.append(val)
        
        id_parts.append(_now().strftime(SCHEMA.DATE_TIME_FILE_FORMAT))
        
        logger.debug(
            '_get_filename: %r, is_for_detail: %r', id_parts, is_for_detail)
        return '_'.join(id_parts)
        
    
    # def _get_filename_old(self, readable_filter_hash, schema, filename=None, **extra):
    #     MAX_VAL_LENGTH = 20
    #     file_elements = [self._meta.resource_name]
    #     if filename is not None:
    #         file_elements.append(filename)
    #     if extra is not None:
    #         for key,val in extra.items():
    #             file_elements.append(str(key))
    #             if val is not None:
    #                 val = default_converter(str(val))
    #                 val = val[:MAX_VAL_LENGTH]
    #                 file_elements.append(val)
    #     for key,val in readable_filter_hash.items():
    #         if key not in schema['id_attribute']:
    #             file_elements.append(str(key))
    #         val = default_converter(str(val))
    #         val = val[:MAX_VAL_LENGTH]
    #         file_elements.append(val)
    #             
    #     # if len(file_elements) > 1:
    #     #     # Add an extra separator for the resource name
    #     #     file_elements.insert(1,'_')
    #         
    #     filename = '_'.join(file_elements)
    #     logger.info('filename: %r', filename)
    #     MAX_FILENAME_LENGTH = 128
    #     filename = filename[:128]
    
    def get_file_id(self, data,  schema, is_for_detail):
        
        resource_name = schema.get('short_key')
        if not resource_name:
            if is_for_detail:
                resource_name = schema.get('title')
                resource_name = default_converter(resource_name)
            else:
                resource_name = schema.get('listing_title')
                if not resource_name:
                    resource_name = schema.get('title')
                resource_name = default_converter(resource_name)
        else:
            if not is_for_detail:
                resource_name += 's'
        id_parts = [resource_name]

        id_attributes = schema.get('file_id_attribute',schema.get('id_attribute'))
        
        logger.debug('id_attributes for file: %r, data: %r', id_attributes, data)
        
        if not id_attributes:
            return id_parts
        
        MAX_VAL_LENGTH = 20
        for id_attribute in id_attributes:
            if id_attribute in data:
                val = data[id_attribute]
                if val:
                    val = default_converter(str(val))
                    val = val[:MAX_VAL_LENGTH]
                    id_parts.append(val)
            elif id_attribute in schema:
                val = schema[id_attribute]
                if val:
                    val = default_converter(str(val))
                    val = val[:MAX_VAL_LENGTH]
                    id_parts.append(val)
        return id_parts
    
    def get_id(self,deserialized, validate=False, schema=None, **kwargs):
        ''' 
        return the full ID for the resource, as defined by the "id_attribute"
        - if validate=True, raises ValidationError if some or all keys are missing
        - otherwise, returns keys that can be found
        '''
        if schema is None:
            raise ProgrammingError
        
        id_attribute = schema['id_attribute']
        fields = schema[RESOURCE.FIELDS]
        logger.debug('get_id: %r, %r', self._meta.resource_name, id_attribute)
        kwargs_for_id = {}
        not_found = []
        for id_field in id_attribute:
            # NOTE: order of priority is kwargs, deserialized[resource_uri], deserialized
            # kwargs are checked first in case deserialized data represents a PATCH
            val = None
            if kwargs and kwargs.get(id_field,None):
                val = parse_val(
                    kwargs.get(
                        id_field,None), id_field,fields[id_field]['data_type']) 
            elif deserialized and deserialized.get(id_field,None):
                val = parse_val(
                    deserialized.get(
                        id_field,None), id_field,fields[id_field]['data_type']) 
            if val is not None:
                # ID values must all be non null/non-empty
                if len(str(val)) > 0:
                    kwargs_for_id[id_field] = val
            else:
                not_found.append(id_field)
            
        if not_found:
            if 'resource_uri' in deserialized:
                return self.find_key_from_resource_uri(
                    deserialized['resource_uri'], schema=schema)
        if not_found:
            if validate:
                logger.error(
                    'not all id fields found: %r: in %r, id_attribute: %r, resource: %r',
                    not_found,deserialized, id_attribute, schema['key'])
                raise ValidationError({
                    k:'required' for k in not_found })
        logger.debug('kwargs_for_id: %r', kwargs_for_id)   
        return kwargs_for_id

    def find_key_from_resource_uri(self,resource_uri, schema=None):
        if schema is None:
            logger.debug('create schema for find_key_from_resource_uri: %r', resource_uri)
            schema = self.build_schema()
        id_attribute = schema['id_attribute']
        resource_name = self._meta.resource_name + '/'
        logger.debug(
            'find_key_from_resource_uri: %r, resource_name: %r, id_attribute: %r', 
            resource_uri, resource_name, id_attribute)
         
        index = resource_uri.find(resource_name)
        if index > -1:
            index = index + len(resource_name)
            keystring = resource_uri[index:]
        else:
            keystring = resource_uri
        keys = keystring.strip('/').split('/')
        logger.debug('keys: %r, id_attribute: %r', keys, id_attribute)
        if len(keys) < len(id_attribute):
            raise BadRequestError({
                'resource uri': '%r does not contain all id attributes: %r'
                % (resource_uri,id_attribute)})
        else:
            return dict(zip(id_attribute,keys))

    def parse(self,deserialized, create=False, fields=None, schema=None):
        ''' parse schema fields from the deserialized dict '''
        DEBUG_PARSE = False or logger.isEnabledFor(logging.DEBUG)
        if DEBUG_PARSE:
            logger.info('parse: %r:%r', self._meta.resource_name, deserialized)

        if not deserialized or not isinstance(deserialized, dict):
            logger.warn('no deserialized data found')
            return {}
               
        errors = {}
        mutable_fields = fields        
        if fields is None:
            
            if schema is None:
                logger.info('build raw schema for parse')
                schema = self.build_schema()
        
            fields = schema[RESOURCE.FIELDS]
            mutable_fields = {}
            for key,field in fields.items():
                editability = field.get('editability', None)
                if editability and (
                    'u' in editability or (create and 'c' in editability )):
                    mutable_fields[key] = field
        if DEBUG_PARSE:
            logger.info('r: %r, mutable fields: %r', self._meta.resource_name, 
            mutable_fields.keys() )
        
        initializer_dict = {}
        for key,field in mutable_fields.items():
            if DEBUG_PARSE:
                logger.info('try field: %r: %r',
                    key, deserialized.get(key))
            alias = field.get('alias')
            
            if key in deserialized or alias in deserialized:            
                _val = deserialized.get(key,deserialized.get(alias,None))
                if _val is not None:
                    try:
                        _val = parse_val(_val,key,field['data_type']) 
                        if DEBUG_PARSE:
                            logger.info('parsed: %r, %r, %r',
                                key, deserialized.get(key,None), _val)
                    except ValidationError, e:
                        logger.info('error: %r', e)
                        errors.update(e.errors)
                
                # TODO: review the workflow intended for applying a default
                # value to a field: should this be handled by the client instead?
                # - are there any good workflows to justify this?
                # * used by the field initialization process
                if _val is None:
                    _val = field.get('default',None)
                    if _val == '': 
                        _val = None
                        
                initializer_dict[key] = _val
        
        if errors:
            raise ValidationError(errors)
        return initializer_dict
    
    @read_authorization
    def search(self, request, **kwargs):
        '''
        Treats a POST request like a GET (with posted raw_search_data):
        - for large raw_search_data (too large for URL param encoding)
        '''
         
        DEBUG_SEARCH = False or logger.isEnabledFor(logging.DEBUG)
         
        if SCHEMA.API_PARAM_COMPLEX_SEARCH_ID not in kwargs:
            raise MissingParam(SCHEMA.API_PARAM_COMPLEX_SEARCH_ID)
        search_ID = kwargs[SCHEMA.API_PARAM_COMPLEX_SEARCH_ID]
         
        all_params = self._convert_request_to_dict(request)
        all_params.update(kwargs)

        raw_search_data = all_params.get(SCHEMA.API_PARAM_SEARCH, None)
         
        if raw_search_data:
            raw_search_data = urllib.unquote(raw_search_data)
            # cache the search data on the session, to support subsequent requests
            # to download or modify
            request.session[search_ID] = raw_search_data  
        else:
            if search_ID in request.session:
                raw_search_data = request.session[search_ID]
            else:
                raise BadRequestError({
                    SCHEMA.API_PARAM_SEARCH: 
                    'required for search: %r, %s/%s' % (search_ID, 
                        self._meta.resource_name, SCHEMA.URI_PATH_COMPLEX_SEARCH)})
        
        if DEBUG_SEARCH:
            logger.info('complex search: %: %r', SCHEMA.API_PARAM_SEARCH, raw_search_data)
        kwargs[SCHEMA.API_PARAM_SEARCH] = raw_search_data
 
        return self.get_list(request,**kwargs)
    
    def _parse_list_ids(self, deserialized, schema):
        id_query_params = defaultdict(set)
        # store ids by row for ValidationError key
        rows_to_ids = defaultdict(dict)
        for _row,_data in enumerate(deserialized):
            try:
                id_kwargs = self.get_id(_data, schema=schema)
            except ValidationError as e:
                # Consider CumulativeError
                e.errors['input_row'] = _row
                raise
            logger.debug('found id_kwargs: %r from %r', id_kwargs, _data)
            if id_kwargs:
                rows_to_ids[_row] = id_kwargs
                for idkey,idval in id_kwargs.items():
                    id_param = '%s__in' % idkey
                    id_query_params[id_param].add(idval)
        return (id_query_params,rows_to_ids)
    
    @write_authorization
    @un_cache
    @transaction.atomic   
    def patch_list(self, request, **kwargs):

        logger.info('patch list, user: %r, resource: %r' 
            % ( request.user.username, self._meta.resource_name))
        logger.debug('patch list: %r' % kwargs)

        deserialize_meta = None
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized, deserialize_meta = self.deserialize(
                request, format=kwargs.get('format', None))

        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
        
        # Unspool the generator in memory
        # Note: Memory performance optimizations should override patch_list;
        # e.g. see WellResource.patch_list
        if isinstance(deserialized, dict):
            deserialized = [deserialized]
        else:
            deserialized = [x for x in deserialized]
            
        if len(deserialized) == 0:
            meta = { 
                SCHEMA.API_MSG_RESULT: {
                    SCHEMA.API_MSG_SUBMIT_COUNT : 0, 
                    SCHEMA.API_MSG_UPDATED: 0, 
                    SCHEMA.API_MSG_CREATED: 0,
                    SCHEMA.API_MSG_UNCHANGED: 0, 
                    SCHEMA.API_MSG_COMMENTS: 'no data patched'
                }
            }
            logger.info('PATCH list: %r', meta)
            return self.build_response(
                request, { API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)

        if len(deserialized) == 1 or isinstance(deserialized, dict):
            # send to patch detail to bypass parent log creation
            if not isinstance(deserialized,dict):
                kwargs['data'] = deserialized[0]
            else:
                kwargs['data'] = deserialized
            return self.patch_detail(request, **kwargs)
        
        logger.debug(
            'Limit the potential candidates for logging to found id_kwargs...')
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        kwargs_for_log = kwargs.copy()
        kwargs_for_log['visibilities'] = ['d','l']
        kwargs_for_log['schema'] = schema

        # 20180227 - set visibilities to detail and list to make up for removing
        # includes='*' from get_list_internal
        includes = set()
        for _data in deserialized:
            includes |= set(_data.keys())
        kwargs_for_log['includes'] = list(includes)
            
        (id_query_params,rows_to_ids) = self._parse_list_ids(deserialized, schema)
        if not id_query_params:
            logger.info('No ids found for PATCH (may be ok if id is generated)')
        kwargs_for_log.update(id_query_params)

        try:
            logger.debug('get original state, for logging... %r', 
                { k:v for k,v in kwargs_for_log.items() if k != 'schema'} )
            original_data = self._get_list_response_internal(**kwargs_for_log)
            logger.info('original state retrieved: %d', len(original_data))
        except Exception as e:
            logger.exception('original state not obtained')
            original_data = []

        if 'parent_log' not in kwargs:
            parent_log = self.make_log(request, schema=schema)
            parent_log.key = self._meta.resource_name
            parent_log.uri = self._meta.resource_name
            parent_log.save()
            kwargs['parent_log'] = parent_log
        parent_log = kwargs['parent_log']    
        logger.info('perform patch_list: %d', len(deserialized))
        for _row,_dict in enumerate(deserialized):
            try:
                self.patch_obj(request, _dict, **kwargs)
            except ValidationError, e:
                # TODO: consider CumulativeError
                e.errors['input_row'] = _row
                e.errors.setdefault('input_id', rows_to_ids[_row])
                raise
        logger.debug('Get new state, for logging: %r...',
            {k:v for k,v in kwargs_for_log.items() if k != 'schema'})
        new_data = self._get_list_response_internal(**kwargs_for_log)
        logger.info('new data: %d, log patches...', len(new_data))
        
        logs = self.log_patches(request, original_data,new_data,schema=schema,**kwargs)
        logger.info('patch logs created: %d', len(logs) if logs else 0 )
        
        patch_count = len(deserialized)
        update_count = len([x for x in logs if x.diffs ])
        logger.debug('updates: %r', [x for x in logs if x.diffs ])
        create_count = len([x for x in logs if x.api_action == API_ACTION.CREATE])
        unchanged_count = patch_count - update_count - create_count
        meta = { 
            SCHEMA.API_MSG_RESULT: {
                SCHEMA.API_MSG_SUBMIT_COUNT : patch_count, 
                SCHEMA.API_MSG_UPDATED: update_count, 
                SCHEMA.API_MSG_CREATED: create_count,
                SCHEMA.API_MSG_UNCHANGED: unchanged_count, 
                SCHEMA.API_MSG_COMMENTS: parent_log.comment
            }
        }
        if deserialize_meta:
            meta.update(deserialize_meta)
        
        param_hash = self._convert_request_to_dict(request)
        if 'test_only' in param_hash:
            logger.info('test_only flag: %r', kwargs.get('test_only'))    
            raise InformationError({
                'test_only': 'successful patch, "test_only" flag is set, rollback...',
                API_RESULT_META: meta 
            })

        if not self._meta.always_return_data:
            return self.build_response(
                request, { API_RESULT_META: meta }, response_class=HttpResponse)
        else:
            kwargs_for_log['limit'] = 0
            logger.debug(
                'return data with post response: %r, kwargs: %r', 
                meta, kwargs_for_log)
            response = self.get_list(request, meta=meta, **kwargs_for_log)             
            response.status_code = 200
            return response
 
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def put_list(self,request, **kwargs):

        # TODO: enforce a policy that either objects are patched or deleted
        # and then posted/patched
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        logger.debug('resource: %r, schema.key: %r', 
            self._meta.resource_name, schema['key'])
        logger.info('put list, user: %r, resource: %r' 
            % ( request.user.username, self._meta.resource_name))
        logger.debug('schema field keys: %r', schema[RESOURCE.FIELDS].keys())
        deserialize_meta = None
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized, deserialize_meta = self.deserialize(
                request, format=kwargs.get('format', None))

        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]

        # Unspool the generator in memory
        # Note: Memory performance optimizations should override put_list;
        # e.g. see WellResource.patch_list
        if isinstance(deserialized, dict):
            deserialized = [deserialized]
        else:
            deserialized = [x for x in deserialized]
            
        kwargs_for_log = kwargs.copy()
        kwargs_for_log['schema'] = schema
        # 20180227 - set visibilities to detail and list to make up for removing
        # includes='*' from get_list_internal
        includes = set()
        for _data in deserialized:
            includes |= set(_data.keys())
        kwargs_for_log['includes'] = list(includes)

        (id_query_params,rows_to_ids) = self._parse_list_ids(deserialized, schema)

        if not id_query_params:
            logger.info('No ids found for PATCH (may be ok if id is generated)')
        kwargs_for_log.update(id_query_params)

        try:
            logger.info('get original state, for logging...')
            logger.debug('kwargs_for_log: %r', kwargs_for_log)
            original_data = self._get_list_response_internal(**kwargs_for_log)
        except Exception as e:
            logger.exception('original state not obtained')
            original_data = []

        logger.info('put list kwargs_for_log: %r', 
            {k:v for k,v in kwargs_for_log.items() if k != 'schema' })
        logger.debug('put list %s, %s',deserialized,kwargs)

        # TODO: review REST actions:
        # PUT deletes the endpoint

        if 'parent_log' not in kwargs:
            parent_log = self.make_log(request, schema=schema)
            parent_log.key = self._meta.resource_name
            parent_log.uri = self._meta.resource_name
            parent_log.save()
            kwargs['parent_log'] = parent_log
        
        logger.info('PUT: delete queryset: %r, %r', 
            request.user.username, self._meta.resource_name)
        self.delete_list(request, **kwargs);
        new_objs = []
        logger.info('PUT: create new objs: %d', len(deserialized))

        for _row,_dict in enumerate(deserialized):
            try:
                new_objs.append(self.put_obj(request, _dict, **kwargs))
            except ValidationError, e:
                # TODO: consider CumulativeError
                e.errors['input_row'] = _row
                e.errors.setdefault('input_id', rows_to_ids[_row])
                raise

        # Get new state, for logging
        # After patch, the id keys must be present
        id_attribute = resource = schema['id_attribute']
        for idkey in id_attribute:
            id_param = '%s__in' % idkey
            ids = set(kwargs_for_log.get(id_param,[]))
            for new_obj in new_objs:
                if hasattr(new_obj, idkey):
                    idval = getattr(new_obj, idkey)
                    ids.add(idval)
            kwargs_for_log[id_param] = ids
        try:
            logger.info('get new state, for logging...')
            logger.debug('kwargs_for_log: %r', kwargs_for_log)
            new_data = self._get_list_response_internal(**kwargs_for_log)
        except Exception as e:
            logger.exception('original state not obtained')
            new_data = []
        
        logger.debug('put list done, new data: %d', len(new_data))
        #self.log_patches(request, original_data,new_data,schema=schema,**kwargs)
        logs = self.log_patches(request, original_data,new_data,schema=schema,**kwargs)
        logger.info('put logs created: %d', len(logs) if logs else 0 )
        put_count = len(deserialized)
        update_count = len([x for x in logs if x.diffs ])
        delete_count = len([x for x in logs if x.api_action == API_ACTION.DELETE])
        meta = { 
            SCHEMA.API_MSG_RESULT: {
                SCHEMA.API_MSG_SUBMIT_COUNT : put_count, 
                SCHEMA.API_MSG_UPDATED: update_count,
                'Deleted': delete_count,
                SCHEMA.API_MSG_COMMENTS: parent_log.comment
            }
        }
        if deserialize_meta:
            meta.update(deserialize_meta)
        
        param_hash = self._convert_request_to_dict(request)
        if 'test_only' in param_hash:
            logger.info('test_only flag: %r', kwargs.get('test_only'))    
            raise InformationError({
                'test_only': 'successful patch, "test_only" flag is set, rollback...',
                API_RESULT_META: meta 
            })
        
        logger.info('put_list done.')
        if not self._meta.always_return_data:
            logger.info('put success, no data')
            return HttpResponse(status=202)
        else:
            kwargs_for_log['limit'] = 0
            response = self.get_list(request, meta=meta, **kwargs_for_log)             
            response.status_code = 200
            return response 

    @write_authorization
    @un_cache        
    @transaction.atomic
    def post_list(self, request, schema=None, **kwargs):
        '''
        POST is used to create or update a resource; not idempotent;
        - The LIMS client will use POST LIST to create and update because
        Django will only "attach" files with POST (not PATCH)
        '''
        logger.info('post list, user: %r, resource: %r' 
            % ( request.user.username, self._meta.resource_name))
        logger.debug('post list: %r' % kwargs)

        deserialize_meta = None
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized, deserialize_meta = self.deserialize(
                request, format=kwargs.get('format', None))

        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]

        # Unspool the generator in memory
        # Note: Memory performance optimizations should override post_list;
        # e.g. see WellResource.patch_list
        if isinstance(deserialized, dict):
            deserialized = [deserialized]
        else:
            deserialized = [x for x in deserialized]
            
        if len(deserialized) == 0:
            meta = { 
                SCHEMA.API_MSG_RESULT: {
                    SCHEMA.API_MSG_SUBMIT_COUNT : 0, 
                    SCHEMA.API_MSG_UPDATED: 0, 
                    SCHEMA.API_MSG_CREATED: 0,
                    SCHEMA.API_MSG_UNCHANGED: 0, 
                    SCHEMA.API_MSG_COMMENTS: 'no data posted'
                }
            }
            return self.build_response(
                request, { API_RESULT_META: meta }, response_class=HttpResponse, **kwargs)
        
        # Post list may actually be a post_detail
        if len(deserialized) == 1 or isinstance(deserialized, dict):
            # send to post detail to bypass parent log creation
            if not isinstance(deserialized,dict):
                kwargs['data'] = deserialized[0]
            else:
                kwargs['data'] = deserialized
            return self.post_detail(request, **kwargs)

        if schema is None:
            raise Exception('schema not initialized')
        id_attribute = schema['id_attribute']
        
        # Limit the potential candidates for logging to found id_kwargs
        kwargs_for_log = kwargs.copy()
        kwargs_for_log['schema'] = schema
        kwargs_for_log['visibilities'] = ['d','l']
        
        # 20180227 - set visibilities to detail and list to make up for removing
        # includes='*' from get_list_internal
        includes = set()
        for _data in deserialized:
            includes |= set(_data.keys())
        kwargs_for_log['includes'] = list(includes)
            
        
        (id_query_params,rows_to_ids) = self._parse_list_ids(deserialized, schema)
        if not id_query_params:
            logger.info('No ids found for PATCH (may be ok if id is generated)')
        kwargs_for_log.update(id_query_params)
        
        try:
            logger.debug('get original state, for logging...')
            logger.debug('kwargs_for_log: %r', kwargs_for_log)
            original_data = self._get_list_response_internal(**kwargs_for_log)
        except Exception as e:
            logger.exception('original state not obtained')
            original_data = []

        if 'parent_log' not in kwargs:
            parent_log = self.make_log(request, schema=schema)
            parent_log.key = self._meta.resource_name
            parent_log.uri = self._meta.resource_name
            parent_log.save()
            kwargs['parent_log'] = parent_log
        
        new_objs = []
        for _row,_dict in enumerate(deserialized):
            try:
                new_objs.append(self.patch_obj(request, _dict, **kwargs))
            except ValidationError, e:
                # TODO: consider CumulativeError
                e.errors['input_row'] = _row
                e.errors.setdefault('input_id', rows_to_ids[_row])
                raise

        # Get new state, for logging
        # After patch, the id keys must be present
        for idkey in id_attribute:
            id_param = '%s__in' % idkey
            extant_ids = set(kwargs_for_log.get(id_param,[]))
            for new_obj in new_objs:
                if hasattr(new_obj, idkey):
                    idval = getattr(new_obj, idkey)
                    extant_ids.add(idval)
            kwargs_for_log[id_param] = extant_ids
        new_data = self._get_list_response_internal(**kwargs_for_log)
        logger.debug('post list done, new data: %d', len(new_data))

        logs = self.log_patches(request, original_data,new_data,schema=schema, **kwargs)
        logger.info('post logs created: %d', len(logs) if logs else 0 )
        patch_count = len(deserialized)
        update_count = len([x for x in logs if x.diffs ])
        create_count = len([x for x in logs if x.api_action == API_ACTION.CREATE])
        unchanged_count = patch_count - update_count - create_count
        meta = { 
            SCHEMA.API_MSG_RESULT: {
                SCHEMA.API_MSG_SUBMIT_COUNT: patch_count, 
                SCHEMA.API_MSG_UPDATED: update_count, 
                SCHEMA.API_MSG_CREATED: create_count,
                SCHEMA.API_MSG_UNCHANGED: unchanged_count, 
                SCHEMA.API_MSG_COMMENTS: parent_log.comment
            }
        }
        if deserialize_meta:
            meta.update(deserialize_meta)
        
        logger.info('POST list: %r', meta)
        if not self._meta.always_return_data:
            return self.build_response(
                request, { API_RESULT_META: meta }, response_class=HttpResponse)
        else:
            logger.info('return data with post response')
            
            response = self.get_list(request, meta=meta, **kwargs)             
            response.status_code = 200
            return response

    @write_authorization
    @un_cache  
    @transaction.atomic      
    def post_detail(self, request, **kwargs):
        '''
        POST is used to create a resource; 
        - The LIMS client will use POST to create exclusively
        - Error if the resource exists
        '''
        deserialize_meta = None
        deserialized = kwargs.pop('data', None)
        if not deserialized:
            deserialized, deserialize_meta = self.deserialize(
                request, format=kwargs.get('format', None))
            
        _data = self.build_post_detail(request, deserialized, **kwargs)
        logger.debug('post detail %s, %r', 
            { k:v for k,v in kwargs.items() if k != 'schema'},deserialized)
        
        
        # TODO: status_code=201
        return self.build_response(
            request, _data, response_class=HttpResponse, **kwargs)
    
    def build_post_detail(self, request, deserialized, log=None, **kwargs):
        ''' Inner post detail method, returns native dict '''

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

#         if not deserialized:
#             return {}
        
        kwargs_for_log = self.get_id(deserialized,validate=False,schema=schema,**kwargs)
        
        id_attribute = schema['id_attribute']
        
        logger.info('post detail: %r, %r, %r', 
            self._meta.resource_name, kwargs_for_log, id_attribute)

        # NOTE: create a log if possible, with id_attribute, for downstream
        log = self.make_log(
            request, kwargs_for_log, id_attribute=id_attribute, schema=schema)
        original_data = None
        if kwargs_for_log and len(kwargs_for_log.items())==len(id_attribute):
            # A full id exists, query for the existing state
            try:
                original_data = self._get_detail_response_internal(**kwargs_for_log)
            except Exception, e: 
                logger.exception('exception when querying for existing obj: %s', 
                    kwargs_for_log)
            
            if original_data is not None and len(original_data) != 0:
                raise ValidationError({ 
                    k: '%r Already exists' % v for k,v in kwargs_for_log.items() })
        logger.debug('patch_obj: %r, %r', deserialized, log)
        logger.debug('patch_obj: %r', kwargs)
        patch_result = self.patch_obj(request, deserialized, log=log, **kwargs)
        if API_RESULT_OBJ in patch_result:
            obj = patch_result[API_RESULT_OBJ]
        else:
            # TODO: 20170109, legacy, convert patch_obj to return:
            # { API_RESULT_OBJ, API_RESULT_META }
            obj = patch_result
        for id_field in id_attribute:
            if id_field not in kwargs_for_log:
                val = getattr(obj, id_field,None)
                if val is not None:
                    kwargs_for_log['%s' % id_field] = val
                else:
                    logger.warn('id field: %r not found in new obj: %r', id_field, obj)
        meta = {}
        if API_RESULT_META in patch_result:
            meta = patch_result[API_RESULT_META]

        # get new state, for logging
        logger.debug('get new state, for logging, kwargs: %r', kwargs_for_log)
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        logger.debug('new post data: %r', new_data)
        if not new_data:
            raise BadRequestError({
                'POST': 'no data found for the new obj created by post: %r' % meta })
        patched_log = self.log_patch(
            request, original_data,new_data,log=log, 
            id_attribute=id_attribute, schema=schema, **kwargs)
        if patched_log:
            patched_log.save()
            logger.debug('post log: %r', patched_log)
            meta[SCHEMA.API_MSG_RESULT] = SCHEMA.API_MSG_SUCCESS
        else:
            logger.info('no post log')
            meta[SCHEMA.API_MSG_WARNING] = 'No Changes were detected'

        param_hash = self._convert_request_to_dict(request)
        if 'test_only' in param_hash:
            logger.info('test_only flag: %r', kwargs.get('test_only'))
            message = {
                'test_only': 'successful patch, "test_only" flag is set, rollback...',
                'patch_log': ApiLog.json_dumps(patched_log)
            }
            meta.update(message)
            raise InformationError({API_RESULT_META: meta})

        # 20170109 - return complex data
        new_data = { 
            API_RESULT_META: meta,
            API_RESULT_DATA: [new_data,] 
        }
        return new_data
    
    @write_authorization
    @un_cache  
    @transaction.atomic      
    def patch_detail(self, request, **kwargs):
        '''
        PATCH is used to create or update a resource; not idempotent
        '''
        
        deserialize_meta = None
        deserialized = kwargs.pop('data', None)
        if not deserialized:
            deserialized, deserialize_meta = self.deserialize(
                request, format=kwargs.get('format', None))
        _data = self.build_patch_detail(request, deserialized, meta=deserialize_meta, **kwargs)
        
        logger.debug('deserialized: %r', deserialized)
        logger.debug('data: %r', _data)    
        logger.debug('kwargs: %r', {k:v for k,v in kwargs.items() if k !='schema'})
        return self.build_response(
            request, _data, response_class=HttpResponse, **kwargs)

    def build_patch_detail(self, request, deserialized, **kwargs):
        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        logger.info('deserialized: %r', deserialized)
        
        id_attribute = schema['id_attribute']
        id_kwargs = self.get_id(
            deserialized, validate=False, schema=schema, **kwargs)

        if not id_kwargs:
            return {}
        
        # NOTE: creating a log, even if no data have changed (may be comment only)
        log = self.make_log(request, schema=schema)

        original_data = self._get_detail_response_internal(**id_kwargs)
        if not log.key:
            self.make_log_key(log, original_data, id_attribute=id_attribute,
                schema=schema, **kwargs)
        logger.debug('original data: %r', original_data)
        log.save()

        if not deserialized:
            return {}
        
        patch_result = self.patch_obj(request, deserialized, log=log, **kwargs)
        meta = {}
        if API_RESULT_META in patch_result:
            meta = patch_result[API_RESULT_META]
        if API_RESULT_OBJ in patch_result:
            obj = patch_result[API_RESULT_OBJ]
        else:
            obj = patch_result
        logger.debug('build patch detail: %r', obj)
        
        for id_field in id_attribute:
            if id_field not in id_kwargs:
                val = getattr(obj, id_field,None)
                if val is not None:
                    id_kwargs['%s' % id_field] = val
        new_data = self._get_detail_response_internal(**id_kwargs)
        logger.debug('new_data: %r', new_data)
        patched_log = self.log_patch(
            request, original_data,new_data,log=log, 
            id_attribute=id_attribute, schema=schema, **kwargs)
        if patched_log:
            patched_log.save()
            meta[SCHEMA.API_MSG_RESULT] = SCHEMA.API_MSG_SUCCESS
            logger.debug('patch log: %r', patched_log)
        else:
            logger.info('no patch log')
            meta[SCHEMA.API_MSG_WARNING] = 'No Changes were detected'
        param_hash = self._convert_request_to_dict(request)
        if 'test_only' in param_hash:
            logger.info('test_only flag: %r', kwargs.get('test_only'))
            message = {
                'test_only': 'successful patch, "test_only" flag is set, rollback...',
                'patch_log': ApiLog.json_dumps(patched_log)
            }
            meta.update(message)
            raise InformationError({ API_RESULT_META: meta })

        # 20170109 - return complex data
        new_data = { 
            API_RESULT_META: meta,
            API_RESULT_DATA: [new_data,], 
            }
            
        return new_data


    @write_authorization
    @un_cache 
    @transaction.atomic       
    def put_detail(self, request, **kwargs):
        '''
        PUT is used to create or overwrite a resource; idempotent
        Note: this version of PUT cannot be used if the resource must create
        the ID keys on create (use POST/PATCH)
        '''
                
        # TODO: enforce a policy that either objects are patched or deleted
        # and then posted/patched
        
        # TODO: if put_detail is used: rework based on post_detail
        raise ApiNotImplemented(self._meta.resource_name, 'put_detail')
    
        # schema = kwargs.pop('schema', None)
        # if not schema:
        #     raise Exception('schema not initialized')
        # 
        # if kwargs.get('data', None):
        #     # allow for internal data to be passed
        #     deserialized = kwargs['data']
        # else:
        #     deserialized = self.deserialize(
        #         request, format=kwargs.get('format', None))
        # 
        # logger.debug('put detail: %r, %r' % (deserialized,kwargs))
        #  
        # # cache state, for logging
        # # Look for id's kwargs, to limit the potential candidates for logging
        # # Note: this version of PUT cannot be used if the resource must create
        # # the ID keys on create (use POST/PATCH)
        #  
        # kwargs_for_log = self.get_id(
        #     deserialized, validate=True, schema=schema, **kwargs)
        # id_attribute = schema['id_attribute']
        # # kwargs_for_log = {}
        # # for id_field in id_attribute:
        # # if deserialized.get(id_field,None):
        # #     kwargs_for_log[id_field] = deserialized[id_field]
        # # elif kwargs.get(id_field,None):
        # #     kwargs_for_log[id_field] = kwargs[id_field]
        # logger.debug('put detail: %s, %s' %(deserialized,kwargs_for_log))
        # try:
        #     logger.info('get original state, for logging...')
        #     logger.debug('kwargs_for_log: %r', kwargs_for_log)
        #     original_data = self._get_list_response_internal(**kwargs_for_log)
        # except Exception as e:
        #     logger.exception('original state not obtained')
        #     original_data = []
        #  
        # try:
        #     logger.debug('call put_obj')
        #     obj = self.put_obj(request, deserialized, **kwargs)
        # except ValidationError as e:
        #     logger.exception('Validation error: %r', e)
        #     raise e
        #          
        # if not kwargs_for_log:
        #     for id_field in id_attribute:
        #         val = getattr(obj, id_field,None)
        #         kwargs_for_log['%s' % id_field] = val
        # 
        # # get new state, for logging
        # new_data = self._get_list_response_internal(**kwargs_for_log)
        # self.log_patches(request, original_data,new_data,**kwargs)
        #  
        # # TODO: add "Result" data to meta section, see patch_list
        #  
        # if not self._meta.always_return_data:
        #     logger.info('put success, no response data')
        #     return HttpResponse(status=200)
        # else:
        #     response.status_code = 200
        #     return response

    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_list(self, request, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'delete_list')

    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_detail(self, request, schema=None, **kwargs):

        logger.debug('delete_detail: %s,  %s' 
            % (self._meta.resource_name, kwargs))

        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        if schema is None:
            raise Exception('schema not initialized')
        id_attribute = schema['id_attribute']
        kwargs_for_log = {}
        for id_field in id_attribute:
            if kwargs.get(id_field,None):
                kwargs_for_log[id_field] = kwargs[id_field]
        logger.debug('delete detail: %s' %(kwargs_for_log))
        if not kwargs_for_log:
            raise Exception('required id keys %s' % id_attribute)
        original_data = self._get_detail_response_internal(**kwargs_for_log)
        log = self.make_log(
            request, attributes=original_data, id_attribute=id_attribute, schema=schema)
        if 'parent_log' in kwargs:
            log.parent_log = kwargs.get('parent_log', None)
        log.api_action = API_ACTION.DELETE
        log.save()
        
        self.delete_obj(request, {}, log=log, **kwargs_for_log)

        # Log
        logger.info('deleted: %s' %kwargs_for_log)
        
        # 20170601 - no diffs for delete
        # log.diffs = { k:[v,None] for k,v in original_data.items()}
        log.save()
        logger.info('delete, api log: %r', log)

        # TODO: return meta information
        return HttpResponse(status=204)

    @write_authorization
#     @un_cache        
    @transaction.atomic    
    def put_obj(self,request, deserialized, **kwargs):
        try:
            self.delete_obj(request, deserialized, **kwargs)
        except ObjectDoesNotExist,e:
            pass 
        
        return self.patch_obj(request, deserialized, **kwargs)            

    @write_authorization
#     @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'delete_obj')
    
    @write_authorization
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'patch_obj')

    def validate(self, _dict, patch=False, schema=None, fields=None):
        '''
        Perform declarative validations according the the field schema:
        @param patch if False then check all fields (for required); not just the 
        patched fields (use if object is being created). When patching, only 
        need to check the fields that are present in the _dict
        
        @return a dict of field_key->[erors] where errors are string messages
        '''
        DEBUG_VALIDATION = False or logger.isEnabledFor(logging.DEBUG)
        if DEBUG_VALIDATION:
            logger.info('validate: %r, %r', patch, schema is not None)
        
        if fields is None:
            if schema is None:
                logger.info('build raw schema for parse')
                schema = self.build_schema()
            fields = schema[RESOURCE.FIELDS]
        
        id_attribute = None
        if schema is not None:
            id_attribute = schema['id_attribute']
        
        # do validations
        errors = {}
        
        for name, field in fields.items():
            if DEBUG_VALIDATION:
                logger.info('validate key: %r, field: %r', name,field)
            if name == 'resource_uri':
                continue
            
            keyerrors = []
            if patch:
                if name not in _dict:
                    continue
                else: 
                    if id_attribute and name in id_attribute:
                        continue
                    editability = field.get('editability',None)
                    if editability and 'u' not in editability:
                        logger.info('field: %r, %r, %r', name, editability, field)
                        errors[name] = 'cannot be changed'
                        continue
                
            value = _dict.get(name,None)
            
            if DEBUG_VALIDATION:
                logger.info('validate: %r:%r',name,value)
                
            if field.get('required', False):
                if value is None:
                     keyerrors.append('required')
                if isinstance(value, basestring):
                    if len(value.strip()) == 0:
                        keyerrors.append('required')
                        
            if not value or isinstance(value, (list, tuple)) and not value[0]:
                if keyerrors:
                    errors[name] = keyerrors            
                continue
            
            ##FIXME: some vocab fields are not choices fields
            logger.debug('field: %s, choices: %r', name, field.get('choices'))
            if 'choices' in field and field['choices']:
                if field['data_type'] != 'list':
                    # note: comparing as string
                    if str(value) not in field['choices']: 
                        keyerrors.append(
                            "'%s' is not one of %r" % (value, field['choices']))
                else:
                    for x in value:
                        # note: comparing as string
                        if str(x) not in field['choices']: 
                            keyerrors.append(
                                '%r are not members of %r' 
                                % (value, field['choices']))

            if 'regex' in field and field['regex']:
                logger.debug('name: %s, value: %s check regex: %s', 
                    name, value, field['regex'] )
                # FIXME validate regex on input
                matcher = re.compile(field['regex'])
                if field['data_type'] != 'list':
                    if not matcher.match(value):
                        msg = field.get('validation_message', None)
                        if not msg:
                            msg = ( "'%s' does not match pattern: '%s'" 
                                % (value, field['regex']))
                        keyerrors.append(msg)
                else:
                    for x in value:
                        if not matcher.match(x):
                            msg = field.get('validation_message', None)
                            if not msg:
                                msg = ( "'%s' does not match pattern: '%s'" 
                                    % (x, field['regex']))
                            keyerrors.append(msg)

            if keyerrors:
                errors[name] = keyerrors

            if DEBUG_VALIDATION:
                logger.info('validate: %r:%r - %r',name,value,keyerrors)
                
        if errors:
            logger.warn('errors in submitted data: %r, errs: %s', _dict, errors)
        return errors

    @staticmethod
    def datainterchange_title_function(field_hash, id_attribute=[]):
        logger.debug(
            'get datainterchange_title_function: id_attribute: %r, keys: %r', 
            id_attribute, field_hash.keys())
        def title_function(key):
            new_title = '%s (not updatable)' % key
            # editable_flags = set(['c','u'])
            logger.debug('key: %r, %r', key, key in field_hash)
            if key in field_hash:
                fi = field_hash[key]
                editability = fi.get('editability',[])
                if ( (editability and 'u' in editability)
                    or key in id_attribute):
                    new_title = key
            logger.debug('new title: %r', new_title)
            return new_title
        return title_function
                 
    @staticmethod
    def get_vocabularies(field_hash):
        vocabularies = {}
        for key, field in field_hash.iteritems():
            if field.get('vocabulary_scope_ref', None):
                scope = field.get('vocabulary_scope_ref')
                vocabularies[key] = \
                    VocabularyResource()._get_vocabularies_by_scope(scope)
                if not vocabularies[key]:
                    logger.warn('no vocabulary found for scope: %r, field: %r', 
                        scope, key)
        return vocabularies

    @staticmethod
    def create_siunit_rowproxy_generator(field_hash, extant_generator):
        '''
        Use schema information to convert raw decimal values to SI Unit values.
        
        e.g. convert ".010" (L) to "10 mL"
        
        Schema data found in the "display_options" of the field.
        Parameters:
        default_unit: (decimal value) values will be scaled either to the 
            "default_unit" or to the nearest SI prefix value.
        symbol: (string) the unit symbol
        multiplier: (decimal value) values are scaled by the multiplier, if given
        decimals: (integer) precision is limited to decimals, if given
        '''
        DEFAULT_PRECISION = Decimal('1e-%d'%9)
        
        siunit_default_units = {}
        for key, field in field_hash.iteritems():
            if field.get('display_type', None) == 'siunit':
                display_options = field.get('display_options', None)
                if display_options is not None:
                    # TODO: move display option parsing to the Field resource
                    try:
                        display_options = display_options.replace(r"'", '"')
                        display_options = json.loads(display_options)
                        logger.debug('key: %r, decoded display_options: %r', 
                            key, display_options)
                        
                        default_unit = display_options.get('defaultUnit')
                        if not default_unit:
                            logger.error(
                                'SIUNIT Field configuration error, '
                                'no "defaultUnit" in display options: '
                                'key: %r, scope: %r, %r', 
                                key, field['scope'], display_options)
                            continue
                        symbol = display_options.get('symbol')
                        if symbol is None:
                            logger.error(
                                'SIUNIT Field configuration error, '
                                'no "symbol" in display options: '
                                'key: %r, scope: %r, %r', 
                                key, field['scope'], display_options)
                            continue
                        # NOTE: convert to string to avoid float numeric errors
                        _dict = { 
                            'default_unit': Decimal(str(default_unit)),
                            'symbol': symbol 
                        }
                        multiplier = display_options.get('multiplier', None)
                        if multiplier:
                            _dict['multiplier'] = Decimal(str(multiplier))
                        decimals = display_options.get('decimals', None)
                        if decimals:
                            _dict['decimals'] = int(decimals)
                        
                        siunit_default_units[key] = _dict
                            
                    except Exception, e:
                        logger.exception(
                            'key: %r, scope: %r, error in display options: %r - %r', 
                            key, field['scope'], display_options, e)
                else:
                    logger.error(
                        'SIUNIT Field configuration error, no display options: '
                        'key: %r, scope: %r', key, field['scope'])
        logger.debug('siunit_default_units: %r', siunit_default_units)
        def siunit_rowproxy_generator(cursor):
            if extant_generator is not None:
                cursor = extant_generator(cursor)
            class Row:
                def __init__(self, row):
                    self.row = row
                def has_key(self, key):
                    return self.row.has_key(key)
                def keys(self):
                    return self.row.keys();
                def __getitem__(self, key):
                    if not row[key]:
                        return row[key]
                    if key in siunit_default_units:

                        options = siunit_default_units[key]

                        raw_val = row[key]
                        val = Decimal(raw_val)
                        
                        if val and 'multiplier' in options:
                            val *= options['multiplier']
                        
                        default_unit = options['default_unit']
                        if val >= default_unit:
                            return '{} {}{}'.format(
                                si_unit.convert_decimal(
                                    val,default_unit, options.get('decimals',3)),
                                si_unit.get_siunit_symbol(options['default_unit']),
                                options['symbol'])
                        else:
                            (symbol,default_unit) = si_unit.get_siunit(val)
                            
                            return '{} {}{}'.format(
                                si_unit.convert_decimal(
                                    val,default_unit, options.get('decimals',3)),
                                symbol, options['symbol'])
                    else:
                        
                        return self.row[key]
            for row in cursor:
                yield Row(row)
        return siunit_rowproxy_generator

    @staticmethod    
    def create_vocabulary_rowproxy_generator(field_hash, extant_generator=None):
        '''
        Create cursor row generator:
        - generator wraps a sqlalchemy.engine.ResultProxy (cursor)
        - yields a wrapper for sqlalchemy.engine.RowProxy on each iteration
        - the wrapper will return vocabulary titles for valid vocabulary values
        in each row[key] for the key columns that are vocabulary columns.
        - returns the regular row[key] value for other columns
        '''
        vocabularies = ApiResource.get_vocabularies(field_hash)

        class Row:
            def __init__(self, row):
                self.row = row
            def has_key(self, key):
                return self.row.has_key(key)
            def keys(self):
                return self.row.keys();
            def __getitem__(self, key):
                val = self.row[key]
                if not val:
                    return val
                if key in vocabularies:
                    val = str(val)
                    if val not in vocabularies[key]:
                        logger.error(
                            ('Unknown vocabulary:'
                             ' scope:%s key:%s val:%r, keys defined: %r'),
                            field_hash[key]['vocabulary_scope_ref'], key, 
                            val,vocabularies[key].keys() )
                        return self.row[key]
                    else:
                        return vocabularies[key][val][SCHEMA.VOCABULARY.TITLE]
                else:
                    return self.row[key]
                
                
#                 
#                 if key in vocabularies:
#                     raw_val = self.row[key]
#                     if raw_val is None or not str(raw_val):
#                         return raw_val
#                     vocab = vocabularies[key].get(raw_val,
#                         vocabularies[key].get(str(raw_val)))
#                     if vocab is None:
#                         logger.error(
#                             ('Unknown vocabulary:'
#                              ' scope:%s key:%s val:%r, keys defined: %r'),
#                             field_hash[key]['vocabulary_scope_ref'], key, 
#                             raw_val,vocabularies[key].keys() )
#                         return raw_val
#                     else:
#                         return vocab[SCHEMA.VOCABULARY.TITLE]
#                 else:
#                     return self.row[key]
        
        def vocabulary_rowproxy_generator(cursor):
            if extant_generator is not None:
                cursor = extant_generator(cursor)
            for row in cursor:
                yield Row(row)
        return vocabulary_rowproxy_generator
    
    def make_child_log(self, parent_log):
        '''
        Create an *unsaved* log using values from the parent_log. 
        Note: client code may use the result of diff to elect to save.
        
        Args:
        :parent_log : will use the action, user and time information
        '''
        
        log =ApiLog()
        log.api_action = parent_log.api_action
        log.ref_resource_name = self._meta.resource_name
        log.username = parent_log.username 
        log.user_id = parent_log.user_id 
        log.date_time = parent_log.date_time
        log.comment = parent_log.comment
        log.parent_log = parent_log
    
        return log
    
    def make_log(self, request, attributes=None, **kwargs):
        ''' 
        Create an *unsaved* log using values from the request and kwargs:
        
        Args:
        :request : the request used for the action, values used:
            :request.method
            :request.user
        '''
        logger.debug(
            'make_log: %r, %r, %r', self._meta.resource_name, attributes, kwargs)
        
        log = ApiLog()
        log.api_action = str((request.method)).upper()
        log.ref_resource_name = self._meta.resource_name
        log.username = request.user.username 
        log.user_id = request.user.id 
        log.date_time = _now()
 
        if HEADER_APILOG_COMMENT in request.META:
            log.comment = request.META[HEADER_APILOG_COMMENT]
            logger.debug('HEADER_APILOG_COMMENT in request.META: %r', log.comment)
        elif kwargs and HEADER_APILOG_COMMENT in kwargs:
            log.comment = kwargs[HEADER_APILOG_COMMENT]
            logger.debug('HEADER_APILOG_COMMENT in kwargs: %r', log.comment)
        
        if kwargs:
            for key, value in kwargs.items():
                if hasattr(log, key):
                    setattr(log, key, value)
        
        if attributes is not None:
            self.make_log_key(log, attributes, **kwargs)
        if DEBUG_PATCH_LOG:
            logger.info('created log: %r', log) 
        return log

    def make_log_key(
            self, log, attributes, id_attribute=None, schema=None, **kwargs):
        
        logger.debug('make_log_key: %r, %r, %r', 
            id_attribute, schema is not None, kwargs)
        
        if id_attribute is None:
            if schema is None:
                logger.debug('build raw schema for make_log_key')
                schema = self.build_schema()
            id_attribute = schema['id_attribute']
        log.key = '/'.join([
            str(attributes[x]) for x in id_attribute if x in attributes])
        log.uri = '/'.join([log.ref_resource_name,log.key])

        logger.debug('make_log_key: %r, %r', log.key, log.uri)
    
    def log_patch(self, request, prev_dict, new_dict, log=None, 
            id_attribute=None, excludes=None, exclude_patterns=None, 
            full_create_log=False, log_empty_diffs=True, **kwargs):
        '''
        @param full_create_log create a diff log showing initial state on create
        @param log_empty_diffs creates a log even if there are no diffs, if there 
            is a comment to record
        NOTE: the log entry is not saved in this function and must be saved by
        the calling function
        '''
        default_exclude_patterns = ['comment_array']
        if exclude_patterns is None:
            exclude_patterns = default_exclude_patterns
            
        if DEBUG_PATCH_LOG:
            logger.info(
                'prev_dict: %s, ======new_dict====: %s', 
                prev_dict, new_dict)

        if log is None:
            log = self.make_log(
                request, attributes=new_dict, id_attribute=id_attribute, **kwargs)
        if not log.comment:
            if HEADER_APILOG_COMMENT in request.META:
                log.comment = request.META[HEADER_APILOG_COMMENT]
            elif kwargs and HEADER_APILOG_COMMENT in kwargs:
                log.comment = kwargs[HEADER_APILOG_COMMENT]
        if log.comment is not None and DEBUG_PATCH_LOG is True:
            logger.info('log comment: %r', log.comment)

        
        if not log.key:
            self.make_log_key(log, new_dict, id_attribute=id_attribute,
                **kwargs)

        log.parent_log = kwargs.get('parent_log', None)
        
        if prev_dict:
            log.diffs = compare_dicts(prev_dict,new_dict, excludes,exclude_patterns)
            if not log.diffs:
                if DEBUG_PATCH_LOG:
                    logger.info('no diffs found: %r, %r' 
                        % (prev_dict,new_dict))
                else:
                    if DEBUG_PATCH_LOG:
                        logger.info('no diffs found: %r', log.uri) 
                if log_empty_diffs is not True:
                    log = None
                elif not log.comment:
                    # In this case, there is nothing to log
                    log = None
            else:
                if DEBUG_PATCH_LOG:
                    logger.info('log PATCH for %r', log.uri) 
                    logger.debug('log diffs for %r: %r', log.uri, log.diffs) 
                
        else: # creating
            log.api_action = API_ACTION.CREATE
            if full_create_log is True:
                log.diffs = { 
                    key: [None,val] for key,val in new_dict.items() 
                        if val is not None }
                if log.diffs and excludes:
                    for key in excludes:
                        if key in log.diffs:
                            del log.diffs[key]
            if DEBUG_PATCH_LOG:
                logger.info('CREATE, api log: %r: %r', log, log.diffs)
            else:
                logger.debug('CREATE: %r', log.uri) 
    
        if DEBUG_PATCH_LOG:
            logger.info('log patch done: %r', log)
        return log

    @transaction.atomic
    def log_patches(self,request, original_data, new_data, schema=None, 
        full_create_log=False, log_empty_diffs=False, **kwargs):
        '''
        log differences between dicts having the same identity in the arrays:
        @param original_data - data from before the API action
        @param new_data - data from after the API action
        - dicts have the same identity if the id_attribute keys have the same
        value.
        '''
        logs = []
        
        if DEBUG_PATCH_LOG:
            logger.info('log patches: %s' %kwargs)
            logger.info('log patches original: %s, =====new data===== %s',
                original_data,new_data)
        
        if schema is None:
            logger.info('log patches-build_schema: %r', self._meta.resource_name)
            schema = self.build_schema(user=request.user)
        id_attribute = schema['id_attribute']
        
        deleted_items = list(original_data)        
        for new_dict in new_data:
            if not new_dict:
                continue
            prev_dict = None
            for c_dict in original_data:
                if c_dict:
                    prev_dict = c_dict
                    for key in id_attribute:
                        if new_dict[key] != c_dict[key]:
                            prev_dict = None
                            break
                    if prev_dict:
                        break # found
            if prev_dict:
                # if found, then it is modified, not deleted
                if DEBUG_PATCH_LOG:
                    logger.info('remove from deleted dict %r, %r',
                        prev_dict, deleted_items)
                deleted_items.remove(prev_dict)
                
            log = self.log_patch(
                request, prev_dict, new_dict, id_attribute=id_attribute, 
                schema=schema, full_create_log=full_create_log, 
                log_empty_diffs=log_empty_diffs, **kwargs)  
            if DEBUG_PATCH_LOG:
                logger.info('patch log: %r', log)          
            if log:
                logs.append(log)   
            
        for deleted_dict in deleted_items:
            
            log = self.make_log(
                request, attributes=deleted_dict, id_attribute=id_attribute,
                schema=schema)
            if 'parent_log' in kwargs:
                log.parent_log = kwargs.get('parent_log', None)

            log.api_action = API_ACTION.DELETE
            # 20170601 - no diffs for delete
            # log.diffs = { key:[val,None] for key,val in deleted_dict.items() }
            logs.append(log)
            if DEBUG_PATCH_LOG:
                logger.info('delete, api log: %r',log)

        logger.debug('logs: %r', logs)
        logs = ApiLog.bulk_create(logs)
        return logs

                
class ApiLogResource(ApiResource):
    
    class Meta:
        queryset = ApiLog.objects.all().order_by(
            'ref_resource_name', 'username','date_time')
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        resource_name='apilog' 
        authorization= ApiLogAuthorization(resource_name)
        ordering = []
        serializer = LimsSerializer()
        excludes = []
        always_return_data = True
        max_limit = 100000
    
    def __init__(self, **kwargs):
        self.scope = 'fields.apilog'
        super(ApiLogResource,self).__init__(**kwargs)

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/clear_cache%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_clear_cache'), name="api_clear_cache"),
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)/children%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_apilog_childview'), 
                name="api_dispatch_apilog_childview"),
            url((r"^(?P<resource_name>%s)/children"
                 r"/(?P<ref_resource_name>[\w\d_.\-:]+)"
                 r"/(?P<key>[\w\d_.\-\+: \/]+)"
                 r"/(?P<date_time>[\w\d_.\-\+:]+)%s$")
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_apilog_childview2'), 
                name="api_dispatch_apilog_childview2"),
            url((r"^(?P<resource_name>%s)/(?P<ref_resource_name>[\w\d_.\-:]+)"
                 r"/(?P<key>[\w\d_.\-\+: \/]+)"
                 r"/(?P<date_time>[\w\d_.\-\+:]+)%s$")
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
        
    def dispatch_clear_cache(self, request, **kwargs):
        self.clear_cache(request, **kwargs)
        return self.build_response(request, 'ok', **kwargs)

    @read_authorization
    def get_detail(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True

        id = kwargs.get('id', None)
        if id:
            return self.build_list_response(request, **kwargs)
            
        ref_resource_name = kwargs.get('ref_resource_name', None)
        if not ref_resource_name:
            raise MissingParam('ref_resource_name')
        
        key = kwargs.get('key', None)
        if not key:
            raise MissingParam('key')
        
        date_time = kwargs.get('date_time', None)
        if not date_time:
            raise MissingParam('date_time')

        return self.build_list_response(request, **kwargs)
        
    @classmethod
    def get_resource_comment_subquery(cls, resource_name):
        bridge = get_tables()
        _apilog = bridge['reports_apilog']
        _logdiff = bridge['reports_logdiff']
        _user_cte = UserResource.get_user_cte().cte('users-'+resource_name)
        _comment_apilogs = (
            select([
                _apilog.c.date_time, 
                _apilog.c.key,
                _apilog.c.username,
                _user_cte.c.name,
                _apilog.c.comment])
            .select_from(_apilog.join(
                _user_cte,_apilog.c.username==_user_cte.c.username))
            .where(_apilog.c.ref_resource_name==resource_name)
            .where(_apilog.c.comment!=None)
            .where(not_(exists(
                select([None]).select_from(_logdiff)
                    .where(_logdiff.c.log_id==_apilog.c.id))))
            .order_by(desc(_apilog.c.date_time)))
        return _comment_apilogs
    
    @read_authorization
    def get_list(self,request,**kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])

        return self.build_list_response(request, **kwargs)
            
            
    @staticmethod
    def decode_before_after(val):
        '''
        Decode list values stored as JSON in the LogDiff.before/after fields
        - convert empty strings to None
        '''
        if val is None:
            return None
        if isinstance(val, six.string_types):
            if not val:
                return None

            # NOTE: all data for diffs are stored as string 
            # - only nested lists require decoding
            if len(val) > 2 and val[0] == '[':
                # Indicates that a json representation of the original list
                # value was stored in the diff
                
                # FIXME: logdiff was using repr to generate before/after;
                # remove the unicode specifier and fix quotes
                # - this was fixed 20180710, so this can be removed
                val = re.sub(r"'",'"',val )
                val = re.sub(r'''u(['"])''',r'\1',val )
                try:
                    # try to decode nested list values
                    val = json.loads(val)
                except Exception, e:
                    logger.exception(
                        'on json loads for list value: %r, %r, %r', 
                        val, type(val), e)
        return val


    @staticmethod    
    def create_apilog_preview_rowproxy_generator(list_preview_fields, extant_generator=None):
        '''
        Create a Row proxy that will decode list values embedded as JSON:
        - generator wraps a sqlalchemy.engine.ResultProxy (cursor)
        - yields a wrapper for sqlalchemy.engine.RowProxy on each iteration
        - the wrapper will decode list values embedded as JSON in the specified 
        list_preview_fields
        - convert empty strings to None
        - returns the regular row[key] value for other columns
        
        '''
        logger.info('create_apilog_preview_rowproxy_generator: %r', list_preview_fields)
        class Row:
            def __init__(self, row):
                self.row = row
            def has_key(self, key):
                return self.row.has_key(key)
            def keys(self):
                return self.row.keys();
            def __getitem__(self, key):
                val = self.row[key]
                if isinstance(val, six.string_types) and not val:
                    return None
                if key in list_preview_fields:
                    val = ApiLogResource.decode_before_after(val)
                return val
        
        def preview_rowproxy_generator(cursor):
            if extant_generator is not None:
                cursor = extant_generator(cursor)
            for row in cursor:
                yield Row(row)
        return preview_rowproxy_generator
    
        
    def build_list_response(self,request, **kwargs):
        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        manual_field_includes = set(param_hash.get('includes', []))
        manual_field_includes.add('id')
        if DEBUG_GET_LIST: 
            logger.info('manual_field_includes: %r', manual_field_includes)

        is_for_detail = kwargs.pop('is_for_detail', False)

        if is_for_detail:
            param_hash['offset'] = 0

        try:
            
            # general setup
          
            (filter_expression, filter_hash, readable_filter_hash) = SqlAlchemyResource.\
                build_sqlalchemy_filters(schema, param_hash=param_hash)
            filename = self._get_filename(
                readable_filter_hash, schema, is_for_detail)

            if filter_expression is None and 'parent_log_id' not in kwargs:
                raise InformationError(
                    key='Filter Data',
                    msg='Please enter a filter expression to begin')
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                                  
            order_params = param_hash.get('order_by',[])
            field_hash = self.get_visible_fields(
                schema[RESOURCE.FIELDS], filter_hash.keys(), manual_field_includes, 
                param_hash.get('visibilities'), 
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.\
                build_sqlalchemy_ordering(order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    ApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            # specific setup 
            _log = self.bridge['reports_apilog']
            _parent_log = self.bridge['reports_apilog']
            _parent_log = _parent_log.alias('parent_log')
            _logdiffs = self.bridge['reports_logdiff']
            
            base_query_tables = ['reports_apilog']
            
            custom_columns = {
                'diff_keys': (
                    select([
                        func.array_to_string(func.array_agg(
                            literal_column('diffkey')),LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([_logdiffs.c.field_key.label('diffkey')])
                        .select_from(_logdiffs)
                        .where(_logdiffs.c.log_id==text('reports_apilog.id'))
                        .order_by(_logdiffs.c.field_key).alias('inner'))
                    ),
                #  create a full ISO-8601 date format
                'parent_log_uri': literal_column(
                    "parent_log.ref_resource_name "
                    "|| '/' || parent_log.key || '/' "
                    "|| to_char(parent_log.date_time, 'YYYY-MM-DD\"T\"HH24:MI:SS.MS') "
                    "|| to_char(extract('timezone_hour' from parent_log.date_time),'S00')" 
                    "||':'" 
                    "|| to_char(extract('timezone_minute' from parent_log.date_time),'FM00')" 
                    ).label('parent_log_uri'),
                'parent_log_id': _parent_log.c.id,
                'child_logs': literal_column(
                    "(select count(*) from reports_apilog ra where ra.parent_log_id=reports_apilog.id)"
                    ).label('child_logs'),
                # NOTE: log_uri does not depend on timezone
                'log_uri': literal_column(
                    "reports_apilog.ref_resource_name "
                    "|| '/' || reports_apilog.key || '/' "
                    "|| to_char(reports_apilog.date_time, 'YYYY-MM-DD\"T\"HH24:MI:SS.MS') ")
            }
            
            if 'date_time' in filter_hash:
                # ISO 8601 only supports millisecond precision, 
                # but postgres supports microsecond
                custom_columns['date_time'] = \
                    literal_column(
                        "date_trunc('millisecond',reports_apilog.date_time)")
                # TODO: 20180530
                # Remove timzone from queries: assume all times are for the current
#                         "date_trunc('millisecond',reports_apilog.date_time) AT TIME ZONE 'UTC' ")
                        # TODO: convert to UTC - SQLAlchemy and raw SQL return different values
                        # for the following (SQLAlchemy generates a timezone, raw sql UTC)
#                         "to_char(reports_apilog.date_time, 'YYYY-MM-DD\"T\"HH24:MI:SS.MS')"
#                         "|| to_char(extract('timezone_hour' from parent_log.date_time),'S00')" 
#                         "||':'" 
#                         "|| to_char(extract('timezone_minute' from parent_log.date_time),'FM00')" )
                # timezone:
                # to_char(reports_apilog.date_time, 'YYYY-MM-DD\"T\"HH24:MI:SS.MS') 
                # - this should be applied to all times
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement

            j = join(
                _log, _parent_log, _log.c.parent_log_id == _parent_log.c.id, isouter=True )
            
            stmt = select(columns.values()).select_from(j)
            
            if 'diff_keys' in filter_hash:
                # manually filter for the diff keys, for performance:
                # FIXME: "search_data" not supported
                expression = None
                for filter_expr,value in param_hash.items():
                    (field_name, filter_type, inverted) = \
                        SqlAlchemyResource.parse_filter(filter_expr)
                    value = SqlAlchemyResource.parse_filter_value(value, filter_type)
                    if field_name == 'diff_keys':
                        expression = SqlAlchemyResource.build_filter(
                            'field_key', 'string', filter_type, inverted, value)
                        stmt = stmt.where(exists(
                            select([None]).select_from(_logdiffs)
                                .where(_logdiffs.c.log_id==_log.c.id)
                                .where(expression)))
                        break
            if 'diffs' in filter_hash:
                raise BadRequestError(key='diffs', msg='filtering is not implemented')
            
                
            # general setup
            stmt = stmt.order_by('ref_resource_name','key', 'date_time')
             
            (stmt,count_stmt) = self.wrap_statement(
                stmt,order_clauses,filter_expression )
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            
            def create_diff_generator(generator):
                bridge = self.bridge
                _apilog = bridge['reports_apilog']
                _logdiff = bridge['reports_logdiff']
                conn = get_engine().connect()
                query = (
                    select([
                        _logdiff.c.field_key,
                        _logdiff.c.before,_logdiff.c.after])
                    .select_from(_logdiff))
                
                def diff_generator(cursor):
                    if generator:
                        cursor = generator(cursor)
                    class Row:
                        def __init__(self, row):
                            self.row = row
                        def has_key(self, key):
                            if key == 'diffs': 
                                return True
                            return self.row.has_key(key)
                        def keys(self):
                            return self.row.keys();
                        def __getitem__(self, key):
                            if key == 'diffs':
                                _diffs = conn.execute(
                                    query.where(_logdiff.c.log_id==row['id']))
                                if _diffs:
                                    val = {}
                                    for x in _diffs:
                                        diffkey = x[0]
                                        before = ApiLogResource.decode_before_after(x[1])
                                        after = ApiLogResource.decode_before_after(x[2])
                                        val[diffkey] = [before, after]
                                    return val
                            elif key == 'json_field':
                                val = self.row[key]
                                if val:
                                    try:
                                        return json.loads(val)
                                    except:
                                        logger.exception('error decoding json from json_field: %r',val)
                                return val
                            else:
                                return self.row[key]
                    try:
                        for row in cursor:
                            yield Row(row)
                    finally:
                        conn.close()
                        
                return diff_generator
            
            diffs = field_hash.get('diffs')
            json_field = field_hash.get('json_field')
            if diffs or json_field:
                rowproxy_generator = create_diff_generator(rowproxy_generator)
            
#             compiled_stmt = str(stmt.compile(
#                 dialect=postgresql.dialect(),
#                 compile_kwargs={"literal_binds": True}))
#             logger.info('compiled_stmt %s', compiled_stmt)
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get(API_RESULT_META, None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get_list')
            raise e  
    
    def build_schema(self, user=None, **kwargs):
        schema = super(ApiLogResource,self).build_schema(user=user, **kwargs)
        temp = [ x.key for x in 
            MetaHash.objects.all().filter(scope='resource').distinct('key')]
        schema['extraSelectorOptions'] = { 
            'label': 'Resource', 
            'searchColumn': 'ref_resource_name', 'options': temp }
        return schema        
    
    def dispatch_apilog_childview(self, request, **kwargs):
        kwargs['parent_log_id'] = kwargs.pop('id')
        logger.info('dispatch_apilog_childview for %r', kwargs['parent_log_id'])
        return ApiLogResource().dispatch('list', request, **kwargs)    

    def dispatch_apilog_childview2(self, request, **kwargs):
        parent_log = self._get_detail_response_internal(**kwargs)
        logger.info('dispatch_apilog_childview2 for parent log: %r', parent_log)
        
        ref_resource_name = kwargs.pop('ref_resource_name')
        key = kwargs.pop('key')
        date_time = kwargs.pop('date_time')
        kwargs['parent_log_id'] = parent_log['id']
        return ApiLogResource().dispatch('list', request, **kwargs)    


class FieldResource(ApiResource):
    
    class Meta:
        
        queryset = MetaHash.objects.filter(
            scope__startswith="fields.").order_by('scope','ordinal','key')
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'field'
        authorization= AllAuthenticatedReadAuthorization(resource_name)
        ordering = []
        filtering = {} 
        serializer = LimsSerializer()
        excludes = [] 
        always_return_data = True 

    def __init__(self, **kwargs):
        super(FieldResource,self).__init__(**kwargs)
        self.resource_resource = None

    def clear_cache(self, request, **kwargs):
        logger.debug('clear_cache: FieldResource...')
        ApiResource.clear_cache(self, request, **kwargs)
        self.get_resource_resource().clear_cache(request, **kwargs)
        
    def get_resource_resource(self):
        if self.resource_resource is None:
            self.resource_resource = ResourceResource()
        return self.resource_resource    
        
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.]+)/(?P<key>[\w\d_]+)%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
    
    def build_schema(self, user=None, **kwargs):
        # start with the default schema for bootstrapping
        default_field = {
            'data_type': 'string',
            'editability': ['c','u'],
            'table': 'reports_metahash',
        }
        
        default_schema = {
            'key':  {
                'key': 'key',
                'scope': 'fields.field',
                'ordinal': 1,
                'json_field_type': '',
                'data_type': 'string',
                'editability': ['c'],
            },
            'scope':  {
                'key': 'scope',
                'scope': 'fields.field',
                'ordinal': 2,
                'json_field_type': '',
                'data_type': 'string',
                'editability': ['c'],
                
            },
            'ordinal':  {
                'key': 'ordinal',
                'scope': 'fields.field',
                'ordinal': 3,
                'json_field_type': '',
                'data_type': 'integer',
                'editability': ['c','u'],
                
            },
            'json_field_type':  {
                'key': 'json_field_type',
                'scope': 'fields.field',
                'ordinal': 4,
                'json_field_type': '',
                'data_type': 'string',
                'editability': ['c','u'],
                
            },
        }
        
        field_schema = deepcopy(
            MetaHash.objects.get_and_parse(
                scope='fields.field', field_definition_scope='fields.field',
                clear=True))
        for key,val in field_schema.items():
            for k,v in default_field.items():
                if k not in val or val.get(k)==None:
                    val[k] = v
        # do not allow the default values to be changed
        for key, val in default_schema.items():
            if key in field_schema:
                field_schema[key].update(default_schema[key])
            else:
                field_schema[key] = default_schema[key]
        
        field_schema['resource_uri'] = { 
                'key': 'resource_uri',
                'scope': 'fields.%s' % self._meta.resource_name,
                'title': 'URI',
                'description': 'URI for the record',
                'data_type': 'string',
                'table': 'None', 
                'visibility':[] 
        }
        
        # TODO: the ResourceResource should create the schema; 
        # provide one here for the bootstrap case
        schema = {
            'content_types': ['json'],
            'description': 'The fields resource',
            'id_attribute': ['scope','key'],
            'key': 'field',
            'scope': 'resource',
            'table': 'metahash',
            'title_attribute': ['scope','key'],
            'ordinal': 0,
            'resource_uri': BASE_URI +'/resource/field',
            'api_name': 'reports',
            'supertype': '',
            RESOURCE.FIELDS: field_schema,
        }
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 
            'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
        
        return schema
    
    @read_authorization
    def get_detail(self, request, **kwargs):
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self,request,**kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self,request, **kwargs):
        
        # NOTE: current implementation creates the fields hash in memory:
        # TODO: sort
        # TODO: limit
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        logger.debug('param_hash: %r', 
            { k:v for k,v in param_hash.items() if k!='schema'})
        filenames = [RESOURCE.FIELDS]
        
        # Construct filters:
        # TODO: Supporting only limited filters for the Field Resource
        scope = param_hash.get('scope', None)
        if not scope:
            scope = param_hash.get('scope__exact', None)
        if not scope:
            scope = param_hash.get('scope__eq', None)
        scopes = param_hash.get('scope__in', None)
        logger.debug('scope: %r, scope_in: %r', scope, scopes)
        key = param_hash.get('key', None)
        if not key:
            key = param_hash.get('key__exact', None)
        if not key:
            key = param_hash.get('key__eq', None)
        
        key_in = param_hash.get('key__in', None)
            
        if not scope and not scopes:
            logger.debug('get all scopes...')
            scopes = MetaHash.objects.all().filter(
                scope__icontains='fields.').values_list('scope',flat=True).distinct()
            if not scopes.exists():
                # bootstrap case
                scopes = ['fields.field',]
            logger.debug('scopes retrieved')
        else:
            if scope:
                filenames=[scope]
                scopes = [scope,]
            if scopes:
                filenames = scopes
        
        logger.debug('scopes: %r', scopes)
        fields = self._build_fields(scopes)
            
        if key_in:
            keys_in = key_in
            if not isinstance(key_in,(list,tuple,set)):
                keys_in = key_in.split(LIST_DELIMITER_URL_PARAM)
            fields = [x for x in fields if x['key'] in keys_in ]
            
        # Implement limited filtering by scope
        
        # logger.info('fields: %r', [(field['key'],field['scope']) for field in fields ])
        response_hash = None
        if scope and key:
            for field in fields:
                if field['key'] == key:
                    response_hash = field
                    break
            if not response_hash:
                logger.info('Field %s/%s not found' % (scope,key))
                raise Http404('Field %s/%s not found' % (scope,key))
        else:    

            limit = param_hash.get('limit', 25)        
            try:
                limit = int(limit)
            except Exception:
                raise BadRequestError({
                    'limit': 'Please provide a positive integer: %r' % limit})
            offset = param_hash.get('offset', 0 )
            try:
                offset = int(offset)
            except Exception:
                raise BadRequestError({
                    'offset': 'Please provide a positive integer: %r' % offset })
            if offset < 0:    
                offset = -offset
            
            if limit > len(fields):
                limit = len(fields)
            if offset > len(fields):
                offset = len(fields)
            
            meta = { 'limit': limit, 'offset': offset, 'total_count': len(fields) }
            
            if offset > 0:
                fields = fields[offset:]
            if limit > 0:
                fields = fields[:limit]
            
            if kwargs.get(API_RESULT_META, None):
                temp = kwargs[API_RESULT_META]
                logger.debug('meta found in kwargs: %r', temp)
                temp.update(meta)
                meta = temp
                logger.debug('meta: %r', meta)
            logger.debug('meta: %r', meta)
            response_hash = { 
                API_RESULT_META: meta, 
                self._meta.collection_name: fields 
            }
            logger.debug('Field resource rebuilt')
        kwargs['filename'] = '_'.join(filenames)
        logger.debug('FieldResource build response: %r, %r, %r', 
            self._meta.resource_name, request, 
            {k:v for k,v in kwargs.items() if k!='schema'})
        return self.build_response(request, response_hash, **kwargs)

    def _build_fields(self, scopes=None):
        ''' Internal callers - build the schema.fields hash
        '''
        
        if not scopes:
            scopes = MetaHash.objects.all()\
                .filter(scope__icontains='fields.')\
                .values_list('scope',flat=True).distinct()
            if not scopes.exists(): 
                # bootstrap case
                scopes = ['fields.field',]

        scopes = set(scopes)
        logger.info('build_fields for scopes: %r', scopes)

        fields = {}
        field_key = '{scope}/{key}'
        for scope in scopes:
            logger.debug('build scope: %r', scope)
            field_hash = deepcopy(
                MetaHash.objects.get_and_parse(
                    scope=scope, field_definition_scope='fields.field', 
                    clear=True))
            for kvalue, field in field_hash.items():
                
                key = field_key.format(**field)
                if key in fields:
                    # FIXME: is this exception needed? seems to mess up program
                    # flow; in particular for the "bootstrap" case of field
                    # initialization, re: fields.field/key
                    # throwing an exception forces _get_list_response (internal)
                    # return an empty response.
                    raise Exception(
                        'field key is already defined: %r: %r' % (key, field))
                fields[key] = field
            
        recursion_test = list()
        def fill_field_refs(key):
            if DEBUG_RESOURCES:
                logger.info('fill field for %r', key)
            
            if key not in fields:
                logger.debug('key: %r not found in %r', key, fields.keys())
            
            else:
                field = fields[key]
                logger.debug('field: %r', key)
                
                ref_field_key = field.get('ref', None)
                if ref_field_key:
                    if key not in recursion_test:
                        recursion_test.append(key)
                    else:
                        raise Exception(
                            'recursive field ref relationship for: %r, parents: %r'
                            % (key, recursion_test))
                    ref_field = fill_field_refs(ref_field_key)
                    logger.debug('ref_key: %r found ref field: %r', 
                                ref_field_key, ref_field)
                    if ref_field:
                        temp = deepcopy(ref_field)
                        for k,v in field.items():
                            if v is not None and v != '':
                                temp[k] = v
                        recursion_test.pop()
                        return temp
                return field
            
        for key,field in fields.items():
            field['1'] = field['scope']
            field['2'] = field['key']
             
            fields[key] = fill_field_refs(key)
# 
#             vocab_scope_ref = field.get('vocabulary_scope_ref')
#             if vocab_scope_ref:
#                 vocab_hash = VocabularyResource()._get_vocabularies_by_scope(vocab_scope_ref)
#                 
#                 if vocab_hash:
#                     field['choices'] = [
#                         key for key,vocab in vocab_hash.items() 
#                             if vocab.get('is_retired') is not True]
#                 else:
#                     logger.info('could not find vocab for %r, %r, ref: %r', 
#                         key, field.get('scope'), vocab_scope_ref)    
        
        decorated = [(x['scope'],x['ordinal'],x['key'], x) for x in fields.values()]
        decorated.sort(key=itemgetter(0,1,2))
        fields = [field for scope,ordinal,key,field in decorated]
        
        return fields
    
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_list(self, request, **kwargs):
        logger.info('delete_list for fields...')
        query = MetaHash.objects.all().filter(scope__icontains='fields.')
        if query.exists():
            query.delete()
        return HttpResponse(status=204)
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_obj(self, request, deserialized, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        logger.info('delete: %r', id_kwargs)
        MetaHash.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):
        
        if DEBUG_FIELDS:
            logger.debug('patch_obj: %r: %r', request.user, self._meta.resource_name)
            logger.info('deserialized: %r', deserialized)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema[RESOURCE.FIELDS]
        
        if not deserialized:
            return {}

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        if DEBUG_FIELDS:
            logger.info('id_kwargs: %r', id_kwargs)
        
        initializer_dict = self.parse(deserialized, schema=schema, create=True)

        if DEBUG_FIELDS:
            logger.info('initializer_dict: %r', initializer_dict)
        
        if not initializer_dict:
            return {}

        field = None
        try:
            field = MetaHash.objects.get(**id_kwargs)
            if DEBUG_FIELDS:
                logger.info('got field: %r: %r', id_kwargs, field)
            errors = self.validate(deserialized, patch=True, schema=schema)
            if errors:
                if DEBUG_FIELDS:
                    logger.info('field validation errors for: %r: %r', 
                        id_kwargs, errors)
                raise ValidationError(errors)
        except ObjectDoesNotExist, e:
            if DEBUG_FIELDS:
                logger.info(
                    'Metahash field %s does not exist, creating', id_kwargs)
            field = MetaHash(**id_kwargs)
            errors = self.validate(initializer_dict, patch=False, schema=schema)
            if errors:
                raise ValidationError(errors)

        for key,val in initializer_dict.items():
            if hasattr(field,key):
                setattr(field,key,val)
        
        if field.json_field:
            json_obj = json.loads(field.json_field)
        else:
            json_obj = {}
        
        for key,val in initializer_dict.items():
            fieldinformation = fields[key]
            if fieldinformation.get('json_field_type', None):
                json_obj[key] = parse_json_field(
                    val, key, fieldinformation['json_field_type'])
                
        field.json_field = json.dumps(json_obj)
        if DEBUG_FIELDS:
            logger.info('save: %r: as %r', 
                id_kwargs, field)
        field.save()
                
        logger.debug('patch_obj done')
        return { API_RESULT_OBJ: field }
            

class ResourceResource(ApiResource):
    
    class Meta:
        queryset = MetaHash.objects.filter(
            scope="resource").order_by('scope','ordinal','key')
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'resource'
        authorization= AllAuthenticatedReadAuthorization(resource_name)
        ordering = []
        filtering = {} 
        serializer = LimsSerializer()
        excludes = [] 
        always_return_data = True 

    def __init__(self, **kwargs):
        super(ResourceResource,self).__init__(**kwargs)
        self.field_resource = None
        
    def get_field_resource(self):
        if not self.field_resource:
            self.field_resource = FieldResource()
        return self.field_resource
    
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/app_data%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_app_data'), name="api_get_app_data"),
            url(r"^(?P<resource_name>%s)/(?P<key>[\w\d_.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
        
    def get_schema(self, request, **kwargs):
        return self.build_response(
            request, self.build_schema(user=request.user, **kwargs), **kwargs)

    def get_app_data(self, request, **kwargs):
        
        app_data = {attr:value for
            attr, value in settings.APP_PUBLIC_DATA.__dict__.iteritems()
                if '__' not in attr }
        
        return self.build_response(
            request, app_data, **kwargs)

    def clear_cache(self, request, **kwargs):
        logger.debug('clear_cache ResourceResource ..')
        ApiResource.clear_cache(self, request, **kwargs)
        caches['resource_cache'].clear()
        
    def _get_resource_schema(self,resource_key, user, **kwargs):
        ''' For internal callers
        '''
        if DEBUG_RESOURCES:
            logger.info('_get_resource_schema: %r %r...', resource_key, user)
        resources = self._build_resources_internal(user)
        
        if resource_key not in resources:
            raise BadRequestError({
                'resource': 'Resource is not initialized: %r' % resource_key })
        
        schema =  resources[resource_key]
        if DEBUG_RESOURCES:
            logger.info('schema fields: %r',
                [(field['key'],field['scope']) 
                    for field in schema[RESOURCE.FIELDS].values()])
        
        full_schema = kwargs.pop('full_schema',None)
        if full_schema:
            vocabularies = self.get_vocabularies(schema[RESOURCE.FIELDS])
            for key,field in schema[RESOURCE.FIELDS].items():
                if key in vocabularies:
                    field['vocabulary'] = vocabularies[key]
        return schema
    
    @read_authorization
    def get_list(self,request,**kwargs):
        '''
        For external callers - Dispatch GET list requests
        '''
        
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)
         
    @read_authorization
    def get_detail(self, request, **kwargs):
        '''
        For external callers - Dispatch GET detail
        '''

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
        return self.build_list_response(request, **kwargs)

    def build_schema(self, user=None, **kwargs):
        '''
        Override resource method - bootstrap the "Resource" resource schema
        '''
        logger.debug('build_schema for %r: %r', self._meta.resource_name, user)
        resource_fields = self.get_field_resource()._get_list_response_internal(
            scope='fields.resource', includes='*')
        field_hash = {}
        for field in resource_fields:
            field_hash[field['key']]=field
        # default schema for bootstrap
        resource_schema = {
            'content_types': ['json'],
            'description': 'The resource resource',
            'id_attribute': ['scope','key'],
            # 'id_attribute': ['scope','key'], # FIXME: 20170926, should be this
            'key': 'resource',
            'scope': 'resource',
            'table': 'metahash',
            'title_attribute': ['key'],
            'ordinal': 0,
            'resource_uri': BASE_URI + '/resource/resource',
            'api_name': 'reports',
            'supertype': '',
            RESOURCE.FIELDS: field_hash
        }
        
        return resource_schema
    
    def build_list_response(self,request, **kwargs):

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        logger.info('calling _build_resources...')
        resources = self._build_resources_internal(user=request.user)
        logger.info('done calling _build_resources...')
                    
        # TODO: pagination, sort, filter

        # Only filter by key and scope at this point
        key = param_hash.get('key', None)
        if key:
            if key not in resources:
                raise Http404('Resource not found: %r' % key)
            response_hash = resources[key]
        else:
            values = resources.values()
            values.sort(key=lambda resource: resource['key'])
            meta = { 'limit': 0, 'offset': 0, 'total_count': len(values) }

            if kwargs.get(API_RESULT_META, None):
                temp = kwargs[API_RESULT_META]
                temp.update(meta)
                meta = temp
            logger.debug('meta: %r', meta)
            
            response_hash = { 
                API_RESULT_META: meta, 
                self._meta.collection_name: values
            }
        
        return self.build_response(request, response_hash, **kwargs)
    
    def _filter_resource(self, schema, user=None):
        '''
        Filter resource based on user authorization
        '''
        logger.debug('filter resource %r: %r', schema['key'], user)
        usergroups = set()
        is_superuser = user is not None and user.is_superuser
            
        if is_superuser:
            return schema
        
        # Qualify user's authorization to filter fields
        if user is not None:
            usergroups = set([x.name for x in user.userprofile.get_all_groups()])
        schema = deepcopy(schema)
        fields = deepcopy(schema[RESOURCE.FIELDS])
        filtered_fields = {}
        disallowed_fields = {}
        for key, field in fields.items():
            include = True
            if key not in schema['id_attribute']:
                view_groups = field.get('view_groups',[])
                if view_groups:
                    if not set(view_groups) & usergroups:
                        include = False
                        if DEBUG_AUTHORIZATION:
                            logger.info(
                                'disallowed field: %r with view_groups: %r', 
                                key, view_groups)
                    else:
                        if DEBUG_AUTHORIZATION:
                            logger.info(
                                'allowed field: %r with view_groups: %r', 
                                key, view_groups)
            if include is True:
                filtered_fields[key] = field
            else:
                disallowed_fields[key] = field
        if disallowed_fields:
            if DEBUG_AUTHORIZATION:
                logger.info('user: %r, resource: %r, disallowed fields: %r', 
                    user, schema['key'], disallowed_fields.keys())
        schema[RESOURCE.FIELDS] = filtered_fields
        
        return schema
    
    def _build_resources_internal(self, user, use_cache=True):
        '''
        Internal callers - return the resource keyed hash, from cache if possible
        '''
        if DEBUG_RESOURCES:
            logger.info('_build_resources: %r: %r', user, use_cache)
        resources = None
        user_resources = None
        resource_cache = caches['resource_cache']
        if (use_cache and self.use_cache 
            and user is not None and user.is_superuser is not True):
            
            user_cache_key = 'resources_%s' % user.username
            user_resources = resource_cache.get(user_cache_key)
        if user_resources:
            logger.debug(
                'user resource retrieved from cache: %r', user_resources.keys())
        else:    
            logger.debug('user resources not cached, building')
            if use_cache is True and self.use_cache is True:
                resources = resource_cache.get('resources')
                
            if resources:
                if DEBUG_RESOURCES is True:
                    logger.info('building user resources from cached resources')
                    for key,r in resources.items():
                        logger.info('r: %r, fields: %r', 
                            key, r[RESOURCE.FIELDS].keys())
                    
            else:
                logger.info('rebuilding resources')
            
                resources = deepcopy(
                    MetaHash.objects.get_and_parse(
                        scope='resource', field_definition_scope='fields.resource', 
                        clear=True))
                if not resources:
                    # If there are no resources, use self to bootstrap
                    logger.info('no resources found, using default resource '
                        'to bootstrap...')
                    resource = self.build_schema(user=user)
                    resources = { resource['key']: resource }
                    
                logger.info('Get the field hash...')
                all_fields = self.get_field_resource()._build_fields()
                field_hash = {}
                # build a hash out of the fields
                for field in all_fields:
                    _fields = field_hash.get(field['scope'],{})
                    _fields[field['key']]=field
                    field_hash[field['scope']] = _fields
                
                recursion_test = list()
                def get_resource(key):
                    if key not in recursion_test:
                        recursion_test.append(key)
                    else:
                        raise Exception(
                            'recursive resource relationship for: %r, parents: %r'
                            % (key, recursion_test))
                    if DEBUG_RESOURCES:
                        logger.info('build resource for %r', key)
                    resource = resources[key]
                    logger.debug('resource: %r', key)
                    resource['1'] = resource['key']
                    resource[RESOURCE.FIELDS] = field_hash.get('fields.%s'%key, {})
                    resource['resource_uri'] = '/'.join([
                        self._meta.resource_name,resource['key']
                    ])
                    
                    for field in resource[RESOURCE.FIELDS].values():
                        if not field.get('table',None):
                            field['table'] = resource.get('table', None)
                            
                    supertype_key = resource.get('supertype')
                    if DEBUG_RESOURCES:
                        logger.info('key: %r, supertype: %r', key, supertype_key)
                    if supertype_key:
                        
                        supertype_resource = get_resource(supertype_key)
                        inherited_fields = deepcopy(supertype_resource[RESOURCE.FIELDS])
                        if DEBUG_RESOURCES:
                            logger.info('%r inherits: %r', key, inherited_fields.keys())
                        inherited_fields.update(resource[RESOURCE.FIELDS])
                        resource[RESOURCE.FIELDS] = inherited_fields

                    # Always support csv (and json)
                    content_types = set(resource.get('content_types',[]))
                    content_types.add('csv')
                    resource['content_types'] = sorted(content_types)
                    
                    # extend resource specific data
                    self.extend_specific_data(resource)
                    
                    recursion_test.pop()
                    
                    return resource

                for key in resources.keys():
                    resource = get_resource(key)
                    
                if use_cache is True and self.use_cache is True:
                    logger.info('caching resources')
                    resource_cache.set('resources', resources)
                if DEBUG_RESOURCES:
                    logger.info('all unfiltered resources: %r', 
                        resources.keys())
            if user and user.is_superuser is True:
                user_resources = resources
            else:
                if DEBUG_RESOURCES:
                    logger.info('filter resources...')        
                user_resources = {}
                for key, resource in resources.items():
                    
                    user_resources[key] = self._filter_resource(resource, user)
                    if DEBUG_RESOURCES:
                        logger.info(
                            'resource: %r, %r', key, resource[RESOURCE.FIELDS].keys())
                if (use_cache is True and self.use_cache is True
                    and user is not None):
                    user_cache_key = 'resources_%s' % user.username
                    resource_cache.set(user_cache_key, user_resources)
                
            logger.debug('return resources')
        return user_resources
 
    def extend_specific_data(self, resource):
        key = resource['key']
        
        if key == 'field':
            temp = [ x.scope for x in MetaHash.objects.all()
                .filter(scope__icontains='fields.').distinct('scope')]
            resource['extraSelectorOptions'] = { 
                'label': 'Resource', 'searchColumn': 'scope', 'options': temp }
        if key == 'vocabulary':
            temp = [ x.scope for x in Vocabulary.objects.all().distinct('scope')]
            resource['extraSelectorOptions'] = { 
                'label': 'Scope', 'searchColumn': 'scope', 'options': temp }
            
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_list(self, request, **kwargs):
        MetaHash.objects.all().filter(scope='resource').delete()
        return HttpResponse(status=204)

    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_obj(self, request, deserialized, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        if DEBUG_RESOURCES:
            logger.info('delete: %r', id_kwargs)
        MetaHash.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @transaction.atomic       
    def patch_obj(self, request, deserialized, **kwargs):
        logger.debug('patch_obj: %r:%r', request.user, self._meta.resource_name)
        logger.debug('patch_obj: %r', deserialized)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema[RESOURCE.FIELDS]
        initializer_dict = {}
        for key in fields.keys():
            if key in deserialized:
                initializer_dict[key] = parse_val(
                    deserialized.get(key,None), key,
                    fields[key].get('data_type','string')) 
        
        id_kwargs = self.get_id(deserialized,schema=schema, **kwargs)
        logger.debug('id_kwargs: %r', id_kwargs)
        try:
            field = None
            try:
                field = MetaHash.objects.get(**id_kwargs)
                errors = self.validate(deserialized, patch=True, schema=schema)
                if errors:
                    raise ValidationError(errors)
            except ObjectDoesNotExist, e:
                if DEBUG_RESOURCES:
                    logger.info('Metahash resource %s does not exist, creating',
                        id_kwargs)
                field = MetaHash(**id_kwargs)
                errors = self.validate(deserialized, patch=False, schema=schema)
                if errors:
                    raise ValidationError(errors)

            for key,val in initializer_dict.items():
                if hasattr(field,key):
                    setattr(field,key,val)
            
            if field.json_field:
                json_obj = json.loads(field.json_field)
            else:
                json_obj = {}
            
            for key,val in initializer_dict.items():
                fieldinformation = fields[key]
                if fieldinformation.get('json_field_type', None):
                    json_obj[key] = parse_json_field(
                        val, key, fieldinformation['json_field_type'])
                    
            field.json_field = json.dumps(json_obj)
            field.save()
            
            if DEBUG_RESOURCES:
                logger.info('patch_obj Resource done')
            return { API_RESULT_OBJ: field }
            
        except Exception, e:
            logger.exception('on patch detail')
            raise e  


class VocabularyResource(ApiResource):
    '''
    '''

    class Meta:
        bootstrap_fields = ['scope', 'key', 'ordinal', 'json_field']
        queryset = Vocabulary.objects.all().order_by(
            'scope', 'ordinal', 'key')
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'vocabulary'
        authorization= AllAuthenticatedReadAuthorization(resource_name)
        serializer = LimsSerializer()
        always_return_data = True

    def __init__(self, **kwargs):
        super(VocabularyResource,self).__init__(**kwargs)
        # for debugging
        self.patch_elapsedtime1 = 0
        self.patch_elapsedtime2 = 0
        self.patch_elapsedtime3 = 0
        self.patch_elapsedtime4 = 0
        self.requried_fields = set(['scope','key','ordinal','title',])

    def get_debug_times(self):
        
        logger.info('vocabulary times: %r, %r, %r, %r', 
            self.patch_elapsedtime1, self.patch_elapsedtime2, 
            self.patch_elapsedtime3, self.patch_elapsedtime4)
        self.patch_elapsedtime1 = 0
        self.patch_elapsedtime2 = 0
        self.patch_elapsedtime3 = 0
        self.patch_elapsedtime4 = 0

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name="api_get_schema"),            
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.\-:]+)/"
                 r"(?P<key>[\w\d_.\-\+:]+)%s$" ) 
                        % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]

    def build_schema(self, user=None, **kwargs):
        schema = super(VocabularyResource,self).build_schema(user=user, **kwargs)
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 
            'label': 'Vocabulary', 'searchColumn': 'scope', 'options': temp }
        return schema
    
    @read_authorization
    def get_detail(self, request, **kwargs):
        key = kwargs.get('key', None)
        if not key:
            raise MissingParam('key')
        
        scope = kwargs.get('scope', None)
        if not scope:
            raise MissingParam('scope')
        
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
        return self.build_list_response(request, **kwargs)
        
        
    def _get_list_response(self, request, key=None, **kwargs):
        ''' For internal callers '''
        
        # FIXME: check request params for caching
        vocabularies = self.get_cache().get('vocabularies')
        if not vocabularies  or not self.use_cache:
            vocabularies =  ApiResource._get_list_response(self, request, **kwargs)
            self.get_cache().set('vocabularies', vocabularies)
        
        if key:
            return [vocab for vocab in vocabularies if vocab['key']==key]
        else:
            return vocabularies

        
    @read_authorization
    def get_list(self,request,**kwargs):
 
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)
         
    def clear_cache(self, request, **kwargs):
        logger.info('clear vocabulary caches')
        super(VocabularyResource,self).clear_cache(request, **kwargs)
        caches['resource_cache'].clear()
        
    def build_list_response(self,request, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            schema = self.build_schema(request.user)
            raise Exception('schema not initialized')

        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        logger.info('VocabularyResource, build list response...')
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)

        key = param_hash.pop('key', None)
        if key:
            param_hash['key__eq'] = key

        scope = param_hash.pop('scope', None)
        if scope:
            param_hash['scope__eq'] = scope
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            if DEBUG_GET_LIST: 
                logger.info('manual_field_includes: %r', manual_field_includes)
  
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(
                readable_filter_hash, schema, is_for_detail)
            
            if DEBUG_GET_LIST: 
                logger.info('filter_fields: %r, kwargs: %r', 
                    filter_hash.keys(),kwargs)
            
            
            original_field_hash = schema[RESOURCE.FIELDS]
            # Add convenience fields "1" and "2", which aid with json viewers
            original_field_hash['1'] = {
                'key': '1',
                'scope': 'fields.vocabulary',
                'data_type': 'string',
                'json_field_type': 'convenience_field',
                'ordering': 'false',
                'visibilities': []
                }
            original_field_hash['2'] = {
                'key': '2',
                'scope': 'fields.vocabulary',
                'data_type': 'string',
                'json_field_type': 'convenience_field',
                'ordering': 'false',
                'visibilities': []
                }
            original_field_hash['resource_uri'] = {
                'key': 'resource_uri',
                'scope': 'fields.vocabulary',
                'data_type': 'string',
                'json_field_type': 'convenience_field',
                'ordering': 'false',
                'visibilities': []
                }
            fields_for_sql = { 
                    key:field for key, field in original_field_hash.items() 
                if not field.get('json_field_type',None) }
            fields_for_json = { 
                    key:field for key, field in original_field_hash.items() 
                if field.get('json_field_type',None) }
            
            order_params = param_hash.get('order_by',[])
            logger.debug('VocabularyResource, build field hash...')
            field_hash = self.get_visible_fields(
                fields_for_sql, filter_hash.keys(), manual_field_includes, 
                param_hash.get('visibilities',[]), 
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            field_hash['json_field'] = {
                'key': 'json_field',
                'scope': 'fields.vocabulary',
                'data_type': 'string',
                'table': 'reports_vocabulary',
                'field': 'json_field',
                'ordering': 'false',
                'visibilities': ['l','d']
                }
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
            
            def json_field_rowproxy_generator(cursor):
                '''
                Wrap connection cursor to fetch fields embedded in the 'json_field'
                '''
                class Row:
                    def __init__(self, row):
                        self.row = row
                        if row.has_key('json_field') and row['json_field']:
                            self.json_content = json.loads(row['json_field'])
                        else:
                            self.json_content = None
                    def has_key(self, key):
                        return (key in fields_for_json or self.row.has_key(key))
                    def keys(self):
                        return self.row.keys() + fields_for_json.keys();
                    def __getitem__(self, key):
                        if key == '1':
                            return row['scope']
                        elif key == '2':
                            return row['key']
                        elif key == 'resource_uri':
                            return '/'.join([
                                'vocabularies', row['scope'], row['key']])
                        elif key not in row:
                            if key in fields_for_json:
                                if self.json_content and key not in self.json_content:
                                    logger.debug(
                                        'key %r not found in json content %r', 
                                        key, self.json_content)
                                    return None
                                elif self.json_content:
                                    return self.json_content[key]
                                else:
                                    return None
                            else:
                                return None
                        else:
                            return self.row[key]
                for row in cursor:
                    yield Row(row)
                    
            # specific setup
            _vocab = self.bridge['reports_vocabulary']
            custom_columns = {
                'json_field' : literal_column('json_field')
                }
            base_query_tables = ['reports_vocabulary'] 
            logger.debug('VocabularyResource, build_sqlalchemy_columns...')
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )
            j = _vocab
            stmt = select(columns.values()).select_from(j)

            # general setup
            (stmt,count_stmt) = self.wrap_statement(
                stmt,order_clauses,filter_expression )
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=original_field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=json_field_rowproxy_generator,
                title_function=title_function, meta=kwargs.get(API_RESULT_META, None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get list')
            raise e  
    
    def _get_vocabularies_by_scope(self, scope):
        ''' Utility method:
        @return a dict of raw_value -> vocab
        
        Note: caches all of the vocabularies in a two level dict
        - keyed by [scope][key]
        '''
        vocabularies = self.get_cache().get('vocabularies');
        if not vocabularies:
            logger.debug('vocabularies not cached, rebuilding')
            vocabularies = {}
            kwargs = {
                'limit': '0',
                'includes': '*'
            }
            _data = self._get_list_response_internal(**kwargs)
            
            for v in _data:
                _scope = v['scope']
                if _scope not in vocabularies:
                     vocabularies[_scope] = {}
                     logger.debug('created vocab by scope: %r', _scope)
                vocabularies[_scope][v['key']] = v
                
#             # Hack: activity.type is serviceactivity.type + activity.class
#             # TODO: reinstate this if needed? 20160510
#             # if 'serviceactivity.type' in vocabularies:
#             #     vocabularies['activity.type'] = \
#             #         deepcopy(vocabularies['serviceactivity.type'])
#             if 'activity.class' in vocabularies:
#                 vocabularies['activity.type'].update(
#                     deepcopy(vocabularies['activity.class']))
            
            self.get_cache().set('vocabularies', vocabularies);
        
        def find_vocab(scope):
            
            if scope[-1] == '*':
                vocab = {}
                search_scope = scope[:-1]
                for current_scope, current_vocab in vocabularies.items():
                    if search_scope in current_scope:
                        vocab.update(current_vocab)
            else:
                return vocabularies.get(scope)
            
        found_vocab = find_vocab(scope)
#         if scope in vocabularies:
#             vocab = vocabularies[scope]
        if found_vocab:
            logger.debug(
                'found vocabulary, scope: %r, result: %r', scope, found_vocab)
            return deepcopy(found_vocab)
        else:
            logger.warn('---unknown vocabulary scope: %r, %r', scope, vocabularies.keys())
            return {}
    
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_list(self, request, **kwargs):
        Vocabulary.objects.all().delete()
        return HttpResponse(status=204)
    
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def put_list(self, request, **kwargs):
        self.suppress_errors_on_bootstrap = True
        result = super(VocabularyResource, self).put_list(request, **kwargs)
        self.suppress_errors_on_bootstrap = False
        return result
    
    @write_authorization
#     @un_cache 
    @transaction.atomic       
    def delete_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        logger.debug('delete: %r', id_kwargs)
        MetaHash.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_list(self, request, **kwargs):
        # TODO: implment custom patch_list for performance
        logger.info('VocabularyResource.patch_list ...')
        self.get_debug_times()
        result = super(VocabularyResource, self).patch_list(request, **kwargs)
        # Clear cache to make sure that resource cache is also cleared
        self.clear_cache(request)
        logger.info('VocabularyResource.patch_list done.')
        self.get_debug_times()
        return result
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def post_list(self, request, schema=None, **kwargs):
        result = ApiResource.post_list(self, request, schema=schema, **kwargs)
        # Clear cache to make sure that resource cache is also cleared
        self.clear_cache(request)
        return result
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def patch_detail(self, request, **kwargs):
        result = ApiResource.patch_detail(self, request, **kwargs)
        # Clear cache to make sure that resource cache is also cleared
        self.clear_cache(request)
        return result
    
    @write_authorization
    @un_cache
    @transaction.atomic
    def post_detail(self, request, **kwargs):
        result = ApiResource.post_detail(self, request, **kwargs)
        # Clear cache to make sure that resource cache is also cleared
        self.clear_cache(request)
        return result

    @write_authorization
    @un_cache
    @transaction.atomic
    def put_detail(self, request, **kwargs):
        result = ApiResource.put_detail(self, request, **kwargs)
        # Clear cache to make sure that resource cache is also cleared
        self.clear_cache(request)
        return result

    @write_authorization
    @un_cache
    @transaction.atomic
    def put_list(self, request, **kwargs):
        result = ApiResource.put_list(self, request, **kwargs)
        # Clear cache to make sure that resource cache is also cleared
        self.clear_cache(request)
        return result
    
    def validate(self, _dict, patch=False, schema=None):
        if (patch is False and not set(_dict.keys()) & self.requried_fields):
            raise ValidationError(key='required_fields', msg='%r' %self.requried_fields)
        return {}
    
    @write_authorization
    @transaction.atomic       
    def patch_obj(self, request, deserialized, **kwargs):
        logger.debug('patch_obj: %r: %r', request.user, self._meta.resource_name)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        logger.debug('patching: %r', deserialized)
        start_time = time.time()

        fields = schema[RESOURCE.FIELDS]
        initializer_dict = {}
        for key in fields.keys():
            if key in deserialized:
                initializer_dict[key] = parse_val(
                    deserialized.get(key,None), key,
                    fields[key].get('data_type','string')) 
        logger.debug('initializer: %r', initializer_dict)
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        self.patch_elapsedtime1 += (time.time() - start_time)
        start_time = time.time()
        
        try:
            vocab = None
            try:
                logger.debug('get: %r', id_kwargs)
                vocab = Vocabulary.objects.get(**id_kwargs)
                self.patch_elapsedtime2 += (time.time() - start_time)
                start_time = time.time()
                errors = self.validate(deserialized, patch=True, schema=schema)
                if errors:
                    raise ValidationError(errors)
            except ObjectDoesNotExist, e:
                logger.debug('Vocab %s does not exist, creating', id_kwargs)
                vocab = Vocabulary(**id_kwargs)
                errors = self.validate(deserialized, patch=False, schema=schema)
                if errors:
                    raise ValidationError(errors)

            self.patch_elapsedtime3 += (time.time() - start_time)
            start_time = time.time()

            for key,val in initializer_dict.items():
                logger.debug('key: %r, val: %r', key, val)
                if hasattr(vocab,key):
                    logger.debug('setting...')
                    setattr(vocab,key,val)

            if vocab.json_field:
                json_obj = json.loads(vocab.json_field)
            else:
                json_obj = {}
            
            for key,val in initializer_dict.items():
                fieldinformation = fields[key]
                if fieldinformation.get('json_field_type', None):
                    json_obj[key] = parse_json_field(
                        val, key, fieldinformation['json_field_type'])
            
            vocab.json_field = json.dumps(json_obj)
            logger.debug('save: %r, as %r', deserialized, vocab)
            vocab.save()
            
            self.patch_elapsedtime4 += (time.time() - start_time)
            logger.debug('patch vocab done: %r', id_kwargs)        
            return { API_RESULT_OBJ: vocab }
            
        except Exception, e:
            logger.exception('on patch detail')
            raise e  

class UserResourceAuthorization(UserGroupAuthorization):
    
    def _is_resource_authorized(
            self, user, permission_type, **kwargs):
        authorized = super(UserResourceAuthorization, self)._is_resource_authorized(
            user, permission_type, **kwargs)
        logger.debug('is UserResource authorized: user: %r: %r', user, authorized)
        if authorized is True:
            return True
        
        return user.is_active

    def is_restricted_view(self, user):
        restricted =  super(UserResourceAuthorization, self).is_restricted_view(user)
        logger.debug('is restricted: %r: %r', user, restricted)
        return restricted 
    
    def get_associated_users(self, user):
        associated_users =  [user]
        logger.info('user: %r, associated users: %r', user, associated_users)
        return associated_users
    
    def get_effective_access_level(self, user, row):
        username = username = row['username']
        
        effective_access_level = None
        
        associated_users = set([user.username for user in self.get_associated_users(user)])
        logger.info('assoc: %r', associated_users)
        if user.username == username:
            effective_access_level = 3
        elif username in associated_users:
            effective_access_level = 1
        else:
            effective_access_level = None

        if DEBUG_AUTHORIZATION:
            logger.info('user: %r, effective_access_level for %r: %r', 
                user.username, username, effective_access_level)
        return effective_access_level
                
    def filter(self, user, filter_expression):
        if self.is_restricted_view(user):
            associated_users = \
                self.get_associated_users(user)
            logger.info('create_authorized_user_filter for %r', user.username)
            logger.info('associated users %r', associated_users)
            auth_filter = column('username').in_(
                [user.username for user in associated_users])
            if filter_expression is not None:
                filter_expression = and_(filter_expression, auth_filter)
            else:
                filter_expression = auth_filter
        return filter_expression

    def get_row_property_generator(self, user, fields, extant_generator):
        '''
        Filter result properties based on authorization rules
        '''
        return self.get_access_level_property_generator(
           user, fields, extant_generator)
        
class UserResource(ApiResource):

    def __init__(self, **kwargs):
        super(UserResource,self).__init__(**kwargs)
        
        self.permission_resource = None
        self.usergroup_resource = None
        
    class Meta:
        queryset = UserProfile.objects.all().order_by('username') 
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        
        # TODO: implement field filtering: users can only see other users details
        # if in readEverythingAdmin group
        authorization = UserResourceAuthorization('user')
        ordering = []
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name = 'user'
    
    def clear_cache(self, request, **kwargs):
        ApiResource.clear_cache(self, request, **kwargs)
#         self.get_resource_resource().clear_cache()
        
    def get_permission_resource(self):
        if not self.permission_resource:
            self.permission_resource = PermissionResource()
        return self.permission_resource
    
    def get_usergroup_resource(self):
        if not self.usergroup_resource:
            self.usergroup_resource = UserGroupResource()
        return self.usergroup_resource
    
    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name="api_get_schema"),            
            
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))/groups%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_user_groupview'), 
                name="api_dispatch_user_groupview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))/permissions%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_user_permissionview'), 
                name="api_dispatch_user_permissionview"),
            ]    

    def dispatch_user_groupview(self, request, **kwargs):
        # signal to include extra column
        return UserGroupResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_permissionview(self, request, **kwargs):
        # signal to include extra column
        return PermissionResource().dispatch('list', request, **kwargs)    


    def build_schema(self, user=None, **kwargs):
        
        schema = super(UserResource,self).build_schema(user=user, **kwargs)
        try:
            if 'usergroups' in schema[RESOURCE.FIELDS]: # may be blank on initiation
                schema[RESOURCE.FIELDS]['usergroups']['choices'] = \
                    [x.name for x in UserGroup.objects.all()]
        except Exception, e:
            logger.exception('on get_schema')
            raise e  
        return schema

    def get_custom_columns(self):

        _up = self.bridge['reports_userprofile']
        _p = self.bridge['reports_permission']
        _ug = self.bridge['reports_usergroup']
        _upp = self.bridge['reports_userprofile_permissions']
        _ugu = self.bridge['reports_usergroup_users']
        
        # Create a recursive CTE to enumerate all groups/supergroups/subgroups
        group_all_supergroups = \
            UserGroupResource.recursive_supergroup_query(self.bridge)

        group_all_permissions = \
            UserGroupResource.recursive_permissions_query(
                self.bridge,group_all_supergroups)
        
        group_all_subgroups = \
            UserGroupResource.recursive_subgroups_query(
                self.bridge,group_all_supergroups)
            
        group_all_users = \
            UserGroupResource.recursive_group_all_users(
                self.bridge,group_all_subgroups)

        user_all_group_permissions = select([
            _ugu.c.userprofile_id,
            func.array_agg(distinct(_p.c.id)).label('all_permissions')]).\
            select_from(
                _ugu.join(group_all_permissions,_ugu.c.usergroup_id
                    ==group_all_permissions.c.usergroup_id)).\
            where(_p.c.id==text('any(gap.permission_ids)')).\
            group_by(_ugu.c.userprofile_id)
        user_all_group_permissions = user_all_group_permissions.cte('uagp') 
        
        user_all_permissions = user_all_group_permissions.union_all(
            select([_upp.c.userprofile_id, func.array_agg(_p.c.id)]).\
            select_from(_p.join(_upp,_upp.c.permission_id==_p.c.id)).\
            group_by(_upp.c.userprofile_id)).alias('uap')
        
        # FIXME: ICCB-L specific data sharing groups
        small_molecule_usergroups = set([
            'smDsl1MutualScreens','smDsl2MutualPositives','smDsl3SharedScreens'])
        rna_usergroups = set([
            'rnaiDsl1MutualScreens','rnaiDsl2MutualPositives','rnaiDsl3SharedScreens'])
                         
        _ugu1=_ugu.alias('ugu1')
        _ugx = _ug.alias('ugx')
        custom_columns = {
            'resource_uri': func.array_to_string(array([
                BASE_URI,'user',text('reports_userprofile.username')]),'/'),
            'permissions': 
                select([func.array_to_string(
                        func.array_agg(text('inner_perms.permission')),
                        LIST_DELIMITER_SQL_ARRAY)]).\
                select_from(
                    select([func.array_to_string(array([
                            _p.c.scope,_p.c.key,_p.c.type]),'/').label('permission')
                        ]).\
                    select_from(_p.join(_upp,_p.c.id==_upp.c.permission_id)).
                    where(text('reports_userprofile.id')==_upp.c.userprofile_id).
                    order_by(text('permission')).alias('inner_perms')),
            'usergroups': 
                select([func.array_to_string(
                        func.array_agg(text('inner_groups.name')), 
                        LIST_DELIMITER_SQL_ARRAY)]).\
                select_from(
                    select([_ugx.c.name]).
                    select_from(
                        _ugx.join(_ugu1,_ugx.c.id==_ugu1.c.usergroup_id)).
                    where(_ugu1.c.userprofile_id==text('reports_userprofile.id')).
                    order_by('name').alias('inner_groups')),
            'data_sharing_levels': 
                select([func.array_to_string(
                        func.array_agg(text('inner_groups.name')), 
                        LIST_DELIMITER_SQL_ARRAY)]).\
                select_from(
                    select([_ugx.c.name]).
                    select_from(
                        _ugx.join(_ugu1,_ugx.c.id==_ugu1.c.usergroup_id)).
                    where(_ugu1.c.userprofile_id==text('reports_userprofile.id')).
                    where(_ugx.c.name.in_(small_molecule_usergroups|rna_usergroups)).
                    order_by('name').alias('inner_groups')),
            'all_permissions':
                select([func.array_to_string(func.array_agg(
                    text('innerp.permission')),LIST_DELIMITER_SQL_ARRAY)]).\
                select_from(
                    select([func.array_to_string(array([
                            _p.c.scope,_p.c.key,_p.c.type]),'/').label('permission')
                        ]).\
                    select_from(user_all_permissions).\
                    where(and_(
                        user_all_permissions.c.userprofile_id
                            ==text('reports_userprofile.id'),
                        _p.c.id==text('any(uap.all_permissions)'))
                        ).\
                    order_by(text('permission')).alias('innerp')),
            }

        return custom_columns

    @classmethod
    def get_user_cte(cls):
        
        bridge = get_tables()
        _up = bridge['reports_userprofile']
        _au = bridge['auth_user']
        
        j = _up
        j = j.join(_au, _au.c.id == _up.c.user_id)
        user_table = (
            select([
                _au.c.username,
                _concat(_au.c.first_name, ' ', _au.c.last_name).label('name'),
                _concat(_au.c.last_name, ', ', _au.c.first_name).label('last_first'),
                _au.c.email
                ])
            .select_from(j))
        return user_table

    @read_authorization
    def get_detail(self, request, **kwargs):

        username = kwargs.get('username', None)
        if not username:
            raise MissingParam('username')

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self,request,**kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])

        return self.build_list_response(request, **kwargs)

    def build_list_response(self,request, **kwargs):

        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        username = param_hash.pop('username', None)
        if username:
            param_hash['username__eq'] = username
        groupname = param_hash.pop('groupname', None)
        if groupname:
            param_hash['usergroups__eq'] = groupname

        try:
            
            # general setup
             
            manual_field_includes = set(param_hash.get('includes', ['username']))
            
            if DEBUG_GET_LIST: 
                logger.info('manual_field_includes: %r', manual_field_includes)
  
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(
                readable_filter_hash, schema, is_for_detail)
            
            filter_expression = self._meta.authorization.filter(
                request.user, filter_expression)
            
            order_params = param_hash.get('order_by',[])
            field_hash = self.get_visible_fields(
                schema[RESOURCE.FIELDS], filter_hash.keys(), manual_field_includes, 
                param_hash.get('visibilities'), 
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_params.append('username')
            order_clauses = \
                SqlAlchemyResource.build_sqlalchemy_ordering(
                    order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    ApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            # specific setup
            custom_columns = {
                'resource_uri': func.array_to_string(array([
                    BASE_URI,'user',text('auth_user.username')]),'/'),
            }
            columns = self.build_sqlalchemy_columns(
                field_hash.values(),custom_columns=custom_columns )

            # build the query statement
            
            _au = get_tables()['auth_user']
            _up = get_tables()['reports_userprofile']

            j = _up
            j = j.join(_au,_up.c.user_id==_au.c.id, isouter=True)
            stmt = select(columns.values()).select_from(j)

            # general setup
             
            (stmt,count_stmt) = \
                self.wrap_statement(stmt,order_clauses,filter_expression )
            
            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get(API_RESULT_META, None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get_list')
            raise e  

    def build_sqlalchemy_columns(self, fields, custom_columns=None ):
        
        if not custom_columns:
            custom_columns = {}
        custom_columns.update(self.get_custom_columns())
        base_query_tables = ['auth_user','reports_userprofile'] 

        return super(UserResource,self).build_sqlalchemy_columns(
            fields,base_query_tables=base_query_tables,custom_columns=custom_columns)
        
    def get_id(self, deserialized, schema=None, **kwargs):
        
        id_kwargs = {}
        username = kwargs.get('username', None)
        if username is None:
            username = deserialized.get('username', None)
        if username:
            username = username.strip()
            if len(username) > 0:
                id_kwargs['username'] = username
                
        if not id_kwargs:
            ecommons_id = kwargs.get('ecommons_id', None)
            if ecommons_id is None:
                ecommons_id = deserialized.get('ecommons_id', None)
            if ecommons_id:
                ecommons_id = ecommons_id.strip()
                if len(ecommons_id) > 0:
                    id_kwargs['ecommons_id'] = ecommons_id
        if not id_kwargs:
            if 'resource_uri' in deserialized:
                return self.find_key_from_resource_uri(
                    deserialized['resource_uri'], schema=schema)
        
        return id_kwargs
    
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        logger.info('delete userprofile: %r', id_kwargs)
        UserProfile.objects.get(**id_kwargs).delete()
    
    @write_authorization
    @transaction.atomic       
    def patch_obj(self, request, deserialized, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
            schema = self.build_schema(request.user)
        fields = schema[RESOURCE.FIELDS]

        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        logger.info('patch user: %r', id_kwargs)

        auth_user_fields = { name:val for name,val in fields.items() 
            if val['table'] and val['table']=='auth_user'}
        userprofile_fields = { name:val for name,val in fields.items() 
            if (val['table'] and val['table']=='reports_userprofile'
                or name in ['username', 'ecommons_id'])}
        is_patch = False
        try:
            # create the auth_user
            username = id_kwargs.get('username', None)
            ecommons_id = id_kwargs.get('ecommons_id', None)
            if not username:
                logger.info('username not specified, setting username to'
                    ' ecommons_id: %s', ecommons_id)
                username = ecommons_id
            
            try:
                user = DjangoUser.objects.get(username=username)
                is_patch = True
                errors = self.validate(deserialized, patch=is_patch, schema=schema)
                if errors:
                    raise ValidationError(errors)
            except ObjectDoesNotExist, e:
                logger.info('User %s does not exist, creating', id_kwargs)
                errors = self.validate(deserialized, patch=is_patch, schema=schema)
                if errors:
                    raise ValidationError(errors)
                user = DjangoUser.objects.create_user(username=username)
                # User is not active by default
                # user.is_active = True
                logger.info('created Auth.User: %s', user)

            initializer_dict = self.parse(
                deserialized, schema=schema, create=not is_patch)

            auth_initializer_dict = {
                k:v for k,v in initializer_dict.items() if k in auth_user_fields}
            if auth_initializer_dict:
                for key,val in auth_initializer_dict.items():
                    if hasattr(user,key):
                        setattr(user,key,val)
                user.save()
                logger.info('== auth user updated, is_patch:%r : %r', 
                    is_patch, user.username)
            else:
                logger.debug('no auth_user fields to update %s', deserialized)
                
            userprofile = None
            try:
                userprofile = UserProfile.objects.get(**id_kwargs)
            except ObjectDoesNotExist, e:
                if hasattr(user, 'userprofile'):
                    raise ValueError(
                        'user already exists: %s: %s' % (user, user.userprofile))
                logger.info('Reports User %s does not exist, creating' % id_kwargs)
                userprofile = UserProfile.objects.create(**id_kwargs)
                logger.info('created UserProfile: %s', userprofile)
            
            userprofile.user = user
            userprofile.save()

            # create the reports userprofile
            userprofile_initializer_dict = {
                k:v for k,v in initializer_dict.items() if k in userprofile_fields}
            if userprofile_initializer_dict:
                logger.info('initializer dict: %r', initializer_dict)
                for key,val in userprofile_initializer_dict.items():
                    logger.debug('set: %s to %r, %s',key,val,hasattr(userprofile, key))
                    
                    if key == 'permissions':
                        # FIXME: first check if permissions have changed
                        userprofile.permissions.clear()
                        if val:
                            # NOTE: staff requirement is enforced in the db application
                            # if user.is_staff != True:
                            #     raise ValidationError(
                            #         key='permissions',
                            #         msg='May only be set for staff users')
                            
                            pr = self.get_permission_resource()
                            for p in val:
                                permission_key = ( 
                                    pr.find_key_from_resource_uri(p))
                                try:
                                    permission = Permission.objects.get(**permission_key)
                                    userprofile.permissions.add(permission)
                                except ObjectDoesNotExist, e:
                                    logger.warn(
                                        'no such permission: %r, %r, %r', 
                                        p, permission_key, initializer_dict)
                                    # if permission does not exist, create it
                                    # TODO: created through the permission resource
                                    permission = Permission.objects.create(
                                        **permission_key)
                                    permission.save()
                                    logger.info('created permission: %r', permission)
                                    userprofile.permissions.add(permission)
                                    userprofile.save()
                    elif key == 'usergroups':
                        # FIXME: first check if groups have changed
                        logger.info('patch usergroups: %r', val)
                        userprofile.usergroup_set.clear()
                        if val:
                            # NOTE: staff requirement is enforced in the db application
                            # if user.is_staff != True:
                            #     raise ValidationError(
                            #         key='usergroups',
                            #         msg='May only be set for staff users')
                            ugr = self.get_usergroup_resource()
                            for g in val:
                                usergroup_key = ugr.find_key_from_resource_uri(g)
                                logger.info('usergroup_key: %r', usergroup_key)
                                try:
                                    usergroup = UserGroup.objects.get(**usergroup_key)
                                    usergroup.users.add(userprofile)
                                    usergroup.save()
                                    logger.debug(
                                        'added user %r, %r to usergroup %r', 
                                        userprofile,userprofile.user, usergroup)
                                except ObjectDoesNotExist as e:
                                    msg = ('no such usergroup: %r, initializer: %r'
                                        % (usergroup_key, initializer_dict))
                                    logger.exception(msg)
                                    raise ValidationError(msg)
                    elif hasattr(userprofile,key):
                        setattr(userprofile,key,val)

                userprofile.save()
                logger.info('created/updated userprofile %r', user.username)
            else:
                logger.info('no reports_userprofile fields to update')

            return { API_RESULT_OBJ: userprofile }
            
        except Exception:
            logger.exception('on patch_obj')
            raise  


class UserGroupResource(ApiResource):
    
    class Meta:
        queryset = UserGroup.objects.all();        
        
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        authorization = UserGroupAuthorization('usergroup') #SuperUserAuthorization()        

        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name='usergroup' 
    
    def __init__(self, **kwargs):
        super(UserGroupResource,self).__init__(**kwargs)
    
        self.permission_resource = None
        self.user_resource = None

    def clear_cache(self, request, **kwargs):
        ApiResource.clear_cache(self, request, **kwargs)
        
    def get_permission_resource(self):
        if not self.permission_resource:
            self.permission_resource = PermissionResource()
        return self.permission_resource
    
    def get_user_resource(self):
        if not self.user_resource:
            self.user_resource = UserResource()
        return self.user_resource

    def find_name(self,deserialized, **kwargs):
        name = kwargs.get('name', None)
        if not name:
            name = deserialized.get('name', None)
        if not name and 'resource_uri' in deserialized:
            keys = self.find_key_from_resource_uri(deserialized['resource_uri'])
            name = keys.get('name', None)
        if not name:
            raise NotImplementedError('must provide a group "name" parameter')
        return name
    
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_detail(self,deserialized, **kwargs):
        deserialize_meta = None
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized, deserialize_meta = self.deserialize(
                request, format=kwargs.get('format', None))
        try:
            self.delete_obj(request, deserialized, **kwargs)
            return HttpResponse(status=204)
        except ObjectDoesNotExist,e:
            return HttpResponse(status=404)
    
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_obj(self, request, deserialized, **kwargs):
        name = self.find_name(deserialized,**kwargs)
        UserGroup.objects.get(name=name).delete()
    
    @write_authorization
    @transaction.atomic       
    def patch_obj(self, request, deserialized, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        name = self.find_name(deserialized,**kwargs)
        fields = schema[RESOURCE.FIELDS]

        group_fields = { name:val for name,val in fields.items() 
            if val['table'] and val['table']=='reports_usergroup'}
        create = False

        usergroup = None
        try:
            usergroup = UserGroup.objects.get(name=name)
            logger.info('modifying usergroup: %r', usergroup)
        except ObjectDoesNotExist, e:
            logger.info('Reports UserGroup %r does not exist, creating',name)
            usergroup = UserGroup.objects.create(name=name)
            create = True
            usergroup.save()

        initializer_dict = self.parse(deserialized, create=create, schema=schema)
        logger.info('initializer_dict: %r', initializer_dict)
        
        for key,val in initializer_dict.items():
            if key == 'permissions':
                logger.info('patching %r: %r', key, val)
                usergroup.permissions.clear()
                pr = self.get_permission_resource()
                for p in val:
                    permission_key = ( 
                        pr.find_key_from_resource_uri(p))
                    try:
                        permission = Permission.objects.get(**permission_key)
                    except ObjectDoesNotExist, e:
                        logger.warn('no such permission %r, %r, %r', 
                            p, permission_key, initializer_dict)
                        # if permission does not exist, create it
                        # TODO: should be created through the permission resource
                        permission = Permission.objects.create(**permission_key)
                        permission.save()
                    usergroup.permissions.add(permission)
                    usergroup.save()
                    logger.debug(
                        'added permission %r to group %r', 
                        permission,usergroup)
            elif key == 'users':
                logger.info('patching %r: %r', key, val)
                usergroup.users.clear()
                ur = self.get_user_resource()
                for u in val:
                    user_key = ur.find_key_from_resource_uri(u)
                    try:
                        user = UserProfile.objects.get(**user_key)
                        usergroup.users.add(user)
                        logger.debug('added user %r to group %r', 
                            user, usergroup)
                    except ObjectDoesNotExist, e:
                        logger.info('no such user: %r, %r, %r', 
                            u, user_key, initializer_dict)
                        raise Exception(
                            'group: %r, no such user: %r', initializer_dict, u)
            elif key == 'super_groups':
                logger.info('patching %r: %r', key, val)
                usergroup.super_groups.clear()
                usergroup.save()
                for ug in val:
                    ug_key = self.find_key_from_resource_uri(ug)
                    try:
                        supergroup = UserGroup.objects.get(**ug_key)
                        usergroup.super_groups.add(supergroup)
                        logger.debug('added supergroup %r to group: %r',
                            supergroup,usergroup)
                    except ObjectDoesNotExist, e:
                        logger.warn(
                            'no such supergroup: %r, initializer: %r',
                            ug_key,initializer_dict)
                        raise Exception(
                            'group: %r, no such supergroup: %r', initializer_dict, ug)
            elif key == 'sub_groups':
                logger.info('patching %r: %r', key, val)
                usergroup.sub_groups.clear()
                usergroup.save()
                for ug in val:
                    ug_key = self.find_key_from_resource_uri(ug)
                    try:
                        subgroup = UserGroup.objects.get(**ug_key)
                        subgroup.super_groups.add(usergroup)
                        subgroup.save()
                        logger.debug('added subgroup %r to group %r', 
                            subgroup, usergroup)
                    except ObjectDoesNotExist, e:
                        logger.warn(
                            'no such subgroup: %r, initializer: %r',
                            ug_key,initializer_dict)
                        raise Exception(
                            'group: %r, no such subgroup: %r', initializer_dict, ug)
            elif key in group_fields and hasattr(usergroup,key):
                setattr(usergroup,key,val)
            else:
                logger.debug(
                    'unknown attribute: %r:%r, usergroup: %r, initializer: %r', 
                    key, val, usergroup,initializer_dict)
        usergroup.save()
        return { API_RESULT_OBJ: usergroup }
            
    @staticmethod    
    def recursive_supergroup_query(bridge):
        '''
        Create a recursive CTE to enumerate all groups/supergroups.
        - For use in building sqlalchemy statements.
        columns: 
        - id: the usergroup id
        - name: the usergroup name
        - sg_ids: supergroup ids (recursive) for the usergroup
        @param bridge an instance of reports.utils.sqlalchemy_bridge.Bridge
        @return: an sqlalchemy statement
        WITH group_super_rpt as (
            WITH RECURSIVE group_supergroups(from_id, sg_ids, cycle) AS 
            (
                SELECT 
                    ugsg.from_usergroup_id,
                    array[ugsg.to_usergroup_id],
                    false as cycle
                from reports_usergroup_super_groups ugsg
                UNION ALL
                SELECT
                    sgs.from_usergroup_id,
                    sgs.to_usergroup_id || g_s.sg_ids as sg_ids,
                    sgs.from_usergroup_id = any(sg_ids)
                from reports_usergroup_super_groups sgs, group_supergroups g_s
                where sgs.to_usergroup_id=g_s.from_id
                and not cycle 
            )
            select 
                ug.id, ug.name,gs.* 
            from reports_usergroup ug 
            left join group_supergroups gs on gs.from_id=ug.id 
            order by name 
        )
        select ug1.id, ug1.name,
        (
            select array_agg(distinct(ug2.id)) 
            from reports_usergroup ug2, group_super_rpt
            where ug2.id=any(group_super_rpt.sg_ids) 
            and group_super_rpt.from_id=ug1.id) as sg_ids 
        from
        reports_usergroup ug1 order by name;
        '''
        
        #Note: using the postgres specific ARRAY and "any" operator
        try:
            _ug = bridge['reports_usergroup']
            _ugsg = bridge['reports_usergroup_super_groups']
    
            ugsg1 = _ugsg.alias('ugsg1')
            group_supergroups = (
                select([
                    ugsg1.c.from_usergroup_id.label('from_id'),
                    literal_column('array[ugsg1.to_usergroup_id]').label('sg_ids'),
                    literal_column('false').label('cycle')
                ])
                .select_from(ugsg1)
                .cte('group_supergroups',recursive=True))
            gsg_alias = group_supergroups.alias('gsg')

            _ugsg_outer = _ugsg.alias('ugsg2')
            group_all_supergroups = gsg_alias.union_all(
                select([
                    _ugsg_outer.c.from_usergroup_id,
                    func.array_append(
                        gsg_alias.c.sg_ids,_ugsg_outer.c.to_usergroup_id),
                    _ugsg_outer.c.from_usergroup_id==text('any(gsg.sg_ids)')
                    ])
                .select_from(gsg_alias)
                .where(and_(
                    _ugsg_outer.c.to_usergroup_id==gsg_alias.c.from_id,
                    gsg_alias.c.cycle==False)))
            group_all_supergroups = group_all_supergroups.alias('gsg_union')
            
            # The query so far returns each path to a supergroup as a separate 
            # row, so the following aggregates all supergroups per item
            _ug1 = _ug.alias('ug1')
            _ug2 = _ug.alias('ug2')
            group_supergroup_rpt = (
                select([
                    _ug2.c.id,
                    _ug2.c.name.label('supergroup_name'),
                    select([
                        func.array_agg(distinct(_ug1.c.id))])
                    .select_from(group_all_supergroups)
                    .where(and_(
                        _ug1.c.id==text('any(gsg_union.sg_ids)'),
                        group_all_supergroups.c.from_id==_ug2.c.id))
                    .label('sg_ids')
                ])
                .select_from(_ug2)
                .order_by(_ug2.c.name))
            return group_supergroup_rpt.cte('group_sg_rpt')
        except Exception, e:
            logger.exception('on recursive_supergroup_query construction')
            raise e  

    @staticmethod
    def recursive_permissions_query(bridge,group_all_supergroups):

        _ugp = bridge['reports_usergroup_permissions']
        
        group_all_permissions = (
            select([
                group_all_supergroups.c.id.label('usergroup_id'),
                func.array_agg(_ugp.c.permission_id).label('permission_ids')])
            .where(or_(
                _ugp.c.usergroup_id==group_all_supergroups.c.id,
                _ugp.c.usergroup_id==text('any(group_sg_rpt.sg_ids)')))
            .group_by(group_all_supergroups.c.id))
        group_all_permissions = group_all_permissions.cte('gap')
        
        return group_all_permissions
    
    @staticmethod    
    def recursive_subgroups_query(bridge, group_all_supergroups):
        _ug = bridge['reports_usergroup']
        group_all_subgroups = (
            select([
                _ug.c.id,
                select([func.array_agg(group_all_supergroups.c.id)])
                    .select_from(group_all_supergroups)
                    .where(_ug.c.id==text('any(group_sg_rpt.sg_ids)'))
                    .label('subgroup_ids')
            ]).select_from(_ug))
        
        return group_all_subgroups.cte('gasubg')

    @staticmethod
    def recursive_group_all_users(bridge,group_all_subgroups):
        _up = bridge['reports_userprofile']
        _ugu = bridge['reports_usergroup_users']
        group_all_users = (
            select([
                group_all_subgroups.c.id,
                select([func.array_agg(_up.c.id)])
                    .select_from(_up.join(_ugu,_up.c.id==_ugu.c.userprofile_id))
                    .where(or_(
                        _ugu.c.usergroup_id==text('any(gasubg.subgroup_ids)'),
                        _ugu.c.usergroup_id==
                            group_all_subgroups.c.id))
                    .label('userprofile_ids')
                ])
            .select_from(group_all_subgroups))

        return group_all_users.cte('gau')
    
    @read_authorization
    def get_detail(self, request, **kwargs):

        name = kwargs.get('name', None)
        if not name:
            raise MissingParam('name')
        
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self,request,**kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])

        return self.build_list_response(request, **kwargs)
        
    def build_list_response(self,request, **kwargs):

        DEBUG_GET_LIST = False or logger.isEnabledFor(logging.DEBUG)

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        
        name = param_hash.pop('name', None)
        if name:
            param_hash['name__eq'] = name
        username = param_hash.pop('username', None)
        if username:
            param_hash['all_users__eq'] = username
        
        sub_group_name = param_hash.pop('sub_groupname',None)
        if sub_group_name:
            param_hash['all_sub_groups__eq']=sub_group_name
        
        super_groupname = param_hash.pop('super_groupname', None)
        if super_groupname:
            param_hash['all_super_groups__eq']=super_groupname    
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            if DEBUG_GET_LIST: 
                logger.info('manual_field_includes: %r', manual_field_includes)
  
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(
                readable_filter_hash, schema, is_for_detail)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
              
            order_params = param_hash.get('order_by',[])
            field_hash = self.get_visible_fields(
                schema[RESOURCE.FIELDS], filter_hash.keys(), manual_field_includes, 
                param_hash.get('visibilities'), 
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    ApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
            
            # specific setup
            _up = self.bridge['reports_userprofile']
            _p = self.bridge['reports_permission']
            _ug = self.bridge['reports_usergroup']
            _ugu = self.bridge['reports_usergroup_users']
            _ugp = self.bridge['reports_usergroup_permissions']
            _ugsg = self.bridge['reports_usergroup_super_groups']
            base_query_tables = ['reports_usergroup'] 
            
            # Create a recursive CTE to enumerate all groups/supergroups/subgroups
            group_all_supergroups = \
                UserGroupResource.recursive_supergroup_query(self.bridge)
            group_all_permissions = \
                UserGroupResource.recursive_permissions_query(
                    self.bridge,group_all_supergroups)
            group_all_subgroups = \
                UserGroupResource.recursive_subgroups_query(
                    self.bridge,group_all_supergroups)
            group_all_users = \
                UserGroupResource.recursive_group_all_users(
                    self.bridge,group_all_subgroups)
                
            _ug1 = _ug.alias('ug1')
            _ug2 = _ug.alias('ug2')
            _ug3 = _ug.alias('ug3')
            _ug4 = _ug.alias('ug4')
            custom_columns = {
                'resource_uri': (
                    func.array_to_string(
                        array([
                            BASE_URI,'usergroup',text('reports_usergroup.name')])
                        ,'/')),
                'permissions': (
                    select([
                        func.array_to_string(
                            func.array_agg(text('innerperm.permission')),
                                LIST_DELIMITER_SQL_ARRAY)
                    ])
                    .select_from(
                        select([
                            func.array_to_string(
                                array([_p.c.scope,_p.c.key,_p.c.type]),'/')
                                    .label('permission')])
                        .select_from(
                            _p.join(_ugp,_p.c.id==_ugp.c.permission_id))
                        .where(_ugp.c.usergroup_id==text('reports_usergroup.id'))
                        .order_by(text('permission'))
                        .alias('innerperm')
                    )),
                'users': (
                    select([
                        func.array_to_string(
                            func.array_agg(text('inner1.username')),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([_up.c.username])
                        .select_from(
                            _up.join(_ugu,_up.c.id==_ugu.c.userprofile_id))
                        .where(_ugu.c.usergroup_id==text('reports_usergroup.id'))
                        .order_by('username')
                        .alias('inner1')
                    )),
                'sub_groups': (
                    select([
                        func.array_to_string(
                            func.array_agg(text('inner1.name')),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([_ug2.c.name.label('name')]).
                        select_from(
                            _ug2.join(_ugsg,_ug2.c.id==_ugsg.c.from_usergroup_id))
                        .where(_ugsg.c.to_usergroup_id==text('reports_usergroup.id'))
                        .order_by('name')
                        .alias('inner1')
                    )),
                'super_groups': (
                    select([
                        func.array_to_string(
                            func.array_agg(text('inner1.name')),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([_ug3.c.name.label('name')])
                        .select_from(
                            _ug3.join(_ugsg,_ug3.c.id==_ugsg.c.to_usergroup_id))
                        .where(
                            _ugsg.c.from_usergroup_id==text('reports_usergroup.id'))
                        .order_by('name')
                        .alias('inner1')
                    )),
                'all_permissions': (
                    select([
                        func.array_to_string(
                            func.array_agg(text('allperm.permission')),                            
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([
                            func.array_to_string(
                                array([_p.c.scope,_p.c.key,_p.c.type])
                                ,'/')
                            .label('permission'),
                            group_all_permissions.c.usergroup_id ])
                        .select_from(group_all_permissions)
                        .where(_p.c.id==text('any(gap.permission_ids)'))
                        .where(
                            group_all_permissions.c.usergroup_id==
                                text('reports_usergroup.id'))
                        .order_by(text('permission'))
                        .alias('allperm')
                    )),
                'all_super_groups': (
                    select([
                        func.array_to_string(
                            func.array_agg(text('supergroup.supergroup_name')),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        select([distinct(_ug1.c.name).label('supergroup_name')])
                        .select_from(group_all_supergroups)
                        .where(and_(
                            _ug1.c.id==text('any(group_sg_rpt.sg_ids)'),
                            group_all_supergroups.c.id
                                ==text('reports_usergroup.id')))
                        .order_by(_ug1.c.name).alias('supergroup')
                    )),
                # NOTE: for orchestra/pgsql 8.4 compatability, MUST use group/agg
                # here, see following notes for this
                'all_sub_groups': (
                    func.array_to_string(
                            func.array_agg(group_all_supergroups.c.supergroup_name),
                            LIST_DELIMITER_SQL_ARRAY)
                            ),
                # NOTE: follows also does not work on orchestra/pg 8.4 ??
                # 'all_sub_groups': (
                #     select([
                #         func.array_to_string(
                #             func.array_agg(text('subgroup.name')),
                #             LIST_DELIMITER_SQL_ARRAY)])
                #     .select_from(
                #         select([distinct(_ug4.c.name)])
                #         .select_from(
                #             _ug4.join(group_all_subgroups,
                #                 _ug4.c.id==func.any(group_all_subgroups.c.subgroup_ids)))
                #         .where(group_all_subgroups.c.id==text('reports_usergroup.id'))
                #         .order_by(_ug4.c.name)
                #         .alias('subgroup')
                #     )),
                # NOTE: following form produces no results on orchestra/pg 8.4 ??
                # 'all_sub_groups': (
                #     select([
                #         func.array_to_string(
                #             func.array_agg(text('subgroup.supergroup_name')),
                #             LIST_DELIMITER_SQL_ARRAY)])
                #     .select_from(
                #         select([distinct(group_all_supergroups.c.supergroup_name)])
                #         .select_from(group_all_supergroups)
                #         .where(text(
                #             'reports_usergroup.id=any(group_sg_rpt.sg_ids)'))
                #         .order_by(group_all_supergroups.c.supergroup_name)
                #         .alias('subgroup')
                #     )),
                'all_users': (
                    select([
                        func.array_to_string(
                            func.array_agg(text('inneruser.username')),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from( 
                        select([distinct(_up.c.username)])
                        .select_from(group_all_users)
                        .where(_up.c.id==text('any(gau.userprofile_ids)'))
                        .where(
                            group_all_users.c.id==text('reports_usergroup.id'))
                        .order_by(_up.c.username).alias('inneruser')
                    )),
                }

            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            # build the query statement
            
            j = _ug.join(group_all_supergroups,
                 _ug.c.id==func.any(group_all_supergroups.c.sg_ids),isouter=True)
            j = j.join(group_all_users,
                _ug.c.id==group_all_users.c.id, isouter=True)
            stmt = select(columns.values()).select_from(j)
            # NOTE: for orchestra/pgsql 8.4 compatability, the all_subgroups 
            # requires the use of group/agg here, see notes above for this
            stmt = stmt.group_by(_ug.c.id, _ug.c.name, _ug.c.description)
            stmt = stmt.order_by('name')
            # general setup
             
            (stmt,count_stmt) = \
                self.wrap_statement(stmt,order_clauses,filter_expression )
            
            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
                        
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get(API_RESULT_META, None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get_list')
            raise e  
        
    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('get_schema'), name="api_get_schema"),            
            
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/users%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_group_userview'), 
                name="api_dispatch_group_userview"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/permissions%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_group_permissionview'), 
                name="api_dispatch_group_permissionview"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/supergroups%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_group_supergroupview'), 
                name="api_dispatch_group_supergroupview"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/subgroups%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_group_subgroupview'), 
                name="api_dispatch_group_subgroupview"),
            ]

    def dispatch_group_userview(self, request, **kwargs):
        # signal to include extra column
        kwargs['groupname'] = kwargs.pop('name')  
        return UserResource().dispatch('list', request, **kwargs)    
    
    def dispatch_group_permissionview(self, request, **kwargs):
        # signal to include extra column
        kwargs['groupname'] = kwargs.pop('name')  
        return PermissionResource().dispatch('list', request, **kwargs)       
   
    def dispatch_group_supergroupview(self, request, **kwargs):
        # signal to include extra column
        kwargs['sub_groupname'] = kwargs.pop('name')  
        return self.dispatch('list', request, **kwargs)       

    def dispatch_group_subgroupview(self, request, **kwargs):
        # signal to include extra column
        kwargs['super_groupname'] = kwargs.pop('name')  
        return self.dispatch('list', request, **kwargs)    

class PermissionResource(ApiResource):

    class Meta:
        queryset = Permission.objects.all().order_by('scope', 'key')
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        authorization= UserGroupAuthorization('permission')         
        object_class = object
        
        ordering = []
        filtering = {}
        serializer = LimsSerializer()
        excludes = [] 
        includes = [] 
        always_return_data = True # this makes Backbone happy
        resource_name='permission' 
    
    def __init__(self, **kwargs):
        super(PermissionResource,self).__init__(**kwargs)
        
        # create all of the permissions on startup
        resources = MetaHash.objects.filter(
            Q(scope='resource')|Q(scope__contains='fields.'))
        query = Permission.objects.all()
        permissionTypes = Vocabulary.objects.all().filter(
            scope='permission.type')
        try:
            for r in resources:
                found = False
                for perm in query:
                    if perm.scope==r.scope and perm.key==r.key:
                        found = True
                if not found:
                    logger.debug('initialize permission: %r:%r'
                        % (r.scope, r.key))
                    for ptype in permissionTypes:
                        p = Permission.objects.create(
                            scope=r.scope, key=r.key, type=ptype.key)
                        p.save()
                        logger.debug('bootstrap created permission %s' % p)
        except Exception, e:
            logger.info('startup exception, %r', e)
            
    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.\-:]+)/"
                 r"(?P<key>[\w\d_.\-\+:]+)/(?P<type>[\w\d_.\-\+:]+)%s$" ) 
                        % (self._meta.resource_name, TRAILING_SLASH), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]
    
    @read_authorization
    def get_detail(self, request, **kwargs):

        if not kwargs.get('scope'):
            raise MissingParam('scope')
        if not kwargs.get('key'):
            raise MissingParam('key')
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
        return self.build_list_response(request, **kwargs)
        
    @read_authorization
    def get_list(self,request,**kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self,request, **kwargs):

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        is_for_detail = kwargs.pop('is_for_detail', False)
        scope = param_hash.pop('scope', None)
        if scope:
            param_hash['scope__eq'] = scope
        key = param_hash.pop('key', None)
        if key:
            param_hash['key__eq'] = key
        
        try:
            
            # general setup
          
            manual_field_includes = set(param_hash.get('includes', []))
            
            (filter_expression, filter_hash, readable_filter_hash) = \
                SqlAlchemyResource.build_sqlalchemy_filters(
                    schema, param_hash=param_hash)
            filename = self._get_filename(
                readable_filter_hash, schema, is_for_detail)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                  
            order_params = param_hash.get('order_by',[])
            field_hash = self.get_visible_fields(
                schema[RESOURCE.FIELDS], filter_hash.keys(), manual_field_includes, 
                param_hash.get('visibilities'), 
                exact_fields=set(param_hash.get('exact_fields', [])),
                order_params=order_params)
            order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
                order_params, field_hash)
             
            rowproxy_generator = None
            if use_vocab is True:
                rowproxy_generator = \
                    ApiResource.create_vocabulary_rowproxy_generator(field_hash)
            rowproxy_generator = \
                self._meta.authorization.get_row_property_generator(
                    request.user, field_hash, rowproxy_generator)
 
            # specific setup
            _p = self.bridge['reports_permission']
            _up = self.bridge['reports_userprofile']
            _upp = self.bridge['reports_userprofile_permissions']
            _ug = self.bridge['reports_usergroup']
            _ugp = self.bridge['reports_usergroup_permissions']
            
            custom_columns = {
                'users': (
                    select([
                        func.array_to_string(
                            func.array_agg(_up.c.username),
                            LIST_DELIMITER_SQL_ARRAY)])
                        .select_from(
                            _up.join(_upp,_up.c.id==_upp.c.userprofile_id))
                        .where(_upp.c.permission_id==_p.c.id)),
                'groups': (
                    select([
                        func.array_to_string(
                            func.array_agg(_ug.c.name),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        _ug.join(_ugp,_ug.c.id==_ugp.c.usergroup_id))
                    .where(_ugp.c.permission_id==_p.c.id)),
                'usergroups': (
                    select([
                        func.array_to_string(
                            func.array_agg(_concat('usergroup/',_ug.c.name)),
                            LIST_DELIMITER_SQL_ARRAY)])
                    .select_from(
                        _ug.join(_ugp,_ug.c.id==_ugp.c.usergroup_id))
                    .where(_ugp.c.permission_id==_p.c.id)),
                'resource_uri':
                    _concat('permission/',_p.c.scope,'/',_p.c.key,'/',_p.c.type),
                        
            }

            base_query_tables = ['reports_permission'] 
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )

            j = _p
            stmt = select(columns.values()).select_from(j)
            # general setup
             
            (stmt,count_stmt) = self.wrap_statement(
                stmt,order_clauses,filter_expression )
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get(API_RESULT_META, None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

    def delete_obj(self, request, deserialized, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'delete_obj')
    
    def patch_obj(self, request, deserialized, **kwargs):
        raise ApiNotImplemented(self._meta.resource_name, 'patch_obj')
    

class JobResourceAuthorization(UserGroupAuthorization):
    
    def _is_resource_authorized(
            self, user, permission_type, **kwargs):
        authorized = super(JobResourceAuthorization, self)._is_resource_authorized(
            user, permission_type, **kwargs)
        logger.debug('is JobResource authorized: user: %r: %r', user, authorized)
        if authorized is True:
            return True
        # Allow active staff users
        return user.is_active and user.is_staff

    def is_restricted_view(self, user):
        restricted =  super(JobResourceAuthorization, self).is_restricted_view(user)
        logger.debug('is restricted: %r: %r', user, restricted)
        return restricted 
    
    def filter(self, user, filter_expression):
        if self.is_restricted_view(user):
            logger.info('create job_authorized_user_filter for %r', user.username)
            auth_filter = column('username') == user.username
            if filter_expression is not None:
                filter_expression = and_(filter_expression, auth_filter)
            else:
                filter_expression = auth_filter
        return filter_expression

    def get_row_property_generator(self, user, fields, extant_generator):
        '''
        Filter result properties based on authorization rules
        '''
        return extant_generator


class JobResource(ApiResource):
    
    class Meta:
        authentication = MultiAuthentication(
            IccblBasicAuthentication(), IccblSessionAuthentication())
        resource_name = 'job'
        authorization= JobResourceAuthorization(resource_name)
        ordering = []
        filtering = {} 
        serializer = LimsSerializer()
        excludes = [] 
        always_return_data = True 

    def __init__(self, **kwargs):
        super(JobResource,self).__init__(**kwargs)
    
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, TRAILING_SLASH),
                self.wrap_view('get_schema'), name="api_get_schema"),
            url((r"^(?P<resource_name>%s)/"
                 r"(?P<id>(\d+))%s$") 
                    % (self._meta.resource_name, TRAILING_SLASH),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/test_job%s$" 
                % (self._meta.resource_name, TRAILING_SLASH),
                self.wrap_view('test_job'), name="api_test_job"),
        ]
    
    @read_authorization
    def get_detail(self, request, **kwargs):
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail'] = True
        return self.build_list_response(request, **kwargs)

    @read_authorization
    def get_list(self, request, **kwargs):

        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)

    def build_list_response(self, request, **kwargs):
        
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')

        param_hash = self._convert_request_to_dict(request)
        param_hash.update(kwargs)
        is_data_interchange = param_hash.get(HTTP_PARAM_DATA_INTERCHANGE, False)
        use_vocab = param_hash.get(HTTP_PARAM_USE_VOCAB, False)
        use_titles = param_hash.get(HTTP_PARAM_USE_TITLES, False)
        if is_data_interchange:
            use_vocab = False
            use_titles = False
        
        logger.info('params: %r', param_hash.keys())
        
        is_for_detail = kwargs.pop('is_for_detail', False)
        screen_facility_id = param_hash.pop('screen_facility_id', None)
        if screen_facility_id:
            param_hash['screen_facility_id__eq'] = screen_facility_id
        publication_id = param_hash.pop('publication_id', None)
        if publication_id:
            param_hash['publication_id__eq'] = publication_id
        
        # general setup
      
        manual_field_includes = set(param_hash.get('includes', []))
        
        (filter_expression, filter_hash, readable_filter_hash) = \
            SqlAlchemyResource.build_sqlalchemy_filters(
                schema, param_hash=param_hash)
        filename = self._get_filename(
            readable_filter_hash, schema, is_for_detail)

        filter_expression = \
            self._meta.authorization.filter(request.user,filter_expression)

        order_params = param_hash.get('order_by', [])
        field_hash = self.get_visible_fields(
            schema['fields'], filter_hash.keys(), manual_field_includes,
            param_hash.get('visibilities'),
            exact_fields=set(param_hash.get('exact_fields', [])),
            order_params=order_params)
        order_clauses = SqlAlchemyResource.build_sqlalchemy_ordering(
            order_params, field_hash)
             
        rowproxy_generator = None
        if use_vocab is True:
            rowproxy_generator = \
                DbApiResource.create_vocabulary_rowproxy_generator(field_hash)
        rowproxy_generator = \
           self._meta.authorization.get_row_property_generator(
               request.user, field_hash, rowproxy_generator)
 
        # specific setup
        _job = self.bridge['reports_job']
        _user = self.bridge['reports_userprofile']
        
        j = _job
        j = j.join(
            _user, _job.c.user_requesting_id==_user.c.id)
            
        custom_columns = {}

        base_query_tables = ['reports_job', 'reports_userprofile' ] 
        columns = self.build_sqlalchemy_columns(
            field_hash.values(), base_query_tables=base_query_tables,
            custom_columns=custom_columns)
            
        stmt = select(columns.values()).select_from(j)
        # general setup
         
        (stmt, count_stmt) = self.wrap_statement(
            stmt, order_clauses, filter_expression)
        if not order_clauses and filter_expression is None:
            _alias = Alias(stmt)
            stmt = select([text('*')]).select_from(_alias)
        stmt = stmt.order_by('-id')

        # compiled_stmt = str(stmt.compile(
        #     dialect=postgresql.dialect(),
        #     compile_kwargs={"literal_binds": True}))
        # logger.info('compiled_stmt %s', compiled_stmt)
            
        title_function = None
        if use_titles is True:
            def title_function(key):
                return field_hash[key]['title']
        if is_data_interchange:
            title_function = DbApiResource.datainterchange_title_function(
                field_hash,schema['id_attribute'])
        
        return self.stream_response_from_statement(
            request, stmt, count_stmt, filename,
            field_hash=field_hash,
            param_hash=param_hash,
            is_for_detail=is_for_detail,
            rowproxy_generator=rowproxy_generator,
            title_function=title_function, meta=kwargs.get('meta', None),
            use_caching=True)
    

    @write_authorization
    @un_cache  
    @transaction.atomic      
    def patch_detail(self, request, **kwargs):
        '''
        '''
        schema = kwargs.pop('schema', None)
        if not schema:
            schema = self.build_schema(user=request.user)
#             raise Exception('schema not initialized')
        deserialized = kwargs.pop('data', None)
        # allow for internal data to be passed
        deserialize_meta = None
        if deserialized is None:
            deserialized, deserialize_meta = self.deserialize(
                request, format=kwargs.get('format', None))
        
        if API_RESULT_DATA in deserialized:
            deserialized = deserialized[API_RESULT_DATA][0]
            
        logger.info('patch detail %s, %s', deserialized,kwargs)

        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(
            deserialized, schema=schema, validate=False,**kwargs)
        original_data = None
        if kwargs_for_log:
            try:
                logger.info('patch job: %r, get original state...', kwargs_for_log)
                original_data = self._get_detail_response_internal(**kwargs_for_log)
                kwargs['original_data'] = original_data
            except Exception, e: 
                logger.exception('exception when querying for existing obj: %s', 
                    kwargs_for_log)
        try:
            parent_log = kwargs.get('parent_log', None)
            log = self.make_log(request)
            log.parent_log = parent_log
            log.save()
            kwargs['parent_log'] = log
            patch_response = self.patch_obj(request, deserialized, **kwargs)
            
            obj = patch_response[API_RESULT_OBJ]
            for id_field in id_attribute:
                if id_field not in kwargs_for_log:
                    val = getattr(obj, id_field,None)
                    if val is not None:
                        kwargs_for_log['%s' % id_field] = val
            
        except ValidationError as e:
            logger.exception('Validation error: %r', e)
            raise e

        # get new state, for logging
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        logger.debug('original: %r, new: %r', original_data, new_data)
        log = self.log_patch(
            request, original_data,new_data,log=log, full_create_log=True, 
            **kwargs_for_log)
        if log:
            log.save()
        meta = {}
        if deserialize_meta:
            meta.update(deserialize_meta)
        if API_RESULT_META in patch_response:
            meta = patch_response[API_RESULT_META]
        return self.build_response(
            request,  { API_RESULT_META: meta, API_RESULT_DATA: [new_data,] }, 
            response_class=HttpResponse, **kwargs)
            
    
    @write_authorization
    @transaction.atomic
    def patch_obj(self, request, deserialized, **kwargs):
        
        logger.info('patch job')
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_kwargs = self.get_id(deserialized, schema=schema, validate=False, **kwargs)
        logger.info('id_kwargs: %r, %r', id_kwargs, kwargs)
        create = not bool(id_kwargs)
        
        initializer_dict = self.parse(deserialized, schema=schema, create=create)
        errors = self.validate(initializer_dict, schema=schema,patch=not create)
        if errors:
            raise ValidationError(errors)
        
        if create is False:
            job = Job.objects.get(**id_kwargs)
        else:
            if SCHEMA.JOB.USERNAME not in initializer_dict:
                raise ValidationError(
                    key=SCHEMA.JOB.USERNAME, msg='required')
            username = initializer_dict[SCHEMA.JOB.USERNAME]
            try:
                user = UserProfile.objects.get(username=username)
                logger.info('Job does not exist, creating: %r',initializer_dict)
                job = Job(user_requesting=user)
            except ObjectDoesNotExist:
                raise ValidationError(
                    key=SCHEMA.JOB.USERNAME, msg='user: %r does not exist' % username)
            
            # Enforce business rule: only one job per URI for mutating methods
            method = initializer_dict[SCHEMA.JOB.METHOD]
            if method != 'GET':
                uri = initializer_dict[SCHEMA.JOB.URI]
                query = Job.objects.all().filter(
                    uri=uri, 
                    state__in=[
                        SCHEMA.VOCAB.job.state.PENDING,
                        SCHEMA.VOCAB.job.state.SUBMITTED,
                        SCHEMA.VOCAB.job.state.PROCESSING,
                        ])
                if query.exists():
                    raise ValidationError(
                        key=SCHEMA.JOB.URI,
                        msg=('Pending jobs exist for method: %s uri: %s: '
                            'job ids: %s' % ( method, uri, 
                                ', '.join(map(str,[job.id for job in query])))))
            
        for key, val in initializer_dict.items():
            if hasattr(job, key):
                setattr(job, key, val)
        job.save()
        
        return { API_RESULT_OBJ: job }
    
    @transaction.atomic
    def create_new_job_from_request(self, request, **kwargs):
        ''' 
        Create a Job instance from the original request that is wrapped
        by the background_job wrapper:
        - Stores the request data to reinvoke later from the 
        background job processing client
        @see reports.utils.background_processor
        '''

        JOB = SCHEMA.JOB
        post_data_directory = \
            settings.BACKGROUND_PROCESSOR['post_data_directory']
        if not os.path.exists(post_data_directory):
            os.makedirs(post_data_directory)

        job_data, raw_post_data = \
            background_processor.get_django_request_job_params(request)
        
        try:
            logger.info('create a new job...')
            job_patch_result = self._patch_detail_internal(job_data)
            logger.info('job patch result: %r', job_patch_result)
            new_job_data = job_patch_result[API_RESULT_DATA][0]
        except ValidationError, e:
            logger.exception('on create job')
            response = self.build_error_response(
                request, { API_RESULT_ERROR: e.errors }, **kwargs)
            raise BackgroundJobImmediateResponse(response)
        except Exception, e:
            logger.exception('Job create error: %r, job: %r', e, job_data)
            raise ProgrammingError('Job create error: %r' % e) 

        if raw_post_data:
            logger.info('store request multipart raw_post_data...')
            logger.info('raw post data length: %d', len(raw_post_data))
            
            post_data_filename = \
                os.path.join(post_data_directory, str(new_job_data[JOB.ID]))
            
            with open(post_data_filename,'w') as job_file:
                logger.info('write request body to file: %r', job_file)
                
                job_file.write(raw_post_data)
                job_file.close()
                logger.info('write completed')
        else:
            logger.info('no raw_post_data/body found for request...')
            
        return new_job_data
    
    def service_job(
        self, job_id, wrapped_resource_instance, wrapped_function, **kwargs):
        '''
        Service the job:
        @param wrapped_resource_instance - Resource instance 
        @param wrapped_function - function on the Resource instance that has 
        been decorated with the background_job
        @param kwargs sent with the wrapped function from the dispatcher: 
        - @see reports.api_base.IccblBaseResource.wrap_view
        '''
        
        JOB = SCHEMA.JOB

        logger.info('1. Service the request for job: %r', job_id)
        try:
            recreated_job_request, job_patch_result = \
                self.setup_job_processing_request(job_id)
        except Exception, e:
            logger.exception(
                'on setting job state to processing: %r', job_id)
            job_patch = {
                JOB.ID: job_id,
                JOB.STATE: SCHEMA.VOCAB.job.state.FAILED,
                JOB.RESPONSE_CONTENT: 'Exception created: %r'% e 
                }
            job_patch_result = self._patch_detail_internal(job_patch)
            logger.info('job setup exception patched: %r', 
                job_patch_result)
            raise
                            
        logger.info('--- process the original job request now ---')
        # Perform the exception handling, 
        # as in reports.api_base.IccblBaseResource.wrap_view
        def _inner(*args, **kwargs):
            # Note: wrap the _func because the exception_handler 
            # wrapper expects the func to be bound to "self" already, 
            # and the _func expects "self" to be passed
            logger.info('wrapped_resource_instance: %r', wrapped_resource_instance)
            logger.info('args: %r', args)
            return wrapped_function(wrapped_resource_instance,*args, **kwargs)
        response = wrapped_resource_instance.exception_handler(_inner)(
            wrapped_resource_instance, recreated_job_request, **kwargs)  
        # response = _func(self, recreated_job_request, **kwargs)
        logger.info('--- processing done: response: %r', 
            response.status_code)
        
        if response.status_code >= 300:
            logger.info('error processing job: %r!', 
                response.status_code)
            
            content = self._meta.serializer.deserialize(
                LimsSerializer.get_content(response), JSON_MIMETYPE)
            logger.warn('fail content: %r', content)
            job_patch = {
                JOB.ID: job_id,
                JOB.STATE: SCHEMA.VOCAB.job.state.FAILED,
                JOB.RESPONSE_CONTENT: json.dumps(content),
                JOB.RESPONSE_STATUS_CODE: response.status_code,
                JOB.DATE_TIME_COMPLETED: _now().isoformat() }
            job_patch_result = self._patch_detail_internal(job_patch)
            logger.info('job patch result: %r', job_patch_result)
            new_job_data = job_patch_result[API_RESULT_DATA][0]
        else:
            logger.info('success processing job: %r!', 
                response.status_code)
            content = self._meta.serializer.deserialize(
                LimsSerializer.get_content(response), JSON_MIMETYPE)
            logger.warn('success content: %r', content)
            job_patch = {
                JOB.ID: job_id,
                JOB.STATE: SCHEMA.VOCAB.job.state.COMPLETED,
                JOB.RESPONSE_CONTENT: json.dumps(content),
                JOB.RESPONSE_STATUS_CODE: response.status_code,
                JOB.DATE_TIME_COMPLETED: _now().isoformat(),
                }
            job_patch_result = self._patch_detail_internal(job_patch)
            logger.info('job patch result: %r', job_patch_result)
            new_job_data = job_patch_result[API_RESULT_DATA][0]
        return new_job_data
    
    @transaction.atomic
    def submit_job(self, job_data):
        '''
        Send a pending job to the background_process_script defined in
        settings.BACKGROUND_PROCESSOR['background_process_script']
        '''
        JOB = SCHEMA.JOB
        job_id = job_data[JOB.ID]
        try:
            use_sbatch = 'sbatch_settings' in settings.BACKGROUND_PROCESSOR
            background_client_module_name = \
                settings.BACKGROUND_PROCESSOR['background_process_script']
            logger.info('import background_process_script: %r ', 
                background_client_module_name)
            background_client_module = \
                importlib.import_module(background_client_module_name)
            output = background_client_module.execute_from_python(
                job_id, sbatch=use_sbatch)
            logger.info('submitted job: %r, output: %r', 
                job_id, output)
            if output:
                job_patch = {
                    JOB.ID: job_id,
                    JOB.STATE: SCHEMA.VOCAB.job.state.SUBMITTED,
                    JOB.PROCESS_ID: output,
                    JOB.DATE_TIME_SUBMITTED: _now().isoformat() }
                job_patch_result = self._patch_detail_internal(job_patch)
                logger.info('submitted job: %r', job_patch_result)
                job_data = job_patch_result[API_RESULT_DATA][0]
                
            return job_data
        except subprocess.CalledProcessError, e:
            logger.exception('unexpected exception on job %r submission: %r', 
                job_id, e)
            job_patch = {
                JOB.ID: job_id,
                JOB.STATE: SCHEMA.VOCAB.job.state.FAILED,
                JOB.DATE_TIME_COMPLETED: _now().isoformat(),
                JOB.PROCESS_MESSAGES: e.output 
            }
            job_patch_result = self._patch_detail_internal(job_patch)
            logger.info('submitted job: %r', job_patch_result)
            job_data = job_patch_result[API_RESULT_DATA][0]
            return job_data
        except Exception, e:
            logger.exception('unexpected exception on job %r submission: %r', 
                job_id, e)
            job_patch = {
                JOB.ID: job_id,
                JOB.STATE: SCHEMA.VOCAB.job.state.FAILED,
                JOB.DATE_TIME_COMPLETED: _now().isoformat(),
                JOB.PROCESS_MESSAGES: str(e) 
            }
            job_patch_result = self._patch_detail_internal(job_patch)
            logger.info('submitted job: %r', job_patch_result)
            job_data = job_patch_result[API_RESULT_DATA][0]
            return job_data
            
    @transaction.atomic    
    def setup_job_processing_request(self, job_id, ):
        ''' Use the job_id to recreate the original request to be sent to the
        original resource.
        
        @return recreated_job_request original request with post data
        @return job_patch_result metadata for the current job after setting to
        "processing" state
        '''
        
        JOB = SCHEMA.JOB
        post_data_directory = \
            settings.BACKGROUND_PROCESSOR['post_data_directory']
        if not os.path.exists(post_data_directory):
            os.makedirs(post_data_directory)
    
        query_params = { JOB.ID: job_id }
        job_data = self._get_detail_response_internal(**query_params)
        logger.info('job_data: %r', job_data)
        if not job_data:
            raise ValidationError(key=JOB.ID, msg='no such job: %r' % job_id)
        
        if job_data[JOB.STATE] != SCHEMA.VOCAB.job.state.PENDING:
            logger.warn('NOTICE: job state is not pending: %r', job_data)
            
            if job_data[JOB.STATE] == SCHEMA.VOCAB.job.state.PROCESSING:
                raise ValidationError(key=JOB.STATE,
                    msg='Job: %r is in state: %r,  processing not allowed'
                        % (job_id, job_data[JOB.STATE]))
            
        raw_data = ''
        post_data_filename = \
            os.path.join(post_data_directory, str(job_data[JOB.ID]))
        logger.info(
            'look for request multipart raw_post_data at %r', 
            post_data_filename)
        if os.path.exists(post_data_filename):
            logger.info(
                'found request multipart raw_post_data at %r', 
                post_data_filename)
            with open(post_data_filename,'r') as job_file:
                raw_data = job_file.read()
        else:
            logger.info('no raw data found at: %r', post_data_filename)
            if MULTIPART_MIMETYPE in job_data[JOB.CONTENT_TYPE]:
                raise ValidationError(
                    key=JOB.CONTENT_TYPE,
                    msg='content_type: %r, requires raw data not found at: %r' 
                        % (job_data[JOB.CONTENT_TYPE], post_data_filename))
        
        recreated_job_request = background_processor.create_request_from_job(
            job_data, raw_data)
        try:
            recreated_job_request.user = \
                DjangoUser.objects.get(username=job_data[JOB.USERNAME])
        except ObjectDoesNotExist:
            logger.error('Job user %s not found: %r', 
                job_data[JOB.USERNAME], job_data)
            raise
        
        process_env = background_processor.get_process_env()
        logger.info('process_env: %r', process_env)
        job_patch = {
            JOB.ID: job_id,
            JOB.STATE: SCHEMA.VOCAB.job.state.PROCESSING,
            JOB.DATE_TIME_PROCESSING: _now().isoformat(),
            JOB.PROCESS_ENV: json.dumps(process_env)
         }
        job_patch_result = self._patch_detail_internal(job_patch)
        logger.debug('done setting job state to processing: %r', 
            job_patch_result)
        return recreated_job_request, job_patch_result
        
    @write_authorization
    @background_job
    @un_cache        
    @transaction.atomic    
    def test_job(self, request, **kwargs):
        '''
        Create a test job as a POST operation:
        '''
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        
        return self.build_response(request, { 'test_job': 'created!' })
        