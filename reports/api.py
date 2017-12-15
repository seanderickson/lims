# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from copy import deepcopy
from collections import defaultdict
from decimal import Decimal
import decimal
from functools import wraps
import json
import logging
from operator import itemgetter
import os
import re
import sys
import time
import urllib

from aldjemy.core import get_tables, get_engine
from django.conf import settings
from django.conf.urls import url
from django.contrib.auth.models import User as DjangoUser
from django.contrib.sessions.models import Session
from django.core.cache import cache, caches
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.core.serializers.json import DjangoJSONEncoder
from django.db import transaction
from django.db.models import Q
from django.db.models.aggregates import Max
from django.forms.models import model_to_dict
from django.http.request import HttpRequest
from django.http.response import HttpResponse, Http404, HttpResponseBase,\
    HttpResponseNotFound
from sqlalchemy import select, asc, text
from sqlalchemy.dialects import postgresql
from sqlalchemy.dialects.postgresql import array
from sqlalchemy.sql import and_, or_, not_, asc, desc, func
from sqlalchemy.sql.elements import literal_column
from sqlalchemy.sql.expression import column, join, distinct, exists
from tastypie.authentication import BasicAuthentication, MultiAuthentication
from tastypie.exceptions import BadRequest
from tastypie.utils.urls import trailing_slash

from reports import LIST_DELIMITER_SQL_ARRAY, LIST_DELIMITER_URL_PARAM, \
    HTTP_PARAM_USE_TITLES, HTTP_PARAM_USE_VOCAB, HEADER_APILOG_COMMENT, \
    HTTP_PARAM_DATA_INTERCHANGE, InformationError
from reports import ValidationError, _now
from reports.api_base import IccblBaseResource, un_cache, Authorization, \
    IccblSessionAuthentication
from reports.models import MetaHash, Vocabulary, ApiLog, LogDiff, Permission, \
                           UserGroup, UserProfile, API_ACTION_DELETE, \
                           API_ACTION_CREATE
from reports.serialize import parse_val, parse_json_field, XLSX_MIMETYPE, \
    SDF_MIMETYPE, XLS_MIMETYPE
from reports.serializers import LimsSerializer
from reports.sqlalchemy_resource import SqlAlchemyResource, _concat


logger = logging.getLogger(__name__)

URI_VERSION = 'v1'
BASE_URI = '/reports/api/' + URI_VERSION

API_RESULT_OBJ = 'object'
API_RESULT_DATA = 'objects'
API_RESULT_META = 'meta'

API_PARAM_OVERRIDE = 'override'
API_MSG_SUBMIT_COUNT = 'Data submitted'
API_MSG_RESULT = 'Result'
API_MSG_WARNING = 'Warning'
API_MSG_NOT_ALLOWED = 'Action not allowed'
API_MSG_UPDATED = 'Updated'
API_MSG_CREATED = 'Created'
API_MSG_UNCHANGED = 'Unchanged'
API_MSG_COMMENTS = 'Comments'
API_MSG_ACTION = 'Action'

DEBUG_RESOURCES = False or logger.isEnabledFor(logging.DEBUG)
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
        self, user, permission_type, **kwargs):
        
        if DEBUG_AUTHORIZATION:
            logger.info("_is_resource_authorized: %s, user: %s, type: %s",
                self.resource_name, user, permission_type)
        scope = 'resource'
        prefix = 'permission'
        uri_separator = '/'
        permission_str =  uri_separator.join([
            prefix,scope,self.resource_name,permission_type])       

        if DEBUG_AUTHORIZATION:
            logger.info('authorization query: %s, user %s, %s' 
                % (permission_str, user, user.is_superuser))
        
        if user.is_superuser:
            if DEBUG_AUTHORIZATION:
                logger.info('%s:%s access allowed for super user: %s' 
                    % (self.resource_name,permission_type,user))
            return True
        
        userprofile = user.userprofile
        permission_types = [permission_type]
        if permission_type == 'read':
            permission_types.append('write')
        query = userprofile.permissions.all().filter(
            scope=scope, key=self.resource_name, type__in=permission_types)
        if query.exists():
            if DEBUG_AUTHORIZATION:
                logger.info(
                    'user %s, auth query: %s, found matching user permissions %s'
                    % (user,permission_str,[str(x) for x in query]))
            logger.info('%s:%s user explicit permission for: %s' 
                % (self.resource_name,permission_type,user))
            return True
        
        if DEBUG_AUTHORIZATION:
            logger.info(
                'user %s, auth query: %s, not found in user permissions %s'
                % (user,permission_str,[str(x) for x in query]))
        
        permissions_group = [ permission 
                for group in userprofile.usergroup_set.all() 
                for permission in group.get_all_permissions(
                    scope=scope, key=self.resource_name, type__in=permission_types)]
        if permissions_group:
            if DEBUG_AUTHORIZATION:
                logger.info(
                    'user: %r, auth query: %r, found usergroup permissions: %r'
                    ,user, permission_str, permissions_group)
            logger.debug('%s:%s usergroup permission for: %s' 
                % (self.resource_name,permission_type,user))
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
            logger.info('write auth failed for: %r, %r', 
                request.user, resource_name)
            msg = 'write auth failed for: %r, %r' % (request.user, resource_name)
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
            logger.info('read auth for: %r using: %r', resource_name, self._meta.authorization)
        if not self._meta.authorization._is_resource_authorized(
                request.user,'read', **kwargs):
            msg = 'read auth failed for: %r, %r' % (request.user, resource_name)
            logger.warn(msg)
            raise PermissionDenied(msg)

        # Use the user specific settings to build the schema
        if resource_name != 'screenresult':
            kwargs['schema'] = self.build_schema(user=request.user)
        
        return _func(self, *args, **kwargs)

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
    

def compare_dicts(dict1, dict2, excludes=None, exclude_patterns=None):
    '''
    @param full (default False) 
    - a full compare shows added keys as well as diff keys
    
    Note: "full" logs would log all of the created data in the resource;
    whereas "not full" only logs the creation; with this strategy, logs 
    must be played backwards to recreate an entity state.
    '''
    if DEBUG_PATCH_LOG:
        logger.info('compare dicts: %r, %r, excludes: %r, %r', 
            dict1, dict2, excludes, exclude_patterns)
    
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
    diffs = {}
    for key in union_keys:
        v1 = dict1.get(key, None)
        v2 = dict2.get(key, None)
        if v1 is not None or v2 is not None:
            if v1 != v2:
                diffs[key] = [v1,v2]
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
        self.resource_resource = None

    def get_resource_resource(self):
        if not self.resource_resource:
            self.resource_resource = ResourceResource()
        return self.resource_resource
    
    def get_schema(self, request, **kwargs):
    
        return self.build_response(request, self.build_schema(user=request.user), **kwargs)

    @read_authorization
    def get_detail(self, request, **kwargs):
        '''
        Override the Tastypie 
        operations will be handled using SqlAlchemy
        '''
        raise NotImplemented(
            'get_detail must be implemented for resource: %r' 
            % self._meta.resource_name)
    
    @read_authorization
    def get_list(self, request, **kwargs):
        '''
        Override the Tastypie 
        operations will be handled using SqlAlchemy
        '''
        raise NotImplemented(
            'get_list must be implemented for resource: %r' 
            % self._meta.resource_name)
        
    def build_list_response(self,request, **kwargs):
        raise NotImplemented(
            'build_list_response must be implemented for the SqlAlchemyResource: %r' 
            % self._meta.resource_name)

    def build_schema(self, user=None):
        
        logger.debug('build schema for: %r: %r', self._meta.resource_name, user)
        schema = self.get_resource_resource()._get_resource_schema(
            self._meta.resource_name, user)
        if DEBUG_RESOURCES:
            logger.info('schema fields: %r', schema['fields'].keys())
        return schema
    
    def get_id(self,deserialized, validate=False, schema=None, **kwargs):
        ''' 
        return the full ID for the resource, as defined by the "id_attribute"
        - if validate=True, raises ValidationError if some or all keys are missing
        - otherwise, returns keys that can be found
        '''
        if schema is None:
            raise ProgrammingError
        
        id_attribute = schema['id_attribute']
        fields = schema['fields']
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
            logger.info('create schema for find_key_from_resource_uri: %r', resource_uri)
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
            raise NotImplementedError(
                'resource uri %r does not contain all id attributes: %r'
                % (resource_uri,id_attribute))
        else:
            return dict(zip(id_attribute,keys))

    def parse(self,deserialized, create=False, fields=None, schema=None):
        ''' parse schema fields from the deserialized dict '''
        DEBUG_PARSE = False or logger.isEnabledFor(logging.DEBUG)
        if DEBUG_PARSE:
            logger.info('parse: %r:%r', self._meta.resource_name, schema is None)
            logger.info('parse: %r', deserialized)
       
        mutable_fields = fields        
        if fields is None:
            
            if schema is None:
                logger.info('build raw schema for parse')
                schema = self.build_schema()
        
            fields = schema['fields']
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
            _val = None
            if key in deserialized:
                _val = parse_val(
                    deserialized.get(key,None), key,field['data_type']) 
                if DEBUG_PARSE:
                    logger.info('parsing: %r, %r, %r',
                        key, deserialized.get(key,None), _val)
                if _val is None:
                    _val = field.get('default',None)
                    if _val == '': 
                        _val = None
                initializer_dict[key] = _val
        return initializer_dict
    
    @read_authorization
    def search(self, request, **kwargs):
        '''
        Treats a POST request like a GET (with posted raw_search_data):
        - for large raw_search_data (too large for URL param encoding)
        '''
         
        DEBUG_SEARCH = False or logger.isEnabledFor(logging.DEBUG)
         
        search_ID = kwargs['search_ID']
         
        all_params = self._convert_request_to_dict(request)
        all_params.update(kwargs)

        # NOTE: remove "nested_search_data" 20170515;
        # - this data must be sent with each request
        # nested_search_data = all_params.get('nested_search_data', None)
        raw_search_data = all_params.get('raw_search_data', None)
         
        if raw_search_data:
            raw_search_data = urllib.unquote(raw_search_data)
            # cache the search data on the session, to support subsequent requests
            # to download or modify
            request.session[search_ID] = raw_search_data  
        else:
            if search_ID in request.session:
                raw_search_data = request.session[search_ID]
            else:
                raise BadRequest(
                    'raw_search_data for id missing: %r, %r.'
                    '.search requires a "raw_search_data" param'
                    % (search_ID, self._meta.resource_name))
        
        if DEBUG_SEARCH:
            logger.info('raw_search_data: %r', raw_search_data)
        kwargs['raw_search_data'] = raw_search_data
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
 
        response = self.build_list_response(request,**kwargs)
        
        # Because this view bypasses the IccblBaseResource.dispatch method, we
        # are implementing the downloadID cookie here, for now.
        downloadID = all_params.get('downloadID', None)
        if downloadID:
            logger.info('set cookie downloadID: %r', downloadID)
            response.set_cookie('downloadID', downloadID)
        else:
            logger.info('no downloadID')
        
        return response
    
    @write_authorization
    @un_cache
    @transaction.atomic   
    def patch_list(self, request, **kwargs):

        logger.info('patch list, user: %r, resource: %r' 
            % ( request.user.username, self._meta.resource_name))
        logger.debug('patch list: %r' % kwargs)

        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
        
        if len(deserialized) == 0:
            meta = { 
                API_MSG_RESULT: {
                    API_MSG_SUBMIT_COUNT : 0, 
                    API_MSG_UPDATED: 0, 
                    API_MSG_CREATED: 0,
                    API_MSG_UNCHANGED: 0, 
                    API_MSG_COMMENTS: 'no data patched'
                }
            }
            logger.info('PATCH list: %r', meta)
            return self.build_response(
                request, { 'meta': meta }, response_class=HttpResponse, **kwargs)

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
        kwargs_for_log['schema'] = schema
#         if 'schema' in kwargs_for_log:
#             del kwargs_for_log['schema']
        for _data in deserialized:
            id_kwargs = self.get_id(_data, schema=schema)
            logger.debug('found id_kwargs: %r from %r', id_kwargs, _data)
            if id_kwargs:
                for idkey,idval in id_kwargs.items():
                    id_param = '%s__in' % idkey
                    id_vals = kwargs_for_log.get(id_param, set())
                    id_vals.add(idval)
                    kwargs_for_log[id_param] = id_vals
        try:
            logger.debug('get original state, for logging... %r', kwargs_for_log)
            original_data = self._get_list_response_internal(**kwargs_for_log)
            logger.info('original state retrieved: %d', len(original_data))
        except Exception as e:
            logger.exception('original state not obtained')
            original_data = []

        try:
            # FIXME: move transaction to method decorator
            with transaction.atomic():
                if 'parent_log' not in kwargs:
                    parent_log = self.make_log(request, schema=schema)
                    parent_log.key = self._meta.resource_name
                    parent_log.uri = self._meta.resource_name
                    parent_log.save()
                    kwargs['parent_log'] = parent_log
                logger.info('perform patch_list: %d', len(deserialized))
                for _dict in deserialized:
                    self.patch_obj(request, _dict, **kwargs)
        except ValidationError as e:
            logger.exception('Validation error: %r', e)
            raise e
            
        logger.debug('Get new state, for logging: %r...',
            {k:v for k,v in kwargs_for_log.items() if k != 'schema'})
        new_data = self._get_list_response_internal(**kwargs_for_log)
        logger.info('new data: %d, log patches...', len(new_data))
        logs = self.log_patches(request, original_data,new_data,schema=schema,**kwargs)
        logger.info('patch logs created: %d', len(logs) if logs else 0 )
        patch_count = len(deserialized)
        update_count = len([x for x in logs if x.diffs ])
        logger.debug('updates: %r', [x for x in logs if x.diffs ])
        create_count = len([x for x in logs if x.api_action == API_ACTION_CREATE])
        unchanged_count = patch_count - update_count
        meta = { 
            API_MSG_RESULT: {
                API_MSG_SUBMIT_COUNT : patch_count, 
                API_MSG_UPDATED: update_count, 
                API_MSG_CREATED: create_count,
                API_MSG_UNCHANGED: unchanged_count, 
                API_MSG_COMMENTS: parent_log.comment
            }
        }
        if not self._meta.always_return_data:
            return self.build_response(
                request, { 'meta': meta }, response_class=HttpResponse)
        else:
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
        logger.debug('schema field keys: %r', schema['fields'].keys())
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]
        
        # Limit the potential candidates for logging to found id_kwargs
        id_attribute = resource = schema['id_attribute']
        kwargs_for_log = kwargs.copy()
        kwargs_for_log['schema'] = schema
        for _data in deserialized:
            id_kwargs = self.get_id(_data, schema=schema)
            logger.debug('found id_kwargs: %r from %r', id_kwargs, _data)
            if id_kwargs:
                for idkey,idval in id_kwargs.items():
                    id_param = '%s__in' % idkey
                    id_vals = kwargs_for_log.get(id_param, [])
                    id_vals.append(idval)
                    kwargs_for_log[id_param] = id_vals
        try:
            logger.info('get original state, for logging...')
            logger.debug('kwargs_for_log: %r', kwargs_for_log)
            original_data = self._get_list_response_internal(**kwargs_for_log)
        except Exception as e:
            logger.exception('original state not obtained')
            original_data = []

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
        self._meta.queryset.delete()
        new_objs = []
        logger.info('PUT: create new objs: %d', len(deserialized))
        for _dict in deserialized:
            new_objs.append(self.put_obj(request, _dict, **kwargs))

        # Get new state, for logging
        # After patch, the id keys must be present
        for idkey in id_attribute:
            id_param = '%s__in' % idkey
            ids = set(kwargs_for_log[id_param])
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
        self.log_patches(request, original_data,new_data,schema=schema,**kwargs)
        
        logger.info('put_list done.')
        if not self._meta.always_return_data:
            logger.info('put success, no data')
            return HttpResponse(status=202)
        else:
            response = self.get_list(request, **kwargs)             
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

        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        if self._meta.collection_name in deserialized:
            deserialized = deserialized[self._meta.collection_name]

        if len(deserialized) == 0:
            meta = { 
                API_MSG_RESULT: {
                    API_MSG_SUBMIT_COUNT : 0, 
                    API_MSG_UPDATED: 0, 
                    API_MSG_CREATED: 0,
                    API_MSG_UNCHANGED: 0, 
                    API_MSG_COMMENTS: 'no data posted'
                }
            }
            return self.build_response(
                request, { 'meta': meta }, response_class=HttpResponse, **kwargs)
        
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
        for _data in deserialized:
            id_kwargs = self.get_id(_data, schema=schema)
            logger.debug('found id_kwargs: %r from %r', id_kwargs, _data)
            if id_kwargs:
                for idkey,idval in id_kwargs.items():
                    id_param = '%s__in' % idkey
                    id_vals = kwargs_for_log.get(id_param, set())
                    id_vals.add(idval)
                    kwargs_for_log[id_param] = id_vals
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
        
        for _dict in deserialized:
            self.patch_obj(request, _dict, **kwargs)
        new_objs = []
        for _dict in deserialized:
            new_objs.append(self.patch_obj(request, _dict))

        # Get new state, for logging
        # After patch, the id keys must be present
        for idkey in id_attribute:
            id_param = '%s__in' % idkey
            ids = set(kwargs_for_log[id_param])
            for new_obj in new_objs:
                if hasattr(new_obj, idkey):
                    idval = getattr(new_obj, idkey)
                    ids.add(idval)
            kwargs_for_log[id_param] = ids
        new_data = self._get_list_response_internal(**kwargs_for_log)
        logger.debug('post list done, new data: %d', len(new_data))

        logs = self.log_patches(request, original_data,new_data,schema=schema, **kwargs)
        logger.info('post logs created: %d', len(logs) if logs else 0 )
        patch_count = len(deserialized)
        update_count = len([x for x in logs if x.diffs ])
        create_count = len([x for x in logs if x.api_action == API_ACTION_CREATE])
        unchanged_count = patch_count - update_count
        meta = { 
            API_MSG_RESULT: {
                API_MSG_SUBMIT_COUNT: patch_count, 
                API_MSG_UPDATED: update_count, 
                API_MSG_CREATED: create_count,
                API_MSG_UNCHANGED: unchanged_count, 
                API_MSG_COMMENTS: parent_log.comment
            }
        }
        logger.info('POST list: %r', meta)
        if not self._meta.always_return_data:
            return self.build_response(
                request, { 'meta': meta }, response_class=HttpResponse)
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
        
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs.pop('data')
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))

        logger.debug('post detail %s, %s', 
            { k:v for k,v in kwargs.items() if k != 'schema'},deserialized)
        
        _data = self.build_post_detail(request, deserialized, **kwargs)
        
        return self.build_response(
            request, _data, response_class=HttpResponse, **kwargs)
    
    def build_post_detail(self, request, deserialized, log=None, **kwargs):
        ''' Inner post detail method, returns native dict '''

        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
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
        logger.info('build post detail: %r', obj)
        for id_field in id_attribute:
            if id_field not in kwargs_for_log:
                val = getattr(obj, id_field,None)
                if val is not None:
                    kwargs_for_log['%s' % id_field] = val
                else:
                    logger.warn('id field: %r not found in new obj: %r', id_field, obj)

        # get new state, for logging
        logger.debug('get new state, for logging, kwargs: %r', kwargs_for_log)
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        logger.debug('new post data: %r', new_data)
        if not new_data:
            raise BadRequest(
                'no data found for the new obj created by post: %r', obj)
        patched_log = self.log_patch(
            request, original_data,new_data,log=log, 
            id_attribute=id_attribute, schema=schema, **kwargs)
        if patched_log:
            patched_log.save()
            logger.debug('post log: %r', patched_log)
        
        # 20170109 - return complex data
        new_data = { API_RESULT_DATA: [new_data,] }
        if API_RESULT_META in patch_result:
            new_data[API_RESULT_META] = patch_result[API_RESULT_META]
        return new_data
    
    @write_authorization
    @un_cache  
    @transaction.atomic      
    def patch_detail(self, request, **kwargs):
        '''
        PATCH is used to create or update a resource; not idempotent
        '''
        
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs.pop('data')
        else:
            deserialized = self.deserialize(
                request, format=kwargs.get('format', None))
        _data = self.build_patch_detail(request, deserialized, **kwargs)

        return self.build_response(
            request, _data, response_class=HttpResponse, **kwargs)

    def build_patch_detail(self, request, deserialized, **kwargs):
        # cache state, for logging
        # Look for id's kwargs, to limit the potential candidates for logging
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        id_attribute = schema['id_attribute']
        kwargs_for_log = self.get_id(
            deserialized, validate=True, schema=schema, **kwargs)

        # NOTE: creating a log, even if no data have changed (may be comment only)
        log = self.make_log(request, schema=schema)
        original_data = None
        if kwargs_for_log:
            try:
                original_data = self._get_detail_response_internal(**kwargs_for_log)
                if not log.key:
                    self.make_log_key(log, original_data, id_attribute=id_attribute,
                        schema=schema, **kwargs)
                logger.debug('original data: %r', original_data)
            except Exception, e: 
                logger.exception('exception when querying for existing obj: %s', 
                    kwargs_for_log)
        log.save()
        patch_result = self.patch_obj(request, deserialized, log=log, **kwargs)
        if API_RESULT_OBJ in patch_result:
            obj = patch_result[API_RESULT_OBJ]
        else:
            # TODO: 20170109, legacy, convert patch_obj to return:
            # { API_RESULT_OBJ, API_RESULT_META }
            obj = patch_result
        logger.debug('build patch detail: %r', obj)
        
        for id_field in id_attribute:
            if id_field not in kwargs_for_log:
                val = getattr(obj, id_field,None)
                if val is not None:
                    kwargs_for_log['%s' % id_field] = val
        new_data = self._get_detail_response_internal(**kwargs_for_log)
        logger.debug('new_data: %r', new_data)
        patched_log = self.log_patch(
            request, original_data,new_data,log=log, 
            id_attribute=id_attribute, schema=schema, **kwargs)
        if patched_log:
            patched_log.save()
            logger.debug('patch log: %r', patched_log)
        else:
            logger.info('no patch log')
        # 20170109 - return complex data
        new_data = { API_RESULT_DATA: [new_data,], }
        if API_RESULT_META in patch_result:
            new_data[API_RESULT_META] = patch_result[API_RESULT_META]
            
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
        raise NotImplementedError('put_detail must be implemented')
    
#         schema = kwargs.pop('schema', None)
#         if not schema:
#             raise Exception('schema not initialized')
# 
#         if kwargs.get('data', None):
#             # allow for internal data to be passed
#             deserialized = kwargs['data']
#         else:
#             deserialized = self.deserialize(
#                 request, format=kwargs.get('format', None))
# 
#         logger.debug('put detail: %r, %r' % (deserialized,kwargs))
#         
#         # cache state, for logging
#         # Look for id's kwargs, to limit the potential candidates for logging
#         # Note: this version of PUT cannot be used if the resource must create
#         # the ID keys on create (use POST/PATCH)
#         
#         kwargs_for_log = self.get_id(
#             deserialized, validate=True, schema=schema, **kwargs)
#         id_attribute = schema['id_attribute']
#         # kwargs_for_log = {}
#         # for id_field in id_attribute:
#         # if deserialized.get(id_field,None):
#         #     kwargs_for_log[id_field] = deserialized[id_field]
#         # elif kwargs.get(id_field,None):
#         #     kwargs_for_log[id_field] = kwargs[id_field]
#         logger.debug('put detail: %s, %s' %(deserialized,kwargs_for_log))
#         try:
#             logger.info('get original state, for logging...')
#             logger.debug('kwargs_for_log: %r', kwargs_for_log)
#             original_data = self._get_list_response_internal(**kwargs_for_log)
#         except Exception as e:
#             logger.exception('original state not obtained')
#             original_data = []
#         
#         try:
#             logger.debug('call put_obj')
#             obj = self.put_obj(request, deserialized, **kwargs)
#         except ValidationError as e:
#             logger.exception('Validation error: %r', e)
#             raise e
#                 
#         if not kwargs_for_log:
#             for id_field in id_attribute:
#                 val = getattr(obj, id_field,None)
#                 kwargs_for_log['%s' % id_field] = val
# 
#         # get new state, for logging
#         new_data = self._get_list_response_internal(**kwargs_for_log)
#         self.log_patches(request, original_data,new_data,**kwargs)
#         
#         # TODO: add "Result" data to meta section, see patch_list
#         
#         if not self._meta.always_return_data:
#             logger.info('put success, no response data')
#             return HttpResponse(status=200)
#         else:
#             response.status_code = 200
#             return response

    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_list(self, request, **kwargs):
        logger.info('delete list...')
        raise NotImplementedError('delete_list is not implemented for %s'
            % self._meta.resource_name )

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
        log.api_action = API_ACTION_DELETE
        log.save()
        
        self.delete_obj(request, {}, log=log, **kwargs_for_log)

        # Log
        logger.info('deleted: %s' %kwargs_for_log)
#         log_comment = None
#         if HEADER_APILOG_COMMENT in request.META:
#             log_comment = request.META[HEADER_APILOG_COMMENT]
        
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
        raise NotImplementedError('delete obj must be implemented')
    
    @write_authorization
    @transaction.atomic    
    def patch_obj(self, request, deserialized, **kwargs):
        raise NotImplementedError('patch obj must be implemented')

    def validate(self, _dict, patch=False, schema=None):
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
        if schema is None:
            logger.info('build raw schema for validate')
            schema = self.build_schema()
        fields = schema['fields']
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
                    if name in id_attribute:
                        continue
                    editability = field.get('editability',None)
                    if not editability or 'u' not in editability:
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
        
        def title_function(key):
            new_title = '%s (not updatable)' % key
            # editable_flags = set(['c','u'])
            logger.debug('key: %r, %r', key, key in field_hash)
            if key in field_hash:
                fi = field_hash[key]
                editability = fi.get('editability',[])
                logger.info('editability: %r', editability)
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
        Adjust fields to their "default_unit" as defined in the "display_options"
        - e.g. if default_unit = 1e-6:
            adjust the raw_val = raw_val.scaleb(6)
        - e.g. if default_unit = 1e6:
            adjust the raw_val = raw_val.scaleb(-6)
        '''
        DEFAULT_PRECISION = Decimal('1e-%d'%9)
        siunit_default_units = {}
        for key, field in field_hash.iteritems():
            if field.get('display_type', None) == 'siunit':
                display_options = field.get('display_options', None)
                try:
                    if display_options is not None:
                        display_options = display_options.replace(r"'", '"')
                        display_options = json.loads(display_options)
                        logger.debug('decoded display_options: %r', display_options)
                    default_unit = display_options.get('defaultUnit',None)
                    multiplier = display_options.get('multiplier', None)
                    decimals = display_options.get('decimals', None)
                    if default_unit is not None:
                        # get the scale (exponent) of the default unit
                        # negate the scale for use with Decimal.scaleb()
                        scale = -Decimal(str(default_unit)).adjusted()
                        if multiplier is not None:
                            # get the exponent of the multiplier for use with scaleb
                            multiplier = Decimal(str(multiplier)).adjusted()
                            if multiplier != 0:
                                scale = scale+multiplier
                        if scale != 0:
                            _dict = { 'scale': scale }
                            if decimals is not None:
                                _dict['decimals'] = Decimal('1e-%d'%int(decimals))
                            else: 
                                _dict['decimals'] = DEFAULT_PRECISION
                            siunit_default_units[key] = _dict
                except Exception, e:
                    logger.exception(
                        'key: %r, error in default unit and multipler: %r', 
                        key, display_options)
        logger.info('siunit_default_units: %r', siunit_default_units)
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
                        raw_val = row[key]
                        try:
                            scale = siunit_default_units[key]['scale']
                            decimals = siunit_default_units[key]['decimals']
                            val = Decimal(raw_val)
                            val = val.scaleb(scale)
                            val = val.quantize(decimals, decimal.ROUND_HALF_UP)
                            logger.debug('converted: %r, %r:%r', key, raw_val, val)
                            return val
                        except Exception, e:
                            logger.exception(
                                'key: %r, error in creating siunit adjusted val for: %r', 
                                key, raw_val)
                            return self.row[key]
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
                if not self.row[key]:
                    return None
                if key in vocabularies:
                    if self.row[key] not in vocabularies[key]:
                        logger.error(
                            ('Unknown vocabulary:'
                             ' scope:%s key:%s val:%r, keys defined: %r'),
                            field_hash[key]['vocabulary_scope_ref'], key, 
                            self.row[key],vocabularies[key].keys() )
                        return self.row[key] 
                    else:
                        return vocabularies[key][self.row[key]]['title']
                else:
                    return self.row[key]
        
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
        
        logger.debug('make_log_key: %r, %r', id_attribute, schema is not None)
        
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
            full_create_log=False, **kwargs):

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
        if HEADER_APILOG_COMMENT in request.META:
            log.comment = request.META[HEADER_APILOG_COMMENT]
        if not log.key:
            self.make_log_key(log, new_dict, id_attribute=id_attribute,
                **kwargs)

        if 'parent_log' in kwargs:
            log.parent_log = kwargs.get('parent_log', None)
        if prev_dict:
            log.diffs = compare_dicts(prev_dict,new_dict, excludes,exclude_patterns)
            if not log.diffs:
                if DEBUG_PATCH_LOG:
                    logger.info('no diffs found: %r, %r' 
                        % (prev_dict,new_dict))
                else:
                    logger.info('no diffs found: %r', log.uri) 
                if log.comment is None:
                    # if no comment, no need to save the log
                    log = None
            else:
                logger.info('PATCH: %r', log.uri) 
                
        else: # creating
            log.api_action = API_ACTION_CREATE
            if full_create_log is True:
                log.diffs = { 
                    key: [None,val] for key,val in new_dict.items() 
                        if val is not None }
            if DEBUG_PATCH_LOG:
                logger.info('create, api log: %r', log)
            else:
                logger.info('CREATE: %r', log.uri) 
    
        if DEBUG_PATCH_LOG:
            logger.info('log patch done: %r', log)
        return log

    @transaction.atomic
    def log_patches(self,request, original_data, new_data, schema=None, **kwargs):
        '''
        log differences between dicts having the same identity in the arrays:
        @param original_data - data from before the API action
        @param new_data - data from after the API action
        - dicts have the same identity if the id_attribute keys have the same
        value.
        '''
        logs = []
        
#         log_comment = None
#         if HEADER_APILOG_COMMENT in request.META:
#             log_comment = request.META[HEADER_APILOG_COMMENT]
        
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
                schema=schema, **kwargs)  
            if DEBUG_PATCH_LOG:
                logger.info('patch log: %r', log)          
            if log:
                logs.append(log)   
            
        for deleted_dict in deleted_items:
            
            log = self.make_log(
                request, attributes=deleted_dict, id_attribute=id_attribute,
                schema=schema)

            # user can specify any valid, escaped json for this field
            # if 'apilog_json_field' in bundle.data:
            #     log.json_field = bundle.data['apilog_json_field']
            
            if 'parent_log' in kwargs:
                log.parent_log = kwargs.get('parent_log', None)

            log.api_action = API_ACTION_DELETE
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
            BasicAuthentication(), IccblSessionAuthentication())
        resource_name='apilog' 
        authorization= ApiLogAuthorization(resource_name) #Authorization()        
        ordering = []
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        max_limit = 100000
    
    def __init__(self, **kwargs):
        self.scope = 'fields.apilog'
        super(ApiLogResource,self).__init__(**kwargs)

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/clear_cache%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_clear_cache'), name="api_clear_cache"),
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)/children%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_apilog_childview'), 
                name="api_dispatch_apilog_childview"),
            url((r"^(?P<resource_name>%s)/children"
                 r"/(?P<ref_resource_name>[\w\d_.\-:]+)"
                 r"/(?P<key>[\w\d_.\-\+: \/]+)"
                 r"/(?P<date_time>[\w\d_.\-\+:]+)%s$")
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_apilog_childview2'), 
                name="api_dispatch_apilog_childview2"),
            url((r"^(?P<resource_name>%s)/(?P<ref_resource_name>[\w\d_.\-:]+)"
                 r"/(?P<key>[\w\d_.\-\+: \/]+)"
                 r"/(?P<date_time>[\w\d_.\-\+:]+)%s$")
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
        
    def dispatch_clear_cache(self, request, **kwargs):
        self.clear_cache()
        return self.build_response(request, 'ok', **kwargs)

    @read_authorization
    def get_detail(self, request, **kwargs):

        id = kwargs.get('id', None)
        if id:
            return self.get_list(request, **kwargs)
            
        ref_resource_name = kwargs.get('ref_resource_name', None)
        if not ref_resource_name:
            logger.info('no ref_resource_name provided')
            raise NotImplementedError(
                'must provide a ref_resource_name parameter')
        
        key = kwargs.get('key', None)
        if not key:
            logger.info('no key provided')
            raise NotImplementedError('must provide a key parameter')
        
        date_time = kwargs.get('date_time', None)
        if not date_time:
            logger.info('no date_time provided')
            raise NotImplementedError('must provide a date_time parameter')

        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
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
            filename = self._get_filename(readable_filter_hash, schema)

            if filter_expression is None and 'parent_log_id' not in kwargs:
                raise InformationError(
                    key='Filter Data',
                    msg='Please enter a filter expression to begin')
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                                  
            order_params = param_hash.get('order_by',[])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes, 
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
                raise NotImplementedError('Diff filtering is not implemented')
                
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
                query = (
                    select([
                        _logdiff.c.field_key,
                        array([_logdiff.c.before,_logdiff.c.after])])
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
                                diff_dict = { x[0]:x[1] for x in _diffs }
                                return json.dumps(diff_dict)
                            else:
                                return self.row[key]
                    conn = get_engine().connect()
                    try:
                        for row in cursor:
                            yield Row(row)
                    finally:
                        conn.close()
                        
                return diff_generator
            
            if 'diffs' in field_hash:
                rowproxy_generator = create_diff_generator(rowproxy_generator)
            
            # compiled_stmt = str(stmt.compile(
            #     dialect=postgresql.dialect(),
            #     compile_kwargs={"literal_binds": True}))
            # logger.info('compiled_stmt %s', compiled_stmt)
            
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get_list')
            raise e  
    
    def build_schema(self, user=None):
        schema = super(ApiLogResource,self).build_schema(user=user)
        temp = [ x.key for x in 
            MetaHash.objects.all().filter(scope='resource').distinct('key')]
        schema['extraSelectorOptions'] = { 
            'label': 'Resource', 
            'searchColumn': 'ref_resource_name', 'options': temp }
        return schema        
    
    def dispatch_apilog_childview(self, request, **kwargs):
        kwargs['parent_log_id'] = kwargs.pop('id')
        return ApiLogResource().dispatch('list', request, **kwargs)    

    def dispatch_apilog_childview2(self, request, **kwargs):
        parent_log = self._get_detail_response_internal(**kwargs)

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
            BasicAuthentication(), IccblSessionAuthentication())
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

    def clear_cache(self):
        ApiResource.clear_cache(self)
        self.get_resource_resource().clear_cache()
        
    def get_resource_resource(self):
        if self.resource_resource is None:
            self.resource_resource = ResourceResource()
        return self.resource_resource    
        
    def prepend_urls(self):

        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.]+)/(?P<key>[\w\d_]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
    
    def create_fields(self):
        pass

    def build_schema(self, user=None):
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
            'fields': field_schema,
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
        filenames = ['fields']
        
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
            scopes = MetaHash.objects.all().filter(
                scope__icontains='fields.').values_list('scope',flat=True).distinct()
            if not scopes.exists():
                # bootstrap case
                scopes = ['fields.field',]
                
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
                raise BadRequest(
                    "Invalid limit '%s' provided. Please provide a positive integer." 
                    % limit)
            offset = param_hash.get('offset', 0 )
            try:
                offset = int(offset)
            except Exception:
                raise BadRequest(
                    "Invalid offset '%s' provided. Please provide a positive integer." 
                    % offset)
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
            
            if kwargs.get('meta', None):
                temp = kwargs['meta']
                logger.debug('meta found in kwargs: %r', temp)
                temp.update(meta)
                meta = temp
                logger.debug('meta: %r', meta)
            logger.debug('meta: %r', meta)
            response_hash = { 
                'meta': meta, 
                self._meta.collection_name: fields 
            }
        kwargs['filename'] = '_'.join(filenames)
        logger.debug('FieldResource build response: %r, %r, %r', 
            self._meta.resource_name, request, 
            {k:v for k,v in kwargs.items() if k!='schema'})

        return self.build_response(request, response_hash, **kwargs)

    def _build_fields(self, scopes=None):
        ''' Internal callers - build the schema.fields hash
        '''
        if not scopes:
            scopes = MetaHash.objects.all().filter(
                scope__icontains='fields.').values_list('scope',flat=True).distinct()
            if not scopes.exists():
                # bootstrap case
                scopes = ['fields.field',]

        fields = []
        for scope in scopes:
            field_hash = deepcopy(
                MetaHash.objects.get_and_parse(
                    scope=scope, field_definition_scope='fields.field', 
                    clear=True))
            fields.extend(field_hash.values())
            
        for field in fields:
            field['1'] = field['scope']
            field['2'] = field['key']
        
        decorated = [(x['scope'],x['ordinal'],x['key'], x) for x in fields]
        decorated.sort(key=itemgetter(0,1,2))
        fields = [field for scope,ordinal,key,field in decorated]

        return fields
    
    
    @write_authorization
    @un_cache        
    @transaction.atomic    
    def delete_list(self, request, **kwargs):
        MetaHash.objects.all().filter(scope__icontains='fields.').delete()

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
        
        logger.debug('patch_obj: %r: %r', request.user, self._meta.resource_name)
        logger.debug('deserialized: %r', deserialized)
        schema = kwargs.pop('schema', None)
        if not schema:
            raise Exception('schema not initialized')
        fields = schema['fields']
        initializer_dict = {}
        for key in fields.keys():
            if key in deserialized:
                initializer_dict[key] = parse_val(
                    deserialized.get(
                        key,None), key,fields[key].get('data_type','string')) 
        
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        
        try:
            field = None
            try:
                field = MetaHash.objects.get(**id_kwargs)
                errors = self.validate(deserialized, patch=True, schema=schema)
                if errors:
                    raise ValidationError(errors)
            except ObjectDoesNotExist, e:
                logger.debug(
                    'Metahash field %s does not exist, creating', id_kwargs)
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
            logger.debug('save: %r, as %r', deserialized, field)
            field.save()
                    
            logger.debug('patch_obj done')
            return { API_RESULT_OBJ: field }
            
        except Exception, e:
            logger.exception('on patch detail')
            raise e  


class ResourceResource(ApiResource):
    
    class Meta:
        queryset = MetaHash.objects.filter(
            scope="resource").order_by('scope','ordinal','key')
        authentication = MultiAuthentication(
            BasicAuthentication(), IccblSessionAuthentication())
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
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),
            url(r"^(?P<resource_name>%s)/app_data%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_app_data'), name="api_get_app_data"),
            url(r"^(?P<resource_name>%s)/(?P<key>[\w\d_.\-\+: ]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
        ]
        
    def create_fields(self):
        pass
    
    def get_schema(self, request, **kwargs):
        return self.build_response(
            request, self.build_schema(user=request.user), **kwargs)

    def get_app_data(self, request, **kwargs):
        
        app_data = {attr:value for
            attr, value in settings.APP_PUBLIC_DATA.__dict__.iteritems()
                if '__' not in attr }
        
        return self.build_response(
            request, app_data, **kwargs)

    def clear_cache(self):
        ApiResource.clear_cache(self)
        caches['resource_cache'].clear()
        
    def _get_resource_schema(self,resource_key, user):
        ''' For internal callers
        '''
        if DEBUG_RESOURCES:
            logger.info('_get_resource_schema: %r %r...', resource_key, user)
        resources = self._build_resources_internal(user)
        
        if resource_key not in resources:
            raise BadRequest('Resource is not initialized: %r', resource_key)
        
        schema =  resources[resource_key]
        if DEBUG_RESOURCES:
            logger.info('schema fields: %r',
                [(field['key'],field['scope']) 
                    for field in schema['fields'].values()])
        
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

    def build_schema(self, user=None):
        '''
        Override resource method - bootstrap the "Resource" resource schema
        '''
        logger.debug('build_schema for %r: %r', self._meta.resource_name, user)
        resource_fields = self.get_field_resource()._get_list_response_internal(
            scope='fields.resource')
        # build a hash out of the fields
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
            'fields': field_hash
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

            if kwargs.get('meta', None):
                temp = kwargs['meta']
                temp.update(meta)
                meta = temp
            logger.debug('meta: %r', meta)
            
            response_hash = { 
                'meta': meta, 
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

        if user is not None:
            usergroups = set([x.name for x in user.userprofile.get_all_groups()])
        schema = deepcopy(schema)
        fields = deepcopy(schema['fields'])
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
        schema['fields'] = filtered_fields
        
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
            logger.debug('user resource retrieved from cache: %r', user_resources)
        if not user_resources:    
            logger.debug('user resources not cached, build resources')
            if use_cache and self.use_cache:
                resources = resource_cache.get('resources')
                
            if resources:
                if DEBUG_RESOURCES:
                    logger.info('using cached resources')
                    for key,r in resources.items():
                        logger.info('r: %r, fields: %r', key, r['fields'].keys())
                
            else:
                logger.info('rebuilding resources')
            
                resources = deepcopy(
                    MetaHash.objects.get_and_parse(
                        scope='resource', field_definition_scope='fields.resource', 
                        clear=True))
                if not resources:
                    # If there are no resources, use self to bootstrap
                    logger.info('no resources found, using default resource to bootstrap...')
                    resource = self.build_schema(user=user)
                    resources = { resource['key']: resource }
                    
                # get all of the fields
                all_fields = self.get_field_resource()._build_fields()
                field_hash = {}
                # build a hash out of the fields
                for field in all_fields:
                    _fields = field_hash.get(field['scope'],{})
                    _fields[field['key']]=field
                    field_hash[field['scope']] = _fields
                
                # For each resource, pull in the fields of the supertype resource
                for key,resource in resources.items():
                    logger.debug('resource: %r', key)
                    resource['1'] = resource['key']
                    resource['fields'] = field_hash.get('fields.%s'%key, {})
                    resource['resource_uri'] = '/'.join([
                        self._meta.resource_name,resource['key']
                    ])
                    
                    # set the field['table'] 
                    for field in resource['fields'].values():
                        if not field.get('table',None):
                            field['table'] = resource.get('table', None)
                    supertype = resource.get('supertype', None)
                    if supertype:
                        if supertype in resources:
                            logger.debug('find the supertype fields: %r for %r', 
                                supertype, key)
                            supertype_fields = field_hash.get(
                                'fields.%s' % supertype, None)
                            if not supertype_fields:
                                # ok, if the supertype is not built yet
                                if DEBUG_RESOURCES:
                                    logger.warning('no fields for supertype: %r, %r', 
                                        supertype, field_hash.keys())
                            else:
                                # explicitly copy out all supertype fields, then update 
                                # with fields from the current resource
                                inherited_fields = {}
                                for field in supertype_fields.values():
                                    inherited_field = deepcopy(field)
                                    inherited_field['scope'] = \
                                        'fields.%s' % resource['key']
                                    if not inherited_field.get('table',None):
                                        inherited_field['table'] = \
                                            resources[supertype].get('table', None)
                                    
                                    inherited_fields[inherited_field['key']] = \
                                        inherited_field
                                inherited_fields.update(resource['fields'])
                                resource['fields'] = inherited_fields
                        else:
                            logger.error(
                                'supertype: %r, not found in resources: %r', 
                                supertype, resources.keys())
                    
                    update_fields = set()
                    create_fields = set()
                    for key,field in resource['fields'].items():
                        editability = field.get('editability', None)
                        if editability:
                            if 'u' in editability:
                                update_fields.add(key)
                                create_fields.add(key)
                            if 'c' in editability:
                                create_fields.add(key)
                    resource['update_fields'] = list(update_fields)
                    resource['create_fields'] = list(create_fields)
                    
                    resource.get('content_types',[]).append('csv')
                    
                    # extend resource specific data
                    self.extend_specific_data(resource)
                    
                if use_cache and self.use_cache:
                    logger.debug('caching resources')
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
                            'resource: %r, %r', key, resource['fields'].keys())
                if user:
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
            
    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_list(self, request, **kwargs):
        MetaHash.objects.all().filter(scope='resource').delete()

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
        fields = schema['fields']
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
            BasicAuthentication(), IccblSessionAuthentication())
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
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),            
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.\-:]+)/"
                 r"(?P<key>[\w\d_.\-\+:]+)%s$" ) 
                        % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]

    def build_schema(self, user=None):
        schema = super(VocabularyResource,self).build_schema(user=user)
        temp = [ x.scope for x in self.Meta.queryset.distinct('scope')]
        schema['extraSelectorOptions'] = { 
            'label': 'Vocabulary', 'searchColumn': 'scope', 'options': temp }
        return schema
    
    @read_authorization
    def get_detail(self, request, **kwargs):
        key = kwargs.get('key', None)
        if not key:
            raise NotImplementedError('must provide a key parameter')
        
        scope = kwargs.get('scope', None)
        if not scope:
            logger.info('no scope provided')
            raise NotImplementedError('must provide a scope parameter')
        
        kwargs['visibilities'] = kwargs.get('visibilities', ['d'])
        kwargs['is_for_detail']=True
        return self.build_list_response(request, **kwargs)
        
        
    def _get_list_response(self, request, key=None, **kwargs):
        ''' For internal callers '''
        
        # FIXME: check request params for caching
        vocabularies = cache.get('vocabularies')
        if not vocabularies  or not self.use_cache:
            vocabularies =  ApiResource._get_list_response(self, request, **kwargs)
            cache.set('vocabularies', vocabularies)
        
        if key:
            return [vocabularies for vocabularies in vocabularies if vocabularies['key']==key]
        else:
            return vocabularies

        
    @read_authorization
    def get_list(self,request,**kwargs):
 
        kwargs['visibilities'] = kwargs.get('visibilities', ['l'])
        return self.build_list_response(request, **kwargs)
         
#         # FIXME - need to check and cache on kwargs
#         # Ok for the time being, because vocabularies are not filtered
#         cached_content = cache.get('vocabulary_listing')
#         if not cached_content or not self.use_cache:
#             response = self.build_list_response(request, format='json',
#                 visibilities=['l'])
#             cached_content = self._meta.serializer.from_json(
#                 self._meta.serializer.get_content(response))
#             cache.set('vocabulary_listing', cached_content)
#         else:
#             logger.info('using cached vocabularies')
#         desired_format = self.get_format(request, **kwargs)
#         return HttpResponse(
#             content=self.serialize(request, cached_content, desired_format), 
#             content_type=build_content_type(desired_format))
            
        
    def clear_cache(self):
        logger.info('clear vocabulary caching...')
        ApiResource.clear_cache(self)
        cache.delete('vocabulary_listing')
        
    def build_list_response(self,request, **kwargs):
        schema = kwargs.pop('schema', None)
        if not schema:
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
            filename = self._get_filename(readable_filter_hash, schema)
            
            if DEBUG_GET_LIST: 
                logger.info('filter_fields: %r, kwargs: %r', 
                    filter_hash.keys(),kwargs)
            
            
            original_field_hash = schema['fields']
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
            logger.info('VocabularyResource, build field hash...')
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
            logger.info('VocabularyResource, build_sqlalchemy_columns...')
            columns = self.build_sqlalchemy_columns(
                field_hash.values(), base_query_tables=base_query_tables,
                custom_columns=custom_columns )
            j = _vocab
            stmt = select(columns.values()).select_from(j)

            logger.info('VocabularyResource, wrap_statement...')
            # general setup
            (stmt,count_stmt) = self.wrap_statement(
                stmt,order_clauses,filter_expression )
            
            title_function = None
            if use_titles is True:
                def title_function(key):
                    return field_hash[key]['title']
            
            logger.info('vocabularies done, stream response...')
            return self.stream_response_from_statement(
                request, stmt, count_stmt, filename, 
                field_hash=original_field_hash, 
                param_hash=param_hash,
                is_for_detail=is_for_detail,
                rowproxy_generator=json_field_rowproxy_generator,
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get list')
            raise e  
    
    def _get_vocabularies_by_scope(self, scope):
        ''' Utility method
        Retrieve and cache all of the vocabularies in a two level dict
        - keyed by [scope][key]
        '''
        vocabularies = cache.get('vocabularies');
        if not vocabularies:
            vocabularies = {}
            kwargs = {
                'limit': '0'
            }
            _data = self._get_list_response_internal(**kwargs)
            for v in _data:
                _scope = v['scope']
                if _scope not in vocabularies:
                     vocabularies[_scope] = {}
                     logger.debug('created vocab by scope: %r', _scope)
                vocabularies[_scope][v['key']] = v
                
            # Hack: activity.type is serviceactivity.type + activity.class
            # TODO: reinstate this if needed? 20160510
            # if 'serviceactivity.type' in vocabularies:
            #     vocabularies['activity.type'] = \
            #         deepcopy(vocabularies['serviceactivity.type'])
            if 'activity.class' in vocabularies:
                vocabularies['activity.type'].update(
                    deepcopy(vocabularies['activity.class']))
            
            cache.set('vocabularies', vocabularies);
        if scope in vocabularies:
            return deepcopy(vocabularies[scope])
        else:
            logger.warn('---unknown vocabulary scope: %r, %r', scope, vocabularies)
            return {}
    
    def clear_cache(self):
        super(VocabularyResource,self).clear_cache()
        cache.delete('vocabularies');

    @write_authorization
    @un_cache 
    @transaction.atomic       
    def delete_list(self, request, **kwargs):
        Vocabulary.objects.all().delete()

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
        logger.info('VocabularyResource.patch_list done.')
        self.get_debug_times()
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

        fields = schema['fields']
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
            BasicAuthentication(), IccblSessionAuthentication())
        
        # TODO: implement field filtering: users can only see other users details
        # if in readEverythingAdmin group
        authorization = UserResourceAuthorization('user')
        ordering = []
        serializer = LimsSerializer()
        excludes = [] #['json_field']
        always_return_data = True # this makes Backbone happy
        resource_name = 'user'
    
    def clear_cache(self):
        ApiResource.clear_cache(self)
        self.get_resource_resource().clear_cache()
        
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
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),            
            
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))/groups%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_user_groupview'), 
                name="api_dispatch_user_groupview"),
            url(r"^(?P<resource_name>%s)/(?P<username>([\w\d_]+))/permissions%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_user_permissionview'), 
                name="api_dispatch_user_permissionview"),
            ]    

    def dispatch_user_groupview(self, request, **kwargs):
        # signal to include extra column
        return UserGroupResource().dispatch('list', request, **kwargs)    
    
    def dispatch_user_permissionview(self, request, **kwargs):
        # signal to include extra column
        return PermissionResource().dispatch('list', request, **kwargs)    


    def build_schema(self, user=None):
        
        schema = super(UserResource,self).build_schema(user=user)
        try:
            if 'usergroups' in schema['fields']: # may be blank on initiation
                schema['fields']['usergroups']['choices'] = \
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
                    order_by('permission').alias('inner_perms')),
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
                    order_by('permission').alias('innerp')),
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
            raise NotImplementedError('must provide a username parameter')

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
            filename = self._get_filename(readable_filter_hash, schema)
            
            filter_expression = self._meta.authorization.filter(
                request.user, filter_expression)
            
            order_params = param_hash.get('order_by',[])
            field_hash = self.get_visible_fields(
                schema['fields'], filter_hash.keys(), manual_field_includes, 
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
                title_function=title_function, meta=kwargs.get('meta', None),
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
        fields = schema['fields']

        logger.info('patch_obj UserResource: %r', deserialized )
        id_kwargs = self.get_id(deserialized, schema=schema, **kwargs)
        logger.info('patch user: %r', id_kwargs)

        auth_user_fields = { name:val for name,val in fields.items() 
            if val['table'] and val['table']=='auth_user'}
        userprofile_fields = { name:val for name,val in fields.items() 
            if (val['table'] and val['table']=='reports_userprofile'
                or name in ['username', 'ecommons_id'])}
        logger.info('userprofile_fields: %r', userprofile_fields.keys())
        is_patch = False
        try:
            # create the auth_user
            username = id_kwargs.get('username', None)
            ecommons_id = id_kwargs.get('ecommons_id', None)
            if not username:
                logger.info(
                    'username not specified, setting username to ecommons_id: %s', 
                    ecommons_id)
                username = ecommons_id
            
            logger.debug('deserialized: %r', deserialized) 
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
                user.is_active = True
                logger.info('created Auth.User: %s', user)

            initializer_dict = self.parse(
                deserialized, schema=schema, create=not is_patch)
            logger.info('initializer_dict: %r', initializer_dict)

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
                logger.info(
                    'no reports_userprofile fields to update %s', deserialized)

            return { API_RESULT_OBJ: userprofile }
            
        except Exception:
            logger.exception('on patch_obj')
            raise  


class UserGroupResource(ApiResource):
    
    class Meta:
        queryset = UserGroup.objects.all();        
        
        authentication = MultiAuthentication(
            BasicAuthentication(), IccblSessionAuthentication())
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

    def clear_cache(self):
        ApiResource.clear_cache(self)
        self.get_resource_resource().clear_cache()
        
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
        if kwargs.get('data', None):
            # allow for internal data to be passed
            deserialized = kwargs['data']
        else:
            deserialized = self.deserialize(
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
        fields = schema['fields']

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
            logger.info('no group name provided')
            raise NotImplementedError('must provide a group name parameter')

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
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
              
            order_params = param_hash.get('order_by',[])
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
                        .order_by('permission')
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
                        .order_by('permission')
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
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get_list')
            raise e  
        
    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/schema%s$" 
                % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('get_schema'), name="api_get_schema"),            
            
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()),
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/users%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_group_userview'), 
                name="api_dispatch_group_userview"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/permissions%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_group_permissionview'), 
                name="api_dispatch_group_permissionview"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/supergroups%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_group_supergroupview'), 
                name="api_dispatch_group_supergroupview"),
            url(r"^(?P<resource_name>%s)/(?P<name>[^/]+)/subgroups%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
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
            BasicAuthentication(), IccblSessionAuthentication())
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
#         query = self._meta.queryset._clone()
        query = Permission.objects.all()
        permissionTypes = Vocabulary.objects.all().filter(
            scope='permission.type')
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

    def prepend_urls(self):
        return [
            url(r"^(?P<resource_name>%s)/(?P<id>[\d]+)%s$" 
                    % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            url((r"^(?P<resource_name>%s)/(?P<scope>[\w\d_.\-:]+)/"
                 r"(?P<key>[\w\d_.\-\+:]+)/(?P<type>[\w\d_.\-\+:]+)%s$" ) 
                        % (self._meta.resource_name, trailing_slash()), 
                self.wrap_view('dispatch_detail'), name="api_dispatch_detail"),
            ]
    
    @read_authorization
    def get_detail(self, request, **kwargs):

        scope = kwargs.get('scope', None)
        if not scope:
            logger.info('no scope provided')
            raise NotImplementedError('must provide a scope parameter')
        key = kwargs.get('key', None)
        if not key:
            logger.info('no key provided')
            raise NotImplementedError('must provide a key parameter')
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
            filename = self._get_filename(readable_filter_hash, schema)
            filter_expression = \
                self._meta.authorization.filter(request.user,filter_expression)
                  
            order_params = param_hash.get('order_by',[])
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
                title_function=title_function, meta=kwargs.get('meta', None),
                use_caching=True  )
             
        except Exception, e:
            logger.exception('on get list')
            raise e  

    def delete_obj(self, request, deserialized, **kwargs):
        raise NotImplementedError('delete obj is not implemented for Permission')
    
    def patch_obj(self, request, deserialized, **kwargs):
        raise NotImplementedError('patch obj is not implemented for Permission')
    
