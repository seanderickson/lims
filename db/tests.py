from __future__ import unicode_literals

from collections import OrderedDict, defaultdict
from decimal import Decimal
import filecmp
import json
import logging
import os
import re
import sys

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.core.urlresolvers import resolve
from django.db import connection
from django.test import TestCase
from django.test.client import MULTIPART_CONTENT
from django.utils.timezone import now
from tastypie.test import ResourceTestCase, TestApiClient
import xlrd

from db.api import API_MSG_SCREENING_PLATES_UPDATED, \
    API_MSG_SCREENING_ADDED_PLATE_COUNT, API_MSG_SCREENING_DELETED_PLATE_COUNT, \
    API_MSG_SCREENING_EXTANT_PLATE_COUNT, API_MSG_SCREENING_TOTAL_PLATE_COUNT, \
    API_RESULT_DATA, API_RESULT_META, API_MSG_LCPS_INSUFFICIENT_VOLUME, \
    API_MSG_SCP_CREATED, API_MSG_SCP_UNSELECTED, API_MSG_SCP_RESELECTED, \
    API_MSG_LCP_CHANGED, API_MSG_LCP_DESELECTED, API_MSG_LCP_SELECTED, \
    API_MSG_LCP_PLATES_ASSIGNED, API_MSG_LCP_ASSAY_PLATES_CREATED, \
    API_MSG_LCPS_VOLUME_OVERRIDDEN, API_MSG_LCPS_CREATED, \
    API_MSG_PLATING_CANCELED, API_MSG_COPYWELLS_DEALLOCATED, \
    API_MSG_CPR_ASSAY_PLATES_REMOVED, API_MSG_LCPS_REMOVED, \
    API_MSG_LCP_MULTIPLE_SELECTIONS_SUBMITTED, \
    API_MSG_LCPS_MUST_BE_DELETED, API_MSG_CPR_PLATES_PLATED, \
    API_MSG_CPR_PLATES_SCREENED, API_MSG_LCPS_UNFULFILLED,\
    API_PARAM_SET_DESELECTED_TO_ZERO, API_MSG_SCPS_DELETED,\
    API_MSG_COPYWELLS_ALLOCATED, API_MSG_COPYWELLS_DEALLOCATED,\
    LibraryCopyPlateResource, LibraryScreeningResource
from db.api import API_PARAM_SHOW_OTHER_REAGENTS, API_PARAM_SHOW_COPY_WELLS, \
    API_PARAM_SHOW_RETIRED_COPY_WELlS, API_PARAM_VOLUME_OVERRIDE
from db.api import VOCAB_LCP_STATUS_SELECTED, VOCAB_LCP_STATUS_UNFULFILLED, \
    VOCAB_LCP_STATUS_PLATED

import db.api
from db.models import Reagent, Substance, Library, ScreensaverUser, \
    UserChecklist, AttachedFile, ServiceActivity, Screen, Well, Publication, \
    PlateLocation, LibraryScreening, LabAffiliation
import db.models
from db.support import lims_utils, screen_result_importer, bin_packer
from db.test.factories import LibraryFactory, ScreenFactory, \
    ScreensaverUserFactory, LabAffiliationFactory
from reports import ValidationError, HEADER_APILOG_COMMENT, _now, \
    API_RESULT_ERROR
from reports.api import API_MSG_COMMENTS, API_MSG_CREATED, \
    API_MSG_SUBMIT_COUNT, API_MSG_UNCHANGED, API_MSG_UPDATED, \
    API_MSG_ACTION, API_MSG_RESULT, API_MSG_WARNING, API_MSG_NOT_ALLOWED, \
    API_RESULT_META, API_PARAM_OVERRIDE
from reports.models import ApiLog, UserProfile, UserGroup, API_ACTION_PATCH, \
    API_ACTION_CREATE, Vocabulary
from reports.serialize import XLSX_MIMETYPE, SDF_MIMETYPE, JSON_MIMETYPE
from reports.serializers import CSVSerializer, XLSSerializer, LimsSerializer, \
    ScreenResultSerializer
from reports.tests import IResourceTestCase, equivocal
from reports.tests import assert_obj1_to_obj2, find_all_obj_in_list, \
    find_obj_in_list, find_in_dict
import unittest
from __builtin__ import False
import decimal


logger = logging.getLogger(__name__)


BASE_URI = '/db/api/v1'
BASE_REPORTS_URI = '/reports/api/v1'
import db; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__))
BASE_URI_DB = '/db/api/v1'

LCP_COPYWELL_KEY = db.api.LabCherryPickResource.LCP_COPYWELL_KEY


class DBResourceTestCase(IResourceTestCase):
    """
    Invoke _bootstrap_init_files on setup
    """
    def __init__(self,*args,**kwargs):
    
        super(DBResourceTestCase, self).__init__(*args,**kwargs)
        self.directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')

    def _bootstrap_init_files(self, reinit_pattern=None):
        logger.info( 'bootstrap reinit_pattern: %r', reinit_pattern)
        super(DBResourceTestCase, self)._bootstrap_init_files(
            reinit_pattern=reinit_pattern)
        self.directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        serializer=CSVSerializer() 
        testApiClient = TestApiClient(serializer=serializer) 
        input_actions_file = os.path.join(
            self.directory, 'api_init_actions.csv')
        logger.info('open input_actions file: %r', input_actions_file)
        self._run_api_init_actions(
            input_actions_file, reinit_pattern=reinit_pattern)
    
    def create_library(self, data=None):
        ''' Create a test Library through the API'''

        input_data = LibraryFactory.attributes()
        if data:
            input_data.update(data)
            
        logger.info('create library: %r', input_data)
        resource_uri = '/'.join([BASE_URI_DB, 'library'])
        test_uri = '/'.join([resource_uri,input_data['short_name']])
        return self._create_resource(input_data,resource_uri,test_uri)

    @staticmethod
    def create_small_molecule_test_well(
        plate, well_index, platesize=384, **kwargs):
        ''' Generate a test well for library initialization'''
        
        library_well_types = [
            'empty','experimental','dmso','library_control' ]
        well_name = lims_utils.well_name_from_index(well_index, platesize)
        library_well_type = kwargs.get(
            'library_well_type',library_well_types[well_index%3])
        _data = {
            'plate_number': plate, 
            'well_name': well_name,
            'well_id': '%s:%s' % (str(plate).zfill(5),well_name),
            'library_well_type' : library_well_type,
        }
        if library_well_type == 'experimental':
            _data['molar_concentration'] = '0.00%d' % (well_index+1)
            _data['vendor_name'] = 'vendorX'
            _data['vendor_identifier'] = 'ID-%d' % well_index
            _data['vendor_batch_id'] = 2
            _data['smiles'] = \
                'H%dC%dN%d' % (well_index%10,well_index%11,well_index%12)
            _data['compound_name'] = \
                ['name-%d'%well_index, 'name-%d'%((well_index+11)%100)]
    
        for k,v in kwargs.items():
            _data[k] = v
        return _data

    @staticmethod
    def create_test_well(
        plate, well_index, platesize=384, **kwargs):
        ''' Generate a test well for library initialization'''
        
        library_well_types = [
            'empty','experimental','dmso','library_control' ]
        well_name = lims_utils.well_name_from_index(well_index, platesize)
        library_well_type = kwargs.get(
            'library_well_type',library_well_types[well_index%3])
        _data = {
            'plate_number': plate, 
            'well_name': well_name,
            'well_id': '%s:%s' % (str(plate).zfill(5),well_name),
            'library_well_type' : library_well_type,
        }
        if library_well_type == 'experimental':
            _data['molar_concentration'] = '0.00%d' % (well_index+1)
            _data['vendor_name'] = 'vendorX'
            _data['vendor_identifier'] = 'ID-%d' % well_index
            _data['vendor_batch_id'] = 2
    
        for k,v in kwargs.items():
            _data[k] = v
        return _data

    def create_screen(self, data=None, uri_params=None, resource_uri=None):
        ''' Create a test Screen through the API'''
        
        
        input_data = ScreenFactory.attributes()
        logger.info('input_data: %r', input_data)
        if data:
            input_data.update(data)
        if 'lab_head_id' not in input_data:
            lab_head = self.create_lab_head()
            input_data['lab_head_id'] = str(lab_head['screensaver_user_id'])
        if 'lead_screener_id' not in input_data:
            lead_screener = self.create_screensaveruser()
            input_data['lead_screener_id'] = str(lead_screener['screensaver_user_id'])
        if 'collaborator_ids' not in input_data:
            collaborator1 = self.create_screensaveruser()
            collaborator2 = self.create_screensaveruser()
            input_data['collaborator_ids'] = [
                collaborator1['screensaver_user_id'], collaborator2['screensaver_user_id']]
        input_data['collaborator_ids'] = [ 
            str(x) for x in input_data['collaborator_ids']]
        if resource_uri is None:
            resource_uri = '/'.join([BASE_URI_DB, 'screen'])
        
        if uri_params is not None:
            resource_uri += '?' + '&'.join(uri_params)
        _data_for_get = { 
            'limit': 0,
            'includes': '*'
        }
        logger.info('screen input_data to create: %r', input_data)
        logger.info('post to %r...', resource_uri)
        resp = self.api_client.post(
            resource_uri, format='json', data=input_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code, self.get_content(resp)))
    
        new_obj = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in new_obj)
        self.assertEqual(len(new_obj[API_RESULT_DATA]), 1)        
        new_obj = new_obj[API_RESULT_DATA][0]
        new_obj = self.get_screen(new_obj['facility_id'])  
        logger.debug('screen created: %r', new_obj)
        result,msg = assert_obj1_to_obj2(input_data,new_obj)
        self.assertTrue(result, msg)
        return new_obj
    
    def create_copy(self, copy_input_data ):   
        resource_uri = BASE_URI_DB + '/librarycopy'
        resource_test_uri = '/'.join([
            resource_uri,copy_input_data['library_short_name'],
            copy_input_data['copy_name']])
        copy_data = self._create_resource(
            copy_input_data, resource_uri, resource_test_uri, 
            excludes=['initial_plate_well_volume','initial_plate_status'])
        logger.info('created library copy: %r', copy_data)
        return copy_data
    
    def get_screen(self, facility_id, data_for_get=None):
        ''' Retrieve a Screen from the API'''
        resource_uri = '/'.join([BASE_URI_DB, 'screen', facility_id])
        return self.get_single_resource(resource_uri, data_for_get)

    def get_schema(self, resource_name, authentication=None, data=None):
        
        data_for_get = { 'includes': '*' }
        if data is not None:
            data_for_get.update(data)
        
        if authentication is None:
            authentication=self.get_credentials()
        
        resource_uri = '/'.join([BASE_URI_DB, resource_name, 'schema'])

        resp = self.api_client.get(
            resource_uri, format='json', authentication=authentication,
            data=data_for_get)
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        return self.deserialize(resp)

    def create_screensaveruser(self, data=None):
        ''' Create a test ScreensaverUser through the API'''
        
        input_data = ScreensaverUserFactory.attributes()
        if data:
            input_data.update(data)
        resource_uri = '/'.join([BASE_URI_DB, 'screensaveruser'])
        test_uri = '/'.join([resource_uri,input_data['username']])
        logger.info('create user: %r', input_data)
        return self._create_resource(input_data,resource_uri,test_uri)
    
    def create_lab_head(self, data=None):

        lab_head_data = self.create_screensaveruser(data=data)
        logger.info('lab_head_data: %r', lab_head_data)
        return self.update_user_to_lab_head(lab_head_data)
    
    def update_user_to_lab_head(self, lab_head_data):
        lab_affiliation = self.create_lab_affiliation()
        
        user_patch_data = {
            'screensaver_user_id': lab_head_data['screensaver_user_id'],
            'classification': 'principal_investigator',
            'lab_affiliation_id': lab_affiliation['lab_affiliation_id']
            }
        
        resource_uri = BASE_URI_DB + '/screensaveruser/'
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=user_patch_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        new_lab_head = _data[API_RESULT_DATA][0]
        logger.info('the new lab affiliation has been set: %r', new_lab_head)
        
        self.assertEqual(
            lab_affiliation['name'], 
            new_lab_head['lab_affiliation_name'])
        
        return new_lab_head
    
    def create_lab_affiliation(self, data=None):
        attributes = LabAffiliationFactory.attributes()
        if data is not None:
            attributes.update(data)
        
        resource_uri = '/'.join([
            BASE_REPORTS_URI, 'vocabulary','labaffiliation.category',attributes['category']])
        lab_affiliation_category = self.get_single_resource(resource_uri)
        if lab_affiliation_category is None:
            lab_affiliation_category = {
                'scope': 'labaffiliation.category',
                'key': attributes['category'],
                'ordinal': attributes['ordinal'],
                'title': 'Lab Affiliation Category ' + attributes['category'],
                'description': 'Lab Affiliation Category desc: ' + attributes['category']}
            
            resp = self.api_client.post(
                resource_uri, 
                format='json', 
                data=lab_affiliation_category, 
                authentication=self.get_credentials())
            self.assertTrue(
                resp.status_code in [200,201,202], 
                (resp.status_code, self.get_content(resp)))
            _data = self.deserialize(resp)
            self.assertTrue(API_RESULT_DATA in _data)
            self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
            new_affiliation_category = _data[API_RESULT_DATA][0]
        
            logger.info('created category: %r', new_affiliation_category)
            for key,val in lab_affiliation_category.items():
                self.assertEqual(
                    lab_affiliation_category[key],new_affiliation_category[key])
        
        # create a lab_affiliation
        resource_uri = BASE_URI_DB + '/labaffiliation/'        
        lab_affiliation = {
            'category': lab_affiliation_category['key'],
            'name': attributes['name'],
        }
        resp = self.api_client.post(
            resource_uri, 
            format='json', 
            data=lab_affiliation, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        new_lab_affiliation =_data[API_RESULT_DATA][0]
        logger.info('created lab: %r', new_lab_affiliation)
        
        self.assertTrue('lab_affiliation_id' in new_lab_affiliation, 
            'no id created for lab_affiliation: %r' % new_lab_affiliation)
        
        for key,val in lab_affiliation.items():
            self.assertEqual(lab_affiliation[key],new_lab_affiliation[key])

        return new_lab_affiliation
    

def setUpModule():

    logger.info('=== setup Module')
    keepdb = False
    reinit_metahash = False
    reinit_pattern = None
    
    if len(sys.argv) > 1:
        for i,arg in enumerate(sys.argv):
            logger.info('arg: %d: %r',i, arg)
            if 'keepdb' in arg:
                keepdb = True
            if 'reinit_metahash' in arg:
                reinit_metahash = True
            if 'reinit_pattern' in arg:
                # grab the next arg
                reinit_pattern = sys.argv[i+1]
    # Set up a superuser
    try:
        logger.info('create/find superuser %s...', IResourceTestCase.username)
        IResourceTestCase.user = User.objects.get(
            username=IResourceTestCase.username)
        logger.info('superuser found: %r', IResourceTestCase.user)
        logger.info('users: %r', [str(u) for u in User.objects.all()])
    except ObjectDoesNotExist:
        logger.info('creating superuser: %s', IResourceTestCase.username)
        IResourceTestCase.user = User.objects.create_superuser(
            IResourceTestCase.username, '1testsuperuser@example.com', 
            IResourceTestCase.password)
        logger.info('superuser created.')

    if reinit_metahash or not keepdb:
        testContext = DBResourceTestCase(methodName='_bootstrap_init_files')
        testContext.setUp()
        testContext._bootstrap_init_files(reinit_pattern=reinit_pattern)
    else:
        print 'skip database metahash initialization when using keepdb'

    try:
        su = ScreensaverUser.objects.get(username=IResourceTestCase.username)
        logger.info('found ss user: %r', su)
        temp_test_case = DBResourceTestCase(methodName='get_from_server')
        resource_uri = \
            BASE_URI_DB + '/screensaveruser/' +IResourceTestCase.username
        DBResourceTestCase.admin_user = \
            temp_test_case.get_from_server(resource_uri)
        logger.info('got admin user: %r', DBResourceTestCase.admin_user)
    except Exception:
        logger.exception('create an admin screensaveruser...')
        temp_test_case = DBResourceTestCase(methodName='create_screensaveruser')
        DBResourceTestCase.admin_user = temp_test_case.create_screensaveruser({ 
            'username': temp_test_case.username,
            'is_superuser': True
        })
        logger.info('admin screensaveruser created')
    
    logger.info('=== setup Module done')

def tearDownModule():

    logger.info('=== teardown Module')
#     Vocabulary.objects.all().filter(scope__contains='labaffiliation.').delete()
#     LabAffiliation.objects.all().delete()

    # remove the admin user
    # ScreensaverUser.objects.all().delete() 
    # UserProfile.objects.all().delete()
    # User.objects.all().delete()

class LibraryResource(DBResourceTestCase):

    def setUp(self):

        super(LibraryResource, self).setUp()

    def tearDown(self):

        DBResourceTestCase.tearDown(self)
        logger.info('delete library resources')
        Library.objects.all().delete()
        # removing the library should remove dependent resources
        # Well.objects.all().delete()
        # Reagent.objects.all().delete()
        PlateLocation.objects.all().delete()
        ApiLog.objects.all().delete()
    
    def test_a_plate_search_parser(self):
        
        test_search_1 = '''
        1000
        1000-2000 A
        B 3000-4000
        5000-6000 A,B,C     
        9000,2000,3000 5000-4000 A,"b-c",D,"Stock A"
        '''
        
        # NOTE that sublists should also be sorted in ascending order,
        # so [A,'b-c',D] becomes [A,D,b-c] (sort capital letters first)
        expected_parsed = [
            { 'plate': 1000 },
            { 'plate_range': [1000,2000], 'copy': 'A' },
            { 'plate_range': [3000,4000], 'copy': 'B' },
            { 'plate_range': [5000,6000], 'copy': ['A','B','C'] },
            { 'plate': [2000,3000,9000], 'plate_range': [4000,5000], 
                'copy': ['A','D', 'b-c','Stock A'] },
        ]
        
        parsed_searches = \
            LibraryCopyPlateResource.parse_plate_copy_search(test_search_1)
        
        # extract single values from lists for convenience
        for parsed_search in parsed_searches:
            for k,v in parsed_search.items():
                if len(v)==1:
                    parsed_search[k] = v[0]
        
        not_found = []
        for expected_search in expected_parsed:
            found = False
            for parsed_search in parsed_searches:
                for k,v in expected_search.items():
                    if parsed_search[k] != v:
                        break
                found = True    
            if found is not True:
                not_found.append(expected_search)
        
        self.assertTrue(len(not_found)==0,
            'expected_searches %r not found in %r' % (not_found, parsed_searches))
        
    def test1_create_library(self):

        logger.info('test1_create_library ...')
        
        resource_uri = BASE_URI_DB + '/library'
        
        library1 = self.create_library({ 
            'short_name': 'testlibrary1','start_plate': '1534', 
            'end_plate': '1534', 'plate_size': '384' })
        library2 = self.create_library({ 
            'short_name': 'testlibrary2','start_plate': '1535', 
            'end_plate': '1537', 'plate_size': '384' })
        
        logger.info('Find the <undefined> library2 wells that were created...')
        resource_uri = '/'.join([
            BASE_URI_DB,'library',library2['short_name'],'well'])
        data_for_get={ 'limit': 0, 'includes': ['*', '-structure_image'] }
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        self.assertTrue(resp['Content-Type'].startswith('application/json'))
        new_obj = self.deserialize(resp)
        expected_count = 384*3
        self.assertEqual(
            len(new_obj[API_RESULT_DATA]), expected_count, 
            'wrong number of wells: %d, expected: %d' 
                % (len(new_obj[API_RESULT_DATA]), expected_count))
        
        index = 0
        platesize = 384
        plate = 1535        
        substance_ids = set()
        # Examine wells - first plate only for speed
        for j in range(384):
            well_name = lims_utils.well_name_from_index(j, platesize)
            well_id = lims_utils.well_id(plate,well_name)
            well_search = {'well_id': well_id}
            result, well = find_obj_in_list(
                well_search, new_obj[API_RESULT_DATA])
            self.assertTrue(result, well)
            
        logger.info('Load wells to the library1...')
        plate = 1534
        input_data = [
            self.create_small_molecule_test_well(plate,i) 
                for i in range(0,384)]
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library1['short_name'],resource_name])
        resp = self.api_client.put(
            resource_uri, format='sdf', data={ 'objects': input_data } , 
            authentication=self.get_credentials(), **data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # Examine Well/Reagents created:
        # NOTE: the well resource is not the linked resource type, and does not
        # have the reagent fields.  
        resource_name = 'reagent'
        reagent_resource_uri = '/'.join([
            BASE_URI_DB,'library', library1['short_name'],resource_name])
        returned_data = \
            self.get_list_resource(reagent_resource_uri, data_for_get)
        expected_count = 384
        self.assertEqual(
            len(returned_data), expected_count, 
            ('library: %r'% library1, 'expected', 
                expected_count, 'found',len(returned_data)))
        
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        self.validate_wells(input_data, returned_data, fields)
        
    def test1a_create_library_comments(self):

        logger.info('test1a_create_library_comments ...')
        
        resource_uri = BASE_URI_DB + '/library'
        
        library1 = self.create_library({ 
            'short_name': 'testlibrary1','start_plate': '1534', 
            'end_plate': '1534', 'plate_size': '384' })
        
        patch_uri = '/'.join([resource_uri,library1['short_name']]) 
        
        test_comment = 'Some test comment 123 xyz'
        
        header_data = { HEADER_APILOG_COMMENT: test_comment}
        
        resp = self.api_client.patch(
            patch_uri, format='json', #data={ 'objects': [ copywell_input] }, 
            authentication=self.get_credentials(), **header_data )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        patch_response = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in patch_response)
        self.assertEqual(len(patch_response[API_RESULT_DATA]), 1)        
        patch_response = patch_response[API_RESULT_DATA][0]
        self.assertTrue('comment_array' in patch_response, 
            'patch_response: %r' % patch_response)
        comment_array = patch_response['comment_array']
        self.assertTrue(len(comment_array),1)
        self.assertTrue(test_comment in patch_response['comment_array'][0], 
            'test_comment: %r not found in library response obj: %r' % (
                test_comment, patch_response))

        test_comment2 = 'another test comment...'
        
        header_data = { HEADER_APILOG_COMMENT: test_comment2}
        
        resp = self.api_client.patch(
            patch_uri, format='json', #data={ 'objects': [ copywell_input] }, 
            authentication=self.get_credentials(), **header_data )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        patch_response = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in patch_response)
        self.assertEqual(len(patch_response[API_RESULT_DATA]), 1)        
        patch_response = patch_response[API_RESULT_DATA][0]
        self.assertTrue('comment_array' in patch_response, 
            'patch_response: %r' % patch_response)
        comment_array = patch_response['comment_array']
        self.assertTrue(len(comment_array),2)
        self.assertTrue(test_comment in comment_array[1], 
            'test_comment: %r not found in library response obj: %r' % (
                test_comment, patch_response))
        self.assertTrue(test_comment2 in comment_array[0], 
            'test_comment2: %r not found in library response obj: %r' % (
                test_comment2, patch_response))
        
    def validate_wells(self, input_data, output_data, fields):
        ''' 
        Validate that the input well/reagents were created in the output_data
        '''
        substance_ids = set()
        for inputobj in input_data:
            # 1. Validate well/reagents exist
            search = { 'well_id': 
                lims_utils.well_id(
                    inputobj['plate_number'],inputobj['well_name']) }
            result, outputobj = find_obj_in_list(
                search,output_data) #, excludes=excludes )
            self.assertTrue(
                result, 
                ('not found', search,outputobj,'=== objects returned ===', 
                      output_data ) ) 
            # 2. Validate well/reagent fields
            expected_data = { 
                key: inputobj[key] for key in fields.keys() if key in inputobj}
            result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
            self.assertTrue(
                result, (msgs, 'input', expected_data, 'output', outputobj))

            substance_id = outputobj['substance_id']
            self.assertTrue(
                substance_id not in substance_ids, 
                ('substance_id not unique', substance_id, substance_ids))
            substance_ids.add(substance_id)
        
    def test10_create_library_copy_specific_wells(self):

        logger.info('test10_create_library_copy ...')
        
        # 1. Create Library
        logger.info('create library a library...')
        start_plate = 1000
        end_plate = 1005
        plate_size = 384
        library_data = self.create_library({
            'start_plate': start_plate, 
            'end_plate': end_plate,
            'plate_size': plate_size,
            'screen_type': 'small_molecule' })
        short_name = library_data['short_name']
        
        # 2. Create Library Wells: wells have various concentrations 
        logger.info('create and load well data, plates: %r-%r...', 
            start_plate, end_plate)
        well_input_data = {}
        for plate in range(start_plate, end_plate+1):
            for i in range(0,plate_size):
                # Create test wells, concentrations will be varied
                well_input = self.create_small_molecule_test_well(
                    plate,i,library_well_type='experimental')
                well_input_data[well_input['well_id']] = well_input
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', short_name, resource_name])
        resp = self.api_client.put(
            resource_uri, format='sdf', 
            data={ 'objects': well_input_data.values()}, 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 3. Create a Library Copy:
        # - will also create Plates: 
        # -- initial_well_volume will be set
        # -- because the well concentrations vary, 
        # Plate.concentration fields are not set
        logger.info('Create library copy...')
        copy_input_data = {
            'library_short_name': library_data['short_name'],
            'copy_name': "A",
            'usage_type': "library_screening_plates",
            'initial_plate_well_volume': Decimal('0.000040')
        
        }        
        resource_uri = BASE_URI_DB + '/librarycopy'
        resource_test_uri = '/'.join([
            resource_uri,copy_input_data['library_short_name'],
            copy_input_data['copy_name']])
        library_copy = self._create_resource(
            copy_input_data, resource_uri, resource_test_uri, 
            excludes=['initial_plate_well_volume'])

        # 4. Verify created Plates
        logger.info('Verify plates created...')
        uri = '/'.join([resource_test_uri,'plate'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        start_plate = int(library_data['start_plate'])
        end_plate = int(library_data['end_plate'])
        number_of_plates = end_plate-start_plate+1
        self.assertEqual(len(new_obj[API_RESULT_DATA]),number_of_plates)
        for plate_data in new_obj[API_RESULT_DATA]:
            self.assertEqual(library_copy['copy_name'],plate_data['copy_name'])
            plate_number = int(plate_data['plate_number'])
            self.assertTrue(
                plate_number>=start_plate and plate_number<=end_plate,
                'wrong plate_number: %r' % plate_data)
            self.assertEqual(
                Decimal(copy_input_data['initial_plate_well_volume']),
                Decimal(plate_data['well_volume']))
            self.assertTrue(plate_data['molar_concentration'] is None)
            self.assertTrue(plate_data['mg_ml_concentration'] is None)
            # Min Molar concentration is set, because copy_wells were created
            self.assertTrue(plate_data['min_molar_concentration'] is not None)
            # Min mg/ml concentration not set, as copy_wells only have molar
            self.assertTrue(plate_data['min_mg_ml_concentration'] is None)
            self.assertEqual(plate_data['status'], 'not_specified')
            
        # 5. Verify created CopyWells
        logger.info('Verify copy_wells created (check varying concentrations)')
        uri = '/'.join([
            BASE_URI_DB,
            'library',
            copy_input_data['library_short_name'],
            'copy',
            copy_input_data['copy_name'],
            'copywell'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        copy_wells_returned = new_obj[API_RESULT_DATA]
        self.assertEqual(len(copy_wells_returned),len(well_input_data))

        for copy_well_data in copy_wells_returned:
            self.assertEqual(
                library_copy['copy_name'],copy_well_data['copy_name'])
            self.assertTrue(
                copy_well_data['plate_number'] >= start_plate
                and copy_well_data['plate_number'] <= end_plate,
                'copy_well returned for wrong plate: %r' % copy_well_data )
            well_input = well_input_data.get(copy_well_data['well_id'], None)
            self.assertTrue(well_input is not None, 
                'copy well not in wells: %r' % copy_well_data)
            if well_input['library_well_type'] == 'experimental':
                self.assertEqual(
                    Decimal(copy_input_data['initial_plate_well_volume']),
                    Decimal(copy_well_data['initial_volume']))
                self.assertEqual(
                    copy_well_data['initial_volume'], 
                    copy_well_data['volume'])
                self.assertEqual(
                    Decimal(well_input['molar_concentration']),
                    Decimal(copy_well_data['molar_concentration']), 
                    'molar concentration mismatch: %r output %r' 
                    % (well_input,copy_well_data))
                self.assertTrue(copy_well_data['mg_ml_concentration'] is None)
            else:
                self.assertTrue(
                    copy_well_data['initial_volume'] is None,
                    'non-experimental well has data: %r' % copy_well_data)
                self.assertTrue(copy_well_data['volume'] is None,
                    'non-experimental well has data: %r' % copy_well_data)
                self.assertTrue(copy_well_data['molar_concentration'] is None,
                    'non-experimental well has data: %r' % copy_well_data)
        return (library_data, library_copy, copy_wells_returned)

    def test10a_create_library_copy_simple_wells(self):

        logger.info('test10a_create_library_copy_simple_wells ...')
        logger.info('creates a library wells with single concentration')
        
        # 1. Create a Library
        logger.info('create library a library...')
        start_plate = 1000
        end_plate = 1005
        plate_size = 384
        library_data = self.create_library({
            'start_plate': start_plate, 
            'end_plate': end_plate,
            'plate_size': plate_size,
            'screen_type': 'small_molecule' })
        short_name = library_data['short_name']
        
        # 2. Create Library Wells: wells all have the same concentration 
        logger.info('create and load well data, plates: %r-%r...', 
            start_plate, end_plate)
        well_input_data = {}
        single_molar_concentration = '0.010'
        for plate in range(start_plate, end_plate+1):
            for i in range(0,plate_size):
                well_input = self.create_small_molecule_test_well(
                    plate,i)
                if well_input['library_well_type'] == 'experimental':
                    well_input['molar_concentration'] = single_molar_concentration
                well_input_data[well_input['well_id']] = well_input
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', short_name, resource_name])
        resp = self.api_client.put(
            resource_uri, format='sdf', 
            data={ 'objects': well_input_data.values()}, 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        # 3. Create a Library Copy        
        logger.info('create library copy...')
        copy_input_data = {
            'library_short_name': library_data['short_name'],
            'copy_name': "A",
            'usage_type': "library_screening_plates",
            'initial_plate_well_volume': '0.000040'
        }        
        resource_uri = BASE_URI_DB + '/librarycopy'
        resource_test_uri = '/'.join([
            resource_uri,copy_input_data['library_short_name'],
            copy_input_data['copy_name']])
        library_copy = self._create_resource(
            copy_input_data, resource_uri, resource_test_uri, 
            excludes=['initial_plate_well_volume'])
        
        # 4. Verify created Plates
        logger.info('Verify plates created...')
        uri = '/'.join([resource_test_uri,'plate'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        start_plate = int(library_data['start_plate'])
        end_plate = int(library_data['end_plate'])
        number_of_plates = end_plate-start_plate+1
        self.assertEqual(len(new_obj[API_RESULT_DATA]),number_of_plates)
        for plate_data in new_obj[API_RESULT_DATA]:
            self.assertEqual(library_copy['copy_name'],plate_data['copy_name'])
            plate_number = int(plate_data['plate_number'])
            self.assertTrue(
                plate_number>=start_plate and plate_number<=end_plate,
                'wrong plate_number: %r' % plate_data)
            self.assertEqual(
                Decimal(copy_input_data['initial_plate_well_volume']),
                Decimal(plate_data['well_volume']))
            self.assertEqual(
                Decimal(plate_data['molar_concentration']),
                Decimal(single_molar_concentration))
            self.assertEqual(
                plate_data['molar_concentration'],
                plate_data['min_molar_concentration'])
            self.assertTrue(plate_data['mg_ml_concentration'] is None)
            self.assertTrue(plate_data['min_mg_ml_concentration'] is None)

        # 5. Verify CopyWell data can be queried, 
        # but that all have single concentration
        logger.info('Verify copy_wells created ()...')
        uri = '/'.join([
            BASE_URI_DB,
            'library',
            copy_input_data['library_short_name'],
            'copy',
            copy_input_data['copy_name'],
            'copywell'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        copy_wells_returned = new_obj[API_RESULT_DATA]
        self.assertEqual(len(copy_wells_returned),len(well_input_data))

        for copy_well_data in copy_wells_returned:
            self.assertEqual(
                library_copy['copy_name'],copy_well_data['copy_name'])
            self.assertTrue(
                copy_well_data['plate_number'] >= start_plate
                and copy_well_data['plate_number'] <= end_plate,
                'copy_well returned for wrong plate: %r' % copy_well_data )
            well_input = well_input_data.get(copy_well_data['well_id'], None)
            self.assertTrue(well_input is not None, 
                'copy well not in wells: %r' % copy_well_data)
            if well_input['library_well_type'] == 'experimental':
                self.assertEqual(
                    Decimal(copy_input_data['initial_plate_well_volume']),
                    Decimal(copy_well_data['initial_volume']))
                self.assertEqual(
                    Decimal(copy_well_data['initial_volume']), 
                    Decimal(copy_well_data['volume']))
                self.assertEqual(
                    Decimal(well_input['molar_concentration']),
                    Decimal(copy_well_data['molar_concentration']), 
                    'molar concentration mismatch: %r output %r' 
                        %(well_input,copy_well_data))
                self.assertTrue(copy_well_data['mg_ml_concentration'] is None)
            else:
                self.assertTrue(
                    copy_well_data['initial_volume'] is None,
                    'non-experimental well has data: %r' % copy_well_data)
                self.assertTrue(copy_well_data['volume'] is None,
                    'non-experimental well has data: %r' % copy_well_data)
                self.assertTrue(copy_well_data['molar_concentration'] is None,
                    'non-experimental well has data: %r' % copy_well_data)
        return (library_data, library_copy, copy_wells_returned)
        
    def test11_create_library_copy_invalids(self):
        # TODO: try to create duplicate copy names
        # TODO: try to create invalid types, plate ranges
        
        pass
    
    def test12_modify_copy_wells(self):

        logger.info('test12_modify_copy_wells ...')
    
        (library_data, copy_data, plate_data) = \
            self.test10_create_library_copy_specific_wells()
        end_plate = int(library_data['end_plate'])
        start_plate = int(library_data['start_plate'])
        short_name = library_data['short_name']
        plate_size = int(library_data['plate_size'])

        logger.info('retrieve the copy_wells...')
        resource_uri = '/'.join([
            BASE_URI_DB,'library',library_data['short_name'],'copy',
            copy_data['copy_name'],'copywell'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        data_for_get['plate_number__eq'] = start_plate
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        copywell_data = new_obj[API_RESULT_DATA]
        expected_wells = plate_size;
        logger.info('returned copywell: %r', copywell_data[0])
        logger.info('returned copywell: %r', copywell_data[-1])
        self.assertEqual(len(copywell_data),expected_wells)
        
        logger.info('patch a copywell...')
        
        copywell_input = copywell_data[0]
        original_volume = Decimal(copywell_input['volume'])
        volume_adjustment = Decimal('0.000008')
        copywell_input['volume'] = original_volume-volume_adjustment
        
        patch_uri = '/'.join([resource_uri,copywell_input['well_id']]) 
        resp = self.api_client.patch(
            patch_uri, format='json', data=copywell_input, 
            authentication=self.get_credentials(), **data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        patch_response = self.deserialize(resp)
        
        self.assertTrue(API_RESULT_DATA in patch_response)
        self.assertEqual(len(patch_response[API_RESULT_DATA]), 1)        
        new_copywell = patch_response[API_RESULT_DATA][0]
        self.assertEqual(
            volume_adjustment, Decimal(new_copywell['consumed_volume']))
        self.assertEqual(
            Decimal(copywell_input['volume']), Decimal(new_copywell['volume']))
    
    def test13_plate_locations(self):

        logger.info('test13_plate_locations ...')
        (library_data, copy_data, plate_data) = \
            self.test10_create_library_copy_specific_wells()
        
        end_plate = library_data['end_plate']
        start_plate = library_data['start_plate']
        short_name = library_data['short_name']
        
        # 1. Simple test
        lps_format = '{library_short_name}:{copy_name}:{start_plate}-{end_plate}'
        copy_plate_ranges = [
            lps_format.format(
                library_short_name=library_data['short_name'],
                copy_name=copy_data['copy_name'],
                start_plate=start_plate,
                end_plate=end_plate,)
        ]
        
        plate_location_input = {
            'room': 'room1', 'freezer': 'freezer1', 'shelf': 'shelf1',
            'bin': 'bin1',
            'copy_plate_ranges': copy_plate_ranges,
        }
        
        resource_uri = BASE_URI_DB + '/platelocation'
        resource_test_uri = resource_uri
        plate_location = self._create_resource(
            plate_location_input, resource_uri, resource_test_uri,)
        self.assertEqual(
            plate_location['copy_plate_ranges'],
            plate_location_input['copy_plate_ranges'],
            'plate ranges not equal: %r - %r' % (
                plate_location['copy_plate_ranges'],
                plate_location_input['copy_plate_ranges']))

        # Verify that plates have the location set as well
        resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'plate_number__range': [start_plate,end_plate],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        plates_data = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        location_fields = ['room','freezer','shelf','bin']
        for plate_data in plates_data:
            for field in location_fields:
                self.assertTrue(
                    plate_data[field]==plate_location_input[field],
                    'plate location: expected: %r, rcvd: %r'
                    % (plate_location_input, plate_data))
        
        # Test apilogs
        resource_uri = BASE_REPORTS_URI + '/apilog'
        
        # Test platelocation apilog
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), 
            data={ 
                'limit': 0, 
                'ref_resource_name': 'platelocation'})
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        expected_count = 1 # create
        self.assertEqual( 
            len(new_obj[API_RESULT_DATA]), expected_count , 
            str((len(new_obj[API_RESULT_DATA]), expected_count, new_obj)))
        
        # Test librarycopyplate apilogs
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), 
            data={ 
                'limit': 0, 
                'ref_resource_name': 'librarycopyplate'})
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        expected_count = 6 # one for each plate in the copy_range
        self.assertEqual( 
            len(new_obj[API_RESULT_DATA]), expected_count , 
            str((len(new_obj[API_RESULT_DATA]), expected_count, new_obj)))
        for logvalue in new_obj[API_RESULT_DATA]:
            diffs = json.loads(logvalue['diffs'])
            self.assertTrue(diffs['bin']==[None, 'bin1'],
                'wrong diff: %r' % diffs )

        # 2. remove plates from the range
        copy_plate_ranges = [
            lps_format.format(
                library_short_name=library_data['short_name'],
                copy_name=copy_data['copy_name'],
                start_plate=start_plate,
                end_plate=end_plate-2,)
        ]
        plate_location_input = {
            'room': 'room1', 'freezer': 'freezer1', 'shelf': 'shelf1',
            'bin': 'bin1',
            'copy_plate_ranges': copy_plate_ranges,
        }
        resource_uri = BASE_URI_DB + '/platelocation'
        
        # TODO: test multipart/form upload, which is what the client will use
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data={'objects': [plate_location_input],}, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        # TODO: check for PATCH meta information
        
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 'includes': ['copy_plate_ranges'],}, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('resp: %r', new_obj)
        new_obj = new_obj[API_RESULT_DATA][0]
        self.assertEqual(
            new_obj['copy_plate_ranges'],
            plate_location_input['copy_plate_ranges'],
            'plate ranges not equal: %r - %r' % (
                new_obj['copy_plate_ranges'],
                plate_location_input['copy_plate_ranges']))
        
        # Verify that removed plates have a location set to null
        resource_uri = BASE_URI_DB + '/librarycopyplate'
        data_for_get={ 
            'plate_number__range': [end_plate-1,end_plate],
            'copy_name__eq': copy_data['copy_name']
            } 
        resp = self.api_client.get(
            resource_uri,format='json', data=data_for_get,
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        self.assertTrue(len(new_obj[API_RESULT_DATA])==2, 
            'could not find all modified plates for %r, %r'
            %(data_for_get, new_obj))
        for plate_data in new_obj[API_RESULT_DATA]:
            for field in location_fields:
                self.assertTrue(
                    plate_data[field]==None,
                    'plate location: %r should be None, %r'
                    % (field, plate_data))
    
    def test17_batch_edit_copyplate_info(self):
        
        logger.info('test17_batch_edit_copyplate_info ...')
        
        (library_data, copy_data, plate_data) = \
            self.test10_create_library_copy_specific_wells()
        end_plate = library_data['end_plate']
        start_plate = library_data['start_plate']
        short_name = library_data['short_name']
        
        # 1. POST new data to the copyplates
        # NOTE: the librarycopyplate batch_edit uses a POST form to send both
        # the search data (3 types of search filter: "nested_search_data",  
        # "raw_search_data", and GET search params),
        # as well as the update data (in the form of 
        # "plate_info" and "plate_location")
        plate_info = {
            'status': 'available',
            'remaining_well_volume': '0.000030', 
            'plate_type': 'nunc_96'
        }
        data = {
            'plate_info': json.dumps(plate_info) , 
            'nested_search_data': json.dumps({
                'library_short_name': library_data['short_name'],
                'copy_name': copy_data['copy_name'] })
            }
        # NOTE: the tastypie client does not correctly encode the multipart
        # form data for a POST, so the django client is used
        data_for_get = {'HTTP_AUTHORIZATION': self.get_credentials()}
        data_for_get['HTTP_ACCEPT'] = JSON_MIMETYPE
        resource_uri = BASE_URI_DB + '/librarycopyplate/batch_edit'
        resp = self.django_client.post(
            resource_uri, content_type=MULTIPART_CONTENT, 
            data=data, 
            **data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        post_response = self.deserialize(resp)
        
        # Inspect meta "Result" section
        self.assertTrue(
            API_RESULT_META in post_response, '%r' % post_response)
        meta = post_response[API_RESULT_META]
        self.assertTrue(API_MSG_RESULT in meta, '%r' % post_response)
        self.assertTrue(
            API_MSG_SUBMIT_COUNT in meta[API_MSG_RESULT], 
            '%r' % post_response)
        self.assertTrue(
            meta[API_MSG_RESULT][API_MSG_SUBMIT_COUNT]==6, 
            'Wrong "%r" count: %r' 
                % (API_MSG_SUBMIT_COUNT, meta))
        logger.info('post_response: %r', post_response)
        
        # 2. Verify plates
        resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'library_short_name__eq': library_data['short_name'],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        plates_data = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        
        plates_data_output = plates_data
        for plate_data_output in plates_data_output:
            for field in plate_info.keys():
                self.assertTrue(
                    equivocal(plate_data_output[field],plate_info[field]),
                    'plate data: expected: %r, rcvd: %r'
                    % (plate_info[field],plate_data_output[field]))
            self.assertEqual(
                _now().date().strftime("%Y-%m-%d"), 
                plate_data_output['date_plated'],
                'expected date_plated: %r, %r' 
                    %(_now().date(), plate_data_output['date_plated']))

        # 3. POST new data - set the status to retired
        plate_info['status'] = 'retired'
        resource_uri = BASE_URI_DB + '/librarycopyplate/batch_edit'
        
        data = {
            'plate_info': json.dumps(plate_info) , 
            'nested_search_data': json.dumps({
                'library_short_name': library_data['short_name'],
                'copy_name': copy_data['copy_name'] })
            }
        data_for_get = {'HTTP_AUTHORIZATION': self.get_credentials()}
        data_for_get['HTTP_ACCEPT'] = JSON_MIMETYPE
        resp = self.django_client.post(
            resource_uri, content_type=MULTIPART_CONTENT, 
            data=data, 
            **data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        post_response = self.deserialize(resp)
        
        # Inspect meta "Result" section
        self.assertTrue(
            API_RESULT_META in post_response, '%r' % post_response)
        meta = post_response[API_RESULT_META]
        self.assertTrue(API_MSG_RESULT in meta, '%r' % post_response)
        self.assertTrue(
            API_MSG_SUBMIT_COUNT in meta[API_MSG_RESULT], '%r' % post_response)
        self.assertTrue(meta[API_MSG_RESULT][API_MSG_SUBMIT_COUNT]==6, 
            'Wrong "%r" count: %r' 
            % (API_MSG_SUBMIT_COUNT, meta))
        logger.info('post_response: %r', post_response)
        
        # 4. Verify plates status changed to "retired"
        resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'library_short_name__eq': library_data['short_name'],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        plates_data = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        
        plates_data_output = plates_data    
        for plate_data_output in plates_data_output:
            self.assertEqual('retired', plate_data_output['status'])
            self.assertEqual(
                _now().date().strftime("%Y-%m-%d"), 
                plate_data_output['date_retired'],
                'expected date_plated: %r, %r' 
                    %(_now().date(), plate_data_output['date_plated']))
    
    def test16_batch_edit_copyplate_location(self):
        
        logger.info('test16_batch_edit_copyplate_location ...')

        (library_data, copy_data, plate_data) = \
            self.test10_create_library_copy_specific_wells()
        end_plate = library_data['end_plate']
        start_plate = library_data['start_plate']
        short_name = library_data['short_name']
        
        plate_location_input = {
            'room': 'room1', 'freezer': 'freezer1', 'shelf': 'shelf1',
            'bin': 'bin1',
        }
        # NOTE: the librarycopyplate batch_edit uses a POST form to send both
        # the search data (3 types of search filter: "nested_search_data",  
        # "raw_search_data", and GET search params),
        # as well as the update data (in the form of 
        # "plate_info" and "plate_location")
        data = {
            'plate_location': json.dumps(plate_location_input) , 
            'nested_search_data': json.dumps({
                'library_short_name': library_data['short_name'],
                'copy_name': copy_data['copy_name'] })
            }

        logger.info('Patch batch_edit: data: %r', data)

        # NOTE: the tastypie client does not correctly encode the multipart
        # form data for a POST, so the django client is used
        data_for_get = {'HTTP_AUTHORIZATION': self.get_credentials()}
        data_for_get['HTTP_ACCEPT'] = JSON_MIMETYPE
        resource_uri = BASE_URI_DB + '/librarycopyplate/batch_edit'
        resp = self.django_client.post(
            resource_uri, content_type=MULTIPART_CONTENT, 
            data=data, 
            **data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        post_response = self.deserialize(resp)
        
        # Inspect meta "Result" section
        self.assertTrue(
            API_RESULT_META in post_response, '%r' % post_response)
        meta = post_response[API_RESULT_META]
        self.assertTrue(API_MSG_RESULT in meta, '%r' % post_response)
        self.assertTrue(
            API_MSG_SUBMIT_COUNT in meta[API_MSG_RESULT], '%r' % post_response)
        self.assertTrue(meta[API_MSG_RESULT][API_MSG_SUBMIT_COUNT]==6, 
            'Wrong %r count: %r' 
            % (API_MSG_SUBMIT_COUNT,meta))
        
        # Get plates as defined
        resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'library_short_name__eq': library_data['short_name'],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        plates_data = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        
        plates_data_input = plates_data
        for plate_data in plates_data_input:
            for field in plate_location_input.keys():
                self.assertTrue(
                    plate_data[field]==plate_location_input[field],
                    'plate location: expected: %r, rcvd: %r'
                    % (plate_location_input[field],plate_data))

        # check ApiLogs
        # one for each plate, one parent_log, one for each location
        resource_uri = BASE_REPORTS_URI + '/apilog'
        
        # Test platelocation apilog
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), 
            data={ 
                'limit': 0, 
                'ref_resource_name': 'platelocation'})
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        expected_count = 1 # created
        self.assertEqual( 
            len(new_obj[API_RESULT_DATA]), expected_count , 
            str((len(new_obj[API_RESULT_DATA]), expected_count, new_obj)))
        
        # Test librarycopyplate apilogs
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), 
            data={ 
                'limit': 0, 
                'ref_resource_name': 'librarycopyplate',
                'includes': ['added_keys']})
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        expected_count = 7 # one for each plate in the copy_range, one parent log
        self.assertEqual( 
            len(new_obj[API_RESULT_DATA]), expected_count , 
            str((len(new_obj[API_RESULT_DATA]), expected_count, new_obj)))
        for logvalue in new_obj[API_RESULT_DATA]:
            if logvalue['parent_log_uri']:
                diffs = json.loads(logvalue['diffs'])
                self.assertTrue(diffs['bin']==[None, 'bin1'],
                    'wrong diff: %r' % diffs )
            else:
                logger.info('parent log: %r', logvalue)
                self.assertTrue(logvalue['key']=='librarycopyplate',
                    'parent_log key should be "librarycopyplate", %r' % logvalue)
    
    def test15_modify_copyplate_info(self):

        logger.info('test15_modify_copyplate_info ...')
        
        (library_data, copy_data, plate_data) = \
            self.test10_create_library_copy_specific_wells()
        end_plate = library_data['end_plate']
        start_plate = library_data['start_plate']
        short_name = library_data['short_name']
        
        new_plate_data = {
            'status': 'available',
            'remaining_well_volume': '0.000030', 
            'plate_type': 'nunc_96'
        }
        # 1. Get plates as defined
        resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'plate_number__range': [start_plate,end_plate],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        plates_data = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        
        plates_data_input = plates_data
        for plate_data in plates_data_input:
            plate_data.update(new_plate_data)

        # 2. Patch the plates
        logger.info('Patch the plates: cred: %r', self.username)
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data={'objects': plates_data_input,}, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        patch_response = self.deserialize(resp)
        self.assertTrue(
            API_RESULT_META in patch_response, '%r' % patch_response)
        meta = patch_response[API_RESULT_META]
        self.assertTrue(API_MSG_RESULT in meta, '%r' % patch_response)
        self.assertTrue(
            API_MSG_SUBMIT_COUNT in meta[API_MSG_RESULT], 
            '%r' % patch_response)
        self.assertEqual(
            meta[API_MSG_RESULT][API_MSG_SUBMIT_COUNT],
            (end_plate-start_plate+1),
            '"%r" : %r, expected: %r' 
            % (API_MSG_SUBMIT_COUNT,
                meta, (end_plate-start_plate+1)))
        
        # 3. Verify that the plates have the expected location
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'plate_number__range': [start_plate,end_plate],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        plates_data_output = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        for plate_data in plates_data_output:
            for field in new_plate_data.keys():
                self.assertTrue(
                    equivocal(plate_data[field],new_plate_data[field]),
                    'plate data: expected: %r, rcvd: %r'
                    % (new_plate_data[field],plate_data[field]))
                self.assertEqual(
                    _now().date().strftime("%Y-%m-%d"), 
                    plate_data['date_plated'],
                    'expected date_plated: %r, %r' 
                        %(_now().date(), plate_data['date_plated']))

        # Test ApiLogs:
        # plate - one for each plate
        # plate_location - one for each plate addition to the range
        
    def test14_modify_copy_plate_locations(self):

        logger.info('test14_modify_copy_plate_locations ...')
        
        (library_data, copy_data, plate_data) = \
            self.test10_create_library_copy_specific_wells()
        end_plate = library_data['end_plate']
        start_plate = library_data['start_plate']
        short_name = library_data['short_name']

        plate_location_input = {
            'room': 'Room1', 'freezer': 'Freezer1', 'shelf': 'Shelf1', 
            'bin': 'Bin1' }

        # Get plates as defined
        resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'plate_number__range': [start_plate,end_plate],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        plates_data = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        
        plates_data_input = plates_data
        for plate_data in plates_data_input:
            for field in plate_location_input.keys():
                self.assertTrue(
                    plate_data[field]==None,
                    'plate location: expected: None, rcvd: %r'
                    % (plate_data))
                plate_data[field] = plate_location_input[field]

        # Patch the plates
        logger.info('Patch the plates: cred: %r', self.username)
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data={'objects': plates_data_input,}, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        patch_response = self.deserialize(resp)
        self.assertTrue(
            API_RESULT_META in patch_response, '%r' % patch_response)
        meta = patch_response[API_RESULT_META]
        self.assertTrue(API_MSG_RESULT in meta, '%r' % patch_response)
        self.assertTrue(
            API_MSG_SUBMIT_COUNT in meta[API_MSG_RESULT], 
            '%r' % patch_response)
        self.assertEqual(
            meta[API_MSG_RESULT][API_MSG_SUBMIT_COUNT],
            (end_plate-start_plate+1),
            '"%r" : %r, expected: %r' 
            % (API_MSG_SUBMIT_COUNT,
                meta, (end_plate-start_plate+1)))
        
        # Verify that the plates have the expected location
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 
                'plate_number__range': [start_plate,end_plate],
                'copy_name__eq': copy_data['copy_name']
                }, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        plates_data = new_obj[API_RESULT_DATA]
        self.assertTrue(len(plates_data)==end_plate-start_plate+1, 
            'could not find all of the plates in the range: %r, %r'
            % ([start_plate,end_plate],plates_data))
        for plate_data in plates_data:
            for field in plate_location_input.keys():
                self.assertTrue(
                    plate_data[field]==plate_location_input[field],
                    'plate location: expected: None, rcvd: %r'
                    % (plate_location_input))

        # Verify that the plate range is set on the location:
        lps_format = '{library_short_name}:{copy_name}:{start_plate}-{end_plate}'
        expected_copy_plate_ranges = [
            lps_format.format(
                library_short_name=library_data['short_name'],
                copy_name=copy_data['copy_name'],
                start_plate=start_plate,
                end_plate=end_plate)]
        
        plate_location_input['copy_plate_ranges'] = expected_copy_plate_ranges
        resource_uri = BASE_URI_DB + '/platelocation'
        resp = self.api_client.get(
            resource_uri,format='json', 
            data={ 'includes': ['copy_plate_ranges'],}, 
            authentication=self.get_credentials(),)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        self.assertTrue(
            len(new_obj[API_RESULT_DATA])==1, 'too many locations were created')
        new_obj = new_obj[API_RESULT_DATA][0]
        self.assertEqual(
            new_obj['copy_plate_ranges'],
            plate_location_input['copy_plate_ranges'],
            'plate ranges not equal: %r - %r' % (
                new_obj['copy_plate_ranges'],
                plate_location_input['copy_plate_ranges']))
        
        # Test ApiLogs:
        # plate - one for each plate
        # plate_location - one for each plate addition to the range
        
    def test4_create_library_invalids(self):
        ''' Test Library schema "required" validations'''
         
        logger.info('test create_library_invalids...')
        library_resource = self.get_resource_from_server('library')
        fields = library_resource['fields']
        resource_uri = BASE_URI_DB + '/library'
        for key,field in fields.items():
            if field.get('required', False):
                logger.debug('testing required field: %r, %r', key, field)
                
                library_item = LibraryFactory.attributes()
                library_item[key] = None
                resp = self.api_client.post(
                    resource_uri, format='json', data=library_item, 
                    authentication=self.get_credentials())
                self.assertTrue(
                    resp.status_code in [400], 
                    (resp.status_code, self.get_content(resp)))
        
                data = self.deserialize(resp)
                logger.info('response: %r', data) 
                data = data[API_RESULT_ERROR]
                self.assertTrue(find_in_dict(key, data), 
                    'Error: response error not found: %r, obj: %r' %(key, data))

        # Test invalid Library name                
        library_item = LibraryFactory.attributes()
        library_item['library_name'] = 'invalid & name'
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        
        logger.debug('response.content.library message: %r', 
                         getattr(resp, 'content'))
        
        key = 'library_name'
        data = self.deserialize(resp)
        logger.info('response: %r', data) 
        data = data[API_RESULT_ERROR]
        self.assertTrue(find_in_dict(key, data), 
            'Error: response error not found: %r, obj: %r' %(key, data))

        # Test invalid Library type
        library_item = LibraryFactory.attributes()
        library_item['library_type'] = 'invalid_type'
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        key = 'library_type'
        data = self.deserialize(resp)
        logger.info('response: %r', data) 
        data = data[API_RESULT_ERROR]
        self.assertTrue(find_in_dict(key, data), 
            'Error: response error not found: %r, obj: %r' %(key, data))

        # TODO: test regex and number: min/max, vocabularies
        
    def test61_small_molecule_delete_compound_name(self):
        # TODO: test that sm.compound names can be removed
        pass
    
    def test6_load_small_molecule_file(self):
        
        logger.info('test6_load_small_molecule_file')
        
        library_item = self.create_library({ 
            'start_plate': '1536', 
            'end_plate': '1536', 
            'plate_size': '384',
            'screen_type': 'small_molecule' })

        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        filename = (
            '%s/db/static/test_data/libraries/clean_data_small_molecule.sdf'
                % APP_ROOT_DIR )

        data_for_get = { 'limit': 0, 'includes': ['*'] }
        data_for_get['HTTP_ACCEPT'] = SDF_MIMETYPE

        logger.info('Open and PUT file: %r', filename)
        with open(filename) as input_file:
            
            input_data = self.serializer.from_sdf(input_file.read())
            input_data = input_data['objects']
            expected_count = 8
            self.assertEqual(
                len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            logger.debug('======Submitting patch file: %r, uri: %r ...', 
                filename, resource_uri)
            resp = self.api_client.put(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
        
        # Examine Well/Reagents created:
        # NOTE: get the library "reagents" instead of the wells: reagents have 
        # one-to-many reln to wells, so the well endpoint returns only first;
        # also, the well resource is not the linked resource type, and does not
        # have the reagent fields.  
        resource_name = 'reagent'
        reagent_resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        resp = self.api_client.get(
            reagent_resource_uri, format='sdf', 
            authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        returned_data = new_obj[API_RESULT_DATA]
        expected_count = 384
        self.assertEqual(
            len(returned_data), expected_count, 
            ('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count,
                'returned_data',returned_data))

        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        self.validate_wells(input_data, returned_data, fields)
                    
        # Test 2: update some wells, check results, and api logs
        filename = ( APP_ROOT_DIR 
            + '/db/static/test_data/libraries/clean_data_small_molecule_update.sdf')

        data_for_get = { 'limit': 0, 'includes': ['*', '-structure_image'] }
        data_for_get['HTTP_ACCEPT'] = SDF_MIMETYPE

        logger.info('Open and PATCH file: %r', filename)
        with open(filename) as input_file:
            
            input_data = self.serializer.from_sdf(input_file.read())
            input_data = input_data[API_RESULT_DATA]
            expected_count = 4
            self.assertEqual(
                len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            logger.info('PATCH: %r ...', resource_uri)
            resp = self.api_client.patch(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
        
            resource_name = 'well'
            resource_uri = '/'.join([
                BASE_URI_DB,'library', library_item['short_name'],
                resource_name])
            resp = self.api_client.get(
                reagent_resource_uri, format='sdf', 
                authentication=self.get_credentials(), 
                data=data_for_get)
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            returned_data = new_obj[API_RESULT_DATA]
            expected_count = 384
            self.assertEqual(
                len(returned_data), expected_count, 
                ('returned_data of ',filename,'found',
                    len(returned_data), 'expected',expected_count,
                    'returned_data',returned_data))
        
            # 1. test the specific wells
            well_ids_to_check = input_data
    
            for update_well in well_ids_to_check:
                search = { 
                    'well_name': update_well['well_name'], 
                    'plate_number': update_well['plate_number']}
                result, outputobj = find_obj_in_list(
                    search,returned_data) #, excludes=excludes )
                self.assertTrue(
                    result, 
                    ('not found', search,outputobj,'=== objects returned ===', 
                          returned_data ) ) 
                result, msgs = assert_obj1_to_obj2(update_well, outputobj)
                self.assertTrue(result, (msgs, update_well, outputobj))
                self.assertTrue(
                    'library_well_type' in outputobj,
                    'library_well_type is missing in %r' % outputobj)
                self.assertTrue(
                    outputobj['library_well_type']=='experimental',
                    'wrong library_well_type: %r' % outputobj)
            
            # 2. check the apilogs - library
            resource_uri = BASE_REPORTS_URI + '/apilog'
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 
                    'limit': 0, 
                    'ref_resource_name': 'library', 
                    'key': library_item['short_name'] })
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            expected_count = 3 # create, post, update
            self.assertEqual( 
                len(new_obj[API_RESULT_DATA]), expected_count , 
                str((len(new_obj[API_RESULT_DATA]), expected_count, new_obj)))

            # 3. check the apilogs - wells
            resource_uri = BASE_REPORTS_URI + '/apilog' 
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 'limit': 0, 'ref_resource_name': 'well' })
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            logs = new_obj[API_RESULT_DATA]
            # (none for create), 4 for update
            expected_count = 4
            self.assertEqual( 
                len(logs), expected_count , 
                str((len(logs), expected_count)))
            
            self.assertEqual(
                len(logs), 4, 
                ('expected %d patch logs, found: %d' %(4,len(logs)), logs))
            
            for log in logs:
                self.assertTrue(
                    ( 'parent_log_uri' in log and
                      library_item['short_name'] in log['parent_log_uri'] ), 
                    'parent_log_uri should contain: %r - %r' 
                    % (library_item['short_name'], log) )
                if log['key'] == '01536:A01':
                    logger.info('log: %r', log)
                    self.assertTrue(
                        set(log['diff_keys'])==
                        set([
                            "pubchem_cid", "vendor_identifier", 
                            "vendor_batch_id", "compound_name", "smiles"]),
                        'diff_keys: %r should equal %r' % (
                            log['diff_keys'],[
                                "pubchem_cid", "vendor_identifier", 
                                "vendor_batch_id", "compound_name", "smiles"]))
            # TODO: check parent_log - library log/ version
    
    def test7_load_sirnai(self):

        logger.info('test7_load_sirnai ...')
        
        filename = ('%s/db/static/test_data/libraries/clean_data_rnai.xlsx'
                    % APP_ROOT_DIR )
        
        library_item = self.create_library({ 
            'start_plate': 50001,  
            'end_plate': 50001, 
            'plate_size': '384',
            'screen_type': 'rnai' })
        resource_uri = BASE_URI_DB + '/library'
        
        self._load_xls_reagent_file(filename,library_item, 5,384 )

    # def test7a_sirna_validations(self):
    #     # TODO: test validations
    #     pass

    def test8_sirna_duplex(self):
        
        logger.info('test8_sirna_duplex ...')
        
        filename = (  
            '%s/db/static/test_data/libraries/clean_rnai_duplex_50440_50443.xlsx'
            % APP_ROOT_DIR)
        
        # create the duplex library
        library_item = self.create_library({ 
            'start_plate': 50440,  
            'end_plate': 50443, 
            'plate_size': '384',
            'screen_type': 'rnai' })
        
        self._load_xls_reagent_file(filename,library_item, 532, 1536 )
        filename = ( '%s/db/static/test_data/libraries/clean_rnai_pool.xlsx'
                     % APP_ROOT_DIR )
        # create the pool library
        library_item = self.create_library({ 
            'start_plate': 50439,  
            'end_plate': 50439, 
            'plate_size': '384',
            'screen_type': 'rnai',
            'is_pool': True })
        
        self._load_xls_reagent_file(filename,library_item, 133, 384 )
        
        ## TODO: test the duplex wells are set on the pool well
        
    def test9_natural_product(self):
        
        logger.info('test9_natural_product ...')
        
        filename = (  
            '%s/db/static/test_data/libraries/clean_data_natural_product.xlsx'
            % APP_ROOT_DIR )
        
        library_item = self.create_library({ 
            'start_plate': 2037,  
            'end_plate': 2037, 
            'plate_size': '384',
            'library_type': 'natural_products' })
        
        self._load_xls_reagent_file(filename,library_item, 352, 384 )

    def _load_xls_reagent_file(
            self,filename,library_item, expected_in, expected_count):
        ''' Test the loading of an xls well/reagent file''' 

        start_plate = library_item['start_plate']
        end_plate = library_item['end_plate']

        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        
        data_for_get = { 'limit': 0, 'includes': ['*', '-structure_image'] }
        data_for_get['HTTP_ACCEPT'] = XLSX_MIMETYPE
        xls_serializer = XLSSerializer()
        
        logger.info('Open and PUT file: %r', filename)
        with open(filename) as input_file:
            ### NOTE: submit the data as an object, because the test framework
            ### will convert it to the target format.
            input_data = self.serializer.from_xlsx(input_file.read())
            input_data = [x for x in input_data[API_RESULT_DATA]]
            self.assertEqual(
                len(input_data), expected_in, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_in)))
            
            resp = self.api_client.put(
                resource_uri, format='xlsx', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(
                resp.status_code in [200, 204], 
                (resp.status_code, 
                 self.get_content(resp)))
        
        resource_name = 'reagent'
        reagent_resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        resp = self.api_client.get(
            reagent_resource_uri, format='xlsx', 
            authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        returned_data = [x for x in new_obj[API_RESULT_DATA]]
        self.assertEqual(
            len(returned_data), expected_count, 
            str(('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count )))

        # 1. test well keys
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        
        self.validate_wells(input_data, returned_data, fields)
        
        # TODO: test substance IDs    
        # substance_ids = set()
        # for inputobj in input_data:
        #     # first just search for the well
        #     search = { 
        #         'well_id': lims_utils.well_id(
        #             inputobj['plate_number'],inputobj['well_name']) }
        #     result, outputobj = find_obj_in_list(search,returned_data)
        #     self.assertTrue(
        #         result, 
        #         ('not found', search,outputobj,'=== objects returned ===', 
        #               returned_data )) 
        # 
        #     # second look at the found item
        #     expected_data = { key: inputobj[key] for key in fields.keys() 
        #         if key in inputobj }
        #     logger.debug('checking: expected: %r', expected_data )
        #     result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
        #     self.assertTrue(result, (msgs, expected_data, outputobj))
        #      
        #     substance_id = outputobj['substance_id']
        #     self.assertTrue(
        #         substance_id not in substance_ids, 
        #         ('substance_id not unique', substance_id, substance_ids))
        #     substance_ids.add(substance_id)
                    

class ScreenResultSerializerTest(TestCase):
    ''' Test serialization of the ScreenResult file'''
    
    def test1_read_write(self):
        
        serializer = ScreenResultSerializer()
        file = 'ScreenResultTest_1_valid.xlsx'
        file_out = 'ScreenResultTest_1_valid_out.xlsx'
        filename = '%s/db/static/test_data/screens/%s' %(APP_ROOT_DIR,file)
        filename_out = (
            '%s/db/static/test_data/screens/%s' %(APP_ROOT_DIR,file_out) )
        
        logger.info('read test file: %r', filename)
        with open(filename, 'rb') as input_file:
            wb = xlrd.open_workbook(file_contents=input_file.read())
            input_data = screen_result_importer.read_workbook(wb)
            logger.info('input_data: keys: %r', input_data.keys())
            
        lookup_key = 'meta'
        self.assertTrue(lookup_key in input_data, 
            '%r not in input_data: %r' % (lookup_key, input_data.keys()))
        meta = input_data[lookup_key]
        lookup_key = 'screen_facility_id'
        self.assertTrue(lookup_key in meta, 
            '%r not in meta: %r' % (lookup_key, meta))
        screen_facility_id = meta[lookup_key]
        
        lookup_key = 'fields'
        self.assertTrue(lookup_key in input_data, 
            '%r not in input_data: %r' % (lookup_key, input_data.keys()))
        data_columns = input_data[lookup_key]
        logger.info('data columns: %r', data_columns.keys())
        # TODO: spot check
        
        lookup_key = API_RESULT_DATA
        self.assertTrue(lookup_key in input_data, 
            '%r not in input_data: %r' % (lookup_key, input_data.keys()))
        result_values = input_data[lookup_key]
        with open(filename_out, 'wb') as output_file:
            logger.info('write back out to %r', filename_out)
            data = screen_result_importer.create_output_data(
                screen_facility_id, data_columns, result_values)
            output_file.write(serializer.to_xlsx(data))
                
        with open(filename_out, 'rb') as output_file_in:
            try:
                logger.info('read back in from output file: %r', filename_out)
                output_data = serializer.from_xlsx(output_file_in.read())
                logger.info('output_data: keys: %r', output_data.keys())
            except ValidationError,e:
                logger.warn('validation error: %r', e.errors)
                raise
        self.validate(self, input_data, output_data)
        
    def test2_errors(self):
        ''' Test ScreenResult file format errors'''
        
        serializer = ScreenResultSerializer()
        file = 'ScreenResultParseErrorsTest.xlsx'
        file_out = 'ScreenResultParseErrorsTest_out.xlsx'
        filename = '%s/db/static/test_data/screens/%s' %(APP_ROOT_DIR,file)
        filename_out = (
            '%s/db/static/test_data/screens/%s' %(APP_ROOT_DIR,file_out) )
        logger.info('check file: %r ', os.path.exists(filename))
        logger.info('Open and validate file: %r', filename)
        with open(filename, 'rb') as input_file:
            try:
                input_data = serializer.from_xlsx(input_file.read())
                
                # Iterate some of the result values to trigger the 
                # validation error
                for rv in input_data['objects']:
                    logger.info('rv: %r', rv)
                    self.fail(
                        'error workbook should generate errors exception: %r'
                        % filename)
            except ValidationError, e:
                logger.info('(Expected) errors reported: %r', e.errors)
                self.assertTrue(
                    len(e.errors.keys())==3, 
                    'should be 3 errors, one for each plate sheet without '
                    'a "plate" field')
                
    @staticmethod
    def validate(testinstance,input_data,output_data):
        ''' Validate ScreenResult serialization round trip '''
        
        logger.info('validate input/output data')
        logger.info('input_data: keys: %r', input_data.keys())
        logger.info('output_data: keys: %r', output_data.keys())
        logger.info(
            'input_data: field keys: %r', input_data['fields'].keys())
        logger.info(
            'output_data: field keys: %r', output_data['fields'].keys())
                    
        # Test DataColumns
        input_fields = input_data['fields']
        output_fields = output_data['fields']
        testinstance.assertTrue(
            len(input_fields)==len(output_fields),
            'input/output fields: %r, %r' 
            % (input_fields.keys(), output_fields.keys()))    
        for colname,input_field in input_fields.items():
            found = None
            for output_field in output_fields.values():
                if input_field['name'] == output_field['name']:
                    found = output_field
                    break
            testinstance.assertTrue(
                found is not None, 
                'input datacolumn not found: %r, output datacolumns: %r'
                % (input_field, output_fields))

            for column_field, val in input_field.items():
                val = input_field[column_field]
                val2 = output_field[column_field]
                if val and column_field == 'derived_from_columns':
                    val = set([x.upper() for x in re.split(r'[,\s]+', val)])
                    val2 = set([x.upper() for x in re.split(r'[,\s]+', val2)])
                
                if column_field == 'data_type' and val == 'numeric':
                    if input_field.get('decimal_places', 0) > 0:
                        testinstance.assertTrue(val2 in ['numeric', 'decimal'],
                            'unrecognized type conversion: %r to %r for %r'
                            % (val, val2, colname))
                        if val2 != 'decimal':
                            logger.warn(
                            'col: %r: data_type: "numeric" with decimal_places > 0 '
                            'should be migrated to '
                            '"decimal": val2: %r',colname, val2)
                    else:
                        
                        testinstance.assertTrue(val2 in ['numeric', 'integer'],
                            'unrecognized type conversion: %r to %r for %r'
                            % (val, val2, colname))
                        if val2 != 'integer':
                            # Note 201707: "numeric" type is migrated to 
                            # "decimal" and "integer" by the api; some tests
                            # do not roundtrip the data to the api
                            logger.warn(
                            'col: %r: data_type: "numeric" with decimal_places not > 0 '
                            'should be migrated to '
                            '"integer": val2: %r', colname, val2)
                elif column_field == 'data_type' and val == 'text':
                    testinstance.assertTrue(val2 in ['text','string'],
                        'unrecognized type conversion: %r to %r for %r'
                        % (val, val2, colname))
                else:
                    result,msg = equivocal(val, val2)
                    testinstance.assertTrue(
                        result,
                        'column not equal: column %r:%r, input: %r, output: %r, %r'
                        % (colname,column_field, val, val2, msg))
        
        # test ResultValues
        items_to_test = 1000
        try:
            # save output data generator to a list
            output_list = [x for x in output_data['objects']]
        except Exception, e:
            logger.exception('%r', e)
            raise
        for i,input_rv in enumerate(input_data['objects']):
            found_row = None
            for output_rv in output_list:
                logger.debug('try: %r - %r', input_rv,output_rv)
                if input_rv['well_id'] == output_rv['well_id']:
                    found_row = output_rv
                    break;
                logger.debug('not matched')
            if not found_row:
                testinstance.fail(
                    'result value well_id not found: %r in %r' 
                    % (input_rv['well_id'], [x for x in output_list]))
            
            logger.debug('input_row: %r', input_rv)
            logger.debug('found_row: %r', found_row)
            for key,val in input_rv.items():
                if val is not None:
                    # NOTE: 20170724 - if input values are null, ResultValues are not created
                    testinstance.assertTrue(
                        key in found_row, 
                        ('key: %r, from input: %r,'
                            ' not found in output input_rv: %r') 
                        % (key, input_rv, found_row) )
                else:
                    if key not in found_row:
                        logger.info(
                            'FYI: Null input %r:%r not found in ouput row: %r',
                            key, val, found_row)
                        continue
                val2 = found_row[key]
#                 if val2 == 'NP':
#                     testinstance.assertTrue(val=='NP' or val==None,
#                         ('partition positive val: %r, key: %r, %r - %r'
#                             % (val, key, input_rv, found_row)))
#                 elif val2 == 'NT':
#                     testinstance.assertTrue(val=='NT' or val==None,
#                         ('confirmed positive val: %r, key: %r, %r - %r'
#                             % (val, key, input_rv, found_row)))
#                 else:
                result,msg = equivocal(val, val2)
                testinstance.assertTrue(
                    result,
                    ('meta field not equal: %r %r != %r, %r'
                         'input: %r, output: %r')
                    % (key, val, val2, msg, input_rv, found_row))
            if i == items_to_test:
                break
                
class ScreenResultResource(DBResourceTestCase):

    def __init__(self,*args,**kwargs):
    
        super(DBResourceTestCase, self).__init__(*args,**kwargs)
        self.sr_serializer = ScreenResultSerializer()
        self.sr_api_client = TestApiClient(serializer=self.sr_serializer)       

    def setUp(self):

        super(ScreenResultResource, self).setUp()

    def tearDown(self):
        
        DBResourceTestCase.tearDown(self)
        logger.info('delete resources')
        with connection.cursor() as cursor:
            try:
                cursor.execute('delete from well_query_index;')
            except Exception as e:
                logger.exception('on delete well_query_index')
            try:
                cursor.execute('delete from well_data_column_positive_index;')
            except Exception as e:
                logger.exception('on delete well_data_column_positive_index')
            try:
                cursor.execute('delete from cached_query;')
            except Exception as e:
                logger.exception('on delete cached_query')
        Screen.objects.all().delete()
        Library.objects.all().delete()
        LabAffiliation.objects.all().delete()
        ApiLog.objects.all().delete()
        ScreensaverUser.objects.all().exclude(username='testsuper').delete()

        # ScreensaverUser.objects.all().filter(username='adminuser').delete()

    def _setup_test_config(self):

        # Setup ScreenResult dependencies
        
        logger.info('create library...')
        library1 = self.create_library({
            'start_plate': 1, 
            'end_plate': 20,
            'screen_type': 'small_molecule' })

        logger.info('create and load well data...')
        plate = 1
        exp_well_count = 20
        input_data = [
            self.create_small_molecule_test_well(
                plate,i,library_well_type='experimental') 
            for i in range(0,exp_well_count-1)]
        
        logger.info('library plate_size: %r', library1)
        control_well_index = exp_well_count
        self.control_well1 = \
            lims_utils.well_name_from_index(
                control_well_index, int(library1['plate_size']))
        logger.info('control well 1: %r', self.control_well1)
        input_data.append(
            self.create_small_molecule_test_well(
                plate,control_well_index,library_well_type='empty') 
            )
        logger.debug('input_data: %r', input_data)
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library1['short_name'],resource_name])
        resp = self.api_client.put(
            resource_uri, format='sdf', data={ 'objects': input_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

    def test1_load_example_file(self):
        
        logger.info('test1_load_example_file...')
        
        default_data_for_get = { 'limit': 0, 'includes': ['*'] }
        default_data_for_get['HTTP_AUTHORIZATION'] = self.get_credentials()
        
        logger.info('create library...')
        library1 = self.create_library({
            'start_plate': 1, 
            'end_plate': 20,
            'screen_type': 'small_molecule' })

        logger.info('create screen...')        
        screen = self.create_screen({ 'screen_type': 'small_molecule' })

        logger.info('create and load well data...')
        input_data = []
        for plate in range(1,4):
            for i in range(0,384):
                
                if plate == 1 and i in [20,21,22]: # A21, 22, 23
                    # setup control wells
                    input_data.append(self.create_small_molecule_test_well(
                        plate,i,library_well_type='empty'))
                else:
                    input_data.append(self.create_small_molecule_test_well(
                        plate,i,library_well_type='experimental'))
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library1['short_name'],resource_name])
        data_for_get = {}
        data_for_get.update(default_data_for_get)
        data_for_get['CONTENT_TYPE'] = JSON_MIMETYPE
        logger.info('PUT library well data ...')
        resp = self.django_client.put(
            resource_uri, data=json.dumps({ 'objects': input_data }) ,
            **data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        data_for_get = {}
        data_for_get.update(default_data_for_get)
        data_for_get['CONTENT_TYPE'] = XLSX_MIMETYPE
        data_for_get['HTTP_ACCEPT'] = XLSX_MIMETYPE
        
        logger.info('load the screen result file...')
        file = 'ScreenResultTest_1_valid.xlsx'
        filename = '%s/db/static/test_data/screens/%s' %(APP_ROOT_DIR,file)
        logger.info('Open and POST file: %r', filename)
        with open(filename) as input_file:
            resource_uri = '/'.join([
                BASE_URI_DB,'screenresult',screen['facility_id']])
            logger.info('POST screen result to the server...')
            # NOTE: content_type arg is req'd with django.test.Client.post
            resp = self.django_client.post(
                resource_uri, content_type=XLSX_MIMETYPE,
                data=input_file.read(), **data_for_get)
            if resp.status_code not in [200, 204]:
                content = self.get_content(resp)
                if content:
                    logger.info('resp: %r', 
                        [[str(y) for y in x] 
                            for x in self.serializer.from_xlsx(content)])
            self.assertTrue(
                resp.status_code in [200, 204], 
                (resp.status_code))
            input_file.seek(0)
            input_data = self.sr_serializer.from_xlsx(input_file.read())            

        logger.info('screen result loaded, refetch from server...')

        # TODO: replace all tastypie.test.TestApiClient clients with the
        # django.test.client.Client instances:
        # Tastypie mucks with the HEADERS and uses non-standard "format" arg:
        # resp = self.sr_api_client.get(
        #     resource_uri, authentication=self.get_credentials(), 
        #     format='xlsx', **data_for_get)
        
        resp = self.django_client.get(
            resource_uri, authentication=self.get_credentials(), 
            **data_for_get)
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        output_data = self.sr_serializer.deserialize(
            self.get_content(resp), XLSX_MIMETYPE)
        
        ScreenResultSerializerTest.validate(self, input_data, output_data)    
    
    def _create_valid_input(self, screen_facility_id):
        # create data as already parsed input
        # TODO: rework the screenresult import into a generic serialization:
        # fields: keyed by (unique) datacolumn name
        # result_values: as an array of dicts (using a generator)
        # result_values: dropt the "numeric_value/value"        
        fields = {
            'E': {
                'ordinal': 0,
                'name': 'Field1',
                'data_worksheet_column': 'E',
                'data_type': 'text', 
                'description': 'field 1 description',
                'replicate_ordinal': 1,
            },
            'F': {
                'ordinal': 1,
                'name': 'Field2',
                'data_worksheet_column': 'F',
                'data_type': 'numeric',
                'decimal_places': 2, 
                'description': 'field 2 description',
                'replicate_ordinal': 1,
                'is_follow_up_data': True,
                'assay_readout_type': 'luminescence',
            },
            'G': {
                'ordinal': 2,
                'name': 'Field3',
                'data_worksheet_column': 'G',
                'data_type': 'numeric',
                'decimal_places': 4, 
                'description': 'field 3 description',
                'replicate_ordinal': 2,
                'is_follow_up_data': True,
                'assay_readout_type': 'flourescence_intensity',
            },
            'H': {
                'ordinal': 3,
                'name': 'Field4',
                'data_worksheet_column': 'H',
                'data_type': 'confirmed_positive_indicator',
                'description': 'field 4 description',
                'how_derived': 'z-score > 1 for Field3'
            },
            'I': {
                'ordinal': 4,
                'name': 'Field5',
                'data_worksheet_column': 'I',
                'data_type': 'partition_positive_indicator',
                'description': 'field 5 description',
                'how_derived': 'ranking of z-score <.5, ,.5<=x<=1, >1'
            },
        }
        result_values = [
            { 
                'well_id': '00001:A01', 
                'E': 'test value 1',
                'F': 91.19 ,
                'G': .0011 ,
                'H': None ,  # should (*NOT 20170724) be interpreted as 'NT' 
                'I': None , # should (*NOT 20170724) be interpreted as 'NP'
            },
            { 
                'well_id': '00001:A02', 
                'E': 'test value 2',
                'F': 0.99 ,
                'G': 1.0331 ,
                'H': 'CP',
                'I': 'W' ,
            },
            { 
                'well_id': '00001:A03', 
                'E': 'test value 2',
                'F': 1.99 ,
                'G': 1.032 ,
                'H': 'I',
                'I': 'M' ,
            },
            { 
                'well_id': '00001:A04', 
                'E': 'test value 2',
                'F': 1.99 ,
                'G': 1.032 ,
                'H': 'NT',
                'I': 'S' ,
            },
            { 
                'well_id': '00001:A05', 
                'E': 'test value 2',
                'F': 1.1 ,
                'G': 1.1 ,
                'H': 'CP',
                'I': 'M' ,
            },
            { 
                'well_id': '00001:A06', 
                'E': 'test value 2',
                'F': 1.1 ,
                'G': 1.1 ,
                'H': 'FP',
                'I': 'M' ,
            },
            { 
                'well_id': '00001:A07',
                'exclude': ['E','F','G','H','I'], 
                'E': 'test value 2',
                'F': 1.1 ,
                'G': 1.1 ,
                'H': 'FP',
                'I': 'M' ,
            },
            { 
                'well_id': '00001:%s' % self.control_well1,
                'assay_well_control_type': 'assay_control', 
                'E': 'test value 2',
                'F': 1.1 ,
                'G': 1.1 ,
                # NOTE: 20170724 - changed importer to raise a ValidationError here
                # 'H': 'FP', # non experimental well should be ignored
                # 'I': 'M' , # non experimental well should be ignored
            },
        ]
        input_data = OrderedDict((
            ('meta', {'screen_facility_id': screen_facility_id } ),
            ('fields', fields),
            ('objects', result_values),
        ))
        return input_data
    
    def test2_load_valid_input(self):
        
        logger.info('test2_load_valid_input...')
        default_data_for_get = { 'limit': 0, 'includes': ['*'] }
        default_data_for_get['HTTP_AUTHORIZATION'] = self.get_credentials()
        
        self._setup_test_config()
        
        logger.info('create screen...')        
        screen = self.create_screen({ 'screen_type': 'small_molecule' })
        self.screen1 = screen
        
        screen_facility_id = screen['facility_id']
        logger.info('created screen %r', screen_facility_id)        
        
        logger.info('create test screen result...')
        input_data = self._create_valid_input(screen_facility_id)

        # The ScreenResultSerializer only recognizes the XLSX format:
        # So serialize the input data into an XLSX file
        input_data_put = screen_result_importer.create_output_data(
            screen_facility_id, 
            input_data['fields'], 
            input_data['objects'] )
        input_data_put = self.sr_serializer.serialize(
            input_data_put, XLSX_MIMETYPE)
        data_for_get = {}
        data_for_get.update(default_data_for_get)
        data_for_get['CONTENT_TYPE'] = XLSX_MIMETYPE
        data_for_get['HTTP_ACCEPT'] = XLSX_MIMETYPE
        screen_facility_id = screen['facility_id']
        resource_name = 'screenresult'
        resource_uri = '/'.join([
            BASE_URI_DB,resource_name,screen_facility_id])
        
        logger.info('PUT screen result to the server...')
        resp = self.django_client.put(
            resource_uri, data=input_data_put, **data_for_get )
        if resp.status_code not in [200, 204]:
            content = self.get_content(resp)
            if content:
                logger.info('resp: %r', 
                    [[str(y) for y in x] 
                        for x in self.serializer.from_xlsx(content)])
        self.assertTrue(
            resp.status_code in [200, 204], resp.status_code)

        logger.info('refetch screen result from server...')
        resp = self.django_client.get(resource_uri, **data_for_get)
        if resp.status_code not in [200, 204]:
            content = self.get_content(resp)
            if content:
                logger.info('resp: %r', 
                    [[str(y) for y in x] 
                        for x in self.serializer.from_xlsx(content)])
        self.assertTrue(resp.status_code in [200,201,202],resp.status_code)
        output_data = self.sr_serializer.deserialize(
            self.get_content(resp), XLSX_MIMETYPE)
        ScreenResultSerializerTest.validate(self, input_data, output_data)
        
        # Test statistics:
        # result value count
        # experimental wells, replicates, etc
        
        screen = self.get_screen(screen['facility_id'])
        
        key = 'experimental_well_count'
        expected_value = 7
        self.assertTrue(screen[key]==expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))
        key = 'library_plates_data_loaded'
        expected_value = 1
        self.assertTrue(screen[key]==expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))
        
        # NOTE: data_loaded statistics have been removed
        # key = 'min_data_loaded_replicate_count'
        # expected_value = 2
        # self.assertTrue(screen[key]==expected_value,
        #     (key,'expected_value',expected_value,
        #         'returned value',screen[key]))
        # key = 'max_data_loaded_replicate_count'        
        # expected_value = 2
        # self.assertTrue(screen[key]==expected_value,
        #     (key,'expected_value',expected_value,
        #         'returned value',screen[key]))

        key = 'assay_readout_types'
        expected_value = ['luminescence', 'flourescence_intensity']
        self.assertTrue(set(screen[key]) <= set(expected_value),
            (key,'expected_value',expected_value,
                'returned value',screen[key]))

        # Test Datacolumn positive counts
        
        resource_name = 'datacolumn'
        resource_uri = '/'.join([BASE_URI_DB,resource_name])
        data_for_get = {
            'screen_facility_id__eq': screen_facility_id,
            'includes': '*',
            'assay_data_type__in': [
              'partition_positive_indicator','boolean_positive_indicator',
              'confirmed_positive_indicator'],
            'order_by': ['ordinal']
        }
        data_for_get.setdefault('limit', 0)
        resp = self.api_client.get(
            resource_uri, authentication=self.get_credentials(), 
            format='json', data=data_for_get)
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        output_data = self.deserialize(resp)
        logger.debug('datacolumn output_data: %r', output_data )
        self.assertTrue(len(output_data['objects'])==2, 
            ('should show two positive indicator columns', output_data))

        confirmed_positive_col = None
        partion_positive_col = None
        self.screen1_positive_columns = []
        for col in output_data['objects']:
            if col['assay_data_type'] == 'confirmed_positive_indicator':
                confirmed_positive_col = col
                self.screen1_positive_columns.append(col['key'])
            if col['assay_data_type'] == 'partition_positive_indicator':
                partion_positive_col = col
                self.screen1_positive_columns.append(col['key'])
        self.assertTrue(confirmed_positive_col is not None, 
            ('confirmed_positive_col not found in %r',output_data))
        self.assertTrue(partion_positive_col is not None, 
            ('partion_positive_col not found in %r',output_data))
        
        key = 'positives_count'
        expected_value=2
        self.assertTrue(confirmed_positive_col[key]==expected_value,
            (key,'expected_value',expected_value,
                'returned value',confirmed_positive_col[key], 
                'col', confirmed_positive_col))
        key = 'positives_count'
        expected_value=5
        self.assertTrue(partion_positive_col[key]==expected_value,
            (key,'expected_value',expected_value,
                'returned value',partion_positive_col[key]))
        key = 'strong_positives_count'
        expected_value=1
        self.assertTrue(partion_positive_col[key]==expected_value,
            (key,'expected_value',expected_value,
                'returned value',partion_positive_col[key]))
        key = 'medium_positives_count'
        expected_value=3
        self.assertTrue(partion_positive_col[key]==expected_value,
            (key,'expected_value',expected_value,
                'returned value',partion_positive_col[key]))
        key = 'weak_positives_count'
        expected_value=1
        self.assertTrue(partion_positive_col[key]==expected_value,
            (key,'expected_value',expected_value,
                'returned value',partion_positive_col[key]))

        logger.info('test2_load_valid_input, done.')        

    def test3_mutual_positives(self):

        logger.info('test3_mutual_positives...')
        default_data_for_get = { 'limit': 0, 'includes': ['*'] }
        default_data_for_get['HTTP_AUTHORIZATION'] = self.get_credentials()
        
        self.test2_load_valid_input()
         
        logger.info(
            'Mutual positives: Create another an overlapping screen result...')
        logger.info('create screen...')        
        screen = self.create_screen({ 'screen_type': 'small_molecule' })
        self.screen2 = screen
        screen_facility_id = screen['facility_id']
        logger.info('created screen %r', screen_facility_id)        
        
        logger.info('create second screen result...')
        input_data = self._create_valid_input(screen_facility_id)

        # The ScreenResultSerializer only recognizes the XLSX format:
        # So serialize the input data into an XLSX file
        input_data_put = screen_result_importer.create_output_data(
            screen_facility_id, 
            input_data['fields'], 
            input_data['objects'] )
        input_data_put = self.sr_serializer.serialize(
            input_data_put, XLSX_MIMETYPE)
        data_for_get = {}
        data_for_get.update(default_data_for_get)
        data_for_get['CONTENT_TYPE'] = XLSX_MIMETYPE
        data_for_get['HTTP_ACCEPT'] = XLSX_MIMETYPE
        logger.info('PUT screen result to the server...')
        screen_facility_id = screen['facility_id']
        resource_name = 'screenresult'
        resource_uri = '/'.join([
            BASE_URI_DB,resource_name,screen_facility_id])
        logger.info('PUT the second screen result to the server...')
        resp = self.django_client.put(
            resource_uri, data=input_data_put, **data_for_get )
        if resp.status_code not in [200, 204]:
            content = self.get_content(resp)
            if content:
                logger.info('resp: %r', 
                    [[str(y) for y in x] 
                        for x in self.serializer.from_xlsx(content)])
        self.assertTrue(
            resp.status_code in [200, 204], resp.status_code)

        logger.info('refetch screen result from server...')
        resp = self.django_client.get(resource_uri, **data_for_get)
        if resp.status_code not in [200, 204]:
            content = self.get_content(resp)
            if content:
                logger.info('resp: %r', 
                    [[str(y) for y in x] 
                        for x in self.serializer.from_xlsx(content)])
        self.assertTrue(resp.status_code in [200,201,202],resp.status_code)
        output_data = self.sr_serializer.deserialize(
            self.get_content(resp), XLSX_MIMETYPE)
        ScreenResultSerializerTest.validate(self, input_data, output_data)
        
        # Verify the "screen.overlapping_positive_screens" field
        self.screen1 = self.get_screen(
            self.screen1['facility_id'], { 'includes': '*' })
        self.screen2 = self.get_screen(
            self.screen2['facility_id'], { 'includes': '*' })
        overlapping1 = self.screen1.get('overlapping_positive_screens', None)
        overlapping2 = self.screen2.get('overlapping_positive_screens', None)
        self.assertTrue(
            self.screen2['facility_id'] in overlapping1,
            '%r not found in %r' % (self.screen2['facility_id'], overlapping1) )
        self.assertTrue(
            self.screen1['facility_id'] in overlapping2,
            '%r not found in %r' % (self.screen1['facility_id'], overlapping2) )
        
        # Verify that the screen1 mutual columns = screen 2 mutual columns
        screen_facility_id = screen['facility_id']
        resource_name = 'screenresult'
        resource_uri = '/'.join([
            BASE_URI_DB,resource_name,screen_facility_id])
         # use JSON to retrieve the data, mutual positives columns would be 
         # ignored when deserializing from XLSX
        data_for_get = {}
        data_for_get.update(default_data_for_get)
        data_for_get['HTTP_ACCEPT'] = JSON_MIMETYPE
        data = {}
        data['includes'] = '*'
        data['show_mutual_positives'] = True
        logger.info('get the screen result: %r', data_for_get)
        resp = self.django_client.get(resource_uri, data=data, **data_for_get)
        if resp.status_code not in [200, 204]:
            content = self.get_content(resp)
            if content:
                logger.info('resp: %r', 
                    [[str(y) for y in x] 
                        for x in self.serializer.from_xlsx(content)])
        self.assertTrue(resp.status_code in [200,201,202],resp.status_code)
        logger.info('deserialize to json and check mutual positives columns...')
        # use generic deserialization to show all of the fields
        content = self.get_content(resp)
        output_data = self.serializer.deserialize(
            content, JSON_MIMETYPE)
        
        try:
            # save output data generator to a list
            output_list = [x for x in output_data['objects']]
        except Exception, e:
            logger.exception('%r', e)
            raise
        self.assertTrue(len(output_list)>0, 
            'no screen results in output: %r' % output_data)
        logger.info((
            'verify that mutual positive columns are in the result value '
            'fields: %r'), output_list[0] )

        mutual_positive_types = ['confirmed_positive_indicator','partition_positive_indicator']
        resource_name = 'screenresult/' + self.screen1['facility_id']
        screen1_schema = self.get_schema(resource_name)
        resource_name = 'screenresult/' + self.screen2['facility_id']
        screen2_schema = self.get_schema(resource_name)
        mutual_positive_cols_screen_1 = [field['key'] 
            for field in screen1_schema['fields'].values()
                if field.get('assay_data_type',None) in mutual_positive_types ]
        mutual_positive_cols_screen_2 = [field['key'] 
            for field in screen2_schema['fields'].values()
                if field.get('assay_data_type',None) in mutual_positive_types ]
        for col in mutual_positive_cols_screen_1:
            self.assertTrue(col in output_list[0].keys(), 
                ('mutual positive column %r is missing in output cols: %r' 
                 % (col, output_list[0].keys())))
        # if a col is positive, check that the mutual pos cols have data as well
        
        for row in output_list:
            if row['is_positive']:
                for col1 in mutual_positive_cols_screen_1:
                    for col2 in mutual_positive_cols_screen_2:
                        self.assertTrue(col2 in row,
                            'could not find %r in row: %r'
                            % (col2, row))
                        if row[col2] is not None:
                            self.assertTrue(row[col1] is not None,
                                ('mutual positive column %r is missing data in row: %r' 
                                    % (col, row)))
            
        # NOTE: Further testing of mutual column field visibilty in the 
        # DataSharingLevel tests
        
    def test4_result_value_errors_from_file(self):
        
        logger.info('test4_result_value_errors_from_file...')
        logger.info('create library...')
        library1 = self.create_library({
            'start_plate': 1, 
            'end_plate': 20,
            'screen_type': 'small_molecule' })

        logger.info('create and load well data...')
        input_data = []
        plate = 1
        for i in range(0,384):
            if i in [0,16,112]: # A01, A04, A07
                # setup control wells
                input_data.append(self.create_small_molecule_test_well(
                    plate,i,library_well_type='empty'))
            else:
                input_data.append(self.create_small_molecule_test_well(
                    plate,i,library_well_type='experimental'))
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library1['short_name'],resource_name])
        logger.info('put library...')
        resp = self.api_client.put(
            resource_uri, format='sdf', data={ 'objects': input_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        logger.info('create screen...')        
        screen = self.create_screen({ 'screen_type': 'small_molecule' })
        
        data_for_get = { 'limit': 0, 'includes': ['*'] }
        data_for_get['CONTENT_TYPE'] = XLSX_MIMETYPE
        data_for_get['HTTP_ACCEPT'] = JSON_MIMETYPE
        data_for_get['HTTP_AUTHORIZATION'] = self.get_credentials()
        
        file = 'ScreenResultRVErrorsTest.xlsx'
        file_out = 'ScreenResultRVErrorsTest_out.xlsx'
        filename = '%s/db/static/test_data/screens/%s' %(APP_ROOT_DIR,file)
        filename_out = (
            '%s/db/static/test_data/screens/%s' %(APP_ROOT_DIR,file_out) )
        logger.info('Open and POST file: %r', filename)
        with open(filename, 'rb') as input_file:

            # NOTE: now posting as xls, so that the parse errors will be 
            # processed along with the validation errors
            
            resource_name = 'screenresult'
            screen_facility_id = screen['facility_id']
            resource_uri = '/'.join([
                BASE_URI_DB,resource_name,screen_facility_id])
            logger.info('PUT screen result to the server... %r', data_for_get)
            # NOTE: content_type arg is req'd with django.test.Client.post
            resp = self.django_client.post(
                resource_uri, content_type=XLSX_MIMETYPE,
                data=input_file.read(), **data_for_get)
            self.assertTrue(
                resp.status_code == 400, resp.status_code)
            logger.info('content-type: %r', resp['Content-Type'])
            content = self.deserialize(resp)
            logger.info('content; %r', content)
            
            # TODO: check for expected errors
            
            key = 'PL_00001'
            self.assertTrue(find_in_dict(key, content), 
                'key: %r should be in errors: %r' % (key, content))
            sheet_errors = find_in_dict(key, content)
            
            key = '00001:A01-G'
            self.assertTrue(key in sheet_errors, 
                'key: %r should be in errors: %r' % (key, sheet_errors))
            self.assertTrue(
                'could not convert' in str(sheet_errors[key]),
                'error should be a conversion error: %r, %r' 
                    % (key,sheet_errors[key]) )
            key = '00001:A03'
            self.assertTrue(key in sheet_errors, 
                'key: %r should be in errors: %r' % (key, sheet_errors))
            self.assertTrue(
                'unknown exclude' in str(sheet_errors[key]),
                'error should be a excluded cols error: %r, %r' 
                    % (key,sheet_errors[key]) )
            key = '00001:A05-E'
            self.assertTrue(key in sheet_errors, 
                'key: %r should be in errors: %r' % (key, sheet_errors))
            self.assertTrue(
                'could not convert' in str(sheet_errors[key]),
                'error should be a conversion error: %r, %r' 
                    % (key,sheet_errors[key]) )

    def test5_data_column_errors(self):
        
        logger.info('test5_data_column_errors...')
        logger.info('create library...')
        library1 = self.create_library({
            'start_plate': 1, 
            'end_plate': 20,
            'screen_type': 'small_molecule' })

        logger.info('create screen...')        
        screen = self.create_screen({ 'screen_type': 'small_molecule' })
        
        fields = {
            'E': {
                'ordinal': 0,
                'name': 'Field1',
                'data_worksheet_column': 'E',
                'data_type': 'text', 
                'description': 'field 1 description',
                'replicate_ordinal': 1,
            },
            'F': {
                'ordinal': 1,
                'name': 'Field2',
                'data_worksheet_column': 'F',
                #'data_type': 'numeric',
                'decimal_places': 2, 
                'description': 'field 2 description',
                'replicate_ordinal': 1,
                'is_follow_up_data': True,
                'assay_readout_type': 'luminescence',
            },
            'G': {
                'ordinal': 2,
                'name': 'Field3',
                'data_worksheet_column': 'G',
                'data_type': 'numeric',
                'decimal_places': 'x', 
                'description': 'field 3 description',
                'replicate_ordinal': 2,
                'is_follow_up_data': True,
                'assay_readout_type': 'flourescence_intensity',
            },
            'H': {
                'ordinal': 3,
                'name': 'Field4',
                'data_worksheet_column': 'H',
                'data_type': 'confirmed_positive_indicator',
                'description': 'field 4 description',
                'how_derived': 'z-score > 1 for Field3',
                'derived_from_columns': 'F,Y'
            },
            'I': {
                'ordinal': 4,
                'name': 'Field5',
                'data_worksheet_column': 'I',
                'data_type': 'partition_positive_indicator',
                'description': 'field 5 description',
                'how_derived': 'ranking of z-score <.5, ,.5<=x<=1, >1'
            },
        }
        result_values = [
            { 
                'well_id': '00001:A01', 
                'E': 'test value 1',
                'F': 91.19,
                'G': .0011,
                'H': None,  # should be interpreted as 'NT'
                'I': None, # should be interpreted as 'NP'
            },
        ]
        input_data = OrderedDict((
            ('meta', {'screen_facility_id': screen['facility_id'] } ),
            ('fields', fields),
            ('objects', result_values),
        ))
        
        # NOTE: now posting directly, w/out the tp client,
        # so that the parse errors will be 
        # processed along with the validation errors
        # Also, HTTP_ACCEPT is set to JSON, simulating the UI
        data_for_get = { 'limit': 0, 'includes': ['*'] }
        data_for_get['CONTENT_TYPE'] = JSON_MIMETYPE
        data_for_get['HTTP_ACCEPT'] = JSON_MIMETYPE
        data_for_get['HTTP_AUTHORIZATION'] = self.get_credentials()
        
        resource_name = 'screenresult'
        screen_facility_id = screen['facility_id']
        resource_uri = '/'.join([
            BASE_URI_DB,resource_name,screen_facility_id])
        logger.info('PUT screen result to the server...')
        # NOTE: content_type arg is req'd with django.test.Client.post
        resp = self.django_client.post(
            resource_uri, content_type='application/json', 
            data=json.dumps(input_data), **data_for_get)
        self.assertTrue(
            resp.status_code == 400, resp.status_code)
        content = self.deserialize(resp)
        logger.info('content; %r', content)
        
        key = 'F' # missing data_type
        self.assertTrue(find_in_dict(key, content), 
            'key: %r should be in errors: %r' % (key, content))
        col_errors = find_in_dict(key, content)
        self.assertTrue('data_type' in str(col_errors),
            'error should be a "data_type" error: %r, %r' %(key,col_errors) )
        
        key = 'G' # decimal_places not number
        self.assertTrue(find_in_dict(key, content), 
            'key: %r should be in errors: %r' % (key, content))
        col_errors = find_in_dict(key, content)
        self.assertTrue('decimal_places' in str(col_errors),
            'error should be a "decimal_places" error: %r, %r' %(key,col_errors) )
        
        key = 'H' # derived from col dne
        self.assertTrue(find_in_dict(key, content), 
            'key: %r should be in errors: %r' % (key, content))
        col_errors = find_in_dict(key, content)
        self.assertTrue(
            'derived_from_columns' in str(col_errors),
            'error should be a "derived_from_columns" error: %r, %r' 
                % (key,col_errors) )
        

class ScreenResource(DBResourceTestCase):
        
    def setUp(self):
        super(ScreenResource, self).setUp()

    def tearDown(self):
        DBResourceTestCase.tearDown(self)
        Screen.objects.all().delete()
        Library.objects.all().delete()
        LibraryScreening.objects.all().delete()
        ApiLog.objects.all().delete()
        ScreensaverUser.objects.all().exclude(username='testsuper').delete()
        LabAffiliation.objects.all().delete()
        
    def test1_create_screen(self):

        logger.info('test1_create_screen...')        
        data = {
            'screen_type': 'small_molecule',
            'cell_lines': ['293_hek_293','colo_858'],
        }
        screen_item = self.create_screen(data=data)
        
        self.assertTrue(
            'facility_id' in screen_item, 
            'the facility_id was not created')
        
        for key, value in data.items():
            self.assertEqual(value, screen_item[key], 
                'key %r, val: %r not expected: %r' 
                    % (key, value, screen_item[key]))
        logger.debug('screen created: %r', screen_item)

    def test1a_create_screen_with_facility_id_override(self):

        logger.info('test1_create_screen...')        
        data = {
            'screen_type': 'small_molecule',
            'cell_lines': ['293_hek_293','colo_858'],
            'species': 'bacteria',
            'facility_id': '10'
        }
        screen_item = self.create_screen(
            data=data, uri_params=['%s=true' % API_PARAM_OVERRIDE,])
        
        self.assertEqual(
            screen_item['facility_id'],data['facility_id'])
        
        for key, value in data.items():
            self.assertEqual(value, screen_item[key], 
                'key %r, val: %r not expected: %r' 
                    % (key, value, screen_item[key]))
        logger.debug('screen created with facility id: %r', screen_item)
        
    def test1b_create_follow_up_screen(self):    
        logger.info('test1_create_screen...')        
        data = {
            'screen_type': 'small_molecule',
            'cell_lines': ['293_hek_293','colo_858'],
            'species': 'bacteria',
        }
        screen_item = self.create_screen(data=data)
        
        self.assertTrue(
            'facility_id' in screen_item, 
            'the facility_id was not created')
        
        for key, value in data.items():
            self.assertEqual(value, screen_item[key], 
                'key %r, val: %r not expected: %r' 
                    % (key, value, screen_item[key]))
        logger.debug('screen created: %r', screen_item)
        
        follow_up_data = {
            'facility_id': screen_item['facility_id'] + 'A',
            'primary_screen': screen_item['facility_id'], 
            'title': screen_item['title']+'-follow up',
            'summary': screen_item['summary'],
            'lead_screener_id': screen_item['lead_screener_id'],
            'lab_head_id': screen_item['lab_head_id'],
            'collaborator_ids': screen_item['collaborator_ids'],
            'data_sharing_level': screen_item['data_sharing_level'],
            'screen_type': screen_item['screen_type']
        }
        
        resource_uri = BASE_URI_DB + '/screen?override=true'
        
        resp = self.api_client.post(
            resource_uri,format='json', 
            data=follow_up_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        data_for_get = {
            'facility_id': follow_up_data['facility_id']}
        follow_up_screen = self.get_single_resource(resource_uri, data_for_get)
        for key, value in follow_up_data.items():
            (result,msgs) = equivocal(value, follow_up_screen[key])
            self.assertTrue(result, 
                'key %r, val: %r not expected: %r' 
                    % (key, value, follow_up_screen[key]))
        
    def test1c_create_study(self):
        logger.info('test1c_create_study...')        

        lab_head = self.create_lab_head()
        lead_screener = self.create_screensaveruser()
        collaborator1 = self.create_screensaveruser()
        collaborator2 = self.create_screensaveruser()
        input_data = ScreenFactory.attributes()
        data = {
            'facility_id': '100000',
            'screen_type': 'small_molecule',
            'study_type': 'in_silico',
            'data_sharing_level': None,
            'lab_head_id': lab_head['screensaver_user_id'],
            'lead_screener_id': lead_screener['screensaver_user_id'],
            'collaborator_ids': [
                str(collaborator1['screensaver_user_id']), 
                str(collaborator2['screensaver_user_id']) ]
            
        }
        input_data.update(data)
        logger.info('input_data: %r', input_data)
        
        resource_uri = '/'.join([BASE_URI_DB, 'study'])
        resp = self.api_client.post(
            resource_uri, format='json', data=input_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code, self.get_content(resp)))
    
        new_obj = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in new_obj)
        self.assertEqual(len(new_obj[API_RESULT_DATA]), 1)        
        new_obj = new_obj[API_RESULT_DATA][0]
        logger.info('created %r', new_obj)

        data_for_get = { 
            'limit': 0,
            'includes': '*'
        }
        test_resource_uri = '/'.join([BASE_URI_DB, 'study', new_obj['facility_id']])
        screen_item =  self.get_single_resource(resource_uri, data_for_get)
        
        for key, value in data.items():
            if key == 'data_sharing_level':
                continue
            self.assertEqual(value, screen_item[key], 
                'key %r, val: %r not expected: %r' 
                    % (key, value, screen_item[key]))
            
        self.assertEqual(screen_item['data_sharing_level'], 0, 
            'if not specified, study data sharing level should be 0: %r' 
            % screen_item['data_sharing_level'])
        logger.debug('study created: %r', screen_item)
        
        # TODO: study result values
        
    def test2a_create_library_screening_cherry_picked_copies(self):
        logger.info('test2a_create_library_screening_cherry_picked_copies...')
                
        logger.info('A. Set up dependencies...')
        logger.info('create users...')
        self.screening_user = self.create_screensaveruser({ 
            'username': 'screening1'
        })
        
        # FIXME: create an "LibraryScreeningPerformers" admin group
        performed_by = self.admin_user
        
        logger.info('create libraries...')
        library1 = self.create_library({
            'start_plate': 1000, 
            'end_plate': 1005,
            'screen_type': 'small_molecule' })

        library2 = self.create_library({
            'start_plate': 2000, 
            'end_plate': 2040,
            'screen_type': 'small_molecule' })

        logger.info('set some experimental wells...')
        plate = 1000
        experimental_well_count = 20
        input_well_data = [
            self.create_small_molecule_test_well(
                plate,i,library_well_type='experimental',
                molar_concentration='0.001') 
            for i in range(0,experimental_well_count)]
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library1['short_name'],resource_name])
        resp = self.api_client.put(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        logger.info('create library copy...')
        library_copy1_input = {
            'library_short_name': library1['short_name'],
            'copy_name': "A",
            'usage_type': "library_screening_plates",
            'initial_plate_well_volume': '0.000010',
            'initial_plate_status': 'available'
        }  
        resource_uri = BASE_URI_DB + '/librarycopy'
        resource_test_uri = '/'.join([
            resource_uri,library_copy1_input['library_short_name'],
            library_copy1_input['copy_name']])
        library_copy1 = self._create_resource(
            library_copy1_input, resource_uri, resource_test_uri, 
            excludes=['initial_plate_well_volume','initial_plate_status'])
        logger.info('created: %r', library_copy1)
 
        library_copy2_input = library_copy1_input.copy()
        library_copy2_input['library_short_name'] = library2['short_name']
        resource_test_uri = '/'.join([
            resource_uri,library_copy2_input['library_short_name'],
            library_copy2_input['copy_name']])
        library_copy2 = self._create_resource(
            library_copy2_input, resource_uri, resource_test_uri,
            excludes=['initial_plate_well_volume','initial_plate_status'])
        logger.info('created: %r', library_copy2)
        
        logger.info('create screen...')        
        screen = self.create_screen({
            'screen_type': 'small_molecule'
            })

        logger.info('1. Modify copy-wells (to be used for screening)')

        plate_to_modify = 1000
        plate_size = 384
        logger.info('retrieve the copy_wells for plate: %d ...', plate_to_modify)
        resource_uri = '/'.join([
            BASE_URI_DB,'library',library1['short_name'],'copy',
            library_copy1_input['copy_name'],'copywell'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        data_for_get['plate_number__eq'] = plate_to_modify
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        copywell_data = new_obj[API_RESULT_DATA]
        expected_wells = plate_size;
        logger.debug('returned copywell: %r', copywell_data[0])
        logger.debug('(last) returned copywell: %r', copywell_data[-1])
        self.assertEqual(len(copywell_data),expected_wells)
        
        # Only modify one well here
        copywell_input = copywell_data[0]
        copywell_id = copywell_input['copywell_id']
        logger.info('copywell input: %r', copywell_input)
        copywell_plate = '%s/%s' % (
            copywell_input['copy_name'],str(copywell_input['plate_number']))
        logger.info('1.A Patch the copy well: %r', copywell_input)
        
        original_volume = Decimal(copywell_input['volume'])
        volume_adjustment = Decimal('0.000004')
        copywell_input['volume'] = original_volume-volume_adjustment
        
        patch_uri = '/'.join([resource_uri,copywell_input['well_id']]) 
        resp = self.api_client.patch(
            patch_uri, format='json', data=copywell_input, 
            authentication=self.get_credentials(), **data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        patch_response = self.deserialize(resp)
        
        self.assertTrue(API_RESULT_DATA in patch_response)
        self.assertEqual(len(patch_response[API_RESULT_DATA]), 1)        
        new_copywell = patch_response[API_RESULT_DATA][0]
        self.assertEqual(
            volume_adjustment, Decimal(new_copywell['consumed_volume']))
        self.assertEqual(
            Decimal(copywell_input['volume']), Decimal(new_copywell['volume']))

        logger.info('2. Create valid library screening input...')

        lps_format = \
            '{library_short_name}:{copy_name}:{{start_plate}}-{{end_plate}}'
        plate_range1 = lps_format.format(**library_copy1).format(**library1)
        plate_range2 = lps_format.format(**library_copy2).format(
                start_plate=library2['start_plate'],
                end_plate=int(library2['start_plate']+10))
        library_plates_screened = [ plate_range1, plate_range2 ]
        volume_to_transfer = "0.000000600"
        library_screening_input = {
            'screen_facility_id': screen['facility_id'],
            'date_of_activity': "2008-01-18",
            # REMOVED: 20170412 per JAS/KR
            # 'assay_protocol':'test assay protocol',
            # 'assay_protocol_type': '',
            'is_for_external_library_plates': False,
            'library_plates_screened': library_plates_screened,
            'number_of_replicates': 2,
            'performed_by_user_id': performed_by['screensaver_user_id'],
            'volume_transferred_per_well_from_library_plates': volume_to_transfer,
            'volume_transferred_per_well_to_assay_plates': "0.000000300"
        }
        resource_uri = BASE_URI_DB + '/libraryscreening'
        resource_test_uri = BASE_URI_DB + '/libraryscreening'

        resp = self.api_client.post(
            resource_uri,format='json', 
            data=library_screening_input, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 1.a Inspect meta "Result" section
        resp = self.deserialize(resp)
        logger.info('resp: %r', resp)
        
        self.assertTrue(API_RESULT_META in resp, '%r' % resp)
        self.assertTrue(API_MSG_RESULT in resp[API_RESULT_META], '%r' % resp)
        
        actual_result_meta = resp[API_RESULT_META][API_MSG_RESULT]
        expected_result_meta = {
            API_MSG_SCREENING_PLATES_UPDATED: 17,
            API_MSG_SCREENING_ADDED_PLATE_COUNT: 17, 
            API_MSG_SCREENING_DELETED_PLATE_COUNT: 0,
            API_MSG_SCREENING_EXTANT_PLATE_COUNT: 0,
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 17,
            API_MSG_COPYWELLS_ALLOCATED: { 
                copywell_plate: 1 }
        }
        for k,v in expected_result_meta.items():
            self.assertEqual(v,actual_result_meta[k], 
                'k: %r, expected: %r, actual: %r' % (k, v, actual_result_meta))
        
        self.assertTrue(len(resp[API_RESULT_DATA]) == 1)
        library_screening_output = resp[API_RESULT_DATA][0]
        
        # 1.b Inspect the returned library screening
        for k,v in library_screening_input.items():
            v2 = library_screening_output[k]
            self.assertTrue(equivocal(v,v2),
                'test key: %r:%r != %r' % (k, v, v2))
        
        # 1.c check the copywell separately
        cw_resource_uri = '/'.join([
            BASE_URI_DB,'library',library1['short_name'],'copy',
            library_copy1_input['copy_name'],'copywell'])
        cw_data_for_get ={ 'limit': 0 }        
        cw_data_for_get.setdefault('includes', ['*'])
        cw_data_for_get['well_id__eq'] = copywell_input['well_id']
        cw_data_for_get['copy_name__eq'] = copywell_input['copy_name']
        resp = self.api_client.get(
            cw_resource_uri, format='json', 
            authentication=self.get_credentials(), 
            data=cw_data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        self.assertTrue(API_RESULT_DATA in new_obj)
        self.assertEqual(len(new_obj[API_RESULT_DATA]), 1)        
        new_copywell_data = new_obj[API_RESULT_DATA][0]
        expected_volume = (
            Decimal(copywell_input['volume']) - Decimal(volume_to_transfer))
        self.assertEqual(
            Decimal(new_copywell_data['volume']),expected_volume,
            'copywell_data: vol expected %r, found: %r' % (
                expected_volume, Decimal(new_copywell_data['volume'])) )

        # 1.D Check logs (17 plate, 1 copywell)
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'ref_resource_name': 'libraryscreening', 
            'key': library_screening_output['activity_id'],
            'api_action': API_ACTION_CREATE
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        logger.debug('logs: %d', len(apilogs))
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('apilog: %r', apilog)
        
        expected_child_logs = 18 #  17 plate log, 1 copywell log
        self.assertTrue(apilog['child_logs'], expected_child_logs)
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'parent_log_id': apilog['id']
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        self.assertEqual(
            len(apilogs),expected_child_logs, 
            'wrong number of child logs returned: %d' % len(apilogs))
        
        for apilog in apilogs:
            logger.info('apilog: %r, %r, %r', 
                apilog['uri'],apilog['ref_resource_name'],apilog['key'])
            if apilog['ref_resource_name'] == 'librarycopyplate':
                logger.info('lcp log: %r', apilog)
                self.assertTrue('remaining_well_volume' in apilog['diff_keys'])
                self.assertTrue('screening_count' in apilog['diff_keys'])
            elif apilog['ref_resource_name'] == 'copywell':
                logger.info('copywell log: %r', apilog)
                self.assertEqual(copywell_id, apilog['key'])
                self.assertTrue('volume' in apilog['diff_keys'])
            else:
                self.fail('unknown log: %r', apilog)

        # 2. deallocate (remove) the screening plate with the copy and check  
        # copywell volume has been adjusted
        plate_range1 = lps_format.format(**library_copy1).format(
            start_plate=library1['start_plate']+1, end_plate=library1['end_plate'])
        plate_removed_key = '/'.join([ library1['short_name'],
            library_copy1['copy_name'], str(library1['start_plate']) ])
        plate_range2 = lps_format.format(**library_copy2).format(
                start_plate=library2['start_plate'],
                end_plate=int(library2['start_plate']+10))
        library_plates_screened = [ plate_range1, plate_range2 ]
        new_library_screening_input = library_screening_output.copy()
        new_library_screening_input['library_plates_screened'] = \
            library_plates_screened
        resource_uri = BASE_URI_DB + '/libraryscreening/' \
            + str(new_library_screening_input['activity_id'])

        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=new_library_screening_input, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 2.a Inspect meta "Result" section
        resp = self.deserialize(resp)
        logger.info('resp: %r', resp)
        
        self.assertTrue(API_RESULT_META in resp, '%r' % resp)
        self.assertTrue(API_MSG_RESULT in resp[API_RESULT_META], '%r' % resp)
        
        actual_result_meta = resp[API_RESULT_META][API_MSG_RESULT]
        expected_result_meta = {
            API_MSG_SCREENING_PLATES_UPDATED: 1,
            API_MSG_SCREENING_ADDED_PLATE_COUNT: 0, 
            API_MSG_SCREENING_DELETED_PLATE_COUNT: 1,
            API_MSG_SCREENING_EXTANT_PLATE_COUNT: 16,
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 16,
            API_MSG_COPYWELLS_DEALLOCATED: { 
                copywell_plate: 1 }
        }
        for k,v in expected_result_meta.items():
            self.assertEqual(v,actual_result_meta[k], 
                'k: %r, expected: %r, actual: %r' % (k, v, actual_result_meta))
        
        self.assertTrue(len(resp[API_RESULT_DATA]) == 1)
        library_screening_output = resp[API_RESULT_DATA][0]
        
        # 2.c check the copywell separately
        resp = self.api_client.get(
            cw_resource_uri, format='json', 
            authentication=self.get_credentials(), 
            data=cw_data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        self.assertTrue(API_RESULT_DATA in new_obj)
        self.assertEqual(len(new_obj[API_RESULT_DATA]), 1)        
        deallocated_copywell_data = new_obj[API_RESULT_DATA][0]
        expected_volume = (
            Decimal(new_copywell_data['volume']) + Decimal(volume_to_transfer))
        self.assertEqual(
            Decimal(deallocated_copywell_data['volume']),expected_volume,
            'copywell_data: vol expected %r, found: %r' % (
                expected_volume, Decimal(deallocated_copywell_data['volume'])))
        
        # 2.D Check logs (1 plate, 1 copywell)
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'ref_resource_name': 'libraryscreening', 
            'key': library_screening_output['activity_id'],
            'api_action': API_ACTION_PATCH
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        logger.debug('logs: %d', len(apilogs))
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('apilog: %r', apilog)
        
        self.assertTrue('library_plates_screened_count' in apilog['diffs'])
        self.assertTrue( '["17", "16"]' in apilog['diffs'])

        expected_child_logs = 2 #  1 plate log, 1 copywell log
        self.assertTrue(apilog['child_logs'], expected_child_logs)
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'parent_log_id': apilog['id']
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        self.assertEqual(
            len(apilogs),expected_child_logs, 
            'wrong number of child logs returned: %d: %r' 
            % (len(apilogs), apilogs))
        
        for apilog in apilogs:
            logger.info('apilog: %r, %r, %r', 
                apilog['uri'],apilog['ref_resource_name'],apilog['key'])
            if apilog['ref_resource_name'] == 'librarycopyplate':
                logger.info('lcp log: %r', apilog)
                self.assertEqual(plate_removed_key, apilog['key'])
            elif apilog['ref_resource_name'] == 'copywell':
                self.assertEqual(copywell_id, apilog['key'])
                logger.info('copywell log: %r', apilog)
            else:
                self.fail('unknown log: %r', apilog)
                
                
    def test2_create_library_screening(self):

        logger.info('test2_create_library_screening...')
                
        # Set up dependencies
        logger.info('create users...')
        self.screening_user = self.create_screensaveruser({ 
            'username': 'screening1'
        })
        
        # FIXME: create an "LibraryScreeningPerformers" admin group
        performed_by_user = self.admin_user
        
        logger.info('create library...')
        library1 = self.create_library({
            'start_plate': 1000, 
            'end_plate': 1005,
            'screen_type': 'small_molecule' })

        library2 = self.create_library({
            'start_plate': 2000, 
            'end_plate': 2040,
            'screen_type': 'small_molecule' })

        logger.info('set some experimental wells...')
        plate = 1000
        experimental_well_count = 20
        input_well_data = [
            self.create_small_molecule_test_well(
                plate,i,library_well_type='experimental',
                molar_concentration='0.001') 
            for i in range(0,experimental_well_count)]
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library1['short_name'],resource_name])
        resp = self.api_client.put(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        logger.info('create library copy...')
        
        # Start with less volume than required
        min_allowed_small_molecule_volume = \
            LibraryScreeningResource.MIN_WELL_VOL_SMALL_MOLECULE
        valid_test_volume = Decimal('0.000000200') # 200nL
        
        library_copy1_input = {
            'library_short_name': library1['short_name'],
            'copy_name': "A",
            'usage_type': "library_screening_plates",
            'initial_plate_well_volume':
                 '{:.9f}'.format(min_allowed_small_molecule_volume + valid_test_volume*5),
            'initial_plate_status': 'available'
        }  
        resource_uri = BASE_URI_DB + '/librarycopy'
        resource_test_uri = '/'.join([
            resource_uri,library_copy1_input['library_short_name'],
            library_copy1_input['copy_name']])
        library_copy1 = self._create_resource(
            library_copy1_input, resource_uri, resource_test_uri, 
            excludes=['initial_plate_well_volume','initial_plate_status'])
        logger.info('created: %r', library_copy1)

        library_copy1_cpp_input = library_copy1_input.copy()
        library_copy1_cpp_input['usage_type'] = 'cherry_pick_source_plates'
        library_copy1_cpp_input['copy_name'] = 'A1_cpsp'
        resource_test_uri = '/'.join([
            resource_uri,library_copy1_cpp_input['library_short_name'],
            library_copy1_cpp_input['copy_name']])
        library_copy1_cpp = self._create_resource(
            library_copy1_cpp_input, resource_uri, resource_test_uri, 
            excludes=['initial_plate_well_volume','initial_plate_status'])
        logger.info('created: %r', library_copy1_cpp)
        
        library_copy2_input = library_copy1_input.copy()
        library_copy2_input['library_short_name'] = library2['short_name']
        resource_test_uri = '/'.join([
            resource_uri,library_copy2_input['library_short_name'],
            library_copy2_input['copy_name']])
        library_copy2 = self._create_resource(
            library_copy2_input, resource_uri, resource_test_uri,
            excludes=['initial_plate_well_volume','initial_plate_status'])
        logger.info('created: %r', library_copy2)
        
        logger.info('create screen...')        
        screen = self.create_screen({
            'screen_type': 'small_molecule'
            })

        lps_format = \
            '{library_short_name}:{copy_name}:{{start_plate}}-{{end_plate}}'
        plate_range1 = lps_format.format(**library_copy1).format(**library1)
        plate_range2 = lps_format.format(**library_copy2).format(
                start_plate=library2['start_plate'],
                end_plate=int(library2['start_plate']+10))
        library_plates_screened = [ plate_range1, plate_range2 ]
        
        library_screening_input = {
            'screen_facility_id': screen['facility_id'],
            'date_of_activity': "2008-01-18",
            # REMOVED: 20170412 per JAS/KR
            # 'assay_protocol':'test assay protocol',
            # 'assay_protocol_type': '',
            'is_for_external_library_plates': False,
            'library_plates_screened': library_plates_screened,
            'number_of_replicates': 2,
            'performed_by_user_id': performed_by_user['screensaver_user_id'],
            'volume_transferred_per_well_from_library_plates': 
                '{0:.9f}'.format(valid_test_volume*2),
            'volume_transferred_per_well_to_assay_plates': 
                '{0:.9f}'.format(valid_test_volume)
        }
        resource_uri = BASE_URI_DB + '/libraryscreening'
        resource_test_uri = BASE_URI_DB + '/libraryscreening'
        data_for_get = {
            'screen_facility_id__eq': screen['facility_id']
        }
        
        # TESTS - Begin with "failing" tests
        
        logger.info('1. Creating (failed) library_screening with invalid inputs...')
        key = 'screen_facility_id' 
        msg = 'input missing a %r should fail' % key
        logger.info('test %r', msg)
        invalid_input1 = library_screening_input.copy()
        del invalid_input1['screen_facility_id'] 
        errors, resp = self._create_resource(
            invalid_input1, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(
            find_in_dict(key, errors), 
            'test: %s, not in errors: %r' %(msg,errors))
        
        key = 'library_plates_screened' 
        logger.info('2. create a library_plate_range with invalid library name')
        # Even though the plate range is correct, this should fail,
        # and the transaction should cancel the creation of the libraryscreening
        value = [ lps_format.format(
            library_short_name='**',copy_name='A').format(**library1)]
        msg = 'invalid format for %r should fail' % key
        logger.info('test %r', msg)
        invalid_input2 = library_screening_input.copy()
        invalid_input2[key] =  value
        errors, resp = self._create_resource(
            invalid_input2, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(find_in_dict(key, errors), 
            'Error: response error not found: %r, obj: %r' %(key, errors))
        
        logger.info('3. test invalid plate range')
        key = 'library_plates_screened' 
        value = [ lps_format.format(**library_copy1).format(**{
            'start_plate': library1['start_plate'],
            'end_plate': int(library1['start_plate'])-1 }) ]
        msg = 'invalid plate range in %r for  %r should fail' % (value,key)
        invalid_input3 = library_screening_input.copy()
        invalid_input3[key] =  value
        errors, resp = self._create_resource(
            invalid_input3, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(find_in_dict(key, errors), 
            'test: %s, not in errors: %r' %(key,errors))

        logger.info('4. test overlapping plate range... ')
        key = 'library_plates_screened' 
        value = [ 
            lps_format.format(**library_copy2).format(**{
                'start_plate': library2['start_plate'],
                'end_plate': int(library2['start_plate'])+2 }),
            lps_format.format(**library_copy2).format(**{
                'start_plate': library2['start_plate']+1,
                'end_plate': int(library2['start_plate'])+4 }),
        ]
        msg = '4 - overlapping plate ranges in %r for  %r should fail' % (value,key)
        logger.info('test %r', msg)
        invalid_input4 = library_screening_input.copy()
        invalid_input4[key] =  value
        errors, resp = self._create_resource(
            invalid_input4, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(find_in_dict(key, errors), 
            'test: %s, not in errors: %r' %(key,errors))

        logger.info('5. test with cherry_pick_source_plates')
        key = 'library_plates_screened' 
        plate_range1cpp = lps_format.format(**library_copy1_cpp).format(**library1)
        value = [ plate_range1cpp, plate_range2 ]
        msg = 'cherry_pick_source_plate copy type should fail %r' % value
        logger.info('test %r', msg)
        invalid_input5 = library_screening_input.copy()
        invalid_input5['library_plates_screened'] = value
        errors, resp = self._create_resource(
            invalid_input5, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(find_in_dict(key, errors), 
            'test: %s, not in errors: %r' %(key,errors))
        
        # 10 - "valid" test
        
        logger.info('10. Create valid library screening input...')
        logger.info('Patch batch_edit: cred: %r', self.username)
        resp = self.api_client.post(
            resource_uri,format='json', 
            data=library_screening_input, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 10.a Inspect meta "Result" section
        resp = self.deserialize(resp)
        logger.info('resp: %r', resp)

        self.assertTrue(API_RESULT_META in resp, '%r' % resp)
        self.assertTrue(API_MSG_RESULT in resp[API_RESULT_META], '%r' % resp)
        
        actual_result_meta = resp[API_RESULT_META][API_MSG_RESULT]
        expected_result_meta = {
            API_MSG_SCREENING_PLATES_UPDATED: 17,
            API_MSG_SCREENING_ADDED_PLATE_COUNT: 17, 
            API_MSG_SCREENING_DELETED_PLATE_COUNT: 0,
            API_MSG_SCREENING_EXTANT_PLATE_COUNT: 0,
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 17
        }
        for k,v in expected_result_meta.items():
            self.assertEqual(v,actual_result_meta[k], 
                'k: %r, expected: %r, actual: %r' % (k, v, actual_result_meta))
        
        self.assertTrue(len(resp[API_RESULT_DATA]) == 1)
        library_screening_output = resp[API_RESULT_DATA][0]
        
        # 10.b Inspect the returned library screening
        for k,v in library_screening_input.items():
            v2 = library_screening_output[k]
            self.assertTrue(equivocal(v,v2),
                'test key: %r:%r != %r' % (k, v, v2))

        # 10.c verify new plate volumes (just the first library)
        # retrieve plates from the server
        _data_for_get = { 
            'plate_number__range': 
                [library1['start_plate'],library1['end_plate']],
            'copy_name__eq': library_copy1['copy_name']
            }
        plate_resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            plate_resource_uri, format='json', 
            authentication=self.get_credentials(), data=_data_for_get)
        new_obj = self.deserialize(resp)
        expected_plates = library1['end_plate']-library1['start_plate']+1
        self.assertTrue(len(new_obj[API_RESULT_DATA]),expected_plates)
        
        expected_remaining_volume = (
            Decimal(library_copy1_input['initial_plate_well_volume']) - 
            Decimal(library_screening_input[
                'volume_transferred_per_well_from_library_plates']))
        for plate_data in new_obj[API_RESULT_DATA]:
            logger.info('inspect plate: %r', plate_data)
            self.assertEqual(
                Decimal(plate_data['remaining_well_volume']),
                expected_remaining_volume,
                'expected: %r, actual: %r' % (
                    expected_remaining_volume, plate_data))

        # 10.d. Test Screen statistics
        screen = self.get_screen(screen['facility_id'])

        key = 'library_screenings'
        expected_value = 1
        self.assertTrue(screen[key] == expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))

        key = 'screened_experimental_well_count'
        expected_value = experimental_well_count
        self.assertTrue(screen[key] == expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))
        key = 'unique_screened_experimental_well_count'
        expected_value = experimental_well_count
        self.assertTrue(screen[key] == expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))
        key = 'library_plates_screened'
        expected_value = 17 # 6 in library1, 11 in library2
        self.assertTrue(screen[key] == expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))
        key = 'assay_plates_screened'
        expected_value = 34 # lps * 2 replicates
        self.assertTrue(screen[key] == expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))

        # 10.e Logs
        
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'ref_resource_name': 'libraryscreening', 
            'key': library_screening_output['activity_id'],
            'api_action': API_ACTION_CREATE
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        logger.debug('logs: %r', apilogs)
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('apilog: %r', apilog)
        self.assertTrue(apilog['diff_keys'] is not None, 'no diff_keys' )
        self.assertTrue('library_plates_screened' in apilog['diff_keys'])
        self.assertTrue(plate_range1 in apilog['diffs'], 
            'plate_range1: %r not in diffs: %r'
            % (plate_range1, apilog['diffs']))
        
        self.assertTrue(apilog['child_logs'] == 17, 
            'wrong child_logs count: %r' % apilog)
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'parent_log_id': apilog['id']
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        self.assertTrue(
            len(apilogs) == 17, 
            'wrong number of child logs returned: %d: %r' 
            % (len(apilogs), apilogs))
        apilog = apilogs[0]
        logger.info('apilog: %r', apilog)
        self.assertTrue(apilog['diff_keys'] is not None, 'no diff_keys' )
        expected_remaining_volume = (
            Decimal(library_copy1_input['initial_plate_well_volume']) - 
            Decimal(library_screening_output[
                'volume_transferred_per_well_from_library_plates']))
        self.assertTrue('remaining_well_volume' in apilog['diff_keys'])
        diffs = json.loads(apilog['diffs'])
        logger.debug('diffs: %r', diffs)
        self.assertEqual(
            expected_remaining_volume, 
            Decimal(diffs['remaining_well_volume'][1]))
        
        # 11 Add a plate range:
        # Inspect PATCH logs after adding a plate range
        
        logger.info('11. Test - add a plate_range...')
        library_screening_input2 = { 
            'activity_id': library_screening_output['activity_id'],
            'library_plates_screened': 
                library_screening_output['library_plates_screened'] 
            }
        
        # add 6 more plates
        added_plate_range_start = library2['start_plate']+15
        added_plate_range_end = library2['start_plate']+20
        added_plate_range = lps_format.format(**library_copy2).format(
                start_plate=added_plate_range_start,
                end_plate=added_plate_range_end )
        logger.info('add plate range: %r', added_plate_range)
        library_screening_input2['library_plates_screened'].append(
            added_plate_range)
        logger.info('add plate range input: %r', library_screening_input2)
        resp = self.api_client.patch(
            resource_uri, 
            format='json', data=library_screening_input2, 
            authentication=self.get_credentials())
        
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        # 11.a Inspect meta "Result" section
        resp = self.deserialize(resp)
        logger.info('resp: %r', resp)

        self.assertTrue(API_RESULT_META in resp, '%r' % resp)
        self.assertTrue(API_MSG_RESULT in resp[API_RESULT_META], '%r' % resp)
        actual_result_meta = resp[API_RESULT_META][API_MSG_RESULT]
        expected_result_meta = {
            API_MSG_SCREENING_PLATES_UPDATED: 6,
            API_MSG_SCREENING_ADDED_PLATE_COUNT: 6, 
            API_MSG_SCREENING_DELETED_PLATE_COUNT: 0,
            API_MSG_SCREENING_EXTANT_PLATE_COUNT: 17,
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 23,
        }
        for k,v in expected_result_meta.items():
            self.assertEqual(v,actual_result_meta[k], 
                'k: %r, expected: %r, actual: %r' % (k, v, actual_result_meta))
        
        self.assertTrue(len(resp[API_RESULT_DATA]) == 1)
        library_screening_output2 = resp[API_RESULT_DATA][0]
        for k,v in library_screening_input2.items():
            v2 = library_screening_output2[k]
            self.assertTrue(equivocal(v,v2),
                'test key: %r:%r != %r' % (k, v, v2))
        
        # 11.c Inspect PATCH logs after adding a plate range
        logger.info('new plate ranges: %r', 
            library_screening_output2['library_plates_screened'])
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'ref_resource_name': 'libraryscreening', 
            'key': library_screening_output2['activity_id'],
            'api_action': API_ACTION_PATCH
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        logger.info('logs: %r', apilogs)
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('apilog: %r', apilog)
        self.assertTrue(apilog['diff_keys'] is not None, 'no diff_keys' )
        self.assertTrue('library_plates_screened' in apilog['diff_keys'])
        self.assertTrue(added_plate_range in apilog['diffs'], 
            'added_plate_range: %r not in diffs: %r'
            % (added_plate_range, apilog['diffs']))
        
        self.assertTrue(apilog['child_logs'] == 6, 
            'wrong child_logs count: %r' % apilog)
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'parent_log_id': apilog['id']
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        self.assertTrue(
            len(apilogs) == 6, 
            'wrong number of child logs returned: %d: %r' 
            % (len(apilogs), apilogs))
        apilog = apilogs[0]
        logger.debug('apilog: %r', apilog)
        self.assertTrue(apilog['diff_keys'] is not None, 'no diff_keys' )
        expected_remaining_volume = (
            Decimal(library_copy1_input['initial_plate_well_volume']) - 
            Decimal(library_screening_output2[
                'volume_transferred_per_well_from_library_plates']))
        self.assertTrue('remaining_well_volume' in apilog['diff_keys'])
        diffs = json.loads(apilog['diffs'])
        logger.info('diffs: %r', diffs)
        self.assertEqual(
            expected_remaining_volume, 
            Decimal(diffs['remaining_well_volume'][1]))
        
        # 11.e check a modified plate
        plate_resource_uri = BASE_URI_DB + '/librarycopyplate'
        modified_plate = self.get_single_resource(plate_resource_uri, {
            'copy_name': library_copy2['copy_name'],
            'plate_number': added_plate_range_start })
        logger.info('modified plate: %r', modified_plate)
        self.assertTrue(modified_plate is not None)
        self.assertEqual(
            expected_remaining_volume, 
            Decimal(modified_plate['remaining_well_volume']))
        
        logger.info('8. Test - delete the first plate_range (6 plates)..')
        logger.info('delete the plate: %r', 
            library_screening_output2['library_plates_screened'][0])
        plate_range_to_delete = \
            library_screening_output2['library_plates_screened'][0]
        library_screening_input3 = { 
            'activity_id': library_screening_output2['activity_id'],
            'library_plates_screened': 
                library_screening_output2['library_plates_screened'][1:] 
            }
        logger.info('input: %r', library_screening_input3)
        resp = self.api_client.patch(
            resource_uri, 
            format='json', data=library_screening_input3, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 11.a Inspect meta "Result" section
        resp = self.deserialize(resp)
        logger.info('resp: %r', resp)

        self.assertTrue(API_RESULT_META in resp, '%r' % resp)
        self.assertTrue(API_MSG_RESULT in resp[API_RESULT_META], '%r' % resp)
        actual_result_meta = resp[API_RESULT_META][API_MSG_RESULT]
        expected_result_meta = {
            API_MSG_SCREENING_PLATES_UPDATED: 6,
            API_MSG_SCREENING_ADDED_PLATE_COUNT: 0, 
            API_MSG_SCREENING_DELETED_PLATE_COUNT: 6,
            API_MSG_SCREENING_EXTANT_PLATE_COUNT: 17,
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 17
        }
        for k,v in expected_result_meta.items():
            self.assertEqual(v,actual_result_meta[k], 
                'k: %r, expected: %r, actual: %r' % (k, v, actual_result_meta))
        
        self.assertTrue(len(resp[API_RESULT_DATA]) == 1)
        library_screening_output3 = resp[API_RESULT_DATA][0]
        for k,v in library_screening_input3.items():
            v2 = library_screening_output3[k]
            self.assertTrue(equivocal(v,v2),
                'test key: %r:%r != %r' % (k, v, v2))

        # verify new plate volumes = initial_plate_well_volume
        # retrieve plates from the server
        _data_for_get = { 
            'plate_number__range': 
                [library1['start_plate'],library1['end_plate']],
            'copy_name__eq': library_copy1['copy_name']
            }
        plate_resource_uri = BASE_URI_DB + '/librarycopyplate'
        resp = self.api_client.get(
            plate_resource_uri, format='json', 
            authentication=self.get_credentials(), data=_data_for_get)
        new_obj = self.deserialize(resp)
        expected_plates = library1['end_plate']-library1['start_plate']+1
        self.assertTrue(len(new_obj[API_RESULT_DATA]),expected_plates)
        
        expected_remaining_volume = \
            Decimal(library_copy1_input['initial_plate_well_volume'])
        for plate_data in new_obj[API_RESULT_DATA]:
            self.assertEqual(
                Decimal(plate_data['remaining_well_volume']),
                expected_remaining_volume,
                'expected: %r, actual: %r' % (
                    expected_remaining_volume, plate_data))
            self.assertTrue(
                Decimal(plate_data['well_volume'])==expected_remaining_volume,
                'initial well volume expected: %r, actual: %r' % (
                    expected_remaining_volume, plate_data))
        
        # 12. test valid input with start_plate==end_plate (no end plate)
        
        logger.info('test valid single plate input...')

        single_plate_lps_format = '{library_short_name}:{copy_name}:{{start_plate}}'
        single_plate_lps_return_format = \
            '{library_short_name}:{copy_name}:{{start_plate}}-{{start_plate}}'
        library_plates_screened = [
            single_plate_lps_format.format(**library_copy1).format(**library1),
            single_plate_lps_format.format(**library_copy2).format(
                start_plate=library2['start_plate'])
        ]
        library_plates_screened_return_formatted = [
            single_plate_lps_return_format.format(
                **library_copy1).format(**library1),
            single_plate_lps_return_format.format(**library_copy2).format(
                start_plate=library2['start_plate'])
        ]

        library_screening_input4 = library_screening_input.copy()
        library_screening_input4['date_of_activity'] = '2016-08-01'
        # 200 nL requested (should be 7.9 uL)
        library_screening_input4['volume_transferred_per_well_from_library_plates'] = \
            '0.000000200'
        library_screening_input4['volume_transferred_per_well_to_assay_plates'] = \
            '0.000000200'
        library_screening_input4['number_of_replicates'] = 1
        library_screening_input4['library_plates_screened'] =  \
            library_plates_screened
        resp = self.api_client.post(
            resource_uri,format='json', 
            data=library_screening_input4, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        # 12.a Inspect meta "Result" section
        resp = self.deserialize(resp)

        self.assertTrue(API_RESULT_META in resp, '%r' % resp)
        self.assertTrue(API_MSG_RESULT in resp[API_RESULT_META], '%r' % resp)
        actual_result_meta = resp[API_RESULT_META][API_MSG_RESULT]
        expected_result_meta = {
            API_MSG_SCREENING_PLATES_UPDATED: 2,
            API_MSG_SCREENING_ADDED_PLATE_COUNT: 2, 
            API_MSG_SCREENING_DELETED_PLATE_COUNT: 0,
            API_MSG_SCREENING_EXTANT_PLATE_COUNT: 0,
            API_MSG_SCREENING_TOTAL_PLATE_COUNT: 2
        }
        for k,v in expected_result_meta.items():
            self.assertEqual(v,actual_result_meta[k], 
                'k: %r, expected: %r, actual: %r' % (k, v, actual_result_meta))

        self.assertTrue(len(resp[API_RESULT_DATA]) == 1)
        library_screening_output4 = resp[API_RESULT_DATA][0]
        for k,v in library_screening_input4.items():
            if k == 'library_plates_screened':
                continue # see next test, below
            v2 = library_screening_output4[k]
            self.assertTrue(equivocal(v,v2),
                'test key: %r:%r != %r' % (k, v, v2))
        logger.info('created library screening, with single-plate range: %r',
            library_screening_output4)
        found_count = 0
        for lps in library_screening_output4['library_plates_screened']:
            for lps_expected in library_plates_screened_return_formatted:
                if lps_expected==lps:
                    found_count += 1
        self.assertTrue(found_count==2, 
            'single plate test, returned ranges, expected: %r, returned: %r'
            %(library_plates_screened_return_formatted,
                library_screening_output4['library_plates_screened']))
        
        # Test Screen statistics

        screen = self.get_screen(screen['facility_id'])

        key = 'library_screenings'
        expected_value = 2
        self.assertTrue(screen[key] == expected_value,
            (key,'expected_value',expected_value,
                'returned value',screen[key]))
        
        # 10. test deletion of a library screening
        
        
        # TODO: test volume warning message:
        # 20170605 - JAS - allow insufficient volume        
        # logger.info('5.a. test insufficient volume...')
        # key = 'library_plates_screened' 
        # test_xfer_vol = \
        #     Decimal(library_copy1_input['initial_plate_well_volume']) - valid_test_volume
        # msg = ('insufficient volume {initial_plate_well_volume} '
        #     'for copy1: {copy_name} should fail').format(**library_copy1_input)
        # logger.info('test %r', msg)
        # invalid_input5 = library_screening_input.copy()
        # invalid_input5['volume_transferred_per_well_from_library_plates'] = \
        #     '{:.9f}'.format(test_xfer_vol) 
        # errors, resp = self._create_resource(
        #     invalid_input5, resource_uri, resource_test_uri, expect_fail=True)
        # self.assertTrue(resp.status_code==400, msg)
        # logger.info('errors: %r', errors)
        # self.assertTrue(find_in_dict(key, errors), 
        #     'test: %s, not in errors: %r' %(key,errors))
        # key2 = API_MSG_LCPS_INSUFFICIENT_VOLUME
        # self.assertTrue(find_in_dict(key2, errors), 
        #     'test: %s, not in errors: %r' %(key2,errors))

        
        
    def test3_create_publication(self):

        logger.info('test3_create_publication ...')
        
        self.screening_user = self.create_screensaveruser({ 
            'username': 'screening1'
        })
        logger.info('create screen...')        
        screen = self.create_screen({
            'screen_type': 'small_molecule'
            })
        
        publication_data = {
            'pubmed_id': 'PM00001',
            'pubmed_central_id': 12121212,
            'authors': 
                "Smith JA, White EA, Sowa ME, Powell ML, Ottinger M, Howley PM",
            'journal': (
                'Proceedings of the National Academy of Sciences of the '
                'United States of America'),
            'pages': "3752-7",
            'screen_facility_id': screen['facility_id'],
            'title': (
                "Genome-wide siRNA screen identifies SMCX, EP400, and Brd4 "
                "as E2-dependent regulators of human papillomavirus "
                "oncogene expression."),
            'volume': "107",
            'year_published': "2010"            
        }
        resource_uri = '/'.join([
            BASE_URI_DB, 'screen',screen['facility_id'],'publications'])
        # NOTE: post data because publications may include attached files
        content_type = MULTIPART_CONTENT
        resource_uri = '/'.join([
            BASE_URI_DB, 'screen',screen['facility_id'],'publications'])
        authentication=self.get_credentials()
        kwargs = {}
        kwargs['HTTP_AUTHORIZATION'] = authentication
        kwargs['HTTP_ACCEPT'] = JSON_MIMETYPE

        # Add an attached file and post the publication
        base_filename = 'iccbl_sm_user_agreement_march2015.pdf'
        filename = ('%s/db/static/test_data/useragreement/%s' 
            %(APP_ROOT_DIR,base_filename))
        logger.info('Open and POST file: %r', filename)
        with open(filename) as input_file:

            logger.info('POST publication with attached_file to the server...')
            publication_data['attached_file'] = input_file
            # NOTE: content_type arg is req'd with django.test.Client.post
            resp = self.django_client.post(
                resource_uri, content_type=content_type, 
                data=publication_data, **kwargs)
            self.assertTrue(
                resp.status_code == 201, 
                (resp.status_code,self.get_content(resp)))
        
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(),
            data={ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('new obj: %s ' % new_obj)
        self.assertEqual(
            len(new_obj[API_RESULT_DATA]),1,
            'wrong number of publications returned')
    
        publication_received = new_obj[API_RESULT_DATA][0]
        result,msgs = assert_obj1_to_obj2(
            publication_data, publication_received, 
            excludes=['attached_file'])
        self.assertTrue(result,msgs)
    
        # Check for the attached file
        uri = ('/db/publication/%s/attached_file' 
            % publication_received['publication_id'])
        try:
            # FIXME: create an admin user with publication/write privs
            admin_user = User.objects.get(username=self.admin_user['username'])
            view, args, kwargs = resolve(uri)
            kwargs['request'] = self.api_client.client.request()
            kwargs['request'].user=admin_user
            result = view(*args, **kwargs)
            self.assertTrue(
                result.status_code == 200, 
                'No attached file found, status_code: %r, %r' 
                    % (result.status_code, result))
            output_filename = '%s.out.%s' % tuple(filename.split('.'))
            logger.info('write %s to %r', filename, output_filename)
            with open(output_filename, 'w') as out_file:
                out_file.write(self.get_content(result))
            self.assertTrue(filecmp.cmp(filename,output_filename), 
                'input file: %r, not equal to output file: %r' 
                % (filename, output_filename))
            os.remove(output_filename)    
        except Exception, e:
            logger.info('no file found at: %r', uri)
            raise
        
        # Test apilog
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'publication', 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('publication log: %r', apilog)
        self.assertTrue(apilog['api_action'] == 'CREATE')
        self.assertTrue('attached_filename' in apilog['diff_keys'])
        self.assertTrue(base_filename in apilog['diffs'])
        
        return publication_received
        
    def test4_delete_publication(self):
        
        logger.info('test4_delete_publication...')
        
        publication_received = self.test3_create_publication()
        resource_uri = '/'.join([
            BASE_URI_DB, 'screen',publication_received['screen_facility_id'],
            'publications', str(publication_received['publication_id'])])
        
        logger.info('submit the publication for deletion...')
        resp = self.api_client.delete(
            resource_uri, authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code == 204, 
            (resp.status_code, self.get_content(resp)))
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(),
            data={ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code == 404, 
            ('error, publication should be deleted', resp.status_code, 
                self.get_content(resp)))

        # Check for the attached file record
        resource_uri = '/'.join([
            BASE_URI_DB, 'attached_file', 
            str(publication_received['attached_file_id'])])
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(),
            data={ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code == 404, 
            (resp.status_code, self.get_content(resp)))
        
        # Check for the attached file
        uri = ('/db/publication/%s/attached_file' 
            % publication_received['publication_id'])
        try:
            # FIXME: create an admin user with publication/write privs
            admin_user = User.objects.get(username=self.admin_user['username'])
            view, args, kwargs = resolve(uri)
            kwargs['request'] = self.api_client.client.request()
            kwargs['request'].user=admin_user
            result = view(*args, **kwargs)
            self.assertTrue(
                result.status_code == 404, 
                'Attached file found after delete, status_code: %r, %r' 
                    % (result.status_code, result))
        except Exception, e:
            logger.exception('exception when trying to locate: %r', uri)
            raise
    
    def test5_pin_transfer_approval(self):
        
        logger.info('test5_pin_transfer_approval...')
        # Create a Screen
        data = {
            'screen_type': 'small_molecule',
        }
        screen_item = self.create_screen(data=data)
        
        self.assertTrue(
            'facility_id' in screen_item, 
            'the facility_id was not created')
        
        for key, value in data.items():
            self.assertEqual(value, screen_item[key], 
                'key %r, val: %r not expected: %r' 
                % (key, value, screen_item[key]))

        logger.info('screen created: %r', screen_item)

        self.screening_user = self.create_screensaveruser({ 
            'username': 'screening1'
        })

        # 1. Set the pin_transfer data
        # FIXME: admin approved pin tranfer user only
        pin_transfer_data_expected = {
            'pin_transfer_approved_by_username': self.screening_user['username'],
            'pin_transfer_date_approved': _now().date().strftime("%Y-%m-%d"),
            'pin_transfer_comments': 'test pin_transfer_comments' }
        
        screen_update_data = {
            'facility_id': screen_item['facility_id']
            }
        screen_update_data.update(pin_transfer_data_expected)
        resource_uri = \
            BASE_URI_DB + '/screen/' + screen_update_data['facility_id']
        resp = self.api_client.patch(
            resource_uri, 
            format='json', data=screen_update_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        logger.info('get the updated pin_transfer screen data')
        new_screen_item = self.get_single_resource(resource_uri)
        logger.info('retrieved: %r', new_screen_item)
        for key,val in pin_transfer_data_expected.items():
            self.assertEqual(
                new_screen_item[key],pin_transfer_data_expected[key],
                'key: %r, %r, %r' 
                % (key,new_screen_item[key],pin_transfer_data_expected[key]))
            
        # 2. Modify the pin_transfer comment
        screen_update_data = {
            'facility_id': screen_item['facility_id']
            }
        screen_update_data['pin_transfer_comments'] = \
            'New test pin transfer comment'
        resp = self.api_client.patch(
            resource_uri, 
            format='json', data=screen_update_data,
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        logger.info('get the updated pin_transfer screen data')
        new_screen_item = self.get_single_resource(resource_uri)
        for key,val in pin_transfer_data_expected.items():
            if key == 'pin_transfer_comments':
                self.assertEqual(new_screen_item[key],screen_update_data[key])
            else:
                self.assertEqual(new_screen_item[key],
                    pin_transfer_data_expected[key],
                    'key: %r, %r, %r' 
                    % (key,new_screen_item[key],
                        pin_transfer_data_expected[key]))
    
    def test6_service_activity(self):
        logger.info('test6_service_activity...')
        # Create a Screen
        data = {
            'screen_type': 'small_molecule',
        }
        screen_item = self.create_screen(data=data)
        
        self.assertTrue(
            'facility_id' in screen_item, 
            'the facility_id was not created')
        
        # FIXME: performed_by_username belongs to ServiceActivityPerformers group
        performed_by_user = self.create_screensaveruser(
            { 'username': 'service_activity_performer'})

        service_activity_post = {
            'screen_facility_id': screen_item['facility_id'],
            'type': "image_analysis",
            'comments': "test",
            'date_of_activity': "2015-10-27",
            'funding_support': "clardy_grants",
            'performed_by_user_id': performed_by_user['screensaver_user_id'],
        }
        
        resource_uri = BASE_URI_DB + '/serviceactivity'
        resp = self.api_client.post(
            resource_uri, 
            format='json', 
            data=service_activity_post, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))

        data_for_get = { 'limit': 0, 'includes': ['*'] }
        resp = self.api_client.get(
            resource_uri,
            format='json', 
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        result,msgs = assert_obj1_to_obj2(
            service_activity_post, new_obj[API_RESULT_DATA][0])
        self.assertTrue(result,msgs)

        # Test apilog
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'serviceactivity', 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.debug('serviceactivity log: %r', apilog)
        self.assertTrue(apilog['api_action'] == 'CREATE')
        
        
class CherryPickRequestResource(DBResourceTestCase):
        
    def __init__(self, *args, **kwargs):
        DBResourceTestCase.__init__(self, *args, **kwargs)
        
    def setUp(self):
        super(CherryPickRequestResource, self).setUp()

    def tearDown(self):
        DBResourceTestCase.tearDown(self)
        logger.info('delete resources')
        Screen.objects.all().delete()
        Library.objects.all().delete()
        LibraryScreening.objects.all().delete()
        LabAffiliation.objects.all().delete()
        ApiLog.objects.all().delete()
        ScreensaverUser.objects.all().exclude(username='testsuper').delete()

    def _setup_duplex_data(self):
        logger.info('create users...')
        self.test_admin_user = self.create_screensaveruser(
            { 'username': 'adminuser',
              'permissions': 'resource/cherrypickrequest/write'  
            })

        logger.info('create library...')
        self.pool_library1 = self.create_library({
            'start_plate': 50000, 
            'end_plate': 50000,
            'screen_type': 'rnai',
            'is_pool': True })
        self.duplex_library1 = self.create_library({
            'start_plate': 50001, 
            'end_plate': 50004,
            'screen_type': 'rnai' })

        logger.info('set some experimental wells...')
        # Create the duplex library wells first
        logger.info('create duplex library %r wells...', 
            self.duplex_library1['short_name'])
        resource_name = 'well'
        input_well_data = []
        for plate in range(
            self.duplex_library1['start_plate'],
            self.duplex_library1['end_plate']+1):
            experimental_well_count = 384
            duplex_well_data = [
                self.create_test_well(
                    plate,i,library_well_type='experimental',
                    molar_concentration='0.001',
                    vendor_name='duplex_vendor2') 
                for i in range(0,experimental_well_count)]
            input_well_data.extend(duplex_well_data)
        logger.info(
            'patch duplex library %r wells...', 
            self.duplex_library1['short_name'])
        resource_uri = '/'.join([
            BASE_URI_DB,'library', 
            self.duplex_library1['short_name'],resource_name])
        resp = self.api_client.patch(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        self.duplex_wells = self.get_list_resource(resource_uri)
        self.duplex_wells = {
            well['well_id']:well for well in self.duplex_wells }

        # Create the pool library wells, using duplex wells
        logger.info(
            'create pool library %r wells...', self.pool_library1['short_name'])
        plate = 50000
        experimental_well_count = 384
        input_well_data = []
        for i in range(0,experimental_well_count):
            input_well = self.create_test_well(
                plate,i,library_well_type='experimental',
                molar_concentration='0.001',
                vendor_name='rna_vendor1') 
            duplex_wells = []
            for duplex_plate in range(
                    self.duplex_library1['start_plate'],
                    self.duplex_library1['end_plate']+1):
                duplex_wells.append(
                    lims_utils.well_id(duplex_plate, input_well['well_name']))
            input_well['duplex_wells'] = duplex_wells
            input_well_data.append(input_well)
        logger.info(
            'patch pool library %r wells...', self.pool_library1['short_name'])
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', 
            self.pool_library1['short_name'],resource_name])
        resp = self.api_client.patch(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        self.pool_wells = self.get_list_resource(resource_uri)
        for pool_well in self.pool_wells:
            self.assertTrue('duplex_wells' in pool_well)
            for duplex_well in pool_well['duplex_wells']:
                self.assertTrue(duplex_well in self.duplex_wells,
                    'duplex well not found: %r in %r' 
                    % (duplex_well, self.duplex_wells.keys()))
        self.pool_wells = { well['well_id']:well for well in self.pool_wells }
        
        logger.info('Create library copies...')
        
        # copy1: cherry_pick_source_plate
        duplex_library_copy1_input = {
            'library_short_name': self.duplex_library1['short_name'],
            'copy_name': "copy1",
            'usage_type': "cherry_pick_source_plates",
            'initial_plate_well_volume': '0.000010',
            'initial_plate_status': 'available'
        }  
        self.duplex_library_copy1 = self.create_copy(duplex_library_copy1_input)
        logger.info(
            'created duplex_library_copy1: %r', self.duplex_library_copy1)
 
        logger.info('create rnai screen...')        
        self.rnai_screen = self.create_screen({
            'screen_type': 'rnai'
            })

    def _setup_data(self):
        # Set up dependencies
        logger.info('create users...')
        self.test_admin_user = self.create_screensaveruser(
            { 'username': 'adminuser',
              'permissions': 'resource/cherrypickrequest/write'  
            })

        logger.info('create library...')
        self.library1 = self.create_library({
            'start_plate': 1000, 
            'end_plate': 1005,
            'screen_type': 'small_molecule' })
        self.library2 = self.create_library({
            'start_plate': 2000, 
            'end_plate': 2040,
            'screen_type': 'small_molecule' })

        # library 3 will contain alternate reagents for library1
        self.library3 = self.create_library({
            'start_plate': 3000, 
            'end_plate': 3005,
            'screen_type': 'small_molecule' })

        # library4: source wells will be forced to another destination plate,
        # if "keep_source_cherry_picks_together" is true 
        self.library4 = self.create_library({
            'start_plate': 4000, 
            'end_plate': 4001,
            'screen_type': 'small_molecule' })
        
        logger.info('set some experimental wells...')
        plate = 1000
        experimental_well_count = 384
        input_well_data = [
            self.create_small_molecule_test_well(
                plate,i,library_well_type='experimental',
                molar_concentration='0.001',
                vendor_name='vendor1') 
            for i in range(0,experimental_well_count)]
        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', self.library1['short_name'],resource_name])
        resp = self.api_client.patch(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        plate = 2001
        experimental_well_count = 384
        input_well_data = [
            self.create_small_molecule_test_well(
                plate,i,library_well_type='experimental',
                molar_concentration='0.001',
                vendor_name='vendor2') 
            for i in range(0,experimental_well_count)]
        resource_uri = '/'.join([
            BASE_URI_DB,'library', self.library2['short_name'],resource_name])
        resp = self.api_client.patch(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        plate = 3000
        experimental_well_count = 384
        input_well_data = [
            self.create_small_molecule_test_well(
                plate,i,library_well_type='experimental',
                molar_concentration='0.001',
                vendor_name='vendor1') 
            for i in range(0,experimental_well_count)]
        resource_uri = '/'.join([
            BASE_URI_DB,'library', self.library3['short_name'],resource_name])
        resp = self.api_client.patch(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        plate = 4000
        experimental_well_count = 384
        input_well_data = [
            self.create_small_molecule_test_well(
                plate,i,library_well_type='experimental',
                molar_concentration='0.001',
                vendor_name='vendorx') 
            for i in range(0,experimental_well_count)]
        resource_uri = '/'.join([
            BASE_URI_DB,'library', self.library4['short_name'],resource_name])
        resp = self.api_client.patch(
            resource_uri, format='sdf', data={ 'objects': input_well_data } , 
            authentication=self.get_credentials(), 
            **{ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        logger.info('Create library copies...')
        
        # copy1: cherry_pick_source_plate
        library_copy1_input = {
            'library_short_name': self.library1['short_name'],
            'copy_name': "copy1",
            'usage_type': "cherry_pick_source_plates",
            'initial_plate_well_volume': '0.000010',
            'initial_plate_status': 'available'
        }  
        self.library1_copy1 = self.create_copy(library_copy1_input)
        logger.info('created library copy1: %r', self.library1_copy1)
 
        # library1 copy1a: cherry_pick_source_plate
        resource_uri = BASE_URI_DB + '/librarycopy'
        library_copy1a_input = library_copy1_input.copy()
        library_copy1a_input['copy_name'] = 'copy1a'
        self.library1_copy1a = self.create_copy(library_copy1a_input)

        # library2, copy1: cherry_pick_source_plate
        library2_copy1_input = library_copy1_input.copy()
        library2_copy1_input['copy_name'] = 'library2copy1'
        library2_copy1_input['library_short_name'] = self.library2['short_name']
        self.library2_copy1 = self.create_copy(library2_copy1_input)
        logger.info('created library2 copy1: %r', self.library2_copy1)
        
        # library4,copy1: source wells are forced to other destination plate,
        # if "keep_source_cherry_picks_together" is true 
        library4_copy1_input = library_copy1_input.copy()
        library4_copy1_input['copy_name'] = 'library4copy1'
        library4_copy1_input['library_short_name'] = self.library4['short_name']
        self.library4_copy1 = self.create_copy(library4_copy1_input)
        logger.info('created library4 copy1: %r', self.library4_copy1)

        # Create extra copies:

        # copy3 - library_screening_plates - available
        # - also make the plates "available", so they won't be used for CPR
        # - "available" "library_screening_plates" should not show up in the
        # API_PARAM_SHOW_RETIRED_COPY_WELlS search 
        library_copy3_input = {
            'library_short_name': self.library1['short_name'],
            'copy_name': "copy3",
            'usage_type': "library_screening_plates",
            'initial_plate_well_volume': '0.000010',
            'initial_plate_status': 'available'
        }  
        self.library1_copy3 = self.create_copy(library_copy3_input)

        # copy4 - library_screening_plates - retired 
        # - also make the plates "retired", so they WILL be used for CPR
        # - "retired" "library_screening_plates" SHOULD show up in the
        # API_PARAM_SHOW_RETIRED_COPY_WELlS search 
        library_copy4_input = {
            'library_short_name': self.library1['short_name'],
            'copy_name': "copy4",
            'usage_type': "library_screening_plates",
            'initial_plate_well_volume': '0.000010',
            'initial_plate_status': 'retired'
        }  
        self.library1_copy4 = self.create_copy(library_copy4_input)

        logger.info('create screen...')        
        self.screen = self.create_screen({
            'screen_type': 'small_molecule'
            })
    
    def _create_cherry_pick_request(self, screen, data=None):
        # 1. Create the cherry pick object
        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        
        # FIXME: test requested by username is in screen
        
        # FIXME: create a "CherryPickRequestVolumeApprovers" group
        volume_approver = self.test_admin_user
        
        new_cpr_data = {
            'screen_facility_id': screen['facility_id'],
            # TODO: use a "CherryPickRequestAdmin"
            'requested_by_id': str(screen['collaborator_ids'][0]), 
            'date_requested': '2016-12-05',
            'transfer_volume_per_well_requested': '0.000000002', 
            'transfer_volume_per_well_approved': '0.000000002',
            'volume_approved_by_username': volume_approver['username'],
            'assay_plate_type': 'eppendorf_384',
            # wells to use: B03,C03,B04,C04,B05,C05
            'wells_to_leave_empty': (
                'Col:1, Col:2, Col:6, Col:7, Col:8, Col:9, Col:10, Col:11, '
                'Col:12, Col:13, Col:14, Col:15, Col:16, Col:17, Col:18, '
                'Col:19, Col:20, Col:21, Col:22, Col:23, Col:24, '
                'Row:A, Row:D, Row:E, Row:F, Row:G, Row:H, Row:I, Row:J, '
                'Row:K, Row:L, Row:M, Row:N, Row:O, Row:P')
            }
        if data is not None:
            new_cpr_data.update(data)
        # NOTE: wells_to_leave_empty leaves all but 6 cells
        resp = self.api_client.post(
            resource_uri,format='json', 
            data=new_cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        _data = self.deserialize(resp)
        new_cpr = _data[API_RESULT_DATA][0]
        logger.info('new cpr: %r', new_cpr)
        
        self.assertTrue('cherry_pick_request_id' in new_cpr, 
            'cherry_pick_request_id field missing: %r' % new_cpr)
        for key in new_cpr_data.keys():
            if key in ['transfer_volume_per_well_requested',
                'transfer_volume_per_well_approved' ]:
                self.assertEqual(
                    Decimal(new_cpr_data[key]), 
                    Decimal(new_cpr[key]))
            else:
                self.assertEqual(
                    str(new_cpr_data[key]), str(new_cpr[key]),
                    '%s %r!=%r' % (key, new_cpr_data[key], new_cpr[key]))
        
        return new_cpr
        
    def _get_scps(self, cpr_id, data_for_get=None):

        _data_for_get={ 'limit': 0 }
        if data_for_get:
            _data_for_get.update(data_for_get)
        resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_id),
            'screener_cherry_pick'])
        return self.get_list_resource(resource_uri, data_for_get)
    
    def _get_lcps(self, cpr_id, data_for_get=None):
        
        _data_for_get={ 'limit': 0, 'includes': '*' }
        if data_for_get:
            _data_for_get.update(data_for_get)
        resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_id),
            'lab_cherry_pick'])
        return self.get_list_resource(resource_uri, data_for_get)
        
    def _get_cpaps(self, cpr_id, data_for_get=None):
        
        _data_for_get={ 'limit': 0, 'includes': '*' }
        if data_for_get:
            _data_for_get.update(data_for_get)
            
        resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_id),
            'cherry_pick_plate'])
        return self.get_list_resource(resource_uri, data_for_get)
    
    def _test1_cherry_pick_request(self):
        
        logger.info('test1_cherry_pick_request')
        self._setup_data()
        return self._create_cherry_pick_request(
            self.screen,
            data={ 'keep_source_plate_cherry_picks_together': False})
    
    def _test2a_set_rnai_screener_cherry_picks(self):
        logger.info('test2a_set_rnai_screener_cherry_picks')
        self._setup_duplex_data()
        cpr_data = self._create_cherry_pick_request(self.rnai_screen)
    
        # 2. Set Screener Cherry Picks
        # TODO: test input by plate-well-range
        
        well_id_template = '{start_plate}:%s'.format(**self.pool_library1)
        screener_cherry_picks = [
            well_id_template % 'A01',
            well_id_template % 'A02',
            well_id_template % 'A03',
            well_id_template % 'A04',
            well_id_template % 'A05',
            well_id_template % 'A06',
            ]
        expected_screener_cherry_picks = screener_cherry_picks[:]
            
        cpr_data['screener_cherry_picks'] = screener_cherry_picks
        
        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 2.a check the Screener Cherry Picks
        scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        self.assertEqual(
            len(expected_screener_cherry_picks), len(scp_well_data))
        
        for well_id in expected_screener_cherry_picks:
            found = False
            for scp_well in scp_well_data:
                if scp_well['screened_well_id'] == well_id:
                    found = True
                    break
            self.assertTrue(
                found, 'well_id: %r not found in %r'
                % (well_id, scp_well_data)) 
        
        return (cpr_data, scp_well_data)
    
    def _test2_set_screener_cherry_picks(self):
        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        cpr_data = self._test1_cherry_pick_request()

        # 2. Set Screener Cherry Picks
        # TODO: test input by plate-well-range
        
        well_id_template = '01000:%s'
        screener_cherry_picks = [
            well_id_template % 'A01',
            well_id_template % 'A02',
            well_id_template % 'A03',
            well_id_template % 'A04',
            well_id_template % 'A05',
            well_id_template % 'A06',
            ]
        expected_screener_cherry_picks = screener_cherry_picks[:]
        screener_cherry_picks.append('2001:P19 P20 P21')
        expected_screener_cherry_picks.extend([
            '02001:P19','02001:P20', '02001:P21'])
        screener_cherry_picks.append('04000 A01,B01,C01,D01,E01,P23')
        expected_screener_cherry_picks.extend([
            '04000:A01','04000:B01','04000:C01','04000:D01','04000:E01',
            '04000:P23',])

        cpr_data['screener_cherry_picks'] = screener_cherry_picks
        
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 2.a check the Screener Cherry Picks
        scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        self.assertEqual(
            len(expected_screener_cherry_picks), len(scp_well_data))
        scp_well_data = {scp['screened_well_id']:scp for scp in scp_well_data}
        for well_id in expected_screener_cherry_picks:
            self.assertTrue(well_id in scp_well_data,
                'well_id: %r not found in scp_well_data: %r'
                % (well_id, scp_well_data.keys()))
        return (
            cpr_data, scp_well_data.values(), expected_screener_cherry_picks)

    def test2a_set_screener_cherry_picks_alternate_searches(self):
        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        cpr_data = self._test1_cherry_pick_request()
        well_id_template = '01000 %s'
        screener_cherry_picks = [
            well_id_template % 'A01',
            well_id_template % 'A02',
            well_id_template % 'A03',
            well_id_template % 'A04',
            well_id_template % 'A05 A06',
            ]
        well_id_template = '01000:%s'
        expected_screener_cherry_picks = [
            well_id_template % 'A01',
            well_id_template % 'A02',
            well_id_template % 'A03',
            well_id_template % 'A04',
            well_id_template % 'A05',
            well_id_template % 'A06',
            ]
        cpr_data['screener_cherry_picks'] = screener_cherry_picks
        
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 2.a check the Screener Cherry Picks
        scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        self.assertEqual(
            len(expected_screener_cherry_picks), len(scp_well_data))
        scp_well_data = {scp['screened_well_id']:scp for scp in scp_well_data}
        for well_id in expected_screener_cherry_picks:
            self.assertTrue(well_id in scp_well_data,
                'well_id: %r not found in scp_well_data: %r'
                % (well_id, scp_well_data.keys()))
            
        # 3. Delete Screener Cherry Picks
        cpr_data['screener_cherry_picks'] = None
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('response: %r', _data)
        _meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_SCPS_DELETED in _meta)
        self.assertEqual(_meta[API_MSG_SCPS_DELETED],len(scp_well_data))
        scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        self.assertEqual(len(scp_well_data),0)
        
        # 4. Alternate patterns:
        cpr_data['screener_cherry_picks'] = '01000:A01 A02 A03 A04 A05 A06'
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 4.a check the Screener Cherry Picks
        scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        self.assertEqual(
            len(expected_screener_cherry_picks), len(scp_well_data))
        scp_well_data = {scp['screened_well_id']:scp for scp in scp_well_data}
        for well_id in expected_screener_cherry_picks:
            self.assertTrue(well_id in scp_well_data,
                'well_id: %r not found in scp_well_data: %r'
                % (well_id, scp_well_data.keys()))
        # 4.b delete    
        cpr_data['screener_cherry_picks'] = None
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        cpr_data['screener_cherry_picks'] = '01000:A01,A02,A03,A04,A05,A06'
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 4.a check the Screener Cherry Picks
        scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        self.assertEqual(
            len(expected_screener_cherry_picks), len(scp_well_data))
        scp_well_data = {scp['screened_well_id']:scp for scp in scp_well_data}
        for well_id in expected_screener_cherry_picks:
            self.assertTrue(well_id in scp_well_data,
                'well_id: %r not found in scp_well_data: %r'
                % (well_id, scp_well_data.keys()))
        # 4.c delete    
        cpr_data['screener_cherry_picks'] = None
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        cpr_data['screener_cherry_picks'] = '01000 A01,A02,A03,A04,A05,A06'
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        # 4.a check the Screener Cherry Picks
        scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        self.assertEqual(
            len(expected_screener_cherry_picks), len(scp_well_data))
        scp_well_data = {scp['screened_well_id']:scp for scp in scp_well_data}
        for well_id in expected_screener_cherry_picks:
            self.assertTrue(well_id in scp_well_data,
                'well_id: %r not found in scp_well_data: %r'
                % (well_id, scp_well_data.keys()))
        
        
    def _modify_copy_well_volume(self, copy_data, volume_adjustment, well_id):
        
        logger.info('_modify_copy_well_volume: %r, adj: %r, %r',
            copy_data, volume_adjustment, well_id)
        cw_resource_uri = '/'.join([
            BASE_URI_DB,'library',copy_data['library_short_name'],'copy',
            copy_data['copy_name'],'copywell',well_id])
        cw = self.get_single_resource(cw_resource_uri)
        original_volume = Decimal(cw['volume'])
        cw['volume'] = original_volume - Decimal(volume_adjustment)
        resp = self.api_client.patch(
            cw_resource_uri, format='json', data=cw, 
            authentication=self.get_credentials() )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        patch_response = self.deserialize(resp)
        
        self.assertTrue(API_RESULT_DATA in patch_response)
        self.assertEqual(len(patch_response[API_RESULT_DATA]), 1)        
        new_copywell = patch_response[API_RESULT_DATA][0]
        logger.info('adjusted copywell: %r', new_copywell)
        self.assertEqual(
            Decimal(volume_adjustment), 
            Decimal(new_copywell['consumed_volume']))
        self.assertEqual(
            Decimal(cw['volume']), Decimal(new_copywell['volume']))
        self.assertEqual(
            original_volume, Decimal(new_copywell['initial_volume']))
        
        return new_copywell
    
    def _submit_screener_cherry_picks(
        self, cpr_id, lcp_target='set_lab_cherry_picks'):
        '''
        Submit the **already created** screener_cherry_picks as lab_cherry_picks
        return the lab_cherry_picks that were created
        '''
        
        lcp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_id),
            lcp_target])
        resp = self.api_client.post(
            lcp_resource_uri, format='json',
            data = {},
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        meta = _data[API_RESULT_META]
        logger.info(
            'response to set_lab_cherry_picks: %r',  meta)
        # TODO: verify the response metadata
        
        # 3.a Get the lab_cherry_picks that were created
        lcp_well_data = self._get_lcps(cpr_id)        

        return (lcp_well_data, meta)
    
    def _test3_set_lab_cherry_picks(self):
        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        (cpr_data, scp_well_data,expected_scps) = \
            self._test2_set_screener_cherry_picks()
        
        # Prep:
        # Modify copy1:A01 volume so that it will not be chosen
        # (copy1a will be chosen instead)
        self._modify_copy_well_volume(
            self.library1_copy1, '0.000010', '01000:A01')

        # 3 Submit the screener_cherry_picks -> lab_cherry_picks
        (lcp_well_data,meta) = \
            self._submit_screener_cherry_picks(
                cpr_data['cherry_pick_request_id'])
        
        self.assertTrue(API_MSG_WARNING not in meta, 
            'meta warning: %r' % meta)
        if API_MSG_LCPS_UNFULFILLED in meta:
            self.assertEqual(meta[API_MSG_LCPS_UNFULFILLED], 0)
        
        # 3.b verify the lab cherry picks
        scp_well_data = { scp['screened_well_id']:scp for scp in scp_well_data}
        # - verify that copy1a was chosen instead of copy1 due to well 01000:A01
        library_copy_map = {
            self.library1['short_name']: self.library1_copy1a['copy_name'],
            self.library2['short_name']: self.library2_copy1['copy_name'],
            self.library4['short_name']: self.library4_copy1['copy_name']
            }
        for lcp_well in lcp_well_data:
            logger.debug('test lcp: %r', lcp_well)
            source_well_id = lcp_well['source_well_id']
            assigned_library = lcp_well['library_short_name']
            assigned_copy = lcp_well['source_copy_name']
            self.assertTrue(source_well_id in scp_well_data,
                'source well: %r not found in %r' 
                    % (source_well_id, scp_well_data.keys()))
            self.assertTrue(assigned_library in library_copy_map,
                'well: %r, assigned library %r not in expected: %r'
                % (source_well_id,assigned_library,library_copy_map))
            self.assertEqual(assigned_copy,library_copy_map[assigned_library])
        
        return (cpr_data, lcp_well_data)
    
    def test_3b_set_duplex_lab_cherry_picks(self):
        (cpr_data, scp_well_data) = \
            self._test2a_set_rnai_screener_cherry_picks()
        
        # 3 Submit the screener_cherry_picks -> lab_cherry_picks
        (lcp_well_data,meta) = \
            self._submit_screener_cherry_picks(
                cpr_data['cherry_pick_request_id'],
                lcp_target='set_duplex_lab_cherry_picks')
        self.assertEqual(len(scp_well_data)*4, len(lcp_well_data))
        
        lcp_well_data = {lcp['source_well_id']:lcp for lcp in lcp_well_data}
        
        # 3.b verify the lab cherry picks
        for scp in scp_well_data:
            well_id = scp['screened_well_id']
            pool_well = self.pool_wells[well_id]
            duplex_wells = pool_well['duplex_wells']
            logger.info('scp pool well: %r, %r', well_id, duplex_wells)
            for duplex_well in duplex_wells:
                self.assertTrue(duplex_well in lcp_well_data)
                lcp_well = lcp_well_data[duplex_well]
                self.assertEqual(
                    lcp_well['library_short_name'], 
                    self.duplex_library1['short_name'])
                self.assertEqual(
                    lcp_well['source_copy_name'], 
                    self.duplex_library_copy1['copy_name'])
                self.assertEqual(lcp_well['screener_well_id'],well_id)
        return (cpr_data, lcp_well_data)

    def test_4b_reserve_map_keep_source_plates_together_false(self):
        (cpr_data, lcp_well_data) = self._test3_set_lab_cherry_picks()
        cpr_id = cpr_data['cherry_pick_request_id']
        # Modify the CPR.keep_source_plates_together
        cpr_data['keep_source_plate_cherry_picks_together'] = False
        del cpr_data['screener_cherry_picks']
        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        resp = self.api_client.patch(
            resource_uri,format='json', 
            data=cpr_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        _data = self.deserialize(resp)
        logger.info('cpr patch result: %r', _data)
        cpr_data = self.get_single_resource(resource_uri)
        logger.info('new cpr: %r', cpr_data)
        self.assertTrue(
            cpr_data['keep_source_plate_cherry_picks_together']==False)
        # Map
        # 4. Map/reserve source copies (that are fulfilled) (allocates volumes)
        resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'reserve_map_lab_cherry_picks'])
        resp = self.api_client.post(
            resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)

        logger.info('data: %r', _data)

        self.assertTrue(API_RESULT_META in _data)
        _meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_LCP_PLATES_ASSIGNED in _meta)
        copy_plate_assigned_msg = _meta[API_MSG_LCP_PLATES_ASSIGNED]
        plate_copy1a = '1000:%s' % self.library1_copy1a['copy_name']
        plate_copy_l2_c1 = '2001:%s' % self.library2_copy1['copy_name']
        plate_copy_l4_c1 = '4000:%s' % self.library4_copy1['copy_name']
        
        for msg in copy_plate_assigned_msg:
            if plate_copy1a in msg:
                self.assertEqual(6, msg[1])
            elif plate_copy_l2_c1 in msg:
                self.assertEqual(3, msg[1])
            elif plate_copy_l4_c1 in msg:
                self.assertEqual(6, msg[1])
            else:
                self.fail(
                    'unknown platecopy has been assigned: %r, '
                    'expected platecopies: %r'
                    % (msg,[plate_copy1a,plate_copy_l2_c1,plate_copy_l4_c1]))
        
        self.assertTrue(API_MSG_LCP_ASSAY_PLATES_CREATED in _meta )
        expected_assay_plates = 3
        self.assertTrue(
            len(_meta[API_MSG_LCP_ASSAY_PLATES_CREATED]),expected_assay_plates)
        
        plated_lab_cherry_picks = self._get_lcps(cpr_id)
        
        # check the assay plate/well assignments
        # confirm new plate mapping
        copyplate_to_assayplates = defaultdict(set)
        for lcp in plated_lab_cherry_picks:
            copy_plate = '{source_copy_name}/{library_plate}'.format(**lcp)
            copyplate_to_assayplates[copy_plate].add(
                lcp['cherry_pick_plate_number'])
        logger.info('copyplate_to_assayplates: %r', copyplate_to_assayplates)
        expected_assignments = {
            '%s/1000'%self.library1_copy1a['copy_name']: [1],
            '%s/2001'%self.library2_copy1['copy_name']: [2],
            # would only be 3 if keep_source_plate_cherry_picks_together == True
            '%s/4000'%self.library4_copy1['copy_name']: [2,3],
            }
        for copy_plate in expected_assignments:
            actual_assignment = copyplate_to_assayplates[copy_plate]
            expected_assignment = set(expected_assignments[copy_plate])
            logger.info(
                'copy_plate: %r, %r, %r', 
                copy_plate, actual_assignment, expected_assignment)
            self.assertEqual(expected_assignment,actual_assignment)

        self._validate_plate_mapping(
            plated_lab_cherry_picks, cpr_data['wells_to_leave_empty'],
            is_random=False)
        
        lcps_by_plate = defaultdict(list)
        for lcp in plated_lab_cherry_picks:
            copy_well = '{source_copy_name}/{source_well_id}'.format(**lcp)
            lcps_by_plate[lcp['cherry_pick_plate_number']].append(
                (copy_well,lcp['destination_well']))
        # TODO: verify the actual well assignments
        for plate,copywells in lcps_by_plate.items():
            logger.info('%r: %r', plate, copywells) 
            
        # TODO: verify random/non-random layout (random at this point)
        
    def test_4_reserve_map_lab_cherry_picks(self):
        '''
        Tests the simple-case:
        - no overrides needed
        - copy plates are available
        '''
        
        (cpr_data, lcp_well_data) = self._test3_set_lab_cherry_picks()
        cpr_id = cpr_data['cherry_pick_request_id']

        previous_lab_cherry_picks = self._get_lcps(cpr_id)
        previous_lab_cherry_picks = { lcp['source_well_id']: lcp 
            for lcp in previous_lab_cherry_picks }
        # 4. Map/reserve source copies (that are fulfilled) (allocates volumes)
        resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'reserve_map_lab_cherry_picks'])
        resp = self.api_client.post(
            resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)

        logger.info('data: %r', _data)

        self.assertTrue(API_RESULT_META in _data)
        _meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_LCP_PLATES_ASSIGNED in _meta)
        copy_plate_assigned_msg = _meta[API_MSG_LCP_PLATES_ASSIGNED]
        plate_copy1a = '1000:%s' % self.library1_copy1a['copy_name']
        plate_copy_l2_c1 = '2001:%s' % self.library2_copy1['copy_name']
        plate_copy_l4_c1 = '4000:%s' % self.library4_copy1['copy_name']
        
        for msg in copy_plate_assigned_msg:
            if plate_copy1a in msg:
                self.assertEqual(6, msg[1])
            elif plate_copy_l2_c1 in msg:
                self.assertEqual(3, msg[1])
            elif plate_copy_l4_c1 in msg:
                self.assertEqual(6, msg[1])
            else:
                self.fail(
                    'unknown platecopy has been assigned: %r, '
                    'expected platecopies: %r'
                    % (msg,[plate_copy1a,plate_copy_l2_c1,plate_copy_l4_c1]))
        
        self.assertTrue(API_MSG_LCP_ASSAY_PLATES_CREATED in _meta )
        expected_assay_plates = 3
        self.assertTrue(
            len(_meta[API_MSG_LCP_ASSAY_PLATES_CREATED]),expected_assay_plates)
        
        plated_lab_cherry_picks = self._get_lcps(cpr_id)
        self._validate_plate_mapping(
            plated_lab_cherry_picks, cpr_data['wells_to_leave_empty'],
            is_random=False)
        
        # 4.A check the assay plate/well assignments
        # - use the lcp.copywell fields to check copywell adjustments
        assay_plates = set()
        for lcp in plated_lab_cherry_picks:
            logger.debug('lcp: %r', lcp)
            assay_plates.add(lcp['cherry_pick_plate_number'])
            # Ensure that each source copy gets a different assay_plate
            # only if keep_source_plate_cherry_picks_together
            # source_copy_id = lcp['source_copy_id']
            # cpp_number = lcp['cherry_pick_plate_number']
            # if source_copy_id not in copy_to_assay_plate:
            #     copy_to_assay_plate[source_copy_id] = cpp_number
            # else:
            #     self.assertEqual(
            #         cpp_number,
            #         copy_to_assay_plate[source_copy_id],
            #         'lcp: %r, copy wells split between plates: %r and %r'
            #         % (lcp, cpp_number, copy_to_assay_plate[source_copy_id]))

            # Check that the copy-wells have had their volumes adjusted
            transfer_volume_per_well_approved = \
                Decimal(cpr_data['transfer_volume_per_well_approved'])
            self.assertEqual(
                Decimal(lcp['source_copy_well_consumed_volume']),
                transfer_volume_per_well_approved, 
                'lcp vol consumed should be %r, %r'
                % ( transfer_volume_per_well_approved, lcp) )
            
            previous_lcp = previous_lab_cherry_picks[lcp['source_well_id']]
            expected_vol = (
                Decimal(previous_lcp['source_copy_well_volume'])
                    -transfer_volume_per_well_approved)
            self.assertEqual(
                expected_vol,
                Decimal(lcp['source_copy_well_volume']),
                'lcp well volume reported should be %r, %r'
                % ( expected_vol, lcp) )
            
        # 4.b Check that copywells have the correct information
        cw_resource_uri = '/'.join([
            BASE_URI_DB, 'copywell'])
        copy_wells_to_get = [lcp['source_copywell_id']
            for lcp in plated_lab_cherry_picks]
        data_for_get = { 'copywell_id__in': copy_wells_to_get }
        cws = self.get_list_resource(cw_resource_uri, data_for_get) 
        self.assertEqual(len(cws), len(plated_lab_cherry_picks))
        for cw in cws:
            self.assertEqual(cw['cherry_pick_screening_count'],1, 
                'cw: {copywell_id}, cpsc: {cherry_pick_screening_count}'.format(**cw))
        
        # 4.c Check that the plates have the correct information
        
        platecopy_resource_uri = '/'.join([BASE_URI_DB, 'librarycopyplate'])
        plates_to_get = set([
            '%s/%s' % (lcp['source_copy_name'],str(lcp['library_plate']))
            for lcp in plated_lab_cherry_picks])
        data_for_get = { 'copyplate_id__in': [x for x in plates_to_get] }
        logger.info('get the plates: %r', plates_to_get)
        plates = self.get_list_resource(platecopy_resource_uri, data_for_get) 
        self.assertEqual(len(plates), len(plates_to_get))
        for plate in plates:
            self.assertEqual(plate['cplt_screening_count'],1, 
                'plate: {copyplate_id}, cpsc: {cplt_screening_count}'.format(**plate))
        
        # Expect 3 assay plates to hold the 6+3+6=15 lcp's with 6 spaces per plate
        self.assertTrue(len(assay_plates) == expected_assay_plates, 
            'wrong number of assay_plates created: %r' % assay_plates)
        
        cpaps = self._get_cpaps(cpr_id)
        logger.info('cpaps created: %r', cpaps)
        self.assertEqual(len(cpaps),expected_assay_plates)
        
        # Test that the cpr.date_volume_reserved has been set
        resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id'])])
        new_cpr_data = self.get_single_resource(resource_uri)
        expected_date = _now().date().strftime("%Y-%m-%d")
        self.assertEqual(new_cpr_data['date_volume_reserved'],expected_date)

        # TODO: Test logs:
        # - copy well logs
        # - copyplate logs
        # - cpr parent log
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'cherrypickrequest', 
            'key': cpr_data['cherry_pick_request_id'],
            'diff_keys': 'date_volume_reserved' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )

        try:
            logger.info('logs: %r', apilogs)
            self.assertTrue(
                len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
            apilog = apilogs[0]
            self.assertTrue('date_volume_reserved' in apilog['diffs'])
            self.assertTrue(expected_date in apilog['diffs'])
            
        except AssertionError:
            logger.exception(
                'logs: %r, expected_date: %r', apilogs, expected_date)
            raise
        expected_child_logs = 18 #  15 copy well logs, 3 plate logs
        self.assertTrue(apilog['child_logs'], expected_child_logs)
        data_for_get={ 
            'limit': 0, 
            'includes': ['*'],
            'parent_log_id': apilog['id']
        }
        apilogs = self.get_list_resource(
            BASE_REPORTS_URI + '/apilog', 
            data_for_get=data_for_get )
        self.assertEqual(
            len(apilogs),expected_child_logs, 
            'wrong number of child logs returned: %d: %r' 
            % (len(apilogs), apilogs))
        
        for apilog in apilogs:
            logger.info('apilog: %r, %r, %r', 
                apilog['uri'],apilog['ref_resource_name'],apilog['key'])
            if apilog['ref_resource_name'] == 'copywell':
                self.assertTrue('volume' in apilog['diffs'], 
                    'apilog.diffs: %r'% apilog)
                # TODO: verify volumes
                self.assertTrue('cherry_pick_screening_count' 
                    in apilog['diffs'], 'apilog.diffs: %r'% apilog)
                
            if apilog['ref_resource_name'] == 'librarycopyplate':
                self.assertTrue('cplt_screening_count' in apilog['diffs'], 
                    'apilog.diffs: %r'% apilog)
        
        return (cpr_data, lcp_well_data)
        
    def _test_2a_change_screener_selections(self):    
        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        (cpr_data, scp_well_data,expected_scps) = \
            self._test2_set_screener_cherry_picks()

        # Test that using the API_PARAM_SHOW_OTHER_REAGENTS includes the (6) 
        # reagents from library3
        data_for_get={ 
            'limit': 0,
            API_PARAM_SHOW_OTHER_REAGENTS: True }
        expanded_scp_well_data = \
            self._get_scps(
                cpr_data['cherry_pick_request_id'],
                data_for_get=data_for_get)
        logger.debug('found expanded_scps: %r', expanded_scp_well_data)
        expected_count = len(scp_well_data) + 6
        self.assertEqual(expected_count, len(expanded_scp_well_data))
        
        for scp in scp_well_data:
            well_id = scp['screened_well_id']
            foundSelected = False
            foundExpanded = False
            for scp_well in expanded_scp_well_data:
                if scp_well['screened_well_id'] == well_id:
                    foundSelected = True
                elif scp_well['searched_well_id'] == well_id:
                    foundExpanded = True
                    self.assertEqual(
                        scp_well['library_short_name'],
                        self.library3['short_name'])
            self.assertTrue(
                foundSelected, 'selected well_id: %r not found in %r'
                % (well_id, scp_well_data)) 
            if scp['library_short_name'] == self.library1['short_name']:
                self.assertTrue(
                    foundExpanded, 'expanded well_id: %r not found in %r'
                    % (well_id, scp_well_data)) 
        
        # Now select two of the expanded wells
        other_reagents_to_post = []
        for scp_well in expanded_scp_well_data:
            if scp_well['library_short_name'] == self.library3['short_name']:
                scp_well['selected'] = True
                other_reagents_to_post.append(scp_well)
            if len(other_reagents_to_post) == 2:
                break

        scp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'screener_cherry_pick'])
        
        resp = self.api_client.patch(
            scp_resource_uri, 
            format='json', 
            data=other_reagents_to_post, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('patch screener_cherry_picks: %r', _data)
        
        _meta = _data[API_RESULT_META]
        created_scps = _meta[API_MSG_SCP_CREATED]
        unselected_scps = _meta[API_MSG_SCP_UNSELECTED]
        for changed_scp in other_reagents_to_post:
            self.assertTrue(changed_scp['screened_well_id'] in created_scps,
                '%r not in %r' %(changed_scp,created_scps))
            self.assertTrue(changed_scp['searched_well_id'] in unselected_scps,
                '%r not in %r' %(changed_scp,unselected_scps))
        
        # Retrieve the new SCP's
        new_scp_well_data = self._get_scps(cpr_data['cherry_pick_request_id'])
        logger.debug('found new_scp_well_data: %r', new_scp_well_data)
        expected_count = len(scp_well_data) + 2 # for the two new selections
        self.assertEqual(expected_count, len(new_scp_well_data))
        
        for scp in new_scp_well_data:
            screened_well_id = scp['screened_well_id']
            if screened_well_id in created_scps:
                self.assertTrue(scp['selected'], 'should be selected: %r'%scp)
            elif screened_well_id in unselected_scps:
                self.assertFalse(scp['selected'], 
                    'should be un selected: %r'%scp)
            else:
                self.assertTrue(scp['selected'], 
                    'all other scps should be selected: %r'%scp)
        
        # Check the cpr stats:
        cpr_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id'])
            ])
        new_cpr_data = self.get_single_resource(cpr_resource_uri)
        logger.info('new cpr data: %r', new_cpr_data)
        expected_scp_count = len(scp_well_data)
        self.assertEqual(
            new_cpr_data['screener_cherry_pick_count'],
            expected_scp_count )
        
        # TODO: check for a CPR.screener_cherry_pick_count log
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'cherrypickrequest', 
            'key': cpr_data['cherry_pick_request_id'],
            'diff_keys': 'screener_cherry_pick_count' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.debug('logs: %r', apilogs)
        expected_count = 3 # create, set SCPs, update SCPs
        self.assertEqual(
            len(apilogs),expected_count, 
            'wrong number of apilogs found: %r' % apilogs)
        apilog = apilogs[1]
        self.assertTrue('screener_cherry_pick_count' in apilog['diffs'], 
            'diffs: %r' % apilog['diffs'])
        apilog = apilogs[2]
        self.assertTrue('screener_cherry_pick_count' in apilog['diffs'], 
            'diffs: %r' % apilog['diffs'])

        return (cpr_data, new_scp_well_data)
    
    def test_2b_change_screener_selections_after_lcp(self):
        ''' 
        - Test that lcp's must be deleted for scp's to be reassigned after 
        lcp's have been selected.
        '''
        
        (cpr_data, scp_well_data) = self._test_2a_change_screener_selections()
        cpr_id = cpr_data['cherry_pick_request_id']

        selected_scp_well_data = self._get_scps(
            cpr_data['cherry_pick_request_id'],
            data_for_get={ 'selected': True })
        expected_count = 6 + 3 + 6
        self.assertEqual(len(selected_scp_well_data), expected_count)
        
        # 1 Submit the screener_cherry_picks -> lab_cherry_picks
        (lcp_well_data,meta) = \
            self._submit_screener_cherry_picks(cpr_id)
        self.assertEqual(len(selected_scp_well_data), len(lcp_well_data))
        
        scp_well_data = {scp['screened_well_id']:scp for scp in scp_well_data}
        
        # 2. Try to switch wells from library3 back to library1
        # - verify error when LCP's are assigned
        other_reagents_to_post = []
        unselected_reagents = []
        for screened_well_id, scp_well in scp_well_data.items():
            if scp_well['library_short_name'] == self.library3['short_name']:
                other_well = scp_well_data[scp_well['searched_well_id']]
                other_well['selected'] = True
                other_reagents_to_post.append(other_well)
                unselected_reagents.append(scp_well)
        scp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'screener_cherry_pick'])
        resp = self.api_client.patch(
            scp_resource_uri, 
            format='json', 
            data=other_reagents_to_post, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('patch screener_cherry_picks error: %r', _data)
        _errors = _data[API_RESULT_ERROR]
        self.assertTrue(API_MSG_NOT_ALLOWED in _errors)
        self.assertTrue(API_MSG_LCPS_MUST_BE_DELETED in _errors)
        
        # 3 delete lab cherry picks
        delete_lcps_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'delete_lab_cherry_picks'])
        resp = self.api_client.post(
            delete_lcps_resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        _meta = _data[API_RESULT_META]
        logger.info('meta: %r', _meta)
        self.assertTrue(API_MSG_LCPS_REMOVED in _meta)
        expected_delete_count = len(lcp_well_data)
        self.assertEqual(
            _meta[API_MSG_LCPS_REMOVED], expected_delete_count)

        # 3.a verify LCPs dne
        lcp_well_data = self._get_lcps(cpr_id)        
        self.assertEqual(len(lcp_well_data), 0)

        # 3.b verify that SCPs still exist
        selected_scp_well_data = self._get_scps(
            cpr_data['cherry_pick_request_id'],
            data_for_get={ 'selected': True })
        self.assertEqual(len(selected_scp_well_data), expected_count)
        
        # 4 verify that scps can be changed (back to original library from test 2)
        resp = self.api_client.patch(
            scp_resource_uri, 
            format='json', 
            data=other_reagents_to_post, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        _meta = _data[API_RESULT_META]
        reselected_scps = _meta[API_MSG_SCP_RESELECTED]
        unselected_scps = _meta[API_MSG_SCP_UNSELECTED]
        for changed_scp in other_reagents_to_post:
            self.assertTrue(changed_scp['screened_well_id'] in reselected_scps,
                '%r not in %r' %(changed_scp,reselected_scps))
        for unselected_scp in unselected_reagents:
            self.assertTrue(
                unselected_scp['screened_well_id'] in unselected_scps,
                '%r not in %r' %(unselected_scp,unselected_scps))
        
        # TODO: test delete screener_cherry_picks
        
    def _test_3a_change_lab_cherry_picks(self):
        ''' 
        Like test_3, 
        1. but now with unfulfillable wells that will have to be resolved,
        2. with extra Copies that can be shown and chosen, override required
        '''
        
        (cpr_data, scp_well_data,expected_scps) = \
            self._test2_set_screener_cherry_picks()
        cpr_id = cpr_data['cherry_pick_request_id']
        
        # Prep:
        # A. Make library1 unavailable by making well 01000:A01 unfulfillable
        # - for all copies:
        # Modify copy1:A01 volume so that it will not be chosen
        # (well will be unfulfilled instead) for both copies
#         self._modify_copy_well_volume(
#             self.library1_copy1, '0.000000010', '01000:A01')
#         self._modify_copy_well_volume(
#             self.library1_copy1a, '0.000000010', '01000:A01')
        self._modify_copy_well_volume(
            self.library1_copy1, '0.000010', '01000:A01')
        self._modify_copy_well_volume(
            self.library1_copy1a, '0.000010', '01000:A01')

        # 1 Submit the screener_cherry_picks -> lab_cherry_picks
        (lcp_well_data,meta) = \
            self._submit_screener_cherry_picks(cpr_id)
        self.assertEqual(len(scp_well_data), len(lcp_well_data))
        
        # 1.A verify the lab cherry picks
        scp_well_data = { scp['screened_well_id']:scp for scp in scp_well_data}
        library_copy_map = {
            self.library1['short_name']: self.library1_copy1['copy_name'],
            self.library2['short_name']: self.library2_copy1['copy_name'],
            self.library4['short_name']: self.library4_copy1['copy_name']
            }
        # - verify well 01000:A01 unavailability causes "unfulfilled"
        lcp_well_1_a1 = None
        lcp_well_1_a2 = None
        for lcp_well in lcp_well_data:
            logger.debug('test lcp: %r', lcp_well)
            source_well_id = lcp_well['source_well_id']
            assigned_library = lcp_well['library_short_name']
            assigned_copy = lcp_well['source_copy_name']
            self.assertTrue(source_well_id in scp_well_data,
                'source well: %r not found in %r' 
                    % (source_well_id, scp_well_data.keys()))
            if lcp_well['library_short_name'] == self.library1['short_name']:
                if source_well_id == '01000:A01':
                    logger.debug(
                        'lcp well should be unfulfilled: %r', lcp_well)
                    self.assertEqual(
                        VOCAB_LCP_STATUS_UNFULFILLED, lcp_well['status'])
                    self.assertEqual(lcp_well['source_copy_name'], None)
                    lcp_well_1_a1 = lcp_well
                    continue
                else:
                    if lcp_well['source_well_id'] == '01000:A02':
                        lcp_well_1_a2 = lcp_well
                    self.assertEqual(
                        VOCAB_LCP_STATUS_SELECTED, lcp_well['status'])
            self.assertTrue(assigned_library in library_copy_map,
                'well: %r, assigned library %r not in expected: %r'
                % (source_well_id,assigned_library,library_copy_map))
            self.assertEqual(assigned_copy,library_copy_map[assigned_library])
        
        # 1.B Verify that the API_PARAM_SHOW_COPY_WELLS settings shows
        # the other "cherry_pick_source_plate" copy- "available" plates when 
        # available (but not the "retired" or "available" 
        # "library_screening_plates" copy-plates)
        data_for_get = { API_PARAM_SHOW_COPY_WELLS: True }
        expanded_lcp_data = self._get_lcps(cpr_id, data_for_get)

        number_wells_copy1 = 6
        number_wells_copy1a = 6
        expected_count = len(scp_well_data) + number_wells_copy1a
        self.assertEqual(len(expanded_lcp_data),expected_count)
        
        # 1.C force selection of 01000:A01 -> copy1
        copy1_name = self.library1_copy1['copy_name']
        copy4_name = self.library1_copy4['copy_name']
        logger.info('Force selection of well 01000:A01 to %r, '
            '(requires override next step, plating)...', copy1_name)
        lcp_well_1_a1['source_copy_name'] = copy1_name
        lcp_well_1_a1['selected'] = True
        lcp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'lab_cherry_pick'])
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a1], 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        
        # Verify meta
        _meta = _data[API_RESULT_META]
        logger.info('lcp patch meta: %r', _meta)

        selected_lcps = _meta[API_MSG_LCP_SELECTED]
        deselected_lcps = _meta[API_MSG_LCP_DESELECTED]
        changed_lcps = _meta[API_MSG_LCP_CHANGED]
        self.assertEqual(len(selected_lcps),1, 'selected: %r' % selected_lcps)
        self.assertEqual(len(deselected_lcps),0)
        self.assertEqual(len(changed_lcps),0)
        
        self.assertTrue(
            LCP_COPYWELL_KEY.format(**lcp_well_1_a1) in selected_lcps)
        
        self.assertTrue(API_MSG_WARNING in _meta, 
            'missing warning: %r' % API_MSG_WARNING)
        self.assertTrue(API_MSG_LCPS_INSUFFICIENT_VOLUME 
            in _meta[API_MSG_WARNING])
        insufficient_volume_cws = \
            _meta[API_MSG_WARNING][API_MSG_LCPS_INSUFFICIENT_VOLUME]
        
        cw_name =  LCP_COPYWELL_KEY.format(**lcp_well_1_a1)
        self.assertTrue(cw_name in insufficient_volume_cws,
            'insufficient cws: %r' % insufficient_volume_cws)
        
        # Verify that all lcps are "selected" now
        new_lcps = self._get_lcps(cpr_id)
        for lcp_well in new_lcps:
            self.assertEqual(
                VOCAB_LCP_STATUS_SELECTED, lcp_well['status'])

        # 2.A Verify that the API_PARAM_SHOW_RETIRED_COPY_WELlS shows
        # (A) plus the "retired" "libary_screening_plates"
        data_for_get = { API_PARAM_SHOW_RETIRED_COPY_WELlS: True }
        expanded_lcp_data2 = self._get_lcps(cpr_id, data_for_get)

        number_wells_copy1 = 6
        number_wells_copy1a = 6
        number_wells_copy4 = 6
        expected_count = \
            len(scp_well_data) + number_wells_copy1a + number_wells_copy4
        self.assertEqual(len(expanded_lcp_data2),expected_count)
        
        # 2.B Switch already selected well 01000:A02 to copy4
        logger.info('Force selection of well 01000:A02 to %r, '
            'library_screening_plates, '
            '(requires override next step, plating)...', copy4_name)
        lcp_well_1_a2['source_copy_name'] = copy4_name
        lcp_well_1_a2['source_copy_id'] = None
        lcp_well_1_a2['selected'] = True
        lcp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'lab_cherry_pick'])
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a2], 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        
        # Verify meta
        _meta = _data[API_RESULT_META]
        logger.info('lcp patch meta: %r', _meta)

        selected_lcps = _meta[API_MSG_LCP_SELECTED]
        deselected_lcps = _meta[API_MSG_LCP_DESELECTED]
        changed_lcps = _meta[API_MSG_LCP_CHANGED]

        self.assertEqual(len(selected_lcps),0)
        self.assertEqual(len(deselected_lcps),0)
        self.assertEqual(len(changed_lcps),1, 'changed: %r' % changed_lcps)
        
        self.assertTrue(
            LCP_COPYWELL_KEY.format(**lcp_well_1_a2) in changed_lcps[0])
        
        self.assertTrue(API_MSG_WARNING in _meta, 
            'missing warning: %r' % API_MSG_WARNING)
        self.assertTrue(
            API_MSG_LCPS_INSUFFICIENT_VOLUME in _meta[API_MSG_WARNING])
        insufficient_volume_cws = \
            _meta[API_MSG_WARNING][API_MSG_LCPS_INSUFFICIENT_VOLUME]
        
        new_lcps = self._get_lcps(cpr_id)
        
        for lcp_well in new_lcps:
            if lcp_well['source_well_id'] == '01000:A01':
                self.assertEqual(lcp_well['source_copy_name'],copy1_name)
            if lcp_well['source_well_id'] == '01000:A02':
                self.assertEqual(lcp_well['source_copy_name'],copy4_name)
        
        # 3. Test fail: try to select two copies for the same source_well_id
        lcp_well_1_a2_second = lcp_well_1_a2.copy()
        lcp_well_1_a2_second['source_copy_name'] = copy1_name 
        lcp_well_1_a2_second['source_copy_id'] = None
        lcp_well_1_a2_second['selected'] = True
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a2,lcp_well_1_a2_second], 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        self.assertTrue(API_RESULT_ERROR in _data)
        _errors = _data[API_RESULT_ERROR]
        self.assertTrue(API_MSG_LCP_MULTIPLE_SELECTIONS_SUBMITTED in _errors)
        _errors = _errors[API_MSG_LCP_MULTIPLE_SELECTIONS_SUBMITTED]
        self.assertEqual(len(_errors),1)
        self.assertTrue(lcp_well_1_a2['source_well_id'] in _errors[0],
             '%r not in errors: %r' 
                % (lcp_well_1_a2['source_well_id'], _errors[0]))
        # Verify meta
        
        
        return (cpr_data, new_lcps)
    
    def _validate_plate_mapping(
        self, lcps, wells_to_leave_empty, is_random=False):
        
        plate_size = 384
        available_wells = lims_utils.assay_plate_available_wells(
            wells_to_leave_empty, plate_size)

        # Check the well assignments (non-random):
        # cherry_pick_plate_number and destination well_name should increment
        # in this ordering (if not keep_source_plate_cherry_picks_together):
        lcps_by_plate_copy_well_dest_well = [
            (lcp['library_plate'],lcp['source_copy_name'],
                lcp['source_well_id'],lcp['cherry_pick_plate_number'],
                lcp['destination_well']) 
                for lcp in lcps]
            
        lcps_by_plate_copy_well_dest_well = sorted(
            lcps_by_plate_copy_well_dest_well)
        
        previous_index = 0
        previous_plate_number = 0
        platesize = 384
        random_assigments = []
        for lcp_plating in lcps_by_plate_copy_well_dest_well:
            logger.info('lcp plating: %r', lcp_plating)
            assay_plate_number = lcp_plating[3]
            destination_well = lcp_plating[4]
            
            self.assertTrue(destination_well in available_wells,
                'destination_well: %r not in available_wells: %r'
                % (destination_well, available_wells))
            
            self.assertTrue(assay_plate_number >= previous_plate_number,
                'assay_plate_number: %d !>= previous: %d'
                % (assay_plate_number, previous_plate_number))
            if assay_plate_number > previous_plate_number:
                previous_plate_number = assay_plate_number
                previous_index = 0
                if is_random is True:
                    self.assertTrue(random_assigments!=sorted(random_assigments),
                        'random assigments are not random: %r' % random_assigments)
                random_assigments = []
                
            current_index = lims_utils.index_from_well_name(
                destination_well, platesize)
            # (knowing the first available well) test start point as well
            if previous_index == 0 and is_random is False:
                self.assertEqual(destination_well,available_wells[0])
            if is_random is False and previous_index > 0:
                self.assertTrue(
                    current_index > previous_index,
                    'previous_index: %d ! > current_index: %d, %r'
                    % (previous_index, current_index, lcp_plating))
            if is_random is True:
                random_assigments.append(current_index)
            previous_index = current_index
        
    
    def test_4a_update_reservation_and_mapping(self):
        
        (cpr_data, current_lcps) = self._test_3a_change_lab_cherry_picks()
        cpr_id = cpr_data['cherry_pick_request_id']
        
        current_lcps = { lcp['source_well_id']:lcp for lcp in current_lcps }
        
        # 1. Reserve and plate
        # 1.A verify override required for 01000:A01 - insufficient volume
        # - was set to 0 in previous test
        # Map/reserve the source copies (that are fulfilled) (allocates volumes)
        plating_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'reserve_map_lab_cherry_picks'])
        resp = self.api_client.post(
            plating_resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        
        self.assertTrue(API_RESULT_ERROR in _data)
        errors = _data[API_RESULT_ERROR]
        self.assertTrue('transfer_volume_per_well_approved' in errors)
        copywell_1000_a01 = LCP_COPYWELL_KEY.format(**current_lcps['01000:A01'])
        self.assertTrue(API_MSG_LCPS_INSUFFICIENT_VOLUME in errors )
        volume_errors = errors[API_MSG_LCPS_INSUFFICIENT_VOLUME]
        logger.info('volume_errors: %r', volume_errors)
        self.assertTrue(len(volume_errors)==1)
        self.assertTrue(
            copywell_1000_a01 in volume_errors[0] )
        self.assertTrue(API_PARAM_VOLUME_OVERRIDE in errors)
        
        # 1.B Repeat volume with override for well A01
        resource_uri = \
            plating_resource_uri + '?' + API_PARAM_VOLUME_OVERRIDE + '=true'
        resp = self.api_client.post(
            resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)

        # 1.B.1 meta data checks
        logger.info('override resp: %r', _data)
        self.assertTrue(API_RESULT_META in _data)
        _meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_LCPS_VOLUME_OVERRIDDEN in _meta)
        volume_overrides = _meta[API_MSG_LCPS_VOLUME_OVERRIDDEN]
        logger.info('volume_overrides: %r', volume_overrides)
        self.assertTrue(len(volume_overrides)==1)
        self.assertTrue(
            copywell_1000_a01 in volume_overrides[0])

        self.assertTrue(API_MSG_LCP_PLATES_ASSIGNED in _meta)
        copy_plate_assigned_msg = _meta[API_MSG_LCP_PLATES_ASSIGNED]
        expected_copyplate_assignments = {
            '1000:%s'%self.library1_copy1['copy_name']: 5,
            '1000:%s'%self.library1_copy4['copy_name']: 1,
            '2001:%s'%self.library2_copy1['copy_name']: 3,
            '4000:%s'%self.library4_copy1['copy_name']: 6
            }
        logger.info(
            'check expected_copyplate_assignments: %r', 
            expected_copyplate_assignments)
        logger.info('actual: %r', copy_plate_assigned_msg)
        for copyname,assigned_count in expected_copyplate_assignments.items():
            found = False
            for msg in copy_plate_assigned_msg:
                if copyname in msg:
                    self.assertEqual(
                        expected_copyplate_assignments[copyname],msg[1])

        self.assertTrue(API_MSG_LCP_ASSAY_PLATES_CREATED in _meta )
        expected_assay_plates = 3
        self.assertTrue(
            len(_meta[API_MSG_LCP_ASSAY_PLATES_CREATED]),expected_assay_plates)

        # 1.C verify the copy-wells have had their volumes adjusted
        lcp_well_data = self._get_lcps(cpr_id)        
        lcp_well_data = {lcp['source_well_id']:lcp for lcp in lcp_well_data}
        transfer_volume_per_well_approved = \
            Decimal(cpr_data['transfer_volume_per_well_approved'])
        for source_well_id, lcp in lcp_well_data.items():
            self.assertEqual(lcp['status'], VOCAB_LCP_STATUS_PLATED)
            if lcp['source_well_id'] != '01000:A01':
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume']),
                    transfer_volume_per_well_approved, 
                    'lcp vol consumed should be %r, %r'
                    % ( transfer_volume_per_well_approved, lcp) )
            else:
                # well 01000:A01: previous consumed = 10mL
                expected_vol = \
                    Decimal('0.000010000') + transfer_volume_per_well_approved
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume']),
                    expected_vol, 
                    'lcp vol consumed should be %r, %r'
                    % ( expected_vol, lcp) )
        
        # 1.D validate plate mapping ordering
        self._validate_plate_mapping(
            lcp_well_data.values(), cpr_data['wells_to_leave_empty'],
            is_random=False)
                
        # 20170227 - Modification:
        # Cannot PASTE selection changes to /lab_cherry_pick after plated,
        # Can PASTE selection changes to /plate_mapping after plated
        # - if API_PARAM_SET_DESELECTED_TO_ZERO==True, 
        # set deselected copy-well vols to zero 
        # - wipe out LCP's: still requires cancel plating

        # 2.A Try to change the LCP assignment after plating with /lab_cherry_pick:
        # - Fails because already plated
        # - Fails because override required, OR, reservation must be canceled first
        
        copy4_name = self.library1_copy4['copy_name']
        logger.info(
            'Try to Force selection of well 01000:A03 to %r', copy4_name)
        lcp_well_1_a3 = current_lcps['01000:A03']
        lcp_well_1_a3['source_copy_name'] = copy4_name
        lcp_well_1_a3['source_copy_id'] = None
        lcp_well_1_a3['selected'] = True
        lcp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'lab_cherry_pick'])
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a3], 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        
        self.assertTrue(API_RESULT_ERROR in _data)
        self.assertTrue(API_PARAM_OVERRIDE in _data[API_RESULT_ERROR])
#         self.assertTrue(API_MSG_LCP_PLATES_ASSIGNED in _data[API_RESULT_ERROR])

        # 2.C cancel plating reservation
        cpaps = self._get_cpaps(cpr_id)
        cancel_reservation_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'cancel_reservation'])
        resp = self.api_client.post(
            cancel_reservation_resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        _meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_COPYWELLS_DEALLOCATED in _meta)
        expected_deallocated_count = len(lcp_well_data)
        self.assertEqual(
            _meta[API_MSG_COPYWELLS_DEALLOCATED], expected_deallocated_count)
        expected_assay_plates_removed = len(cpaps)
        self.assertEqual(
            _meta[API_MSG_CPR_ASSAY_PLATES_REMOVED], 
            expected_assay_plates_removed)

        # 2.D Verify copy well volumes are reset
        # - but all the copy assignments are still set
        new_lab_cherry_picks = self._get_lcps(cpr_id)
        new_lab_cherry_picks = {
            LCP_COPYWELL_KEY.format(**lcp):lcp 
            for lcp in new_lab_cherry_picks}
        for source_well_id, lcp in new_lab_cherry_picks.items():
            logger.info('lcp: %r', source_well_id)
            self.assertTrue(lcp['source_copy_name'] != None, 'lcp: %r' % lcp)
            self.assertTrue(
                lcp.get('cherry_pick_plate_number',"")==None,
                'cherry_pick_plate_number: %r' % lcp)
            self.assertTrue(lcp.get('destination_well',None)==None,
                'destination well: %r' % lcp)
            # Test that the copy-wells have had their volumes deallocated
            if lcp['source_well_id'] != '01000:A01':
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume']),0,
                    'lcp vol consumed should be 0, %r'% lcp)
                self.assertEqual(
                    Decimal(lcp['source_copy_well_initial_volume']),
                    Decimal(lcp['source_copy_well_volume']),
                    'lcp source_copy_well_initial_volume '
                    '!= source_copy_well_volume: %r'
                    % lcp )
            else:
                # 01000:A01 was set to 0 in test3 prep
                self.assertEqual(
                    Decimal(lcp['source_copy_well_volume']), 0, 
                    'source_copy_well_volume 0, %r, %r'
                    % (source_well_id, lcp))
                # well 01000:A01: previous consumed = 10mL
                expected_vol = Decimal('0.000010000')
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume']),
                    expected_vol, 
                    'lcp vol consumed should be %r, %r'
                    % (expected_vol,lcp))
                
        # 2.D1 Check that copywells have the correct information
        cw_resource_uri = '/'.join([
            BASE_URI_DB, 'copywell'])
        copy_wells_to_get = [lcp['source_copywell_id']
            for lcp in new_lab_cherry_picks.values()]
        data_for_get = { 'copywell_id__in': copy_wells_to_get }
        cws = self.get_list_resource(cw_resource_uri, data_for_get) 
        self.assertEqual(len(cws), len(new_lab_cherry_picks))
        for cw in cws:
            self.assertEqual(cw['cherry_pick_screening_count'],0, 
                'cw: {copywell_id}, cpsc: {cherry_pick_screening_count}'.format(**cw))
        
        # 2.D2 Check that the plates have the correct information
        
        platecopy_resource_uri = '/'.join([BASE_URI_DB, 'librarycopyplate'])
        plates_to_get = set([
            '%s/%s' % (lcp['source_copy_name'],str(lcp['library_plate']))
            for lcp in new_lab_cherry_picks.values()])
        data_for_get = { 'copyplate_id__in': [x for x in plates_to_get] }
        logger.info('get the plates: %r', plates_to_get)
        plates = self.get_list_resource(platecopy_resource_uri, data_for_get) 
        self.assertEqual(len(plates), len(plates_to_get))
        for plate in plates:
            self.assertEqual(plate['cplt_screening_count'],0, 
                'plate: {copyplate_id}, cpsc: {cplt_screening_count}'.format(**plate))
                
        # 2.E Verify copy well logs 
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'cherrypickrequest', 
            'key': cpr_id,
            'diff_keys': 'date_volume_reserved' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.debug('date_volume_reserved logs: %r', apilogs)
        
        # Find the deallocate log:
        deallocate_log = None
        for apilog in apilogs:
            diffs = json.loads(apilog['diffs'])
            logger.info('diffs: %r', diffs)
            if diffs['date_volume_reserved'][1] is None:
                deallocate_log = apilog
        data_for_get={ 
            'limit': 0, 
            'parent_log_id': deallocate_log['id']
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.debug('deallocate child logs: %r', apilogs)
        expected_deallocate_logs = 19 # 6+3+6 copywell, 4 plate
        self.assertEqual(expected_deallocate_logs, len(apilogs))
        for apilog in apilogs:
            if apilog['ref_resource_name'] == 'copywell':
                self.assertTrue('volume' in apilog['diffs'], 
                    'apilog.diffs: %r'% apilog)
                # TODO: verify volumes
                self.assertTrue('cherry_pick_screening_count' 
                    in apilog['diffs'], 'apilog.diffs: %r'% apilog)
                
                self.assertTrue(API_MSG_PLATING_CANCELED in apilog['comment'])
            if apilog['ref_resource_name'] == 'librarycopyplate':
                self.assertTrue('cplt_screening_count' in apilog['diffs'], 
                    'apilog.diffs: %r'% apilog)

        # 2.C Verify the LCP's can be altered now
        # NOTE: override not req'd until plating
        logger.info('Try 2 Force selection of well 01000:A03 to copy4,')
        lcp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'lab_cherry_pick'])
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a3], 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        # TODO: verify the response metadata
        # Verify warning about insufficient well vol A01
        
        # Verify selection
        new_lab_cherry_picks = self._get_lcps(cpr_id)
        new_lab_cherry_picks = {
            lcp['source_well_id']:lcp for lcp in new_lab_cherry_picks}
        test_well_id = lcp_well_1_a3['source_well_id']
        self.assertTrue(test_well_id in new_lab_cherry_picks)
        self.assertEqual(
            lcp_well_1_a3['source_copy_name'], 
            new_lab_cherry_picks[test_well_id]['source_copy_name'])
        
        # 2.E Try to de-select an LCP
        logger.info('2.E Try to de-select an LCP')
        test_data = lcp_well_1_a3.copy()
        test_data['selected'] = False
        lcp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'lab_cherry_pick'])
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[test_data], 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.debug('data: %r', _data)
        new_lab_cherry_picks = self._get_lcps(cpr_id)
        new_lab_cherry_picks = {
            lcp['source_well_id']:lcp for lcp in new_lab_cherry_picks}
        test_well_id = lcp_well_1_a3['source_well_id']
        self.assertTrue(test_well_id in new_lab_cherry_picks)
        logger.info('new lcp: %r', new_lab_cherry_picks[test_well_id])
        self.assertTrue(
            new_lab_cherry_picks[test_well_id]['source_copy_name']==None)
        
        # 3.E1 reselect
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a3], 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        # Verify selection
        new_lab_cherry_picks = self._get_lcps(cpr_id)
        new_lab_cherry_picks = {
            lcp['source_well_id']:lcp for lcp in new_lab_cherry_picks}
        test_well_id = lcp_well_1_a3['source_well_id']
        self.assertTrue(test_well_id in new_lab_cherry_picks)
        self.assertEqual(
            lcp_well_1_a3['source_copy_name'], 
            new_lab_cherry_picks[test_well_id]['source_copy_name'])

        
        # Part 3:
        # Plate again and then use "cancel_reservation"
        # (Note: 01000:A01 volume override needed)
        # (Note: 01000:A03 now also copy4)
        resource_uri = \
            plating_resource_uri + '?' + API_PARAM_VOLUME_OVERRIDE + '=true'
        resp = self.api_client.post(
            resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        self.assertTrue(API_RESULT_META in _data)
        _meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_LCP_PLATES_ASSIGNED in _meta)
        copy_plate_assigned_msg = _meta[API_MSG_LCP_PLATES_ASSIGNED]
        expected_copyplate_assignments = {
            '1000:%s'%self.library1_copy1['copy_name']: 4,
            '1000:%s'%self.library1_copy4['copy_name']: 2,
            '2001:%s'%self.library2_copy1['copy_name']: 3,
            '4000:%s'%self.library4_copy1['copy_name']: 6
            }
        logger.info(
            'check expected_copyplate_assignments: %r', 
            expected_copyplate_assignments)
        logger.info('actual: %r', copy_plate_assigned_msg)
        for copyname,assigned_count in expected_copyplate_assignments.items():
            found = False
            for msg in copy_plate_assigned_msg:
                if copyname in msg:
                    expected_count = expected_copyplate_assignments[copyname]
                    self.assertEqual(expected_count,msg[1],
                        'copy: %r expected: %r, %r' 
                            % (copyname, expected_count, msg))
        self.assertTrue(API_MSG_LCP_ASSAY_PLATES_CREATED in _meta )
        expected_assay_plates = 3
        self.assertTrue(
            len(_meta[API_MSG_LCP_ASSAY_PLATES_CREATED]),expected_assay_plates)

        # 3.A1 verify the copy-wells have had their volumes adjusted
        lcp_well_data = self._get_lcps(cpr_id)        
        lcp_well_data = {lcp['source_well_id']:lcp for lcp in lcp_well_data}
        transfer_volume_per_well_approved = \
            Decimal(cpr_data['transfer_volume_per_well_approved'])
        for source_well_id, lcp in lcp_well_data.items():
            self.assertEqual(lcp['status'], VOCAB_LCP_STATUS_PLATED)
            if lcp['source_well_id'] != '01000:A01':
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume']),
                    transfer_volume_per_well_approved, 
                    'lcp vol consumed should be %r, %r'
                    % ( transfer_volume_per_well_approved, lcp) )
            else:
                # well 01000:A01: previous consumed = 10mL
                expected_vol = \
                    Decimal('0.000010000') + transfer_volume_per_well_approved
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume']),
                    expected_vol, 
                    'lcp vol consumed should be %r, %r'
                    % ( expected_vol, lcp) )
        
        # 3.A2 validate plate mapping ordering
        self._validate_plate_mapping(
            lcp_well_data.values(), cpr_data['wells_to_leave_empty'],
            is_random=False)

        # 3.B cancel plating resrvation
        cancel_reservation_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'cancel_reservation'])
        resp = self.api_client.post(
            cancel_reservation_resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        _meta = _data[API_RESULT_META]
        logger.info('meta: %r', _meta)
        self.assertTrue(API_MSG_COPYWELLS_DEALLOCATED in _meta)
        expected_deallocated_count = len(lcp_well_data)
        self.assertEqual(
            _meta[API_MSG_COPYWELLS_DEALLOCATED], expected_deallocated_count)
        expected_assay_plates_removed = 3
        self.assertEqual(
            _meta[API_MSG_CPR_ASSAY_PLATES_REMOVED], 
            expected_assay_plates_removed)
        
        # 3.B1 verify copy-wells are restored
        new_lab_cherry_picks = self._get_lcps(cpr_id)
        new_lab_cherry_picks = {
            lcp['source_well_id']:lcp for lcp in new_lab_cherry_picks}
        for source_well_id, lcp in new_lab_cherry_picks.items():
            logger.debug('lcp: %r', source_well_id)
            self.assertEqual(lcp['status'], VOCAB_LCP_STATUS_SELECTED)
            self.assertTrue(lcp['source_copy_name'] != None, 'lcp: %r' % lcp)
            self.assertTrue(
                lcp.get('cherry_pick_plate_number',"")==None,
                'cherry_pick_plate_number: %r' % lcp)
            # Test that the copy-wells have had their volumes deallocated
            if lcp['source_well_id'] != '01000:A01':
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume'] or 0),
                    0, 
                    'lcp vol consumed should be 0, %r'% lcp)
                self.assertEqual(
                    Decimal(lcp['source_copy_well_initial_volume'] or 0),
                    Decimal(lcp['source_copy_well_volume'] or 0),
                    'lcp source_copy_well_initial_volume '
                    '!= source_copy_well_volume: %r'
                    % lcp )
            else:
                self.assertEqual(
                    Decimal(lcp['source_copy_well_volume'] ), 0,
                    'lcp source_copy_well_volume: %r'
                    % lcp )
                # well 01000:A01: previous consumed = 10mL
                expected_vol = Decimal('0.000010000')
                self.assertEqual(
                    Decimal(lcp['source_copy_well_consumed_volume']),
                    expected_vol, 
                    'lcp vol consumed should be %r, %r'
                    % ( expected_vol, lcp) )
                
        # TODO: verify plate mapping:
        # - cpr.keep_source_plate_cherry_picks_together
        # - cpr.random
        # - cpr.wells_to_leave_empty
        
        # 5. Verify updates are disallowed after plating
        # - SCP selection
        # - LCP assignment
        
    def test_4c_update_reservation_keep_plating_zero_deselected(self):
        '''
        Update the plating and mapping of the LCP's, but this time, do not 
        cancel the reservation, instead, update one of the copy selections:
        - deallocate the previous copy-well
        - allocate the new copy-well
        - if API_PARAM_SET_DESELECTED_TO_ZERO==True, then zero out the 
        previously selected copy-well volume
        '''
        
        (cpr_data, current_lcps) = self.test_4_reserve_map_lab_cherry_picks()
        cpr_id = cpr_data['cherry_pick_request_id']
        
        current_lcps = { lcp['source_well_id']:lcp for lcp in current_lcps }
        
        # 20170227 - Modification:
        # Cannot PASTE selection changes to /lab_cherry_pick after plated,
        # unless API_PARAM_OVERRIDE=True
        # - if API_PARAM_SET_DESELECTED_TO_ZERO==True, 
        # set deselected copy-well vols to zero 
        # - wipe out LCP's: still requires cancel plating

        # 2.A Try to change the LCP assignment after plating with /lab_cherry_pick:
        # - Fails because already plated
        # - Fails because override required, OR, reservation must be canceled first
        
        copy4_name = self.library1_copy4['copy_name']
        test_wellid_to_change = '01000:A03'
        test_patch_comment = 'Update lab cherry pick selections after plating'
        header_data = { HEADER_APILOG_COMMENT: test_patch_comment}
        logger.info(
            'Try to Force selection of well %r to %r', test_wellid_to_change,copy4_name)
        lcp_well_1_a3 = current_lcps[test_wellid_to_change]
        copywell_to_deselect = lcp_well_1_a3['source_copywell_id']
        lcp_well_1_a3['source_copy_name'] = copy4_name
        lcp_well_1_a3['source_copy_id'] = None
        lcp_well_1_a3['selected'] = True
        lcp_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', 
            str(cpr_data['cherry_pick_request_id']),
            'lab_cherry_pick_plating'])
        resp = self.api_client.patch(
            lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a3], 
            authentication=self.get_credentials(),
            **header_data)
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        
        self.assertTrue(API_RESULT_ERROR in _data)
        self.assertTrue(API_PARAM_OVERRIDE in _data[API_RESULT_ERROR])
#         self.assertTrue(API_MSG_LCP_PLATES_ASSIGNED in _data[API_RESULT_ERROR])

        # 2.B Try to change the LCP assignment after plating,
        # now using the API_PARAM_OVERRIDE==True
        # part B - if API_PARAM_SET_DESELECTED_TO_ZERO==True, 
        # deselected copy-well vols are set to zero
        
        override_lcp_resource_uri = \
            lcp_resource_uri + '?' + API_PARAM_OVERRIDE + '=true' \
            + '&' + API_PARAM_SET_DESELECTED_TO_ZERO + '=true'
        resp = self.api_client.patch(
            override_lcp_resource_uri, 
            format='json', 
            data=[lcp_well_1_a3], 
            authentication=self.get_credentials(),
            **header_data)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        _meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_COPYWELLS_DEALLOCATED in _meta,
            '%r not in meta: %r' 
            % (API_MSG_COPYWELLS_DEALLOCATED, _meta))
        
        # 2.D Verify copy-well volumes:
        # previous well deallocated
        # new well allocated
        # all other LCP's are unchanged
        new_lab_cherry_picks = self._get_lcps(cpr_id)
        new_lab_cherry_picks = { lcp['source_well_id']:lcp 
            for lcp in new_lab_cherry_picks }
        
        self.assertTrue(test_wellid_to_change in new_lab_cherry_picks,
            '%r not found in %r' %(test_wellid_to_change, new_lab_cherry_picks))
        new_lcp = new_lab_cherry_picks[test_wellid_to_change]
        self.assertEqual(new_lcp['source_copy_name'], copy4_name)
        
        transfer_volume_per_well_approved = \
            Decimal(cpr_data['transfer_volume_per_well_approved'])
        self.assertEqual(
            Decimal(new_lcp['source_copy_well_consumed_volume']),
            transfer_volume_per_well_approved, 
            'lcp vol consumed should be %r, %r'
            % ( transfer_volume_per_well_approved, new_lcp) )
            
        expected_vol = (
            Decimal(lcp_well_1_a3['source_copy_well_initial_volume'])
                -transfer_volume_per_well_approved)
        self.assertEqual(
            expected_vol,
            Decimal(new_lcp['source_copy_well_volume']),
            'lcp well volume reported should be %r, %r'
            % ( expected_vol, new_lcp) )
        
        # Check all other wells are unchanged against "current" (previous) lcps
        for source_well_id, lcp in current_lcps.items():
            self.assertTrue(source_well_id in new_lab_cherry_picks,
                'source well %r not found in new: %r' 
                % (source_well_id, new_lab_cherry_picks.keys()))
            if source_well_id != test_wellid_to_change:
                prev_lcp = current_lcps[source_well_id]
                logger.debug('new_lcp: %r ',lcp)
                logger.debug('prev lcp: %r', prev_lcp)
                for key,val in prev_lcp.items():
                    self.assertEqual(val,lcp[key], 
                        'key: %r, prev: %r, new: %r' % (key, lcp[key],val))
            
        # 2.D1 Check that copywells (from the copywell resource) have the correct information
        cw_resource_uri = '/'.join([
            BASE_URI_DB, 'copywell'])
        copy_wells_to_get = [lcp['source_copywell_id']
            for lcp in new_lab_cherry_picks.values()]
        data_for_get = { 'copywell_id__in': copy_wells_to_get }
        cws = self.get_list_resource(cw_resource_uri, data_for_get) 
        self.assertEqual(len(cws), len(new_lab_cherry_picks))
        for cw in cws:
            self.assertEqual(cw['cherry_pick_screening_count'],1, 
                'cw: {copywell_id}, cpsc: {cherry_pick_screening_count}'.format(**cw))
            self.assertEqual(Decimal(cw['volume']), expected_vol,
                ('cw: {copywell_id}, volume: {volume}'.format(**cw), expected_vol))
        
        # Check the deselected copywell
        
        data_for_get = { 'copywell_id': copywell_to_deselect }
        cws = self.get_list_resource(cw_resource_uri, data_for_get) 
        self.assertEqual(len(cws), 1, 
            'on getting: %r: returned: %r' % (copywell_to_deselect,cws))
        deselected_cw = cws[0]
        self.assertEqual(
            Decimal(deselected_cw['volume']), Decimal(0), 
            'deselected_cw: %r'% deselected_cw)
        
        # 2.D2 Check that the plates have the correct information
        
        current_copies = defaultdict(list)
        for lcp in new_lab_cherry_picks.values():
            current_copies[lcp['source_copy_name']].append(lcp['source_well_id'])
        for copy_name,lcps in current_copies.items():
            logger.info('copy: %r, lcps: %r', copy_name, lcps)
            
        platecopy_resource_uri = '/'.join([BASE_URI_DB, 'librarycopyplate'])
        plates_to_get = set([
            '%s/%s' % (lcp['source_copy_name'],str(lcp['library_plate']))
            for lcp in new_lab_cherry_picks.values()])
        data_for_get = { 'copyplate_id__in': [x for x in plates_to_get] }
        logger.info('get the plates: %r', plates_to_get)
        plates = self.get_list_resource(platecopy_resource_uri, data_for_get) 
        self.assertEqual(len(plates), len(plates_to_get))
        for plate in plates:
            self.assertEqual(plate['cplt_screening_count'],1, 
                'plate: {copyplate_id}, cpsc: {cplt_screening_count}'.format(**plate))
                
        # 2.E Verify copy well logs 
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'cherrypickrequest', 
            'key': cpr_id,
            'diff_keys': 'date_volume_reserved' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('lab_cherry_pick_updates apilogs for cpr: %r', cpr_id)
        selection_update_log = None
        for apilog in apilogs:
            if test_wellid_to_change in apilog['json_field']:
                selection_update_log = apilog
                self.assertTrue(test_patch_comment in apilog['comment'],
                    '%r' % apilog)
        self.assertTrue(selection_update_log is not None)
        logger.info('found selection update log: %r', selection_update_log)
        
        data_for_get={ 
            'limit': 0, 
            'parent_log_id': selection_update_log['id']
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('selection_update_log child logs: %r', apilogs)
        expected_deallocate_logs = 3 # 2 copywell, 1 plate
        self.assertEqual(expected_deallocate_logs, len(apilogs))
        for apilog in apilogs:
            if apilog['ref_resource_name'] == 'copywell':
                self.assertTrue('volume' in apilog['diffs'], 
                    'apilog.diffs: %r'% apilog)
                self.assertTrue(test_patch_comment in apilog['comment'],
                    '%r' % apilog)
                if not copywell_to_deselect in apilog['key']:
                    # NOTE: screening count is not updated for the deselected well
                    self.assertTrue('cherry_pick_screening_count' 
                        in apilog['diffs'], 'apilog.diffs: %r'% apilog)
                    
                self.assertTrue(test_patch_comment in apilog['comment'])
            if apilog['ref_resource_name'] == 'librarycopyplate':
                self.assertTrue(test_patch_comment in apilog['comment'],
                    '%r' % apilog)
                self.assertTrue('cplt_screening_count' in apilog['diffs'], 
                    'apilog.diffs: %r'% apilog)

        
    def test_5_set_plated_screened(self):
        
        (cpr_data, lcps) = self.test_4_reserve_map_lab_cherry_picks()
        
        cpr_id = cpr_data['cherry_pick_request_id']
        cpaps = self._get_cpaps(cpr_id)

        # 1. Patch the plating_date        
        collaborators = self.screen['collaborator_ids']
        collaborators.append(self.screen['lead_screener_id'])
        plates = [cpap['plate_ordinal'] for cpap in cpaps]
        plating_data = {
            'plating_date': '2017-02-10',
            'plated_by_id': collaborators[0],
            'comments': 'test comment for plating...' }
        logger.info('plating data: %r', plating_data)
        
        patch_data = []
        for plate in plates:
            patch_dict = plating_data.copy()
            patch_dict['plate_ordinal'] = plate
            patch_data.append(patch_dict)

        header_data = { HEADER_APILOG_COMMENT: plating_data['comments']}

        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        patch_uri = '/'.join([resource_uri,str(cpr_id),'cherry_pick_plate']) 
        resp = self.api_client.patch(
            patch_uri, format='json', 
            data={ 'objects': patch_data }, 
            authentication=self.get_credentials(), 
            **header_data )
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        #1.a Verify meta from patch
        self.assertTrue(API_RESULT_META in _data)
        self.assertTrue(API_MSG_CPR_PLATES_PLATED in _data[API_RESULT_META])
        self.assertEqual(_data[API_RESULT_META][API_MSG_CPR_PLATES_PLATED],3)
        
        #1.b Retrieve and verify cp assay plates
        cpaps = self._get_cpaps(cpr_id)
        logger.info('cpaps: %r', cpaps)
        
        for cpap in cpaps:
            self.assertEqual(plating_data['plating_date'],cpap['plating_date'])
        
        #1.c Verify plating_date logs
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'cherrypickrequest', 
            'key': cpr_id,
            'diff_keys': 'last_plating_activity_date' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('cpr plate logs: %r', apilogs)
        self.assertTrue(
            len(apilogs)==1, 'wrong number apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        self.assertTrue(plating_data['plating_date'] in apilog['diffs'],
            'wrong diffs: %r' % apilog)
        self.assertEqual(plating_data['comments'], apilog['comment'])

        data_for_get={ 
            'limit': 0, 
            'parent_log_id': apilog['id']
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('cpap plate logs: %r', apilogs)
        self.assertEqual(len(apilogs),len(cpaps), 
            'wrong number cpap apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        self.assertTrue(plating_data['plating_date'] in apilog['diffs'],
            'wrong diffs: %r' % apilog)
        self.assertEqual(plating_data['comments'], apilog['comment'])
        
        #2. Patch the screening date
        screening_data = {
            'screening_date': '2017-02-14',
            'screened_by_id': collaborators[1],
            'comments': 'test comment for screening...' }
        logger.info('screening data: %r', screening_data)
        
        patch_data = []
        for plate in plates:
            patch_dict = screening_data.copy()
            patch_dict['plate_ordinal'] = plate
            patch_data.append(patch_dict)

        header_data = { HEADER_APILOG_COMMENT: screening_data['comments']}

        resource_uri = BASE_URI_DB + '/cherrypickrequest'
        patch_uri = '/'.join([resource_uri,str(cpr_id),'cherry_pick_plate']) 
        resp = self.api_client.patch(
            patch_uri, format='json', 
            data={ 'objects': patch_data }, 
            authentication=self.get_credentials(), 
            **header_data )
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('data: %r', _data)
        #1.a Verify meta from patch
        self.assertTrue(API_RESULT_META in _data)
        meta = _data[API_RESULT_META]
        self.assertTrue(API_MSG_CPR_PLATES_SCREENED in meta)
        expected_cpaps = len(cpaps)
        self.assertEqual(
            meta[API_MSG_CPR_PLATES_SCREENED],expected_cpaps)
        
        #2.b Retrieve and verify cp assay plates
        cpaps = self._get_cpaps(cpr_id)
        logger.info('cpaps: %r', cpaps)
        
        for cpap in cpaps:
            self.assertEqual(
                screening_data['screening_date'],cpap['screening_date'])
        
        #2.c Verify screening_date logs
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'cherrypickrequest', 
            'key': cpr_id,
            'diff_keys': 'last_screening_activity_date' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('cpr plate logs: %r', apilogs)
        self.assertTrue(len(apilogs)==1, 
            'wrong number apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        self.assertTrue(screening_data['screening_date'] in apilog['diffs'],
            'wrong diffs: %r' % apilog)
        self.assertEqual(screening_data['comments'], apilog['comment'])

        data_for_get={ 
            'limit': 0, 
            'parent_log_id': apilog['id']
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('cpap plate logs: %r', apilogs)
        self.assertEqual(len(apilogs),len(cpaps), 
            'wrong number cpap apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        self.assertTrue(screening_data['screening_date'] in apilog['diffs'],
            'wrong diffs: %r' % apilog)
        self.assertEqual(screening_data['comments'], apilog['comment'])
        
        # 3.b Verify delete_lab_cherry_picks disallowed if plating date is set
        delete_lcps_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'delete_lab_cherry_picks'])
        resp = self.api_client.post(
            delete_lcps_resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('delete lcps response: %r', _data)
        self.assertTrue(API_RESULT_ERROR in _data)
        self.assertTrue(API_MSG_NOT_ALLOWED in _data[API_RESULT_ERROR])
        
        # 3.c Verify that cancel_reservation disallowed if plating date is set
        cancel_reservation_resource_uri = '/'.join([
            BASE_URI_DB, 'cherrypickrequest', str(cpr_id), 
            'cancel_reservation'])
        resp = self.api_client.post(
            cancel_reservation_resource_uri,format='json', 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        logger.info('cancel_reservation response: %r', _data)
        self.assertTrue(API_RESULT_ERROR in _data)
        self.assertTrue(API_MSG_NOT_ALLOWED in _data[API_RESULT_ERROR])
    
    def test_b_minimal_set_finder(self):
        
        complete_set = [1,2,3,4,5,6,7]
        
        instance_sets = [
            [2,6],
            [1,2,7],
            [3,4,7],
            [6,7],
            [1,2,3]
        ]
        
        chosen_minimal = lims_utils.find_minimal_satisfying_set(
            complete_set, instance_sets)
        
        logger.info('chosen: %r', chosen_minimal)
        expected = [2,7]
        # Other:
        self.assertEqual(set(expected),set(chosen_minimal))
        
    def test_a_bin_packer(self):
        '''
        Test for bin_packer:
        For the Cherry Pick use case:
        "package" is the collection of all picks from one source plate, 
        "size" is the number of picks for that source plate.
        "bin" is a destination cherry pick assay plate
        "capacity" is the number of wells available on the assay plate 
        (plate size - number of wells to leave empty)
        '''
        capacity = 6
        bins = [2,3,6,7,10]
        package_array = [{'name':x, 'size': x} for x in bins]
        # Note: bins are sorted internally by name
        expected_packed_bins = [
            [{'name': 6, 'size': 6}], 
            [{'name': 2, 'size': 2},{'name': 3, 'size': 3}, 
                {'name': 7, 'size': 1}], 
            [{'name': 10, 'size': 6}], 
            [{'name': 7, 'size': 6}], 
            [{'name': 10, 'size': 4}]]
        
        packed_bins = bin_packer.pack_bins(capacity, package_array)
        logger.info('input bins: %r', bins)
        logger.info('packed bins: %r', packed_bins)
        
        for expected_bin in expected_packed_bins:
            self.assertTrue(expected_bin in packed_bins,
                'expected bin: %r, not found in packed bins: %r'
                %(expected_bin, packed_bins))

        capacity = 8
        expected_packed_bins = [
            [{'name': 7, 'size': 7}], 
            [{'name': 2, 'size': 2},{'name': 6, 'size': 6},], 
            [{'name': 3, 'size': 3}, {'name': 10, 'size': 2}], 
            [{'name': 10, 'size': 8}]]
        
        packed_bins = bin_packer.pack_bins(capacity, package_array)
        logger.info('packed bins: %r', packed_bins)
        for expected_bin in expected_packed_bins:
            self.assertTrue(expected_bin in packed_bins,
                'expected bin: %r, not found in packed bins: %r'
                %(expected_bin, packed_bins))
        
        capacity = 10
        expected_packed_bins = [
            [{'name': 10, 'size': 10}], 
            [{'name': 3, 'size': 3},{'name': 7, 'size': 7},], 
            [{'name': 2, 'size': 2},{'name': 6, 'size': 6},]]
        packed_bins = bin_packer.pack_bins(capacity, package_array)
        logger.info('packed bins: %r', packed_bins)
        for expected_bin in expected_packed_bins:
            self.assertTrue(expected_bin in packed_bins,
                'expected bin: %r, not found in packed bins: %r'
                %(expected_bin, packed_bins))
            
        capacity = 6
        bins = [6,4,2,2]
        package_array = [{'name':x, 'size': x} for x in bins]
        expected_packed_bins = [
            [{'name': 6, 'size': 6}], 
            [{'name': 2, 'size': 2},{'name': 4, 'size': 4}], 
            [{'name': 2, 'size': 2}],
        ]
        
        packed_bins = bin_packer.pack_bins(capacity, package_array)
        logger.info('input bins: %r', bins)
        logger.info('packed bins: %r', packed_bins)
        
        for expected_bin in expected_packed_bins:
            self.assertTrue(expected_bin in packed_bins,
                'expected bin: %r, not found in packed bins: %r'
                %(expected_bin, packed_bins))

        capacity = 6
        bins = [5,4,3,7]
        package_array = [{'name':x, 'size': x} for x in bins]
        expected_packed_bins = [
            [{'name': 5, 'size': 5}, {'name': 7, 'size': 1}], 
            [{'name': 4, 'size': 4}], 
            [{'name': 3, 'size': 3}], 
            [{'name': 7, 'size': 6}]]
        
        packed_bins = bin_packer.pack_bins(capacity, package_array)
        logger.info('input bins: %r', bins)
        
        for expected_bin in expected_packed_bins:
            self.assertTrue(expected_bin in packed_bins,
                'expected bin: %r, not found in packed bins: %r'
                %(expected_bin, packed_bins))
            
        expected_packed_bins = [
            [{'name': 5, 'size': 5}, {'name': 7, 'size': 1}], 
            [{'name': 4, 'size': 4}], 
            [{'name': 3, 'size': 3}], 
            [{'name': 7, 'size': 6}]]
        
        packed_bins = bin_packer.pack_bins(capacity, package_array)
        logger.info('input bins: %r', bins)
        
        for expected_bin in expected_packed_bins:
            self.assertTrue(expected_bin in packed_bins,
                'expected bin: %r, not found in packed bins: %r'
                %(expected_bin, packed_bins))
    
#     def test_A_bin_packer_find_two_bins(self):
#         '''
#         Test for keep_source_plates_together=False:
#         - fit packages from unfilled bins into other unfilled bins;
#         - never split a package over more than two bins 
#         '''
#         
#         capacity = 6
#         package = { 'size': 5, 'name': 5 }
#         already_packed_bins = [
#             [{'name': 2, 'size': 2},{'name': 2, 'size': 2},], # available 2
#             [{'name': 3, 'size': 3},{'name': 2, 'size': 2},], # available 1
#             [{'name': 1, 'size': 1},{'name': 2, 'size': 2},], # available 3
#         ]
#         expected_available = [
#             [{'name': 2, 'size': 2},{'name': 2, 'size': 2},], # available 2
#             [{'name': 1, 'size': 1},{'name': 2, 'size': 2},], # available 3
#         ]
#         available_bins = bin_packer.find_two_bins_for_package(
#             capacity, package, already_packed_bins)
#         logger.info('found available_bins: %r', available_bins)
#         self.assertTrue(len(available_bins)>0)
#         for available_bin in available_bins:
#             self.assertTrue(available_bin in expected_available,
#                 'expected bin: %r, not found in packed bins: %r'
#                 %(available_bin, expected_available))
# 
#         package = { 'size': 5, 'name': 5 }
#         already_packed_bins = [
#             [{'name': 2, 'size': 2},{'name': 4, 'size': 4},], # available 0
#             [{'name': 3, 'size': 3},{'name': 2, 'size': 2},], # available 1
#             [{'name': 2, 'size': 2},{'name': 2, 'size': 2},], # available 2 **
#             [{'name': 1, 'size': 1},{'name': 2, 'size': 2},], # available 3 **
#             [{'name': 2, 'size': 2},], # available 4
#         ]
#         # ** don't pick these because only 1 space needed after picking 4
#         expected_available = [
#             [{'name': 2, 'size': 2},], # available 4
#             [{'name': 3, 'size': 3},{'name': 2, 'size': 2},], # available 1
#         ]
#         available_bins = bin_packer.find_two_bins_for_package(
#             capacity, package, already_packed_bins)
#         logger.info('found available_bins: %r', available_bins)
#         self.assertTrue(len(available_bins)>0)
#         for available_bin in available_bins:
#             self.assertTrue(available_bin in expected_available,
#                 'expected bin: %r, not found in packed bins: %r'
#                 %(available_bin, expected_available))
# 
#         package = { 'size': 5, 'name': 5 }
#         already_packed_bins = [
#             [{'name': 2, 'size': 2},{'name': 4, 'size': 4},], # available 0
#             [{'name': 3, 'size': 3},{'name': 2, 'size': 2},], # available 1
#             [{'name': 2, 'size': 2},{'name': 2, 'size': 2},], # available 2
#             [{'name': 2, 'size': 2},{'name': 2, 'size': 2},], # available 2
#         ]
#         available_bins = bin_packer.find_two_bins_for_package(
#             capacity, package, already_packed_bins)
#         logger.info('found available_bins: %r', available_bins)
#         # Nothing picked because it would require splitting over 3
#         self.assertTrue(len(available_bins)==0)



class MutualScreensTest(DBResourceTestCase,ResourceTestCase):
    
    def test_mutual_positives_to_screen(self):
        pass
        # create two screens
        # create two screen results
        # set data sharing levels
        # find assay wells that overlap
        # verify that members of screen1 can see mutual positives for screen2
        
class ScreensaverUserResource(DBResourceTestCase):
        
    def setUp(self):
        super(ScreensaverUserResource, self).setUp()

    def tearDown(self):
        DBResourceTestCase.tearDown(self)
        
        logger.info('delete resources')
        UserChecklist.objects.all().delete()
        AttachedFile.objects.all().delete()
        ServiceActivity.objects.all().delete()
        logger.info('delete users, including: %r', self.username)
        ScreensaverUser.objects.all().exclude(username=self.username).delete()
         
        UserGroup.objects.all().delete()
        UserProfile.objects.all().exclude(username=self.username).delete()
        User.objects.all().exclude(username=self.username).delete()
        ApiLog.objects.all().delete()
        
        # removed: should not be nec
        # Vocabulary.objects.all().filter(scope__contains='labaffiliation.').delete()
        LabAffiliation.objects.all().delete()
        
    def test01_create_user_iccbl(self):

        logger.info('test01_create_user_iccbl...')
        # create users using only ecommons, username
        simple_user_input = { 
            'ecommons_id': 'tester01c',
            'first_name': 'FirstName01c',
            'last_name': 'LastName01c',    
            'email': 'tester01c@limstest.com',    
            'harvard_id': '332122',
            'harvard_id_expiration_date': '2018-05-01',
        }
        resource_uri = BASE_URI_DB + '/screensaveruser'
        resource_test_uri = '/'.join([
            resource_uri,simple_user_input['ecommons_id']])
        created_user = self._create_resource(
            simple_user_input, resource_uri, resource_test_uri)
        self.assertTrue(
            simple_user_input['ecommons_id']==created_user['username'],
            'username should equal the ecommons id if only ecommons is'
            ' provided: %r, %r' % (simple_user_input,created_user))
        
        # Verify that the ecommons & username cannot be changed
        
        user_update = {'ecommons_id': 'testerxxxx'}
        resource_uri = '/'.join([
            BASE_URI_DB,'screensaveruser',simple_user_input['ecommons_id']])
        resp = self.api_client.patch(
            resource_uri, 
            format='json', data=user_update, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        resp_data = self.deserialize(resp)
        logger.info('(expected) error response: %r', resp_data)
        
        self.assertTrue(API_RESULT_ERROR in resp_data)
        self.assertTrue('username' in resp_data[API_RESULT_ERROR])
        
    def test0_create_user(self):
        
        logger.info('test0_create_user...')
        self.user1 = self.create_screensaveruser({ 
            'username': 'st1',
            'is_active': False
        })

        logger.info('test0_create_user 2...')
        self.user2 = self.create_screensaveruser({ 
            'username': 'st2',
            'is_active': False
            })
        logger.info('test0_create_user screening user...')
        self.screening_user = self.create_screensaveruser(
            { 'username': 'screening1', 'is_active': True })
        
        self.assertTrue(self.screening_user['is_active'])
        
        # FIXME: test more specific admin user permissions
        logger.info('test0_create_user admin user...')
        self.test_admin_user = self.create_screensaveruser(
            { 'username': 'adminuser'})
        logger.info('test0_create_user admin user 2...')
        self.test_admin_user2 = self.create_screensaveruser(
            { 'username': 'adminuser2'})
        # create an admin
        patch_obj = { 'objects': [
            {
                'username': 'adminuser',
                'is_superuser': True,
                'is_active': True
            },
            {
                'username': 'adminuser2',
                'is_superuser': True,
                'is_active': True
            },
        ]}
        resource_uri = BASE_URI_DB + '/screensaveruser'

        logger.info('test0_create_user patch admin users...')
        resp = self.api_client.patch(
            resource_uri, 
            format='json', data=patch_obj, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))

    def test02_create_user_without_username(self):
    
        # 1. create a user with no username (first, last must be unique)
        logger.info('test02_create_user_without_username...')
        # create a non-login user with no username
        user1_input_data = { 
            'first_name': 'FirstNameUniq1',
            'last_name': 'LastNameUniq1',    
        }
        resource_uri = BASE_URI_DB + '/screensaveruser'
        
        _data_for_get = { 
            'limit': 0,
            'includes': '*',
            'HTTP_ACCEPT': 'application/json'
        }
        resp = self.api_client.post(
            resource_uri, format='json', data=user1_input_data, 
            authentication=self.get_credentials(), **_data_for_get)
        logger.info('resp: %r', resp)
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in new_obj)
        self.assertEqual(len(new_obj[API_RESULT_DATA]),1,
            'more than one object returned for: %r, returns: %r'
            % (resource_uri,new_obj))
        user1_output_data = new_obj[API_RESULT_DATA][0]
        logger.info('post create user (no username): %r', user1_output_data)
        # 1.a verify that a user with no username has been created
        self.assertIsNone(user1_output_data.get('username'))
        self.assertFalse(user1_output_data.get('is_active'))
        
        # 1.b verify that a second attempt fails with the same first/last name
        resp = self.api_client.post(
            resource_uri, format='json', data=user1_input_data, 
            authentication=self.get_credentials(), **_data_for_get)
        
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        self.assertTrue(API_RESULT_ERROR in new_obj)
        # NOTE: should be an error for a non-unique user
        errors = new_obj[API_RESULT_ERROR]
        self.assertTrue('first_name' in errors)
        self.assertTrue('last_name' in errors)
        logger.info('second attempt (expected) errors reported: %r', errors)
        
        # 1.c create another
        user2_input_data = { 
            'first_name': 'FirstNameUniq2',
            'last_name': 'LastNameUniq2',    
        }
        resp = self.api_client.post(
            resource_uri, format='json', data=user2_input_data, 
            authentication=self.get_credentials(), **_data_for_get)
        
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in new_obj)
        self.assertEqual(len(new_obj[API_RESULT_DATA]),1,
            'more than one object returned for: %r, returns: %r'
            % (resource_uri,new_obj))
        user2_output_data = new_obj[API_RESULT_DATA][0]
        logger.info('post create user (no username): %r', user2_output_data)

        # 2. set the user as a lab_head for users and screens
        # 2.a convert the user to a lab head
        
        user1_data_as_lab_head = self.update_user_to_lab_head(user1_output_data)
        self.assertEqual(
            user1_data_as_lab_head['screensaver_user_id'], 
            user1_output_data['screensaver_user_id'])
        self.assertEqual(
            user1_data_as_lab_head['classification'], 'principal_investigator')
        
        user2_data_patch2 = {
            'screensaver_user_id': user2_output_data['screensaver_user_id'],
            'lab_head_id': user1_data_as_lab_head['screensaver_user_id'] }
        
        resource_uri = BASE_URI_DB + '/screensaveruser/'
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=user2_data_patch2, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        user2_output_data2 = _data[API_RESULT_DATA][0]
        
        self.assertTrue(
            user2_output_data2['lab_head_id'], 
            user1_data_as_lab_head['screensaver_user_id'])
        
        # 3. Set the username for user1
        # - this creates a reports_userprofile for the user
        user1_input_data2 = {
            'username': 'usr1',
        }
        resource_uri = '/'.join([
            BASE_URI_DB,'screensaveruser', str(user1_output_data['screensaver_user_id'])])
        logger.info('patch in the username for the user: %r, %r', 
            resource_uri, user1_input_data2)
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=user1_input_data2, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        user1_output_data2 = _data[API_RESULT_DATA][0]
        
        # 3.a Verify username is set
        self.assertEqual(
            user1_input_data2['username'], user1_output_data2['username'])

    def test03_update_lab_affiliation_name(self): 
        ''' 
        Simple test for the updating the name of a Lab Affiliation
        '''
        
        lab_affiliation = self.create_lab_affiliation()
        
        lab_affiliation['name'] = 'Test new Lab Affiliation Name'
        resource_uri = BASE_URI_DB + '/labaffiliation/'        
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=lab_affiliation, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        new_lab_affiliation =_data[API_RESULT_DATA][0]
        logger.info('created lab: %r', new_lab_affiliation)
        
        for key,val in lab_affiliation.items():
            self.assertEqual(lab_affiliation[key],new_lab_affiliation[key])
        
    def test01_create_lab_head(self):

        logger.info('test01_create_lab_head...')
        
        logger.info('1. Create the lab head...')
        lab_head = self.create_lab_head()
        logger.info('lab_head created: %r', lab_head)
        self.assertTrue(
            'lab_affiliation_id' in lab_head, 
            'Lab head does not contain "lab_affiliation_id": %r' % lab_head)
        self.assertTrue(
            'lab_head_id' in lab_head, 
            'Lab head does not contain "lab_head_id": %r' % lab_head)
        self.assertEqual(lab_head['screensaver_user_id'], lab_head['lab_head_id'])
        self.assertEqual(lab_head['classification'], 'principal_investigator')
        
        logger.info('2. Assign a user to the Lab Head...')
        
        user = self.create_screensaveruser({ 'username': 'st1'})
        user_data = {
            'username': user['username'],
            'lab_head_id': lab_head['screensaver_user_id']
            }
        resource_uri = '/'.join([BASE_URI_DB, 'screensaveruser'])
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=user_data, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        updated_user = _data[API_RESULT_DATA][0]

        logger.info('User: %r (with lab head set)', updated_user)
        
        self.assertEqual(
            updated_user['lab_head_id'], lab_head['screensaver_user_id'])
        self.assertEqual(updated_user['lab_name'], lab_head['lab_name'])
        self.assertEqual(
            updated_user['lab_affiliation_name'], 
            lab_head['lab_affiliation_name'])
        self.assertEqual(
            updated_user['lab_affiliation_category'], 
            lab_head['lab_affiliation_category'])
        
        logger.info('3. Change the Lab Head classification to "unassigned"...') 
        # - removes the lab members
        
        lab_head_update = {
            'username': lab_head['username'],
            'classification': 'unassigned'
        }
        logger.info('lab_head update: %r', lab_head_update)
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=lab_head_update, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        updated_lab_head = _data[API_RESULT_DATA][0]

        logger.info(
            'updated user; no longer lab_head: %r (classification unset)', 
            updated_lab_head)
        
        self.assertEqual(
            updated_lab_head['lab_affiliation_id'],None, 
            'updated_lab_head should not contain "lab_affiliation_id": %r' 
                % updated_lab_head)
        self.assertEqual(
            updated_lab_head['lab_head_id'],None, 
            'updated_lab_head should not contain "lab_head_id": %r' 
                % updated_lab_head)
        self.assertEqual(updated_lab_head['lab_name'], None)
        self.assertEqual(
            updated_lab_head['lab_affiliation_name'], None)
        self.assertEqual(
            updated_lab_head['lab_affiliation_category'], None)
        self.assertEqual(updated_lab_head['classification'], 'unassigned')
        
        logger.info('3.A check that the user no longer has the lab head...')
        
        updated_user = self.get_single_resource(
            resource_uri, {'username': user['username']})
        
        self.assertEqual(updated_user['lab_head_id'], None)
        self.assertEqual(updated_user['lab_name'], None)
        self.assertEqual(
            updated_user['lab_affiliation_name'], None)
        self.assertEqual(
            updated_user['lab_affiliation_category'], None)
        
        logger.info('4. Update the lab members directly on the lab head...')
        
        logger.info('4a. reset the user to a lab_head...')
        lab_head_update2 = {
            'username': lab_head['username'],
            'classification': 'principal_investigator',
            'lab_affiliation_id': lab_head['lab_affiliation_id']
        }
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=lab_head_update2, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        updated_lab_head2 = _data[API_RESULT_DATA][0]
        logger.info(
            'updated_lab_head: %r (classification reset to pi)', 
            updated_lab_head2)
        self.assertEqual(
            updated_lab_head2['classification'], 'principal_investigator')
        
        logger.info('4.b set lab member...')
        
        lab_head_update3 = {
            'username': lab_head['username'],
            'lab_member_ids': [str(user['screensaver_user_id']),]
        }
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=lab_head_update3, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        _data = self.deserialize(resp)
        self.assertTrue(API_RESULT_DATA in _data)
        self.assertEqual(len(_data[API_RESULT_DATA]), 1)        
        updated_lab_head3 = _data[API_RESULT_DATA][0]
        logger.info('updated_lab_head3: %r (add member)', updated_lab_head3)
        self.assertEqual(
            updated_lab_head3['lab_member_ids'], 
            lab_head_update3['lab_member_ids'])
        
        user_after_update3 = self.get_single_resource(
            resource_uri, {'username': user['username']})
        logger.info('User: %r (lab head set)', user_after_update3)
        
        self.assertEqual(
            user_after_update3['lab_head_id'], lab_head['screensaver_user_id'])
        self.assertEqual(user_after_update3['lab_name'], lab_head['lab_name'])
        self.assertEqual(
            user_after_update3['lab_affiliation_name'], 
            lab_head['lab_affiliation_name'])
        self.assertEqual(
            user_after_update3['lab_affiliation_category'], 
            lab_head['lab_affiliation_category'])
        
    def test01a_update_lab_head_dsl(self):
        
        # Verify that the lab member dsl's are updated on lab head update.
        
        
        
        
        pass
        
        
    def test1_patch_usergroups(self):
        ''' 
        Verify the reports.UserResource usergroup functionality from the 
        ScreensaverUserResource
        '''
        
        logger.info('test1_patch_usergroups...')
        self.test0_create_user();
        group_patch = { 'objects': [
            { 
                'name': 'usergroup1'
            },
            { 
                'name': 'usergroup2'
            },
            { 
                'name': 'usergroup3'
            }
        ]};
        
        try:       
            resource_uri = BASE_REPORTS_URI + '/usergroup/'
            resp = self.api_client.put(resource_uri, 
                format='json', data=group_patch, 
                authentication=self.get_credentials())
            self.assertTrue(
                resp.status_code in [200,201,202], 
                (resp.status_code, self.get_content(resp)))

            resp = self.api_client.get(resource_uri, format='json', 
                authentication=self.get_credentials(), data={ 'limit': 999 })
            self.assertTrue(
                resp.status_code in [200,201,202], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            self.assertEqual(
                len(new_obj[API_RESULT_DATA]), 
                len(group_patch[API_RESULT_DATA]), new_obj)
            
            for i,item in enumerate(group_patch[API_RESULT_DATA]):
                result, obj = find_obj_in_list(item, new_obj[API_RESULT_DATA])
                self.assertTrue(
                    result, 
                    ('bootstrap item not found', item, 
                        new_obj[API_RESULT_DATA]))
                logger.info('item found: %r', obj)        
        except Exception, e:
            logger.exception('on group_patch: %r', group_patch)
            raise

        userpatch = [   
            {
                'username': self.user1['username'],
                'usergroups': ['usergroup1',]
            },
            {
                'username': self.screening_user['username'],
                'usergroups': ['usergroup2',]
            },
            {
                'username': self.test_admin_user['username'],
                'usergroups': ['usergroup3',]
            },
        ]
        resource_uri = BASE_URI_DB + '/screensaveruser'
        resp = self.api_client.patch(resource_uri, 
            format='json', data={ 'objects': userpatch }, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))

        for patchobj in userpatch:
            logger.info('test patchobj: %r', patchobj)
            groupname = patchobj['usergroups'][0]
            username = patchobj['username']
            data_for_get = { 
                'limit': 0, 
                'includes': ['*'],
                'usergroups__eq': groupname
            }
            resp = self.api_client.get(resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get )
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            result = self.deserialize(resp)
            new_objs = result[API_RESULT_DATA]
            self.assertEqual(
                len(new_objs), 1,
                'wrong number of users found for group: %r, %r: %r'
                % (len(new_objs), groupname, [x['username'] for x in new_objs] ))
            self.assertEqual(username, new_objs[0]['username'])

    def test2_user_checklist(self):
        
        logger.info('test2_user_checklist...')
        self.test0_create_user();
        
        # Note "get_credentials" returns the superuser
        admin_performing_operation = self.username
        test_su_id = str(self.user1['screensaver_user_id'])
        
        # FIXME: create a "ChecklistAdmin" usergroup
        # Note admin assigned to the checklist operation
        checklist_admin = self.test_admin_user
        checklist_admin2 = self.test_admin_user2

        # 1. Create a UserChecklist item        
        checklist_patch = {
            'admin_username': checklist_admin['username'], 
            'name': "added_to_iccb_l_users_email_list",
            'is_checked': True,
            'date_effective': "2015-09-02",
            'screensaver_user_id': test_su_id
        }
        
        test_comment = 'Some test comment 123 xyz'
        
        header_data = { HEADER_APILOG_COMMENT: test_comment}
        header_data['HTTP_ACCEPT'] = JSON_MIMETYPE
        
        patch_uri = '/'.join([BASE_URI_DB,'userchecklist',test_su_id])
        resp = self.api_client.patch(
            patch_uri, 
            format = 'json', 
            data = { API_RESULT_DATA: checklist_patch}, 
            authentication=self.get_credentials(),
            **header_data)
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))

        data_for_get = { 'limit': 0, 'includes': ['*'] }
        checklist_uri = '/'.join([
            BASE_URI_DB,'userchecklist',test_su_id, checklist_patch['name']])
        resp = self.api_client.get(
            checklist_uri,
            format='json', 
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('UserChecklist created: %r', new_obj)
        result,msgs = assert_obj1_to_obj2(
            checklist_patch, new_obj)
        self.assertTrue(result,msgs)
        self.assertEqual(new_obj['status'],'activated')
        
        # 1.A checklistitem logs
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'userchecklist', 
            'key__contains': test_su_id + '/' + new_obj['name']
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        self.assertEqual(
            len(apilogs),1, 'wrong number of apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('log: %r', apilog)
        self.assertTrue(apilog['comment']==test_comment,
            'comment %r should be: %r' % (apilog['comment'], test_comment))
        self.assertTrue('status' in apilog['diff_keys'])
        diffs = json.loads(apilog['diffs'])
        self.assertTrue('status' in diffs)
        self.assertTrue('activated' in diffs['status'])
        self.assertEqual(admin_performing_operation, apilog['username'],
            'wrong admin username recorded: %r'% apilog)
        
        # 2. Deactivate
        checklist_patch2 = {
            'admin_username': checklist_admin2['username'], 
            'name': "added_to_iccb_l_users_email_list",
            'is_checked': False,
            'date_effective': "2016-09-02",
            'screensaver_user_id': test_su_id
        }
        
        test_comment = 'Some test comment 123 xyz'
        
        header_data = { HEADER_APILOG_COMMENT: test_comment}
        header_data['HTTP_ACCEPT'] = JSON_MIMETYPE
        
        resp = self.api_client.patch(
            patch_uri, 
            format = 'json', 
            data = { API_RESULT_DATA: checklist_patch2}, 
            authentication=self.get_credentials(),
            **header_data)
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))

        data_for_get = { 'limit': 0, 'includes': ['*'] }
        resp = self.api_client.get(
            checklist_uri,
            format='json', 
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('UserChecklist created: %r', new_obj)
        result,msgs = assert_obj1_to_obj2(
            checklist_patch2, new_obj)
        self.assertTrue(result,msgs)
        self.assertEqual(new_obj['status'],'deactivated')
        
        # 2.A check for deactivate log
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'userchecklist', 
            'key__contains': test_su_id + '/' + new_obj['name']
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        self.assertTrue(
            len(apilogs) == 2, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[1]
        logger.info('log: %r', apilog)
        self.assertTrue(apilog['comment']==test_comment,
            'comment %r should be: %r' % (apilog['comment'], test_comment))
        self.assertTrue('status' in apilog['diff_keys'])
        diffs = json.loads(apilog['diffs'])
        self.assertTrue('status' in diffs)
        self.assertTrue('deactivated' in diffs['status'])
        self.assertEqual(admin_performing_operation, apilog['username'],
            'wrong admin username recorded: %r'% apilog)
        
        # 3. modify date only
        
    def test3_attached_files(self):
        
        logger.info('test3_attached_files...')
        self.test0_create_user();

        # Test using embedded "contents" field               
        test_su_id = str(self.user1['screensaver_user_id'])
        # FIXME: create an admin with ScreensaverUser/write permission
        attached_file_admin = self.test_admin_user
        
        attachedfile_item_post = {
            'created_by_username': attached_file_admin['username'], 
            'type': '2009_iccb_l_nsrb_small_molecule_user_agreement', 
            'filename': "test_pasted_text.txt",
            'contents': "This is a test of pasted text\n1\n2\n3\n\n end\n",
            'file_date': '2015-10-10'
            }

        content_type = MULTIPART_CONTENT
        resource_uri = \
            BASE_URI_DB + '/screensaveruser/%s/attachedfiles/' % test_su_id
        
        authentication=self.get_credentials()
        kwargs = {}
        kwargs['HTTP_AUTHORIZATION'] = authentication
        kwargs['HTTP_ACCEPT'] = JSON_MIMETYPE
        
        logger.info('Post attached file item: %r', attachedfile_item_post)
        # NOTE: content_type arg is req'd with django.test.Client.post
        # NOTE: content_type defaults to MULTIPART_CONTENT
        resp = self.django_client.post(
            resource_uri, content_type=content_type, 
            data=attachedfile_item_post, **kwargs)
        self.assertTrue(
            resp.status_code in [201, 202], 
            (resp.status_code,self.get_content(resp)))
        
        data_for_get = { 'limit': 0, 'includes': ['*'] }
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('new obj: %s ' % new_obj)
        result,msgs = assert_obj1_to_obj2(
            attachedfile_item_post, new_obj[API_RESULT_DATA][0], 
            excludes=['contents'])
        self.assertTrue(result,msgs)
        af = new_obj[API_RESULT_DATA][0]
        uri = '/db/attachedfile/%s/content' % af['attached_file_id']
        try:
            admin_user = User.objects.get(username=attached_file_admin['username'])
            view, args, kwargs = resolve(uri)
            kwargs['request'] = self.api_client.client.request()
            kwargs['request'].user=admin_user
            result = view(*args, **kwargs)
            logger.info('attached_file request view result: %r',result)
            self.assertEqual(
                result.content, attachedfile_item_post['contents'], 
                'download file view returns wrong contents: %r' 
                    % result.content)
        except Exception, e:
            logger.info('no file found at: %r', uri)
            raise
    
    def test3a_attached_file_filesystem(self):
        
        logger.info('test3a_attached_file_filesystem...')
        
        self.test0_create_user();

        # Test using embedded "contents" field               
        test_su_id = str(self.user1['screensaver_user_id'])
        admin_username = self.test_admin_user['username']
        attachedfile_item_post = {
            'created_by_username': admin_username, 
            'type': '2009_iccb_l_nsrb_small_molecule_user_agreement', 
        }

        content_type = MULTIPART_CONTENT
        resource_uri = \
            BASE_URI_DB + '/screensaveruser/%s/attachedfiles/' % test_su_id
        authentication=self.get_credentials()
        kwargs = {}
        kwargs['HTTP_AUTHORIZATION'] = authentication
        kwargs['HTTP_ACCEPT'] = JSON_MIMETYPE
        file = 'iccbl_sm_user_agreement_march2015.pdf'
        filename = \
            '%s/db/static/test_data/useragreement/%s' %(APP_ROOT_DIR,file)
        logger.info('Open and POST file: %r', filename)
        with open(filename) as input_file:

            logger.info('POST with attached_file to the server')
            attachedfile_item_post['attached_file'] = input_file

            logger.info('Post attached file %r', filename)
            # NOTE: content_type arg is req'd with django.test.Client.post
            # NOTE: content_type defaults to MULTIPART_CONTENT
            resp = self.django_client.post(
                resource_uri, content_type=content_type, 
                data=attachedfile_item_post, **kwargs)
            self.assertTrue(
                resp.status_code in [201, 202], 
                (resp.status_code,self.get_content(resp)))
        
        data_for_get = { 'limit': 0, 'includes': ['*'] }
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('new obj: %s ' % new_obj)
        result,msgs = assert_obj1_to_obj2(
            attachedfile_item_post, new_obj[API_RESULT_DATA][0], 
            excludes=['attached_file'])
        self.assertTrue(result,msgs)
        af = new_obj[API_RESULT_DATA][0]
        uri = '/db/attachedfile/%s/content' % af['attached_file_id']
        try:
            admin_user = User.objects.get(username=admin_username)
            view, args, kwargs = resolve(uri)
            kwargs['request'] = self.api_client.client.request()
            kwargs['request'].user=admin_user
            result = view(*args, **kwargs)
            logger.info('attached_file request view result: %r',result)
            output_filename = '%s.out.%s' % tuple(filename.split('.'))
            logger.info('write %s to %r', filename, output_filename)
            with open(output_filename, 'w') as out_file:
                out_file.write(self.get_content(result))
            self.assertTrue(filecmp.cmp(filename,output_filename), 
                'input file: %r, not equal to output file: %r' 
                % (filename, output_filename))
            os.remove(output_filename)    
        except Exception, e:
            logger.info('no file found at: %r', uri)
            raise
        
        # TODO: delete attached file
        # TODO: attachedfile logs
    
    def test4_user_agreement_updator(self):
        
        logger.info('test4_user_agreement_updator...')
        self.test0_create_user();
        test_su_id = str(self.user1['screensaver_user_id'])
        # FIXME: create an admin with ScreensaverUser/write priv
        admin_username = self.test_admin_user['username']
        
        # 1. small molecule data sharing
        
        useragreement_item_post = {
            'admin_user': admin_username,
            'created_by_username': admin_username, 
            'type': 'sm',
            'data_sharing_level': 2 
            }
        test_comment = 'test update comment for user agreement'
        content_type = MULTIPART_CONTENT
        resource_uri = \
            BASE_URI_DB + '/screensaveruser/%s/useragreement/' % test_su_id
        
        authentication=self.get_credentials()
        kwargs = { 'limit': 0, 'includes': ['*'] }
        kwargs['HTTP_AUTHORIZATION'] = authentication
        kwargs[HEADER_APILOG_COMMENT] = test_comment
        kwargs['HTTP_ACCEPT'] = JSON_MIMETYPE
        
        file = 'iccbl_sm_user_agreement_march2015.pdf'
        filename = \
            '%s/db/static/test_data/useragreement/%s' %(APP_ROOT_DIR,file)
        logger.info('Open and POST file: %r', filename)
        with open(filename) as input_file:

            logger.info('POST user agreement to the server...')
            useragreement_item_post['attached_file'] = input_file
            
            resp = self.django_client.post(
                resource_uri, content_type=MULTIPART_CONTENT, 
                data=useragreement_item_post, **kwargs)
            if resp.status_code not in [200]:
                logger.info(
                    'resp code: %d, resp: %r, content: %r', 
                    resp.status_code, resp, resp.content)
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code))
        
        # Tests: 
        # 1.a check that the user agreement is an attached file to the user
        attached_type = db.api.UserAgreementResource.FILE_TYPE_SM_UA
        data_for_get = { 
            'limit': 0, 'includes': ['*'],
            'type__eq': attached_type
        }
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('new obj: %s ' % new_obj)
        self.assertTrue(API_RESULT_DATA in new_obj)
        af = new_obj[API_RESULT_DATA][0]
        uri = '/db/attachedfile/%s/content' % af['attached_file_id']
        try:
            admin_user = User.objects.get(username=admin_username)
            view, args, kwargs = resolve(uri)
            kwargs['request'] = self.api_client.client.request()
            kwargs['request'].user=admin_user
            result = view(*args, **kwargs)
            output_filename = '%s.out.%s' % tuple(filename.split('.'))
            logger.info('write %s to %r', filename, output_filename)
            with open(output_filename, 'w') as out_file:
                out_file.write(self.get_content(result))
            self.assertTrue(filecmp.cmp(filename,output_filename), 
                'input file: %r, not equal to output file: %r' 
                % (filename, output_filename))    
        except Exception, e:
            logger.exception('no file found at: %r', uri)
            raise
        
        # 1.b check that the data sharing level for the type is assigned
        
        resource_uri = BASE_URI_DB + '/screensaveruser'
        resource_uri = '/'.join([resource_uri,self.user1['username']])
        user_data = self.get_single_resource(resource_uri)
        
        self.assertEqual(
            user_data['sm_data_sharing_level'],
            useragreement_item_post['data_sharing_level'])
        self.assertEqual(
            user_data['is_active'], True)

        # 1.c check that a checklist item has been created for the user agreement
        
        resource_uri = BASE_URI_DB + '/userchecklist'
        resource_uri = '/'.join([resource_uri,test_su_id])
        checklist_items = self.get_list_resource(
            resource_uri, {'status__eq': 'activated'})
        logger.info('checklist_items: %r', checklist_items)
        self.assertTrue(len(checklist_items)==1)
        val = db.api.UserAgreementResource.CHK_TYPE_SM_UA_ACTIVE
        self.assertTrue(checklist_items[0]['name'] == val,
            'wrong checklist item - expected name: %r, %r'
            %(val, checklist_items[0]))
        
        # 1.d check logs
        
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'screensaveruser', 
            'key': test_su_id,
            'diff_keys__contains': 'sm_data_sharing_level' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('logs: %r', apilogs)
        self.assertEqual(
            len(apilogs),1, 'wrong apilog count: %r' % apilogs)
        apilog = apilogs[0]
        self.assertTrue(apilog['comment']==test_comment,
            'comment %r should be: %r' % (apilog['comment'], test_comment))
        self.assertTrue('sm_data_sharing_level' in apilog['diff_keys'])

        # 2 rnai data sharing
        
        useragreement_item_post = {
            'admin_user': admin_username,
            'created_by_username': admin_username, 
            'type': 'rnai',
            'data_sharing_level': 1 
            }
        test_comment = 'test update rna comment for user agreement'
        content_type = MULTIPART_CONTENT
        resource_uri = \
            BASE_URI_DB + '/screensaveruser/%s/useragreement/' % test_su_id
        
        authentication=self.get_credentials()
        kwargs = { 'limit': 0, 'includes': ['*'] }
        kwargs['HTTP_AUTHORIZATION'] = authentication
        kwargs[HEADER_APILOG_COMMENT] = test_comment
        kwargs['HTTP_ACCEPT'] = JSON_MIMETYPE
        
        file = 'iccbl_rnai_ua_march2015.pdf'
        filename = \
            '%s/db/static/test_data/useragreement/%s' %(APP_ROOT_DIR,file)
        logger.info('Open and POST file: %r', filename)
        with open(filename) as input_file:

            logger.info('PUT user agreement to the server...')
            useragreement_item_post['attached_file'] = input_file
            
            resp = self.django_client.post(
                resource_uri, content_type=MULTIPART_CONTENT, 
                data=useragreement_item_post, **kwargs)
            if resp.status_code not in [200]:
                logger.info(
                    'resp code: %d, resp: %r, content: %r', 
                    resp.status_code, resp, resp.content)
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code))
        
        # 2.a check that the user agreement is an attached file to the user
        attached_type = db.api.UserAgreementResource.FILE_TYPE_RNA_UA
        data_for_get = { 
            'limit': 0, 'includes': ['*'], 
            'type__eq': attached_type 
        }
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        self.assertEqual(
            len(new_obj[API_RESULT_DATA]), 1, 
            'wrong count of user agreements returned: %r' % new_obj)
        
        af = new_obj[API_RESULT_DATA][0]
        uri = '/db/attachedfile/%s/content' % af['attached_file_id']
        try:
            admin_user = User.objects.get(username=admin_username)
            view, args, kwargs = resolve(uri)
            kwargs['request'] = self.api_client.client.request()
            kwargs['request'].user=admin_user
            result = view(*args, **kwargs)
            output_filename = '%s.out.%s' % tuple(filename.split('.'))
            logger.info('write %s to %r', filename, output_filename)
            with open(output_filename, 'w') as out_file:
                out_file.write(self.get_content(result))
            self.assertTrue(filecmp.cmp(filename,output_filename), 
                'input file: %r, not equal to output file: %r' 
                % (filename, output_filename))    
        except Exception, e:
            logger.exception('no file found at: %r', uri)
            raise

        # 2.b check that the data sharing level for the type is assigned
        
        resource_uri = BASE_URI_DB + '/screensaveruser'
        resource_uri = '/'.join([resource_uri,test_su_id])
        user_data = self.get_single_resource(resource_uri)
        
        self.assertEqual(
            user_data['rnai_data_sharing_level'],
            useragreement_item_post['data_sharing_level'])
        self.assertEqual(
            user_data['is_active'], True)

        # 2.c check that checklist item has been created for user agreement
        
        resource_uri = BASE_URI_DB + '/userchecklist'
        resource_uri = '/'.join([resource_uri,test_su_id])
        checklist_items = self.get_list_resource(resource_uri, 
            {'status__eq': 'activated',
             'name__eq': u'current_rnai_user_agreement_active'})
        self.assertEqual(len(checklist_items),1)
        val = db.api.UserAgreementResource.CHK_TYPE_RNA_UA_ACTIVE
        self.assertTrue(checklist_items[0]['name'] == val,
            'wrong checklist item - expected name: %r, %r'
            %(val, checklist_items[0]))
        
        # 2.d check logs
        
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'screensaveruser', 
            'key': test_su_id,
            'diff_keys__contains': 'rnai_data_sharing_level',
            'order_by': 'date_created' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        logger.info('logs: %r', apilogs)
        self.assertEqual(
            len(apilogs),1, 'wrong apilogs count: %r' % apilogs)
        apilog = apilogs[0]
        self.assertTrue(apilog['comment']==test_comment,
            'comment %r should be: %r' % (apilog['comment'], test_comment))
        self.assertTrue('rnai_data_sharing_level' in apilog['diff_keys'])
    
    # TODO: test expire dsl: create a "UserAgreementResource.expire"
        
    def test5_service_activity(self):
        
        logger.info('test5_service_activity...')
        self.test0_create_user();
        
        serviced_user = self.user1
        
        # FIXME: make sure performed_by_username belongs to 
        # ServiceActivityPerformers group
        performed_by_user = self.create_screensaveruser(
            { 'username': 'service_activity_performer'})
        performed_by_user2 = self.create_screensaveruser(
            { 'username': 'service_activity_performer2'})

        service_activity_post = {
            'serviced_user_id': serviced_user['screensaver_user_id'],
            'type': "image_analysis",
            'comments': "test",
            'date_of_activity': "2015-10-27",
            'funding_support': "clardy_grants",
            'performed_by_user_id': performed_by_user['screensaver_user_id'],
        }

        # 1. Create        
        resource_uri = BASE_URI_DB + '/serviceactivity'
        resp = self.api_client.post(
            resource_uri, 
            format='json', 
            data=service_activity_post, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))

        data_for_get = { 'limit': 0, 'includes': ['*'] }
        resp = self.api_client.get(
            resource_uri,
            format='json', 
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('new service activity result: %r', new_obj)
        self.assertTrue(API_RESULT_DATA in new_obj)
        new_obj = new_obj[API_RESULT_DATA][0]
        result,msgs = assert_obj1_to_obj2(
            service_activity_post, new_obj)
        self.assertTrue(result,msgs)

        # 1a. Test apilog
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'serviceactivity', 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.debug('serviceactivity log: %r', apilog)
        self.assertTrue(apilog['api_action'] == 'CREATE')
        self.assertEquals(
            apilog['uri'], 
            'screensaveruser/{serviced_user_id}/serviceactivity/{activity_id}'
                .format(**new_obj))
        
        # TODO: test with a serviced screen
        
        # 2. patch
        service_activity_post = {
            'activity_id': new_obj['activity_id'],
            'performed_by_user_id': performed_by_user2['screensaver_user_id']}

        resource_uri = BASE_URI_DB + '/serviceactivity'
        resp = self.api_client.patch(
            resource_uri, 
            format='json', 
            data=service_activity_post, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))

        data_for_get = { 'limit': 0, 'includes': ['*'] }
        resp = self.api_client.get(
            resource_uri,
            format='json', 
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('new service activity result: %r', new_obj)
        self.assertTrue(API_RESULT_DATA in new_obj)
        self.assertEquals(1, len(new_obj[API_RESULT_DATA]))
        new_obj = new_obj[API_RESULT_DATA][0]
        result,msgs = assert_obj1_to_obj2(
            service_activity_post, new_obj)
        self.assertTrue(result,msgs)
        
        # 2.a Test apilog
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'serviceactivity',
            'diff_keys': 'performed_by_user_id' 
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        self.assertEqual(len(apilogs),1, 
            'wrong number of apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('serviceactivity log: %r', apilog)
        self.assertTrue(apilog['api_action'] == 'PATCH')
        self.assertEquals(
            apilog['uri'], 
            'screensaveruser/{serviced_user_id}/serviceactivity/{activity_id}'
                .format(**new_obj))
        
        # 3 delete serviceactivity
        resource_uri = '/'.join([
            BASE_URI_DB,'serviceactivity',str(new_obj['activity_id'])])
        resp = self.api_client.delete(
            resource_uri, authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code == 204, 
            (resp.status_code, self.get_content(resp)))
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(),
            data={ 'limit': 0, 'includes': '*'} )
        self.assertTrue(
            resp.status_code == 404, 
            ('error, publication should be deleted', resp.status_code, 
                self.get_content(resp)))

        # 3.a Test delete apilog
        resource_uri = BASE_REPORTS_URI + '/apilog'
        data_for_get={ 
            'limit': 0, 
            'ref_resource_name': 'serviceactivity',
            'api_action': 'DELETE'
        }
        apilogs = self.get_list_resource(
            resource_uri, data_for_get=data_for_get )
        self.assertTrue(
            len(apilogs) == 1, 'too many apilogs found: %r' % apilogs)
        apilog = apilogs[0]
        logger.info('serviceactivity log: %r', apilog)
        self.assertTrue(apilog['api_action'] == 'DELETE')
        self.assertEquals(
            apilog['uri'], 
            'screensaveruser/{serviced_user_id}/serviceactivity/{activity_id}'
                .format(**new_obj))


class DataSharingLevel(DBResourceTestCase):
    
    def __init__(self, *args, **kwargs):
        DBResourceTestCase.__init__(self, *args, **kwargs)
        self.sr_serializer = ScreenResultSerializer()
        self.general_user_password = 'testpass1'

    def tearDown(self):
        DBResourceTestCase.tearDown(self)
        # NOTE: tearDown may be eliminated for iterative testing:
        # the setup_data method will check for existence of data before creating
        
        logger.info('tearDown...')
        with connection.cursor() as cursor:
            try:
                cursor.execute('delete from well_query_index;')
            except Exception as e:
                logger.exception('on delete well_query_index')
            try:
                cursor.execute('delete from well_data_column_positive_index;')
            except Exception as e:
                logger.exception('on delete well_data_column_positive_index')
            try:
                cursor.execute('delete from cached_query;')
            except Exception as e:
                logger.exception('on delete cached_query')
        Screen.objects.all().delete()
        ScreensaverUser.objects.all().exclude(username=self.username).delete()
        LabAffiliation.objects.all().delete()
        Library.objects.all().delete()
        ApiLog.objects.all().delete()

    
    def setup_data(self):
        # NOTE: tearDown may be eliminated for iterative testing:
        # the setup_data method will check for existence of data before creating
        self.setup_library()
        # Create 2 labs for each level
        def set_user_password(user):
            # assign password to the test user
            # NOTE: may only be done through the Django Model, for now
            # TODO: superuser should be able to assign password through secure connection
            userObj = User.objects.get(username=user['username'])
            userObj.set_password(self.general_user_password)
            userObj.save()
            return user

        # Store extant users so that they will not have to be recreated if 
        # already created
        reference_users = self.get_list_resource(BASE_URI + '/screensaveruser') 
        reference_users = { 
            user['username']: user 
                for user in reference_users}
        logger.info('starting reference users: %r', reference_users.keys())
        reference_screens = self.get_list_resource(BASE_URI + '/screen') 
        reference_screens = { 
            screen['facility_id']: screen 
                for screen in reference_screens}
        logger.info('starting reference screens: %r', reference_screens.keys())
        user_data = {
        }
        for dsl in range(1,4):
            for lab in ['a','b']:
                _data = dict(
                    user_data,
                    sm_data_sharing_level=dsl)
                lab_head_username = 'lab_head%d%s' % (dsl,lab)
                lab_head = reference_users.get(lab_head_username, None)
                if lab_head is None:
                    _data = dict(_data,username=lab_head_username)
                    logger.info('setup: user not found; creating %r...', _data)
                    lab_head = set_user_password(self.create_lab_head(_data))
                
                _data = dict(_data,lab_head_id=lab_head['screensaver_user_id'])
                lead_screener_username = 'lead_screener%d%s' % (dsl,lab)
                lead_screener = reference_users.get(lead_screener_username, None)
                if lead_screener is None:
                    _data = dict(_data, username=lead_screener_username)
                    logger.info('setup: user not found; creating %r...', _data)
                    lead_screener = set_user_password(self.create_screensaveruser(_data))
                    
                collaborator_username = 'collaborator%d%s' % (dsl,lab)
                collaborator = reference_users.get(collaborator_username, None)
                if collaborator is None:
                    _data = dict(_data, username=collaborator_username)
                    logger.info('setup: user not found; creating %r...', _data)
                    collaborator = set_user_password(self.create_screensaveruser(_data))
                
                # Create screens for each lab, by level
                screen_data = {
                    'screen_type': 'small_molecule',
                    'lab_head_id': lab_head['screensaver_user_id'],
                    'lead_screener_id': lead_screener['screensaver_user_id'],
                    'collaborator_ids': [collaborator['screensaver_user_id'],],
                }
                
                for screen_dsl in range(0,4):
                    facility_id = '%d%s%d' % (dsl,lab,screen_dsl)
                    if facility_id not in reference_screens:
                        _data = dict(screen_data,
                            facility_id = '%d%s%d' % (dsl,lab,screen_dsl),
                            title = '%d%s%d' % (dsl,lab,screen_dsl),
                            data_sharing_level = screen_dsl )
                        logger.info('setup: screen not found; creating %r...', _data)
                        screen = self.create_screen(
                            data=_data, uri_params=['override=true',])

        reference_users = self.get_list_resource(BASE_URI + '/screensaveruser') 
        reference_users_by_id = {
            user['screensaver_user_id']: user for user in reference_users}
        reference_users = { 
            user['username']: user for user in reference_users}
        self.reference_users = reference_users
        logger.info('reference_users: %r', reference_users.keys())
        self.users_by_level = defaultdict(set)
        for username,user in reference_users.items():
            self.users_by_level[user['sm_data_sharing_level']].add(username)
        logger.info('users_by_level: %r', self.users_by_level)
        
        reference_screens = self.get_list_resource(BASE_URI + '/screen') 
        reference_screens = { 
            screen['facility_id']: screen 
                for screen in reference_screens}
        self.reference_screens = reference_screens
        logger.info('reference_screens: %r', reference_screens.keys())
        screens_by_lead = defaultdict(set)
        for facility_id,screen in reference_screens.items():
            lead_screener = reference_users_by_id[screen['lead_screener_id']]
            screens_by_lead[lead_screener['username']].add(facility_id)
        # organize by level
        self.screens_by_lead = {}
        for username, screens in screens_by_lead.items():
            self.screens_by_lead[username] = sorted(
                screens, key=lambda facility_id: reference_screens[facility_id]['data_sharing_level'])
        logger.info('screens_by_lead: %r', self.screens_by_lead)
        self.screens_by_level = defaultdict(set)
        for facility_id,screen in reference_screens.items():
            self.screens_by_level[screen['data_sharing_level']].add(facility_id)
        logger.info('screens_by_level: %r', self.screens_by_level)
        
        # IF setup data were preserved; delete the screen results
        self.clear_out_screenresults()
        
    def clear_out_screenresults(self):        
        logger.info(
            'remove all screen results before tests')
        reference_screens = self.get_list_resource(BASE_URI + '/screen') 
        reference_screens = { 
            screen['facility_id']: screen 
                for screen in reference_screens}
        for screen_facility_id,screen in reference_screens.items():
            if screen['has_screen_result'] in [1,2]:
                logger.info('delete screen result: %r', screen_facility_id)
                delete_url = '/'.join([
                    BASE_URI, 'screenresult',screen_facility_id])
                self.api_client.delete(delete_url, authentication=self.get_credentials())
        
    def get_entity(self, resource_uri, as_username=None, data_for_get=None):
        
        _data_for_get = { 'limit': 0, 'includes': ['*'] }
        if data_for_get is not None:
            _data_for_get.update(data_for_get)
        
        if as_username is not None:
            authentication = self.create_basic(
                as_username, self.general_user_password)
        else:
            authentication = self.get_credentials()
            
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=authentication, 
            data=_data_for_get)
        if resp.status_code == 200:
            new_obj = self.deserialize(resp)
            if API_RESULT_DATA in new_obj:
                self.assertEqual(len(new_obj[API_RESULT_DATA]),1,
                    'more than one object returned for: %r, returns: %r'
                    % (resource_uri,new_obj))
                new_obj = new_obj[API_RESULT_DATA][0]
            logger.debug('obj: %r', new_obj)
            return new_obj
        elif resp.status_code == 404:
            return None
        else:
            self.fail((resp.status_code,self.get_content(resp)))
    
    def get_user(self, username, as_username=None, data_for_get=None):
        
        resource_uri = '/'.join([BASE_URI_DB, 'screensaveruser', username])
        return self.get_entity(resource_uri, as_username, data_for_get)
    
    def get_screen(self, facility_id, as_username=None, data_for_get=None):
        
        resource_uri = '/'.join([BASE_URI_DB, 'screen', facility_id])
        return self.get_entity(resource_uri, as_username, data_for_get)
    
    def setup_library(self):
        
        short_name = 'library1'
        resource_uri = '/'.join([BASE_URI_DB, 'library',short_name])
        self.library1 = self.get_entity(resource_uri)
        if self.library1 is None:
            logger.info('library not found, creating: %s', short_name)
            self.library1 = self.create_library({
                'short_name': short_name,
                'start_plate': 1, 
                'end_plate': 20,
                'screen_type': 'small_molecule' })
    
            logger.info('create wells...')
            plate = 1
            input_data = [
                self.create_small_molecule_test_well(
                    plate,i,library_well_type='experimental') 
                for i in range(0,20)]
            # setup one control well: i:20 - E02
            input_data.append(
                self.create_small_molecule_test_well(
                    plate,20,library_well_type='empty') 
                )
            resource_name = 'well'
            resource_uri = '/'.join([
                BASE_URI_DB,'library', self.library1['short_name'],resource_name])
            resp = self.api_client.put(
                resource_uri, format='sdf', data={ 'objects': input_data } , 
                authentication=self.get_credentials(), 
                **{ 'limit': 0, 'includes': '*'} )
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))

    @staticmethod
    def _create_screen_result_test_data(
        screen_facility_id, well_ids, 
        confirmed_positive_wells=None, partition_positive_wells=None):
        
        confirmed_positive_wells = confirmed_positive_wells or []
        partition_positive_wells = partition_positive_wells or []   
        
        fields = {
            'E': {
                'name': 'Field1_non_positive_indicator',
                'data_worksheet_column': 'E',
                'data_type': 'numeric',
                'decimal_places': 4, 
                'description': 'Non-positive field',
                'replicate_ordinal': 2,
                'is_follow_up_data': True,
                'assay_readout_type': 'flourescence_intensity',
            },
            'F': {
                'name': 'Field2_positive_indicator',
                'data_worksheet_column': 'F',
                'data_type': 'confirmed_positive_indicator',
                'description': 'positive indicator field',
                'how_derived': 'z-score > 1 for Field3'
            },
            'G': {
                'name': 'Field3_partitioned_positive',
                'data_worksheet_column': 'G',
                'data_type': 'partition_positive_indicator',
                'description': 'partition positive indicator',
                'how_derived': 'ranking of z-score <.5, ,.5<=x<=1, >1'
            },
        }
        
        result_values = []
        for i, well_id in enumerate(well_ids):
            
            result_value = {
                'well_id': well_id,
                'E': i+1
            }
            if well_id in confirmed_positive_wells:
                result_value['F'] = 'CP'
            if well_id in partition_positive_wells:
                result_value['G'] = 'W'
            result_values.append(result_value)
        
        input_data = OrderedDict((
            ('meta', {'screen_facility_id': screen_facility_id } ),
            ('fields', fields),
            ('objects', result_values),
        ))
        return input_data


    def create_screen_result_for_test(
            self, screen_facility_id,well_ids, 
            confirmed_positive_wells=None, partition_positive_wells=None):
        
        input_data = self._create_screen_result_test_data(
            screen_facility_id,
            well_ids, confirmed_positive_wells, partition_positive_wells)

        # The ScreenResultSerializer only recognizes the XLSX format:
        # So serialize the input data into an XLSX file
        input_data_put = screen_result_importer.create_output_data(
            screen_facility_id, 
            input_data['fields'], 
            input_data['objects'] )
        input_data_put = self.sr_serializer.serialize(
            input_data_put, XLSX_MIMETYPE)

        data_for_get = { 'limit': 0, 'includes': ['*'] }
        data_for_get['HTTP_AUTHORIZATION'] = self.get_credentials()
        data_for_get['CONTENT_TYPE'] = XLSX_MIMETYPE
        data_for_get['HTTP_ACCEPT'] = XLSX_MIMETYPE
        logger.info('PUT screen result to the server...')

        resource_name = 'screenresult'
        resource_uri = '/'.join([
            BASE_URI_DB,resource_name,screen_facility_id])
        resp = self.django_client.put(
            resource_uri, data=input_data_put, **data_for_get )
        if resp.status_code not in [200, 204]:
            content = self.get_content(resp)
            if content:
                logger.info('resp: %r', 
                    [[str(y) for y in x] 
                        for x in self.serializer.from_xlsx(content)])
        self.assertTrue(
            resp.status_code in [200, 204], resp.status_code)

        output_data = self.get_screenresult(screen_facility_id)
        
        ScreenResultSerializerTest.validate(self, input_data, output_data)

    def get_screenresult_json(self, screen_facility_id, username=None, data=None):
        data_for_get = { 'includes': '*' }
        if data is not None:
            data_for_get.update(data)

        if username:
            authentication=self.create_basic(username, self.general_user_password)
        else:
            authentication=self.get_credentials()

        resource_uri = BASE_URI + '/screenresult/' + screen_facility_id
        resp = self.api_client.get(
            resource_uri, format='json', data=data_for_get, 
            authentication=authentication )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        new_data = new_obj[API_RESULT_DATA]
        return new_data

    def get_screenresult(self,screen_facility_id, username=None, data=None):
        
        if username:
            authentication=self.create_basic(username, self.general_user_password)
        else:
            authentication=self.get_credentials()
        
        data_for_get = { 'limit': 0, 'includes': ['*'] }
        data_for_get['HTTP_AUTHORIZATION'] = authentication
        data_for_get['CONTENT_TYPE'] = XLSX_MIMETYPE
        data_for_get['HTTP_ACCEPT'] = XLSX_MIMETYPE
        resource_name = 'screenresult'
        resource_uri = '/'.join([
            BASE_URI_DB,resource_name,screen_facility_id])
        
        logger.info('GET screen result %r from the server, %r', 
            screen_facility_id, data_for_get)
        resp = self.django_client.get(resource_uri, data=data, **data_for_get)
        if resp.status_code not in [200, 204]:
            if resp.status_code == 403:
                logger.info('resp: %r',self.get_content(resp))
            content = self.get_content(resp)
            if content:
                logger.info('resp: %r', 
                    [[str(y) for y in x] 
                        for x in self.serializer.from_xlsx(content)])
        self.assertTrue(resp.status_code in [200,201,202],resp.status_code)
        output_data = self.sr_serializer.deserialize(
            self.get_content(resp), XLSX_MIMETYPE)
        return output_data

    def get_screens(self, username, data=None):
        
        data_for_get = { 'includes': '*' }
        if data is not None:
            data_for_get.update(data)
        
        # Get a general screen listing
        resource_uri = BASE_URI + '/screen'
        resp = self.api_client.get(
            resource_uri, format='json', data=data_for_get, 
            authentication=self.create_basic(username, self.general_user_password) )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        new_data = new_obj[API_RESULT_DATA]
        return new_data
        
    def get_datacolumns(self, username=None, data=None):
        
        data_for_get = { 'includes': '*' }
        if data is not None:
            data_for_get.update(data)
        
        if username:
            authentication=self.create_basic(username, self.general_user_password)
        else:
            authentication=self.get_credentials()
        
        resource_uri = BASE_URI + '/datacolumn'
        resp = self.api_client.get(
            resource_uri, format='json', data=data_for_get, 
            authentication=authentication )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        new_data = new_obj[API_RESULT_DATA]
        return new_data
        
    def get_schema(self, resource_name, username=None, data=None):
        
        data_for_get = { 'includes': '*' }
        if data is not None:
            data_for_get.update(data)
        
        if username:
            authentication=self.create_basic(username, self.general_user_password)
        else:
            authentication=self.get_credentials()
        
        resource_uri = '/'.join([BASE_URI_DB, resource_name, 'schema'])

        resp = self.api_client.get(
            resource_uri, format='json', authentication=authentication,
            data=data_for_get)
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        return self.deserialize(resp)
        
#     def test1_schema(self):
#         
#         self.setup_users()
# 
#         # Get the Screen schema for the user
#         screen_schema = self.get_schema(
#             'screen', authentication=self.create_basic(
#                 self.lab_head1['username'], self.general_user_password))
#         logger.debug('screen schema: %r', screen_schema)
#         self.assertIsNotNone(screen_schema)
#         
#         fields = screen_schema['fields']
#         
#         # Verify ordering/filtering turned off for level1, 2, 3 fields
#         restricted_fields = set()
#         fields_by_level = defaultdict(set)
#         for key,field in screen_schema['fields'].items():
#             for level in range(0,4):
#                 if field['data_access_level'] == level:
#                     fields_by_level[level].add(key)
#                     if level > 0:
#                         restricted_fields.add(key)
#                 elif field['data_access_level'] == 0 and field['viewGroups'] is None:
#                     fields_by_level[0].add(key)
#         logger.info('fields by level: %r', fields_by_level)
#         
#         # Test some fields (not comprehensive)
#         self.assertTrue(set(['summary','publishable_protocol']) & set(fields_by_level[1]))
#         self.assertTrue(set(['positives_summary','assay_readout_types']) & set(fields_by_level[2]))
#         self.assertTrue(set(['activity_count','cherry_pick_screenings']) & set(fields_by_level[3]))
#         
#         for key in restricted_fields:
#             self.assertTrue(
#                 fields[key]['filtering'] is not True, 
#                 'field filtering: %r: %r' % (key, fields[key]))
#             self.assertTrue(
#                 fields[key]['ordering'] is not True, 
#                 'field ordering: %r: %r' % (key, fields[key]))
        
    def test2_overlapping_screens_property(self):
         
        self.setup_data()
        facility1 = self.screens_by_lead['lead_screener2a'][2]
        facility2 = self.screens_by_lead['lead_screener2b'][2]
        screen1_2 = self.get_screen(facility1)
        screen2_2 = self.get_screen(facility2)
 
        # 1.A before data upload, no overlapping of result values
        self.assertIsNone(
            screen1_2['overlapping_positive_screens'],
            screen1_2['overlapping_positive_screens'])
        self.assertIsNone(
            screen2_2['overlapping_positive_screens'],
            screen2_2['overlapping_positive_screens'])
         
        # 2. create overlapping data
        well_ids_to_create = ['00001:A0%d'%i for i in range(1,6)]
        self.create_screen_result_for_test(
            facility1, well_ids_to_create,
            confirmed_positive_wells=['00001:A02'])
        self.create_screen_result_for_test(
            facility2, well_ids_to_create,
            confirmed_positive_wells=['00001:A02'])
 
        screen1_2 = self.get_screen(facility1)
        screen2_2 = self.get_screen(facility2)
         
        self.assertTrue(
            screen2_2['facility_id'] in screen1_2['overlapping_positive_screens'],
            '%r not found in %r' 
                % (screen2_2['facility_id'],screen1_2['overlapping_positive_screens']))
        self.assertTrue(
            screen1_2['facility_id'] in screen2_2['overlapping_positive_screens'],
            '%r not found in %r' 
                % (screen1_2['facility_id'],screen2_2['overlapping_positive_screens']))
        # 3. Remove screen results
        delete_url = '/'.join([BASE_URI, 'screenresult/%s'])
        self.api_client.delete(delete_url % facility1, authentication=self.get_credentials())
        self.api_client.delete(delete_url % facility2, authentication=self.get_credentials())
        
        screen1_2 = self.get_screen(facility1)
        screen2_2 = self.get_screen(facility2)

        self.assertIsNone(
            screen1_2['overlapping_positive_screens'],
            screen1_2['overlapping_positive_screens'])
        self.assertIsNone(
            screen2_2['overlapping_positive_screens'],
            screen2_2['overlapping_positive_screens'])
        
    def test3_visibility(self):
        '''
        Test User data visibility (access level) as per the ICCBL User Agreements:
        
        General Screen visibility Rules:
        
        before deposit: Own Screens and Level 0 Screens
        after deposit: Own Screens and Levels 0,1,2 Screens
        (level 3 users may only ever see Own Screens and Level 0 Screens)
        
        Field Access Level Rules:
        
        User-Screen Access Level 0:
        - for any Screen the User is authorized to see
        User-Screen Access Level 1:
        - Level 1 or 2 Users viewing Level 0 Screens
        - Level 1 Users viewing Level 1 Screens
        - Level 1 Users viewing overlapping level 2 Screens
        - Level 2 Users viewing overlapping level 1 or 2 Screens
        User-Screen Access Level 2 :
        - Level 1 or 2 Users viewing Level 0 Screens
        - Level 1 Users viewing Level 1 Screens
        User-Screen Access Level 3:
        - When viewing own screens
        
        '''
        self.setup_data()
        self.do_level_1_visibility()
        self.do_level_2_visibility()
        
    def do_visibility_test(
            self, fields, user, expected_screens_access_levels):

        logger.info('do_visibility_test for: %r', user['username'])

        # Retrieve the screens as superuser for reference
        reference_screens = self.get_list_resource(BASE_URI + '/screen') 
        reference_screens = { screen['facility_id']: screen 
            for screen in reference_screens}

        # Retrieve screens as the user        
        reported_screens = self.get_screens(user['username'])
        reported_screens = { screen['facility_id']:screen 
            for screen in reported_screens}

        reported_screens_access_levels = defaultdict(set)
        for facility_id, screen in reported_screens.items():
            reported_level = screen['user_access_level_granted']
            reported_screens_access_levels[reported_level].add(facility_id)
        logger.info('expected_screens_access_levels: %r', 
            expected_screens_access_levels)
        logger.info('reported_screens_access_levels: %r', 
            reported_screens_access_levels)
        # Note filter default values added to the defaultdict inadvertantly
        self.assertEqual(
            set([key for key,value in expected_screens_access_levels.items() if value ]),
            set(reported_screens_access_levels.keys()))

        for access_level_expected, ids_expected in expected_screens_access_levels.items():
            self.assertEqual(
                ids_expected, 
                reported_screens_access_levels[access_level_expected])
            for facility_id in ids_expected:
                logger.info('do visibility test, screen: %r', facility_id)
                logger.info(
                    'test screen: %r, access level: %r', 
                    facility_id, access_level_expected)
                reported_screen = reported_screens[facility_id]
                reference_screen = reference_screens[facility_id]
                self.assertEqual(
                    access_level_expected, 
                    reported_screen['user_access_level_granted'])
                
                for key, value in reported_screen.items():
                    if key == 'user_access_level_granted':
                        continue
                    if key not in fields:
                        if value is None:
                            continue
                        logger.warn('fields missing key: %r, reported: %r',
                            key, value)
                    
                    field = fields[key]
                    reference_value = reference_screen[key]
                    if reference_value is None:
                        # skip null fields for now
                        continue
                    logger.debug('screen: %r, test value: %r: %r', 
                        facility_id, key, reference_value)
                    if field['data_access_level'] > access_level_expected:
                        self.assertIsNone(
                            value, '%r - %r:%r' % (facility_id, key,value))
                        continue
                    if field.get('view_groups', None):
                        usergroups = user.get('usergroups')
                        
                        if ( not usergroups
                             or not set(field['view_groups']) & set(usergroups)):
                            logger.info('view group restricted field: %r: %r', 
                                key, user['username'])
                            logger.debug(
                                'view group restricted field: '
                                '%r, user %r, usergroups: %r, view_groups: %r', 
                                key, user, user.get('usergroups'), field['view_groups'])
                            continue
                    if key == 'has_screen_result':
                        if reference_value == 1:
                            if access_level_expected in [2,3]:
                                self.assertEqual(value, 1)
                            else:
                                self.assertEqual(value, 2) # "not shared"
                        else:
                            self.assertEqual(value, 3) # "none"
                    elif key == 'overlapping_positive_screens':
                        if access_level_expected == 2:
                            expected_overlapping = set(reference_value)
                            expected_overlapping &= (
                                expected_screens_access_levels[2]
                                | expected_screens_access_levels[3])
                            self.assertEqual(expected_overlapping, set(value))
                        elif access_level_expected == 3:
                            expected_overlapping = set(reference_value)
                            expected_overlapping &= (
                                expected_screens_access_levels[1]
                                | expected_screens_access_levels[2]
                                | expected_screens_access_levels[3])
                            self.assertEqual(expected_overlapping, set(value))
                        else:
                            self.assertIsNone(value)
                    else:
                        self.assertEqual(reference_value, value,
                            '%r:%r: %r != %r' % (
                                facility_id,key,reference_value,value))

    def do_screen_result_visible_test(self, user, screen_facility_id):
        reference_screen_result = self.get_screenresult(screen_facility_id)
        reported_screen_results = self.get_screenresult(
            screen_facility_id, user['username'])
        
        ScreenResultSerializerTest.validate(
            self, reference_screen_result, reported_screen_results)
    
    def do_data_columns_test(self, user, expected_screens_access_levels,
        expected_positive_count_by_col_title=None):
        
        # Retrieve the screens as superuser for reference
        reference_screens = self.get_list_resource(BASE_URI + '/screen') 
        reference_screens = { screen['facility_id']: screen 
            for screen in reference_screens}

        reference_data_columns = self.get_datacolumns()
        logger.info('reference datacolumns: %r', 
            [(dc['screen_facility_id'],dc['name'],dc['user_access_level_granted'])
                for dc in reference_data_columns])
        reference_data_columns_by_title = { dc['title']:dc for dc in reference_data_columns }
        if expected_positive_count_by_col_title:
            for title, count in expected_positive_count_by_col_title.items():
                logger.info('test pos count for: %r, %r, %r',
                    title, dc['key'], dc['positives_count'])
                dc = reference_data_columns_by_title[title]
                self.assertEqual(dc['positives_count'], count)
        
        # 1. Collate the datacolumns reported by expected access level

        dc_by_screen = defaultdict(list)
        for dc in reference_data_columns:
            dc_by_screen[dc['screen_facility_id']].append(dc)
        if not dc_by_screen:
            self.fail('no reference datacolumns found')
        for facility_id,dcs in dc_by_screen.items():
            logger.info('dcs reported: %r: %r', facility_id,
                [dc['name'] for dc in dcs])
            
        expected_cols_access_levels = defaultdict(set)
        for access_level,facility_ids in expected_screens_access_levels.items():
            if access_level == 0:
                continue
            for facility_id in facility_ids:
                screen = reference_screens.get(facility_id) 
                if screen['has_screen_result'] == 3:
                    continue
                dcs = dc_by_screen.get(facility_id,None)
                if not dcs:
                    self.fail('no datacolumns found: %r' % facility_id)
                for dc in dcs:
                    logger.debug(
                        'access_level: %r, dc: %r', 
                        access_level, dc['data_column_id'])
                    if access_level > 1:
                        expected_cols_access_levels[access_level].add(
                            dc['data_column_id'])
                    elif access_level == 1:
                        # must be a positives column with positive values
                        if dc.get('positives_count',0) > 0:
                            logger.info('level 1 pos column: %r, %r, %r', 
                                dc['key'],dc['title'], dc['positives_count'])
                            expected_cols_access_levels[access_level].add(
                                dc['data_column_id'])
                        else:
                            logger.info(
                                'level 1 column not allowed, no positives: %r, %r',
                                dc['key'], dc['title'])
                    else:
                        logger.info(
                            'data column not allowed: %r:%r:%r',
                            dc['screen_facility_id'], dc['data_column_id'],
                            dc['name'])

        # 2. compare expected datacolumns by level to reported
        
        reported_datacolumns = self.get_datacolumns(username=user['username'])
        
        dc_by_screen = defaultdict(list)
        for dc in reported_datacolumns:
            dc_by_screen[dc['screen_facility_id']].append(dc['name'])
        logger.info('reported dc_by_screen: %r', dc_by_screen)
        
        reported_dc_by_access_level = defaultdict(set)
        for dc in reported_datacolumns:
            reported_dc_by_access_level[dc['user_access_level_granted']]\
                .add(dc['data_column_id'])
        logger.info('reported_dc_by_access_level: %r', reported_dc_by_access_level) 

        for access_level,col_ids in expected_cols_access_levels.items():
            logger.info('access level: %r, expected: %r', access_level, col_ids)
            reported_col_ids = reported_dc_by_access_level.get(access_level, None)
            self.assertEqual(col_ids, reported_col_ids)        

    def do_screenresult_test(self, user, screen_facility_id, 
        expected_access_level, mutual_wells=None):
        '''
        Test the screen results; after the do_data_columns_test has been run to 
        verify visibility of columns
        '''
        logger.info('do_screenresult_test: %r, %r', user, screen_facility_id)
        
        
        reference_screen = self.get_screen(screen_facility_id)
        
        reported_screen = self.get_screen(screen_facility_id, user['username'])
        
        user_access_level_granted = reported_screen['user_access_level_granted']
        logger.info('screen: %r, user accessing: %r, access_level: %r',
            screen_facility_id, user['username'], user_access_level_granted )
        self.assertEqual(user_access_level_granted, expected_access_level)
        
        key = 'has_screen_result'
        reference_value = reference_screen[key]
        if reference_value == 3:
            self.fail('do_screenresult_test: screen: %r reports no screen result: %r'
                % (screen_facility_id, reference_screen))
        self.assertEqual(1, reference_value, 
            'admin user should never see has_screen_result other than (1,2): %r'
            % reference_value)
        if user_access_level_granted < 2:
            self.assertEqual(2, reported_screen[key])
            
            # TODO: verify that the screen result can not be directly accessed
        else:
            self.assertEqual(1, reported_screen[key])
            
            data = { 'screen_type': reference_screen['screen_type'] }
            reference_datacolumns = self.get_datacolumns(data=data)
            reference_datacolumns_by_key = {
                dc['key']: dc for dc in reference_datacolumns }
            reference_datacolumns = { dc['data_column_id']:dc 
                for dc in reference_datacolumns}
            reported_datacolumns = self.get_datacolumns(
                username=user['username'], data=data)
            reported_datacolumns_by_key = {
                dc['key']: dc for dc in reported_datacolumns }
            reported_datacolumns = { dc['data_column_id']:dc 
                for dc in reported_datacolumns }
            dc_by_screen = defaultdict(list)
            for dc_id, dc in reported_datacolumns.items():
                dc_by_screen[dc['screen_facility_id']].append(dc_id)
            logger.info('reported dc_by_screen: %r', dc_by_screen)
            
            reported_dc_by_access_level = defaultdict(set)
            for dc_id,dc in reported_datacolumns.items():
                reported_dc_by_access_level[dc['user_access_level_granted']]\
                    .add(dc_id)
            logger.info('reported_dc_by_access_level: %r', reported_dc_by_access_level) 
            
            # Must check the schema first, and build the expected fields from that
            
            sr_schema_resource = '/'.join(['screenresult',screen_facility_id])
            
            # grab the "maximal" screen result: all available columns from all screens added
            dc_ids = reference_datacolumns.keys()
            
            reference_expected_dc_keys = set([reference_datacolumns[dc_id]['key'] 
                for dc_id in dc_ids ])
            reference_screen_result_schema = self.get_schema(
                sr_schema_resource, data={ 'dc_ids':dc_ids })
            logger.info('verify that all of the datacolumns are present in the reference')
            logger.info('reported reference fields: %r', 
                reference_screen_result_schema['fields'].keys())
            reference_dc_keys = set([field['key'] 
                for field in reference_screen_result_schema['fields'].values()
                    if field.get('is_datacolumn',False)])
            self.assertEqual(reference_expected_dc_keys, reference_dc_keys)
            
            reported_screen_result_schema = self.get_schema(
                sr_schema_resource, username=user['username'], data={ 'dc_ids':dc_ids })
            self.assertIsNotNone(reported_screen_result_schema)
             
            expected_dc_keys = set()
            for dc_id, dc in reported_datacolumns.items():
                # Note: positives_count is a level 2 field, so must use ref col
                reference_dc = reference_datacolumns[dc_id]
                if dc_id in reported_dc_by_access_level[3]:
                    expected_dc_keys.add(dc['key'])
                elif dc_id in reported_dc_by_access_level[2]:
                    expected_dc_keys.add(dc['key'])
                elif dc_id in reported_dc_by_access_level[1]:
                    logger.info('level 1 column: %r, %r, %r, %r, %r', 
                        dc['key'],dc['title'], reference_dc['positives_count'], 
                        dc['screen_facility_id'],
                        reported_screen['overlapping_positive_screens'])
                    if ( dc['screen_facility_id'] 
                            in reported_screen['overlapping_positive_screens']):
                        if reference_dc['positives_count'] > 0:
                            logger.info('allowed')
                            expected_dc_keys.add(dc['key'])    
            reported_dc_keys = set([field['key'] 
                for field in reported_screen_result_schema['fields'].values()
                    if field.get('is_datacolumn',False)])
            self.assertEqual(expected_dc_keys, reported_dc_keys)
                        
            # Test the results
            
            reference_screenresult = self.get_screenresult_json(
                screen_facility_id, data={ 'dc_ids': dc_ids} )
            
            for rv in reference_screenresult:
                self.assertTrue(reference_expected_dc_keys <= set(rv.keys()),
                    'keys not found %r' 
                        % (reference_expected_dc_keys-set(rv.keys())))
            reference_result_values = {
                rv['well_id']:rv for rv in reference_screenresult }
            
            reported_screenresult = self.get_screenresult_json(
                screen_facility_id, username=user['username'], 
                data={ 'dc_ids': dc_ids} )
            for rv in reported_screenresult:
                well_id = rv['well_id']
                self.assertTrue(expected_dc_keys <= set(rv.keys()),
                    '%r, keys not found %r' 
                        % (well_id, expected_dc_keys-set(rv.keys())))
                disallowed_keys = reference_expected_dc_keys - reported_dc_keys
                logger.info('disallowed dc keys: %r', disallowed_keys)
                self.assertTrue(disallowed_keys.isdisjoint(set(rv.keys())),
                    '%r, disallowed keys found %r' 
                        % (well_id, disallowed_keys & set(rv.keys())))
                
                reference_rv = reference_result_values[well_id]
                
                for key,value in rv.items():
                    reference_value = reference_rv[key]
                    dc = reported_datacolumns_by_key.get(key,None)
                    if dc is None:
                        # non-datacolumn fields should be equal
                        self.assertEqual(
                            reference_value, value,
                            '%r, key: %r, %r != %r' 
                                % (well_id, key, reference_value, value))
                    else:
                        access_level = dc['user_access_level_granted']
                        logger.info('test: %r, %r:%r, %r, %r', 
                            well_id, key, reference_value, value, rv['is_positive'])
                        if access_level > 1:
                            self.assertEqual(
                                reference_value, value,
                                '%r, key: %r, %r != %r' 
                                    % (well_id, key, reference_value, value))
                        else:
                            if reference_value is not None:
                                if rv['is_positive'] is True:
                                    self.assertEqual(
                                        reference_value, value,
                                        '%r: %r, %r != %r' 
                                        % (well_id, key, reference_value, value))
                                else:
                                    logger.info('non overlapping: %r: %r: %r', 
                                        rv['well_id'], key, value)
                                    self.assertIsNone(value,
                                        '%r: %r, non-overlapping should be none: %r: %r' 
                                        % (well_id, key, reference_value, value))
                            else:
                                logger.info('reference value is none: %r, %r:%r,%r',
                                    well_id, key, reference_value, value)
    
    def get_fields_by_level(self, fields):
        fields_by_level = defaultdict(set)
        restricted_fields = set()
        view_group_restricted_fields = set()
        for key,field in fields.items():
            if field.get('view_groups'):
               view_group_restricted_fields.add(key)
               continue 
            for level in range(0,4):
                if field['data_access_level'] == level:
                    fields_by_level[level].add(key)
                elif not field['data_access_level']:
                    fields_by_level[0].add(key)
                else:
                    restricted_fields.add(key)
        logger.info('view group restricted fields: %r', view_group_restricted_fields)
        logger.info('restricted fields: %r', restricted_fields)
        logger.info('fields by level: %r', fields_by_level)
        return fields_by_level
    
    def do_level_1_visibility(self):
        ''' Test level 1 user data visibility '''

        user = self.reference_users['lead_screener1a']
        own_screens = self.screens_by_lead[user['username']]
        logger.info('screens by lead: %r', self.screens_by_lead)

        logger.info('user: %r, own screens: %r', user['username'], own_screens)
        
        # Get the Screen schema for the user
        screen_schema = self.get_schema('screen', user['username'])
        self.assertIsNotNone(screen_schema)
        
        fields_by_level = self.get_fields_by_level(screen_schema['fields'])
        
        logger.info(
            '1. Visibility before data deposit (own & level 0 screens only)')
        
        expected_screens_access_levels = defaultdict(set)
        expected_screens_access_levels[2] = \
            set(self.screens_by_level[0]) - set(own_screens)
        expected_screens_access_levels[3] = set(own_screens)
        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)

        logger.info(
            '2. After data deposit')
        well_ids_to_create = ['00001:A0%d'%i for i in range(1,6)]
         
        logger.info('2.a non-qualifying deposit (level 2 screen for level 1 user)')
        # - Screen visibility should be unchanged
        screen_facility_id = own_screens[2]
        self.create_screen_result_for_test(
            screen_facility_id, well_ids_to_create,
            confirmed_positive_wells=['00001:A02'])
        expected_positive_count_by_col_title = {
            '%s [%s]' % ('Field2_positive_indicator',screen_facility_id): 1 }
        # test that user can get own screen result
        
        screen_result = self.get_screenresult(
            screen_facility_id, user['username'])
 
        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)
        
        logger.info(
            '2.b qualifying deposit, can now see data for level 1,2 screens')
        screen_facility_id = own_screens[1]
        self.create_screen_result_for_test(
            screen_facility_id, well_ids_to_create,
            confirmed_positive_wells=['00001:A01','00001:A02'])
        expected_positive_count_by_col_title[
            '%s [%s]' % ('Field2_positive_indicator',screen_facility_id)] = 2

        expected_screens_access_levels = defaultdict(set)
        expected_screens_access_levels[0] = \
            set(self.screens_by_level[2]) - set(own_screens)
        expected_screens_access_levels[2] = set(self.screens_by_level[0])
        expected_screens_access_levels[2].update(self.screens_by_level[1])
        expected_screens_access_levels[2] -= set(own_screens)
        expected_screens_access_levels[3] = set(own_screens)
        expected_screens_access_levels_after_deposit_no_overlap = \
            expected_screens_access_levels
        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)
        self.do_data_columns_test(
            user, expected_screens_access_levels,
            expected_positive_count_by_col_title=expected_positive_count_by_col_title)
        
        logger.info(
            '2.c Create overlapping data in a level 2 Screen')
        screen_facility_id = self.screens_by_lead['lead_screener2a'][2]
        self.create_screen_result_for_test(
            screen_facility_id, well_ids_to_create,
            confirmed_positive_wells=['00001:A02','00001:A03'])
        expected_positive_count_by_col_title[
            '%s [%s]' % ('Field2_positive_indicator',screen_facility_id)] = 2

        expected_screens_access_levels = defaultdict(set)
        expected_screens_access_levels[0] = \
            set(self.screens_by_level[2]) - set(own_screens)
        expected_screens_access_levels[0].remove(screen_facility_id)
        expected_screens_access_levels[1] = set([screen_facility_id])
        expected_screens_access_levels[2] = set(self.screens_by_level[0])
        expected_screens_access_levels[2].update(self.screens_by_level[1])
        expected_screens_access_levels[2] -= set(own_screens)
        expected_screens_access_levels[3] = set(own_screens)
        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)

        self.do_data_columns_test(
            user, expected_screens_access_levels,
            expected_positive_count_by_col_title=expected_positive_count_by_col_title)
        self.do_screenresult_test(user, screen_facility_id, 1) # should be level 1 access
        logger.info('do level 1 user viewing level 2 mutual wells test')
        self.do_screenresult_test(user, own_screens[1], 3,
            mutual_wells=['00001:A02'])
#         self.do_screen_result_overlapping_test(
#             user, self.screens1[1]['facility_id'])

        logger.info(
            '3 remove screen result')
        delete_url = '/'.join([
            BASE_URI, 'screenresult',screen_facility_id])
        self.api_client.delete(delete_url, authentication=self.get_credentials())

        self.do_visibility_test(
            screen_schema['fields'], user, 
            expected_screens_access_levels_after_deposit_no_overlap)
        self.do_data_columns_test(user, expected_screens_access_levels_after_deposit_no_overlap)

        self.clear_out_screenresults()
        
    def do_level_2_visibility(self):
        ''' Test level 2 user data visibility '''

        user = self.reference_users['lead_screener2a']
        own_screens = self.screens_by_lead[user['username']]
        logger.info('screens by lead: %r', self.screens_by_lead)
        logger.info('user: %r, own screens: %r', user['username'], own_screens)
        
        # Get the Screen schema for the user
        screen_schema = self.get_schema('screen', user['username'])
        self.assertIsNotNone(screen_schema)

        fields_by_level = self.get_fields_by_level(screen_schema['fields'])
        
        logger.info(
            '1. Visibility before data deposit (own & level 0 screens only)')
        
        expected_screens_access_levels = defaultdict(set)
        expected_screens_access_levels[2] = \
            set(self.screens_by_level[0]) - set(own_screens)
        expected_screens_access_levels[3] = set(own_screens)
        logger.info('expected_screens_access_levels: %r',
            expected_screens_access_levels)

        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)

        logger.info('2. After data deposit')
        well_ids_to_create = ['00001:A0%d'%i for i in range(1,6)]
         
        logger.info(
            '2.a non-qualifying deposit (level 1 screen for level 2 user)')
        # - Screen visibility should be unchanged
        screen_facility_id = own_screens[1]
        self.create_screen_result_for_test(
            screen_facility_id, well_ids_to_create,
            confirmed_positive_wells=['00001:A03'])

        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)
        
        logger.info(
            '2.b qualifying deposit')
        # - can now see data for level 1,2 screens
        screen_facility_id = own_screens[2]
        self.create_screen_result_for_test(
            screen_facility_id, well_ids_to_create,
            confirmed_positive_wells=['00001:A03'])
        
        expected_screens_access_levels = defaultdict(set)
        expected_screens_access_levels[0] = set(self.screens_by_level[2]) 
        expected_screens_access_levels[0].update(self.screens_by_level[1])
        expected_screens_access_levels[0] -= set(own_screens)
        expected_screens_access_levels[2] = set(self.screens_by_level[0])
        expected_screens_access_levels[2] -= set(own_screens)
        expected_screens_access_levels[3] = set(own_screens)
        expected_screens_access_levels_after_deposit_no_overlap = \
            expected_screens_access_levels
        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)
        
        logger.info(
            '2.c Create overlapping data in a level 2 Screen')
        screen_facility_id = self.screens_by_lead['lead_screener2b'][2]
        self.create_screen_result_for_test(
            screen_facility_id, well_ids_to_create,
            confirmed_positive_wells=['00001:A03'])

        expected_screens_access_levels = defaultdict(set)
        expected_screens_access_levels[0] = set(self.screens_by_level[2]) 
        expected_screens_access_levels[0].update(self.screens_by_level[1])
        expected_screens_access_levels[0] -= set(own_screens)
        expected_screens_access_levels[0].remove(screen_facility_id)
        expected_screens_access_levels[1] = set([screen_facility_id])
        expected_screens_access_levels[2] = set(self.screens_by_level[0])
        expected_screens_access_levels[2] -= set(own_screens)
        expected_screens_access_levels[3] = set(own_screens)
        self.do_visibility_test(
            screen_schema['fields'], user, expected_screens_access_levels)

        self.do_data_columns_test(user, expected_screens_access_levels)

        logger.info(
            '3 remove screen result')
        delete_url = '/'.join([
            BASE_URI, 'screenresult',screen_facility_id])
        self.api_client.delete(
            delete_url, authentication=self.get_credentials())
        self.do_visibility_test(
            screen_schema['fields'], user, 
            expected_screens_access_levels_after_deposit_no_overlap)
    
        self.clear_out_screenresults()
        
