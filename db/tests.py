"""
This file demonstrates writing tests using the unittest module. These will pass
when you run "manage.py test".

Replace this with more appropriate tests for your application.
"""

from django.test import TestCase
from reports.tests import MetaHashResourceBootstrap
from tastypie.test import ResourceTestCase, TestApiClient
import os
import logging
from lims.tests import assert_obj1_to_obj2, find_all_obj_in_list, find_obj_in_list
import json
import factory
import db.models

from test.factories import *


# from tests.factories import ScreenFactory

logger = logging.getLogger(__name__)

BASE_URI_DB = '/db/api/v1'


class LibraryContentLoadTest(MetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        print '============== LibraryContentLoadTest setup ============'
        super(LibraryTest, self).setUp()
        super(LibraryTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, and the resource definitions
        super(LibraryTest, self)._bootstrap_init_files()
        print '============== LibraryContentLoadTest setup: begin ============'
        
    def test_load_sdf(self):
        pass;
        

class LibraryTest(MetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        print '============== LibraryTest setup ============'
        super(LibraryTest, self).setUp()
        super(LibraryTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, and the resource definitions
        super(LibraryTest, self)._bootstrap_init_files()
        print '============== LibraryTest setup: begin ============'
        
        testApiClient = TestApiClient(serializer=self.csv_serializer) 

        filename = os.path.join(self.directory,'metahash_fields_library.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields:library'})

        print '============== LibraryTest setup: done ============'
        

    def test1_create_library(self):
        logger.info(str(('==== test1_create_library =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        
        logger.info(str((library_item)))
        resp = self.api_client.post(resource_uri, 
            format='json', data=library_item, authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        # create a second library
        library_item = LibraryFactory.attributes()
        
        logger.info(str(('item', library_item)))
        resp = self.api_client.post(resource_uri, 
            format='json', data=library_item, authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        resp = self.api_client.get( resource_uri, 
            format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(library_item, new_obj['objects'])
        self.assertTrue(result, str(('bootstrap item not found', obj, library_item, new_obj['objects'])))
        logger.info(str(('item found', obj)))


class ScreensTest(MetaHashResourceBootstrap,ResourceTestCase):
    
    def setUp(self):
        print '============== ScreensTest setup ============'
        super(ScreensTest, self).setUp()
        super(ScreensTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, and the resource definitions
        super(ScreensTest, self)._bootstrap_init_files()
        print '============== ScreensTest setup: begin ============'
        
        testApiClient = TestApiClient(serializer=self.csv_serializer) 

        filename = os.path.join(self.directory,'metahash_fields_screen.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields:screen'})

        print '============== ScreensTest setup: done ============'
        

    def test1_create_screen(self):
        logger.info(str(('==== test1_create_screen =====')))
        
        resource_uri = BASE_URI_DB + '/screen'
        
        screen_item = ScreenFactory.attributes()
        
        logger.info(str((screen_item)))
        resp = self.api_client.post(resource_uri, 
            format='json', data=screen_item, authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        screen_item = ScreenFactory.attributes()
        
        logger.info(str(('item', screen_item)))
        resp = self.api_client.post(resource_uri, 
            format='json', data=screen_item, authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        resp = self.api_client.get( resource_uri, 
            format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(screen_item, new_obj['objects'])
        self.assertTrue(result, str(('bootstrap item not found', obj, screen_item, new_obj['objects'])))
        self.assertTrue('facility_id' in obj, 'the facility_id was not created')
        logger.info(str(('item found', obj)))

    def test0_create_screen(self):
        logger.info(str(('==== test_create_screen =====')))
        
        # the simplest of tests, create some simple screens
        self.bootstrap_items = [   
            {
                'facility_id': "1",
                'project_phase': "Primary Screen",
                'screen_type': "Small Molecule",
                'title': "Test screen 1",
                'data_sharing_level': 1,
                'total_plated_lab_cherry_picks': 0,
            },
            {
                'facility_id': "2",
                'project_phase': "Primary Screen",
                'screen_type': "RNAi",
                'title': "Test Screen 2",
                'data_sharing_level': 0,
                'total_plated_lab_cherry_picks': 0,
            },
        ]
        
        resource_uri = BASE_URI_DB + '/screen'
        logger.info(str(('--posting to:', resource_uri)))
        for i,item in enumerate(self.bootstrap_items):
            logger.info(str(('item', item)))         
            resp = self.api_client.post(resource_uri, 
                format='json', data=item, authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
#             self.assertHttpCreated(resp)
            
        logger.info('created items, now get them')
        resp = self.api_client.get( resource_uri, 
            format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        for i,item in enumerate(self.bootstrap_items):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(result, str(('bootstrap item not found', item, new_obj['objects'])))
            self.assertTrue('facility_id' in obj, 'the facility_id was not created')
            self.assertEqual(str((i+1)), obj['facility_id'], str(('expected the facility_id returned to be incremented to ', i+1, obj)) )
            logger.info(str(('item found', obj)))
            assert_obj1_to_obj2(item, obj)

        logger.info(str(('==== test_create_screen done =====')))
        

class MutualScreensTest(MetaHashResourceBootstrap,ResourceTestCase):
    
    def test_mutual_positives_to_screen(self):
        pass
        # create two screens
        # create two screen results
        # set data sharing levels
        # find assay wells that overlap
        
