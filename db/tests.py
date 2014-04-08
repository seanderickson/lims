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
import reports.tests

from test.factories import *
from lims.api import CSVSerializer


logger = logging.getLogger(__name__)


BASE_URI = '/db/api/v1'
BASE_REPORTS_URI = '/reports/api/v1'
import db; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__))
BASE_URI_DB = '/db/api/v1'


class TestApiInit(reports.tests.TestApiInit):
    
    def setUp(self):
        super(TestApiInit, self).setUp()
#         super(TestApiInit, self)._setUp()
#         super(TestApiInit, self).test_0bootstrap_metahash()
        # NOTE: run reports tests to set up the metahash resource
        
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')

        
    def test_2api_init(self):
        
        print '***================ super: db test_2api_init =============== '

        super(TestApiInit, self).test_2api_init()
        
        print '***================ local: db test_2api_init =============== '
        
        
        serializer=CSVSerializer() 
        # todo: doesn't work for post, see TestApiClient.post() method, it is 
        # incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=serializer) 
        
        filename = os.path.join(self.db_directory,'api_init_actions.csv')
        with open(filename) as input_file:
            api_init_actions = serializer.from_csv(input_file.read(), root=None)
            
            bootstrap_files = [
                'metahash_resource_data.csv',
                'vocabularies_data.csv']
            for action in api_init_actions:
                
                print '\n++++=========== processing action', json.dumps(action)
                command = action['command'].lower() 
                resource = action['resource'].lower()
                resource_uri = BASE_REPORTS_URI + '/' + resource
                
                if command == 'delete':
                    resp = testApiClient.delete(
                        resource_uri, authentication=self.get_credentials())
                    self.assertHttpAccepted(resp)
                
                else:
                    filename = os.path.join(self.db_directory,action['file'])
                    search_excludes = []
                    # exclude 'resource_uri' from equivalency check during 
                    # bootstrap, because resource table needs to be loaded for
                    # the uri generation
                    if action['file'] in bootstrap_files:
                        search_excludes = ['resource_uri'] 
                    logger.info(str(('+++++++++++processing file', filename)))
                    with open(filename) as data_file:
                        input_data = serializer.from_csv(data_file.read())
                        
                        if command == 'put':
                            resp = testApiClient.put(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
                            logger.debug(str(('action: ', json.dumps(action), 
                                              'response: ' , resp.status_code)))
                            self.assertTrue(
                                resp.status_code in [200], str((resp.status_code, resp)))
                            
                            # now see if we can get these objects back
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), data={ 'limit': 999 })
                            self.assertTrue(
                                resp.status_code in [200], str((resp.status_code, resp)))
                            #   self.assertValidJSONResponse(resp)
                            new_obj = self.deserialize(resp)
                            result, msgs = find_all_obj_in_list(
                                input_data['objects'], new_obj['objects'], 
                                excludes=search_excludes)
                            self.assertTrue(
                                result, str((command, 'input file', filename, 
                                             msgs, new_obj['objects'])) )
                        
                        elif command == 'patch':
                            resp = testApiClient.patch(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
#                             self.assertHttpAccepted(resp)
                            self.assertTrue(resp.status_code in [202, 204], str((
                                'response not accepted, resource_uri:', resource_uri, 
                                'response', resp)))
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), data={ 'limit': 999 } )
                            self.assertTrue(
                                resp.status_code in [200], str((resp.status_code, resp)))
                            #                             self.assertValidJSONResponse(resp)
                            new_obj = self.deserialize(resp)
                            with open(filename) as f2:
                                input_data2 = serializer.from_csv(f2.read())
                                result, msgs = find_all_obj_in_list(
                                    input_data2['objects'], new_obj['objects'], 
                                    excludes=search_excludes)
                                self.assertTrue(
                                    result, str(( command, 'input file', filename, msgs )) )
                        
                        elif command == 'post':
                            self.fail((
                                'resource entry: ' + json.dumps(action) + '; '
                                'cannot POST multiple objects to tastypie; '
                                'therefore the "post" command is invalid with '
                                'the initialization scripts'))
                        else:
                            self.fail('unknown command: ' + command + ', ' + json.dumps(action))


class LibraryContentLoadTest(MetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        print '============== LibraryContentLoadTest setup ============'
        super(LibraryContentLoadTest, self).setUp()
        super(LibraryContentLoadTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(LibraryContentLoadTest, self)._bootstrap_init_files()
        print '============== LibraryContentLoadTest setup: begin ============'
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
    def test_load_sdf(self):
        pass;
        

class LibraryTest(MetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        print '============== LibraryTest setup ============'
        super(LibraryTest, self).setUp()
        super(LibraryTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(LibraryTest, self)._bootstrap_init_files()
        print '============== LibraryTest setup: begin ============'
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
        testApiClient = TestApiClient(serializer=self.csv_serializer) 

        filename = os.path.join(self.db_directory,'metahash_fields_library.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.library'})

        print '============== LibraryTest setup: done ============'
        

    def test1_create_library(self):
        logger.info(str(('==== test1_create_library =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        
        logger.info(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        # create a second library
        library_item = LibraryFactory.attributes()
        
        logger.info(str(('item', library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(library_item, new_obj['objects'])
        self.assertTrue(
            result, str(('bootstrap item not found', obj, 
                         library_item, new_obj['objects'])))
        logger.info(str(('item found', obj)))


    def test2_create_library_invalid_library_type(self):
        logger.info(str(('==== test2_create_library_invalid_library_type =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        library_item['library_type'] = 'invalid_type'
        
        logger.info(str(('try to create an invalid library_type:', library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        
        from reports.dump_obj import dumpObj
        logger.info(str(('response', dumpObj(resp))))
        
        self.assertTrue(resp.status_code in [400], str((resp.status_code, resp)))
        
        logger.info(str(('response.content.library message', 
                         getattr(resp, 'content'))))
        
        obj = json.loads(getattr(resp, 'content'))
        logger.info(str(('dump', dumpObj(obj))))

class ScreensTest(MetaHashResourceBootstrap,ResourceTestCase):
    
    def setUp(self):
        print '============== ScreensTest setup ============'
        super(ScreensTest, self).setUp()
        super(ScreensTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(ScreensTest, self)._bootstrap_init_files()
        print '============== ScreensTest setup: begin ============'
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
        testApiClient = TestApiClient(serializer=self.csv_serializer) 

        filename = os.path.join(self.db_directory,'metahash_fields_screen.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.screen'})

        print '============== ScreensTest setup: done ============'
        

    def test1_create_screen(self):
        logger.info(str(('==== test1_create_screen =====')))
        
        resource_uri = BASE_URI_DB + '/screen'
        
        screen_item = ScreenFactory.attributes()
        
        logger.info(str(('screen_item',screen_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        screen_item = ScreenFactory.attributes()
        
        logger.info(str(('item', screen_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(screen_item, new_obj['objects'])
        self.assertTrue(
            result, str(('bootstrap item not found', obj, 
                         screen_item, new_obj['objects'])))
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
            resp = self.api_client.post(
                resource_uri, format='json', data=item, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
#             self.assertHttpCreated(resp)
            
        logger.info('created items, now get them')
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        for i,item in enumerate(self.bootstrap_items):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(
                result, str(('bootstrap item not found', item, new_obj['objects'])))
            self.assertTrue(
                'facility_id' in obj, 'the facility_id was not created')
            self.assertEqual(
                str((i+1)), obj['facility_id'], 
                str(('expected the facility_id returned to be incremented to ',
                     i+1, obj)) )
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
        
