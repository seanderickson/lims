from django.test import TestCase
import json
import factory
import os
import logging
from tastypie.test import ResourceTestCase, TestApiClient

import db.models
from test.factories import *
from reports.tests import assert_obj1_to_obj2, find_all_obj_in_list, find_obj_in_list
from reports.serializers import CSVSerializer
from reports.tests import MetaHashResourceBootstrap
import reports.tests


logger = logging.getLogger(__name__)


BASE_URI = '/db/api/v1'
# import reports.BASE_URI as BASE_REPORTS_URI
BASE_REPORTS_URI = '/reports/api/v1'
import db; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__))
BASE_URI_DB = '/db/api/v1'


# class TestApiInit(reports.tests.TestApiInit):
#     
#     def setUp(self):
#         super(TestApiInit, self).setUp()
# #         super(TestApiInit, self)._setUp()
# #         super(TestApiInit, self).test_0bootstrap_metahash()
#         # NOTE: run reports tests to set up the metahash resource
#         
#         self.db_resource_uri = BASE_URI + '/metahash'
#         self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
# 
#         
#     def test_2api_init(self):
#         
#         logger.debug('***================ super: db test_2api_init =============== ')
# 
#         super(TestApiInit, self).test_2api_init()
#         
#         logger.debug('***================ local: db test_2api_init =============== ')
#         
#         
#         serializer=CSVSerializer() 
#         # todo: doesn't work for post, see TestApiClient.post() method, it is 
#         # incorrectly "serializing" the data before posting
#         testApiClient = TestApiClient(serializer=serializer) 
#         
#         filename = os.path.join(self.db_directory,'api_init_actions.csv')
#         with open(filename) as input_file:
#             api_init_actions = serializer.from_csv(input_file.read(), root=None)
#             
#             bootstrap_files = [
#                 'metahash_resource_data.csv',
#                 'vocabularies_data.csv']
#             for action in api_init_actions:
#                 
#                 logger.debug('\n++++=========== processing action', json.dumps(action))
#                 command = action['command'].lower() 
#                 resource = action['resource'].lower()
#                 resource_uri = BASE_REPORTS_URI + '/' + resource
#                 
#                 if command == 'delete':
#                     resp = testApiClient.delete(
#                         resource_uri, authentication=self.get_credentials())
#                     self.assertHttpAccepted(resp)
#                 
#                 else:
#                     filename = os.path.join(self.db_directory,action['file'])
#                     search_excludes = []
#                     # exclude 'resource_uri' from equivalency check during 
#                     # bootstrap, because resource table needs to be loaded for
#                     # the uri generation
#                     if action['file'] in bootstrap_files:
#                         search_excludes = ['resource_uri'] 
#                     logger.debug(str(('+++++++++++processing file', filename)))
#                     with open(filename) as data_file:
#                         input_data = serializer.from_csv(data_file.read())
#                         
#                         if command == 'put':
#                             resp = testApiClient.put(
#                                 resource_uri, format='csv', data=input_data, 
#                                 authentication=self.get_credentials() )
#                             logger.debug(str(('action: ', json.dumps(action), 
#                                               'response: ' , resp.status_code)))
#                             self.assertTrue(
#                                 resp.status_code in [200], str((resp.status_code, resp)))
#                             
#                             # now see if we can get these objects back
#                             resp = testApiClient.get(
#                                 resource_uri, format='json', 
#                                 authentication=self.get_credentials(), data={ 'limit': 999 })
#                             self.assertTrue(
#                                 resp.status_code in [200], str((resp.status_code, resp)))
#                             #   self.assertValidJSONResponse(resp)
#                             new_obj = self.deserialize(resp)
#                             result, msgs = find_all_obj_in_list(
#                                 input_data['objects'], new_obj['objects'], 
#                                 excludes=search_excludes)
#                             self.assertTrue(
#                                 result, str((command, 'input file', filename, 
#                                              msgs, new_obj['objects'])) )
#                         
#                         elif command == 'patch':
#                             resp = testApiClient.patch(
#                                 resource_uri, format='csv', data=input_data, 
#                                 authentication=self.get_credentials() )
# #                             self.assertHttpAccepted(resp)
#                             self.assertTrue(resp.status_code in [202, 204], str((
#                                 'response not accepted, resource_uri:', resource_uri, 
#                                 'response', resp)))
#                             resp = testApiClient.get(
#                                 resource_uri, format='json', 
#                                 authentication=self.get_credentials(), data={ 'limit': 999 } )
#                             self.assertTrue(
#                                 resp.status_code in [200], str((resp.status_code, resp)))
#                             #                             self.assertValidJSONResponse(resp)
#                             new_obj = self.deserialize(resp)
#                             with open(filename) as f2:
#                                 input_data2 = serializer.from_csv(f2.read())
#                                 result, msgs = find_all_obj_in_list(
#                                     input_data2['objects'], new_obj['objects'], 
#                                     excludes=search_excludes)
#                                 self.assertTrue(
#                                     result, str(( command, 'input file', filename, msgs )) )
#                         
#                         elif command == 'post':
#                             self.fail((
#                                 'resource entry: ' + json.dumps(action) + '; '
#                                 'cannot POST multiple objects to tastypie; '
#                                 'therefore the "post" command is invalid with '
#                                 'the initialization scripts'))
#                         else:
#                             self.fail('unknown command: ' + command + ', ' + json.dumps(action))


class LibraryContentLoadTest(MetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        logger.debug('============== LibraryContentLoadTest setup ============')
        super(LibraryContentLoadTest, self).setUp()
        super(LibraryContentLoadTest, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(LibraryContentLoadTest, self)._bootstrap_init_files()
        logger.debug('============== LibraryContentLoadTest setup: begin ============')
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
    def test_load_sdf(self):
        pass;

class DBMetaHashResourceBootstrap(MetaHashResourceBootstrap):

    def _bootstrap_init_files(self):
        '''
        test loads the essential files of the api initialization, the 'bootstrap':
        - PUT metahash_fields_initial.csv
        - PATCH metahash_fields_initial_patch.csv
        - PATCH metahash_fields_resource.csv
        - PATCH metahash_vocabularies.csv
        - PUT vocabularies_data.csv
        - PUT metahash_resource_data.csv
        '''
        logger.debug('------------- DBMetaHashResourceBootstrap _bootstrap_init_files -----------------')        
        super(DBMetaHashResourceBootstrap, self)._bootstrap_init_files()

        
        serializer=CSVSerializer() 
        resource_uri = BASE_URI + '/metahash'
        directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')

        # Note, once the resources are loaded, can start checking the 
        # resource_uri that is returned
        filename = os.path.join(directory, 'metahash_resource_data.csv')
        (input, output) = self._patch_test('resource', filename, id_keys_to_check=['key'])
        filename = os.path.join(directory, 'vocabularies_data.csv')
        self._patch_test('vocabularies', filename, id_keys_to_check=['key'])

        
        logger.debug('------------- Done: DBMetaHashResourceBootstrap _bootstrap_init_files -----------------')        
        
        

class LibraryResource(DBMetaHashResourceBootstrap,ResourceTestCase):

    def setUp(self):
        logger.debug('============== LibraryResource setup ============')
        super(LibraryResource, self).setUp()
        super(LibraryResource, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(LibraryResource, self)._bootstrap_init_files()
        logger.debug('============== LibraryResource setup: begin ============')
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
        testApiClient = TestApiClient(serializer=reports.serializers.LimsSerializer) 

        filename = os.path.join(self.db_directory,'metahash_fields_library.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.library'})

        logger.debug('============== LibraryResource setup: done ============')
        

    def test1_create_library(self):
        logger.debug(str(('==== test1_create_library =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        # create a second library
        library_item = LibraryFactory.attributes()
        
        logger.debug(str(('item', library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(library_item, new_obj['objects'])
        self.assertTrue(
            result, str(('bootstrap item not found', obj, 
                         library_item, new_obj['objects'])))
        logger.debug(str(('item found', obj)))


    def test2_create_library_invalid_library_type(self):
        logger.debug(str(('==== test2_create_library_invalid_library_type =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        library_item['library_type'] = 'invalid_type'
        
        logger.debug(str(('try to create an invalid library_type:', library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        
        self.assertTrue(resp.status_code in [400], str((resp.status_code, resp)))
        
        logger.debug(str(('response.content.library message', 
                         getattr(resp, 'content'))))
        
        obj = json.loads(getattr(resp, 'content'))
        logger.debug(str(('content', obj)))    
        self.assertTrue('library_type' in obj['library'], 
            str(('content should have an entry for the faulty "name"', obj)))
        logger.debug(str(('==== done: test2_create_library_invalid_library_type =====')))

    def test3_create_library_invalid_library_name(self):
        logger.debug(str(('==== test3_create_library_invalid_library_name =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        library_item['library_name'] = 'invalid & name'
        
        logger.debug(str(('try to create an invalid library name:', library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        
        self.assertTrue(resp.status_code in [400], str((resp.status_code, resp)))
        
        logger.debug(str(('response.content.library message', 
                         getattr(resp, 'content'))))
        
        obj = json.loads(getattr(resp, 'content'))
        logger.debug(str(('content', obj)))    
        self.assertTrue('library_name' in obj['library'], 
            str(('content should have an entry for the faulty "name"', obj)))
        
        # TODO: test the error message
        
        logger.debug(str(('==== done: test3_create_library_invalid_library_name =====')))

    def test4_create_library_invalids(self):
        '''
        Test the schema "required" validations 
        
        - todo: this can be a template for other endpoints
        ''' 
        logger.debug(str(('==== test4_create_library_invalids =====')))

        resource_uri = BASE_REPORTS_URI + '/resource/library'
        logger.debug(str(('Get the library schema', resource_uri )))
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        self.assertTrue(resp.status_code in [200], 
                        str((resp.status_code, resp.serialize())))
        new_obj = self.deserialize(resp)
        fields = new_obj['schema']['fields']
        logger.debug(str(('=== field keys', fields.keys())))
        resource_uri = BASE_URI_DB + '/library'
        
        # make sure the default works
        library_item = LibraryFactory.attributes()
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        for key,field in fields.items():
            logger.debug(str(('key, field', key,field)))
            if field.get('required', False):
                logger.debug(str(('testing required field', key, field)))
                
                library_item = LibraryFactory.attributes()
                library_item[key] = None
                resp = self.api_client.post(
                    resource_uri, format='json', data=library_item, 
                    authentication=self.get_credentials())
                self.assertTrue(resp.status_code in [400], str((resp.status_code, resp)))
                logger.debug(str(('response.content.library message', 
                                 getattr(resp, 'content'))))
        
                obj = json.loads(getattr(resp, 'content'))
                logger.debug(str(('==========content', obj)))
                self.assertTrue(key in obj['library'], 
                    str(('error content should have an entry for the faulty key:',key, obj)))
                
                
        # TODO: test regex and number: min/max
        
        
        logger.debug(str(('==== done: test4_create_library_invalids =====')))

class ScreenResource(DBMetaHashResourceBootstrap,ResourceTestCase):
    
        
    def setUp(self):
        logger.debug('============== ScreenResource setup ============')
        super(ScreenResource, self).setUp()
        super(ScreenResource, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(ScreenResource, self)._bootstrap_init_files()
        logger.debug('============== ScreenResource setup: begin ============')
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
#         testApiClient = TestApiClient(serializer=self.csv_serializer) 
        testApiClient = TestApiClient(serializer=reports.serializers.LimsSerializer) 

        filename = os.path.join(self.db_directory,'metahash_fields_screen.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.screen'})

        logger.debug('============== ScreenResource setup: done ============')
        

    def test1_create_screen(self):
        logger.debug(str(('==== test1_create_screen =====')))
        
        resource_uri = BASE_URI_DB + '/screen'
        
        screen_item = ScreenFactory.attributes()
        
        logger.debug(str(('screen_item',screen_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        screen_item = ScreenFactory.attributes()
        
        logger.debug(str(('item', screen_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(screen_item, new_obj['objects'])
        self.assertTrue(
            result, str(('bootstrap item not found', obj, 
                         screen_item, new_obj['objects'])))
        self.assertTrue('facility_id' in obj, 'the facility_id was not created')
        logger.debug(str(('item found', obj)))

    def test0_create_screen(self):
        logger.debug(str(('==== test_create_screen =====')))
        
        # the simplest of tests, create some simple screens
        self.bootstrap_items = [   
            {
                'facility_id': "1",
                'project_phase': "primary_screen",
                'screen_type': "small_molecule",
                'title': "Test screen 1",
                'data_sharing_level': 1,
                'total_plated_lab_cherry_picks': 0,
            },
            {
                'facility_id': "2",
                'project_phase': "primary_screen",
                'screen_type': "rnai",
                'title': "Test Screen 2",
                'data_sharing_level': 0,
                'total_plated_lab_cherry_picks': 0,
            },
        ]
        
        resource_uri = BASE_URI_DB + '/screen'
        logger.debug(str(('--posting to:', resource_uri)))
        for i,item in enumerate(self.bootstrap_items):
            logger.debug(str(('item', item)))         
            resp = self.api_client.post(
                resource_uri, format='json', data=item, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
#             self.assertHttpCreated(resp)
            
        logger.debug('created items, now get them')
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
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
            logger.debug(str(('item found', obj)))
            assert_obj1_to_obj2(item, obj)

        logger.debug(str(('==== test_create_screen done =====')))
        

class MutualScreensTest(DBMetaHashResourceBootstrap,ResourceTestCase):
    
    def test_mutual_positives_to_screen(self):
        pass
        # create two screens
        # create two screen results
        # set data sharing levels
        # find assay wells that overlap
        
