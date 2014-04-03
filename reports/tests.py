
"""
About testing:

From: https://docs.djangoproject.com/en/dev/topics/testing/overview/s
The test database
Tests that require a database (namely, model tests) will not use your "real"
(production) database. Separate, blank databases are created for the tests.

Regardless of whether the tests pass or fail, the test databases are destroyed 
when all the tests have been executed.

By default the test databases get their names by prepending test_ to the value 
of the NAME settings for the databases defined in DATABASES. When using the 
SQLite database engine the tests will by default use an in-memory database 
(i.e., the database will be created in memory, bypassing the filesystem entirely!). 
If you want to use a different database name, specify TEST_NAME in the dictionary 
for any given database in DATABASES.

Aside from using a separate database, the test runner will otherwise use all of
the same database settings you have in your settings file: ENGINE, USER, HOST, 
etc. The test database is created by the user specified by USER, so you'll need
to make sure that the given user account has sufficient privileges to create a
new database on the system.

NOTE: to run using sqlite in memory, use a test_settings file with the database as:
DATABASES['default'] = {'ENGINE': 'django.db.backends.sqlite3'}
and run like:
$ ./manage.py test --settings=lims.test_settings

"""

import os
import logging
import json

from django.test import TestCase
from django.contrib.auth.models import User
from tastypie.test import ResourceTestCase, TestApiClient
from reports.models import Vocabularies
from tastypie.fields import BooleanField

from lims.api import CSVSerializer
from lims.api import CsvBooleanField
from lims.tests import assert_obj1_to_obj2, find_all_obj_in_list, find_obj_in_list


logger = logging.getLogger(__name__)

BASE_URI = '/reports/api/v1'
import reports; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(reports.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(reports.__path__))
    
class SerializerTest(TestCase):

    def test_csv(self):
        logger.info(str(('======== test_csv =========')))
        directory = APP_ROOT_DIR
        serializer = CSVSerializer() 
        
        input = [['one','two', 'three', 'four', 'five'],
                ['uno', '2', 'true', 'false', '' ]]
        
        
        input_data = serializer.from_csv_iterate(input, root=None)
        for obj in input_data:
            self.assertTrue(obj['one']=='uno')
            self.assertTrue(obj['two']=='2')
            self.assertTrue(obj['three']=='true')
            self.assertTrue(obj['four']=='false')
            self.assertTrue(obj['five']=='')
        logger.info(str(('input to object', input_data)))
        csv_data = serializer.to_csv(input_data, root=None)
        logger.info(str(('back to csv', csv_data)))

        with open(directory + '/reports/test/test_csv_.csv', 'w') as f:
            f.write(csv_data)
            
        with open(directory + '/reports/test/test_csv_.csv') as fin:    
            final_data = serializer.from_csv(fin.read(), root=None)
            logger.info(str(('final_data', final_data)))
            for obj in final_data:
                self.assertTrue(obj['one']=='uno')
                self.assertTrue(obj['two']=='2')
                self.assertTrue(obj['three']=='true')
                self.assertTrue(obj['four']=='false')
                self.assertTrue(obj['five']=='')
        
        # TODO: delete the file
        logger.info(str(('======== test_csv done =========')))

class HydrationTest(TestCase):
    
    def test_hydrate_boolean_tastypie(self):
        logger.info(str(('======== test_hydrate_boolean =========')))
        
        field = BooleanField()
        
        test_data = [ True, False, 'true', 'TRUE', 'FALSE', 'false', '' ]
        expected  = [ True, False, True, True, True, True, False ]
        
        for i, item in enumerate(test_data):
            result = field.convert(item)
            self.assertEqual(result, expected[i], 
                             str((i,' is not equal', item, result, expected[i])))
        logger.info(str(('======== test_hydrate_boolean done =========')))
            
    def test_hydrate_boolean_lims(self):
        field = CsvBooleanField()
        
        test_data = [ True, False, 'true', 'TRUE', 'FALSE', 'false', '' ]
        expected  = [ True, False, True, True, False, False, False ]
        
        for i, item in enumerate(test_data):
            result = field.convert(item)
            self.assertEqual(result, expected[i], 
                             str((i,' is not equal', item, result, expected[i])))
  
    
class MetaHashResourceBootstrap():
    
    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)

    def _setUp(self):
        # Create a user.
        self.username = 'daniel'
        self.password = 'pass'
        self.user = User.objects.create_user(
            self.username, 'daniel@example.com', self.password)
        
        self.resource_uri = BASE_URI + '/metahash'
        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        self.csv_serializer=CSVSerializer() 
        # todo: doesn't work for post, see TestApiClient.post() method,
        # it is incorrectly "serializing" the data before posting
        self.testApiClient = TestApiClient(serializer=self.csv_serializer) 
   
    def _patch_test(self,resource_name, filename, keys_not_to_check=[], 
                    id_keys_to_check=[], data_for_get={}):
        '''
        data_for_get - dict of extra header information to send with the GET request
        '''
        data_for_get.setdefault('limit', 999 )
        resource_uri = BASE_URI + '/' + resource_name
        
        logger.info(str(('===resource_uri', resource_uri)))
        with open(filename) as bootstrap_file:
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())

            logger.info(str(('Submitting patch...', bootstrap_file)))
            resp = self.testApiClient.patch(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials() )
            logger.info(str(('Response: ' , resp.status_code)))
#            self.assertHttpAccepted(resp)
            self.assertTrue(resp.status_code in [202, 204], str((resp)))
            
            logger.info(str(('check patched data for',resource_name,
                             'execute get on:',resource_uri)))
            resp = self.api_client.get(
                resource_uri, format='json', authentication=self.get_credentials(), 
                data=data_for_get)
            logger.info(str(('--------resp to get:', resp.status_code)))
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
            new_obj = self.deserialize(resp)
            
            for inputobj in input_data['objects']:
                result, outputobj = find_obj_in_list(
                    inputobj,new_obj['objects'], excludes=keys_not_to_check )
                self.assertTrue(
                    result, 
                    str(('not found', outputobj,'=== objects returned ===', 
                         new_obj['objects'] )) ) 
                self.assertTrue(
                    resource_name in outputobj['resource_uri'], 
                    str(('wrong resource_uri returned:', outputobj,
                         'should contain', resource_name)))
                for id_key in id_keys_to_check:
                    self.assertTrue(
                        inputobj[id_key] in outputobj['resource_uri'], 
                        str(('wrong resource_uri returned:', outputobj,
                             'should contain id key', id_key, 'val', inputobj[id_key])))

            # return both collections for further testing
            return (input_data['objects'], new_obj['objects']) 

    def _put_test(self, resource_name, filename, keys_not_to_check=[], 
                  id_keys_to_check=[], data_for_get={}):
        '''
        id_keys_to_check if the resource data has been loaded, 
            then these are id keys to check to see if they are being used in 
            the returned resource_uri field
        '''
        data_for_get.setdefault('limit', 999 )
        resource_uri = BASE_URI + '/' + resource_name

        with open(filename) as bootstrap_file:
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())
            
            logger.info(str(('Submitting put...', bootstrap_file)))
            resp = self.testApiClient.put(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials() )
            logger.info(str(('Response: ' , resp.status_code)))
#            self.assertHttpAccepted(resp)
            self.assertTrue(resp.status_code in [200, 202, 204], str((resp)))
    
            logger.info(str(('check put data for',resource_name,
                             'execute get on:',resource_uri)))
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get)
            logger.info(str(('--------resp to get:', resp.status_code)))
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
            new_obj = self.deserialize(resp)
            # do a length check, since put will delete existing resources
            self.assertEqual(len(new_obj['objects']), len(input_data['objects']), 
                             'input length != output length: ' + str((new_obj)))
            
            for inputobj in input_data['objects']:
                result, outputobj = find_obj_in_list(
                    inputobj,new_obj['objects'], excludes=keys_not_to_check)
                self.assertTrue(result, str(('not found', outputobj, 
                                             new_obj['objects'] )) )
                self.assertTrue(resource_name in outputobj['resource_uri'], 
                    str(('wrong resource_uri returned:', outputobj,
                         'should contain', resource_name)))
                for id_key in id_keys_to_check:
                    self.assertTrue(inputobj[id_key] in outputobj['resource_uri'], 
                        str(('wrong resource_uri returned:', outputobj,
                             'should contain id key', id_key, 'val', inputobj[id_key])))
                
            # return both collections for further testing
            return (input_data['objects'], new_obj['objects']) 
                    
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
        print '------------- _bootstrap_init_files -----------------'
        serializer=CSVSerializer() 
        resource_uri = BASE_URI + '/metahash'
        # todo: doesn't work for post, see TestApiClient.post() method, it is 
        # incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=serializer) 

        filename = os.path.join(self.directory, 'metahash_fields_initial.csv')
        self._put_test('metahash', filename, keys_not_to_check=['resource_uri'])

        filename = os.path.join(self.directory, 'metahash_fields_initial_patch.csv')
        self._patch_test('metahash', filename, keys_not_to_check=['resource_uri'])

        filename = os.path.join(self.directory, 'metahash_fields_resource.csv')
        self._patch_test('metahash', filename, keys_not_to_check=['resource_uri'], 
                         data_for_get={ 'scope':'fields.resource' })
                        
        filename = os.path.join(self.directory,'metahash_fields_vocabularies.csv')
        self._patch_test('metahash', filename, keys_not_to_check=['resource_uri'], 
                         data_for_get={ 'scope':'fields.vocabularies' })

        # Note, once the resources are loaded, can start checking the 
        # resource_uri that is returned
        filename = os.path.join(self.directory, 'metahash_resource_data.csv')
        (input, output) = self._put_test('resource', filename, id_keys_to_check=['key'])
        filename = os.path.join(self.directory, 'vocabularies_data.csv')
        self._put_test('vocabularies', filename, id_keys_to_check=['key'])

  
        print '------------- done _bootstrap_init_files -----------------'

        
class TestApiInit(MetaHashResourceBootstrap,ResourceTestCase):
    
    def setUp(self):
        super(TestApiInit, self).setUp()
        super(TestApiInit, self)._setUp()

    def test_0bootstrap_metahash(self):
        
        print '================ test_0bootstrap_metahash =============== '
        # in order for the metahash resource to work, the metahash itself must 
        # be "bootstrapped":
        # the metahash must be filled with the fields that describe itself
        # these are the "bootstrap" fields
        bootstrap_items = [   
            {
                'key': 'scope',
                'scope': 'fields.metahash',
                'ordinal': 0    
            },
            {
                'key': 'key',
                'scope': 'fields.metahash',
                'ordinal': 1   
            },
            {
                'key': 'ordinal',
                'scope': 'fields.metahash',
                'ordinal': 2   
            },
            {
                'key': 'json_field_type',
                'scope': 'fields.metahash',
                'ordinal': 3    
            }
        ]
        
        for item in bootstrap_items:         
            resp = self.api_client.post(
                self.resource_uri, format='json', data=item, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
            #             self.assertHttpCreated(resp)
            
        logger.info('created items, now get them')
        resp = self.api_client.get(
            self.resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp, resp.status_code)))
        new_obj = self.deserialize(resp)
        logger.info(str(('deserialized object:', json.dumps(new_obj))))
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        self.assertEqual(len(new_obj['objects']), 4, str((new_obj)))
        
        for inputobj in bootstrap_items:
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, str(('not found', outputobj )) )
        
        print '================ test_0bootstrap_metahash done ========== '
    
    def test_1bootstrap_init(self):
        print '================ test_1bootstrap_init ================'
        self._bootstrap_init_files()
        print '================ test_1bootstrap_init done ================'
        
    def test_2api_init(self):
        
        print '***================ test_2api_init =============== '
        serializer=CSVSerializer() 
        # todo: doesn't work for post, see TestApiClient.post() method, it is 
        # incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=serializer) 
        
        filename = os.path.join(self.directory,'api_init_actions.csv')
        with open(filename) as input_file:
            api_init_actions = serializer.from_csv(input_file.read(), root=None)
            
            bootstrap_files = [
                'metahash_fields_initial.csv',
                'metahash_fields_initial_patch.csv',
                'metahash_fields_resource.csv',
                'metahash_resource_data.csv',
                'metahash_fields_vocabularies.csv',
                'vocabularies_data.csv']
            for action in api_init_actions:
                
                print '\n++++=========== processing action', json.dumps(action)
                command = action['command'].lower() 
                resource = action['resource'].lower()
                resource_uri = BASE_URI + '/' + resource
                
                if command == 'delete':
                    resp = testApiClient.delete(
                        resource_uri, authentication=self.get_credentials())
                    self.assertHttpAccepted(resp)
                
                else:
                    filename = os.path.join(self.directory,action['file'])
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
#                             self.assertHttpAccepted(resp)
                            
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
                            self.assertHttpAccepted(resp)
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
                    

class UserResource(MetaHashResourceBootstrap,ResourceTestCase):
    
    def setUp(self):
        print '============== User setup ============'
        super(UserResource, self).setUp()
        super(UserResource, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(UserResource, self)._bootstrap_init_files()
        print '============== User setup: begin ============'
        
        # Create the User resource field entries
        testApiClient = TestApiClient(serializer=self.csv_serializer) 

        filename = os.path.join(self.directory,'metahash_fields_user.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.user'})

        print '============== User setup: done ============'
    
    def _test0_create_user(self):
        logger.info(str(('==== test_create_user =====')))
        
        # the simplest of tests, create some simple users
        self.bootstrap_items = [   
            {
                'ecommons_id': 'st1',
                'first_name': 'Sally',
                'last_name': 'Tester', 
                'email': 'sally.tester@limstest.com',    
            },
            {
                'login_id': 'jt1',
                'first_name': 'Joe',
                'last_name': 'Tester',    
                'email': 'joe.tester@limstest.com',    
            },
            {
                'login_id': 'bt1',
                'first_name': 'Bad',
                'last_name': 'TestsALot',    
                'email': 'bad.fester@slimstest.com',    
            },
        ]

        for i,item in enumerate(self.bootstrap_items):         
            resp = self.api_client.post(self.resource_uri, 
                format='json', data=item, authentication=self.get_credentials())
            self.assertHttpCreated(resp)
            
        logger.info('created items, now get them')
        resp = self.api_client.get(self.resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        self.assertEqual(len(new_obj['objects']), 3, str((new_obj)))
        
        for i,item in enumerate(self.bootstrap_items):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(
                result, str(('bootstrap item not found', item, new_obj['objects'])))
            self.assertTrue(
                'screensaver_user_id' in obj, 'the primary key was not created')
            self.assertEqual(
                i+1, obj['screensaver_user_id'], 
                str(('expected the screensaver_user_id returned to be incremented to ', i+1, obj)) )
            logger.info(str(('item found', obj)))

        logger.info(str(('==== test_create_user done =====')))

    def _test1_user_test_data(self):
        logger.info(str(('==== test1_user_test_data =====')))
        
        filename = os.path.join(self.directory,'test_data/users1.csv')
        self._put_test('user', filename)
        
        logger.info(str(('==== test1_user_test_data done =====')))
    
    def test2_user_permissions(self):
        logger.info(str(('==== test2_user_permissions =====')))
        self._test1_user_test_data()
        
        logger.info(str(('==== test2_user_permissions start =====')))
        filename = os.path.join(self.directory,'test_data/users2_patch.csv')
        self._patch_test('user', filename)
        
        logger.info(str(('==== test2_user_permissions done =====')))

class UserGroupResource(UserResource):
    
    def setUp(self):
        logger.info(str(( '============== UserGroup setup ============')))
        super(UserGroupResource, self).setUp()
        
        logger.info(str(( '============== UserGroup setup: begin ============')))
        meta_resource_uri = BASE_URI + '/metahash'
        
        # Create the User resource field entries
        # todo: doesn't work for post, see TestApiClient.post() method, 
        # it is incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=self.csv_serializer) 
        
        filename = os.path.join(self.directory,'metahash_fields_usergroup.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.usergroup'})
        
        filename = os.path.join(self.directory,'metahash_fields_permission.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.permission'})

        logger.info(str(( '============== UserGroup setup done ============')))

    def _test2_usergroup(self):
        logger.info(str(('==== test2_usergroup =====')))
        #create users
        self._test1_user_test_data()
        
        logger.info(str(('----- test2_usergroup =====')))

        filename = os.path.join(self.directory,'test_data/usergroups1.csv')
        self._put_test('usergroup', filename)


        logger.info(str(('==== test2_usergroup done =====')))


    def test3_user_groups(self):
        logger.info(str(('==== test3_user_groups =====')))
        self._test2_usergroup()
        logger.info(str(('==== test3_user_groups start =====')))
        filename = os.path.join(self.directory,'test_data/users3_groups_patch.csv')
        self._patch_test('user', filename)
      
        logger.info(str(('==== test3_user_groups done =====')))
