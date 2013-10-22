
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

from django.test import TestCase

import datetime
import csv
import StringIO
from django.contrib.auth.models import User
from django.test.client import Client
from tastypie.test import ResourceTestCase, TestApiClient
from reports.models import MetaHash, Vocabularies
from lims.api import CSVSerializer

import logging
import json
from tastypie.fields import BooleanField
from reports.api import CsvBooleanField

logger = logging.getLogger(__name__)

BASE_URI = '/reports/api/v1'

def assert_obj1_to_obj2(obj1, obj2, keys=['key', 'scope', 'ordinal', 'username', 'name'], excludes=['resource_uri']):
    original_keys = set(obj1.keys())
    updated_keys = set(obj2.keys())
        
    intersect_keys = original_keys.intersection(updated_keys)
    if intersect_keys != original_keys:
        return False, ('keys missing', original_keys-intersect_keys)
#    logger.info(str(('intersect_keys', intersect_keys)))
    for key in keys:
        if key not in obj1:
            continue
        if key not in obj2:
            return False, None
        if str(obj1[key]) != str(obj2[key]):
            return False, None
    keys_to_search = set(obj1.keys()) - set(keys) - set(excludes)
#    logger.info(str(('equal so far, keys to search', keys_to_search)))
    csvBooleanField = CsvBooleanField()
    
    for key in keys_to_search:
        if isinstance(obj1[key], basestring):
            val1 = str(obj1[key])
            if val1 != str(obj2[key]):
    #            logger.info('not equal')
                if csvBooleanField.convert(val1) != csvBooleanField.convert(obj1[key]):
                    return False, ('key not equal', key, obj1, obj2, 'val obj1', str(obj1[key]), 'val 2', str(obj2[key]))
    #            logger.info(str((key, 'boolean equal')))
    #        logger.info(str((key, 'equal', obj1, obj2)))
        else: # better be a list
            for v in obj1[key]:
                found = False
                for v2 in obj2[key]:
                    if str(v2) == str(v): # have to do this because return values are being unicoded
                        found = True
                if not found:
                    return False, ('list key not equal', key, obj1, obj2, 'val obj1', str(obj1[key]), 'val 2', str(obj2[key]))

    return True, (key, obj1, obj2)


def find_obj_in_list(obj, item_list):
    list_msgs = []
    for item in item_list:
        result, msgs = assert_obj1_to_obj2(obj, item)
        if result:
            logger.debug(str(('found', obj)))
            return True, (item)
        else:
            list_msgs.append(msgs)
#        elif msgs:
#            return result, msgs
    return False, ('obj not found in list', obj, list_msgs)

def find_all_obj_in_list(list1, list2):
    for item in list1:
        result, msgs = find_obj_in_list(item, list2)
        if not result:
            logger.debug(str(('-----not found', item, list2, msgs)))
            return False, msgs
    return True, msgs


class BaselineTest(TestCase):
    
    def test_assert_obj1_to_obj2(self):
        # all of obj1 in obj2 (may be more in obj2
        
        obj1 = { 'one': 1, 'two': 'two', 'three':''}
        obj2 = { 'one': 1, 'two': 'two', 'three':'', 'four': 'blah'}
        result, msgs = assert_obj1_to_obj2(obj1, obj2)
        self.assertTrue(result, msgs)
        result, msgs = assert_obj1_to_obj2(obj2, obj1)
        self.assertFalse(result, msgs)

    def test_find_obj_in_list(self):
        
        obj1 = { 'one': 1, 'two': 'two', 'three':''}
        obj1a = { 'one': 2, 'two': 'two', 'three':''}
        obj2 = { 'one': 1, 'two': 'two', 'three':'', 'four': 'blah'}
        obj3 = { 'xxx': 1, 'two': 'two', 'three':'', 'four': 'blah'}
        
        obj_list = [ obj2, obj1, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertTrue(result, msgs)

        obj_list = [ obj2, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertTrue(result, msgs)

        obj_list = [ obj3, obj1a ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        logger.info(str((result,msgs)))
        self.assertFalse(result, str((msgs)))

        obj_list = [ obj1a, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertFalse(result, msgs)

    def test_find_all_obj_in_list(self):
        
        obj1 = { 'one': 1, 'two': 'two', 'three':''}
        obj1a = { 'one': 2, 'two': 'two', 'three':''}
        obj2 = { 'one': 1, 'two': 'two', 'three':'', 'four': 'blah'}
        obj3 = { 'xxx': 1, 'two': 'two', 'three':'', 'four': 'blah'}
        
        obj_list = [ obj1a, obj1 ]
        obj_list2 = [ obj2, obj1a, obj3 ]
        result, msgs = find_all_obj_in_list(obj_list, obj_list2)
        self.assertTrue(result, msgs)

        obj_list2 = [ obj2, obj3 ]
        result, msgs = find_all_obj_in_list(obj_list, obj_list2)
        self.assertFalse(result, msgs)


class SerializerTest(TestCase):

    def test_csv(self):
        directory = './reports/static/api_init'
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

        with open(directory + '/test_csv_.csv', 'w') as f:
            f.write(csv_data)
            
        with open(directory + '/test_csv_.csv') as fin:    
            final_data = serializer.from_csv(fin.read(), root=None)
            logger.info(str(('final_data', final_data)))
            for obj in final_data:
                self.assertTrue(obj['one']=='uno')
                self.assertTrue(obj['two']=='2')
                self.assertTrue(obj['three']=='true')
                self.assertTrue(obj['four']=='false')
                self.assertTrue(obj['five']=='')


class HydrationTest(TestCase):
    
    def test_hydrate_boolean_tastypie(self):
        
        field = BooleanField()
        
        test_data = [ True, False, 'true', 'TRUE', 'FALSE', 'false', '' ]
        expected  = [ True, False, True, True, True, True, False ]
        
        for i, item in enumerate(test_data):
            result = field.convert(item)
            self.assertEqual(result, expected[i], str((i,' is not equal', item, result, expected[i])))
            
    def test_hydrate_boolean_lims(self):
        field = CsvBooleanField()
        
        test_data = [ True, False, 'true', 'TRUE', 'FALSE', 'false', '' ]
        expected  = [ True, False, True, True, False, False ]
        
        for i, item in enumerate(test_data):
            result = field.convert(item)
            self.assertEqual(result, expected[i], str((i,' is not equal', item, result, expected[i])))
            
    
class MetaHashTest(ResourceTestCase):

    def test_get_and_parse(self):
        print '================ test_get_and_parse ========== '
        initializer = {
                       'key': 'key',
                       'scope': 'fields:metahash',
                       'ordinal': 0    }
        MetaHash.objects.create(**initializer)
        hash = MetaHash.objects.get_and_parse(scope='fields:metahash')
        self.assertTrue(len(hash)==1, str((hash)))
        print '================ test_get_and_parse done ========== '


class MetaHashResourceBootstrap(ResourceTestCase):
    
    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)

    def setUp(self):
        super(MetaHashResourceBootstrap, self).setUp()

        # Create a user.
        self.username = 'daniel'
        self.password = 'pass'
        self.user = User.objects.create_user(self.username, 'daniel@example.com', self.password)
        
        self.resource_uri = BASE_URI + '/metahash'
                
    def test_bootstrap_metahash(self):
        
        print '================ test_bootstrap_metahash =============== '
        # in order for the metahash resource to work, the metahash itself must be "bootstrapped":
        # the metahash must be filled with the fields that describe itself
        # these are the "bootstrap" fields
        bootstrap_items = [   
            {
                'key': 'scope',
                'scope': 'fields:metahash',
                'ordinal': 0    
            },
            {
                'key': 'key',
                'scope': 'fields:metahash',
                'ordinal': 1   
            },
            {
                'key': 'ordinal',
                'scope': 'fields:metahash',
                'ordinal': 2   
            },
            {
                'key': 'json_field_type',
                'scope': 'fields:metahash',
                'ordinal': 3    
            }
        ]
        
        for item in bootstrap_items:         
            resp = self.api_client.post(self.resource_uri, format='json', data=item, authentication=self.get_credentials())
            print resp
            self.assertHttpCreated(resp)
            
        logger.info('created items, now get them')
        resp = self.api_client.get(self.resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp, resp.status_code)))
        new_obj = self.deserialize(resp)
        logger.info(str(('deserialized object:', json.dumps(new_obj))))
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 4, str((new_obj)))
        
        for item in bootstrap_items:
            self.assertTrue(find_obj_in_list(item, new_obj['objects']), str(('bootstrap item not found', item, new_obj['objects'])))
        
        print '================ test_bootstrap_metahash done ========== '
        
    def test_bootstrap_init_files(self):
        
        print '================ test_bootstrap_init_files =============== '
        serializer=CSVSerializer() 
        directory = '/home/sde4/workspace/iccbl-lims/'
        resource_uri = BASE_URI + '/metahash'
        testApiClient = TestApiClient(serializer=serializer) # todo: doesn't work for post, see TestApiClient.post() method, it is incorrectly "serializing" the data before posting
        with open(directory + 'reports/static/api_init/metahash_fields_initial.csv') as bootstrap_file:
            input_data = serializer.from_csv(bootstrap_file.read())

            resp = testApiClient.put(resource_uri, format='csv', data=input_data, authentication=self.get_credentials() )
            print 'Response: ' , resp.status_code
            self.assertHttpAccepted(resp)
    
            logger.info('===== get1:'+ resource_uri)
            resp = self.api_client.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
            logger.info(str(('--------resp to get:', resp, resp.status_code)))
            new_obj = self.deserialize(resp)
            self.assertValidJSONResponse(resp)
            self.assertEqual(len(new_obj['objects']), len(input_data['objects']), 'input length != output length: ' + str((new_obj)))
            
            for inputobj in input_data['objects']:
                self.assertTrue(find_obj_in_list(inputobj,new_obj['objects']), str(('not found', inputobj, new_obj['objects'])) )


        with open(directory + 'reports/static/api_init/metahash_fields_resource1.csv') as bootstrap_file1:
            input_data = serializer.from_csv(bootstrap_file1.read())
            print 'Submitting...'
            resp = testApiClient.patch(resource_uri, format='csv', data=input_data, authentication=self.get_credentials() )
            print 'Response: ' , resp.status_code
            self.assertHttpAccepted(resp)
        
            logger.info('===== get1:'+ resource_uri)
            resp = self.api_client.get(resource_uri, format='json', authentication=self.get_credentials())
            logger.info(str(('--------resp to get:', resp, resp.status_code)))
            new_obj = self.deserialize(resp)
            self.assertValidJSONResponse(resp)
#            self.assertEqual(len(new_obj['objects']), len(input_data['objects']), 'input length != output length: ' + str((new_obj)))
#            self.assertEqual(len(new_obj['objects']), len(input_data['objects']), 'input length != output length: ' + str((new_obj)))
            
            for inputobj in input_data['objects']:
                self.assertTrue(find_obj_in_list(inputobj,new_obj['objects']), str(('not found', inputobj, new_obj['objects'])) )
            
#        todo: look at TestCase to find out how to reference files as  "fixtures"
#        see: https://docs.djangoproject.com/en/dev/ref/django-admin/#loaddata-fixture-fixture
        print '================ done test_bootstrap_init_files =============== '

        
    def test_api_init(self):
        
        print '***================ test_api_init =============== '
        serializer=CSVSerializer() 
        testApiClient = TestApiClient(serializer=serializer) # todo: doesn't work for post, see TestApiClient.post() method, it is incorrectly "serializing" the data before posting
        directory = '/home/sde4/workspace/iccbl-lims/reports/static/api_init/'
        
        with open(directory + '/api_init_actions.csv') as input_file:
            api_init_actions = serializer.from_csv(input_file.read(), root=None)
            
            for action in api_init_actions:
                
                print '\n++++=========== processing action', json.dumps(action)
                command = action['command'].lower() 
                resource = action['resource'].lower()
                resource_uri = BASE_URI + '/' + resource
                
                if command == 'delete':
                    resp = testApiClient.delete(resource_uri, authentication=self.get_credentials())
                    self.assertHttpAccepted(resp)
                
                else:
                    input_file = action['file']
                    logger.info(str(('+++++++++++processing file', input_file)))
                    with open(directory + '/' + input_file ) as data_file:
                        input_data = serializer.from_csv(data_file.read())
                        
                        if command == 'put':
                            resp = testApiClient.put(resource_uri, format='csv', data=input_data, authentication=self.get_credentials() )
                            logger.debug(str(('action: ', json.dumps(action), 'response: ' , resp.status_code)))
                            self.assertHttpAccepted(resp)
                            
                            # now see if we can get these objects back
                            resp = testApiClient.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
                            self.assertValidJSONResponse(resp)
                            new_obj = self.deserialize(resp)
                            result, msgs = find_all_obj_in_list(input_data['objects'], new_obj['objects'])
                            self.assertTrue(result, str(( command, 'input file', input_file, msgs, new_obj['objects'])) )
                        
                        elif command == 'patch':
                            resp = testApiClient.patch(resource_uri, format='csv', data=input_data, authentication=self.get_credentials() )
#                            print 'action: ', json.dumps(action), 'response: ' , resp.status_code, resp
                            self.assertHttpAccepted(resp)
                            resp = testApiClient.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 } )
#                            print '----- response', resp
                            self.assertValidJSONResponse(resp)
                            new_obj = self.deserialize(resp)
                            with open(directory + '/' + input_file) as f2:
                                input_data2 = serializer.from_csv(f2.read())
                                print 'test input against returned data'
                                result, msgs = find_all_obj_in_list(input_data2['objects'], new_obj['objects'])
                                self.assertTrue(result, str(( command, 'input file', input_file, msgs )) )
                        
                        elif command == 'post':
                            self.fail('resource entry: ' + json.dumps(action) 
                                + '; cannot POST multiple objects to tastypie; therefore the "post" command is invalid with the initialization scripts')
                        else:
                            self.fail('unknown command: ' + command + ', ' + json.dumps(action))
        
class MetaResourceBase(ResourceTestCase):

    def setUp(self):
        super(MetaResourceBase, self).setUp()

        # Create a user.
        self.username = 'daniel'
        self.password = 'pass'
        self.user = User.objects.create_user(self.username, 'daniel@example.com', self.password)
        
        # "Bootstrap" the metahash table.  (See the bootstrapping test above)
        bootstrap_items = [   
            {
                'key': 'scope',
                'scope': 'fields:metahash',
                'ordinal': 0    
            },
            {
                'key': 'key',
                'scope': 'fields:metahash',
                'ordinal': 1   
            },
            {
                'key': 'ordinal',
                'scope': 'fields:metahash',
                'ordinal': 2   
            },
            {
                'key': 'json_field_type',
                'scope': 'fields:metahash',
                'ordinal': 3    
            }
        ]
        for item in bootstrap_items:
            MetaHash.objects.create(**item)

    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)


class ResourceResourceTest(MetaResourceBase):
    def setUp(self):
        
        logger.info(str(('---setUp----')))
        super(ResourceResourceTest, self).setUp()

        self.resource_uri = BASE_URI + '/resource'

        bootstrap_items = [   
            {
                'key': 'key',
                'scope': 'fields:resource',
                'ordinal': 0    
            },
            {
                'key': 'scope',
                'scope': 'fields:resource',
                'ordinal': 1    
            },
            {
                'key': 'ordinal',
                'scope': 'fields:resource',
                'ordinal': 2    
            },
            {
                'key': 'id_attribute',
                'scope': 'fields:resource',
                'ordinal': 3,
                'json_field_type': 'fields.ListField'    
            }
        ]
        for item in bootstrap_items:
            MetaHash.objects.create(**item)
        logger.info(str(('test resource fields created')))   
        for field in MetaHash.objects.all():
            print '===== metahash field', field

        
    def test_resource_create(self):
        self._do_test_resource_create()
        
    def _do_test_resource_create(self):
        
        bootstrap_data = { 'objects': [
            {
                'key': 'resource',
                'scope': 'resource',
                'id_attribute': ['scope','key'],
                'ordinal': 0
            },
            {
                'key': 'usergroup',
                'scope': 'resource',
                'id_attribute': ['name'],
                'ordinal': 1
            },
            {
                'key': 'user',
                'scope': 'resource',
                'id_attribute': ['username'],
                'ordinal': 2
            }
        ]}
        resource_resource_uri = BASE_URI + '/resource'
        logger.info('===== post1:'+ resource_resource_uri )
        resp = self.api_client.put(resource_resource_uri, format='json', data=bootstrap_data, authentication=self.get_credentials())
        logger.info(str(('--------resp to post:', resp, resp.status_code)))
        new_obj = self.deserialize(resp)
        logger.info(str(('deserialized object:', json.dumps(new_obj))))
        
        for oldobj in bootstrap_data['objects']:
            self.assertTrue(find_obj_in_list(oldobj,new_obj['objects'] ), str(('not found', oldobj, '   new obj:', new_obj['objects'])))
            
        
class UserGroupResource(ResourceResourceTest):
    
    def xsetUp(self):
        
        logger.info(str(('---setUp----')))
        super(UserGroupResource, self).setUp()

        self.resource_uri = BASE_URI + '/usergroup'
        
        # create some more users 
        self.superusername = 'superuser'
        password = 'pass'
        self.superuser = User.objects.create_superuser(self.superusername, 'superuser@example.com', password)
        
        self.testuser1 = User.objects.create_user('testuser1', 'testuser1@example.com', 'pass')
        self.testuser2 = User.objects.create_user('testuser2', 'testuser2@example.com', 'pass')
        
        logger.info(str(('test users created', self.testuser1, self.testuser2)))
       # "Bootstrap" the usergroup table.  (See the bootstrapping test above)
        bootstrap_items = [   
            {
                'key': 'name',
                'scope': 'fields:usergroup',
                'ordinal': 0    
            },
            {
                'key': 'users',
                'scope': 'fields:usergroup',
                'ordinal': 1   
            },
            {
                'key': 'permissions',
                'scope': 'fields:usergroup',
                'ordinal': 2   
            },
            {
                'key': 'user_list',
                'scope': 'fields:usergroup',
                'ordinal': 3    
            }
        ]
        for item in bootstrap_items:
            MetaHash.objects.create(**item)
        logger.info(str(('test fields created')))   

        self._do_test_resource_create();
        
        for field in MetaHash.objects.all().filter(scope='fields:resource'):
            print '===== metahash field', field


             
    def xtest_create_group(self):
        
        logger.info(str(('============= test_create_group========', self.testuser1 )))
        initializer = {
            'name': 'test_group1',
            'users': [ BASE_URI + '/user/' + self.testuser1.username ]
            }

        logger.info('===== post1:'+ self.resource_uri)
        resp = self.api_client.post(self.resource_uri, format='json', data=initializer, authentication=self.get_credentials())
        logger.info(str(('--------resp to post:', resp, resp.status_code)))
        new_obj = self.deserialize(resp)
        logger.info(str(('deserialized object:', json.dumps(new_obj))))
        self.assertHttpCreated(resp)       
        
        logger.info('===== get1:'+ self.resource_uri)
        resp = self.api_client.get(self.resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.info(str(('--------resp to get:', resp, resp.status_code)))
        new_obj = self.deserialize(resp)
        logger.info(str(('deserialized object:', json.dumps(new_obj))))
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 1, str((new_obj)))
        self.assertTrue(find_obj_in_list(initializer,new_obj['objects']), str(('not found', initializer, new_obj['objects'])))
        
        logger.info(str(('============= test_create_group done ========' )))
        
#
#class MetaHashResourceUsage(MetaResourceBase):
#    # Use ``fixtures`` & ``urls`` as normal. See Django's ``TestCase``
#    # documentation for the gory details.
#    fixtures = ['test_entries.json']
#
#    def setUp(self):
#        super(MetaHashResourceUsage, self).setUp()
#    
#    def test_get_list_unauthorzied(self):
#        self.assertHttpUnauthorized(self.api_client.get('/reports/api/v1/metahash/', format='json'))
#
#    def test_get_list_json(self):
#        print '=============== test get_list ==================='
#
#        resp = self.api_client.get('/reports/api/v1/metahash/', format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
#        print resp
#        self.assertValidJSONResponse(resp)
#        
#    def test_put_json(self):
#        print '=============== test put_json ==================='
#        resource_uri = '/reports/api/v1/metahash'
#        
#        # first get, see what's there
#        logger.info('=====get prelim:'+ resource_uri)
#        resp = self.api_client.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
#        logger.info(str(('--------resp to get:', resp, resp.status_code)))
#        new_obj = self.deserialize(resp)
#        
#        extant_object_count = len(new_obj['objects'])
#        
#        # first create a "json" field by specifying the "json_type" of the field 
#        # (if this weren't specified, it would throw an error, since there is no "real" field for this (todo: test))
#        initializer = {
#               'key': 'test_field',
#               'scope': 'fields:metahash',
#               'ordinal': 0, 
#               'json_field_type': 'fields.CharField'    }
#
#        logger.info('===== post1:'+ resource_uri)
#        resp = self.api_client.post(resource_uri, format='json', data=initializer, authentication=self.get_credentials())
#        logger.info(str(('--------resp to post:', resp, resp.status_code)))
#        new_obj = self.deserialize(resp)
#        logger.info(str(('deserialized object:', json.dumps(new_obj))))
#        self.assertHttpCreated(resp)
#
#        logger.info('===== get1:'+ resource_uri)
#        resp = self.api_client.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
#        logger.info(str(('--------resp to get:', resp, resp.status_code)))
#        new_obj = self.deserialize(resp)
#        logger.info(str(('deserialized object:', json.dumps(new_obj))))
#        self.assertValidJSONResponse(resp)
#        self.assertEqual(len(new_obj['objects']), extant_object_count+1, str((new_obj)))
#        self.assertTrue(find_obj_in_list(initializer,new_obj['objects'] ), str(('not found', initializer, new_obj['objects'])))
#        
#        # now update with some field instance data
#        initializer = {
#               'key': 'test_field',
#               'scope': 'fields:metahash',
#               'ordinal': 0, 
#               'json_field_type': 'fields.CharField', 
#               'test_field': 'foo and bar!'    }
#        uri = resource_uri + '/' + initializer['scope']+ '/' + initializer['key'] + '/'
#        
#        logger.info('===== put1:'+ uri)
#        resp = self.api_client.put(uri, data=initializer, authentication=self.get_credentials())
#        logger.info(str(('--------resp to put:', resp, resp.status_code)))
#        new_obj = self.deserialize(resp)
#        logger.info(str(('deserialized object:', json.dumps(new_obj))))
#        self.assertHttpAccepted(resp)
#        self.assertEqual(new_obj['test_field'], 'foo and bar!', 'unexpected result: ' +new_obj['test_field'])
#        self.assertTrue(assert_obj1_to_obj2(initializer, new_obj), str(('assert_obj1_to_obj2', initializer, new_obj)))
#
#        logger.info('===== get final:'+ resource_uri)
#        resp = self.api_client.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
#        logger.info(str(('--------resp to get:', resp, resp.status_code)))
#        new_obj = self.deserialize(resp)
#        logger.info(str(('deserialized object:', json.dumps(new_obj))))
#        self.assertValidJSONResponse(resp)
#        self.assertEqual(len(new_obj['objects']), extant_object_count+1, str((new_obj)))
#
#        print '===================  finished: test put_json modification ==================='
#        
#    # this should work, but the tastypie.ResourceTestCase methods around post for my custom serializer are broken
#    # todo: just test the api directly through the requests api
#    def xxxx_test_csv_serialization(self):
#        print '=============== test_csv_serialization1 ==================='
#        ''' 
#        Note: we are only verifying simple property lists, i.e. key-value lists; 
#        where values must be: either a string or a list of strings
#        '''
#        serializer=CSVSerializer() 
#        self.api_client = TestApiClient(serializer=serializer) # todo: doesn't work for post, see TestApiClient.post() method, it is incorrectly "serializing" the data before posting
#        posting_client = Client()
#        
#        header = ['key', 'scope', 'ordinal', 'json_field_type']
#        vals = ['test_field', 'fields:metahash', 0, 'fields.CharField']
#        
#        raw_data = StringIO.StringIO()
#        writer = csv.writer(raw_data)
#        writer.writerow(header)
#        writer.writerow(vals)
#        
#        kwargs = {  'content_type': 'text/csv',
#                    'HTTP_AUTHORIZATION': self.get_credentials()
#                }
#
#        resp = posting_client.post('/reports/api/v1/metahash/', data=raw_data.getvalue(), **kwargs)
#        new_obj = self.deserialize(resp)
#        print '----- resp:' , resp, '---', new_obj
#        logger.info(str(( '----- resp:' , resp, '---', new_obj)))
#        self.assertHttpCreated(resp)
#
#        print '-------------------  test put_json modification ==================='
#        # now create an some field instance data
#        header = ['key', 'scope', 'test_field']
#        vals = ['test_field', 'fields:metahash', 'foo and bar!']
#        
#        raw_data = StringIO.StringIO()
#        writer = csv.writer(raw_data)
#        writer.writerow(header)
#        writer.writerow(vals)
#
#        uri = '/reports/api/v1/metahash/'+ str(new_obj.get('id')) + '/'
#        print 'uri:', uri
#        resp = self.api_client.put(uri, format='csv', serializer=serializer, data=raw_data.getvalue(), authentication=self.get_credentials())
#        print '--------resp to put:', resp
#        self.assertHttpAccepted(resp)
#        new_obj = self.deserialize(resp)
#        self.assertEqual(new_obj['test_field'], 'foo and bar!', 'unexpected result: ' +new_obj['test_field'])
#
#        resp = self.api_client.get('/reports/api/v1/metahash/', format='json', authentication=self.get_credentials())
#        self.assertValidJSONResponse(resp)
#        self.assertEqual(len(self.deserialize(resp)['objects']), 2)
#        logger.info(str(('---- resp:' , resp)))
        
        
        
#        # Scope out the data for correctness.
#        self.assertEqual(len(self.deserialize(resp)['objects']), 12)
#        # Here, we're checking an entire structure for the expected data.
#        self.assertEqual(self.deserialize(resp)['objects'][0], {
#            'pk': str(self.entry_1.pk),
#            'user': '/api/v1/user/{0}/'.format(self.user.pk),
#            'title': 'First post',
#            'slug': 'first-post',
#            'created': '2012-05-01T19:13:42',
#            'resource_uri': '/api/v1/entry/{0}/'.format(self.entry_1.pk)
#        })

#    def test_get_list_xml(self):
#        self.assertValidXMLResponse(self.api_client.get('/api/v1/entries/', format='xml', authentication=self.get_credentials()))

#    def test_get_detail_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.get(self.detail_url, format='json'))
#
#    def test_get_detail_json(self):
#        resp = self.api_client.get(self.detail_url, format='json', authentication=self.get_credentials())
#        self.assertValidJSONResponse(resp)
#
#        # We use ``assertKeys`` here to just verify the keys, not all the data.
#        self.assertKeys(self.deserialize(resp), ['created', 'slug', 'title', 'user'])
#        self.assertEqual(self.deserialize(resp)['name'], 'First post')
#
##    def test_get_detail_xml(self):
##        self.assertValidXMLResponse(self.api_client.get(self.detail_url, format='xml', authentication=self.get_credentials()))
#
#    def test_post_list_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.post('/api/v1/entries/', format='json', data=self.post_data))
#
#    def test_post_list(self):
#        # Check how many are there first.
#        self.assertEqual(Entry.objects.count(), 5)
#        self.assertHttpCreated(self.api_client.post('/api/v1/entries/', format='json', data=self.post_data, authentication=self.get_credentials()))
#        # Verify a new one has been added.
#        self.assertEqual(Entry.objects.count(), 6)
#
#    def test_put_detail_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.put(self.detail_url, format='json', data={}))
#
#    def test_put_detail(self):
#        # Grab the current data & modify it slightly.
#        original_data = self.deserialize(self.api_client.get(self.detail_url, format='json', authentication=self.get_credentials()))
#        new_data = original_data.copy()
#        new_data['title'] = 'Updated: First Post'
#        new_data['created'] = '2012-05-01T20:06:12'
#
#        self.assertEqual(Entry.objects.count(), 5)
#        self.assertHttpAccepted(self.api_client.put(self.detail_url, format='json', data=new_data, authentication=self.get_credentials()))
#        # Make sure the count hasn't changed & we did an update.
#        self.assertEqual(Entry.objects.count(), 5)
#        # Check for updated data.
#        self.assertEqual(Entry.objects.get(pk=25).title, 'Updated: First Post')
#        self.assertEqual(Entry.objects.get(pk=25).slug, 'first-post')
#        self.assertEqual(Entry.objects.get(pk=25).created, datetime.datetime(2012, 3, 1, 13, 6, 12))
#
#    def test_delete_detail_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.delete(self.detail_url, format='json'))
#
#    def test_delete_detail(self):
#        self.assertEqual(Entry.objects.count(), 5)
#        self.assertHttpAccepted(self.api_client.delete(self.detail_url, format='json', authentication=self.get_credentials()))
#        self.assertEqual(Entry.objects.count(), 4)