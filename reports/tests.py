
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
# from reports.models import Vocabularies
from tastypie.fields import BooleanField

from reports.serializers import CsvBooleanField,CSVSerializer, SDFSerializer
from reports.utils.sdf2py import MOLDATAKEY
from tastypie.utils.dict import dict_strip_unicode_keys


logger = logging.getLogger(__name__)

BASE_URI = '/reports/api/v1'
import reports; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(reports.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(reports.__path__))


def is_boolean(field):
    if type(field) == bool: return True
    
    if isinstance(field, basestring):
        val = field.lower()
        if val=='true' or val=='false':
            return True
    return False
    
csvBooleanField = CsvBooleanField()

####
# NOTE: equivocal, and other equivalency methods are for end-to-end testing 
# of the API with values read from csv or json.
# Values are submitted, then read back as json from the API.
####

def equivocal(val1, val2):
    '''
    For testing equality from the CSV input to the JSON response of what is 
    created in the DB
    obj1 input
    obj2 output 
    '''
    if val1 == val2:
        return True, ('val1', val1, 'val2', val2 )
    
    if not val1:
        if val2:
            if val2 or len(val2) > 0:
                return False, ('val1', val1, 'val2', val2 )
        return True, ('val1', val1, 'val2', val2 )

    if isinstance(val1, (int, long, float, complex)):
        val1 = str(val1)
        val2 = str(val2)
        if val1 != val2:
            return False, ('val1', val1, 'val2', val2 )
        
    
    if isinstance(val1, basestring):
        val1 = str(val1)
        val2 = str(val2)
        if val1 != val2:
            if (is_boolean(val1) and 
                    csvBooleanField.convert(val1)==csvBooleanField.convert(val2)):
                return True, ('val1', val1, 'val2', val2 )
            return False, ('val1', val1, 'val2', val2 )
    else: # better be a list
        assert isinstance(val1, list) and isinstance(val2, list), \
               str(('Must be a list if not a string', 'val1', val1, 'val2', val2))
        for v in val1:
            if not v: 
                if val2 and len(val2) > 0:
                    return False, ('val1', val1, 'list val2 not empty', val2 )
            if v not in val2:
                if v not in [str(v2) for v2 in val2]:
                    return False, ('val1', val1, 'val2', val2 )
    return True, ('val1', val1, 'val2', val2 )
    
#NOTE only works for String comparisons
def assert_obj1_to_obj2(
        obj1, obj2, 
        keys=[],
#         keys=['key', 'scope', 'ordinal', 'username', 'name'], 
        excludes=['resource_uri']):
    '''
    For testing equality from the CSV input to the JSON response of what is 
    created in the DB
    obj1 input
    obj2 output 
    '''
    original_keys = set(obj1.keys())
    updated_keys = set(obj2.keys())
    
    intersect_keys = original_keys.intersection(updated_keys)
    if intersect_keys != original_keys:
        return False, ('keys missing', original_keys-intersect_keys)
    for key in keys:
        if key not in obj1:
            continue
        if key not in obj2:
            return False, ('key not found',key)  
        result, msgs =  equivocal(obj1[key], obj2[key])
        if not result:
            # don't report for this section, just move along to the next items to test
            return False, () 
        #        if str(obj1[key]) != str(obj2[key]):
        #            return False, None  # don't report for this section
        #    logger.debug(str(('equal so far, keys to search', keys_to_search)))

    fkey = 'resource_uri'
    if fkey not in excludes:
        excludes.append(fkey)
        if fkey in obj1:
            if fkey not in obj2:
                return False, ('no resource uri in obj2', obj1, obj2)
            val1 = str(obj1[fkey])
            val2 = str(obj2[fkey])
            
            if val1 != val2:
                if val1.endswith('/'): val1 = val1[0:len(val1)-1]
                if val2.endswith('/'): val2 = val2[0:len(val2)-1]
                if val1 != val2:
                    if val1 not in val2:
                        return False, ('resource uri', val1, val2)
                    logger.warn(str((
                        'imprecise uri matching, equivocating:', val1, val2)))
    
    keys_to_search = set(obj1.keys()) - set(keys) - set(excludes)

    csvBooleanField = CsvBooleanField()
    
    for key in keys_to_search:
        result, msgs =  equivocal(obj1[key], obj2[key])
        if not result:
            return False, ('key not equal', key, 'obj1', obj1, 'obj2', obj2, 'msgs', msgs)
            
    return True, ('obj1:', obj1, 'obj2:', obj2)


def find_obj_in_list(obj, item_list, **kwargs):
    list_msgs = []
    for item in item_list:
        result, msgs = assert_obj1_to_obj2(obj, item, **kwargs)
        if result:
            logger.debug(str(('found', obj, item)))
            return True, (item)
        else:
            if not msgs in list_msgs:
                list_msgs.append(msgs)
    return False, ('obj not found in list', obj, list_msgs)

def find_all_obj_in_list(list1, list2, **kwargs):
    msgs = ['not run yet']
    for item in list1:
        result, msgs = find_obj_in_list(item, list2, **kwargs)
        if not result:
            logger.debug(str(('-----not found', item, list2, msgs)))
            return False, msgs
    return True, msgs

class SDFSerializerTest(TestCase):
    
    def test1_read(self):
        logger.debug('=== test1 SDF read')
        
        records = [{
            'smiles': 'Cl',
            'chemical_name': 'HCl (hydrochloric acid)',
            'molecular_mass': '35.976677742',
            MOLDATAKEY: '''
  SciTegic12121315112D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
M  END
'''
            },
            {
            'smiles': 'C([C@@H](C(=O)O)O)C(=O)O',
            'inchi': 'InChI=1S/C4H6O5/c5-2(4(8)9)1-3(6)7/h2,5H,1H2,(H,6,7)(H,8,9)/t2-/m0/s1',
            'chemical_name': 'malate',
            'molecular_mass': '134.021523292',
            MOLDATAKEY: '''
  SciTegic12121315112D

  9  8  0  0  1  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  3  5  1  0  0  0  0
  2  6  1  1  0  0  0
  1  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
M  END
'''                
            },{
            'smiles': 'O.O.Cl',
            'inchi': 'InChI=1S/ClH.2H2O/h1H;2*1H2',
            'chemical_name': 'hydrochloride hydrate',
            'molecular_mass': '71.99780711',
            MOLDATAKEY: '''
  SciTegic12121315112D

  3  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
M  END
'''
            }]
        
        serializer = SDFSerializer()
        
        sdf_data = serializer.to_sdf(records)
        with open(APP_ROOT_DIR + '/lims/test/test1.sdf', 'w') as fout:    
            fout.write(sdf_data)
            
        fout.close()
    
        with open(APP_ROOT_DIR + '/lims/test/test1.sdf') as fin:    
            _data = serializer.from_sdf(fin.read(), root=None)
            final_data = _data['objects']
            logger.debug(str(('final_data', final_data)))
            
            self.assertTrue(
                len(final_data)==len(records), 
                str(('len is', len(final_data),len(records))))
            for obj in final_data:
                logger.debug(str(('object: ', obj)))
                
                for record in records:
                    if record['smiles'] == obj['smiles']:
                        for k,v in record.items():
                            self.assertTrue(k in obj, str(('obj does not contain',k)))
                            self.assertTrue(
                                obj[k] == v, 
                                str(('values not equal', k, record[k], '\n read: ', v)))

class BaselineTest(TestCase):

    def test0_equivocal(self):
        logger.debug(str(('>> test0_equivocal <<')))
        
        test_true = [
            ['',''],
            ['', None],
            ['1','1'],
            ['1',1], # integer strings can come back as integers (if integerfield is defined)
            ['sal','sal'],
            ['True','TRUE'],
            ['true','True'],
            ['TRUE',True], # boolean strings can come back as booleans (if booleanfield is defined)
            ['false','False'],
            ['FALSE','false'],
            ['False',False],
            ]
            ## TODO: test date strings 
        
        for [a,b] in test_true:
            logger.debug(str(('csv serialization equivocal truthy test', a, b)))
            result, msgs = equivocal(a,b)
            self.assertTrue(result, msgs)

        test_false = [
            ['1',''],
            ['1',2],    # although string integers may be converted, still want them to be equivocal
            ['','1'],
            ['sal','val'],
            ['True','1'],
            ['true',''],
            ['TRUE',False],
            ['false','true'],
            ['False',True],
            ]
        for [a,b] in test_false:
            logger.debug(str(('csv serialization equivocal falsy test', a, b)))
            result, msgs = equivocal(a,b)
            self.assertFalse(result, msgs)
        logger.debug(str(('>> test0_equivocal done <<')))
    
    def test1_assert_obj1_to_obj2(self):
        logger.debug(str(('====== test1_assert_obj1_to_obj2 ====')))
        # all of obj1 in obj2 (may be more in obj2
        
        obj1 = { 'one': '1', 'two': 'two', 'three':''}
        obj2 = { 'one': '1', 'two': 'two', 'three':'', 'four': 'blah'}
        result, msgs = assert_obj1_to_obj2(obj1, obj2)
        self.assertTrue(result, msgs)
        result, msgs = assert_obj1_to_obj2(obj2, obj1)
        self.assertFalse(result, msgs)
        obj3 = { 'one':'1', 'two':'2', 'three':'' }
        result, msgs = assert_obj1_to_obj2(obj3, obj1)
        self.assertFalse(result, str((result,msgs)))

        obj4 = { 'one':'1', 'two':'2', 'three':'' }
        obj5 = { 'one':'1', 'two':'2', 'three':None }
        result, msgs = assert_obj1_to_obj2(obj5, obj4)
        self.assertTrue(result, str((result,msgs)))
        logger.debug(str(('====== test1_assert_obj1_to_obj2 done ====')))

    def test2_find_obj_in_list(self):
        logger.debug(str(('====== test2_find_obj_in_list ====')))
        
        obj1 = { 'one': '1', 'two': 'two', 'three':''}
        obj1a = { 'one': '2', 'two': 'two', 'three':''}
        obj2 = { 'one': '1', 'two': 'two', 'three':'', 'four': 'blah'}
        obj3 = { 'xxx': '1', 'two': 'two', 'three':'', 'four': 'blah'}
        
        obj_list = [ obj2, obj1, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertTrue(result, msgs)

        obj_list = [ obj2, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertTrue(result, msgs)

        obj_list = [ obj3, obj1a ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        logger.debug(str((result,msgs)))
        self.assertFalse(result, str((msgs)))

        obj_list = [ obj1a, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertFalse(result, msgs)
        logger.debug(str(('====== test2_find_obj_in_list done ====')))

    def test3_find_all_obj_in_list(self):
        logger.debug(str(('====== test3_find_all_obj_in_list ====')))
        
        obj1 = { 'one': '1', 'two': 'two', 'three':''}
        obj1a = { 'one': '2', 'two': 'two', 'three':''}
        obj2 = { 'one': '1', 'two': 'two', 'three':'', 'four': 'blah'}
        obj3 = { 'xxx': '1', 'two': 'two', 'three':'', 'four': 'blah'}
        
        obj_list = [ obj1a, obj1 ]
        obj_list2 = [ obj2, obj1a, obj3 ]
        result, msgs = find_all_obj_in_list(obj_list, obj_list2)
        self.assertTrue(result, msgs)

        obj_list2 = [ obj2, obj3 ]
        result, msgs = find_all_obj_in_list(obj_list, obj_list2)
        self.assertFalse(result, str((result, msgs  )))
        logger.debug(str(('====== test3_find_all_obj_in_list done ====')))

    def test4_user_example(self):
        logger.debug(str(('====== test4_user_example ====')))
        bootstrap_items = [   
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
        item = {'login_id': 'jt1', 'first_name': 'Joe', 'last_name': 'Tester', 
                'email': 'joe.tester@limstest.com'}
        result, obj = find_obj_in_list(item, bootstrap_items)
        logger.debug(str((result, obj)))
        
        self.assertTrue(obj['email'] == 'joe.tester@limstest.com')
        
        logger.debug(str(('====== test4_user_example done ===='))) 
    
    
import reports.utils.serialize
    
class SerializerTest(TestCase):

    def test_csv(self):
        logger.debug(str(('======== test_csv =========')))
        directory = APP_ROOT_DIR
        serializer = CSVSerializer() 
        
        input = [['one','two', 'three', 'four', 'five'],
                ['uno', '2', 'true', 'false', '' ]]
        
        
        input_data = reports.utils.serialize.from_csv_iterate(input)
        for obj in input_data:
            self.assertTrue(obj['one']=='uno')
            self.assertTrue(obj['two']=='2')
            self.assertTrue(obj['three']=='true')
            self.assertTrue(obj['four']=='false')
            self.assertTrue(obj['five']=='')
        logger.debug(str(('input to object', input_data)))
        csv_data = serializer.to_csv(input_data, root=None)
        logger.debug(str(('back to csv', csv_data)))

        with open(directory + '/reports/test/test_csv_.csv', 'w') as f:
            f.write(csv_data)
            
        with open(directory + '/reports/test/test_csv_.csv') as fin:    
            final_data = serializer.from_csv(fin.read(), root=None)
            logger.debug(str(('final_data', final_data)))
            for obj in final_data:
                self.assertTrue(obj['one']=='uno')
                self.assertTrue(obj['two']=='2')
                self.assertTrue(obj['three']=='true')
                self.assertTrue(obj['four']=='false')
                self.assertTrue(obj['five']=='')
        
        # TODO: delete the file
        logger.debug(str(('======== test_csv done =========')))

class HydrationTest(TestCase):
    
    def test_hydrate_boolean_tastypie(self):
        logger.debug(str(('======== test_hydrate_boolean =========')))
        
        field = BooleanField()
        
        test_data = [ True, False, 'true', 'TRUE', 'FALSE', 'false', '' ]
        expected  = [ True, False, True, True, True, True, False ]
        
        for i, item in enumerate(test_data):
            result = field.convert(item)
            self.assertEqual(result, expected[i], 
                             str((i,' is not equal', item, result, expected[i])))
        logger.debug(str(('======== test_hydrate_boolean done =========')))
            
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
        self.username = 'testsuper'
        self.password = 'pass'
        self.user = User.objects.create_superuser(
            self.username, 'testsuperuser@example.com', self.password)
        
        self.resource_uri = BASE_URI + '/metahash'
        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        self.csv_serializer=CSVSerializer() 
        # todo: doesn't work for post, see TestApiClient.post() method,
        # it is incorrectly "serializing" the data before posting
        self.testApiClient = TestApiClient(serializer=self.csv_serializer) 
   
    def _patch_test(self,resource_name, filename, keys_not_to_check=[], 
                    id_keys_to_check=[], data_for_get={}):
        '''
        data_for_get - dict of extra header krmation to send with the GET request
        '''
        data_for_get.setdefault('limit', 999 )
        data_for_get.setdefault('HTTP_APILOG_COMMENT', 'patch_test: %s' % filename )
        resource_uri = BASE_URI + '/' + resource_name
        
        logger.debug(str(('===resource_uri', resource_uri)))
        with open(filename) as bootstrap_file:
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())

            logger.debug(str(('Submitting patch...', bootstrap_file)))
            resp = self.testApiClient.patch(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            logger.debug(str(('Response: ' , resp.status_code)))
            #            self.assertHttpAccepted(resp)
            self.assertTrue(resp.status_code in [202, 204], str((resp)))
            
            logger.debug(str(('check patched data for',resource_name,
                             'execute get on:',resource_uri)))
            resp = self.api_client.get(
                resource_uri, format='json', authentication=self.get_credentials(), 
                data=data_for_get)
            logger.debug(str(('--------resp to get:', resp.status_code)))
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
            #TODO: GET the apilogs expected and test them
            
            
            
            
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
        data_for_get.setdefault('HTTP_APILOG_COMMENT', 'put_test: %s' % filename )
        resource_uri = BASE_URI + '/' + resource_name

        with open(filename) as bootstrap_file:
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())
            
            logger.debug(str(('Submitting put...', bootstrap_file)))
            resp = self.testApiClient.put(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            logger.debug(str(('Response: ' , resp.status_code)))
#            self.assertHttpAccepted(resp)
            self.assertTrue(resp.status_code in [200, 202, 204], 
                            str((resp.status_code, resp.serialize() )) )
    
            logger.debug(str(('check put data for',resource_name,
                             'execute get on:',resource_uri)))
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get)
            logger.debug(str(('--------resp to get:', resp.status_code)))
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp.serialize())))
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
                
            #TODO: GET the apilogs expected and test them
            
            
            
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
        logger.debug('------------- _bootstrap_init_files -----------------')
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

  
        logger.debug('------------- done _bootstrap_init_files -----------------')

        
class TestApiInit(MetaHashResourceBootstrap,ResourceTestCase):
    
    def setUp(self):
        super(TestApiInit, self).setUp()
        super(TestApiInit, self)._setUp()

    def test0_bootstrap_metahash(self):
        
        logger.debug('================ reports test0_bootstrap_metahash =============== ')
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
            self.assertTrue(resp.status_code in [201], 
                str((resp.status_code, resp, 'resource_uri', self.resource_uri,
                    'item', item)))
            #             self.assertHttpCreated(resp)
            
        logger.debug('created items, now get them')
        resp = self.api_client.get(
            self.resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp, resp.status_code)))
        new_obj = self.deserialize(resp)
        logger.debug(str(('deserialized object:', json.dumps(new_obj))))
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        self.assertEqual(len(new_obj['objects']), 4, str((new_obj)))
        
        for inputobj in bootstrap_items:
            logger.debug(str(('testing:', inputobj)))
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, str(('not found', inputobj, outputobj )) )
        
        logger.debug('================ test0_bootstrap_metahash done ========== ')
    
    def test1_bootstrap_init(self):
        logger.debug('================ reports test1_bootstrap_init ================')
        self._bootstrap_init_files()
        logger.debug('================ test1_bootstrap_init done ================')
        
    def test2_api_init(self):
        
        logger.debug('***================ reports test2_api_init =============== ')
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
                
                logger.debug('\n++++=========== processing action', json.dumps(action))
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
                    logger.debug(str(('+++++++++++processing file', filename)))
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
        logger.debug('============== User setup ============')
        super(UserResource, self).setUp()
        super(UserResource, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
        super(UserResource, self)._bootstrap_init_files()
        logger.debug('============== User setup: begin ============')
        
        # Create the User resource field entries
        testApiClient = TestApiClient(serializer=self.csv_serializer) 
        filename = os.path.join(self.directory,'metahash_fields_user.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.user'})
        
        self.resource_uri = BASE_URI + '/user'
        logger.debug('============== User setup: done ============')
    
    def test0_create_user(self):
        logger.debug(str(('==== test_create_user =====')))
        
        # the simplest of tests, create some simple users
        self.bootstrap_items = [   
            {
                'ecommons_id': 'st1',
                'first_name': 'Sally',
                'last_name': 'Tester', 
                'email': 'sally.tester@limstest.com',    
            },
            {
                'username': 'jt1',
                'first_name': 'Joe',
                'last_name': 'Tester',    
                'email': 'joe.tester@limstest.com',    
            },
            {
                'username': 'bt1',
                'first_name': 'Bad',
                'last_name': 'TestsALot',    
                'email': 'bad.tester@slimstest.com',    
            },
        ]

        for i,item in enumerate(self.bootstrap_items):  
            try:       
                resp = self.api_client.post(self.resource_uri, 
                    format='json', data=item, authentication=self.get_credentials())
                self.assertTrue(resp.status_code in [201], str((resp.status_code, resp.serialize())))
                #                 self.assertHttpCreated(resp)
            except Exception, e:
                logger.error(str(('on creating', item, '==ex==', e)))
                raise
            
        logger.debug('created items, now get them')
        resp = self.api_client.get(self.resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        self.assertEqual(len(new_obj['objects']), 3, str((new_obj)))
        
        for i,item in enumerate(self.bootstrap_items):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(
                result, str(('bootstrap item not found', item, new_obj['objects'])))
            logger.debug(str(('item found', obj)))

        logger.debug(str(('==== test_create_user done =====')))

    def test1_create_user_with_permissions(self):
        logger.debug(str(('==== test1_user_test_data =====')))
        
        filename = os.path.join(self.directory,'test_data/users1.csv')
        # NOTE: verify the username manually, because the username will
        # be set to the ecommons id if not set
        self._put_test('user', filename, keys_not_to_check=['username'])
        
        logger.debug(str(('==== test1_user_test_data done =====')))
    
    def test2_patch_user_permissions(self):
        logger.debug(str(('==== test2_patch_user_permissions =====')))
        self.test1_create_user_with_permissions()
        
        logger.debug(str(('==== test2_patch_user_permissions start =====')))
        filename = os.path.join(self.directory,'test_data/users2_patch.csv')
        self._patch_test('user', filename)
        
        logger.debug(str(('==== test2_patch_user_permissions done =====')))
        
    def test3_user_read_permissions(self):
        '''
        Try to do something we don't have permissions for - 
        read a resource (metahash) 
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        '''
        
        logger.debug(str(('==== test3_user_read_permissions =====')))
        self.test2_patch_user_permissions()
                
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # Try to do some unauthorized actions

        resource_uri = BASE_URI + '/metahash'
        resp = self.testApiClient.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [403], str((resp.status_code, resp)))
        
        # Now add the needed permission
        
        user_patch = {
            'resource': 'user/' + username,
            'permissions': ['permission/resource/metahash/read'] };

        uri = self.resource_uri + '/' + username
        logger.debug(str(('=====now add the permission needed to this user:', user_patch, uri)))
        resp = self.api_client.patch( uri, 
                    format='json', data=user_patch, 
                    authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        str((resp.status_code, resp.serialize())))  
        
        #         # is it set?
        #         resp = self.api_client.get(
        #                 uri, format='json', authentication=self.get_credentials())
        #         logger.warn(str(('response: ' , self.deserialize(resp) )))

        # now try again as the updated user:
        
        resp = self.testApiClient.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [200], 
                        str((resp.status_code, resp.serialize())))
       
        logger.debug(str(('==== test3_user_read_permissions done =====')))

    def test4_user_write_permissions(self):
        '''
        Try to do something we don't have permissions for - 
        write a resource (user),
        then give that permission to the user, and try again 
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        '''
        logger.debug(str(('==== test4_user_write_permissions =====')))

        self.test2_patch_user_permissions()
                
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # Try to do some unauthorized actions

        resource_uri = BASE_URI + '/usergroup'
        
        # just adding a group, for instance
        
        json_data = { 'objects': [{ 'name': 'test_group_x' }] }
        
        resp = self.testApiClient.patch(
            resource_uri, format='json', data=json_data, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [403], str((resp.status_code, resp.serialize())))
        
        # Now add the needed permission
        
        user_patch = {
            'resource': 'user/' + username,
            'permissions': ['permission/resource/usergroup/write'] };

        logger.debug(str(('now add the permission needed to this user:', user_patch)))
        uri = self.resource_uri + '/' + username
        resp = self.api_client.patch( uri, 
                    format='json', data=user_patch, 
                    authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        str((resp.status_code, resp.serialize())))  
        
        # now try again as the updated user:
        
        resp = self.testApiClient.patch(
            resource_uri, format='json', data=json_data, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [202, 204], 
                        str((resp.status_code, resp.serialize())))

        #         # is it set?
        #         resp = self.api_client.get(
        #                 uri, format='json', authentication=self.get_credentials())
        #         logger.warn(str(('response: ' , self.deserialize(resp) )))

        logger.debug(str(('==== test4_user_write_permissions done =====')))
        
class UserGroupResource(UserResource):
    
    def setUp(self):
        logger.debug(str(( '============== UserGroup setup ============')))
        super(UserGroupResource, self).setUp()
        
        logger.debug(str(( '============== UserGroup setup: begin ============')))
        meta_resource_uri = BASE_URI + '/metahash'
        
        # Create the User resource field entries
        # todo: doesn't work for post, see TestApiClient.post() method, 
        # it is incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=self.csv_serializer) 
        
        filename = os.path.join(self.directory,'metahash_fields_usergroup.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.usergroup'})
        
        filename = os.path.join(self.directory,'metahash_fields_permission.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.permission'})

#         self.resource_uri = BASE_URI + '/usergroup'
        logger.debug(str(( '============== UserGroup setup done ============')))

    def test2_create_usergroup_with_permissions(self):
        logger.debug(str(('==== test2_usergroup =====')))
        #create users
        self.test1_create_user_with_permissions()
        
#         
#         data_for_get = {}
#         data_for_get.setdefault('limit', 999 )
#         resp = self.api_client.get(
#                 '/reports/api/v1/user', format='json', 
#                 authentication=self.get_credentials(), data=data_for_get)
#         logger.debug(str(('--------resp to get:', resp.status_code)))
#         self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
#         new_obj = self.deserialize(resp)

        logger.debug(str(('----- test2_usergroup =====')))

        filename = os.path.join(self.directory,'test_data/usergroups1.csv')
        # note: excluding sub_groups here because the one sub_group is set when 
        # "testGroupX" sets super_groups=['testGroup3']; and thereby testGroup3 
        # gets sub_groups=['testGroupX']; even though that's not in the input file.
        self._put_test('usergroup', filename, keys_not_to_check=['sub_groups'])

        logger.debug(str(('==== test2_usergroup done =====')))

    def test3_patch_users_groups(self):
        logger.debug(str(('==== test3_patch_users_groups =====')))
        self.test2_create_usergroup_with_permissions()
        
        logger.debug(str(('==== test3_patch_users_groups start =====')))
        filename = os.path.join(self.directory,'test_data/users3_groups_patch.csv')
        self._patch_test('user', filename) #, keys_not_to_check=['groups'])
      
        logger.debug(str(('==== test3_patch_users_groups done =====')))

#     def test3a_patch_usersgroups_groups(self):
#         logger.debug(str(('==== test3a_patch_usersgroups_groups =====')))
#         self.test2_create_usergroup_with_permissions()
#         
#         logger.debug(str(('==== test3a_patch_usersgroups_groups start =====')))
#         filename = os.path.join(self.directory,'test_data/users3_groups_patch.csv')
#         self._patch_test('user', filename, keys_not_to_check=['groups'])
#       
#         logger.debug(str(('==== test3a_patch_usersgroups_groups done =====')))


    def test4_user_group_permissions(self):
        '''
        Test of a group permission -
        first test that the user in the group has the permission,
        then remove this user from this group and try again
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.debug(str(('==== test4_user_group_permissions =====')))
        
        self.test2_create_usergroup_with_permissions()
                        
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # 1 read test - should have permission through group
        resource_uri = BASE_URI + '/metahash'
        resp = self.testApiClient.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp.serialize())))
        
        #         # is it set?
        #         uri = self.resource_uri + '/' + username
        #         resp = self.api_client.get(
        #                 uri, format='json', authentication=self.get_credentials())
        #         logger.info(str(('is it set? response: ' , self.deserialize(resp) )))

        # now patch this user's usergroups, removing the user from the group 'testgroup1'
        # which will remove the permissions as well 
        user_patch = {
            'usergroups': ['usergroup/testGroup3'] };

        logger.debug(str(('now reset this users groups and remove testGroup1:', user_patch)))
        uri = BASE_URI + '/user' + '/' + username
        resp = self.api_client.patch(uri, format='json', data=user_patch, 
                                     authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        str((resp.status_code, resp.serialize())))  
        
        #         # is it set?
        #         resp = self.api_client.get(
        #                 uri, format='json', authentication=self.get_credentials())
        #         logger.info(str(('is it set? response: ' , self.deserialize(resp) )))

        # now try again as the updated user:
        resp = self.testApiClient.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [403], 
                        str((resp.status_code, resp.serialize())))
       
        logger.debug(str(('==== test4_user_group_permissions done =====')))


    def test5_usergroup_can_contain_group_permissions(self):
        '''
        Test of an inherited group permission  -
        first test that the user in the group doesn't have the permission,
        then add the user's usergroup to a group with the permission and try again
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.debug(str(('==== test5_usergroup_can_contain_group_permissions =====')))
        
        self.test2_create_usergroup_with_permissions()
                        
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # 1 read test - user, user's group don't have the permission
        resource_uri = BASE_URI + '/vocabularies'
        resp = self.testApiClient.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(resp.status_code in [403], str((resp.status_code, resp.serialize())))
        
        # now create a new group, with previous user's group as a member,
        # then add permissions to this new group to read (vocabularies)
        # note: double nest the groups also as a test
        usergroup_patch = { 'objects': [
            {
            'name': 'testGroup5',
            'super_groups': ['usergroup/testGroup3'] },
            {
            'name': 'testGroup6',
            'users': ['user/sde4'],
            'super_groups': ['usergroup/testGroup5'] },
        ]}
        
        logger.debug(str(('now set the new group:', usergroup_patch)))
#         uri = self.resource_uri + '/'
        uri = BASE_URI + '/usergroup'
        
        resp = self.api_client.patch(uri, format='json', 
            data=usergroup_patch, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        str((resp.status_code, resp.serialize())))  
        
        # is it set?
        resp = self.api_client.get(
                uri, format='json', authentication=self.get_credentials())
        new_obj = self.deserialize(resp)
        result, outputobj = find_all_obj_in_list(
            usergroup_patch['objects'],new_obj['objects']) #, excludes=keys_not_to_check )
        self.assertTrue(
            result, 
            str(('not found', outputobj,'=== objects returned ===', 
                 new_obj['objects'] )) ) 
        
        # is it set-2, does group inherit the permissions?
        for obj in new_obj['objects']:
            if obj['name'] == 'testGroup6':
                testGroup6 = obj
            if obj['name'] == 'testGroup3':
                testGroup3 = obj
                
        self.assertTrue(testGroup6 and testGroup3)
        for permission in testGroup3['all_permissions']:
            logger.info(str(('find permission', permission)))
            self.assertTrue(permission in testGroup6['all_permissions'], 
                str(('could not find permission', permission, 
                     'in testGroup6 permissions', testGroup6['all_permissions'])))
        
        # 2 read test - user has permissions through inherited permissions,
        resource_uri = BASE_URI + '/vocabularies'
        resp = self.testApiClient.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp.serialize())))

        logger.debug(str(('==== Done: test5_usergroup_can_contain_group_permissions =====')))

    
    def test6_usergroup_can_contain_group_users(self):
        '''
        Test that group "contains" the subgroup users (through "all_users")
        first test that the group doesn't have the user, then add the subgroup
        with the user, test that the group has the user in "all_users
        - done in prev test: create a new user (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.debug(str(('==== test6_usergroup_can_contain_group_users =====')))
        
        self.test2_create_usergroup_with_permissions()
                        
        # now create a new group, with a previous group as a sub_group
        usergroup_patch = { 'objects': [
            {
            'name': 'testGroup5',
            'sub_groups': ['usergroup/testGroup2'] },
            {
            'name': 'testGroup6',
            'sub_groups': ['usergroup/testGroup5'] },
        ]}
        
        logger.debug(str(('now set the new groups:', usergroup_patch)))
#         uri = self.resource_uri + '/'
        uri = BASE_URI + '/usergroup'
        resp = self.api_client.patch(uri, format='json', 
            data=usergroup_patch, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        str((resp.status_code, resp.serialize())))  
        
        # is it set?
        resp = self.api_client.get(
                uri + '/testGroup6', format='json', authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [200], 
                        str((resp.status_code, resp.serialize())))
        new_obj = self.deserialize(resp)
        # TODO: review: TP returns all dict/keys as unicode
#         new_obj = dict_strip_unicode_keys(new_obj)
#         logger.info(str(('results from get', new_obj)))
        
        self.assertTrue(new_obj['all_users'])
        self.assertTrue('user/sde4' in new_obj['all_users'])
        
        # TODO: could also test that testGroup2 now has super_group=testGroup5
        
        logger.debug(str(('==== Done: test6_usergroup_can_contain_group_users =====')))
