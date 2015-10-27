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

from __future__ import unicode_literals

import cStringIO
import json
import logging
import os
import unittest
from unittest.util import safe_repr

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.db.utils import ProgrammingError
from django.http.response import StreamingHttpResponse
from django.test import TestCase
from django.test.simple import DjangoTestSuiteRunner
from django.test.testcases import SimpleTestCase
from django.utils.encoding import force_text
from tastypie import fields
from tastypie.fields import BooleanField
from tastypie.test import ResourceTestCase, TestApiClient
from tastypie.utils.dict import dict_strip_unicode_keys

from reports import dump_obj, HEADER_APILOG_COMMENT
from reports.api import compare_dicts
from reports.dump_obj import dumpObj
from reports.models import API_ACTION_CREATE
from reports.serializers import CsvBooleanField, CSVSerializer, SDFSerializer, \
    LimsSerializer, XLSSerializer
from reports.sqlalchemy_resource import SqlAlchemyResource
import reports.utils.log_utils
from reports.utils.sdf2py import MOLDATAKEY
import reports.utils.serialize


_ = reports.utils.log_utils.LogMessageFormatter   # optional, to improve readability
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

def numerical_equivalency(val, val2):
    try:
        if float(val) == float(val2):
            return True
    except:
        pass
    return False

####
# NOTE: equivocal, and other equivalency methods are for end-to-end testing, 
# using string equals, of the API with values read from csv or json.
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
        if numerical_equivalency(val1, val2):
            return True, ('val1', val1, 'val2', val2 )

        val1 = str(val1)
        val2 = str(val2)
        if val1 != val2:
            return False, ('val1', val1, 'val2', val2 )
        
    
    if isinstance(val1, basestring):
        val1 = str(val1)
        if hasattr(val2, "__getitem__") or hasattr(val2, "__iter__"): 
            # allow single item lists to be equal to string
            if len(val2) == 1:
                val2 = val2[0]
        val2 = str(val2)
        if val1 != val2:
            if (is_boolean(val1) and 
                    csvBooleanField.convert(val1)==csvBooleanField.convert(val2)):
                return True, ('val1', val1, 'val2', val2 )
            elif numerical_equivalency(val1,val2):
                return True, ('val1', val1, 'val2', val2 )
                
            return False, ('val1', val1, 'val2', val2 )
    else: # better be a list
        if not isinstance(val1, list) and isinstance(val2, list):
            return False, ('Must be a list if not a string', 'val1', val1, 'val2', val2)
        if not val2:
            return False, ('val2 is None', 'val1', val1 )
        for v in val1:
            if not v: 
                if val2 and len(val2) > 0:
                    return False, ('val1', val1, 'list val2 not empty', val2 )
            if v not in val2:
                if v not in [str(v2) for v2 in val2]:
                    return False, ('val1', val1, 'val2', val2 )
    return True, ('val1', val1, 'val2', val2 )
    
def assert_obj1_to_obj2( obj1, obj2, keys=[], excludes=['resource_uri']):
    '''
    For testing equality of the (CSV) input to the response of what is 
    created and deserialized from the DB
    - NOTE only works for String comparisons
    @param obj1 input
    @param obj2 output 
    '''
    original_keys = set(obj1.keys())
    original_keys = original_keys.difference(excludes)
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
        #    logger.debug(_('equal so far, keys to search', keys_to_search))

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
                        return False, ('resource uri not equal', val1, val2, obj2)
                    logger.warn(_(
                        'imprecise uri matching, equivocating:', val1, val2))
    
    keys_to_search = set(obj1.keys()) - set(keys) - set(excludes)
    # logger.debug(_('keys to search', keys_to_search,excludes))
    for key in keys_to_search:
        result, msgs =  equivocal(obj1[key], obj2[key])
        if not result:
            return False, ('key not equal', key, obj1[key], obj2[key], 
                'obj1', obj1, 'obj2', obj2, 'msgs', msgs)
            
    return True, ('obj1:', obj1, 'obj2:', obj2)


def find_obj_in_list(obj, item_list, id_keys_to_check=[], **kwargs):
    '''
    @param id_keys_to_check is specified, then "find" the first obj in the list
        based only on these id keys (a secondary check can be done for complete 
        equality
    '''
    list_msgs = []
    for item in item_list:
        if id_keys_to_check:
            found = item
            for key in id_keys_to_check:
                if key not in obj:
                    raise AssertionError('id key is not in obj', (key, obj))
                if key in item and obj[key] != item[key]:
                    found = None
                    break
            if found:
                return True, (item)    
        result, msgs = assert_obj1_to_obj2(obj, item, **kwargs)
        if result:
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
            return False, msgs
    return True, msgs

# To Run tests without a database:
# Example:
# ./manage.py test reports.SDFSerializerTest.test2_clean_data_sdf \
#      --settings=lims.settings_testing_debug --verbosity=2 \
#      --testrunner=reports.tests.NoDbTestRunner
# NOTE: when using this testrunner, the class finder is different; note that the 
# path excludes the module, so 
# "reports.SDFSerializerTest.test2_clean_data_sdf"
# not 
# "reports.test.SDFSerializerTest.test2_clean_data_sdf"
class NoDbTestRunner(DjangoTestSuiteRunner):
  """ A test runner to test without database creation """

  def setup_databases(self, **kwargs):
    """ Override the database creation defined in parent class """
    pass

  def teardown_databases(self, old_config, **kwargs):
    """ Override the database teardown defined in parent class """
    pass


# Create a TestRunner that will swallow exceptions on teardown of the test
# database.
# - motivation: use this runner for Travis tests so that the teardown error is
# not reported as a testing failure.
# 
# FIXME for django 1.7 - use django.test.runner.DiscoverRunner
# FIXME: override of testrunner class to catch teardown error:
# - caused by holding another db connection for the sqlalchemy bridge
# NOTE: when using this testrunner, the class finder is different; note that the 
# path excludes the module, so 
# "reports.SDFSerializerTest.test2_clean_data_sdf"
# not 
# "reports.test.SDFSerializerTest.test2_clean_data_sdf"
class IccblTestRunner(DjangoTestSuiteRunner):
    
    def teardown_databases(self, old_config, **kwargs):
        try:
            DjangoTestSuiteRunner.teardown_databases(self, old_config, **kwargs)
        except Exception, e:
            logger.exception('on teardown')
            
def print_find_errors(outputobj):
    '''
    print nested error found in an end-to-end equivalency test.
    '''
    for x in outputobj:
        if not isinstance(x, (basestring, dict)):
            for y in x:
                if not isinstance(x, (basestring, dict)):
                    for z in x:
                        if isinstance(z, (list,tuple)):
                            for zz in z:
                                print 'zz', zz, '\n'
                        print 'a', z, '\n'
                else:
                    if isinstance(x, basestring):
                        for z in x.split('msgs'), '\n':
                            print 'b',z, '\n'
                    else:
                        print 'c',z, '\n'
        else:
            if isinstance(x, basestring):
                for z in x.split('msgs'), '\n':
                    print 'd',z, '\n'
            else:
                print 'e',z, '\n'
        
class BaselineTest(TestCase):
    '''
    Test the equivalency methods.
    '''
    
    def test0_equivocal(self):
        logger.debug(_('>> test0_equivocal <<'))
        
        test_true = [
            ['',''],
            ['', None],
            ['1','1'],
            # Test that integer strings can come back as integers 
            # (if integerfield is defined)
            ['1',1], 
            ['sal','sal'],
            ['True','TRUE'],
            ['true','True'],
            # Test that boolean strings can come back as booleans 
            # (if booleanfield is defined)
            ['TRUE',True], 
            ['false','False'],
            ['FALSE','false'],
            ['False',False],
            ]
            ## TODO: test date strings 
        
        for [a,b] in test_true:
            logger.warn(_('csv serialization equivocal truthy test', a, b))
            result, msgs = equivocal(a,b)
            self.assertTrue(result, msgs)

        test_false = [
            ['1',''],
            # Test that string is not equal to the wrong integer after conversion
            ['1',2],    
            ['','1'],
            ['sal','val'],
            ['True','1'],
            ['true',''],
            ['TRUE',False],
            ['false','true'],
            ['False',True],
            ]
        for [a,b] in test_false:
            logger.debug(_('csv serialization equivocal falsy test', a, b))
            result, msgs = equivocal(a,b)
            self.assertFalse(result, msgs)
        logger.debug(_('>> test0_equivocal done <<'))
    
    def test1_assert_obj1_to_obj2(self):
        logger.debug(_('====== test1_assert_obj1_to_obj2 ===='))

        # test that all of obj1 in obj2 (may be more in obj2)
        obj1 = { 'one': '1', 'two': 'two', 'three':''}
        obj2 = { 'one': '1', 'two': 'two', 'three':'', 'four': 'blah'}
        result, msgs = assert_obj1_to_obj2(obj1, obj2)
        self.assertTrue(result, msgs)
        result, msgs = assert_obj1_to_obj2(obj2, obj1)
        self.assertFalse(result, msgs)
        obj3 = { 'one':'1', 'two':'2', 'three':'' }
        result, msgs = assert_obj1_to_obj2(obj3, obj1)
        self.assertFalse(result, _(result,msgs))

        obj4 = { 'one':'1', 'two':'2', 'three':'' }
        obj5 = { 'one':'1', 'two':'2', 'three':None }
        result, msgs = assert_obj1_to_obj2(obj5, obj4)
        self.assertTrue(result, _(result,msgs))
        logger.debug(_('====== test1_assert_obj1_to_obj2 done ===='))

    def test2_find_obj_in_list(self):
        logger.debug(_('====== test2_find_obj_in_list ===='))
        
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
        logger.debug(_(result,msgs))
        self.assertFalse(result, msgs)

        obj_list = [ obj1a, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertFalse(result, msgs)
        logger.debug(_('====== test2_find_obj_in_list done ===='))

    def test3_find_all_obj_in_list(self):
        logger.debug(_('====== test3_find_all_obj_in_list ===='))
        
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
        self.assertFalse(result, (result, msgs  ))
        logger.debug(_('====== test3_find_all_obj_in_list done ===='))

    def test4_user_example(self):
        logger.debug(_('====== test4_user_example ===='))
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
                'login_id': 'ut1',
                'first_name': 'userx',
                'last_name': 'Tester1',    
                'email': 'user.tester1@slimstest.com',    
            },
        ]        
        item = {'login_id': 'jt1', 'first_name': 'Joe', 'last_name': 'Tester', 
                'email': 'joe.tester@limstest.com'}
        result, obj = find_obj_in_list(item, bootstrap_items)
        logger.debug(_(result, obj))
        
        self.assertTrue(obj['email'] == 'joe.tester@limstest.com')
        
        logger.debug(_('====== test4_user_example done ====')) 
    
    
class SerializerTest(TestCase):

    def test_csv(self):
        logger.debug(_('======== test_csv ========='))
        directory = APP_ROOT_DIR
        serializer = CSVSerializer() 
        
        input = [['one','two', 'three', 'four', 'five','six'],
                ['uno', '2', 'true', 'false', '',['a','b','c']]]
        
        input_data = reports.utils.serialize.from_csv_iterate(input)
        for obj in input_data:
            self.assertTrue(obj['one']=='uno')
            self.assertTrue(obj['two']=='2')
            self.assertTrue(obj['three']=='true')
            self.assertTrue(obj['four']=='false')
            self.assertTrue(obj['five']=='')
            self.assertTrue(obj['six']==['a','b','c'])
        logger.debug(_('input to object', input_data))
        csv_data = serializer.to_csv(input_data, root=None)
        logger.debug(_('back to csv', csv_data))

        with open(directory + '/reports/test/test_csv_.csv', 'w') as f:
            f.write(csv_data)
            
        with open(directory + '/reports/test/test_csv_.csv') as fin:    
            final_data = serializer.from_csv(fin.read(), root=None)
            logger.debug(_('final_data', final_data))
            for obj in final_data:
                self.assertTrue(obj['one']=='uno')
                self.assertTrue(obj['two']=='2')
                self.assertTrue(obj['three']=='true')
                self.assertTrue(obj['four']=='false')
                self.assertTrue(obj['five']=='')
                self.assertTrue(obj['six']==['a','b','c'])
        
        # TODO: delete the file
        logger.debug(_('======== test_csv done ========='))

class SDFSerializerTest(SimpleTestCase):
    
    def test1_read(self):
        logger.debug('=== test1 SDF read')
        
        records = [{
            'smiles': 'Cl',
            'compound_name': 'HCl (hydrochloric acid)',
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
            'compound_name': 'malate',
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
            'compound_name': 'hydrochloride hydrate',
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
        filename = "lims/static/test_data/test1_output.sdf"
        sdf_data = serializer.to_sdf(records)
        with open(os.path.join(APP_ROOT_DIR, filename), 'w') as fout:    
            fout.write(sdf_data)
        fout.close()
    
        with open(os.path.join(APP_ROOT_DIR, filename)) as fin:    
            _data = serializer.from_sdf(fin.read(), root=None)
            final_data = _data
            logger.debug(_('final_data', final_data))
            
            self.assertTrue(
                len(final_data)==len(records), 
                str(('len is', len(final_data),len(records))))
            for obj in final_data:
                logger.debug(_('object: ', obj))
                
                for record in records:
                    if record['smiles'] == obj['smiles']:
                        for k,v in record.items():
                            self.assertTrue(
                                obj[k] == v.strip(), 
                                ('values not equal', k, record[k], 'read:', v))

    def test2_clean_data_sdf(self):
        ''' more extensive, read/write test '''
        logger.debug(_('==== test2_clean_data_sdf ====='))

        record_one = {
            r'vendor': r'Biomol-TimTec',
            r'vendor_reagent_id': r'SPL000058',
            r'vendor_batch_id': r'HM-001_TM-20090805',
            r'plate_number': r'1536',
            r'well_name': r'A01',
            r'library_well_type': r'experimental',
            r'facility_reagent_id': r'ICCB-00589081',
            r'facility_batch_id':r'008',
            r'compound_name': [r'fake compound name 1',r'fake compound name 2'],
            r'smiles':r'O=C1CC(C)(C)CC(=O)C1C(c1ccccc1)C1=C(O)CC(C)(C)CC1=O',
            r'inchi':r'InChI=1/C23H28O4/c1-22(2)10-15(24)20(16(25)11-22)19(14-8-6-5-7-9-14)21-17(26)12-23(3,4)13-18(21)27/h5-9,19-20,26H,10-13H2,1-4H3',
            r'pubchem_cid': [r'558309',r'7335957'],
            r'chembank_id':[r'1665724',r'6066882'],
            r'mg_ml_concentration':r'.111',
            r'chembl_id':[r'100001',r'100002',r'111102'],
            r'pubmed_id':[r'20653109',r'20653081'] }

        record_two = {
            r'Library': r'Biomol-TimTec1',
            r'Source': r'Biomol-TimTec',
            r'vendor': r'Biomol-TimTec',
            r'vendor_reagent_id': r'ST001215',
            r'plate_number': r'1536',
            r'well_name': r'A02',
            r'library_well_type': r'experimental',
            r'Row': r'1',
            r'Col':r'2',
            r'facility_reagent_id': r'ICCB-00589082',
            r'molar_concentration': r'.005',
            r'compound_name':r'fake compound name 1',
            r'CAS_Number': r'fake cas number 1',
            r'smiles':r'Clc1ccc(\\C=C/c2c(C)n(C)n(c3ccccc3)c2=O)c(Cl)c1',
            r'inchi':r'InChI=1/C23H28O4/c1-22(2)10-15(24)20(16(25)11-22)19(14-8-6-5-7-9-14)21-17(26)12-23(3,4)13-18(21)27/h5-9,19-20,26H,10-13H2,1-4H3',
            MOLDATAKEY:
r'''Structure89
csChFnd70/04290511482D

 16 17  0  0  0  0  0  0  0  0999 V2000
    4.8373    2.0813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8373    0.7093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6256    2.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0491    2.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4179    2.0813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4179    0.7093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.0813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2308    6.2888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2308    4.9141    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2076    2.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2076    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0491    4.2034    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2308    2.1223    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0491    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.6256    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  1  2  0  0  0  0
 16  2  1  0  0  0  0
  5  1  1  0  0  0  0
  6  3  1  0  0  0  0
  7  6  2  0  0  0  0
 15  2  2  0  0  0  0
 14  5  2  0  0  0  0
 13  5  1  0  0  0  0
 11  6  1  0  0  0  0
 12  7  1  0  0  0  0
 10 13  1  0  0  0  0
  9 10  1  0  0  0  0
  8 11  2  0  0  0  0
  4  8  1  0  0  0  0
 16  7  1  0  0  0  0
  4 12  2  0  0  0  0
M  END'''            }
        last_record = {
            'Library': 'Biomol-TimTec1',
            'Source': 'Biomol-TimTec',
            'vendor': 'Biomol-TimTec',
            'vendor_reagent_id': '',
            'plate_number': '1536',
            'well_name': 'A09',
            'library_well_type': 'empty',
            'facility_reagent_id': '',
            'compound_name': '',
            'smiles': '',
            'inchi': '',
            }
        
        
        serializer = SDFSerializer()
        filename = ( APP_ROOT_DIR 
            + '/db/static/test_data/libraries/clean_data_small_molecule.sdf' )

        with open(os.path.join(APP_ROOT_DIR, filename)) as fin:    
            _data = serializer.from_sdf(fin.read(), root=None)
            input_data = _data
            
            if logger.isEnabledFor(logging.DEBUG):
                print 'data read in'
                for x in input_data:
                    print x, '\n'
            
            expected_count = 8
            self.assertEqual(len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            result, msgs = assert_obj1_to_obj2(record_one, input_data[0])
            self.assertTrue(result, msgs)
 
            if logger.isEnabledFor(logging.DEBUG):
                print 'record 1 expected:'
                print input_data[1][MOLDATAKEY]
                print '====='
                print 'record 2 read:'
                print record_two[MOLDATAKEY]
            
                for i,c in enumerate(input_data[1][MOLDATAKEY]):
                    if record_two[MOLDATAKEY][i] != c:
                        print 'i', i, c,record_two[MOLDATAKEY][i]
                        break 

            result, msgs = assert_obj1_to_obj2(record_two, input_data[1])
            self.assertTrue(result, msgs)
 
            result, msgs = assert_obj1_to_obj2(last_record, input_data[-1])
            self.assertTrue(result, msgs)
 
            # Now test the whole system by writing back out and reading back in
            
            out_filename = os.path.join(APP_ROOT_DIR, filename + 'out')
            sdf_data = serializer.to_sdf(input_data)
            with open(out_filename, 'w') as fout:    
                fout.write(sdf_data)
            fout.close()

            with open(out_filename) as fin:    
                _fdata = serializer.from_sdf(fin.read(), root=None)
                final_data = _fdata
                
                self.assertEqual(len(input_data), len(final_data), 
                    str(('initial serialization of ',out_filename,'found',
                        len(final_data), 'expected',len(input_data),
                        'final_data',final_data)))
                
                keys_not_to_check=[]
                for i,inputobj in enumerate(input_data):
                    result, outputobj = find_obj_in_list(
                        inputobj,final_data, excludes=keys_not_to_check )
                    if not result:
                        print 'input obj not found'
                        print inputobj, '\n'
                        print 'final data read in:'
                        for x in final_data:
                            print x , '\n'
                        print 'messages'
                        print_find_errors(outputobj)
                        
                        self.fail('input object not found')


class XLSSerializerTest(SimpleTestCase):
    
    def test1_read(self):
        logger.debug('==== test1 XLS read')
        
        records = [{
            'key1': 'aval1',
            'key2': 'aval2',
            'key3': 'aval3',
            'key4': ['a','b','c']           
            },{
            'key1': 'bval1',
            'key2': 'bval2',
            'key3': 'bval3'            
            },
            ]

        serializer = XLSSerializer()
        filename = "lims/static/test_data/test1_output.xls"
        _data = serializer.to_xls(records)
        with open(os.path.join(APP_ROOT_DIR, filename), 'w') as fout:    
            fout.write(_data)
            
        fout.close()
    
        with open(os.path.join(APP_ROOT_DIR, filename)) as fin:    
            _data = serializer.from_xls(fin.read(), root=None)
            final_data = _data
            logger.debug(_('final_data', final_data))
            
            self.assertTrue(
                len(final_data)==len(records), 
                str(('len is', len(final_data),len(records))))
            for obj in final_data:
                logger.debug(_('object: ', obj))
                
                for record in records:
                    if record['key1'] == obj['key1']:
                        for k,v in record.items():
                            self.assertTrue(
                                obj[k] == v, 
                                ('values not equal', k, record[k], 'read:', v))

    def test2_clean_data(self):

        serializer = XLSSerializer()
        
        test_input_data = [
            {
                'plate_number': '50001', 
                'well_name': 'A05', 
                'library_well_type': 'experimental', 
                'molar_concentration': '0.0001', 
                'mg_ml_concentration': '.115', 
                'vendor': 'vendorX', 
                'vendor_reagent_id': 'M-005300-00', 
                'facility_reagent_id': 'F-005300-00', 
                'silencing_reagent_type': 'sirna', 
                'sequence': 'GACAUGCACUGCCUAAUUA;GUACAGAACUCUCCCAUUC;GAUGAAAUGUGCCUUGAAA;GAAGGUGGAUUUGCUAUUG', 
                'anti_sense_sequence': 'GACAUGCACUGCCUAAUUA;GUACAGAACUCUCCCAUUC;GAUGAAAUGUGCCUUGAAA;GAAGGUGGAUUUGCUAUUA', 
                'vendor_entrezgene_id': '22848', 
                'vendor_entrezgene_symbols': ['AAK1','AAK2'], 
                'vendor_gene_name': 'VendorGeneNameX', 
                'vendor_genbank_accession_numbers': ['NM_014911','NM_014912'], 
                'vendor_species': 'VendorSpeciesX', 
                'facility_entrezgene_id': '1111', 
                'facility_entrezgene_symbols': ['AAK3','AAK4'], 
                'facility_gene_name': 'FacilityGeneNameX',
                'facility_genbank_accession_numbers': ['F_014911','F_014914'], 
                'facility_species': 'FacilitySpeciesX', 
                },{
                'plate_number': '50001', 
                'well_name': 'A07', 
                'library_well_type': 'library_control', 
                'vendor': 'vendorX', 
                'vendor_reagent_id': 'M-000000-00', 
                'silencing_reagent_type': 'sirna', 
                'sequence': 'GUACAGAGAGGACUACUUC;GGUACGAGGUGAUGCAGUU;UCAGUGGCCUCAACGAGAA;GCAAGUACAGAGAGGACUA', 
                'anti_sense_sequence': 'GUACAGAGAGGACUACUUC;GGUACGAGGUGAUGCAGUU;UCAGUGGCCUCAACGAGAA;GCAAGUACAGAGAGGACUG', 
                'vendor_entrezgene_id': '9625', 
                'vendor_entrezgene_symbols': ['AATK'], 
                'vendor_genbank_accession_numbers': ['XM_375495'], 
                },            
            ]
        # TODO: referencing a "db" subproject from the "reports" subproject
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_rnai.xlsx'
        with open(os.path.join(APP_ROOT_DIR, filename)) as fin:    
            _data = serializer.from_xls(fin.read(), root=None)
            logger.debug(_('final_data', _data))
            
            expected_count = 5
            self.assertTrue(
                len(_data)==expected_count, 
                str(('len is', len(_data),expected_count)))

            keys_not_to_check=[]
            for i,inputobj in enumerate(test_input_data):
                result, outputobj = find_obj_in_list(
                    inputobj,_data, excludes=keys_not_to_check )
                if not result:
                    print '========input obj not found==========='
                    print inputobj, '\n'
                    print 'messages'
                    print outputobj
                    print 'final data read in:'
                    for x in _data:
                        print x , '\n'
                    print_find_errors(outputobj)
                    
                    self.fail('input object not found')

        # Test XLSX
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_rnai.xlsx'
        with open(os.path.join(APP_ROOT_DIR, filename)) as fin:    
            _data = serializer.from_xlsx(fin.read(), root=None)
            logger.debug(_('xlsx final_data', _data))
            
            expected_count = 5
            self.assertTrue(
                len(_data)==expected_count, 
                str(('len is', len(_data),expected_count)))

            keys_not_to_check=[]
            for i,inputobj in enumerate(test_input_data):
                result, outputobj = find_obj_in_list(
                    inputobj,_data, excludes=keys_not_to_check )
                if not result:
                    print '========xlsx input obj not found==========='
                    print inputobj, '\n'
                    print 'messages'
                    print outputobj
                    print 'final data read in:'
                    for x in _data:
                        print x , '\n'
                    print_find_errors(outputobj)
                    
                    self.fail('xlsx input object not found')


class HydrationTest(TestCase):
    
    def test_hydrate_boolean_tastypie(self):
        logger.debug(_('======== test_hydrate_boolean ========='))
        
        field = BooleanField()
        
        test_data = [ True, False, 'true', 'TRUE', 'FALSE', 'false', '' ]
        expected  = [ True, False, True, True, True, True, False ]
        
        for i, item in enumerate(test_data):
            result = field.convert(item)
            self.assertEqual(result, expected[i], 
                (i,' is not equal', item, result, expected[i]))
        logger.debug(_('======== test_hydrate_boolean done ========='))
            
    def test_hydrate_boolean_lims(self):
        field = CsvBooleanField()
        
        test_data = [ True, False, 'true', 'TRUE', 'FALSE', 'false', '' ]
        expected  = [ True, False, True, True, False, False, False ]
        
        for i, item in enumerate(test_data):
            result = field.convert(item)
            self.assertEqual(result, expected[i],
                (i,' is not equal', item, result, expected[i]))
    def test_hydrate_list(self):
        test_data = ['one','two','three']
        field = fields.ListField()
        
        result = field.convert(test_data)
        
        self.assertEqual(test_data, result)
  
class LogCompareTest(TestCase):
    
    def test_compare_dicts(self):
        
        dict1 = {
            'one': None,
            'two': 'value2a',
            'four': 'value4',
            'five': 'unchanged',
            'six': 'removed',
        }
        
        dict2 = {
            'one': 'value1',
            'two': 'value2b',
            'three': 'added',
            'four': None,
            'five': 'unchanged',
             }
        
        diff_dict = compare_dicts(dict1, dict2)
        
        self.assertTrue(diff_dict['added_keys'] == ['three'], diff_dict['added_keys'])
        self.assertTrue(diff_dict['removed_keys'] == ['six'], diff_dict['removed_keys'])
        self.assertTrue(set(diff_dict['diff_keys']) 
            == set(['one','two','four']), diff_dict['diff_keys'])

        # because it's not "full", diff_dict doesn't show diffs added
        self.assertTrue('three' not in diff_dict['diffs'], diff_dict['diffs'])
        self.assertTrue(diff_dict['diffs']['two']==['value2a', 'value2b'])
        

class IResourceTestCase(SimpleTestCase):
    """
    Override the Django SimpleTestCase, not using the TransactionTestCase
    necessary so that the SqlAlchemy connection can see the same database as the 
    django test code (Django ORM).

    Serves as a base class for tests, and, through _bootstrap_init_files, 
    as the setup for this module.

    TODO/FIXME: initialize the SqlAlchemyResource inside the transactions so that 
    the normal TransactionTestCase can be used.
    """
    def __init__(self,*args,**kwargs):
    
        super(IResourceTestCase, self).__init__(*args,**kwargs)
        # Create a user.
        self.username = 'testsuper'
        self.password = 'pass'
        try:
            self.user = User.objects.get(username=self.username)
        except ObjectDoesNotExist:
            self.user = User.objects.create_superuser(
                self.username, 'testsuperuser@example.com', self.password)
        
        self.resource_uri = BASE_URI + '/metahash'
        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        self.csv_serializer=CSVSerializer() 
        
        self.serializer = LimsSerializer()
        self.api_client = TestApiClient(serializer=self.serializer)       
    
    def setUp(self):
        super(IResourceTestCase, self).setUp()
 
    def create_basic(self, username, password):
        """
        Creates & returns the HTTP ``Authorization`` header for use with BASIC
        Auth.
        """
        import base64
        return 'Basic %s' % base64.b64encode(
            ':'.join([username, password]).encode('utf-8')).decode('utf-8')

    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)
        
    def get_resource_from_server(self, resource_name):
        '''
        Utility to get a resource description from the server
        '''
        resource_uri = BASE_URI + '/resource/' + resource_name
        logger.debug(_('Get the schema', resource_uri ))
        return self.get_from_server(resource_uri)
        
    def get_from_server(self, resource_uri):
        logger.debug(_('get_from_server', resource_uri))
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(_('--------resp to get:', resp.status_code))
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, self.serialize(resp)))
        return self.deserialize(resp)
    
    def _patch_test(self,resource_name, filename, keys_not_to_check=[], 
                    id_keys_to_check=[], data_for_get={}):
        '''
        @param data_for_get - dict of extra header information to send with the 
        GET request
        '''
        data_for_get.setdefault('limit', 999 )
        data_for_get.setdefault('includes', '*' )
        data_for_get.setdefault( HEADER_APILOG_COMMENT, 'patch_test: %s' % filename )
        resource_uri = BASE_URI + '/' + resource_name
        
        logger.debug(_('===resource_uri', resource_uri))
        with open(filename) as bootstrap_file:
            # NOTE / TODO: we have to deserialize the input, because the TP test method 
            # will expect a python data object, which it will serialize!
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())
            
            logger.debug(_('Submitting patch...', filename))
            resp = self.api_client.patch(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [202, 204], 
                (resp.status_code, self.deserialize(resp)))
            
            logger.debug(_('check patched data for',resource_name,
                             'execute get on:',resource_uri))
            resp1 = self.api_client.get(
                resource_uri, format='json', authentication=self.get_credentials(), 
                data=data_for_get)
            logger.debug(_('--------resp to get:', resp.status_code))
            new_obj = self.deserialize(resp1)
            self.assertTrue(resp1.status_code in [200], 
                (resp1.status_code, new_obj))
            
            for inputobj in input_data['objects']:
                # use id_keys_to_check to perform a search only on those keys
                result, outputobj = find_obj_in_list(
                    inputobj,new_obj['objects'], id_keys_to_check=id_keys_to_check, 
                    excludes=keys_not_to_check )
                self.assertTrue(
                    result, (
                    'not found', outputobj,'=== objects returned ===', 
                         new_obj['objects'])) 
                # once found, perform equality based on all keys (in obj1)
                result, msg = assert_obj1_to_obj2(inputobj, outputobj,
                    excludes=keys_not_to_check)
                self.assertTrue(result,
                    ('not equal', msg, inputobj, outputobj))
                self.assertTrue(
                    resource_name in outputobj['resource_uri'], 
                    ('wrong resource_uri returned:', filename, outputobj['resource_uri'],
                         'should contain', resource_name))
                for id_key in id_keys_to_check:
                    self.assertTrue(
                        inputobj[id_key] in outputobj['resource_uri'], 
                        ('wrong resource_uri returned:', filename, outputobj['resource_uri'],
                             'should contain id key', id_key, 'val', inputobj[id_key]))
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
        data_for_get.setdefault('includes', '*' )
        data_for_get.setdefault( HEADER_APILOG_COMMENT, 'put_test: %s' % filename )
        resource_uri = BASE_URI + '/' + resource_name

        with open(filename) as bootstrap_file:
            # NOTE / TODO: we have to deserialize the input, because the TP test method 
            # will expect a python data object, which it will serialize!
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())
            
            logger.debug(_('Submitting put...', filename))
            resp = self.api_client.put(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            logger.debug(_('Response: ' , resp.status_code))
            self.assertTrue(resp.status_code in [200, 202, 204], 
                            str((resp.status_code, self.serialize(resp) )) )
    
            logger.debug(_('check put data for',resource_name,
                             'execute get on:',resource_uri))
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get)
            logger.debug(_('--------resp to get:', resp.status_code))
            self.assertTrue(resp.status_code in [200], 
                (resp.status_code, self.serialize(resp)))
            new_obj = self.deserialize(resp)
            # do a length check, since put will delete existing resources
            self.assertEqual(len(new_obj['objects']), len(input_data['objects']), 
                str(('input length != output length: ',
                    len(new_obj['objects']), len(input_data['objects']),
                    input_data,'\n\n', new_obj)))
            
            for inputobj in input_data['objects']:
                result, outputobj = find_obj_in_list(
                    inputobj,new_obj['objects'],
                    id_keys_to_check=id_keys_to_check,
                    excludes=keys_not_to_check)
                self.assertTrue(result, ('not found', outputobj, 
                                             new_obj['objects'] ) )
                
                # once found, perform equality based on all keys (in obj1)
                result, msg = assert_obj1_to_obj2(inputobj, outputobj,
                    excludes=keys_not_to_check)
                self.assertTrue(result,
                    ('not equal', msg, inputobj, outputobj))

                logger.debug(_('outputobj[resource_uri]', outputobj['resource_uri']))
                self.assertTrue(resource_name in outputobj['resource_uri'], 
                    ('wrong resource_uri returned:', filename, outputobj['resource_uri'],
                         'should contain', resource_name))
                for id_key in id_keys_to_check:
                    self.assertTrue(inputobj[id_key] in outputobj['resource_uri'], 
                        ('wrong resource_uri returned:', outputobj,
                             'should contain id key', id_key, 'val', inputobj[id_key]))
            #TODO: GET the apilogs expected and test them

            # return both collections for further testing
            return (input_data['objects'], new_obj['objects']) 

    def _bootstrap_init_files(self):

        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        serializer=CSVSerializer() 
        testApiClient = TestApiClient(serializer=serializer) 
        input_actions_file = os.path.join(self.directory, 'api_init_actions.csv')
        logger.info(_('open input_actions file', input_actions_file))
        
        self._run_api_init_actions(input_actions_file)

    def _run_api_init_actions(self, input_actions_file):
        
        if not os.path.exists(input_actions_file):
            raise AssertionError('no such file %s' % input_actions_file)
        with open(input_actions_file) as input_file:
            api_init_actions = reports.utils.serialize.from_csv(input_file)
            if not api_init_actions:
                raise AssertionError('no actions read from file: %s' % input_actions_file)
            for action in api_init_actions:
                logger.info(_('process action', action))
                command = action['command'].lower() 
                resource = action['resource'].lower()

                if command == 'delete':
                    # no need to delete with the empty test database
                    pass
                else:
                    data_file = os.path.join(self.directory,action['file'])
                    logger.info(_('Data File: ', data_file))
                    if command == 'put':
                        self._put_test(resource, data_file, keys_not_to_check=['resource_uri'])
                    elif command == 'patch':
                        self._patch_test(resource, data_file, keys_not_to_check=['resource_uri'])
                    else:
                        raise AssertionError('Unknown API command %s' % command)
                    
    def get_content(self, resp):
        return self.serializer.get_content(resp);
#         if isinstance(resp, StreamingHttpResponse):
#             buffer = cStringIO.StringIO()
#             for line in resp.streaming_content:
#                 buffer.write(line)
#             logger.info(str(('get_content returns: ',buffer.getvalue())))
#             return buffer.getvalue()
#         else:
#             return resp.content
    
    def deserialize(self, resp):
        """
        Given a ``HttpResponse`` coming back from using the ``client``, this method
        checks the ``Content-Type`` header & attempts to deserialize the data based on
        that.

        It returns a Python datastructure (typically a ``dict``) of the serialized data.
        """
        # Override to allow for use of the StreamingHttpResponse or the HttpResponse
        return self.serializer.deserialize(self.get_content(resp), format=resp['Content-Type'])

    def serialize(self, data, format='application/json'):
        """
        Given a Python datastructure (typically a ``dict``) & a desired content-type,
        this method will return a serialized string of that data.
        """
        return self.serializer.serialize(data, format=format)


def setUpModule():
    logger.info(_('=== setup Module'))

    # FIXME: running the bootstrap method as a test suite setup:
    # TODO: create a local TestRunner,TestSuite,
    # so that this can be run before the suite
    testContext = IResourceTestCase(methodName='_bootstrap_init_files')
    testContext.setUp()
    testContext._bootstrap_init_files()

    logger.info(_('=== setup Module done'))

def tearDownModule():
    logger.info(_('=== teardown Module'))


@unittest.skip("Only run this test if testing the api initialization")        
class TestApiInit(IResourceTestCase):
    
    def setUp(self):
        # Create a user.
        self.resource_uri = BASE_URI + '/metahash'
        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        self.csv_serializer=CSVSerializer() 
        
        self.serializer = LimsSerializer()
        self.api_client = TestApiClient(serializer=self.serializer)
        

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
                (resp.status_code, resp, 'resource_uri', self.resource_uri,
                    'item', item))
            
        logger.debug('created items, now get them')
        resp = self.api_client.get(
            self.resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.debug(_('--------resp to get:', resp, resp.status_code))
        new_obj = self.deserialize(resp)
        logger.debug(_('deserialized object:', json.dumps(new_obj)))
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        self.assertEqual(len(new_obj['objects']), 4, (new_obj))
        
        for inputobj in bootstrap_items:
            logger.debug(_('testing:', inputobj))
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, ('not found', inputobj, outputobj ) )
        
        logger.debug('================ test0_bootstrap_metahash done ========== ')
    
    def test2_api_init(self):
        
        logger.debug('***================ reports test2_api_init =============== ')
        serializer=CSVSerializer() 
        # todo: doesn't work for post, see TestApiClient.post() method, it is 
        # incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=serializer) 
        
        filename = os.path.join(self.directory,'api_init_actions.csv')
        with open(filename) as input_file:
            # NOTE / TODO: we have to deserialize the input, because the TP test method 
            # will expect a python data object, which it will serialize!
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
                    logger.debug(_('+++++++++++processing file', filename))
                    with open(filename) as data_file:
                        # NOTE / TODO: we have to deserialize the input, because the TP test method 
                        # will expect a python data object, which it will serialize!
                        input_data = serializer.from_csv(data_file.read())
                        
                        if command == 'put':
                            resp = testApiClient.put(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
                            logger.debug(str(('action: ', json.dumps(action), 
                                              'response: ' , resp.status_code)))
                            self.assertTrue(
                                resp.status_code in [200], (resp.status_code, resp))
                            
                            # now see if we can get these objects back
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), data={ 'limit': 999 })
                            self.assertTrue(
                                resp.status_code in [200], (resp.status_code))
                            new_obj = self.deserialize(resp)
                            result, msgs = find_all_obj_in_list(
                                input_data['objects'], new_obj['objects'], 
                                excludes=search_excludes)
                            self.assertTrue(
                                result, (command, 'input file', filename, 
                                             msgs, new_obj['objects']) )
                        
                        elif command == 'patch':
                            resp = testApiClient.patch(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
                            self.assertHttpAccepted(resp)
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), data={ 'limit': 999 } )
                            self.assertTrue(
                                resp.status_code in [200], (resp.status_code))
                            new_obj = self.deserialize(resp)
                            with open(filename) as f2:
                                input_data2 = serializer.from_csv(f2.read())
                                result, msgs = find_all_obj_in_list(
                                    input_data2['objects'], new_obj['objects'], 
                                    excludes=search_excludes)
                                self.assertTrue(
                                    result, ( command, 'input file', filename, msgs ) )
                        
                        elif command == 'post':
                            self.fail((
                                'resource entry: ' + json.dumps(action) + '; '
                                'cannot POST multiple objects to tastypie; '
                                'therefore the "post" command is invalid with '
                                'the initialization scripts'))
                        else:
                            self.fail('unknown command: ' + command + ', ' + json.dumps(action))
        

class UserUsergroupSharedTest(object):
            
    def test1_create_user_with_permissions(self, test_log=False):
        logger.info(_('==== test1_create_user_with_permissions =====', test_log))
        
        filename = os.path.join(self.directory,'test_data/users1.csv')
        from datetime import datetime
        data_for_get = { HEADER_APILOG_COMMENT: 
            'patch_test: file: %s, %s' % (filename, datetime.now().isoformat() ) }
        input_data,output_data = self._put_test(
            'user', filename, id_keys_to_check=['username'],
            data_for_get=data_for_get)
        
        # test the logs
        # FIXME: create separate tests for the apilogs (cannot ensure state, 
        # as there are multiple callers here )
        if test_log:
            resource_uri = BASE_URI + '/apilog' #?ref_resource_name=record'
            logger.info(_('get', resource_uri, data_for_get))
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 'limit': 999, 'ref_resource_name': 'user',
                    'comment__eq': data_for_get[HEADER_APILOG_COMMENT]}) #,
#                     'api_action__in': 'CREATE' })
            self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
            new_obj = self.deserialize(resp)
            objects = new_obj['objects']
            logger.debug(_('===apilogs:', json.dumps(new_obj)))
            # look for 6 "CREATE" logs
            self.assertEqual( len(objects), 6, 
                str((6,len(objects), 'wrong # of api logs', objects)))
            for obj in objects:
                self.assertTrue(obj['api_action'] == API_ACTION_CREATE, 
                    ('action should be create', obj))
                self.assertTrue( obj['comment'] == data_for_get[HEADER_APILOG_COMMENT], 
                    ('comment not set', data_for_get[HEADER_APILOG_COMMENT],obj))
        
        logger.debug(_('==== test1_create_user_with_permissions ====='))
                    
class UserResource(IResourceTestCase, UserUsergroupSharedTest):

    def setUp(self):
        super(UserResource, self).setUp()
    
    def test0_create_user(self):
        logger.debug(_('==== test_create_user ====='))
        
        # the simplest of tests, create some simple users
        self.bootstrap_items = { 'objects': [   
            {
                'username': 'st1',
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
                'username': 'ut1',
                'first_name': 'userx',
                'last_name': 'Tester1',    
                'email': 'user.tester1@slimstest.com',    
            },
        ]}
        try:       
            uri = BASE_URI + '/user'
            resp = self.api_client.put(uri, 
                format='json', data=self.bootstrap_items, authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [200,201,202], 
                (resp.status_code, self.serialize(resp)))
            #                 self.assertHttpCreated(resp)
        except Exception, e:
            logger.error(_('on creating', self.bootstrap_items, '==ex==', e))
            raise

        logger.debug('created items, now get them')
        resp = self.api_client.get(uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.debug(_('--------resp to get:', resp.status_code))
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        self.assertEqual(len(new_obj['objects']), 3, (new_obj))
        
        for i,item in enumerate(self.bootstrap_items['objects']):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(
                result, ('bootstrap item not found', item, new_obj['objects']))
            logger.debug(_('item found', obj))

        logger.debug(_('==== test_create_user done ====='))

    def test2_patch_user_permissions(self,test_log=False):
        '''
        @param test_log set to false so that only the first time tests the logs
        '''
        
        logger.info(_('==== test2_patch_user_permissions =====', test_log))
        self.test1_create_user_with_permissions(test_log=test_log)
        
        logger.debug(_('==== test2_patch_user_permissions start ====='))
        filename = os.path.join(self.directory,'test_data/users2_patch.csv')
        data_for_get = { HEADER_APILOG_COMMENT: 'patch_test: %s' % filename }

        input_data,output_data = self._patch_test(
            'user', filename, data_for_get=data_for_get)
        
        # # test the logs
        # resource_uri = BASE_URI + '/apilog' #?ref_resource_name=record'
        # logger.debug(_('get', resource_uri))
        # resp = self.api_client.get(
        #     resource_uri, format='json', 
        #     authentication=self.get_credentials(), 
        #     data={ 'limit': 999, 'ref_resource_name': 'user',
        #         'uri__in': ','.join([x['resource_uri'][x['resource_uri'].find('user'):] for x in input_data]),
        #         'api_action__eq': 'PATCH' })
        # self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        # new_obj = self.deserialize(resp)
        # objects = new_obj['objects']
        # logger.debug(_('===apilogs:', json.dumps(new_obj)))
        # # look for 3 logs: 3 for patch list
        # self.assertEqual( len(objects), 3, 
        #     str((len(objects), 'wrong # of api logs', objects)))
        # for obj in objects:
        #     self.assertTrue('permissions' in json.dumps(obj['diff_keys']), 
        #         ('no "permissions" key in the diff keys', obj))
        #     self.assertTrue( obj['comment'] == data_for_get[HEADER_APILOG_COMMENT], 
        #         ('comment not set', data_for_get[HEADER_APILOG_COMMENT],obj))
        
        logger.debug(_('==== test2_patch_user_permissions done ====='))
        
    def test3_user_read_permissions(self):
        '''
        Try to do something we don't have permissions for - 
        read a resource (metahash) 
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        '''
        
        logger.debug(_('==== test3_user_read_permissions ====='))
        self.test2_patch_user_permissions(test_log=False)
                
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # Try to do some unauthorized actions

        resource_uri = BASE_URI + '/metahash'
        resp = self.api_client.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [403], (resp.status_code, resp))
        
        # Now add the needed permission
        
        user_patch = {
            'resource_uri': 'user/' + username,
            'permissions': ['resource/metahash/read'] };

        uri = BASE_URI + '/user/' + username
        logger.debug(_('=====now add the permission needed to this user:', user_patch, uri))
        resp = self.api_client.patch( uri, 
                    format='json', data=user_patch, 
                    authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        (resp.status_code, self.serialize(resp)))  
        
        # now try again as the updated user:
        
        resp = self.api_client.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [200], 
                        (resp.status_code, self.serialize(resp)))
       
        logger.debug(_('==== test3_user_read_permissions done ====='))

    def test4_user_write_permissions(self):
        '''
        Try to do something we don't have permissions for - 
        write a resource (user),
        then give that permission to the user, and try again 
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        '''
        logger.debug(_('==== test4_user_write_permissions ====='))

        self.test2_patch_user_permissions(test_log=False)
                
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
        
        resp = self.api_client.patch(
            resource_uri, format='json', data=json_data, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [403], (resp.status_code, self.serialize(resp)))
        
        # Now add the needed permission
        
        user_patch = {
            'resource_uri': 'user/' + username,
            'permissions': ['resource/usergroup/write'] };

        logger.debug(_('now add the permission needed to this user:', user_patch))
        uri = BASE_URI + '/user/' + username
        resp = self.api_client.patch( uri, 
                    format='json', data=user_patch, 
                    authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        (resp.status_code, self.serialize(resp)))  
        
        # now try again as the updated user:
        
        resp = self.api_client.patch(
            resource_uri, format='json', data=json_data, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [202, 204], 
                        (resp.status_code, self.serialize(resp)))

        #         # is it set?
        #         resp = self.api_client.get(
        #                 uri, format='json', authentication=self.get_credentials())
        #         logger.warn(str(('response: ' , self.deserialize(resp) )))

        logger.debug(_('==== test4_user_write_permissions done ====='))
        
class UserGroupResource(IResourceTestCase, UserUsergroupSharedTest):

    def setUp(self):
        super(UserGroupResource, self).setUp()

    def test2_create_usergroup_with_permissions(self):
        logger.debug(_('==== test2_usergroup ====='))
        #create users
        self.test1_create_user_with_permissions(test_log=False)
        
        filename = os.path.join(self.directory,'test_data/usergroups1.csv')
        # note: excluding sub_groups here because the one sub_group is set when 
        # "testGroupX" sets super_groups=['testGroup3']; and thereby testGroup3 
        # gets sub_groups=['testGroupX']; even though that's not in the input file.
        self._put_test('usergroup', filename, id_keys_to_check=['name'],
            keys_not_to_check=['sub_groups'])

        logger.debug(_('==== test2_usergroup done ====='))

    def test3_patch_users_groups(self):
        logger.debug(_('==== test3_patch_users_groups ====='))
        self.test2_create_usergroup_with_permissions()
        
        logger.debug(_('==== test3_patch_users_groups start ====='))
        filename = os.path.join(self.directory,'test_data/users3_groups_patch.csv')
        # don't check permissions, because the patch file is setting groups with permissions too
        self._patch_test('user', filename, id_keys_to_check=['username'],
            keys_not_to_check=['permissions']) 
      
        logger.debug(_('==== test3_patch_users_groups done ====='))


    def test4_user_group_permissions(self):
        '''
        Test of a group permission -
        first test that the user in the group has the permission,
        then remove this user from this group and try again
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.debug(_('==== test4_user_group_permissions ====='))
        
        self.test2_create_usergroup_with_permissions()
                        
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # 1 read test - should have permission through group
        resource_uri = BASE_URI + '/metahash'
        resp = self.api_client.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(resp.status_code in [200], (resp.status_code, self.serialize(resp)))
        
        # now patch this user's usergroups, removing the user from the group 'testgroup1'
        # which will remove the permissions as well 
        user_patch = {
            'usergroups': ['testGroup3'] };

        logger.debug(_('now reset this users groups and remove testGroup1:', user_patch))
        uri = BASE_URI + '/user' + '/' + username
        resp = self.api_client.patch(uri, format='json', data=user_patch, 
                                     authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        (resp.status_code, self.serialize(resp)))  
        
        # check the user settings/groups
        uri = BASE_URI + '/user/' + username
        resp = self.api_client.get(uri,format='json', data={'includes': '*'}
            , authentication=self.get_credentials() )
        logger.debug(_('--------resp to get:', resp.status_code))
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, new_obj))
        logger.info(_('user', new_obj))
        self.assertTrue(new_obj['usergroups'] == ['testGroup3'],
            ('wrong usergroups', new_obj))
        
        # now try again as the updated user:
        resp = self.api_client.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(resp.status_code in [403], 
                        (resp.status_code, self.serialize(resp)))
       
        logger.debug(_('==== test4_user_group_permissions done ====='))


    def test5_usergroup_can_contain_group_permissions(self):
        '''
        Test of an inherited group permission  -
        first test that the user in the group doesn't have the permission,
        then add the user's usergroup to a group with the permission and try again
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.debug(_('==== test5_usergroup_can_contain_group_permissions ====='))
        
        self.test2_create_usergroup_with_permissions()
                        
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # 1 read test - user, user's group don't have the permission
        resource_uri = BASE_URI + '/vocabularies'
        resp = self.api_client.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(resp.status_code in [403], 
            (resp.status_code, self.serialize(resp)))
        
        # now create a new group, with previous user's group as a member,
        # then add permissions to this new group to read (vocabularies)
        # note: double nest the groups also as a test
        usergroup_patch = { 'objects': [
            {
            'name': 'testGroup5',
            'super_groups': ['testGroup3'] },
            {
            'name': 'testGroup6',
            'users': ['sde4'],
            'super_groups': ['testGroup5'] },
        ]}
        
        logger.debug(_('now set the new group:', usergroup_patch))
        uri = BASE_URI + '/usergroup'
        
        resp = self.api_client.patch(uri, format='json', 
            data=usergroup_patch, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
            (resp.status_code, self.serialize(resp)))  
        
        # is it set?
        resp = self.api_client.get(
                uri, format='json', authentication=self.get_credentials())
        new_obj = self.deserialize(resp)
        result, outputobj = find_all_obj_in_list(
            usergroup_patch['objects'],new_obj['objects'],
            id_keys_to_check=['name']) #, excludes=keys_not_to_check )
        self.assertTrue(
            result, 
            ('not found', outputobj,'=== objects returned ===', 
                 new_obj['objects'] ) ) 
        
        # is it set-2, does group inherit the permissions?
        for obj in new_obj['objects']:
            if obj['name'] == 'testGroup6':
                testGroup6 = obj
            if obj['name'] == 'testGroup3':
                testGroup3 = obj
                
        self.assertTrue(testGroup6 and testGroup3)
        for permission in testGroup3['all_permissions']:
            logger.debug(_('find permission', permission))
            self.assertTrue(permission in testGroup6['all_permissions'], 
                ('could not find permission', permission, 
                     'in testGroup6 permissions', testGroup6['all_permissions']))
        
        # 2 read test - user has permissions through inherited permissions,
        resource_uri = BASE_URI + '/vocabularies'
        resp = self.api_client.get(
            resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, self.serialize(resp)))

        logger.debug(_('==== Done: test5_usergroup_can_contain_group_permissions ====='))

    
    def test6_usergroup_can_contain_group_users(self):
        '''
        Test that group "contains" the subgroup users (through "all_users")
        first test that the group doesn't have the user, then add the subgroup
        with the user, test that the group has the user in "all_users
        - done in prev test: create a new user (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.debug(_('==== test6_usergroup_can_contain_group_users ====='))
        
        self.test2_create_usergroup_with_permissions()
                        
        # now create a new group, with a previous group as a sub_group
        usergroup_patch = { 'objects': [
            {
            'resource_uri': 'usergroup/testGroup5',
            'name': 'testGroup5',
            'sub_groups': ['testGroup2'] },
            {
            'resource_uri': 'usergroup/testGroup6',
            'name': 'testGroup6',
            'sub_groups': ['testGroup5'] },
        ]}
        
        logger.debug(_('now set the new groups:', usergroup_patch))
#         uri = self.resource_uri + '/'
        uri = BASE_URI + '/usergroup'
        resp = self.api_client.patch(uri, format='json', 
            data=usergroup_patch, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202, 204], 
                        (resp.status_code, self.serialize(resp)))  
        
        # is it set?
        data_for_get={ 'includes': '*'}
        resp = self.api_client.get(
                uri + '/testGroup6', 
                format='json', 
                authentication=self.get_credentials(),
                data=data_for_get)
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, self.serialize(resp)))
        new_obj = self.deserialize(resp)
        
        self.assertTrue(new_obj['all_users'])
        self.assertTrue('sde4' in new_obj['all_users'],
            ('could not find user', 'sde4', new_obj['all_users']))
        
        # TODO: could also test that testGroup2 now has super_group=testGroup5
        
        logger.debug(_('==== Done: test6_usergroup_can_contain_group_users ====='))
        
@unittest.skip('TODO: remove or restore')
class RecordResource(IResourceTestCase):
    
    def setUp(self):
        super(RecordResource, self).setUp()
        
    def test0_bootstrap_metahash(self):
        
        logger.debug('================ reports test0_bootstrap_metahash =============== ')
        
        # TODO: fields for the "metahash:fields.metahash" definitions input file
        bootstrap_items = [   
            {
                'key': 'linked_field_type',
                'scope': 'fields.metahash',
                'json_field_type': 'fields.CharField',
                'ordinal': 4,
                'description': 'The tastypie field used to serialize'    
            },
            {
                'key': 'linked_field_module',
                'scope': 'fields.metahash',
                'json_field_type': 'fields.CharField',
                'ordinal': 6,
                'description': 'The table model used to hold the linked field',
                'comment': 'Use the full Python style path, i.e. <module.module.Class>'
            },
            {
                'key': 'linked_field_value_field',
                'scope': 'fields.metahash',
                'json_field_type': 'fields.CharField',
                'ordinal': 7,
                'description': 'The field name in the linked field module that holds the value'
            },
            {
                'key': 'linked_field_parent',
                'scope': 'fields.metahash',
                'json_field_type': 'fields.CharField',
                'ordinal': 6,
                'description': 'The name of the field linking back to the parent table/resource'
            },
        ]
        resource_uri = BASE_URI + '/metahash'
         
        resp = self.api_client.patch(
            resource_uri, format='json', data={ 'objects': bootstrap_items }, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202,204], 
            (resp.status_code, resp, 'resource_uri', self.resource_uri,
                'item', bootstrap_items))
             
        logger.debug('created items, now get them')
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.debug(_('--------resp to get:', resp, resp.status_code))
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
         
        for inputobj in bootstrap_items:
            logger.debug(_('testing:', inputobj))
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, ('not found', inputobj, outputobj ) )

        # TODO: add this to the "metahash:fields.resource" file 
        # for the complex linked table instance (SMR, RNAi)
        # Can use either the resource-level "linked_table_module" or the 
        # field level "linked_field_module"
        bootstrap_items = [   
            {
                'key': 'linked_table_module',
                'scope': 'fields.resource',
                'ordinal': 4,    
                'json_field_type': 'fields.CharField'    
            },
        ]
        resource_uri = BASE_URI + '/metahash'
        resp = self.api_client.patch(
            resource_uri, format='json', data={ 'objects': bootstrap_items }, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202,204], 
            (resp.status_code, resp, 'resource_uri', self.resource_uri,
                'item', bootstrap_items))

        # create a resource for the RecordResource
        resource_data = [
            {   'key':  'record',
                'scope': 'resource',
                'ordinal': '0',
                'api_name': 'reports',
                # TODO: use either this or the "linked_field_module" value on the field defs
                # TODO: this is not used at the moment
                'linked_table_module': 'reports.models.RecordValueComplex', 
                'id_attribute': ['id'],
            },
        ]
        resource_uri = BASE_URI + '/resource'
        
        for item in resource_data:         
            resp = self.api_client.post(
                resource_uri, format='json', data=item, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], 
                (resp.status_code, resp, 'resource_uri', self.resource_uri,
                    'item', item))

        logger.debug('created items, now get them')
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.debug(_('--------resp to get:', resp, resp.status_code))
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))

        for inputobj in resource_data:
            logger.debug(_('testing:', inputobj))
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, ('not found', inputobj, outputobj ) )
 
        logger.debug('================ test0_bootstrap_metahash done ========== ')
    
    def test1_persist_to_store(self):
    
        logger.debug('================ setup =============== ')
        self.test0_bootstrap_metahash()
        
        logger.debug('================ reports test1_persist_to_store =============== ')
        
        # simplest test, create the data manually, then try to retrieve it.
        
        # create field definitions in the metahash for the record resource
        resource_uri = BASE_URI + '/metahash'
        record_fields = [
            {   'key':  'id',
                'scope': 'fields.record',
                'ordinal': '0',
            },
            {   'key':  'scope',
                'scope': 'fields.record',
                'ordinal': '1',
                'filtering': 'true',
            },
            {   'key':  'base_value1',
                'scope': 'fields.record',
                'ordinal': '2',
                'filtering': 'true',
            },
            {   'key':  'field1',
                'scope': 'fields.record',
                'ordinal': '2',
                'linked_field_type': 'fields.CharField',
                'linked_field_module': 'reports.models.RecordValue',
                'linked_field_value_field': 'value',
                'linked_field_parent': 'parent',
                'linked_field_meta_field': 'field_meta',
            },
            {   'key':  'field2',
                'scope': 'fields.record',
                'ordinal': '2',
                'linked_field_type': 'fields.CharField',
                'linked_field_module': 'reports.models.RecordValue',
                'linked_field_value_field': 'value',
                'linked_field_parent': 'parent',
                'linked_field_meta_field': 'field_meta',
            },
            {   'key':  'field3',
                'scope': 'fields.record',
                'ordinal': '3',
                'linked_field_type': 'fields.ListField',
                'linked_field_module': 'reports.models.RecordMultiValue',
                'linked_field_value_field': 'value',
                'linked_field_parent': 'parent',
                'linked_field_meta_field': 'field_meta',
            },
            # complex type storing both field4 and field5
            {   'key':  'field4',
                'scope': 'fields.record',
                'ordinal': '4',
                'linked_field_type': 'fields.CharField',
                # Omit the linked_field_modele, so that this will use the 
                # "linked_table_module" value from the resource definition
                # 'linked_field_module': 'reports.models.RecordMultiValue',
                # also omit the 'linked_field_meta_field': 'field_meta', because
                # the complex type will have only one entry per parent
                'linked_field_value_field': 'value1',
                'linked_field_parent': 'parent'
            },
            {   'key':  'field5',
                'scope': 'fields.record',
                'ordinal': '5',
                'linked_field_type': 'fields.CharField',
                # Omit the linked_field_module, so that this will use the 
                # "linked_table_module" value from the resource definition
                # 'linked_field_module': 'reports.models.RecordMultiValue',
                # also omit the 'linked_field_meta_field': 'field_meta', because
                # the complex type will have only one entry per parent
                'linked_field_value_field': 'value2',
                'linked_field_parent': 'parent'
            },
        ]

        resp = self.api_client.patch(
            resource_uri, format='json', data={ 'objects': record_fields}, 
            authentication=self.get_credentials())
        logger.debug(str(('===resp',self.deserialize(resp) )))
        self.assertTrue(resp.status_code in [202,204], 
            (resp.status_code, resp, 'resource_uri', self.resource_uri,
                'item', record_fields))
            #             self.assertHttpCreated(resp)
            
        logger.debug('created items, now get them')
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), 
                data={ 'limit': 999, 'scope':'fields.record' })
        logger.debug(_('--------resp to get:', resp, resp.status_code))
        new_obj = self.deserialize(resp)
        logger.debug(_('deserialized object:', json.dumps(new_obj)))
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        self.assertEqual(len(new_obj['objects']), len(record_fields), 
            str((len(new_obj['objects']), new_obj)))
        logger.debug(_('=== returned objs', new_obj['objects']))
        
        for inputobj in record_fields:
            logger.debug(_('=====testing:', inputobj))
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, ('not found', inputobj, outputobj ) )
        
        logger.debug(_('==== now create datapoints in the record table'))
        
        datapoints = [
            {   'scope': 'record',
                'base_value1': 'base value 1 1',
                'field1': 'test value to store/retrieve',
                'field2': '2nd field test value to store/retrieve',
                'field3': ['valueC','ValueE','ValueA'],
                'field4': '1st recordvaluecomplex',
                'field5': '2nd recordvaluecomplex',
            },
            {   'scope': 'record',  # vanilla 'record'
                'field1': 'test value to store/retrieve',
                'field2': '2b field test value to store/retrieve',
                'field3': ['valueA','ValueB','ValueC'],
                'field4': '1st recordvaluecomplex 2',
                'field5': '2nd recordvaluecomplex 2',
            }
        ]
        self.datapoints = datapoints
        resource_uri = BASE_URI + '/record'
        
        for item in datapoints:         
            resp = self.api_client.post(
                resource_uri, format='json', data=item, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], 
                (resp.status_code, resp, 'resource_uri', self.resource_uri,
                    'item', item))

        logger.debug('=== get the datapoints created in the record table')
        
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999, }) #'scope':'record' })
        
        logger.debug(_('--------resp to get:', resp, resp.status_code))
        
        new_obj = self.deserialize(resp)
        
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        self.assertEqual(len(new_obj['objects']), len(datapoints))
        
        for inputobj in datapoints:
            logger.debug(_('=====testing:', inputobj))
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, ('not found', inputobj, outputobj ) )
            logger.debug(_('===found: ', outputobj))
        
        logger.debug('================ done: reports test1_persist_to_store =============== ')

        #### NOTE: combining the two tests into one:
        #### although the django.test.TransactionTestCase isolates the tests in transactions,
        #### it does not rollback the sequence values; so the ID for the first and 
        #### the second record are "1" and "2" on the first run, if test2 re-runs
        #### test1, then the ID's for the second run will be "3" and "4"
        #### since we don't know if test1 is run before test2, we'll need to find
        #### the records before we can update them, unless we create a natural key for them.

        logger.debug('================ reports test2_update =============== ')
            
        logger.debug(_('==== now update datapoints in the record table'))
        
        datapoints = [
            {   'resource_uri':'record/2',
                'base_value1': 'updated base value 1 1',
                'field1': 'xxx updated',
                'field3': ['valueA','ValueE','ValueC'],
                'field4': 'xxx updated field 4 1',
            },
            {   'resource_uri':'record/1',
                'field2': 'xxx updated',
                'field3': ['valueA','ValueB','ValueC'],
                'field5': 'xxx updated field 5 2',
            }
        ]
        
        resource_uri = BASE_URI + '/record'
        resp = self.api_client.patch(
            resource_uri, format='json', data={ 'objects': datapoints }, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [202,204], 
            (resp.status_code, resp, 'resource_uri', self.resource_uri))

        logger.debug('=== get the datapoints patched in the record table')
        
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999, }) #'scope':'record' })
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        self.assertEqual(len(new_obj['objects']), len(datapoints), 
            str((len(new_obj['objects']), len(datapoints),new_obj)))
        
        for inputobj in datapoints:
            logger.debug(_('=====testing:', inputobj))
            result, outputobj = find_obj_in_list(inputobj,new_obj['objects'])
            self.assertTrue(result, (outputobj,new_obj['objects'] ) )
            logger.debug(_('=== found: ', outputobj))
            self.assertTrue(inputobj['resource_uri'] in outputobj['resource_uri'],
                (inputobj['resource_uri'] ,outputobj['resource_uri']))
        
        # test the logs
        resource_uri = BASE_URI + '/apilog' #?ref_resource_name=record'
        logger.debug(_('get', resource_uri))
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), 
            data={ 'limit': 999, 'ref_resource_name': 'record' })
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        new_obj = self.deserialize(resp)
        logger.debug(_('===apilogs:', json.dumps(new_obj)))
        
        # look for 5 logs; 2 for create, two for update, one for patch list
        self.assertEqual( len(new_obj['objects']), len(datapoints)*2+1, 
            str((len(new_obj['objects']), len(datapoints)*2+1,new_obj)))
        
        logger.debug('================ done: reports test2_update =============== ')
            
