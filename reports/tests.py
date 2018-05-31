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

IccblTestRunner: set in the settings.py:
run as:
./manage.py test db.tests.TestName --settings=lims.settings_testing --keepdb
Options:
--keepdb:
    if not set, the database is always reinitialized
    if set, then the database is not destroyed after tests are run
--reinit_metahash
    if not set, then the metahash that would be left over from running with 
    --keepdb will be reused
    otherwise, all of the "_bootstrap_init_files" procedures are run


To Run tests without a database:
Example:
./manage.py test reports.SDFSerializerTest.test2_clean_data_sdf \
     --settings=lims.settings_testing_debug --verbosity=2 \
     --testrunner=reports.tests.NoDbTestRunner
NOTE: when using this testrunner, the class finder is different; note that the 
path excludes the module, so 
"reports.SDFSerializerTest.test2_clean_data_sdf"
not 
"reports.test.SDFSerializerTest.test2_clean_data_sdf"

NOTE: to run using sqlite in memory, use a test_settings file with the database as:
DATABASES['default'] = {'ENGINE': 'django.db.backends.sqlite3'}
and run like:
$ ./manage.py test --settings=lims.test_settings

"""

from __future__ import unicode_literals

from decimal import Decimal
import json
import logging
import os
import re
import sys
import unittest
import urlparse

import dateutil.parser
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.db.utils import ProgrammingError
from django.http.response import StreamingHttpResponse
from django.test import TestCase
from django.test.client import Client, FakePayload
from django.test.runner import DiscoverRunner
from django.test.testcases import SimpleTestCase

from reports import HEADER_APILOG_COMMENT
from reports.api import compare_dicts, API_RESULT_DATA, API_RESULT_META
from reports.dump_obj import dumpObj
from reports.models import MetaHash, UserGroup, \
    UserProfile, ApiLog, Permission, Job
import reports.schema as SCHEMA
from reports.serialize import parse_val, JSON_MIMETYPE
import reports.serialize.csvutils as csvutils
from reports.serialize.sdfutils import MOLDATAKEY
from reports.serializers import CSVSerializer, SDFSerializer, \
    LimsSerializer, XLSSerializer
import reports.utils.background_processor
import reports.utils.log_utils


logger = logging.getLogger(__name__)

BASE_URI = '/reports/api/v1'

# NOTE: the simulated django client requests expect the HTTP Header "Accept" to 
# be stored in the variable "HTTP_ACCEPT"
DJANGO_ACCEPT_PARAM = 'HTTP_ACCEPT'

# Required for non-staff users to log in
settings.IS_PRODUCTION_READY = True

import reports; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(reports.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(reports.__path__))

DEBUG = False

def find_in_dict(key,content):
    ''' 
    @return a value nested in a key of a dictionary - tree
    '''
    if isinstance(content, dict):
        if key in content:
            return content[key]
        else:
            for value in content.values():
                val = find_in_dict(key,value)
                if val:
                    return val
    return None

def is_boolean(field):
    if type(field) == bool: 
        return True
    
    if isinstance(field, basestring):
        val = field.lower()
        if val=='true' or val=='false':
            return True
    return False
    

def numerical_equivalency(val1, val2):
    try:
        if DEBUG:
            logger.info('numerical equivalency: %r to %r, %r to %r', 
                val1, val2, Decimal(val1), Decimal(val2))
        if Decimal(val1) == Decimal(val2):
            return True
    except:
        pass
    return False

####
# Equivocal, and other equivalency methods are for end-to-end testing, 
# using string equals, of the API with values read from csv or json.
# Values are submitted, then read back as json from the API.
####

def decode_from_utf8(val):
    '''
    Convert all values to unicode for comparison
    - NOTE: csv package does not support unicode, it returns bytes, so decode
    byte strings for comparison.
    '''
    try:
        if isinstance(val, str):
            # if string, then decode back to unicode, assuming utf-8 encoding
            val = val.decode('utf-8')
        return val
    except Exception, e:
        logger.exception('decode: %r, %r', val, type(val))
        raise
    
def equivocal(val1, val2, skip_null_values=False, key=None):
    '''
    Test for equivalence between submitted and returned (csv) values:
    Compares string,number,boolean,or list where:
    - either object has been converted to a string,
    - with lists, ordering has changed, or members have been converted to 
    their string representation.
    NOTE: skip_null_values added to address fields with parent "ref" set
    '''
    if DEBUG:
        logger.info('equivocal: %r, %r', val1, val2)

    if val1 is None or val1 == '':
        if skip_null_values is True:
            logger.debug('key: %r, skip nulls: %r',key, skip_null_values)
            return True, ('skip_null_values')
    if val1 == val2:
        if DEBUG:
            logger.info('1equivocal: %r, %r', val1, val2)
        return True, ('val1', val1, 'val2', val2 )
    
    if (is_boolean(val1) is True or is_boolean(val2) is True ):
        if DEBUG:
            logger.info('boolean equivocal: %r, %r', val1, val2)
        if (parse_val(val1,'testval1','boolean')
            == parse_val(val2,'testval2','boolean')):
            if DEBUG:
                logger.info('boolean equivalent: %r, %r',
                    parse_val(val1,'testval1','boolean'),
                    parse_val(val2,'testval2','boolean') )
            return True, ('val1', val1, 'val2', val2 )
        else:
            return False, ('boolean not equivalent: val1', val1, 'val2', val2 )
        
    if ( isinstance(val1, (int, long, float, complex, Decimal))
        or isinstance(val2, (int, long, float, complex, Decimal))):
        if DEBUG:
            logger.info('numerical equivocal: %r, %r', val1, val2)
        if numerical_equivalency(val1, val2):
            return True, ('val1', val1, 'val2', val2 )

        val1 = str(val1)
        val2 = str(val2)
        if val1 != val2:
            if DEBUG:
                logger.info('fail final numerical test: %r to %r', val1, val2)
            return False, ('val1', val1, 'val2', val2 )
        
    
    if isinstance(val1, basestring):
        if DEBUG:
            logger.info('string equivocal: %r, %r', val1, val2)

        # TODO: rework API: equates empty string to "None"
        if not val1:
            if not val2:
                return True, ('val1', val1, 'val2', val2 )
            else:
                return False, ('val1', val1, 'val2', val2 )
        if not val2:
            return False, ('val2 is empty', 'val1', val1 )
            
        # if string value, then decode back to unicode
        if DEBUG:
            logger.info('val1 %r, decode: %r', val1, decode_from_utf8(val1))
        val1 = decode_from_utf8(val1)
        if hasattr(val2, "__getitem__") or hasattr(val2, "__iter__"): 
            # allow single item lists to be equal to string
            if len(val2) == 1:
                val2 = val2[0]
        if DEBUG:
            logger.info('val2 %r, decode: %r', val2, decode_from_utf8(val2))
        val2 = decode_from_utf8(val2)
        if DEBUG:
            logger.info('string equivalency: %r to %r', val1, val2)

        if val1 != val2:
            if ((is_boolean(val1) or is_boolean(val2) ) 
                and parse_val(val1,'testval1','boolean')
                    == parse_val(val2,'testval2','boolean')):
                if DEBUG:
                    logger.info('boolean true')
                return True, ('val1', val1, 'val2', val2 )
            elif numerical_equivalency(val1,val2):
                return True, ('val1', val1, 'val2', val2 )
                
            return False, ('val1', val1, 'val2', val2 )
        else:
            if DEBUG:
                logger.info('String true: %r: %r', val1, val2)
            return True, ('val1', val1, 'val2', val2    )
    else: # better be a list
        if DEBUG:
            logger.info('list equivocal: %r, %r', val1, val2)
        # TODO: rework API: equates empty list to "None"
        if not val1:
            if not val2:
                return True, ('val1', val1, 'val2', val2 )
            else:
                if skip_null_values is True:
                    logger.info('skip None: %r, %r', val1, val2)
                    return True, ('skip_null_values is True',)
                else:
                    return False, ('val1', val1, 'val2', val2 )
        if not val2:
            return False, ('val2 is empty', 'val1', val1 )
        if not isinstance(val1, list) and isinstance(val2, list):
            return False, (
                'Must be a list if not a string', 'val1', val1, 'val2', val2)
        for v in val1:
            if not v: 
                if val2 and len(val2) > 0:
                    return False, ('val1', val1, 'list val2 not empty', val2 )
            if v not in val2:
                v = decode_from_utf8(v)
                if v not in [decode_from_utf8(v2) for v2 in val2]:
                    return False, ('val1', val1, 'val2', val2 )
        if DEBUG:
            logger.info('list equivocal True: %r, %r', val1, val2)
        return True, ('val1', val1, 'val2', val2 )
    
    logger.error(
        'equivocal: val not recognized (simple, boolean, numerical, string or list)'
        ', val1: %r:%r, val2: %r:%r', type(val1),val1, type(val2), val2)
    return False, ('val1', val1, 'val2', val2)
    
def assert_obj1_to_obj2( obj1, obj2, excludes=['resource_uri'], skip_null_values=False):
    '''
    For testing equality of the (CSV) input to API returned values
    @param obj1 input
    @param obj2 output 
    NOTE: skip_null_values added to address fields with parent "ref" set
    '''
    if obj1 is None:
        return False, ('obj1 is None')
    if obj2 is None:
        return False, ('obj2 is None')
    
    original_keys = set([k for k in obj1.keys() if k])
    original_keys = original_keys.difference(excludes)
    updated_keys = set(obj2.keys())
    intersect_keys = original_keys.intersection(updated_keys)
    if intersect_keys != original_keys:
        logger.info('missing: %r', sorted(original_keys-intersect_keys))
        return False, ('keys missing:', 
            original_keys-intersect_keys, 
            'original_keys', sorted(original_keys), 
            'updated_keys', sorted(updated_keys))

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
                        return False, (
                            'resource uri not equal', val1, val2, obj2)
                    logger.warn(
                        'imprecise uri matching, equivocating %r to %r:', 
                        val1, val2)
    
    keys_to_search = original_keys
    logger.debug('keys to test: %r: %r', keys_to_search, skip_null_values)
    for key in keys_to_search:
        result, msgs =  equivocal(obj1[key], obj2[key], skip_null_values=skip_null_values, key=key)
        if not result:
            return False, ('msgs', msgs,'key', key, obj1[key], obj2[key], 
                'obj1', obj1, 'obj2', obj2)
            
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
        
        else:    
                
            result, msgs = assert_obj1_to_obj2(obj, item, **kwargs)
            if result:
                return True, (item)
            else:
                if not msgs in list_msgs:
                    list_msgs.append(msgs)
    if id_keys_to_check:
        return False, ('id_keys_to_check: %r, item: %r, not found in %r',
            id_keys_to_check,
            [ obj[key] for key in id_keys_to_check],
            [ [obj[key] for key in id_keys_to_check] for obj in item_list ])
    else:
        return False, ('obj not found in list', obj, list_msgs)

def find_all_obj_in_list(list1, list2, **kwargs):
    msgs = ['not run yet']
    for item in list1:
        result, msgs = find_obj_in_list(item, list2, **kwargs)
        if not result:
            return False, msgs
    return True, msgs


# Run tests without a database:
# Example:
# ./manage.py test reports.SDFSerializerTest.test2_clean_data_sdf \
#      --settings=lims.settings_testing_debug --verbosity=2 \
#      --testrunner=reports.tests.NoDbTestRunner
# NOTE: when using this testrunner, the class finder is different; note that the 
# path excludes the module, so 
# "reports.SDFSerializerTest.test2_clean_data_sdf"
# not 
# "reports.test.SDFSerializerTest.test2_clean_data_sdf"
class NoDbTestRunner(DiscoverRunner):
  """ A test runner to test without database creation """

  def setup_databases(self, **kwargs):
    """ Override the database creation defined in parent class """
    pass

  def teardown_databases(self, old_config, **kwargs):
    """ Override the database teardown defined in parent class """
    pass


# NOTE:
# Create a TestRunner that will swallow exceptions on teardown of the test
# database.
# - motivation: use this runner for Travis tests so that the teardown error is
# not reported as a testing failure.
# 
# FIXME: override of testrunner class to catch teardown error:
# - caused by holding another db connection for the sqlalchemy bridge
# NOTE-(verify in DiscoverRunner): 
# when using this testrunner, the class finder is different; note that the 
# path excludes the module, so 
# "reports.SDFSerializerTest.test2_clean_data_sdf"
# not 
# "reports.test.SDFSerializerTest.test2_clean_data_sdf"
class IccblTestRunner(DiscoverRunner):
    
    @classmethod
    def add_arguments(cls, parser):
        super(IccblTestRunner, cls).add_arguments(parser)
        
        parser.add_argument(
            '--reinit_metahash',
            action='store_true',
            help=(
                'if set (with --keepdb), the _bootstrap_init_files method '
                'will be called (without --keepdb, _bootstrap_init_files'
                ' is always called)'))
        parser.add_argument(
            '--reinit_pattern [pattern]',
            action='store',
            help=(
                'if set (with --reinit_metahash), will only process input'
                'file actions with "resource" or "file" that matches the '
                'pattern given'))
        
    def teardown_databases(self, old_config, **kwargs):
        try:
            DiscoverRunner.teardown_databases(self, old_config, **kwargs)
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
        logger.info('>> test0_equivocal <<')
        
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
            [ True, 'True'],
            ['TRUE',1], 
            ['false','False'],
            ['FALSE','false'],
            ['False',False],
            # test that bytestrings encoded as UTF-8 are recognized
            # [u'\u03bc','\xce\xbc'] 
            ]
            ## TODO: test date strings 
        
        for [a,b] in test_true:
            result, msgs = equivocal(a,b)
            logger.info(
                'equivocal true %r to %r: %r: %r', 
                a, b, result, msgs)
            self.assertTrue(result, msgs)

        test_false = [
            ['1',''],
            ['1',2],    
            ['','1'],
            ['sal','val'],
            ['True','2'],
            ['True','0'],
            ['true',''],
            ['TRUE',False],
            ['false','true'],
            ['False',True],
            ]
        for [a,b] in test_false:
            result, msgs = equivocal(a,b)
            logger.info(
                'equivocal false: %r to %r: %r: %r', 
                a, b, result, msgs)
            self.assertFalse(result, '%r != %r, %r' % (a, b, msgs))
    
    def test1_assert_obj1_to_obj2(self):
        logger.info('test1_assert_obj1_to_obj2...')
        # test that all of obj1 in obj2 (may be more in obj2)
        obj1 = { 'one': '1', 'two': 'two', 'three':''}
        obj2 = { 'one': '1', 'two': 'two', 'three':'', 'four': 'blah'}
        result, msgs = assert_obj1_to_obj2(obj1, obj2)
        self.assertTrue(result, msgs)
        result, msgs = assert_obj1_to_obj2(obj2, obj1)
        self.assertFalse(result, msgs)
        obj3 = { 'one':'1', 'two':'2', 'three':'' }
        result, msgs = assert_obj1_to_obj2(obj3, obj1)
        self.assertFalse(result, (result,msgs))

        obj4 = { 'one':'1', 'two':'2', 'three':'' }
        obj5 = { 'one':'1', 'two':'2', 'three':None }
        result, msgs = assert_obj1_to_obj2(obj5, obj4)
        self.assertTrue(result, (result,msgs))
        logger.debug('====== test1_assert_obj1_to_obj2 done ====')

    def test2_find_obj_in_list(self):
        logger.info('====== test2_find_obj_in_list ====')
        
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
        logger.debug('%r, %r', result,msgs)
        self.assertFalse(result, msgs)

        obj_list = [ obj1a, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertFalse(result, msgs)
        logger.debug('====== test2_find_obj_in_list done ====')

    def test3_find_all_obj_in_list(self):
        logger.info('====== test3_find_all_obj_in_list ====')
        
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
        logger.debug('====== test3_find_all_obj_in_list done ====')

    def test4_user_example(self):
        logger.info('====== test4_user_example ====')
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
        
        self.assertTrue(obj['email'] == 'joe.tester@limstest.com')
        
        logger.debug('====== test4_user_example done ====') 
    
    
class SerializerTest(TestCase):

    def test_parse_val(self):
        
        test_array = [
            (('1', 'integer_test', 'integer'),1),
            (('true', 'boolean_test', 'boolean'),True),
            (('TRUE', 'boolean_test', 'boolean'),True),
            (('FALSE', 'boolean_test', 'boolean'),False),
            (('0', 'boolean_test', 'boolean'),False),
            (('1', 'boolean_test', 'boolean'),True),
            (('10', 'boolean_test', 'boolean'),False),
            (('x', 'boolean_test', 'boolean'),False),
            ((None, 'boolean_test', 'boolean'),None),
            (('2017-04-25', 'date_test', 'date'),
                dateutil.parser.parse('2017-04-25').date()),
        ]
        for test_data in test_array:
            result = parse_val(*test_data[0])
            self.assertEqual(result, test_data[1],
                '%r != %r' % (result,test_data))
        
        
    def test_csv(self):
        # NOTE this tests the obsoleted (non-streaming) Serializers;
        # new serializer uses csv.writer
        logger.debug('======== test_csv =========')
        directory = APP_ROOT_DIR
        serializer = CSVSerializer() 
        
        input = [['one','two', 'three', 'four', 'five','six','seven'    ],
                ['uno', '2', 'true', 'false', '',['a','b','c'], u'\u03bc']]
        
        input_data = csvutils.from_csv_iterate(input)
        for obj in input_data:
            self.assertTrue(obj['one']=='uno')
            self.assertTrue(obj['two']=='2')
            self.assertTrue(obj['three']=='true')
            self.assertTrue(obj['four']=='false')
            self.assertTrue(obj['five']=='')
            self.assertTrue(obj['six']==['a','b','c'])
            self.assertTrue(obj['seven']==u'\u03bc')
        csv_data = serializer.to_csv(input_data, root=None)

        with open(directory + '/reports/test_csv_.csv', 'w') as f:
            f.write(csv_data)
            
        with open(directory + '/reports/test_csv_.csv') as fin:    
            final_data = serializer.from_csv(fin.read(), root=None)
            for obj in final_data:
                self.assertTrue(obj['one']=='uno')
                self.assertTrue(obj['two']=='2')
                self.assertTrue(obj['three']=='true')
                self.assertTrue(obj['four']=='false')
                self.assertTrue(obj['five']=='')
                self.assertTrue(obj['six']==['a','b','c'])
                self.assertTrue(obj['seven']==u'\u03bc')
        
        # TODO: delete the file

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
            
            self.assertTrue(
                len(final_data)==len(records), 
                str(('len is', len(final_data),len(records))))
            for obj in final_data:
                for record in records:
                    if record['smiles'] == obj['smiles']:
                        for k,v in record.items():
                            self.assertTrue(
                                obj[k] == v.strip(), 
                                ('values not equal', k, record[k], 'read:', v))

    def test2_clean_data_sdf(self):
        ''' more extensive, read/write test '''

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


class XlsSerializerTest(SimpleTestCase):
    
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
            
            self.assertTrue(
                len(final_data)==len(records), 
                str(('len is', len(final_data),len(records))))
            for obj in final_data:
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
        
        self.assertTrue(set(diff_dict.keys()) 
            == set(['one','two','three','four','six']), diff_dict.keys())

        self.assertTrue('three' in diff_dict, diff_dict)
        self.assertTrue(diff_dict['two']==['value2a', 'value2b'])
        

class IResourceTestCase(SimpleTestCase):
   
    username = 'testsuper'
    password = 'pass'
    general_user_password = 'testpass1'
    
    
    """
    Override the Django SimpleTestCase, not using the TransactionTestCase
    necessary so that the SqlAlchemy connection can see the same database as the 
    django test code (Django ORM).

    Serves as a base class for tests, and, through _bootstrap_init_files, 
    as the setup for this module.

    FIXME: initialize the SqlAlchemyResource inside the transactions so that 
    the normal TransactionTestCase can be used.
    """
    def __init__(self,*args,**kwargs):
    
        super(IResourceTestCase, self).__init__(*args,**kwargs)
        
        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        self.csv_serializer=CSVSerializer() 
        self.serializer = LimsSerializer()
        self.api_client = TestApiClient(serializer=self.serializer)       

        # TODO: replace all tastypie.test.TestApiClient clients with the
        # django.test.client.Client instances:
        # Tastypie mucks with the HEADERS and uses the non-standard "format" arg:
        # resp = self.sr_api_client.get(
        #     resource_uri, authentication=self.get_credentials(), 
        #     format='xlsx', **data_for_get)
        # Tastypie PUT/POST requires that the data be serialized before posting,
        # and does not create the multipart/form-data header
        self.django_client = Client()
        settings.BACKGROUND_PROCESSING = False
            
    def setUp(self):
        super(IResourceTestCase, self).setUp()

    def _bootstrap_init_files(self, reinit_pattern=None):
        
        logger.info('_bootstrap_init_files...: %r', reinit_pattern)
        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        serializer=CSVSerializer() 
        testApiClient = TestApiClient(serializer=serializer) 
        input_actions_file = os.path.join(
            self.directory, 'api_init_actions.csv')
        
        self._run_api_init_actions(
            input_actions_file, reinit_pattern=reinit_pattern)

    def create_basic(self, username, password):
        """
        Creates & returns the HTTP ``Authorization`` header for use with BASIC
        Auth.
        """
        import base64
        return 'Basic %s' % base64.b64encode(
            ':'.join([username, password]).encode('utf-8')).decode('utf-8')

    def set_user_password(self, username, password):
        # assign password to the test user
        # NOTE: may only be done through the Django Model, for now
        # TODO: superuser should be able to assign password through secure connection
        userObj = User.objects.get(username=username)
        userObj.set_password(password)
        userObj.save()
        return userObj

    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)
        
    def _create_resource(
            self,input_data,resource_uri,resource_test_uri, 
            data_for_get= None, expect_fail=False, excludes=[]):
        
        _data_for_get = { 
            'limit': 0,
            'includes': '*',
            DJANGO_ACCEPT_PARAM: 'application/json'
        }
        if data_for_get:
            _data_for_get.update(data_for_get)
            
        logger.info('resource: %r, resource_test_uri: %r', 
            resource_uri, resource_test_uri)
        resp = self.api_client.post(
            resource_uri, format='json', data=input_data, 
            authentication=self.get_credentials(), **_data_for_get)
        # FIXME: tp uses status code 200 instead of 201
        if expect_fail:
            self.assertFalse(
                resp.status_code in [200,201], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)

            resp1 = self.api_client.get(
                resource_test_uri, format='json', 
                authentication=self.get_credentials(), data=_data_for_get)
            new_obj1 = self.deserialize(resp1)
            self.assertTrue(
                API_RESULT_DATA not in new_obj1 
                    or len(new_obj1[API_RESULT_DATA])==0, 
                'Error: failed create returns objects: %r' % new_obj1)
            return (new_obj,resp)
        else:
            self.assertTrue(
                resp.status_code in [200,201], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            logger.debug('post response: %r', new_obj)
            
        new_obj = self.get_single_resource(
            resource_test_uri, data_for_get=_data_for_get)
        result,msg = assert_obj1_to_obj2(input_data,new_obj, excludes=excludes)
        self.assertTrue(result, msg)
        logger.debug('item created: %r', new_obj)
        return new_obj
    
    def get_list_resource(self, resource_uri, data_for_get=None):
        _data_for_get = { 
            'limit': 0,
            'includes': '*'
        }
        if data_for_get:
            _data_for_get.update(data_for_get)
        logger.info('get from %r... %r', resource_uri, _data_for_get)
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), data=_data_for_get)
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        if API_RESULT_DATA in new_obj:
            new_obj = new_obj[API_RESULT_DATA]
        return new_obj
    
    def get_single_resource(self, resource_uri, data_for_get=None):
        '''
        Retrieve a single item from the resource_uri
        -- assertion failure if unsuccessful
        '''
        _data_for_get = { 
            'limit': 0,
            'includes': '*'
        }
        if data_for_get:
            _data_for_get.update(data_for_get)
        logger.info('get from %r... %r', resource_uri, _data_for_get)
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), data=_data_for_get)
        if resp.status_code == 404:
            logger.info('resp: %r', resp)
            logger.info('no resource found: %r, %r', resource_uri, data_for_get)
            return None
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        # NOTE: not all responses have data nested in API_RESULT_DATA
        if API_RESULT_DATA in new_obj:
            self.assertEqual(len(new_obj[API_RESULT_DATA]),1,
                'more than one object returned for: %r, returns: %r'
                % (resource_uri,new_obj))
            new_obj = new_obj[API_RESULT_DATA][0]
        logger.debug('obj: %r', new_obj)
        return new_obj
    
    def get_resource_from_server(self, resource_name):
        '''
        Utility to get a resource description (schema) from the server
        '''
        resource_uri = BASE_URI + '/resource/' + resource_name
        logger.info('Get the resource schema: %r', resource_uri )
        return self.get_single_resource(resource_uri)
        
    def _patch_test(
        self,resource_name, filename, keys_not_to_check=['resource_uri'], 
        id_keys_to_check=[], data_for_get={}):
        '''
        @param data_for_get - dict of extra header information to send with the 
        GET request
        '''
        logger.info('_patch test: %r, %r', resource_name, filename)
        data_for_get.setdefault('limit', 999 )
        data_for_get.setdefault('includes', '*' )
        data_for_get[HEADER_APILOG_COMMENT] =  'patch_test: %s' % filename
        resource_uri = BASE_URI + '/' + resource_name
        
        logger.debug('===resource_uri: %r', resource_uri)
        with open(filename) as bootstrap_file:
            # Note: we have to deserialize the input, because the TP test method 
            # will expect a python data object, which it will serialize!
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())
            
            logger.info('Submitting patch... %r: %r', filename, resource_uri)
            resp = self.api_client.patch(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(
                resp.status_code <= 204, 
                (resp.status_code, self.get_content(resp)))
            logger.info('get: %r,%r, %r', resource_uri, data_for_get, id_keys_to_check)
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data=data_for_get)
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            
            new_obj = self.deserialize(resp)
            for inputobj in input_data[API_RESULT_DATA]:
                # use id_keys_to_check to perform a search only on those keys
                result, outputobj = find_obj_in_list(
                    inputobj,new_obj[API_RESULT_DATA], 
                    id_keys_to_check=id_keys_to_check, 
                    excludes=keys_not_to_check )
                logger.debug('objects returned: %r', new_obj[API_RESULT_DATA])
                self.assertTrue(
                    result, 
                    'not found: %r, msg: %r' % (inputobj, outputobj))
                # once found, perform equality based on all keys (in obj1)
                logger.debug('found: %r: %r', inputobj.get('scope'), inputobj.get('key'))

                # For fields patches only:
                # NOTE: skip_null_values added to address fields with parent "ref" set
                skip_null_values = False
                if 'scope' in inputobj and 'fields.' in inputobj.get('scope'):
                    ref = inputobj.get('ref')
                    if ref is not None:
                        logger.info('%r: ref: %r',inputobj.get('key'), ref)
                        skip_null_values = True
                result, msg = assert_obj1_to_obj2(inputobj, outputobj,
                    excludes=keys_not_to_check, skip_null_values=skip_null_values)
                self.assertTrue(result,
                    'not equal: %r: %r - %r' % ( msg, inputobj, outputobj))
            
            # return both collections for further testing
            return (input_data[API_RESULT_DATA], new_obj[API_RESULT_DATA]) 

    def _put_test(
        self, resource_name, filename, keys_not_to_check=['resource_uri'], 
        id_keys_to_check=[], data_for_get={}):
        '''
        id_keys_to_check if the resource data has been loaded, 
            then these are id keys to check to see if they are being used in 
            the returned resource_uri field
        '''
        data_for_get.setdefault('limit', 0 )
        data_for_get.setdefault('includes', '*' )
        data_for_get.setdefault( 
            HEADER_APILOG_COMMENT, 'put_test: %s' % filename )
        resource_uri = BASE_URI + '/' + resource_name

        with open(filename) as bootstrap_file:
            # NOTE: have to deserialize the input, because the TP test method 
            # will expect a python data object, which it will serialize!
            input_data = self.csv_serializer.from_csv(bootstrap_file.read())
            
            resp = self.api_client.put(
                resource_uri, format='csv', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(
                resp.status_code <= 204, 
                '%r: %r' % (resp.status_code, self.get_content(resp)))
            logger.info('get: %r', resource_uri)
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get)
            self.assertTrue(
                resp.status_code in [200], 
                '%r: %r' % (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            # do a length check, since put will delete existing resources
            self.assertEqual(
                len(new_obj[API_RESULT_DATA]), len(input_data[API_RESULT_DATA]), 
                str(('input length != output length: ',
                    len(new_obj[API_RESULT_DATA]), 
                    len(input_data[API_RESULT_DATA]),
                    input_data,'\n\n', new_obj)))
            
            for inputobj in input_data[API_RESULT_DATA]:
                result, outputobj = find_obj_in_list(
                    inputobj,new_obj[API_RESULT_DATA],
                    id_keys_to_check=id_keys_to_check,
                    excludes=keys_not_to_check)
                self.assertTrue(
                    result, 
                    'not found: %r: %r' % (outputobj, new_obj[API_RESULT_DATA]))
                # once found, perform equality based on all keys (in obj1)
                result, msg = assert_obj1_to_obj2(inputobj, outputobj,
                    excludes=keys_not_to_check)
                self.assertTrue(
                    result,
                    'not equal: %r: %r - %r' % ( msg, inputobj, outputobj))
            #TODO: GET the apilogs expected and test them
            return (input_data[API_RESULT_DATA], new_obj[API_RESULT_DATA]) 

    def _run_api_init_actions(self, input_actions_file, reinit_pattern=None):
        logger.info( '_run_api_init_actions: %r, reinit_pattern: %r', 
            input_actions_file, reinit_pattern)
        
        if not os.path.exists(input_actions_file):
            raise AssertionError('no such file %s' % input_actions_file)
        with open(input_actions_file) as input_file:
            api_init_actions = csvutils.from_csv(input_file)
            if not api_init_actions:
                raise AssertionError(
                    'no actions read from file: %s' % input_actions_file)
            for action in api_init_actions:
                command = action['command'].lower() 
                resource = action['resource'].lower()
                filename = action.get('file',None)
                
                def processAction():
                    logger.warn('process action: %r', action)
                    if command == 'delete':
                        # no need to delete with the empty test database
                        pass
                    else:
                        if filename is None:
                            raise AssertionError(
                                '"put" and "patch" require a filename argument: %r'
                                    % action)
                        data_file = os.path.join(self.directory,filename)
                        logger.info('Data File: %r', data_file)
                        if command == 'put':
                            self._put_test(
                                resource, data_file, 
                                keys_not_to_check=['resource_uri'])
                        elif command == 'patch':
                            id_keys_to_check = []
                            if resource == 'field':
                                id_keys_to_check = ['key','scope']
                            self._patch_test(
                                resource, data_file, 
                                id_keys_to_check=id_keys_to_check,
                                keys_not_to_check=['resource_uri'])
                        else:
                            raise AssertionError(
                                'Unknown API command %r' % action)
                        logger.info('Completed, data File: %r', data_file)
           
                if reinit_pattern is not None:
                    pattern = re.compile(
                        r'%s'%reinit_pattern,flags=re.IGNORECASE)
                    if pattern.search(resource):
                        logger.info(
                            'reinit_pattern resource match: %r', resource)
                        processAction()
                    else:
                        logger.info(
                            'reinit_pattern: %r doesnt match resource: %r',
                            reinit_pattern, resource)
                    if pattern.search(filename):
                        logger.info(
                            'reinit_pattern filename match: %r', filename)
                        processAction()
                    else:
                        logger.info(
                            'reinit_pattern: %r doesnt match filename: %r',
                            reinit_pattern, filename)
                else:
                    processAction()
                    
        logger.info('Completed, input actions file: %r', input_actions_file)
        
    def get_content(self, resp):
        
        return self.serializer.get_content(resp);
    
    def deserialize(self, resp):

        return self.serializer.deserialize(
            self.get_content(resp), resp['Content-Type'])

    def serialize(self, data, format):
        content_type = self.serializer.get_content_type_for_format(format)
        return self.serializer.serialize(data, content_type)

runTestApiInit = [False]
def setUpModule():
    logger.info('=== setup Module')

    # FIXME: running the bootstrap method as a test suite setup:
    # TODO: create a local TestRunner,TestSuite,
    # so that this can be run before the suite
    global runTestApiInit
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
            if 'TestApiInit' in arg:
                runTestApiInit[0] = True
        if keepdb and runTestApiInit[0]:
            raise Exception(
                'The TestApiInit test cannot be run with an existing test database')
    
    logger.info(
        'init vars: keepdb: %r, reinit_metahash: %r,reinit_pattern: %r, '
        'TestApiInit: %r ', 
        keepdb,reinit_metahash,reinit_pattern,runTestApiInit[0])

    # Set up a superuser
    print 'create a superuser...'
    try:
        logger.info('create/find superuser %s...', IResourceTestCase.username)
        IResourceTestCase.user = User.objects.get(
            username=IResourceTestCase.username)
        logger.warn('superuser found: %r', IResourceTestCase.user)
        logger.warn('users: %r', [str(u) for u in User.objects.all()])
    except ObjectDoesNotExist:
        logger.warn('creating superuser: %s', IResourceTestCase.username)
        IResourceTestCase.user = User.objects.create_superuser(
            IResourceTestCase.username, '1testsuperuser@example.com', 
            IResourceTestCase.password)
    print 'superuser created.'
    logger.info('reinit_pattern: %r', reinit_pattern)
    if (reinit_metahash or not keepdb) and not runTestApiInit[0]:
        testContext = IResourceTestCase(methodName='_bootstrap_init_files')
        testContext.setUp()
        testContext._bootstrap_init_files(reinit_pattern=reinit_pattern)
        logger.info('database initialization finished')
    else:
        print 'skip database metahash initialization when using keepdb'


    logger.info('=== setup Module done')

def tearDownModule():
    logger.info('=== teardown Module')


class TestApiClient(object):


    def __init__(self, serializer=None):
        """
        """
        self.client = Client()
        self.serializer = serializer

        if not self.serializer:
            self.serializer = LimsSerializer()
        
        super(TestApiClient, self).__init__()
        
        
    def get(self, uri, format='json', data=None, authentication=None, **kwargs):
        """
        Performs a simulated ``GET`` request to the provided URI.
        """
        content_type = self.serializer.get_content_type_for_format(format)
        if DJANGO_ACCEPT_PARAM not in kwargs:
            kwargs[DJANGO_ACCEPT_PARAM] = content_type

        # GET & DELETE are the only times we don't serialize the data.
        if data is not None:
            kwargs['data'] = data

        if authentication is not None:
            kwargs['HTTP_AUTHORIZATION'] = authentication

        return self.client.get(uri, **kwargs)

    def post(
        self, uri, format='json', data=None, authentication=None, **kwargs):
        """
        Performs a simulated ``POST`` request to the provided URI.

        Optionally accepts a ``data`` kwarg. 
        **Unlike** ``GET``, in ``POST`` the ``data`` gets serialized & sent 
        as the body instead of becoming part of the URI.
        """
        content_type = self.serializer.get_content_type_for_format(format)
        logger.info('content_type: %r', content_type)
        if 'content_type' not in kwargs:
            kwargs['content_type'] = content_type
        if DJANGO_ACCEPT_PARAM not in kwargs:
            kwargs[DJANGO_ACCEPT_PARAM] = content_type

        if data is not None:
            kwargs['data'] = self.serializer.serialize(data, content_type)

        if authentication is not None:
            kwargs['HTTP_AUTHORIZATION'] = authentication

        return self.client.post(uri, **kwargs)

    def put(self, uri, format='json', data=None, authentication=None, **kwargs):
        """
        Performs a simulated ``PUT`` request to the provided URI.

        Optionally accepts a ``data`` kwarg. 
        **Unlike** ``GET``, in ``PUT`` the ``data`` gets serialized & sent as
        the body instead of becoming part of the URI.
        """
        content_type = self.serializer.get_content_type_for_format(format)
        logger.info('content_type: %r', content_type)
        if 'content_type' not in kwargs:
            kwargs['content_type'] = content_type
        if DJANGO_ACCEPT_PARAM not in kwargs:
            kwargs[DJANGO_ACCEPT_PARAM] = content_type

        if data is not None:
            kwargs['data'] = self.serializer.serialize(data, content_type)

        if authentication is not None:
            kwargs['HTTP_AUTHORIZATION'] = authentication

        return self.client.put(uri, **kwargs)

    def patch(
        self, uri, format='json', data=None, authentication=None, **kwargs):
        """
        Performs a simulated ``PATCH`` request to the provided URI.

        Optionally accepts a ``data`` kwarg. 
        **Unlike** ``GET``, in ``PATCH`` the ``data`` gets serialized & sent 
        as the body instead of becoming part of the URI.
        """
        content_type = self.serializer.get_content_type_for_format(format)
        if 'content_type' not in kwargs:
            kwargs['content_type'] = content_type
        if DJANGO_ACCEPT_PARAM not in kwargs:
            kwargs[DJANGO_ACCEPT_PARAM] = content_type

        if data is not None:
            kwargs['data'] = self.serializer.serialize(data, content_type)

        if authentication is not None:
            kwargs['HTTP_AUTHORIZATION'] = authentication

        # Django doesn't support PATCH natively.
        parsed = urlparse.urlparse(uri)
        r = {
            'CONTENT_TYPE': content_type,
            'PATH_INFO': self.client._get_path(parsed),
            'QUERY_STRING': parsed[4],
            'REQUEST_METHOD': 'PATCH',
        }
        if data:
            r['CONTENT_LENGTH'] = len(kwargs['data'])
            r['wsgi.input'] = FakePayload(kwargs['data'])
        else:
            r['CONTENT_LENGTH'] = 0
        
        r.update(kwargs)
        return self.client.request(**r)

    def delete(
        self, uri, format='json', data=None, authentication=None, **kwargs):
        """
        Performs a simulated ``DELETE`` request to the provided URI.

        """
        content_type = self.serializer.get_content_type_for_format(format)
        logger.info('content_type: %r', content_type)
        if 'content_type' not in kwargs:
            kwargs['content_type'] = content_type
        if DJANGO_ACCEPT_PARAM not in kwargs:
            kwargs[DJANGO_ACCEPT_PARAM] = content_type

        # GET & DELETE are the only times we don't serialize the data.
        if data is not None:
            kwargs['data'] = data

        if authentication is not None:
            kwargs['HTTP_AUTHORIZATION'] = authentication

        return self.client.delete(uri, **kwargs)
    

@unittest.skipIf(runTestApiInit[0]!=True, 
                 "Only run this test if testing the api initialization")   
class TestApiInit(IResourceTestCase):
    
    def setUp(self):
        
        super(TestApiInit, self).setUp()

        self.directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        self.csv_serializer=CSVSerializer() 
        
    def tearDown(self):
        super(TestApiInit, self).tearDown()
        
        logger.info('delete resources')
        MetaHash.objects.all().delete()

    def test0_bootstrap_metahash(self):
        
        logger.info('run test0_bootstrap_metahash....')
        
        resource_uri = BASE_URI + '/field'
        # NOTE:
        # In order for the metahash resource to work, the metahash itself must 
        # be "bootstrapped":
        # The metahash must be filled with the fields that describe itself
        # these are the "bootstrap" fields.
        bootstrap_items =  [   
            {
                'key': 'scope',
                'scope': 'fields.field',
                'ordinal': 0,
            },
            {
                'key': 'key',
                'scope': 'fields.field',
                'ordinal': 1,   
            },
            {
                'key': 'ordinal',
                'scope': 'fields.field',
                'ordinal': 2,   
            },
            {
                'key': 'json_field_type',
                'scope': 'fields.field',
                'ordinal': 3,    
            },
            {
                'key': 'data_type',
                'scope': 'fields.field',
                'ordinal': 4,
                'json_field_type': 'fields.CharField',    
            },
            {
                'key': 'visibility',
                'scope': 'fields.field',
                'ordinal': 5,    
                'json_field_type': 'fields.ListField',    
            },
            {
                'key': 'editability',
                'scope': 'fields.field',
                'ordinal': 6,    
                'json_field_type': 'fields.ListField',    
            },
            {
                'key': 'table',
                'scope': 'fields.field',
                'ordinal': 7,    
                'json_field_type': 'fields.CharField',    
            },
            {
                'key': 'field',
                'scope': 'fields.field',
                'ordinal': 8,    
                'json_field_type': 'fields.CharField',    
            },
            {
                'key': 'regex',
                'scope': 'fields.field',
                'ordinal': 9,    
                'json_field_type': 'fields.CharField',    
            },
            {
                'key': 'vocabulary_scope_ref',
                'scope': 'fields.field',
                'ordinal': 10,    
                'json_field_type': 'fields.CharField',    
            },
            {
                'key': 'display_options',
                'scope': 'fields.field',
                'ordinal': 11,    
                'json_field_type': 'fields.CharField',    
            },
            {
                'key': 'value_template',
                'scope': 'fields.field',
                'ordinal': 12,    
                'json_field_type': 'fields.CharField',    
            },
        ]    
        
        # Note: the initial bootstrap fields _must_ be patched as a list:
        # The server requires these fields to bootstrap          
        resp = self.api_client.patch(
            resource_uri, format='json', 
            data={ API_RESULT_DATA: bootstrap_items }, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
            
        logger.debug('created items, now get them')
        data_for_get = {
             'limit': '*',
             'scope': 'fields.field'
        }
        resp = self.api_client.get(
            resource_uri, format='json', 
            authentication=self.get_credentials(), )
        new_obj = self.deserialize(resp)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        new_objects = new_obj[API_RESULT_DATA]
        self.assertEqual(
            len(new_objects), len(bootstrap_items), 
            (len(new_objects), 'expected', len(bootstrap_items), new_objects))
        for inputobj in bootstrap_items:
            result, outputobj = find_obj_in_list(inputobj,new_objects)
            self.assertTrue(result, ('not found', inputobj, outputobj ) )
    
    def test2_create_resource_fields(self):
        
        self.test0_bootstrap_metahash()

        logger.info('run test2_create_resource_fields....')
        
        resource_uri = BASE_URI + '/field'
        bootstrap_items =  [   
            {
                'key': 'scope',
                'scope': 'fields.resource',
                'ordinal': 0,
                'data_type': 'string',
            },
            {
                'key': 'key',
                'scope': 'fields.resource',
                'ordinal': 1,   
                'data_type': 'string',
            },
            {
                'key': 'ordinal',
                'scope': 'fields.resource',
                'ordinal': 2,   
                'data_type': 'integer',
            },
            {
                'key': 'title',
                'scope': 'fields.resource',
                'ordinal': 3,    
                'json_field_type': 'fields.CharField',    
                'data_type': 'string',
            },
            {
                'key': 'id_attribute',
                'scope': 'fields.resource',
                'ordinal': 4,    
                'json_field_type': 'fields.ListField',    
                'data_type': 'list',
            },
            {
                'key': 'api_name',
                'scope': 'fields.resource',
                'ordinal': 5,    
                'json_field_type': 'fields.CharField',    
                'data_type': 'string',
            },
            {
                'key': 'table',
                'scope': 'fields.resource',
                'ordinal': 6,    
                'json_field_type': 'fields.CharField',    
                'data_type': 'string',
            },
            {
                'key': 'supertype',
                'scope': 'fields.resource',
                'ordinal': 7,    
                'json_field_type': 'fields.CharField',    
                'data_type': 'string',
            },
        ]    
        logger.info('patch %r....',resource_uri)
        resp = self.api_client.patch(
            resource_uri, format='json', 
            data={ API_RESULT_DATA: bootstrap_items }, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
            
        logger.debug('created items, now get them')
        data_for_get = {
             'limit': '*',
             'scope': 'fields.resource'
        }
        resp = self.api_client.get(
            resource_uri, format='json',
            data = data_for_get, 
            authentication=self.get_credentials(), )
        new_obj = self.deserialize(resp)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        new_objects = new_obj[API_RESULT_DATA]
        self.assertEqual(
            len(new_objects), len(bootstrap_items), (len(new_objects)))

        for inputobj in bootstrap_items:
            result, outputobj = find_obj_in_list(inputobj,new_objects)
            self.assertTrue(result, ('not found', inputobj, outputobj ) )
    
    def test3_create_resources(self):

        self.test2_create_resource_fields()

        logger.info('create the fields for each resource...')
        resource_uri = BASE_URI + '/field'
        bootstrap_items =  [   
            {
                'key': 'key1',
                'scope': 'fields.resource1',
                'ordinal': 0,
                'data_type': 'string',
                'visibility': ['l','d'],
            },
            {
                'key': 'key2',
                'scope': 'fields.resource1',
                'ordinal': 1,   
                'data_type': 'string',
            },
            {
                'key': 'keya',
                'scope': 'fields.resource2',
                'ordinal': 0,   
                'data_type': 'string',
                'visibility': ['l','d'],
                # create a dependency  on the keyb field
                'display_options': '{ hrefTemplate: "resource2/{keyb}" }'
            },
            {
                'key': 'keyb',
                'scope': 'fields.resource2',
                'ordinal': 1,    
                'data_type': 'string',
            },
        ]
        logger.info('patch %r....',resource_uri)
        resp = self.api_client.patch(
            resource_uri, format='json', 
            data={ API_RESULT_DATA: bootstrap_items }, 
            authentication=self.get_credentials())
            
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        
        logger.info('create the resource table entries...')
        resource_uri = BASE_URI + '/resource'
        bootstrap_items =  [   
            {
                'key': 'resource1',
                'scope': 'resource',
                'ordinal': 0,
                'title': "Resource 1",
                'id_attribute': ['key1','key2'],
                'api_name': 'reports'
            },
            {
                'key': 'resource2',
                'scope': 'resource',
                'ordinal': 0,
                'title': "Resource 2",
                'id_attribute': ['keya','keyb'],
                'api_name': 'reports',
                'supertype': 'resource1'
            },
        ]    
        logger.info('patch %r....',resource_uri)
        resp = self.api_client.patch(
            resource_uri, format='json', 
            data={ API_RESULT_DATA: bootstrap_items }, 
            authentication=self.get_credentials())
            
        self.assertTrue(resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
            
        logger.debug('created items, now get them')
        data_for_get = {
             'limit': 0,
        }
        resp = self.api_client.get(
            resource_uri, format='json',
            data = data_for_get, 
            authentication=self.get_credentials(), )
        new_obj = self.deserialize(resp)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

        new_objects = new_obj[API_RESULT_DATA]
        self.assertEqual(
            len(new_objects), len(bootstrap_items), (len(new_objects)))

        for inputobj in bootstrap_items:
            result, outputobj = find_obj_in_list(inputobj,new_objects)
            self.assertTrue(result, ('not found', inputobj, outputobj ) )
            if outputobj['key'] == 'resource1':
                resource1 = outputobj
            if outputobj['key'] == 'resource2':
                resource2 = outputobj
                
        logger.info('for resource2, verify supertype fields are inherited...')
        resource1_fields = resource1['fields'].values()
        resource2_fields = resource2['fields'].values()
        logger.info('fields: %r', resource2_fields)
        self.assertTrue(
            len(resource2_fields)==4, 
            'Should have %d fields but found %d' %(4,len(resource2_fields)))
        for resource1_field in resource1_fields:
            resource1_field['scope'] = 'fields.resource2'
            result, found_field = find_obj_in_list(
                resource1_field, resource2_fields)
            self.assertTrue(
                result, 'cannot find supertype field: %r' % found_field )
            logger.info('found field: %r', found_field)

    def test99_api_init(self):
        
        logger.debug('================ reports test2_api_init =============== ')
        serializer=CSVSerializer() 
        # Note: doesn't work for post, see TestApiClient.post() method, it is 
        # incorrectly "serializing" the data before posting
        testApiClient = TestApiClient(serializer=serializer) 
        
        filename = os.path.join(self.directory,'api_init_actions.csv')
        with open(filename) as input_file:
            # NOTE: we have to deserialize the input, because the TP test method 
            # will expect a python data object, which it will serialize!
            api_init_actions = serializer.from_csv(input_file.read(), root=None)
            
            bootstrap_files = [
                'metahash_fields_initial.csv',
                'metahash_fields_initial_patch.csv',
                'metahash_fields_resource.csv',
                'metahash_resource_data.csv',
                'metahash_fields_vocabulary.csv',
                'vocabulary_data.csv']
            for action in api_init_actions:
                
                logger.debug('\n+++===== processing action', json.dumps(action))
                command = action['command'].lower() 
                resource = action['resource'].lower()
                resource_uri = BASE_URI + '/' + resource
                
                if command == 'delete':
                    resp = testApiClient.delete(
                        resource_uri, authentication=self.get_credentials())
                    self.assertTrue(
                        resp.status_code in [202, 204], 
                        (resp.status_code,self.get_content(resp)))
                
                else:
                    filename = os.path.join(self.directory,action['file'])
                    search_excludes = []
                    # exclude 'resource_uri' from equivalency check during 
                    # bootstrap, because resource table needs to be loaded for
                    # the uri generation
                    if action['file'] in bootstrap_files:
                        search_excludes = ['resource_uri'] 
                    logger.debug('+++++++++++processing file %r', filename)
                    with open(filename) as data_file:
                        # NOTE: Must deserialize the input, because the TP test
                        # method will expect a python data object, 
                        # which it will serialize!
                        input_data = serializer.from_csv(data_file.read())
                        
                        if command == 'put':
                            resp = testApiClient.put(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
                            self.assertTrue(
                                resp.status_code <= 204, 
                                (resp.status_code, self.get_content(resp)))
                            
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), 
                                data={ 'limit': 0 })
                            self.assertTrue(
                                resp.status_code in [200], 
                                (resp.status_code,self.get_content(resp)))
                            new_obj = self.deserialize(resp)
                            result, msgs = find_all_obj_in_list(
                                input_data[API_RESULT_DATA], 
                                new_obj[API_RESULT_DATA], 
                                excludes=search_excludes)
                            self.assertTrue(
                                result, (command, 'input file', filename, 
                                             msgs, new_obj[API_RESULT_DATA]) )
                        
                        elif command == 'patch':
                            resp = testApiClient.patch(
                                resource_uri, format='csv', data=input_data, 
                                authentication=self.get_credentials() )
                            self.assertTrue(
                                resp.status_code in [200], 
                                (resp.status_code, self.get_content(resp)))
                            resp = testApiClient.get(
                                resource_uri, format='json', 
                                authentication=self.get_credentials(), 
                                data={ 'limit': 0 } )
                            self.assertTrue(
                                resp.status_code in [200], 
                                (resp.status_code, self.get_content(resp)))
                            new_obj = self.deserialize(resp)
                            with open(filename) as f2:
                                input_data2 = serializer.from_csv(f2.read())
                                result, msgs = find_all_obj_in_list(
                                    input_data2[API_RESULT_DATA], 
                                    new_obj[API_RESULT_DATA], 
                                    excludes=search_excludes)
                                self.assertTrue(
                                    result, 
                                    ( command, 'input file', filename, msgs ) )
                        
                        elif command == 'post':
                            self.fail((
                                'resource entry: ' + json.dumps(action) + '; '
                                'cannot POST multiple objects to tastypie; '
                                'therefore the "post" command is invalid with '
                                'the initialization scripts'))
                        else:
                            self.fail('unknown command: %r, %r',command,
                                      json.dumps(action))
        
class VocabularyResource(IResourceTestCase):
    
    def test1_create_read(self):
        
        test_vocabs = [
            {'scope': 'test.vocab', 'key': 'test1', 'ordinal': 1, 
             'description': 'test1 vocab', 'title': 'Test 1' },
            {'scope': 'test.vocab', 'key': 'test2', 'ordinal': 2, 
             'description': 'test2 vocab', 'title': 'Test 2' },
            {'scope': 'test.vocab', 'key': 'test3', 'ordinal': 3, 
             'description': 'test3 vocab', 'title': 'Test 3' },
        ]
        try:       
            uri = BASE_URI + '/vocabulary'
            resp = self.api_client.patch(uri, 
                format='json', data={ API_RESULT_DATA: test_vocabs }, 
                authentication=self.get_credentials())
            self.assertTrue(
                resp.status_code <= 204, 
                (resp.status_code, self.get_content(resp)))
        except Exception, e:
            logger.exception('on creating: %r', test_vocabs)
            raise

        logger.debug('created vocab items, now get them')
        resp = self.api_client.get(uri, format='json', 
            authentication=self.get_credentials(), 
            data={ 
                'limit': 0,
                'scope__eq': 'test.vocab' })
        new_obj = self.deserialize(resp)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        self.assertEqual(
            len(new_obj[API_RESULT_DATA]), len(test_vocabs), (new_obj))
        
        for i,item in enumerate(test_vocabs):
            result, obj = find_obj_in_list(item, new_obj[API_RESULT_DATA])
            self.assertTrue(
                result, 
                ('vocab item not found', item, new_obj[API_RESULT_DATA]))


class UserUsergroupSharedTest(object):
            
    def test1_create_user_with_permissions(self, test_log=False):
        logger.info('test1_create_user_with_permissions...')
        
        filename = os.path.join(self.directory,'test_data/users1.csv')
        from datetime import datetime
        data_for_get = { HEADER_APILOG_COMMENT: 
            'patch_test: file: %s, %s' % (filename, datetime.now().isoformat()),
            'includes': '*' }
        input_data,output_data = self._patch_test(
            'user', filename, id_keys_to_check=['username'],
            data_for_get=data_for_get)
        
        # Test the logs
        # FIXME: create separate tests for the apilogs
        if test_log:
            resource_uri = BASE_URI + '/apilog' #?ref_resource_name=record'
            logger.info('get: %r, %r', resource_uri, data_for_get)
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 'limit': 999, 'ref_resource_name': 'user',
                    'comment__eq': data_for_get[HEADER_APILOG_COMMENT]}) 
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            objects = new_obj[API_RESULT_DATA]
            # look for 6 "CREATE" logs
            self.assertEqual( 
                len(objects), 6, 
                str((6,len(objects), 'wrong # of api logs', objects)))
            for obj in objects:
                self.assertTrue(
                    obj['api_action'] == VOCAB.apilog.api_action.CREATE, 
                    ('action should be create', obj))
                self.assertTrue( 
                    obj['comment'] == data_for_get[HEADER_APILOG_COMMENT], 
                    ('comment not set', data_for_get[HEADER_APILOG_COMMENT],obj))
        

class UserResource(IResourceTestCase, UserUsergroupSharedTest):
    '''
    NOTE: User/Group permissions are set declaritively, granting blanket
    "read" or "write" permissions for a resource to administrative users.
    See the "DataSharingLevel" tests in the "db" module for LIMS (Screensaver
    Screening User) read permissions testing.
    '''

    def setUp(self):
        super(UserResource, self).setUp()

    def tearDown(self):
        IResourceTestCase.tearDown(self)
        
        logger.info('delete users, including: %r', self.username)
        UserGroup.objects.all().delete()
        UserProfile.objects.all().exclude(username=self.username).delete()
        User.objects.all().exclude(username=self.username).delete()
        ApiLog.objects.all().delete()
    
    def test0_create_user(self):
        
        logger.info('test0_create_user...')
        
        # create some simple users
        bootstrap_items = [   
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
        ]
        uri = BASE_URI + '/user'
        resp = self.api_client.patch(uri, 
            format='json', data={ API_RESULT_DATA: bootstrap_items}, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code <= 204, 
            (resp.status_code, self.get_content(resp)))

        logger.debug('created users, now GET them')
        data_for_get = { 'limit': 0 }
        data_for_get['username__in'] = [u['username'] for u in bootstrap_items]
        
        resp = self.api_client.get(uri, format='json', 
            authentication=self.get_credentials(), data=data_for_get)
        new_obj = self.deserialize(resp)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        self.assertEqual(
            len(new_obj[API_RESULT_DATA]), 3, (new_obj))
        
        for i,item in enumerate(bootstrap_items):
            result, obj = find_obj_in_list(item, new_obj[API_RESULT_DATA])
            self.assertTrue(
                result, 
                ('bootstrap item not found', item, new_obj[API_RESULT_DATA]))

    def test2_patch_user_permissions(self,test_log=False):
        '''
        @param test_log set to false so that only the first time tests the logs
        '''
        
        logger.info('test2_patch_user_permissions...')
        self.test1_create_user_with_permissions(test_log=test_log)
        
        logger.debug('==== test2_patch_user_permissions start =====')
        filename = os.path.join(self.directory,'test_data/users2_patch.csv')
        data_for_get = { HEADER_APILOG_COMMENT: 'patch_test: %s' % filename }

        input_data,output_data = self._patch_test(
            'user', filename, data_for_get=data_for_get)
        
        logger.info('patched output: %r', output_data)
        # # test the logs
        
    def test3_user_read_permissions(self):
        '''
        Try to do something we don't have permissions for - 
        read a resource (metahash) 
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        '''
        logger.info('test3_user_read_permissions...')
        self.test2_patch_user_permissions(test_log=False)
                
        # assign password to the test user
        username = 'sde1'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # Try to do some unauthorized actions

        test_resource_uri = '/'.join([BASE_URI, 'apilog'])
        resp = self.api_client.get(
            test_resource_uri, format='json', data={
                'ref_resource_name': 'user'}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(
            resp.status_code in [403], 
            (resp.status_code, self.get_content(resp)))
        
        # Now add the needed permission
        
        user_patch = {
            'resource_uri': 'user/' + username,
            'permissions': ['resource/apilog/read'] };

        uri = BASE_URI + '/user/' + username
        logger.debug('add permission to user: %r: %r', user_patch, uri)
        resp = self.api_client.patch( uri, 
                    format='json', data=user_patch, 
                    authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code <= 204,
            (resp.status_code, self.get_content(resp)))  
        
        # now try again as the updated user:
        
        resp = self.api_client.get(
            test_resource_uri, format='json', data={
                'ref_resource_name': 'user'}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
       
    def test4_user_write_permissions(self):
        '''
        Try to do something we don't have permissions for - 
        write a resource (user),
        then give that permission to the user, and try again 
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        '''
        logger.info('test4_user_write_permissions...')
        self.test2_patch_user_permissions(test_log=False)
                
        # assign password to the test user
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # Try to do some unauthorized actions

        resource_uri = BASE_URI + '/usergroup'
        
        logger.info('try unauthorized action: adding a usergroup...')
        
        json_data = { API_RESULT_DATA: [{ 'name': 'test_group_x' }] }
        
        resp = self.api_client.patch(
            resource_uri, format='json', data=json_data, 
            authentication=self.create_basic(username, password) )
        logger.info('resp: %r, %r', resp.status_code, self.get_content(resp))
        self.assertTrue(
            resp.status_code in [403], 
            (resp.status_code, self.get_content(resp)))
        
        # Now add the needed permission
        
        user_patch = {
            'resource_uri': 'user/' + username,
            'permissions': ['resource/usergroup/write'] };

        logger.debug(
            'now add the permission needed to this user: %r', user_patch)
        uri = BASE_URI + '/user/' + username
        resp = self.api_client.patch( uri, 
                    format='json', data=user_patch, 
                    authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code <= 204, 
            (resp.status_code, self.get_content(resp)))  
        
        logger.info(
            'action now authorized: adding a usergroup: %r', resource_uri)
        
        resp = self.api_client.patch(
            resource_uri, format='json', data=json_data, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(
            resp.status_code <= 204, 
            (resp.status_code, self.get_content(resp)))

        logger.debug('==== test4_user_write_permissions done =====')
        

class UserGroupResource(IResourceTestCase, UserUsergroupSharedTest):

    def setUp(self):
        super(UserGroupResource, self).setUp()

    def tearDown(self):
        IResourceTestCase.tearDown(self)
        
        logger.info('delete users, including: %r', self.username)
        UserGroup.objects.all().delete()
        UserProfile.objects.all().exclude(username=self.username).delete()
        User.objects.all().exclude(username=self.username).delete()
        ApiLog.objects.all().delete()

    def test2_create_usergroup_with_permissions(self):
        logger.info('test2_create_usergroup_with_permissions...')
        #create users
        self.test1_create_user_with_permissions(test_log=False)
        
        filename = os.path.join(self.directory,'test_data/usergroups1.csv')
        # Note: Excluding sub_groups here because the one sub_group is set when 
        # "testGroupX" sets super_groups=['testGroup3']; and thereby testGroup3 
        # gets sub_groups=['testGroupX']; although that's not in the input file.
        self._patch_test('usergroup', filename, id_keys_to_check=['name'],
            keys_not_to_check=['sub_groups'])

    def test3_patch_users_groups(self):
        logger.info('test3_patch_users_groups...')
        self.test2_create_usergroup_with_permissions()
        
        filename = os.path.join(self.directory,'test_data/users3_groups_patch.csv')
        # don't check permissions, 
        # because the patch file is setting groups with permissions too
        self._patch_test('user', filename, id_keys_to_check=['username'],
            keys_not_to_check=['permissions']) 
      
    def test4_user_group_permissions(self):
        '''
        Test of a group permission -
        first test that the user in the group has the permission,
        then remove this user from this group and try again
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.info('test4_user_group_permissions...')
        
        self.test2_create_usergroup_with_permissions()
                        
        # Manually assign password to the test user: API does not support 
        # passwords at this time
        username = 'sde4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # 1 read test - should have permission through group 3
        test_resource_uri = '/'.join([BASE_URI, 'apilog'])
        resp = self.api_client.get(
            test_resource_uri, format='json', data={
                'ref_resource_name': 'user'}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        logger.debug('api logs: %r', self.get_content(resp))
        # now patch this user's usergroups, 
        # removing the user from the group 'testgroup3,2'
        # which will remove the permissions for the "apilog" resource as well 
        user_patch = {
            'usergroups': ['testGroup1'] };

        logger.debug(
            'reset this users groups and remove testGroup1: %r', user_patch)
        uri = BASE_URI + '/user' + '/' + username
        resp = self.api_client.patch(
            uri, format='json', data=user_patch, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code <= 204, 
            (resp.status_code, self.get_content(resp)))  
        
        # check the user settings/groups
        uri = BASE_URI + '/user/' + username
        resp = self.api_client.get(
            uri,format='json', data={'includes': '*'}, 
            authentication=self.get_credentials() )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        self.assertTrue(
            new_obj['usergroups'] == ['testGroup1'],
            ('wrong usergroups', new_obj))
        
        # now try again as the updated user:
        test_resource_uri = '/'.join([BASE_URI, 'apilog'])
        logger.info('test user: %r, resource: %r',
            username, test_resource_uri)
        resp = self.api_client.get(
            test_resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password) )
        self.assertTrue(
            resp.status_code in [403], 
            (resp.status_code, self.get_content(resp)))
       
    def test5_usergroup_can_contain_group_permissions(self):
        '''
        Test of an inherited group permission  -
        first test that the user in the group doesn't have the permission,
        then add the user's usergroup to a group with the permission retry.
        - done in prev test: create a new user (as superuser)
        - done in prev test: assign some permissions (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.info('test5_usergroup_can_contain_group_permissions...')
        
        self.test2_create_usergroup_with_permissions()
                        
        # assign password to the test user
        username = 'tester4'
        password = 'testpass1'
        user = User.objects.get(username=username)
        user.set_password(password)
        user.save()

        # 1 read test - user, user's group don't have the permission
        resource_uri = BASE_URI + '/usergroup'
        test_resource_uri = '/'.join([BASE_URI, 'usergroup'])
        logger.info('get: %r', test_resource_uri)
        resp = self.api_client.get(
            test_resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(
            resp.status_code in [403], 
            (resp.status_code, self.get_content(resp)))
        
        # now create a new group, with previous user's group as a member,
        # then add permissions to this new group to read (apilogs)
        # note: double nest the groups also as a test
        usergroup_patch = { API_RESULT_DATA: [
            {
            'name': 'testGroup5',
            'super_groups': ['testGroup3'] },
            {
            'name': 'testGroup6',
            'users': [username],
            'super_groups': ['testGroup5'] },
        ]}
        
        logger.debug('now set the new group: %r', usergroup_patch)
        uri = BASE_URI + '/usergroup'
        
        resp = self.api_client.patch(uri, format='json', 
            data=usergroup_patch, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code <= 204, 
            (resp.status_code, self.get_content(resp)))  
        
        # is it set?
        resp = self.api_client.get(
                uri, format='json', authentication=self.get_credentials())
        new_obj = self.deserialize(resp)
        result, outputobj = find_all_obj_in_list(
            usergroup_patch[API_RESULT_DATA],new_obj[API_RESULT_DATA],
            id_keys_to_check=['name']) #, excludes=keys_not_to_check )
        self.assertTrue(
            result, 
            ('not found', outputobj,'=== objects returned ===', 
                new_obj[API_RESULT_DATA])) 
        
        # is it set-2, does group inherit the permissions?
        for obj in new_obj[API_RESULT_DATA]:
            if obj['name'] == 'testGroup6':
                testGroup6 = obj
            if obj['name'] == 'testGroup3':
                testGroup3 = obj
                
        self.assertTrue(testGroup6 and testGroup3)
        for permission in testGroup3['all_permissions']:
            logger.debug('find permission: %r', permission)
            self.assertTrue(
                permission in testGroup6['all_permissions'], 
                ('could not find permission', permission, 
                     'in testGroup6 permissions', 
                     testGroup6['all_permissions']))
        
        # 2 read test - user has permissions through inherited permissions,
        resp = self.api_client.get(
            test_resource_uri, format='json', data={}, 
            authentication=self.create_basic(username, password ))
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))

    def test6_usergroup_can_contain_group_users(self):
        '''
        Test that group "contains" the subgroup users (through "all_users")
        first test that the group doesn't have the user, then add the subgroup
        with the user, test that the group has the user in "all_users
        - done in prev test: create a new user (as superuser)
        - add/remove users to groups (as superuser)
        '''
        logger.info('test6_usergroup_can_contain_group_users...')
        
        username = 'sde4'
        self.test2_create_usergroup_with_permissions()
                        
        # now create a new group, with a previous group as a sub_group
        usergroup_patch = { API_RESULT_DATA: [
            {
            'resource_uri': 'usergroup/testGroup5',
            'name': 'testGroup5',
            'sub_groups': ['testGroup2'] },
            {
            'resource_uri': 'usergroup/testGroup6',
            'name': 'testGroup6',
            'sub_groups': ['testGroup5'] },
        ]}
        
        logger.debug('now set the new groups: %r', usergroup_patch)
        uri = BASE_URI + '/usergroup'
        resp = self.api_client.patch(uri, format='json', 
            data=usergroup_patch, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code <= 204, 
            (resp.status_code, self.get_content(resp)))  
        
        # is it set?
        data_for_get={ 'includes': '*'}
        resp = self.api_client.get(
                uri + '/testGroup6', 
                format='json', 
                authentication=self.get_credentials(),
                data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        
        self.assertTrue(new_obj['all_users'])
        self.assertTrue(
            username in new_obj['all_users'],
            ('could not find user', username, new_obj['all_users']))
        
        # TODO: could also test that testGroup2 now has super_group=testGroup5
        
        
class JobResource(IResourceTestCase):
    
    def setUp(self):
        super(JobResource, self).setUp()
        settings.BACKGROUND_PROCESSING = True
        
    def tearDown(self):
        IResourceTestCase.tearDown(self)
        
        logger.info('delete users, including: %r', self.username)
        Job.objects.all().delete()
        UserGroup.objects.all().delete()
        UserProfile.objects.all().exclude(username=self.username).delete()
        User.objects.all().exclude(username=self.username).delete()
        ApiLog.objects.all().delete()
        settings.BACKGROUND_PROCESSING = False
    
    def test1_test_job(self):
        '''
        Basic test of the background_job wrapper and the Job Resource
        '''
        logger.info('test1_test_job...')
        JOB = SCHEMA.JOB
        USER = SCHEMA.USER
        JOB_STATE = SCHEMA.VOCAB.job.state
        # Setup: create user
        
        username = 'st1'
        patch_obj = {
            USER.USERNAME: username,
            USER.FIRST_NAME: 'Sally',
            USER.LAST_NAME: 'Tester', 
            USER.EMAIL: 'sally.tester@limstest.com',    
            USER.IS_STAFF: True,
        }
        resource_uri = '/'.join([BASE_URI, USER.resource_name])
        test_uri = '/'.join([resource_uri,patch_obj['username']])
        user_obj = self._create_resource(patch_obj, resource_uri, test_uri)
        logger.info('created user: %r', user_obj)
        
        self.set_user_password(username, self.general_user_password)
        # Now add the needed permission
        user_patch = {
            'resource_uri': 'user/' + username,
            'permissions': ['resource/job/write'] };
        uri = BASE_URI + '/user/' + username
        logger.debug('add permission to user: %r: %r', user_patch, uri)
        resp = self.api_client.patch( uri, 
                    format='json', data=user_patch, 
                    authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code <= 204,
            (resp.status_code, self.get_content(resp)))  
        
        # 1. Invoke the "test_job" from the client
        job_resource_uri = '/'.join([BASE_URI, JOB.resource_name])
        test_background_job_decorated_uri = '/'.join(
            [BASE_URI, JOB.resource_name, 'test_job'])
        patch_obj = {
            'foo': 'bar',
        }
        resp = self.api_client.patch(
            test_background_job_decorated_uri, 
            format='json', data={API_RESULT_DATA: [patch_obj] }, 
            authentication=self.create_basic(username, self.general_user_password))
        self.assertTrue(
            resp.status_code in [200,201,202], 
            (resp.status_code, self.get_content(resp)))
        job_response = self.deserialize(resp)
        logger.info('job created: %r', job_response)
        job_obj = job_response[API_RESULT_META][JOB.resource_name]
        logger.info('job created: %r', job_obj)
        self.assertEqual(JOB_STATE.SUBMITTED,job_obj[JOB.STATE])
        job_id = job_obj[JOB.ID]
        
        # 2. Validate state (pending)
        params = { JOB.ID: job_id }
        job_data = self.get_single_resource(job_resource_uri, params)
        
        self.assertEqual(JOB_STATE.SUBMITTED,job_data[JOB.STATE])

        # 2. Invoke the job using the background processor
        api_client = \
            reports.utils.background_processor.ApiClient(self.username, self.password)
        background_client = \
            reports.utils.background_processor.BackgroundClient(api_client)
        job_service_response = background_client.service(job_id)
        
        logger.info('job_service_response: %r', job_service_response)
        self.assertTrue(API_RESULT_META in job_service_response)
        new_job_data = job_service_response[API_RESULT_META]
        self.assertTrue(JOB.resource_name in new_job_data)
        new_job_data = new_job_data[JOB.resource_name]
        
        self.assertEqual(new_job_data[JOB.RESPONSE_STATUS_CODE], 200)
        self.assertTrue(new_job_data[JOB.STATE], JOB_STATE.COMPLETED)
        self.assertIsNotNone(new_job_data[JOB.DATE_TIME_PROCESSING])
        self.assertIsNotNone(new_job_data[JOB.DATE_TIME_COMPLETED])
        self.assertTrue(JOB.RESPONSE_CONTENT in new_job_data)
        response_content = json.loads(new_job_data[JOB.RESPONSE_CONTENT])
        self.assertEquals(response_content, { 'test_job': 'created!'})
        
    
    