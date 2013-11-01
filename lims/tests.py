from django.test.testcases import TestCase
from lims.api import CSVSerializer, CsvBooleanField

import logging

logger = logging.getLogger(__name__)


def is_boolean(field):
    if type(field) == bool: return True
    
    if isinstance(field, basestring):
        val = field.lower()
        if val=='true' or val=='false':
            return True
    return False
    
csvBooleanField = CsvBooleanField()

# equivocal, and other equivalency methods are for end-to-end testing of the API with
# values read from csv or json, submitted, then read back as json from the API.
def equivocal(val1, val2):
    '''
    For testing equality from the CSV input to the JSON response of what is created in the DB
    obj1 input
    obj2 output 
    '''
    if val1 == val2:
        return True, ('val1', val1, 'val2', val2 )
    
    if not val1:
        if val2:
            if len(val2) > 0:
                return False, ('val1', val1, 'val2', val2 )
        return True, ('val1', val1, 'val2', val2 )
    
    if isinstance(val1, basestring):
        val1 = str(val1)
        val2 = str(val2)
        if val1 != val2:
            if is_boolean(val1) and csvBooleanField.convert(val1) == csvBooleanField.convert(val2):
                return True, ('val1', val1, 'val2', val2 )
            return False, ('val1', val1, 'val2', val2 )
                
    else: # better be a list
        assert isinstance(val1, list) and isinstance(val2, list), str(('Must be a list if not a string', 'val1', val1, 'val2', val2))
        for v in val1:
            if not v: 
                if val2 and len(val2) > 0:
                    return False, ('val1', val1, 'list val2 not empty', val2 )
            if v not in val2:
                if v not in [str(v2) for v2 in val2]:
                    return False, ('val1', val1, 'val2', val2 )
    return True, ('val1', val1, 'val2', val2 )
    
#NOTE only works for String comparisons
def assert_obj1_to_obj2(obj1, obj2, keys=['key', 'scope', 'ordinal', 'username', 'name'], excludes=['resource_uri']):
    '''
    For testing equality from the CSV input to the JSON response of what is created in the DB
    obj1 input
    obj2 output 
    '''
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
            return False, ('key not found',key)  
        result, msgs =  equivocal(obj1[key], obj2[key])
        if not result:
            return False, () # don't report for this section, just move along to the next items to test
#        if str(obj1[key]) != str(obj2[key]):
#            return False, None  # don't report for this section
#    logger.info(str(('equal so far, keys to search', keys_to_search)))

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
                    logger.warn(str(('using imprecise uri matching, equivocating: ', val1, val2)))
    
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


class BaselineTest(TestCase):

    def test_0equivocal(self):
        logger.info(str(('====== test_0equivocal ====')))
        
        test_true = [
            ['',''],
            ['', None],
            ['1','1'],
            ['1',1],        # integer strings can come back as integers (if integerfield is defined)
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
            logger.info(str(('csv serialization equivocal truthy test', a, b)))
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
            logger.info(str(('csv serialization equivocal falsy test', a, b)))
            result, msgs = equivocal(a,b)
            self.assertFalse(result, msgs)
        logger.info(str(('====== test_0equivocal done ====')))
    
    def test_1assert_obj1_to_obj2(self):
        logger.info(str(('====== test_1assert_obj1_to_obj2 ====')))
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
        logger.info(str(('====== test_1assert_obj1_to_obj2 done ====')))

    def test_2find_obj_in_list(self):
        logger.info(str(('====== test_2find_obj_in_list ====')))
        
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
        logger.info(str((result,msgs)))
        self.assertFalse(result, str((msgs)))

        obj_list = [ obj1a, obj3 ]
        result, msgs = find_obj_in_list(obj1, obj_list)
        self.assertFalse(result, msgs)
        logger.info(str(('====== test_2find_obj_in_list done ====')))

    def test_3find_all_obj_in_list(self):
        logger.info(str(('====== test_3find_all_obj_in_list ====')))
        
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
        logger.info(str(('====== test_3find_all_obj_in_list done ====')))

    def test_4user_example(self):
        logger.info(str(('====== test_4user_example ====')))
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
        item = {'login_id': 'jt1', 'first_name': 'Joe', 'last_name': 'Tester', 'email': 'joe.tester@limstest.com'}
        result, obj = find_obj_in_list(item, bootstrap_items)
        logger.info(str((result, obj)))
        
        self.assertTrue(obj['email'] == 'joe.tester@limstest.com')
        
        logger.info(str(('====== test_4user_example done ===='))) # if this passes, then why does line 431 fail in reports.tests.py?
        
        
        
        
        
        
        
        
        
        
        
        
        