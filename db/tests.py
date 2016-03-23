from __future__ import unicode_literals

import cStringIO
import json
import logging
import os
import sys

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.core.urlresolvers import resolve
from django.http.response import StreamingHttpResponse
from django.test import TestCase
from django.test.client import encode_multipart, BOUNDARY, MULTIPART_CONTENT
from django.utils.timezone import now
from tastypie.test import ResourceTestCase, TestApiClient

from db.models import Reagent, Substance, Library, ScreensaverUser, \
    UserChecklistItem, AttachedFile, ServiceActivity, Screen, Well
import db.models
from db.support import lims_utils
from db.support import lims_utils
from db.test.factories import LibraryFactory, ScreenFactory, ScreensaverUserFactory
from reports.dump_obj import dumpObj
from reports.models import ApiLog, UserProfile, UserGroup
from reports.serializers import CSVSerializer, XLSSerializer, LimsSerializer
from reports.sqlalchemy_resource import SqlAlchemyResource
from reports.tests import IResourceTestCase, equivocal
from reports.tests import assert_obj1_to_obj2, find_all_obj_in_list, find_obj_in_list
import reports.tests
import reports.utils.log_utils

logger = logging.getLogger(__name__)


BASE_URI = '/db/api/v1'
BASE_REPORTS_URI = '/reports/api/v1'
import db; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__))
BASE_URI_DB = '/db/api/v1'

class DBResourceTestCase(IResourceTestCase):
    """
    Serves as a base class for tests, and, through _bootstrap_init_files, 
    as the setup for this module.
    """
    def __init__(self,*args,**kwargs):
    
        super(DBResourceTestCase, self).__init__(*args,**kwargs)
        self.directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')

    def _bootstrap_init_files(self):
        
        super(DBResourceTestCase, self)._bootstrap_init_files()
        self.directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        serializer=CSVSerializer() 
        testApiClient = TestApiClient(serializer=serializer) 
        input_actions_file = os.path.join(self.directory, 'api_init_actions.csv')
        logger.info('open input_actions file: %r', input_actions_file)
        self._run_api_init_actions(input_actions_file)
    
    def create_library(self, data=None):

        input_data = LibraryFactory.attributes()
        if data:
            input_data.update(data)
        resource_uri = '/'.join([BASE_URI_DB, 'library'])
        test_uri = '/'.join([resource_uri,input_data['short_name']])
        
        return self._create_resource(input_data,resource_uri,test_uri)

    def create_screen(self, data=None):
        
        input_data = ScreenFactory.attributes()
        if data:
            input_data.update(data)
        resource_uri = '/'.join([BASE_URI_DB, 'screen'])
        _data_for_get = { 
            'limit': 0,
            'includes': '*'
        }
        logger.info('screen input_data to create: %r', input_data)

        logger.info('post to %r...', resource_uri)
        resp = self.api_client.post(
            resource_uri, format='json', data=input_data, 
            authentication=self.get_credentials())
        new_obj = self.deserialize(resp)
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code, self.get_content(resp)))
        
        resource_test_uri = resource_uri + '/' + new_obj['facility_id']        
        logger.info('get from %r... %r', resource_test_uri, _data_for_get)
        resp = self.api_client.get(
            resource_test_uri, format='json', 
            authentication=self.get_credentials(), data=_data_for_get)
        self.assertTrue(
            resp.status_code in [200,201], 
            (resp.status_code,self.get_content(resp)))
        new_obj = self.deserialize(resp)
        if 'objects' in new_obj:
            self.assertEqual(len(new_obj['objects']),1,
                'created resource test_uri: %r, returns multiple objects: %r'
                % (resource_test_uri,new_obj))
            new_obj = new_obj['objects'][0]
        logger.debug('created obj: %r', new_obj)
        result,msg = assert_obj1_to_obj2(input_data,new_obj)
        self.assertTrue(result, msg)
        logger.debug('item created: %r', new_obj)
        return new_obj

    def create_screensaveruser(self, data=None):
        
        input_data = ScreensaverUserFactory.attributes()
        if data:
            input_data.update(data)
        resource_uri = '/'.join([BASE_URI_DB, 'screensaveruser'])
        test_uri = '/'.join([resource_uri,input_data['username']])
        
        return self._create_resource(input_data,resource_uri,test_uri)
        
        
def setUpModule():

    logger.info('=== setup Module')
    
    keepdb = False
    reinit_metahash = False
    if len(sys.argv) > 1:
        for arg in sys.argv:
            print 'arg: ', arg
            if 'keepdb' in arg:
                keepdb = True
            if 'reinit_metahash' in arg:
                reinit_metahash = True
    if reinit_metahash or not keepdb:
        testContext = DBResourceTestCase(methodName='_bootstrap_init_files')
        testContext.setUp()
        testContext._bootstrap_init_files()
    else:
        print 'skip database metahash initialization when using keepdb'

    logger.info('=== setup Module done')

def tearDownModule():

    logger.info('=== teardown Module')
    
    # FIXME: close the sqlalchemy bridge connection on tearDown
    # - the right solution probably requires a custom TestRunner
    # This does not work:
    # bridge.get_engine().connect().close()
    # bridge.get_engine().dispose()
    # bridge = None
    logger.info('=== teardown Module done')

class LibraryResource(DBResourceTestCase):

    def setUp(self):

        super(LibraryResource, self).setUp()

    def tearDown(self):

        DBResourceTestCase.tearDown(self)
        logger.info('delete library resources')
        Library.objects.all().delete()
        Well.objects.all().delete()
        ApiLog.objects.all().delete()

    def test1_create_library(self):
        
        resource_uri = BASE_URI_DB + '/library'
        
        library1 = self.create_library({ 
            'short_name': 'testlibrary1','start_plate': '1534', 
            'end_plate': '1534', 'plate_size': '384' })
        library2 = self.create_library({ 
            'short_name': 'testlibrary2','start_plate': '1535', 
            'end_plate': '1537', 'plate_size': '384' })
        
        # now find the library wells
        resource_uri = '/'.join([
            BASE_URI_DB,'library',library2['short_name'],'well'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
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
            len(new_obj['objects']), expected_count, 
            str(('wrong number of wells',len(new_obj['objects']), expected_count)))
        
        logger.info('examine the wells created...')
        index = 0
        platesize = 384
        plate = 1535        
        substance_ids = set()
        # Examine first plate only for speed
        for j in range(384):
            well_name = lims_utils.well_name_from_index(j, platesize)
            well_id = lims_utils.well_id(plate,well_name)
            well_search = {'well_id': well_id}
            result, well = find_obj_in_list(well_search, new_obj['objects'])
            self.assertTrue(result, well)
                
    def test10_create_library_copy(self):
        
        library = self.create_library({
            'start_plate': 1000, 
            'end_plate': 1005})

        logger.info('create library copy...')
        input_data = {
            'library_short_name': library['short_name'],
            'name': "A",
            'usage_type': "library_screening_plates"
        }        
        resource_uri = BASE_URI_DB + '/librarycopy'
        resource_test_uri = '/'.join([
            resource_uri,input_data['library_short_name'],input_data['name']])
        library_copy = self._create_resource(
            input_data, resource_uri, resource_test_uri)
        
        logger.info('verify plates created...')
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
        start_plate = int(library['start_plate'])
        end_plate = int(library['end_plate'])
        number_of_plates = end_plate-start_plate+1
        self.assertEqual(len(new_obj['objects']),number_of_plates)

        for obj in new_obj['objects']:
            self.assertEqual(library_copy['name'],obj['copy_name'])
            plate_number = int(obj['plate_number'])
            self.assertTrue(
                plate_number>=start_plate and plate_number<=end_plate,
                'wrong plate_number: %r' % obj)

    def test4_create_library_invalids(self):
        '''
        Test the schema "required" validations 
        ''' 
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
        
                obj = self.deserialize(resp)
                self.assertTrue(
                    key in obj, 
                    'Error: response error not found: %r, obj: %r' %(key, obj))

        # test 2                
        library_item = LibraryFactory.attributes()
        library_item['library_name'] = 'invalid & name'
        
        logger.debug('try to create an invalid library name: %r', library_item)
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        
        logger.debug('response.content.library message: %r', 
                         getattr(resp, 'content'))
        
        obj = self.deserialize(resp)
        self.assertTrue(
            'library_name' in obj, 
            'content should have an entry for the faulty "name": %r' % obj)

        # test 3
        library_item = LibraryFactory.attributes()
        library_item['library_type'] = 'invalid_type'
        
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(
            resp.status_code in [400], 
            (resp.status_code, self.get_content(resp)))
        obj = self.deserialize(resp)
        self.assertTrue(
            'library_type' in obj, 
            'content should have an entry for the faulty "name": %r' % obj)

        # TODO: test regex and number: min/max, vocabularies
        
    def test61_small_molecule_delete_compound_name(self):
        # TODO: test that sm.compound names can be removed
        pass
    
    def test6_load_small_molecule_file(self):
        
        library_item = self.create_library({ 
            'start_plate': '1536', 
            'end_plate': '1536', 
            'plate_size': '384',
            'screen_type': 'small_molecule' })

        resource_name = 'well'
        resource_uri = '/'.join([BASE_URI_DB,'library', library_item['short_name'],resource_name])
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_small_molecule.sdf'

        data_for_get={}
        data_for_get.setdefault('limit', 0)
        data_for_get.setdefault('includes', ['*'])
        data_for_get.setdefault('HTTP_ACCEPT', 'chemical/x-mdl-sdfile' )

        with open(filename) as input_file:
            
            input_data = self.serializer.from_sdf(input_file.read())
            input_data = input_data['objects']
            expected_count = 8
            self.assertEqual(
                len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            logger.debug('======Submitting patch file: %r, uri: %r', 
                filename, resource_uri)
            resp = self.api_client.put(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
        
        # NOTE: get the library "reagents" instead of the wells: reagents have a 
        # one-to-many reln to wells, so the well endpoint returns only first;
        # also, the well resource is not the linked resource type, and does not have
        # the reagent fields.  
        resource_name = 'reagent'
        reagent_resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        resp = self.api_client.get(
            reagent_resource_uri, format='sdf', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        returned_data = new_obj['objects']
        expected_count = 384
        self.assertEqual(
            len(returned_data), expected_count, 
            ('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count,
                'returned_data',returned_data))
        
        # test 1: test well keys
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        
        substance_ids = set()
        for inputobj in input_data:
            # first just search for the well
            search = { 'well_id': 
                lims_utils.well_id(
                    inputobj['plate_number'],inputobj['well_name']) }
            result, outputobj = find_obj_in_list(
                search,returned_data) #, excludes=excludes )
            self.assertTrue(
                result, 
                ('not found', search,outputobj,'=== objects returned ===', 
                      returned_data ) ) 
            # second look at the found item
            expected_data = { 
                key: inputobj[key] for key in fields.keys() if key in inputobj }
            result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
            self.assertTrue(
                result, (msgs, 'input', expected_data, 'output', outputobj))

            substance_id = outputobj['substance_id']
            self.assertTrue(
                substance_id not in substance_ids, 
                ('substance_id not unique', substance_id, substance_ids))
            substance_ids.add(substance_id)
                    
        logger.debug('==== done: test6_load_small_molecule_file: create =====')
    
        # test 2: update some wells, check results, and api logs
        filename = ( APP_ROOT_DIR 
            + '/db/static/test_data/libraries/clean_data_small_molecule_update.sdf')

        data_for_get={}
        data_for_get.setdefault('limit', 0)
        data_for_get.setdefault('includes', ['*'])
        data_for_get.setdefault('HTTP_ACCEPT', 'chemical/x-mdl-sdfile' )

        with open(filename) as input_file:
            
            input_data = self.serializer.from_sdf(input_file.read())
            input_data = input_data['objects']
            expected_count = 4
            self.assertEqual(
                len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
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
            returned_data = new_obj['objects']
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
                self.assertTrue('library_well_type' in outputobj,
                    'library_well_type is missing in %r' % outputobj)
                self.assertTrue(outputobj['library_well_type']=='experimental',
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
                len(new_obj['objects']), expected_count , 
                str((len(new_obj['objects']), expected_count, new_obj)))

            # 2. check the apilogs - wells
            resource_uri = BASE_REPORTS_URI + '/apilog' 
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 'limit': 0, 'ref_resource_name': 'well' })
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            logs = new_obj['objects']
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
                    self.assertTrue(
                        set(json.loads(log['diff_keys']))==
                        set([
                            "pubchem_cid", "vendor_identifier", 
                            "vendor_batch_id", "compound_name", "smiles"]),
                        'diff_keys: %r should equal %r' % (
                            log['diff_keys'],[
                                "pubchem_cid", "vendor_identifier", 
                                "vendor_batch_id", "compound_name", "smiles"]))
            # TODO: check parent_log - library log/ version
    
    def test7_load_sirnai(self):

        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_rnai.xlsx'
        
        library_item = self.create_library({ 
            'start_plate': 50001,  
            'end_plate': 50001, 
            'plate_size': '384',
            'screen_type': 'rnai' })
        resource_uri = BASE_URI_DB + '/library'
        
        self._load_xls_reagent_file(filename,library_item, 5,384 )

    #     def test7a_sirna_validations(self):
    #         # TODO: test validations
    #         pass
    #     

    def test8_sirna_duplex(self):
        
        filename = ( APP_ROOT_DIR 
            + '/db/static/test_data/libraries/clean_rnai_duplex_50440_50443.xlsx')
        
        # create the duplex library
        library_item = self.create_library({ 
            'start_plate': 50440,  
            'end_plate': 50443, 
            'plate_size': '384',
            'screen_type': 'rnai' })
        
        self._load_xls_reagent_file(filename,library_item, 532, 1536 )
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_rnai_pool.xlsx'
        
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
        
        filename = ( APP_ROOT_DIR 
            + '/db/static/test_data/libraries/clean_data_natural_product.xlsx')
        
        library_item = self.create_library({ 
            'start_plate': 2037,  
            'end_plate': 2037, 
            'plate_size': '384',
            'library_type': 'natural_products' })
        
        self._load_xls_reagent_file(filename,library_item, 352, 384 )

    def _load_xls_reagent_file(
            self,filename,library_item, expected_in, expected_count):

        start_plate = library_item['start_plate']
        end_plate = library_item['end_plate']

        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        
        data_for_get={}
        data_for_get.setdefault('limit', 0)
        data_for_get.setdefault('includes', ['*'])
        data_for_get.setdefault('HTTP_ACCEPT', 
            'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' )
        xls_serializer = XLSSerializer()
        
        with open(filename) as input_file:
            ### NOTE: submit the data as an object, because the test framework
            ### will convert it to the target format.
            input_data = self.serializer.from_xlsx(input_file.read())
            input_data = input_data['objects']
        
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
            reagent_resource_uri, format='xlsx', authentication=self.get_credentials(), 
            data=data_for_get)
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        returned_data = new_obj['objects']
        self.assertEqual(
            len(returned_data), expected_count, 
            str(('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count )))

        # 1. test well keys
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        
        substance_ids = set()
        for inputobj in input_data:
            # first just search for the well
            search = { 
                'well_id': lims_utils.well_id(
                    inputobj['plate_number'],inputobj['well_name']) }
            result, outputobj = find_obj_in_list(search,returned_data)
            self.assertTrue(
                result, 
                ('not found', search,outputobj,'=== objects returned ===', 
                      returned_data )) 

            # second look at the found item
            expected_data = { key: inputobj[key] for key in fields.keys() 
                if key in inputobj }
            logger.debug('checking: expected: %r', expected_data )
            result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
            self.assertTrue(result, (msgs, expected_data, outputobj))
            
            substance_id = outputobj['substance_id']
            self.assertTrue(
                substance_id not in substance_ids, 
                ('substance_id not unique', substance_id, substance_ids))
            substance_ids.add(substance_id)
                    

class ScreenResource(DBResourceTestCase):
        
    def setUp(self):
        logger.debug('============== ScreenResource setup ============')
        super(ScreenResource, self).setUp()

    def tearDown(self):
        DBResourceTestCase.tearDown(self)
        
        logger.info('delete resources')
        Screen.objects.all().delete()
        Library.objects.all().delete()
        ApiLog.objects.all().delete()
        
    def test2_create_library_screening(self):
        
        logger.info('create users...')
        self.screening_user = self.create_screensaveruser({ 
            'username': 'screening1'
        })
        self.admin_user = self.create_screensaveruser({ 
            'username': 'adminuser',
            'is_superuser': True
        })
        
        logger.info('create library...')
        library1 = self.create_library({
            'start_plate': 1000, 
            'end_plate': 1005,
            'screen_type': 'small_molecule' })

        library2 = self.create_library({
            'start_plate': 2000, 
            'end_plate': 2040,
            'screen_type': 'small_molecule' })

        logger.info('create library copy...')
        input_data = {
            'library_short_name': library1['short_name'],
            'name': "A",
            'usage_type': "library_screening_plates"
        }        
        resource_uri = BASE_URI_DB + '/librarycopy'
        resource_test_uri = '/'.join([
            resource_uri,input_data['library_short_name'],input_data['name']])
        library_copy1 = self._create_resource(
            input_data, resource_uri, resource_test_uri)
        logger.info('created: %r', library_copy1)
 
        input_data['library_short_name'] = library2['short_name']
        resource_test_uri = '/'.join([
            resource_uri,input_data['library_short_name'],input_data['name']])
        library_copy2 = self._create_resource(
            input_data, resource_uri, resource_test_uri)
        logger.info('created: %r', library_copy2)
        logger.info('create screen...')        
        screen = self.create_screen({
            'screen_type': 'small_molecule'
            })

        lps_format = '{library_short_name}:{name}:{{start_plate}}-{{end_plate}}'
        library_plates_screened = [
            lps_format.format(**library_copy1).format(**library1),
            lps_format.format(**library_copy2).format(
                start_plate=library2['start_plate'],
                end_plate=int(library2['start_plate']+10))
        ]
        
        library_screening_input = {
            'screen_facility_id': screen['facility_id'],
            'date_of_activity': "2008-01-18",
            'assay_protocol':'test assay protocol',
            'assay_protocol_type': '',
            'is_for_external_library_plates': False,
            'library_plates_screened': library_plates_screened,
            'number_of_replicates': 2,
            'performed_by_username': self.admin_user['username'],
            'volume_transferred_per_well_from_library_plates': "0.000006000",
            'volume_transferred_per_well_to_assay_plates': "0.000002000"
        }
        resource_uri = BASE_URI_DB + '/libraryscreening'
        resource_test_uri = BASE_URI_DB + '/libraryscreening'
        data_for_get = {
            'screen_facility_id__eq': screen['facility_id']
        }

        # test validations
        key = 'screen_facility_id' 
        msg = 'input missing a %r should fail' % key
        logger.info('test %r', msg)
        input = library_screening_input.copy()
        del input['screen_facility_id'] 
        errors, resp = self._create_resource(
            input, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(
            key in errors, 'test: %s, not in errors: %r' %(msg,errors))
        
        key = 'library_plates_screened' 
#         value = [ 'short_name_good:copy_name_good:1000-a3000' ]
        # create a library_plate_range with invalid library name:
        # even though the plate range is correct, this should fail,
        # and the transaction should cancel the creation of the libraryscreening
        value = [lps_format.format(library_short_name='**',name='A').format(**library1)]
        msg = 'invalid format for %r should fail' % key
        logger.info('test %r', msg)
        input = library_screening_input.copy()
        input[key] =  value
        errors, resp = self._create_resource(
            input, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(
            key in errors, 'test: %s, not in errors: %r' %(msg,errors))

        key = 'library_plates_screened' 
        value = [ lps_format.format(**library_copy1).format(**{
            'start_plate': library1['start_plate'],
            'end_plate': int(library1['start_plate'])-1 }) ]
        msg = 'invalid plate range in %r for  %r should fail' % (value,key)
        logger.info('test %r', msg)
        input = library_screening_input.copy()
        input[key] =  value
        errors, resp = self._create_resource(
            input, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(
            key in errors, 'test: %s, not in errors: %r' %(msg,errors))

        key = 'library_plates_screened' 
        value = [ 
            lps_format.format(**library_copy1).format(**{
                'start_plate': library2['start_plate'],
                'end_plate': int(library2['start_plate'])+2 }),
            lps_format.format(**library_copy1).format(**{
                'start_plate': library2['start_plate']+1,
                'end_plate': int(library2['start_plate'])+4 }),
        ]
        msg = 'overlapping plate ranges in %r for  %r should fail' % (value,key)
        logger.info('test %r', msg)
        input = library_screening_input.copy()
        input[key] =  value
        errors, resp = self._create_resource(
            input, resource_uri, resource_test_uri, expect_fail=True)
        self.assertTrue(resp.status_code==400, msg)
        self.assertTrue(
            key in errors, 'test: %s, not in errors: %r' %(msg,errors))

        logger.info('test valid input...')
        library_screening = self._create_resource(
            library_screening_input, resource_uri, resource_test_uri)

        # test update
        # add a plate_range
        input = { 
            'activity_id': library_screening['activity_id'],
            'library_plates_screened': library_screening['library_plates_screened'] 
            }
        input['library_plates_screened'].append(
            lps_format.format(**library_copy2).format(
                start_plate=library2['start_plate']+15,
                end_plate=int(library2['start_plate']+20)))
        library_screening = self._create_resource(
            input, resource_uri, resource_test_uri)

        # test update
        # delete a plate_range
        input = { 
            'activity_id': library_screening['activity_id'],
            'library_plates_screened': library_screening['library_plates_screened'][:-1] 
            }
        logger.info('input: %r', input)
        library_screening = self._create_resource(
            input, resource_uri, resource_test_uri)
        
    
    def test1_create_screen(self):
        logger.info('==== test1_create_screen =====')
        
        screen_item = self.create_screen()
        
        self.assertTrue(
            'facility_id' in screen_item, 
            'the facility_id was not created')
        

class MutualScreensTest(DBResourceTestCase,ResourceTestCase):
    
    def test_mutual_positives_to_screen(self):
        pass
        # create two screens
        # create two screen results
        # set data sharing levels
        # find assay wells that overlap
        
class ScreensaverUserResource(DBResourceTestCase):
        
    def setUp(self):
        logger.debug('============== ScreensaverUserResource setup ============')
        super(ScreensaverUserResource, self).setUp()

    def tearDown(self):
        DBResourceTestCase.tearDown(self)
        
        logger.info('delete resources')
        UserChecklistItem.objects.all().delete()
        AttachedFile.objects.all().delete()
        ServiceActivity.objects.all().delete()
        ScreensaverUser.objects.all().delete()
        
        UserGroup.objects.all().delete()
        UserProfile.objects.all().delete()
        User.objects.all().delete()
        
        ApiLog.objects.all().delete()

    def test01_create_user_iccbl(self):
        logger.debug('=== test01_create_user_iccbl')

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
        resource_test_uri = '/'.join([resource_uri,simple_user_input['ecommons_id']])
        created_user = self._create_resource(
            simple_user_input, resource_uri, resource_test_uri)
        self.assertTrue(
            simple_user_input['ecommons_id']==created_user['username'],
            'username should equal the ecommons id if only ecommons is provided: %r, %r' 
                % (simple_user_input,created_user))
        
    def test0_create_user(self):
        
        self.user1 = self.create_screensaveruser({ 'username': 'st1'})
        self.screening_user = self.create_screensaveruser({ 'username': 'screening1'})
        self.admin_user = self.create_screensaveruser({ 'username': 'adminuser'})

        # create an admin
        
        patch_obj = { 'objects': [
            {
                'username': 'adminuser',
                'is_superuser': True
            }]
        }
        resource_uri = BASE_REPORTS_URI + '/user'

        try:       
            resp = self.api_client.patch(resource_uri, 
                format='json', data=patch_obj, 
                authentication=self.get_credentials())
            self.assertTrue(
                resp.status_code in [200,201,202], 
                (resp.status_code, self.get_content(resp)))
        except Exception, e:
            logger.exception('on patching adminuser %s' % patch_obj)
            raise

    def test1_patch_usergroups(self):
        
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
                format='json', data=group_patch, authentication=self.get_credentials())
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
                len(new_obj['objects']), len(group_patch['objects']), new_obj)
            
            for i,item in enumerate(group_patch['objects']):
                result, obj = find_obj_in_list(item, new_obj['objects'])
                self.assertTrue(
                    result, 
                    ('bootstrap item not found', item, new_obj['objects']))
                logger.info('item found: %r', obj)        
        except Exception, e:
            logger.exception('on group_patch: %r', group_patch)
            raise

        userpatch = { 'objects': [   
            {
                'username': self.user1['username'],
                'usergroups': ['usergroup1']
            },
            {
                'username': self.screening_user['username'],
                'usergroups': ['usergroup2']
            },
            {
                'username': self.admin_user['username'],
                'usergroups': ['usergroup3']
            },
        ]};
        try:       
            resource_uri = BASE_URI_DB + '/screensaveruser'
            resp = self.api_client.patch(resource_uri, 
                format='json', data=userpatch, authentication=self.get_credentials())
            self.assertTrue(
                resp.status_code in [200,201,202], 
                (resp.status_code, self.get_content(resp)))

            data_for_get={ 'limit': 0 }        
            data_for_get.setdefault('includes', ['*'])
            resp = self.api_client.get(resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get )
            new_obj = self.deserialize(resp)
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            self.assertEqual(len(new_obj['objects']), 3, (new_obj))
            
            for i,item in enumerate(userpatch['objects']):
                result, obj = find_obj_in_list(item, new_obj['objects'])
                self.assertTrue(
                    result, 
                    ('bootstrap item not found', item, new_obj['objects']))
                logger.debug('item found: %r', obj)        
        
        except Exception, e:
            logger.exception('on userpatch %r, %e', userpatch, e)
            raise

    def test2_user_checklist_items(self):
        self.test0_create_user();
        
        test_username = self.user1['username']
        checklist_item_patch = {
            'objects': [
                { 
                    'admin_username': self.admin_user['username'], 
                    'item_group': "mailing_lists_wikis",
                    'item_name': "added_to_iccb_l_nsrb_email_list",
                    'status': "activated",
                    'status_date': "2015-09-02",
                    'username': test_username
                }
            ]}
        
        try:       
            resource_uri = BASE_URI_DB + '/userchecklistitem/%s' % test_username
            resp = self.api_client.patch(resource_uri, 
                format='json', 
                data=checklist_item_patch, authentication=self.get_credentials())
            self.assertTrue(
                resp.status_code in [200,201,202], 
                (resp.status_code, self.get_content(resp)))

            data_for_get={ 'limit': 0 }        
            data_for_get.setdefault('includes', ['*'])
            resp = self.api_client.get(
                resource_uri + '/mailing_lists_wikis/added_to_iccb_l_nsrb_email_list',
                format='json', 
                authentication=self.get_credentials(), data=data_for_get )
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            result,msgs = assert_obj1_to_obj2(
                checklist_item_patch['objects'][0], new_obj)
            self.assertTrue(result,msgs)
        
        except Exception, e:
            logger.exception('on userchecklist')
            raise e
        
        # TODO: checklistitem logs
                    
    def test3_attached_files(self):
        '''
        '''
        self.test0_create_user();

        # 1. test using embedded "contents" field               
        test_username = self.user1['username']
        admin_username = self.admin_user['username']
        attachedfile_item_post = {
            'created_by_username': admin_username, 
            'type': '2009_iccb_l_nsrb_small_molecule_user_agreement', 
            'filename': "test_pasted_text.txt",
            'contents': "This is a test of pasted text\n1\n2\n3\n\n end\n",
            'file_date': '2015-10-10'
            }

        content_type = MULTIPART_CONTENT
        resource_uri = BASE_URI_DB + '/screensaveruser/%s/attachedfiles/' % test_username
        
        authentication=self.get_credentials()
        kwargs = {}
        kwargs['HTTP_AUTHORIZATION'] = authentication
        
        django_test_client = self.api_client.client
        resp = django_test_client.post(resource_uri, 
            content_type=content_type, 
            data=attachedfile_item_post, **kwargs)
        self.assertTrue(
            resp.status_code in [201, 202], 
            (resp.status_code,self.get_content(resp)))
        
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            resource_uri,
            authentication=self.get_credentials(), data=data_for_get )
        self.assertTrue(
            resp.status_code in [200], 
            (resp.status_code, self.get_content(resp)))
        new_obj = self.deserialize(resp)
        logger.info('new obj: %s ' % new_obj)
        result,msgs = assert_obj1_to_obj2(
            attachedfile_item_post, new_obj['objects'][0], 
            excludes=['contents'])
        self.assertTrue(result,msgs)
        af = new_obj['objects'][0]
        uri = '/db/attachedfile/content/%s' % af['attached_file_id']
        try:
            admin_user = User.objects.get(username=admin_username)
            view, args, kwargs = resolve(uri)
            kwargs['request'] = self.api_client.client.request()
            kwargs['request'].user=admin_user
            result = view(*args, **kwargs)
            logger.info('attached_file request view result: %r',result)
            self.assertEqual(
                result.content, attachedfile_item_post['contents'], 
                'download file view returns wrong contents: %r' % result.content)
        except Exception, e:
            logger.info('no file found at: %r', uri)
            raise
    

        # TODO: 2. test attached file from file system
        # TODO: delete attached file
        # TODO: attachedfile logs

    def test4_service_activity(self):
        
        self.test0_create_user();
        
        test_username = self.user1['username']
        admin_username = self.admin_user['username']
        service_activity_post = {
            'serviced_username': test_username,
            'type': "image_analysis",
            'comments': "test",
            'date_of_activity': "2015-10-27",
            'funding_support': "clardy_grants",
            'performed_by_username': admin_username,
        }

        
        try:       
            resource_uri = BASE_URI_DB + '/serviceactivity'
            resp = self.api_client.post(resource_uri, 
                format='json', 
                data=service_activity_post, authentication=self.get_credentials())
            self.assertTrue(
                resp.status_code in [200,201,202], 
                (resp.status_code, self.get_content(resp)))

            data_for_get={ 'limit': 0 }        
            data_for_get.setdefault('includes', ['*'])
            resp = self.api_client.get(
                resource_uri,
                format='json', 
                authentication=self.get_credentials(), data=data_for_get )
            self.assertTrue(
                resp.status_code in [200], 
                (resp.status_code, self.get_content(resp)))
            new_obj = self.deserialize(resp)
            result,msgs = assert_obj1_to_obj2(
                service_activity_post, new_obj['objects'][0])
            self.assertTrue(result,msgs)
        
        except Exception, e:
            logger.exception('on serviceactivity test')
            raise e
        
        # TODO: delete serivceactivity
        # TODO: serviceactivity logs
