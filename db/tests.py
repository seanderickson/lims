import cStringIO
import json
import logging
import os

from django.http.response import StreamingHttpResponse
from django.test import TestCase
from django.utils.timezone import now
import factory
from south.migration import get_migrator
from south.migration.base import Migrations
from south.migration.migrators import FakeMigrator
from tastypie.test import ResourceTestCase, TestApiClient

from db.models import Reagent, Substance, Library
import db.models
from db.support import lims_utils
from db.test.factories import LibraryFactory
from reports.dump_obj import dumpObj
from reports.serializers import CSVSerializer, XLSSerializer
from reports.tests import MetaHashResourceBootstrap
from reports.tests import assert_obj1_to_obj2, find_all_obj_in_list, find_obj_in_list
import reports.tests
from test.factories import *


logger = logging.getLogger(__name__)


BASE_URI = '/db/api/v1'
BASE_REPORTS_URI = '/reports/api/v1'
import db; 
try:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__[0]))
except:
    APP_ROOT_DIR = os.path.abspath(os.path.dirname(db.__path__))
BASE_URI_DB = '/db/api/v1'



class DBMetaHashResourceBootstrap(MetaHashResourceBootstrap):
    
    def setUp(self):
        super(DBMetaHashResourceBootstrap, self).setUp()
        self._bootstrap_init_files()

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
        self._patch_test('vocabularies', filename, id_keys_to_check=['key','scope'])

        
        logger.debug('------------- Done: DBMetaHashResourceBootstrap _bootstrap_init_files -----------------')        
                  

class LibraryResource(DBMetaHashResourceBootstrap):

    def setUp(self):
        logger.debug('============== LibraryResource setup ============')
        
        super(LibraryResource, self).setUp()
        logger.debug('============== LibraryResource setup: begin ============')
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        filename = os.path.join(self.db_directory,'metahash_fields_library.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.library'})
        filename = os.path.join(self.db_directory,'metahash_fields_well.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.well'})

        filename = os.path.join(self.db_directory,'metahash_fields_reagent.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.reagent'})

        filename = os.path.join(self.db_directory,'metahash_fields_smallmoleculereagent.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.smallmoleculereagent'})

        filename = os.path.join(self.db_directory,'metahash_fields_naturalproductreagent.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.naturalproductreagent'})

        filename = os.path.join(self.db_directory,'metahash_fields_silencingreagent.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.silencingreagent'})
        
        logger.debug('============== LibraryResource setup: done ============')

    def test1_create_library(self):
        
        logger.debug(str(('==== test1_create_library =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 'start_plate': '1534', 'end_plate': '1534', 'plate_size': '384' })
        
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))

        # create a second library
        library_item = LibraryFactory.attributes()
        library_item.update({ 'start_plate': '1535', 'end_plate': '1537', 'plate_size': '384' })
         
        logger.debug(str(('item', library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
         
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.info(str(('--------resp to get:', resp.status_code)))
        # self.assertValidJSONResponse(resp)
        self.assertTrue(resp.status_code in [200], str((resp.status_code)))
        new_obj = self.deserialize(resp)
        logger.info(str(('resp', new_obj)))
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(library_item, new_obj['objects'])
        self.assertTrue(
            result, str(('library_item', obj, 
                         library_item, new_obj['objects'])))
        logger.debug(str(('item found', obj)))
        
        # now find the library wells
        resource_uri = '/'.join([BASE_URI_DB,'library',library_item['short_name'],'well'])
        logger.debug(str(('GET', resource_uri)))
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.debug(str(('--------resp to get:', resp.status_code, resp['Content-Type'], resp)))
        self.assertTrue(resp.status_code in [200], str((resp.status_code)))
        self.assertTrue(resp['Content-Type'].startswith('application/json'))
        new_obj = self.deserialize(resp)
        expected_count = 384*3
        self.assertEqual(len(new_obj['objects']), expected_count, 
            str(('wrong number of wells',len(new_obj['objects']), expected_count)))
        
        # now examine the wells created
        index = 0
        platesize = 384
        plate = 1535        
        substance_ids = set()
        from db.support import lims_utils
        for i in range(3):
            for j in range(384):
                _plate = plate + i
                well_name = lims_utils.well_name_from_index(j, platesize)
                well_id = lims_utils.well_id(_plate,well_name)
                well_search = {'well_id': well_id}
                result, well = find_obj_in_list(well_search, new_obj['objects'])
                self.assertTrue(result, well)
                
        logger.debug(str(('==== done: test1_create_library =====')))

    def test4_create_library_invalids(self):
        '''
        Test the schema "required" validations 
        
        - todo: this can be a template for other endpoints
        ''' 
        logger.debug(str(('==== test4_create_library_invalids =====')))

        library_resource = self.get_resource_from_server('library')
        fields = library_resource['schema']['fields']

        # make sure the default works
        resource_uri = BASE_URI_DB + '/library'
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

        # another test                
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

        # another test
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

        # TODO: test regex and number: min/max, vocabularies
        logger.debug(str(('==== done: test4_create_library_invalids =====')))
        
    def test6_load_small_molecule_file(self):
        logger.debug(str(('==== test6_load_small_molecule_file =====')))
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 'start_plate': '1536', 'end_plate': '1536', 'plate_size': '384' })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))

        logger.debug(str(('==== start: test6_load_small_molecule_file =====')))

        resource_name = 'well'
        resource_uri = '/'.join([BASE_URI_DB,'library', library_item['short_name'],resource_name])
        
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_small_molecule.sdf'

        data_for_get={}
        data_for_get.setdefault('limit', 0)
        data_for_get.setdefault('includes', ['*'])
        data_for_get.setdefault('HTTP_ACCEPT', 'chemical/x-mdl-sdfile' )

        with open(filename) as input_file:
            
            ### NOTE: submit the data as an object, because the test framework
            ### will convert it to the target format.
            input_data = self.serializer.from_sdf(input_file.read())
            input_data = input_data['objects']
            #logger.info(str(('===== input data', input_data)))
        
            expected_count = 8
            self.assertEqual(len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            logger.debug(str(('======Submitting patch...', filename, resource_uri)))
        
            resp = self.api_client.put(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [200, 204], 
                str((resp.status_code, self.deserialize(resp))))
        
        logger.debug(str(('check patched data for',resource_name,
            'execute get on:',resource_uri)))

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

        logger.debug(str(('--------resp to get:', resp.status_code)))
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        new_obj = self.deserialize(resp)
        returned_data = new_obj['objects']
        expected_count = 384
        self.assertEqual(len(returned_data), expected_count, 
            str(('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count,
                'returned_data',returned_data)))
        
        # 1. test well keys
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        logger.debug(str(('=== well fields', fields.keys())))
        
        substance_ids = set()
        for inputobj in input_data:
            # first just search for the well
            search = { 'well_id': lims_utils.well_id(inputobj['plate_number'],inputobj['well_name']) }
            result, outputobj = find_obj_in_list(
                search,returned_data) #, excludes=excludes )
            self.assertTrue(result, 
                str(('not found', search,outputobj,'=== objects returned ===', 
                      returned_data )) ) 
            logger.debug(str(('found', search)))
            # second look at the found item
            expected_data = { key: inputobj[key] for key in fields.keys() if key in inputobj }
            result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
            logger.debug(str((result, msgs)))
            self.assertTrue(result, str((msgs, 'input', expected_data, 'output', outputobj)))

            substance_id = outputobj['substance_id']
            self.assertTrue(substance_id not in substance_ids, 
                str(('substance_id not unique', substance_id, substance_ids)))
            substance_ids.add(substance_id)
                    
        logger.debug(str(('==== done: test6_load_small_molecule_file =====')))
    
        ########
        # Next: update some wells, check results, and api logs
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_small_molecule_update.sdf'

        data_for_get={}
        data_for_get.setdefault('limit', 0)
        data_for_get.setdefault('includes', ['*'])
        data_for_get.setdefault('HTTP_ACCEPT', 'chemical/x-mdl-sdfile' )

        with open(filename) as input_file:
            
            input_data = self.serializer.from_sdf(input_file.read())
            input_data = input_data['objects']
            #logger.info(str(('===== input data', input_data)))
        
            expected_count = 4
            self.assertEqual(len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            logger.debug(str(('======Submitting patch...', filename, resource_uri)))
        
            resp = self.api_client.put(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [200, 204], 
                str((resp.status_code, self.deserialize(resp))))
        
            logger.debug(str(('check updated/patched data for',resource_name,
                'execute get on:',resource_uri)))
    
            resource_name = 'well'
            resource_uri = '/'.join([BASE_URI_DB,'library', library_item['short_name'],resource_name])
            resp = self.api_client.get(
                reagent_resource_uri, format='sdf', authentication=self.get_credentials(), 
                data=data_for_get)
    
            logger.debug(str(('--------resp to get:', resp.status_code)))
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
            new_obj = self.deserialize(resp)
            returned_data = new_obj['objects']
            expected_count = 384
            self.assertEqual(len(returned_data), expected_count, 
                str(('returned_data of ',filename,'found',
                    len(returned_data), 'expected',expected_count,
                    'returned_data',returned_data)))
        
            # 1. test the specific wells
            well_ids_to_check = input_data
    
            for update_well in well_ids_to_check:
                search = { 'well_name': update_well['well_name'], 'plate_number': update_well['plate_number']}
                result, outputobj = find_obj_in_list(
                    search,returned_data) #, excludes=excludes )
                self.assertTrue(result, 
                    str(('not found', search,outputobj,'=== objects returned ===', 
                          returned_data )) ) 
                logger.debug(str(('found', search)))
                
                result, msgs = assert_obj1_to_obj2(update_well, outputobj)
                logger.debug(str((result, msgs)))
                self.assertTrue(result, str((msgs, update_well, outputobj)))

            # 2. check the apilogs - library
            resource_uri = BASE_REPORTS_URI + '/apilog' #?ref_resource_name=record'
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 'limit': 0, 'ref_resource_name': 'library', 'key': library_item['short_name'] })
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
            new_obj = self.deserialize(resp)
            logger.debug(str(('===library apilogs:', json.dumps(new_obj))))
            
            expected_count = 3 # create, post, update
            self.assertEqual( len(new_obj['objects']), expected_count , 
                str((len(new_obj['objects']), expected_count, new_obj)))

            
            # 2. check the apilogs
            resource_uri = BASE_REPORTS_URI + '/apilog' #?ref_resource_name=record'
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 'limit': 0, 'ref_resource_name': 'well' })
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
            new_obj = self.deserialize(resp)
            logger.debug(str(('===apilogs:', json.dumps(new_obj))))
            
            # look for 12 logs; 8 for create, 4 for update
            expected_count = 12
            self.assertEqual( len(new_obj['objects']), expected_count , 
                str((len(new_obj['objects']), expected_count)))
        
    def test7_load_sirnai(self):

        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_rnai.xlsx'
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': 50001,  
            'end_plate': 50001, 
            'plate_size': '384',
            'screen_type': 'rnai' })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        self._load_xls_reagent_file(filename,library_item, 5,384 )

#     def test7a_sirna_validations(self):
#         # TODO: test validations
#         pass
#     
    def test8_sirna_duplex(self):
        
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_rnai_duplex_50440_50443.xlsx'
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': 50440,  
            'end_plate': 50443, 
            'plate_size': '384',
            'screen_type': 'rnai' })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        self._load_xls_reagent_file(filename,library_item, 532, 1536 )

        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_rnai_pool.xlsx'

        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': 50439,  
            'end_plate': 50439, 
            'plate_size': '384',
            'screen_type': 'rnai',
            'is_pool': True })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        self._load_xls_reagent_file(filename,library_item, 133, 384 )
        
        ## TODO: test the duplex wells are set on the pool well
        
    def test9_natural_product(self):
        
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_natural_product.xlsx'
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': 2037,  
            'end_plate': 2037, 
            'plate_size': '384',
            'library_type': 'natural_products' })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        self._load_xls_reagent_file(filename,library_item, 352, 384 )


    def _load_xls_reagent_file(self, filename,library_item, expected_in, expected_count):
        logger.debug(str(('==== test_load_xls_reagent_file =====')))
        
        start_plate = library_item['start_plate']
        end_plate = library_item['end_plate']

        logger.debug(str(('==== start: _load_xls_reagent_file =====')))

        resource_name = 'well'
        resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        
        data_for_get={}
        data_for_get.setdefault('limit', 0)
        data_for_get.setdefault('includes', ['*'])
        data_for_get.setdefault('HTTP_ACCEPT', 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' )
        xls_serializer = XLSSerializer()
        with open(filename) as input_file:
            
            ### NOTE: submit the data as an object, because the test framework
            ### will convert it to the target format.
            input_data = self.serializer.from_xlsx(input_file.read())
            input_data = input_data['objects']
        
            self.assertEqual(len(input_data), expected_in, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_in)))
            
            logger.debug(str(('======Submitting patch...', filename, resource_uri)))
        
            resp = self.api_client.put(
                resource_uri, format='xlsx', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [200, 204], 
                str((resp.status_code, xls_serializer.from_xlsx(resp.content))))
        
        logger.debug(str(('check patched data for',resource_name,
            'execute get on:',resource_uri)))

        resource_name = 'reagent'
        reagent_resource_uri = '/'.join([BASE_URI_DB,'library', library_item['short_name'],resource_name])
        resp = self.api_client.get(
            reagent_resource_uri, format='xlsx', authentication=self.get_credentials(), 
            data=data_for_get)

        logger.debug(str(('--------resp to get:', resp.status_code)))
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        new_obj = self.deserialize(resp)
        returned_data = new_obj['objects']
        self.assertEqual(len(returned_data), expected_count, 
            str(('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count )))

        # 1. test well keys
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        logger.debug(str(('=== well fields', fields.keys())))
        
        substance_ids = set()
        for inputobj in input_data:
            # first just search for the well
            search = { 'well_id': lims_utils.well_id(inputobj['plate_number'],inputobj['well_name']) }
            result, outputobj = find_obj_in_list(search,returned_data)
            self.assertTrue(result, 
                str(('not found', search,outputobj,'=== objects returned ===', 
                      returned_data )) ) 
            logger.debug(str(('found', search)))
            # second look at the found item
            
            expected_data = { key: inputobj[key] for key in fields.keys() 
                if key in inputobj }
            result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
            logger.debug(str((result, msgs)))
            self.assertTrue(result, str((msgs, expected_data, outputobj)))
            
            substance_id = outputobj['substance_id']
            self.assertTrue(substance_id not in substance_ids, 
                str(('substance_id not unique', substance_id, substance_ids)))
            substance_ids.add(substance_id)
                    

class ScreenResource(DBMetaHashResourceBootstrap):
        
    def setUp(self):
        logger.debug('============== ScreenResource setup ============')
        super(ScreenResource, self).setUp()
        logger.debug('============== ScreenResource setup: begin ============')
        self.db_resource_uri = BASE_URI + '/metahash'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        
        testApiClient = TestApiClient(serializer=reports.serializers.LimsSerializer) 

        filename = os.path.join(self.db_directory,'metahash_fields_screen.csv')
        self._patch_test(
            'metahash', filename, id_keys_to_check=['scope','key'], 
            data_for_get={ 'scope':'fields.screen'})

        logger.debug('============== ScreenResource setup: done ============')        

    def test1_create_screen(self):
        logger.debug(str(('==== test1_create_screen =====')))
        
        resource_uri = BASE_URI_DB + '/screen'
        
        screen_item = ScreenFactory.attributes()
        
        logger.info(str(('screen_item',screen_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], 
            str((resp.status_code, self.deserialize(resp))))
        
        screen_item = ScreenFactory.attributes()
        
        logger.info(str(('item', screen_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))
        
        data_for_get={ 'limit': 0, 'includes': ['*'] }        

        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.debug(str(('--------resp to get:', resp.status_code)))
        # self.assertValidJSONResponse(resp)
        self.assertTrue(resp.status_code in [200], str((resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(screen_item, new_obj['objects'])
        self.assertTrue(
            result, str(('bootstrap item not found', obj, 
                         screen_item, new_obj['objects'])))
        self.assertTrue('facility_id' in obj, 'the facility_id was not created')
        logger.debug(str(('item found', obj)))

    # FIXME: DISABLE this test - does not use Factory
    def _test0_create_screen(self):
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
            logger.info(str(('post item', item)))         
            resp = self.api_client.post(
                resource_uri, format='json', data=item, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], 
                str((resp.status_code, self.deserialize(resp))))
            logger.info(str(('resp to create', self.deserialize(resp))))
            
        data_for_get={ 'limit': 0, 'includes': ['*'] }        

        logger.debug('created items, now get them')
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.info(str(('--------resp to get:', resp.status_code)))
        # self.assertValidJSONResponse(resp)
        self.assertTrue(resp.status_code in [200], str((resp.status_code)))
        new_obj = self.deserialize(resp)
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
        
class ScreensaverUserResource(DBMetaHashResourceBootstrap):
        
    def setUp(self):
        logger.debug('============== ScreensaverUserResource setup ============')
        super(ScreensaverUserResource, self).setUp()
        logger.debug('============== ScreensaverUserResource setup: begin ============')
        self.db_resource_uri = BASE_URI + '/screensaveruser'
        self.db_directory = os.path.join(APP_ROOT_DIR, 'db/static/api_init')
        self.reports_directory = os.path.join(APP_ROOT_DIR, 'reports/static/api_init')
        
        testApiClient = TestApiClient(serializer=reports.serializers.LimsSerializer) 

        filename = os.path.join(self.reports_directory,'metahash_fields_user.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.user' })
        
        filename = os.path.join(self.directory,'metahash_fields_usergroup.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.usergroup'})
        
        filename = os.path.join(self.directory,'metahash_fields_permission.csv')
        self._patch_test('metahash', filename, data_for_get={ 'scope':'fields.permission'})

        filename = os.path.join(self.db_directory,'metahash_fields_screensaveruser.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.screensaveruser'},
            id_keys_to_check=['key','scope'])

        logger.debug('============== ScreensaverUserResource setup: done ============')        

    def test0_create_user(self):
        logger.debug(str(('==== test_create_user =====')))
        
        # the simplest of tests, create some simple users
        self.bootstrap_items = { 'objects': [   
            {
                'username': 'st1',
                'ecommons_id': 'st1',
                'first_name': 'Sally',
                'last_name': 'Tester', 
                'email': 'sally.tester@limstest.com',  
                'harvard_id': '123456',
                'harvard_id_expiration_date': '2015-05-01'  
            },
            {
                'username': 'jt1',
                'first_name': 'Joe',
                'last_name': 'Tester',    
                'email': 'joe.tester@limstest.com',    
                'harvard_id': 'A4321',
                'harvard_id_expiration_date': '2015-05-02'  
            },
            {
                'username': 'bt1',
                'first_name': 'Bad',
                'last_name': 'TestsALot',    
                'email': 'bad.tester@slimstest.com',    
                'harvard_id': '332122',
                'harvard_id_expiration_date': '2018-05-01'  
            },
        ]}
        resource_uri = BASE_URI_DB + '/screensaveruser'

        try:       
            resp = self.api_client.put(resource_uri, 
                format='json', data=self.bootstrap_items, authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [200,201,202], 
                str((resp.status_code, self.serialize(resp))))
            #                 self.assertHttpCreated(resp)
        except Exception, e:
            logger.error(str(('on creating', self.bootstrap_items, '==ex==', e)))
            raise

        logger.debug('created items, now get them')
        resp = self.api_client.get(resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        self.assertEqual(len(new_obj['objects']), 3, str((new_obj)))
        
        for i,item in enumerate(self.bootstrap_items['objects']):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(
                result, str(('bootstrap item not found', item, new_obj['objects'])))
            logger.debug(str(('item found', obj)))

        logger.debug(str(('==== test_create_user done =====')))

    def test2_patch_usergroups(self):
        
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
            self.assertTrue(resp.status_code in [200,201,202], 
                str((resp.status_code, self.serialize(resp))))

            resp = self.api_client.get(resource_uri, format='json', 
                authentication=self.get_credentials(), data={ 'limit': 999 })
            logger.debug(str(('--------resp to get:', resp.status_code)))
            new_obj = self.deserialize(resp)
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
            self.assertEqual(len(new_obj['objects']), len(group_patch['objects']), str((new_obj)))
            
            for i,item in enumerate(group_patch['objects']):
                result, obj = find_obj_in_list(item, new_obj['objects'])
                self.assertTrue(
                    result, str(('bootstrap item not found', item, new_obj['objects'])))
                logger.info(str(('item found', obj)))        
        except Exception, e:
            logger.error(str(('on group_patch', group_patch, '==ex==', e)))
            raise

        userpatch = { 'objects': [   
            {
                'username': 'st1',
                'usergroups': ['usergroup1']
            },
            {
                'username': 'jt1',
                'usergroups': ['usergroup2']
            },
            {
                'username': 'bt1',
                'usergroups': ['usergroup3']
            },
        ]};
        try:       
            resource_uri = BASE_URI_DB + '/screensaveruser'
            resp = self.api_client.put(resource_uri, 
                format='json', data=userpatch, authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [200,201,202], 
                str((resp.status_code, self.serialize(resp))))

            data_for_get={ 'limit': 0 }        
            data_for_get.setdefault('includes', ['*'])
            resp = self.api_client.get(resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get )
            logger.debug(str(('--------resp to get:', resp.status_code)))
            new_obj = self.deserialize(resp)
            self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
            self.assertEqual(len(new_obj['objects']), 3, str((new_obj)))
            
            for i,item in enumerate(userpatch['objects']):
                result, obj = find_obj_in_list(item, new_obj['objects'])
                self.assertTrue(
                    result, str(('bootstrap item not found', item, new_obj['objects'])))
                logger.debug(str(('item found', obj)))        
        
        except Exception, e:
            logger.error(str(('on userpatch ', userpatch, '==ex==', e)))
            raise
