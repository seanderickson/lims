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
from db.models import Reagent, Substance, Library
from south.migration.base import Migrations
from south.migration.migrators import FakeMigrator
from south.migration import get_migrator
from db.test.factories import LibraryFactory
from django.utils.timezone import now
from reports.dump_obj import dumpObj


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


# FIXME: reinstate this - following the version in reports
# class TestApiInit(MetaHashResourceBootstrap,ResourceTestCase):
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
        self._patch_test('vocabularies', filename, id_keys_to_check=['key'])

        
        logger.debug('------------- Done: DBMetaHashResourceBootstrap _bootstrap_init_files -----------------')        
        
          

class LibraryResource(DBMetaHashResourceBootstrap,ResourceTestCase):

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

#         self.testApiClient = TestApiClient(serializer=reports.serializers.SmallMoleculeSerializer()) 

        logger.debug('============== LibraryResource setup: done ============')
        

    def test1_create_library(self):
        logger.debug(str(('==== test1_create_library =====')))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 'start_plate': '1534', 'end_plate': '1534', 'plate_size': '384' })
        
        self.library1 = library_item # store for later use
        
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
         
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
        
        result, obj = find_obj_in_list(library_item, new_obj['objects'])
        self.assertTrue(
            result, str(('library_item', obj, 
                         library_item, new_obj['objects'])))
        logger.debug(str(('item found', obj)))
        
        # now find the library wells
        
        resource_uri = '/'.join([BASE_URI_DB,'library',self.library1['short_name'],'well'])
        logger.info(str(('GET', resource_uri)))
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 384, str((new_obj)))
        
        # now examine the wells created
        index = 0
        platesize = 384
        plate = 1534
        from db.support import lims_utils
        for well in new_obj['objects']:
            logger.info(str(('testing well', well)))
            well_name = lims_utils.well_name_from_index(index, platesize)
            well_id = lims_utils.well_id(plate,well_name)
            self.assertEqual(well_id, well['well_id'], 
                str(('not equal',index, well_id, well['well_id'])))
            index += 1
            
        logger.debug(str(('==== done: test1_create_library =====')))

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
                
                
        # TODO: test regex and number: min/max
        logger.debug(str(('==== done: test4_create_library_invalids =====')))

    def test6_load_small_molecule_file(self):
        logger.debug(str(('==== test6_load_small_molecule_file =====')))
        
        filename = os.path.join(self.db_directory,'metahash_fields_smallmoleculereagent.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.smallmoleculereagent'})

        filename = os.path.join(self.db_directory,'metahash_fields_naturalproductreagent.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.naturalproductreagent'})

        filename = os.path.join(self.db_directory,'metahash_fields_silencingreagent.csv')
        self._patch_test(
            'metahash', filename, data_for_get={ 'scope':'fields.silencingreagent'})

        logger.debug(str(('==== start: test6_load_small_molecule_file =====')))
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 'start_plate': '1534', 'end_plate': '1534', 'plate_size': '384' })
        
        self.library1 = library_item # store for later use
        
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))

        resource_name = 'reagent'
        resource_uri = '/'.join([BASE_URI_DB,'library', self.library1['short_name'],resource_name])
        wells_uri = resource_uri
        
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_data_small_molecule.sdf'

        data_for_get={}
        data_for_get.setdefault('limit', 999)
        data_for_get.setdefault('HTTP_ACCEPT', 'chemical/x-mdl-sdfile' )

        with open(filename) as input_file:
            
            input_data = self.serializer.from_sdf(input_file.read())
            input_data = input_data['objects']
            
            # "convert" the well data from SS1 keys/values
            ## Discuss: the "old" SS1 format has different field naming conventions, 
            # so we'll need to have a built in functionality to map them ##
            from db.support.data_converter import convert_well_data
            input_data = [convert_well_data(x) for x in input_data]
        
            expected_count = 8
            self.assertEqual(len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            logger.info(str(('======Submitting patch...', filename, resource_uri)))
        
            resp = self.api_client.put(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [200, 204], self.deserialize(resp))
        
        resource_uri = BASE_URI_DB + '/reagent'
        logger.debug(str(('check patched data for',resource_name,
            'execute get on:',resource_uri)))
 
        resp = self.api_client.get(
            resource_uri, format='sdf', authentication=self.get_credentials(), 
            data=data_for_get)
 
        logger.debug(str(('--------resp to get:', resp.status_code, self.deserialize(resp))))
        
        resource_uri = BASE_URI_DB + '/well'
        logger.debug(str(('check patched data for',resource_name,
            'execute get on:',resource_uri)))

        resp = self.api_client.get(
            resource_uri, format='sdf', authentication=self.get_credentials(), 
            data=data_for_get)

        logger.debug(str(('--------resp to get:', resp.status_code)))
        self.assertTrue(resp.status_code in [200], str((resp.status_code, resp)))
        new_obj = self.deserialize(resp)
        returned_data = new_obj['objects']
        self.assertEqual(len(returned_data), 384, 
            str(('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count,
                'returned_data',returned_data)))
        
        # test returned data:
        # NOTE: not all keys will be returned, since the clean_data file contains unused 
        # keys.  Rather, use the schema to get the well, reagent, and 
        # small_molecule_reagent fields and check only them
        
        # 1. test well keys
        specific_schema = self.get_from_server(wells_uri + '/schema')
#         well_resource = self.get_resource_from_server('well')
        fields = specific_schema['fields']
        logger.debug(str(('=== well fields', fields.keys())))
        
        excludes=['resource_uri','molecular_weight','molecular_mass'] # TODO: fixme mol wt
        for inputobj in input_data:
            well_data = { key: inputobj[key] for key in fields.keys() if key in inputobj }
            result, outputobj = find_obj_in_list(
                well_data,returned_data, excludes=excludes )
            self.assertTrue(
                result, 
                str(('not found', outputobj,'=== objects returned ===', 
                     returned_data )) ) 
            if result:
                logger.info(str(('found', inputobj, 'outputobj', outputobj)))
        logger.debug(str(('==== done: test6_load_small_molecule_file =====')))

          
class ReagentResource(DBMetaHashResourceBootstrap):
            
    def setUp(self):
        logger.debug('============== ReagentResource setup ============')
        
        super(ReagentResource, self).setUp()

        logger.debug('============== ReagentResource setup: begin ============')
        
        
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

        library_item = LibraryFactory.attributes()
        library_item.update({ 'start_plate': '1534', 'end_plate': '1534', 'plate_size': '384' })
        
        self.library1 = library_item # store for later use
        
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(str((library_item)))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], str((resp.status_code, resp)))

        
        logger.debug('============== done: ReagentResource setup ============')
        
    def test1_create_reagent(self):
        logger.debug(str(('==== test1_create_reagent =====')))
        
        resource_uri = '/'.join([BASE_URI_DB,'library', self.library1['short_name'], 'reagent'])
        reagent_item1 = ReagentFactory.attributes()
        reagent_item1['well_id'] = '1534:A01'
        logger.debug(str((reagent_item1)))
        resp = self.api_client.post(
            resource_uri, format='json', data=reagent_item1, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], 
            str((resp.status_code, self.deserialize(resp))))
         
        # create a second
        reagent_item2 = ReagentFactory.attributes()
        reagent_item2['well_id'] = '1534:A02'
         
        logger.debug(str((reagent_item1)))
        resp = self.api_client.post(
            resource_uri, format='json', data=reagent_item2, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], 
            str((resp.status_code, self.deserialize(resp))))
         
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data={ 'limit': 999 })
        logger.debug(str(('--------resp to get:', resp.status_code)))
        new_obj = self.deserialize(resp)
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(new_obj['objects']), 2, str((new_obj)))
         
        result, obj = find_obj_in_list(reagent_item1, new_obj['objects'])
        self.assertTrue(
            result, str(('item not found', obj, 
                         reagent_item1, new_obj['objects'])))
        logger.debug(str(('item found', obj)))
 
        result, obj = find_obj_in_list(reagent_item2, new_obj['objects'])
        self.assertTrue(
            result, str(('bootstrap item2 not found', obj, 
                         reagent_item2, new_obj['objects'])))
        logger.debug(str(('item2 found', obj)))
 
        logger.debug(str(('==== done: test1_create_reagent =====')))




class ScreenResource(DBMetaHashResourceBootstrap,ResourceTestCase):
        
    def setUp(self):
        logger.debug('============== ScreenResource setup ============')
        super(ScreenResource, self).setUp()
#         super(ScreenResource, self)._setUp()
        # load the bootstrap files, which will load the metahash fields, 
        # and the resource definitions
#         super(ScreenResource, self)._bootstrap_init_files()
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
        
