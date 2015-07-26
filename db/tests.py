import cStringIO
import json
import logging
import os

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.http.response import StreamingHttpResponse
from django.test import TestCase
from django.utils.timezone import now
from tastypie.test import ResourceTestCase, TestApiClient

from db.models import Reagent, Substance, Library
import db.models
from db.support import lims_utils
from db.test.factories import LibraryFactory, ScreenFactory
from reports.dump_obj import dumpObj
from reports.serializers import CSVSerializer, XLSSerializer, LimsSerializer
from reports.sqlalchemy_resource import SqlAlchemyResource
from reports.tests import IResourceTestCase
from reports.tests import assert_obj1_to_obj2, find_all_obj_in_list, find_obj_in_list
import reports.tests
import reports.utils.log_utils


_ = reports.utils.log_utils.LogMessageFormatter
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
        logger.info(_('open input_actions file', input_actions_file))
        
        self._run_api_init_actions(input_actions_file)


def setUpModule():

    logger.info(_('=== setup Module'))
    
    # FIXME: running the bootstrap method as a test suite setup:
    # TODO: create a local TestRunner,TestSuite,
    # so that this can be run before the suite
    testContext = DBResourceTestCase(methodName='_bootstrap_init_files')
    testContext.setUp()
    testContext._bootstrap_init_files()

    logger.info(_('=== setup Module done'))

def tearDownModule():

    logger.info(_('=== teardown Module'))
    
    # FIXME: close the sqlalchemy bridge connection on tearDown
    # - the right solution probably requires a custom TestRunner
    # This does not work:
    # bridge.get_engine().connect().close()
    # bridge.get_engine().dispose()
    # bridge = None
    logger.info(_('=== teardown Module done'))

class LibraryResource(DBResourceTestCase):

    def setUp(self):
        logger.debug('============== LibraryResource setup ============')
        
        super(LibraryResource, self).setUp()

    def test1_create_library(self):
        
        logger.debug(_('==== test1_create_library ====='))
        
        resource_uri = BASE_URI_DB + '/library'
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': '1534', 'end_plate': '1534', 'plate_size': '384' })
        
        logger.debug(_(library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))

        # create a second library
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': '1535', 'end_plate': '1537', 'plate_size': '384' })
         
        logger.debug(_('item', library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))
         
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.info(_('--------resp to get:', resp.status_code))
        self.assertTrue(resp.status_code in [200], (resp.status_code))
        new_obj = self.deserialize(resp)
        logger.info(_('resp', new_obj))
        self.assertEqual(len(new_obj['objects']), 2, (new_obj))
        
        result, obj = find_obj_in_list(library_item, new_obj['objects'])
        self.assertTrue(
            result, _('library_item', obj, 
                         library_item, new_obj['objects']))
        logger.debug(_('item found', obj))
        
        # now find the library wells
        resource_uri = '/'.join([
            BASE_URI_DB,'library',library_item['short_name'],'well'])
        data_for_get={ 'limit': 0 }        
        data_for_get.setdefault('includes', ['*'])
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.debug(_('--------resp to get:', 
            resp.status_code, resp['Content-Type'], resp))
        self.assertTrue(resp.status_code in [200], (resp.status_code))
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
                
        logger.debug(_('==== done: test1_create_library ====='))

    def test4_create_library_invalids(self):
        '''
        Test the schema "required" validations 
        
        - todo: this can be a template for other endpoints
        ''' 
        logger.debug(_('==== test4_create_library_invalids ====='))

        library_resource = self.get_resource_from_server('library')
        fields = library_resource['schema']['fields']

        # make sure the default works
        resource_uri = BASE_URI_DB + '/library'
        library_item = LibraryFactory.attributes()
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))
        
        for key,field in fields.items():
            if field.get('required', False):
                logger.debug(_('testing required field', key, field))
                
                library_item = LibraryFactory.attributes()
                library_item[key] = None
                resp = self.api_client.post(
                    resource_uri, format='json', data=library_item, 
                    authentication=self.get_credentials())
                self.assertTrue(resp.status_code in [400], (resp.status_code, resp))
                logger.debug(_('response.content.library message', 
                                 getattr(resp, 'content')))
        
                obj = json.loads(getattr(resp, 'content'))
                logger.debug(_('==========content', obj))
                self.assertTrue(key in obj['library'], 
                    _('error content should have an entry for the faulty key:',key, obj))

        # test 2                
        library_item = LibraryFactory.attributes()
        library_item['library_name'] = 'invalid & name'
        
        logger.debug(_('try to create an invalid library name:', library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        
        self.assertTrue(resp.status_code in [400], (resp.status_code, resp))
        
        logger.debug(_('response.content.library message', 
                         getattr(resp, 'content')))
        
        obj = json.loads(getattr(resp, 'content'))
        self.assertTrue('library_name' in obj['library'], 
            _('content should have an entry for the faulty "name"', obj))

        # test 3
        library_item = LibraryFactory.attributes()
        library_item['library_type'] = 'invalid_type'
        
        logger.debug(_('try to create an invalid library_type:', library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [400], (resp.status_code, resp))
        logger.debug(_('response.content.library message', 
                         getattr(resp, 'content')))
        obj = json.loads(getattr(resp, 'content'))
        self.assertTrue('library_type' in obj['library'], 
            _('content should have an entry for the faulty "name"', obj))

        # TODO: test regex and number: min/max, vocabularies

        logger.debug(_('==== done: test4_create_library_invalids ====='))
        
    def test6_load_small_molecule_file(self):
        
        logger.debug(_('==== test6_load_small_molecule_file ====='))
        
        # setup library
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': '1536', 'end_plate': '1536', 'plate_size': '384' })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(_(library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))

        logger.debug(_('==== start: test6_load_small_molecule_file ====='))

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
            expected_count = 8
            self.assertEqual(len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            logger.debug(_('======Submitting patch...', filename, resource_uri))
            resp = self.api_client.put(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [200, 204], 
                _(resp.status_code, self.deserialize(resp)))
        
        logger.debug(_('check patched data for',resource_name,
            'execute get on:',resource_uri))

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
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        new_obj = self.deserialize(resp)
        returned_data = new_obj['objects']
        expected_count = 384
        self.assertEqual(len(returned_data), expected_count, 
            ('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count,
                'returned_data',returned_data))
        
        # test 1: test well keys
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        logger.debug(_('=== well fields', fields.keys()))
        
        substance_ids = set()
        for inputobj in input_data:
            # first just search for the well
            search = { 'well_id': 
                lims_utils.well_id(
                    inputobj['plate_number'],inputobj['well_name']) }
            result, outputobj = find_obj_in_list(
                search,returned_data) #, excludes=excludes )
            self.assertTrue(result, 
                _('not found', search,outputobj,'=== objects returned ===', 
                      returned_data ) ) 
            logger.debug(_('found', search))
            # second look at the found item
            expected_data = { 
                key: inputobj[key] for key in fields.keys() if key in inputobj }
            result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
            logger.debug(_(result, msgs))
            self.assertTrue(result, (msgs, 'input', expected_data, 'output', outputobj))

            substance_id = outputobj['substance_id']
            self.assertTrue(substance_id not in substance_ids, 
                _('substance_id not unique', substance_id, substance_ids))
            substance_ids.add(substance_id)
                    
        logger.debug(_('==== done: test6_load_small_molecule_file ====='))
    
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
            self.assertEqual(len(input_data), expected_count, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_count,
                    'input_data',input_data)))
            
            logger.debug(_('======Submitting patch...', filename, resource_uri))
            resp = self.api_client.put(
                resource_uri, format='sdf', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [200, 204], 
                _(resp.status_code, self.deserialize(resp)))
        
            logger.debug(_('check updated/patched data for',resource_name,
                'execute get on:',resource_uri))
            resource_name = 'well'
            resource_uri = '/'.join([BASE_URI_DB,'library', library_item['short_name'],resource_name])
            resp = self.api_client.get(
                reagent_resource_uri, format='sdf', authentication=self.get_credentials(), 
                data=data_for_get)
            self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
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
                search = { 
                    'well_name': update_well['well_name'], 
                    'plate_number': update_well['plate_number']}
                result, outputobj = find_obj_in_list(
                    search,returned_data) #, excludes=excludes )
                self.assertTrue(result, 
                    _('not found', search,outputobj,'=== objects returned ===', 
                          returned_data ) ) 
                logger.debug(_('found', search))
                result, msgs = assert_obj1_to_obj2(update_well, outputobj)
                logger.debug(_(result, msgs))
                self.assertTrue(result, (msgs, update_well, outputobj))

            # 2. check the apilogs - library
            resource_uri = BASE_REPORTS_URI + '/apilog' #?ref_resource_name=record'
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 
                    'limit': 0, 
                    'ref_resource_name': 'library', 
                    'key': library_item['short_name'] })
            self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
            new_obj = self.deserialize(resp)
            logger.debug(_('===library apilogs:', json.dumps(new_obj)))
            expected_count = 3 # create, post, update
            self.assertEqual( len(new_obj['objects']), expected_count , 
                str((len(new_obj['objects']), expected_count, new_obj)))

            # 2. check the apilogs - wells
            resource_uri = BASE_REPORTS_URI + '/apilog' 
            resp = self.api_client.get(
                resource_uri, format='json', 
                authentication=self.get_credentials(), 
                data={ 'limit': 0, 'ref_resource_name': 'well' })
            self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
            new_obj = self.deserialize(resp)
            logger.debug(_('===apilogs:', json.dumps(new_obj)))
            
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
        logger.debug(_(library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))
        
        self._load_xls_reagent_file(filename,library_item, 5,384 )

    #     def test7a_sirna_validations(self):
    #         # TODO: test validations
    #         pass
    #     

    def test8_sirna_duplex(self):
        
        filename = ( APP_ROOT_DIR 
            + '/db/static/test_data/libraries/clean_rnai_duplex_50440_50443.xlsx')
        
        # create the duplex library
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': 50440,  
            'end_plate': 50443, 
            'plate_size': '384',
            'screen_type': 'rnai' })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(_(library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))
        
        self._load_xls_reagent_file(filename,library_item, 532, 1536 )
        filename = APP_ROOT_DIR + '/db/static/test_data/libraries/clean_rnai_pool.xlsx'
        
        # create the pool library
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': 50439,  
            'end_plate': 50439, 
            'plate_size': '384',
            'screen_type': 'rnai',
            'is_pool': True })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(_(library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))
        
        self._load_xls_reagent_file(filename,library_item, 133, 384 )
        
        ## TODO: test the duplex wells are set on the pool well
        
    def test9_natural_product(self):
        
        filename = ( APP_ROOT_DIR 
            + '/db/static/test_data/libraries/clean_data_natural_product.xlsx')
        
        library_item = LibraryFactory.attributes()
        library_item.update({ 
            'start_plate': 2037,  
            'end_plate': 2037, 
            'plate_size': '384',
            'library_type': 'natural_products' })
        resource_uri = BASE_URI_DB + '/library'
        logger.debug(_(library_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=library_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))
        
        self._load_xls_reagent_file(filename,library_item, 352, 384 )

    def _load_xls_reagent_file(self,
            filename,library_item, expected_in, expected_count):

        logger.debug(_('==== test_load_xls_reagent_file ====='))
        
        start_plate = library_item['start_plate']
        end_plate = library_item['end_plate']

        logger.debug(_('==== start: _load_xls_reagent_file ====='))

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
        
            self.assertEqual(len(input_data), expected_in, 
                str(('initial serialization of ',filename,'found',
                    len(input_data), 'expected',expected_in)))
            
            logger.debug(_('======Submitting patch...', filename, resource_uri))
        
            resp = self.api_client.put(
                resource_uri, format='xlsx', data=input_data, 
                authentication=self.get_credentials(), **data_for_get )
            self.assertTrue(resp.status_code in [200, 204], 
                _(resp.status_code, xls_serializer.from_xlsx(resp.content)))
        
        logger.debug(_('check patched data for',resource_name,
            'execute get on:',resource_uri))

        resource_name = 'reagent'
        reagent_resource_uri = '/'.join([
            BASE_URI_DB,'library', library_item['short_name'],resource_name])
        resp = self.api_client.get(
            reagent_resource_uri, format='xlsx', authentication=self.get_credentials(), 
            data=data_for_get)

        logger.debug(_('--------resp to get:', resp.status_code))
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        new_obj = self.deserialize(resp)
        returned_data = new_obj['objects']
        self.assertEqual(len(returned_data), expected_count, 
            str(('returned_data of ',filename,'found',
                len(returned_data), 'expected',expected_count )))

        # 1. test well keys
        specific_schema = self.get_from_server(resource_uri + '/schema')
        fields = specific_schema['fields']
        logger.debug(_('=== well fields', fields.keys()))
        
        substance_ids = set()
        for inputobj in input_data:
            # first just search for the well
            search = { 
                'well_id': lims_utils.well_id(
                    inputobj['plate_number'],inputobj['well_name']) }
            result, outputobj = find_obj_in_list(search,returned_data)
            self.assertTrue(result, 
                _('not found', search,outputobj,'=== objects returned ===', 
                      returned_data ) ) 
            logger.debug(_('found', search))
            # second look at the found item
            
            expected_data = { key: inputobj[key] for key in fields.keys() 
                if key in inputobj }
            result, msgs = assert_obj1_to_obj2(expected_data, outputobj)
            logger.debug(_(result, msgs))
            self.assertTrue(result, (msgs, expected_data, outputobj))
            
            substance_id = outputobj['substance_id']
            self.assertTrue(substance_id not in substance_ids, 
                _('substance_id not unique', substance_id, substance_ids))
            substance_ids.add(substance_id)
                    

class ScreenResource(DBResourceTestCase):
        
    def setUp(self):
        logger.debug('============== ScreenResource setup ============')
        super(ScreenResource, self).setUp()

    def test1_create_screen(self):
        logger.debug(_('==== test1_create_screen ====='))
        
        resource_uri = BASE_URI_DB + '/screen'
        
        # test 1: create a screen
        screen_item = ScreenFactory.attributes()
        logger.info(_('screen_item',screen_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], 
            _(resp.status_code, self.deserialize(resp)))
        
        # create another screen
        screen_item = ScreenFactory.attributes()
        logger.info(_('item', screen_item))
        resp = self.api_client.post(
            resource_uri, format='json', data=screen_item, 
            authentication=self.get_credentials())
        self.assertTrue(resp.status_code in [201], (resp.status_code, resp))

        # find the screens        
        data_for_get={ 'limit': 0, 'includes': ['*'] }        
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.debug(_('--------resp to get:', resp.status_code))
        self.assertTrue(resp.status_code in [200], (resp.status_code))
        new_obj = self.deserialize(resp)
        self.assertEqual(len(new_obj['objects']), 2, (new_obj))
        
        result, obj = find_obj_in_list(screen_item, new_obj['objects'])
        self.assertTrue(
            result, _('bootstrap item not found', obj, 
                         screen_item, new_obj['objects']))
        self.assertTrue('facility_id' in obj, 'the facility_id was not created')
        logger.debug(_('item found', obj))

    def _test0_create_screen(self):
        logger.debug(_('==== test_create_screen ====='))
        
        # FIXME: should use the Factory
        # the simplest of tests, create some simple screens
        simple_screen_input = [   
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
        logger.debug(_('--posting to:', resource_uri))
        for i,item in enumerate(simple_screen_input):
            logger.info(_('post item', item))         
            resp = self.api_client.post(
                resource_uri, format='json', data=item, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [201], 
                _(resp.status_code, self.deserialize(resp)))
            logger.info(_('resp to create', self.deserialize(resp)))
            
        data_for_get={ 'limit': 0, 'includes': ['*'] }        

        logger.debug('created items, now get them')
        resp = self.api_client.get(
            resource_uri, format='json', authentication=self.get_credentials(), 
            data=data_for_get)
        logger.info(_('--------resp to get:', resp.status_code))
        self.assertTrue(resp.status_code in [200], (resp.status_code))
        new_obj = self.deserialize(resp)
        self.assertEqual(len(new_obj['objects']), 2, new_obj)
        
        for i,item in enumerate(simple_screen_input):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(
                result, _('bootstrap item not found', item, new_obj['objects']))
            self.assertTrue(
                'facility_id' in obj, 'the facility_id was not created')
            self.assertEqual(
                str(i+1), obj['facility_id'], 
                ('expected the facility_id returned to be incremented to ',
                 i+1, obj) )
            logger.debug(_('item found', obj))
            assert_obj1_to_obj2(item, obj)

        logger.debug(_('==== test_create_screen done ====='))
        

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

    def test0_create_user(self):
        logger.debug(_('==== test_create_user ====='))
        
        # the simplest of tests, create some simple users
        simple_user_input = { 'objects': [   
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
                format='json', data=simple_user_input, 
                authentication=self.get_credentials())
            self.assertTrue(resp.status_code in [200,201,202], 
                _(resp.status_code, self.serialize(resp)))
        except Exception, e:
            logger.exception(_('on creating', simple_user_input), e)
            raise

        logger.debug('created items, now get them')
        resp = self.api_client.get(resource_uri, format='json', 
            authentication=self.get_credentials(), data={ 'limit': 999 })
        new_obj = self.deserialize(resp)
        self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
        self.assertEqual(len(new_obj['objects']), 3, (new_obj))
        
        for i,item in enumerate(simple_user_input['objects']):
            result, obj = find_obj_in_list(item, new_obj['objects'])
            self.assertTrue(
                result, _('bootstrap item not found', item, new_obj['objects']))
            logger.debug(_('item found', obj))

        logger.debug(_('==== test_create_user done ====='))

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
                _(resp.status_code, self.serialize(resp)))

            resp = self.api_client.get(resource_uri, format='json', 
                authentication=self.get_credentials(), data={ 'limit': 999 })
            new_obj = self.deserialize(resp)
            self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
            self.assertEqual(len(new_obj['objects']), len(group_patch['objects']), new_obj)
            
            for i,item in enumerate(group_patch['objects']):
                result, obj = find_obj_in_list(item, new_obj['objects'])
                self.assertTrue(result, 
                    _('bootstrap item not found', item, new_obj['objects']))
                logger.info(_('item found', obj))        
        except Exception, e:
            logger.exception(_('on group_patch', group_patch), e)
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
                _(resp.status_code, self.serialize(resp)))

            data_for_get={ 'limit': 0 }        
            data_for_get.setdefault('includes', ['*'])
            resp = self.api_client.get(resource_uri, format='json', 
                authentication=self.get_credentials(), data=data_for_get )
            new_obj = self.deserialize(resp)
            self.assertTrue(resp.status_code in [200], (resp.status_code, resp))
            self.assertEqual(len(new_obj['objects']), 3, (new_obj))
            
            for i,item in enumerate(userpatch['objects']):
                result, obj = find_obj_in_list(item, new_obj['objects'])
                self.assertTrue(
                    result, _('bootstrap item not found', item, new_obj['objects']))
                logger.debug(_('item found', obj))        
        
        except Exception, e:
            logger.exception(_('on userpatch ', userpatch), e)
            raise
