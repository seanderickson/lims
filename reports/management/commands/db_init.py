import os
import logging
import requests
import sys
import os
import traceback

from subprocess import Popen, PIPE
from lims.api import CSVSerializer
import json
from django.core.management.base import BaseCommand, CommandError
from optparse import make_option

logger = logging.getLogger(__name__)

class ApiError(Exception):
    
    def __init__(self, url, action, result):
        err_msg = ''
        try:
            json_status = result.json()
            if json_status and 'error_message' in json_status:
                err_msg = json_status['error_message']
        except ValueError,e:
            logger.info('There is no json in the response')
            logger.debug(str(('-----raw response text-------', result.text)) )
        self.message = str((url,'action',action, result.reason, result.status_code, err_msg ))
        if(logger.isEnabledFor(logging.DEBUG)):
            self.message = str((url,'action',action, result.reason, result.status_code, result ))

    def __str__(self):
        return repr(self.message)

class Command(BaseCommand):
    ''' 
    Initialize the reports/lims database:
    - with metahash (field, resource) information
    - data load files 
    '''
    help = '''Initialize the reports/lims database:
            - with metahash (field, resource) information
            - data load files 
            '''
    option_list = BaseCommand.option_list + (
        make_option('-d', '--inputDirectory', action='store', dest='inputDirectory',
                    metavar='DIRECTORY', 
                    help='input file path (where the input actions files can be found)'),
        make_option('-f', '--input_actions_file', action='store', dest='inputActionsFile',
                    metavar='INPUT_ACTIONS', 
                    help='''input actions file.  
                            An input action is a [command,resource,file]; 
                            see the reports/static/api_init/api_init_actions.csv file''' ),
        make_option('-u', '--url', action='store', dest='url',
                    metavar='URL', 
                    help='api url, e.g.: "http://localhost:8000/reports/api/v1"'),
        make_option('-a', '--auth', action='store', dest='auth',
                    metavar='AUTH',
                    help='api authentication to use <user:password>'),
        )

    args = '<' + ','.join([str(x) for x in option_list]) + '>'

    def handle(self, *args, **options):
        logger.info(str(('db_init, args', args, 'options', options)))
        
        inputDir = options['inputDirectory']
        if not inputDir:
            self.print_help('db_init', 'require inputDirectory option')
            raise CommandError("Option `--inputDirectory=...` must be specified.")
            return
        
        inputActionsFile = options['inputActionsFile']
        if not inputActionsFile:
            self.print_help('db_init', 'require inputActionsFile option')
            raise CommandError("Option `--inputActionsFile=...` must be specified.")
            return
        url = options['url']
        if not url:
            self.print_help('db_init', 'require url option')
            raise CommandError("Option `--url=...` must be specified.")
            return
        auth = options['auth']        
        if not auth:
            self.print_help('db_init', 'require auth option')
            raise CommandError("Option `--auth=...` must be specified.")
            return
        
        authentication = None
        if(auth):
            authentication = tuple(auth.split(':'))
        headers =   { 
            'json': {'content-type': 'application/json'},
            'csv':  {'content-type': 'text/csv'},
            } 
                    
        serializer=CSVSerializer() 
    
    #    # Identify some files to treat specially during verification
    #    # TODO: make verification optional
    #    bootstrap_files = ['metahash_fields_initial.csv','metahash_fields_initial_patch.csv','metahash_fields_resource.csv','metahash_resource_data.csv']
    #    # spefically, what not to verify for these files (Note: the resource table has to be filled before uri's will be correct)
    #    excludes = ['resource_uri']
    
#        filename = os.path.join(inputDir,inputActionsFile)
        with open(inputActionsFile) as input_file:
            api_init_actions = serializer.from_csv(input_file.read(), root=None)
    
            for action in api_init_actions:
                
                print '\n++++=========== processing action', json.dumps(action)
                command = action['command'].lower() 
                resource = action['resource'].lower()
                resource_uri = url + '/' + resource
                
                if command == 'delete':
                    header = headers['json']
                    delete(resource_uri, header, authentication)
                
                else:
                    data_file = os.path.join(inputDir,action['file'])
                    extension = os.path.splitext(data_file)[1].strip('.')
                    header = headers[extension]
    
    #                search_excludes = excludes
    #                if filename in bootstrap_files:
    #                    search_excludes = []
                    logger.info(str(('+++++++++++processing file', data_file)))
    #                with open(filename) as data_file:
    #                    input_data = serializer.deserialize(data_file.read(), header['content-type'])
                        
                    if command == 'put':
                        put(data_file, resource_uri, header, authentication)
                            
    #                    # now see if we can get these objects back
    #                    resp = testApiClient.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 })
    #                    self.assertValidJSONResponse(resp)
    #                    new_obj = self.deserialize(resp)
    #                    result, msgs = find_all_obj_in_list(input_data['objects'], new_obj['objects'], excludes=search_excludes)
    #                    self.assertTrue(result, str(( command, 'input file', filename, msgs, new_obj['objects'])) )
                        
                    elif command == 'patch':
                        patch(data_file, resource_uri, header, authentication)
      
    #                    resp = testApiClient.patch(resource_uri, format='csv', data=input_data, authentication=self.get_credentials() )
    ##                            print 'action: ', json.dumps(action), 'response: ' , resp.status_code, resp
    #                    self.assertHttpAccepted(resp)
    #                    resp = testApiClient.get(resource_uri, format='json', authentication=self.get_credentials(), data={ 'limit': 999 } )
    ##                            print '----- response', resp
    #                    self.assertValidJSONResponse(resp)
    #                    new_obj = self.deserialize(resp)
    #                    with open(filename) as f2:
    #                        input_data2 = serializer.from_csv(f2.read())
    #                        result, msgs = find_all_obj_in_list(input_data2['objects'], new_obj['objects'], excludes=search_excludes)
    #                        self.assertTrue(result, str(( command, 'input file', filename, msgs )) )
                        
                    else:
                        raise Exception(
                            str((   'Only [ PUT, PATCH, DELETE ] are supported.',
                                    '- cannot POST multiple objects to tastypie; therefore the "post" command is invalid with the initialization scripts',
                                    'resource entry: ',json.dumps(action),
                                    )) )


def delete(obj_url, headers, authentication):
    try:
        r = requests.delete(obj_url, auth=authentication, headers=headers)
        if(r.status_code != 204):
            print "ERROR", r
            raise ApiError(obj_url,'DELETE',r)
        print 'DELETE: ', obj_url, ' ,response:', r.status_code
        logger.info(str(('DELETE', obj_url)))
    except Exception, e:
        logger.error(str(('exception recorded while contacting server', e)))
        raise e

def put(input_file, obj_url,headers,authentication):
    try:
        with open(input_file) as f:
            r = requests.put(obj_url, auth=authentication, headers=headers, data=f.read() )
            if(r.status_code != 200):
                raise ApiError(obj_url,'PUT',r)
            print 'PUT: ' , input_file, 'to ', obj_url,' ,response:', r.status_code, ', count: ',len(r.json()['objects'])
            if(logger.isEnabledFor(logging.DEBUG)):
                logger.debug('--- PUT objects:')
                try:
                    for obj in r.json()['objects']:
                        logger.debug(str((obj)))
                except ValueError,e:
                    logger.debug('----no json object to report')
                    logger.debug(str(('text response', r.text)))
    except Exception, e:
        extype, ex, tb = sys.exc_info()
        logger.warn(str(('throw', e, tb.tb_frame.f_code.co_filename, 'error line', tb.tb_lineno, extype, ex)))
        logger.error(str(('exception recorded while contacting server', e)))
        raise e
    
def patch(patch_file, obj_url,headers,authentication):
    try:
        print 'PATCH: ' , patch_file, 'to ', obj_url
        with open(patch_file) as f:
            r = requests.patch(obj_url, auth=authentication, headers=headers, data=f.read() )
            if(r.status_code not in [200,202,204]):
                raise ApiError(obj_url,'PATCH',r)
            print 'PATCH: ', patch_file, ', to: ',obj_url,' ,response:', r.status_code,', count: ',len(r.json()['objects'])
            if(logger.isEnabledFor(logging.DEBUG)):
                logger.debug('--- PATCHED objects:')
                try:
                    for obj in r.json()['objects']:
                        logger.debug(str((obj)))
                except ValueError,e:
                    logger.debug('----no json object to report')
                    logger.debug(str(('text response', r.text)))
    except Exception, e:
        logger.error(str(('exception recorded while contacting server', e)))
        raise e
    