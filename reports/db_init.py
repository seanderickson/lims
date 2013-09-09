import argparse
import api_importer
from api_importer import ApiError
import os
import logging
import requests

from subprocess import Popen, PIPE

logger = logging.getLogger(__name__)

def main(inputDir, url, auth):
    
    url_base = url
    
    authentication = None
    if(auth):
        authentication = tuple(auth.split(':'))
    headers = {'content-type': 'application/json'}
    input_file_dir = inputDir
        
    # delete the metahash records
    obj_url = url_base + '/metahash/'
    delete(obj_url, headers, authentication)

    # fyi, system way to do this, has been replaced above    
    #pipe = Popen(["curl", "-H 'Content-Type: application/json' --user sde4:3xc3l3nc -X DELETE  http://localhost:8000/reports/api/v1/metahash/",
    #             shell=True, stdin=PIPE).stdin
        
    # initial import
#    api_importer.main(input_file_dir + '/metahash_fields_initial.json', url_base + '/metahash/',  'POST', auth=':'.join(authentication))
#    obj_url = url_base + '/metahash/'
#    headers = {'content-type': 'application/json'}
#    input_file = input_file_dir + '/metahash_fields_initial.json'
#    # NOTE: cannot POST multiple objects to tastypie at this time.  PUT a list will delete also
#    put(input_file, obj_url, headers, authentication)
    
    
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_initial.csv'
    # NOTE: cannot POST multiple objects to tastypie at this time.  PUT a list will delete also
    put(input_file, obj_url, headers, authentication)
    
    # now that the field defs are created, including the json fields, we can populate these json field values
    
    # patch the metahash records
    
    # Note: one-at-a-time api_importer version shown:
    # NOTE: using only "key" here, although the real key is the composite key (scope, key).  this works here
    # because the table is only loaded with the "metahash:fields" scope in the beginning
    #    api_importer.main(input_file_dir + '/metahash_fields_initial_patch.json',url_base + '/metahash/',  'PATCH', auth=':'.join(authentication)) #,  id_key='key')
    
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_initial_patch.csv'
    patch(input_file, obj_url, headers, authentication)    
    
    # define resource fields
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_resource.csv'
    patch(input_file, obj_url, headers, authentication)
    
    # import resource defs
    obj_url = url_base + '/resource/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_resource_data.csv'
    put(input_file, obj_url, headers, authentication)
            
    # define vocabularies fields
#    api_importer.main(input_file_dir + '/metahash_fields_vocabularies.json',url_base + '/metahash/',  'POST', auth=':'.join(authentication))
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_vocabularies.csv'
    patch(input_file, obj_url, headers, authentication)
    
    # define screensaveruser
#    api_importer.main(input_file_dir + '/metahash_fields_screensaveruser.json',url_base + '/metahash/',  'POST', auth=':'.join(authentication))
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_screensaveruser.csv'
    patch(input_file, obj_url, headers, authentication)

    # define screen
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_screen.csv'
    patch(input_file, obj_url, headers, authentication)

    # define screensummary
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_screensummary.csv'
    patch(input_file, obj_url, headers, authentication)

    # define screenresult (non-dynamic fields)
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_screenresult.csv'
    patch(input_file, obj_url, headers, authentication)

    # define datacolumn
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_datacolumn.csv'
    patch(input_file, obj_url, headers, authentication)

    # define apilog
    obj_url = url_base + '/metahash/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/metahash_fields_apilog.csv'
    patch(input_file, obj_url, headers, authentication)

    # delete the vocabularies records
    obj_url = url_base + '/vocabularies/'
    headers = {'content-type': 'application/json'}
    delete(obj_url, headers, authentication)

    # import vocabularies
#    api_importer.main(input_file_dir + '/vocabularies_data.json',url_base + '/vocabularies/',  'POST', auth=':'.join(authentication))
    obj_url = url_base + '/vocabularies/'
    headers = {'content-type': 'text/csv'}
    input_file = input_file_dir + '/vocabularies_data.csv'
    put(input_file, obj_url, headers, authentication)

def delete(obj_url, headers, authentication):
    try:
        r = requests.delete(obj_url, auth=authentication, headers=headers)
        if(r.status_code != 204):
            raise ApiError(obj_url,'DELETE',r)
        print 'DELETE: ', obj_url
        logger.info(str(('DELETE', obj_url)))
    except Exception, e:
        logger.error(str(('exception recorded while contacting server', e)))
        raise e

def put(input_file, obj_url,headers,authentication):
    try:
        with open(input_file) as f:
            r = requests.put(obj_url, auth=authentication, headers=headers, data=f.read() )
            if(r.status_code != 202):
                raise ApiError(obj_url,'PUT',r)
            print 'PUT: ' , input_file, 'to ', obj_url, ', count: ',len(r.json()['objects'])
            if(logger.isEnabledFor(logging.DEBUG)):
                logger.debug('--- PUT objects:')
                try:
                    for obj in r.json()['objects']:
                        logger.debug(str((obj)))
                except ValueError,e:
                    logger.debug('----no json object to report')
                    logger.debug(str(('text response', r.text)))
    except Exception, e:
        logger.error(str(('exception recorded while contacting server', e)))
        raise e
    
def patch(patch_file, obj_url,headers,authentication):
    try:
        with open(patch_file) as f:
            r = requests.patch(obj_url, auth=authentication, headers=headers, data=f.read() )
            if(r.status_code != 202):
                raise ApiError(obj_url,'PATCH',r)
            print 'PATCH: ' , patch_file, 'to ', obj_url, ', count: ',len(r.json()['objects'])
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
    
    
parser = argparse.ArgumentParser(description='Import file')

parser.add_argument('-v', '--verbose', dest='verbose', action='count',
                help="Increase verbosity (specify multiple times for more)")    
parser.add_argument('-d', action='store', dest='inputDirectory',
                    metavar='DIRECTORY', required=True,
                    help='input file path')
parser.add_argument('-u', action='store', dest='url',
                    metavar='URL', required=True,
                    help='api url, e.g.: "http://localhost:8000/reports/api/v1"')
parser.add_argument('-auth', action='store', dest='auth',
                    metavar='AUTH', required=True,
                    help='api authentication to use <user:password>')

if __name__ == "__main__":
    args = parser.parse_args()
 
    log_level = logging.WARNING # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    if args.verbose:
        # NOTE: when running with the django settings file, the logging configured there will augment this, and 
        # cause double logging. So this will manually override that.
        # Probably a better solution would be to configure this utility as a "management command"
        # and then let manage.py run it.  see: https://docs.djangoproject.com/en/1.4/howto/custom-management-commands/
        logging.basicConfig(level=log_level, format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')   
        logger.setLevel(log_level)
        
    main(args.inputDirectory, args.url, args.auth)