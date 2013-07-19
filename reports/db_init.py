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
    try:
        r = requests.delete(obj_url, auth=authentication, headers=headers)
        print 'Deleted metahash table'
    except Exception, e:
        logger.error(str(('exception recorded while contacting server', e)))
        raise e
    if(r.status_code != 204):
        raise ApiError('DELETE')

    # fyi, system way to do this, has been replaced above    
    #pipe = Popen(["curl", "-H 'Content-Type: application/json' --user sde4:3xc3l3nc -X DELETE  http://localhost:8000/reports/api/v1/metahash/",
    #             shell=True, stdin=PIPE).stdin
        
    # initial import
    api_importer.main(input_file_dir + '/metahash_fields_initial.json', url_base + '/metahash/',  'POST', auth=':'.join(authentication))

    # now that the field defs are created, including the json fields, we can populate these json field values
    api_importer.main(input_file_dir + '/metahash_fields_initial_patch.json',url_base + '/metahash/',  'PATCH', auth=':'.join(authentication),  id_key='key')

    # define vocabularies
    api_importer.main(input_file_dir + '/metahash_fields_vocabularies.json',url_base + '/metahash/',  'POST', auth=':'.join(authentication))

    # define screensaveruser
    api_importer.main(input_file_dir + '/metahash_fields_screensaveruser.json',url_base + '/metahash/',  'POST', auth=':'.join(authentication))

    # delete the vocabularies records
    obj_url = url_base + '/vocabularies/'
    try:
        r = requests.delete(obj_url, auth=authentication, headers=headers)
        print 'Deleted vocabularies table'
    except Exception, e:
        logger.error(str(('exception recorded while contacting server', e)))
        raise e
    if(r.status_code != 204):
        raise ApiError(r,'DELETE')

    # import vocabularies
    api_importer.main(input_file_dir + '/vocabularies_data.json',url_base + '/vocabularies/',  'POST', auth=':'.join(authentication))

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