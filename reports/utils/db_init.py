import json
import logging
import requests
import sys, os
from urlparse import urlparse
import argparse
import getpass

from django_requests import get_logged_in_session

import reports.utils.serialize

logger = logging.getLogger(__name__)

# TODO: this replaces /reports/management/commands/db_init.py - sde4 - 201404


class ApiError(Exception):
    
    def __init__(self, url, action, result):
        err_msg = ''
        try:
            json_status = result.json()
            if json_status and 'error_message' in json_status:
                err_msg = json_status['error_message']
        except ValueError,e:
            logger.warn('There is no json in the response')
            logger.warn(str(('-----raw response text-------', result.text)) )
            err_msg = result.content

        self.message = str((
            url,'action',action, result.reason, result.status_code, err_msg )) \
            .replace('\\n','') \
            .replace('\\','')
        if(logger.isEnabledFor(logging.DEBUG)):
            self.message = str((url,'action',action, result.reason, 
                                result.status_code, result.content )).replace('\\','')

    def __str__(self):
        return str(self.message)


def delete(obj_url, headers, session=None, authentication=None):
    try:
        if session:
            r = session.delete(obj_url, headers=headers,verify=False)
        elif authentication:
            r = requests.delete(obj_url, auth=authentication, headers=headers,verify=False)
        
        if(r.status_code != 204):
            print "DELETE ERROR", r, r.text
            raise ApiError(obj_url,'DELETE',r)
        print 'DELETE: ', obj_url, ' ,response:', r.status_code
        logger.info(str(('DELETE', obj_url)))
    except Exception, e:
        logger.error(str(('delete', obj_url, 'exception recorded while contacting server', e)))
        raise e

def put(input_file, obj_url,headers, session=None, authentication=None):
    try:
        with open(input_file) as f:
            if session:
                r = session.put(
                    obj_url, headers=headers, data=f.read(),verify=False)
            elif authentication:
                r = requests.put(
                    obj_url, auth=authentication, headers=headers, data=f.read(),verify=False)
            if(r.status_code != 200):
                raise ApiError(obj_url,'PUT',r)
            print ('PUT: ' , input_file, 'to ', obj_url,' ,response:', 
                   r.status_code, ', count: ',len(r.json()['objects']))
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
        logger.warn(str((
            'throw', e, tb.tb_frame.f_code.co_filename, 'error line', 
            tb.tb_lineno, extype, ex)))
        logger.error(str(('put', obj_url, 'exception recorded while contacting server', e)))
        raise e
    
def patch(patch_file, obj_url,headers, session=None, authentication=None):
    try:
        print 'PATCH: ' , patch_file, 'to ', obj_url
        with open(patch_file) as f:
            if session:
                r = session.patch(
                    obj_url, headers=headers, data=f.read(),verify=False)
            elif authentication:
                r = requests.patch(
                    obj_url, auth=authentication, headers=headers, data=f.read(),verify=False)
            if(r.status_code not in [200,202,204]):
                
                raise ApiError(obj_url,'PATCH',r)
            print ('PATCH: ', patch_file, ', to: ',obj_url,' ,response:', 
                    r.status_code,', count: ',len(r.json()['objects']))
            if(logger.isEnabledFor(logging.DEBUG)):
                logger.debug('--- PATCHED objects:')
                try:
                    for obj in r.json()['objects']:
                        logger.debug(str((obj)))
                except ValueError,e:
                    logger.debug('----no json object to report')
                    logger.debug(str(('text response', r.text)))
    except Exception, e:
        logger.error(str(('patch', obj_url, 'exception recorded while contacting server', e)))
        raise e


parser = argparse.ArgumentParser(description='url')
parser.add_argument(
    '-u', '--url', required=True,
    help='api url, e.g.: "http://localhost:8000/reports/api/v1')
parser.add_argument(
    '-d', '--input_dir', required=True,
    help='input file path (where the input actions files can be found)')
parser.add_argument(
    '-f', '--input_actions_file', required=True,
    help='''input actions file.
            An input action is a [command,resource,file]; 
            see the reports/static/api_init/api_init_actions.csv file''')
parser.add_argument(
    '-U', '--username', required=True,
    help='username for the api authentication')
parser.add_argument(
    '-p', '--password',
    help='user password for the api authentication')
parser.add_argument(
    '-v', '--verbose', dest='verbose', action='count',
    help="Increase verbosity (specify multiple times for more)")    

if __name__ == "__main__":
    args = parser.parse_args()
    log_level = logging.WARNING # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
        DEBUG=True
	from reports.utils.django_requests import DEBUG as DJANGO_REQUESTS_DEBUG
	DJANGO_REQUESTS_DEBUG=True
    logging.basicConfig(
        level=log_level, 
        format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')        

    CONTENT_TYPES =   { 
        'json': {'content-type': 'application/json'},
        'csv':  {'content-type': 'text/csv'},
        } 
                
    url = args.url
    u = urlparse(url)
    base_url = '%s://%s' % (u.scheme,u.netloc)

    password = args.password
    if not password:
        password = getpass.getpass()
    
    headers ={}

    #### log in using django form-based auth, and keep the sesion
    session = get_logged_in_session(
        args.username, password, base_url)
    # django session based auth requires a csrf token
    headers['X-CSRFToken'] = session.cookies['csrftoken']
    
    
    with open(args.input_actions_file) as input_file:
        api_init_actions = reports.utils.serialize.from_csv(input_file)

        for action in api_init_actions:
            
            print '\n++++=========== processing action', json.dumps(action)
            command = action['command'].lower() 
            resource = action['resource'].lower()
            resource_uri = url + '/' + resource
            
            # tastypie session based auth requires the referer to be set
            headers['Referer']=resource_uri
            
            if command == 'delete':
                headers.update(CONTENT_TYPES['json'])
                delete(resource_uri, headers, session=session)
            
            else:
                data_file = os.path.join(args.input_dir,action['file'])
                print 'Data File: ', data_file
                extension = os.path.splitext(data_file)[1].strip('.')
                headers.update(CONTENT_TYPES[extension])
                        
                if command == 'put':
                    put(data_file, resource_uri, headers, session=session)
                        
                elif command == 'patch':
                    patch(data_file, resource_uri, headers, session=session)
  
                else:
                    raise Exception(str((
                        'Only [ PUT, PATCH, DELETE ] are supported.',
                        '- cannot POST multiple objects to tastypie; ',
                        'therefore the "post" command is invalid with the initialization scripts',
                        'resource entry: ',json.dumps(action) )) )


