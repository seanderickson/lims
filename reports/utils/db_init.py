'''
utils.db_init:

This module implements API requests using the requests package:
 - actions specified in an api_init_actions" csv file.
 
@see static/api_init/api_init_actions.csv

'''
from __future__ import unicode_literals

import argparse
import getpass
import json
import logging
import sys, os
from urlparse import urlparse

import requests

import django_requests
import reports.schema as SCHEMA
import reports.serialize.csvutils as csvutils
from reports.serializers import CONTENT_TYPES
from reports.utils import parse_credentials


logger = logging.getLogger(__name__)

class ApiError(Exception):
    
    def __init__(self, url, action, result):
        err_msg = ''
        try:
            err_msg = result.json()
        except ValueError,e:
            err_msg = str(result.content)

        self.message = {
            'url': url,
            'action': action, 
            'reason': result.reason,
            'status_code': result.status_code,
            'message': err_msg
            }

    def __str__(self):
        return str(self.message)

def _print_result(r):
    result = r.json()
    result = result.get(SCHEMA.API_RESULT_META, result)
    return result.get(SCHEMA.API_MSG_RESULT, result)
    

def delete(obj_url, headers, session=None, authentication=None):
    if session:
        r = session.delete(obj_url, headers=headers,verify=False)
    elif authentication:
        r = requests.delete(
            obj_url, auth=authentication, headers=headers,verify=False)
    if not r.status_code < 300:
        raise ApiError(obj_url,'DELETE',r)

    print 'DELETE: {}, response: {}'.format(
        obj_url, r.status_code)
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug('result: %r', r)

def put(input_file, obj_url,headers, session=None, authentication=None):
    with open(input_file) as f:
        if session:
            r = session.put(
                obj_url, headers=headers, data=f.read(),verify=False)
        elif authentication:
            r = requests.put(
                obj_url, auth=authentication, 
                headers=headers, data=f.read(),verify=False)
        if not r.status_code < 300:
            raise ApiError(obj_url,'PUT',r)
        print 'PUT: {} to {}, response: {}, count: {}, result: {} '.format(
            input_file, obj_url, r.status_code, 
            len(r.json().get(SCHEMA.API_RESULT_DATA, [])),
            _print_result(r))

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug('result: %r', r.json())

def patch(input_file, obj_url,headers, session=None, authentication=None):
    with open(input_file) as f:
        if session:
            r = session.patch(
                obj_url, headers=headers, data=f.read(),verify=False)
        elif authentication:
            r = requests.patch(
                obj_url, auth=authentication, headers=headers, 
                data=f.read(),verify=False)
        if not r.status_code < 300:
            raise ApiError(obj_url,'PATCH',r)
        
        print 'PATCH: {} to {}, response: {}, result: {} '.format(
            input_file, obj_url, r.status_code, 
            _print_result(r))
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug('result: %r', r.json())


def post(input_file, extension, obj_url,headers, session=None, authentication=None):
    with open(input_file) as f:
        logger.info('headers: %r', headers)
        # NOTE: the extension as the files key is a ICCBL LIMS specific
        files = { extension: f }
        
        if session:
            print 'POST with session auth'
            r = session.post(
                obj_url, headers=headers, files=files,verify=False)
        elif authentication:
            print 'POST with basic auth'
            r = requests.post(
                obj_url, auth=authentication, headers=headers, 
                files=files,verify=False)
        if not r.status_code < 300:
            raise ApiError(obj_url,'PATCH',r)
        
        print 'POST: {} to {}, response: {}, result: {} '.format(
            input_file, obj_url, r.status_code, 
            _print_result(r))
        
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug('result: %r', r.json())


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
    '-U', '--username', 
    help='username for the api authentication')
parser.add_argument(
    '-p', '--password',
    help='user password for the api authentication')
parser.add_argument(
    '-c', '--credential_file',
    help = 'credential file containing the username:password')
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
        django_requests.DEBUG = True
    logging.basicConfig(
        level=log_level, 
        format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')        
                
    url = args.url
    u = urlparse(url)
    base_url = '%s://%s' % (u.scheme,u.netloc)
    
    username = None
    if args.credential_file:
        username,password = parse_credentials(args.credential_file)
    if username is None:
        username = args.username
        if username is None:
            parser.error(
                'username is required if not specifying the credential_file')
        password = args.password
        if not password:
            password = getpass.getpass()
    
    session_headers ={}

    logger.info('begin processing file: %r', args.input_actions_file)
    session = django_requests.get_logged_in_session(
        username, password, base_url)
    # Django session based auth requires a csrf token
    session_headers['X-CSRFToken'] = session.cookies['csrftoken']
    # Always accept json for debugging the returned values
    session_headers['Accept'] = 'application/json'
    
    with open(args.input_actions_file) as input_file:
        api_init_actions = csvutils.from_csv(input_file)

        for action in api_init_actions:
            
            print '\n++++=========== processing action', json.dumps(action)
            command = action['command'].lower() 
            resource = action['resource'].lower()
            resource_uri = url + '/' + resource
            
            headers = session_headers.copy()
            headers['Referer']=resource_uri
            
            if command == 'delete':
                delete(resource_uri, headers, session=session)
            else:
                data_file = os.path.join(args.input_dir,action['file'])
                print 'Data File: ', data_file
                extension = os.path.splitext(data_file)[1].strip('.')
                
                if extension not in CONTENT_TYPES:
                    raise Exception('Unknown file extension: %r, types: %r'
                        % (extension, CONTENT_TYPES))
                        
                if command == 'put':
                    headers.update({ 'content-type': CONTENT_TYPES[extension] })
                    put(data_file, resource_uri, headers, session=session)
                elif command == 'patch':
                    headers.update({ 'content-type': CONTENT_TYPES[extension] })
                    patch(data_file, resource_uri, headers, session=session)
                elif command == 'post':
                    post(data_file, extension, resource_uri, headers, session=session)
                else:
                    raise Exception(str((
                        'Only [ PATCH, POST, PUT, DELETE ] are supported.',
                        'resource entry: ',json.dumps(action) )) )


