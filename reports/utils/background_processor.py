import argparse
import datetime
import getpass
import json
import logging
from multiprocessing.process import Process
import os
import socket
import time
import urllib

import django
from django.conf import settings
from django.conf.global_settings import DEFAULT_CHARSET
from django.test.client import RequestFactory
import django.test.client
import django.utils.http

from reports import HEADER_APILOG_COMMENT
from reports.schema import API_RESULT_DATA
import reports.schema as SCHEMA
from reports.serialize import MULTIPART_MIMETYPE, JSON_MIMETYPE
from reports.serializers import DJANGO_ACCEPT_PARAM, LimsSerializer
from reports.utils import parse_credentials


logger = logging.getLogger(__name__)

JOB = SCHEMA.JOB

DEBUG_BACKGROUND = True or logger.isEnabledFor(logging.DEBUG)

def create_request_from_job(job_data, raw_data=''):
    
    if DEBUG_BACKGROUND is True:
        logger.info('create_request_from_job: %r', job_data)
    
    request_factory = RequestFactory()
    
    path = job_data[JOB.URI]
    method = job_data[JOB.METHOD]
    encoding = job_data[JOB.ENCODING]
    content_type = job_data[JOB.CONTENT_TYPE]
    # For job processing, force response to JSON
    #     accept = urllib.unquote(job_data[JOB.HTTP_ACCEPT])
    accept = JSON_MIMETYPE
    comment = job_data[JOB.COMMENT]
    params = job_data[JOB.PARAMS]
    if params:
        params = json.loads(params)
    else:
        params = {}
    
    if raw_data:
        if DEBUG_BACKGROUND is True:
            logger.info('add raw data: %d', len(raw_data))
        if MULTIPART_MIMETYPE not in content_type:
            msg = 'content type must contain %r for raw data post: found: %r'\
                 % (MULTIPART_MIMETYPE, content_type)
            logger.error(msg)
            raise ValidationError(key=JOB.CONTENT_TYPE, msg=msg)
        if method != 'POST':
            errmsg = 'method %r is not %r, required for raw data post' % (
                method, 'POST')
            raise ValidationError(key=JOB.METHOD, msg=errmsg)
    else:
        if DEBUG_BACKGROUND is True:
            logger.info('no raw data to add')
        
    if DEBUG_BACKGROUND is True:
        logger.info('create_request_from_job content type: %r', content_type)
    request = request_factory.generic(
        method, path, data=raw_data, HTTP_ACCEPT=accept,
        content_type=content_type,**params )
    
    if DEBUG_BACKGROUND is True:
        logger.info('create_request_from_job: META: %r', request.META )
        logger.info('create_request_from_job: FILES: %r', request.FILES )
    
    return request
    

def get_django_request_job_params(request):
    '''
    Parse Job parameters from the request:
    @param request django.http.request.HttpRequest
    NOTE: 
    - because authentication is run before this method, Django has 
    processed the request already: the request object has members:
    request.GET, request.POST, request.FILES
    as well as the raw_post_data: _body
    '''
    
    path = request.path
    method = request.method
    encoding = request.encoding or DEFAULT_CHARSET
    content_type = request.META.get('CONTENT_TYPE', '')
    accept = request.META.get(DJANGO_ACCEPT_PARAM, '*/*')
    username = request.user.username
    comment = request.META.get(HEADER_APILOG_COMMENT, '')
    context_data = {}
    # Note: cookies are not stored -- yet 20180402
    
    # query_string
    query_params = {}
    if request.GET:
        for k,v in request.GET.items():
            query_params[k] = v
    
    if DEBUG_BACKGROUND is True:
        logger.info('user: %r, path: %r, method: %r, encoding: %r, '
            'content_type: %r, accept: %r, params: %r',
            username, path, method, encoding, content_type, 
            accept, query_params)
    
    post_data = {}
    if request.POST:
        if MULTIPART_MIMETYPE not in content_type:
            logger.warn('background processor POST data, content_type '
                'does not contain %r, %r', MULTIPART_MIMETYPE, content_type)
        for key,val in request.POST.items():
            logger.info('found post param: %r: %r', key, val)
            post_data[key] = val
    if request.FILES:
        context_data['filenames'] = []
        if MULTIPART_MIMETYPE not in content_type:
            logger.warn('background processor POST data, content_type '
                'does not contain %r, %r', MULTIPART_MIMETYPE, content_type)
        for key,_file in request.FILES.items():
            logger.info('set post form file: %r, %r', key, _file)
            post_data[key] = _file
            context_data['filenames'].append(_file.name)

    raw_post_data = None
    if post_data:
        logger.info('encode multipart post data: %r', post_data.keys())
        content_type = 'multipart/form-data; boundary=%s' % django.test.client.BOUNDARY
        # 'multipart/form-data; boundary=BoUnDaRyStRiNg'
        raw_post_data = django.test.client.encode_multipart(
            django.test.client.BOUNDARY, post_data)

    job_data = {
        JOB.USERNAME: username,
        JOB.URI: path,
        JOB.METHOD: method,
        JOB.ENCODING: encoding,
        JOB.CONTENT_TYPE: content_type,
        JOB.HTTP_ACCEPT:  urllib.quote(accept),
        JOB.PARAMS: json.dumps(query_params),
        JOB.COMMENT: comment,
        JOB.CONTEXT_DATA: json.dumps(context_data)
    }
    
    return job_data, raw_post_data

def get_process_env():
    ''' Information about the current process environment '''
    
    data = {}
    try:
        os_props = ['pid', 'uid', 'uname']
        for prop in os_props:
            func = getattr(os, 'get%s'%prop, None)
            logger.info('get os prop: %s, %r', prop, func)
            if func is not None:
                data[prop] = func()
        data['user'] = getpass.getuser()
        data['hostname'] = socket.gethostname()
    except Exception, e:
        logger.exception('on get_process_env')
        
    return data

class ApiClient(django.test.client.Client):
    
    def __init__(self, username, password, serializer=None):
        self.serializer = serializer
        self.username = username
        self.password = password
        
        if not self.serializer:
            self.serializer = LimsSerializer()

        super(ApiClient, self).__init__(
            enforce_csrf_checks=False,
            SERVER_NAME='localhost')

    def deserialize(self, resp):

        return self.serializer.deserialize(
            self.serializer.get_content(resp), resp['Content-Type'])

    def create_basic(self, username, password):
        """
        Creates & returns the HTTP ``Authorization`` header for use with BASIC
        Auth.
        """
        import base64
        return 'Basic %s' % base64.b64encode(
            ':'.join([username, password]).encode('utf-8')).decode('utf-8')
    
    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)

    def request(self, **request):
        
        if DJANGO_ACCEPT_PARAM not in request:
            request[DJANGO_ACCEPT_PARAM] = JSON_MIMETYPE
        
        if 'HTTP_AUTHORIZATION' not in request:
            request['HTTP_AUTHORIZATION'] = self.get_credentials()
        
        return django.test.client.Client.request(self, **request)

class BackgroundClient(object):
    
    def __init__(self, api_client):
        self.api_client = api_client
    
    def service(self, job_id):
        self._service(job_id)
        
    def _service(self, job_id):
        
        uri = '/' + '/'.join(
            [SCHEMA.REPORTS_API_URI, JOB.resource_name,str(job_id)])
        data = {}
        kwargs = {}
        logger.info('get the job: %r', uri)
        resp = self.api_client.get(
            uri, data=data, **kwargs)
        logger.info('get job response status: %r', resp.status_code)
        job_data = self.api_client.deserialize(resp)
        if API_RESULT_DATA in job_data:
            job_data = job_data[API_RESULT_DATA]
        
        logger.info('got job data: %r', job_data)
        method = job_data[JOB.METHOD]
        uri = job_data[JOB.URI]
        job_id = job_data[JOB.ID]
        query_params = { JOB.JOB_PROCESSING_FLAG: job_id }

        # Create ApiClient instance for each job, so that the 
        # user assigned to the job is used
        as_user = '%s:%s' % (job_data[JOB.USERNAME], self.api_client.username)
        logger.info('process job: %r for user: %s', job_id, as_user )
        user_api_client = ApiClient(as_user, self.api_client.password)
        params = {
            'QUERY_STRING': django.utils.http.urlencode(query_params, doseq=True)
        }
        params['HTTP_AUTHORIZATION'] = self.api_client.get_credentials()
        # data will be filled in on the server
        # content type will be modified as needed on the server
        content_type = JSON_MIMETYPE
        logger.info('re-process the original request...')
        resp = user_api_client.generic(method, uri, content_type, **params)
        logger.info('job service response status: %r', resp.status_code)
        job_service_response = self.api_client.deserialize(resp)
        logger.info('job_service_response: %r', job_service_response)
        
        return job_service_response
        
    
class BackgroundProcessor(object):
    ''' 
    Process the the queue of background jobs:
    - fire of a BackgroundClient job for "pending" jobs
    - manage state of jobs
    '''
    
    def __init__(self, api_client):
        self.api_client = api_client
        self.client_processor = BackgroundClient(api_client)
    
    def service_background_jobs(self):
        logger.info('service_background_jobs')
        
        # NOTE: paths must begin with a "/", indicating that the first part of 
        # the URI is a script name (which each app, i.e "reports" serves as).
        # see django.core.handlers.wsgi.__init__
        uri = '/' + '/'.join([SCHEMA.REPORTS_API_URI, JOB.resource_name ])
        data = {
            JOB.STATE: SCHEMA.VOCAB.job.state.PENDING }
        kwargs = {}

        logger.info('get jobs: %r', uri)
        resp = self.api_client.get(
            uri, data=data, **kwargs)
        
        job_listing = self.api_client.deserialize(resp)
        if API_RESULT_DATA in job_listing:
            job_listing = job_listing[API_RESULT_DATA]
        
        for job in job_listing:
            logger.info('found job: %r', job)
            job_id =job[JOB.ID]
            logger.info('Process the job: %r', job_id)
            p = Process(target=self.client_processor.service,args=(job_id,) )
            # make the parent process wait: 
            # p.daemon = True
            # if set to true, then the parent process won't wait.  
            logger.info('start')
            p.start();
            logger.info('started...')
            
        logger.debug('servicing completed')

    def service(self, loop_time_seconds=30, time_to_run=60):
        logger.info('service: loop_time_seconds: %d, time_to_run: %d', 
            loop_time_seconds, time_to_run)
        time_start = time.time()
        while time.time()-time_start < time_to_run:
            self.service_background_jobs()
            time.sleep(loop_time_seconds)
        logger.info('service loop exit: run time (s): %d', (time.time()-time_start))    
      

parser = argparse.ArgumentParser(description=
    'Process and update background jobs')
# Admin credentials
parser.add_argument(
    '-U', '--username',
    help='username for the api authentication')
parser.add_argument(
    '-p', '--password',
    help='user password for the api authentication')
parser.add_argument(
    '-c', '--credential_file',
    help='credential file containing the username:password for api '
    'authentication')

parser.add_argument(
    '-lt', '--service_loop_time_s', required=True, type=int, 
    help='Number of seconds between each check for new pending jobs in the database')
parser.add_argument(
    '-rt', '--run_time_s', required=True, type=int, 
    help='Number of seconds to run and then exit')

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
    logging.basicConfig(
        level=log_level, 
        format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')        

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

    print 'run background processing service', datetime.datetime.now()
    try:
        django.setup()
        
        client = ApiClient(username, password)
        background_processor = BackgroundProcessor(client)
        background_processor.service(loop_time_seconds=args.service_loop_time_s, 
            time_to_run=args.run_time_s);  
        
    except Exception, e:
        logger.exception('in background service method')
        raise e
    print 'exit background processing service', datetime.datetime.now()
    