import requests
from requests import Request
import argparse
import logging
from urlparse import urlparse
import sys
import getpass

DEBUG=False
LOGIN_FORM = '/db/login/'
ACCOUNTS_LOGIN_FORM = '/accounts/login/'
CONTENTTYPE_HEADERS =   { 
    'json': {'content-type': 'application/json'},
    'csv':  {'content-type': 'text/csv'},
    } 

logger = logging.getLogger(__name__)


def get_logged_in_session(username, password, base_url, 
                          login_form=None,
                          headers={}):
    ''' 
    Preferred form, first gets the csrf token from the account form, then uses
    this to perform subsequent operations.
    '''
    if not login_form:
        login_form=ACCOUNTS_LOGIN_FORM
    
    if DEBUG:
        print 'username:%s:%s,%s,%s' % (username,password,base_url,login_form)
    s = requests.Session()
    # first, get the page, to induce django to generate the csrf token
    r = s.get('%s%s' % (base_url, login_form), verify=False)
    
    if not r.status_code in [200]:
        raise Exception("Error getting login form: %s, response: %s, %s" 
                        % (login_form, r.status_code, r.content))        
    # verify 
    if not 'csrftoken' in s.cookies:
        raise Exception('failed to get a csrf token from form %s\n%s' 
                        % (login_form, r.content))
    if DEBUG:
        sys.stderr.write('2. csrftoken: %s\n' % s.cookies['csrftoken'])
    
    # with the csrftoken, we will be allowed to post    
    r = Request('POST', '%s%s' % (base_url, login_form), 
                data={'username':username,'password':password },
                headers={'Referer': '%s%s' % (base_url, login_form),
                         'X-CSRFToken': s.cookies['csrftoken'] })

    if DEBUG:
        sys.stderr.write('request: %s\n' % str((r)) )
    
    prepped = s.prepare_request(r)
    
    # Post the login request to the url
    logger.info(str(('post the login request', r.url, r.method)))
    r = s.send(prepped)
    logger.info(str(('response', r.status_code)))
    if DEBUG:
        sys.stderr.write('login response: %s\n' % r.content )
    
    if r.status_code in [200]:
        if 'logged in as:' not in r.content:
            if DEBUG:
                sys.stderr.write('login failed:\n%s\n' % r.content )
            raise Exception("Authentication failed, username/password rejected.")
        return s
    else:
        raise Exception("Authentication error: status: %s, %s" 
                        % (r.status_code, r.content))


def get_session(username=None, password=None, base_url=None, login_form=LOGIN_FORM):
    '''
    deprecated - simplified version, if allowed to make the post to the 
    login form without a
    csrftoken (not the default, in this case, can be done with a custom 
    view.loging, with csrf_exempt decorator.
    '''
    s = requests.Session()
    r = s.post('%s%s' % (base_url, login_form), 
               {'username': username,'password':password } , verify=False)
    
    if DEBUG:
        sys.stderr.write('status: %s\n' %r.status_code)
    if not 'csrftoken' in s.cookies:
        raise Exception('failed to get a csrf token from form %s' 
                        % (login_form, r.content))
    if DEBUG:
        sys.stderr.write('1. csrftoken: %s\n' % s.cookies['csrftoken'])
    
    if r.status_code in [200]:
        return s
    else:
        raise Exception("Authentication error: status: %s, %s" 
                        % (r.status_code, r.content))


def delete(url, request_or_session, headers):
    r = request_or_session.delete(url,headers=headers);
    
    if r.status_code not in [204]:
        raise Exception("Error: status: %s, %s" 
                        % (r.status_code, r.content))
    logger.info(str(('DELETE', url, r.status_code)))
    return r


def get(url,request_or_session,headers): 
    logger.info(str(('GET', url)))
    r = request_or_session.get(url,headers=headers)
    
    if r.status_code not in [200]:
        raise Exception("Error: status: %s, %s" 
                        % (r.status_code, r.content))
    return r
       
    
def put(url,request_or_session,headers ):
    with open(file) as f:
        
        headers['Referer']=url
        headers['X-CSRFToken'] = request_or_session.cookies['csrftoken']
        if DEBUG:
            sys.stderr.write('csrftoken: %s\n' % request_or_session.cookies['csrftoken'])
            sys.stderr.write('headers: %s\n' % str((headers)))
        
        r = Request('PUT', url,
                    headers=headers,
                    data=f.read())
        
        prepped = request_or_session.prepare_request(r)
    
        r = request_or_session.send(prepped)        
            
        if r.status_code not in [200,202]:
            raise Exception("Error: status: %s\n%s\n%s" 
                            % (r.status_code, r.headers, r.content))
        return r
    
def patch(url,request_or_session,file, headers={} ):
    with open(file) as f:
        
        headers['Referer']=url
        headers['X-CSRFToken'] = request_or_session.cookies['csrftoken']
        if DEBUG:
            sys.stderr.write('csrftoken: %s\n' % request_or_session.cookies['csrftoken'])
            sys.stderr.write('headers: %s\n' % str((headers)))
        
        r = Request('PATCH', url,
                    headers=headers,
                    data=f.read())
        
        prepped = request_or_session.prepare_request(r)
    
        r = request_or_session.send(prepped)        
            
        if r.status_code not in [200,202]:
            raise Exception("Error: status: %s\n%s\n%s" 
                            % (r.status_code, r.headers, r.content))
        return r

    
    
parser = argparse.ArgumentParser(description='url')
parser.add_argument('url', help='url to connect to')
parser.add_argument('-u', '--username', help='username', required=True)
parser.add_argument('-p', '--password', help='password', required=False)
parser.add_argument('-f', '--file', help='file', required=False)
parser.add_argument('-a', '--action', help='HTTP action', required=True,
                    choices=['GET','PUT','PATCH','DELETE'])
parser.add_argument('--header', action='append')
parser.add_argument('-v', '--verbose', dest='verbose', action='count',
                help="Increase verbosity (specify multiple times for more)")    
parser.add_argument('--login_form', help='login form to use, i.e. /accounts/login/', 
                    required=False)

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

    password = args.password
    if not password:
        password = getpass.getpass()
    
    # Default headers
    headers = {'content-type': 'application/json',
               'Accept': 'application/json'}

    if args.header:
        for x in args.header:
            if ':' not in x:
                raise Exception("--header must be of the form 'header:value' ")
            header = x.split(":")
            headers[header[0]]=header[1]
    
    logger.info(str(('headers: ', headers)))

    url = args.url
    u = urlparse(url)
    base_url = '%s://%s' % (u.scheme,u.netloc)
    s = get_logged_in_session(args.username,password,base_url,login_form=args.login_form)
    
    action = args.action
    
    if action == 'GET':
        r = get(url, s, headers)
        print r.content
    elif action == 'DELETE':
        delete(url, s, headers)
        print 'DELETE: ', url, ' ,response:', r.status_code
    elif action == 'PATCH':
        r = patch(url,s,args.file,headers)
        print r.content
    elif action == 'PUT':
        r = put(url,s,args.file,headers)        
        print r.content
    else:
        raise Exception("unknown action %s" % action)
    
    
