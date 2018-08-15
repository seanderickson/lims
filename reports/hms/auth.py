# HMS ICCB-L Ecommons Authentication


import argparse
import logging
import requests
import re
from datetime import timedelta
from datetime import datetime
from urllib import unquote_plus
import getpass

import ldap3

logger = logging.getLogger(__name__)


LDAP_SERVER = 'gal.med.harvard.edu'
LDAP_USER_PREFIX = 'med\\'

def ldap_authenticate(ecommons_id, ecommons_password):

    try:
        logger.info('attempt to use ldap server %r to authenticate: %r', 
                    LDAP_SERVER, ecommons_id)
        server = ldap3.Server(LDAP_SERVER, use_ssl=True)
        conn = ldap3.Connection(
            server, 
            user=LDAP_USER_PREFIX + ecommons_id, 
            password=ecommons_password, auto_bind=True)
        logger.info('connection success: who_am_i: %r', 
                    conn.extend.standard.who_am_i())
        return True
    except ldap3.LDAPBindError, e:
        logger.info('ldap3.LDAPBindError: %r', e)
        raise
    except Exception, e:
        logger.exception('LDAP exception: %r', e)
        raise

def authenticate(ecommons_id, ecommons_password):
    '''
    Authenticate with the eCommons Server using the credentials
    '''
    return ldap_authenticate(ecommons_id, ecommons_password)

def authenticateOld(ecommons_id, ecommons_password):
    """
    Performs authentication of a user using their eCommons ID and password. 

    Sends an HTTP request to an HMS-hosted eCommons authentication server, 
    using a SAML [1] inspired  format. 
    See also the HMS RITG wiki [2] for more details.

    This code must be executed on a host with an IP address that is explicitly
    registered with the eCommons authentication server (orchestra and trumpet are
    approved).  You will also need a valid "issuer ID".
        
    [1] http://www.oasis-open.org/committees/tc_home.php?wg_abbrev=security
    [2] https://wiki.med.harvard.edu/IT/RITG/AuthenticationWebService
        TODO: this link is stale, get a new link from RITG
    [3] https://ddh4.dev.cbmi.med.harvard.edu/java/ecommons-authen/0.0.1/
    [4] https://ddh4.dev.cbmi.med.harvard.edu/java/ecommons-authen/0.0.1/doc/HMSAuthenticationWebService.20090429.doc
    """
    
    PATTERN_STATUS_CODE = r'.*<StatusCode>(\d+)</StatusCode'
    PATTERN_STATUS_CODE_CATEGORY = r'.*<StatusCodeCategory>(\w+)</StatusCodeCategory'
    PATTERN_STATUS_MESSAGE = r'.*<StatusMessage>([^<]+)</StatusMessage'

    url = 'https://authenticate.med.harvard.edu/wsAuthenticate.asp'
            
    authreq = """<AuthNRequest>
<RequestId>{9BA34ADF-287E-47CE-8F61-8795943FEC1C}</RequestId>
<Issuer>Orchestra_ATTR_CLIENT</Issuer>
<Signature>ExampleRequestSignature</Signature>
<IssueInstant>%s</IssueInstant>
<ValidityInterval>
<NotBefore>%s</NotBefore>
<NotAfter>%s</NotAfter>
</ValidityInterval>
<RequestApp>HMSICCBL</RequestApp>
<AuthNData>
<Id>%s</Id>
<Password>%s</Password>
</AuthNData>
</AuthNRequest>"""

    logger.info( 'Authenticating %s ...' % ecommons_id) 
    
    headers = {'Content-Type':'text/xml',
               'SOAPMethodName':'urn:myserver:AuthenticationReply#GetXIDAuthenticateResponse'}
    issue_instant = datetime.now()
    not_before = issue_instant
    not_after = issue_instant + timedelta(days=7)
    prepared_req = authreq % (
        issue_instant, not_before, not_after, ecommons_id, ecommons_password)
    r = requests.post(url, data=prepared_req, headers=headers, timeout=10 )
    if r.status_code != 200: 
        raise Exception('HTTP response code: %r, %r' % (r.status_code, r.content()))

    matchObject = re.match(PATTERN_STATUS_CODE, r.text)
    if matchObject:
        status=matchObject.group(1)
        if status=='2000':
            logger.info('successful authentication for %r',ecommons_id)
            return True
        else:
            logger.warn(str((
                'failed authentication', status,'ecommons_id', ecommons_id,
                 re.match(PATTERN_STATUS_CODE_CATEGORY,r.text).group(1),
                 unquote_plus(str(re.match(PATTERN_STATUS_MESSAGE,r.text).group(1))) )))
            return False
    else:
        msg =str((
            'ecommons authentication request not recognized (cannot find the'
            ' status, looking for the pattern)', r.text, PATTERN_STATUS_CODE))
        logger.error(msg)
        raise Exception(msg)

parser = argparse.ArgumentParser(description='Authenticate User')
parser.add_argument('-i', action='store', dest='id',
                    metavar='ID', required=True,
                    help='eCommons ID')
parser.add_argument('-p', action='store', dest='pw',
                    metavar='PW',
                    help='eCommons password')
parser.add_argument('-v', '--verbose', dest='verbose', action='count',
                help="Increase verbosity (specify multiple times for more)")    

if __name__ == "__main__":
    args = parser.parse_args()

    log_level = logging.WARNING # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(
        level=log_level, 
        format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')        
    logger.setLevel(log_level)

    password = args.pw
    if not password:
        password = getpass.getpass()
    
    authenticate(args.id,password)
