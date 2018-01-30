
import logging

logger = logging.getLogger(__name__)

def parse_credentials(credential_file):
    username, password = None,None
    
    with open(credential_file) as f:
        for line in f:
            if username is not None:
                raise ArgumentError(credential_file, 
                    'must contain a single "username:password"')
            (username,password) = line.strip().split(':')
            logger.debug('read: %r:has_password:%r', 
                username, password is not None)
    return username, password
