from __future__ import unicode_literals 
import re
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

# Sorting based on https://nedbatchelder.com/blog/200712/human_sorting.html
def tryfloat(s):
    try:
        return float(s)
    except:
        return s
     
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryfloat(c) for c in re.split('(-?[0-9\.]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    return sorted(l, key=alphanum_key)
