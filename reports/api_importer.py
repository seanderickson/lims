import argparse
import logging
import requests
import json

from dump_obj import dumpObj

logger = logging.getLogger(__name__)

#DEFAULT_ID_KEY='id'

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
    

def main_csv(inputFilePath, url, action, auth=None ):
    pass
                    
def main(inputFilePath, url, action, auth=None): #, id_key=DEFAULT_ID_KEY ):
    logger.info('input path:' + inputFilePath + ', ' + url)
    '''
    Tool for submitting batches of json 'objects' to the REST api, one at a time.
    - inputFile is a json file, containing an array of objects, or a nested 'objects' array.
    '''
    
    authentication = None
    if(auth):
        authentication = auth.split(':')
        
    headers = {'content-type': 'application/json'}

    with open(inputFilePath, 'r') as fh:

        data = fh.read() #.decode(DEFAULT_ENCODING)
        objects = json.loads(data)
        
        if 'objects' in objects:
            objects = objects['objects']
        
        obj_url = url
        if( obj_url[len(obj_url)-1:] != '/'):
            obj_url = obj_url + '/'
        i = 0
        for obj in objects:
            if action == 'PUT':
                obj_url = url # + str(obj[id_key]) + '/'
                logger.info(str(('PUT (modify) the object: ', obj, 'url', obj_url)) )
                try:
                    # TODO: since not using the id in the url, may have to use the data = { objects: [obj]} wrapper as in patch
                    r = requests.put(obj_url, data=json.dumps(obj), auth=tuple(authentication), headers=headers)
                except Exception, e:
                    logger.error(str(('exception recorded while contacting server', e)))
                    raise e
                if(r.status_code != 202):
                    raise ApiError(obj_url,action, r) 

            elif action == 'PATCH':
                # NOTE: patch to the resource_uri value inside the obj
                obj_url = url # + str(obj[id_key]) + '/'
                logger.info(str(('PATCH (partial modify) the object: ', obj, 'url', obj_url)) )
                data = { "objects": [obj]} 
                try:
                    r = requests.patch(obj_url, data=json.dumps(data), auth=tuple(authentication), headers=headers)
                except Exception, e:
                    logger.error(str(('exception recorded while contacting server', e)))
                    raise e
                if(r.status_code != 202): 
                    raise ApiError(obj_url, action, r) 

            elif action == 'POST':
                logger.info(str(('POST (create) the object: ', obj, 'url', obj_url)) )
                try:
                    r = requests.post(obj_url, data=json.dumps(obj), auth=tuple(authentication), headers=headers)
                except Exception, e:
                    logger.error(str(('exception recorded while contacting server', e)))
                    raise e
                if(r.status_code != 201): 
                    raise ApiError(obj_url,action, r) 
            else:
                raise Exception('unknown action: "' + action + '"')
            
            logger.info(str(('result', r.status_code, r.text)))
            i += 1
        
        print action, i , ' objects to', url
            
                
parser = argparse.ArgumentParser(description='Import file')
parser.add_argument('-f', action='store', dest='inputFile',
                    metavar='FILE', required=True,
                    help='input file path')
parser.add_argument('-u', action='store', dest='url',
                    metavar='URL', required=True,
                    help='api url')
parser.add_argument('-i', action='store', dest='id_key',
                    metavar='ID_KEY', required=False,
                    help='specify the key for the ID in each record; note that this key must be specified in the API urls in prepend_urls ')
parser.add_argument('-auth', action='store', dest='auth',
                    metavar='AUTH', required=True,
                    help='api authentication to use <user:password>')
parser.add_argument('-a', '--action', dest='action', choices=['PUT', 'POST', 'PATCH'],
                help="choose PUT:for modifying records (each record must have an ID), POST: for creating records, PATCH: for modifying with partial update (supply the id_key)",
                required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='count',
                help="Increase verbosity (specify multiple times for more)")    

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

    print 'importing ', args.inputFile, ', to: ', args.url
    
    kwargs = {'auth':args.auth}
#    if args.id_key:
#        kwargs['id_key'] = args.id_key
    
    main(args.inputFile, args.url, args.action, **kwargs)
    
