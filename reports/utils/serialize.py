import csv
import logging

logger = logging.getLogger(__name__)


def from_csv(csvfile):
    '''
    NOTES: 
    - nested lists are denoted by brackets, i.e. '[]',
    - to escape use '\[...' (i.e. when embedding a regex expression)
    TODO: version 2 - read from a stream
    '''
    reader = csv.reader(csvfile)
    return from_csv_iterate(reader)
    
def from_csv_iterate(iterable):
    data_result = []
    i = 0
    keys = []
    list_keys = [] # cache
    for row in iterable:
        if i == 0:
            keys = row
        else:
            item = dict(zip(keys,row))
            for key in item.keys():
                val = item[key]
                if val and len(val)> 1:
                    if val[0] == '\\' and val[1] == '[':
                        # this could denote an escaped bracket, i.e. for a regex
                        item[key] = val[1:]
                    elif key in list_keys or val[0]=='[':
                        # due to the simplicity of the serializer, above, any 
                        # quoted string is a nested list
                        list_keys.append(key)
                        item[key] = val.strip('"[]').split(',')
            data_result.append(item)
        i += 1
    logger.debug('read in data, count: ' + str(len(data_result)) )   
    return data_result
