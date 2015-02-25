import csv
import logging

logger = logging.getLogger(__name__)


def from_csv(csvfile, list_delimiter=';', list_keys=None):
    '''
    @param list_keys overrides nested list eval for column keys; no brackets '[]' are 
        needed to denote these columns as list columns - however, to comply with 
        the csv standard, they still have to be quoted (if list_delimiter=csv_delimiter)
    NOTES: 
    - nested lists are denoted by brackets, i.e. '[]',
    - to escape use '\[...' (i.e. when embedding a regex expression)
    TODO: version 2 - read from a stream
    '''
    reader = csv.reader(csvfile)
    return from_csv_iterate(reader, list_delimiter=list_delimiter, list_keys=list_keys)
    
def from_csv_iterate(iterable, list_delimiter=';', list_keys=None):
    list_keys = list_keys or []
    data_result = []
    i = 0
    keys = []
    list_keys = list(list_keys) 
    for row in iterable:
        if i == 0:
            keys = [x for x in row]
        else:
            item = dict(zip(keys,row))
            for key in item.keys():
                val = item[key]
#                 logger.info(str(('val', key, val)))
                if val and len(val)> 1:
                    if val[0] == '\\' and val[1] == '[':
                        # this could denote an escaped bracket, i.e. for a regex
                        item[key] = val[1:]
                    elif key in list_keys or val[0]=='[':
                        # due to the simplicity of the serializer, above, any 
                        # quoted string is a nested list
                        list_keys.append(key)
                        item[key] = [x.strip() for x in val.strip('"[]').split(list_delimiter)]
#                         logger.info(str(('item[key]',item[key], type(item[key]))))
            data_result.append(item)
        i += 1
    logger.debug('read in data, count: ' + str(len(data_result)) )   
    return data_result

