# Utilities for generating gray_code values for integers in a given base

# Gray code is a form of binary encoding where transitions between consecutive 
# numbers differ by only one bit.
# see http://rosettacode.org/wiki/Gray_code#Python

import string
import logging

logger = logging.getLogger(__name__)

def int2bin(n):
    'From positive integer to list of binary bits, msb at index 0'
    if n:
        bits = []
        while n:
            n,remainder = divmod(n, 2)
            bits.insert(0, remainder)
        return bits
    else: return [0]
    
def bin2int(bits):
    'From binary bits, msb at index 0 to integer'
    i = 0
    for bit in bits:
        i = i * 2 + bit
    return i    

def bin2gray(bits):
    return bits[:1] + [i ^ ishift for i, ishift in zip(bits[:-1], bits[1:])]

def gray2bin(bits):
    b = [bits[0]]
    for nextb in bits[1:]: b.append(b[-1] ^ nextb)
    return b

def xnor(a, b):
    if a==b:
        return 1
    else:
        return 0
    
def bin2igray(bits):
    if bits[len(bits)-1] == 0: # even
        return bits[:1] + [i ^ ishift for i, ishift in zip(bits[:-1], bits[1:])]
    else:    
        return bits[:1] + [xnor(i,ishift) for i, ishift in zip(bits[:-1], bits[1:])]

def igray2bin(bits):
    b = [bits[0]]
    if bits[len(bits)-1] == 0: # even
        for nextb in bits[1:]: b.append(b[-1] ^ nextb)
    else:
        for nextb in bits[1:]: b.append(xnor(b[-1], nextb))
        
    return b


digs = string.digits + string.uppercase
digs = digs.translate(None,'01OIS') # remove unwanted chars

# http://stackoverflow.com/questions/2267362/convert-integer-to-a-string-in-a-given-numeric-base-in-python
def int2base(x,symbols):
    '''
    @param x number to convert to base x where x = (len(symbols))
    @param symbols - non repeating char sequence
    '''
    base = len(symbols)
    if x < 0: sign = -1
    elif x==0: return symbols[x]
    else: sign = 1
    x *= sign
    digits = []
    while x:
        digits.append(symbols[x % base])
        x /= base
    if sign < 0:
        digits.append('-')
    digits.reverse()
    return ''.join(digits)

def int2gray_base(x, symbols):
    '''
    @param x number to convert to gray code in base x where x = (len(symbols))
    @param symbols - non repeating char sequence
    '''
    return int2base(bin2int(bin2gray(int2bin(x))),symbols)

def int2igray_base(x, symbols):
    '''
    @param x number to convert to gray code in base x where x = (len(symbols))
    @param symbols - non repeating char sequence
    '''
    return int2base(bin2int(bin2igray(int2bin(x))),symbols)

start = len(digs)**8 # start with a full 8 digit register       

def create_substance_id(_id):
    if(logger.isEnabledFor(logging.DEBUG)):
        logger.debug(str(('create substance id from', start, _id)))
    return int2gray_base(start + _id, digs);
    
if __name__ == '__main__':
#     start = len(digs)**7 # start with a full register
#     start = 32
#     start = 0
    size = 35
    base = len(digs)
    

    print 'digs', digs, 'len', len(digs), 'start', start
    _set = set()
    for i in range(100):
#         print ('input: %s,  base %s: substance_id %s' 
#             % ( i, base, create_substance_id(i)) )
        x = i + start
#         print('int:%2i -> bin:%12r -> igray:%12r -> bin:%12r -> int:%2i' %
#               ( x,
#                 int2bin(x),
#                 bin2igray(int2bin(x)),
#                 gray2bin(bin2gray(int2bin(x))),
#                 bin2int(igray2bin(bin2igray(int2bin(x))))
#               ))
#         print ('input: %s,  base %s: igray coded %s' 
#             % ( x, base, int2igray_base(x, digs)) )

#         val = int2igray_base(x, digs)
#         if val not in _set:
#             _set.add(val)
#         else:
#             raise Exception(str(('duplicate value', val, _set)))
#         if x%100000 == 0:
        print ('input: %s, base %s, \n_bin: %s,\nibin: %s, igray: %s' 
            % ( x, base, int2bin(x),bin2gray(int2bin(x)), int2gray_base(x, digs) ))

