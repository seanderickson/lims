from __future__ import unicode_literals

import logging
import argparse
from itertools import combinations

logger = logging.getLogger(__name__)

def sum_bin(bin):
    return sum([p['size'] for p in bin])

def pack_bins_largest_first(capacity, package_array, packed_bins=None):
    '''
    Fit the packages in the package_array into the packed_bins.
    - Fit the largest packages first 
    '''
    
    package_array = sorted(package_array, reverse=True, key=lambda x: x['size'])
    logger.info('package array: %r, existing bins: %r', package_array, packed_bins)
    
    if packed_bins is None:
        packed_bins = []

    while package_array:
        
        placed = None
        for bin in packed_bins:
            space = capacity-sum_bin(bin)
            if space == 0: 
                continue
            placed = []
            logger.info('find package to fit in bin: %r, space: %r', bin, space)
            for i,package in enumerate(package_array):
                space = capacity-sum_bin(bin)
                if space >= package['size']:
                    bin.append(package)
                    placed.append(i)
                    logger.info('place package: %r, in bin: %r', package, bin)
            
            placed = sorted(placed,reverse=True)
            for i in placed:
                logger.info('remove placed package: %r', package_array[i])
                package_array.pop(i)
            if not package_array:
                break
        if not placed and package_array:
            logger.info('create new bin for: %r', package_array[0])
            packed_bins.append([package_array.pop(0)])
            
    return packed_bins
    
def pack_bins(capacity, package_array):
    '''
    Pack the packages in to bins having a capacity given.

    Pack all of the packages that fit first, then split packages that don't 
    fit, packing in available spaces, or creating new bins to fit them.
    
    Args:
        :param: package_array: a list of packages
        :type package_array: list[ dict[ name, size ] ]
        :param capacity: the size of the bins to pack
        :param split_bins: if True, then bin packing will include a final pass
        to reassign packages from unfilled bins to compress. Note that no 
        package will be assigned to more than 2 bins (plates)
    :return list[ dict[ name, size_packed ] ] 
    '''
    # First, pack all of the packages that fit, then run again to pack the 
    # packages that must be split
    
    packages_that_fit = [
        package for package in package_array if package['size']<=capacity]
    oversize_packages = [
        package for package in package_array if package['size']>capacity]
    
    packed_bins = pack_bins_largest_first(capacity, packages_that_fit)
    
    if oversize_packages:

        oversize_packages = sorted(oversize_packages, reverse=True, key=lambda x: x['size'])
        available_packed_bins = [bin for bin in packed_bins if capacity-sum_bin(bin)>0]
        
        # First, try to fit the oversized packages into (2) existing bins
        placed = []
        logger.info('Try to fit oversize_packages: %r', oversize_packages)
        for i,package in enumerate(oversize_packages):
            
            required_size = package['size']
            available_bins = [bin for bin in available_packed_bins
                if capacity-sum_bin(bin)>0]
            
            available_bin_pairs = combinations(available_bins, 2)
            
            for bin_pair in available_bin_pairs:
                
                capacity1 = capacity - sum_bin(bin_pair[0])
                capacity2 = capacity - sum_bin(bin_pair[1])
                
                if capacity1+capacity2 >= required_size:
                    bin_pair[0].append({
                        'name': package['name'], 'size': package['size']-capacity1 })
                    bin_pair[1].append({
                        'name': package['name'], 'size': package['size']-capacity2 })
                    placed.append(i)
                    
                    # Remove the bins from available bins,
                    # because there should only ever be one split package in a bin
                    available_packed_bins.remove(bin_pair[0])
                    available_packed_bins.remove(bin_pair[1])
                    
        placed = sorted(placed,reverse=True)
        for i in placed:
            logger.info('remove %d from %r', i, oversize_packages)
            oversize_packages.pop(i)
        
        # Second, split the oversize_packages, and fit the remainders
        remainder_packages = []
        for package in oversize_packages:
            
            full_bins = int(package['size']/capacity)
            for i in range(full_bins):
                packed_bins.append([{
                    'name': package['name'], 'size': capacity }])
            
            remainder_size = package['size'] % capacity
            logger.info('split package: %r into %d full bins, remainder_size: %d',
                package, full_bins, remainder_size)
            if remainder_size:
                remainder_packages.append({ 
                    'name': package['name'], 'size':remainder_size })
            
        if remainder_packages:
            logger.info('fit remainder packages: %r', remainder_packages)
            packed_bins = pack_bins_largest_first(
                capacity, remainder_packages, packed_bins=packed_bins)
    
#         if split_bins is True:
#             refilled_bins = []
#             unfilled_bins = []
#             # reassign packages from unfilled bins to compress. Note that no 
#             # package will be assigned to more than 2 bins (plates)
#             for bin in packed_bins:
#                 available = capacity-sum_bin(bin)
#                 if available > 0:
#                     bin['available': available]
#                     unfilled_bins.append(bin)
#             if len(unfilled_bins) > 1:
#                 # sort by available size, desc
#                 unfilled_bins = sorted(unfilled_bins, key=lambda x: x['available'])
#                 logger.info('unfilled bins: %r', unfilled_bins)
#                 while unfilled_bins:
#                     unfilled_bin = unfilled_bins.pop(0)
#                             
    # Sort within the bins by name, per convention
    packed_and_internally_sorted_bins = []
    for bin in packed_bins:
        packed_and_internally_sorted_bins.append(
            sorted(bin,key=lambda package: package['name']))
        
    return packed_and_internally_sorted_bins

def find_two_bins_for_package(capacity, package, packed_bins):
    two_bins = []
    required_space = package['size']
    # Sort by available space desc
    packed_bins = sorted(packed_bins, key=lambda x: capacity-sum_bin(x),reverse=True)
    for bin in packed_bins:
        available = capacity - sum_bin(bin)
        logger.info('bin: %r, available: %d', bin, available)
        if available < 1:
            continue
        if len(two_bins)>0:
            if required_space > available:
                continue
        two_bins.append(bin)
        required_space -= available
        logger.info('two bins: %r, reqd space: %d', two_bins, required_space)
        if required_space <= 0:
            break
        if len(two_bins) == 2:
            break
    if required_space > 0:
        return []
    if len(two_bins) != 2:
        return []
    return two_bins
    
parser = argparse.ArgumentParser(
    description='Pack the packages of varying sizes in bins having a fixed capacity')

parser.add_argument(
    '-c', '--capacity', required=True, type=int, 
    help='bin capacity')
parser.add_argument(
    '-p', '--package-array', required=True, type=int, nargs='+',
    help='comma separated list of package sizes')
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


    # for testing purposes, name of the package is its size
    packages = [ { 'name': p, 'size': p } for p in args.package_array ]
    packed_bins = pack_bins(args.capacity, packages)
    
    print 'packed_bins', packed_bins