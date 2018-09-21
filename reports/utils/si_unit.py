# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from decimal import Decimal
import decimal
import logging


logger = logging.getLogger(__name__)

def get_siunit_symbol(test_value):
    
    return get_siunit(test_value)[0]

def get_siunit(test_value):
    '''
    Return the best match (SI Unit symbol, value) for the given
    test_value, such that:
    test_value can be represented a number between 1 and 1000;
    (best_match_symbol_val)<=test_value<(next_higher_symbol_val)
    '''
    siunits = [
      ['T', 1e12],
      ['G', 1e9],
      ['M', 1e6],
      ['k', 1e3],
      ['', 1],
      ['m', 1e-3,],
      ['μ', 1e-6,],
      ['n', 1e-9 ],
      ['p', 1e-12 ]
      ]
    for symbol,val in siunits:
        if val <= abs(test_value):
            return (symbol,val)
    return None

def convert_decimal(
    raw_val, default_unit=1e-6, decimals=1, multiplier=None, 
    track_significance=False):
    '''
    Convert a decimal by scaling to the default unit and rounding to the 
    given decimals (decimal digits), optionally multiplying by a multiplier.
    
    @param default_unit adjust raw_val to the "default_unit" 
    (as defined in the "display_options")
    - e.g. if default_unit = 1e-6:
        adjust the raw_val = raw_val.scaleb(6)
    - e.g. if default_unit = 1e6:
        adjust the raw_val = raw_val.scaleb(-6)
    
    @param decimals digits of precision to apply, rounding using ROUND_HALF_UP
    @param multiplier (Note: only the exponent of the multiplier is used, so
    only powers of 10 may be used)
    @param track_significance if False, then trailing zeros (after rounding
    specified by "decimals") are dropped 
    '''
    assert decimals >= 0, 'decimals must be >= 0'
    logger.debug('convert_decimal: %r, %r, %r, %r, %r',
        raw_val, default_unit, decimals, multiplier, track_significance)
    # get the scale (exponent) of the default unit
    # negate the scale for use with Decimal.scaleb()
    scale = -Decimal(str(default_unit)).adjusted()
    if multiplier is not None:
        # get the scale (exponent) of the multiplier
        multiplier = Decimal(str(multiplier)).adjusted()
        if multiplier != 0:
            scale = scale+multiplier
    val = Decimal(raw_val)
    if scale != 0:
        val = val.scaleb(scale)
    decimals = Decimal('1e-%d'%int(decimals))
    val = val.quantize(decimals, decimal.ROUND_HALF_UP)
        
    def remove_exponent(d):
        return d.quantize(Decimal(1)) if d == d.to_integral() else d.normalize()
    if track_significance is False:
        val = remove_exponent(val)
    
    return val
