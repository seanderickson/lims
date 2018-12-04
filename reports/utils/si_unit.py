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
      ['T', Decimal('1e12')],
      ['G', Decimal('1e9')],
      ['M', Decimal('1e6')],
      ['k', Decimal('1e3')],
      ['', Decimal(1)],
      ['m', Decimal('1e-3')],
      ['μ', Decimal('1e-6')],
      ['n', Decimal('1e-9')],
      ['p', Decimal('1e-12')]
      ]
    for symbol,val in siunits:
        if val <= abs(test_value):
            return (symbol,val)

    return None

def remove_exponent(d):
    '''Remove exponent and trailing zeros (see python decimal doc)
    '''
    return d.quantize(Decimal(1)) if d == d.to_integral() else d.normalize()


def round_decimal(raw_val, decimals=1, track_significance=False):
    '''
    Round to the specified decimals
    
    @param track_significance if True then remove trailing zeros
    
    @return Decimal representing the rounded raw_val as a decimal
    '''
    assert decimals >= 0, 'decimals must be >= 0'
    
    if raw_val is None: 
        return None
    
    val = Decimal(raw_val)

    logger.debug('convert_decimal: %r, %r, %r, %r',
        raw_val, decimals, track_significance)
    
    # Convert the decimal argument to a Decimal value
    decimals = Decimal('1e-%d'%int(decimals))
    # Quantize: Return a value equal to the first operand after rounding and 
    # having the exponent of the second operand. (see python doc)
    val = val.quantize(decimals, decimal.ROUND_HALF_UP)
        
    if track_significance is False:
        val = remove_exponent(val)

    return val

def convert_decimal(
    raw_val, default_unit='1e-6', decimals=1, multiplier=None, 
    track_significance=False):
    '''
    Convert raw_val as a Decimal by:
    - scaling to the default unit (e.g. 1e-6, convert 0.0000015 to 1.5)
    - rounding to the given decimals (decimal digits), 
    - optionally multiplying by a multiplier.
    
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
    
    @return a Decimal representation of the raw_val
    '''
    assert decimals >= 0, 'decimals must be >= 0'
    
    if raw_val is None:
        return None
    
    logger.debug('convert_decimal: %r, %r, %r, %r, %r',
        raw_val, default_unit, decimals, multiplier, track_significance)

    val = Decimal(raw_val)

    # - Get the scale (exponent) of the default unit
    #     DOC: "the adjusted exponent after shifting out the coefficient’s 
    #           rightmost digits until only the lead digit remain"
    # - Negate the scale for use with Decimal.scaleb()
    scale = -Decimal(str(default_unit)).adjusted()
    if multiplier is not None:
        # get the scale (exponent) of the multiplier
        multiplier = Decimal(str(multiplier)).adjusted()
        if multiplier != 0:
            scale = scale+multiplier

    # Convert the value using the determined scaling factor:
    if scale != 0:
        val = val.scaleb(scale)

    val = round_decimal(val, decimals=decimals, track_significance=track_significance)    
    
#     decimals = Decimal('1e-%d'%int(decimals))
#     val = val.quantize(decimals, decimal.ROUND_HALF_UP)
#         
#     def remove_exponent(d):
#         return d.quantize(Decimal(1)) if d == d.to_integral() else d.normalize()
#     if track_significance is False:
#         val = remove_exponent(val)
    
    return val
