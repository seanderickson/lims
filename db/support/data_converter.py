import re
# MIGRATE "OLD" SS1 fields to new fields

def convert_well_data(
        data,
        value_converters={
            'library_well_type': lambda x: default_converter(x) },
        field_mapping = {
            'plate': 'plate_number',
            'well': 'well_name',
            'well_type': 'library_well_type',
            'facility_reagent_id': 'facility_id',
            'vendor':'vendor_name',
            'vendor_reagent_id': 'vendor_identifier',
        } ):
    '''
    MIGRATE "OLD" SS1 fields to new fields
    @param data dict of data, converted in place
    '''
    for oldname,newname in field_mapping.items():
        if newname not in data:
            data[newname] = next(
                (val for (key,val) in data.items() if key.lower()==oldname),
                [None])
            if newname in value_converters:
                data[newname] = value_converters[newname](data[newname])
    return data

DEFAULT_CONVERTER=re.compile(r'[\W]+')
def default_converter(original_text):
        temp = DEFAULT_CONVERTER.sub(' ', original_text)
        return '_'.join(temp.lower().split())      
  