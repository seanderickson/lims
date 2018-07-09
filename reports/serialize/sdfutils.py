from __future__ import unicode_literals

import cStringIO
import logging
import re
import six
from reports.serialize import INPUT_FILE_DESERIALIZE_LINE_NUMBER_KEY

logger = logging.getLogger(__name__)

VERBOSE = False,
ENCODING = u'utf8',
MOLDATAKEY = u'molfile'
COMMENTKEY = u'comment'
COMMENTTAG = u'comment'



_moldata_re = re.compile(ur'M\s+END')
_dos2unix = re.compile(ur'\r\n')
_dos_unix_le = re.compile(ur'[\r\n]{1,2}')

def parse_mol(data, _delimre=re.compile(ur'(>\s+<[^>]+>[^\r^\n]*[\r\n]*)')):

    #original: def parse_mol(data, _delimre=re.compile(ur'((?<=\n)>\s+<[^>]+>[^\S\n]*\n)')): 
    #dos2unix: def parse_mol(data, _delimre=re.compile(ur'((?<=\n)>\s+<[^>]+>[^\r^\n]*[\r\n]*)')):
    # also doesn't work with all line endings in the clean_data_small_molecule.sdf file
    
    # TODO: review usage of the look-behind regex which should eliminate strip()
    parts = _delimre.split(data.strip())
    
    if not parts[0]:
        # if the separator matched at the start of the string, then there is an
        # empty value in the first position
        parts.pop(0)
    
    if _moldata_re.search(parts[0]):
        yield (MOLDATAKEY, _dos2unix.sub(ur'\n',parts.pop(0).strip()))
    
    comments = {}
        
    last_comment = None
    for tag, val in zip(*[iter(parts)] * 2):
        key = tag[tag.find(u'<')+1:tag.rfind(u'>')]
        v = val.strip()

        v = _dos_unix_le.split(v)
        if len(v) == 1: v = v[0]
        if len(v) == 0:
            v = None
        
        # TODO: review spec for comments:
        # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
        if key == COMMENTKEY:
            last_comment = v
            continue

        if last_comment:
            comments[key] = last_comment
            last_comment = None
        yield key, v

    assert not last_comment
    
    if comments:
        yield (COMMENTTAG, comments if comments else None)

def first_nonempty_line(preamble):
    for txt in [line.strip() for line in preamble.splitlines()]:
        if len(txt):
            return txt
                                        
def parse_sdf(data, _delimre=re.compile(ur'(?<=\n)\$\$\$\$')):
    """
    for mol_record in sdf2py.parse_sdf(sdf_data_as_a_string):
        preamble = mol_record.pop(_params.MOLDATAKEY)
        title = first_nonempty_line(preamble) or 'WARNING: NO TITLE FOUND'
        print title
        for tag, value in mol_record.items():
            print (u'%s: %s' % (tag, value)).encode(_params.ENCODING)
        print
        
    # TODO: implement this as a generator
    """
    result = []
    if isinstance(data, six.string_types):
        data = data.strip()
        data = data.strip(u'$')
        
        cumulative_lines = 1

        mols = _delimre.split(data)
        for mol in mols:
            mol_lines = len(_dos_unix_le.split(mol))
            x = dict(parse_mol(mol))
            x[INPUT_FILE_DESERIALIZE_LINE_NUMBER_KEY] = cumulative_lines
            cumulative_lines += mol_lines-1
            
            result.append(x)
    else: # treat the data as an iterable
        buffer = cStringIO.StringIO()
        linecount = 0
        record_line = 0
        for line in data:
            if _delimre.match(line):
                x = dict(parse_mol(buffer.read()))
                x[INPUT_FILE_DESERIALIZE_LINE_NUMBER_KEY] = record_line
                record_line = 0
                result.append(x)
                buffer = cStringIO.StringIO()
            else:
                if record_line == 0:
                    record_line = linecount
                buffer.write(line)
            linecount += 1
    return tuple(result)


def to_sdf(data,output):
    '''
    @param data an iterable containing dict objects to be written
    @param output a file like object
    '''
    
    if isinstance(data,dict):
        data = [data]
            
    for d in data:
        
        if d.get(MOLDATAKEY, None):
            output.write(str(d[MOLDATAKEY]))
            output.write('\n') 
            # because we've not copied the data, don't delete it
            # future optimize: implement data as iterable
            #                 del d[MOLDATAKEY]
        for k,v in d.items():
            if k == MOLDATAKEY: 
                continue
            output.write('> <%s>\n' % k)
            # according to 
            # http://download.accelrys.com/freeware/ctfile-formats/ctfile-formats.zip
            # "only one blank line should terminate a data item"
            if v is not None:
                # find lists, but not strings (or dicts)
                # Note: a dict here will be non-standard; probably an error 
                # report, so just stringify dicts as is.
                if not hasattr(v, "strip") and isinstance(v, (list,tuple)): 
                    for x in v:
                        # DB should be UTF-8, so this should not be necessary,
                        # however, it appears we have legacy non-utf data in 
                        # some tables (i.e. small_molecule_compound_name 193090
                        output.write(unicode.encode(x,'utf-8'))
#                             output.write(str(x))
                        output.write('\n')
                else:
                    output.write(str(v))
                    output.write('\n')

            output.write('\n')
        output.write('$$$$\n')

        

if __name__ == '__main__':
    import sys
    for p in sys.argv[1:]:
        print '%s:' % p
        with open(p) as fh:
            data = fh.read().decode(ENCODING)

        for mol_record in parse_sdf(data):
            preamble = mol_record.pop(MOLDATAKEY)
            title = first_nonempty_line(preamble) or 'WARNING: NO TITLE FOUND'
            print title
            for tag, value in mol_record.items():
                print (u'%s: %s' % (tag, value)).encode(ENCODING)
            print
