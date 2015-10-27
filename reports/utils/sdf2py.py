from __future__ import unicode_literals

import cStringIO
import logging
import re

import six


# Imports an SDF/molfile
# Original version by Gabriel Barriz
logger = logging.getLogger(__name__)
# ---------------------------------------------------------------------------


VERBOSE = False,
ENCODING = u'utf8',
MOLDATAKEY = u'molfile'
COMMENTKEY = u'comment'
COMMENTTAG = u'comment'

# ---------------------------------------------------------------------------
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

    yield (COMMENTTAG, comments if comments else None)

def first_nonempty_line(preamble):
    for txt in [line.strip() for line in preamble.splitlines()]:
        if len(txt):
            return txt

def parse_sdf(data, _delimre=re.compile(ur'(?<=\n)\$\$\$\$')): #ur'(?<=\n)\$\$\$\$\n')):
    """
    for mol_record in sdf2py.parse_sdf(sdf_data_as_a_string):
        preamble = mol_record.pop(_params.MOLDATAKEY)
        title = first_nonempty_line(preamble) or 'WARNING: NO TITLE FOUND'
        print title
        for tag, value in mol_record.items():
            print (u'%s: %s' % (tag, value)).encode(_params.ENCODING)
        print
    """
    #     data = data.strip(u'\r\n')
    #     data = data.strip(u'\n')
    result = []
    if isinstance(data, six.string_types):
        data = data.strip()
        data = data.strip(u'$')

        mols = _delimre.split(data)
        # original
        #     return tuple(dict(parse_mol(mol)) for mol in mols)
    
        for mol in mols:
            x = dict(parse_mol(mol))
            result.append(x)
    else: # treat the data as an iterable
        
        buffer = cStringIO.StringIO()
        for line in data:
            if _delimre.match(line):
                x = dict(parse_mol(buffer.read()))
                result.append(x)
                buffer = cStringIO.StringIO()
            else:
                buffer.write(line)
    
    return tuple(result)
        

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
