from __future__ import unicode_literals

import re

WELL_ID_PATTERN = re.compile(r'^(\d{1,5}):(([a-zA-Z]{1,2})(\d{1,2}))$')
WELL_NAME_PATTERN = re.compile(r'^([a-zA-Z]{1,2})(\d{1,2})$')
PLATE_PATTERN = re.compile(r'^(\d{1,5})$')
PLATE_RANGE_PATTERN = re.compile(r'^(\d{1,5})-(\d{1,5})$')
# Copy name pattern can be any string, starting with an alph char
# FIXME/TODO: constrain this pattern
COPY_NAME_PATTERN = re.compile(r'^[A-Za-z]+.*')
