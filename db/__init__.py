from __future__ import unicode_literals

import re

WELL_ID_PATTERN = re.compile(r'^(\d+):?(([a-zA-Z]{1,2})(\d{1,2}))$')
COPYWELL_ID_PATTERN = re.compile(r'^(.*\/)?([A-Za-z]+[\w\- :]*)\/(\d+):?(([a-zA-Z]{1,2})(\d{1,2}))$')
WELL_NAME_PATTERN = re.compile(r'^([a-zA-Z]{1,2})(\d{1,2})$')
PLATE_PATTERN = re.compile(r'^(\d+)$')
PLATE_RANGE_PATTERN = re.compile(r'^(\d+)\s*-\s*(\d+)$$')
# Copy name pattern:
# - must start with an ascii alphabetic character
# - may contain spaces,dashes, and colons
COPY_NAME_PATTERN = re.compile(r'^[A-Za-z]+[\w\- :]*$')
