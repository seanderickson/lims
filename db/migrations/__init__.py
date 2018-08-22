from __future__ import unicode_literals 
import json
import logging

import datetime
from django.db import migrations, models
from pytz import timezone
import pytz
from reports.models import ApiLog

logger = logging.getLogger(__name__)

# Hash of log times: 
# Workaround, because some activities have identical date_of_activity
times_seen = set()

def create_log_time(key,input_date):
    ''' Create a log timestamp that is unique for the key 
    ("plate_number", "screen_facility_id", etc):
    - Some activities in SS1 have identical date_of_activity
    '''
    
    date_time = pytz.timezone('US/Eastern').localize(
        datetime.datetime.combine(input_date, datetime.datetime.min.time()))
    i = 0
    
    timekey = '%r:%r'
    _timekey = timekey % (key,date_time)
    # if time has been seen, increment by a second, until time is unique
    # (to avoid api log key collisions)
    while _timekey in times_seen:
        i += 1
        logger.debug('%d, key: %s, adjust time: %s to %s', i, key,
            date_time.isoformat(), (date_time + datetime.timedelta(0,1)))
        date_time += datetime.timedelta(0,1)
        _timekey = timekey % (key,date_time)
    times_seen.add(_timekey)
    return date_time

def _child_log_from(parent_log):
    child_log = ApiLog()
    child_log.parent_log = parent_log
    child_log.username = parent_log.username
    child_log.user_id = parent_log.user_id
    child_log.date_time = parent_log.date_time
    child_log.api_action = parent_log.api_action
    child_log.comment = parent_log.comment
    return child_log
