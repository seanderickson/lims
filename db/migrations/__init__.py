
import json
import logging

# from datetime import datetime, timedelta
import datetime
from django.db import migrations, models
from pytz import timezone
import pytz
from reports.models import ApiLog

logger = logging.getLogger(__name__)

# # unique offset for the logs in this migration to avoid collisions
# plate_vol_time_offset = 1111
# def create_log_time(input_date):
#     date_time = pytz.timezone('US/Eastern').localize(
#         datetime.datetime.combine(input_date, datetime.datetime.min.time()))
#     i = 0
#     date_time += datetime.timedelta(0,plate_vol_time_offset)
#     while date_time in times_seen:
#         i += 1
#         date_time += datetime.timedelta(0,i)
#     times_seen.add(date_time)
#     return date_time

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
    while _timekey in times_seen:
        i += 1
        logger.info('key: %s, adjust time: %s to %s', 
            key,
            date_time.isoformat(), (date_time + datetime.timedelta(0,i)))
        date_time += datetime.timedelta(0,i)
        _timekey = timekey % (key,date_time)
    times_seen.add(_timekey)
    return date_time


def _create_generic_log(activity):
    
    log = ApiLog()

    log.comment = activity.comments
    log.date_time = create_log_time(activity.date_of_activity)
    activity_user = activity.created_by
    # NOTE: if activity.date_performed < 2010-01-15, the "created_by" field
    # is empty
    if activity_user is None:
        activity_user = activity.performed_by
    log.username = activity_user.username
    if not log.username:
        log.username = '%s: %s' % (
            activity_user.screensaver_user_id,
            activity_user.email)
    
    log.user_id = activity_user.screensaver_user_id
    
    return log

def _child_log_from(parent_log):
    child_log = ApiLog()
    child_log.parent_log = parent_log
    child_log.username = parent_log.username
    child_log.user_id = parent_log.user_id
    child_log.date_time = parent_log.date_time
    child_log.api_action = parent_log.api_action
    child_log.comment = parent_log.comment
    return child_log
