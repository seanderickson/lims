
import datetime
import json
import logging

from django.db import migrations, models
from pytz import timezone
import pytz
from reports.models import ApiLog


# this is a workaround, because some activities have identical date_of_activity
times_seen = set()
# unique offset for the logs in this migration to avoid collisions
plate_vol_time_offset = 1111
def create_log_time(input_date):
    date_time = pytz.timezone('US/Eastern').localize(
        datetime.datetime.combine(input_date, datetime.datetime.min.time()))
    i = 0
    date_time += datetime.timedelta(0,plate_vol_time_offset)
    while date_time in times_seen:
        i += 1
        date_time += datetime.timedelta(0,i)
    times_seen.add(date_time)
    return date_time

def _create_generic_log(activity):
    
    log = ApiLog()

    log.comment = activity.comments
    log.date_time = create_log_time(activity.date_of_activity)
    log.username = activity.performed_by.username
    if not log.username:
        log.username = '%s: %s' % (
            activity.performed_by.screensaver_user_id,
            activity.performed_by.email)
    
    log.user_id = activity.performed_by.screensaver_user_id
    
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
