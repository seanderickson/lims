# -*- coding: utf-8 -*-    
from __future__ import unicode_literals

from datetime import datetime, time, date, timedelta
import json
import logging

from django.db import migrations, models, transaction
from django.db.utils import IntegrityError
from pytz import timezone
import pytz

from db.support.data_converter import default_converter
from reports.models import ApiLog


logger = logging.getLogger(__name__)

DB_API_URI = '/db/api/v1'

times_seen = set()
def create_log_time(key,input_date):
    date_time = pytz.timezone('US/Eastern').localize(
        datetime.combine(input_date, datetime.min.time()))
    i = 0
    
    timekey = '%r:%r'
    _timekey = timekey % (key,date_time)
    while _timekey in times_seen:
        i += 1
        logger.info('key: %s, adjust time: %s to %s', 
            key,
            date_time.isoformat(), (date_time + timedelta(0,i)))
        date_time += timedelta(0,i)
        _timekey = timekey % (key,date_time)
    times_seen.add(_timekey)
    return date_time

def make_log(
    apps, input_date, ref_resource_name, key, 
    diffs=None, comment=None, username=None):

    apilog_model = apps.get_model('reports','apilog')
    collision_counter=0

    if diffs is None:
        diffs = {}
    
    log = ApiLog()
    log.date_time = create_log_time(key,input_date) 
    log.user_id = 1
    log.username = username
    log.ref_resource_name = ref_resource_name
    log.key = key
    log.uri = '/'.join([DB_API_URI, ref_resource_name, key])
    log.diffs = diffs
    log.comment = comment
    try:
        # check for log key (date_time) collisions; this shouldn't 
        # happen with the "create_log_time()", but, in case it does
        with transaction.atomic():
            log.save()
    except IntegrityError as e:
        q = apilog_model.objects.filter(
                ref_resource_name=ref_resource_name,
                key = key).order_by('-date_time')        
        if q.exists():    
            max_datetime = ( q.values_list('date_time', flat=True))[0]
        else:
            max_datetime = log.date_time
        logger.info('log time collision: %s, adjust log time from : %s to %s', 
            e, max_datetime.isoformat(), 
            (max_datetime + timedelta(0,collision_counter)))
        max_datetime += timedelta(0,collision_counter)
        times_seen.add(max_datetime)
        log.date_time = max_datetime
        collision_counter = collision_counter + 1
        
    return log

def migrate_pin_transfer_approval(apps, schema_editor):
    '''
    Migrate the pin_transfer_approval activities to the
    "pin_transfer_approved_by" field, and recored the date and comment
    '''
    Screen = apps.get_model('db', 'Screen')
    AdministrativeActivityModel = apps.get_model('db', 'AdministrativeActivity')
    count = 0
    for s in ( Screen.objects.all()
        .filter(pin_transfer_admin_activity__isnull=False)):
        activity = s.pin_transfer_admin_activity;
        
        logger.info('migrate pin transfer activity for screen: %r', s)
        logger.debug('process pin transfer activity: %r', activity)
        
        # Create an ApiLog - note, no de-duplication needed; for the current 
        # database, the pin transfer approval has only been set once
        diffs = {
            'pin_transfer_approved_by_username': 
                [None, activity.performed_by.username],
            'pin_transfer_date_approved': 
                [None, activity.date_of_activity.strftime("%Y-%m-%d")],
            'pin_transfer_comments': [None, activity.comments],
        }
        
        log = make_log(
            apps,
            activity.date_of_activity, 'screen', 
            s.facility_id, 
            diffs=diffs, comment=activity.comments,
            username=activity.created_by.username)
        count = count + 1
    
    logger.info('migrated %d pin_transfer_admin_activity logs', count)

def migrate_screen_status(apps,schema_editor):
    '''
    Migrate the screen_status_item table to the "status" field of the screen
    object and create the needed ApiLog entries to record the history
    NOTE: manual migration 0002 must be run first:
    - adds an "id" field to screen_status_item
    - copies the latest status to the screen.status field    
    '''
        
    # Create a history log for all of the status's for each screen, 
    # and store the _current/latest_ status on the new screen.status field
    count=0
    ScreenStatusItem = apps.get_model('db','ScreenStatusItem')
    for screen in (
            apps.get_model('db','Screen').objects.all().order_by('facility_id')):
        
        logger.info('process screen: %s', screen.facility_id)
        if screen.status:
            # Clean up the vocabulary used in the status_item table
            # NOTE: migration 0003 shows a newer way of generating the vocabs
            # - this is ok for these
            screen.status = default_converter(screen.status)
            screen.save()
        # now scan the screen_status_items to recreate logs
        prev_item = None
        for status in ( ScreenStatusItem.objects.filter(screen=screen)
                .order_by('status_date')):
            prev_status = default_converter(prev_item.status)
            new_status = default_converter(status.status)
            diffs = {}
            if prev_item:
                diffs['status'] = [prev_status, new_status]
            else:
                diffs['status'] = [None, new_status]

            log = make_log(
                apps,status.status_date, 'screen', screen.facility_id, 
                diffs, username='sde')
            logger.debug('created log: %d: %r', count, log)
            
            prev_item = status
            count = count + 1
    
    logger.info('updated: %d screen status entries', count)


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0005_usergroups'),
    ]

    operations = [
        migrations.RunPython(migrate_screen_status),
        migrations.RunPython(migrate_pin_transfer_approval),

    ]
