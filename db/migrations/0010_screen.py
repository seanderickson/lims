# -*- coding: utf-8 -*-    
from __future__ import unicode_literals

import datetime
import json
import logging

from django.db import migrations, models, transaction
from django.db.models import F
from django.db.utils import IntegrityError
from pytz import timezone
import pytz

from db.migrations import create_log_time, times_seen
from reports.utils import default_converter
from reports.models import ApiLog
from reports.schema import DATE_FORMAT


logger = logging.getLogger(__name__)

DB_API_URI = '/db/api/v1'

# Default Admin User for log records (should be updated for other facilities)
USERNAME_ADMIN = 'sde_EDIT'
USER_ID_ADMIN = 7

def make_log(status_item):

    logger.debug('make status log: %r', status_item)
    collision_counter=1

    log = ApiLog()
    log.key = status_item['screen_facility_id']
    log.date_time = create_log_time(log.key,status_item['date']) 
    log.username = status_item['username']
    log.user_id = status_item['user_id']
    log.ref_resource_name = 'screen'
    log.uri = '/'.join([DB_API_URI, log.ref_resource_name, log.key])
    log.diffs = status_item.get('diffs')
    log.comment = status_item.get('comments')
    try:
        # check for log key (date_time) collisions; this shouldn't 
        # happen with the "create_log_time()", but, in case it does
        with transaction.atomic():
            log.save()
    except IntegrityError as e:
        q = ApiLog.objects.filter(
                ref_resource_name=log.ref_resource_name, key = log.key).order_by('-date_time')        
        if q.exists():    
            max_datetime = ( q.values_list('date_time', flat=True))[0]
        else:
            max_datetime = log.date_time
        logger.info('log time collision: %s, adjust log time from : %s to %s', 
            e, max_datetime.isoformat(), 
            (max_datetime + datetime.timedelta(0,collision_counter)))
        max_datetime += datetime.timedelta(0,collision_counter)
        times_seen.add(max_datetime)
        log.date_time = max_datetime
        collision_counter = collision_counter + 1

def migrate_screen_status(apps,schema_editor):
    '''
    Migrate the screen_status_item table to the "status" field of the screen
    object and create the needed ApiLog entries to record the history
    NOTE: manual migration 0002 must be run first:
    - adds an "id" field to screen_status_item
    - copies the latest status to the screen.status field    
    '''
    Screen = apps.get_model('db', 'Screen')
    ScreenStatusItem = apps.get_model('db','ScreenStatusItem')
        
    # Create a history log for all of the status's for each screen, 
    # and store the _current/latest_ status on the new screen.status field
    count=0
    for screen in Screen.objects.all().order_by('facility_id'):
        
        logger.info('process screen: %s', screen.facility_id)
        if screen.status:
            # Note: screen.status is set in manual/0002
            # Clean up the vocabulary used in the status_item table
            # NOTE: migration 0003 shows a newer way of generating the vocabs
            # - this is ok for these
            screen.status = default_converter(screen.status)
            screen.save()
            
        # Now scan both the screen.status_items and the screen.pin_transfer_activities
        # to concoct a status history table
        
        status_items= []
        
        for status in ScreenStatusItem.objects\
            .filter(screen=screen).order_by('status_date'):
            
            status_item = {
                'date': status.status_date,
                'status': default_converter(status.status),
                'username': USERNAME_ADMIN, # NOTE: no logs available in SS1
                'user_id': USER_ID_ADMIN,
                'screen_facility_id': str(screen.facility_id)
                }
            status_items.append(status_item)

        for s in Screen.objects\
            .filter(pin_transfer_admin_activity_id__isnull=False):
        
            activity = s.pin_transfer_admin_activity;

            status_item = {
                'date': activity.date_of_activity, 
                'status': 'transfer_approved',
                'username': activity.created_by.username,
                'user_id': activity.created_by.user_id,
                'screen_facility_id': str(screen.facility_id),
                'comment': activity.comments
                }
            status_items.append(status_item)
        
        # sort by date to get the right status history
        prev_item = None
        status_items = sorted(status_items, key=lambda x: x['date'])
        for status_item in status_items:
            diffs = {
                'status': [None, status_item['status']],
                'date': [None, status_item['date'].strftime(DATE_FORMAT)]
                }
            if prev_item:
                diffs['status'][0] = prev_item['status']
                diffs['date'][0] = prev_item['date'].strftime(DATE_FORMAT)

            status_item['diffs'] = diffs
            prev_item = status_item
            
        for status_item in status_items:
            make_log(status_item)
            count = count+1

    logger.info('updated: %d screen status entries', count)


def migrate_screen_project_phase(apps,schema_editor):
    '''
    project_phase is used to distinguish primary screens, follow up screens, 
    and studies.
    - migration:
    - populate screen.parent_screen
    - if project_phase != annotation, set the screen.study_type = null
    '''
    
    # 1. Identify and set the parent screens for follow up screens
    Screen = apps.get_model('db','Screen')
    for screen in (Screen.objects.all()
        .filter(project_id__isnull=False)
        .exclude(project_id__exact='')
        .exclude(project_id__exact=F('facility_id'))):
        logger.info('screen: %r find parent screen: %r', 
            screen.facility_id, screen.project_id)
        parent_screen = Screen.objects.get(facility_id=screen.project_id)
        screen.parent_screen = parent_screen
        logger.info(
            'screen: %s, parent_screen: %s', 
            screen.facility_id, parent_screen.facility_id)
        screen.save()
        logger.info('%r, %r', screen.facility_id, screen.parent_screen.facility_id)

#     # 2. verify that all "follow up" screens have a parent screen
#     query = ( Screen.objects.all()
#         .filter(project_phase__exact='follow_up_screen')
#         .exclude(parent_screen__isnull=True))
#     if query.exists():
#         raise Exception('not all follow_up_screens were converted: %r',
#             [x.facility_id for x in query.all()])
# 
    # 3. Set the "study_type" to null for all non-study screens
    query = ( Screen.objects.all().exclude(project_phase__icontains='annotation'))
    count = query.update(study_type=None)
    logger.info('study_type set to null for screens (count): %d', count)

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0008_usergroups'),
        ('reports', '0001_initial'), 
    ]

    operations = [
        migrations.RunPython(migrate_screen_status),
        migrations.RunPython(migrate_screen_project_phase),
    ]
