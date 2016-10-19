# -*- coding: utf-8 -*-    
from __future__ import unicode_literals

import json
import logging

from django.db import migrations, models, transaction
from django.utils import timezone
from django.db.utils import IntegrityError

from db.support.data_converter import default_converter
from reports.models import ApiLog

from datetime import datetime, time, date, timedelta
from pytz import timezone
import pytz

logger = logging.getLogger(__name__)

times_seen = set()
def create_log_time(input_date):
    date_time = pytz.timezone('US/Eastern').localize(
        datetime.combine(input_date, datetime.min.time()))
    i = 0
    while date_time in times_seen:
        i += 1
        date_time += timedelta(0,i)
        logger.debug('adjust time: %s', date_time.isoformat())
    times_seen.add(date_time)
    return date_time

def migrate_screen_status(apps,schema_editor):
    '''
    Migrate the screen_status_item table to the "status" field of the screen
    object and create the needed ApiLog entries to record the history
    NOTE: manual migration 0002 must be run first:
    - adds an "id" field to screen_status_item
    - copies the latest status to the screen.status field    
    '''
        
    # first, clean up the vocabulary used in the status_item table
    # NOTE: migration 0003 shows a newer way of generating the vocabs
    # - this is ok for these
    
    for obj in apps.get_model('db','ScreenStatusItem').objects.all():
        attr = 'status'
        temp = getattr(obj, attr)
        if temp:
            temp2 = default_converter(temp)
            setattr(obj, attr, temp2)
            logger.debug(str(( 'screen', attr, temp,temp2))) 
        obj.save()
        
    # Now create a history log for all of the status's for each screen, 
    # and store the _current/latest_ status on the new screen.status field
    j=0
    ScreenStatusItem = apps.get_model('db','ScreenStatusItem')
    for screen in apps.get_model('db','Screen').objects.all():
        logger.info('process screen: %s', screen.facility_id)
        i=0
        if screen.status:
            screen.status = default_converter(screen.status)
            screen.save()
        # now scan the screen_status_items to recreate logs
        prev_item = None
        for status in ( ScreenStatusItem.objects.filter(screen=screen)
                .order_by('status_date')):
            log = ApiLog()
            log.date_time = create_log_time(status.status_date) 
            log.user_id = 1
            log.username = ''
            log.ref_resource_name = 'screen'
            log.key = screen.facility_id
            log.uri = '/db/api/v1/screen/' + screen.facility_id
            log.diffs = {}
            if prev_item:
                log.diffs['status'] = [prev_item.status, status.status]
            else:
                log.diffs['status'] = [None, status.status]
            logger.info('create log: %s: %r' , j, log)
            try:
                # use a nested atomic block to delimit rollback (the entire
                # migration is is the outer atomic block)
                with transaction.atomic():
                    log.save()
            except IntegrityError as e:
                apilog_model = apps.get_model('reports','apilog')
                max_datetime = ( 
                    apilog_model.objects.filter(
                        ref_resource_name='screen',
                        key = log.key)
                    .order_by('-date_time')
                    .values_list('date_time', flat=True))[0]
                max_datetime += timedelta(0,i)
                times_seen.add(max_datetime)
                logger.info('new log time: %s', max_datetime.isoformat())
                log.date_time = max_datetime

            prev_item = status
            i = i + 1
            j = j + 1
    
    logger.info(str(( 'updated ', j, 'statuses')))


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0005_usergroups'),
    ]

    operations = [
        migrations.RunPython(migrate_screen_status),
    ]
