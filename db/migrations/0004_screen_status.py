# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import json
import logging

from django.db import migrations, models
from django.utils import timezone

from db.support.data_converter import default_converter
from reports.models import ApiLog


logger = logging.getLogger(__name__)

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
        i=0
        if screen.status:
            screen.status = default_converter(screen.status)
            screen.save()
        # no scan the screen_status_items to recreate logs
        prev_item = None
        for status in ScreenStatusItem.objects.filter(screen=screen):
            log = ApiLog()
            log.date_time = timezone.now() 
            #status.status_date + datetime.timedelta(0,i) 
            # hack add 1 sec to avoid duplicate key error #timezone.now() 
            log.user_id = 1
            log.username = 'sde4'
            log.ref_resource_name = 'screen'
            log.key = screen.facility_id
            log.uri = '/db/api/v1/screen/' + screen.facility_id
            log.diff_keys = '["status","status_date"]'
            diffs = {}
            if prev_item:
                diffs['status'] = [prev_item.status, status.status]
                diffs['status_date'] = [unicode(prev_item.status_date), unicode(status.status_date)]
            else:
                diffs['status'] = [None, status.status]
                diffs['status_date'] = [None, unicode(status.status_date)]
            log.diffs = json.dumps(diffs)
            logger.debug(str(( 'create log: ' , j, log)))
            log.save()
             
            prev_item = status
            i = i + 1
            j = j + 1
    
    logger.info(str(( 'updated ', j, 'statuses')))


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0003_db_migration_prep'),
    ]

    operations = [
        migrations.RunPython(migrate_screen_status),
    ]
