# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

import logging

from reports.models import ApiLog
from db.migrations import create_log_time

logger = logging.getLogger(__name__)


def migrate_well_deprecations(apps, schema_editor):
    
    logger.info('migrate_well_deprecations')
    
    Well = apps.get_model('db','Well')

    def create_log(well, activity):
        log = ApiLog()
        log.ref_resource_name = 'well'
        log.key = well.well_id
        log.uri = '/'.join(['db/api/v1',log.ref_resource_name, log.key])
        log.username = activity.performed_by.username
        log.user_id = activity.performed_by.screensaver_user_id
        log.api_action = 'PATCH'
        if log.username is None:
            log.username = 'sde_EDIT'
        log.date_time = create_log_time(log.key,activity.date_of_activity) 
        log.diffs = {
            'is_deprecated': [False,True] }
        log.comment = activity.comments
        log.save()
        return log
    
    count = 0
    for well in Well.objects.all().filter(is_deprecated=True):
        a = well.deprecation_admin_activity.activity
        well.deprecation_reason = a.comments
        well.save()
        create_log(well, a)
        
        count += 1
        if count % 500 == 0:
            logger.info('converted %d wells', count)
    
    logger.info('done: converted %d wells', count)

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0019_raw_data_transformer'),
    ]

    operations = [
        migrations.AddField(
            model_name='well',
            name='deprecation_reason',
            field=models.TextField(null=True),
        ),

        migrations.RunPython(migrate_well_deprecations),
        migrations.RemoveField(
            model_name='well',
            name='deprecation_admin_activity'
        ),
        
    ]
