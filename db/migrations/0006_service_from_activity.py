# -*- coding: utf-8 -*-    
from __future__ import unicode_literals

from django.db import migrations, models
import django.utils.timezone

import logging
logger = logging.getLogger(__name__)

# from datetime import datetime, time, date, timedelta
from pytz import timezone
import pytz
# times_seen = set()
# def create_log_time(input_date):
#     date_time = pytz.timezone('US/Eastern').localize(
#         datetime.combine(input_date, datetime.min.time()))
#     i = 0
#     while date_time in times_seen:
#         i += 1
#         date_time += timedelta(0,i)
#         logger.info('adjust time: %s', date_time.isoformat())
#     times_seen.add(date_time)
#     return date_time

def migrate_service_lab_activities(apps,schema_editor):
     
    ServiceModel = apps.get_model('db', 'Service')
    
    def convert_date_time(date_time):
        if date_time:
            return pytz.timezone('US/Eastern').localize(date_time)    
        else:
            return None
        
    def set_activity_fields(a,s):
        s.date_time_created = convert_date_time(a.date_created)
        s.date_performed = a.date_of_activity
        s.created_by = a.created_by
        s.performed_by = a.performed_by
        s.comments = a.comments
         
        s.date_loaded = convert_date_time(a.date_loaded)
        s.date_publicly_available = convert_date_time(a.date_publicly_available)
    
    service_count = 0 
    for sa in apps.get_model('db','ServiceActivity').objects.all():
  
        sm = ServiceModel()
          
        set_activity_fields(sa.activity, sm)
         
        sm.type = sa.service_activity_type
        sm.screen = sa.serviced_screen
        sm.user = sa.serviced_user
        sm.funding_support = sa.funding_support
         
        sm.save()
        service_count += 1
    
    la_count = 0    
    for la in apps.get_model('db','LabActivity').objects.all():
         
        sm = ServiceModel()
 
        set_activity_fields(la.activity, sm)
         
        if hasattr(la,'screening'):
            if hasattr(la.screening, 'libraryscreening'):
                if la.screening.libraryscreening.is_for_external_library_plates:
                    sm.type = 'external_library_screening'
                else:
                    sm.type = 'library_screening'
            else:
                sm.type = 'cherry_pick_screening'
        else:
            sm.type = 'cherry_pick_liquid_transfer'
         
        sm.comments = None
        sm.screen = la.screen
        sm.user = None
        sm.funding_support = None
         
        sm.save()
        la_count += 1

    logger.info('migrated %d service and %d lab activities', service_count, la_count)
class Migration(migrations.Migration):

    dependencies = [
        ('db', '0005_screen_status'),
    ]

    operations = [
        migrations.CreateModel(
            name='Service',
            fields=[
                ('id', 
                    models.AutoField(verbose_name='ID', serialize=False, 
                        auto_created=True, primary_key=True)),
                ('date_time_created', 
                    models.DateTimeField(default=django.utils.timezone.now)),
                ('date_performed', models.DateField()),
                ('type', models.TextField()),
                ('comments', models.TextField(null=True, blank=True)),
                ('funding_support', models.TextField(null=True)),
                ('date_loaded', 
                    models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', 
                    models.DateTimeField(null=True, blank=True)),
                ('created_by', 
                    models.ForeignKey(related_name='services_created', 
                        blank=True, to='db.ScreensaverUser', null=True)),
                ('performed_by', 
                    models.ForeignKey(related_name='services_performed', 
                        to='db.ScreensaverUser')),
                ('screen', models.ForeignKey(to='db.Screen', null=True)),
                ('user', models.ForeignKey(to='db.ScreensaverUser', null=True)),
            ],
            options={
                'db_table': 'service',
            },
        ),
        migrations.RunPython(migrate_service_lab_activities),
    ]
