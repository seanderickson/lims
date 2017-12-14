# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
# from datetime import timedelta
# import datetime
import json
import logging
import os

from django.contrib.auth.models import User, UserManager
from django.db import migrations, models
import pytz

from db.support.data_converter import default_converter
from lims.base_settings import PROJECT_ROOT


logger = logging.getLogger(__name__)


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0003_db_migration_prep'),
        ('auth', '0001_initial')
    ]

    operations = [
        # NOTE: placing these schema migrations here for convenience (so that
        # they do not conflict with data migrations in 0003, 0004, 0007)
        migrations.AddField(
            model_name='screensaveruser',
            name='lab_head',
            field=models.ForeignKey(
                related_name='lab_members', to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='screensaveruser',
            name='lab_affiliation',
            field=models.ForeignKey(related_name='lab_heads', 
                to='db.LabAffiliation', null=True),
        ),
#         migrations.AddField(
#             model_name='screensaveruser',
#             name='rnai_data_sharing_level',
#             field=models.IntegerField(null=True),
#         ),
#         migrations.AddField(
#             model_name='screensaveruser',
#             name='sm_data_sharing_level',
#             field=models.IntegerField(null=True),
#         ),
        # 20170918; stashing here, to avoid sql pending trigger error in 0003
        migrations.RemoveField(
            model_name='screen',
            name='transfection_agent',
        ),
        migrations.RenameField(
            model_name='screen', 
            old_name='transfection_agent_text', 
            new_name='transfection_agent'
        ),
        
        
    ]
