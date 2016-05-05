# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import logging
import os

from django.db import migrations, models

from lims.base_settings import PROJECT_ROOT
from reports.utils.gray_codes import create_substance_id
from db.models import create_id

logger = logging.getLogger(__name__)

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0020_current_work'),
    ]

    operations = [
        migrations.AlterField(
            model_name='attachedfile',
            name='type',
            field=models.TextField()),
        migrations.RemoveField(
            model_name='attachedfile',name='attached_file_type_id'),
        
#         migrations.RunSQL('ALTER TABLE screen DROP COLUMN cell_line_id ; '),
#         migrations.RunSQL('DROP TABLE screen_cell_line; '),
        migrations.RunSQL('ALTER TABLE screen DROP COLUMN transfection_agent_id ; '),
        migrations.RunSQL('DROP TABLE transfection_agent; '),
        migrations.RunSQL('DROP TABLE service_activity; '),
        migrations.RunSQL('DROP TABLE screen_funding_support_link; '),
        migrations.RunSQL('DROP TABLE funding_support; '),

#         migrations.DeleteModel('ScreeningRoomUser'),
#         migrations.DeleteModel('LabHead'),
    ]
