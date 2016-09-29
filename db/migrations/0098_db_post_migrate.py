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
            model_name='attachedfile',
            name='attached_file_type',
        ),

        # TODO Reinstate: after migration finished
        # migrations.RunSQL('ALTER TABLE screen DROP COLUMN cell_line_id ; '),
        # migrations.RunSQL('DROP TABLE screen_cell_line; '),
        migrations.RunSQL('ALTER TABLE screen DROP COLUMN transfection_agent_id ; '),
        migrations.RunSQL('DROP TABLE transfection_agent; '),
        migrations.RunSQL(
            'ALTER TABLE service_activity DROP COLUMN funding_support_id; '),
        migrations.RemoveField(
            model_name='screenfundingsupportlink',
            name='funding_support',
        ),
        migrations.RemoveField(
            model_name='screenfundingsupportlink',
            name='screen',
        ),
        migrations.DeleteModel(
            name='ScreenFundingSupportLink',
        ),

        # TODO: service_activity depends on funding support
        migrations.RunSQL('DROP TABLE funding_support; '),
        
        migrations.RemoveField(
            model_name='datacolumnderivedfromlink',
            name='derived_data_column',
        ),
        migrations.RemoveField(
            model_name='datacolumnderivedfromlink',
            name='derived_from_data_column',
        ),
        migrations.DeleteModel(
            name='DataColumnDerivedFromLink',
        ),

        # Operations already handled in migration 0002
        migrations.AlterField(
            model_name='screensaveruser',
            name='lab_head',
            field=models.ForeignKey(related_name='lab_member', to='db.ScreensaverUser', null=True),
        ),
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='volume_approved_by',
            field=models.ForeignKey(
                related_name='approved_cherry_pick', to='db.ScreensaverUser', null=True),
        ),
        
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='requested_by',
            field=models.ForeignKey(related_name='requested_cherry_pick', to='db.ScreensaverUser'),
        ),
        migrations.AlterField(
            model_name='screensaveruser',
            name='lab_head',
            field=models.ForeignKey(related_name='lab_member', blank=True, to='db.ScreensaverUser', null=True),
        ),
        migrations.RemoveField(
            model_name='well',
            name='latest_released_reagent',
        ),
        migrations.DeleteModel(
            name='AnnotationValue',
        ),
        migrations.DeleteModel(
            name='AnnotationType',
        ),
        migrations.DeleteModel(
            name='StudyReagentLink',
        ),

        
#         migrations.DeleteModel('ScreeningRoomUser'),
#         migrations.DeleteModel('LabHead'),
    ]
