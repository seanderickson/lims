# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import logging
import os

from django.db import migrations, models


logger = logging.getLogger(__name__)

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0030_current_work'),
    ]

    operations = [
        
        migrations.CreateModel(
            name='WellQueryIndex',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('query', models.ForeignKey(to='db.CachedQuery')),
                ('well', models.ForeignKey(to='db.Well')),
            ],
            options={
                'db_table': 'well_query_index',
            },
        ),
        
        migrations.AlterField(
            model_name='attachedfile',
            name='type',
            field=models.TextField()),
        migrations.RemoveField(
            model_name='attachedfile',
            name='attached_file_type',
        ),
        
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
        migrations.RemoveField(
            model_name='serviceactivity',
            name='funding_support_link',
        ),
        # Note: funding_support_link is an alias for funding_support_id
        # migrations.RunSQL(
        #     'ALTER TABLE service_activity DROP COLUMN funding_support_id; '),
        migrations.DeleteModel(
            name='FundingSupport',
        ),
        
        migrations.DeleteModel(
            name='CellLine'),
            
        migrations.DeleteModel(
            name='ScreenStatusItem'),
            
        
        migrations.DeleteModel(
            name='TransfectionAgent'),

#         migrations.RemoveField(
#             model_name='ScreeningRoomUser',
#             name='last_notified_smua_checklist_item_event'),
# 
#         migrations.RemoveField(
#             model_name='ScreeningRoomUser',
#             name='last_notified_rnaiua_checklist_item_event'),

        migrations.DeleteModel(
            name='AttachedFileType',
        ),

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
       
#         migrations.AlterField(
#             model_name='cherrypickrequest',
#             name='volume_approved_by',
#             field=models.ForeignKey(
#                 related_name='approved_cherry_pick', to='db.ScreensaverUser', 
#                 null=True),
#         ),
#         migrations.AlterField(
#             model_name='cherrypickrequest',
#             name='requested_by',
#             field=models.ForeignKey(
#                 related_name='requested_cherry_pick', to='db.ScreensaverUser'),
#         ),
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
        
        # TODO: migration may not remove if already removed from the model        
        migrations.RemoveField(
            model_name='copy',
            name='well_concentration_dilution_factor',
        ),
        
        migrations.RemoveField(
            model_name='plate',
            name='plated_activity_id',
        ),
        migrations.RemoveField(
            model_name='plate',
            name='retired_activity_id',
        ),

#         migrations.RemoveField(
#             model_name='screenresult',
#             name='experimental_well_count',
#         ),


        migrations.DeleteModel(
            name='ChecklistItemEvent',
        ),
        
        migrations.DeleteModel(
            name='ChecklistItem',
        ),
        
        migrations.DeleteModel(
            name='AdministratorUser',
        ),
        migrations.DeleteModel(
            name='ScreeningRoomUser'),
                  
        migrations.DeleteModel(
            name='LabHead',
        ),
        
        migrations.DeleteModel(
            name='ScreensaverUserRole',
        ),

        # 20180920 - leaving this here for archival purposes
        # migrations.RunSQL('''
        #     DROP TABLE screening_room_user_facility_usage_role;
        # '''),

        
        # FIXME: not working on orchestra: moved to manual migration 0002
        # Keep here to convince makemigrations that this is done
        migrations.AlterUniqueTogether(
            name='assayplate',
            unique_together=set([('library_screening', 'plate', 'replicate_ordinal')]),
        ),       
        
    ]
