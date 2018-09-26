# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import logging
import os

from django.db import migrations, models


logger = logging.getLogger(__name__)

def temp_migrate_breaker(apps,schema_editor):
    raise Exception('xxx stop migration xxx')

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0030_current_work'),
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

        migrations.RunSQL(
            'ALTER TABLE service_activity DROP COLUMN funding_support_id; '),
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
#         migrations.RunPython(temp_migrate_breaker),

# NOTE: 20180926 does not work, presumably because the fk for this has been altered in 0002
#         migrations.RemoveField(
#             model_name='screen',
#             name='pin_transfer_admin_activity',
#         ),
        

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
