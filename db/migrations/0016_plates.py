
# -*- coding: utf-8 -*-
from __future__ import unicode_literals


from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0015_create_lib_content_diffs_rnai'),
    ]

    # Post-manual migration updates
    operations = [
        migrations.RemoveField(model_name='plate',name='min_molar_concentration'),
        migrations.RemoveField(model_name='plate',name='max_molar_concentration'),
        migrations.RemoveField(model_name='plate',name='min_mg_ml_concentration'),
        migrations.RemoveField(model_name='plate',name='max_mg_ml_concentration'),
        migrations.RemoveField(model_name='plate',name='primary_well_molar_concentration'),
        migrations.RemoveField(model_name='plate',name='primary_well_mg_ml_concentration'),

        migrations.RemoveField(model_name='copy',name='min_molar_concentration'),
        migrations.RemoveField(model_name='copy',name='max_molar_concentration'),
        migrations.RemoveField(model_name='copy',name='min_mg_ml_concentration'),
        migrations.RemoveField(model_name='copy',name='max_mg_ml_concentration'),
        migrations.RemoveField(model_name='copy',name='primary_well_molar_concentration'),
        migrations.RemoveField(model_name='copy',name='primary_well_mg_ml_concentration'),
        
        # These fields are replaced with sql queries in the api
        migrations.RemoveField(
            model_name='copy',
            name='primary_plate_status',
        ),
        migrations.RemoveField(
            model_name='copy',
            name='primary_plate_location_id',
        ),
        migrations.RemoveField(
            model_name='copy',
            name='plate_locations_count',
        ),
        migrations.RemoveField(
            model_name='copy',
            name='plates_available',
        ),
        
        
        migrations.AddField(
            model_name='labcherrypick',
            name='is_manually_selected',
            field=models.NullBooleanField(),
        ),
#         migrations.AddField(
#             model_name='plate',
#             name='key',
#             field=models.TextField(null=True),
#         ),
        
    ]
    
