
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
        
    ]
    
