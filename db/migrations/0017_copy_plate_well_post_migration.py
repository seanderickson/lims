# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models



class Migration(migrations.Migration):
    # Post-migration after 0016:
    # remove old concentration fields from copy and plate
    
    dependencies = [
        ('db', '0016_plates'),
    ]

    
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
