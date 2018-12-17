# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

def migrate_extant_libraries(apps, schema_editor):

    Library = apps.get_model('db','Library')
    
    Library.objects.all().update(is_released=True)
    
    Library.objects.all()\
        .filter(library_type__in=['dos','annotation','nci'])\
        .update(is_archived=True)
        
class Migration(migrations.Migration):

    dependencies = [
        ('db', '0020_well_deprecation'),
    ]

    operations = [
        migrations.AddField(
            model_name='library',
            name='is_released',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='library',
            name='is_archived',
            field=models.BooleanField(default=False),
        ),
        migrations.RunPython(migrate_extant_libraries),
    ]
