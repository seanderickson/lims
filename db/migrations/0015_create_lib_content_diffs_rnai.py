# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging
from django.db import migrations, models
from reports.models import ApiLog

logger = logging.getLogger(__name__)

def create_rnai_diffs(apps,schema_editor):
    from db.support.library_content_migrator import Migrator
    Migrator().do_migration(apps, schema_editor, screen_type='rnai')

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0014_create_lib_content_diffs_smallmolecule'),
    ]

    operations = [
        migrations.RunPython(create_rnai_diffs),
    ]
