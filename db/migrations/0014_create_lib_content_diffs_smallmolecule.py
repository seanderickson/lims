# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging
from django.db import migrations, models

logger = logging.getLogger(__name__)

def create_sm_diffs(apps,schema_editor):
    from db.support.library_content_migrator import Migrator
    Migrator().do_migration(apps, schema_editor, screen_type='small_molecule')

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0013_library_comments'),
    ]

    operations = [
        migrations.RunPython(create_sm_diffs),
    ]
