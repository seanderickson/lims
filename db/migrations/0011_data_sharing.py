# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging
import os
import re

from django.db import migrations, models

from db.support.data_converter import default_converter
import lims
from lims.base_settings import PROJECT_ROOT
import unicodecsv as csv

logger = logging.getLogger(__name__)


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0010_lab_affiliation'),
    ]

    operations = [
#         migrations.AddField(
#             model_name='screensaveruser',
#             name='rnai_data_sharing_level',
#             field=models.IntegerField(null=True),
#         ),
#         migrations.AddField(
#             model_name='screensaveruser',
#             name='sm_data_sharing_level',
#             field=models.IntegerField(null=True),
#         ),
    ]
