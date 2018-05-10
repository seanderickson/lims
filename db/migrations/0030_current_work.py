# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import logging
import os
import re

from django.db import migrations, models

from db.support.data_converter import default_converter
import lims


logger = logging.getLogger(__name__)


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0020_well_deprecation'),
    ]

    operations = [
    ]
