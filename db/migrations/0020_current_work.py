# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import logging
import os

from django.db import migrations, models

import lims
from lims.base_settings import PROJECT_ROOT
from db.support.data_converter import default_converter


logger = logging.getLogger(__name__)
import re


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0016_plate_volume'),
    ]

    operations = [
    ]
