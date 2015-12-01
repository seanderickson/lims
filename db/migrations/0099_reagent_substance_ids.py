# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import logging
import os

from django.db import migrations, models

from lims.base_settings import PROJECT_ROOT
from reports.utils.gray_codes import create_substance_id
from db.models import create_id

logger = logging.getLogger(__name__)

def create_substance_ids1(apps, schema_editor):
    # Retired 20151029
    # Legacy method - write out to a csv file, then use manual SQL to read in.
    # SQL method - takes 14 min on laptop
    # so 5x faster

    filename = os.path.join(PROJECT_ROOT, '..','new_substance_ids.csv')
    import csv
    
    _count = 0
    with open(filename, 'wb') as f:
        writer = csv.writer(f)
        # writer.writeheader(['reagent_id', 'substance_id'])
        for r_id in ( apps.get_model('db','Reagent').objects.all()
                .order_by('reagent_id') #.filter(reagent_id__gt=first)
                .values_list('reagent_id') ):
            reagent_id = r_id[0]
            # TODO: replace with db.models.create_id - so that the substance_id_seq is used
            writer.writerow([reagent_id, create_id()])
            _count += 1
            if _count %10000 == 0:
                logger.info('processed: %d' % _count )

    logger.info('wrote : %d substance_ids' % _count )

def create_substance_ids(apps, schema_editor):
    connection = schema_editor.connection
    cursor = connection.cursor()
    count = 0
    batch_size = 1000
    logger.info('begin create_substance_ids')
    cursor.execute(
        'create temp table reagent_sub_ids as ' 
            "select nextval('substance_id_seq'), reagent_id " 
            "from reagent;")
    cursor.execute(
        'create temp table reagent_sub_ids1 '
            '(reagent_id integer, substance_id text )')
    cursor1 = connection.cursor()
    sql1 = 'insert into reagent_sub_ids1 values (%s, %s)'
    cursor.execute('Select * from reagent_sub_ids')
    cursor.arraysize = batch_size
    try:
        while True:
            rows = cursor.fetchmany()
            if not rows:
                break
            newrows = []
            for row in rows:
                newrows.append([row[1], create_substance_id(row[0])])
                count += 1
            cursor1.executemany(sql1, newrows)
            if count % 100000 == 0:
                logger.info('processed %d rows', count)

        logger.info('processed %d rows, now update reagent table', count)
        cursor.execute(
            'update reagent ' 
                'set substance_id=a.substance_id ' 
                'from reagent_sub_ids1 a '
                'where a.reagent_id=reagent.reagent_id')
    
    except Exception as e:
        logger.exception('create substance ids, ex after %d rows',count)
        raise
    
    logger.info('create substance ids, done after %d rows',count)


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0098_db_post_migrate'),
    ]

    operations = [
        migrations.RunPython(create_substance_ids),
        migrations.AlterField(
            model_name='reagent',
            name='substance_id',
            field=models.CharField(unique=True, max_length=8),
            ),
    ]
