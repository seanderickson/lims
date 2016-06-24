# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging
from django.db import migrations, models
from db.support.data_converter import default_converter
from reports.models import ApiLog

logger = logging.getLogger(__name__)

def library_screen_vocab_conversion(apps,schema_editor):
    '''
    Performs the vocabulary conversion for the library and screen entities.
    - we are converting from the stored, display values, to a key value;
    this key will reference the entry for the respective vocab in the 
    reports_vocabularies table.  (Referencing by scope, key).
    Note: some conversions are also being done in SQL (manual scripts)
    '''
    # NOTE: migration 0003 shows a newer way of generating the vocabs
    # - this is ok for these
    
    library_vocabs = ['library_type', 'solvent', 'screen_type', 'screening_status']
    
    count = 0
    for obj in apps.get_model('db','Library').objects.all():
        for attr in library_vocabs:
            temp = getattr(obj, attr)
            if temp:
                temp2 = default_converter(temp)
                setattr(obj, attr, temp2)
                if(logger.isEnabledFor(logging.DEBUG)):
                    logger.debug(str(('convert library vocab', attr, temp,'to',temp2)))

        obj.save()
        count = count +1
    logger.info(str(('converted vocabs on', count, 'libraries')))
    
    
    screen_vocabs = ['screen_type', 'project_phase' ] 
    # Note: data_sharing_levels: using int values instead of names
    
    count = 0
    for obj in apps.get_model('db','Screen').objects.all():
        for attr in screen_vocabs:
            temp = getattr(obj, attr)
            if temp:
                temp2 = default_converter(temp)
                setattr(obj, attr, temp2)
                if(logger.isEnabledFor(logging.DEBUG)):
                    logger.debug(str(('convert screen vocab ', attr, temp,'to', temp2))) 
        obj.save()
        count += 1
    
    logger.info(str(('converted vocabs on', count, 'screens')))


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0006_screen_status'),
    ]

    operations = [
        migrations.RunPython(library_screen_vocab_conversion),
    ]
