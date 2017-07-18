# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models

# def create_lab_affiliation_vocab(apps):
#  
#     # populate the title field, change the name field to a key
#     replace_phrases = [
#         ['harvard medical school', 'hms'],
#         ['harvard university', 'harvard'],
#         ['European Molecular Biology Laboratory', 'embl'],
#         ['[embl]',''],
#         ['Dana Farber Cancer Institute', 'dfci'],
#         ['University of California', 'uc'],
#         ['University of Massachusetts', 'umass'],
#         ['Institute of Chemistry and Cell Biology', 'iccb'],
#         ['Beth Israel Deaconess Medical Center', 'bidmc'],
#         ['Tufts University', 'tufts'],
#         ['University of California', 'UC'],
#         ['University of Massachusetts', 'umass'],
#         ['[NYU]', ''],
#         ['the',''],
#         ['of',''],
#         ['in',''],
#         ["women's", 'womens'],
#         ["children's", 'childrens']
#     ]
#     replace_phrases = [[re.compile(r'\b%s\b' % x, re.IGNORECASE),y] 
#         for [x,y] in replace_phrases ]
#     vocab_file = os.path.join(
#         lims.settings.PROJECT_ROOT, '..',
#         'db', 'static', 'api_init', 'vocabulary_lab_affiliations_data.csv')
#     logger.info('write vocabularies to %s' % vocab_file)
# #     resource_uri = 'vocabulary/%s/%s/'
#     with open(vocab_file, 'a+') as _file:
#          
#         header = ['key', 'scope', 'ordinal', 'title', 'comment'] 
#         reader = csv.reader(_file)
#         vocab_writer = csv.writer(_file)
#         defined_vocabs = {}
#         try:
#             header = reader.next()
#             logger.info('read header: %s', header)
#             for row in reader:
#                 defined_vocabs[row[4]] = dict(zip(header,row))
#         except StopIteration as e:
#             logger.info('no entries in %s, writing a new file',vocab_file)
#             vocab_writer.writerow(header)
#  
#         scope = 'labaffiliation.category'
#         default_order = ['hms','hms_affiliated_hospital','hsph',
#             'harvard_fas','broad_icg','other']
#         for i,la in enumerate(apps.get_model('db', 'LabAffiliation')
#                 .objects.all().distinct('affiliation_category')):
#             if la.affiliation_category not in defined_vocabs:
#                 
#                 key = default_converter(la.affiliation_category)
#                 if key in default_order:
#                     ordinal = default_order.index(key)
#                 else:
#                     ordinal = i + len(default_order)
#                 title = la.affiliation_category
#                 row = [key, scope, ordinal, title, la.affiliation_category ]
#                 vocab_writer.writerow(row)
#                 defined_vocabs[la.affiliation_category] = dict(zip(header,row))
#                 logger.debug('created %s', row)
#             else:
#                 logger.info('vocabulary already exists: %s - %s', 
#                     la.affiliation_category, defined_vocabs[la.affiliation_category])
#  
#         _scope = 'labaffiliation.category.%s'
#         count_updated = 0
#         for i,la in enumerate(apps.get_model('db', 'LabAffiliation')
#                 .objects.all()
#                 .order_by('affiliation_category','affiliation_name')):
#             if la.affiliation_name not in defined_vocabs:
#                 name = la.affiliation_name.lower()
#                 for replacer,replacement in replace_phrases:
#                     logger.info(
#                         'replacer: %s, replacement: %s, name: %s', 
#                         str(replacer),replacement,name)
#                     name = replacer.sub(replacement.lower(),name)
#                     logger.info('new name: %s', name)    
#                  
#                 title = la.affiliation_name
#                 key = default_converter(name)
#                 scope = _scope % default_converter(la.affiliation_category)
#                 ordinal = len(defined_vocabs) + 1
#                 row = [key, scope, ordinal, title, la.affiliation_category ]
#                 defined_vocabs[la.affiliation_name] = dict(zip(header,row))
#                 vocab_writer.writerow(row)
#                 logger.debug('created row: %s', row)
#             else:
#                 logger.info('vocabulary already exists: %s - %s', 
#                     la.affiliation_name, defined_vocabs[la.affiliation_name])
#                  
#             # now set the screensaveruser field
#             ScreensaverUser = apps.get_model('db','ScreensaverUser')
#             if la.labhead_set.exists():
#                 for lh in la.labhead_set.all():
#                     su = ScreensaverUser.objects.get(
#                         screensaver_user_id=lh.screensaver_user_id)
#                     new_value = defined_vocabs[la.affiliation_name]['key']
#                     logger.debug(
#                         'updating user %s, lab_affiliation: %s', 
#                         su.username,new_value )
#                     su.lab_head_affiliation = new_value;
#                     su.save()
#                     count_updated += 1
#                  
#                  
#     logger.info('labaffiliation vocabulary creation done, updated: %s users, %s vocabs',
#          count_updated, len(defined_vocabs))

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0006_screen'),
    ]

    operations = [

        # NOTE: see 0003 for schema field migrations for lab head
        # - because of transactions, schema migration must be elsewhere        
        migrations.RunSQL(
            'UPDATE screensaver_user su '
            ' set lab_head_id=sru.lab_head_id '
            ' from  screening_room_user sru '
            ' where sru.screensaver_user_id=su.screensaver_user_id'),
        migrations.RunSQL(
            'UPDATE screensaver_user su '
            ' set lab_head_id=lh.screensaver_user_id '
            ' from  lab_head lh '
            ' where lh.screensaver_user_id=su.screensaver_user_id'),
         
        migrations.RunSQL(
            'UPDATE screensaver_user su '
            ' set lab_affiliation_id=lh.lab_affiliation_id '
            ' from  lab_head lh '
            ' where lh.screensaver_user_id=su.screensaver_user_id'),
         
        # TODO: drop the LabHead and ScreeningRoomUser models - 20170705
        migrations.RunSQL(
            'ALTER TABLE lab_head '
            'DROP CONSTRAINT fk_lab_head_to_screening_room_user'),
        migrations.RunSQL(
            'ALTER TABLE screening_room_user '
            'DROP CONSTRAINT fk_screening_room_user_to_lab_head'),
        
        
    ]
