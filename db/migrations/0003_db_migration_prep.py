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

#####
# Migration prep:
# Perform migrations required so that Django may use the legacy schema for 
# subsequent actions
# * This migration must be performed before the migration bootstrap step
# * see migration.sh
#####

vocab_replace_map = {
    'ICCB-L/NSRB staff': 'ICCB-L staff',
    }

vocab_key_replace_map = {
    'iccb_l_staff': 'staff',
    'hms_affiliate_with_hms_quad_appointment': 'hms_affiliate_hms_quad',
    }
vocab_ignore_map = {
    'user.classification': ['unassigned',]}

def create_vocab(vocab_writer, attr, scope, query, write_to_file=True):
    logger.info('create simple vocab: %s, %s', attr,scope)
    vocabs = []
    for ordinal, attr_value in (
            enumerate(query.values_list(attr, flat=True)
                .distinct(attr).order_by(attr))):
        if not attr_value: continue
        
        title = attr_value
        if attr_value in vocab_replace_map:
            title = vocab_replace_map[attr_value]
        key = default_converter(title)
        
        if scope in vocab_ignore_map:
            if key in vocab_ignore_map[scope]:
                continue
            
        if key in vocab_key_replace_map:
            key = vocab_key_replace_map[key]

        query.filter(**{ '%s__exact' % attr: attr_value }).update(**{ attr: key })
        
        vocabs.append([key, scope, ordinal, title])

    for row in vocabs:
        if write_to_file is True:
            vocab_writer.writerow(row)
        logger.info('updated vocab: %r' % row)

def create_serviceactivity_vocab(vocab_writer, attr, scope, query):
    logger.info('create simple vocab: %s, %s', attr,scope)
    vocabs = []
    for ordinal, attr_value in (enumerate(
        query.values_list(attr, flat=True)
            .distinct(attr).order_by(attr))):
        if not attr_value: continue
        key = default_converter(attr_value)
        title = attr_value
        vocabs.append([key, scope, ordinal, title])
    for row in vocabs:
        title = row[3]
        key = row[0]
        # NOTE: do not run update; this is the second run of the service activity
        #         query.filter(**{ '%s__exact' % attr: title }).update(**{ attr: key })
        vocab_writer.writerow(row)
        logger.info('updated vocab: %r' % row)


def create_simple_vocabularies(apps):

    # simple vocab cases: update without linked tables
    
    vocab_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'vocabulary_data_generated.csv')

    logger.info('write vocabularies to %s' % vocab_file)
    
    with open(vocab_file, 'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['key', 'scope', 'ordinal', 'title'] 
        vocab_writer.writerow(header)
        # Run it twice for service activities, so that they can be separated
        # from the vanilla activities; 
        create_serviceactivity_vocab(
            vocab_writer, 'service_activity_type', 'serviceactivity.type',
            apps.get_model('db', 'ServiceActivity').objects.all())
        
        input_args = [
                ['species', 'screen.species',
                    apps.get_model('db', 'Screen').objects.all()],
                ['assay_type', 'screen.assay_type',
                    apps.get_model('db', 'Screen').objects.all()],
                ['study_type', 'study.type',
                    apps.get_model('db', 'Screen').objects.all()],
                ['assay_readout_type', 'datacolumn.assay_readout_type',
                    apps.get_model('db', 'DataColumn').objects.all()],
                ['value', 'funding_support',
                    apps.get_model('db', 'FundingSupport').objects.all()],
                ['service_activity_type', 'activity.type',
                    apps.get_model('db', 'ServiceActivity').objects.all()],
                # dep: 0020(provisional) -> moves to 0002: move classification column
                ['classification', 'user.classification',
                    apps.get_model('db', 'ScreensaverUser').objects.all()],
                ['lab_head_appointment_category', 'lab_head.appointment_category',
                    apps.get_model('db', 'ScreensaverUser').objects.all()],
                ['lab_head_appointment_department', 'lab_head.appointment_department',
                    apps.get_model('db', 'ScreensaverUser').objects.all()],
                ['value', 'cell_line',
                    apps.get_model('db', 'CellLine').objects.all()],
                ['value', 'transfection_agent',
                    apps.get_model('db', 'TransfectionAgent').objects.all()],
                ['assay_protocol_type', 'assayprotocol.type',
                    apps.get_model('db', 'Screening').objects.all()],

                ['screen_type', 'screen.type',
                    apps.get_model('db', 'Screen').objects.all()],
#                 ['project_phase', 'screen.project_phase',
#                     apps.get_model('db', 'Screen').objects.all()],
                
                ['screening_status', 'library.screening_status',
                    apps.get_model('db', 'Library').objects.all()],
                ['library_type', 'library.type',
                    apps.get_model('db', 'Library').objects.all()],
                ['solvent', 'library.solvent',
                    apps.get_model('db', 'Library').objects.all()],
                # NOTE: library.screen_type is the same as screen.type
                ['screen_type', 'library.screen_type',
                    apps.get_model('db', 'Library').objects.all(),False],

                ['affiliation_category', 'labaffiliation.category',
                    apps.get_model('db', 'LabAffiliation').objects.all()],

                # NOTE: plate.type vocabulary is being handled in the manual migration
                # ['plate_type', 'plate.type', 
                #     apps.get_model('db', 'Plate'),objects.all(), False]
                # ['assay_plate_type', 'plate.type', 
                #     apps.get_model('db', 'CherryPickRequest'),objects.all(), False]
                # ['assay_plate_type', 'plate.type', 
                #     apps.get_model('db', 'CherryPickAssayPlate'),objects.all(), False]

            ]            
        for arg_list in input_args:
            create_vocab(vocab_writer, *arg_list)

    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file, api_init_actions_file))
    with open(api_init_actions_file, 'a+') as _file:
        new_row = ['patch', 'vocabulary', os.path.basename(vocab_file)]
        reader = csv.reader(_file)
        found = False
        for row in reader:
            if row == new_row:
                found = True
                break
        if not found:
            writer = csv.writer(_file)
            writer.writerow(new_row)
        else:
            logger.info('api_init entry for row already created: %r' % row)

    # update the new service_activity.funding_support
    # see below for manual sql for screen funding supports
    for fs in apps.get_model('db', 'FundingSupport').objects.all():
        apps.get_model('db', 'ServiceActivity').objects.all().filter(
            funding_support_link=fs).update(funding_support=fs.value)
    
    
def create_attached_file_type_vocab(apps):
    
    vocab_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'vocabulary_attachedfiletype_data.csv')
    logger.info('write vocabularies to %s' % vocab_file)

    replace_map = {
        'ICCB-L/NSRB RNAi User Agreement': 'ICCB-L RNAi User Agreement',
        '2009 ICCB-L/NSRB Small Molecule User Agreement': 
            '2009 ICCB-L Small Molecule User Agreement',
        '2010 ICCB-L/NSRB Small Molecule User Agreement': 
            'ICCB-L Small Molecule User Agreement',
        'ICCB-L/NSRB Application (user)': 'ICCB-L Application (user)'
    }
    
    with open(vocab_file, 'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['key', 'scope', 'ordinal', 'title', 'is_retired'] 
        vocab_writer.writerow(header)

        _scope = 'attachedfiletype.%s'
        _ordinal = 1
        for_entity_type = None
        for i, obj in enumerate(apps.get_model('db', 'AttachedFileType')
                .objects.all().order_by('for_entity_type','value')):
            
            if for_entity_type != obj.for_entity_type:
                for_entity_type = obj.for_entity_type
                _ordinal = 1
            
            value = obj.value
            
            if value in replace_map:
                value = replace_map[value]
            key = default_converter(value.lower())
            
            scope = _scope % obj.for_entity_type
            title = value
            ordinal = _ordinal
            _ordinal += 1
            is_retired = key in [
                'iccb_l_nsrb_application',
                '2009_iccb_l_small_molecule_user_agreement',
                'marcus_application','miare_document',
                'nerce_screener_supplies_list']
            row = [key, scope, ordinal, title, is_retired] 
            vocab_writer.writerow(row)
            logger.info('created: %r', row)
            
            (apps.get_model('db', 'AttachedFile').objects
                .filter(attached_file_type=obj).update(type=key))

    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file, api_init_actions_file))
    with open(api_init_actions_file, 'a+') as _file:
        new_row = ['patch', 'vocabulary', os.path.basename(vocab_file)]
        reader = csv.reader(_file)
        found = False
        for row in reader:
            if row == new_row:
                found = True
                break
        if not found:
            writer = csv.writer(_file)
            writer.writerow(new_row)
        else:
            logger.info('api_init entry for row already created: %r' % row)

    logger.info('done')


def create_checklist_vocabularies(apps):
    
    vocab_file = os.path.join(
        PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'vocabulary_checklists_data.csv')
    logger.info('write vocabularies to %s' % vocab_file)
    with open(vocab_file, 'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['key', 'scope', 'ordinal', 'title','comment'] 
        vocab_writer.writerow(header)

        scope = 'userchecklist.name'
        for i,obj in enumerate(apps.get_model('db', 'ChecklistItem')
                .objects.all()
                .order_by('checklist_item_group','order_statistic')):
            item_name = obj.item_name
            key = default_converter(item_name)
            comment = 'group: %s' % obj.checklist_item_group
            title = item_name
            ordinal = i
            row = [key, scope, ordinal, title, comment ]
            vocab_writer.writerow(row)
            logger.info('created: %r', row)
    
    
def create_vocabularies(apps, schema_editor):
    
    create_simple_vocabularies(apps)
    create_attached_file_type_vocab(apps)
    create_checklist_vocabularies(apps)
    

def update_facility_usage_roles(apps, schema_editor):
    UserFacilityUsageRole = apps.get_model('db', 'UserFacilityUsageRole')
    
    role_updates = (
#         ('smallMoleculeScreener', 'small_molecule_screener'),
#         ('rnaiScreener', 'rnai_screener'),
        ('iccblProjectUser', 'iccbl_project_user'),
        ('analyticalChemistryUser', 'analytical_chemistry_user'),
        ('nonScreeningUser', 'non_screening_user'),
        ('mouseImagerUser', 'mouse_image_user'),
        ('qpcrUser', 'qpcr_user'),
        ('medicinalChemistUser', 'medicinal_chemist_user'),
    )
    
    for ru in role_updates:
        (UserFacilityUsageRole.objects.all()
            .filter(facility_usage_role=ru[0])
            .update(facility_usage_role=ru[1]))

    UserFacilityUsageRole.objects.all()\
        .filter(facility_usage_role='smallMoleculeScreener').delete()
    UserFacilityUsageRole.objects.all()\
        .filter(facility_usage_role='rnaiScreener').delete()

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0002_db_migration_prep'),
    ]

    operations = [
        migrations.RunPython(create_vocabularies),
        migrations.RunSQL('delete from screen_funding_supports;'),
        migrations.RunSQL('''
            insert into screen_funding_supports (id,screen_id,funding_support) (
              select nextval('screen_funding_supports_id_seq'), screen.screen_id, value
              from screen
              join screen_funding_support_link fsl on(screen.screen_id=fsl.screen_id)
              join funding_support fs on(fsl.funding_support_id=fs.funding_support_id) 
              order by screen.screen_id,value
              );
        '''.strip()),
        migrations.RunSQL('delete from user_facility_usage_role;'),
        migrations.RunSQL('''
            insert into user_facility_usage_role (id,screensaver_user_id,facility_usage_role) (
              select nextval('user_facility_usage_role_id_seq'), su.screensaver_user_id, facility_usage_role
              from screening_room_user su
              join screening_room_user_facility_usage_role roles on(su.screensaver_user_id=roles.screening_room_user_id)
              order by su.screensaver_user_id,roles.facility_usage_role
              );
        '''.strip()),
        # Transfer the cell_line vocab from the "screen_cell_line" table to the 
        # "screen_cell_lines" vocab table
        migrations.RunSQL('delete from screen_cell_lines;'),
        migrations.RunSQL('''
            insert into screen_cell_lines (id,screen_id,cell_line) (
              select nextval('screen_cell_lines_id_seq'), screen.screen_id, value
              from screen
              join screen_cell_line scl on(screen.screen_id=scl.screen_id)
              join cell_line cl on(scl.cell_line_id=cl.cell_line_id) 
              order by screen.screen_id,value
              );
        '''.strip()),
        # Transfer the transfection_agent vocab from the transfection_agent
        # table to the screen.transfection_agent field
        migrations.RunSQL('''
            update screen set transfection_agent_text = value 
            from transfection_agent ta 
            where ta.transfection_agent_id=screen.transfection_agent_id;
        '''.strip()),
           
#         migrations.RemoveField(
#             model_name='screen',
#             name='transfection_agent',
#         ),
#         migrations.RenameField(
#             model_name='screen', 
#             old_name='transfection_agent_text', 
#             new_name='transfection_agent'
#         ),
        
        # Lab affiliation migration prep: data migration in 0007
        migrations.RenameField(
            model_name='labaffiliation',
            old_name='affiliation_category',
            new_name='category',
        ),
        migrations.RenameField(
            model_name='labaffiliation',
            old_name='affiliation_name',
            new_name='name',
        ),
        
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


        migrations.RunPython(update_facility_usage_roles),
    ]
