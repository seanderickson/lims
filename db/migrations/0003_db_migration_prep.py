# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import logging
import os
import re

from django.db import migrations, models

import lims
from lims.base_settings import PROJECT_ROOT
from db.support.data_converter import default_converter


logger = logging.getLogger(__name__)


def create_vocab(vocab_writer, attr, scope, query):
    resource_uri = '/reports/api/v1/vocabularies/%s/%s/'
    logger.info('create simple vocab: %s, %s', attr,scope)
    vocabs = []
    for ordinal, attr_value in (enumerate(
        query.values_list(attr, flat=True)
            .distinct(attr).order_by(attr))):
        if not attr_value: continue
        key = default_converter(attr_value)
        title = attr_value
        _resource_uri = resource_uri % (scope, key)
        vocabs.append([_resource_uri, key, scope, ordinal, title])
    for row in vocabs:
        title = row[4]
        key = row[1]
        query.filter(**{ '%s__exact' % attr: title }).update(**{ attr: key })
        vocab_writer.writerow(row)
        logger.info('updated vocab: %r' % row)


def create_simple_vocabularies(apps):

    # simple vocab cases: update without linked tables
    
    vocab_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'vocabularies_data_generated.csv')

    logger.info('write vocabularies to %s' % vocab_file)
    
    with open(vocab_file, 'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['resource_uri', 'key', 'scope', 'ordinal', 'title'] 
        vocab_writer.writerow(header)
        
        input_args = [
                ['species', 'screen.species',
                    apps.get_model('db', 'Screen').objects.all()],
                ['assay_type', 'screen.assay_type',
                    apps.get_model('db', 'Screen').objects.all()],
                ['assay_readout_type', 'datacolumn.assay_readout_type',
                    apps.get_model('db', 'DataColumn').objects.all()],
                ['value', 'funding_support',
                    apps.get_model('db', 'FundingSupport').objects.all()],
                ['service_activity_type', 'serviceactivity.type',
                    apps.get_model('db', 'ServiceActivity').objects.all()],
                # dep: 0020(provisional) -> moves to 0002: move classification column
                ['classification', 'user.classification',
                    apps.get_model('db', 'ScreensaverUser').objects.all()],
                ['value', 'cell_line',
                    apps.get_model('db', 'CellLine').objects.all()],
                ['value', 'transfection_agent',
                    apps.get_model('db', 'TransfectionAgent').objects.all()],
            ]            
        for arg_list in input_args:
            create_vocab(vocab_writer, *arg_list)

    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file, api_init_actions_file))
    with open(api_init_actions_file, 'a+') as _file:
        new_row = ['patch', 'vocabularies', os.path.basename(vocab_file)]
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
    
def create_lab_affiliation_vocab(apps):
 
    # populate the title field, change the name field to a key
    replace_phrases = [
        ['harvard medical school', 'hms'],
        ['harvard university', 'harvard'],
        ['European Molecular Biology Laboratory', 'embl'],
        ['[embl]',''],
        ['Dana Farber Cancer Institute', 'dfci'],
        ['University of California', 'uc'],
        ['University of Massachusetts', 'umass'],
        ['Institute of Chemistry and Cell Biology', 'iccb'],
        ['Beth Israel Deaconess Medical Center', 'bidmc'],
        ['Tufts University', 'tufts'],
        ['University of California', 'UC'],
        ['University of Massachusetts', 'umass'],
        ['[NYU]', ''],
        ['the',''],
        ['of',''],
        ['in',''],
        ["women's", 'womens'],
        ["children's", 'childrens']
    ]
    replace_phrases = [[re.compile(r'\b%s\b' % x, re.IGNORECASE),y] 
        for [x,y] in replace_phrases ]
    vocab_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'vocabularies_lab_affiliations_data.csv')
    logger.info('write vocabularies to %s' % vocab_file)
    resource_uri = 'vocabularies/%s/%s/'
    with open(vocab_file, 'a+') as _file:
         
        header = ['resource_uri', 'key', 'scope', 'ordinal', 'title', 'comment'] 
        reader = csv.reader(_file)
        vocab_writer = csv.writer(_file)
        defined_vocabs = {}
        try:
            header = reader.next()
            logger.info('read header: %s', header)
            for row in reader:
                defined_vocabs[row[4]] = dict(zip(header,row))
        except StopIteration as e:
            logger.info('no entries in %s, writing a new file',vocab_file)
            vocab_writer.writerow(header)
 
        scope = 'labaffiliation.category'
        default_order = ['hms','hms_affiliated_hospital','hsph',
            'harvard_fas','broad_icg','other']
        for i,la in enumerate(apps.get_model('db', 'LabAffiliation')
                .objects.all().distinct('affiliation_category')):
            if la.affiliation_category not in defined_vocabs:
                
                key = default_converter(la.affiliation_category)
                if key in default_order:
                    ordinal = default_order.index(key)
                else:
                    ordinal = i + len(default_order)
                title = la.affiliation_category
                row = [resource_uri % (scope, key), key, scope, ordinal, title, 
                    la.affiliation_category ]
                vocab_writer.writerow(row)
                defined_vocabs[la.affiliation_category] = dict(zip(header,row))
                logger.debug('created %s', row)
            else:
                logger.info('vocabulary already exists: %s - %s', 
                    la.affiliation_category, defined_vocabs[la.affiliation_category])
 
        _scope = 'labaffiliation.category.%s'
        count_updated = 0
        for i,la in enumerate(apps.get_model('db', 'LabAffiliation')
                .objects.all()
                .order_by('affiliation_category','affiliation_name')):
            if la.affiliation_name not in defined_vocabs:
                name = la.affiliation_name.lower()
                for replacer,replacement in replace_phrases:
                    logger.info('replacer: %s, replacement: %s, name: %s', str(replacer),replacement,name)
                    name = replacer.sub(replacement.lower(),name)
                    logger.info('new name: %s', name)    
                 
                title = la.affiliation_name
                key = default_converter(name)
                scope = _scope % default_converter(la.affiliation_category)
                ordinal = len(defined_vocabs) + 1
                row = [resource_uri % (scope, key), key, scope, ordinal, title, 
                    la.affiliation_category ]
                defined_vocabs[la.affiliation_name] = dict(zip(header,row))
                vocab_writer.writerow(row)
                logger.debug('created row: %s', row)
            else:
                logger.info('vocabulary already exists: %s - %s', 
                    la.affiliation_name, defined_vocabs[la.affiliation_name])
                 
            # now set the screensaveruser field
            ScreensaverUser = apps.get_model('db','ScreensaverUser')
            if la.labhead_set.exists():
                for lh in la.labhead_set.all():
                    su = ScreensaverUser.objects.get(screensaver_user_id=lh.screensaver_user_id)
                    new_value = defined_vocabs[la.affiliation_name]['key']
                    logger.debug('updating user %s, lab_affiliation: %s', su.username,new_value )
                    su.lab_head_affiliation = new_value;
                    su.save()
                    count_updated += 1
                 
                 
    logger.info('labaffiliation vocabulary creation done, updated: %s users, %s vocabs',
         count_updated, len(defined_vocabs))
    
def create_attached_file_type_vocab(apps):
    
    vocab_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'vocabularies_attachedfiletype_data.csv')
    logger.info('write vocabularies to %s' % vocab_file)
    resource_uri = '/reports/api/v1/vocabularies/%s/%s/'
    with open(vocab_file, 'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['resource_uri', 'key', 'scope', 'ordinal', 'title'] 
        vocab_writer.writerow(header)

        _scope = 'attachedfiletype.%s'
        for i, obj in enumerate(apps.get_model('db', 'AttachedFileType')
                .objects.all().order_by('value')):
            key = default_converter(obj.value)
            scope = _scope % obj.for_entity_type
            title = obj.value
            ordinal = i
            row = [resource_uri % (scope, key), key, scope, ordinal, title] 
            vocab_writer.writerow(row)
            logger.info(str(('created', row)))
            
            (apps.get_model('db', 'AttachedFile').objects
                .filter(attached_file_type=obj).update(type=key))

    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file, api_init_actions_file))
    with open(api_init_actions_file, 'a+') as _file:
        new_row = ['patch', 'vocabularies', os.path.basename(vocab_file)]
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
    # create a separate vocab file: checklist_item_vocab, add to api_init.csv
    # output vocabs into a vocabulary patch file
    vocab_file = os.path.join(
        PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'vocabularies_checklists_data.csv')
    logger.info('write vocabularies to %s' % vocab_file)
    resource_uri = '/reports/api/v1/vocabularies/%s/%s/'
    with open(vocab_file, 'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['resource_uri', 'key', 'scope', 'ordinal', 'title', 
            'expire_interval_days', 'expire_notifiy_days', ''] 
        vocab_writer.writerow(header)

        # ci_group_map = {}
        scope = 'checklistitem.group'
        default_order = ['mailing', 'forms', 'non-harvard', 'imaging', 'legacy']
        for obj in (apps.get_model('db', 'ChecklistItem')
                .objects.all().distinct('checklist_item_group')):
            key = default_converter(obj.checklist_item_group)
            ordinal = 0
            for i, x in enumerate(default_order):
                if x in obj.checklist_item_group.lower():
                    ordinal = i 
                    break
                
            title = obj.checklist_item_group
            row = [resource_uri % (scope, key), key, scope, ordinal, title, None, None]
            vocab_writer.writerow(row)
            # ci_group_map[obj.checklist_item_group] = key 
            logger.info(str(('created', row)))

        _scope = 'checklistitem.%s.name'
        # ci_name_map = {}
        
        
        for obj in (apps.get_model('db', 'ChecklistItem')
                .objects.all().distinct('item_name')):
            key = default_converter(obj.item_name)
            scope = _scope % default_converter(obj.checklist_item_group)
            
            # NOTE: fore user_checklist_item overload of vocabularies:
            if key in ('current_rnai_user_agreement_active',
                'current_small_molecule_user_agreement_active'):
                expire_interval_days = 720 
                expire_notifiy_days = 60
            else:
                expire_interval_days = None 
                expire_notifiy_days = None
                
            title = obj.item_name
            ordinal = obj.order_statistic
            row = [resource_uri % (scope, key), key, scope, ordinal, title, 
                expire_interval_days, expire_notifiy_days ]
            vocab_writer.writerow(row)
            # ci_name_map[obj.item_name] = key
            logger.info(str(('created', row)))
            
        scope = 'checklistitem.status'
        status_values = [
            {
                'key': 'not_completed',
                'title': 'Not Completed',
                'ordinal': 0
            },
            {
                'key': 'activated',
                'title': 'Activated',
                'ordinal': 1
            },
            {
                'key': 'deactivated',
                'title': 'Deactivated',
                'ordinal': 2
            },
            {
                'key': 'na',
                'title': 'N/A',
                'ordinal': 3
            },
            {
                'key': 'completed',
                'title': 'Completed',
                'ordinal': 4
            },
        ]
        
        for _dict in status_values:
            _dict['scope'] = scope
            row = [resource_uri % (_dict['scope'], _dict['key']),
                _dict['key'], _dict['scope'], _dict['ordinal'], _dict['title'], None, None]
            vocab_writer.writerow(row)
            logger.info(str(('created', row)))
    
    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db', 'static', 'api_init', 'api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file, api_init_actions_file))
    with open(api_init_actions_file, 'a+') as _file:
        new_row = ['patch', 'vocabularies', os.path.basename(vocab_file)]
        reader = csv.reader(_file)
        found = False
        for row in reader:
            if row == new_row:
                found = True
                break
        if not found:
            logger.info('write api_init row: %s' % new_row)
            writer = csv.writer(_file)
            writer.writerow(new_row)
        else:
            logger.info('api_init entry for row already created: %r' % row)


    logger.info('vocabulary creation done')


def create_vocabularies(apps, schema_editor):
    
    create_simple_vocabularies(apps)
    create_attached_file_type_vocab(apps)
    create_checklist_vocabularies(apps)

    create_lab_affiliation_vocab(apps)
    
    

def update_facility_usage_roles(apps, schema_editor):
    UserFacilityUsageRole = apps.get_model('db', 'UserFacilityUsageRole')
    
    role_updates = (
        ('smallMoleculeScreener', 'small_molecule_screener'),
        ('rnaiScreener', 'rnai_screener'),
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
            update screen set transfection_agent = value 
            from transfection_agent ta 
            where ta.transfection_agent_id=screen.transfection_agent_id;
        '''.strip()),
        migrations.RunPython(update_facility_usage_roles),
    ]
