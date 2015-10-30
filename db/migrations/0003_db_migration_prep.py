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


def create_vocab(vocab_writer,attr,scope,query):
    resource_uri = '/reports/api/v1/vocabularies/%s/%s/'

    vocabs = []
    for ordinal,attr_value in (enumerate( 
        query.values_list(attr, flat=True)
            .distinct(attr).order_by(attr)) ):
        if not attr_value: continue
        key = default_converter(attr_value)
        title = attr_value
        _resource_uri = resource_uri % (scope,key)
        vocabs.append([_resource_uri,key,scope,ordinal,title])
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
        'db','static','api_init','vocabularies_data_generated.csv')

    logger.info('write vocabularies to %s' % vocab_file)
    
    with open(vocab_file,'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['resource_uri','key','scope','ordinal','title'] 
        vocab_writer.writerow(header)
        
        input_args = [
                ['species', 'screen.species',
                    apps.get_model('db','Screen').objects.all()],
                ['assay_type', 'screen.assay_type',
                    apps.get_model('db','Screen').objects.all()],
                ['assay_readout_type', 'datacolumn.assay_readout_type',
                    apps.get_model('db','DataColumn').objects.all()],
                ['value', 'funding_support',
                    apps.get_model('db','FundingSupport').objects.all()],
                ['species', 'screen.species',
                    apps.get_model('db','Screen').objects.all()],
                ['service_activity_type', 'serviceactivity.type',
                    apps.get_model('db','ServiceActivity').objects.all()],
            ]            
        for arg_list in input_args:
            create_vocab(vocab_writer,*arg_list)

    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db','static','api_init','api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file,api_init_actions_file))
    with open(api_init_actions_file,'a+') as _file:
        new_row = ['patch','vocabularies',os.path.basename(vocab_file)]
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
    for fs in apps.get_model('db','FundingSupport').objects.all():
        apps.get_model('db','ServiceActivity').objects.all().filter(
            funding_support_link=fs).update(funding_support=fs.value)
    
    
def create_attached_file_type_vocab(apps):
    
    vocab_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db','static','api_init','vocabularies_attachedfiletype_data.csv')
    logger.info('write vocabularies to %s' % vocab_file)
    resource_uri = '/reports/api/v1/vocabularies/%s/%s/'
    with open(vocab_file,'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['resource_uri','key','scope','ordinal','title'] 
        vocab_writer.writerow(header)

        _scope = 'attachedfiletype.%s'
        for i,obj in enumerate(apps.get_model('db','AttachedFileType')
                .objects.all().order_by('value')):
            key = default_converter(obj.value)
            scope = _scope % obj.for_entity_type
            title = obj.value
            ordinal = i
            row = [resource_uri % (scope,key),key,scope,ordinal,title] 
            vocab_writer.writerow(row)
            logger.info(str(('created', row)))
            
            (apps.get_model('db','AttachedFile').objects
                .filter(attached_file_type=obj).update(type=key))

    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db','static','api_init','api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file,api_init_actions_file))
    with open(api_init_actions_file,'a+') as _file:
        new_row = ['patch','vocabularies',os.path.basename(vocab_file)]
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
        'db','static','api_init','vocabularies_checklists_data.csv')
    logger.info('write vocabularies to %s' % vocab_file)
    resource_uri = '/reports/api/v1/vocabularies/%s/%s/'
    with open(vocab_file,'w') as _file:
        vocab_writer = csv.writer(_file)
        header = ['resource_uri','key','scope','ordinal','title'] 
        vocab_writer.writerow(header)

        # ci_group_map = {}
        scope = 'checklistitem.group'
        default_order = ['mailing', 'forms', 'non-harvard', 'imaging','legacy']
        for obj in (apps.get_model('db','ChecklistItem')
                .objects.all().distinct('checklist_item_group')):
            key = default_converter(obj.checklist_item_group)
            ordinal = 0
            for i,x in enumerate(default_order):
                if x in obj.checklist_item_group.lower():
                    ordinal=i 
                    break
                
            title = obj.checklist_item_group
            row = [resource_uri % (scope,key),key,scope,ordinal,title]
            vocab_writer.writerow(row)
            #ci_group_map[obj.checklist_item_group] = key 
            logger.info(str(('created', row)))

        _scope = 'checklistitem.%s.name'
        # ci_name_map = {}
        for obj in (apps.get_model('db','ChecklistItem')
                .objects.all().distinct('item_name')):
            key = default_converter(obj.item_name)
            scope = _scope % default_converter(obj.checklist_item_group)
            
            title = obj.item_name
            ordinal = obj.order_statistic
            row = [resource_uri % (scope,key),key,scope,ordinal,title]
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
            _dict['scope']=scope
            row = [resource_uri % (_dict['scope'],_dict['key']),
                _dict['key'],_dict['scope'],_dict['ordinal'],_dict['title']]
            vocab_writer.writerow(row)
            logger.info(str(('created', row)))
    
    api_init_actions_file = os.path.join(
        lims.settings.PROJECT_ROOT, '..',
        'db','static','api_init','api_init_actions.csv')
    logger.info('write %s entry to %s' % (vocab_file,api_init_actions_file))
    with open(api_init_actions_file,'a+') as _file:
        new_row = ['patch','vocabularies',os.path.basename(vocab_file)]
        reader = csv.reader(_file)
        found = False
        for row in reader:
            if row == new_row:
                found = True
                break
        if not found:
            logger.info('write api_init row: %s' % new_row )
            writer = csv.writer(_file)
            writer.writerow(new_row)
        else:
            logger.info('api_init entry for row already created: %r' % row)


    logger.info('vocabulary creation done')


def create_vocabularies(apps, schema_editor):
    
    create_simple_vocabularies(apps)
    create_attached_file_type_vocab(apps)
    create_checklist_vocabularies(apps)


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0002_db_migration_prep'),
    ]

    operations = [
        migrations.RunPython(create_vocabularies),
        migrations.AlterField(
            model_name='attachedfile',
            name='type',
            field=models.TextField()),
        migrations.RemoveField(
            model_name='attachedfile',name='attached_file_type_id'),
    ]
