# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from collections import OrderedDict, defaultdict

from django.db import migrations, models

import json
import logging
import re

logger = logging.getLogger(__name__)

screen_log_cols = OrderedDict({
    'activity_id': 'a.activity_id',
    'username': 'username',
    'email': 'email',
    'screensaver_user_id': 'screensaver_user_id',
    'date_of_activity': 'date_of_activity',
    'comments': 'a.comments',
    'screen_facility_id': 'screen.facility_id',
    'screen_id': 'screen.screen_id',
    })
cpr_log_cols = screen_log_cols.copy()
cpr_log_cols.update({'cherry_pick_request_id': 'cpr.cherry_pick_request_id'})

screen_log_sql =\
''' 
    from activity a 
    join screensaver_user on(performed_by_id=screensaver_user_id) 
    join administrative_activity using(activity_id) 
    join screen_update_activity on(update_activity_id=activity_id)
    join screen using(screen_id)  
    where administrative_activity_type='Plate Raw Data Transformation' 
    order by screen.screen_id '''
cpr_log_sql =\
''' 
    from activity a 
    join screensaver_user on(performed_by_id=screensaver_user_id) 
    join administrative_activity using(activity_id) 
    join cherry_pick_request_update_activity on(update_activity_id=activity_id)
    join cherry_pick_request cpr using(cherry_pick_request_id)
    join screen using(screen_id)  
    where administrative_activity_type='Plate Raw Data Transformation' 
    order by cpr.cherry_pick_request_id '''



def migrate_screen_transformer_logs(apps, schema_editor):
    '''
    Migrate JSON data from the 
    administrative_activity_type='Plate Raw Data Transformation' logs for 
    Screens;
    NOTE: join tables (screen_update_activity) are not valid Django models, so 
    migrate using raw sql.
    '''
    cols = screen_log_cols
        
    colkeys = cols.keys()
    _cols = ', '.join([ '%s as %s' % (value,key) for key, value in cols.items() ])
    sql = 'select ' + _cols + screen_log_sql
    
    migrate_transformer_logs(apps, schema_editor, sql, colkeys)
    
def migrate_cpr_transformer_logs(apps, schema_editor):
    '''
    Migrate JSON data from the 
    administrative_activity_type='Plate Raw Data Transformation' logs for 
    CherryPickRequests;
    NOTE: join tables (cherry_pick_request_update_activity) are not valid Django models, so 
    migrate using raw sql.
    '''
    cols = cpr_log_cols
        
    colkeys = cols.keys()
    _cols = ', '.join([ '%s as %s' % (value,key) for key, value in cols.items() ])
    sql = 'select ' + _cols + cpr_log_sql
    
    migrate_transformer_logs(apps, schema_editor, sql, colkeys)
    
    
def migrate_transformer_logs(apps, schema_editor, sql, colkeys):
        
    RawDataTransform = apps.get_model('db','RawDataTransform')
    RawDataInputFile = apps.get_model('db','RawDataInputFile')
    
    connection = schema_editor.connection
    cursor = connection.cursor()
        
    # create some stats:
    facility_ids = set()
    cpr_ids = set()
    transform_fields = set()
    input_file_fields = set()
    output_types = set()
    assay_plate_sizes = set()
    readout_types = set()
    readouts = set()
    conditions = set()
    collation_orderings = set()
    replicates = set()
    
    transform_field_map = {
        'outputFormatSelection': 'output_sheet_option',
        'outputFormat': 'output_sheet_option',
        'outputFileName': 'output_filename',
        'assayPlateSize': 'assay_plate_size',
        'negativeControls': 'assay_negative_controls',
        'assayNegativeControls': 'assay_negative_controls',
        'positiveControls': 'assay_positive_controls',
        'assayPositiveControls': 'assay_positive_controls',
        'assayOtherControls': 'assay_other_controls',
        'otherControls': 'assay_other_controls',
        'libraryControls': 'library_controls', 
        'plates': 'plate_ranges',        
    }
    transform_vocabs = {
        'output_sheet_option': {
            'PlatePerWorksheet': 'plate_per_worksheet',
            'AllPlatesInSingleWorksheet':'all_plates_in_single_worksheet'
        },
        'assay_plate_size': {
            'WELLS_1536': 1536, 
            'WELLS_384': 384,
            'WELLS_96': 96 }
    }
    
    input_file_field_map = {
        'readouts': 'readouts',
        'replicates': 'replicates',
        'readoutTypeSelection': 'readout_type',
        'conditions': 'conditions',
        'collationOrderOrdering': 'collation_order', 
        'uploadedFilename': 'filename',
        }
    input_file_vocabs = {
        'collation_order': {
            'P,R,C,R': 'PQ_REP_C_READ',
            'P,Q,R,C,R': 'PQ_REP_C_READ',
            'P,C,R,R': 'PQ_C_REP_READ', 
            'P,Q,C,R,R': 'PQ_C_REP_READ',
            'C,P,Q,R,R': 'C_PQ_REP_READ', 
            'C,P,R,R': 'C_PQ_REP_READ',  
            'C,R,P,R': 'C_REP_PQ_READ', 
            'C,R,P,Q,R': 'C_REP_PQ_READ', 
            'R,P,Q,C,R': 'REP_PQ_C_READ', 
            'R,C,P,Q,R': 'REP_C_PQ_READ', 
            }
        
        }
    cursor.execute(sql)
    _list = cursor.fetchall()
    if len(_list) == 0:
        raise Exception('No screen_update_activities found with sql: %r' % sql)
    input_file_count = 0
    
    # Since only the last RawDataTransform activity matters, overwrite previous
    previous_rdt = None
    
    for i,_data in enumerate(_list):
        _activity = dict(zip(colkeys, _data))
        
        screen_facility_id = _activity.get('screen_facility_id', None)
        screen_id = _activity.get('screen_id', None)
        cherry_pick_request_id = _activity.get('cherry_pick_request_id', None)
        if screen_facility_id is None and cherry_pick_request_id is None:
            raise Exception('no screen or cpr: %r', _activity)
        
        facility_ids.add(screen_facility_id)
        cpr_ids.add(cherry_pick_request_id)

        
        raw_parts = _activity['comments']
        if raw_parts is not None:
            # semicolon separates RawDataTransform from RawDataInputFile
            raw_parts = raw_parts.split(';')
        else:
            logger.warn('no comments: %r', _activity)
            continue
        
        transform_part = json.loads(raw_parts[0])
        
        output_types.add(transform_part.get('outputFormatSelection',None))
        assay_plate_sizes.add(transform_part.get('assayPlateSize',None))
        logger.debug(
            'parsed for screen: %r: cpr: %r, %r', 
            screen_facility_id,cherry_pick_request_id, transform_part)
        
        # Since only the last RawDataTransform activity matters, overwrite previous
        rdt = None
        if previous_rdt is not None:
            if cherry_pick_request_id is not None:
                if cherry_pick_request_id == previous_rdt.cherry_pick_request_id:
                    rdt = previous_rdt
                    rdt.rawdatainputfile_set.all().delete()
            elif screen_id is not None:
                if screen_id == previous_rdt.screen_id:
                    rdt = previous_rdt
                    rdt.rawdatainputfile_set.all().delete()
                        
        if rdt is None:
            if cherry_pick_request_id is not None:
                rdt = RawDataTransform(
                    cherry_pick_request_id=cherry_pick_request_id)
            elif screen_id is not None:
                rdt = RawDataTransform(
                    screen_id=screen_id)
                
        for k,v in transform_part.items():
            if k in transform_field_map:
                new_key = transform_field_map[k]
                if v is not None and new_key in transform_vocabs:
                    new_val = transform_vocabs[new_key][v]
                    setattr(rdt,new_key,new_val)
                else:
                    setattr(rdt,new_key,v)
        rdt.save()
        previous_rdt = rdt
        
        input_file_parts = []
        if len(raw_parts) > 1:
            for j,raw_part in enumerate(raw_parts[1:]):
                input_file_count += 1
                input_file_part = json.loads(raw_part)
                input_file_fields.update(input_file_part.keys())

                rdif = RawDataInputFile(
                    raw_data_transform=rdt,
                    ordinal=j)
                
                readout_type = input_file_part.get('readoutTypeSelection',None)
                readout_types.add(readout_type)
                rdif.readout_type = readout_type.lower()
                
                collation_order = ','.join(
                    [ x[0] for x in input_file_part.get('collationOrderOrdering',None)])
                collation_orderings.add(collation_order)
                rdif.collation_order = input_file_vocabs['collation_order'][collation_order]
                
                readout = input_file_part.get('readouts', None)
                if readout:
                    readout = ','.join(re.split('[\s,]+', readout))
                readouts.add(readout)
                rdif.readouts = readout
                
                replicate = input_file_part.get('replicates', None)
                replicates.add(replicate)
                rdif.replicates = replicate
                
                condition = input_file_part.get('conditions',None)
                if condition:
                    condition = ','.join(re.split('[\s,]+',condition))
                conditions.add(condition)
                rdif.conditions = condition
                
                rdif.filename = input_file_part.get('uploadedFilename',None)
                
                rdif.save()
                
                input_file_parts.append(input_file_part)

        
        logger.debug('parsed input file parts: %r', input_file_parts)
        
        transform_fields.update(transform_part.keys())
        
        # Example:
        #  {u'outputFormatSelection': u'PlatePerWorksheet', 
        # u'negativeControls': u'1, 2', u'positiveControls': u'24', 
        # u'outputFileName': u'1110_2089-2090', u'otherControls': u'23', 
        # u'plates': u'2089-2090'}
        # Example:
        # {u'readouts': u'', u'readoutTypeSelection': u'FLUORESCENCE_INTENSITY', 
        # u'replicates': 2, u'conditions': u'', 
        # u'collationOrderOrdering': [u'Plates', u'Conditions', u'Replicates', u'Readouts'], 
        # u'uploadedFilename': u'1110-raw.txt'}
        
        # transform_fields: 
        # set([u'outputFormatSelection', u'assayPlateSize', u'negativeControls', 
        # u'positiveControls', u'assayNegativeControls', u'outputFormat', 
        # u'outputFileName', u'assayPositiveControls', u'libraryControls', 
        # u'assayOtherControls', u'otherControls', u'plates'])
        # input_file_fields: 
        # set([u'readouts', u'replicates', u'readoutTypeSelection', 
        # u'conditions', u'collationOrderOrdering', u'uploadedFilename'])
        
        if i % 100 == 0:
            logger.info('processed %d activities, %d input_files', i, input_file_count)
    
    logger.info('stats: %d activities, input_files: %d, %d screens, %d cprs', 
        i, input_file_count, len(facility_ids),len(cpr_ids) )


    logger.info('transform_fields: %r', transform_fields)
    logger.info('input_file_fields: %r', input_file_fields)
    logger.info('output_types: %r', output_types)
    logger.info('assay_plate_sizes: %r', assay_plate_sizes)
    logger.info('readout_types:%r', readout_types)
    # readout_types:set([u'FP', u'LUMINESCENCE', u'UNSPECIFIED', u'IMAGING', 
    # u'FLUORESCENCE_INTENSITY', u'ABSORBANCE', u'PHOTOMETRY'])
    logger.info('collation_orderings: %r', collation_orderings)
    # set([u'P,R,C,R', u'C,P,Q,R,R', u'R,P,Q,C,R', u'P,Q,R,C,R', u'R,C,P,Q,R', 
    # u'P,C,R,R', u'C,R,P,R', u'C,R,P,Q,R', u'C,P,R,R', u'P,Q,C,R,R'])
    logger.info('readouts: %r', readouts)
    # set([u'', u'CrosstalkCorr, Lumin', u'CrossTalkCorr,Lumin', u'x,y,z', 
    # u'CrossTalkCorr\r\nL', u'CrossTalkCorr,L', u'Crosstalk, Lumin', 
    # u'Crosstalk_Corr, Lumin', u'FP, PChannel, SChannel', u'CrossTalkCorr,lumin', 
    # u'CrossTalkCorr,\r\nL', u'CrosstalkCorr, L', u'FP, P channel, S channel'])
    logger.info('replicates: %r', replicates)
    # set([32, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 24, 36, 30])
    logger.info('conditions: %r', conditions)
#     raise Exception('xxxx')
    
class Migration(migrations.Migration):

    dependencies = [
        ('db', '0018_library_screening_and_cpr_logs'),
    ]

    operations = [
        migrations.CreateModel(
            name='RawDataInputFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('ordinal', models.IntegerField()),
                ('collation_order', models.TextField(null=True)),
                ('conditions', models.TextField(null=True)),
                ('readouts', models.TextField(null=True)),
                ('readout_type', models.TextField(null=True)),
                ('replicates', models.IntegerField(null=True)),
                ('filename', models.TextField(null=True)),
            ],
            options={
                'db_table': 'raw_data_input_file',
            },
        ),
        migrations.CreateModel(
            name='RawDataTransform',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('plate_ranges', models.TextField(null=True)),
                ('output_filename', models.TextField(null=True)),
                ('temp_output_filename', models.TextField(null=True)),
                ('output_sheet_option', models.TextField(null=True)),
                ('assay_plate_size', models.TextField(null=True)),
                ('assay_positive_controls', models.TextField(null=True)),
                ('assay_negative_controls', models.TextField(null=True)),
                ('assay_other_controls', models.TextField(null=True)),
                ('library_controls', models.TextField(null=True)),
                ('library_plate_size', models.TextField(null=True)),
                ('comments', models.TextField(null=True)),
                ('cherry_pick_request', models.OneToOneField(null=True, to='db.CherryPickRequest')),
                ('screen', models.OneToOneField(null=True, to='db.Screen')),
            ],
            options={
                'db_table': 'raw_data_transform',
            },
        ),
        migrations.AddField(
            model_name='rawdatainputfile',
            name='raw_data_transform',
            field=models.ForeignKey(to='db.RawDataTransform'),
        ),

        migrations.RunPython(migrate_screen_transformer_logs),
        migrations.RunPython(migrate_cpr_transformer_logs),
    ]
