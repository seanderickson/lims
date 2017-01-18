# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from db import WELL_NAME_PATTERN
from db.models import CherryPickRequest
from db.support import lims_utils
import logging

logger = logging.getLogger(__name__)


def migrate_cherry_pick_request_empty_wells(apps,schema_editor):
    logger.info('migrate cherry_pick_request_empty_wells')

    CherryPickRequestModel = apps.get_model('db', 'CherryPickRequest')
    
    connection = schema_editor.connection
    
    sql = (
        'select well_name from cherry_pick_request_empty_well where '
        'cherry_pick_request_id=%s' )
    
    
    for cpr in CherryPickRequest.objects.all(): #.filter(cherry_pick_request_id=43882): 
        
        plate_size = cpr.assay_plate_size
        nrows = lims_utils.get_rows(plate_size)
        ncols = lims_utils.get_cols(plate_size)
        logger.info('cpr: %d, plate_size: %d nrows: %d, ncols: %d', 
            cpr.cherry_pick_request_id, plate_size, nrows, ncols)
        
        cursor = connection.cursor()
        cursor.execute(sql, [cpr.cherry_pick_request_id])
        _list = cursor.fetchall()
        
        well_names = set(x[0] for x in _list)
        wells_to_leave_empty = set()
        rowhash = {}
        colhash = {}
        wells = set()
        for wellname in well_names:
            match = WELL_NAME_PATTERN.match(wellname)
            if not match:
                logger.error('well name pattern is not recognized: %r', wellname)
                raise Exception('cpr %d', cpr.cherry_pick_request_id)
                continue
            wells.add(wellname)
            wellrow = match.group(1)
            wellcol = int(match.group(2))
            _row_wells = rowhash.get(wellrow,set())
            _row_wells.add(wellname)
            rowhash[wellrow] = _row_wells
            _col_wells = colhash.get(wellcol,set())
            _col_wells.add(wellname)
            colhash[wellcol] = _col_wells
        for _col, colwells in colhash.items():
            # Note: using >= nrows because of an error in extant data for 96 well plate types
            if len(colwells) >= nrows:
                wells_to_leave_empty.add('Col:%d' % _col)
        for _row, rowwells in rowhash.items():
            if len(rowwells) >= ncols:
                wells_to_leave_empty.add('Row:%s' % _row) 
                
#         logger.info('wells: %r', sorted(wells))
#         logger.info('well row cols: %r', wells_to_leave_empty)
#         logger.info('col hash: %r', colhash)
#         logger.info('row_hash: %r', rowhash)
        for wellname in wells:
            match = WELL_NAME_PATTERN.match(wellname)
            wellrow = match.group(1)
            wellcol = int(match.group(2))
            
            if 'Col:%d' % wellcol in wells_to_leave_empty:
                logger.debug('col wellname found: %r', wellname)
                continue
            if 'Row:%s' % wellrow in wells_to_leave_empty:
                logger.debug('row wellname found: %r', wellname)
                continue
            wells_to_leave_empty.add(wellname)
        logger.debug('wells_to_leave_empty: %r', wells_to_leave_empty)
        decorated = []
        for wellname in wells_to_leave_empty:
            if 'Col:' in wellname:
                decorated.append((1,int(wellname.split(':')[1]),wellname))
            elif 'Row:' in wellname:
                decorated.append((2,wellname.split(':')[1],wellname))
            else:
                match = WELL_NAME_PATTERN.match(wellname)
                decorated.append((match.group(1),match.group(1), wellname))
        wells_to_leave_empty = [x[2] for x in sorted(decorated)]
        
        cpr.wells_to_leave_empty = ','.join(wells_to_leave_empty)
        cpr.save()
        
        logger.info('cpr: %d, wells_to_leave_empty: %s', 
            cpr.cherry_pick_request_id, cpr.wells_to_leave_empty)



class Migration(migrations.Migration):
    # Post-migration after 0016:
    # remove old concentration fields from copy and plate
    
    dependencies = [
        ('db', '0017_copy_plate_well_post_migration'),
    ]

    
    operations = [
        migrations.AddField(
            model_name='cherrypickrequest',
            name='wells_to_leave_empty',
            field=models.TextField(null=True),
        ),
        migrations.RunPython(migrate_cherry_pick_request_empty_wells),
        
        # Set LabCherryPick.copy using well volume adjustments
        migrations.RunSQL(
            'update lab_cherry_pick set copy_id=c.copy_id '
            'from copy c join well_volume_adjustment wva using(copy_id) '
            'where lab_cherry_pick.lab_cherry_pick_id=wva.lab_cherry_pick_id;'),
            
    ]
