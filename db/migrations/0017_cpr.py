
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging

from django.db import migrations, models

from db import WELL_NAME_PATTERN
from db.support import lims_utils


logger = logging.getLogger(__name__)

def migrate_cherry_pick_request_empty_wells(apps,schema_editor):
    
    logger.info('migrate cherry_pick_request_empty_wells')

    CherryPickRequest = apps.get_model('db', 'CherryPickRequest')
    
    connection = schema_editor.connection
    
    sql = (
        'select well_name from cherry_pick_request_empty_well where '
        'cherry_pick_request_id=%s' )
    
    count = 0
    for cpr in CherryPickRequest.objects.all(): #.filter(cherry_pick_request_id=43882): 
        
        plate_size = cpr.assay_plate_size
        nrows = lims_utils.get_rows(plate_size)
        ncols = lims_utils.get_cols(plate_size)
        logger.debug('cpr: %d, plate_size: %d nrows: %d, ncols: %d', 
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
        
        logger.debug('cpr: %d, wells_to_leave_empty: %s', 
            cpr.cherry_pick_request_id, cpr.wells_to_leave_empty)
        count += 1
        
    logger.info('migrated %d cpr.wells_to_leave_empty', count)

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0016_plates'),
    ]

    operations = [
        migrations.AddField(
            model_name='screenercherrypick',
            name='selected',
            field=models.NullBooleanField(),
        ),
        migrations.AddField(
            model_name='screenercherrypick',
            name='searched_well',
            field=models.ForeignKey(
                related_name='searched_screener_cherry_pick', 
                to='db.Well', null=True),
        ),
        migrations.AddField(
            model_name='labcherrypick',
            name='copy',
            field=models.ForeignKey(related_name='copy_lab_cherry_picks', 
                to='db.Copy', null=True),
        ),
        
        migrations.RunPython(migrate_cherry_pick_request_empty_wells),
        
        # Set LabCherryPick.copy using well volume adjustments
        migrations.RunSQL(
            'update lab_cherry_pick '
            'set copy_id=c.copy_id '
            'from copy c join well_volume_adjustment wva using(copy_id) '
            'where lab_cherry_pick.lab_cherry_pick_id=wva.lab_cherry_pick_id;'),
        migrations.RunSQL(
            'update screener_cherry_pick '
            'set searched_well_id=screened_well_id; '),
        migrations.RunSQL(
            'update screener_cherry_pick '
            'set selected=true; '),
            
        migrations.RunSQL(
            'update cherry_pick_request '
            'set date_volume_reserved = cplts.date_volume_reserved '
            'from ( '
            ' select cherry_pick_request_id, max(date_of_activity) as date_volume_reserved '
            ' from cherry_pick_assay_plate '
            ' join cherry_pick_liquid_transfer '
            ' on(activity_id=cherry_pick_liquid_transfer_id) '
            ' join activity using(activity_id) '
            ' group by cherry_pick_request_id) cplts '
            ' where cplts.cherry_pick_request_id '
            '   = cherry_pick_request.cherry_pick_request_id;'),
        # set null=Fals after update
        migrations.AlterField(
            model_name='screenercherrypick',
            name='searched_well',
            field=models.ForeignKey(
                related_name='searched_screener_cherry_pick', 
                to='db.Well', null=False),
        ),
        migrations.AddField(
            model_name='cherrypickassayplate',
            name='plating_date',
            field=models.DateField(null=True),
        ),
        migrations.AddField(
            model_name='cherrypickassayplate',
            name='plated_by',
            field=models.ForeignKey(
                related_name='plated_cherry_pick_plates', 
                to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='cherrypickassayplate',
            name='screened_by',
            field=models.ForeignKey(
                related_name='screened_cherry_pick_plates', 
                to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='cherrypickassayplate',
            name='screening_date',
            field=models.DateField(null=True),
        ),
        # deprecate created_by - replaced by logs
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='created_by',
            field=models.ForeignKey(
                related_name='created_cherry_pick', 
                to='db.ScreensaverUser', null=True),
        ),
        
        migrations.RunSQL(
            'update cherry_pick_assay_plate '
            'set plating_date=a.date_of_activity '
            'from activity a '
            'where a.activity_id=cherry_pick_assay_plate.cherry_pick_liquid_transfer_id;'),
        # NOTE: CherryPickPlating: use "activity.created_by" for the log user,
        # and for the cpap.performed_by field; "performed_by" should an admin,
        # but the legacy activities show the screener user as the "activity.performed_by"
        migrations.RunSQL('''
            update cherry_pick_assay_plate
            set plated_by_id=coalesce(a.created_by_id,a.performed_by_id) 
            from activity a 
            where a.activity_id=cherry_pick_assay_plate.cherry_pick_liquid_transfer_id;
        '''),
            
        # Adjust the plate_ordinal to be one's based
        migrations.RunSQL(
            'update cherry_pick_assay_plate '
            'set plate_ordinal = cherry_pick_assay_plate.plate_ordinal+1 '
            'from  '
            '( '
            '    select distinct(cherry_pick_request_id) '
            '    from cherry_pick_assay_plate '
            '    where plate_ordinal = 0 '
            '    order by cherry_pick_request_id ' 
            ') cprs join ( '
            '    select '
            '    cherry_pick_assay_plate_id, cherry_pick_request_id, plate_ordinal  '
            '    from cherry_pick_assay_plate  '
            '    order by cherry_pick_request_id, plate_ordinal desc '
            ') ids using(cherry_pick_request_id) '
            'where cherry_pick_assay_plate.cherry_pick_assay_plate_id '
            '= ids.cherry_pick_assay_plate_id;'),
        # NOTE: the legacy system allows for multiple screening activities;
        # use the last one to set the screening_date
        migrations.RunSQL('''
            update cherry_pick_assay_plate 
            set screening_date=(select max(date_of_activity) 
              from activity a 
              join cherry_pick_assay_plate_screening_link cpapsl 
                on(cherry_pick_screening_id=a.activity_id) 
              where cpapsl.cherry_pick_assay_plate_id
                =cherry_pick_assay_plate.cherry_pick_assay_plate_id);
        '''),
        # NOTE: for the Cherry Pick Screening, the "performed_by" field is a 
        # screen member
        migrations.RunSQL('''
            update cherry_pick_assay_plate 
            set screened_by_id=(
              select a.performed_by_id 
              from activity a 
              join cherry_pick_assay_plate_screening_link cpapsl
                on(cherry_pick_screening_id=a.activity_id) 
              where cpapsl.cherry_pick_assay_plate_id
                =cherry_pick_assay_plate.cherry_pick_assay_plate_id
              order by date_of_activity desc limit 1);
        '''),
    ]
    
