# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
from db import WELL_NAME_PATTERN
from db.models import CherryPickRequest
from db.support import lims_utils
import logging

from db.migrations import create_log_time, _create_generic_log,\
    _child_log_from

logger = logging.getLogger(__name__)

base_uri = '/db/api/v1'
cpap_resource_name = 'cherrypickassayplate'
cpr_resource_name = 'cherrypickrequest'
screen_resource_name = 'screen'
cpap_uri = '/db/api/v1/' + cpap_resource_name

def _create_cpr_log(cpr, activity):
    
    cpr_log = _create_generic_log(activity)
    cpr_log.ref_resource_name = cpr_resource_name
    cpr_log.key = str(cpr.cherry_pick_request_id)
    cpr_log.uri = '/'.join([
        base_uri,screen_resource_name,cpr.screen.facility_id,
        cpr_resource_name, cpr_log.key])
    cpr_log.api_action = 'PATCH'

    return cpr_log

def migrate_cherry_pick_request_empty_wells(apps,schema_editor):
    
    logger.info('migrate cherry_pick_request_empty_wells')

    CherryPickRequestModel = apps.get_model('db', 'CherryPickRequest')
    
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
        
        logger.debug('cpr: %d, wells_to_leave_empty: %s', 
            cpr.cherry_pick_request_id, cpr.wells_to_leave_empty)
        count += 1
        
    logger.info('migrated %d cpr.wells_to_leave_empty', count)

def migrate_plating_activity_logs(apps,schema_editor):
    logger.info('migrate_plating_activity_date')
    # TODO: replaces migration 0016
    # for each cpap, grab the cplt activity and create a new log for this
    
    # CPAP logs
    # CPR log
    
    # TODO
    
    
def create_lab_cherry_pick_logs(apps, schema_editor):
    
    # a. create logs for all the cplt's, and child cpap's, then finally for the
    # wvas on each cpap
    # strategy:
    # - iterate through all cplts, create uber-parent log.
    # - iterate cpaps, creating parent log,
    # - iterate wva's on cpap, create copywell volume change logs,
    # then look up the wva-lcp-cpap-cplt and find the log for it & set it as the parent log
    # for that cw log, lcp log
    # TODO: 
    #   - create a state variable on CPR "latest cplt activity" that 
    #   can be changed as part of the uber-parent log.  Note, not necessary 
    #   but may improve usability.  ( Otherwise, uber-parent log can just have
    #   comments, but no diff keys.)
    # Alternative - create a log for "plated plates" or "number of wells plated",
    #     Or, use the "last_plating_activity_date"
    # - create a state variable on cpap, change the state with a log (date plated)
    # 
    
    cpap_parent_logs = {}
    cpr_last_date_plated = {}
    
    CherryPickLiquidTransfer = apps.get_model('db', 'CherryPickLiquidTransfer')
    
    i = 0
    liquid_transfers = CherryPickLiquidTransfer.objects.all().order_by()
    logger.info('create logs for %d liquid_transfers', len(liquid_transfers))
    for cplt in liquid_transfers:
        lab_activity = cplt.labactivitylink
        activity = lab_activity.activitylink
        cpr = cplt.cherrypickassayplate_set.all()[0].cherry_pick_request
        screen_facility_id = lab_activity.screen.facility_id

        # 1. create a parent log on the cherry pick
        extra_information = {}
        cpr_log = _create_cpr_log(cpr, activity)
        
        # Note: "plating_activity" is a pseudo-key: this belies the need for 
        # an "activity" controlleed vocabulary for batch activities.
        
        # TODO: check: previous state will always be "not plated"
        # set the new cpap state
        previous_state = 'not_plated'
        cpap_state = 'not_plated'
        if cplt.status == 'Successful':
            cpap_state = 'plated'
        elif cplt.status == 'Failed':
            cpap_state = 'failed'
        elif cplt.status == 'Canceled':
            cpap_state = 'canceled'
        else:
            logger.warn('unknown state: %r, %r, %r', cplt.status,cpr,cplt)
        
        cpr_log.diffs = {}
#         if cpap_state == 'plated':
        prev_last_plating_activity_date = \
            cpr_last_date_plated.get(cpr.cherry_pick_request_id, None)
        cpr_log.diffs.update({
            'last_plating_activity_date': [
                prev_last_plating_activity_date, activity.date_of_activity]
            })
        prev_last_plating_activity_date[cpr.cherry_pick_request_id] = \
            activity.date_of_activity
#         else:
#             if cpap_state == 'failed':
#                   cpr_log.diffs.update({ 'status': cpap_state })  
#         cpr_log.diffs = {
#             'plating_activity': [previous_state, cpap_state]
#             }
        
        # store other activity information, in case needed
#         extra_information['date_of_activity'] = str(activity.date_of_activity)
#         extra_information['created_by_id'] = activity.created_by_id
        # n/a for cherrypick plates
        # extra_information['volume_transferred_per_well_from_library_plates'] = \
        #    lab_activity.volume_transferred_per_well_from_library_plates
#         extra_information['screen'] = screen_facility_id
#         cpr_log.json_field = json.dumps(extra_information)
        cpr_log.save()
        
        j = 0
        # TODO: should organize the logs by Copy Plate, not assay plate
        # i.e. store plate_number: 
        copy_plate_screening = {}
        
        for cpap in cplt.cherrypickassayplate_set.all():
            cpap_log = _child_log_from(cpr_log)
            
            cpap_log.ref_resource_name = "cherrypickassayplate"
            cpap_log.key = '/'.join(str(x) for x in [
                cpap.cherry_pick_request_id, 
                cpap.plate_ordinal, 
                cpap.attempt_ordinal ])
            cpap_log.uri = '/'.join([base_uri,cpap_log.ref_resource_name,cpap_log.key])
            
            # Note: previous state will always be "not plated"
            cpap_log.diffs = { 
                'plating_date': [None, activity.date_of_activity],
                'liquid_transfer_status': [previous_state,cpap_state] 
                }
            
            cpap_log.save()
            
            cpap_parent_logs[cpap.cherry_pick_assay_plate_id] = cpap_log
            j += 1
        logger.info('cpr: %r, cpaps processed %d',cpr.cherry_pick_request_id,j)
        
    logger.info(
        'finished step 1: %d, cpap parent logs: %d',len(liquid_transfers),
        len(cpap_parent_logs))
    create_lcp_wva_logs(apps,cpap_parent_logs, schema_editor)

def create_lcp_wva_logs(apps,cpap_parent_logs, schema_editor):
    
    logger.info('now create the child logs for all of the '
        'wvas on the cpap on the cplts' )
    
    CherryPickLiquidTransfer = apps.get_model('db', 'CherryPickLiquidTransfer')
    connection = schema_editor.connection
    cursor = connection.cursor()

    # create apilogs for cherry_pick_liquid_transfers
    # 1. well volume adjustments (copy_id, well_id)
    #    a. well volume correction activities
    #    - get date time from correction activity, who, comment
    #    b. lab_cherry_pick -> cherry_pick_assay_plate -> cherry_pick_liquid_transfer -> lab_activity
    #     or cherry_pick_screening -> cpap's -> lcp -> wva
    cols = OrderedDict({
        'well_id': 'well_id' ,
        'library_short_name': 'l.short_name',
        'copy_name':'c.name',
        'plate_number':'p.plate_number',
        'volume_adjustment': 'wva.volume',
        'initial_volume': 'p.well_volume',
        'comments': 'a.comments',
        'ecommons_id': 'u.ecommons_id',
        'email': 'u.email',
        'date_created': 'a.date_created',
        'date_of_activity': 'a.date_of_activity',
        'performed_by_id': 'performed_by_id',
        'login_id': 'u.login_id',
        'legacy_plate_name': 'cpap.legacy_plate_name',
        'cpap_id': 'cpap.cherry_pick_assay_plate_id'
        })
    _cols = ', '.join([ '%s as %s' % (value,key) for key, value in cols.items() ])
    # TODO: how to do a stored procedure here?
    query_sql = '\n'.join([
        'select ',
        _cols ,
        'from well_volume_adjustment wva ' ,
        'join copy c using(copy_id) ',
        'join library l using(library_id) '
        'join well w using(well_id) ',
        'join plate p on(wva.copy_id=p.copy_id and w.plate_number=p.plate_number) ',
        'join lab_cherry_pick lcp on (wva.lab_cherry_pick_id=lcp.lab_cherry_pick_id) ',
        'join cherry_pick_assay_plate cpap ',
        '  on (cpap.cherry_pick_assay_plate_id = lcp.cherry_pick_assay_plate_id)',
        'join cherry_pick_liquid_transfer cplt ', 
        '    on(cpap.cherry_pick_liquid_transfer_id=cplt.activity_id) '
        'join activity a on(a.activity_id = cplt.activity_id)',
        'join screensaver_user u on(u.screensaver_user_id=a.performed_by_id)' 
        ' where a.activity_id=%s ',
        'AND p.well_volume is not null ',
        'order by c.name, well_id,wva.well_volume_adjustment_id ',
        ])        
    logger.info('query_sql: %r', query_sql)

    liquid_transfers = CherryPickLiquidTransfer.objects.all()
    logger.info('create logs for %d liquid_transfers', len(liquid_transfers))
    
    count = 0
    parent_log_count = 0
    
    copywells = {}
    prev_logs = set()
    
    for xfer in liquid_transfers:
        cpap = xfer.cherrypickassayplate_set.all()[0] 
        
        if cpap.cherry_pick_assay_plate_id not in cpap_parent_logs:
            raise Exception('could not find a parent log for cpap %r, xfer %r' 
                % (cpap,xfer))
            
        parent_log = cpap_parent_logs[cpap.cherry_pick_assay_plate_id]
        
        parent_log_count += 1
        
        colkeys = cols.keys()
        # FIXME: create a log for each activity/cpap
        cursor.execute(query_sql, [xfer.labactivitylink.activitylink.activity_id])
        _list = cursor.fetchall()
        if len(_list) == 0:
            logger.info('no adjustments for %r, %r, %r', parent_log,xfer,xfer.status )
            continue;
    
        i = 0
        prev_volume = None
        for adjustment in _list:
            
            adj = dict(zip(colkeys, adjustment))
            
            log = _child_log_from(parent_log)
            try:
                
                log.ref_resource_name = copywell_resource_name
                log.api_action = 'PATCH'

                log.key = '/'.join([
                    adj['library_short_name'], adj['copy_name'], adj['well_id']])
                log.uri = copywell_uri + '/' + log.key
                log.json_field = str(adj['volume_adjustment'])
                
                # Short cut, just record the volume adjustment
                log.diffs = {'volume': ['',str(adj['volume_adjustment']) ]}
                # temporarily store the adjustment in the json field
                log.json_field = str(adj['volume_adjustment'])

                # Fixed: Manual migration 0016 uses all of the WVA's to set the copy well volume
                #  
                # # FIXME: because this list may not contain all wva's for the well,
                # # we must do a different, final pass over all wva's for the well 
                # # to construct the vol change.
                # log.diff_keys = json.dumps(['volume'])
                # if log.key in copywells:
                #     prev_volume = copywells[log.key]
                # else:
                #     prev_volume = round(adj['initial_volume'], 10)
                #     copywells[log.key] = prev_volume
                #  
                # new_volume = round(prev_volume + float(adj['volume_adjustment']),10)
                # log.diffs = json.dumps({ 
                #     'volume':[prev_volume, new_volume ]})
                # copywells[log.key] = new_volume
                
                log.comment = adj['comments'] or ''
                if  adj['legacy_plate_name']:
                    log.comment = log.comment + '. ' + adj['legacy_plate_name']

                if (log.ref_resource_name,log.key,log.date_time) in prev_logs :
                    logger.warn('log key already exists! %r', log)
                    log.date_time = log.date_time + timedelta(0,i+parent_log_count) # hack, add second
                
                log.save()
                prev_logs.add((log.ref_resource_name,log.key,log.date_time))
                count += 1
                i += 1
            except Exception, e:
                logger.exception('on lcp_wva log migrate')
                raise e
            if count % 10000 == 0:
                logger.info('parent_log_count: %d logs created: %d', 
                    parent_log_count,count )
#             if parent_log_count > 10: break
            
        logger.info('done, parent_log_count: %d, cplt count %d total: %d',
            parent_log_count, i, count)
        
    logger.info('total: %d, parent_log_count: %d',count,parent_log_count)

def create_well_correction_logs(apps, schema_editor):
    '''
    Scan all WVA's:
    - if WVA.comment matches a Cherry Pick ID, assign the WVA to the CP
    
    '''
    
    logger.info('create well correction logs...')
    
    cp_pattern = re.compile(r'.*cp(\d+)',re.I)
    
    i = 0
    total_corrections = 0
    query = apps.get_model('db', 'WellVolumeCorrectionActivity').objects.all()\
        .order_by('activity__activity__date_of_activity')
    for wvac in query:
        activity = wvac.activity.activity
        matched = False
        
        # TODO: can look at activity.lab_cherry_pick.cherry_pick
        
        if activity.comments:
            matched = cp_pattern.match(activity.comments)
            if matched:
                # for a cpr
                cpr_id = matched.group(1)
                logger.debug('locate cpr: %r', cpr_id)
                
                parent_log = None
                try:
                    cpr = apps.get_model('db', 'CherryPickRequest').objects.get(pk=cpr_id)  
                    logger.info('process correction activity for cpr: %r',cpr_id)
                    parent_log = _create_cpr_log(cpr, activity)
                    
                    parent_log.save()
                except Exception, e:
                    logger.info('could not find cpr: %r', cpr_id)
                    parent_log = _create_generic_log(activity)
                    
                    # FIXME: need a parent resource?
                    parent_log.ref_resource_name = copywell_resource_name
                    parent_log.key = copywell_resource_name
                    parent_log.uri = '/'.join([
                        base_uri, parent_log.ref_resource_name, parent_log.key])
                    
                    parent_log.save()
                
                # create wva logs
                j = 0
                for wva in wvac.wellvolumeadjustment_set.all():
                    log = _child_log_from(parent_log)

                    log.ref_resource_name = copywell_resource_name
                    log.key = '/'.join([
                        wva.copy.library.short_name,
                        wva.copy.name,
                        wva.well_id])
                    log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
                    
                    # Short cut, just record the volume adjustment
                    log.diffs = {'volume': [ None,str(wva.volume) ]}

                    # temporarily store the adjustment in the json field
                    log.json_field = str(wva.volume)
                    
                    log.save()
                    j += 1
                logger.info('processed %d', j)
                total_corrections += j
                
        if not matched:
            # this is a manual adjustment, create a generic parent log
            parent_log = _create_generic_log(activity)
            
            # FIXME: need a parent resource?
            parent_log.ref_resource_name = copywell_resource_name
            parent_log.key = copywell_resource_name
            parent_log.uri = '/'.join([
                base_uri, parent_log.ref_resource_name, parent_log.key])
            
            parent_log.save()
            logger.info('create log for manual wvca: %r', parent_log)
            j = 0
            for wva in wvac.wellvolumeadjustment_set.all():
                log = _child_log_from(parent_log)

                log.ref_resource_name = copywell_resource_name
                log.key = '/'.join([
                    wva.copy.library.short_name,
                    wva.copy.name,
                    wva.well_id])
                log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
            
                # Short cut, just record the volume adjustment
                log.diffs = {'volume': [ None,str(wva.volume) ]}

                # temporarily store the adjustment in the json field
                log.json_field = str(wva.volume)
                log .save()
                j += 1
            logger.info('processed %d', j)
            total_corrections += j
        i += 1       
        
#             if i>10: break
        
    # finally, case c. where the wva has no parent:
    # we know these are all for one copy: 
    # copy_id = 664, name = 'C', library_short_name = 'Human4 Duplexes'
    
    copy = apps.get_model('db', 'Copy').objects.get(pk=664)
    copy1 = apps.get_model('db', 'Copy').objects.get(pk=659)
    
    parent_log = ApiLog()
    parent_log.date_time = create_log_time(datetime.date(2000, 1, 1))
    parent_log.ref_resource_name = copywell_resource_name
    parent_log.key = copywell_resource_name
    parent_log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
    
    parent_log1 = ApiLog()
    parent_log1.date_time = create_log_time(datetime.date(2000, 1, 2))
    parent_log1.ref_resource_name = copywell_resource_name
    parent_log1.key = copywell_resource_name
    parent_log1.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
    
    # assign to Stewart Rudnicki
    parent_log.user_id = 767
    parent_log.username = 'sr50'
    parent_log.comment = 'Manual well volume correction activity with no log information'
    parent_log.save()
    
    parent_log1.user_id = 767
    parent_log1.username = 'sr50'
    parent_log1.comment = 'Manual well volume correction activity with no log information'
    parent_log1.save()
    
    query = apps.get_model('db', 'WellVolumeAdjustment').objects\
        .filter(lab_cherry_pick__isnull=True)\
        .filter(well_volume_correction_activity__isnull=True)\
        .order_by('well__well_id')
    j = 0
    for wva in query:
        if wva.copy_id not in [659,664]:
            raise Exception('manual wva for unknown copy: %r', wva.copy_id,[659,664])
        
        if wva.copy_id == copy.copy_id:
            log = _child_log_from(parent_log)
        else:
            log = _child_log_from(parent_log1)
            
        log.ref_resource_name = copywell_resource_name
        log.key = '/'.join([
            wva.copy.library.short_name,
            wva.copy.name,
            wva.well_id])
        log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
        
        # Short cut, just record the volume adjustment
        log.diffs = {'volume': [ None,str(wva.volume) ]}

        # temporarily store the adjustment in the json field
        log.json_field = str(wva.volume)
        log .save()
        j += 1
    logger.info('orphaned wvas processed: %d', j)
    total_corrections += j
    logger.info('done, processed activities %d, corrections %d', i, total_corrections)


class Migration(migrations.Migration):
    
    dependencies = [
        ('db', '0017_copy_plate_well_post_migration'),
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
            '= ids.cherry_pick_assay_plate_id;') 
#         migrations.RunPython(create_well_correction_logs),
#         migrations.RunPython(create_lab_cherry_pick_logs),
            
    ]
