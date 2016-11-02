
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import datetime
import json
import logging

from django.db import migrations, models
from pytz import timezone
import pytz

from db.support.data_converter import default_converter
from reports.models import ApiLog
from collections import OrderedDict
import re
from decimal import Decimal


logger = logging.getLogger(__name__)

base_uri = '/db/api/v1'
copywell_resource_name = 'copywell'
copywell_uri = '/db/api/v1/' + copywell_resource_name
cpap_resource_name = 'cherrypickassayplate'
cpap_uri = '/db/api/v1/' + cpap_resource_name
librarycopyplate_resource_name = 'librarycopyplate'

# this is a workaround, because some activities have identical date_of_activity
times_seen = set()
# unique offset for the logs in this migration to avoid collisions
plate_vol_time_offset = 1111
def create_log_time(input_date):
    date_time = pytz.timezone('US/Eastern').localize(
        datetime.datetime.combine(input_date, datetime.datetime.min.time()))
    i = 0
    date_time += datetime.timedelta(0,plate_vol_time_offset)
    while date_time in times_seen:
        i += 1
        date_time += datetime.timedelta(0,i)
    times_seen.add(date_time)
    return date_time

def _create_generic_log(activity):
    
    log = ApiLog()

    log.comment = activity.comments
    log.date_time = create_log_time(activity.date_of_activity)
    log.username = activity.performed_by.username
    if not log.username:
        log.username = '%s: %s' % (
            activity.performed_by.screensaver_user_id,
            activity.performed_by.email)
    
    log.user_id = activity.performed_by.screensaver_user_id
    
    return log

def _create_cpr_log(cpr, activity):
    
    cpr_log = _create_generic_log(activity)
    cpr_log.ref_resource_name = 'cherrypickrequest'
    cpr_log.key = str(cpr.cherry_pick_request_id)
    cpr_log.uri = '/'.join([base_uri,cpr_log.ref_resource_name,cpr_log.key])
    cpr_log.api_action = 'PATCH'

    return cpr_log

def _child_log_from(parent_log):
    child_log = ApiLog()
    child_log.parent_log = parent_log
    child_log.username = parent_log.username
    child_log.user_id = parent_log.user_id
    child_log.date_time = parent_log.date_time
    child_log.api_action = parent_log.api_action
    child_log.comment = parent_log.comment
    return child_log

def create_library_screening_logs(apps, schema_editor):
    
    logger.info('create library screening logs')
    
    LibraryScreening = apps.get_model('db', 'LibraryScreening')
    screen_to_screening_count = {}
    copyplate_to_screening_count = {}
    copyplate_to_volume = {}
    i = 0 
    total_plate_logs = 0
    for screening in LibraryScreening.objects.all().order_by(
            'screeninglink__labactivitylink__activitylink__date_of_activity'):
        lab_activity = screening.screeninglink.labactivitylink
        
        screen = lab_activity.screen
        activity = lab_activity.activitylink
        # create parent logs:
        logger.debug('for screen: %r', screen.facility_id)
        
        # screen log for library screening
        screen_log = _create_generic_log(activity)
        screen_log.ref_resource_name = 'screen'
        screen_log.key = screen.facility_id
        screen_log.uri = '/'.join([base_uri,screen_log.ref_resource_name,screen_log.key])
         # TODO: the edit UI will set a "screenings" count
        screen_count = screen_to_screening_count.get(screen.facility_id, 0)
        screen_log.diffs = {'screening_count': [screen_count, screen_count+1] }
        screen_count +=1
        screen_to_screening_count[screen.facility_id] = screen_count
        # TODO: actually create the screening_count var on screen
        try:
            screen_log.save()
        except Exception,e :
            logger.exception(
                'exception on save of %r, %s', screen_log, screen_log.ref_resource_name)
            raise e
        j = 0
        cp_logs = []
        for assay_plate in screening.assayplate_set.all().filter(replicate_ordinal=0):
            plate = assay_plate.plate
            # copyplate log for library screening
            cp_log = _child_log_from(screen_log)
            cp_log.ref_resource_name = librarycopyplate_resource_name 
            cp_log.key = '/'.join([
                plate.copy.library.short_name, plate.copy.name, 
                str(plate.plate_number).zfill(5)])
            cp_log.uri = '/'.join([base_uri,cp_log.ref_resource_name,cp_log.key])
             # TODO: the edit UI will set a "screenings" count
            screen_count = copyplate_to_screening_count.get(cp_log.key, 0)
            old_volume = copyplate_to_volume.get(cp_log.key, plate.well_volume )
            if not old_volume:
                # not sure what these library plates with no initial well volume are -
                # but cannot compute if not known
                old_volume = 0
            adjustment = lab_activity.volume_transferred_per_well_from_library_plates
            if not adjustment:
                 # not sure what these "library screenings" with no volume are:
                 # -- some are external library plates, presumably
                 # -- some have no AP's and are "z prime" logs?
                 # -- some are legacy records from before ss 2.0
                adjustment = 0
            new_volume = old_volume-adjustment
            cp_log.diffs = {
                'screening_count': [screen_count, screen_count+1],
                'remaining_volume': [str(old_volume),str(new_volume)]
                 }
            screen_count +=1
            copyplate_to_screening_count[cp_log.key] = screen_count
            copyplate_to_volume[cp_log.key] = new_volume
            cp_log.json_field = json.dumps({
                'volume_transferred_per_well_from_library_plates': str(adjustment)
                })
            cp_logs.append(cp_log)
            j += 1
            # TODO: possibly create logs for the assay plates created (assaywells)
        # todo, use bulk create
        for cp_log in cp_logs:
            cp_log.save()

        i += 1
        total_plate_logs += j
        if total_plate_logs % 1000  == 0:
            logger.info('created screen logs: %d, plate logs: %d',i, total_plate_logs)
    logger.info('created screen logs: %d, plate logs: %d',i, total_plate_logs)
  
def create_plate_activity_logs(apps, schema_editor):  

    logger.info('create plate activity logs')

    Activity = apps.get_model('db', 'Activity')

    cols = OrderedDict({
        'activity_id': 'a.activity_id',
        'username': 'username',
        'screensaver_user_id': 'screensaver_user_id',
        'date_of_activity': 'date_of_activity',
        'comments': 'a.comments',
        'plate_number': 'plate_number',
        'copy_name': 'copy.name',
        'library_short_name': 'library.short_name',
        })
    colkeys = cols.keys()
    _cols = ', '.join([ '%s as %s' % (value,key) for key, value in cols.items() ])
    sql = (
        'select ' + _cols + 
    ''' from activity a join screensaver_user on(performed_by_id=screensaver_user_id) 
    join plate on (activity_id=plate.plated_activity_id)
    join copy using(copy_id)
    join library using(library_id); '''
    )
    
    connection = schema_editor.connection
    cursor = connection.cursor()
        
    cursor.execute(sql)
    _list = cursor.fetchall()
    if len(_list) == 0:
        raise Exception('No plate plated_activities found with sql: %r' % sql)
    for i,_data in enumerate(_list):
        _activity = dict(zip(colkeys, _data))
        log = ApiLog()
        log.ref_resource_name = librarycopyplate_resource_name
        log.key = '/'.join([
            _activity['library_short_name'],_activity['copy_name'],
            str(int(_activity['plate_number'])).zfill(5)])
        log.uri = '/'.join([base_uri,log.ref_resource_name,log.key])
        log.comment = _activity['comments']
        log.date_time = create_log_time(_activity['date_of_activity'])
        log.username = _activity['username']
        log.user_id = _activity['screensaver_user_id']
        if "'available'" in log.comment.lower():
            log.diffs = { 'status': ['not_specied','available']}
        elif "'not available'" in log.comment.lower():
            log.diffs = { 'status': ['not_specied','not_available']}
        else:
            raise Exception('unknown plate.plated_activity comment: %r', _activity)
        log.save()
        
        if i % 1000 == 0:
            logger.info('processed %d plate plated activity logs', i)
    logger.info('processed %d plate plated activity logs', i)

    sql = (
        'select ' + _cols + 
    ''' from activity a join screensaver_user on(performed_by_id=screensaver_user_id) 
    join plate on (activity_id=plate.retired_activity_id)
    join copy using(copy_id)
    join library using(library_id); '''
    )

    cursor.execute(sql)
    _list = cursor.fetchall()
    status_change_pattern = re.compile(r".*from '([^\']+)'.*to '([^\']+)'.*")
    if len(_list) == 0:
        raise Exception('No plate retired_activities found with sql: %r' % sql)
    
    status_terms_recognized = set()
    for i,_data in enumerate(_list):
        _activity = dict(zip(colkeys, _data))
        log = ApiLog()
        log.ref_resource_name = librarycopyplate_resource_name
        log.key = '/'.join([
            _activity['library_short_name'],_activity['copy_name'],
            str(int(_activity['plate_number'])).zfill(5)])
        log.uri = '/'.join([base_uri,log.ref_resource_name,log.key])
        log.comment = _activity['comments']
        log.date_time = create_log_time(_activity['date_of_activity'])
        log.username = _activity['username']
        log.user_id = _activity['screensaver_user_id']
        
        match = status_change_pattern.match(log.comment)
        if not match:
            raise Exception('unknown plate.retired_activity comment: %r', _activity)
        log.diffs = {'status': [
            default_converter(match.group(1)),
            default_converter(match.group(2))]}
        log.save()
        status_terms_recognized.add(default_converter(match.group(1)))
        status_terms_recognized.add(default_converter(match.group(2)))
        
        if i % 1000 == 0:
            logger.info('processed %d plate retired activity logs', i)
    logger.info('processed %d plate retired activity logs', i)
    logger.info('status terms recognized: %r', status_terms_recognized)


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
    # - create a state variable on cpap, change the state with a log
    # 
    
    cpap_parent_logs = {}
    
    CherryPickLiquidTransfer = apps.get_model('db', 'CherryPickLiquidTransfer')
    
    i = 0
    liquid_transfers = CherryPickLiquidTransfer.objects.all()
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
            
        cpr_log.diffs = {
            'plating_activity': [previous_state, cpap_state]
            }
        
        # store other activity information, in case needed
        extra_information['date_of_activity'] = str(activity.date_of_activity)
        extra_information['created_by_id'] = activity.created_by_id
        # n/a for cherrypick plates
        # extra_information['volume_transferred_per_well_from_library_plates'] = \
        #    lab_activity.volume_transferred_per_well_from_library_plates
        extra_information['screen'] = screen_facility_id
        cpr_log.json_field = json.dumps(extra_information)
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
            
            cpap_log.diffs = { 'state': [previous_state,cpap_state] }
            
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
    
    logger.info('create well correction logs...')
    
    cp_pattern = re.compile(r'.*cp(\d+)',re.I)
    
    i = 0
    total_corrections = 0
    query = apps.get_model('db', 'WellVolumeCorrectionActivity').objects.all()\
        .order_by('activity__activity__date_of_activity')
    for wvac in query:
        activity = wvac.activity.activity
        matched = False
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


#  This would be a method to make complete diff logs
##
## TODO: this is tricky - need to verify that the copy_well volumes (final)
## match the copy_well volumes from the manual sql migration 0016
##
def create_copywell_adjustments(apps, schema_editor):
    # Now go back over all of the corrections and construct the diffs
     
    # first get the plate volume = initial volume
    logger.info('create copywell adjustments...')
    librarycopyplate_pattern = re.compile(r'^[^\/]+\/[^\/]+\/\d+$')
    
    copy_plate_initial_volumes = {}
    for plate in apps.get_model('db', 'Plate').objects.all():
        key = ( '%s/%s/%s' 
            % ( plate.copy.library.short_name, 
                plate.copy.name,
                str(plate.plate_number).zfill(5)))
        copy_plate_initial_volumes[key] = plate.well_volume
         
    logger.info(
        'plate volume map built, iterate through copywells for %d plates...',
        len(copy_plate_initial_volumes))
    prev_wellcopy_key = None
    prev_plate_key = None
    initial_plate_vol = None
    current_wellcopy_volume = None
    j = 0
    for log in ( apps.get_model('reports', 'ApiLog').objects.all()
            .filter(ref_resource_name='copywell')
#             .filter(diff_keys__contains='volume')
            .order_by('key','date_time')):
        plate_key = log.key.split(':')[0]
        if not librarycopyplate_pattern.match(plate_key):
            # skip these logs, may be "parent" logs
            logger.info('non-copyplate log_key: %r', log.key)
            continue
        
        if plate_key != prev_plate_key:
            prev_plate_key = plate_key
            initial_plate_vol = copy_plate_initial_volumes[plate_key]
         
        if log.key != prev_wellcopy_key:
            current_wellcopy_volume = initial_plate_vol
            prev_wellcopy_key = log.key
        # TODO: would be better to make a hash than to store in the json_field
        new_volume = current_wellcopy_volume + Decimal(log.json_field)
        log.diffs = {
            'volume': [str(current_wellcopy_volume),
                       str(new_volume)],
           }
        current_wellcopy_volume = new_volume
        log.json_field = json.dumps({ 'volume_adjustment': log.json_field })
        log.save()
        j += 1
        if j % 10000 == 0:
            logger.info('processed: %d',j)
    logger.info('volume apilogs adjusted: %d', j)

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0015_create_lib_content_diffs_rnai'),
    ]

    # 1 (done). assay_plate / library screenings xfers:
    #    - assay_plate -> library_screening ->> lab_activity
    # library_screening -> assay_plate -> plate
    #                                    -> screen
    # Follows are preliminary-not yet implemented
    # # 2. create ApiLogs for WVA's: temporarily store volume adjusment on the log
    # 
    # # a. do cplt's
    # self.create_lab_cherry_pick_logs(orm) 
    #         
    # # b. do wvca's
    # # c. do orphaned wvas (no cplt or wvca attached)
    # self.create_well_correction_logs(orm)
    #  
    # # d. do a final pass, over the logs sorted by date, in order to 
    # #    construct the previous_volume, new volume in the log
    # self.create_copywell_adjustments(orm)
    
    operations = [
        migrations.RunPython(create_library_screening_logs),
        migrations.RunPython(create_lab_cherry_pick_logs),
        migrations.RunPython(create_plate_activity_logs),
        migrations.RunPython(create_well_correction_logs),
        migrations.RunPython(create_copywell_adjustments),
    ]
