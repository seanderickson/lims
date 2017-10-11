# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from collections import OrderedDict
from collections import defaultdict
import datetime
from decimal import Decimal
import json
import logging
import re

from django.db import migrations, models

from db.migrations import create_log_time, _create_generic_log, \
    _child_log_from
from db.support.data_converter import default_converter
from reports.models import ApiLog


logger = logging.getLogger(__name__)

base_uri = '/db/api/v1'
cpap_resource_name = 'cherrypickassayplate'
cpr_resource_name = 'cherrypickrequest'
plate_resource_name = 'librarycopyplate'
copywell_resource_name = 'copywell'
screen_resource_name = 'screen'
cpap_uri = '/db/api/v1/' + cpap_resource_name

copywell_resource_name = 'copywell'
copywell_uri = '/db/api/v1/' + copywell_resource_name
librarycopyplate_resource_name = 'librarycopyplate'
library_screening_resource_name = 'libraryscreening'


extant_plate_logs = defaultdict(set)
def find_extant_plate_logs(apps, schema_editor):
    # Build a hash of the current librarycopyplate/key to date_times;
    # Hack to avoid collisions caused by repeating times (between plate logs of 
    # different migrations
    resource_name = 'librarycopyplate'
    sql = '''
    select key, date_time
    from reports_apilog
    where ref_resource_name = %s
    '''
    with schema_editor.connection.cursor() as c:
        c.execute(sql,[resource_name])
        _row = c.fetchone()
        while _row:
            extant_plate_logs[_row[0]].add(_row[1])
            _row = c.fetchone()
            
        logger.info('extant_plate_logs: %d', len(extant_plate_logs))

extant_cpr_logs = defaultdict(set)
def find_extant_cpr_logs(apps, schema_editor):
    # Build a hash of the current librarycopyplate/key to date_times;
    # Hack to avoid collisions caused by repeating times (between plate logs of 
    # different migrations
    resource_name = 'cherrypickrequest'
    sql = '''
    select key, date_time
    from reports_apilog
    where ref_resource_name = %s
    '''
    
    with schema_editor.connection.cursor() as c:
        c.execute(sql,[resource_name])
        _row = c.fetchone()
        while _row:
            extant_cpr_logs[_row[0]].add(_row[1])
            _row = c.fetchone()
            
        logger.info('extant_cpr_logs: %d', len(extant_cpr_logs))

def _create_cpr_log(
    date_of_activity=None, cpr_id=None, screen_facility_id=None,
    username=None, email=None, performed_by_id=None, comments=None, **kwargs ):

    log = ApiLog()
    log.date_time = create_log_time(date_of_activity)
    if username:
        log.username = username
    else:
        if email:
            log.username = email
        else:
            logger.info(
                'cpr log w/o username or email: %r: %r', cpr_id,performed_by_id)
            log.username = performed_by_id
    log.user_id = performed_by_id
    if comments is not None:
        log.comment = comments
    else:
        log.comment = 'Cherry Pick Request reservation log (migration)'
    log.ref_resource_name = cpr_resource_name
    log.key = str(cpr_id)
    
    # Hack: to avoid integrity collisions between test migrations
    if log.date_time in extant_cpr_logs[log.key]:
        log.date_time = create_log_time(log.date_time)
    extant_cpr_logs[log.key].add(log.date_time)
        
    log.uri = '/'.join([
        'screen', screen_facility_id, log.ref_resource_name,log.key])
    log.api_action = 'PATCH'
    
    return log

def _create_wvac_log(
    date_of_activity=None, 
    username=None, email=None, performed_by_id=None, comments=None, **kwargs ):

    log = ApiLog()
    log.date_time = create_log_time(date_of_activity)
    if username:
        log.username = username
    else:
        log.username = email
    log.user_id = performed_by_id
    if comments is not None:
        log.comment = comments
    else:
        log.comment = 'Manual Well Volume Correction (migration)'
    log.ref_resource_name = copywell_resource_name
    log.key = copywell_resource_name
    log.uri = copywell_resource_name
    log.api_action = 'PATCH'
    
    return log

def _create_ls_log(
    date_of_activity=None, activity_id=None, screen_facility_id=None,
    performed_by_username=None, performed_by_id=None, 
    created_by_username=None, created_by_id=None, 
    comments=None, **kwargs ):

    log = ApiLog()
    log.ref_resource_name = library_screening_resource_name
    log.key = str(activity_id)
    log.uri = '/'.join([
        'screen', screen_facility_id, log.ref_resource_name,log.key])
    log.api_action = 'PATCH'
    log.date_time = create_log_time(date_of_activity)
    if created_by_username:
        log.username = created_by_username
    else:
        log.username = 'sr50'
    if created_by_id:
        log.user_id = created_by_id
    else:
        log.user_id = 767
    if comments is not None:
        log.comment = comments
    else:
        log.comment = 'Library Screening log (migration)'
    
    return log

def _create_plate_activity_log(activity_dict):
    log = ApiLog()
    log.ref_resource_name = librarycopyplate_resource_name
    log.api_action = 'PATCH'
    log.key = '/'.join([
        activity_dict['library_short_name'],activity_dict['copy_name'],
        str(int(activity_dict['plate_number']))])
#         str(int(activity_dict['plate_number'])).zfill(5)])
    log.uri = '/'.join([base_uri,log.ref_resource_name,log.key])
    log.comment = activity_dict['comments']
    log.date_time = create_log_time(activity_dict['date_of_activity'])
    if log.date_time in extant_plate_logs[log.key]:
        log.date_time = create_log_time(log.date_time)
    extant_plate_logs[log.key].add(log.date_time)    
    log.username = activity_dict['username']
    if log.username is None:
        log.username = activity_dict['email']
    log.user_id = activity_dict['screensaver_user_id']
    return log

def create_copy_comments_logs(apps, schema_editor):
    
    #TODO convert copy.comments to apilogs
    
    pass
    
def create_screen_comments_logs(apps, schema_editor):
    
    #TODO convert screen.comments to apilogs
    
    pass
    
def create_user_comments_logs(apps, schema_editor):
    
    #TODO convert screen.comments to apilogs
    
    pass
    

def create_plate_generic_logs(apps, schema_editor):
    logger.info('create plate generic activity logs')

    Activity = apps.get_model('db', 'Activity')

    cols = OrderedDict({
        'activity_id': 'a.activity_id',
        'username': 'username',
        'email': 'email',
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
        '''
    from activity a 
    join screensaver_user on(performed_by_id=screensaver_user_id)
    join plate_update_activity on (activity_id=update_activity_id)
    join administrative_activity aa using(activity_id)  
    join plate using(plate_id)
    join copy using(copy_id)
    join library using(library_id)
    where
    aa.administrative_activity_type = 'Comment' 
    and not exists (select null from plate where plated_activity_id = activity_id) 
    and not exists (select null from plate where retired_activity_id = activity_id) 
    and a.comments not ilike '%plate type changed%' 
    and a.comments not ilike '%Freezer 1 inventory updates for Richard Siu and Rachel Warden%'
    '''
    )
#     and a.comments not ilike '%activity created by database migration%'
#     and a.comments not ilike '%well volume changed%' 
#     and a.comments not ilike '%concentration changed from%' 
#     and a.comments not ilike '%concentration (molar) changed from%' 
#     and a.comments not ilike '%concentration (mg/ml) changed from%' 
    
    connection = schema_editor.connection
    cursor = connection.cursor()
    cursor.execute(sql)
    _list = cursor.fetchall()
    if len(_list) == 0:
        raise Exception('No plate generic activities found with sql: %r' % sql)
    for i,_data in enumerate(_list):
        _activity = dict(zip(colkeys, _data))
        log = _create_plate_activity_log(_activity)
        log.json_field = { 
            'migration': 'LibraryCopyPlate (comment log)',
            'data': { 'plate_update_activity.activity_id': 
                _activity['activity_id'] }
        }
        log.save()
        if i % 1000 == 0:
            logger.info('plate generic log: %r, %r', log, log.log_uri)
            logger.info('processed %d plate generic activity logs', i)
    
    logger.info('processed %d plate generic activity logs', i)
    
def create_plate_plated_and_retired_logs(apps, schema_editor):  

    logger.info('create plate activity logs')

    Activity = apps.get_model('db', 'Activity')

    cols = OrderedDict({
        'activity_id': 'a.activity_id',
        'username': 'username',
        'email': 'email',
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
    ''' 
        from activity a join screensaver_user on(performed_by_id=screensaver_user_id) 
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
        log = _create_plate_activity_log(_activity)
        if "'available'" in log.comment.lower():
            log.diffs['status'] = ['not_specied','available']
        elif "'not available'" in log.comment.lower():
            log.diffs['status'] = ['not_specied','not_available']
        else:
            raise Exception('unknown plate.plated_activity comment: %r', _activity)
        log.diffs['date_plated'] = [None, _activity['date_of_activity']]
        log.json_field = { 
            'migration': 'LibraryCopyPlate plated',
            'data': { 
                'plate_update_activity.activity_id': _activity['activity_id'],
                'comment': log.comment
            }
        }
        log.comment = None # update activity comment is unneeded
        log.save()
        
        if i % 1000 == 0:
            logger.info('plate plated log: %r, %r', log, log.log_uri)
            logger.info('processed %d plate plated activity logs', i)
    logger.info('processed %d plate plated activity logs', i)

    sql = (
        'select ' + _cols + 
    ''' 
        from activity a join screensaver_user on(performed_by_id=screensaver_user_id) 
        join plate on (activity_id=plate.retired_activity_id)
        join copy using(copy_id)
        join library using(library_id); '''
    )

    cursor.execute(sql)
    _list = cursor.fetchall()
    if len(_list) == 0:
        raise Exception('No plate retired_activities found with sql: %r' % sql)
    
    status_change_pattern = re.compile(r".*from '([^\']+)'.*to '([^\']+)'.*")
    status_terms_recognized = set()
    for i,_data in enumerate(_list):
        _activity = dict(zip(colkeys, _data))
        log = _create_plate_activity_log(_activity)
        match = status_change_pattern.match(log.comment)
        if not match:
            raise Exception('unknown plate.retired_activity comment: %r', _activity)
        log.diffs['status'] = [
            default_converter(match.group(1)),
            default_converter(match.group(2))]
        log.diffs['date_retired'] = [None, _activity['date_of_activity']]
        log.json_field = { 
            'migration': 'LibraryCopyPlate retired',
            'data': { 
                'plate_update_activity.activity_id': _activity['activity_id'],
                'comment': log.comment
            }
        }
        log.comment = None # retired activity comment is unneeded
        log.save()
        status_terms_recognized.add(default_converter(match.group(1)))
        status_terms_recognized.add(default_converter(match.group(2)))
        
        if i % 1000 == 0:
            logger.info('plate retired log: %r, %r', log, log.log_uri)
            logger.info('processed %d plate retired activity logs', i)
    logger.info('processed %d plate retired activity logs', i)
    logger.info('status terms recognized: %r', status_terms_recognized)

def create_library_screening_logs(apps, schema_editor):

    ls_cols = [
        'activity_id',
        'screen_facility_id',
        'plate_id',
        'plate_number',
        'initial_plate_well_volume',
        'copy_name',
        'library_short_name',
        'comments',
        'date_of_activity',
        'performed_by_username',
        'performed_by_id',
        'created_by_username',
        'created_by_id',
        'screened_experimental_well_count',
        'libraries_screened_count',
        'library_plates_screened_count',
        'volume_transferred_per_well_from_library_plates'
    ]
    ls_query_sql = '''
    select
    ls.activity_id,
    s.facility_id as screen_facility_id,
    p.plate_id,
    p.plate_number,
    p.well_volume as initial_plate_well_volume,
    c.name as copy_name,
    l.short_name as library_short_name,
    a.comments,
    a.date_of_activity,
    performed_by.username as performed_by_username,
    performed_by.screensaver_user_id as performed_by_id,
    created_by.username as created_by_username,
    created_by.screensaver_user_id as created_by_id,
    ls.screened_experimental_well_count,
    ls.libraries_screened_count,
    ls.library_plates_screened_count,
    la.volume_transferred_per_well_from_library_plates
    from assay_plate ap 
    join plate p using(plate_id)
    join copy c using(copy_id)
    join library l using(library_id)
    join library_screening ls on(library_screening_id=ls.activity_id)
    join lab_activity la on(ls.activity_id=la.activity_id)
    join screen s on (la.screen_id=s.screen_id)
    join activity a on(ls.activity_id=a.activity_id)
    left join screensaver_user performed_by on (a.performed_by_id=performed_by.screensaver_user_id)
    left join screensaver_user created_by on (a.created_by_id=created_by.screensaver_user_id)
    where ap.replicate_ordinal = 0
    order by date_of_activity, activity_id, plate_number
    '''
    library_screening_to_log = {}
    with schema_editor.connection.cursor() as c:
        c.execute(ls_query_sql)
        _row = c.fetchone()
        
        copyplate_to_screening_count = {}
        copyplate_to_volume = {}
        count = 1
        while _row:
            ls = dict(zip(ls_cols, _row))
            activity_id = ls['activity_id']
            screen_facility_id = ls['screen_facility_id']
            
            if activity_id not in library_screening_to_log:
                ls_log = _create_ls_log(**ls)
                ls_log.json_field = {
                    'migration': 'LibraryScreening',
                    'data': { 'library_screening.activity_id': activity_id }
                }
                ls_log.save()
                library_screening_to_log[activity_id] = ls_log
            # Create plate log
            cp_log = _child_log_from(ls_log)
            cp_log.comment = None
            cp_log.ref_resource_name = librarycopyplate_resource_name 
            cp_log.key = '/'.join([
                ls['library_short_name'], ls['copy_name'], 
                str(ls['plate_number'])])
#                 str(ls['plate_number']).zfill(5)])
            cp_log.uri = '/'.join([base_uri,cp_log.ref_resource_name,cp_log.key])
            if cp_log.date_time in extant_plate_logs[cp_log.key]:
                cp_log.date_time = create_log_time(cp_log.date_time)
            extant_plate_logs[cp_log.key].add(cp_log.date_time)    
            previous_screening_count = copyplate_to_screening_count.get(
                ls['plate_id'], 0)
            previous_volume = copyplate_to_volume.get(
                ls['plate_id'], ls['initial_plate_well_volume'] )
            if previous_volume is not None:
                previous_volume = Decimal(previous_volume)
            adjustment = ls['volume_transferred_per_well_from_library_plates']
            if not adjustment:
                 # not sure what these "library screenings" with no volume are:
                 # -- some are external library plates, presumably
                 # -- some have no AP's and are "z prime" logs?
                 # -- some are legacy records from before ss 2.0
                adjustment = Decimal(0)
            else:
                adjustment = Decimal(adjustment)
            
            cp_log.diffs['screening_count'] = [
                previous_screening_count, previous_screening_count+1]
            copyplate_to_screening_count[ls['plate_id']] = previous_screening_count +1

            # NOTE: some library plates have no well volume set: skip these
            if previous_volume is not None and ls['initial_plate_well_volume'] != 0:
                new_volume = previous_volume-adjustment
                cp_log.diffs['remaining_volume'] = [
                    str(previous_volume),str(new_volume)]
                copyplate_to_volume[ls['plate_id']] = new_volume
            
            cp_log.json_field = {
                'migration': 'LibraryScreening',
                'data': {
                    'volume_transferred_per_well_from_library_plates': 
                        str(adjustment)
                }
            }
            cp_log.save()
            
            _row = c.fetchone()
            count += 1
            if count % 10000  == 0:
                logger.info('library screening log: %r, %r', 
                    ls_log, ls_log.log_uri)
                logger.info('cp_log: %r, %r', cp_log, cp_log.log_uri)
                logger.info('created library screening logs: %d', count)
            
        logger.info('created library screening logs: %d', count)

    LibraryScreening = apps.get_model('db', 'LibraryScreening')
     
    for i,screening in enumerate( LibraryScreening.objects.all()
        .order_by('screeninglink__labactivitylink__activitylink__date_of_activity')):
        if not screening.assayplate_set.exists():
            continue
        lab_activity = screening.screeninglink.labactivitylink
        activity = lab_activity.activitylink
         
        if activity.activity_id not in library_screening_to_log:
            raise Exception(
                'library_screening not found: %r' % activity.activity_id)
        ls_log = library_screening_to_log[activity.activity_id]
        activity.apilog_uri = ls_log.log_uri
        activity.save()
         
        if i % 10000 == 0:
            logger.info('updated %d ls_log URIs', i)
         
    logger.info('updated %d ls_log URIs', i)

def create_well_volume_adjustment_logs(apps, schema_editor):
    '''
    Migrate logs for well volume adjustments:
    1. For lab cherry picks
    - new API CherryPickRequest.reserve logs:
    - CherryPickRequest (parent) log: date_volume_reserved, number (assay) plates
    - CopyWell logs (cherry_pick_screening_count, volume)
    - CopyPlate logs (cplt_screening_count)
    - TODO: Cherry Pick Assay Plate created logs? (TODO for API also)
    2. For manual well volume correction activities
    '''
    CherryPickLiquidTransfer = apps.get_model('db', 'CherryPickLiquidTransfer')
    CherryPickRequest = apps.get_model('db', 'CherryPickRequest')

    # Strategy:
    # build hash cpr-> plates
    # build hash copywell -> last_volume hash
    # build hash copywell -> cherry_pick_screening_count hash
    connection = schema_editor.connection

    plates_query_sql = '''
        with cpr_plates as (
        select
        lcp.cherry_pick_request_id,
        plate_id
        from lab_cherry_pick lcp
        join well on(well_id=source_well_id)
        join plate on(well.plate_number=plate.plate_number and plate.copy_id=lcp.copy_id)
        group by lcp.cherry_pick_request_id, plate_id
        )
        select
        cherry_pick_request_id,
        l.short_name,
        c.name,
        p.plate_number
        from 
        cpr_plates
        join plate p using (plate_id)
        join copy c using (copy_id)
        join library l using(library_id)
        order by cherry_pick_request_id, l.short_name, c.name,p.plate_number;
    '''
    cpr_plate_hash = defaultdict(list)
    
    with connection.cursor() as c:
        c.execute(plates_query_sql)
        _row = c.fetchone()
        while _row:
            cpr_plate_hash[_row[0]].append(_row[1:])
            _row = c.fetchone()
    if len(cpr_plate_hash) == 0:
        raise Exception('no cpr plates found')
    logger.info('cpr plates hash[0]: %r', cpr_plate_hash.items()[0])
    
    wva_cols = [
        'well_volume_adjustment_id',
        'well_id',
        'library_short_name',
        'copy_name',
        'plate_number',
        'volume_adjustment',
        'plate_initial_volume',
        'comments',
        'username',
        'email',
        'performed_by_id',
        'date_of_activity',
        'cpr_id',
        'cplt_id',
        'wvac_id',
        'screen_facility_id',
        'assay_plate_count']
    
    wva_query_sql = '''
        select * from (
        select
            well_volume_adjustment_id,
            well_id,
            l.short_name as library_short_name,
            c.name as copy_name,
            p.plate_number,
            wva.volume as volume_adjustment,
            p.well_volume as plate_initial_volume,
            a.comments,
            u.username,
            u.email,
            a.performed_by_id,
            a.date_of_activity,
            cpap.cherry_pick_request_id as cpr_id,
            cplt.activity_id as cplt_id,
            null as wvac_id,
            screen.facility_id as screen_facility_id,
            (select count(*) from 
                cherry_pick_assay_plate cpap1 
                where cpap1.cherry_pick_request_id=cpr.cherry_pick_request_id)::INTEGER as assay_plate_count
            from well_volume_adjustment wva
            join copy c using(copy_id) 
            join library l using(library_id) 
            join well w using(well_id) 
            join plate p on(wva.copy_id=p.copy_id and w.plate_number=p.plate_number)
            join lab_cherry_pick lcp on (wva.lab_cherry_pick_id=lcp.lab_cherry_pick_id)
            join cherry_pick_assay_plate cpap
                on (cpap.cherry_pick_assay_plate_id = lcp.cherry_pick_assay_plate_id)
            join cherry_pick_request cpr on(cpap.cherry_pick_request_id=cpr.cherry_pick_request_id)
            join screen on (cpr.screen_id=screen.screen_id)
            join cherry_pick_liquid_transfer cplt
                on(cpap.cherry_pick_liquid_transfer_id=cplt.activity_id) 
            join activity a on(a.activity_id = cplt.activity_id)
            join screensaver_user u on(u.screensaver_user_id=a.performed_by_id) 
        union
        select
            well_volume_adjustment_id,
            well_id,
            l.short_name as library_short_name,
            c.name as copy_name,
            p.plate_number,
            wva.volume as volume_adjustment,
            p.well_volume as plate_initial_volume,
            a.comments,
            u.username,
            u.email,
            a.performed_by_id,
            a.date_of_activity,
            null as cpr_id,
            null as cplt_id,
            well_volume_correction_activity_id as wvac_id,
            null as screen_facility_id,
            null as assay_plate_count
            from well_volume_adjustment wva
            join copy c using(copy_id) 
            join library l using(library_id) 
            join well w using(well_id) 
            join plate p on(wva.copy_id=p.copy_id and w.plate_number=p.plate_number)
            join activity a on(a.activity_id = well_volume_correction_activity_id)
            join screensaver_user u on(u.screensaver_user_id=a.performed_by_id)
        union
        select
            well_volume_adjustment_id,
            well_id,
            l.short_name as library_short_name,
            c.name as copy_name,
            p.plate_number,
            wva.volume as volume_adjustment,
            p.well_volume as plate_initial_volume,
            null as comments,
            null as username,
            null as email,
            null as performed_by_id,
            null as date_of_activity,
            null as cpr_id,
            null as cplt_id,
            null as wvac_id,
            null as screen_facility_id,
            null as assay_plate_count
            from well_volume_adjustment wva
            join copy c using(copy_id) 
            join library l using(library_id) 
            join well w using(well_id) 
            join plate p on(wva.copy_id=p.copy_id and w.plate_number=p.plate_number)
            where well_volume_correction_activity_id is null
            and lab_cherry_pick_id is null
        ) as a
        order by date_of_activity asc nulls first, well_volume_adjustment_id            
    '''
    cpr_logs = {}
    cpr_plate_logs = {}
    wva_correction_logs = {}
    manual_wvac_parent_log = _create_wvac_log(
        date_of_activity=datetime.date(2000, 1, 1), 
        username='sr50', 
        performed_by_id=767, 
        comments=('Migration: manual well volume correction '
            'activity with no log information'))
    manual_wvac_parent_log.json_field = { 
        'migration': 'WellVolumeCorrectionActivity (manual)(orphaned)' }
    manual_wvac_parent_log.save()
    cw_prev_data = {}
    plate_prev_data = {}
    cw_reserve_count = 0
    wvac_count = 0
    with connection.cursor() as c:
        c.execute(wva_query_sql)
        _row = c.fetchone()
        while _row:
            wva = dict(zip(wva_cols, _row))
            cpr_id = wva['cpr_id']
            wvac_id = wva['wvac_id']
            if cpr_id is not None:
                cpr_parent_log = cpr_logs.get(cpr_id, None)
                if cpr_parent_log is None:
                    cpr_parent_log = _create_cpr_log(**wva)
                    cpr_parent_log.diffs = {
                        'date_volume_reserved': 
                            [None, wva['date_of_activity'] ],
                        'number_plates': [0, wva['assay_plate_count']]
                    }
                    cpr_parent_log.json_field = { 
                        'migration': 'CherryPickRequest.reservation',
                        'data': { 'cherry_pick_liquid_transfer.activity_id': 
                            wva['cplt_id'] }
                    }
                    cpr_parent_log.save()
                    cpr_logs[cpr_id] = cpr_parent_log
    
                plate_key = '/'.join([
                        wva['library_short_name'], wva['copy_name'],str(wva['plate_number'])])
                cpr_plate_key = '%s/%s' % (cpr_id,plate_key)
                plate_log = cpr_plate_logs.get(cpr_plate_key, None) 
                if plate_log is None:               
                    # Create plate log
                    # Note: only one plate log per cpr/activity/plate
                    plate_log = _child_log_from(cpr_parent_log)

                    # Hack: to avoid integrity collisions between plate migrations
                    if plate_log.date_time in extant_plate_logs[plate_log.key]:
                        plate_log.date_time = create_log_time(plate_log.date_time)
                    extant_plate_logs[plate_log.key].add(plate_log.date_time)
                    plate_log.comment = None
                    plate_log.ref_resource_name = plate_resource_name
                    plate_log.key = plate_key
                    plate_log.uri = '/'.join([plate_log.ref_resource_name,plate_log.key])
                    prev_cplt_screening_count = plate_prev_data.get(plate_log.key, None)
                    if prev_cplt_screening_count is None:
                        prev_cplt_screening_count = 0
                    new_plate_cplt_screening_count = prev_cplt_screening_count+1
                    plate_prev_data[plate_log.key] = new_plate_cplt_screening_count
                    plate_log.diffs = {
                        'cplt_screening_count': [
                            prev_cplt_screening_count,
                            new_plate_cplt_screening_count ],
                        }
                    plate_log.save()
                    logger.debug('plate log: %r, %r', plate_log, plate_log.diffs)
                    cpr_plate_logs[cpr_plate_key] = plate_log
                cw_log = _child_log_from(cpr_parent_log)
                # NOTE: prevent integrity collisions with wvac_logs:
                cw_log.date_time = create_log_time(wva['date_of_activity'])
                cw_log.comment = None
                cw_log.ref_resource_name = copywell_resource_name
                cw_log.key = '/'.join([
                    wva['library_short_name'], wva['copy_name'], wva['well_id']])
                cw_log.uri = cw_log.ref_resource_name + '/' + cw_log.key
                cw_log.api_action = 'PATCH'
                prev_data = cw_prev_data.get(cw_log.key,None)
                if prev_data is None:
                    prev_data = [wva['plate_initial_volume'], 0]
                prev_volume = prev_data[0]
                new_volume = Decimal(prev_volume)+Decimal(wva['volume_adjustment'])
                prev_cp_screening_count = prev_data[1]
                new_cp_screening_count = prev_cp_screening_count+1
                cw_prev_data[cw_log.key] = [new_volume, new_cp_screening_count]
                
                cw_log.diffs = {
                    'volume': [str(prev_volume),str(new_volume)],
                    'cherry_pick_screening_count': [
                        prev_cp_screening_count,new_cp_screening_count]}
                # temporarily store the adjustment in the json field
                cw_log.json_field = { 
                    'migration': 'Cherry Pick Reservation (migration)',
                    'data': {'volume_adjustment': str(wva['volume_adjustment']) }
                }
                cw_log.comment = None
                cw_log.save()
                cw_reserve_count += 1
                if cw_reserve_count % 1000 == 0:
                    logger.info('cpr_parent_log: %r', cpr_parent_log.log_uri)
                    logger.info('cw_reserve_count %d', cw_reserve_count)
                    
            elif wvac_id is not None:
                logger.debug('wva: %r', wva)
                wvac_parent_log = wva_correction_logs.get(wvac_id, None)
                if wvac_parent_log is None:
                    wvac_parent_log = _create_wvac_log(**wva)
                    wvac_parent_log.json_field = { 
                        'migration': 'WellVolumeAdjustment (manual)',
                        'data': { 'well_volume_correction.activity_id': wvac_id }
                    }
                    wvac_parent_log.save()
                    wva_correction_logs[wvac_id] = wvac_parent_log
            
                cw_log = _child_log_from(wvac_parent_log)
                # NOTE: prevent integrity collisions with cpr_logs:
                cw_log.date_time = create_log_time(wva['date_of_activity'])
                cw_log.ref_resource_name = copywell_resource_name
                cw_log.key = '/'.join([
                    wva['library_short_name'], wva['copy_name'], wva['well_id']])
                cw_log.uri = cw_log.ref_resource_name + '/' + cw_log.key
                cw_log.api_action = 'PATCH'
                prev_data = cw_prev_data.get(cw_log.key,None)
                if prev_data is None:
                    prev_data = [wva['plate_initial_volume'], 0]
                prev_volume = prev_data[0]
                new_volume = Decimal(prev_volume)+Decimal(wva['volume_adjustment'])
                cw_prev_data[cw_log.key] = [new_volume, prev_data[1]]
                
                cw_log.diffs = {
                    'volume': [str(prev_volume),str(new_volume)] }
                # temporarily store the adjustment in the json field
                cw_log.json_field = { 
                    'migration': 'WellVolumeAdjustment (manual)',
                }
                cw_log.save()
                wvac_count += 1
                if wvac_count % 1000 == 0:
                    logger.info('wvac_parent_log: %r', wvac_parent_log.log_uri)
                    logger.info('wvac_count %d', wvac_count)
            else:
                logger.info('no cpr or wvac found: %r', wva)
                known_orphaned_wva_copies = ['Human4 Duplexes/C', 'Human1 Duplexes/D']                
                wva_copy = '{library_short_name}/{copy_name}'.format(**wva)
                if wva_copy not in known_orphaned_wva_copies:
                    raise Exception('unknown wva: %r, copy not recognized: %r')

                cw_log = _child_log_from(manual_wvac_parent_log)
                # NOTE: prevent integrity collisions with cpr_logs:
#                 cw_log.date_time = create_log_time(cw_log.date_time)
                cw_log.ref_resource_name = copywell_resource_name
                cw_log.key = '/'.join([
                    wva['library_short_name'], wva['copy_name'], wva['well_id']])
                cw_log.uri = cw_log.ref_resource_name + '/' + cw_log.key
                cw_log.api_action = 'PATCH'
                prev_data = cw_prev_data.get(cw_log.key,None)
                if prev_data is None:
                    prev_data = [wva['plate_initial_volume'], 0]
                prev_volume = prev_data[0]
                new_volume = Decimal(prev_volume)+Decimal(wva['volume_adjustment'])
                cw_prev_data[cw_log.key] = [new_volume, prev_data[1]]
                
                cw_log.diffs = {
                    'volume': [str(prev_volume),str(new_volume)] }
                # temporarily store the adjustment in the json field
                cw_log.json_field = { 'volume_adjustment': str(wva['volume_adjustment']) }
                cw_log.save()
                wvac_count += 1
                
            _row = c.fetchone()
        logger.info('cw_reserve_count %d', cw_reserve_count)
        logger.info('wvac_count %d', wvac_count)
                
def create_cherry_pick_plating_logs(apps, schema_editor):
    '''
    Create logs for the plating (liquid transfers)
    - parent_log: Cherry Pick Request: 
        last_plating_activity_date
    - child logs: for each Cherry Pick Assay Plate
    NOTE: keep CherryPickLiquidTransfer entries -> LabActivity (for now)
   '''
    CherryPickLiquidTransfer = apps.get_model('db', 'CherryPickLiquidTransfer')
    
    i = 0
    liquid_transfers = CherryPickLiquidTransfer.objects.all().order_by()
    logger.info('create logs for %d liquid_transfers', len(liquid_transfers))
    cpr_parent_logs = {}
    for cplt in liquid_transfers:
        lab_activity = cplt.labactivitylink
        activity = lab_activity.activitylink
        cpr = cplt.cherrypickassayplate_set.all()[0].cherry_pick_request
        screen_facility_id = lab_activity.screen.facility_id
        # NOTE: CherryPickPlating: use "activity.created_by" for the log user,
        # and for the cpap.performed_by field; "performed_by" should an admin,
        # but the legacy activities show the screener user as the "activity.performed_by"
        screensaver_user_legacy_performed_by = activity.performed_by
        screensaver_user_admin = activity.created_by
        
        # NOTE: if activity.date_performed < 2010-01-15, the "created_by" field
        # is empty
        if screensaver_user_admin is None:
            screensaver_user_admin = screensaver_user_legacy_performed_by
            logger.debug('Activity has no "created_by" id %r', activity.activity_id)
        
        cpr_parent_log = cpr_parent_logs.get(cpr.cherry_pick_request_id, None)
        if cpr_parent_log is None:
            cpr_parent_log = _create_cpr_log(
                date_of_activity=activity.date_of_activity, 
                cpr_id=cpr.cherry_pick_request_id, 
                screen_facility_id=screen_facility_id, 
                username=screensaver_user_admin.username, 
                email=screensaver_user_admin.email, 
                performed_by_id=screensaver_user_admin.screensaver_user_id, 
                comments=activity.comments)
            cpr_parent_log.diffs = {
                'last_plating_activity_date': 
                    [None, activity.date_of_activity ],
            }
            cpr_parent_log.json_field = { 
                'migration': 'CherryPickAssayPlate.plating_date',
                'data': { 
                    'cherry_pick_liquid_transfer.activity_id': activity.activity_id,
                    'status': cplt.status 
                }
            }
            cpr_parent_log.save()
            cpr_parent_logs[cpr.cherry_pick_request_id] = cpr_parent_log
            
        activity.apilog_uri = cpr_parent_log.log_uri
        activity.save()
            
        # NOTE: hack: use the json_field as a flag to find the cplt's that were
        # with more than one entry for a cpr: in this case the status has changed
        # (e.g. from "Failed" to "Canceled" or "Success")
        if cpr_parent_log.json_field is not None:
            prev_status = cpr_parent_log.json_field['data']['status']
            if prev_status != cplt.status:
                cpr_parent_log = _create_cpr_log(
                    date_of_activity=activity.date_of_activity, 
                    cpr_id=cpr.cherry_pick_request_id, 
                    screen_facility_id=screen_facility_id, 
                    username=screensaver_user_admin.username, 
                    email=screensaver_user_admin.email, 
                    performed_by_id=screensaver_user_admin.screensaver_user_id, 
                    comments=activity.comments)
                cpr_parent_log.diffs = {
                    'last_plating_activity_date': 
                        [None, activity.date_of_activity ],
                }
                cpr_parent_log.json_field = { 
                    'migration': 'CherryPickAssayPlate.plating_date',
                    'data': { 
                        'cherry_pick_liquid_transfer.activity_id': activity.activity_id,
                        'status': cplt.status 
                    }
                }
                cpr_parent_log.save()
                cpr_parent_logs[cpr.cherry_pick_request_id] = cpr_parent_log
                    
        for cpap in cplt.cherrypickassayplate_set.all():
            logger.debug('cpap: %r, %r, %r, %r, %r', 
                cpap.cherry_pick_request_id, cpap.plate_ordinal, cplt.status,
                activity.activity_id,
                cpr_parent_log.log_uri)  
            cpap_log = _child_log_from(cpr_parent_log)
            cpap_log.ref_resource_name = cpap_resource_name
            cpap_log.key = '/'.join(str(x) for x in [
                cpap.cherry_pick_request_id, 
                cpap.plate_ordinal ])
            cpap_log.uri = '/'.join([base_uri,cpap_log.ref_resource_name,cpap_log.key])
            # Note: previous state will always be "not plated"
            cpap_log.comment = (
                'Cherry Pick Plating (migration); performed_by: %s %s (%s)' 
                    % (screensaver_user_legacy_performed_by.first_name,
                       screensaver_user_legacy_performed_by.last_name,
                       screensaver_user_legacy_performed_by.username))
            cpap_log.diffs = { 
                'plating_date': [None, activity.date_of_activity],
#                 'liquid_transfer_status': [None, cplt.status],
                'plated_by_username': [None, screensaver_user_admin.username],
                'plated_by_name': [
                    None, '%s %s' 
                        % (screensaver_user_admin.first_name, 
                           screensaver_user_admin.last_name)],
                'plating_comments': [None, activity.comments]
                }
            
            cpap_log.save()
            i += 1
            
            if i % 100 == 0:
                logger.info('cpap: %r, %r, %r, %r, %r', 
                    cpap.cherry_pick_request_id, cpap.plate_ordinal, cplt.status,
                    activity.activity_id,
                    cpr_parent_log.log_uri)  
                logger.info('created %d cherry pick assay plate plating logs', i)
    logger.info('created %d cherry pick assay plate plating logs', i)
    
def create_cherry_pick_screening_logs(apps, schema_editor):
    '''
    Create logs for the screening
    - parent_log: Cherry Pick Request: 
        last_screening_activity_date, is_completed?, 
    - child logs: for each Cherry Pick Assay Plate
    Note: keep CherryPickScreening entries (for now)
   '''
    CherryPickScreening = apps.get_model('db', 'CherryPickScreening')
    CherryPickAssayPlate = apps.get_model('db', 'CherryPickAssayPlate')
    cpap_to_cps_sql = '''
        select cherry_pick_assay_plate_id
        from cherry_pick_assay_plate_screening_link
        where cherry_pick_screening_id = %s '''
    
    cursor = schema_editor.connection.cursor()
    
    i = 0
    cp_screenings = CherryPickScreening.objects.all().order_by()
    logger.info('create logs for %d CherryPickScreening', len(cp_screenings))
    cpr_parent_logs = {}
    for cp_screening in cp_screenings:
        
        lab_activity = cp_screening.screeninglink.labactivitylink
        activity = lab_activity.activitylink
        cpr = cp_screening.cherry_pick_request
        screen_facility_id = lab_activity.screen.facility_id
        screensaver_user = activity.performed_by
        screensaver_user_admin = activity.created_by
        
        # NOTE: if activity.date_performed < 2010-01-15, the "created_by" field
        # is empty
        if screensaver_user_admin is None:
            screensaver_user_admin = screensaver_user
            logger.debug('Activity has no "created_by" id %r', activity.activity_id)
        
        cpr_parent_log = cpr_parent_logs.get(cpr.cherry_pick_request_id, None)
        if cpr_parent_log is None:
            # NOTE: for the Cherry Pick Screening, the "performed_by" field is a 
            # screen member; but the ApiLog user will be the admin
            cpr_parent_log = _create_cpr_log(
                date_of_activity=activity.date_of_activity, 
                cpr_id=cpr.cherry_pick_request_id, 
                screen_facility_id=screen_facility_id, 
                username=screensaver_user_admin.username, 
                email=screensaver_user_admin.email, 
                performed_by_id=screensaver_user_admin.screensaver_user_id, 
                comments=activity.comments)
            cpr_parent_log.diffs = {
                'last_screening_activity_date': 
                    [None, activity.date_of_activity ],
            }
            cpr_parent_log.json_field = { 
                'migration': 'CherryPickAssayPlate.screening_date',
                'data': { 
                    'cherry_pick_screening.activity_id': activity.activity_id,
                }
            }
            
            cpr_parent_log.save()
            cpr_parent_logs[cpr.cherry_pick_request_id] = cpr_parent_log

        # NOTE: hack: use the json_field to store the previous activity_id:
        # In the legacy system, more than one screening activity may be recorded
        # for an assay plate (also in the new, but handled differently)
        if cpr_parent_log.json_field is not None:
            prev_activity_id = cpr_parent_log.json_field['data']['cherry_pick_screening.activity_id']
            if prev_activity_id != activity.activity_id:
                # NOTE: for the Cherry Pick Screening, the "performed_by" field is a 
                # screen member; but the ApiLog user will be the admin
                cpr_parent_log = _create_cpr_log(
                    date_of_activity=activity.date_of_activity, 
                    cpr_id=cpr.cherry_pick_request_id, 
                    screen_facility_id=screen_facility_id, 
                    username=screensaver_user_admin.username, 
                    email=screensaver_user_admin.email, 
                    performed_by_id=screensaver_user_admin.screensaver_user_id, 
                    comments=activity.comments)
                cpr_parent_log.diffs = {
                    'last_screening_activity_date': 
                        [None, activity.date_of_activity ],
                }
                cpr_parent_log.json_field = { 
                    'migration': 'CherryPickAssayPlate.screening_date',
                    'data': { 
                        'cherry_pick_screening.activity_id': activity.activity_id,
                    }
                }
                cpr_parent_log.save()
                cpr_parent_logs[cpr.cherry_pick_request_id] = cpr_parent_log

        activity.apilog_uri = cpr_parent_log.log_uri
        activity.save()
        # NOTE: because the cherry_pick_assay_plate_screening_link table (and
        # the cpap->cps many-to-many link) has not been set up, must join manually
        cursor.execute(cpap_to_cps_sql, [activity.activity_id])
        cpap_ids = [x[0] for x in cursor.fetchall()]
        logger.debug('activity: %r, cpap_ids: %r, cpr_parent_log.log_uri: %r', 
            activity.activity_id, cpap_ids,cpr_parent_log.log_uri)
        
        for cpap in CherryPickAssayPlate.objects.all().filter(
            cherry_pick_assay_plate_id__in=cpap_ids):
        # for cpap in cp_screening.cherrypickassayplatescreeninglink_set.all():  
            cpap_log = _child_log_from(cpr_parent_log)
            cpap_log.date_time = create_log_time(cpr_parent_log.date_time)
            cpap_log.ref_resource_name = cpap_resource_name
            cpap_log.key = '/'.join(str(x) for x in [
                cpap.cherry_pick_request_id, 
                cpap.plate_ordinal ])
            cpap_log.uri = '/'.join([base_uri,cpap_log.ref_resource_name,cpap_log.key])
            # Note: previous state will always be "not plated"
            # NOTE: for the Cherry Pick Screening, the "performed_by" field is a 
            # screen member, not the admin
            cpap_log.diffs = { 
                'screening_date': [None, activity.date_of_activity],
                'screened_by_username': [None, screensaver_user.username],
                'screened_by_name': [
                    None, '%s %s' 
                        % (screensaver_user.first_name, screensaver_user.last_name)],
                'screening_comments': [None, activity.comments]
                }
            
            cpap_log.save()
            
            if activity.activity_id == 1049111:
                logger.info('cpap_log saved: %r', cpap_log)
            
            i += 1
            
            if i % 100 == 0:
                logger.info('cpr_parent_log: %r', cpr_parent_log)
                logger.info('created %d cherry pick assay plate screening logs', i)
    logger.info('created %d cherry pick assay plate screening logs', i)
    
class Migration(migrations.Migration):
    
    dependencies = [
        ('db', '0017_cpr'),
    ]
    
    operations = [
        migrations.RunPython(create_well_volume_adjustment_logs),
        migrations.RunPython(create_cherry_pick_plating_logs),
        migrations.RunPython(create_cherry_pick_screening_logs),

        migrations.RunPython(find_extant_plate_logs),
        migrations.RunPython(find_extant_cpr_logs),
        migrations.RunPython(create_library_screening_logs),
        migrations.RunPython(create_plate_plated_and_retired_logs),
        migrations.RunPython(create_plate_generic_logs),
    ]

