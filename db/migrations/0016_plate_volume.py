# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from datetime import datetime, time, date, timedelta
import json
import logging

from django.db import migrations, models
from pytz import timezone
import pytz

from db.support.data_converter import default_converter
from reports.models import ApiLog
from django.db.utils import IntegrityError


logger = logging.getLogger(__name__)

base_uri = '/db/api/v1'
copywell_resource_name = 'copywell'
copywell_uri = '/db/api/v1/' + copywell_resource_name
cpap_resource_name = 'cherrypickassayplate'
cpap_uri = '/db/api/v1/' + cpap_resource_name
plate_resource_name = 'plate'
plate_uri = '/db/api/v1/' + plate_resource_name

# this is a workaround, because some activities have identical date_of_activity


times_seen = set()
# unique offset for the logs in this migration to avoid collisions
plate_vol_time_offset = 1111
def create_log_time(input_date):
    date_time = pytz.timezone('US/Eastern').localize(
        datetime.combine(input_date, datetime.min.time()))
    i = 0
    date_time += timedelta(0,plate_vol_time_offset)
    while date_time in times_seen:
        i += 1
        date_time += timedelta(0,i)
    times_seen.add(date_time)
    return date_time

def _create_generic_log(activity):
    
    log = ApiLog()

    log.comment = activity.comments
    log.date_time = create_log_time(activity.date_of_activity)
    log.username = activity.performed_by.ecommons_id
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
    
    logger.info(str(('create library screening logs')))
    
    LibraryScreening = apps.get_model('db', 'LibraryScreening')
    screen_to_screening_count = {}
    copyplate_to_screening_count = {}
    copyplate_to_volume = {}
    i = 0 
    total_plate_logs = 0
    for screening in LibraryScreening.objects.all().order_by('activity__activity__activity__date_of_activity'):
        lab_activity = screening.activity.activity
        screen = lab_activity.screen
        activity = lab_activity.activity
        # create parent logs:
        logger.debug(str(('for screen', screen.facility_id)))
        
        # screen log for library screening
        screen_log = _create_generic_log(activity)
        screen_log.ref_resource_name = 'screen'
        screen_log.key = screen.facility_id
        screen_log.uri = '/'.join([base_uri,screen_log.ref_resource_name,screen_log.key])
         # TODO: the edit UI will set a "screenings" count
        screen_log.diff_keys = json.dumps(['screening_count'])
        screen_count = screen_to_screening_count.get(screen.facility_id, 0)
        screen_log.diff_keys = json.dumps(['screening_count'])
        screen_log.diffs = json.dumps({'screening_count': [screen_count, screen_count+1] })
        screen_count +=1
        screen_to_screening_count[screen.facility_id] = screen_count
        # TODO: actually create the screening_count var on screen
        try:
            screen_log.save()
        except Exception,e :
            logger.warn(str(('exception on save', screen_log, screen_log.ref_resource_name, e)))
            raise e
        j = 0
        cp_logs = []
        for assay_plate in screening.assayplate_set.all().filter(replicate_ordinal=0):
            plate = assay_plate.plate
            # copyplate log for library screening
            cp_log = _child_log_from(screen_log)
            cp_log.ref_resource_name = 'librarycopyplate' #TODO: "plate" or "copyplate"
            cp_log.key = '/'.join([plate.copy.name, str(plate.plate_number).zfill(5)])
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
            cp_log.diff_keys = json.dumps(['screening_count','remaining_volume'])
            cp_log.diffs = json.dumps({
                'screening_count': [screen_count, screen_count+1],
                'remaining_volume': [str(old_volume),str(new_volume)]
                 })
            screen_count +=1
            copyplate_to_screening_count[cp_log.key] = screen_count
            copyplate_to_volume[cp_log.key] = new_volume
            cp_log.json_field = json.dumps({
                'volume_transferred_per_well_from_library_plates': str(adjustment)
                })
#                 logger.info(str(('plate', cp_log.key, 'assay_plate', assay_plate)))
#                 logger.info(str(('saving', cp_log)))
            cp_logs.append(cp_log)
            j += 1
            # TODO: possibly create logs for the assay plates created (assaywells)
        # todo, use bulk create
        for cp_log in cp_logs:
            cp_log.save()

        i += 1
        total_plate_logs += j
        if total_plate_logs - (total_plate_logs/1000)*1000  == 0:
            logger.info(str(('created screen logs:', i, 'plate logs', total_plate_logs)))
    logger.info(str(('created screen logs:', i, 'plate logs', total_plate_logs)))
        
def create_lcp_logs(apps):
    
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
    
    i = 0
    liquid_transfers = CherryPickLiquidTransfer.objects.all()
    logger.info(str(('create logs for ', len(liquid_transfers), 'liquid_transfers')))
    for cplt in liquid_transfers:
        lab_activity = cplt.activity
        activity = lab_activity.activity
        cpr = cplt.cherrypickassayplate_set.all()[0].cherry_pick_request
        screen_facility_id = lab_activity.screen.facility_id

        # 1. create a parent log on the cherry pick
        extra_information = {}
        cpr_log = _create_cpr_log(cpr, activity)
        
        # Note: "plating_activity" is a pseudo-key: this belies the need for 
        # an "activity" controlleed vocabulary for batch activities.
        cpr_log.diff_keys = json.dumps(['plating_activity'])
        
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
            logger.warn(str(('unknown state', cplt.status,cpr,cplt)))
            
        cpr_log.diffs = json.dumps({
            'plating_activity': [previous_state, cpap_state]
            })
        
        # store other activity information, in case needed
        extra_information['date_of_activity'] = str(activity.date_of_activity)
        extra_information['created_by_id'] = activity.created_by_id
        # n/a for cherrypick plates
        # extra_information['volume_transferred_per_well_from_library_plates'] = \
        #    lab_activity.volume_transferred_per_well_from_library_plates
        extra_information['screen'] = screen_facility_id
        cpr_log.json_field = json.dumps(extra_information)
#             logger.info(str(('create cpr_log', cpr_log)))
        cpr_log.save()
        
        # b. cpap child logs        
        j = 0
        for cpap in cplt.cherrypickassayplate_set.all():
            cpap_log = _child_log_from(cpr_log)
            
            cpap_log.ref_resource_name = "cherrypickassayplate"
            cpap_log.key = '/'.join(str(x) for x in [
                cpap.cherry_pick_request_id, 
                cpap.plate_ordinal, 
                cpap.attempt_ordinal ])
            cpap_log.uri = '/'.join([base_uri,cpap_log.ref_resource_name,cpap_log.key])
            
            cpap_log.diff_keys = json.dumps(['state'])
            cpap_log.diffs = json.dumps([previous_state,cpap_state])
            
            cpap_log.save()
            
            cpap_parent_logs[cpap.cherry_pick_assay_plate_id] = cpap_log
            j += 1
        logger.info(str(('cpr',cpr.cherry_pick_request_id,'cpaps processed',j)))
            # lcp/cpap logs
        
#             if i % 100 == 0:
#                 break
        
    logger.info(str((
        'finished step 1:',len(liquid_transfers), 'cpap parent logs',
        len(cpap_parent_logs))))        
    create_lcp_wva_logs(apps,cpap_parent_logs)

    
def create_lcp_wva_logs(apps,cpap_parent_logs):
    
    logger.info(str(('now create the child logs for all of the '
        'wvas on the cpap on the cplts' )))
    
    # create apilogs for cherry_pick_liquid_transfers
    # 1. well volume adjustments (copy_id, well_id)
    #    a. well volume correction activities
    #    - get date time from correction activity, who, comment
    #    b. lab_cherry_pick -> cherry_pick_assay_plate -> cherry_pick_liquid_transfer -> lab_activity
    #     or cherry_pick_screening -> cpap's -> lcp -> wva
    cols = OrderedDict({
        'well_id': 'well_id' ,
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
    logger.info(str(('query_sql', query_sql)))

    liquid_transfers = CherryPickLiquidTransfer.objects.all()
    logger.info(str(('create logs for ', len(liquid_transfers), 'liquid_transfers')))
    
    count = 0
    parent_log_count = 0
    
    copywells = {}
    prev_logs = set()
    
    for xfer in liquid_transfers:
        cpap = xfer.cherrypickassayplate_set.all()[0] 
        
        if cpap.cherry_pick_assay_plate_id not in cpap_parent_logs:
            raise Exception(str(('could not find a parent log for cpap', cpap,xfer)))
            
        parent_log = cpap_parent_logs[cpap.cherry_pick_assay_plate_id]
        
        parent_log_count += 1
        
        colkeys = cols.keys()
        # FIXME: create a log for each activity/cpap
        _list = db.execute(query_sql, [xfer.activity_id])

        if len(_list) == 0:
            logger.info(str(('no adjustments for ', parent_log,xfer,xfer.status )))
            continue;
    
        i = 0
        prev_volume = None
        for adjustment in _list:
            
            adj = dict(zip(colkeys, adjustment))
            
            log = _child_log_from(parent_log)
            try:
                
                log.ref_resource_name = copywell_resource_name
                log.api_action = 'PATCH'

                log.key = adj['copy_name'] + '/'+ adj['well_id']
                log.uri = copywell_uri + '/' + log.key
                log.json_field = str(adj['volume_adjustment'])

#                     
#                     # FIXME: because this list may not contain all wva's for the well,
#                     # we must do a different, final pass over all wva's for the well 
#                     # to construct the vol change.
#                     log.diff_keys = json.dumps(['volume'])
#                     if log.key in copywells:
#                         prev_volume = copywells[log.key]
#                     else:
#                         prev_volume = round(adj['initial_volume'], 10)
#                         copywells[log.key] = prev_volume
#                     
#                     new_volume = round(prev_volume + float(adj['volume_adjustment']),10)
#                     log.diffs = json.dumps({ 
#                         'volume':[prev_volume, new_volume ]})
#                     copywells[log.key] = new_volume
                
                log.comment = adj['comments'] or ''
                if  adj['legacy_plate_name']:
                    log.comment = log.comment + '. ' + adj['legacy_plate_name']

#                     log.json_field = json.dumps({
#                         'volume_adjustment': round(float(adj['volume_adjustment']),10) })
                if (log.ref_resource_name,log.key,log.date_time) in prev_logs :
                    logger.warn(str(('log key already exists!', log)))
                    log.date_time = log.date_time + timedelta(0,i+parent_log_count) # hack, add second
                
                log.save()
                prev_logs.add((log.ref_resource_name,log.key,log.date_time))
                count += 1
                i += 1
            except Exception, e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
                msg = str((e,activity.activity_id, screen_facility_id ))
                
                logger.warn(str(('on migrate', 
                     msg, exc_type, fname, exc_tb.tb_lineno)))
                raise e
            if count % 10000 == 0:
                logger.info(str(( 'parent_log_count: ', parent_log_count,
                    ' logs created', count )))
#             if parent_log_count > 10: break
            
        logger.info(str(('done, parent_log_count',parent_log_count,'cplt count', i, 'total',count)))
        
    logger.info(str(('total', count, 'parent_log_count',parent_log_count)))   
        
    
def create_well_correction_logs(apps):
    
    logger.info(str(('create well correction logs...')))
    
    cp_pattern = re.compile(r'.*cp(\d+)',re.I)
    
    i = 0
    total_corrections = 0
    query = apps.get_mdodel('db', 'WellVolumeCorrectionActivity').objects.all()\
        .order_by('activity__activity__date_of_activity')
    for wvac in query:
        activity = wvac.activity.activity
        matched = False
        if activity.comments:
            match = cp_pattern.match(activity.comments)
            if match:
                matched = True
                # for a cpr
                cpr_id = match.group(1)
                logger.debug(str(('locate cpr', cpr_id)))
                
                parent_log = None
                try:
                    cpr = apps.get_mdodel('db', 'CherryPickRequest').objects.get(pk=cpr_id)  
                    logger.info(str(('process correction activity for cpr',cpr_id)))
                    parent_log = _create_cpr_log(cpr, activity)
                    
                    parent_log.save()
                except Exception, e:
                    logger.info(str(('could not find cpr', cpr_id)))
                    parent_log = _create_generic_log(activity)
                    
                    # FIXME: need a parent resource?
                    parent_log.ref_resource_name = 'wellvolumecorrectionactivity'
                    parent_log.key = str(activity.activity_id)
                    parent_log.uri = '/'.join([
                        base_uri, parent_log.ref_resource_name, parent_log.key])
                    
                    parent_log.save()
                
                # create wva logs
                j = 0
                for wva in wvac.wellvolumeadjustment_set.all():
                    log = _child_log_from(parent_log)

                    log.ref_resource_name = copywell_resource_name
                    log.key = wva.copy.name + '/'+ wva.well_id
                    log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
                    
                    # temporarily store the adjustment in the json field
                    log.json_field = str(wva.volume)
                    
                    log.save()
                    j += 1
                logger.info(str(('processed', j)))
                total_corrections += j
                
        if not matched:
            # this is a manual adjustment, create a generic parent log
            parent_log = _create_generic_log(activity)
            
            # FIXME: need a parent resource?
            parent_log.ref_resource_name = 'wellvolumecorrectionactivity'
            parent_log.key = str(activity.activity_id)
            parent_log.uri = '/'.join([
                base_uri, parent_log.ref_resource_name, parent_log.key])
            
            parent_log.save()
            logger.info(str(('create generic wvca', parent_log)))
            j = 0
            for wva in wvac.wellvolumeadjustment_set.all():
                log = _child_log_from(parent_log)

                log.ref_resource_name = copywell_resource_name
                log.key = wva.copy.name + '/'+ wva.well_id
                log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
                
                # temporarily store the adjustment in the json field
                log.json_field = str(wva.volume)
                log .save()
                j += 1
            logger.info(str(('processed', j)))
            total_corrections += j
        i += 1       
        
#             if i>10: break
        
    # finally, case c. where the wva has no parent:
    # we know these are all for one copy: 
    # copy_id = 664, name = 'C', library_short_name = 'Human4 Duplexes'
    
    copy = apps.get_mdodel('db', 'Copy').objects.get(pk=664)
    copy1 = apps.get_mdodel('db', 'Copy').objects.get(pk=659)
    
    parent_log = ApiLog()
    parent_log.date_time = datetime.date(2000, 1, 1)
    parent_log.ref_resource_name = 'wellvolumecorrectionactivity'
    parent_log.key = 'unknown'
    parent_log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
    
    parent_log1 = ApiLog()
    parent_log1.date_time = datetime.date(2000, 1, 2)
    parent_log1.ref_resource_name = 'wellvolumecorrectionactivity'
    parent_log1.key = 'unknown'
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
    
    query = apps.get_mdodel('db', 'WellVolumeAdjustment').objects\
        .filter(lab_cherry_pick__isnull=True)\
        .filter(well_volume_correction_activity__isnull=True)\
        .order_by('well__well_id')
    j = 0
    for wva in query:
        if wva.copy_id not in [659,664]:
            raise Exception(str(('manual wva for unknown copy', wva.copy_id,[659,664])))
        
        if wva.copy_id == copy.copy_id:
            log = _child_log_from(parent_log)
        else:
            log = _child_log_from(parent_log1)
            
        log.ref_resource_name = copywell_resource_name
        log.key = wva.copy.name + '/'+ wva.well_id
        log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
        
        # temporarily store the adjustment in the json field
        log.json_field = str(wva.volume,10)
        log .save()
        j += 1
    logger.info(str(('orphaned wvas processed', j)))
    total_corrections += j
    logger.info(str(('done, processed activities', i, 'corrections', total_corrections)))        


def create_copywell_adjustments(apps):
    # Now go back over all of the corrections and construct the diffs
    
    # first get the plate volume = initial volume
    logger.info(str(('create copywell adjustments...')))
    copy_plate_initial_volumes = {}
#         for plate in orm.Plate.objects.all()\
#                 .filter(copy__usage_type='cherry_pick_source_plates'):
    for plate in apps.get_mdodel('db', 'Plate').objects.all():
        key = '%s/%s' % (plate.copy.name,str(plate.plate_number).zfill(5))
        copy_plate_initial_volumes[key] = plate.well_volume
        
    logger.info(str(('plate volume map built, iterate through copywells...')))
    prev_wellcopy_key = None
    prev_plate_key = None
    initial_plate_vol = None
    current_wellcopy_volume = None
    diff_keys = json.dumps(['volume'])
    j = 0
    for log in apps.get_mdodel('reports', 'ApiLog').objects.all()\
            .filter(ref_resource_name='copywell')\
            .order_by('key','date_time'):
        plate_key = log.key.split(':')[0]
        if plate_key != prev_plate_key:
            prev_plate_key = plate_key
            initial_plate_vol = copy_plate_initial_volumes[plate_key]
        
        if log.key != prev_wellcopy_key:
            current_wellcopy_volume = initial_plate_vol
            prev_wellcopy_key = log.key
        # TODO: would be better to make a hash than to store in the json_field
        new_volume = Decimal(log.json_field)
        log.diff_keys = diff_keys
        log.diffs = json.dumps({
            'volume': [str(current_wellcopy_volume),
                       str(new_volume)],
           })
        current_wellcopy_volume = new_volume
        log.json_field = json.dumps({ 'volume_adjustment': log.json_field })
        log.save()
        j += 1
        if j % 10000 == 0:
            logger.info(str(('processed',j)))
    logger.info(str(('volume apilogs adjusted', j)))

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
    # self.create_lcp_logs(orm) 
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
    ]
