# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import DataMigration
from django.utils import timezone, tzinfo
from django.db import models
from reports.models import ApiLog
import json
from django.utils.timezone import make_aware
import datetime
import logging
from collections import OrderedDict
from db.models import CherryPickScreening, CherryPickLiquidTransfer,\
    LibraryScreening
from reports.models import ApiLog
import sys
import os
import re
logger = logging.getLogger(__name__)


base_uri = '/db/api/v1'
copywell_resource_name = 'copywell'
copywell_uri = '/db/api/v1/' + copywell_resource_name
cpap_resource_name = 'cherrypickassayplate'
cpap_uri = '/db/api/v1/' + cpap_resource_name
plate_resource_name = 'plate'
plate_uri = '/db/api/v1/' + plate_resource_name

class Migration(DataMigration):

    def _create_generic_log(self, activity):
        
        log = ApiLog()

        log.comment = activity.comments
        log.date_time = datetime.datetime.combine(
            activity.date_of_activity, datetime.time())
#             log.date_time = make_aware(
#                 activity.date_of_activity, timezone.get_default_timezone())
        if log.date_time not in self.times_seen:
            self.times_seen.add(log.date_time)
        else:
            i = 0
            while log.date_time in self.times_seen:
                i += 1
                log.date_time += datetime.timedelta(0,i)
            self.times_seen.add(log.date_time)
        log.date_time = make_aware(log.date_time, timezone.get_default_timezone())
        
        
        log.username = activity.performed_by.ecommons_id
        if not log.username:
            log.username = '%s: %s' % (
                activity.performed_by.screensaver_user_id,
                activity.performed_by.email)
        
        log.user_id = activity.performed_by.screensaver_user_id
        
        return log

    def _create_cpr_log(self, cpr, activity):
        
        cpr_log = self._create_generic_log(activity)
        cpr_log.ref_resource_name = 'cherrypickrequest'
        cpr_log.key = str(cpr.cherry_pick_request_id)
        cpr_log.uri = '/'.join([base_uri,cpr_log.ref_resource_name,cpr_log.key])
        cpr_log.api_action = 'PATCH'
    
        return cpr_log
    
    def _child_log_from(self,parent_log):
        child_log = ApiLog()
        child_log.parent_log = parent_log
        child_log.username = parent_log.username
        child_log.user_id = parent_log.user_id
        child_log.date_time = parent_log.date_time
        child_log.api_action = parent_log.api_action
        child_log.comment = parent_log.comment
        return child_log

    def forwards(self,orm):
        self.times_seen = set() # hack, because some activities have identical date_of_activity

        # 1. assay_plate / library screenings xfers:
        #    - assay_plate -> library_screening ->> lab_activity
        # library_screening -> assay_plate -> plate
        #                                    -> screen
        self.create_library_screening_logs(orm)

#         # 2. create ApiLogs for WVA's: temporarily store volume adjusment on the log
# 
#         # a. do cplt's
#         self.create_lcp_logs(orm) 
#                
#         # b. do wvca's
#         # c. do orphaned wvas (no cplt or wvca attached)
#         self.create_well_correction_logs(orm)
#         
#         # d. do a final pass, over the logs sorted by date, in order to 
#         #    construct the previous_volume, new volume in the log
#         self.create_copywell_adjustments(orm)
        
        
    def create_library_screening_logs(self,orm):
        
        logger.info(str(('create library screening logs')))
    
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
            screen_log = self._create_generic_log(activity)
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
                cp_log = self._child_log_from(screen_log)
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
                     # not sure what these "library screenings" with no volume are still:
                     # -- some are external library plates, presumably
                     # -- some have not AP's and are "z prime" logs?
                     # -- some are legacy records from before ss 2.0
                    adjustment = 0
                adjustment = round(float(adjustment),10)
                new_volume = round(old_volume-adjustment,10)
                cp_log.diff_keys = json.dumps(['screening_count','remaining_volume'])
                cp_log.diffs = json.dumps({
                    'screening_count': [screen_count, screen_count+1],
                    'remaining_volume': [old_volume,new_volume]
                     })
                screen_count +=1
                copyplate_to_screening_count[cp_log.key] = screen_count
                copyplate_to_volume[cp_log.key] = new_volume
                cp_log.json_field = json.dumps({
                    'volume_transferred_per_well_from_library_plates': adjustment
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
            
    def create_lcp_logs(self, orm):
        
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
            cpr_log = self._create_cpr_log(cpr, activity)
            
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
                cpap_log = self._child_log_from(cpr_log)
                
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
        self.create_lcp_wva_logs(orm,cpap_parent_logs)
    
        
    def create_lcp_wva_logs(self, orm,cpap_parent_logs):
        
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
                
                log = self._child_log_from(parent_log)
                try:
                    
                    log.ref_resource_name = copywell_resource_name
                    log.api_action = 'PATCH'
    
                    log.key = adj['copy_name'] + '/'+ adj['well_id']
                    log.uri = copywell_uri + '/' + log.key
                    log.json_field = str(round(float(adj['volume_adjustment']),10))

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
                        log.date_time = log.date_time + datetime.timedelta(0,i+parent_log_count) # hack, add second
                    
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
            
        
    def create_well_correction_logs(self, orm):
        
        logger.info(str(('create well correction logs...')))
        
        cp_pattern = re.compile(r'.*cp(\d+)',re.I)
        
        i = 0
        total_corrections = 0
        query = orm.WellVolumeCorrectionActivity.objects.all()\
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
                        cpr = orm.CherryPickRequest.objects.get(pk=cpr_id)  
                        logger.info(str(('process correction activity for cpr',cpr_id)))
                        parent_log = self._create_cpr_log(cpr, activity)
                        
                        parent_log.save()
                    except Exception, e:
                        logger.info(str(('could not find cpr', cpr_id)))
                        parent_log = self._create_generic_log(activity)
                        
                        # FIXME: need a parent resource?
                        parent_log.ref_resource_name = 'wellvolumecorrectionactivity'
                        parent_log.key = str(activity.activity_id)
                        parent_log.uri = '/'.join([
                            base_uri, parent_log.ref_resource_name, parent_log.key])
                        
                        parent_log.save()
                    
                    # create wva logs
                    j = 0
                    for wva in wvac.wellvolumeadjustment_set.all():
                        log = self._child_log_from(parent_log)

                        log.ref_resource_name = copywell_resource_name
                        log.key = wva.copy.name + '/'+ wva.well_id
                        log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
                        
                        # temporarily store the adjustment in the json field
                        log.json_field = str(round(wva.volume,10))
                        
                        log.save()
                        j += 1
                    logger.info(str(('processed', j)))
                    total_corrections += j
                    
            if not matched:
                # this is a manual adjustment, create a generic parent log
                parent_log = self._create_generic_log(activity)
                
                # FIXME: need a parent resource?
                parent_log.ref_resource_name = 'wellvolumecorrectionactivity'
                parent_log.key = str(activity.activity_id)
                parent_log.uri = '/'.join([
                    base_uri, parent_log.ref_resource_name, parent_log.key])
                
                parent_log.save()
                logger.info(str(('create generic wvca', parent_log)))
                j = 0
                for wva in wvac.wellvolumeadjustment_set.all():
                    log = self._child_log_from(parent_log)

                    log.ref_resource_name = copywell_resource_name
                    log.key = wva.copy.name + '/'+ wva.well_id
                    log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
                    
                    # temporarily store the adjustment in the json field
                    log.json_field = str(round(wva.volume,10))
                    log .save()
                    j += 1
                logger.info(str(('processed', j)))
                total_corrections += j
            i += 1       
            
#             if i>10: break
            
        # finally, case c. where the wva has no parent:
        # we know these are all for one copy: 
        # copy_id = 664, name = 'C', library_short_name = 'Human4 Duplexes'
        
        copy = orm.Copy.objects.get(pk=664)
        copy1 = orm.Copy.objects.get(pk=659)
        
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
        
        query = orm.WellVolumeAdjustment.objects\
            .filter(lab_cherry_pick__isnull=True)\
            .filter(well_volume_correction_activity__isnull=True)\
            .order_by('well__well_id')
        j = 0
        for wva in query:
            if wva.copy_id not in [659,664]:
                raise Exception(str(('manual wva for unknown copy', wva.copy_id,[659,664])))
            
            if wva.copy_id == copy.copy_id:
                log = self._child_log_from(parent_log)
            else:
                log = self._child_log_from(parent_log1)
                
            log.ref_resource_name = copywell_resource_name
            log.key = wva.copy.name + '/'+ wva.well_id
            log.uri = '/'.join([base_uri, log.ref_resource_name, log.key])
            
            # temporarily store the adjustment in the json field
            log.json_field = str(round(wva.volume,10))
            log .save()
            j += 1
        logger.info(str(('orphaned wvas processed', j)))
        total_corrections += j
        logger.info(str(('done, processed activities', i, 'corrections', total_corrections)))        

    
    def create_copywell_adjustments(self, orm):
        # Now go back over all of the corrections and construct the diffs
        
        # first get the plate volume = initial volume
        logger.info(str(('create copywell adjustments...')))
        copy_plate_initial_volumes = {}
#         for plate in orm.Plate.objects.all()\
#                 .filter(copy__usage_type='cherry_pick_source_plates'):
        for plate in orm.Plate.objects.all():
            key = '%s/%s' % (plate.copy.name,str(plate.plate_number).zfill(5))
            copy_plate_initial_volumes[key] = plate.well_volume
            
        logger.info(str(('plate volume map built, iterate through copywells...')))
        prev_wellcopy_key = None
        prev_plate_key = None
        initial_plate_vol = None
        current_wellcopy_volume = None
        diff_keys = json.dumps(['volume'])
        j = 0
        for log in orm['reports.ApiLog'].objects.all()\
                .filter(ref_resource_name='copywell')\
                .order_by('key','date_time'):
            plate_key = log.key.split(':')[0]
            if plate_key != prev_plate_key:
                prev_plate_key = plate_key
                initial_plate_vol = copy_plate_initial_volumes[plate_key]
            
            if log.key != prev_wellcopy_key:
                current_wellcopy_volume = round(float(initial_plate_vol),10)
                prev_wellcopy_key = log.key
            new_volume = round(current_wellcopy_volume+float(log.json_field),10)
            log.diff_keys = diff_keys
            log.diffs = json.dumps({
                'volume': [str(current_wellcopy_volume),
                           str(new_volume)],
               })
            current_wellcopy_volume = new_volume
            log.json_field = json.dumps({ 'volume_adjustment': log.json_field })
            log.save()
#             logger.info(str(('log', log.key, current_wellcopy_volume, 
#                 float(log.json_field), log.diffs)))
            j += 1
            if j % 10000 == 0:
                logger.info(str(('processed',j)))
        logger.info(str(('volume apilogs adjusted', j)))
        
#     def create_logs(self, copy):
#         
#         cols = OrderedDict({
#             'well_id': 'well_id' ,
#             'copy_name':'c.name',
#             'plate_number':'p.plate_number',
#             'volume_adjustment': 'wva.volume',
#             'initial_volume': 'p.well_volume',
#             'comments': 'a.comments',
#             'ecommons_id': 'u.ecommons_id',
#             'email': 'u.email',
#             'date_created': 'a.date_created',
#             'date_of_activity': 'a.date_of_activity',
#             'performed_by_id': 'performed_by_id',
#             'login_id': 'u.login_id'
#             })
#         _cols = ', '.join([ '%s as %s' % (value,key) for key, value in cols.items() ])
# 
#         # TODO: re-organize this by 
#         # 1. cpap, as the parent log
#         # 1.a or cpap->cherry_pick_screening
#         #    cherry pick assay plate states: 
#         #    - not plated: 
#         #    - plated: _cherryPickLiquidTransfer != null && _cherryPickLiquidTransfer.isSuccessful();
#         #    - failed: cherryPickLiquidTransfer != null && _cherryPickLiquidTransfer.isFailed();
#         #    - cancelled: _cherryPickLiquidTransfer != null && _cherryPickLiquidTransfer.isCancelled();
#         #    - plated & screened: !_cherryPickScreenings.isEmpty();
#         # ** wva occurs even if the cplt is failed, see
#         # LibrariesDao.findRemainingVolumesInWellCopies
#         
#         #    lab cherry pick states:
#         #    - unfulfilled 
#         #    - allocated: wva's > 0
#         #    - mapped: CherryPickAssayPlate != null
#         #    - mapped+unallocated
#         #    - mapped+allocated
#         #    - failed
#         #    - canceled: assayPlate.isCancelled
#         #    - plated: wva's>0, assayPlate.isPlated: cherryPickLiquidTransfer != null
#         # 2. cherry_pick_request as the parent, parent
#         adjustment_queries = {
#             'corrections': '\n'.join([
#                 'select ',
#                 _cols ,
#                 'from well_volume_adjustment wva ' ,
#                 'join copy c using(copy_id) ',
#                 'join well w using(well_id) ',
#                 'join plate p on(wva.copy_id=p.copy_id and w.plate_number=p.plate_number) ',
#                 'join activity a on(wva.well_volume_correction_activity_id = a.activity_id) ', 
#                 'join screensaver_user u on(u.screensaver_user_id=a.performed_by_id)' 
#                 ' where c.copy_id=%s ',
#                 'order by c.name, well_id,wva.well_volume_adjustment_id ',
#                 ]),
#             'cherry_picks': '\n'.join([
#                 'select ',
#                 _cols ,
#                 ', cpap.legacy_plate_name as legacy_plate_name',
#                 'from well_volume_adjustment wva ' ,
#                 'join copy c using(copy_id) ',
#                 'join well w using(well_id) ',
#                 'join plate p on(wva.copy_id=p.copy_id and w.plate_number=p.plate_number) ',
#                 'join lab_cherry_pick lcp on (wva.lab_cherry_pick_id=lcp.lab_cherry_pick_id) ',
#                 'join cherry_pick_assay_plate cpap on (cpap.cherry_pick_assay_plate_id = lcp.cherry_pick_assay_plate_id)',
#                 'join activity a on(a.activity_id = cpap.cherry_pick_liquid_transfer_id)',
#                 'join screensaver_user u on(u.screensaver_user_id=a.performed_by_id)' 
#                 ' where c.copy_id=%s ',
#                 'AND p.well_volume is not null ',
#                 'order by c.name, well_id,wva.well_volume_adjustment_id ',
#                 ]),
#             
#         }
#             
#         copywell_resource_name = 'copywell'
#         copywell_uri = '/db/api/v1/' + copywell_resource_name
#         plate_resource_name = 'plate'
#         plate_uri = '/db/api/v1/' + plate_resource_name
#         
#         # 1. well volume adjustments (copy_id, well_id)
#         count = 0
#         prev_volume = None
#         prev_well_id = None
#         
#         prev_times = set() # this one is a hack,because some of the activity times are not uniq
#         current_plate_id = None
#         current_volume = None
#         
#         prev_logs = set()
#         
#         for key,query_sql in adjustment_queries.items():
# #             logger.info(str(('sql', key, query_sql)))
# #             logger.info(str(('key', key,copy.copy_id, copy.library.short_name,copy.name,query_sql )))
#             colkeys = cols.keys()
#             if key == 'cherry_picks':
#                 colkeys.append('legacy_plate_name')
#             _list = db.execute(query_sql, [copy.copy_id])
# 
#             if len(_list) == 0:
#                 logger.info(str(('no adjustments for ', key, copy.library.short_name, copy.name )))
#             
#             i = 0;
#             prev_log = None
#             for adjustment in _list:
#                 
#                 adj = dict(zip(colkeys, adjustment))
# #                 logger.info(str(('adj', adj)))
#                 
#                 if adj['well_id'] != prev_well_id:
#                     prev_well_id = adj['well_id']
#                     plate_id = (adj['copy_name'],adj['plate_number'])
#                     if ( not current_plate_id or
#                             plate_id != current_plate_id):
#                         current_plate_id = plate_id
#                         prev_volume = round(adj['initial_volume'], 10)
#                     prev_times = set()
#                 log = ApiLog()
# 
#                 try:
#                     log.username = adj['ecommons_id']
#                     if not log.username:
#                         # TODO: construct a username
#                         log.username = adj['login_id']
#                     
#                     if not log.username:
#                         log.username = adj['performed_by_id']
#     #                     logger.debug(str(('no username found: ', copy.copy_id, adj)))
#                     
#                     # log.user_id = getattr(activity.performed_by.user, 'id', log.username)
#                     log.user_id = 1    
#                     if 'performed_by_id' in adj:
#                         log.user_id = adj['performed_by_id']
#                         
#                     log.ref_resource_name = 'copywell'
#                     log.api_action = 'PATCH'
#     
#                     log.key = adj['copy_name'] + '/'+ adj['well_id']
#                     log.uri = copywell_uri + '/' + log.key
#                     log.diff_keys = json.dumps(['volume'])
#                     
#                     new_volume = round(prev_volume + float(adj['volume_adjustment']),10)
#                     log.diffs = json.dumps([
#                         prev_volume, new_volume ])
#                     prev_volume = new_volume
#                     
#                     log.comment = adj['comments']
#                     if key == 'cherry_picks' and adj['legacy_plate_name']:
#                         if log.comment:
#                             log.comment = log.comment + '. ' + adj['legacy_plate_name']
#                         else:
#                             log.comment = adj['legacy_plate_name']
#     #                 logger.info(str(('created log', log)))
#     
#                     log.date_time = make_aware(
#                         adj['date_created'], timezone.get_default_timezone())
#                     
#                     if (log.ref_resource_name,log.key,log.date_time) in prev_logs :
#                         log.date_time = log.date_time + datetime.timedelta(0,count) # hack, add second
#     
#                     # TODO: create a parent log for the lab_cherry_pick, and/or for
#                     # the cherry_pick_assay_plate creation
#                     # log.parent_id = cherry_pick_assay_plate_log_id
#                 
#                     # finally, dump everything known into the json field
#                     log.json_field = json.dumps({
#                         k:str(v) for k,v in adj.items() if k not in [
#                             'comments','date_created','plate_number','copy_name','well_id']})
# 
#                     log.save()
#                     prev_logs.add((log.ref_resource_name,log.key,log.date_time))
#                     count += 1
#                     i += 1
#                 except Exception, e:
#                     msg = str(('exception on save: ', log,', sql', query_sql, copy, adj, e))
#                     logger.info(msg)
#                     raise e
#                 if count % 10000 == 0:
#                     logger.info(str(( 'Count: ', count, ' logs created', datetime.datetime.now() )))
#             logger.info(str((key, ',', i)))
# 
#         # TODO: reorganize this using the library screening, screen as the parent logs
#         logger.info(str(('2. assay_plate / library screenings xfers for ', copy.name)))
#         cols = OrderedDict({
#             'copy_name':'c.name',
#             'plate_number':'p.plate_number',
#             'volume_adjustment': '-la.volume_transferred_per_well_from_library_plates',
#             'initial_volume': 'p.well_volume',
#             'comments': 'a.comments',
#             'ecommons_id': 'u.ecommons_id',
#             'email': 'u.email',
#             'date_created': 'a.date_created',
#             'date_of_activity': 'a.date_of_activity',
#             'performed_by_id': 'performed_by_id',
#             'login_id': 'u.login_id'
#             })
#         colkeys = cols.keys()
#         _cols = ', '.join([ '%s as %s' % (value,key) 
#             for key, value in cols.items() ])
#         query_sql = '\n'.join([
#             'select',
#             _cols ,
#             'from ',
#             'plate p ',
#             'join copy c on(p.copy_id=c.copy_id) ',
#             'join assay_plate ap using(plate_id) ',
#             'join screening ls on(ls.activity_id = ap.library_screening_id) ',
#             'join lab_activity la using(activity_id) ',
#             'join activity a using(activity_id) ',
#             'join screensaver_user u on(u.screensaver_user_id=a.performed_by_id)' 
#             'where ap.replicate_ordinal = 0 ',
#             'AND la.volume_transferred_per_well_from_library_plates is not null '
#             'AND p.well_volume is not null ',
#             'AND c.copy_id=%s ',
#             'order by p.plate_number, c.name,a.activity_id '
#             ])
#         _list = db.execute(query_sql, [copy.copy_id])
# 
#         if len(_list) == 0:
#             logger.info(str(('no library screening adjustments for ', 
#                 key, copy.library.short_name, copy.name )))
#         i = 0;
#         prev_log = None
#         prev_resource_id = None
#         prev_volume = None
#         for adjustment in _list:
#             
#             adj = dict(zip(colkeys, adjustment))
# #             logger.info(str((adj)))
#             resource_id = '%s/%s' % (adj['plate_number'],adj['copy_name'])
#             
#             if resource_id != prev_resource_id:
#                 prev_resource_id = resource_id
#                 prev_times = set()
#                 prev_volume = adj['initial_volume']
#     
#             log = ApiLog()
#             
#             try:
#     
#                 log.username = adj['ecommons_id']
#                 if not log.username:
#                     log.username = adj['login_id']
#                 if not log.username:
#                     log.username = adj['performed_by_id']
#                 
#                 log.user_id = 1    
#                 if 'performed_by_id' in adj:
#                     log.user_id = adj['performed_by_id']
#                     
#                 log.ref_resource_name = plate_resource_name
#                 log.api_action = 'PATCH'
#     
#                 log.key = resource_id
#                 log.uri = plate_uri  + '/' + log.key
#                 log.diff_keys = json.dumps(['remaining_volume'])
#                 
#                 new_volume = round(prev_volume + float(adj['volume_adjustment']),10)
#                 log.diffs = json.dumps([
#                     prev_volume, new_volume ])
#                 prev_volume = new_volume
#                 
#                 log.comment = adj['comments']
#     
#                 log.date_time = make_aware(
#                     adj['date_created'], timezone.get_default_timezone())
#                 
#                 if (log.ref_resource_name,log.key,log.date_time) in prev_logs :
#                     log.date_time = log.date_time + datetime.timedelta(0,count) # hack, add second
#     
#                 # TODO: create a parent log for the library screening, and/or for
#                 # the assay_plate creation
#                 # log.parent_id = assay_plate_log_id
#                 
#                 # finally, dump everything known into the json field
#                 log.json_field = json.dumps({
#                     k:str(v) for k,v in adj.items() if k not in ['comments','date_created','plate_number','copy_name']})
# 
#                 log.save()
#                 prev_logs.add((log.ref_resource_name,log.key,log.date_time))
#                 count += 1
#                 i = i+1
#             except Exception, e:
#                 msg = str(('exception on save,', log, ', sql', query_sql, copy, adj, e))
#                 logger.info(msg)
#                 raise e
#         
#         logger.info(str(('library_screening', i)))
#         logger.info(str(('Completed', copy.library.short_name, copy.name, count, ' logs')))
        

    def backwards(self, orm):
        # Deleting field 'Plate.remaining_volume'
        db.delete_column(u'plate', 'remaining_volume')


        # Changing field 'Plate.well_volume'
        db.alter_column(u'plate', 'well_volume', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9))

    models = {
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'db.abasetestset': {
            'Meta': {'object_name': 'AbaseTestset', 'db_table': "u'abase_testset'"},
            'abase_testset_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'testset_date': ('django.db.models.fields.DateField', [], {}),
            'testset_name': ('django.db.models.fields.TextField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.activity': {
            'Meta': {'object_name': 'Activity', 'db_table': "u'activity'"},
            'activity_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'activities_created'", 'null': 'True', 'to': u"orm['db.ScreensaverUser']"}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_of_activity': ('django.db.models.fields.DateField', [], {}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'performed_by': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'activities_performed'", 'to': u"orm['db.ScreensaverUser']"}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.activityupdateactivity': {
            'Meta': {'object_name': 'ActivityUpdateActivity', 'db_table': "u'activity_update_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.administrativeactivity': {
            'Meta': {'object_name': 'AdministrativeActivity', 'db_table': "u'administrative_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']", 'primary_key': 'True'}),
            'administrative_activity_type': ('django.db.models.fields.TextField', [], {})
        },
        u'db.administratoruser': {
            'Meta': {'object_name': 'AdministratorUser', 'db_table': "u'administrator_user'"},
            'screensaver_user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.ScreensaverUser']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'db.annotationtype': {
            'Meta': {'object_name': 'AnnotationType', 'db_table': "u'annotation_type'"},
            'annotation_type_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_numeric': ('django.db.models.fields.BooleanField', [], {}),
            'name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'study': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.annotationvalue': {
            'Meta': {'object_name': 'AnnotationValue', 'db_table': "u'annotation_value'"},
            'annotation_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AnnotationType']", 'null': 'True', 'blank': 'True'}),
            'annotation_value_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'numeric_value': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']", 'null': 'True', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.assayplate': {
            'Meta': {'object_name': 'AssayPlate', 'db_table': "u'assay_plate'"},
            'assay_plate_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'library_screening': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LibraryScreening']", 'null': 'True', 'blank': 'True'}),
            'plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Plate']", 'null': 'True', 'blank': 'True'}),
            'plate_number': ('django.db.models.fields.IntegerField', [], {}),
            'replicate_ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'screen_result_data_loading': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.assaywell': {
            'Meta': {'object_name': 'AssayWell', 'db_table': "u'assay_well'"},
            'assay_well_control_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_well_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'confirmed_positive_value': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_positive': ('django.db.models.fields.BooleanField', [], {}),
            'screen_result': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenResult']"}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"})
        },
        u'db.attachedfile': {
            'Meta': {'object_name': 'AttachedFile', 'db_table': "u'attached_file'"},
            'attached_file_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'attached_file_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AttachedFileType']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'file_contents': ('django.db.models.fields.TextField', [], {}),
            'file_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'filename': ('django.db.models.fields.TextField', [], {}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']", 'null': 'True', 'blank': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']", 'null': 'True', 'blank': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']", 'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.attachedfiletype': {
            'Meta': {'object_name': 'AttachedFileType', 'db_table': "u'attached_file_type'"},
            'attached_file_type_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'for_entity_type': ('django.db.models.fields.CharField', [], {'max_length': '31'}),
            'value': ('django.db.models.fields.TextField', [], {})
        },
        u'db.attachedfileupdateactivity': {
            'Meta': {'object_name': 'AttachedFileUpdateActivity', 'db_table': "u'attached_file_update_activity'"},
            'attached_file': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AttachedFile']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.cachedquery': {
            'Meta': {'object_name': 'CachedQuery', 'db_table': "u'cached_query'"},
            'count': ('django.db.models.fields.IntegerField', [], {'null': 'True'}),
            'datetime': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'key': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'params': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'sql': ('django.db.models.fields.TextField', [], {}),
            'uri': ('django.db.models.fields.TextField', [], {}),
            'username': ('django.db.models.fields.CharField', [], {'max_length': '128'})
        },
        u'db.cellline': {
            'Meta': {'object_name': 'CellLine', 'db_table': "u'cell_line'"},
            'cell_line_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.checklistitem': {
            'Meta': {'object_name': 'ChecklistItem', 'db_table': "u'checklist_item'"},
            'checklist_item_group': ('django.db.models.fields.TextField', [], {}),
            'checklist_item_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'is_expirable': ('django.db.models.fields.BooleanField', [], {}),
            'item_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'order_statistic': ('django.db.models.fields.IntegerField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.checklistitemevent': {
            'Meta': {'object_name': 'ChecklistItemEvent', 'db_table': "u'checklist_item_event'"},
            'checklist_item_event_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'checklist_item_id': ('django.db.models.fields.IntegerField', [], {}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_performed': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'is_expiration': ('django.db.models.fields.BooleanField', [], {}),
            'is_not_applicable': ('django.db.models.fields.BooleanField', [], {}),
            'screening_room_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"})
        },
        u'db.checklistitemeventupdateactivity': {
            'Meta': {'object_name': 'ChecklistItemEventUpdateActivity', 'db_table': "u'checklist_item_event_update_activity'"},
            'checklist_item_event': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ChecklistItemEvent']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.cherrypickassayplate': {
            'Meta': {'unique_together': "((u'cherry_pick_request', u'plate_ordinal', u'attempt_ordinal'),)", 'object_name': 'CherryPickAssayPlate', 'db_table': "u'cherry_pick_assay_plate'"},
            'assay_plate_type': ('django.db.models.fields.TextField', [], {}),
            'attempt_ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'cherry_pick_assay_plate_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'cherry_pick_assay_plate_type': ('django.db.models.fields.CharField', [], {'max_length': '31'}),
            'cherry_pick_liquid_transfer': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickLiquidTransfer']", 'null': 'True', 'blank': 'True'}),
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            'legacy_plate_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'plate_ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'status': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.cherrypickassayplatescreeninglink': {
            'Meta': {'object_name': 'CherryPickAssayPlateScreeningLink', 'db_table': "u'cherry_pick_assay_plate_screening_link'"},
            'cherry_pick_assay_plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickAssayPlate']"}),
            'cherry_pick_screening': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickScreening']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.cherrypickliquidtransfer': {
            'Meta': {'object_name': 'CherryPickLiquidTransfer', 'db_table': "u'cherry_pick_liquid_transfer'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabActivity']", 'primary_key': 'True'}),
            'status': ('django.db.models.fields.TextField', [], {})
        },
        u'db.cherrypickrequest': {
            'Meta': {'object_name': 'CherryPickRequest', 'db_table': "u'cherry_pick_request'"},
            'assay_plate_type': ('django.db.models.fields.TextField', [], {}),
            'assay_protocol_comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cherry_pick_assay_protocols_followed': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cherry_pick_followup_results_status': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cherry_pick_request_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_requested': ('django.db.models.fields.DateField', [], {}),
            'date_volume_approved': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'is_randomized_assay_plate_layout': ('django.db.models.fields.BooleanField', [], {}),
            'keep_source_plate_cherry_picks_together': ('django.db.models.fields.BooleanField', [], {}),
            'legacy_cherry_pick_request_number': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'max_skipped_wells_per_plate': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'number_unfulfilled_lab_cherry_picks': ('django.db.models.fields.IntegerField', [], {}),
            'requested_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'transfer_volume_per_well_approved': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'transfer_volume_per_well_requested': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'volume_approved_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministratorUser']", 'null': 'True', 'blank': 'True'})
        },
        u'db.cherrypickrequestemptywell': {
            'Meta': {'object_name': 'CherryPickRequestEmptyWell', 'db_table': "u'cherry_pick_request_empty_well'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'well_name': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'})
        },
        u'db.cherrypickrequestupdateactivity': {
            'Meta': {'object_name': 'CherryPickRequestUpdateActivity', 'db_table': "u'cherry_pick_request_update_activity'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'unique': 'True'})
        },
        u'db.cherrypickscreening': {
            'Meta': {'object_name': 'CherryPickScreening', 'db_table': "u'cherry_pick_screening'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screening']", 'primary_key': 'True'}),
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"})
        },
        u'db.collaboratorlink': {
            'Meta': {'object_name': 'CollaboratorLink', 'db_table': "u'collaborator_link'"},
            'collaborator': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.copy': {
            'Meta': {'object_name': 'Copy', 'db_table': "u'copy'"},
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'copy_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_plated': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'max_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'max_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'min_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'min_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'plate_locations_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'plates_available': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'primary_plate_location_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'primary_plate_status': ('django.db.models.fields.TextField', [], {}),
            'primary_well_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'primary_well_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'usage_type': ('django.db.models.fields.TextField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well_concentration_dilution_factor': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '8', 'decimal_places': '2', 'blank': 'True'})
        },
        u'db.copyupdateactivity': {
            'Meta': {'object_name': 'CopyUpdateActivity', 'db_table': "u'copy_update_activity'"},
            'copy_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'})
        },
        u'db.copywell': {
            'Meta': {'object_name': 'CopyWell', 'db_table': "u'copy_well'"},
            'adjustments': ('django.db.models.fields.IntegerField', [], {}),
            'copy': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Copy']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'initial_volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
#             'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Plate']"}),
            'plate_number': ('django.db.models.fields.IntegerField', [], {}),
            'volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"}),
        },
        u'db.datacolumn': {
            'Meta': {'object_name': 'DataColumn', 'db_table': "u'data_column'"},
            'assay_phenotype': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_readout_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'channel': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'data_column_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'data_type': ('django.db.models.fields.TextField', [], {}),
            'decimal_places': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'how_derived': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_derived': ('django.db.models.fields.BooleanField', [], {}),
            'is_follow_up_data': ('django.db.models.fields.BooleanField', [], {}),
            'medium_positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'replicate_ordinal': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'screen_result': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenResult']"}),
            'strong_positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'time_point': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'time_point_ordinal': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'weak_positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'zdepth_ordinal': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'})
        },
        u'db.datacolumnderivedfromlink': {
            'Meta': {'object_name': 'DataColumnDerivedFromLink', 'db_table': "u'data_column_derived_from_link'"},
            'derived_data_column': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.DataColumn']"}),
            'derived_from_data_column': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'derived_from'", 'to': u"orm['db.DataColumn']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.equipmentused': {
            'Meta': {'object_name': 'EquipmentUsed', 'db_table': "u'equipment_used'"},
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'equipment': ('django.db.models.fields.TextField', [], {}),
            'equipment_used_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'lab_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabActivity']"}),
            'protocol': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.fundingsupport': {
            'Meta': {'object_name': 'FundingSupport', 'db_table': "u'funding_support'"},
            'funding_support_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'unique': 'True', 'blank': 'True'})
        },
        u'db.gene': {
            'Meta': {'object_name': 'Gene', 'db_table': "u'gene'"},
            'entrezgene_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'gene_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'gene_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'species_name': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.genegenbankaccessionnumber': {
            'Meta': {'unique_together': "((u'gene', u'genbank_accession_number'),)", 'object_name': 'GeneGenbankAccessionNumber', 'db_table': "u'gene_genbank_accession_number'"},
            'genbank_accession_number': ('django.db.models.fields.TextField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.genesymbol': {
            'Meta': {'unique_together': "((u'gene', u'ordinal'),)", 'object_name': 'GeneSymbol', 'db_table': "u'gene_symbol'"},
            'entrezgene_symbol': ('django.db.models.fields.TextField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.labactivity': {
            'Meta': {'object_name': 'LabActivity', 'db_table': "u'lab_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']", 'primary_key': 'True'}),
            'molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'volume_transferred_per_well_from_library_plates': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'})
        },
        u'db.labaffiliation': {
            'Meta': {'object_name': 'LabAffiliation', 'db_table': "u'lab_affiliation'"},
            'affiliation_category': ('django.db.models.fields.TextField', [], {}),
            'affiliation_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'lab_affiliation_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.labcherrypick': {
            'Meta': {'object_name': 'LabCherryPick', 'db_table': "u'lab_cherry_pick'"},
            'assay_plate_column': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'assay_plate_row': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'cherry_pick_assay_plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickAssayPlate']", 'null': 'True', 'blank': 'True'}),
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            'lab_cherry_pick_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'screener_cherry_pick': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenerCherryPick']"}),
            'source_well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.labhead': {
            'Meta': {'object_name': 'LabHead', 'db_table': "u'lab_head'"},
            'lab_affiliation': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabAffiliation']", 'null': 'True', 'blank': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.ScreeningRoomUser']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'db.legacysmallmoleculecasnumber': {
            'Meta': {'object_name': 'LegacySmallMoleculeCasNumber', 'db_table': "u'_legacy_small_molecule_cas_number'"},
            'cas_number': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'smiles': ('django.db.models.fields.CharField', [], {'max_length': '2047'})
        },
        u'db.library': {
            'Meta': {'object_name': 'Library', 'db_table': "u'library'"},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_received': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_screenable': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'end_plate': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'experimental_well_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'is_pool': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'latest_released_contents_version_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'library_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'library_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'library_type': ('django.db.models.fields.TextField', [], {}),
            'loaded_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'libraries_loaded'", 'null': 'True', 'to': u"orm['db.ScreensaverUser']"}),
            'owner_screener': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']", 'null': 'True', 'blank': 'True'}),
            'plate_size': ('django.db.models.fields.TextField', [], {}),
            'provider': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'screen_type': ('django.db.models.fields.TextField', [], {}),
            'screening_status': ('django.db.models.fields.TextField', [], {}),
            'short_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'solvent': ('django.db.models.fields.TextField', [], {}),
            'start_plate': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'version_number': ('django.db.models.fields.IntegerField', [], {'default': '0'})
        },
        u'db.librarycontentsversion': {
            'Meta': {'object_name': 'LibraryContentsVersion', 'db_table': "u'library_contents_version'"},
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'library_contents_loading_activity': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'lcv_load'", 'to': u"orm['db.AdministrativeActivity']"}),
            'library_contents_release_activity': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'lcv_release'", 'null': 'True', 'to': u"orm['db.AdministrativeActivity']"}),
            'library_contents_version_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'version_number': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.libraryscreening': {
            'Meta': {'object_name': 'LibraryScreening', 'db_table': "u'library_screening'"},
            'abase_testset_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screening']", 'primary_key': 'True'}),
            'is_for_external_library_plates': ('django.db.models.fields.BooleanField', [], {}),
            'libraries_screened_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'library_plates_screened_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'screened_experimental_well_count': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.libraryupdateactivity': {
            'Meta': {'object_name': 'LibraryUpdateActivity', 'db_table': "u'library_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.molfile': {
            'Meta': {'object_name': 'Molfile', 'db_table': "u'molfile'"},
            'molfile': ('django.db.models.fields.TextField', [], {}),
            'reagent': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.Reagent']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'db.naturalproductreagent': {
            'Meta': {'object_name': 'NaturalProductReagent', 'db_table': "u'natural_product_reagent'"},
            'reagent': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.Reagent']", 'unique': 'True', 'primary_key': 'True'})
        },
        u'db.plate': {
            'Meta': {'object_name': 'Plate', 'db_table': "u'plate'"},
            'avg_remaining_volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'copy': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Copy']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'facility_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'max_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'max_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'max_remaining_volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'min_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'min_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'min_remaining_volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plate_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'plate_location': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.PlateLocation']", 'null': 'True', 'blank': 'True'}),
            'plate_number': ('django.db.models.fields.IntegerField', [], {}),
            'plate_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'plated_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True', 'null': 'True', 'blank': 'True'}),
            'primary_well_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'primary_well_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'quadrant': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'remaining_volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'retired_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True', 'null': 'True', 'blank': 'True'}),
            'screening_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.TextField', [], {}),
            'stock_plate_number': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well_volume': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'})
        },
        u'db.platelocation': {
            'Meta': {'object_name': 'PlateLocation', 'db_table': "u'plate_location'"},
            'bin': ('django.db.models.fields.TextField', [], {}),
            'freezer': ('django.db.models.fields.TextField', [], {}),
            'plate_location_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'room': ('django.db.models.fields.TextField', [], {}),
            'shelf': ('django.db.models.fields.TextField', [], {})
        },
        u'db.plateupdateactivity': {
            'Meta': {'object_name': 'PlateUpdateActivity', 'db_table': "u'plate_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'plate_id': ('django.db.models.fields.IntegerField', [], {}),
            'update_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'})
        },
        u'db.publication': {
            'Meta': {'object_name': 'Publication', 'db_table': "u'publication'"},
            'attached_file': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AttachedFile']", 'unique': 'True', 'null': 'True', 'blank': 'True'}),
            'authors': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'journal': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'pages': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'publication_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'pubmed_central_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pubmed_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'title': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'volume': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'year_published': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.reagent': {
            'Meta': {'object_name': 'Reagent', 'db_table': "u'reagent'"},
            'library_contents_version': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LibraryContentsVersion']", 'null': 'True'}),
            'reagent_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'substance_id': ('django.db.models.fields.CharField', [], {'default': "'UWPRH2ZJ'", 'unique': 'True', 'max_length': '8'}),
            'vendor_batch_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'vendor_identifier': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'vendor_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']", 'null': 'True'})
        },
        u'db.reagentpublicationlink': {
            'Meta': {'object_name': 'ReagentPublicationLink', 'db_table': "u'reagent_publication_link'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'publication_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"})
        },
        u'db.resultvalue': {
            'Meta': {'object_name': 'ResultValue', 'db_table': "u'result_value'"},
            'assay_well_control_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'data_column': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.DataColumn']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_exclude': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'is_positive': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'numeric_value': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'result_value_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']", 'null': 'True', 'blank': 'True'})
        },
        u'db.rnaicherrypickrequest': {
            'Meta': {'object_name': 'RnaiCherryPickRequest', 'db_table': "u'rnai_cherry_pick_request'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']", 'primary_key': 'True'})
        },
        u'db.schemahistory': {
            'Meta': {'object_name': 'SchemaHistory', 'db_table': "u'schema_history'"},
            'comment': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'date_updated': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'screensaver_revision': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'})
        },
        u'db.screen': {
            'Meta': {'object_name': 'Screen', 'db_table': "u'screen'"},
            'abase_protocol_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'abase_study_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amount_to_be_charged_for_screen': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '9', 'decimal_places': '2', 'blank': 'True'}),
            'assay_plates_screened_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'assay_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'billing_comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'billing_info_return_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'coms_approval_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'coms_registration_number': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'data_meeting_complete': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_meeting_scheduled': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_privacy_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_privacy_expiration_notified_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_sharing_level': ('django.db.models.fields.IntegerField', [], {}),
            'date_charged': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_completed5kcompounds': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_faxed_to_billing_department': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_of_application': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'facilities_and_administration_charge': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '9', 'decimal_places': '2', 'blank': 'True'}),
            'facility_id': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'fee_form_requested_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'fee_form_requested_initials': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'image_url': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_billing_for_supplies_only': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_fee_form_on_file': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'lab_head': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabHead']", 'null': 'True', 'blank': 'True'}),
            'lead_screener': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']", 'null': 'True', 'blank': 'True'}),
            'libraries_screened_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'library_plates_data_analyzed_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'library_plates_data_loaded_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'library_plates_screened_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'max_allowed_data_privacy_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'max_data_loaded_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'max_screened_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'min_allowed_data_privacy_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'min_data_loaded_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'min_screened_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'perturbagen_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'perturbagen_ug_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'pin_transfer_admin_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'null': 'True', 'blank': 'True'}),
            'project_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'project_phase': ('django.db.models.fields.TextField', [], {}),
            'pubchem_assay_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pubchem_deposited_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'publishable_protocol': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'publishable_protocol_comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'publishable_protocol_date_entered': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'publishable_protocol_entered_by': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'screen_id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen_type': ('django.db.models.fields.TextField', [], {}),
            'screened_experimental_well_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'see_comments': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'species': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'status': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'status_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'study_type': ('django.db.models.fields.TextField', [], {}),
            'summary': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'title': ('django.db.models.fields.TextField', [], {}),
            'to_be_requested': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'total_plated_lab_cherry_picks': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'transfection_agent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.TransfectionAgent']", 'null': 'True', 'blank': 'True'}),
            'unique_screened_experimental_well_count': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'url': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well_studied': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']", 'null': 'True', 'blank': 'True'})
        },
        u'db.screenbillingitem': {
            'Meta': {'object_name': 'ScreenBillingItem', 'db_table': "u'screen_billing_item'"},
            'amount': ('django.db.models.fields.DecimalField', [], {'max_digits': '9', 'decimal_places': '2'}),
            'date_sent_for_billing': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'item_to_be_charged': ('django.db.models.fields.TextField', [], {}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screenercherrypick': {
            'Meta': {'object_name': 'ScreenerCherryPick', 'db_table': "u'screener_cherry_pick'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            'screened_well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"}),
            'screener_cherry_pick_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.screenfundingsupportlink': {
            'Meta': {'object_name': 'ScreenFundingSupportLink', 'db_table': "u'screen_funding_support_link'"},
            'funding_support': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.FundingSupport']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screening': {
            'Meta': {'object_name': 'Screening', 'db_table': "u'screening'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabActivity']", 'primary_key': 'True'}),
            'assay_protocol': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_protocol_last_modified_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'assay_protocol_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_well_volume': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'number_of_replicates': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'volume_transferred_per_well_to_assay_plates': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'})
        },
        u'db.screeningroomuser': {
            'Meta': {'object_name': 'ScreeningRoomUser', 'db_table': "u'screening_room_user'"},
            'coms_crhba_permit_number': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'coms_crhba_permit_principal_investigator': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'lab_head': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabHead']", 'null': 'True', 'blank': 'True'}),
            'last_notified_rnaiua_checklist_item_event': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'rnai_ua_user'", 'null': 'True', 'to': u"orm['db.ChecklistItemEvent']"}),
            'last_notified_smua_checklist_item_event': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'smua_user'", 'null': 'True', 'to': u"orm['db.ChecklistItemEvent']"}),
            'screensaver_user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.ScreensaverUser']", 'unique': 'True', 'primary_key': 'True'}),
            'user_classification': ('django.db.models.fields.TextField', [], {})
        },
        u'db.screeningroomuserfacilityusagerole': {
            'Meta': {'object_name': 'ScreeningRoomUserFacilityUsageRole', 'db_table': "u'screening_room_user_facility_usage_role'"},
            'facility_usage_role': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screening_room_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"})
        },
        u'db.screenkeyword': {
            'Meta': {'object_name': 'ScreenKeyword', 'db_table': "u'screen_keyword'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'keyword': ('django.db.models.fields.TextField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screenpublicationlink': {
            'Meta': {'object_name': 'ScreenPublicationLink', 'db_table': "u'screen_publication_link'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'publication_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screenresult': {
            'Meta': {'object_name': 'ScreenResult', 'db_table': "u'screen_result'"},
            'channel_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'experimental_well_count': ('django.db.models.fields.IntegerField', [], {}),
            'replicate_count': ('django.db.models.fields.IntegerField', [], {}),
            'screen': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.Screen']", 'unique': 'True'}),
            'screen_result_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.screenresultupdateactivity': {
            'Meta': {'object_name': 'ScreenResultUpdateActivity', 'db_table': "u'screen_result_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen_result': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenResult']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.screensaveruser': {
            'Meta': {'object_name': 'ScreensaverUser', 'db_table': "u'screensaver_user'"},
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'digested_password': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'ecommons_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'email': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'first_name': ('django.db.models.fields.TextField', [], {}),
            'harvard_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'harvard_id_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'harvard_id_requested_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'last_name': ('django.db.models.fields.TextField', [], {}),
            'login_id': ('django.db.models.fields.TextField', [], {'unique': 'True', 'blank': 'True'}),
            'mailing_address': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'phone': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'screensaver_user_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']", 'null': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {'default': '1', 'blank': 'True'})
        },
        u'db.screensaveruserrole': {
            'Meta': {'object_name': 'ScreensaverUserRole', 'db_table': "u'screensaver_user_role'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']"}),
            'screensaver_user_role': ('django.db.models.fields.TextField', [], {})
        },
        u'db.screensaveruserupdateactivity': {
            'Meta': {'object_name': 'ScreensaverUserUpdateActivity', 'db_table': "u'screensaver_user_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.screenstatusitem': {
            'Meta': {'object_name': 'ScreenStatusItem', 'db_table': "u'screen_status_item'", 'index_together': "((u'screen', u'status', u'status_date'),)"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'status': ('django.db.models.fields.TextField', [], {}),
            'status_date': ('django.db.models.fields.DateField', [], {})
        },
        u'db.screenupdateactivity': {
            'Meta': {'object_name': 'ScreenUpdateActivity', 'db_table': "u'screen_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.serviceactivity': {
            'Meta': {'object_name': 'ServiceActivity', 'db_table': "u'service_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']", 'primary_key': 'True'}),
            'service_activity_type': ('django.db.models.fields.TextField', [], {}),
            'serviced_screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']", 'null': 'True', 'blank': 'True'}),
            'serviced_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"})
        },
        u'db.silencingreagent': {
            'Meta': {'object_name': 'SilencingReagent', 'db_table': "u'silencing_reagent'"},
            'anti_sense_sequence': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'duplex_wells': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['db.Well']", 'symmetrical': 'False'}),
            'facility_gene': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'facility_reagent'", 'unique': 'True', 'null': 'True', 'to': u"orm['db.Gene']"}),
            'is_restricted_sequence': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'reagent': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.Reagent']", 'unique': 'True', 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'silencing_reagent_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'vendor_gene': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'vendor_reagent'", 'unique': 'True', 'null': 'True', 'to': u"orm['db.Gene']"})
        },
        u'db.silencingreagentduplexwells': {
            'Meta': {'unique_together': "((u'silencing_reagent', u'well'),)", 'object_name': 'SilencingReagentDuplexWells', 'db_table': "u'silencing_reagent_duplex_wells'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'silencing_reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.SilencingReagent']"}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"})
        },
        u'db.smallmoleculechembankid': {
            'Meta': {'object_name': 'SmallMoleculeChembankId', 'db_table': "u'small_molecule_chembank_id'"},
            'chembank_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"})
        },
        u'db.smallmoleculechemblid': {
            'Meta': {'object_name': 'SmallMoleculeChemblId', 'db_table': "u'small_molecule_chembl_id'"},
            'chembl_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"})
        },
        u'db.smallmoleculecherrypickrequest': {
            'Meta': {'object_name': 'SmallMoleculeCherryPickRequest', 'db_table': "u'small_molecule_cherry_pick_request'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']", 'primary_key': 'True'})
        },
        u'db.smallmoleculecompoundname': {
            'Meta': {'object_name': 'SmallMoleculeCompoundName', 'db_table': "u'small_molecule_compound_name'"},
            'compound_name': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"})
        },
        u'db.smallmoleculepubchemcid': {
            'Meta': {'object_name': 'SmallMoleculePubchemCid', 'db_table': "u'small_molecule_pubchem_cid'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'pubchem_cid': ('django.db.models.fields.IntegerField', [], {}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"})
        },
        u'db.smallmoleculereagent': {
            'Meta': {'object_name': 'SmallMoleculeReagent', 'db_table': "u'small_molecule_reagent'"},
            'inchi': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_restricted_structure': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'molecular_formula': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'molecular_mass': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'molecular_weight': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'reagent': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['db.Reagent']", 'unique': 'True', 'primary_key': 'True'}),
            'smiles': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.studyreagentlink': {
            'Meta': {'object_name': 'StudyReagentLink', 'db_table': "u'study_reagent_link'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"}),
            'study': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.substance': {
            'Meta': {'object_name': 'Substance'},
            'comment': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.transfectionagent': {
            'Meta': {'object_name': 'TransfectionAgent', 'db_table': "u'transfection_agent'"},
            'transfection_agent_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.well': {
            'Meta': {'object_name': 'Well', 'db_table': "u'well'"},
            'deprecation_admin_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'null': 'True', 'blank': 'True'}),
            'facility_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_deprecated': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'latest_released_reagent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'reagent_well'", 'null': 'True', 'to': u"orm['db.Reagent']"}),
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'library_well_type': ('django.db.models.fields.TextField', [], {}),
            'mg_ml_concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'molar_concentration': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'plate_number': ('django.db.models.fields.IntegerField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'well_id': ('django.db.models.fields.TextField', [], {'primary_key': 'True'}),
            'well_name': ('django.db.models.fields.TextField', [], {})
        },
        u'db.wellvolumeadjustment': {
            'Meta': {'object_name': 'WellVolumeAdjustment', 'db_table': "u'well_volume_adjustment'"},
            'copy': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Copy']"}),
            'lab_cherry_pick': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabCherryPick']", 'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'volume': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"}),
            'well_volume_adjustment_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'well_volume_correction_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.WellVolumeCorrectionActivity']", 'null': 'True', 'blank': 'True'})
        },
        u'db.wellvolumecorrectionactivity': {
            'Meta': {'object_name': 'WellVolumeCorrectionActivity', 'db_table': "u'well_volume_correction_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'primary_key': 'True'})
        },
        u'reports.apilog': {
            'Meta': {'unique_together': "(('ref_resource_name', 'key', 'date_time'),)", 'object_name': 'ApiLog'},
            'added_keys': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'api_action': ('django.db.models.fields.CharField', [], {'max_length': '10'}),
            'comment': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'date_time': ('django.db.models.fields.DateTimeField', [], {}),
            'diff_keys': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'diffs': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json_field': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'parent_log': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'child_logs'", 'null': 'True', 'to': u"orm['reports.ApiLog']"}),
            'ref_resource_name': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'removed_keys': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'uri': ('django.db.models.fields.TextField', [], {}),
            'user_id': ('django.db.models.fields.IntegerField', [], {}),
            'username': ('django.db.models.fields.CharField', [], {'max_length': '128'})
        },
        u'reports.job': {
            'Meta': {'object_name': 'Job'},
            'comment': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'date_time_fullfilled': ('django.db.models.fields.DateTimeField', [], {'null': 'True'}),
            'date_time_processing': ('django.db.models.fields.DateTimeField', [], {'null': 'True'}),
            'date_time_requested': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'input_filename': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'path_info': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'remote_addr': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'request_method': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'response_code': ('django.db.models.fields.IntegerField', [], {}),
            'response_content': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'response_filename': ('django.db.models.fields.TextField', [], {'null': 'True'})
        },
        u'reports.listlog': {
            'Meta': {'unique_together': "(('apilog', 'ref_resource_name', 'key', 'uri'),)", 'object_name': 'ListLog'},
            'apilog': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['reports.ApiLog']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'ref_resource_name': ('django.db.models.fields.CharField', [], {'max_length': '35'}),
            'uri': ('django.db.models.fields.TextField', [], {})
        },
        u'reports.metahash': {
            'Meta': {'unique_together': "(('scope', 'key'),)", 'object_name': 'MetaHash'},
            'alias': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json_field': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'json_field_type': ('django.db.models.fields.CharField', [], {'max_length': '128', 'null': 'True', 'blank': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            'linked_field_type': ('django.db.models.fields.CharField', [], {'max_length': '128', 'null': 'True', 'blank': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'scope': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'})
        },
        u'reports.permission': {
            'Meta': {'unique_together': "(('scope', 'key', 'type'),)", 'object_name': 'Permission'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            'scope': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '15'})
        },
        u'reports.record': {
            'Meta': {'object_name': 'Record'},
            'base_value1': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'scope': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'})
        },
        u'reports.recordmultivalue': {
            'Meta': {'unique_together': "(('field_meta', 'parent', 'ordinal'),)", 'object_name': 'RecordMultiValue'},
            'field_meta': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['reports.MetaHash']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['reports.Record']"}),
            'value': ('django.db.models.fields.TextField', [], {})
        },
        u'reports.recordvalue': {
            'Meta': {'object_name': 'RecordValue'},
            'field_meta': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['reports.MetaHash']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['reports.Record']"}),
            'value': ('django.db.models.fields.TextField', [], {'null': 'True'})
        },
        u'reports.recordvaluecomplex': {
            'Meta': {'object_name': 'RecordValueComplex'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['reports.Record']", 'unique': 'True'}),
            'value1': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'value2': ('django.db.models.fields.TextField', [], {'null': 'True'})
        },
        u'reports.usergroup': {
            'Meta': {'object_name': 'UserGroup'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['reports.Permission']", 'symmetrical': 'False'}),
            'super_groups': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'sub_groups'", 'symmetrical': 'False', 'to': u"orm['reports.UserGroup']"}),
            'users': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['reports.UserProfile']", 'symmetrical': 'False'})
        },
        u'reports.userprofile': {
            'Meta': {'object_name': 'UserProfile'},
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'ecommons_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'gender': ('django.db.models.fields.CharField', [], {'max_length': '15', 'null': 'True'}),
            'harvard_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'harvard_id_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'harvard_id_requested_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json_field': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'json_field_type': ('django.db.models.fields.CharField', [], {'max_length': '128', 'null': 'True', 'blank': 'True'}),
            'mailing_address': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['reports.Permission']", 'symmetrical': 'False'}),
            'phone': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'}),
            'username': ('django.db.models.fields.TextField', [], {})
        },
        u'reports.vocabularies': {
            'Meta': {'unique_together': "(('scope', 'key'),)", 'object_name': 'Vocabularies'},
            'alias': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json_field': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '128', 'blank': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'scope': ('django.db.models.fields.CharField', [], {'max_length': '128', 'blank': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '512', 'blank': 'True'})
        }
        
    }

    complete_apps = ['reports','db']
    
    