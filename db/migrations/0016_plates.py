
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

# import datetime
import json
import logging

from django.db import migrations, models
# from pytz import timezone
# import pytz

from db.support.data_converter import default_converter
from reports.models import ApiLog
from collections import OrderedDict
import re
from decimal import Decimal

from db.migrations import create_log_time, _create_generic_log,\
    _child_log_from

logger = logging.getLogger(__name__)

base_uri = '/db/api/v1'
copywell_resource_name = 'copywell'
copywell_uri = '/db/api/v1/' + copywell_resource_name
librarycopyplate_resource_name = 'librarycopyplate'


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
    
    # TODO: track the plate.cplt_screening_count
    
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
        migrations.RunPython(create_plate_activity_logs),
        migrations.RunPython(create_copywell_adjustments),
    ]
