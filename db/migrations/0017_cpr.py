
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from collections import OrderedDict
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
copywell_resource_name = 'copywell'
copywell_uri = '/db/api/v1/' + copywell_resource_name
librarycopyplate_resource_name = 'librarycopyplate'


def _create_plate_activity_log(activity_dict):
    log = ApiLog()
    log.ref_resource_name = librarycopyplate_resource_name
    log.key = '/'.join([
        activity_dict['library_short_name'],activity_dict['copy_name'],
        str(int(activity_dict['plate_number'])).zfill(5)])
    log.uri = '/'.join([base_uri,log.ref_resource_name,log.key])
    log.comment = activity_dict['comments']
    log.date_time = create_log_time(activity_dict['date_of_activity'])
    log.username = activity_dict['username']
    log.user_id = activity_dict['screensaver_user_id']
    return log

def migrate_plate_generic_logs(apps, schema_editor):
    logger.info('create plate generic activity logs')

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
        '''
    from activity a 
    join screensaver_user on(performed_by_id=screensaver_user_id)
    join plate_update_activity on (activity_id=update_activity_id) 
    join plate using(plate_id)
    join copy using(copy_id)
    join library using(library_id)
    where
    not exists (select null from plate where plated_activity_id = activity_id) 
    and not exists (select null from plate where retired_activity_id = activity_id) 
    and a.comments not ilike '%Facility ID changed%' 
    and a.comments not ilike '%location changed from%' 
    and a.comments not ilike '%well volume changed%' 
    and a.comments not ilike '%plate type changed%' 
    and a.comments not ilike '%concentration changed from%' 
    and a.comments not ilike '%concentration (molar) changed from%' 
    and a.comments not ilike '%concentration (mg/ml) changed from%' 
    and a.comments not ilike '%activity created by database migration%'
    and a.comments not ilike '%Freezer 1 inventory updates for Richard Siu and Rachel Warden%'
    '''
    )
    
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
            logger.info('plate generic log: %r, %r', log, log.comment)
            logger.info('processed %d plate generic activity logs', i)
    
    
def create_plate_plated_and_retired_logs(apps, schema_editor):  

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
        log = _create_plate_activity_log(_activity)
        if "'available'" in log.comment.lower():
            log.diffs['status'] = ['not_specied','available']
        elif "'not available'" in log.comment.lower():
            log.diffs['status'] = ['not_specied','not_available']
        else:
            raise Exception('unknown plate.plated_activity comment: %r', _activity)
        log.diffs['date_plated'] = [None, _activity['date_of_activity']]
        log.json_field = { 
            'migration': 'LibraryCopyPlate plating',
            'data': { 
                'plate_update_activity.activity_id': _activity['activity_id'] }
        }
        log.save()
        
        if i % 1000 == 0:
            logger.info('plate plated log: %r, %r', log, log.comment)
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
                'plate_update_activity.activity_id': _activity['activity_id'] }
        }
        log.save()
        status_terms_recognized.add(default_converter(match.group(1)))
        status_terms_recognized.add(default_converter(match.group(2)))
        
        if i % 1000 == 0:
            logger.info('plate retired log: %r, %r', log, log.comment)
            logger.info('processed %d plate retired activity logs', i)
    logger.info('processed %d plate retired activity logs', i)
    logger.info('status terms recognized: %r', status_terms_recognized)

def create_library_screening_logs(apps, schema_editor):
     
    logger.info('create library screening logs')
     
    LibraryScreening = apps.get_model('db', 'LibraryScreening')
    ScreeensaverUser = apps.get_model('db', 'ScreensaverUser')
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
        screensaver_user_performed_by = activity.performed_by
        screensaver_user_created_by = activity.created_by
        # NOTE: if activity.date_performed < 2010-01-15, the "created_by" field
        # is empty
        default_admin_user = ScreeensaverUser.objects.get(username='sr50')
        if screensaver_user_created_by is None:
            screensaver_user_created_by = default_admin_user
        
        # NOTE: not creating a screen "parent" log;
        # libraryScreening activities will stand on their own
        # Create a LibraryScreening log
        library_screening_log = _create_generic_log(activity)
        library_screening_log.ref_resource_name = 'libraryscreening'
        library_screening_log.key = str(activity.activity_id)
        library_screening_log.uri = '/'.join([
            'screen', screen.facility_id, 
            library_screening_log.ref_resource_name,library_screening_log.key])
        library_screening_log.diffs = {
            'screened_experimental_well_count': 
                [None,screening.screened_experimental_well_count],
            'libraries_screened_count': [None, screening.libraries_screened_count],
            'library_plates_screened_count': [None, screening.library_plates_screened_count],
            'is_for_external_library_plates': [None, screening.is_for_external_library_plates]
            }
        library_screening_log.json_field = {
            'migration': 'LibraryScreening',
            'data': { 'library_screening.activity_id': activity.activity_id 
            }
        }
        library_screening_log.save()
         
        j = 0
        cp_logs = []
        
        for assay_plate in screening.assayplate_set.all().filter(replicate_ordinal=0):
            plate = assay_plate.plate
            # copyplate log for library screening
            cp_log = _child_log_from(library_screening_log)
            cp_log.comment = None
            cp_log.ref_resource_name = librarycopyplate_resource_name 
            cp_log.key = '/'.join([
                plate.copy.library.short_name, plate.copy.name, 
                str(plate.plate_number).zfill(5)])
            cp_log.uri = '/'.join([base_uri,cp_log.ref_resource_name,cp_log.key])
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
            cp_log.json_field = {
                'migration': 'LibraryScreening',
                'data': { 'library_screening.activity_id': activity.activity_id 
                }
            }
            screen_count +=1
            copyplate_to_screening_count[cp_log.key] = screen_count
            copyplate_to_volume[cp_log.key] = new_volume
            cp_log.json_field = json.dumps({
                'migration: volume_transferred_per_well_from_library_plates': str(adjustment)
                })
            cp_logs.append(cp_log)
            j += 1
 
        # TODO use bulk create
        for cp_log in cp_logs:
            cp_log.save()
 
        i += 1
        total_plate_logs += j
        if total_plate_logs % 1000  == 0:
            logger.info('library screening log: %r, %r, %d', 
                library_screening_log, library_screening_log.comment, len(cp_logs))
            logger.info('created screen logs: %d, plate logs: %d',i, total_plate_logs)
    logger.info('created screen logs: %d, plate logs: %d',i, total_plate_logs)


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0016_plates'),
    ]

    operations = [
        migrations.RunPython(create_library_screening_logs),
        migrations.RunPython(create_plate_plated_and_retired_logs),
        migrations.RunPython(migrate_plate_generic_logs),
    ]
    
