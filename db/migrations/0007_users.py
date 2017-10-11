# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
from datetime import timedelta
import datetime
import json
import logging
import os

from django.contrib.auth.models import User, UserManager
from django.db import migrations, models
import pytz

from db.support.data_converter import default_converter
from lims.base_settings import PROJECT_ROOT
from reports.models import ApiLog
from db.api import _now


logger = logging.getLogger(__name__)

    
def create_screensaver_users(apps, schema_editor):
    '''
    Create entries in auth_user, reports.UserProfile for each valid entry 
    in the screensaver_user table.
    - 'valid' entries have [ecommons_id or login_id] and email addresses.
    - Duplicated ecommons/and/or/login ID's are invalid; the second
    (duplicate) entry will be ignored.
    Not idempotent - can be re-run by deleting the reports -
    Usergroup_Users, and UserProfile tables, and also the Auth_User
    entries - (make sure not to delete the superuser).
    '''
    
    count=0
    
    AuthUserClass = apps.get_model('auth','User')
    UserProfileClass = apps.get_model('reports','UserProfile')
    ScreensaverUser = apps.get_model('db', 'ScreensaverUser')
    ScreensaverUserRole = apps.get_model('db', 'ScreensaverUserRole')
    auth_user_username_limit = 30 # field size limit for auth_user
    
    for su in ScreensaverUser.objects.all():
        logger.info('processing: %r: %r', su, su.screensaver_user_id)
        if su.screensaver_user_id == 4712:
            logger.info('skip erroneous duplicate user account: %r', su)
            continue
        if su.screensaver_user_id == 830:
            logger.info('skip Ruchir Shah rcs12 old account: %r', su)
            continue
        
        logger.info("processing %r, ecommons: %r, login_id: %r, email: %r, %s,%s",
            su.screensaver_user_id, su.ecommons_id, su.login_id, su.email, 
            su.first_name, su.last_name)
        au = None
        up = None
        username = None
        if su.ecommons_id: 
            # convert in case it has an error
            username = default_converter(str(su.ecommons_id)) 
            logger.info('username: converted ecommons: %r to %r', 
                su.ecommons_id, username)
        elif su.login_id:
            username = su.login_id
        else:
            continue
        
        su.username = username
        logger.info('save new su.username: %r, %r', su.screensaver_user_id, su.username)
        su.save()
        
        has_login = ScreensaverUserRole.objects.all()\
            .filter(screensaver_user_id=su.screensaver_user_id)\
            .filter(screensaver_user_role='screensaverUser').exists()
        
        if has_login is not True:
            continue
        
        try:
            au = AuthUserClass.objects.get(username=username)
            logger.info('found auth_user: %s', username)
        except Exception, e:
            au = AuthUserClass(
                username = username, 
                email = su.email if su.email else 'none', 
                first_name = su.first_name, 
                last_name = su.last_name,
                date_joined = datetime.datetime.now().replace(
                    tzinfo=pytz.timezone('US/Eastern')),
                is_active=False,
                is_staff=False)

        if hasattr(su, 'administratoruser'):
            au.is_staff = True
        au.save()
        
        try:
            up = UserProfileClass.objects.get(username=username)
            logger.info('found userprofile: %s', username)
        except Exception, e:
            if su.created_by: 
                created_by_username = su.created_by.username
            else:
                created_by_username = 'sde_EDIT'
                
            up = UserProfileClass()
            
            up.username = username 
            up.phone = su.phone
            up.mailing_address = su.mailing_address
            up.comments = su.comments
            up.ecommons_id = su.ecommons_id
            up.harvard_id = su.harvard_id
            up.harvard_id_expiration_date = \
                su.harvard_id_expiration_date
            up.harvard_id_requested_expiration_date = \
                su.harvard_id_requested_expiration_date
            up.created_by_username = created_by_username            
            
            up.user = au
            up.save()

        su.user = up
        logger.info('su: (%r, %r, %r, %r, %r, %r), up: %r, au: %r, username: %r',
            su.screensaver_user_id, su.ecommons_id, su.login_id, su.email, 
            su.first_name,su.last_name, up.id, au.id, au.username)
        su.save()
        
        logger.info('user created: %r', au.username)
        count += 1
    logger.info('Converted %d users with login access', count)
    
def create_user_checklists(apps, schema_editor):
    '''
    Convert ChecklistItemEvent entries into UserChecklist
    - create ApiLogs to track status changes
    - also track the status_notified
    
    - not idempotent, can be re-run by deleting UCIs and Apilogs
    '''
    # prerequisites: 
    # - convert checklist_item / checklist_item_event entries into into 
    # checklistitem.* vocabularies (migration 0003)
    # - create the user_checklist_item table (0002)

    ChecklistItem = apps.get_model('db','ChecklistItem')
    UserChecklist = apps.get_model('db','UserChecklist')
    
    uc_name_map = {}
    for obj in ChecklistItem.objects.all().distinct('item_name'):
        key = default_converter(obj.item_name)
        uc_name_map[obj.item_name] = key

    # create entries in the user_checklist table
    # note: status values are hard-coded to correspond to the vocabulary
    # keys (created in migration 0003)
    sql_keys = [
        'checklist_item_event_id', 'suid','cigroup','ciname',
        'su_username','admin_username','admin_suid','admin_upid',
        'date_performed', 'date_created','status','is_notified'
        ]
    sql = '''
select
cie.checklist_item_event_id,
screening_room_user_id,
ci.checklist_item_group,
ci.item_name,
su.username su_username,
admin.username admin_username,
admin.screensaver_user_id admin_suid,
up.id admin_upid,
cie.date_performed,
cie.date_created,
case 
    when cie.is_not_applicable then 'n_a'
    when ci.is_expirable and cie.date_performed is not null then
        case 
            when cie.is_expiration then 'deactivated' 
            else 'activated' end
    when cie.date_performed is not null then 'completed'     
    else 'not_completed'
end as status,
( 
   select 1 from screening_room_user sru 
    where sru.last_notified_smua_checklist_item_event_id = cie.checklist_item_event_id
       UNION    
   select 1 from screening_room_user sru 
    where sru.last_notified_rnaiua_checklist_item_event_id = cie.checklist_item_event_id
) as is_notified
from checklist_item ci
join checklist_item_event cie using(checklist_item_id)
join screensaver_user su on screening_room_user_id=su.screensaver_user_id
join screensaver_user admin on cie.created_by_id=admin.screensaver_user_id
left join reports_userprofile up on up.id=admin.user_id
order by screening_room_user_id, checklist_item_group, item_name, cie.date_performed asc;
'''
    now = datetime.datetime.now()
    now = pytz.timezone('US/Eastern').localize(now)
    earliest_allowed_time = datetime.datetime(2000, 1, 1)
    earliest_allowed_time = pytz.timezone('US/Eastern').localize(earliest_allowed_time)
    
    connection = schema_editor.connection
    cursor = connection.cursor()
    resource_name = 'userchecklist'
    _dict = None
    log = None
    ucl_hash = {}
    notified_ucl_hash = {}
    unique_log_keys = set()

    try:
        cursor.execute(sql)
        i = 0
        
        # Iterate through the ChecklistItemEvents:
        # - Ordered by date_performed asc
        # - create a UserChecklist - key [username,checklist_item_name] 
        # - first occurrence creates a new UCL
        # - subsequent occurrences represent updates
        # - keep track of UCLs in hash
        # - look for previous UCL 
        
        for row in cursor:
            _dict = dict(zip(sql_keys,row))
            ucl = None
            checklist_name = uc_name_map[_dict['ciname']]
            key = '/'.join([str(_dict['suid']),checklist_name])
            previous_dict = ucl_hash.get(key, None)
            notified_previous_dict = notified_ucl_hash.get(key, None)
            logger.debug('previous_dict: %s:%s' % (key,previous_dict))

            date_created = _dict['date_created']
            if date_created.tzinfo is None:
                date_created = pytz.timezone('US/Eastern').localize(date_created)
            date_time = date_created
            if date_created.date() != _dict['date_performed']:
                # only use the less accurate date_performed date if that date
                # is not equal to the date_created date
                date_performed = pytz.timezone('US/Eastern').localize(
                    datetime.datetime.combine(
                        _dict['date_performed'],
                        datetime.datetime.min.time()))
                if date_performed < now:
                    if date_performed > earliest_allowed_time:
                        date_time = date_performed
                else:
                    logger.info('date performed is erroneous: > than today: %r, %r',
                        date_performed, _dict)
            if date_time < earliest_allowed_time:
                logger.info('erroneous times: date_created: %r, date_performed: %r', 
                    date_created, date_performed)
                date_time = earliest_allowed_time
            if date_time > now:
                logger.error('log time > now: %r, %r, %r', date_time, now, _dict)
                raise Exception('log time > now')
            # Set up Log
            # create the apilog for this item
            log = ApiLog()
            log.ref_resource_name = resource_name
            log.username = _dict['admin_username']
            log.user_id = _dict['suid']
            log.date_time = date_time
            log.api_action = 'PATCH'
            log.key = key
            log.uri = '/'.join([log.ref_resource_name,log.key])
            log.json_field = { 
                'migration': 'ChecklistItemEvent',
                'data': { 'checklist_item_event_id': 
                    _dict['checklist_item_event_id'] }          
                }
            if log.username is None:
                log.username = 'sde_EDIT'
            logger.debug('creating log: %r', log)
            # For logging: is the key (date_time, actually) unique?
            full_key = '/'.join([log.ref_resource_name,log.key,str(log.date_time)])
            while full_key in unique_log_keys:
                # add a second to make it unique; because date performed is a date,
                logger.info('time collision for: %r',full_key)
                log.date_time = log.date_time  + datetime.timedelta(milliseconds=1)
                logger.info('new log time: %r', log.date_time)
                full_key = '/'.join([log.ref_resource_name,log.key,str(log.date_time)])
            unique_log_keys.add(full_key)
                
            if previous_dict:
                # NOTE: for SMUA, there may be multiple events before "deactivation" event;
                # - so the previous dict will be updated multiple times here,
                # - 60 days from the last
                ucl = previous_dict['obj']
                ucl.admin_user_id = int(_dict['admin_suid'])
                previous_state = ucl.is_checked
                ucl.is_checked = False
                if _dict['status'] in ['activated', 'completed']:
                    ucl.is_checked = True
                ucl.date_effective = _dict['date_performed']

                if _dict['status'] == 'deactivated':
                    if notified_previous_dict:
                        ucl.date_notified = (
                            _dict['date_performed'] - datetime.timedelta(days=60))
                
                logger.debug(
                    'saving, dict: %s, prev_dict: %s, date_effective %s, date_notified: %s', 
                    _dict, previous_dict, ucl.date_effective, ucl.date_notified)
                ucl.save()
                
                # Make a log
                log.diffs['is_checked'] = [previous_state,ucl.is_checked]
                if previous_dict['admin_username'] != _dict['admin_username']:    
                    log.diffs['admin_username'] = \
                        [previous_dict['admin_username'], _dict['admin_username']]
                log.diffs['status_date'] = [
                    previous_dict['date_performed'].isoformat(), 
                    _dict['date_performed'].isoformat()]
                log.diffs['status'] = [previous_dict['status'],_dict['status']]
            else:
                # create
                ucl_hash[key] = _dict
                if _dict['is_notified']:
                    notified_ucl_hash[key] = _dict
                logger.debug('create user checklist item: %r, %r', _dict, 
                    _dict['date_performed'].isoformat())

                is_checked = False
                if _dict['status'] in ['activated', 'completed']:
                    is_checked = True
                ucl = UserChecklist.objects.create(
                    screensaver_user_id = int(_dict['suid']),
                    admin_user_id = int(_dict['admin_suid']),
                    name = checklist_name,
                    is_checked = is_checked,
                    date_effective = _dict['date_performed'])
                ucl.save()
                _dict['obj'] = ucl
                
                # Make a log
                log.api_action = 'CREATE'
                log.diffs['is_checked'] = [None,ucl.is_checked]
                log.diffs['date_effective'] = [
                    None,_dict['date_performed'].isoformat()]
                
                logger.debug('created ucl: %r', ucl)
                i += 1

            log.save()

            logger.debug('created log: %r, %r', log, log.diffs )
            if i%1000 == 0:
                logger.info('created %d logs', i)
                logger.info('key: %r', key)
                logger.info('created log: %r, %r', log, log.diffs )
    
    except Exception, e:
        logger.exception('migration exc')
        raise e  
    
    logger.info('created %d user_checklist_items', i)
    
    
class Migration(migrations.Migration):

    dependencies = [
        ('db', '0004_postprep'),
        ('auth', '0001_initial')
    ]

    operations = [
        
        
        # Note: for users migration;
        # - the fields last_notified_smua_checklist_item_event, (also rnai)
        # will be moved to the userchecklistitem status_notified_date;
        # And nightly batch scripts will have to be modified accordingly
        # (property needed on the vocabulary/checklistitem: 
        # - "expiration_interval","expire_notify_days" properties,
        # - state machine here: i.e. status ordering, 
        # --allowed status transitions??  
        # --(note also, should have an "available_status_types" to restrict the 
        # status-states that the item can move through, using this to run the 
        # expiration reports for nightly automatic expirations.
        #     not_completed
        #     activated
        #     deactivated 
        #     na
        #     completed

        migrations.RunPython(create_screensaver_users),
        migrations.RunPython(create_user_checklists),
    ]
