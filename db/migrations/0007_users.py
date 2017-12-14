# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
# from datetime import timedelta
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
from django.db.utils import ProgrammingError
import re

RESOURCE_USER_CHECKLIST = 'userchecklist'
RESOURCE_USER_AGREEMENT = 'useragreement'

CHECKLIST_NAME_RNA_UA = 'current_rnai_user_agreement_active'
CHECKLIST_NAME_SM_UA = 'current_small_molecule_user_agreement_active'
checklist_to_agreement_map = {
    CHECKLIST_NAME_RNA_UA: 'rnai',
    CHECKLIST_NAME_SM_UA: 'sm',
}

logger = logging.getLogger(__name__)


unique_log_keys = set()
now = datetime.datetime.now()
now = pytz.timezone('US/Eastern').localize(now)
earliest_allowed_time = datetime.datetime(2000, 1, 1)
earliest_allowed_time = pytz.timezone('US/Eastern').localize(earliest_allowed_time)
    
def set_log_time(log, date_time):
    ''' Ensure that the full log key (date_time) is unique
    '''

    log.date_time = date_time
    full_key = '/'.join([log.ref_resource_name,log.key,str(log.date_time)])
    while full_key in unique_log_keys:
        # add a second to make it unique; because date performed is a date,
        logger.debug('time collision for: %r',full_key)
        log.date_time = log.date_time  + datetime.timedelta(milliseconds=1)
        logger.debug('new log time: %r', log.date_time)
        full_key = '/'.join([log.ref_resource_name,log.key,str(log.date_time)])
    unique_log_keys.add(full_key)

def get_activity_date_time(date_time_created, date_performed):
    if date_time_created.tzinfo is None:
        date_time_created = pytz.timezone('US/Eastern').localize(date_time_created)
    date_time = date_time_created
    if date_time_created.date() != date_performed:
        # only use the less accurate date_performed date if that date
        # is not equal to the date_time_created date
        date_performed = pytz.timezone('US/Eastern').localize(
            datetime.datetime.combine(
                date_performed,
                datetime.datetime.min.time()))
        if date_performed < now:
            if date_performed > earliest_allowed_time:
                date_time = date_performed
        else:
            logger.warn('date performed %r is erroneous: > than today, date_time_created %r',
                date_performed, date_time_created)
    if date_time < earliest_allowed_time:
        logger.warn('erroneous times: date_time_created: %r, date_performed: %r', 
            date_time_created, date_performed)
        date_time = earliest_allowed_time
    if date_time > now:
        logger.error('log time > now: %r, %r, %r, %r', date_time, now, 
            date_time_created, date_performed)
        raise Exception
    
    return date_time

    
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
    
def create_user_checklist_from_checklist_item_events(apps, schema_editor):
    '''
    Convert ChecklistItemEvent entries into UserChecklist
    - create ApiLogs to track status changes
    - also track the status_notified
    
    - not idempotent, can be re-run by deleting user_checklist, user_agreement
    and reports_apilog/reports_logdiff;
        /* clear for new migration 0007 */
        
        delete from reports_logdiff where exists( select null from reports_apilog where ref_resource_name = 'userchecklist' and log_id=id);
        delete from reports_apilog where ref_resource_name = 'userchecklist';
        delete from user_checklist ;
        delete from reports_logdiff where exists( select null from reports_apilog where ref_resource_name = 'useragreement' and log_id=id);
        delete from reports_apilog where ref_resource_name = 'useragreement';
        delete from user_agreement ;
    
    '''
    # prerequisites: 
    # - convert checklist_item / checklist_item_event entries into into 
    # checklistitem.* vocabularies (migration 0003)
    # - create the user_checklist_item table (0002)

    ChecklistItem = apps.get_model('db','ChecklistItem')
    UserChecklist = apps.get_model('db','UserChecklist')
    UserAgreement = apps.get_model('db','UserAgreement')
    
    # Create a map from ci names to new names:
    uc_name_map = {}
    for obj in ChecklistItem.objects.all().distinct('item_name'):
        key = default_converter(obj.item_name)
        uc_name_map[obj.item_name] = key
    logger.info('uc_name_map: %r', uc_name_map)
    
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
left outer join screensaver_user admin on cie.created_by_id=admin.screensaver_user_id
left outer join reports_userprofile up on up.id=admin.user_id
order by screening_room_user_id, checklist_item_group, item_name, checklist_item_event_id asc;
'''
    connection = schema_editor.connection
    cursor = connection.cursor()
    resource_name = 'userchecklist'
    _dict = None
    log = None
    ucl_hash = {}
    notified_ucl_hash = {}

    cursor.execute(sql)
    i = 0
    
    # Iterate through the ChecklistItemEvents:
    # - Ordered by date_performed asc
    # - create a UserChecklist - key [username,checklist_item_name] 
    # - first occurrence creates a new UserChecklist or UserAgreement
    # - subsequent occurrences represent updates
    # - keep track of UCLs in hash
    # - look for previous UCL 
    
    for row in cursor:
        
        _dict = dict(zip(sql_keys,row))
        checklist_name = uc_name_map[_dict['ciname']]
        if checklist_name in checklist_to_agreement_map:
            resource_name = RESOURCE_USER_AGREEMENT
            key = '/'.join([
                str(_dict['suid']),
                checklist_to_agreement_map[checklist_name] ])
        else:
            resource_name = RESOURCE_USER_CHECKLIST
            key = '/'.join([str(_dict['suid']),checklist_name])
        previous_dict = ucl_hash.get(key, None)
        notified_previous_dict = notified_ucl_hash.get(key, None)

        logger.debug('previous_dict: %s:%s' % (key,previous_dict))
        
        # Create a log for every event
        log = ApiLog()
        log.ref_resource_name = resource_name
        log.username = _dict['admin_username']
        log.user_id = _dict['suid']
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
        
        date_time = get_activity_date_time(_dict['date_created'],_dict['date_performed'])
        set_log_time(log, date_time)
            
        if previous_dict:
            # NOTE: for SMUA, there may be multiple events before the
            # "deactivation" event;
            # - so the previous dict will be updated multiple times here,
            # - 60 days from the last
            if 'ucl' in previous_dict:
                ucl = previous_dict['ucl']
                ucl.admin_user_id = int(_dict['admin_suid'])
                previous_state = ucl.is_checked
                ucl.is_checked = False
                if _dict['status'] in ['activated', 'completed']:
                    ucl.is_checked = True
                elif _dict['status'] == 'deactivated':
                    if notified_previous_dict:
                        ucl.date_notified = (
                            _dict['date_performed'] - datetime.timedelta(days=60))
                        log.diffs['date_notified'] = [
                            None, ucl.date_notified.isoformat()]
                ucl.date_effective = _dict['date_performed']
                
                logger.debug(
                    'dict: %s, prev_dict: %s, date_effective %s, date_notified: %s', 
                    _dict, previous_dict, ucl.date_effective, ucl.date_notified)
                
                # Make a log
                log.diffs['is_checked'] = [previous_state,ucl.is_checked]
                if previous_dict['admin_username'] != _dict['admin_username']:    
                    log.diffs['admin_username'] = \
                        [previous_dict['admin_username'], _dict['admin_username']]
                log.diffs['status_date'] = [
                    previous_dict['date_performed'].isoformat(), 
                    _dict['date_performed'].isoformat()]
                log.diffs['status'] = [previous_dict['status'],_dict['status']]
                ucl.save()
            elif 'ua' in previous_dict:
                user_agreement = previous_dict['ua']
                if _dict['status'] in ['activated', 'completed']:
                    user_agreement.date_active = _dict['date_performed']
                    previous_expired = user_agreement.date_expired
                    user_agreement.date_expired = None
                    previous_notified = user_agreement.date_notified
                    user_agreement.date_notified = None
                    # NOTE: implied that UA has been nulled out before reactivating
                    log.diffs['date_active'] = \
                        [None,_dict['date_performed'].isoformat()]
                    if previous_expired is not None:
                        log.diffs['date_expired'] = \
                            [previous_expired.isoformat(), None]
                    if previous_notified is not None:
                        log.diffs['date_notified'] = \
                            [previous_notified.isoformat(), None]
                        
                if _dict['status'] == 'deactivated':
                    user_agreement.date_expired = _dict['date_performed']
                    # NOTE: implied that UA has been nulled out 
                    # before reactivating/deactivating
                    log.diffs['date_expired'] = \
                        [None,_dict['date_performed'].isoformat()]
                    if notified_previous_dict:
                        user_agreement.date_notified = (
                            _dict['date_performed'] - datetime.timedelta(days=60))
                        log.diffs['date_notified'] = [
                            None, user_agreement.date_notified.isoformat()]
                user_agreement.save()
                
            else:
                logger.error(
                    'no obj found in prev dict: %r', previous_dict)
                raise ProgrammingError
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
                
            if checklist_name in checklist_to_agreement_map:
                user_agreement = UserAgreement.objects.create(
                    screensaver_user_id = int(_dict['suid']),
                    type=checklist_to_agreement_map[checklist_name],
                )
                log.api_action = 'CREATE'
                
                if is_checked is True:
                    user_agreement.date_active=_dict['date_performed']
                    log.diffs['date_active'] = [
                        None,_dict['date_performed'].isoformat()]
                else:
                    logger.warn('first ua record is expiration: %r', _dict)
                    user_agreement.date_expired=_dict['date_performed']
                    log.diffs['date_expired'] = [
                        None,_dict['date_performed'].isoformat()]
                
                user_agreement.save()
                _dict['ua'] = user_agreement

                logger.debug('created ua: %r', user_agreement)
                
            else:    
                ucl = UserChecklist.objects.create(
                    screensaver_user_id = int(_dict['suid']),
                    admin_user_id = int(_dict['admin_suid']),
                    name = checklist_name,
                    is_checked = is_checked,
                    date_effective = _dict['date_performed'])
                ucl.save()
                _dict['ucl'] = ucl
                
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
            logger.info('created %d checklists & user_agreements', i)
            logger.info('key: %r', key)
            logger.info('created log: %r, %r', log, log.diffs )
    
    logger.info('created %d user_checklist and user_agreement instances', i)

#     set_user_data_sharing_level_attached_file(apps, schema_editor)

def set_user_data_sharing_level_and_attached_file(apps, schema_editor):
    '''
    Purpose:
    - set the data sharing level for the extant user agreements
    -- NOTE: will not set dsl if no UA found
    - set the attached file for extant user agreements
    MUST be run after create_user_checklist_from_checklist_item_events


/**
  To verify with current roles:

select * from (
with roles as (
  select screensaver_user_id, first_name, last_name,
  (select screensaver_user_role from screensaver_user_role sru
    where screensaver_user_role ~ 'smDsl'
    and sru.screensaver_user_id = su.screensaver_user_id
    order by screensaver_user_role asc
    limit 1 ) user_role,
  (exists(select null from screensaver_user_role sru 
    where sru.screensaver_user_id = su.screensaver_user_id
    and screensaver_user_role = 'screensaverUser') ) can_login
  from screensaver_user su 
)
select roles.*,
  user_role ~ data_sharing_level::text as match,  
  ua.*
  from roles 
  left outer join user_agreement ua on(roles.screensaver_user_id=ua.screensaver_user_id and ua.type='sm')
  where roles.user_role is not null
) a where match is not true;

select * from (
with roles as (
  select screensaver_user_id, first_name, last_name,
  (select screensaver_user_role from screensaver_user_role sru
    where screensaver_user_role ~ 'rnaiDsl'
    and sru.screensaver_user_id = su.screensaver_user_id
    order by screensaver_user_role asc
    limit 1 ) user_role,
  (exists(select null from screensaver_user_role sru 
    where sru.screensaver_user_id = su.screensaver_user_id
    and screensaver_user_role = 'screensaverUser') ) can_login
  from screensaver_user su 
)
select roles.*,
  user_role ~ data_sharing_level::text as match,  
  ua.*
  from roles 
  left outer join user_agreement ua on(roles.screensaver_user_id=ua.screensaver_user_id and ua.type='rnai')
  where roles.user_role is not null
) a where match is not true;

**/

    '''
    
    UserAgreement = apps.get_model('db','UserAgreement')
    connection = schema_editor.connection
    cursor = connection.cursor()

    sql = '''
    with agreement_query as (
      select
      screensaver_user_id, first_name, last_name,
      (select checklist_item_event_id from (
              select * from checklist_item_event cie 
              join checklist_item using(checklist_item_id) 
              where item_name ='{checklist_item_name}' 
              and cie.screening_room_user_id = su.screensaver_user_id
              order by checklist_item_event_id desc, date_performed desc limit 1
              ) as cie_last
         ) checklist_item_event_id,
      (select screensaver_user_role from screensaver_user_role sr
        where sr.screensaver_user_id = su.screensaver_user_id
        and screensaver_user_role ~ '{data_sharing_role_pattern}'
        order by screensaver_user_role limit 1) dsl_role,
      (select attached_file_id 
        from attached_file af 
        where af.screensaver_user_id = su.screensaver_user_id
        and type='{attached_file_type}'
        order by date_created desc
        limit 1) attached_file_id
      from screensaver_user su
    )
    select
    agreement_query.*,
    extract(epoch from (af.date_created-ua.date_active)) date_difference,
    af.date_created af_date_created,
    af.filename,
    cie.is_expiration,
    cie.date_performed cie_date_performed,
    ua.date_active ua_date_active,
    ua.date_expired ua_date_expired,
    ua.id as user_agreement_id
    from 
    agreement_query
    left outer join attached_file af using(attached_file_id)
    left outer join checklist_item_event cie using(checklist_item_event_id)
    left outer join user_agreement ua on(
        agreement_query.screensaver_user_id=ua.screensaver_user_id 
        and ua.type='{user_agreement_type}' );
'''    
    sql_keys = [
        'screensaver_user_id', 'first_name','last_name',
        'checklist_item_event_id','dsl_role','attached_file_id','date_difference',
        'af_date_created','filename','is_expiration','cie_date_performed',
        'ua_date_active','ua_date_expired','user_agreement_id'
        ]
   
    def update_user_agreement(_dict):
        user_agreement = UserAgreement.objects.get(pk=_dict['user_agreement_id'])
        
        user_role = _dict['dsl_role']
        
        if user_role:
            if user_agreement.date_expired is None:
                if '1' in user_role:
                    user_agreement.data_sharing_level = 1
                elif '2' in user_role:
                    user_agreement.data_sharing_level = 2
                elif '3' in user_role:
                    user_agreement.data_sharing_level = 3
                else:
                    logging.error(
                        'unknown user role: %r, %r'
                        % (user_role, _dict))
                    raise ProgrammingError
            else:
                logger.warn('Not setting data_sharing_role, UA is expired on %r for %r',
                    user_agreement.date_expired, _dict)
        
        if _dict['is_expiration'] is True:
            if _dict['cie_date_performed'] != _dict['ua_date_expired']:
                logger.warn('expiration date error: %r, %r', _dict, user_agreement)
                user_agreement.date_expired = _dict['cie_date_performed']
        
        elif _dict['is_expiration'] is False:
            if _dict['cie_date_performed'] != _dict['ua_date_active']:
                logger.warn('active date error: %r, %r', _dict, user_agreement)
                user_agreement.date_active = _dict['cie_date_performed']
            if user_agreement.date_expired is not None:
                logger.warn('user agreement error, active but has date_expired: %r, %r',
                    _dict, user_agreement.date_expired)
                user_agreement.date_expired = None
            if user_agreement.date_notified is not None:
                logger.warn('user agreement error, active but has date_notified: %r, %r',
                    _dict, user_agreement.date_notified)
                user_agreement.date_notified = None
        else:
            logging.error('is_expiration is not set: %r', _dict)
            raise ProgrammingError
                
        if _dict['attached_file_id'] is not None:
            if _dict['date_difference'] > 0:
                user_agreement.file_id = _dict['attached_file_id']
            else:
                if _dict['date_difference'] is None:
                    logging.error('no date difference %r' % _dict)
                    raise ProgrammingError
                if (_dict['date_difference']/(-60*60*24)) > 365:
                    logger.warn(
                        'File date is more than 1 year older than the '
                        'date_performed(date_active), not using: %r', _dict)
                    if _dict['is_expiration'] is False:
                        logger.warn('Active UA has no file assigned')
                else:
                    logger.info(
                        'setting attached file, even though difference is negative: %r',
                        _dict)
                    user_agreement.file_id = _dict['attached_file_id']
        user_agreement.save()
        
    logger.info('Process Small Molecule agreements....')
    cursor.execute(
        sql.format(
            checklist_item_name='Current Small Molecule User Agreement active',
            data_sharing_role_pattern='smDsl',
            attached_file_type='iccb_l_small_molecule_user_agreement',
            user_agreement_type='sm',
        ))
    i = 0
    for row in cursor:
        _dict = dict(zip(sql_keys,row))
        if _dict['user_agreement_id']:
            update_user_agreement(_dict)
            i += 1
            if i % 100 == 0:
                logger.info('Updated %d sm user_agreement records', i)
    logger.info('Updated %d sm user_agreement records', i)

    logger.info('Process RNAi agreements....')
    cursor.execute(
        sql.format(
            checklist_item_name='Current RNAi User Agreement active',
            data_sharing_role_pattern='rnaiDsl',
            attached_file_type='iccb_l_rnai_user_agreement',
            user_agreement_type='rnai',
        ))
    i = 0
    for row in cursor:
        _dict = dict(zip(sql_keys,row))
        if _dict['user_agreement_id']:
            update_user_agreement(_dict)
            i += 1
            if i % 100 == 0:
                logger.info('Updated %d rnai user_agreement records', i)
        
    logger.info('Updated %d rnai user_agreement records', i)


def create_data_sharing_level_update_logs(apps, schema_editor):
    
    sql = ''' 
    select 
        a.activity_id,
        a.performed_by_id,
        a.date_of_activity,
        a.date_created,
        a.comments,
        admin.username as admin_username,
        su.screensaver_user_id, 
        su.first_name, 
        su.last_name
    from activity a
    join screensaver_user admin on(performed_by_id=admin.screensaver_user_id)
    left outer join screensaver_user_update_activity sua on (activity_id=update_activity_id)
    join screensaver_user su on(sua.screensaver_user_id=su.screensaver_user_id) 
    where a.comments ~* 'screens level'
    order by su.screensaver_user_id, date_of_activity asc;
    '''
    sql_keys = [
        'activity_id', 'performed_by_id', 'date_of_activity','date_created',
        'comments','admin_username', 'screensaver_user_id', 'first_name', 
        'last_name'
    ]
    def create_log(_dict):
        log = ApiLog()
        log.ref_resource_name = 'useragreement'
        log.username = _dict['admin_username']
        log.user_id = _dict['performed_by_id']
        log.api_action = 'PATCH'
        log.json_field = { 
            'migration': 'User Data Sharing Roles',
            'data':  {'screensaver_user_update_activity':_dict['activity_id'] } 
            }
        if log.username is None:
            log.username = 'sde_EDIT'
        return log
    
    connection = schema_editor.connection
    cursor = connection.cursor()

    cursor.execute(sql)

    sm_pattern = re.compile(r"Small Molecule Screens Level (\d)")
    sm_from_pattern = re.compile(r"from 'Small Molecule Screens Level (\d)'")
    sm_to_pattern = re.compile(r"to 'Small Molecule Screens Level (\d)'")
    rnai_pattern = re.compile(r"RNAi Screens Level (\d)")
    rnai_from_pattern = re.compile(r"from 'RNAi Screens Level (\d)'")
    rnai_to_pattern = re.compile(r"to 'RNAi Screens Level (\d)'")
    
    i = 0
    for row in cursor:
        _dict = dict(zip(sql_keys,row))
        log = create_log(_dict)
        comments = _dict['comments']
        log.comment = comments
        if 'added' in comments:
            match = sm_pattern.search(comments)
            if match:
                log.key = '/'.join([str(_dict['screensaver_user_id']), 'sm'])
                log.diffs['data_sharing_level'] = [None,match.group(1)]
                logger.debug('1a: %r, %r', log.key, _dict)

            match = rnai_pattern.search(comments)
            if match:
                if log.key:
                    logger.error('activity matches both: %r, %r',log.key, _dict)
                    raise ProgrammingError
                log.key = '/'.join([str(_dict['screensaver_user_id']), 'rnai'])
                log.diffs['data_sharing_level'] = [None,match.group(1)]
                logger.debug('1b: %r, %r', log.key, _dict)
                
        elif 'updated' in comments:
            
            match2 = sm_to_pattern.search(comments)
            if match2:
                val2 = match2.group(1)
                val1 = None
                match1 = sm_from_pattern.search(comments)
                if match1:
                    val1 = match1.group(1)
                
                log.key = '/'.join([str(_dict['screensaver_user_id']), 'sm'])
                log.diffs['data_sharing_level'] = [val1,val2]
                logger.debug('2a: %r %r', log.key, _dict)
            
            match2 = rnai_to_pattern.search(comments)
            if match2:
                if log.key:
                    logger.error('activity matches both: %r, %r',log.key, _dict)
                    raise ProgrammingError
                val2 = match2.group(1)
                val1 = None
                match1 = rnai_from_pattern.search(comments)
                if match1:
                    val1 = match1.group(1)
                
                log.key = '/'.join([str(_dict['screensaver_user_id']), 'rnai'])
                log.diffs['data_sharing_level'] = [val1,val2]
                logger.debug('2b: %r %r', log.key, _dict)

        elif 'removed' in comments:
            match = sm_pattern.search(comments)
            if match:
                log.key = '/'.join([str(_dict['screensaver_user_id']), 'sm'])
                log.diffs['data_sharing_level'] = [match.group(1),None]
                logger.debug('3a: %r, %r', log.key, _dict)

            match = rnai_pattern.search(comments)
            if match:
                if log.key:
                    logger.error('activity matches both: %r, %r',log.key, _dict)
                    raise ProgrammingError
                log.key = '/'.join([str(_dict['screensaver_user_id']), 'rnai'])
                log.diffs['data_sharing_level'] = [match.group(1),None]
                logger.debug('3b: %r, %r', log.key, _dict)
            
            if 'expiring' in comments:
                log.comment = 'Automated User Agreement expiration'
        else:
            logger.error('unrecognized: %r', _dict)
            raise ProgrammingError
 
        if not log.key:
            logger.error('Actvity not recognized: %r', _dict)
            raise ProgrammingError
        
        log.uri = '/'.join([log.ref_resource_name,log.key])
        
        
        set_log_time(log, get_activity_date_time(
            _dict['date_created'], _dict['date_of_activity']))
        log.save()
        i += 1
        
        if i % 100 == 0:
            logger.info('processed %d data sharing level update logs', i)
            logger.info('latest log: %r', log)
            logger.info('%r', _dict)
    logger.info('processed %d data sharing level update logs: %r', i, log)
    
class Migration(migrations.Migration):

    dependencies = [
        ('db', '0004_postprep'),
        ('auth', '0001_initial')
    ]

    operations = [
        # Note: for users migration;
        # - the fields last_notified_smua_checklist_item_event, (also rnai)
        # will be moved to the userchecklistitem status_notified_date;
        # UPDATE: 20171025 - moved to the user_agreement.date_notified
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
        migrations.RunPython(
            create_user_checklist_from_checklist_item_events),
        migrations.RunPython(set_user_data_sharing_level_and_attached_file),
        migrations.RunPython(create_data_sharing_level_update_logs),
    ]
