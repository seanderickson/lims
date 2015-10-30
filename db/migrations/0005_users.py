# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
import datetime
import json
import logging
import os

from django.contrib.auth.models import User, UserManager
from django.db import migrations, models
from django.utils.timezone import utc, make_aware
import pytz

from db.support.data_converter import default_converter
from lims.base_settings import PROJECT_ROOT
from reports.models import Vocabularies, ApiLog


# from reports.utils.sqlalchemy_bridge import Bridge
logger = logging.getLogger(__name__)

def create_screensaver_users(apps, schema_editor):
    '''
    Create entries in auth_user, reports.UserProfile for each valid entry 
    in the screensaver_user table.
    - 'valid' entries have [ecommons_id or login_id] and email addresses.
    - Duplicated ecommons/and/or/login ID's are invalid; the second
    (duplicate) entry will be ignored.
    '''
    
    i=0
    skip_count = 0
    
    AuthUserClass = apps.get_model('auth','User')
    UserProfileClass = apps.get_model('reports','UserProfile')
    ScreensaverUser = apps.get_model('db', 'ScreensaverUser')
    auth_user_username_limit = 30
    
    # delete this inconsistent user:
    #  866 |      14 | 2007-05-24 00:00:00-04 | Ernebjerg  | Morten    | morten_ernebjerg@hms.harvard.edu | 617-432-6392 |                 | using wellmate                           |          |                   | me44        | 70572885   |                            |                                      |               |             |                         |        |         | me44
    ssuid = 866
    try:
        obj = ScreensaverUser.objects.get(screensaver_user_id=ssuid)
        obj.delete()
    except Exception,e:
        logger.error(str(('cannot find/delete screensaver_user_id', ssuid, e)))
    # remove the second jgq10 erroneous account
    ssuid = 3166
    try:
        su = ScreensaverUser.objects.get(screensaver_user_id=ssuid)
        username = '%s_%s' % (su.first_name, su.last_name)
        username = default_converter(username)[:auth_user_username_limit]
        su.username = username
        su.save()
    except Exception,e:
        logger.error(str(('cannot find/delete screensaver_user_id', ssuid, e)))

    ssuid = 830
    # for ruchir shahs old acct
    try:
        su = ScreensaverUser.objects.get(screensaver_user_id=ssuid)
        username = '%s_%s' % (su.first_name, su.last_name)
        username = default_converter(username)[:auth_user_username_limit]
        su.username = username
        su.save()
    except Exception,e:
        logger.error(str(('cannot find/delete screensaver_user_id', ssuid, e)))
    ssuid = 3945
    # for min-joon han dupl
    try:
        su = ScreensaverUser.objects.get(screensaver_user_id=ssuid)
        username = '%s_%s' % (su.first_name, su.last_name)
        username = default_converter(username)[:auth_user_username_limit]
        su.username = username
        su.save()
    except Exception,e:
        logger.error(str(('cannot find/delete screensaver_user_id', ssuid, e)))
    ssuid = 129
    # for maria chmura
    try:
        su = ScreensaverUser.objects.get(screensaver_user_id=ssuid)
        username = '%s_%s' % (su.first_name, su.last_name)
        username = default_converter(username)[:auth_user_username_limit]
        su.username = username
        su.save()
    except Exception,e:
        logger.error(str(('cannot find/delete screensaver_user_id', ssuid, e)))
    
    # sean johnston second account
    ssuid = 3758 
    try:
        su = ScreensaverUser.objects.get(screensaver_user_id=ssuid)
        username = '%s_%s_2' % (su.first_name, su.last_name)
        username = default_converter(username)[:auth_user_username_limit]
        su.username = username
        su.save()
    except Exception,e:
        logger.error(str(('cannot find screensaver_user_id', ssuid, e)))
    
    # to be deleted, duplicate account for Zecai      | Liang     | zl59 
    ssuid = 4505
    try:
        su = ScreensaverUser.objects.get(screensaver_user_id=ssuid)
        username = '%s_%s_to_be_deleted' % (su.first_name, su.last_name)
        username = default_converter(username)[:auth_user_username_limit]
        su.username = username
        su.save()
    except Exception,e:
        logger.error(str(('cannot find screensaver_user_id', ssuid, e)))
    
    for su in ScreensaverUser.objects.all():
        logger.debug("processing ecommons: %r, login_id: %r, email: %r, %s,%s"
            % (su.ecommons_id, su.login_id, su.email, su.first_name, su.last_name))
        au = None
        up = None
        if not su.username:
            username = None
            if su.ecommons_id: 
                username = default_converter(str(su.ecommons_id)) # convert in case it has an error
            elif su.login_id:
                username = default_converter(str(su.login_id))
            elif su.first_name is not None and su.last_name is not None:
                username = '%s_%s' % (su.first_name, su.last_name)
                username = default_converter(username)[:auth_user_username_limit]
            elif su.email:
                username = default_converter(su.email)[:auth_user_username_limit]
            else:
                msg = (str((
                     'Cannot create a login account, does not have id information', 
                    su.screensaver_user_id, ',e', su.ecommons_id, ',l', 
                    su.login_id, ',', su.email, su.first_name, su.last_name)) )
                raise Exception(msg)
        
            su.username = username

        username = su.username
        # find or create the auth_user
        try:
            au = AuthUserClass.objects.get(username=username)
            logger.info(str(('found auth_user', username)))
        except Exception, e:
            pass;
            
        if not au:
            try:
                au = AuthUserClass(
                    username = username, 
                    email = su.email if su.email else 'none', 
                    first_name = su.first_name, 
                    last_name = su.last_name,
                    date_joined = datetime.datetime.utcnow().replace(tzinfo=utc),
                    is_active=False,
                    is_staff=False)

                au.save()
                logger.debug('created user from ecommons: %s, %r' 
                    % (su.ecommons_id, au))
            except Exception, e:
                logger.error(str(('cannot create user ',username,e)))
                raise
#                     continue;
        # find or create the userprofile
        try:
            up = UserProfileClass.objects.get(username=username)
            logger.info(str(('found userprofile', username)))
        except Exception, e:
            logger.info(str(('no userprofile', username)))
#                 created_by_username = None
#                 if su.created_by:
#                     up
#                     if su.created_by.ecommons_id:
#                         created_by_username = su.created_by.ecommons_id
#                         if not created_by_username:
#                             created_by_username = su.created_by.login_id
#                             if not created_by_username:
#                                 created_by_username = su.created_by_id
            if su.created_by: 
                created_by_username = su.created_by.username
            else:
                created_by_username = 'sde_EDIT'
                
            up = UserProfileClass()
            
            up.username = username 
            up.email = su.email
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
                
#                 if orm.ScreensaverUser.objects.all().filter(
#                         user=up).exists():
#                     msg = ('==== error: duplicate user found: ',
#                         username, su,
#                         orm.ScreensaverUser.objects.all().filter(user=up))
#                     print msg
#                     logger.error(str((msg, 'skipping')))
#                     continue
#                 
        su.user = up
        su.save()
        logger.info('saved %r, %s' % (up,up.username))
        i += 1
        
    logger.info(str(( 'Converted ', i , ' users, skipped: ', skip_count)))
    
# FIXME: 20150722 - NOT FINISHED ROLES
def create_roles(apps, schema_editor):
    '''
    Create the needed Usergroups and Permissions to implement the legacy
    ScreensaverUserRoles
    '''
    
    AuthUserClass = apps.get_model('auth', 'User')
    UserProfileClass = apps.get_model('reports','UserProfile')
    UserGroupClass = apps.get_model('reports','UserGroup')
    ScreensaverUser = apps.get_model('db', 'ScreensaverUser')
    role_group_map = {}
    
    role_group_map['screensaverUser'] = \
        UserGroupClass.objects.get_or_create(name='screensaverUser')[0]
    role_group_map['smDsl1MutualScreens'] = \
        UserGroupClass.objects.get_or_create(name='smallMoleculeLevel1')[0]
    role_group_map['smDsl2MutualPositives'] = \
        UserGroupClass.objects.get_or_create(name='smallMoleculeLevel2')[0]
    role_group_map['smDsl3SharedScreens'] = \
        UserGroupClass.objects.get_or_create(name='smallMoleculeLevel3')[0]
    role_group_map['rnaiDsl1MutualScreens'] = \
        UserGroupClass.objects.get_or_create(name='rnaiDsl1MutualScreens')[0]
    role_group_map['rnaiDsl2MutualPositives'] = \
        UserGroupClass.objects.get_or_create(name='rnaiDsl2MutualPositives')[0]
    role_group_map['rnaiDsl3SharedScreens'] = \
        UserGroupClass.objects.get_or_create(name='rnaiDsl3SharedScreens')[0]
#      screensaver_user_role     
# -------------------------------
#  billingAdmin
#  cherryPickRequestsAdmin
#  developer
#  labHeadsAdmin
#  librariesAdmin
#  libraryCopiesAdmin
#  marcusAdmin
#  readEverythingAdmin
#  rnaiDsl1MutualScreens
#  rnaiDsl2MutualPositives
#  rnaiDsl3SharedScreens
#  screenDataSharingLevelsAdmin
#  screenDslExpirationNotify
#  screenResultsAdmin
#  screensAdmin
#  screensaverUser
#  serviceActivityAdmin
#  smDsl1MutualScreens
#  smDsl2MutualPositives
#  smDsl3SharedScreens
#  userAgreementExpirationNotify
#  userChecklistItemsAdmin
#  userRolesAdmin
#  usersAdmin

    i = j = 0
    roles_assigned = 0
    for up in UserProfileClass.objects.all():
        su_query = ScreensaverUser.objects.all().filter(user=up)
        if su_query.exists():
            su = su_query[0]
            
            for role in su.screensaveruserrole_set.all():
                 logger.debug(str(( 'user', up.username , 'found role', role.screensaver_user_role)))
                 if j==0: i += 1
                 j += 1
                 if role.screensaver_user_role in role_group_map:
                     ug = role_group_map[role.screensaver_user_role]
                     ug.users.add(up)
                     ug.save()
                 else:
                    logger.error(str(('unknown group',role.screensaver_user_role)) )
            roles_assigned += j
            j = 0
    logger.info(str(('created',roles_assigned,'roles for',i,'users')))


def create_user_checklist_items(apps, schema_editor):
    
    # prerequisites: 
    # - convert checklist_item / checklist_item_event entries into into 
    # checklistitem.* vocabularies (migration 0002)
    # - create the user_checklist_item table (0002)

    ChecklistItem = apps.get_model('db','ChecklistItem')
    UserChecklistItem = apps.get_model('db','UserChecklistItem')
    ci_group_map = {}
    for obj in ChecklistItem.objects.all().distinct('checklist_item_group'):
        key = default_converter(obj.checklist_item_group)
        ci_group_map[obj.checklist_item_group] = key 
    
    ci_name_map = {}
    for obj in ChecklistItem.objects.all().distinct('item_name'):
        key = default_converter(obj.item_name)
        ci_name_map[obj.item_name] = key

    # create entries in the user_checklist_item table
    # note: status values are hard-coded to correspond to the vocabulary
    # keys (created in migration 0002)
    sql_keys = [
        'suid','cigroup','ciname',
        'su_username','admin_username','admin_suid','admin_upid',
        'date_performed', 'date_created','status',
        ]
    sql = '''
select
screening_room_user_id,
ci.checklist_item_group,
ci.item_name,
su.username su_username,
admin.username admin_username,
admin.screensaver_user_id admin_suid,
up.id admin_upid,
cie.date_performed,
cie.date_created,
case when cie.is_not_applicable then 'n_a'
 when ci.is_expirable and cie.date_performed is not null then
case when cie.is_expiration then 'deactivated' else 'activated' end
 when cie.date_performed is not null then 'completed'     
 else 'not_completed'
 end as status
from checklist_item ci
join checklist_item_event cie using(checklist_item_id)
join screensaver_user su on screening_room_user_id=su.screensaver_user_id
join screensaver_user admin on cie.created_by_id=admin.screensaver_user_id
left join reports_userprofile up on up.id=admin.user_id
order by screening_room_user_id, checklist_item_group, item_name, cie.date_performed asc;
'''
    connection = schema_editor.connection
    cursor = connection.cursor()

#     bridge = Bridge()
#     conn = bridge.get_engine().connect()

    log_ref_resource_name = 'userchecklistitem'
    
    _dict = None
    log = None
    
    uci_hash = {}
    unique_log_keys = set()
    try:
        cursor.execute(sql)
        i = 0
        for row in cursor:
            _dict = dict(zip(sql_keys,row))
            
            key = '/'.join([str(_dict['suid']),_dict['cigroup'],_dict['ciname']])
            previous_dict = uci_hash.get(key)
            logger.debug('prev_dict: %s:%s' % (key,previous_dict))
            if previous_dict:
                uci = previous_dict['obj']
                uci.admin_user_id = int(_dict['admin_suid'])
                uci.status = _dict['status']
                uci.status_date = _dict['date_performed']
                uci.save()
                logger.debug('updated: %r' % uci)
                
            else:
                uci_hash[key] = _dict
                logger.debug(str(('create user checklist item', _dict, 
                    _dict['date_performed'].isoformat())))
                uci = UserChecklistItem.objects.create(
                    screensaver_user_id = int(_dict['suid']),
                    admin_user_id = int(_dict['admin_suid']),
                    item_group = ci_group_map[_dict['cigroup']],
                    item_name = ci_name_map[_dict['ciname']],
                    status = _dict['status'],
                    status_date = _dict['date_performed'])
                uci.save()
                _dict['obj'] = uci
                
                logger.debug('created: %r' % (uci))
                i += 1

            date_time = pytz.utc.localize(_dict['date_created'])                
            if date_time.date() != _dict['date_performed']:
                # only use the less accurate date_performed date if that date
                # is not equal to the date_created date
#                     date_time = _dict['date_performed']
                date_time = datetime.datetime.combine(
                    _dict['date_performed'],
                    datetime.datetime.min.time())
            # create the apilog for this item
            log = ApiLog()
            log.ref_resource_name = log_ref_resource_name
            log.key = '/'.join([_dict['su_username'],uci.item_group,uci.item_name])
            log.username = _dict['admin_username']
            log.user_id = _dict['admin_upid']
            log.date_time = date_time
            log.api_action = 'PATCH'
            log.uri = '/'.join([log.ref_resource_name,log.key])
            log.comment = 'status=%s' % _dict['status']
            
            # is the key (date_time, actually) unique?
            full_key = '/'.join([log.ref_resource_name,log.key,str(log.date_time)])
            while full_key in unique_log_keys:
                # add a second to make it unique; because date performed is a date,
                logger.info(str(('time collision for: ',full_key)))
                log.date_time = log.date_time  + datetime.timedelta(0,1)
                full_key = '/'.join([log.ref_resource_name,log.key,str(log.date_time)])
                
            unique_log_keys.add(full_key)
            if previous_dict:
                diff_keys = ['status']
                diffs = {}
                logger.debug(str(('found previous_dict', previous_dict)))
                diff_keys.append('admin_username')
                diffs['admin_username'] = [previous_dict['admin_username'], _dict['admin_username']]
                
                diff_keys.append('status_date')
                diffs['status_date'] = [
                    previous_dict['date_performed'].isoformat(), 
                    _dict['date_performed'].isoformat()]
                
                diffs['status'] = [previous_dict['status'],_dict['status']]
            
                log.diff_keys = json.dumps(diff_keys)
                log.diffs = json.dumps(diffs)
 
            logger.debug(str(('create log', log)))
            
            log.save()
            log = None
            if i%1000 == 0:
                logger.info(str(('created', i, 'logs')))
    except Exception, e:
        logger.exception('migration exc')
        raise e  

    print 'created %d user_checklist_items' % i


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0004_screen_status'),
        ('auth', '0001_initial')
    ]

    operations = [
        migrations.RunPython(create_screensaver_users),
        migrations.RunPython(create_roles),
        migrations.RunPython(create_user_checklist_items),
    ]
