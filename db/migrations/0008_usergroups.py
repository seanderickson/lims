# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import csv
# from datetime import timedelta
# import datetime
import json
import logging
import os

from django.contrib.auth.models import User, UserManager
from django.db import migrations, models
import pytz

from db.support.data_converter import default_converter
from lims.base_settings import PROJECT_ROOT

logger = logging.getLogger(__name__)

def create_roles(apps, schema_editor):
    '''
    Create the needed Usergroups and Permissions to implement the legacy
    ScreensaverUserRoles
    '''
    
    AuthUserClass = apps.get_model('auth', 'User')
    UserProfileClass = apps.get_model('reports','UserProfile')
    UserGroupClass = apps.get_model('reports','UserGroup')
    ScreensaverUser = apps.get_model('db', 'ScreensaverUser')
    ScreensaverUserRole = apps.get_model('db', 'ScreensaverUserRole')
    role_group_map = {}

    # These are the roles that will be migrated
    user_group_new_names = {
        'billingAdmin': 'billingAdmin',
        'cherryPickRequestsAdmin': 'cherryPickRequestAdmin',
        'labHeadsAdmin': 'labHeadAdmin',
        'librariesAdmin': 'libraryAdmin',
        'libraryCopiesAdmin': 'libraryCopyAdmin',
        'screenResultsAdmin': 'screenResultAdmin',
        'screensAdmin': 'screenAdmin',
        'screenDslExpirationNotify': 'screenDslExpirationNotify',
        'serviceActivityAdmin': 'serviceActivityAdmin',
        'userAgreementExpirationNotify': 'userAgreementExpirationNotify',
        'userChecklistItemsAdmin': 'userChecklistAdmin',
        'userRolesAdmin': 'userGroupAdmin',
        'usersAdmin': 'screensaverUserAdmin',
        'readEverythingAdmin': 'readEverythingAdmin'
    }
    
    # Create UserGroups
    # - Group permissions are set using 
    # /lims/static/production_data/screensaver_usergroups-prod.csv   
    for ss_role_name in ( ScreensaverUserRole.objects.all()
        .distinct('screensaver_user_role')
        .values_list('screensaver_user_role', flat=True) ):
        if ss_role_name in user_group_new_names:
            role_name = user_group_new_names[ss_role_name]
            group,created = UserGroupClass.objects.get_or_create(
                name=role_name)
            group.description = 'Created in migration'
            role_group_map[ss_role_name] = group
            logger.info('created user group: %r, %s',created, role_name)
    
    roles_assigned = 0
    users_assigned = 0
    active_user_count = 0
    
    for su in ScreensaverUser.objects.all():
        logger.info('processing: %r: %r', su.screensaver_user_id, su.username)
        if su.user_id is None:
            logger.warn('skipping user account: %r, no user profile exists', 
                su.screensaver_user_id)
            continue
        up = su.user
        auth_user = up.user
        
        if su.screensaveruserrole_set.exists(): 
            users_assigned += 1
        
        roles = set(
            su.screensaveruserrole_set.all()
                .values_list('screensaver_user_role', flat=True))
        logger.info('user: %r, roles: %r', su.username, roles)
#         if 'rnaiDsl1MutualScreens' in roles:
#             logger.info('set rnai dsl 1')
#             su.rnai_data_sharing_level = 1
#         elif 'rnaiDsl2MutualPositives' in roles:
#             logger.info('set rnai dsl 2')
#             su.rnai_data_sharing_level = 2
#         elif 'rnaiDsl3SharedScreens' in roles:
#             logger.info('set rnai dsl 3')
#             su.rnai_data_sharing_level = 3
#         if 'smDsl1MutualScreens' in roles:
#             logger.info('set sm dsl 1')
#             su.sm_data_sharing_level = 1
#         elif 'smDsl2MutualPositives' in roles:
#             logger.info('set sm dsl 2')
#             su.sm_data_sharing_level = 2
#         elif 'smDsl3SharedScreens' in roles:
#             logger.info('set sm dsl 3')
#             su.sm_data_sharing_level = 3
        
        if 'readEverythingAdmin' in roles or su.classification == 'staff':
            auth_user.is_staff = True
        
        if 'screensaverUser' in roles:
            logger.info('setting user to active: %r', su.username )
            auth_user.is_active = True
            active_user_count += 1
        
        for ss_role_name in roles:
            if ss_role_name in role_group_map:
                ug = role_group_map[ss_role_name]
                ug.users.add(up)
                ug.save()
                roles_assigned += 1
        
        auth_user.save()
        su.save()
        logger.info('su: %r',su.username)
        
    logger.info('created %d roles for %d users, %d active', 
        roles_assigned, users_assigned, active_user_count)
    

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0007_users'),
        ('auth', '0001_initial')
    ]

    operations = [
#         
#         migrations.AddField(
#             model_name='screensaveruser',
#             name='rnai_data_sharing_level',
#             field=models.IntegerField(null=True),
#         ),
#         migrations.AddField(
#             model_name='screensaveruser',
#             name='sm_data_sharing_level',
#             field=models.IntegerField(null=True),
#         ),
        
        migrations.RunPython(create_roles),
    ]
