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


# IF false, do not migrate user's roles
IS_PRODUCTION_READY = False

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

    user_group_new_names = {
        'cherryPickRequestsAdmin': 'cherryPickRequestAdmin',
        'labHeadsAdmin': 'labHeadAdmin',
        'librariesAdmin': 'libraryAdmin',
        'libraryCopiesAdmin': 'libraryCopyAdmin',
        'screenResultsAdmin': 'screenResultAdmin',
        'screensAdmin': 'screenAdmin',
        'userChecklistItemsAdmin': 'userChecklistItemAdmin',
        'userRolesAdmin': 'userGroupAdmin',
        'usersAdmin': 'userAdmin'
    }
    
    # Create UserGroups
    # - Group permissions are set using 
    # /lims/static/production_data/screensaver_usergroups-prod.csv    
    for ssrole in ( ScreensaverUserRole.objects.all()
        .distinct('screensaver_user_role')
        .values_list('screensaver_user_role', flat=True) ):
        role = user_group_new_names.get(ssrole,ssrole)
        group,created = UserGroupClass.objects.get_or_create(name=role)
        role_group_map[role] = group
        logger.info('created user group: %r, %s',created, role)
    
    roles_assigned = 0
    users_assigned = 0
    for su in ScreensaverUser.objects.all():
        up = su.user
        if su.screensaveruserrole_set.exists(): 
            users_assigned += 1
        
        roles = set(
            su.screensaveruserrole_set.all()
                .values_list('screensaver_user_role', flat=True))
        logger.debug('user: %r, roles: %r', su.username, roles)
        if 'rnaiDsl1MutualScreens' in roles:
            roles = roles - set(('rnaiDsl3SharedScreens','rnaiDsl2MutualPositives'))
        if 'rnaiDsl2MutualPositives' in roles:
            roles = roles - set(('rnaiDsl3SharedScreens',))
        if 'smDsl1MutualScreens' in roles:
            roles = roles - set(('smDsl3SharedScreens','smDsl2MutualPositives'))
        if 'smDsl2MutualPositives' in roles:
            roles = roles - set(('smDsl3SharedScreens',))
            
        for ssrole in roles:
            role = user_group_new_names.get(ssrole,ssrole)
            if role in role_group_map:
                ug = role_group_map[role]
                if IS_PRODUCTION_READY:
                    ug.users.add(up)
                    ug.save()
                else:
                    logger.info(
                        'role: %r, not assigned, setting '
                        'IS_PRODUCTION_READY is False', role)
                roles_assigned += 1
            else:
                logger.error('unknown group: %s',role)
            
    logger.info(str(('created',roles_assigned,'roles for',users_assigned,'users')))
    

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0004_users'),
        ('auth', '0001_initial')
    ]

    operations = [
        migrations.RunPython(create_roles),
    ]
