# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging
from django.db import migrations, models
from reports.models import ApiLog

logger = logging.getLogger(__name__)

def create_library_comments(apps,schema_editor):
    
    sql_keys = [
        'activity_id', 'library_id', 'short_name', 'date_created', 'comments', 
        'username', 'performed_by_id']
    sql = '''
        select 
        a.activity_id,
        l.library_id, 
        l.short_name, 
        a.date_created,
        a.comments,
        su.username,
        a.performed_by_id
        from activity a
        join screensaver_user su on(performed_by_id=screensaver_user_id) 
        join library_update_activity lua on(activity_id=lua.update_activity_id) 
        join administrative_activity aa using(activity_id) 
        join library l using(library_id) 
        where administrative_activity_type='Comment'
        order by l.library_id asc, a.date_created asc;    
    '''
    connection = schema_editor.connection
    cursor = connection.cursor()
    try:
        cursor.execute(sql)
        i = 0
        for row in cursor:
            _dict = dict(zip(sql_keys,row))
            
#             logger.info('processing: %r',_dict)
            
            log = ApiLog()
            # Note: as long as 0004_users migration has been completed, all
            # user accounts will have a "username"
            log.username = _dict['username'] 
            log.user_id = _dict['performed_by_id'] 
            log.date_time = _dict['date_created']
            log.api_action = 'PATCH'
            log.ref_resource_name = 'library'
            log.key = _dict['short_name']
            log.uri = '/'.join([log.ref_resource_name,log.key])
            log.comment = _dict['comments']
            log.save()
            
            i += 1

    except Exception, e:
        logger.exception('migration exc')
        raise e  

    print 'created %d library comments' % i


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0009_convert_studies_to_screenresult'),
        ('reports', '0001_initial')
    ]

    operations = [
        migrations.RunPython(create_library_comments),
    ]
