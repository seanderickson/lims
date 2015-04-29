# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'ApiLog'
        db.create_table(u'reports_apilog', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user_id', self.gf('django.db.models.fields.IntegerField')()),
            ('username', self.gf('django.db.models.fields.CharField')(max_length=128)),
            ('ref_resource_name', self.gf('django.db.models.fields.CharField')(max_length=128)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=128)),
            ('uri', self.gf('django.db.models.fields.TextField')()),
            ('date_time', self.gf('django.db.models.fields.DateTimeField')()),
            ('api_action', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('added_keys', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('removed_keys', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('diff_keys', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('diffs', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('parent_log', self.gf('django.db.models.fields.related.ForeignKey')(related_name='child_logs', null=True, to=orm['reports.ApiLog'])),
            ('json_field', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'reports', ['ApiLog'])

        # Adding unique constraint on 'ApiLog', fields ['ref_resource_name', 'key', 'date_time']
        db.create_unique(u'reports_apilog', ['ref_resource_name', 'key', 'date_time'])

        # Adding model 'ListLog'
        db.create_table(u'reports_listlog', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('apilog', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['reports.ApiLog'])),
            ('ref_resource_name', self.gf('django.db.models.fields.CharField')(max_length=35)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=128)),
            ('uri', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'reports', ['ListLog'])

        # Adding unique constraint on 'ListLog', fields ['apilog', 'ref_resource_name', 'key', 'uri']
        db.create_unique(u'reports_listlog', ['apilog_id', 'ref_resource_name', 'key', 'uri'])

        # Adding model 'MetaHash'
        db.create_table(u'reports_metahash', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('scope', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('alias', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('json_field_type', self.gf('django.db.models.fields.CharField')(max_length=128, null=True, blank=True)),
            ('json_field', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('linked_field_type', self.gf('django.db.models.fields.CharField')(max_length=128, null=True, blank=True)),
        ))
        db.send_create_signal(u'reports', ['MetaHash'])

        # Adding unique constraint on 'MetaHash', fields ['scope', 'key']
        db.create_unique(u'reports_metahash', ['scope', 'key'])

        # Adding model 'Vocabularies'
        db.create_table(u'reports_vocabularies', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('scope', self.gf('django.db.models.fields.CharField')(max_length=128, blank=True)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=128, blank=True)),
            ('alias', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('title', self.gf('django.db.models.fields.CharField')(max_length=512, blank=True)),
            ('json_field', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal(u'reports', ['Vocabularies'])

        # Adding unique constraint on 'Vocabularies', fields ['scope', 'key']
        db.create_unique(u'reports_vocabularies', ['scope', 'key'])

        # Adding model 'Permission'
        db.create_table(u'reports_permission', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('scope', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=15)),
        ))
        db.send_create_signal(u'reports', ['Permission'])

        # Adding unique constraint on 'Permission', fields ['scope', 'key', 'type']
        db.create_unique(u'reports_permission', ['scope', 'key', 'type'])

        # Adding model 'UserGroup'
        db.create_table(u'reports_usergroup', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.TextField')(unique=True)),
        ))
        db.send_create_signal(u'reports', ['UserGroup'])

        # Adding M2M table for field users on 'UserGroup'
        m2m_table_name = db.shorten_name(u'reports_usergroup_users')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('usergroup', models.ForeignKey(orm[u'reports.usergroup'], null=False)),
            ('userprofile', models.ForeignKey(orm[u'reports.userprofile'], null=False))
        ))
        db.create_unique(m2m_table_name, ['usergroup_id', 'userprofile_id'])

        # Adding M2M table for field permissions on 'UserGroup'
        m2m_table_name = db.shorten_name(u'reports_usergroup_permissions')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('usergroup', models.ForeignKey(orm[u'reports.usergroup'], null=False)),
            ('permission', models.ForeignKey(orm[u'reports.permission'], null=False))
        ))
        db.create_unique(m2m_table_name, ['usergroup_id', 'permission_id'])

        # Adding M2M table for field super_groups on 'UserGroup'
        m2m_table_name = db.shorten_name(u'reports_usergroup_super_groups')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('from_usergroup', models.ForeignKey(orm[u'reports.usergroup'], null=False)),
            ('to_usergroup', models.ForeignKey(orm[u'reports.usergroup'], null=False))
        ))
        db.create_unique(m2m_table_name, ['from_usergroup_id', 'to_usergroup_id'])

        # Adding model 'UserProfile'
        db.create_table(u'reports_userprofile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['auth.User'], unique=True)),
            ('username', self.gf('django.db.models.fields.TextField')()),
            ('gender', self.gf('django.db.models.fields.CharField')(max_length=15, null=True)),
            ('phone', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('mailing_address', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('ecommons_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('harvard_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('harvard_id_expiration_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('harvard_id_requested_expiration_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('json_field_type', self.gf('django.db.models.fields.CharField')(max_length=128, null=True, blank=True)),
            ('json_field', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal(u'reports', ['UserProfile'])

        # Adding M2M table for field permissions on 'UserProfile'
        m2m_table_name = db.shorten_name(u'reports_userprofile_permissions')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('userprofile', models.ForeignKey(orm[u'reports.userprofile'], null=False)),
            ('permission', models.ForeignKey(orm[u'reports.permission'], null=False))
        ))
        db.create_unique(m2m_table_name, ['userprofile_id', 'permission_id'])

        # Adding model 'Record'
        db.create_table(u'reports_record', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('base_value1', self.gf('django.db.models.fields.TextField')()),
            ('scope', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
        ))
        db.send_create_signal(u'reports', ['Record'])

        # Adding model 'RecordValue'
        db.create_table(u'reports_recordvalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['reports.Record'])),
            ('field_meta', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['reports.MetaHash'])),
            ('value', self.gf('django.db.models.fields.TextField')(null=True)),
        ))
        db.send_create_signal(u'reports', ['RecordValue'])

        # Adding model 'RecordMultiValue'
        db.create_table(u'reports_recordmultivalue', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['reports.Record'])),
            ('field_meta', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['reports.MetaHash'])),
            ('value', self.gf('django.db.models.fields.TextField')()),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'reports', ['RecordMultiValue'])

        # Adding unique constraint on 'RecordMultiValue', fields ['field_meta', 'parent', 'ordinal']
        db.create_unique(u'reports_recordmultivalue', ['field_meta_id', 'parent_id', 'ordinal'])

        # Adding model 'RecordValueComplex'
        db.create_table(u'reports_recordvaluecomplex', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['reports.Record'], unique=True)),
            ('value1', self.gf('django.db.models.fields.TextField')(null=True)),
            ('value2', self.gf('django.db.models.fields.TextField')(null=True)),
        ))
        db.send_create_signal(u'reports', ['RecordValueComplex'])

        # Adding model 'Job'
        db.create_table(u'reports_job', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('remote_addr', self.gf('django.db.models.fields.TextField')(null=True)),
            ('request_method', self.gf('django.db.models.fields.TextField')(null=True)),
            ('path_info', self.gf('django.db.models.fields.TextField')(null=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True)),
            ('input_filename', self.gf('django.db.models.fields.TextField')(null=True)),
            ('date_time_fullfilled', self.gf('django.db.models.fields.DateTimeField')(null=True)),
            ('date_time_processing', self.gf('django.db.models.fields.DateTimeField')(null=True)),
            ('date_time_requested', self.gf('django.db.models.fields.DateTimeField')(default=datetime.datetime.now)),
            ('response_filename', self.gf('django.db.models.fields.TextField')(null=True)),
            ('response_content', self.gf('django.db.models.fields.TextField')(null=True)),
            ('response_code', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'reports', ['Job'])


    def backwards(self, orm):
        # Removing unique constraint on 'RecordMultiValue', fields ['field_meta', 'parent', 'ordinal']
        db.delete_unique(u'reports_recordmultivalue', ['field_meta_id', 'parent_id', 'ordinal'])

        # Removing unique constraint on 'Permission', fields ['scope', 'key', 'type']
        db.delete_unique(u'reports_permission', ['scope', 'key', 'type'])

        # Removing unique constraint on 'Vocabularies', fields ['scope', 'key']
        db.delete_unique(u'reports_vocabularies', ['scope', 'key'])

        # Removing unique constraint on 'MetaHash', fields ['scope', 'key']
        db.delete_unique(u'reports_metahash', ['scope', 'key'])

        # Removing unique constraint on 'ListLog', fields ['apilog', 'ref_resource_name', 'key', 'uri']
        db.delete_unique(u'reports_listlog', ['apilog_id', 'ref_resource_name', 'key', 'uri'])

        # Removing unique constraint on 'ApiLog', fields ['ref_resource_name', 'key', 'date_time']
        db.delete_unique(u'reports_apilog', ['ref_resource_name', 'key', 'date_time'])

        # Deleting model 'ApiLog'
        db.delete_table(u'reports_apilog')

        # Deleting model 'ListLog'
        db.delete_table(u'reports_listlog')

        # Deleting model 'MetaHash'
        db.delete_table(u'reports_metahash')

        # Deleting model 'Vocabularies'
        db.delete_table(u'reports_vocabularies')

        # Deleting model 'Permission'
        db.delete_table(u'reports_permission')

        # Deleting model 'UserGroup'
        db.delete_table(u'reports_usergroup')

        # Removing M2M table for field users on 'UserGroup'
        db.delete_table(db.shorten_name(u'reports_usergroup_users'))

        # Removing M2M table for field permissions on 'UserGroup'
        db.delete_table(db.shorten_name(u'reports_usergroup_permissions'))

        # Removing M2M table for field super_groups on 'UserGroup'
        db.delete_table(db.shorten_name(u'reports_usergroup_super_groups'))

        # Deleting model 'UserProfile'
        db.delete_table(u'reports_userprofile')

        # Removing M2M table for field permissions on 'UserProfile'
        db.delete_table(db.shorten_name(u'reports_userprofile_permissions'))

        # Deleting model 'Record'
        db.delete_table(u'reports_record')

        # Deleting model 'RecordValue'
        db.delete_table(u'reports_recordvalue')

        # Deleting model 'RecordMultiValue'
        db.delete_table(u'reports_recordmultivalue')

        # Deleting model 'RecordValueComplex'
        db.delete_table(u'reports_recordvaluecomplex')

        # Deleting model 'Job'
        db.delete_table(u'reports_job')


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

    complete_apps = ['reports']