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
            ('username', self.gf('django.db.models.fields.CharField')(max_length=35)),
            ('ref_resource_name', self.gf('django.db.models.fields.CharField')(max_length=35)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=128)),
            ('uri', self.gf('django.db.models.fields.TextField')()),
            ('date_time', self.gf('django.db.models.fields.DateTimeField')()),
            ('api_action', self.gf('django.db.models.fields.CharField')(max_length=10)),
            ('added_keys', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('removed_keys', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('diff_keys', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('diffs', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
            ('json_field', self.gf('django.db.models.fields.TextField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'reports', ['ApiLog'])

        # Adding unique constraint on 'ApiLog', fields ['ref_resource_name', 'key', 'date_time']
        db.create_unique(u'reports_apilog', ['ref_resource_name', 'key', 'date_time'])

        # Adding model 'MetaHash'
        db.create_table(u'reports_metahash', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('scope', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('alias', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('json_field_type', self.gf('django.db.models.fields.CharField')(max_length=128, null=True, blank=True)),
            ('json_field', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal(u'reports', ['MetaHash'])

        # Adding unique constraint on 'MetaHash', fields ['scope', 'key']
        db.create_unique(u'reports_metahash', ['scope', 'key'])

        # Adding model 'Vocabularies'
        db.create_table(u'reports_vocabularies', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('scope', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('key', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('alias', self.gf('django.db.models.fields.CharField')(max_length=35, blank=True)),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('title', self.gf('django.db.models.fields.CharField')(max_length=100, blank=True)),
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

        # Adding M2M table for field permissions on 'UserGroup'
        m2m_table_name = db.shorten_name(u'reports_usergroup_permissions')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('usergroup', models.ForeignKey(orm[u'reports.usergroup'], null=False)),
            ('permission', models.ForeignKey(orm[u'reports.permission'], null=False))
        ))
        db.create_unique(m2m_table_name, ['usergroup_id', 'permission_id'])


    def backwards(self, orm):
        # Removing unique constraint on 'Permission', fields ['scope', 'key', 'type']
        db.delete_unique(u'reports_permission', ['scope', 'key', 'type'])

        # Removing unique constraint on 'Vocabularies', fields ['scope', 'key']
        db.delete_unique(u'reports_vocabularies', ['scope', 'key'])

        # Removing unique constraint on 'MetaHash', fields ['scope', 'key']
        db.delete_unique(u'reports_metahash', ['scope', 'key'])

        # Removing unique constraint on 'ApiLog', fields ['ref_resource_name', 'key', 'date_time']
        db.delete_unique(u'reports_apilog', ['ref_resource_name', 'key', 'date_time'])

        # Deleting model 'ApiLog'
        db.delete_table(u'reports_apilog')

        # Deleting model 'MetaHash'
        db.delete_table(u'reports_metahash')

        # Deleting model 'Vocabularies'
        db.delete_table(u'reports_vocabularies')

        # Deleting model 'Permission'
        db.delete_table(u'reports_permission')

        # Deleting model 'UserGroup'
        db.delete_table(u'reports_usergroup')

        # Removing M2M table for field permissions on 'UserGroup'
        db.delete_table(db.shorten_name(u'reports_usergroup_permissions'))


    models = {
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
            'ref_resource_name': ('django.db.models.fields.CharField', [], {'max_length': '35'}),
            'removed_keys': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'uri': ('django.db.models.fields.TextField', [], {}),
            'user_id': ('django.db.models.fields.IntegerField', [], {}),
            'username': ('django.db.models.fields.CharField', [], {'max_length': '35'})
        },
        u'reports.metahash': {
            'Meta': {'unique_together': "(('scope', 'key'),)", 'object_name': 'MetaHash'},
            'alias': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json_field': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'json_field_type': ('django.db.models.fields.CharField', [], {'max_length': '128', 'null': 'True', 'blank': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
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
        u'reports.usergroup': {
            'Meta': {'object_name': 'UserGroup'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['reports.Permission']", 'symmetrical': 'False'})
        },
        u'reports.vocabularies': {
            'Meta': {'unique_together': "(('scope', 'key'),)", 'object_name': 'Vocabularies'},
            'alias': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'json_field': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'key': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'scope': ('django.db.models.fields.CharField', [], {'max_length': '35', 'blank': 'True'}),
            'title': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'})
        }
    }

    complete_apps = ['reports']