# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.utils.timezone
import django.db.models.deletion
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='ApiLog',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('user_id', models.IntegerField()),
                ('username', models.CharField(max_length=128)),
                ('ref_resource_name', models.CharField(max_length=128, db_index=True)),
                ('key', models.CharField(max_length=128, db_index=True)),
                ('uri', models.TextField()),
                ('date_time', models.DateTimeField()),
                ('api_action', models.CharField(max_length=10, choices=[('POST', 'POST'), ('PUT', 'PUT'), ('CREATE', 'CREATE'), ('PATCH', 'PATCH'), ('DELETE', 'DELETE')])),
                ('added_keys', models.TextField(null=True, blank=True)),
                ('removed_keys', models.TextField(null=True, blank=True)),
                ('diff_keys', models.TextField(null=True, blank=True)),
                ('diffs', models.TextField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
                ('json_field', models.TextField(null=True)),
                ('parent_log', models.ForeignKey(related_name='child_logs', to='reports.ApiLog', null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Job',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('remote_addr', models.TextField(null=True)),
                ('request_method', models.TextField(null=True)),
                ('path_info', models.TextField(null=True)),
                ('comment', models.TextField(null=True)),
                ('input_filename', models.TextField(null=True)),
                ('date_time_fullfilled', models.DateTimeField(null=True)),
                ('date_time_processing', models.DateTimeField(null=True)),
                ('date_time_requested', models.DateTimeField(default=django.utils.timezone.now)),
                ('response_filename', models.TextField(null=True)),
                ('response_content', models.TextField(null=True)),
                ('response_code', models.IntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='ListLog',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('ref_resource_name', models.CharField(max_length=64)),
                ('key', models.CharField(max_length=128)),
                ('uri', models.TextField()),
                ('apilog', models.ForeignKey(to='reports.ApiLog')),
            ],
        ),
        migrations.CreateModel(
            name='MetaHash',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('scope', models.CharField(max_length=64, blank=True)),
                ('key', models.CharField(max_length=64, blank=True)),
                ('alias', models.CharField(max_length=64, blank=True)),
                ('ordinal', models.IntegerField()),
                ('json_field_type', models.CharField(max_length=128, null=True)),
                ('json_field', models.TextField(null=True)),
                ('linked_field_type', models.CharField(max_length=128, null=True, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='Permission',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('scope', models.CharField(max_length=64, blank=True)),
                ('key', models.CharField(max_length=64, blank=True)),
                ('type', models.CharField(max_length=35)),
            ],
        ),
        migrations.CreateModel(
            name='Record',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('base_value1', models.TextField()),
                ('scope', models.CharField(max_length=64, blank=True)),
            ],
        ),
        migrations.CreateModel(
            name='RecordMultiValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.TextField()),
                ('ordinal', models.IntegerField()),
                ('field_meta', models.ForeignKey(to='reports.MetaHash')),
                ('parent', models.ForeignKey(to='reports.Record')),
            ],
        ),
        migrations.CreateModel(
            name='RecordValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value', models.TextField(null=True)),
                ('field_meta', models.ForeignKey(to='reports.MetaHash')),
                ('parent', models.ForeignKey(to='reports.Record')),
            ],
        ),
        migrations.CreateModel(
            name='RecordValueComplex',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('value1', models.TextField(null=True)),
                ('value2', models.TextField(null=True)),
                ('parent', models.OneToOneField(to='reports.Record', unique=True)),
            ],
        ),
        migrations.CreateModel(
            name='UserGroup',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.TextField(unique=True)),
                ('title',models.TextField(unique=True, null=True, blank=True)),
                ('permissions', models.ManyToManyField(to='reports.Permission')),
                ('super_groups', models.ManyToManyField(related_name='sub_groups', to='reports.UserGroup')),
            ],
        ),
        migrations.CreateModel(
            name='UserProfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('username', models.TextField(unique=True)),
                ('phone', models.TextField(null=True)),
                ('mailing_address', models.TextField(null=True)),
                ('comments', models.TextField(null=True)),
                ('ecommons_id', models.TextField(null=True)),
                ('harvard_id', models.TextField(null=True)),
                ('harvard_id_expiration_date', models.DateField(null=True)),
                ('harvard_id_requested_expiration_date', models.DateField(null=True)),
                ('created_by_username', models.TextField(null=True)),
                ('gender', models.CharField(max_length=15, null=True)),
                ('email', models.TextField(null=True)),
                ('json_field_type', models.CharField(max_length=128, null=True)),
                ('json_field', models.TextField(null=True)),
                ('permissions', models.ManyToManyField(to='reports.Permission')),
                ('user', models.OneToOneField(null=True, on_delete=django.db.models.deletion.SET_NULL, to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='Vocabularies',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('scope', models.CharField(max_length=128, blank=True)),
                ('key', models.CharField(max_length=128, blank=True)),
                ('alias', models.CharField(max_length=64, blank=True)),
                ('ordinal', models.IntegerField()),
                ('title', models.CharField(max_length=512, blank=True)),
                ('json_field', models.TextField(null=True)),
            ],
        ),
        migrations.AlterUniqueTogether(
            name='vocabularies',
            unique_together=set([('scope', 'key')]),
        ),
        migrations.AddField(
            model_name='usergroup',
            name='users',
            field=models.ManyToManyField(to='reports.UserProfile'),
        ),
        migrations.AlterUniqueTogether(
            name='permission',
            unique_together=set([('scope', 'key', 'type')]),
        ),
        migrations.AlterUniqueTogether(
            name='metahash',
            unique_together=set([('scope', 'key')]),
        ),
        migrations.AlterUniqueTogether(
            name='recordmultivalue',
            unique_together=set([('field_meta', 'parent', 'ordinal')]),
        ),
        migrations.AlterUniqueTogether(
            name='listlog',
            unique_together=set([('apilog', 'ref_resource_name', 'key', 'uri')]),
        ),
        migrations.AlterUniqueTogether(
            name='apilog',
            unique_together=set([('ref_resource_name', 'key', 'date_time')]),
        ),
    ]
