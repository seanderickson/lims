# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.conf import settings
from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='ApiLog',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('user_id', models.IntegerField()),
                ('username', models.CharField(max_length=128)),
                ('ref_resource_name', models.CharField(max_length=128, 
                    db_index=True)),
                ('key', models.CharField(max_length=128, db_index=True)),
                ('uri', models.TextField()),
                ('date_time', models.DateTimeField()),
                ('api_action', models.CharField(max_length=10, 
                    choices=[('POST', 'POST'), ('PUT', 'PUT'), 
                        ('CREATE', 'CREATE'), ('PATCH', 'PATCH'), 
                        ('DELETE', 'DELETE')])),
                ('comment', models.TextField(null=True)),
                ('json_field', models.TextField(null=True)),
                ('parent_log', models.ForeignKey(
                    related_name='child_logs', to='reports.ApiLog', null=True)),
            ],
        ),
        migrations.CreateModel(
            name='LogDiff',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('field_key', models.TextField()),
                ('field_scope', models.TextField()),
                ('before', models.TextField(null=True)),
                ('after', models.TextField(null=True)),
                ('log', models.ForeignKey(to='reports.ApiLog')),
            ],
        ),
        migrations.CreateModel(
            name='MetaHash',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('scope', models.CharField(max_length=64)),
                ('key', models.CharField(max_length=64)),
                ('alias', models.CharField(max_length=64)),
                ('ordinal', models.IntegerField()),
                ('json_field_type', 
                    models.CharField(max_length=128, null=True)),
                ('json_field', models.TextField(null=True)),
                ('linked_field_type', 
                    models.CharField(max_length=128, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Permission',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('scope', models.CharField(max_length=64)),
                ('key', models.CharField(max_length=64)),
                ('type', models.CharField(max_length=35)),
            ],
        ),
        migrations.CreateModel(
            name='UserGroup',
            fields=[
                ('id', models.AutoField(verbose_name='ID', 
                    serialize=False, auto_created=True, primary_key=True)),
                ('name', models.TextField(unique=True)),
                ('title', models.TextField(unique=True, null=True)),
                ('permissions', 
                    models.ManyToManyField(to='reports.Permission')),
                ('super_groups', models.ManyToManyField(
                    related_name='sub_groups', to='reports.UserGroup')),
                ('description', models.TextField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='UserProfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', 
                    serialize=False, auto_created=True, primary_key=True)),
                ('username', models.TextField(unique=True)),
                ('phone', models.TextField(null=True)),
                ('mailing_address', models.TextField(null=True)),
                ('comments', models.TextField(null=True)),
                ('ecommons_id', models.TextField(null=True)),
                ('harvard_id', models.TextField(null=True)),
                ('harvard_id_expiration_date', models.DateField(null=True)),
                ('harvard_id_requested_expiration_date', 
                    models.DateField(null=True)),
                ('created_by_username', models.TextField(null=True)),
                ('gender', models.CharField(max_length=15, null=True)),
                ('permissions', 
                    models.ManyToManyField(to='reports.Permission')),
                ('user', models.OneToOneField(
                    null=True, to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='Vocabulary',
            fields=[
                ('id', models.AutoField(verbose_name='ID', 
                    serialize=False, auto_created=True, primary_key=True)),
                ('scope', models.CharField(max_length=128)),
                ('key', models.CharField(max_length=128)),
                ('alias', models.CharField(max_length=64)),
                ('ordinal', models.IntegerField()),
                ('title', models.CharField(max_length=512)),
                ('is_retired', models.NullBooleanField()),
                ('comment', models.TextField(null=True)),
                ('description', models.TextField(null=True)),
                ('expire_interval_days', models.IntegerField(null=True)),
                ('expire_notifiy_days', models.IntegerField(null=True)),
                ('json_field', models.TextField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Job',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('uri', models.TextField()),
                ('method', models.TextField()),
                ('encoding', models.TextField()),
                ('content_type', models.TextField()),
                ('http_accept', models.TextField()),
                ('params', models.TextField()),
                ('comment', models.TextField(null=True)),
                ('context_data', models.TextField(null=True)),
                ('process_id', models.TextField(null=True)),
                ('state', models.TextField(default=b'pending', 
                    choices=[
                        (b'COMPLETED', b'completed'), 
                        (b'FAILED', b'failed'), 
                        (b'PENDING', b'pending'),
                        (b'PROCESSING', b'processing'),
                        (b'SUBMITTED', b'submitted'),
                        ])),
                ('date_time_requested', 
                    models.DateTimeField(default=django.utils.timezone.now)),
                ('date_time_submitted', models.DateTimeField(null=True)),
                ('date_time_processing', models.DateTimeField(null=True)),
                ('date_time_completed', models.DateTimeField(null=True)),
                ('response_content', models.TextField(null=True)),
                ('response_status_code', models.IntegerField(null=True)),
                ('user_requesting', models.ForeignKey(to='reports.UserProfile')),
            ],
        ),
        
        
        migrations.AlterUniqueTogether(
            name='vocabulary',
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
            name='logdiff',
            unique_together=set([('log', 'field_key', 'field_scope')]),
        ),
        migrations.AlterUniqueTogether(
            name='apilog',
            unique_together=set([('ref_resource_name', 'key', 'date_time')]),
        ),
    ]
