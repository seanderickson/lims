# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging

from django.db import migrations, models

import db.models
import django.db.models.deletion


#####
# Migration prep:
# Perform migrations required so that Django may use the legacy schema for 
# subsequent actions
# * This migration must be performed before the migration bootstrap step
# * see migration.sh
#####

logger = logging.getLogger(__name__)

def _update_table_autofield(db, table, column):
    
    ###
    # Converting to Autofields in Django
    # Creating primary key "AutoFields" for the old java/hibernate 
    # "GenericGenerator" class
    #
    # example:
    # Changing field 'Screen.screen_id' to auto field
    # *NOTE: the following does not work with Postgres, for an already 
    # existing field:
    # db.alter_column(u'screen', 'screen_id', self.gf(
    #     'django.db.models.fields.AutoField')(primary_key=True))
    # ( Postgres can create a field of type 'serial', but it is not a real type, 
    # so postgres will not convert a field to 'serial'; as would be needed to work
    # with the Django Autofield;
    # see: http://www.postgresql.org/docs/8.3/interactive/datatype-numeric.html#DATATYPE-SERIAL
    # see: http://south.aeracode.org/ticket/407, 
    # Fix is as follows:
    
    # Note: we don't need to create the sequence; just re-associate the old one
    # Note - 20140826:
    # -- moved out of 0012_create_reagent_autofield
    # because the migration bootstrap step needs these fields to be here.
    # Also: note: we have not modified the "models" section below, because some 
    # we would also have to change all the other migrations for consistency;
    # we'll just keep the very last schema migration up-to-date.
    
    # Note: we don't need to create the sequence; just re-associate the old one

    # first change the name to use the standard django naming
    db.execute(
        "ALTER SEQUENCE {column}_seq RENAME to {table}_{column}_seq".format(
            table=table, column=column))
    db.execute(
        ( "ALTER TABLE {table} ALTER COLUMN {column} "
          "SET DEFAULT nextval('{table}_{column}_seq'::regclass)").format(
              table=table, column=column))
    db.execute(
        "ALTER SEQUENCE {table}_{column}_seq OWNED BY {table}.{column}".format(
            table=table, column=column))

def convert_django_autofields(apps, schema_editor):
    
    table_autofields = (
        ('result_value', 'result_value_id'),
        ('reagent', 'reagent_id'),
        ('screen', 'screen_id'),
        ('library', 'library_id'),
        ('gene', 'gene_id'),
        ('screensaver_user', 'screensaver_user_id'),
        ('attached_file', 'attached_file_id'),
        ('activity', 'activity_id'),
        ('equipment_used','equipment_used_id'),
        ('annotation_type','annotation_type_id'),
        ('assay_plate','assay_plate_id'),
        ('assay_well','assay_well_id'),
        ('cherry_pick_request','cherry_pick_request_id'),
        ('lab_cherry_pick','lab_cherry_pick_id'),
        ('cherry_pick_assay_plate','cherry_pick_assay_plate_id'),
        ('lab_affiliation','lab_affiliation_id'),
        ('publication','publication_id'),
        ('screen_result','screen_result_id'),
        ('data_column','data_column_id'),
        ('screener_cherry_pick','screener_cherry_pick_id'),
        ('copy','copy_id'),
        ('plate','plate_id'),
        ('plate_location','plate_location_id'),
        ('well_volume_adjustment','well_volume_adjustment_id'),
        )
    for (table, key_field) in table_autofields:
        logger.info(str(('_update_table_autofield', table, key_field)))
        _update_table_autofield(schema_editor,table, key_field)


def _alter_table_reference(db, sub_table, fk_column, 
                        new_parent, new_parent_column, new_primary_key=None,
                        null_ok=True ):
    '''
    Purpose, flatten the table hierarchy for specified tables.
    '''
    # alter table molfile rename column reagent_id to smr_reagent_id;
    # alter table molfile add column reagent_id integer;
    # update molfile set reagent_id = smr_reagent_id;
    # alter table molfile add constraint reagent_fk FOREIGN KEY (reagent_id) 
    #     REFERENCES reagent (reagent_id);
    # alter table molfile alter column reagent_id set NOT NULL;
    # alter table molfile drop column smr_reagent_id;
    # ALTER TABLE small_molecule_compound_name ADD PRIMARY KEY (reagent_id, ordinal);

    ## NOTE: we are copying/deleting/making new foreign key because it is 
    ## proving difficult to find the constraint to drop for the extant foreign key
    db.execute(
        ( "ALTER TABLE {table} RENAME COLUMN {column} to tmp_{column}").format(
              table=sub_table, column=fk_column))
    db.execute(
        ( "ALTER TABLE {table} ADD COLUMN {column} integer").format(
              table=sub_table, column=fk_column))
    db.execute(
        ( "update {table} set {column} = tmp_{column}").format(
              table=sub_table, column=fk_column))
    db.execute(
        ("ALTER TABLE {table} ADD CONSTRAINT fk_{column} "
            "FOREIGN KEY ({column}) "
            "REFERENCES {other_table} ({other_column}) ").format(
                table=sub_table, column=fk_column, 
                other_table=new_parent, other_column=new_parent_column))
    db.execute(
        ( "ALTER TABLE {table} DROP COLUMN tmp_{column} ").format(
              table=sub_table, column=fk_column))

    if new_primary_key or not null_ok:
        db.execute(
            ( "ALTER TABLE {table} ALTER COLUMN {column} set NOT NULL").format(
                  table=sub_table, column=fk_column))
    if new_primary_key:
        db.execute(
            ( "ALTER TABLE {table} ADD PRIMARY KEY ({primary_key})").format(
                  table=sub_table, primary_key=new_primary_key))
    

def alter_table_parents(apps, schema_editor):
    '''
    Flatten the inheritance hierarchy for reagent tables
    '''
    
    db = schema_editor
    _alter_table_reference(
        db,'molfile','reagent_id', 'reagent', 'reagent_id', 'reagent_id')
    _alter_table_reference(
        db,'small_molecule_compound_name', 
        'reagent_id', 'reagent', 'reagent_id', 'reagent_id,ordinal')
    _alter_table_reference(
        db,'small_molecule_pubchem_cid',
        'reagent_id', 'reagent', 'reagent_id', 'reagent_id,pubchem_cid')
    _alter_table_reference(
        db,'small_molecule_chembank_id',
        'reagent_id', 'reagent', 'reagent_id', 'reagent_id,chembank_id')
    _alter_table_reference(
        db,'small_molecule_chembl_id',
        'reagent_id', 'reagent', 'reagent_id', 'reagent_id,chembl_id')


def alter_table_references(apps, schema_editor):
    '''
    Flatten the inheritance hierarchy for these tables:
    - move ScreeningRoomUser to ScreensaverUser
    - move AdministrativeActivity to Activity
    '''
    db = schema_editor
    
    # move foreign key from screening_room_user to screensaver_user
    _alter_table_reference(
        db,'attached_file', 'screensaver_user_id', 'screensaver_user', 
        'screensaver_user_id', null_ok=True)
    _alter_table_reference(
        db,'library', 'owner_screener_id','screensaver_user', 
        'screensaver_user_id')
    _alter_table_reference(
        db,'service_activity', 'serviced_user_id','screensaver_user', 
        'screensaver_user_id', null_ok=False)
    _alter_table_reference(
        db,'cherry_pick_request', 'requested_by_id','screensaver_user', 
        'screensaver_user_id', null_ok=False)
    _alter_table_reference(
        db,'cherry_pick_request', 'volume_approved_by_id','screensaver_user', 
        'screensaver_user_id', null_ok=True)
    _alter_table_reference(
        db,'screen', 'lead_screener_id','screensaver_user', 
        'screensaver_user_id', null_ok=True)
    _alter_table_reference(
        db,'screen', 'lab_head_id','screensaver_user', 
        'screensaver_user_id', null_ok=True)

    _alter_table_reference(
        db,'screen', 'pin_transfer_admin_activity_id','activity', 
        'activity_id', null_ok=True)


def add_timezone_to_timestamp_fields(apps, schema_editor):
    
    # Change timestamp fields to use timezone information:
    # From the docs:
    # The PostgreSQL backend stores datetimes as timestamp with time zone. 
    # In practice, this means it converts datetimes from the connection’s time zone 
    # to UTC on storage, and from UTC to the connection’s time zone on retrieval.
    # 
    # As a consequence, if you’re using PostgreSQL, you can switch between 
    # USE_TZ = False and USE_TZ = True freely. The database connection’s time zone 
    # will be set to TIME_ZONE or UTC respectively, so that Django obtains correct 
    # datetimes in all cases. You don’t need to perform any data conversions.
    # 
    # Note: it appears that legacy Screensaver dates were stored without timezone
    # information. This modification will apply the database default timezone 
    # (America/New_York for orchestra), which will be correct.
    
    table_columns = (
        ('attached_file', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('activity', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('checklist_item_event', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('cherry_pick_request', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('copy', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('library', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('plate', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('screen', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('screensaver_user', ('date_created', 'date_loaded', 'date_publicly_available')),
        ('screen_result', ('date_created', 'date_loaded', 'date_publicly_available')),
        )
    
    for table, columns in table_columns:
        for column in columns:
            schema_editor.execute(
                ('ALTER TABLE {table} ALTER COLUMN {column} '
                    'SET DATA TYPE timestamp with time zone').format(
                      table=table, column=column))


def create_reagent_ids(apps, schema_editor):
    
    for reagent in apps.get_model('db','Reagent').objects.all():
        reagent.substance_id = db.models.create_id()
        reagent.save()


class Migration(migrations.Migration):

    dependencies = [
        ('db', '0001_initial'),
        ('reports', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Substance',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, 
                    auto_created=True, primary_key=True)),
                ('comment', models.TextField(null=True)),
            ],
            options={
                'db_table': 'substance',
            },
        ),
        
        migrations.AddField(
            model_name='library',
            name='version_number',
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name='reagent',
            name='substance_id',
            field=models.CharField(null=True,max_length=8),
        ),
        migrations.AlterField(
            model_name='silencingreagent',
            name='is_restricted_sequence',
            field=models.NullBooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='smallmoleculereagent',
            name='is_restricted_structure',
            field=models.NullBooleanField(default=False),
        ),
        
        migrations.AddField(
            model_name='screen',
            name='status_date',
            field=models.DateField(null=True),
        ),
        migrations.AddField(
            model_name='screen',
            name='status',
            field=models.TextField(null=True),
        ),
        
        migrations.AddField(
            model_name='screen',
            name='parent_screen',
            field=models.ForeignKey(
                related_name='follow_up_screen', to='db.Screen', null=True),
        ),
        migrations.AlterField(
            model_name='screen',
            name='study_type',
            field=models.TextField(null=True),
        ),
        
        migrations.AddField(
            model_name='silencingreagent',
            name='vendor_gene',
            field=models.OneToOneField(
                related_name='vendor_reagent', null=True, 
                to='db.Gene', unique=True),
        ),
        migrations.AddField(
            model_name='silencingreagent',
            name='facility_gene',
            field=models.OneToOneField(
                related_name='facility_reagent', null=True, 
                to='db.Gene', unique=True),
        ),
        migrations.AddField(
            model_name='library',
            name='loaded_by',
            field=models.ForeignKey(related_name='libraries_loaded', 
                to='db.ScreensaverUser', null=True),
        ),
        migrations.RunPython(convert_django_autofields),
        migrations.RemoveField(model_name='reagent',name='facility_batch_id'),
        migrations.RemoveField(model_name='smallmoleculereagent',name='salt_form_id'),
        migrations.RemoveField(model_name='well',name='version'),
        migrations.RemoveField(model_name='library',name='version'),
        migrations.RemoveField(model_name='screensaveruser',name='version'),
        migrations.RemoveField(model_name='activity',name='version'),
        migrations.RemoveField(model_name='attachedfile',name='version'),
        migrations.RemoveField(model_name='abasetestset',name='version'),
        migrations.RemoveField(model_name='equipmentused',name='version'),
        migrations.RemoveField(model_name='annotationtype',name='version'),
        migrations.RemoveField(model_name='assayplate',name='version'),
        migrations.RemoveField(model_name='assaywell',name='version'),
        migrations.RemoveField(model_name='cellline',name='version'),
        migrations.RemoveField(model_name='checklistitem',name='version'),
        migrations.RemoveField(model_name='cherrypickrequest',name='version'),
        migrations.RemoveField(model_name='labcherrypick',name='version'),
        migrations.RemoveField(model_name='cherrypickassayplate',name='version'),
        migrations.RemoveField(model_name='labaffiliation',name='version'),
        migrations.RemoveField(model_name='publication',name='version'),
        migrations.RemoveField(model_name='screen',name='version'),
        migrations.RemoveField(model_name='screenresult',name='version'),
        migrations.RemoveField(model_name='datacolumn',name='version'),
        migrations.RemoveField(model_name='screenercherrypick',name='version'),
        migrations.RemoveField(model_name='transfectionagent',name='version'),
        migrations.RemoveField(model_name='librarycontentsversion',name='version'),
        migrations.RemoveField(model_name='copy',name='version'),
        migrations.RemoveField(model_name='plate',name='version'),
        migrations.RemoveField(model_name='wellvolumeadjustment',name='version'),

        # screening_room_user decommission
        migrations.AddField(
            model_name='screensaveruser',
            name='classification', 
            field=models.TextField(null=True)),
        migrations.RunSQL(
            'UPDATE screensaver_user su '
            ' set classification=user_classification '
            ' from  screening_room_user sru '
            ' where sru.screensaver_user_id=su.screensaver_user_id'),

        migrations.RunPython(alter_table_references),
        migrations.RunPython(add_timezone_to_timestamp_fields),
        
        migrations.RunPython(alter_table_parents),
        migrations.RunSQL(("ALTER TABLE {table} DROP COLUMN {column} ").format(
                  table='molfile', column='ordinal')),
        migrations.RunSQL((
            "ALTER TABLE {table} ADD CONSTRAINT "
            "    {table}_{column}_unique UNIQUE({column})").format(
            table='screen_result', column='screen_id')),
        migrations.AddField(
            model_name='screensaveruser',
            name='user',
            field=models.OneToOneField(
                on_delete=django.db.models.deletion.SET_NULL, 
                to='reports.UserProfile', null=True),
        ),
        migrations.AddField(
            model_name='screensaveruser',
            name='username', 
            field=models.TextField(unique=True, null=True)),
        migrations.AlterField(
            model_name='screensaveruser',
            name='login_id',
            field=models.TextField(unique=True, null=True)),
        
        # IN PROGRESS - well volume/concentration related fields
        
        migrations.AddField(
            model_name='plate',
            name='remaining_well_volume', 
            field=models.DecimalField(
                null=True, max_digits=10, decimal_places=9)),
        migrations.AddField(
            model_name='plate',
            name='experimental_well_count', 
            field=models.IntegerField(null=True)),
        migrations.AddField(
            model_name='plate',
            name='screening_count', 
            field=models.IntegerField(default=0, null=True)),
        migrations.AddField(
            model_name='plate',
            name='cplt_screening_count', 
            field=models.IntegerField(default=0, null=True)),
        migrations.AddField(
            model_name='plate',
            name='mg_ml_concentration',
            field=models.DecimalField(null=True, max_digits=5, decimal_places=3),
        ),
        migrations.AddField(
            model_name='plate',
            name='molar_concentration',
            field=models.DecimalField(null=True, max_digits=13, decimal_places=12),
        ),
        
        migrations.AddField(
            model_name='plate',
            name='date_plated',
            field=models.DateField(null=True),
        ),
        
        migrations.AddField(
            model_name='plate',
            name='date_retired',
            field=models.DateField(null=True),
        ),
        # 20170407 - added for screening room staff
        migrations.AddField(
            model_name='plate',
            name='is_active',
            field=models.NullBooleanField(),
        ),        
        
        migrations.CreateModel(
            name='CopyWell',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, 
                    auto_created=True, primary_key=True)),
                # ('plate_number', models.IntegerField()),
                ('volume', models.DecimalField(null=True, max_digits=10, decimal_places=9)),
                ('initial_volume', models.DecimalField(null=True, max_digits=10, decimal_places=9)),
                ('cherry_pick_screening_count', models.IntegerField(null=True)),
                ('copy', models.ForeignKey(to='db.Copy')),
                ('mg_ml_concentration',
                    models.DecimalField(null=True, max_digits=5, decimal_places=3)),
                ('molar_concentration',
                    models.DecimalField(null=True, max_digits=13, decimal_places=12)),
            ],
            options={
                'db_table': 'copy_well',
            },
        ),
        migrations.AddField(
            model_name='copywell',
            name='plate',
            field=models.ForeignKey(to='db.Plate'),
        ),
        migrations.AddField(
            model_name='copywell',
            name='well',
            field=models.ForeignKey(to='db.Well'),
        ),
        migrations.AlterUniqueTogether(
            name='copywell',
            unique_together=set([('copy', 'plate', 'well')]),
        ),
        # TODO: moved to 102, not working
        # migrations.AlterUniqueTogether(
        #     name='assayplate',
        #     unique_together=set([('library_screening', 'plate', 'replicate_ordinal')]),
        # ),
        migrations.AlterField(
            model_name='assayplate',
            name='plate_number',
            field=models.IntegerField(db_index=True),
        ),
        migrations.AlterField(
            model_name='assayplate',
            name='replicate_ordinal',
            field=models.IntegerField(db_index=True),
        ),

        #####
        
        migrations.CreateModel(
            name='CachedQuery',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, 
                    auto_created=True, primary_key=True)),
                ('key', models.TextField(unique=True)),
                ('sql', models.TextField()),
                ('uri', models.TextField()),
                ('params', models.TextField(null=True)),
                ('datetime', models.DateTimeField(default=django.utils.timezone.now)),
                ('username', models.CharField(max_length=128)),
                ('count', models.IntegerField(null=True)),
            ],
            options={
                'db_table': 'cached_query',
            },
        ),
        migrations.CreateModel(
            name='WellQueryIndex',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('query', models.ForeignKey(to='db.CachedQuery')),
                ('well', models.ForeignKey(to='db.Well')),
            ],
            options={
                'db_table': 'well_query_index',
            },
        ),
        migrations.CreateModel(
            name='UserChecklist',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('name', models.TextField()),
                ('is_checked', models.BooleanField()),
                ('date_effective', models.DateField()),
                ('date_notified', models.DateField(null=True)),
                ('admin_user', models.ForeignKey(
                    related_name='userchecklistitems_created', 
                    to='db.ScreensaverUser')),
                ('screensaver_user', models.ForeignKey(to='db.ScreensaverUser')),
            ],
            options={
                'db_table': 'user_checklist',
            },
        ),
        migrations.AlterUniqueTogether(
            name='userchecklist',
            unique_together=set([('screensaver_user', 'name')]),
        ),
        migrations.AddField(
            model_name='attachedfile',
            name='type',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='attachedfile',
            name='attached_file_type',
            field=models.ForeignKey(to='db.AttachedFileType', null=True)),
        migrations.AddField(
            model_name='serviceactivity',
            name='funding_support',
            field=models.TextField(null=True),
        ),
        # 20170524
        # Allow serviced user to be null; either serviced_user or 
        # serviced_screen required (TODO: verify workflow with JAS/KR)
        migrations.AlterField(
            model_name='serviceactivity',
            name='serviced_user',
            field=models.ForeignKey(to='db.ScreensaverUser', null=True),
        ),
        
        migrations.CreateModel(
            name='ScreenFundingSupports',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, 
                    auto_created=True, primary_key=True)),
                ('funding_support', models.TextField()),
                ('screen', models.ForeignKey(
                    related_name='fundingsupports', to='db.Screen')),
            ],
            options={
                'db_table': 'screen_funding_supports',
            },
        ),
        migrations.AlterUniqueTogether(
            name='screenfundingsupports',
            unique_together=set([('screen', 'funding_support')]),
        ),
        migrations.CreateModel(
            name='ScreenCellLines',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('cell_line', models.TextField()),
                ('screen', models.ForeignKey(
                    related_name='celllines', to='db.Screen')),
            ],
            options={
                'db_table': 'screen_cell_lines',
            },
        ),
        migrations.AlterUniqueTogether(
            name='screencelllines',
            unique_together=set([('screen', 'cell_line')]),
        ),
        
        # NOTE: temporary for migration
        migrations.AddField(
            model_name='screen',
            name='transfection_agent_text',
            field=models.TextField(null=True),
        ),
        
        migrations.CreateModel(
            name='UserFacilityUsageRole',
            fields=[
                ('id', models.AutoField(
                    verbose_name='ID', serialize=False, auto_created=True, 
                    primary_key=True)),
                ('facility_usage_role', models.TextField()),
                ('screensaver_user', models.ForeignKey(to='db.ScreensaverUser')),
            ],
            options={
                'db_table': 'user_facility_usage_role',
            },
        ),
        migrations.AlterUniqueTogether(
            name='userfacilityusagerole',
            unique_together=set([('screensaver_user', 'facility_usage_role')]),
        ),
        migrations.RunSQL('alter table reagent alter column library_contents_version_id drop not null'),

        #  Update assay_well with the plate_number to expedite plate data loading stats 
        migrations.AddField(
            model_name='assaywell',
            name='plate_number', 
            field=models.IntegerField(null=True)),
        migrations.RunSQL(
            'update assay_well '
            'set plate_number = substring(well_id from 1 for 5 )::integer;'),
        migrations.AlterField(
            model_name='assaywell',
            name='plate_number',
            field=models.IntegerField(null=False)),
        migrations.AlterField(
            model_name='screensaveruser',
            name='ecommons_id',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='screensaveruser',
            name='email',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='screenresult',
            name='experimental_well_count',
            field=models.IntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='screenresult',
            name='replicate_count',
            field=models.IntegerField(default=0),
        ),
        migrations.AlterField(
            model_name='screenresult',
            name='channel_count',
            field=models.IntegerField(default=0, null=True),
        ),
        migrations.AlterField(
            model_name='assaywell',
            name='is_positive',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='datacolumn',
            name='is_follow_up_data',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='datacolumn',
            name='is_derived',
            field=models.BooleanField(default=False),
        ),
        
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='assay_protocol_comments',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='cherry_pick_assay_protocols_followed',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='cherry_pick_followup_results_status',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='comments',
            field=models.TextField(null=True),
        ),
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='is_randomized_assay_plate_layout',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='keep_source_plate_cherry_picks_together',
            field=models.BooleanField(default=True),
        ),
        migrations.AlterField(
            model_name='cherrypickrequest',
            name='number_unfulfilled_lab_cherry_picks',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='cherrypickrequest',
            name='wells_to_leave_empty',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='cherrypickrequest',
            name='date_volume_reserved',
            field=models.DateField(null=True),
        ),
        migrations.AddField(
            model_name='activity',
            name='apilog_uri',
            field=models.TextField(null=True),
        ),
        
        migrations.AddField(
            model_name='screensaveruser',
            name='lab_head_appointment_category', 
            field=models.TextField(null=True)),

        migrations.AddField(
            model_name='screensaveruser',
            name='lab_head_appointment_department', 
            field=models.TextField(null=True)),

        migrations.AddField(
            model_name='screensaveruser',
            name='lab_head_appointment_update_date', 
            field=models.DateField(null=True)),

        migrations.RunSQL(
            'UPDATE screensaver_user su '
            ' set lab_head_appointment_category=lh.lab_head_appointment_category '
            ' from  lab_head lh '
            ' where lh.screensaver_user_id=su.screensaver_user_id'),
        migrations.RunSQL(
            'UPDATE screensaver_user su '
            ' set lab_head_appointment_department=lh.lab_head_appointment_department '
            ' from  lab_head lh '
            ' where lh.screensaver_user_id=su.screensaver_user_id'),
        migrations.RunSQL(
            'UPDATE screensaver_user su '
            ' set lab_head_appointment_update_date=lh.lab_head_appointment_update_date '
            ' from  lab_head lh '
            ' where lh.screensaver_user_id=su.screensaver_user_id'),
        # end TODO: 20170918 ======

        migrations.CreateModel(
            name='UserAgreement',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.TextField()),
                ('data_sharing_level', models.IntegerField(null=True)),
                ('date_active', models.DateField(null=True)),
                ('date_expired', models.DateField(null=True)),
                ('date_notified', models.DateField(null=True)),
                ('file', models.ForeignKey(to='db.AttachedFile', null=True)),
            ],
            options={
                'db_table': 'user_agreement',
            },
        ),
        migrations.AddField(
            model_name='useragreement',
            name='screensaver_user',
            field=models.ForeignKey(to='db.ScreensaverUser'),
        ),
        migrations.AlterUniqueTogether(
            name='useragreement',
            unique_together=set([('screensaver_user', 'type')]),
        ),

        
    ]
