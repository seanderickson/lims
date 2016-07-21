# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import django.utils.timezone
# import db.models
# import django.db.models.deletion

## Note: existing migrations are intended for migrating an existing screensaver 
## database to the new database schema.
## Therefore, these should not be used to create a new instance.
## instead, the migrations directory should be removed and recreated with 
## $> manage.py makemigrations db
## which will generate an initial migration by inspecting the 
## models.py file, and then run:
## $> manage.py migrate --fake-initial

## Developer note:
## New models and fields to be added should not be included in the 0001 file, 
## they should be applied in subsequent files, as migration 0001 is "faked" in
## the migration.sh script.

class Migration(migrations.Migration):

    dependencies = [
        ('reports', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='AbaseTestset',
            fields=[
                ('abase_testset_id', models.IntegerField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('testset_name', models.TextField()),
                ('comments', models.TextField()),
                ('testset_date', models.DateField()),
            ],
            options={
                'db_table': 'abase_testset',
            },
        ),
        migrations.CreateModel(
            name='Activity',
            fields=[
                ('activity_id', models.AutoField(serialize=False, primary_key=True)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('comments', models.TextField(blank=True)),
                ('date_of_activity', models.DateField()),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
                ('version', models.IntegerField()),
            ],
            options={
                'db_table': 'activity',
            },
        ),
        migrations.CreateModel(
            name='ActivityUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'db_table': 'activity_update_activity',
            },
        ),
        migrations.CreateModel(
            name='AnnotationType',
            fields=[
                ('annotation_type_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('name', models.TextField(blank=True)),
                ('description', models.TextField(blank=True)),
                ('ordinal', models.IntegerField()),
                ('is_numeric', models.BooleanField()),
            ],
            options={
                'db_table': 'annotation_type',
            },
        ),
        migrations.CreateModel(
            name='AnnotationValue',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('annotation_value_id', models.IntegerField(null=True, blank=True)),
                ('numeric_value', models.FloatField(null=True, blank=True)),
                ('value', models.TextField(blank=True)),
                ('annotation_type', models.ForeignKey(blank=True, to='db.AnnotationType', null=True)),
            ],
            options={
                'db_table': 'annotation_value',
            },
        ),
        migrations.CreateModel(
            name='AssayPlate',
            fields=[
                ('assay_plate_id', models.AutoField(serialize=False, primary_key=True)),
                ('replicate_ordinal', models.IntegerField()),
                ('version', models.IntegerField()),
                ('plate_number', models.IntegerField()),
            ],
            options={
                'db_table': 'assay_plate',
            },
        ),
        migrations.CreateModel(
            name='AssayWell',
            fields=[
                ('assay_well_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('assay_well_control_type', models.TextField(null=True,blank=True)),
                ('is_positive', models.BooleanField()),
                ('confirmed_positive_value', models.TextField(null=True,blank=True)),
                ('plate_number', models.IntegerField()),
            ],
            options={
                'db_table': 'assay_well',
            },
        ),
        migrations.CreateModel(
            name='AttachedFile',
            fields=[
                ('attached_file_id', models.AutoField(serialize=False, primary_key=True)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('contents', models.BinaryField()),
                ('filename', models.TextField()),
#                 ('type', models.TextField()),
                ('file_date', models.DateField(null=True, blank=True)),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
                ('version', models.IntegerField()),
            ],
            options={
                'db_table': 'attached_file',
            },
        ),
        migrations.CreateModel(
            name='AttachedFileType',
            fields=[
                ('for_entity_type', models.CharField(max_length=31)),
                ('attached_file_type_id', models.IntegerField(serialize=False, primary_key=True)),
                ('value', models.TextField()),
            ],
            options={
                'db_table': 'attached_file_type',
            },
        ),
        migrations.CreateModel(
            name='AttachedFileUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('attached_file', models.ForeignKey(to='db.AttachedFile')),
            ],
            options={
                'db_table': 'attached_file_update_activity',
            },
        ),
        
        
#         migrations.CreateModel(
#             name='CachedQuery',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#                 ('key', models.TextField(unique=True)),
#                 ('sql', models.TextField()),
#                 ('uri', models.TextField()),
#                 ('params', models.TextField(null=True)),
#                 ('datetime', models.DateTimeField(default=django.utils.timezone.now)),
#                 ('username', models.CharField(max_length=128)),
#                 ('count', models.IntegerField(null=True)),
#             ],
#             options={
#                 'db_table': 'cached_query',
#             },
#         ),
        migrations.CreateModel(
            name='CellLine',
            fields=[
                ('cell_line_id', models.IntegerField(serialize=False, primary_key=True)),
                ('value', models.TextField(unique=True)),
                ('version', models.IntegerField()),
            ],
            options={
                'db_table': 'cell_line',
            },
        ),
        migrations.CreateModel(
            name='ChecklistItem',
            fields=[
                ('checklist_item_id', models.IntegerField(serialize=False, primary_key=True)),
                ('checklist_item_group', models.TextField()),
                ('is_expirable', models.BooleanField()),
                ('item_name', models.TextField(unique=True)),
                ('order_statistic', models.IntegerField()),
                ('version', models.IntegerField()),
            ],
            options={
                'db_table': 'checklist_item',
            },
        ),
        migrations.CreateModel(
            name='ChecklistItemEvent',
            fields=[
                ('checklist_item_event_id', models.IntegerField(serialize=False, primary_key=True)),
                ('date_performed', models.DateField(null=True, blank=True)),
                ('is_expiration', models.BooleanField()),
                ('checklist_item_id', models.IntegerField()),
                ('is_not_applicable', models.BooleanField()),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
            ],
            options={
                'db_table': 'checklist_item_event',
            },
        ),
        migrations.CreateModel(
            name='ChecklistItemEventUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('checklist_item_event', models.ForeignKey(to='db.ChecklistItemEvent')),
            ],
            options={
                'db_table': 'checklist_item_event_update_activity',
            },
        ),
        migrations.CreateModel(
            name='CherryPickAssayPlate',
            fields=[
                ('cherry_pick_assay_plate_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('plate_ordinal', models.IntegerField()),
                ('attempt_ordinal', models.IntegerField()),
                ('assay_plate_type', models.TextField()),
                ('legacy_plate_name', models.TextField(blank=True)),
                ('cherry_pick_assay_plate_type', models.CharField(max_length=31)),
#                 ('status', models.TextField(null=True)),
            ],
            options={
                'db_table': 'cherry_pick_assay_plate',
            },
        ),
        migrations.CreateModel(
            name='CherryPickAssayPlateScreeningLink',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('cherry_pick_assay_plate', models.ForeignKey(to='db.CherryPickAssayPlate')),
            ],
            options={
                'db_table': 'cherry_pick_assay_plate_screening_link',
            },
        ),
        migrations.CreateModel(
            name='CherryPickRequest',
            fields=[
                ('cherry_pick_request_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('comments', models.TextField(blank=True)),
                ('is_randomized_assay_plate_layout', models.BooleanField()),
                ('legacy_cherry_pick_request_number', models.IntegerField(null=True, blank=True)),
                ('number_unfulfilled_lab_cherry_picks', models.IntegerField()),
                ('assay_plate_type', models.TextField()),
                ('transfer_volume_per_well_approved', models.DecimalField(null=True, max_digits=10, decimal_places=9, blank=True)),
                ('transfer_volume_per_well_requested', models.DecimalField(null=True, max_digits=10, decimal_places=9, blank=True)),
                ('date_requested', models.DateField()),
                ('date_volume_approved', models.DateField(null=True, blank=True)),
                ('assay_protocol_comments', models.TextField(blank=True)),
                ('cherry_pick_assay_protocols_followed', models.TextField(blank=True)),
                ('cherry_pick_followup_results_status', models.TextField(blank=True)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('keep_source_plate_cherry_picks_together', models.BooleanField()),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
                ('max_skipped_wells_per_plate', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'db_table': 'cherry_pick_request',
            },
        ),
        migrations.CreateModel(
            name='CherryPickRequestEmptyWell',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('well_name', models.CharField(max_length=255, blank=True)),
            ],
            options={
                'db_table': 'cherry_pick_request_empty_well',
            },
        ),
        migrations.CreateModel(
            name='CherryPickRequestUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'db_table': 'cherry_pick_request_update_activity',
            },
        ),
#         migrations.CreateModel(
#             name='CollaboratorLink',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#             ],
#             options={
#                 'db_table': 'collaborator_link',
#             },
#         ),
        migrations.CreateModel(
            name='Copy',
            fields=[
                ('version', models.IntegerField()),
                ('usage_type', models.TextField()),
                ('name', models.TextField()),
                ('copy_id', models.AutoField(serialize=False, primary_key=True)),
                ('comments', models.TextField(blank=True)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('date_plated', models.DateField(null=True, blank=True)),
                ('primary_plate_status', models.TextField()),
                ('primary_plate_location_id', models.IntegerField(null=True, blank=True)),
                ('plates_available', models.IntegerField(null=True, blank=True)),
                ('plate_locations_count', models.IntegerField(null=True, blank=True)),
                ('primary_well_mg_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
                ('primary_well_molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('min_molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('max_molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('min_mg_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
                ('max_mg_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
                ('well_concentration_dilution_factor', models.DecimalField(null=True, max_digits=8, decimal_places=2, blank=True)),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
            ],
            options={
                'db_table': 'copy',
            },
        ),
        migrations.CreateModel(
            name='CopyUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('copy_id', models.IntegerField()),
                ('update_activity_id', models.IntegerField(unique=True)),
            ],
            options={
                'db_table': 'copy_update_activity',
            },
        ),
#         migrations.CreateModel(
#             name='CopyWell',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#                 ('plate_number', models.IntegerField()),
#                 ('volume', models.FloatField(null=True, blank=True)),
#                 ('initial_volume', models.FloatField(null=True, blank=True)),
#                 ('adjustments', models.IntegerField()),
#                 ('copy', models.ForeignKey(to='db.Copy')),
#             ],
#             options={
#                 'db_table': 'copy_well',
#             },
#         ),
        migrations.CreateModel(
            name='DataColumn',
            fields=[
                ('data_column_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('ordinal', models.IntegerField()),
                ('replicate_ordinal', models.IntegerField(null=True, blank=True)),
                ('assay_phenotype', models.TextField(blank=True)),
                ('assay_readout_type', models.TextField(blank=True)),
                ('comments', models.TextField(blank=True)),
                ('description', models.TextField(blank=True)),
                ('how_derived', models.TextField(blank=True)),
                ('is_follow_up_data', models.BooleanField()),
                ('name', models.TextField()),
                ('time_point', models.TextField(blank=True)),
                ('is_derived', models.BooleanField()),
                ('positives_count', models.IntegerField(null=True, blank=True)),
                ('channel', models.IntegerField(null=True, blank=True)),
                ('time_point_ordinal', models.IntegerField(null=True, blank=True)),
                ('zdepth_ordinal', models.IntegerField(null=True, blank=True)),
                ('data_type', models.TextField()),
                ('decimal_places', models.IntegerField(null=True, blank=True)),
                ('strong_positives_count', models.IntegerField(null=True, blank=True)),
                ('medium_positives_count', models.IntegerField(null=True, blank=True)),
                ('weak_positives_count', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'db_table': 'data_column',
            },
        ),
        migrations.CreateModel(
            name='DataColumnDerivedFromLink',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('derived_data_column', models.ForeignKey(to='db.DataColumn')),
                ('derived_from_data_column', models.ForeignKey(related_name='derived_from', to='db.DataColumn')),
            ],
            options={
                'db_table': 'data_column_derived_from_link',
            },
        ),
        migrations.CreateModel(
            name='EquipmentUsed',
            fields=[
                ('equipment_used_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('protocol', models.TextField(blank=True)),
                ('description', models.TextField(blank=True)),
                ('equipment', models.TextField()),
            ],
            options={
                'db_table': 'equipment_used',
            },
        ),
        migrations.CreateModel(
            name='FundingSupport',
            fields=[
                ('funding_support_id', models.IntegerField(serialize=False, primary_key=True)),
                ('value', models.TextField(unique=True, blank=True)),
            ],
            options={
                'db_table': 'funding_support',
            },
        ),
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('gene_id', models.AutoField(serialize=False, primary_key=True)),
                ('entrezgene_id', models.IntegerField(null=True, blank=True)),
                ('gene_name', models.TextField(blank=True)),
                ('species_name', models.TextField(blank=True)),
            ],
            options={
                'db_table': 'gene',
            },
        ),
        migrations.CreateModel(
            name='GeneGenbankAccessionNumber',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('genbank_accession_number', models.TextField()),
                ('gene', models.ForeignKey(to='db.Gene')),
            ],
            options={
                'db_table': 'gene_genbank_accession_number',
            },
        ),
        migrations.CreateModel(
            name='GeneSymbol',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('entrezgene_symbol', models.TextField()),
                ('ordinal', models.IntegerField()),
                ('gene', models.ForeignKey(to='db.Gene')),
            ],
            options={
                'db_table': 'gene_symbol',
            },
        ),
        migrations.CreateModel(
            name='LabAffiliation',
            fields=[
                ('version', models.IntegerField()),
                ('affiliation_name', models.TextField(unique=True)),
                ('affiliation_category', models.TextField()),
                ('lab_affiliation_id', models.AutoField(serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'lab_affiliation',
            },
        ),
        migrations.CreateModel(
            name='LabCherryPick',
            fields=[
                ('lab_cherry_pick_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('assay_plate_row', models.IntegerField(null=True, blank=True)),
                ('assay_plate_column', models.IntegerField(null=True, blank=True)),
                ('cherry_pick_assay_plate', models.ForeignKey(blank=True, to='db.CherryPickAssayPlate', null=True)),
            ],
            options={
                'db_table': 'lab_cherry_pick',
            },
        ),
        migrations.CreateModel(
            name='LegacySmallMoleculeCasNumber',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('smiles', models.CharField(max_length=2047)),
                ('cas_number', models.TextField()),
            ],
            options={
                'db_table': '_legacy_small_molecule_cas_number',
            },
        ),
        migrations.CreateModel(
            name='Library',
            fields=[
                ('library_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField(null=True, blank=True)),
                ('library_name', models.TextField(unique=True)),
                ('short_name', models.TextField(unique=True)),
                ('description', models.TextField(blank=True)),
                ('provider', models.TextField(blank=True)),
                ('screen_type', models.TextField()),
                ('library_type', models.TextField()),
                ('start_plate', models.IntegerField(unique=True)),
                ('end_plate', models.IntegerField(unique=True)),
                ('screening_status', models.TextField()),
                ('date_received', models.DateField(null=True, blank=True)),
                ('date_screenable', models.DateField(null=True, blank=True)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('plate_size', models.TextField()),
                ('latest_released_contents_version_id', models.IntegerField(null=True, blank=True)),
                ('experimental_well_count', models.IntegerField(null=True, blank=True)),
                ('is_pool', models.NullBooleanField()),
                ('solvent', models.TextField()),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
#                 ('version_number', models.IntegerField(default=0)),
            ],
            options={
                'db_table': 'library',
            },
        ),
        migrations.CreateModel(
            name='LibraryContentsVersion',
            fields=[
                ('library_contents_version_id', models.IntegerField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('version_number', models.IntegerField()),
                ('library', models.ForeignKey(to='db.Library')),
            ],
            options={
                'db_table': 'library_contents_version',
            },
        ),
        migrations.CreateModel(
            name='LibraryUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('library', models.ForeignKey(to='db.Library')),
            ],
            options={
                'db_table': 'library_update_activity',
            },
        ),
        migrations.CreateModel(
            name='Plate',
            fields=[
                ('plate_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('plate_type', models.TextField(blank=True)),
                ('plate_number', models.IntegerField()),
                ('well_volume', models.DecimalField(null=True, max_digits=10, decimal_places=9, blank=True)),
#                 ('remaining_volume', models.FloatField(null=True, blank=True)),
#                 ('avg_remaining_volume', models.FloatField(null=True, blank=True)),
#                 ('min_remaining_volume', models.FloatField(null=True, blank=True)),
#                 ('max_remaining_volume', models.FloatField(null=True, blank=True)),
#                 ('screening_count', models.IntegerField(null=True, blank=True)),
                ('facility_id', models.TextField(blank=True)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('status', models.TextField()),
                ('retired_activity_id', models.IntegerField(unique=True, null=True, blank=True)),
                ('plated_activity_id', models.IntegerField(unique=True, null=True, blank=True)),
                ('stock_plate_number', models.IntegerField(null=True, blank=True)),
                ('quadrant', models.IntegerField(null=True, blank=True)),
                ('min_molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('max_molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('min_mg_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
                ('max_mg_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
                ('primary_well_molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('primary_well_mg_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
                ('copy', models.ForeignKey(to='db.Copy')),
            ],
            options={
                'db_table': 'plate',
            },
        ),
        migrations.CreateModel(
            name='PlateLocation',
            fields=[
                ('plate_location_id', models.AutoField(serialize=False, primary_key=True)),
                ('bin', models.TextField()),
                ('freezer', models.TextField()),
                ('room', models.TextField()),
                ('shelf', models.TextField()),
            ],
            options={
                'db_table': 'plate_location',
            },
        ),
        migrations.CreateModel(
            name='PlateUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('plate_id', models.IntegerField()),
                ('update_activity_id', models.IntegerField(unique=True)),
            ],
            options={
                'db_table': 'plate_update_activity',
            },
        ),
        migrations.CreateModel(
            name='Publication',
            fields=[
                ('publication_id', models.AutoField(serialize=False, primary_key=True)),
                ('authors', models.TextField(blank=True)),
                ('journal', models.TextField(blank=True)),
                ('pages', models.TextField(blank=True)),
                ('pubmed_id', models.TextField(blank=True)),
                ('title', models.TextField(blank=True)),
                ('version', models.IntegerField()),
                ('volume', models.TextField(blank=True)),
                ('year_published', models.TextField(blank=True)),
                ('pubmed_central_id', models.IntegerField(null=True, blank=True)),
#                 ('attached_file', models.OneToOneField(
#                     null=True, to='db.AttachedFile', unique=True, on_delete=models.CASCADE)),
            ],
            options={
                'db_table': 'publication',
            },
        ),
        migrations.CreateModel(
            name='Reagent',
            fields=[
                ('reagent_id', models.AutoField(serialize=False, primary_key=True)),
#                 ('substance_id', models.CharField(default=db.models.create_id, unique=True, max_length=8)),
                ('vendor_identifier', models.TextField(blank=True)),
                ('vendor_name', models.TextField(blank=True)),
                ('vendor_batch_id', models.TextField(blank=True)),
                ('facility_batch_id', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'db_table': 'reagent',
            },
        ),
#         migrations.CreateModel(
#             name='ReagentPublicationLink',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#                 ('publication_id', models.IntegerField(unique=True)),
#             ],
#             options={
#                 'db_table': 'reagent_publication_link',
#             },
#         ),
        migrations.CreateModel(
            name='ResultValue',
            fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('result_value_id', models.AutoField(serialize=False, primary_key=True)),
                ('assay_well_control_type', models.TextField(null=True, blank=True)),
                ('is_exclude', models.NullBooleanField()),
                ('numeric_value', models.FloatField(null=True, blank=True)),
                ('is_positive', models.NullBooleanField()),
                ('value', models.TextField(null=True,blank=True)),
                ('data_column', models.ForeignKey(blank=True, to='db.DataColumn', null=True)),
            ],
            options={
                'db_table': 'result_value',
            },
        ),
        migrations.CreateModel(
            name='SchemaHistory',
            fields=[
                ('screensaver_revision', models.IntegerField(serialize=False, primary_key=True)),
                ('date_updated', models.DateTimeField(null=True, blank=True)),
                ('comment', models.TextField(blank=True)),
            ],
            options={
                'db_table': 'schema_history',
            },
        ),
        migrations.CreateModel(
            name='Screen',
            fields=[
                ('screen_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('screen_type', models.TextField()),
                ('title', models.TextField()),
                ('summary', models.TextField(blank=True)),
                ('comments', models.TextField(blank=True)),
                ('abase_study_id', models.TextField(blank=True)),
                ('abase_protocol_id', models.TextField(blank=True)),
                ('publishable_protocol', models.TextField(blank=True)),
                ('publishable_protocol_entered_by', models.TextField(blank=True)),
                ('publishable_protocol_comments', models.TextField(blank=True)),
                ('study_type', models.TextField()),
                ('url', models.TextField(blank=True)),
                ('data_meeting_complete', models.DateField(null=True, blank=True)),
                ('data_meeting_scheduled', models.DateField(null=True, blank=True)),
                ('date_of_application', models.DateField(null=True, blank=True)),
                ('publishable_protocol_date_entered', models.DateField(null=True, blank=True)),
                ('amount_to_be_charged_for_screen', models.DecimalField(null=True, max_digits=9, decimal_places=2, blank=True)),
                ('billing_comments', models.TextField(blank=True)),
                ('is_billing_for_supplies_only', models.BooleanField(default=False)),
                ('billing_info_return_date', models.DateField(null=True, blank=True)),
                ('date_charged', models.DateField(null=True, blank=True)),
                ('date_completed5kcompounds', models.DateField(null=True, blank=True)),
                ('date_faxed_to_billing_department', models.DateField(null=True, blank=True)),
                ('facilities_and_administration_charge', models.DecimalField(null=True, max_digits=9, decimal_places=2, blank=True)),
                ('is_fee_form_on_file', models.BooleanField(default=False)),
                ('fee_form_requested_date', models.DateField(null=True, blank=True)),
                ('fee_form_requested_initials', models.TextField(blank=True)),
                ('see_comments', models.BooleanField(default=False)),
                ('to_be_requested', models.BooleanField(default=False)),
                ('coms_registration_number', models.TextField(blank=True)),
                ('coms_approval_date', models.DateField(null=True, blank=True)),
                ('data_sharing_level', models.IntegerField()),
                ('data_privacy_expiration_date', models.DateField(null=True, blank=True)),
                ('max_allowed_data_privacy_expiration_date', models.DateField(null=True, blank=True)),
                ('min_allowed_data_privacy_expiration_date', models.DateField(null=True, blank=True)),
                ('data_privacy_expiration_notified_date', models.DateField(null=True, blank=True)),
                ('screened_experimental_well_count', models.IntegerField(default=0)),
                ('unique_screened_experimental_well_count', models.IntegerField(default=0)),
                ('total_plated_lab_cherry_picks', models.IntegerField(default=0)),
                ('assay_plates_screened_count', models.IntegerField(default=0)),
                ('library_plates_screened_count', models.IntegerField(default=0)),
                ('library_plates_data_loaded_count', models.IntegerField(default=0)),
                ('library_plates_data_analyzed_count', models.IntegerField(default=0)),
                ('min_screened_replicate_count', models.IntegerField(null=True, blank=True)),
                ('max_screened_replicate_count', models.IntegerField(null=True, blank=True)),
                ('min_data_loaded_replicate_count', models.IntegerField(null=True, blank=True)),
                ('max_data_loaded_replicate_count', models.IntegerField(null=True, blank=True)),
                ('libraries_screened_count', models.IntegerField(null=True, blank=True)),
                ('facility_id', models.TextField(unique=True)),
                ('project_phase', models.TextField()),
                ('project_id', models.TextField(blank=True)),
                ('pubchem_assay_id', models.IntegerField(null=True, blank=True)),
                ('pubchem_deposited_date', models.DateField(null=True, blank=True)),
                ('image_url', models.TextField(blank=True)),
                ('species', models.TextField(blank=True)),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
                ('perturbagen_molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('perturbagen_ug_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
#                 ('status', models.TextField(null=True, blank=True)),
#                 ('status_date', models.DateField(null=True, blank=True)),
                ('assay_type', models.TextField(null=True, blank=True)),
            ],
            options={
                'db_table': 'screen',
            },
        ),
        migrations.CreateModel(
            name='ScreenBillingItem',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('amount', models.DecimalField(max_digits=9, decimal_places=2)),
                ('date_sent_for_billing', models.DateField(null=True, blank=True)),
                ('item_to_be_charged', models.TextField()),
                ('ordinal', models.IntegerField()),
                ('screen', models.ForeignKey(to='db.Screen')),
            ],
            options={
                'db_table': 'screen_billing_item',
            },
        ),
        migrations.CreateModel(
            name='ScreenerCherryPick',
            fields=[
                ('screener_cherry_pick_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
            ],
            options={
                'db_table': 'screener_cherry_pick',
            },
        ),
        migrations.CreateModel(
            name='ScreenFundingSupportLink',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('funding_support', models.ForeignKey(to='db.FundingSupport')),
                ('screen', models.ForeignKey(to='db.Screen')),
            ],
            options={
                'db_table': 'screen_funding_support_link',
            },
        ),
#         migrations.CreateModel(
#             name='ScreenFundingSupports',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#                 ('funding_support', models.TextField()),
#                 ('screen', models.ForeignKey(to='db.Screen')),
#             ],
#             options={
#                 'db_table': 'screen_funding_supports',
#             },
#         ),
        migrations.CreateModel(
            name='ScreeningRoomUserFacilityUsageRole',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('facility_usage_role', models.TextField()),
            ],
            options={
                'db_table': 'screening_room_user_facility_usage_role',
            },
        ),
        migrations.CreateModel(
            name='ScreenKeyword',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('keyword', models.TextField()),
                ('screen', models.ForeignKey(to='db.Screen', related_name='keywords')),
            ],
            options={
                'db_table': 'screen_keyword',
            },
        ),
#         migrations.CreateModel(
#             name='ScreenPublicationLink',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#                 ('publication_id', models.IntegerField(unique=True)),
#                 ('screen', models.ForeignKey(to='db.Screen')),
#             ],
#             options={
#                 'db_table': 'screen_publication_link',
#             },
#         ),
        migrations.CreateModel(
            name='ScreenResult',
            fields=[
                ('screen_result_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('replicate_count', models.IntegerField()),
                ('experimental_well_count', models.IntegerField()),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('channel_count', models.IntegerField(null=True, blank=True)),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
            ],
            options={
                'db_table': 'screen_result',
            },
        ),
        migrations.CreateModel(
            name='ScreenResultUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('screen_result', models.ForeignKey(to='db.ScreenResult')),
            ],
            options={
                'db_table': 'screen_result_update_activity',
            },
        ),
        migrations.CreateModel(
            name='ScreensaverUser',
            fields=[
                ('screensaver_user_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField(default=1, blank=True)),
                ('date_created', models.DateTimeField(default=django.utils.timezone.now)),
                ('harvard_id', models.TextField(blank=True)),
                ('harvard_id_expiration_date', models.DateField(null=True, blank=True)),
                ('harvard_id_requested_expiration_date', models.DateField(null=True, blank=True)),
                ('first_name', models.TextField()),
                ('last_name', models.TextField()),
                ('email', models.TextField(blank=True)),
                ('phone', models.TextField(blank=True)),
                ('mailing_address', models.TextField(blank=True)),
                ('ecommons_id', models.TextField(blank=True)),
                ('date_loaded', models.DateTimeField(null=True, blank=True)),
                ('date_publicly_available', models.DateTimeField(null=True, blank=True)),
                ('login_id', models.TextField(null=True)),
                ('digested_password', models.TextField(null=True)),
                ('comments', models.TextField(blank=True)),
#                 ('username', models.TextField(unique=True, null=True)),
            ],
            options={
                'db_table': 'screensaver_user',
            },
        ),
        migrations.CreateModel(
            name='ScreensaverUserRole',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('screensaver_user_role', models.TextField()),
            ],
            options={
                'db_table': 'screensaver_user_role',
            },
        ),
        migrations.CreateModel(
            name='ScreensaverUserUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'db_table': 'screensaver_user_update_activity',
            },
        ),
        migrations.CreateModel(
            name='ScreenStatusItem',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('status', models.TextField()),
                ('status_date', models.DateField()),
                ('screen', models.ForeignKey(to='db.Screen')),
            ],
            options={
                'db_table': 'screen_status_item',
            },
        ),
        migrations.CreateModel(
            name='ScreenUpdateActivity',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('screen', models.ForeignKey(to='db.Screen')),
            ],
            options={
                'db_table': 'screen_update_activity',
            },
        ),
        migrations.CreateModel(
            name='SmallMoleculeChembankId',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chembank_id', models.IntegerField()),
            ],
            options={
                'db_table': 'small_molecule_chembank_id',
            },
        ),
        migrations.CreateModel(
            name='SmallMoleculeChemblId',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('chembl_id', models.IntegerField()),
            ],
            options={
                'db_table': 'small_molecule_chembl_id',
            },
        ),
        migrations.CreateModel(
            name='SmallMoleculeCompoundName',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('compound_name', models.TextField()),
                ('ordinal', models.IntegerField()),
            ],
            options={
                'db_table': 'small_molecule_compound_name',
            },
        ),
        migrations.CreateModel(
            name='SmallMoleculePubchemCid',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('pubchem_cid', models.IntegerField()),
            ],
            options={
                'db_table': 'small_molecule_pubchem_cid',
            },
        ),
        migrations.CreateModel(
            name='StudyReagentLink',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
            ],
            options={
                'db_table': 'study_reagent_link',
            },
        ),
#         migrations.CreateModel(
#             name='Substance',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#                 ('comment', models.TextField()),
#             ],
#         ),
        migrations.CreateModel(
            name='TransfectionAgent',
            fields=[
                ('transfection_agent_id', models.IntegerField(serialize=False, primary_key=True)),
                ('value', models.TextField(unique=True)),
                ('version', models.IntegerField()),
            ],
            options={
                'db_table': 'transfection_agent',
            },
        ),
#         migrations.CreateModel(
#             name='UserChecklistItem',
#             fields=[
#                 ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
#                 ('item_group', models.TextField()),
#                 ('item_name', models.TextField()),
#                 ('status', models.TextField()),
#                 ('status_date', models.DateField()),
#             ],
#             options={
#                 'db_table': 'user_checklist_item',
#             },
#         ),
        migrations.CreateModel(
            name='Well',
            fields=[
                ('well_id', models.TextField(serialize=False, primary_key=True)),
                ('version', models.IntegerField(null=True, blank=True)),
                ('plate_number', models.IntegerField()),
                ('well_name', models.TextField()),
                ('facility_id', models.TextField(blank=True)),
                ('library_well_type', models.TextField()),
                ('is_deprecated', models.BooleanField(default=False)),
                ('molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
                ('mg_ml_concentration', models.DecimalField(null=True, max_digits=5, decimal_places=3, blank=True)),
                ('barcode', models.TextField(unique=True, null=True)),
            ],
            options={
                'db_table': 'well',
            },
        ),
        migrations.CreateModel(
            name='WellVolumeAdjustment',
            fields=[
                ('well_volume_adjustment_id', models.AutoField(serialize=False, primary_key=True)),
                ('version', models.IntegerField()),
                ('volume', models.DecimalField(null=True, max_digits=10, decimal_places=9, blank=True)),
                ('copy', models.ForeignKey(to='db.Copy')),
                ('lab_cherry_pick', models.ForeignKey(blank=True, to='db.LabCherryPick', null=True)),
                ('well', models.ForeignKey(to='db.Well')),
            ],
            options={
                'db_table': 'well_volume_adjustment',
            },
        ),
        migrations.CreateModel(
            name='AdministrativeActivity',
            fields=[
                ('activity', models.OneToOneField(primary_key=True, serialize=False, to='db.Activity')),
                ('administrative_activity_type', models.TextField()),
            ],
            options={
                'db_table': 'administrative_activity',
            },
        ),
        migrations.CreateModel(
            name='AdministratorUser',
            fields=[
                ('screensaver_user', models.OneToOneField(
                    primary_key=True, serialize=False, to='db.ScreensaverUser')),
            ],
            options={
                'db_table': 'administrator_user',
            },
        ),
        migrations.CreateModel(
            name='LabActivity',
            fields=[
                ('activitylink', models.OneToOneField(
                    primary_key=True, parent_link=True, db_column='activity_id', 
                    serialize=False, to='db.Activity')),
                ('volume_transferred_per_well_from_library_plates', models.DecimalField(null=True, max_digits=10, decimal_places=9, blank=True)),
                ('molar_concentration', models.DecimalField(null=True, max_digits=13, decimal_places=12, blank=True)),
            ],
            options={
                'db_table': 'lab_activity',
            },
        ),
        migrations.CreateModel(
            name='Molfile',
            fields=[
                ('molfile', models.TextField()),
                ('reagent', models.OneToOneField(primary_key=True, serialize=False, to='db.Reagent')),
            ],
            options={
                'db_table': 'molfile',
            },
        ),
        migrations.CreateModel(
            name='RnaiCherryPickRequest',
            fields=[
                ('cherry_pick_request', models.OneToOneField(primary_key=True, serialize=False, to='db.CherryPickRequest')),
            ],
            options={
                'db_table': 'rnai_cherry_pick_request',
            },
        ),
        migrations.CreateModel(
            name='ScreeningRoomUser',
            fields=[
                ('screensaver_user', models.OneToOneField(primary_key=True, 
                    serialize=False, to='db.ScreensaverUser')),
                ('user_classification', models.TextField()),
                ('coms_crhba_permit_number', models.TextField(blank=True)),
                ('coms_crhba_permit_principal_investigator', models.TextField(blank=True)),
            ],
            options={
                'db_table': 'screening_room_user',
            },
        ),
        migrations.CreateModel(
            name='ServiceActivity',
            fields=[
                ('service_activity_type', models.TextField()),
                ('activity', models.OneToOneField(primary_key=True, serialize=False, to='db.Activity')),
#                 ('funding_support', models.TextField(null=True)),
                ('funding_support_link', models.ForeignKey(db_column='funding_support_id', to='db.FundingSupport', null=True)),
            ],
            options={
                'db_table': 'service_activity',
            },
        ),
        migrations.CreateModel(
            name='SilencingReagent',
            fields=[
                ('reagentlink', models.OneToOneField(
                    primary_key=True, parent_link=True,db_column='reagent_id',
                    serialize=False, to='db.Reagent')),
#                 ('reagent', models.OneToOneField(primary_key=True, serialize=False, to='db.Reagent')),
                ('sequence', models.TextField(blank=True)),
                ('anti_sense_sequence', models.TextField(blank=True)),
                ('silencing_reagent_type', models.TextField(blank=True)),
                ('is_restricted_sequence', models.BooleanField(default=False)),
            ],
            options={
                'db_table': 'silencing_reagent',
            },
        ),
        migrations.CreateModel(
            name='SmallMoleculeReagent',
            fields=[
                ('reagentlink', models.OneToOneField(
                    primary_key=True, parent_link=True,db_column='reagent_id',
                    serialize=False, to='db.Reagent')),
#                 ('reagent', models.OneToOneField(primary_key=True, serialize=False, to='db.Reagent')),
                ('inchi', models.TextField(blank=True)),
                ('molecular_formula', models.TextField(blank=True)),
                ('molecular_mass', models.DecimalField(null=True, max_digits=15, decimal_places=9, blank=True)),
                ('molecular_weight', models.DecimalField(null=True, max_digits=15, decimal_places=9, blank=True)),
                ('smiles', models.TextField(blank=True)),
                ('is_restricted_structure', models.BooleanField(default=False)),
                ('salt_form_id', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'db_table': 'small_molecule_reagent',
            },
        ),
        migrations.CreateModel(
            name='NaturalProductReagent',
            fields=[
                ('reagentlink', models.OneToOneField(
                    primary_key=True, parent_link=True,db_column='reagent_id',
                    serialize=False, to='db.Reagent')),
#                 ('reagent', models.OneToOneField(primary_key=True, serialize=False, to='db.Reagent')),
            ],
            options={
                'db_table': 'natural_product_reagent',
            },
        ),
        migrations.CreateModel(
            name='SmallMoleculeCherryPickRequest',
            fields=[
                ('cherry_pick_request', models.OneToOneField(primary_key=True, serialize=False, to='db.CherryPickRequest')),
            ],
            options={
                'db_table': 'small_molecule_cherry_pick_request',
            },
        ),
        # New datacolumn_derived_from_columns table, see manual/0002 sql for creation
        migrations.AddField(
            model_name='datacolumn',
            name='derived_from_columns',
            field=models.ManyToManyField(related_name='derived_columns', to='db.DataColumn'),
        ),
        migrations.AddField(
            model_name='well',
            name='latest_released_reagent',
            field=models.ForeignKey(related_name='reagent_well', blank=True, to='db.Reagent', null=True),
        ),
        migrations.AddField(
            model_name='publication',
            name='screen',
            field=models.ForeignKey(to='db.Screen', null=True, on_delete=models.CASCADE),
        ),
        migrations.AddField(
            model_name='publication',
            name='reagent',
            field=models.ForeignKey(to='db.Reagent', null=True, on_delete=models.CASCADE),
        ),
        migrations.AddField(
            model_name='well',
            name='library',
            field=models.ForeignKey(to='db.Library'),
        ),
#         migrations.AddField(
#             model_name='userchecklistitem',
#             name='admin_user',
#             field=models.ForeignKey(related_name='userchecklistitems_created', to='db.ScreensaverUser'),
#         ),
#         migrations.AddField(
#             model_name='userchecklistitem',
#             name='screensaver_user',
#             field=models.ForeignKey(to='db.ScreensaverUser'),
#         ),
        migrations.AddField(
            model_name='studyreagentlink',
            name='reagent',
            field=models.ForeignKey(to='db.Reagent'),
        ),
        migrations.AddField(
            model_name='studyreagentlink',
            name='study',
            field=models.ForeignKey(to='db.Screen'),
        ),
        migrations.AddField(
            model_name='smallmoleculepubchemcid',
            name='reagent',
            field=models.ForeignKey(to='db.Reagent'),
        ),
        migrations.AddField(
            model_name='smallmoleculecompoundname',
            name='reagent',
            field=models.ForeignKey(to='db.Reagent'),
        ),
        migrations.AddField(
            model_name='smallmoleculechemblid',
            name='reagent',
            field=models.ForeignKey(to='db.Reagent'),
        ),
        migrations.AddField(
            model_name='smallmoleculechembankid',
            name='reagent',
            field=models.ForeignKey(to='db.Reagent'),
        ),
        migrations.AddField(
            model_name='screensaveruserupdateactivity',
            name='screensaver_user',
            field=models.ForeignKey(to='db.ScreensaverUser'),
        ),
        migrations.AddField(
            model_name='screensaveruserrole',
            name='screensaver_user',
            field=models.ForeignKey(to='db.ScreensaverUser'),
        ),
        migrations.AddField(
            model_name='screensaveruser',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', 
                null=True, related_name='created_user'),
        ),
#         migrations.AddField(
#             model_name='screensaveruser',
#             name='user',
#             field=models.ForeignKey(on_delete=django.db.models.deletion.SET_NULL, to='reports.UserProfile', null=True),
#         ),
        migrations.AddField(
            model_name='screenresult',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='screenresult',
            name='screen',
            field=models.OneToOneField(to='db.Screen'),
        ),
        migrations.AddField(
            model_name='screenercherrypick',
            name='cherry_pick_request',
            field=models.ForeignKey(to='db.CherryPickRequest'),
        ),
        migrations.AddField(
            model_name='screenercherrypick',
            name='screened_well',
            field=models.ForeignKey(to='db.Well'),
        ),
        migrations.AddField(
            model_name='screen',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='screen',
            name='transfection_agent',
            field=models.ForeignKey(blank=True, to='db.TransfectionAgent', null=True),
        ),
        migrations.AddField(
            model_name='screen',
            name='well_studied',
            field=models.ForeignKey(blank=True, to='db.Well', null=True),
        ),
        migrations.AddField(
            model_name='resultvalue',
            name='well',
            field=models.ForeignKey(blank=True, to='db.Well', null=True),
        ),
#         migrations.AddField(
#             model_name='reagentpublicationlink',
#             name='reagent',
#             field=models.ForeignKey(to='db.Reagent'),
#         ),
        migrations.AddField(
            model_name='reagent',
            name='library_contents_version',
            field=models.ForeignKey(to='db.LibraryContentsVersion', null=True),
        ),
        migrations.AddField(
            model_name='reagent',
            name='well',
            field=models.ForeignKey(to='db.Well', null=True),
        ),
        migrations.AddField(
            model_name='plate',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='plate',
            name='plate_location',
            field=models.ForeignKey(blank=True, to='db.PlateLocation', null=True),
        ),
        migrations.AddField(
            model_name='library',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', null=True),
        ),
#         migrations.AddField(
#             model_name='library',
#             name='loaded_by',
#             field=models.ForeignKey(related_name='libraries_loaded', blank=True, to='db.ScreensaverUser', null=True),
#         ),
        migrations.AddField(
            model_name='labcherrypick',
            name='cherry_pick_request',
            field=models.ForeignKey(to='db.CherryPickRequest'),
        ),
        migrations.AddField(
            model_name='labcherrypick',
            name='screener_cherry_pick',
            field=models.ForeignKey(to='db.ScreenerCherryPick'),
        ),
        migrations.AddField(
            model_name='labcherrypick',
            name='source_well',
            field=models.ForeignKey(to='db.Well'),
        ),
        migrations.AddField(
            model_name='datacolumn',
            name='screen_result',
            field=models.ForeignKey(to='db.ScreenResult'),
        ),
#         migrations.AddField(
#             model_name='copywell',
#             name='plate',
#             field=models.ForeignKey(to='db.Plate'),
#         ),
#         migrations.AddField(
#             model_name='copywell',
#             name='well',
#             field=models.ForeignKey(to='db.Well'),
#         ),
        migrations.AddField(
            model_name='copy',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='copy',
            name='library',
            field=models.ForeignKey(to='db.Library'),
        ),
#         migrations.AddField(
#             model_name='collaboratorlink',
#             name='screen',
#             field=models.ForeignKey(to='db.Screen'),
#         ),
        migrations.AddField(
            model_name='cherrypickrequestupdateactivity',
            name='cherry_pick_request',
            field=models.ForeignKey(to='db.CherryPickRequest'),
        ),
        migrations.AddField(
            model_name='cherrypickrequestemptywell',
            name='cherry_pick_request',
            field=models.ForeignKey(to='db.CherryPickRequest'),
        ),
        migrations.AddField(
            model_name='cherrypickrequest',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', 
                null=True, related_name='created_cherry_pick'),
        ),
        migrations.AddField(
            model_name='cherrypickrequest',
            name='screen',
            field=models.ForeignKey(to='db.Screen'),
        ),
        migrations.AddField(
            model_name='cherrypickassayplate',
            name='cherry_pick_request',
            field=models.ForeignKey(to='db.CherryPickRequest'),
        ),
        migrations.AddField(
            model_name='checklistitemevent',
            name='created_by',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='attachedfile',
            name='created_by',
            field=models.ForeignKey(related_name='attachedfilecreated', to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='attachedfile',
            name='reagent',
            field=models.ForeignKey(to='db.Reagent', null=True),
        ),
        migrations.AddField(
            model_name='attachedfile',
            name='screen',
            field=models.ForeignKey(to='db.Screen', null=True),
        ),
        migrations.AddField(
            model_name='attachedfile',
            name='screensaver_user',
            field=models.ForeignKey(to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='attachedfile',
            name='publication',
            field=models.OneToOneField(to='db.Publication', null=True),
        ),
        migrations.AddField(
            model_name='attachedfile',
            name='attached_file_type',
            field=models.ForeignKey(blank=True, to='db.AttachedFileType', null=False),
        ),
        migrations.AddField(
            model_name='assaywell',
            name='screen_result',
            field=models.ForeignKey(to='db.ScreenResult'),
        ),
        migrations.AddField(
            model_name='assaywell',
            name='well',
            field=models.ForeignKey(to='db.Well'),
        ),
        migrations.AddField(
            model_name='assayplate',
            name='plate',
            field=models.ForeignKey(blank=True, to='db.Plate', null=True),
        ),
        migrations.AddField(
            model_name='assayplate',
            name='screen',
            field=models.ForeignKey(to='db.Screen'),
        ),
        migrations.AddField(
            model_name='annotationvalue',
            name='reagent',
            field=models.ForeignKey(blank=True, to='db.Reagent', null=True),
        ),
        migrations.AddField(
            model_name='annotationtype',
            name='study',
            field=models.ForeignKey(to='db.Screen'),
        ),
        migrations.AddField(
            model_name='activityupdateactivity',
            name='activity',
            field=models.ForeignKey(to='db.Activity'),
        ),
        migrations.AddField(
            model_name='activity',
            name='created_by',
            field=models.ForeignKey(related_name='activities_created', blank=True, to='db.ScreensaverUser', null=True),
        ),
        migrations.AddField(
            model_name='activity',
            name='performed_by',
            field=models.ForeignKey(related_name='activities_performed', to='db.ScreensaverUser'),
        ),
        migrations.AddField(
            model_name='abasetestset',
            name='screen',
            field=models.ForeignKey(to='db.Screen'),
        ),
        migrations.CreateModel(
            name='CherryPickLiquidTransfer',
            fields=[
                ('status', models.TextField()),
                ('labactivitylink', models.OneToOneField(
                    primary_key=True, parent_link=True,db_column='activity_id',
                    serialize=False, to='db.LabActivity')),
            ],
            options={
                'db_table': 'cherry_pick_liquid_transfer',
            },
        ),
        migrations.CreateModel(
            name='LabHead',
            fields=[
                ('screensaver_user', models.OneToOneField(primary_key=True, serialize=False, to='db.ScreeningRoomUser')),
                ('lab_affiliation', models.ForeignKey(blank=True, 
                    to='db.LabAffiliation', null=True)),
            ],
            options={
                'db_table': 'lab_head',
            },
        ),
        migrations.CreateModel(
            name='Screening',
            fields=[
                ('labactivitylink', models.OneToOneField(
                    primary_key=True, serialize=False, parent_link=True,
                    db_column='activity_id', to='db.LabActivity')),
                ('assay_protocol', models.TextField(blank=True)),
                ('number_of_replicates', models.IntegerField(null=True, blank=True)),
                ('assay_protocol_type', models.TextField(blank=True)),
                ('assay_protocol_last_modified_date', models.DateField(null=True, blank=True)),
                ('assay_well_volume', models.DecimalField(null=True, max_digits=10, decimal_places=9, blank=True)),
                ('volume_transferred_per_well_to_assay_plates', models.DecimalField(null=True, max_digits=10, decimal_places=9, blank=True)),
            ],
            options={
                'db_table': 'screening',
            },
        ),
        migrations.CreateModel(
            name='WellVolumeCorrectionActivity',
            fields=[
                ('activity', models.OneToOneField(primary_key=True, serialize=False, to='db.AdministrativeActivity')),
            ],
            options={
                'db_table': 'well_volume_correction_activity',
            },
        ),
        migrations.AddField(
            model_name='well',
            name='deprecation_admin_activity',
            field=models.ForeignKey(blank=True, to='db.AdministrativeActivity', null=True),
        ),
#         migrations.AlterUniqueTogether(
#             name='userchecklistitem',
#             unique_together=set([('screensaver_user', 'item_group', 'item_name')]),
#         ),
        migrations.AddField(
            model_name='silencingreagent',
            name='duplex_wells',
            field=models.ManyToManyField(to='db.Well'),
        ),
#         migrations.AddField(
#             model_name='silencingreagent',
#             name='facility_gene',
#             field=models.ForeignKey(related_name='facility_reagent', null=True, blank=True, to='db.Gene', unique=True),
#         ),
#         migrations.AddField(
#             model_name='silencingreagent',
#             name='vendor_gene',
#             field=models.ForeignKey(related_name='vendor_reagent', null=True, blank=True, to='db.Gene', unique=True),
#         ),
        migrations.AddField(
            model_name='serviceactivity',
            name='serviced_screen',
            field=models.ForeignKey(blank=True, to='db.Screen', null=True),
        ),
        migrations.AddField(
            model_name='serviceactivity',
            name='serviced_user',
            field=models.ForeignKey(to='db.ScreensaverUser'),
        ),
        migrations.AddField(
            model_name='screenupdateactivity',
            name='update_activity',
            field=models.ForeignKey(to='db.AdministrativeActivity'),
        ),
        migrations.AlterIndexTogether(
            name='screenstatusitem',
            index_together=set([('screen', 'status', 'status_date')]),
        ),
        migrations.AddField(
            model_name='screensaveruserupdateactivity',
            name='update_activity',
            field=models.ForeignKey(to='db.AdministrativeActivity'),
        ),
        migrations.AddField(
            model_name='screenresultupdateactivity',
            name='update_activity',
            field=models.ForeignKey(to='db.AdministrativeActivity'),
        ),
        migrations.AddField(
            model_name='screeningroomuserfacilityusagerole',
            name='screening_room_user',
            field=models.ForeignKey(to='db.ScreeningRoomUser'),
        ),
        migrations.AddField(
            model_name='screeningroomuser',
            name='last_notified_rnaiua_checklist_item_event',
            field=models.ForeignKey(related_name='rnai_ua_user', blank=True, to='db.ChecklistItemEvent', null=True),
        ),
        migrations.AddField(
            model_name='screeningroomuser',
            name='last_notified_smua_checklist_item_event',
            field=models.ForeignKey(related_name='smua_user', blank=True, to='db.ChecklistItemEvent', null=True),
        ),
#         migrations.AlterUniqueTogether(
#             name='screenfundingsupports',
#             unique_together=set([('screen', 'funding_support')]),
#         ),
        migrations.AddField(
            model_name='screen',
            name='lead_screener',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', 
                null=True, related_name='led_screen'),
        ),
        migrations.AddField(
            model_name='screen',
            name='pin_transfer_admin_activity',
            field=models.ForeignKey(blank=True, to='db.AdministrativeActivity', null=True),
        ),
        migrations.AddField(
            model_name='libraryupdateactivity',
            name='update_activity',
            field=models.ForeignKey(to='db.AdministrativeActivity'),
        ),
        migrations.AddField(
            model_name='librarycontentsversion',
            name='library_contents_loading_activity',
            field=models.ForeignKey(related_name='lcv_load', to='db.AdministrativeActivity'),
        ),
        migrations.AddField(
            model_name='librarycontentsversion',
            name='library_contents_release_activity',
            field=models.ForeignKey(related_name='lcv_release', blank=True, to='db.AdministrativeActivity', null=True),
        ),
        migrations.AddField(
            model_name='library',
            name='owner_screener',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', 
                null=True, related_name='owned_library'),
        ),
        migrations.AddField(
            model_name='labactivity',
            name='screen',
            field=models.ForeignKey(to='db.Screen'),
        ),
        migrations.AlterUniqueTogether(
            name='genesymbol',
            unique_together=set([('gene', 'ordinal')]),
        ),
        migrations.AlterUniqueTogether(
            name='genegenbankaccessionnumber',
            unique_together=set([('gene', 'genbank_accession_number')]),
        ),
        migrations.AlterUniqueTogether(
            name='screenkeyword',
            unique_together=set([('screen', 'keyword')]),
        ),
        
        migrations.AddField(
            model_name='equipmentused',
            name='lab_activity',
            field=models.ForeignKey(to='db.LabActivity'),
        ),
#         migrations.AddField(
#             model_name='collaboratorlink',
#             name='collaborator',
#             field=models.ForeignKey(to='db.ScreeningRoomUser'),
#         ),
        
        # New screen_collaborators table, see manual/0002 sql for creation
        migrations.AddField(
            model_name='screen',
            name='collaborators',
            field=models.ManyToManyField(to='db.ScreensaverUser',
                related_name='collaborating_screens'),
        ),
        
        migrations.AddField(
            model_name='cherrypickrequestupdateactivity',
            name='update_activity',
            field=models.OneToOneField(to='db.AdministrativeActivity', unique=True),
        ),
        migrations.AddField(
            model_name='cherrypickrequest',
            name='requested_by',
            field=models.ForeignKey(to='db.ScreeningRoomUser',
                                    related_name='requested_cherry_pick'),
        ),
        migrations.AddField(
            model_name='cherrypickrequest',
            name='volume_approved_by',
            field=models.ForeignKey(blank=True, to='db.AdministratorUser', 
                null=True, related_name='approved_cherry_pick'),
        ),
        migrations.AddField(
            model_name='checklistitemeventupdateactivity',
            name='update_activity',
            field=models.ForeignKey(to='db.AdministrativeActivity'),
        ),
        migrations.AddField(
            model_name='checklistitemevent',
            name='screening_room_user',
            field=models.ForeignKey(to='db.ScreeningRoomUser'),
        ),
        migrations.AddField(
            model_name='attachedfileupdateactivity',
            name='update_activity',
            field=models.ForeignKey(to='db.AdministrativeActivity'),
        ),
        migrations.AddField(
            model_name='assayplate',
            name='screen_result_data_loading',
            field=models.ForeignKey(blank=True, to='db.AdministrativeActivity', null=True),
        ),
        migrations.AddField(
            model_name='activityupdateactivity',
            name='update_activity',
            field=models.ForeignKey(to='db.AdministrativeActivity'),
        ),
        migrations.CreateModel(
            name='CherryPickScreening',
            fields=[
                ('screeninglink', models.OneToOneField(
                    primary_key=True, parent_link=True, db_column='activity_id', 
                    serialize=False, to='db.Screening')),
                ('cherry_pick_request', models.ForeignKey(to='db.CherryPickRequest')),
            ],
            options={
                'db_table': 'cherry_pick_screening',
            },
        ),
        migrations.CreateModel(
            name='LibraryScreening',
            fields=[
                ('screeninglink', models.OneToOneField(
                    primary_key=True, parent_link=True, db_column='activity_id',
                    serialize=False, to='db.Screening')),
                ('abase_testset_id', models.TextField(blank=True)),
                ('is_for_external_library_plates', models.BooleanField()),
                ('screened_experimental_well_count', models.IntegerField(default=0)),
                ('libraries_screened_count', models.IntegerField(null=True, blank=True)),
                ('library_plates_screened_count', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'db_table': 'library_screening',
            },
        ),
        migrations.AddField(
            model_name='wellvolumeadjustment',
            name='well_volume_correction_activity',
            field=models.ForeignKey(blank=True, to='db.WellVolumeCorrectionActivity', null=True),
        ),
        migrations.AddField(
            model_name='screeningroomuser',
            name='lab_head',
            field=models.ForeignKey(blank=True, to='db.LabHead', null=True),
        ),
        migrations.AddField(
            model_name='screen',
            name='lab_head',
            field=models.ForeignKey(blank=True, to='db.ScreensaverUser', 
                related_name='lab_head_screen', null=True),
        ),
        migrations.AddField(
            model_name='cherrypickassayplate',
            name='cherry_pick_liquid_transfer',
            field=models.ForeignKey(blank=True, to='db.CherryPickLiquidTransfer', null=True),
        ),
        migrations.AddField(
            model_name='cherrypickassayplatescreeninglink',
            name='cherry_pick_screening',
            field=models.ForeignKey(to='db.CherryPickScreening'),
        ),
        migrations.AlterUniqueTogether(
            name='cherrypickassayplate',
            unique_together=set([('cherry_pick_request', 'plate_ordinal', 'attempt_ordinal')]),
        ),
        migrations.AddField(
            model_name='assayplate',
            name='library_screening',
            field=models.ForeignKey(blank=True, to='db.LibraryScreening', null=True),
        ),
    ]
