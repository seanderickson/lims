# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'LegacySmallMoleculeCasNumber'
        db.create_table(u'_legacy_small_molecule_cas_number', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('smiles', self.gf('django.db.models.fields.CharField')(max_length=2047)),
            ('cas_number', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['LegacySmallMoleculeCasNumber'])

        # Adding model 'AbaseTestset'
        db.create_table(u'abase_testset', (
            ('abase_testset_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('testset_name', self.gf('django.db.models.fields.TextField')()),
            ('comments', self.gf('django.db.models.fields.TextField')()),
            ('testset_date', self.gf('django.db.models.fields.DateField')()),
        ))
        db.send_create_signal(u'db', ['AbaseTestset'])

        # Adding model 'Activity'
        db.create_table(u'activity', (
            ('activity_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('performed_by', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'activities_performed', to=orm['db.ScreensaverUser'])),
            ('date_of_activity', self.gf('django.db.models.fields.DateField')()),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'activities_created', null=True, to=orm['db.ScreensaverUser'])),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['Activity'])

        # Adding model 'ActivityUpdateActivity'
        db.create_table(u'activity_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Activity'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['ActivityUpdateActivity'])

        # Adding model 'AdministrativeActivity'
        db.create_table(u'administrative_activity', (
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Activity'], primary_key=True)),
            ('administrative_activity_type', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['AdministrativeActivity'])

        # Adding model 'AdministratorUser'
        db.create_table(u'administrator_user', (
            ('screensaver_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['AdministratorUser'])

        # Adding model 'AnnotationType'
        db.create_table(u'annotation_type', (
            ('annotation_type_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('study', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('name', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('is_numeric', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'db', ['AnnotationType'])

        # Adding model 'AnnotationValue'
        db.create_table(u'annotation_value', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('annotation_value_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('numeric_value', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('value', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('annotation_type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AnnotationType'], null=True, blank=True)),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Reagent'], null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['AnnotationValue'])

        # Adding model 'AssayPlate'
        db.create_table(u'assay_plate', (
            ('assay_plate_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('replicate_ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('library_screening', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LibraryScreening'], null=True, blank=True)),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('plate', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Plate'], null=True, blank=True)),
            ('plate_number', self.gf('django.db.models.fields.IntegerField')()),
            ('screen_result_data_loading', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'], null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['AssayPlate'])

        # Adding model 'AssayWell'
        db.create_table(u'assay_well', (
            ('assay_well_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('assay_well_control_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('is_positive', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('screen_result', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreenResult'])),
            ('well', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'])),
            ('confirmed_positive_value', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal(u'db', ['AssayWell'])

        # Adding model 'AttachedFile'
        db.create_table(u'attached_file', (
            ('attached_file_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('file_contents', self.gf('django.db.models.fields.TextField')()),
            ('filename', self.gf('django.db.models.fields.TextField')()),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'], null=True, blank=True)),
            ('screensaver_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'], null=True, blank=True)),
            ('attached_file_type', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AttachedFileType'])),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Reagent'], null=True, blank=True)),
            ('file_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['AttachedFile'])

        # Adding model 'AttachedFileType'
        db.create_table(u'attached_file_type', (
            ('for_entity_type', self.gf('django.db.models.fields.CharField')(max_length=31)),
            ('attached_file_type_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['AttachedFileType'])

        # Adding model 'AttachedFileUpdateActivity'
        db.create_table(u'attached_file_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('attached_file', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AttachedFile'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['AttachedFileUpdateActivity'])

        # Adding model 'Cell'
        db.create_table(u'cell', (
            ('cell_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('alternate_id', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('alternate_name', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('batch_id', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('cell_type', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('cell_type_detail', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('center_name', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('center_specific_id', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('clo_id', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('disease', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('disease_detail', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('facility_id', self.gf('django.db.models.fields.CharField')(unique=True, max_length=255)),
            ('genetic_modification', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('mutations_explicit', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('mutations_reference', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('organ', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('organism', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('organism_gender', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('recommended_culture_conditions', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('tissue', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('vendor', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('vendor_catalog_id', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('verification', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('verification_reference_profile', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['Cell'])

        # Adding model 'CellGrowthProperties'
        db.create_table(u'cell_growth_properties', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Cell'])),
            ('growth_property', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['CellGrowthProperties'])

        # Adding model 'CellLine'
        db.create_table(u'cell_line', (
            ('cell_line_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.TextField')(unique=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['CellLine'])

        # Adding model 'CellLineage'
        db.create_table(u'cell_lineage', (
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Cell'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['CellLineage'])

        # Adding model 'CellMarkers'
        db.create_table(u'cell_markers', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.PrimaryCell'])),
            ('cell_markers', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['CellMarkers'])

        # Adding model 'CellRelatedProjects'
        db.create_table(u'cell_related_projects', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Cell'])),
            ('related_project', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['CellRelatedProjects'])

        # Adding model 'CellUpdateActivity'
        db.create_table(u'cell_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Cell'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'], unique=True)),
        ))
        db.send_create_signal(u'db', ['CellUpdateActivity'])

        # Adding model 'ChecklistItem'
        db.create_table(u'checklist_item', (
            ('checklist_item_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('checklist_item_group', self.gf('django.db.models.fields.TextField')()),
            ('is_expirable', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('item_name', self.gf('django.db.models.fields.TextField')(unique=True)),
            ('order_statistic', self.gf('django.db.models.fields.IntegerField')()),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['ChecklistItem'])

        # Adding model 'ChecklistItemEvent'
        db.create_table(u'checklist_item_event', (
            ('checklist_item_event_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('date_performed', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('is_expiration', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('checklist_item_id', self.gf('django.db.models.fields.IntegerField')()),
            ('screening_room_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'])),
            ('is_not_applicable', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['ChecklistItemEvent'])

        # Adding model 'ChecklistItemEventUpdateActivity'
        db.create_table(u'checklist_item_event_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('checklist_item_event', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ChecklistItemEvent'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['ChecklistItemEventUpdateActivity'])

        # Adding model 'CherryPickAssayPlate'
        db.create_table(u'cherry_pick_assay_plate', (
            ('cherry_pick_assay_plate_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'])),
            ('plate_ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('attempt_ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('cherry_pick_liquid_transfer', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickLiquidTransfer'], null=True, blank=True)),
            ('assay_plate_type', self.gf('django.db.models.fields.TextField')()),
            ('legacy_plate_name', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cherry_pick_assay_plate_type', self.gf('django.db.models.fields.CharField')(max_length=31)),
        ))
        db.send_create_signal(u'db', ['CherryPickAssayPlate'])

        # Adding model 'CherryPickAssayPlateScreeningLink'
        db.create_table(u'cherry_pick_assay_plate_screening_link', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cherry_pick_assay_plate', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickAssayPlate'])),
            ('cherry_pick_screening', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickScreening'])),
        ))
        db.send_create_signal(u'db', ['CherryPickAssayPlateScreeningLink'])

        # Adding model 'CherryPickLiquidTransfer'
        db.create_table(u'cherry_pick_liquid_transfer', (
            ('status', self.gf('django.db.models.fields.TextField')()),
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LabActivity'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['CherryPickLiquidTransfer'])

        # Adding model 'CherryPickRequest'
        db.create_table(u'cherry_pick_request', (
            ('cherry_pick_request_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('requested_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'])),
            ('is_randomized_assay_plate_layout', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('legacy_cherry_pick_request_number', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('volume_approved_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministratorUser'], null=True, blank=True)),
            ('number_unfulfilled_lab_cherry_picks', self.gf('django.db.models.fields.IntegerField')()),
            ('assay_plate_type', self.gf('django.db.models.fields.TextField')()),
            ('transfer_volume_per_well_approved', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9, blank=True)),
            ('transfer_volume_per_well_requested', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9, blank=True)),
            ('date_requested', self.gf('django.db.models.fields.DateField')()),
            ('date_volume_approved', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('assay_protocol_comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cherry_pick_assay_protocols_followed', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cherry_pick_followup_results_status', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('keep_source_plate_cherry_picks_together', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('max_skipped_wells_per_plate', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['CherryPickRequest'])

        # Adding model 'CherryPickRequestEmptyWell'
        db.create_table(u'cherry_pick_request_empty_well', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'])),
            ('well_name', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
        ))
        db.send_create_signal(u'db', ['CherryPickRequestEmptyWell'])

        # Adding model 'CherryPickRequestUpdateActivity'
        db.create_table(u'cherry_pick_request_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'], unique=True)),
        ))
        db.send_create_signal(u'db', ['CherryPickRequestUpdateActivity'])

        # Adding model 'CherryPickScreening'
        db.create_table(u'cherry_pick_screening', (
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screening'], primary_key=True)),
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'])),
        ))
        db.send_create_signal(u'db', ['CherryPickScreening'])

        # Adding model 'CollaboratorLink'
        db.create_table(u'collaborator_link', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('collaborator', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'])),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
        ))
        db.send_create_signal(u'db', ['CollaboratorLink'])

        # Adding model 'Copy'
        db.create_table(u'copy', (
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('usage_type', self.gf('django.db.models.fields.TextField')()),
            ('library', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Library'])),
            ('name', self.gf('django.db.models.fields.TextField')()),
            ('copy_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('date_plated', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('primary_plate_status', self.gf('django.db.models.fields.TextField')()),
            ('primary_plate_location_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('plates_available', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('plate_locations_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('primary_well_mg_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
            ('primary_well_molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('min_molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('max_molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('min_mg_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
            ('max_mg_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
            ('well_concentration_dilution_factor', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=8, decimal_places=2, blank=True)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['Copy'])

        # Adding model 'CopyUpdateActivity'
        db.create_table(u'copy_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('copy_id', self.gf('django.db.models.fields.IntegerField')()),
            ('update_activity_id', self.gf('django.db.models.fields.IntegerField')(unique=True)),
        ))
        db.send_create_signal(u'db', ['CopyUpdateActivity'])

        # Adding model 'DataColumn'
        db.create_table(u'data_column', (
            ('data_column_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('replicate_ordinal', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('assay_phenotype', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('assay_readout_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('how_derived', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('is_follow_up_data', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('name', self.gf('django.db.models.fields.TextField')()),
            ('time_point', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('is_derived', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('positives_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('screen_result', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreenResult'])),
            ('channel', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('time_point_ordinal', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('zdepth_ordinal', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('data_type', self.gf('django.db.models.fields.TextField')()),
            ('decimal_places', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('strong_positives_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('medium_positives_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('weak_positives_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['DataColumn'])

        # Adding model 'DataColumnDerivedFromLink'
        db.create_table(u'data_column_derived_from_link', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('derived_data_column', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.DataColumn'])),
            ('derived_from_data_column', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'derived_from', to=orm['db.DataColumn'])),
        ))
        db.send_create_signal(u'db', ['DataColumnDerivedFromLink'])

        # Adding model 'EquipmentUsed'
        db.create_table(u'equipment_used', (
            ('equipment_used_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('protocol', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('equipment', self.gf('django.db.models.fields.TextField')()),
            ('lab_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LabActivity'])),
        ))
        db.send_create_signal(u'db', ['EquipmentUsed'])

        # Adding model 'ExperimentalCellInformation'
        db.create_table(u'experimental_cell_information', (
            ('experimental_cell_information_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Cell'])),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
        ))
        db.send_create_signal(u'db', ['ExperimentalCellInformation'])

        # Adding model 'FundingSupport'
        db.create_table(u'funding_support', (
            ('funding_support_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.TextField')(unique=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['FundingSupport'])

        # Adding model 'Gene'
        db.create_table(u'gene', (
            ('gene_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('entrezgene_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('gene_name', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('species_name', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal(u'db', ['Gene'])

        # Adding model 'GeneGenbankAccessionNumber'
        db.create_table(u'gene_genbank_accession_number', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('gene_id', self.gf('django.db.models.fields.IntegerField')()),
            ('genbank_accession_number', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['GeneGenbankAccessionNumber'])

        # Adding model 'GeneOldEntrezgeneId'
        db.create_table(u'gene_old_entrezgene_id', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('old_entrezgene_id', self.gf('django.db.models.fields.IntegerField')()),
            ('gene_id', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['GeneOldEntrezgeneId'])

        # Adding model 'GeneOldEntrezgeneSymbol'
        db.create_table(u'gene_old_entrezgene_symbol', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('old_entrezgene_symbol', self.gf('django.db.models.fields.TextField')()),
            ('gene_id', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['GeneOldEntrezgeneSymbol'])

        # Adding model 'GeneSymbol'
        db.create_table(u'gene_symbol', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('gene', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Gene'])),
            ('entrezgene_symbol', self.gf('django.db.models.fields.TextField')()),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['GeneSymbol'])

        # Adding model 'LabActivity'
        db.create_table(u'lab_activity', (
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Activity'], primary_key=True)),
            ('volume_transferred_per_well_from_library_plates', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9, blank=True)),
            ('molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
        ))
        db.send_create_signal(u'db', ['LabActivity'])

        # Adding model 'LabAffiliation'
        db.create_table(u'lab_affiliation', (
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('affiliation_name', self.gf('django.db.models.fields.TextField')(unique=True)),
            ('affiliation_category', self.gf('django.db.models.fields.TextField')()),
            ('lab_affiliation_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
        ))
        db.send_create_signal(u'db', ['LabAffiliation'])

        # Adding model 'LabCherryPick'
        db.create_table(u'lab_cherry_pick', (
            ('lab_cherry_pick_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'])),
            ('screener_cherry_pick', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreenerCherryPick'])),
            ('source_well', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'])),
            ('cherry_pick_assay_plate', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickAssayPlate'], null=True, blank=True)),
            ('assay_plate_row', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('assay_plate_column', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['LabCherryPick'])

        # Adding model 'Library'
        db.create_table(u'library', (
            ('library_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('library_name', self.gf('django.db.models.fields.TextField')(unique=True)),
            ('short_name', self.gf('django.db.models.fields.TextField')(unique=True)),
            ('description', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('provider', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('screen_type', self.gf('django.db.models.fields.TextField')()),
            ('library_type', self.gf('django.db.models.fields.TextField')()),
            ('start_plate', self.gf('django.db.models.fields.IntegerField')(unique=True)),
            ('end_plate', self.gf('django.db.models.fields.IntegerField')(unique=True)),
            ('screening_status', self.gf('django.db.models.fields.TextField')()),
            ('date_received', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('date_screenable', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('plate_size', self.gf('django.db.models.fields.TextField')()),
            ('latest_released_contents_version_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('experimental_well_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('is_pool', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('owner_screener', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'], null=True, blank=True)),
            ('solvent', self.gf('django.db.models.fields.TextField')()),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['Library'])

        # Adding model 'LibraryContentsVersion'
        db.create_table(u'library_contents_version', (
            ('library_contents_version_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('version_number', self.gf('django.db.models.fields.IntegerField')()),
            ('library', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Library'])),
            ('library_contents_loading_activity', self.gf('django.db.models.fields.related.ForeignKey')(related_name=u'lcv_load', to=orm['db.AdministrativeActivity'])),
            ('library_contents_release_activity', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'lcv_release', null=True, to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['LibraryContentsVersion'])

        # Adding model 'LibraryScreening'
        db.create_table(u'library_screening', (
            ('abase_testset_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('is_for_external_library_plates', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screening'], primary_key=True)),
            ('screened_experimental_well_count', self.gf('django.db.models.fields.IntegerField')()),
            ('libraries_screened_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('library_plates_screened_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['LibraryScreening'])

        # Adding model 'LibraryUpdateActivity'
        db.create_table(u'library_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('library', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Library'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['LibraryUpdateActivity'])

        # Adding model 'Molfile'
        db.create_table(u'molfile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('molfile', self.gf('django.db.models.fields.TextField')()),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.SmallMoleculeReagent'], unique=True)),
        ))
        db.send_create_signal(u'db', ['Molfile'])

        # Adding model 'NaturalProductReagent'
        db.create_table(u'natural_product_reagent', (
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Reagent'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['NaturalProductReagent'])

        # Adding model 'Plate'
        db.create_table(u'plate', (
            ('plate_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('plate_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('plate_number', self.gf('django.db.models.fields.IntegerField')()),
            ('well_volume', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9, blank=True)),
            ('copy', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Copy'])),
            ('facility_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('plate_location', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.PlateLocation'], null=True, blank=True)),
            ('status', self.gf('django.db.models.fields.TextField')()),
            ('retired_activity_id', self.gf('django.db.models.fields.IntegerField')(unique=True, null=True, blank=True)),
            ('plated_activity_id', self.gf('django.db.models.fields.IntegerField')(unique=True, null=True, blank=True)),
            ('stock_plate_number', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('quadrant', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('min_molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('max_molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('min_mg_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
            ('max_mg_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
            ('primary_well_molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('primary_well_mg_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['Plate'])

        # Adding model 'PlateLocation'
        db.create_table(u'plate_location', (
            ('plate_location_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('bin', self.gf('django.db.models.fields.TextField')()),
            ('freezer', self.gf('django.db.models.fields.TextField')()),
            ('room', self.gf('django.db.models.fields.TextField')()),
            ('shelf', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['PlateLocation'])

        # Adding model 'PlateUpdateActivity'
        db.create_table(u'plate_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('plate_id', self.gf('django.db.models.fields.IntegerField')()),
            ('update_activity_id', self.gf('django.db.models.fields.IntegerField')(unique=True)),
        ))
        db.send_create_signal(u'db', ['PlateUpdateActivity'])

        # Adding model 'PrimaryCell'
        db.create_table(u'primary_cell', (
            ('age_in_years', self.gf('django.db.models.fields.IntegerField')()),
            ('donor_ethnicity', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('donor_health_status', self.gf('django.db.models.fields.CharField')(max_length=255, blank=True)),
            ('passage_number', self.gf('django.db.models.fields.IntegerField')()),
            ('cell', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Cell'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['PrimaryCell'])

        # Adding model 'Publication'
        db.create_table(u'publication', (
            ('publication_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('authors', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('journal', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('pages', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('pubmed_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('title', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('volume', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('year_published', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('attached_file', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AttachedFile'], unique=True, null=True, blank=True)),
            ('pubmed_central_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['Publication'])

        # Adding model 'Reagent'
        db.create_table(u'reagent', (
            ('reagent_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('vendor_identifier', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('vendor_name', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('library_contents_version', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LibraryContentsVersion'])),
            ('well', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'])),
            ('facility_batch_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('vendor_batch_id', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal(u'db', ['Reagent'])

        # Adding model 'ReagentPublicationLink'
        db.create_table(u'reagent_publication_link', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Reagent'])),
            ('publication_id', self.gf('django.db.models.fields.IntegerField')(unique=True)),
        ))
        db.send_create_signal(u'db', ['ReagentPublicationLink'])

        # Adding model 'ResultValue'
        db.create_table(u'result_value', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('result_value_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('assay_well_control_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('is_exclude', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('numeric_value', self.gf('django.db.models.fields.FloatField')(null=True, blank=True)),
            ('is_positive', self.gf('django.db.models.fields.NullBooleanField')(null=True, blank=True)),
            ('value', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('data_column', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.DataColumn'], null=True, blank=True)),
            ('well', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'], null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['ResultValue'])

        # Adding model 'RnaiCherryPickRequest'
        db.create_table(u'rnai_cherry_pick_request', (
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['RnaiCherryPickRequest'])

        # Adding model 'SchemaHistory'
        db.create_table(u'schema_history', (
            ('screensaver_revision', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('date_updated', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('comment', self.gf('django.db.models.fields.TextField')(blank=True)),
        ))
        db.send_create_signal(u'db', ['SchemaHistory'])

        # Adding model 'Screen'
        db.create_table(u'screen', (
            ('screen_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('screen_type', self.gf('django.db.models.fields.TextField')()),
            ('title', self.gf('django.db.models.fields.TextField')()),
            ('summary', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('abase_study_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('abase_protocol_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('publishable_protocol', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('lead_screener', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'], null=True, blank=True)),
            ('lab_head', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LabHead'], null=True, blank=True)),
            ('publishable_protocol_entered_by', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('publishable_protocol_comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('study_type', self.gf('django.db.models.fields.TextField')()),
            ('url', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('data_meeting_complete', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('data_meeting_scheduled', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('date_of_application', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('publishable_protocol_date_entered', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('amount_to_be_charged_for_screen', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=9, decimal_places=2, blank=True)),
            ('billing_comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('is_billing_for_supplies_only', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('billing_info_return_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('date_charged', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('date_completed5kcompounds', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('date_faxed_to_billing_department', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('facilities_and_administration_charge', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=9, decimal_places=2, blank=True)),
            ('is_fee_form_on_file', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('fee_form_requested_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('fee_form_requested_initials', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('see_comments', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('to_be_requested', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('coms_registration_number', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('coms_approval_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('pin_transfer_admin_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'], null=True, blank=True)),
            ('data_sharing_level', self.gf('django.db.models.fields.IntegerField')()),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('data_privacy_expiration_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('max_allowed_data_privacy_expiration_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('min_allowed_data_privacy_expiration_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('data_privacy_expiration_notified_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('screened_experimental_well_count', self.gf('django.db.models.fields.IntegerField')()),
            ('unique_screened_experimental_well_count', self.gf('django.db.models.fields.IntegerField')()),
            ('total_plated_lab_cherry_picks', self.gf('django.db.models.fields.IntegerField')()),
            ('assay_plates_screened_count', self.gf('django.db.models.fields.IntegerField')()),
            ('library_plates_screened_count', self.gf('django.db.models.fields.IntegerField')()),
            ('library_plates_data_loaded_count', self.gf('django.db.models.fields.IntegerField')()),
            ('library_plates_data_analyzed_count', self.gf('django.db.models.fields.IntegerField')()),
            ('min_screened_replicate_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('max_screened_replicate_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('min_data_loaded_replicate_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('max_data_loaded_replicate_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('libraries_screened_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('facility_id', self.gf('django.db.models.fields.TextField')(unique=True)),
            ('project_phase', self.gf('django.db.models.fields.TextField')()),
            ('project_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('pubchem_assay_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('pubchem_deposited_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('image_url', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('well_studied', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'], null=True, blank=True)),
            ('species', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('cell_line', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CellLine'], null=True, blank=True)),
            ('transfection_agent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.TransfectionAgent'], null=True, blank=True)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('perturbagen_molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('perturbagen_ug_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
        ))
        db.send_create_signal(u'db', ['Screen'])

        # Adding model 'ScreenBillingItem'
        db.create_table(u'screen_billing_item', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('amount', self.gf('django.db.models.fields.DecimalField')(max_digits=9, decimal_places=2)),
            ('date_sent_for_billing', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('item_to_be_charged', self.gf('django.db.models.fields.TextField')()),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['ScreenBillingItem'])

        # Adding model 'ScreenFundingSupportLink'
        db.create_table(u'screen_funding_support_link', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('funding_support', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.FundingSupport'])),
        ))
        db.send_create_signal(u'db', ['ScreenFundingSupportLink'])

        # Adding model 'ScreenKeyword'
        db.create_table(u'screen_keyword', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('keyword', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['ScreenKeyword'])

        # Adding model 'ScreenPublicationLink'
        db.create_table(u'screen_publication_link', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('publication_id', self.gf('django.db.models.fields.IntegerField')(unique=True)),
        ))
        db.send_create_signal(u'db', ['ScreenPublicationLink'])

        # Adding model 'ScreenResult'
        db.create_table(u'screen_result', (
            ('screen_result_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('replicate_count', self.gf('django.db.models.fields.IntegerField')()),
            ('experimental_well_count', self.gf('django.db.models.fields.IntegerField')()),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('channel_count', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['ScreenResult'])

        # Adding model 'ScreenResultUpdateActivity'
        db.create_table(u'screen_result_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screen_result', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreenResult'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['ScreenResultUpdateActivity'])

        # Adding model 'ScreenStatusItem'
        db.create_table(u'screen_status_item', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('status', self.gf('django.db.models.fields.TextField')()),
            ('status_date', self.gf('django.db.models.fields.DateField')()),
        ))
        db.send_create_signal(u'db', ['ScreenStatusItem'])

        # Adding index on 'ScreenStatusItem', fields ['screen', 'status', 'status_date']
        db.create_index(u'screen_status_item', ['screen_id', 'status', 'status_date'])

        # Adding model 'ScreenUpdateActivity'
        db.create_table(u'screen_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['ScreenUpdateActivity'])

        # Adding model 'ScreenerCherryPick'
        db.create_table(u'screener_cherry_pick', (
            ('screener_cherry_pick_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'])),
            ('screened_well', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'])),
        ))
        db.send_create_signal(u'db', ['ScreenerCherryPick'])

        # Adding model 'Screening'
        db.create_table(u'screening', (
            ('assay_protocol', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('number_of_replicates', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('assay_protocol_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LabActivity'], primary_key=True)),
            ('assay_protocol_last_modified_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('assay_well_volume', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9, blank=True)),
            ('volume_transferred_per_well_to_assay_plates', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9, blank=True)),
        ))
        db.send_create_signal(u'db', ['Screening'])

        # Adding model 'ScreeningRoomUser'
        db.create_table(u'screening_room_user', (
            ('screensaver_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], primary_key=True)),
            ('user_classification', self.gf('django.db.models.fields.TextField')()),
            ('lab_head', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LabHead'], null=True, blank=True)),
            ('coms_crhba_permit_number', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('coms_crhba_permit_principal_investigator', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('last_notified_smua_checklist_item_event', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'smua_user', null=True, to=orm['db.ChecklistItemEvent'])),
            ('last_notified_rnaiua_checklist_item_event', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'rnai_ua_user', null=True, to=orm['db.ChecklistItemEvent'])),
        ))
        db.send_create_signal(u'db', ['ScreeningRoomUser'])

        # Adding model 'LabHead'
        db.create_table(u'lab_head', (
            ('screensaver_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'], primary_key=True)),
            ('lab_affiliation', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LabAffiliation'], null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['LabHead'])

        # Adding model 'ScreeningRoomUserFacilityUsageRole'
        db.create_table(u'screening_room_user_facility_usage_role', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screening_room_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'])),
            ('facility_usage_role', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['ScreeningRoomUserFacilityUsageRole'])

        # Adding model 'ScreensaverUser'
        db.create_table(u'screensaver_user', (
            ('screensaver_user_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')(blank=True)),
            ('date_created', self.gf('django.db.models.fields.DateTimeField')()),
            ('first_name', self.gf('django.db.models.fields.TextField')()),
            ('last_name', self.gf('django.db.models.fields.TextField')()),
            ('email', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('phone', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('mailing_address', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('comments', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('login_id', self.gf('django.db.models.fields.TextField')(unique=True, blank=True)),
            ('digested_password', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('ecommons_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('harvard_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('harvard_id_expiration_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('harvard_id_requested_expiration_date', self.gf('django.db.models.fields.DateField')(null=True, blank=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'], null=True, blank=True)),
            ('date_loaded', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
            ('date_publicly_available', self.gf('django.db.models.fields.DateTimeField')(null=True, blank=True)),
        ))
        db.send_create_signal(u'db', ['ScreensaverUser'])

        # Adding model 'ScreensaverUserRole'
        db.create_table(u'screensaver_user_role', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screensaver_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'])),
            ('screensaver_user_role', self.gf('django.db.models.fields.TextField')()),
        ))
        db.send_create_signal(u'db', ['ScreensaverUserRole'])

        # Adding model 'ScreensaverUserUpdateActivity'
        db.create_table(u'screensaver_user_update_activity', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('screensaver_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreensaverUser'])),
            ('update_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'])),
        ))
        db.send_create_signal(u'db', ['ScreensaverUserUpdateActivity'])

        # Adding model 'ServiceActivity'
        db.create_table(u'service_activity', (
            ('service_activity_type', self.gf('django.db.models.fields.TextField')()),
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Activity'], primary_key=True)),
            ('serviced_screen', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'], null=True, blank=True)),
            ('serviced_user', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.ScreeningRoomUser'])),
        ))
        db.send_create_signal(u'db', ['ServiceActivity'])

        # Adding model 'SilencingReagent'
        db.create_table(u'silencing_reagent', (
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Reagent'], primary_key=True)),
            ('sequence', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('silencing_reagent_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('vendor_gene', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'vendor_reagent', unique=True, null=True, to=orm['db.Gene'])),
            ('facility_gene', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'facility_reagent', unique=True, null=True, to=orm['db.Gene'])),
            ('is_restricted_sequence', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'db', ['SilencingReagent'])

        # Adding model 'SilencingReagentDuplexWells'
        db.create_table(u'silencing_reagent_duplex_wells', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('silencing_reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.SilencingReagent'])),
            ('well', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'])),
        ))
        db.send_create_signal(u'db', ['SilencingReagentDuplexWells'])

        # Adding model 'SmallMoleculeChembankId'
        db.create_table(u'small_molecule_chembank_id', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.SmallMoleculeReagent'])),
            ('chembank_id', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['SmallMoleculeChembankId'])

        # Adding model 'SmallMoleculeChemblId'
        db.create_table(u'small_molecule_chembl_id', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.SmallMoleculeReagent'])),
            ('chembl_id', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['SmallMoleculeChemblId'])

        # Adding model 'SmallMoleculeCherryPickRequest'
        db.create_table(u'small_molecule_cherry_pick_request', (
            ('cherry_pick_request', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.CherryPickRequest'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['SmallMoleculeCherryPickRequest'])

        # Adding model 'SmallMoleculeCompoundName'
        db.create_table(u'small_molecule_compound_name', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.SmallMoleculeReagent'])),
            ('compound_name', self.gf('django.db.models.fields.TextField')()),
            ('ordinal', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['SmallMoleculeCompoundName'])

        # Adding model 'SmallMoleculePubchemCid'
        db.create_table(u'small_molecule_pubchem_cid', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.SmallMoleculeReagent'])),
            ('pubchem_cid', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['SmallMoleculePubchemCid'])

        # Adding model 'SmallMoleculeReagent'
        db.create_table(u'small_molecule_reagent', (
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Reagent'], primary_key=True)),
            ('inchi', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('molecular_formula', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('molecular_mass', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=15, decimal_places=9, blank=True)),
            ('molecular_weight', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=15, decimal_places=9, blank=True)),
            ('smiles', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('salt_form_id', self.gf('django.db.models.fields.IntegerField')(null=True, blank=True)),
            ('is_restricted_structure', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'db', ['SmallMoleculeReagent'])

        # Adding model 'StudyReagentLink'
        db.create_table(u'study_reagent_link', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('study', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Screen'])),
            ('reagent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Reagent'])),
        ))
        db.send_create_signal(u'db', ['StudyReagentLink'])

        # Adding model 'TransfectionAgent'
        db.create_table(u'transfection_agent', (
            ('transfection_agent_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('value', self.gf('django.db.models.fields.TextField')(unique=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'db', ['TransfectionAgent'])

        # Adding model 'Well'
        db.create_table(u'well', (
            ('well_id', self.gf('django.db.models.fields.TextField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('plate_number', self.gf('django.db.models.fields.IntegerField')()),
            ('well_name', self.gf('django.db.models.fields.TextField')()),
            ('facility_id', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('library_well_type', self.gf('django.db.models.fields.TextField')()),
            ('library', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Library'])),
            ('deprecation_admin_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'], null=True, blank=True)),
            ('is_deprecated', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('latest_released_reagent', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name=u'reagent_well', null=True, to=orm['db.Reagent'])),
            ('molar_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=13, decimal_places=12, blank=True)),
            ('mg_ml_concentration', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=5, decimal_places=3, blank=True)),
        ))
        db.send_create_signal(u'db', ['Well'])

        # Adding model 'WellVolumeAdjustment'
        db.create_table(u'well_volume_adjustment', (
            ('well_volume_adjustment_id', self.gf('django.db.models.fields.IntegerField')(primary_key=True)),
            ('version', self.gf('django.db.models.fields.IntegerField')()),
            ('well', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Well'])),
            ('lab_cherry_pick', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.LabCherryPick'], null=True, blank=True)),
            ('well_volume_correction_activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.WellVolumeCorrectionActivity'], null=True, blank=True)),
            ('volume', self.gf('django.db.models.fields.DecimalField')(null=True, max_digits=10, decimal_places=9, blank=True)),
            ('copy', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.Copy'])),
        ))
        db.send_create_signal(u'db', ['WellVolumeAdjustment'])

        # Adding model 'WellVolumeCorrectionActivity'
        db.create_table(u'well_volume_correction_activity', (
            ('activity', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['db.AdministrativeActivity'], primary_key=True)),
        ))
        db.send_create_signal(u'db', ['WellVolumeCorrectionActivity'])


    def backwards(self, orm):
        # Removing index on 'ScreenStatusItem', fields ['screen', 'status', 'status_date']
        db.delete_index(u'screen_status_item', ['screen_id', 'status', 'status_date'])

        # Deleting model 'LegacySmallMoleculeCasNumber'
        db.delete_table(u'_legacy_small_molecule_cas_number')

        # Deleting model 'AbaseTestset'
        db.delete_table(u'abase_testset')

        # Deleting model 'Activity'
        db.delete_table(u'activity')

        # Deleting model 'ActivityUpdateActivity'
        db.delete_table(u'activity_update_activity')

        # Deleting model 'AdministrativeActivity'
        db.delete_table(u'administrative_activity')

        # Deleting model 'AdministratorUser'
        db.delete_table(u'administrator_user')

        # Deleting model 'AnnotationType'
        db.delete_table(u'annotation_type')

        # Deleting model 'AnnotationValue'
        db.delete_table(u'annotation_value')

        # Deleting model 'AssayPlate'
        db.delete_table(u'assay_plate')

        # Deleting model 'AssayWell'
        db.delete_table(u'assay_well')

        # Deleting model 'AttachedFile'
        db.delete_table(u'attached_file')

        # Deleting model 'AttachedFileType'
        db.delete_table(u'attached_file_type')

        # Deleting model 'AttachedFileUpdateActivity'
        db.delete_table(u'attached_file_update_activity')

        # Deleting model 'Cell'
        db.delete_table(u'cell')

        # Deleting model 'CellGrowthProperties'
        db.delete_table(u'cell_growth_properties')

        # Deleting model 'CellLine'
        db.delete_table(u'cell_line')

        # Deleting model 'CellLineage'
        db.delete_table(u'cell_lineage')

        # Deleting model 'CellMarkers'
        db.delete_table(u'cell_markers')

        # Deleting model 'CellRelatedProjects'
        db.delete_table(u'cell_related_projects')

        # Deleting model 'CellUpdateActivity'
        db.delete_table(u'cell_update_activity')

        # Deleting model 'ChecklistItem'
        db.delete_table(u'checklist_item')

        # Deleting model 'ChecklistItemEvent'
        db.delete_table(u'checklist_item_event')

        # Deleting model 'ChecklistItemEventUpdateActivity'
        db.delete_table(u'checklist_item_event_update_activity')

        # Deleting model 'CherryPickAssayPlate'
        db.delete_table(u'cherry_pick_assay_plate')

        # Deleting model 'CherryPickAssayPlateScreeningLink'
        db.delete_table(u'cherry_pick_assay_plate_screening_link')

        # Deleting model 'CherryPickLiquidTransfer'
        db.delete_table(u'cherry_pick_liquid_transfer')

        # Deleting model 'CherryPickRequest'
        db.delete_table(u'cherry_pick_request')

        # Deleting model 'CherryPickRequestEmptyWell'
        db.delete_table(u'cherry_pick_request_empty_well')

        # Deleting model 'CherryPickRequestUpdateActivity'
        db.delete_table(u'cherry_pick_request_update_activity')

        # Deleting model 'CherryPickScreening'
        db.delete_table(u'cherry_pick_screening')

        # Deleting model 'CollaboratorLink'
        db.delete_table(u'collaborator_link')

        # Deleting model 'Copy'
        db.delete_table(u'copy')

        # Deleting model 'CopyUpdateActivity'
        db.delete_table(u'copy_update_activity')

        # Deleting model 'DataColumn'
        db.delete_table(u'data_column')

        # Deleting model 'DataColumnDerivedFromLink'
        db.delete_table(u'data_column_derived_from_link')

        # Deleting model 'EquipmentUsed'
        db.delete_table(u'equipment_used')

        # Deleting model 'ExperimentalCellInformation'
        db.delete_table(u'experimental_cell_information')

        # Deleting model 'FundingSupport'
        db.delete_table(u'funding_support')

        # Deleting model 'Gene'
        db.delete_table(u'gene')

        # Deleting model 'GeneGenbankAccessionNumber'
        db.delete_table(u'gene_genbank_accession_number')

        # Deleting model 'GeneOldEntrezgeneId'
        db.delete_table(u'gene_old_entrezgene_id')

        # Deleting model 'GeneOldEntrezgeneSymbol'
        db.delete_table(u'gene_old_entrezgene_symbol')

        # Deleting model 'GeneSymbol'
        db.delete_table(u'gene_symbol')

        # Deleting model 'LabActivity'
        db.delete_table(u'lab_activity')

        # Deleting model 'LabAffiliation'
        db.delete_table(u'lab_affiliation')

        # Deleting model 'LabCherryPick'
        db.delete_table(u'lab_cherry_pick')

        # Deleting model 'Library'
        db.delete_table(u'library')

        # Deleting model 'LibraryContentsVersion'
        db.delete_table(u'library_contents_version')

        # Deleting model 'LibraryScreening'
        db.delete_table(u'library_screening')

        # Deleting model 'LibraryUpdateActivity'
        db.delete_table(u'library_update_activity')

        # Deleting model 'Molfile'
        db.delete_table(u'molfile')

        # Deleting model 'NaturalProductReagent'
        db.delete_table(u'natural_product_reagent')

        # Deleting model 'Plate'
        db.delete_table(u'plate')

        # Deleting model 'PlateLocation'
        db.delete_table(u'plate_location')

        # Deleting model 'PlateUpdateActivity'
        db.delete_table(u'plate_update_activity')

        # Deleting model 'PrimaryCell'
        db.delete_table(u'primary_cell')

        # Deleting model 'Publication'
        db.delete_table(u'publication')

        # Deleting model 'Reagent'
        db.delete_table(u'reagent')

        # Deleting model 'ReagentPublicationLink'
        db.delete_table(u'reagent_publication_link')

        # Deleting model 'ResultValue'
        db.delete_table(u'result_value')

        # Deleting model 'RnaiCherryPickRequest'
        db.delete_table(u'rnai_cherry_pick_request')

        # Deleting model 'SchemaHistory'
        db.delete_table(u'schema_history')

        # Deleting model 'Screen'
        db.delete_table(u'screen')

        # Deleting model 'ScreenBillingItem'
        db.delete_table(u'screen_billing_item')

        # Deleting model 'ScreenFundingSupportLink'
        db.delete_table(u'screen_funding_support_link')

        # Deleting model 'ScreenKeyword'
        db.delete_table(u'screen_keyword')

        # Deleting model 'ScreenPublicationLink'
        db.delete_table(u'screen_publication_link')

        # Deleting model 'ScreenResult'
        db.delete_table(u'screen_result')

        # Deleting model 'ScreenResultUpdateActivity'
        db.delete_table(u'screen_result_update_activity')

        # Deleting model 'ScreenStatusItem'
        db.delete_table(u'screen_status_item')

        # Deleting model 'ScreenUpdateActivity'
        db.delete_table(u'screen_update_activity')

        # Deleting model 'ScreenerCherryPick'
        db.delete_table(u'screener_cherry_pick')

        # Deleting model 'Screening'
        db.delete_table(u'screening')

        # Deleting model 'ScreeningRoomUser'
        db.delete_table(u'screening_room_user')

        # Deleting model 'LabHead'
        db.delete_table(u'lab_head')

        # Deleting model 'ScreeningRoomUserFacilityUsageRole'
        db.delete_table(u'screening_room_user_facility_usage_role')

        # Deleting model 'ScreensaverUser'
        db.delete_table(u'screensaver_user')

        # Deleting model 'ScreensaverUserRole'
        db.delete_table(u'screensaver_user_role')

        # Deleting model 'ScreensaverUserUpdateActivity'
        db.delete_table(u'screensaver_user_update_activity')

        # Deleting model 'ServiceActivity'
        db.delete_table(u'service_activity')

        # Deleting model 'SilencingReagent'
        db.delete_table(u'silencing_reagent')

        # Deleting model 'SilencingReagentDuplexWells'
        db.delete_table(u'silencing_reagent_duplex_wells')

        # Deleting model 'SmallMoleculeChembankId'
        db.delete_table(u'small_molecule_chembank_id')

        # Deleting model 'SmallMoleculeChemblId'
        db.delete_table(u'small_molecule_chembl_id')

        # Deleting model 'SmallMoleculeCherryPickRequest'
        db.delete_table(u'small_molecule_cherry_pick_request')

        # Deleting model 'SmallMoleculeCompoundName'
        db.delete_table(u'small_molecule_compound_name')

        # Deleting model 'SmallMoleculePubchemCid'
        db.delete_table(u'small_molecule_pubchem_cid')

        # Deleting model 'SmallMoleculeReagent'
        db.delete_table(u'small_molecule_reagent')

        # Deleting model 'StudyReagentLink'
        db.delete_table(u'study_reagent_link')

        # Deleting model 'TransfectionAgent'
        db.delete_table(u'transfection_agent')

        # Deleting model 'Well'
        db.delete_table(u'well')

        # Deleting model 'WellVolumeAdjustment'
        db.delete_table(u'well_volume_adjustment')

        # Deleting model 'WellVolumeCorrectionActivity'
        db.delete_table(u'well_volume_correction_activity')


    models = {
        u'db.abasetestset': {
            'Meta': {'object_name': 'AbaseTestset', 'db_table': "u'abase_testset'"},
            'abase_testset_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'testset_date': ('django.db.models.fields.DateField', [], {}),
            'testset_name': ('django.db.models.fields.TextField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.activity': {
            'Meta': {'object_name': 'Activity', 'db_table': "u'activity'"},
            'activity_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'activities_created'", 'null': 'True', 'to': u"orm['db.ScreensaverUser']"}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_of_activity': ('django.db.models.fields.DateField', [], {}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'performed_by': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'activities_performed'", 'to': u"orm['db.ScreensaverUser']"}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.activityupdateactivity': {
            'Meta': {'object_name': 'ActivityUpdateActivity', 'db_table': "u'activity_update_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.administrativeactivity': {
            'Meta': {'object_name': 'AdministrativeActivity', 'db_table': "u'administrative_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']", 'primary_key': 'True'}),
            'administrative_activity_type': ('django.db.models.fields.TextField', [], {})
        },
        u'db.administratoruser': {
            'Meta': {'object_name': 'AdministratorUser', 'db_table': "u'administrator_user'"},
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'primary_key': 'True'})
        },
        u'db.annotationtype': {
            'Meta': {'object_name': 'AnnotationType', 'db_table': "u'annotation_type'"},
            'annotation_type_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_numeric': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'study': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.annotationvalue': {
            'Meta': {'object_name': 'AnnotationValue', 'db_table': "u'annotation_value'"},
            'annotation_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AnnotationType']", 'null': 'True', 'blank': 'True'}),
            'annotation_value_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'numeric_value': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']", 'null': 'True', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.assayplate': {
            'Meta': {'object_name': 'AssayPlate', 'db_table': "u'assay_plate'"},
            'assay_plate_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'library_screening': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LibraryScreening']", 'null': 'True', 'blank': 'True'}),
            'plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Plate']", 'null': 'True', 'blank': 'True'}),
            'plate_number': ('django.db.models.fields.IntegerField', [], {}),
            'replicate_ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'screen_result_data_loading': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.assaywell': {
            'Meta': {'object_name': 'AssayWell', 'db_table': "u'assay_well'"},
            'assay_well_control_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_well_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'confirmed_positive_value': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_positive': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'screen_result': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenResult']"}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"})
        },
        u'db.attachedfile': {
            'Meta': {'object_name': 'AttachedFile', 'db_table': "u'attached_file'"},
            'attached_file_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'attached_file_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AttachedFileType']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'file_contents': ('django.db.models.fields.TextField', [], {}),
            'file_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'filename': ('django.db.models.fields.TextField', [], {}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']", 'null': 'True', 'blank': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']", 'null': 'True', 'blank': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']", 'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.attachedfiletype': {
            'Meta': {'object_name': 'AttachedFileType', 'db_table': "u'attached_file_type'"},
            'attached_file_type_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'for_entity_type': ('django.db.models.fields.CharField', [], {'max_length': '31'}),
            'value': ('django.db.models.fields.TextField', [], {})
        },
        u'db.attachedfileupdateactivity': {
            'Meta': {'object_name': 'AttachedFileUpdateActivity', 'db_table': "u'attached_file_update_activity'"},
            'attached_file': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AttachedFile']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.cell': {
            'Meta': {'object_name': 'Cell', 'db_table': "u'cell'"},
            'alternate_id': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'alternate_name': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'batch_id': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'cell_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'cell_type': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'cell_type_detail': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'center_name': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'center_specific_id': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'clo_id': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'disease': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'disease_detail': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'facility_id': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '255'}),
            'genetic_modification': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'mutations_explicit': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'mutations_reference': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'organ': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'organism': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'organism_gender': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'recommended_culture_conditions': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'tissue': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'vendor': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'vendor_catalog_id': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'verification': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'verification_reference_profile': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.cellgrowthproperties': {
            'Meta': {'object_name': 'CellGrowthProperties', 'db_table': "u'cell_growth_properties'"},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Cell']"}),
            'growth_property': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.cellline': {
            'Meta': {'object_name': 'CellLine', 'db_table': "u'cell_line'"},
            'cell_line_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.celllineage': {
            'Meta': {'object_name': 'CellLineage', 'db_table': "u'cell_lineage'"},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Cell']", 'primary_key': 'True'})
        },
        u'db.cellmarkers': {
            'Meta': {'object_name': 'CellMarkers', 'db_table': "u'cell_markers'"},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.PrimaryCell']"}),
            'cell_markers': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.cellrelatedprojects': {
            'Meta': {'object_name': 'CellRelatedProjects', 'db_table': "u'cell_related_projects'"},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Cell']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'related_project': ('django.db.models.fields.TextField', [], {})
        },
        u'db.cellupdateactivity': {
            'Meta': {'object_name': 'CellUpdateActivity', 'db_table': "u'cell_update_activity'"},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Cell']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'unique': 'True'})
        },
        u'db.checklistitem': {
            'Meta': {'object_name': 'ChecklistItem', 'db_table': "u'checklist_item'"},
            'checklist_item_group': ('django.db.models.fields.TextField', [], {}),
            'checklist_item_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'is_expirable': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'item_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'order_statistic': ('django.db.models.fields.IntegerField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.checklistitemevent': {
            'Meta': {'object_name': 'ChecklistItemEvent', 'db_table': "u'checklist_item_event'"},
            'checklist_item_event_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'checklist_item_id': ('django.db.models.fields.IntegerField', [], {}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_performed': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'is_expiration': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_not_applicable': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'screening_room_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"})
        },
        u'db.checklistitemeventupdateactivity': {
            'Meta': {'object_name': 'ChecklistItemEventUpdateActivity', 'db_table': "u'checklist_item_event_update_activity'"},
            'checklist_item_event': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ChecklistItemEvent']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.cherrypickassayplate': {
            'Meta': {'object_name': 'CherryPickAssayPlate', 'db_table': "u'cherry_pick_assay_plate'"},
            'assay_plate_type': ('django.db.models.fields.TextField', [], {}),
            'attempt_ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'cherry_pick_assay_plate_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'cherry_pick_assay_plate_type': ('django.db.models.fields.CharField', [], {'max_length': '31'}),
            'cherry_pick_liquid_transfer': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickLiquidTransfer']", 'null': 'True', 'blank': 'True'}),
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            'legacy_plate_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'plate_ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.cherrypickassayplatescreeninglink': {
            'Meta': {'object_name': 'CherryPickAssayPlateScreeningLink', 'db_table': "u'cherry_pick_assay_plate_screening_link'"},
            'cherry_pick_assay_plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickAssayPlate']"}),
            'cherry_pick_screening': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickScreening']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.cherrypickliquidtransfer': {
            'Meta': {'object_name': 'CherryPickLiquidTransfer', 'db_table': "u'cherry_pick_liquid_transfer'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabActivity']", 'primary_key': 'True'}),
            'status': ('django.db.models.fields.TextField', [], {})
        },
        u'db.cherrypickrequest': {
            'Meta': {'object_name': 'CherryPickRequest', 'db_table': "u'cherry_pick_request'"},
            'assay_plate_type': ('django.db.models.fields.TextField', [], {}),
            'assay_protocol_comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cherry_pick_assay_protocols_followed': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cherry_pick_followup_results_status': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'cherry_pick_request_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_requested': ('django.db.models.fields.DateField', [], {}),
            'date_volume_approved': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'is_randomized_assay_plate_layout': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'keep_source_plate_cherry_picks_together': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'legacy_cherry_pick_request_number': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'max_skipped_wells_per_plate': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'number_unfulfilled_lab_cherry_picks': ('django.db.models.fields.IntegerField', [], {}),
            'requested_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'transfer_volume_per_well_approved': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'transfer_volume_per_well_requested': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'volume_approved_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministratorUser']", 'null': 'True', 'blank': 'True'})
        },
        u'db.cherrypickrequestemptywell': {
            'Meta': {'object_name': 'CherryPickRequestEmptyWell', 'db_table': "u'cherry_pick_request_empty_well'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'well_name': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'})
        },
        u'db.cherrypickrequestupdateactivity': {
            'Meta': {'object_name': 'CherryPickRequestUpdateActivity', 'db_table': "u'cherry_pick_request_update_activity'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'unique': 'True'})
        },
        u'db.cherrypickscreening': {
            'Meta': {'object_name': 'CherryPickScreening', 'db_table': "u'cherry_pick_screening'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screening']", 'primary_key': 'True'}),
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"})
        },
        u'db.collaboratorlink': {
            'Meta': {'object_name': 'CollaboratorLink', 'db_table': "u'collaborator_link'"},
            'collaborator': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.copy': {
            'Meta': {'object_name': 'Copy', 'db_table': "u'copy'"},
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'copy_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_plated': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'max_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'max_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'min_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'min_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'plate_locations_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'plates_available': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'primary_plate_location_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'primary_plate_status': ('django.db.models.fields.TextField', [], {}),
            'primary_well_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'primary_well_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'usage_type': ('django.db.models.fields.TextField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well_concentration_dilution_factor': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '8', 'decimal_places': '2', 'blank': 'True'})
        },
        u'db.copyupdateactivity': {
            'Meta': {'object_name': 'CopyUpdateActivity', 'db_table': "u'copy_update_activity'"},
            'copy_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'update_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'})
        },
        u'db.datacolumn': {
            'Meta': {'object_name': 'DataColumn', 'db_table': "u'data_column'"},
            'assay_phenotype': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_readout_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'channel': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'data_column_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'data_type': ('django.db.models.fields.TextField', [], {}),
            'decimal_places': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'how_derived': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_derived': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_follow_up_data': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'medium_positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'name': ('django.db.models.fields.TextField', [], {}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'replicate_ordinal': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'screen_result': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenResult']"}),
            'strong_positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'time_point': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'time_point_ordinal': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'weak_positives_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'zdepth_ordinal': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'})
        },
        u'db.datacolumnderivedfromlink': {
            'Meta': {'object_name': 'DataColumnDerivedFromLink', 'db_table': "u'data_column_derived_from_link'"},
            'derived_data_column': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.DataColumn']"}),
            'derived_from_data_column': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'derived_from'", 'to': u"orm['db.DataColumn']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.equipmentused': {
            'Meta': {'object_name': 'EquipmentUsed', 'db_table': "u'equipment_used'"},
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'equipment': ('django.db.models.fields.TextField', [], {}),
            'equipment_used_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'lab_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabActivity']"}),
            'protocol': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.experimentalcellinformation': {
            'Meta': {'object_name': 'ExperimentalCellInformation', 'db_table': "u'experimental_cell_information'"},
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Cell']"}),
            'experimental_cell_information_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.fundingsupport': {
            'Meta': {'object_name': 'FundingSupport', 'db_table': "u'funding_support'"},
            'funding_support_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'unique': 'True', 'blank': 'True'})
        },
        u'db.gene': {
            'Meta': {'object_name': 'Gene', 'db_table': "u'gene'"},
            'entrezgene_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'gene_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'gene_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'species_name': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.genegenbankaccessionnumber': {
            'Meta': {'object_name': 'GeneGenbankAccessionNumber', 'db_table': "u'gene_genbank_accession_number'"},
            'genbank_accession_number': ('django.db.models.fields.TextField', [], {}),
            'gene_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'})
        },
        u'db.geneoldentrezgeneid': {
            'Meta': {'object_name': 'GeneOldEntrezgeneId', 'db_table': "u'gene_old_entrezgene_id'"},
            'gene_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'old_entrezgene_id': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.geneoldentrezgenesymbol': {
            'Meta': {'object_name': 'GeneOldEntrezgeneSymbol', 'db_table': "u'gene_old_entrezgene_symbol'"},
            'gene_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'old_entrezgene_symbol': ('django.db.models.fields.TextField', [], {})
        },
        u'db.genesymbol': {
            'Meta': {'object_name': 'GeneSymbol', 'db_table': "u'gene_symbol'"},
            'entrezgene_symbol': ('django.db.models.fields.TextField', [], {}),
            'gene': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Gene']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.labactivity': {
            'Meta': {'object_name': 'LabActivity', 'db_table': "u'lab_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']", 'primary_key': 'True'}),
            'molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'volume_transferred_per_well_from_library_plates': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'})
        },
        u'db.labaffiliation': {
            'Meta': {'object_name': 'LabAffiliation', 'db_table': "u'lab_affiliation'"},
            'affiliation_category': ('django.db.models.fields.TextField', [], {}),
            'affiliation_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'lab_affiliation_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.labcherrypick': {
            'Meta': {'object_name': 'LabCherryPick', 'db_table': "u'lab_cherry_pick'"},
            'assay_plate_column': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'assay_plate_row': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'cherry_pick_assay_plate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickAssayPlate']", 'null': 'True', 'blank': 'True'}),
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            'lab_cherry_pick_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'screener_cherry_pick': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenerCherryPick']"}),
            'source_well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.labhead': {
            'Meta': {'object_name': 'LabHead', 'db_table': "u'lab_head'"},
            'lab_affiliation': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabAffiliation']", 'null': 'True', 'blank': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']", 'primary_key': 'True'})
        },
        u'db.legacysmallmoleculecasnumber': {
            'Meta': {'object_name': 'LegacySmallMoleculeCasNumber', 'db_table': "u'_legacy_small_molecule_cas_number'"},
            'cas_number': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'smiles': ('django.db.models.fields.CharField', [], {'max_length': '2047'})
        },
        u'db.library': {
            'Meta': {'object_name': 'Library', 'db_table': "u'library'"},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_received': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_screenable': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'description': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'end_plate': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'experimental_well_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'is_pool': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'latest_released_contents_version_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'library_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'library_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'library_type': ('django.db.models.fields.TextField', [], {}),
            'owner_screener': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']", 'null': 'True', 'blank': 'True'}),
            'plate_size': ('django.db.models.fields.TextField', [], {}),
            'provider': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'screen_type': ('django.db.models.fields.TextField', [], {}),
            'screening_status': ('django.db.models.fields.TextField', [], {}),
            'short_name': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'solvent': ('django.db.models.fields.TextField', [], {}),
            'start_plate': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.librarycontentsversion': {
            'Meta': {'object_name': 'LibraryContentsVersion', 'db_table': "u'library_contents_version'"},
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'library_contents_loading_activity': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "u'lcv_load'", 'to': u"orm['db.AdministrativeActivity']"}),
            'library_contents_release_activity': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'lcv_release'", 'null': 'True', 'to': u"orm['db.AdministrativeActivity']"}),
            'library_contents_version_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'version_number': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.libraryscreening': {
            'Meta': {'object_name': 'LibraryScreening', 'db_table': "u'library_screening'"},
            'abase_testset_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screening']", 'primary_key': 'True'}),
            'is_for_external_library_plates': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'libraries_screened_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'library_plates_screened_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'screened_experimental_well_count': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.libraryupdateactivity': {
            'Meta': {'object_name': 'LibraryUpdateActivity', 'db_table': "u'library_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.molfile': {
            'Meta': {'object_name': 'Molfile', 'db_table': "u'molfile'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'molfile': ('django.db.models.fields.TextField', [], {}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.SmallMoleculeReagent']", 'unique': 'True'})
        },
        u'db.naturalproductreagent': {
            'Meta': {'object_name': 'NaturalProductReagent', 'db_table': "u'natural_product_reagent'"},
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']", 'primary_key': 'True'})
        },
        u'db.plate': {
            'Meta': {'object_name': 'Plate', 'db_table': "u'plate'"},
            'copy': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Copy']"}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'facility_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'max_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'max_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'min_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'min_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'plate_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'plate_location': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.PlateLocation']", 'null': 'True', 'blank': 'True'}),
            'plate_number': ('django.db.models.fields.IntegerField', [], {}),
            'plate_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'plated_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True', 'null': 'True', 'blank': 'True'}),
            'primary_well_mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'primary_well_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'quadrant': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'retired_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True', 'null': 'True', 'blank': 'True'}),
            'status': ('django.db.models.fields.TextField', [], {}),
            'stock_plate_number': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well_volume': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'})
        },
        u'db.platelocation': {
            'Meta': {'object_name': 'PlateLocation', 'db_table': "u'plate_location'"},
            'bin': ('django.db.models.fields.TextField', [], {}),
            'freezer': ('django.db.models.fields.TextField', [], {}),
            'plate_location_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'room': ('django.db.models.fields.TextField', [], {}),
            'shelf': ('django.db.models.fields.TextField', [], {})
        },
        u'db.plateupdateactivity': {
            'Meta': {'object_name': 'PlateUpdateActivity', 'db_table': "u'plate_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'plate_id': ('django.db.models.fields.IntegerField', [], {}),
            'update_activity_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'})
        },
        u'db.primarycell': {
            'Meta': {'object_name': 'PrimaryCell', 'db_table': "u'primary_cell'"},
            'age_in_years': ('django.db.models.fields.IntegerField', [], {}),
            'cell': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Cell']", 'primary_key': 'True'}),
            'donor_ethnicity': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'donor_health_status': ('django.db.models.fields.CharField', [], {'max_length': '255', 'blank': 'True'}),
            'passage_number': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.publication': {
            'Meta': {'object_name': 'Publication', 'db_table': "u'publication'"},
            'attached_file': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AttachedFile']", 'unique': 'True', 'null': 'True', 'blank': 'True'}),
            'authors': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'journal': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'pages': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'publication_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'pubmed_central_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pubmed_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'title': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'volume': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'year_published': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.reagent': {
            'Meta': {'object_name': 'Reagent', 'db_table': "u'reagent'"},
            'facility_batch_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'library_contents_version': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LibraryContentsVersion']"}),
            'reagent_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'vendor_batch_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'vendor_identifier': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'vendor_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"})
        },
        u'db.reagentpublicationlink': {
            'Meta': {'object_name': 'ReagentPublicationLink', 'db_table': "u'reagent_publication_link'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'publication_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"})
        },
        u'db.resultvalue': {
            'Meta': {'object_name': 'ResultValue', 'db_table': "u'result_value'"},
            'assay_well_control_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'data_column': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.DataColumn']", 'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_exclude': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'is_positive': ('django.db.models.fields.NullBooleanField', [], {'null': 'True', 'blank': 'True'}),
            'numeric_value': ('django.db.models.fields.FloatField', [], {'null': 'True', 'blank': 'True'}),
            'result_value_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']", 'null': 'True', 'blank': 'True'})
        },
        u'db.rnaicherrypickrequest': {
            'Meta': {'object_name': 'RnaiCherryPickRequest', 'db_table': "u'rnai_cherry_pick_request'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']", 'primary_key': 'True'})
        },
        u'db.schemahistory': {
            'Meta': {'object_name': 'SchemaHistory', 'db_table': "u'schema_history'"},
            'comment': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'date_updated': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'screensaver_revision': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'})
        },
        u'db.screen': {
            'Meta': {'object_name': 'Screen', 'db_table': "u'screen'"},
            'abase_protocol_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'abase_study_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'amount_to_be_charged_for_screen': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '9', 'decimal_places': '2', 'blank': 'True'}),
            'assay_plates_screened_count': ('django.db.models.fields.IntegerField', [], {}),
            'billing_comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'billing_info_return_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'cell_line': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CellLine']", 'null': 'True', 'blank': 'True'}),
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'coms_approval_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'coms_registration_number': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'data_meeting_complete': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_meeting_scheduled': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_privacy_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_privacy_expiration_notified_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'data_sharing_level': ('django.db.models.fields.IntegerField', [], {}),
            'date_charged': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_completed5kcompounds': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_faxed_to_billing_department': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_of_application': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'facilities_and_administration_charge': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '9', 'decimal_places': '2', 'blank': 'True'}),
            'facility_id': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'fee_form_requested_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'fee_form_requested_initials': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'image_url': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_billing_for_supplies_only': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_fee_form_on_file': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'lab_head': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabHead']", 'null': 'True', 'blank': 'True'}),
            'lead_screener': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']", 'null': 'True', 'blank': 'True'}),
            'libraries_screened_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'library_plates_data_analyzed_count': ('django.db.models.fields.IntegerField', [], {}),
            'library_plates_data_loaded_count': ('django.db.models.fields.IntegerField', [], {}),
            'library_plates_screened_count': ('django.db.models.fields.IntegerField', [], {}),
            'max_allowed_data_privacy_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'max_data_loaded_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'max_screened_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'min_allowed_data_privacy_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'min_data_loaded_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'min_screened_replicate_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'perturbagen_molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'perturbagen_ug_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'pin_transfer_admin_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'null': 'True', 'blank': 'True'}),
            'project_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'project_phase': ('django.db.models.fields.TextField', [], {}),
            'pubchem_assay_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'pubchem_deposited_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'publishable_protocol': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'publishable_protocol_comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'publishable_protocol_date_entered': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'publishable_protocol_entered_by': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'screen_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'screen_type': ('django.db.models.fields.TextField', [], {}),
            'screened_experimental_well_count': ('django.db.models.fields.IntegerField', [], {}),
            'see_comments': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'species': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'study_type': ('django.db.models.fields.TextField', [], {}),
            'summary': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'title': ('django.db.models.fields.TextField', [], {}),
            'to_be_requested': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'total_plated_lab_cherry_picks': ('django.db.models.fields.IntegerField', [], {}),
            'transfection_agent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.TransfectionAgent']", 'null': 'True', 'blank': 'True'}),
            'unique_screened_experimental_well_count': ('django.db.models.fields.IntegerField', [], {}),
            'url': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well_studied': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']", 'null': 'True', 'blank': 'True'})
        },
        u'db.screenbillingitem': {
            'Meta': {'object_name': 'ScreenBillingItem', 'db_table': "u'screen_billing_item'"},
            'amount': ('django.db.models.fields.DecimalField', [], {'max_digits': '9', 'decimal_places': '2'}),
            'date_sent_for_billing': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'item_to_be_charged': ('django.db.models.fields.TextField', [], {}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screenercherrypick': {
            'Meta': {'object_name': 'ScreenerCherryPick', 'db_table': "u'screener_cherry_pick'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']"}),
            'screened_well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"}),
            'screener_cherry_pick_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.screenfundingsupportlink': {
            'Meta': {'object_name': 'ScreenFundingSupportLink', 'db_table': "u'screen_funding_support_link'"},
            'funding_support': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.FundingSupport']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screening': {
            'Meta': {'object_name': 'Screening', 'db_table': "u'screening'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabActivity']", 'primary_key': 'True'}),
            'assay_protocol': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_protocol_last_modified_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'assay_protocol_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'assay_well_volume': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'number_of_replicates': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'volume_transferred_per_well_to_assay_plates': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'})
        },
        u'db.screeningroomuser': {
            'Meta': {'object_name': 'ScreeningRoomUser', 'db_table': "u'screening_room_user'"},
            'coms_crhba_permit_number': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'coms_crhba_permit_principal_investigator': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'lab_head': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabHead']", 'null': 'True', 'blank': 'True'}),
            'last_notified_rnaiua_checklist_item_event': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'rnai_ua_user'", 'null': 'True', 'to': u"orm['db.ChecklistItemEvent']"}),
            'last_notified_smua_checklist_item_event': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'smua_user'", 'null': 'True', 'to': u"orm['db.ChecklistItemEvent']"}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'primary_key': 'True'}),
            'user_classification': ('django.db.models.fields.TextField', [], {})
        },
        u'db.screeningroomuserfacilityusagerole': {
            'Meta': {'object_name': 'ScreeningRoomUserFacilityUsageRole', 'db_table': "u'screening_room_user_facility_usage_role'"},
            'facility_usage_role': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screening_room_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"})
        },
        u'db.screenkeyword': {
            'Meta': {'object_name': 'ScreenKeyword', 'db_table': "u'screen_keyword'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'keyword': ('django.db.models.fields.TextField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screenpublicationlink': {
            'Meta': {'object_name': 'ScreenPublicationLink', 'db_table': "u'screen_publication_link'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'publication_id': ('django.db.models.fields.IntegerField', [], {'unique': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.screenresult': {
            'Meta': {'object_name': 'ScreenResult', 'db_table': "u'screen_result'"},
            'channel_count': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'experimental_well_count': ('django.db.models.fields.IntegerField', [], {}),
            'replicate_count': ('django.db.models.fields.IntegerField', [], {}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'screen_result_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.screenresultupdateactivity': {
            'Meta': {'object_name': 'ScreenResultUpdateActivity', 'db_table': "u'screen_result_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen_result': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreenResult']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.screensaveruser': {
            'Meta': {'object_name': 'ScreensaverUser', 'db_table': "u'screensaver_user'"},
            'comments': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']", 'null': 'True', 'blank': 'True'}),
            'date_created': ('django.db.models.fields.DateTimeField', [], {}),
            'date_loaded': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'date_publicly_available': ('django.db.models.fields.DateTimeField', [], {'null': 'True', 'blank': 'True'}),
            'digested_password': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'ecommons_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'email': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'first_name': ('django.db.models.fields.TextField', [], {}),
            'harvard_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'harvard_id_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'harvard_id_requested_expiration_date': ('django.db.models.fields.DateField', [], {'null': 'True', 'blank': 'True'}),
            'last_name': ('django.db.models.fields.TextField', [], {}),
            'login_id': ('django.db.models.fields.TextField', [], {'unique': 'True', 'blank': 'True'}),
            'mailing_address': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'phone': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'screensaver_user_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {'blank': 'True'})
        },
        u'db.screensaveruserrole': {
            'Meta': {'object_name': 'ScreensaverUserRole', 'db_table': "u'screensaver_user_role'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']"}),
            'screensaver_user_role': ('django.db.models.fields.TextField', [], {})
        },
        u'db.screensaveruserupdateactivity': {
            'Meta': {'object_name': 'ScreensaverUserUpdateActivity', 'db_table': "u'screensaver_user_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screensaver_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreensaverUser']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.screenstatusitem': {
            'Meta': {'object_name': 'ScreenStatusItem', 'db_table': "u'screen_status_item'", 'index_together': "((u'screen', u'status', u'status_date'),)"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'status': ('django.db.models.fields.TextField', [], {}),
            'status_date': ('django.db.models.fields.DateField', [], {})
        },
        u'db.screenupdateactivity': {
            'Meta': {'object_name': 'ScreenUpdateActivity', 'db_table': "u'screen_update_activity'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"}),
            'update_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']"})
        },
        u'db.serviceactivity': {
            'Meta': {'object_name': 'ServiceActivity', 'db_table': "u'service_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Activity']", 'primary_key': 'True'}),
            'service_activity_type': ('django.db.models.fields.TextField', [], {}),
            'serviced_screen': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']", 'null': 'True', 'blank': 'True'}),
            'serviced_user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.ScreeningRoomUser']"})
        },
        u'db.silencingreagent': {
            'Meta': {'object_name': 'SilencingReagent', 'db_table': "u'silencing_reagent'"},
            'facility_gene': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'facility_reagent'", 'unique': 'True', 'null': 'True', 'to': u"orm['db.Gene']"}),
            'is_restricted_sequence': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']", 'primary_key': 'True'}),
            'sequence': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'silencing_reagent_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'vendor_gene': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'vendor_reagent'", 'unique': 'True', 'null': 'True', 'to': u"orm['db.Gene']"})
        },
        u'db.silencingreagentduplexwells': {
            'Meta': {'object_name': 'SilencingReagentDuplexWells', 'db_table': "u'silencing_reagent_duplex_wells'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'silencing_reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.SilencingReagent']"}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"})
        },
        u'db.smallmoleculechembankid': {
            'Meta': {'object_name': 'SmallMoleculeChembankId', 'db_table': "u'small_molecule_chembank_id'"},
            'chembank_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.SmallMoleculeReagent']"})
        },
        u'db.smallmoleculechemblid': {
            'Meta': {'object_name': 'SmallMoleculeChemblId', 'db_table': "u'small_molecule_chembl_id'"},
            'chembl_id': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.SmallMoleculeReagent']"})
        },
        u'db.smallmoleculecherrypickrequest': {
            'Meta': {'object_name': 'SmallMoleculeCherryPickRequest', 'db_table': "u'small_molecule_cherry_pick_request'"},
            'cherry_pick_request': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.CherryPickRequest']", 'primary_key': 'True'})
        },
        u'db.smallmoleculecompoundname': {
            'Meta': {'object_name': 'SmallMoleculeCompoundName', 'db_table': "u'small_molecule_compound_name'"},
            'compound_name': ('django.db.models.fields.TextField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'ordinal': ('django.db.models.fields.IntegerField', [], {}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.SmallMoleculeReagent']"})
        },
        u'db.smallmoleculepubchemcid': {
            'Meta': {'object_name': 'SmallMoleculePubchemCid', 'db_table': "u'small_molecule_pubchem_cid'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'pubchem_cid': ('django.db.models.fields.IntegerField', [], {}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.SmallMoleculeReagent']"})
        },
        u'db.smallmoleculereagent': {
            'Meta': {'object_name': 'SmallMoleculeReagent', 'db_table': "u'small_molecule_reagent'"},
            'inchi': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_restricted_structure': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'molecular_formula': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'molecular_mass': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '15', 'decimal_places': '9', 'blank': 'True'}),
            'molecular_weight': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '15', 'decimal_places': '9', 'blank': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']", 'primary_key': 'True'}),
            'salt_form_id': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'smiles': ('django.db.models.fields.TextField', [], {'blank': 'True'})
        },
        u'db.studyreagentlink': {
            'Meta': {'object_name': 'StudyReagentLink', 'db_table': "u'study_reagent_link'"},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'reagent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Reagent']"}),
            'study': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Screen']"})
        },
        u'db.transfectionagent': {
            'Meta': {'object_name': 'TransfectionAgent', 'db_table': "u'transfection_agent'"},
            'transfection_agent_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'value': ('django.db.models.fields.TextField', [], {'unique': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {})
        },
        u'db.well': {
            'Meta': {'object_name': 'Well', 'db_table': "u'well'"},
            'deprecation_admin_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'null': 'True', 'blank': 'True'}),
            'facility_id': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'is_deprecated': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'latest_released_reagent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "u'reagent_well'", 'null': 'True', 'to': u"orm['db.Reagent']"}),
            'library': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Library']"}),
            'library_well_type': ('django.db.models.fields.TextField', [], {}),
            'mg_ml_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '5', 'decimal_places': '3', 'blank': 'True'}),
            'molar_concentration': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '13', 'decimal_places': '12', 'blank': 'True'}),
            'plate_number': ('django.db.models.fields.IntegerField', [], {}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'well_id': ('django.db.models.fields.TextField', [], {'primary_key': 'True'}),
            'well_name': ('django.db.models.fields.TextField', [], {})
        },
        u'db.wellvolumeadjustment': {
            'Meta': {'object_name': 'WellVolumeAdjustment', 'db_table': "u'well_volume_adjustment'"},
            'copy': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Copy']"}),
            'lab_cherry_pick': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.LabCherryPick']", 'null': 'True', 'blank': 'True'}),
            'version': ('django.db.models.fields.IntegerField', [], {}),
            'volume': ('django.db.models.fields.DecimalField', [], {'null': 'True', 'max_digits': '10', 'decimal_places': '9', 'blank': 'True'}),
            'well': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.Well']"}),
            'well_volume_adjustment_id': ('django.db.models.fields.IntegerField', [], {'primary_key': 'True'}),
            'well_volume_correction_activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.WellVolumeCorrectionActivity']", 'null': 'True', 'blank': 'True'})
        },
        u'db.wellvolumecorrectionactivity': {
            'Meta': {'object_name': 'WellVolumeCorrectionActivity', 'db_table': "u'well_volume_correction_activity'"},
            'activity': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['db.AdministrativeActivity']", 'primary_key': 'True'})
        }
    }

    complete_apps = ['db']