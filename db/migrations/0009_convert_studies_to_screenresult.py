# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import logging
from django.db import migrations, models

logger = logging.getLogger(__name__)

def convert_studies(apps,schema_editor):
    
    Screen = apps.get_model('db','Screen')
    ScreenResult = apps.get_model('db','ScreenResult')
    AnnotationType = apps.get_model('db', 'AnnotationType')
    AnnotationValue = apps.get_model('db', 'AnnotationValue')
    DataColumn = apps.get_model('db', 'DataColumn')
    
    # 1. Create placeholder ScreenResults for all studies
    for screen in ( Screen.objects.all()
        .filter(project_phase__exact='annotation')
        .order_by('facility_id')):
        
        logger.info('processing: %s', screen.facility_id)
        screen_result = ScreenResult.objects.create(
            screen=screen,
            date_created=screen.date_created,
            date_loaded=screen.date_loaded,
            date_publicly_available=screen.date_publicly_available,
            created_by=screen.created_by)
    
        # 2. Create AssayWells
        schema_editor.execute(
            'insert into assay_well (well_id, plate_number, screen_result_id) '
            'select distinct(reagent.well_id), plate_number, screen_result_id  '
            'from annotation_value av '
            'join annotation_type at using(annotation_type_id) '
            'join screen on(study_id=screen_id) '
            'join screen_result using(screen_id) '
            'join reagent using(reagent_id) '
            'join well on(reagent.well_id=well.well_id) '
            'where screen.facility_id = %s order by well_id;',
            params=[screen.facility_id])
    
        # 3. Migrate AnnotationTypes to DataColumns
        for at in AnnotationType.objects.all().filter(study=screen):
            
            logger.info(
                'processing annotation type: %s, %s', at.annotation_type_id, 
                at.name)
            data_type = 'text'
            if at.is_numeric:
                data_type = 'numeric'
            decimal_places = None
            if at.annotation_type_id == 22588:
                # Study 200003 - weighted average
                decimal_places = 2
            
            data_column = DataColumn.objects.create(
                name=at.name,
                description=at.description,
                ordinal=at.ordinal,
                screen_result=screen_result,
                data_type=data_type,
                decimal_places=decimal_places)

            logger.info('Create result_values for data_column: %d', 
                data_column.data_column_id)
            # 4. Migrate AnnotationValues to ResultValues
            value_column = 'value'
            if data_type == 'numeric':
                value_column = 'numeric_value'
            sql = (
                'insert into result_value (data_column_id, well_id, {value_column}) '
                "select '{data_column_id}', reagent.well_id, {value_column} "
                "from annotation_value av "
                "join reagent using(reagent_id) "
                "where av.annotation_type_id=%s "
                "order by reagent.well_id; ").format(
                    data_column_id=data_column.data_column_id,
                    value_column=value_column),
            logger.info('sql: %s', sql)
            schema_editor.execute(sql[0],params=[at.annotation_type_id])
    
    
    # 5. Test: 
    #     - importer
    # 6. Cleanup: remove AT, AV

    

class Migration(migrations.Migration):

    dependencies = [
        ('db', '0006_screen'),
    ]

    operations = [
        migrations.RunPython(convert_studies),
        
    ]
