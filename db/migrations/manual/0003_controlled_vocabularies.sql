BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20140003,
current_timestamp,
'migrate different controlled vocabularies';

/**
  Purpose: convert the old controlled vocabulary types to newer names 
  * newer names have no spaces, and only alphanumeric and "_","." chars
  * so that they can be used as URI values







  * TODO: Dynamically generate or perform these updates (see migration 0003)
**/

/**
 library_well_type 
-------------------
 experimental
 <undefined>
 empty
 RNAi buffer
 DMSO
 library control
**/
update well set library_well_type='library_control' where library_well_type='library control';
update well set library_well_type='undefined' where library_well_type='<undefined>';
update well set library_well_type='rnai_buffer' where library_well_type='RNAi buffer';
update well set library_well_type='dmso' where library_well_type='DMSO';

/**
     plate_type      
---------------------
 
 Eppendorf 384 CB PP
 ABgene 384 CB PP
 Eppendorf 96 CB PP
 Nunc 96 VB PS
 Marsh 384 VB PP
 Costar 96 RB PS
 Genetix 384 CB PP
 
**/
update plate set plate_type = 'eppendorf_384' where plate_type='Eppendorf 384 CB PP';
update plate set plate_type = 'eppendorf_96' where plate_type='Eppendorf 96 CB PP';
update plate set plate_type = 'abgene_384' where plate_type='ABgene 384 CB PP';
update plate set plate_type = 'nunc_96' where plate_type='Nunc 96 VB PS';
update plate set plate_type = 'marsh_384' where plate_type='Marsh 384 VB PP';
update plate set plate_type = 'costar_96' where plate_type='Costar 96 RB PS';
update plate set plate_type = 'genetix_384' where plate_type='Genetix 384 CB PP';

/**
             status             
--------------------------------
 Retired
 Lost
 Given Away
 Discarded (volume transferred)
 Not specified
 Discarded
 Not available
 Not created
 Available
**/

update plate set status = 'available' where status='Available';
update plate set status = 'discarded' where status='Discarded';
update plate set status = 'discarded_volume_transferred' where status='Discarded (volume transferred)';
update plate set status = 'given_away' where status='Given Away';
update plate set status = 'lost' where status='Lost';
update plate set status = 'not_available' where status='Not available';
update plate set status = 'not_created' where status='Not created';
update plate set status = 'not_specified' where status='Not specified';
update plate set status = 'retired' where status='Retired';


/**
        usage_type         
---------------------------
 Cherry Pick Source Plates
 Library Screening Plates
 Stock Plates
 96 Stock Plates
**/

update copy set usage_type = 'cherry_pick_source_plates' where usage_type = 'Cherry Pick Source Plates';
update copy set usage_type = 'library_screening_plates' where usage_type = 'Library Screening Plates';
update copy set usage_type = 'stock_plates' where usage_type = 'Stock Plates';
update copy set usage_type = '96_stock_plates' where usage_type = '96 Stock Plates';


/** 
 assay_well_control_type 
-------------------------
 
 assay control
 experimental
 empty
 assay positive control
 buffer
 other
 library control

**/

update assay_well set assay_well_control_type = 'assay_control' where assay_well_control_type = 'assay control';
update assay_well set assay_well_control_type = 'assay_positive_control' where assay_well_control_type = 'assay positive control';
update assay_well set assay_well_control_type = 'library_control' where assay_well_control_type = 'library control';

/**
          data_type           
------------------------------
 Numeric
 Partition Positive Indicator
 Text
 Confirmed Positive Indicator
 Boolean Positive Indicator
**/

update data_column set data_type = 'numeric' where data_type = 'Numeric';
update data_column set data_type = 'partition_positive_indicator' where data_type = 'Partition Positive Indicator';
update data_column set data_type = 'text' where data_type = 'Text';
update data_column set data_type = 'confirmed_positive_indicator' where data_type = 'Confirmed Positive Indicator';
update data_column set data_type = 'boolean_positive_indicator' where data_type = 'Boolean Positive Indicator';

/** 
-- moved to python migration 0003

insert into screen_funding_supports (id,screen_id,funding_support) (
  select nextval('screen_funding_supports_id_seq'), screen.screen_id, value
  from screen
  join screen_funding_support_link fsl on(screen.screen_id=fsl.screen_id)
  join funding_support fs on(fsl.funding_support_id=fs.funding_support_id) 
  order by screen.screen_id,value
  );
**/
COMMIT;
