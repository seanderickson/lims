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
  * TODO: Dynamically generate or perform these updates (in code)
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

COMMIT;
