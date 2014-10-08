BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20140003,
current_timestamp,
'migrate different controlled vocabularies';

/**
  Purpose: convert the old controlled vocabulary types to newer names 
  * newer names have no spaces, and only alphanumeric and "_","." chars
  * for simplicity in URI construction & parsing
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

COMMIT;
