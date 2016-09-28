BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20140013,
current_timestamp,
'copy the new substance ids into the database';

/**
  Purpose: copy the new substance ids into the database efficiently:
  - first write out the substance id table to csv.
  - copy the csv file to a temp table
  - update reagent.substance_id from the table
  
  NOTE: 20151029 - retired; this functionality is now run from the substance
  ID migration script.
**/


CREATE TEMP TABLE new_substance_ids ( 
  reagent_id integer, 
  substance_id text,
  primary key (reagent_id)
);

/** check the current path**/
\! pwd

\copy new_substance_ids from 'new_substance_ids.csv' with CSV

update reagent 
  set substance_id=a.substance_id 
  from new_substance_ids a
  where a.reagent_id=reagent.reagent_id;

COMMIT;
