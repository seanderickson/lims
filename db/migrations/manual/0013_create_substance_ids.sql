BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20140013,
current_timestamp,
'copy the new substance ids into the database';

/**
  Purpose: copy the new substance ids into the database efficiently
**/


CREATE TEMP TABLE new_substance_ids ( 
  reagent_id integer, 
  substance_id text
  primary key (reagent_id) );

/** check the current path**/
\! pwd

\copy new_substance_ids from 'new_substance_ids.csv' with CSV

update reagent 
  set substance_id=a.substance_id 
  from new_substance_ids a
  where a.reagent_id=reagent.reagent_id;

COMMIT;
