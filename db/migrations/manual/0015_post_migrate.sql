/**
  Purpose: clean up all of the old reagents that are no longer needed because
  they aren't the "latest_released_reagent":
  
screensaverlims=# select count(*) from reagent;
  count  
---------
 1868502
(1 row)

screensaverlims=# select count(r.*) from reagent r join well w using(well_id) where w.latest_released_reagent_id != r.reagent_id;
 count  
--------
 635410
(1 row)

  ** NOTE: too big for a transaction **
  
  **/

delete from reagent r 
  using well w 
  where r.well_id=r.well_id 
  and w.latest_released_reagent_id != r.reagent_id;


