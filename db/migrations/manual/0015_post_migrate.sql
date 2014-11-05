/**
  Purpose: clean up all of the old reagents that are no longer needed because
  they aren't the "latest_released_reagent":
  
screensaverlims=# select count(*) from reagent;
  count  
---------
 1868502
(1 row)

screensaverlims=# select count(r.*) from reagent r join well w using(well_id) 
  where w.latest_released_reagent_id != r.reagent_id;
 count  
--------
 635410
(1 row)

  ** NOTE: too big for a transaction **
  
  **/

create temp table temp_reagents_to_delete 
  as select reagent_id 
  from reagent r join well w using(well_id) 
  where w.latest_released_reagent_id != r.reagent_id;

delete from small_molecule_reagent r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from molfile r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from small_molecule_compound_name r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from small_molecule_pubchem_cid r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from small_molecule_chembank_id r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from small_molecule_chembl_id r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from silencing_reagent_duplex_wells dw
  using temp_reagents_to_delete d
  where dw.silencingreagent_id=d.reagent_id; 

delete from silencing_reagent r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from natural_product_reagent r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

/** 
 FIXME: these will all be broken:
 study_reagent_link
 reagent_publication_link
 attached_file
 annotation_value
**/
delete from annotation_value r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from study_reagent_link r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 




delete from reagent r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 


