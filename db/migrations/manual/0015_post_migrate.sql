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

create temp table temp_reagents_to_delete (
  reagent_id integer not null unique
  ); 
  
insert into temp_reagents_to_delete
  ( select reagent_id from reagent r join well w using(well_id) 
      where w.latest_released_reagent_id != r.reagent_id );


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

/** REMOVED: 20150114 - table deleted in migration
delete from silencing_reagent_duplex_wells dw
  using temp_reagents_to_delete d
  where dw.silencingreagent_id=d.reagent_id; 
**/
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
 study_reagent_link
 publication
 attached_file
 annotation_value
**/
delete from study_reagent_link r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from publication p 
  using temp_reagents_to_delete d
  where p.reagent_id=d.reagent_id; 

delete from attached_file r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

delete from annotation_value r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 


/** Remove indexes **/
alter table annotation_value  drop constraint fk_annotation_value_to_reagent ;

alter table attached_file  drop constraint fk_attached_file_to_reagent ;

alter table publication  drop constraint publication_reagent_id_fk_reagent_reagent_id ;

alter table small_molecule_reagent drop CONSTRAINT "fkf5a7b431161ea629" ;

alter table natural_product_reagent drop CONSTRAINT fkc0f2d4c161ea629 ;

alter table silencing_reagent drop CONSTRAINT fkba0f3291161ea629 ;

alter table study_reagent_link drop  CONSTRAINT fk_reagent_to_study ;

alter table well drop CONSTRAINT fk37a0ce68c7e7c9 ;

delete from reagent r 
  using temp_reagents_to_delete d
  where r.reagent_id=d.reagent_id; 

alter table annotation_value add constraint fk_annotation_value_to_reagent 
  FOREIGN KEY (reagent_id) REFERENCES reagent(reagent_id);       

alter table attached_file add constraint fk_attached_file_to_reagent
  FOREIGN KEY (reagent_id) REFERENCES reagent(reagent_id);       

alter table publication  add constraint publication_reagent_id_fk_reagent_reagent_id
  FOREIGN KEY (reagent_id) REFERENCES reagent(reagent_id);       
  
alter table small_molecule_reagent add CONSTRAINT fkf5a7b431161ea629
  FOREIGN KEY (reagent_id) REFERENCES reagent(reagent_id);       
  
alter table natural_product_reagent add CONSTRAINT fkc0f2d4c161ea629
  FOREIGN KEY (reagent_id) REFERENCES reagent(reagent_id);       
  
alter table silencing_reagent add CONSTRAINT fkba0f3291161ea629
  FOREIGN KEY (reagent_id) REFERENCES reagent(reagent_id);       

alter table study_reagent_link add CONSTRAINT fk_reagent_to_study
  FOREIGN KEY (reagent_id) REFERENCES reagent(reagent_id);       

/**
alter table well add CONSTRAINT fk37a0ce68c7e7c9
  FOREIGN KEY (latest_released_reagent_id) REFERENCES reagent(reagent_id);       
alter table well drop column latest_released_reagent;
**/

