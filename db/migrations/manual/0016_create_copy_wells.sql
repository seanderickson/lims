BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20161012,
current_timestamp,
'create new copywell data using assay plate information';

select current_time;


/** 
  1. Create copy_wells for well_volume_adjustments:
  - set copy_well.volume (remaining_volume) = current plate.well_volume - sum(wva's)
  - set # of adjustments = count(wva's)
  Note: all wva's should be used, not just some of successful, failed, canceled, 
  because wva is created only when liquid xfer is done.
  -- 9 canceled, 1 failed
  -- Note no need to decrement the cherry_pick_screening_count for the (9) canceled
  ---- because the WVA's are removed for these by SS1 process.
**/
create temp table temp_well_volume_adjustment_summary as (
  select 
  c.copy_id,
  w.well_id,
  p.plate_id,
  p.well_volume,
  p.well_volume + sum(wva.volume) as well_remaining_volume,
  count(wva) as cherry_pick_screening_count
  from well w
  join plate p using(plate_number)
  join copy c using(copy_id)
  join library l on(c.library_id=l.library_id),
  well_volume_adjustment wva
  where wva.copy_id=c.copy_id and wva.well_id = w.well_id
  group by c.copy_id, p.plate_id, w.well_id, p.well_volume
  order by w.well_id);

insert into copy_well
    ( id, copy_id, plate_id, well_id, initial_volume,volume, cherry_pick_screening_count)  
select 
  nextval('copy_well_id_seq'), 
  copy_id, 
  plate_id,
  well_id, 
  well_volume,
  well_remaining_volume,
  cherry_pick_screening_count
from temp_well_volume_adjustment_summary;

/*
  2. Adjust library well canonical concentrations to the current screening copy
    concentrations.
  NOTE: ICCBL Specific: 
    - Set the concentration for screening as the "canonical" concentration 
      on the wells.   
    - Stock copies will be adjusted by reversing any "dilution_factor" ( e.g.
      if the dilution was 2, then the new stock concentration will be the 
      original well concentration / 2 ).
    - Screening and Cherry Pick Request copies will have no concentration set.
      (and dilution factor is not used any more).
    - Any library with heterogeneous well concentrations AND dilution factors
      will need to have copy wells created for the stock copies that were 
      diluted (only XuTan7).
*/

/** 
  Store the original well concentrations -
**/
create table original_well_concentrations as (
  select library_id, well_id, mg_ml_concentration, molar_concentration 
  from well
  where (mg_ml_concentration is not null or molar_concentration is not null) );
/** 
  Store the original copy dilution factors -  
**/
create table original_copy_dilution_factors as (
  select copy_id, library_id, name, well_concentration_dilution_factor from copy );
  
/**
  2.a Use the (min) screening copy dilution factor 
    to adjust the well canonical concentration to that of the (max) screening 
    concentration.
*/
/* 
  For reference: copies with multiple screening concentrations:
    - these should be cleaned up before final migration, for now choose the 
    max concentration (minimum dilution) and use that as the new copy screening
    concentration.
    20180430:
 short_name  |  screen_type   |   df   |           copies           |        usage_type         |     molar      | molar_scaled | mgml  | mgml_scaled 
-------------+----------------+--------+----------------------------+---------------------------+----------------+--------------+-------+-------------
 ICCBBio1    | small_molecule |  1.000 | N;A;O;E;B;D;T;C            | library_screening_plates  |                |              | 5.000 |        5.00
 ICCBBio1    | small_molecule |  4.550 | H;P, R;P;G;F;Q;I           | library_screening_plates  |                |              | 5.000 |        1.10
 ICCBBio1    | small_molecule | 16.670 | J;L;M;K                    | library_screening_plates  |                |              | 5.000 |        0.30
 ICCBBio1    | small_molecule | 20.000 | R;S                        | library_screening_plates  |                |              | 5.000 |        0.25
 BiomolICCB1 | small_molecule |  1.000 | HA;C;A                     | cherry_pick_source_plates |                |              | 5.000 |        5.00
 BiomolICCB1 | small_molecule |  4.550 | MA                         | cherry_pick_source_plates |                |              | 5.000 |        1.10
 BiomolICCB1 | small_molecule |  1.000 | HF;HI;HC;HH;HG;HB;HE;HD;HJ | library_screening_plates  |                |              | 5.000 |        5.00
 BiomolICCB1 | small_molecule |  4.550 | MF;MG;MD;ME;MB;MC          | library_screening_plates  |                |              | 5.000 |        1.10
 BiomolICCB1 | small_molecule |  1.000 | H St A                     | stock_plates              |                |              | 5.000 |        5.00
 BiomolICCB1 | small_molecule |  4.550 | MSA                        | stock_plates              |                |              | 5.000 |        1.10
 HumanDUbs1  | rnai           |  1.000 | K;L                        | library_screening_plates  | 0.000010000000 |    0.0000100 |       |            
 HumanDUbs1  | rnai           | 10.000 | J;D;E;F;I;H;G;C            | library_screening_plates  | 0.000010000000 |    0.0000010 |       |            
 HumanDUbs1  | rnai           |  1.000 | Stock A;Stock B            | stock_plates              | 0.000010000000 |    0.0000100 |       |            
(13 rows)
    
**/
create table copies_with_multiple_screening_concentrations as (
  select
  short_name, screen_type,
  well_concentration_dilution_factor df,
  array_to_string(array_agg(name),';') as copies, 
  usage_type,
  (select array_to_string(array_agg(distinct(well.molar_concentration)),';') from well 
      where library_id = library.library_id ) molar,
  (select round((well.molar_concentration/well_concentration_dilution_factor),7) from well 
      where molar_concentration is not null and library.library_id = library_id limit 1) molar_scaled,
  (select array_to_string(array_agg(distinct(well.mg_ml_concentration)),';') from well 
      where library_id = library.library_id) mgml,
  (select round(well.mg_ml_concentration/well_concentration_dilution_factor,2) from well 
      where mg_ml_concentration is not null and library.library_id = library_id limit 1) mgml_scaled
  from copy
  join library using(library_id)
  join (
      select
      library_id
      from copy c join library using(library_id)
      where
      c.usage_type in ('cherry_pick_source_plates','library_screening_plates')
      group by library_id
      having count(distinct(well_concentration_dilution_factor))>1 
  ) lcs using(library_id)
  group by library.library_id, short_name, screen_type, usage_type, well_concentration_dilution_factor
);

/*
  Store the libraries with concentrations
*/
create temp table temp_libraries_with_concentrations as (
  select 
    library_id, short_name, 
    count(distinct(molar_concentration)) count_molar,
    count(distinct(mg_ml_concentration)) count_mg_ml
  from library join well using (library_id) 
  group by library_id, short_name 
  having count(distinct(molar_concentration)) > 1 or count(distinct(mg_ml_concentration)) > 1 
    or (count(distinct(molar_concentration)) > 0 and count(distinct(mg_ml_concentration)) > 0 )
  order by count_molar desc
);

/* 
  Save the concentration for homogeneous libraries 
**/

create temp table homogeneous_library_original_concentrations as (
  select
    library_id,
    max(mg_ml_concentration) as mg_ml_concentration,
    max(molar_concentration) as molar_concentration
  from library
  join well using (library_id)
  where library_id not in (
    select library_id from temp_libraries_with_concentrations)
  group by library_id);



/* 
  2.a Perform the canonical well update to the (max) screening concentration:
  - using min dilution because some libraries have multiple screening 
    concentrations (to be removed by prod time).
*/

create temp table screening_copy_min_dilutions as (
  select 
    library_id, min(well_concentration_dilution_factor) min_screening_dilution
  from copy
  where usage_type in ('library_screening_plates','cherry_pick_source_plates')
  group by library_id
);
create temp table screening_copy_max_dilutions as (
  select 
    library_id, max(well_concentration_dilution_factor) max_screening_dilution
  from copy
  where usage_type in ('library_screening_plates','cherry_pick_source_plates')
  group by library_id
);

/**
  NOTE: showing only 3 libraries where max/min are different (will be cleaned up)

select 
  short_name, 
  screen_type, 
  max_screening_dilution, min_screening_dilution,
  mg_ml_concentration, molar_concentration 
from homogeneous_library_original_concentrations 
join library using(library_id) 
join screening_copy_max_dilutions using(library_id)
join screening_copy_min_dilutions using(library_id)
where max_screening_dilution != min_screening_dilution;
 short_name  |  screen_type   | max_screening_dilution | min_screening_dilution | mg_ml_concentration | molar_concentration 
-------------+----------------+------------------------+------------------------+---------------------+---------------------
 HumanDUbs1  | rnai           |                 10.000 |                  1.000 |                     |      0.000010000000
 ICCBBio1    | small_molecule |                 20.000 |                  1.000 |               5.000 |                    
 BiomolICCB1 | small_molecule |                  4.550 |                  1.000 |               5.000 |                    
(3 rows)


  Other tests:

** shows 99 copies with min molar adjustments
  
select 
  short_name, 
  screen_type, 
  min_screening_dilution, 
  mg_ml_concentration, molar_concentration 
from homogeneous_library_original_concentrations 
join library using(library_id) 
join screening_copy_min_dilutions using(library_id) 
where molar_concentration is not null;

** shows 99 with max molar adjustments **

select 
  short_name, 
  screen_type, 
  max_screening_dilution, 
  mg_ml_concentration, molar_concentration 
from homogeneous_library_original_concentrations 
join library using(library_id) 
join screening_copy_max_dilutions using(library_id) 
where molar_concentration is not null;

** shows no records for (min) mg_ml concentration adjustments 
select 
  short_name, 
  screen_type, 
  min_screening_dilution, 
  mg_ml_concentration, molar_concentration 
from homogeneous_library_original_concentrations 
join library using(library_id) 
join screening_copy_min_dilutions using(library_id) 
where mg_ml_concentration is not null;

** shows only one for (max) mg_ml concentration adjustments - ICCBBio1 **
select 
  short_name, 
  screen_type, 
  max_screening_dilution, 
  mg_ml_concentration, molar_concentration 
from homogeneous_library_original_concentrations 
join library using(library_id) 
join screening_copy_max_dilutions using(library_id) 
where mg_ml_concentration is not null;


**/

create temp table new_well_concentrations as (
  select
    library_id, 
    well_id, 
    round(mg_ml_concentration/min_screening_dilution, 3) as mg_ml_concentration,
    round(molar_concentration/min_screening_dilution, 9) as molar_concentration
  from original_well_concentrations
  join screening_copy_min_dilutions using(library_id)
);

/** shows that no mg_ml adjustments will be made **/
select 
  count(*)
  from original_well_concentrations o 
  join new_well_concentrations n using(well_id) 
  where o.mg_ml_concentration is not null 
  and o.mg_ml_concentration != n.mg_ml_concentration;

/** shows 311,101 molar well adjustments (out of 629,232) **/
select 
  count(*)
  from original_well_concentrations o 
  join new_well_concentrations n using(well_id) 
  where o.molar_concentration is not null 
  and o.molar_concentration != n.molar_concentration;


update well set
    mg_ml_concentration = round(mg_ml_concentration/min_screening_dilution, 3),
    molar_concentration = round(molar_concentration/min_screening_dilution, 9)
  from screening_copy_min_dilutions 
  where screening_copy_min_dilutions.library_id = well.library_id
  and (well.mg_ml_concentration is not null or well.molar_concentration is not null);

/** check counts again **/
select 
  count(*)
  from original_well_concentrations o 
  join well n using(well_id) 
  where o.mg_ml_concentration is not null 
  and o.mg_ml_concentration != n.mg_ml_concentration;
select 
  count(*)
  from original_well_concentrations o 
  join well n using(well_id) 
  where o.molar_concentration is not null 
  and o.molar_concentration != n.molar_concentration;
/**
  2.b Create stock copy concentrations: 
  Find the stock copies with dilution factors and reverse this dilution factor
  so that the stock copies reflect the well concentration (for homogeneous libraries).
*/

/*
  2.b-1 Find the heterogeneous libraries with unique concentrations for wells:
*/

/*
  2.b-2 Set the new stock copy concentrations (homogeneous) for diluted stock
*/

create temp table new_stock_concentrations1 as (
  select library_id, copy_id,
    round(mg_ml_concentration/well_concentration_dilution_factor, 3) as mg_ml_concentration,
    round(molar_concentration/well_concentration_dilution_factor, 9) as molar_concentration
  from copy
  join homogeneous_library_original_concentrations using (library_id)
  where usage_type in ('stock_plates', '96_stock_plates')
  and well_concentration_dilution_factor !=1 and well_concentration_dilution_factor is not null
);

/** Shows (20180430) 17 copies/82 plates (one mg_ml) **/
select 
  short_name,
  hl.mg_ml_concentration, ns.mg_ml_concentration,
  hl.molar_concentration, ns.molar_concentration
  from homogeneous_library_original_concentrations hl 
  join library using (library_id)
  join new_stock_concentrations1 ns using(library_id);

update plate set
  mg_ml_concentration = nsc.mg_ml_concentration,
  molar_concentration = nsc.molar_concentration
  from new_stock_concentrations1 nsc
  where nsc.copy_id=plate.copy_id;

/**
  2.b-3 Set the new stock copy concentrations (homogeneous) for non-diluted stock
  - iff the screening_copy_min_dilution is != 1
**/
create temp table new_stock_concentrations2 as (
  select library_id, copy_id,
  hl.mg_ml_concentration,
  hl.molar_concentration
  from copy
  join screening_copy_min_dilutions scmin using(library_id)
  join homogeneous_library_original_concentrations hl using (library_id)
  where usage_type in ('stock_plates', '96_stock_plates')
  and scmin.min_screening_dilution != 1
  and well_concentration_dilution_factor=1 or well_concentration_dilution_factor is null
);

/** using SS1 data (imported), this shows (correctly) that there will be no 
  change from the ss1 values for these after migration ***

select
  copy.library_id,
  c.copy_id,
  copy.name,
  short_name,
  c.primary_well_mg_ml_concentration/copy.well_concentration_dilution_factor, ns.mg_ml_concentration,
  c.primary_well_molar_concentration/copy.well_concentration_dilution_factor, ns.molar_concentration
  from original_copy_concentrations c
  join copy using(copy_id) 
  join new_stock_concentrations1 ns using(copy_id)
  join library on copy.library_id=library.library_id
  where ( 
    round(primary_well_molar_concentration/copy.well_concentration_dilution_factor,9) != ns.molar_concentration
  or round(primary_well_mg_ml_concentration/copy.well_concentration_dilution_factor,3) != ns.mg_ml_concentration )
  and copy.usage_type in ('stock_plates', '96_stock_plates');
 
select 
  short_name,
  c.primary_well_mg_ml_concentration, ns.mg_ml_concentration,
  c.primary_well_molar_concentration, ns.molar_concentration
  from original_copy_concentrations c
  join copy using(copy_id) 
  join library using (library_id)
  join new_stock_concentrations2 ns using(library_id)
  where primary_well_mg_ml_concentration != ns.mg_ml_concentration;

select 
  short_name,
  c.primary_well_mg_ml_concentration, ns.mg_ml_concentration,
  c.primary_well_molar_concentration, ns.molar_concentration
  from original_copy_concentrations c
  join copy using(copy_id) 
  join library using (library_id)
  join new_stock_concentrations2 ns using(library_id)
  where primary_well_molar_concentration != ns.molar_concentration;


** are there any plates left after this that should have concentration set?
** this (correctly) shows only XuTan7; 

select 
  short_name,copy.name,
  copy_id, well_concentration_dilution_factor
  from copy
  join library using (library_id)
  where copy_id not in (select copy_id from new_stock_concentrations1) 
  and copy_id not in (select copy_id from new_stock_concentrations2)
  and well_concentration_dilution_factor != 1
  and usage_type in ('stock_plates', '96_stock_plates');

**/
/** Shows (20180430) 79 copies/ 2326 plates (no mg_ml) **/
select 
  short_name,
  hl.mg_ml_concentration, ns.mg_ml_concentration,
  hl.molar_concentration, ns.molar_concentration
  from homogeneous_library_original_concentrations hl 
  join library using (library_id)
  join new_stock_concentrations2 ns using(library_id);

update plate set
  mg_ml_concentration = nsc.mg_ml_concentration,
  molar_concentration = nsc.molar_concentration
  from new_stock_concentrations2 nsc
  where nsc.copy_id=plate.copy_id;
  
/**
  2.d Set the stock copy-well concentrations for non-homogeneous libraries:
    - create copy wells on the stock copies that match the original well / dilution
    concentration: NOTE: this is only XuTan7 for ICCB-L
**/

create temp table temp_copywell_concentrations as (
  select
    copy_id, plate_id, well_id,
    round(oc.mg_ml_concentration/well_concentration_dilution_factor, 3) as mg_ml_concentration,
    round(oc.molar_concentration/well_concentration_dilution_factor, 9) as molar_concentration
    from original_well_concentrations oc
    join temp_libraries_with_concentrations using(library_id)
    join copy using(library_id) 
    join plate using (copy_id)
    where well_concentration_dilution_factor != 1
    and (oc.mg_ml_concentration is not null or oc.molar_concentration is not null)
);

insert into copy_well
  ( id, copy_id, plate_id, well_id, mg_ml_concentration, molar_concentration )
  select 
    nextval('copy_well_id_seq'), 
    copy_id, 
    plate_id,
    well_id, 
    mg_ml_concentration, 
    molar_concentration
  from temp_copywell_concentrations tw
  where not exists (
    select null from temp_well_volume_adjustment_summary twvas
    where twvas.copy_id=tw.copy_id and twvas.well_id=tw.well_id);

/**
  2b. Update copy_well.mg_ml_concentration, copy_well.molar_concentration:
  - with temp_well_concentrations if they do exist.
**/

update copy_well
  set mg_ml_concentration = tw.mg_ml_concentration,
      molar_concentration = tw.molar_concentration
  from temp_copywell_concentrations tw 
  where copy_well.copy_id = tw.copy_id
  and copy_well.plate_id = tw.plate_id
  and copy_well.well_id = tw.well_id;

/*
  3a. Set plate.remaining_well_volume:
  - plate.remaining_well_volume == plate.well_volume - sum(library_screening.vol_xfer)
*/

select 'set plate remaining volume for library screening plates' as action;

create temp table temp_plate_remaining_volume as (
  select
    p.copy_id, p.plate_id,
    p.well_volume - sum(la.volume_transferred_per_well_from_library_plates) as plate_remaining_volume,
    count(la.*) as adjustments /** note this should be the same as screening_count **/
  from plate p
  join copy c using(copy_id)
  left join assay_plate ap using(plate_id)
  left join screening ls on(ls.activity_id = ap.library_screening_id)
  left join lab_activity la using(activity_id)
  where ap.replicate_ordinal = 0
  group by p.plate_id, p.copy_id, p.well_volume order by p.plate_id, copy_id 
);

update plate
  set remaining_well_volume = temp_plate_remaining_volume.plate_remaining_volume
  from temp_plate_remaining_volume 
  where temp_plate_remaining_volume.plate_id = plate.plate_id;


/** 
  3.c Set plate.remaining_well_volume = plate.well_volume (initial volume) for anything left
-- this shouldn't happen?
**/

select 'plates with remaining vol unset', count(*) 
  from plate 
  where remaining_well_volume is null;
/**  and avg_remaining_volume is null; **/

update plate set remaining_well_volume = well_volume 
  where remaining_well_volume is null;
/**  and avg_remaining_volume is null; **/

/** 
  3.d Set plate.experimental well count 
  - will be used to calculate avg volumes 
**/
update plate 
  set experimental_well_count = (
    select count(*) from well 
    where well.plate_number = plate.plate_number 
    and well.library_well_type='experimental');

/** 
  3.e Set plate.screening_count
**/

select 'default plate screening statistics ' as action;

create temp table plate_screening_count as (
  select 
  p.plate_id, 
  p.plate_number, 
  count(distinct(ls))  
  from plate p 
  join assay_plate ap using(plate_id) 
  join library_screening ls on(ap.library_screening_id=ls.activity_id) 
  where ap.replicate_ordinal=0 
  group by p.plate_id, p.plate_number 
  order by plate_id );

update plate as p 
  set screening_count = psc.count
  from plate_screening_count psc 
  where p.plate_id=psc.plate_id;
  
/**
  3.f Add to screening_count for cherry_pick_screenings
  
  - find wvas for cp
  - group by plate
  - count distinct cplt's per plate
  NOTE: could also count the distinct cherry_pick_requests, but this is wrong
  because multiple liquid transfer actions may be performed.
  
  
**/  
select 'Update cplt_based_screening_count' as action;

create temp table cplt_based_screening_count as
  select
  count(distinct(cplt.*)) cplt_screening_count,
  plate.plate_number,
  plate.plate_id,
  copy.name
  from 
  well_volume_adjustment wva
  join copy using(copy_id)
  join well using(well_id)
  join plate on(well.plate_number=plate.plate_number and plate.copy_id=wva.copy_id)
  join lab_cherry_pick using(lab_cherry_pick_id)
  join cherry_pick_assay_plate using(cherry_pick_assay_plate_id)
  join cherry_pick_liquid_transfer cplt on (cherry_pick_liquid_transfer_id=cplt.activity_id)
  where cplt.status != 'Canceled' 
  group by plate.plate_number, plate.plate_id, copy.name;

update plate
  set cplt_screening_count = cplt_based_screening_count.cplt_screening_count
  from cplt_based_screening_count
  where cplt_based_screening_count.plate_id = plate.plate_id;

select 'Update plate.date_plated, date_retired' as action;
update plate
  set date_plated=activity.date_of_activity
  from activity 
  where plate.plated_activity_id=activity.activity_id;

update plate
  set date_retired=activity.date_of_activity
  from activity 
  where plate.retired_activity_id=activity.activity_id;

/** 
 Try to close the transaction, because otherwise postgres gives errors:
 cannot ALTER TABLE "copy_well" because it has pending trigger events
**/
commit;

select 'create indexes and cleanup ' as action;

/** create indexes not already there **/
create index assay_plate_plate_id on assay_plate (plate_id);
create index copy_library_id on copy (library_id);
create index plate_copy_id on plate (copy_id);
create index well_library_id on well (library_id);

/** indexes for the result_value and assay_well queries 
**/
create index assay_well_sr_is_positive_idx2 on assay_well(screen_result_id, is_positive, well_id);

create index result_value_is_excluded_index on result_value(data_column_id, is_exclude);

/** add back in fk constraint **/

alter table copy_well 
  add constraint well_id_refs_well_id 
  foreign key (well_id) 
  references well;

/**
 TODO: 20180404 - Verify that all indexes match to 0001_initial and models.py
alter table assay_well
  add constraint fk_assay_well_to_well
  foreign key (well_id)
  references well;
alter table assay_well
  add constraint fk_assay_well_screen_result
  foreign key (screen_result_id)
  references screen_result;

create index assay_well_well_id on assay_well(well_id);

create index assay_well_confirmed_positives_data_only_index
  btree (confirmed_positive_value) 
  WHERE confirmed_positive_value IS NOT NULL;

alter table result_value
  add constraint fk_result_value_to_data_column
  FOREIGN KEY (data_column_id) REFERENCES data_column(data_column_id);

alter table result_value
  add constraint fk_result_value_to_well
  FOREIGN KEY (well_id) REFERENCES well(well_id);

## not shown to help: create index well_id_idx on result_value(well_id);
  
**/  

/**
vacuum analyze;
**/
select current_time;

