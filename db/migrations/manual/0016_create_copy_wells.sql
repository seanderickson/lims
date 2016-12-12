BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20161012,
current_timestamp,
'create new copywell data using assay plate information';

select current_time;

/**
  Create copywell entries as needed for:
  - well volume adjustments
**/

/** 
  1. Create copy_wells for well_volume_adjustments:
  - set copy_well.volume (remaining_volume) = current plate.well_volume - sum(wva's)
  - set # of adjustments = count(wva's)
  Note: all wva's should be used, not just some of successful, failed, canceled, 
  because wva is created only when liquid xfer is done.
**/
create temp table temp_well_volume_adjustment_summary as (
  select 
  c.copy_id,
  w.well_id,
  p.plate_id,
  p.well_volume,
  p.well_volume + sum(wva.volume) as well_remaining_volume,
  count(wva) as adjustments
  from well w
  join plate p using(plate_number)
  join copy c using(copy_id)
  join library l on(c.library_id=l.library_id),
  well_volume_adjustment wva
  where wva.copy_id=c.copy_id and wva.well_id = w.well_id
  group by c.copy_id, p.plate_id, w.well_id, p.well_volume
  order by w.well_id);

insert into copy_well
    ( id, copy_id, plate_id, well_id, initial_volume,volume)  
select 
  nextval('copy_well_id_seq'), 
  copy_id, 
  plate_id,
  well_id, 
  well_volume,
  well_remaining_volume
from temp_well_volume_adjustment_summary;

/*
  2. create copy_wells for wells having specific concentration values:
  - create a temp table with wells with unique concentration values
*/
create temp table temp_libraries_with_concentrations as (
  select 
    library_id, short_name, 
    count(distinct(molar_concentration)) 
  from library join well using (library_id) 
  group by library_id, short_name 
  having count(distinct(molar_concentration)) > 1
  order by count desc
);

create temp table temp_well_concentrations as (
  select
    copy_id,
    plate_id,
    well.plate_number,
    well_id,
    well_volume,
    well_volume as well_remaining_volume,
    mg_ml_concentration,
    molar_concentration
  from well 
  join temp_libraries_with_concentrations tl using(library_id)
  join (
  select plate_id, copy_id, library_id, plate_number, well_volume
  from plate join copy using(copy_id) ) cp on(well.plate_number=cp.plate_number and cp.library_id=tl.library_id)
  where well.library_well_type = 'experimental');

/*
  2a. Create/Set copy_well.mg_ml_concentration, copy_well.molar_concentration:
  - for concentration-specific wells if they don't exist -
  (leaving volume null, since the volume will default to the plate volume if
   not set),
*/
insert into copy_well
    ( id, copy_id, plate_id, well_id, mg_ml_concentration, molar_concentration)  
  select 
    nextval('copy_well_id_seq'), 
    copy_id, 
    plate_id,
    well_id, 
    mg_ml_concentration, 
    molar_concentration
  from temp_well_concentrations tw
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
  from temp_well_concentrations tw
  where tw.copy_id=copy_well.copy_id and tw.well_id=copy_well.well_id;

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

/*
  3b. Set plate concentration fields:
  - for all libraries having a single concentration value.
  -- NOTE: the min/max/primary concentration values will differ because they do 
  -- not reflect the multiplication by the well_concentration_dilution_factor.
  -- (verified: this operation matches current values).
*/

create temp table current_library_mg_ml_concentrations as (
  select 
  library_id, max(mg_ml_concentration) as mg_ml_concentration 
  from well 
  where library_well_type = 'experimental' 
  and mg_ml_concentration is not null 
  and not exists (
    select null from temp_libraries_with_concentrations 
    where library_id = well.library_id) 
  group by library_id );

create temp table current_library_molar_concentrations as (
  select 
  library_id, max(molar_concentration) as molar_concentration 
  from well 
  where library_well_type = 'experimental' 
  and molar_concentration is not null 
  and not exists (
    select null from temp_libraries_with_concentrations 
    where library_id = well.library_id) 
  group by library_id );

update plate
  set mg_ml_concentration 
    = round(l.mg_ml_concentration/copy.well_concentration_dilution_factor,9) 
    /** 9 dec places = nano liter **/
  from current_library_mg_ml_concentrations l
  join copy using(library_id)
  where plate.copy_id=copy.copy_id;

update plate
  set molar_concentration 
    = round(l.molar_concentration/copy.well_concentration_dilution_factor,9) 
    /** 9 dec places = nano liter **/
  from current_library_molar_concentrations l
  join copy using(library_id)
  where plate.copy_id=copy.copy_id;
 

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


/** add back in fk constraint **/

alter table copy_well 
  add constraint well_id_refs_well_id 
  foreign key (well_id) 
  references well;

/**
vacuum analyze;
**/
select current_time;

