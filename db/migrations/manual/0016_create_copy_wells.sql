BEGIN;


INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20161012,
current_timestamp,
'create new copywell data using assay plate information';

/**
  Create copywell entries as needed for:
  - well volume adjustments
**/

/** 
- 1. create a copy well for all well volume adjustments
- current plate.well_volume - adjustments
- # of adjustments
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
    ( id, copy_id, plate_id, well_id, initial_volume,volume,adjustments)  
select 
  nextval('copy_well_id_seq'), 
  copy_id, 
  plate_id,
  well_id, 
  well_volume,
  well_remaining_volume,
  adjustments
from temp_well_volume_adjustment_summary;

/*
- 2. create a temp table with wells with unique concentration values
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
- 2a. insert temp_well_concentrations into copy_well if they dne (leaving volume null, 
  since the volume will default to the plate volume if not set, adjustments to the screening_count),
- 2b. update copy_well with temp_well_concentrations if the do exist.
*/
insert into copy_well
    ( id, copy_id, plate_id, well_id, mg_ml_concentration, molar_concentration, adjustments)  
select 
  nextval('copy_well_id_seq'), 
  copy_id, 
  plate_id,
  well_id, 
  mg_ml_concentration, 
  molar_concentration,
  0
from temp_well_concentrations tw
where not exists (
  select null from temp_well_volume_adjustment_summary twvas
  where twvas.copy_id=tw.copy_id and twvas.well_id=tw.well_id);

update copy_well
  set mg_ml_concentration = tw.mg_ml_concentration,
  molar_concentration = tw.molar_concentration
from temp_well_concentrations tw
where tw.copy_id=copy_well.copy_id and tw.well_id=copy_well.well_id;

/*
- 3a. create plate initial and remaining volumes for all plates
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
  set remaining_volume = pv.plate_remaining_volume,
  avg_remaining_volume = pv.plate_remaining_volume,
  min_remaining_volume = pv.plate_remaining_volume,
  max_remaining_volume = pv.plate_remaining_volume
  from temp_plate_remaining_volume pv 
  where pv.plate_id = plate.plate_id;

/*
- 3b. create plate concentrations for all libraries not covered by the wells in 2
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
- 4. set the plate statistics 
  (**TODO: this goes away if stats can be done on the fly**)
  (**Use as intermediate check step for 3.c, below)
**/
select 'create plate volume statistics for copy_wells for cherry_picks' as action;

create temp table plate_volume2 as (
  select
  cw.plate_id,
  round(avg(cw.volume)::NUMERIC,9) avg_remaining_volume, /** 9 dec pl = nano liter **/
  min(cw.volume) min_remaining_volume,
  max(cw.volume) max_remaining_volume
  from copy_well cw
  group by cw.plate_id );

update plate set 
  avg_remaining_volume = pv.avg_remaining_volume, 
  min_remaining_volume = pv.min_remaining_volume, 
  max_remaining_volume = pv.max_remaining_volume
  from plate_volume2 pv 
  where plate.plate_id = pv.plate_id;

/** 
- 3.c Set remaining_volume = well_volume (initial volume) for anything left
-- this shouldn't happen?
**/

select 'plates with remaining vol unset', count(*) 
  from plate 
  where remaining_volume is null 
  and avg_remaining_volume is null;

update plate set remaining_volume = well_volume 
  where remaining_volume is null 
  and avg_remaining_volume is null;

/** 
  Done - volume statistics
  Next - plate screening statistics
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

/** TODO: copy screening statistics - may wish to store on copy 

**/
select 'copy screening statistics ' as action;

insert into copy_screening_statistics
( SELECT c.copy_id,
    c.name,
    l.short_name
    ,count(distinct(ls)) screening_count 
    ,count(distinct(ap)) ap_count 
    ,count(distinct(srdl)) dl_count
    ,min(srdl.date_of_activity) first_date_data_loaded 
    ,max(srdl.date_of_activity) last_date_data_loaded 
    ,min(lsa.date_of_activity) first_date_screened 
    ,max(lsa.date_of_activity) last_date_screened 
    FROM copy c 
    join plate p using(copy_id) 
    join library l using(library_id) 
    join assay_plate ap on(ap.plate_id=p.plate_id)  
        left outer join library_screening ls on(library_screening_id=ls.activity_id)
        left outer join activity lsa on(ls.activity_id=lsa.activity_id) 
    left outer join (
      select a.activity_id, a.date_of_activity 
      from activity a 
      join administrative_activity using(activity_id) ) srdl 
          on(screen_result_data_loading_id=srdl.activity_id)  
    group by c.copy_id, c.name, l.short_name);
/**
 ** REMOVED - calculate dynamically
insert into plate_screening_statistics
( SELECT
    p.plate_id, 
    p.plate_number,
    c.copy_id,
    c.name,
    l.short_name,
    l.library_id,
    count(distinct(ls)) screening_count, 
    count(distinct(ap)) ap_count ,
    count(distinct(srdl)) dl_count,
    min(srdl.date_of_activity) first_date_data_loaded, 
    max(srdl.date_of_activity) last_date_data_loaded, 
    min(lsa.date_of_activity) first_date_screened, 
    max(lsa.date_of_activity) last_date_screened 
    FROM copy c 
    join plate p using(copy_id) 
    join library l using(library_id) 
    join assay_plate ap on(ap.plate_id=p.plate_id)  
      left outer join library_screening ls on(library_screening_id=ls.activity_id)
            left outer join activity lsa on(ls.activity_id=lsa.activity_id) 
    left outer join (
      select a.activity_id, a.date_of_activity 
      from activity a 
      join administrative_activity using(activity_id) ) srdl 
    on(screen_result_data_loading_id=srdl.activity_id)  
    group by p.plate_id, p.plate_number, c.copy_id, c.name, l.library_id, l.short_name );
**
**/


select 'create indexes and cleanup ' as action;

/** create indexes not already there **/
create index assay_plate_plate_id on assay_plate (plate_id);
create index copy_library_id on copy (library_id);
create index plate_copy_id on plate (copy_id);
create index well_library_id on well (library_id);

/** indexes for the result_value and assay_well queries **/
create index assay_well_sr_is_positive_idx2 on assay_well(screen_result_id, is_positive, well_id);


/** add back in fk constraint 
alter table copy_well 
  add constraint well_id_refs_well_id 
  foreign key (well_id) 
  references well;
**/

COMMIT;

/**
vacuum analyze;
**/
select current_time;

