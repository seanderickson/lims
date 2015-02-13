BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20150016,
current_timestamp,
'create new copywell data using assay plate information';

/**
  Purpose: create new copywell data using assay plate information
  *** preliminary ***
**/


select current_time;

insert into copy_well
    ( id, copy_id, plate_id, well_id, well_name, plate_number)  
    select nextval('copy_well_id_seq'), 
    c.copy_id, 
    p.plate_id,
    w.well_id, 
    w.well_name, 
    p.plate_number 
    from well w 
    left join plate p using (plate_number) 
    left join copy c using(copy_id)
    join library l on(c.library_id=l.library_id) 
    order by w.well_id;

create temp table plate_volume as (
    select
    p.copy_id, p.plate_id,
    p.well_volume - sum(la.volume_transferred_per_well_from_library_plates) as plate_remaining_volume
    from plate p
    join copy c using(copy_id)
    left join assay_plate ap using(plate_id)
    left join screening ls on(ls.activity_id = ap.library_screening_id)
    left join lab_activity la using(activity_id)
    where ap.replicate_ordinal = 0
    group by p.plate_id, p.copy_id, p.well_volume order by p.plate_id, copy_id );

update copy_well
set volume = pv.plate_remaining_volume
from plate_volume as pv
where pv.copy_id = copy_well.copy_id and pv.plate_id=copy_well.plate_id;

create temp table well_volume_adjustments as (
    select 
    c.copy_id,
    w.well_id,
    p.well_volume,
    p.well_volume + sum(wva.volume) as well_remaining_volume
    from well w
    join plate p using(plate_number)
    join copy c using(copy_id)
    join library l on(c.library_id=l.library_id),
    well_volume_adjustment wva
    where wva.copy_id=c.copy_id and wva.well_id = w.well_id
    group by c.copy_id, w.well_id, p.well_volume
    order by w.well_id);

update copy_well as cw
set volume = wva.well_remaining_volume
from well_volume_adjustments as wva
where cw.copy_id = wva.copy_id and cw.well_id = wva.well_id;

update copy_well cw set volume = p.well_volume 
from plate p
where p.copy_id = cw.copy_id and p.plate_id=cw.plate_id
and volume is null;

/** set the plate statistics **/

/** OLD ****
    with plate_volume as (
      select 
      p.plate_id,
      p.well_volume - sum(la.volume_transferred_per_well_from_library_plates) as plate_remaining_volume
      from plate p
      left join assay_plate ap using(plate_id) 
      left join screening ls on(ls.activity_id = ap.library_screening_id) 
      left join lab_activity la using(activity_id) 
      where ap.replicate_ordinal = 0 
      group by p.plate_id, p.well_volume order by p.plate_id )
    update plate as p
    set remaining_volume = pv.plate_remaining_volume
    from plate_volume as pv
    where pv.plate_id = p.plate_id;


**/


create temp table plate_volume2 as (
  select
  cw.plate_id,
  round(avg(cw.volume),9) avg_remaining_volume, /** 9 dec pl = nano liter **/
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


/** update plate set remaining_volume = well_volume where remaining_volume is null;
**/

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

/** TODO: copy screening statistics - may wish to store on copy **/

create table copy_screening_statistics ( 
    copy_id integer, 
    name text,
    short_name text,
    screening_count integer,
    ap_count integer,
    dl_count integer,
    first_date_data_loaded date,
    last_date_data_loaded date,
    first_date_screened date,
    last_date_screened date,
    primary key (copy_id)
);

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

create table plate_screening_statistics ( 
    plate_id integer,
    plate_number integer,
    copy_id integer,
    copy_name text,
    library_short_name text,
    library_id integer,
    screening_count integer,
    ap_count integer,
    dl_count integer,
    first_date_data_loaded date,
    last_date_data_loaded date,
    first_date_screened date,
    last_date_screened date,
    primary key (plate_id)
);

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



create index assay_plate_plate_id on assay_plate (plate_id);
create index copy_library_id on copy (library_id);
create index plate_copy_id on plate (copy_id);


COMMIT;

/**
vacuum analyze;
**/
select current_time;
