/**
  Clean up the result value table:
  
  - remove duplicate result values located for certain data columns
  - unclear why these exist
  - Screensaver 1 arbitrarily uses the last value from these duplicates, when
  using the default sort order of the result_value table. Presumably these last
  values are the last values that were entered into screensaver for these records. 

**/

/**
 NOTE: result_value_duplicates has already been created on ss1 database (20180529) 
 
drop table result_value_duplicates;
create table result_value_duplicates as 
  select data_column_id, well_id, count(*)
  from result_value
  group by data_column_id, well_id
  HAVING count(*) > 1;
**/


drop table result_value_duplicates2;
create table result_value_duplicates2 as 
  select rv.*  
  from result_value rv 
  join result_value_duplicates rvd on (rv.well_id=rvd.well_id and rvd.data_column_id=rv.data_column_id); 

create sequence result_value_duplicates_seq;

alter table result_value_duplicates2 add column sequence integer;

update result_value_duplicates2 set sequence = nextval('result_value_duplicates_seq');


delete from result_value where result_value_id in (
  select result_value_id from result_value_duplicates2 where sequence not in (
    select seq from (
      select max(sequence) as seq, well_id, data_column_id
      from result_value_duplicates2
      group by well_id, data_column_id 
    ) a 
  )
);

alter table result_value 
  ADD constraint data_column_well_id unique(data_column_id, well_id);

/**
drop index result_value_data_column_and_well_index;
**/