/**
  Refactor the result_value table to be more efficient:
  
  - divide small_molecule and rnai 

  # PRELIMINARY 2017-12-18 - May not be necessary with result_value_cleanup.sql
**/

create table result_value_sm as (
  select * from result_value rv
  where exists(select null from data_column
    join screen_result using(screen_result_id)
    join screen using(screen_id)
    where data_column.data_column_id = rv.data_column_id
    and screen.screen_type = 'small_molecule') );

alter table result_value_sm add column id integer;

create sequence result_value_sm_sequence;

update result_value_sm set id=nextval('result_value_sm_sequence')
  order by result_value_sm.result_value_id asc;
  
create index result_value_sm_data_column_and_numeric_value_index 
  on result_value_sm(data_column_id, numeric_value);
create index result_value_sm_data_column_and_value_index 
  on result_value_sm(data_column_id, value);
create index result_value_sm_data_column_and_well_index 
  on result_value_sm(data_column_id, well_id);
/** new **/
create index result_value_sm_is_excluded 
  on result_value_sm(well_id, data_column_id, is_exclude);

alter table result_value_sm 
  alter column id 
  set default nextval('result_value_sm_sequence'::regclass);

alter table result_value_sm drop column result_value_id;

/**
 RNAi
**/
create table result_value_rnai as (
  select * from result_value rv
  where exists(select null from data_column
    join screen_result using(screen_result_id)
    join screen using(screen_id)
    where data_column.data_column_id = rv.data_column_id
    and screen.screen_type = 'rnai') );

alter table result_value_rnai add column id integer;

create sequence result_value_rnai_sequence;

update result_value_rnai set id=nextval('result_value_rnai_sequence')
  order by result_value_rnai.result_value_id asc;
    
create index result_value_rnai_data_column_and_numeric_value_index 
  on result_value_rnai(data_column_id, numeric_value);
create index result_value_rnai_data_column_and_value_index 
  on result_value_rnai(data_column_id, value);
create index result_value_rnai_data_column_and_well_index 
  on result_value_rnai(data_column_id, well_id);
/** new **/
create index result_value_rnai_is_excluded 
  on result_value_rnai(well_id, data_column_id, is_exclude);

