BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
7155,
current_timestamp,
'migrate screen_status_item into the api_history and screen';

/* add an id to capture the natural ordering of the table (missing an id field, this also will allow the Django ORM to read it */

alter table screen_status_item add column id integer;
create sequence screen_status_item_sequence;
update screen_status_item set id=nextval('screen_status_item_sequence');

update screen set status = si.status 
    from screen_status_item si 
    join screen s on(si.screen_id=s.screen_id) 
    where  (si.id = (select max(si1.id) from screen_status_item si1 where si1.screen_id = s.screen_id limit 1) );

update screen set status_date = si.status_date 
    from screen_status_item si 
    join screen s on(si.screen_id=s.screen_id) 
    where  (si.id = (select max(si1.id) from screen_status_item si1 where si1.screen_id = s.screen_id limit 1) );

/** TODO this migration must be run 1st, then the following is performed with South **/
/** insert into reports_apilog ( user_id, username, ref_resource_name, key, uri, date_time, api_action, diff_keys, diffs) **/
/** select 1, 'sde4', 'screen', s.facility_id, '/db/api/v1/screen/'||s.facility_id, '["status","status_date"]', ' **/


COMMIT;
