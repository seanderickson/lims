BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
9990004,
current_timestamp,
'migrate screen_update_activity';

/* add an id to capture the natural ordering of the table (missing an id field, this also will allow the Django ORM to read it */
/* because the table has no ID, Django ORM doesn't know what to do with it - so, ironically? the South migration cannot do it either */

alter table screen_update_activity add column id integer;
create sequence screen_update_activity_sequence;
update screen_update_activity set id=nextval('screen_update_activity_sequence');

COMMIT;
