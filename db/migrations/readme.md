######## Migrations for an existing Screensaver installation  #####

meta:

-- Note: existing migrations are intended for migrating an existing screensaver 
-- database to the new database schema.
-- Therefore, these should not be used to create a new instance.
-- instead, the migrations directory should be removed and recreated with 
-- makemigrations, which will generate an initial migration by inspecting the 
-- models.py file, and then manage.py migrate --fake-initial should be run.


-- Recreating the migration initial file (for clean installs of the app).
-- New Django >= 1.7

./manage.py makemigrations db
./manage.py migrate --fake-initial
-- OLD: Django < 1.7 
-- ./manage.py schemamigration db --initial # regenerate 001 script
-- # now if the db has changed, can generate a migration script for it:
-- ./manage.py schemamigration db --auto


0001 - initial

- should be faked, if the datbase is already existent, or,
-- if new project, run only this migration. 
-- Note: this migration represents the legacy Screensaver schema; all updates
-- to the schema are contained in subsequent migrations.
./manage.py migrate db 0001 --fake  


0002_screen_status:

Creates a status field on on Screen; this will replace the screen_status table
./manage.py migrate db 0002

* After performing the schema migration, use the manual script to populate the 
Screen table with the latest status.

* Purpose: add an id to capture the natural ordering of the status table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.


psql -Uuser database -f ./migrations/manual/0003_screen_status.sql

0003_screen_status:

**NOTE: "manual/0003_screen_status.sql" must be run before 
"0003_screen_status.py" is run, in order to convert the table to have an id 
field so that the Django ORM can read it.


migrate the screen_status_item table to the "status" field of the screen object 
and create the needed ApiLog entries to record the history

0004_auto__chg_field_administratoruser_screensaver_user__chg_field_labhead

* change screensaver_user fields on administrator_user, lab_head, 
screening_room_user to onetoone
* change screen field on screen_result to OneToOne
* add screensaver_user field to django auth.user

0005_create_django_users

* create an entry in the django auth.user table for every screensaver_user that
has email and:
1) ecommons_id: username=ecommons
2) login_id: username=login_id
3) otherwise, skip

0006_screensaver_user_permissions:

* create a m2m table joining permissions to screensaver_user

0008_library_controlled_vocabularies

Performs the vocabulary conversion for the library and screen entities.
        - we are converting from the stored, display values, to a key value;
        this key will reference the entry for the respective vocab in the 
        reports_vocabularies table.  (Referencing by scope, key).
        Note: there will be a conversion script for all vocabularies, eventually.
        
* convert the following vocabs on library: 
['library_type', 'solvent', 'screen_type', 'screening_status']
from title form to key form (replaces non-characters with underscores, 
convert the whole string to lowercase).

* convert the following vocabs on Screen: 
['screen_type', 'project_phase', 'species' ]

* note, this conversion should be updated if new vocabs need to be converted

0009_create_autofields_screen_library

Re-associate the sequences for the screen_id and library_id fields.
        - set the default nextval on the column
        - set ownership on the sequence to the table
        
        Note: we will be doing this for _every_ sequence in the database.
        (Old method in Hibernate used a Java based Generator class methodology
        to associate).

        
0010_auto__add_field_reagent_substance_id:

* create reagent.substance_id 

0011_reagent_substance_ids:

* set reagent substance IDs

0012_create_reagent_autofield

* setup autofield reagent_id

0013_create_library_content_diffs

** NOTE: before running this migration the manual migration 
0003_screen_status.sql must be run first.  This script converts the gene_symbol
table to have an internal pk ID field; this makes the table readable by Django.

This migration converts the data tracked in the "library_content_version" (LCV)
table into equivalent report_apilog table entries.  For each well in a library,
each reagent version is compared; diff_keys and diffs entries are created on 
the Apilog.
