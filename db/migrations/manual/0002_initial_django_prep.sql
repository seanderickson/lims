BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20140002,
current_timestamp,
'Initial Django prep: 
  migrate screen_status_item into the api_history and screen, 
  migrate gene_symbol table';

/**
  *** Initial manual migration ***
  
  Depends on migration script 0002_initial_django_prep:
  - the migration script creates fields used here, such as:
  status, status_date fields on Screen
  ** This script must be run before the migrations/0004_screen_status.py is run.
**/


/**
  *** Add id field to screen_status ***
  Purpose: add an id to capture the natural ordering of the status table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/

alter table screen_status_item add column id integer;
create sequence screen_status_item_sequence;
update screen_status_item set id=nextval('screen_status_item_sequence');

update screen set status = si.status 
    from screen_status_item si 
    join screen s on(si.screen_id=s.screen_id) 
    where  (si.id = (select max(si1.id) from screen_status_item si1 where si1.screen_id = screen.screen_id limit 1) );
    

update screen set status_date = si.status_date 
    from screen_status_item si 
    join screen s on(si.screen_id=s.screen_id) 
    where  (si.id = (select max(si1.id) from screen_status_item si1 where si1.screen_id = screen.screen_id limit 1) );

/** TODO this migration must be run 1st, then the following is performed with South **/
/** insert into reports_apilog ( user_id, username, ref_resource_name, key, uri, date_time, api_action, diff_keys, diffs) **/
/** select 1, 'sde4', 'screen', s.facility_id, '/db/api/v1/screen/'||s.facility_id, '["status","status_date"]', ' **/

/**
  *** Add id field to cherry_pick_request_empty_well  ***
  Purpose: add an id to capture the natural ordering of the cherry_pick_request_empty_well table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/


/**
  *** Add id field to gene_symbol ***
  Purpose: add an id to capture the natural ordering of the gene_symbol table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/

alter table gene_symbol add column id integer;
create sequence gene_symbol_sequence;
update gene_symbol set id=nextval('gene_symbol_sequence');

/** since this table will still be used, make the changes for the autofield **/
alter table gene_symbol drop CONSTRAINT gene_symbol_pkey;
alter table gene_symbol add primary key(id);
alter table gene_symbol ADD constraint gene_symbol_natural_key unique(gene_id,ordinal);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table gene_symbol alter column id set default nextval('gene_symbol_sequence'::regclass);
alter sequence gene_symbol_sequence owned by gene_symbol.id;

/**
  *** Add id field to gene_genbank_accession_number ***
  Purpose: add an id to gene_genbank_accession_number table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/

alter table gene_genbank_accession_number add column id integer;
create sequence gene_genbank_accession_number_sequence;
update gene_genbank_accession_number set id=nextval('gene_genbank_accession_number_sequence');
/** since this table will still be used, make the changes for the autofield **/
alter table gene_genbank_accession_number drop CONSTRAINT gene_genbank_accession_number_pkey1;
alter table gene_genbank_accession_number add primary key(id);
alter table gene_genbank_accession_number 
    ADD constraint gene_genbank_accession_number_natural_key unique(gene_id,genbank_accession_number);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table gene_genbank_accession_number alter column id set default nextval('gene_genbank_accession_number_sequence'::regclass);
alter sequence gene_genbank_accession_number_sequence owned by gene_genbank_accession_number.id;


/**
  *** Add id field to small_molecule sub-tables***
  Purpose: add an id;
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/

\set _table small_molecule_pubchem_cid
\set _unique_constraint reagent_id,pubchem_cid

\set _unique_constraint_name :_table _unique
\set _table_sequence :_table _sequence
\set _table_pkey :_table _pkey
\set _quoted_sequence '\'' :_table_sequence '\''

alter table :_table drop constraint :_table_pkey;
alter table :_table add constraint :_unique_constraint_name UNIQUE(:_unique_constraint);
alter table :_table add column id integer;
create sequence :_table_sequence;
update :_table set id=nextval(:_quoted_sequence);
/** since this table will still be used, make the changes for the autofield **/
alter table :_table add primary key(id);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table :_table alter column id set default nextval(:_quoted_sequence::regclass);
alter sequence :_table_sequence owned by :_table.id;


\set _table small_molecule_chembl_id
\set _unique_constraint reagent_id,chembl_id

\set _unique_constraint_name :_table _unique
\set _table_sequence :_table _sequence
\set _table_pkey :_table _pkey
\set _quoted_sequence '\'' :_table_sequence '\''

alter table :_table drop constraint :_table_pkey;
alter table :_table add constraint :_unique_constraint_name UNIQUE(:_unique_constraint);
alter table :_table add column id integer;
create sequence :_table_sequence;
update :_table set id=nextval(:_quoted_sequence);
/** since this table will still be used, make the changes for the autofield **/
alter table :_table add primary key(id);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table :_table alter column id set default nextval(:_quoted_sequence::regclass);
alter sequence :_table_sequence owned by :_table.id;


\set _table small_molecule_chembank_id
\set _unique_constraint reagent_id,chembank_id

\set _unique_constraint_name :_table _unique
\set _table_sequence :_table _sequence
\set _table_pkey :_table _pkey
\set _quoted_sequence '\'' :_table_sequence '\''

alter table :_table drop constraint :_table_pkey;
alter table :_table add constraint :_unique_constraint_name UNIQUE(:_unique_constraint);
alter table :_table add column id integer;
create sequence :_table_sequence;
update :_table set id=nextval(:_quoted_sequence);
/** since this table will still be used, make the changes for the autofield **/
alter table :_table add primary key(id);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table :_table alter column id set default nextval(:_quoted_sequence::regclass);
alter sequence :_table_sequence owned by :_table.id;


\set _table small_molecule_compound_name
\set _unique_constraint reagent_id,ordinal

\set _unique_constraint_name :_table _unique
\set _table_sequence :_table _sequence
\set _table_pkey :_table _pkey
\set _quoted_sequence '\'' :_table_sequence '\''

alter table :_table drop constraint :_table_pkey;
alter table :_table add constraint :_unique_constraint_name UNIQUE(:_unique_constraint);
alter table :_table add column id integer;
create sequence :_table_sequence;
update :_table set id=nextval(:_quoted_sequence);
/** since this table will still be used, make the changes for the autofield **/
alter table :_table add primary key(id);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table :_table alter column id set default nextval(:_quoted_sequence::regclass);
alter sequence :_table_sequence owned by :_table.id;



\set _table screen_keyword
\set _unique_constraint screen_id,keyword

\set _unique_constraint_name :_table _unique
\set _table_sequence :_table _sequence
\set _table_pkey :_table _pkey
\set _quoted_sequence '\'' :_table_sequence '\''

alter table :_table drop constraint :_table_pkey;
alter table :_table add constraint :_unique_constraint_name UNIQUE(:_unique_constraint);
alter table :_table add column id integer;
create sequence :_table_sequence;
update :_table set id=nextval(:_quoted_sequence);
/** since this table will still be used, make the changes for the autofield **/
alter table :_table add primary key(id);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table :_table alter column id set default nextval(:_quoted_sequence::regclass);
alter sequence :_table_sequence owned by :_table.id;



/**
  *** Add id field to silencing_reagent_duplex_wells ***
  Purpose: add an id to silencing_reagent_duplex_wells table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/

alter table silencing_reagent_duplex_wells add column id integer;
create sequence silencing_reagent_duplex_wells_sequence;
update silencing_reagent_duplex_wells 
    set id=nextval('silencing_reagent_duplex_wells_sequence');
/** since this table will still be used, make the changes for the autofield **/
alter table silencing_reagent_duplex_wells 
    drop CONSTRAINT silencing_reagent_duplex_wells_pkey;
alter table silencing_reagent_duplex_wells add primary key(id);
alter table silencing_reagent_duplex_wells 
    ADD constraint silencing_reagent_duplex_wells_natural_key 
        unique(silencing_reagent_id,well_id);
/** Set the default nextval for Django ORM (AutoField) usage. **/
alter table silencing_reagent_duplex_wells 
    alter column id set default nextval('silencing_reagent_duplex_wells_sequence'::regclass);
alter sequence silencing_reagent_duplex_wells_sequence 
    owned by silencing_reagent_duplex_wells.id;


alter table silencing_reagent_duplex_wells 
    RENAME COLUMN silencing_reagent_id to silencingreagent_id;

/** 
  *** Copy the silencing reagent gene information from the link tables
  These tables were created to allow multiple genes to be associated with the
  siRNAs, for hairpin libraries.  This feature is not used in SS, and we'll 
  review genes at some point in the future. 
**/

update silencing_reagent si 
    set facility_gene_id = fvg.gene_id 
    from reagent_facility_genes fvg where fvg.reagent_id=si.reagent_id;
update silencing_reagent si 
    set vendor_gene_id = rvg.gene_id 
    from reagent_vendor_genes rvg where rvg.reagent_id=si.reagent_id;

drop table reagent_facility_genes cascade;
drop table reagent_vendor_genes cascade;

/**
  Make the well concentration fields units specific
**/

/** discontinue: 20150216
update well set molar_concentration = molar_concentration*1000000;
alter table well rename COLUMN molar_concentration to micro_molar_concentration;
**/

/** 
  Change timestamp fields to use timezone information:
  From the docs:
  The PostgreSQL backend stores datetimes as timestamp with time zone. 
  In practice, this means it converts datetimes from the connection’s time zone 
  to UTC on storage, and from UTC to the connection’s time zone on retrieval.
  
  As a consequence, if you’re using PostgreSQL, you can switch between 
  USE_TZ = False and USE_TZ = True freely. The database connection’s time zone 
  will be set to TIME_ZONE or UTC respectively, so that Django obtains correct 
  datetimes in all cases. You don’t need to perform any data conversions.
  
  Note: it appears that legacy Screensaver dates were stored without timezone
  information. This modification will apply the database default timezone 
  (America/New_York for orchestra), which will be correct.

** MOVED TO MIGRATION 0002 ***

alter table library alter column date_created SET DATA TYPE timestamp with time zone;
alter table library alter column date_loaded SET DATA TYPE timestamp with time zone; 
alter table library alter column date_publicly_available SET DATA TYPE timestamp with time zone;

alter table screen alter column date_created SET DATA TYPE timestamp with time zone;
alter table screen alter column date_loaded SET DATA TYPE timestamp with time zone; 
alter table screen alter column date_publicly_available SET DATA TYPE timestamp with time zone;

alter table screensaver_user alter column date_created SET DATA TYPE timestamp with time zone;
alter table screensaver_user alter column date_loaded SET DATA TYPE timestamp with time zone; 
alter table screensaver_user alter column date_publicly_available SET DATA TYPE timestamp with time zone;

alter table screen_result alter column date_loaded SET DATA TYPE timestamp with time zone;
alter table screen_result alter column date_publicly_available SET DATA TYPE timestamp with time zone;
alter table screen_result alter column date_created SET DATA TYPE timestamp with time zone;
 **/


/**
  Add id field to screensaver_user_role:
  Purpose: add an id to capture the natural ordering of the screensaver_user_role table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/

alter table screensaver_user_role add column id integer;
create sequence screensaver_user_role_sequence;
update screensaver_user_role set id=nextval('screensaver_user_role_sequence');


/**
  Create a many-to-many join table for the screen.collaborators field, then
  populate it using the legacy table (todo: remove the legacy table)
**/
CREATE TABLE "screen_collaborators" (
  "id" serial NOT NULL PRIMARY KEY, 
  "screen_id" integer NOT NULL, 
  "screensaveruser_id" integer NOT NULL, 
  UNIQUE ("screen_id", "screensaveruser_id")
  );
INSERT into screen_collaborators 
  ( select nextval('screen_collaborators_id_seq'), screen_id, collaborator_id 
      from collaborator_link );
        
ALTER TABLE "screen_collaborators" ADD CONSTRAINT "screen_collaborators_screen_id_fk" 
  FOREIGN KEY ("screen_id") 
  REFERENCES "screen" ("screen_id") DEFERRABLE INITIALLY DEFERRED;
ALTER TABLE "screen_collaborators" ADD CONSTRAINT "screen_collaborators_screensaver_user_id_fk" 
  FOREIGN KEY ("screensaveruser_id") 
  REFERENCES "screensaver_user" ("screensaver_user_id") DEFERRABLE INITIALLY DEFERRED;
CREATE INDEX "screen_collaborators_screen_idx" ON "screen_collaborators" ("screen_id");
CREATE INDEX "screen_collaborators_screensaver_user_idx" ON "screen_collaborators" ("screensaveruser_id");
DROP TABLE collaborator_link;

/** done - collaborator_link table **/

/**
  Create a foreignkey for the publication.screen field, then
  populate it using the legacy table (todo: remove the legacy table)
**/

ALTER TABLE "publication" ADD COLUMN "screen_id" integer NULL;
ALTER TABLE "publication" ADD CONSTRAINT "publication_screen_id_fk_screen_screen_id" 
  FOREIGN KEY ("screen_id") REFERENCES "screen" ("screen_id") DEFERRABLE INITIALLY DEFERRED;
UPDATE publication set screen_id = spl.screen_id
  FROM screen_publication_link spl
  where publication.publication_id=spl.publication_id;
CREATE INDEX "publication_screen_id_index" ON "publication" ("screen_id");
DROP TABLE screen_publication_link;

/** Must commit publication changes **/
COMMIT;
/** done - screen_publication_link table **/

BEGIN;
/** Reverse the publication->attached_file link, so that on_delete CASCADE works
    when deleting the publication 
**/
ALTER TABLE "attached_file" ADD COLUMN "publication_id" integer NULL UNIQUE;
CREATE INDEX "attached_file_publication" ON "attached_file" ("publication_id");
ALTER TABLE "attached_file" ADD CONSTRAINT "attached_file_publication_id_fk_publication_publication_id" 
  FOREIGN KEY ("publication_id") 
  REFERENCES "publication" ("publication_id") DEFERRABLE INITIALLY DEFERRED;
UPDATE attached_file set publication_id = p.publication_id
  FROM publication p
  where p.attached_file_id = attached_file.attached_file_id;
ALTER TABLE publication DROP COLUMN attached_file_id;

/** Must commit attached_file changes **/
COMMIT;
BEGIN;

/** convert attached_file file_contents into a "contents" field of type bytea **/

ALTER TABLE attached_file ADD COLUMN contents bytea;

UPDATE attached_file SET contents = loread(lo_open(file_contents, 262144), 10000000) 
  WHERE attached_file.file_contents IS NOT NULL;
COMMIT;
BEGIN;
  
ALTER TABLE attached_file ALTER COLUMN file_contents DROP NOT NULL;

/**
  TODO:
  use "\lo_list"
  and
  "\lo_unlink" to drop largeobjects

 *** TODO: cannot do following:
 *** pg_largeobject is owned by user "postgres" on orchestra and cannot 
 *** be deleted.
DELETE FROM pg_largeobject USING attached_file WHERE loid=file_contents;
ALTER TABLE attached_file DROP COLUMN file_contents;
**/

/**
  Create a foreignkey for the publication.reagent field, then
  populate it using the legacy table.
  TODO: this appears to be unused?
**/

ALTER TABLE "publication" ADD COLUMN "reagent_id" integer NULL;
ALTER TABLE "publication" ADD CONSTRAINT "publication_reagent_id_fk_reagent_reagent_id" 
  FOREIGN KEY ("reagent_id") REFERENCES "reagent" ("reagent_id") DEFERRABLE INITIALLY DEFERRED;
UPDATE publication set reagent_id = rpl.reagent_id
  FROM reagent_publication_link rpl
  where publication.publication_id=rpl.publication_id;
CREATE INDEX "publication_reagent_id_index" ON "publication" ("reagent_id");
DROP TABLE reagent_publication_link;

/** done - reagent_publication_link table **/


/**
  Create a many-to-many join table for the datacolumn.derived_from field, then
  populate it using the legacy table (todo: remove the legacy table)
**/
CREATE TABLE "data_column_derived_from_columns" (
  "id" serial NOT NULL PRIMARY KEY, 
  "from_datacolumn_id" integer NOT NULL, 
  "to_datacolumn_id" integer NOT NULL, 
  UNIQUE ("from_datacolumn_id", "to_datacolumn_id")
);
INSERT into data_column_derived_from_columns 
  ( select nextval('data_column_derived_from_columns_id_seq'), 
      derived_from_data_column_id, derived_data_column_id 
      from data_column_derived_from_link );

ALTER TABLE "data_column_derived_from_columns" 
  ADD CONSTRAINT "fk_from_datacolumn_id_data_column_id" 
  FOREIGN KEY ("from_datacolumn_id") 
  REFERENCES "data_column" ("data_column_id") DEFERRABLE INITIALLY DEFERRED;
ALTER TABLE "data_column_derived_from_columns" 
  ADD CONSTRAINT "fk_to_datacolumn_id_data_column_id" 
  FOREIGN KEY ("to_datacolumn_id") 
  REFERENCES "data_column" ("data_column_id") DEFERRABLE INITIALLY DEFERRED;
CREATE INDEX "data_column_derived_from_columns_from_dc_idx" 
  ON "data_column_derived_from_columns" ("from_datacolumn_id");
CREATE INDEX "data_column_derived_from_columns_to_dc_idx" 
  ON "data_column_derived_from_columns" ("to_datacolumn_id");

/**
20160408
TODO: verify that the ldld here is the same as the sr.date_loaded
# 'last_data_loading_date': literal_column(
#     '( select activity.date_created '
#     '  from activity '
#     '  join administrative_activity aa using(activity_id) '
#     '  join screen_update_activity on update_activity_id=activity_id  '
#     "  where administrative_activity_type = 'Screen Result Data Loading' " 
#     '  and screen_id=screen.screen_id '
#     '  order by date_created desc limit 1 )'
update screen_result set date_loaded = TODO: 
**/

ALTER TABLE assay_well ALTER COLUMN is_positive SET default false;

/** 
20170413 Delete "data_loading" assay plates:
- these plates were created exclusively to track data loading replicate counts,
- per JAS we will no longer track these values
- accounting algorithm in SS1 is inaccurate
**/
delete from assay_plate 
  where screen_result_data_loading_id is not null 
  and library_screening_id is null;

DROP TABLE cell_lineage;
DROP TABLE cell_related_projects;
DROP TABLE cell_markers;
DROP TABLE cell_growth_properties;
DROP TABLE experimental_cell_information;
DROP TABLE cell_update_activity;
DROP TABLE primary_cell;
DROP TABLE cell;
DROP TABLE gene_old_entrezgene_id;
DROP TABLE gene_old_entrezgene_symbol;
DROP TABLE silencing_reagent_non_targetted_genbank_accession_number;


COMMIT;
