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

update well set molar_concentration = molar_concentration*1000000;
alter table well rename COLUMN molar_concentration to micro_molar_concentration;


COMMIT;
