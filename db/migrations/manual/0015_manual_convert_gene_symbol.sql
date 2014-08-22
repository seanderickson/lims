BEGIN;

INSERT INTO schema_history (screensaver_revision, date_updated, comment)
SELECT
20140013,
current_timestamp,
'migrate gene_symbol';

/**
  Purpose: add an id to capture the natural ordering of the gene_symbol table; 
  because the table has no ID, Django ORM doesn't know what to do with it -
  and the South migration cannot do it either.
**/

alter table gene_symbol add column id integer;
create sequence gene_symbol_sequence;
update gene_symbol set id=nextval('gene_symbol_sequence');
alter table gene_symbol drop CONSTRAINT gene_symbol_pkey;
alter table gene_symbol add primary key(id);
alter table gene_symbol ADD constraint gene_symbol_natural_key unique(gene_id,ordinal);
# Set the default nextval for Django ORM (AutoField) usage.
alter table gene_symbol alter column id set default nextval('gene_symbol_sequence'::regclass);
# TODO: this fails:
#alter sequence gene_symbol_sequence owned by gene_symbol.id;
COMMIT;
