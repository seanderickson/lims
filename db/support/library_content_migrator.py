# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import json
import logging

from django.core.exceptions import ObjectDoesNotExist
from django.db import models
from django.db.models import Count
from django.forms.models import model_to_dict
from django.utils import timezone

from db.api import SmallMoleculeReagentResource, WellResource, \
    SilencingReagentResource, NaturalProductReagentResource, LibraryResource
from reports.api import compare_dicts, is_empty_diff
from reports.models import ApiLog
from db.migrations import create_log_time


logger = logging.getLogger(__name__)

class Migrator:
    '''
    purpose: support class for the library content version migrations
    '''

    smrResource = SmallMoleculeReagentResource()
    silencingReagentResource = SilencingReagentResource()
    naturalProductResource = NaturalProductReagentResource()
    wellResource = WellResource()
    libraryResource = LibraryResource()

    rnai_keys = [
                'well_id', 'vendor_identifier', 'vendor_name', 'vendor_batch_id',
                'vendor_name_synonym', #'substance_id',
                'sequence', 'silencing_reagent_type',

                'vendor_entrezgene_id',
                'vendor_entrezgene_symbols',
                'vendor_gene_name',
                'vendor_gene_species',
                'vendor_genbank_accession_numbers',

                'facility_entrezgene_id',
                'facility_entrezgene_symbols',
                'facility_gene_name',
                'facility_gene_species',
                'facility_genbank_accession_numbers',

                'duplex_wells'
                ]
            
    rnai_sql = '''select 
well_id, vendor_identifier, vendor_name, vendor_batch_id, vendor_name_synonym,
sequence, silencing_reagent_type,

vg.entrezgene_id as vendor_entrezgene_id,
(select '["' || array_to_string(array_agg(entrezgene_symbol), '","') || '"]'
    from (select entrezgene_symbol from gene_symbol gs 
          where gs.gene_id = sr.vendor_gene_id 
          order by ordinal) a ) as vendor_entrezgene_symbols,
vg.gene_name as vendor_gene_name,
vg.species_name as vendor_gene_species,
(select '["' || array_to_string(array_agg(genbank_accession_number), '","') || '"]'
    from (select genbank_accession_number from gene_genbank_accession_number ga
    where ga.gene_id = sr.vendor_gene_id order by ga.id ) a ) as vendor_genbank_accession_numbers,

fg.entrezgene_id as facility_entrezgene_id,
(select '["' || array_to_string(array_agg(entrezgene_symbol), '","') || '"]'
    from (select entrezgene_symbol from gene_symbol gs 
          where gs.gene_id = sr.facility_gene_id 
          order by ordinal) a ) as facility_entrezgene_symbols,
fg.gene_name as facility_gene_name,
fg.species_name as facility_gene_species,
(select '["' || array_to_string(array_agg(genbank_accession_number), '","') || '"]'
    from (select genbank_accession_number from gene_genbank_accession_number ga
    where ga.gene_id = sr.facility_gene_id order by ga.id) a ) as facility_genbank_accession_numbers,

(select '["' || array_to_string(array_agg(well_id), '","') || '"]' 
    from (select well_id from silencing_reagent_duplex_wells dw 
    where dw.silencingreagent_id=r.reagent_id order by dw.well_id) a ) as duplex_wells

from reagent r join silencing_reagent sr using(reagent_id)
left join gene vg on(vendor_gene_id=vg.gene_id)
left join gene fg on(facility_gene_id=fg.gene_id)
where r.library_contents_version_id=%s order by well_id;
'''
    sm_keys = [
                'well_id', 'vendor_identifier', 'vendor_name', 'vendor_batch_id',
                'vendor_name_synonym',# 'substance_id',
                'inchi', 'smiles', 
                'molecular_formula', 'molecular_mass', 'molecular_weight',
                'compound_name', 'pubchem_cid', 'chembl_id', 'chembank_id',
                ]
            
    sm_sql = '''select 
well_id, vendor_identifier, vendor_name, vendor_batch_id, vendor_name_synonym,
inchi, smiles, 
molecular_formula, molecular_mass, molecular_weight,
(select '["' || array_to_string(array_agg(compound_name), '","') || '"]' 
    from (select compound_name from small_molecule_compound_name smr 
    where smr.reagent_id=r.reagent_id order by ordinal) a) as compound_name,
(select '["' || array_to_string(array_agg(pubchem_cid), '","') || '"]' 
    from ( select pubchem_cid from small_molecule_pubchem_cid p 
           where p.reagent_id=r.reagent_id order by id ) a ) as pubchem_cid,
(select '["' || array_to_string(array_agg(chembl_id), '","') || '"]' 
    from (select chembl_id from small_molecule_chembl_id cb 
          where cb.reagent_id=r.reagent_id order by id ) a )   as chembl_id,
(select '["' || array_to_string(array_agg(chembank_id), '","') || '"]' 
    from (select chembank_id from small_molecule_chembank_id cbk 
          where cbk.reagent_id=r.reagent_id order by id ) a ) as chembank_id 
from reagent r join small_molecule_reagent using(reagent_id)
where r.library_contents_version_id=%s order by well_id;
'''
    
    np_keys = [ 
                'well_id', 'vendor_identifier', 'vendor_name', 'vendor_batch_id',
                'vendor_name_synonym',# 'substance_id',
                ]
    np_sql = '''select 
well_id, vendor_identifier, vendor_name, vendor_batch_id, vendor_name_synonym
from reagent r join natural_product_reagent using(reagent_id)
where r.library_contents_version_id=%s order by well_id;
'''
          
    def do_migration(self, apps, schema_editor, screen_type=None):                
        i=0
        
        query = apps.get_model('db','LibraryContentsVersion').objects.all()
        if screen_type:
            query = (query.filter(library__screen_type=screen_type))
            #.exclude(library__library_type='natural_products'))
        library_ids = [x['library'] for x in (query
                .values('library')  # actually, library id's here
                .order_by('library') )]
        logger.info('libraries to consider: %r', library_ids)
    
        for library in (apps.get_model('db','Library').objects.all()
                        .filter(library_id__in=library_ids)):
            
            prev_version = None
            for version in (library.librarycontentsversion_set.all()
                            .order_by('version_number')):
                # create an apilog for the library
                activity = (version.library_contents_loading_activity.activity)
                log = ApiLog()
                if getattr(activity.performed_by, 'ecommons_id', None):
                    log.username = activity.performed_by.ecommons_id
                    
                    if log.username == 'dwrobel':
                        log.username = 'djw11'
                        log.user_id = 761
                    
                if getattr(activity.performed_by, 'user', None):
                    log.user_id = getattr(
                        activity.performed_by.user, 'id', log.username)
                if not log.user_id:
                    log.user_id = 1
                log.ref_resource_name = self.libraryResource._meta.resource_name
                log.api_action = 'PATCH'
                log.json_field = {
                    'migration': 'Library (contents)',
                    'data': { 
                        'library_contents_version.activity_id': 
                            activity.activity_id,
                    }
                }
                log.key = library.short_name
                log.uri = '/'.join(['library',log.key])
                log.date_time = create_log_time(log.key, activity.date_of_activity)
                log.diffs = {
                    'version_number': [
                        prev_version.version_number if prev_version else 0, 
                        version.version_number]
                        }
                log.comment = activity.comments
                log.save()
                if prev_version:
                    self.diff_library_wells(
                        schema_editor,library, prev_version, version, log)            

                prev_version = version
                
                # add version to library
                library.version_number = version.version_number
                library.loaded_by = activity.performed_by
            
            library.save()
            i=i+1

            ## TODO: 20140826
            ## - set all the reagent.library values
            ## - prune out all the old reagents

    
        print 'processed: ', i, 'libraries'

    def diff_library_wells(
        self, schema_editor,library, version1, version2, parent_log):
        
        connection = schema_editor.connection
        
        screen_type = library.screen_type;
        logger.info(
            'processing library: %r, screen type: %r, versions: %r, exp wells: %d ',
            library.short_name, library.screen_type, [version1, version2],
            library.experimental_well_count)
        i = 0
        keys = self.sm_keys
        sql = self.sm_sql
        if library.library_type == 'natural_products':
            keys = self.np_keys
            sql = self.np_sql
        if screen_type == 'rnai':
            keys = self.rnai_keys
            sql = self.rnai_sql
            
        prev_well_map = None
        for version in [version1, version2]:
            logger.debug('version: %s, execute sql %r', 
                version.library_contents_version_id, sql)
            cursor = connection.cursor()
            cursor.execute(sql, [version.library_contents_version_id])
            _list = cursor.fetchall()
            logger.debug('executed, list: %r', _list)
            
            if len(_list) == 0:
                logger.error('no wells for %r', library.short_name)
                continue

            cumulative_diff_keys = set()
            
            well_map = dict((x[0], x) for x in _list)
            if prev_well_map:
                for key,row in prev_well_map.items():
                    new_row = well_map.get(key, None)
                    if new_row:
                        prev_well = dict(zip(keys,row))
                        new_well = dict(zip(keys,new_row))
                        log = self.create_well_log(
                            version, prev_well, new_well, parent_log)
                        if log:
                            if log.diffs:
                                cumulative_diff_keys.update(
                                    set(log.diffs.keys()))
                            i += 1
                    else:
                        logger.error(
                            'no new well/reagent entry found for %r', key)
                    if i>0 and i % 1000 == 0:
                        logger.info(
                            '===created logs for %s: %r-%r, %d', 
                            library.short_name, version1, version2, i)
            else:
                prev_well_map = well_map
            
            parent_log.save()
        logger.info('===created %d logs for %r, %r', 
            i, library.short_name, parent_log.comment)
    
    def create_well_log(self, version, prev_dict, current_dict, parent_log):
        
        difflog = compare_dicts(
            prev_dict, current_dict,
            excludes=['reagent_id', 'resource_uri'])
        if is_empty_diff(difflog):
            return None
        
        log = ApiLog()
        
        if parent_log:
            log.username = parent_log.username
            log.user_id = parent_log.user_id
            log.comment = parent_log.comment
            log.date_time = parent_log.date_time
        else:
            activity = version.library_contents_loading_activity.activity
            if getattr(activity.performed_by, 'ecommons_id', None):
                log.username = activity.performed_by.ecommons_id
            else:
                log.username = 'sde'
                
            if getattr(activity.performed_by, 'login_id', None):
                log.username = activity.performed_by.login_id
            # FIXME
            log.user_id = 1
        
            log.date_time = activity.date_of_activity

        log.ref_resource_name = self.wellResource._meta.resource_name
        log.api_action = 'PATCH'
        log.key = prev_dict['well_id']
        log.uri = '/db/api/v1/well/'+log.key 
        log.diffs = difflog
        log.json_field = json.dumps({
            'version': version.version_number })
        log.parent_log = parent_log
        log.save()
        return log
    
