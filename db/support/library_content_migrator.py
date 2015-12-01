# -*- coding: utf-8 -*-

import json
import logging

from django.core.exceptions import ObjectDoesNotExist
from django.db import models
from django.db.models import Count
from django.utils import timezone, tzinfo
from django.utils.timezone import make_aware
from tastypie.bundle import Bundle
from tastypie.resources import Resource

from db.api import SmallMoleculeReagentResource, WellResource, \
    SilencingReagentResource, NaturalProductReagentResource, LibraryResource
from db.models import SmallMoleculeReagent
from reports.api import compare_dicts, is_empty_diff
from reports.models import ApiLog


# from south.utils import datetime_utils as datetime
# from south.db import db
# from south.v2 import DataMigration
logger = logging.getLogger(__name__)

class Migrator:
    '''
    purpose: support class for the 0013,0014 library content version migrations
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
#                           .exclude(library__library_type='natural_products'))
        library_ids = [x['library'] for x in (query
                .values('library')  # actually, library id's here
#                 .annotate(count=Count('library'))
#                 .filter(count__gt=1)
                .order_by('library') )]
        logger.info(str(('libraries to consider', library_ids)))
    
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
                if getattr(activity.performed_by, 'user', None):
                    log.user_id = getattr(activity.performed_by.user, 'id', log.username)
                if not log.user_id:
#                     logger.info(str(("can't find a user id for version", version, activity)))
                    log.user_id = 1
                log.date_time = make_aware(
                    activity.date_created,timezone.get_default_timezone())
                log.ref_resource_name = self.libraryResource._meta.resource_name
                # TODO: what action? could also be a REST specifier, i.e. 'PATCH'
                log.api_action = 'PUT'
                # store the old type in the catch-all field
                log.json_field = json.dumps( {
                    'administrative_activity_type': 
                    version.library_contents_loading_activity.administrative_activity_type
                    })
                log.uri = self.libraryResource.get_resource_uri(bundle_or_obj=library)
                log.key = '/'.join([str(x) for x in (
                    self.libraryResource.detail_uri_kwargs(library).values()) ])
                log.diff_keys = json.dumps(['version_number'])
                log.diffs = json.dumps([
                    prev_version.version_number if prev_version else 0, version.version_number])
                log.comment = activity.comments
                log.save()
                logger.debug(str(('created log', log)))
                if prev_version:
                    self.diff_library_wells(schema_editor,library, prev_version, version, log)            

                
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

    def diff_library_wells(self, schema_editor,library, version1, version2, parent_log):
        
        connection = schema_editor.connection
        
        screen_type = library.screen_type;
        logger.info(str(('processing: ', library.short_name, 
           'type', library.screen_type, library.library_type, 
           'versions',  [version1, version2],
           'experimental wells', library.experimental_well_count)) )
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
                logger.error(str(('no wells for ', library.short_name)))
                continue

            cumulative_diff_keys = set()
            
            well_map = dict((x[0], x) for x in _list)
            if prev_well_map:
                for key,row in prev_well_map.items():
                    new_row = well_map.get(key, None)
                    if new_row:
                        prev_well = dict(zip(keys,row))
                        new_well = dict(zip(keys,new_row))
                        logger.debug(str(('===create log for ', prev_well['well_id'])))
                        log = self.create_well_log(version, prev_well, new_well, parent_log)
                        if log:
                            if log.diff_keys:
                                cumulative_diff_keys.update(set(json.loads(log.diff_keys)))
                            i += 1
                    else:
                        logger.error(str(('no new well/reagent entry found for', key)))
                    if i>0 and i % 1000 == 0:
                        logger.info(str(('===created logs for ', library.short_name, i)))
            else:
                prev_well_map = well_map
            
            comment = parent_log.comment
            if not comment or len(comment)==0:
                comment = ''
            else:
                comment += ': '
            parent_log.comment = comment + json.dumps(list(cumulative_diff_keys))
            parent_log.save()
        logger.info(str(('===created logs for ', library.short_name, parent_log.comment, 'count', i)))
    
    def create_well_log(self, version, prev_dict, current_dict, parent_log):
        
        difflog = compare_dicts(
            prev_dict, current_dict,
            excludes=['reagent_id', 'resource_uri'])
        if is_empty_diff(difflog):
            return None
        
        activity = version.library_contents_loading_activity.activity
        log = ApiLog()
        #             log.username = activity.performed_by.ecommons_id
        #             log.user_id = activity.performed_by.user.id 
        
        if getattr(activity.performed_by, 'ecommons_id', None):
            log.username = activity.performed_by.ecommons_id
        else:
            log.username = 'sde'
            
        if getattr(activity.performed_by, 'login_id', None):
            log.username = activity.performed_by.login_id
        # FIXME
        log.user_id = 1
        
        log.date_time = make_aware(activity.date_created,timezone.get_default_timezone())
        log.ref_resource_name = self.wellResource._meta.resource_name
        # TODO: what types here? could also be a REST specifier, i.e. 'PATCH'
        log.api_action = 'MIGRATION'
        log.key = prev_dict['well_id']
        log.uri = '/db/api/v1/well/'+log.key # self.wellResource.get_resource_uri(bundle_or_obj=well)
        log.diff_dict_to_api_log(difflog)
        
        log.json_field = json.dumps({
            'version': version.version_number })
        log.parent_log = parent_log
        log.save()
        return log
    
#     def diff_library_wells_old(self, library, screen_type):
#         
#         versions = [x for x in (library.librarycontentsversion_set.all()
#                 .order_by('version_number')) ] 
#         if len(versions) < 2:
#             return
# 
#         logger.info(str(('processing: ', library.short_name, 
#            'type', library.screen_type, library.library_type, 
#            'versions',  [x.version_number for x in versions],
#            'experimental wells', library.experimental_well_count)) )
#         
#         if screen_type == 'rnai':
#             keys = [
#                 'well_id', 'vendor_identifier', 'vendor_name', 'vendor_batch_id',
#                 'vendor_name_synonym','substance_id',
#                 'sequence', 'silencing_reagent_type',
#                 'vendor_entrezgene_id',
#                 'vendor_entrezgene_symbols',
#                 'vendor_gene_name',
#                 'vendor_gene_species',
#                 'vendor_genbank_accession_numbers',
#                 'facility_entrezgene_id',
#                 'facility_entrezgene_symbols',
#                 'facility_gene_name',
#                 'facility_gene_species',
#                 'facility_genbank_accession_numbers',
#                 'silencing_reagent_duplex_wells'
#                 ]
#             
#             sql = '''select 
# well_id, vendor_identifier, vendor_name, vendor_batch_id, vendor_name_synonym,substance_id,
# vg.entrezgene_id as vendor_entrezgene_id,
# vg.gene_name as vendor_gene_name,
# vg.species_name as vendor_gene_species,
# (select '["' || array_to_string(array_agg(entrezgene_symbol), '","') || '"]'
#     from (select entrezgene_symbol from gene_symbol gs 
#           where gs.gene_id = sr.vendor_gene_id 
#           order by ordinal) a ) as vendor_entrezgene_symbols,
# (select '["' || array_to_string(array_agg(genbank_accession_number), '","') || '"]'
#     from gene_genbank_accession_number ga
#     where ga.gene_id = sr.vendor_gene_id ) as vendor_genbank_accession_numbers,
# fg.entrezgene_id as facility_entrezgene_id,
# fg.gene_name as facility_gene_name,
# fg.species_name as facility_gene_species,
# (select '["' || array_to_string(array_agg(entrezgene_symbol), '","') || '"]'
#     from (select entrezgene_symbol from gene_symbol gs 
#           where gs.gene_id = sr.facility_gene_id 
#           order by ordinal) a ) as facility_entrezgene_symbols,
# (select '["' || array_to_string(array_agg(genbank_accession_number), '","') || '"]'
#     from gene_genbank_accession_number ga
#     where ga.gene_id = sr.facility_gene_id ) as facility_genbank_accession_numbers,
# (select '["' || array_to_string(array_agg(well_id), '","') || '"]' 
#     from silencing_reagent_duplex_wells dw where dw.silencingreagent_id=r.reagent_id) as duplex_wells
# from reagent r join silencing_reagent sr using(reagent_id)
# left join gene vg on(vendor_gene_id=vg.gene_id)
# left join gene fg on(facility_gene_id=fg.gene_id)
# where r.library_contents_version_id=%s order by well_id;
# '''
#         else:
#             # create hash: version->well->smr
#             keys = [
#                 'well_id', 'vendor_identifier', 'vendor_name', 'vendor_batch_id',
#                 'vendor_name_synonym','substance_id',
#                 'inchi', 'smiles', 
#                 'molecular_formula', 'molecular_mass', 'molecular_weight',
#                 'compound_name', 'pubchem_cid', 'chembl_id', 'chembank_id',
#                 ]
#             
#             sql = '''select 
# well_id, vendor_identifier, vendor_name, vendor_batch_id, vendor_name_synonym,substance_id,
# inchi, smiles, 
# molecular_formula, molecular_mass, molecular_weight,
# (select '["' || array_to_string(array_agg(compound_name), '","') || '"]' 
#     from (select compound_name from small_molecule_compound_name smr 
#     where smr.reagent_id=r.reagent_id order by ordinal) a) as compound_name,
# (select '["' || array_to_string(array_agg(pubchem_cid), '","') || '"]' 
#     from small_molecule_pubchem_cid p where p.reagent_id=r.reagent_id) as pubchem_cid,
# (select '["' || array_to_string(array_agg(chembl_id), '","') || '"]' 
#     from small_molecule_chembl_id cb where cb.reagent_id=r.reagent_id)   as chembl_id,
# (select '["' || array_to_string(array_agg(chembank_id), '","') || '"]' 
#     from small_molecule_chembank_id cbk where cbk.reagent_id=r.reagent_id)  as chembank_id 
# from reagent r join small_molecule_reagent using(reagent_id)
# where r.library_contents_version_id=%s order by well_id;
# '''
#         prev_well_map = None
#         i = 0
#         for version in versions:
#             _list = db.execute(sql, [version.library_contents_version_id])
#             
#             if len(_list) == 0:
#                 logger.error(str(('no wells for ', library.short_name, version.version)))
#                 continue
# 
#             well_map = dict((x[0], x) for x in _list)
#             if prev_well_map:
#                 for key,row in prev_well_map.items():
#                     new_row = well_map.get(key, None)
#                     if new_row:
#                         prev_well = dict(zip(keys,row))
#                         new_well = dict(zip(keys,new_row))
#                         log = self.create_well_log(version, prev_well, new_well)
#                         if log:
#                             i += 1
#                     else:
#                         logger.error(str(('no new well/reagent entry found for', key)))
#             else:
#                 prev_well_map = well_map
#                 
#         logger.info(str(('===created logs for ', library.short_name, i)))
    
    
    ########### OLD ORM METHODOLOGY: (not performant)        
#     def diff_smr(self,r1,r2):    
#         
#         bundle1 = self.smrResource.full_dehydrate(Bundle(obj=r1))
#         bundle2 = self.smrResource.full_dehydrate(Bundle(obj=r2))
#         
#         diff_log = compare_dicts(
#             bundle1.data, bundle2.data,
#             excludes=['reagent_id', 'resource_uri'])
#         return diff_log
# 
#     def diff_rnai(self,r1,r2):  
#         bundle1 = self.silencingReagentResource.full_dehydrate(Bundle(obj=r1))
#         bundle2 = self.silencingReagentResource.full_dehydrate(Bundle(obj=r2))
#         
#         diff_log = compare_dicts(
#             bundle1.data, bundle2.data,
#             excludes=['reagent_id', 'resource_uri'])
#         return diff_log
# 
#     def diff_natural_product(self,r1,r2):    
#         bundle1 = self.naturalProductResource.full_dehydrate(Bundle(obj=r1))
#         bundle2 = self.naturalProductResource.full_dehydrate(Bundle(obj=r2))
#         
#         diff_log = compare_dicts(
#             bundle1.data, bundle2.data,
#             excludes=['reagent_id', 'resource_uri'])
#         return diff_log
# 
#     def diff_library_wells_old_orm_method(self, library):
#         
#         versions = [x.version_number 
#             for x in (library.librarycontentsversion_set.all()
#                 .order_by('version_number')) ] 
#         if len(versions) < 2:
#             return
# 
#         logger.info(str(('processing: ', library.short_name, 
#                'type', library.screen_type, library.library_type, 
#                'versions',  versions,
#                'experimental wells', library.experimental_well_count)) )
#         #build a hash of well->[reagents by version]
#             
#         base_query = None
#         diff_function = None      
#         if library.screen_type == 'rnai':
#             base_query = orm.SilencingReagent.objects.all()
#             diff_function = self.diff_rnai
#             attribute = 'silencingreagent'
#         else:
#             if library.library_type == 'natural_products':
#                 base_query = orm.NaturalProductReagent.objects.all()
#                 diff_function = self.diff_natural_product
#                 attribute = 'naturalproductreagent'
#             else:
#                 base_query = orm.SmallMoleculeReagent.objects.all()
#                 diff_function = self.diff_smr
#                 attribute = 'smallmoleculereagent'
#         logs_created = 0
#        
#         for well in library.well_set.all():
#             
#             prev_version_reagent = None
#         
#             for reagent in (
#                 base_query
#                 .filter(reagent__well=well)
#                 .order_by('reagent__library_contents_version__version_number')):
#         
#                 if prev_version_reagent:
#                     if self.create_diff_log(well, 
#                                             prev_version_reagent, reagent,
#                                             diff_function):
#                         logs_created +=1
#                         
#                 prev_version_reagent = reagent
#                 
#             if logs_created and logs_created % 1000 == 0:
#                 logger.info(str(('created ',logs_created, ' logs for ',library.short_name )))                        
#         logger.info(str(('library', library.short_name, 
#                          'logs_created', logs_created)) )
# 
# 
# 
#     def create_diff_log_old_orm(self,well,a,b,diff_function):
#         '''
#         @param a reagent a
#         @param b reagent b
#         '''
#         difflog = diff_function(a,b)
#         if ( 'added_keys' in difflog or 
#              'removed_keys' in difflog or
#              'diff_keys' in difflog ):
#             diff_keys = difflog['diff_keys'] or []
#             diffs = difflog['diffs'] or {}
#             diff_keys.append('substance_id')
#             diffs['substance_id'] = [a.reagent.substance_id, 
#                                      b.reagent.substance_id]
#             difflog['diff_keys'] = diff_keys
#             difflog['diffs'] = diffs
#             
#             activity = (b.reagent.library_contents_version
#                          .library_contents_loading_activity.activity)
#             log = ApiLog()
# 
#             if getattr(activity.performed_by, 'ecommons_id', None):
#                 log.username = activity.performed_by.ecommons_id
#             if getattr(activity.performed_by, 'login_id', None):
#                 log.username = activity.performed_by.login_id
#             if not log.username:
#                 logger.warn(str(("can't find alog.username, library, ss user id ", 
#                     b.reagent.library_contents_version.library.short_name, activity.performed_by.id)))
#                 log.username = 'sde'
# 
#             # FIXME: to make this work, we need to copy all of the screensaver_user 
#             # table to our user table
#             log.user_id = 1
#             ### end FIXME ###    
#             
#             log.date_time = make_aware(activity.date_created,timezone.get_default_timezone())
#             log.ref_resource_name = self.wellResource._meta.resource_name
#             # TODO: what types here? could also be a REST specifier, i.e. 'PATCH'
#             log.api_action = 'MIGRATION'
#             log.uri = self.wellResource.get_resource_uri(bundle_or_obj=well)
#             log.key = '/'.join([str(x) for x in (
#                 self.wellResource.detail_uri_kwargs(well).values()) ])
#             log.diff_dict_to_api_log(difflog)
#             
#             log.json_field = json.dumps({
#                 'version': b.reagent.library_contents_version.version_number })
#             log.save()
#             return log
#         return None
#     
