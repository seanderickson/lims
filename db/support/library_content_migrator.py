# -*- coding: utf-8 -*-

import json

from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import DataMigration
from django.db import models

from db.api import SmallMoleculeReagentResource, WellResource,\
    SilencingReagentResource, NaturalProductReagentResource, LibraryResource
from tastypie.resources import Resource
from tastypie.bundle import Bundle
from reports.models import ApiLog
from django.utils import timezone, tzinfo

import logging
from django.utils.timezone import make_aware, UTC

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
            
    def diff_smr(self,r1,r2):    
        bundle1 = self.smrResource.full_dehydrate(Bundle(obj=r1))
        bundle2 = self.smrResource.full_dehydrate(Bundle(obj=r2))
        
        diff_log = self.smrResource.compare_dicts(
            bundle1.data, bundle2.data,
            excludes=['reagent_id', 'resource_uri'])
        return diff_log

    def diff_rnai(self,r1,r2):  
#         print 'considering', r1.reagent.substance_id  
        bundle1 = self.silencingReagentResource.full_dehydrate(Bundle(obj=r1))
        bundle2 = self.silencingReagentResource.full_dehydrate(Bundle(obj=r2))
        
        diff_log = self.silencingReagentResource.compare_dicts(
            bundle1.data, bundle2.data,
            excludes=['reagent_id', 'resource_uri'])
        return diff_log

    def diff_natural_product(self,r1,r2):    
        bundle1 = self.naturalProductResource.full_dehydrate(Bundle(obj=r1))
        bundle2 = self.naturalProductResource.full_dehydrate(Bundle(obj=r2))
        
        diff_log = self.naturalProductResource.compare_dicts(
            bundle1.data, bundle2.data,
            excludes=['reagent_id', 'resource_uri'])
        return diff_log

    def do_migration(self, orm, screen_type=None):                
        i=0
        
        # First, get only the libraries having > 1 version
        from django.db.models import Count
        query = orm.LibraryContentsVersion.objects.all()
        if screen_type:
            query = query.filter(library__screen_type=screen_type)
        library_ids = [x['library'] for x in ( 
            query
                .values('library')  # actually, library id's here
                .annotate(count=Count('library'))
                .filter(count__gt=1)
                .order_by('library') )]
        logger.info(str(('libraries to consider', library_ids)))
    
        for library in (orm.Library.objects.all()
                        .filter(library_id__in=library_ids)):
#                         .filter(screen_type='small_molecule')
#                         .exclude(library_type='natural_products')):
            versions = [x.version_number for x in 
                            (library.librarycontentsversion_set.all()
                             .order_by('version_number')) ] 
            if len(versions) < 2:
                continue

            logger.info(str(('processing: ', library.short_name, 
                   'type', library.screen_type, library.library_type, 
                   'versions',  versions,
                   'experimental wells', library.experimental_well_count)) )
            
            #build a hash of well->[reagents by version]
            
            base_query = None
            diff_function = None      
            if library.screen_type == 'rnai':
                base_query = orm.SilencingReagent.objects.all()
                diff_function = self.diff_rnai
            else:
                if library.library_type == 'natural_products':
                    base_query = orm.NaturalProductReagent.objects.all()
                    diff_function = self.diff_natural_product
                else:
                    base_query = orm.SmallMoleculeReagent.objects.all()
                    diff_function = self.diff_smr
            
            logs_created = 0
            for well in library.well_set.all():
                prev_version_reagent = None
                for reagent in (
                    base_query
                    .filter(reagent__well=well)
                    .order_by('reagent__library_contents_version__version_number')):

                    if prev_version_reagent:
                        if self.create_diff_log(well, 
                                                prev_version_reagent, reagent,
                                                diff_function):
                            logs_created +=1
                    prev_version_reagent = reagent
            logger.info(str(('library', library.short_name, 
                             'logs_created', logs_created)) )
            
            prev_version = None
            for version in (library.librarycontentsversion_set.all()
                            .order_by('version_number')):
                # create an apilog for the library
                activity = (version.library_contents_loading_activity.activity)
                log = ApiLog()
                log.username = activity.performed_by.ecommons_id
                log.user_id = activity.performed_by.user.id 
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
                log.diffs = json.dumps([prev_version, version.version_number])
                log.comment = activity.comments
                log.save()
                logger.info(str(('created log', log)))
                prev_version = version.version_number
                
                # add version to library
                library.version_number = version.version_number
                library.loaded_by = activity.performed_by
            library.save()
            i=i+1
            
#             wells = {}
#             prev_version = None
#             for version in (library.librarycontentsversion_set.all()
#                             .order_by('version_number')):
#                 
#                 for reagent in (base_query
#                             .filter(reagent__well__library=library)
#                             .filter(reagent__library_contents_version=version)
#                             .order_by('reagent__well__well_id')):
#                     well = reagent.reagent.well
#                     if well in wells:
#                         wells[well].append(reagent)
#                     else:
#                         wells[well] = [reagent]
# 
#                 # create an apilog for the library
#                 activity = (version.library_contents_loading_activity.activity)
#                 log = ApiLog()
#                 log.username = activity.performed_by.ecommons_id
#                 log.user_id = activity.performed_by.user.id 
#                 log.date_time = activity.date_created
#                 log.ref_resource_name = self.libraryResource._meta.resource_name
#                 # TODO: what types here? could also be a REST specifier, i.e. 'PATCH'
#                 log.api_action = 'PUT'
#                 log.json_field = json.dumps( {
#                     'administrative_activity_type': 
#                     version.library_contents_loading_activity.administrative_activity_type
#                     })
#                 log.uri = self.libraryResource.get_resource_uri(bundle_or_obj=library)
#                 log.key = '/'.join([str(x) for x in (
#                     self.libraryResource.detail_uri_kwargs(library).values()) ])
#                 log.diff_keys = json.dumps(['version'])
#                 log.diffs = json.dumps([prev_version, version.version_number])
#                 log.comment = activity.comments
#                 log.save()
#                 print 'log', log
#                 prev_version = version.version_number
#             
#             self.diff_versions(wells,diff_function)
    
        print 'processed: ', i, 'libraries'

    def create_diff_log(self,well,a,b,diff_function):
        difflog = diff_function(a,b)
        if ( 'added_keys' in difflog or 
             'removed_keys' in difflog or
             'diff_keys' in difflog ):
            diff_keys = difflog['diff_keys'] or []
            diffs = difflog['diffs'] or {}
            diff_keys.append('substance_id')
            diffs['substance_id'] = [a.reagent.substance_id, 
                                     b.reagent.substance_id]
            difflog['diff_keys'] = diff_keys
            difflog['diffs'] = diffs
            
            activity = (b.reagent.library_contents_version
                         .library_contents_loading_activity.activity)
            log = ApiLog()
            log.username = activity.performed_by.ecommons_id
            log.user_id = activity.performed_by.user.id 
            log.date_time = make_aware(activity.date_created,timezone.get_default_timezone())
            log.ref_resource_name = self.wellResource._meta.resource_name
            # TODO: what types here? could also be a REST specifier, i.e. 'PATCH'
            log.api_action = 'MIGRATION'
            log.uri = self.wellResource.get_resource_uri(bundle_or_obj=well)
            log.key = '/'.join([str(x) for x in (
                self.wellResource.detail_uri_kwargs(well).values()) ])
            log.diff_dict_to_api_log(difflog)
            
            log.json_field = json.dumps({
                'version': b.reagent.library_contents_version.version_number })
            log.save()
            return log
        return None
    
#     def diff_versions(self, wells, diff_function):
#         print 'diff versions by well: ', len(wells)
#         
#         logs_created = 0
#         skipped_log = 0
#         for well, smrs in wells.items():
#             b = None
#             a = None
#             for smr in smrs:
#                 b = a
#                 a = smr
#                 if not b:
#                     continue
#                 difflog = diff_function(a,b)
#                 if ( 'added_keys' in difflog or 
#                      'removed_keys' in difflog or
#                      'diff_keys' in difflog ):
#                     diff_keys = difflog['diff_keys'] or []
#                     diffs = difflog['diffs'] or {}
#                     diff_keys.append('substance_id')
#                     diffs['substance_id'] = [a.reagent.substance_id, 
#                                              b.reagent.substance_id]
#                     difflog['diff_keys'] = diff_keys
#                     difflog['diffs'] = diffs
#                     
#                     activity = (b.reagent.library_contents_version
#                                  .library_contents_loading_activity.activity)
#                     log = ApiLog()
#                     log.username = activity.performed_by.ecommons_id
#                     log.user_id = activity.performed_by.user.id 
#                     log.date_time = activity.date_created
#                     log.ref_resource_name = self.wellResource._meta.resource_name
#                     # TODO: what types here? could also be a REST specifier, i.e. 'PATCH'
#                     log.api_action = 'MIGRATION'
#                     log.uri = self.wellResource.get_resource_uri(bundle_or_obj=well)
#                     log.key = '/'.join([str(x) for x in (
#                         self.wellResource.detail_uri_kwargs(well).values()) ])
#                     log.diff_dict_to_api_log(difflog)
#                     
#                     log.json_field = json.dumps({
#                         'version': b.reagent.library_contents_version.version_number })
#                     log.save()
#                     logs_created += 1
#                 else:
#                     skipped_log += 1
#         print 'logs created: ', logs_created, 'skipped_log', skipped_log
            
