''' Schema utils:

Define schema constants
- Field names used by the API and intended as constant over versions
- Vocabulary values used by the API and intended as constant over versions

Utilities for parsing the schema
 '''

import logging
from django.conf import settings

import reports.schema
from reports.schema import *

logger = logging.getLogger(__name__)

DB_API_URI = 'db/api/%s' % VERSION

class SCREEN(schema_obj):
    resource_name = 'screen'
    
    FACILITY_ID = 'facility_id'
    SCREEN_TYPE = 'screen_type'
    SCREEN_RESULT_AVAILABILITY = 'has_screen_result'
    STATUS = 'status'
    STATUS_DATE = 'status_date'
    STUDY_TYPE = 'study_type'
    TITLE = 'title'
    LEAD_SCREENER_ID = 'lead_screener_id'
    LAB_HEAD_ID = 'lab_head_id'
    COLLABORATORS_ID = 'collaborator_ids'
    LAB_NAME = 'lab_name'
    DATA_SHARING_LEVEL = 'data_sharing_level'
    DATA_PRIVACY_EXPIRATION_DATE = 'data_privacy_expiration_date'
    DATA_PRIVACY_EXPIRATION_NOTIFIED_DATE = 'data_privacy_expiration_notified_date'
    MIN_DATA_PRIVACY_EXPIRATION_DATE = 'min_allowed_data_privacy_expiration_date'
    MAX_DATA_PRIVACY_EXPIRATION_DATE = 'max_allowed_data_privacy_expiration_date'
    LAST_LIBRARY_SCREENING_DATE = 'date_of_last_library_screening'
    PUBLICATIONS = 'publications'
    PUBLICATION_IDS = 'publication_ids'

class DATA_COLUMN(schema_obj):
    resource_name = 'datacolumn'
    
    NAME = 'name'
    DESCRIPTION = 'description'
    ASSAY_DATA_TYPE = 'assay_data_type'
    DATA_TYPE = 'data_type'
    DECIMAL_PLACES = 'decimal_places'
    REPLICATE_ORDINAL = 'replicate_ordinal'
    TIME_POINT = 'time_point'
    CHANNEL = 'channel'
    ZDEPTH_ORDINAL = 'zdepth_ordinal'
    ASSAY_READOUT_TYPE = 'assay_readout_type'
    DERIVED_FROM_COLUMNS = 'derived_from_columns'
    HOW_DERIVED = 'how_derived'
    IS_FOLLOW_UP_DATA = 'is_follow_up_data'
    IS_DERIVED = 'is_derived'
    ASSAY_PHENOTYPE = 'assay_phenotype'
    POSITIVES_COUNT = 'positives_count'
    STRONG_POSITIVES_COUNT = 'strong_positives_count'
    MEDIUM_POSITIVES_COUNT = 'medium_positives_count'
    WEAK__POSITIVES_COUNT = 'weak_positives_count'
    ORDINAL = 'ordinal'
    TIME_POINT_ORDINAL = 'time_point_ordinal'
    SCREEN_FACILITY_ID = 'screen_facility_id'
    COMMENTS = 'comments'
    TITLE = 'title'
    KEY = 'key'

class SCREEN_RESULT(schema_obj):
    resource_name = 'screenresult'
    
    WELL_ID = 'well_id'
    PLATE_NUMBER = 'plate_number'
    WELL_NAME = 'well_name'
    LIBRARY_SHORT_NAME = 'short_name'
    LIBRARY_WELL_TYPE = 'library_well_type'
    VENDOR_NAME = 'vendor_name'
    VENDOR_ID = 'vendor_identifier'
    ASSAY_CONTROL_TYPE = 'assay_well_control_type'
    SCREEN_FACILITY_ID = 'screen_facility_id'
    #     'screen_title'
    #     'mouseover'
    IS_POSITIVE = 'is_positive'
    CONFIRMED_POSITIVE_VALUE = 'confirmed_positive_value'
    EXCLUDE = 'exclude'


class SCREENSAVER_USER(schema_obj):
    resource_name = 'screensaveruser'
    
    SCREENSAVER_USER_ID = 'screensaver_user_id'
    FIRST_NAME = 'first_name'
    LAST_NAME = 'last_name'
    EMAIL = 'email'

class USER_AGREEMENT(schema_obj):
    resource_name = 'useragreement'
    
    SCREENSAVER_USER_ID = 'screensaver_user_id'
    USERNAME = 'username'
    USER_EMAIL = 'user_email'
    USER_NAME = 'user_name'
    USER_FIRST_NAME = 'user_first_name'
    USER_LAST_NAME = 'user_last_name'
    TYPE = 'type'
    DATA_SHARING_LEVEL = 'data_sharing_level'
    STATUS = 'status'
    DATE_ACTIVE = 'date_active'
    DATE_EXPIRED = 'date_expired'
    DATE_NOTIFIED = 'date_notified'
    
class PUBLICATION(schema_obj):
    resource_name = 'publication'
    
    PUBLICATION_ID = 'publication_id'
    PUBMED_ID = 'pubmed_id'
    PUBMED_CENTRAL_ID = 'pubmed_central_id'
    AUTHORS = 'authors'
    TITLE = 'title'
    JOURNAL = 'journal'
    VOLUME = 'volume'
    YEAR_PUBLISHED = 'year_published'
    PAGES = 'pages'
    ATTACHED_FILE_ID = 'attached_file_id'
    ATTACHED_FILENAME = 'attached_filename'
    SCREEN_FACILITY_ID = 'screen_facility_id'

    @classmethod
    def format_publication(cls,pub):
        txt = ''
        val = pub.get(cls.AUTHORS)
        if val:
            txt += '%s ' % val
        val = pub.get(cls.YEAR_PUBLISHED)
        if val:
            txt += '(%s) ' % val
        val = pub.get(cls.TITLE)
        if val:
            if val[-1] != '.':
                val += '.'
            txt += '%s ' % val
        val = pub.get(cls.JOURNAL)
        if val:
            txt += '%s ' % val
        val = pub.get(cls.VOLUME)
        if val:
            txt += '%s, ' % val
        val = pub.get(cls.PAGES)
        if val:
            txt += '%s.' % val
            
        return txt    


class VOCAB(reports.schema.VOCAB):
    
    class datacolumn(schema_obj):
        class data_type(schema_obj):
            BOOLEAN_POSITIVE = 'boolean_positive_indicator'
            BOOLEAN = 'boolean'
            CONFIRMED_POSITIVE = 'confirmed_positive_indicator'
            DECIMAL = 'decimal'
            PARTITIONED_POSITIVE = 'partition_positive_indicator'
            INTEGER = 'integer'
            NUMERIC = 'numeric'
            STRING = 'string'
            TEXT = 'text'
            
            numeric_types = (INTEGER, DECIMAL, NUMERIC)
    
    class resultvalue(schema_obj):
        class partitioned_positive(schema_obj):
            NP = 0
            W = 1
            M = 2
            S = 3
        class confirmed_positive(schema_obj):
            NT = 0
            I = 1
            FP = 2
            CP = 3
            
    class assaywell(schema_obj):

        class control_type(schema_obj):
            BUFFER = 'buffer'
            EMPTY = 'empty'
            EXPERIMENTAL = 'experimental'
            ASSAY_CONTROL = 'assay_control'
            ASSAY_POSITIVE_CONTROL = 'assay_positive_control'
            LIBRARY_CONTROL = 'library_control'
            OTHER = 'other'
            
    
    class screen(schema_obj):
        # data_sharing_level = namedtuple('data_sharing_level', 
        #     ['SHARED','MUTUAL', 'MUTUAL_POSITIVES','PRIVATE'])\
        #     (0,1,2,3)
        class data_sharing_level(schema_obj):
            SHARED = 0
            MUTUAL = 1
            MUTUAL_POSITIVES = 2
            PRIVATE = 3
        class user_role(schema_obj):
            PRINCIPAL_INVESTIGATOR = 'principal_investigator'
            LEAD_SCREENER = 'lead_screener'
            COLLABORATOR = 'collaborator'
        class screen_result_availability(schema_obj):
            AVAILABLE = 1
            NOT_SHARED = 2
            NONE = 3
    
    class user_agreement(schema_obj):
        
        class type(schema_obj):
            SM = 'sm'
            RNAI = 'rnai'
        
        class status(schema_obj):
            ACTIVE = 'active'
            EXPIRED = 'expired'
            INACTIVE = 'inactive'
