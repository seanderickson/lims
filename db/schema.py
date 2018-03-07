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

