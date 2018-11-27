from __future__ import unicode_literals 


class APP_PUBLIC_DATA:
    ''' Include public data about the application'''

    # NOTE: Override this file in production, and include in the settings.py 

    @classmethod
    def as_dict(cls):
        _vars = vars(APP_PUBLIC_DATA)
        _vars.update(vars(cls))
        _vars = { k:v for k, v in _vars.iteritems() 
            if '__' not in k 
                and hasattr(cls, k)
                and callable(getattr(cls, k)) is not True }
        return _vars

    app_name = 'Screensaver'
    facility_name = 'Screening Facility'
    facility_url = 'http://facility.website.address.here/'
    site_url = 'https://website.address.here/lims'
    contact_feedback_name = 'Feedback'
    contact_feedback_email = 'feedback@email.address.here'
    contact_facility_director_name = 'Facility Director'
    contact_facility_director_email = 'facility_director@email.address.here'
    contact_site_admin_name = 'Site Admin'
    contact_site_admin_email = 'site_admin@email.address.here'
    contact_informatics_name = 'Informatics Contact'
    contact_informatics_email = 'informatics_contact@email.address.here'

    software_development_facility = 'ICCB-Longwood Screening Facility'
    software_development_facility_url = 'http://iccb.med.harvard.edu/'
    software_repository_url = 'https://github.com/hmsiccbl/lims'
    
    page_title = 'Screensaver'
    login_help = 'If you are an approved Screensaver user, you may use your username and password to log in.'
    login_contact_help = '<a href="mailto:%s">Account request, log in help, and Screensaver queries</a>' % (contact_feedback_email)
    site_disclaimer = \
        '''
        Information in this database is confidential and is to be shared
        only among members of the [facility_name] screening
        community. By logging into this database, I am agreeing to hold in
        confidence all information that I learn, download, or print until
        the information is publicly available. Thus, deposition of
        information into this database does not constitute a public
        disclosure and those who deposit information, including myself, can
        preserve their ability to publish and patent the results of their
        work if they so choose.'''

    small_molecule_cherry_pick_ratio_allowed = 0.003
    rnai_cherry_pick_ratio_allowed = 0.0165

    # Interval between user activity and logout for the UI
    BROWSER_IDLE_TIMEOUT_SECONDS = 2 * 60 * 60  # 2 hours
