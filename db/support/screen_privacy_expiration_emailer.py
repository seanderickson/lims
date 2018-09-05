from __future__ import unicode_literals

import argparse
from collections import OrderedDict
import datetime
import json
import logging
import os.path
from urlparse import urlparse

import dateutil.parser
from django.conf import settings
from requests.packages import urllib3

from db.schema import get_href, get_vocab_title, replace_vocabularies, \
    replace_html_values, DB_API_URI, DATE_FORMAT
import db.schema as SCHEMA
from reports import HEADER_APILOG_COMMENT_CLIENT
from reports.serialize import LimsJSONEncoder, parse_val
from reports.utils import parse_credentials, sort_nicely
from reports.utils.admin_emailer import Emailer, read_email_template, \
    create_prettytable, validate_email
from reports.utils.django_requests import get_resource, get_resource_listing, \
    get_resource_schema
import reports.utils.django_requests as django_requests


logger = logging.getLogger(__name__)

SCREEN = SCHEMA.SCREEN
USER = SCHEMA.SCREENSAVER_USER
PUB = SCHEMA.PUBLICATION
FIELD_DSL = SCREEN.DATA_SHARING_LEVEL
FIELD_DPED = SCREEN.DATA_PRIVACY_EXPIRATION_DATE
FIELD_NOTIFICATION_DATE = SCREEN.DATA_PRIVACY_EXPIRATION_NOTIFIED_DATE
FIELD_MIN_DPED = SCREEN.MIN_DATA_PRIVACY_EXPIRATION_DATE
FIELD_MAX_DPED = SCREEN.MAX_DATA_PRIVACY_EXPIRATION_DATE

SCREEN_STATUS_IGNORE =['dropped_technical','dropped_resources',
    'transferred_to_broad_institute']

VOCAB = SCHEMA.VOCAB
VOCAB_DSL = VOCAB.screen.data_sharing_level
VOCAB_USER_ROLE = VOCAB.screen.user_role
VOCAB_SCREEN_RESULT_AVAILABILITY = VOCAB.screen.screen_result_availability

DEFAULT_APILOG_COMMENT = 'Automated user agreement action'

EMAIL_MESSAGE_TEMPLATES = {
    'html_wrapper': 'email_template.html',
    'msg_admin_adjust_data_privacy_expiration_date':
        'msg_admin_adjust_data_privacy_expiration_date.txt',
    'msg_admin_adjust_data_privacy_expiration_no_action':
        'msg_admin_adjust_data_privacy_expiration_no_action.txt',
    'msg_screen_member_notification': 'msg_screen_member_notification.txt',
    'msg_admin_screens_notified': 'msg_admin_screens_notified.txt',
    'msg_admin_no_notifications': 'msg_admin_no_notifications.txt',
    'msg_admin_screens_expired': 'msg_admin_screens_expired.txt',
    'msg_admin_no_expirations': 'msg_admin_no_expirations.txt',
    'msg_admin_no_published_not_shared':'msg_admin_no_published_not_shared.txt',
    'msg_admin_published_not_shared':'msg_admin_published_not_shared.txt',
    }

parser = argparse.ArgumentParser(description=
    'Expire or Notify users with user agreements older than the cutoff dates')

# Email setup options
parser.add_argument(
    '-email_message_directory', '--email_message_directory', required=True,
    help='Email message directory, contains: %r' 
        % [str(s) for s in EMAIL_MESSAGE_TEMPLATES.values()] )
parser.add_argument(
    '-email_log_filename', required=False, 
    help='print email messages to log file before sending')
parser.add_argument(
    '-no_email', '--no_email', action='store_true',
    help='do not send any email notifications (to users or admins)')

parser.add_argument(
    '-admin_from_email', required=False,
    help='If specified, used for the "From:" header; otherwise, defaults to the '
    '"admin_email_address specified or from the admin\'s user acct.')
parser.add_argument(
    '-admin_email', required=False,
    help='If specified, used for the "To:" header; otherwise, defaults to the '
    '"admin_email_address specified or from the admin\'s user acct.')
parser.add_argument(
    '-mail_recipient_list', nargs='+',
    help='additional admin recipients, comma separated list of email addresses')
        

parser.add_argument('-contact_info', '--contact_info', required=True, 
    help='Contact information, e.g. "Administrator user (admin@lims.edu)')

# Admin credentials
parser.add_argument(
    '-U', '--username',
    help='username for the api authentication')
parser.add_argument(
    '-p', '--password',
    help='user password for the api authentication')
parser.add_argument(
    '-c', '--credential_file',
    help = 'credential file containing the username:password for api authentication')

# Testing flow options
parser.add_argument(
    '-test_only', '--test_only', action='store_true',
    help='test run only: no user emails are sent, database changes are rolled back')
parser.add_argument(
    '-admin_email_only', '--admin_email_only', action='store_true',
    help='perform actions and send email to the admin account only, no user emails are sent.')

# Query/Expire/Notify Options
parser.add_argument(
    '-screen_type', '--screen_type', required=True,
    choices=('sm', 'rnai'),
    help='User Agreement for Screen Type (sm="Small Molecule", rnai="RNAi"')
parser.add_argument(
    '-adjust_expiration_days_from_activity', type=int,
    help='Set the screen data_privacy_expiration_date by the given number of '
        '(days) from the date_of_last_library_screening')
parser.add_argument(
    '-days_ahead_to_notify', '--days_ahead_to_notify', type=int,
    help='number of days ahead of expiration to notify users; '
    'Note: this option is exclusive: if this'
    'option is used, notification only will be performed')
parser.add_argument(
    '-expire', action='store_true',
    help='find and expire screens having a data_privacy_expiration_date on or '
    'before the current date. Set: data_sharing_level to 1. '
    'Note: this option is exclusive: if this'
    'option is used, expire only will be performed')
parser.add_argument(
    '-notifyofpublications', action='store_true',
    help='Find private screens that have been published.'
    'Send a notification message to the admins')

parser.add_argument(
    '-u', '--url', required=True,
    help='server base url, e.g.: "http://localhost:8000')


parser.add_argument(
    '-comment', '--comment', 
    help='Comment to annotate server logs for this action', 
    default=DEFAULT_APILOG_COMMENT)

parser.add_argument(
    '-v', '--verbose', dest='verbose', action='count',
    help="Increase verbosity (specify multiple times for more)")    


if __name__ == "__main__":
    args = parser.parse_args()
    log_level = logging.WARNING # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
        DEBUG=True
    logging.basicConfig(
        level=log_level, 
        format='%(msecs)d:%(module)s:%(lineno)d:%(levelname)s: %(message)s')        

    exclusive_args = [ 
        'adjust_expiration_days_from_activity',
        'days_ahead_to_notify', 'expire', 'notifyofpublications' ]
    found_args = [x for x in exclusive_args if getattr(args,x)]
    if not found_args:
        parser.error('Must specify an action: (%s)' % ', '.join(exclusive_args))
    if len(found_args) > 1:
        parser.error('May only specify one of: (%s) at a time' % ', '.join(found_args))
                
    username = None
    password = None
    
    if args.credential_file:
        username,password = parse_credentials(args.credential_file)
    if username is None:
        username = args.username
        if username is None:
            parser.error(
                'username is required if not specifying the credential_file')
        password = args.password
        if not password:
            password = getpass.getpass()

    # Retrieve the email message templates
    html_email_wrapper = None
    for k,v in EMAIL_MESSAGE_TEMPLATES.items():
        _file = os.path.join(args.email_message_directory,v)
        with open(_file) as f:
            if k == 'html_wrapper':
                html_email_wrapper = f.read()
            else:
                EMAIL_MESSAGE_TEMPLATES[k] = read_email_template(f)

    # Log in using django form-based auth, and keep the session
    u = urlparse(args.url)
    base_url = '%s://%s' % (u.scheme,u.netloc)
    # FIXME: suppressing ssl warnings:
    urllib3.disable_warnings()
    session = django_requests.get_logged_in_session(
        username, password, base_url)
    
    # Get admin user data
    headers ={}
    headers['Accept'] = 'application/json'
    admin_user_url = '/'.join([base_url,DB_API_URI, USER.resource_name,username])
    admin_user_data = get_resource(admin_user_url, session, headers)
    logger.info('admin user: %r', admin_user_data)
    
    admin_email = admin_user_data[USER.EMAIL]
    if args.admin_email:
        admin_email = args.admin_email
    
    admin_from_email = admin_email
    if args.admin_from_email:
        admin_from_email = args.admin_from_email

    admin_email_list = [admin_email]
    if args.mail_recipient_list:
        admin_email_list.extend(args.mail_recipient_list)

    # Set up the emailer        
    emailer = Emailer(admin_email, admin_from_email_address=admin_from_email,
        no_email=args.no_email, admin_email_only=args.admin_email_only, 
        test_only=args.test_only, html_email_wrapper=html_email_wrapper,
        email_log_filename=args.email_log_filename)
    
    # Get metadata - schemas
    # Get the screen schema
    schema_url = '/'.join([base_url, DB_API_URI, SCREEN.resource_name, 'schema'])
    screen_schema = get_resource_schema(schema_url, session, headers)
    
    # Get the user schema
    schema_url = '/'.join([base_url, DB_API_URI, USER.resource_name, 'schema'])
    user_schema = get_resource_schema(schema_url, session, headers)
    
    # Specify some non-schema fields for the reports
    FIELD_NEW_DPED = 'new_%s' % FIELD_DPED
    PSEUDO_FIELD_MEMBERS = 'members'
        
    screen_type = args.screen_type
    if args.screen_type == 'sm':
        screen_type = 'small_molecule'
    
    # Screen Schema utils
    def get_title(key, schema):
        
        if key == FIELD_NEW_DPED: 
            return 'New %s' % schema['fields'][FIELD_DPED]['title']
        elif key == PSEUDO_FIELD_MEMBERS:
            return 'Members'
        else:
            return SCHEMA.get_title(key, schema)
    
    def get_role(user_id,screen):
        role_field = 'screensaver_user_role'
        if user_id == screen[SCREEN.LAB_HEAD_ID]:
            return get_vocab_title(
                role_field, VOCAB_USER_ROLE.PRINCIPAL_INVESTIGATOR, screen_schema)
        elif user_id == screen[SCREEN.LEAD_SCREENER_ID]:
            return get_vocab_title(
                role_field, VOCAB_USER_ROLE.LEAD_SCREENER, screen_schema)
        else:
            return get_vocab_title(
                role_field, VOCAB_USER_ROLE.COLLABORATOR, screen_schema)

    def get_member_ids(screen):
        members = [screen[SCREEN.LEAD_SCREENER_ID]]
        members.append(screen[SCREEN.LAB_HEAD_ID])
        if screen.get(SCREEN.COLLABORATORS_ID, None) is not None:
            # FIXME: collaborator_ids should be a list of integers
            members.extend([int(x) for x in screen[SCREEN.COLLABORATORS_ID]])
        return members
    
    def get_users_for_screens(screens):

        # Retrieve all of the user information
        user_ids_to_retrieve = set()
        for screen in screens.values():
            user_ids_to_retrieve.update(get_member_ids(screen))
        logger.info('user_ids_to_retrieve: %d, %r', 
            len(user_ids_to_retrieve), sorted(user_ids_to_retrieve))

        user_url = '/'.join([base_url,DB_API_URI,USER.resource_name])
        params = {
            '%s__in'%USER.SCREENSAVER_USER_ID : user_ids_to_retrieve,
            'includes': '*', 
            'limit': 0
            }
        users,meta = get_resource_listing(user_url, session, headers, params)
        users = { u[USER.SCREENSAVER_USER_ID]:u for u in users }

        return users

    def get_publications_for_screens(screens):
        pub_ids_to_retrieve = set()
        for screen in screens.values():
            pub_ids_to_retrieve.update(screen.get(SCREEN.PUBLICATION_IDS,[]))
        logger.info('pub_ids_to_retrieve: %d, %r', 
            len(pub_ids_to_retrieve), sorted(pub_ids_to_retrieve))

        pub_url = '/'.join([base_url,DB_API_URI,PUB.resource_name])
        params = {
            '%s__in' % PUB.PUBLICATION_ID : pub_ids_to_retrieve,
            'includes': '*', 
            'limit': 0
            }
        pubs,meta = get_resource_listing(pub_url, session, headers, params)
        pubs = { p[PUB.PUBLICATION_ID]:p for p in pubs }

        return pubs
        
        
    users_text_format = '({role}) #{%s} {%s} {%s} <{%s}>' % (
        USER.SCREENSAVER_USER_ID, USER.FIRST_NAME,
        USER.LAST_NAME, USER.EMAIL)
    users_html_format = '<b>{role}</b><br/>{%s} {%s} {%s}<br/>{%s}' % (
        USER.SCREENSAVER_USER_ID, USER.FIRST_NAME,
        USER.LAST_NAME, USER.EMAIL)
    
    # Default args for all queries:
    # - exclude studies
    # - have screen results
    # - DSL is 2 (overlapping) or 3 (private)    
    # - exclude dropped and transferred screens
    default_report_args = {
        '%s__is_null' % SCREEN.STUDY_TYPE : True,
        SCREEN.SCREEN_TYPE: screen_type,
        '%s__in' % SCREEN.SCREEN_RESULT_AVAILABILITY: 
            VOCAB.screen.screen_result_availability.AVAILABLE,
        '%s__in' % FIELD_DSL:  
            [VOCAB_DSL.MUTUAL_POSITIVES,VOCAB_DSL.PRIVATE],
        '-%s__in' % SCREEN.STATUS: SCREEN_STATUS_IGNORE,
        'includes': [
            FIELD_MIN_DPED,
            FIELD_MAX_DPED,
            FIELD_DPED,
            FIELD_NOTIFICATION_DATE,
            SCREEN.LAST_LIBRARY_SCREENING_DATE,
            SCREEN.COLLABORATORS_ID,
            '-%s' % SCREEN.STATUS_DATE ],
        'limit': 0
        }

    screen_url = '/'.join([base_url,DB_API_URI, SCREEN.resource_name])
    current_time = datetime.datetime.now()
    # NOTE: time/time zone will not be transmitted by the script (dates only):
    # all date calculations are in the local time zone; assumed to be the same
    # as the server
    # current_time = pytz.timezone('US/Eastern').localize(current_time)
    
    if args.notifyofpublications:

        report_args = default_report_args.copy()
        # Query for the screens with attached publications that are listed as private
        report_args['%s__is_null' % SCREEN.PUBLICATION_IDS] = False
        # Find DSL > 0
        report_args['%s__in' % FIELD_DSL] = \
            [VOCAB_DSL.MUTUAL_POSITIVES,VOCAB_DSL.PRIVATE, VOCAB_DSL.MUTUAL]

        screens,meta = get_resource_listing(screen_url, session, headers, report_args)
        screens = { s[SCREEN.FACILITY_ID]:s for s in screens }
        if not screens:
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES[
                    'msg_admin_no_published_not_shared']
            msg_subject = msg_subject.format(screen_type=
                get_vocab_title(SCREEN.SCREEN_TYPE, screen_type, screen_schema))
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines)
        else:
            
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_published_not_shared']
            
            msg_params = {
                'count': len(screens), 
                'screen_type': get_vocab_title(
                    SCREEN.SCREEN_TYPE, screen_type, screen_schema),    
                }
            msg_subject = msg_subject.format(**msg_params)
            msg_body_lines = [txt.format(**msg_params) for txt in msg_body_lines]
            txt_msg_body_lines = msg_body_lines
            html_msg_body_lines = msg_body_lines[:]
            
            # Create the report URL
            report_title = '%s Screen Publication Report for %s'\
                % (get_vocab_title(SCREEN.SCREEN_TYPE, screen_type, screen_schema), 
                    current_time.strftime(DATE_FORMAT))
            txt_msg_body_lines.append(report_title)
            report_resource_subtype = { 
                'small_molecule': 'small_molecule_screens',
                'rnai': 'rnai_screens' }[screen_type]
            report_url = '/'.join([
                settings.APP_PUBLIC_DATA.site_url,'#%s'%report_resource_subtype])
            
            # TODO: refactor report url generation
            report_url += '/includes/' + ','.join([
                SCREEN.PUBLICATIONS,SCREEN.LAST_LIBRARY_SCREENING_DATE,
                SCREEN.SCREEN_RESULT_AVAILABILITY,
                ])
            search_args = []
            for k,v in report_args.items():
                if k in ['includes','limit']:
                    continue
                if isinstance(v,(list,tuple)):
                    search_args.append('%s=%s' % (k,','.join(map(str,v))))
                else:
                    search_args.append('%s=%s' % (k,str(v)))
                    
            search_args = ';'.join(search_args)
            report_url += '/search/' + search_args
            txt_msg_body_lines.append(report_url)
            html_msg_body_lines.append('<a href="%s">%s</a>' % (
                report_url, report_title ))
            
            report_fields = [
                SCREEN.FACILITY_ID, SCREEN.TITLE, FIELD_DSL, FIELD_DPED,
                SCREEN.LAST_LIBRARY_SCREENING_DATE, SCREEN.PUBLICATIONS,
                PSEUDO_FIELD_MEMBERS]
            table_fields = OrderedDict(
                (key,get_title(key,screen_schema)) for key in report_fields)
            users = get_users_for_screens(screens)
            publications = get_publications_for_screens(screens)
            # Create text and html versions of the screens
            def replace_for_text(screen):
                screen = replace_vocabularies(screen, screen_schema)
                members = []
                for user_id in get_member_ids(screen):
                    user = users[user_id].copy()
                    user['role'] = get_role(user[USER.SCREENSAVER_USER_ID],screen)
                    text = users_text_format.format(**user)
                    members.append(text)
                screen[PSEUDO_FIELD_MEMBERS] = '\n'.join(members)
                pubs_text = []
                for pub_id in screen.get(SCREEN.PUBLICATION_IDS):
                    pub = publications[int(pub_id)]
                    pubs_text.append(PUB.format_publication(pub))
                screen[SCREEN.PUBLICATIONS] = '\n'.join(pubs_text)
                
                return screen
    
            def replace_for_html(screen):
                screen = replace_html_values(screen, screen_schema)
                members = []
                for user_id in get_member_ids(screen):
                    user = users[user_id].copy()
                    user['role'] = get_role(user[USER.SCREENSAVER_USER_ID],screen)
                    user[USER.SCREENSAVER_USER_ID] = get_href(
                        '#{}'.format(user[USER.SCREENSAVER_USER_ID]), 
                        USER.SCREENSAVER_USER_ID, user_schema,user)
                    text = users_html_format.format(**user)
                    members.append(text)
                screen[PSEUDO_FIELD_MEMBERS] = '<br/>'.join(members)
                pubs_text = []
                for pub_id in screen.get(SCREEN.PUBLICATION_IDS):
                    pub = publications[int(pub_id)]
                    pubs_text.append(PUB.format_publication(pub))
                screen[SCREEN.PUBLICATIONS] = '<br/>'.join(pubs_text)
                return screen
            
            txt_screens = []
            html_screens = []
            for facility_id in sort_nicely(screens.keys()):
                screen = screens[facility_id]
                txt_screens.append(
                    replace_for_text(screen))
                html_screens.append(
                    replace_for_html(screen))
                
            pt = create_prettytable(txt_screens, table_fields)
            txt_msg_body_lines.append(pt.get_string(padding_width=5))
            pt = create_prettytable(html_screens, table_fields)
            html_msg_body_lines.append(pt.get_html_string(
                attributes = {"class": "content-table"})
                .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
    
            emailer.send_email(admin_email_list, msg_subject, txt_msg_body_lines,
                html_msg_body_lines=html_msg_body_lines)
            
        
    
    if args.adjust_expiration_days_from_activity:
        
        report_args = default_report_args.copy()
        
        # 20180419 - fixed; consider all screens for adjustment; not just:
        # Query for the screens with a data_privacy_expiration date less than
        # the cutoff specified by days from last library screening
        # date_of_last_activity_cutoff = current_time + datetime.timedelta(
        #    days=-args.adjust_expiration_days_from_activity)
        # report_args['%s__lt' % SCREEN.LAST_LIBRARY_SCREENING_DATE] = \
        #     date_of_last_activity_cutoff.strftime(DATE_FORMAT)

        screens,meta = get_resource_listing(screen_url, session, headers, report_args)
        if not screens:
            # TODO: this condition should not happen; consider checking for 
            # no patches and no overrides to send this message
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES[
                    'msg_admin_adjust_data_privacy_expiration_no_action']
            msg_subject = \
                msg_subject.format(days_ahead_to_notify=args.days_ahead_to_notify)
            msg_body_lines = [
                txt.format(days_ahead_to_notify=args.days_ahead_to_notify) 
                    for txt in msg_body_lines]
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines)
        else:

            screen_patches = []
            screens_unchanged = []
            screens_set_as_usual = []
            screens_set_to_min = []
            screens_set_to_max = []
            
            # Iterate the screens and create patches to the dped
            for screen in screens:
                facility_id = screen[SCREEN.FACILITY_ID]
                screen[FIELD_NEW_DPED] = None
                screen_patch = { SCREEN.FACILITY_ID: facility_id }
                date_of_last_activity = screen[SCREEN.LAST_LIBRARY_SCREENING_DATE]
                if date_of_last_activity:
                    date_of_last_activity = dateutil.parser.parse(date_of_last_activity)
                current_dped = screen[FIELD_DPED]
                if current_dped:
                    current_dped = dateutil.parser.parse(current_dped)
                max_allowed_dped = screen[FIELD_MAX_DPED]
                if max_allowed_dped:
                    max_allowed_dped = dateutil.parser.parse(max_allowed_dped)
                min_allowed_dped = screen[FIELD_MIN_DPED]
                if min_allowed_dped:
                    min_allowed_dped = dateutil.parser.parse(min_allowed_dped)
                
                new_dped = date_of_last_activity + datetime.timedelta(
                    days=args.adjust_expiration_days_from_activity)
                
                if min_allowed_dped and new_dped < min_allowed_dped:
                    new_dped = min_allowed_dped
                    screens_set_to_min.append(screen)
                elif max_allowed_dped and new_dped > max_allowed_dped:
                    new_dped = max_allowed_dped
                    screens_set_to_max.append(screen)
                elif new_dped != current_dped:
                    screens_set_as_usual.append(screen)
                            
                if new_dped != current_dped:
                    screen_patch[FIELD_DPED] = new_dped.strftime(DATE_FORMAT)
                    screen[FIELD_NEW_DPED] = new_dped.strftime(DATE_FORMAT)
                    screen_patch[FIELD_NOTIFICATION_DATE] = None
                    screen_patches.append(screen_patch)
                    logger.info('%s, patch: %s to %s', facility_id, 
                        current_dped.strftime(DATE_FORMAT), 
                        new_dped.strftime(DATE_FORMAT))
                else:
                    screens_unchanged.append(screen)
                    logger.info('unchanged: %s', facility_id)
                    
            if screen_patches:
                logger.info('patch %d screens', len(screen_patches))
                _url = screen_url
                if args.test_only is True:
                    _url += '?test_only=true'
                headers = {
                    'Content-Type': 'application/json',
                    'Accept': 'application/json',
                    HEADER_APILOG_COMMENT_CLIENT: args.comment
                }
                data = json.dumps(
                    screen_patches, cls=LimsJSONEncoder,
                    encoding="utf-8")
                r = django_requests.patch(
                    _url, session,
                    data = data,
                    headers = headers)
                if r.status_code not in [200]:
                    content = json.loads(r.content)
                    if args.test_only is True:
                        logger.info('"test_only" run: response: %r', content)
                    else:
                        raise Exception("Error: status: %s, %s" 
                            % (r.status_code, content))
                content = json.loads(r.content)
                logger.info('PATCH result: %r', content)
                logger.info('content: %r', content.keys())
                logger.info('meta: %r', content.get(SCHEMA.API_RESULT_META,None))

            # send admin email
            def fill_parms(txt):
                ''' fill text format parameters for messages '''
                return txt.format(
                    adjust_expiration_days_from_activity = \
                        args.adjust_expiration_days_from_activity,
                    screens_set_as_usual = len(screens_set_as_usual),
                    screens_set_to_min = len(screens_set_to_min),
                    screens_set_to_max = len(screens_set_to_max),
                    patch_count = len(screen_patches),
                    unchanged_count = len(screens_unchanged),
                    screens_patched = \
                        ', '.join([s[SCREEN.FACILITY_ID] for s in screen_patches]),
                    screens_unchanged = \
                        ', '.join([s[SCREEN.FACILITY_ID] for s in screens_unchanged])
                    )
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_adjust_data_privacy_expiration_date']
            msg_subject = fill_parms(msg_subject)
            msg_body_lines = [fill_parms(txt) for txt in msg_body_lines]
            html_msg_body_lines = msg_body_lines[:]
    
            # Create the report URL
            report_title = '%s Screen DSL Adjustment Report for %s'\
                % (get_vocab_title(SCREEN.SCREEN_TYPE, screen_type, screen_schema), 
                    current_time.strftime(DATE_FORMAT))
            msg_body_lines.append(report_title)
            report_url = '/'.join([
                settings.APP_PUBLIC_DATA.site_url,'#screen'])
            report_url += '/includes/' + ','.join(report_args.get('includes',''))
            report_screens = set([screen[SCREEN.FACILITY_ID] for screen 
                in screen_patches + screens_set_as_usual + screens_set_to_max
                    + screens_set_to_min])
            report_url += '/search/facility_id__in='\
                 + ','.join(report_screens)
            msg_body_lines.append(report_url)
            html_msg_body_lines.append('<a href="%s">%s</a>' % (
                report_url, report_title ))
                
            # Create the report summary tables
            
            if screens_set_as_usual:
                title = 'Adjusted to %d days from the %s: ' % (
                    args.adjust_expiration_days_from_activity, 
                    get_title(SCREEN.LAST_LIBRARY_SCREENING_DATE, screen_schema))
                msg_body_lines.append('')
                msg_body_lines.append(title)
                html_msg_body_lines.append('<hr>')
                html_msg_body_lines.append('<h4>%s</h4><br/>'%title)
    
                report_fields = [
                    SCREEN.FACILITY_ID, SCREEN.TITLE, FIELD_DSL, FIELD_DPED,
                    FIELD_NEW_DPED, SCREEN.LAST_LIBRARY_SCREENING_DATE]
                table_fields = OrderedDict(
                    (key,get_title(key, screen_schema)) for key in report_fields)
                report_screens = [ replace_vocabularies(screen,screen_schema) 
                    for screen in screens_set_as_usual]
                pt = create_prettytable(report_screens, table_fields)
                msg_body_lines.append(pt.get_string(padding_width=5))
                report_screens = [ replace_html_values(screen,screen_schema) 
                    for screen in screens_set_as_usual]
                pt = create_prettytable(report_screens, table_fields)
                html_msg_body_lines.append(pt.get_html_string(
                    attributes = {"class": "content-table"})
                    .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
    
            if screens_set_to_min:
                title = 'Overridden, set to the %s: ' \
                    % get_title(FIELD_MIN_DPED, screen_schema)
                msg_body_lines.append('')
                msg_body_lines.append(title)
                html_msg_body_lines.append('<hr>')
                html_msg_body_lines.append('<h4>%s</h4><br/>'%title)
    
                report_fields = [
                    SCREEN.FACILITY_ID, SCREEN.TITLE, FIELD_DSL, FIELD_DPED,
                    FIELD_NEW_DPED, FIELD_MIN_DPED, 
                    FIELD_NOTIFICATION_DATE, SCREEN.LAST_LIBRARY_SCREENING_DATE]
                table_fields = OrderedDict(
                    (key,get_title(key,screen_schema)) for key in report_fields)
                report_screens = [ replace_vocabularies(screen,screen_schema) 
                    for screen in screens_set_to_min]
                pt = create_prettytable(report_screens, table_fields)
                msg_body_lines.append(pt.get_string(padding_width=5))
                report_screens = [ replace_html_values(screen,screen_schema) 
                    for screen in screens_set_to_min]
                pt = create_prettytable(report_screens, table_fields)
                html_msg_body_lines.append(pt.get_html_string(
                    attributes = {"class": "content-table"})
                    .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
                
            if screens_set_to_max:
    
                title = 'Overridden, set to the %s:' \
                    % get_title(FIELD_MIN_DPED,screen_schema)
                msg_body_lines.append('')
                msg_body_lines.append(title)
                html_msg_body_lines.append('<hr>')
                html_msg_body_lines.append('<h4>%s</h4><br/>'%title)
                
                report_fields = [
                    SCREEN.FACILITY_ID, SCREEN.TITLE, FIELD_DSL, FIELD_DPED,
                    FIELD_NEW_DPED, FIELD_MAX_DPED, 
                    FIELD_NOTIFICATION_DATE, SCREEN.LAST_LIBRARY_SCREENING_DATE]
                table_fields = OrderedDict(
                    (key,get_title(key,screen_schema)) for key in report_fields)
                report_screens = [ replace_vocabularies(screen,screen_schema) 
                    for screen in screens_set_to_max]
                pt = create_prettytable(report_screens, table_fields)
                msg_body_lines.append(pt.get_string(padding_width=5))
                report_screens = [ replace_html_values(screen,screen_schema) 
                    for screen in screens_set_to_max]
                pt = create_prettytable(report_screens, table_fields)
                html_msg_body_lines.append(pt.get_html_string(
                    attributes = {"class": "content-table"})
                    .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
                
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines,
                html_msg_body_lines=html_msg_body_lines )
        
    if args.expire:

        report_args = default_report_args.copy()
        report_args['%s__lte' % FIELD_DPED] = \
            current_time.strftime(DATE_FORMAT)
        screens_to_expire,meta = get_resource_listing(
            screen_url, session, headers, report_args)

        if not screens_to_expire:
            logger.info('No screens found to expire, send email for no expirations...')
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_no_expirations']
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines)
            
        else:
        
            screens_to_expire = { s[SCREEN.FACILITY_ID]:s 
                for s in screens_to_expire }
            logger.info('retrieved %d screens with %s <= %s', 
                len(screens_to_expire), FIELD_DPED, 
                current_time.strftime(DATE_FORMAT))
            
            screen_patches = []
            for facility_id, screen in screens_to_expire.items():
                screen_patch = { SCREEN.FACILITY_ID: facility_id }
                screen_patch[FIELD_DSL] = VOCAB_DSL.MUTUAL
                
                logger.info(('expire: {%s}: {%s}'
                    %(SCREEN.FACILITY_ID,FIELD_DPED)).format(**screen))
                
                screen_patches.append(screen_patch)
                
            logger.info('expire patches: %r', screen_patches)    
            _url = screen_url
            if args.test_only is True:
                _url += '?test_only=true'
            headers = {
                'Content-Type': 'application/json',
                'Accept': 'application/json',
                HEADER_APILOG_COMMENT_CLIENT: args.comment
            }
            data = json.dumps(
                screen_patches, cls=LimsJSONEncoder,
                encoding="utf-8")
            r = django_requests.patch(
                _url, session,
                data = data,
                headers = headers)
            if r.status_code not in [200]:
                content = json.loads(r.content)
                if args.test_only is True:
                    logger.info('"test_only" run: response: %r', content)
                else:
                    raise Exception("Error: status: %s, %s" 
                        % (r.status_code, content))
            content = json.loads(r.content)
            logger.info('PATCH result: %r', content)
            logger.info('content: %r', content.keys())
            logger.info('meta: %r', content.get(SCHEMA.API_RESULT_META,None))
                        
            # Create the admin message
    
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_screens_expired']
            msg_params = {
                'current_date': current_time.strftime(DATE_FORMAT),
                'screen_count': len(screen_patches),
                'data_sharing_level_1': get_vocab_title(
                    SCREEN.DATA_SHARING_LEVEL, VOCAB_DSL.MUTUAL, screen_schema)
            }
            msg_subject = msg_subject.format(**msg_params)
            msg_body_lines = [txt.format(**msg_params) for txt in msg_body_lines]
            
            report_fields = [
                SCREEN.FACILITY_ID, SCREEN.TITLE, FIELD_DSL, FIELD_DPED,
                SCREEN.LAST_LIBRARY_SCREENING_DATE, PSEUDO_FIELD_MEMBERS]
            table_fields = OrderedDict(
                (key,get_title(key,screen_schema)) for key in report_fields)
            users = get_users_for_screens(screens_to_expire)
            # Create text and html versions of the screens
            def replace_for_text(screen):
                screen = replace_vocabularies(screen, screen_schema)
                members = []
                for user_id in get_member_ids(screen):
                    user = users[user_id].copy()
                    user['role'] = get_role(user[USER.SCREENSAVER_USER_ID],screen)
                    text = users_text_format.format(**user)
                    members.append(text)
                screen[PSEUDO_FIELD_MEMBERS] = '\n'.join(members)
                return screen
    
            def replace_for_html(screen):
                screen = replace_html_values(screen, screen_schema)
                members = []
                for user_id in get_member_ids(screen):
                    user = users[user_id].copy()
                    user['role'] = get_role(user[USER.SCREENSAVER_USER_ID],screen)
                    user[USER.SCREENSAVER_USER_ID] = get_href(
                        '#{}'.format(user[USER.SCREENSAVER_USER_ID]), 
                        USER.SCREENSAVER_USER_ID, user_schema,user)
                    text = users_html_format.format(**user)
                    members.append(text)
                screen[PSEUDO_FIELD_MEMBERS] = '<br/>'.join(members)
                return screen
            
            txt_notified_screens = []
            html_notified_screens = []
            for facility_id in sort_nicely(screens_to_expire.keys()):
                screen = screens_to_expire[facility_id]
                txt_notified_screens.append(
                    replace_for_text(screen))
                html_notified_screens.append(
                    replace_for_html(screen))
                
            txt_msg_body_lines = msg_body_lines
            html_msg_body_lines = msg_body_lines[:]
            pt = create_prettytable(txt_notified_screens, table_fields)
            txt_msg_body_lines.append(pt.get_string(padding_width=5))
            pt = create_prettytable(html_notified_screens, table_fields)
            html_msg_body_lines.append(pt.get_html_string(
                attributes = {"class": "content-table"})
                .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
    
            emailer.send_email(admin_email_list, msg_subject, txt_msg_body_lines,
                html_msg_body_lines=html_msg_body_lines)
            
    if args.days_ahead_to_notify:

        report_args = default_report_args.copy()

        # Query for the screens with data_privacy_expiration date
        # DPED < current_date + (days_ahead_to_notify), 
        # and that have not been notified
        notification_cutoff = \
            current_time + datetime.timedelta(days=args.days_ahead_to_notify)
        logger.info('days: %d, notification cutoff date: %s',
            args.days_ahead_to_notify, notification_cutoff.strftime(DATE_FORMAT))
        report_args['%s__lt' % FIELD_DPED] = \
            notification_cutoff.strftime(DATE_FORMAT)
        report_args['%s__is_null' % FIELD_NOTIFICATION_DATE] = True
        
        screens,meta = get_resource_listing(screen_url, session, headers, report_args)
        
        if not screens:
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_no_notifications']
            msg_subject = msg_subject.format(days_ahead_to_notify=args.days_ahead_to_notify)
            msg_body_lines = [
                txt.format(days_ahead_to_notify=args.days_ahead_to_notify) 
                    for txt in msg_body_lines]
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines)
        else:
            screens = {s[SCREEN.FACILITY_ID]:s for s in screens}
            users = get_users_for_screens(screens)
            
            screens_notified = set()
            screens_not_notified = set()
            user_notification_results = dict()
    
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_screen_member_notification']
            args.sample_message_subject = None
            args.sample_message_lines = None
            args.sample_html_message_lines = None
            args.sample_recipient = None
            for facility_id, screen in screens.items():
                logger.info('attempt to notify users for screen: %s', facility_id)
                
                screen_title_format = '#{facility_id}: "{title}"'
                msg_parms = {
                    'screen_title_text': screen_title_format.format(**screen),
                    'data_sharing_level_title': get_vocab_title(
                        FIELD_DSL, screen[FIELD_DSL], screen_schema),
                    'data_sharing_level_1_title': 
                        get_vocab_title(FIELD_DSL, VOCAB_DSL.MUTUAL, screen_schema),
                    'data_privacy_expiration_date': screen[FIELD_DPED],
                    'contact_info': args.contact_info,
                }
                msg_parms_html = msg_parms.copy()
                msg_parms_html['screen_title_text'] = \
                    screen_title_format.format(
                        facility_id = \
                            get_href(screen[SCREEN.FACILITY_ID],SCREEN.FACILITY_ID,
                                screen_schema, screen),
                        title = screen[SCREEN.TITLE])

                def notify_screen(screen,user):
                    ''' send email, True if success'''
                    valid_email = validate_email(user.get(USER.EMAIL, None))
                    if not valid_email:
                        return False, 'Invalid email address'
                    
                    text_msg_body_lines = \
                        [txt.format(**msg_parms) for txt in msg_body_lines]
                    html_msg_body_lines = \
                        [txt.format(**msg_parms_html) 
                            for txt in msg_body_lines]
                    if args.sample_message_subject is None:
                        args.sample_message_subject = 'Subject: %s' % msg_subject
                        args.sample_message_lines = text_msg_body_lines
                        args.sample_html_message_lines = html_msg_body_lines[:]
                        args.sample_recipient = valid_email
                    emailer.send_email([valid_email], msg_subject, 
                        text_msg_body_lines,
                        html_msg_body_lines=html_msg_body_lines)
                    
                    return True, None
                
                # Screen is not "notified" unless the PI (lab head) 
                # or Lead Screener is notified
                # 20180312 -verified with JAS60
                user = users[screen[SCREEN.LAB_HEAD_ID]]
                is_notified, msg = notify_screen(screen, user)
                user_notification_results[user[USER.SCREENSAVER_USER_ID]] = msg
                if is_notified is True:
                    screens_notified.add(facility_id)
                else:
                    logger.info('Could not notify Lab Head for screen #%s: %r',
                        facility_id, user)
                user = users[screen[SCREEN.LEAD_SCREENER_ID]]
                result, msg = notify_screen(screen, user)
                user_notification_results[user[USER.SCREENSAVER_USER_ID]] = msg
                if result is True:
                    screens_notified.add(facility_id)
                    is_notified = result
                else:
                    logger.info('Could not notify Lead Screener for screen #%s: %r',
                        facility_id, user)
                if is_notified is not True:
                    screens_not_notified.add(facility_id)

                if screen.get(SCREEN.COLLABORATORS_ID, None):
                    for user_id in screen[SCREEN.COLLABORATORS_ID]:
                        # FIXME: collaborator_ids is a list of strings, should be ints
                        user = users[int(user_id)]
                        result, msg = notify_screen(screen, user)
                        user_notification_results[user[USER.SCREENSAVER_USER_ID]] = msg

            # TODO: patch the notification date
            if screens_notified:
                screen_patches = []
                for facility_id in screens_notified:
                    screen_patches.append({
                        SCREEN.FACILITY_ID: facility_id,
                        SCREEN.DATA_PRIVACY_EXPIRATION_NOTIFIED_DATE: 
                            current_time.strftime(DATE_FORMAT) })
                logger.info('patch notification date for %d screens', 
                    len(screen_patches))
                _url = screen_url
                if args.test_only is True:
                    _url += '?test_only=true'
                headers = {
                    'Content-Type': 'application/json',
                    'Accept': 'application/json',
                    HEADER_APILOG_COMMENT_CLIENT: args.comment
                }
                data = json.dumps(
                    screen_patches, cls=LimsJSONEncoder,
                    encoding="utf-8")
                r = django_requests.patch(
                    _url, session,
                    data = data,
                    headers = headers)
                if r.status_code not in [200]:
                    content = json.loads(r.content)
                    if args.test_only is True:
                        logger.info('"test_only" run: response: %r', content)
                    else:
                        raise Exception("Error: status: %s, %s" 
                            % (r.status_code, content))
                content = json.loads(r.content)
                logger.info('PATCH result: %r', content)
                logger.info('content: %r', content.keys())
                logger.info('meta: %r', content.get(SCHEMA.API_RESULT_META,None))
                            
            # Create the admin message
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_screens_notified']

            # Create the notification summary tables
            report_fields = [
                SCREEN.FACILITY_ID, SCREEN.TITLE, FIELD_DSL, FIELD_DPED,
                SCREEN.LAST_LIBRARY_SCREENING_DATE, PSEUDO_FIELD_MEMBERS]
            table_fields = OrderedDict(
                (key,get_title(key,screen_schema)) for key in report_fields)
            
            # Create a summary table for screens notified
            
            # Create a copy of the screens in the list, insert display 
            # formatted values (links and vocabularies)
            def replace_for_text(screen):
                # Create a copy of the screen and format values as needed
                screen = replace_vocabularies(screen, screen_schema)
                member_ids = get_member_ids(screen)
                members =[]
                for user_id in member_ids:
                    if user_id not in users:
                        logger.warn('user info not retrieved: %s', user_id)
                        continue
                    user = users[user_id].copy()
                    msg = user_notification_results[user_id]
                    user['role'] = get_role(user[USER.SCREENSAVER_USER_ID],screen)
                    text = users_text_format.format(**user)
                    if msg is not None:
                        text += ', %s' % msg
                    members.append(text)
                screen[PSEUDO_FIELD_MEMBERS] = '\n'.join(members)
                return screen
            
            def replace_for_html(screen):
                # Create a copy of the screen and format values as needed
                member_ids = get_member_ids(screen)
                members =[]
                for user_id in member_ids:
                    if user_id not in users:
                        logger.warn('user info not retrieved: %s', user_id)
                        continue
                    user = users[user_id].copy()
                    msg = user_notification_results[user_id]
                    user['role'] = get_role(user[USER.SCREENSAVER_USER_ID],screen)
                    user[USER.SCREENSAVER_USER_ID] = get_href(
                        '#{}'.format(user[USER.SCREENSAVER_USER_ID]), 
                        USER.SCREENSAVER_USER_ID, user_schema,user)
                    text = users_html_format.format(**user)
                    if msg is not None:
                        text += ', %s' % msg
                    members.append(text)
                screen[PSEUDO_FIELD_MEMBERS] = '<br/>'.join(members)
                screen = replace_html_values(screen, screen_schema)
                return screen
            
            txt_notified_screens = []
            html_notified_screens = []
            for facility_id in sort_nicely(screens_notified):
                screen = screens[facility_id]
                txt_notified_screens.append(replace_for_text(screen))
                html_notified_screens.append(replace_for_html(screen))
                
            txt_not_notified_screens = []
            html_not_notified_screens = []
            for facility_id in sort_nicely(screens_not_notified):
                screen = screens[facility_id]
                txt_not_notified_screens.append(replace_for_text(screen))
                html_not_notified_screens.append(replace_for_html(screen))
            
            msg_parms = {
                'screen_count': len(screens),
                'days_ahead_to_notify': args.days_ahead_to_notify,
                'screens_notified_count': len(screens_notified),
                'screens_not_notified_count': len(screens_not_notified),
            }
            msg_subject = msg_subject.format(**msg_parms)
            txt_msg_body_lines = []
            html_msg_body_lines = []
            for line in msg_body_lines:
                txt_msg_body_lines.append(line.format(**msg_parms))
                html_msg_body_lines.append(line.format(**msg_parms))
                if 'screens_notified_count' in line and txt_notified_screens:
                    pt = create_prettytable(txt_notified_screens, table_fields)
                    txt_msg_body_lines.append(pt.get_string(padding_width=5))
                    pt = create_prettytable(html_notified_screens, table_fields)
                    html_msg_body_lines.append(pt.get_html_string(
                        attributes = {"class": "content-table"})
                        .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
            
                elif 'screens_not_notified_count' in line and txt_not_notified_screens:
                    pt = create_prettytable(txt_not_notified_screens, table_fields)
                    txt_msg_body_lines.append(pt.get_string(padding_width=5))
                    pt = create_prettytable(html_not_notified_screens, table_fields)
                    html_msg_body_lines.append(pt.get_html_string(
                        attributes = {"class": "content-table"})
                        .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
            
            if args.sample_message_subject:
                txt_msg_body_lines.append('[example email]')
                txt_msg_body_lines.append('To: %s' % args.sample_recipient)
                txt_msg_body_lines.append(args.sample_message_subject)
                txt_msg_body_lines.extend(args.sample_message_lines)
    
                html_msg_body_lines.append('<hr>[example email]<hr/>')
                html_msg_body_lines.append('To: %s' % args.sample_recipient)
                html_msg_body_lines.append('')
                html_msg_body_lines.append(args.sample_message_subject)
                html_msg_body_lines.append('')
                html_msg_body_lines.extend(args.sample_html_message_lines)
            
            emailer.send_email(admin_email_list, msg_subject, txt_msg_body_lines,
                html_msg_body_lines=html_msg_body_lines)
            
            
