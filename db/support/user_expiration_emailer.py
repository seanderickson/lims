from __future__ import unicode_literals

from argparse import ArgumentError
import argparse
from collections import defaultdict, OrderedDict
import datetime
import json
import logging
import os.path
from urlparse import urlparse

from django.conf import settings
from requests.packages import urllib3

from reports import InformationError, HEADER_APILOG_COMMENT_CLIENT
from reports.api import API_RESULT_DATA, API_RESULT_META
from reports.utils import parse_credentials
from reports.utils.admin_emailer import Emailer, read_email_template, \
    create_prettytable, validate_email
import reports.utils.django_requests as django_requests


logger = logging.getLogger(__name__)

BASE_URI = 'db/api/v1'
USERAGREEMENT_RESOURCE = 'useragreement'
USER_RESOURCE = 'screensaveruser'

TESTING_MODE = True

FIELD_USER_NAME = 'user_name'
FIELD_USER_EMAIL = 'user_email'
FIELD_AGREEMENT_TYPE = 'type'
SCREEN_TYPES = {
    'sm': 'Small Molecule', 'rnai': 'RNAi' }  

DEFAULT_APILOG_COMMENT = 'Automated user agreement action'

EMAIL_MESSAGE_TEMPLATES = {
    'html_wrapper': 'email_template.html',
    'msg_admin_notification':'msg_admin_user_notification.txt',
    'msg_admin_expiration': 'msg_admin_user_expiration.txt',
    'msg_admin_notify_no_actions':'msg_admin_notify_no_actions.txt',
    'msg_admin_expire_no_actions':'msg_admin_expire_no_actions.txt',
    'msg_user_expiration': 'msg_user_expiration.txt'
    }

parser = argparse.ArgumentParser(
    description='Expire or Notify users with user agreements older than the cutoff dates')

# Email setup options
parser.add_argument(
    '-email_message_directory', '--email_message_directory', required=True,
    help='Email message directory, contains: %r' 
        % [str(s) for s in EMAIL_MESSAGE_TEMPLATES.values()] )

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
parser.add_argument(
    '-mr', '--mail_recipient_list', nargs='+',
    help='additional admin recipients, comma separated list of email addresses')

# Testing flow options
parser.add_argument(
    '-test_only', '--test_only', action='store_true',
    help='test run only: no user emails are sent, database changes are rolled back')
parser.add_argument(
    '-no_email', '--no_email', action='store_true',
    help='do not send any email notifications (to users or admins)')
parser.add_argument(
    '-admin_email_only', '--admin_email_only', action='store_true',
    help='perform actions and send email to the admin account only, no user emails are sent.')

# Query/Expire/Notify Options
parser.add_argument(
    '-screen_type', '--screen_type', required=True,
    choices=('sm', 'rnai'),
    help='User Agreement for Screen Type (sm="Small Molecule", rnai="SiRNA"')
parser.add_argument(
    '-days_to_expire', '--days_to_expire', type=int, required=True,
    help='number of days for a user agreement to expire')
parser.add_argument(
    '-days_ahead_to_notify', '--days_ahead_to_notify', type=int,
    help='number of days ahead of expiration to notify users; if this'
    'option is used, notification only will be performed')

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
    CONTENT_TYPES =   { 
        'json': {'content-type': 'application/json'},
        'csv':  {'content-type': 'text/csv'},
        } 
                
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

    #### log in using django form-based auth, and keep the session
    url = args.url
    u = urlparse(url)
    base_url = '%s://%s' % (u.scheme,u.netloc)

    # FIXME: suppressing ssl warnings:
    urllib3.disable_warnings()
    
    session = django_requests.get_logged_in_session(
        username, password, base_url)
    # django session based auth requires a csrf token
    headers ={}
    headers['X-CSRFToken'] = session.cookies['csrftoken']
    headers['Accept'] = 'application/json'
    logger.debug('successfully acquired token: %r', headers)
    
    # Get admin user email
    admin_user_url = '/'.join([url,BASE_URI, USER_RESOURCE,username])
    logger.info('GET: %r', admin_user_url)
    r = session.get(admin_user_url,headers=headers)
    if r.status_code not in [200]:
        raise Exception("Error for %r, status: %s, %s" 
                        % (admin_user_url, r.status_code, r.content))
    
    admin_user_data = json.loads(r.content)
    logger.info('admin_user_data: %r', admin_user_data)
    admin_from_email = admin_user_data['email']
    logger.info('admin_from_email: %r', admin_from_email)
    
    admin_email_list = [admin_from_email]
    if args.mail_recipient_list:
        admin_email_list.extend(args.mail_recipient_list)
        
    emailer = Emailer(admin_from_email, no_email=args.no_email, 
        admin_email_only=args.admin_email_only, 
        test_only=args.test_only, html_email_wrapper=html_email_wrapper)
        
    # Begin processing user agreements
    # Query for the user agreements
    
    current_time = datetime.datetime.now()
    date_active_cutoff = \
        current_time + datetime.timedelta(days=-(args.days_to_expire))
    if args.days_ahead_to_notify is not None:
        date_active_cutoff = \
            date_active_cutoff + datetime.timedelta(
                days=args.days_ahead_to_notify)
    
    date_active_cutoff = date_active_cutoff.strftime("%Y-%m-%d")
    logger.info(
        'days_to_expire: %d, days_ahead_to_notify: %r, date_active_cutoff: %s', 
        args.days_to_expire, args.days_ahead_to_notify, date_active_cutoff)
    report_args = {
        'status__in': 'active',
        FIELD_AGREEMENT_TYPE: args.screen_type,
        'date_active__lte': date_active_cutoff,
        'includes': [FIELD_USER_EMAIL],
        'limit': 0
        }

    if args.days_ahead_to_notify is not None:
        report_args['date_notified__is_null'] = True
    
    useragreement_url = '/'.join([url,BASE_URI, USERAGREEMENT_RESOURCE])
    logger.info('GET: %r', useragreement_url)
    r = session.get(useragreement_url,headers=headers,params=report_args)
    if r.status_code not in [200]:
        raise Exception("Error: status: %s, %s" % (r.status_code, r.content))
    content = json.loads(r.content)
    logger.info('url: %r, meta: %r', useragreement_url, content[API_RESULT_META])
    if API_RESULT_DATA not in content:
        raise InformationError(
            'Content not recognized, no %r found' % API_RESULT_DATA)
    user_agreements = content[API_RESULT_DATA]

    def send_user_notification(user_agreement, date_to_expire):
        ''' 
        Send user notifications: 
        - make sure admin_email_only, test_only checks are met
        '''
        (msg_user_expiration_subject, msg_user_expiration_lines) = \
            EMAIL_MESSAGE_TEMPLATES['msg_user_expiration']
        msg_user_expiration_lines = msg_user_expiration_lines[:]
        valid_email = validate_email(ua.get(FIELD_USER_EMAIL, None))
        if not valid_email:
            return False, 'Invalid email address'
        
        msg_user_expiration_subject = msg_user_expiration_subject.format(
            screen_type=SCREEN_TYPES[user_agreement['type']])
        for i,line in enumerate(msg_user_expiration_lines):
            msg_user_expiration_lines[i] = line.format(
                screen_type=SCREEN_TYPES[user_agreement['type']],
                date_to_expire=date_to_expire,
                data_sharing_level=user_agreement['data_sharing_level'],
                contact_info=args.contact_info )
        
        emailer.send_email([valid_email], msg_user_expiration_subject, 
            msg_user_expiration_lines)
        
        return True, 'Success'

    def fill_parms(txt):
        ''' fill text format parameters for messages '''
        return txt.format(
            count=len(user_agreements), 
            screen_type=SCREEN_TYPES[args.screen_type],
            date_active_cutoff=date_active_cutoff,
            current_date=current_time.strftime("%Y-%m-%d"),
            days_to_expire=args.days_to_expire,
            days_ahead_to_notify=args.days_ahead_to_notify)

    if args.days_ahead_to_notify:
        
        if len(user_agreements) == 0:
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_notify_no_actions']
            msg_subject = fill_parms(msg_subject)
            msg_body_lines = [fill_parms(txt) for txt in msg_body_lines]
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines)
        
        else:
            # Notify the users and update the date_notified, send admin email
            
            date_to_expire = current_time + datetime.timedelta(
                days=args.days_ahead_to_notify)
            date_to_expire = date_to_expire.strftime("%Y-%m-%d")
            
            user_agreement_notify_success = defaultdict(list)
            user_agreement_notify_fail = defaultdict(list)
            for ua in user_agreements:
                result, msg = send_user_notification(ua, date_to_expire)
                if result is True:
                    user_agreement_notify_success[msg].append(ua)
                else:
                    user_agreement_notify_fail[msg].append(ua)
            # PATCH the date_notified 
            ua_updates = []
            for msg, uas in user_agreement_notify_success.items():
                logger.info('create updates for %r', msg)
                for ua in uas:
                    ua_updates.append({
                        'screensaver_user_id': ua['screensaver_user_id'],
                        FIELD_AGREEMENT_TYPE: ua[FIELD_AGREEMENT_TYPE],
                        'date_notified': current_time.strftime("%Y-%m-%d") 
                    })
            if ua_updates:

                logger.info('PATCH: user agreement updates: %r', ua_updates)
                
                headers = {
                    'Content-Type': 'application/json',
                    HEADER_APILOG_COMMENT_CLIENT: args.comment
                }
                logger.info('headers: %r', headers)
                _url = useragreement_url
                if args.test_only is True:
                    _url += '?test_only=true'
                r = django_requests.patch(
                    _url, session,
                    data = json.dumps(ua_updates),
                    headers = headers)
                if r.status_code not in [200]:
                    content = json.loads(r.content)
                    if args.test_only is True:
                        logger.info('"test_only" run: response: %r', content)
                    else:
                        raise Exception("Error: status: %s, %s" 
                            % (r.status_code, content))
                else:
                    content = json.loads(r.content)
                    logger.info('PATCH result: %r', content)
                    logger.info('content: %r', content.keys())
                    logger.info('meta: %r', content[API_RESULT_META])
                    
            # Send the Admin email
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_notification']
            msg_subject = fill_parms(msg_subject)
            msg_body_lines = [fill_parms(txt) for txt in msg_body_lines]
            html_msg_body_lines = ['%s<br>' % line for line in msg_body_lines]

            # Create the report URL
            msg_body_lines.append('User Agreement Report:')
            html_msg_body_lines.append('<hr>')
            report_url = '/'.join([
                settings.APP_PUBLIC_DATA.site_url,'#useragreement'])
            report_url += '/includes/' + ','.join(report_args.get('includes',''))
            del report_args['includes']
            del report_args['date_notified__is_null']
            report_url += '/search/' + ';'.join([
                '%s=%s'%(k,v) for k,v in report_args.items()])
            msg_body_lines.append(report_url)
            html_msg_body_lines.append('<a href="%s">%s</a>' % (
                report_url, 
                '%s User Agreement Report' % SCREEN_TYPES[args.screen_type] ))
            
            # Create the report summary tables
            user_msg_fields = OrderedDict((
                (FIELD_USER_NAME, 'User'),
                ('screensaver_user_id', 'ID'),
                (FIELD_USER_EMAIL,'Email'),
                ('data_sharing_level','Data Sharing Level'),
                ('date_active','Date Active')
                ))
            html_msg_body_lines.append('<hr>')
            if user_agreement_notify_success:
                msg_body_lines.append('')
                msg_body_lines.append('Notification Summary:')
                html_msg_body_lines.append('')
                html_msg_body_lines.append('Notification Summary:')
                html_msg_body_lines.append('<hr>')
                html_msg_body_lines.append('')
                for msg, uas in user_agreement_notify_success.items():
                    msg_body_lines.append(msg)
                    html_msg_body_lines.append('<h4>%s:</h4>' % msg)
                    pt = create_prettytable(uas, user_msg_fields)
                    msg_body_lines.append(pt.get_string(padding_width=5))
                    html_msg_body_lines.append(pt.get_html_string(
                        attributes = {"class": "content-table"}))
                    
            if user_agreement_notify_fail:
                msg_body_lines.append('')
                html_msg_body_lines.append('')
                html_msg_body_lines.append('<hr>')
                for msg, uas in user_agreement_notify_fail.items():
                    msg_body_lines.append(msg)
                    html_msg_body_lines.append('<h4>%s:</h4>' % msg)
                    pt = create_prettytable(uas, user_msg_fields)
                    msg_body_lines.append(pt.get_string(padding_width=5))
                    html_msg_body_lines.append(pt.get_html_string(
                        attributes = {"class": "content-table"}))
            
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines,
                html_msg_body_lines=html_msg_body_lines )
        
    else:
        # Expire: "days_ahead_to_notify" is not set - expire the user accounts
        if len(user_agreements) == 0:
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_expire_no_actions']
            msg_subject = fill_parms(msg_subject)
            msg_body_lines = [fill_parms(txt) for txt in msg_body_lines]
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines)
        
        else:
            # Prepare PATCH - status 'expired'
            ua_updates = []
            for ua in user_agreements:
                ua_updates.append({
                    'screensaver_user_id': ua['screensaver_user_id'],
                    FIELD_AGREEMENT_TYPE: ua[FIELD_AGREEMENT_TYPE],
                    'status': 'expired'
                })

            logger.info('PATCH: user agreement updates: %r', ua_updates)
            
            headers = {
                'Content-Type': 'application/json',
                HEADER_APILOG_COMMENT_CLIENT: args.comment
            }
            
            _url = useragreement_url
            if args.test_only is True:
                _url += '?test_only=true'
            r = django_requests.patch(
                _url, session,
                data = json.dumps(ua_updates),
                headers = headers)
            if r.status_code not in [200]:
                content = json.loads(r.content)
                if args.test_only is True:
                    logger.info('"test_only" run: response: %r', content)
                else:
                    raise Exception("Error: status: %s, %s" 
                        % (r.status_code, content))
            else:
                content = json.loads(r.content)
                logger.info('PATCH result: %r', content)
                logger.info('content: %r', content.keys())
                logger.info('meta: %r', content[API_RESULT_META])
            
            # Send the Admin email
            (msg_subject, msg_body_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_admin_expiration']
            msg_subject = fill_parms(msg_subject)
            msg_body_lines = [fill_parms(txt) for txt in msg_body_lines]
            html_msg_body_lines = ['%s<br>' % line for line in msg_body_lines]
            
            # Create the report URL
            msg_body_lines.append('User Agreement Report:')
            html_msg_body_lines.append('<hr>')
            report_url = '/'.join([
                settings.APP_PUBLIC_DATA.site_url,'#useragreement'])
            report_url += '/includes/' + ','.join(report_args.get('includes',''))
            del report_args['includes']
            report_url += '/search/' + ';'.join([
                '%s=%s'%(k,v) for k,v in report_args.items()])
            msg_body_lines.append(report_url)
            html_msg_body_lines.append('<a href="%s">%s</a>' % (
                report_url, 
                '%s User Agreement Report' % SCREEN_TYPES[args.screen_type] ))
            
            # Create the report summary tables
            user_msg_fields = OrderedDict((
                (FIELD_USER_NAME, 'User'),
                ('screensaver_user_id', 'ID'),
                (FIELD_USER_EMAIL,'Email'),
                ('data_sharing_level','Data Sharing Level'),
                ('date_active','Date Active')
                ))
            html_msg_body_lines.append('<hr>')
            msg_body_lines.append('')
            pt = create_prettytable(user_agreements, user_msg_fields)
            msg_body_lines.append(pt.get_string(padding_width=5))
            html_msg_body_lines.append(
                pt.get_html_string(attributes = {"class": "content-table"}))

            emailer.send_email(admin_email_list, msg_subject, msg_body_lines,
                html_msg_body_lines=html_msg_body_lines )
                    

