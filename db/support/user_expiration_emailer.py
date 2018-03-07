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

from db.schema import get_href, get_title, get_vocab_title, \
    replace_vocabularies, replace_html_values, DB_API_URI, DATE_FORMAT
import db.schema
from reports import InformationError, HEADER_APILOG_COMMENT_CLIENT
from reports.api import API_RESULT_DATA, API_RESULT_META
from reports.utils import parse_credentials
from reports.utils.admin_emailer import Emailer, read_email_template, \
    create_prettytable, validate_email
from reports.utils.django_requests import get_resource_listing, get_resource, \
    get_resource_schema
import reports.utils.django_requests as django_requests


logger = logging.getLogger(__name__)

USER = db.schema.SCREENSAVER_USER
UA = db.schema.USER_AGREEMENT
VOCAB = db.schema.VOCAB

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
    '-ua_type', '--useragreement_type', required=True,
    choices=VOCAB.user_agreement.type.get_dict().values(),
    help='User Agreement for Screen Type (sm="Small Molecule", rnai="SiRNA"')
parser.add_argument(
    '-days_to_expire', '--days_to_expire', type=int, required=True,
    help='number of days for a user agreement to expire. If set without the'
    '"days_ahead_to_notify"; will expire the user agreements.')
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
    
    # Get the user agreement schema
    schema_url = '/'.join([base_url, DB_API_URI, UA.resource_name, 'schema'])
    ua_schema = get_resource_schema(schema_url, session, headers)
    ua_type=get_vocab_title(UA.TYPE, args.useragreement_type,ua_schema)
    
    # Begin processing user agreements
    # Query for the user agreements
    
    current_time = datetime.datetime.now()
    date_active_cutoff = \
        current_time + datetime.timedelta(days=-(args.days_to_expire))
    if args.days_ahead_to_notify is not None:
        date_active_cutoff = \
            date_active_cutoff + datetime.timedelta(
                days=args.days_ahead_to_notify)
    
    date_active_cutoff = date_active_cutoff.strftime(DATE_FORMAT)
    logger.info(
        'days_to_expire: %d, days_ahead_to_notify: %r, date_active_cutoff: %s', 
        args.days_to_expire, args.days_ahead_to_notify, date_active_cutoff)
    report_args = {
        '%s__in' % UA.STATUS: VOCAB.user_agreement.status.ACTIVE,
        UA.TYPE: args.useragreement_type,
        '%s__lte' % UA.DATE_ACTIVE: date_active_cutoff,
        'includes': [UA.USER_EMAIL],
        'limit': 0
        }

    if args.days_ahead_to_notify is not None:
        report_args['%s__is_null' % UA.DATE_NOTIFIED] = True
    
    useragreement_url = '/'.join([base_url,DB_API_URI, UA.resource_name])
    
    user_agreements,meta = \
        get_resource_listing(useragreement_url, session, headers, report_args)
    

    def fill_parms(txt):
        ''' fill text format parameters for messages '''
        return txt.format(
            count=len(user_agreements), 
            screen_type=ua_type,
            date_active_cutoff=date_active_cutoff,
            current_date=current_time.strftime(DATE_FORMAT),
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
            date_to_expire = date_to_expire.strftime(DATE_FORMAT)
            
            user_agreement_notify_success = defaultdict(list)
            user_agreement_notify_fail = defaultdict(list)
            
            args.sample_message_subject = None
            args.sample_message_lines = None
            args.sample_recipient = None
            (msg_user_expiration_subject, msg_user_expiration_lines) = \
                EMAIL_MESSAGE_TEMPLATES['msg_user_expiration']
            def send_user_notification(user_agreement, date_to_expire):
                ''' 
                Send user notifications: 
                - make sure admin_email_only, test_only checks are met
                '''
                text_msg_body_lines = msg_user_expiration_lines[:]
                valid_email = validate_email(ua.get(UA.USER_EMAIL, None))
                if not valid_email:
                    return False, 'Invalid email address'
                
                msg_subject = msg_user_expiration_subject.format(
                    screen_type=ua_type)
                for i,line in enumerate(text_msg_body_lines):
                    text_msg_body_lines[i] = line.format(
                        screen_type=ua_type,
                        date_to_expire=date_to_expire,
                        data_sharing_level=get_vocab_title(
                            UA.DATA_SHARING_LEVEL, 
                            user_agreement[UA.DATA_SHARING_LEVEL], ua_schema),
                        contact_info=args.contact_info )

                if args.sample_message_subject is None:
                    args.sample_message_subject = 'Subject: %s' % msg_subject
                    args.sample_message_lines = text_msg_body_lines[:]
                    args.sample_recipient = valid_email
                                    
                emailer.send_email([valid_email], msg_subject, 
                    text_msg_body_lines)
                
                return True, 'Success'
            
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
                        UA.SCREENSAVER_USER_ID: 
                            ua[UA.SCREENSAVER_USER_ID],
                        UA.TYPE: ua[UA.TYPE],
                        UA.DATE_NOTIFIED: current_time.strftime(DATE_FORMAT) 
                    })
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
            del report_args['%s__is_null' % UA.DATE_NOTIFIED]
            report_url += '/search/' + ';'.join([
                '%s=%s'%(k,v) for k,v in report_args.items()])
            msg_body_lines.append(report_url)
            html_msg_body_lines.append('<a href="%s">%s</a>' % (
                report_url, 
                '%s User Agreement Report' % ua_type ))
            
            # Create the report summary tables
            report_fields = [
                UA.USER_NAME, UA.SCREENSAVER_USER_ID,
                UA.USER_EMAIL,UA.DATA_SHARING_LEVEL,
                UA.DATE_ACTIVE]
            table_fields = OrderedDict(
                (key,get_title(key, ua_schema)) for key in report_fields)
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
                    pt = create_prettytable(
                        [replace_vocabularies(ua,ua_schema) for ua in uas], table_fields)
                    msg_body_lines.append(pt.get_string(padding_width=5))

                    html_pt = create_prettytable(
                        [replace_html_values(ua,ua_schema) for ua in uas], 
                        table_fields)
                    html_msg_body_lines.append(html_pt.get_html_string(
                        attributes = {"class": "content-table"})
                        .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
                    
            if user_agreement_notify_fail:
                msg_body_lines.append('')
                html_msg_body_lines.append('')
                html_msg_body_lines.append('<hr>')
                for msg, uas in user_agreement_notify_fail.items():
                    msg_body_lines.append(msg)
                    html_msg_body_lines.append('<h4>%s:</h4>' % msg)
                    pt = create_prettytable(
                        [replace_vocabularies(ua,ua_schema) for ua in uas], table_fields)
                    msg_body_lines.append(pt.get_string(padding_width=5))

                    html_pt = create_prettytable(
                        [replace_html_values(ua,ua_schema) for ua in uas], 
                        table_fields)
                    html_msg_body_lines.append(html_pt.get_html_string(
                        attributes = {"class": "content-table"})
                        .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))

            # Append the sample message
            if args.sample_message_subject:
                msg_body_lines.append('[example email]')
                msg_body_lines.append('To: %s' % args.sample_recipient)
                msg_body_lines.append(args.sample_message_subject)
                msg_body_lines.extend(args.sample_message_lines)
    
                html_msg_body_lines.append('<hr>[example email]<hr/>')
                html_msg_body_lines.append('To: %s' % args.sample_recipient)
                html_msg_body_lines.append('')
                html_msg_body_lines.append(args.sample_message_subject)
                html_msg_body_lines.append('')
                html_msg_body_lines.extend(args.sample_message_lines)
            
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
                    UA.SCREENSAVER_USER_ID: 
                        ua[UA.SCREENSAVER_USER_ID],
                    UA.TYPE: ua[UA.TYPE],
                    UA.STATUS: VOCAB.user_agreement.status.EXPIRED
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
                '%s User Agreement Report' % ua_type ))
            
            # Create the report summary tables
            report_fields = [
                UA.USER_NAME, UA.SCREENSAVER_USER_ID,
                UA.USER_EMAIL,UA.DATA_SHARING_LEVEL,
                UA.DATE_ACTIVE]
            table_fields = OrderedDict(
                (key,get_title(key, ua_schema)) for key in report_fields)
            html_msg_body_lines.append('<hr>')
            msg_body_lines.append('')
            pt = create_prettytable(
                [replace_vocabularies(ua,ua_schema) for ua in user_agreements], 
                table_fields)
            msg_body_lines.append(pt.get_string(padding_width=5))

            html_pt = create_prettytable(
                [replace_html_values(ua,ua_schema) for ua in user_agreements], 
                table_fields)
            html_msg_body_lines.append(html_pt.get_html_string(
                attributes = {"class": "content-table"})
                .replace('&lt;','<').replace('&gt;','>').replace('&amp;','&'))
            
            emailer.send_email(admin_email_list, msg_subject, msg_body_lines,
                html_msg_body_lines=html_msg_body_lines )
                    

