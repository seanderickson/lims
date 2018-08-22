from __future__ import unicode_literals

from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText        
import logging
import re
import smtplib
import os.path


from django.conf import settings
import prettytable
import datetime
import base64
from logging.handlers import RotatingFileHandler

MAX_ARCHIVE_FILE_SIZE=(2**20*20) # 20 Mb
MAX_ARCHIVE_FILE_BACKUPS=20 # save up to 20 files in the rotating list

logger = logging.getLogger(__name__)

# If set, admin email only will be sent
TESTING_MODE = True

def create_prettytable(list_of_dicts, keys_and_titles ):

    _table = prettytable.PrettyTable(keys_and_titles.values())
    for _dict in list_of_dicts:
        _table.add_row([_dict[key] for key in keys_and_titles.keys()])    
    return _table

def read_email_template(f):
    '''
    read the email template:
    - first line is the subject
    - other lines are the message
    @return (subject, [message_lines])
    '''
    
    subject = None
    message_lines = []
    for i,line in enumerate(f):
        if i == 0:
            subject = line.strip()
        else:
            message_lines.append(line.strip())
    return (subject,message_lines)

def validate_email(email_addr):
    if not email_addr:
        return False
    email_addr = email_addr.strip()
    EMAIL_REGEX = re.compile(r'[^@\s]+@[^@\s]+\.[^@\s]+')
    
    if not EMAIL_REGEX.match(email_addr):
        return False
    if 'has.missing.email' in email_addr:
        return False
    return email_addr

def fill_email_wrapper(
    email_wrapper_as_text, email_title, email_content, email_footer,
    preheader_text=None):
    ''' 
    Fill the HTML email wrapper content holes:
    <![CDATA[preheader_text]]>
    <![CDATA[email_title]]>
    <![CDATA[email_content]]>
    <![CDATA[email_footer]]>
    '''
    if '<![CDATA[email_title]]>' not in email_wrapper_as_text:
        raise Exception('no title...')
    logger.info('email_title: %r', email_title)
    if preheader_text is None:
        preheader_text = email_title
    new_text = email_wrapper_as_text.replace(
        '<![CDATA[email_title]]>', email_title)
    new_text = new_text.replace(
        '<![CDATA[preheader_text]]>', preheader_text)
    new_text = new_text.replace(
        '<![CDATA[email_content]]>', email_content)
    new_text = new_text.replace(
        '<![CDATA[email_footer]]>', email_footer)
    return new_text

class Emailer(object):
    '''
    Sends admin emails; with stateful flags for testing modes
    '''
    
    def __init__(self, admin_email_address, admin_from_email_address=None,
        no_email=False, admin_email_only=False, test_only=False,
        html_email_wrapper=None, email_log_filename=None):
        '''
        @param admin_email_address - email header "from"
        @param no_email - print all email to stdout only
        @param admin_email_only - send all email to the admin_email_address only
        @param test_only - prepend "Test" to messages; admin email only
        @param html_email_wrapper - HTML template for email
        @param email_log_filename - (optional) - print email to log file given
        '''
        self.admin_email_address = admin_email_address
        self.admin_from_email_address = admin_email_address
        if admin_from_email_address is not None:
            self.admin_from_email_address = admin_from_email_address
         
        self.no_email = no_email
        self.admin_email_only = admin_email_only
        self.test_only = test_only
        if test_only is True:
            self.admin_email_only = True
        self.html_email_wrapper = html_email_wrapper
        if email_log_filename:
            if not os.path.exists(
                os.path.dirname(email_log_filename)):
                raise Exception(
                    'specified email_log_filename : %r path does not exist %r'
                    % (email_log_filename, os.path.dirname(email_log_filename)))
            else:
                # Use the logging RotatingFileHandler as a convenient wrapper
                self.email_log_filename = email_log_filename
                self.archiving_logger = logging.getLogger('email_log_filename')
                self.archiving_logger.propagate = 0
                self.archiving_logger.setLevel(logging.DEBUG)
                self.archiving_handler = RotatingFileHandler(
                    email_log_filename, maxBytes=MAX_ARCHIVE_FILE_SIZE, 
                    backupCount=MAX_ARCHIVE_FILE_BACKUPS)
                self.archiving_logger.addHandler(self.archiving_handler)                
        
    def send_email(self,
        address_to_list, msg_subject, msg_body_lines,
        html_msg_body_lines=None, msg_footer=None):
        ''' 
        send an email from the admin to the address list
        - if self.html_email_wrapper has been set also send as html, look for
        html_msg_body_lines
        '''

        if self.test_only is True:
            msg_subject = '(Test Mode) ' + msg_subject
            
        if self.admin_email_only is True or TESTING_MODE is True:
            if address_to_list != [self.admin_email_address]:
                msg_body_lines.insert(0,'')
                msg_body_lines.insert(0,
                    '(Admin email only) Message for user: %s' 
                        % ', '.join(address_to_list))
                if html_msg_body_lines:
                    html_msg_body_lines.insert(0,'<hr>')
                    html_msg_body_lines.insert(0,
                        '(Admin email only) - Message for user: %s' 
                            % ', '.join(address_to_list))
                address_to_list = [self.admin_email_address]

        if not msg_footer:
            msg_footer = ' | '.join([
                '<a href="%s">%s</a>' % (
                    settings.APP_PUBLIC_DATA.site_url,
                    settings.APP_PUBLIC_DATA.app_name ),
                '<a href="%s">%s</a>' % (
                    settings.APP_PUBLIC_DATA.facility_url,
                    settings.APP_PUBLIC_DATA.facility_name ),
                 '<a href="mailto:%s">%s</a>' % (
                     settings.APP_PUBLIC_DATA.contact_feedback_email,
                     settings.APP_PUBLIC_DATA.contact_feedback_name ),
                ])
        
        # Create message container
        msg = MIMEMultipart('alternative')
        msg['Subject'] = msg_subject
        msg['From'] = self.admin_from_email_address
        msg['To'] = ', '.join(address_to_list)

        # Record the MIME types of both parts - text/plain and text/html.
        # NOTE: ensure utf-8 encoding
        plaintext = '\n'.join(msg_body_lines).encode('utf-8')
        part1 = MIMEText(plaintext, 'plain',  'utf-8')
        
        if self.html_email_wrapper:
            if not html_msg_body_lines:
                html_msg_body_lines = msg_body_lines[:]
            lines = ['<p>']
            for line in html_msg_body_lines:
                if line.strip():
                    lines.append(' %s' % line.strip())
                else:
                    lines.append('</p><p>')
            lines.append('</p>')
            email_content = ''.join(lines) #.encode('utf-8')
            html_email_body = fill_email_wrapper(
                self.html_email_wrapper, 
                msg_subject, email_content, msg_footer)
            html_email_body = html_email_body.encode('utf-8')
            part2 = MIMEText(html_email_body, 'html', 'utf-8')
        else:
            email_content = '<br/>'.join(msg_body_lines).encode('utf-8')
            part2 = MIMEText(email_content, 'html', 'utf-8')
        
        logger.debug('send message to: %r msg_subject: %r, msg_body_lines: %r', 
            address_to_list, msg_subject, msg_body_lines)
        # Attach parts into message container.
        # According to RFC 2046, the last part of a multipart message, in this case
        # the HTML message, is best and preferred.
        msg.attach(part1)
        msg.attach(part2)        

        if self.no_email is not True:
            # Send the message via our own SMTP server, but don't include the
            # envelope header.
            s = smtplib.SMTP('localhost')
            send_result = s.sendmail(
                self.admin_from_email_address, address_to_list, msg.as_string())
            logger.info('sendmail result: %r', send_result)
            s.quit()        
        else:
            logger.info('"no_email" flag is set')
            print 'From:', self.admin_from_email_address
            print 'To:', address_to_list
            print 'Subject:', msg_subject
            print '\n'.join(msg_body_lines)
        if self.archiving_logger is not None:
            logger.info('log to file: %r', self.email_log_filename)
            timestamp = str(datetime.datetime.now())
            self.archiving_logger.debug('\n\n=== begin: %s ===\n' % timestamp)
            self.archiving_logger.debug(msg.as_string())
            self.archiving_logger.debug('\n\n=== end: %s ===\n' % timestamp)

