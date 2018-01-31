from __future__ import unicode_literals

from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText        
import logging
import re
import smtplib

from django.conf import settings
import prettytable


logger = logging.getLogger(__name__)


def create_prettytable(list_of_dicts, keys_and_titles ):

    _table = prettytable.PrettyTable(keys_and_titles.values())
    for _dict in list_of_dicts:
        _table.add_row([_dict[key] for key in keys_and_titles.keys()])    
    return _table

def read_email_template(f):
    '''
    read the email template:
    - first line is the subject
    - other lines are the messagex
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
    if preheader_text is None:
        preheader_text = email_title
    new_text = email_wrapper_as_text.replace(
        '<![CDATA[email_title]]>', email_title)
    new_text = email_wrapper_as_text.replace(
        '<![CDATA[preheader_text]]>', preheader_text)
    new_text = new_text.replace(
        '<![CDATA[email_content]]>', email_content)
    new_text = new_text.replace(
        '<![CDATA[email_footer]]>', email_footer)
    return new_text

class Emailer(object):
    '''
    Sends admin emails
    '''
    
    def __init__(self, admin_email_address, 
        no_email=False, admin_email_only=False, test_only=False,
        html_email_wrapper=None ):
        '''
        @param admin_email_address - email header "from"
        @param no_email - print all email to stdout only
        @param admin_email_only - send all email to the admin_email_address only
        @param test_only - prepend "Test" to messages; admin email only
        @param html_email_wrapper - HTML template for email
        '''
        
        self.admin_email_address = admin_email_address
        self.no_email = no_email
        self.admin_email_only = admin_email_only
        self.test_only = test_only
        if test_only is True:
            self.admin_email_only = True
        self.html_email_wrapper = html_email_wrapper
        
        
    def send_email(self,
        address_to_list, msg_subject, msg_body_lines,
        html_msg_body_lines=None, msg_footer=None):

        if self.no_email is True:
            logger.info('"no_email" flag is set')
            
            print '"no_email" mode: email to be sent:'
            print 'From:', self.admin_email_address
            print 'To:', address_to_list
            print 'Subject:', msg_subject
            print '\n'.join(msg_body_lines)
            return
        else:
            if not html_msg_body_lines:
                html_msg_body_lines = ['%s<br>' % line for line in msg_body_lines]
                
            if self.test_only is True:
                msg_subject = '(Test Mode) ' + msg_subject
                
            if self.admin_email_only is True or TESTING_MODE is True:
                if address_to_list != [self.admin_email_address]:
                    msg_body_lines.insert(0,'')
                    msg_body_lines.insert(0,
                        '(Admin email only) Message for user: %r' % address_to_list)
                    html_msg_body_lines.insert(0,'<hr>')
                    html_msg_body_lines.insert(0,
                        '(Admin email only) - Message for user: %r' % address_to_list)
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
            msg['From'] = self.admin_email_address
            msg['To'] = ', '.join(address_to_list)
    
            # Record the MIME types of both parts - text/plain and text/html.
            part1 = MIMEText('\n'.join(msg_body_lines), 'plain')
            
            if self.html_email_wrapper:
                html_email_body = fill_email_wrapper(
                    self.html_email_wrapper, 
                    msg_subject, ' '.join(html_msg_body_lines), msg_footer)
                part2 = MIMEText(html_email_body, 'html')
            else:
                part2 = MIMEText('<br/>'.join(msg_body_lines), 'html')
            
            logger.info('send message to: %r msg_subject: %r, msg_body_lines: %r', 
                address_to_list, msg_subject, msg_body_lines)
            # Attach parts into message container.
            # According to RFC 2046, the last part of a multipart message, in this case
            # the HTML message, is best and preferred.
            msg.attach(part1)
            msg.attach(part2)        
            # Send the message via our own SMTP server, but don't include the
            # envelope header.
            s = smtplib.SMTP('localhost')
            send_result = s.sendmail(
                self.admin_email_address, address_to_list, msg.as_string())
            logger.info('sendmail result: %r', send_result)
            s.quit()        

