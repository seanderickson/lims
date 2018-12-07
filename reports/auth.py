from __future__ import unicode_literals
import logging
import re
from django.conf import settings
from django.contrib.auth.models import User

from reports import LoginFailedException
import reports.hms.auth


logger = logging.getLogger(__name__)

DEBUG_AUTHENTICATION = False

USER_PROXY_LOGIN_PATTERN = re.compile(r'(\w+)\:(\w+)')
USER_PROXY_ADMIN_GROUP = 'userView'

class CustomAuthenticationBackend():
    """
    Authenticate users - if auth_user has a valid password, use that, otherwise,
    use eCommons.
    """

    def authenticate(self, request, username=None, password=None):
        '''
        Authenticate the user given by the username and password.
        
        @return the user object.
        
        Note: Special superuser login as another user:
        if the syntax <some_superuser-username>:<username> is used, 
        then log in using the <superuser-username> and superuser password,
        but return the user for the subordinate <username> given.
        '''
        username = username.lower()
        logger.info('find and authenticate the user: %s', username)

        # Proxy Login 
        matchObject = USER_PROXY_LOGIN_PATTERN.match(username)
        if(matchObject):
            staff_user = matchObject.group(1)
            logged_in_as = matchObject.group(2)
            s_user = self._inner_authenticate(staff_user, password)
            if s_user:
                is_allowed = s_user.is_staff
                if not is_allowed:
                    try:
                        is_allowed = (
                            s_user.userprofile.usergroup_set
                                .filter(name=USER_PROXY_ADMIN_GROUP).exists())
                        if is_allowed:
                            logger.info('user: %r is in the group: %r',
                                staff_user,USER_PROXY_ADMIN_GROUP)
                    except Exception, e:
                        logger.exception(
                            'Unexpected error querying user groups: %r', e)
                if is_allowed:
                    try:
                        user = User.objects.get(username=logged_in_as)
                        logger.info('logged in staff user %r as %r',
                            staff_user, logged_in_as)
                        if user.is_superuser:
                            if not s_user.is_superuser:
                                raise PermissionDenied(
                                    '"login as" may not be used to access a superuser account')
                        if user.is_staff:
                            if not s_user.is_superuser:
                                raise PermissionDenied(
                                    '"login as" may only be used by super users to access another staff account')
                        return user
                    except User.DoesNotExist, e:
                        msg = 'no such user with the id: %r' % logged_in_as
                        logger.warn(msg)
                        raise LoginFailedException(msg)
                else:
                    msg = 'Log in as user: %r fails for user %r: '\
                        'Must be a superuser, or member of group: %r'
                    msg = msg % (logged_in_as, staff_user, USER_PROXY_ADMIN_GROUP)
                    logger.warn(msg)
                    raise LoginFailedException(msg)
            return None
        else:
            user = self._inner_authenticate(username, password)
            if settings.IS_PRODUCTION_READY is not True:
                if not any([user.is_staff,user.is_superuser]):
                    msg = ('login not allowed for non-staff users when '
                        '"settings.IS_PRODUCTION_READY" is not set; user: %r'
                        % user)
                    logger.warn(msg)
                    raise LoginFailedException(msg)
            logger.info('User logged in: %r', username)
            return user

    def _inner_authenticate(self, username, password):
        username = username.lower()
        if DEBUG_AUTHENTICATION:
            logger.info('inner_authenticate: %r', username)
        try:
            user = User.objects.get(username=username)
            if user.has_usable_password():
                if(user.check_password(password)):
                    logger.info('logged in user %r using password',username)
                    return user
                else:
                    msg = 'User password authentication failed: "%s"' % username
                    logger.info(msg)
                    raise LoginFailedException(msg)
            elif reports.hms.auth.authenticate(username, password):
                logger.info(
                    'eCommons authentication succeeded: %r', username)
                if(user.is_active):
                    return user
                else:
                    msg = 'User authenticated, but is not active: "%s"' % username
                    logger.warn(msg)
                    raise LoginFailedException(msg)
            else:
                # Note: reports.hms.auth.authenticate raises an Exception instead
                msg = 'User not authenticated with the ecommons server: %s' % username
                logger.warn(msg)
                raise LoginFailedException(msg)
        except User.DoesNotExist, e:
            msg = 'No user found for: "%s"' % username
            logger.warn(msg)
            raise LoginFailedException(msg)
        except Exception, e:
            logger.warn('auth ex: %r', e)
            raise LoginFailedException(e)

    def get_user(self, user_id):
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None