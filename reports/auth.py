from __future__ import unicode_literals
import logging
import re
from django.conf import settings
from django.contrib.auth.models import User
from reports.hms.auth import authenticate
from django.core.exceptions import PermissionDenied

logger = logging.getLogger(__name__)

IS_PRODUCTION_READY = settings.IS_PRODUCTION_READY or False

USER_PROXY_LOGIN_PATTERN = re.compile(r'(\w+)\:(\w+)')
USER_PROXY_ADMIN_GROUP = 'screensaverUserAdmin'

class CustomAuthenticationBackend():
    """
    Authenticate users - if auth_user has a valid password, use that, otherwise,
    use eCommons.
    """

    def authenticate(self, username=None, password=None):
        logger.info('find and authenticate the user: %s', username)

        # Proxy Login 
        matchObject = USER_PROXY_LOGIN_PATTERN.match(username)
        if(matchObject):
            superuser = matchObject.group(1)
            logged_in_as = matchObject.group(2)
            s_user = self._inner_authenticate(superuser, password)
            if s_user:
                is_allowed = s_user.is_superuser
                if not is_allowed:
                    try:
                        is_allowed = (
                            s_user.userprofile.usergroup_set
                                .filter(name=USER_PROXY_ADMIN_GROUP).exists())
                        if is_allowed:
                            logger.info('user: %r is in the group: %r',
                                superuser,USER_PROXY_ADMIN_GROUP)
                    except Exception:
                        logger.info('user %r is not in the %r group',
                            superuser,USER_PROXY_ADMIN_GROUP)
                if is_allowed:
                    try:
                        user = User.objects.get(username=logged_in_as)
                        logger.info('logged in super user %r as %r',
                            superuser, logged_in_as)
                        return user
                    except User.DoesNotExist, e:
                        msg = 'no such user with the id: %r' % username
                        logger.warn(msg)
                        raise PermissionDenied(msg)
                else:
                    msg = (
                        'logging in as another user requires superuser privileges'
                        ', user: %r') % superuser
                    logger.warn(msg)
                    raise PermissionDenied(msg)
            return None
        else:
            user = self._inner_authenticate(username, password)
            if IS_PRODUCTION_READY is not True:
                if user.is_staff is not True:
                    logger.warn('login not allowed for non-staff users when '
                        '"settings.IS_PRODUCTION_READY" is set True; user: %r', 
                        user)
                    raise PermissionDenied
            return user

    def _inner_authenticate(self, username=None, password=None):
        if username is None:
            raise PermissionDenied('username not set')
        username = username.lower()
        logger.debug('inner_authenticate: %r', username)
        try:
            user = User.objects.get(username=username)
            if user.has_usable_password():
                if(user.check_password(password)):
                    logger.info('logged in user %r using password',username)
                    return user
                else:
                    msg = 'user password authentication failed: %r' % username
                    logger.info(msg)
                    raise PermissionDenied(msg)
            if(authenticate(username, password)):
                logger.info('user %r authenticated with the ecommons server', user)
                if(user.is_active):
                    return user
                else:
                    msg = 'user authenticated, but is not active: %r' % user
                    logger.warn(msg)
                    raise PermissionDenied(msg)
            else:
                msg = 'user not authenticated with the ecommons server: %r' % user
                logger.warn(msg)
                raise PermissionDenied(msg)
        except User.DoesNotExist, e:
            msg = 'no such user with the id: %r' % username
            logger.warn(msg)
            raise PermissionDenied(msg)

    def get_user(self, user_id):
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None