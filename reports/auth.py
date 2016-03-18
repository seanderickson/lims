from __future__ import unicode_literals
import logging
import re
from django.contrib.auth.models import User
from reports.hms.auth import authenticate
from django.core.exceptions import ValidationError

logger = logging.getLogger(__name__)

class CustomAuthenticationBackend():
    """
    Authenticate users - if auth_user has a valid password, use that, otherwise,
    use eCommons.
    """

    def authenticate(self, username=None, password=None):
        logger.info('find and authenticate the user: %s', username)

        # Test for superuser logging in as another user
        matchObject = re.match(r'(\w+)\:(\w+)', username)
        if(matchObject):
            superuser = matchObject.group(1)
            logged_in_as = matchObject.group(2)
            s_user = self._inner_authenticate(superuser, password)
            if s_user:
                if s_user.is_superuser:
                    try:
                        user = User.objects.get(username=logged_in_as)
                        if(user.is_active):
                            return user
                        else:
                            msg = 'user authenticated, but is not active: %r' % user
                            logger.warn(msg)
                            raise ValidationError(msg)
                    except User.DoesNotExist, e:
                        msg = 'no such user with the id: %r' % username
                        logger.warn(msg)
                        raise ValidationError(msg)
                else:
                    msg = 'user: %r does not have superuser privileges' % superuser
                    logger.warn(msg)
                    raise ValidationError(msg)
            return None
        else:
            return self._inner_authenticate(username, password)

    def _inner_authenticate(self, username=None, password=None):
        logger.info('innner_authenticate: %r', username)
        try:
            user = User.objects.get(username=username)
            if user.has_usable_password():
                if(user.check_password(password)):
                    return user
                else:
                    msg = 'user password authentication failed: %r' % username
                    logger.info(msg)
                    raise ValidationError(msg)
            logger.info("no password set, trying to authenticate with ecommons...")
            if(authenticate(username, password)):
                logger.info('user %r authenticated with the ecommons server', user)
                if(user.is_active):
                    return user
                else:
                    msg = 'user authenticated, but is not active: %r' % user
                    logger.warn(msg)
                    raise ValidationError(msg)
            else:
                msg = 'user not authenticated with the ecommons server: %r' % user
                logger.warn(msg)
                raise ValidationError(msg)
        except User.DoesNotExist, e:
            msg = 'no such user with the id: %r' % username
            logger.warn(msg)
            raise ValidationError(msg)

    def get_user(self, user_id):
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None