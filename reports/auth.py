from __future__ import unicode_literals
import logging
import re
from django.contrib.auth.models import User
from reports.hms.auth import authenticate

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
                            logger.warn(str(('user authenticated, but is not active',user)))
                    except User.DoesNotExist, e:
                        logger.error(str(('no such user with the id', username)))
                    except Exception, e:
                        logger.error(str(('failed to authenticate', username, e)))
                else:
                    logger.warn(str(('user:', superuser, 'does not have superuser privileges')))
            return None
        else:
            return self._inner_authenticate(username, password)

    def _inner_authenticate(self, username=None, password=None):
        
        try:
            user = User.objects.get(username=username)
            if user.has_usable_password():
                if(user.check_password(password)):
                    return user
                else:
                    logger.info('user password authentication failed')
                    return None
            logger.info("no password set, trying to authenticate with ecommons...")
            if(authenticate(username, password)):
                logger.info(str(('user authenticated with the ecommons server', user)))
                if(user.is_active):
                    return user
                else:
                    logger.warn(str(('user authenticated, but is not active',user)))
            else:
                logger.warn(str(('user not authenticated with the ecommons server', user)))
        except User.DoesNotExist, e:
            logger.error(str(('no such user with the id', username)))
        except Exception, e:
            logger.error(str(('failed to authenticate', username, e)))

    def get_user(self, user_id):
        logger.debug(str(('get_user',user_id)))
        try:
            return User.objects.get(pk=user_id)
        except User.DoesNotExist:
            return None