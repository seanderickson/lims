# Django settings for lims project

try:
    from base_settings import *
except ImportError:
    import sys
    print >>sys.stderr, '''Base Settings not defined.  Please configure a version of
    base_settings.py for this site.'''
    del sys
    
import os.path

# NOTE: the parent settings file defines the PROJECT_ROOT
PROJECT_ROOT = '.'
print 'PROJECT_ROOT: ', PROJECT_ROOT

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    ('Sean D Erickson', 'sean_erickson@hms.harvard.edu'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql', 
        'NAME': 'postgres',
        # The following settings are not used with sqlite3:
        'USER': 'postgres',
        'PASSWORD': '',
        'HOST': '',                      
        'PORT': '',                      
    }
}

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['localhost']

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'US/Eastern'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = ''

# URL prefix for static files.
# Example: "http://example.com/static/", "http://static.example.com/"
STATIC_URL = '/static/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
    os.path.join(PROJECT_ROOT, 'lims', 'static'),
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

#AUTHENTICATION_BACKENDS = ('reports.auth.CustomAuthenticationBackend',)

# Make this unique, and don't share it with anybody.
SECRET_KEY = 'x-=2g@c)#_z_1xxn0fpxj+y)v%n7&e#xpm+coo0(5^*@l___%8i*0travis-only'


BACKGROUND_PROCESSING = False

# set SQLALCHEMY_POOL_CLASS=sqlalchemy.pool.NullPool for testing
# environments, so that the test database can be destroyed
import sqlalchemy.pool
SQLALCHEMY_POOL_CLASS = sqlalchemy.pool.NullPool

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s:%(lineno)d %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(asctime)s: %(name)s:%(lineno)d: %(message)s'
        },
        'simple1': {
            'format': '%(levelname)s %(msecs)s: %(name)s:%(lineno)d: %(message)s'
        },
    },
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'console':{
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },  
    },
    'loggers': {
        'django.request': {
            'handlers': ['console'],
            'level': 'WARN',
            'propagate': False,
        },
        'db': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'WARN',
        },        
        'db.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'lims': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'WARN',
        },               
        'reports': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'WARN',
        },        
        'django.db': {  # for SQL
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'utils': {  # 
            'handlers': ['console'],
            'propagate': True,
            'level': 'WARN',
        },
    }
}

TEST_RUNNER='reports.tests.IccblTestRunner'

# disable migrations while testing
# see http://stackoverflow.com/questions/25161425/disable-migrations-when-running-unit-tests-in-django-1-7
class DisableMigrations(object):

    def __contains__(self, item):
        return True

    def __getitem__(self, item):
        return None

if 'test' in sys.argv[1:] or 'travis' in sys.argv[1:]:
    MIGRATION_MODULES = DisableMigrations()