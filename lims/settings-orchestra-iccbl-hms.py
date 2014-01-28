# Django settings for lims project on the HMS Orchestra system
# Copy this file to settings.py and change local values for your installation

from os import environ

try:
    from base_settings import *
except ImportError:
    import sys
    print >>sys.stderr, '''Base Settings not defined.  Please configure a version of
    base_settings.py for this site.'''
    del sys
    
    
DEBUG = True
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    ('Admin', 'admin@email.com'),
)

MANAGERS = ADMINS


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', 
        'NAME': 'devscreensaver2', 
        'USER': 'devscreensaver2web',
        'PASSWORD': '',
        'HOST': 'dev.pgsql.orchestra:',
        'PORT': '',                      # Set to empty string for default.
    }
}

# Note that the SCREENSAVER_PGSQL variables can be found in the appropriate file at:
# /opt/apache/conf/auth/[server_address]
# to access these variables from the commmand line, see "/setenv_and_run.sh"  
_dbdefault = DATABASES['default']
if 'SCREENSAVER_PGSQL_SERVER' in environ:
    # explicit db configuration for lincs site in environment variables
    _dbdefault['NAME'] = environ['SCREENSAVER_PGSQL_DB']
    _dbdefault['HOST'] = environ['SCREENSAVER_PGSQL_SERVER']
    _dbdefault['USER'] = environ['SCREENSAVER_PGSQL_USER']
    _dbdefault['PASSWORD'] = environ['SCREENSAVER_PGSQL_PASSWORD']



# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = []

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
STATIC_ROOT = os.path.join(BASE_DIR, '..', '..', '..', 'docroot', '_static')

# URL prefix for static files.
# Example: "http://example.com/static/", "http://static.example.com/"
STATIC_URL = '/_static/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = 'x1gg450))screensaver2-lims2014001212045'

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s:%(lineno)d %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(asctime)s: %(pathname)s:%(lineno)d:%(levelname)s: %(message)s'
        },
    },
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'logfile': {
            'level':'DEBUG',
            'class':'logging.handlers.RotatingFileHandler',
#            'filename': os.path.join(PROJECT_ROOT, '..') +  "/logs/screensaver2.log",
            'filename': "/www/dev.screensaver2.med.harvard.edu/support/logs/screensaver2.log",
            'maxBytes': 5000000,
            'backupCount': 2,
            'formatter': 'simple',
        },
        'console':{
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },  
    },
    'loggers': {
        'django.request': {
            'handlers': ['logfile'],
            'level': 'ERROR',
            'propagate': False,
        },
        'db': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
        'lims': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },
        'db.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'lims.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'django': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
        'utils': {  # for SQL
            'handlers': ['logfile'],
            'propagate': True,
            'level': 'INFO',
        },        
        'tastypie': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
    }
}

