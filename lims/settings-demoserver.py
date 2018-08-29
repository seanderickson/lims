# Django settings for the Screensaver Demo Server
# Copy this file to settings.py and change local values for your installation

from os import environ

try:
    from base_settings import *
except ImportError:
    import sys
    print >>sys.stderr, '''Base Settings not defined.  Please configure a version of
    base_settings.py for this site.'''
    del sys
    
# NOTE: the parent settings file defines the PROJECT_ROOT
print 'PROJECT_ROOT: ', PROJECT_ROOT
    
# NOTE THAT DEBUG SHOULD NEVER BE True FOR A PUBLIC FACING INSTALLATION
# - If debugging is required on the Orchestra server, first disable non-HMS
# access through the docroot/.htaccess file.
# (LEAKS environment variables, i.e. database password)
DEBUG = True

# If not True, then only staff may log in to the system
# see reports/auth.py
IS_PRODUCTION_READY = True

# NOTE: SSL may only be enforced on the production server
# NOTE: the migration app uses insecure HTTP to initialize, 
# TODO: enable these settings when in production
# if IS_PRODUCTION_READY is True:
#     SECURE_SSL_REDIRECT = True
#     SESSION_COOKIE_SECURE = True
#     CSRF_COOKIE_SECURE = True

ADMINS = (
    ('Admin', 'admin@email.com'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql', 
        'NAME': 'demoscreensaver', 
        'USER': 'demoscreensaver',
        'PASSWORD': '',
        'HOST': 'localhost',
        'PORT': '',                      # Set to empty string for default.
    },
}

# Note that the SCREENSAVER_PGSQL variables can be found in the appropriate file at:
# /opt/apache/conf/auth/[server_address]
# to access these variables from the commmand line, see "/setenv_and_run.sh"  
_dbdefault = DATABASES['default']
if 'DEMOSCREENSAVER_PGSQL_SERVER' in environ:
    # explicit db configuration for lincs site in environment variables
    _dbdefault['NAME'] = environ['DEMOSCREENSAVER_PGSQL_DB']
    _dbdefault['HOST'] = environ['DEMOSCREENSAVER_PGSQL_SERVER']
    _dbdefault['USER'] = environ['DEMOSCREENSAVER_PGSQL_USER']
    _dbdefault['PASSWORD'] = environ['DEMOSCREENSAVER_PGSQL_PASSWORD']



# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
# NOTE that 'dev.screensaver2.med.harvard.edu' is an alias for
# 'dev.orchestraweb.med.harvard.edu'
ALLOWED_HOSTS = [
    '127.0.0.1',
    'localhost',
    'dev.orchestraweb.med.harvard.edu', 
    'demo.screensaver.med.harvard.edu'
]

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'US/Eastern'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

AUTHENTICATION_BACKENDS = ('reports.auth.CustomAuthenticationBackend',)

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = ''
# STATIC_ROOT = os.path.join(PROJECT_ROOT, '..', '..', 'docroot', '_static')

# URL prefix for static files.
# Example: "http://example.com/static/", "http://static.example.com/"
STATIC_URL = '/_static/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
    os.path.join(PROJECT_ROOT, "lims", "static"),
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = 'tell_no_one1_demo2##$%^^!xxx'

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'django_cache',
        'TIMEOUT': None,
        'OPTIONS': {
            'MAX_ENTRIES': 50000 
        },
    },
    'reports_cache': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'reports_cache',
        'TIMEOUT': None,
        'OPTIONS': {
            'MAX_ENTRIES': 100000,
        }
    },
    'resource_cache': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'resource_cache',
        'TIMEOUT': None,
        'OPTIONS': {
            'MAX_ENTRIES': 100000,
        }
    },
    'db_cache': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'db_cache',
        'TIMEOUT': None,
        'OPTIONS': {
            'MAX_ENTRIES': 100000,
        }
    },
    'screen_cache': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'screen_cache',
        'TIMEOUT': None,
        'OPTIONS': {
            'MAX_ENTRIES': 100000,
        }
    },
}

# if structure image cache directory is available.  see db.api for details.
WELL_STRUCTURE_IMAGE_DIR='/home/sde4/docs/work/images/structure-images'

BACKGROUND_PROCESSING = True
BACKGROUND_PROCESSOR = {
    'post_data_directory': 
        os.path.join(PROJECT_ROOT,'background','post_data'),
    'job_output_directory': 
        os.path.join(PROJECT_ROOT,'background','job_output'),
    'credential_file': 
        os.path.join(PROJECT_ROOT, 'demoadmin_credentials.txt'),
    'python_environ_script':
        os.path.join(PROJECT_ROOT, 'run_dev.sh'),

    'sbatch_setings': {
        'partition': 'short',
        'time': '00:02:00',
        'mem': '16G',
        },
    }

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
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        },
        'console':{
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },  
        'logfile': {
            'level':'DEBUG',
            'class':'logging.handlers.RotatingFileHandler',
            'filename': os.path.join(PROJECT_ROOT, 'logs', 'iccbl-lims.log'),
            'maxBytes': 5000000,
            'backupCount': 4,
            'formatter': 'simple',
        },
    },
    'loggers': {
        'django': {
            'handlers': ['mail_admins','console','logfile'],
            'level': 'WARN',
            'propagate': False,
        },
        'db': {  # set a default handler
            'handlers': ['console','logfile'],
            'propagate': True,
            'level': 'INFO',
        },        
        'lims': {  # set a default handler
            'handlers': ['console','logfile'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports': {  # set a default handler
            'handlers': ['console','logfile'],
            'propagate': True,
            'level': 'INFO',
        },        
        # suppress streaming errors (i.e. image not found)
        'reports.serialize.streaming_serializers': {  
            'handlers': ['console','logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
#         'db.views': {  
#             'handlers': ['logfile'],
#             'propagate': True,
#             'level': 'WARN',
#         },        
        'django.db': {  # for SQL
            'handlers': ['console','logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
        'utils': {  # 
            'handlers': ['console','logfile'],
            'propagate': True,
            'level': 'INFO',
        },
    }
}
