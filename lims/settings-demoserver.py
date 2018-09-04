# Django settings for the Screensaver Demo Server
# Copy this file to settings.py and change local values for your installation

from os import environ
import sys

try:
    from base_settings import *
except ImportError:
    print >>sys.stderr, '''Base Settings not defined.  Please configure a version of
    base_settings.py for this site.'''

print 'PROJECT_ROOT: ', PROJECT_ROOT

try:
    from app_data_iccbl import APP_PUBLIC_DATA
except ImportError:
    print >>sys.stderr, '''APP_PUBLIC_DATA not defined.  Please configure a version of
    app_data.py for this site.'''

# NOTE THAT DEBUG SHOULD NEVER BE True FOR A PUBLIC FACING INSTALLATION
DEBUG = False 

# If not True, then only staff may log in to the system
# see reports/auth.py
IS_PRODUCTION_READY = True

# NOTE: SSL may only be enforced on the production server
# NOTE: the migration app uses insecure HTTP to initialize, 
# TODO: enable these settings when in production
#if IS_PRODUCTION_READY is True:
#    SECURE_SSL_REDIRECT = True
#    SESSION_COOKIE_SECURE = True
#    CSRF_COOKIE_SECURE = True

ADMINS = (
    ('Admin', 'demoadmin@lims.com'),
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

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = [
    '127.0.0.1',
    'localhost',
    ]

TIME_ZONE = 'US/Eastern'

LANGUAGE_CODE = 'en-us'

SITE_ID = 1

AUTHENTICATION_BACKENDS = ('reports.auth.CustomAuthenticationBackend',)

STATIC_ROOT = os.path.join(PROJECT_ROOT, '..', '..', 'docroot', '_static')

STATIC_URL = '/_static/'

STATICFILES_DIRS = (
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = '##%%@@unique-key-demo-server11xx!'

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
WELL_STRUCTURE_IMAGE_DIR='/structure-images'

BACKGROUND_PROCESSING = True
APP_PUBLIC_DATA.BACKGROUND_PROCESSING = BACKGROUND_PROCESSING

BACKGROUND_PROCESSOR = {
    'post_data_directory':
        os.path.join(PROJECT_ROOT,'..','logs','background','post_data'),
    'job_output_directory':
        os.path.join(PROJECT_ROOT,'..','logs','background','job_output'),
    'credential_file':
        os.path.join(PROJECT_ROOT, 'db','static','demo_data','demosuperuser_credentials.txt'),
    'python_environ_script':
        os.path.join(PROJECT_ROOT, 'run.sh'),
    'background_process_script': 'reports.utils.background_client_util',
    # if "sbatch_settings" is set, use SLURM/sbatch to process
    'sbatch_settingsx': {
        'partition': 'short',
        'time': '00:01:00',
        'mem': '16G',
        'open-mode':'append',
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
            'handlers': ['logfile'],
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
