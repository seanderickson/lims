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

# NOTE: the parent settings file defines the PROJECT_ROOT
print 'PROJECT_ROOT: ', PROJECT_ROOT
    
try:
    from app_data_iccbl import APP_PUBLIC_DATA
    APP_PUBLIC_DATA.SESSION_COOKIE_AGE = SESSION_COOKIE_AGE
except ImportError:
    print >>sys.stderr, '''APP_PUBLIC_DATA not defined.  Please configure a version of
    app_data.py for this site.'''
    
# NOTE THAT DEBUG SHOULD NEVER BE True FOR A PUBLIC FACING INSTALLATION
# - If debugging is required on the Orchestra server, first disable non-HMS
# access through the docroot/.htaccess file.
# (LEAKS environment variables, i.e. database password)
DEBUG = False
TEMPLATE_DEBUG = DEBUG

# If not True, then only staff may log in to the system
# see reports/auth.py
IS_PRODUCTION_READY = False

# NOTE: SSL may only be enforced on the production server
# NOTE: the migration app uses insecure HTTP to initialize, 
# TODO: enable these settings when in production
if IS_PRODUCTION_READY is True:
    SECURE_SSL_REDIRECT = True
    SESSION_COOKIE_SECURE = True
    CSRF_COOKIE_SECURE = True

ADMINS = (
    ('Admin', 'admin@email.com'),
)

MANAGERS = ADMINS


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql', 
        'NAME': 'devscreensaver2', 
        'USER': 'devscreensaver2web',
        'PASSWORD': '',
        'HOST': 'dev.pgsql96.orchestra',
        'PORT': '',                      # Set to empty string for default.
    },
}

# Note that the SCREENSAVER_PGSQL variables can be found in the appropriate file at:
# /opt/apache/conf/auth/[server_address]
# to access these variables from the commmand line, see "/setenv_and_run.sh"  
_dbdefault = DATABASES['default']
if 'SCREENSAVER_PGSQL_SERVER' in environ:
    # explicit db configuration for lincs site in environment variables
    _dbdefault['NAME'] = environ['SCREENSAVER_PGSQL_DB']
    _dbdefault['HOST'] = environ['SCREENSAVER_PGSQL96_SERVER']
    _dbdefault['USER'] = environ['SCREENSAVER_PGSQL_USER']
    _dbdefault['PASSWORD'] = environ['SCREENSAVER_PGSQL_PASSWORD']



# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
# NOTE that 'dev.screensaver2.med.harvard.edu' is an alias for
# 'dev.orchestraweb.med.harvard.edu'
ALLOWED_HOSTS = [
    '127.0.0.1',
    'localhost',
    'dev.orchestraweb.med.harvard.edu', 
    'dev.screensaver2.med.harvard.edu']

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'US/Eastern'
APP_PUBLIC_DATA.TIME_ZONE = TIME_ZONE

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

AUTHENTICATION_BACKENDS = ('reports.auth.CustomAuthenticationBackend',)

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = os.path.join(PROJECT_ROOT, '..', '..', 'docroot', '_static')

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
)

# Make this unique, and don't share it with anybody.
SECRET_KEY = 'tell_no_one1_xxxw##!!!xsls%%#)*@'

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
WELL_STRUCTURE_IMAGE_DIR='/groups/screensaver/image_directory/structure-images'

BACKGROUND_PROCESSING = True
APP_PUBLIC_DATA.BACKGROUND_PROCESSING = BACKGROUND_PROCESSING

# NOTE PROJECT_ROOT is abs path, so remove the first slash
drive, path_and_file = os.path.splitdrive(PROJECT_ROOT)
def get_path_parts(path):
    folders = []
    while 1:
        path, folder = os.path.split(path)
        if folder != "":
            folders.append(folder)
        else:
            if path != "":
                folders.append(path)
            break
    folders.reverse()
    if folders[0] == os.path.sep:
        folders = folders[1:]
    print 'path', path, 'folders', folders
    return folders
def fix_path_for_o2(path):
    parts = get_path_parts(path)
    if parts[0] != 'n':
        parts.insert(0,'n')
    return os.path.join(os.path.sep,*parts) 
#O2_PROJECT_ROOT=os.path.join(os.path.sep,'n',*get_path_parts(PROJECT_ROOT))
O2_PROJECT_ROOT=fix_path_for_o2(PROJECT_ROOT)
print 'new O2 project root', O2_PROJECT_ROOT
# NOTE: if running on orchestra, should use PROJECT_ROOT here
# NOTE: 20180410 - /n/www, /n/groups *are* currently mounted on orchestra
BACKGROUND_PROCESSOR = {
    'post_data_directory': 
        os.path.join(O2_PROJECT_ROOT,'..','logs','background','post_data'),
    'job_output_directory': 
        os.path.join(O2_PROJECT_ROOT,'..','logs','background','job_output'),
    'credential_file': 
        os.path.join(O2_PROJECT_ROOT, '..','production_data','sde_credentials.1.txt'),
    'python_environ_script':
        os.path.join(O2_PROJECT_ROOT, 'run_prod.webconf02.sh'),
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
        'logfile': {
            'level':'DEBUG',
            'class':'logging.handlers.RotatingFileHandler',
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
    }
}

