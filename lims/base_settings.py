# Django settings for lims project

import sys
import os
import django.template

try:
    from app_data import APP_PUBLIC_DATA
except ImportError:
    print >>sys.stderr, '''app_data.py file not found.'''

PROJECT_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),'..'))

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    ('Site Admin', 'site_admin@email.com'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'sqlite3', 
        'NAME': os.path.join(PROJECT_ROOT, 'project.db'),       
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    }
}

# Hosts/domain names that are valid for this site; required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = []

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# FIXME: set to the local time zone for the installation
TIME_ZONE = 'UTC'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/var/www/example.com/media/"
MEDIA_ROOT = ''

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://example.com/media/", "http://media.example.com/"
MEDIA_URL = ''

# /accounts/login is the default
LOGIN_URL = '/accounts/login/'

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself; store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = ''

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

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

TEMPLATE_CONTEXT_PROCESSORS = (
    "django.contrib.auth.context_processors.auth",
    "django.core.context_processors.debug",
    "django.core.context_processors.i18n",
    "django.core.context_processors.media",
    "django.core.context_processors.static",
    "django.contrib.messages.context_processors.messages",
    "django.core.context_processors.request",
    "lims.webpack_bundle_hash_name_processor.bundle_context_processor",    
#     "lims.context_processors.login_url_with_redirect",
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    # Uncomment the next line for simple clickjacking protection:
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'lims.urls'

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'lims.wsgi.application'

SESSION_EXPIRE_AT_BROWSER_CLOSE = True

TEMPLATE_DIRS = [
    os.path.join(PROJECT_ROOT, 'lims','templates'),
    os.path.join(PROJECT_ROOT, 'reports','templates'),
    ]

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.admin',
    'django.contrib.admindocs',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'tastypie', # manual says this is "not necessary, but useful"
    'reports',
    'lims',
    # NOTE: initial reports migration: may require removal of "db" to solve
    # circular dependency
    'db',
)

# for tastypie: will evaluate resource URIs the same with or without the trailing slash
APPEND_SLASH=True
TASTYPIE_ALLOW_MISSING_SLASH=True

# turn off migrations during testing (just make the database from models.py)
SOUTH_TESTS_MIGRATE = False

# Default if "next" is not given as a request param
LOGIN_REDIRECT_URL='/lims'

# directory for temp files created on download
TEMP_FILE_DIR='/tmp'

# base path for profiling
PROFILE_LOG_BASE='/tmp'

# if structure image cache directory is available.  see db.views for details.
WELL_STRUCTURE_IMAGE_DIR=''

# maximum rows to cache in the database table well_query_index
# see db/api.ScreenResultResource
MAX_WELL_INDEXES_TO_CACHE=3e+08

# maximum rows to cache for cached_resultproxy:
# see reports/sqlalchemy_resource
# (max rows per query by hash)
MAX_ROWS_FOR_CACHE_RESULTPROXY=1e3

# minimum wells for insertion into the well_query_index before dropping indexes
# - for performance tuning on screen result / well queries
MIN_WELLS_TO_CLEAR_INDEXES = 3e5

# If not True, then only staff may log in to the system
# see reports/auth.py
IS_PRODUCTION_READY = False

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'unique-snowflake',        
    },
    'reports_cache': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'reports_cache'
    },
    'resource_cache': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'resource_cache'
    },
    'db_cache': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'db_cache'
    },
    'screen_cache': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
        'LOCATION': 'screen_cache'
    },
}

# set SQLALCHEMY_POOL_CLASS=sqlalchemy.pool.NullPool for testing
# environments, so that the test database can be destroyed
# import sqlalchemy.pool
# SQLALCHEMY_POOL_CLASS = sqlalchemy.pool.NullPool

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
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
        }
    },
    'loggers': {
        'django': {
            'handlers': ['mail_admins'],
            'level': 'WARN',
            'propagate': True,
        },
    }
}

