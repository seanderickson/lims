from __future__ import unicode_literals 

# Django settings for lims project

import sys
import os

try:
    from app_data import APP_PUBLIC_DATA
except ImportError:
    print >>sys.stderr, '''app_data.py file not found.'''

PROJECT_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),'..'))

DEBUG = True

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


# Hosts/domain names that are valid for this site.
# INSTALL TODO: required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = []

# NOTE: Django internationalization is not being used
USE_I18N = False
# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
# NOTE: Django internationalization is not being used; so this has no effect.
LANGUAGE_CODE = 'en-us'
# Django format localization is not being used
USE_L10N = False


# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# INSTALL TODO: set to the local time zone for the installation
TIME_ZONE = 'UTC'
# Use timezone aware datetime information when handling datetimes in Django:
# NOTE: postgresql (stores) values internally as UTC, and the connection 
# converts using the TIME_ZONE setting .
USE_TZ = True

# NOTE: /accounts/login is the default
LOGIN_URL = '/accounts/login/'
# Default if "next" is not given as a request param
LOGIN_REDIRECT_URL='/lims/'
LOGOUT_REDIRECT_URL='/lims/'

# Timeout, in seconds; will cause user logout: 8 hours
# NOTE: this will not log the browser out until a request is made.
SESSION_COOKIE_AGE = 60*60*8
SESSION_EXPIRE_AT_BROWSER_CLOSE = True
# NOTE: SSL may only be enforced on the production server
# SECURE_SSL_REDIRECT = True
# SESSION_COOKIE_SECURE = True
# CSRF_COOKIE_SECURE = True

STATIC_ROOT = ''
STATIC_URL = '/_static/'
STATICFILES_DIRS = ()
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)
# Make this unique, and don't share it with anybody.
SECRET_KEY = 'tell_no_1'
SITE_ID = 1

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.contrib.messages.context_processors.messages',
                "lims.webpack_bundle_hash_name_processor.bundle_context_processor",    
            ],
        },
    },
]


MIDDLEWARE = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    # Uncomment the next line for simple clickjacking protection:
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'lims.urls'

WSGI_APPLICATION = 'lims.wsgi.application'

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.admin',
    'django.contrib.admindocs',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
#     'aldjemy',
    'reports',
    'lims',
    'db',
)

# Turn off migrations during testing (just make the database from models.py)
SOUTH_TESTS_MIGRATE = False

# ICCBL-Setting: Directory for temp files created on download
TEMP_FILE_DIR='/tmp'

# Base path for profiling
PROFILE_LOG_BASE='/tmp'

# ICCBL-Setting: structure image cache directory if available.  
# @see db.views for details
WELL_STRUCTURE_IMAGE_DIR=''

# ICCBL-Setting: Maximum rows to cache in the database table "well_query_index"
# @see db.api.ScreenResultResource
MAX_WELL_INDEXES_TO_CACHE=3e+08

# ICCBL-Setting: Maximum rows to cache per query for cached_resultproxy:
# @see reports.sqlalchemy_resource
MAX_ROWS_FOR_CACHE_RESULTPROXY=1e4

# ICCBL-Setting: Minimum wells for insertion into the well_query_index before 
# clearing older indexes; for performance tuning on screen result / well queries.
# @see db.api.ScreenResultResource
MIN_WELLS_TO_CLEAR_INDEXES = 3e5

# If not True, then only staff may log in to the system
# see reports/auth.py
IS_PRODUCTION_READY = False

# NOTE: SSL may only be enforced on the production server
# NOTE: do not use with the migration app, as it uses insecure HTTP to initialize, 
# TODO: enable these settings when in production
if IS_PRODUCTION_READY is True:
    SECURE_SSL_REDIRECT = True
    SESSION_COOKIE_SECURE = True
    CSRF_COOKIE_SECURE = True

# ICCBL-Setting: For use when authenticating
# @see reports.api_base
BASIC_AUTH_REALM='screensaver'

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

