# Django settings for 1km project
import sys

try:
    from settings import *
except ImportError:
    print >>sys.stderr, '''Settings not defined.  Please configure a version of
    settings.py for this site.'''
    
import os.path


print 'PROJECT_ROOT: ', PROJECT_ROOT, ', ' , os.path.join(PROJECT_ROOT, '..')

# set SQLALCHEMY_POOL_CLASS=sqlalchemy.pool.NullPool for testing
# environments, so that the test database can be destroyed
import sqlalchemy.pool
SQLALCHEMY_POOL_CLASS = sqlalchemy.pool.NullPool

# FIXME: sqllite3 db does not work - errors on "DISTINCT ON" clause
# DATABASES['default'] = {'ENGINE': 'django.db.backends.sqlite3',
#                         'NAME': ':memory'}

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s:%(lineno)d %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(asctime)s: %(name)s:%(funcName)s:%(lineno)d: %(message)s'
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
            'filename': os.path.join(PROJECT_ROOT, '..') +  "/logs/iccbl-testing.log",
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
#         '': {
#             'handlers': ['console'],
#             'propagate': False,
#             'level': 'INFO',
#         },
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
        'db.api': {  # set a default handler
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
        'reports.api': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },
        # suppress streaming errors (i.e. image not found)
        'reports.serialize.streaming_serializers': {  
            'handlers': ['logfile'],
            'propagate': True,
            'level': 'WARN',
        },        
        'db.views': {  
            'handlers': ['logfile'],
            'propagate': True,
            'level': 'WARN',
        },        
        'db.tests': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
        'lims.tests': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
        },        
        'django': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
        'django.db.backends': {  # for SQL
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

TEST_RUNNER='reports.tests.IccblTestRunner'

# disable migrations while testing
# see http://stackoverflow.com/questions/25161425/disable-migrations-when-running-unit-tests-in-django-1-7
class DisableMigrations(object):

    def __contains__(self, item):
        return True

    def __getitem__(self, item):
        return "notmigrations"

if 'test' in sys.argv[1:] or 'travis' in sys.argv[1:]:
    print 'tests in progres, no migrations...'
    MIGRATION_MODULES = DisableMigrations()