# Django settings to run tests with no db
# NOTE: the test runner finds classes differently with this runner:
# usually use:
# ./manage.py test reports.tests.SDFSerializerTest.test2_clean_data_sdf --settings=lims.settings_testing_debug
# but now use:
# ./manage.py test reports.SDFSerializerTest.test2_clean_data_sdf --settings=lims.settings_testing_debug
# >NOTE< the second form does not include the module "test"
#
#

try:
    from settings import *
except ImportError:
    import sys
    print >>sys.stderr, '''Settings not defined.  Please configure a version of
    settings.py for this site.'''
    del sys
    
import os.path


print 'PROJECT_ROOT: ', PROJECT_ROOT, ', ' , os.path.join(PROJECT_ROOT, '..')


# make tests faster
# use from the command line with testing like
# ./manage.py test --settings=lims.testing-settings
SOUTH_TESTS_MIGRATE = False

# set SQLALCHEMY_POOL_CLASS=sqlalchemy.pool.NullPool for testing
# environments, so that the test database can be destroyed
import sqlalchemy.pool
SQLALCHEMY_POOL_CLASS = sqlalchemy.pool.NullPool

# import reports.tests.NoDbTestRunner
TEST_RUNNER='reports.tests.NoDbTestRunner'

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
            'format': '%(levelname)s %(msecs)s: %(name)s:%(funcName)s:%(lineno)d: %(message)s'
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
            'handlers': ['console'],
            'level': 'ERROR',
            'propagate': False,
        },
        'db': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'lims': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports': {  # set a default handler
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },
        'reports.api': {  # set a default handler
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },
        'db.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
        },        
        'lims.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
        },        
        'django': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'utils': {  # for SQL
            'handlers': ['console'],
            'propagate': True,
            'level': 'DEBUG',
        },        
        'tastypie': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
    }
}
