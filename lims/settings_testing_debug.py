# Django settings for 1km project

try:
    from settings_testing import *
except ImportError:
    import sys
    print >>sys.stderr, '''Settings not defined.  Please configure a version of
    settings.py for this site.'''
    del sys
    
import os.path

print 'PROJECT_ROOT: ', PROJECT_ROOT, ', ' , os.path.join(PROJECT_ROOT, '..')


### uncomment to turn datetime warnings into a runtime exception with traceback
# import warnings
# warnings.filterwarnings(
#         'error', r"DateTimeField .* received a naive datetime",
#         RuntimeWarning, r'django\.db\.models\.fields')



# todo: does this exist?
# TEST_RUNNER='test_utils.test_runners.keep_database'

# just found out that this will force a full stack trace (with DEBUG=False?)
TASTYPIE_FULL_DEBUG=True

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
            'handlers': ['console'],
            'level': 'ERROR',
            'propagate': False,
        },
        'db': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        }, 
        'db.api': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
        },        
        'lims': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },               
        'reports': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },
        'reports.api': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
        },
        'reports.api.detail_uri_kwargs': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'WARN',
        },
        'reports.tests': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
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
        'django': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'django.db.backends': {  # for SQL
            'handlers': ['console'],
            'propagate': False,
            'level': 'DEBUG',
        },        
        'tastypie': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
    }
}
