# overlay settings for running scripts from the command line on the server
# set's up a different log file so that the server logs aren't overwritten, 
# and don't get the permissions changed.

from os import environ

try:
    from settings import *
except ImportError:
    import sys
    print >>sys.stderr, '''Settings not defined.  Please configure a version of
    settings.py for this site.'''
    del sys
    
    
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
        'console':{
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple'
        },  
    },
    'loggers': {
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
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
        'utils': {  # for SQL
            'handlers': ['console'],
            'propagate': True,
            'level': 'INFO',
        },        
        'tastypie': {  # set a default handler
            'handlers': ['console'],
            'propagate': False,
            'level': 'INFO',
        },        
    }
}

