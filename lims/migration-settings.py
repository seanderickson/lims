# Django settings for lims project

try:
    from settings import *
except ImportError:
    import sys
    print >>sys.stderr, '''Settings not defined.  Please configure a version of
    settings.py for this site.'''
    del sys
    
import os.path


print 'PROJECT_ROOT: ', PROJECT_ROOT

# Turn off these settings so bootstrapping to a local non-HTTPS server will work
BACKGROUND_PROCESSING = False
IS_PRODUCTION_READY = False
SECURE_SSL_REDIRECT = False
SESSION_COOKIE_SECURE = False
CSRF_COOKIE_SECURE = False

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
            'filename': 
                os.path.join(PROJECT_ROOT, 'logs', 'screensaver-migration.log'),
            'maxBytes': 5000000,
            'backupCount': 2,
            'formatter': 'simple',
        },
    },
    'loggers': {
        '': {  # set a default handler
            'handlers': ['logfile'],
            'propagate': True,
            'level': 'INFO',
        },        
        'django.request': {
            'handlers': ['mail_admins'],
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
        'django.db': {  # for SQL
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
        'django.db.backends': {  # for SQL
            'handlers': ['logfile'],
            'propagate': False,
            'level': 'INFO',
        },        
        'utils': {  # 
            'handlers': ['logfile'],
            'propagate': True,
            'level': 'INFO',
        },        
    }
}
