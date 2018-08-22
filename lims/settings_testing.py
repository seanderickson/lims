# Django settings for 1km project
import sys

try:
    from settings import *
except ImportError:
    print >>sys.stderr, '''Settings not defined.  Please configure a version of
    settings.py for this site.'''
    
import os.path


print 'PROJECT_ROOT: ', PROJECT_ROOT

# set SQLALCHEMY_POOL_CLASS=sqlalchemy.pool.NullPool for testing
# environments, so that the test database can be destroyed
import sqlalchemy.pool
SQLALCHEMY_POOL_CLASS = sqlalchemy.pool.NullPool

# FIXME: sqllite3 db does not work - errors on "DISTINCT ON" clause
# DATABASES['default'] = {'ENGINE': 'django.db.backends.sqlite3',
#                         'NAME': ':memory'}
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': 'test_screensaverlims',
        'USER': '',
        'PASSWORD': '',
        'HOST': '',
        'PORT': '',
    }
}

BACKGROUND_PROCESSING = False

TEST_RUNNER='reports.tests.IccblTestRunner'

# disable migrations while testing
# see http://stackoverflow.com/questions/25161425/disable-migrations-when-running-unit-tests-in-django-1-7
class DisableMigrations(object):

    def __contains__(self, item):
        return True

    def __getitem__(self, item):
        return None

if 'test' in sys.argv[1:] or 'travis' in sys.argv[1:]:
    print 'tests in progres, no migrations...'
    MIGRATION_MODULES = DisableMigrations()
