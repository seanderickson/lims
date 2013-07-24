from settings import *

# make tests faster
# use from the command line with testing like
# ./manage.py test --settings=lims.test_settings
SOUTH_TESTS_MIGRATE = False
DATABASES['default'] = {'ENGINE': 'django.db.backends.sqlite3',
                        'NAME': ':memory'}
