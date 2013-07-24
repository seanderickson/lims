
"""
About testing:

From: https://docs.djangoproject.com/en/dev/topics/testing/overview/s
The test database
Tests that require a database (namely, model tests) will not use your "real"
(production) database. Separate, blank databases are created for the tests.

Regardless of whether the tests pass or fail, the test databases are destroyed 
when all the tests have been executed.

By default the test databases get their names by prepending test_ to the value 
of the NAME settings for the databases defined in DATABASES. When using the 
SQLite database engine the tests will by default use an in-memory database 
(i.e., the database will be created in memory, bypassing the filesystem entirely!). 
If you want to use a different database name, specify TEST_NAME in the dictionary 
for any given database in DATABASES.

Aside from using a separate database, the test runner will otherwise use all of
the same database settings you have in your settings file: ENGINE, USER, HOST, 
etc. The test database is created by the user specified by USER, so you'll need
to make sure that the given user account has sufficient privileges to create a
new database on the system.

NOTE: to run using sqlite in memory, use a test_settings file with the database as:
DATABASES['default'] = {'ENGINE': 'django.db.backends.sqlite3'}
and run like:
$ ./manage.py test --settings=lims.test_settings

"""

from django.test import TestCase

import datetime
import csv
import StringIO
from django.contrib.auth.models import User
from django.test.client import Client
from tastypie.test import ResourceTestCase, TestApiClient
from reports.models import MetaHash, Vocabularies
from lims.api import CSVSerializer

import logging

logger = logging.getLogger(__name__)

class CreateMetaHash(ResourceTestCase):
    # Use ``fixtures`` & ``urls`` as normal. See Django's ``TestCase``
    # documentation for the gory details.
    fixtures = ['test_entries.json']

    def setUp(self):
        super(CreateMetaHash, self).setUp()

        # Create a user.
        self.username = 'daniel'
        self.password = 'pass'
        self.user = User.objects.create_user(self.username, 'daniel@example.com', self.password)
        
        initializer = {
                       'key': 'key',
                       'scope': 'metahash:fields',
                       'ordinal': 0    }
        MetaHash.objects.create(**initializer)
        
#        # Fetch the ``Entry`` object we'll use in testing.
#        # Note that we aren't using PKs because they can change depending
#        # on what other tests are running.
#        self.entry_1 = Entry.objects.get(slug='first-post')
#
#        # We also build a detail URI, since we will be using it all over.
#        self.detail_url = '/api/v1/entry/{0}/'.format(self.entry_1.pk)
#
#        # The data we'll send on POST requests. Again, because we'll use it
#        # frequently (enough).
#        self.post_data = {
#            'user': '/api/v1/user/{0}/'.format(self.user.pk),
#            'title': 'Second Post!',
#            'slug': 'second-post',
#            'created': '2012-05-01T22:05:12'
#        }

    def get_credentials(self):
        return self.create_basic(username=self.username, password=self.password)

    def test_get_list_unauthorzied(self):
        self.assertHttpUnauthorized(self.api_client.get('/reports/api/v1/metahash/', format='json'))

    def test_get_list_json(self):
        print '=============== test get_list ==================='

        resp = self.api_client.get('/reports/api/v1/metahash/', format='json', authentication=self.get_credentials())
        print resp
        self.assertValidJSONResponse(resp)

    def test_post_simple(self):
        print '=============== test test_post_simple ==================='
        
        initializer = {
               'key': 'scope',
               'scope': 'metahash:fields',
               'ordinal': 0    }

        resp = self.api_client.post('/reports/api/v1/metahash/', format='json', data=initializer, authentication=self.get_credentials())
        print resp
        self.assertHttpCreated(resp)

        resp = self.api_client.get('/reports/api/v1/metahash/', format='json', authentication=self.get_credentials())
        print resp
        self.assertValidJSONResponse(resp)
        
    def test_put_json(self):
        print '=============== test put_json ==================='
        # first create a "json" field by specifying the "json_type" of the field 
        # (if this weren't specified, it would throw an error, since there is no "real" field for this (todo: test))
        initializer = {
               'key': 'test_field',
               'scope': 'metahash:fields',
               'ordinal': 0, 
               'json_field_type': 'fields.CharField'    }

        resp = self.api_client.post('/reports/api/v1/metahash/', format='json', data=initializer, authentication=self.get_credentials())
        new_obj = self.deserialize(resp)
        print '----- resp:' , resp, '---', new_obj
        self.assertHttpCreated(resp)

        print '-------------------  test put_json modification ==================='
        # now create an some field instance data
        initializer = {
               'key': 'test_field',
               'scope': 'metahash:fields',
               'ordinal': 0, 
               'json_field_type': 'fields.CharField', 
               'test_field': 'foo and bar!'    }

        uri = '/reports/api/v1/metahash/'+ str(new_obj.get('id')) + '/'
        print 'uri:', uri
        resp = self.api_client.put(uri, data=initializer, authentication=self.get_credentials())
        print '--------resp to put:', resp
        self.assertHttpAccepted(resp)
        new_obj = self.deserialize(resp)
        self.assertEqual(new_obj['test_field'], 'foo and bar!', 'unexpected result: ' +new_obj['test_field'])

        resp = self.api_client.get('/reports/api/v1/metahash/', format='json', authentication=self.get_credentials())
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(self.deserialize(resp)['objects']), 2)
        logger.info(str(('---- resp:' , resp)))
        

    # this should work, but the tastypie.ResourceTestCase methods around post for my custom serializer are broken
    # todo: just test the api directly through the requests api
    def xxxx_test_csv_serialization(self):
        print '=============== test_csv_serialization1 ==================='
        ''' 
        Note: we are only verifying simple property lists, i.e. key-value lists; 
        where values must be: either a string or a list of strings
        '''
        serializer=CSVSerializer() 
        self.api_client = TestApiClient(serializer=serializer) # todo: doesn't work for post, see TestApiClient.post() method, it is incorrectly "serializing" the data before posting
        posting_client = Client()
        
        header = ['key', 'scope', 'ordinal', 'json_field_type']
        vals = ['test_field', 'metahash:fields', 0, 'fields.CharField']
        
        raw_data = StringIO.StringIO()
        writer = csv.writer(raw_data)
        writer.writerow(header)
        writer.writerow(vals)
        
        kwargs = {  'content_type': 'text/csv',
                    'HTTP_AUTHORIZATION': self.get_credentials()
                }

        resp = posting_client.post('/reports/api/v1/metahash/', data=raw_data.getvalue(), **kwargs)
        new_obj = self.deserialize(resp)
        print '----- resp:' , resp, '---', new_obj
        logger.info(str(( '----- resp:' , resp, '---', new_obj)))
        self.assertHttpCreated(resp)

        print '-------------------  test put_json modification ==================='
        # now create an some field instance data
        header = ['key', 'scope', 'test_field']
        vals = ['test_field', 'metahash:fields', 'foo and bar!']
        
        raw_data = StringIO.StringIO()
        writer = csv.writer(raw_data)
        writer.writerow(header)
        writer.writerow(vals)

        uri = '/reports/api/v1/metahash/'+ str(new_obj.get('id')) + '/'
        print 'uri:', uri
        resp = self.api_client.put(uri, format='csv', serializer=serializer, data=raw_data.getvalue(), authentication=self.get_credentials())
        print '--------resp to put:', resp
        self.assertHttpAccepted(resp)
        new_obj = self.deserialize(resp)
        self.assertEqual(new_obj['test_field'], 'foo and bar!', 'unexpected result: ' +new_obj['test_field'])

        resp = self.api_client.get('/reports/api/v1/metahash/', format='json', authentication=self.get_credentials())
        self.assertValidJSONResponse(resp)
        self.assertEqual(len(self.deserialize(resp)['objects']), 2)
        logger.info(str(('---- resp:' , resp)))
        
        
        
#        # Scope out the data for correctness.
#        self.assertEqual(len(self.deserialize(resp)['objects']), 12)
#        # Here, we're checking an entire structure for the expected data.
#        self.assertEqual(self.deserialize(resp)['objects'][0], {
#            'pk': str(self.entry_1.pk),
#            'user': '/api/v1/user/{0}/'.format(self.user.pk),
#            'title': 'First post',
#            'slug': 'first-post',
#            'created': '2012-05-01T19:13:42',
#            'resource_uri': '/api/v1/entry/{0}/'.format(self.entry_1.pk)
#        })

#    def test_get_list_xml(self):
#        self.assertValidXMLResponse(self.api_client.get('/api/v1/entries/', format='xml', authentication=self.get_credentials()))

#    def test_get_detail_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.get(self.detail_url, format='json'))
#
#    def test_get_detail_json(self):
#        resp = self.api_client.get(self.detail_url, format='json', authentication=self.get_credentials())
#        self.assertValidJSONResponse(resp)
#
#        # We use ``assertKeys`` here to just verify the keys, not all the data.
#        self.assertKeys(self.deserialize(resp), ['created', 'slug', 'title', 'user'])
#        self.assertEqual(self.deserialize(resp)['name'], 'First post')
#
##    def test_get_detail_xml(self):
##        self.assertValidXMLResponse(self.api_client.get(self.detail_url, format='xml', authentication=self.get_credentials()))
#
#    def test_post_list_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.post('/api/v1/entries/', format='json', data=self.post_data))
#
#    def test_post_list(self):
#        # Check how many are there first.
#        self.assertEqual(Entry.objects.count(), 5)
#        self.assertHttpCreated(self.api_client.post('/api/v1/entries/', format='json', data=self.post_data, authentication=self.get_credentials()))
#        # Verify a new one has been added.
#        self.assertEqual(Entry.objects.count(), 6)
#
#    def test_put_detail_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.put(self.detail_url, format='json', data={}))
#
#    def test_put_detail(self):
#        # Grab the current data & modify it slightly.
#        original_data = self.deserialize(self.api_client.get(self.detail_url, format='json', authentication=self.get_credentials()))
#        new_data = original_data.copy()
#        new_data['title'] = 'Updated: First Post'
#        new_data['created'] = '2012-05-01T20:06:12'
#
#        self.assertEqual(Entry.objects.count(), 5)
#        self.assertHttpAccepted(self.api_client.put(self.detail_url, format='json', data=new_data, authentication=self.get_credentials()))
#        # Make sure the count hasn't changed & we did an update.
#        self.assertEqual(Entry.objects.count(), 5)
#        # Check for updated data.
#        self.assertEqual(Entry.objects.get(pk=25).title, 'Updated: First Post')
#        self.assertEqual(Entry.objects.get(pk=25).slug, 'first-post')
#        self.assertEqual(Entry.objects.get(pk=25).created, datetime.datetime(2012, 3, 1, 13, 6, 12))
#
#    def test_delete_detail_unauthenticated(self):
#        self.assertHttpUnauthorized(self.api_client.delete(self.detail_url, format='json'))
#
#    def test_delete_detail(self):
#        self.assertEqual(Entry.objects.count(), 5)
#        self.assertHttpAccepted(self.api_client.delete(self.detail_url, format='json', authentication=self.get_credentials()))
#        self.assertEqual(Entry.objects.count(), 4)