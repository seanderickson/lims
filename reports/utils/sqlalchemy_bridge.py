'''
Based on "django-sabridge" at https://github.com/johnpaulett/django-sabridge
Copyright (c) 2011, John Paulett
All rights reserved.
'''
from sqlalchemy import create_engine, MetaData, Table
from django.db import connection
from django.conf import settings
import sqlalchemy.pool
import logging

logger = logging.getLogger(__name__)

class Bridge(object):
    def __init__(self):
        self._meta = None
        self._tables = {}

    def connection_url(self):
        """Build a URL for :py:func:`sqlalchemy.create_engine`
        based upon the database defined by :py:attr:`django.db.connection`
        """
        # we lazily import connection since it requires a DJANGO_SETTINGS_MODULE
        
        return urlbuild(
            scheme=connection.vendor,
            path=connection.settings_dict['NAME'],
            username=connection.settings_dict['USER'],
            password=connection.settings_dict['PASSWORD'],
            hostname=connection.settings_dict['HOST'],
            port=connection.settings_dict['PORT']
        )
    
    def __getitem__(self, db_table):
        # not thread-safe, but at worst, we re-autoload the table
        if db_table not in self._tables:
            self._tables[db_table] = self._map_model(db_table)
        return self._tables[db_table]

    def _map_model(self, db_table):
        return Table(db_table, self.meta, autoload=True)
    
    def get_engine(self):
        return self.meta.bind
    
    @property
    def meta(self):
        """:py:class:`sqlalchemy.schema.MetaData` instance bound to the current
        Django database connection."""
        # Lazily connect to the database. By waiting until meta is
        # needed, we prevent capturing stale connection info from Django
        # (e.g. during a TestCase)
        if self._meta is None:
            self._meta = MetaData()
            # set SQLALCHEMY_POOL_CLASS == sqlalchemy.pool.NullPool for testing
            # environments, so that the test database can be destroyed
            if getattr(settings, 'SQLALCHEMY_POOL_CLASS', None):
                self._meta.bind = create_engine(self.connection_url(), poolclass=settings.SQLALCHEMY_POOL_CLASS)
            else:
                self._meta.bind = create_engine(self.connection_url())
        return self._meta


def urlbuild(scheme, path, username=None, password=None, hostname=None, port=None):
    """Builds a URL to pass into :py:func:`sqlalchemy.create_engine`

    dialect+driver://username:password@host:port/database
    """
    # we could use urlparse.ParseResult, but we would still need to manually
    # build the whole netloc, so we don't gain much
    
    # build the netloc depending on what parts we are given,
    # later join the netloc together
    netloc = []
    if username:
        netloc.append(username)
    if password:
        netloc.append(':')
        netloc.append(password)
    if username or password:
        netloc.append('@')
    if hostname:
        netloc.append(hostname)
    if port:
        netloc.append(':')
        netloc.append(str(port))

    # put it all together
    return '{scheme}://{netloc}/{path}'.format(scheme=scheme,
                                               netloc=''.join(netloc),
                                               path=path)
