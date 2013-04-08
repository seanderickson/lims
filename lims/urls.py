from django.conf.urls import patterns, include, url
from lims.views import *

from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^admin/', include(admin.site.urls)),

    # Login / logout.
    # Note: the name "login_url" name is set to the request by the registered hmslincs.context_procesor.login_url_with_redirect
    (r'^db/login/$', 'django.contrib.auth.views.login', {'template_name': 'db/login.html'}),
    url(r'^db/logout/$', logout_page, name='logout'),
    
    url(r'^db/', include('db.urls')),
)

