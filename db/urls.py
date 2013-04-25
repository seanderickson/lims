from django.conf.urls import patterns, url, include
from db import views
from tastypie.api import Api

from db.api import ScreensaverUserResource

v1_api = Api(api_name='v1')
v1_api.register(ScreensaverUserResource())

urlpatterns = patterns('',
    url(r'^$', views.main, name="home"),
    url(r'^screeners/$', views.screeners, name="screeners"),
    url(r'^screeners1/$', views.screeners1, name="screeners1"),
    url(r'^staff/$', views.staff, name="staff"),
    url(r'^screeners_datatables/$', views.screeners_datatables, name="screeners_datatables"),
    (r'^api/', include(v1_api.urls)),
)