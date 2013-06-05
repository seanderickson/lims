from django.conf.urls import patterns, url, include
from db import views
from tastypie.api import Api

from db.api import ScreensaverUserResource, ScreenResource, LabHeadResource, LabAffiliationResource, ScreeningRoomUserResource

v1_api = Api(api_name='v1')
v1_api.register(ScreensaverUserResource())
v1_api.register(ScreenResource())
v1_api.register(LabHeadResource())
v1_api.register(LabAffiliationResource())
v1_api.register(ScreeningRoomUserResource())

urlpatterns = patterns('',
    url(r'^$', views.main, name="home"),
    url(r'^screeners/$', views.screeners, name="screeners"),
    url(r'^screeners1/$', views.screeners1, name="screeners1"),
    url(r'^screeners2/$', views.screeners2, name="screeners2"),
    url(r'^staff/$', views.staff, name="staff"),
    url(r'^screeners_datatables/$', views.screeners_datatables, name="screeners_datatables"),
    url(r'^screens_sm/$', views.screens_sm, name="screens_sm"),
    (r'^api/', include(v1_api.urls)),
)