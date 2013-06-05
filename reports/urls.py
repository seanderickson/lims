from django.conf.urls import patterns, url, include
from db import views
from tastypie.api import Api

from reports.api import FieldInformationResource


v1_api = Api(api_name='v1')
v1_api.register(FieldInformationResource())

urlpatterns = patterns('',
#    url(r'^$', views.main, name="home"),
    (r'^api/', include(v1_api.urls)),
)