from django.conf.urls import patterns, url, include
from reports import views
from tastypie.api import Api

from reports.api import MetaHashResource, VocabulariesResource


v1_api = Api(api_name='v1')
v1_api.register(MetaHashResource())
v1_api.register(VocabulariesResource())

urlpatterns = patterns('',
    url(r'^$', views.main, name="home"),
    (r'^api/', include(v1_api.urls)),
)