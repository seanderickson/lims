from __future__ import unicode_literals
from django.conf.urls import url, include

from reports.api_base import Api
from reports.api import VocabularyResource, ResourceResource,\
    ApiLogResource, UserResource, UserGroupResource, PermissionResource, \
    FieldResource, JobResource

v1_api = Api(api_name='v1')
v1_api.register(VocabularyResource())
v1_api.register(ResourceResource())
v1_api.register(ApiLogResource())
v1_api.register(UserResource())
v1_api.register(UserGroupResource())
v1_api.register(PermissionResource())
v1_api.register(JobResource())

v1_api.register(FieldResource())

# #     url(r'^$', views.main, name="home"),
urlpatterns = [
    url(r'^api/', include(v1_api.urls)),
]