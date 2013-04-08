from django.conf.urls import patterns, url, include
from db import views

urlpatterns = patterns('',
    url(r'^$', views.main, name="home"),
)