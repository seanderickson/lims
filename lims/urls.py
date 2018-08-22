from __future__ import unicode_literals
from django.conf.urls import include, url
from lims import views

from django.contrib import admin
import django.contrib.auth.views

admin.autodiscover()

urlpatterns = [
    url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    url(r'^admin/', admin.site.urls),

    # Login / logout.
    # Note: login.html is actually served by the reports project:
    # that is, reports/templates/login.html version; 
    # this is done because at this time only the reports project has
    # all of the necessary css and javascript installed
    
    ## Django > 1.11
    url(r'^accounts/login/$',
        django.contrib.auth.views.LoginView.as_view(
            template_name='login.html',
            redirect_authenticated_user=True),
        name='login'),
    url(r'^accounts/logout/$', 
        django.contrib.auth.views.LogoutView.as_view(), name='logout'),
        
    ## Django < 1.11 
    # url(r'^accounts/login/$', django.contrib.auth.views.login, 
    #     {'template_name': 'login.html'}, name='login'),
    # url(r'^accounts/logout/$', views.logout_page, name='logout'),
    # 
    # url(r'^lims/$', views.main, name="home"),
    # url(r'^db/', include('db.urls')),
    # url(r'^reports/', include('reports.urls')),
]

