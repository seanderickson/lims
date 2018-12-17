from __future__ import unicode_literals 

import logging

from django.conf import settings
from django.contrib.auth import logout
from django.http import HttpResponseRedirect
from django.shortcuts import render


logger = logging.getLogger(__name__)


def main(request):
    
    public_data = settings.APP_PUBLIC_DATA.as_dict()
    return render(request, 'lims/index.html', public_data)

def logout_page(request):
    """
    Log users out and re-direct them to the main page.
    """
    logout(request)
    return HttpResponseRedirect('/lims') # todo: use system property

