import logging
from django.shortcuts import render
from django.contrib.auth import logout
from django.http import HttpResponseRedirect

logger = logging.getLogger(__name__)


def main(request):
    logger.info('retrieve home for user: %r: %r', 
        request.user, request.user.is_authenticated)
    return render(request, 'lims/index_require.html', {})

def logout_page(request):
    """
    Log users out and re-direct them to the main page.
    """
    logout(request)
    return HttpResponseRedirect('/lims') # todo: use system property

