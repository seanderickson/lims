import logging
from django.shortcuts import render
from django.http import HttpResponse
import json
from django.utils.encoding import smart_str

logger = logging.getLogger(__name__)


def main(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'reports/index_require.html', {'search': search})

#def fieldinformation(request):
#    search = request.GET.get('search', '')
#    logger.debug(str(('fieldinformation')))
#    
#    # TODO: all these urls should be located using reverse
#    url_schema = '/reports/api/v1/fieldinformation/schema/' # TODO: how to use django url tag here
#    url = '/reports/api/v1/fieldinformation/?format=json' # TODO: how to use django url tag here
#
#    return render(request, 'reports/index_backbone.html', {'search': search, 'api_url': url, 'api_url_schema': url_schema })
