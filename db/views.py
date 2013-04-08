import logging
from django.shortcuts import render

logger = logging.getLogger(__name__)


def main(request):
    search = request.GET.get('search','')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index.html', {'search':search })
