import logging
from django.shortcuts import render

logger = logging.getLogger(__name__)


def main(request):
    search = request.GET.get('search','')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index.html', {'search':search })


def smallmolecule(request):
    search = request.GET.get('search','')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index_jqgrid.html', {'search':search })

def smallmolecule1(request):
    search = request.GET.get('search','')
    logger.debug(str(('index_jqgrid_notemplate  search: ', search)))
    return render(request, 'db/index_jqgrid_notemplate.html', {'search':search })

def smallmoleculelisting(request):
    search = request.GET.get('search','')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index_jqgrid.html', {'search':search })

