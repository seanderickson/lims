from __future__ import unicode_literals

import json
import logging

from django.http import HttpResponse
from django.shortcuts import render
from django.utils.encoding import smart_str


logger = logging.getLogger(__name__)

#  
# def main(request):
#     search = request.GET.get('search', '')
#     logger.debug(str(('main search: ', search)))
#     return render(request, 'reports/index_require.html', {'search': search})

