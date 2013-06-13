import logging
from django.shortcuts import render
from db.models import ScreensaverUser
from django.http import HttpResponse
import json
from django.utils.encoding import smart_str

logger = logging.getLogger(__name__)


def main(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index.html', {'search': search})

def screeners(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index_jqgrid.html', {'search': search})

def screeners1(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index_datatables.html', {'search': search})

def screeners2(request):
    search = request.GET.get('search', '')
    logger.debug(str(('screeners2')))
    
    # TODO: all these urls should be located using reverse
    #    root_url = '/db/screeners2/'
    url_schema = '/db/api/v1/screensaveruser/schema/' # TODO: how to use django url tag here
    url = '/db/api/v1/screensaveruser/?format=json' # TODO: how to use django url tag here

    return render(request, 'db/index_backbone.html', {'search': search, 'api_url': url, 'api_url_schema': url_schema })

def staff(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index_jqgrid.html', {'search': search})

def screens_sm(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))

    # TODO: all these urls should be located using reverse
    root_url = '/db/screens_sm/'
    url_schema = '/db/api/v1/screen/schema/' # TODO: how to use django url tag here
    url = '/db/api/v1/screen/?format=json' # TODO: how to use django url tag here

    return render(request, 'db/index_backbone.html', {'search': search, 'root_url': root_url, 'api_url': url, 'api_url_schema': url_schema })

from django.forms.models import model_to_dict

def screeners_datatables(request):
#    logger.info(str(('screeners_datatables:', request)))
    sEcho = 0
    if (request.GET.get('sEcho')):
        sEcho = int(request.GET.get('sEcho'))
    output1 = {
        "sEcho": sEcho,
        "iTotalRecords": 0,
        "iTotalDisplayRecords": 0,
        "aaData": [],
        "aoColumns": [],
    }
    fields = get_fields(ScreensaverUser)
    logger.info(str(('fields', len(fields))))
    for prop in fields:
        output1['aoColumns'].append({'sTitle': prop})

    total = request.GET.get('iDisplayLength')
    if total:
        total = int(total)
    else:
        total = 40
    logger.info(str(('total', total)))
    i=0
    for obj in ScreensaverUser.objects.all():
        output1['aaData'].append([smart_str(getattr(obj, field), 'utf-8', errors='ignore') for field in fields])
        i=i+1
        if i > total:
            break

    output1['iTotalRecords'] = len(ScreensaverUser.objects.all())
    output1['iTotalDisplayRecords'] = total

    json_return = json.dumps(output1)
    logger.info(str(('return', json_return)))
    return HttpResponse(json_return, mimetype="application/x-javascript")


def get_fields(model):
    fields = {}
    
    for field in model._meta.fields:
        if(field.name in fields):
            break
        else:
            fields.setdefault(field.name,field.name)
    return fields.keys()


import types

def get_properties(obj):
    """
    Custom method to grab all of the fields _and_ properties on model instances
    """

    props = []
    MethodWrapperType = type(object().__hash__)

    for slot in dir(obj):
        try:
            attr = getattr(obj, slot)

            if (slot.find('_') == 0 or slot == '__dict__' or slot == '__class__'
                or slot == '__doc__'
                or slot == '__module__'
                or (isinstance(attr, types.BuiltinMethodType)
                    or isinstance(attr, MethodWrapperType))
                or (isinstance(attr, types.MethodType)
                    or isinstance(attr, types.FunctionType))
                or isinstance(attr, types.TypeType)):
                    continue
            else:
                props.append(slot)
        except Exception, e:
            logger.debug(str(('can not introspect', e)))
    return props

