import logging
import os.path
import sys
from django.shortcuts import render
from db.models import ScreensaverUser, Reagent
from django.http import HttpResponse
import json
from django.utils.encoding import smart_str
from django.conf import settings
from django.http.response import Http404

logger = logging.getLogger(__name__)


def main(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index.html', {'search': search})

def well_image(request, well_id):
    
    if not request.user.is_authenticated():
        raise ImmediateHttpResponse(response=HttpForbidden(
            str(('user',request.user ,'not authorized for smiles_image view'))))

    reagent = Reagent.objects.get(well_id=well_id)
    response = HttpResponse(mimetype="image/png")
    if reagent.smallmoleculereagent:
        well = reagent.well
        _plate = '{:0>5d}'.format(well.plate_number)
        _name = '%s%s.png' %(_plate, well.well_name)
        try:
            structure_image_dir = os.path.abspath(settings.WELL_STRUCTURE_IMAGE_DIR)
            structure_image_path = os.path.join(
                structure_image_dir, _plate, _name)
            
            if os.path.exists(structure_image_path):
                from PIL import Image
                image = Image.open(structure_image_path)
                image.save(response, "PNG")
                return response
            else:
                msg = str(('well_image not found', well_id, structure_image_path))
                logger.info(msg)
                raise Http404(msg)

        except Exception, e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]      
            msg = str(e)
            logger.warn(str(('well_image', well_id, 
                msg, exc_type, fname, exc_tb.tb_lineno)))
            raise e   

def smiles_image(request, well_id):

    if not request.user.is_authenticated():
        raise ImmediateHttpResponse(response=HttpForbidden(
            str(('user',request.user ,'not authorized for smiles_image view'))))
        
    import rdkit.Chem
    import rdkit.Chem.AllChem
    import rdkit.Chem.Draw
#         import matplotlib
# 
# 
    reagent = Reagent.objects.get(well_id=well_id)
    smiles = reagent.smallmoleculereagent.smiles
    logger.info(str((well_id, str(smiles) )))
    m = rdkit.Chem.MolFromSmiles(str(smiles))
#         m = reagent.smallmoleculereagent.molfile.molfile
#         logger.info(str(('molfile', m)))
    rdkit.Chem.AllChem.Compute2DCoords(m)
    im = rdkit.Chem.Draw.MolToImage(m)
#         return im
# #         matplotlib.pyplot.imshow(im)
#         
    response = HttpResponse(mimetype="image/png")
    im.save(response, "PNG")
    return response
    



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

