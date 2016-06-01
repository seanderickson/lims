from __future__ import unicode_literals

import json
import logging
import os.path
import sys
import types
from wsgiref.util import FileWrapper

from PIL import Image
from django.conf import settings
from django.forms.models import model_to_dict
from django.http import HttpResponse
from django.http.response import Http404, HttpResponseServerError
from django.shortcuts import render

from db import WELL_ID_PATTERN
from db.models import ScreensaverUser, Reagent, AttachedFile
from reports.api import UserGroupAuthorization


logger = logging.getLogger(__name__)


def main(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index.html', {'search': search})

def well_image(request, well_id):
    
    if not request.user.is_authenticated():
        raise ImmediateHttpResponse(response=HttpForbidden(
            'user %r is not authorized for well_image view' % request.user ))

    match = WELL_ID_PATTERN.match(well_id)
    if not match:
        logger.warn('invalid well_id format: %d, pattern: %s' 
            % (well_id,WELL_ID_PATTERN))
        raise Http404('invalid well id format: %s' % well_id)
    else:
        _plate = match.group(1)
        _well_name = match.group(2)
        _name = '%s%s.png' % (_plate,_well_name)
        structure_image_dir = os.path.abspath(settings.WELL_STRUCTURE_IMAGE_DIR)
        structure_image_path = os.path.join(
            structure_image_dir, _plate, _name)

        if os.path.exists(structure_image_path):
            try:
                image = Image.open(structure_image_path)
                response = HttpResponse(content_type="image/png")
                image.save(response, "PNG")
                return response
            except Exception as e:
                logger.exception('well_image exception for %r, %r' 
                    % (well_id, e))
                raise HttpResponseServerError
        else:
            logger.info('well_image for %s not found at %s', 
                well_id, structure_image_path)
            raise Http404


def smiles_image(request, well_id):
    # TODO: not tested
    if not request.user.is_authenticated():
        raise ImmediateHttpResponse(response=HttpForbidden(
            'user: %r not authorized for smiles_image view' % request.user))
        
    import rdkit.Chem
    import rdkit.Chem.AllChem
    import rdkit.Chem.Draw
    reagent = Reagent.objects.get(well_id=well_id)
    smiles = reagent.smallmoleculereagent.smiles
    logger.info(str((well_id, str(smiles) )))
    m = rdkit.Chem.MolFromSmiles(str(smiles))
    rdkit.Chem.AllChem.Compute2DCoords(m)
    im = rdkit.Chem.Draw.MolToImage(m)
    response = HttpResponse(mimetype="image/png")
    im.save(response, "PNG")
    return response

    
def attached_file(request, attached_file_id):
    af = None
    try:
        af = AttachedFile.objects.get(attached_file_id=attached_file_id)
    except Exception,e:
        msg = 'could not find attached file object for id: %r' % attached_file_id
        logger.exception(msg)
        raise Http404(msg)
    return _download_file(request,af)

        
def _download_file(request, attached_file):   
    """                                                                         
    Send a file through Django without loading the whole file into              
    memory at once. The FileWrapper will turn the file object into an           
    iterator for chunks of 8KB.                                                 
    """
    logger.debug('download attached file: %s', attached_file)
    try:     
        if(not request.user.is_authenticated()):
            logger.warn('access to restricted: user: %r, file: %r',
                request.user,attached_file) 
            return HttpResponse('Log in required.', status=401)
        
        userprofile = request.user.userprofile
        if attached_file.screensaver_user.user != userprofile:
            authorization = UserGroupAuthorization()
            if not authorization._is_resource_authorized(
                    'attachedfiles',request.user,'read'):
                msg = ('%r UserGroupAuthorization needed for user: %r'
                    % ('attachedfiles/read',request.user))
                logger.warn(msg)
                return HttpResponse(msg, status=401)
        else:
            logger.info('User allowed attached file access to own files %s',
                request.user)
        
        response = HttpResponse(
            str(attached_file.contents), content_type='application/octet-stream') 
        response['Content-Disposition'] = \
            'attachment; filename=%s' % unicode(attached_file.filename)
        response['Content-Length'] = len(attached_file.contents)
        return response
    except Exception,e:
        logger.exception('on accessing attached file %s' % attached_file)
        raise

