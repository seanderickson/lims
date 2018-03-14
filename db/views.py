# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import json
import logging
import os.path
import sys
import types

from PIL import Image
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponse
from django.http.response import Http404, HttpResponseServerError
from django.shortcuts import render

from db import WELL_ID_PATTERN
from db.api import AttachedFileAuthorization, PublicationAuthorization
from db.models import ScreensaverUser, Reagent, AttachedFile, Publication, \
    RawDataTransform, SmallMoleculeReagent
from reports.api import UserGroupAuthorization
from reports.serialize import XLSX_MIMETYPE


logger = logging.getLogger(__name__)


def main(request):
    search = request.GET.get('search', '')
    logger.debug(str(('main search: ', search)))
    return render(request, 'db/index.html', {'search': search})

@login_required
def well_image(request, well_id):

    # TODO: use group authorization - not required at ICCBL
    auth = UserGroupAuthorization('reagent')
    # if auth.has_read_authorization(
    #     request.user, well_id) is not True:
    #     return HttpResponse(status=403)
    match = WELL_ID_PATTERN.match(well_id)
    if not match:
        logger.warn('invalid well_id format: %d, pattern: %s' 
            % (well_id,WELL_ID_PATTERN))
        raise Http404('invalid well id format: %s' % well_id)
    else:
        _plate = match.group(1)
        _well_name = match.group(2)
        well_id = '%s:%s' % (_plate,_well_name)
        if auth._is_resource_authorized(
            request.user, 'read', resource_name='reagent') is not True:
            smr = SmallMoleculeReagent.objects.get(
                well_id=well_id)
            if smr.is_restricted_structure is True:
                logger.warn('structure is restricted: %s', well_id)
                raise Http404
        
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
            logger.info('well_image for %s not found at "%s"', 
                well_id, structure_image_path)
            raise Http404

@login_required
def smiles_image(request, well_id):
    # TODO: not tested
        
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


@login_required
def publication_attached_file(request, publication_id):
    response = None
    try:
        publication = Publication.objects.get(publication_id=publication_id)
        auth = PublicationAuthorization('publication')
        if auth.has_publication_read_authorization(
            request.user, publication_id) is not True:
            response = HttpResponse(status=403)
        else:
            response = _download_file(request,publication.attachedfile)
    except ObjectDoesNotExist,e:
        msg = 'could not find publication object for id: %r' % publication_id
        logger.exception(msg)
        response = HttpResponse(status=404)
    return process_response(request, response)

@login_required
def attached_file(request, attached_file_id):
    logger.info('request: %r', request)
    response = None
    try:
        logger.info('get attached file: %r, %r', request, attached_file_id)
        af = AttachedFile.objects.get(attached_file_id=attached_file_id)
        logger.info('found: %r', af)
        auth = AttachedFileAuthorization('attachedfile')
        if auth.has_file_read_authorization(
                request.user, attached_file_id) is not True:
            logger.info('no read access...')
            response = HttpResponse(status=403)
        else:
            logger.info('creating response...')
            response = _download_file(request,af)
            logger.info('response: %r', response)
    except ObjectDoesNotExist,e:
        msg = 'AttachedFile not found: %r' % attached_file_id
        logger.exception(msg)
        response = HttpResponse(status=404)

    return process_response(request, response)
        
@login_required
def _download_file(request, attached_file):   
    """                                                                         
    TODO: Send a file through Django without loading the whole file into              
    memory at once. The FileWrapper will turn the file object into an           
    iterator for chunks of 8KB.                                                 
    """
    logger.debug('download attached file: %s', attached_file)
    try:     
        response = HttpResponse(
            str(attached_file.contents), 
            content_type='application/octet-stream') 
        response['Content-Disposition'] = \
            'attachment; filename=%s' % unicode(attached_file.filename)
        response['Content-Length'] = len(attached_file.contents)
        logger.info('download file response %r', response)
        return response
    except Exception,e:
        logger.exception('on accessing attached file %s' % attached_file)
        raise

@login_required
def screen_raw_data_transform(
        request,screen_facility_id):
    logger.info('download raw_data_transform result for screen %r', 
        screen_facility_id)
    
    response = None
    try:
        rdt = RawDataTransform.objects.get(
            screen__facility_id=screen_facility_id)
        if not os.path.exists(rdt.temp_output_filename):
            logger.info('raw data file for %d not found at "%s"', 
                rdt, rdt.temp_output_filename)
            response = HttpResponse(status=404)
        else:
            logger.info('attempt to stream file: %r', rdt.temp_output_filename)
            response = stream_raw_data_file(
                request, rdt.temp_output_filename,rdt.output_filename)
    except ObjectDoesNotExist:
        response = HttpResponse(status=404)
        
    return process_response(request, response)

@login_required
def cpr_raw_data_transform(
        request,cherry_pick_request_id):
    logger.info('download raw_data_transform result for cpr %r', 
        cherry_pick_request_id)

    response = None
    try:
        rdt = RawDataTransform.objects.get(
            cherry_pick_request__cherry_pick_request_id=cherry_pick_request_id)
        if not os.path.exists(rdt.temp_output_filename):
            logger.info('raw data file for %d not found at "%s"', 
                rdt, rdt.temp_output_filename)
            response = HttpResponse(status=404)
        else:
            logger.info('attempt to stream file: %r', rdt.temp_output_filename)
            response = stream_raw_data_file(
                request, rdt.temp_output_filename,rdt.output_filename)
    except ObjectDoesNotExist:
        response = HttpResponse(status=404)
    
    return process_response(request, response)

def stream_raw_data_file(request, filename, output_filename):

    with open(filename) as temp_file:
        temp_file.seek(0, os.SEEK_END)
        size = temp_file.tell()
        temp_file.seek(0)   
        response = HttpResponse(temp_file.read())
        response['Content-Length'] = size
        response['Content-Type'] = XLSX_MIMETYPE
        response['Content-Disposition'] = \
            'attachment; filename=%s.xlsx' % unicode(output_filename)
        return response
    
def process_response(request, response):

    if hasattr(request, 'GET'):    
        downloadID = request.GET.get('downloadID', None)
        if downloadID:
            logger.info('set cookie "downloadID" %r', downloadID )
            response.set_cookie('downloadID', downloadID)
        else:
            logger.debug('no downloadID: %s' % request.GET )
    
    return response
        
    
