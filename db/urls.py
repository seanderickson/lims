from __future__ import unicode_literals

from django.conf.urls import patterns, url, include
from tastypie.api import Api

from db import views
from db.api import ScreensaverUserResource, ScreenResource, \
    ScreenResultResource, StudyResource, \
    DataColumnResource, LibraryResource, \
    LibraryCopyResource, LibraryCopyPlateResource, \
    WellResource, ActivityResource, ReagentResource, \
    SmallMoleculeReagentResource, SilencingReagentResource, NaturalProductReagentResource, \
    CopyWellResource, UserChecklistItemResource, \
    CherryPickRequestResource, CherryPickPlateResource, \
    AttachedFileResource, ServiceActivityResource, LibraryScreeningResource,\
    CherryPickLiquidTransferResource, CherryPickScreeningResource, \
    UserAgreementResource, PublicationResource,PlateLocationResource

import db.api


v1_api = Api(api_name='v1')
v1_api.register(ScreensaverUserResource())
v1_api.register(ScreenResource())
v1_api.register(StudyResource())
v1_api.register(ScreenResultResource())
v1_api.register(DataColumnResource())
v1_api.register(LibraryResource())
v1_api.register(LibraryCopyResource())
v1_api.register(LibraryCopyPlateResource())
v1_api.register(PlateLocationResource())
v1_api.register(WellResource())
v1_api.register(ActivityResource())
v1_api.register(ServiceActivityResource())
v1_api.register(ReagentResource())
# v1_api.register(CherryPickRequest2Resource())
v1_api.register(CherryPickRequestResource())
# v1_api.register(CherryPickWellResource())
v1_api.register(CherryPickPlateResource())
v1_api.register(UserChecklistItemResource())
v1_api.register(AttachedFileResource())
v1_api.register(SmallMoleculeReagentResource())
v1_api.register(SilencingReagentResource())
v1_api.register(NaturalProductReagentResource())
v1_api.register(CopyWellResource())
v1_api.register(LibraryScreeningResource())
v1_api.register(CherryPickLiquidTransferResource())
v1_api.register(CherryPickScreeningResource())
v1_api.register(UserAgreementResource())
v1_api.register(PublicationResource())
v1_api.register(db.api.ResourceResource())

urlpatterns = patterns('',
    url(r'^$', views.main, name="home"),
    url(r'^smiles_image/(?P<well_id>\S+)$','db.views.smiles_image', name="smiles_image" ),
    url(r'^well_image/(?P<well_id>\S+)$','db.views.well_image', name="well_image" ),
    url(r'^attachedfile/(?P<attached_file_id>\d+)/content$',
        'db.views.attached_file', name="attached_file" ),
    url(r'^publication/(?P<publication_id>\d+)/attached_file$',
        'db.views.publication_attached_file', name="publication_attached_file" ),
    (r'^api/', include(v1_api.urls)),
)