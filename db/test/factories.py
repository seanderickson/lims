import factory
import db.models
from django.utils.timezone import now
from factory.fuzzy import FuzzyChoice, FuzzyInteger

class ScreenFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = db.models.Screen
    data_sharing_level = 1
    screen_type = "small_molecule"
    title = factory.Sequence(lambda n: 'Screen Title_'+ str(n) )
    summary = factory.Sequence(lambda n: 'Screen summary No_'+ str(n) )
    publishable_protocol = factory.Sequence(lambda n: 'Screen protocol No_'+ str(n) )

class ScreensaverUserFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = db.models.ScreensaverUser

    username = factory.Sequence(lambda n: 'username_'+ str(n) )
    first_name = factory.Sequence(lambda n: 'first_'+ str(n) )
    last_name = factory.Sequence(lambda n: 'last_'+ str(n) )
    email = factory.Sequence(lambda n: 'testemail_%d@testemail.com' % n )
    is_active = True
    
class LibraryFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = db.models.Library
    
    short_name = factory.Sequence(lambda n: 'library_'+ str(n) )
    library_name = factory.Sequence(lambda n: 'library_'+ str(n) + '_long' )
    screen_type = 'small_molecule'
    screening_status = 'allowed'
    provider = 'test provider'
    solvent = 'dmso'
    library_type = 'commercial'
    start_plate = factory.Sequence(lambda n: str(n*10) )
    end_plate = factory.Sequence(lambda n: str(n*10+9) )
    plate_size = str(384)
#     date_created = now()
        
class WellFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = db.models.Well
    plate_number = '{:0>5d}'.format(FuzzyInteger(0, 55000).fuzz())
    col = FuzzyInteger(0, 23).fuzz()
    row = FuzzyInteger(0, 16).fuzz()
    well_name = chr(ord('A')+row) + '{:0>2d}'.format(col)
    well_id = '%s:%s' % (plate_number, well_name)

    facility_id = 'HMSL{:0>5d}'.format(FuzzyInteger(0, 16).fuzz())
    library_well_type = FuzzyChoice(['experimental','empty','dmso','library_control','rnai_buffer'])
    
        
class ReagentFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = db.models.Reagent
    vendor_identifier = factory.Sequence(lambda n: 'rvi_'+ str(n) )
    vendor_name = factory.Sequence(lambda n: 'vendor_'+ str(n) )
    vendor_batch_id = factory.Sequence(lambda n: 'vendor_batch_'+ str(n) + '_long' )
        
class SmallMoleculeReagentFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = db.models.SmallMoleculeReagent
    smiles = factory.Sequence(lambda n: 'CC[C@H](CO)NC1=NC2=C(C(=N1)NCC3=CC=CC=C3)N=CN2C(C)C%d' % n )
    moldata = r'''Structure89
csChFnd70/04290511482D

 16 17  0  0  0  0  0  0  0  0999 V2000
    4.8373    2.0813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8373    0.7093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6256    2.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0491    2.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4179    2.0813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4179    0.7093    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.0813    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2308    6.2888    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2308    4.9141    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2076    2.7850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2076    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0491    4.2034    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.2308    2.1223    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.0491    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.6256    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0  0  0  0
  3  1  2  0  0  0  0
 16  2  1  0  0  0  0
  5  1  1  0  0  0  0
  6  3  1  0  0  0  0
  7  6  2  0  0  0  0
 15  2  2  0  0  0  0
 14  5  2  0  0  0  0
 13  5  1  0  0  0  0
 11  6  1  0  0  0  0
 12  7  1  0  0  0  0
 10 13  1  0  0  0  0
  9 10  1  0  0  0  0
  8 11  2  0  0  0  0
  4  8  1  0  0  0  0
 16  7  1  0  0  0  0
  4 12  2  0  0  0  0
M  END'''         
    
class SubstanceFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = db.models.Substance

class LabAffiliationFactory(factory.Factory):

    category = factory.Sequence(lambda n: 'LAC_'+ str(n) )
    name = factory.Sequence(lambda n: 'Lab Title '+ str(n) )    
    ordinal = factory.Sequence(lambda n: n)
    
    