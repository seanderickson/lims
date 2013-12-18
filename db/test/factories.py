import factory
import db.models

class ScreenFactory(factory.django.DjangoModelFactory):
    FACTORY_FOR=db.models.Screen

    data_sharing_level = 1
    facility_id = factory.Sequence(lambda n: str(n))
    project_phase = "Primary Screen"
    screen_type = "Small Molecule"
    title = "Test screen 1"


class LibraryFactory(factory.django.DjangoModelFactory):
    FACTORY_FOR=db.models.Library
    
    short_name = factory.Sequence(lambda n: 'library_'+ str(n) )
    library_name = factory.Sequence(lambda n: 'library_'+ str(n) + '_long' )
    screen_type = 'smallmolecule'
    library_type = 'smallmolecule'
    start_plate = factory.Sequence(lambda n: str(n*10) )
    end_plate = factory.Sequence(lambda n: str(n*10+9) )
    plate_size = str(384)
        