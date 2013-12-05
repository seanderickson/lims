import factory
import db.models

class ScreenFactory(factory.django.DjangoModelFactory):
    FACTORY_FOR=db.models.Screen
    data_sharing_level = 1
    facility_id = factory.Sequence(lambda n: str(n))
    project_phase = "Primary Screen"
    screen_type = "Small Molecule"
    title = "Test screen 1"
