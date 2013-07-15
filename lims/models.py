from django.db import models


class GetOrNoneManager(models.Manager):
    """Adds get_or_none method to objects
    """
    def get_or_none(self, function=None, **kwargs):
        try:
            x = self.get(**kwargs)
            if x and function:
                return function(x)
            else:
                return x
        except self.model.DoesNotExist:
            return None
