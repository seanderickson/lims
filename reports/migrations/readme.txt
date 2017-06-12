201706012 - To make initial migration:
- avoid circular dependency between reports and db apps:
	i.e. errors like: ValueError: Lookup failed for model referenced by field db.ScreensaverUser.user: reports.UserProfile
- solution: temporarily remove "db" from lims/base_settings.py "INSTALLED_APPS"
