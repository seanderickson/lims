# ICCB-L (HMS) specific instructions
[![Build Status](https://travis-ci.org/hmsiccbl/lims.svg?branch=master)](https://travis-ci.org/hmsiccbl/lims.svg?branch=master)
# To run manage.py on the orchestra file system:
- PGSQL_SCREENSAVER_DB, _USER, _HOST, _PASSWORD variables must be in the environment; the settings.py looks for them
- copy lims/settings-orchestra-iccbl-hms.py to lims/settings.py
- then, the cgi-bin/django.cgi script will be run with these variables set by the orchestra system
- to run manage.py from the command line, use the "setenv_and_run.sh" script to grab the vars from the apache/conf/auth file:
./setenv_and_run.sh /opt/apache/conf/auth/dev.screensaver2.med.harvard.edu ./manage.py <other args>

# To reset/regenerate the migrations (useful in reports, before final deploy, to reset)
mv reports/migrations reports/migrations.bak
./manage.py schemamigration reports --initial
./manage.py migrate reports 0001 --fake --delete-ghost-migrations
