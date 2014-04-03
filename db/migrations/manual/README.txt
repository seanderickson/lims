NOTES: 

#TODO: can the django custom sql hooks?

- strategy will be to have a migrations directory *only* for migrating from 
Screensaver 1.
- For clean installations, we will create a new 0001_initial script just before 
release of the first production version.
- controlled vocabularies are being reworked: 
thus all instances of controlled vocabularies will require a conversion script.    
- activities are being reworked: 
thus all instances of activities will require conversion scripts (see the 
screen status conversion scripts).
    
For de-novo databases:
- we will create a final "0001_initial" migration script; it will have the most current (at release time)
representation of the schema.
- installers on de-novo installations will run:
./manage.py migrate myapp 0001 --fake

For transformation installations (Screensaver 1 to 2):

- the schema will still
The stock "initial" should be re-created at release time and will not need these scripts at all.
    
Other:
* These migrations must be run before the corresponding South migrations are run *
So "0003_screen_status.sql" is run before "0003_screen_status.py"
