The  api_init files are used to "bootstrap" the database with meta information
that defines the various Resources available through the API.

- The metahash_fields_resource.csv file defines what a "Resource" record looks
like.
- metahash_resource.data.csv defines actual instances of Resources.

- All of the "fields" files define the fields in the various other resources
that are available through the API.

- Of particular interest is the metahash_fields_initial.csv and the 
metahash_fields_initial.patch.csv.  Together, these files define the fields of 
a "Metahash" record; these Metahash records define what a "field" is.  There
are the two files; the "initial" file is loaded first; it provides the
definition of all the fields, which, once defined, are then used by the
"initial_patch" file. So the first file defines the fields used to define
_every_ resource, even itself.

- api_init_actions.csv is a recipe file used by the "db_init" script, telling
it what order to load each of the api_init files in the bootstrapping process. 
It is run like this:
run the server in a localhost port:
(virtualenv)$./manage.py runserver 55001
run the bootloader script:
(virtualenv)$ ./manage.py db_init  --inputDirectory=./reports/static/api_init/ \
  -f ./reports/static/api_init/api_init_actions.csv \
  -u http://localhost:55001/reports/api/v1 -a <user:password>
