#!/usr/bin/env bash

##
# Build a demo server
#
# Prerequisites
# - Python 2.7.x, pip, virtualenv
# - Node, NPM
# - git
# - postgres > 9.4
##

function error () {
  warn "$@"
  exit 1
}

WARNINGS=''
function warn () {
  WARNINGS = "${WARNINGS}, $@"
}

REALPATH=${REALPATH:-"$(which realpath 2>/dev/null)"}
if [[ -z $REALPATH ]]; then
  PROG=$( basename $0 )
  error 'cannot find realpath'
fi

SCRIPTPATH="$($REALPATH $0)"
SAVEPATH="$SCRIPTPATH.save"
BASEDIR=${BASEDIR:-"$(dirname $SCRIPTPATH)"}
SUPPORTDIR=${SUPPORTDIR-"$(dirname $BASEDIR)"}
LOGFILE=${LOGFILE:-${BASEDIR:+"$BASEDIR/create_demo.log"}}

if [[ $# -lt 2 ]]
then
  echo "Usage: $0 [branch] [repository] [properties_file]"
  exit $WRONG_ARGS
fi 

BRANCH=$1
REMOTE=$2
MIGRATION_PROPERTIES_FILE=${3:-${MIGRATION_PROPERTIES_FILE:-"./create_demo.properties"}}

source $MIGRATION_PROPERTIES_FILE

DBUSER=${DBUSER:-"demoscreensaver"}  
DB=${DB:-"demoscreensaver"}  
DBHOST=${DBHOST:-"localhost"}  
DBPASSWORD=${DBPASSWORD:-""}
#SETENV_SCRIPT=${SETENV_SCRIPT:-""}
DEBUG=${DEBUG:-false}
RUN_DB_TESTS=${RUN_DB_TESTS:-false}

# PATH TO node, npm, leave blank if part of the env path already
NODE_PATH=${NODE_PATH:-""}
if [[ -n "$NODE_PATH" ]]; then
  export PATH=${PATH}:$NODE_PATH
fi
GIT_PATH=${GIT_PATH:-""}
if [[ -n "$GIT_PATH" ]]; then
  export PATH=$GIT_PATH:${PATH}
fi

if $DEBUG; then
  echo "set logfile to tty"
  LOGFILE=$(tty)
fi

VENV=${VENV:-${SUPPORTDIR:+"$SUPPORTDIR/virtualenv"}}
if [[ -z $VENV ]]; then
  error 'no virtualenv available'
fi

DJANGO_CMD=./manage.py


source ./utils.sh

#if [[ -n "$SETENV_SCRIPT" ]]; then
#  AUTH_FILE=${AUTH_FILE:-"/opt/apache/conf/auth/dev.screensaver2.med.harvard.edu" }
#
#  DJANGO_CMD="$SETENV_SCRIPT $AUTH_FILE ./manage.py"
#fi

if [[ -n $AUTH_FILE ]]; then
  read_auth_file $AUTH_FILE
fi

# set an overlay settings.py (overrides logging settings, so we don't 
# overwrite/or change the perms for the server logs)
export DJANGO_SETTINGS_MODULE=lims.settings-server-commandline

DJANGO_CMD="./manage.py"

function ts {
  date +%Y%m%dT%H%M%S%z
}

function _debug {
  if $DEBUG; then echo "$@"; fi
}

function set_settings {

  echo "set_settings: $(ts) ..." >> "$LOGFILE"

  mkdir logs

  if [[ ! -e lims/settings.py ]]; then
    error "lims.settings.py must be created for the demo server; for example, see lims/settings-demoserver.py"
  fi
  
  if [[ ! -n "$APP_DATA_FILE" ]]; then
    error "APP_DATA_FILE property must be set"
  fi
  if [[ -e "$APP_DATA_FILE" ]]; then
    cp $APP_DATA_FILE lims/
  else
    error "no $APP_DATA_FILE found for APP_DATA_FILE migration.properties setting"
  fi

#  login_template="reports/templates/login.html"
#  demo_login_template="reports/templates/login_demo.html"
#  if [[ -e $demo_login_template ]]; then
#    echo "copy $demo_login_template to $login_template"
#    cp $demo_login_template $login_template
#  fi

  echo "set_settings done: $(ts) ..." >> "$LOGFILE"
}

function gitpull {
  
  _debug -n "pulling branch $BRANCH from $REMOTE... "
  git fetch $REMOTE >>"$LOGFILE" 2>&1 || error "git-fetch failed: $?"
  git checkout $BRANCH >>"$LOGFILE" 2>&1 || error "git-checkout failed: $?"

  # Forgo the heavy-handed approach  
  #  git reset --hard $REMOTE/$BRANCH >> "$LOGFILE" 2>&1 || error "git hard reset failed: $?"
  git checkout $REMOTE/$BRANCH ./db/static/api_init/*.csv
  git pull --ff-only $REMOTE $BRANCH >>"$LOGFILE" 2>&1 || error "git-pull failed: $?"

  _debug 'done'
  return 0
}

function create_demo_db {

  echo "create_demo_db $(ts) ..." >> "$LOGFILE"

  # drop DB if exists
  if [[ -x "$DROP_DB_COMMAND" ]]; then
    echo "execute DROP_DB_COMMAND: $DROP_DB_COMMAND at $(ts) ..." >> "$LOGFILE"
    $DROP_DB_COMMAND $DB $DBUSER >>"$LOGFILE" 2>&1 || error "dropdb failed: $?"
  else
    # test if the db exists
    psql -h $DBHOST -U $DBUSER -lqt | cut -d \| -f 1 | grep -w $DB
    if [[ "$?" -eq 0 ]]; then
      echo "execute dropdb $(ts) ..." >> "$LOGFILE"
      dropdb -U $DBUSER $DB -h $DBHOST >>"$LOGFILE" 2>&1 || error "dropdb failed: $?"
    fi
  fi
  
  if [[ $CREATE_DB -ne 0 ]]; then
    echo "execute createdb $(ts) ..." >> "$LOGFILE"
    createdb -U $DBUSER $DB -h $DBHOST --encoding=UTF8 >>"$LOGFILE" 2>&1 || error "createdb fails with status $?" 
  fi

  echo "create_demo_db done: $(ts) ..." >> "$LOGFILE"
}


function django_migrate {

  echo "django_migrate: $(ts) ..." >> "$LOGFILE"
  
  # Migrate all of the apps used by the LIMS: 
  # note; some of the system migrations (auth, sites, admin, contenttypes ) will
  # generate errors due to interdependencies; these appear to resolve themselves - 20151030
  
  # back up the migrations directory (currently only used for migrating existing databases)  
  if [[ ! -e db/migrations_new_install ]]; then
    error "demo server creation requires db/migrations_new_install directory"
  else
    echo "back up the existing migrations directory..."
    mv db/migrations db/migrations_bak
    mv db/migrations_new_install db/migrations
  fi
  
  
  for x in sites auth contenttypes admin sessions reports db;
  do
    echo "migrate app: $x ..." >> "$LOGFILE"
    $DJANGO_CMD migrate $x; #  --fake-initial; 
    if [[ $? -ne 0 ]]; then
      error "migrate app $x failed: $?" >> "$LOGFILE"
    fi
    $DJANGO_CMD showmigrations $x; #  --fake-initial; 
  done
  
  echo "- Create the superuser: $adminuser" >> "$LOGFILE"
  # update, as of Django >1.7, initial_data is no longer used to initialize superuser
  # try this method instead
  
  _username=`awk -F ':' '{print $1}' $credential_file`
  _username="'"$_username"'"
  _userpass=`awk -F ':' '{print $2}' $credential_file`
  _userpass="'"$_userpass"'"
  _useremail="'"$adminemail"'"
  python -c "import django; django.setup(); \
   from django.contrib.auth import get_user_model; \
   get_user_model()._default_manager.db_manager().create_superuser( \
   username=$_username, \
   email=$_useremail, \
   password=$_userpass)"  
  if [[ $? -ne 0 ]]; then
    error "create superuser failed: $?" >> "$LOGFILE"
  fi
  
  echo "- Create the demoadmin..." >> "$LOGFILE"
  _username="'demoadmin'"
  _userpass="'demoadmin'"
  _useremail="'demoadmin@iccbl.hms.harvard.edu'"
  python -c "import django; django.setup(); \
   from django.contrib.auth import get_user_model; \
   get_user_model()._default_manager.db_manager().create_superuser( \
   username=$_username, \
   email=$_useremail, \
   password=$_userpass)"  
  if [[ $? -ne 0 ]]; then
    error "create demoadmin failed: $?" >> "$LOGFILE"
  fi
  
  echo "- Create the demolabhead..." >> "$LOGFILE"
  _username="'demolabhead'"
  _userpass="'demolabhead'"
  _useremail="'demolabhead@iccbl.hms.harvard.edu'"
  python -c "import django; django.setup(); \
   from django.contrib.auth import get_user_model; \
   get_user_model()._default_manager.db_manager().create_user( \
   username=$_username, \
   email=$_useremail, \
   password=$_userpass, is_staff=True)"  
  if [[ $? -ne 0 ]]; then
    error "create demolabhead failed: $?" >> "$LOGFILE"
  fi
  
  echo "- Create the demoscreener1..." >> "$LOGFILE"
  _username="'demoscreener1'"
  _userpass="'demoscreener1'"
  _useremail="'demoscreener1@iccbl.hms.harvard.edu'"
  python -c "import django; django.setup(); \
   from django.contrib.auth import get_user_model; \
   get_user_model()._default_manager.db_manager().create_user( \
   username=$_username, \
   email=$_useremail, \
   password=$_userpass)"  
  
  echo "- Create the demoscreener2..." >> "$LOGFILE"
  _username="'demoscreener2'"
  _userpass="'demoscreener2'"
  _useremail="'demoscreener2@iccbl.hms.harvard.edu'"
  python -c "import django; django.setup(); \
   from django.contrib.auth import get_user_model; \
   get_user_model()._default_manager.db_manager().create_user( \
   username=$_username, \
   email=$_useremail, \
   password=$_userpass)"  
  
  # TODO: remove if memcached is installed
  $DJANGO_CMD createcachetable

  #  echo "restore migrations_existing_database and migrations_new_install directories..."
  #  mv db/migrations db/migrations_new_install
  #  mv db/migrations_existing_database db/migrations
  if [[ -e db/migrations_bak ]]; then
    echo "restore migrations directory..."
    mv db/migrations db/migrations_new_install
    mv db/migrations_bak db/migrations
  fi
  
  
  echo "django_migrate done: $(ts) ..." >> "$LOGFILE"
}

function load_demo_data {

  echo "load_demo_data $(ts) ..." >> "$LOGFILE"

  echo "run a local dev server on port $BOOTSTRAP_PORT..." >> "$LOGFILE"
  
  nohup $DJANGO_CMD runserver --settings=lims.migration-settings --nothreading \
  --noreload $BOOTSTRAP_PORT  &
  server_pid=$!
  if [[ "$?" -ne 0 ]]; then
    runserver_status =$?
    echo "bootstrap error, dev runserver status: $runserver_status" >> "$LOGFILE"
    exit $runserver_status
  fi

  echo "pause for server to start: $server_pid" >> "$LOGFILE"
  sleep 30

  echo "load_demo_data: $(ts) ...">> "$LOGFILE"
  
  echo "load demo reports data..." >> "$LOGFILE"
  PYTHONPATH=. python reports/utils/db_init.py \
    --input_dir=./db/static/demo_data/ \
    -f ./db/static/demo_data/demo_init_actions_reports.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/reports/api/v1 -c ${credential_file} >>"$LOGFILE" 2>&1 
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "load demo reports data failed: $?" >> "$LOGFILE"
  fi

  echo "load demo db data..." >> "$LOGFILE"
  PYTHONPATH=. python reports/utils/db_init.py  \
    --input_dir=./db/static/demo_data/ \
    -f ./db/static/demo_data/demo_init_actions.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/db/api/v1 -U demoadmin -p demoadmin >>"$LOGFILE" 2>&1 
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "load demo db data failed: $?" >> "$LOGFILE"
  fi
    
  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill final_server_pid: $final_server_pid" >> "$LOGFILE"
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  sleep 3

  echo "load_demo_data done: $(ts)" >> "$LOGFILE"

}

function bootstrap {
  echo "Bootstrapping the web server: $(ts) ...">> "$LOGFILE"
  
  BOOTSTRAP_PORT=${BOOTSTRAP_PORT:-55999}
  
  echo "run a local dev server on port $BOOTSTRAP_PORT..." >> "$LOGFILE"
  
  nohup $DJANGO_CMD runserver --settings=lims.migration-settings --nothreading \
  --noreload $BOOTSTRAP_PORT  &
  
  server_pid=$!
  if [[ "$?" -ne 0 ]]; then
    runserver_status =$?
    echo "bootstrap error, dev runserver status: $runserver_status" >> "$LOGFILE"
    exit $runserver_status
  fi
  echo "pause for server to start: $server_pid" >> "$LOGFILE"
  sleep 10
  
  echo "bootstrap the metahash data..." >> "$LOGFILE"
  PYTHONPATH=. python reports/utils/db_init.py  \
    --input_dir=./reports/static/api_init/ \
    -f ./reports/static/api_init/api_init_actions.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/reports/api/v1 -c ${credential_file} >>"$LOGFILE" 2>&1 
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "bootstrap reports failed: $?"
  fi
    
  PYTHONPATH=. python reports/utils/db_init.py  \
    --input_dir=./db/static/api_init/ \
    -f ./db/static/api_init/api_init_actions.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/reports/api/v1 -c ${credential_file} >>"$LOGFILE" 2>&1
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "bootstrap db failed: $?"
  fi

  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill final_server_pid: $final_server_pid" >> "$LOGFILE"
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  sleep 3
  echo "Bootstrapping the web server done: $(ts)" >> "$LOGFILE"
}

function frontend_setup {
  echo "frontend_setup: $(ts) ..." >> "$LOGFILE"

  cd reports/static >>"$LOGFILE" 2>&1
  npm --version 2>&1 || error "npm not found: $?"
  rm -rf ./node_modules 
  npm install >>"$LOGFILE" 2>&1 || error "npm install failed: $?"
  cd ../..
  
  echo "frontend_setup done: $(ts)" >> "$LOGFILE"

}

function frontend_deploy {

  echo "frontend_deploy: $(ts) ..." >> "$LOGFILE"

  cd reports/static >>"$LOGFILE" 2>&1
  rm bundle.*.js
  rm 1.bundle.*.js
  npm run build 2>&1 || error "npm run build failed: $?"
  
  # TODO: frontend tests
  
  cd ../..
  
  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    $DJANGO_CMD collectstatic --noinput --clear \
        --ignore="*node_modules*" \
        --ignore="*bower_components*" \
        --ignore="*test_data*" \
        --ignore="*.json" \
        --ignore="Gruntfile.js" \
        --ignore="*api_init*" \
        --ignore="js*" || error "collectstatic failed: $?"
    # FIXME: image and css files must be copied still
    # --ignore="css*"
    # --ignroe="images*"      
        
  fi
  
  if [ -e ../wsgi/app.wsgi ]; then
    touch ../wsgi/app.wsgi
  fi

  echo "frontend_deploy done: $(ts)" >> "$LOGFILE"

}


function main {
  # Steps:
  
  gitpull
  
  create_demo_db
  
  maybe_activate_virtualenv
  
  pip install -r requirements.txt >>"$LOGFILE" 2>&1

  set_settings

  django_migrate

  frontend_setup

  frontend_deploy
  
  bootstrap

  load_demo_data
  

}


echo "start create demo server script: $(ts) ..." >>"$LOGFILE"

 main "$@"

#  gitpull
  
#  create_demo_db
  
#  maybe_activate_virtualenv
  
#  pip install -r requirements.txt >>"$LOGFILE" 2>&1

#  set_settings

#  django_migrate

#  frontend_setup

#  frontend_deploy
  
#  bootstrap

#  load_demo_data
  


echo "demo server creation script finished: $(ts)" >>"$LOGFILE"

