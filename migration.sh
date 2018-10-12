#!/usr/bin/env bash

##
# HMS - ICCBL specific:
# 
# Migrate the Screensaver 1 database to the Screensaver 2 database
#
# Prerequisites
# - Python 2.7.x, pip, virtualenv
# - Node, NPM
# - git
# - postgres
##

if [[ $# -lt 2 ]]
then
  echo "Usage: $0 [branch] [repository] [migration_file]"
  exit $WRONG_ARGS
fi 

BRANCH=$1
REMOTE=$2
MIGRATION_PROPERTIES_FILE=${3:-${MIGRATION_PROPERTIES_FILE:-"./migration.properties"}}

source ./utils.sh

REALPATH=${REALPATH:-"$(which realpath 2>/dev/null)"}
if [[ -z $REALPATH ]]; then
  PROG=$( basename $0 )
  error 'cannot find realpath'
fi

SCRIPTPATH="$($REALPATH $0)"
SAVEPATH="$SCRIPTPATH.save"
BASEDIR=${BASEDIR:-"$(dirname $SCRIPTPATH)"}
SUPPORTDIR=${SUPPORTDIR-"$(dirname $BASEDIR)"}
LOGFILE=${LOGFILE:-${BASEDIR:+"$BASEDIR/migration.log"}}

source $MIGRATION_PROPERTIES_FILE

DBUSER=${DBUSER:-"screensaver-lims"}  
DB=${DB:-"screensaver-lims"}  
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

VENV=${VENV:-${SUPPORTDIR:+"$SUPPORTDIR/virtualenv_o2"}}
if [[ -z $VENV ]]; then
  error 'no virtualenv available'
fi

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

function gitpull {
  
  _debug -n "pulling branch $BRANCH from $REMOTE... "
  git fetch $REMOTE >>"$LOGFILE" 2>&1 || error "git-fetch failed: $?"
  git checkout $BRANCH >>"$LOGFILE" 2>&1 || error "git-checkout failed: $?"

  # reset api_init files that may have been modified during the previous migration
  if [[ -e db/static/api_init/vocabulary_data_generated.csv ]]; then
    # save the old version, might be useful
    mv db/static/api_init/vocabulary_data_generated.csv db/static/api_init/vocabulary_data_generated.old
  fi
  # Forgo the heavy-handed approach  
  #  git reset --hard $REMOTE/$BRANCH >> "$LOGFILE" 2>&1 || error "git hard reset failed: $?"
  git checkout $REMOTE/$BRANCH ./db/static/api_init/*.csv
  git pull --ff-only $REMOTE $BRANCH >>"$LOGFILE" 2>&1 || error "git-pull failed: $?"

  _debug 'done'
  return 0
}

function restoredb {
  echo "restoredb: DB_LOAD_SCHEMA_ONLY:$DB_LOAD_SCHEMA_ONLY  $(ts) ..." >> "$LOGFILE"

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
  
  D=${PG_RESTORE_DIR:-.}
  filespec=${PG_RESTORE_FILESPEC:-''}
  
  # NOTE: restore generates some errors - manual verification is required.

  if [[ $DB_LOAD_SCHEMA_ONLY -eq 0 ]]; then
    echo "execute pg_restore: ${D}/screensaver*${filespec}.excl_big_data.pg_dump - $(ts) ..." >> "$LOGFILE"
    pg_restore -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
     `ls -1 ${D}/screensaver*${filespec}.excl_big_data.pg_dump` >>"$LOGFILE" 2>&1 
  else
    # For quick testing
    echo "execute pg_restore: ${D}/screensaver*${filespec}.schema.pg_dump - $(ts) ..." >> "$LOGFILE"
    pg_restore -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
      `ls -1 ${D}/screensaver*${filespec}.schema.pg_dump` >>"$LOGFILE" 2>&1 
  fi

  # clear out the job directory on db rebuild
  rm -rf ${SUPPORTDIR}/logs/background
  
  echo "restoredb completed: $(ts) " >> "$LOGFILE"

}

function restoredb_data {
  echo "restoredb_data: DB_LOAD_SCHEMA_ONLY: $DB_LOAD_SCHEMA_ONLY $(ts) ..." >> "$LOGFILE"
  echo "restoredb_data: DB_SKIP_BIG_FILES: $DB_SKIP_BIG_FILES $(ts) ..." >> "$LOGFILE"

  D=${PG_RESTORE_DIR:-.}
  filespec=${PG_RESTORE_FILESPEC:-''}

  echo "'pg_restore' ${D}/screensaver*${filespec}.attached_file_schema.pg_dump $(ts) ..." >> "$LOGFILE"
  pg_restore -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
    `ls -1 ${D}/screensaver*${filespec}.attached_file_schema.pg_dump`  >>"$LOGFILE" 2>&1 
  echo "'pg_restore' ${D}/screensaver*${filespec}.result_data_schema.pg_dump ..." >> "$LOGFILE"
  pg_restore -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
    `ls -1 ${D}/screensaver*${filespec}.result_data_schema.pg_dump`  >>"$LOGFILE" 2>&1 

  if [[ ( $DB_LOAD_SCHEMA_ONLY -eq 0 && $DB_SKIP_BIG_FILES -eq 0 ) ]]; then
    echo "+++ LOADING attached file and result data ... "
    echo "+++ LOADING attached file data: ${D}/screensaver*${filespec}.attached_file.pg_dump $(ts)" >> "$LOGFILE"

    PG_RESTORE_CMD="pg_restore"
    if [[ $IS_DEV_SERVER -eq 1 ]]; then
      PG_RESTORE_CMD="nice -5 pg_restore"
    fi

    $PG_RESTORE_CMD -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
      `ls -1 ${D}/screensaver*${filespec}.attached_file.pg_dump`  >>"$LOGFILE" 2>&1 

    echo "+++ LOADING result data: ${D}/screensaver*${filespec}.result_data.pg_dump $(ts)" >> "$LOGFILE"
    $PG_RESTORE_CMD -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
      `ls -1 ${D}/screensaver*${filespec}.result_data.pg_dump`  >>"$LOGFILE" 2>&1 
  fi
  
  if [[ ( $DB_LOAD_TESTING_DATA -eq 1 ) ]]; then
    for x in `ls -1 ${TESTING_DATA_DIR}/result_value*`; do
      echo "importing: $x $(ts)"
      psql -h $DBHOST -d $DB -U $DBUSER -c "\copy result_value from $x"
    done          
    for x in `ls -1 ${TESTING_DATA_DIR}/assay_well*`; do
      echo "importing: $x $(ts)"
      psql -h $DBHOST -d $DB -U $DBUSER -c "\copy assay_well from $x"
    done          
  fi  
  
  echo "vacuum analyze $(ts)"
  psql -h $DBHOST -d $DB -U $DBUSER -c "vacuum analyze;"
  
  # TODO: create a check to validate db imports
  
  echo "restoredb_data completed: $(ts) " >> "$LOGFILE"
  
  return 0
}

function django_syncdb {

  echo "django_syncdb: $(ts) ..." >> "$LOGFILE"
  
  # migrate replaces syncdb; 
  # note; some of the system migrations (auth, sites, admin, contenttypes ) will
  # generate errors due to interdependencies; these appear to resolve themselves - 20151030
  
  for x in sites auth contenttypes admin sessions reports;
  do
    echo "migrate app: $x ..." >> "$LOGFILE"
    $DJANGO_CMD migrate $x; #  --fake-initial; 
  done
  
  echo "- Create the adminuser: $adminuser"
  # update, as of Django >1.7, initial_data is no longer used to initialize superuser
  # try this method instead
  _adminuser="'"$adminuser"'"
  _adminpass="'"$adminpass"'"
  _adminemail="'"$adminemail"'"
  python -c "import django; django.setup(); \
   from django.contrib.auth import get_user_model; \
   get_user_model()._default_manager.db_manager().create_superuser( \
   username=$_adminuser, \
   email=$_adminemail, \
   password=$_adminpass)"  
  

  # TODO: remove if memcached is installed
  $DJANGO_CMD createcachetable

  echo "django_syncdb done: $(ts) ..." >> "$LOGFILE"
}

function premigratedb {
  echo "pre migrations: $(ts) ..." >> "$LOGFILE"

  # these migrations can be run before bootstrapping
  
  # check which migrations are completed 
  # - if we've skipped db restore, then only apply latest
  completed_migrations=$($DJANGO_CMD showmigrations db | grep '[X]' | awk '{print $2}')
  echo "completed migrations: $completed_migrations" >> "$LOGFILE"
  
  migration='0001'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration --fake >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
  fi
  migration='0002'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"

    psql -U $DBUSER $DB -h $DBHOST -a -v ON_ERROR_STOP=1 \
      -f ./db/migrations/manual/0002_initial_django_prep.sql \
      >>"$LOGFILE" 2>&1 || error "manual script 0002 failed: $?"

  fi
    
  migration='0003'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"

    psql -U $DBUSER $DB -h $DBHOST -a -v ON_ERROR_STOP=1 \
      -f ./db/migrations/manual/0003_controlled_vocabularies.sql \
      >>"$LOGFILE" 2>&1 || error "manual script 0003 failed: $?"
  fi
  migration='0004'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
  fi
    
  echo "pre migrations completed: $(ts) " >> "$LOGFILE"
 
}

function migratedb {
  echo "running migrations: $(ts) ..." >> "$LOGFILE"

  # check which migrations are completed - if we've skipped db restore, then only apply latest
  completed_migrations=$($DJANGO_CMD showmigrations db | grep '\[X\]' | awk '{print $2}')
  echo "completed migrations: $completed_migrations" >> "$LOGFILE"
  
  migration='0007'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0008'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0010'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0011'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0012'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0013'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0014'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  migration='0015'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    psql -U $DBUSER $DB -h $DBHOST -a -v ON_ERROR_STOP=1 \
        -f ./db/migrations/manual/0015_post_migrate.sql >>"$LOGFILE" 2>&1 || error "manual script 0015 failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  migration='0016'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    psql -U $DBUSER $DB -h $DBHOST -a -v ON_ERROR_STOP=1 \
        -f ./db/migrations/manual/0016_create_copy_wells.sql >>"$LOGFILE" 2>&1 || error "manual script 0016 failed: $?"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0017'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0018'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0019'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0020'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  

  migration='0021'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
    
  migration='0022'
  if [[ ! $completed_migrations =~ $migration ]]; then
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
    

# TEMP: 20170614 disable post migrations; leaves vestigal fields/tables in place TODO: reinstate    
  migration='0098' 
  if [[ ! $completed_migrations =~ $migration ]]; then
   echo "migration $migration: $(ts) ..." >> "$LOGFILE"
   $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
   echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi

  # substance ID generation left until last
#    migration='0099' 
#    if [[ ! $completed_migrations =~ $migration ]]; then
#      echo "migration $migration: $(ts) ..." >> "$LOGFILE"
#      $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
#      echo "migration $migration complete: $(ts)" >> "$LOGFILE"
#    fi
    
  echo "migrations completed: $(ts) " >> "$LOGFILE"

}

function migrate_result_values {
  psql -U $DBUSER $DB -h $DBHOST -a -v ON_ERROR_STOP=1 \
    -f ./db/migrations/manual/result_value_migration.sql >>"$LOGFILE" 2>&1 \
    || error "manual migrate_result_values failed: $?"
}

function result_value_cleanup {
  echo "Result value cleanup: $(ts) ...">> "$LOGFILE"
  psql -U $DBUSER $DB -h $DBHOST -a \
    -f ./db/migrations/manual/result_value_cleanup.sql >>"$LOGFILE" 2>&1 \
    || error "manual migrate_result_cleanup failed: $?"
  echo "Result value cleanup completed: $(ts) " >> "$LOGFILE"
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
  echo "wait for server process: $server_pid to start..." >>"$LOGFILE" 
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

  echo "PATCH the db/static/api_init/labaffiliation_updates.csv file..."  >> "$LOGFILE"
  curl -v  --dump-header - -H "Content-Type: text/csv" --user sde:${adminpass} \
    -X PATCH http://localhost:${BOOTSTRAP_PORT}/db/api/v1/labaffiliation \
    --data-binary @db/static/api_init/labaffiliation_updates.csv >>"$LOGFILE" 2>&1
  
  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill $final_server_pid ..." >>"$LOGFILE" 
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  # kill $server_pid

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
  rm bundle.*
  rm 1.bundle.*
  rm 2.bundle.*
  rm 3.bundle.*
  
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

function setup_production_users {

  echo "setup_production_users: $(ts) ..." >> "$LOGFILE"

  echo "setup_production_users using a local dev server on port $BOOTSTRAP_PORT..." >> "$LOGFILE"
  
  # FIXME: using the local server may not be necessary as the server has been bootstrapped
  nohup $DJANGO_CMD runserver --settings=lims.migration-settings --nothreading \
  --noreload $BOOTSTRAP_PORT  &
  
  server_pid=$!
  if [[ "$?" -ne 0 ]]; then
    runserver_status =$?
    echo "setup test data error, dev runserver status: $runserver_status" >> "$LOGFILE"
    exit $runserver_status
  fi
  echo "wait for server process: $server_pid to start..." >>"$LOGFILE"  
  sleep 10
  echo "run api init ..." >>"$LOGFILE" 
  PYTHONPATH=. python reports/utils/db_init.py  \
    --input_dir=$BOOTSTRAP_PRODUCTION_DIR \
    -f ${BOOTSTRAP_PRODUCTION_DIR}/api_init_actions_patch.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/reports/api/v1 -c ${credential_file} >>"$LOGFILE" 2>&1
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "bootstrap production data failed: $?"
  fi

  echo "Add user 'sde' to the screensaver_users table..." >> "$LOGFILE"
  curl -v  --dump-header - -H "Content-Type: text/csv" --user sde:${adminpass} \
    -X PATCH http://localhost:${BOOTSTRAP_PORT}/db/api/v1/screensaveruser/ \
    --data-binary @${BOOTSTRAP_PRODUCTION_DIR}/screensaver_users-db-prod.csv  >>"$LOGFILE" 2>&1

  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill $final_server_pid"
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  # kill $server_pid

  echo "setup_production_users done: $(ts)" >> "$LOGFILE"

}

function run_expiration_scripts {
  # Run user and screen privacy expiration scripts
  echo "run user and screen privacy expiration scripts: $(ts) ..." >> "$LOGFILE"

  # TODO: move these to cron jobs on deployment
  
  echo "1.a Send SM user data privacy expiration notifications..." >>"$LOGFILE"
  PYTHONPATH=. python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type sm -days_to_expire 730 -days_ahead_to_notify 14 \
  -email_message_directory db/static/user_agreement/ \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_log_filename ../logs/mail_user_agreement_notification.log \
  -v -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.b Expire SM user agreements ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=. python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type sm -days_to_expire 730 \
  -email_message_directory db/static/user_agreement/ \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_log_filename ../logs/mail_user_agreement_expiration.log \
  -v -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.c Send RNAi user data privacy expiration notifications..." >>"$LOGFILE"
  PYTHONPATH=. python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type rnai -days_to_expire 730 -days_ahead_to_notify 14 \
  -email_message_directory db/static/user_agreement/ \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_log_filename ../logs/mail_user_agreement_notification.log \
  -v -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.d Expire RNAi user agreements ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=. python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type rnai -days_to_expire 730 \
  -email_message_directory db/static/user_agreement/ \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_log_filename ../logs/mail_user_agreement_expiration.log \
  -v -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.a Adjust screen Data Privacy Expiration Dates ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=. python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -adjust_expiration_days_from_activity 790 \
  -email_log_filename ../logs/mail_screen_dped_adjust.log \
  -v -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.b Notify of screen privacy expirations ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=. python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -days_ahead_to_notify 60 \
  -email_log_filename ../logs/mail_screen_dped_notification.log \
  -v -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.c Expire screen data sharing levels ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=. python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -expire \
  -email_log_filename ../logs/mail_screen_dped_expiration.log \
  -v -test_only -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.d Notify admins for screen publications ..." >>"$LOGFILE" 2>&1
  echo "Using pass file: ${credential_file}" >>"$LOGFILE"
  PYTHONPATH=. python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info 'Jen Smith (jennifer_smith@hms.harvard.edu)' \
  -admin_from_email screensaver-feedback@hms.harvard.edu \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -notifyofpublications \
  -email_log_filename ../logs/mail_screen_notifyofpublications.log \
  -v -test_only -admin_email_only >>"$LOGFILE" 2>&1


  echo "Done: user and screen privacy expiration scripts: $(ts)" >> "$LOGFILE"
}

function create_studies {

  # Create in_silico statistical studies
  echo "create in_silico statistical studies: $(ts) ..." >> "$LOGFILE"

  echo "create in_silico statistical studies using a local dev server on port $BOOTSTRAP_PORT..." >>"$LOGFILE" 2>&1
  
  # FIXME: using the local server may not be necessary as the server has been bootstrapped
   
  nohup $DJANGO_CMD runserver --settings=lims.migration-settings --nothreading \
  --noreload $BOOTSTRAP_PORT  &
  
  server_pid=$!
  if [[ "$?" -ne 0 ]]; then
    runserver_status =$?
    echo "setup test data error, dev runserver status: $runserver_status"
    exit $runserver_status
  fi
  
  echo "wait for server process: $server_pid to start..." >>"$LOGFILE" 
  sleep 10
  echo "create study $study_id, $study_file ..." >>"$LOGFILE" 
  study_id=200001
  study_file=docs/studies/study_${study_id}.json
  # lead_screener=sde_edit
  PYTHONPATH=. python reports/utils/django_requests.py -c ${credential_file} \
    -a POST http://localhost:${BOOTSTRAP_PORT}/db/api/v1/study/create_screened_count_study \
    --header "Content-type: application/json" \
    --header "HTTP-Accept: application/json" \
    -f ${study_file} >>"$LOGFILE" 2>&1

  study_id=200002
  study_file=docs/studies/study_${study_id}.json
  echo "create study $study_id, $study_file ..." >>"$LOGFILE" 
  # lead_screener=sde_edit
  PYTHONPATH=. python reports/utils/django_requests.py -c ${credential_file} \
    -a POST http://localhost:${BOOTSTRAP_PORT}/db/api/v1/study/create_screened_count_study \
    --header "Content-type: application/json" \
    --header "HTTP-Accept: application/json" \
    -f ${study_file} >>"$LOGFILE" 2>&1

  study_id=200003
  study_file=docs/studies/study_${study_id}.json
  echo "create study $study_id, $study_file ..." >>"$LOGFILE" 
  # lead_screener=sde_edit
  PYTHONPATH=. python reports/utils/django_requests.py -c ${credential_file} \
    -a POST http://localhost:${BOOTSTRAP_PORT}/db/api/v1/study/create_confirmed_positive_study \
    --header "Content-type: application/json" \
    --header "HTTP-Accept: application/json" \
    -f ${study_file} >>"$LOGFILE" 2>&1

  echo "ping the studies to test..." >>"$LOGFILE" 
  study_id=200001
  PYTHONPATH=. python reports/utils/django_requests.py -c ${credential_file} \
    -a GET http://localhost:${BOOTSTRAP_PORT}/db/api/v1/screenresult/${study_id}?limit=25 \
    --header "HTTP-Accept: application/json" \
    | mail -s "Study data ${study_id}" sean.erickson.hms@gmail.com 
  study_id=200002
  PYTHONPATH=. python reports/utils/django_requests.py -c ${credential_file} \
    -a GET http://localhost:${BOOTSTRAP_PORT}/db/api/v1/screenresult/${study_id}?limit=25 \
    --header "HTTP-Accept: application/json" \
    | mail -s "Study data ${study_id}" sean.erickson.hms@gmail.com
  study_id=200003
  PYTHONPATH=. python reports/utils/django_requests.py -c ${credential_file} \
    -a GET http://localhost:${BOOTSTRAP_PORT}/db/api/v1/screenresult/${study_id}?limit=25 \
    --header "HTTP-Accept: application/json" \
    | mail -s "Study data ${study_id}" sean.erickson.hms@gmail.com

  ####
  echo "create in_silico statistical studies finished, stop server ..." >>"$LOGFILE" 
  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill $final_server_pid" >>"$LOGFILE" 
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  # kill $server_pid

  echo "create in_silico statistical studies done: $(ts)" >> "$LOGFILE"

}


function main {
  # Steps:
  
  gitpull
  
  restoredb
  
  restoredb_data
  
  maybe_activate_virtualenv
  
  pip install -r requirements.txt >>"$LOGFILE" 2>&1

# NOTE: 20170726 - removed from the migration;
# - a data exposure ocurred due to the settings "DEBUG" flag being reset by 
# this action. The settings file should be manually updated in the future.
#  if [[ -n "$SETTINGS_FILE" ]]; then
#    cp $SETTINGS_FILE lims/settings.py
#  else
#    echo "no SETTINGS_FILE migration.properties setting"
#    # cp lims/settings-dist.py lims/settings.py
#  fi
  
  mkdir logs
  
  if [[ $RUN_DB_TESTS -ne 0 ]] ; then
    $DJANGO_CMD test --verbosity=2 --settings=lims.settings_testing  \
      >> "$LOGFILE" 2>&1 || error "django tests failed: $?"
  fi
  
  frontend_setup
  
  frontend_deploy  
  
  django_syncdb

  premigratedb

  bootstrap
  
  # the later migrations require the bootstrapped data
  migratedb
  
  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    tail -400 migration.log | mail -s "Migration completed $(ts)" sean.erickson.hms@gmail.com
  fi  
  
  setup_production_users
  
  create_studies

  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    run_expiration_scripts
  fi
  # put this here to see if LSF will start reporting results
  # exit 0
    
  # Integration test: run some grunt-mocha chai-jquery to test the UI?  
  
  # get the largest dataset, prime the well query/mutual positives 
  
  # wget https://dev.screensaver2.med.harvard.edu/db/api/v1/screenresult/1158?page=1&limit=25&offset=0&library_well_type__eq=experimental

  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    PYTHONPATH=. python reports/utils/django_requests.py  -c ${credential_file} \
      -a GET "https://dev.screensaver2.med.harvard.edu/db/api/v1/screenresult/1158?page=1&limit=25&offset=0&library_well_type__eq=experimental"
  fi
  
  result_value_cleanup
  
  if [[ $MIGRATE_RESULT_VALUE_TABLE -ne 0 ]]; then
    migrate_result_values
  fi

}


echo "start migration: $(ts) ..." >>"$LOGFILE"

main "$@"

#  gitpull
  
#  restoredb
  
#  restoredb_data
    
#  maybe_activate_virtualenv
  
#  pip install -r requirements.txt >>"$LOGFILE" 2>&1

#  frontend_setup

#  frontend_deploy
  
#  django_syncdb

#  premigratedb  

#  bootstrap
  
# the later migrations require the bootstrapped data
#  migratedb
  
#  setup_production_users
  
#  create_studies

# run_expiration_scripts

#  result_value_cleanup




echo "migration finished: $(ts)" >>"$LOGFILE"

