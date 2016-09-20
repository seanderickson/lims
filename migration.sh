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
LOGFILE=${LOGFILE:-${BASEDIR:+"$BASEDIR/migration.log"}}

if [[ $# -lt 2 ]]
then
  echo "Usage: $0 [branch] [repository] [migration_file]"
  exit $WRONG_ARGS
fi 

BRANCH=$1
REMOTE=$2
#MIGRATION_PROPERTIES_FILE=${3:-${MIGRATION_PROPERTIES_FILE:-${BASE_DIR:+"$BASE_DIR/migration.properties"}}}
MIGRATION_PROPERTIES_FILE=${3:-${MIGRATION_PROPERTIES_FILE:-"./migration.properties"}}

source $MIGRATION_PROPERTIES_FILE

DBUSER=${DBUSER:-"screensaver-lims"}  
DB=${DB:-"screensaver-lims"}  
DBHOST=${DBHOST:-"localhost"}  
DBPASSWORD=${DBPASSWORD:-""}
SETENV_SCRIPT=${SETENV_SCRIPT:-""}
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

if [[ -n "$SETENV_SCRIPT" ]]; then
  AUTH_FILE=${AUTH_FILE:-"/opt/apache/conf/auth/dev.screensaver2.med.harvard.edu" }

  DJANGO_CMD="$SETENV_SCRIPT $AUTH_FILE ./manage.py"
fi


source ./utils.sh

function ts {
  date +%Y%m%dT%H%M%S%z
}

function _debug {
  if $DEBUG; then echo "$@"; fi
}

function gitpull {
  # git fetch --all
  # git pull
  # git merge seanderickson/master
  
  _debug -n "pulling branch $BRANCH from $REMOTE... "
  #  cp -a $SCRIPTPATH $SAVEPATH
  git fetch $REMOTE >>"$LOGFILE" 2>&1 || error "git-fetch failed: $?"
  git checkout $BRANCH >>"$LOGFILE" 2>&1 || error "git-checkout failed: $?"
  
  # reset the api_init_actions file, modified during migration
  git checkout ./db/static/api_init/api_init_actions.csv

  #  git checkout -- $SCRIPTPATH >>"$LOGFILE" 2>&1 || error "git-checkout $SCRIPTPATH failed: $?"

  if [[ -e db/static/api_init/vocabulary_data_generated.csv ]]; then
    mv db/static/api_init/vocabulary_data_generated.csv db/static/api_init/vocabulary_data_generated.old.csv
  fi
  
  mkdir bak
  mv db/static/api_init/vocabulary_* bak  
  
  git pull --ff-only $REMOTE $BRANCH >>"$LOGFILE" 2>&1 || error "git-pull failed: $?"

  #  mv -f $SAVEPATH $SCRIPTPATH
  #  update_deploy_info "$STATUS"

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

    pg_restore -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
      `ls -1 ${D}/screensaver*${filespec}.attached_file.pg_dump`  >>"$LOGFILE" 2>&1 

    echo "+++ LOADING result data: ${D}/screensaver*${filespec}.result_data.pg_dump $(ts)" >> "$LOGFILE"
    pg_restore -Fc --no-owner -h $DBHOST -d $DB -U $DBUSER \
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
  
  # migrate replaces syncdb; 
  # note; some of the system migrations (auth, sites, admin, contenttypes ) will
  # generate errors due to interdependencies; these appear to resolve themselves - 20151030
  
  for x in sites auth contenttypes admin sessions tastypie reports;
  do
    echo "migrate app: $x ..." >> "$LOGFILE"
    $DJANGO_CMD migrate $x; #  --fake-initial; 
  done
  
  echo "- Create the adminuser: $adminuser"
  # update, as of Django >1.7, initial_data is no longer used to initialize superuser
  # try this method instead
  
  # FIXME: echo password vulnerability
  
  _adminuser="'"$adminuser"'"
  _adminpass="'"$adminpass"'"
  _adminemail="'"$adminemail"'"
  echo "from django.contrib.auth.models import User; User.objects.create_superuser($_adminuser, $_adminemail , $_adminpass)" | $DJANGO_CMD shell

  # TODO: remove if memcached is installed
  $DJANGO_CMD createcachetable
}

function premigratedb {
  echo "pre migrations: $(ts) ..." >> "$LOGFILE"

  # these migrations can be run before bootstrapping
  
  # check which migrations are completed 
  # - if we've skipped db restore, then only apply latest
  completed_migrations=$($DJANGO_CMD migrate db --list | grep '[X]' | awk '{print $2}')
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
    
  echo "pre migrations completed: $(ts) " >> "$LOGFILE"
 
}

function migratedb {
  echo "running migrations: $(ts) ..." >> "$LOGFILE"

  # check which migrations are completed - if we've skipped db restore, then only apply latest
  completed_migrations=$($DJANGO_CMD migrate db --list | grep '\[X\]' | awk '{print $2}')
  echo "completed migrations: $completed_migrations" >> "$LOGFILE"
  
  migration='0004'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0005'
  if [[ ! $completed_migrations =~ $migration ]]; then
    echo "migration $migration: $(ts) ..." >> "$LOGFILE"
    $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
    echo "migration $migration complete: $(ts)" >> "$LOGFILE"
  fi
  
  migration='0006'
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
  if [[ $DB_FULL_MIGRATION -eq 1 ]]; then
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
    
    migration='0098' 
    if [[ ! $completed_migrations =~ $migration ]]; then
      echo "migration $migration: $(ts) ..." >> "$LOGFILE"
      $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
      echo "migration $migration complete: $(ts)" >> "$LOGFILE"
    fi

    # substance ID generation left until last
    migration='0099' 
    if [[ ! $completed_migrations =~ $migration ]]; then
      echo "migration $migration: $(ts) ..." >> "$LOGFILE"
      $DJANGO_CMD migrate db $migration >>"$LOGFILE" 2>&1 || error "db $migration failed: $?"
      echo "migration $migration complete: $(ts)" >> "$LOGFILE"
    fi
  fi
    
  echo "migrations completed: $(ts) " >> "$LOGFILE"

}

function bootstrap {
  echo "Bootstrapping the web server: $(ts) ...">> "$LOGFILE"
  
  BOOTSTRAP_PORT=${BOOTSTRAP_PORT:-55999}
  
  echo "run a local dev server on port $BOOTSTRAP_PORT..."
  nohup $DJANGO_CMD runserver --settings=lims.migration-settings --nothreading --noreload $BOOTSTRAP_PORT  &
  server_pid=$!
  if [[ "$?" -ne 0 ]]; then
    runserver_status =$?
    echo "bootstrap error, dev runserver status: $runserver_status"
    exit $runserver_status
  fi
#  echo "wait for server process: ($!) to start..."
#  wait $server_pid
  sleep 3
  
  echo "bootstrap the metahash data..."
  PYTHONPATH=. python reports/utils/db_init.py  \
    --input_dir=./reports/static/api_init/ \
    -f ./reports/static/api_init/api_init_actions.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/reports/api/v1 -U sde -p ${adminpass} >>"$LOGFILE" 2>&1 
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "bootstrap reports failed: $?"
  fi
    
  PYTHONPATH=. python reports/utils/db_init.py  \
    --input_dir=./db/static/api_init/ \
    -f ./db/static/api_init/api_init_actions.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/reports/api/v1 -U sde -p ${adminpass} >>"$LOGFILE" 2>&1
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "bootstrap db failed: $?"
  fi
  
  echo "bootstrap the server production data..."
  
  PYTHONPATH=. python reports/utils/db_init.py  \
    --input_dir=$BOOTSTRAP_PRODUCTION_DIR \
    -f ${BOOTSTRAP_PRODUCTION_DIR}/api_init_actions.csv \
    -u http://localhost:${BOOTSTRAP_PORT}/reports/api/v1 -U sde -p ${adminpass} >>"$LOGFILE" 2>&1
  if [[ $? -ne 0 ]]; then
    kill $server_pid
    error "bootstrap production data failed: $?"
  fi

  # add user "sde" to the screensaver_users table 20150831
  curl -v  --dump-header - -H "Content-Type: text/csv" --user sde:${adminpass} \
    -X PATCH http://localhost:${BOOTSTRAP_PORT}/db/api/v1/screensaveruser/ \
    --data-binary @${BOOTSTRAP_PRODUCTION_DIR}/screensaver_users-db-prod.csv

  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill $final_server_pid"
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
  
  npm run build 2>&1 || error "npm run build failed: $?"
  
  # TODO: frontend tests
  
  # add bootstrap to external folder - needed for the login template 
  mkdir css/external
  cp node_modules/bootstrap/dist/css/bootstrap.min.css css/external
  
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
    
  echo "frontend_setup done: $(ts)" >> "$LOGFILE"

}

function main {
  # Steps:
  
  gitpull
  
  restoredb
  
  restoredb_data
    
  maybe_activate_virtualenv
  
  pip install -r requirements.txt >>"$LOGFILE" 2>&1

  if [[ -n "$SETTINGS_FILE" ]]; then
    cp $SETTINGS_FILE lims/settings.py
  else
    echo "no SETTINGS_FILE migration.properties setting"
    # cp lims/settings-dist.py lims/settings.py
  fi
  
  mkdir logs
  
  if [[ $RUN_DB_TESTS -ne 0 ]] ; then
    $DJANGO_CMD test --verbosity=2 --settings=lims.settings_testing  \
      >> "$LOGFILE" 2>&1 || error "django tests failed: $?"
  fi
  
  frontend_setup
  
  django_syncdb

  premigratedb  

  bootstrap
  
  # the later migrations require the bootstrapped data
  migratedb
  
  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    tail -400 migration.log | mail -s "Migration completed $(ts)" sean.erickson.hms@gmail.com
  fi  
  # put this here to see if LSF will start reporting results
  exit 0
    
  # Integration test: run some grunt-mocha chai-jquery to test the UI?  
  
  # get the largest dataset, prime the well query/mutual positives 
  
  # wget https://dev.screensaver2.med.harvard.edu/db/api/v1/screenresult/1158?page=1&limit=25&offset=0&library_well_type__eq=experimental

  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    PYTHONPATH=. reports/utils/django_requests.py -u sde  \
      -a GET "https://dev.screensaver2.med.harvard.edu/db/api/v1/screenresult/1158?page=1&limit=25&offset=0&library_well_type__eq=experimental"
  fi

}

function code_bootstrap {
  # performs a full update, minus the database restoration
  
  gitpull
  
#  restoredb
  
#  restoredb_data
    
  maybe_activate_virtualenv
  
  pip install -r requirements.txt >>"$LOGFILE" 2>&1

  if [[ -n "$SETTINGS_FILE" ]]; then
    cp $SETTINGS_FILE lims/settings.py
  else
    echo "no SETTINGS_FILE migration.properties setting"
    # cp lims/settings-dist.py lims/settings.py
  fi
  
  mkdir logs
  
  if [[ $RUN_DB_TESTS -ne 0 ]] ; then
    $DJANGO_CMD test --verbosity=2 --settings=lims.settings_testing  \
      >> "$LOGFILE" 2>&1 || error "django tests failed: $?"
  fi
  
  frontend_setup
  
  django_syncdb

  premigratedb  

  bootstrap
  
  # the later migrations require the bootstrapped data
  migratedb
  
  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    if [[ $WARNINGS -ne '' ]]; then
      msg = "${WARNINGS} \n `tail -400 migration.log`"
      echo $msg | mail -s "Migration completed $(ts), with warning" sean.erickson.hms@gmail.com
    else
      tail -400 migration.log | mail -s "Migration completed $(ts)" sean.erickson.hms@gmail.com
    fi
  fi  
  # put this here to see if LSF will start reporting results
  exit 0

  
}  

echo "start migration: $(ts) ..."

main "$@"

#  restoredb
  
#  restoredb_data
    
#  maybe_activate_virtualenv
  
  
#  django_syncdb

#  premigratedb  

#  bootstrap
  
  # the later migrations require the bootstrapped data
#  migratedb







echo "migration finished: $(ts)"


function frontend_setup_grunt {
  echo "frontend_setup: $(ts) ..." >> "$LOGFILE"

  cd reports/static >>"$LOGFILE" 2>&1
  
  npm --version 2>&1 || error "npm not found: $?"
  
  rm -rf ./node_modules # untested 20150923
  npm install >>"$LOGFILE" 2>&1 || error "npm install failed: $?"
  npm install bower
  ./node_modules/.bin/bower cache clean
  rm -rf ./bower_components # untested
  ./node_modules/.bin/grunt bowercopy >>"$LOGFILE" 2>&1 || error "grunt bowercopy failed: $?"
  
  ./node_modules/.bin/grunt test >>"$LOGFILE" 2>&1 || warn "grunt test failed, see logfile: $?"
  
  cd ../..
  
  if [[ $IS_DEV_SERVER -ne 1 ]]; then
    $DJANGO_CMD collectstatic --noinput --ignore="*node_modules*" \
        --ignore="*bower_components*" --ignore="*api_init*" || error "collectstatic failed: $?"
  fi

  if [ -e ../wsgi/app.wsgi ]; then
    touch ../wsgi/app.wsgi
  fi
    
  echo "frontend_setup done: $(ts)" >> "$LOGFILE"
  
}
