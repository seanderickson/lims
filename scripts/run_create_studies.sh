#!/usr/bin/env bash

##
# Re-generate the in-silico studies (o2):
# - check
# Prerequisites
# - Server is configured: Python 2.7.x, virtualenv
#
##

if [[ $# -eq 1 ]]; then
  OPTION=$1
  if [[ "$OPTION" != "real_run" ]]; then
    echo "Usage: $0 [real_run] (optional, if not set, do a trial run, with no email and no db changes)"
    exit 0
  fi
fi


### MODIFY AS NEEDED ###

credential_file=${credential_file:-"$SUPPORTDIR/production_data/sde_credentials.txt"}
BOOTSTRAP_PORT=53999

########################

REALPATH=${REALPATH:-"$(which realpath 2>/dev/null)"}
if [[ -z $REALPATH ]]; then
  PROG=$( basename $0 )
  error 'cannot find realpath'
fi

DEBUG=true

SCRIPTPATH="$($REALPATH $0)"
SCRIPTPATH="$(dirname $SCRIPTPATH)"
# strip the "scripts" directory to get the real BASEDIR
BASEDIR=${BASEDIR:-"$(dirname $SCRIPTPATH)"}
SUPPORTDIR=${SUPPORTDIR-"$(dirname $BASEDIR)"}
LOGFILE=${LOGFILE:-${BASEDIR:+"$BASEDIR/run_expirations.log"}}

if $DEBUG; then
  echo "set logfile to tty"
  LOGFILE=$(tty)
fi

VENV=${VENV:-${SUPPORTDIR:+"$SUPPORTDIR/virtualenv_o2"}}
if [[ -z $VENV ]]; then
  error 'no virtualenv available'
fi

source ${SCRIPTPATH}/utils.sh

DJANGO_CMD=$BASEDIR/manage.py
export DJANGO_SETTINGS_MODULE=lims.settings-server-commandline

cd $BASEDIR

function test_setup {

  maybe_activate_virtualenv
  
#  $DJANGO_CMD showmigrations db  
  
  echo "REAL RUN: $REAL_RUN"
  TEST_RUN_SETTINGS=" -admin_email_only -test_only " 
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"
  if [[ $REAL_RUN -eq 1 ]]; then
    TEST_RUN_SETTINGS=""
  fi
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"
  
  

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

  study_id=200001
  study_file=docs/studies/study_${study_id}.json
  echo "create study $study_id, $study_file ..." >>"$LOGFILE" 
  # lead_screener=sde_edit
  PYTHONPATH=${BASEDIR} python reports/utils/django_requests.py -c ${credential_file} \
    -a POST http://localhost:${BOOTSTRAP_PORT}/db/api/v1/study/create_screened_count_study \
    --header "Content-type: application/json" \
    --header "HTTP-Accept: application/json" \
    -f ${study_file} >>"$LOGFILE" 2>&1

  study_id=200002
  study_file=docs/studies/study_${study_id}.json
  echo "create study $study_id, $study_file ..." >>"$LOGFILE" 
  # lead_screener=sde_edit
  PYTHONPATH=${BASEDIR} python reports/utils/django_requests.py -c ${credential_file} \
    -a POST http://localhost:${BOOTSTRAP_PORT}/db/api/v1/study/create_screened_count_study \
    --header "Content-type: application/json" \
    --header "HTTP-Accept: application/json" \
    -f ${study_file} >>"$LOGFILE" 2>&1

  study_id=200003
  study_file=docs/studies/study_${study_id}.json
  echo "create study $study_id, $study_file ..." >>"$LOGFILE" 
  # lead_screener=sde_edit
  PYTHONPATH=${BASEDIR} python reports/utils/django_requests.py -c ${credential_file} \
    -a POST http://localhost:${BOOTSTRAP_PORT}/db/api/v1/study/create_confirmed_positive_study \
    --header "Content-type: application/json" \
    --header "HTTP-Accept: application/json" \
    -f ${study_file} >>"$LOGFILE" 2>&1

  echo "ping the studies to test..." >>"$LOGFILE" 
  study_id=200001
  PYTHONPATH=${BASEDIR} python reports/utils/django_requests.py -c ${credential_file} \
    -a GET http://localhost:${BOOTSTRAP_PORT}/db/api/v1/screenresult/${study_id}?limit=25 \
    --header "HTTP-Accept: application/json" \
    | mail -s "Study data ${study_id}" sean.erickson.hms@gmail.com 
  study_id=200002
  PYTHONPATH=${BASEDIR} python reports/utils/django_requests.py -c ${credential_file} \
    -a GET http://localhost:${BOOTSTRAP_PORT}/db/api/v1/screenresult/${study_id}?limit=25 \
    --header "HTTP-Accept: application/json" \
    | mail -s "Study data ${study_id}" sean.erickson.hms@gmail.com
  study_id=200003
  PYTHONPATH=${BASEDIR} python reports/utils/django_requests.py -c ${credential_file} \
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

if [[ $OPTION == "real_run" ]]; then
  REAL_RUN=1
fi

maybe_activate_virtualenv

create_studies

#test_setup
  

