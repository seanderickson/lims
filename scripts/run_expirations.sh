#!/usr/bin/env bash

##
# Run the expiration scripts (o2):
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

DJANGO_CMD=${BASEDIR}/manage.py
export DJANGO_SETTINGS_MODULE=lims.settings-server-commandline

### MODIFY AS NEEDED ###

credential_file=${credential_file:-"$SUPPORTDIR/production_data/sde_credentials.txt"}
BOOTSTRAP_PORT=53999
MAIL_RECIPIENT_LIST='sean.erickson.hms@gmail.com'

# TODO: move these to settings.py

USER_DSL_DAYS_TO_EXPIRE=730
USER_DSL_DAYS_AHEAD_TO_NOTIFY=14
SCREEN_DSL_EXPIRE_DAYS_FROM_ACTIVITY=790
SCREEN_DSL_DAYS_AHEAD_TO_NOTIFY=60

ADMIN_FROM_EMAIL='screensaver-feedback@hms.harvard.edu'

CONTACT_PRINTED_INFO='Jen Smith (jennifer_smith@hms.harvard.edu)'

########################

cd $BASEDIR

function test_setup {

  echo "REAL RUN: $REAL_RUN"
  TEST_RUN_SETTINGS=" -admin_email_only -test_only " 
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"
  if [[ $REAL_RUN -eq 1 ]]; then
    TEST_RUN_SETTINGS=""
  fi
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"

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

  SERVER_URL=http://localhost:${BOOTSTRAP_PORT}

  echo "2.c Expire screen data sharing levels ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -mail_recipient_list ${MAIL_RECIPIENT_LIST} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -expire \
  -email_log_filename ../logs/mail_screen_dped_expiration.log \
  -v $TEST_RUN_SETTINGS -test_only -admin_email_only >>"$LOGFILE" 2>&1

  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill $final_server_pid ..." >>"$LOGFILE" 
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  # kill $server_pid

}


function run_expiration_scripts {
  # Run user and screen privacy expiration scripts
  echo "run user and screen privacy expiration scripts: $(ts) ..." >> "$LOGFILE"

  echo "REAL RUN: $REAL_RUN"
  TEST_RUN_SETTINGS=" -admin_email_only -test_only " 
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"
  if [[ $REAL_RUN -eq 1 ]]; then
    TEST_RUN_SETTINGS=""
  fi
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"

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

  SERVER_URL=http://localhost:${BOOTSTRAP_PORT}

  echo "1.a Send SM user data privacy expiration notifications..." >>"$LOGFILE"
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type sm -days_to_expire $USER_DSL_DAYS_TO_EXPIRE -days_ahead_to_notify $USER_DSL_DAYS_AHEAD_TO_NOTIFY \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_log_filename ../logs/mail_user_agreement_notification.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.b Expire SM user agreements ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type sm -days_to_expire $USER_DSL_DAYS_TO_EXPIRE \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -mail_recipient_list ${MAIL_RECIPIENT_LIST} \
  -email_log_filename ../logs/mail_user_agreement_expiration.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.c Send RNAi user data privacy expiration notifications..." >>"$LOGFILE"
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type rnai -days_to_expire $USER_DSL_DAYS_TO_EXPIRE -days_ahead_to_notify $USER_DSL_DAYS_AHEAD_TO_NOTIFY \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_log_filename ../logs/mail_user_agreement_notification.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.d Expire RNAi user agreements ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type rnai -days_to_expire $USER_DSL_DAYS_TO_EXPIRE \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -mail_recipient_list ${MAIL_RECIPIENT_LIST} \
  -email_log_filename ../logs/mail_user_agreement_expiration.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.a Adjust screen Data Privacy Expiration Dates ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -mail_recipient_list ${MAIL_RECIPIENT_LIST} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -adjust_expiration_days_from_activity $SCREEN_DSL_EXPIRE_DAYS_FROM_ACTIVITY \
  -email_log_filename ../logs/mail_screen_dped_adjust.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.b Notify of screen privacy expirations ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -days_ahead_to_notify $SCREEN_DSL_DAYS_AHEAD_TO_NOTIFY \
  -email_log_filename ../logs/mail_screen_dped_notification.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.c Expire screen data sharing levels ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -mail_recipient_list ${MAIL_RECIPIENT_LIST} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -expire \
  -email_log_filename ../logs/mail_screen_dped_expiration.log \
  -v $TEST_RUN_SETTINGS -test_only -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.d Notify admins for screen publications ..." >>"$LOGFILE" 2>&1
  echo "Using pass file: ${credential_file}" >>"$LOGFILE"
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -mail_recipient_list ${MAIL_RECIPIENT_LIST} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -notifyofpublications \
  -email_log_filename ../logs/mail_screen_notifyofpublications.log \
  -v $TEST_RUN_SETTINGS -test_only -admin_email_only >>"$LOGFILE" 2>&1

  ####
  echo "run_expirations finished, stop server ..." >>"$LOGFILE" 
  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill $final_server_pid ..." >>"$LOGFILE" 
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  # kill $server_pid

  echo "Done: user and screen privacy expiration scripts: $(ts)" >> "$LOGFILE"
}

if [[ $OPTION == "real_run" ]]; then
  REAL_RUN=1
fi

maybe_activate_virtualenv

run_expiration_scripts
#test_setup
  





function archive_run_expiration_scripts_using_public_server {
  # Run user and screen privacy expiration scripts
  echo "run user and screen privacy expiration scripts: $(ts) ..." >> "$LOGFILE"

  echo "REAL RUN: $REAL_RUN"
  TEST_RUN_SETTINGS=" -admin_email_only -test_only " 
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"
  if [[ $REAL_RUN -eq 1 ]]; then
    TEST_RUN_SETTINGS=""
  fi
  echo "TEST_RUN_SETTINGS: $TEST_RUN_SETTINGS"
  
  echo "1.a Send SM user data privacy expiration notifications..." >>"$LOGFILE"
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type sm -days_to_expire $USER_DSL_DAYS_TO_EXPIRE -days_ahead_to_notify $USER_DSL_DAYS_AHEAD_TO_NOTIFY \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_log_filename ../logs/mail_user_agreement_notification.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.b Expire SM user agreements ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type sm -days_to_expire $USER_DSL_DAYS_TO_EXPIRE \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_log_filename ../logs/mail_user_agreement_expiration.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.c Send RNAi user data privacy expiration notifications..." >>"$LOGFILE"
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type rnai -days_to_expire $USER_DSL_DAYS_TO_EXPIRE -days_ahead_to_notify $USER_DSL_DAYS_AHEAD_TO_NOTIFY \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_log_filename ../logs/mail_user_agreement_notification.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "1.d Expire RNAi user agreements ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/user_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -ua_type rnai -days_to_expire $USER_DSL_DAYS_TO_EXPIRE \
  -email_message_directory db/static/user_agreement/ \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_log_filename ../logs/mail_user_agreement_expiration.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.a Adjust screen Data Privacy Expiration Dates ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -adjust_expiration_days_from_activity $SCREEN_DSL_EXPIRE_DAYS_FROM_ACTIVITY \
  -email_log_filename ../logs/mail_screen_dped_adjust.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.b Notify of screen privacy expirations ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -days_ahead_to_notify $SCREEN_DSL_DAYS_AHEAD_TO_NOTIFY \
  -email_log_filename ../logs/mail_screen_dped_notification.log \
  -v $TEST_RUN_SETTINGS -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.c Expire screen data sharing levels ..." >>"$LOGFILE" 2>&1
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -expire \
  -email_log_filename ../logs/mail_screen_dped_expiration.log \
  -v $TEST_RUN_SETTINGS -test_only -admin_email_only >>"$LOGFILE" 2>&1

  echo "2.d Notify admins for screen publications ..." >>"$LOGFILE" 2>&1
  echo "Using pass file: ${credential_file}" >>"$LOGFILE"
  PYTHONPATH=${BASEDIR} python \
  db/support/screen_privacy_expiration_emailer.py \
  -c ${credential_file} \
  -u ${SERVER_URL} \
  -contact_info "${CONTACT_PRINTED_INFO}" \
  -admin_from_email ${ADMIN_FROM_EMAIL} \
  -email_message_directory db/static/screen_privacy/ \
  -screen_type sm -notifyofpublications \
  -email_log_filename ../logs/mail_screen_notifyofpublications.log \
  -v $TEST_RUN_SETTINGS -test_only -admin_email_only >>"$LOGFILE" 2>&1

  ####
  echo "run_expirations finished, stop server ..." >>"$LOGFILE" 
  final_server_pid=$(ps aux |grep runserver| grep ${BOOTSTRAP_PORT} | awk '{print $2}')
  echo "kill $final_server_pid ..." >>"$LOGFILE" 
  kill $final_server_pid || error "kill server $final_server_pid failed with $?"
  # kill $server_pid


  echo "Done: user and screen privacy expiration scripts: $(ts)" >> "$LOGFILE"
}
