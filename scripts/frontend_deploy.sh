#!/usr/bin/env bash

##
# Rebuild the frontend (o2):
# - check
# Prerequisites
# - Server is configured: Python 2.7.x, virtualenv
# - Node, NPM
#
# NOTE: this script requires database credentials in .pgpass
##

if [[ $# -eq 1 ]]; then
  OPTION=$1
  if [[ "$OPTION" != "setup" ]]; then
    echo "Usage: $0 [setup] (optional, to run npm install)"
    exit 0
  fi
fi


REALPATH=${REALPATH:-"$(which realpath 2>/dev/null)"}
if [[ -z $REALPATH ]]; then
  PROG=$( basename $0 )
  error 'cannot find realpath'
fi

DEBUG=false

SCRIPTPATH="$($REALPATH $0)"
SCRIPTPATH="$(dirname $SCRIPTPATH)"
# strip the "scripts" directory to get the real BASEDIR
BASEDIR=${BASEDIR:-"$(dirname $SCRIPTPATH)"}
SUPPORTDIR=${SUPPORTDIR-"$(dirname $BASEDIR)"}
LOGFILE=${LOGFILE:-${BASEDIR:+"$BASEDIR/frontend_deploy.log"}}

if [[ ! $(command -v npm) ]]; then
  NODE_PATH=${NODE_PATH:-${SUPPORTDIR:+"$SUPPORTDIR/node/bin"}}
  if [[ -n "$NODE_PATH" ]]; then
    export PATH=${PATH}:$NODE_PATH
  fi
  if [[ ! $(command -v npm) ]]; then
    error "npm command not found"
  fi
fi

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

cd $BASEDIR

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

  maybe_activate_virtualenv
  
  $DJANGO_CMD collectstatic --noinput --clear \
      --ignore="*node_modules*" \
      --ignore="*demo_data*" \
      --ignore="*screen_privacy*" \
      --ignore="*user_agreement*" \
      --ignore="*test_data*" \
      --ignore="*.json" \
      --ignore="Gruntfile.js" \
      --ignore="*api_init*" \
      --ignore="js*" || error "collectstatic failed: $?"
  # FIXME: image and css files must be copied still
  # --ignore="css*"
  # --ignroe="images*"      
        
  if [ -e ../wsgi/app.wsgi ]; then
    touch ../wsgi/app.wsgi
  fi

  echo "frontend_deploy done: $(ts)" >> "$LOGFILE"

}

if [[ $OPTION == "setup" ]]; then
  echo "running frontend setup..."
  frontend_setup
fi
frontend_deploy
  

