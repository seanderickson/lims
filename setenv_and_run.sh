#!/usr/bin/env bash

##
# HMS - ICCBL specific:
# 
# This script will grab the postgres connection variables from the appropriate apache conf file
# and set them to the environment.  This will make them available to settings.py
# when run (via manage.py).
# - This script makes it possible to run manage.py easily from the commmand line.
#
# NOTES:
# 
# - we are using the lims/settings-server-commandline.py settings file to use ../logs/screensaver2-cli.log
# the advantage of this is that we won't be overwriting or changing perms of the server logs
# (** this will cause a server fault)
#
# Orchestra apache auth conf files are located at:
# /opt/apache/conf/auth/[server_address]
#
# Examples:
#
# ./setenv_and_run.sh /opt/apache/conf/auth/dev.screensaver2.med.harvard.edu ./manage.py shell
# 
# System bootstrapping:
# ./setenv_and_run.sh /opt/apache/conf/auth/dev.screensaver2.med.harvard.edu ./manage.py db_init  --inputDirectory=./reports/static/api_init/ -f ./reports/static/api_init/api_init_actions.csv -u http://localhost:55001/reports/api/v1 -a 
#
##

function warn () {
  echo "$PROG: $@" >&2
}

function error () {
  warn "$@"
  exit 1
}

REALPATH=${REALPATH:-"$(which realpath 2>/dev/null)"}
if [[ -z $REALPATH ]]; then
  PROG=$( basename $0 )
  error 'cannot find realpath'
fi
SCRIPTPATH="$($REALPATH $0)"
PROG=$(basename $SCRIPTPATH)
BASEDIR=${BASEDIR:-"$(dirname $SCRIPTPATH)"}
SUPPORTDIR=${SUPPORTDIR-"$(dirname $BASEDIR)"}

## 
# Activate the virtualenv
##

VENV=${VENV:-${SUPPORTDIR:+"$SUPPORTDIR/virtualenv"}}
if [[ -z $VENV ]]; then
  error 'no virtualenv available'
fi

function invenv {
  local pyex="$VENV/bin/python"
  if [[ -z $VENV || ! -x $pyex ]]; then
    error "VENV ($VENV) does not point to a valid virtualenv directory"
  fi

  [[ $VIRTUAL_ENV == $VENV && $(which python) == $pyex ]] && \
    python -c 'import sys; sys.real_prefix' 2>/dev/null
}

function maybe_activate_virtualenv {
  if [[ -z $VENV ]]; then
    error 'assertion error'
  fi

  if [[ ! -e $VENV ]]; then
    virtualenv --no-site-packages --distribute $VENV
  fi
  VENV="$($REALPATH $VENV)"

  if invenv; then return; fi

  source "$VENV/bin/activate" || error "failed to activate virtualenv: $?"

  if invenv; then return; fi

  error 'failed to activate virtualenv (reason unknown)'
}
#######

if [[ $# -lt 1 ]]
then
  error "Usage: $0 [path_to_apache_conf_auth_file]"
fi

AUTH_FILE=$1
shift

if [[ ! -e $AUTH_FILE ]]; then
  error "apache auth file not found at: \"$AUTH_FILE\""
fi 
 
#set -x

for x in `sed  's/SetEnv\s\+\(\S\+\)\s\+\(.*\)$/\1=\2/;tx;d;:x' $AUTH_FILE `; do 
  export $x  
done

maybe_activate_virtualenv


# set an overlay settings.py (overrides logging settings, so we don't 
# overwrite/or change the perms for the server logs)
export DJANGO_SETTINGS_MODULE=lims.settings-server-commandline

if [[ ! -z $@ ]]; then
  exec $@
fi
