#!/usr/bin/env bash


function warn () {
  echo "$PROG: $@" >&2
}

function error () {
  warn "$@"
  exit 1
}

function read_auth_file() {
  auth_file=$1
  if [[ ! -e $auth_file ]]; then
     error "apache auth file not found at: \"$auth_file\""
  else 

#set -x

  for x in `sed  's/SetEnv\s\+\(\S\+\)\s\+\(.*\)$/\1=\2/;tx;d;:x' $auth_file `; do
    export $x  
  done
  fi
}

function invenv {
  if [[ -z "${VIRTUAL_ENV}" ]] 
  then
    local pyex="$VENV/bin/python"
    if [[ -z $VENV || ! -x $pyex ]]; then
      error "VENV ($VENV) does not point to a valid virtualenv directory"
    fi
  
    [[ $VIRTUAL_ENV == $VENV && $(which python) == $pyex ]] && \
      python -c 'import sys; sys.real_prefix' 2>/dev/null
  else
    echo "Using already active virtual environment: \"${VIRTUAL_ENV}\" "
  fi
      
}

function maybe_activate_virtualenv {

  if invenv; then return; fi

  if [[ -z $VENV ]]; then
    error "no VENV environment variable set"
  fi

  if [[ ! -e $VENV ]]; then
    error "VENV: \"$VENV\" does not exist"
  fi
  VENV="$($REALPATH $VENV)"

  if type "module" > /dev/null; then
    module load gcc/6.2.0
    module load python/2.7.12
    module load perl
  fi

  source "$VENV/bin/activate" || error "failed to activate virtualenv: $?"

  if invenv; then return; fi

  error 'failed to activate virtualenv (reason unknown)'
}

function maybe_activate_virtualenv_webconf {
  # Activate the virtualenv in the webconf environment:
  # modules do not need to be loaded

  if [[ -z $VENV ]]; then
    error 'assertion error'
  fi

  VENV="$($REALPATH $VENV)"

  if invenv; then return; fi

#  module load gcc/6.2.0
#  module load python/2.7.12
#  module load perl

  source "$VENV/bin/activate" || error "failed to activate virtualenv: $?"

  if invenv; then return; fi

  error 'failed to activate virtualenv (reason unknown)'
}

function ts {
  date +%Y%m%dT%H%M%S%z
}

