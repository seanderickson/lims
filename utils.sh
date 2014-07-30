#!/usr/bin/env bash

function warn () {
  echo "$PROG: $@" >&2
}

function error () {
  warn "$@"
  exit 1
}

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
