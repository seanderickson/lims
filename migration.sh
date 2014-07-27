#!/usr/bin/env bash

##
# HMS - ICCBL specific:
# 
# Migrate the Screensaver 1 database to the Screensaver 2 database
#
##

source ./setenv_and_run.sh

# Prerequisites
# - Python 2.7.x, pip, virtualenv
# - Node, NPM
# - git

SCRIPTPATH="$($REALPATH $0)"


if [[ $# -lt 2 ]]
then
  echo "Usage: $0 [branch] "
  exit $WRONG_ARGS
fi

BRANCH=$1
REMOTE=$2

LOGFILE=./migration.log

function _debug {
  if $DEBUG; then echo "$@"; fi
}

function gitpull {
  _debug -n "pulling branch $BRANCH from $REMOTE... "
  cp -a $SCRIPTPATH $SAVEPATH
  git fetch $REMOTE >>"$LOGFILE" 2>&1 || error "git-fetch failed: $?"
  git checkout $BRANCH >>"$LOGFILE" 2>&1 || error "git-checkout failed: $?"
  git checkout -- $SCRIPTPATH >>"$LOGFILE" 2>&1 || error "git-checkout $SCRIPTPATH failed: $?"
  git pull --ff-only $REMOTE $BRANCH >>"$LOGFILE" 2>&1 || error "git-pull failed: $?"
  mv -f $SAVEPATH $SCRIPTPATH
#  update_deploy_info "$STATUS"
  _debug 'done'
  return 0
}

# Steps:

# Git pull
# git fetch --all
# git merge seanderickson/master

gitpull

# start virtualenv

maybe_activate_virtualenv

pip install -r requirements.txt

# django test
./manage.py test


# cd reports/static
# npm install

# grunt bowercopy
# grunt test

# drop SS2 database
# import Screensaver database
# pg_restore -Fc --no-owner -h localhost -d screensaver_test -Uscreensaver_test screensaver.2014-01-10.excl_big_data.pg_dump 
# bsub -q short -W 4:0 ../../scripts/pg_restore.sh 2014-01-10 devscreensaver2 devscreensaver2web dev.pgsql.orchestra

# Migrate
./manage.py migrate reports

./manage.py migrate tastypie

./manage.py migrate db 0001 --fake  

./manage.py migrate db 0002

psql -Uuser database -f ./migrations/manual/0003_screen_status.sql

# run the rest of the migrations
./manage.py migrate db

# Bootstrapping the web server

# run a local server
./setenv_and_run.sh/opt/apache/conf/auth/dev.screensaver2.med.harvard.edu ./manage.py runserver 55001

# bootstrap the metahash 
PYTHONPATH=. python reports/utils/db_init.py  \
  --input_dir=./reports/static/api_init/ \
  -f ./reports/static/api_init/api_init_actions.csv \
  -u http://localhost:8000/reports/api/v1 -U sde
  
PYTHONPATH=. python reports/utils/db_init.py  \
  --input_dir=./db/static/api_init/ \
  -f ./db/static/api_init/api_init_actions.csv \
  -u http://localhost:8000/reports/api/v1 -U sde

# Setup the server:

PYTHONPATH=. python reports/utils/db_init.py  \
  --input_dir=./lims/static/production_data/ \
  -f ./lims/static/production_data/api_init_actions.csv \
  -u http://localhost:8000/reports/api/v1 -U sde



# run some grunt-mocha chai-jquery to test the UI?  