#!/usr/bin/env bash

# Runs the UI and Server unit test suites
function error () {
  warn "$@"
  exit 1
}



# frontend tests

  cd reports/static 2>&1
  pwd
  npm --version 2>&1 || error "npm not found: $?"
  
  npm install  2>&1 || error "npm install failed: $?"
  npm install bower
  ./node_modules/.bin/bower cache clean
  ./node_modules/.bin/grunt bowercopy 2>&1 || error "grunt bowercopy failed: $?"
  
  ./node_modules/.bin/grunt test 2>&1 || error "grunt test failed: $?"
  
  cd ../..
  pwd
