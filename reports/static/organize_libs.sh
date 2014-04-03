#!/bin/bash

if [[ $# -lt 1 ]]
then 
  echo "Usage: $0 [bower_dir] [dest_dir]"
  echo "bower_dir - directory containing /bower_components dir"
  echo "dest_dir - directory containing /js and /css dirs"
  exit $WRONG_ARGS
fi

# unfancy script that copies my js/css libs from the bower components directory 
# to the deploy location

BOWER_DIR=$1

# TODO: grab the minified versions
# TODO: concat, minify, testing...

PACKAGES=(  'requirejs/require.js'
            'jquery/jquery.js'
            'underscore/underscore.js' 
            'backbone-amd/backbone.js' 
            'backbone-pageable/lib/backbone-pageable.js'
            'backgrid/lib/backgrid.js'
            'backgrid-filter/backgrid-filter.js'
            'backgrid-paginator/backgrid-paginator.js'
            'backgrid-select-all/backgrid-select-all.js'
            'backbone.stickit/backbone.stickit.js'
            'Backbone.ModelBinder/Backbone.ModelBinder.js'
            'backbone-forms/distribution/backbone-forms.js'
            'layoutmanager/backbone.layoutmanager.js'
            'lunr.js/lunr.js'
            'requirejs-text/text.js' 
            'bootstrap/dist/js/bootstrap.js'         
          )
          
CSS=( 'Bootstrap-Admin-Theme/bootstrap/css/bootstrap.css'
      'Bootstrap-Admin-Theme/bootstrap/css/bootstrap-responsive.css'
      'backgrid/lib/backgrid.css'
      'backgrid-paginator/backgrid-paginator.css'
      'backgrid-filter/backgrid-filter.css'
    )

DEST_DIR=$2/js/libs

rm "${DEST_DIR}/*.js"

for package in ${PACKAGES[*]}; do
  command="cp $BOWER_DIR/bower_components/$package $DEST_DIR"
  echo "executing: $command"
  $($command)
done 

DEST_DIR=$2/css/external

rm "${DEST_DIR}/*.css"

for css in ${CSS[*]}; do
  command="cp $BOWER_DIR/bower_components/$css $DEST_DIR"
  echo "executing: $command"
  $($command)
done 