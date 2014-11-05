/*global module*/
module.exports = function (grunt) {
  'use strict';
  
  var gruntConfig = {};
  
  // specify a local PhantomJS server, see
  // https://github.com/gruntjs/grunt-contrib-qunit & 
  // https://github.com/gruntjs/grunt-contrib-connect
  gruntConfig.connect = {
    server: {
      options: {
        port: 8008,
        base: '.',
        debug: true
      }
    }
  };

  gruntConfig.mocha = {
    test: {
      options: {
        urls: ['http://localhost:8008/js/test/index_tests_mocha.html'],
        //        reporter: 'Spec',
        //        reporter: 'Min',
        reporter: 'JSON' ,
        //        reporter: 'List',
        run: false,
        timeout: 25000,
        log: true
      }
    }
    
  };
  
  gruntConfig.bowercopy = {
    options: {
        srcPrefix: 'bower_components',
        clean: true
    },
    scripts: {
        options: {
            destPrefix: 'js/libs'
        },
        files: {
          'require.js' : 'requirejs/require.js',
          'jquery.js' : 'jquery/dist/jquery.js',
          'underscore.js' : 'underscore/underscore.js', 
          'backbone.js' : 'backbone-amd/backbone.js', 
//          'backbone-pageable.js' : 'backbone-pageable/lib/backbone-pageable.js',
          'backbone.paginator': 'backbone.paginator/lib/backbone.paginator.js',
          'backgrid.js' : 'backgrid/lib/backgrid.js',
          'backgrid-filter.js' : 'backgrid-filter/backgrid-filter.js',
          'backgrid-paginator.js' : 'backgrid-paginator/backgrid-paginator.js',
          'backgrid-select-all.js' : 'backgrid-select-all/backgrid-select-all.js',
          'backbone.stickit.js' : 'backbone.stickit/backbone.stickit.js',
          'Backbone.ModelBinder.js' : 'Backbone.ModelBinder/Backbone.ModelBinder.js',
          'backbone-forms.js' : 'backbone-forms/distribution/backbone-forms.js',
          'backbone.layoutmanager.js' : 'layoutmanager/backbone.layoutmanager.js',
          'lunr.js' : 'lunr.js/lunr.js',
          'text.js' : 'requirejs-text/text.js', 
          'bootstrap.js' : 'bootstrap/dist/js/bootstrap.js'
        }
    }, 
    css: {
      options: {
        destPrefix: 'css/external'
      },
      files: {
        'bootstrap.css': 'bootstrap/dist/css/bootstrap.css',
        'backgrid.css': 'backgrid/lib/backgrid.css',
        'backgrid-paginator.css': 'backgrid-paginator/backgrid-paginator.css',
        'backgrid-filter.css': 'backgrid-filter/backgrid-filter.css'
      }
    },
    fonts: {
      files: {
        'css/fonts': 'bootstrap/dist/fonts/*'
      }
    }, 
    
    devfiles: {
      options: {
        destPrefix: 'js/test/resources'
      },
      files: {
        'mocha.js': 'mocha/mocha.js',
        'mocha.css': 'mocha/mocha.css',
        'chai.js': 'chai/chai.js',
        'chai-jquery.js': 'chai-jquery/chai-jquery.js',
        'sinon.js': 'sinon/lib/sinon.js',
        'sinon': 'sinon/lib/sinon/*'
      }
    }
  };

  grunt.initConfig(gruntConfig);

  grunt.loadNpmTasks('grunt-contrib-jshint');
  grunt.loadNpmTasks('grunt-bowercopy');
  grunt.loadNpmTasks('grunt-contrib-connect');
  grunt.loadNpmTasks('grunt-mocha');
  grunt.registerTask('test', ['connect', 'mocha']);
  
};
