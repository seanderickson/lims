/*global module*/
module.exports = function (grunt) {
  'use strict';
 
//  grunt.event.on('qunit.spawn', function (url) {
//    grunt.log.ok("qunit.spawn - Running test: " + url);
//  });
  grunt.event.on('fail.timeout', function (url) {
    grunt.log.ok("qunit.fail.timeout: " + url);
  });
//grunt.event.on('qunit.log', function (result, actual, expected, message, source) {
//  grunt.log.ok("qunit.fail.timeout: " + result +',' 
//      + actual +',' + expected +',' +  message + ',' + source);
//});
  
  var gruntConfig = {};

//  gruntConfig.qunit = {
//    all: {
//      options: {
//        urls: [
//          'http://localhost:8008/js/test/index_tests.html'
//        ],
//        run: true,
//        ignoreLeaks: false
//      }
//    }
//  };
  
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

  gruntConfig.jshint = {
      options: { bitwise: true, camelcase: true, curly: true, eqeqeq: true, forin: true, immed: true,
          indent: 2, latedef: true, newcap: true, noarg: true, noempty: true, nonew: true, plusplus: true,
          quotmark: true, regexp: true, undef: true, unused: true, strict: true, trailing: true,
          maxparams: 3, maxdepth: 2, maxstatements: 500},
      all: [
          'Gruntfile.js',
          'js/models/**/*.js',
          'js/views/**/*.js',
          'js/*.js'
      ]
  };
  
  
  gruntConfig.mocha = {
    test: {
      options: {
        urls: ['http://localhost:8008/js/test/index_tests_mocha.html'],
        run: false,
      }
    }
    
  };
  
  gruntConfig.bowercopy = {
    options: {
        srcPrefix: 'bower_components'
    },
    scripts: {
        options: {
            destPrefix: 'js/libs'
        },
        files: {
          'require.js' : 'requirejs/require.js',
          'jquery.js' : 'jquery/jquery.js',
          'underscore.js' : 'underscore/underscore.js', 
          'backbone.js' : 'backbone-amd/backbone.js', 
          'backbone-pageable.js' : 'backbone-pageable/lib/backbone-pageable.js',
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
//          'bootstrap.css.map': 'bootstrap/dist/css/bootstrap.css.map',
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
//  grunt.loadNpmTasks('grunt-contrib-qunit');
  grunt.loadNpmTasks('grunt-contrib-connect');
  grunt.loadNpmTasks('grunt-mocha');
//grunt.registerTask('test', ['connect', 'qunit']);
//  grunt.registerTask('test', ['connect:server', 'mocha:test']);
  grunt.registerTask('test', ['connect', 'mocha']);
  

  
};