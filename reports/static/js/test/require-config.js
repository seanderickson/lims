require.config({
  baseUrl: '../',
  
  paths: {
//    'mocha': 'test/resources/mocha',
    'chai': 'test/resources/chai',
    'chai-jquery': 'test/resources/chai-jquery',
    'sinon': 'test/resources/sinon',
    'xhr': 'test/resources/sinon/util/fake_xml_http_request',
    'fakeServer': 'test/resources/sinon/util/fake_server',
    'event': 'test/resources/sinon/util/event',
    'tests': 'test/tests/'
  },
  shim: {
    'chai-jquery': ['jquery', 'chai'],
//     'QUnit': {
//         exports: 'QUnit',
//         init: function() {
//             QUnit.config.autoload = false;
//             QUnit.config.autostart = false;
//         }
//     },
     'sinon': {
       exports: 'sinon'
     },
     'xhr': {
       deps: ['sinon','event'],
       exports: 'xhr'
     },
     'fakeServer': {
       deps: ['sinon','xhr'],
       exports: 'fakeServer'
     }
  }
//  ,
//  urlArgs: 'bust=' + (new Date()).getTime()
});

// Chai
var should = chai.should();

//Protect from barfs
console = window.console || function() {};
 

