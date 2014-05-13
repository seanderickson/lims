require.config({
  baseUrl: '../',
  
  paths: {
    'QUnit': 'test/resources/qunit',
    'sinon': 'test/resources/sinon',
    'xhr': 'test/resources/fake_xml_http_request',
    'fakeServer': 'test/resources/fake_server',
    'event': 'test/resources/event',
    'tests': 'test/tests/'
  },
  shim: {
     'QUnit': {
         exports: 'QUnit',
         init: function() {
             QUnit.config.autoload = false;
             QUnit.config.autostart = false;
         }
     },
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
});