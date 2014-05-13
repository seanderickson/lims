/**
 * Application tests loading script for the Iccbl-lims app.
 */

requirejs(['../require-config'], function() { // First, load application require-config
  //second, load the testing require-config, which sets the baseUrl to ".." :
  requirejs(['require-config'], function() {
    // load test code
    requirejs(['QUnit'], function(QUnit){
      requirejs(['tests/dummyModelTest'], function(dummyModelTest) {
        dummyModelTest.run();
      });
      requirejs(['tests/testAppState'], function(testAppState) {
        testAppState.run();
      });
        
      // start QUnit.
      QUnit.load();
      QUnit.start();  
    });
  });
});
