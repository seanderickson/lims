define(['models/app_state',
        'text!test/models/ui_resources_fixture.json',
        'text!test/models/resource_from_server.js',
        'text!test/models/menu_fixture.json'
        ], function(appModel, ui_resources_raw, resource_raw, menu_raw) {
    
  var run = function() {
    module("Resource from server test", {
      setup: function () {
          console.log('setup the fake server: ' + resource_raw);
          this.server = sinon.fakeServer.create();
          this.server.respondWith(
              "GET", "/reports/api/v1/resource?limit=0", 
              [200, { "Content-Type": "application/json" }, '{ "objects": [' + resource_raw + "] }" ]);
      },
      teardown: function () {
          this.server.restore();
      }
    });  
    
    test( "a basic appModel test", function() {
      equal(appModel.get('api_root_url'), '/reports/api/v1', 
            "Expect value to be '/reports/api/v1'" );
    });
    
    test( "set/get resource", function() {
      appModel.set('ui_resources', JSON.parse(ui_resources_raw));
      appModel.set('menu', JSON.parse(menu_raw));
      
      var resources = appModel.get('ui_resources');
      ok(resources, 'resources loaded check');
      ok(resources['admin'], 'resources admin check');
    });
    
    test( "tests that server resource is merged to ui resource", function() {
      console.log('run async test');
      expect(2);
      appModel.set('ui_resources', JSON.parse(ui_resources_raw));
      appModel.set('menu', JSON.parse(menu_raw));
      stop();
      appModel.start(function(){
        console.log('finished appModel start');
        var resource = appModel.getResource('resource')

        ok(resource, 'resource loaded check');
        ok(!_.isUndefined(resource['schema']), 
            'resource has server (schema) information added in');
        
        start();
        console.log('finished async test');
      });
      this.server.respond();
      
    });
  };
  
  // TODO: should we just implicitly run this?
  return {run: run};
});