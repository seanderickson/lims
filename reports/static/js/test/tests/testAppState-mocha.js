define([
        'sinon', 'fakeServer', 'chai',
        'models/app_state',
        'test/models/ui_resources_fixture.json',
        'test/models/resource_from_server.json',
        'test/models/vocabulary_from_server.json',
        'test/models/user_resource_from_server.json',
        'test/models/test_user.json',
        'test/models/menu_fixture.json'
        ], function(sinon, fakeServer, chai, 
            appModel, ui_resources_raw, resource_raw, vocabulary_raw, 
            user_resource_raw, test_user_raw, menu_raw) {
  var expect = chai.expect;
  console.log('start testAppState');
  
  window.user = 'testuser';  
  
  describe("Resource from server test:", function(){
    before(function(){
        console.log('-> Resource from server test');
        this.server = sinon.fakeServer.create();
        this.server.autoRespond = true;

        // Note: load the server with two responses:
        // First for the getResources call
        this.server.respondWith(
            "GET", "/reports/api/v1/resource?limit=0", 
            [200, { "Content-Type": "application/json" }, 
                  '{ "objects": [' + resource_raw + "] }" ]);
        this.server.respondWith(
            "GET", "/reports/api/v1/vocabulary?limit=0", 
            [200, { "Content-Type": "application/json" }, 
                  '{ "objects": [' + vocabulary_raw + "] }" ]);
        // Second, for the getUser call
        this.server.respondWith(
            "GET", "reports/api/v1/user/testuser?includes=*", 
            [200, { "Content-Type": "application/json" }, 
                  '{ "objects": [' + test_user_raw + "] }" ]);
    });
    after(function(){
        this.server.restore();
    });
  
    it( "should report the correct api_root_url", function() {
      expect(appModel.get('api_root_url')).to.equal('/reports/api/v1');
    });
  
    it( "should set/get resource", function() {
      appModel.set('ui_resources', JSON.parse(ui_resources_raw));
      appModel.set('menu', JSON.parse(menu_raw));
      
      var resources = appModel.get('ui_resources');
      expect(resources);
      expect(resources['admin'], 'resources admin check');
    });
  
    it("server resource should be merged to ui resource", function(done) {
      console.log('run async test');
      appModel.set('ui_resources', JSON.parse(ui_resources_raw));
      appModel.set('menu', JSON.parse(menu_raw));
      appModel.start(function(){
        console.log('finished appModel start');
        var resource = appModel.getResource('resource')
        expect(resource, 'resource loaded check');
        expect(resource['schema'], 
            'resource has server (schema) information added in');
        console.log('finished async test');
        if(done){
          done();
        }
      });

      // Note: no need to manually respond as the autoRespond flag = true
      // respond once for the getResources call in appModel
      //      this.server.respond();
      //      // again for the user request that is in the appModel.start
      //      this.server.respond();
    });
  });  
});