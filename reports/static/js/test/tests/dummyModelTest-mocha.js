define([
        'sinon', 'fakeServer', 'chai',
        'underscore',
        'backbone',
        'models/app_state',
        'text!test/models/ui_resources_fixture.json',
        'text!test/models/resource_from_server.js'
        ], 
      function(sinon, fakeServer, chai, _, Backbone, 
          appModel, ui_resources_raw, resource_raw) {
  var testData = {           
      object1: {
          property1: 'val1',
          property2: 'val2'
      }
  };
  var expect = chai.expect;
  
  describe('Reference ajax/backbone/fakeserver test: ', function() {
    before(function(){
      this.server = sinon.fakeServer.create();
      this.server.autoRespond = true;
      this.server.respondWith(
          "GET", "/api/testmodel/1", 
          [200, { "Content-Type": "application/json" }, 
          JSON.stringify(testData)]);
      
      
    });
    after(function(){
      this.server.restore();
    });    
    
    it('should retrieve testdata on fetch"', function(done) {
      var TestModel = Backbone.Model.extend({
        url: function() { return '/api/testmodel' + (this.id ? '/'+this.id : ''); }
      });
      var model = new TestModel({ id: 1 });
      model.fetch({
        success: function (m, resp) {
          console.log('success: %o, %o', m, resp);
          expect(model.get('object1')).eql(testData['object1']); 
          model.get('object1').should.eql(testData['object1']); 
          if(done){
            console.log('done called');
            done();
          }
        }, 
        error: function (m, resp) {
          console.log('error: %o, %o', m, resp);
        }
      });
      //      console.log('waiting...');
      //      this.server.respond();
      //      console.log('done waiting');
      
    });
 
  });
  
  
});