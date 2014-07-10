define([
        'sinon', 'fakeServer', 
        'underscore',
        'backbone',
        'models/app_state'], 
      function(sinon, fakeServer, _, Backbone, appModel) {
  var testData = {           
      object1: {
          property1: 'val1',
          property2: 'val2'
      }
  };
  
  describe('Reference ajax/backbone/fakeserver test', function() {
    before(function(){
      this.server = sinon.fakeServer.create();
      this.server.respondWith(
          "GET", "/api/testmodel/1", 
          [200, { "Content-Type": "application/json" }, 
          JSON.stringify(testData)]);
      
      
    });
    after(function(){
      this.server.restore();
    });    
    
    describe('simple get model', function(done) {
      it('should retrieve testdata on fetch"', function() {
        var TestModel = Backbone.Model.extend({
          url: function() { return '/api/testmodel' + (this.id ? '/'+this.id : ''); }
        });
        var model = new TestModel({ id: 1 });
        model.fetch({
          success: function (m, resp) {
            console.log('success: %o, %o', m, resp);
            if(done) done();
          }, 
          error: function (m, resp) {
            console.log('error: %o, %o', m, resp);
          }
        });
        this.server.respond();
        model.get('object1').should.eql(testData['object1']); 
      });
    });
 
  });
  
  
});