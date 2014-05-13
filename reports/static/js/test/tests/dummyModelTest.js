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
  var run = function() {
    
    module("Reference ajax/backbone test", {
      setup: function () {
          this.server = sinon.fakeServer.create();
          this.server.respondWith(
              "GET", "/api/testmodel/1", 
              [200, { "Content-Type": "application/json" }, 
              JSON.stringify(testData)]);
      },
      teardown: function () {
          this.server.restore();
      }
    });
   
    test( "dummy Backbone model test", function() {
      var TestModel = Backbone.Model.extend({
        url: function() { return '/api/testmodel' + (this.id ? '/'+this.id : ''); }
      });
      var model = new TestModel({ id: 1 });
      model.fetch({
        success: function (m, resp) {
          console.log('success: %o, %o', m, resp);
        }, 
        error: function (m, resp) {
          console.log('error: %o, %o', m, resp);
        }
      });
      this.server.respond();
   
      // TODO: equal should be working here, not sure why not
      propEqual(model.get('object1'),testData['object1'] , 
            'fetched model has expected attributes: ' + JSON.stringify(model.attributes));
    });
  };

  return {run: run}
});