define(function(require) {
  var models = require('test/models/models');
 
  describe('Reference Backbone model test', function() {
    beforeEach(function(){
      console.log('running outer');
    });
    
    it('should default "urlRoot" property to "/api/samples"', function() {
      var sample = new models.Sample();
      sample.urlRoot.should.equal('/api/samples');
    });
 
  });
 
});