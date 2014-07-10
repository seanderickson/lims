define(function(require) {
  var models = require('test/models/models');
 
  describe('Models', function() {
    beforeEach(function(){
      console.log('running outer');
    });
    
    describe('Sample Model', function() {
      it('should default "urlRoot" property to "/api/samples"', function() {
        var sample = new models.Sample();
        sample.urlRoot.should.equal('/api/samples');
      });
    });
 
  });
 
});