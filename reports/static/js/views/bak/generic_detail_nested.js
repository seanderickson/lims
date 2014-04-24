define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_modelbinder',
    'iccbl_backgrid',
    'models/app_state',
    'text!templates/generic-detail-modelbinder.html'

], function( $, _, Backbone, modelbinder, Iccbl, appModel, detailTemplate) {
	
  //FIXME: 20140402 - convert this to backbone-stickit:
  // prefer it's binding mechanism: xpath->to->model, 
  // rather than model->to->xpath like backbone_modelbinder
	var DetailView = Backbone.Layout.extend({
    
	  initialize: function(args) {
	    var self = this;
      var schema = this.model.resource.schema;
      
      var keys = this.detailKeys = schema.detailKeys(this.model);
      
      //var enhancedModel = this.enhancedModel = new Backbone.Model({});
      var descriptions = this.descriptions = new Backbone.Model.extend({});
      _.each(keys, function(key) {
      });
      
      var reagent = this.reagent = new Backbone.Model(this.model.get('reagent'));
      this.model.set('reagent', reagent );
	    
	  },

    afterRender : function() {
      var binder = new Backbone.ModelBinder();
      binder.bind(this.model, this.el);
      return this;
    },

    serialize: function() {
      var schema = this.model.resource.schema;
    
      return {
        'title': Iccbl.getTitleFromTitleAttribute(this.model, schema),
        'keys': _.chain(this.detailKeys),
      };      
    },    
    
    template: _.template(detailTemplate)

	});

	return DetailView;
});