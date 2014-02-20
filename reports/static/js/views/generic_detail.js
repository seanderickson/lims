define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'backbone_forms',
    'iccbl_backgrid',
    'text!templates/generic-detail.html'

], function( $, _, Backbone, stickit, backbone_forms, Iccbl, detailTemplate) {
	
	var DetailView = Backbone.Layout.extend({
    
    serialize: function() {
    
      var schema = this.model.resourceSchema;
      var keys = Iccbl.sortOnOrdinal(
          _.keys(this.model.attributes), schema.fields)
      var detailKeys = _(keys).filter(function(key){
        return _.has(schema.fields, key) && 
            _.has(schema.fields[key], 'visibility') && 
            _.contains(schema.fields[key]['visibility'], 'detail');
      });
      return {
        'fieldDefinitions': schema.fields,
        'title': Iccbl.getTitleFromTitleAttribute(this.model, schema),
        'keys': _.chain(detailKeys),
        'object': this.model.attributes
      };      
    },    
    
    afterRender: function() { 
      this.stickit(this.model);
    },
  
    template: _.template(detailTemplate)

	});

	return DetailView;
});