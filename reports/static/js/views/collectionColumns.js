
define([
    'jquery',
    'underscore',
    'backbone',
    'text!templates/generic-collection-columns.html',
], function($, _, Backbone, collectionColumnsTemplate) {

  /**
   * Show an in memory collection as a series of columns rather than rows.
   */
  var CollectionView = Backbone.View.extend({

    initialize : function(attributes, options) {
      console.log('initialize CollectionColumns view');
      var self = this;
      this.options = options;
      Iccbl.assert(
          !_.isUndefined(options.schemaResult), 
          'collection view requires a schemaResult struct');
      Iccbl.assert(
          !_.isUndefined(options.router), 'collection view requires a router');
      Iccbl.assert(
          !_.isUndefined(options.collection), 'collection view requires a collection');

      headers = [];
      rows = [];
      this.attributeKeys = Iccbl.sortOnOrdinal(
          _.keys(options.schemaResult.fields), options.schemaResult.fields);
      var resource_definition = self.options.schemaResult;
      options.collection.each(function(model){
          headers.push(
              Iccbl.getTitleFromTitleAttribute(model, self.options.schemaResult));
      });
      headers.unshift(
          self.options.schemaResult['title']);
      _.each(
          this.attributeKeys, 
            function(key){
              var attributeTitle = self.options.schemaResult.fields[key]['title'];
              if(_.isUndefined(attributeTitle)) attributeTitle = key;
              row = [{ value: attributeTitle} ]
              self.options.collection.each(function(model){
                  row.push( {value: model.get(key)} );
              });
              rows.push(row);
            });

      var data = {
          headers: headers, rows: rows,
          caption: resource_definition['description']
          };

      var compiledTemplate = this.compiledTemplate = _.template( collectionColumnsTemplate, data );
    },

    render : function() {
      var self = this;
      this.$el.html(this.compiledTemplate);
      return this;
    },
  });

  return CollectionView;
});