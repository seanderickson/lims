// provisional...
define([
    'jquery',
    'underscore',
    'backbone',
    'views/detail_stickit',
    'text!templates/generic-collection-columns.html',
], function($, _, Backbone, DetailView, collectionColumnsTemplate) {
    var CollectionView = Backbone.View.extend({

        initialize : function(attributes, options) {
            console.log('initialize CollectionColumns view');
            var self = this;
            this.options = options;
            Iccbl.assert( !_.isUndefined(options.schemaResult), 'collection view requires a schemaResult struct');
            Iccbl.assert( !_.isUndefined(options.schemaResult['resource_definition']), 'collection view schemaResult requires a resource_definition');
            Iccbl.assert( !_.isUndefined(options.router), 'collection view requires a router');
            Iccbl.assert( !_.isUndefined(options.collection), 'collection view requires a collection');

            // Headers
            headers = [];
            rows = [];
            this.attributeKeys = Iccbl.sortOnOrdinal(_.keys(options.schemaResult.fields), options.schemaResult.fields);
            console.log('attribute keys: ' + JSON.stringify(this.attributeKeys));
            var resource_definition = self.options.schemaResult['resource_definition'];
            options.collection.each(function(model){
                var header = Iccbl.getTitleFromTitleAttribute(model, self.options.schemaResult);
                 // _.reduce(resource_definition['title_attribute'],
                    // function(memo, item){
                        // if( model.has(item) ) memo += model.get(item)
                        // else memo += item
                        // return memo ;
                    // }, '');

                headers.push(header);
            });
            headers.unshift(self.options.schemaResult.fields[resource_definition['title_attribute']]['title']);
            console.log('attributeKeys: ' + JSON.stringify(this.attributeKeys));
            _.each(this.attributeKeys, function(key){
                var attributeTitle = self.options.schemaResult.fields[key]['title'];
                if(_.isUndefined(attributeTitle)) attributeTitle = key;
                row = [{ value: attributeTitle} ]
                self.options.collection.each(function(model){
                    row.push( { value: model.get(key)} );
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