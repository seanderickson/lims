define([
  'jquery',
  'underscore',
  'backbone',
  'backbone_pageable',
  'backgrid',
  'iccbl_backgrid',
  'text!templates/simple-list.html',
  'text!templates/modal_ok_cancel.html',
  'text!templates/generic-selector.html',
], function($, _, Backbone, BackbonePageableCollection, Backgrid, Iccbl, listTemplate, modalTemplate, genericSelectorTemplate ){

    // for compatibility with require.js, attach PageableCollection in the right place on the Backbone object
    // see https://github.com/wyuenho/backbone-pageable/issues/62
    Backbone.PageableCollection = BackbonePageableCollection;

    var ajaxStart = function(){
        $('#loading').fadeIn({duration:100});
    };
    var ajaxComplete = function(){
        $('#loading').fadeOut({duration:100});
    };
    $(document).bind("ajaxComplete", function(){
        ajaxComplete(); // TODO: bind this closer to the collection
    });

    var ListView = Backbone.View.extend({

        initialize : function(attributes, options) {
            var self = this;
            Backgrid.requireOptions(options, ["columns", "collection", "schemaResult"]);

            this.options = options;
            var schemaResult = this.schemaResult = options.schemaResult;

            this.objects_to_destroy = _([]);

            var data = { message: '' };
            var compiledTemplate = this.compiledTemplate = _.template( listTemplate, data );

            var collection = this.collection = options.collection;

            var grid = this.grid = new Backgrid.Grid({
              columns: this.options.columns,
              collection: this.options.collection,
            });

            this.objects_to_destroy.push(grid);

            // Paginator
            var paginator = self.paginator = new Backgrid.Extension.Paginator({
              collection: self.collection
            });
            this.objects_to_destroy.push(paginator);

            // Extraselector
            if( _.has(schemaResult, 'extraSelectorOptions')){
                var searchHash = self.model.get('search');
                console.log('extraselector init: searchTerms: ' + JSON.stringify(searchHash));

                var extraSelectorModel = new Backbone.Model({ selection: '' });
                var extraSelectorKey = schemaResult.extraSelectorOptions.searchColumn;
                _.each(_.keys(searchHash), function(key){
                    console.log('key: ' + key + ', extrSelectorKey: ' + extraSelectorKey);
                    if( key == extraSelectorKey){
                        extraSelectorModel.set({ selection: searchHash[key] });
                    }
                });
                var extraSelectorInstance = self.extraSelectorInstance =
                    new Iccbl.GenericSelector({ model: extraSelectorModel } , schemaResult.extraSelectorOptions );
                this.objects_to_destroy.push(extraSelectorInstance);

                this.listenTo(this.model, 'change: search', function(){
                    var searchHash = self.model.get('search');
                    console.log('extraselector, search changed: ' + JSON.stringify(searchHash));
                    _.each(_.keys(searchHash), function(key){
                        console.log('key: ' + key + ', extrSelectorKey: ' + extraSelectorKey);
                        if( key == extraSelectorKey){
                            extraSelectorModel.set({ selection: searchHash[key] });
                        }
                    });
                });
                this.listenTo(extraSelectorModel, 'change', function() {
                    console.log('===--- extraSelectorModel change');
                    var searchHash = _.clone(self.model.get('search'));
                    var value = extraSelectorModel.get('selection');
                    searchHash[extraSelectorKey] = value;
                    self.model.set('search', searchHash);
                });
            }

            // TODO: work out the specifics of communication complete event
            this.listenTo(self.collection, 'request', ajaxStart);
            this.listenTo(self.collection, 'error', ajaxComplete);
            this.listenTo(self.collection, 'complete', ajaxComplete);

            console.log('list view initialized');
        },


        remove: function(){
            console.log('ListView remove called');
            Backbone.View.prototype.remove.apply(this);
        },

        onClose: function(){
            console.log('Extra onclose method called');
            if(_.isObject(this.objects_to_destroy)){
                this.objects_to_destroy.each(function(view_obj){
                    view_obj.remove();
                    view_obj.unbind();
                    view_obj.stopListening();
                });
            }
        },

        render: function(){
            var self = this;
            this.$el.html(this.compiledTemplate);
            this.$("#example-table").append(this.grid.render().$el);
            self.$("#paginator-div").append(self.paginator.render().$el);
            if(!_.isUndefined(self.extraSelectorInstance)){
                self.$("#extra-selector-div").append(self.extraSelectorInstance.render().$el);
            }

            var page = this.model.get('page');
            if(!_.isUndefined(page)){
                this.grid.collection.setPage(page);
            }
            this.delegateEvents();
            return this;
        }

    });

  return ListView;
});