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
            this.options = options;
            Iccbl.assert( !_.isUndefined(options.columns), 'simple list view options.columns is required');
            Iccbl.assert( !_.isUndefined(options.collection), 'simple list view requires a collection');

            this.objects_to_destroy = _([]);

            var data = { message: '' };
            var compiledTemplate = this.compiledTemplate = _.template( listTemplate, data );

            var grid = this.grid = new Backgrid.Grid({
              columns: this.options.columns,
              collection: this.options.collection,
            });

            this.objects_to_destroy.push(grid);

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
            this.$el.html(this.compiledTemplate);
            this.$("#example-table").append(this.grid.render().$el);
            this.delegateEvents();
            return this;
        }

    });

  return ListView;
});