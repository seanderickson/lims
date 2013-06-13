define([
  'jquery',
  'underscore',
  'backbone',
  'backbone_pageable',
  'backgrid',
  'iccbl_backgrid',
  'text!templates/rows-per-page.html',
  'text!templates/list.html'
], function($, _, Backbone, BackbonePageableCollection, Backgrid, App, rowsPerPageTemplate, listTemplate ){

    // for compatibility with require.js, attach PageableCollection in the right place on the Backbone object
    // see https://github.com/wyuenho/backbone-pageable/issues/62
    Backbone.PageableCollection = BackbonePageableCollection;

    var ajaxStart = function(){
        console.log("start spinner...");
        $('#loading').fadeIn({duration:100});
    };
    var ajaxComplete = function(){
        console.log("stop spinner...");
        $('#loading').fadeOut({duration:100});
    };

  var ListView = Backbone.View.extend({

    options: { // TODO: needs to be set from menu changes/app model
        url_schema : '/reports/api/v1/fieldinformation/schema/',
        url : '/reports/api/v1/fieldinformation/?format=json',
        router: {
            navigate: function(route){
                console.log('todo: implement the router event notifications: ' + route);
            }
        }
    },

    initialize : function() {
        console.log('initialize ListView');
    },

    reset_grid : function(options){
        var self = this;
        // Call the schema from the server to get the field definitions
        $.ajax({
            type: "GET",
            url: options.url_schema,
            data: "",
            dataType: "json",
            success: function(result) {
                self.buildGrid(result.fields, options.url, options);

            }, // end success outer ajax call
            error: function(x, e) {
                alert(x.readyState + " "+ x.status +" "+ e.msg);
            }
        });
    },

    buildGrid : function(fieldDefinitions, _url, options) {
        console.log('buildGrid...');
        // TODO: make sure req'd options are present (i.e. router, type)
        var collection = new App.MyCollection({ 'url': _url, router: options.router, type: options.type });

        var columns = App.createBackgridColModel(fieldDefinitions, App.MyHeaderCell, collection);

        var grid = new Backgrid.Grid({
          columns: columns,
          collection: collection
        });

        var selector = new App.ItemsPerPageSelector(
            { 'selections': ['25','50','200','1000'], 'template': rowsPerPageTemplate },
            collection);

        var paginator = new Backgrid.Extension.Paginator({
          collection: collection
        });

        var data = { message: 'title, and other information' }; //'hello world!' };
        var compiledTemplate = _.template( listTemplate, data );
        // Append our compiled template to this Views "el"
        //$("#container").append( compiledTemplate );
//
        // $("#container").append(paginator.render().$el);
        // $("#container").append(selector.render().$el);
        // $("#container").append(grid.render().$el);

        this.el.innerHTML = compiledTemplate;
        $("#paginator-div").append(paginator.render().$el);
        $("#rows-selector-div").append(selector.render().$el);
        $("#example-table").append(grid.render().$el);

        var allEvent = function(event){
            console.log("event fired: " + event);
        };
        collection.on({
          "all": allEvent,  // TODO: this is debug code
        });

        // TODO: tie these events into the collection override?
        //        collection.bind('request', window.App.ajaxStart, this);
        collection.on('request', ajaxStart); // NOTE: can use bind or on
        collection.bind('sync', ajaxComplete, this);
        collection.bind('sync', selector.updateSelection, selector);  // TODO: selector.listenTo(collection, 'sync'...

        if (options.page){
            collection.state.currentPage = options.page;
        }
        if (options.pageSize){
            collection.state.pageSize = options.pageSize;
        }
        if (options.sortKey){
            collection.state.sortKey = options.sortKey;
            collection.state.order = options.order;
            // Notify header cells
            collection.trigger("backgrid:sort", options.sortKey, options.direction, null, collection);
        }

        if (options.searchBy){
            collection.searchBy = options.searchBy;
            console.log('trigger search');
            collection.trigger("MyServerSideFilter:search", options.searchColumn, options.searchTerm, collection);
            console.log('done: trigger search');
        }
        collection.fetch({ reset: true } );

        console.log('list view initialized');
    },

    render: function(){
        console.log('render listView');

        this.reset_grid(this.options);

        return this;
    }
  });

  return ListView;
});