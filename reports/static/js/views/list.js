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
        $('#loading').fadeIn({duration:100});
    };
    var ajaxComplete = function(){
        $('#loading').fadeOut({duration:100});
    };
    $(document).bind("ajaxComplete", function(){
        ajaxComplete(); // TODO: bind this closer to the collection
    });

    var ListView = Backbone.View.extend({

        initialize : function() {
                console.log('initialize ListView');
        },

        setOptions : function(options){
            this.options = options;
        },

        remove: function(){
            console.log('ListView remove called');
            Backbone.View.prototype.remove.apply(this);
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

            this.objects_to_destroy = _([]);
            // TODO: make sure req'd options are present (i.e. router, type)
            var collection = this.collection = new App.MyCollection({ 'url': _url, router: options.router, type: options.type });
            this.objects_to_destroy.push(collection);

            var columns = App.createBackgridColModel(fieldDefinitions, App.MyHeaderCell, collection);

            var grid = this.grid = new Backgrid.Grid({
              columns: columns,
              collection: collection
            });
            this.objects_to_destroy.push(grid);

            var selector = new App.ItemsPerPageSelector(
                { 'selections': ['25','50','200','1000'], 'template': rowsPerPageTemplate },
                collection);
            this.objects_to_destroy.push(selector);

            var paginator = new Backgrid.Extension.Paginator({
              collection: collection
            });
            this.objects_to_destroy.push(paginator);


            // var allEvent = function(event, options){
                // console.log("event fired: " + event);
            // };
            // collection.on({
              // "all": allEvent,  // TODO: this is debug code
            // });

            var data = { message: '' };
            if (this.options.header_message){
                data.message = this.options.header_message; //'hello world!' };
            }
            var compiledTemplate = _.template( listTemplate, data );
            this.el.innerHTML = compiledTemplate;
            $("#paginator-div").append(paginator.render().$el);
            $("#rows-selector-div").append(selector.render().$el);
            $("#example-table").append(grid.render().$el);

            // TODO: tie these events into the collection override?
            //        collection.bind('request', window.App.ajaxStart, this);
            //collection.on('request', ajaxStart); // NOTE: can use bind or on
            this.listenTo(collection, 'request', ajaxStart);
            this.listenTo(collection, 'MyCollection:setRoute', this.setRoute);
            //collection.bind('sync', ajaxComplete, this);
            this.listenTo(collection, 'error', ajaxComplete); // TODO: not tested!
            this.listenTo(collection, 'complete', ajaxComplete); // TODO: not tested!

            //collection.bind('sync', selector.updateSelection, selector);  // TODO: selector.listenTo(collection, 'sync'...
            this.listenTo(collection, 'sync', selector.updateSelection );

            if (options.page){
                collection.state.currentPage = options.page;
            }
            if (options.pageSize){
                collection.state.pageSize = options.pageSize;
            }
            // if (options.sortKey){
                // collection.state.sortKey = options.sortKey;
                // collection.state.order = options.order;
                // // Notify header cells
                // collection.trigger("backgrid:sort", options.sortKey, options.direction, null, collection);
            // }
            if(typeof options.orderBy !== 'undefined' && options.orderBy !== null ){
                var direction = 'ascending';
                var order = -1; // according to the docs, -1 == ascending
                var sortKey = options.orderBy;
                if(sortKey.charAt(0) === '-'){
                    order = 1; // according to the docs, 1 == descending
                    direction = 'descending';
                    sortKey = sortKey.substring(1);
                }
                collection.state.sortKey = sortKey;
                collection.state.order = order;
                // Notify header cells
                collection.trigger("backgrid:sort", sortKey, direction, null, collection);
            }
            if(typeof options.searchBy !== 'undefined' && options.searchBy !== null){
                // TODO: only can search one term at a time
                var p = /([^=]+)=([^=]+)/
                var match = p.exec(options.searchBy);
                if (match){
                    // data[match[1] + '__contains'] = match[2];
                    // console.log('parsed search: ' + JSON.stringify(data));
                    var searchColumn = match[1];
                    var searchTerm = match[2];
                    console.log('parsed search: ' + searchColumn + ', ' + searchTerm );

                    collection.searchBy = options.searchBy;
                    console.log('trigger search');
                    collection.trigger("MyServerSideFilter:search", searchColumn, searchTerm, collection);
                    console.log('done: trigger search');
                 }else{
                     console.log('unrecognized searchBy: ' + options.searchBy);
                 }
            }

            collection.fetch({ reset: true } );

            console.log('list view initialized');
        },

        setRoute: function(route){
            console.log('setRoute triggered: ' + route);
            this.model.set({ route: 'list/' + this.options.type + '/' + route } );
        },

        onClose: function(){
            console.log('Extra onclose method called');
            this.objects_to_destroy.each(function(view_obj){
                view_obj.remove();
                view_obj.unbind();
                view_obj.stopListening();
            });
        },

        render: function(){ // TODO: should the build grid be called from here?
            console.log('render listView');
            this.reset_grid(this.options);
            return this;
        }
    });

  return ListView;
});