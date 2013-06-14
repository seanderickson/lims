define([
  'jquery',
  'underscore',
  'backbone',
  'views/menu',
  'views/list',
], function($, _, Backbone, MenuView, ListView){

  var AppRouter = Backbone.Router.extend({

    initialize: function(options){
        console.log('initialize router...' );
        // since the router isn't a view, it doesn't have a model, creating one here
        this.model = options.model;

        _.bindAll(this,'change_route'); // binds all of the objects function properties to this instance
        this.model.bind('change:route', this.change_route );
        console.log('router initialized...');
    },

    change_route: function(){
        var newRoute = this.model.get('route');
        this.model.set({ content_options: {} }); // unset specific content options
        console.log('change route to: ' + newRoute);
        this.navigate( newRoute, { trigger: false, replace: true } );
    },

    routes: {
        '': 'index',
        // 'list/:type': 'toList',
        // 'list/:type(/page/:page)': 'toPage',
        // 'rpp/:rpp(/page/:page)': 'toRowsPerPage',
        // 'order_by/:orderBy(/rpp/:rpp)(/page/:page)': 'toOrderedToPage',
        'list/:type(/search/:searchBy)(/order_by/:orderBy)(/rpp/:rpp)(/page/:page)': 'toListSearchOrderedToPage',
        '*unknownAction': 'unknownAction',
    },

    unknownAction: function(unknownAction){
        alert('Unknown action entered: ' + unknownAction);
    },

    index: function(){
        console.log("Index route has been called..");
        this.model.set({ menu_item:'home'});
    },

    toListSearchOrderedToPage: function(type,searchBy, orderBy,rpp, page){
        console.log("toSearchOrderedToPage route: searchBy: " + searchBy
            + ", order: "+  orderBy + ", rpp: " + rpp + ", page: " + page + ', type: ' + type);

        options = { type: type };

        options.page = null;
        if(typeof page !== 'undefined' && page !== null ){
            options.page = parseInt(page);
        }

        options.pageSize = null;
        if(typeof rpp !== 'undefined' && rpp !== null ){
            options.pageSize = parseInt(rpp);
        }

        options.orderBy = orderBy;
        // if(typeof orderBy !== 'undefined' && orderBy !== null ){
            // var direction = 'ascending';
            // var order = -1;
            // if(orderBy.charAt(0) === '-'){
                // order = 1; // according to the docs, 1 == descending
                // direction = 'descending';
                // orderBy = orderBy.substring(1);
            // }
            // options.sortKey = orderBy;
            // options.order = order;
            // options.direction = direction;
        // }

        options.searchBy = searchBy;
        // if(typeof searchBy !== 'undefined' && searchBy !== null){
            // // TODO: only can search one term at a time
            // var p = /([^=]+)=([^=]+)/
            // var match = p.exec(searchBy);
            // if (match){
                // // data[match[1] + '__contains'] = match[2];
                // // console.log('parsed search: ' + JSON.stringify(data));
                // var searchColumn = match[1];
                // var searchTerm = match[2];
                // console.log('parsed search: ' + searchColumn + ', ' + searchTerm );
                // // notify search listeners
                // options.searchBy = searchBy;
                // options.searchColumn = searchColumn;
                // options.searchTerm = searchTerm;
             // }
        // }

        this.model.set({
            content_options: options,
            menu_item: 'list_' + options.type
            // TODO: this still feels a little hackish, we're encoding the list/type in the menu item
            // perhaps a controller passed in to the router is a better option
            // right now, this is an aggressive use of the application model change event system
        });

    },

  });

  return AppRouter;
});