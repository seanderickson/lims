define([
  'jquery',
  'underscore',
  'backbone',
  'views/menu',
  'views/list',
], function($, _, Backbone, MenuView, ListView){

  var AppRouter = Backbone.Router.extend({

    initialize: function(options){
        console.log('initialize router...' + JSON.stringify(options));
        // since the router isn't a view, it doesn't have a model, creating one here
        this.model = options.model;

        Backbone.history.start();
        _.bindAll(this,'change_route'); // binds all of the objects function properties to this instance
        this.model.bind('change:route', this.change_route );
        console.log('router initialized...');
    },

    change_route: function(){
        var newRoute = this.model.get('route');
        console.log('change route to: ' + newRoute);
        this.navigate( newRoute );
    },

    routes: {
        '': 'index',
        // 'list/:type': 'toList',
        // 'list/:type(/page/:page)': 'toPage',
        // 'rpp/:rpp(/page/:page)': 'toRowsPerPage',
        // 'order_by/:orderBy(/rpp/:rpp)(/page/:page)': 'toOrderedToPage',
        'list/:type(/search/:searchBy)(/order_by/:orderBy)(/rpp/:rpp)(/page/:page)': 'toSearchOrderedToPage',
        '*unknownAction': 'unknownAction',
    },

    unknownAction: function(unknownAction){
        alert('Unknown action entered: ' + unknownAction);
    },

    index: function(){
        console.log("Index route has been called..");
        // TODO: now set the app_state.menu_item, and presumably it will fire menu change event
        this.model.set({ menu_item:'home'});
    },

    toSearchOrderedToPage: function(type,searchBy, orderBy,rpp, page){
        console.log("toSearchOrderedToPage route: searchBy: " + searchBy
            + ", order: "+  orderBy + ", rpp: " + rpp + ", page: " + page + ', type: ' + type);

        options = { router: this, type: type };

        options.page = null;
        if(typeof page !== 'undefined' && page !== null ){
            console.log('page: ' + page);
            options.page = parseInt(page);
            // // set the state instead, will be caught and used with a refresh                //              App._collection.getPage(page);
            // App._collection.state.currentPage = parseInt(page);
        }

        options.pageSize = null;
        if(typeof rpp !== 'undefined' && rpp !== null ){
            console.log('rpp: ' + rpp);
            options.pageSize = parseInt(rpp);
            // set the state instead, will be caught and used with a refresh
            //              App._collection.getPage(page);
            // App._collection.state.pageSize = parseInt(rpp);
            // // TODO: catch when pageSize is too large
        }

        options.sortKey = null;
        if(typeof orderBy !== 'undefined' && orderBy !== null ){
            var direction = 'ascending';
            var order = -1;
            if(orderBy.charAt(0) === '-'){
                order = 1; // according to the docs, 1 == descending
                direction = 'descending';
                orderBy = orderBy.substring(1);
            }
            options.sortKey = orderBy;
            options.order = order;
            options.direction = direction;

            // App._collection.state.sortKey = orderBy;
            // App._collection.state.order = order;
            // // Notify header cells
            // App._collection.trigger("backgrid:sort", orderBy, direction, null, App._collection);
        }

        options.searchBy = null;
        if(typeof searchBy !== 'undefined' && searchBy !== null){
            // TODO: only can search one term at a time
            var p = /([^=]+)=([^=]+)/
            var match = p.exec(searchBy);
            if (match){
                // data[match[1] + '__contains'] = match[2];
                // console.log('parsed search: ' + JSON.stringify(data));
                var searchColumn = match[1];
                var searchTerm = match[2];
                console.log('parsed search: ' + searchColumn + ', ' + searchTerm );
                // notify search listeners
                options.searchBy = searchBy;
                options.searchColumn = searchColumn;
                options.searchTerm = searchTerm;

                // App._collection.searchBy = searchBy;
                // App._collection.trigger("MyServerSideFilter:search", searchColumn, searchTerm, App._collection);
             }
        }

        // TODO: switch list with the different types
        if (type ==='fieldinformation'){
            options.url_schema = '/reports/api/v1/fieldinformation/schema/'; // TODO: how to use django url tag here
            options.url = '/reports/api/v1/fieldinformation/?format=json'; // TODO: how to use django url tag here
        }else if (type ==='screensaveruser'){
            options.url_schema = '/db/api/v1/screensaveruser/schema/'; // TODO: how to use django url tag here
            options.url = '/db/api/v1/screensaveruser/?format=json'; // TODO: how to use django url tag here
        }else{
            window.alert('unknown type: ' + type);
        }
        ListView.initialize( options ); // TODO: this moves to the app_state, and then set the menu at the same time
    },

  });

  return AppRouter;
});