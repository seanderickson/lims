define([
  'jquery',
  'underscore',
  'backbone',
  'views/menu',
  'views/list',
], function($, _, Backbone, MenuView, ListView){

  var AppRouter = Backbone.Router.extend({

    routes: {
        '': 'index',
        'list/:type(/search/:searchBy)(/order_by/:orderBy)(/rpp/:rpp)(/page/:page)': 'toList',
        '*unknownAction': 'unknownAction',
    },

    initialize: function(options){
        console.log('initialize router...' );
        // since the router isn't a view, it doesn't have a model, creating one here
        this.model = options.model;

        // _.bindAll(this,'change_route'); // make this available to change_route
        // this.model.bind('change:route', this.change_route );
        this.listenTo(this.model, 'change:route', this.change_route);
        console.log('router initialized...');
    },

    change_route: function(){
        var newRoute = this.model.get('route');
        this.model.set({ content_options: {} }); // unset specific content options
        console.log('change route to: ' + newRoute);
        this.navigate( newRoute, { trigger: false, replace: true } );
    },

    unknownAction: function(unknownAction){
        alert('Unknown action entered: ' + unknownAction);
    },

    index: function(){
        console.log("Index route has been called..");
        this.model.set({ menu_item:'home'});
    },

    toList: function(type,searchBy, orderBy,rpp, page){
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
        options.searchBy = searchBy;

        this.model.set({
            content_options: options,
            current_submodel: options.type
            // TODO: this still feels a little hackish, we're encoding the list/type in the menu item
            // perhaps a controller passed in to the router is a better option
            // right now, this is an aggressive use of the application model change event system
        });
    },
  });

  return AppRouter;
});