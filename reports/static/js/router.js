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
        'detail/:type/:id(/:id2)(/)': 'toDetail',
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
        console.log('------change route to: ' + newRoute );
        //this.navigate( newRoute, { trigger: false, replace: true } );
        options = { trigger: false };
        if(this.model.get('routing_options')){
            console.log('routing options: ' + JSON.stringify(this.model.get('routing_options')));
            _.extend(options, this.model.get('routing_options'));
            console.log('routing options: ' + JSON.stringify(options));
        }
        this.navigate( newRoute, options);
        this.model.set({'routing_options': {}} );
    },

    unknownAction: function(unknownAction){
        alert('Unknown action entered: ' + unknownAction);
    },

    index: function(){
        console.log("Index route has been called..");
        this.model.set({ menu_item:'home', view: 'home_view' });
    },

    toList: function(type,searchBy, orderBy,rpp, page){
        console.log("toSearchOrderedToPage route: searchBy: " + searchBy
            + ", order: "+  orderBy + ", rpp: " + rpp + ", page: " + page + ', type: ' + type);

        var _content_options = { type: type, 'view': 'list_view' };

        _content_options.page = null;
        if(typeof page !== 'undefined' && page !== null ){
            _content_options.page = parseInt(page);
        }

        _content_options.pageSize = null;
        if(typeof rpp !== 'undefined' && rpp !== null ){
            _content_options.pageSize = parseInt(rpp);
        }

        _content_options.orderBy = orderBy;
        _content_options.searchBy = searchBy;

        this.model.set({
            content_options: _content_options,
            current_submodel_id: _content_options.type
            // TODO: this still feels a little hackish, we're encoding the list/type in the menu item
            // perhaps a controller passed in to the router is a better option
            // right now, this is an aggressive use of the application model change event system
        });
    },

    toDetail: function(type, id, id2){
        console.log('to detail page, type: ' + type + ', ' + id + ', ' + id2);
        var _content_options = { type: type, id: id, view: 'detail_view'};
        if(!_.isUndefined(id2)){
            _content_options['id'] = [id,id2]; // allow for composite ids
        }
        this.model.set({
            content_options: _content_options,
            current_submodel_id: _content_options.type
        });
    },

  });

  return AppRouter;
});