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
        'list/:ui_resource_id(/search/:searchBy)(/order_by/:orderBy)(/rpp/:rpp)(/page/:page)': 'toList',
        'detail/:ui_resource_id/:key(/:key2)(/)': 'toDetail',
        '*unknownAction': 'unknownAction',
    },

    initialize: function(attributes, options){
        console.log('initialize router...' );
        // since the router isn't a view, it doesn't have a model, creating one here
        this.model = attributes.model;
        this.routesHit = 0;

        // _.bindAll(this,'change_route'); // make this available to change_route
        // this.model.bind('change:route', this.change_route );
        // this.listenTo(this.model, 'change:route', this.change_route);

        Backbone.history.on('route', function(router, route, params) {
                this.routesHit++;
                //console.log('detected route: ' + route + ', params: ' + JSON.stringify(params) + ', routesHit:' + this.routesHit);
             }, this);

        this.listenTo(this.model, 'change:current_view', this.model_set_route);
        this.listenTo(this.model, 'change:current_ui_resource_id', this.model_set_route);
        this.listenTo(this.model, 'change:current_route_options', this.model_set_route);
        this.listenTo(this.model, 'change:current_route_update', this.model_update_route);


        console.log('router initialized...');
    },

    back: function() {  // TODO: not used yet, from example, how to have a safe back action
        if(this.routesHit > 1) {
          //more than one route hit -> user did not land to current page directly
          this.routesHit--;
          window.history.back();
        } else {
          //otherwise go to the home page. Use replaceState if available so
          //the navigation doesn't create an extra history entry
          this.navigate('/', {trigger:true, replace:true});
        }
      },

    // change_route: function(){
        // var newRoute = this.model.get('route');
        // this.model.set({ content_options: {} }); // unset specific content options
        // //this.navigate( newRoute, { trigger: false, replace: true } );
        // options = { trigger: false };
        // if(this.model.get('routing_options')){
            // console.log('routing options: ' + JSON.stringify(this.model.get('routing_options')));
            // _.extend(options, this.model.get('routing_options'));
            // console.log('routing options: ' + JSON.stringify(options));
        // }
        // this.navigate( newRoute, options);
        // this.model.set({'routing_options': {}} );
    // },

    unknownAction: function(unknownAction){
        alert('Unknown action entered: ' + unknownAction);
    },

    get_route: function(){
        var current_view = this.model.get('current_view');
        Iccbl.assert( !_.isUndefined(current_view), 'current_view is not defined');
        var current_ui_resource_id = this.model.get('current_ui_resource_id');
        Iccbl.assert( !_.isUndefined(current_ui_resource_id), 'current_ui_resource_id is not defined');
        var _route = current_view + '/' + current_ui_resource_id;
        var current_route_options = this.model.get('current_route_options');
        Iccbl.assert( !_.isUndefined(current_route_options), 'current_route_options');
        if(_.isString(current_route_options)){
            _route += '/' + current_route_options;
        }else if(_.isArray(current_route_options)){
            _route = _.reduce(_.each(current_route_options), function(route, option){
                return route + '/' + option;
                }, _route );
        }else if(_.isArray(current_route_options)){
            _route = _.reduce(_.pairs(current_route_options), function(route, pair){
                return route + '/' + pair[0] + '/' + pair[1];
                }, _route );
        }
        console.log('get_route: ' + _route);
        return _route;
    },

    model_set_route: function(){
        var options = { trigger: false };

        this.navigate( this.get_route(), options );
    },

    model_update_route: function(){
        var current_route_update = this.model.get('current_route_update');
        Iccbl.assert( !_.isUndefined(current_route_update), 'current_route_update');

        var options = { trigger: false, replace: true };
        console.log('update route: ' + current_route_update );
        this.navigate( this.get_route() + '/' + current_route_update, options );
    },

    index: function(){
        console.log("Index route has been called..");
        this.model.set({ menu_item:'home', view: 'home' });
    },

    toList: function(ui_resource_id,searchBy, orderBy,rpp, page){
        console.log("toSearchOrderedToPage route: searchBy: " + searchBy
            + ", order: "+  orderBy + ", rpp: " + rpp + ", page: " + page + ', ui_resource_id: ' + ui_resource_id);

        var _content_options = { ui_resource_id: ui_resource_id, 'view': 'list' };

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

        // TODO: all these content options should really be going into the current_route_options,
        // and the list view should understand the mappings (i.e. rpp->pageSize)

        this.model.set({ current_view: {}, current_route_options: {} }, {silent:true});

        this.model.set({
            content_options: _content_options,
            current_view: 'list',
            current_ui_resource_id: _content_options.ui_resource_id
            // TODO: this still feels a little hackish, we're encoding the list/ui_resource_id in the menu item
            // perhaps a controller passed in to the router is a better option
            // right now, this is an aggressive use of the application model change event system
        });
    },


    toDetail: function(ui_resource_id, key, key2){
        console.log('to detail page, ui_resource_id: ' + ui_resource_id + ', ' + key + ', ' + key2);
        var _content_options = { ui_resource_id: ui_resource_id, key: key, view: 'detail'};
        if(!_.isUndefined(key2)){
            _content_options['key'] = [key,key2]; // allow for composite ids
        }
        // TODO: cleanup; content options are not used by conent view for detail?
        this.model.set({
            content_options: _content_options,
            current_view: 'detail',
            current_ui_resource_id: _content_options.ui_resource_id,
            current_route_options: content_options['key']
        });
    },

  });

  return AppRouter;
});