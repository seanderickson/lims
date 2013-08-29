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
        'home(/home)': 'toHome',
        'list/:ui_resource_id(/rpp/:rpp)(/page/:page)(/order_by/:orderBy)(/search/:searchBy)': 'toList',
        'list/:ui_resource_id(/search/:searchBy)(/order_by/:orderBy)(/page/:page)(/rpp/:rpp)': 'toList1',
        'list/:ui_resource_id(/search/:searchBy)(/order_by/:orderBy)(/rpp/:rpp)(/page/:page)': 'toList2',
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
        this.listenTo(this.model, 'change:current_resource_id', this.model_set_route);
        this.listenTo(this.model, 'change:current_options', this.model_set_route);

        console.log('router initialized...');
    },

    back: function() {  // TODO: not used yet, from example, how to have a safe back action
        if(this.routesHit > 1) {
          console.log('back, routesHit: ' + this.routesHit);
          //more than one route hit -> user did not land to current page directly
          this.routesHit--;
          window.history.back();
        } else {
          console.log('first route in site, back to home...');
          //otherwise go to the home page. Use replaceState if available so
          //the navigation doesn't create an extra history entry
          this.navigate('/', {trigger:true, replace:true});
        }
      },

    unknownAction: function(unknownAction){
        alert('Unknown action entered: ' + unknownAction);
    },

    get_route: function(){
        var current_view = this.model.get('current_view');
        Iccbl.assert( !_.isUndefined(current_view), 'current_view is not defined');

        var current_resource_id = this.model.get('current_resource_id');
        Iccbl.assert( !_.isUndefined(current_resource_id), 'current_resource_id is not defined');

        var current_options = this.model.get('current_options');
        Iccbl.assert( !_.isUndefined(current_options), 'router.get_route: current_options');

        var _route = current_view + '/' + current_resource_id;

        console.log('getting route: ' + current_resource_id + ', ' + JSON.stringify(current_options));
        var route_fragment = '';

        if( current_view == 'list' && ! _.isEmpty(current_options) ) { // in this case, parse the set of list options in the order needed for the route parsing
            var default_arg_order = ['rpp', 'page','order_by','search'];
            _.each(default_arg_order, function(option){
                if(_.has(current_options, option)){
                    route_fragment += '/' + option + '/' + current_options[option];
                }
            });
        }else{
            if(_.isString(current_options)){
                route_fragment += '/' + current_options;
            }else if(_.isArray(current_options)){
                route_fragment = '/' + _.reduce(current_options, function(route, option){
                    if(!_.isNull(option)) route += option + '/';
                    return  route; // use a url-encoded slash between keys
                    }, route_fragment );
            }else if(_.isObject(current_options)){  // generic, option order not defined
                route_fragment = _.reduce(_.pairs(current_options), function(route, pair){
                    return route + '/' + pair[0] + '/' + pair[1];
                    }, route_fragment );
            }
        }


        _route += route_fragment;
        //console.log('get_route: ' + _route);
        return _route;
    },

    model_set_route: function(){
        // trigger false to suppress further parsing, replace false (default) to create browser history
        var options = { trigger: false }; // , replace:false
        var routing_options = this.model.get('routing_options');
        this.model.set({ routing_options: {} });
        if(!_.isUndefined(routing_options)){
            options = _.extend(options, routing_options);
        }

        var route = this.get_route();
        console.log('--- model_set_route: ' + route + ', ' + JSON.stringify(options) + ', ' + JSON.stringify(routing_options) );
        this.navigate( route, options );
    },

    // model_update_route: function(){
        // var current_route_update = this.model.get('current_route_update');
        // Iccbl.assert( !_.isUndefined(current_route_update), 'current_route_update');
        // console.log('model_update_route: ' + this.get_route() + " , " +  current_route_update);
//
        // // trigger false to suppress further parsing, replace to modify in place w/out creating browser history
        // var options = { trigger: false, replace: true };
        // console.log('update route: ' + current_route_update );
        // this.navigate( this.get_route() + '/' + current_route_update, options );
    // },

    index: function(){
        console.log("Index route has been called..");
        this.model.set({ menu_item:'home', view: 'home' });
    },

    toList1: function(ui_resource_id,searchBy, orderBy,page, rpp){  // kinda crappy that they can't figure this out in either order
        this.toList(ui_resource_id, rpp, page, orderBy, searchBy);
    },

    toList2: function(ui_resource_id,searchBy, orderBy,rpp, page){
    },

    toList: function(ui_resource_id, rpp, page, orderBy, searchBy ){
        console.log("toList: searchBy: " + searchBy
            + ", order: "+  orderBy + ", rpp: " + rpp + ", page: " + page + ', ui_resource_id: ' + ui_resource_id);

        //var _content_options = { ui_resource_id: ui_resource_id, view: 'list' };
        var _content_options = {};

        if( _.isString(page)){
            _content_options.page = parseInt(page);
        }

        if( _.isString(rpp) ){
            _content_options.rpp = parseInt(rpp);
        }


        if( _.isString(orderBy) ){
            _content_options.order = orderBy;
        }

        if( _.isString(searchBy)){
            _content_options.search = searchBy;
        }

        this.model.set({ current_view: {}, current_options: {} }, {silent:true});

        // console.log('toList model.set: ' + ui_resource_id + ', ' + JSON.stringify(_content_options));
        this.model.set({
            current_options: _content_options,
            current_view: 'list',
            current_resource_id: ui_resource_id,
            routing_options: { trigger: false, replace: true } // TODO: necessary?
        });
    },


    toDetail: function(ui_resource_id, key, key2){
        console.log('to detail page, ui_resource_id: ' + ui_resource_id + ', ' + key + ', ' + key2);
        var _current_options = key;
        if(!_.isUndefined(key2)){
            _current_options = [key,key2]; // allow for composite ids
        }
        this.model.set({
            current_view: 'detail',
            current_resource_id: ui_resource_id,
            current_options: _current_options,
            routing_options: { trigger: false, replace: true } // TODO: necessary?
        });
    },

    toHome: function(){
        this.model.set({
            current_view: 'home',
            current_resource_id: 'home',
            current_options: {}
        });
    },

    // toList1: function(ui_resource_id,searchBy, orderBy,rpp, page){
        // console.log("toList: searchBy: " + searchBy
            // + ", order: "+  orderBy + ", rpp: " + rpp + ", page: " + page + ', ui_resource_id: ' + ui_resource_id);
//
        // var _content_options = { ui_resource_id: ui_resource_id, 'view': 'list' };
//
        // var _route_options = {};
        // _content_options.page = null;
        // if(typeof page !== 'undefined' && page !== null ){
            // _content_options.page = parseInt(page);
            // _route_options.page = parseInt(page);
        // }
//
        // _content_options.pageSize = null;
        // if(typeof rpp !== 'undefined' && rpp !== null ){
            // _content_options.pageSize = parseInt(rpp);
            // _route_options.pageSize = parseInt(rpp);
        // }
//
//
        // if( _.isString(orderBy) ){
            // _content_options.orderBy = orderBy;
            // _route_options.orderBy = orderBy;
        // }
//
        // if( _.isString(searchBy)){
            // _content_options.searchBy = searchBy;
            // _route_options.searchBy = searchBy;
        // }
//
        // // TODO: all these content options should really be going into the current_route_options,
        // // and the list view should understand the mappings (i.e. rpp->pageSize)
//
        // this.model.set({ current_view: {}, current_route_options: {} }, {silent:true});
//
        // console.log('toList model.set: ' + JSON.stringify(_route_options));
        // this.model.set({
            // content_options: _content_options,
            // current_view: 'list',
            // current_resource_id: _content_options.ui_resource_id,
            // current_route_options: _route_options
            // // TODO: this still feels a little hackish, we're encoding the list/ui_resource_id in the menu item
            // // perhaps a controller passed in to the router is a better option
            // // right now, this is an aggressive use of the application model change event system
        // });
    // },

  });

  return AppRouter;
});