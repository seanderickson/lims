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
        'list/:ui_resource_id(/rpp/:rpp)(/page/:page)(/order/:orderBy)(/search/:searchBy)(/)': 'toList', // note: can tolerate missing routes, but not out of order
        'detail/:ui_resource_id/:key(/:key2)(/tab/:tab)(/rpp/:rpp)(/page/:page)(/order/:orderBy)(/search/:searchBy)(/)': 'toDetail',
        'menu/:action' : 'toMenu',
        '*unknownAction': 'unknownAction',
    },
    LIST_ROUTE_ORDER: ['rpp', 'page','order','search'],

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
    get_list_route: function(current_options){
        console.log('get list route: ' + JSON.stringify(current_options));
        var route_fragment = '';
        _.each(this.LIST_ROUTE_ORDER, function(option){

            if(_.has(current_options, option)){
                var optionValue = current_options[option]
                if(!_.isUndefined(optionValue) && !_.isNull(optionValue)){
                    if(route_fragment.length > 0) route_fragment += '/';
                    route_fragment += option + '/';
                    if(option == 'rpp' || option == 'page' ){
                        route_fragment += optionValue;

                    }else if(option == 'order' ){
                        // TODO: use the type of search comparator as well
                        var frag = _.reduce(_.pairs(optionValue), function(frag, pair){
                            if(frag.length > 0 ) frag += ',';
                            if(pair[1] == '-' ) frag += '-';
                            return frag + pair[0];
                        }, '' );
                        route_fragment += frag;
                    }else if(option == 'search'){
                        // TODO: use the type of search comparator as well
                        var frag = _.reduce(_.pairs(optionValue), function(frag, pair){
                            if(frag.length > 0 ) frag += ',';

                            // // // TODO: discuss whether full encoding of the search fragment is necessary.
                            // // // to-date, we know that the forward slash messes up the backbone router parsing, but other URL chars do not,
                            // // // and full encoding reduces the usability of the URL for the end user
                            // // //_history_search = encodeURIComponent(_history_search);
                            // // new_options['search'] = this.searchBy.replace('\/','%2F');
                            var value = pair[1].replace('\/','%2F');
                            return frag + pair[0] + '=' + value;
                        }, '' );
                        route_fragment += frag;
                    }
                }
            }
        });
        return route_fragment;
    },

    get_route: function(){
        var current_view = this.model.get('current_view');
        Iccbl.assert( !_.isUndefined(current_view), 'router: current_view is not defined');

        var current_resource_id = this.model.get('current_resource_id');
        Iccbl.assert( !_.isUndefined(current_resource_id), 'router: current_resource_id is not defined');

        var current_options = this.model.get('current_options');
        Iccbl.assert( !_.isUndefined(current_options), 'router.get_route: current_options');

        var _route = current_view + '/' + current_resource_id;

        console.log('getting route: ' + current_resource_id + ', ' + JSON.stringify(current_options));
        var route_fragment = '';

        if( ! _.isEmpty(current_options) ) {
            if( current_view == 'list' ) { // in this case, parse the set of list options in the order needed for the route parsing
                route_fragment += '/' + this.get_list_route(current_options);

            }else if( current_view == 'detail' ){
                route_fragment = '/';
                var key = Iccbl.getKey(current_options);
                if(!_.isEmpty(key)) route_fragment += key;
                if(key.charAt(key.length-1) != '/' ) route_fragment += '/'; // TODO: cleanup code so no trailing slashes
                if(!_.isEmpty(current_options['tab'])) route_fragment += 'tab/' + current_options.tab;


                var list_route = this.get_list_route(current_options);
                if(!_.isEmpty(list_route)) route_fragment += '/';
                route_fragment += list_route;
                // _.each(this.LIST_ROUTE_ORDER, function(option){
                    // if(_.has(current_options, option)){
                        // route_fragment += '/' + option + '/' + current_options[option];
                    // }
                // });
            }else{
                if(_.isString(current_options)){
                    route_fragment += '/' + current_options;
                }else if(_.isArray(current_options)){ // an array is just a set of keys, to be separated by slashes
                    route_fragment = '/' + Iccbl.getKey(current_options);
                }else if(_.isObject(current_options)){  // generic, option order not defined
                    route_fragment = _.reduce(_.pairs(current_options), function(route, pair){

                        return route + '/' + pair[0] + '/' + pair[1];
                    }, route_fragment );
                }
            }
        }

        _route += route_fragment;
        //console.log('get_route: ' + _route);
        return _route;
    },

    model_set_route: function(){
        console.log('model_set_route');
        // trigger false to suppress further parsing, replace false (default) to create browser history
        var options = { trigger: false }; // , replace:false
        var routing_options = this.model.get('routing_options');
        this.model.set({ routing_options: {} });
        if(!_.isUndefined(routing_options)){
            options = _.extend(options, routing_options);
        }

        var route = this.get_route();
        console.log('--- model_set_route: ' + route + ', ' + JSON.stringify(options) );
        this.navigate( route, options );
    },

    index: function(){
        console.log("Index route has been called..");
        this.model.set({ menu_item:'home', view: 'home' });
    },

    toList: function(ui_resource_id, rpp, page, orderBy, searchBy ){
        var self = this;
        console.log("toList: searchBy: " + searchBy
            + ", order: "+  orderBy + ", rpp: " + rpp + ", page: " + page + ', ui_resource_id: ' + ui_resource_id);

        var _content_options = {};

        if( _.isString(page)){
            _content_options.page = parseInt(page);
        }

        if( _.isString(rpp) ){
            _content_options.rpp = parseInt(rpp);
        }


        if( _.isString(orderBy) ){
            _content_options.order = self.parseOrder(orderBy);
        }

        if( _.isString(searchBy)){

            self.parseSearch(searchBy, {
                success: function(search_full_string_encoded, searchExpressionHash){
                    _content_options.search = searchExpressionHash;
                },
                error: function(search_full_string_encoded, msg){
                    console.log("ERROR: " + search_full_string_encoded + ', parse failed with message: ' + msg);
                }
            });
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

    // Custom search param parsing method
    // successCallback: function(search_full_string_encoded,searchExpressionsHash)
    // errorFunction: function( search_full_string_encoded, errorMessage);
    parseSearch: function(search_full_string_encoded, options){
        var searchExpressions = {};
        // TODO: query terms that tastypie will understand.  these are to be set on the MyHeaderCell
        // QUERY_TERMS = set([
            // 'exact', 'iexact', 'contains', 'icontains', 'gt', 'gte', 'lt', 'lte', 'in',
            // 'startswith', 'istartswith', 'endswith', 'iendswith', 'range', 'year',
            // 'month', 'day', 'week_day', 'isnull', 'search', 'regex', 'iregex',
        // ])

        _(search_full_string_encoded.split(',')).each(function(searchItem){
            var searchExpression = searchItem.split('=');
            if(searchExpression.length != 2 ){
                if(_.has(options, 'error')){
                    options.error(search_full_string_encoded, 'Warning: invalid search item: ' + searchItem );
                }else{
                    console.log('Error: invalid search item: ' + searchItem + ', in: ' + search_full_string_encoded);
                }
                return;
            }else{
                searchExpressions[searchExpression[0]] = searchExpression[1];
            }
        });

        if(!_.isEmpty(searchExpressions)){
            if(_.has(options, 'success')){
                options.success(search_full_string_encoded, searchExpressions);
            }else{
                return searchExpressions;
            }
        }else{
            var msg = 'Error: no search expressions found';
            if(_.has(options, 'error')){
                errorFunction(search_full_string_encoded, msg );
            }else{
                console.log(msg + ', in: ' + search_full_string_encoded);
            }
        }
    },

    parseOrder: function(orderString){
        var orderHash = {};
        if(!_.isEmpty(orderString)){
            _.each(orderString.split(','), function(item){
                var key = item;
                var dir = '+';
                if(item.charAt(0) == '-' ) {
                    key = item.substring(1);
                    dir = '-';
                }
                orderHash[key]=dir;
           });
        }
        return orderHash;
    },


    toDetail: function(ui_resource_id, key, key2, _tab, rpp, page, orderBy, searchBy){
        console.log('to detail page, ui_resource_id: ' + ui_resource_id + ', ' + key + ', ' + key2 + ', ' + rpp + ', ' + page + ', ' + orderBy + ', ' + searchBy);
        var _current_options = { 'key': key };
        if(!_.isUndefined(key2) && !_.isNull(key2) ){
            _current_options['key'] = [key,key2]; // allow for composite ids
        }
        if(!_.isEmpty(_tab)){
            _current_options['tab'] = _tab;
        }
        if( _.isString(page)){
            _current_options.page = parseInt(page);
        }

        if( _.isString(rpp) ){
            _current_options.rpp = parseInt(rpp);
        }


        if( _.isString(orderBy) ){
            _current_options.order = this.parseOrder(orderBy);
        }

        if( _.isString(searchBy)){
            this.parseSearch(searchBy, {
                success: function(search_full_string_encoded, searchExpressionHash){
                    _content_options.search = searchExpressionHash;
                },
                error: function(search_full_string_encoded, msg){
                    console.log("ERROR: " + search_full_string_encoded + ', parse failed with message: ' + msg);
                }
            });
        }

        this.model.set({
            current_view: 'detail',
            current_resource_id: ui_resource_id,
            current_options: _current_options,
            current_detail: _current_options,
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

    toMenu: function(action){
        this.model.set({
            current_view: 'menu',
            current_resource_id: action,
            current_options: {}
        });
    },


  });

  return AppRouter;
});