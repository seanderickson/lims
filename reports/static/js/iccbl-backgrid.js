/**
 * @summary     ICCBL-lims Backgrid Extension Functions
 * @description Utility Functions for the iccbl-lims
 * @version     0.1
 * @file        iccbl-backgrid.js
 * @author      Sean Erickson
 * @contact     sean_erickson “AT” hms.harvard.edu
 *
 * @copyright Copyright 2013 Harvard University, all rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This source file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license files for details.
 *
 * For details please refer to: https://github.com/hmsiccbl/lims
 *
 **/

define([
  'jquery',
  'underscore',
  'backbone',
  'backbone_pageable',
  'backgrid',
  'backgrid_filter',
  'backgrid_paginator',
  'lunr', //TODO: lunr should not be a requirement - it is for client side filtering, and backgrid_filter is requiring it.
  'text!templates/generic-selector.html'
], function($, _, Backbone, BackbonePageableCollection, Backgrid, BackgridFilter, BackgridPaginator, lunr, genericSelectorTemplate){

    // for compatibility with require.js, attach PageableCollection in the right place on the Backbone object
    // see https://github.com/wyuenho/backbone-pageable/issues/62
    Backbone.PageableCollection = BackbonePageableCollection;

    var root = window;
    // create a root namespace (note: copying this strategy from Backgrid source
    var Iccbl = root.Iccbl = {

      VERSION: "0.0.1",

      // Extension: {},
//
      // requireOptions: function (options, requireOptionKeys) {
        // for (var i = 0; i < requireOptionKeys.length; i++) {
          // var key = requireOptionKeys[i];
          // if (_.isUndefined(options[key])) {
            // throw new TypeError("'" + key  + "' is required");
          // }
        // }
      // },
//
      // resolveNameToClass: function (name) {
        // if (_.isString(name)) {
          // // var key = _.map(name.split('-'), function (e) { return capitalize(e); }).join('') + suffix;
          // var klass = Iccbl[key]; // || Backgrid.Extension[key];
          // if (_.isUndefined(klass)) {
            // throw new ReferenceError("Class '" + key + "' not found");
          // }
          // return klass;
        // }
//
        // return name;
      // }
    };
    //var Iccbl = {};

    var assertIccbl = Iccbl.assert = function (condition, message) {
        if (!condition) {
            throw message || "Assertion failed";
        }
    };

    var stringToFunction = Iccbl.stringToFunction = function(str) {
      var arr = str.split(".");

      var fn = (window || this);
      for (var i = 0, len = arr.length; i < len; i++) {
        fn = fn[arr[i]];
      }

      if (typeof fn !== "function") {
        throw new ReferenceError("function not found: " + str);
      }

      return  fn;
    };

    var sortOnOrdinal = Iccbl.sortOnOrdinal = function(keys, fieldHash){
        var sorted = _(keys).sort(function(a,b){
            console.log('comparing: ' + a + ', to ' + b );
            order_a = fieldHash[a]['ordinal'];  // TODO: need an edit order by
            order_b = fieldHash[b]['ordinal'];
            if(_.isNumber(order_a) && _.isNumber(order_b)){
                return order_a - order_b;
            }else if(_.isNumber(order_a)){
                return -1;
            }else if(_.isNumber(order_b)){
                return 1;
            }else{
                return 0;
            }
        });
        return sorted;
    };

    var getKey = Iccbl.getKey = function( options ){
        var route_fragment = '';
        if(_.isString(options)){
            route_fragment += options;
        }else if( _.isArray(options)) {// an array is just a set of keys, to be separated by slashes
            route_fragment = _.reduce(options, function(route, option){
                    if(!_.isNull(option)){
                        if(!_.isEmpty(route)) route += '/';
                        route += option;
                    }
                    return  route;
                }, route_fragment );
        }else if(_.isObject(options)){  // generic, option order not defined
            if( _.has(options, 'key')){
                return Iccbl.getKey(options.key);
            }
        }

        return route_fragment;
    };

    // var getKey = Iccbl.getKey = function( key_array ){
        // var _key = key_array;
        // Iccbl.assert( !_.isEmpty(_key), 'content:detail: current_options must be defined (as the key), if not schemaResult, model supplied');
        // // handle composite keys
        // if(_.isArray(_key)){
            // _key = _.reduce(_key, function(memo, item){
                // if(!_.isNull(item)) memo += item + '/';
                // return memo;
            // }, '');
        // }
        // return _key;
    // };

    var getSchema = Iccbl.getSchema = function (schema_url, callback) {
        $.ajax({
            type: "GET",
            url: schema_url, //options.url_schema,
            data: "",
            dataType: "json",
            success: function(schemaResult) {
                callback(schemaResult);
            }, // end success outer ajax call
            error: function(x, e) {
                alert(x.readyState + " "+ x.status +" "+ e.msg); // TODO: use error div in Bootstrap
            }
        });
    };

    var getModel = Iccbl.getModel = function(schemaResult, url, callback) {
        var ModelClass = Backbone.Model.extend({url: url, defaults: {} });
        var instance = new ModelClass();
        instance.fetch({
            success: function(model){
                callback(schemaResult, model);
            },
            error: function(model, response, options){
                //console.log('error fetching the model: '+ model + ', response: ' + JSON.stringify(response));
                var msg = 'Error locating resource: ' + url;
                var sep = '\n';
                if(!_.isUndefined(response.status)) msg += sep + response.status;
                if(!_.isUndefined(response.statusText)) msg += sep+ response.statusText;
                if(!_.isEmpty(response.responseText)) msg += sep+ response.responseText;
                window.alert(msg); // TODO: use Bootstrap inscreen alert classed message div
            }
        });
    };

    var getCollection = Iccbl.getCollection = function(schemaResult, url, callback) {
        var CollectionClass = Iccbl.CollectionOnClient.extend({url: url, defaults: {} });
        var instance = new CollectionClass();
        instance.fetch({
            success: function(collection){
                callback(schemaResult, collection);
            },
            error: function(model, response, options){
                //console.log('error fetching the model: '+ model + ', response: ' + JSON.stringify(response));
                var msg = 'Error locating resource: ' + url;
                var sep = '\n';
                if(!_.isUndefined(response.status)) msg += sep + response.status;
                if(!_.isUndefined(response.statusText)) msg += sep+ response.statusText;
                if(!_.isEmpty(response.responseText)) msg += sep+ response.responseText;
                window.alert(msg); // TODO: use Bootstrap inscreen alert classed message div
            }
        });
    };

    var getCollection2 = Iccbl.getCollection2 = function(schemaResult, url, callback) {
        var CollectionClass = Iccbl.CollectionInColumns.extend({url: url, defaults: {} });
        var instance = new CollectionClass();
        instance.fetch({
            success: function(collection){
                callback(schemaResult, collection);
            },
            error: function(model, response, options){
                //console.log('error fetching the model: '+ model + ', response: ' + JSON.stringify(response));
                var msg = 'Error locating resource: ' + url;
                var sep = '\n';
                if(!_.isUndefined(response.status)) msg += sep + response.status;
                if(!_.isUndefined(response.statusText)) msg += sep+ response.statusText;
                if(!_.isEmpty(response.responseText)) msg += sep+ response.responseText;
                window.alert(msg); // TODO: use Bootstrap inscreen alert classed message div
            }
        });
    };

    var getSchemaAndModel = Iccbl.getSchemaAndModel = function(schema_url, url, callback){
        Iccbl.getSchema(schema_url, function(schemaResult) {
            console.log('schemaResult callback: ' + schemaResult + ', ' + url);
            Iccbl.getModel(schemaResult, url, callback);
        });
    };


    var getSchemaAndCollection = Iccbl.getSchemaAndCollection = function(schema_url, url, callback){
        Iccbl.getSchema(schema_url, function(schemaResult) {
            console.log('schemaResult callback: ' + schemaResult + ', ' + url);
            Iccbl.getCollection(schemaResult, url, callback);
        });
    };

    var getSchemaAndCollection2 = Iccbl.getSchemaAndCollection2 = function(schema_url, url, callback){
        Iccbl.getSchema(schema_url, function(schemaResult) {
            console.log('schemaResult callback: ' + schemaResult + ', ' + url);
            Iccbl.getCollection2(schemaResult, url, callback);
        });
    };

    var MyModel = Iccbl.MyModel = Backbone.Model.extend({
        // TODO: we want to make sure there is a trailing slash, or tastypie doesn't work.
        url : function(){
            var url = Backbone.Model.prototype.url.call(this);
            // console.log('---- url1: ' + url);
            return url + (url.charAt(url.length - 1) === '/' ? '' : '/') ;
        },

        initialize: function () {
            Backbone.Model.prototype.initialize.apply(this, arguments);
            // console.log('x--- urlRoot: ' + this.urlRoot + ", " + this.id + ', ' + this.collection.url);
            var self = this;
            // we want to make sure there is a trailing slash, or tastypie doesnt work.
            // TODO: not sure why we have to override url function like this
            // this.url = function(){
              // var url = Backbone.Model.prototype.url.call(self);
              // console.log('---- url1: ' + url);
              // return url + (url.charAt(url.length - 1) === '/' ? '' : '/') ;
            // };
            // definition above should work, but doesn't.
            // however, when overriding url like above only, the function has to be attached to the prototype manually here.  why?
            this.url = MyModel.prototype.url;
        },
    });


    var GenericSelector = Iccbl.GenericSelector = Backbone.View.extend({

        initialize: function(attributes, options){

            this.listenTo(this.model, 'change', this.changeNotified);
            this._options = options;
            _.bindAll(this, 'changeNotified');
        },

        events: {
            "change #generic_selector": "selectorChanged"
        },

        selectorChanged: function(e){
            e.preventDefault();
            var option = e.currentTarget.value;
            this.model.set({'selection': option });
        },

        changeNotified: function(){
            var selection = this.model.get('selection');
            console.log('change notified for ' + this._options.label + ', selection: ' + selection + ', options: ' + JSON.stringify(this._options.options) );
            this.$('#generic_selector').val(String(selection));
        },

        render: function(){
            this.$el.empty();
            if(!_.contains(this._options.options, '' )){
                this._options.options.unshift(''); // create a blank entry
            }
            this.$el.append(_.template( genericSelectorTemplate,
                { label: this._options.label,
                  'options': _(this._options.options ) }));  // TODO: this should come from the metahash schema
            if(!_.isUndefined(this._options.selectorClass)){
                this.$('#generic_selector').removeClass().addClass(this._options.selectorClass);
            };
            this.changeNotified();
            this.delegateEvents();
            return this;
        },
    });


    var LinkCell = Iccbl.LinkCell = Backgrid.Cell.extend({
        className: "link-cell",
        events : {
            "click #link" : "toLink",
        },

        initialize: function(options){
            Backgrid.Cell.prototype.initialize.apply(this, arguments);
        },

        render: function () {
            this.$el.empty();
            var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
            this.$el.append($("<a id='link' >", {
                tabIndex : -1,
                href : '',
                title : formattedValue,
                //target : "_blank"
            }).text(formattedValue));

            this.delegateEvents();
            return this;
        },

        toLink: function(e){
            e.preventDefault();
            this.model.collection.trigger("MyCollection:link", this.model, this.column.get("name") );
        },
    });

    var EditCell = Iccbl.EditCell = Backgrid.Cell.extend({
        className: "detail-cell",
        events : {
            "click #edit" : "editDetail",
        },

        initialize: function(options){
            Backgrid.Cell.prototype.initialize.apply(this, arguments);
        },

        render: function () {
            this.$el.empty();
            var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
            this.$el.append($("<a id='edit' >", {
                tabIndex : -1,
                href : '',
                title : formattedValue,
                //target : "_blank"
            }).text(formattedValue));

            this.delegateEvents();
            return this;
        },

        editDetail: function(e){
            e.preventDefault();
            this.model.collection.trigger("MyCollection:edit", this.model);
        },
    });

    /**
     * uses the options.attributes.label
     */
    var DeleteCell = Iccbl.DeleteCell = Backgrid.Cell.extend({
        className: "delete-cell",
        events : {
            "click #delete" : "delete"
        },

        initialize: function(options){
            Backgrid.Cell.prototype.initialize.apply(this, arguments);
        },

        render: function () {
            this.$el.empty();

            this.$el.append("&nbsp;");
            this.$el.append($("<a id='delete' >", {
                tabIndex : -1,
                href : '',
            }).text(this.options.column.attributes['text']));
            this.delegateEvents();
            return this;
        },

        delete: function(e){
            e.preventDefault();
            this.model.collection.trigger("MyCollection:delete", this.model);
        }
    });

    var CollectionOnClient = Iccbl.CollectionOnClient = Backbone.Collection.extend({
        /**
         *  Override collection parse method:
         *      Parse server response data.
         */
        parse: function(response) {
            console.log('Collection on client, parse called');
            // hack the response for tastypie:
            return response.objects;
        },
    });
    var CollectionInColumns = Iccbl.CollectionInColumns = Backbone.Collection.extend({
        /**
         *  Override collection parse method:
         *      Parse server response data.
         * untested
         */
        parse: function(response) {
            console.log('Collection on client, parse called');
            // hack the response for tastypie:
            var pivoted = {};
            var i = 0;
            _.each(response.objects, function(obj){
                _.pairs(obj, function(pair){
                    if(_.has(pivoted,pair[0])){
                        pivoted[pair[0]] = {};
                    }
                    pivoted[pair[0]][i] = pair[1];
                });
                i++;
            });
            return _.values(pivoted);
        },
    });

    var MyCollection = Iccbl.MyCollection = Backbone.PageableCollection.extend({

        initialize: function(options){
            var self = this;
            if (typeof(options) !== 'undefined'){
                //console.log('options: ' + JSON.stringify(options));
            }else{
                window.alert('no options defined');
            }
            // TODO: require these options
            this.url = options.url;

            Iccbl.assert(!_.isUndefined(options.listModel), 'Iccbl.Collection requires options.listModel');

            this.listModel = options.listModel;

            this.listenTo(this.listModel, 'change:search', function(){
                console.log('===--- collection detect: listModel change:search');
                var searchHash = self.listModel.get('search');
                self.setSearch(searchHash);
            });

            this.listenTo(this.listModel, 'change:order', function(){
                console.log('===--- listModel change:order');
                var orderHash = self.listModel.get('order');
                this.setOrder(orderHash);
            });

            this.listenTo(this.listModel, 'change:rpp', function(){
                var pageSize = parseInt(self.listModel.get('rpp'));
                console.log('===--- listModel change:rpp: ' + pageSize);
                self.setPageSize(pageSize);
            });

            this.listenTo(this.listModel, 'change:page', function(){
                console.log('===--- listModel change:page');
                var page = parseInt(self.listModel.get('page'));
                self.setPage(page);
                // self.state.currentPage = page;
                // self.fetch();
            });


            Backbone.PageableCollection.prototype.initialize.apply(this, options); // call super constructor
        },

        url: function() {
            return this.url;
        },

        // searchBy: null,
        searchHash: {},
        model: MyModel,
        state: {
            pageSize: 25,  // TODO: probably not necessary
        },
        queryParams: {
            // adjust the query params for tastypie
            pageSize: 'limit',
            offset: function(){
                return (this.state.currentPage-1) * this.state.pageSize;
            },
            totalRecords: null, // unset for tastypie
            totalPages: null, // unset for tastypie
            sortKey: "order_by", // modified for tastypie
            order: null, // unset for tastypie
            order_by: function() { // modified for tastypie: use only "order_by=(-)field_name"
                if (typeof this.state !== 'undefined' && this.state.sortKey && this.state.order) {
                    var dir = "-";
                    if (this.state.order<0){  // according to docs, -1 == ascending
                        dir = "";
                    }
                    return dir + this.state.sortKey;
                }
            },
            directions: {
                "-1": "asc",
                "1": "desc"
            }

        },

        /**
         *  Override pageable collection parse method:
         *      Parse server response data.
         */
        parse: function(response) {
            // hack the response for tastypie:
            // note, this is because the pageable collection doesn't work with the backbone-tastypie.js fix
            var state = _.clone(this.state);
            state.totalRecords = response.meta.total_count;
            if(Math.ceil(state.totalRecords/state.pageSize) < state.currentPage ){
                console.log('adjust currentPage');
                state.currentPage=1;
            }
            this.state = this._checkState(state); // recalculate the state and do sanity checks.
            //if(! _(state).isEqual(this.state)){
                // TODO: call setOptions explicitly where needed, i.e. when setting searchBy
                console.log('parse: new state: ' + JSON.stringify(this.state));
                // trigger so our app knows about the route change, replace to modify in place w/out creating browser history
                this.setOptions({trigger: false, replace: true });
            //}
            return response.objects;
        },

        /**
           @property {-1|0|1} [state.order=-1] The order to use for sorting. Specify
           -1 for ascending order or 1 for descending order. If 0, no client side
           sorting will be done and the order query parameter will not be sent to
         */
        setOptions: function(route_options) {
            console.log('- setOptions: ' //+ this.searchBy + ', '
                + this.state.sortKey + ', ' + this.state.order + ', '+ this.state.pageSize
                + ', ', + this.state.currentPage + ', options: ' + JSON.stringify(route_options));
            var new_options = {};
            // if(this.searchBy !== null){
//
                // // TODO: discuss whether full encoding of the search fragment is necessary.
                // // to-date, we know that the forward slash messes up the backbone router parsing, but other URL chars do not,
                // // and full encoding reduces the usability of the URL for the end user
                // //_history_search = encodeURIComponent(_history_search);
                // new_options['search'] = this.searchBy.replace('\/','%2F');
                // // route += 'search/'+ encodeURIComponent(this.searchBy) ;
            // }
            // if(typeof this.state.sortKey !== 'undefined' && this.state.sortKey !== null){
                // var sortKey = this.state.sortKey;
                // if(this.state.order > 0){
                    // sortKey = '-' + sortKey;
                // }
                // new_options['order_by'] = sortKey;
            // }

            //if (!this.state.sortKey) this.clearOrder();

            new_options['page'] = this.state.currentPage;
            new_options['rpp'] =  this.state.pageSize;

            // signal to the parent view that the options have changed
            this.trigger( "MyCollection:changeOptions", new_options, route_options); // TODO: need to set replace=false when modifying to add the page, rpp for the first time: { replace: true } );
        },

        // Override
        // called from HeaderCells
        setSorting: function(sortKey,order,options) { // override and hack in sorting URL support
            console.log('setSorting called: ' + sortKey+ ', order_by: ' + order + ', ' + typeof(order) + ', options: ' + options );
            var dir = '-'; // desc
            if(typeof order !== 'undefined' && order < 0){
                dir = '';  // asc
            }
            var obj = Backbone.PageableCollection.prototype.setSorting.call(this,sortKey, order);

            if(!_.isEmpty(sortKey)){
                var orderHash = {};
                orderHash[sortKey] = dir;
                this.addOrder(orderHash);
            }else{
                this.clearOrder();
            }
            return obj;
        },

        // Custom - called from external
        setOrder: function(orderHash){
            console.log('setOrder called: ' + JSON.stringify(orderHash));
            var self = this;

            // TODO: only one order supported on the Pageable collection!!
            _.each(_.keys(orderHash), function(key){
                var dir = orderHash[key];

                var direction = 'ascending';
                var order = -1; // according to the docs, -1 == ascending

                if(dir === '-'){
                    order = 1; // according to the docs, 1 == descending
                    direction = 'descending';
                }

                console.log('setting order: ' + key);
                self.state.sortKey = key;
                self.state.order = order;
                // Notify header cells
                self.trigger("backgrid:sort", key, direction, null, self);
            });
        },

        // Custom search method
        setSearch: function(searchHash){
            var self = this;
            // this.searchBy = search_full_string_encoded;
            var searchHash = _.clone(searchHash);
            console.log('collection.setSearch: trigger search to headers:' + JSON.stringify(searchHash) );
            this.trigger("MyServerSideFilter:search", searchHash, this);
            // console.log('done: trigger search');

            // if the search key is not in the queryParams, then it is not a column search.
            // this will add it manually to the queryParams (which are serialized in the fetch to the server)
            _.each(_.keys(searchHash), function(key){
                // check if param dne, or if param exists and is a value to be set
                // the reason for the "isFunction" check is that the Backgrid-filter defined params are function calls to get the current value in the searchbox
                if(!_.has(self.queryParams, key) || !_.isFunction(self.queryParams[key])){

                    self.queryParams[key] = searchHash[key];  // NOTE: setting the key on the queryParams means that it will be sent to the server on the next fetch
                }
            });
            this.fetch();
        },

        // Custom
        clearSearch: function(searchKeys){
            console.log('clearsearch: ' + JSON.stringify(searchKeys));
            var self = this;
            var searchHash = {};
            if(!_.isUndefined(searchKeys)){
                searchHash = _.clone(self.listModel.get('search'));
                _.each(searchKeys, function(searchKey){
                    delete searchHash[searchKey];
                });
            }
            self.listModel.set({'search': searchHash });

            // this.trigger("MyServerSideFilter:removeSearch", searchKeys, this);
            this.state.currentPage=1;  // if filtering, always set the page to 1

            //this.fetch();
        },

        // Custom
        clearOrder: function(orderKeys){
            var self = this;
            if(!_.isUndefined(orderKeys)){
                var orderHash = _.clone(self.listModel.get('order'));
                _.each(orderKeys, function(orderKey){
                    delete orderHash[orderKey];
                    // self.sort(columnName, null);
                });
                self.listModel.set({'order': orderHash });
            }else{
                self.listModel.set({'order': {}});
            }

            //this.trigger("MyServerSideFilter:removeOrder", orderKeys, this);

            // this.state.currentPage=1;  // if filtering, always set the page to 1
            // this.fetch();
        },


        // Custom
        addOrder: function(orderHash){
            var self = this;
            self.listModel.set({'order': _.clone(orderHash) });
            // TODO: multisort: for the time being, hash only contains one item.  so no need for extend
            // var oldorderHash = _.clone(self.listModel.get('order'));
            // oldorderHash = _.extend(oldorderHash, orderHash);

            //console.log('trigger listModel:change:order: ' + JSON.stringify(oldorderHash) );
            // self.listModel.set({'order': oldorderHash });

            //this.trigger("MyServerSideFilter:addOrder", orderHash, this);


            //self.fetch();
        },

        // Custom
        addSearch: function(searchHash){
            var self = this;
            var oldsearchHash = _.clone(self.listModel.get('search'));
            console.log('collection addSearch: current: ' + JSON.stringify(oldsearchHash) + ', adding: ' + JSON.stringify(searchHash));
            oldsearchHash = _.extend(oldsearchHash, searchHash);
            self.listModel.set({'search': oldsearchHash });

            //this.trigger("MyServerSideFilter:addSearch", searchHash, this);


            self.state.currentPage=1;  // if filtering, always set the page to 1
            // not needed because already called with backgrid-filter from header cells
            // self.fetch();
        },

        // Override
        getPage: function(page) { // override and hack in paging URL support
            console.log('getPage: ' + page);
            this.listModel.set({'page': page});

            var obj = Backbone.PageableCollection.prototype.getPage.call(this,page);
            return obj;
        },

        setPage: function(page){
            var page = parseInt(page);
            this.state.currentPage = page;
            this.fetch();
        },

        // TODO: verify this has been replaced by PageableCollection.setPageSize(pageSize, options) ...
        // setPage: function(page){
            // var state = _.clone(this.state);
            // state.pageSize = parseInt(page);
            // if(Math.ceil(state.totalRecords/state.pageSize) < state.currentPage ){
                // console.log('adjust currentPage');
                // state.currentPage=Math.ceil(state.totalRecords/state.pageSize);
            // }
            // this.collection.state = this.collection._checkState(state); // recalculate the state and do sanity checks.
//
            // this.collection.fetch({reset:true});
        // },

        // fetch: function(options) {
            // model = this;
            // if (typeof options !== 'undefined') { // TODO: review; find a better way to _always_ stop the ajax spinner
                // options.error = function(resp){
                    // console.log('error retrieving collection');
                // };
            // }
            // return Backbone.PageableCollection.prototype.fetch.call(this,options);
        // }
    });// end definition of Collection extension


    // NOTE: Backgrid instantiates the HeaderCell, so we don't have the listModel here explicitly.
    // Rather, we delegate up to the collection for actions, and the collection can set the listModel accordingly
    var MyServerSideFilter = Iccbl.MyServerSideFilter = Backgrid.Extension.ServerSideFilter.extend({
        // override so that we can keep a handle to the containing column name.
        // TODO: can handle this with events instead (so that the filter notifies the containing headercell?)
        columnName : null,  // TODO: use "name"

        initialize: function (options) {
            this.columnName = options.columnName;
            Backgrid.Extension.ServerSideFilter.prototype.initialize.apply(this, [options])

        },
    });

    // TODO: this should be named "contains" header cell, because it searches using "__contains"
    var MyHeaderCell = Iccbl.MyHeaderCell = Backgrid.HeaderCell.extend({

        _serverSideFilter : null,

        initialize: function (options) {

            Backgrid.HeaderCell.prototype.initialize.apply(this, [options]);

            this._serverSideFilter = new MyServerSideFilter({
              collection: this.collection, // TODO: Try to remove this: the collection should be passed as an option
              name: this.column.get("name")+"__contains", // the name of the URL query parameter for tastypie/django TODO: support other filters
              placeholder: "Search "+ this.column.get("label"), // HTML5 placeholder for the search box
              columnName: this.column.get("name"),
            });

            this.listenTo(this.collection, "MyServerSideFilter:search", this._search);

            this._serverSideFilter['events']['click .close'] = function(e) {
                if (e) e.preventDefault();

                this.remove(); // this is the filter
                this.clear();
                this.collection.clearSearch([this.name]);
                this.collection.fetch();
            };

            this._serverSideFilter['events']['submit'] = function(e) {
                var searchHash = {};
                searchHash[this.name] = this.$el.find("input[type=text]").val();
                console.log('server side filter add search: ' + JSON.stringify(searchHash));
                this.collection.addSearch(searchHash);
                this.search(e);
            };
        },

        remove: function(){
            console.log('headercell remove called');
            this._serverSideFilter.remove();
            this._serverSideFilter.unbind();
            this._serverSideFilter.collection = null;
            this.unbind();

            Backgrid.HeaderCell.prototype.remove.apply(this);
        },

        // function to listen for router generated custom search event MyServerSideFilter:search
        _search: function(searchHash, collection){
            var self = this;
            // console.log('Header cell respond to MyServerSideFilter:search trigger: ' + self._serverSideFilter.name + ', ' + JSON.stringify(searchHash));
            if (collection == this.collection) {
                var found = false;
                _.each(_(searchHash).pairs(), function(pair){
                    var key = pair[0];
                    var val = pair[1];
                    if( self._serverSideFilter.name == key){
                        console.log('--found search: ' + key + '=' + val + ', on: ' + self.column.get('name'));
                        found = true;
                        self.$el.append(self._serverSideFilter.render().el); // create the DOM element
                        self._serverSideFilter.$el.find("input[type=text]").val(val) // set the search term
                    }
                });

                if(!found){
                    if( !_.isEmpty(this.$el.find("input[type=text]").val())){
                        this._serverSideFilter.remove();
                        // this.collection.clearSearch([self._serverSideFilter.name])
                    }
                }
            }
        },

      /**
         Renders a header cell with a sorter and a label.
       */
      render: function () {
        Backgrid.HeaderCell.prototype.render.apply(this);
        var _handle = this;
        var filterIcon = $('<i class="icon-search"></i>');
        filterIcon.click(function(e) {
            _handle.$el.append(_handle._serverSideFilter.render().el);
        });
        this.$el.prepend(filterIcon);

        this.$el.prop('title', this.options['column']['attributes']["description"]);
        return this;
      }
    });

    return Iccbl;
});
