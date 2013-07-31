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
], function($, _, Backbone, BackbonePageableCollection, Backgrid, BackgridFilter, BackgridPaginator, lunr){

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

    // // Given an array property on an object, determine if it contains the value
    // var safeArrayPropertyContains = Iccbl.safeArrayPropertyContains = function(obj, property, value) {
        // if(!_.isArray(obj[property])){
            // var option = option[property];
            // return option.contains(value);
        // }
//
    // };

    // attach some objs to function scope //TODO: is this needed/correct?
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

    var MyServerSideFilter = Iccbl.MyServerSideFilter = Backgrid.Extension.ServerSideFilter.extend({
        // override so that we can keep a handle to the containing column name.
        // TODO: can handle this with events instead (so that the filter notifies the containing headercell?)
        columnName : null,

        initialize: function (options) {
            this.columnName = options.columnName;
            Backgrid.Extension.ServerSideFilter.prototype.initialize.apply(this, [options])

        },
    });

//    return {

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


    var MyCollection = Iccbl.MyCollection = Backbone.PageableCollection.extend({

        initialize: function(options){
            if (typeof(options) !== 'undefined'){
                //console.log('options: ' + JSON.stringify(options));
            }else{
                window.alert('no options defined');
            }
            // TODO: require these options
            this.url = options.url;
            //this.router = options.router; // TODO: set the AppModel.route property instead
            //this.type = options.type;

            // if(!_.isUndefined(options.searchBy)) this.searchBy = options.searchBy;

            Backbone.PageableCollection.prototype.initialize.apply(this, options); // call super constructor

        },

        url: function() {
            return this.url;
        },

        searchBy: null,
        model: MyModel,
        state: {
            pageSize: 25,
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
        parse: function(response) {
            // hack the response for tastypie:
            // note, this is because the pageable collection doesn't work with the backbone-tastypie.js fix
            //this.state.totalRecords = response.meta.total_count;
            var state = _.clone(this.state);
            state.totalRecords = response.meta.total_count;
            if(Math.ceil(state.totalRecords/state.pageSize) < state.currentPage ){
                console.log('adjust currentPage');
                state.currentPage=1;
            }
            this.state = this._checkState(state); // recalculate the state and do sanity checks.
            console.log('new state: ' + JSON.stringify(this.state));
            this.setRoutes();
            return response.objects;
        },
        /**
           @property {-1|0|1} [state.order=-1] The order to use for sorting. Specify
           -1 for ascending order or 1 for descending order. If 0, no client side
           sorting will be done and the order query parameter will not be sent to
         */

        setRoutes: function() {
            console.log('setRoutes: ' + this.searchBy + ', '
                + this.state.sortKey + ', ' + this.state.order + ', '+ this.state.pageSize + ', ', + this.state.currentPage);
            //var route = 'list/' + this.type;
            var route = '';
            if(this.searchBy !== null){
                if(route.length > 0) route += '/';
                route += 'search/'+this.searchBy ;
            }
            if(typeof this.state.sortKey !== 'undefined' && this.state.sortKey !== null){
                var sortKey = this.state.sortKey;
                if(this.state.order > 0){
                    sortKey = '-' + sortKey;
                }
                if(route.length > 0) route += '/';

                route += 'order_by/' + sortKey;
                //router.navigate('order_by/' + sortKey +'/page/'+ this.state.currentPage);
            }
            if(route.length > 0) route += '/';

            route += 'rpp/' + this.state.pageSize + '/page/' + this.state.currentPage;

            //this.router.navigate(route); // TODO: set AppModel.route instead
            this.trigger( "MyCollection:setRoute", route );
        },

        // Override
        setSorting: function(sortKey,order,options) { // override and hack in sorting URL support
            console.log('setSorting called: ' + sortKey+ ', order: ' + order + typeof(order) + ', options: ' + options );
            console.log('searchBy: ' + JSON.stringify(this.data));
            var dir = '-';
            var direction = 'descending';
            if(typeof order !== 'undefined' && order < 0){
                dir = '';
                direction = 'ascending';
            }
            var obj = Backbone.PageableCollection.prototype.setSorting.call(this,sortKey, order);

            return obj;
        },

        // Override
        getPage: function(page) { // override and hack in paging URL support
            var obj = Backbone.PageableCollection.prototype.getPage.call(this,page);
            //this.setRoutes();
            return obj;
        },

        fetch: function(options) {
            model = this;
            if (typeof options !== 'undefined') { // TODO: review; find a better way to _always_ stop the ajax spinner
                options.error = function(resp){
                    console.log('error retrieving collection');
                };
            }
            return Backbone.PageableCollection.prototype.fetch.call(this,options);
        }
    });// end definition of Collection extension


    var MyHeaderCell = Iccbl.MyHeaderCell = Backgrid.HeaderCell.extend({
        _serverSideFilter : null,
        initialize: function (options) {
            //this.constructor.__super__.initialize.apply(this, [options]);
            Backgrid.HeaderCell.prototype.initialize.apply(this, [options]);

            this._serverSideFilter = new MyServerSideFilter({
              collection: this.collection, // TODO: Try to remove this: the collection should be passed as an option

              name: this.column.get("name")+"__contains", // the name of the URL query parameter for tastypie/django TODO: support other filters
              placeholder: "Search "+ this.column.get("label"), // HTML5 placeholder for the search box
              columnName: this.column.get("name"),
            });

            this.listenTo(this.collection, "MyServerSideFilter:search", this._search);  // TODO: research use of "bindTo" instead
            //this.collection.bind("MyServerSideFilter:search", this._search, this);

            this._serverSideFilter['events']['click .close'] = function(e) {
                if (e) e.preventDefault();

                this.remove(); // this is the filter
                this.clear();
                this.collection.searchBy = null;
                this.collection.trigger('MyServerSideFilter:clearSearch');
            };

            // listen for search submission by the user and set the routes
            this._serverSideFilter['events']['submit'] = function(e) {
                // TODO: this should be handled by a backbone event
                this.collection.searchBy = this.columnName + '=' + this.$el.find("input[type=text]").val();
                this.collection.state.currentPage=1;  // if filtering, always set the page to 1
                //this.collection.setRoutes();
                this.search(e);
            };
            // _.bindAll(this, 'render');
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
        _search: function(searchColumn, searchTerm, collection){
            console.log('_search excuting from header cell');
            if (collection == this.collection) {
                if (searchColumn !== this.column.get("name")){
                    this._serverSideFilter.remove();
                }else{
                    console.log('on MyServerSideFilter:search: ' + searchColumn );
                    this.$el.append(this._serverSideFilter.render().el); // create the DOM element
                    this._serverSideFilter.$el.find("input[type=text]").val(searchTerm) // set the search term

                    // it works to just call fetch, I think, because the super constructor has
                    // the following, which is run when server-side data is retreived.
                      // if (Backbone.PageableCollection &&
                          // collection instanceof Backbone.PageableCollection &&
                          // collection.mode == "server") {
                        // collection.queryParams[this.name] = function () {
                          // return self.$el.find("input[type=text]").val();
                        // };
                      // }
                    // TODO: it might be better to explicitly call this.search, (but how?) (i.e. some event...)

                    this.collection.fetch({ reset: true });
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


    var ItemsPerPageSelector = Iccbl.ItemsPerPageSelector = Backbone.View.extend({

        template: function(){ return _.template(this._template); }, // TODO: this appears to be unused (see render)

        events: {
            "change #rowsperpage_selector": "setItemsPerPage"
        },

        initialize: function(options, collection){
            this.options = _(options.selections);
            this.collection = collection;
            this._template = options.template;
            this.render();
        },
        setItemsPerPage: function(e){
            e.preventDefault();
            var option = e.currentTarget.value;
            console.log('option: ' + option + ', clicked');
            // this.collection.state.pageSize = parseInt(option);
            var state = _.clone(this.collection.state);
            state.pageSize = parseInt(option);
            if(Math.ceil(state.totalRecords/state.pageSize) < state.currentPage ){
                console.log('adjust currentPage');
                state.currentPage=Math.ceil(state.totalRecords/state.pageSize);
            }
            this.collection.state = this.collection._checkState(state); // recalculate the state and do sanity checks.

            this.collection.fetch({reset:true});
            //this.collection.setRoutes();
        },
        render: function(){
            var selector_template = _.template(this._template, {
                options: this.options,
            });
            this.$el.html(selector_template);
            this.delegateEvents();
            return this;
        },
        updateSelection: function(){
            $('#rowsperpage_selector').val(String(this.collection.state.pageSize));
        },
    });
//    };
    // TODO: I don't think it's standard to return the main object this way?
    return Iccbl;
});
