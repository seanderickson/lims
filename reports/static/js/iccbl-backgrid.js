/**
 * @summary     ICCBL-lims Utility Functions
 * @description Utility Functions for the iccbl-lims
 * @version     0.1
 * @file        iccbl-lims.js
 * @author      Sean Erickson
 * @contact     sean_erickson “AT”hms.harvard.edu
 *
 * @copyright Copyright 2013 Harvard University, all rights reserved.
 *
 * This source file is free software, under either the GPL v2 license
 *
 * This source file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license files for details.
 *
 * For details please refer to: https://github.com/hmsiccbl/lims
 *
 * Depends on iccbl-lims.js
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

    MyModel = Backbone.Model.extend({
      initialize: function () {
        Backbone.Model.prototype.initialize.apply(this, arguments);
        var self = this;
        this.on("change", function (model, options) {
          if (options && options.save === false) return;
          // TODO: not sure what the best way to "fix" the URL for the model is yet, but here is a way -sde4
          self._url = model.id;
          self.url = function() {
              return this._url;
          }
          model.save();
        });
      },
    });

    MyServerSideFilter = Backgrid.Extension.ServerSideFilter.extend({
        // override so that we can keep a handle to the containing column name.
        // TODO: can handle this with events instead (so that the filter notifies the containing headercell?)
        columnName : null,

        initialize: function (options) {
            this.columnName = options.columnName;
            Backgrid.Extension.ServerSideFilter.prototype.initialize.apply(this, [options])

        },

        // remove: function() {
          // this.model.unbind('change', this.render);
          // return Backbone.View.prototype.remove.apply(this, arguments);
        // },
    });

    return {

        createBackgridColModel: function(fields_from_rest, optionalHeaderCell, collection) {
            var colModel = [];
            var i = 0;
            var _total_count = 0;
            for (var field in fields_from_rest){
                if (fields_from_rest.hasOwnProperty(field)) { // filter
                    var prop = fields_from_rest[field];
                    colModel[i] = { 'name':field, 'label':prop['name'], cell: 'string', order: prop['order']};
                    if (optionalHeaderCell){
                        colModel[i]['headerCell'] = optionalHeaderCell;
                    }
                    i++;
                }
            }

            colModel.sort(function(a,b){
                //console.log('sort: ' + a['order'] +',' + b['order']);
                if(typeof a['order'] !== 'undefined' && typeof b['order'] !== 'undefined'){
                    return a['order']-b['order'];
                }else if(typeof a['order'] !== 'undefined'){
                    return -1;
                }else if(typeof b['order'] !== 'undefined'){
                    return 1;
                }else{
                    return 0;
                }
            });
            //console.log('colModel: ' + JSON.stringify(colModel));
            //var _colWidth = 1/i * _width;
            //console.log('colWidth:' + _colWidth);
            return colModel;
        },


        MyCollection: Backbone.PageableCollection.extend({

            initialize: function(options){
                if (typeof(options) !== 'undefined'){
                    //console.log('options: ' + JSON.stringify(options));
                }else{
                    window.alert('no options defined');
                }
                // TODO: require these options
                this.url = options.url;
                //this.router = options.router; // TODO: set the AppModel.route property instead
                this.type = options.type;

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
                },

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
//                        window.App.ajaxComplete();
                        console.log('error retrieving collection');
                        // TODO: bind this to the ajaxComplete handler
                    };
                }
                return Backbone.PageableCollection.prototype.fetch.call(this,options);
            }
        }),// end definition of Collection extension


        MyHeaderCell: Backgrid.HeaderCell.extend({
            _serverSideFilter : null,
            initialize: function (options) {
                // console.log('MyHeaderCell initialize, options: ' + JSON.stringify(options) ); // TODO: looking for the collection instance here
                this.constructor.__super__.initialize.apply(this, [options])


                this._serverSideFilter = new MyServerSideFilter({
                  collection: this.collection, // TODO: Try to remove this: the collection should be passed as an option

                  name: this.column.get("name")+"__contains", // the name of the URL query parameter for tastypie/django TODO: support other filters
                  placeholder: "Search "+ this.column.get("label"), // HTML5 placeholder for the search box
                  columnName: this.column.get("name"),
                });

                this.listenTo(this.collection, "MyServerSideFilter:search", this._search);  // TODO: research use of "bindTo" instead

                this._serverSideFilter['events']['click .close'] = function(e) {
                    if (e) e.preventDefault();

                    this.remove(); // this is the filter
                    this.clear();
                    this.collection.searchBy = null;
                };

                // listen for search submission by the user and set the routes
                this._serverSideFilter['events']['submit'] = function(e) {
                    // TODO: this should be handled by a backbone event
                    this.collection.searchBy = this.columnName + '=' + this.$el.find("input[type=text]").val();
                    this.collection.state.currentPage=1;  // if filtering, always set the page to 1
                    //this.collection.setRoutes();
                    this.search(e);
                };
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
            return this;
          }
        }),


        ItemsPerPageSelector: Backbone.View.extend({
        //        el: $("#selector-div"), // view must have the element to get events for
            // template: _.template($("#rows-per-page-template").html()),
            template: function(){ return _.template(this._template); },

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
            updateSelection: function(e){
                $('#rowsperpage_selector').val(String(this.collection.state.pageSize));
            },
        }),


    };
});
