define([
  'jquery',
  'underscore',
  'backbone',
  'backbone_pageable',
  'backgrid',
  'iccbl_backgrid',
  'views/detail_stickit',
  'text!templates/rows-per-page.html',
  'text!templates/list.html',
  'text!templates/modal_ok_cancel.html',
  'text!templates/generic-selector.html',
], function($, _, Backbone, BackbonePageableCollection, Backgrid,  Iccbl, DetailView, rowsPerPageTemplate, listTemplate, modalTemplate, genericSelectorTemplate ){

    // for compatibility with require.js, attach PageableCollection in the right place on the Backbone object
    // see https://github.com/wyuenho/backbone-pageable/issues/62
    Backbone.PageableCollection = BackbonePageableCollection;

    var ajaxStart = function(){
        $('#loading').fadeIn({duration:100});
    };
    var ajaxComplete = function(){
        $('#loading').fadeOut({duration:100});
    };
    $(document).bind("ajaxComplete", function(){
        ajaxComplete(); // TODO: bind this closer to the collection
    });

    var ListView = Backbone.View.extend({

        initialize : function(attributes, options) {
            console.log('initialize ListView');
            this.router = options.router;
            this._options = options;
        },

        // setOptions : function(options){
            // this.options = options;
        // },

        remove: function(){
            console.log('ListView remove called');
            Backbone.View.prototype.remove.apply(this);
        },

        reset_grid : function(){
            var self = this;
            // Call the schema from the server to get the field definitions
            $.ajax({
                type: "GET",
                url: this._options.url_schema,
                data: "",
                dataType: "json",
                success: function(result) {
                    self.buildGrid(result);
                }, // end success outer ajax call
                error: function(x, e) {
                    alert(x.readyState + " "+ x.status +" "+ e.msg);
                }
            });
        },

        buildGrid : function(schemaResult) {

            console.log('buildGrid...');

            var self = this; // todo: "this" is the parent ListView, all of this should be handled by events
            this.objects_to_destroy = _([]);
            // TODO: make sure req'd options are present (i.e. router, type)
            // type == table or query or api resource type
            var collection = this.collection = new Iccbl.MyCollection({ 'url': this._options.url }); //, router: this._options.router  });

            this.objects_to_destroy.push(collection);

            var columns = this.createBackgridColModel(schemaResult.fields, Iccbl.MyHeaderCell);//, col_options );
            //columns.unshift({ name: 'deletor', label: 'Delete', text:'X', description: 'delete record', cell: Iccbl.DeleteCell, sortable: false });
            var grid = this.grid = new Backgrid.Grid({
              columns: columns,
              collection: collection,
            });
            this.objects_to_destroy.push(grid);

            self.listenTo(collection, "MyCollection:link", function (model, column) {
                console.log('---- process link for '+ column);

                var fieldDef = schemaResult.fields[column];
                if( _.has(fieldDef,'backgrid_cell_options')) {
                    // NOTE: format for backgrid cell options is "/{attribute_key}/"
                    backgrid_cell_options = fieldDef['backgrid_cell_options'];
                    console.log('backgrid_cell_options: ' + backgrid_cell_options);

                    _route = backgrid_cell_options.replace(/{([^}]+)}/g, function (match) {
                        console.log('matched: ' + match + ', model: ' + model);
                        match = match.replace(/[{}]/g,'');
                        console.log('matched: ' + match + ', model: ' + model.get(match));
                        return typeof model.get(match) != "undefined" ? model.get(match) : match;
                    });
                    console.log('route: ' + _route);
                    this.router.navigate(_route, {trigger:true});
                }else{
                    console.log('no options defined for link cell');
                }
            });


            self.listenTo(collection, "MyCollection:edit", function (model) {
                console.log('---- create detail view for '+ this._options.ui_resource_id);
                // TODO: Get the title, and other items from the meta data
                // var detailView = new DetailView({ model: model},
                    // { title: "Details", fields:schemaResult.fields,
                      // app_model: this.model, resource_definition: schemaResult['resource_definition'],
                      // router: self.router  } ); // todo get title
//
                // $('#list-container').hide();

                //var _routeFrom = this.model.get('route');

                // Note: some links must use composite keys - because the composite key is the public key
                // (don't want to expose the private, possibly transient key)
                var id = '/' + model.get('id');
                if(_.has(schemaResult['resource_definition'], 'id_attribute')){
                    console.log('create id from ' + schemaResult['resource_definition']['id_attribute']);
                    id = _.reduce(schemaResult['resource_definition']['id_attribute'],
                            function(memo, item){ return memo += model.get(item) + '/';}, '');
                }else{
                    console.log('Warn: schema for this type has no resource_definition,id_attribute; type: ' + this._options.ui_resource_id);
                }
                console.log('id: ' + id);
                // TODO: Move route setting to the parent (contentview controller/view)
                var _route = 'detail/' + this._options.ui_resource_id + '/' + id;


                // // TODO: if we can make the "back" button do a navigate back, then all we have to do is set the route...
                // this.model.set({ route: _route } );

                console.log('-- set route: ' + _route);
                //this.router.navigate(_route, {trigger: true} );  // trigger false since don't want route actions firing
                this.model.set({    content_options: { schemaResult: schemaResult, model: model} ,
                                    current_view: 'detail',
                                    current_route_options: id } ); // signal to the app_model that the current view has changed // todo: separate out app_model from list_model
                //this.model.set({ current_view: 'detail_view'}, { silent: true } ); // silent: true since we don't want content view reacting to change event

                // NOTE: having self bind to the detailView like this:
                // self.listenTo(detailView, 'remove', function(){
                // causes the detailView to hang around in memory until self is closed
                // so either:
                // detailView.on('remove', function(){
                // self.listenToOnce(detailView, 'remove', function(){
                    // //self.collection.fetch({reset:true});
                    // $('#list-container').show();
                    // detailView.close();
                    // //this.model.set({ route: _routeFrom } );
                // });
//
                // $('#detail-container').append(detailView.render().$el);
            });

            self.listenTo(collection, "MyCollection:delete", function (model) {
                var modalDialog = new Backbone.View({
                    el: _.template(modalTemplate, { body: "Please confirm deletion of record: '" + model.get('toString') + "'", title: "Please confirm deletion" } ),
                    events: {
                        'click #modal-cancel':function(event) {
                            console.log('cancel button click event, '); // + JSON.stringify(fieldDefinitions));
                            event.preventDefault();
                            $('#modal').modal('hide'); // TODO: read-up on modal!  this is not ideal with the reference to template elements!
                        },
                        'click #modal-ok':function(event) {
                            console.log('ok button click event, '); // + JSON.stringify(fieldDefinitions));
                            event.preventDefault();
                            model.destroy();
                            $('#modal').modal('hide');
                        }
                    },
                });
                modalDialog.render();
                $('#modal').empty();
                $('#modal').html(modalDialog.$el);
                $('#modal').modal();
                console.log("removing model: " + JSON.stringify(model));
                console.log('----delete resource_uri: ' + model.get('resource_uri') );
                //model.destroy();
            });


            // TODO: is there a way to create this without extending Backbone.View (like footer below)?
            var extraSelector = Backbone.View.extend({
                initialize: function(options){
                    // console.log('extraselector init: '  + options);
                    this._options = options;

                },
                events: {
                    "change #generic_selector": "selectorChanged"
                },
                selectorChanged: function(e){
                    e.preventDefault();
                    var option = e.currentTarget.value;
                    var searchTerm, searchColumn, searchExpression;
                    searchTerm = option;
                    searchColumn = this._options.searchColumn
                    searchExpression = searchColumn + '=' + searchTerm;
                    collection.searchBy = searchExpression;
                    collection.trigger("MyServerSideFilter:search", { searchColumn: searchTerm }, collection);
                },
                updateSelection: function( searchItems ){
                    //console.log('-extraselector updateSelection: ' + JSON.stringify(searchItems) );
                    if( _(searchItems).has(this._options.searchColumn)){
                        $('#generic_selector').val(searchItems[this._options.searchColumn]);
                    }else{
                        console.log("extra selector set for column: " + this._options.searchColumn + ", not in searchItems: " + JSON.stringify(searchItems));
                    }
                },
                render: function(){
                    // console.log('===============render extra selector' + this)
                    this.delegateEvents();
                    this._options.options.unshift(' '); // create a blank entry

                    this.$el.html(_.template( genericSelectorTemplate,
                        { label: this._options.label,
                          'options': _(this._options.options ) }));  // TODO: this should come from the metahash schema
                    return this;
                },
            });

            var selector = new Iccbl.ItemsPerPageSelector(
                { 'selections': ['25','50','200','1000'], 'template': rowsPerPageTemplate },
                collection);
            this.objects_to_destroy.push(selector);

            var paginator = new Backgrid.Extension.Paginator({
              collection: collection
            });
            this.objects_to_destroy.push(paginator);

            var data = { message: '' };
            if (this._options.header_message){
                data.title = this._options.title;
                data.message = this._options.header_message; //'hello world!' };
            }
            var compiledTemplate = _.template( listTemplate, data );
            this.el.innerHTML = compiledTemplate;
            if( _.has(schemaResult, 'extraSelectorOptions')){
                console.log('searchTerms: ' + JSON.stringify(schemaResult.searchTerms));
                var extraSelectorInstance = new extraSelector( schemaResult.extraSelectorOptions ); // HMMm: need to understand the "delegate events" better here
                this.objects_to_destroy.push(extraSelectorInstance);
                $("#extra-selector-div").append(extraSelectorInstance.render().$el);
                extraSelectorInstance.listenTo(collection, 'MyServerSideFilter:search', extraSelectorInstance.updateSelection);
                extraSelectorInstance.listenTo(collection, 'MyServerSideFilter:clearSearch', extraSelectorInstance.updateSelection);

            }
            $("#paginator-div").append(paginator.render().$el);
            $("#rows-selector-div").append(selector.render().$el);
            $("#example-table").append(grid.render().$el);

            // encapsulate the footer in a view, help grab button click
            var footer = new Backbone.View({
                el: $("<form><button type='button' id='addRecord'>Add</button></form>"),
                events: {
                    'click button':function(event) {
                        console.log('button click event, '); // + JSON.stringify(fieldDefinitions));
                        event.preventDefault();
                        // TODO: set the defaults, also determine if should be set on create, from the Meta Hash
                        var defaults = {};
                        _.each(schemaResult.fields, function(value, key){
                            if (key == 'resource_uri') {
                                defaults[key] = self._options.url;
                            } else if (key == 'id'){ // nop // TODO: using the meta-hash, always exclude the primary key from create
                            } else {
                                 defaults[key] = '';
                            }
                        });
                        var NewModel = Backbone.Model.extend({urlRoot: self._options.url, defaults: defaults });
                        var detailView = new DetailView({ model: new NewModel}, { isEditMode: true, title: "Add new record", fields:schemaResult.fields}); // TODO: get the model "edit title" from the metainformation_hash

                        $('#list-container').hide();
                        // NOTE: having self bind to the detailView like this:
                        // self.listenTo(detailView, 'remove', function(){
                        // causes the detailView to hang around in memory until self is closed
                        // detailView.on('remove', function(){
                        self.listenToOnce(detailView, 'remove', function(){
                            console.log('... remove detailView event');
                            self.collection.fetch({reset:true});
                            $('#list-container').show();
                            detailView.close();
                        });
                        $('#detail-container').append(detailView.render().$el);

                    },
                },

            });
            $("#table-footer-div").append(footer.$el);
            this.objects_to_destroy.push(footer);


            // Note on event subscriptions: prefer listenTo over "on" (alias for _.bind/model.bind) as this
            // will allow the object to unbind all observers at once.
            //collection.on('request', ajaxStart); // NOTE: can use bind or on
            //collection.bind('sync', ajaxComplete, this);

            this.listenTo(collection, 'request', ajaxStart);
            this.listenTo(collection, 'MyCollection:setRoute', this.setRoute);
            this.listenTo(collection, 'sync', selector.updateSelection );

            // TODO: work out the specifics of communication complete event.  the following are superceded by the global handler for "ajaxComplete"
            this.listenTo(collection, 'error', ajaxComplete);
            this.listenTo(collection, 'complete', ajaxComplete);

            // var allEvent = function(event){
                // console.log("event fired: " + event);
            // };
            // self.on({
              // "all": allEvent,  // TODO: this is debug code
            // });

            if (this._options.page){
                collection.state.currentPage = this._options.page;
            }
            if (this._options.pageSize){
                collection.state.pageSize = this._options.pageSize;
            }
            if(typeof this._options.orderBy !== 'undefined' && this._options.orderBy !== null ){
                var direction = 'ascending';
                var order = -1; // according to the docs, -1 == ascending
                var sortKey = options.orderBy;
                if(sortKey.charAt(0) === '-'){
                    order = 1; // according to the docs, 1 == descending
                    direction = 'descending';
                    sortKey = sortKey.substring(1);
                }
                collection.state.sortKey = sortKey;
                collection.state.order = order;
                // Notify header cells
                collection.trigger("backgrid:sort", sortKey, direction, null, collection);
            }
            if(typeof this._options.searchBy !== 'undefined' && this._options.searchBy !== null){

                var searchExpressions = {};
                // TODO: query terms that tastypie will understand.  these are to be set on the MyHeaderCell
                // QUERY_TERMS = set([
                    // 'exact', 'iexact', 'contains', 'icontains', 'gt', 'gte', 'lt', 'lte', 'in',
                    // 'startswith', 'istartswith', 'endswith', 'iendswith', 'range', 'year',
                    // 'month', 'day', 'week_day', 'isnull', 'search', 'regex', 'iregex',
                // ])

                _(this._options.searchBy.split(',')).each(function(searchItem){
                    var searchExpression = searchItem.split('=');
                    if(searchExpression.length != 2 ){
                        console.log('Warning: invalid search item: ' + searchItem + ', in: ' + this_.options.searchBy);
                    }else{
                        searchExpressions[searchExpression[0]] = searchExpression[1];
                    }
                    // var p = /([^=]+)=([^=]+)/
                    // var match = p.exec(searchItem);
                    // if (match){
                        // // data[match[1] + '__contains'] = match[2];
                        // // console.log('parsed search: ' + JSON.stringify(data));
                        // var searchColumn = match[1];
                        // var searchTerm = match[2];
                });

                if(!_.isEmpty(searchExpressions)){
                    collection.searchBy = this._options.searchBy;
                    // console.log('trigger search:' + this._options.searchBy );
                    collection.trigger("MyServerSideFilter:search", searchExpressions, collection);
                    // console.log('done: trigger search');
                }else{
                    console.log('Warn: no search expressions found in: '  + this._options.searchBy );
                }
            }
            console.log('collection fetch trigger');
            collection.fetch({ reset: true } );

            console.log('list view initialized');
        },

        /**
         *
         * @param {Object} fields_from_rest - hash of fields for the current dataset:
         *      field properties { visibility: [array of strings], title: a label for the field, order: display order of the field }
         * @param {Object} optionalHeaderCell - a Backgrid.HeaderCell to use for each column
         * @param {Object} options - a hash of { fieldKey: [custom cell: extend Backgrid.Cell] } to map custom cell implementations to fields
         */
        createBackgridColModel: function(restFields, optionalHeaderCell) {
            console.log('--createBackgridColModel'); //: restFields: ' + JSON.stringify(restFields));
            var colModel = [];
            var i = 0;
            var _total_count = 0;
            _.each(_.pairs(restFields), function(pair){
                var key = pair[0];
                var prop = pair[1];

                var visible = _.has(pair[1], 'visibility') && _.contains(pair[1]['visibility'], 'list');
                if(visible){

                    var backgridCellType = 'string';
                    if( !_.isEmpty(prop['backgrid_cell_type'])){
                        backgridCellType = prop['backgrid_cell_type'];
                        try{
//                            console.log('look for ' + key + ', ' + prop['backgrid_cell_type']);
                            var klass = Iccbl.stringToFunction(prop['backgrid_cell_type']);
//                            console.log('got  ' + klass);
                            if(!_.isUndefined(klass)){
//                                console.log('----- cell found: ' + klass);
                                backgridCellType = klass;
                            }
                        }catch(ex){
                            var msg = '----Warn: field: ' + key + ', no Iccbl class found for type: ' + prop['backgrid_cell_type'];
                            console.log(msg + ': ' + JSON.stringify(ex));
                        }
                    }
                    colModel[i] = {
                        'name':key,
                        'label':prop['title'],
                        'description':prop['description'],
                        cell: backgridCellType,
                        order: prop['ordinal'],
                        editable: false,
                    };
                    if (optionalHeaderCell){
                        colModel[i]['headerCell'] = optionalHeaderCell;
                    }
                    i++;
                }else{
                    //console.log('field not visible in list view: ' + key)
                }
            });


            //console.log('colModel: ' + JSON.stringify(colModel));
            colModel.sort(function(a,b){
                if(_.isNumber(a['order']) && _.isNumber( b['order'])){
                    return a['order']-b['order'];
                }else if(_.isNumber( a['order'])){
                    return -1;
                }else if(_.isNumber(b['order'])){
                    return 1;
                }else{
                    return 0;
                }
            });
            return colModel;
        },

        setRoute: function(route, options){
            // var _route = 'list/' + this._options.ui_resource_id + '/' + route;
            // console.log('setRoute triggered: ' + _route + ', ' + JSON.stringify(options) );
            // this.router.navigate(_route, options);

            this.model.set({ current_route_update: route });
        },

        onClose: function(){
            console.log('Extra onclose method called');
            if(_.isObject(this.objects_to_destroy)){
                this.objects_to_destroy.each(function(view_obj){
                    view_obj.remove();
                    view_obj.unbind();
                    view_obj.stopListening();
                });
            }
        },

        render: function(){
            console.log('render listView');
            this.reset_grid();
            return this;
        }
    });

  return ListView;
});