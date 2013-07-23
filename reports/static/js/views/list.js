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

        initialize : function() {
            console.log('initialize ListView');
        },

        setOptions : function(options){
            this.options = options;
        },

        remove: function(){
            console.log('ListView remove called');
            Backbone.View.prototype.remove.apply(this);
        },

        reset_grid : function(options){
            var self = this;
            // Call the schema from the server to get the field definitions
            $.ajax({
                type: "GET",
                url: options.url_schema,
                data: "",
                dataType: "json",
                success: function(result) {
                    self.buildGrid(result, options.url, options);
                }, // end success outer ajax call
                error: function(x, e) {
                    alert(x.readyState + " "+ x.status +" "+ e.msg);
                }
            });
        },

        wrapAsUnderscore: function(obj) {
            var self = this;
            var newObj = {};
            _.each(_.pairs(obj), function(pair){
                if(_.isArray(pair[1])){
                    newObj[pair[0]] = _(pair[1]);
                }else if(_.isObject(pair[1])){
                    newObj[pair[0]] = self.wrapAsUnderscore(pair[1]);
                }else{
                    newObj[pair[0]] = pair[1];
                }
            });
            return newObj;

        },

        buildGrid : function(schemaResult, _url, options) {

            var fieldDefinitions = schemaResult.fields;

            var field_defs = this.wrapAsUnderscore(fieldDefinitions);
            // console.log('field_defs: ' + JSON.stringify(field_defs));

            console.log('buildGrid...');

            var self = this; // todo: "this" is the parent ListView, all of this should be handled by events
            this.objects_to_destroy = _([]);
            // TODO: make sure req'd options are present (i.e. router, type)
            // type == table or query or api resource type
            var collection = this.collection = new Iccbl.MyCollection({ 'url': _url, router: options.router, type: options.type  });

            this.objects_to_destroy.push(collection);

            // var ClickableRow = Backgrid.Row.extend({
              // events: {
                // "click": "onClick"
              // },
              // onClick: function () {
                // self.trigger("MyCollection:edit", this.model);
              // }
            // });


            // TODO: p-o-c of how to do a Link-to-Detail/Edit: should be controlled by the field type
            //var col_options = { customCells: { 'title': Iccbl.EditCell, 'screensaver_user_id': Iccbl.EditCell } };
            var columns = this.createBackgridColModel(fieldDefinitions, Iccbl.MyHeaderCell);//, col_options );
            columns.unshift({ name: 'deletor', label: 'Delete', text:'X', description: 'delete record', cell: Iccbl.DeleteCell, sortable: false });
            var grid = this.grid = new Backgrid.Grid({
              columns: columns,
              collection: collection,
              // row: ClickableRow,
            });
            this.objects_to_destroy.push(grid);

            self.listenTo(collection, "MyCollection:edit", function (model) {
                console.log('edit:', model);
                // TODO: Get the title, and other items from the meta data
                var detailView = new DetailView({ model: model}, {title: "Details", fields:field_defs} ); // TODO: get the model "edit title" from the metainformation_hash

                $('#list-container').hide();
                // NOTE: having self bind to the detailView like this:
                // self.listenTo(detailView, 'remove', function(){
                // causes the detailView to hang around in memory until self is closed
                // so either:
                // detailView.on('remove', function(){
                self.listenToOnce(detailView, 'remove', function(){
                    self.collection.fetch({reset:true});
                    $('#list-container').show();
                    detailView.close();
                });

                $('#detail-container').append(detailView.render().$el);
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
                    console.log('init: '  + options);
                    this.label = options.label ;
                    this._options = options.options;
                    this._options.unshift(' ');

                },
                events: {
                    "change #generic_selector": "selectorChanged"
                },
                selectorChanged: function(e){
                    e.preventDefault();
                    var option = e.currentTarget.value;
                    var searchTerm, searchColumn, searchExpression;
                    searchTerm = option;
                    searchColumn = 'scope'
                    searchExpression = searchColumn + '=' + searchTerm;
                    collection.searchBy = searchExpression;
                    collection.trigger("MyServerSideFilter:search", searchColumn, searchTerm, collection);
                },
                updateSelection: function(searchColumn, searchTerm){
                    console.log('updateSelection: ' + searchColumn + ', ' + searchTerm);
                    $('#generic_selector').val(searchTerm);
                },
                render: function(){
                    console.log('===============render extra selector' + this)
                    this.delegateEvents();
                    this.$el.html(_.template( genericSelectorTemplate,
                        { label: this.label,
                          'options': _(this._options) }));  // TODO: this should come from the metahash schema
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
            if (this.options.header_message){
                data.title = this.options.title;
                data.message = this.options.header_message; //'hello world!' };
            }
            var compiledTemplate = _.template( listTemplate, data );
            this.el.innerHTML = compiledTemplate;
            if( _.has(schemaResult, 'searchTerms')){
                console.log('searchTerms: ' + JSON.stringify(schemaResult.searchTerms));
                var extraSelectorInstance = new extraSelector({ label:'Scope: ', options: schemaResult.searchTerms}); // HMMm: need to understand the "delegate events" better here
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
                        _.each(fieldDefinitions, function(value, key, list){
                            if (key == 'resource_uri') {
                                defaults[key] = self.options.url;
                            } else if (key == 'id'){ // nop // TODO: using the meta-hash, always exclude the primary key from create
                            } else {
                                 defaults[key] = '';
                            }
                        });
                        var NewModel = Backbone.Model.extend({urlRoot: self.options.url, defaults: defaults });
                        var detailView = new DetailView({ model: new NewModel}, { isEditMode: true, title: "Add new record", fields:field_defs}); // TODO: get the model "edit title" from the metainformation_hash

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

            if (options.page){
                collection.state.currentPage = options.page;
            }
            if (options.pageSize){
                collection.state.pageSize = options.pageSize;
            }
            if(typeof options.orderBy !== 'undefined' && options.orderBy !== null ){
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
            if(typeof options.searchBy !== 'undefined' && options.searchBy !== null){
                // TODO: only can search one term at a time
                var p = /([^=]+)=([^=]+)/
                var match = p.exec(options.searchBy);
                if (match){
                    // data[match[1] + '__contains'] = match[2];
                    // console.log('parsed search: ' + JSON.stringify(data));
                    var searchColumn = match[1];
                    var searchTerm = match[2];
                    console.log('parsed search: ' + searchColumn + ', ' + searchTerm );

                    collection.searchBy = options.searchBy;
                    console.log('trigger search:' + searchColumn + ', ' + searchTerm );
                    collection.trigger("MyServerSideFilter:search", searchColumn, searchTerm, collection);
                    console.log('done: trigger search');

                 }else{
                     console.log('unrecognized searchBy: ' + options.searchBy);
                 }
            }

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
            // console.log('createBackgridColModel: restFields: ' + JSON.stringify(restFields));
            var colModel = [];
            var i = 0;
            var _total_count = 0;
            _.each(_.pairs(restFields), function(pair){
                var key = pair[0];
                var prop = pair[1];

                var visible = _.has(pair[1], 'visibility') && _.contains(pair[1]['visibility'], 'list');
                //visible = visible || options.customCells.hasOwnProperty(key);
//                console.log('field: ' + key + ', ' + options.customCells.hasOwnProperty(key) + ', ' + visible);
                if(visible){
                    // console.log('process list field: ' + pair[0] + ', ' + JSON.stringify(pair[1]));
                    //var cell = 'string';
                    // var cell = prop['ui_type']
//
                    // if(cell == 'choice' || cell == 'multiselect'){ // TODO: should have two types, "gridCellType", and "editFormType"?
                        // cell = 'string';
                    // }
                    var backgridCellType = 'string';
                    if( !_.isEmpty(prop['backgrid_cell_type'])){
                        backgridCellType = prop['backgrid_cell_type'];
                        try{
                            console.log('look for ' + key + ', ' + prop['backgrid_cell_type']);
                            var klass = Iccbl.stringToFunction(prop['backgrid_cell_type']);
                            console.log('got  ' + klass);
                            if(!_.isUndefined(klass)){
                                console.log('----- cell found: ' + klass);
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
                    console.log('field not visible in list view: ' + key)
                }
            });


            console.log('colModel: ' + JSON.stringify(colModel));
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








        //////////TODO: DELETE!




        //////////TODO: DELETE!
        createBackgridColModel_orig: function(restFields, optionalHeaderCell, options) {
            // console.log('createBackgridColModel: restFields: ' + JSON.stringify(restFields));
            var colModel = [];
            var i = 0;
            var _total_count = 0;
            _.each(_.pairs(restFields), function(pair){
                var key = pair[0];
                var prop = pair[1];

                var visible = _.has(pair[1], 'visibility') && _.contains(pair[1]['visibility'], 'list');
                visible = visible || options.customCells.hasOwnProperty(key);
                console.log('field: ' + key + ', ' + options.customCells.hasOwnProperty(key) + ', ' + visible);
                if(visible){
                    // console.log('process list field: ' + pair[0] + ', ' + JSON.stringify(pair[1]));
                    //var cell = 'string';
                    var cell = prop['ui_type']

                    if(cell == 'choice' || cell == 'multiselect'){ // TODO: should have two types, "gridCellType", and "editFormType"?
                        cell = 'string';
                    }

                    if(!_.isUndefined(options['customCells'][key])){ // TODO this should go in the metahash
                        cell = options['customCells'][key];
                    }
                    colModel[i] = {
                        'name':key,
                        'label':prop['title'],
                        'description':prop['description'],
                        cell: cell,
                        order: prop['ordinal'],
                        editable: false,
                    };
                    if (optionalHeaderCell){
                        colModel[i]['headerCell'] = optionalHeaderCell;
                    }
                    i++;
                }else{
                    console.log('field not visible in list view: ' + key)
                }
            });


            // console.log('colModel: ' + JSON.stringify(colModel));
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


        setRoute: function(route){
            var _route = 'list/' + this.options.type + '/' + route;
            console.log('setRoute triggered: ' + _route);
            this.model.set({ route: _route } );
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

        render: function(){ // TODO: should the build grid be called from here?
            console.log('render listView');
            this.reset_grid(this.options);
            return this;
        }
    });

  return ListView;
});