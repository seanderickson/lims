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
            console.log('initialize ListView: ');
            var self = this;
            Iccbl.assert( !_.isUndefined(options.ui_resource_id), 'listView options: ui_resource_id is required');
            Iccbl.assert( !_.isUndefined(options.router), 'listView options: router is required');
            Iccbl.assert( !_.isUndefined(options.url), 'listView options: url is required');
            Iccbl.assert( !_.isUndefined(options.schemaResult), 'listView options: schemaResult is required');
            // optional: Iccbl.assert( !_.isUndefined(options.url_schema), 'listView options: url_schema is required');
            Iccbl.assert( !_.isUndefined(options.header_message), 'listView options: header_message is required');
            Iccbl.assert( !_.isUndefined(options.title), 'listView options: title is required');

            this.router = options.router;
            this._options = options;

            var ListModel = Backbone.Model.extend({ defaults: { rpp: 25, page: 1, order: {}, search: {}} });
            var listModel = this.listModel = new ListModel(_.clone(this.model.get('current_options')));
            this.listenTo(this.listModel, 'change:search', function(){
                console.log('===--- listModel change:search');

                var current_options = _.clone(self.model.get('current_options'));
                current_options.search = self.listModel.get('search');
                self.model.set({current_options: current_options });
            });
            this.listenTo(this.listModel, 'change:order', function(){
                console.log('===--- listModel change:order');

                var orderHash = self.listModel.get('order');
                var current_options = _.clone(self.model.get('current_options'));
                current_options.order = orderHash;
                self.model.set({current_options: current_options });
            });

            this.listenTo(this.listModel, 'change:rpp', function(){
                console.log('===--- listModel change:rpp');
                var pageSize = parseInt(self.listModel.get('rpp'));

                // self.collection.setPageSize(pageSize);
                // // Triggers a fetch

                var current_options = _.clone(self.model.get('current_options'));
                current_options.rpp = pageSize;
                self.model.set({current_options: current_options });
            });

            this.listenTo(this.listModel, 'change:page', function(){
                console.log('===--- listModel change:page');
                var page = parseInt(self.listModel.get('page'));

                // Not needed since page changes come from the paginator->pageablecollection
                // self.collection.state.currentPage = page;
                // self.collection.fetch();

                var current_options = _.clone(self.model.get('current_options'));
                current_options.page = page;
                self.model.set({current_options: current_options });
            });

            var data = { message: '' };
            if (this._options.header_message){
                data.title = this._options.title;
                data.message = this._options.header_message; //'hello world!' };
            }
            var compiledTemplate = this.compiledTemplate = _.template( listTemplate, data );
            this.objects_to_destroy = _([]);

            var collection  = this.collection = new Iccbl.MyCollection({
                    'url': this._options.url,
                    currentPage: parseInt(self.listModel.get('page')),
                    pageSize: parseInt(self.listModel.get('rpp')),
                    listModel: listModel
                });
            this.objects_to_destroy.push(collection);

            this.buildGrid(this._options.schemaResult);
        },

        buildGrid : function(schemaResult) {

            console.log('buildGrid...');

            var self = this; // todo: "this" is the parent ListView, all of this should be handled by events

            self.listenTo(self.collection, "MyCollection:link", function (model, column) {
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

            self.listenTo(self.collection, "MyCollection:edit", function (model) {
                console.log('---- create detail view for '+ this._options.ui_resource_id);
                // Note: some links must use composite keys - because the composite key is the public key
                // (don't want to expose the private, possibly transient key)
                var id = '/' + model.get('id');
                if(_.has(schemaResult['resource_definition'], 'id_attribute')){
                    console.log('create id from ' + schemaResult['resource_definition']['id_attribute']);
                    id = _.reduce(schemaResult['resource_definition']['id_attribute'],
                            function(memo, item){
                                if(!_.isEmpty(memo)) memo += '/';
                                return memo += model.get(item);
                            }, '');
                }else{
                    console.log('Warn: schema for this type has no resource_definition,id_attribute; type: ' + this._options.ui_resource_id);
                }
                console.log('id: ' + id);
                this.model.set({    current_scratch: { schemaResult: schemaResult, model: model} ,
                                    current_view: 'detail',
                                    current_options: { key: id },
                                    routing_options: {trigger: false, replace: false}
                               }); // signal to the app_model that the current view has changed // todo: separate out app_model from list_model
            });

            self.listenTo(self.collection, "MyCollection:delete", function (model) {
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

            // Rows-per-page selector
            var rppModel = new Backbone.Model({ selection: this.listModel.get('rpp') });
            var rppSelectorInstance = self.rppSelectorInstance = new Iccbl.GenericSelector(
                { model: rppModel }, {label: 'Rows per page:', options: ['25','50','200','1000']} );
            this.objects_to_destroy.push(rppSelectorInstance);
            self.$("#extra-selector-div").append(rppSelectorInstance.render().$el);
            this.listenTo(this.listModel, 'change: rpp', function(){
                rppModel.set({ selection: self.listModel.get('rpp') });
            });
            this.listenTo(rppModel, 'change', function() {
                console.log('===--- rppModel change');
                self.listModel.set('rpp', rppModel.get('selection'));
            });

            // Paginator
            var paginator = self.paginator = new Backgrid.Extension.Paginator({
              collection: self.collection
            });
            this.objects_to_destroy.push(paginator);

            // Extraselector
            if( _.has(schemaResult, 'extraSelectorOptions')){
                console.log('searchTerms: ' + JSON.stringify(schemaResult.searchTerms));

                var extraSelectorModel = new Backbone.Model({ selection: '' });
                var extraSelectorKey = schemaResult.extraSelectorOptions.searchColumn;
                var extraSelectorInstance = self.extraSelectorInstance = new Iccbl.GenericSelector({ model: extraSelectorModel} , schemaResult.extraSelectorOptions ); // HMMm: need to understand the "delegate events" better here
                this.objects_to_destroy.push(extraSelectorInstance);
                self.$("#extra-selector-div").append(extraSelectorInstance.render().$el);

                this.listenTo(this.listModel, 'change: search', function(){
                    var searchHash = self.listModel.get('search');
                    _.each(_.keys(searchHash), function(key){
                        if( key == extraSelectorKey){
                            extraSelectorModel.set({ selection: searchHash[key] });
                        }
                    });
                });
                this.listenTo(extraSelectorModel, 'change', function() {
                    console.log('===--- extraSelectorModel change');
                    var searchHash = _.clone(self.listModel.get('search'));
                    var value = extraSelectorModel.get('selection');
                    searchHash[extraSelectorKey] = value;
                    self.listModel.set('search', searchHash);
                });

                // extraSelectorInstance.listenTo(self.collection, 'MyServerSideFilter:search', extraSelectorInstance.updateSelection);
                // extraSelectorInstance.listenTo(self.collection, 'MyServerSideFilter:clearSearch', extraSelectorInstance.updateSelection);

            }
            // self.$("#paginator-div").append(paginator.render().$el);
            // self.$("#rows-selector-div").append(selector.render().$el);
            var columns = this.createBackgridColModel(this._options.schemaResult.fields, Iccbl.MyHeaderCell);//, col_options );
            //columns.unshift({ name: 'deletor', label: 'Delete', text:'X', description: 'delete record', cell: Iccbl.DeleteCell, sortable: false });
            var grid = this.grid = new Backgrid.Grid({
              columns: columns,
              collection: self.collection,
            });
            self.$("#example-table").append(grid.render().$el);
            this.objects_to_destroy.push(grid);

            // encapsulate the footer in a view, help grab button click
            var footer = self.footer = new Backbone.View({
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
            self.$("#table-footer-div").append(footer.$el);
            this.objects_to_destroy.push(footer);


            // Note on event subscriptions: prefer listenTo over "on" (alias for _.bind/model.bind) as this
            // will allow the object to unbind all observers at once.
            //collection.on('request', ajaxStart); // NOTE: can use bind or on
            //collection.bind('sync', ajaxComplete, this);

            this.listenTo(self.collection, 'request', ajaxStart);
            // this.listenTo(collection, 'MyCollection:setRoute', this.setRoute);
            this.listenTo(self.collection, 'MyCollection:changeOptions', this.change_options);
//            this.listenTo(self.collection, 'sync', this.selector.updateSelection );

            // TODO: work out the specifics of communication complete event.  the following are superceded by the global handler for "ajaxComplete"
            this.listenTo(self.collection, 'error', ajaxComplete);
            this.listenTo(self.collection, 'complete', ajaxComplete);

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
                                console.log('----- class found: ' + klass);
                                backgridCellType = klass;
                            }
                        }catch(ex){
                            var msg = '----for: field: ' + key + ', no Iccbl class found for type: ' + prop['backgrid_cell_type'] + ', this may be a Backgrid cell type';
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

        change_options: function(new_options, _routing_options ){
            console.log('changeOptions triggered: ' + JSON.stringify(new_options) + ' , ' + JSON.stringify(_routing_options) );

            var updated_current_options = _.extend({}, this.model.get('current_options'), new_options);
            this.model.set({ routing_options: _routing_options,
                             current_options: updated_current_options });
        },

        remove: function(){
            console.log('ListView remove called');
            Backbone.View.prototype.remove.apply(this);
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
            var self = this;
            this.$el.html(this.compiledTemplate);
            self.$("#example-table").append(this.grid.render().$el);
            self.$("#paginator-div").append(self.paginator.render().$el);
            self.$("#rows-selector-div").append(self.rppSelectorInstance.render().$el);

            if(!_.isUndefined(self.extraSelectorInstance)){
                self.$("#extra-selector-div").append(self.extraSelectorInstance.render().$el);
            }else{
                console.log('wherez theextra selector');
            }


            // // TODO: can set the fetchOptions directly as well
            // if (!_.isEmpty(self.listModel.get('page'))){
                // self.collection.state.currentPage = this._options.page;
            // }
            // if (!_.isEmpty(self.listModel.get('rpp'))){
                // // self.collection.state.pageSize = this._options.rpp;
                // self.collection.setPageSize(parseInt(self.listModel.get('rpp')));
            // }
            // if(_.isString(this._options.order) ){
                // var direction = 'ascending';
                // var order = -1; // according to the docs, -1 == ascending
                // var sortKey = this._options.order;
                // if(sortKey.charAt(0) === '-'){
                    // order = 1; // according to the docs, 1 == descending
                    // direction = 'descending';
                    // sortKey = sortKey.substring(1);
                // }
                // self.collection.state.sortKey = sortKey;
                // self.collection.state.order = order;
                // // Notify header cells
                // self.collection.trigger("backgrid:sort", sortKey, direction, null, self.collection);
            // }
            var orderHash = self.listModel.get('order');
            if(!_.isEmpty(orderHash)){

                self.collection.setOrder(orderHash);

            }

            var searchHash = self.listModel.get('search');
            if(!_.isEmpty(searchHash)){

                self.collection.setSearch(searchHash);

            }

            if(_.isString(this._options.search)){
                // var searchExpressions = {};
                // // TODO: query terms that tastypie will understand.  these are to be set on the MyHeaderCell
                // // QUERY_TERMS = set([
                    // // 'exact', 'iexact', 'contains', 'icontains', 'gt', 'gte', 'lt', 'lte', 'in',
                    // // 'startswith', 'istartswith', 'endswith', 'iendswith', 'range', 'year',
                    // // 'month', 'day', 'week_day', 'isnull', 'search', 'regex', 'iregex',
                // // ])
//
                // _(this._options.search.split(',')).each(function(searchItem){
                    // var searchExpression = searchItem.split('=');
                    // if(searchExpression.length != 2 ){
                        // console.log('Warning: invalid search item: ' + searchItem + ', in: ' + this._options.search);
                    // }else{
                        // searchExpressions[searchExpression[0]] = searchExpression[1];
                    // }
               // });

                // self.collection.parseSearch(this._options.search, {
                    // success: function(search_full_string_encoded, searchHash){
                        // self.collection.setSearch(search_full_string_encoded, searchHash);
                        // //fetchOptions['data'] = searchHash;
                    // },
                    // error: function(search_full_string_encoded, msg){
                        // console.log("ERROR: " + search_full_string_encoded + ', parse failed with message: ' + msg);
                    // }
                // });


                // if(!_.isEmpty(searchExpressions)){
                    // self.collection.setSearch(this.
                    // self.collection.searchBy = this._options.search;
                    // console.log('trigger search:' + this._options.search + ', ' + JSON.stringify(searchExpressions) );
                    // self.collection.trigger("MyServerSideFilter:search", searchExpressions, self.collection);
                    // // console.log('done: trigger search');
                // }else{
                    // console.log('Warn: no search expressions found in: '  + this._options.search );
                // }
            }
            console.log('collection fetch trigger');
            var fetchOptions = { reset: true };
            self.collection.fetch(fetchOptions);



            return this;
        }

    });

  return ListView;
});