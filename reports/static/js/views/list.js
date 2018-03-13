define([
  'jquery',
  'underscore',
  'backbone',
//  'backbone_pageable',
  'backgrid',
  'iccbl_backgrid',
  'templates/rows-per-page.html',
  'templates/list.html',
  'templates/modal_ok_cancel.html',
], function($, _, Backbone, Backgrid,  Iccbl, 
	rowsPerPageTemplate, listTemplate, modalTemplate ){

//  // for compatibility with require.js, attach PageableCollection in the right 
//// place on the Backbone object
//  // see https://github.com/wyuenho/backbone-pageable/issues/62
//  Backbone.PageableCollection = BackbonePageableCollection;

  var ajaxStart = function(){
      $('#loading').fadeIn({duration:100});
  };
//  var ajaxComplete = function(){
//      $('#loading').fadeOut({duration:100});
//  };
//  $(document).bind("ajaxComplete", function(){
//      ajaxComplete(); // TODO: bind this closer to the collection
//  });

  var ListView = Backbone.View.extend({

    initialize : function(attributes, options) {
      console.log('initialize ListView: ');
      var self = this;
      Iccbl.assert( !_.isUndefined(options.resource), 
      		'listView options: "resource" is required');

      Iccbl.assert( !_.isUndefined(options.router), 
      		'listView options: router is required');
      Iccbl.assert( !_.isUndefined(options.url), 
      		'listView options: url is required');
      Iccbl.assert( !_.isUndefined(options.header_message), 
      		'listView options: header_message is required');
      Iccbl.assert( !_.isUndefined(options.title), 
      		'listView options: title is required');

      this.router = options.router;
      this._options = options;

      var ListModel = Backbone.Model.extend({
        defaults: {
            rpp: 25,
            page: 1,
            order: {},
            search: {}}
        });
      var current_options = _.clone(this.model.get('current_options'));
      var listModel = this.listModel = new ListModel({
        rpp: current_options.rpp,
        page: current_options.page,
        order: _.clone(current_options.order),
        search: _.clone(current_options.search)
      });

      this.objects_to_destroy = _([]);

      var _state = {
        currentPage: parseInt(self.listModel.get('page')),
        pageSize: parseInt(self.listModel.get('rpp'))
      };

      var orderHash = self.listModel.get('order');
      if(!_.isEmpty(orderHash)){
        _.each(_.keys(orderHash), function(key) {
            var dir = orderHash[key];
            var direction = 'ascending';
            var order = -1;
            // according to the docs, -1 == ascending
            if (dir === '-') {
                // according to the docs, 1 == descending
                direction = 'descending';
                order = 1;
            }
            _state['sortKey'] = key;
            _state['order'] = order;
        });
      }

      var Collection = Iccbl.MyCollection.extend({
      	state: _state,
      	url: this._options.url,
      	listModel: listModel
      });
      var collection = self.collection = new Collection();
//      var collection = self.collection = new Collection({
//        'url': this._options.url,
//        listModel: listModel
//      });
      this.objects_to_destroy.push(collection);

      this.listenTo(this.listModel, 'change:search', function(){
        var searchHash = self.listModel.get('search')
        var current_options = _.clone(self.model.get('current_options'));
        console.log('===--- list detect: listModel change:search old: ' 
        		+ JSON.stringify(current_options.search) 
        		+ ', ' + JSON.stringify(searchHash));
        current_options.search = searchHash;
        self.model.set({current_options: current_options });
      });

      this.listenTo(this.listModel, 'change:order', function(){
        var orderHash = self.listModel.get('order');
        console.log('===--- list detect: listModel change:order:' +
        		JSON.stringify(orderHash));
        var current_options = _.clone(self.model.get('current_options'));
        current_options.order = orderHash;
        self.model.set({current_options: current_options });
      });

      this.listenTo(this.listModel, 'change:rpp', function(){
        console.log('===--- list detect: listModel change:rpp');
        var pageSize = parseInt(self.listModel.get('rpp'));
        var current_options = _.clone(self.model.get('current_options'));
        current_options.rpp = pageSize;
        self.model.set({current_options: current_options });
      });

      this.listenTo(this.listModel, 'change:page', function(){
        console.log('===--- list detect: listModel change:page');
        var page = parseInt(self.listModel.get('page'));
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

      this.buildGrid(this._options.resource);
    },

    buildGrid : function(resource) {

      console.log('buildGrid...');
      var self = this;

      self.listenTo(self.collection, "MyCollection:link", 
  		function (model, column) {
          console.log('---- process link for '+ column);

          var fieldDef = resource.fields[column];
          if( _.has(fieldDef,'backgrid_cell_options')) {
              // NOTE: format for backgrid cell options is "/{attribute_key}/"
              backgrid_cell_options = fieldDef['backgrid_cell_options'];
              console.log('backgrid_cell_options: ' + backgrid_cell_options);

              _route = backgrid_cell_options.replace(/{([^}]+)}/g, 
          		function (match) {
                  	console.log('matched: ' + match + ', model: ' + model);
                  	match = match.replace(/[{}]/g,'');
                  	console.log('matched: ' + match + ', model: ' 
                  			+ model.get(match));
                  	return !_.isUndefined(model.get(match)) ? model.get(match) : match;
//                        	return typeof model.get(match) != "undefined" ? model.get(match) : match;
              	});
              console.log('route: ' + _route);
              this.router.navigate(_route, {trigger:true});
          }else{
              console.log('no options defined for link cell');
          }
      });

      self.listenTo(
        self.collection, "MyCollection:edit", 
    		function (model) {
          var id = Iccbl.getIdFromIdAttribute( model, resource );
          // signal to the app_model that the current view has changed 
          // todo: separate out app_model from list_model
          this.model.set({    current_scratch: { resource: resource, model: model} ,
                              current_view: 'edit',
                              current_options: { key: id },
                              routing_options: {trigger: false, replace: false}
                         }); 
      });

      self.listenTo(self.collection, "MyCollection:detail", function (model) {
        var id = Iccbl.getIdFromIdAttribute( model, resource );

          this.model.set({    current_scratch: { resource: resource, model: model} ,
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
      var rppModel = self.rppModel = new Backbone.Model({ selection: String(self.listModel.get('rpp')) });
      var rppSelectorInstance = self.rppSelectorInstance = new Iccbl.GenericSelector(
          { model: rppModel }, {label: 'Rows per page:', options: ['', '25','50','200','1000'], selectorClass: 'input-small' } );
      this.objects_to_destroy.push(rppSelectorInstance);
      this.listenTo(this.listModel, 'change: rpp', function(){
          rppModel.set({ selection: String(self.listModel.get('rpp')) });
      });
      this.listenTo(rppModel, 'change', function() {
          console.log('===--- rppModel change');
          self.listModel.set('rpp', String(rppModel.get('selection')));
      });

      var paginator = self.paginator = new Backgrid.Extension.Paginator({

      	  // If you anticipate a large number of pages, you can adjust
      	  // the number of page handles to show. The sliding window
      	  // will automatically show the next set of page handles when
      	  // you click next at the end of a window.
      	  // windowSize: 20, // Default is 10

      	  // Used to multiple windowSize to yield a number of pages to slide,
      	  // in the case the number is 5
      	  //slideScale: 0.25, // Default is 0.5

      	  // Whether sorting should go back to the first page
      	  // from https://github.com/wyuenho/backgrid/issues/432
      	  goBackFirstOnSort: false, // Default is true

      	  collection: self.collection
      	});            
      this.objects_to_destroy.push(paginator);

      // Extraselector
      if( _.has(resource, 'extraSelectorOptions')){
          var searchHash = self.listModel.get('search');
          console.log('extraselector init: searchTerms: ' + JSON.stringify(searchHash));

          var extraSelectorModel = new Backbone.Model({ selection: '' });
          var extraSelectorKey = resource.extraSelectorOptions.searchColumn;
          _.each(_.keys(searchHash), function(key){
              console.log('key: ' + key + ', extrSelectorKey: ' + extraSelectorKey);
              if( key == extraSelectorKey || key  === extraSelectorKey+ '__exact'){
                  extraSelectorModel.set({ selection: searchHash[key] });
              }
          });
          var extraSelectorInstance = self.extraSelectorInstance =
              new Iccbl.GenericSelector({ model: extraSelectorModel }, 
                                          resource.extraSelectorOptions );
          this.objects_to_destroy.push(extraSelectorInstance);

          this.listenTo(this.listModel, 'change: search', function(){
              var searchHash = self.listModel.get('search');
              console.log('extraselector, search changed: ' + JSON.stringify(searchHash));
              _.each(_.keys(searchHash), function(key){
                  console.log('key: ' + key + ', extrSelectorKey: ' + extraSelectorKey);
                  if( key === extraSelectorKey || key  === extraSelectorKey+ '__exact'){
                      extraSelectorModel.set({ selection: searchHash[key] });
                  }
              });
          });
          this.listenTo(extraSelectorModel, 'change', function() {
              console.log('===--- extraSelectorModel change');
              var searchHash = _.clone(self.listModel.get('search'));
              var value = extraSelectorModel.get('selection');
              searchHash[extraSelectorKey + '__exact'] = value;
              self.listModel.set('search', searchHash);
              self.collection.setSearch(searchHash);
          });
      }

      var columns = Iccbl.createBackgridColModel(
        this._options.resource.fields, [],{}, Iccbl.MyHeaderCell);

      var grid = this.grid = new Backgrid.Grid({
        columns: columns,
        collection: self.collection,
      });

      this.objects_to_destroy.push(grid);

      // encapsulate the footer in a view, help grab button click
      var footer = self.footer = new Backbone.View({
          el: $("<form><button type='button' id='addRecord'>Add</button></form>"),
          events: {
              'click button':function(event) {
                  console.log('button click event, '); 
                  event.preventDefault();
                  // TODO: set the defaults, also determine if should be 
                  // set on create, from the Meta Hash
                  var defaults = {};

                  id_attributes = self._options.resource['id_attribute']
                  _.each(resource.fields, function(value, key){
                      if (key == 'resource_uri') {
                          defaults[key] = self._options.url;
                      } else if (key == 'id'){
                          // nop // TODO: using the meta-hash, always exclude the primary key from create
                      // } else if (_.contains(id_attributes,key)){
                          // // nop // TODO: using the meta-hash, always exclude the primary key from create
                      } else {
                           defaults[key] = '';
                      }
                  });
                  var NewModel = Backbone.Model.extend({urlRoot: self._options.url, defaults: defaults });

                  self.collection.trigger('MyCollection:edit', new NewModel());

                  // // TODO: get the model "edit title" from the metainformation_hash
                  // var detailView = new DetailView(
                      // { model: new NewModel},
                      // { isEditMode: true, title: "Add new record", resource:resource, router:self._options.router});
//
                  // $('#list-container').hide();
                  // // NOTE: having self bind to the detailView like this:
                  // // self.listenTo(detailView, 'remove', function(){
                  // // causes the detailView to hang around in memory until self is closed
                  // // detailView.on('remove', function(){
                  // self.listenToOnce(detailView, 'remove', function(){
                      // console.log('... remove detailView event');
                      // self.collection.fetch({reset:true});
                      // $('#list-container').show();
                      // detailView.close();
                  // });
                  // $('#detail-container').append(detailView.render().$el);

              },
          },

      });
      this.objects_to_destroy.push(footer);

      // Note: prefer listenTo over "on" (alias for _.bind/model.bind) as this
      // will allow the object to unbind all observers at once.
      //collection.on('request', ajaxStart); // NOTE: can use bind or on
      //collection.bind('sync', ajaxComplete, this);

      this.listenTo(self.collection, 'request', ajaxStart);

      // TODO: work out the specifics of communication complete event.  
      // the following are superceded by the global handler for "ajaxComplete"
//      this.listenTo(self.collection, 'error', ajaxComplete);
//      this.listenTo(self.collection, 'complete', ajaxComplete);
      console.log('list view initialized');
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
      console.log('--render start');
      var self = this;
      self.listenTo(self.collection, "add", self.reportState);
      self.listenTo(self.collection, "remove", self.reportState);
      self.listenTo(self.collection, "reset", self.reportState);


      this.$el.html(this.compiledTemplate);
      var finalGrid = self.finalGrid = this.grid.render();
      self.$("#example-table").append(finalGrid.$el);
      self.$("#paginator-div").append(self.paginator.render().$el);
      self.$("#rows-selector-div").append(
      		self.rppSelectorInstance.render().$el);
      if(!_.isUndefined(self.extraSelectorInstance)){
        self.$("#extra-selector-div").append(
        		self.extraSelectorInstance.render().$el);
      }
      self.$("#table-footer-div").append(self.footer.$el);

      this.delegateEvents();
      
      var fetched = false;
      
      var searchHash = self.listModel.get('search');
      if(!_.isEmpty(searchHash)){
        self.collection.setSearch(searchHash);
        fetched = true;
      }

      var orderHash = self.listModel.get('order');
      if(!_.isEmpty(orderHash)){
        _.each(_.keys(orderHash), function(key) {
          var dir = orderHash[key];
          var direction = 'ascending';
          var order = -1;
          // according to the docs, -1 == ascending
          if (dir === '-') {
              // according to the docs, 1 == descending
              direction = 'descending';
              order = 1;
          }
          finalGrid.sort(key, direction);
          fetched = true;
        });
      }

      this.listenTo(self.collection, 'sync', function(event){
        self.$('#header_message').html(self._options.header_message + 
        		", total records: " + self.collection.state.totalRecords);
      });
      
      if ( !fetched ) {
        var fetchOptions = { reset: false };
        self.collection.fetch(fetchOptions);
      }

      console.log('rendered');
      return this;
    },
    
    reportState: function(){
    	var self = this;
    	var state = self.collection.state;
      var currentPage = Math.max(state.currentPage, state.firstPage);

      // Order: note, single sort only at this time
      var orderHash = {};
      if(state.order && state.order != 0 && state.sortKey ){
          // Note: 
      	// backbone-pageable: state.order: ascending=-1, descending=1
      	// tastypie: "-"=descending, ""=ascending (not specified==ascending)
      	orderHash[state.sortKey] = state.order == -1 ? '' : '-';
      }
      
      self.listModel.set({ 'rpp': state.pageSize, 
  						 'page': currentPage,
  						 'order': orderHash });
    }

  });

  return ListView;
});