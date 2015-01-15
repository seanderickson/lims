// REFACTOR OF list.js to use layoutmanager
define([
  'jquery',
  'underscore',
  'backbone',
//  'backbone_pageable',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_selector',
  'text!templates/rows-per-page.html',
  'text!templates/list2.html',
  'text!templates/modal_ok_cancel.html'
], function(
      $, _, Backbone, Backgrid,  
      Iccbl, appModel, genericSelector,
      rowsPerPageTemplate, listTemplate, modalTemplate) {

//  var ajaxStart = function(){
//      $('#loading').fadeIn({duration:100});
//  };
//  var ajaxComplete = function(){
//      $('#loading').fadeOut({duration:100});
//  };
//  $(document).bind("ajaxComplete", function(){
//      ajaxComplete(); 
//  });

  var ListView = Backbone.View.extend({
    LIST_ROUTE_ORDER: ['rpp', 'page','order','search'],
    
    initialize : function(args) {
      console.log('initialize ListView: ');
      var self = this;
      var _options = self._options = args.options;
      var resource = self._options.resource;
      
      var ListModel = Backbone.Model.extend({
        defaults: {
            rpp: 25,
            page: 1,
            order: {},
            search: {}}
        });
      

      // grab any preset searches
      var preset_searches = {};
      if(!_.isUndefined(resource.options)
          && ! _.isUndefined(resource.options.search)){
        _.each(_.keys(resource.options.search), function(key){
          preset_searches[key] = resource.options.search[key];
        });
      }

      // convert the uriStack into the listmodel
      var urlSuffix = self.urlSuffix = "";
      var listInitial = {};
      var searchHash = {};
      var orderings = []
      if(_.has(self._options,'uriStack')){
        var stack = self._options.uriStack;
        for (var i=0; i<stack.length; i++){
          var key = stack[i];
          i = i+1;
          if (i==stack.length) continue;
          var value = stack[i];
          if (!value || _.isEmpty(value) ){
            continue;
          }
          
          if(key == 'log'){
            self.urlSuffix = key + '/' + value;
            continue;
          }
          if(key == 'children'){
            // This is a hack to show the children of an apilog, see 
            // reports.api.ApiLogResource.prepend_urls for further details
            var substack = _.rest(stack,i)
            var substack_consumed = []
            var _key = Iccbl.popKeyFromStack(resource, substack, substack_consumed);
            i += substack_consumed.length;
            self.urlSuffix = key + '/' + _key;
            console.log('urlSuffix: ' + self.urlSuffix);
          }
          
          if(_.contains(this.LIST_ROUTE_ORDER, key)){
            
            if(key === 'search') {
              var searches = value.split(',');
              _.each(searches, function(search){
                var parts = search.split('=');
                if (!parts || parts.length!=2) {
                  var msg = 'invalid search parts: ' + search;
                  console.log(msg);
                  appModel.error(msg);
                } else if (_.isEmpty(parts[1])) {
                  // pass, TODO: prevent empty searches from notifying
                } else {
                  searchHash[parts[0]] = parts[1];
                }
              });
            } else if (key === 'order') {
              orderings = value.split(',');
              //listInitial[key] = value.split(',');        
            }else {
              listInitial[key] = value;
            }
          }
        }
      }
      // TODO: do the same thing here to implement preset orderings
      searchHash = _.extend({}, preset_searches,searchHash);
      if(!_.isEmpty(searchHash)){
        listInitial['search'] = searchHash;
      }
      
      if(!_.isUndefined(resource.options)
          && ! _.isUndefined(resource.options.order)){
        if(!_.isEmpty(resource.options.order)){
          // TODO: to complicated: should deal with orderings as a hash, to simplify
          orderings = _.union(orderings,resource.options.order);
          var orderHash = {};
          _.each(orderings, function(item){
            var dir = item.substring(0,1);
            var fieldname = item;
            if(dir == '-'){
              fieldname = item.substring(1);
            }else{
              dir = '';
            }
            if(!_.has(orderHash,fieldname)){
              orderHash[fieldname] = dir;
            }else{
              orderings = _.without(orderings, item);
            }
          });
        }
      }
      listInitial['order'] = orderings;
      
      var listModel = this.listModel = new ListModel(listInitial);

      this.objects_to_destroy = _([]);

      var _state = {
        currentPage: parseInt(self.listModel.get('page')),
        pageSize: parseInt(self.listModel.get('rpp'))
      };

      var orderStack = self.listModel.get('order') || [];
      _state.orderStack = _.clone(orderStack);

      var url = self._options.resource.apiUri + '/' + self.urlSuffix;
      if (_.has(self._options, 'url')) {
        url = self._options.url;
      } else {
        self._options.url = url;   // TODO: cleanup messy initializations
      }

      var collection;
      if( !_options.collection){  
        var Collection = Iccbl.MyCollection.extend({
          state: _state  
        });
        
        collection = self.collection = new Collection({
          'url': url,
          listModel: listModel
        });
        this.objects_to_destroy.push(collection);
      }else{
        collection = self.collection = _options.collection;
        collection.listModel = listModel;
        collection.state = _state;
      }
      
      var columns;
      if(!_options.columns){
        columns = Iccbl.createBackgridColModel(
            this._options.schemaResult.fields, Iccbl.MyHeaderCell, orderStack);
      }else{
        columns = _options.columns;
      }
      
      // ==== testing - images
      // FIXME: FOR SMR only - poc - 
      // first find the image field
      _.each(_.keys(this._options.schemaResult.fields), function(key){
        var field = self._options.schemaResult.fields[key];
        if (field['ui_type'] == 'image'){
          var col = {
              'name' : field['key'],
              'label' : field['title'],
              'description' : field['description'],
              'backgrid_cell_options': field['backgrid_cell_options'],
              cell : Iccbl.ImageCell,
              order : field['ordinal'],
              editable : false,
          };
          console.log('adding image col: ' + field['key'] + ', ' );
          columns.unshift(col);
        }
      });
      
      // end - FIXME: FOR SMR only - poc - 
      //////////
      

      this.listenTo(this.listModel, 'change', this.reportState );

      var compiledTemplate = this.compiledTemplate = _.template(listTemplate);

      this.buildGrid( columns, self._options.schemaResult );
    },
    
    getCollectionUrl: function(limit) {
      var self = this;
      var urlparams = '';
      _.each(self.LIST_ROUTE_ORDER, function(route){
        var value = self.listModel.get(route);
        if ( (!_.isObject(value) && value ) || 
             ( _.isObject(value) && !_.isEmpty(value))) {
          
          if (route === 'search') {
            var val = '';
            _.each(value, function(v,k){
              if (val !== '' ) val += '&';
              
              val += k + "=" + v;
            });
            self.trigger('update_title', 'Search: ' + val);
            if(!_.isEmpty(urlparams)) urlparams += '&';
            urlparams += val;
          }else if (route === 'order') {
            urlparams += 'order_by=' + value.join('&order_by=');
          }
        }
      });
      if(_.isUndefined(limit)){
        var limit = self.collection.state.totalRecords;
        if(!_.isNumber(limit)){
          limit = 0;
        }
      }
      
      var url = self.collection.url + '?' + 
                '&limit=' + limit;
      if(!_.isEmpty(urlparams)) url += '&' + urlparams;
      console.log('collection url: ' + url)
      return url;
      
    },
    
    reportState: function() {
      var self = this;
      var newStack = [];
      
      // If a suffix was consumed, then put it back
      if(self.urlSuffix != ""){
        newStack = self.urlSuffix.split('/');
      }
      
      _.each(self.LIST_ROUTE_ORDER, function(route){
        var value = self.listModel.get(route);
        if ( (!_.isObject(value) && value ) || 
             ( _.isObject(value) && !_.isEmpty(value))) {
          newStack.push(route);
          
          if (route === 'search') {
            var val = '';
            _.each(value, function(v,k){
              if (val !== '' ) val += ',';
              
              val += k + "=" + v;
            });
            
            newStack.push(val);
            
          }else if (route === 'order') {
            newStack.push(value);
          } else {
            newStack.push(value);
          }
        }
      });
      self.trigger('uriStack:change', newStack );

      // FIXME: TODO: see reports.ManagedResource.create_response:
      // We need to set the "Content-Disposition" header to trigger the server to 
      // bounce it back 
      var limitForDownload = 0;
      $('#download_link').attr('href', self.getCollectionUrl(limitForDownload) + '&format=csv');
    },
    
    checkState: function(){
      var self = this;
      var state = self.collection.state;
      var currentPage = Math.max(state.currentPage, state.firstPage);

      
      // search: set in Iccbl Collection
      
      self.listModel.set({ 
        'rpp': state.pageSize, 
        'page': currentPage,
        'order': _.clone(state.orderStack) 
      });
      
    },

    buildGrid: function( columns, schemaResult ) {

      console.log( 'buildGrid...');
      var self = this;

      self.listenTo(self.collection, "MyCollection:link", 
    		function(model, column) {
          console.log('---- process link for '+ column);
  
          var fieldDef = schemaResult.fields[column];
          if( _.has(fieldDef,'backgrid_cell_options')) {
            // NOTE: format for backgrid cell options is "/{attribute_key}/"
            var backgrid_cell_options = fieldDef['backgrid_cell_options'];
            var _route = Iccbl.replaceTokens(model,backgrid_cell_options);
            console.log('route: ' + _route);
            appModel.router.navigate(_route, {trigger:true});
          }else{
            console.log('no options defined for link cell');
          }
        });

      self.listenTo(
          self.collection, "MyCollection:detail", 
          function (model) {
            console.log('detail handler for' + model.get('toString'));
            if(!_.isUndefined(self._options.detailHandler)){
              self._options.detailHandler(model);
            }else{
//            self.trigger('detail', model);
              var id = Iccbl.getIdFromIdAttribute( model, schemaResult );
              model.resource = self._options.resource;
              appModel.router.navigate(model.resource.key + '/' + id, {trigger:true});
            }
          });

      // TODO: not used
      self.listenTo(self.collection, "MyCollection:delete", function (model) {
          var modalDialog = new Backbone.View({
              el: _.template(modalTemplate, { 
                  body: "Please confirm deletion of record: '" + model.get('toString') + 
                        "'", title: "Please confirm deletion" } ),
              events: {
                  'click #modal-cancel':function(event) {
                      console.log('cancel button click event, '); 
                      event.preventDefault();
                      // TODO: read-up on modal!  this is not ideal with the 
                      // reference to template elements!
                      $('#modal').modal('hide'); 
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
      var rppModel = self.rppModel = new Backbone.Model({ 
          selection: String(self.listModel.get('rpp')) 
        });
      var rppSelectorInstance = self.rppSelectorInstance = new genericSelector(
          { model: rppModel }, 
          { label: 'Rows', 
            options: ['25','50','200','1000'] });
      this.objects_to_destroy.push(rppSelectorInstance);
      this.listenTo(this.listModel, 'change: rpp', function(){
          rppModel.set({ selection: String(self.listModel.get('rpp')) });
      });
      this.listenTo(rppModel, 'change', function() {
          var rpp = parseInt(rppModel.get('selection'));
          self.listModel.set('rpp', String(rpp));

          console.log('===--- rppModel change: ' + rpp );
          
          self.collection.setPageSize(rpp);

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

    	  collection: self.collection,
    	  
    	});            
      this.objects_to_destroy.push(paginator);

      // Downloadselector
      var downloadSelectorModel = new Backbone.Model({ selection: '' });
      var downloadSelectorOptions = []
      if(_.has(self._options.resource, 'content_types')){
        // exclude JSON for downloads
        downloadSelectorOptions = _.without(self._options.resource.content_types, 'json');
      }
      if (!_.contains(downloadSelectorOptions, '')) {
        downloadSelectorOptions.unshift('');
        // create a blank entry
      }
      var downloadSelectorInstance = self.downloadSelectorInstance =
        new genericSelector({ model: downloadSelectorModel }, { options: downloadSelectorOptions } );
      this.objects_to_destroy.push(downloadSelectorInstance);
      this.listenTo(downloadSelectorModel, 'change', function() {
        var val = downloadSelectorModel.get('selection');
        var limitForDownload = 0;
        $('#download_link').attr('href', self.getCollectionUrl(limitForDownload) + '&format=' + val);
      });
      
      // Extraselector
      if( _.has(schemaResult, 'extraSelectorOptions')){
        var searchHash = self.listModel.get('search');
        console.log('extraselector init: searchTerms: ' + JSON.stringify(searchHash));

        var extraSelectorModel = new Backbone.Model({ selection: '' });
        var extraSelectorKey = schemaResult.extraSelectorOptions.searchColumn;
        
        if ( !_.isEmpty(searchHash)){
          _.each(_.keys(searchHash), function(key){
              console.log('key: ' + key + ', extrSelectorKey: ' + extraSelectorKey);
              if( key == extraSelectorKey 
                  || key  === extraSelectorKey+ '__exact'
                  || key  === extraSelectorKey+ '__ne'){
                  extraSelectorModel.set({ selection: searchHash[key] });
              }
          });
        }
        var selectorOptions = _.clone(schemaResult.extraSelectorOptions);
        if (!_.contains(selectorOptions.options, '')) {
          selectorOptions.options.unshift('');
          // create a blank entry
        }

        var extraSelectorInstance = self.extraSelectorInstance =
            new genericSelector({ model: extraSelectorModel }, selectorOptions );
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
            var searchHash = _.clone(self.listModel.get('search'));
            var val = extraSelectorModel.get('selection');
            var value =  _.isUndefined(val.value) ? val: val.value ;

            delete searchHash[extraSelectorKey + '__exact']
            delete searchHash[extraSelectorKey + '__ne']
            delete searchHash[extraSelectorKey]
            
            if(!_.isEmpty(value) && !_.isEmpty(value.trim())){
              var field = self._options.schemaResult.fields[extraSelectorKey];
              if(field.ui_type=='boolean'){
                searchHash[extraSelectorKey] = value;
              }else{
                // Convert the !negated values to value__ne 
                // (this is better for the URL)
                if(value.indexOf('!')==0){
                  searchHash[extraSelectorKey + '__ne'] = value.substring(1);
                }else{
                  searchHash[extraSelectorKey + '__exact'] = value;
                }
              }
            }
            self.listModel.set('search', searchHash);
            self.collection.setSearch(searchHash);
        });
      }

      var grid = this.grid = new Backgrid.Grid({
        columns: columns,
        body: Iccbl.MultiSortBody,
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
                  appModel.router.navigate(self._options.resource.api_resource + "/+add", {trigger:true});
                  return;
              },
          },
      });
      this.objects_to_destroy.push(footer);

      // Note: prefer listenTo over "on" (alias for _.bind/model.bind) as this
      // will allow the object to unbind all observers at once.
      //collection.on('request', ajaxStart); // NOTE: can use bind or on
      //collection.bind('sync', ajaxComplete, this);

//      this.listenTo(self.collection, 'request', ajaxStart);

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
    

    beforeRender: function(){
      var self = this;
      self.listenTo(self.collection, "add", self.checkState);
      self.listenTo(self.collection, "remove", self.checkState);
      self.listenTo(self.collection, "reset", self.checkState);
      self.listenTo(self.collection, "sort", self.checkState);

      this.$el.html(this.compiledTemplate);
      var finalGrid = self.finalGrid = this.grid.render();
      self.$("#example-table").append(finalGrid.$el);
      // FIXME: move css classes out to templates
      finalGrid.$el.addClass("col-sm-12 table-striped table-condensed table-hover");

      self.$("#paginator-div").append(self.paginator.render().$el);
//      self.$('#list-header').append(
//          '<div class="pull-right pull-down"><a class="btn btn-default btn-sm pull-down" role="button" id="download_link" href="' + 
//          self.getCollectionUrl() +
//          '">download</a></div>');
      self.$("#rppselector").html(
      		self.rppSelectorInstance.render().$el);
      if(!_.isUndefined(self.extraSelectorInstance)){
        self.$("#extraselector").html(
            self.extraSelectorInstance.render().$el);
      }
      if(!_.isUndefined(self.downloadSelectorInstance)){
        self.$("#downloadselector").html(
            self.downloadSelectorInstance.render().$el);
      }
              
//      // FIXME: "add" feature should be enabled declaratively, by user/group status:
//      if(appModel.getCurrentUser().is_superuser
//          && _.contains(self._options.resource.visibility, 'add')){
//        self.$("#table-footer-div").append(self.footer.$el);
//      }
      
      if(appModel.hasPermission(self._options.resource.key, 'write')){
        self.$("#table-footer-div").append(self.footer.$el);
      }

      this.delegateEvents();
      
      var fetched = false;
      
      var searchHash = self.listModel.get('search');
      if(!_.isEmpty(searchHash)){
        self.collection.setSearch(searchHash);
        fetched = true;
      }

      var orderStack = self.listModel.get('order') || [];
      self.collection.state.orderStack = _.clone(orderStack);
      // TODO: test further: added 20150114
      if(!_.isEmpty(orderStack)){
        self.collection.setSorting();
        //        fetched = true;
      }

      this.listenTo(self.collection, 'sync', function(event){
        var msg = ''; 
        if (self._options.header_message) {
          if (msg) msg += ', ';
          msg += self._options.header_message;
        }
        if (msg) msg += ', ';
        msg += 'Page ' + self.collection.state.currentPage + 
               ' of ' + ( self.collection.state.lastPage ? self.collection.state.lastPage:1 ) + 
               ' pages, ' + self.collection.state.totalRecords + 
               ' ' + self._options.resource.title  + ' records';
        self.$('#header_message').html(msg);
      });
      
      if ( !fetched ) {
        var fetchOptions = { reset: false, error: appModel.jqXHRerror };
        self.collection.fetch(fetchOptions);
      }
      
      // Note: replace: true - to suppress router history:
      // at this point, reportState is modifying the URL to show rpp, pagesSize, etc.
      appModel.set('routing_options', {replace: true});  
      return this;
    },
    
    afterRender: function() {
      this.reportState();
    },
    
    

  });

  return ListView;
});