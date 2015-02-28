define([
  'jquery',
  'underscore',
  'backbone',
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

  
  var ListView = Backbone.View.extend({

    LIST_ROUTE_ORDER: ['rpp', 'page', 'includes', 'order','search'],
    SEARCH_DELIMITER: ';',
    
    events: {
      'click .btn#select_columns': 'select_columns',
      'click .btn#download_link': 'download'
    },

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
      var orderings = [];
      var includes = [];
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
              var searches = value.split(this.SEARCH_DELIMITER);
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
            } else if (key === 'includes') {
              includes = value.split(',');
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
      
      if(!_.isUndefined(resource.options)
          && ! _.isUndefined(resource.options.includes)){
        if(!_.isEmpty(resource.options.order)){
          includes = _.union(includes, resource.options.includes);
        }
      }
      
      listInitial['includes'] = includes;

      if(!_.isUndefined(resource.options)
          && ! _.isUndefined(resource.options.rpp)){
        listInitial['rpp'] = resource.options.rpp;
      }

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
            this._options.schemaResult.fields, 
            Iccbl.MyHeaderCell, 
            orderStack,
            includes);
      }else{
        columns = _options.columns;
      }
      
      this.listenTo(this.listModel, 'change', this.reportState );
      this.buildGrid( columns, self._options.schemaResult );
    },
    
    getCollectionUrl: function(limit) {
      var self = this;
      var urlparams = '';
      self.trigger('update_title', '');
      
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
            if(!_.isEmpty(urlparams)) urlparams += '&';
            urlparams += 'order_by=' + value.join('&order_by=');
          }else if (route === 'includes') {
            if(!_.isEmpty(urlparams)) urlparams += '&';
            urlparams += 'includes=' + value.join('&includes=');
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
              if (val !== '' ) val += self.SEARCH_DELIMITER;
              
              val += k + "=" + v;
            });
            
            newStack.push(val);
            
          }else if (route === 'order') {
            newStack.push(value);
          }else if (route === 'includes') {
            newStack.push(value.join(','));
          } else {
            newStack.push(value);
          }
        }
      });
      console.log('newStack: ' + JSON.stringify(newStack));
      self.trigger('uriStack:change', newStack );
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
      
      var rpp_selections = ['25','50','200','500','1000'];
      if(!_.isUndefined(self._options.resource.options)
          && ! _.isUndefined(self._options.resource.options.rpp_selections)){
        rpp_selections = self._options.resource.options.rpp_selections;
      }
      var rppSelectorInstance = self.rppSelectorInstance = new genericSelector(
          { model: rppModel }, 
          { label: 'Rows', 
            options: rpp_selections });
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
            self.collection.fetch();
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
      console.log('list view initialized');
    },


    remove: function(){
      console.log('ListView remove called');
      Backbone.View.prototype.remove.apply(this,arguments);
    },

    /** Backbone.layoutmanager callback **/
    cleanup: function(){
      this.onClose();
    },
    onClose: function(){
      this.$el.empty();
      
      if(_.isObject(this.objects_to_destroy)){
          this.objects_to_destroy.each(function(view_obj){
            console.log('destroy: ' + view_obj);
              view_obj.remove();
              view_obj.off();
              view_obj.stopListening();
              view_obj = null;
          });
      }

      this.objects_to_destroy = null; // critical for memory mgt
      this.paginator = null; // critical for memory mgt

      this.finalGrid = null;
      this.grid = null;
      this.collection = null;
      this.extraSelectorInstance = null;
      this.rppSelectorInstance = null;
      this.footer = null;
      this.listModel = null;
      this._options = null;
      this.off();
    },
    

    afterRender: function(){
      var self = this;
      self.listenTo(self.collection, "add", self.checkState);
      self.listenTo(self.collection, "remove", self.checkState);
      self.listenTo(self.collection, "reset", self.checkState);
      self.listenTo(self.collection, "sort", self.checkState);

      this.$el.html(_.template(listTemplate));
      var finalGrid = self.finalGrid = this.grid.render();
      self.objects_to_destroy.push(finalGrid);
      self.$("#example-table").append(finalGrid.el);
      // FIXME: move css classes out to templates
      finalGrid.$el.addClass("col-sm-12 table-striped table-condensed table-hover");

      self.$("#paginator-div").append(self.paginator.render().$el);
      self.$("#rppselector").html(
      		self.rppSelectorInstance.render().$el);
      if(!_.isUndefined(self.extraSelectorInstance)){
        self.$("#extraselector").html(
            self.extraSelectorInstance.render().$el);
      }

      var clearSearchesButton = $([
        '<div class="pull-right pull-down">',
        '<div class="col-xs-1">',
        '<a id="clear-searches" >clear searches</a>',
        '</div>',
        '</div>'
        ].join(''));
      clearSearchesButton.click(function(e){
        e.preventDefault();
        self.clear_searches();
      });
      self.$('#list-header').append(clearSearchesButton);
      
      var clearSortsButton = $([
        '<div class="pull-right pull-down">',
        '<div class="col-xs-1">',
        '<a id="clear-sorts" >clear sorts</a>',
        '</div>',
        '</div>'
        ].join(''));
      clearSortsButton.click(function(e){
        e.preventDefault();
        self.clear_sorts();
      });
      self.$('#list-header').append(clearSortsButton);
                               
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
      this.reportState();
      return this;
    },
    
    clear_searches: function(){
      this.collection.trigger("Iccbl:clearSearches");
      this.collection.getFirstPage({reset: true, fetch: true});
    },
    
    clear_sorts: function(){
      this.collection.trigger("Iccbl:clearSearches");
      this.listModel.unset('order');
      this.collection.state.orderStack = [];
      this.collection.getFirstPage({reset: true, fetch: true});
    },
    
    download: function(e){
      var self = this;
      e.preventDefault();
      var limitForDownload = 0;
      appModel.download(self.getCollectionUrl(limitForDownload), self._options.resource);
    },
     
    select_columns: function(event){
      var self = this;
      var includes = self.listModel.get('includes') || [];
      var altCheckboxTemplate =  _.template('\
              <div class="form-group" style="margin-bottom: 0px;" > \
                <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
                  <label for="<%= editorId %>"><span data-editor\><%= title %></label>\
                </div>\
              </div>\
            ');

      var formSchema = {};
      _.each(_.pairs(this._options.schemaResult.fields), function(pair){
        var prop = pair[1];
        var key = prop['key'];
        if(key == 'resource_uri' || key == 'id') return;

        var title = prop['title'];
        formSchema[key] = { 
            title: title, 
            ordinal: prop['ordinal'],
            key:  key, 
            type: 'Checkbox',
            template: altCheckboxTemplate };
      });
      var FormFields = Backbone.Model.extend({
        schema: formSchema
      });
      
      var formFields = new FormFields();
      var already_visible = {};
      var default_visible = {};
      _.each(_.pairs(this._options.schemaResult.fields), function(pair){
        var key = pair[1]['key'];
        var prop = pair[1];
        if(key == 'resource_uri' || key == 'id') return;
        
        var _visible = (_.has(prop, 'visibility') && 
            _.contains(prop['visibility'], 'list'));
        default_visible[key] = _visible;
        _visible = _visible || _.contains(includes, key);
        _visible = _visible && !_.contains(includes, '-'+key);
        
        already_visible[key] = _visible;
        formFields.set( key, _visible);
      });
      
      var orderedFieldKeys = _.sortBy(formFields.keys(), function(key){
        console.log('find:' + key + ',' + formSchema[key]);
        return formSchema[key]['ordinal'];
      });
      console.log('orderedFields: ' + JSON.stringify(orderedFieldKeys));
  
      var form = new Backbone.Form({
          model: formFields,
          fields: orderedFieldKeys,
          template: _.template([
            "<form data-fieldsets class='form-horizontal container' >",
            "<div class='btn btn-default btn-sm ' id='select-all' >select all</div>",
            "<div class='btn btn-default btn-sm ' id='clear-all' >clear all</div>",
            "</form>"
            ].join(''))
      });
      
      form.events = {
        'click .btn#select-all': function(){
          $("form input:checkbox").each(function(){
            $(this).prop("checked",true);
          });
        },
        'click .btn#clear-all': function(){
          $("form input:checkbox").each(function(){
            $(this).prop("checked",false);
          });
        }
      };

      appModel.showModal({
        ok: function(){
          form.commit();
          if(_.isUndefined(
              _.find(formFields.values(),function(val){ return val==true; }))){
            // TODO: display "nothing selected" error
            console.log('error: nothing selected');
            self.select_columns();
            return;
          }
          
          var new_includes = [];
          console.log('formFields: ' + JSON.stringify(formFields.toJSON()));
          _.each(formFields.keys(), function(key){
            var value = formFields.get(key);
            
            if(value && !already_visible[key] ){
              self.grid.insertColumn(
                  Iccbl.createBackgridColumn(
                      key,self._options.schemaResult.fields[key],
                      Iccbl.MyHeaderCell, self.collection.state.orderStack));
            }
            if(!value && default_visible[key]){
              new_includes.unshift('-' + key);
              column =  self.grid.columns.find(function(column){
                if(column.get('name') == key){
                  self.grid.removeColumn(column);
                  return true;
                }
              });
            }else if(!value){
              column =  self.grid.columns.find(function(column){
                if(column.get('name') == key){
                  self.grid.removeColumn(column);
                  return true;
                }
              });
            }
            if(value && !default_visible[key]){
              new_includes.unshift(key);
            }
          });
          
          console.log('new_includes: ' + JSON.stringify(new_includes));
          self.listModel.set({'includes': new_includes });
          self.collection.fetch();
          self.listModel.set({'includes': new_includes });
        },
        view: form.render().el,
        title: 'Select columns'  
      });
      
    },
    
  });
  

  return ListView;
});
