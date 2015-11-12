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

    //LIST_ROUTE_ORDER: ['rpp', 'page', 'includes', 'order','search'],
    LIST_ROUTE_ORDER: appModel.LIST_ARGS,
    SEARCH_DELIMITER: appModel.SEARCH_DELIMITER,
    
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
      
      var urlSuffix = self.urlSuffix = "";
      var listInitial = {};

      // first set presets
      
      if(!_.isUndefined(resource.options)
          && ! _.isUndefined(resource.options.rpp)){
        listInitial['rpp'] = resource.options.rpp;
      }
      if(resource.options && !_.isEmpty(resource.options.order)){
        listInitial['order'] = _.clone(resource.options.order);
      }
      if(resource.options && !_.isEmpty(resource.options.includes)){
        listInitial['includes'] = _.clone(resource.options.includes);
      }
      if(resource.options && !_.isEmpty(resource.options.search)){
        listInitial['search'] = _.clone(resource.options.search);
      }
      
      // set with the uriStack values
      if(_.has(self._options,'uriStack')){
        var stack = self._options.uriStack;
        console.log('initialize ListView: ',stack);
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
              var searchHash = {};
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
              listInitial[key] = searchHash;
            } else if (key === 'order') {
              listInitial[key] = value.split(',');        
            } else if (key === 'includes') {
              listInitial[key] = value.split(',');        
            }else {
              listInitial[key] = value;
            }
          }
        }
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
            orderStack,
            listModel.get('includes'));
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
      
      var search_title_val = '';
      
      var _translate_sql_specifier = function(k){
        if(! k || _.isEmpty(k)) return k;
        
        var parts = k.split('__');
        var name = k;
        if(parts.length == 2){
          var field = self._options.schemaResult.fields[parts[0]];
          if(field){
            name = field.title;
          }
          name += ' ' + parts[1];
        }
        return name;
      };
      
      _.each(self.LIST_ROUTE_ORDER, function(route){
        var value = self.listModel.get(route);
        if ( (!_.isObject(value) && value ) || 
             ( _.isObject(value) && !_.isEmpty(value))) {

          if (route === 'search') {
            
            var val = '';
            if( route == 'search'){
              _.each(value, function(v,k){
                if(val !== '' ) val += '&';
                if(search_title_val !== '' ) search_title_val += ' AND ';
                val += k + "=" + v;
                search_title_val += _translate_sql_specifier(k) + " (" + v + ") ";
              });
            }
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
      
      if(_.has(self._options, 'search_data')){
        console.log('search data found on the self._options: ', self._options.search_data);
        var or_clauses = '';
        var search_list = self._options.search_data;
        _.each(search_list, function(hash){
          var temp = '';
          _.each(_.keys(hash), function(k){
            var v = hash[k];
            if (temp !== '' ) temp += ' AND ';
            temp += _translate_sql_specifier(k) + " (" + v + ") ";
          });
          if(or_clauses != '') or_clauses += "<br>OR ";
          or_clauses += temp
        });
        if(search_title_val !== '') search_title_val += '<br>AND ';
        search_title_val += '[ ' + or_clauses + ' ]';
      }            
      if(search_title_val !== ''){
        self.trigger('update_title', 'Search: ' + search_title_val);
      }
      
      return url;
    },
    
    reportState: function(args) {
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
      var previousStack = self.currentStack;
      self.currentStack = newStack;
      if(previousStack && _.isEqual(previousStack,newStack)){
        console.log('no new stack updates');
      }else{
//        var route_update = false;
//        var changedAttributes = self.listModel.changedAttributes();
//        var previousAttributes = self.listModel.previousAttributes();
//        _.each(_.keys(changedAttributes), function(key){
//          if(_.has(previousAttributes, key)){ 
//            route_update = true;
//            console.log('valid changed attribute', key, changedAttributes, previousAttributes);
//          }
//        });
//        
////        if(route_update){
////            // replace false (default) to create browser history
////          appModel.set('routing_options', {replace: true});  
////        }
//
        // calling this to update the title
        self.getCollectionUrl(0);
        self.trigger('uriStack:change', newStack );
      }
    },
    
    checkState: function(){
      var self = this;
      var state = self.collection.state;
      var currentPage = 1;
      if(state.currentPage){
        if(state.firstPage){
          currentPage = Math.max(state.currentPage, state.firstPage);    
        }else{
          currentPage = state.currentPage;
        }
      }
      
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

      
      self.listenTo(self.collection, "MyCollection:delete", function (model) {
        var form_template = [
           "<div class='form-horizontal container' id='delete_form' >",
           "<form data-fieldsets class='form-horizontal container' >",
           "</form>",
           "</div>"].join('');      
        var formSchema = {};
        formSchema['comments'] = {
          title: 'Comments',
          key: 'comments',
          validators: ['required'],
          type: 'TextArea'
        };
        var FormFields = Backbone.Model.extend({
          schema: formSchema,
        });
        var formFields = new FormFields();
        var form = new Backbone.Form({
          model: formFields,
          template: _.template(form_template)
        });
        var _form_el = form.render().el;
        var title = Iccbl.getTitleFromTitleAttribute(model,self._options.schemaResult);
        var dialog = appModel.showModal({
            okText: 'confirm',
            ok: function(e){
              var errors = form.commit({ validate: true }); // runs schema and model validation
              console.log('errors', errors);
              if(!_.isEmpty(errors) ){
                console.log('form errors, abort submit: ',errors);
                return false;
              }else{
                var values = form.getValue();
                console.log('form values', values);
                var comments = values['comments'];
                var headers = {};
                headers[appModel.HEADER_APILOG_COMMENT] = comments;
              
                e.preventDefault();
                model.collection = self.collection;
                // Backbone will only send DELETE if the model has an id
                model.set('id', Iccbl.getIdFromIdAttribute(model,self._options.schemaResult));
                model.destroy({
                  wait: true,
                  headers: headers,
                  success: function(model,response){
                    console.log('model removed successfully', model, response);
                  },
                  error: appModel.backboneFetchError
                });
              }
            },
            view: _form_el,
            title: 'Confirm deletion of "' + title + '"'  
        });
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
          self.listModel.set('page',1);
          self.collection.state.currentPage = 1; // set this because of how checkstate is triggered
          console.log('===--- rppModel change: ' + rpp );
          appModel.set('routing_options', {replace: false});  
          self.collection.setPageSize(rpp, { first: true });
      });

      if(self.collection instanceof Backbone.PageableCollection){
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
      }
      
      self.listenTo(self.listModel, 'change:search', function(){
        // TODO: this listener should be set in the collection initializer
        var searchHash = _.clone(self.listModel.get('search'));
        self.collection.setSearch(searchHash);

//        self.collection.fetch();
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

        self.listenTo(self.listModel, 'change:search', function(){
            var searchHash = _.clone(self.listModel.get('search'));
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
            
            console.log('extra selector change:' + value);

            delete searchHash[extraSelectorKey + '__exact']
            delete searchHash[extraSelectorKey + '__ne']
            delete searchHash[extraSelectorKey]
            
            if(!_.isEmpty(value) && !_.isEmpty(value.trim())){
              var field = self._options.schemaResult.fields[extraSelectorKey];
              if(field.data_type=='boolean'){
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
        });
      }
      var grid = this.grid = new Backgrid.Grid({
        columns: columns,
        body: Iccbl.MultiSortBody,
        collection: self.collection,
      });
      this.objects_to_destroy.push(grid);

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
      this.listModel = null;
      this._options = null;
      this.off();
    },
    
    beforeRender: function(){
      this.$el.html(_.template(listTemplate));
    },
      
    afterRender: function(){
      var self = this;
      self.listenTo(self.collection, "add", self.checkState);
      self.listenTo(self.collection, "remove", self.checkState);
      self.listenTo(self.collection, "reset", self.checkState);
      self.listenTo(self.collection, "sort", self.checkState);

      var finalGrid = self.finalGrid = this.grid.render();
      self.objects_to_destroy.push(finalGrid);
      self.$("#example-table").append(finalGrid.el);
      // FIXME: move css classes out to templates
      finalGrid.$el.addClass("col-sm-12 table-striped table-condensed table-hover");

      if(self.collection instanceof Backbone.PageableCollection){
        self.$("#paginator-div").append(self.paginator.render().$el);
        self.$("#rppselector").html(
            self.rppSelectorInstance.render().$el);
      }
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
                               
      if(_.has(self._options,'extraControls')){
        _.each(self._options.extraControls, function(control){
          self.$('#extra_controls').append(control);
        });
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
        var fetchOptions = { reset: false, error: appModel.backboneFetchError };
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
      this.collection.trigger("Iccbl:clearSorts");
      this.listModel.unset('order');
      this.collection.state.orderStack = [];
      this.collection.getFirstPage({reset: true, fetch: true});
    },
    
    download: function(e){
      var self = this;
      var limitForDownload = 0;

      e.preventDefault();
      
      // find any extra search data
      if(_.has(self._options, 'search_data')){
        console.log('search data found on the list._options: ', self._options.search_data);
        // endcode for the post_data arg; post data elements are sent 
        // as urlencoded values via a POST form for simplicity
        appModel.download(
            self.getCollectionUrl(limitForDownload), 
            self._options.resource, 
            { search_data: JSON.stringify(self._options.search_data) } ); 
      }else{
        appModel.download(
            self.getCollectionUrl(limitForDownload), 
            self._options.resource);
      }
    },
     
    /** Build the select columns dialog **/
    select_columns: function(event){
      
      var self = this;
      var form_template = [
        "<form  class='form-horizontal container' >",
        "<div class='btn btn-default btn-sm ' id='select-all' >select all</div>",
        "<div class='btn btn-default btn-sm ' id='clear-all' >clear all</div>"
      ];
      var field_template = '<div data-fields="<%= name %>" ></div>';
      var sub_field_template = '<div data-fields="<%= name %>" >   </div>';
      var header_template = [
        '<div class="form-group" >',
        '<input class="selection-group" type="checkbox" id="<%= id %>-checkbox"> </input>',
        '<label class="selection-group .h5 " id="<%= id %>" title="<%= help %>" ><%= name %> columns</label>',
        '</div>'
        ].join('');
      var altCheckboxTemplate =  _.template('\
          <div class="form-group" style="margin-bottom: 0px;" > \
            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
              <label title="<%= help %>" for="<%= editorId %>"><span data-editor\><%= title %></label>\
            </div>\
          </div>\
        ');      
      var altSubCheckboxTemplate =  _.template('\
          <div class="form-group sub-resource-field" style="margin-bottom: 0px;" > \
            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
            <label for="<%= editorId %>" > - </label>\
              <label for="<%= editorId %>"><span data-editor\><%= title %></label>\
            </div>\
          </div>\
        ');
      
      var includes = self.listModel.get('includes') || [];
      var _fields = this._options.schemaResult.fields;
      var _extra_scopes = [];
      var formSchema= {};
      
      // Create the schema fields:
      // - fields for the current resource are shown as normal
      // - if the scope is not the resource scope, create an extra_scope entry,
      // these items are indented.
      _.each(_.pairs(_fields), function(pair){

        var prop = pair[1];
        var key = prop['key'];
        var scope = prop['scope'];
        var fieldType = scope.split('.')[0]
        var fieldResource = scope.split('.')[1];
        
        if(fieldResource != self._options.resource.key){
          console.log('sub resource: ' + fieldResource + ',' + key);
          if(!_.has(_extra_scopes, scope)){
            _extra_scopes[scope] = [];
          }
          _extra_scopes[scope].push(key);
          formSchema[key] = {
            title: prop['title'],
            key: key,
            type: 'Checkbox',
            help: prop['description'],
            template: altSubCheckboxTemplate
          }
        }else{
          formSchema[key] = { 
            title: prop['title'], 
            key:  key, 
            type: 'Checkbox',
            help: prop['description'],
            template: altCheckboxTemplate };
        }
      });
      
      // default checkbox states
      var already_visible = {};
      var default_visible = {};

      // Build the form model
      var FormFields = Backbone.Model.extend({
        schema: formSchema
      });
      var formFields = new FormFields();
      _.each(_.pairs(_fields), function(pair){
        var key = pair[1]['key'];
        var prop = pair[1];
        var _visible = (_.has(prop, 'visibility') && 
            _.contains(prop['visibility'], 'l'));
        default_visible[key] = _visible;
        _visible = _visible || _.contains(includes, key);
        _visible = _visible && !_.contains(includes, '-'+key);
        already_visible[key] = _visible;
        formFields.set( key, _visible);
      });
      
      // Build the form template: build and append a field template for each field
      // "main scope" first - the fields for this resource
      var main_scope = 'fields.' + self._options.resource.key;
      var main_keys = _.filter(_.keys(_fields), function(key) {
        return _fields[key]['scope'] == main_scope;
      });
      main_keys = _.sortBy(main_keys, function(key){
        return _fields[key]['ordinal'];
      });
      _.each(main_keys, function(key){
        form_template.push( _.template(field_template, { 
          editorId: key+'-id', 
          title: _fields[key]['title'],
          name: key 
        }));
      });
      // second, any fields from other scopes/resources
      var _extra_scopes_shown = [];
      _.each(_.keys(_extra_scopes), function(scope){
        console.log('scope',scope)
        var fieldType = scope.split('.')[0]
        var fieldResource = scope.split('.')[1];
        var sub_resource;
        if(fieldType == 'datacolumn'){
          sub_resource = {
            title: fieldResource,
            help: 'Screen result data column field'
          }
        }else{
          sub_resource = appModel.getResource(fieldResource);
        }
        form_template.push(
          _.template(header_template, 
            {
              id: scope,
              name: sub_resource.title,
              help: sub_resource.description
            }));
        
        _.each(_extra_scopes[scope], function(sub_key){
          form_template.push( _.template(field_template, { name: sub_key }) );
          if(formFields.get(sub_key)){
            _extra_scopes_shown.push(scope);
          }
        });
        
      });
      form_template.push('</form>');
      
      console.log('extra_scopes_shown: ' + _extra_scopes_shown);
      var form = new Backbone.Form({
        model: formFields,
        template: _.template(form_template.join(''))
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
        },
        'click label.selection-group': function(e){
          //e.preventDefault();
          //e.stopPropagation();
          var debug_el = form.$el.find('input').each(function(){
            var key = $(this).attr('name');
            if(_.has(_fields, key) && _fields[key]['scope'] == e.target.id){
              $(this).closest('.form-group').toggle();
            }
          })
        },
        'click input.selection-group': function(e){
          e.preventDefault();
          //e.stopPropagation();
          var id = e.target.id;
          var scope = id.split('-')[0];
          console.log('id: ' + id + ', ' + scope + ', ' + e.target.checked );
          var debug_el = form.$el.find('input').each(function(){
            var key = $(this).attr('name');
            if(_.has(_fields, key) && _fields[key]['scope'] == scope ){
              $(this).closest('.form-group').show();
              form.setValue(key, e.target.checked );
            }
          });
        },
      };
      
      var _form_el = form.render().el;
      
      $(_form_el).find('.sub-resource-field').find('input').each(function(){
        var key = $(this).attr('name');
        var shown = !_.isUndefined(_.find(_extra_scopes_shown, function(scope){
          return scope == _fields[key]['scope']; }));
        $(this).closest('.form-group').toggle(shown);
      });

      _.each(_extra_scopes_shown, function(scope){
        if(_.every(_extra_scopes[scope], function(sub_key){
          return formFields.get(sub_key);
        })){
          $(_form_el).find(
              '#' + scope.split('.').join('\\.') + '-checkbox').prop("checked", true);
        }
      });
      
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
                      key,_fields[key],
                      self.collection.state.orderStack));
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
          
          // trigger an event to notify new header forms to self-display
          self.collection.trigger("MyServerSideFilter:search", 
            self.listModel.get('search'), self.collection);
        },
        view: _form_el,
        title: 'Select columns'  
      });
      
    },

    /**
     * Special function for screen result lists
     */
    show_mutual_positives: function(show_mutual_positives){
      var self = this;
      var _fields = this._options.schemaResult.fields;

      if (show_mutual_positives){
        
        _.each(_.pairs(_fields), function(pair){
          var key = pair[1]['key'];
          var prop = pair[1];
          var fieldType = prop['scope'].split('.')[0]
          if(fieldType == 'datacolumn'){
            var currentColumn = self.grid.columns.findWhere(
              { name: key });
            if(! currentColumn){
              self.grid.insertColumn(
                Iccbl.createBackgridColumn(
                    key,prop,
                    self.collection.state.orderStack));
            }
          }
        });
        
        var searchHash = _.clone(self.listModel.get('search'));
        searchHash['show_mutual_positives'] = 'true';
        self.listModel.set('search',searchHash);
      }else{
        _.each(_.pairs(_fields), function(pair){
          var key = pair[1]['key'];
          var prop = pair[1];
          var fieldType = prop['scope'].split('.')[0]
          if(fieldType == 'datacolumn'){
            var currentColumn = self.grid.columns.findWhere(
              { name: key });
            if(currentColumn){
              self.grid.removeColumn(currentColumn);
            }
          }
        });
        var searchHash = _.clone(view.listModel.get('search'));
        if(_.has(searchHash,'show_mutual_positives')){
          delete searchHash['show_mutual_positives'];
          view.listModel.set('search',searchHash);
        }        
      }
    },
    
  });
  

  return ListView;
});
