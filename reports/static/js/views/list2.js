define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_selector',
  'templates/rows-per-page.html',
  'templates/list2.html',
  'templates/modal_ok_cancel.html'
], function(
      $, _, Backbone, Backgrid,  
      Iccbl, appModel, genericSelector,
      rowsPerPageTemplate, listTemplate, modalTemplate) {

  
  var ListView = Backbone.View.extend({

    LIST_ROUTE_ORDER: appModel.LIST_ARGS,
    SEARCH_DELIMITER: appModel.SEARCH_DELIMITER,
    
    events: {
      'click button#select_columns': 'select_columns',
      'click button#download_link': 'download',
      'click button#clear_sorts': 'clear_sorts',
      'click button#clear_searches': 'clear_searches'
    },

    initialize : function(args) {
      
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

      // Presets
      
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
      
      // Update with the uriStack values
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
              var searchHash = self.searchHash = {};
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
            searchHash,
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
        if (!previousStack){
          // Note: replace: true - to suppress router history:
          // at this point, reportState is modifying the URL to show rpp, pagesSize, etc.
          appModel.set('routing_options', {replace: true});  
        }
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
          var fieldDef = schemaResult.fields[column];
          if( _.has(fieldDef,'backgrid_cell_options')) {
            var backgrid_cell_options = fieldDef['backgrid_cell_options'];
            var _route = Iccbl.formatString(backgrid_cell_options, model);
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
                  }
                }).fail(function(){ appModel.jqXHRfail.apply(this,arguments); });      
              }
            },
            view: _form_el,
            title: 'Confirm deletion of "' + title + '"'  
        });
      });
      
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
      	  goBackFirstOnSort: false, // Default is true
      	  collection: self.collection      	  
      	});            
        this.objects_to_destroy.push(paginator);
      }
      
      self.listenTo(self.listModel, 'change:search', function(){
        // TODO: this listener should be set in the collection initializer
        var searchHash = _.clone(self.listModel.get('search'));
        self.collection.setSearch(searchHash);
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
        className: "backgrid col-sm-12 table-striped table-condensed table-hover"
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
      self.$("#table-div").append(finalGrid.el);

      if(self.collection instanceof Backbone.PageableCollection){
        self.$("#paginator-div").append(self.paginator.render().$el);
        self.$("#rppselector").html(
            self.rppSelectorInstance.render().$el);
      }
      if(!_.isUndefined(self.extraSelectorInstance)){
        self.$("#extraselector").html(
            self.extraSelectorInstance.render().$el);
      }

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
      if(!_.isEmpty(orderStack)){
        self.collection.setSorting();
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
        self.$('#pagination_message').html(msg);
      });
      
      if ( !fetched ) {
        var fetchOptions = { reset: false };
        self.collection.fetch(
          fetchOptions
        ).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); }
        ).done(function(){ });      
      }
      
      // Note: replace: true - to suppress router history:
      // at this point, reportState is modifying the URL to show rpp, pagesSize, etc.
      appModel.set('routing_options', {replace: true});  
      this.reportState();
      return this;
    },
    
    clear_searches: function(){
      this.collection.trigger("Iccbl:clearSearches");
      // FIXME: any call to getFirstPage results in a fetch, disabled for now
      //this.collection.getFirstPage({reset: true, fetch: true});
      //this.collection.getFirstPage();
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
        "<div class='btn btn-default btn-sm ' id='select-all' >Select All</div>",
        "<div class='btn btn-default btn-sm ' id='clear-all' >Clear All</div>",
        "<button class='btn btn-default btn-sm ' id='modal-cancel'>Cancel</button>",
        "<button class='btn btn-default btn-sm ' id='modal-ok'>Ok</button>"
      ];
      var field_template = '<div data-fields="<%= name %>" ></div>';
      var optgroupSelectionTemplate = [
        '<div class="form-group" >',
        '<input class="selection-group" type="checkbox" id="<%= id %>-checkbox"> </input>',
        '<label class="selection-group .h5 " id="<%= id %>" title="<%= help %>" ><%= name %> columns</label>',
        '</div>'
        ].join('');
      var fieldCheckboxTemplate =  _.template('\
          <div class="form-group" style="margin-bottom: 0px;" > \
            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
              <label title="<%= help %>" for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>\
            </div>\
          </div>\
        ');      
      var optgroupFieldCheckboxTemplate =  _.template('\
          <div class="form-group sub-resource-field" style="margin-bottom: 0px;" > \
            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
            <label for="<%= editorId %>" > - </label>\
              <label for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>\
            </div>\
          </div>\
        ');
      
      var includes = self.listModel.get('includes') || [];
      var _fields = this._options.schemaResult.fields;
      
      // Create the (two level) tree of fields:
      // - fields for the current resource are shown as normal
      // - if the field information specifies an optgroup
      // - or if the field scope is not the current resource scope, create an optGroup entry,
      // these items are indented and nested in an optgroup form-group.
      var _optGroups = {};
      var already_visible = {};
      var default_visible = {};
      var _optgroups_shown = [];
      var defaultScope = 'fields.' + self._options.resource.key;
      var orderedKeys = _.sortBy(_.keys(_fields), function(key){
        return _fields[key]['ordinal'];
      });
      _optGroups[defaultScope] = {
        title: defaultScope,
        keys: []
      }
      
      _.each(orderedKeys, function(key){

        var prop = _fields[key];
        var optGroup = defaultScope;
        var title, fieldType, fieldResource, subResource;
        
        if (prop.display_options && prop.display_options.optgroup ){
          optGroup = prop.display_options.optgroup;
          title = optGroup.charAt(0).toUpperCase() + optGroup.slice(1);
          if(!_.has(_optGroups, optGroup)){
            _optGroups[optGroup] = 
              {
                title: title,
                help: title,
                keys: []
              };
          }
        }
        else if (prop['scope'] != defaultScope ){
          optGroup = prop['scope'];
          fieldType = optGroup.split('.')[0]
          fieldResource = optGroup.split('.')[1];
          if(fieldType == 'datacolumn'){
            title = fieldResource.charAt(0).toUpperCase() + fieldResource.slice(1);
            if(!_.has(_optGroups, optGroup)){
              _optGroups[optGroup] = 
                {
                  title: title,
                  help: 'Screen result data column field',
                  keys: []
                }
            }
          }else{
            if(!_.has(_optGroups, optGroup)){
              subResource = appModel.getResource(fieldResource);
              _optGroups[optGroup] = 
                {
                  title: subResource.title,
                  help: subResource.description,
                  keys: []
                }
            }
          }
        }
        _optGroups[optGroup].keys.push(key);

        var _visible = (_.has(prop, 'visibility') && 
            _.contains(prop['visibility'], 'l'));
        default_visible[key] = _visible;
        _visible = _visible || _.contains(includes, key);
        _visible = _visible && !_.contains(includes, '-'+key);
        already_visible[key] = _visible;
        
        if(_visible){
          _optgroups_shown[optGroup] = true;
        }
      });

      // Build the form model
      // Build the form template; manually lay out the fields/optgroups
      var formSchema= {};
      var FormFields = Backbone.Model.extend({
        schema: formSchema
      });
      var formFields = new FormFields();
      
      // first, the fields not in an optGroup
      _.each(_optGroups[defaultScope].keys,function(key){
        formSchema[key] = formFieldSchema = { 
          title: _fields[key]['title'], 
          key:  key, 
          type: 'Checkbox',
          help: _fields[key]['description'],
          template: fieldCheckboxTemplate 
        };

        form_template.push( 
          _.template(field_template)({ 
              editorId: key+'-id', 
              title: _fields[key]['title'],
              name: key 
            })
        );
        if(already_visible[key]){
          formFields.set( key, true);
        }
      });
      _.each(_.keys(_optGroups),function(optGroup){
        if (optGroup == defaultScope ) return;
        
        form_template.push(
          _.template(optgroupSelectionTemplate)( 
            {
              id: optGroup,
              name: _optGroups[optGroup].title,
              help: _optGroups[optGroup].description
            })
        );
        _.each(_optGroups[optGroup].keys,function(key){
          formSchema[key] = formFieldSchema = { 
            title: _fields[key]['title'], 
            key:  key, 
            type: 'Checkbox',
            help: _fields[key]['description'],
            template: optgroupFieldCheckboxTemplate 
          };
        
          form_template.push( 
            _.template(field_template)({ name: key }) );
          if(already_visible[key]){
            _optgroups_shown.push(optGroup);
            formFields.set( key, true);
          }
        });
      });
        
      form_template.push('</form>');

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
        // click on the optgroup text expands
        'click label.selection-group': function(e){
          _.each(_optGroups[e.target.id].keys, function(key){
            form.$el.find('[name='+key+']').closest('.form-group').toggle();
          });
        },
        // click on the optgroup checkbox checks all sub boxes
        'click input.selection-group': function(e){
          var optGroup = e.target.id.split('-')[0];
          console.log('id: ' + e.target.id + ', ' + optGroup + ', ' + e.target.checked );
          _.each(_optGroups[optGroup].keys, function(key){
            form.$el.find('[name='+key+']').closest('.form-group').show();
            form.setValue(key, e.target.checked );
          });
        },
      };
      
      var _form_el = form.render().el;
      
      // if any fields in the optGroup are shown, toggle the optGroup
      _.each(_optgroups_shown,function(optGroup){
        _.each(_optGroups[optGroup].keys,function(key){
          form.$el.find('[name='+key+']').closest('.form-group').toggle();
        });
      });
      
      // if all fields in the optGroup are shown, check the optGroup
      _.each(_optgroups_shown, function(optGroup){
        if(_.every(_optGroups[optGroup].keys, function(sub_key){
          return formFields.get(sub_key);
        })){
          $(_form_el).find(
              '#' + optGroup + '-checkbox').prop("checked", true);
        }
      });
      
      appModel.showModal({
        cancelText: 'Cancel',
        okText: 'Ok',
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
    show_mutual_positives: function(screen_facility_id, show_mutual_positives){
      var self = this;
      var _fields = this._options.schemaResult.fields;

      if (show_mutual_positives){
        
        _.each(_.pairs(_fields), function(pair){
          var key = pair[1]['key'];
          var prop = pair[1];
          var fieldType = prop['scope'].split('.')[0];
          var field_screen_facility_id = _.result(prop,'screen_facility_id', '');
          
          if(fieldType == 'datacolumn' 
            && field_screen_facility_id != screen_facility_id ){
            prop['visibility'].push('l');
            var currentColumn = self.grid.columns.findWhere({ name: key });
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
        var searchHash = _.clone(self.listModel.get('search'));
        var orderStack = self.listModel.get('order') || [];
        _.each(_.pairs(_fields), function(pair){
          var key = pair[1]['key'];
          var prop = pair[1];
          var fieldType = prop['scope'].split('.')[0]
          var field_screen_facility_id = _.result(prop,'screen_facility_id', '');
          // Note: if filtering/ordering on one of the mutual positive columns, 
          // do not remove it here.
          if(fieldType == 'datacolumn'
              && field_screen_facility_id != screen_facility_id 
              && ! _.findKey(searchHash, function(val,hashkey){
                return hashkey.indexOf(key) > -1
              })
              && ! _.find(orderStack, function(orderkey){
                return orderkey.indexOf(key) > -1
              })
            ){
            console.log('remove: ', key);
            prop['visibility'] = _.without(prop['visibility'], 'l');
            var currentColumn = self.grid.columns.findWhere({ name: key });
            if(currentColumn){
              self.grid.removeColumn(currentColumn);
            }
          }
        });
        if(_.has(searchHash,'show_mutual_positives')){
          delete searchHash['show_mutual_positives'];
          self.listModel.set('search',searchHash);
        }        
      }
    },
    
  });
  

  return ListView;
});
