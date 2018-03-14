define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_selector',
  'utils/treeSelector',
  'templates/rows-per-page.html',
  'templates/list2.html',
  'templates/modal_ok_cancel.html'
], function(
      $, _, Backbone, Backgrid,  
      Iccbl, appModel, genericSelector, TreeSelector,
      rowsPerPageTemplate, listTemplate, modalTemplate) {

  var ListView = Backbone.View.extend({

    LIST_MODEL_ROUTE_KEYS: appModel.LIST_ARGS,
    SEARCH_DELIMITER: appModel.SEARCH_DELIMITER,
    
    events: {
      'click button#select_columns': 'select_columns',
      'click button#download_link': 'download',
      'click button#clear_sorts': 'clear_sorts',
      'click button#clear_searches': function(e){
        e.preventDefault();
        this.clear_searches();
      }
    },
    
    /**
     * appModel.LIST_ARGS values may be sent as initial values in the args hash
     */
    initialize : function(args) {
      
      var self = this;
      var _options = self._options = args;
      var resource = self.resource = args.resource;
      self._classname = 'List2 - ' + resource.key;
      var urlSuffix = self.urlSuffix = "";
      var uriStack = args.uriStack || [];

      if (!self._options.extraControls){
        self._options.extraControls = [];
      }

      var ListModel = Backbone.Model.extend({
        defaults: {
            rpp: 25,
            page: 1,
            order: [],
            search: {},
            includes: [],
            esearch: '',
            csearch: '' }
        });

      var listInitial = this.parseUriStack(uriStack, _options.resource.options);
      var listModel = this.listModel = new ListModel(listInitial);

      this.objects_to_destroy = _([]);

      var _state = {
        currentPage: parseInt(self.listModel.get('page')),
        pageSize: parseInt(self.listModel.get('rpp'))
      };

      if (!_.has(self._options, 'url')) {
        self._options.url = self._options.resource.apiUri + '/' + self.urlSuffix;
      }
      
      var collection;
      if( !_options.collection){  
        var Collection = Iccbl.MyCollection.extend({
          state: _state,  
          modelId: function(attrs) {
            return Iccbl.getIdFromIdAttribute( attrs, resource);
          },
          'url': self._options.url,
          listModel: listModel,
          extraIncludes: _.result(args,'extraIncludes')
        });
        
        collection = self.collection = new Collection();
        this.objects_to_destroy.push(collection);
      }else{
        collection = self.collection = _options.collection;
        collection.listModel = listModel;
        collection.state = _state;
        collection.url = self._options.url;
      }

      var $modifySearch = $([
        '<a class="btn btn-default btn-sm pull-down" ',
        'role="button" id="modify_search_button" ',
        'href="#">',
        'Modify Search</a>'
      ].join(''));
      $modifySearch.click(function(e){
        e.preventDefault();
        self.modifySearch();
      });
      
      function setComplexSearch(){
        var searchId = self.listModel.get(appModel.URI_PATH_COMPLEX_SEARCH);
        if (!_.isUndefined(searchId)){
          console.log('change: ' + appModel.URI_PATH_COMPLEX_SEARCH, searchId);
          collection.type = 'POST';
          var searchData = appModel.getSearch(searchId);
          var new_url = self._options.url + '/' + appModel.URI_PATH_COMPLEX_SEARCH + '/' + searchId;
          collection.url = new_url;
          collection.fetch = function(options){
            var options = options || {};
            options.data = _.extend({}, options.data);
            options.data[appModel.API_PARAM_SEARCH] = searchData;
            options.type = 'POST';
            return Iccbl.MyCollection.prototype.fetch.call(this, options);
          };
        } else {
          console.log('no change/unset: ' + appModel.URI_PATH_COMPLEX_SEARCH, null);
        }
      };
      function setEncodedSearch(){
        var searchData = self.listModel.get(appModel.URI_PATH_ENCODED_SEARCH);
        if (!_.isUndefined(searchData)){
          console.log('change: ' + appModel.URI_PATH_ENCODED_SEARCH, searchData);
          collection.url = self._options.url ;
          collection.fetch = function(options){
            var options = options || {};
            options.data = _.extend({}, options.data);
            options.data[appModel.API_PARAM_SEARCH] = searchData;
            options.type = 'GET';
            return Iccbl.MyCollection.prototype.fetch.call(this, options);
          };
        } else {
          console.log('no change/unset: ' + appModel.URI_PATH_ENCODED_SEARCH, null);
        }
      };
      self.listenTo(self.listModel, 'change:'+appModel.URI_PATH_COMPLEX_SEARCH,setComplexSearch);
      self.listenTo(self.listModel, 'change:'+appModel.URI_PATH_ENCODED_SEARCH,setEncodedSearch);
      
      if (!_.isEmpty(listModel.get(appModel.URI_PATH_COMPLEX_SEARCH))){
        setComplexSearch();
        self._options.extraControls.unshift($modifySearch);
      }else if (! _.isEmpty(listModel.get(appModel.URI_PATH_ENCODED_SEARCH))){
        setEncodedSearch();
        self._options.extraControls.unshift($modifySearch);
      }
      
      var columns;
      if(!_options.columns){
        columns = Iccbl.createBackgridColModel(
            this._options.resource.fields, 
            listModel.get('order'),
            listModel.get(appModel.URI_PATH_SEARCH),
            listModel.get('includes'));
      }else{
        columns = _options.columns;
      }
      
      this.listenTo(this.listModel, 'change', this.reportState );
      this.buildGrid( columns, self._options.resource );
    },
    
    /**
     * Get the current extra search data for the list view:
     * appModel.URI_PATH_COMPLEX_SEARCH
     * appModel.URI_PATH_ENCODED_SEARCH
     */
    getSearchData: function() {
      var self = this;
      var search_data;
      var searchId = self.listModel.get(appModel.URI_PATH_COMPLEX_SEARCH);
      if (!_.isUndefined(searchId)  && !_.isEmpty(searchId)){
        search_data = appModel.getSearch(searchId);
      }else if (! _.isEmpty(self.listModel.get(appModel.URI_PATH_ENCODED_SEARCH))){
        search_data = self.listModel.get(appModel.URI_PATH_ENCODED_SEARCH);
      }
      return search_data;
    },
    
    /** 
     * Show/modify the reagent search list view 
     **/
    modifySearch: function() {
      var self = this;
      var resource = self.resource;
      var search_data;
      var search_path;
      var searchId = self.listModel.get(appModel.URI_PATH_COMPLEX_SEARCH);
      if (!_.isUndefined(searchId)  && !_.isEmpty(searchId)){
        search_data = appModel.getSearch(searchId);
        search_path = appModel.URI_PATH_COMPLEX_SEARCH;
      }else if (! _.isEmpty(self.listModel.get(appModel.URI_PATH_ENCODED_SEARCH))){
        search_data = self.listModel.get(appModel.URI_PATH_ENCODED_SEARCH);
        search_path = appModel.URI_PATH_ENCODED_SEARCH;
      } else {
        throw "no search found in the listmodel: " + self.listModel.toJSON();
      }
      
      if (_.isEmpty(search_data)){
        throw "no search found in the listmodel: " + self.listModel.toJSON();
      }
      console.log('modify search:', resource.key, search_data);
      
      
      function parse(value, errors){
        var parsedData;
        console.log('parse', value);
        if(resource.key == 'well' || resource.key == 'reagent'){
          parsedData = Iccbl.parseRawWellSearch(value, errors);
        }else if (resource.key == 'compound_search'){
          parsedData = Iccbl.parseCompoundVendorIDSearch(value,errors);
        }else if (resource.key == 'librarycopyplate'){
          parsedData = Iccbl.parseRawPlateSearch(value,errors);
        }else{
          throw 'Unknown search resource: ' + resource.key;
        }
          
        if (_.isEmpty(parsedData)){
          errors.push('no values found for input');
        } else {
          console.log('parsedData', parsedData);
        }
        return parsedData;
      };
      
      function validateSearch(value,formValues){
        var errors = [];
        var parsedData = parse(value, errors);
        if (!_.isEmpty(errors)){
          return {
            type: 'searchVal',
            message: errors.join('; ')
          };
        }
      };
      function processSearch(text_to_search){
        var errors =[];
        var parsedData = parse(text_to_search, errors);
        if (!_.isEmpty(errors)){
          throw 'Unexpected error on submit: ' + errors.join(', ');
        }
        
        if (parsedData.length <= appModel.MAX_RAW_SEARCHES_IN_URL){
          var search_text = _.map(parsedData, function(parsedLine){
            if (_.has(parsedLine, 'combined')){
              return parsedLine.combined.join(' ');
            }else{
              return parsedLine;
            }
          }).join(appModel.UI_PARAM_RAW_SEARCH_LINE_ENCODE);
          if (search_path != appModel.URI_PATH_ENCODED_SEARCH){
            appModel.set('routing_options', {replace: false});  
          }
          self.listModel.unset(appModel.URI_PATH_COMPLEX_SEARCH, { silent: true});
          self.listModel.set(appModel.URI_PATH_ENCODED_SEARCH,search_text);
          
        }else{
          // Send complex search data as a POST
          // must change the route, and create a post
          // TODO: key the search data using the searchId: 
          // this allows for browser "back" in the session
          // will also need to listen to URIStack changes and grab the data
          // from the search ID
          var searchId = '' + ( new Date() ).getTime();
          appModel.setSearch(searchId,text_to_search);
          
          if (search_path != appModel.URI_PATH_COMPLEX_SEARCH){
            appModel.set('routing_options', {replace: false});  
          }
          self.listModel.unset(appModel.URI_PATH_ENCODED_SEARCH, { silent: true});
          self.listModel.set(appModel.URI_PATH_COMPLEX_SEARCH,searchId);
          
        }
        self.collection.getFirstPage({reset: true, fetch: true});
        
      };
      
      function searchToText(parsedData){
        search_text = _.map(parsedData, function(parsedLine){
          if (_.has(parsedLine, 'combined')){
            return parsedLine.combined.join(' ');
          }else{
            return parsedLine;
          }
        }).join('\n');
        return search_text;
      };
      var schema = {
        searchVal: {
          title: 'Search Data',
          key: 'searchVal',
          type: Backbone.Form.editors.TextArea.extend({
              initialize: function() {
                Backbone.Form.editors.TextArea.prototype.initialize.apply(this, arguments);
                this.$el.attr('rows',12);
              }
          }),
          editorClass: 'input-full',
          validators: ['required', validateSearch],
          template: appModel._field_template //altFieldTemplate
        }
      };

      var errors = [];
      var parsedData = parse(search_data, errors);
      if (!_.isEmpty(errors)){
        throw 'errors in parsed search:' + errors.join(', ');
      } else {
        search_data = searchToText(parsedData);
      }
      var FormFields = Backbone.Model.extend({
        schema: schema
      });
      var formFields = new FormFields({
        searchVal: search_data
      });
      var form = new Backbone.Form({
        model: formFields,
        template: appModel._form_template //_.template(form_template.join(''))
      });
      var _form_el = form.render().el;
      var title = 'Search for ' + resource.title;
      
      appModel.showModal({
        title: title,
        view: _form_el,
        ok: function() {
          var errors = form.commit();
          if(!_.isEmpty(errors)){
            console.log('form errors, abort submit: ' + JSON.stringify(errors));
            return false;
          }
          
          processSearch(form.getValue()['searchVal']);
        } 
      });
    },
    
    /**
     * Parse the uriStack and set the listModel state
     */
    parseUriStack: function(uriStack, initialOptions ){
      console.log('parseUriStack', uriStack, initialOptions);
      var self = this;
      var listInitial = initialOptions || {};
      listInitial = _.extend({},_.pick(listInitial,self.LIST_MODEL_ROUTE_KEYS));
      // Update with the uriStack values
      for (var i=0; i<uriStack.length; i++){
        var key = uriStack[i];
        i = i+1;
        if (i==uriStack.length) continue;
        var value = uriStack[i];
        if (!value || _.isEmpty(value) ){
          continue;
        }
        
        if(key == 'log'){
          self.urlSuffix = key + '/' + value;
          continue;
        }
        
        if(_.contains(self.LIST_MODEL_ROUTE_KEYS, key)){
          
          if(key === appModel.URI_PATH_SEARCH) {
            var searchHash = {};
            var searches = value.split(self.SEARCH_DELIMITER);
            _.each(searches, function(search){
              var parts = search.split('=');
              if (!parts || parts.length!=2) {
                var msg = 'invalid search parts: ' + search;
                console.log(msg);
                appModel.error(msg);
              } else if (_.isEmpty(parts[1])) {
              } else {
                searchHash[parts[0]] = decodeURIComponent(parts[1]);
              }
            });
            console.log('searchHash', searchHash);
            
            listInitial[key] = searchHash;
          } else if (key == appModel.URI_PATH_ENCODED_SEARCH) {
            listInitial[key] = value;
          } else if (key == appModel.URI_PATH_COMPLEX_SEARCH){
            var searchId = '' + value;
            if (_.isEmpty(searchId)){
              var msg = 'no complex search id found:' + uriStack.join('/');
              console.log(msg);
              appModel.error(msg);
            }
            if (!_.isNaN(parseInt(searchId))) {
              var search_data = appModel.getSearch(searchId);
              if(_.isUndefined(search_data) || _.isEmpty(search_data)){
                var msg = 'Browser state no longer contains searchId:'+searchId 
                  +'", search must be performed again';
                console.log(msg);
                appModel.error(msg);
              }else{
                listInitial[key] = searchId;
              }  
            }else{
              var msg = 'search ID is not a valid number: ' + searchId;
              console.log(msg);
              appModel.error(msg);
            }
          } else if (key === 'order') {
            listInitial[key] = value.split(',');        
          } else if (key === 'includes') {
            listInitial[key] = value.split(',');        
          }else {
            listInitial[key] = value;
          }
        }
      }
      console.log('parseUriStack, listInitial', listInitial);
      return listInitial;
    },
    
    /**
     * Get the current URL for the collection:
     * - for downloads
     * - create a display title
     */
    getCollectionUrl: function(limit) {
      var self = this;
      var url = self.collection.url;
      var urlparams = '';
      var search_title_val = '';
      
      function printFieldAndOperator(k){
        if(! k || _.isEmpty(k)) return k;
        var parts = k.split('__');
        var name = k;
        if(parts.length == 2){
          var field = self._options.resource.fields[parts[0]];
          if(field){
            name = field.title;
          }
          name += ' ' + parts[1];
        }
        return name;
      };
      
      _.each(self.LIST_MODEL_ROUTE_KEYS, function(routeKey){
        var routeEntry = self.listModel.get(routeKey);
        if ( ! _.isEmpty(routeEntry)) {
          if (routeKey === appModel.URI_PATH_SEARCH) {
            if(!_.isEmpty(urlparams)) urlparams += '&';
            urlparams += _.map(
              _.zip(_.keys(routeEntry),_.values(routeEntry)), 
              function(kv){
                // NOTE: 201705 support for the raw search data
                kv[1] = encodeURIComponent(kv[1]);
                return kv.join('=');
              }).join('&');
            routeEntry = _.extend({},routeEntry);
            // TODO: dc_ids used for reagent/screen_result views, should not be in search title
            delete routeEntry['dc_ids'];
            if(search_title_val !== '') search_title_val += '<br>AND ';
            search_title_val += _.map(
              _.zip(_.keys(routeEntry),_.values(routeEntry)), 
              function(kv){
                return [printFieldAndOperator(kv[0]),kv[1]].join('=');
              }).join(' AND ');
          } else if (routeKey == appModel.URI_PATH_ENCODED_SEARCH) {
            if(!_.isEmpty(urlparams)) urlparams += '&';
            urlparams += appModel.API_PARAM_SEARCH + '=' + routeEntry;
            if(search_title_val !== '') search_title_val += '<br>AND ';
            search_title_val += '[' + routeEntry + ']';
          } else if (routeKey == appModel.URI_PATH_COMPLEX_SEARCH){
            // nop
          }else if (routeKey === 'order') {
            if(!_.isEmpty(urlparams)) urlparams += '&';
            urlparams += 'order_by=' + routeEntry.join('&order_by=');
          }else if (routeKey === 'includes') {
            if(!_.isEmpty(urlparams)) urlparams += '&';
            urlparams += 'includes=' + routeEntry.join('&includes=');
          }
        }
      });
      
//      if (_.has(self._options, appModel.API_PARAM_ENCODED_SEARCH)){
//        if(!_.isEmpty(urlparams)) urlparams += '&';
//        urlparams += appModel.API_PARAM_SEARCH + '=' 
//          + self._options[appModel.API_PARAM_ENCODED_SEARCH];
//        if(search_title_val !== '') search_title_val += '<br>AND ';
//        search_title_val += '[' + self._options[appModel.API_PARAM_ENCODED_SEARCH] + ']';
//      }
      
      
      if(_.isUndefined(limit)){
        var limit = self.collection.state.totalRecords;
        if(!_.isNumber(limit)){
          limit = 0;
        }
      }
      url += '?limit=' + limit;
      if(!_.isEmpty(urlparams)){
        url += '&' + urlparams;
      }
      console.log('collection url: ' + url)
      
      if(_.has(self._options, appModel.API_PARAM_NESTED_SEARCH)){
        console.log('nested_search_data found on the self._options: ', 
          self._options[appModel.API_PARAM_NESTED_SEARCH]);
        var or_clauses = '';
        var search_list = self._options[appModel.API_PARAM_NESTED_SEARCH];
        _.each(search_list, function(hash){
          var temp = '';
          _.each(_.keys(hash), function(k){
            var v = hash[k];
            if (temp !== '' ) temp += ' AND ';
            temp += printFieldAndOperator(k) + " (" + v + ") ";
          });
          if(or_clauses != '') or_clauses += "<br>OR ";
          or_clauses += temp
        });
        if (!_.isEmpty(or_clauses)){
          if(search_title_val !== '') search_title_val += '<br>AND ';
          search_title_val += '[ ' + or_clauses + ' ]';
        }
      }            

      if(search_title_val !== ''){
        self.trigger('update_title', 'Search: ' + search_title_val);
      } else {
        self.trigger('update_title', '');
      }
      
      return url;
    },
    
    /**
     * Report the listmodel "state" to the application using the uriStack
     */
    reportState: function(args) {
      console.log('report state', args);
      var self = this;
      var newStack = [];
      var previousStack = self.currentStack;
      
      // Note: this is repeated in the specific change listeners
      $('#clear_sorts').toggle(!_.isEmpty(self.listModel.get('order')));
      $('#clear_searches').toggle(!_.isEmpty(self.listModel.get(appModel.URI_PATH_SEARCH)));
     
      // If a suffix was consumed, then put it back
      if(self.urlSuffix != ""){
        newStack = self.urlSuffix.split('/');
      }
      _.each(self.LIST_MODEL_ROUTE_KEYS, function(routeKey){
        var routeEntry = self.listModel.get(routeKey);
        if ( (!_.isObject(routeEntry) && routeEntry ) || 
             ( _.isObject(routeEntry) && !_.isEmpty(routeEntry))) {
          newStack.push(routeKey);
          if (routeKey === appModel.URI_PATH_SEARCH) {
            newStack.push(_.map(
              _.zip(_.keys(routeEntry),_.values(routeEntry)), 
              function(kv){
                // NOTE: 201705 support for the raw search data
                // kv[1] = encodeURIComponent(kv[1]);
                return kv.join('=');
              }).join(self.SEARCH_DELIMITER));
          } else if (routeKey == appModel.URI_PATH_ENCODED_SEARCH) {
            newStack.push(routeEntry);
          } else if (routeKey == appModel.URI_PATH_COMPLEX_SEARCH){
            newStack.push(routeEntry);
          } else if (routeKey === 'order') {
            newStack.push(routeEntry);
          }else if (routeKey === 'includes') {
            newStack.push(routeEntry.join(','));
          } else {
            newStack.push(routeEntry);
          }
        }
      });
      self.currentStack = newStack;
      console.log('newStack: ' + JSON.stringify(newStack));
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
      
      self.listModel.set({ 
        'rpp': state.pageSize, 
        'page': currentPage,
      }, {silent: true});
      
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
          editorClass: 'input-full',
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
        var title = Iccbl.getTitleFromTitleAttribute(model,self._options.resource);
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
                model.set('id', Iccbl.getIdFromIdAttribute(model,self._options.resource));
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
      this.listenTo(this.listModel, 'change:rpp', function(){
          rppModel.set({ selection: String(self.listModel.get('rpp')) });
      });
      this.listenTo(rppModel, 'change', function() {
          var rpp = parseInt(rppModel.get('selection'));
          console.log('===--- rppModel change: ' + rpp );
          self.listModel.set('rpp', String(rpp));
          self.listModel.set('page',1);
          // set this because of how checkstate is triggered
          self.collection.state.currentPage = 1; 
          appModel.set('routing_options', {replace: false});  
          self.collection.setPageSize(rpp, { first: true });
      });
      
      this.listenTo(this.listModel, 'change:page', function(){
        var page = self.listModel.get('page');
        self.collection.getPage(page, {reset: true });
        self.reportState();
      });
      this.listenTo(this.listModel, 'change:'+ appModel.URI_PATH_SEARCH, function(){
        // TODO: this listener should be set in the initializer
        var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
        // Note: this is repeated in reportState
        $('#clear_searches').toggle(!_.isEmpty(searchHash));
        if (self.listModel.get('page') !== 1){
          self.collection.setSearch(searchHash, { fetch: false });
          self.listModel.set('page', 1, {silent: true});
          self.collection.getPage(1);
        } else {
          self.collection.setSearch(searchHash);
        }
      });
      
      this.listenTo(this.listModel, 'change:order', function(){
        // Note: this is repeated in reportState
        $('#clear_sorts').toggle(!_.isEmpty(self.listModel.get('order')));
        
      });
      
      this.listenTo(this.listModel, 'change:includes', 
        function(model, changed, options){
          
          var reset = true;
          if (options && options.reset === false){
            reset = options.reset;
          }
          
          var fields = self._options.resource.fields;
          var toAdd = [];
          var toRemove = [];
          var previous = self.listModel.previous('includes');
          var current = self.listModel.get('includes');
          console.log('previous includes:', previous, 'current', current);
          
          _.each(_.difference(previous,current), function(premoved){
            var negate = false;
            if(premoved.charAt(0)=='-'){
              premoved = premoved.slice(1);
              negate = true;
            }
            if (!_.has(fields, premoved)){
              console.log('unknown column', premoved);
              return;
            }
            var defaultVisible = _.contains(fields[premoved]['visibility'], 'l');
            if (defaultVisible==true){
              if (negate==true){
                toAdd.push(premoved);
              }else{
                //pass
              }
            } else {
              if (negate != true){
                toRemove.push(premoved);
              } else {
                //pass
              }
            }
          });
          _.each(_.difference(current,previous), function(cadded){
            var negate = false;
            if(cadded.charAt(0)=='-'){
              cadded = cadded.slice(1);
              negate = true;
            }
            if (!_.has(fields, cadded)){
              console.log('unknown column', cadded);
              return;
            }

            var defaultVisible = _.contains(fields[cadded]['visibility'], 'l');
            if (defaultVisible==true){
              if (negate==true){
                toRemove.push(cadded);
              }else{
                //pass
              }
            } else {
              if (negate != true){
                toAdd.push(cadded);
              } else {
                //pass
              }
            }
          });
          
          console.log('toRemove:', toRemove, 'toAdd', toAdd);
          
          _.each(toAdd, function(key){
            var field = fields[key];
            console.log('add column', key, field['ordinal']);
            var column = self.grid.columns.findWhere({ name: key });
            if (!column){
              // find out where it goes
              var ordinal = field['ordinal'];
              var index = 0;
              self.grid.columns.find(function(column){
                var colKey = column.get('name');
                var colField = fields[colKey];
                var colOrdinal = colField['ordinal'];
                if(colOrdinal>ordinal){
                  console.log('add col', key, ordinal, 'before col', colKey,colOrdinal)
                  return true;
                }
                index += 1;
              });
              
              self.grid.insertColumn(
                Iccbl.createBackgridColumn(
                    key,field,
                    self.listModel.get('order')),
                    { at: index});
            } else {
              console.log('column already included', key)
            }
          });
          _.each(toRemove, function(key){
            var column = self.grid.columns.findWhere({ name: key });
            if (!column){
              console.log('column already not present', key)
            } else {
              self.grid.removeColumn(column);
            }
          });
          
          if (!_.isEmpty(toRemove)){
            self.clear_searches(toRemove);
          }
          
          if(reset){
            self.collection.fetch();
          }
          
          // trigger an event to notify new header forms to self-display
          self.collection.trigger("MyServerSideFilter:search", 
            self.listModel.get(appModel.URI_PATH_SEARCH), self.collection);
      });

      if(self.collection instanceof Backbone.PageableCollection){
        var paginator = self.paginator = new Backgrid.Extension.Paginator({
      	  goBackFirstOnSort: false, // Default is true
      	  collection: self.collection      	  
      	});            
        this.objects_to_destroy.push(paginator);
      }
      
      // Extraselector
      if( _.has(schemaResult, 'extraSelectorOptions')){
        var searchHash = self.listModel.get(appModel.URI_PATH_SEARCH);
        console.log('extraselector init: searchTerms: ' + JSON.stringify(searchHash));

        var extraSelectorModel = new Backbone.Model({ selection: '' });
        var extraSelectorKey = schemaResult.extraSelectorOptions.searchColumn;
        
        if ( !_.isEmpty(searchHash)){
          _.each(_.keys(searchHash), function(key){
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
          var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
          _.each(_.keys(searchHash), function(key){
            if( key === extraSelectorKey || key  === extraSelectorKey+ '__exact'){
                extraSelectorModel.set({ selection: searchHash[key] });
            }
          });
        });
        this.listenTo(extraSelectorModel, 'change', function() {
          var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
          var val = extraSelectorModel.get('selection');
          var value =  _.isUndefined(val.value) ? val: val.value ;
          delete searchHash[extraSelectorKey + '__exact']
          delete searchHash[extraSelectorKey + '__ne']
          delete searchHash[extraSelectorKey]
          
          if(!_.isEmpty(value) && !_.isEmpty(value.trim())){
            var field = self._options.resource.fields[extraSelectorKey];
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
          self.collection.setSearch(searchHash);
        });
      }
      
      var body = Iccbl.MultiSortBody;
      if (self._options.body) {
        body = body.extend(self._options.body);
      }
      
      var grid = this.grid = new Backgrid.Grid({
        columns: columns,
        body: body,
        row: self._options.row,
        collection: self.collection,
        className: "backgrid left-align col-sm-12 table-striped table-condensed table-hover"
      });
      this.objects_to_destroy.push(grid);
      //var gridHeader = this.gridHeader = new Backgrid.Grid({
      //  columns: columns,
      //  body: body,
      //  row: self._options.row,
      //  collection: self.collection,
      //  className: "backgrid col-sm-12 table-striped table-condensed table-hover"
      //});
      //this.objects_to_destroy.push(gridHeader);

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
      // TODO: use the backbone template and serialize attributes
      var title = _.result(self._options, 'title', '');
      this.$el.html(_.template(listTemplate)(
          { title: _.result(self._options, 'title', '')})
      );
    },
      
//      Aborted attempt to follow:
//      https://github.com/santomegonzalo/backgrid.fixedheader
//    _scrollingContainer: function(evt) {
//      var self = this;
//      var $element = $(evt.target);
//      var top = $element.scrollTop();
//      var left = $element.scrollLeft();
//      
//      self.container.find('.backgrid-fixed-header').css('margin-left', -Math.abs(left));
//      
//      if( self.container.find('.backgrid-fixed-content').scrollTop() 
//          >= self.container.find('.backgrid-fixed-content').find('.backgrid').height() 
//            - self.container.find('.backgrid-fixed-content').height() - 100) {
//        $.event.trigger('BackgridContainer:scrollEnd');
//      }
//    },
//    _resize: function(evt) {
//      var self = this;
//      var width = self.container.find('.backgrid-fixed-content').find('.backgrid').width();
//      self.container.find('.backgrid-fixed-header').width(width);
//    },

    afterRender: function(){
      var self = this;
      var fetched = false;
      
      self.listenTo(self.collection, "add", self.checkState);
      self.listenTo(self.collection, "remove", self.checkState);
      self.listenTo(self.collection, "reset", self.checkState);
      self.listenTo(self.collection, "sort", self.checkState);
      self.listenTo(self.collection, "pageable:state:change", function(state){
        if (state && state.currentPage){
          self.listModel.set('page', state.currentPage, {silent: true});
          self.reportState();
        }
      });
      var finalGrid = self.finalGrid = this.grid.render();
      self.objects_to_destroy.push(finalGrid);
      self.$("#table-div").append(finalGrid.el);
      
      //Aborted attempt to follow:
      //https://github.com/santomegonzalo/backgrid.fixedheader
      //var container = self.container = self.$("#table-div");
      //
      //$(window).bind('resize', self._resize);
      //container.append('<div class="backgrid-fixed-header"></div>');
      //container.append('<div class="backgrid-fixed-content"></div>');      
      //container.find('.backgrid-fixed-content').bind('scroll', self._scrollingContainer);      
      //
      //this.grid.collection.bind("add", function(){
      //  self._resize();
      //});      
      //
      //self.grid.render();
      //self.gridHeader.render();
      //
      //container.find('.backgrid-fixed-header').append(self.gridHeader.$el);
      //container.find('.backgrid-fixed-content').append(self.grid.$el);
      //
      ////      self.gridHeader.$el.find('tbody').hide(); //addClass('remove_body');
      ////      self.grid.$el.find('thead').hide(); // addClass('remove_header');
      
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
        self.$('#extra_controls').append(
          '<div id="extra_controls_div" class="panel"></div>');
        _.each(self._options.extraControls, function(control){
          // Adjust the checkbox types, so that the first also has a margin
          // otherwise, wrapped checkboxes are offset
          control.has('input[type="checkbox"]').css('margin-left','10px');
          console.log('append extra control: ', control);
          self.$('#extra_controls_div').append(control);
        });
      }

      this.delegateEvents();
      
      var searchHash = self.listModel.get(appModel.URI_PATH_SEARCH);
      if(!_.isEmpty(searchHash)){
        self.collection.setSearch(searchHash, {reset: false});
        fetched = true;
      }

      var orderStack = self.listModel.get('order') || [];
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
        //{ reset: false } (default) - uses set to (intelligently) merge the fetched 
        //models ("add" events are fired),
        //{reset: true}, in which case the collection will be (efficiently) reset 
        //(no "add" events will be fired)
        self.collection.fetch({ reset: false }
        ).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); }
        ).done(function(){ });      
      }
      
      // Note: replace: true - to suppress router history:
      // at this point, reportState is modifying the URL to show rpp, pagesSize, etc.
      appModel.set('routing_options', {replace: true});  
      this.reportState();
      return this;
    },
    
    clear_searches: function(fields_to_clear){
      var self = this;
      // notify the column headers
      if (fields_to_clear){
        this.collection.trigger(
          "Iccbl:clearSearches", {fields_to_clear: fields_to_clear});
      } else {
        this.collection.trigger("Iccbl:clearSearches");
      }

      var searchHash = _.clone(this.listModel.get(appModel.URI_PATH_SEARCH));
      console.log('clear hash', searchHash, fields_to_clear);
      var fields = this._options.resource.fields;
      if (fields_to_clear){
        fields = _.pick(fields,fields_to_clear);
      }
      _.each(_.keys(searchHash), function(key){
        var originalKey = key;
        if (key.indexOf('__')>0){
          key = key.slice(0,key.indexOf('__'));
        }
        if (key.charAt(0)=='-'){
          key = key.slice(1);
        }
        if (_.has(fields, key)){
          console.log('clearing key', originalKey);
          delete searchHash[originalKey];
        } else {
          console.log('Not clearing search term:', originalKey);
        }
      });
      console.log('cleared hash', searchHash);
      this.listModel.set(appModel.URI_PATH_SEARCH, searchHash);
      // FIXME: any call to getFirstPage results in a fetch, disabled for now
      //this.collection.getFirstPage({reset: true, fetch: true});
      //this.collection.getFirstPage();
    },
    
    clear_sorts: function(){
      this.collection.trigger("Iccbl:clearSorts");
      this.listModel.unset('order');
      this.collection.getFirstPage({reset: true, fetch: true});
    },
    
    download: function(e){
      var self = this;
      var limitForDownload = 0;

      e.preventDefault();
      
      // find any extra search data
      if(_.has(self._options, appModel.API_PARAM_SEARCH)){
        // endcode for the post_data arg; post data elements are sent 
        // as urlencoded values via a POST form for simplicity
        var data = {};
        data[appModel.API_PARAM_SEARCH] = 
                self._options[appModel.API_PARAM_SEARCH];
        appModel.download(
            self.getCollectionUrl(limitForDownload), 
            self._options.resource, data); 
      }else if(_.has(self._options, appModel.API_PARAM_NESTED_SEARCH)){
        // endcode for the post_data arg; post data elements are sent 
        // as urlencoded values via a POST form for simplicity
        var data = {};
        data[appModel.API_PARAM_NESTED_SEARCH] = 
                self._options[appModel.API_PARAM_NESTED_SEARCH];
        appModel.download(
            self.getCollectionUrl(limitForDownload), 
            self._options.resource,data); 
      }else{
        appModel.download(
            self.getCollectionUrl(limitForDownload), 
            self._options.resource);
      }
    },
    
    /**
     * Select columns to display
     * - using TreeSelector widget
     */
    select_columns: function(event){
      console.log('select_columns...');
      var self = this;
      
      // 1. build a columns collection
      var includes = self.listModel.get('includes') || [];
      var _fields = this._options.resource.fields;
      
      _fields = _.omit(_fields, function(value, key, object) {
        return _.contains(value.visibility, 'api');
      });
      
      var screenModel = _.result(self._options,'screen');
      if (screenModel){
        // Special case "datacolumn" screen result fields: 
        // - only select data columns for this screen; other screen data
        // data columns are selected through the special "other screen columns"
        // selector.
        _fields = _.omit(_fields, function(value, key, object) {
          if (value.is_datacolumn === true){
            if (screenModel.get('facility_id') 
                == value.screen_facility_id ){
              return false;
            } else {
              return true;
            }
          } else {
            return false;
          }
        });
      }      
      var ColumnCollection = Backbone.Collection.extend({
        modelId: function(attrs){
          var id = Iccbl.getIdKeys(attrs,{ id_attribute: ['scope','key']}).join('-');
          id = id.replace('.','_');
          return id;
        }
      });
      
      var columnCollection = new ColumnCollection();
      _.each(_fields, function(field){
        var key = field['key'];
        var model = new Backbone.Model(field);
        if (model.has('scope')){
          var scope = model.get('scope');
          var fieldType = scope.split('.')[0]
          var fieldResource = scope.split('.')[1];
          fieldResource = fieldResource.charAt(0).toUpperCase() + fieldResource.slice(1);
          
          model.set('resource', fieldResource);
        }

        var _visible = (_.has(field, 'visibility') && 
            _.contains(field['visibility'], 'l'));
        _visible = _visible || _.contains(includes, key);
        
        if (_.contains(includes, '-'+key)){
          _visible = false;
        }
        
        if (_visible == true){
          model.set('checked', true);
        }
        columnCollection.add(model);
      });
      console.log('columnCollection', columnCollection);

      function showColumns(collection){
        var new_includes = [];
        collection.each(function(model){
          var key = model.get('key');
          var isDefaultVisible = _.contains(model.get('visibility'),'l');
          if (model.get('checked') == true){
            if (isDefaultVisible !== true){
              new_includes.unshift(key);
            }
          } else {
            if (isDefaultVisible === true){
              new_includes.unshift('-'+key);
            }
          }
        });
        console.log('new includes', new_includes);
        self.listModel.set({'includes': new_includes });
      };
      
      // 2. initialize the tree selector
       var dcView = new TreeSelector({
        collection: columnCollection,
        treeAttributes: ['resource', 'title'],
        extraControls: [],
        startExpanded: true
       });
       Backbone.Layout.setupView(dcView);
  
        var el = dcView.render().el;
        var dialog = appModel.showModal({
          buttons_on_top: true,
          css: { 
              display: 'table',
              height: '500px',
              width: '80%'
            },
          css_modal_content: {
            overflow: 'hidden'
          },
          ok: function(){
            showColumns(columnCollection);
          },
          view: el,
          title: 'Select Columns to display'
        });
    },
    
//    /** Build the select columns dialog **/
//    /** FIXME: 20170809 - replace with the TreeSelector **/
//    select_columns1: function(event){
//
//      console.log('select_columns...');
//      
//      var self = this;
//      var form_template = [
//        "<form  class='form-horizontal container' >",
//        "<div class='btn btn-default btn-sm ' id='select-all' >Select All</div>",
//        "<div class='btn btn-default btn-sm ' id='clear-all' >Clear All</div>",
//        "<button class='btn btn-default btn-sm ' id='modal-cancel'>Cancel</button>",
//        "<button class='btn btn-default btn-sm ' id='modal-ok'>Ok</button>"
//      ];
//      var field_template = '<div data-fields="<%= name %>" ></div>';
//      var optgroupSelectionTemplate = [
//        '<div class="form-group" >',
//        '<input class="selection-group" type="checkbox" id="<%= id %>-checkbox"> </input>',
//        '<label class="selection-group .h5 " id="<%= id %>" title="<%= help %>" ><%= name %> columns</label>',
//        '</div>'
//        ].join('');
//      var fieldCheckboxTemplate =  _.template('\
//          <div class="form-group" style="margin-bottom: 0px;" > \
//            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
//              <label title="<%= help %>" for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>\
//            </div>\
//          </div>\
//        ');      
//      var optgroupFieldCheckboxTemplate =  _.template('\
//          <div class="form-group sub-resource-field" style="margin-bottom: 0px;" > \
//            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
//            <label for="<%= editorId %>" > - </label>\
//              <label title="<%= help %>" for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>\
//            </div>\
//          </div>\
//        ');
//      var optgroupFieldCheckboxTemplate1 =  _.template('\
//          <div class="form-group sub-resource-field" style="margin-bottom: 0px;" > \
//            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
//              <label title="<%= help %>" for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>\
//            </div>\
//          </div>\
//        ');
//      
//      var includes = self.listModel.get('includes') || [];
//      var _fields = this._options.resource.fields;
//      
//      _fields = _.omit(_fields, function(value, key, object) {
//        return _.contains(value.visibility, 'api');
//      });
//      
//      // Create the (two level) tree of fields:
//      // - fields for the current resource are shown as normal
//      // - if the field information specifies an optgroup
//      // - or if the field scope is not the current resource scope, create an optGroup entry,
//      // these items are indented and nested in an optgroup form-group.
//      var _optGroups = {};
//      var already_visible = {};
//      var default_visible = {};
//      var _optgroups_shown = [];
//      var defaultScope = 'fields.' + self._options.resource.key;
//      console.log('defaultScope: ' , defaultScope);
//      var orderedKeys = _.sortBy(_.keys(_fields), function(key){
//        return _fields[key]['ordinal'];
//      });
//      _optGroups[defaultScope] = {
//        title: defaultScope,
//        keys: []
//      }
//      
//      _.each(orderedKeys, function(key){
//
//        var prop = _fields[key];
//        var optGroup = defaultScope;
//        var title, fieldType, fieldResource, subResource;
//        
//        if (prop.display_options && prop.display_options.optgroup ){
//          optGroup = prop.display_options.optgroup;
//          title = optGroup.charAt(0).toUpperCase() + optGroup.slice(1);
//          if(!_.has(_optGroups, optGroup)){
//            _optGroups[optGroup] = 
//              {
//                title: title,
//                help: title,
//                keys: []
//              };
//          }
//        }
//        else if (prop['scope'] == 'otherscreen.datacolumns'){
//          optGroup = 'otherscreen';
//          if(!_.has(_optGroups, optGroup)){
//            _optGroups[optGroup] = 
//              {
//                title: 'Other Screen',
//                help: 'Other screen data column fields',
//                keys: []
//              }
//          }
//          
//        }
//        else if (prop['scope'] != defaultScope ){
//          optGroup = prop['scope'];
//          fieldType = optGroup.split('.')[0]
//          fieldResource = optGroup.split('.')[1];
//          if(fieldType == 'datacolumn'){
//            title = fieldResource.charAt(0).toUpperCase() + fieldResource.slice(1);
//            if(!_.has(_optGroups, optGroup)){
//              _optGroups[optGroup] = 
//                {
//                  title: title,
//                  help: 'Screen result data column field',
//                  keys: []
//                }
//            }
//          }else{
//            if(!_.has(_optGroups, optGroup)){
//              subResource = appModel.getResource(fieldResource);
//              _optGroups[optGroup] = 
//                {
//                  title: subResource.title,
//                  help: subResource.description,
//                  keys: []
//                }
//            }
//          }
//        }
//        
//        _optGroups[optGroup].keys.push(key);
//
//        var _visible = (_.has(prop, 'visibility') && 
//            _.contains(prop['visibility'], 'l'));
//        default_visible[key] = _visible;
//        _visible = _visible || _.contains(includes, key);
//        _visible = _visible && !_.contains(includes, '-'+key);
//        already_visible[key] = _visible;
//
//        console.log('col: ', key, optGroup, _visible);
//        
//        if(_visible){
//          _optgroups_shown[optGroup] = true;
//        }
//      });
//
//      // Build the form model
//      // Build the form template; manually lay out the fields/optgroups
//      var formSchema= {};
//      var FormFields = Backbone.Model.extend({
//        schema: formSchema
//      });
//      var formFields = new FormFields();
//      
//      // first, the fields not in an optGroup
//      _.each(_optGroups[defaultScope].keys,function(key){
//        formSchema[key] = formFieldSchema = { 
//          title: _fields[key]['title'], 
//          key:  key, 
//          type: 'Checkbox',
//          help: _fields[key]['description'],
//          template: fieldCheckboxTemplate 
//        };
//
//        form_template.push( 
//          _.template(field_template)({ 
//              editorId: key+'-id', 
//              title: _fields[key]['title'],
//              name: key 
//            })
//        );
//        if(already_visible[key]){
//          formFields.set( key, true);
//        }
//      });
//      _.each(_.keys(_optGroups),function(optGroup){
//        if (optGroup == defaultScope ) return;
//        
//        form_template.push(
//          _.template(optgroupSelectionTemplate)( 
//            {
//              id: optGroup,
//              name: _optGroups[optGroup].title,
//              help: _optGroups[optGroup].description
//            })
//        );
//        _.each(_optGroups[optGroup].keys,function(key){
//          formSchema[key] = formFieldSchema = { 
//            title: _fields[key]['title'], 
//            key:  key, 
//            type: 'Checkbox',
//            help: _fields[key]['description'],
//            template: optgroupFieldCheckboxTemplate 
//          };
//          if (optGroup == 'otherscreen'){
//            formSchema[key]['template'] = optgroupFieldCheckboxTemplate1;
//          }
//          
//          form_template.push( 
//            _.template(field_template)({ name: key }) );
//          if(already_visible[key]){
//            _optgroups_shown.push(optGroup);
//            formFields.set( key, true);
//          }
//        });
//      });
//        
//      form_template.push('</form>');
//
//      var form = new Backbone.Form({
//        model: formFields,
//        template: _.template(form_template.join(''))
//      });
//
//      form.events = {
//        'click .btn#select-all': function(){
//          $("form input:checkbox").each(function(){
//            $(this).prop("checked",true);
//          });
//        },
//        'click .btn#clear-all': function(){
//          $("form input:checkbox").each(function(){
//            $(this).prop("checked",false);
//          });
//        },
//        // click on the optgroup text expands
//        'click label.selection-group': function(e){
//          _.each(_optGroups[e.target.id].keys, function(key){
//            form.$el.find('[name='+key+']').closest('.form-group').toggle();
//          });
//        },
//        // click on the optgroup checkbox checks all sub boxes
//        'click input.selection-group': function(e){
//          var optGroup = e.target.id.split('-')[0];
//          console.log('id: ' + e.target.id + ', ' + optGroup + ', ' + e.target.checked );
//          _.each(_optGroups[optGroup].keys, function(key){
//            form.$el.find('[name='+key+']').closest('.form-group').show();
//            form.setValue(key, e.target.checked );
//          });
//        },
//      };
//      
//      var _form_el = form.render().el;
//      
//      // if any fields in the optGroup are shown, toggle the optGroup
//      _.each(_optgroups_shown,function(optGroup){
//        _.each(_optGroups[optGroup].keys,function(key){
//          form.$el.find('[name='+key+']').closest('.form-group').toggle();
//        });
//      });
//      
//      // if all fields in the optGroup are shown, check the optGroup
//      _.each(_optgroups_shown, function(optGroup){
//        if(_.every(_optGroups[optGroup].keys, function(sub_key){
//          return formFields.get(sub_key);
//        })){
//          $(_form_el).find(
//              '#' + optGroup + '-checkbox').prop("checked", true);
//        }
//      });
//      
//      appModel.showModal({
//        cancelText: 'Cancel',
//        okText: 'Ok',
//        ok: function(){
//
//          var new_includes = [];
//          var other_screens =[];
//        
//          form.commit();
//          if(_.isUndefined(
//              _.find(formFields.values(),function(val){ return val==true; }))){
//            // TODO: display "nothing selected" error
//            console.log('error: nothing selected');
//            self.select_columns();
//            return;
//          }
//          
//          _.each(formFields.keys(), function(key){
//            var value = formFields.get(key);
//            var max_ordinal = _.max(_fields, function(field){
//              return field['ordinal'];
//            });
//            var field = _fields[key];
//            if(value && !already_visible[key] ){
//              if (_.result(field, 'is_screen_column')===true){
//                new_includes.unshift(key); 
//                other_screens.push(field['screen_facility_id']);
//                value = false; // stop processing this value
//              }else{
//                // 20161213 - Moved to listModel listener
////                self.grid.insertColumn(
////                    Iccbl.createBackgridColumn(
////                        key,field,
////                        self.collection.state.orderStack));
//              }
//            }
//            if(!value && default_visible[key]){
//              new_includes.unshift('-' + key);
//              
//              // 20161213 - Moved to listModel listener
////              column =  self.grid.columns.find(function(column){
////                if(column.get('name') == key){
////                  self.grid.removeColumn(column);
////                  return true;
////                }
////              });
//            }else if(!value){
//              // 20161213 - Moved to listModel listener
////              column =  self.grid.columns.find(function(column){
////                if(column.get('name') == key){
////                  self.grid.removeColumn(column);
////                  return true;
////                }
////              });
//            }
//            if(value && !default_visible[key]){
//              new_includes.unshift(key);
//            }
//          });
//          
//          if (!_.isEmpty(other_screens)){
//            self.show_other_screens(other_screens);
//          }else{
//            self.listModel.set({'includes': new_includes });
//            // New 20161213 - use listmodel event listener instead            
//            //self.collection.fetch();
//            //
//            //// trigger an event to notify new header forms to self-display
//            //self.collection.trigger("MyServerSideFilter:search", 
//            //  self.listModel.get(appModel.URI_PATH_SEARCH), self.collection);
//          }
//        },
//        view: _form_el,
//        title: 'Select columns'  
//      });
//      
//    },
    
//    show_other_screens: function(other_screens){
//      console.log('show_other_screens', other_screens);
//      var self = this;
//      var schemaUrl = [self._options.url,'schema'].join('/');
//      var orderStack = self.listModel.get('order') || [];
//      appModel.getResourceFromUrl(schemaUrl, function(newSchema){
//        var count = 0;
//        _.each(_.values(newSchema['fields']),function(newField){
//          if (_.result(newField,'is_datacolumn')===true){
//            if (other_screens === newField['screen_facility_id']
//                || _.contains(other_screens, newField['screen_facility_id'])){
//              self.grid.insertColumn(
//                  Iccbl.createBackgridColumn(
//                      newField['key'],newField,
//                      orderStack));
//            }
//          }
//        });
//
//        self._options.resource = newSchema;                
//        var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
//        searchHash['other_screens'] = other_screens;
//        self.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
//        // New 20161213 see listmodel change listener
//        //self.collection.fetch();
//        //
//        //// trigger an event to notify new header forms to self-display
//        //self.collection.trigger("MyServerSideFilter:search", 
//        //  self.listModel.get(appModel.URI_PATH_SEARCH), self.collection);
//        
//      },
//      { other_screens: other_screens });
//    },

//    /**
//     * Special function for screen result lists
//     * 20170905 - Deprecated: other screen result columns are now selected using
//     * the treeSelector
//     */
//    show_mutual_positives: function(screen_facility_id, show_mutual_positives){
//      var self = this;
//      var _fields = this._options.resource.fields;
//      var orderStack = self.listModel.get('order') || [];
//
//      if (show_mutual_positives){
//
//        var schemaUrl = [self._options.url,
//                         'schema'].join('/');
//        appModel.getResourceFromUrl(schemaUrl, function(newSchema){
//          var count = 0;
//          _.each(_.values(newSchema['fields']),function(newField){
//            if (_.result(newField,'is_datacolumn')===true){
//              if (screen_facility_id != newField['screen_facility_id']){
//                console.log('insert: ', newField);
//                self.grid.insertColumn(
//                    Iccbl.createBackgridColumn(
//                        newField['key'],newField,
//                        orderStack));
//              }
//            }
//          });
//  
//          self._options.resource = newSchema;                
//          var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
//          searchHash['show_mutual_positives'] = 'true';
//          self.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
//
//          // New 20161213 - use listmodel event listener instead            
//          //self.collection.fetch();
//          //
//          //// trigger an event to notify new header forms to self-display
//          //self.collection.trigger("MyServerSideFilter:search", 
//          //  self.listModel.get(appModel.URI_PATH_SEARCH), self.collection);
//          
//        },
//        { show_mutual_positives: true});
//      
//      
//      
//      }else{
//        var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
//        _.each(_.pairs(_fields), function(pair){
//          var key = pair[1]['key'];
//          var prop = pair[1];
//          var fieldType = prop['scope'].split('.')[0]
//          var field_screen_facility_id = _.result(prop,'screen_facility_id', '');
//          // Note: if filtering/ordering on one of the mutual positive columns, 
//          // do not remove it here.
//          if(fieldType == 'datacolumn'
//              && field_screen_facility_id != screen_facility_id 
//              && ! _.findKey(searchHash, function(val,hashkey){
//                return hashkey.indexOf(key) > -1
//              })
//              && ! _.find(orderStack, function(orderkey){
//                return orderkey.indexOf(key) > -1
//              })
//            ){
//            console.log('remove: ', key);
//            prop['visibility'] = _.without(prop['visibility'], 'l');
//            var currentColumn = self.grid.columns.findWhere({ name: key });
//            if(currentColumn){
//              self.grid.removeColumn(currentColumn);
//            }
//          }
//        });
//        if(_.has(searchHash,'show_mutual_positives')){
//          delete searchHash['show_mutual_positives'];
//          self.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
//        }        
//        
//      }    
//    },
    
//    /**
//     * Special function for screen result lists
//     */
//    show_mutual_positives_bak: function(screen_facility_id, show_mutual_positives){
//      var self = this;
//      var _fields = this._options.resource.fields;
//      var orderStack = self.listModel.get('order') || [];
//
//      if (show_mutual_positives){
//        
//        _.each(_.pairs(_fields), function(pair){
//          var key = pair[1]['key'];
//          var prop = pair[1];
//          var fieldType = prop['scope'].split('.')[0];
//          var field_screen_facility_id = _.result(prop,'screen_facility_id', '');
//          
//          if(fieldType == 'datacolumn' 
//            && field_screen_facility_id != screen_facility_id ){
//            prop['visibility'].push('l');
//            var currentColumn = self.grid.columns.findWhere({ name: key });
//            if(! currentColumn){
//              self.grid.insertColumn(
//                Iccbl.createBackgridColumn(
//                    key,prop,
//                    orderStack));
//            }
//          }
//        });
//        
//        var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
//        searchHash['show_mutual_positives'] = 'true';
//        self.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
//      }else{
//        var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
//        _.each(_.pairs(_fields), function(pair){
//          var key = pair[1]['key'];
//          var prop = pair[1];
//          var fieldType = prop['scope'].split('.')[0]
//          var field_screen_facility_id = _.result(prop,'screen_facility_id', '');
//          // Note: if filtering/ordering on one of the mutual positive columns, 
//          // do not remove it here.
//          if(fieldType == 'datacolumn'
//              && field_screen_facility_id != screen_facility_id 
//              && ! _.findKey(searchHash, function(val,hashkey){
//                return hashkey.indexOf(key) > -1
//              })
//              && ! _.find(orderStack, function(orderkey){
//                return orderkey.indexOf(key) > -1
//              })
//            ){
//            console.log('remove: ', key);
//            prop['visibility'] = _.without(prop['visibility'], 'l');
//            var currentColumn = self.grid.columns.findWhere({ name: key });
//            if(currentColumn){
//              self.grid.removeColumn(currentColumn);
//            }
//          }
//        });
//        if(_.has(searchHash,'show_mutual_positives')){
//          delete searchHash['show_mutual_positives'];
//          self.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
//        }        
//      }
//    },
    
  });
  

  return ListView;
});
