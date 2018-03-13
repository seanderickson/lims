define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_detail_layout',
  'views/generic_edit',
  'views/library/library',
  'views/library/libraryCopy',
  'views/library/libraryCopyPlate',
  'views/library/libraryWells',
  'views/screen/screen',
  'views/screen/libraryScreening',
  'views/screen/cherryPickRequest',
  'views/plateLocation',
  'views/user/user2',
  'views/user/screensaveruser',
  'views/usergroup/usergroup2',
  'views/activityListView',
  'utils/uploadDataForm',
  'test/detailTest',
  'utils/wellSelector',
  'views/search_box',
  'templates/content.html',
  'templates/welcome.html',
  'templates/welcome-screener.html',
  'templates/about.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         EditView, LibraryView, LibraryCopyView, LibraryCopyPlateView,
         LibraryWellsView,
         ScreenView, LibraryScreeningView, CherryPickRequestView,
         PlateLocationView, UserAdminView, 
         UserView, UserGroupAdminView, ActivityListView, UploadDataForm, DetailTestView, 
         WellSelectorView, SearchBox, layout, 
         welcomeLayout, welcomeScreenerLayout, aboutLayout) {
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout, 
    'LibraryView': LibraryView,
    'LibraryCopyView': LibraryCopyView, 
    'LibraryCopyPlateView': LibraryCopyPlateView,
    'ScreenView': ScreenView,
    'LibraryScreeningView': LibraryScreeningView,
    'CherryPickRequestView': CherryPickRequestView,
    'PlateLocationView': PlateLocationView,
    'UserView': UserView,
    'UserAdminView': UserAdminView,
    'UserGroupAdminView': UserGroupAdminView,
    'DetailTestView': DetailTestView,
    'WellSelectorView': WellSelectorView,
    'ActivityListView': ActivityListView,
    'LibraryWellsView': LibraryWellsView
  };
  
  var ContentView = Iccbl.UriContainerView.extend({
    
    template: _.template(layout),
    
    initialize: function() {
      console.log('initialize content.js', arguments);
      Iccbl.UriContainerView.prototype.initialize.apply(this,arguments);
      this.title = null;
      this.consumedStack = [];
      this.objects_to_destroy = [];
      _.bindAll(this, 'cleanup');
    },
        
    // Create a new model for editing
    showAdd: function(resource, uriStack){
      
      var self = this;
      var viewClass = DetailLayout;
      var view;
      
      var newModel = appModel.createNewModel(resource.key);
      newModel.resource = resource;
      this.$('#content_title').html('Create a new ' + resource.title );
      this.$('#content_title_row').show();
      
      if (_.has(resource, 'detailView')){
        if (_.has(VIEWS, resource['detailView'])) {
          viewClass = VIEWS[resource['detailView']];
          console.log('found view ' + resource['detailView']);
        }else{
          var msg = 'detailView class specified could not be found: ' + 
              resource['detailView'];
          console.log(msg);
          throw msg;
        }
      }
      view = new viewClass({
        model: newModel, 
        uriStack: uriStack,
        isCreate: true,
      });

      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
      self.objects_to_destroy.push(view);
    },
    
    // Show a mock detail view to test UI components
    showDetailTest: function(uriStack){
      var self = this;
      this.$('#content_title').html("Detail Test");
      this.$('#content_title_row').show();
      var view = new DetailTestView({uriStack: uriStack});
      self.setView('#content', view).render();
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.objects_to_destroy.push(view);
    },
    
    // Show a mock well selector view to test UI components
    showWellSelectorTest: function(uriStack){
      var self = this;
      this.$('#content_title').html("Well Selector Test");
      this.$('#content_title_row').show();
      var view = new WellSelectorView({
        uriStack: uriStack
      });
      self.setView('#content', view).render();
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.objects_to_destroy.push(view);
    },
    
    // Show the standard view of a resource
    showDetail: function(uriStack, model){
      
      console.log('showDetail', uriStack, model.resource.key);
      
      var self = this;
      var uriStack = _.clone(uriStack);
      var resource = model.resource;
      var viewClass = DetailLayout;

      if (_.has(resource, 'detailView')){
        if (_.has(VIEWS, resource['detailView'])) {
          viewClass = VIEWS[resource['detailView']];
          console.log('found view ' + resource['detailView']);
        }else{
          var msg = 'detailView class specified could not be found: ' + 
              resource['detailView'];
          console.log(msg);
          throw msg;
        }
      }
      
      var view = new viewClass({ model: model, uriStack: uriStack});
      var titleFunc = function setContentTitle(){
        var title = model.resource.title + ': ' + 
            Iccbl.getTitleFromTitleAttribute(model,model.resource);          
        if (_.isFunction(view.getTitle)){
          title = view.getTitle();
        }
        self.$('#content_title').html(title);
        self.$('#content_title_row').show();
      }
      titleFunc();
      model.on('sync', function(model){
        // TODO: it would be better to watch just the title attribute
        titleFunc();
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
      self.objects_to_destroy.push(view);
    
    }, // end showDetail
    
    // Show a listing of a resource, get the parameters from the uriStack
    showList: function(resource, uriStack, complex_search) {
      
      console.log('showList: resource:', resource, 'uriStack', uriStack, 
        'complex search', complex_search);

      var self = this;
      var uriStack = _.clone(uriStack);
      var viewClass = ListView;
      var view;
      var url = resource.apiUri;
      
      if (_.has(resource, 'listView')){
        if (_.has(VIEWS, resource['listView'])) {
          viewClass = VIEWS[resource['listView']];
          console.log('found view ' + resource['listView']);
        }else{
          var msg = 'listView class specified could not be found: ' + 
              resource['listView'];
          window.alert(msg);
          throw msg;
        }
      }
      
      if(uriStack.length > 1 && uriStack[0] == 'children'){
        
        self.consumedStack.push(uriStack.shift());
        
        var _key = Iccbl.popKeyFromStack(resource, uriStack, self.consumedStack);
        url += '/children/' + _key;
        
        this.$('#content_title').html('Child logs for: ' + _key);
        this.$('#content_title_row').show();
      }else{
        this.$('#content_title').html(resource.title + ' listing');
        this.$('#content_title_row').show();
      }
    
      if (!_.isEmpty(complex_search)){
        view = self.createListSearchView(resource, uriStack, viewClass, complex_search);
      
      }else{ // normal list view

        var Collection = Iccbl.MyCollection.extend({
          url: url
        });
        collection = self.collection = new Collection();
        var extraControls = [];
        
        if (false){
          // FIXME: not all types can have an "add"
          //if (appModel.hasPermission(resource.key, 'write')){
          
          if (_.contains(resource.visibility, 'c')){
            var showAddButton = $([
              '<a class="btn btn-default btn-sm pull-down" ',
                'role="button" id="add_resource" href="#">',
                'Add</a>'
              ].join(''));   
            showAddButton.click(function(e){
              e.preventDefault();
              var route = resource.key + '/+add';
              appModel.router.navigate(route, {trigger: true});
            });
            extraControls.push(showAddButton);
          }
          if (_.contains(resource.visibility, 'e')){
            var showUploadButton = $([
              '<a class="btn btn-default btn-sm pull-down" ',
                'role="button" id="patch_resource" href="#">',
                'Upload data</a>'
              ].join(''));   
            showUploadButton.click(function(e){
              e.preventDefault();
              UploadDataForm.postUploadFileDialog(
                collection.url, resource.content_types)
                .done(function(){
                  collection.fetch({ reset: true });
                })
                .fail(function(){
                  appModel.jqXHRfail.apply(this,arguments); 
                });
            });
            extraControls.push(showUploadButton);
          }
        }
        view = new viewClass({ 
            uriStack: uriStack,
            resource: resource,
            collection: collection,
            extraControls: extraControls
          });
      }
      
      // FIXME: createListSearchView does not return a view; should return a deferred
      if (!_.isUndefined(view)){
        self.listenTo(view, 'update_title', function(val){
          if(val) {
            this.$('#content_title').html(resource.title + ': <small>' + val + '</small>');
            this.$('#content_title_row').show();
          }else{
            this.$('#content_title').html(resource.title);
            this.$('#content_title_row').show();
          }
        });
      
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView('#content', view ).render();
        self.objects_to_destroy.push(view);
      }

    }, // end showList

    /**
     * Prepare the list view to send API_PARAM_SEARCH search data that has 
     * been encoded in the URI stack, or that is referenced by the 
     * API_PARAM_SEARCH_ID in localstorage
     * TODO: 20180312 - refactor this to the list view class.
     */
    createListSearchView: function(resource, uriStack, viewClass, complex_search){
      console.log('createListSearchView...')
      var self = this;
      var resources_supporting_raw_search = [
        'well','reagent','compound_search','librarycopyplate'];
      if (!_.contains(resources_supporting_raw_search,resource.key)){
        throw "search data not supported for resource: " + resource.key + 
          ', must be one of: ' + resources_supporting_raw_search.join(',');
      }
      
      function showView(schemaResult){
        // Update the current resource and fields:
        // TODO: modify the api: sirna, rna resource schema with default visibilities
        // - current resource prop and field definitions override the new 
        // definitions, if defined, but the new fields are used if not.
        var newFields = _.extend({}, schemaResult['fields'], resource['fields'])
        var newResource = _.extend({}, schemaResult, resource );
        newResource['fields'] = newFields;
        
        var url = newResource.apiUri;
        
        var includes = [];
        if (newResource.key == 'compound_search'){
          url += '/compound_search';
          // FIXME: hack to add columns; fix is to implement sirna/smr resource
          // schema as superset of reagent schema
          newResource['options']['includes'] = [
           'inchi','smiles','structure_image','molecular_formula',
           'molecular_mass','molecular_weight','pubchem_cid','chembank_id',
           'chembl_id'];
        }
        
        var viewOptions = {
          url: url,
          resource: newResource, 
          extraControls: extraControls
        };
        
        // If Complex Search, then search data are sent as a post header
        if (_.has(complex_search, appModel.API_PARAM_SEARCH_ID)){
          var searchId = complex_search[appModel.API_PARAM_SEARCH_ID];
          self.consumedStack.push(appModel.URI_PATH_COMPLEX_SEARCH);
          self.consumedStack.push(searchId);
          url += '/' + [appModel.URI_PATH_COMPLEX_SEARCH, searchId].join('/');
          
          uriStack = _.result(complex_search,'uriStack',[]);
          viewOptions['uriStack'] = uriStack;
          // NOTE: setting API_PARAM_SEARCH on the list options:
          // - for libraryCopyPlate.js
          // TODO: 20180312 - refactor list.js to process complex search
          viewOptions[appModel.API_PARAM_SEARCH] = 
            complex_search[appModel.API_PARAM_SEARCH];
          // Set up the collection add extra API_PARAM_SEARCH to POST data
          var Collection = Iccbl.MyCollection.extend({
            url: url,
            fetch: function(options){
              var options = options || {};
              options.data = _.extend({}, options.data);
              options.data[appModel.API_PARAM_SEARCH] = 
                complex_search[appModel.API_PARAM_SEARCH];
              options.type = 'POST';
              return Iccbl.MyCollection.prototype.fetch.call(this, options);
            }
          });
          var collection = new Collection();
          viewOptions['collection'] = collection;
        } else {
          self.consumedStack.push(appModel.URI_PATH_ENCODED_SEARCH);
          self.consumedStack.push(complex_search[appModel.API_PARAM_SEARCH]);
          viewOptions[appModel.API_PARAM_ENCODED_SEARCH] = 
            complex_search[appModel.API_PARAM_SEARCH];
          fetchType = 'GET';
          uriStack = _.result(complex_search,'uriStack',[]);
          viewOptions['uriStack'] = uriStack;
          // Set up the collection to add extra API_PARAM_SEARCH to GET data
          var Collection = Iccbl.MyCollection.extend({
            url: url,
            fetch: function(options){
              var options = options || {};
              options.data = _.extend({}, options.data);
              options.data[appModel.API_PARAM_SEARCH] = 
                complex_search[appModel.API_PARAM_SEARCH];
              options.type = 'GET';
              return Iccbl.MyCollection.prototype.fetch.call(this, options);
            }
          });
          var collection = new Collection();
          viewOptions['collection'] = collection;
        }
  
        
        var extraControls = [];
        var $modifySearch = $([
          '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="modify_search_button" ',
          'href="#">',
          'Modify Search</a>'
        ].join(''));
        extraControls.push($modifySearch);
        $modifySearch.click(function(e){
          e.preventDefault();
          self.showSearch(newResource, complex_search);
        });
        viewOptions['extraControls'] = extraControls;
        var view = new viewClass(viewOptions);
        
        self.$('#content_title_row').show();
          
        self.listenTo(view, 'update_title', function(val){
          if(val) {
            self.$('#content_title').html(newResource.title + ': <small>' + val + '</small>');
          }else{
            self.$('#content_title').html(newResource.title);
          }
        });
        
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView('#content', view ).render();
        self.objects_to_destroy.push(view);
      }; // showView
      
      var schemaUrl = resource.apiUri + '/schema';
      appModel.getResourceFromUrl(schemaUrl, showView);
      
      // FIXME: should return a deferred call
      //      return null;
    },
    
    /** 
     * Show/modify the reagent search list view 
     **/
    showSearch: function(resource, complex_search) {
      console.log('modify search:', resource.key, complex_search);

      function parse(value, errors){
        var parsedData;
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
        var errors;
        if(resource.key == 'well' || resource.key == 'reagent'){
          errors = SearchBox.prototype.processWellSearch.call(this, text_to_search);
        }else if (resource.key == 'compound_search'){
          errors = SearchBox.prototype.processCompoundSearch.call(this, text_to_search);
        }else if (resource.key == 'librarycopyplate'){
          errors = SearchBox.prototype.processCopyPlateSearch.call(this,text_to_search);
        }else{
          throw 'Unknown search resource: ' + resource.key;
        }
        if (errors){
          throw 'Unexpected error on submit: ' + errors.join(', ');
        }
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

      var search_text = complex_search[appModel.API_PARAM_SEARCH];
      var errors = [];
      var parsedData = parse(search_text, errors);
      if (!_.isEmpty(errors)){
        throw 'errors in parsed search:' + errors.join(', ');
      } else {
        search_text = _.map(parsedData, function(parsedLine){
          if (_.has(parsedLine, 'combined')){
            return parsedLine.combined.join(' ');
          }else{
            return parsedLine;
          }
        }).join('\n');
      }
      var FormFields = Backbone.Model.extend({
        schema: schema
      });
      var formFields = new FormFields({
        searchVal: search_text
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
    
    // Backbone layoutmanager callback
    cleanup: function(){
      console.log('cleanup:', this.objects_to_destroy.length);

      var oldView = this.getView('#content');
      if(!_.isUndefined(oldView)){
        oldView.off();
      }
      this.removeView('#content');
      this.off(); //this.unbind(); 
      _.each(this.objects_to_destroy, function(obj){
        obj.remove();
        obj.off();
      });
      this.$('#content_title_message').empty();
      this.$('#content_title').empty();
      this.objects_to_destroy = [];
    },

    // Main view control method
    // delegated from UriContainerView on appModel change:uriStack
    changeUri: function(uriStack) {
      
      console.log('changeUri:', uriStack);
      
      var self = this;
      var consumedStack = this.consumedStack = [];
      var uiResourceId;
      var resource;
      
      self.cleanup();
      
      if (!_.isEmpty(uriStack)){
        if(uriStack[0] == 'detail_test'){
          this.consumedStack.push(uriStack.shift());
          self.showDetailTest(uriStack);
          return;
        }else if (uriStack[0] == 'wellselectortest'){
          this.consumedStack.push(uriStack.shift());
          self.showWellSelectorTest(uriStack);
        }
        uiResourceId = uriStack.shift();
        this.consumedStack.push(uiResourceId);
      }else{
        uiResourceId = 'home';
      }
      
      resource = appModel.getResource(uiResourceId);
      
      if (uiResourceId == 'home'){
        
        var currentUser = appModel.getCurrentUser();
        var template = welcomeLayout;
        if(!currentUser.is_staff && !currentUser.is_superuser ){
          template = welcomeScreenerLayout;
        }
        
        var WelcomeView = Backbone.Layout.extend({
          template: _.template(template), 
          serialize: function(){
            var appData = appModel.getAppData();
            return {
              user_fullname: currentUser.first_name + ' ' + currentUser.last_name,
              app_data: appData.toJSON()
            };
          }
        });
        $('#navbar').children().removeClass('active');
        $('#navbar').children('#home').addClass('active');
        var view = new WelcomeView();
        self.setView('#content', view).render();
        return;
      }
      
      if (uiResourceId == 'about'){
        var WelcomeView = Backbone.Layout.extend({
          template: _.template(aboutLayout)
        });
        $('#navbar').children().removeClass('active');
        $('#navbar').children('#about').addClass('active');
        var view = new WelcomeView();
        self.setView('#content', view).render();
        return;
      }
      var complex_search = appModel.findComplexSearch(uriStack);
      if (_.has(complex_search,'errors')){
        var errMsg = complex_search['errors'];
        console.log(errMsg);
        appModel.error(errMsg);
        complex_search = null;
      }
      
      if (!_.isEmpty(uriStack) 
            && !_.contains(appModel.LIST_ARGS, uriStack[0]) 
            && _.isEmpty(complex_search)) {
        // DETAIL VIEW
        
        if(uriStack[0] == '+add'){
          self.showAdd(resource, uriStack);
        }else{ 
          try{
            var _key = Iccbl.popKeyFromStack(resource, uriStack, consumedStack );
            var options = {};
            if (uiResourceId == 'screen'){
              // Use the special "ui" url for screen
              options.url = [resource.apiUri, _key, 'ui'].join('/');
            }
            appModel.getModel(
              uiResourceId, _key, 
              function(model){
                model.resource = resource;
                self.showDetail(uriStack, model);
              }, 
              options);
          }catch(e){
            var msg = 'Unable to display resource: ' + uiResourceId;
            console.log(msg,e);
            appModel.error(msg);
          }
        }
      } else {
        // LIST VIEW
        self.showList(resource, uriStack, complex_search);
      }
    }
  });

  return ContentView;
});