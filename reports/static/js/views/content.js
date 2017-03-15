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
  'views/screen/screen',
  'views/screen/libraryScreening',
  'views/screen/cherryPickRequest',
  'views/plateLocation',
  'views/user/user2',
  'views/user/screensaveruser',
  'views/usergroup/usergroup2',
  'utils/uploadDataForm',
  'test/detailTest',
  'utils/wellSelector',
  'templates/content.html',
  'templates/welcome.html',
  'templates/about.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         EditView, LibraryView, LibraryCopyView, LibraryCopyPlateView, 
         ScreenView, LibraryScreeningView, CherryPickRequestView,
         PlateLocationView, UserAdminView, 
         UserView, UserGroupAdminView, UploadDataForm, DetailTestView, 
         WellSelectorView, layout, 
         welcomeLayout, aboutLayout) {
  
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
    'WellSelectorView': WellSelectorView
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
      var newModel = appModel.createNewModel(resource.key);
      var viewClass = DetailLayout;
      var view;
      
      newModel.resource = resource;
      this.$('#content_title').html('Create a new ' + resource.title );
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
      function saveSuccessCallBack(model){
        console.log('new model created: ', model);
        if (resource.key == 'screen'){
          appModel.unset('screens');
        }
        else if (resource.key == 'library'){
          appModel.unset('libraries');
        }
        else if (resource.key == 'screensaveruser'){
          appModel.unset('users');
        }
        model = new Backbone.Model(model);
        var key = Iccbl.getIdFromIdAttribute( model,resource );
        model.key = resource.key + '/' + key;
        appModel.router.navigate(resource.key + '/' + key, {trigger:true});
      };
      
      view = new viewClass({
        model: newModel, 
        uriStack: uriStack,
        isCreate: true,
        saveSuccessCallBack: saveSuccessCallBack
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
      var view = new DetailTestView({uriStack: uriStack});
      self.setView('#content', view).render();
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.objects_to_destroy.push(view);
    },
    
    // Show a mock well selector view to test UI components
    showWellSelectorTest: function(uriStack){
      var self = this;
      this.$('#content_title').html("Well Selector Test");
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
      
      var titleFunc = function setContentTitle(val){
        self.$('#content_title').html(val);
      }
      
      titleFunc(model.resource.title + ': ' + 
        Iccbl.getTitleFromTitleAttribute(model,model.resource) );
      
      var view = new viewClass({ model: model, uriStack: uriStack});
      model.on('sync', function(model){
        // TODO: it would be better to watch just the title attribute
        titleFunc(model.resource.title + ': ' + 
          Iccbl.getTitleFromTitleAttribute(model,model.resource) );
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
      self.objects_to_destroy.push(view);
    
    }, // end showDetail
    
    // Show a listing of a resource, get the parameters from the uriStack
    showList: function(resource, uriStack) {
      
      console.log('showList: uriStack', uriStack);

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
      }else{
        this.$('#content_title').html(resource.title + ' listing');
      }
    
      if(uriStack.length > 1 && uriStack[0] == 'search' 
          && !_.isNaN(parseInt(uriStack[1]))){

        view = self.createListSearchView(resource, viewClass, uriStack);
      
      }else{ // normal list view

        var Collection = Iccbl.MyCollection.extend({
          url: url
        });
        collection = self.collection = new Collection();
//        
//        collection = self.collection = new Collection({
//          'url': url
//        });
        
        var extraControls = [];
        if (appModel.hasPermission(resource.key, 'write')){
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
//            model: appModel, 
            uriStack: uriStack,
            schemaResult: resource, 
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
          }else{
            this.$('#content_title').html(resource.title);
          }
        });
      
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView('#content', view ).render();
        self.objects_to_destroy.push(view);
      }

    }, // end showList
    
    /**
     * Create a list view using the stateful searchID given on the uriStack. 
     */
    createListSearchView: function(resource, viewClass, uriStack){
      var self = this;
      console.log('createListSearchView')
      var url;
      var searchID;
      var search_data;

      console.log('create a collection with search data: ', uriStack);
      
      this.consumedStack.push(uriStack.shift());
      searchID = uriStack.shift();
      this.consumedStack.push(searchID);
      search_data = appModel.getSearch(searchID);
      
      if(_.isUndefined(search_data) || _.isEmpty(search_data)){
        var msg = 'Content List search requires a "search_data:'
          +searchID +'" in current browser state';
        console.log(msg);
        appModel.error(msg);
        return;
      }
      url = resource.apiUri + '/search/' + searchID;
      var Collection = Iccbl.MyCollection.extend({
        url: url,
        fetch: function(options){
          var options = options || {};
          // Encode for the post_data arg; post data elements are sent 
          // as urlencoded values via a POST form for simplicity
          options.data = _.extend(
            { search_data: JSON.stringify(search_data) }, 
            options.data);
          options.type = 'POST';
          return Iccbl.MyCollection.prototype.fetch.call(this, options);
        }
      });
      var collection = new Collection();
      
      function showView(newResource){
        console.log('got new resource:', newResource, resource);
        
        var extraControls = [];
        
        // Special button to modifiy the reagent search
        // TODO: implement this for other searches, if search_data is complex/large
        if (resource.key == 'reagent'){
          var $modifySearch = $([
            '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="modify_search_button" ',
            'href="#">',
            'Modify Search</a>'
          ].join(''));
          extraControls.push($modifySearch);
          $modifySearch.click(function(e){
            e.preventDefault();
            self.showReagentSearch(resource, search_data);
          });
        }
        
        var view = new viewClass({ 
            uriStack: uriStack,
            schemaResult: newResource, 
            resource: newResource, 
            collection: collection, 
            search_data: search_data,
            extraControls: extraControls
          });
        this.$('#content_title').html(resource.title + ' search');
        
        
        self.listenTo(view, 'update_title', function(val){
          if(val) {
            this.$('#content_title').html(resource.title + ': <small>' + val + '</small>');
          }else{
            this.$('#content_title').html(resource.title);
          }
        });
      
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView('#content', view ).render();
        self.objects_to_destroy.push(view);
        
      }
      var schemaUrl = resource.apiUri + '/schema?search='+ searchID;
      appModel.getResourceFromUrl(schemaUrl, showView);

      // FIXME: should return a deferred call
      //      return null;
    },
    
    /** 
     * Show/modify the reagent search list view 
     **/
    showReagentSearch: function(resource, search_data) {
      var schema = {
        search_data: {
          title: 'Well Search',
          key: 'search_data',
          type: Backbone.Form.editors.TextArea.extend({
              initialize: function() {
                Backbone.Form.editors.TextArea.prototype.initialize.apply(this, arguments);
                this.$el.attr('rows',12);
              }
          }),
          editorClass: 'input-full',
          validators: ['required'],
          template: appModel._field_template //altFieldTemplate
        }
      };
      var FormFields = Backbone.Model.extend({
        schema: schema
      });
      var formFields = new FormFields({
        search_data: search_data
      });
      var form = new Backbone.Form({
        model: formFields,
        template: self._form_template //_.template(form_template.join(''))
      });
      var _form_el = form.render().el;
      appModel.showModal({
        title: 'Search for Wells',
        view: _form_el,
        ok: function() {
          var errors = form.commit();
          if(!_.isEmpty(errors)){
            console.log('form errors, abort submit: ' + JSON.stringify(errors));
            return false;
          }
          var searchID = ( new Date() ).getTime();
          appModel.setSearch(searchID,form.getValue()['search_data']);
          this.searchID = searchID;
          appModel.set('routing_options', {replace: false});  
          var _route = resource.key + '/search/'+ searchID;
          appModel.router.navigate(_route, {trigger:true});
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
        var WelcomeView = Backbone.Layout.extend({
          template: _.template(welcomeLayout)
        });
        $('#navbar').children().removeClass('active');
        $('#navbar').children('#home').addClass('active');
        var view = new WelcomeView();
        self.setView('#content', view).render();
        this.$('#content_title').html(resource.title);
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
        this.$('#content_title').html(resource.title);
        return;
      }
      
      // Test for list args, if not found, then it's a detail view
      if (!_.isEmpty(uriStack) 
            && !_.contains(appModel.LIST_ARGS, uriStack[0]) 
            && uriStack[0] != 'search') {
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
        self.showList(resource, uriStack);
      }
    }
  });

  return ContentView;
});