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
  'views/library/libraryWell',
  'views/screen/screen',
  'views/screen/libraryScreening',
  'views/screen/cherryPickRequest',
  'views/plateLocation',
  'views/user/user2',
  'views/user/screensaveruser',
  'views/usergroup/usergroup2',
  'views/activityListView',
  'views/apilogView',
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
         LibraryWellsView, LibraryWellView,
         ScreenView, LibraryScreeningView, CherryPickRequestView,
         PlateLocationView, UserAdminView, UserView, UserGroupAdminView, 
         ActivityListView, ApilogView, UploadDataForm, DetailTestView, 
         WellSelectorView, SearchBox, layout, welcomeLayout, 
         welcomeScreenerLayout, aboutLayout) {
  
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
    'LibraryWellsView': LibraryWellsView,
    'LibraryWellView': LibraryWellView,
    'ApilogView': ApilogView
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
      this.$('#content_title').html('Create ' + resource.title );
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
    showList: function(resource, uriStack) {
      
      console.log('showList: resource:', resource, 'uriStack', uriStack); 
      var self = this;
      var uriStack = _.clone(uriStack);
      var viewClass = ListView;
      var view;

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
      var url = resource.apiUri;
      if(uriStack.length > 1 && uriStack[0] == 'children'){
        
        self.consumedStack.push(uriStack.shift());
        
        var _key = Iccbl.popKeyFromStack(resource, uriStack, self.consumedStack);
//        url += '/children/' + _key;
        url = [resource.apiUri, 'children', _key].join('/');
        
        
        this.$('#content_title').html('Child logs for: ' + _key);
        this.$('#content_title_row').show();
      }
      else{
//        if (!_.isEmpty(_.intersection(uriStack, 
//              [appModel.URI_PATH_COMPLEX_SEARCH,appModel.URI_PATH_ENCODED_SEARCH]))){
//          this.$('#content_title').html(resource.title + ' listing');
//          this.$('#content_title_row').show();
//        }
      }
      
      console.log('collection url', url);
      var Collection = Iccbl.MyCollection.extend({
        url: url
      });
      var collection = new Collection();
      var extraControls = [];
      
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
        if (resource.key in ['reagent','well','copywell',
                             'librarycopy','librarycopyplate']){
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
      
      self.listenTo(view, 'update_title', function(val){
        if(val) {
          this.$('#content_title').html(resource.title + ': <small>' + val + '</small>');
          this.$('#content_title_row').show();
        }else{
          if (!_.contains(['well', 'smallmoleculereagent','silencingreagent'],resource.key)){
            this.$('#content_title').html(resource.title);
            this.$('#content_title_row').show();
          }
        }
      });
    
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#content', view ).render();
      self.objects_to_destroy.push(view);

    }, // end showList

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
      this.off();
      appModel.clearErrors();

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
      if (!_.isEmpty(uriStack) 
            && !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
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