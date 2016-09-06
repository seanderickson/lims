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
  'views/plateLocation',
  'views/user/user2',
  'views/user/screensaveruser',
  'views/usergroup/usergroup2',
  'test/detailTest',
  'templates/content.html',
  'templates/welcome.html',
  'templates/about.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         EditView, LibraryView, LibraryCopyView, LibraryCopyPlateView, 
         ScreenView, LibraryScreeningView, PlateLocationView, UserAdminView, 
         UserView, UserGroupAdminView, DetailTestView, layout, welcomeLayout, 
         aboutLayout) {
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout, 
    'LibraryView': LibraryView,
    'LibraryCopyView': LibraryCopyView, 
    'LibraryCopyPlateView': LibraryCopyPlateView,
    'ScreenView': ScreenView,
    'LibraryScreeningView': LibraryScreeningView,
    'PlateLocationView': PlateLocationView,
    'UserView': UserView,
    'UserAdminView': UserAdminView,
    'UserGroupAdminView': UserGroupAdminView,
    'DetailTestView': DetailTestView
  };
    
  var ContentView = Iccbl.UriContainerView.extend({
    
    template: _.template(layout),
    
    initialize: function() {
      console.log('initialize content.js', arguments);
      Iccbl.UriContainerView.prototype.initialize.apply(this,arguments);
      this.title = null;
      this.consumedStack = [];
      this.objects_to_destroy = [];
    },
        
    // Create a new model for editing
    showAdd: function(resource, uriStack){
      
      var self = this;
      var newModel = appModel.createNewModel(resource.key);
      var viewClass = DetailLayout;
      var view;
      
      newModel.resource = resource;
      this.$('#content_title').html(resource.title + ': Add' );
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
      view = new viewClass({ model: newModel, uriStack: uriStack});

      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
    },
    
    // Show a mock detail view to test UI components
    showDetailTest: function(uriStack){
      var self = this;
      this.$('#content_title').html("Detail Test");
      var view = new DetailTestView({uriStack: uriStack});
      self.setView('#content', view).render();
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
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
    
    }, // end showDetail
    
    // Show a listing of a resource, get the parameters from the uriStack
    showList: function(resource, uriStack, schemaResult) {
      
      console.log('showList: uriStack', uriStack);

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
      
      if(uriStack.length > 1 && uriStack[0] == 'children'){
        var substack = _.rest(uriStack,1);
        var _key = Iccbl.popKeyFromStack(resource, substack, []);
        this.$('#content_title').html('Child logs for: ' + _key);
      }else{
        this.$('#content_title').html(resource.title + ' listing');
      }
    
      if(uriStack.length > 1 && uriStack[0] == 'search' 
          && !_.isNaN(parseInt(uriStack[1]))){

        view = self.createListSearchView(resource, schemaResult, viewClass, uriStack);
      
      }else{ // normal list view
      
        view = new viewClass({ 
            model: appModel, 
            uriStack: uriStack,
            schemaResult: schemaResult, 
            resource: resource
          });
      }
    
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

    }, // end showList
    
    // Create a list view using the stateful searchID given on the uriStack
    createListSearchView: function(resource, schemaResult, viewClass, uriStack){
      
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
        fetch: function(options){
          // Endcode for the post_data arg; post data elements are sent 
          // as urlencoded values via a POST form for simplicity
          options.data = _.extend(
            { search_data: JSON.stringify(search_data) }, 
            options.data);
          options.type = 'POST';
          return Iccbl.MyCollection.prototype.fetch.call(this, options);
        }
      });
      var collection = new Collection({
        url: url,
      });
      var view = new viewClass({ 
          model: appModel, 
          uriStack: uriStack,
          schemaResult: schemaResult, 
          resource: resource, 
          collection: collection, 
          search_data: search_data
        });
      this.$('#content_title').html(resource.title + ' search');
      return view;
    },
    
    // Backbone layoutmanager callback
    cleanup: function(){
      console.log('cleanup');
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
    },

    // Main view control method
    // delegated from UriContainerView on appModel change:uriStack
    changeUri: function(uriStack) {
      
      console.log('changeUri:', uriStack);
      
      var self = this;
      var consumedStack = this.consumedStack = [];
      var uiResourceId;
      var resource;
      
      self.removeView('#content');
      self.cleanup();
      self.off();
      
      if (!_.isEmpty(uriStack)){
        if(uriStack[0] == 'detail_test'){
          this.consumedStack.push(uriStack.shift());
          self.showDetailTest(uriStack);
          return;
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
        self.showList(resource, uriStack,resource);
      }
    }
  });

  return ContentView;
});