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
  'views/screen/screen',
  'views/screen/libraryScreening',
  'views/user/user2',
  'views/user/screensaveruser',
  'views/usergroup/usergroup2',
  'test/detailTest',
  'text!templates/content.html',
  'text!templates/welcome.html',
  'text!templates/about.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         EditView, LibraryView, ScreenView, LibraryScreeningView, UserAdminView, 
         UserView, UserGroupAdminView, DetailTestView, layout, welcomeLayout, 
         aboutLayout) {
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout, 
    'LibraryView': LibraryView,
    'ScreenView': ScreenView,
    'LibraryScreeningView': LibraryScreeningView,
    'UserView': UserView,
    'UserAdminView': UserAdminView,
    'UserGroupAdminView': UserGroupAdminView,
    'DetailTestView': DetailTestView
  };
    
  var ContentView = Iccbl.UriContainerView.extend({
    
    template: _.template(layout),
    className: "col-sm-10 col-md-10 col-lg-10",
    
    initialize: function() {
      console.log('initialize content.js', arguments);
      Iccbl.UriContainerView.prototype.initialize.apply(this,arguments);
      this.title = null;
      this.consumedStack = [];
      this.objects_to_destroy = [];
    },
        
    showAdd: function(resource, uriStack){
      var self = this;
      var newModel = appModel.createNewModel(resource.key);
      newModel.resource = resource;
      this.$('#content_title').html(resource.title + ': Add' );
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
      var view = new viewClass({ model: newModel, uriStack: uriStack});

      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
    },
    
    showDetailTest: function(uriStack){
      var self = this;
      self.removeView('#content');
      self.cleanup();
      self.off();
      this.$('#content_title').html("Detail Test");
      var view = new DetailTestView({uriStack: uriStack});
      self.setView('#content', view).render();
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
    },
    
    showDetail: function(uriStack, model){
      var self = this;
			// TODO: validate cleanup operations
      self.removeView('#content');
      self.cleanup();
      self.off();
      //
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
        this.$('#content_title').html(val);
      }
      
      titleFunc(model.resource.title + ': ' + 
        Iccbl.getTitleFromTitleAttribute(model,model.resource) );
      
      var view = new viewClass({ model: model, uriStack: uriStack});
      self.listenTo(view, 'update_title', titleFunc );
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
    
    }, // end showDetail
    
    showList: function(resource, uriStack, schemaResult) {
      
      var self = this;
			// TODO: validate cleanup operations
      self.removeView('#content');
      self.cleanup();
      self.off();
      
      var addNewResourceButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="add_button" href="#">',
          'Add</a>'
        ].join(''));
      addNewResourceButton.click(function(e){
        e.preventDefault();
        newUriStack = [resource.key,'+add'];
        var route = newUriStack.join('/');
        console.log('history route: ' + route);
        appModel.router.navigate(route, {trigger: true});
      });
      var extraControls = [];
      if (appModel.hasPermission(resource, 'edit')){
        extraControls.push(addNewResourceButton);
      }
      
      var uriStack = _.clone(uriStack);
      var viewClass = ListView;
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
        // search view
        console.log('create a collection with search data');
        // IF search, override collection fetch; POST search data
        // TODO: if search data contains a simple AND list, then just pass these
        // params on to the list and conduct the searches from the header field mechanism
        
        this.consumedStack.push(uriStack.shift());
        // TODO: searchID not used - will be needed to persist searches on server
        var searchID = uriStack.shift();

        this.consumedStack.push(searchID);
        this.$('#content_title').html(resource.title + ' search');
        
        var url = resource.apiUri + '/search/' + searchID;
        
        var search_data = appModel.getSearch(searchID);
        if(_.isUndefined(search_data) || _.isEmpty(search_data)){
          var msg = 'Content List search requires a "search_data:'+searchID +'" in current browser state';
          console.log(msg);
          appModel.error(msg);
          return;
        }
        var Collection = Iccbl.MyCollection.extend({
          url: url,
          fetch: function(options){
            // endcode for the post_data arg; post data elements are sent 
            // as urlencoded values via a POST form for simplicity
            console.log('execute POST for collection form data: ', search_data);
            options = options || {};
            var data = options.data || {};
            options.data = _.extend({ search_data: JSON.stringify(search_data) }, data);
            options.type = 'POST';
            Collection.__super__.fetch.call(this, options);
          }
        });
        var collection = new Collection({
          url: url,
        });

        var view = new viewClass({ model: appModel, 
          options: { 
            uriStack: uriStack,
            schemaResult: schemaResult, 
            resource: resource, 
            collection: collection, 
            search_data: search_data,
            extraControls: extraControls
          }
        });
      }else{ // normal list view
        var view = new viewClass({ model: appModel, 
          options: { 
            uriStack: uriStack,
            schemaResult: schemaResult, 
            resource: resource,
            extraControls: extraControls
          }
        });
      }
    
      self.listenTo(view, 'update_title', function(val){
        if(val) {
          this.$('#content_title').html(
              resource.title + 
              ': <small>' + val + '</small>');
        }else{
          this.$('#content_title').html(resource.title);
        }
      });
    
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      
      Backbone.Layout.setupView(view);
      self.setView('#content', view ).render();
      self.objects_to_destroy.push(view);

    }, // end showList
    
    /** Backbone.layoutmanager callback **/
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

    // delegated from UriContainerView on appModel change:uriStack
    changeUri: function(uriStack) {
      console.log('changeUri:', uriStack);
      var self = this;
      var consumedStack = this.consumedStack = [];

      if (!_.isEmpty(uriStack)){
        if(uriStack[0] == 'detail_test'){
          this.consumedStack.push(uriStack.shift());
          self.showDetailTest(uriStack);
          return;
        }
        var uiResourceId = uriStack.shift();
        this.consumedStack.push(uiResourceId);
      }else{
        uiResourceId = 'home';
      }
      
      var resource = appModel.getResource(uiResourceId);
      
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
      
      if (_.isUndefined(resource) || _.isUndefined(resource)){
        var msg = "Resource schema not defined: " + uiResourceId;
        appModel.error(msg)
        throw msg;
      }

      // Re-fetch the specific resource schema from the Resource endpoint:
      // - if the Resource has customizations of the schema
      // - i.e. "extraSelectorOptions
      appModel.getResourceFromUrl(uiResourceId, resource.apiUri + '/schema', function(resource){
        // Test for list args, if not found, then it's a detail view
        if (!_.isEmpty(uriStack) && 
              !_.contains(appModel.LIST_ARGS, uriStack[0]) &&
              uriStack[0] != 'search') {
          // DETAIL VIEW
          
          if(uriStack[0] == '+add'){
            self.showAdd(resource, uriStack);
          }else{ 
            var _key = Iccbl.popKeyFromStack(resource, uriStack, consumedStack );
            appModel.getModel(uiResourceId, _key, function(model){
              model.resource = resource;
              self.showDetail(uriStack, model);
            });
          }
        } else {
          // LIST VIEW
          self.showList(resource, uriStack,resource);
        }
        
      });

    }
  });

  return ContentView;
});