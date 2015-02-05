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
  'views/user/user2',
  'views/usergroup/usergroup2',
  'text!templates/content.html',
  'text!templates/welcome.html',
  'text!templates/about.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         EditView, LibraryView, UserAdminView, UserGroupAdminView, layout, welcomeLayout, 
         aboutLayout) {
  
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout, 
    'LibraryView': LibraryView,
    'UserAdminView': UserAdminView,
    'UserGroupAdminView': UserGroupAdminView
  };
    
  var ContentView = Iccbl.UriContainerView.extend({
    
    template: _.template(layout),
    className: "col-sm-10 col-md-10 col-lg-10",
    
    initialize: function() {
      console.log('initialize content.js');
      Iccbl.UriContainerView.prototype.initialize.apply(this,arguments);
      this.title = null;
      this.consumedStack = [];
    },
    
//    /**
//     * Child view bubble up URI stack change event
//     */
//    reportUriStack: function(reportedUriStack) {
//      var consumedStack = this.consumedStack || [];
//      var actualStack = consumedStack.concat(reportedUriStack);
//      appModel.set({'uriStack': actualStack}, {source: this});     
//    },
//    
//    /**
//     * Backbone.Model change event handler
//     * @param options.source = the event source triggering view
//     */
//    uriStackChange: function(model, val, options) {
//      if(options.source === this){
//        console.log('self generated uristack change');
//        return;
//      }else{
//        this.changeUri();
//      }
//    },
    
//    serialize: function(){
//      return {
//        'title': this.title
//      };
//    },
    
    showAdd: function(resource, uriStack){
      var self = this;

      // Note: have to fill up the default fields so that the edit form will
      // show those fields
      var defaults = {};
      _.each(_.keys(resource.schema.fields), function(key){
        
        var field = resource.schema.fields[key];
        if(_.isArray(field.visibility)  
            && _.contains(field.visibility, 'edit')){
          console.log('add field: ' + field.key );
        
          if (key == 'resource_uri') {
              defaults[key] = self.options.url;
          } else if (key == 'id'){
            // nop 
            // always exclude the id field to signal create case to the server
          } else {
            if(!_.isEmpty(field.default)){
              defaults[key] = field.default;
            }else{
              defaults[key] = ''; // fill the rest of the fields with blanks
            }
          }
        }
      });

      var NewModel = Backbone.Model.extend({
        urlRoot: resource.apiUri , defaults: defaults 
      });
      var newModel = new NewModel();
      newModel.resource = resource;
      
      this.$('#content_title').html(resource.title + ': Add' );
      var viewClass = DetailLayout;
      var view = new viewClass({ model: newModel, uriStack: uriStack});
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
    },
    
    showDetail: function(uriStack, model){
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
      
      this.$('#content_title').html(
          model.resource.title + 
          ': ' + Iccbl.getTitleFromTitleAttribute(model,model.resource.schema) + '' 
          );
      

      var view = new viewClass({ model: model, uriStack: uriStack});
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
    },
    
    
    showList: function(resource, uriStack, schemaResult) {
      var self = this;
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
//        var subresource = appModel.getResource(uriStack[1]);
        var substack = _.rest(uriStack,1);
        var _key = Iccbl.popKeyFromStack(resource, substack, []);
        this.$('#content_title').html('Child logs for: ' + _key);
      }else{
        this.$('#content_title').html(resource.title + ' listing');
      }
      
      view = new viewClass({ model: appModel, 
        options: { 
          uriStack: uriStack,
          schemaResult: schemaResult, 
          resource: resource
        }
      });
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
    },

    changeUri: function(uriStack) {
      var self = this;
      var consumedStack = this.consumedStack = [];

      if (!_.isEmpty(uriStack)){
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
      
      if (_.isUndefined(resource) || _.isUndefined(resource.schema)){
        var msg = "Resource schema not defined: " + uiResourceId;
        appModel.error(msg)
        throw msg;
      }

      // Test for list args, if not found, then it's a detail view
      if (!_.isEmpty(uriStack) && 
            !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        // DETAIL VIEW
        
        if(uriStack[0] == '+add'){
//        consumedStack = uriStack.unshift();
        self.showAdd(resource, uriStack);
        }else{ 
          var _key = Iccbl.popKeyFromStack(resource, uriStack, consumedStack );
          appModel.getModel(uiResourceId, _key, function(model){
            self.showDetail(uriStack, model);
          });
        }
      } else {
        // LIST VIEW
        self.showList(resource, uriStack,resource.schema);
      }
    }
  });

  return ContentView;
});