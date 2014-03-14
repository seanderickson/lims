define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_detail_layout',
  'views/library',
  'views/user',
  'text!templates/content.html',
  'text!templates/welcome.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         LibraryView, UserAdminView, layout, welcomeLayout) {
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout, 
    'LibraryView': LibraryView,
    'UserAdminView': UserAdminView
  };
    
  var ContentView = Iccbl.UriContainerView.extend({
    
    template: _.template(layout),
    
    initialize: function() {
      console.log('initialize content.js');
      Iccbl.UriContainerView.prototype.initialize.apply(this,arguments);
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
    
    showDetail: function(resource, uriStack, model){
      var self = this;
      var uriStack = _.clone(uriStack);
      var viewClass = DetailLayout;
      if (_.has(resource, 'detailView')){
        if (_.has(VIEWS, resource['detailView'])) {
          viewClass = VIEWS[resource['detailView']];
          console.log('found view ' + resource['detailView']);
        }else{
          var msg = 'detailView class specified could not be found: ' + 
              resource['detailView'];
          
          window.alert(msg);
          throw msg;
        }
      }
      var view = new viewClass({ model: model, uriStack: uriStack });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
    },
    
    showList: function(resource, uriStack, schemaResult) {
      var self = this;
      var uriStack = _.clone(uriStack);
      var view = new ListView({ model: appModel, 
        options: { 
          uriStack: uriStack,
          schemaResult: schemaResult, 
          resource: resource
        }
      });

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.listenTo(view, 'detail', function(model){
        var key = Iccbl.getIdFromIdAttribute(model,schemaResult);
        var keysToReport = Iccbl.getIdKeys(model,schemaResult);
        
        // NOTE: TODO: have to call this explicitly because this is the
        // top level of the view hierarchy.  Could/would delegate, but then
        // need to have the child view consume the keys from the stack
        self.reportUriStack(keysToReport, {source: view});
        var newUriStack = [];
        model.resource = resource;
        model.key = key;
//        appModel.updateModel(resource.key,key,model,function(model) {
          self.showDetail(resource, newUriStack, model);          
//        });
      });
      Backbone.Layout.setupView(view);
      self.setView('#content', view ).render();
    },

    changeUri: function(uriStack) {
      var self = this;
      var consumedStack = this.consumedStack = [];
      
      if (_.isEmpty(uriStack)){
        var WelcomeView = Backbone.Layout.extend({
          template: _.template(welcomeLayout)
        });
        var view = new WelcomeView();
        self.setView('#content', view).render();
        return;
      }
      
      // Test for a view ID
      var uiResourceId = uriStack.shift();
      this.consumedStack.push(uiResourceId);
      var resource = appModel.getResource(uiResourceId);

      // Test for list args, if not found, then it's a detail view
      if (!_.isEmpty(uriStack) && 
            !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        var _key = Iccbl.popKeyFromStack(resource, uriStack, consumedStack );
        appModel.getModel(uiResourceId, _key, function(model){
          self.showDetail(resource, uriStack, model);
        });
      } else {
        
        self.showList(
            resource, uriStack, 
            resource.schema);
//        appModel.getSchema(uiResourceId, function(schemaResult) {
//          self.showList(resource, uriStack, schemaResult);
//        });
      }
    }
    
  });

  return ContentView;
});