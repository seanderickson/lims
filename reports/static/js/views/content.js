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
  'text!templates/content.html',
  'text!templates/welcome.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         LibraryView, layout, welcomeLayout) {
  
  // put the view args in a keyed hash for lookup by name
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout, 
    'LibraryView': LibraryView
  };
    
  var ContentView = Iccbl.UriContainerView.extend({
    
    template: _.template(layout),
    
    initialize: function() {
      console.log('ContentView initialize');
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
      // get the view class
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

    changeUri: function(uriStack) {
    
      var self = this;

      console.log('got uristack: ' + JSON.stringify(uriStack));
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
          
        // DETAIL VIEW
        // The "scratch" variable holds the some of the prior state; for instance,
        // if "edit" view is called, it may contain the model to edit
//        var current_scratch = appModel.get('current_scratch');
//        appModel.set({ current_scratch: {} });
          
  
        var _key = Iccbl.popKeyFromStack(resource, uriStack, consumedStack );
        appModel.getModel(uiResourceId, _key, function(model){
          self.showDetail(resource, uriStack, model);
        });
        
//        if(_.isUndefined(current_scratch.schemaResult) || 
//            _.isUndefined(current_scratch.model)){  
//          appModel.getModel(uiResourceId, _key, self.showDetail );
//
//        }else{
//          appModel.updateModel(uiResourceId, _key, current_scratch.model, createDetail );
//        }
      } else {
        // LIST VIEW
        this.$('#title').html(resource.title);
        
        var createList = function(schemaResult) {
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
//              self.consumedStack = keysToReport;
              var newUriStack = [];
              appModel.updateModel(uiResourceId,key,model,function(model) {
                self.showDetail(resource, newUriStack, model);          
              });

//              appModel.updateModel(resourceId, key, model, self.showDetail);
//              var idList = Iccbl.getIdKeys(model,schemaResult);
//              // this is the standard action, will re-evaluate to the detail view
//              // NOTE: by setting the source to the view, signal to self to process
////              self.trigger('uriStack:change', idList, {source: view}); 
//              self.reportUriStack(idList, {source: view});
            });
            Backbone.Layout.setupView(view);
            self.setView('#content', view ).render();
        };
        appModel.getSchema(uiResourceId, createList);
        
      }
    }

  });

  return ContentView;
});