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
  'text!templates/content.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         LibraryView, layout) {
  
  // put the view args in a keyed hash for lookup by name
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout, 
    'LibraryView': LibraryView
  };
    
  var ContentView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function() {
      console.log('ContentView initialize');
      this.listenTo(appModel, 'change:uriStack', this.uriStackChange );
    },
    
    reportUriStack: function(reportedUriStack) {
      var actualStack = this.consumedStack.concat(reportedUriStack);
      appModel.set({'uriStack': actualStack}, {source: this});     
    },
    
    uriStackChange: function(model, val, options) {
      if(options.source === this){
        console.log('self generated uristack change');
        return;
      }else{
        this.changeUri();
      }
    },
    
    changeUri: function() {
    
      var self = this;

      var uriStack = _.clone(appModel.get('uriStack'));
      console.log('got uristack: ' + JSON.stringify(uriStack));
      var consumedStack = this.consumedStack = [];
      
      if (_.isEmpty(uriStack)){
        // TODO: go to home
        return;
      }
      
      // Test for a view ID
      var uiResourceId = uriStack.shift();
      this.consumedStack.push(uiResourceId);
      var resource = appModel.getResource(uiResourceId);
      
      if (!_.isEmpty(uriStack) && 
          ( !_.contains(appModel.LIST_ARGS, uriStack[uriStack.length-1]) || 
            uriStack[uriStack.length-1] === 'key' ) ) {
          
        // DETAIL VIEW
        // The "scratch" variable holds the some of the prior state; for instance,
        // if "edit" view is called, it may contain the model to edit
        var current_scratch = appModel.get('current_scratch');
        appModel.set({ current_scratch: {} });
          
        // get the view class
        var viewClass = DetailLayout;
        if (_.has(resource, 'detailView')){
          if (_.has(VIEWS, resource['detailView'])) {
            viewClass = VIEWS[resource['detailView']];
            console.log('found view ' + resource['detailView']);
          }else{
            window.alert('detailView class specified could not be found: ' + 
                resource['detailView']);
          }
        }
        var createDetail = function(model){
          var view = new viewClass({ model: model, uriStack: uriStack });
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView('#content', view).render();
        }
  
        var _key = Iccbl.popKeyFromStack(resource, uriStack, consumedStack );
        if(_.isUndefined(current_scratch.schemaResult) || 
            _.isUndefined(current_scratch.model)){  
          appModel.getModel(uiResourceId, _key, createDetail );

        }else{
          appModel.updateModel(uiResourceId, _key, current_scratch.model, createDetail );
        }
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
            Backbone.Layout.setupView(view);
            self.setView('#content', view ).render();
        };
        appModel.getSchema(uiResourceId, createList);
        
      }
    }

  });

  return ContentView;
});