define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_detail_layout',
  'views/library/library',
  'text!templates/genericResource.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         LibraryView, layout) {
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout
  };
    
  var LibraryWellView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      this.uriStack = args.uriStack;
      this.library = args.library;
      this.consumedStack = [];
      _.bindAll(this, 'showDetail');
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },

    // layoutmanager hook
    afterRender: function() {
      var self = this;
      
      var uriStack = this.uriStack;
      var library = this.library;

      // NOTE: library/well view is actually a "reagents" view
      // - use the well specific schema:
      var url = library.resource.apiUri +'/' + library.key + '/well';
      // but get the reagent list from the link:
      var reagents_url = library.resource.apiUri +'/' + library.key + '/reagent';
      
      var setupFunction = function(resource){
        console.log('setupFunction: ' + resource);
        // Test for list args, if not found, then it's a detail view
        if (!_.isEmpty(uriStack) && !_.isEmpty(uriStack[0]) &&
                !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
          var _key = Iccbl.popKeyFromStack(resource,uriStack,self.consumedStack);
  
          appModel.getModel(resource.key, _key, function(model){
            model.resource = resource;
            self.showDetail(model);
          } );
        } else {
          self.consumedStack = [];
          self.showList(resource, reagents_url);
        }
        
      };
      
      // get the library well specific resource
      // TODO: the well view should also get it itself so that the well view looks correct?
      // -or- the libary-specific well view has it's own schema and url: i.e.:
      // library/short_name/well/schema
      // library/short_name/well/well_id
      var wellSpecificResource;
      var wellResourceId = library.key + '-wells';
      try{
        // try to get cached resource
        setupFunction(appModel.getResource(wellResourceId));
      } catch(e) {
        console.log('get resource from server ' + appModel);
        appModel.getResourceFromUrl(wellResourceId, url+'/schema', setupFunction );
      }
      
//      var ModelClass = Backbone.Model.extend({
//        url : url + 'schema',
//        defaults : {}
//      });
//      var instance = new ModelClass();
//      instance.fetch({
//          success : function(model) {
//            
//            var resourceId = 'well';
//            var resource = appModel.getResource(resourceId);
//            
//            resource.schema = model.toJSON();
//            resource.schema = _.extend(resource.schema, appModel.schemaClass);
//            
//            // Test for list args, if not found, then it's a detail view
//            if (!_.isEmpty(uriStack) && !_.isEmpty(uriStack[0]) &&
//                    !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
//              var _key = Iccbl.popKeyFromStack(resource,uriStack,self.consumedStack);
//      
//              appModel.getModel(resourceId, _key, function(model){
//                model.resource = resource;
//                self.showDetail(model);
//              } );
//            } else {
//              self.consumedStack = [];
//              self.showList(resource, url);
//            }
//          },
//          error: appModel.jqXHRerror
//      });      
          
    },    
    
    showDetail: function(model) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      // get the view class
      var viewClass = DetailLayout;
      var view = new viewClass({ model: model, uriStack: uriStack });

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#content', view).render();
      
      
      self.$('#content').append(
          '<img style="position: absolute; top: 8em; right: 3em; height: 16em;" src="/db/well_image/' 
          + model.get('well_id') + '" alt="image" />')
    },
    
    showList: function(resource, url) {

      var self = this;
      var uriStack = _.clone(this.uriStack);
      var detailHandler = function(model) {
        var key = Iccbl.getIdFromIdAttribute(model,resource.schema);
        model.resource = resource;
        model.key = key;
        var keysToReport = Iccbl.getIdKeys(model,resource.schema);
        self.consumedStack = keysToReport;
        self.showDetail(model);
      };
      var view = new ListView({ options: {
        uriStack: uriStack,
        schemaResult: resource.schema,
        resource: resource,
        url: url, 
        detailHandler: detailHandler
      }});
    
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#content', view ).render();
    }
  });

  return LibraryWellView;
});