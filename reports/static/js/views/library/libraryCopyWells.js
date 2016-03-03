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
    
  var LibraryCopyWellView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      this.uriStack = args.uriStack;
      this.library = args.library;
      this.copy = args.copy;
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
    afterRender: function(){
      var uriStack = this.uriStack;
      var library = this.library;
      var copy = this.copy;
      var copykey = this.copy.key.split('/')[1]; // in this context, only need the copy name

      var url = library.resource.apiUri +'/' + library.key + '/copy/' + copykey + '/copywell';
      var resourceId = 'copywell';
      var resource = appModel.getResource(resourceId);

      // Test for list args, if not found, then it's a detail view
      if (!_.isEmpty(uriStack) && !_.isEmpty(uriStack[0]) &&
              !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {

        // Determine the key from the stack:
        // Note: because this is in the context of a library/copy, the key for the 
        // librarycopywell can be simpler, ignoring the library__short_name, copy_name
        var _key;
        var stackItem = uriStack.shift();
        if(stackItem !== copy.key) {
          // assume that it is the well_id
        }else{
          stackItem = uriStack.shift();
        }
        this.consumedStack = [stackItem];
        
        _key = library.key + '/' + copy.key + '/' + stackItem;

        appModel.getModel(resourceId, _key, this.showDetail );
      } else {
        this.consumedStack = [];
        this.showList(resource, url);
      }
    },     
    
    showDetail: function(model) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      // get the view class
      var viewClass = DetailLayout;

      var view = new viewClass({ model: model, uriStack: uriStack });

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.$('#content_title').html("");
      self.setView('#content', view).render();
      
    },
    
    showList: function(resource, url) {

      var self = this;
      var uriStack = _.clone(this.uriStack);

      var view = new ListView({ options: {
        uriStack: uriStack,
        schemaResult: resource,
        resource: resource,
        url: url
      }});
      self.listenTo(view, 'update_title', function(val){
        if(val) {
          this.$('#content_title').html('<small>' + val + '</small>');
        }else{
          this.$('#content_title').html("");
        }
      });
    
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#content', view ).render();
    }
  });

  return LibraryCopyWellView;
});