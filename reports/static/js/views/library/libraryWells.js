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
      var reagents_url = library.resource.apiUri +'/' + library.key + '/well';
      var url = reagents_url;
      var resource = appModel.getResource('reagent');
      if (!_.isEmpty(uriStack) && !_.isEmpty(uriStack[0]) &&
              !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        var _key = Iccbl.popKeyFromStack(resource,uriStack,self.consumedStack);

        appModel.getModel('well', _key, function(model){
          model.resource = resource;
          self.showDetail(model);
        } );
      } else {
        self.consumedStack = [];
        self.showList(resource, reagents_url);
      }
        
    },    
    
    showDetail: function(model) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      // get the view class
      var viewClass = DetailLayout;

      // prevent normal display of structure image value
      var si_field = model.resource.schema.fields['structure_image']
      if(si_field){
        si_field['visibility'] = ['api'];
      }
      var view = new viewClass({ model: model, uriStack: uriStack });

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.$('#content_title').html("");
      self.setView('#content', view).render();
      
      // TODO: support for generic images
      if(model.has('structure_image')){
        self.$('#content').append(
            '<img style="position: absolute; top: 8em; right: 3em; height: 16em;" src="' 
            + model.get('structure_image') + '" alt="image" />')
      }
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

  return LibraryWellView;
});