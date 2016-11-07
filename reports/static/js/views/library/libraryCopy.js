define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/generic_detail_layout',
  'views/generic_edit',
  'views/list2',
  'views/library/libraryCopyWells',
  'views/library/libraryCopyPlate',
  'utils/tabbedController'
], function($, _, Backbone, Iccbl, layoutmanager, appModel, 
            DetailLayout, EditView, ListView, LibraryCopyWellsView, 
            LibraryCopyPlateView, TabbedController ) {
  
  var LibraryCopyView = TabbedController.extend({

    initialize: function(args) {
      var self = this;
      this._classname = 'libraryCopy';

      this.library = args.library;
      TabbedController.prototype.initialize.apply(this,arguments);
    },

    tabbed_resources: {
      detail: { 
        description: 'Copy Details', 
        title: 'Copy Details', 
        invoke: 'setDetail'
      },
      plate: { 
        description: 'Copy Plates', title: 'Copy Plates', invoke: 'setCopyPlates',
        resource: 'librarycopyplates'
      },
      copywell: { 
        description: 'Copy Wells', title: 'Copy Wells', invoke: 'setCopyWells',
        resource: 'copywell'
      },
      platelocation: { 
        description: 'Copy Plate Locations', title: 'Copy Plate Locations', invoke: 'setCopyPlateLocations',
        resource: 'platelocation'
      }
    },      
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      return {
        'base_url': [
           self.library.resource.key,self.library.key,'copy',
           self.model.key].join('/'),
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    setDetail: function(delegateStack) {
      var key = 'detail';
      
      var view = new DetailLayout({ 
        model: this.model,
        uriStack: delegateStack
      });
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    setCopyPlateLocations: function(delegateStack){
      var self = this;
      var key = 'platelocation';
      var resource = appModel.getResource(key);
      var url = [resource.apiUri, 
                 '?library_copy=' + 
                 self.library.get('short_name') + '/' 
                 + self.model.get('copy_name')].join('/')
      // use a child container to manage batch edit functionality
      var view = new ListView({
        uriStack: delegateStack,
        url: url,
        resource: resource
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      
    },
    
    setCopyPlates: function(delegateStack) {
      var self = this;
      // use a child container to manage batch edit functionality
      var view = new LibraryCopyPlateView({
        library: self.library,
        copy: self.model,
        uriStack: delegateStack,
        resource: appModel.getResource('librarycopyplate')
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    setCopyWells: function(delegateStack) {
      var self = this;
      var key = 'copywell';
      // TODO: review the way this is done in screen.js with sub-tabs
      var view = new LibraryCopyWellsView({
        library: self.library,
        copy: self.model,
        uriStack: delegateStack
      });
      
      view.on('showCopy', function(){
        self.change_to_tab('detail');
      });
      view.on('showPlate', function(plate_number){
        self.uriStack = [plate_number];
        self.change_to_tab('plate');
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    }

  });
  
  return LibraryCopyView;
});