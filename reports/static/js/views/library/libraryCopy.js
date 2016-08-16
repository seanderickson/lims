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
  'templates/generic-tabbed.html'
], function($, _, Backbone, Iccbl, layoutmanager, appModel, 
            DetailLayout, EditView, ListView, LibraryCopyWellsView, 
            libraryCopyTemplate ) {
  
  var LibraryCopyView = Backbone.Layout.extend({
    
    template: _.template(libraryCopyTemplate),

    initialize: function(args) {
      var self = this;
      this.library = args.library;
      this.tabViews = {}; 
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      
      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail'){
          var permission = self.tabbed_resources[key].permission;
          if (_.isUndefined(permission)){
            permission = self.tabbed_resources[key].resource;
          }
          if (!appModel.hasPermission(permission)){
            delete self.tabbed_resources[key];
          }
        }
      });
      _.bindAll(this, 'click_tab');
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
      }
    },      
    
    events: {
        'click ul.nav-tabs >li': 'click_tab',
    },

    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      return {
        'base_url': [
           self.library.resource.key,self.library.key,self.model.resource.key,
           self.model.key].join('/'),
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if (viewId == '+add') {
          this.uriStack.unshift(viewId);
          this.showAdd();
          return;
        }else if (viewId == 'edit'){
          this.uriStack.unshift(viewId);
          this.showEdit();
          return;
        }

        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },
    
    click_tab : function(event){
      var self = this;
      event.preventDefault();
      event.stopPropagation();
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      appModel.requestPageChange({
        ok: function(){
          self.change_to_tab(key);
        }
      });
      
    },

    change_to_tab: function(key){
      if(_.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active');
        this.$('#' + key).addClass('active');
        if(key !== 'detail'){
          this.consumedStack = [key];
        }else{
          this.consumedStack = [];
        }
        var delegateStack = _.clone(this.uriStack);
        this.uriStack = [];
        var method = this[this.tabbed_resources[key]['invoke']];
        if (_.isFunction(method)) {
          method.apply(this,[delegateStack]);
        } else {
          throw "Tabbed resources refers to a non-function: " + this.tabbed_resources[key]['invoke']
        }
      }else{
        var msg = 'Unknown tab: ' + key;
        appModel.error(msg);
        throw msg;
      }
    },
    
    setDetail: function(delegateStack) {
      var key = 'detail';
      
      var view = this.tabViews[key];
      if ( !view ) {
        view = new DetailLayout({ 
          model: this.model,
          uriStack: delegateStack, 
          buttons: ['download', 'history'] });
        this.tabViews[key] = view;
      } 
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    showEdit: function() {
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack, 
        buttons: ['download', 'upload','history']
      });
      Backbone.Layout.setupView(view);

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },
    
    setCopyPlates: function(delegateStack){
      var self = this;
      var plateResource = appModel.getResource('librarycopyplate'); 

      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        // Detail view
        var plate_number = delegateStack.shift();
        self.consumedStack.push(plate_number);
        var _key = [self.model.get('name'),plate_number].join('/');
        appModel.getModel(plateResource.key, _key, function(model){
          var view = new DetailLayout({ 
            model: model,
            uriStack: delegateStack, 
            buttons: ['download','history'] });
          Backbone.Layout.setupView(view);
          self.listenTo(view ,'uriStack:change',self.reportUriStack);
          self.setView("#tab_container",view ).render();
        });        
        return;
      }else{
        // List view
        // special url because list is a child view for librarycopy
        var url = [self.model.resource.apiUri, 
                   self.model.key,
                   'plate'].join('/');
        view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: plateResource,
          resource: plateResource,
          url: url
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        this.setView("#tab_container", view ).render();
      }
    },
    
    setCopyPlatesOld: function(delegateStack){
      var self = this;
      var key = 'plate';
      var view = this.tabViews[key];
      if ( !view ) {
        var view = new LibraryCopyPlatesView({
          library: self.library,
          copy: self.model,
          uriStack: delegateStack
        });
        self.tabViews[key] = view;
        Backbone.Layout.setupView(view);
      } else {
        self.reportUriStack([]);
      }
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();      
    },

    setCopyWells: function(delegateStack) {
      var self = this;
      var key = 'copywell';
      var view = this.tabViews[key];
      if ( !view ) {
        // TODO: review the way this is done in screen.js with sub-tabs
        var view = new LibraryCopyWellsView({
          library: self.library,
          copy: self.model,
          uriStack: delegateStack
        });
        self.tabViews[key] = view;
        Backbone.Layout.setupView(view);
      } else {
        self.reportUriStack([]);
      }
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }
  
  });
  
  return LibraryCopyView;
});