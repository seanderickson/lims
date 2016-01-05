define([
    'jquery',
    'underscore',
    'backbone',
//    'backbone_pageable',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_layout',
    'views/list2',
//    'views/simple-list',
//    'views/usergroup/users',
//    'views/usergroup/permissions',
//    'views/usergroup/super_groups',
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, layout) {

  // TODO: create a genericTabbedLayout base class?
  var UserView = Backbone.Layout.extend({

    initialize: function(args) {
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema);
      
      _.bindAll(this, 'click_tab');
    },
    
    template: _.template(layout),
    
    tabbed_resources: {
        detail: { 
          description: 'Group Details', 
          title: 'Group Details', 
          invoke: 'setDetail' },
        users: { 
          description: 'Group Users', 
          title: 'Group Users', 
          invoke: 'setUsers' },
        supergroups: { 
          description: 'Groups that this group inherits permissions from', 
          title: 'Permission Super Groups', 
          invoke: 'setPermissionGroups' },
        subgroups: { 
          description: 'Groups that this group inherits users from', 
          title: 'User Subgroups', 
          invoke: 'setSubgroups' },
    },
    
    events: {
      'click ul.nav-tabs >li': 'click_tab',
    },

    /**
     * Layoutmanager hook
     */
    serialize: function() {
      return {
        'base_url': self.model.resource.key + '/' + self.model.key,
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
        if(viewId == '+add'){
          this.uriStack.unshift(viewId);
          viewId = 'detail';
        }
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          console.log('Error: ' + msg);
          appModel.error(msg);
          return;
        }
      }
      this.change_to_tab(viewId);
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    click_tab : function(event){
      var self = this;
      event.preventDefault();
      // Block clicks from the wrong elements
      // TODO: how to make this specific to this view? (it is also catching
      //clicks on the table paginator)
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      
      if(this.key && this.key === key){
        return;
      }
      
      appModel.requestPageChange({
        ok: function(){
          self.change_to_tab(key);
        }
      });
    },

    change_to_tab: function(key){

      if(_.has(this.tabbed_resources, key)){
        var delegateStack = _.clone(this.uriStack);
//        delegateStack.push(key);
        if(!_.isUndefined(this.tabbed_resources[key].alias)){
          key = this.tabbed_resources[key].alias;
        }
        if(this.key && this.key === key){
          return;
        }else{
          this.key = key;
        }        
        
        this.$('li').removeClass('active'); // TODO: use bootstrap tabs
        this.$('#' + key).addClass('active');
        
        this.uriStack = [];
        var method = this[this.tabbed_resources[key]['invoke']];
        if (_.isFunction(method)) {
          method.apply(this,[delegateStack]);
        } else {
          throw "Tabbed resources refers to a non-function: " + this.tabbed_resources[key]['invoke']
        }
      }else{
        var msg = 'Unknown tab: ' + key;
        window.alert(msg);
        throw msg;
      }
    },
    
    setDetail: function(delegateStack){
      var self = this;
      var key = 'detail';
      
      var view = this.tabViews[key];
      
      this.model.resource.schema.fields['permissions']['choices'] = appModel.get('permissionOptions');
      
      var onEditCallBack = function(displayFunction){
        appModel.getUserGroupOptions(function(options){
          self.model.resource.schema.fields['super_groups']['choices'] = options;
          self.model.resource.schema.fields['sub_groups']['choices'] = options;
          appModel.getUserOptions(function(options){
            self.model.resource.schema.fields['users']['choices'] = options;
            
            displayFunction();
          });
        });
      };
      
      if ( !view ) {
        view = new DetailLayout({ 
          model: this.model, 
          uriStack: delegateStack,
          onEditCallBack: onEditCallBack 
        });
        this.tabViews[key] = view;
      }
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // Note: since detail_layout reports the tab, the consumedStack is empty here
      this.consumedStack = []; 
      this.setView("#tab_container", view ).render();
      return view;
    },
    
    setUsers: function(delegateStack){
      var self = this;
      var key = 'users';
      var resource = appModel.getResource('user');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'users'].join('/');
      view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource.schema,
        resource: resource,
        url: url
      }});
      Backbone.Layout.setupView(view);
      this.consumedStack = [key]; 
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    //    setUsers: function(delegateStack){
    //      var self = this;
    //      var resource = appModel.getResource('user');
    //      var key = resource.key;
    //      var view = this.tabViews[key];
    //      if ( !view ) {      
    //        view = new UsersView({ model: this.model, uriStack: delegateStack });
    //        self.tabViews[key] = view;
    //      }
    //      this.consumedStack = [key]; 
    //      self.listenTo(view , 'uriStack:change', self.reportUriStack);
    //      self.setView("#tab_container", view ).render();
    //    },

    setSubgroups: function(delegateStack){
      var self = this;
      var key = 'subgroups';
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'subgroups'].join('/');
      view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: self.model.resource.schema,
        resource: self.model.resource,
        url: url
      }});
      Backbone.Layout.setupView(view);
      this.consumedStack = [key]; 
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
//      }
      
    },
    
//    setPermissionGroups: function(delegateStack){
//      var self = this;
//      var key = 'supergroups';
//      var url = [self.model.resource.apiUri, 
//                 self.model.key,
//                 'supergroups'].join('/');
//      view = new ListView({ options: {
//        uriStack: _.clone(delegateStack),
//        schemaResult: self.model.resource.schema,
//        resource: self.model.resource,
//        url: url
//      }});
//      Backbone.Layout.setupView(view);
//      this.consumedStack = [key]; 
//      self.reportUriStack([]);
//      self.listenTo(view , 'uriStack:change', self.reportUriStack);
//      this.setView("#tab_container", view ).render();
//    },
//    setPermissionGroups: function(delegateStack){
//      var self = this;
//      var key = 'supergroups';
//      var view = this.tabViews[key];
//      if ( !view ) {      
//        view = new SuperGroupsView({ model: this.model, uriStack: delegateStack });
//        self.tabViews[key] = view;
//      }
//      this.consumedStack = [key]; 
//      self.listenTo(view , 'uriStack:change', self.reportUriStack);
//      self.setView("#tab_container", view ).render();
//      
//    },
    
    //    setPermissions: function(delegateStack){
    //      var self = this;
    //      var key = 'permission';
    //      var view = this.tabViews[key];
    //      if ( !view ) {      
    //        view = new PermissionsView({ model: this.model, uriStack: delegateStack });
    //        self.tabViews[key] = view;
    //      }
    //      this.consumedStack = [key]; 
    //      self.listenTo(view , 'uriStack:change', self.reportUriStack);
    //      self.setView("#tab_container", view ).render();
    //    },
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });

  return UserView;
});