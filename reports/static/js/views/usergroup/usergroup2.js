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
    'templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout, EditView,
            ListView, layout) {

  var UserView = Backbone.Layout.extend({

    initialize: function(args) {
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource);
      
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
      var self = this;
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
        if(viewId == 'edit'){
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
      
      this.model.resource.fields['permissions']['choices'] = appModel.get('permissionOptions');
      
      var editView = EditView.extend({});
      
      if ( !view ) {
        view = new DetailLayout({ 
          model: this.model, 
          uriStack: delegateStack,
          EditView: editView
        });
        view.showEdit = function(){
          appModel.initializeAdminMode(function(){
            var fields = self.model.resource.fields;
            var options = appModel.getUserGroupOptions();
            fields['super_groups']['choices'] = options;
            fields['sub_groups']['choices'] = options;
            fields['users']['choices'] = appModel.getUsernameOptions();
            fields['permissions']['choices'] = appModel.getPermissionsOptions();
            DetailLayout.prototype.showEdit.call(view,arguments);
          });  
        };
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
      view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url
      });
      Backbone.Layout.setupView(view);
      this.consumedStack = [key]; 
//      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    setSubgroups: function(delegateStack){
      var self = this;
      var key = 'subgroups';
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'subgroups'].join('/');
      view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: self.model.resource,
        resource: self.model.resource,
        url: url
      });
      Backbone.Layout.setupView(view);
      this.consumedStack = [key]; 
//      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      
    },
    
    setPermissionGroups: function(delegateStack){
      var self = this;
      var key = 'supergroups';
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'supergroups'].join('/');
      view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: self.model.resource,
        resource: self.model.resource,
        url: url
      }});
      Backbone.Layout.setupView(view);
      this.consumedStack = [key]; 
//      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    onClose: function() {
      this.tabViews = {};
      this.remove();
    }

  });

  return UserView;
});