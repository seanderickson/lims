define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/user/user_group_permissions',
    'views/generic_detail_layout',
    'views/list2',
    'templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, UserGroupPermissionView, DetailLayout, 
            ListView, layout) {

  var UserView = Backbone.Layout.extend({

    initialize: function(args) {
      var self = this;
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource);

      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail' && !appModel.hasPermission(
            self.tabbed_resources[key].resource,'read')){
          delete self.tabbed_resources[key];
        }
      });
      
      _.bindAll(this, 'click_tab');
    },
    
    template: _.template(layout),
    
    tabbed_resources: {
        detail: { 
          description: 'User Details', 
          title: 'User Details', 
          invoke: 'setDetail' },
//        usergroup: { 
//          description: 'User Groups', 
//          title: 'User Groups', 
//          invoke: 'setGroups', 
//          resource: 'usergroup' },
        usergrouppermissions: { 
          description: 'User Groups and Permissions', 
          title: 'User Groups and Permissions', 
          invoke: 'setGroupsPermissions'
        },
//        permission: { 
//          description: 'User Permissions', 
//          title: 'User Permissions', 
//          invoke: 'setPermissions',
//          resource: 'permission' },
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
          this.$('ul.nav-tabs > li').addClass('disabled');
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }
        if(viewId == 'edit'){
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          window.alert(msg);
          throw msg;
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
      
      console.log('++++hack remove disabled state');
      this.$('ul.nav-tabs > li').removeClass('disabled');
    },
    
    click_tab : function(event){
      var self = this;
      event.preventDefault();
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      if(this.$('#'+key).hasClass('disabled')){
        return;
      }
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
        if(key !== 'detail'){
          this.consumedStack = [key];
        }else{
          this.consumedStack = [];
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
      var key = 'detail';
 
      this.model.resource.fields['permissions']['choices'] = appModel.get('permissionOptions');
      
      console.log('setDetail: delegateStack: ', delegateStack);
        var view = new DetailLayout({ 
          model: this.model, 
          uriStack: delegateStack
        });
        this.tabViews[key] = view;
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // Note: since detail_layout reports the tab, the consumedStack is empty here
      this.consumedStack = []; 
      this.setView("#tab_container", view ).render();
      return view;
    },
    
    setGroupsPermissions(delegateStack){
      var self = this;

      var pUserModel = this.model.clone();
      pUserModel.key = this.model.key;
      var resource = _.extend({},this.model.resource);
      resource = _.extend({}, this.model.resource );
      resource.fields = _.pick(
        this.model.resource.fields,
        ['username','first_name','last_name','usergroups','permissions']);
      resource.fields['first_name']['visibility'] = [];
      resource.fields['last_name']['visibility'] = [];
      resource.fields['usergroups']['visibility'] = [];
      resource.fields['usergroups']['editability'] = ['u'];
      resource.fields['permissions']['editability'] = ['u'];
      resource.fields['first_name']['editability'] = [];
      resource.fields['last_name']['editability'] = [];
      pUserModel.resource = resource;
      pUserModel.resource.fields['permissions']['choices'] = (
          appModel.get('permissionOptions'));
      
      var view = new DetailLayout({ 
        model: pUserModel, 
        uriStack: delegateStack,
        DetailView: UserGroupPermissionView
      });
      
      Backbone.Layout.setupView(view);
      this.consumedStack = ['usergrouppermissions']; 
      //      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

//    setGroups: function(delegateStack){
//      var self = this;
//      var key = 'usergroup';
//      var resource = appModel.getResource('usergroup');
//      var url = [self.model.resource.apiUri, 
//                 self.model.key,
//                 'groups'].join('/');
//      var view = new ListView({ options: {
//        uriStack: _.clone(delegateStack),
//        schemaResult: resource,
//        resource: resource,
//        url: url
//      }});
//      Backbone.Layout.setupView(view);
//      this.consumedStack = [key]; 
//      self.reportUriStack([]);
//      self.listenTo(view , 'uriStack:change', self.reportUriStack);
//      this.setView("#tab_container", view ).render();
//    },    
    
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });

  return UserView;
});