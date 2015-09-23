define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_layout',
    'views/list2',
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, layout) {

  var UserView = Backbone.Layout.extend({

    initialize: function(args) {
      var self = this;
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema);

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
        usergroup: { 
          description: 'User Groups', 
          title: 'User Groups', 
          invoke: 'setGroups', 
          resource: 'usergroup' },
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
      return {
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
      var key = 'detail';
 
      this.model.resource.schema.fields['permissions']['choices'] = appModel.get('permissionOptions');
      
//      var onEditCallBack = function(displayFunction){
//        appModel.getUserGroupOptions(function(options){
//          self.model.resource.schema.fields['super_groups']['choices'] = options;
//          self.model.resource.schema.fields['sub_groups']['choices'] = options;
//          appModel.getUserOptions(function(options){
//            self.model.resource.schema.fields['users']['choices'] = options;
//            
//            displayFunction();
//          });
//        });
//      };
      
      var view = this.tabViews[key];
      if ( !view ) {
        view = new DetailLayout({ 
          model: this.model, 
          uriStack: delegateStack
//          onEditCallBack: onEditCallBack 
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

    setGroups: function(delegateStack){
      var self = this;
      var key = 'usergroup';
      var resource = appModel.getResource('usergroup');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'groups'].join('/');
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
    
//    setGroups: function(delegateStack){
//      var self = this;
//      var key = 'usergroup';
//      var view = this.tabViews[key];
//      if ( !view ) {      
//        view = new GroupsView({ model: this.model, uriStack: delegateStack });
//        self.tabViews[key] = view;
//      }
//      this.consumedStack = [key]; 
//      self.listenTo(view , 'uriStack:change', self.reportUriStack);
//      self.setView("#tab_container", view ).render();
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