define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_pageable',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_layout',
    'views/list',
    'views/simple-list',
    'views/user/groups',
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, BackbonePageableCollection, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, SimpleListView, GroupsView, layout) {

  // for compatibility with require.js, attach PageableCollection in the 
  // right place on the Backbone object
  // see https://github.com/wyuenho/backbone-pageable/issues/62
  Backbone.PageableCollection = BackbonePageableCollection;

  // TODO: create a genericTabbedLayout base class?
  var UserView = Backbone.Layout.extend({

    initialize: function(args) {
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema);
      this.model.id = Iccbl.getIdFromIdAttribute(this.model, this.model.resource.schema);
      
      _.bindAll(this, 'click_tab');
    },
    
    template: _.template(layout),
    
    tabbed_resources: {
        detail: { 
          description: 'User Details', 
          title: 'User Details', 
          invoke: 'setDetail' },
        groups: { 
          description: 'User Groups', 
          title: 'User Groups', 
          invoke: 'setGroups' },
        permissions: { 
          description: 'User Permissions', 
          title: 'User Permissions', 
          invoke: 'setPermissions' },
    },
    
    events: {
      'click li': 'click_tab',
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
      return {
//        'title': Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema),
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
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          window.alert(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
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
      
      if(appModel.isPagePending()){
        appModel.showModal(function(){
          self.change_to_tab(key);
        })
      }else{
        self.change_to_tab(key);
      }
    },

    change_to_tab: function(key){
      if(this.key && this.key === key){
        return;
      }else{
        this.key = key;
      }
      if(_.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active'); // TODO: use bootstrap tabs
        this.$('#' + key).addClass('active');
        
        this.consumedStack = [key];
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
        window.alert(msg);
        throw msg;
      }
    },
    
    setDetail: function(){
      var key = 'detail';
      
      var view = this.tabViews[key];
      if ( !view ) {
        view = new DetailLayout({ model: this.model});
        this.tabViews[key] = view;
      }
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // NOTE: since subview doesn't report stack, report it here
      //      this.reportUriStack([]);
      this.setView("#tab_container", view ).render();
    },

    setGroups: function(delegateStack){
      var self = this;
      var key = 'groups';
      var view = this.tabViews[key];
      if ( !view ) {      
        view = new GroupsView({ model: this.model, uriStack: delegateStack });
        self.tabViews[key] = view;
      }
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
//      self.reportUriStack([]);
    },
    
    setPermissions: function(){
      var self = this;
      var key = 'permissions';
      
      var view = self.tabViews[key];
      if ( !view ) {
        
        // TODO: make a userPermissions subview
        var resourceId = 'permission';
        var resource = appModel.getResource(resourceId);

        // get all of the permissions first, then checkoff the ones that are for this user
        var columns = Iccbl.createBackgridColModel(
            resource.schema.fields, Iccbl.MyHeaderCell);

        // make a checkbox column that will be used to manage the user->permission relationship
        columns.unshift({
            column: "isForUser",
            name: "isForUser",
            label: "Select permissions for " + self.title,
            cell: "boolean"
        });

        var ListModel = Backbone.Model.extend({
            defaults: {
                rpp: 25,
                page: 1,
                order: {},
                search: {}}
            });
        var listModel = new ListModel();

        var _state = {
            currentPage: parseInt(listModel.get('page')),
            pageSize: parseInt(listModel.get('rpp'))
          };

        var Collection = Iccbl.MyCollection.extend({
          state: _state
        });
        var collection  = new Collection({
          'url': resource.apiUri,
          listModel: listModel
        });

        // update the collection models with a new attribute tied to the checkbox column
        collection.bind("sync", function(){
          collection.each(function(model){
            var currentPerms = self.model.get('permissions');
            var currentPermissionUri = model.get('resource_uri');
            model.set('isForUser', Iccbl.containsByMatch(currentPerms, currentPermissionUri));

            // when the checkbox is selected, use a custom model event to 
            // update the group model relationship and refetch the user model
            model.on('change:isForUser', function(model){
                var currentPerms = self.model.get('permissions');
                var permissionUri = model.get('resource_uri');
                var containsByMatch = Iccbl.containsByMatch(currentPerms, permissionUri);
                var changed = false;
                if(model.get('isForUser') && ! containsByMatch ){
                  currentPerms.push(permissionUri);
                  changed = true;
                }else if(!model.get('isForUser')  && containsByMatch ) {
                  currentPerms = _.without(currentPerms, permissionUri);
                  changed = true
                }
                if(changed){
                  self.model.save({'permissions': currentPerms },{
                    patch: true,
                    success: function(){
                      console.log('refetch the permission model');
                      model.fetch();
                    }
                  });
                }
              });
            });
        });

        collection.fetch();

        var view = new SimpleListView(
            { model: listModel },
            { columns: columns,collection: collection, schemaResult: resource.schema});
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.tabViews[key] = view;
        self.setView("#tab_container", view ).render();
            
      } else {
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
      }
    },

    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });

  return UserView;
});