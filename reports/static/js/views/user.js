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
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, BackbonePageableCollection, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, SimpleListView, layout) {

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
        'title': Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema),
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
      event.preventDefault();
      // Block clicks from the wrong elements
      // TODO: how to make this specific to this view? (it is also catching
      //clicks on the table paginator)
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      this.change_to_tab(key);
    },

    change_to_tab: function(key){
      if(_.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active');
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

    setGroups: function(){
      var self = this;
      var key = 'groups';
      
      var view = this.tabViews[key];
      if ( !view ) {
        
        // TODO: make a userGroups subview
        var resourceId = 'usergroup';
        var resource = appModel.getResource(resourceId);

//        var createGroups = function(schemaResult){
//          var resource = schemaResult.resource;
          var columns = Iccbl.createBackgridColModel(
              resource.schema.fields, Iccbl.MyHeaderCell);

          // make a checkbox column that will be used to manage the group->user relationship
          columns.unshift({
              column: "isForUser",
              name: "isForUser",
              label: "Select groups for " + self.title,
              cell: "boolean"
          });

          var ListModel = Backbone.Model.extend({
              defaults: {
                  rpp: 25,
                  page: 1,
                  order: {},
                  search: {}}
              });
          var listModel = new ListModel({});
          
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

          // collection.bind("change reset add remove", function(){
          collection.bind("sync", function(){
            collection.each(function(model){
              var currentUsers = model.get('users');
              var currentUserUri = self.model.get('resource_uri');
              // Note; the resource_uri for the current user may be full, and the
              // uri's from the server may be relative (so not including "api/v1")
              // model.set('isForUser', _.find(currentUsers, function(user_uri) {
                // return  currentUserUri.indexOf(user_uri) != -1;
              // }));

              model.set('isForUser', Iccbl.containsByMatch(currentUsers, currentUserUri) );

              // when the checkbox is selected, use a custom model event to 
              // update the group model relationship and refetch the user model
              // reln is group->users, so update it that way
              model.on('change:isForUser', function(model){
                console.log('model change: ' + JSON.stringify(model));
                var currentUsers = model.get('users');
                var currentUserUri = self.model.get('resource_uri');
                var changed = false;
                var alreadyContained = Iccbl.containsByMatch(currentUsers, currentUserUri);
                if(model.get('isForUser') && ! alreadyContained) { 
                  currentUsers.unshift(currentUserUri);
                  changed = true;
                }else if(!model.get('isForUser') && alreadyContained ){
                  currentUsers = _.without(currentUsers, currentUserUri);
                  changed = true;
                }
                if(changed){
                  model.save({'users':currentUsers},{
                    patch: true,
                    success: function(){
                      self.model.fetch();
                    }
                  });
                }
              });
            });
          });
          

          var view = new SimpleListView(
              { model: listModel },
              { columns: columns,collection: collection, schemaResult: resource.schema });
          Backbone.Layout.setupView(view);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.tabViews[key] = view;
          self.setView("#tab_container", view ).render();
//        };
//  
//        appModel.getSchema(resourceId, createGroups);
      } else {
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
      }
    },

    setPermissions: function(){
      var self = this;
      var key = 'permissions';
      
      var view = self.tabViews[key];
      if ( !view ) {
        
        // TODO: make a userPermissions subview
        var resourceId = 'permission';
        var resource = appModel.getResource(resourceId);

//        // get all of the permissions first, then checkoff the ones that are for this user
//        var createPermissions = function(schemaResult){
//          var resource = schemaResult.resource;
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
            
//        };
//
//        appModel.getSchema(resourceId, createPermissions);
      } else {
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
      }
    },


//    initialize : function(attributes, options) {
//      this.options = options;
//      var self = this;
//      Iccbl.requireOptions(options, ['url_root', 'current_options', 
//                                     'user_model', 'user_schema', 'router']);
//
//      console.log('initialize UserView: current options: ' + JSON.stringify(options.current_options));
//
//      _.bindAll(this, 'setDetail', 'setGroups', 'setPermissions');
//
//      this._tabbed_resources = {
//          detail: { description: 'User Details', title: 'User Details', invoke : this.setDetail },
//          groups: { description: 'User Groups', title: 'User Groups', invoke : this.setGroups },
//          permissions: { description: 'User Permissions', title: 'User Permissions', invoke : this.setPermissions },
//      };
//
//      var compiledTemplate = _.template( tabbedTemplate, { tab_resources: this._tabbed_resources, });
//
//      this.$el.html(compiledTemplate);
//
//      var tab = options.current_options['tab'];
//      if(_.isUndefined(tab) || _.isEmpty(tab) ) tab = 'detail';
//
//      this.listenTo(this.model, 'change:current_detail', this.change_current_detail);
//
//      self.user = self.options.user_model;
//      self.userSchema = self.options.user_schema;
//
//      // set the id specifically on the model: backbone requires this to determine whether a "POST" or "PATCH" will be used
//      self.user.id = Iccbl.getIdFromIdAttribute(self.user, self.userSchema);
//      self.title = Iccbl.getTitleFromTitleAttribute(self.user, self.userSchema);
//
//      self._tabbed_resources[tab].invoke();
//      self.$('#' + tab).addClass('active');
//
//      console.log('initialized....');
//    },
//
//    change_current_detail: function(){
//        var _current_options = this.model.get('current_detail');
//        console.log('change_current_detail: ' + JSON.stringify(_current_options) );
//
//        //TODO: will there be a case where we have to reset all the options?
//        var tab = _current_options['tab'];
//        if(this.options.current_options['tab'] != tab)
//            if(!_.isUndefined(tab)) this.change_to_tab(tab);
//    },
//    render : function() {
//        window.scrollTo(0, 0);
//        return this;
//    }

    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });

  return UserView;
});