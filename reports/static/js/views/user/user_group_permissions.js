define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'iccbl_backgrid',
    'models/app_state',
    'views/generic_detail_stickit',
    'views/generic_edit',
    'views/list2',
    'views/simple-list',
    'templates/user-groups-permissions-layout.html',

], function($, _, Backbone, stickit, Iccbl, appModel,
            DetailView, EditView, ListView, SimpleListView, layoutTemplate ) {

  var View = Backbone.Layout.extend({
    
    initialize: function(args) {
      var self = this;
      console.log('---- initialize user_group_permissions');
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.subviews = {};

      // NOTE: edit of groups/permissions is handled on the main user page
      var buttons = this.buttons = ['download','history','back'];
    },

    afterRender: function(){
      var delegateStack = _.clone(this.uriStack);
      this.showDetail(delegateStack);
    },
    
    showDetail: function(delegateStack){
      var self = this;
      
      // create the usergroup table
      var resource = appModel.getResource('usergroup');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'groups'].join('/');
      var ListModel = Backbone.Model.extend({
        defaults: {
            rpp: 50,
            page: 1,
            order: {},
            search: {}}
        });
      var listModel = new ListModel();
      var Collection = Iccbl.MyCollection.extend({
        url: url,
        listModel: listModel
      });
      var collection = self.collection = new Collection();
      var groupsView = new SimpleListView({
        uriStack: _.clone(delegateStack),
        collection: collection,
        resource: resource,
        includes: ['-all_users','-users','-permissions','-sub_groups','-super_groups']
      });
      Backbone.Layout.setupView(groupsView);
      this.setView("#group_content", groupsView ).render();
      
      var permissionsView = new DetailView({
        model: self.model,
        buttons: []
      });
      Backbone.Layout.setupView(permissionsView);
      this.setView("#permission_content", permissionsView ).render();
      self.$("#edit_content").hide();
      
    },
    
//    showEdit: function() {
//      console.log('showEdit...');
//      var self = this;
//
//      self.model.resource.fields['permissions']['choices'] = (
//          appModel.get('permissionOptions'));
//
//      var view = new EditView({ 
//        model: self.model, 
//        uriStack: [] 
//      });
//      Backbone.Layout.setupView(view);
//      self.listenTo(view,'remove',function(){
//        self.removeView(view);
//        self.$("#edit_content").hide();
//        self.$("#detail_content").show();
//      });
//      self.$("#detail_content").hide();
//      view = self.setView("#edit_content", view ).render();
//      self.$("#edit_content").show();
//      return view;
//    },    

    serialize: function() {
      return {
        'buttons': _.chain(this.buttons),
        'title': Iccbl.getTitleFromTitleAttribute(
            this.model,
            this.model.resource),
      }      
    },    
  
    template: _.template(layoutTemplate),
    
    onClose: function() {
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    }
    
  });
  
  return View;
});
