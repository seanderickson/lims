define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_edit',
  'text!templates/list-with-selection-form.html',
  'text!templates/modal_ok_cancel.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, 
         EditForm, layout, modalOkCancel) {
  
  var ListSelectionView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    events: {
      'click button#save': 'save',
      'click button#cancel': 'cancel'
    },
    
    initialize: function(args) {

      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema);
      this.model.id = Iccbl.getIdFromIdAttribute(this.model, this.model.resource.schema);
      
      var groupResourceId = 'usergroup';
      var groupResource = this.groupResource = appModel.getResource(groupResourceId);      
      
      this.currentGroups = _.clone(this.model.get('usergroups'));
      
      this.addGroups = new Backbone.Collection();
      this.removeGroups = new Backbone.Collection();
      
      // dirty flag
      self.changed = false;
    },
    
    changeGroups: function() {
      var self = this;
      
      this.groupsCollection.each(function(model){
        if(model.hasChanged() && model.get('isForUser') ) { 
          if(!Iccbl.containsByMatch(self.currentGroups, model.get('resource_uri'))){
            self.addGroups.add(model);
          }
          self.removeGroups.remove(model);
        }else if(model.hasChanged() && !model.get('isForUser')){
          if(Iccbl.containsByMatch(self.currentGroups, model.get('resource_uri'))){
            self.removeGroups.add(model);
          }
          self.addGroups.remove(model);
        }        
      });
      
      if(!self.addGroups.isEmpty() || !self.removeGroups.isEmpty() ){
        self.changed = true;
      }else{
        self.changed = false;
      }
      if(self.changed){
        this.$('#save').removeClass('disabled');
        this.$('#cancel').removeClass('disabled');
        this.$('#add-selection-list').html(self.addGroups.reduce(function(memo,model){
          if(memo != "") memo += ", ";
          memo += Iccbl.getTitleFromTitleAttribute(model, self.groupResource.schema);
          return memo;
        },""));
        this.$('#remove-selection-list').html(self.removeGroups.reduce(function(memo,model){
          if(memo != "") memo += ", ";
          memo += Iccbl.getTitleFromTitleAttribute(model, self.groupResource.schema);
          return memo;
        },""));
        if(self.changed){
          console.log('changed');
          
          appModel.setPagePending(true);

          this.$('#collapseTwo').collapse('show');
          //          this.$('#collapseOne').collapse('show');
        }
      }else{
        appModel.setPagePending(false);

        this.$('#save').addClass('disabled');
        this.$('#cancel').addClass('disabled');
      }
    },
    
    save: function(){
      var self = this;
      self.addGroups.each(function(model){
        self.currentGroups.push(model.get('resource_uri'));
      });
      self.removeGroups.each(function(model){
        self.currentGroups = _.without(self.currentGroups, model.get('resource_uri'));
      });
      
      self.model.save({'usergroups':self.currentGroups},{patch: true});
      appModel.setPagePending(false);
      this.reset();
    },
    
    cancel: function(){
      appModel.setPagePending(false);
      this.reset();
    },
    
    reset: function(){
      var self=this;
      self.addGroups.reset();
      self.removeGroups.reset();
      this.$('#add-selection-list').html("");
      this.$('#remove-selection-list').html("");
      this.$('#save').addClass('disabled');
      this.$('#cancel').addClass('disabled');
      this.groupsCollection.fetch();
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
    afterRender: function(){

      var self = this;

      var resourceId = 'usergroup';
      var resource = appModel.getResource(resourceId);

      var columns = Iccbl.createBackgridColModel(
          resource.schema.fields, Iccbl.MyHeaderCell);

      // make a checkbox column that will be used to manage the group->user relationship
      columns.unshift({
          column: "isForUser",
          name: "isForUser",
          label: "",
          cell: "boolean"
      });

      var collection  = this.groupsCollection = new Iccbl.MyCollection({
        'url': resource.apiUri
      });

      // collection.bind("change reset add remove", function(){
      collection.fetch({ success: function(){
        collection.each(function(model){
          var currentUsers = model.get('users');
          var currentUserUri = self.model.get('resource_uri');
          // Note; the resource_uri for the current user may be full, and the
          // uri's from the server may be relative (so not including "api/v1")
          model.set('isForUser', Iccbl.containsByMatch(currentUsers, currentUserUri));
        });
      }});
      self.listenTo(self.groupsCollection, 'change', self.changeGroups );
      
      var uriStack = _.clone(this.uriStack);
      
      var view = new ListView({ options: {
        uriStack: uriStack,
        schemaResult: resource.schema,
        resource: resource,
        url: resource.apiUri,
        collection: collection,
        columns: columns
      }});
      Backbone.Layout.setupView(view);
      
      this.setView("#list-with-selection-table", view).render();
      
    },    
    cleanup: function(event){
      console.log('cleanup');
    },
    
    remove: function(args){
      $(this).remove();
    }
    
    
  });
  
  return ListSelectionView;
});