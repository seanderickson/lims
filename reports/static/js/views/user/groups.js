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
//      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema);
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
      this.currentGroups = _.clone(this.model.get('usergroups'));
      
      this.groupsCollection.each(function(model){
        if(model.hasChanged() && model.get('is_for_user') ) { 
          if(!Iccbl.containsByMatch(self.currentGroups, model.get('resource_uri'))){
            self.addGroups.add(model);
          }
          self.removeGroups.remove(model);
        }else if(model.hasChanged() && !model.get('is_for_user')){
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
        
        appModel.setPagePending(true);
        this.$('#collapseTwo').collapse('show');
      }else{
        appModel.setPagePending(false);
        this.$('#save').addClass('disabled');
        this.$('#cancel').addClass('disabled');
      }
    },
    
    save: function(){
      var self = this;
      
      if(!self.changed){
        appModel.error('nothing changed');
        return;
      }
      self.addGroups.each(function(model){
        self.currentGroups.push(model.get('resource_uri'));
      });
      self.removeGroups.each(function(model){
        self.currentGroups = _.without(self.currentGroups, model.get('resource_uri'));
      });
      
      var comment = this.$('#comment').val();
      
      if(_.isEmpty(comment)){
        this.$('#comment').parent().prepend('Required');
        this.$('#comment').parent().addClass('text-danger has-error');
        return;
      }
      
      self.model.save(
        {'usergroups':self.currentGroups},
        { patch: true,
          headers: {
            'APILOG_COMMENT': comment
          }, 
          success: function(){
            self.groupView.render();
            self.reset();
          }
        });
    },
    
    cancel: function(){
      this.reset();
    },
    
    reset: function(){
      var self=this;
      appModel.setPagePending(false);
      self.changed = false;
      self.addGroups.reset();
      self.removeGroups.reset();
      this.$('#add-selection-list').html("");
      this.$('#remove-selection-list').html("");
      this.$('#save').addClass('disabled');
      this.$('#cancel').addClass('disabled');
      this.$('#comment').parent().removeClass('text-danger has-error');
      this.$('#collapseTwo').collapse('hide');
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
      var resource = _.extend({}, appModel.getResource(resourceId));
      
      resource.schema.extraSelectorOptions = { 
          'label': 'For user', 'searchColumn': 'is_for_user', 
          'options': [ {value:'',label:'-'},
                       {value: self.model.id,label:'yes'},
                       {value:'false', label:'no'} ]}; //FIXME: need a not_eq operator
      
      var fields = resource.schema.fields;
      fields['name']['backgrid_cell_type'] = 'Iccbl.LinkCell';
      fields['name']['backgrid_cell_options'] = '/groups/{name}'
      var columns = Iccbl.createBackgridColModel(
          fields, Iccbl.MyHeaderCell);
      
      // make the (normally not visible) checkbox column
      var editable = appModel.getCurrentUser().is_superuser;
      columns.unshift({
          column: "is_for_user",
          name: "is_for_user",
          label: "Member",
          cell: "boolean",
          editable: editable
      });

      var collection  = this.groupsCollection = new Iccbl.MyCollection({
        'url': self.model.resource.apiUri + '/' + self.model.id + '/groups'
      });

      var uriStack = _.clone(this.uriStack);
      
      var view = this.groupView = new ListView({ options: {
        uriStack: uriStack,
        schemaResult: resource.schema,
        resource: resource,
        collection: collection,
        columns: columns
      }});
      Backbone.Layout.setupView(view);

      self.listenTo(self.groupsCollection, 'change', self.changeGroups );
      self.listenTo(view , 'uriStack:change', self.reportUriStack);

      this.setView("#list-with-selection-table", view).render();
      
    }    
    
  });
  
  return ListSelectionView;
});