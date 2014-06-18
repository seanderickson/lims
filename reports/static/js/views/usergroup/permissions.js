define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_edit',
  'text!templates/list-with-selection-form.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, 
         EditForm, layout) {
  
  var ListSelectionView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    events: {
      'click button#save': 'save',
      'click button#cancel': 'cancel'
    },

    initialize: function(args) {

      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.model.id = Iccbl.getIdFromIdAttribute(this.model, this.model.resource.schema);
      
      var permissionResourceId = 'permission';
      var permissionResource = this.permissionResource = 
          _.extend({}, appModel.getResource(permissionResourceId)); 
      
      this.currentPermissions = _.clone(this.model.get('permissions'));
      
      this.addPermissions = new Backbone.Collection();
      this.removePermissions = new Backbone.Collection();
      
      self.changed = false;
      _.bindAll(this,'changePermissions');
    },
    
    changePermissions: function() {
      var self = this;
      this.currentPermissions = _.clone(this.model.get('permissions'));
      
      this.permissionCollection.each(function(model){
        if(model.hasChanged() && model.get('is_for_group') ) { 
          if(!Iccbl.containsByMatch(self.currentPermissions, model.get('resource_uri'))){
            self.addPermissions.add(model);
          }
          self.removePermissions.remove(model);
        }else if(model.hasChanged() && !model.get('is_for_group')){
          if(Iccbl.containsByMatch(self.currentPermissions, model.get('resource_uri'))){
            self.removePermissions.add(model);
          }
          self.addPermissions.remove(model);
        }        
      });
      
      if(!self.addPermissions.isEmpty() || !self.removePermissions.isEmpty() ){
        self.changed = true;
      }else{
        self.changed = false;
      }
      
      if(self.changed){
        this.$('#save').removeClass('disabled');
        this.$('#cancel').removeClass('disabled');
        
        this.$('#add-selection-list').html(self.addPermissions.reduce(function(memo,model){
          if(memo != "") memo += ", ";
          memo += Iccbl.getTitleFromTitleAttribute(model, self.permissionResource.schema);
          return memo;
        },""));
        this.$('#remove-selection-list').html(self.removePermissions.reduce(function(memo,model){
          if(memo != "") memo += ", ";
          memo += Iccbl.getTitleFromTitleAttribute(model, self.permissionResource.schema);
          return memo;
        },""));
        appModel.setPagePending(function(){
          // function to run if user requests another page then cancels that request.
          self.$('#collapseTwo').collapse('show');
          self.$('#collapseOne').collapse('hide');
        });
        this.$('#collapseTwo').collapse('show');
      }else{
        appModel.clearPagePending();
        this.$('#save').addClass('disabled');
        this.$('#cancel').addClass('disabled');
        self.reset();
      }
    },
    
    save: function(){
      var self = this;
      
      if(!self.changed){
        appModel.error('nothing changed');
        self.reset();
        return;
      }

      self.addPermissions.each(function(model){
        self.currentPermissions.push(model.get('resource_uri'));
      });
      self.removePermissions.each(function(model){
        self.currentPermissions = _.without(self.currentPermissons, model.get('resource_uri'));
      });
      
      var comment = this.$('#comment').val();
      
      if(_.isEmpty(comment)){
        this.$('#comment').parent().prepend('Required');
        this.$('#comment').parent().addClass('text-danger has-error');
        return;
      }
      
      self.model.save(
        {'permissions':self.currentPermissions},
        { patch: true,
          headers: {
            'APILOG_COMMENT': comment
          }, 
          success: function(){
            self.reset();
            self.permissionView.render();
          }
        });
    },
    
    cancel: function(){
      this.reset();
    },
    
    reset: function(){
      var self=this;
      appModel.clearPagePending();
      self.changed = false;
      self.addPermissions.reset();
      self.removePermissions.reset();
      this.$('#add-selection-list').html("");
      this.$('#remove-selection-list').html("");
      this.$('#save').addClass('disabled');
      this.$('#cancel').addClass('disabled');
      this.$('#comment').parent().removeClass('text-danger has-error');
      this.$('#comment').val("");
      this.$('#collapseOne').collapse('show');
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

      var resource = self.permissionResource;
      resource.schema.extraSelectorOptions = { 
          'label': 'For user', 'searchColumn': 'groups', 
          'options': [ {value: self.model.id,label:'yes'},
                       {value:'!'+self.model.id, label:'no'} ]}; //FIXME: need a not_eq operator
      
      var fields = resource.schema.fields;
//      // FIXME: is this wired?
//      fields['key']['backgrid_cell_type'] = 'Iccbl.LinkCell';
//      fields['key']['backgrid_cell_options'] = '/permissions/{scope}/{key}'; 
      var columns = Iccbl.createBackgridColModel(
          fields, Iccbl.MyHeaderCell);
      
      // make the (normally not visible) checkbox column
      var editable = appModel.getCurrentUser().is_superuser;
      columns.unshift({
          column: "is_for_group",
          name: "is_for_group",
          label: "Member",
          cell: "boolean",
          editable: editable
      });

      var collection  = this.permissionCollection = new Iccbl.MyCollection({
        'url': self.model.resource.apiUri + '/' + self.model.id + '/permissions'
      });

      var uriStack = _.clone(this.uriStack);
      
      var view = this.permissionView = new ListView({ options: {
        uriStack: uriStack,
        schemaResult: resource.schema,
        resource: resource,
        collection: collection,
        columns: columns
      }});
      Backbone.Layout.setupView(view);

      self.listenTo(self.permissionCollection, 'sync', function(){
        if(self.changed){
          self.permissionCollection.each(function(model){
            if(Iccbl.containsByMatch(
                self.addPermissions.pluck("resource_uri"), model.get('resource_uri'))){
              model.set({'is_for_group': true});
            }
            if(Iccbl.containsByMatch(
                self.removePermissions.pluck("resource_uri"), model.get('resource_uri'))){
              model.set({'is_for_group': false});
            }
          });
        }
      });

      // React to changes in the form:
      // NOTE: only listen -after- is editing, otherwise, change is called when syncing too
      self.listenTo(self.permissionCollection, 'backgrid:edit', function(){
        self.permissionCollection.once('change:is_for_group', self.changePermissions);
      })
      self.listenTo(view , 'uriStack:change', self.reportUriStack);

      this.setView("#list-with-selection-table", view).render();
    }    
    
  });
  
  return ListSelectionView;
});