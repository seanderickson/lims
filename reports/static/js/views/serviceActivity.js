define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_detail_stickit',
  'views/generic_edit',
  'views/generic_detail_layout'
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, 
         DetailView, EditView, DetailLayout) {

  var ServiceActivityView = DetailLayout.extend({
    
    initialize: function(args) {
      
      console.log('serviceActivityView', args);
      var self = this;
      self._args = args;
      self.screen = args.screen;
      self.user = args.user;
      
      var fields = self.model.resource.fields;
      
      var detailView = DetailView.extend( {
        
        afterRender: function(){
          
          DetailView.prototype.afterRender.apply(this,arguments);

        }
      }, self._args);
      
      var editView = EditView.extend({
        
        save_success: function(data, textStatus, jqXHR){
          var meta = _.result(data, 'meta', null);
          if (meta) {
            appModel.showJsonMessages(meta);
          }
          var urlPath = [];
          if (self.user) {
            urlPath = [self.user.resource.key,self.user.key,'serviceactivity'];
          } else if (self.screen){
            urlPath = [self.screen.resource.key,self.screen.key,'activities'];
          }
          if (! this.model.isNew()){
            urlPath.push(Iccbl.getIdFromIdAttribute( this.model,this.model.resource ));
          }           
          appModel.clearPagePending();
          appModel.router.navigate(urlPath.join('/'),{trigger:true});

        },
        
        afterRender: function(){

          EditView.prototype.afterRender.apply(this,arguments);

        }//editView.afterRender
      }, self._args);//editView
      args.EditView = editView;
      args.DetailView = detailView;
      
      DetailLayout.prototype.initialize.call(this,args);

    },
    
    afterRender: function() {
      DetailLayout.prototype.afterRender.apply(this,arguments);
    },    
    
    showEdit: function() {
      var self = this;

      function getUser(callback){
        if (!self.user && !_.isEmpty(self.model.get('serviced_username'))){
          appModel.getModel(
            'screensaveruser', self.model.get('serviced_username'), 
            function(servicedUserModel){
              self.user = servicedUserModel;
              callback();
            });
        } else {
          callback();
        }
      };
      
        
      function setupEditView(callback){
        var resource = self.model.resource;
        resource.fields['performed_by_username']['choices'] = 
          appModel.getAdminUserOptions();

        if (self.user ){
          var screen_facility_field = resource.fields['screen_facility_id'];
          screen_facility_field['edit_type'] = 'select';
          var choices = appModel._get_user_screen_choices(self.user);
          choices.unshift({ 'val': '', 'label': 'not selected' });
          screen_facility_field['choices'] = choices;
        }
        if (self.screen){
          
          var serviced_username_field = resource.fields['serviced_username'];
          serviced_username_field['edit_type'] = 'select';
          serviced_username_field['choices'] = 
            appModel._get_screen_member_choices(self.screen);
        }
        view = DetailLayout.prototype.showEdit.apply(self,arguments);
      };
      
      $(this).queue([
       getUser,
       appModel.getAdminUserOptions,
       setupEditView]);
      
    }
  });

  return ServiceActivityView;
});