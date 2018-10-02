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
      
      var detailView = DetailView.extend( {
        
        afterRender: function(){
          
          DetailView.prototype.afterRender.apply(this,arguments);

          if (appModel.hasPermission('activity', 'write')){
            var add_another_button = $([
              '<button class="btn btn-default btn-sm " role="button" ',
              'id="add_another_sa_button" title="Add another Service Activity" >',
              'Add another Service Activity',
              '</button>'
              ].join(''));
            $('#generic-detail-buttonpanel-left').append(add_another_button);
            add_another_button.click(function(e){
              if (self.screen){
                var uriStack = ['screen', self.screen.key,
                                'activities','+add'];
                console.log('route: ', uriStack);
                appModel.setUriStack(uriStack);
              }else{
                var uriStack = ['screensaveruser', self.user.key,
                                'activity','+add'];
                console.log('route: ', uriStack);
                appModel.setUriStack(uriStack);
              }
            });
          }

        }
      }, self._args);
      
      var ActivityEditView = EditView.extend({
        
        save_success: function(data, textStatus, jqXHR){
          // 20180912 - no need to display (kls4)
          //var meta = _.result(data, 'meta', null);
          //if (meta) {
          //  appModel.showJsonMessages(meta);
          //}
          var urlPath = [];
          if (self.user) {
            urlPath = [self.user.resource.key,self.user.key,'activity'];
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
          var editFormSelf = this;
          EditView.prototype.afterRender.apply(this,arguments);
          if(appModel.hasPermission('usergroup', 'edit')){
            var addServiceActivityPerformerButton = $([
              '<button type="submit" ',
              'class="btn btn-default btn-sm" " >Add Service Activity Performer</input>',
              ].join(''));
            this.$el.find('[key="form-group-performed_by_username"]').append(
              addServiceActivityPerformerButton);
            addServiceActivityPerformerButton.click(function(e){
              e.preventDefault();
              var urlPath = ['usergroup','serviceActivityPerformers','edit'];
              appModel.router.navigate(urlPath.join('/'),{trigger:true});
            });
          }
          
          self.listenTo(editFormSelf, "classification:change", function(e){
            var classification = editFormSelf.getValue('classification');
            console.log('classification change: ', classification);
            var type_vocab_ref = self.model.resource.fields['type']['vocabulary_scope_ref'];
            type_vocab_ref = type_vocab_ref.replace('*', classification);
            var type_vocab = Iccbl.appModel.getVocabularySelectOptions(type_vocab_ref);
            var chosenSelectEl = editFormSelf.$el.find('[key="type"]')
              .find('.chosen-select');
            chosenSelectEl.empty();
            
            _.each(type_vocab, function(v){
              chosenSelectEl.append($('<option>',{
                value: v.val }).text(v.label));
            });
            chosenSelectEl.trigger('chosen:updated');
          });
          
          
        }//editView.afterRender
      }, self._args);//editView
      args.EditView = ActivityEditView;
      args.DetailView = detailView;
      
      DetailLayout.prototype.initialize.call(this,args);

    },
    
    afterRender: function() {
      DetailLayout.prototype.afterRender.apply(this,arguments);
    },    
    
    showEdit: function() {
      var self = this;

      function getUser(callback){
        if (!self.user && !_.isEmpty(self.model.get('serviced_user_id'))){
          appModel.getModel(
            'screensaveruser', self.model.get('serviced_user_id'), 
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
        
        // Only use performed_by_userid from the UI
        self.model.resource.fields['performed_by_username']['editability'] = [];
        appModel.getUserIdsInGroupOptions('serviceActivityPerformers', function(options){
          resource.fields['performed_by_user_id'].choiceHash = options;
  
          if (self.user ){
            var screen_facility_field = resource.fields['screen_facility_id'];
            screen_facility_field['edit_type'] = 'select';
            var choices = appModel._get_user_screen_choices(self.user);
            choices.unshift({ 'val': '', 'label': 'not selected' });
            screen_facility_field.choiceHash = choices;
          }
          if (self.screen){
            
            var serviced_user_field = resource.fields['serviced_user_id'];
            serviced_user_field['edit_type'] = 'select';
            serviced_user_field.choiceHash = 
              appModel._get_screen_member_choices(self.screen);
          }
          view = DetailLayout.prototype.showEdit.apply(self,arguments);
        });
      };
      
      $(this).queue([
       getUser,
       appModel.getAdminUserOptions,
       setupEditView]);
      
    }
  });

  return ServiceActivityView;
});