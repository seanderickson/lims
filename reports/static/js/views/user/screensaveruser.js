define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_layout',
    'views/list2',
    'views/generic_edit',
    'views/generic_detail_stickit',
    'views/user/user2',
    'views/screen/screen',
    'views/activityListView',
    'views/serviceActivity',
    'utils/uploadDataForm',
    'templates/generic-tabbed.html',
    'bootstrap-datepicker'
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, EditView, DetailView, ReportsUserView, ScreenView,
            ActivityListView, ServiceActivityView,
            UploadDataForm, layout) {

  var UserView = ReportsUserView.extend({

    screensaver_tabbed_resources: {
      screens: {
        description: "Screens associated with this user",
        title: "Screens",
        invoke: "setScreens",
        resource: 'screen'
      },
      userchecklist: {
        description: "User Checklist",
        title: "User Checklist",
        invoke: "setUserChecklist",
        resource: 'userchecklist'
      },
      attachedfile: {
        description: "Attached Files",
        title: "Attached Files",
        invoke: "setAttachedFiles",
        resource: 'attachedfile'
      },
      serviceactivity: {
        description: "Service Activities",
        title: "Service Activities",
        invoke: "setServiceActivities",
        resource: 'serviceactivity'
      }
    },
    
    initialize: function(args) {
      UserView.__super__.initialize.apply(this, arguments);      
      var self = this;
      // Merge the parent user.js tabs and modify the tab order
      var keyOrder = ['detail'].concat(_.keys(this.screensaver_tabbed_resources));
      if (_.has(this.tabbed_resources, 'usergrouppermissions')){
        keyOrder.push('usergrouppermissions');
      }
      var tempTabs = _.extend({},
        this.tabbed_resources, this.screensaver_tabbed_resources);
      var orderedTabs = {};
      _.each(keyOrder, function(key){
        orderedTabs[key] = tempTabs[key];
      });
      this.tabbed_resources = orderedTabs;
      
      if (appModel.getCurrentUser().username != this.model.get('username')){
        _.each(_.keys(this.tabbed_resources), function(key){
          if(key !== 'detail' && !appModel.hasPermission(
              self.tabbed_resources[key].resource,'read')){
            delete self.tabbed_resources[key];
          }
        });
      }
    },

    setDetail: function(delegateStack){
      
      console.log('setDetail: ', delegateStack);

      var key = 'detail';
      var self,outerSelf;
      outerSelf = self = this;
      var resource = self.model.resource;
      
      if (_.has(resource.fields,'permissions')){
        resource.fields['permissions']['editability'] = [];
        // TODO: usergroupview permission fields are being altered; should not 
        // affect the resource here.
        resource.fields['permissions']['visibility'] = [];
        resource.fields['usergroups']['editability'] = [];
      }
      
      // Manage classification changes
      var originalClassification = self.model.get('classification');
      if ( originalClassification == 'principal_investigator'){
        resource.fields['lab_member_usernames']['visibility'] = ['d'];
        // set up a custom vocabulary that joins username to name; will be 
        // used as the text of the linklist
        resource.fields['lab_member_usernames'].vocabulary = (
            _.object(self.model.get('lab_member_usernames'),
              self.model.get('lab_member_names')));
      } else {
        resource.fields['lab_member_usernames']['editability'] = [];
      }

      // Create model validation rules based on classification
      // FIXME: attach model validation on showAdd as well
      this.model.validate = function(attrs) {
        var errs = {};
        if ( attrs.classification == 'principal_investigator' &&
            _.isEmpty(attrs.lab_affiliation_id) ){
          errs.lab_affiliation_id = 'Required if PI';
        }
        if (!attrs.username && !attrs.ecommons_id){
          errs.username = errs.ecommons_id = 'Either Ecommons or Username must be entered'
        }
        if (!_.isEmpty(errs)) return errs;
      };
      
      var editView = EditView.extend({
        save_success: function(data, textStatus, jqXHR){
          EditView.prototype.save_success.apply(this,arguments);
          console.log('user saved. unset users...');
          appModel.unset('users');
          console.log('unset done');
          appModel.getUsers();
        },
        
        afterRender: function(){
          var self = this;
          console.log('override afterRender');
          
          // Set up the dynamic add lab_affiliation form
          // - add listener to update options dynamically
          // - attach the "add lab affiliation" button to the view
          
          var addLabAffiliationButton = $([
            '<a class="btn btn-default btn-sm" ',
              'role="button" id="add_button" href="#">',
              'Add</a>'
            ].join(''));
          
          addLabAffiliationButton.click(function(event){
            event.preventDefault();
            outerSelf.addLabAffiliation(self);
          });

          // Render the editForm; then add the add lab affiliation button
          self.$el.find('div[key="lab_affiliation_id"]')
            .append(addLabAffiliationButton);
          
          // Set up lab_head_affiliation availability based on classification
          if(self.getValue('classification') != 'principal_investigator'){
              self.$el.find('[key="form-group-lab_affiliation_id"]').hide();
          }else{
              self.$el.find('[key="form-group-lab_head_username"]').hide();
          }
          // attach classification change listener
          self.listenTo(this, "classification:change", function(e){
            var classification = self.getValue('classification');
            console.log('classification:' + classification)
            if(classification == 'principal_investigator'){
              self.$el.find('[key="form-group-lab_affiliation_id"]').show();
              self.$el.find('[key="form-group-lab_member_usernames"]').show();
              self.$el.find('[key="form-group-lab_head_username"]').hide();
              self.setValue('lab_head_username', '');
            } else {
              self.setValue('lab_affiliation_id','');
              self.setValue('lab_member_usernames','');
              if (originalClassification == 'principal_investigator'){
                console.log('unset lab head username')
                self.setValue('lab_head_username','');
                self.$el.find('[key="lab_head_username"]')
                  .find('.chosen-select').trigger("chosen:updated");
              }
              self.$el.find('[key="form-group-lab_affiliation_id"]').hide();
              self.$el.find('[key="form-group-lab_member_usernames"]').hide();
              self.$el.find('[key="form-group-lab_head_username"]').show();
            }
          });
        
          EditView.prototype.afterRender.call(this,arguments);
        }
      });
      
      var detailView = DetailView.extend({
        afterRender: function(){
          DetailView.prototype.afterRender.call(this,arguments);
          if(appModel.hasPermission('screen', 'write')){
            var addSMScreenButton = $([
              '<a class="btn btn-default btn-sm pull-down" ',
                'role="button" id="add_smscreen" href="#">',
                'Add Small Molecule Screen</a>'
              ].join(''));
            addSMScreenButton.click(function(e){
              e.preventDefault();
              self.addScreen({ screen_type: 'small_molecule' });
            });
            var addRnaiScreenButton = $([
              '<a class="btn btn-default btn-sm pull-down" ',
                'role="button" id="add_rnaiscreen" href="#">',
                'Add RNAi Screen</a>'
              ].join(''));
            addRnaiScreenButton.click(function(e){
              e.preventDefault();
              self.addScreen({ screen_type: 'rnai' });
            });
            var updateSMUAButton = $([
              '<a class="btn btn-default btn-sm pull-down" ',
                'role="button" id="update_smua" href="#">',
                'Update SM User Agreement</a>'
              ].join(''));
            updateSMUAButton.click(function(e){
              e.preventDefault();
              self.updateUserAgreement('sm');
            });
            var updateRNAUAButton = $([
              '<a class="btn btn-default btn-sm pull-down" ',
                'role="button" id="update_rnaua" href="#">',
                'Update RNAi User Agreement</a>'
              ].join(''));
            updateRNAUAButton.click(function(e){
              e.preventDefault();
              self.updateUserAgreement('rnai');
            });
            this.$el.find('#generic-detail-buttonpanel').append([
               addSMScreenButton,addRnaiScreenButton,updateRNAUAButton,updateSMUAButton]);
          }
        }
      });
 
      var ScreenDetailLayout = DetailLayout.extend({
        showEdit: function() {
          var self = this;
          // Lazy load choices
          appModel.initializeAdminMode(function(){
            appModel.getLabAffiliationOptions(function(){
              internalShowEdit();
            });
          });
            
          function internalShowEdit(){
            var fields = self.model.resource.fields;

            fields['lab_head_username']['choices'] = 
              appModel.getPrincipalInvestigatorOptions();
            fields['lab_member_usernames']['choices'] = 
              appModel.getUserOptions();
            fields['lab_member_usernames']['title'] = 'Lab Members';
            fields['lab_affiliation_id']['choices'] = 
              appModel.getLabAffiliationOptions();
  
            // - add listener to update view options when classification changes
            //    - note we want to replace this with model-driven events (backbone.stickit)
            self.model.on('sync', function(){
              // TODO: should only need to do this if the classification has changed
              // to "PI"; but the changedAttributes are unreliable for detecting this
              if(self.model.get('classification')=='principal_investigator'){
                appModel.unset('principal_investigators');
              }
            });
            
            ScreenDetailLayout.__super__.showEdit.apply(self,arguments);
          };
        }
      });
      
      view = new ScreenDetailLayout({ 
        model: this.model, 
        uriStack: delegateStack,
        EditView: editView,
        DetailView: detailView
      });
      
      this.tabViews[key] = view;
      
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // Note: since detail_layout reports the tab, the consumedStack is empty here
      this.consumedStack = []; 
      view = this.setView("#tab_container", view ).render();
      return view;
    },

    addLabAffiliation: function(editForm){
      var self = this;
      var choiceHash = {}
      var currentAffiliationNames, vocabulary;
      var vocabulary_scope_ref = 'labaffiliation.category';
      try{
        vocabulary = Iccbl.appModel.getVocabulary(vocabulary_scope_ref);
          _.each(_.keys(vocabulary),function(choice){
            if (vocabulary[choice].is_retired) {
              console.log('skipping retired vocab: ',choice,vocabulary[choice].title );
            } else {
              choiceHash[choice] = vocabulary[choice].title;
            }
          });
        currentAffiliationNames = Iccbl.appModel.getVocabulary('labaffiliation.category.*');
      }catch(e){
        console.log('on get vocabulary', e);
        self.appModel.error('Error locating vocabulary: ' + vocabulary_scope_ref);
      }
      
      var formSchema = {};
      formSchema['category'] = {
        title: 'Affiliation Category',
        key: 'category',
        type: 'Select',
        validators: ['required'],
        options: choiceHash,
        editorClass: 'form-control',
        template: appModel._field_template
      };
      formSchema['name'] = {
        title: 'Affiliation Name',
        key: 'name',
        type: 'Text',
        validators: ['required'],
        editorClass: 'form-control',
        template: appModel._field_template
      };
      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          console.log('form validate', attrs);
          var errs = {};
            if(_.has(currentAffiliationNames,attrs['name'])){
              errs['name'] = '"'+ attrs['name'] + '" is already used';
            }
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: appModel._form_template
      });
      var _form_el = form.render().el;

      var dialog = appModel.showModal({
        okText: 'create',
        view: _form_el,
        title: 'Create a new Lab Affiliation',
        ok: function(e){
          e.preventDefault();
          var errors = form.commit({ validate: true }); // runs schema and model validation
          if(!_.isEmpty(errors) ){
            console.log('form errors, abort submit: ',errors);
            _.each(_.keys(errors), function(key){
              $('[name="'+key +'"').parents('.form-group').addClass('has-error');
            });
            return false;
          }else{
            var values = form.getValue();
            var resource = appModel.getResource('labaffiliation');
            var data = {
              'category': values['category'],
              'name': values['name'],
            };
            
            $.ajax({
              url: resource.apiUri,    
              data: JSON.stringify(data),
              contentType: 'application/json',
              method: 'POST',
              success: function(data){
                data = data[appModel.API_RESULT_DATA];
                appModel.unset('labAffiliations')
                appModel.unset('labAffiliationOptions')
                appModel.showModalMessage({
                  title: 'Lab Affiliation Created',
                  okText: 'ok',
                  body: '"' + values['name'] + '"',
                  ok: function(e){
                    e.preventDefault();
                    // manually add the affiliations to the multiselects
                    // TODO: possibly use backbone.stickit here to update the options 
                    // on model events.
                    var category = appModel.getVocabularyTitle(
                      'labaffiliation.category', data['category']);
                    editForm.$el.find('[key="lab_affiliation_id"]')
                      .find('.chosen-select').append($('<option>',{
                        value: data['lab_affiliation_id']
                      }).text(category + ' - ' + data['name']));
                    editForm.$el.find('[key="lab_affiliation_id"]')
                      .find('.chosen-select').trigger("chosen:updated");
                  }
                });
              },
              done: function(model, resp){
                // TODO: done replaces success as of jq 1.8
                console.log('done');
              }
            }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });
          
            return true;
          }
        }
      });
    },
    
    addScreen: function(extra_defaults){
      var self = this;
      console.log('add screen for ls: ', self.model.get('username'));
      var defaults = _.extend({
        lead_screener_username: self.model.get('username'),
        lab_head_username: self.model.get('lab_head_username')
      }, extra_defaults);
      var newModel = appModel.createNewModel('screen', defaults);

      var view = new ScreenView({ 
        model: newModel, 
        uriStack: ['+add']
      });
      
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      this.consumedStack = ['screen','+add'];
      self.reportUriStack([]);
    },
        
    setScreens: function(delegateStack) {
      var self = this;
      var addSMScreenButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="add_smscreen" href="#">',
          'Add Small Molecule Screen</a>'
        ].join(''));
      addSMScreenButton.click(function(e){
        e.preventDefault();
        self.addScreen({ screen_type: 'small_molecule' });
      });
      var addRnaiScreenButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="add_rnaiscreen" href="#">',
          'Add RNAi Screen</a>'
        ].join(''));
      addRnaiScreenButton.click(function(e){
        e.preventDefault();
        self.addScreen({ screen_type: 'rnai' });
      });
      var key = 'screens';
      var originalResource = appModel.getResource('screen');
      var resource = _.extend({},originalResource);
      resource.fields = _.pick(
        originalResource.fields,
        ['facility_id','title','screen_type','screensaver_user_role', 'status','status_date',
         'date_of_last_activity','date_created']);
      resource.fields['screen_type']['visibility'] = ['l'];
      resource.fields['screensaver_user_role']['visibility'] = ['l'];
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'screens'].join('/');
      var view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: resource,
          resource: resource,
          url: url,
          extraControls: [addRnaiScreenButton,addSMScreenButton]
      });
      Backbone.Layout.setupView(view);
      self.consumedStack = [key]; 
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();      
    },

    getTitle: function() {
      var self = this;
      var title = [
         '<H4 id="title">',
         '<a href="#' + self.model.resource.key,
         '/{username}" >{first_name} {last_name} ({username})</a>'];
      title.push('</H4>');
      var titleDiv = $(Iccbl.formatString(title.join(''), this.model));
      if (appModel.hasGroup('readEverythingAdmin')){
        if (!_.isEmpty(self.model.get('comments'))){
          titleDiv.append(
            Iccbl.createCommentIcon([self.model.get('comments')],'Commments'));
        }
      }
      return titleDiv;
    },
    
    setServiceActivities: function(delegateStack) {
      var self = this;
      var saResource = appModel.getResource('serviceactivity');
      saResource.fields['serviced_username']['editability'] = [];
      saResource.fields['activity_class']['visibility'] = [];
      
      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        
        if (delegateStack[0]!='+add') {
          var activity_id = delegateStack.shift();
          self.consumedStack.push(activity_id);
          var _key = activity_id
          appModel.getModelFromResource(saResource, _key, function(model){
            view = new ServiceActivityView({
              model: model, 
              user: self.model,
              uriStack: _.clone(delegateStack)
            });
            Backbone.Layout.setupView(view);
            self.listenTo(view , 'uriStack:change', self.reportUriStack);
            self.setView("#tab_container", view ).render();
          });        
          return;
          
        } else {
          
          // Do not allow return to +add screen
          delegateStack.shift();
          self.setServiceActivities(delegateStack);
          return;
        }

      }else{
        // List view
        // FIXME: refactor with ActivityListView
        (function listView(){
          var view, url;
          var extraControls = [];
          var addServiceActivityButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
              'role="button" id="add_button" href="#">',
              'Add</a>'
            ].join(''));
          var showDeleteButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
              'role="button" id="showDeleteButton" href="#">',
              'Delete</a>'
            ].join(''));
          var showHistoryButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
              'role="button" id="showHistoryButton" href="#">',
              'History</a>'
            ].join(''));
          
          addServiceActivityButton.click(function(e){
            e.preventDefault();
            self.addServiceActivity();
          });
          showHistoryButton.click(function(e){
            e.preventDefault();
            var newUriStack = ['apilog','order','-date_time', 'search'];
            var search = {};
            search['ref_resource_name'] = 'serviceactivity';
            search['uri__contains'] = 'screensaveruser/' + self.model.get('username');
            newUriStack.push(appModel.createSearchString(search));
            var route = newUriStack.join('/');
            console.log('history route: ' + route);
            appModel.router.navigate(route, {trigger: true});
            self.remove();
          });
          if(appModel.hasPermission(saResource.key, 'edit')){
            extraControls.unshift(addServiceActivityButton);
          }
          if(appModel.hasPermission(saResource.key, 'edit')){
            extraControls.unshift(showDeleteButton);
          }
          extraControls.unshift(showHistoryButton);
          console.log('extraControls',extraControls);
          url = [self.model.resource.apiUri, 
                     self.model.key,
                     'serviceactivities'].join('/');
          view = new ActivityListView({ 
            uriStack: _.clone(delegateStack),
            schemaResult: saResource,
            resource: saResource,
            url: url,
            extraControls: extraControls
          });
          showDeleteButton.click(function(e){
            e.preventDefault();
            if (! view.grid.columns.findWhere({name: 'deletor'})){
              view.grid.columns.unshift({ 
                name: 'deletor', label: 'Delete', text:'X', 
                description: 'delete record', 
                cell: Iccbl.DeleteCell, sortable: false });
            }
          });
          Backbone.Layout.setupView(view);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
        })();
      }
    },
    
    addServiceActivity: function() {
      var self = this;
      
      var saResource = Iccbl.appModel.getResource('serviceactivity');
      saResource.fields['serviced_username']['editability'] = [];
      
      var defaults = {
        serviced_username: self.model.get('username'),
        serviced_user: self.model.get('name'),
        performed_by_username: appModel.getCurrentUser().username
      };
      var newModel = appModel.newModelFromResource(saResource, defaults);
      var view = new ServiceActivityView({
        model: newModel,
        user: self.model,
        uriStack: ['+add']
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();

      self.consumedStack = ['serviceactivity'];
      self.reportUriStack([]);
      view.reportUriStack(['+add']);

    },
    
    setAttachedFiles: function(delegateStack) {
      var self = this;
      var key = 'attachedfile';
      var resource = appModel.getResource('attachedfile');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'attachedfiles'].join('/');
      var uploadAttachedFileButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="save_button" href="#">',
          'Add</a>'
        ].join(''));
      var showDeleteButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="showDeleteButton" href="#">',
            'Delete</a>'
          ].join(''));
      
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: [uploadAttachedFileButton, showDeleteButton]
      });
      uploadAttachedFileButton.click(function(e){
        e.preventDefault();
        UploadDataForm.uploadAttachedFileDialog(
          url, 'attachedfiletype.user'
        ).done(function(){
          view.collection.fetch({ reset: true });
        }).fail(function(){
          appModel.jqXHRfail.apply(this,arguments); 
        });
      });
      showDeleteButton.click(function(e){
        e.preventDefault();
        if (! view.grid.columns.findWhere({name: 'deletor'})){
          view.grid.columns.unshift({ 
            name: 'deletor', label: 'Delete', text:'X', 
            description: 'delete record', 
            cell: Iccbl.DeleteCell, sortable: false });
        }
      });
      Backbone.Layout.setupView(view);
      self.consumedStack = [key]; 
//      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      
    },
    
    updateUserAgreement: function(type){
      var self = this;
      // TODO: generate these from the appModel, should be vocabulary like
      var choiceHashes = {
        'sm': {
          'smDsl1MutualScreens': 'Small Molecule Screens Level 1',
          'smDsl2MutualPositives': 'Small Molecule Screens Level 2',
          'smDsl3SharedScreens': 'Small Molecule Screens Level 3'
        },
        'rnai': {
          'rnaiDsl1MutualScreens': 'RNAi Screens Level 1',
          //'rnaiDsl2MutualPositives': 'RNAi Screens Level 2',
          'rnaiDsl3SharedScreens': 'RNAi Screens Level 3'
        }
      };
      var choiceHash = choiceHashes[type];
      
      var form_template = [
         "<div class='form-horizontal container' id='uploadUAButton_form' >",
         "<form class='form-horizontal container' >",
         "<div data-fieldsets ></div>",
//         "<div class='form-group' ><input type='file' name='fileInput' /></div>",
         "</form>",
         "</div>"].join('');      
      var fieldTemplate = _.template([
        '<div class="form-group" >',
        '    <label class="control-label col-sm-6" for="<%= editorId %>"><%= title %></label>',
        '    <div class="col-sm-6" >',
        '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
        '      <div data-error class="text-danger" ></div>',
        '      <div><%= help %></div>',
        '    </div>',
        '  </div>',
      ].join(''));
      
      var DisabledField = EditView.DisabledField.extend({
        tagName: 'p',
        className: 'form-control-static'
      });
      var formSchema = {};
      formSchema['user_fullname'] = {
        title: 'User',
        key: 'user_fullname',
        type: DisabledField, 
        template: fieldTemplate
      };
      formSchema['current_value'] = {
        title: 'Current Data Sharing Level',
        key: 'current_value',
        type: DisabledField, 
        template: fieldTemplate
      };
      formSchema['lab_head_dsl'] = {
        title: 'Lab Head Data Sharing Level',
        key: 'lab_head_dsl',
        type: DisabledField, 
        template: fieldTemplate
      };
      formSchema['user_has_login'] = {
        title: 'User has login privileges',
        key: 'user_has_login',
        type: DisabledField, 
        template: fieldTemplate
      };
     formSchema['usergroup'] = {
        title: 'New Data Sharing Level',
        key: 'usergroup',
        type: Backbone.Form.editors.Select.extend({
            className: 'form-control'
          }),
        options: choiceHash,
        template: fieldTemplate,
        validators: ['required']
      };
      formSchema['fileInputPlaceholder'] = {
        title: 'File',
        key: 'file_input',
        type: DisabledField, 
        template: fieldTemplate
      };
      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          console.log('form validate', attrs);
          var errs = {};
          var file = $('input[name="fileInput"]')[0].files[0]; 
          if (file) {
          } else {
            errs.fileInputPlaceholder = 'Must specify a User Agreement File to upload';
          }
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      
      if (type == 'sm'){
        formFields.set('current_value', appModel.getSMUALevel(self.model));
        formFields.set('usergroup', appModel.getSMUALevel(self.model));
      } else {
        formFields.set('current_value', appModel.getRNAUALevel(self.model));
        formFields.set('usergroup', appModel.getRNAUALevel(self.model));
      }
      formFields.set(
        'user_fullname',self.model.get('last_name') + ', ' 
        + self.model.get('first_name'));
      
      if(_.contains(self.model.get('usergroups'), 'screensaverUser')){
        formFields.set('user_has_login', 'yes');
      } else {
        formFields.set('user_has_login', 'no');
      }
      
      function buildForm(){
        var form = new Backbone.Form({
          model: formFields,
          template: _.template(form_template)
        });
        var _form_el = form.render().el;
      
        form.$el.find("[name='fileInputPlaceholder']").append(
          '<label class="btn btn-default btn-file">' + 
          'Browse<input type="file" name="fileInput" style="display: none;"></label>'+ 
          '<p id="filechosen" class="form-control-static" ></p>');
        form.$el.on('change', ':file', function() {
          var input = $(this),
              numFiles = input.get(0).files ? input.get(0).files.length : 1,
              label = input.val().replace(/\\/g, '/').replace(/.*\//, '');
          form.$el.find('#filechosen').html(label);
        });
      
        var dialog = appModel.showModal({
          okText: 'upload',
          ok: function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }); // runs schema and model validation
            if(!_.isEmpty(errors) ){
              console.log('form errors, abort submit: ',errors);
              _.each(_.keys(errors), function(key){
                $('[name="'+key +'"').parents('.form-group').addClass('has-error');
              });
              if (_.has(errors,'_others')){
                $form = $('#uploadUAButton_form');
                $errorDiv = $('<div class="panel text-danger" />');
                _.each(errors['_others'], function(otherError){
                  _.each(_.values(otherError), function(errMsg){
                    $errorDiv.append('Error: ' + errMsg );
                  });
                });
                $form.append($errorDiv);
              }
              return false;
            }else{
              var values = form.getValue();
              var comments = values['comments'];
              var headers = {};
              headers[appModel.HEADER_APILOG_COMMENT] = comments;
              
              var data = new FormData();
              _.each(_.keys(values), function(key){
                data.append(key,values[key])
              });
              data.append('type', type);
              // NOTE: admin user defaults to the current logged in user:
              // - attached_file.created_by_username
              // - userchecklistitem.admin_user
              // - apilog.username
              
              var file = $('input[name="fileInput"]')[0].files[0];
              data.append('attached_file',file);

              var url = [self.model.resource.apiUri, 
                         self.model.key,
                         'useragreement'].join('/');
              $.ajax({
                url: url,    
                data: data,
                cache: false,
                contentType: false,
                processData: false,
                type: 'POST',
                headers: headers, 
                success: function(data){
                  self.model.fetch({ reset: true });
                  appModel.showModalMessage({
                    title: 'Updated User Agreement',
                    okText: 'ok',
                    body: 'data sharing level: ' + values['usergroup'] +
                      '<br>file: ' + file.name + ' '
                  });
                },
                done: function(model, resp){
                  // TODO: done replaces success as of jq 1.8
                  console.log('done');
                }
              }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });
              return true;
            }
          },
          view: _form_el,
          title: 'Update User Agreement'  });
      };
      
      var lab_head_username = self.model.get('lab_head_username');
      if (lab_head_username == self.model.get('username')){
        formFields.set('lab_head_dsl', '&lt;user is lab head&gt;');
        buildForm();
      }else{
        appModel.getModel('screensaveruser',lab_head_username,
          function(model){
            console.log('got model: ', model);
            if (type == 'sm'){
              formFields.set('lab_head_dsl', appModel.getSMUALevel(model));
            } else {
              formFields.set('lab_head_dsl', appModel.getRNAUALevel(model));
            }
            buildForm();
          }
        );
      }
    },

    setUserChecklist: function(delegateStack) {
      var self = this;
      var key = 'userchecklist';
      var resource = appModel.getResource('userchecklist');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'checklist'].join('/');

      var showSaveButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="save_button" style="display: none; href="#">',
          'save</a>'
        ].join(''));
      var showHistoryButton = $([
      '<a class="btn btn-default btn-sm pull-down" ',
        'role="button" id="showHistoryButton" href="#">',
        'History</a>'
      ].join(''));
      showHistoryButton.click(function(e){
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', 'search'];
        var search = {};
        search['ref_resource_name'] = 'userchecklist';
        search['key__contains'] = self.model.get('username') + '/';
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        console.log('history route: ' + route);
        appModel.router.navigate(route, {trigger: true});
        self.remove();
      });
      var Collection = Iccbl.MyCollection.extend({
        url: url
      });
      var fetchCollection = new Collection();

      
      if (appModel.hasPermission('userchecklist', 'write')){
        
        resource.fields['date_effective']['backgridCellType'] = Iccbl.DateCell2.extend({
          
          isEditable: function(){
            if (this.model.edited == true){
              return true;
            } else {
              if (this.model.get('status') != 'not_completed'){
                return true;
              }
            }
            return false;
          },
          
          initialize: function(){
            Iccbl.DateCell2.__super__.initialize.apply(this, arguments);
            var self = this;
            self.model.on('change:date_effective' , function(){
              if (self.model.edited){
                self.$el.addClass('edited');
              } else {
                self.$el.parent().children().removeClass('edited');
              }
            });
          }, 
        });
        resource.fields['status']['editability'] = ['l'];
        resource.fields['status']['backgridCellType'] = Iccbl.BooleanCell.extend({
          initialize: function(){
            Iccbl.BooleanCell.__super__.initialize.apply(this, arguments);
            var self = this;
            this.column.editable = function(){ return false; }
            var model = this.model;
            model.on('change:is_checked' , function(){
              if (model.edited){
                self.$el.addClass('edited');
              } else {
                self.$el.parent().children().removeClass('edited');
              }
            });
          },
          
          cellClicked: function(){
            var model = this.model, column = this.column;
            var checked = !model.get('is_checked');
            model.set('is_checked', checked);
            this.render();
          },

          render: function(){
            var self = this;
            Iccbl.BooleanCell.prototype.render.apply(this, arguments);
            
            var text;
            var rawValue = this.model.get("is_checked");
            var statusValue = this.model.get('status');
            
            if (rawValue != true){
              if(_.isEmpty(this.model.get('date_effective'))){
                text = 'Not Completed';
              } else { 
                text = 'Deactivated';
              }
            } else {
              text = 'Activated';
            }
            this.$el.text(text);
            this.delegateEvents();
            
            return this;
          }
        });
    
        // Set up a "post" collection to track changed models
        var PostCollection = Backbone.Collection.extend({
          url: url,
          // explicitly define the id so that collection compare & equals work
          modelId: function(attrs) {
            return Iccbl.getIdFromIdAttribute( attrs, resource);
          }
        });
        var changedCollection = new PostCollection();
        var UserChecklistModel = Backbone.Model.extend({
          url: url,
          edited: false,
          initialize : function() {
            Backbone.Model.prototype.initialize.call(this,arguments);
            var originalModel = new Backbone.Model(this.attributes);
            this.on('change:is_checked', function(model, options) {
              this.edited = true;
              var checked = model.get('is_checked');
              if (originalModel.get('status') == 'activated'){
                if (checked){
                  this.edited = false;
                }
              }else{
                if (!checked){
                  this.edited = false;
                }
              }
              if (this.edited){
                model.set('date_effective', Iccbl.getISODateString(new Date()));
                var adminUser = appModel.getCurrentUser();
                model.set('admin_username', adminUser.username);
                model.set('admin_name', adminUser.first_name + ' ' + adminUser.last_name);
                changedCollection.add(model);
              } else {
                model.set(originalModel.attributes);
                changedCollection.remove(model);
              }
              
              if (!changedCollection.isEmpty()){
                showSaveButton.show();
                appModel.setPagePending();
              }else{
                showSaveButton.hide();
                appModel.clearPagePending();
              }
            });
            this.on('change:date_effective', function(model, options){
              // if the date is changed separately, then only do not toggle other values
              if (originalModel.get('date_effective') != model.get('date_effective')) {
                this.edited = true;
                changedCollection.add(model);
              } else {
                if (this.edited != true){
                  changedCollection.remove(model);
                }
              }
              if (!changedCollection.isEmpty()){
                showSaveButton.show();
                appModel.setPagePending();
              }
            });
          },
        });
        var Collection = Iccbl.MyCollection.extend({
          url: url
        });
        var fetchCollection = new Collection();
        fetchCollection.model = UserChecklistModel;
    
        showSaveButton.click(function(e){
          
          e.preventDefault();
          console.log('changed collection', changedCollection,changedCollection.url);
          
          if(changedCollection.isEmpty()){
            appModel.error('No changes to save');
            return;
          }
          var headers = {};
          Backbone.sync(
            "patch",
            changedCollection,
            {
              headers: headers,
              error: function(){
                appModel.jqXHRfail.apply(this,arguments);
                console.log('error, refetch', arguments);
                changedCollection.reset();
                fetchCollection.fetch({ reset: true });
              },
              success: function(){
                console.log('success');
                fetchCollection.fetch();
              },
              done: function(){
                console.log('done');
              }
            }
          );
          // 20170615 - removed per JAS/KR; comment not necessary
          //appModel.showOkCommentForm({
          //  ok: function(formValues){
          //    console.log('form values', formValues);
          //    var comments = formValues['comments'];
          //    var headers = {};
          //    headers[appModel.HEADER_APILOG_COMMENT] = comments;
          //    
          //    Backbone.sync(
          //      "patch",
          //      changedCollection,
          //      {
          //        headers: headers,
          //        error: function(){
          //          appModel.jqXHRfail.apply(this,arguments);
          //          console.log('error, refetch', arguments);
          //          changedCollection.reset();
          //          collection.fetch({ reset: true });
          //        },
          //        success: function(){
          //          console.log('success');
          //          collection.fetch();
          //        },
          //        done: function(){
          //          console.log('done');
          //        }
          //      }
          //    );
          //  }
          //}); // showOkCommentForm
        }); // save button click handler
        displayUserChecklist(fetchCollection);
      } else {
        displayUserChecklist(fetchCollection);
      }
      
      function displayUserChecklist(collection) {
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: resource,
          resource: resource,
          url: url,
          collection: collection,
          extraControls: [showSaveButton, showHistoryButton],
          tableClass: 'left-align'
        });
        view.grid.columns.on('update', function(){
          view.$el.find('td').removeClass('edited');
        });
        collection.on('sync', function(){
          console.log('synced');
          view.$el.find('td').removeClass('edited');
          appModel.clearPagePending();
        });
        
        Backbone.Layout.setupView(view);
        self.consumedStack = [key]; 
  //      self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
      };
        
    }
    
  });

  return UserView;
});