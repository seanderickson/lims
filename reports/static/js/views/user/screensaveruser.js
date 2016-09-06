define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_layout',
    'views/list2',
    'views/user/user2',
    'views/screen/screen',
    'views/generic_edit',
    'views/generic_detail_stickit',
    'templates/generic-tabbed.html',
    'bootstrap-datepicker'
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, ReportsUserView, ScreenView, EditView, DetailView, layout) {

  var UserView = ReportsUserView.extend({

    screensaver_tabbed_resources: {
      screens: {
        description: "Screens associated with this user",
        title: "Screens",
        invoke: "setScreens",
        resource: 'screen'
      },
      userchecklistitem: {
        description: "User Checklist Items",
        title: "User Checklist Items",
        invoke: "setUserChecklistItems",
        resource: 'userchecklistitem'
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
      this.tabbed_resources = _.extend({},
        this.tabbed_resources, this.screensaver_tabbed_resources);
      
      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail' && !appModel.hasPermission(
            self.tabbed_resources[key].resource,'read')){
          delete self.tabbed_resources[key];
        }
      });
    },

    setDetail: function(delegateStack){
      var key = 'detail';
      var self,outerSelf;
      outerSelf = self = this;
      console.log('setDetail: ', delegateStack);
        
      // Create model validation rules based on classification
      // FIXME: attach model validation on showAdd as well
      this.model.validate = function(attrs) {
        var errs = {};
        if ( attrs.classification == 'principal_investigator' &&
            _.isEmpty(attrs.lab_head_affiliation) ){
          errs.lab_head_affiliation = 'Required if PI';
        }
        if (!attrs.username && !attrs.ecommons_id){
          errs.username = errs.ecommons_id = 'Either Ecommons or Username must be entered'
        }
        if (!_.isEmpty(errs)) return errs;
      };
      
      var editView = EditView.extend({
        
        afterRender: function(){
          var self = this;
          console.log('override afterRender');
          
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
          self.$el.find('div[key="lab_head_affiliation"]').append(addLabAffiliationButton);
          
          // Set up lab_head_affiliation availability based on classification
          if(self.getValue('classification') != 'principal_investigator'){
              self.$el.find('[key="form-group-lab_head_affiliation"]').hide();
          }
          // attach classification change listener
          self.listenTo(this, "classification:change", function(e){
            var classification = self.getValue('classification');
            console.log('classification:' + classification)
            if(classification == 'principal_investigator'){
              self.$el.find('[key="form-group-lab_head_affiliation"]').show();
              self.setValue('lab_head_username',self.model.get('username'));
              self.$el.find('[key="lab_head_username"]').find('.chosen-select').trigger("chosen:updated");
            } else {
              self.setValue('lab_head_affiliation','');
              //editForm.setValue('lab_head_username','');
              self.$el.find('[key="lab_head_username"]').find('.chosen-select').trigger("chosen:updated");
              self.$el.find('[key="form-group-lab_head_affiliation"]').hide();
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
            this.$el.find('#generic-form-top').append([
               addSMScreenButton,addRnaiScreenButton,updateRNAUAButton,updateSMUAButton]);
          }
        }
      });
      
      view = new DetailLayout({ 
        model: this.model, 
        uriStack: delegateStack,
        EditView: editView,
        DetailView: detailView
      });
      view.showEdit = function(){
        appModel.initializeAdminMode(function(){
          
          self.model.resource.fields['permissions']['choices'] = 
            appModel.getPermissionsOptions();
          self.model.resource.fields['usergroups']['choices'] = 
            appModel.getUserGroupOptions();
          self.model.resource.fields['lab_head_username']['choices'] = 
            appModel.getPrincipalInvestigatorOptions();

          // Set up the dynamic add lab_affiliation form
          // with the edit view available, set up the lab_head_affiliation rules
          // - add listener to update options dynamically
          // - attach the "add lab affiliation" button to the view
          
          // - add listener to update view options when classification changes
          //    - note we want to replace this with model-driven events (backbone.stickit)
          self.model.on('sync', function(){
            // TODO: should only need to do this if the classification has changed
            // to "PI"; but the changedAttributes are unreliable for detecting this
            if(self.model.get('classification')=='principal_investigator'){
              appModel.unset('principal_investigators');
            }
          });
          
          DetailLayout.prototype.showEdit.call(view,arguments);
        });
      };
//      view.showEdit = function(){
//        this.model.resource.fields['permissions']['choices'] = (
//            appModel.get('permissionOptions'));
//        appModel.getPrincipalInvestigatorOptions(function(options){
//
//          self.model.resource.fields['lab_head_username']['choices'] = options;
//          
//          // Set up the dynamic add lab_affiliation form
//          // with the edit view available, set up the lab_head_affiliation rules
//          // - add listener to update options dynamically
//          // - attach the "add lab affiliation" button to the view
//          
//          // - add listener to update view options when classification changes
//          //    - note we want to replace this with model-driven events (backbone.stickit)
//          self.model.on('sync', function(){
//            // TODO: should only need to do this if the classification has changed
//            // to "PI"; but the changedAttributes are unreliable for detecting this
//            if(self.model.get('classification')=='principal_investigator'){
//              appModel.unset('principal_investigators');
//            }
//          });
//          DetailLayout.prototype.showEdit.call(view,arguments);
//        });  
//      };
      
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
            choiceHash[choice] = vocabulary[choice].title;
          });
        currentAffiliationNames = Iccbl.appModel.getVocabulary('labaffiliation.category.*');
      }catch(e){
        console.log('on get vocabulary', e);
        self.appModel.error('Error locating vocabulary: ' + vocabulary_scope_ref);
      }
      var form_template = [
         "<div class='form-horizontal container' id='addLabAffiliationForm' >",
         "<form data-fieldsets class='form-horizontal container' >",
         "</form>",
         "</div>"].join('');      
      var fieldTemplate = _.template([
        '<div class="form-group" >',
        '    <label class="control-label " for="<%= editorId %>"><%= title %></label>',
        '    <div class="" >',
        '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
        '      <div data-error class="text-danger" ></div>',
        '      <div><%= help %></div>',
        '    </div>',
        '  </div>',
      ].join(''));
      
      var formSchema = {};
      formSchema['affiliation_category'] = {
        title: 'Affiliation Category',
        key: 'affiliation_category',
        type: 'Select',
        validators: ['required'],
        options: choiceHash,
        template: fieldTemplate
      };
      formSchema['affiliation_name'] = {
        title: 'Affiliation Name',
        key: 'affiliation_name',
        type: 'Text',
        validators: ['required'],
        template: fieldTemplate
      };
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        type: 'TextArea',
        template: fieldTemplate
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          console.log('form validate', attrs);
          var errs = {};
          var newVal = attrs['affiliation_name'];
          if (newVal){
            newVal = newVal.toLowerCase().replace(/\W+/g, '_');
            if(_.has(currentAffiliationNames,newVal)){
              errs['affiliation_name'] = '"'+ attrs['affiliation_name'] + '" is already used';
            }
          }
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: _.template(form_template)
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
            var resource = appModel.getResource('vocabulary');
            var key = values['affiliation_name'].toLowerCase().replace(/\W+/g, '_');
            var ordinal = currentAffiliationNames.length + 1;
            var max_item = _.max(currentAffiliationNames, function(affil){ return affil.ordinal });
            if (max_item){
              ordinal = max_item.ordinal + 1;
            }
            var data = {
              'scope': 'labaffiliation.category.' + values['affiliation_category'],
              'key': key,
              'title': values['affiliation_name'],
              'description': values['affiliation_name'],
              'ordinal': ordinal,
              'comment': values['comment']
            };
            
            $.ajax({
              url: resource.apiUri,    
              data: JSON.stringify(data),
              contentType: 'application/json',
              method: 'POST',
              success: function(data){
                appModel.getVocabularies(function(vocabularies){
                  appModel.set('vocabularies', vocabularies);
                });
                appModel.showModalMessage({
                  title: 'Lab Affiliation Created',
                  okText: 'ok',
                  body: '"' + values['affiliation_name'] + '"',
                  ok: function(e){
                    e.preventDefault();
                    // manually add the new lha and lab head username to the multiselects
                    // TODO: possibly use backbone.stickit here to update the options 
                    // on model events.
                    
                    editForm.$el.find('[key="lab_head_affiliation"]')
                      .find('.chosen-select').append($('<option>',{
                        value: data['key']
                      }).text(data['title']));
                    editForm.$el.find('[key="lab_head_affiliation"]')
                      .find('.chosen-select').trigger("chosen:updated");

                    appModel.unset('principal_investigators');
                    appModel.unset('vocabularies');
                    appModel.getVocabularies(function(vocabularies){
                      appModel.set('vocabularies', vocabularies);
                      appModel.getPrincipalInvestigatorOptions(function(vocabulary){
                        var choiceHash = []
                        _.each(_.keys(vocabulary),function(choice){
                          if(vocabulary[choice].is_retired){
                            console.log('skipping retired vocab: ',choice,vocabulary[choice].title );
                          }else{
                            choiceHash.push({ val: choice, label: vocabulary[choice].title });
                          }
                        });
                        
                        editForm.fields['lab_head_username'].editor.setOptions(choiceHash);
                        editForm.$el.find('[key="lab_head_username"]')
                          .find('.chosen-select').trigger("chosen:updated");
                      });
                    });
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
    
    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      UserView.__super__.afterRender.apply(this, arguments);      
    },
    
    click_tab : function(event){
      UserView.__super__.click_tab.apply(this, arguments);      
    },

    change_to_tab: function(key){
      UserView.__super__.change_to_tab.apply(this, arguments);      
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
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();      
    },
    
    setServiceActivities: function(delegateStack) {
      var self = this;
      var key = 'serviceactivity';
      var resource = appModel.getResource('serviceactivity');

      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        // Detail view
        var activity_id = delegateStack.shift();
        self.consumedStack.push(activity_id);
        var _key = activity_id
        appModel.getModel(resource.key, _key, function(model){
          view = new DetailLayout({
            model: model, 
            uriStack: _.clone(delegateStack)
          });
          Backbone.Layout.setupView(view);
          //self.tabViews[key] = view;

          // NOTE: have to re-listen after removing a view
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
        });        
        return;
      }else{
        // List view
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
            self.addServiceActivity(delegateStack);
          });
          showHistoryButton.click(function(e){
            e.preventDefault();
            var newUriStack = ['apilog','order','-date_time', 'search'];
            var search = {};
            search['ref_resource_name'] = 'serviceactivity';
            search['changes__icontains'] = '"serviced_username": "' + self.model.get('username') + '"';
            newUriStack.push(appModel.createSearchString(search));
            var route = newUriStack.join('/');
            console.log('history route: ' + route);
            appModel.router.navigate(route, {trigger: true});
            self.remove();
          });
          if(appModel.hasPermission(self.model.resource.key, 'edit')){
            extraControls.unshift(addServiceActivityButton);
          }
          if(appModel.hasPermission(self.model.resource.key, 'edit')){
            extraControls.unshift(showDeleteButton);
          }
          extraControls.unshift(showHistoryButton);
          console.log('extraControls',extraControls);
          url = [self.model.resource.apiUri, 
                     self.model.key,
                     'serviceactivities'].join('/');
          view = new ListView({ 
            uriStack: _.clone(delegateStack),
            schemaResult: resource,
            resource: resource,
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
          self.consumedStack = [key]; 
          self.reportUriStack([]);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
        })();
      }
    },
    
    addServiceActivity: function(delegateStack) {
      var self = this;
      
      var resource = Iccbl.appModel.getResource('serviceactivity');
      var defaults = {};
      appModel.initializeAdminMode(function(){
        resource.fields['performed_by_username']['choices'] = 
          appModel.getAdminUserOptions();

        _.each(resource.fields, function(value, key){
            if (key == 'resource_uri') {
              defaults[key] = resource.resource_uri;
            } else if (key == 'id'){
            } else {
              defaults[key] = '';
            }
        });
        
        defaults['serviced_username'] = self.model.get('username');
        defaults['serviced_user'] = self.model.get('name');
        defaults['performed_by_username'] = appModel.getCurrentUser().username;
        
        delegateStack.unshift('+add');
        var NewModel = Backbone.Model.extend({
          resource: resource,
          urlRoot: resource.apiUri, 
          defaults: defaults
        });
        var view = new DetailLayout({
          model: new NewModel(), 
          uriStack: _.clone(delegateStack)
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view,'remove',function(){
          self.setServiceActivities([]);
          self.removeView(view); // todo - test
        });
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();

      });
      
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
        EditView.uploadAttachedFileDialog(url, view.collection, 'attachedfiletype.user');
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
      self.reportUriStack([]);
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

    setUserChecklistItems: function(delegateStack) {
      var self = this;
      var key = 'userchecklistitem';
      var resource = appModel.getResource('userchecklistitem');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'checklistitems'].join('/');

      var showSaveButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="save_button" href="#">',
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
        search['ref_resource_name'] = 'userchecklistitem';
        search['changes__icontains'] = '"username": "' + self.model.get('username') + '"';
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        console.log('history route: ' + route);
        appModel.router.navigate(route, {trigger: true});
        self.remove();
      });

      var PostCollection = Backbone.Collection.extend({
        url: url,
        toJSON: function(){
          return {
            objects: Collection.__super__.toJSON.apply(this) 
          };
        }
      });
      var changedCollection = new PostCollection();
      var UserChecklistModel = Backbone.Model.extend({
        url: url,
        initialize : function() {
          this.on('change', function(model, options) {
            // Prevent save on update
            if (options.save === false)
                return;
            model.url = url;
            if(_.isEmpty(model.get('status_date'))){
              model.set('status_date', Iccbl.getISODateString(new Date()));
            }
            if(_.isEmpty(model.get('admin_username'))){
              model.set('admin_username', appModel.getCurrentUser().username);
            }
            changedCollection.add(model);
            appModel.setPagePending();
          });
        },
      });

      var Collection = Iccbl.MyCollection.extend({
        url: url
      });
      collection = new Collection({
        url: url,
      });
      collection.model = UserChecklistModel;

      showSaveButton.click(function(e){
        
        e.preventDefault();
        console.log('changed collection', changedCollection,changedCollection.url);
        
        if(changedCollection.isEmpty()){
          appModel.error('No changes to save');
          return;
        }
        
        appModel.showSaveWithComments(function(formValues){
          console.log('form values', formValues);
          var comments = formValues['comments'];
          var headers = {};
          headers[appModel.HEADER_APILOG_COMMENT] = comments;
          
          Backbone.sync("patch",changedCollection,
            {
              headers: headers,
              error: function(){
                appModel.jqXHRfail.apply(this,arguments);
                console.log('error, refetch', arguments);
                changedCollection.reset();
                collection.fetch({ reset: true });
              },
            }
          );
        });
      });
        
      view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        collection: collection,
        extraControls: [showSaveButton, showHistoryButton]
      });
      Backbone.Layout.setupView(view);
      self.consumedStack = [key]; 
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
    }
    
  });

  return UserView;
});