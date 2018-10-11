define([
    'jquery',
    'underscore',
    'backbone',
    'backgrid',
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
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, 
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
      activity: {
        description: "Activities",
        title: "Activities",
        invoke: "setActivities",
        resource: 'activity'
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
      self.tabbed_resources = orderedTabs;
      
      if (appModel.getCurrentUser().username != this.model.get('username')){
        _.each(_.keys(this.tabbed_resources), function(key){
          if(key !== 'detail' && !appModel.hasPermission(
              self.tabbed_resources[key].resource,'read')){
            delete self.tabbed_resources[key];
          }
        });
      }
      
      if (this.model.get('is_staff') != true){
        delete self.tabbed_resources['usergrouppermissions'];
      }
      
      if (!this.model.has('screens')){
        delete self.tabbed_resources['screens'];
      }
      _.bindAll(this,'addServiceActivity');
    },

    setLabHeadFieldVisibility: function(classification, originalClassification){
      var self = this;
      if(classification == appModel.VOCAB_USER_CLASSIFICATION_PI){
        self.$el.find('[key="form-group-lab_affiliation_id"]').show();
        self.$el.find('[key="form-group-lab_member_ids"]').show();
        self.$el.find('[key="form-group-lab_head_id"]').hide();
        self.$el.find('[key="form-group-lab_head_appointment_category"]').show();
        self.$el.find('[key="form-group-lab_head_appointment_department"]').show();
        self.$el.find('[key="form-group-lab_head_appointment_update_date"]').show();
      } else {
        self.$el.find('[key="form-group-lab_affiliation_id"]').hide();
        self.$el.find('[key="form-group-lab_member_ids"]').hide();
        self.$el.find('[key="form-group-lab_head_id"]').show();
        self.$el.find('[key="form-group-lab_head_appointment_category"]').hide();
        self.$el.find('[key="form-group-lab_head_appointment_department"]').hide();
        self.$el.find('[key="form-group-lab_head_appointment_update_date"]').hide();
      }
    },
    
    setDetail: function(delegateStack){
      
      var key = 'detail';
      var self,outerSelf;
      outerSelf = self = this;
      var resource = appModel.getResource(self.model.resource.key);
      if (_.has(resource.fields,'permissions')){
        resource.fields['permissions']['editability'] = [];
        resource.fields['permissions']['visibility'] = [];
        resource.fields['usergroups']['editability'] = [];
      }
      
      if (!_.isEmpty(self.model.get('ecommons_id'))){
        resource.fields['ecommons_id']['visibility'].push('e');
        resource.fields['ecommons_id']['editability']= [];
        resource.fields['username']['visibility']= [];
        resource.fields['username']['editability']= [];
      }
      else if (!_.isEmpty(self.model.get('username'))){
        // NOTE: only checking username; if ecommons id is set, username should
        // also be set
        resource.fields['username']['visibility'].push('e');
        resource.fields['username']['editability']= [];
        resource.fields['ecommons_id']['editability']= [];
      }
      
      var vocab_scope_ref = 
        self.model.resource.fields['rnai_data_sharing_level']['vocabulary_scope_ref'];
      var vocab_dsl = Iccbl.appModel.getVocabularySelectOptions(vocab_scope_ref);
      vocab_dsl = _.reject(vocab_dsl, function(obj){
        return obj.val == 2;
      });
      resource.fields['rnai_data_sharing_level'].choiceHash = vocab_dsl;
      
      var originalClassification = self.model.get('classification');
      if ( originalClassification == appModel.VOCAB_USER_CLASSIFICATION_PI){
        resource.fields['classification']['editability'] = [];
        resource.fields['classification']['visibility'] = ['d','e']
        resource.fields['lab_head_appointment_category']['visibility'] = ['d','e']
        resource.fields['lab_head_appointment_department']['visibility'] = ['d','e']
        resource.fields['lab_head_appointment_update_date']['visibility'] = ['d','e']
        resource.fields['classification']['visibility'] = ['d','e']
        resource.fields['lab_member_ids']['visibility'] = ['d'];
        resource.fields['lab_member_ids'].vocabulary = {};
        if (self.model.has('lab_member_ids')){
          var lab_member_ids = self.model.get('lab_member_ids');
          var lab_member_names = self.model.get('lab_member_names');
          var lab_member_emails = self.model.get('lab_member_emails');
          for (var i=0; i<lab_member_ids.length; i++){
            var name = lab_member_names[i];
            var email = lab_member_emails[i];
            if (!_.isEmpty(email) && email != 'null'){
              name += ' &lt;' + email + '&gt;';
            }
            resource.fields['lab_member_ids'].vocabulary[lab_member_ids[i]] = name;
          }
        }
      } else {
        resource.fields['lab_member_ids']['editability'] = [];
      }

      this.model.validate = function(attrs) {
        var errs = {};
        if ( attrs.classification == appModel.VOCAB_USER_CLASSIFICATION_PI &&
            !_.isNumber(attrs.lab_affiliation_id) ){
          errs.lab_affiliation_id = 'Required if PI';
        }
        if (!attrs.username && !attrs.ecommons_id){
          if (attrs.is_active == true){
            errs.username = errs.ecommons_id = 
              'Either Ecommons or Username must be entered for active users'
          }
        }
        if (!_.isEmpty(errs)) return errs;
      };
      
      var editView = EditView.extend({
        save_success: function(data, textStatus, jqXHR){
          EditView.prototype.save_success.apply(this,arguments);
          if (appModel.DEBUG){
            console.log('user saved. unset users...');
          }
          appModel.unset('users');
          appModel.getUsers();
        },
        
        afterRender: function(){
          var self = this;
          
          var addLabAffiliationButton = $([
            '<a class="btn btn-default btn-sm" ',
              'role="button" id="add_button" href="#">',
              'Add</a>'
            ].join(''));
          
          addLabAffiliationButton.click(function(event){
            event.preventDefault();
            outerSelf.addLabAffiliation(self);
          });

          self.$el.find('div[key="lab_affiliation_id"]')
            .append(addLabAffiliationButton);
          
          if(self.getValue('classification') != appModel.VOCAB_USER_CLASSIFICATION_PI){
              self.$el.find('[key="form-group-lab_affiliation_id"]').hide();
          }else{
              self.$el.find('[key="form-group-lab_head_id"]').hide();
          }
          outerSelf.setLabHeadFieldVisibility(originalClassification, originalClassification);
          
          self.listenTo(this, "classification:change", function(e){
            var classification = self.getValue('classification');
            outerSelf.setLabHeadFieldVisibility(classification, originalClassification);
            if(classification == appModel.VOCAB_USER_CLASSIFICATION_PI){
              self.setValue('lab_head_id', '');
            } else {
              self.setValue('lab_affiliation_id','');
              self.setValue('lab_member_ids','');
              if (originalClassification == appModel.VOCAB_USER_CLASSIFICATION_PI){
                self.setValue('lab_head_id','');
                self.$el.find('[key="lab_head_id"]')
                  .find('.chosen-select').trigger("chosen:updated");
              }
            }
          });
          
          self.$el.find('input[name="is_active"]').click(function(e){
            if (outerSelf.model.get('is_staff') != true){
              if (e.target.checked){
                e.preventDefault();
                appModel.showModal({
                  okText: 'Allow User Login',
                  title: 'Perform administrative override of the login setting?',
                  body: 'Use the "Update User Agreement" buttons to enable login ' +
                  'for screeners.<br/> Selecting this option will manually override the ' + 
                  'User Agreement requirement for login.'
                  ,
                  ok: function(){
                    self.setValue('is_active', true);
                  }
                });
              }
            }
          });
        
          EditView.prototype.afterRender.apply(this,arguments);
        }
      });
      
      var detailView = DetailView.extend({
        afterRender: function(){
          DetailView.prototype.afterRender.apply(this,arguments);
          if(appModel.hasPermission('screensaveruser', 'write')){
            
            var addLabMemberButton = $([
              '<div class="pull-down">',
              '<a class="btn btn-default btn-sm" ',
                'title="Create a new user as a member of this lab" ',
                'role="button" id="add_lab_member_button" href="#">',
                'Add Lab Member</a>',
                '</div>'
              ].join(''));
            this.$el.find('#lab_member_ids').append(addLabMemberButton);
            addLabMemberButton.click(function(e){
              e.preventDefault();
              self.addLabMember();
            });
            
          }
          
          self.showUserAgreements();
          
        }
      });
 
      var ScreensaverUserDetailLayout = DetailLayout.extend({
        showEdit: function() {
          var self = this;
          // Lazy load choices
          appModel.initializeAdminMode(function(){
            appModel.getLabAffiliationOptions(function(){
              internalShowEdit();
            });
          });
          outerSelf.$('ul.nav-tabs > li').addClass('disabled');
            
          function internalShowEdit(){
            var fields = resource.fields;

            fields['lab_head_id'].choiceHash = 
              appModel.getPrincipalInvestigatorOptions();
            fields['lab_member_ids']['title'] = 'Lab Members';
            fields['lab_affiliation_id'].choiceHash = 
              appModel.getLabAffiliationOptions();
  
            // - add listener to update view options when classification changes
            //    - note we want to replace this with model-driven events (backbone.stickit)
            self.model.on('sync', function(){
              // TODO: should only need to do this if the classification has changed
              // to "PI"; but the changedAttributes are unreliable for detecting this
              if(self.model.get('classification')==appModel.VOCAB_USER_CLASSIFICATION_PI){
                appModel.unset('principal_investigators');
              }
            });
            
            ScreensaverUserDetailLayout.__super__.showEdit.apply(self,arguments);
          };
        }
      });
      
      view = new ScreensaverUserDetailLayout({ 
        model: this.model, 
        resource: resource,
        uriStack: delegateStack,
        EditView: editView,
        DetailView: detailView,
        table_class: 'col-sm-12',
        label_col_class: 'col-xs-2',
        value_col_class: 'col-xs-10',
        buttons: ['edit', 'history','download']
      });
      
      this.tabViews[key] = view;
      
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.consumedStack = []; 
      view = this.setView("#tab_container", view ).render();
      return view;
    },
    
    showUserAgreements: function(){
      var self = this;
      if (!appModel.hasPermission('useragreement')){
        console.log('user does not have permission to query "useragreement"');
        return;
      }
      var hasWritePermission = appModel.hasPermission('useragreement','write');

      var uaResource = appModel.getResource('useragreement');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: uaResource.apiUri 
      });
      
      var ua_collection = new CollectionClass();
      ua_collection.fetch({
        data: { 
          limit: 0,
          screensaver_user_id: self.model.get('screensaver_user_id'),
          order_by: ['type'],
          includes: '*'
        },
        success: build_table
      }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      
      function build_table(collection) {
        if (collection.isEmpty()) {
          return;
        }
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell.extend({
            render: function () {
              this.$el.empty();
              var column = this.column;
              var label = Iccbl.createLabel(column.get("label"), 10);
              this.$el.append(label);
              this.$el.prop('title', column.get("description"))
              return this;
            }
          })
        };
        
        var typeVocabRef = uaResource.fields['type']['vocabulary_scope_ref'];
        var dslVocabRef = uaResource.fields['data_sharing_level']['vocabulary_scope_ref'];
        var statusVocabRef = uaResource.fields['status']['vocabulary_scope_ref'];
        var filenameHref = uaResource.fields['filename']['display_options']['hrefTemplate']
        var columns = [];
        var key = 'type';
        columns.push(_.extend({},colTemplate,{
            'name' : key,
            'label' : uaResource.fields[key]['title'],
            'description' : uaResource.fields[key]['description'],
            'order': 1,
            'sortable': true,
            'cell': Iccbl.SelectCell.extend({
              optionValues:
                Iccbl.appModel.getVocabularySelectCellArray(
                  typeVocabRef),
              vocabulary_scope_ref: typeVocabRef 
            })
          }));
        key = 'data_sharing_level';
        columns.push(_.extend({},colTemplate,{
            'name' : key,
            'label' : uaResource.fields[key]['title'],
            'description' : uaResource.fields[key]['description'],
            'order': 1,
            'sortable': true,
            'cell': Iccbl.SelectCell.extend({
              optionValues:
                Iccbl.appModel.getVocabularySelectCellArray(
                  dslVocabRef),
              vocabulary_scope_ref: dslVocabRef 
            })
          }));
        key = 'status';
        columns.push(_.extend({},colTemplate,{
            'name' : key,
            'label' : uaResource.fields[key]['title'],
            'description' : uaResource.fields[key]['description'],
            'order': 1,
            'sortable': true,
            'cell': Iccbl.SelectCell.extend({
              optionValues:
                Iccbl.appModel.getVocabularySelectCellArray(
                  statusVocabRef),
              vocabulary_scope_ref: statusVocabRef 
            })
          }));
        key = 'date_active';
        columns.push(_.extend({},colTemplate,{
            'name' : key,
            'label' : uaResource.fields[key]['title'],
            'description' : uaResource.fields[key]['description'],
            'order': 1,
            'sortable': true,
            'cell': 'Date'
          }));
        key = 'date_notified';
        columns.push(_.extend({},colTemplate,{
            'name' : key,
            'label' : uaResource.fields[key]['title'],
            'description' : uaResource.fields[key]['description'],
            'order': 1,
            'sortable': true,
            'cell': 'Date'
          }));
        key = 'date_expired';
        columns.push(_.extend({},colTemplate,{
            'name' : key,
            'label' : uaResource.fields[key]['title'],
            'description' : uaResource.fields[key]['description'],
            'order': 1,
            'sortable': true,
            'cell': 'Date'
          }));
        key = 'filename';
        columns.push(_.extend({},colTemplate,{
            'name' : key,
            'label' : uaResource.fields[key]['title'],
            'description' : uaResource.fields[key]['description'],
            'order': 1,
            'sortable': true,
            'cell': Iccbl.LinkCell.extend({
              hrefTemplate: filenameHref                
            })
          }));
        if (hasWritePermission) {
          columns.push(
            _.extend({},colTemplate,{
              'name' : 'action',
              'label' : 'Action',
              'description' : 'Update/Reset',
              'order': 1,
              'sortable': false,
              'cell': Iccbl.TextWrapCell.extend({
                className: 'text-wrap-cell-extra-narrow',
                render : function() {
                  var innerself = this;
                  this.$el.empty();
                  this.$el.css('text-align','center');
                  this.$el.append(
                    $('<input type="button" class="btn btn-default btn-sm" value="Update" />')
                      .on('click', function(){
                        self.updateUserAgreement(innerself.model); 
                      })
                    );
                  this.$el.append(
                    $('<input type="button" class="btn btn-default btn-sm" value="History" />')
                      .on('click', function(){
                        var newUriStack = ['apilog','order','-date_time', appModel.URI_PATH_SEARCH];
                        var search = {};
                        search['ref_resource_name'] = uaResource.key;
                  
                        //search['key'] = this.model.key;
                        search['key'] = encodeURIComponent(self.model.key 
                          + '/' + innerself.model.get('type'));
                        
                        newUriStack.push(appModel.createSearchString(search));
                        var route = newUriStack.join('/');
                        appModel.router.navigate(route, {trigger: true});
                      })
                    );

                  return this;
                }
              })
            })
          );
        }
        
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();

        // Walk the DOM to find/replace the DSL fields with the 
        // User Agreement table
        // NOTE: if not found, may already exist
        $('#sm_data_sharing_level').closest('tr').remove();
        var $target_el = $('#rnai_data_sharing_level');
        var $container;
        if ($target_el.length ){
          // First time, create a container for the table
          $target_el = $target_el.closest('tr');
          $target_el.empty();
          $target_el = $target_el.parent().closest('td');
          $container = $('<div id="agreement_table"></div>');
          $target_el.prepend($container);
        } else {
          // Look for the container
          $container = $('#agreement_table');
          if ($container.length){
            $container.empty();
          } else {
            appModel.error('Programming Error: Agreement Table container not found');
          }
        }
        var ua_grid = new Backgrid.Grid({
          columns: colModel,
          collection: collection,
          className: 'backgrid table-striped table-condensed table-hover '
        });
        $container.html(ua_grid.render().$el);
        
      }; // build_table
      
    },

    updateUserAgreement: function(uaModel){
      
      var self = this;
      var uaResource = appModel.getResource('useragreement');
      var type = uaModel.get('type');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'useragreement'].join('/');
      var form_template = [
         "<div class='form-horizontal container' id='uploadUAButton_form' >",
         "<form class='form-horizontal container' >",
         "<div data-fieldsets ></div>",
         "</form>",
         "</div>"].join('');      
      var fieldTemplate = appModel._horizontal_form_field_template;
      var DisabledField = EditView.DisabledField.extend({
        tagName: 'div',
        className: 'form-control-static'
      });
      var formSchema = {};

      var lab_head_id = self.model.get('lab_head_id');
      function showLabHeadWarning(msg){
        appModel.showModal({
          body: msg,
          title: 'Proceed?',
          ok: function(){
            buildForm();
            return false;
          }
        });
      };
      
      if (uaModel.get('status') == 'active'){
        formSchema['deactivate'] = {
          title: 'Deactivate',
          key: 'deactivate',
          type: 'Checkbox', 
          template: fieldTemplate,
          help: 'Select deactivate to null out the previous agreement settings, '
            + 'and remove login capability'
        };
        // 20171101 - expire not shown in UI - JAS
        //formSchema['expire'] = {
        //  title: 'Expire',
        //  key: 'expire',
        //  type: 'Checkbox', 
        //  template: fieldTemplate
        //};
      }
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
        title: 'User currently has login privileges',
        key: 'user_has_login',
        type: DisabledField, 
        template: fieldTemplate
      };
      formSchema['date_active'] = {
        title: 'Date Active',
        key: 'date_active',
        type: EditView.DatePicker, 
        template: fieldTemplate
      };
      
      var vocab_scope_ref = 
        self.model.resource.fields['sm_data_sharing_level']['vocabulary_scope_ref'];
      var vocab_dsl = Iccbl.appModel.getVocabularySelectOptions(vocab_scope_ref);
      if (type=='rnai'){
        vocab_dsl = _.reject(vocab_dsl, function(obj){
          return obj.val == 2;
        });
      }
      formSchema['data_sharing_level'] = {
        title: 'New Data Sharing Level',
        key: 'data_sharing_level',
        type: Backbone.Form.editors.Select.extend({
            className: 'form-control'
          }),
        options: vocab_dsl,
        template: fieldTemplate,
        validators: []
      };
      formSchema['fileInputPlaceholder'] = {
        title: 'File',
        key: 'file_input',
        type: EditView.DisabledField.extend({
          tagName: 'div',
          className: ''
        }), 
        template: fieldTemplate
      };
      formSchema['comments'] = {
        title: 'Comments (optional)',
        key: 'comments',
        type: 'TextArea',
        editorClass: 'input-full',
        template: fieldTemplate
      };
      
      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          var errs = {};
          
          var dateVal = new Date(attrs['date_active']);
          dateVal = new Date(
            dateVal.getFullYear(),dateVal.getMonth(),dateVal.getDate());
          var currentDate = new Date();
          currentDate = new Date(
            currentDate.getFullYear(),currentDate.getMonth(),currentDate.getDate());
          if (dateVal > currentDate){
            errs['date_active'] = 'May not be in the future';
          }
          
          var dsl = attrs['data_sharing_level'];
          if (_.isNumber(uaModel.get('lab_head_data_sharing_level'))){
            if (dsl != uaModel.get('lab_head_data_sharing_level')){
              if (self.model.get('classification') != appModel.VOCAB_USER_CLASSIFICATION_PI){
                errs['data_sharing_level'] = 'Must match the Lab Head';
              }
            }
          }
          
          if (attrs['expire'] != true && attrs['deactivate'] != true){
            // Only require DSL if not expiring/deactivating
            var new_dsl = attrs['data_sharing_level'];
            if ( new_dsl === null || new_dsl == ''){
              errs['data_sharing_level'] = 'Required'
            }
            if (!uaModel.has('filename') || uaModel.has('date_expired')){
              var file = $('input[name="fileInput"]')[0].files[0];
              if (!file){
                errs['fileInputPlaceholder'] = 'Attached file is required to activate';
              }
            }
          }
          
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var typeLabel = appModel.getVocabularyTitle(
        uaResource.fields['type']['vocabulary_scope_ref'], type);
      formTitle = Iccbl.formatString(
        'Update {type} User Agreement', { type: typeLabel});
      formFields.set(
        'user_fullname',self.model.get('last_name') + ', ' 
        + self.model.get('first_name'));
      
      formFields.set('user_has_login', 'false');
      if (self.model.get('is_active')==true){
        formFields.set('user_has_login', 'true');
      }
      
      if (uaModel.get('status')=='active'){
        // NOT RESET CASE - allow patching
        formFields.set('date_active', uaModel.get('date_active'));
        formFields.set('data_sharing_level', uaModel.get('data_sharing_level'));
        formFields.set('current_value', appModel.getVocabularyTitle(
            uaResource.fields['data_sharing_level']['vocabulary_scope_ref'], 
            uaModel.get('data_sharing_level')));
      } else {
        formFields.set('date_active', new Date());
        var default_dsl = uaModel.get('lab_head_data_sharing_level');
        if (_.isEmpty(default_dsl)){
          default_dsl = 3;
        }
        formFields.set('data_sharing_level', default_dsl);
      }
      
      if (self.model.get('classification') == appModel.VOCAB_USER_CLASSIFICATION_PI){
        formFields.set('lab_head_dsl', '&lt;User is Lab Head&gt;');
      } else {
        if (!_.isNumber(uaModel.get('lab_head_data_sharing_level'))){
          showLabHeadWarning(
            'Lab Head has no active ' + typeLabel + ' - Data Sharing Level');
        } else {
          formFields.set('lab_head_dsl', appModel.getVocabularyTitle(
            uaResource.fields['data_sharing_level']['vocabulary_scope_ref'], 
            uaModel.get('lab_head_data_sharing_level')));
        }
      }

      function buildForm(){
        var form = new Backbone.Form({
          model: formFields,
          template: _.template(form_template)
        });
        var _form_el = form.render().el;
      
        var file_input_txt = '<label id="file-button" class="btn btn-default btn-file">' + 
          'Browse<input type="file" name="fileInput" style="display: none;"></label>'; 
        if (uaModel.has('filename')) {
          formFields.set('file_input', '(' + uaModel.get('filename') + ')');
          file_input_txt += 
            '<p id="filechosen" class="form-control-static" >&lt;' + 
            uaModel.get('filename') + '&gt;</p>';
        } else {
          file_input_txt += '<p id="filechosen" class="form-control-static" ></p>';
        }
        form.$el.find("[name='fileInputPlaceholder']").append(file_input_txt);
        form.$el.on('change', ':file', function() {
          var input = $(this),
              numFiles = input.get(0).files ? input.get(0).files.length : 1,
              label = input.val().replace(/\\/g, '/').replace(/.*\//, '');
          form.$el.find('#filechosen').html(label);
        });
        form.listenTo(form, "deactivate:change", function(e){
          var value = form.getValue('deactivate');
          if(value == true){
            form.$el.find('input').not('[name="deactivate"]').prop('disabled', true);
            form.$el.find('select').prop('disabled', true);
            form.$el.find('#file-button').addClass('disabled');
          } else {
            form.$el.find('input').prop('disabled', false);
            form.$el.find('select').prop('disabled', false);
            form.$el.find('#file-button').removeClass('disabled');
          }
        });
        form.listenTo(form, "expire:change", function(e){
          var value = form.getValue('expire');
          if(value == true){
            form.$el.find('input').not('[name="expire"]').prop('disabled', true);
            form.$el.find('select').prop('disabled', true);
            form.$el.find('#file-button').addClass('disabled');
          } else {
            form.$el.find('input').prop('disabled', false);
            form.$el.find('select').prop('disabled', false);
            form.$el.find('#file-button').removeClass('disabled');
          }
        });
        
        var dialog = appModel.showModal({
          okText: 'ok',
          view: _form_el,
          title: formTitle,
          ok: function(e){
            e.preventDefault();
            $('#general_error').remove();
            var errors = form.commit({ validate: true }); // runs schema and model validation
            if(!_.isEmpty(errors) ){
              _.each(_.keys(errors), function(key){
                $('[name="'+key +'"]').parents('.form-group').addClass('has-error');
              });
              if (_.has(errors,'_others')){
                $form = $('#uploadUAButton_form');
                $errorDiv = $('<div id="general_error" class="panel text-danger" />');
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

              var headers = {};
              var comments = values['comments'];
              headers[appModel.HEADER_APILOG_COMMENT] = comments;
              
              if (_.result(values,'expire')){
                appModel.showModalError('Manual Expiration not allowed');
                return 
              }
              else if (_.result(values,'deactivate')){
                var contentType = 'application/json';
                var data = JSON.stringify({ 
                  status: 'inactive', type: type });
                
                $.ajax({
                  url: url,    
                  data: data,
                  cache: false,
                  contentType: 'application/json',
                  dataType: 'json', // what is expected back from the server
                  processData: false,
                  type: 'PATCH',
                  headers: headers, 
                  success: function(data){
                    self.showUserAgreements();
                    self.model.fetch({ reset: true }).done(function(){
                      self.buildMessages();
                    });
                    appModel.showModalMessage({
                      title: formTitle + ' - complete',
                      okText: 'ok',
                      body: 'deactivated'
                    });
                  }
                }).fail(function(jqXHR, textStatus, errorThrown){ 
                  Iccbl.appModel.jqXHRfail.apply(this,arguments); 
                }).done(function(model, resp){
                    // TODO: done replaces success as of jq 1.8
                });
                return true;
              }

              // Only process form if not expiring/deactivating
              var mutableFields = _.filter(_.keys(uaResource.fields), function(key){
                return _.contains(uaResource.fields[key]['editability'],'u');
              });
              var lab_head_dsl = values['lab_head_dsl'];
              var dsl = values['data_sharing_level'];
              if (_.isNumber(lab_head_dsl)){
                if( lab_head_dsl != dsl){
                  appModel.showModal({
                    title: 'Continue?',
                    body: 'Selected Data Sharing Level is different than Lab Head DSL',
                    
                  });
                }
              }
              // TODO: Should admin be required to override if DSL is less 
              // restrictive than screen?              
              // (Show a warning if DSL does not match one of user's screens)
              //var warnScreens = _.filter(self.getUserScreens(), function(screen){
              //  return ( screen.get('screen_type') == type 
              //      && screen.get('data_sharing_level') > dsl);
              //});
              //if (!_.isEmpty(warnScreens)){
              //  $errorDiv = $('<div id="general_error" class="panel text-danger" />');
              //  $errorDiv.append('Error: User screens are more restrictive: (' + 
              //    _.pluck(warnScreens, 'facility_id').join(', ') + ')');
              //  $form = $('#uploadUAButton_form');
              //  $form.append($errorDiv);
              //}
              var httpType = 'PATCH';
              var contentType = 'application/json';
              var file = $('input[name="fileInput"]')[0].files[0];
              var data = {};
              
              if (file) {
                // Send as multipart/form-data
                httpType = 'POST'; 
                contentType = false;
                data = new FormData();
                data.append('type', type);
                data.append('attached_file',file);
                _.each(_.keys(values), function(key){
                  if (_.contains(mutableFields,key)){
                    data.append(key,values[key])
                  }
                });
              } else {
                _.each(_.keys(values), function(key){
                  if (_.contains(mutableFields,key)){
                    data[key] = values[key];
                  }
                });
                data['type'] = type;
                data = JSON.stringify(data);
              }
              
              function processPostPatch(patchUrl){
                $.ajax({
                  url: patchUrl,
                  data: data,
                  cache: false,
                  contentType: contentType,
                  dataType: 'json', // what is expected back from the server
                  processData: false,
                  type: httpType,
                  headers: headers, 
                  success: function(data){
                    self.showUserAgreements();
                    self.model.fetch({ reset: true }).done(function(){
                      self.buildMessages();
                    });
                    
                    var messages = {};
                    messages['data_sharing_level'] =  values['data_sharing_level'];
                    if (file){
                      messages['filename'] = file.name;
                    }
                    var meta = _.result(data, 'meta', null);
                    if (meta) {
                      _.extend(messages, meta);
                    }
                    appModel.showJsonMessages(messages, {
                      title: 'Update Success'
                    });
                  }
                })
                .fail(processError)
                .done(function(model, resp){
                    // TODO: done replaces success as of jq 1.8
                });
              };
              
              function processError(jqXHR, textStatus, errorThrown){
                var jsonError = _.result(jqXHR, 'responseJSON');
                if (!_.isUndefined(jsonError)){
                  var error = _.result(jsonError, 'errors');
                  var overrideFlag = _.result(error,appModel.API_PARAM_OVERRIDE)
                  if (!_.isUndefined(overrideFlag)){
                    var bodyMsg = appModel.print_json(error);
                    body = $('<textarea class="input-full" rows=' 
                      + appModel.MAX_ROWS_IN_DIALOG_MSG + ' ></textarea>');
                    body.val(bodyMsg);

                    appModel.showModal({
                      okText: 'Confirm Override',
                      title: 'Override Required: New Data Sharing Level will invalidate Lab Member Agreements:',
                      view: body,
                      ok: function(){
                        processPostPatch(url + '?' + appModel.API_PARAM_OVERRIDE + '=true');
                      }
                    });
                    
                  }
                } else {
                  Iccbl.appModel.jqXHRfail.apply(this,arguments); 
                }                
              };

              processPostPatch(url);
              return true;
            }
          }
        });
      };
      buildForm();

    },
    
    addLabAffiliation: function(editForm){
      var self = this;
      var choiceHash = {}
      var currentAffiliationNames, vocabulary;
      var vocabulary_scope_ref = 'labaffiliation.category';
      try{
        vocabulary = Iccbl.appModel.getVocabulary(vocabulary_scope_ref);
          _.each(_.keys(vocabulary),function(choice){
            if (vocabulary[choice].is_retired !== true ) {
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
                if (_.isArray(data)){
                  data = data[0];
                }
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
              }
            }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });
          
            return true;
          }
        }
      });
    },
    
    addScreen: function(extra_defaults){
      var self = this;
      var defaults = _.extend({
        lead_screener_id: self.model.get('screensaver_user_id'),
        lab_head_id: self.model.get('lab_head_id')
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
        
    addLabMember: function(extra_defaults){
      var self = this;
      self.$('ul.nav-tabs > li').addClass('disabled');

      if (appModel.DEBUG){
        console.log('add lab member for lab head: ', self.model.get('screensaver_user_id'));
      }
      var defaults = _.extend({
        lab_head_id: self.model.get('lab_head_id'),
        lab_name: self.model.get('lab_name')
      }, extra_defaults);
      var newModel = appModel.createNewModel('screensaveruser', defaults);
      var resource = appModel.getResource('screensaveruser');
      newModel.resource = resource;
      resource.fields['lab_name']['visibility'] = ['e'];
      
      var vocab_scope_ref = 
        resource.fields['classification']['vocabulary_scope_ref'];
      var vocab_cls = Iccbl.appModel.getVocabularySelectOptions(vocab_scope_ref);
      vocab_cls = _.reject(vocab_cls, function(obj){
        return obj.val == appModel.VOCAB_USER_CLASSIFICATION_PI;
      });
      resource.fields['classification'].choiceHash = vocab_cls;
      resource.fields['classification']['vocabulary_scope_ref'] = null;
      var hideForLabMember = [
        'lab_head_id','lab_affiliation_id','lab_member_ids',
        'lab_head_appointment_category','lab_head_appointment_department',
        'lab_head_appointment_update_date','permissions','usergroups',
        'is_staff', 'is_superuser','is_active'];
      _.each(_.pick(resource.fields,hideForLabMember), function(field){
        field['editability'] = [];
      });
      
      var NewLabUserEditView = EditView.extend({
        save_success: function(data, textStatus, jqXHR){
          appModel.getModel('screensaveruser', self.model.key, 
            function(model){
              self.model = model;
              self.change_to_tab('detail');
              self.$('ul.nav-tabs > li').removeClass('disabled');
            });
        },
        cancel: function(e) {
          e.preventDefault();
          appModel.clearPagePending();
          this.remove();
          self.change_to_tab('detail');
          self.$el.find('#tab_container-title').empty();
          self.$('ul.nav-tabs > li').removeClass('disabled');
        }
      });
      var view = new NewLabUserEditView({ 
        model: newModel,
        resource: resource,
        uriStack: ['+add'],
      });
      this.$el.find('#tab_container-title').append('Add Lab Member...');
      Backbone.Layout.setupView(view);
      self.setView("#tab_container", view ).render();
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
      var resource = appModel.getResource('screen');
      var visibleFields =  [
          'facility_id','title','screen_type','screensaver_user_role', 
          'status','status_date','date_of_last_activity','date_created'];
      var extraIncludes = [];
      _.each(_.values(resource.fields), function(field){
        if (_.contains(visibleFields, field.key)){
          if (!_.contains(field['visibility'],'l')){
            extraIncludes.push(field.key);
          }
        } else {
          field['visibility'] = [];
        }
      });
      resource.options['includes'] = extraIncludes;
      resource.options['order'] = ['facility_id'];
      var url = [self.model.resource.apiUri, 
                 self.model.get('screensaver_user_id'),
                 'screens'].join('/');
      var Collection = Iccbl.MyCollection.extend({
        url: url
      });
      var collection = self.collection = new Collection();
      var view = new ListView({ 
          uriStack: _.clone(delegateStack),
          resource: resource,
          url: url,
          collection: collection,
          extraControls: [addRnaiScreenButton,addSMScreenButton]
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.consumedStack = ['screens']; 
      self.setView("#tab_container", view ).render();      
    },

    getTitle: function() {
      var self = this;
      var title = [
         '<H4 id="title">',
         '<a href="#' + self.model.resource.key,
         '/{screensaver_user_id}" >{first_name} {last_name} ({screensaver_user_id})</a>'];
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
    
    setActivities: function(delegateStack) {
      // Note: modified 20171108 to include all activities associated with the
      // user, either by:
      // - serviced user == user
      // - screen == user's screens
      // - performed by == user
      
      var self = this;
      var activityResource = appModel.getResource('activity');
      activityResource.fields['serviced_user_id']['editability'] = [];

      if (!_.isEmpty(delegateStack) && delegateStack[0]=='+add') {
          self.addServiceActivity();
      }else if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        
        var activity_id = delegateStack.shift();
        self.consumedStack.push(activity_id);
        var _key = activity_id
        appModel.getModelFromResource(activityResource, _key, function(model){
          view = new ServiceActivityView({
            model: model, 
            user: self.model,
            uriStack: _.clone(delegateStack)
          });
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          Backbone.Layout.setupView(view);
          self.setView("#tab_container", view ).render();
        });        
        return;
          
      }else{
        // List view
        var view, url;
        var extraControls = [];
        var addServiceActivityButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="add_button" href="#" ',
            'title="Add a service activity for the user" > ',
            'Add Service Activity</a>'
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
        
        addServiceActivityButton.click(self.addServiceActivity);
        
        // TODO: Shows histor only for service activities; could show for all
        // activities?
        showHistoryButton.click(function(e){
          e.preventDefault();
          var newUriStack = ['apilog','order','-date_time', appModel.URI_PATH_SEARCH];
          var search = {};
          search['ref_resource_name'] = 'activity';
          search['uri__contains'] = 'screensaveruser/' + self.model.get('screensaver_user_id');
          newUriStack.push(appModel.createSearchString(search));
          var route = newUriStack.join('/');
          appModel.router.navigate(route, {trigger: true});
          self.remove();
        });
        if(appModel.hasPermission(activityResource.key, 'edit')){
          extraControls.unshift(addServiceActivityButton);
          extraControls.unshift(showDeleteButton);
        }
        extraControls.unshift(showHistoryButton);
        url = [self.model.resource.apiUri, 
                   self.model.key,
                   'activities'].join('/');
        view = new ActivityListView({ 
          uriStack: _.clone(delegateStack),
          resource: activityResource,
          url: url,
          extraControls: extraControls,
          user: self.model
        });
        showDeleteButton.click(function(e){
          e.preventDefault();
          if (! view.grid.columns.findWhere({name: 'deletor'})){
            view.grid.columns.unshift({ 
              name: 'deletor', label: 'Delete', text:'X', 
              description: 'delete record', 
              sortable: false,
              cell: Iccbl.DeleteCell.extend({
                render: function(){
                  if (this.model.get('classification') == 'activity'){
                    Iccbl.DeleteCell.prototype.render.apply(this, arguments);
                  }
                  return this;
                }
              }), 
            });
          }
        });
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView("#tab_container", view ).render();
      }
    },
    
    addServiceActivity: function(e) {
      
      if (e){
        e.preventDefault();
      }
      var self = this;
      
      var activityResource = appModel.getResource('activity');
      activityResource.fields['serviced_user_id']['editability'] = [];
      activityResource.fields['serviced_username']['editability'] = [];
      activityResource.fields['serviced_user']['visibility'] = ['l','d','e'];
      
      var vocab_scope_ref = 
        activityResource.fields['classification']['vocabulary_scope_ref'];
      var vocab_classification = Iccbl.appModel.getVocabularySelectOptions(vocab_scope_ref);
      vocab_type = _.reject(vocab_classification, function(obj){
        return obj.val == 'screening';
      });
      activityResource.fields['classification'].choiceHash = vocab_type;
      
      var defaults = {
        serviced_user_id: self.model.get('screensaver_user_id'),
        serviced_user: self.model.get('name'),
      };
      var newModel = appModel.newModelFromResource(activityResource, defaults);
      var view = new ServiceActivityView({
        model: newModel,
        user: self.model,
        uriStack: ['+add']
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();

      self.consumedStack = ['activity'];
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
        resource: resource,
        url: url,
        extraControls: [uploadAttachedFileButton, showDeleteButton]
      });
      uploadAttachedFileButton.click(function(e){
        e.preventDefault();
        var type_choices = appModel.getVocabularySelectOptions('attachedfiletype.user');
        // NOTE: do not show user agreement types in generic uploads (JAS)
        type_choices = _.reject(type_choices, function(choice){
          return choice.val.match(/user_agreement/gi);
        });
        UploadDataForm.uploadAttachedFileDialog(
          url, type_choices
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
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      
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
        var newUriStack = ['apilog','order','-date_time', appModel.URI_PATH_SEARCH];
        var search = {};
        search['ref_resource_name'] = 'userchecklist';
        search['key__contains'] = self.model.get('screensaver_user_id') + '/';
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
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
            Backbone.Model.prototype.initialize.apply(this,arguments);
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
          
          if(changedCollection.isEmpty()){
            appModel.showModalError('No changes to save');
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
                fetchCollection.fetch();
              },
              done: function(){
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
          view.$el.find('td').removeClass('edited');
          appModel.clearPagePending();
        });
        
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.consumedStack = [key]; 
        self.setView("#tab_container", view ).render();
      };
        
    },
    
    getUserScreens: function() {
      var self = this;
      var allScreens = appModel.getScreens();
      var screens = allScreens.filter(function(model){
        return _.contains(self.model.get('screens'),model.get('facility_id'));
      });
      return screens;
    },

    buildMessages: function() {
      var self = this;
      
      if (!appModel.hasGroup('readEverythingAdmin')){
        return;
      }
   
      $('#content_title_message').find('#user_status_message').remove();
      var screens = self.getUserScreens();
      if (_.isEmpty(screens)) {
        return;
      }
      var error_screens_sm = [];
      var error_screens_rnai = [];
      var userDslSm = self.model.get('sm_data_sharing_level');
      var userDslRnai = self.model.get('rnai_data_sharing_level');
      _.each(screens, function(model){
        var screenDsl = model.get('data_sharing_level');
        if (model.get('screen_type')=='small_molecule'){
          if (_.isNumber(userDslSm)){
            if (userDslSm < screenDsl){
              error_screens_sm.push(
                model.get('facility_id') + ': ' + 
                appModel.getVocabularyTitle(
                  'screen.data_sharing_level',screenDsl));
            }
          }
        }else if (model.get('screen_type')=='rnai') {
          if (_.isNumber(userDslRnai)){
            if (userDslRnai < screenDsl){
              error_screens_rnai.push(
                model.get('facility_id') + ': ' + 
                appModel.getVocabularyTitle(
                  'screen.data_sharing_level',screenDsl));
            }
          }
        }else{
          console.log('uknown screen type', model.get('screen_type'), model);
        }
      });
      if (!_.isEmpty(error_screens_sm)){
        $('#content_title_message').append(
          $('<div id="user_status_message" class="alert alert-danger"></div>').html(
            "User Small Molecule Sharing (" 
            + appModel.getVocabularyTitle(
              'useragreement.data_sharing_level',userDslSm) + ')'
            + " is less restrictive than user's Screens: (" 
            + error_screens_sm.join(',') + ')'
            + ' - (violates policy)'));
      }
      if (!_.isEmpty(error_screens_rnai)){
        $('#content_title_message').append(
          $('<div id="user_status_message" class="alert alert-danger"></div>').html(
            "User RNAi Sharing (" 
            + appModel.getVocabularyTitle(
              'useragreement.data_sharing_level',userDslRnai) + ')'
            + " is less restrictive than user's Screens: (" 
            + error_screens_rnai.join(',') + ')'
            + ' - (violates policy)'));
      }
    },

    afterRender: function(){
      var self = this;
      ReportsUserView.prototype.afterRender.apply(this,arguments);
      $('#content_title_message').find('#user_status_message').remove();
      
      // check for screen dsl mismatches
      $(this).queue([
        appModel.getScreenOptions,
        self.buildMessages ]);
      
    }
    
  });

  return UserView;
});