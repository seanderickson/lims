/**
 * Screen form/view
 */
define([
  'jquery', 
  'underscore', 
  'backbone', 
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/screen/screenSummary', 
  'views/screen/screenData',
  'views/screen/cherryPickRequest',
  'views/generic_detail_layout', 
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'views/list2', 
  'utils/uploadDataForm',
  'utils/tabbedController',
  'templates/generic-detail-screen.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, 
            ScreenSummaryView, ScreenDataView, CherryPickRequestView, 
            DetailLayout, DetailView, EditView, 
            ListView, UploadDataForm, TabbedController, screenTemplate) {

  var ScreenView = TabbedController.extend({
    isAdmin: false,
    OK_STATUSES: ['accepted','ongoing'],
    COMPLETED_STATUSES: ['completed','completed_duplicate_with_ongoing' ],
    
    initialize: function(args) {
      var self = this;
      self.args = args;
      this._classname = 'ScreenView';

      this.tabbed_resources = this.screen_tabbed_resources;
      if (self.model.get('project_phase') == 'annotation') {
        this.tabbed_resources = this.study_tabbed_resources;
      } 
      
      TabbedController.prototype.initialize.apply(this,arguments);
    },

    study_tabbed_resources: {
      detail: {
        description: 'Study Details',
        title: 'Study Details',
        invoke: 'setDetail'
      },
      data : {
        description : 'Data',
        title : 'Data',
        invoke : 'setData',
        permission: 'screenresult'
      }
    },
    
    screen_tabbed_resources: {
      detail : {
        description : 'Screen Details',
        title : 'Screen Details',
        invoke : 'setDetail'
      },
      summary : {
        description : 'Screening Summary',
        title : 'Screening Summary',
        invoke : 'setSummary',
        permission: 'screensummary'
      },
      data : {
        description : 'Data',
        title : 'Data',
        invoke : 'setData',
        permission: 'screenresult'
      },
      cherrypickrequest: {
        description : 'Cherry Pick Requests',
        title : 'Cherry Picks',
        invoke : 'setCherryPicks',
        permission: 'cherrypickrequest'
      },
      activities: {
        description : 'Activities',
        title : 'Activities',
        invoke : 'setActivities',
        permission: 'activity'
      },
      billingItems: {
        description : 'Billing information',
        title : 'Billing',
        invoke : 'setBilling',
        permission: 'screenbilling'
      },
    },

    render : function() {
      window.scrollTo(0, 0);
      return this;
    },

    getTitle: function() {
      var self = this;
      var title = [
         '<H4 id="title">',
         '<a href="#screen/{facility_id}" >{facility_id}</a>: '];
      if (this.model.get('screen_type') == 'small_molecule'){
        title.push('<small>A screen for compounds that... </small>');
      }
      title.push('{title}');
      title.push('</H4>');
      var titleDiv = $(Iccbl.formatString(title.join(''), this.model));
      if (appModel.hasGroup('readEverythingAdmin')){
        titleDiv.append(
          Iccbl.createCommentIcon([self.model.get('comments')],'Commments'));
      }
      return titleDiv;
    },
        
    setDetail: function(delegateStack) {
      var self,outerSelf = self = this;
      var key = 'detail';
      var fields = self.model.resource.fields;
      // set up a custom vocabulary that joins username to name; will be 
      // used as the text of the linklist
      if (this.model.get('project_phase')!= 'annotation') {
        fields['collaborator_usernames'].vocabulary = (
            _.object(this.model.get('collaborator_usernames'),
              this.model.get('collaborator_names')));
      }
      var editView = EditView.extend({
        save_success: function(data, textStatus, jqXHR){
          EditView.prototype.save_success.apply(this,arguments);
          appModel.unset('users');
        },
        afterRender: function() {
          var form = this;
          outerSelf._addVocabularyButton(
            this, 'cell_lines', 'cell_line', 'Cell Line', { description: 'ATCC Designation' });
         
          outerSelf._addVocabularyButton(
            this, 'transfection_agent', 'transfection_agent', 'Transfection Agent');
         
          // 2016-11-15 removed per JAS
          //outerSelf._addVocabularyButton(
          //  this, 'species', 'screen.species', 'Screened Species');

          form.$el.find('div[key="title"]').parent().prepend(
            '<span id="title-sm-screen">A screen for compounds that...</span>');
          function screenTypeSettings(screen_type){
            // 20161128 - adjust visible fields based on the screen type
            if (screen_type == 'small_molecule'){
              form.$el.find('#title-sm-screen').show();
              form.$el.find('div[data-fields="transfection_agent"]').hide();
              form.$el.find('div[data-fields="perturbagen_ug_ml_concentration"]').hide();
              form.$el.find('div[data-fields="perturbagen_molar_concentration"]').hide();
            }else{
              form.$el.find('#title-sm-screen').hide();
              form.$el.find('div[data-fields="transfection_agent"]').show();
              form.$el.find('div[data-fields="perturbagen_ug_ml_concentration"]').hide();
              form.$el.find('div[data-fields="perturbagen_molar_concentration"]').show();
            }
          };

          screenTypeSettings(this.model.get('screen_type'));
          // Note the "screen_type:change" syntax is backwards from Backbone
          form.on('screen_type:change', function(e){
            screenTypeSettings(form.getValue('screen_type'));
          });
          form.on('status:change', function(e){
            form.setValue('status_date', new Date());
          });
          EditView.prototype.afterRender.apply(this,arguments);
        }
      });
      
      var view = this.tabViews[key];
      if (view) {
        this.removeView(this.tabViews[key]);
      }
      var detailView = DetailView.extend({
        afterRender: function() {
          var dview = DetailView.prototype.afterRender.apply(this,arguments);
          
          if (self.model.get('project_phase')=='annotation') {
            // do nothing
          } else {
            self.createStatusHistoryTable($('#status_table'));
            self.createActivitySummary($('#activity_summary'));
            self.createCprTable($('#cpr_table'));
            self.createPublicationTable(this.$el.find('#publications'));
            self.createAttachedFileTable(this.$el.find('#attached_files'));
            
          }
          
          if (appModel.hasGroup('readEverythingAdmin')) {
            var adminControl = $('<a id="admin-control"></a>');
            if (self.isAdmin){
              adminControl.append('admin&nbsp;&lt;&lt;')
            }else{
              adminControl.append('admin&nbsp;&gt;&gt;')
            }
            $('#generic-detail-buttonpanel').append(adminControl);
            adminControl.click(function(e){
              e.preventDefault();
              self.isAdmin = !self.isAdmin;
              self.render();
            });
          }
        },
        
        serialize: function() {
          // special handling of the grouped keys for the screen template

          var data = DetailView.prototype.serialize.apply(this,arguments);
          var informationKeys = [];
          var groupedKeys = [];
          data['groupedKeys'].each(function(groupKey) {
            if (_.result(groupKey,'title',null) == 'Information') {
              informationKeys = informationKeys.concat(groupKey.fields);
            } else {
              groupedKeys.push(groupKey);
            }
          });
          data['groupedKeys'] = _.chain(groupedKeys);
          data['informationKeys'] = _.chain(informationKeys);
          return data;
        },
        
        template: _.template(screenTemplate)
      });

      var temp_validate = this.model.validate;
      this.model.validate = function(attrs, options) {
        errs = temp_validate();
        if (!_.isEmpty(_.result(attrs,'data_privacy_expiration_notified_date'))) {
          if (!_.isEmpty(_.result(attrs,'max_allowed_data_privacy_expiration_date'))) {
            errs['max_allowed_data_privacy_expiration_date'] = (
              'can not be set if the expiration notified date is set');
          }
          if (!_.isEmpty(_.result(attrs,'min_allowed_data_privacy_expiration_date'))) {
            errs['min_allowed_data_privacy_expiration_date'] = (
              'can not be set if the expiration notified date is set');
          }
        }
        return _.isEmpty(errs) ? null : errs;
      };
        
      var editableKeys = self.model.resource.updateKeys();
      if (self.model.isNew()) {
        editableKeys = _.union(editableKeys,self.model.resource.createKeys());
      }
      editableKeys = _.filter(editableKeys, function(key){
        return ! _.contains(fields[key].visibility, 'billing');
      });
      var editVisibleKeys = self.model.resource.allEditVisibleKeys();
      editVisibleKeys = _.filter(editVisibleKeys, function(key){
        return ! _.contains(fields[key].visibility, 'billing');
      });
      var detailKeys = self.model.resource.detailKeys();
      var adminKeys = self.model.resource.adminKeys();
      console.log('adminKeys', adminKeys);
      if (!self.isAdmin){
        detailKeys = _.difference(detailKeys, adminKeys);
      }
      
      /////
      // pick just the non-billing fields: prevent backbone save from sending
      // uninitialized billing fields on create
      var modelKeys = detailKeys.concat(editVisibleKeys);
      model = new Backbone.Model(this.model.pick(modelKeys));
      model.set('id', this.model.get('id'));
      model.resource = this.model.resource;
      model.urlRoot = this.model.resource.apiUri;
      /////
      
      view = new DetailLayout(_.extend(self.args, { 
        model: model, 
        uriStack: delegateStack,
        EditView: editView,
        editableKeys: editableKeys,
        detailKeys: detailKeys,
        editVisibleKeys: editVisibleKeys,
        DetailView: detailView
      }));
      view.showEdit = function() {
        appModel.initializeAdminMode(function() {
          var userOptions = appModel.getUserOptions();
          fields['collaborator_usernames']['choices'] = userOptions;
          fields['lead_screener_username']['choices'] = (
              [{ val: '', label: ''}].concat(userOptions));
          fields['lab_head_username']['choices'] = (
              appModel.getPrincipalInvestigatorOptions() );
          fields['publishable_protocol_entered_by']['choices'] = (
              appModel.getAdminUserOptions() );
          fields['pin_transfer_approved_by_username']['choices'] = (
              appModel.getAdminUserOptions() );
          // TODO: this should default only after pp value is entered
          //view.model.set('publishable_protocol_entered_by',
          //    appModel.getCurrentUser().username);
          DetailLayout.prototype.showEdit.apply(view,arguments);
        });  
      };
      this.$("#tab_container-title").html("");

      this.tabViews[key] = view;
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.consumedStack = []; 
      this.setView("#tab_container", view ).render();
      
      return view;
    },
    
    createActivitySummary: function($target_el) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'activities'].join('/');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: url
      });
      $target_el.empty();
      
      function build_table(collection) {
          if (collection.isEmpty()) {
            return;
          }
          var activity = collection.at(0);
          
          $target_el.append($([
            '<div class="col-xs-12"><strong>Activity Summary</strong></div>',
            '<div id="" class="col-xs-12" >',
            '<table id="activity_summary_table" class="table-condensed data-list">',
            '<tr>',
            '<td class="dl-title small">Activities</td>', 
            '<td class="dl-data small">',
            Iccbl.formatString(
              '<a href="#screen/{facility_id}/activities">{activity_count}</a>', 
              self.model),
            '</td>', 
            '</tr>',
            '<tr>',
            '<td class="dl-title small">Last Activity</td>', 
            '<td class="dl-data small">',
            '<table id="last_activity" class="table-condensed data-list">',
            '<tr>',
            '<td class="dl-title ">Date</td>', 
            '<td class="dl-data ">' + activity.get('date_of_activity') + '</td>',
            '</tr>',
            '<tr>',
            '<td class="dl-title ">Activity Type</td>', 
            '<td class="dl-data ">',
            appModel.getVocabularyTitle('activity.type',activity.get('type')),
            '</td>',
            '</tr>',
            '<tr>',
            '<td class="dl-title col-xs-6">Performed By</td>', 
            '<td class="dl-data ">' + activity.get('performed_by_name') + '</td>',
            '</tr>',
            '</table>',
            '</td>',
            '</tr>',
            '</table>',
            '</div>'
            ].join('')));

          $('#activity_count').closest('tr').remove();
      }      
      if (self.model.has('latest_activities_data')) {
        build_table(new CollectionClass(self.model.get('latest_activities_data')));
      } else {
        var cpr_collection = new CollectionClass();
        cpr_collection.fetch({
          data: { 
            limit: 1,
            order_by: ['-date_of_activity']
          },
          success: build_table
        }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      }
    },

    addPublicationDialog: function(collection) {
      var self = this;
      var url = [self.model.resource.apiUri, 
         self.model.key,
         'publications'].join('/');
      var form_template = [
         "<div class='form-horizontal container' id='publication_form' >",
         "<form data-fieldsets class='form-horizontal container' >",
         //"<div class='form-group' ><input type='file' name='fileInput' /></div>",
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
      var DisabledField = EditView.DisabledField.extend({
        tagName: 'p',
        className: 'form-control-static'
      });
      
      var formSchema = {};
      formSchema['pubmed_id'] = {
        title: 'Pubmed ID',
        key: 'pubmed_id',
        type: 'Text',
        template: fieldTemplate
      };
      formSchema['pubmed_central_id'] = {
        title: 'Pubmed Central ID',
        key: 'pubmed_central_id',
        type: 'Text',
        template: fieldTemplate
      };
      formSchema['title'] = {
        title: 'Title',
        key: 'title',
        type: 'TextArea',
        editorClass: 'input-full',
        validators: ['required'],
        template: fieldTemplate
      };
      formSchema['authors'] = {
        title: 'Authors',
        key: 'authors',
        type: 'TextArea',
        editorClass: 'input-full',
        validators: ['required'],
        template: fieldTemplate
      };
      formSchema['journal'] = {
        title: 'Journal',
        key: 'journal',
        type: 'TextArea',
        editorClass: 'input-full',
        validators: ['required'],
        template: fieldTemplate
      };
      formSchema['volume'] = {
        title: 'Volume',
        key: 'volume',
        type: 'Text',
        template: fieldTemplate
      };
      formSchema['year_published'] = {
        title: 'Year Published',
        key: 'year_published',
        type: 'Text',
        template: fieldTemplate
      };
      formSchema['pages'] = {
        title: 'Pages',
        key: 'pages',
        type: 'Text',
        template: fieldTemplate
      };
      formSchema['fileInputPlaceholder'] = {
        title: 'Upload File',
        key: 'file_input',
        type: DisabledField, 
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
        validate: function(attrs) {
          var errs = {};
          if (_.isEmpty(attrs.pubmed_central_id) 
              && _.isEmpty(attrs.pubmed_id)) {
            var msg = 'must specify either "Pubmed ID" or "Pubmed CID"';
            errs['pubmed_id'] = msg;
            errs['pubmed_central_id'] = msg;
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
        okText: 'ok',
        ok: function(e) {
          e.preventDefault();
          var errors = form.commit({ validate: true }); // runs schema and model validation
          if (!_.isEmpty(errors)) {
            console.log('form errors, abort submit: ',errors);
            _.each(_.keys(errors), function(key) {
              $('[name="'+key +'"]').parents('.form-group').addClass('has-error');
            });
            return false;
          } else {
            var values = form.getValue();
            var comments = _.result(values,'comments','');
            var publication_id = _.result(values,'pubmed_id', _.result(values,'pubmed_central_id'));
            var headers = {};
            headers[appModel.HEADER_APILOG_COMMENT] = comments;
            
            var data = new FormData();
            _.each(_.keys(values), function(key) {
              if (values[key]) {
                data.append(key,values[key]);
              }
            });
            
            var file = $('input[name="fileInput"]')[0].files[0];
            var filename;
            if (file) {
              data.append('attached_file',file);
              filename = file.name;
              if (!_.isEmpty(values['filename'])) {
                filename = values['filename'];
              }
              if (_.isEmpty(values['file_date'])) {
                data.append(
                  'file_date',
                  Iccbl.getISODateString(_.result(file,'lastModifiedDate', null)));
              }
            } else {
              filename = values['filename'];
            }
            
            $.ajax({
              url: url,    
              data: data,
              cache: false,
              contentType: false,
              processData: false,
              dataType: 'json', // what is expected back from the server
              type: 'POST',
              headers: headers
            }).fail(function() {
              appModel.jqXHRfail.apply(this,arguments); 
            }).done(function(data, textStatus, jqXHR) {
              collection.fetch({ reset: true });
              var msg = 'Success';
              if (data) {
                msg = data['title'];
              }
              appModel.showModalMessage({
                title: 'Publication added',
                okText: 'ok',
                body: msg
              });
            });
          
            return true;
          }
        },
        view: _form_el,
        title: 'Upload an Attached File'  
      });
    },
    
    createPublicationTable: function($target_el) {
      var self = this;
      var key = 'publication';
      var resource = appModel.getResource(key);
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'publications'].join('/');
      
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: url
      });
      
      var collection = new CollectionClass();
      var colModel = Iccbl.createBackgridColModel(resource.fields); 
      var grid = new Backgrid.Grid({
        columns: colModel,
        collection: collection,
        className: 'backgrid table-striped table-condensed table-hover'
      });
      $target_el.empty();
      var cell = $('<div>',{ class: 'col-xs-4' });
      var grid_el = grid.render().$el;
      cell.html(grid_el);

      collection.on('all', function(){
        if (collection.isEmpty()) {
          grid_el.hide();
        } else {
          grid_el.show();
        }
      });
      if (appModel.hasPermission('screen','write')) {
        var showDeleteButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
              'role="button" id="showDeleteButton" href="#">',
              'Delete</a>'
            ].join(''));
        
        var showAddButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
              'role="button" id="showAddButton" href="#">',
              'Add</a>'
            ].join(''));
        
        cell.append(showAddButton);
        cell.append(showDeleteButton);
        showDeleteButton.click(function(e) {
          e.preventDefault();
          if (! grid.columns.findWhere({name: 'deletor'})) {
            grid.columns.unshift({ 
              name: 'deletor', label: 'Delete', text:'X', 
              description: 'delete record', 
              cell: Iccbl.DeleteCell, sortable: false });
          }
        });
        self.listenTo(collection, "MyCollection:delete", function (model) {
          
          var title = 'Confirm deletion of publication: ' +
            Iccbl.getTitleFromTitleAttribute(model, resource);
          appModel.showOkCommentForm({ 
            title: title, 
            ok: function(values){
              appModel.clearPagePending();
              var headers = {};
              headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
              
              model.collection = collection;
              var patchUrl = [resource.apiUri, 
                         Iccbl.getIdFromIdAttribute(model, resource),
                         ].join('/');
              model.url = patchUrl;
              
              // Backbone will only send DELETE if the model has an id
              model.set('id', Iccbl.getIdFromIdAttribute(model,resource));
              model.destroy({
                wait: true,
                headers: headers,
                success: function(model,response){
                  console.log('model removed successfully', model, response);
                }
              }).fail(function(){ appModel.jqXHRfail.apply(this,arguments); });      
            }
          });
        });            
        
        showAddButton.click(function(e) {
          e.preventDefault();
          self.addPublicationDialog(collection);
        });
      } 
      $target_el.append(cell);
      
      collection.fetch({
        data: { 
          limit: 0
        },
        success: function(collection, response) {
          console.log('publications fetched.');
        }
      }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });
    },

    createAttachedFileTable: function($target_el) {
      var self = this;
      var key = 'attachedfile';
      var resource = appModel.getResource(key);
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'attachedfiles'].join('/');
      
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: url
      });
      
      var collection = new CollectionClass();
      var colModel = Iccbl.createBackgridColModel(resource.fields); 
      var grid = new Backgrid.Grid({
        columns: colModel,
        collection: collection,
        className: 'backgrid table-striped table-condensed table-hover'
      });
      $target_el.empty();
      var cell = $('<div>',{ class: 'col-xs-4' });
      var grid_el = grid.render().$el;
      cell.html(grid_el);

      collection.on('all', function(){
        if (collection.isEmpty()) {
          grid_el.hide();
        } else {
          grid_el.show();
        }
      });
      if (appModel.hasPermission('screen','write')) {
        var showDeleteButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
              'role="button" id="showDeleteButton" href="#">',
              'Delete</a>'
            ].join(''));
        
        var showAddButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
              'role="button" id="showAddButton" href="#">',
              'Add</a>'
            ].join(''));
        
        cell.append(showAddButton);
        cell.append(showDeleteButton);
        showDeleteButton.click(function(e) {
          e.preventDefault();
          if (! grid.columns.findWhere({name: 'deletor'})) {
            grid.columns.unshift({ 
              name: 'deletor', label: 'Delete', text:'X', 
              description: 'delete record', 
              cell: Iccbl.DeleteCell, sortable: false });
          }
        });
        self.listenTo(collection, "MyCollection:delete", function (model) {
          
          var title = 'Confirm deletion of attached file: ' + 
            Iccbl.getTitleFromTitleAttribute(model, resource);
          appModel.showOkCommentForm({
            title: title, 
            ok: function(values){
              appModel.clearPagePending();
              var headers = {};
              headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
              
              model.collection = collection;
              var patchUrl = [resource.apiUri, 
                         Iccbl.getIdFromIdAttribute(model, resource),
                         ].join('/');
              model.url = patchUrl;
              
              // Backbone will only send DELETE if the model has an id
              model.set('id', Iccbl.getIdFromIdAttribute(model,resource));
              model.destroy({
                wait: true,
                headers: headers,
                success: function(model,response){
                  console.log('model removed successfully', model, response);
                }
              }).fail(function(){ appModel.jqXHRfail.apply(this,arguments); });
            }
          });
        });            
        
        showAddButton.click(function(e) {
          e.preventDefault();
          UploadDataForm.uploadAttachedFileDialog(
            url, 'attachedfiletype.screen'
          ).done(function(data, textStatus, jqXHR){
            var msg = 'Success';
            if (data) {
              msg = data['filename'];
            }
            appModel.showModalMessage({
              title: 'Attached File uploaded',
              okText: 'ok',
              body: msg
            });
            collection.fetch({ reset: true });
          }).fail(function(){
            appModel.jqXHRfail.apply(this,arguments); 
          });
          
        });
      } 
      $target_el.append(cell);
      
      collection.fetch({
        data: { 
          limit: 0
        },
        success: function(collection, response) {
          console.log('attached files fetched.');
        }
      }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });
    },
    
    createCprTable: function($target_el) {
      var self = this;
      var cprResource = appModel.getResource('cherrypickrequest');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: cprResource.apiUri 
      });
      $target_el.empty();
      
      function build_table(collection) {
        if (collection.isEmpty()) {
          return;
        }
        collection.each(function(model) {
        });
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell
        };
        var columns = [
          _.extend({},colTemplate,{
            'name' : 'cherry_pick_request_id',
            'label' : '#',
            'description' : 'Cherry Pick Request ID',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.LinkCell.extend({
              hrefTemplate: 
                '#screen/{screen_facility_id}/cherrypickrequest/{cherry_pick_request_id}'
            })
          }),
          _.extend({},colTemplate,{
            'name' : 'date_requested',
            'label' : 'Date Requested',
            'description' : 'Date',
            'order': 1,
            'sortable': true,
            'cell': 'Date'
          }),
          _.extend({},colTemplate,{
            'name' : 'requested_by_name',
            'label' : 'Requested By',
            'description' : 'Requested By',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.TextWrapCell.extend({
              className: 'text-wrap-cell-extra-narrow'
            })
          })
        ];
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();

        $target_el.append($([
          '<div class="col-xs-12"><strong>',
          'Recent <a href="#screen/' + self.model.get('facility_id'),
          '/cherrypick">Cherry Pick Requests</a></strong></div>',
          '<div class="col-xs-12" id="cprs"/>'].join('')));
        
        var _grid = new Backgrid.Grid({
          columns: colModel,
          collection: collection,
          className: 'backgrid table-striped table-condensed table-hover '
        });
        $target_el.find('#cprs').html(_grid.render().$el);
        
      }
      
      if (self.model.has('cherry_pick_request_data')) {
        build_table(new CollectionClass(
          self.model.get('cherry_pick_request_data')));
      } else {
        var cpr_collection = new CollectionClass();
        cpr_collection.fetch({
          data: { 
            limit: 10,
            screen_facility_id__eq: self.model.key,
            order_by: ['-date_requested']
          },
          success: build_table
        }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      }
    },
    
    /**
     * Update the screen status with a status history table: populate
     * using the apilog history of the status attribute
     **/
    createStatusHistoryTable: function($target_el) {
      console.log('createStatusHistoryTable');
      var self = this;
      var apilogResource = appModel.getResource('apilog');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: apilogResource.apiUri 
      });
      $target_el.empty();
      
      function build_table(collection) {
        console.log('build status history table', collection);
        if (collection.isEmpty()) {
          return;
        }
        collection.each(function(model) {
          var diffs = JSON.parse(model.get('diffs'));
          if (!_.isEmpty(diffs.status[1])){
            model.set(
              'status', 
              appModel.getVocabularyTitle('screen.status', diffs.status[1]));
          }else{
            collection.remove(model);
          }
        });
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell
        };
        var columns = [
            _.extend({},colTemplate,{
              'name' : 'status',
              'label' : 'Status',
              'description' : 'Screen status',
              'order': 1,
              'sortable': true,
              'cell': Iccbl.TextWrapCell.extend({
                className: 'text-wrap-cell-extra-narrow'
              })
            }),
            _.extend({},colTemplate,{
              'name' : 'date_time',
              'label' : 'Date',
              'description' : 'Date',
              'order': 1,
              'sortable': true,
              'cell': 'Date'
            })];
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();

        $('#status').closest('tr').remove();
        $('#status_date').closest('tr').remove();
        $target_el.append($([
          '<div class="col-xs-12"><strong>Status Items</strong></div>',
          '<div class="col-xs-12" id="status_items"/>'].join('')));
        
        var status_grid = new Backgrid.Grid({
          columns: colModel,
          collection: collection,
          className: 'backgrid table-striped table-condensed table-hover '
        });
        $target_el.find('#status_items').html(status_grid.render().$el);
        
      }
      
      if (self.model.has('status_data')) {
        build_table(new CollectionClass(self.model.get('status_data')));
      } else {
        var status_collection = new CollectionClass();
        status_collection.fetch({
          data: { 
            limit: 0,
            key: self.model.get('facility_id'),
            ref_resource_name: self.model.resource.key,
            diff_keys__icontains: '"status"',
            order_by: ['-date_time']
          },
          success: build_table
        }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      }
      
    },
    
    _addVocabularyButton: function(
        editForm, fieldKey, vocabulary_scope_ref, vocabulary_name, options) {
      var options = options || {};
      var addButton = $([
        '<a class="btn btn-default btn-sm" ',
          'role="button" id="add_' + fieldKey + '_button" href="#">',
          'Add a vocabulary item...</a>'
        ].join(''));
      addButton.click(function(event) {
        event.preventDefault();
        appModel.addVocabularyItemDialog(
          vocabulary_scope_ref, vocabulary_name,
          function(new_vocab_item) {
            // manually add the new vocabulary options to the multiselect
            editForm.$el.find('[key="' + fieldKey + '"]')
              .find('.chosen-select').append($('<option>',{
                value: new_vocab_item['key']
              }).text(new_vocab_item['title']));
            editForm.$el.find('[key="' + fieldKey + '"]')
              .find('.chosen-select').trigger("chosen:updated");
          }, options);
      });
      editForm.$el.find('div[key="form-group-' + fieldKey + '"]').append(addButton);
    },
    
    setCherryPicks: function(delegateStack) {
      var self = this;
      var cpResource = appModel.getResource('cherrypickrequest'); 
      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        // Detail view
        var cherryPickRequestId = delegateStack.shift();
        if (cherryPickRequestId == '+add'){
           appModel.initializeAdminMode(function() {
             self.showAddCherryPick();
           });
        } else {
          self.consumedStack.push(cherryPickRequestId);
          var _key = cherryPickRequestId;
          appModel.getModel(cpResource.key, _key, function(model){
            view = new CherryPickRequestView({
              model: model, 
              uriStack: _.clone(delegateStack),
              screen: self.model
            });
            Backbone.Layout.setupView(view);
            self.listenTo(view , 'uriStack:change', self.reportUriStack);
            self.setView("#tab_container", view ).render();
            
            console.log('title: ', 
              Iccbl.getTitleFromTitleAttribute(model, model.resource));
            self.$("#tab_container-title").html(view.getTitle());
            self.$("#tab_container-title").show();
          });        
        }
        return;
      }else{
        // List view
        var extraControls = [];
        var url = [self.model.resource.apiUri, 
                   self.model.key,
                   'cherrypicks'].join('/');
        var Collection = Iccbl.MyCollection.extend({
          url: url
        });
        collection = new Collection();
        if (appModel.hasPermission(cpResource.key, 'write')){
          var showAddButton = $([
             '<a class="btn btn-default btn-sm pull-down" ',
               'role="button" id="add_resource" href="#">',
               'Add</a>'
             ].join(''));   
           showAddButton.click(function(e){
             e.preventDefault();

             appModel.initializeAdminMode(function() {
               self.showAddCherryPick();
             });
           });
           extraControls.push(showAddButton);
        }    
        
        cpResource.fields['cherry_pick_request_id'].backgridCellType = 
          Iccbl.LinkCell.extend(_.extend({},
            cpResource.fields['cherry_pick_request_id'].display_options,
            {
              linkCallback: function(e){
                e.preventDefault();
                // re-fetch the full model
                var _key = this.model.get('cherry_pick_request_id');
                appModel.getModel(cpResource.key, _key, function(model){
                  view = new CherryPickRequestView({
                    model: model, 
                    uriStack: [],
                    screen: self.model
                  });
                  Backbone.Layout.setupView(view);
                  self.consumedStack = ['cherrypickrequest',_key];
                  self.listenTo(view , 'uriStack:change', self.reportUriStack);
                  self.setView("#tab_container", view ).render();
                  self.$("#tab_container-title").html(view.getTitle());
                  self.$("#tab_container-title").show();
                });        
                return;
              }
            }));
        
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: cpResource,
          resource: cpResource,
          url: url,
          collection: collection,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.$("#tab_container-title").empty();
        self.$("#tab_container-title").hide();
        self.setView("#tab_container", view ).render();
      }
    },
    
    showAddCherryPick: function() {
      console.log('add cherry pick request');
      var self = this;

      var defaults = _.extend(
        {}, 
        _.pick(
          self.model.attributes, 
          ['screen_type','project_id', 'project_phase', 'lab_head_username', 
           'lab_name', 'lead_screener_username', 'lead_screener_name']),
        {
          screen_facility_id: self.model.get('facility_id'),
          screen_title: self.model.get('title'),
          date_requested: Iccbl.getISODateString(new Date()),
          requested_by_username: self.model.get('lead_screener_username')
      });
      
      if (self.model.get('screen_type') == 'small_molecule'){
        defaults['transfer_volume_per_well_requested'] = '0.0000016';
        defaults['transfer_volume_per_well_approved'] = '0.0000016';
      }
      
      var newModel = appModel.createNewModel('cherrypickrequest', defaults);

      newModel.resource.fields['requested_by_username'].choices = 
        appModel._get_screen_members(self.model);
      newModel.resource.fields['volume_approved_by_username'].choices = 
        appModel.getAdminUserOptions();
      
      var view = new CherryPickRequestView({ 
        model: newModel, 
        screen: self.model,
        uriStack: ['+add']
      });
      
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();

      
      this.consumedStack = ['cherrypickrequest','+add'];
      self.reportUriStack([]);
    },
    
    setActivities: function(delegateStack) {
      var self = this;
      var key = 'activities';
      var view = this.tabViews[key];
      
      if (!view) {
        var self = this;
        var url = [self.model.resource.apiUri,self.model.key,'activities'].join('/');
        var resource = appModel.getResource('activity');
        
        var sa_vocab = appModel.getVocabulary('activity.type');
        resource.fields['type'].vocabulary = 
          _.map(sa_vocab, function(v) {
            return [v.title,v.key];
          }); // TODO: app model method for this
        
        console.log('combined vocab', resource.fields['type'].vocabulary );
        
        var view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: resource,
          resource: resource,
          url: url,
          extraControls: []
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        this.$('li').removeClass('active');
        this.$('#'+key).addClass('active');

      } else {
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
      }
      self.$("#tab_container-title").hide();
    },
    
    setData: function(delegateStack) {
      var self = this;
      var key = 'data';
      view = new ScreenDataView({
        model: self.model, 
        uriStack: _.clone(delegateStack)
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.$("#tab_container-title").hide();
    },

    setSummary: function(delegateStack) {
      var self = this;
      var key = 'summary';
      view = new ScreenSummaryView({
        model: self.model, 
        uriStack: _.clone(delegateStack)
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.$("#tab_container-title").hide();
    },

    setBilling: function(delegateStack) {
      var self = this;
      var key = 'billing';
      var view = this.tabViews[key];
      if (!view) {
        
        var billingKeys = self.model.resource.filterKeys('visibility', 'billing');
        var editableKeys = _.intersection(
          self.model.resource.updateKeys(), billingKeys);
        var editVisibleKeys = _.intersection(
          self.model.resource.allEditVisibleKeys(),billingKeys);

        var buttons = ['download','history'];
        
        if (appModel.hasPermission('screen', 'write')) {
          buttons.push('edit');
        }
        
        appModel.getModel(
          self.model.resource.key, self.model.key, 
          function(model) {
            view = new DetailLayout({ 
              model: model, 
              uriStack: delegateStack,
              detailKeys: billingKeys,
              editVisibleKeys: editVisibleKeys,
              editableKeys: editableKeys,
              buttons: buttons
            });
            self.tabViews[key] = view;
            
            self.listenTo(view , 'uriStack:change', this.reportUriStack);
            self.setView("#tab_container", view ).render();
          },{ data_for_get: { visibilities: ['billing']}});
      } else {
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
      }
      self.$("#tab_container-title").hide();
    },

    /**
     * Library Copy Plates Loaded view is a sub-view of Summary
     */
    showCopyPlatesLoaded: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'copyplatesloaded'].join('/');
      var resource = appModel.getResource('librarycopyplate');
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: []
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Library Copy Plates loaded for Screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');
      self.$("#tab_container-title").hide();
    },
    
    afterRender: function(){
      var self = this;
      TabbedController.prototype.afterRender.apply(this,arguments);
      $('#content_title_message').empty();
      if (self.model.has('status')) {
        if (_.contains(self.COMPLETED_STATUSES, self.model.get('status'))){
          $('#content_title_message').append($('<div class="alert alert-success"></div>').html(
            'Screen Status: ' + appModel.getVocabularyTitle(
              'screen.status',self.model.get('status'))));
        }
        else if (!_.contains(self.OK_STATUSES, self.model.get('status'))){
          $('#content_title_message').prepend($('<div class="alert alert-danger"></div>').html(
            'Screen Status: ' + appModel.getVocabularyTitle(
              'screen.status',self.model.get('status'))));
        }
      }
    }

  });

  return ScreenView;
});