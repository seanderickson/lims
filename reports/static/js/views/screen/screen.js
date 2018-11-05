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
  'views/activityListView',
  'views/serviceActivity',
  'views/screen/libraryScreening',
  'views/generic_detail_layout', 
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'views/list2', 
  'utils/uploadDataForm',
  'utils/tabbedController',
  'templates/generic-detail-stickit.html',
  'templates/generic-detail-stickit-one-column.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, 
            ScreenSummaryView, ScreenDataView, CherryPickRequestView, 
            ActivityListView, ServiceActivityView, LibraryScreeningView, 
            DetailLayout, DetailView, EditView, ListView, UploadDataForm, 
            TabbedController, detailTemplate, detailOneColTemplate) {

  var ScreenView = TabbedController.extend({
    isAdmin: false,
    OK_STATUSES: ['accepted','ongoing'],
    COMPLETED_STATUSES: ['completed','completed_duplicate_with_ongoing' ],
    
    initialize: function(args) {
      var self = this;
      self.args = args;
      this._classname = 'ScreenView';

      this.tabbed_resources = this.screen_tabbed_resources;
      if (! _.isEmpty(self.model.get('study_type'))) {
        this.tabbed_resources = this.study_tabbed_resources;
      } 
      
      TabbedController.prototype.initialize.apply(this,arguments);
      
      var access_level = this.model.get('user_access_level_granted');
      if (access_level > 1 && _.isEmpty(self.model.get('study_type'))) {
        this.tabbed_resources['summary'] = this.screen_tabbed_resources['summary'];
        this.tabbed_resources['protocol'] = this.screen_tabbed_resources['protocol'];
      } 
      if (access_level > 2 && _.isEmpty(self.model.get('study_type'))){
        this.tabbed_resources['cherrypickrequest'] = this.screen_tabbed_resources['cherrypickrequest'];
        this.tabbed_resources['activities'] = this.screen_tabbed_resources['activities'];
      }
      if (access_level >= 2 && self.model.get('has_screen_result') == 1){
        this.tabbed_resources['data'] = this.screen_tabbed_resources['data'];
      }
      _.bindAll(this, 'addLibraryScreening','addServiceActivity');
      
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
      protocol : {
        description : 'Screen Protocol Information',
        title : 'Protocol',
        invoke : 'setProtocol',
        permission: 'screen'
      },
      summary : {
        description : 'Summary of Library Screening transfer activities',
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
      if (this.model.get('screen_type') == 'small_molecule'
        && _.isEmpty(this.model.get('study_type'))){
        title.push('<small>A screen for compounds that... </small>');
      }
      title.push('{title}');
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
       
    /**
     * TODO: create a ScreenDetail class
     */
    setDetail: function(delegateStack) {
      
      var self,outerSelf = self = this;
      var key = 'detail';
      var resource = self.model.resource;
      var fields = resource.fields;
      var model = this.model;
      
      if (! model.isNew() 
          && !_.isEmpty(delegateStack)  && delegateStack[0] == '+add'){
        // do not allow return to +add
        console.log('cannot return to screen/id/+add...');
        delegateStack.shift();
      }
      
      // set up a custom vocabulary that joins username to name; will be 
      // used as the text of the linklist
      fields['collaborator_ids'].vocabulary = (
          _.object(model.get('collaborator_ids'),
            model.get('collaborator_names')));

      if (!_.isEmpty(model.get('primary_screen'))){
        fields['primary_screen'].visibility=['d'];
      }
      if (!_.isEmpty(model.get('reconfirmation_screens'))){
        fields['reconfirmation_screens'].visibility=['d'];
      }
      
      //// Hack to change group name:
      //_.each(fields, function(field){
      //  if (field.display_options && 
      //      field.display_options.group == 'Pin Transfer / RNAi Transfection'){
      //    if (self.model.get('screen_type')=='small_molecule'){
      //      field.display_options.group = 'Pin Transfer';
      //    } else {
      //      field.display_options.group = 'RNAi Transfection';
      //    }
      //  }
      //});
      
      // Manage visible fields; admin fields
      var editableKeys = model.resource.updateKeys();
      if (model.isNew()) {
        editableKeys = _.union(editableKeys,model.resource.createKeys());
      }
      editableKeys = _.filter(editableKeys, function(key){
        return ! (   _.contains(fields[key].visibility, 'billing')
                  || _.contains(fields[key].visibility, 'protocol'));
      });
      editableKeys = _.without(editableKeys, 'data_privacy_expiration_date');
      var editVisibleKeys = model.resource.allEditVisibleKeys();
      editVisibleKeys = _.filter(editVisibleKeys, function(key){
        return ! (   _.contains(fields[key].visibility, 'billing')
                  || _.contains(fields[key].visibility, 'protocol'));
      });
      var detailKeys = model.resource.detailKeys();
      var adminKeys = model.resource.adminKeys();
      if (!self.isAdmin){
        detailKeys = _.difference(detailKeys, adminKeys);
      }
      
      var editView = EditView.extend({

        initialize: function() {
          var self = this;
          EditView.prototype.initialize.apply(this,arguments);
        }, 
        
        save_success: function(data, textStatus, jqXHR){
          EditView.prototype.save_success.apply(this,arguments);
          appModel.unset('screens', {silent: true});
          appModel.getScreens();
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

          screenTypeSettings(model.get('screen_type'));
          // Note the forms "screen_type:change" syntax is backwards from Backbone
          form.on('screen_type:change', function(e){
            screenTypeSettings(form.getValue('screen_type'));
          });
          form.on('status:change', function(e){
            form.setValue('status_date', new Date());
          });
          
          EditView.prototype.afterRender.apply(this,arguments);
        },

        validate: function(attrs) {
          //// FIXME: 20170605 refactor: this must be tested to verify it works
          //var errs = EditView.prototype.validate.apply(this,arguments);
          //errors = errs || {};
          //console.log('extra validation...')
          //if (!_.isEmpty(this.model.get('data_privacy_expiration_notified_date'))) {
          //  if (!_.isEmpty(this.model.get('max_allowed_data_privacy_expiration_date'))) {
          //    errors['max_allowed_data_privacy_expiration_date'] = (
          //      'can not be set if the expiration notified date is set');
          //  }
          //  if (!_.isEmpty(this.model.get('min_allowed_data_privacy_expiration_date'))) {
          //    errors['min_allowed_data_privacy_expiration_date'] = (
          //      'can not be set if the expiration notified date is set');
          //  }
          //}
          //
          //if (!_.isEmpty(errors)) return errors;
          return null;
        }
      });
      
      // Custom showEdit function allows lazy loading of user choices fields
      // FIXME: it would be better to extend DetailLayout.showEdit
      var showEdit = function showEdit(updateModel) {
        var self = this;
        var model = updateModel || self.model;
        var fields = model.resource.fields;
        
        function setupEdit() {
          // FIXME: if "edit" page is reloaded, option lists here may not be loaded
          // - caused by search_box dynamic form ajax request contention
          // - result is that the choices are not populated here
          // - fix lazy load search box options; or force showEdit to block 
          // until loaded
          var userOptions = appModel.getUserOptions();
          fields['collaborator_ids'].choiceHash = userOptions;
          fields['lead_screener_id'].choiceHash = (
              [{ val: '', label: ''}].concat(userOptions));
          fields['lab_head_id'].choiceHash = (
              appModel.getPrincipalInvestigatorOptions() );
          //if (_.has(fields,'pin_transfer_approved_by_username')){
          //  fields['pin_transfer_approved_by_username'].choiceHash = 
          //    appModel.getUsernamesInGroupOptions('pinTransferApprovers');
          //}
          
          // pick just the non-billing fields: prevent backbone save from sending
          // uninitialized billing fields on create
          var modelKeys = detailKeys.concat(editVisibleKeys);
          editModel = new Backbone.Model(model.pick(modelKeys));
          editModel.set('id', model.get('id'));
          editModel.resource = model.resource;
          editModel.urlRoot = model.resource.apiUri;
          editModel.url = model.url;
          editModel.parse = model.parse;
          editModel.isNew = model.isNew;

          var editViewInstance = new editView(_.extend({}, self.args, 
            { 
              model: editModel,
              uriStack: self.uriStack,
              isCreate: true,
              editableKeys: editableKeys,
              detailKeys: detailKeys,
              editVisibleKeys: editVisibleKeys
            }));
          Backbone.Layout.setupView(editViewInstance);
          view.listenTo(editViewInstance,'remove',function(){
            view.removeView(editViewInstance);
            view.showDetail();
          });
          view.listenTo(editViewInstance, 'uriStack:change', self.reportUriStack);
          view.setView("#detail_content", editViewInstance).render();
          appModel.setPagePending();
          return editViewInstance;
        };  
        $(this).queue([
          appModel.getPrincipalInvestigatorOptions,
          appModel.getUserOptions,
          function(callback){
            appModel.getUsernamesInGroupOptions(
              'pinTransferApprovers', callback);
          }, 
          setupEdit
        ]);
      };
      
      var view = this.tabViews[key];
      if (view) {
        this.removeView(this.tabViews[key]);
      }
      var detailView = DetailView.extend({
        afterRender: function() {
          var dview = DetailView.prototype.afterRender.apply(this,arguments);

          $('#comments').empty();
          $('#comments').append(Iccbl.collapsibleText(self.model.get('comments'), 240));
          
          if (_.isEmpty(model.get('study_type'))) {
            if(self.model.get('user_access_level_granted') == 3 ){
              self.createStatusHistoryTable($('#detail_extra_information'));
              self.createActivitySummary($('#detail_extra_information'));
              self.createCprTable($('#detail_extra_information'));
            }
            if (appModel.hasGroup('readEverythingAdmin')) {
              self.createAttachedFileTable(this.$el.find('#attached_files'));
            }
          }
          
          if (_.isEmpty(self.model.get('study_type'))
              && _.isEmpty(model.get('primary_screen'))){
            
            // Reconfirmation up screens can be added, if the current screen is not a 
            // reconfirmation up screen itself.
            
            var addReconfirmationScreenControl = $([
              '<a class="btn btn-default btn-sm pull-right" ',
                'role="button" id="ScreenButton" href="#">',
                'Add a Reconfirmation Screen</a>'
              ].join(''));
            if (appModel.hasPermission('screen','write')) {
              $('#generic-detail-buttonpanel-right').append(addReconfirmationScreenControl);
            }
            addReconfirmationScreenControl.click(function(e){
              e.preventDefault();
              
              var newFacilityId = self.model.get('facility_id') + 'A';
              if (!_.isEmpty(model.get('reconfirmation_screens'))){
                var other = self.model.get('reconfirmation_screens').sort();
                var last = other[other.length-1];
                var newLetter =  String.fromCharCode(last.charCodeAt(last.length-1)+1);
                newFacilityId =  self.model.get('facility_id') + newLetter;
              }
              
              var defaults = self.model.pick([
                'screen_type', 'lab_head_id', 'lead_screener_id', 
                'data_sharing_level',
                'collaborator_ids', 'title', 'summary', 'species', 
                'cell_lines', 'transfection_agent', 
                'perturbagen_molar_concentration','perturbagen_ug_ml_concentration'
              ]);
              defaults['facility_id'] = newFacilityId;
              defaults['primary_screen'] = self.model.get('facility_id');
              var newModel = appModel.newModelFromResource(resource, defaults);
              newModel.url = function(){
                return resource.apiUri + '?override=true';
              }
              
              console.log('new reconfirmation screen:', newModel);
              self.consumedStack = ['+add'];
              editVisibleKeys.push('primary_screen');
              showEdit(newModel);
              self.$('ul.nav-tabs > li').addClass('disabled');
              self.reportUriStack([]);
              var titleDiv = $('#detail_container').find('#content_title')
              titleDiv.append(
                '<H4 id="title">Reconfirmation Up Screen for ' + 
                self.model.get('facility_id') + '</H4>');
              titleDiv.show();
            });
                                              
          }
          if (_.isEmpty(self.model.get('study_type'))
              && appModel.hasGroup('readEverythingAdmin')) {
            var adminControl = $('<a id="admin-control"></a>');
            if (self.isAdmin){
              adminControl.append('&nbsp;admin&nbsp;&lt;&lt;&nbsp;')
            }else{
              adminControl.append('&nbsp;admin&nbsp;&gt;&gt;&nbsp;')
            }
            $('#generic-detail-buttonpanel-left').append(adminControl);
            adminControl.click(function(e){
              e.preventDefault();
              outerSelf.isAdmin = !outerSelf.isAdmin;
              outerSelf.setDetail(delegateStack);
            });
          }
        },
        
        template: _.template(detailTemplate)
      });

      ScreenDetailLayout = DetailLayout.extend({
        showEdit: showEdit
      })
      // FIXME: 20170519: it would be better to pick needed values only from the args 
      view = new ScreenDetailLayout(_.extend(self.args, { 
        model: model, 
        uriStack: delegateStack,
        EditView: editView,
        editableKeys: editableKeys,
        detailKeys: detailKeys,
        editVisibleKeys: editVisibleKeys,
        DetailView: detailView
      }));
        
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
            appModel.getVocabularyTitle('activity.type.*',activity.get('type')),
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
      
      var base_fields = ['pubmed_id','pubmed_central_id','attached_filename'];
      _.each(resource.fields,function(field){
        if (_.contains(base_fields,field.key)){
          field['visibility'] = ['l','d'];
        } else {
          field['visibility'] = ['d'];
        }
      });
      var colModel = Iccbl.createBackgridColModel(resource.fields);
      var lastCol = colModel.pop();
      colModel.push({
        'name' : 'citation',
        'label' : 'Citation',
        'sortable': false,
        'searchable': false,
        'editable' : false,
        'visible': true,
        'headerCell': Backgrid.HeaderCell,
        'cell': Iccbl.LinkCell.extend({
          hrefTemplate: '#publication/{publication_id}',
          render: function(){
            var model = this.model;
            this.$el.empty();
            var formattedValue = '';
            if (this.model.has('authors')){
              formattedValue += model.get('authors') + ' '; 
            }
            if (model.has('year_published')){
              formattedValue += '(' + model.get('year_published') + '). ';
            }
            if (model.has('title')){
              var title = model.get('title').trim();
              if (title.charAt(title.length-1)!='.'){
                title += '.';
              }
              formattedValue += title + ' ';
            }
            if (model.has('journal')){
              formattedValue += model.get('journal') + ' '
            }
            if (model.has('volume')){
              formattedValue += model.get('volume') + ', ';
            }
            if (model.has('pages')){
              formattedValue += model.get('pages') + '.';
            }
            console.log('formattedValue', formattedValue);
            var interpolatedVal = Iccbl.formatString(this.hrefTemplate,model);
            this.$el.append($('<a>', {
              tabIndex : -1,
              href : interpolatedVal,
              target : this.target,
              title: 'Citation text for this publication'
            }).text(formattedValue));
            return this;
          }
        })
      });
      colModel.push(lastCol);
      var grid = new Backgrid.Grid({
        columns: colModel,
        collection: collection,
        className: 'backgrid table-striped table-condensed table-hover nested-table'
      });
      $target_el.empty();
      var cell = $('<div>',{ class: '' });
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
      resource.fields['filename']['backgridCellType'] = Iccbl.LinkCell;
      
      var collection = new CollectionClass();
      var colModel = Iccbl.createBackgridColModel(resource.fields); 
      var grid = new Backgrid.Grid({
        columns: colModel,
        collection: collection,
        className: 'backgrid table-striped table-condensed table-hover nested-table'
      });
      $target_el.empty();
      var cell = $('<div>',{ class: '' });
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
          
          var type_choices = appModel.getVocabularySelectOptions('attachedfiletype.screen');
          UploadDataForm.uploadAttachedFileDialog(
            url, type_choices
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
      
      function build_table(collection) {
        if (collection.isEmpty()) {
          return;
        }
        collection.each(function(model) {
        });
        var originalLength = collection.length;
        collection = new Backbone.Collection(collection.slice(0,10));
        
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
            'cell': Iccbl.StringCell
          })
        ];
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();

        $target_el.append($([
          '<div class="col-xs-12"><strong>',
          'Recent Cherry Pick Requests ',
          '<a href="#screen/' + self.model.get('facility_id'),
          '/cherrypickrequest">(Total: ' + originalLength + ')</a></strong></div>',
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
            limit: 0,
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
      var self = this;
      console.log('createStatusHistoryTable, ', self.model.key, self.model);
      if (!appModel.hasPermission('apilog')){
        console.log('user does not have permission to query "apilog"');
        return;
      }
      var apilogResource = appModel.getResource('apilog');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: apilogResource.apiUri 
      });
      
      function build_table(collection) {
        console.log('build status history table', collection);
        if (collection.isEmpty()) {
          return;
        }
        collection.each(function(model) {
          var diffs = model.get('diffs');
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
              'cell': Iccbl.StringCell
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
          'Add a ' + vocabulary_name + ' vocabulary item...</a>'
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
            self.listenTo(view , 'uriStack:change', self.reportUriStack);
            Backbone.Layout.setupView(view);
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
        cpResource.fields['screen_facility_id']['visibility'] = [];
        cpResource.fields['screen_type']['visibility'] = [];
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          resource: cpResource,
          url: url,
          collection: collection,
          extraControls: extraControls
        });
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
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
          ['screen_type','project_id', 'project_phase', 'lab_head_id', 
           'lab_name', 'lead_screener_id', 'lead_screener_name']),
        {
          screen_facility_id: self.model.get('facility_id'),
          screen_title: self.model.get('title'),
          date_requested: Iccbl.getISODateString(new Date()),
          requested_by_id: self.model.get('lead_screener_id')
      });
      
      if (self.model.get('screen_type') == 'small_molecule'){
        defaults['transfer_volume_per_well_requested'] = '0.0000016';
        defaults['transfer_volume_per_well_approved'] = '0.0000016';
      }else{
        defaults['transfer_volume_per_well_requested'] = null;
        defaults['transfer_volume_per_well_approved'] = null;
      }
      
      var newModel = appModel.createNewModel('cherrypickrequest', defaults);

      newModel.resource.fields['requested_by_id'].choiceHash = 
        appModel._get_screen_member_choices(self.model);
      newModel.resource.fields['volume_approved_by_username'].choiceHash = 
        appModel.getAdminUserOptions();
      
      var view = new CherryPickRequestView({ 
        model: newModel, 
        screen: self.model,
        uriStack: ['+add']
      });
      
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView("#tab_container", view ).render();

      
      this.consumedStack = ['cherrypickrequest','+add'];
      self.reportUriStack([]);
    },
    
    /**
     * 
     */
    setActivities: function(delegateStack) {
      console.log('set activities');
      var self = this;
      var activityResource = appModel.getResource('activity');
      
      _.each(_.values(activityResource.fields), function(field){
        if(_.result(field.display_options,'optgroup')=='screen'){
          field['visibility'] = [];
        }
      });
      if (_.has(activityResource.fields, 'funding_support')){
        activityResource.fields['funding_support']['visibility'] = [];
      }
      
      if (!_.isEmpty(delegateStack) && delegateStack[0]=='+add') {
          self.addServiceActivity();
      }
      else if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        // Detail View
        
        activityResource.fields['screen_facility_id'].editability = [];
        activityResource.fields['funding_support']['editability'] = [];
        
        var activity_id = delegateStack.shift();
        self.consumedStack.push(activity_id);
        var _key = activity_id
        appModel.getModelFromResource(activityResource, _key, function(model){
          view = new ServiceActivityView({
            model: model, 
            screen: self.model,
            uriStack: _.clone(delegateStack)
          });
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          Backbone.Layout.setupView(view);
          self.setView("#tab_container", view ).render();
        });        
        return;
          
      }else{
        // List view
        
        var self = this;
        var url = [self.model.resource.apiUri,self.model.key,'activities'].join('/');
        var extraControls = [];
        var addServiceActivityButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="add_button" href="#">',
            'Add Service Activity</a>'
          ].join(''));
        addServiceActivityButton.click(self.addServiceActivity);
        if(appModel.hasPermission(activityResource.key, 'edit')){
          extraControls.unshift(addServiceActivityButton);
        }
        var addLibraryScreeningButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="add_library_screening_button" href="#">',
            'Add Library Screening Visit</a>'
          ].join(''));
        addLibraryScreeningButton.click(self.addLibraryScreening);
        if (appModel.hasPermission('libraryscreening','write')) {
          extraControls.unshift(addLibraryScreeningButton);
        }

        activityResource.fields['screen_facility_id'].visibility = [];
        
        var view = new ActivityListView({ 
          uriStack: _.clone(delegateStack),
          resource: activityResource,
          url: url,
          extraControls: extraControls,
          screen: self.model
        });
        
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView("#tab_container", view ).render();

      }
      self.$("#tab_container-title").hide();
    },

    addServiceActivity: function(e){
      if (e) {
        e.preventDefault();
      }

      var self = this;
      var activityResource = appModel.getResource('activity');
      var defaults = {
        screen_facility_id: self.model.get('facility_id')
      };
      
      activityResource.fields['serviced_user']['visibility'] = [];
      activityResource.fields['serviced_user_id'].choiceHash = 
        appModel._get_screen_member_choices(self.model);
        
      activityResource.fields['serviced_username']['editability'] = [];
      activityResource.fields['screen_facility_id']['editability'] = [];
      // 20170605 - JAS - no funding support if attached to a screen
      activityResource.fields['funding_support']['editability'] = [];
      activityResource.fields['funding_support']['visibility'] = [];

      var vocab_scope_ref = 
        activityResource.fields['classification']['vocabulary_scope_ref'];
      var vocab_classification = Iccbl.appModel.getVocabularySelectOptions(vocab_scope_ref);
      vocab_type = _.reject(vocab_classification, function(obj){
        return obj.val == 'screening';
      });
      activityResource.fields['classification'].choiceHash = vocab_type;
      
      // NOTE: funding support removed from screen.service_activities: redundant
      //
      //// Funding supports: populate select with screen funding supports;
      //// - if only one screen funding support; set it as default
      //var funding_supports = self.model.get('funding_supports');
      //var funding_support_field = saResource.fields['funding_support']
      //if (_.isEmpty(funding_supports)){
      //  appModel.showModalMessage('Screen must have funding supports entered');
      //  appModel.showModalMessage({
      //    title: 'Screen must have funding supports entered',
      //    body: 'Enter screen funding supports before creating service activities'
      //  });
      //} else {
      //  var vocabulary = appModel.getVocabularySelectOptions(
      //    funding_support_field.vocabulary_scope_ref);
      //  vocabulary = _.filter(vocabulary, function(item){
      //    return _.contains(funding_supports,item.val);
      //  });
      //  funding_support_field.choiceHash = vocabulary;
      //  funding_support_field.vocabulary_scope_ref = '';
      //  if (funding_supports.length == 1){
      //    defaults['funding_support'] = funding_supports[0];
      //  }
      //}
      
      var newModel = appModel.newModelFromResource(activityResource, defaults);
      var view = new ServiceActivityView({
        model: newModel,
        screen: self.model,
        uriStack: ['+add']
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView("#tab_container", view ).render();

      self.consumedStack = ['activities'];
      self.reportUriStack([]);
      view.reportUriStack(['+add']);
    },

    addLibraryScreening: function(e) {
      e.preventDefault();
      this.change_to_tab('summary', ['libraryscreening','+add']);
//      var uriStack = ['screen', self.screen.key,
//                      'summary','libraryscreening','+add'];
//      console.log('route: ', uriStack);
//      appModel.setUriStack(uriStack);
//      
//      console.log('add library screening visit');
//      var self = this;
//      var defaults = {
//        screen_facility_id: self.model.get('facility_id'),
//        screen_type: self.model.get('screen_type')
//      };
//      var newModel = appModel.createNewModel('libraryscreening', defaults);
//
//      var view = new LibraryScreeningView({ 
//        model: newModel, 
//        screen: self.model,
//        uriStack: ['+add']
//      });
//      
//      Backbone.Layout.setupView(view);
//      self.listenTo(view , 'uriStack:change', self.reportUriStack);
//      self.setView("#tab_container", view ).render();
//
//      $title = self.$el.find('#tab_container-title');
//      $title.html('<H5 id="title">Add Library Screening Visit</H5>');
//      $title.show();
//      
//      // TODO: do not set for library screening added in the screen view?
//      this.consumedStack = ['activities'];
//      
//      self.reportUriStack([]);
//      view.reportUriStack(['+add']);
    },    
    
    
    setData: function(delegateStack) {
      var self = this;
      var key = 'data';
      view = new ScreenDataView({
        model: self.model, 
        uriStack: _.clone(delegateStack)
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView("#tab_container", view ).render();
      self.$("#tab_container-title").hide();
    },
    
    setProtocol: function(delegateStack) {
      var self = this;
      var key = 'protocol';
      var view = this.tabViews[key];
      if (!view) {
        
        var protocolKeys = self.model.resource.filterKeys('visibility', 'protocol');
        var editableKeys = _.intersection(
          self.model.resource.updateKeys(), protocolKeys);
        var editVisibleKeys = _.intersection(
          self.model.resource.allEditVisibleKeys(),protocolKeys);
        console.log('protocol keys', protocolKeys);
        var buttons = ['download','history'];
        
        if (appModel.hasPermission('screen', 'write')) {
          buttons.push('edit');
        }
        var detailView = DetailView.extend({
          template: _.template(detailOneColTemplate),
          afterRender: function() {
            var dview = DetailView.prototype.afterRender.apply(this,arguments);
            //$('#publishable_protocol').empty();
            //$('#publishable_protocol').append(
            //  Iccbl.collapsibleText(self.model.get('publishable_protocol'), 350));
          
            self.createPublicationTable(this.$el.find('#publications'));
          }
        });
        
        var ProtocolDetailLayout = DetailLayout.extend({
          showEdit: function(){
            var pSelf = this;
            appModel.getAdminUserOptions(function(options){
              pSelf.model.resource.fields['publishable_protocol_entered_by'].choiceHash = options;
              return DetailLayout.prototype.showEdit.apply(pSelf, arguments);
            });
          }
        });
        
        view = new ProtocolDetailLayout({ 
          model: self.model, 
          uriStack: delegateStack,
          detailKeys: protocolKeys,
          editVisibleKeys: editVisibleKeys,
          editableKeys: editableKeys,
          buttons: buttons,
          DetailView: detailView
        });
        self.tabViews[key] = view;
        
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
      } else {
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
      }
      self.$("#tab_container-title").hide();
    },

    setSummary: function(delegateStack) {
      var self = this;
      var key = 'summary';
      view = new ScreenSummaryView({
        model: self.model, 
        uriStack: _.clone(delegateStack)
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
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
            
            self.listenTo(view , 'uriStack:change', self.reportUriStack);
            self.setView("#tab_container", view ).render();
          },{ data_for_get: { visibilities: ['billing']}});
      } else {
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
      self.$("#tab_container-title").hide();
    },

    afterRender: function(){
      var self = this;
      TabbedController.prototype.afterRender.apply(this,arguments);
      $('#content_title_message').find('#screen_status_message').remove();
      if (!_.isEmpty(self.model.get('status'))) {
        if (_.contains(self.COMPLETED_STATUSES, self.model.get('status'))){
          $('#content_title_message').append(
            $('<div id="screen_status_message" class="alert alert-success"></div>').html(
              'Screen Status: ' + appModel.getVocabularyTitle(
                'screen.status',self.model.get('status'))));
        }
        else if (!_.contains(self.OK_STATUSES, self.model.get('status'))){
          $('#content_title_message').prepend(
            $('<div id="screen_status_message" class="alert alert-danger"></div>').html(
              'Screen Status: ' + appModel.getVocabularyTitle(
                'screen.status',self.model.get('status'))));
        }
      }
      $(this).queue([
        appModel.getUserOptions,
        self.showUserDslWarnings]);
    },
    
    showUserDslWarnings: function() {
      var self = this;
      $('#content_title_message').find('#screen_member_dsl_message').remove();
      
      if (!appModel.hasGroup('readEverythingAdmin')){
        return;
      }

      console.log('showUserDslWarnings');
      var users = appModel.getUsers();
      var screenMembers = users.filter(function(model){
        var userId = model.get('screensaver_user_id');
        return ( userId == self.model.get('lead_screener_id')
            || userId == self.model.get('lab_head_id')
            || _.contains(self.model.get('collaborator_ids'),userId));
      });
      
      var warnUsers = _.filter(screenMembers, function(member){
        var userDsl;
        if (self.model.get('screen_type')=='small_molecule'){
          userDsl = member.get('sm_data_sharing_level');
        }
        else if (self.model.get('screen_type')=='rnai'){
          userDsl = member.get('rnai_data_sharing_level');
        } else {
          console.log('unknown screen type!', self.model.get('screen_type'));
        }
        return _.isNumber(userDsl) && userDsl < self.model.get('data_sharing_level');
      });
      if (!_.isEmpty(warnUsers)){
        var userDslProp = 'sm_data_sharing_level';
        if (self.model.get('screen_type') == 'rnai'){
          userDslProp = 'rnai_data_sharing_level';
        }
        $('#content_title_message').prepend(
          $('<div id="screen_member_dsl_message" class="alert alert-danger"></div>').html(
            'Screen Data Sharing Level: ' 
            + appModel.getVocabularyTitle(
              'screen.data_sharing_level',self.model.get('data_sharing_level'))
            + ', is more restrictive than User Data Sharing Level: ' 
            + _.map(warnUsers, function(user){
              return '( #' + user.get('screensaver_user_id') + ' - '
                + user.get('name') + ' - '
                + appModel.getVocabularyTitle(
                  'useragreement.data_sharing_level',user.get(userDslProp)) + ')';
            }).join(', ')
            + ' (violates policy)'));
      }
    }
    
  });

  return ScreenView;
});