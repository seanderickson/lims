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
  'views/screen/libraryScreening', 
  'views/generic_detail_layout', 
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'views/list2', 
  'views/collectionColumns',
  'templates/generic-tabbed.html',
  'templates/generic-detail-screen.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, 
            LibraryScreening, DetailLayout, DetailView, EditView, ListView, 
            CollectionColumnView, tabbedTemplate, screenTemplate){

  var ScreenView = Backbone.Layout.extend({

    initialize: function(args) {
      var self = this;
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.tabbed_resources = _.extend(
        {}, _.mapObject(this.tabbed_resources_template, function(val,key){
          return _.clone(val);
        }));
      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail'){
          var permission = self.tabbed_resources[key].permission;
          if (_.isUndefined(permission)){
            permission = self.tabbed_resources[key].resource;
          }
          if (!appModel.hasPermission(permission)){
            delete self.tabbed_resources[key];
          }
        }
      });
      if (self.model.get('project_phase') == 'annotation'){
        delete self.tabbed_resources['summary'];
        delete self.tabbed_resources['billingItems'];
        delete self.tabbed_resources['cherrypicks'];
        delete self.tabbed_resources['activities'];
        self.tabbed_resources['detail'] = self.study_tab;
        self.tabbed_resources['results'].title = 'Reagents';
        self.tabbed_resources['results'].description = 'Reagents studied';
      }

      
      self.$loadScreenResultsButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
        id="loadScreenResults" >Load Data</button>');
      self.$deleteScreenResultsButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
        id="deleteScreenResults" >Delete Data</button>');

      _.bindAll(this, 'click_tab');
    },

    template: _.template(tabbedTemplate),
    study_tab: {
      description: 'Study Details',
      title: 'Study Details',
      invoke: 'setDetail'
    },
    tabbed_resources_template: {
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
      attachedfile: {
        description: "Attached Files",
        title: "Attached Files",
        invoke: "setAttachedFiles",
        resource: 'attachedfile'
      },      
      publication: {
        description: "Publications",
        title: "Publications",
        invoke: "setPublications",
        resource: 'publication'
      },      
      billingItems: {
        description : 'Billing information',
        title : 'Billing',
        invoke : 'setBilling',
        permission: 'screenbilling'
      },
      datacolumns : {
        description : 'Data Columns',
        title : 'Data Columns',
        invoke : 'setDatacolumns',
        permission: 'screenresult'
      },
      results : {
        description : 'Screen Results',
        title : 'Screen Results',
        invoke : 'setResults',
        permission: 'screenresult'
      },
      cherrypicks: {
        description : 'Cherry Pick Requests',
        title : 'Cherry Picks',
        invoke : 'setCherryPicks',
        permission: 'cherrypick'
      },
      activities: {
        description : 'Activities',
        title : 'Activities',
        invoke : 'setActivities',
        permission: 'activity'
      },
    },

    events: {
      'click ul.nav-tabs >li': 'click_tab',
      'click button#loadScreenResults': 'loadScreenResults',
      'click button#deleteScreenResults': 'deleteScreenResults',
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      var displayed_tabbed_resources = _.extend({},this.tabbed_resources);
      if (! self.model.get('has_screen_result')){
        delete displayed_tabbed_resources['results'];
        delete displayed_tabbed_resources['datacolumns'];
      }
      return {
        'base_url': self.model.resource.key + '/' + self.model.key,
        'tab_resources': displayed_tabbed_resources
      }      
    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      var self = this;
      console.log('afterRender called...');
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if (viewId == '+add') {
          this.$('ul.nav-tabs > li').addClass('disabled');
          this.uriStack.unshift(viewId);
          viewId = 'detail';
          delete self.tabbed_resources['summary'];
        }else if (viewId == 'edit'){
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }else if (_.contains(['libraries'],viewId)){
          this.consumedStack = [viewId];
          this.showLibraries(this.uriStack);
          return;
        }else if (_.contains(['copyplates'],viewId)){
          this.consumedStack = [viewId];
          this.showCopyPlates(this.uriStack);
          return;
        }else if (_.contains(['copyplatesloaded'],viewId)){
          this.consumedStack = [viewId];
          this.showCopyPlatesLoaded(this.uriStack);
          return;
        }else if (viewId == 'libraryscreening'){
          if(_.contains(this.uriStack, '+add') ){
            this.addLibraryScreening();
            return;
          }else{
            this.consumedStack = [viewId];
            this.showLibraryScreenings(this.uriStack); 
            return;
          }
        }
        
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
      console.log('afterRender, done.');
    },

    click_tab : function(event){
      var self = this;
      event.preventDefault();
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      if(this.$('#'+key).hasClass('disabled')){
        return;
      }
      if(this.key && this.key === key){
        return;
      }
      appModel.requestPageChange({
        ok: function(){
          self.change_to_tab(key);
        }
      });
    },

    change_to_tab: function(key){
      console.log('change_to_tab', key);
      
      if(_.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active');
        this.$('#' + key).addClass('active');
        if(key !== 'detail'){
          this.consumedStack = [key];
        }else{
          this.consumedStack = [];
        }
        var delegateStack = _.clone(this.uriStack);
        this.uriStack = [];
        var method = this[this.tabbed_resources[key]['invoke']];
        if (_.isFunction(method)) {
          method.apply(this,[delegateStack]);
        } else {
          throw "Tabbed resources refers to a non-function: " + this.tabbed_resources[key]['invoke']
        }
      }else{
        var msg = 'Unknown tab: ' + key;
        appModel.error(msg);
        throw msg;
      }
    },

    render : function(){
      window.scrollTo(0, 0);
      return this;
    },

    setDetail: function(delegateStack){
      var self,outerSelf = self = this;
      var key = 'detail';
      var fields = self.model.resource.fields;
      // set up a custom vocabulary that joins username to name; will be 
      // used as the text of the linklist
      if (this.model.get('project_phase')!= 'annotation'){
        fields['collaborator_usernames'].vocabulary = (
            _.object(this.model.get('collaborator_usernames'),
              this.model.get('collaborator_names')));
      }
      var editView = EditView.extend({
        afterRender: function(){
          outerSelf._addVocabularyButton(
            this, 'cell_lines', 'cell_line', 'Cell Line', { description: 'ATCC Designation' });
         
          outerSelf._addVocabularyButton(
            this, 'transfection_agent', 'transfection_agent', 'Transfection Agent');
         
          outerSelf._addVocabularyButton(
            this, 'species', 'screen.species', 'Screened Species');
          EditView.prototype.afterRender.apply(this,arguments);
        }
      });
      
      var view = this.tabViews[key];
      if (view) {
        this.removeView(this.tabViews[key]);
      }
      
      var detailView = DetailView.extend({
        afterRender: function(){
          DetailView.prototype.afterRender.apply(this,arguments);
          
          if(self.model.get('project_phase')=='annotation'){
            if (appModel.hasPermission('screenresult','write')){
              if (self.model.get('has_screen_result')){
                this.$el.prepend(self.$deleteScreenResultsButton);
              }
              this.$el.prepend(self.$loadScreenResultsButton);
            }
          }else{
            self.createStatusHistoryTable($('#screen_extra_information'));
            self.createActivitySummary($('#screen_extra_information'));
            self.createCprTable($('#screen_extra_information'));
            
            // TODO: create a metadata setting for "show if present"
            if (!self.model.has('perturbagen_molar_concentration')){
              $('#perturbagen_molar_concentration').closest('tr').remove();
            }
            if (!self.model.has('perturbagen_ug_ml_concentration')){
              $('#perturbagen_ug_ml_concentration').closest('tr').remove();
            }
            if (!self.model.has('transfection_agent')){
              $('#transfection_agent').closest('tr').remove();
            }
          }
          
        },
        
        serialize: function() {
          console.log('serialize...');
          var data = DetailView.prototype.serialize.apply(this,arguments);
          var informationKeys = [];
          var groupedKeys = [];
          data['groupedKeys'].each(function(groupKey){
            if(_.result(groupKey,'title',null) == 'Information'){
              informationKeys = informationKeys.concat(groupKey.fields);
            }else{
              groupedKeys.push(groupKey);
            }
          });
          console.log('groupedKeys', groupedKeys, informationKeys);
          data['groupedKeys'] = _.chain(groupedKeys);
          data['informationKeys'] = _.chain(informationKeys);
          return data;
        },
        
        template: _.template(screenTemplate)
        
      });
      
      var ScreenModel = Backbone.Model.extend({
        key: self.model.key,
        resource: self.model.resource,
        url: self.model.url,
        parse: self.model.parse,
        validate: function(attrs, options) {
          errs = {};
          if (!_.isEmpty(_.result(attrs,'data_privacy_expiration_notified_date'))){
            if (!_.isEmpty(_.result(attrs,'max_allowed_data_privacy_expiration_date'))){
              errs['max_allowed_data_privacy_expiration_date'] = (
                'can not be set if the expiration notified date is set');
            }
            if (!_.isEmpty(_.result(attrs,'min_allowed_data_privacy_expiration_date'))){
              errs['min_allowed_data_privacy_expiration_date'] = (
                'can not be set if the expiration notified date is set');
            }
          }
          return _.isEmpty(errs) ? null : errs;
        }
      });
      this.model = new ScreenModel(this.model.attributes);
      
      view = new DetailLayout({ 
        model: this.model, 
        uriStack: delegateStack,
        EditView: editView,
        DetailView: detailView
      });
      view.showEdit = function(){
        appModel.initializeAdminMode(function(){
          var userOptions = appModel.getUserOptions();
          fields['collaborator_usernames']['choices'] = userOptions;
          fields['lead_screener_username']['choices'] = (
              [{ val: '', label: ''}].concat(userOptions));
          fields['lab_head_username']['choices'] = (
              appModel.getPrincipalInvestigatorOptions() );
          fields['publishable_protocol_entered_by']['choices'] = (
              appModel.getAdminUserOptions() );
          // TODO: this should default only after pp value is entered
          //view.model.set('publishable_protocol_entered_by',
          //    appModel.getCurrentUser().username);
          DetailLayout.prototype.showEdit.apply(view,arguments);
        });  
      };

      this.tabViews[key] = view;
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.consumedStack = []; 
      this.setView("#tab_container", view ).render();
      return view;
    },

    /**
     * Libraries view is a sub-view of Summary
     */
    showLibraries: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'libraries'].join('/');
      var resource = appModel.getResource('library');
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: []
      }});
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event){
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Libraries for screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');

    },
    
    /**
     * Library Copy Plates view is a sub-view of Summary
     */
    showCopyPlates: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'copyplates'].join('/');
      var resource = appModel.getResource('librarycopyplate');
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: []
      }});
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event){
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Library Copy Plates for Screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');

    },
    
    /**
     * Library Copy Plates Loaded view is a sub-view of Summary
     */
    showCopyPlatesLoaded: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'copyplatesloaded'].join('/');
      var resource = appModel.getResource('librarycopyplate');
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: []
      }});
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event){
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Library Copy Plates loaded for Screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');

    },    
    
    setAttachedFiles: function(delegateStack) {
      var self = this;
      var key = 'attachedfile';
      var resource = appModel.getResource('attachedfile');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'attachedfiles'].join('/');
      extraControls = []
      if (appModel.hasPermission('screen','write')){
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
        uploadAttachedFileButton.click(function(e){
          e.preventDefault();
          var url = [self.model.resource.apiUri, 
             self.model.key,
             'attachedfiles'].join('/');
          EditView.uploadAttachedFileDialog(url, view.collection, 'attachedfiletype.screen');
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
        extraControls = extraControls.concat(uploadAttachedFileButton, showDeleteButton);
      }      
      
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: extraControls
      }});
      Backbone.Layout.setupView(view);
      self.consumedStack = [key]; 
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      
    },
        
    setPublications: function(delegateStack) {
      var self = this;
      var key = 'publication';
      var resource = appModel.getResource(key);
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'publications'].join('/');
      extraControls = []
      if (appModel.hasPermission('screen','write')){
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
        showDeleteButton.click(function(e){
          e.preventDefault();
          if (! view.grid.columns.findWhere({name: 'deletor'})){
            view.grid.columns.unshift({ 
              name: 'deletor', label: 'Delete', text:'X', 
              description: 'delete record', 
              cell: Iccbl.DeleteCell, sortable: false });
          }
        });
        showAddButton.click(function(e){
          e.preventDefault();
          self.addPublicationDialog(view.collection);
        });
        extraControls = extraControls.concat(showAddButton, showDeleteButton)
      }      
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: extraControls
      }});
      Backbone.Layout.setupView(view);
      self.consumedStack = [key]; 
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      
    },
    
    addPublicationDialog: function(collection){
      var self = this;
      var url = [self.model.resource.apiUri, 
         self.model.key,
         'publications'].join('/');
      var form_template = [
         "<div class='form-horizontal container' id='uploadAttachedFileButton_form' >",
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
        validators: ['required'],
        template: fieldTemplate
      };
      formSchema['authors'] = {
        title: 'Authors',
        key: 'authors',
        type: 'TextArea',
        validators: ['required'],
        template: fieldTemplate
      };
      formSchema['journal'] = {
        title: 'Journal',
        key: 'journal',
        type: 'TextArea',
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

      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          var errs = {};
          if (_.isEmpty(attrs.pubmed_central_id) 
              && _.isEmpty(attrs.pubmed_id)){
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
              var comments = _.result(values,'comments','');
              var publication_id = _.result(values,'pubmed_id', _.result(values,'pubmed_central_id'));
              var headers = {};
              headers[appModel.HEADER_APILOG_COMMENT] = comments;
              
              var data = new FormData();
              _.each(_.keys(values), function(key){
                if(values[key]){
                  data.append(key,values[key]);
                }
              });
              
              var file = $('input[name="fileInput"]')[0].files[0];
              var filename;
              if(file){
                data.append('attached_file',file);
                filename = file.name;
                if (!_.isEmpty(values['filename'])){
                  filename = values['filename'];
                }
                if (_.isEmpty(values['file_date'])){
                  data.append(
                    'file_date',
                    Iccbl.getISODateString(_.result(file,'lastModifiedDate', null)));
                }
              }else{
                filename = values['filename'];
              }
              
              $.ajax({
                url: url,    
                data: data,
                cache: false,
                contentType: false,
                processData: false,
                type: 'POST',
                headers: headers, 
                success: function(data){
                  collection.fetch({ reset: true });
                  appModel.showModalMessage({
                    title: 'Publication added',
                    okText: 'ok',
                    body: '"' + publication_id + '"'
                  });
                },
                done: function(model, resp){
                  // TODO: done replaces success as of jq 1.8
                  console.log('done');
                }
              }).fail(function(){ appModel.jqXHRfail.apply(this,arguments); });
            
              return true;
            }
          },
          view: _form_el,
          title: 'Upload an Attached File'  });
      
    },
    
    createActivitySummary: function($target_el){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'activities'].join('/');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: url
      });
      var cell = $('<div id="activity_summary" class="row"/>');
      $target_el.append(cell);

      function build_table(collection){
          if (collection.isEmpty()){
            return;
          }
          var activity = collection.at(0);
          
          cell.append($([
            '<div class="col-xs-12"><strong>Activity Summary</strong></div>',
            '<div id="" class="col-xs-12" >',
            '<table id="activity_summary" class="table-condensed data-list">',
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
      if(self.model.has('latest_activities_data')){
        build_table(new CollectionClass(self.model.get('latest_activities_data')));
      }else{
        var cpr_collection = new CollectionClass();
        cpr_collection.fetch({
          data: { 
            limit: 1,
            order_by: ['-date_of_activity']
          },
          success: build_table
        }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      }
    },
    
    createCprTable: function($target_el){
      var self = this;
      var cprResource = appModel.getResource('cherrypickrequest');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: cprResource.apiUri 
      });
      var cell = $('<div id="cpr_table" class="row"/>');
      $target_el.append(cell);
      
      function build_table(collection){
        if (collection.isEmpty()){
          return;
        }
        collection.each(function(model){
        });
        var TextWrapCell = Backgrid.Cell.extend({
          className: 'text-wrap-cell'
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
              hrefTemplate: '#cherrypickrequest/{cherry_pick_request_id}'
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
            'cell': TextWrapCell
          })
        ];
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();

        cell.append($([
          '<div class="col-xs-12"><strong>',
          'Recent <a href="#screen/' + self.model.get('facility_id'),
          '/cherrypicks">Cherry Pick Requests</a></strong></div>',
          '<div class="col-xs-12" id="cprs"/>'].join('')));
        
        var _grid = new Backgrid.Grid({
          columns: colModel,
          collection: collection,
          className: 'backgrid table-striped table-condensed table-hover '
        });
        cell.find('#cprs').html(_grid.render().$el);
        
      }
      
      if(self.model.has('cherry_pick_request_data')){
        build_table(new CollectionClass(self.model.get('cherry_pick_request_data')));
      }else{
        var cpr_collection = new CollectionClass();
        cpr_collection.fetch({
          data: { 
            limit: 10,
            screen_facility_id__eq: self.model.key,
            order_by: ['-date_requested']
          },
          success: build_table
        }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      }
    },
    
    /**
     * Update the screen status with a status history table: populate
     * using the apilog history of the status attribute
     **/
    createStatusHistoryTable: function($target_el){
      var self = this;
      var apilogResource = appModel.getResource('apilog');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: apilogResource.apiUri 
      });
      var cell = $('<div id="status_table" class="row"/>');
      $target_el.append(cell);
      
      function build_table(collection){
        if (collection.isEmpty()){
          return;
        }
        collection.each(function(model){
          var diffs = JSON.parse(model.get('diffs'));
          model.set('status', appModel.getVocabularyTitle('screen.status', diffs.status[1]));
        });
        var TextWrapCell = Backgrid.Cell.extend({
          className: 'text-wrap-cell'
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
              'cell': TextWrapCell
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
        cell.append($([
          '<div class="col-xs-12"><strong>Status Items</strong></div>',
          '<div class="col-xs-12" id="status_items"/>'].join('')));
        
        var status_grid = new Backgrid.Grid({
          columns: colModel,
          collection: collection,
          className: 'backgrid table-striped table-condensed table-hover '
        });
        cell.find('#status_items').html(status_grid.render().$el);
        
      }
      
      if(self.model.has('status_data')){
        build_table(new CollectionClass(self.model.get('status_data')));
      }else{
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
        }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      }
      
    },
    
    _addVocabularyButton: function(
        editForm, fieldKey, vocabulary_scope_ref, vocabulary_name, options){
      var options = options || {};
      var addButton = $([
        '<a class="btn btn-default btn-sm" ',
          'role="button" id="add_' + fieldKey + '_button" href="#">',
          'Add</a>'
        ].join(''));
      addButton.click(function(event){
        event.preventDefault();
        appModel.addVocabularyItemDialog(
          vocabulary_scope_ref, vocabulary_name,
          function(new_vocab_item){
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
    
    setCherryPicks: function(delegateStack){
      var self = this;
      var key = 'cherryPicks';
      var view = this.tabViews[key];
      
      if (!view){
        var self = this;
        var url = [self.model.resource.apiUri,self.model.key,'cherrypicks'].join('/');
        var resource = appModel.getResource('cherrypickrequest');
        var view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: resource,
          resource: resource,
          url: url,
          extraControls: []
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        this.$('li').removeClass('active');
        this.$('#'+key).addClass('active');

      }else{
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
    },
      
    setActivities: function(delegateStack){
      var self = this;
      var key = 'activities';
      var view = this.tabViews[key];
      
      if (!view){
        var self = this;
        var url = [self.model.resource.apiUri,self.model.key,'activities'].join('/');
        var resource = appModel.getResource('activity');
        
        var sa_vocab = appModel.getVocabulary('activity.type');
        resource.fields['type'].vocabulary = 
          _.map(sa_vocab, function(v){
            return [v.title,v.key];
          }); // TODO: app model method for this
        
        console.log('combined vocab', resource.fields['type'].vocabulary );
        
        var view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: resource,
          resource: resource,
          url: url,
          extraControls: []
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        this.$('li').removeClass('active');
        this.$('#'+key).addClass('active');

      }else{
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
    },
    
    /**
     * Loads the screen results using ajax to post the file.
     * - because ajax cannot handle post response attachments (easily), set
     * the response type to 'application/json' and display the errors in a 
     * modal dialog.
     * 
     * NOTE: a "POST" form cannot be used - there is no standard way to signal
     * the result of the post to JavaScript.
     */
    loadScreenResults: function(){
      var self = this;
      var isStudy = self.model.get('project_phase') == 'annotation';
      var form_template = [
         "<div class='form-horizontal container' id='screenresult_form' >",
         "<form data-fieldsets class='form-horizontal container' >",
         "<div class='form-group' ><input type='file' name='fileInput' /></div>",
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
          var errs = {};
          var file = $('input[name="fileInput"]')[0].files[0]; 
          if (! file) {
            errs.fileInput = 'Must specify a file';
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

      function saveFile() {

        var errors = form.commit({ validate: true }); // runs schema and model validation
        if(!_.isEmpty(errors) ){
          console.log('form errors, abort submit: ',errors);
          _.each(_.keys(errors), function(key){
            $('[name="'+key +'"').parents('.form-group').addClass('has-error');
            if( key == '_others'){
              _.each(errors[key], function(otherError){
                // TODO: display file req'd msg
              });
            }
          });
          return false;
        }else{
          var url = [self.model.resource.apiUri,self.model.key,'screenresult'].join('/');
          var values = form.getValue();
          var file = $('input[name="fileInput"]')[0].files[0]; 
          var data = new FormData();
          // NOTE: 'xls' is sent as the key to the FILE object in the upload.
          // Use this as a non-standard way to signal the upload type to the 
          // serializer. 
          data.append('xls', file);
          var headers = {
            'Accept': 'application/json'
          };
          headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
          $.ajax({
            url: url,    
            data: data,
            cache: false,
            contentType: false,
            processData: false,
            type: 'POST',
            headers: headers
          })
          .fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); })
          .done(function(model, resp){
            // FIXME: should be showing a regular message
            appModel.error('success');
            // Must refetch and render, once for the containing tabbed layout,
            // and then (during render), will refetch again for the summary view
            self.model.fetch({
              success: function() {
                if (!isStudy){
                  self.uriStack = ['summary'];
                  // remove the child view before calling render, to prevent
                  // it from being rendered twice, and calling afterRender twice
                  self.removeView('#tab_container');
                }
                self.render();
              }
            }).fail(function(){ 
              Iccbl.appModel.jqXHRfail.apply(this,arguments); 
            });
          });
        }        
      }
      
      var dialog = appModel.showModal({
          ok: saveFile,
          view: _form_el,
          title: 'Select a Screen Result (xlsx workbook) file to upload'
      });
    },
    
    deleteScreenResults: function(){
      var self = this;
      var isStudy = self.model.get('project_phase') == 'annotation';
      var formSchema = {};
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        type: 'TextArea'
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields
      });
      var _form_el = form.render().el;
      var dialog = appModel.showModal({
          ok: function(){
            var errors = form.commit({ validate: true }); 
            if(!_.isEmpty(errors) ){
              _.each(_.keys(errors), function(key){
                $('[name="'+key +'"').parents('.form-group').addClass('has-error');
              });
              return false;
            }else{
              var url = [self.model.resource.apiUri,self.model.key,'screenresult'].join('/');
              var values = form.getValue();
              var headers = {
                'Accept': 'application/json'
              };
              headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
              $.ajax({
                url: url,    
                cache: false,
                contentType: false,
                processData: false,
                type: 'DELETE',
                headers: headers
              })
              .fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); })
              .done(function(model, resp){
                appModel.error('success');
                // Must refetch and render, once for the containing tabbed layout,
                // and then (during render), will refetch again for the summary view
                self.model.fetch({
                  success: function() {
                    if (!isStudy){
                      self.uriStack = ['summary'];
                      // remove the child view before calling render, to prevent
                      // it from being rendered twice
                      self.removeView('#tab_container');
                    }
                    self.render();
                  }
                }).fail(function(){ 
                  Iccbl.appModel.jqXHRfail.apply(this,arguments); 
                });
              });
            }        
          },
          view: _form_el,
          title: 'Delete screen results for ' + self.model.get('facility_id')
      });
      
    },
    
    addLibraryScreening: function(){
      console.log('add library screening');
      var self = this;
      var defaults = {
        screen_facility_id: self.model.get('facility_id'),
        screen_type: self.model.get('screen_type')
      };
      var newModel = appModel.createNewModel('libraryscreening', defaults);

      var view = new LibraryScreening({ 
        model: newModel, 
        uriStack: ['+add']
      });
      
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();

      
      this.consumedStack = ['libraryscreening','+add'];
      self.reportUriStack([]);
    },
    
    showLibraryScreenings: function(delegateStack){
      var self = this;
      var lsResource = appModel.getResource('libraryscreening'); 

      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        var activityId = delegateStack.shift();
        self.consumedStack.push(activityId);
        var _key = [self.model.key,'libraryscreening',activityId ].join('/');
        appModel.getModel(lsResource.key, activityId, function(model){

          var view = new LibraryScreening({ 
            model: model, 
            uriStack: _.clone(delegateStack),
          });
          
          Backbone.Layout.setupView(view);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
        });        
        return;
      }else{
        // List view
        // special url because list is a child view for library
        var url = [self.model.resource.apiUri, 
                   self.model.key,
                   'libraryscreening'].join('/');
        view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: lsResource,
          resource: lsResource,
          url: url
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        this.setView("#tab_container", view ).render();
      }
    },
    
    setSummary : function(delegateStack){
      console.log('setSummary...', delegateStack);
      var self = this;
      var key = 'summary';
      var $addLibraryScreeningButton = $(
        '<a class="btn btn-default btn-sm" role="button" \
        id="addLibraryScreening" href="#">Add Library Screening</a>');
      $addLibraryScreeningButton.click(function(e){
        e.preventDefault();
        self.addLibraryScreening();
      });
      
      function viewLoadHistory(e){
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', 'search'];
        var search = {};
        search['ref_resource_name'] = 'screenresult';
        search['key'] = self.model.key;
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        appModel.router.navigate(route, {trigger: true});
        self.remove();
      }
      
      var summaryKeys = self.model.resource.filterKeys('visibility', 'summary');
      var summaryModel = appModel.getModel(
        self.model.resource.key, self.model.key, 
        function(model){
          view = new DetailView({ 
            events: {
              'click button#history': viewLoadHistory
            },
            model: model, 
            uriStack: _.clone(delegateStack),
            detailKeys: summaryKeys,
            buttons: ['history'],
            afterRender: function(){
              DetailView.prototype.afterRender.apply(this);
              this.$el.find('#libraries_screened_count').click(function(e){
                e.preventDefault();
                self.consumedStack = ['libraries'];
                self.showLibraries(delegateStack);
              });
              this.$el.find('#library_plates_screened').click(function(e){
                e.preventDefault();
                self.consumedStack = ['copyplates'];
                self.showCopyPlates(delegateStack);
              });
              this.$el.find('#library_plates_data_loaded').click(function(e){
                e.preventDefault();
                self.consumedStack = ['copyplates'];
                self.showCopyPlatesLoaded(delegateStack);
              });

              if (self.model.get('has_screen_result')){
                self._createPositivesSummary(this.$el.find('#positives_summary'));
              }
              if (appModel.hasPermission('screenresult','write')){
                if (self.model.get('has_screen_result')){
                  this.$el.prepend(self.$deleteScreenResultsButton);
                }
                this.$el.prepend(self.$loadScreenResultsButton);
              }
              if (appModel.hasPermission('libraryscreening','write')){
                this.$el.prepend($addLibraryScreeningButton);
              }
            }
          });
          self.tabViews[key] = view;
          
          self.setView("#tab_container", view ).render();
          
          self.reportUriStack([]);
          
        },{ data_for_get: { visibilities: ['summary']} }
      );
    },
    
    _createPositivesSummary: function($target_el){
      var self = this;

      var experimental_wells_loaded = self.model.get('experimental_well_count');
      if(!_.isNumber(experimental_wells_loaded) ||
          experimental_wells_loaded < 1 ){
        return;
      }
      function createPositiveStat(raw_value){
        var formatter = new Iccbl.DecimalFormatter({ decimals: 2 });
        if (!raw_value) raw_value = 0;
        if( raw_value === 0 || !_.isNumber(experimental_wells_loaded) ||
            experimental_wells_loaded < 1 ){
          return ''
        }
        return Iccbl.formatString(
            '{val} ({percent}%)', {
              val: raw_value,
              percent: formatter.fromRaw(
                100.0*raw_value/experimental_wells_loaded )
            });
      }
      
      var dcResource = appModel.getResource('datacolumn');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: dcResource.apiUri 
      });
      
      var dcCollection = new CollectionClass();
      dcCollection.fetch({
        data: { 
          limit: 0,
          screen_facility_id__eq: self.model.get('facility_id'),
          data_type__in: [
            'partition_positive_indicator','boolean_positive_indicator',
            'confirmed_positive_indicator'],
          includes: ['positives_count', 'strong_positives_count',
                     'medium_positives_count','weak_positives_count'],
          order_by: ['ordinal']
        },
        success: function(collection, response) {
          if (!collection || collection.isEmpty()){
            return;
          }
          collection.each(function(dc){
            dc.set('total_positives', createPositiveStat(dc.get('positives_count')));
            dc.set('strong_positives', createPositiveStat(dc.get('strong_positives_count')));
            dc.set('medium_positives', createPositiveStat(dc.get('medium_positives_count')));
            dc.set('weak_positives', createPositiveStat(dc.get('weak_positives_count')));
          });
          
          var TextWrapCell = Backgrid.Cell.extend({
            className: 'text-wrap-cell'
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
                'name' : 'name',
                'label' : 'Data Column',
                'description' : 'Data Column',
                'order': 1,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'total_positives',
                'label' : 'Total Positives',
                'description' : 'Total Positives',
                'order': 2,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'strong_positives',
                'label' : 'Strong Positives',
                'description' : 'Strong Positives',
                'order': 3,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'medium_positives',
                'label' : 'Medium Positives',
                'description' : 'Medium Positives',
                'order': 4,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'weak_positives',
                'label' : 'Weak Positives',
                'description' : 'Weak Positives',
                'order': 5,
                'sortable': true,
                'cell': TextWrapCell
              }),
          ];
          var colModel = new Backgrid.Columns(columns);
          colModel.comparator = 'order';
          colModel.sort();
          var positives_grid = new Backgrid.Grid({
            columns: colModel,
            collection: collection,
            className: 'backgrid table-striped table-condensed table-hover'
          });
          
          $target_el.empty();
          var cell = $('<div>',{ class: 'col-xs-4' });
          cell.html(positives_grid.render().$el);
          $target_el.append(cell);
        }
      }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });
      
    },

    setBilling: function(delegateStack){
      var self = this;
      var key = 'billing';
      var view = this.tabViews[key];
      if (!view){
        var billingKeys = self.model.resource.filterKeys('visibility', 'billing');
        var summaryModel = appModel.getModel(
          self.model.resource.key, self.model.key, 
          function(model){
            view = new DetailLayout({ 
              model: model, 
              uriStack: delegateStack,
              detailKeys: billingKeys,
              editKeys: billingKeys,
              editableKeys: billingKeys
            });
            self.tabViews[key] = view;
            
            self.listenTo(view , 'uriStack:change', this.reportUriStack);
            self.setView("#tab_container", view ).render();
            self.reportUriStack([]);
          },{ data_for_get: { visibilities: ['billing']}});
      }else{
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
    },

    setDatacolumns: function(delegateStack){
      
      var self = this;
      var datacolumnResource = appModel.getResource('datacolumn'); 
      var url = [self.model.resource.apiUri,self.model.key,'datacolumns'].join('/');
      
      // construct a collection-in-columns
      // pivot the datacolumn list
      Iccbl.getCollectionOnClient(url, function(collection){
        // create a colModel for the list
        var columns = [];
        var TextWrapCell = Backgrid.Cell.extend({
          className: 'text-wrap-cell'
        })
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell.extend({
            render: function(){
              this.$el.empty();
              var label = Iccbl.createLabel(
                this.column.get("label"), 10);
              this.$el.append(label);
              return this;
            }
          })
        };
        // setup the first column - schema field titles (that are visible)
        var col0 = _.extend({},colTemplate,{
          'name' : 'fieldname',
          'label' : '',
          'description' : 'Datacolumn field',
          'cell': TextWrapCell.extend({
            className: '150_width'
          })
        });
        columns.push(col0);
        
        collection.each(function(datacolumn){
          var col = _.extend({},colTemplate,{
            'name' : datacolumn.get('key'),
            'label' : datacolumn.get('title'),
            'description' : datacolumn.get('description'),
            'order': datacolumn.get('ordinal'),
            'sortable': true,
            'cell': TextWrapCell
          });
          columns.push(col);
          if (datacolumn.has('derived_from_columns')){
            datacolumn.set(
              'derived_from_columns', 
              datacolumn.get('derived_from_columns').join(', '));
          }
        });
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();
        
        // create the collection by pivoting
        orderedFields = _.sortBy(datacolumnResource.fields,'ordinal');
        var pivotCollection = new Backbone.Collection();
        _.each(orderedFields, function(field){
          if(_.contains(field.visibility, 'l') 
              && !_.contains(['name','title','key'],field.key )){
            var row = {'key': field.key, 'fieldname': field.title };
            collection.each(function(datacolumn){
              row[datacolumn.get('key')] = datacolumn.get(field.key);
            });
            pivotCollection.push(row);
          }
        });
        
        // now show as a list view
        var view = new Backgrid.Grid({
          columns: colModel,
          collection: pivotCollection,
          className: 'backgrid table-striped table-condensed table-hover'
        });
        var view = view.render();
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);

        self.setView("#tab_container", view ).render();
     });
      
    },
    
    setResults : function(delegateStack){
      var self = this;
      var screenResultResource = appModel.getResource('screenresult');
      var schemaUrl = [appModel.dbApiUri,
                       'screenresult',
                       self.model.key,
                       'schema'].join('/');
      var url = [appModel.dbApiUri,
                 'screenresult',
                 self.model.key].join('/');

      var _id = self.model.key;
      var screen_facility_id = self.model.get('facility_id');
      
      var createResults = function(schemaResult){
        var extraControls = [];
        var show_positives_control = $([
          '<label class="checkbox-inline">',
          '  <input type="checkbox">positives',
          '</label>'
          ].join(''));
        var show_mutual_positives_control = $([
          '<label class="checkbox-inline">',
           '  <input type="checkbox">mutual positives',
           '</label>'
           ].join(''));
        if (self.model.get('project_phase') != 'annotation'){
          extraControls = extraControls.concat(
            show_positives_control, show_mutual_positives_control);
        }
        // create an option vocab for the exclued col, if needed
        if (_.has(schemaResult['fields'], 'excluded')){
          var options = [];
          _.each(schemaResult['fields'],function(field){
            var dc_col_key = 'dc_' + self.model.key + '_';
            if (field['key'].indexOf(dc_col_key)>-1){
              var option_name = field['key'].substring(dc_col_key.length);
              options.unshift([option_name, option_name]);
            }
          });
          schemaResult['fields']['excluded']['vocabulary'] = options;
          console.log('custom exclude options', options);
        }
        var initialSearchHash;
        view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: schemaResult,
          resource: screenResultResource,
          url: url,
          extraControls: extraControls
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        initialSearchHash = view.listModel.get('search');
        
        if(_.has(initialSearchHash, 'is_positive__eq')
            && initialSearchHash.is_positive__eq.toLowerCase()=='true'){
          show_positives_control.find('input[type="checkbox"]').prop('checked',true);
        }
        if(_.has(initialSearchHash, 'show_mutual_positives')
            && initialSearchHash.show_mutual_positives.toLowerCase()=='true'){
          show_mutual_positives_control.find('input[type="checkbox"]').prop('checked',true);
          view.show_mutual_positives(screen_facility_id, true);
        }
        if(_.has(initialSearchHash, 'other_screens')){
          view.show_other_screens(initialSearchHash['other_screens']);
        }
        show_positives_control.click(function(e){
          if(e.target.checked){
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['is_positive__eq'] = 'true';
            view.listModel.set('search',searchHash);
          }else{
            // make sure unset
            var searchHash = _.clone(view.listModel.get('search'));
            if(_.has(searchHash,'is_positive__eq')){
              delete searchHash['is_positive__eq'];
              view.listModel.set('search',searchHash);
            }
          }
        });
        show_mutual_positives_control.click(function(e){
          if(e.target.checked){
            window.setTimeout(function(){
              view.show_mutual_positives(screen_facility_id, true);
            });
          }else{
            window.setTimeout(function(){
              view.show_mutual_positives(screen_facility_id, false);
            });
          }
        });
        
      };
      appModel.getResourceFromUrl(schemaUrl, createResults);
    },
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });

  return ScreenView;
});