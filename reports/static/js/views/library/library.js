/**
 * Library form/view
 */
define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/library/libraryCopy', 
  'views/library/libraryWell', 
  'views/library/libraryWells', 
  'views/generic_detail_layout',
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'views/list2',
  'utils/uploadDataForm',
  'utils/tabbedController'
  
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, LibraryCopyView, 
            LibraryWellView, LibraryWellsView, DetailLayout, DetailView, EditView, ListView, 
            UploadDataForm, TabbedController ) {

  var LibraryView = TabbedController.extend({
    
    initialize: function(args) {
      var self = this;
      this._classname = 'LibraryView';
      var tabbed_resources = this.tabbed_resources;
      TabbedController.prototype.initialize.apply(this,arguments);
      
      // add back in resources viewable for screeners
      if (appModel.getCurrentUser().is_staff !== true){
        this.tabbed_resources['well'] = tabbed_resources['well'];
      }

      _.bindAll(this, 'upload','download_contents');
      console.log('LibraryView initialized');
    },
    
    tabbed_resources: {
      detail: { 
        description: 'Library Details', 
        title: 'Library Details', 
        invoke: 'setDetail'
      },
      copy: { 
        description: 'Copies', title: 'Copies', invoke: 'setCopies',
        resource: 'librarycopy'
      },
      plate: { 
        description: 'Plates', title: 'Plates', invoke: 'setPlates' ,
        resource: 'librarycopyplate'
      },
      well: { 
        description: 'Well based view of the library contents', 
        title: 'Wells', invoke: 'setWells' ,
        resource: 'well'
      }
    },      
    
    events: {
      'click ul.nav-tabs >li': 'click_tab',
      'click button#upload': 'upload',
    },
    
    afterRender: function(){
      
      console.log('LibraryView.afterRender called...')
      
      var self = this;
      TabbedController.prototype.afterRender.apply(this,arguments);

      function showReleaseMessage(){
        $('#content_title_message').find('#release_message').remove();
        
        if (self.model.isNew() || self.model.get('is_released') == true){
          return;
        }

        var helpMsg = [
              'Library has not been set to "released":',
              '- data will not be available to screeners,',
              '- Contents will be loaded without a preview,',
              '- All other functions and logging are unaffected.'];
        var releaseMessage = $(
          '<div id="release_message" class="alert alert-info"></div>')
            .html('Library has not been');
        releaseMessage.attr('title', helpMsg.join('\n'));
        
        var releaseLibraryInfo = $('<a>', {
          tabIndex : -1,
          href : "#",
          class: 'alert-link',
          title: "Info"
        }).text('released');
        releaseMessage.append('&nbsp;');
        releaseMessage.append(releaseLibraryInfo);
        
        var releaseLibraryLink = $('<a>', {
          tabIndex : -1,
          href : "#",
          class: 'alert-link',
          title: "Set the library to released"
        }).text('[release now]');
        releaseMessage.append('&nbsp;');
        releaseMessage.append(releaseLibraryLink);
        
        $('#content_title_message').append(releaseMessage);
        
        releaseLibraryInfo.click(function(e){
          e.preventDefault();
          e.stopPropagation();
          
          appModel.showModalMessage({
            'title': 'Library is not released',
            'body': helpMsg.join('<br/>')
          });
        });
        
        releaseLibraryLink.click(function(e){
          e.preventDefault();
          e.stopPropagation();

          appModel.showOkCommentForm({
            ok: function(formValues){
              var options = {
                patch: true
              };
              var comment = formValues['comment'];
              if (comment){
                options['headers'] = {};
                options['headers'][appModel.HEADER_APILOG_COMMENT] = comment;
              }
              self.model.set('is_released',true);
              self.model.save({'is_released': true}, options)
                .done(function(data, textStatus, jqXHR){ 
                  console.log('done releasing...');
                  self.setDetail([]);
                }).fail(function() { 
                  Iccbl.appModel.jqXHRfail.apply(this,arguments); 
                });
            },
            okText: 'Release',
            title: 'Release the library: "'
              + Iccbl.getTitleFromTitleAttribute(self.model,self.model.resource) + '" ?' 
          });

        });
        
      };
      
      function showPreviewMessage(){
        $('#content_title_message').find('#preview_message').remove();
        var preview_log_id = self.model.get('preview_log_id');
        
        if (!_.isNumber(preview_log_id)){
          return;
        }
        var previewMessage = $(
          '<div id="preview_message" class="alert alert-info"></div>')
            .html('Library has a preview version of well data loaded. ');
        var showPreviewLink = $('<a>', {
          tabIndex : -1,
          href : "#",
          class: 'alert-link',
          title: "Show the well view of the preview data"
        }).text('[view]');
        previewMessage.append(showPreviewLink);
        
        var previewUriStack= ['#apilog', self.model.get('preview_log_key')]
        var showPreviewDetails = $('<a>', {
          tabIndex : -1,
          href : previewUriStack.join('/'),
          class: 'alert-link',
          target: '_blank',
          title: "Display the API log for the preview load operation"
        }).text('[information]');
        previewMessage.append('&nbsp;');
        previewMessage.append(showPreviewDetails);
        $('#content_title_message').append(previewMessage);

        showPreviewLink.click(function(e){
          e.preventDefault();
          e.stopPropagation();
          self.change_to_tab('well', ['search','show_preview=true']);
        });
        
        appModel.getModel('apilog', self.model.get('preview_log_key'), 
          function(preview_log){
            showPreviewDetails.attr('title', 
              appModel.print_json(preview_log.pick(
                ['date_time', 'username','user_id','key', 'api_action',
                  'child_logs','comment', 'diff','json_field'])));
          }
        );
        
      };

      showPreviewMessage();
      showReleaseMessage();
      
      self.model.on('sync', showPreviewMessage);
      self.model.on('sync', showReleaseMessage);
    },

    /**
     * Upload Library Contents:
     * (Override)
     * 
     * 
     * NOTE: For a "released" library, all Library Reagent uploads will be 
     * committed to a "preview". 
     * In a "preview" upload:
     * - the server will save ApiLog "preview" logs only;
     * - (all normal updates are rolled back)
     * - patch logs can later be confirmed and applied by the admin.

     * 
     * A preview contains all the upload diff data on 
     * the server and may be viewed, applied ("released"), or deleted by the 
     * Administrator users.
     * 
     * @see appModel.API_PARAM_PATCH_PREVIEW_MODE
     */
    upload: function(e){
      e.preventDefault();
      e.stopPropagation();
      var self = this;
      var url = _.result(this.model, 'url') + '/well';
      var content_types = self.model.resource.content_types; // use standard library content types
      if( this.model.get('screen_type') == 'small_molecule') {
        content_types.push('sdf');
      }

      // Only one upload preview may be stored on the server at one time
      if (_.isNumber(self.model.get('preview_log_id'))){
        appModel.showModalMessage({
          title: 'Library has a pending preview release',
          body: ['Library has a pending preview loaded:',
            'The current preview must be either released or deleted before ',
            'other well data may be loaded'].join('<br/>')
        });
        return;
      }
      
      // Set the param appModel.API_PARAM_PATCH_PREVIEW_MODE
      // as extra post_data to be sent along with the upload.
      // API_PARAM_PATCH_PREVIEW_MODE signals the WellResource to patch as 
      // a preview.
      var options = {};
      if (self.model.get('is_released') == true){
        options[appModel.API_PARAM_PATCH_PREVIEW_MODE] = true;
      }
      
      UploadDataForm.postUploadFileDialog(url, content_types, options)
        .done(function(data, textStatus, jqXHR){
          
          function refreshOnUpdate(){
            // Must refetch and render, once for the containing tabbed layout,
            // and then (during render), will refetch again for the table view
            self.model.fetch({
              reset: true,
            }).done(function() {
              console.log('re-rendering library view...');
              self.setDetail([]);
            }).fail(function() { 
              Iccbl.appModel.jqXHRfail.apply(this,arguments); 
            });
          };
          
          if (appModel.getAppData().get('BACKGROUND_PROCESSING') == true){
            data = _.result(data,appModel.API_RESULT_META, data);
            if (!_.has(data, 'job')){
              appModel.error('Job response does not contain job data');
              return;
            }
            var jobModel = appModel.addBackgroundJob(data['job']);
            
            self.listenTo(jobModel, 'change:state', jobStateChanged);
            
            // TODO: Create a listener to update library from any browser instance
            function jobStateChanged(){
              console.log('jobStateChanged:', arguments);
              if (jobModel.get('state') != 'completed'){
                return;
              }
              if (self.isClosed == true){
                console.log('library view has been closed already');
                return;
              }
              refreshOnUpdate();
            };
          } else {
            refreshOnUpdate();
          }
          
        })
        .fail(function(){
          appModel.jqXHRfail.apply(this,arguments); 
        }
      );
    },

    /**
     * Note: for library, upload and download operate on the library contents, 
     * (reagents) and not on the library entity itself.
     */
    download_contents: function(e){
      e.preventDefault();
      e.stopPropagation();
      var self = this;
      var resource;
      if( this.model.get('screen_type') == 'rnai') {
        resource = appModel.getResource('silencingreagent');
      } else {
        resource = appModel.getResource('smallmoleculereagent');
      }
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'reagent'].join('/');
      url += '?limit=0&includes=*';
      
      appModel.download(url, resource);
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      return {
        'base_url': this.model.resource.key + '/' + this.model.key,
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    setDetail: function(delegateStack) {
      console.log('detail view');
      
      var self = this;
      var key = 'detail';
      var buttons = ['download_contents'];
      if (appModel.hasPermission('library', 'write')){
        buttons = buttons.concat(['download','upload','history','edit']);
      }
      
      // Custom library model validation: plate range
      // FIXME: 20180710 - not working
      this.model.validate = function(attrs) {
        var errs = {};
        console.log('validating: ', attrs); 
        if (!_.isNumber(attrs.start_plate)){
          errs.start_plate = 'number required';
        }else{
          var ranges = _.result(self.model.resource,'library_plate_ranges',[]);
          var rangesAsString = _.map(_.zip(
            _.filter(ranges, function(val, index){ return index%2==0; }),
            _.filter(ranges, function(val, index){ return index%2==1; })),
            function(pair){ return pair.join('-'); } ).join(', ');
          for (var i=0; i<ranges.length-1; i++){
            var start = ranges[i];
            var end = ranges[i+1];
            if (attrs.start_plate >= start && attrs.start_plate <= end){
              errs.start_plate = 'range already used: ' + rangesAsString;
              break;
            }
            if (attrs.end_plate >= start && attrs.end_plate <= end){
              errs.end_plate = 'range already used: ' + rangesAsString;
              break;
            }
            i++;
          }
        }
        if (!_.isEmpty(errs)) return errs;
      };

      var view = this.tabViews[key];
      if (view) {
        // remove the view to refresh the page form
        this.removeView(this.tabViews[key]);
      }
      var detailView = DetailView.extend({
        afterRender: function(){
          DetailView.prototype.afterRender.apply(this,arguments);
          var search_data = {
            key: self.model.key,
            ref_resource_name: self.model.resource.key,
            comment__is_blank: false
          };
          appModel.createCommentTable(self.model,search_data, $('#comment_table'));
          $('#download_contents').click(self.download_contents)
          
        }
      });;

      view = new DetailLayout({ 
        model: this.model,
        uriStack: delegateStack, 
        DetailView: detailView,
        EditView: EditView.extend({
          save_success: function(data, textStatus, jqXHR){
            EditView.prototype.save_success.apply(this,arguments);
            appModel.unset('libraries', {silent: true});
            appModel.getLibraries();
          },
          afterRender: function(){
            EditView.prototype.afterRender.apply(this,arguments);
            var form = this;
            
            form.on('screen_type:change', function(e){
              var screen_type = form.getValue('screen_type');
              if (screen_type == 'small_molecule'){
                form.$el.find('div[data-fields="is_pool"]').hide();
              }else{
                form.$el.find('div[data-fields="is_pool"]').show();
              }
            });
            
          }          
        }),
        buttons: buttons 
      });
      this.tabViews[key] = view;

      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    setWells: function(delegateStack) {
      var self = this;
      
      var resource;
      if( this.model.get('screen_type') == 'rnai') {
        resource = appModel.getResource('silencingreagent');
      } else {
        resource = appModel.getResource('smallmoleculereagent');
      }
      
      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        // Detail view
        // TODO: will do a double fetch of the well here and in wellview
        var well_id = delegateStack.shift();
        appModel.getModel(resource.key, well_id, function(model){
          model.resource = resource;
          view = new LibraryWellView({
            model: model, 
            uriStack: _.clone(delegateStack),
            library: self.model
          });
          self.consumedStack.push(well_id);
          Backbone.Layout.setupView(view);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
          self.$("#tab_container-title").html('Well ' + model.get('well_id'));
    
        });        
        return;
      }else{
        // List view
        var url = [self.model.resource.apiUri, 
                   self.model.key,
                   'reagent'].join('/');
        var download_contents_button = $([
           '<a class="btn btn-default btn-sm pull-right" ',
             'title="Download the full library reagent table (with images, if Excel)" ',
             'role="button" id="download_contents_button" >',
             'Full Download</a>'
           ].join(''));
        download_contents_button.click(self.download_contents)
        view = new LibraryWellsView({ 
          uriStack: _.clone(delegateStack),
          resource: resource,
          url: url,
          library: self.model,
          extraListControls: [download_contents_button]
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
      }
    },
        
    setPlates: function(delegateStack) {
      var self = this;
      var copyPlateResource = appModel.getResource('librarycopyplate'); 

      // List view
      // special url because list is a child view for library
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'plate'].join('/');
      view = new ListView({ 
        uriStack: _.clone(delegateStack),
        resource: copyPlateResource,
        url: url
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    setCopies: function(delegateStack) {
      var self = this;
      var copyResource = appModel.getResource('librarycopy'); 
      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        // Detail view
        var copyname = delegateStack.shift();
        self.consumedStack.push(copyname);
        var _key = [self.model.key,copyname].join('/');
        appModel.getModel(copyResource.key, _key, function(model){
          view = new LibraryCopyView({
            model: model, 
            uriStack: _.clone(delegateStack),
            library: self.model
          });
          Backbone.Layout.setupView(view);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
          
          console.log('title: ', Iccbl.getTitleFromTitleAttribute(model, model.resource));
          self.$("#tab_container-title").html('Copy ' + model.get('copy_name'));
        });        
        return;
      }else{
        // List view
        var extraControls = [];
        var url = [self.model.resource.apiUri, 
                   self.model.key,
                   'copy'].join('/');
        var Collection = Iccbl.MyCollection.extend({
          url: url
        });
        collection = new Collection();
        
        if (appModel.hasPermission(copyResource.key, 'write')){
          var showAddButton = $([
             '<a class="btn btn-default btn-sm pull-down" ',
               'role="button" id="add_resource" href="#">',
               'Add</a>'
             ].join(''));   
           showAddButton.click(function(e){
             e.preventDefault();
             self.showAddCopy(collection);
           });
           extraControls.push(showAddButton);
        }    
        
        var concentration_types = self.model.get('concentration_types');
        if (concentration_types) {
          copyResource.fields['min_molar_concentration']['visibility'] = [];
          copyResource.fields['max_molar_concentration']['visibility'] = [];
          copyResource.fields['min_mg_ml_concentration']['visibility'] = [];
          copyResource.fields['max_mg_ml_concentration']['visibility'] = [];
          if (_.contains(concentration_types, 'mg_ml')){
            copyResource.fields['min_mg_ml_concentration']['visibility'] = ['l','d'];
            copyResource.fields['max_mg_ml_concentration']['visibility'] = ['l','d'];
          }        
          if (_.contains(concentration_types, 'molar')){
            copyResource.fields['min_molar_concentration']['visibility'] = ['l','d'];
            copyResource.fields['max_molar_concentration']['visibility'] = ['l','d'];
          }
        }
        
        copyResource.fields['copy_name'].backgridCellType = 
          Iccbl.LinkCell.extend(_.extend({},
            copyResource.fields['copy_name'].display_options,
            {
              linkCallback: function(e){
                e.preventDefault();
                // re-fetch the full model
                var copyname = this.model.get('copy_name');
                var _key = [self.model.key,copyname].join('/');
                appModel.getModel(copyResource.key, _key, function(model){
                  view = new LibraryCopyView({
                    model: model, 
                    uriStack: [],
                    library: self.model
                  });
                  Backbone.Layout.setupView(view);
                  self.consumedStack = ['copy',copyname];
                  self.listenTo(view , 'uriStack:change', self.reportUriStack);
                  self.setView("#tab_container", view ).render();
                  
                  console.log('title: ', 
                    Iccbl.getTitleFromTitleAttribute(model, model.resource));
                  self.$("#tab_container-title")
                    .html('Copy ' + model.get('copy_name'));
                });        
                return;
              }
            }));
        
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          resource: copyResource,
          url: url,
          collection: collection,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.$("#tab_container-title").empty();
        self.setView("#tab_container", view ).render();
      }
    },

    /**
     * Manually create the copy dialog, to set initial volume/concentration
     */
    showAddCopy: function(collection){
      var self = this;
      var formSchema = {};
      var fieldTemplate = appModel._field_template;
      var formTemplate = appModel._form_template;
      var copyResource = appModel.getResource('librarycopy');
      var copyNameField = _.result(copyResource['fields'], 'copy_name', {});
      var copyUsageTypeField = _.result(copyResource['fields'],'usage_type',{});
      var plateResource = appModel.getResource('librarycopyplate');
      var plateStatusField = _.result(plateResource['fields'],'status',{});
      var plateWellVolumeField = _.result(plateResource['fields'], 'well_volume', {});
      var plateMgMlConcentrationField = _.result(
        plateResource['fields'], 'mg_ml_concentration', {});
      var plateMolarConcentrationField = _.result(
        plateResource['fields'], 'molar_concentration', {});
      var copyUrl = copyResource.apiUri;
        
      formSchema['library_short_name'] = {
        title: 'Library',
        key: 'library_short_name',
        type: EditView.DisabledField,
        template: fieldTemplate 
      };
      
      formSchema['copy_name'] = {
        title: 'Name',
        key: 'copy_name',
        type: Backbone.Form.editors.Text,
        editorClass: 'form-control',
        validators: [
          'required',
          function checkRange(value, formValues) {
            var extantCopies = self.model.get('copies');
            if(extantCopies && _.contains(extantCopies,value)){
              return {
                type: 'Unique',
                message: 'value: ' + value 
                  + ', is already used: copies: ' + extantCopies.join(',')
              };
            }
          }
        ],
        template: fieldTemplate 
      };
      if(_.has(copyNameField,'regex') ){
        formSchema['copy_name']['validators'].push(
          { type: 'regexp', 
            regexp: new RegExp(copyNameField['regex']),
            message: _.result(
              copyNameField,'validation_message','Name pattern is incorrect' )
          }
        );
      }
      
      formSchema['usage_type'] = {
        title: 'Usage Type',
        key: 'usage_type',
        type: EditView.ChosenSelect,
        editorClass: 'chosen-select',
        editorAttrs: { widthClass: 'col-sm-5'},
        validators: ['required'],
        options: appModel.getVocabularySelectOptions(
          copyUsageTypeField.vocabulary_scope_ref),
        template: fieldTemplate 
      };
      
      formSchema['initial_plate_status'] = {
        title: 'Initial Plate Status',
        help: 'Initial status for the plates of this copy (optional)',
        key: 'initial_plate_status',
        type: EditView.ChosenSelect,
        editorClass: 'chosen-select',
        editorAttrs: { widthClass: 'col-sm-5'},
        validators: [],
        options: appModel.getVocabularySelectOptions(
          plateStatusField.vocabulary_scope_ref),
        template: fieldTemplate 
      };
      
      formSchema['initial_plate_well_volume'] = {
        title: 'Initial Plate Well Volume',
        key: 'initial_plate_well_volume',
        validators: ['required',EditView.CheckPositiveNonZeroValidator],
        type: EditView.SIunitEditor,
        template: fieldTemplate 
      };
      _.extend(
        formSchema['initial_plate_well_volume'],
        plateWellVolumeField['display_options']);

      formSchema['set_initial_plate_concentration'] = {
        title: '<strong>Set initial plate concentrations:</strong>',
        help: 'Select to override the default library well concentrations',
        key: 'set_initial_plate_concentration',
        type: 'Checkbox',
        template: appModel._alt_checkbox_template
      };
      formSchema['initial_plate_mg_ml_concentration'] = {
        title: 'Initial Plate mg/ml Concentration',
        help: 'If set, this value overrides the default library well concentration',
        key: 'initial_plate_mg_ml_concentration',
        type: Backbone.Form.editors.Number,
        validators: [EditView.CheckPositiveNonZeroValidator],
        editorClass: 'form-control',
        template: fieldTemplate 
      };
      _.extend(
        formSchema['initial_plate_mg_ml_concentration'],
        plateMgMlConcentrationField['display_options']);

      formSchema['initial_plate_molar_concentration'] = {
        title: 'Initial Plate Molar Concentration',
        help: 'If set, this value overrides the default library well concentration',
        key: 'initial_plate_molar_concentration',
        type: EditView.SIunitEditor,
        validators: [EditView.CheckPositiveNonZeroValidator],
        editorClass: 'form-control form-control-zero-padding',
        template: fieldTemplate 
      };
      _.extend(
        formSchema['initial_plate_molar_concentration'],
        plateMolarConcentrationField['display_options']);
      
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        type: 'TextArea',
        editorClass: 'input-full',
        template: fieldTemplate
      };
      
      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs) {
          var errors = {};
          if (!_.isNull(attrs['initial_plate_mg_ml_concentration'])
              && !_.isNull(attrs['initial_plate_molar_concentration'])){
            var msg = 'Must enter either (mg/ml) or (molar)';
            errors['initial_plate_mg_ml_concentration'] = msg;
            errors['initial_plate_molar_concentration'] = msg;
            
          }
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var formFields = new FormFields({
        'initial_plate_mg_ml_concentration': null, 
        'initial_plate_molar_concentration': null, 
        'initial_plate_well_volume': null
      });
      
      var form = new Backbone.Form({
        model: formFields,
        template: formTemplate
      });
      form.setValue('library_short_name', self.model.get('short_name'));
      var formview = form.render();
      var _form_el = formview.el;
      var $form = formview.$el;

      $form.find('[name="initial_plate_mg_ml_concentration"]')
        .prop('disabled', true);
      $form.find('[name="initial_plate_molar_concentration"]')
        .find('input, select').prop('disabled', true);

      form.listenTo(form, "change", function(e){
        var set_initial_plate_concentration = 
            form.getValue('set_initial_plate_concentration');
        if(set_initial_plate_concentration){
          $form.find('[name="initial_plate_mg_ml_concentration"]')
            .prop('disabled', false);
          $form.find('[name="initial_plate_molar_concentration"]')
            .find('input, select').prop('disabled', false);
        } else {
          $form.find('[name="initial_plate_mg_ml_concentration"]')
            .prop('disabled', true);
          $form.find('[name="initial_plate_molar_concentration"]')
            .find('input, select').prop('disabled', true);
        }
      });
      
      $form.find('.chosen-select').chosen({
        disable_search_threshold: 3,
        width: '100%',
        allow_single_deselect: true,
        search_contains: true
        });

      var dialog = appModel.showModal({
        okText: 'Submit',
        view: _form_el,
        title: 'Create a new copy for library: ' + self.model.get('library_name'),
        ok: function(e) {
          e.preventDefault();
          var errors = form.commit({ validate: true }) || {};
          var values = form.getValue();
          
          form.$el.find('.form-group').removeClass('has-error');
          if (!_.isEmpty(errors) ) {
            _.each(_.keys(errors), function(key) {
              form.$el.find('[name="'+key +'"]')
                .parents('.form-group').addClass('has-error');
            });
            return false;
          }            
          
          if (! values['set_initial_plate_concentration']) {
            delete values['initial_plate_mg_ml_concentration'];
            delete values['initial_plate_molar_concentration'];
          } else {
            if (! values['initial_plate_mg_ml_concentration']) {
              delete values['initial_plate_mg_ml_concentration'];
            }
            if (!values['initial_plate_molar_concentration']) {
              delete values['initial_plate_molar_concentration'];
            }
          }
          
          $.ajax({
            url: copyUrl,    
            data: JSON.stringify(values),
            cache: false,
            contentType: 'application/json',
            dataType: 'json', // what is expected back from the server
            processData: false, // do not process data being sent
            type: 'PATCH'
          })
          .fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); })
          .done(function(data) {
            self.model.fetch();
            appModel.showConnectionResult(data, {
              title: 'Copy Created'
            });

            collection.fetch();
          }); // ajax
        } // ok
      }); // dialog
      
    }
  
  });
  
  return LibraryView;
});