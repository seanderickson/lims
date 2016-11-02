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
  'views/generic_detail_layout',
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'views/list2',
  'utils/uploadDataForm',
  'utils/tabbedController'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, LibraryCopyView, 
            LibraryWellView, DetailLayout, DetailView, EditView, ListView, 
            UploadDataForm, TabbedController ) {

  var LibraryView = TabbedController.extend({
    
    initialize: function(args) {
      var self = this;
      this._classname = 'LibraryView';
      _.bindAll(this, 'upload');
      TabbedController.prototype.initialize.apply(this,arguments);
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
      'click button#upload': 'upload'        
    },

    upload: function(event){
      var self = this;
      event.preventDefault();
      var url = _.result(this.model, 'url') + '/well';
      var content_types = self.model.resource.content_types; // use standard library content types
      if( this.model.get('screen_type') == 'small_molecule') {
        content_types.push('sdf');
      }
      UploadDataForm.postUploadFileDialog(url, content_types)
        .done(function(){
          self.model.fetch();
        })
        .fail(function(){
          appModel.jqXHRfail.apply(this,arguments); 
        }
      );
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
      var buttons = ['download'];
      if (appModel.hasPermission('library', 'write')){
        buttons = buttons.concat(['upload','history','edit']);
      }
      
      // Custom library model validation: plate range
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
          self.createCommentHistoryTable($('#comments'), this);
        },
      });

      view = new DetailLayout({ 
        model: this.model,
        uriStack: delegateStack, 
        DetailView: detailView,
        buttons: buttons 
      });
      this.tabViews[key] = view;

      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    /**
     * Update the libryary with a comment table: populate
     * using the apilog history of the library
     **/
    createCommentHistoryTable: function($target_el, view){
      console.log('create the comment history table', $target_el);
      var self = this;
      var apilogResource = appModel.getResource('apilog');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: apilogResource.apiUri 
      });
      var cell = $([
        '<div id="comment_table_div" class="row">',
        '</div>'
        ].join(''));
      
      var $addCommentButton= $(
        '<div class="col-xs-12"><button class="btn btn-default btn-sm" role="button" \
        id="addCommentButton" >Add comment</button></div>');

      if (appModel.hasPermission('library', 'write')){
        cell.append($addCommentButton);
      }
      
      function build_table(collection){
        if (collection.isEmpty()){
          return;
        }
        $target_el.html(cell);
        
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
              'name' : 'comment',
              'label' : 'Comment',
              'description' : 'Comment',
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
            }),
            _.extend({},colTemplate,{
              'name' : 'username',
              'label' : 'Username',
              'description' : 'Username',
              'order': 1,
              'sortable': true,
              'cell': Iccbl.LinkCell.extend({
                'hrefTemplate': '#user/{username}'
              })
            })
        ];
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();

        cell.append($('<div class="col-xs-12" id="comment_items"/>'));
        
        var comment_grid= new Backgrid.Grid({
          columns: colModel,
          collection: collection,
          className: 'backgrid table-striped table-condensed table-hover '
        });
        cell.find('#comment_items').html(comment_grid.render().$el);
      }
      
      var comment_collection = new CollectionClass();
      comment_collection.fetch({
        data: { 
          limit: 0,
          key: self.model.get('short_name'),
          ref_resource_name: self.model.resource.key,
          comment__is_null: false,
          diff_keys__is_null: true,
          order_by: ['-date_time']
        },
        success: build_table
      }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      
      $addCommentButton.click(function(){

        var FormFields = Backbone.Model.extend({
          schema: {
            comment: {
              title: 'Comment',
              key: 'comment',
              type: 'TextArea',
              editorClass: 'input-full',
              validators: ['required'], 
              template: appModel._field_template
            }
          }
        });
        var formFields = new FormFields();
        var form = new Backbone.Form({
          model: formFields,
          template: appModel._form_template
        });
        var _form_el = form.render().el;

        appModel.showModal({
          okText: 'ok',
          ok: function(e){
            e.preventDefault();
            
            var errors = form.commit();
            if(!_.isEmpty(errors)){
              console.log('form errors, abort submit: ' + JSON.stringify(errors));
              return false;
            }
            var values = form.getValue();
            
            // TODO: submit a PATCH with only a header comment
            var comment = values['comment'];
            var headers = {};
            headers[appModel.HEADER_APILOG_COMMENT] = comment;
            $.ajax({
              url: [self.model.resource.apiUri, self.model.key].join('/'),    
              // PATCH with no data; 
              // data: JSON.stringify({ comment: comment }),  
              cache: false,
              contentType: 'application/json',
              dataType: 'json', // what is expected back from the server
              processData: false, // do not process data being sent
              type: 'PATCH',
              headers: headers, 
              success: function(data){
                console.log('success', data);
                self.model.fetch({ reset: true });
                view.render();
              },
              done: function(model, resp){
                // TODO: done replaces success as of jq 1.8
                console.log('done');
              }
            }).fail(function(){ appModel.jqXHRfail.apply(this,arguments); });
          
            return true;
            
          },
          view: _form_el,
          title: 'Add comment'  
        });
        
      });
    },
    
    setWells: function(delegateStack) {
      var self = this;
      var schemaUrl = [self.model.resource.apiUri,self.model.key,'well',
                       'schema'].join('/');
      var wellResource = appModel.getResource('well'); 

      function createWellView(schemaResult){
        if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
            !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
          // Detail view
          var well_id = delegateStack.shift();
          self.consumedStack.push(well_id);
          var _key = [self.model.key,well_id].join('/');
          appModel.getModel(wellResource.key, well_id, function(model){
            model.resource = schemaResult;
            view = new LibraryWellView({
              model: model, 
              uriStack: _.clone(delegateStack),
              library: self.model
            });
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
                     'well'].join('/');
          view = new ListView({ 
            uriStack: _.clone(delegateStack),
            schemaResult: schemaResult,
            resource: schemaResult,
            url: url
          });
          Backbone.Layout.setupView(view);
          self.reportUriStack([]);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
        }
        
      }
      appModel.getResourceFromUrl(schemaUrl, createWellView);

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
        schemaResult: copyPlateResource,
        resource: copyPlateResource,
        url: url
      });
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
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
          self.$("#tab_container-title").html('Copy ' + model.get('name'));
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
        collection = new Collection({
          url: url,
        });
        
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
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: copyResource,
          resource: copyResource,
          url: url,
          collection: collection,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
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
      var copyNameField = _.result(copyResource['fields'], 'name', {});
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
      
      formSchema['name'] = {
        title: 'Name',
        key: 'name',
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
        formSchema['name']['validators'].push(
          { type: 'regexp', 
            regexp: new RegExp(copyNameField['regex']),
            message: _.result(copyNameField,'validation_message','Name pattern is incorrect' )
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
        options: appModel.getVocabularySelectOptions(copyUsageTypeField.vocabulary_scope_ref),
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
        options: appModel.getVocabularySelectOptions(plateStatusField.vocabulary_scope_ref),
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
      
// TODO:
//      formSchema['initial_plate_status'] = {
//        title: 'Initial Plate Status',
//        key: 'initial_plate_status',
//        type: EditView.ChosenSelect,
//        editorClass: 'chosen-select',
//        editorAttrs: { widthClass: 'col-sm-5'},
//        options: appModel.getVocabularySelectOptions(plateStatusField.vocabulary_scope_ref),
//        template: fieldTemplate 
//      };
      
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

      $form.find('[name="initial_plate_mg_ml_concentration"]').prop('disabled', true);
      $form.find('[name="initial_plate_molar_concentration"]').find('input, select').prop('disabled', true);

      form.listenTo(form, "change", function(e){
        console.log('change');
        var set_initial_plate_concentration = form.getValue('set_initial_plate_concentration');
        if(set_initial_plate_concentration){
          $form.find('[name="initial_plate_mg_ml_concentration"]').prop('disabled', false);
          $form.find('[name="initial_plate_molar_concentration"]').find('input, select').prop('disabled', false);
        } else {
          $form.find('[name="initial_plate_mg_ml_concentration"]').prop('disabled', true);
          $form.find('[name="initial_plate_molar_concentration"]').find('input, select').prop('disabled', true);
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
            collection.fetch();
            
            // TODO: refactor - display response
            if (_.isObject(data) && !_.isString(data)) {
              data = _.result(_.result(data,'meta',data),'Result',data);
              var msg_rows = appModel.dict_to_rows(data);
              var bodyMsg = msg_rows;
              if (_.isArray(msg_rows) && msg_rows.length > 1) {
                bodyMsg = _.map(msg_rows, function(msg_row) {
                  return msg_row.join(': ');
                }).join('<br>');
              }
              var title = 'Copy created';
              appModel.showModalMessage({
                body: bodyMsg,
                title: title  
              });
            } else {
              console.log('ajax should have been parsed data as json', data);
              appModel.showModalMessage({
                title: 'Copy created',
                okText: 'ok',
                body: 'Copy name: ' + values['name']
              });
            }
          }); // ajax
        } // ok
      }); // dialog
      
    }
  
  });
  
  return LibraryView;
});