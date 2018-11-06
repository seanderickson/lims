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
  'utils/wellSelector2',
  'templates/genericResource.html',
  'google_palette'
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, 
         DetailView, EditView, WellSelector, genericLayout, palette) {
  
  /**
   * Transform plate matrix data into a plate-well spreadsheet
   */
  var RawDataTransformView = Backbone.Layout.extend({
    
    template: _.template(genericLayout),

    ERR_MSG_: 'Err',
    
    // all - red, positive_controls -red, negative -royal blue, other -mint green, library -yellow
    // '#FF0000',
    CONTROL_WELLS_COLOR_PALETTE: ['#e31a1c', '#2B60DE', '#98FF98', '#FFFF00'],
    
    initialize: function(args) {
      console.log('initialize', args, arguments);
      var self = this;
      self.resource = self.model.resource;
      self.screen = args.screen;
      self.cherry_pick_request = args.cherry_pick_request;
      self.uriStack = args.uriStack;
      self.consumedStack = [];
      _.bindAll(this, 'submitUpload');
      
      self.show_library_controls = true;
      if (self.cherry_pick_request){
        self.show_library_controls = false;
      }
      if (self.screen && self.screen.get('screen_type') == 'small_molecule'){
        self.show_library_controls = false;
      }
      
      // Track the input file metadata
      var InputFileCollection = Backbone.Collection.extend({
        // explicitly define the id so that collection compare & equals work
        modelId: function(attrs) {
          return _.result( attrs, 'ordinal');
        },
        comparator: 'ordinal'
      })
      self.inputFileCollection = new InputFileCollection(this.model.get('input_files'));
      
      if (self.inputFileCollection.isEmpty()){
        self.inputFileCollection.add(
          new Backbone.Model(_.extend(
              appModel.get_field_defaults(self.resource.input_file_fields),
              { ordinal: 0 }))
        );
      }
      self.downloadButton = $([
        '<button type="button" class="btn btn-default btn-xs" ',
        ' title="Download the file with the transformed data, if upload has been completed" ',
        'id="download_button" >download</button>',
      ].join(''));
          
      self.downloadButton.click(function(e){
        e.preventDefault();
        console.log('download transformed file...');
        if (self.screen){
          appModel.downloadUrl(
            '/db/screen_raw_data_transform/' + self.screen.get('facility_id'));
        }else if (self.cherry_pick_request){
          appModel.downloadUrl(
            '/db/cpr_raw_data_transform/' 
              + self.cherry_pick_request.get('cherry_pick_request_id'));
        }
      });
      self.uploadButton = $([
        '<button type="button" class="btn btn-default btn-xs" ',
        ' title="Upload and transform the set of input files that have been specified below" ',
        'id="upload_button" >upload</button>',
      ].join(''));
      self.uploadButton.click(self.submitUpload);
      self.addInputFileButton = $([
        '<button type="button" class="btn btn-default btn-xs controls-right" ',
        'id="add_input_file" >add input file</button>',
      ].join(''));
      
    },

    /**
     * Submit the raw data tranform data:
     * NOTE: the result of the submission may contain error messages; for this
     * reason, the transformed file is not returned as a result of the POST.
     * Ajax cannot easily decode the xlsx file response and/or process error
     * messages at the same time.
     * However, if the upload is successful, a download is triggered with
     * the receipt of the success response.
     */
    submitUpload: function() {
      var self = this;
      console.log('Transform form uploadData...');
      
      var url;
      if (self.screen){
        url = [self.resource.apiUri, self.screen.key].join('/');
      }else if (self.cherry_pick_request){
        url = [self.resource.apiUri, 
               self.cherry_pick_request.get('screen_facility_id'),
               self.cherry_pick_request.get('cherry_pick_request_id')].join('/');
      }

      // find the input files
      var input_errors = [];
      var formData = new FormData();
      self.inputFileCollection.each(function(input_file_model){
        console.log('validate input_file_model', input_file_model);

        var form = input_file_model.form;
        
        var errors = form.commit({ validate: true }); // runs schema and model validation
        console.log('form commit returns ', errors);
        if(!_.isEmpty(errors) ){
          console.log('form errors, abort submit: ',errors);
          _.each(_.keys(errors), function(key){
            form.$el.find('[name="'+key +'"]').parents('.form-group').addClass('has-error');
          });
          if (_.has(errors,'_others')){
            var $errorDiv = $('<div id="general_error" class="panel text-danger" />');
            _.each(errors['_others'], function(otherError){
              _.each(_.values(otherError), function(errMsg){
                $errorDiv.append('Error: ' + errMsg );
              });
            });
            form.$el.append($errorDiv);
          }
          input_errors = input_errors.concat(errors)
          return false;
        }
        var file = form.$el.find('input[name="file_input"]')[0].files[0];
        if (!file){
          throw 'No file for form ' + ordinal;
        }
        var input_file_prefix = 'input_file_'+ input_file_model.get('ordinal');
        formData.append(input_file_prefix,file);
        var values = form.getValue();
        _.each(_.keys(values), function(key){
            formData.append(input_file_prefix + '_' + key,values[key])
        });
      });
      
      var errs = self.view.commit({ validate: true });
      console.log('input_errors', input_errors);
      if (!_.isEmpty(input_errors)){
        return false;
      }
      if (!_.isEmpty(errs)){
        return false;
      }
      
      var values = self.view.getValue();
      _.each(_.keys(values), function(key){
        formData.append(key,values[key])
      });
      for(var pair of formData.entries()) {
        console.log('form value:' , pair[0], pair[1]); 
      }
      $.ajax({
        url:  url,
        data: formData,
        cache: false,
        contentType: false,
        processData: false,
        dataType: 'json', // what is expected back from the server
        type: 'POST',
      }).done(function(data, textStatus, jqXHR){
        console.log('success', data);
        var meta = _.result(data,appModel.API_RESULT_META, {});
        appModel.showJsonMessages(meta);
        self.downloadButton.click();
      }).fail(function(jqXHR, textStatus, errorThrown) { 
        console.log('errors', arguments);
        self.view.save_fail.apply(self.view, arguments);
      });
    },
    
    afterRender: function(){
      var self = this;
      console.log('edit view');
      
      var vocab_scope_ref = 
        self.model.resource.fields['assay_plate_size']['vocabulary_scope_ref'];
      var vocab = Iccbl.appModel.getVocabularySelectOptions(vocab_scope_ref);
      vocab = _.reject(vocab, function(obj){
        return obj.val == 1;
      });
      this.model.resource.fields['assay_plate_size'].choiceHash = vocab;
      vocab_scope_ref = 
        self.model.resource.fields['library_plate_size']['vocabulary_scope_ref'];
      vocab = Iccbl.appModel.getVocabularySelectOptions(vocab_scope_ref);
      vocab = _.reject(vocab, function(obj){
        return obj.val != 384;
      });
      this.model.resource.fields['library_plate_size'].choiceHash = vocab;

      if (self.show_library_controls == false){
        self.model.resource.fields['library_controls']['editability'] = [];
      }
      
      // TODO: Create a separate TransformerEditView class
      var TransformerEditView = EditView.extend({
        
        afterRender: function(){
          var editForm = this;
          console.log('after render...');
          editForm.$el.find('#generic-edit-buttons-top').empty();
          editForm.$el.find('#generic-edit-buttons-bottom').empty();

          editForm.$el.find('#generic-edit-buttons-top').append(self.uploadButton);
          editForm.$el.find('#generic-edit-buttons-top').append(self.downloadButton);
          
          var inputFilesDiv = $([
            '<div class="panel panel-default">',
            '<div class="panel-heading">Input Files: </div>',
            '<div class="panel-body">',
            '<div class="" id="input_file_form_area" ></div>',
            '</div>',
            '</div>'].join(''));
          inputFilesDiv.find('.panel-heading').append(self.addInputFileButton);
          var inputFileArea = inputFilesDiv.find('#input_file_form_area');
          editForm.$el.find('#generic-edit-buttons-bottom').append(inputFilesDiv);
          
          self.inputFileCollection.each(function(input_file_model){
            console.log('input_file', input_file_model);
            self._create_input_file_form(inputFileArea, input_file_model);
          });
          
          self.addInputFileButton.click(function(){
            var ordinal = self.inputFileCollection.size();
            var input_file_model = new Backbone.Model(_.extend(
              appModel.get_field_defaults(self.resource.input_file_fields),
              { ordinal: ordinal }));
            self.inputFileCollection.add(input_file_model);
            self._create_input_file_form(inputFileArea, input_file_model);
          });
          
          self.setupWellControlSelectionButtons(editForm);
          
          editForm.$el.find('[name="plate_ranges"]').change(function(){
            var errors = editForm.commit({ validate: true }); 
            // runs schema and model validation, errors not needed
            var plate_ranges = editForm.getValue('plate_ranges');
            if (_.isEmpty(plate_ranges)) return;
            var plate_range_string = plate_ranges.split(/\s*,\s*/).join('_');
            if (self.screen){
                editForm.setValue('output_filename', 
                  self.screen.get('facility_id') + '_' + plate_range_string);
            }else if (self.cherry_pick_request){
                editForm.setValue('output_filename', 
                  'cpr' + self.cherry_pick_request.get('cherry_pick_request_id')
                  + '_' + plate_range_string);
            }
            // run validation again, to clear filename error
            var errors = editForm.commit({ validate: true }); 
          });
      
          EditView.prototype.afterRender.apply(this,arguments);

          if (!self.model.isNew()){
            console.log('perform initial validation...');
            var errs = this.commit({ validate: true });
            console.log('errs', errs);
          }
    
        },
        commit: function(options) {
          self.inputFileCollection.each(function(input_file_model){
            var $errorDiv = input_file_model.form.$el.find('#general_errors');
            $errorDiv.empty();
            $errorDiv.hide();
          });
          return TransformerEditView.__super__.commit.apply(this, arguments);
        },

        process_errors: function(errors) {
          
          var non_field_errors = TransformerEditView.__super__.process_errors.call(this, errors);
          var new_non_field_errors = {};
          _.each(_.keys(non_field_errors), function(key) {
            var error = errors[key];
            
            var match = key.match(/input_file_(\d)+/);
            
            if (match){
              console.log('errors for file', key, error);
              var ordinal = parseInt(match[1]);
              var input_file_model = self.inputFileCollection.findWhere({'ordinal': ordinal});
              
              if (input_file_model){
                var $errorDiv = input_file_model.form.$el.find('#general_errors');
                $errorDiv.empty();
                $errorDiv.show();
                $errorDiv.append(appModel.print_dict(error));
              }else{
                new_non_field_errors[key] = error;
              }
            }else{
              new_non_field_errors[key] = error;
            }
          });
          return new_non_field_errors;
        }
      }); // TransformerEditView definition
      
      function assayControlWellOverlapValidator(value, formValues){
        console.log('assayControlWellOverlapValidator', value, formValues);
        var assayPlateSize = formValues['assay_plate_size'];
        var wellHash = {};
  
        var parseErrors = [];
        
        var currentNamedRanges = Iccbl.parseNamedWellRanges(value,assayPlateSize,parseErrors);
        var currentWells = _.flatten(_.map(_.values(currentNamedRanges), function(nwr){
          return nwr['wells'];
        }));
        var allAssayControlWellRanges = self.getAllAssayControlWellRanges(formValues, parseErrors);
        // parse errors handled in validate
        var duplicates = 
          Iccbl.find_duplicates(_.map(_.values(allAssayControlWellRanges),
            function(nwr){
              return nwr['wells'];
            }
          ));
        if (!_.isEmpty(_.intersection(currentWells, duplicates))){
          console.log('duplicates error...', value, currentNamedRanges, duplicates);
          return{
            type: 'duplicates',
            message: 'wells are shared with other control types: ' + duplicates.join(', ')
          };
        }
      };
      
      function libraryControlWellOverlapValidator(value, formValues){
        console.log('assayControlWellOverlapValidator', value, formValues);
        var assayPlateSize = formValues['assay_plate_size'];
        var libraryPlateSize = formValues['library_plate_size'];
        var wellHash = {};
  
        var parseErrors = [];
        
        var libraryControlNamedRanges = 
          Iccbl.parseNamedWellRanges(value,libraryPlateSize,parseErrors);
        var currentWells = _.flatten(_.map(_.values(libraryControlNamedRanges), 
          function(nwr){
            return nwr['wells'];
          }));
        var allAssayControlWellRanges = 
          self.getAllAssayControlWellRanges(formValues, parseErrors);
        var assayControlWells = _.flatten(_.map(_.values(allAssayControlWellRanges), 
          function(nwr){
            return nwr['wells'];
          }));
        if (assayPlateSize < libraryPlateSize){
          assayControlWells = Iccbl.convoluteWells(
            assayPlateSize,libraryPlateSize,assayControlWells);
        } else if (assayPlateSize > libraryPlateSize) {
          var assayControlWellsByQuadrant = 
            Iccbl.deconvoluteWells(assayPlateSize,libraryPlateSize,assayControlWells);
          assayControlWells = _.unique(_.flatten(_.values(assayControlWellsByQuadrant)));
        }
        var duplicates = _.intersection(currentWells, assayControlWells);
        if (!_.isEmpty(duplicates)){
          console.log('duplicates error...', value, currentWells, assayControlWells, duplicates);
          return{
            type: 'duplicates',
            message: 'wells are shared with other control types: ' + duplicates.join(', ')
          };
        }
      };
      function assayWellSelectionValidator(value, formValues){
        var errors = [];
        var assayPlateSize = formValues['assay_plate_size'];
        console.log('perform well selection validation', value, assayPlateSize);
        Iccbl.parseNamedWellRanges(value,assayPlateSize,errors);
        if (!_.isEmpty(errors)){
          console.log('wsv error', errors);
          return {
            type: 'parse',
            message: errors.join('; ')
          };
        }
      };
      function libraryWellSelectionValidator(value, formValues){
        var errors = [];
        var libraryPlateSize = formValues['library_plate_size'];
        console.log('perform well selection validation', value, libraryPlateSize);
        Iccbl.parseNamedWellRanges(value,libraryPlateSize,errors);
        if (!_.isEmpty(errors)){
          console.log('wsv error', errors);
          return {
            type: 'parse',
            message: errors.join('; ')
          };
        }
      };
      function plateRangesValidator(value, formValues){
        var errors = [];
        var platesAndRanges = value.split(/\s*,\s*/);
        var plates = [];
        _.each(platesAndRanges, function(pr){
          if (Iccbl.PLATE_RANGE_PATTERN.test(pr)){
            var rangeParts = Iccbl.PLATE_RANGE_PATTERN.exec(pr);
            plates = plates.concat([rangeParts[1],rangeParts[2]]);
          }else{
            plates.push(pr);
          }
        })
        _.each(plates,function(plate){
          if (!Iccbl.PLATE_PATTERN.exec(plate)){
            errors.push('Not a plate number: ' + plate);
          }else if (self.screen){
            // TODO: check plate numbers?
          }else if (self.cherry_pick_request){
            if (parseInt(plate) > self.cherry_pick_request.get('number_plates_completed')){
              errors.push('Plate is out of range: ' + plate 
                + ', max: ' + self.cherry_pick_request.get('number_plates_completed'));
            }
          }
        })
        if (!_.isEmpty(errors)){
          console.log('plates error', errors);
          return {
            type: 'parse',
            message: errors.join('; ')
          };
        }
        
      }
      
      var view = self.view = new TransformerEditView({ 
        model: this.model,
        resource: this.model.resource,
        uriStack: self.uriStack,
        editSchemaOverrides: {
          'assay_positive_controls': {
            validators: [assayWellSelectionValidator, assayControlWellOverlapValidator],
            editorAttrs: { widthClass: 'col-sm-6'}
          },
          'assay_negative_controls': {
            validators: [assayWellSelectionValidator, assayControlWellOverlapValidator],
            editorAttrs: { widthClass: 'col-sm-6'}
          },
          'assay_other_controls': {
            validators: [assayWellSelectionValidator, assayControlWellOverlapValidator],
            editorAttrs: { widthClass: 'col-sm-6'}
          },
          'library_controls': {
            validators: [libraryWellSelectionValidator, libraryControlWellOverlapValidator],
            editorAttrs: { widthClass: 'col-sm-6'}
          },
          'plate_ranges': {
            validators: ['required', plateRangesValidator]
          }
        }
      });
      Backbone.Layout.setupView(view);
      editform_el = this.setView("#resource_content", view ).render().$el;
      self.reportUriStack([]);
      
    }, // afterRender

    getAllAssayControlWellRanges: function(formValues, errors){
      var self = this;
      var controlWellNamedRanges = {};
      var assayPlateSize = formValues['assay_plate_size'];
      var ordinal = 1;
      _.each(['assay_positive_controls','assay_negative_controls',
              'assay_other_controls'], function(control_type){
        var value = formValues[control_type];
        var namedRanges = Iccbl.parseNamedWellRanges(value,assayPlateSize,errors);
        console.log('label', control_type,  namedRanges);
        var wells = [];
        _.each(_.values(namedRanges), function(namedRange){
          wells = wells.concat(namedRange['wells']);
        });
        var label = self.resource.fields[control_type].title;
        controlWellNamedRanges[label] = { 
          ordinal: ordinal, 
          wells: wells,
          label: label
        }
        ordinal += 1;
      });
      
      console.log('controlWellNamedRanges', controlWellNamedRanges);
      return controlWellNamedRanges;
    },
    
    /** Setup the assay controls well selector buttons **/
    setupWellControlSelectionButtons: function(editForm){
      var self = this;
      var select_assay_positives_button = $([
        '<a class="btn btn-default btn-sm" ',
          'role="button" id="assay_positive_controls_button" >',
          'select</a>'
        ].join(''));
      var select_assay_negatives_button = $([
        '<a class="btn btn-default btn-sm" ',
          'role="button" id="assay_negative_controls_button" >',
          'select</a>'
        ].join(''));
      var select_assay_other_button = $([
        '<a class="btn btn-default btn-sm" ',
          'role="button" id="assay_other_controls_button" >',
          'select</a>'
        ].join(''));
      var select_library_controls_button = $([
        '<a class="btn btn-default btn-sm" ',
          'role="button" id="library_controls_button" >',
          'select</a>'
        ].join(''));
      
      function getCurrentAssayPlateSize(){
        var assayPlateSize = editForm.$el.find('[name="assay_plate_size"]').val();
        if (_.isEmpty(assayPlateSize)){
          assayPlateSize = self.model.get('assay_plate_size');
        }
        if (_.isEmpty(assayPlateSize)){
          appModel.showModalError('must define the assay_plate_size');
          return;
        } else {
          return assayPlateSize;
        }
      };
      function getCurrentLibraryPlateSize(){
        var libraryPlateSize = editForm.$el.find('[name="library_plate_size"]').val();
        if (_.isEmpty(libraryPlateSize)){
          libraryPlateSize = self.model.get('library_plate_size');
        }
        if (_.isEmpty(libraryPlateSize)){
          appModel.showModalError('must define the library plate size');
          return;
        } else {
          return libraryPlateSize;
        }
      };
      function assayControlsButtonHandler(fieldKey){
        var assayPlateSize = getCurrentAssayPlateSize();
        if (!assayPlateSize) return;
        var label = self.resource.fields[fieldKey].title;
        var currentValue = editForm.$el.find('[name="'+ fieldKey + '"]').val();
        self.showWellSelector(
          currentValue, assayPlateSize,'Select ' + label, true,
          function(newValue){
            editForm.setValue(fieldKey, newValue);
            self.model.set(fieldKey, newValue);
            editForm.commit({ validate: true });
          });
      };
      function libraryControlsButtonHandler(fieldKey){
        var libraryPlateSize = getCurrentLibraryPlateSize();
        if (!libraryPlateSize) return;
        var label = self.resource.fields[fieldKey].title;
        var currentValue = editForm.$el.find('[name="'+ fieldKey + '"]').val();
        self.showWellSelector(
          currentValue, libraryPlateSize,'Select ' + label, true,
          function(newValue){
            editForm.setValue(fieldKey, newValue);
            self.model.set(fieldKey, newValue);
            editForm.commit({ validate: true }); 
          });
      };
      
      select_assay_positives_button.click(
        _.partial(assayControlsButtonHandler,'assay_positive_controls'));
      editForm.$el.find('div[data-fields="assay_positive_controls"]')
        .find('.form-group').append(select_assay_positives_button);
      editForm.$el.find('[name="assay_positive_controls"]').change(function(){
        // runs schema and model validation, errors not needed
        var errors = editForm.commit({ validate: true }); 
      });
      
      select_assay_negatives_button.click(
        _.partial(assayControlsButtonHandler,'assay_negative_controls'));
      editForm.$el.find('div[data-fields="assay_negative_controls"]')
        .find('.form-group').append(select_assay_negatives_button);
      editForm.$el.find('[name="assay_negative_controls"]').change(function(){
        // runs schema and model validation, errors not needed
        var errors = editForm.commit({ validate: true }); 
      });

      select_assay_other_button.click(
        _.partial(assayControlsButtonHandler,'assay_other_controls'));
      editForm.$el.find('div[data-fields="assay_other_controls"]')
        .find('.form-group').append(select_assay_other_button);
      editForm.$el.find('[name="assay_other_controls"]').change(function(){
        // runs schema and model validation, errors not needed
        var errors = editForm.commit({ validate: true }); 
      });
      
      if (self.show_library_controls === true){
        select_library_controls_button.click(
          _.partial(libraryControlsButtonHandler,'library_controls'));
        editForm.$el.find('div[data-fields="library_controls"]')
          .find('.form-group').append(select_library_controls_button);
        editForm.$el.find('[name="library_controls"]').change(function(){
          // runs schema and model validation, errors not needed
          var errors = editForm.commit({ validate: true }); 
          // errors not needed; commit will cause to display
        });
      }
      editForm.$el.find('[name="assay_plate_size"]').change(function(){
        // runs schema and model validation, errors not needed
        var errors = editForm.commit({ validate: true }); 
      });
      editForm.$el.find('[name="library_plate_size"]').change(function(){
        // runs schema and model validation, errors not needed
        var errors = editForm.commit({ validate: true }); 
      });

      // Show combined assay controls
      var show_all_assay_controls_button = $([
        '<a class="btn btn-default btn-sm" ',
          ' title="Show a combined view of assay positive, negative, and other controls" ',
          'role="button" id="show_all_assay_controls_button" >',
          'Show All Assay Controls</a>'
        ].join(''));
      
      show_all_assay_controls_button.click(function(e){
        e.preventDefault();
        editForm.commit({ validate: true }); 
        var values = editForm.getValue();
        var assayPlateSize = values['assay_plate_size'];

        var errors = []; // already processed in validate
        var controlWellNamedRanges = self.getAllAssayControlWellRanges(values, errors);

        self.showAllAssayControlsWellSelector(assayPlateSize,controlWellNamedRanges)
      
      });

      editForm.$el.find('div[data-fields="assay_other_controls"]')
        .find('.form-group').append(show_all_assay_controls_button);
      
      // Show combined assay and library controls

      var show_all_assay_and_library_controls_button = $([
        '<a class="btn btn-default btn-sm" ',
          ' title="Show a combined view of assay and library controls in the specified library plate size" ',
          'role="button" id="show_all_assay_and_library_controls_button" >',
          'Show All Assay and Library Controls</a>'
        ].join(''));
      if (self.show_library_controls === false){
        show_all_assay_and_library_controls_button = $([
          '<a class="btn btn-default btn-sm" ',
          ' title="SHow a all of the assay controls in the specified library plate size" ',
            'role="button" id="show_all_assay_and_library_controls_button" >',
            'Show All Controls in Library Plate Format</a>'
          ].join(''));
      }      
      show_all_assay_and_library_controls_button.click(function(e){
        e.preventDefault();
        editForm.commit({ validate: true }); 
        var values = editForm.getValue();
        var assayPlateSize = values['assay_plate_size'];
        var libraryPlateSize = values['library_plate_size'];
        var errors = []; // already processed in validate
        var controlWellNamedRanges = self.getAllAssayControlWellRanges(values, errors);
        var libraryControlWells = [];
        
        if (self.show_library_controls === true){
          var libraryControlNamedRanges = Iccbl.parseNamedWellRanges(
            values['library_controls'],libraryPlateSize,errors);
          _.each(_.values(libraryControlNamedRanges), function(namedRange){
            libraryControlWells = libraryControlWells.concat(namedRange['wells']);
          });
        }
        console.log('show_all_assay_and_library_controls', 
          assayPlateSize, libraryPlateSize, errors);
        
        self.showAllControlsWellSelector(assayPlateSize, libraryPlateSize, 
          controlWellNamedRanges, libraryControlWells);
      }); // show_all_assay_and_library_controls_button button handler

      if (self.show_library_controls === true){
        editForm.$el.find('div[data-fields="library_controls"]')
          .find('.form-group').append(show_all_assay_and_library_controls_button);
      }else{
        editForm.$el.find('div[data-fields="assay_other_controls"]')
          .find('.form-group').append(show_all_assay_and_library_controls_button);
      }
    },
    
    showAllAssayControlsWellSelector: function(assayPlateSize, controlWellNamedRanges){
      var self = this;
      var errors = [];
      var duplicates = 
        Iccbl.find_duplicates(_.map(_.values(controlWellNamedRanges),
          function(nwr){
            return nwr['wells'];
          }
        ));
      if (!_.isEmpty(duplicates)){
        errors.push('wells are shared between ranges: ' + duplicates.join(', '));
      }
        
      var wellSelector = new WellSelector({
        plateSize: assayPlateSize,
        namedWellRangeInput: controlWellNamedRanges,
        colorPalette: self.CONTROL_WELLS_COLOR_PALETTE,
        useNamedRanges: true,
        errors: errors,
        editable: false
      });
      var el = wellSelector.render().el;

      var modalDialog = appModel.showModal({
        css: { 
            display: 'table',
            width: 'auto'
          },
        css_modal_content: {
          overflow: 'hidden'
        },
        okText: 'ok',
        view: el,
        title: 'All Assay Control Wells',
        ok: function(e) {
          e.preventDefault();
          if (!_.isEmpty(wellSelector.getErrors())){
            appModel.showJsonMessages({
              'Well Selector Errors': wellSelector.getErrors()
            });
          } else {
            // nop
          }
          wellSelector.remove();
        }              
      });
    },
    
    showAllControlsWellSelector: function(
      assayPlateSize, libraryPlateSize, assayControlNamedRanges, libraryControlWells){
      var self = this;
      var libraryControlWellRange;
      if (!_.isEmpty(libraryControlWells)){
        libraryControlWellRange = {
          wells: libraryControlWells,
          label: self.resource.fields['library_controls'].title
        };
      }
      
      var colorPalette = self.CONTROL_WELLS_COLOR_PALETTE;
      var quadLabels = ['QA1','QA2','QB1','QB2'];
      var controlWellNamedRanges = assayControlNamedRanges;
      
      if (assayPlateSize < libraryPlateSize){
        console.log('convolute assay wells from 96->384...');
        var convolutedControlRanges = {};
        var newOrdinal = 1;
        _.each(_.values(controlWellNamedRanges), function(namedWellRange){
          
          var wells_by_quadrant = [[],[],[],[]];
          // preserve the quadrant that the well came from
          _.each(namedWellRange['wells'], function(wellName){
            var convolutedWells = Iccbl.convoluteWell(assayPlateSize, libraryPlateSize, wellName);
            for(var quadrant=0;quadrant<4;quadrant++){
              wells_by_quadrant[quadrant].push(convolutedWells[quadrant]);
            }
          });
          var baseLabel = namedWellRange['label'];
          for(var quadrant=0;quadrant<4;quadrant++){
            var newLabel = baseLabel + '-' + quadLabels[quadrant];
            convolutedControlRanges[newLabel] = {
              label: newLabel,
              ordinal: newOrdinal + quadrant,
              wells: wells_by_quadrant[quadrant]
            };
          }
          newOrdinal += 4;
        });
        // Each control type is split into it's 4 quadrants:
        // Color palette, create 4 colors for each control type, using the
        // base color for the original:
        // positive_controls -red, negative -blue, other -green, library -yellow
        // see: http://google.github.io/palette.js/ for palettes chosen
        var positivePalette = palette('cb-YlOrRd',8).slice(3,7);
        var negativePalette = palette('cb-PuBu',8).slice(2,6);
        var otherPalette = palette('cb-RdYlGn',8).slice(4);
        var newColorPalette = [];
        newColorPalette = newColorPalette.concat(positivePalette);
        newColorPalette = newColorPalette.concat(negativePalette);
        newColorPalette = newColorPalette.concat(otherPalette)
        // add the library controls color (3)
        newColorPalette.push(colorPalette[3]);
        colorPalette = newColorPalette;
        console.log('new color palette', colorPalette);
        controlWellNamedRanges = convolutedControlRanges;
      }
      
      if (assayPlateSize <= libraryPlateSize){
        console.log('assayPlateSize <= libraryPlateSize...');
        if (libraryControlWellRange){
          controlWellNamedRanges['library_controls'] = _.extend({},
            libraryControlWellRange, { ordinal: _.size(controlWellNamedRanges)+1 });
        }
        var duplicates = 
          Iccbl.find_duplicates(_.map(_.values(controlWellNamedRanges),
            function(nwr){
              return nwr['wells'];
            }
          ));
        var errors = [];
        if (!_.isEmpty(duplicates)){
          errors.push('wells are shared between ranges: ' + duplicates.join(', '));
        }
        
        var wellSelector = new WellSelector({
          plateSize: libraryPlateSize,
          namedWellRangeInput: controlWellNamedRanges,
          colorPalette: colorPalette,
          useNamedRanges: true,
          errors: errors,
          editable: false
        });
        var el = wellSelector.render().el;

        var modalDialog = appModel.showModal({
          css: { display: 'table', width: 'auto' },
          css_modal_content: { overflow: 'hidden'},
          okText: 'ok',
          view: el,
          title: 'All Control Wells (Library plate format)',
          ok: function(e) {
            e.preventDefault();
            if (!_.isEmpty(wellSelector.getErrors())){
              appModel.showJsonMessages({
                'Well Selector Errors': wellSelector.getErrors()
              });
            } else {
              // nop
            }
            wellSelector.remove();
          }              
        });
      } else { // assayPlateSize > libraryPlateSize - show 4 quadrant panels
        console.log('assayPlateSize > libraryPlateSize - show 4 quadrant panels');
        
        var controlWellNamedRangesByQuadrant = [{},{},{},{}];
        
        _.each(_.values(controlWellNamedRanges), function(namedWellRange){
          var wellsByQuadrant = Iccbl.deconvoluteWells(
            assayPlateSize, libraryPlateSize,namedWellRange['wells']);
          console.log('wellsByQuadrant', namedWellRange, wellsByQuadrant)
          for (var quadrant=0;quadrant<4;quadrant++ ){
            controlWellNamedRangesByQuadrant[quadrant][namedWellRange['label']] = 
              _.extend({}, namedWellRange, {
                wells: wellsByQuadrant[quadrant]
              });
            if (libraryControlWellRange){
              controlWellNamedRangesByQuadrant[quadrant][libraryControlWellRange.label] =
                _.extend({},libraryControlWellRange, { ordinal: 4 });
            }
          }
        });
        console.log('controlWellNamedRangesByQuadrant', controlWellNamedRangesByQuadrant);  
        var $template = $([
          '<div id="wellselector_wrapper" >',
          '<div id="button_area"></div>',
          '<div id="grid_area"></div>',
          '</div>'].join(''));
        
        var $nextGridButton = $([
          '<button type="button" class="btn btn-default btn-xs" ',
          'id="next_quadrant_button" >>> show next quadrant</button>',
        ].join(''));
        $template.find('#button_area').append($nextGridButton);
        
        var title = '';
        if (libraryControlWellRange){
          title += 'Library Controls and ';
        }
        var currentQuadrant = 0;
        function showNextGrid(){
          var currentControlWellNamedRanges = controlWellNamedRangesByQuadrant[currentQuadrant];
          var duplicates = 
            Iccbl.find_duplicates(_.map(_.values(currentControlWellNamedRanges),
              function(nwr){
                return nwr['wells'];
              }
            ));
          var errors = [];
          if (!_.isEmpty(duplicates)){
            errors.push('wells are shared between ranges: ' + duplicates.join(', '));
          }
          console.log('showNextGrid', currentQuadrant, currentControlWellNamedRanges);
          var wellSelector = new WellSelector({
            plateSize: libraryPlateSize,
            namedWellRangeInput: currentControlWellNamedRanges,
            colorPalette: colorPalette,
            useNamedRanges: true,
            errors: errors,
            editable: false
          });
          var el = wellSelector.render().el;
          $('#modal_title').html(
            title + 'Assay Controls (Quadrant: '+quadLabels[currentQuadrant] + ')');
          console.log('set el', el);
          $template.find('#grid_area').empty().append(el);
          currentQuadrant +=1;
          if (currentQuadrant > 3) currentQuadrant = 0;
        }; // shwoNextGrid
        
        $nextGridButton.click(showNextGrid);
        showNextGrid();
        var modalDialog = appModel.showModal({
          css: { display: 'table', width: 'auto' },
          css_modal_content: { overflow: 'hidden'},
          okText: 'ok',
          view: $template,
          title: title + 'Assay Controls (Quadrant: '+quadLabels[0] + ')',
          ok: function(e) {
            // TODO: cleanup
          }              
        });
      }
    },

    showWellSelector: function(currentValue, plateSize, title, editable, callBack){
      console.log('showWellSelector...');
      var self = this;
      
      var errors = [];
      var namedWellRangeInput = Iccbl.parseNamedWellRanges(currentValue, plateSize, errors);
      var wellSelector = new WellSelector({
        plateSize: plateSize,
        useNamedRanges: true,
        namedWellRangeInput: namedWellRangeInput,
        errors: errors,
        editable: editable
      });
      var el = wellSelector.render().el;

      var modalDialog = appModel.showModal({
        css: { 
            display: 'table',
            width: 'auto'
          },
        css_modal_content: {
          overflow: 'hidden'
        },
        okText: 'ok',
        view: el,
        title: title,
        ok: function(e) {
          e.preventDefault();
          if (!_.isEmpty(wellSelector.getErrors())){
            appModel.showJsonMessages({
              'Well Selector Errors': wellSelector.getErrors()
            });
          } else {
            var namedWellRanges = wellSelector.getSelectedWells();
            var wellRangeString = Iccbl.generateNamedWellBlockString(
              namedWellRanges,plateSize);
            callBack(wellRangeString);
          }
          wellSelector.remove();
        }              
      });
    }, // end showWellsToLeaveEmptyDialog
    
    _create_input_file_form: function($el, input_file_model) {
      
      var self = this;
      console.log('create_input_file form', $el, input_file_model);

      var fields = self.resource.input_file_fields;
      console.log('input_file_fields', self.resource.input_file_fields);
      var ordinal = input_file_model.get('ordinal');
      var file_input_txt = 
        '<label id="file-button" class="btn btn-default btn-file">' + 
          'Browse<input type="file" name="file_input" style="display: none;"></label>'; 
      if (input_file_model.has('filename')) {
        file_input_txt += 
          '<p id="filechosen" class="form-control-static" >&lt;' + 
          input_file_model.get('filename') + '&gt;</p>';
        input_file_model.unset('filename');
      } else {
        file_input_txt += '<p id="filechosen" class="form-control-static" ></p>';
      }
      var form_template = [
         "<div class='form-horizontal panel panel-info' >",
         " <div class='panel-heading'>",
         "  <h4 class='panel-title' >File: {display_ordinal}",
         '  <button type="button" class="btn btn-default btn-xs controls-right" ',
         '    id="remove_{ordinal}" >remove</button></h4>',
         " </div>",
         " <div class='panel-body'>",
         "  <div class='row'>",
         "   <div id='general_errors' class='alert alert-danger' style='display: none;' ></div>",
         "  </div>",
         "  <form id='input_file_form_{ordinal}' class='form-horizontal col-sm-10' >",
         "   <div data-fieldsets ></div>",
         "  </form>",
         " </div>",
         "</div>"].join('');
      form_template = Iccbl.formatString(
        form_template, 
        { ordinal: ordinal,
          display_ordinal: ordinal+1
        })
      var fieldTemplate = appModel._horizontal_form_2col_field_template;
      var DisabledField = EditView.DisabledField.extend({
        tagName: 'div',
        className: 'form-control-static'
      });

      function CheckFileExistsValidator(value, formValues) {
        var file = $('input[name="file_input"]')[0].files[0];
        console.log('found file:', file);
        if (!file){
          return {
            type: 'file',
            message: 'Input file is required'
          };
        }
      };

      var formSchema = {};
      formSchema['filename'] = {
        title: 'File',
        key: 'filename',
        type: EditView.DisabledField.extend({
          tagName: 'div',
          className: ''
        }), 
        validators: [CheckFileExistsValidator],
        template: fieldTemplate
      };
      formSchema['collation_order'] = {
        title: fields['collation_order'].title,
        help: fields['collation_order'].description,
        key: 'collation_order',
        type: Backbone.Form.editors.Select.extend({
            className: 'form-control'
          }),
        options: Iccbl.appModel.getVocabularySelectOptions(
          fields['collation_order']['vocabulary_scope_ref']),
        template: fieldTemplate,
        validators: ['required']
      };
      formSchema['readout_type'] = {
        title: 'Readout Type',
        key: 'readout_type',
        type: Backbone.Form.editors.Select.extend({
            className: 'form-control'
          }),
        options: Iccbl.appModel.getVocabularySelectOptions(
          fields['readout_type']['vocabulary_scope_ref']),
        template: fieldTemplate,
        validators: ['required']
      };
      formSchema['conditions'] = {
        title: 'Conditions',
        key: 'conditions',
        type: EditView.TextArea2,
        editorClass: 'form-control input-full',
        template: fieldTemplate
      };
      formSchema['replicates'] = {
        title: '# Replicates',
        key: 'replicates',
        validators: [EditView.CheckPositiveNonZeroValidator],
        type: Backbone.Form.editors.Number,
        editorClass: 'form-control',
        template: fieldTemplate
      };
      formSchema['readouts'] = {
        title: 'Readout Names',
        key: 'readouts',
        type: EditView.TextArea2,
        editorClass: 'form-control input-full',
        template: fieldTemplate
      };
      
      input_file_model.schema = formSchema;
      
      var form = new Backbone.Form({
        model: input_file_model,
        template: _.template(form_template),
      });
      input_file_model.form = form;
      var _form_el = form.render().el;

      form.$el.find("[name='filename']").append(file_input_txt);
      form.$el.on('change', ':file', function() {
        var input = $(this),
            numFiles = input.get(0).files ? input.get(0).files.length : 1,
            label = input.val().replace(/\\/g, '/').replace(/.*\//, '');
        form.$el.find('#filechosen').html(label);
        form.setValue('filename', label);
      });
      if (ordinal == 0 ){
        form.$el.find('#remove_'+ordinal).hide();
      }
      form.$el.find('#remove_'+ordinal).click(function(){
        if (ordinal == 0){
          console.log('cannot remove the only file');
        } else {
          form.$el.remove();
          self.inputFileCollection.remove(input_file_model);
        }
      });
      $el.append(_form_el);
    },
    
    /** Backbone.layoutmanager callback **/
    cleanup: function(){
      console.log('cleanup called...');
      this.model = null;
      this.screen = null;
      this.cherry_pick_request = null;
      this.args = null;
    },
    
    /**
     * Child view bubble up URI stack change event
     * TODO: not used here
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    }
  
  });

  return RawDataTransformView;
});