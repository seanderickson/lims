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
  
  var nested_library_plate_pattern = '{library}:{copy}:{start_plate}-{end_plate}';
  
  var LibraryScreeningView = DetailLayout.extend({

    /**
     * Parse a single entry in the list of libray_copy_plate_range
     */
    _parse_library_plate_entry: function(entry){
      var self = this;
      if(self.library_plates_screened_regex){
        var match = self.library_plates_screened_regex.exec(entry);
        if(match !== null){
          return {
            library: match[1],
            copy: match[2],
            start_plate: match[3],
            end_plate: match[4]
          };
        }else{
          appModel.error('library_plates_screened entry does not match pattern: ' + regex );
        }
      }else{
        return entry;
      }
    },
    
    _createPlateRangeTable: function(plate_collection, $target_el, editable){
      var self = this;
      $target_el.empty();
      
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
          'name' : 'library',
          'label' : 'Library',
          'description' : 'Library Short Name',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library}'
          })
        }),
        _.extend({},colTemplate,{
          'name' : 'copy',
          'label' : 'Copy',
          'description' : 'Copy Name',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library}/copy/{copy}'
          })
        }),
        _.extend({},colTemplate,{
          'name' : 'start_plate',
          'label' : 'Start Plate',
          'description' : 'Start Plate',
          'order': 1,
          'sortable': true,
          'cell': TextWrapCell
        }),
        _.extend({},colTemplate,{
          'name' : 'end_plate',
          'label' : 'End Plate',
          'description' : 'End Plate',
          'order': 1,
          'sortable': true,
          'cell': TextWrapCell
        })
      ];
      
      if(editable && editable === true ){
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'delete',
            'label' : '',
            'description' : 'delete',
            'text': 'X',
            'order': 1,
            'sortable': false,
            'cell': Iccbl.DeleteCell
          }));
        self.listenTo(plate_collection, "MyCollection:delete", function (model) {
          var entry = Iccbl.formatString(nested_library_plate_pattern,model);
          plate_collection.remove(model);
        });
      }
      var colModel = new Backgrid.Columns(columns);
      colModel.comparator = 'order';
      colModel.sort();

      var cell = $('<div>',{ class: 'col-sm-10' });
      
      var plate_range_grid = new Backgrid.Grid({
        columns: colModel,
        collection: plate_collection,
        className: 'backgrid table-striped table-condensed table-hover'
      });
      cell.append(plate_range_grid.render().$el);
      $target_el.append(cell);
      
    }, // _createPlateRangeTable

    _addPlateRangeDialog: function(plate_collection){
      
      var self = this;
      var options = options || {};
      var description = 'Enter a plate range';
      var formSchema = {};
      var fieldTemplate = appModel._field_template;
      var formTemplate = appModel._form_template;
      var initfun = function(){
        console.log('initfun...');
        libraryOptions = appModel.getScreeningLibraryOptions(self.model.get('screen_type'));
        libraryOptions.unshift({val:'',label:''});
        
        formSchema['library'] = {
          title: 'Library',
          key: 'library',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          options: libraryOptions,
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['start_plate'] = {
          title: 'Start Plate',
          key: 'start_plate',
          validators: ['required'],
          type: EditView.ChosenSelect,
          options: [{val:'',label:''}],
          editorClass: 'chosen-select',
          template: fieldTemplate
        };
        formSchema['end_plate'] = {
          title: 'End Plate',
          key: 'end_plate',
          validators: ['required'],
          type: EditView.ChosenSelect,
          options: [{val:'',label:''}],
          editorClass: 'chosen-select',
          template: fieldTemplate
        };
        formSchema['copy'] = {
          title: 'Copy',
          key: 'copy',
          validators: ['required'],
          type: EditView.ChosenSelect,
          options: [{val:'',label:''}],
          editorClass: 'chosen-select',
          template: fieldTemplate
        };
  
        var FormFields = Backbone.Model.extend({
          schema: formSchema,
          validate: function(attrs){
            var errs = {};
            var newEntry = Iccbl.formatString(nested_library_plate_pattern, attrs );
            var start2 = parseInt(attrs['start_plate']);
            var end2 = parseInt(attrs['end_plate']);
            if (start2>end2){
              var temp = end2;
              end2 = start2;
              attrs['end_plate'] = end2;
              start2 = temp;
              attrs['start_plate'] = temp;
            }
            plate_collection.each(function(model){
              var plate_string = Iccbl.formatString(nested_library_plate_pattern, model);
              if (newEntry==plate_string){
                errs['library'] = 'Entry exists: ' + newEntry;
              }
              if(model.get('library')==attrs['library']){
                if(model.get('copy')==attrs['copy']){
                  var start1 = parseInt(model.get('start_plate'));
                  var end1 = parseInt(model.get('end_plate'));
                  if(start2<end1 && end2 > start1 ){
                    errs['start_plate'] = 'overlapping range';
                  }
                }
              }
              
            });
            if (!_.isEmpty(errs)) return errs;
          }
        });
        var formFields = new FormFields();
        
        var form = new Backbone.Form({
          model: formFields,
          template: formTemplate
        });
        var formview = form.render();
        var _form_el = formview.el;
        formview.$el.find('.chosen-select').chosen({
          disable_search_threshold: 3,
          width: '100%',
          allow_single_deselect: true,
          search_contains: true
        });

        form.on("library:change", function(e){
          var library = form.getValue('library');
          console.log('change:library_name ' + library);
          var libraries = appModel.getLibraries();
          var libraryRecord = libraries.find(function(model){
            return model.get('short_name') == library;
          });
          var start=parseInt(libraryRecord.get('start_plate'));
          var end=parseInt(libraryRecord.get('end_plate'));
          var fieldKey = 'start_plate';
          var $chosen = form.$el.find('[name="' + fieldKey + '"]').parent()
              .find('.chosen-select');
          $chosen.empty();
          for(var i = start; i<=end; i++){
            $chosen.append($('<option>',{
                value: i
              }).text(i));
          }; 
          $chosen.trigger("chosen:updated");
          fieldKey = 'end_plate';
          $chosen = form.$el.find('[name="' + fieldKey + '"]').parent()
              .find('.chosen-select');
          $chosen.empty();
          for(var i = start; i<=end; i++){
            $chosen.append($('<option>',{
                value: i
              }).text(i));
          }; 
          $chosen.trigger("chosen:updated");

          fieldKey = 'copy';
          $chosen = form.$el.find('[name="' + fieldKey + '"]').parent()
              .find('.chosen-select');
          $chosen.empty();
          _.each(libraryRecord.get('screening_copies'),function(copy){
            $chosen.append($('<option>',{
                value: copy
              }).text(copy));
          });
          $chosen.trigger("chosen:updated");
        });
        form.on("start_plate:change", function(e){
          // TODO: set the end plate > start plate
        });
          
        var dialog = appModel.showModal({
          okText: 'Create',
          view: _form_el,
          title: 'Create a new plate range',
          ok: function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }) || {}; 
            if(_.isEmpty(errors) ){
              plate_collection.add(form.getValue());
              return true;
            }else{
              _.each(_.keys(errors), function(key){
                form.$el.find('[name="'+key +'"]').parents('.form-group').addClass('has-error');
              });
            }            
            return false;
          }
        });
        
      };
      $(this).queue([appModel.getLibraries,initfun]);
    }, //_addPlateRangeDialog
    
    initialize: function(args) {
      
      var self = this;
      
      var detailView = DetailView.extend({
        
        afterRender: function(){
          
          DetailView.prototype.afterRender.apply(this,arguments);

          var $target_el = this.$el.find('#library_plates_screened');
          var plate_collection = new Backbone.Collection;
          plate_collection.comparator = 'start_plate';
          plate_collection.set(
            _.map(self.model.get('library_plates_screened'),
                  self._parse_library_plate_entry));
          self._createPlateRangeTable(plate_collection, $target_el);
          
          if(!_.isEmpty(self.model.get('assay_protocol'))){
            $target_el = this.$el.find('#assay_protocol_last_modified_date');
            var apilogResource = appModel.getResource('apilog');
            var options = {
              data_for_get: {
                limit: 1,
                key: self.model.get('activity_id'),
                ref_resource_name: self.model.resource.key,
                diff_keys__icontains: '"assay_protocol"',
                order_by: ['-date_time']
              }
            };
            Iccbl.getCollectionOnClient(apilogResource.apiUri,function(collection){
              if(collection && !collection.isEmpty()){
                var model = collection.at(0);
                console.log('assay_protocol_last_modified_date',model.get('date_time'));
                self.model.set('assay_protocol_last_modified_date',model.get('date_time'));
              }else{
                // no apilog found, set it to the creation date
                self.model.set('assay_protocol_last_modified_date',
                  self.model.get('date_created'));
              }
            },options);
          }
        
        }
      });
      
      var editView = EditView.extend({
        
        afterRender: function(){

          EditView.prototype.afterRender.apply(this,arguments);

          var self_editform = this;
          var $target_el = this.$el.find('[name="library_plates_screened"]');
          var plate_collection = new Backbone.Collection;
          plate_collection.comparator = 'start_plate';
          self.listenTo(plate_collection, "MyCollection:delete", function (model) {
            var entry = Iccbl.formatString(nested_library_plate_pattern,model);
            plate_collection.remove(model);
          });

          self_editform.listenTo(plate_collection, 'update', function (plate_collection) {
            self_editform.setValue('library_plates_screened', plate_collection.map(function(model){
              return Iccbl.formatString(nested_library_plate_pattern, model);
            }));
          });
          var lps = self.model.get('library_plates_screened');
          if(_.isEmpty(lps)) lps = [];
          plate_collection.set(_.map(lps,self._parse_library_plate_entry));
          self._createPlateRangeTable(plate_collection, $target_el, true);

          // library_copy_plates button
          var addButton = $([
            '<button class="btn btn-default btn-sm" ',
              'role="button" id="add_library_plates_screened_button" href="#">',
              'Add</button>'
            ].join(''));
          addButton.click(function(event){
            event.preventDefault();
            self._addPlateRangeDialog(plate_collection);
          });
          $target_el = $target_el.parent();
          $target_el.append(addButton);
          
          // attach is_for_external_library_plates change listener
          $target_el = self_editform.$el.find('[key=form-group-library_plates_screened]');
          this.listenTo(this, "is_for_external_library_plates:change", function(e){
            var val = self_editform.getValue('is_for_external_library_plates');
            console.log('is_for_external_library_plates:',val)
            if(val){
              $target_el.hide();
              self_editform.setValue('library_plates_screened',null);
            } else {
              $target_el.show();
            }
          });
          // attach listener volume calculator
          function calculateVolFromLibraryPlates(){
            var num_replicates = self_editform.getValue('number_of_replicates');
            var vol_to_assayplates  = self_editform.getValue('volume_transferred_per_well_to_assay_plates');
            if(_.isNumber(num_replicates) && _.isNumber(vol_to_assayplates)){
              if(num_replicates > 0 && vol_to_assayplates > 0){
                var calculated = num_replicates * vol_to_assayplates;
                self_editform.setValue('volume_transferred_per_well_from_library_plates',calculated);
              }
            } else {
              console.log('cannot calculate volume until form values are entered');
            }
          };
          var calcVolButton = $([
            '<button class="btn btn-default btn-sm" ',
              'role="button" id="calc_volume_transferred_per_well_from_library_plates" href="#">',
              'calculate</button>'
            ].join(''));
          $target_el = self_editform.$el.find('[key=form-group-volume_transferred_per_well_from_library_plates]');
          $target_el.append(calcVolButton);
          calcVolButton.click(function(event){
            event.preventDefault();
            calculateVolFromLibraryPlates();
          });
          this.listenTo(this, "number_of_replicates:change", calculateVolFromLibraryPlates);
          this.listenTo(this, "volume_transferred_per_well_to_assay_plates:change", calculateVolFromLibraryPlates);
        }
      });
      args.EditView = editView;
      args.DetailView = detailView;
      
      DetailLayout.prototype.initialize.call(this,args);

      if (!_.has(self.modelFields,'library_plates_screened')
          || !_.has(self.modelFields['library_plates_screened'],'regex')){
        console.log('error....field schema missing/regex for library_plates_screened',
          self, self.modelFields);
        appModel.error('field schema missing/regex for library_plates_screened');
      }
      var regex_string = self.modelFields['library_plates_screened']['regex'];
      try{
        self.library_plates_screened_regex = RegExp(regex_string);
      }catch(e){
        appModel.error('regex misconfigured for "library_plates_screened" in metadata: ' + regex_string);
      }
      
      _.bindAll(this, '_parse_library_plate_entry');
    },
    
    afterRender: function() {
      DetailLayout.prototype.afterRender.apply(this,arguments);
      $title = this.$el.find('#content_title');
      $title.html('Library Screening');
      $title.parent().show();
    },    
    
    showEdit: function() {
      var self = this;
      appModel.initializeAdminMode(function(){
        var userOptions = appModel.getAdminUserOptions();
        var fields = self.model.resource.fields;
        fields['performed_by_username']['choices'] = (
            [{ val: '', label: ''}].concat(userOptions));
        DetailLayout.prototype.showEdit.apply(self,arguments);
      });  
    }
    
  });

  return LibraryScreeningView;
});