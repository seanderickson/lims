define([
  'jquery',
  'underscore',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_edit'  
], 
function($, _, Backgrid, Iccbl, appModel, EditView) {

  var PlateRangeTablePrototype = {

    /**
     * Parse a single entry in the list of libray_copy_plate_range
     */
    _parse_library_plate_entry: function(library_copy_plate_range_regex, entry){
      var self = this;
      // E.G. Regex: /(([^:]+):)?(\w+):(\d+)(-(\d+))?/
      if(library_copy_plate_range_regex){
        var match = library_copy_plate_range_regex.exec(entry);
        if(match !== null){
          temp = {
            library: match[2], // optional
            copy: match[3],
            start_plate: match[4],
            end_plate: match[6]
          };
          return temp;
        }else{
          appModel.error('entry does not match pattern: ' 
            + entry + ', ' + library_copy_plate_range_regex );
        }
      }else{
        return entry;
      }
    },
    
    _createPlateRangeTable: function(plate_collection, $target_el, nested_library_plate_pattern, editable){
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
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library}/copy/{copy}/plate/{start_plate}'
          })
        }),
        _.extend({},colTemplate,{
          'name' : 'end_plate',
          'label' : 'End Plate',
          'description' : 'End Plate',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library}/copy/{copy}/plate/{end_plate}'
          })
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
        // Removed - see libraryScreeningView, already handled there
        //self.listenTo(plate_collection, "MyCollection:delete", function (model) {
        //  var entry = Iccbl.formatString(nested_library_plate_pattern,model);
        //  plate_collection.remove(model);
        //});
        
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

    _addPlateRangeDialog: function(
        plate_collection, libraryOptions, nested_library_plate_pattern, for_screening){
      
      var self = this;
      var options = options || {};
      var description = 'Enter a plate range';
      var formSchema = {};
      var fieldTemplate = appModel._field_template;
      var formTemplate = appModel._form_template;
      var initfun = function(){
        console.log('initfun...');
//        libraryOptions = appModel.getScreeningLibraryOptions(self.model.get('screen_type'));
//        libraryOptions.unshift({val:'',label:''});
        
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
          if (for_screening===true){
            _.each(libraryRecord.get('screening_copies'),function(copy){
              $chosen.append($('<option>',{
                  value: copy
                }).text(copy));
            });
          }else{
            _.each(libraryRecord.get('copies'),function(copy){
              $chosen.append($('<option>',{
                  value: copy
                }).text(copy));
            });
          }
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
    } //_addPlateRangeDialog
  };
  
  return PlateRangeTablePrototype;
});



