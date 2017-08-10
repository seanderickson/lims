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
            library_short_name: match[2], // optional
            copy_name: match[3],
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
    
    _createPlateRangeTable: function(
        plate_collection, $target_el, editable, extra_cols, screen_facility_id){
      var self = this;
      $target_el.empty();
      
      var colTemplate = {
        'cell' : 'string',
        'order' : -1,
        'sortable': false,
        'searchable': false,
        'editable' : false,
        'visible': true,
        'headerCell': Backgrid.HeaderCell
      };
      var columns = [];
      
      var libraryColumn = _.extend({},colTemplate,{
          'name' : 'library_short_name',
          'label' : 'Library',
          'description' : 'Library Short Name',
          'order': 1,
          'sortable': true,
          'cell':
            Iccbl.CommentArrayLinkCell.extend({
              hrefTemplate: '#library/{library_short_name}',
              comment_attribute: 'library_comment_array',
              title_function: function(model){
                return 'Comments for library: ' + model.get('library_short_name');
              }
            })
          });
      columns.push(libraryColumn);
      
      // TODO: rework copy cell when copy comments are implemented in apilogs
      columns.push(          
        _.extend({},colTemplate,{
          'name' : 'copy_name',
          'label' : 'Copy',
          'description' : 'Copy Name',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library_short_name}/copy/{copy_name}',
            render: function(){
              var self = this;
              Iccbl.LinkCell.prototype.render.apply(this, arguments);
              var comments = this.model.get('copy_comments');
              if (!_.isEmpty(comments)){
                this.$el.attr('title', comments);
                this.$el.append(Iccbl.createCommentIcon(
                  comments,
                  'Comments for Copy: ' 
                    + self.model.get('library_short_name') + '/'
                    + self.model.get('source_copy_name')
                  ));
              }
              return this;
            }
          })
        }));
      columns.push(          
        _.extend({},colTemplate,{
          'name' : 'plate_range',
          'label' : 'Plate Range',
          'description' : 'Copy Plate range',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': 
              '#library/{library_short_name}/copy/{copy_name}/plate' +
              '/search/plate_number__range={start_plate},{end_plate}',
            render : function() {
              var self = this;
              this.$el.empty();
              var start_plate = self.model.get('start_plate');
              var end_plate = self.model.get('end_plate');
              var formattedValue = start_plate;
              if (start_plate != end_plate){
                formattedValue += '-' + end_plate;
              }
              var interpolatedVal = Iccbl.formatString(self.hrefTemplate,self.model);
              var $link = $('<a>', {
                  tabIndex : -1,
                  href : interpolatedVal,
                  target : self.target,
                  title: self.title
                }).text(formattedValue);

              var comments = self.model.get('plate_comment_array');
              if (!_.isEmpty(comments)){
                comments = Iccbl.parseComments(comments);
                $link.attr('title',comments);
                $link.append(Iccbl.createCommentIcon(
                  comments,
                  'Comments for ' + formattedValue ));
              }
              self.$el.append($link);
              
              return this;
            },
          })
        }));
      
      if(_.contains(extra_cols, '-library_screening_id')){
        // omit this column
      }
      else if(_.contains(extra_cols,'library_screening_id') ||
          (!plate_collection.isEmpty() 
              && plate_collection.at(0).has('library_screening_id'))){
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'library_screening_id',
            'label' : 'Screening ID',
            'description' : 'Library Screening ID',
            'order': -1,
            'sortable': true,
            'formatter': _.extend({}, Iccbl.StringFormatter.prototype, {
              fromRaw(rawValue){
                var val = Iccbl.StringFormatter.prototype.fromRaw.call(this,rawValue);
                if (val=='0') return '';
                return val;
              }
            }),
            'cell': Iccbl.LinkCell.extend({
              'hrefTemplate': '#screen/'+ screen_facility_id 
              +'/summmary/libraryscreening/{library_screening_id}'
            })
          })
        );
      }
      if(_.contains(extra_cols, '-plate_locations')){
        // omit this column
      }
      else if(_.contains(extra_cols,'plate_locations') ||
          (!plate_collection.isEmpty() 
              && plate_collection.at(0).has('plate_locations'))){
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'plate_locations',
            'label' : 'Locations',
            'description' : 'Plate Locations',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.TextWrapCell
          })
        );
      }      
      
      if(_.contains(extra_cols,'library_screening_status') ||
          (!plate_collection.isEmpty() 
              && plate_collection.at(0).has('library_screening_status'))){
        var optionValues = [];
        var vocabulary_scope_ref = 'library.screening_status';
        try{
          var vocabulary = Iccbl.appModel.getVocabulary(vocabulary_scope_ref);
            _.each(_.keys(vocabulary),function(choice){
              optionValues.push([vocabulary[choice].title,choice]);
            });
        }catch(e){
          console.log('build column errorr: e',e);
        }
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'library_screening_status',
            'label' : 'Library Screening Status',
            'description' : 'Library Screening Status',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.SelectCell.extend({
              optionValues: optionValues,
              vocabulary_scope_ref: vocabulary_scope_ref
            })
          })
        );
      }
      var tableClasses = 'backgrid table-striped table-condensed table-hover'
      var rowClass = null;
      var StatusColorRow = Backgrid.Row.extend({
        _setStyle: function() {
          
          if (!_.isEmpty(this.model.get('errors'))) {
            this.$el.addClass('danger');
          }
          else if (!_.isEmpty(this.model.get('warnings'))) {
            this.$el.addClass('warning');
          }
        },
        render: function() {
          StatusColorRow.__super__.render.apply(this, arguments);
          this._setStyle();
          return this;
        }
      });
      
      if(_.contains(extra_cols,'warnings') ||
          (!plate_collection.isEmpty() && plate_collection.at(0).has('warnings'))){
        rowClass = StatusColorRow;
        tableClasses = 'backgrid table-condensed table';
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'warnings',
            'label' : 'Warnings',
            'description' : 'Warnings',
            'order': 1,
            'sortable': true,
            'formatter': _.extend({}, Iccbl.StringFormatter.prototype, {
              fromRaw(rawValue){
                if (!_.isEmpty(rawValue)){
                  return rawValue.join('<br>');
                }
                return '';
              }
            }),
            'cell': Iccbl.TextWrapCell.extend({
              className: 'text-wrap-cell-narrow'
            })
          })
        );
      }
      
      if(_.contains(extra_cols,'errors') ||
          (!plate_collection.isEmpty() && plate_collection.at(0).has('errors'))){
        rowClass = StatusColorRow;
        tableClasses = 'backgrid table-condensed table';
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'errors',
            'label' : 'Errors',
            'description' : 'Errors: screening not possible',
            'order': 1,
            'sortable': true,
            'formatter': _.extend({}, Iccbl.StringFormatter.prototype, {
              fromRaw(rawValue){
                if (!_.isEmpty(rawValue)){
                  return rawValue.join('<br>');
                }
                return '';
              }
            }),
            'cell': Iccbl.TextWrapCell.extend({
              className: 'text-wrap-cell-narrow'
            })
          })
        );
      }
      
      if(editable && editable === true ){
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'delete',
            'label' : 'Remove',
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

      var cell = $('<div>',{ class: '' });
      
      var plate_range_grid = new Backgrid.Grid({
        row: rowClass,
        columns: colModel,
        collection: plate_collection,
        className: tableClasses
      });
//      this.listenTo(plate_collection,'backgrid:sort', function(col,direction){
//        if (col.get('name')=='plate_range'){
//          // TODO
//        }
//      });

      cell.append(plate_range_grid.render().$el);
      $target_el.append(cell);
      
    }, // _createPlateRangeTable


    // TODO: 20170412 - rework this to use the plate search functionality:
    // see libraryscreening.js
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
    } //_addPlateRangeDialogOld
  };
  
  return PlateRangeTablePrototype;
});



