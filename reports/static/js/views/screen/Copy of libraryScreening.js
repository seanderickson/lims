define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_detail_stickit',
  'views/generic_edit',
  'views/generic_detail_layout'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, 
         DetailView, EditView, DetailLayout) {
  
  var nested_library_plate_pattern = '{library}:{copy}:{start_plate}-{end_plate}';
  
  var parse_library_plate_entry = function(entry){
    var vals = entry.split(':');
    return {
      library: vals[0],
      copy: vals[1],
      start_plate: vals[2].split('-')[0],
      end_plate: vals[2].split('-')[1]
    };
  };
  var LibraryScreeningView = DetailLayout.extend({

    _createPlateRangeTable: function($target_el, editable){
      var self = this;
      var collection = new Backbone.Collection();
      _.each(self.model.get('library_plates_screened'),function(entry){
        collection.add(new Backbone.Model(parse_library_plate_entry(entry)));
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
            'name' : 'library',
            'label' : 'Library',
            'description' : 'Library Short Name',
            'order': 1,
            'sortable': true,
            'cell': TextWrapCell
          }),
          _.extend({},colTemplate,{
            'name' : 'copy',
            'label' : 'Copy',
            'description' : 'Copy Name',
            'order': 1,
            'sortable': true,
            'cell': TextWrapCell
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
        self.listenTo(collection, "MyCollection:delete", function (model) {
          var entry = Iccbl.replaceTokens(model,nested_library_plate_pattern);
          console.log('delete: ' , model, entry);
          var lps = self.model.get('library_plates_screened');
          lps = _.without(lps,_.find(lps,function(val){
            return val === entry;
          }));
          console.log('new lps', lps);
          collection.remove(model);
          self.model.set('library_plates_screened',lps);
        });
        
      }
      var colModel = new Backgrid.Columns(columns);
      colModel.comparator = 'order';
      colModel.sort();

      $target_el.empty();
      var cell = $('<div>',{ class: 'col-sm-4' });
      
      var plate_range_grid = new Backgrid.Grid({
        columns: colModel,
        collection: collection,
        className: 'backgrid table-striped table-condensed table-hover'
      });
      cell.html(plate_range_grid.render().$el);
      $target_el.append(cell);
      
    }, // _createPlateRangeTable

    _addPlateRangeDialog2: function(callback){
      
    },
    
    _addPlateRangeDialog: function(callback){
      
      var self = this;
      var options = options || {};
      var description = 'Enter a plate range';
      var formSchema = {};
      var fieldTemplate = appModel._field_template;
      var formTemplate = appModel._form_template;
      var initfun = function(){
        console.log('initfun...');
        libraryOptions = appModel.getLibraryOptions();
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
            var library_plates = self.model.get('library_plates_screened');
            _.each(library_plates, function(entry){
              if (newEntry==entry){
                errs['library'] = 'Entry exists: ' + newEntry;
              }
              var range = parse_library_plate_entry(entry);
              
              if(range['library']==attrs['library']){
                if(range['copy']==attrs['copy']){
                  var start1 = parseInt(range['start_plate']);
                  var end1 = parseInt(range['end_plate']);
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
          _.each(libraryRecord.get('copies'),function(copy){
            $chosen.append($('<option>',{
                value: copy
              }).text(copy));
          });
          $chosen.trigger("chosen:updated");
        });
        form.on("start_plate:change", function(e){
          
        });
          
        var dialog = appModel.showModal({
          okText: 'Create',
          view: _form_el,
          title: 'Create a new plate range',
          ok: function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }) || {}; // runs schema and model validation
            if(_.isEmpty(errors) ){
              var library_plates = self.model.get('library_plates_screened');
              library_plates.push(newEntry);
              self.model.set('library_plates_screened',library_plates);
              $target_el = self.$el.find('[key="library_plates_screened"]');
              self._createPlateRangeTable($target_el,true);
              // TODO: mark the model/page as needing to be saved
              // implement on back end
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
          $target_el = this.$el.find('#library_plates_screened');
          self._createPlateRangeTable($target_el);
        }
      });
      
      var editView = EditView.extend({
        afterRender: function(){
          EditView.prototype.afterRender.apply(this,arguments);

          $target_el = this.$el.find('[key="library_plates_screened"]');
          self._createPlateRangeTable($target_el,true);

          // create a button and a form to add library copy plate ranges to the model
          var addButton = $([
            '<button class="btn btn-default btn-sm" ',
              'role="button" id="add_library_plates_screened_button" href="#">',
              'Add</button>'
            ].join(''));
          addButton.click(function(event){
            event.preventDefault();
            self._addPlateRangeDialog();
          });
          $target_el = $target_el.parent();
          $target_el.append(addButton);
        }
      });
      args.EditView = editView;
      args.DetailView = detailView;
      DetailLayout.prototype.initialize.call(this,args);
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
        var fields = self.model.resource.schema.fields;
        fields['performed_by_username']['choices'] = (
            [{ val: '', label: ''}].concat(userOptions));
        DetailLayout.prototype.showEdit.apply(self,arguments);
      });  
    }
    
  });

  return LibraryScreeningView;
});