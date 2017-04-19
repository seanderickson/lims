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
        plate_collection, $target_el, editable, extra_cols){
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
      var columns = [
        _.extend({},colTemplate,{
          'name' : 'library_short_name',
          'label' : 'Library',
          'description' : 'Library Short Name',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library_short_name}'
          })
        }),
        _.extend({},colTemplate,{
          'name' : 'copy_name',
          'label' : 'Copy',
          'description' : 'Copy Name',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library_short_name}/copy/{copy_name}'
          })
        }),
        _.extend({},colTemplate,{
          'name' : 'start_plate',
          'label' : 'Start Plate',
          'description' : 'Start Plate',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library_short_name}/copy/{copy_name}/plate/{start_plate}'
          })
        }),
        _.extend({},colTemplate,{
          'name' : 'end_plate',
          'label' : 'End Plate',
          'description' : 'End Plate',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.LinkCell.extend({
            'hrefTemplate': '#library/{library_short_name}/copy/{copy_name}/plate/{end_plate}'
          })
        })
        
      ];
      
      if(_.contains(extra_cols,'library_screening_status') ||
          (!plate_collection.isEmpty() && plate_collection.at(0).has('library_screening_status'))){
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
            'label' : 'Screening Status',
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
      //if(_.contains(extra_cols,'plate_statuses') ||
      //    (!plate_collection.isEmpty() && plate_collection.at(0).has('plate_statuses'))){
      //  columns.push(          
      //    _.extend({},colTemplate,{
      //      'name' : 'plate_statuses',
      //      'label' : 'Plate Status',
      //      'description' : 'Plate Statuses',
      //      'order': 1,
      //      'sortable': true,
      //      'cell': Iccbl.StringCell.extend({
      //        'className': 'select-cell',
      //        formatter: _.extend({}, Iccbl.StringFormatter.prototype, {
      //          fromRaw(rawValue){
      //            if(!_.isEmpty(rawValue)&&_.isArray(rawValue)){
      //              return _.map(rawValue,function(status){
      //                return status.charAt(0).toUpperCase() + status.slice(1);
      //              }).join(', ');
      //            }
      //            return Iccbl.StringFormatter.prototype.fromRaw.call(this,rawValue);
      //          }
      //        })
      //      })
      //    })
      //  );
      //}
      if(_.contains(extra_cols,'warnings') ||
          (!plate_collection.isEmpty() && plate_collection.at(0).has('warnings'))){
        columns.push(          
          _.extend({},colTemplate,{
            'name' : 'warnings',
            'label' : 'Warnings',
            'description' : 'Warnings',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.TextWrapCell
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
        columns: colModel,
        collection: plate_collection,
        className: 'backgrid table-striped table-condensed table-hover'
      });
      cell.append(plate_range_grid.render().$el);
      $target_el.append(cell);
      
    }, // _createPlateRangeTable

//    _addPlateRangeDialog2: function(
//        plate_collection, libraryOptions, nested_library_plate_pattern, for_screening){
//      
//      // create the search form
//      var TextArea2 = Backbone.Form.editors.TextArea.extend({
//        render: function() {
//          TextArea2.__super__.render.apply(this,arguments);
//          this.$el.attr('placeholder', this.schema.placeholder);
//          return this;
//        },        
//      });
//      var formSchema = {};
//      formSchema['plate_search'] = {
//        title: 'Plate Search',
//        key: 'plate_search',
//        editorClass: 'input-full form-control',
//        validators: ['required'],
//        type: TextArea2,
//        template: appModel._field_template
//      };
//      
//      var FormFields = Backbone.Model.extend({
//        schema: formSchema,
//        validate: function(attrs) {
//          var errs = {};
//          if (!_.isEmpty(errs)) return errs;
//        }
//      });
//      var formFields = new FormFields();
//      var form = new Backbone.Form({
//        model: formFields,
//        template: appModel._form_template
//      });
//      
//      
//      
//      var View = Backbone.Layout.extend({
//        template: _.template(genericLayout),
//        afterRender: function(){
//          $('#resource_content').html(form.render().el);
//          form.$el.append([
//            '<button type="submit" class="btn btn-default btn-xs" ',
//            'style="width: 3em;">ok</input>',
//          ].join(''));
//          
//          var helpLink = $('<a >&nbsp;?</a>');
//          form.$el.find('[key="form-group-well_search"]').find('label').append(helpLink);
//          helpLink.click(function(e){
//            e.preventDefault();
//            var bodyMessage = [
//              'By plate number range, e.g.:',
//              '1000',
//              '1000-2000 A',
//              'B 3000-4000',
//              '5000-6000 A,B,C     ',
//              '9000,2000,3000 5000-4000 A,"b-c",D',
//            ];
//            appModel.showModalMessage({
//              title: 'Search for plate ranges using patterns',
//              body: bodyMessage.join('<br/>')
//            });
//          });
//          form.$el.find('[ type="submit" ]').click(function(e){
//            e.preventDefault();
//            var errors = form.commit({ validate: true }); 
//            if(!_.isEmpty(errors)){
//              form.$el.find('#well_search').addClass("has-error");
//              return;
//            }else{
//              form.$el.find('#well_search').removeClass("has-error");
//            }
//            
//            var search_value = form.getValue('well_search');
//            console.log('search_value', search_value);
//            search_value = search_value.split(/\n+/);
//            console.log('search_value', search_value);
//            
//            
//            var url = [self.model.resource.apiUri,self.model.key].join('/');
//            
//            function submit(override){
//              if (! _.isUndefined(override)){
//                url += '?' + appModel.API_PARAM_OVERRIDE + '=true';
//              }
//              data = { 'screener_cherry_picks': search_value };
//              var headers = {}; // can be used to send a comment
//              $.ajax({
//                url: url,     
//                cache: false,
//                contentType: 'application/json', 
//                processData: false,
//                dataType: 'json', // what is expected back from the server
//                data: JSON.stringify(data),
//                type: 'PATCH',
//                headers: headers
//              }).done(function(data, textStatus, jqXHR){
//                console.log('success', data);
//                appModel.showConnectionResult(data, {
//                  title: 'Create Screener Cherry Picks'
//                });
//                self.model.fetch({ reset: true }).done(function(){
//                  self.uriStack = ['screenercherrypicks'];
//                  // remove the child view before calling render, to prevent
//                  // it from being rendered twice, and calling afterRender twice
//                  self.removeView('#tab_container');
//                  self.render();
//                });
//              }).fail(function(jqXHR, textStatus, errorThrown) { 
//                console.log('errors', arguments);
//                
//                var jsonError = _.result(jqXHR, 'responseJSON');
//                if (!_.isUndefined(jsonError)){
//                  var error = _.result(jsonError, 'errors');
//                  var overrideFlag = _.result(error,appModel.API_PARAM_OVERRIDE);
//                  var errorMsg = _.result(error, 'screener_cherry_picks'); 
//                  var librariesToOverride = _.result(error, 'Libraries');
//                  if (!_.isUndefined(librariesToOverride)){
//                    errorMsg += '<br/>Libraries:<br/>';
//                    errorMsg += librariesToOverride.join('<br/>');
//                  }
//                  if (!_.isUndefined(overrideFlag)){
//                    appModel.showModal({
//                      title: 'Override Required',
//                      body: errorMsg,
//                      okText: 'Confirm Override',
//                      ok: function() {
//                        var override = true;
//                        submit(override);
//                      },
//                      cancel: function() {
//                        // nop
//                      }
//                    });
//                  } else {
//                    appModel.jqXHRfail.apply(this,arguments); 
//                  }
//                } else {
//                  appModel.jqXHRfail.apply(this,arguments); 
//                }
//              });
//            };
//            
//            submit();
//          });
//        }
//      });
//      var view = new View();
//
//    
//    
//    }, //_addPlateRangeDialog
    
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
    } //_addPlateRangeDialogOld
  };
  
  return PlateRangeTablePrototype;
});



