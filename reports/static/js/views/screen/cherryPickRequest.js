define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/generic_detail_layout', 
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'views/list2',
  'utils/tabbedController',
  'utils/wellSelector',
  'templates/genericResource.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, 
            DetailLayout, DetailView, EditView,
            ListView, TabbedController, WellSelector, genericLayout) {
  
  var CherryPickView = TabbedController.extend({
    
    initialize: function(args) {
      var self = this;
      self.args = args;
      this._classname = 'CherryPickRequestView';
      
//      this.tabbed_resources = _.extend({},this.cherry_pick_tabbed_resources);
//      if (! this.model.get('number_plates')>0){
//        delete this.tabbed_resources['cherrypickplates'];
//      }
//      if(self.model.get('total_number_lcps') == 0){
//        delete this.tabbed_resources['labcherrypicks'];
//      }      
      TabbedController.prototype.initialize.apply(this,arguments);
      
      if (_.isUndefined(this.screen)){
        appModel.getModel('screen', this.model.get('screen_facility_id'), 
          function(screen){
            self.screen = screen;
          });
      } else {
        if (this.screen.get('facility_id') != this.model.get('screen_facility_id')){
          throw "Wrong screen: " + this.screen.get('facility_id') + " for this Cherry Pick Request";
        }
      }
      
      _.bindAll(this, 'createScpView','createLcpView',
                'showScpSearchForm','showWellsToLeaveEmptyDialog');
      
    },

    cherry_pick_tabbed_resources: {
      detail: { 
        description: 'Details', 
        title: 'Details', 
        invoke: 'setDetail'
      },
      screenercherrypicks: {
        description : 'Screener Cherry Picks',
        title : 'Screener Cherry Picks',
        invoke : 'setScreenerCherryPicks',
        permission: 'cherrypickrequest'
      },
      labcherrypicks: {
        description : 'Lab Cherry Picks',
        title : 'Lab Cherry Picks',
        invoke : 'setLabCherryPicks',
        permission: 'cherrypickrequest'
      },
//      platemapping: {
//        description : 'Plate Mapping',
//        title : 'Plate Mapping',
//        invoke : 'setPlateMapping',
//        permission: 'cherrypickrequest'
//      },
      sourceplates: {
        description : 'Source Plates',
        title : 'Source Plates',
        invoke : 'setSourcePlates',
        permission: 'cherrypickrequest'
      },
      cherrypickplates: {
        description : 'Cherry Pick Plates',
        title : 'Cherry Pick Plates',
        invoke : 'setCherryPickPlates',
        permission: 'cherrypickrequest'
      },
            
    },      
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;

      this.tabbed_resources = _.extend({},this.cherry_pick_tabbed_resources);
      if (! this.model.get('number_plates')>0){
        delete this.tabbed_resources['cherrypickplates'];
      }
      if(self.model.get('total_number_lcps') == 0){
        delete this.tabbed_resources['labcherrypicks'];
        delete this.tabbed_resources['sourceplates'];
      }      
      
      return {
        'base_url': self.model.resource.key + '/' + self.model.key,
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    setDetail: function(delegateStack) {
      console.log('detail view');
      
      var self = this;
      var key = 'detail';
      var buttons = ['download'];
      if (appModel.hasPermission('cherrypickrequest', 'write')){
        buttons = buttons.concat(['history','edit']);
      }
      
      var editView = EditView.extend({
        afterRender: function(){
          var editForm = this;
          console.log('after render...');

          $('[name="assay_plate_type"]').change(function(){
            if($('[name="assay_plate_type"]').val().includes('96')){
              editForm.setValue('wells_to_leave_empty', null);
            }
          });
          var number_plates = self.model.get('number_plates');
          if ( !_.isNumber(number_plates) || number_plates == 0){
            var fieldKey = 'wells_to_leave_empty';
            var wellsToLeaveEmptyButton = $([
              '<a class="btn btn-default btn-sm" ',
                'role="button" id="wellsToLeaveEmpty_button" href="#">',
                'Well picker</a>'
              ].join(''));
            wellsToLeaveEmptyButton.click(function(event) {
              event.preventDefault();
              var currentValue = $('[name="'+ fieldKey + '"]').val();
              console.log('wells_to_leave_empty current value', currentValue);
              self.showWellsToLeaveEmptyDialog(editForm, currentValue, true);
            });
            this.$el.find('div[data-fields="' + fieldKey + '"]')
              .find('div[data-editor]').append(wellsToLeaveEmptyButton);
          }          
          EditView.prototype.afterRender.apply(this,arguments);
        }
      });
      var saveSuccessCallBack = function(response){
        // Override so that the CPR can be displayed in the screen resource:
        // TODO: modify generic_edit to handle custom-nested URLs for edit success URL
        var meta = _.result(response, appModel.API_RESULT_META, null);
        
        if (!_.isEmpty(meta)) {
          appModel.showJsonMessages(meta);
        }
        if (_.isUndefined(self.model) || _.isUndefined(self.model.key)){
          var model = _.result(response, appModel.API_RESULT_DATA, null);
          if (!_.isEmpty(model)){
            model = new Backbone.Model(model);
            model.key = Iccbl.getIdFromIdAttribute( model,self.model.resource );
            appModel.router.navigate([
              self.screen.resource.key,self.screen.key,'cherrypickrequest',
              model.key].join('/'), 
              {trigger:true});
          } else { 
            appModel.router.navigate([
              self.screen.resource.key,self.screen.key,'cherrypickrequest'
              ].join('/'), 
              {trigger:true});
          }
        }else{
          appModel.router.navigate([
            self.screen.resource.key,self.screen.key,'cherrypickrequest',
            self.model.key].join('/'), 
            {trigger:true});
        }
      };
      detailView = DetailView.extend({
        afterRender: function() {
          DetailView.prototype.afterRender.apply(this,arguments);
          var detail_self = this;
          var empty_wells_link = $(
            '<a href="#">' + self.model.get('wells_to_leave_empty') + '</a>');
          empty_wells_link.click(function(e){
            e.preventDefault();
            self.showWellsToLeaveEmptyDialog();
          });
          detail_self.$el.find('#wells_to_leave_empty').html(empty_wells_link);
        }
      });
      view = new DetailLayout(_.extend(self.args, { 
        model: this.model,
        uriStack: delegateStack, 
        buttons: buttons,
        saveSuccessCallBack: saveSuccessCallBack,
        EditView: editView,
        DetailView: detailView
      }));
      view.showEdit = function() {
        var innerself = this;
        appModel.initializeAdminMode(function(){
          var fields = self.model.resource.fields;
          fields['requested_by_username'].choices = 
            appModel._get_screen_members(self.screen);
          // TODO: resource/cherrypickrequest/write
          appModel.getAdminUserOptions(function(options){
            fields['volume_approved_by_username'].choices = options;
          },'cherrypickrequest');
            
          var number_plates = self.model.get('number_plates');
          if ( _.isNumber(number_plates) && number_plates != 0){
            _.each([
              'assay_plate_type',
              'keep_source_plate_cherry_picks_together',
              'is_randomized_assay_plate_layout',
              'wells_to_leave_empty',
              'transfer_volume_per_well_requested',
              'transfer_volume_per_well_approved'],
              function(disallowed_field){
                fields[disallowed_field]['editability'] = [];
              }
            );
          }
          
          DetailLayout.prototype.showEdit.apply(innerself,arguments);
        });  
      };
      this.tabViews[key] = view;
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
      
    }, // end setDetail

    showWellsToLeaveEmptyDialog: function(editForm, currentValue, editable){
      var self = this;
      var rowCollection = new Backbone.Collection();
      if (_.isUndefined(editable)){
        editable = false;
      }
      if (_.isUndefined(currentValue)){
        currentValue = self.model.get('wells_to_leave_empty');
      }
      var assay_plate_type = $('[name="assay_plate_type"]').val();
      if (_.isEmpty(assay_plate_type)){
        assay_plate_type = self.model.get('assay_plate_type');
        if (_.isEmpty(assay_plate_type)){
          appModel.error('must define the assay plate type');
          return;
        }
      }
      var plateSize = 384;
      if (assay_plate_type.includes('96')){
        plateSize = 96;
      }
      try {
        var wellSelector = new WellSelector({
          plateSize: plateSize,
          wellSelections: currentValue,
          editable: editable
        });
        var el = wellSelector.render().el;

        var title = 'Select wells to leave empty';
        if (!editable){
          title = "Wells to leave empty";
        }

        var modalDialog = appModel.showModal({
          ok_only: !editable,
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
            
            if (!editable) return;
            
            console.log('selected wells: ', wellSelector.getSelectedWells());
            editForm.setValue('wells_to_leave_empty', 
              wellSelector.getSelectedWells().join(', '));
            self.model.set(
              'wells_to_leave_empty', 
              wellSelector.getSelectedWells().join(','));
          }              
        });
        if (! editable){
          modalDialog.$el.find('#modal-cancel').hide();
        }
        
      } catch(e) {
        console.log('uncaught error: ' + e);
        appModel.error('Error with selection: ' + e);
        return;
      }
    }, // end showWellsToLeaveEmptyDialog
    
    setSourcePlates: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'source_plate'].join('/');
      var resource = appModel.getResource('librarycopyplate')
      var fields = resource.fields;
      var includes = [];
      _.each([
        'status_date',
        'mg_ml_concentration',
        'molar_concentration',
        'min_mg_ml_concentration',
        'min_molar_concentration',
        'max_mg_ml_concentration',
        'max_molar_concentration',
        'first_date_screened',
        'last_date_screened'],
        function(hidden_field){
          fields[hidden_field]['visibility'] = [];
          includes.unshift('-' + hidden_field);
        }
      );
      _.each([
        'copy_usage_type'],
        function(shown_field){
          fields[shown_field]['visibility'] = ['l','d'];
          includes.unshift(shown_field);
        }
      );
      resource.options.includes = includes;
      
      var extraControls = [];
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        resource: resource,
        includes: includes,
        url: url,
        extraControls: extraControls
      });
      Backbone.Layout.setupView(view);
//      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Source Plates for CPR: ' + self.model.key + '</H4>');
      });
    },
    
    showPlateMappingGrid: function(url){
      var self = this;
      var title = Iccbl.formatString(
        'Screen: {screen_facility_id}, ' +
        'CPR: {cherry_pick_request_id}, Cherry Pick Assay Plate:',
        self.model);
      var assay_plate_type = self.model.get('assay_plate_type');
      var plateSize = 384;
      if (assay_plate_type.includes('96')){
        plateSize = 96;
      }
      Iccbl.getCollectionOnClient(url,
        function(collection){
          if(collection && !collection.isEmpty()){
          } else {
            appModel.error('No plating found for CPR: ' + self.model.key);
          }
          
          var wellSelectionMap = collection.groupBy(function(model){
            return model.get('cherry_pick_plate_number');
          });
          var number_of_plates = _.keys(wellSelectionMap).length;
          function getWellIds(cherry_pick_plate_number){
            if (!_.has(wellSelectionMap,cherry_pick_plate_number)){
              appModel.error(
                'Plate not found: ' + cherry_pick_plate_number +
                ', available: ' + _.keys(wellSelectionMap).join(','));
              return null;
            }
            var lcps = wellSelectionMap[cherry_pick_plate_number];
            well_list = _.map(lcps, function(lcp){
              return lcp.get('destination_well');
            });
            console.log('well list:', well_list);
            return well_list.join(',');
          };
          
          var current_plate_number = 1;
          var wellSelector = new WellSelector({
            plateSize: plateSize,
            wellSelections: getWellIds(current_plate_number),
            editable: false
          });
          var el = wellSelector.render().el;
          var modalDialog = appModel.showModal({
            buttons_on_top: true,
            css: { 
                display: 'table',
                width: 'auto'
              },
            css_modal_content: {
              overflow: 'hidden'
            },
            css_modal_body: {
              'padding-top': '0px'
            },
            okText: 'next',
            cancelText: 'done',
            view: el,
            title: title + ' 1',
            ok: function(e) {
              e.preventDefault();
              
              current_plate_number += 1;
              if (current_plate_number > number_of_plates){
                current_plate_number = 1;
              }
              console.log('next...', current_plate_number);
              
              modalDialog.$el.find('#modal_title').html(
                title + ' ' + current_plate_number);
              var wellSelector = new WellSelector({
                plateSize: plateSize,
                wellSelections: getWellIds(current_plate_number),
                editable: false
              });
              var el = wellSelector.render().el;
              modalDialog.$el.find('.modal-body').empty();
              modalDialog.$el.find('.modal-body').append(el);
              return false;
            }              
          });
        }
      );
      
    }, // end showPlateMappingGrid
    
    setCherryPickPlates: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,
        'cherry_pick_plate'].join('/');
      var resource = appModel.getResource('cherrypickassayplate');
      var selectionCollection = null;

      if (appModel.hasPermission('cherrypickrequest', 'write')){

        // Set up the UI to select for plating and screening
          
        var selectAllControl = $([
          '<div class="checkbox">',
          '<label  ',
          'title="Select all plates (on all pages)" >',
          'Select All<br/>',
          '<input id="selectAll" type="checkbox">',
          '</label><br>&nbsp;',
          '</div>'
        ].join(''));
        resource['fields']['selected'] = {
          ordinal: 0,
          visibility: ['l'],
          editability: ['l'],
          data_type: 'boolean',
          description: 'Selected', 
          title: 'Selected',
          key: 'selected',
          headerCell: Backgrid.HeaderCell.extend({
            render: function() {
              this.$el.empty();
              this.$el.append(selectAllControl);
              return this;
            }
          })
        };
      
        var Collection = Iccbl.CollectionOnClient.extend({
          // explicitly define the id so that collection compare & equals work
          modelId: function(attrs) {
            return Iccbl.getIdFromIdAttribute( attrs, resource);
          },
          url: url 
        });
        
        selectionCollection = new Collection();
        
        // Fetch the entire cpap collection to use for selection tracking
        selectionCollection.fetch({

          data: { limit: 0 },

          success: function(collection, response) {
            
            var extraControls = [];
            var setPlatedDiv = $(
              '<span title="Set the plating activity date for the selected plates"/>'
            );
            var setPlatedButton = $([
              '<a class="btn btn-default btn-sm pull-down disabled" ',
              'role="button" id="set_plated_button" ',
              'href="#">',
              'Set Plated</a>'
            ].join(''));
            setPlatedDiv.append(setPlatedButton);
            extraControls.push(setPlatedDiv);
            var setScreenedDiv = $(
              '<span title="Set the plating activity date for the selected plates"/>'
            );
            var setScreenedButton = $([
              '<a class="btn btn-default btn-sm pull-down disabled" ',
              'role="button" id="set_screened_button" ',
              'title="set the screening activity date for the selected plates" href="#">',
              'Set Screened</a>'
            ].join(''));
            extraControls.push(setScreenedButton);
            
            setPlatedButton.click(function(e){
              e.preventDefault();
              self.showPlatedDateDialog(selectionCollection);
            });
            setScreenedButton.click(function(e){
              e.preventDefault();
              self.showScreenedDateDialog(selectionCollection);
            });
            
            view = createCherryPickAssayPlateView(extraControls);
            
            view.collection.on('change', function(model){
              // record user selections in the selectionCollection
              selectionCollection.add(model, {silent: true, merge: true});
              if (model.get('selected')==false){
                selectAllControl.find('input[type="checkbox"]').prop('checked',false);
              }
              if(!_.isEmpty(selectionCollection.where({selected: true }))){
                setPlatedButton.removeClass('disabled');                
                setScreenedButton.removeClass('disabled');                
              } else {
                setPlatedButton.addClass('disabled');                
                setScreenedButton.addClass('disabled');                
              }
              
            });
            
            selectionCollection.on('change', function(model){
              // catch changes to the selectionCollection by "selectAll"
              view.collection.get(model).set('selected', model.get('selected'));
              
              if(!_.isEmpty(selectionCollection.where({selected: true }))){ 
                setPlatedButton.removeClass('disabled');                
                setScreenedButton.removeClass('disabled');                
              } else {
                setPlatedButton.addClass('disabled');                
                setScreenedButton.addClass('disabled');                
              }
            });

            // Make sure that on reset actions (page changes), selections are persisted
            view.collection.on('reset', function(){
              selectionCollection.each(function(model){
                var retrievedModel = view.collection.get(model);
                if (!_.isUndefined(retrievedModel)){
                    retrievedModel.set('selected', model.get('selected'));
                }
             });
            });
            
            selectAllControl.click(function(e){
              if (e.target.checked) {
                selectionCollection.each(function(model){
                  model.set('selected', true);
                });
              } else {
                selectionCollection.each(function(model){
                  model.set('selected', false);
                });
              }
            });
          },
          always: function(){
            console.log('done: ');
          }
        }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      

        
      } else {
        createCherryPickAssayPlateView();
      }        
      
      function createCherryPickAssayPlateView(extraControls){
        if (_.isUndefined(extraControls)){
          extraControls = [];
        }
        var downloadPlateMappingButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="plate_mapping_file" href="#">',
            'Download Plate Mapping files</a>'
          ].join(''));
        extraControls.push(downloadPlateMappingButton);
        var showPlateMappingButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
             'title="Show the plate mapping in a grid (available if plates are assigned)"',
            'role="button" id="plate_mapping_grid" href="#">',
            'Show Plate Mapping</a>'
          ].join(''));
        extraControls.push(showPlateMappingButton);

        var view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: resource,
          resource: resource,
          url: url,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
//        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.listenTo(view, 'afterRender', function(event) {
          view.$el.find('#list-title').show().append(
            '<H4 id="title">Assay Plates for CPR: ' + self.model.key + '</H4>');
        });
        
        downloadPlateMappingButton.click(function(e){
          e.preventDefault();
          console.log('download plate mapping file');
          var url = [self.model.resource.apiUri,self.model.key,
                     'plate_mapping_file'].join('/');
          appModel.downloadUrl(self.$el,url);
        });
        showPlateMappingButton.click(function(e){
          e.preventDefault();
          var url = [self.model.resource.apiUri,self.model.key,
            'lab_cherry_pick'].join('/');
          self.showPlateMappingGrid(url);
        });
        return view;
      };
    }, // end setCherryPickPlates
    
    showScreenedDateDialog: function(selectionCollection){
      var self = this;
      
      var selectedEntries = selectionCollection.where(
        {selected: true }); //, plating_date: null 
      if (_.isEmpty(selectedEntries)){
        console.log('no selections found');
        return;
      } 
      var notPlatedEntries = _.filter(selectedEntries,
        function(model){
          return !model.has('plating_date');
        });
      if (!_.isEmpty(notPlatedEntries)){
        appModel.showModalMessage({
          title: 'Plates must be plated before being screened',
          body: 
            'Selected plates have not yet been plated: ' +
            _.map(notPlatedEntries,
              function(model){ 
              return model.get('plate_ordinal'); 
            }).join(','),
        });
        return;
      }
      var alreadyScreenedEntries = _.filter(selectedEntries,
        function(model){
          return model.has('screening_date');
        });
      if (!_.isEmpty(alreadyScreenedEntries)){
        appModel.showModal({
          title: 'Already screened',
          body: 
            'Selected plates are already screened: ' +
            _.map(alreadyScreenedEntries,
              function(model){ 
              return model.get('plate_ordinal'); 
            }).join(',') +
            ', change the last date screened?',
          ok: function(event){
            $('#modal').modal('hide');
            dialogFunction();
          }
        })
      }else{
        dialogFunction();
      }
      
      function dialogFunction(){
        appModel.initializeAdminMode(function(){
              
          // Build the form model
          var FormFields = Backbone.Model.extend({
            schema: {
              screening_date: {
                title: 'Date',
                key: 'screening_date',
                type: EditView.DatePicker,
                validators: ['required'], 
                template: appModel._field_template
              },
              screened_by_username: {
                title: 'Screened By',
                key: 'screened_by_username',
                type: EditView.ChosenSelect,
                editorClass: 'chosen-select',
                options: appModel._get_screen_members(self.screen),
                validators: ['required'], 
                template: appModel._field_template
              },
              comments: {
                title: 'Comments',
                key: 'comments',
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
          var _form = form.render();
          _form.$el.find('.chosen-select').chosen({
            disable_search_threshold: 3,
            width: '100%',
            allow_single_deselect: true,
            search_contains: true
            });
    
          appModel.showModal({
            okText: 'ok',
            ok: function(e){
              e.preventDefault();
              var errors = form.commit({ validate: true }); // runs schema and model validation
              if(!_.isEmpty(errors) ){
                console.log('form errors, abort submit: ',errors);
                _.each(_.keys(errors), function(key){
                  $('[name="'+key +'"').parents('.form-group').addClass('has-error');
                });
                if (_.has(errors,'_others')){
                  $errorDiv = $('<div class="panel text-danger" />');
                  _.each(errors['_others'], function(otherError){
                    _.each(_.values(otherError), function(errMsg){
                      $errorDiv.append('Error: ' + errMsg );
                    });
                  });
                  _form.$el.append($errorDiv);
                }
                return false;
              }            
              var values = _form.getValue();
              console.log('form submitted', values, selectedEntries);
              
              var plate_updates = [];
                    
              var SubmitCollection = Backbone.Collection.extend({
                url: selectionCollection.url
              });
              var submitCollection = new SubmitCollection(selectedEntries);
              submitCollection.each(function(model){
                 model.set('screening_date', values['screening_date']);
                 model.set('screened_by_username', values['screened_by_username']);
              });
  
              var headers = {};
              headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
              
              submitCollection.sync(
                'patch', submitCollection, { headers: headers })
                .done(function(data, textStatus, jqXHR){
                  appModel.showConnectionResult(data, {
                    title: 'Screening'
                  });
                  self.model.fetch({ reset: true }).done(function(){
                    self.uriStack = ['cherrypickplates'];
                    // Remove the child view before calling render, to prevent
                    // it from being rendered twice, and calling afterRender twice
                    self.removeView('#tab_container');
                    self.render();
                  });
                }).fail(function(jqXHR, textStatus, errorThrown){
                  console.log('fail', arguments);
                  appModel.jqXHRfail.apply(this,arguments); 
                });
              
            },
            view: _form.el,
            title: 'Set Screened Date for selected plates'
          });
        });      
        
      };
    }, // end showScreenedDateDialog
    
    showPlatedDateDialog: function(selectionCollection){
      var self = this;
      
      var selectedEntries = selectionCollection.where(
        {selected: true }); //, plating_date: null 
      if (_.isEmpty(selectedEntries)){
        console.log('no selections found');
        return;
      } 
      var alreadyPlatedEntries = _.filter(selectedEntries,
        function(model){
          return model.has('plating_date');
        });
      if (!_.isEmpty(alreadyPlatedEntries)){
        appModel.showModal({
          title: 'Already plated',
          body: 
            'Selected plates are already plated: ' +
            _.map(alreadyPlatedEntries,
              function(model){ 
              return model.get('plate_ordinal'); 
            }).join(',') +
            ', change the plating date?',
          ok: function(event){
            $('#modal').modal('hide');
            appModel.getAdminUserOptions(dialogFunction, 'cherrypickrequest');    
          }
        })
      }else{
        appModel.getAdminUserOptions(dialogFunction, 'cherrypickrequest');    
      }
      
      function dialogFunction(cherryPickUserOptions){
        // Build the form model
        var FormFields = Backbone.Model.extend({
          schema: {
            plating_date: {
              title: 'Date',
              key: 'plating_date',
              type: EditView.DatePicker,
              validators: ['required'], 
              template: appModel._field_template
            },
            plated_by_username: {
              title: 'Plated By',
              key: 'plated_by_username',
              type: EditView.ChosenSelect,
              editorClass: 'chosen-select',
              options: cherryPickUserOptions,
              validators: ['required'], 
              template: appModel._field_template
            },
            comments: {
              title: 'Comments',
              key: 'comments',
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
        var _form = form.render();
        _form.$el.find('.chosen-select').chosen({
          disable_search_threshold: 3,
          width: '100%',
          allow_single_deselect: true,
          search_contains: true
          });
  
        appModel.showModal({
          okText: 'ok',
          ok: function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }); // runs schema and model validation
            if(!_.isEmpty(errors) ){
              console.log('form errors, abort submit: ',errors);
              _.each(_.keys(errors), function(key){
                $('[name="'+key +'"').parents('.form-group').addClass('has-error');
              });
              if (_.has(errors,'_others')){
                $errorDiv = $('<div class="panel text-danger" />');
                _.each(errors['_others'], function(otherError){
                  _.each(_.values(otherError), function(errMsg){
                    $errorDiv.append('Error: ' + errMsg );
                  });
                });
                _form.$el.append($errorDiv);
              }
              return false;
            }            
            var values = _form.getValue();
            console.log('form submitted', values, selectedEntries);
            
            var plate_updates = [];
                  
            var SubmitCollection = Backbone.Collection.extend({
              url: selectionCollection.url
            });
            var submitCollection = new SubmitCollection(selectedEntries);
            submitCollection.each(function(model){
               model.set('plating_date', values['plating_date']);
               model.set('plated_by_username', values['plated_by_username']);
            });

            var headers = {};
            headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
            
            submitCollection.sync(
              'patch', submitCollection, { headers: headers })
              .done(function(data, textStatus, jqXHR){
                appModel.showConnectionResult(data, {
                  title: 'Plating'
                });
                self.model.fetch({ reset: true }).done(function(){
                  self.uriStack = ['cherrypickplates'];
                  // Remove the child view before calling render, to prevent
                  // it from being rendered twice, and calling afterRender twice
                  self.removeView('#tab_container');
                  self.render();
                });
              }).fail(function(jqXHR, textStatus, errorThrown){
                console.log('fail', arguments);
                appModel.jqXHRfail.apply(this,arguments); 
              });
            
          },
          view: _form.el,
          title: 'Set Plated Date for selected plates'
        });      
      };
    }, // end showPlatedDateDialog

    
//    /**
//     * LCP Plate Mapping (and QC) view
//     */
//    setPlateMapping: function(delegateStack) {
//      var self = this;
//      var url = [self.model.resource.apiUri, 
//                 self.model.key,
//                 'lab_cherry_pick_plating'].join('/');
//      var extraControls = [];
//      var downloadPlateMappingButton = $([
//        '<a class="btn btn-default btn-sm pull-down" ',
//           'title="Download plate mapping file (available if plates are assigned)"',
//          'role="button" id="plate_mapping_file" href="#">',
//          'Download Plate Mapping files</a>'
//        ].join(''));
//      extraControls.push(downloadPlateMappingButton);
//      var showPlateMappingButton = $([
//        '<a class="btn btn-default btn-sm pull-down" ',
//           'title="Show the plate mapping in a grid (available if plates are assigned)"',
//          'role="button" id="plate_mapping_grid" href="#">',
//          'Show Plate Mapping</a>'
//        ].join(''));
//      extraControls.push(showPlateMappingButton);
//
//      function createPlateMappingView(resource) {
//        
//        view = self.createLcpView(delegateStack, resource,url,extraControls)
//        
//        self.listenTo(view, 'afterRender', function(event) {
//          view.$el.find('#list-title').show().append(
//            '<H4 id="title">Plate Mapping for CPR: ' + self.model.key + '</H4>');
//        });
//        downloadPlateMappingButton.click(function(e){
//          e.preventDefault();
//          console.log('download plate mapping file');
//          var url = [self.model.resource.apiUri,self.model.key,
//                     'plate_mapping_file'].join('/');
//          appModel.downloadUrl(self.$el,url);
//        });
//        
//        showPlateMappingButton.click(function(e){
//          e.preventDefault();
//          self.showPlateMappingGrid(url);
//        });
//      };
//      var schemaUrl = url + '/schema';
//      appModel.getResourceFromUrl(schemaUrl, createPlateMappingView);
//      
//    }, // end setPlateMapping
    
    
    /**
     * Lab Cherry Picks view
     */
    setLabCherryPicks: function(delegateStack) {

      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,
        'lab_cherry_pick'].join('/');
      var schemaUrl = [
        self.model.resource.apiUri,self.model.key,'lab_cherry_pick',
        'schema'].join('/');
      
      var plate_mapping_controls = [];
      
      if (self.model.get('number_plates') != 0){
        url = [self.model.resource.apiUri, 
                   self.model.key,
                   'lab_cherry_pick_plating'].join('/');
        // NOTE: 201703 - display showPlateMapping
        var downloadPlateMappingButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
             'title="Download plate mapping file (available if plates are assigned)"',
            'role="button" id="plate_mapping_file" href="#">',
            'Download Plate Mapping files</a>'
          ].join(''));
        plate_mapping_controls.push(downloadPlateMappingButton);
        var showPlateMappingButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
             'title="Show the plate mapping in a grid (available if plates are assigned)"',
            'role="button" id="plate_mapping_grid" href="#">',
            'Show Plate Mapping</a>'
          ].join(''));
        plate_mapping_controls.push(showPlateMappingButton);
        downloadPlateMappingButton.click(function(e){
          e.preventDefault();
          console.log('download plate mapping file');
          var url = [self.model.resource.apiUri,self.model.key,
                     'plate_mapping_file'].join('/');
          appModel.downloadUrl(self.$el,url);
        });
        
        showPlateMappingButton.click(function(e){
          e.preventDefault();
          self.showPlateMappingGrid(url);
        });
        schemaUrl = url + '/schema';
      }      

      function createView(resource){
        view = self.createLcpView(delegateStack,resource, url);
        self.listenTo(view, 'afterRender', function(event) {
          
          if (!_.isEmpty(plate_mapping_controls)){
            view.$el.find('#list-container').prepend(plate_mapping_controls);
          }
        });
      };
      appModel.getResourceFromUrl(schemaUrl, createView);
    }, // end setLabCherryPicks
    

    createLcpView: function (delegateStack, resource, url, extraControls){
      var self = this;
      if (_.isUndefined(extraControls)){
        extraControls = [];
      }

      var checkboxDiv = $([
          '<div id="show_input_group" class="input-group"></div>'
        ].join(''));
      var showCopyWellsControl = $([
          '<label class="checkbox-inline" ', 
          ' style="margin-left: 10px;" ',
          'title="Show all available wells from Cherry Pick Source Plate copies" >',
          '  <input id="show_copy_wells" type="checkbox">Available',
          '</label>'
        ].join(''));
      var showAllCopyWellsControl = $([
          '<label class="checkbox-inline" ',
          'title="Show all available and retired wells ',
          'from Cherry Pick Source and Library Screening copies">',
          '  <input id="show_available_and_retired_copy_wells" ',
          '     type="checkbox">Available and Retired',
          '</label>'
        ].join(''));
      var showUnfulfilledWellsControl = $([
          '<label class="checkbox-inline" ',
          ' title="Show rows for unfulfilled picks only" >',
          '  <input id="showUnfulfilled" type="checkbox">Unfulfilled',
          '</label>'
        ].join(''));
      var showInsufficientWellsControl = $([
          '<label class="checkbox-inline" ',
          ' title="Show rows for picks with insufficient volume only" >',
          '  <input id="showInsufficient" type="checkbox">Insufficient Volume',
          '</label>'
        ].join(''));
      var showManuallySelectedWellsControl = $([
          '<label class="checkbox-inline" ',
          ' title="Show rows for manually selected picks only" >',
          '  <input id="showManuallySelected" type="checkbox">Manually Selected',
          '</label>'
        ].join(''));
      checkboxDiv.append(showCopyWellsControl);
      checkboxDiv.append(showAllCopyWellsControl);
      checkboxDiv.append(showUnfulfilledWellsControl);
      checkboxDiv.append(showInsufficientWellsControl);
      checkboxDiv.append(showManuallySelectedWellsControl);
      checkboxDiv.prepend('<label for="show_input_group">show</label>');
      extraControls.push(checkboxDiv);
      
      var setSelectedLcpButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="save_button_lcp_selected" href="#">',
          'Save Selections</a>'
        ].join(''));
      extraControls.push(setSelectedLcpButton);
      var updateSelectedLcpButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="update_selected_button" href="#">',
          'Update Selections</a>'
        ].join(''));
      extraControls.push(updateSelectedLcpButton);
      var cancelSelectedButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="cancel_selected_button" href="#">',
          'Cancel Selections</a>'
        ].join(''));
      extraControls.push(cancelSelectedButton);
      var reserveAndMapSelectedButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="reserve_map_selected_button" href="#">',
          'Reserve Selections and Map to Plates</a>'
        ].join(''));
      extraControls.push(reserveAndMapSelectedButton);
      var deleteLcpsButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="deleteLcpsButton" href="#">',
          'Delete Lab Cherry Picks</a>'
        ].join(''));
      extraControls.push(deleteLcpsButton);
      var cancelReservation = $([
          '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="cancel_reservation" href="#">',
          'Cancel Reservation and Delete Plating assignments</a>'
        ].join(''));
      extraControls.push(cancelReservation);

      if (self.model.get('number_plates') == 0){
        deleteLcpsButton.show();
        reserveAndMapSelectedButton.show();
      } else {
        cancelReservation.show();
      }

      ///// Library and Plate comments /////
      
      function parseComments(comment_array){
        return _.map(
          comment_array,
          function(comment){
            // FIXME: make this error tolerant
            comment_array = comment.split('$');
            return '(' + comment_array[0] + ') ' +
              Iccbl.getDateString(comment_array[1]) + 
              ': ' + comment_array[2];
          }).join('\n');
      };
      function createCommentIcon(comments, title){
        var comment_icon = $(
          '<span class="glyphicon glyphicon-comment" ' +
          'style="color: lightgray; " ></span>');
        
        var rows = 10;
        var buttons_on_top = false;
        if (comments.length > 300){
          rows = 30;
          buttons_on_top = true;
        }
        comment_icon.click(function(e){
          e.preventDefault();
          var body = $('<textarea class="input-full" rows=' + rows + ' ></textarea>');
          body.val(comments);
          appModel.showModalMessage({
            title: title,
            view: body,
            buttons_on_top: buttons_on_top
          });
        });
        return comment_icon;
      };
      resource.fields['library_plate']['backgridCellType'] = 
        Iccbl.LinkCell.extend({
          render: function(){
            var self = this;
            Iccbl.LinkCell.prototype.render.apply(this, arguments);
            var comments = this.model.get('library_plate_comment_array');
            if (!_.isEmpty(comments)){
              comments = parseComments(comments);
              this.$el.attr('title',comments);
              this.$el.append(createCommentIcon(
                comments,
                'Comments for Plate: ' 
                  + self.model.get('library_short_name') + '/' 
                  + self.model.get('source_copy_name')  + '/'
                  + self.model.get('library_plate')));
            }
            return this;
          }
        },resource.fields['library_plate'].display_options);

      resource.fields['source_copy_name']['backgridCellType'] = 
        Iccbl.LinkCell.extend({
          render: function(){
            var self = this;
            Iccbl.LinkCell.prototype.render.apply(this, arguments);
            var comments = this.model.get('source_copy_comments');
            if (!_.isEmpty(comments)){
              this.$el.attr('title', comments);
              this.$el.append(createCommentIcon(
                comments,
                'Comments for Copy: ' 
                  + self.model.get('library_short_name') + '/'
                  + self.model.get('source_copy_name')
                ));
            }
            return this;
          }
        },resource.fields['source_copy_name'].display_options);

      var library_link_cell = Iccbl.LinkCell.extend({
          render: function(){
            var self = this;
            Iccbl.LinkCell.prototype.render.apply(this, arguments);
            var comments = this.model.get('library_comment_array');
            if (!_.isEmpty(comments)){
              comments = parseComments(comments);
              this.$el.attr('title',comments);
              this.$el.append(createCommentIcon(
                comments,
                'Comments for library: ' + self.model.get('library_short_name')));
            }
            return this;
          }
        });
      resource.fields['library_name']['backgridCellType'] =
        library_link_cell.extend(resource.fields['library_name'].display_options);
      resource.fields['library_short_name']['backgridCellType'] =
        library_link_cell.extend(resource.fields['library_short_name'].display_options);
      ///// end Library and Plate comments /////

      if(self.model.get('has_pool_screener_cherry_picks') === true){
        resource.fields['pool_reagent_vendor_id']['visibility'] = ['l','d'];
      }
      
      // Track user LCP selections
      var Collection = Backbone.Collection.extend({
        // explicitly define the id so that collection compare & equals work
        modelId: function(attrs) {
          return Iccbl.getIdFromIdAttribute( attrs, resource);
        }
      })
      var lcpSelectionUpdateCollection = new Collection();
      
      // Custom LCP Selection Row with styling
      // Reference to previous row to use for row styling
      var previousSourceWell = null;
      var SelectedLcpRow = Backgrid.Row.extend({

        initialize: function () {
          var self = this;
          SelectedLcpRow.__super__.initialize.apply(this, arguments);
          this.listenTo(
              this.model.collection, 
              'show_available_and_retired_copy_wells', 
              function(){
                self._setStyle();
              }
          );
          this.listenTo(
              this.model.collection, 
              'show_copy_wells', 
              function(){
                self._setStyle();
              }
          );
        },
        
        _setStyle: function(){
          var searchHash = this.model.collection.listModel.get('search');
          if (_.result(searchHash,'show_copy_wells')=='true'
            || _.result(searchHash,'show_available_and_retired_copy_wells')=='true' ){
            if (this.model.has('selected_copy_name')) {
              if (this.model.get('selected_copy_name')==this.model.get('source_copy_name')) {
                this.$el.addClass('selected');
              }
            }
            if (this.model.get('source_well_id')!==previousSourceWell){
              previousSourceWell = this.model.get('source_well_id');
              this.$el.closest('tr').addClass('selector_row_group');
            }
          } else {
            this.$el.closest('tr').removeClass('selector_row_group');
          }
        },
        
        render: function() {
          SelectedLcpRow.__super__.render.apply(this, arguments);
          this._setStyle();
          return this;
        }
      });
      // Set LCP List View
      var ListViewSelect = ListView.extend({
        afterRender: function(){
          ListViewSelect.__super__.afterRender.apply(this, arguments);
          // Backgrid specific: all column headers are given a class for their
          // column name: in this case the "selected" class has other meaning, 
          // so remove it
          this.$('th').removeClass('selected');
          this.grid.$el.addClass('selector_table');
          return this;
        }
      });

      var view = new ListViewSelect({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        row: SelectedLcpRow,
        extraControls: extraControls
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Lab Cherry Picks for : ' + self.model.key + '</H4>');
        view.$el.find('#extra_controls').removeClass('col-sm-5');
        view.$el.find('#extra_controls').addClass('col-sm-8');
        view.$el.find('#list_controls').removeClass('col-sm-7');
        view.$el.find('#list_controls').addClass('col-sm-4');
      });
    
      var initialSearchHash = view.listModel.get('search');
      if (_.has(initialSearchHash, 'show_copy_wells')
          && initialSearchHash.show_copy_wells.toLowerCase()=='true') {
        showCopyWellsControl.find('input[type="checkbox"]').prop('checked',true);
      }
      if (_.has(initialSearchHash, 'show_available_and_retired_copy_wells')
          && initialSearchHash.show_available_and_retired_copy_wells.toLowerCase()=='true') {
        showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',true);
        showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
      }
      if (_.has(initialSearchHash, 'show_unfulfilled')
          && initialSearchHash.show_unfulfilled.toLowerCase()=='true') {
        showUnfulfilledWellsControl.find('input[type="checkbox"]').prop('checked',true);
      }
      if (_.has(initialSearchHash, 'show_insufficient')
          && initialSearchHash.show_insufficient.toLowerCase()=='true') {
        showInsufficientWellsControl.find('input[type="checkbox"]').prop('checked',true);
      }
      if (_.has(initialSearchHash, 'show_manual')
          && initialSearchHash.show_manual.toLowerCase()=='true') {
        showManuallySelectedWellsControl.find('input[type="checkbox"]').prop('checked',true);
      }

      // Manage LCP selection updates
      view.collection.on('add', function(model){
        // cache the 'selected' property for update management
        model.set(
          { selected_on_server: model.get('selected') },
          { silent: true }
        );
      });
      view.collection.on('change', function(model){
        console.log('collection changed', arguments);
        if (!model.has('selected_on_server')){
          // Block changes caused by changing the field schema with other actions 
          // that happen before the selected_on_server flag is set
          return;
        }
//        if (!(showCopyWellsControl.find('input[type="checkbox"]').prop('checked')
//            || showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked'))){
//          return;
//        }

        var current_copy_name = model.get('source_copy_name');
        var source_well_id = model.get('source_well_id');
        var source_copy_id = model.get('source_copy_id');
        var selected_copy_name = model.get('selected_copy_name');

        if (model.get('selected') != model.get('selected_on_server')) {
          if (self.model.get('number_plates') != 0){
            if (_.isEmpty(selected_copy_name)){
              appModel.showModalMessage({
                title: 'Note:',
                body: 'New selections can not be made unless plating '+
                  'assignments are deallocated'
              });
              model.set('selected',false);
            } else {
              lcpSelectionUpdateCollection.add(model);
            }
          } else {
            lcpSelectionUpdateCollection.add(model);
          }
        } else {
          lcpSelectionUpdateCollection.remove(model);
        }
        if (model.get('selected')) {
          // unselect others in the group
          _.each(view.collection.where({'source_well_id':source_well_id}),
            function(collection_model){
              if (collection_model.get('source_copy_name')!==current_copy_name){
                collection_model.set({'selected': false});
              }
          });
        }
      });
      
      lcpSelectionUpdateCollection.on('update reset', function(){
        if (!lcpSelectionUpdateCollection.isEmpty()){
          appModel.setPagePending(null, 'Selection updates are unsaved, continue anyway?');
          cancelSelectedButton.show();
          if (self.model.get('number_plates') == 0){
            setSelectedLcpButton.show();
          }else{
            updateSelectedLcpButton.show();
          }
          reserveAndMapSelectedButton.hide();
          deleteLcpsButton.hide();
          cancelReservation.hide();
        } else {
          setSelectedLcpButton.hide();
          updateSelectedLcpButton.hide();
          cancelSelectedButton.hide();
          if (self.model.get('number_plates') == 0){
            deleteLcpsButton.show();
            reserveAndMapSelectedButton.show();
          } else {
            cancelReservation.show();
          }
          appModel.clearPagePending();
        }
      });
      
      // Make sure that on reset actions (page changes), selections are persisted
      view.collection.on('reset', function(){
        // Note: on "reset" the "add" methods aren't being called
        view.collection.each(function(model){
          model.set(
            { selected_on_server: model.get('selected') },
            { silent: true }
          );
        });
        lcpSelectionUpdateCollection.each(function(model){
          var retrievedModel = view.collection.get(model);
          if (!_.isUndefined(retrievedModel)){
              retrievedModel.set('selected', model.get('selected'));
          }
          if (model.get('selected')) {
            // unselect others in the group
            _.each(view.collection.where({'source_well_id':model.get('source_well_id')}),
              function(collection_model){
                if (collection_model.get('source_copy_name')!==model.get('source_copy_name')){
                  collection_model.set({'selected': false});
                }
            });
          }
       });
      });
      
      setSelectedLcpButton.click(function(e){
        e.preventDefault();
        console.log('selection updates', lcpSelectionUpdateCollection);
        if (lcpSelectionUpdateCollection.isEmpty()) {
          appModel.error('No changes to save');
          return;
        }
        function processClick(formValues){
          console.log('form values', formValues);
          var headers = {};
          headers[appModel.HEADER_APILOG_COMMENT] = formValues['comments'];
          lcpSelectionUpdateCollection.url = view.collection.url;
          lcpSelectionUpdateCollection.sync(
            'patch', lcpSelectionUpdateCollection, { headers: headers })
            .done(function(data, textStatus, jqXHR){
              appModel.showConnectionResult(data, {
                title: 'Lab Cherry Pick copy selection updates'
              });
              setSelectedLcpButton.hide();
              lcpSelectionUpdateCollection.reset(null); // clear
              // On success, clear all the buttons
              var originalSearchHash = _.clone(view.listModel.get('search'));
              var searchHash = _.clone(view.listModel.get('search'));
              delete searchHash['show_copy_wells'];
              delete searchHash['show_available_and_retired_copy_wells'];
              delete searchHash['show_unfulfilled'];
              delete searchHash['show_insufficient'];
              searchHash['show_manual'] = 'true';
              showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showUnfulfilledWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showInsufficientWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showManuallySelectedWellsControl.find('input[type="checkbox"]').prop('checked',true);
              view.listModel.set('search',searchHash);
              if (_.isEqual(originalSearchHash,searchHash)){
                view.collection.fetch({ reset: true });
              }
    
            }).fail(function(jqXHR, textStatus, errorThrown){
              console.log('fail', arguments);
              appModel.jqXHRfail.apply(this,arguments); 
            });
        };
        
        appModel.showOkCommentForm({
          ok: processClick,
          title: 'Save Lab Cherry Pick selections?'
        });
      });
      
      cancelSelectedButton.click(function(e){
        e.preventDefault();
        console.log('cancel selection updates', lcpSelectionUpdateCollection);
        lcpSelectionUpdateCollection.each(function(model){
          var retrievedModel = view.collection.get(model);
          if (!_.isUndefined(retrievedModel)){
              retrievedModel.set('selected', !model.get('selected'));
          }
        });
        lcpSelectionUpdateCollection.reset(null);
        view.$el.find('td').removeClass('edited');
        appModel.clearPagePending();
      });
      
      updateSelectedLcpButton.click(function(e){
        e.preventDefault();
        console.log('selection updates', lcpSelectionUpdateCollection);
        if (lcpSelectionUpdateCollection.isEmpty()) {
          appModel.error('No changes to save');
          return;
        }

        function processClick(formValues){
          console.log('form values', formValues);
          var headers = {};
          headers[appModel.HEADER_APILOG_COMMENT] = formValues['comments'];
          lcpSelectionUpdateCollection.url = function(){
            var url = view.collection.url;
            // NOTE: the comment dialog implies override is confirmed
            url += '?' + appModel.API_PARAM_OVERRIDE + '=true';
            
            if (formValues['set_deselected_to_zero'] === true){
              url += '&' + appModel.API_PARAM_SET_DESELECTED_TO_ZERO +'=true'
            }
            return url;
          };
          
          lcpSelectionUpdateCollection.sync(
            'patch', lcpSelectionUpdateCollection, { headers: headers})
            .done(function(data, textStatus, jqXHR){
              appModel.showConnectionResult(data, {
                title: 'Lab Cherry Pick copy selection updates'
              });
              setSelectedLcpButton.hide();
              lcpSelectionUpdateCollection.reset(null); // clear

              // On success, clear all the buttons
              var searchHash = _.clone(view.listModel.get('search'));
              delete searchHash['show_copy_wells'];
              delete searchHash['show_available_and_retired_copy_wells'];
              delete searchHash['show_unfulfilled'];
              delete searchHash['show_insufficient'];
              searchHash['show_manual'] = 'true';
              showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showUnfulfilledWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showInsufficientWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showManuallySelectedWellsControl.find('input[type="checkbox"]').prop('checked',true);
              view.listModel.set('search',searchHash);
              if (_.isEqual(originalSearchHash,searchHash)){
                view.collection.fetch({ reset: true });
              }
    
            }).fail(function(jqXHR, textStatus, errorThrown){
              console.log('fail', arguments);
              appModel.jqXHRfail.apply(this,arguments); 
            });
        };

        var options = {
          okText: 'Ok',
          cancelText: 'Cancel and return to page',
          title: 'Update Selections?'
        };
        var form_template = appModel._form_template;
        var FormFields = Backbone.Model.extend({
          schema: {
            comments: {
              title: 'Comments',
              type: 'TextArea',
              editorClass: 'input-full',
              validators: ['required'], 
              template: appModel._field_template
            },
            set_deselected_to_zero: {
              title: 'Set Deselected Copy-Well volumes to zero',
              help: 'If selected, deselected Copy Well volumes will be set ' +
                'to zero as part of this operation',
              type: 'Checkbox',
              template: appModel._alt_checkbox_template
            }
          }
        });
        var formFields = new FormFields();
        var form = new Backbone.Form({
          model: formFields,
          template: appModel._form_template
        });
        var _form_el = form.render().el;
        options.view = _form_el;
        options.ok = function(e){
          appModel.clearPagePending();
          var errors = form.commit();
          if(!_.isEmpty(errors)){
            console.log('form errors, abort submit: ' + JSON.stringify(errors));
            return false;
          }else{
            processClick(form.getValue());            
          }
        };        
        appModel.showModal(options);
      });      

      deleteLcpsButton.click(function(e){
        e.preventDefault();
        function processClick(formValues){
          var delete_lab_cherry_picks_url = [
            self.model.resource.apiUri,self.model.key, 
            'delete_lab_cherry_picks'].join('/');
          var headers = {}; // can be used to send a comment
          headers[appModel.HEADER_APILOG_COMMENT] = formValues['comments'];
          $.ajax({
            url: delete_lab_cherry_picks_url,     
            cache: false,
            contentType: 'application/json', 
            dataType: 'json', // what is expected back from the server
            type: 'POST',
            headers: headers
          }).done(function(data, textStatus, jqXHR){
            appModel.showConnectionResult(data, {
              title: 'Delete Lab Cherry Picks'
            });
            self.model.fetch({ reset: true }).done(function(){
              self.uriStack = ['screenercherrypicks'];
              // Remove the child view before calling render, to prevent
              // it from being rendered twice, and calling afterRender twice
              self.removeView('#tab_container');
              self.render();
            });
          }).fail(function(jqXHR, textStatus, errorThrown){
            appModel.jqXHRfail.apply(this,arguments); 
          });
        };
        var msg = 'Delete Lab Cherry Pick selections (return to Screener Cherry Pick View)';
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: function(){
              appModel.showOkCommentForm({
                ok: processClick,
                title: msg
              });
            }
          });
        }else{
          appModel.showOkCommentForm({
            ok: processClick,
            title: msg
          })
        }
      });
      
      reserveAndMapSelectedButton.click(function(e){
        
        e.preventDefault();
        console.log('submit selections for plating');
        var plate_lab_cherrypicks_url = [
          self.model.resource.apiUri,self.model.key, 
          'reserve_map_lab_cherry_picks'].join('/');

        function processClick(overrideInsufficient, comments){
          
          var data = new FormData();
          // send API_PARAM_VOLUME_OVERRIDE = 'volume_override'
          data.append(appModel.API_PARAM_VOLUME_OVERRIDE, overrideInsufficient);
          var headers = {}; // can be used to send a comment
          headers[appModel.HEADER_APILOG_COMMENT] = comments;
          $.ajax({
            url: plate_lab_cherrypicks_url,     
            cache: false,
            contentType: false,
            processData: false,
            dataType: 'json', // what is expected back from the server
            data: data,
            type: 'POST',
            headers: headers
          }).done(function(data, textStatus, jqXHR){
            appModel.showConnectionResult(data, {
              title: 'Lab Cherry Pick Plating result'
            });
            self.model.fetch({ reset: true }).done(function(){
              self.uriStack = ['labcherrypicks'];
              // Remove the child view before calling render, to prevent
              // it from being rendered twice, and calling afterRender twice
              self.removeView('#tab_container');
              self.render();
            });
          }).fail(function(jqXHR, textStatus, errorThrown){
            console.log('errors', arguments);
            
            var jsonError = _.result(jqXHR, 'responseJSON');
            if (!_.isUndefined(jsonError)){
              var error = _.result(jsonError, 'errors');
              var errorFlag = _.result(error,appModel.API_PARAM_VOLUME_OVERRIDE);
              var errorWells = _.result(error, appModel.API_MSG_LCPS_INSUFFICIENT_VOLUME); 
              if (!_.isUndefined(errorFlag )){
                appModel.showOkCommentForm({
                  title: 'Some Copy Wells have insufficient volume, Confirm override?',
                  body: 'Copy Wells: \n' + errorWells.join(', '),
                  okText: 'Override',
                  ok: function(formValues) {
                    var overrideInsufficient = true;
                    processClick(overrideInsufficient, formValues['comments']);
                  }
                });
              } else {
                appModel.jqXHRfail.apply(this,arguments); 
              }
            } else {
              appModel.jqXHRfail.apply(this,arguments); 
            }
          });
        };
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });

      cancelReservation.click(function(e){
        e.preventDefault();
        
        function processClick(formValues){
          var cancel_reservation_url = [
            self.model.resource.apiUri,self.model.key, 
            'cancel_reservation'].join('/');
          var headers = {}; 
          headers[appModel.HEADER_APILOG_COMMENT] = formValues['comments'];
          $.ajax({
            url: cancel_reservation_url,     
            cache: false,
            contentType: 'application/json', 
            dataType: 'json', // what is expected back from the server
            type: 'POST',
            headers: headers
          }).done(function(data, textStatus, jqXHR){
            appModel.showConnectionResult(data, {
              title: 'Cancel Reservation and Remove Plating'
            });
            self.model.fetch({ reset: true }).done(function(){
              self.uriStack = ['labcherrypicks'];
              // Remove the child view before calling render, to prevent
              // it from being rendered twice, and calling afterRender twice
              self.removeView('#tab_container');
              self.render();
            });
          }).fail(function(jqXHR, textStatus, errorThrown){
            appModel.jqXHRfail.apply(this,arguments); 
          });
        };
        var cancelMsg = 'Cancel reserved, allocated copy well volumes and ' +
          'remove plating assignments?';
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: function(){
              appModel.showOkCommentForm({
                ok: processClick,
                title: cancelMsg,
                cancelText: 'No',
                okText: 'Yes'
              });
            }
          });
        }else{
          appModel.showOkCommentForm({
            ok: processClick,
            title: cancelMsg,
            cancelText: 'No',
            okText: 'Yes'
          });
        }
      });
      
      var extra_columns_for_selection = [
        'selected', 'source_copy_well_volume','volume_approved',
        'source_copy_usage_type','source_plate_status',
        'source_plate_date_retired', 'source_plate_screening_count'];
      showUnfulfilledWellsControl.click(function(e){
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(extra_columns_for_selection,includes);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['show_unfulfilled'] = 'true';
            view.listModel.set('search',searchHash);
          } else {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes, extra_columns_for_selection);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            delete searchHash['show_unfulfilled'];
            view.listModel.set('search',searchHash);
          }
          view.$('th').removeClass('selected');
          view.$('tr').removeClass('selected');
        };
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });
      showInsufficientWellsControl.click(function(e){
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(extra_columns_for_selection,includes);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['show_insufficient'] = 'true';
            view.listModel.set('search',searchHash);
          } else {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes, extra_columns_for_selection);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            delete searchHash['show_insufficient'];
            view.listModel.set('search',searchHash);
          }
          view.$('th').removeClass('selected');
          view.$('tr').removeClass('selected');
        };
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });
      showManuallySelectedWellsControl.click(function(e){
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(extra_columns_for_selection,includes);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['show_manual'] = 'true';
            view.listModel.set('search',searchHash);
          } else {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes, extra_columns_for_selection);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            delete searchHash['show_manual'];
            view.listModel.set('search',searchHash);
          }
          view.$('th').removeClass('selected');
          view.$('tr').removeClass('selected');
          
        };
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });
     showCopyWellsControl.click(function(e) {
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(extra_columns_for_selection,includes);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['show_copy_wells'] = 'true';

            showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
            delete searchHash['show_available_and_retired_copy_wells'];
            view.listModel.set('search',searchHash);
            
//            if (self.model.get('number_plates') != 0){
//              appModel.showModalMessage({
//                title: 'Note:',
//                body: 'Selections can not be changed unless plating '+
//                  'assignments are deallocated'
//              });
//            }
          } else {
            // make sure unset
            var searchHash = _.clone(view.listModel.get('search'));
            if (!_.has(searchHash,'show_copy_wells')) {
              return;
            }
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes, extra_columns_for_selection);
            view.listModel.set({ includes: includes}, {reset: false});
            if (_.has(searchHash,'show_copy_wells')) {
              delete searchHash['show_copy_wells'];
              view.listModel.set('search',searchHash);
            }
          }
          // See note above about removing the 'selected' backgrid th class
          view.$('th').removeClass('selected');
          view.$('tr').removeClass('selected');
          view.collection.trigger('show_copy_wells');
        };
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });

      showAllCopyWellsControl.click(function(e) {
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(extra_columns_for_selection,includes);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['show_available_and_retired_copy_wells'] = 'true';
            
            showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
            delete searchHash['show_copy_wells'];
            view.listModel.set('search',searchHash);
//            if (self.model.get('number_plates') != 0){
//              appModel.showModalMessage({
//                title: 'Note:',
//                body: 'Selections can not be changed unless plating '+
//                  'assignments are deallocated'
//              });
//            }
          } else {
            var searchHash = _.clone(view.listModel.get('search'));
            if (!_.has(searchHash,'show_available_and_retired_copy_wells')) {
              return;
            }
            // Make sure unset
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes,extra_columns_for_selection);
            view.listModel.set({ includes: includes}, {reset: false});
            if (_.has(searchHash,'show_available_and_retired_copy_wells')) {
              delete searchHash['show_available_and_retired_copy_wells'];
              view.listModel.set('search',searchHash);
            }
          }
          // See note above about removing the 'selected' backgrid th class
          view.$('th').removeClass('selected');
          view.$('tr').removeClass('selected');
          view.collection.trigger('show_available_and_retired_copy_wells');
          
        };
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });

      return view;
    }, // createLcpView
    
    
    /**
     * Screener Cherry Picks view
     */
    setScreenerCherryPicks: function(delegateStack) {
      var self = this;
      
      if (self.model.get('screener_cherry_pick_count') == 0) {
        self.showScpSearchForm();
      } else { // Screener cherry picks have been entered
        var schemaUrl = [
           self.model.resource.apiUri,self.model.key,'screener_cherry_pick',
           'schema'].join('/');
        appModel.getResourceFromUrl(schemaUrl, function(schemaResult){
          self.createScpView(schemaResult, delegateStack);
        });
      }
    },
    
    showScpSearchForm: function(){
      var self = this;
      // create the search form
      var TextArea2 = Backbone.Form.editors.TextArea.extend({
        render: function() {
          TextArea2.__super__.render.apply(this,arguments);
          this.$el.attr('placeholder', this.schema.placeholder);
          return this;
        },        
      });
      var formSchema = {};
      formSchema['well_search'] = {
        title: 'Well Search',
        key: 'well_search',
        editorClass: 'input-full form-control',
        validators: ['required'],
        type: TextArea2,
        template: appModel._field_template
      };
      
      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs) {
          var errs = {};
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: appModel._form_template
      });
      
      var View = Backbone.Layout.extend({
        template: _.template(genericLayout),
        afterRender: function(){
          $('#resource_content').html(form.render().el);
          form.$el.append([
            '<button type="submit" class="btn btn-default btn-xs" ',
            'style="width: 3em;">ok</input>',
          ].join(''));
          
          var helpLink = $('<a >&nbsp;?</a>');
          form.$el.find('[key="form-group-well_search"]').find('label').append(helpLink);
          helpLink.click(function(e){
            e.preventDefault();
            var bodyMessage = [
              'By Well ID, e.g.:',
              '02091:B15',
              '2091:B15',
              '',
              'Enter multiple wells, separated by newlines, or by commas:',
              '02091:B15',
              '02091:B16',
              '02091:L19',
              'or',
              '02091:B15, 02091:B16, 02091:L19',
              '',
              'By Plate Number, followed by Well Name(s), e.g.:',
              '2091 B15',
              '2091 B16',
              '2091 L19',
              'or',
              '2091 B15 B16 L19',
              'or',
              '2091 B15,B16,L19',
              'or',
              '2091:B15 B16 L19',
              'or',
              '2091:B15,B16,L19'
            ];
            appModel.showModalMessage({
              title: 'Search for wells using patterns',
              body: bodyMessage.join('<br/>')
            });
          });
          form.$el.find('[ type="submit" ]').click(function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }); 
            if(!_.isEmpty(errors)){
              form.$el.find('#well_search').addClass(self.errorClass);
              return;
            }else{
              form.$el.find('#well_search').removeClass(self.errorClass);
            }
            
            var search_value = form.getValue('well_search');
            console.log('search_value', search_value);
            search_value = search_value.split(/\n+/);
            console.log('search_value', search_value);
            
            
            var url = [self.model.resource.apiUri,self.model.key].join('/');
            
            function submit(override){
              if (! _.isUndefined(override)){
                url += '?' + appModel.API_PARAM_OVERRIDE + '=true';
              }
              data = { 'screener_cherry_picks': search_value };
              var headers = {}; // can be used to send a comment
              $.ajax({
                url: url,     
                cache: false,
                contentType: 'application/json', 
                processData: false,
                dataType: 'json', // what is expected back from the server
                data: JSON.stringify(data),
                type: 'PATCH',
                headers: headers
              }).done(function(data, textStatus, jqXHR){
                console.log('success', data);
                appModel.showConnectionResult(data, {
                  title: 'Create Screener Cherry Picks'
                });
                self.model.fetch({ reset: true }).done(function(){
                  self.uriStack = ['screenercherrypicks'];
                  // remove the child view before calling render, to prevent
                  // it from being rendered twice, and calling afterRender twice
                  self.removeView('#tab_container');
                  self.render();
                });
              }).fail(function(jqXHR, textStatus, errorThrown) { 
                console.log('errors', arguments);
                
                var jsonError = _.result(jqXHR, 'responseJSON');
                if (!_.isUndefined(jsonError)){
                  var error = _.result(jsonError, 'errors');
                  var overrideFlag = _.result(error,appModel.API_PARAM_OVERRIDE);
                  var errorMsg = _.result(error, 'screener_cherry_picks'); 
                  var librariesToOverride = _.result(error, 'Libraries');
                  if (!_.isUndefined(librariesToOverride)){
                    errorMsg += '<br/>Libraries:<br/>';
                    errorMsg += librariesToOverride.join('<br/>');
                  }
                  if (!_.isUndefined(overrideFlag)){
                    appModel.showModal({
                      title: 'Override Required',
                      body: errorMsg,
                      okText: 'Confirm Override',
                      ok: function() {
                        var override = true;
                        submit(override);
                      },
                      cancel: function() {
                        // nop
                      }
                    });
                  } else {
                    appModel.jqXHRfail.apply(this,arguments); 
                  }
                } else {
                  appModel.jqXHRfail.apply(this,arguments); 
                }
              });
            };
            
            submit();
          });
        }
      });
      var view = new View();
      self.setView('#tab_container', view).render();
//      self.reportUriStack([]);
    }, // end showScpSearchForm
    
    createScpView: function(schemaResult, delegateStack){
      var self = this;
      
      var url = [self.model.resource.apiUri,self.model.key,
                 'screener_cherry_pick'].join('/');
      
      var extraControls = [];
      var checkboxDiv = $([
          '<div id="show_input_group" class="input-group"></div>'
        ].join(''));
      var showOtherReagentsControl = $([
          '<label class="checkbox-inline" ',
          ' style="margin-left: 10px;" ',
          ' title="Show other reagents with the same vendor ID" >',
          '  <input type="checkbox">Other Reagents',
          '</label>'
        ].join(''));
      checkboxDiv.append(showOtherReagentsControl);
      var showAlternateSelectionsControl = $([
          '<label class="checkbox-inline" ',
          ' style="margin-left: 10px;" ',
          ' title="Show selections of alternate reagents with the same vendor ID as the searched well IDs" >',
          '  <input type="checkbox">Manually Selected',
          '</label>'
        ].join(''));
      checkboxDiv.append(showAlternateSelectionsControl);
      checkboxDiv.prepend('<label for="show_input_group">show</label>');
      extraControls.push(checkboxDiv);

      var deleteScpsButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="deleteScpsButton" href="#">',
          'Delete Screener Cherry Picks</a>'
        ].join(''));
      extraControls.push(deleteScpsButton);

      var setLcpsButtonTitle = 'Set Lab Cherry Picks';
      if(self.model.get('has_pool_screener_cherry_picks') === true){
        setLcpsButtonTitle = 'Set Pool Lab Cherry Picks';
      }
      var setLcpsButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="setLcpsButton" href="#">',
          setLcpsButtonTitle + '</a>'
        ].join(''));
      extraControls.push(setLcpsButton);
      var setDuplexLcpsButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="setDuplexLcpsButton" href="#">',
          'Set Duplex Lab Cherry Picks</a>'
        ].join(''));
      extraControls.push(setDuplexLcpsButton);
      // Set up the grid to record edits of the "selected" column
      var setSelectedButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="save_button_selected" href="#">',
          'Save Selections</a>'
        ].join(''));
      extraControls.push(setSelectedButton);
      var cancelSelectedButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'style="display: none; " ',
          'role="button" id="cancel_selected_button" href="#">',
          'Cancel Selections</a>'
        ].join(''));
      extraControls.push(cancelSelectedButton);
      
      if(self.model.get('total_number_lcps') == 0){
        deleteScpsButton.show();
        if(self.model.get('has_pool_screener_cherry_picks') === true){
          //schemaResult.fields['pool_reagent_vendor_id']['visibility'] = ['l','d'];
          setDuplexLcpsButton.show();
          setLcpsButton.show();
        } else {
          setLcpsButton.show();
        }
      } else {
      }
      if(self.model.get('has_alternate_screener_cherry_pick_selections') === true){
        schemaResult['fields']['searched_well_id']['visibility'] = ['l','d'];
//        if(self.model.get('total_number_lcps') == 0){
//          schemaResult['fields']['selected']['visibility'] = ['l','d'];
//        }
        showAlternateSelectionsControl.show();
      } else {
        showAlternateSelectionsControl.hide();
      }

      // Track user SCP selections
      var Collection = Backbone.Collection.extend({
        // explicitly define the id so that collection compare & equals work
        modelId: function(attrs) {
          return Iccbl.getIdFromIdAttribute( attrs, schemaResult);
        }
      })
      var selectionUpdateCollection = new Collection();

      // Custom SCP Selection Row with styling
      var SelectedScpRow = Backgrid.Row.extend({
        initialize: function () {
          var self = this;
          SelectedScpRow.__super__.initialize.apply(this, arguments);
          
          this.listenTo(this.model.collection, appModel.API_PARAM_SHOW_OTHER_REAGENTS, function(){
            console.log('event:', arguments);
            self._setStyle();
          });
          this.listenTo(this.model.collection, appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS, function(){
            console.log('event:', arguments);
            self._setStyle();
          });
        },
        
        _setStyle: function(){
          var searchHash = this.model.collection.listModel.get('search');
          if( _.result(searchHash,appModel.API_PARAM_SHOW_OTHER_REAGENTS)=='true'
              || _.result(searchHash,appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS)=='true')
          {
            if (this.model.get('searched_well_id')
                ==this.model.get('screened_well_id')) {
              this.$el.closest('tr').addClass('selector_row_group');
            }
          } else {
            this.$el.closest('tr').removeClass('selector_row_group');
          }
          if (_.result(searchHash,appModel.API_PARAM_SHOW_OTHER_REAGENTS)=='true'
            || _.result(searchHash,appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS)=='true' ){
            if (this.model.has('searched_well_id')) {
              if (this.model.get('selected')) {
                this.$el.addClass('selected');
              }
            } 
          } else {
            this.$el.removeClass('selected');
          }
        },
        
        render: function() {
          SelectedScpRow.__super__.render.apply(this, arguments);
          this._setStyle();
          return this;
        }
      });
      // Set SCP List View
      var ListViewSelect = ListView.extend({
        afterRender: function(){
          ListViewSelect.__super__.afterRender.apply(this, arguments);
          // Backgrid specific: all column headers are given a class for their
          // column name: in this case the "selected" class has other meaning, 
          // so remove it
          this.$('th').removeClass('selected');
          this.grid.$el.addClass('selector_table');
          
          return this;
        }
      });
      
      var view = new ListViewSelect({ 
        uriStack: _.clone(delegateStack),
        schemaResult: schemaResult,
        resource: schemaResult,
        url: url,
        extraControls: extraControls,
        row: SelectedScpRow
      });
      Backbone.Layout.setupView(view);
      self.reportUriStack(['screenercherrypicks']);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Screener Cherry Picks for : ' + self.model.key + '</H4>');
      });
      
      selectionUpdateCollection.on('update reset', function(){
        if (!selectionUpdateCollection.isEmpty()){
          appModel.setPagePending(null, 'Selection updates are unsaved, continue anyway?');
          cancelSelectedButton.show();
          setSelectedButton.show();
          deleteScpsButton.hide();
          setLcpsButton.hide();
          setDuplexLcpsButton.hide();
        } else {
          cancelSelectedButton.hide();
          setSelectedButton.hide();

          if(self.model.get('total_number_lcps') == 0){
            deleteScpsButton.show();
            if(self.model.get('has_pool_screener_cherry_picks') === true){
              setDuplexLcpsButton.show();
            } else {
              setLcpsButton.show();
            }
          }
          view.$el.find('td').removeClass('edited');
          appModel.clearPagePending();
        }
      });

      cancelSelectedButton.click(function(e){
        e.preventDefault();
        console.log('cancel selection updates', selectionUpdateCollection);
        selectionUpdateCollection.each(function(model){
          var retrievedModel = view.collection.get(model);
          if (!_.isUndefined(retrievedModel)){
              retrievedModel.set('selected', !model.get('selected'));
          }
        });
        selectionUpdateCollection.reset(null);
        view.$el.find('td').removeClass('edited');
        appModel.clearPagePending();
      });

      deleteScpsButton.click(function(e) {
        e.preventDefault();
        var modelForSubmit = new Backbone.Model(
          self.model.pick(['cherry_pick_request_id'])
        );
        modelForSubmit.set('id', self.model.get('id'));
        modelForSubmit.resource = self.model.resource;
        modelForSubmit.urlRoot = self.model.resource.apiUri;
        modelForSubmit.set('screener_cherry_picks', []);
        modelForSubmit.save(null, {
          patch: true
        }).done(function(data, textStatus, jqXHR) {
          appModel.showConnectionResult(data, {
            title: 'Delete Screener Cherry Picks'
          });
          self.model.fetch({ reset: true }).done(function(){
            self.uriStack = ['screenercherrypicks'];
            // remove the child view before calling render, to prevent
            // it from being rendered twice, and calling afterRender twice
            self.removeView('#tab_container');
            self.render();
          });
        }).fail(function(jqXHR, textStatus, errorThrown) { 
          appModel.jqXHRfail.apply(this,arguments); 
        });
      });

      function processLabCherryPicks(set_lab_cherry_picks_url){
        var headers = {}; // can be used to send a comment
        $.ajax({
          url: set_lab_cherry_picks_url,     
          cache: false,
          contentType: 'application/json', 
          dataType: 'json', // what is expected back from the server
          type: 'POST',
          headers: headers
        }).done(function(data, textStatus, jqXHR){
          console.log('submitted', arguments);
          appModel.showConnectionResult(data, {
            title: 'Set Lab Cherry Picks'
          });
          self.model.fetch({ reset: true }).done(function(){
            self.uriStack = ['labcherrypicks'];
            // remove the child view before calling render, to prevent
            // it from being rendered twice, and calling afterRender twice
            self.removeView('#tab_container');
            self.render();
          });
        }).fail(function(jqXHR, textStatus, errorThrown){
          appModel.jqXHRfail.apply(this,arguments); 
        });
      };
      // submit as Lab Cherry Picks
      setLcpsButton.click(function(e){
        e.preventDefault();
        var set_lab_cherry_picks_url = [
          self.model.resource.apiUri,self.model.key, 'set_lab_cherry_picks'].join('/');
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: function(){
              processLabCherryPicks(set_lab_cherry_picks_url);
            }
          });
        }else{
          processLabCherryPicks(set_lab_cherry_picks_url);
        }
      });

      // submit as Lab Cherry Picks
      setDuplexLcpsButton.click(function(e){
        e.preventDefault();
        var set_lab_cherry_picks_url = [
          self.model.resource.apiUri,self.model.key, 
          'set_duplex_lab_cherry_picks'].join('/');
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: function(){
              processLabCherryPicks(set_lab_cherry_picks_url);
            }
          });
        }else{
          processLabCherryPicks(set_lab_cherry_picks_url);
        }
      });

      // manage SCP selections if lcps not set
      if(self.model.get('total_number_lcps') == 0){
      
        // Manage SCP selection updates
        view.collection.on('add', function(model){
          // cache the 'selected' property for update management
          model.set(
            { selected_on_server: model.get('selected') },
            { silent: true }
          );
        });
        view.collection.on('change', function(model){
          if (!model.has('selected_on_server')){
            // block changes caused by changing fields or other actions that 
            // happen before the selected_on_server flag is set
            return;
          }
          var current_well_search = model.get('searched_well_id');
          var current_well = model.get('screened_well_id');
          
          if (model.get('selected') != model.get('selected_on_server')) {
            selectionUpdateCollection.add(model);
          } else {
            selectionUpdateCollection.remove(model);
          }
          if (model.get('selected')) {
            // unselect others in the group
            _.each(view.collection.where({'searched_well_id':current_well_search}),
              function(model){
                if (model.get('screened_well_id')!==current_well){
                  model.set('selected', false);
                }
            });
          }
          if (!selectionUpdateCollection.isEmpty()){
            appModel.setPagePending(null, 'Selection updates are unsaved, continue anyway?');
            setSelectedButton.show();
          } else {
            setSelectedButton.hide();
          }
        });
        
        // Make sure that on reset actions (page changes), that selections are persisted
        view.collection.on('reset', function(){
          // Note: on "reset" the "add" methods aren't being called
          view.collection.each(function(model){
            model.set(
              { selected_on_server: model.get('selected') },
              { silent: true }
            );
          });
          selectionUpdateCollection.each(function(model){
            // TODO: if the id_attribute is set properly, can use collection.get(model)
            var retrievedModel = view.collection.get(model);
            if (!_.isUndefined(retrievedModel)){
                retrievedModel.set('selected', model.get('selected'));
            }
         });
        });
      } // manage selections if lcps not set
      
      setSelectedButton.click(function(e){
        e.preventDefault();
        console.log('selection updates', selectionUpdateCollection);
        if (selectionUpdateCollection.isEmpty()) {
          appModel.error('No changes to save');
          return;
        }
        function processClick(formValues){
          console.log('form values', formValues);
          var headers = {};
          headers[appModel.HEADER_APILOG_COMMENT] = formValues['comments'];
          selectionUpdateCollection.url = view.collection.url;
          selectionUpdateCollection.sync(
            'patch', selectionUpdateCollection, { headers: headers })
            .done(function(data, textStatus, jqXHR){
              appModel.showConnectionResult(data, {
                title: 'Screener Cherry Pick selection updates'
              });
              setSelectedButton.hide();
              selectionUpdateCollection.reset(null); // clear

              var originalSearchHash = _.clone(view.listModel.get('search'));
              var searchHash = _.clone(view.listModel.get('search'));
              delete searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS];
              searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS] = 'true';
              showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',false);
              showAlternateSelectionsControl.find('input[type="checkbox"]').prop('checked',true);

              view.listModel.set('search',searchHash);
              if (_.isEqual(originalSearchHash,searchHash)){
                view.collection.fetch({ reset: true });
              }
  
            }).fail(function(jqXHR, textStatus, errorThrown){
              appModel.jqXHRfail.apply(this,arguments); 
            });
        };
        appModel.showOkCommentForm({
          ok: processClick,
          title: 'Save Screener Cherry Pick selections?'
        });
        
      });
      
      var extra_colums = [
        'selected', 'searched_well_id','mg_ml_concentration','molar_concentration'];
      showAlternateSelectionsControl.click(function(e) {
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(extra_colums,includes);
            if(self.model.get('total_number_lcps') != 0){
              appModel.showModalMessage({
                title: 'Note:',
                body: 'Selections can not be changed after Lab Cherry Picks have been created'
              });
              includes = _.without(includes, 'selected');
            }
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS] = 'true';

            showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',false);
            delete searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS];
            
            view.listModel.set('search',searchHash);
          } else {
            // make sure unset
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes,extra_colums);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            if (_.has(searchHash,appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS)) {
              delete searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS];
              view.listModel.set('search',searchHash);
            }
          }
          // see note above about removing the 'selected' backgrid th class
          view.$('th').removeClass('selected');
          view.$('tr').removeClass('selected');
          view.collection.trigger(appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS);
        }
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });
      
      // Display or hide other reagents
      showOtherReagentsControl.click(function(e) {
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(extra_colums,includes);
            if(self.model.get('total_number_lcps') != 0){
              appModel.showModalMessage({
                title: 'Note:',
                body: 'Selections can not be changed after Lab Cherry Picks have been created'
              });
              includes = _.without(includes, 'selected');
            }
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS] = 'true';

            showAlternateSelectionsControl.find('input[type="checkbox"]').prop('checked',false);
            delete searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS];
            view.listModel.set('search',searchHash);
            
            
            view.listModel.set('search',searchHash);
          } else {
            // make sure unset
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes,extra_colums);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            if (_.has(searchHash,appModel.API_PARAM_SHOW_OTHER_REAGENTS)) {
              delete searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS];
              view.listModel.set('search',searchHash);
            }
          }
          // see note above about removing the 'selected' backgrid th class
          view.$('th').removeClass('selected');
          view.$('tr').removeClass('selected');
          view.collection.trigger(appModel.API_PARAM_SHOW_OTHER_REAGENTS);
          
        };
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            ok: processClick
          });
        }else{
          processClick();
        }
      });
      var initialSearchHash = view.listModel.get('search');
      if (_.has(initialSearchHash, appModel.API_PARAM_SHOW_OTHER_REAGENTS)
          && initialSearchHash.show_other_reagents.toLowerCase()=='true') {
        showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',true);
      }
      if (_.has(initialSearchHash, appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS)
          && initialSearchHash.show_alternates.toLowerCase()=='true') {
        showAlternateSelectionsControl.find('input[type="checkbox"]').prop('checked',true);
        showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',false);
      }
    }, // end createScpView
    
    /** Backbone.layoutmanager callback **/
    cleanup: function(){
      console.log('cleanup called...');
      this.model = null;
      this.screen = null;
      this.args = null;
      this.tabbed_resources = null;
    }    
  
  });
  
  return CherryPickView;
});