define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/screen/rawDataTransformer',
  'views/generic_detail_layout', 
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'views/list2',
  'utils/tabbedController',
  'utils/wellSelector',
  'templates/genericResource.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, 
            RawDataTransformer, DetailLayout, DetailView, EditView,
            ListView, TabbedController, WellSelector, genericLayout) {
  
  var CherryPickView = TabbedController.extend({
    
    initialize: function(args) {
      
      console.log('initialize cpr view...');
      
      var self = this;
      self.args = args;
      this._classname = 'CherryPickRequestView';
      
      // FIXME: wells_to_leave_empty is a string; should be a list
      // - This is a workaround to produces proper wrapping in the display
      var wellsToLeaveEmpty = self.model.get('wells_to_leave_empty');
      if (!_.isEmpty(wellsToLeaveEmpty)){
        wellsToLeaveEmpty = wellsToLeaveEmpty.split(',').join(', ');
        self.model.set('wells_to_leave_empty', wellsToLeaveEmpty);
      }

      var tabbed_resources_initial = this.tabbed_resources;
      
      TabbedController.prototype.initialize.apply(this,arguments);

      // add back in tabs viewable for screeners
      var screener_tabs = [
        'screenercherrypicks','labcherrypicks'];
      _.each(screener_tabs, function(tab){
        self.tabbed_resources[tab] = tabbed_resources_initial[tab];
      });
      
      if (! this.model.get('number_plates')>0){
        delete this.tabbed_resources['cherrypickplates'];
      }
      if(self.model.get('total_number_lcps') == 0){
        delete this.tabbed_resources['labcherrypicks'];
        delete this.tabbed_resources['sourceplates'];
      }      
      if(self.model.get('number_plates_completed') == 0){
        delete this.tabbed_resources['transformer'];
      }
      _.bindAll(this, 'createScpView','createLcpView',
                'showScpSearchForm','showWellsToLeaveEmptyDialog');
    },

    tabbed_resources: {
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
      //platemapping: {
      //  description : 'Plate Mapping',
      //  title : 'Plate Mapping',
      //  invoke : 'setPlateMapping',
      //  permission: 'cherrypickrequest'
      //},
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
      transformer: {
        title : 'Transform Raw Data',
        description : 'Convert raw data input matrices into a plate/well matrix',
        invoke : 'setRawDataTransformer',
        permission: 'rawdatatransform'
      },
            
    },      

    getTitle: function() {
      return Iccbl.formatString(
        '<H4 id="title">Cherry Pick Request: <a href="#screen/{screen_facility_id}/' + 
        'cherrypickrequest/{cherry_pick_request_id}" >{cherry_pick_request_id}</a>' +
        '</H4>',
        this.model);
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;

//      this.tabbed_resources = _.extend({},this.cherry_pick_tabbed_resources);
//      if (! this.model.get('number_plates')>0){
//        delete this.tabbed_resources['cherrypickplates'];
//      }
//      if(self.model.get('total_number_lcps') == 0){
//        delete this.tabbed_resources['labcherrypicks'];
//        delete this.tabbed_resources['sourceplates'];
//      }      
//      if(self.model.get('number_plates_completed') == 0){
//        delete this.tabbed_resources['transformer'];
//      }
      
      return {
        'base_url': self.model.resource.key + '/' + self.model.key,
        'tab_resources': this.tabbed_resources
      }      
    },
    
    setWarnings: function() {
      var self = this;
      if (!appModel.hasGroup('readEverythingAdmin')){
        // Warnings are only for admin users
        return;
      }
      if (self.model.isNew()){
        return;
      }
      if (self.warningIsProcessing == true){
        return;
      }
      else {
        self.warningIsProcessing = true;
      }
      
      console.log('setWarnings...')
      var $target_el = $('#tab_container-title');
      var mouseTitle = 'Administrative warnings have been detected for the screener cherry picks';
      var warning_div = $('<span id="warning-div" class=""></span>');
      var loading_div = $('<span id="warning-loading-div" class="loading-ellipsis"></span>');
      var warningsButton = $(
        '<button id="cpr-warnings" type="button" class="btn btn-sm btn-danger" '
        + ' title="' + mouseTitle + '" >'
        +'<span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span>'
        +' Warnings</button>');
      
      $target_el.find('#warning-div').remove();
      $target_el.find('#title').append(warning_div);
      warning_div.append('&nbsp;&nbsp;')
      warning_div.append(loading_div);
      
      var url = [self.model.resource.apiUri,self.model.key,'warnings'].join('/');
      var ModelClass = Backbone.Model.extend({
        url : url,
        defaults : { },
        parse: function(resp, options){
          return _.result(resp, Iccbl.appModel.API_RESULT_META, resp);
        }
      });
      var instance = new ModelClass();
      instance.fetch({
        // prevent triggering global event handler: ajaxStart loading gif
        global: false
      }).done(function(){
        
        warning_div.empty();
        $target_el.find('#title').append(warning_div);
        
        var fields_to_show = [
          'cherry_pick_allowance_warning','restricted_libraries',
          'concentration_warnings','copyplates_not_found','screener_picks_not_screened',
          'duplicate_screener_cherry_picks'
        ];
        var shownFields = instance.pick(fields_to_show);
        if (!_.isEmpty(shownFields)){
          
          var comments = [];
          var body = $('<div class="panel"></div>');
          
          _.each(_.keys(shownFields), function(key){
            var val = shownFields[key];
            if (!val || _.isEmpty(val)) return;
            if (key=='cherry_pick_allowance_warning'){
              comments.push(val);
              body.append(val + '<br>');
            }
            else if (key=='restricted_libraries'){
              // TODO: convert library.screening_status vocab
              comments.push('Restricted Libraries:');
              body.append(comments[comments.length-1]+'<br>');
              comments = comments.concat(val);
              body.append(val.join('<br>'));
              body.append('<br>');
            }
            else if (key=='concentration_warnings'){
              comments.push('The following picked wells are not at the maximal '
                + 'concentration found for the reagent:');
              body.append(comments[comments.length-1]+'<br>');
              _.each(_.keys(val), function(well_id){
                console.log('formatting', well_id, val[well_id]);
                var message = 
                  'Selected Well: {well_id} ({selected_concentration}) '
                  + 'Alternate: {alternate_well} ({alternate_concentration})';
                comments.push(Iccbl.formatString(message, 
                  _.extend({well_id: well_id}, val[well_id])));
                var selectedWellLink = $('<a>', {
                  tabIndex : -1,
                  href : Iccbl.formatString('#library/{library_short_name}/well/{well_id}',
                    { library_short_name: val[well_id]['library_short_name'], well_id: well_id }),
                  target : '_blank',
                }).text(well_id);
                var alternateWellLink = $('<a>', {
                  tabIndex : -1,
                  href : Iccbl.formatString('#library/{library_short_name}/well/{well_id}',
                    { library_short_name: val[well_id]['library_short_name'], 
                      well_id: val[well_id]['alternate_well'] }),
                  target : '_blank',
                }).text(val[well_id]['alternate_well']);
                body.append([
                  'Selected: ', selectedWellLink, ' (', val[well_id]['selected_concentration'] ,'), ', 
                  'Alternate: ', alternateWellLink, ' (', val[well_id]['alternate_concentration'] ,'), '
                ]);
                body.append('<br>');
                
              });
            }
            else if (key=='copyplates_not_found'){
              
              var platesEsearch = _.keys(val).join(appModel.SEARCH_DELIMITER);
              var link  = $('<a>', {
                tabIndex : -1,
                href : ['#librarycopyplate', appModel.API_PARAM_ENCODED_SEARCH,
                  platesEsearch].join('/'),
                target : '_blank',
              }).text('No (available) Cherry Pick Screening copies found for the plates:');
              body.append(link);
              body.append('<br>');
              comments.push('No (available) Cherry Pick Screening plates found:');
              _.each(_.keys(val), function(plate){
                comments.push(plate + ': ' + val[plate].join(', '));
                body.append(comments[comments.length-1] + '<br>');
              });
              
            }
            else if (key=='screener_picks_not_screened'){
              comments.push('Screener Picks that have not been screened:');
              body.append(comments[comments.length-1] + '<br>');
              comments.push(val.join(', '));
              body.append(val.join(', ') + '<br>');
            }
            else if (key=='duplicate_screener_cherry_picks'){
              comments.push(
                'The following wells have already been cherry picked '
                + 'for this Screen, which is against policy:');
              body.append(comments[comments.length-1] + '<br>');
              comments.push(val.join(', '));
              body.append(val.join(', ')+ '<br>');
            }
          });
          
          if (!_.isEmpty(comments)){
            warningsButton.attr({'title': comments.join('\n')});
            warningsButton.click(function(e){
              e.preventDefault();
              Iccbl.appModel.showModalMessage({
                title: 'Administrative warnings for Screener Cherry Picks',
                view: body,
                buttons_on_top: true,
              });
            });
            warning_div.append('&nbsp;&nbsp;')
            warning_div.append(warningsButton);
          }
        }
        self.warningIsProcessing = false;
      }).fail(function(jqXHR, textStatus, errorThrown){
        console.log('fail', arguments);
        appModel.jqXHRfail.apply(this,arguments); 
        $target_el.find('#warning-div').remove();
        self.warningIsProcessing = false;
      });
      
    },
    
    afterRender: function(){
      var self = this;
      TabbedController.prototype.afterRender.apply(this,arguments);
      if (_.isUndefined(this.screen)
          ||  this.screen.get('facility_id') != this.model.get('screen_facility_id') ){
        appModel.getModel('screen', this.model.get('screen_facility_id'), 
          function(screen){
            self.screen = screen;
            self.setWarnings();
          });
      } else {
        self.setWarnings();
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
        
        save_success: function(data, textStatus, jqXHR){
          console.log('success');

          // Override so that the CPR can be displayed in the screen resource:
          var meta = _.result(data, appModel.API_RESULT_META, null);
          
          if (!_.isEmpty(meta)) {
            appModel.showJsonMessages(meta);
          }
          if (_.isUndefined(this.model) || _.isUndefined(this.model.key)){
            if (! this.model.isNew()){
              this.model.key = Iccbl.getIdFromIdAttribute( this.model,this.model.resource );
              appModel.router.navigate([
                self.screen.resource.key,self.screen.key,'cherrypickrequest',
                this.model.key].join('/'), 
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
        },
        
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
        EditView: editView,
        DetailView: detailView
      }));
      view.showEdit = function() {
        var innerself = this;
        appModel.initializeAdminMode(function(){
          var fields = self.model.resource.fields;
          fields['requested_by_id'].choiceHash = 
            appModel._get_screen_member_choices(self.screen);
          
          // TODO: resource/cherrypickrequest/write
          appModel.getAdminUserOptions(function(options){
            fields['volume_approved_by_username'].choiceHash = options;
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
          appModel.showModalError('must define the assay plate type');
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
      var resource = appModel.getResource('librarycopyplate');
      console.log('resource', resource);
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
          if (_.has(fields, hidden_field)){
            fields[hidden_field]['visibility'] = [];
            includes.unshift('-' + hidden_field);
          }
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
    
    setRawDataTransformer: function(delegateStack) {
      console.log('setRawDataTransformer', delegateStack);
      var self = this;
      
      var rawdataResource = appModel.getResource('rawdatatransform');
      var schemaUrl = [rawdataResource.apiUri,
                       'schema'].join('/');
      function showTransformer(schemaResult) {
        var cpr_id = self.model.get('cherry_pick_request_id');
        var options = {
          failCallback: function(){
            var newModel = appModel.newModelFromResource(schemaResult);
            newModel.set('cherry_pick_request_id', cpr_id );
            var plate_numbers = _.times(
              self.model.get('number_plates_completed'), 
              function(i){ return i+1;});
            newModel.set('plate_ranges', plate_numbers.join(','));
            newModel.set('output_filename', 'cpr' + cpr_id + '_' + plate_numbers.join('_'));
            newModel.resource = schemaResult;
            showView(newModel);
          }
        };
        var rdtKey = [self.model.get('screen_facility_id'),cpr_id].join('/');
        appModel.getModel(
          rawdataResource.key, rdtKey, 
          function(model){
            model.resource = schemaResult;
            if (!model.has('library_plate_size')){
              model.set('library_plate_size', schemaResult.fields['library_plate_size'].default );
            }
            showView(model);
          }, 
          options);
      };
      function showView(model){
        var view = new RawDataTransformer({
          cherry_pick_request: self.model,
          model: model,
          uriStack: delegateStack
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        this.consumedStack = ['transformer'];
        
      }
      
      appModel.getResourceFromUrl(schemaUrl, showTransformer);
      
    },    
    setCherryPickPlates: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,
        'cherry_pick_plate'].join('/');
      var resource = appModel.getResource('cherrypickassayplate');
      var selectionCollection = null;

      var ListCollectionClass = Iccbl.MyCollection.extend({
        /** explicitly define the id so that collection compare & equals work **/
        modelId: function(attrs) {
          return Iccbl.getIdFromIdAttribute( attrs, resource);
        },
        url: url,
        /** Override parse to set selected default value **/
        parseRecords: function(resp,options){
          var modelArray = ListCollectionClass.__super__.parseRecords.apply(this, arguments);
          _.each(modelArray, function(model){
            model['selected'] = false;
          });
          return modelArray;
        }
      });
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
      
        var SelectionCollectionClass = Iccbl.CollectionOnClient.extend({
          /** explicitly define the id so that collection compare & equals work **/
          modelId: function(attrs) {
            return Iccbl.getIdFromIdAttribute( attrs, resource);
          },
          url: url,
          /** Override parse to set selected default value **/
          parse: function(resp){
            var modelArray = SelectionCollectionClass.__super__.parse.apply(this, arguments);
            _.each(modelArray, function(model){
              model['selected'] = false;
            });
            return modelArray;
          }
        });
        
        selectionCollection = new SelectionCollectionClass();
        
        // Fetch the entire cpap collection to use for selection tracking
        selectionCollection.fetch({

          data: { limit: 0 },

          success: function(collection, response) {
            
            var extraControls = [];
            var setPlatedDiv = $(
              '<span title="Set the plating activity date for the selected plates"/>'
            );
            var setPlatedButton = $([
              '<a class="btn btn-default btn-sm disabled" ',
              'role="button" id="set_plated_button" ',
              'href="#">',
              'Set plated</a>'
            ].join(''));
            setPlatedDiv.append(setPlatedButton);
            extraControls.push(setPlatedDiv);
            var setScreenedDiv = $(
              '<span title="Set the plating activity date for the selected plates"/>'
            );
            var setScreenedButton = $([
              '<a class="btn btn-default btn-sm disabled" ',
              'role="button" id="set_screened_button" ',
              'title="set the screening activity date for the selected plates" href="#">',
              'Set screened</a>'
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
          '<a class="btn btn-default btn-sm" ',
            'role="button" id="plate_mapping_file" href="#">',
            'Download plate mapping files</a>'
          ].join(''));
        extraControls.push(downloadPlateMappingButton);
        var showPlateMappingButton = $([
          '<a class="btn btn-default btn-sm" ',
             'title="Show the plate mapping in a grid (available if plates are assigned)"',
            'role="button" id="plate_mapping_grid" href="#">',
            'Show plate mapping</a>'
          ].join(''));
        extraControls.push(showPlateMappingButton);

        collection = new ListCollectionClass();
        var view = new ListView({ 
          collection: collection,
          uriStack: _.clone(delegateStack),
          resource: resource,
          url: url,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
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
          appModel.downloadUrl(url);
        });
        showPlateMappingButton.click(function(e){
          e.preventDefault();
          var url = [self.model.resource.apiUri,self.model.key,
            'lab_cherry_pick_plating'].join('/');
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
              screened_by_id: {
                title: 'Screened By',
                key: 'screened_by_id',
                type: EditView.ChosenSelect,
                editorClass: 'chosen-select',
                options: appModel._get_screen_member_choices(self.screen),
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
              
              var SubmitCollection = Backbone.Collection.extend({
                url: selectionCollection.url
              });
              var submitCollection = new SubmitCollection(selectedEntries);
              submitCollection.each(function(model){
                 model.set('screening_date', values['screening_date']);
                 model.set('screened_by_id', values['screened_by_id']);
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
            // FIXME: known bug, if user cancels the plating dialog, and then
            // tries to return to it, it does not render
            
            //$('#modal').one('hidden.bs.modal',  function () {
            //  console.log('hidden.bs.modal', arguments);
            //  appModel.getAdminUserOptions(dialogFunction, 'cherrypickrequest');    
            //});
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
            
            var SubmitCollection = Backbone.Collection.extend({
              url: selectionCollection.url
            });
            var submitCollection = new SubmitCollection(selectedEntries);
            submitCollection.each(function(model){
               model.set('plating_date', values['plating_date']);
               model.set('plated_by_username', values['plated_by_username']);
               model.unset('plated_by_id'); // only use username
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
        var downloadPlateMappingButton = $([
          '<a class="btn btn-default btn-sm" ',
             'title="Download plate mapping file (available if plates are assigned)"',
            'role="button" id="plate_mapping_file" href="#">',
            'Download plate mapping files</a>'
          ].join(''));
        var showPlateMappingButton = $([
          '<a class="btn btn-default btn-sm" ',
             'title="Show the plate mapping in a grid (available if plates are assigned)"',
            'role="button" id="plate_mapping_grid" href="#">',
            'Show plate mapping</a>'
          ].join(''));
        downloadPlateMappingButton.click(function(e){
          e.preventDefault();
          console.log('download plate mapping file');
          var url = [self.model.resource.apiUri,self.model.key,
                     'plate_mapping_file'].join('/');
          appModel.downloadUrl(url);
        });
        
        showPlateMappingButton.click(function(e){
          e.preventDefault();
          self.showPlateMappingGrid(url);
        });
        schemaUrl = url + '/schema';

        if (appModel.hasGroup('readEverythingAdmin')){
          plate_mapping_controls.push(downloadPlateMappingButton);
        }      
        plate_mapping_controls.push(showPlateMappingButton);
      }

      function createView(resource){
        view = self.createLcpView(delegateStack,resource, url);
        self.listenTo(view, 'afterRender', function(event) {
          
          if (!_.isEmpty(plate_mapping_controls)){
            view.$el.find('#list_controls').prepend(plate_mapping_controls);
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
          '<div id="show_input_group" class="input-group pull-down pull-left"></div>'
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
      var setSelectedLcpButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="save_button_lcp_selected" href="#">',
          'Save selections</a>'
        ].join(''));
      var updateSelectedLcpButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="update_selected_button" href="#">',
          'Update selections</a>'
        ].join(''));
      var cancelSelectedButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="cancel_selected_button" href="#">',
          'Cancel selections</a>'
        ].join(''));
      var reserveAndMapSelectedButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="reserve_map_selected_button" href="#">',
          'Reserve selections and map to plates</a>'
        ].join(''));
      var deleteLcpsButton = $([
          '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="deleteLcpsButton" href="#">',
          'Delete lab cherry picks</a>'
        ].join(''));
      var cancelReservation = $([
          '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="cancel_reservation" href="#">',
          'Cancel reservation and delete plating assignments</a>'
        ].join(''));
      var extraListControls = [];
      if (appModel.hasGroup('readEverythingAdmin')){
        checkboxDiv.append(showCopyWellsControl);
        checkboxDiv.append(showAllCopyWellsControl);
        if (self.model.get('number_unfulfilled_lab_cherry_picks') > 0){
          checkboxDiv.append(showUnfulfilledWellsControl);
        }
        checkboxDiv.append(showInsufficientWellsControl);
        checkboxDiv.append(showManuallySelectedWellsControl);
        checkboxDiv.prepend('<label for="show_input_group">show</label>');
        extraListControls.push(checkboxDiv);
        extraControls.push(setSelectedLcpButton);
        extraControls.push(updateSelectedLcpButton);
        extraControls.push(cancelSelectedButton);
        if(appModel.hasPermission('labcherrypick','write')){
          extraControls.push(reserveAndMapSelectedButton);
        }
        if(appModel.hasPermission('labcherrypick','write')){
          extraControls.push(deleteLcpsButton);
        }
        if(appModel.hasPermission('labcherrypick','write')){
          extraControls.push(cancelReservation);
        }
      }
      

      if(appModel.hasPermission('labcherrypick','write')){
        if (self.model.get('number_plates') == 0){
          deleteLcpsButton.show();
          reserveAndMapSelectedButton.show();
        } else {
          cancelReservation.show();
        }
      }
      
      if (appModel.hasGroup('readEverythingAdmin')){

        ///// Library and Plate comments /////
        resource.fields['library_plate']['backgridCellType'] = 
          Iccbl.CommentArrayLinkCell.extend({
            comment_attribute: 'library_plate_comment_array',
            title_function: function(model){
              return 'Comments for Plate: ' 
                + model.get('library_short_name') + '/' 
                + model.get('source_copy_name')  + '/'
                + model.get('library_plate');
            }
          });

        resource.fields['source_copy_name']['backgridCellType'] = 
          Iccbl.LinkCell.extend({
            render: function(){
              var self = this;
              Iccbl.LinkCell.prototype.render.apply(this, arguments);
              var comments = this.model.get('source_copy_comments');
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
          },resource.fields['source_copy_name'].display_options);

        resource.fields['library_name']['backgridCellType'] =
          Iccbl.CommentArrayLinkCell.extend({
            comment_attribute: 'library_comment_array',
            title_function: function(model){
              return 'Comments for library: ' + model.get('library_short_name');
            }
          });
      
        resource.fields['library_short_name']['backgridCellType'] =
          Iccbl.CommentArrayLinkCell.extend({
            comment_attribute: 'library_comment_array',
            title_function: function(model){
              return 'Comments for library: ' + model.get('library_short_name');
            }
          });
      
        ///// end Library and Plate comments /////
        if(self.model.get('has_pool_screener_cherry_picks') === true){
          resource.fields['pool_reagent_vendor_id']['visibility'] = ['l','d'];
        }
      
      } // end Admin controls

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
          var searchHash = this.model.collection.listModel.get(appModel.URI_PATH_SEARCH);
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
        resource: resource,
        url: url,
        row: SelectedLcpRow,
        extraControls: extraControls,
        extraListControls: extraListControls
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Lab Cherry Picks for : ' + self.model.key + '</H4>');
      });
    
      if (appModel.hasGroup('readEverythingAdmin')){
        // setup Admin-functions
        
        var initialSearchHash = view.listModel.get(appModel.URI_PATH_SEARCH);
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
        
        view.grid.columns.on('update', function(){
          view.$el.find('td').removeClass('edited');
        });
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
          //if (!(showCopyWellsControl.find('input[type="checkbox"]').prop('checked')
          //    || showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked'))){
          //  return;
          //}
  
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
            appModel.showModalError('No changes to save');
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
                // On success, clear all the buttons
                var originalSearchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
                var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
                delete searchHash['show_copy_wells'];
                delete searchHash['show_available_and_retired_copy_wells'];
                delete searchHash['show_unfulfilled'];
                delete searchHash['show_insufficient'];
                var includes = _.clone(view.listModel.get('includes'));
                includes = _.without(includes, 'selected');
                view.listModel.set({ includes: includes}, {reset: false});
                showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
                showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
                
                showUnfulfilledWellsControl.find('input[type="checkbox"]').prop('checked',false);
                showInsufficientWellsControl.find('input[type="checkbox"]').prop('checked',false);
  
                showManuallySelectedWellsControl.find('input[type="checkbox"]').prop('checked',true);
                searchHash['show_manual'] = 'true';
                view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
                //              if (_.isEqual(originalSearchHash,searchHash)){
                //                view.collection.fetch({ reset: true });
                //              }
                // Note: this may trigger another a superfluous fetch after the listmodel
                // update. Necessary to finally clear the lcpSelectionCollection
                view.collection.fetch({ reset: true }).done(function(){
                  // have to wait for the fetch operation to reset, due to the 
                  // asynchronous event handling
                  lcpSelectionUpdateCollection.reset(null); // clear
                });
      
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
            appModel.showModalError('No changes to save');
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
  
                // On success, clear all the buttons
                var originalSearchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
                var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
                delete searchHash['show_copy_wells'];
                delete searchHash['show_available_and_retired_copy_wells'];
                delete searchHash['show_unfulfilled'];
                delete searchHash['show_insufficient'];
                var includes = _.clone(view.listModel.get('includes'));
                includes = _.without(includes, 'selected');
                view.listModel.set({ includes: includes}, {reset: false});
                
                showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
                showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
  
                showUnfulfilledWellsControl.find('input[type="checkbox"]').prop('checked',false);
                searchHash['show_manual'] = 'true';
  
                showInsufficientWellsControl.find('input[type="checkbox"]').prop('checked',false);
                showManuallySelectedWellsControl.find('input[type="checkbox"]').prop('checked',true);
                view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
                //if (_.isEqual(originalSearchHash,searchHash)){
                //  view.collection.fetch({ reset: true });
                //}
                view.collection.fetch({ reset: true }).done(function(){
                  // have to wait for the fetch operation to reset, due to the 
                  // asynchronous event handling
                  lcpSelectionUpdateCollection.reset(null); // clear
                });
      
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
                    body: 'Copy Wells: ' + errorWells.join(', '),
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
        
        showUnfulfilledWellsControl.find('input[type="checkbox"]').change(function(e){
          function processClick(){
            if (e.target.checked) {
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              showInsufficientWellsControl.find('input[type="checkbox"]').prop('checked',false);
              showManuallySelectedWellsControl.find('input[type="checkbox"]').prop('checked',false);
              delete searchHash['show_manual'];
              delete searchHash['show_insufficient'];
              searchHash['show_unfulfilled'] = 'true';
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
              
            } else {
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              delete searchHash['show_unfulfilled'];
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
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
        showInsufficientWellsControl.find('input[type="checkbox"]').change(function(e){
          function processClick(){
            if (e.target.checked) {
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              showUnfulfilledWellsControl.find('input[type="checkbox"]').prop('checked',false);
              delete searchHash['show_unfulfilled'];
              searchHash['show_insufficient'] = 'true';
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
            } else {
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              delete searchHash['show_insufficient'];
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
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
        showManuallySelectedWellsControl.find('input[type="checkbox"]').change(function(e){
          function processClick(){
            if (e.target.checked) {
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              showUnfulfilledWellsControl.find('input[type="checkbox"]').prop('checked',false);
              delete searchHash['show_unfulfilled'];
              searchHash['show_manual'] = 'true';
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
            } else {
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              delete searchHash['show_manual'];
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
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
        var extra_columns_for_selection = [
          'selected', 'source_copy_well_volume','volume_approved',
          'source_copy_usage_type','source_plate_status',
          'source_plate_date_retired', 'source_plate_screening_count',
          'source_plate_cp_screening_count'];
       showCopyWellsControl.click(function(e) {
          function processClick(){
            var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
            if (e.target.checked) {
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.union(extra_columns_for_selection,includes);
              view.listModel.set({ includes: includes}, {reset: false});
              searchHash['show_copy_wells'] = 'true';
  
              showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              delete searchHash['show_available_and_retired_copy_wells'];
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
            } else {
              // make sure unset
              if (!_.has(searchHash,'show_copy_wells')) {
                return;
              }
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.difference(includes, extra_columns_for_selection);
              view.listModel.set({ includes: includes}, {reset: false});
              if (_.has(searchHash,'show_copy_wells')) {
                delete searchHash['show_copy_wells'];
                view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
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
  
        showAllCopyWellsControl.find('input[type="checkbox"]').change(function(e) {
          function processClick(){
            if (e.target.checked) {
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.union(extra_columns_for_selection,includes);
              view.listModel.set({ includes: includes}, {reset: false});
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              searchHash['show_available_and_retired_copy_wells'] = 'true';
              
              showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              delete searchHash['show_copy_wells'];
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
            } else {
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              if (!_.has(searchHash,'show_available_and_retired_copy_wells')) {
                return;
              }
              // Make sure unset
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.difference(includes,extra_columns_for_selection);
              view.listModel.set({ includes: includes}, {reset: false});
              if (_.has(searchHash,'show_available_and_retired_copy_wells')) {
                delete searchHash['show_available_and_retired_copy_wells'];
                view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
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
      } // END Admin-functions
      
      
      return view;
    }, // createLcpView
    
    /**
     * Screener Cherry Picks view
     */
    setScreenerCherryPicks: function(delegateStack) {
      var self = this;
      
      if (self.model.get('screener_cherry_pick_count') == 0) {
        self.showScpSearchForm();
        self.reportUriStack([]);
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
      
      var formSchema = {};
      formSchema['well_search'] = {
        title: 'Well Search',
        key: 'well_search',
        editorClass: 'input-full form-control',
        validators: ['required'],
        type: EditView.TextArea2,
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
            'style="width: 3em;">Ok</input>',
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
              form.$el.find('#well_search').addClass("has-error");
              return;
            }else{
              form.$el.find('#well_search').removeClass("has-error");
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
              var data = { 'screener_cherry_picks': search_value };
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
                self.setWarnings();
                
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
      Backbone.Layout.setupView(view);
      self.setView('#tab_container', view).render();
    }, // end showScpSearchForm
    
    createScpView: function(schemaResult, delegateStack){
      var self = this;
      
      var url = [self.model.resource.apiUri,self.model.key,
                 'screener_cherry_pick'].join('/');
      
      var extraControls = [];
      var showOtherReagentsControl = $([
          '<label class="checkbox-inline" ',
          ' style="margin-left: 10px;" ',
          ' title="Show other reagents with the same vendor ID" >',
          '  <input type="checkbox">Other Reagents',
          '</label>'
        ].join(''));
      var showAlternateSelectionsControl = $([
          '<label class="checkbox-inline" ',
          ' style="margin-left: 10px;" ',
          ' title="Show selections of alternate reagents with the same vendor ID as the searched well IDs" >',
          '  <input type="checkbox">Alternate Selections',
          '</label>'
        ].join(''));
      var deleteScpsButton = $([
          '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="deleteScpsButton" href="#">',
          'Delete screener cherry picks</a>'
        ].join(''));
      var setLcpsButtonTitle = 'Set lab cherry picks';
      if(self.model.get('has_pool_screener_cherry_picks') === true){
        setLcpsButtonTitle = 'Set pool lab cherry picks';
      }
      var setLcpsButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="setLcpsButton" href="#">',
          setLcpsButtonTitle + '</a>'
        ].join(''));
      var setDuplexLcpsButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="setDuplexLcpsButton" href="#">',
          'Set duplex lab cherry picks</a>'
        ].join(''));
      // Set up the grid to record edits of the "selected" column
      var setSelectedButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="save_button_selected" href="#">',
          'Save selections</a>'
        ].join(''));
      var cancelSelectedButton = $([
        '<a class="btn btn-default btn-sm" ',
          'style="display: none; " ',
          'role="button" id="cancel_selected_button" href="#">',
          'Cancel selections</a>'
        ].join(''));

      if (appModel.hasGroup('readEverythingAdmin')){
        var checkboxDiv = $([
            '<div id="show_input_group" class="input-group pull-down pull-left"></div>'
          ].join(''));
        checkboxDiv.append(showOtherReagentsControl);
        checkboxDiv.append(showAlternateSelectionsControl);
        checkboxDiv.prepend('<label for="show_input_group">show</label>');
        extraControls.push(checkboxDiv);
        extraControls.push(deleteScpsButton);
        extraControls.push(setLcpsButton);
        extraControls.push(setDuplexLcpsButton);
        extraControls.push(setSelectedButton);
        extraControls.push(cancelSelectedButton);
      }

      
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
          var searchHash = this.model.collection.listModel.get(appModel.URI_PATH_SEARCH);
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
        resource: schemaResult,
        url: url,
        extraControls: extraControls,
        row: SelectedScpRow
      });
      Backbone.Layout.setupView(view);
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
            self.setWarnings();
            
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
          appModel.showConnectionResult(data, {
            title: 'Set Lab Cherry Picks'
          });
          self.model.fetch({ reset: true }).done(function(){
            self.uriStack = ['labcherrypicks'];
            // remove the child view before calling render, to prevent
            // it from being rendered twice, and calling afterRender twice
            self.removeView('#tab_container');
            self.render();
            self.setWarnings();
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
          appModel.showModalError('No changes to save');
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

              var originalSearchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
              delete searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS];
              searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS] = 'true';
              showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',false);
              showAlternateSelectionsControl.show()
              showAlternateSelectionsControl.find('input[type="checkbox"]').prop('checked',true);
              
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
              view.collection.fetch({ reset: true }).done(function(){
                // have to wait for the fetch operation to reset, due to the 
                // asynchronous event handling
                selectionUpdateCollection.reset(null); // clear
              });
              self.setWarnings();
  
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
      showAlternateSelectionsControl.find('input[type="checkbox"]').change(function(e) {
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
            var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
            searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS] = 'true';

            showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',false);
            delete searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS];
            
            view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
          } else {
            // make sure unset
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes,extra_colums);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
            if (_.has(searchHash,appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS)) {
              delete searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS];
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
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
      showOtherReagentsControl.find('input[type="checkbox"]').change(function(e) {
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
            var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
            searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS] = 'true';

            showAlternateSelectionsControl.find('input[type="checkbox"]').prop('checked',false);
            delete searchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS];
            view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
            
            
            view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
          } else {
            // make sure unset
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.difference(includes,extra_colums);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get(appModel.URI_PATH_SEARCH));
            if (_.has(searchHash,appModel.API_PARAM_SHOW_OTHER_REAGENTS)) {
              delete searchHash[appModel.API_PARAM_SHOW_OTHER_REAGENTS];
              view.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
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
      var initialSearchHash = view.listModel.get(appModel.URI_PATH_SEARCH);
      if (_.has(initialSearchHash, appModel.API_PARAM_SHOW_OTHER_REAGENTS)
          && initialSearchHash.show_other_reagents.toLowerCase()=='true') {
        showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',true);
      }
      if (_.has(initialSearchHash, appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS)
          && initialSearchHash[appModel.API_PARAM_SHOW_ALTERNATE_SELECTIONS].toLowerCase()=='true') {
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