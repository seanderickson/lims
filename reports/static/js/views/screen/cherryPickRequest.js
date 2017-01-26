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
      
      this.tabbed_resources = _.extend({},this.cherry_pick_tabbed_resources);
      if (! this.model.get('number_plates')>0){
        delete this.tabbed_resources['cherrypickplates'];
      }
      
      TabbedController.prototype.initialize.apply(this,arguments);
      
      if (_.isUndefined(this.screen)){
        appModel.getModel('screen', this.model.get('screen_facility_id'), 
          function(screen){
            self.screen = screen;
          });
      }
      
      _.bindAll(this, 'createScpView','showScpSearchForm','showWellsToLeaveEmptyDialog');
      
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
          
          EditView.prototype.afterRender.apply(this,arguments);
        }
      });
      var saveSuccessCallBack = function(response){
        // TODO: modify generic_edit to handle custom-nested URLs for edit success URL
        var meta = _.result(response, 'meta', null);
        if (meta) {
          appModel.showJsonMessages(meta);
        }
        appModel.router.navigate([
          self.screen.resource.key,self.screen.key,'cherrypickrequest',self.model.key].join('/'), 
          {trigger:true});
      };
      detailView = DetailView.extend({
        afterRender: function() {
          console.log('afterrender...');
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
          fields['volume_approved_by_username'].choices = 
            appModel.getAdminUserOptions();
          
          DetailLayout.prototype.showEdit.apply(innerself,arguments);
        });  
      };
      this.tabViews[key] = view;

      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
      
    },

    showWellsToLeaveEmptyDialog: function(editForm, currentValue, editable){
      var self = this;
      var rowCollection = new Backbone.Collection();
      if (_.isUndefined(editable)){
        editable = false;
      }
      if (_.isUndefined(currentValue)){
        currentValue = self.model.get('wells_to_leave_empty');
      }
      console.log('current_value', currentValue);
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
    },
    
    setCherryPickPlates: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,
        'cherry_pick_plate'].join('/');
      var resource = appModel.getResource('cherrypickassayplate');
      var extraControls = [];
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: extraControls
      });
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Assay Plates for CPR: ' + self.model.key + '</H4>');
      });
      
      
    },
    
    setLabCherryPicks: function(delegateStack) {
      var self = this;

      function createLcpView(schemaResult){
        var lcp_selection_updates = {};
        var selectionUpdateCollection = new Backbone.Collection();
        var url = [self.model.resource.apiUri,self.model.key,
          'lab_cherry_pick'].join('/');

        if(self.model.get('has_pool_screener_cherry_picks') === true){
          schemaResult.fields['pool_reagent_vendor_id']['visibility'] = ['l','d'];
        }
        
        var extraControls = [];
        var showCopyWellsControl = $([
            '<label class="checkbox-inline" ',
            'title="Show all available Cherry Pick Source Plate copies" >',
            '  <input id="showCopyWells" type="checkbox">show Available Cherry Pick Source Wells',
            '</label>'
          ].join(''));
        extraControls.push(showCopyWellsControl);
        var showAllCopyWellsControl = $([
            '<label class="checkbox-inline" ',
            'title="Show all available and retired Cherry Pick Source and Library Screening copies">',
            '  <input id="showAllCopyWells" type="checkbox">show Retired Library Screening Wells',
            '</label>'
          ].join(''));
        extraControls.push(showAllCopyWellsControl);
        var showUnfulfilledWellsControl = $([
            '<label class="checkbox-inline">',
            '  <input id="showUnfulfilled" type="checkbox">show Unfulfilled',
            '</label>'
          ].join(''));
        extraControls.push(showUnfulfilledWellsControl);

        var setSelectedLcpButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'style="display: none; " ',
            'role="button" id="save_button_lcp_selected" href="#">',
            'Save Selections</a>'
          ].join(''));
        extraControls.push(setSelectedLcpButton);
        var reserveAndMapSelectedButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
//            'style="display: none; " ',
            'role="button" id="reserve_map_selected_button" href="#">',
            'Reserve Selections and Map to Plates</a>'
          ].join(''));
        extraControls.push(reserveAndMapSelectedButton);
//        var overrideWellVolumeControl = $([
//            '<label class="checkbox-inline" ',
//            'title="If chosen wellcopy volumes are insufficient, select to override">',
//            '  <input id="overrideWellVolume" type="checkbox">override Well Volumes',
//            '</label>'
//          ].join(''));
        
        
        var previousSourceWell = null;
        var SelectedLcpRow = Backgrid.Row.extend({
          render: function() {
            console.log('render');
            
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
                this.$el.closest('tr').css('border-top','2px solid black');
              }
            }
            return SelectedLcpRow.__super__.render.apply(this, arguments);
          }
        });
        
        var view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: schemaResult,
          resource: schemaResult,
          url: url,
          row: SelectedLcpRow,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.listenTo(view, 'afterRender', function(event) {
          view.$el.find('#list-title').show().append(
            '<H4 id="title">Lab Cherry Picks for : ' + self.model.key + '</H4>');
        });
        
        showUnfulfilledWellsControl.click(function(e){
          function processClick(){
            if (e.target.checked) {
              var searchHash = _.clone(view.listModel.get('search'));
              searchHash['show_unfulfilled'] = 'true';
              view.listModel.set('search',searchHash);
            } else {
              var searchHash = _.clone(view.listModel.get('search'));
              delete searchHash['show_unfulfilled'];
              view.listModel.set('search',searchHash);
            }
          };
          if(appModel.isPagePending()){
            appModel.requestPageChange({
              ok: processClick
            });
          }else{
            processClick();
          }
        });
        
        // Display or hide copy wells
        showCopyWellsControl.click(function(e) {
          function processClick(){
            if (e.target.checked) {
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.union(
                ['selected', 'source_copy_well_volume','volume_approved',
                 'source_copy_usage_type','source_plate_status'],includes);
              view.listModel.set({ includes: includes}, {reset: false});
              var searchHash = _.clone(view.listModel.get('search'));
              searchHash['show_copy_wells'] = 'true';

              showAllCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              delete searchHash['show_available_and_retired_copy_wells'];
              view.listModel.set('search',searchHash);
            } else {
              // make sure unset
              var searchHash = _.clone(view.listModel.get('search'));
              if (!_.has(searchHash,'show_copy_wells')) {
                return;
              }
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.without(
                includes,
                'selected','source_copy_well_volume','volume_approved',
                'source_copy_usage_type','source_plate_status');
              view.listModel.set({ includes: includes}, {reset: false});
              if (_.has(searchHash,'show_copy_wells')) {
                delete searchHash['show_copy_wells'];
                view.listModel.set('search',searchHash);
              }
            }
          };
          if(appModel.isPagePending()){
            appModel.requestPageChange({
              ok: processClick
            });
          }else{
            processClick();
          }
        });

        // Display or hide copy wells
        showAllCopyWellsControl.click(function(e) {
          function processClick(){
            if (e.target.checked) {
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.union(
                ['selected', 'source_copy_well_volume','volume_approved',
                 'source_copy_usage_type','source_plate_status'],includes);
              view.listModel.set({ includes: includes}, {reset: false});
              var searchHash = _.clone(view.listModel.get('search'));
              searchHash['show_available_and_retired_copy_wells'] = 'true';
              
              showCopyWellsControl.find('input[type="checkbox"]').prop('checked',false);
              delete searchHash['show_copy_wells'];
              view.listModel.set('search',searchHash);
            } else {
              var searchHash = _.clone(view.listModel.get('search'));
              if (!_.has(searchHash,'show_available_and_retired_copy_wells')) {
                return;
              }
              // make sure unset
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.without(
                includes,
                'selected','source_copy_well_volume','volume_approved',
                'source_copy_usage_type','source_plate_status');
              view.listModel.set({ includes: includes}, {reset: false});
              if (_.has(searchHash,'show_available_and_retired_copy_wells')) {
                delete searchHash['show_available_and_retired_copy_wells'];
                view.listModel.set('search',searchHash);
              }
            }
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

        // Manage selection updates
        
        
        
        
        
        
        
        /////////////////// FIXME REDO like SCP's
        
//        
//        view.collection.on('change', function(model){
//          console.log('collection changed', arguments);
//          selectionUpdateCollection
        
        
        
        
        
        
        
        view.collection.on('change', function(model){
          console.log('collection changed', arguments);
//          var copy_well = Iccbl.formatString('{source_well_id}/{source_copy_name}', model);
          var current_copy_name = model.get('source_copy_name');
          var source_well_id = model.get('source_well_id');
          var source_copy_id = model.get('source_copy_id');
          var selected_copy_name = model.get('selected_copy_name');
          if (model.get('selected')) {
            
            if (current_copy_name === selected_copy_name ){
              delete lcp_selection_updates[source_well_id];
            }else{
              lcp_selection_updates[source_well_id] = source_copy_id;
            }
            
            // unselect others in the group
            _.each(view.collection.where({'source_well_id':source_well_id}),
              function(model){
                if (model.get('source_copy_name')!==current_copy_name){
                  model.set({'selected': false});
                }
            });
          } else {
            // deselection
            // if it was in the updates, then remove it
            if (_.has(lcp_selection_updates, source_well_id)) {
              if (lcp_selection_updates[source_well_id]==source_copy_id){
                delete lcp_selection_updates[source_well_id];
              }
            } else {
              // otherwise, was not already in the updates, so it must be a deselection
              lcp_selection_updates[source_well_id] = null;
            }
          } 
          if (!_.isEmpty(lcp_selection_updates) ) {
            appModel.setPagePending(null, 'Selection updates are unsaved, continue anyway?');
            setSelectedLcpButton.show();
          } else {
            setSelectedLcpButton.hide();
          }
        });
        // Make sure that on reset actions (page changes), that selections are persisted
        view.collection.on('reset', function(){
          console.log('reset event...');
          if (!_.isEmpty(lcp_selection_updates)){
            view.collection.each(function(model){
              var source_well_id = model.get('source_well_id');
              if (_.has(lcp_selection_updates, source_well_id)){
                if (model.get('source_copy_id')!=lcp_selection_updates[source_well_id]){
                  model.set('selected', false);
                } else {
                  model.set('selected', true);
                }
              }
            });
          }
        });
        
        reserveAndMapSelectedButton.click(function(e){
          
          e.preventDefault();
          console.log('submit selections for plating');
          var plate_lab_cherrypicks_url = [
            self.model.resource.apiUri,self.model.key, 
            'reserve_map_lab_cherry_picks'].join('/');
          var headers = {}; // can be used to send a comment

          function processClick(overrideInsufficient){
            
            var data = new FormData();
            // send API_PARAM_VOLUME_OVERRIDE = 'volume_override'
            data.append('volume_override', overrideInsufficient);

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
            
              console.log('submitted', arguments);
              // Show result message
              console.log('done', data);
              // TODO: refactor success message parsing
              if (_.isObject(data) && !_.isString(data)){
                data = _.result(_.result(data,'meta',data),'Result',data);
                var msg_rows = appModel.dict_to_rows(data);
                var bodyMsg = msg_rows;
                if (_.isArray(msg_rows) && msg_rows.length > 1){
                  bodyMsg = _.map(msg_rows, function(msg_row){
                    return msg_row.join(': ');
                  }).join('<br>');
                }
                var title = 'Lab Cherry Pick Plating result';
                appModel.showModalMessage({
                  body: bodyMsg,
                  title: title  
                });
              }else{
                console.log('data should have been parsed as json', data);
                appModel.showModalMessage({
                  title: 'Lab Cherry Pick Plating result',
                  okText: 'ok',
                  body: 'success'
                });
              }
              console.log('reload...');

              // Show the cherry pick plate tab
              self.tabbed_resources['cherrypickplates'] = 
                  self.cherry_pick_tabbed_resources['cherrypickplates'];
              self.uriStack = ['labcherrypicks'];
              // remove the child view before calling render, to prevent
              // it from being rendered twice, and calling afterRender twice
              self.removeView('#tab_container');
              self.render();
              
            }).fail(function(jqXHR, textStatus, errorThrown){
              console.log('errors', arguments);
              
              var overridden = false;
              var jsonError = _.result(jqXHR, 'responseJSON');
              if (!_.isUndefined(jsonError)){
                var error = _.result(jsonError, 'errors');
                error = _.result(error,'transfer_volume_per_well_approved');
                if (!_.isUndefined(error)){
                  // Check for "volume_override" required message
                  if (error.length==1 && error[0].includes('volume_override')){
                    overridden = true;
                    appModel.showModal({
                      title: 'Insufficient volume for copy wells in lab cherry picks',
                      body: error,
                      okText: 'Override',
                      ok: function() {
                        var overrideInsufficient = true;
                        processClick(overrideInsufficient);
                      },
                      cancel: function() {
                        appModel.jqXHRfail.apply(jqXHR,textStatus,errorThrown); 
                      }
                    });
                  }
                }
              }
              if (overridden == false){
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
        
        setSelectedLcpButton.click(function(e){
          e.preventDefault();
          console.log('selected', lcp_selection_updates);
          if (_.isEmpty(lcp_selection_updates)) {
            appModel.error('No changes to save');
            return;
          }
          self.saveCopyWellsSelected(lcp_selection_updates)
            .done(function(data, textStatus, jqXHR){
              console.log('success', data);
              // TODO: refactor success message parsing
              if (_.isObject(data) && !_.isString(data)){
                data = _.result(_.result(data,'meta',data),'Result',data);
                var msg_rows = appModel.dict_to_rows(data);
                var bodyMsg = msg_rows;
                if (_.isArray(msg_rows) && msg_rows.length > 1){
                  bodyMsg = _.map(msg_rows, function(msg_row){
                    return msg_row.join(': ');
                  }).join('<br>');
                }
                var title = 'Lab Cherry Pick Copy Selection updates';
                appModel.showModalMessage({
                  body: bodyMsg,
                  title: title  
                });
              }else{
                console.log('data should have been parsed as json', data);
                appModel.showModalMessage({
                  title: 'Lab Cherry Pick Copy Selection updates',
                  okText: 'ok',
                  body: 'success'
                });
              }
              lcp_selection_updates = {};
              setSelectedLcpButton.hide();
              view.collection.fetch({ reset: true });
            })
            .fail(function(jqXHR, textStatus, errorThrown){
              appModel.jqXHRfail.apply(this,arguments); 
            });
        });
      };
      
      var schemaUrl = [
        self.model.resource.apiUri,self.model.key,'lab_cherry_pick',
        'schema'].join('/');
      appModel.getResourceFromUrl(schemaUrl, createLcpView);
    },
    
    saveCopyWellsSelected: function(lcp_selection_updates) {
      var self = this;
      var promise = $.Deferred();
      
      appModel.showSaveWithComments(function(formValues){
        console.log('form values', formValues);
        var comments = formValues['comments'];
        var headers = {};
        headers[appModel.HEADER_APILOG_COMMENT] = comments;
        
        var cp_patch_url = [
          self.model.resource.apiUri,self.model.key].join('/');
        $.ajax({
          url: cp_patch_url,    
          data: JSON.stringify({
            lcp_selection_updates: lcp_selection_updates
          }),
          cache: false,
          contentType: 'application/json', 
          dataType: 'json', // what is expected back from the server
          type: 'PATCH',
          headers: headers
        }).done(function(data, textStatus, jqXHR){ 
          promise.resolveWith(this,arguments); 
        }).fail(function(jqXHR, textStatus, errorThrown){
          promise.rejectWith(this,arguments);
        })
      });
      // return a restricted Promise
      return promise.promise();
    },
    
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
      console.log('no screener cherry picks entered...');
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
          form.$el.find('[ type="submit" ]').click(function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }); 
            if(!_.isEmpty(errors)){
              form.$el.find('#well_search').addClass(self.errorClass);
              return;
            }else{
              form.$el.find('#well_search').removeClass(self.errorClass);
            }
            
            var modelForSubmit = new Backbone.Model(
              self.model.pick(['cherry_pick_request_id'])
            );
            modelForSubmit.set('id', self.model.get('id'));
            modelForSubmit.resource = self.model.resource;
            modelForSubmit.urlRoot = self.model.resource.apiUri;
            modelForSubmit.set('screener_cherry_picks', form.getValue('well_search'));
            
            modelForSubmit.save(null, {
              patch: true
            }).done(function(data, textStatus, jqXHR){ 
              console.log('success', data);
              if (_.isObject(data) && !_.isString(data)){
                data = _.result(_.result(data,'meta',data),'Result',data);
                var msg_rows = appModel.dict_to_rows(data);
                var bodyMsg = msg_rows;
                if (_.isArray(msg_rows) && msg_rows.length > 1){
                  bodyMsg = _.map(msg_rows, function(msg_row){
                    return msg_row.join(': ');
                  }).join('<br>');
                }
                var title = 'Screener Cherry Picks';
                appModel.showModalMessage({
                  body: bodyMsg,
                  title: title  
                });
              }else{
                console.log('data should have been parsed as json', data);
                appModel.showModalMessage({
                  title: 'Screener Cherry Picks',
                  okText: 'ok',
                  body: 'uploaded...'
                });
              }
              self.model.fetch({ reset: true }).done(function(){
                self.removeView('#tab_container');
                self.change_to_tab('screenercherrypicks');
              });
            }).fail(function(jqXHR, textStatus, errorThrown) { 
              appModel.jqXHRfail.apply(this,arguments); 
            });
            
          });
        }
      });
      var view = new View();
      self.setView('#tab_container', view).render();
      self.reportUriStack([]);
    },
    
    // Construct the ScreenerCherryPick tab/grid
    createScpView: function(schemaResult, delegateStack){
      var self = this;
      
      var url = [self.model.resource.apiUri,self.model.key,
                 'screener_cherry_pick'].join('/');
      
      var selection_updates = [];
      var deselect_updates = [];
      var extraControls = [];
      
      var Collection = Backbone.Collection.extend({
          modelId: function(attrs) {
            return Iccbl.getIdFromIdAttribute( attrs, schemaResult);
          }
      })
      var selectionUpdateCollection = new Collection();
      
      // Show the delete button
      // TODO: When should delete be disabled?
      if (self.model.get('number_plates') == 0){
        console.log('add delete button...');
        var deleteScpsButton = $([
            '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="deleteScpsButton" href="#">',
            'Delete Screener Cherry Picks</a>'
          ].join(''));
        extraControls.push(deleteScpsButton);
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
          }).done(function(model, resp) {
            // TODO: done replaces success as of jq 1.8
            console.log('delete SCPs completed');
            self.model.fetch({ reset: true }).done(function(){
              self.removeView('#tab_container');
              self.change_to_tab('screenercherrypicks');
            });
          }).fail(function(jqXHR, textStatus, errorThrown) { 
            appModel.jqXHRfail.apply(this,arguments); 
          });
        });
      }
      var showOtherReagentsControl = $([
          '<label class="checkbox-inline">',
          '  <input type="checkbox">showOtherReagents',
          '</label>'
        ].join(''));
      extraControls.push(showOtherReagentsControl);

      var setLcpsButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="setLcpsButton" href="#">',
          'Set Lab Cherry Picks</a>'
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
      
      if(self.model.get('has_pool_screener_cherry_picks') === true){
        //schemaResult.fields['pool_reagent_vendor_id']['visibility'] = ['l','d'];
        setDuplexLcpsButton.show();
      }
      if(self.model.get('has_alternate_screener_cherry_pick_selections') === true){
        schemaResult['fields']['searched_well_id']['visibility'] = ['l','d'];
        schemaResult['fields']['selected']['visibility'] = ['l','d'];
      }      

      var SelectedScpRow = Backgrid.Row.extend({
        initialize: function () {
          var self = this;
          SelectedScpRow.__super__.initialize.apply(this, arguments);
          
          this.listenTo(this.model.collection, 'show_other_reagents', function(){
            console.log('event:', arguments);
            self._setStyle();
          });
        },
        
        _setStyle: function(){
          var searchHash = this.model.collection.listModel.get('search');
          if(self.model.get('has_alternate_screener_cherry_pick_selections') === true
              || _.result(searchHash,'show_other_reagents')=='true')
          {
            if (this.model.get('searched_well_id')
                ==this.model.get('screened_well_id')) {
              this.$el.closest('tr').addClass('selector_row_group');
            }
          } else {
            this.$el.closest('tr').removeClass('selector_row_group');
          }
          if (_.result(searchHash,'show_other_reagents')=='true'){
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

      var ListViewSelect = ListView.extend({
        afterRender: function(){
          ListViewSelect.__super__.afterRender.apply(this, arguments);
          // Backgrid specific: all column headers are given a class for their
          // column name: in this case the "selected" class has other meaning, 
          // so remove it
          this.$('th').removeClass('selected');
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
      
      // submit as Lab Cherry Picks
      setLcpsButton.click(function(e){
        e.preventDefault();

        function processClick(){
          
          var set_lab_cherry_picks_url = [
            self.model.resource.apiUri,self.model.key, 'set_lab_cherry_picks'].join('/');
          var headers = {}; // can be used to send a comment
          $.ajax({
            url: set_lab_cherry_picks_url,     
  //          data: JSON.stringify({
  //            selection_updates: selection_updates,
  //            deselect_updates: deselect_updates
  //          }),
            cache: false,
            contentType: 'application/json', 
            dataType: 'json', // what is expected back from the server
            type: 'POST',
            headers: headers
          }).done(function(data, textStatus, jqXHR){
            console.log('submitted', arguments);
            self.change_to_tab('labcherrypicks');
          }).fail(function(jqXHR, textStatus, errorThrown){
            appModel.jqXHRfail.apply(this,arguments); 
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

      // submit as Lab Cherry Picks
      setDuplexLcpsButton.click(function(e){
        e.preventDefault();
        var set_lab_cherry_picks_url = [
          self.model.resource.apiUri,self.model.key, 
          'set_duplex_lab_cherrypicks'].join('/');
        var headers = {}; // can be used to send a comment
        $.ajax({
          url: set_lab_cherry_picks_url,     
//          data: JSON.stringify({
//            selection_updates: selection_updates,
//            deselect_updates: deselect_updates
//          }),
          cache: false,
          contentType: 'application/json', 
          dataType: 'json', // what is expected back from the server
          type: 'POST',
          headers: headers
        }).done(function(data, textStatus, jqXHR){
          console.log('submitted', arguments);
          self.change_to_tab('labcherrypicks');
        }).fail(function(jqXHR, textStatus, errorThrown){
          appModel.jqXHRfail.apply(this,arguments); 
        });
      });

      // Manage selection updates
      view.collection.on('add', function(model){
        // cache the 'selected' property for update management
        model.set(
          { selected_on_server: model.get('selected') },
          { silent: true }
        );
      });
      view.collection.on('change', function(model){
        console.log('collection changed', arguments);
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
        console.log('reset event...');
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
      setSelectedButton.click(function(e){
        e.preventDefault();
        console.log('selection updates', selectionUpdateCollection);
        if (selectionUpdateCollection.isEmpty()) {
          appModel.error('No changes to save');
          return;
        }
        selectionUpdateCollection.url = view.collection.url;
        selectionUpdateCollection.sync('patch', selectionUpdateCollection)
          .done(function(data, textStatus, jqXHR){
            console.log('success', data);
            if (_.isObject(data) && !_.isString(data)){
              data = _.result(_.result(data,'meta',data),'Result',data);
              var msg_rows = appModel.dict_to_rows(data);
              var bodyMsg = msg_rows;
              if (_.isArray(msg_rows) && msg_rows.length > 1){
                bodyMsg = _.map(msg_rows, function(msg_row){
                  return msg_row.join(': ');
                }).join('<br>');
              }
              var title = 'Screener Cherry Picks';
              appModel.showModalMessage({
                body: bodyMsg,
                title: title  
              });
            }else{
              console.log('data should have been parsed as json', data);
              appModel.showModalMessage({
                title: 'Screener Cherry Picks',
                okText: 'ok',
                body: 'uploaded...'
              });
            }
            
            setSelectedButton.hide();
            selectionUpdateCollection.reset(null); // clear
            view.collection.fetch({ reset: true });

          }).fail(function(jqXHR, textStatus, errorThrown){
            console.log('fail', arguments);
            appModel.jqXHRfail.apply(this,arguments); 
          });
      });
      
      // Display or hide other reagents
      showOtherReagentsControl.click(function(e) {
        function processClick(){
          if (e.target.checked) {
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.union(
                ['selected', 'searched_well_id','mg_ml_concentration',
                 'molar_concentration'],includes);
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['show_other_reagents'] = 'true';
            view.listModel.set('search',searchHash);
          } else {
            // make sure unset
            var includes = _.clone(view.listModel.get('includes'));
            includes = _.without(
                includes,
                'selected', 'searched_well_id','mg_ml_concentration',
                'molar_concentration');
            view.listModel.set({ includes: includes}, {reset: false});
            var searchHash = _.clone(view.listModel.get('search'));
            if (_.has(searchHash,'show_other_reagents')) {
              delete searchHash['show_other_reagents'];
              view.listModel.set('search',searchHash);
            }
          }
          // see note above about removing the 'selected' backgrid th class
          view.$('th').removeClass('selected');
          view.collection.trigger('show_other_reagents');
          
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
      if (_.has(initialSearchHash, 'show_other_reagents')
          && initialSearchHash.show_other_reagents.toLowerCase()=='true') {
        showOtherReagentsControl.find('input[type="checkbox"]').prop('checked',true);
      }
    },
  
  });
  
  return CherryPickView;
});