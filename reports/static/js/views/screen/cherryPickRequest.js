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
      this._classname = 'CherryPickRequest';
      TabbedController.prototype.initialize.apply(this,arguments);
      _.bindAll(this, 'createScpView','showScpSearchForm','showWellsToLeaveEmptyDialog');
      
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

          var wellsToLeaveEmptyButton = $([
            '<a class="btn btn-default btn-sm" ',
              'role="button" id="wellsToLeaveEmpty_button" href="#">',
              'Well picker</a>'
            ].join(''));
          wellsToLeaveEmptyButton.click(function(event) {
            event.preventDefault();
            self.showWellsToLeaveEmptyDialog(editForm);
          });
          var fieldKey = 'wells_to_leave_empty';
          this.$el.find('div[data-fields="' + fieldKey + '"]')
            .find('div[data-editor]').append(wellsToLeaveEmptyButton);
          
          EditView.prototype.afterRender.apply(this,arguments);
        }
      });
      var saveSuccessCallBack = function(model){
        // TODO: modify generic_edit to handle custom-nested URLs for edit success URL
        var meta = _.result(model, 'meta', null);
        if (meta) {
          appModel.showJsonMessages(meta);
        }
        model = new Backbone.Model(model);
        var key = Iccbl.getIdFromIdAttribute( model,self.model.resource );
        model.key = self.model.resource.key + '/' + key;
        appModel.router.navigate([
          self.screen.resource.key,self.screen.key,'cherrypickrequest',key].join('/'), 
          {trigger:true});
      };
      
      view = new DetailLayout(_.extend(self.args, { 
        model: this.model,
        uriStack: delegateStack, 
        buttons: buttons,
        saveSuccessCallBack: saveSuccessCallBack,
        EditView: editView
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

    showWellsToLeaveEmptyDialog: function(editForm){
      var self = this;
      var rowCollection = new Backbone.Collection();
      
      var wellSelector = new WellSelector({
        wellSelections: self.model.get('wells_to_leave_empty')
      });
      var el = wellSelector.render().el;
      appModel.showModal({
        width: 900,
        okText: 'ok',
        view: el,
        title: 'Select wells to leave empty',
        ok: function(e) {
          e.preventDefault();
          console.log('selected wells: ', wellSelector.getSelectedWells());
          editForm.setValue('wells_to_leave_empty', 
            wellSelector.getSelectedWells().join(', '));
          self.model.set(
            'wells_to_leave_empty', 
            wellSelector.getSelectedWells().join(','));
        }              
      });
    },
    
    setLabCherryPicks: function(delegateStack) {
      var self = this;

      function createLcpView(schemaResult){
        var url = [self.model.resource.apiUri,self.model.key,
          'lab_cherry_pick'].join('/');

        if(self.model.get('has_pool_screener_cherry_picks') === true){
          schemaResult.fields['pool_reagent_vendor_id']['visibility'] = ['l','d'];
        }
        
        var extraControls = [];
        var setSelectedButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'style="display: none; " ',
            'role="button" id="save_button_selected" href="#">',
            'Save Selections</a>'
          ].join(''));
        extraControls.push(setSelectedButton);
        var showCopyWellsControl = $([
            '<label class="checkbox-inline">',
            '  <input type="checkbox">showCopyWells',
            '</label>'
          ].join(''));
        extraControls.push(showCopyWellsControl);

        var SelectedLcpRow = Backgrid.Row.extend({
          render: function() {
            console.log('render');
            if (this.model.has('selected_copy_name')) {
              if (this.model.get('selected_copy_name')==this.model.get('source_copy_name')) {
                this.$el.addClass('selected');
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
              view.listModel.set('search',searchHash);
            } else {
              // make sure unset
              var includes = _.clone(view.listModel.get('includes'));
              includes = _.without(
                includes,
                'selected','source_copy_well_volume','volume_approved',
                'source_copy_usage_type','source_plate_status');
              view.listModel.set({ includes: includes}, {reset: false});
              var searchHash = _.clone(view.listModel.get('search'));
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
        var initialSearchHash = view.listModel.get('search');
        if (_.has(initialSearchHash, 'show_copy_wells')
            && initialSearchHash.show_copy_wells.toLowerCase()=='true') {
          showCopyWellsControl.find('input[type="checkbox"]').prop('checked',true);
        }

        // Manage selection updates
        var selection_updates = {};
        view.collection.on('change', function(model){
          console.log('collection changed', arguments);
          var copy_well = Iccbl.formatString('{source_well_id}/{source_copy_name}', model);
          if (model.get('selected')) {
            var current_copy_name = model.get('source_copy_name');
            var source_well_id = model.get('source_well_id');
            var selected_copy_name = model.get('selected_copy_name');
            
            if (current_copy_name === selected_copy_name ){
              delete selection_updates[source_well_id];
            }else{
              selection_updates[source_well_id] = current_copy_name;
            }
            
            // unselect others in the group
            _.each(view.collection.where({'source_well_id':source_well_id}),
              function(model){
                if (model.get('source_copy_name')!==current_copy_name){
                  model.set('selected', false);
                }
            });
          } 
          if (!_.isEmpty(selection_updates) ) {
            appModel.setPagePending();
            setSelectedButton.show();
          } else {
            setSelectedButton.hide();
          }
        });
        // Make sure that on reset actions (page changes), that selections are persisted
        view.collection.on('reset', function(){
          console.log('reset event...');
          if (!_.isEmpty(selection_updates)){
            view.collection.each(function(model){
              var source_well_id = model.get('source_well_id');
              if (_.has(selection_updates, source_well_id)){
                if (model.get('source_copy_name')!=selection_updates[source_well_id]){
                  model.set('selected', false);
                } else {
                  model.set('selected', true);
                }
              }
            });
          }
        });
        setSelectedButton.click(function(e){
          e.preventDefault();
          console.log('selected', selection_updates,deselect_updates);
          if (_.isEmpty(selection_updates) && _.isEmpty(deselect_updates)) {
            appModel.error('No changes to save');
            return;
          }
          self.saveCopyWellsSelected(selection_updates,deselect_updates)
            .done(function(){
              selection_updates = [];
              deselect_updates = [];
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
      
      var extraControls = [];
      
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
        schemaResult.fields['pool_reagent_vendor_id']['visibility'] = ['l','d'];
        setDuplexLcpsButton.show();
      }
      
      var SelectedScpRow = Backgrid.Row.extend({
        render: function() {
          console.log('render');
          if (this.model.has('searched_well_id')) {
            if (this.model.get('searched_well_id')==this.model.get('screened_well_id')) {
              this.$el.addClass('selected');
            }
          }
          return SelectedScpRow.__super__.render.apply(this, arguments);
        }
      });
      
      var view = new ListView({ 
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
        var set_lab_cherrypicks_url = [
          self.model.resource.apiUri,self.model.key, 'set_lab_cherrypicks'].join('/');
        var headers = {}; // can be used to send a comment
        $.ajax({
          url: set_lab_cherrypicks_url,     
          data: JSON.stringify({
            selection_updates: selection_updates,
            deselect_updates: deselect_updates
          }),
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
        })
      });

      // Manage selection updates
      var selection_updates = {}; // TODO
      var selection_updates = [];
      var deselect_updates = [];
      view.collection.on('change', function(model){
        console.log('collection changed', arguments);
        if (model.get('selected')) {
          var current_well_search = model.get('searched_well_id');
          var current_well = model.get('screened_well_id');
          selection_updates.push(model.get('screened_well_id'));
          deselect_updates = _.without(
            deselect_updates, 
            model.get('screened_well_id'));

          // unselect others in the group
          _.each(view.collection.where({'searched_well_id':current_well_search}),
            function(model){
              if (model.get('screened_well_id')!==current_well){
                model.set('selected', false);
              }
          });
          
        } else {
          selection_updates = _.without(
            selection_updates, 
            model.get('screened_well_id'));
          deselect_updates.push(model.get('screened_well_id'));
        }
        if (!_.isEmpty(selection_updates) || !_.isEmpty(deselect_updates)) {
          appModel.setPagePending();
          setSelectedButton.show();
        } else {
          setSelectedButton.hide();
        }
      });
      
      // Make sure that on reset actions (page changes), that selections are persisted
      view.collection.on('reset', function(){
        console.log('reset event...');
        if (!_.isEmpty(selection_updates)){
          view.collection.each(function(model){
            if (_.contains(selection_updates, model.get('screened_well_id'))){
              model.set('selected', true);
            }
          });
        }
        if (!_.isEmpty(deselect_updates)) {
          view.collection.each(function(model){
            if (_.contains(deselect_updates, model.get('screened_well_id'))){
              model.set('selected', false);
            }
          });
        }
      });
      setSelectedButton.click(function(e){
        e.preventDefault();
        console.log('selected', selection_updates,deselect_updates);
        if (_.isEmpty(selection_updates) && _.isEmpty(deselect_updates)) {
          appModel.error('No changes to save');
          return;
        }
        
        self.saveScreenerCherryPicksSelected(selection_updates,deselect_updates)
          .done(function(){
            selection_updates = [];
            deselect_updates = [];
            view.collection.fetch({ reset: true });
          })
          .fail(function(jqXHR, textStatus, errorThrown){
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
    
    saveScreenerCherryPicksSelected: function(
        selection_updates, deselect_updates) {
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
            selection_updates: selection_updates,
            deselect_updates: deselect_updates
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
    }
    
  });
  
  return CherryPickView;
});