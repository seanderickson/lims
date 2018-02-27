define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/list2',
  'utils/treeSelector',
  'utils/tabbedController'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, 
            ListView, TreeSelector, TabbedController) {
  
  var ScreenDataView = TabbedController.extend({
    
    initialize: function(args) {
      var self = this;
      this._classname = 'ScreenDataView';
      
      this.tabbed_resources = this.screen_tabbed_resources;
      if (! _.isEmpty(this.model.get('study_type'))) {
        this.tabbed_resources = this.study_tabbed_resources;
      } 
      
      TabbedController.prototype.initialize.apply(this,arguments);
      
      if (this.model.get('user_access_level_granted') > 1){
        if (! _.isEmpty(this.model.get('study_type'))) {
          this.tabbed_resources['datacolumns'] = this.study_tabbed_resources['datacolumns']
        } else {
          this.tabbed_resources['datacolumns'] = this.screen_tabbed_resources['datacolumns']
        }
      }
    },

    screen_tabbed_resources: {
      detail: { 
        description: 'Result Data', 
        title: 'Screen Results', 
        invoke: 'setDetail'
      },
      datacolumns: {
        description : 'Data Columns',
        title : 'Data Columns',
        invoke : 'setDatacolumns',
        permission: 'screenresult'
      }
    },
    study_tabbed_resources: {
      detail: { 
        description: 'Reagent Annotations', 
        title: 'Reagent Annotations', 
        invoke: 'setDetail'
      },
      datacolumns: {
        description : 'Data Columns',
        title : 'Data Columns',
        invoke : 'setDatacolumns',
        permission: 'screenresult'
      }
    },
    
    events: {
      'click ul.nav-tabs >li': 'click_tab',
      'click button#loadScreenResults': 'loadScreenResults',
      'click button#deleteScreenResults': 'deleteScreenResults',
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      return {
        'base_url': self.model.resource.key + '/' + self.model.key + '/data',
        'tab_resources': self.tabbed_resources
      }      
    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function() {
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)) {
        viewId = this.uriStack.shift();

        if (_.contains(appModel.LIST_ARGS, viewId) ) {
          this.uriStack.unshift(viewId);
          viewId = 'detail'
        }

        if (!_.has(this.tabbed_resources, viewId)) {
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },    
    
    setDetail: function(delegateStack) {

      console.log('set screen result detail...', delegateStack);
      
      var self = this;
      var key = 'detail';

      var $loadScreenResultsButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
        id="loadScreenResults" >Load Data</button>');
      var $deleteScreenResultsButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
        id="deleteScreenResults" >Delete Data</button>');
      
      
      var screenResultResource = appModel.getResource('screenresult');
//      if (! _.isEmpty(this.model.get('study_type'))) {
//        screenResultResource['fields']['assay_well_control_type']['visibility'] = [];
//      }

      var schemaUrl = [appModel.dbApiUri,
                       'screenresult',
                       self.model.key,
                       'schema'].join('/');
      var url = [appModel.dbApiUri,
                 'screenresult',
                 self.model.key].join('/');

      var _id = self.model.key;
      var screen_facility_id = self.model.get('facility_id');
      
      var createResults = function(schemaResult) {
        var extraControls = [];
        var show_positives_control = $([
          '<label class="checkbox-inline">',
          '  <input type="checkbox">show positive rows only',
          '</label>'
          ].join(''));
        var show_mutual_positives_control = $([
          '<label class="checkbox-inline">',
           '  <input type="checkbox">mutual positives',
           '</label>'
           ].join(''));
//        var show_mutual_positives_control = $([
//          '<button class="btn btn-default btn-sm" role="button" ',
//          'id="showMutualPositives" title="Show mutual positive columns" >',
//          'Show mutual positive columns',
//          '</button>'
//           ].join(''));
        var show_other_screen_columns_button = $([
          '<button class="btn btn-default btn-sm pull-right" role="button" ',
          'id="showOtherScreenColumns" title="Show other screen columns" >',
          'Other screen columns',
          '</button>'
          ].join(''))
        if (_.isEmpty(self.model.get('study_type'))) {
          extraControls = extraControls.concat(
            show_positives_control, show_mutual_positives_control);
          extraControls.push($('<span>&nbsp;</span>'));
        }
        if (appModel.hasPermission('screenresult','write')) {
          extraControls = extraControls.concat(
            $deleteScreenResultsButton,$loadScreenResultsButton);
        }
        
        // create an option vocab for the excluded col, if needed
        if (_.has(schemaResult['fields'], 'excluded')) {
          var options = [];
          _.each(schemaResult['fields'],function(field) {
            var dc_col_key = 'dc_' + self.model.key + '_';
            if (field['key'].indexOf(dc_col_key)>-1) {
              var option_name = field['key'].substring(dc_col_key.length);
              options.unshift([option_name, option_name]);
            }
          });
          schemaResult['fields']['excluded']['vocabulary'] = options;
          console.log('custom exclude options', options);
        }
        var initialSearchHash;
        var SRListView = ListView.extend({
          afterRender: function(){
            ListView.prototype.afterRender.apply(this,arguments);
            this.$('#list_controls').append(show_other_screen_columns_button);
            return this;
          }
        });
        
        view = new SRListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: schemaResult,
          resource: screenResultResource,
          url: url,
          extraControls: extraControls,
          screen: self.model
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        
        initialSearchHash = view.listModel.get('search');
        
        if (_.has(initialSearchHash, 'is_positive__eq')
            && initialSearchHash.is_positive__eq.toLowerCase()=='true') {
          show_positives_control.find('input[type="checkbox"]').prop('checked',true);
        }
        // 20170905 - removed show_mutual_positives - other screen columns will
        // be added explicitly using treeSelector
        // 20180220 - reinstated, determine state by included cols
        if (_.has(initialSearchHash, 'show_mutual_positives')
            && initialSearchHash.show_mutual_positives.toLowerCase()=='true') {
          show_mutual_positives_control.find('input[type="checkbox"]').prop('checked',true);
          self.showMutualPositiveColumns(view, true);
        }

        show_other_screen_columns_button.click(function(e){
          self.showOtherScreenColumnsDialog(view);
        });
        show_positives_control.find('input[type="checkbox"]').change(function(e) {
          if (this.checked) {
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['is_positive__eq'] = 'true';
            view.listModel.set('search',searchHash);
          } else {
            // make sure unset
            var searchHash = _.clone(view.listModel.get('search'));
            if (_.has(searchHash,'is_positive__eq')) {
              delete searchHash['is_positive__eq'];
              view.listModel.set('search',searchHash);
            }
          }
        });
        // 20170905 - removed: mutual pos columns will be selected explicitly
        // 20180220 - reinstated, determine state by included cols
        show_mutual_positives_control.find('input[type="checkbox"]').change(function(e) {
          if (this.checked) {
            window.setTimeout(function() {
              self.showMutualPositiveColumns(view, true);
            });
          } else {
            window.setTimeout(function() {
              self.showMutualPositiveColumns(view, false);
            });
          }
        });
//        show_mutual_positives_control.click(function(e) {
//          var searchHash = _.clone(view.listModel.get('search'));
//          searchHash['is_positive__eq'] = 'true';
//          view.listModel.set('search',searchHash);
//          show_positives_control.find('input[type="checkbox"]').prop('checked',true);
//          self.showMutualPositiveColumns(view);
//        });
        
      }; // createResults
      
      if (self.model.get('has_screen_result')) {
        var listOptions = ListView.prototype.parseUriStack(delegateStack);
        console.log('parsed listOptions', listOptions);
        var dc_ids = [];
        if (_.has(listOptions, 'search')){
          if (_.has(listOptions.search,'dc_ids')){
            dc_ids = listOptions.search.dc_ids;
            if (dc_ids.length > 0 && _.isString(dc_ids)){
              dc_ids = dc_ids.split(',');
            }
          }
        }
        dc_ids = _.map(dc_ids, function(dc_id){ return parseInt(dc_id); });
        console.log('dc_ids', dc_ids);
        var options = {};
        if (!_.isEmpty(dc_ids)){
          options['dc_ids'] = dc_ids;
        }
        appModel.getResourceFromUrl(schemaUrl, createResults, options);
      } else {
        if (appModel.hasPermission('screenresult','write')) {
          this.$('#tab_container').append("No data loaded ");
          this.$('#tab_container').append($loadScreenResultsButton);
          self.reportUriStack([]);
        }
      }
    },
    
    showMutualPositiveColumns: function(listView, shown){
      var self = this;

//      var searchHash = listView.collection.listModel.get('search');
//      var dc_ids = _.result(searchHash,'dc_ids', '');
//      if (!_.isArray(dc_ids)){
//        dc_ids = dc_ids.split(',')
//      }
//      dc_ids = _.map(dc_ids, function(dc_id){ return parseInt(dc_id); });
//      console.log('showMutualPositiveColumns', dc_ids);
      
      // Retrieve all columns visible for this screen
      // - check the positive overlapping columns
      var srResource = self.model.resource;
      var resource = appModel.getResource('datacolumn');
      var url = [resource.apiUri, 
                 'for_screen', self.model.get('facility_id')].join('/');
      var data_for_get = { 
        limit: 0,
        includes: [
          'screen_facility_id','screen_title','name','description',
          'assay_data_type','ordinal','vocabulary_scope_ref','edit_type'],
        order_by: ['-study_type','screen_facility_id', 'ordinal'],
        use_vocabularies: false
      };
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: url,
        modelId: function(attrs) {
          return Iccbl.getIdFromIdAttribute( attrs, resource);
        }
      });
      var collection = new CollectionClass();
      collection.fetch({
        data: data_for_get,
        success: function(collection, response) {
          if (shown==true){
            collection.each(function(model){
  //            if(_.contains(dc_ids, model.get('data_column_id'))){
  //              model.set('checked', true);
  //            }
              if (_.contains(self.model.get('overlapping_positive_screens'),
                    model.get('screen_facility_id'))){
                if (model.get('positives_count')>0){
                  model.set('checked',true);
                }
              }
              if (self.model.get('screen_type') == 'rnai'){
                if (model.get('screen_facility_id') == '200001'){
                  // Include the study: Reagent counts for Small Molecule screens
                  model.set('checked', true);
                }
                else if (model.get('screen_facility_id') == '200002'){
                  // Include the study: Reagent counts for RNAi screens
                  model.set('checked', true);
                }
              }
            });
            var dcs_selected = collection.where({checked: true});
            _.each(dcs_selected, function(dc){
              var key = dc.get('key');
              var name = dc.get('title');
              var dc_id = dc.get('data_column_id');
              var field = dc.attributes;
              field.description = Iccbl.formatString(
                '{screen_facility_id}: {screen_title} ({name} - {description})',
                field);
              var grid = listView.grid;
              
              if (dc.get('user_access_level_granted') > 1){
                field['filtering'] = true;
                field['ordering'] = true;
              }
              if (!_.has(srResource.fields, key)){
                srResource.fields[key] = field;
              }
    
              var column = grid.columns.findWhere({ name: key });
              if (!column){
                var index = grid.columns.size();
                grid.insertColumn(
                  Iccbl.createBackgridColumn(
                      key,field,
                      listView.listModel.get('order')),
                      { at: index});
              } else {
                console.log('column already included', key)
              }
            });
//          var dc_ids_selected = _.map(dcs_selected, function(dcmodel){
//            return dcmodel.get('data_column_id');
//          });
          
            var searchHash = _.clone(
              listView.collection.listModel.get('search'));
  //          searchHash['dc_ids'] = dc_ids_selected;
            searchHash['show_mutual_positives'] = true;
            listView.collection.listModel.set('search', searchHash);
          } else { // not shown
            collection.each(function(model){
              if (_.contains(self.model.get('overlapping_positive_screens'),
                    model.get('screen_facility_id'))){
                if (model.get('positives_count')>0){
                  var key = model.get('key');
                  if (_.has(srResource.fields, key)){
                    delete srResource.fields[key];
                    var column = listView.grid.columns.findWhere({ name: key });
                    if (!column){
                      console.log('column already not present', key)
                    } else {
                      listView.grid.removeColumn(column);
                    }
                  }
                }
              }
              
            });
            var searchHash = _.clone(
              listView.collection.listModel.get('search'));
            searchHash['show_mutual_positives'] = false;
            listView.collection.listModel.set('search', searchHash);
          }          
        },
        always: function(){
          console.log('done: ');
        }
      }).fail(function(){ 
        Iccbl.appModel.jqXHRfail.apply(this,arguments); 
      });        
    },
    
//    showMutualPositiveColumns_bak_by_dc_ids: function(listView){
//      var self = this;
//
//      var searchHash = listView.collection.listModel.get('search');
//      var dc_ids = _.result(searchHash,'dc_ids', '');
//      if (!_.isArray(dc_ids)){
//        dc_ids = dc_ids.split(',')
//      }
//      dc_ids = _.map(dc_ids, function(dc_id){ return parseInt(dc_id); });
//      console.log('showMutualPositiveColumns', dc_ids);
//      var srResource = self.model.resource;
//      var resource = appModel.getResource('datacolumn');
//      var url = [resource.apiUri, 
//                 'for_screen', self.model.get('facility_id')].join('/');
//      var data_for_get = { 
//        limit: 0,
//        includes: [
//          'screen_facility_id','screen_title','name','description',
//          'assay_data_type','ordinal'],
//        order_by: ['screen_facility_id', 'ordinal'],
//        use_vocabularies: true
//      };
//      var CollectionClass = Iccbl.CollectionOnClient.extend({
//        url: url,
//        modelId: function(attrs) {
//          return Iccbl.getIdFromIdAttribute( attrs, resource);
//        }
//      });
//      var collection = new CollectionClass();
//      collection.fetch({
//        data: data_for_get,
//        success: function(collection, response) {
//          collection.each(function(model){
//            if(_.contains(dc_ids, model.get('data_column_id'))){
//              model.set('checked', true);
//            }
//            if (_.contains(self.model.get('overlapping_positive_screens'),
//                  model.get('screen_facility_id'))){
//              if (model.get('positives_count')>0){
//                model.set('checked',true);
//              }
//            }
//          });
//          var dcs_selected = collection.where({checked: true});
//          _.each(dcs_selected, function(dc){
//            var key = dc.get('key');
//            var name = dc.get('title');
//            var dc_id = dc.get('data_column_id');
//            var field = dc.attributes;
//            var grid = listView.grid;
//            
//            if (dc.get('user_access_level_granted') > 1){
//              field['filtering'] = true;
//              field['ordering'] = true;
//            }
//            if (!_.has(srResource.fields, key)){
//              srResource.fields[key] = field;
//            }
//  
//            var column = grid.columns.findWhere({ name: key });
//            if (!column){
//              var index = grid.columns.size();
//              grid.insertColumn(
//                Iccbl.createBackgridColumn(
//                    key,field,
//                    listView.listModel.get('order')),
//                    { at: index});
//            } else {
//              console.log('column already included', key)
//            }
//          
//          });
//          var dc_ids_selected = _.map(dcs_selected, function(dcmodel){
//            return dcmodel.get('data_column_id');
//          });
//          
//          var searchHash = _.clone(
//            listView.collection.listModel.get('search'));
//          searchHash['dc_ids'] = dc_ids_selected;
//          listView.collection.listModel.set('search', searchHash);
//          
//        },
//        always: function(){
//          console.log('done: ');
//        }
//      }).fail(function(){ 
//        Iccbl.appModel.jqXHRfail.apply(this,arguments); 
//      });        
//    },

    showOtherScreenColumnsDialog: function(listView){
      var self = this;
      
      var searchHash = listView.collection.listModel.get('search');
      var dc_ids = _.result(searchHash,'dc_ids', '');
      if (!_.isArray(dc_ids)){
        dc_ids = dc_ids.split(',')
      }
      dc_ids = _.map(dc_ids, function(dc_id){ return parseInt(dc_id); });
      console.log('showOtherScreenColumnsDialog', dc_ids);
      var resource = appModel.getResource('datacolumn');
      var url = [resource.apiUri, 
                 'for_screen', self.model.get('facility_id')].join('/');
      var data_for_get = { 
        limit: 0,
        includes: [
          'screen_facility_id','screen_title','name','description',
          'assay_data_type','ordinal','study_type'],
        order_by: ['screen_facility_id', 'ordinal'],
        use_vocabularies: false
      };
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: url,
        modelId: function(attrs) {
          return Iccbl.getIdFromIdAttribute( attrs, resource);
        }
      });
      var collection = new CollectionClass();
      collection.fetch({
        data: data_for_get,
        success: function(collection, response) {
          showTreeSelector(collection);
        },
        always: function(){
          console.log('done: ');
        }
      }).fail(function(){ 
        Iccbl.appModel.jqXHRfail.apply(this,arguments); 
      });        

      function showColumns(dcs_selected) {
        var srResource = self.model.resource;
        var fields = srResource.fields;
        
        console.log('show columns', dcs_selected);
        
        _.each(dcs_selected, function(dc){
          var key = dc.get('key');
          var name = dc.get('title');
          var dc_id = dc.get('data_column_id');
          var field = dc.attributes;
          var grid = listView.grid;
          
          if (dc.get('user_access_level_granted') > 1){
            field['filtering'] = true;
            field['ordering'] = true;
          }
          
          if (!_.has(srResource.fields, key)){
            srResource.fields[key] = field;
          }

          var column = grid.columns.findWhere({ name: key });
          if (!column){
            var index = grid.columns.size();
            grid.insertColumn(
              Iccbl.createBackgridColumn(
                  key,field,
                  listView.listModel.get('order')),
                  { at: index});
          } else {
            console.log('column already included', key)
          }
        
        });
        
        var dc_ids_selected = _.map(dcs_selected, function(dcmodel){
          return dcmodel.get('data_column_id');
        });
        
        var searchHash = _.clone(
          listView.collection.listModel.get('search'));
        searchHash['dc_ids'] = dc_ids_selected;
        listView.collection.listModel.set('search', searchHash);
        
      }; // showColumns
      
      
      function showTreeSelector(collection){
        console.log('showTreeSelector', dc_ids);
        collection.each(function(model){
          var dc_id = model.get('data_column_id');
          if (_.contains(dc_ids, dc_id)){
            model.set('checked', true);
          }
          model.set('screen_info',
            model.get('screen_facility_id') + ' - ' + model.get('screen_title'));
        });
        
        var show_positives_control = $([
          '<label class="checkbox-inline pull-left" ',
          '   title="Show positive columns only" >',
          '  <input type="checkbox">positives',
          '</label>'
          ].join(''));
        
        var show_mutual_screens_control = $([
          '<label class="checkbox-inline pull-left" ',
          '   title="Show mutual screens having overlapping positives" >',
          '  <input type="checkbox">mutual screens',
          '</label>'
          ].join(''));
        
        var show_studies_control = $([
          '<label class="checkbox-inline pull-left" ',
          '   title="Show study annotations" >',
          '  <input type="checkbox">show studies',
          '</label>'
          ].join(''));
        
        var OtherColumnsTreeSelector = TreeSelector.extend({
          search: function(){
            var show_pos_only = 
              show_positives_control.find('input[type="checkbox"]').prop('checked');
            var show_mutual_only = 
              show_mutual_screens_control.find('input[type="checkbox"]').prop('checked');
            var show_studies = 
              show_studies_control.find('input[type="checkbox"]').prop('checked');
            var searchedModels;
            if (!_.isEmpty(this.getSearchVal())){
              searchedModels = OtherColumnsTreeSelector.__super__.search.apply(this,arguments);
            } else if (show_pos_only || show_mutual_only || show_studies ) {
              searchedModels = this.collection.models;
            } else{
              this.clearSearch();
              return;
            }
            searchedModels = _.filter(searchedModels, function(model){
              var shown = true;
              if (shown && show_pos_only){
                shown = model.get('positives_count')>0;
              }
              if (shown && show_mutual_only){
                shown = _.contains(
                  self.model.get('overlapping_positive_screens'),
                  model.get('screen_facility_id'));
              }
              if (shown && show_studies){
                shown = !_.isEmpty(model.get('study_type'));
              }
              return shown;
            });
            return searchedModels;
          }
        });
        
        var dcView = new OtherColumnsTreeSelector({
          collection: collection,
          treeAttributes: ['screen_info', 'title'],
          treeAttributesForTypeAhead: ['screen_info'],
          extraControls: [show_positives_control, show_mutual_screens_control, 
                          show_studies_control]
        });

        show_positives_control.find('input[type="checkbox"]').change(function(e) {
          var searchedModels = dcView.search();
          if (this.checked || !_.isEmpty(searchedModels)) {
            collection.trigger('searchChange', searchedModels);
          }
        });
        show_mutual_screens_control.find('input[type="checkbox"]').change(function(e){
          var searchedModels = dcView.search();
          if (this.checked || !_.isEmpty(searchedModels)) {
            collection.trigger('searchChange', searchedModels);
          }
        });
        show_studies_control.find('input[type="checkbox"]').change(function(e){
          var searchedModels = dcView.search();
          if (this.checked || !_.isEmpty(searchedModels)) {
            collection.trigger('searchChange', searchedModels);
          }
        });
        
        Backbone.Layout.setupView(dcView);
  
        var el = dcView.render().el;
        var dialog = appModel.showModal({
            buttons_on_top: true,
            css: { 
                display: 'table',
                height: '500px',
                width: '80%'
              },
            css_modal_content: {
              overflow: 'hidden'
            },
            ok: function(){
              showColumns(dcView.getSelected());
            },
            view: el,
            title: 'Select Other Screen Columns to display'
        });
      }; // showTreeSelector
      
      
    },
    
    /**
     * Loads the screen results using ajax to post the file.
     * - because ajax cannot handle post response attachments (easily), set
     * the response type to 'application/json' and display the errors in a 
     * modal dialog.
     * 
     * NOTE: a "POST" form cannot be used - there is no standard way to signal
     * the result of the post to JavaScript.
     */
    loadScreenResults: function() {
      var self = this;
      var isStudy = self.model.get('project_phase') == 'annotation';
      var form_template = [
         "<div class='form-horizontal container' id='screenresult_form' >",
         "<form data-fieldsets class='form-horizontal container' >",
         "<div class='form-group' ><input type='file' name='fileInput' /></div>",
         "</form>",
         "</div>"].join('');      
      
      var fieldTemplate = _.template([
        '<div class="form-group" >',
        '    <label class="control-label " for="<%= editorId %>"><%= title %></label>',
        '    <div class="" >',
        '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
        '      <div data-error class="text-danger" ></div>',
        '      <div><%= help %></div>',
        '    </div>',
        '  </div>',
      ].join(''));
      
      var formSchema = {};
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        editorClass: 'input-full',
        type: 'TextArea',
        template: fieldTemplate
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs) {
          var errs = {};
          var file = $('input[name="fileInput"]')[0].files[0]; 
          if (! file) {
            errs.fileInput = 'Must specify a file';
          }
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: _.template(form_template)
      });
      var _form_el = form.render().el;

      function saveFile() {

        var errors = form.commit({ validate: true }); // runs schema and model validation
        if (!_.isEmpty(errors)) {
          console.log('form errors, abort submit: ',errors);
          _.each(_.keys(errors), function(key) {
            $('[name="'+key +'"').parents('.form-group').addClass('has-error');
            if ( key == '_others') {
              _.each(errors[key], function(otherError) {
                // TODO: display file req'd msg
              });
            }
          });
          return false;
        } else {
          var url = [self.model.resource.apiUri,self.model.key,'screenresult'].join('/');
          var values = form.getValue();
          var file = $('input[name="fileInput"]')[0].files[0]; 
          var data = new FormData();
          // NOTE: 'xls' is sent as the key to the FILE object in the upload.
          // Use this as a non-standard way to signal the upload type to the 
          // serializer. 
          data.append('xls', file);
          var headers = {
            'Accept': 'application/json'
          };
          headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
          $.ajax({
            url: url,    
            data: data,
            cache: false,
            contentType: false,
            processData: false,
            type: 'POST',
            headers: headers
          })
          .fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); })
          .done(function(model, resp) {
            // FIXME: should be showing a regular message
            appModel.error('success');
            // Must refetch and render, once for the containing tabbed layout,
            // and then (during render), will refetch again for the table view
            self.model.fetch({
              success: function() {
                if (!isStudy) {
                  self.uriStack = ['detail'];
                  // remove the child view before calling render, to prevent
                  // it from being rendered twice, and calling afterRender twice
                  self.removeView('#tab_container');
                }
                self.render();
              }
            }).fail(function() { 
              Iccbl.appModel.jqXHRfail.apply(this,arguments); 
            });
          });
        }        
      }
      
      var dialog = appModel.showModal({
          ok: saveFile,
          view: _form_el,
          title: 'Select a Screen Result (xlsx workbook) file to upload'
      });
    },
    
    deleteScreenResults: function() {
      var self = this;
      var isStudy = self.model.get('project_phase') == 'annotation';
      var formSchema = {};
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        editorClass: 'input-full',
        validators: ['required'],
        type: 'TextArea'
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields
      });
      var _form_el = form.render().el;
      var dialog = appModel.showModal({
          ok: function() {
            var errors = form.commit({ validate: true }); 
            if (!_.isEmpty(errors)) {
              _.each(_.keys(errors), function(key) {
                $('[name="'+key +'"').parents('.form-group').addClass('has-error');
              });
              return false;
            } else {
              var url = [self.model.resource.apiUri,self.model.key,'screenresult'].join('/');
              var values = form.getValue();
              var headers = {
                'Accept': 'application/json'
              };
              headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
              $.ajax({
                url: url,    
                cache: false,
                contentType: false,
                processData: false,
                type: 'DELETE',
                headers: headers
              })
              .fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); })
              .done(function(model, resp) {
                appModel.error('success');
                // Must refetch and render, once for the containing tabbed layout,
                // and then (during render), will refetch again for the detail view
                self.model.fetch({
                  success: function() {
                    if (!isStudy) {
                      self.uriStack = ['detail'];
                      // remove the child view before calling render, to prevent
                      // it from being rendered twice
                      self.removeView('#tab_container');
                    }
                    self.render();
                  }
                }).fail(function() { 
                  Iccbl.appModel.jqXHRfail.apply(this,arguments); 
                });
              });
            }        
          },
          view: _form_el,
          title: 'Delete screen results for ' + self.model.get('facility_id')
      });
      
    },
    
    setDatacolumns: function(delegateStack) {
      
      var self = this;
      var datacolumnResource = appModel.getResource('datacolumn'); 
      var url = [self.model.resource.apiUri,self.model.key,'datacolumns'].join('/');
      
      // construct a collection-in-columns
      // pivot the datacolumn list
      Iccbl.getCollectionOnClient(url, function(collection) {
        // create a colModel for the list
        var columns = [];
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell.extend({
            render: function() {
              this.$el.empty();
              var label = Iccbl.createLabel(
                this.column.get("label"), 10);
              this.$el.append(label);
              return this;
            }
          })
        };
        // setup the first column - schema field titles (that are visible)
        var col0 = _.extend({},colTemplate,{
          'name' : 'fieldname',
          'label' : '',
          'description' : 'Datacolumn field',
          'cell': Iccbl.TextWrapCell.extend({
            className: '150_width'
          })
        });
        columns.push(col0);
        
        collection.each(function(datacolumn) {
          var col = _.extend({},colTemplate,{
            'name' : datacolumn.get('key'),
            'label' : datacolumn.get('title'),
            'description' : datacolumn.get('description'),
            'order': datacolumn.get('ordinal'),
            'sortable': true,
            'cell': Iccbl.TextWrapCell
          });
          columns.push(col);
          if (datacolumn.has('derived_from_columns')) {
            datacolumn.set(
              'derived_from_columns', 
              datacolumn.get('derived_from_columns').join(', '));
          }
        });
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();
        
        // create the collection by pivoting
        orderedFields = _.sortBy(datacolumnResource.fields,'ordinal');
        var pivotCollection = new Backbone.Collection();
        _.each(orderedFields, function(field) {
          if (_.contains(field.visibility, 'l') 
              && !_.contains(['name','title','key'],field.key )) {
            var row = {'key': field.key, 'fieldname': field.title };
            collection.each(function(datacolumn) {
              row[datacolumn.get('key')] = datacolumn.get(field.key);
            });
            pivotCollection.push(row);
          }
        });
        
        // now show as a list view
        var view = new Backgrid.Grid({
          columns: colModel,
          collection: pivotCollection,
          className: 'backgrid table-striped table-condensed table-hover'
        });
        var view = view.render();
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);

        self.setView("#tab_container", view ).render();
     });
      
    }    

  });
  
  return ScreenDataView;
});