define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_detail_layout',
  'views/library/library',
  'utils/treeSelector',
  'templates/genericResource.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         LibraryView, TreeSelector, layout) {
  
  var LibraryWellsView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      this.args = args;
      this.resource = appModel.cloneResource(args.resource);
      this.uriStack = args.uriStack;
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    afterRender: function(){
      this.showList(this.uriStack);
    },     
    
    showDataColumnsDialog: function(listView){
      var self = this;
      
      var searchHash = listView.collection.listModel.get(appModel.URI_PATH_SEARCH);
      var dc_ids = _.result(searchHash,'dc_ids', '');
      if (!_.isArray(dc_ids)){
        dc_ids = dc_ids.split(',')
      }
      dc_ids = _.map(dc_ids, function(dc_id){ return parseInt(dc_id); });
      console.log('showDataColumnsDialog', dc_ids);
      var resource = appModel.getResource('datacolumn');
      var data_for_get = { 
        limit: 0,
        // access level_1 - only granted for level 2 screens that are overlapping,
        // - or for level 2 users viewing overlapping level 1 screens; this 
        // class of access is not available on well search because must only be
        // used in the context of a screen result
        // access level 2 - shared mutually, and public screens
        // access level 3 - own screens, or if current user is an admin
        user_access_level_granted__gt: 1, 
        includes: [
          'screen_facility_id','screen_title','name','description',
          'assay_data_type','ordinal','study_type'],
        order_by: ['screen_facility_id', 'ordinal'],
        use_vocabularies: false
      };
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: resource.apiUri,
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
        var fields = self.resource.fields;
        var grid = listView.grid;
        
        console.log('show columns', dcs_selected);
        
        _.each(dcs_selected, function(dc){
          var key = dc.get('key');
          var name = dc.get('title');
          var dc_id = dc.get('data_column_id');
          var field = dc.attributes;
          
          if (dc.get('user_access_level_granted') > 1){
            field['filtering'] = true;
            field['ordering'] = true;
          } else {
            throw 'data column should not appear in this search: ' + key
              + ': user_access_level_granted:' + dc.get('user_access_level_granted'); 
          }
          
          if (!_.has(fields, key)){
            fields[key] = field;
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
        // remove unselected
        _.each(dc_ids, function(former_dc_id){
          if (!_.contains(dc_ids_selected,former_dc_id)){
            var column = grid.columns.find(function(gridCol){
              var fi = gridCol.get('fieldinformation');
              if (fi){
                return fi['data_column_id'] == former_dc_id;
              }
            });
            if (!column){
              console.log('column already not present', former_dc_id)
            } else {
              grid.removeColumn(column);
            }
          }
        });
        var searchHash = _.clone(
          listView.collection.listModel.get(appModel.URI_PATH_SEARCH));
        searchHash['dc_ids'] = dc_ids_selected;
        listView.collection.listModel.set(appModel.URI_PATH_SEARCH, searchHash);
        
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
            var show_studies = 
              show_studies_control.find('input[type="checkbox"]').prop('checked');
            var searchedModels;
            if (!_.isEmpty(this.getSearchVal())){
              searchedModels = OtherColumnsTreeSelector.__super__.search.apply(this,arguments);
            } else if (show_pos_only || show_studies ) {
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
          extraControls: [show_positives_control, show_studies_control]
        });

        show_positives_control.find('input[type="checkbox"]').change(function(e) {
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
            title: 'Select Screen and Study Columns to display'
        });
      }; // showTreeSelector
      
      
    },

    showList: function(delegateStack) {
      var self = this;

      function createReagentView(schemaResult) {
        // Update the current resource and fields:
        // - current resource prop and field definitions override the new 
        // definitions, if defined, but the new fields are used if not.
        var newFields = _.extend({}, schemaResult['fields'], self.resource['fields'])
        var newResource = _.extend({}, schemaResult, self.resource );
        
        newResource['fields'] = newFields;
        var url = newResource.apiUri;
        if (newResource.key == 'compound_search'){
          url += '/compound_search';
          // FIXME: hack to add columns; fix is to implement sirna/smr resource
          // schema as superset of reagent schema
          newResource['options']['includes'] = [
           'inchi','smiles','structure_image','molecular_formula',
           'molecular_mass','molecular_weight','pubchem_cid','chembank_id',
           'chembl_id'];
        }
        
        var extraControls = [];
        var show_data_columns_control = $([
          '<button class="btn btn-default btn-sm pull-right" role="button" ',
          'id="show_data_columns_control" title="Show study and screen data columns" >',
          'Select Study and Screen Columns',
          '</button>'
          ].join(''));
        show_data_columns_control.click(function(e){
          self.showDataColumnsDialog(view);
        });
        
        extraControls.push(show_data_columns_control);
        var WellListView = ListView.extend({
          afterRender: function(){
            ListView.prototype.afterRender.apply(this,arguments);
            this.$('#list_controls').append(show_data_columns_control);
            return this;
          }
        });
        
        var viewArgs = _.extend({},self.args,{
          resource: newResource,
          url: url});

        var view = new WellListView (viewArgs);
        self.listenTo(view, 'update_title', function(val) {
          if(val) {
            self.$('#content_title').html(newResource.title + ': <small>' + val + '</small>');
            $('#content_title_row').show();
          } else {
            self.$('#content_title').html(newResource.title);
            $('#content_title_row').show();
          }
        });
          
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView('#resource_content', view ).render();
        
      }
      
      // Parse the delegateStack to determine if extra data columns are specified
      var listOptions = ListView.prototype.parseUriStack(delegateStack);
      console.log('parsed listOptions', listOptions);
      var dc_ids = [];
      if (_.has(listOptions, appModel.URI_PATH_SEARCH)){
        if (_.has(listOptions.search,'dc_ids')){
          dc_ids = listOptions.search.dc_ids;
          if (dc_ids.length > 0 && _.isString(dc_ids)){
            dc_ids = dc_ids.split(',');
          }
        }
      }
      dc_ids = _.map(dc_ids, function(dc_id){ return parseInt(dc_id); });
      console.log('dc_ids', dc_ids);
      if (!_.isEmpty(dc_ids)){
        // Re-fetch the schema if extra data columns are specified
        var options = {};
        options['dc_ids'] = dc_ids;
        var schemaUrl = [self.resource.apiUri,'schema'].join('/');
        appModel.getResourceFromUrl(schemaUrl, createReagentView, options);
      } else {
        var schemaUrl = [self.resource.apiUri,'schema'].join('/');

        appModel.getResourceFromUrl(schemaUrl, createReagentView, options);
//        createReagentView(self.resource);
      }
      
    }
    
  });

  return LibraryWellsView;
});