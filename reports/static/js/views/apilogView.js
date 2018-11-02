define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_detail_stickit',
  'views/generic_edit',
  'views/generic_detail_layout',
  'templates/generic-detail-stickit.html',
  'templates/generic-detail-stickit-one-column.html',
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, 
         DetailView, EditView, DetailLayout, detailTemplate, detailTemplateOneColumn) {

  var ApilogView = DetailLayout.extend({
    
    initialize: function(args) {
      
      console.log('ApilogView', args);
      var self = this;
      self._args = args;
      
      var detailView = DetailView.extend( {
        
        afterRender: function(){
          DetailView.prototype.afterRender.apply(this,arguments);
          
          // Retrieve the schema for the ref_resource_name
          self.get_reference_schema(self.model, function(refSchema){
            var fields = refSchema.fields;
            
            // Replace the diffs element with a diff table
            var $target_el = $('#diffs');
            self.createDiffTable($target_el, refSchema);
            
            var diff_keys = self.model.get('diff_keys');
            diff_keys = _.map(diff_keys, function(key){
              return fields[key]['title'];
            }).join(', ');
            $('#diff_keys').html(diff_keys);
            
            // Format the json_field
            var json_field = self.model.get('json_field');
            if (json_field){
              json_field = appModel.print_json(json_field);
              // insert line breaks
              json_field = json_field.replace(/^(.?)/gm,'<br/>$1');
              // insert hard spaces
              json_field = json_field.replace(/\s/gm,'&nbsp;');
              console.log('json_field', json_field);
              $('#json_field').html(json_field);
            }
            
          });
        },
        template: _.template(detailTemplateOneColumn)
        
      }, self._args);
      
      args.DetailView = detailView;
      
      DetailLayout.prototype.initialize.call(this,args);
      
      _.bindAll(this,'createDiffTable');
    },

    get_reference_schema: function(model, callBack){
      var self = this;
      var ref_resource_name = self.model.get('ref_resource_name');
      
      // If viewing a well based resource, fetch the specific schema
      if (_.contains(['well','reagent'], ref_resource_name)){
        var schemaUrl;
        var parent_log_uri = model.get('parent_log_uri');
        console.log('parent_log_uri', parent_log_uri);
        if (parent_log_uri.indexOf('library') > -1){
          var library_short_name = parent_log_uri.split('/')[1];
          schemaUrl = [appModel.dbApiUri, ref_resource_name, 'schema'].join('/');
          schemaUrl += '?library_short_name=' + library_short_name;
        }
        if (schemaUrl){
          var options = {};
          appModel.getResourceFromUrl(schemaUrl, function(newSchema){
            callBack(newSchema);
          }, options);
        } else {
          // fall back to stock resource schema
          callBack(appModel.getResource(ref_resource_name));
        }
      } else {
        callBack(appModel.getResource(ref_resource_name));
      }
      
    },
    
    /**
     * Update the diffs display with a table
     **/
    createDiffTable: function($target_el, refSchema) {
      var self = this;
      var fields = refSchema.fields;
      var diffs = self.model.get('diffs');
      
      $target_el.empty();
      
      if (_.isEmpty(diffs)){
        return;
      }
      // build a collection out of the diffs
      var collection = new Backbone.Collection();
      _.each(_.keys(diffs), function(key){
        var fieldName = key;
        if (_.has(fields, key)){
          fieldName = fields[key]['title'];
        }
        var diffValues = diffs[key];
        collection.add(new Backbone.Model({
          key: fieldName,
          before: diffValues[0],
          after: diffValues[1]
        }));
      });
    
      console.log('build diff table', collection);
      if (collection.isEmpty()) {
        return;
      }
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
          'name' : 'key',
          'label' : 'Field',
          'description' : 'Field Name',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.StringCell
        }),
        _.extend({},colTemplate,{
          'name' : 'before',
          'label' : 'Before',
          'description' : 'Value before update',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.StringCell
        }),
        _.extend({},colTemplate,{
          'name' : 'after',
          'label' : 'After',
          'description' : 'Value after update',
          'order': 1,
          'sortable': true,
          'cell': Iccbl.StringCell
        })
      ];
      var colModel = new Backgrid.Columns(columns);
      colModel.comparator = 'ordinal';
      colModel.sort();

      var diff_grid = new Backgrid.Grid({
        columns: colModel,
        collection: collection,
        className: 'backgrid table-striped table-condensed table-hover '
      });
      $target_el.html(diff_grid.render().$el);
      
    },
    
    
    afterRender: function() {
      DetailLayout.prototype.afterRender.apply(this,arguments);
    },    
    
  });

  return ApilogView;
});