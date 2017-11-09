define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'layoutmanager',
  'bootstrap-3-typeahead',
  'iccbl_backgrid',
  'models/app_state',
  'utils/plateRangeTable',
  'views/generic_detail_stickit',
  'views/generic_edit',
  'views/generic_detail_layout'
], 
function($, _, Backbone, Backgrid, layoutmanager, typeahead, Iccbl, appModel, 
         PlateRangePrototype, DetailView, EditView, DetailLayout) {

  var nested_library_plate_pattern = '{library}:{copy}:{start_plate}-{end_plate}';
  
  var Prototype = DetailLayout.extend(PlateRangePrototype);
  var PlateLocationView = Prototype.extend({
    
    initialize: function(args) {
      
      var self = this;
      var fields = self.model.resource.fields;
      var library_plate_parser = null; 
      if (!_.has(fields,'copy_plate_ranges')
          || !_.has(fields['copy_plate_ranges'],'regex')){
        console.log('error....field schema missing/regex for copy_plate_ranges',
          self, fields);
        appModel.error('field schema missing/regex for copy_plate_ranges');
      }else{
        var lcp_regex_string = fields['copy_plate_ranges']['regex'];
        try{
          lcp_regex = RegExp(lcp_regex_string);
          library_plate_parser = _.partial(self._parse_library_plate_entry,lcp_regex);
        }catch(e){
          appModel.error(
            'regex misconfigured for "copy_plate_ranges" in metadata: ' 
            + lcp_regex_string);
        }
      }
      
      var detailView = DetailView.extend( {
        
        afterRender: function(){
          
          DetailView.prototype.afterRender.apply(this,arguments);

          var $target_el = this.$el.find('#copy_plate_ranges');
          var plate_collection = new Backbone.Collection;
          plate_collection.comparator = function(platerange){
            return parseInt(platerange.get("start_plate"));
          };
          plate_collection.set(
            _.map(self.model.get('copy_plate_ranges'),library_plate_parser));
          self._createPlateRangeTable(
            plate_collection, $target_el, nested_library_plate_pattern);
          
          $target_el = this.$el.find('#bin_last_modified_date');
          var apilogResource = appModel.getResource('apilog');
          var options = {
            data_for_get: {
              limit: 1,
              key: Iccbl.getIdFromIdAttribute(self.model, self.model.resource),
              ref_resource_name: self.model.resource.key,
              diff_keys__icontains: '"copy_plate_ranges"',
              order_by: ['-date_time']
            }
          };
          Iccbl.getCollectionOnClient(
            apilogResource.apiUri,
            function(collection){
              if(collection && !collection.isEmpty()){
                var model = collection.at(0);
                // TODO: bin last modified date?
              }else{
                // no apilog found, set it to the creation date
                //self.model.set('assay_protocol_last_modified_date',
                //  self.model.get('date_created'));
              }
            },options);
        }
      });
      
      var editView = EditView.extend({
        
        afterRender: function(){

          EditView.prototype.afterRender.apply(this,arguments);

          var self_editform = this;
          var $target_el = this.$el.find('[name="copy_plate_ranges"]');
          var plate_collection = new Backbone.Collection;
          plate_collection.comparator = function(platerange){
            return parseInt(platerange.get("start_plate"));
          };
          self.listenTo(plate_collection, "MyCollection:delete", function (model) {
            var entry = Iccbl.formatString(nested_library_plate_pattern,model);
            plate_collection.remove(model);
          });

          self_editform.listenTo(plate_collection, 'update', function (plate_collection) {
            self_editform.setValue('copy_plate_ranges', plate_collection.map(function(model){
              return Iccbl.formatString(nested_library_plate_pattern, model);
            }));
          });
          var lps = self.model.get('copy_plate_ranges');
          if(_.isEmpty(lps)) lps = [];
          plate_collection.set(_.map(lps,library_plate_parser));
          self._createPlateRangeTable(
            plate_collection, $target_el, nested_library_plate_pattern, true);

          // library_copy_plates button
          var addButton = $([
            '<button class="btn btn-default btn-sm" ',
              'role="button" id="add_copy_plate_ranges_button" href="#">',
              'Add</button>'
            ].join(''));
          addButton.click(function(event){
            event.preventDefault();
            libraryOptions = appModel.getLibraryOptions();
            libraryOptions.unshift({val:'',label:''});
            self._addPlateRangeDialog(
               plate_collection,libraryOptions,nested_library_plate_pattern);
          });
          $target_el = $target_el.parent();
          $target_el.append(addButton);
          
          var plateLocationTree = appModel.getPlateLocationTree();
          console.log('construct the platelocations, ', plateLocationTree );

          var subKey = 'room';
          $('[name="'+subKey +'"]').typeahead({
            autoSelect: false,
            delay: 1,
            minLength: 0,
            items: 'all',
            source: _.keys(plateLocationTree),
            afterSelect: function(val){
              var subKey = 'freezer';
              view.setValue(subKey,null);
              view.setValue('shelf',null);
              view.setValue('bin',null);
              $('[name="'+subKey +'"]').typeahead('destroy')
              $('[name="'+subKey +'"]').typeahead({
                autoSelect: false,
                delay: 1,
                minLength: 0,
                items: 'all',
                source: _.keys(plateLocationTree[val]),
                afterSelect: function(freezerVal){
                  var subKey = 'shelf';
                  view.setValue(subKey,null);
                  view.setValue('bin',null);
                  $('[name="'+subKey +'"]').typeahead('destroy')
                  $('[name="'+subKey +'"]').typeahead({
                    autoSelect: false,
                    delay: 1,
                    minLength: 0,
                    items: 'all',
                    source: _.keys(plateLocationTree[val][freezerVal]),
                    afterSelect: function(shelfVal){
                      var subKey = 'bin';
                      view.setValue(subKey,null);
                      $('[name="'+subKey +'"]').typeahead('destroy')
                      $('[name="'+subKey +'"]').typeahead({
                        autoSelect: false,
                        delay: 1,
                        minLength: 0,
                        items: 'all',
                        source: _.keys(plateLocationTree[val][freezerVal][shelfVal])
                      });// bin
                    }
                  });// shelf
                }
              });//freezer
            }
          });//room
        }//editView.afterRender
      });//editView
      args.EditView = editView;
      args.DetailView = detailView;
      
      DetailLayout.prototype.initialize.apply(this,args);

    },
    
    afterRender: function() {
      DetailLayout.prototype.afterRender.apply(this,arguments);
      //$title = this.$el.find('#content_title');
      //$title.html('Plate Location');
      //$title.parent().show();
    },    
    
    showEdit: function() {
      var self = this;
      var fields = self.model.resource.fields;
      var initfun = function(){
        view = DetailLayout.prototype.showEdit.apply(self,arguments);
      };
      $(this).queue([appModel.getPlateLocationTree,initfun]);
    }
  });

  return PlateLocationView;
});