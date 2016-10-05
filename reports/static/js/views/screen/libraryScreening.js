define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'utils/plateRangeTable',
  'views/generic_detail_stickit',
  'views/generic_edit',
  'views/generic_detail_layout'
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, 
         PlateRangePrototype, DetailView, EditView, DetailLayout) {
  
  var nested_library_plate_pattern = '{library}:{copy}:{start_plate}-{end_plate}';
  
  var Prototype = DetailLayout.extend(PlateRangePrototype);
  
  var LibraryScreeningView = Prototype.extend({

    initialize: function(args) {
      
      var self = this;
      this.screen = args.screen;
      
      var fields = self.model.resource.fields;
      var library_plate_parser = null; 
      if (!_.has(fields,'library_plates_screened')
          || !_.has(fields['library_plates_screened'],'regex')){
        console.log('error....field schema missing/regex for library_plates_screened',
          self, fields);
        appModel.error('field schema missing/regex for library_plates_screened');
      }else{
        var lcp_regex_string = fields['library_plates_screened']['regex'];
        try{
          lcp_regex = RegExp(lcp_regex_string);
          library_plate_parser = _.partial(self._parse_library_plate_entry,lcp_regex);
        }catch(e){
          appModel.error('regex misconfigured for "library_plates_screened" in metadata: ' + lcp_regex_string);
        }
      }
     
      var detailView = DetailView.extend({
        
        afterRender: function(){
          
          DetailView.prototype.afterRender.apply(this,arguments);
          
          var $target_el = this.$el.find('#library_plates_screened');
          var plate_collection = new Backbone.Collection;
          plate_collection.comparator = 'start_plate';
          plate_collection.set(
            _.map(self.model.get('library_plates_screened'),library_plate_parser));
          self._createPlateRangeTable(plate_collection, $target_el, nested_library_plate_pattern);
          
          if(!_.isEmpty(self.model.get('assay_protocol'))){
            $target_el = this.$el.find('#assay_protocol_last_modified_date');
            var apilogResource = appModel.getResource('apilog');
            var options = {
              data_for_get: {
                limit: 1,
                key: self.model.get('activity_id'),
                ref_resource_name: self.model.resource.key,
                diff_keys__icontains: '"assay_protocol"',
                order_by: ['-date_time']
              }
            };
            Iccbl.getCollectionOnClient(apilogResource.apiUri,function(collection){
              if(collection && !collection.isEmpty()){
                var model = collection.at(0);
                console.log('assay_protocol_last_modified_date',model.get('date_time'));
                self.model.set('assay_protocol_last_modified_date',model.get('date_time'));
              }else{
                // no apilog found, set it to the creation date
                self.model.set('assay_protocol_last_modified_date',
                  self.model.get('date_created'));
              }
            },options);
          }
        
        }
      });
      
      var editView = EditView.extend({
        
        afterRender: function(){

          EditView.prototype.afterRender.apply(this,arguments);

          var self_editform = this;
          var $target_el = this.$el.find('[name="library_plates_screened"]');
          var plate_collection = new Backbone.Collection;
          plate_collection.comparator = 'start_plate';
          self.listenTo(plate_collection, "MyCollection:delete", function (model) {
            var entry = Iccbl.formatString(nested_library_plate_pattern,model);
            plate_collection.remove(model);
          });

          self_editform.listenTo(plate_collection, 'update', function (plate_collection) {
            self_editform.setValue('library_plates_screened', plate_collection.map(function(model){
              return Iccbl.formatString(nested_library_plate_pattern, model);
            }));
          });
          var lps = self.model.get('library_plates_screened');
          if(_.isEmpty(lps)) lps = [];
          plate_collection.set(_.map(lps,library_plate_parser));
          self._createPlateRangeTable(
            plate_collection, $target_el, nested_library_plate_pattern, true);

          // library_copy_plates button
          var addButton = $([
            '<button class="btn btn-default btn-sm" ',
              'role="button" id="add_library_plates_screened_button" href="#">',
              'Add</button>'
            ].join(''));
          addButton.click(function(event){
            event.preventDefault();
            libraryOptions = appModel.getScreeningLibraryOptions(self.model.get('screen_type'));
            libraryOptions.unshift({val:'',label:''});
            self._addPlateRangeDialog(
              plate_collection,libraryOptions,nested_library_plate_pattern,
              for_screening=true);
          });
          $target_el = $target_el.parent();
          $target_el.append(addButton);
          
          // attach is_for_external_library_plates change listener
          $target_el = self_editform.$el.find('[key=form-group-library_plates_screened]');
          this.listenTo(this, "is_for_external_library_plates:change", function(e){
            var val = self_editform.getValue('is_for_external_library_plates');
            console.log('is_for_external_library_plates:',val)
            if(val){
              $target_el.hide();
              self_editform.setValue('library_plates_screened',null);
            } else {
              $target_el.show();
            }
          });
          // attach listener volume calculator
          function calculateVolFromLibraryPlates(){
            var num_replicates = self_editform.getValue('number_of_replicates');
            var vol_to_assayplates  = self_editform.getValue('volume_transferred_per_well_to_assay_plates');
            if(_.isNumber(num_replicates) && _.isNumber(vol_to_assayplates)){
              if(num_replicates > 0 && vol_to_assayplates > 0){
                var calculated = num_replicates * vol_to_assayplates;
                self_editform.setValue('volume_transferred_per_well_from_library_plates',calculated);
              }
            } else {
              console.log('cannot calculate volume until form values are entered');
            }
          };
          var calcVolButton = $([
            '<button class="btn btn-default btn-sm" ',
              'role="button" id="calc_volume_transferred_per_well_from_library_plates" href="#">',
              'calculate</button>'
            ].join(''));
          $target_el = self_editform.$el.find('[key=form-group-volume_transferred_per_well_from_library_plates]');
          $target_el.append(calcVolButton);
          calcVolButton.click(function(event){
            event.preventDefault();
            calculateVolFromLibraryPlates();
          });
          this.listenTo(this, "number_of_replicates:change", calculateVolFromLibraryPlates);
          this.listenTo(this, "volume_transferred_per_well_to_assay_plates:change", calculateVolFromLibraryPlates);
        }
      });
      args.EditView = editView;
      args.DetailView = detailView;
      args.saveSuccessCallBack = function(model){
        model = new Backbone.Model(model);
        var key = Iccbl.getIdFromIdAttribute( model,self.model.resource );
        model.key = self.model.resource.key + '/' + key;
        appModel.router.navigate([
          self.screen.resource.key,self.screen.key,'summary/libraryscreening',key].join('/'), 
          {trigger:true});
      };
      
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
        var fields = self.model.resource.fields;
        fields['performed_by_username']['choices'] = (
            [{ val: '', label: ''}].concat(userOptions));
        fields['screen_facility_id']['editability'] = [];
        DetailLayout.prototype.showEdit.apply(self,arguments);
      });  
    }
    
  });

  return LibraryScreeningView;
});