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
      var url = _.result(self.model,'url');
     
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
        
        overrides_granted: {},

        save_success: function(data, textStatus, jqXHR){
          var inner_self = this;
          var meta = _.result(data, 'meta', null);
          if (meta) {
            appModel.showJsonMessages(meta);
          }
          var objects = _.result(data, 'objects', null)
          if (objects && objects.length == 1) {
            model = new Backbone.Model(objects[0]);
            var key = Iccbl.getIdFromIdAttribute( model,self.model.resource );
            appModel.router.navigate([
              self.screen.resource.key,self.screen.key,'summary/libraryscreening',key].join('/'), 
              {trigger:true});
          } else {
            console.log('no objects in the response', data);
            appModel.error('Could not display the server response');
          }
        },
        
        save: function(changedAttributes, options){
          var inner_self = this;
          var options = options || { headers: {} };
          var comment = _.result(options.headers,appModel.HEADER_APILOG_COMMENT);
          this.model.save(changedAttributes, options)
            .done(function(data, textStatus, jqXHR){ 
              inner_self.save_success.apply(this,arguments);
            })
            .fail(function(jqXHR, textStatus, errorThrown) { 
              /// "not allowed" libraries require an "override" param
              /// see cherryPickRequest.js
              console.log('errors', arguments);
              
              var jsonError = _.result(jqXHR, 'responseJSON');
              if (!_.isUndefined(jsonError)){
                var error = _.result(jsonError, 'errors');
                var libraryErrorFlag = _.result(error,appModel.API_PARAM_OVERRIDE);
                var volumeErrorFlag = _.result(error,appModel.API_PARAM_VOLUME_OVERRIDE)
                if (!_.isUndefined(libraryErrorFlag)){
                  var title = _.result(error, 'library_plates_screened');
                  if (!title){
                    console.log('Error: expecting "library_plates_screened" error message');
                    title = 'Override required for libraries that are not "allowed"';
                  }
                  var bodyMsg = 'Not allowed libraries';
                  var errors = _.result(error, 'Libraries');
                  if (!errors){
                    console.log('Libraries should be reported in the response');
                  }else{
                    bodyMsg = 'Libraries: ' + errors.join(', ');
                  }
                  appModel.showOkCommentForm({
                    title: title,
                    body: bodyMsg,
                    okText: 'Confirm Override',
                    ok: function(formValues) {
                      inner_self.overrides_granted[appModel.API_PARAM_OVERRIDE] = 'true';
                      self.model.url = function(){
                        var new_url = url + '?' + _.map(_.pairs(inner_self.overrides_granted),
                          function(keyval){
                            return keyval.join('=');
                          }).join('&')
                        console.log('new url', new_url);
                        return new_url;
                      };
                      var new_comment = formValues['comments']
                      if (new_comment){
                        if (comment){
                          comment += '; ' + new_comment;
                        } else {
                          comment = new_comment
                        }
                        options.headers[appModel.HEADER_APILOG_COMMENT] = comment;
                        console.log('new comments', comment);
                      }
                      inner_self.save(changedAttributes, options);
                    }
                  });
                } else if (!_.isUndefined(volumeErrorFlag)){
                  var title = _.result(error, 'library_plates_screened');
                  if (!title){
                    console.log('Error: expecting "library_plates_screened" error message');
                    title = 'Override required to screen plates with insufficient volume';
                  }
                  var bodyMsg = 'Insufficient volume for plates (unspecified)';
                  var errors = _.result(error, 'Plates');
                  if (!errors){
                    console.log('Plates should be reported in the response');
                  }else{
                    bodyMsg = 'Plates: ' + errors.join(', ');
                  }
                  appModel.showOkCommentForm({
                    title: title,
                    body: bodyMsg,
                    okText: 'Confirm Override',
                    ok: function(formValues) {
                      inner_self.overrides_granted[appModel.API_PARAM_VOLUME_OVERRIDE] = 'true';
                      self.model.url = function(){
                        var new_url = url + '?' + _.map(_.pairs(inner_self.overrides_granted),
                          function(keyval){
                            return keyval.join('=');
                          }).join('&')
                        console.log('new url', new_url);
                        return new_url;
                      };
                      var new_comment = formValues['comments']
                      if (new_comment){
                        if (comment){
                          comment += '; ' + new_comment;
                        } else {
                          comment = new_comment
                        }
                        options.headers[appModel.HEADER_APILOG_COMMENT] = comment;
                        console.log('new comments', comment);
                      }
                      inner_self.save(changedAttributes, options);
                    }
                  });
                } else {
                  inner_self.save_fail.apply(this,arguments);
                }
              } else {
                inner_self.save_fail.apply(this,arguments);
              }
            });
        },
        
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
            appModel.getScreeningLibraryOptions(
              self.model.get('screen_type'),
              function(libraryOptions){
                libraryOptions.unshift({val:'',label:''});
                self._addPlateRangeDialog(
                  plate_collection,libraryOptions,nested_library_plate_pattern,
                  for_screening=true);
              });
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
            var vol_to_assayplates  = self_editform.getValue(
              'volume_transferred_per_well_to_assay_plates');
            if(_.isNumber(num_replicates) && _.isNumber(vol_to_assayplates)){
              if(num_replicates > 0 && vol_to_assayplates > 0){
                var calculated = num_replicates * vol_to_assayplates;
                self_editform.setValue(
                  'volume_transferred_per_well_from_library_plates',calculated);
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
          $target_el = self_editform.$el.find(
            '[key=form-group-volume_transferred_per_well_from_library_plates]');
          $target_el.append(calcVolButton);
          calcVolButton.click(function(event){
            event.preventDefault();
            calculateVolFromLibraryPlates();
          });
          this.listenTo(
            this, "number_of_replicates:change", calculateVolFromLibraryPlates);
          this.listenTo(
            this, "volume_transferred_per_well_to_assay_plates:change", 
            calculateVolFromLibraryPlates);
        }
      });
      
      args.EditView = editView;
      // force all fields to be sent to server (including library_plates_screened)
      // even if unchanged
      args.fullSaveOnEdit = true; 
      args.DetailView = detailView;
//      args.saveSuccessCallBack = function(model){
//        var meta = _.result(model, 'meta', null);
//        if (meta) {
//          appModel.showJsonMessages(meta);
//        }
//        var objects = _.result(model, 'objects', null)
//        if (objects && objects.length == 1) {
//          model = new Backbone.Model(objects[0]);
//          var key = Iccbl.getIdFromIdAttribute( model,self.model.resource );
//          model.key = self.model.resource.key + '/' + key;
//          appModel.router.navigate([
//            self.screen.resource.key,self.screen.key,'summary/libraryscreening',key].join('/'), 
//            {trigger:true});
//        } else {
//          console.log('no objects in the response', model);
//          appModel.error('Could not display the server response');
//        }
//      };
      
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
//        var userOptions = appModel.getAdminUserOptions();
//        fields['performed_by_username']['choices'] = (
//            [{ val: '', label: ''}].concat(userOptions));
        
        // TODO: replace with appModel.get_screen_members
        var fields = self.model.resource.fields;
        var users = appModel.get('users');
        var userOptions = [
          { val: self.screen.get('lead_screener_username'), 
            label: self.screen.get('lead_screener_name') },
          { val: self.screen.get('lab_head_username'), 
            label: self.screen.get('lab_name') },
        ];
        userOptions = userOptions.concat(_.map(_.zip(
          self.screen.get('collaborator_usernames'),
          self.screen.get('collaborator_names')),
          function(pair){
          return { val: pair[0], label: pair[1] };
        }));
        fields['performed_by_username']['choices'] = (
            [{ val: '', label: ''}].concat(userOptions));
        
        fields['screen_facility_id']['editability'] = [];
        DetailLayout.prototype.showEdit.apply(self,arguments);
      });  
    }
    
  });

  return LibraryScreeningView;
});