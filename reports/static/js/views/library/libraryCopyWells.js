define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_detail_layout',
  'views/generic_edit',
  'views/library/library',
  'utils/uploadDataForm',
  'templates/genericResource.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         EditView, LibraryView, UploadDataForm, layout) {
  
  var LibraryCopyWellView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      this.uriStack = args.uriStack;
      this.library = args.library;
      this.copy = args.copy;
      this.consumedStack = [];
      this._classname = 'libraryCopyWells';
      _.bindAll(this, 'showDetail');
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
      var self = this;
      var uriStack = this.uriStack;
      var resourceId = 'copywell';

      if (!_.isEmpty(uriStack) &&
              !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        // Detail view
        if (self.library && self.copy){
          // If in the context of libary/copy, then just the well id
          var well_id = uriStack.shift();
          this.consumedStack = [well_id];
          var _key = self.library.key + '/' + self.copy.get('copy_name')+ '/' + well_id;
          appModel.getModel(resourceId, _key, this.showDetail );
        } else {
          var copy_well_id = uriStack.shift();
          appModel.getModel(resourceId, copy_well_id, this.showDetail );
        }
      } else {
        if (_.isEmpty(uriStack)) {
          self.consumedStack = [];
          self.reportUriStack([]);
          self.showCopyWellSearchForm();
        } else {
          this.$("#tab_container-title").empty();
          this.consumedStack = [];
          this.showList();
        }
      }
    },
    
    showCopyWellSearchForm(){
      var self = this;
      var resource = appModel.getResource('copywell');
      ///// Copy Well search
      var copyWellsSchema = {};
      
      function validateCopyWellSearch(value,formValues){
        var errors = [];
        var parsedData = Iccbl.parseRawCopyWellSearch(value, errors);
        if (_.isEmpty(parsedData)){
          errors.push('no values found for input');
        } else {
          console.log('parsedData', parsedData);
        }
        if (!_.isEmpty(errors)){
          return {
            type: 'searchVal',
            message: errors.join('; ')
          };
        }
      };
      
      copyWellsSchema['searchVal'] = {
        title: 'Find Copy Wells',
        key: 'searchVal',
        placeholder: 'Copy Wells...',
        validators: ['required',validateCopyWellSearch],
        type: EditView.TextArea2,
        template: appModel._fieldTemplate,
        editorClass: 'form-control'
      };
      var CopyWellsFormModel = Backbone.Model.extend({
        schema: copyWellsSchema,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var copyWellsFormModel = new CopyWellsFormModel({
      });
      
      var copyWellsForm = new Backbone.Form({
        model: copyWellsFormModel,
        template: appModel._form_template,
      });

      var View = Backbone.Layout.extend({
        template: _.template(layout),
        afterRender: function(){
          $('#resource_content').html(copyWellsForm.render().el);
          
          copyWellsForm.$el.find('textarea').attr('rows',20);
          
          var helpLink = $('<a >&nbsp;?</a>');
          copyWellsForm.$el.find('.form-group').find('label').append(helpLink);
          helpLink.click(function(e){
            e.preventDefault();
            var bodyMessage = [
              '<b>By well ID and copy:</b>',
              '02091:B15 A',
              'A 2091:B15',
              '',
              '<b>Enter multiple wells, separated by newlines, or by commas:</b>',
              '02091:B15 C',
              '02091:B16 D',
              '02091:L19 E',
              'or',
              '02091:B15, 02091:B16, 02091:L19 C',
              '',
              '<b>By plate number, followed by well name(s):</b>',
              '2091 B15 C',
              '2091 B16 D',
              '2091 L19 E',
              'or',
              '2091 B15 B16 L19 C',
              '',
              '<b>Multiple wells and copies</b>',
              '2091 B15 B16 L19 C',
              '',
              '<b>Search for copies and plate range:</b>',
              '50001-50005 A',
              '50001-50005 A B C D',
              '50001,50002,50004 A B C D',
              '',
              '<b>Quoted copy name with a space in it:</b>',
              '50001-50005 "Stock A"',
            ];
            appModel.showModalMessage({
              title: 'Search for copy wells using plate, well, and copy patterns',
              body: bodyMessage.join('<br/>')
            });
          });

          copyWellsForm.$el.append([
            '<button type="submit" class="btn btn-default btn-xs" ',
            'style="width: 3em;">ok</input>',
          ].join(''));
          
    
          copyWellsForm.$el.find('[ type="submit" ]').click(function(e){
            e.preventDefault();
            copyWellsForm.$el.find('[data-error]').empty();
            var errors = copyWellsForm.commit({ validate: true }); 
            if(!_.isEmpty(errors)){
              copyWellsForm.$el.find('[name="searchVal"]').addClass(self.errorClass);
              return;
            }else{
              copyWellsForm.$el.find('[name="searchVal"]').removeClass(self.errorClass);
            }
            
            var text_to_search = copyWellsForm.getValue('searchVal');
            
            console.log('processCopyWellSearch: ', text_to_search);
            errors = [];
            var parsedSearchArray = Iccbl.parseRawCopyWellSearch(text_to_search,errors);
            if (!_.isEmpty(errors)){
              copyWellsForm.$el.find('[name="searchVal"]').addClass(self.errorClass);
              copyWellsForm.$el.find('[data-error]').html(errors.join('<br/>'));
              return;
            }
      
            
            if (parsedSearchArray.length <= appModel.MAX_RAW_SEARCHES_IN_URL){
              // encode simple searches as a URL param
              var encodedSearches = [];
              _.each(parsedSearchArray, function(parsedSearch){
                encodedSearches.push(_.map(
                  parsedSearch.combined, encodeURIComponent).join(','));
              });
              var uriStack = [resource.key, appModel.URI_PATH_ENCODED_SEARCH, 
                encodedSearches.join(appModel.UI_PARAM_RAW_SEARCH_LINE_ENCODE)];
              appModel.setUriStack(uriStack);
            }else{
              // must change the route, and create a post
              var searchId = ( new Date() ).getTime();
              appModel.setSearch(searchId,text_to_search);
              this.searchId = searchId;
              var uriStack = [
                resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
                searchId];
              appModel.setUriStack(uriStack);
            }
            
          }); // submit

        }
      });
      var view = new View();
      Backbone.Layout.setupView(view);
      self.setView('#tab_container', view).render();
      
    },
    
    showDetail: function(model) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      var view = new DetailLayout({ 
        model: model, 
        uriStack: uriStack,
        buttons: ['download','history']
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#resource_content', view).render();
      
      // FIXME: not working
      self.$("#tab_container-title").html('Well ' + model.get('well_id'));
    },

    showList: function() {
      console.log('show list...');
      var self = this;
      var uriStack = _.clone(this.uriStack);
      var extraControls = [];

      var resource = appModel.getResource('copywell');
      var url = [resource.apiUri];
      if (self.library && self.copy){
        url = [ 
          self.library.resource.apiUri,self.library.key,'copy',self.copy.get('copy_name'),
          'copywell'];
      }
      url = url.join('/');
      
      var Collection = Iccbl.MyCollection.extend({
        url: url
      });
      collection = new Collection();

      resource = appModel.cloneResource(resource);
      if(self.copy) {
        resource.fields['copy_usage_type']['visibility'] = ['d'];
        resource.fields['library_short_name']['visibility'] = ['d'];
      }
      
      // Set concentration type visibility based concentrations found in library wells
      resource.fields['mg_ml_concentration']['visibility'] = [];
      resource.fields['molar_concentration']['visibility'] = [];
      if (self.library){
        var concentration_types = self.library.get('concentration_types');
        if (concentration_types) {
          if (_.contains(concentration_types, 'mg_ml')){
            resource.fields['mg_ml_concentration']['visibility'] = ['l'];
          }
          if (_.contains(concentration_types, 'molar')){
            resource.fields['molar_concentration']['visibility'] = ['l'];
          }
        }
      }

      if (appModel.hasPermission(resource.key, 'write')){
        
        _.each(_.pick(resource['fields'], 
          ['volume', 'mg_ml_concentration', 'molar_concentration', 'cherry_pick_screening_count']), 
          function(field){
            field['editability'] = ['l','u'];
          });
        
        // Set up the grid to record edits of the "well_volume" column
        var showSaveButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'style="display: none; " ',
            'role="button" id="save_button_wells" href="#">',
            'save</a>'
          ].join(''));
        extraControls.push(showSaveButton);
        var showUploadButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="patch_resource" href="#">',
            'Upload data</a>'
          ].join(''));   
        extraControls.push(showUploadButton);
        
        showUploadButton.click(function(e){
          e.preventDefault();
          UploadDataForm.postUploadFileDialog(
            resource.apiUri, resource.content_types)
            .done(function(data, textStatus, jqXHR){
              appModel.showConnectionResult(data);
              collection.fetch({ reset: true });
            })
            .fail(function(){
              appModel.jqXHRfail.apply(this,arguments); 
            });
        });

        // Create a new collection of just the changed items 
        // (so that multi page changes can be collected)
        var PostCollection = Backbone.Collection.extend({
          url: url,
          toJSON: function(){
            return {
              objects: Collection.__super__.toJSON.apply(this) 
            };
          }
        });
        var changedCollection = new PostCollection();
        var CopyWellModel = Backbone.Model.extend({
          url: url,
          initialize : function() {
            this.on('change', function(model, options) {
              // Prevent save on update
              if (options.save === false)
                  return;
              if (!model.hasChanged()) {
                console.log('no changes');
                return;
              }
              var saveModel = new Backbone.Model(model.pick(resource.id_attribute));
              // Save the copy_usage_type to perform warning later if not CPR usage
              saveModel.set('copy_usage_type', model.get('copy_usage_type'));
              // Because the siUnit cell converts to float, must filter spurious changes
              var has_changes = false;
              
              // Check volume
              var modelField = 'volume';
              var newValue = parseFloat(model.get(modelField));
              var prevValue = parseFloat(model.previous(modelField));

              // Round the raw previous value
              var displayOptions = resource.fields[modelField]['display_options'];
              var defaultUnit = _.result(displayOptions,'defaultUnit');
              var precision = _.result(displayOptions, 'decimals');
              if (_.isNumber(defaultUnit) && _.isNumber(precision)){
                prevValue = Iccbl.roundForDefaultUnit(prevValue, precision, defaultUnit);
              }
              
              // NaN !== NaN
              if (!_.isNaN(newValue) && newValue !== prevValue){
                has_changes = true;
                console.log('field value change:', modelField, prevValue, newValue);
                saveModel.set(modelField, model.get(modelField));
              } else {
                self.$el.find('td.'+modelField).removeClass('edited');
                console.log('no change in field: ', modelField);
              }

              // Check the concentrations
              var modelField = 'mg_ml_concentration';
              var newValue = parseFloat(model.get(modelField));
              var prevValue = parseFloat(model.previous(modelField));

              // Round the raw previous value
              var displayOptions = resource.fields[modelField]['display_options'];
              var defaultUnit = _.result(displayOptions,'defaultUnit');
              var precision = _.result(displayOptions, 'decimals');
              if (_.isNumber(defaultUnit) && _.isNumber(precision)){
                prevValue = Iccbl.roundForDefaultUnit(prevValue, precision, defaultUnit);
              }
              
              // NaN !== NaN
              if (!_.isNaN(newValue) && newValue !== prevValue){
                has_changes = true;
                console.log('field value change:', modelField, prevValue, newValue);
                saveModel.set(modelField, model.get(modelField));
              } else {
                self.$el.find('td.'+modelField).removeClass('edited');
                console.log('no change in field: ', modelField);
              }

              var modelField = 'molar_concentration';
              var newValue = parseFloat(model.get(modelField));
              var prevValue = parseFloat(model.previous(modelField));

              // Round the raw previous value
              var displayOptions = resource.fields[modelField]['display_options'];
              var defaultUnit = _.result(displayOptions,'defaultUnit');
              var precision = _.result(displayOptions, 'decimals');
              if (_.isNumber(defaultUnit) && _.isNumber(precision)){
                prevValue = Iccbl.roundForDefaultUnit(prevValue, precision, defaultUnit);
              }
              
              // NaN !== NaN
              if (!_.isNaN(newValue) && newValue !== prevValue){
                has_changes = true;
                console.log('field value change:', modelField, prevValue, newValue);
                saveModel.set(modelField, model.get(modelField));
              } else {
                self.$el.find('td.'+modelField).removeClass('edited');
                console.log('no change in field: ', modelField);
              }
              
              var modelField = 'cherry_pick_screening_count';
              var newValue = model.get(modelField);
              var prevValue = model.previous(modelField);
              if (newValue != prevValue){
                has_changes = true;
                console.log('field value change:', modelField, prevValue, newValue);
                saveModel.set(modelField, model.get(modelField));
              } else {
                self.$el.find('td.'+modelField).removeClass('edited');
                console.log('no change in field: ', modelField);
              }
              
              if (!has_changes) {
                console.log('no changes');
                return;
              }

              changedCollection.add(saveModel);
              showSaveButton.show();
              appModel.setPagePending();
            });
          },
        });
        
        CopyWellModel._label = 'CopyWellModel';
        collection.model = CopyWellModel;
        
        collection.on('backgrid:error', function(model, column, value){
          console.log('backgrid error', arguments);
          // TODO: indicate errors in the cell
        });

        showSaveButton.click(function(e){
          e.preventDefault();
          console.log('changed collection', changedCollection,changedCollection.url);
          if(changedCollection.isEmpty()){
            appModel.showModalError('No changes to save');
            return;
          }
          
          var copy_usage_types = {};
          changedCollection.each(function(model){
            var type = model.get('copy_usage_type');
            console.log('checking', model.toJSON(), type);
            if (type !== 'cherry_pick_source_plates'){
              var typeCount = _.result(copy_usage_types,type, 0);
              copy_usage_types[type] = typeCount + 1;
            }
          });
          if (!_.isEmpty(copy_usage_types)){
            appModel.showModal({
              title: 'Copy types are not "cherry_pick_source_plates"',
              body: appModel.print_json(copy_usage_types),
              ok: function(){
                console.log('ok...');
              }
            });
          }
          var saveFunction = function() {
            appModel.showOkCommentForm({
              ok: function(formValues){
                console.log('form values', formValues);
                var comments = formValues['comments'];
                var headers = {};
                headers[appModel.HEADER_APILOG_COMMENT] = comments;
                
                Backbone.sync("patch",changedCollection,
                  {
                    headers: headers,
                    error: function(){
                      appModel.jqXHRfail.apply(this,arguments);
                      console.log('error, refetch', arguments);
                      changedCollection.reset();
                      collection.fetch({ reset: true });
                    },
                    success: function(){
                      changedCollection.reset();
                      collection.fetch({ reset: true });
                    }
                  }
                );
              }
            });
          
          };
          saveFunction();
          
        });
      
      } else { // User does not have write permission
        _.each(_.values(resource['fields']), function(field){
          field['editability'] = [];
        });
      }
      
      var showHistoryButton = $([
      '<a class="btn btn-default btn-sm pull-down" ',
        'role="button" id="showHistoryButton_wells" href="#">',
        'History</a>'
      ].join(''));
      
      
      if (self.library && self.copy){
        // Only show history if in the context of a library copy
        extraControls.push(showHistoryButton);
        showHistoryButton.click(function(e){
          e.preventDefault();
          var newUriStack = ['apilog','order','-date_time', appModel.URI_PATH_SEARCH];
          var search = {};
          search['ref_resource_name'] = 'copywell';
          search['key__icontains'] = [
            self.library.get('short_name'),
            self.copy.get('copy_name')].join('/');
          newUriStack.push(appModel.createSearchString(search));
          var route = newUriStack.join('/');
          console.log('history route: ' + route);
          appModel.router.navigate(route, {trigger: true});
          self.remove();
        });
      }
      

      resource.fields['copy_name'].backgridCellType = 
        Iccbl.LinkCell.extend(_.extend({},
          resource.fields['copy_name'].display_options,
          {
            linkCallback: function(e){
              e.preventDefault();
              self.trigger('showCopy');
              return;
            }
          }));
      resource.fields['plate_number'].backgridCellType = 
        Iccbl.LinkCell.extend(_.extend({},
          resource.fields['plate_number'].display_options,
          {
            linkCallback: function(e){
              e.preventDefault();
              console.log('link clicked: ', this.model.get('plate_number'));
              self.trigger('showPlate', this.model.get('plate_number'));
              return;
            }
          }));
            
      resource.fields['well_id'].backgridCellType = 
        Iccbl.LinkCell.extend(_.extend({},
          resource.fields['well_id'].display_options,
          {
            linkCallback: function(e){
              e.preventDefault();
              var well_id = this.model.get('well_id');
              self.consumedStack = [well_id];
              var _key = [this.model.get('copy_name'), well_id];
              if (self.library){
                _key.unshift(self.library.key);
              }
              appModel.getModel('copywell', _key.join('/'), self.showDetail );
            }
          }));
            
      
      var view = new ListView({ 
        uriStack: uriStack,
        resource: resource,
        url: url,
        collection: collection,
        extraControls: extraControls
      });
        
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#resource_content', view ).render();
    }
    
  });

  return LibraryCopyWellView;
});