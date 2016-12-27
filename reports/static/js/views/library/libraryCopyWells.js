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
  'utils/uploadDataForm',
  'templates/genericResource.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         LibraryView, UploadDataForm, layout) {
  
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
      
      var uriStack = this.uriStack;
      var library = this.library;
      var copy = this.copy;
      var url = [ 
          library.resource.apiUri,library.key,'copy',copy.get('copy_name'),
          'copywell'].join('/');
      var resourceId = 'copywell';
      var resource = appModel.getResource(resourceId);

      if (!_.isEmpty(uriStack) &&
              !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        // Detail view
        var well_id = uriStack.shift();
        this.consumedStack = [well_id];
        _key = library.key + '/' + copy.get('copy_name')+ '/' + well_id;
        appModel.getModel(resourceId, _key, this.showDetail );
      } else {
        // List view
        this.consumedStack = [];
        this.showList(resource, url);
        this.$("#tab_container-title").empty();
      }
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

    showList: function(resource, url) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      var extraControls = [];
      var Collection = Iccbl.MyCollection.extend({
        url: url
      });
      collection = new Collection({
        url: url,
      });

      resource = appModel.cloneResource(resource);
      if(self.copy) {
        resource.fields['copy_usage_type']['visibility'] = ['d'];
        resource.fields['library_short_name']['visibility'] = ['d'];
      }
      
      // Set concentration type visibility based concentrations found in library wells
      var concentration_types = self.library.get('concentration_types');
      resource.fields['mg_ml_concentration']['visibility'] = [];
      resource.fields['molar_concentration']['visibility'] = [];
      if (concentration_types) {
        if (_.contains(concentration_types, 'mg_ml')){
          resource.fields['mg_ml_concentration']['visibility'] = ['l'];
        }
        if (_.contains(concentration_types, 'molar')){
          resource.fields['molar_concentration']['visibility'] = ['l'];
        }
      }

      if (appModel.hasPermission(resource.key, 'write')){
        
        // Set field visibility based on copy type
        if(self.copy && self.copy.get('usage_type') == 'cherry_pick_source_plates'){
          _.each(_.pick(resource['fields'], 
            ['volume']), 
            function(field){
              field['editability'] = ['l','u'];
            });
        }
        _.each(_.pick(resource['fields'], 
          ['mg_ml_concentration','molar_concentration']), 
          function(field){
            field['editability'] = ['l','d'];
          });
        
        var showUploadButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="patch_resource" href="#">',
            'Upload data</a>'
          ].join(''));   
        extraControls.push(showUploadButton);
        
        showUploadButton.click(function(e){
          e.preventDefault();
          function callbackSuccess() {
          };
          
          UploadDataForm.postUploadFileDialog(
            collection.url, resource.content_types)
            .done(function(){
              collection.fetch({ reset: true });
            })
            .fail(function(){
              appModel.jqXHRfail.apply(this,arguments); 
            });
        });

        // Set up the grid to record edits of the "well_volume" column
        var showSaveButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'style="display: none; " ',
            'role="button" id="save_button_wells" href="#">',
            'save</a>'
          ].join(''));
        extraControls.push(showSaveButton);
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

              // Because the siUnit cell converts to float, must filter spurious changes
              var has_changes = false;
              var newValue = parseFloat(model.get("volume"));
              var prevValue = parseFloat(model.previous("volume"));
              // NaN !== NaN
              if (!_.isNaN(newValue) && newValue !== prevValue){
                has_changes = true;
                saveModel.set('volume', model.get('volume'));
              } else {
                console.log('no change in volume');
              }
              var newValue = parseFloat(model.get("mg_ml_concentration"));
              var prevValue = parseFloat(model.previous("mg_ml_concentration"));
              if (!_.isNaN(newValue) && newValue !== prevValue){
                has_changes = true;
                saveModel.set('mg_ml_concentration', model.get('mg_ml_concentration'));
              } else {
                console.log('no change in mg_ml_concentration');
              }
              var newValue = parseFloat(model.get("molar_concentration"));
              var prevValue = parseFloat(model.previous("molar_concentration"));
              if (!_.isNaN(newValue) && newValue !== prevValue){
                has_changes = true;
                saveModel.set('molar_concentration', model.get('molar_concentration'));
              } else {
                console.log('no change in molar_concentration');
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
            appModel.error('No changes to save');
            return;
          }
          
          var saveFunction = function() {
            appModel.showSaveWithComments(function(formValues){
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
      extraControls.push(showHistoryButton);
      showHistoryButton.click(function(e){
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', 'search'];
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
              _key = [self.library.key, self.copy.get('copy_name'), well_id].join('/');
              appModel.getModel('copywell', _key, self.showDetail );
            }
          }));
            
      
      var view = new ListView({ 
        _name: 'WellsListView',
        uriStack: uriStack,
        schemaResult: resource,
        resource: resource,
        url: url,
        collection: collection,
        extraControls: extraControls
      });
      view._instanceName = 'WellsListView_instance';
        
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#resource_content', view ).render();
    }
    
  });

  return LibraryCopyWellView;
});