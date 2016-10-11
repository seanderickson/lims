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
          library.resource.apiUri,library.key,'copy',copy.get('name'),
          'copywell'].join('/');
      var resourceId = 'copywell';
      var resource = appModel.getResource(resourceId);

      if (!_.isEmpty(uriStack) && !_.isEmpty(uriStack[0]) &&
              !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        // Detail view
        var well_id = uriStack.shift();
        this.consumedStack = [well_id];
        _key = library.key + '/' + copy.get('name')+ '/' + well_id;
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
      
      if (appModel.hasPermission(resource.key, 'write')){
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
              var newValue = parseFloat(model.get("volume"));
              var prevValue = parseFloat(model.previous("volume"));
              if (newValue == prevValue){
                console.log('no change');
                return;
              }
              model.url = url;
              changedCollection.add(model);
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
        });
      
      }else {
        // Turn off editability on the "volume" field
        var volume_field = _.result(resource['fields'], 'volume', null);
        if (volume_field) {
          volume_field['editability'] = [];
        }
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
          self.copy.get('name')].join('/');
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        console.log('history route: ' + route);
        appModel.router.navigate(route, {trigger: true});
        self.remove();
      });

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