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
  'text!templates/genericResource.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         LibraryView, layout) {
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout
  };
    
  var LibraryWellView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      this.uriStack = args.uriStack;
      this.library = args.library;
      this.consumedStack = [];
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

    // layoutmanager hook
    afterRender: function() {
      var self = this;
      
      var uriStack = this.uriStack;
      var library = this.library;

      // NOTE: library/well view is actually a "reagents" view
      var reagents_url = library.resource.apiUri +'/' + library.key + '/well';
      var url = reagents_url;

      appModel.getResourceFromUrl('well', reagents_url + '/schema', function(resource){

        var displayed_scopes = ['fields.well', 'fields.reagent'];
        if(self.library.get('screen_type') == 'rnai'){
          displayed_scopes.push('fields.silencingreagent');
        }else if(self.library.get('library_type') == 'natural_products'){
          displayed_scopes.push('fields.naturalproductreagent');
        }else{
          displayed_scopes.push('fields.smallmoleculereagent');
        }
        var specific_schema = _.clone(resource);
        specific_schema.fields = {}
        _.each(_.keys(resource.fields), 
            function(key){
              var field = resource.fields[key];
              if(_.contains(displayed_scopes, field['scope'])){
                specific_schema.fields[key] = field;
              }
            });
        
        if (!_.isEmpty(uriStack) && !_.isEmpty(uriStack[0]) &&
                !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
          var _key = Iccbl.popKeyFromStack(resource,uriStack,self.consumedStack);
  
          appModel.getModel('well', _key, function(model){
            model.resource = resource;
            self.showDetail(model, specific_schema);
          } );
        } else {
          self.consumedStack = [];
          self.showList(resource, reagents_url, specific_schema);
        }
      });        
    },    
    
    showDetail: function(model, specific_schema) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      // get the view class
      var viewClass = DetailLayout;

      // prevent normal display of structure image value
      var si_field = model.resource.fields['structure_image']
      if(si_field){
        si_field['visibility'] = [];
      }
      
      var view = new viewClass({ schema: specific_schema, model: model, uriStack: uriStack });

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.$('#content_title').html("");
      self.setView('#content', view).render();
      
      // TODO: support for generic images
      if(model.has('structure_image')){
        self.$('#content').append(
            '<img style="position: absolute; top: 8em; right: 3em; height: 16em;" src="' 
            + model.get('structure_image') + '" alt="image" />')
      }
    },
    
    showList: function(resource, url, specific_schema) {

      var self = this;
      var uriStack = _.clone(this.uriStack);
      
      // hack - should fit the reagent_type default columns into the metadata
      var library_type = this.library.get('library_type');
      console.log('library_type: ' + library_type);
      if( this.library.get('screen_type') == 'rnai'){
        uriStack.push('includes');
        uriStack.push('vendor_entrezgene_symbols,vendor_entrezgene_id');
      }else if(library_type != 'natural_products'){
        uriStack.push('includes');
        uriStack.push('compound_name');
      }
      
      
      // FIXME: deprecate this for the LinkCell...
      var detailHandler = function(model) {
        var key = Iccbl.getIdFromIdAttribute(model,resource);
        model.resource = resource;
        model.key = key;
        var keysToReport = Iccbl.getIdKeys(model,resource);
        self.consumedStack = keysToReport;
        self.showDetail(model);
      };
      var view = new ListView({ options: {
        uriStack: uriStack,
        schemaResult: specific_schema,
        resource: resource,
        url: url, 
        detailHandler: detailHandler
      }});
      self.listenTo(view, 'update_title', function(val){
        if(val) {
          this.$('#content_title').html('<small>' + val + '</small>');
        }else{
          this.$('#content_title').html("");
        }
      });
    
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#content', view ).render();
    }
  });

  return LibraryWellView;
});