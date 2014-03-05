/**
 * Screen form/view
 *
 */
define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/libraryCopy', 
  'views/generic_detail_layout',
  'views/list2',
  'text!templates/library.html'
], function($, _, Backbone, Iccbl, layoutmanager, appModel, LibraryCopyView, 
            DetailLayout, ListView, libraryTemplate) {
  
  var LibraryView = Backbone.Layout.extend({
    
    initialize: function(args) {
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      _.bindAll(this, 'click_tab');
      _.bindAll(this, 'setDetail');
    },
    
    tabbed_resources: {
      detail: { description: 'Library Details', title: 'Library Details', invoke: 'setDetail' },
      copy: { description: 'Copies', title: 'Copies', invoke: 'setCopies' },
      contents: { description: 'Library Contents', title: 'Library Contents', invoke: 'setContents' },
    },      
    
    events: {
      // TODO: how to make this specific to this view? 
      // (it is also catching clicks on the table paginator)
      //          'click .tabbable li': 'click_tab', 
        'click li': 'click_tab',
    },
    
    click_tab : function(event){
      event.preventDefault();

      // Block clicks from the wrong elements
      // TODO: how to make this specific to this view? (it is also catching
      //clicks on the table paginator)
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;

      this.change_to_tab(key);
    },

    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    
    change_to_tab: function(key){
      if( _.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active');
        this.$('#' + key).addClass('active');
        this.consumedStack = [key];
        _.result(this, this.tabbed_resources[key]['invoke']);
        //this.trigger('uriStack:change', [key] );
      }else{
        var msg = 'Unknown tab: ' + key;
        window.alert(msg);
        throw msg;
      }
    },
    
    setDetail: function() {
      var key = 'detail';
      
      var view = this.tabViews[key];
      if ( !view ) {
        view = new DetailLayout({ model: this.model});
        this.tabViews[key] = view;
      }
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // NOTE: since subview doesn't report stack, report it here
//      this.reportUriStack([]);
      this.setView("#tab_container", view ).render();
    },
    
    setCopies: function() {
      var self = this;
      var key = 'copy';
      var view = this.tabViews[key];
      if ( !view ) {
        var view = new LibraryCopyView({
          library: self.model,
          uriStack: _.clone(self.uriStack)
        });
        self.tabViews[key] = view;
        Backbone.Layout.setupView(view);
//        var apiResourceId = 'librarycopy';
//        appModel.getSchema(apiResourceId, function(schema) {
//          
//          var lcUrl = this.model.resource.apiUri +'/' + this.model.key + '/copy/'
//          var lcResource = appModel.getResource(libraryCopyResourceId);
//          var view = new ListView({ 
//            options: {
//              resource: lcResource,
//              url: lcUrl , 
//              schemaResult: schema,
//              uriStack: _.clone(self.uriStack) 
//            } 
//          });
//          self.listenTo(view , 'uriStack:change', self.reportUriStack);
//          self.listenTo(view, 'detail', function(model){
//
//            var key = Iccbl.getIdFromIdAttribute(model,schema);
//            var keysToReport = Iccbl.getIdKeys(model,schema);
//            if(keysToReport[0] = this.model.key){
//              keysToReport.shift(); // get rid of the library key part
//            }
//            
//            appModel.updateModel(apiResourceId,key,model,function(model){
//              // TODO: this will be a "LibraryCopy" view
//              var newView = new DetailLayout({model: model, uriStack:[] });
//              self.reportUriStack(keysToReport);
//              self.listenTo(newView , 'uriStack:change', function(uriStack) {
//                // This is proof, despite previous failures at LibraryCopyView
//                // of the need to nest this section in a view
//                self.reportUriStack(keysToReport.concat(uriStack));
//              });
//
//              self.setView("#tab_container", newView).render();
//            });
//          });
          


//          var libraryKey = self.model.key;
//          var _url = self.model.resource.apiUri +'/' + libraryKey + '/copy/'
//          var libraryCopyResource = appModel.getResource(apiResourceId);
//          
//          
//          var copyModel = new Backbone.Model({
//            resource: libraryCopyResource,
//            listUri: _url,
//            schemaResult: schema,
//            library: self.model,
//            uriStack: _.clone(self.uriStack)
//          });
//          var view = new LibraryCopyView({ 
//            model: copyModel
//          });
//          
//          copyModel.on('change:uriStack', function(model, val, options){
//            self.reportUriStack(model.get('uriStack'));
//          });
//        });
      }
      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    setContents: function() {
      // for testing: 
      var key = 'contents';
      var view = this.tabViews[key];
      if ( !view ) {
        view = new DetailLayout({ model: this.model});
        this.tabViews[key] = view;
      }
      // testing only
      this.model.set({'library_name': this.model.get('library_name') + ", tab: " + key }, {silent: true});
      this.setView("#tab_container", new DetailLayout({ model: this.model})).render();
    },
    
    // layoutmanager hook
    serialize: function() {
      return {
        'title': Iccbl.getTitleFromTitleAttribute(this.model, this.model.resourceSchema),
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    // layoutmanager hook
    afterRender: function(){
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
//        this.consumedStack = [viewId];
        
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          window.alert(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    },
  
    template: _.template(libraryTemplate)
  });
  
  return LibraryView;
});