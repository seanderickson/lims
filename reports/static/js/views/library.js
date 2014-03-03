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
  'views/generic_detail_layout',
  'views/list2',
  'text!templates/library.html'
], function($, _, Backbone, Iccbl, layoutmanager, appModel, DetailLayout, ListView, libraryTemplate) {
  
  var LibraryView = Backbone.Layout.extend({
    
    initialize: function(args) {
      this.views = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      console.log('uriStack: ' + JSON.stringify(this.uriStack));
      _.bindAll(this, 'click_tab');
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

    reportUriStack: function(reportedUriStack) {
      var actualStack = this.consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    
    change_to_tab: function(key){
      this.$('li').removeClass('active');
      this.$('#' + key).addClass('active');

      if( _.has(this.tabbed_resources, key)){
        _.result(this, this.tabbed_resources[key]['invoke']);
        this.consumedStack = [key];
        this.trigger('uriStack:change', [key] );
      }else{
        window.alert('Unknown tab: ' + key);
      }
    },
    
    setDetail: function() {
      // for testing: 
      var key = 'detail';
      
      var view = this.views[key];
      if ( !view ) {
        view = new DetailLayout({ model: this.model});
        this.views[key] = view;
      }
      this.setView("#tab_container", view ).render();
    },
    
    setCopies: function() {
      var self = this;
      var key = 'copy';
      this.consumedStack.push(key);
      
      var view = this.views[key];
      if ( !view ) {
        var libraryKey = self.model.key;
        var _url = self.model.resource.apiUri +'/' + libraryKey + '/copy/'
        var apiResourceId = 'librarycopy';
        var libraryResource = appModel.getResource(apiResourceId);
        
        appModel.getSchema(apiResourceId, function(schema) {
          
          var view = new ListView({ 
            model: appModel, 
            options: {
              resource: libraryResource,
              url: _url , 
              schemaResult: schema,
              uriStack: _.clone(self.uriStack) 
            } 
          });
          self.listenTo(view , 'uriStack:change', self.reportUriStack)
          Backbone.Layout.setupView(view);
          
          self.views[key] = view;
          Backbone.Layout.setupView(view);
          
          self.setView("#tab_container", view ).render();
        });
      } else {
        this.setView("#tab_container", view ).render();
      }
    },
    
    setContents: function() {
      // for testing: 
      var key = 'contents';
      var view = this.views[key];
      if ( !view ) {
        view = new DetailLayout({ model: this.model});
        this.views[key] = view;
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
        this.consumedStack = [viewId];
        
        if (!_.has(this.tabbed_resources, viewId)){
          window.alert('could not find the tabbed resource: ' + viewId);
          return;
        }
      }
      this.change_to_tab(viewId);
    },
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.views = {};
      this.remove();
    },
  
    template: _.template(libraryTemplate)
  });
  
  return LibraryView;
});