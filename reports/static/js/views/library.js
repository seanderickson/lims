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
  'text!templates/library.html'
], function($, _, Backbone, Iccbl, layoutmanager, appModel, DetailLayout, libraryTemplate) {
  
  var LibraryView = Backbone.Layout.extend({
    
    initialize: function() {
      this.views = {}; // view cache
    },
    
    tabbed_resources: {
      detail: { description: 'Library Details', title: 'Library Details', invoke: 'setDetail' },
      copies: { description: 'Copies', title: 'Copies', invoke: 'setCopies' },
      plates: { description: 'Plates', title: 'Plates', invoke: 'setPlates' },
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

      // signal to the app_model that the current view has changed 
      var current_options = {'tab': key, 'key':Iccbl.getKey(appModel.get('current_options') ) };
      appModel.set({ current_options: {}, current_detail: {} },{ silent: true });
      appModel.set({
          current_options: current_options,
          routing_options: {trigger: false, replace: false}
      }); 
    },
    
    change_to_tab: function(key){
      this.$('li').removeClass('active');
      this.$('#' + key).addClass('active');

      if( _.has(this.tabbed_resources, key)){
        _.result(this, this.tabbed_resources[key]['invoke']);
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
    
    setPlates: function() {
      // for testing: 
      var key = 'plates';
      var view = this.views[key];
      if ( !view ) {
        view = new DetailLayout({ model: this.model});
        this.views[key] = view;
      }
      // testing only
      this.model.set({'library_name': this.model.get('library_name') + ", tab: " + key }, {silent: true});
      
      this.setView("#tab_container", view ).render();
      
    },    
    
    setCopies: function() {
      // for testing: 
      var key = 'copies';
      var view = this.views[key];
      if ( !view ) {
        view = new DetailLayout({ model: this.model});
        this.views[key] = view;
      }
      // testing only
      this.model.set({'library_name': this.model.get('library_name') + ", tab: " + key }, {silent: true});
      
      this.setView("#tab_container", view ).render();
      
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
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    // layoutmanager hook
    afterRender: function(){
      this.change_to_tab('detail');
    },
    
    onClose: function() {
      console.log('-------------- onclose ---------------------');
      this.views = {};
      this.remove();
    },
  
    template: _.template(libraryTemplate)
  });
  
  return LibraryView;
});