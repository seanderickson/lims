/**
 * Library form/view
 */
define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/library/libraryCopies', 
  'views/library/libraryCopyPlates', 
  'views/library/libraryWells', 
  'views/library/libraryVersions',
  'views/generic_detail_layout',
  'views/generic_edit',
  'views/list2',
  'text!templates/library.html'
], function($, _, Backbone, Iccbl, layoutmanager, appModel, LibraryCopiesView, 
            LibraryCopyPlatesView, LibraryWellsView, LibraryVersionsView,
            DetailLayout, EditView, ListView, 
            libraryTemplate) {
  
  var LibraryView = Backbone.Layout.extend({
    
    initialize: function(args) {
      var self = this;
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      
      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail' && !appModel.hasPermission(self.tabbed_resources[key].resource)){
          delete self.tabbed_resources[key];
        }
      });
      console.log('check1');
      _.bindAll(this, 'click_tab');
      _.bindAll(this, 'saveFile');
      
  
    },
    
    template: _.template(libraryTemplate),

    tabbed_resources: {
      detail: { 
        description: 'Library Details', 
        title: 'Library Details', 
        invoke: 'setDetail'
      },
      copy: { 
        description: 'Copies', title: 'Copies', invoke: 'setCopies',
        resource: 'librarycopy'
      },
      plate: { 
        description: 'Plates', title: 'Plates', invoke: 'setPlates' ,
        resource: 'librarycopyplate'
      },
      well: { 
        description: 'Well based view of the library contents', 
        title: 'Wells', invoke: 'setWells' ,
        resource: 'well'
      },
      version: { 
        description: 'Library contents version', 
        title: 'Versions', invoke: 'setVersions' ,
        resource: 'librarycontentsversion'
      }
    },      
    
    events: {
      // TODO: how to make this specific to this view? 
      // (it is also catching clicks on the table paginator)
      //          'click .tabbable li': 'click_tab', 
        'click li': 'click_tab',
        'click button#upload': 'upload'        
    },
    
    upload: function(event){
      appModel.showModal({
          ok: this.saveFile,
          body: '<input type="file" name="fileInput" />',
          title: 'Select a SDF file to upload'  });
    },
    
    saveFile: function() {
      var file = $('input[name="fileInput"]')[0].files[0]; 
      var data = new FormData();
      var url = _.result(this.model, 'url') + '/well';
      // NOTE: 'sdf' is sent as the key to the FILE object in the upload.
      // We are using this as a non-standard way to signal the upload type to the 
      // serializer. (TP doesn't support mulitpart uploads, so this is patched in).
      data.append('sdf', file);
//      data.append('file', file, {'Content-type': 'chemical/x-mdl-sdfile'});
//      data.append('Content-type','chemical/x-mdl-sdfile');
      $.ajax({
        url: url,    
        data: data,
        cache: false,
        contentType: false,
        processData: false,
        type: 'PUT',
        headers: {
          'APILOG_COMMENT': 'TODO: comment for library content load'
        },        
        success: function(data){
          // FIXME: should be showing a regular message
          appModel.error('success');
          // FIXME: should reload the library view to reflect new comments/etc.
        },
        done: function(model, resp){
          // TODO: done replaces success as of jq 1.8
          console.log('done');
        },
        error: appModel.jqXHRFail
//        error: function(data){
//          alert('no upload');
//          $('#loadingModal').modal('hide');
//        }
      });
    },
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      return {
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if (viewId == '+add') {
          this.uriStack.unshift(viewId);
          this.showAdd();
          return;
        }else if (viewId == 'edit'){
          this.uriStack.unshift(viewId);
          this.showEdit();
          return;
        }

        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },
    
    click_tab : function(event){
      event.preventDefault();
      // Block clicks from the wrong elements
      // TODO: how to make this specific to this view? (is it still also catching
      // clicks on the table paginator?)
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      this.change_to_tab(key);
    },

    showAdd: function() {
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack
      });
      Backbone.Layout.setupView(view);

      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },
    
    showEdit: function() {
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack, 
        buttons: ['download', 'upload']
      });
      Backbone.Layout.setupView(view);

      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },
    
    change_to_tab: function(key){
      if(_.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active');
        this.$('#' + key).addClass('active');
        if(key !== 'detail'){
          this.consumedStack = [key];
        }else{
          this.consumedStack = [];
        }
        var delegateStack = _.clone(this.uriStack);
        this.uriStack = [];
        var method = this[this.tabbed_resources[key]['invoke']];
        if (_.isFunction(method)) {
          method.apply(this,[delegateStack]);
        } else {
          throw "Tabbed resources refers to a non-function: " + this.tabbed_resources[key]['invoke']
        }
      }else{
        var msg = 'Unknown tab: ' + key;
        appModel.error(msg);
        throw msg;
      }
    },
    
    setDetail: function(delegateStack) {
      var key = 'detail';
      
      var view = this.tabViews[key];
      if ( !view ) {
        view = new DetailLayout({ 
          model: this.model,
          uriStack: delegateStack, 
          buttons: ['download', 'upload'] });
        this.tabViews[key] = view;
      } 
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // NOTE: if subview doesn't report stack, report it here
      //      this.reportUriStack([]);
      this.setView("#tab_container", view ).render();
    },

    setVersions: function(delegateStack) {
      var self = this;
      var key = 'version';
      var view = this.tabViews[key];
      if ( !view ) {
        var view = new LibraryVersionsView({
          library: self.model,
          uriStack: delegateStack
        });
        self.tabViews[key] = view;
        Backbone.Layout.setupView(view);
      } else {
        self.reportUriStack([]);
      }
      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    setWells: function(delegateStack) {
      var self = this;
      var key = 'well';
      var view = this.tabViews[key];
      if ( !view ) {
        var view = new LibraryWellsView({
          library: self.model,
          uriStack: delegateStack
        });
        self.tabViews[key] = view;
        Backbone.Layout.setupView(view);
      } else {
        self.reportUriStack([]);
      }
      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    setPlates: function(delegateStack) {
      var self = this;
      var key = 'plate';
      var view = this.tabViews[key];
      if ( !view ) {
        var view = new LibraryCopyPlatesView({
          library: self.model,
          uriStack: delegateStack
        });
        self.tabViews[key] = view;
        Backbone.Layout.setupView(view);
      } else {
        self.reportUriStack([]);
      }
      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    setCopies: function(delegateStack) {
      var self = this;
      var key = 'copy';
      var view = this.tabViews[key];
      if ( !view ) {
        var view = new LibraryCopiesView({
          library: self.model,
          uriStack: delegateStack
        });
        self.tabViews[key] = view;
        Backbone.Layout.setupView(view);
      } else {
        self.reportUriStack([]);
      }
      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }
  
  });
  
  return LibraryView;
});