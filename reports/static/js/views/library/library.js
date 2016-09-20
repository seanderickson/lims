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
  'views/library/libraryCopy', 
  'views/library/libraryWell', 
  'views/generic_detail_layout',
  'views/generic_edit',
  'views/list2',
  'templates/generic-tabbed.html'
], function($, _, Backbone, Iccbl, layoutmanager, appModel, LibraryCopyView, 
            LibraryWellView, DetailLayout, EditView, ListView, 
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
      }
    },      
    
    events: {
        'click ul.nav-tabs >li': 'click_tab',
        'click button#upload': 'upload'        
    },
    
    upload: function(event){
      appModel.showModal({
          ok: this.saveFile,
          body: '<input type="file" name="fileInput" />',
          title: 'Select a SDF file to upload'  });
    },
    
    saveFile: function() {
      var self=this;
      var file = $('input[name="fileInput"]')[0].files[0]; 
      var data = new FormData();
      var url = _.result(this.model, 'url') + '/well';
      // NOTE: 'sdf' is sent as the key to the FILE object in the upload.
      // We are using this as a non-standard way to signal the upload type to the 
      // serializer. (TP doesn't support mulitpart uploads, so this is patched in).
      data.append('sdf', file);
      var headers = {};
      headers[appModel.HEADER_APILOG_COMMENT] = self.model.get('comment');
      $.ajax({
        url: url,    
        data: data,
        cache: false,
        contentType: false,
        processData: false,
        type: 'POST',
        headers: headers
      })
      .fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); })
      .done(function(model, resp){
          // FIXME: should be showing a regular message
          appModel.error('success');
          self.model.fetch();
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
        'base_url': this.model.resource.key + '/' + this.model.key,
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
        if(viewId == '+add'){
          this.$('ul.nav-tabs > li').addClass('disabled');
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }
        if(viewId == 'edit'){
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          window.alert(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },    
    
    click_tab : function(event){
      var self = this;
      event.preventDefault();
      event.stopPropagation();
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      appModel.requestPageChange({
        ok: function(){
          self.change_to_tab(key);
        }
      });
    },

    change_to_tab: function(key){
      if(_.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active');
        this.$('#' + key).addClass('active');
        this.$("#tab_container-title").empty();
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
    
    showAdd: function() {
      console.log('add view');
      
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack
      });
      Backbone.Layout.setupView(view);

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

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },
    
    setDetail: function(delegateStack) {
      console.log('detail view');
      
      var self = this;
      var key = 'detail';
      var buttons = ['download'];
      if (appModel.hasPermission('library', 'write')){
        buttons = buttons.concat(['upload','history','edit']);
      }
      
      // Custom library model validation: plate range
      this.model.validate = function(attrs) {
        var errs = {};
        console.log('validating: ', attrs); 
        if (!_.isNumber(attrs.start_plate)){
          errs.start_plate = 'number required';
        }else{
          var ranges = _.result(self.model.resource,'library_plate_ranges',[]);
          var rangesAsString = _.map(_.zip(
            _.filter(ranges, function(val, index){ return index%2==0; }),
            _.filter(ranges, function(val, index){ return index%2==1; })),
            function(pair){ return pair.join('-'); } ).join(', ');
          for (var i=0; i<ranges.length-1; i++){
            var start = ranges[i];
            var end = ranges[i+1];
            if (attrs.start_plate >= start && attrs.start_plate <= end){
              errs.start_plate = 'range already used: ' + rangesAsString;
              break;
            }
            if (attrs.end_plate >= start && attrs.end_plate <= end){
              errs.start_plate = 'range already used: ' + rangesAsString;
              break;
            }
            i++;
          }
        }
        if (!_.isEmpty(errs)) return errs;
      };

      var view = this.tabViews[key];
      if (view) {
        // remove the view to refresh the page form
        this.removeView(this.tabViews[key]);
      }
      
      view = new DetailLayout({ 
        model: this.model,
        uriStack: delegateStack, 
        buttons: buttons });
      this.tabViews[key] = view;

      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    setWells: function(delegateStack) {
      var self = this;
      var schemaUrl = [self.model.resource.apiUri,self.model.key,'well',
                       'schema'].join('/');
      var wellResource = appModel.getResource('well'); 

      function createWellView(schemaResult){
        if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
            !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
          // Detail view
          var well_id = delegateStack.shift();
          self.consumedStack.push(well_id);
          var _key = [self.model.key,well_id].join('/');
          appModel.getModel(wellResource.key, well_id, function(model){
            model.resource = schemaResult;
            view = new LibraryWellView({
              model: model, 
              uriStack: _.clone(delegateStack),
              library: self.model
            });
            Backbone.Layout.setupView(view);
            self.listenTo(view , 'uriStack:change', self.reportUriStack);
            self.setView("#tab_container", view ).render();
            self.$("#tab_container-title").html('Well ' + model.get('well_id'));
      
          });        
          return;
        }else{
          // List view
          var url = [self.model.resource.apiUri, 
                     self.model.key,
                     'well'].join('/');
          view = new ListView({ 
            uriStack: _.clone(delegateStack),
            schemaResult: schemaResult,
            resource: schemaResult,
            url: url
          });
          Backbone.Layout.setupView(view);
          self.reportUriStack([]);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
        }
        
      }
      appModel.getResourceFromUrl(schemaUrl, createWellView);

    },
        
    setPlates: function(delegateStack) {
      var self = this;
      var copyPlateResource = appModel.getResource('librarycopyplate'); 

      // List view
      // special url because list is a child view for library
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'plate'].join('/');
      view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: copyPlateResource,
        resource: copyPlateResource,
        url: url
      });
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },

    setCopies: function(delegateStack) {
      var self = this;
      var copyResource = appModel.getResource('librarycopy'); 
      if(!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0]) ){
        // Detail view
        var copyname = delegateStack.shift();
        self.consumedStack.push(copyname);
        var _key = [self.model.key,copyname].join('/');
        appModel.getModel(copyResource.key, _key, function(model){
          view = new LibraryCopyView({
            model: model, 
            uriStack: _.clone(delegateStack),
            library: self.model
          });
          Backbone.Layout.setupView(view);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
          
          console.log('title: ', Iccbl.getTitleFromTitleAttribute(model, model.resource));
          self.$("#tab_container-title").html('Copy ' + model.get('name'));
        });        
        return;
      }else{
        // List view
        var url = [self.model.resource.apiUri, 
                   self.model.key,
                   'copy'].join('/');
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: copyResource,
          resource: copyResource,
          url: url
        });
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.$("#tab_container-title").empty();
        self.setView("#tab_container", view ).render();
      }
    },

    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }
  
  });
  
  return LibraryView;
});