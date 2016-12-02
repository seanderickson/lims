define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_layout',
    'views/generic_edit',
    'views/list2',
    'templates/generic-tabbed.html',
    'test/models/detail_test_resource.json', 
    'test/models/detail_test_model.json',
    'test/models/detail_test_vocabularies.json',    
    'test/models/detail_test_users.json'    
    
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout,EditView, 
            ListView, layout, 
            resource_raw, test_model_raw, test_vocabularies_raw, test_users_raw ) {

  var DetailView = Backbone.Layout.extend({

    initialize: function(args) {
      console.log('1detail test initialize');
      var self = this;
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      try{
        this.resource = JSON.parse(resource_raw);
        _.each(_.values(this.resource.schema.fields), appModel.parseSchemaField );
      }catch(e){
        console.log('error reading resource,',e);
        throw e;
      }
      if(_.isUndefined(this.model)){
        try{
          
          this.model = new Backbone.Model(JSON.parse(test_model_raw));
        }catch(e){
          console.log('error reading model,',e);
          throw e;
        }
      }
      this.model.url = 'testingUrl';
      this.model.resource = this.resource;
      this.model.resource.schema = _.extend(this.model.resource.schema, new Iccbl.SchemaClass());
      this.title = Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema);
      try{
        _.extend(appModel.get("vocabularies"), JSON.parse(test_vocabularies_raw));
      }catch(e){
        console.log('error reading vocabularies,',e);
        throw e;
      }
      try{
        if(appModel.has('users')){
          appModel.set("users", new Backbone.Collection(
            _.union(appModel.get("users"), JSON.parse(test_users_raw))));
        }else{
          appModel.set("users", new Backbone.Collection(
            JSON.parse(test_users_raw)));
        }
        appModel.getUserOptions(function(options){
          self.model.resource.schema.fields['field13']['choices'] = options;
        });
      }catch(e){
        console.log('error reading users,',e);
        throw e;
      }
      
      _.bindAll(this, 'click_tab');
      console.log('detail test initialize done');
    },
    
    template: _.template(layout),
    
    tabbed_resources: {
        detail: { 
          description: 'Object Details', 
          title: 'Object Details', 
          invoke: 'setDetail' },
        usergroup: { 
          description: 'Object List', 
          title: 'Object List', 
          invoke: 'setList', 
          resource: 'test_resource' }
    },
    
    events: {
      'click ul.nav-tabs >li': 'click_tab',
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
      console.log('afterRender...', this.uriStack);
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if(viewId == '+add'){
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
      console.log('afterRender: viewId: ', viewId);
      this.change_to_tab(viewId);
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    click_tab : function(event){
      var self = this;
      event.preventDefault();
      // Block clicks from the wrong elements
      // TODO: how to make this specific to this view? (it is also catching
      //clicks on the table paginator)
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      
      if(this.key && this.key === key){
        return;
      }
      
      appModel.requestPageChange({
        ok: function(){
          self.change_to_tab(key);
        }
      });
    },

    change_to_tab: function(key){
      console.log('change to tab: ', key, this.uriStack)
      if(_.has(this.tabbed_resources, key)){
        var delegateStack = _.clone(this.uriStack);
        if(!_.isUndefined(this.tabbed_resources[key].alias)){
          key = this.tabbed_resources[key].alias;
        }
        if(this.key && this.key === key){
          return;
        }else{
          this.key = key;
        }        
        
        this.$('li').removeClass('active'); // TODO: use bootstrap tabs
        this.$('#' + key).addClass('active');
        
        this.uriStack = [];
        var method = this[this.tabbed_resources[key]['invoke']];
        if (_.isFunction(method)) {
          method.apply(this,[delegateStack]);
        } else {
          throw "Tabbed resources refers to a non-function: " + this.tabbed_resources[key]['invoke']
        }
      }else{
        var msg = 'Unknown tab: ' + key;
        window.alert(msg);
        throw msg;
      }
    },
    
    setDetail: function(delegateStack){
      console.log('setDetail', delegateStack);
      // TODO: create/or/import a model & schema that exercises test cases
      // TODO: write mocha test cases
      
      var key = 'detail';
 
      var view = this.tabViews[key];
      var saveCallBack = function(model){
        alert('save model: ' + JSON.stringify(model.toJSON()));
      };
      
      if ( !view ) {
        view = new DetailLayout({ 
          model: this.model, 
          uriStack: delegateStack,
          saveCallBack: saveCallBack
        });
        this.tabViews[key] = view;
      }
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // Note: since detail_layout reports the tab, the consumedStack is empty here
      this.consumedStack = []; 
      this.setView("#tab_container", view ).render();
      return view;
    },

    setList: function(delegateStack){

      // TODO: FOR TESTING: create or import a list of test objects
      
      
      var self = this;
      var key = 'test_resource';
      var resource = appModel.getResource('usergroup');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'groups'].join('/');
      view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource.schema,
        resource: resource,
        url: url
      });
      Backbone.Layout.setupView(view);
      this.consumedStack = [key]; 
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
    },    
    
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });

  return DetailView;
});