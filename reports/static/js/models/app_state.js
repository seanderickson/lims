define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid',
  'text!templates/modal_ok_cancel.html'    
], function($, _, Backbone, Iccbl, modalOkCancelTemplate ){
  
  var API_VERSION = 'api/v1';
  var REPORTS_API_URI = '/reports/' + API_VERSION;
  var DB_API_URI = '/db/' + API_VERSION;
  var DEBUG = true;
  

  var SchemaClass = function() {};
  SchemaClass.prototype.detailKeys = function()
  {
    var self = this;
    
    var keys = Iccbl.sortOnOrdinal(
      _.keys(self.fields), self.fields)
    var detailKeys = _(keys).filter(function(key){
      return _.has(self.fields, key) && 
          _.has(self.fields[key], 'visibility') && 
          _.contains(self.fields[key]['visibility'], 'detail');
    });
    return detailKeys;

  };  

  var AppState = Backbone.Model.extend({
    
    defaults: {

      // TODO: deprecate these variables
      // use the REPORTS_API_URI, and DB_API_URI defined below
      root_url: '/reports',  // used for the backbone history
      api_root_url: '/reports/api/v1',

      path: '',
      actualStack: [],
      
      // NOTE: in backbone, change notifications are only detected for the 
      // object identity.  For this reason, each of the following state items is
      // kept at the top level.
      current_view: 'home',
      current_resource_id: 'home',
      current_options: {},
      routing_options: {},
      current_scratch: {},
      current_details: {},
      login_information: {},


      list_defaults: {
          page: 1,
          rpp: 25,
          order: null,
          search: null,
      },
      detail_defaults: {
      },

    },

    initialize : function() {
      _.bindAll(this,'jqXHRerror', 'error');
    },
        
    start: function(callBack) {
      var self = this;
      console.log('start app_state.js');
      this.getResources(function(){
        self.setCurrentUser(callBack);
      });
    },
    
    setCurrentUser: function(callBack) {
      var self = this;
      this.getModel('user', window.user, function(model){
        self.currentUser = model.toJSON();
        console.log('Current user: ' + JSON.stringify(this.currentUser));
        
        if(callBack) callBack(); 
      });
    },
    
    getCurrentUser: function(){
      return this.currentUser;
    },

    getResources: function(callBack){
      console.log('- getResources from the server...')
      var self = this;
      // Retrieve the resource definitions from the server
      var resourceUrl = self.reportsApiUri + '/resource'
      Iccbl.getCollectionOnClient(resourceUrl, function(collection){
        
        // Store the URI for each resource.
        // TODO: store the apiUri on the resource in the server
        var schemaClass = new SchemaClass();
        var resourcesCollection = collection.toJSON();
        _.each(resourcesCollection,function(resource){
          resource.apiUri = '/' + resource.api_name + '/' + 
            self.apiVersion + '/' + resource.key;
          resource.schema = _.extend(resource.schema, schemaClass);
        });

        // 1. Create a resource hash, keyed by the resource id key
        // 2. Augment uiResources with the (api)resources
        var ui_resources = self.get('ui_resources');
        var resources = {};
        _.each(resourcesCollection, function(resource){
          resources[resource.key] = resource;
          if (_.has(ui_resources, resource.key)) {
            ui_resources[resource.key] = _.extend(
                {}, ui_resources[resource.key], resource);
          } else {
            ui_resources[resource.key] = _.extend({},resource);
          }
        });

        // 3. for "virtual" uiResources, find the underlying apiResource, 
        // and augment the uiResource with it
        _.each(_.keys(ui_resources), function(key){
          var ui_resource = ui_resources[key];
          // then augment the ui resources with their api_resource, if different
          if ( key !== ui_resource.api_resource ) {
            ui_resources[key] = _.extend(
                {}, resources[ui_resource.api_resource], ui_resource);
          }
        });
        if(callBack) callBack();                
      });
      
      console.log('finished getResources')
    },

    
    /**
     * Note that schema now comes from resource.schema
     * 
     * A schema: (for ResourceResource)
     * 
     * schema: {
        extraSelectorOptions: {},
        fields: {
          api_name,
          comment,
          description,
          id,
          id_attribute,
          is_restricted,
          key,
          ordinal,
          resource_uri,
          scope,
          title,
          title_attribute,
          visibility
        },
        resource_definition: {}
        },
     * 
     */
    getSchema: function(resourceId) {
      return this.getResource(resourceId).schema;
    },
    
    /**
     * Get a model from the server
     */
    getModel: function(resourceId, key, callBack) {
      var self = this;
      var resource = this.getResource(resourceId);
      var url = resource.apiUri + '/' + key;

      var ModelClass = Backbone.Model.extend({
        url : url,
        defaults : {}
      });
      var instance = new ModelClass();
      instance.fetch({
          success : function(model) {
            model.resource = resource;
            model.key = key;
            callBack(model);
          },
          error: this.jqXHRerror

      });
    },
    
    
    getMenu: function(){
      
      var menu = this.get('menu');
      
      if(!this.getCurrentUser().is_superuser){
        menu.submenus = _.omit(menu.submenus, 'admin');
      }
      return menu;
    },
    
    /**
     * A resource (ResourceResource):
     * {
          key: "resource",
          scope: "resource",
          title: "Resources",
          api_name: "reports",
          comment: null,
          description: "Resource describes tables, queries that are available through the API",
          id: 680,
          id_attribute: [
            "key"
          ],
          is_restricted: null,
          ordinal: 0,
          resource_uri: "/reports/api/v1/resource/resource",
          schema: {},
          title_attribute: [
            "title"
          ],
          visibility: [
            "list",
            "detail"
          ]
        }
     *    
     */
    getResource: function(resourceId){
      var uiResources = this.get('ui_resources');
      var resources = this.get('resources');
      if(!_.has(uiResources, resourceId)) {
        if(!_.has(resources, resourceId)) {
          this.error('Unknown resource: ' + resourceId);
          throw "Unknown resource: " + resourceId;
        }
        return resources[resourceId];
      }
      return uiResources[resourceId];
    },
    
    jqXHRerror: function(model, response, options) {
      //console.log('error fetching the model: '+ model + ', response:
      // ' + JSON.stringify(response));
      var url = options.url || model.url || '';
      var msg = 'Error locating resource: ' + url;
      this.error(msg);
      
      var sep = '\n';
      if (!_.isUndefined(response.status))
          msg += sep + response.status;
      if (!_.isUndefined(response.statusText))
          msg += sep + response.statusText;
      if (!_.isEmpty(response.responseText))
          msg += sep + response.responseText;
      
      if(DEBUG) window.alert(msg);
      else console.log(msg);
      // TODO: use Bootstrap inscreen alert classed message div
    },
    
    error: function(msg){
      var msgs = this.get('messages');
      if (msgs && msgs.length > 0) {
        msgs.push(msg);
        
        if(msgs.length > 5){
          msgs = msgs.splice(4, msgs.length-5);
        }
        this.set('messages', msgs);
      }else{
        this.set('messages', [msg]);
      }
    }, 
    
    setUriStack: function(value){
      var self = this;
      var _deferred = function(){
       self.set({ uriStack: value });
      };
      if(self.isPagePending()){
        self.showModal(_deferred);
      }else{
        _deferred();
      }
    },
    
    /**
     * Set flag to signal that the current page has pending changes;
     * (see setUriStack; modal confirm dialog will be triggered).
     */
    setPagePending: function(value){
      this.set({'pagePending': value});
    },
    
    isPagePending: function(){
      return this.has('pagePending') && this.get('pagePending');
    },
    
    showModal: function(callback){
      var self = this;
      console.log('showModal: ' + JSON.stringify(this.model));
      var modalDialog = new Backbone.View({
          el: _.template(modalOkCancelTemplate, { 
            body: "Pending changes in the page: continue anyway?",
            title: "Please confirm" } ),
          events: {
              'click #modal-cancel':function(event) {
                  console.log('cancel button click event, '); 
                  event.preventDefault();
                  $('#modal').modal('hide'); 
              },
              'click #modal-ok':function(event) {
                  console.log('ok button click event, '); 
                  event.preventDefault();
                  $('#modal').modal('hide');
                  self.setPagePending(false);
                  callback();
              }
          },
      });
      modalDialog.render();
      modalDialog.$el.find('#modal-cancel').html('Cancel and return to page');
      modalDialog.$el.find('#modal-ok').html('Continue');
      $('#modal').empty();
      $('#modal').html(modalDialog.$el);
      $('#modal').modal({show:true, backdrop: 'static'});
    }   
  
  });

  var appState = new AppState();
  
  appState.resources = {};   // TO be retrieved from the server 
  
  appState.apiVersion = API_VERSION;
  appState.reportsApiUri = REPORTS_API_URI;
  appState.dbApiUri = DB_API_URI;
  appState.LIST_ARGS = ['page','rpp','order','search','log'];      
  
  
  return appState;
});