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
      _.bindAll(this,'jqXHRerror', 'jqXHRFail', 'error');
    },
        
    start: function(callBack) {
      var self = this;
      console.log('start app_state.js');
      this.getResources(function(){
        self.getVocabularies(function(vocabularies){
          self.set('vocabularies', vocabularies);
          self.setCurrentUser(callBack);
        });
      });
    },
    
    setCurrentUser: function(callBack) {
      var self = this;
      this.getModel('user', window.user, function(model){
        self.currentUser = model.toJSON();
        console.log('Current user: ' + JSON.stringify(self.currentUser));
        
        if(callBack) callBack(); 
      });
    },
    
    getCurrentUser: function(){
      return this.currentUser;
    },
    
    getVocabularies: function(callBack){
      console.log('getVocabularies from the server...');
      var self = this;
      
      var resourceUrl = self.reportsApiUri + '/vocabularies'
      Iccbl.getCollectionOnClient(resourceUrl, function(collection){
        var vocabularies = {};
        collection.each(function(vModel){
          var v = vModel.toJSON();
          if(!_.has(vocabularies, v.scope)){
            vocabularies[v.scope] = {};
          }
          vocabularies[v.scope][v.key]=v;
        });
        callBack(vocabularies);
      });
    },
    
    getVocabulary: function(scope){
      if(!this.has('vocabularies')){
        throw "Vocabularies aren't initialized";
      }
      var vocabularies = this.get('vocabularies');
      if(!_.has(vocabularies, scope)) {
          throw "Unknown vocabulary: " + scope;
      }
      return vocabularies[scope];
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
            try{
              callBack(model);
            } catch (e) {
              console.log('uncaught error: ' + e);
              self.error('error displaying model: ' + model + ', '+ e);
            }
          },
          error: this.jqXHRerror

      });
    },
    
    
    getMenu: function(){
      var self = this;
      var currentUser = this.getCurrentUser();
      var menu = this.get('menu');
      
      if(!currentUser.is_superuser){
        // TODO: use permissions here
        menu.submenus = _.omit(menu.submenus, 'admin');
        
        // TODO: iterate over each menu: if user doesn't have read perm for 
        // resource, omit
        var new_submenus = {}
        _.each(_.keys(menu.submenus), function(key){
          if(self.hasPermission(key)){
            new_submenus[key] = menu.submenus[key];
          }else{
            console.log('user: ' + currentUser.username 
                + ', doesnt have permission to view the menu: ' + key );
          }
//          var r_perm = 'permission/resource/'+key
//          var match = false;
//          _.each(self.getCurrentUser().all_permissions, function(permission){
//            // NOTE: allow match of [either] 'read' or 'write'
//            if(permission.indexOf(r_perm) > -1 ) {
//              match = true;
//            }
//          })
//          if(match){
//            new_submenus[key] = menu.submenus[key];
//          }else{
//            console.log('user: ' + currentUser.username 
//                + ', doesnt have permission to view the menu: ' + key );
//          }
        });
        menu.submenus = new_submenus;
      }
      return menu;
    },
    
//    hasPermission: function(p_check){
//      var self = this;
//      var current_user = self.getCurrentUser();
//      if (current_user.is_superuser){
//        return true;
//      }
//      var result = _.find(current_user.all_permissions, function(permission){
//        console.log('check for perm: ' + permission);
//        return permission.indexOf(p_check) > -1;
//      });
//      console.log("perm check: " + p_check + ',' + result);
//      return result;
//    },
    
    /**
     * Test if the current user has the resource/permission - 
     * - if permission is unset, 
     * will check if the user has *any* permission on the resource.
     */
    hasPermission: function(resource, permission){
      
      var self = this;
      if(self.getCurrentUser().is_superuser) return true;
      
      var r_perm = 'permission/resource/'+ resource;
      if(!_.isUndefined(permission)){
        r_perm += '/'+ permission;
      }// otherwise, will return true if user has either permission
      var match = _.find(
          self.getCurrentUser().all_permissions, 
          function(p){
            if(p.indexOf(r_perm) > -1 ) {
              return true;
            }
          });
      return !_.isUndefined(match);
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
      if(!_.has(uiResources, resourceId)) {
          throw "Unknown resource: " + resourceId;
      }
      return uiResources[resourceId];
    },
        
    getResourceFromUrl: function(resourceId, schemaUrl, callback){
      
      var uiResources = this.get('ui_resources');
      var ui_resource = {};
      if(_.has(uiResources, resourceId)) {
        ui_resource = uiResources[resourceId];
      }
      
      var ModelClass = Backbone.Model.extend({
        url : schemaUrl,
        defaults : {}
      });
      var instance = new ModelClass();
      instance.fetch({
          success : function(model) {
            console.log('resource schema model', model.toJSON());
            schema = model.toJSON();
            ui_resource = _.extend({}, ui_resource, schema.resource_definition);
            var schemaClass = new SchemaClass();
            ui_resource.schema = _.extend(schema, schemaClass);
            
            callback(ui_resource)
          },
          error: this.jqXHRerror
      });      
      
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
    
    jqXHRFail: function(xhr, text, message){
      var self = this;
      
      var msg = message;
      if ( _.has(xhr,'responseJSON') && !_.isEmpty(xhr.responseJSON) ) {
        
        if ( _.has( xhr.responseJSON, 'error_message') ) {
          msg += ': ' + xhr.responseJSON.error_message;
        }
        else if ( _.has( xhr.responseJSON, 'error') ) {
          msg += ': ' + xhr.responseJSON.error;
        }else{
          msg += xhr.responseText;
        }
        
      } else {
        
        var re = /([\s\S]*)Request Method/;
        var match = re.exec(xhr.responseText);
        
        if (match) {
          msg += ': ' + match[1] + ': ' + xhr.status + ':' + xhr.statusText;
        } else {
        
          try{
            var rtext = JSON.parse(xhr.responseText);
            msg += ':' + rtext.error_message;
          } catch (e) {
            console.log('couldnt parse responseText: ' + xhr.responseText)
            msg += ': ' +  xhr.status + ':' + xhr.statusText + ', ' + xhr.responseText;
          }
        }              
        
      }
      self.error(msg);
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
      this.set({ uriStack: value });
    },

//    setUriStack: function(value){
//      var self = this;
//      var _deferred = function(){
//       self.set({ uriStack: value });
//      };
//      
//      // TODO: push the requestPageChange method up the call stack to the client code:
//      // - this will allow us to prevent other actions on the client side.
//      
//      if(self.isPagePending()){
//        self.showModal({
//          ok: _deferred
//        });
//      }else{
//        _deferred();
//      }
//    },
    
    /**
     * Set flag to signal that the current page has pending changes;
     * (see setUriStack; modal confirm dialog will be triggered).
     */
    setPagePending: function(callback){
      this.set({'pagePending': callback});
    },
    clearPagePending: function(){
      this.unset('pagePending');
    },    
    isPagePending: function(){
      return this.has('pagePending');
    },
    
    /**
     * options.ok = ok function
     * options.cancel = cancel function
     */
    requestPageChange: function(options){
      var self = this;
      var callbackOk = options.ok;
      options.ok = function(){
        if(callbackOk) callbackOk();
        self.clearPagePending();
      };
      

      if(! self.isPagePending()){
        options.ok();
      }else{
        var pendingFunction = this.get('pagePending');
        if(_.isFunction(pendingFunction)){
          options.title = 'Please confirm';
          options.body = "Pending changes in the page: continue anyway?";
          options.cancel = pendingFunction;
        }
        self.showModal(options);
      }
    },
    
    /**
     * options.ok = ok function
     * options.cancel = cancel function
     * options.body
     * options.title
     */
    showModal: function(options){
      
      var self = this;
      var callbackOk = (options && options.ok)? options.ok : function(){};
      var callbackCancel = (options && options.cancel)? options.cancel: function(){};
          
      console.log('showModal: ' + options);
      var modalDialog = new Backbone.View({
          el: _.template(modalOkCancelTemplate, { 
            body: options.body,
            title: options.title } ),
          events: {
              'click #modal-cancel':function(event) {
                  console.log('cancel button click event, '); 
                  event.preventDefault();
                  $('#modal').modal('hide'); 
                  callbackCancel();
              },
              'click #modal-ok':function(event) {
                  console.log('ok button click event, '); 
                  event.preventDefault();
                  $('#modal').modal('hide');
                  self.clearPagePending();
                  callbackOk();
              }
          },
      });
      modalDialog.render();
      if(!_.isUndefined(options.view)){
        modalDialog.$el.find('.modal-body').append(options.view.render().el);
      }
      modalDialog.$el.find('#modal-cancel').html('Cancel and return to page');
      modalDialog.$el.find('#modal-ok').html('Continue');
      $('#modal').empty();
      $('#modal').html(modalDialog.$el);
      $('#modal').modal({show:true, backdrop: 'static'});
    }   
  
  });

  var appState = new AppState();
  
  appState.schemaClass = new SchemaClass(); // make accessible to outside world
  
  appState.resources = {};   // TO be retrieved from the server 
  
  appState.apiVersion = API_VERSION;
  appState.reportsApiUri = REPORTS_API_URI;
  appState.dbApiUri = DB_API_URI;
  appState.LIST_ARGS = ['page','rpp','includes','order','search','log', 'children'];      
  
  Iccbl.appModel = appState;
  
  return appState;
});