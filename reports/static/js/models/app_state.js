define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid'
], function($, _, Backbone, Iccbl ){
  
  var API_VERSION = 'api/v1';
  var REPORTS_API_URI = '/reports/' + API_VERSION;
  var DB_API_URI = '/db/' + API_VERSION;
  
  // TODO: feed AppState from the server
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

      // Define the menu data structure:
      // menu's can only be two levels deep (todo: recursive cool menu thing)
      // - the menu tree here is composed of ui_resources, defined below.
      // - the topographical structure here represents the visual structure in the app.
      menu: {
        view: 'home',
        submenus:{
          'library': {
            expanded: false,
            view: 'list',
            submenus: {
              'smallmoleculelibrary':{
                  view: 'list'
              },
              'rnalibrary':{
                  view: 'list'
              },
            }
          },
          'screensaveruser': {
            expanded: false,
            view: 'list',
            submenus: {
              'screeners':{
                view: 'list'
              },
              'staff':{
                view: 'list'
              },
            }
          },
          'screen': {
            view: 'list',
            expanded: false,
            submenus: {
              'small_molecule_screens':{
                  view: 'list',
              },
              'rnai_screens':{
                  view: 'list',
              },
            }
          },
          'admin': {
            expanded: false,
            // TODO: menu view not implemented; would show children menu items as links in the page
            view: 'menu',  
            submenus: {
              'metahash':{
                  view: 'list',
              },
              'resource':{
                  view: 'list',
              },
              'vocabularies':{
                  view: 'list',
              },
              'apilog':{
                  view: 'list',
              },
              'users':{
                  view: 'list',
              },
              'groups':{
                  view: 'list',
              },
              'permissions':{
                  view: 'list',
              },
            }
          },
        }
      },
      
      // will search ui_resources for each of the items identified as resources on the stack
      // TODO: this should also come from the API?
      ui_resources: {   
          admin: {
              title: 'Admin',
              route: '',
              view: 'HomeView',
              content_header: 'Welcome',
              description: 'Perform admin activities'
          },
          home: {
              title: 'Screensaver Reporting',
              route: '/',
              view: 'HomeView',
              content_header: 'Welcome',
              description: 'Menu starting point'
          },

          metahash: {
              header_message: 'Define fields for display on detail and list views',
              title: 'Field Information',
              route: 'list/metahash',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'metahash',
              url_root: '/reports/api/v1',
              options: { search: {'scope__exact': 'fields.metahash'}},
              description: 'Control field settings'
          },

          resource: {
              header_message: 'Define resources for display in the reporting interface',
              title: 'Resource Information',
              route: 'list/resource',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'resource',
              url_root: '/reports/api/v1',
              description: 'Control resource information'
          },

          vocabularies: {
              header_message: 'Define controlled vocabularies',
              title: 'Application Vocabularies',
              route: 'list/vocabularies',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'vocabularies',
              url_root: '/reports/api/v1',
              description: 'Enter controlled vocabularies'
          },

          apilog: {
              header_message: 'View change logs',
              title: 'Change logs',
              route: 'list/apilog',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'apilog',
              url_root: '/reports/api/v1',
              description: 'Change logs'
          },
          users: {
              header_message: 'Django users',
              title: 'Users',
              route: 'list/users',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'user',
              url_root: '/reports/api/v1',
              description: 'Django user'
          },
          groups: {
              header_message: 'User Groups',
              title: 'User Groups',
              route: 'list/groups',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'usergroup',
              url_root: '/reports/api/v1',
              description: 'User Groups'
          },
          permissions: {
              header_message: 'Django permissions',
              title: 'Permissions',
              route: 'list/permissions',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'permission',
              url_root: '/reports/api/v1',
              description: 'Permissions'
          },

          screensaveruser: {
              header_message: 'All users (Screeners and Staff)',
              title: 'Screensaver Users',
              route: 'list/screensaveruser',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'screensaveruser',
              url_root: '/db/api/v1',
              description: 'View user information'
          },

          screeners: {
              header_message: 'Screening Users',
              title: 'Screeners',
              route: 'list/screeners',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'screensaveruser',
              url_root: '/db/api/v1',
              description: 'View user information',
              options: { search: {'screeningroomuser__isnull': 'False'} }
          },

          staff: {
              header_message: 'Staff',
              title: 'Staff Users',
              route: 'list/staff',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'screensaveruser',
              url_root: '/db/api/v1',
              description: 'View user information',
              options: { search: {'administratoruser__isnull': 'False'} }
          },
          screen: {
              header_message: 'All screens (Small Molecule and RNAi)',
              title: 'Screens',
              route: 'list/screen',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'screen',
              url_root: '/db/api/v1',
              description: 'View screen information',
              options: { order_by: { facility_id:'-'} }  //TODO: instructions like this don't work
          },
          small_molecule_screens: {
              header_message: 'Small Molecule Screens',
              title: 'Small Molecule',
              route: 'list/small_molecule_screens',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'screen',
              url_root: '/db/api/v1',
              description: 'View small molecule screen information',
              options: { search: { screen_type: 'small_molecule'} }
          },
          rnai_screens: {
              header_message: 'All screens (Small Molecule and RNAi)',
              title: 'RNAi',
              route: 'list/rnai_screens',
              list_view: 'ListView',
              detailView: 'DetailView',
              api_resource: 'screen',
              url_root: '/db/api/v1',
              description: 'View rnai screen information',
              options: { search: { screen_type: 'rnai'} }
          },
          library: {
              header_message: 'All libraries (Small Molecule and RNAi)',
              title: 'Libraries',
              route: 'list/library',
              list_view: 'ListView',
              detailView: 'LibraryView',
              api_resource: 'library',
              url_root: '/db/api/v1',
              description: 'View library information'
          },
          smallmoleculelibrary: {
              header_message: 'Small Molecule Libraries',
              title: 'Small Molecule',
              route: 'list/smallmoleculelibrary',
              list_view: 'ListView',
              detailView: 'LibraryView',
              api_resource: 'library',
              url_root: '/db/api/v1',
              description: 'View Small Molecule Library information',
              options: { search: { screen_type: 'small_molecule'} }
          },
          rnalibrary: {
              header_message: 'RNAi Libraries',
              title: 'RNAi',
              route: 'list/rnalibrary',
              list_view: 'ListView',
              detailView: 'LibraryView',
              api_resource: 'library',
              url_root: '/db/api/v1',
              description: 'View RNAi library information',
              options: { search: { screen_type: 'rnai'} }
          }
      },
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
      this.set('schemas', {});
      console.log('AppState initialized');
    },
        
    start: function(callBack){
      
      this.getResources(callBack);
    },

    getResources: function(callBack){
      var self = this;
      // Retrieve the resource definitions from the server
      var resourceUrl = self.reportsApiUri + '/resource'
      Iccbl.getCollectionOnClient(resourceUrl, function(collection){
        console.log('got resources: ' + resources);
        
        // TODO: store the apiUri on the resource in the server
        var resourcesCollection = collection.toJSON();
        _.each(resourcesCollection,function(resource){
          resource.apiUri = '/' + resource.api_name + '/' + 
            self.apiVersion + '/' + resource.key;
        });

        // make a hash out of the resources, keyed by the resource id key
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

        // now combine that with the ui_resources hash
        _.each(_.keys(ui_resources), function(key){
          var ui_resource = ui_resources[key];
          // then augment the ui resources with their api_resource, if different
          if ( key !== ui_resource.api_resource ) {
            ui_resources[key] = _.extend(
                {}, resources[ui_resource.api_resource], ui_resource);
          }
        });
                
        self.set('resources', resources);
        self.set('resources', ui_resources);
        callBack();
      });
    },

    getSchema: function(resourceId, callBack) {
      var self = this;
      var schemas = this.get('schemas');
      if(_.has(schemas, resourceId)){
        callBack(schemas[resourceId]);
      } else {
        var resource = self.getResource(resourceId);
        
        var schema_url = resource.apiUri + '/schema'
        Iccbl.getSchema(schema_url, function(schema){
          schemas[resourceId] = schema;
          callBack(schema);
        });
      }
    },
    
    /**
     * Update a model with the member attributes needed for the system:
     * - resource
     * - resourceSchema
     * - key
     */
    updateModel: function(resourceId, key, model, callBack) {
      var resource = this.getResource(resourceId);
      this.getSchema(resourceId, function(schema){
        model.key = key;
        model.resource = resource;
        model.resourceSchema = schema;
        callBack(model);
      });      
    },
    
    /**
     * Get a model from the server
     */
    getModel: function(resourceId, key, callBack) {
      
      var resource = this.getResource(resourceId);
      this.getSchema(resourceId, function(schema){
        var url = resource.apiUri + '/' + key;

        var ModelClass = Backbone.Model.extend({
          url : url,
          defaults : {}
        });
        var instance = new ModelClass();
        instance.fetch({
            success : function(model) {
              model.resourceSchema = schema;
              model.resource = resource;
              model.key = key;
              callBack(model);
            },
            error : function(model, response, options) {
                //console.log('error fetching the model: '+ model + ', response:
                // ' + JSON.stringify(response));
                var msg = 'Error locating resource: ' + url;
                var sep = '\n';
                if (!_.isUndefined(response.status))
                    msg += sep + response.status;
                if (!_.isUndefined(response.statusText))
                    msg += sep + response.statusText;
                if (!_.isEmpty(response.responseText))
                    msg += sep + response.responseText;
                window.alert(msg);
                // TODO: use Bootstrap inscreen alert classed message div
            }
        });
      });
    },
    
    
    getResource: function(resourceId){
      var uiResources = this.get('ui_resources');
      var resources = this.get('resources');
      if(!_.has(uiResources, resourceId)) {
        if(!_.has(resources, resourceId)) {
          window.alert('Unknown resource: ' + resourceId);
          throw "Unknown resource: " + resourceId;
        }
        return resources[resourceId];
      }
      var uiResource = uiResources[resourceId];
      return uiResource;
//      // test for nested resource
//      if (_.has(uiResource, 'api_resource') 
//          && resourceId !=== uiResource.api_resource ){
//        
//      }
//      return this.getResource(uiResources[resourceId].api_resource);
//      if (!_.has(resources, resourceId)){
//      }
//      return ;
    }
  });

//  return AppState;
  var appState = new AppState();
  
  
  

  appState.resources = {};   // TO be retrieved from the server  
  
  appState.apiVersion = API_VERSION;
  appState.reportsApiUri = REPORTS_API_URI;
  appState.dbApiUri = DB_API_URI;
  appState.LIST_ARGS = ['page','rpp','order','search'];      
  
  
  return appState;
});