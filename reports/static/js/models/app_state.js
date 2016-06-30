define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid',
  'templates/modal_ok_cancel.html'    
], function($, _, Backbone, Iccbl, modalOkCancelTemplate ){
  
  var API_VERSION = 'api/v1';
  var REPORTS_API_URI = '/reports/' + API_VERSION;
  var DB_API_URI = '/db/' + API_VERSION;
  var SEARCH_DELIMITER = ';';
  var DEBUG = false;
  
  var SchemaClass = Iccbl.SchemaClass = function() {};
  SchemaClass.prototype.detailKeys = function()
  {
    return this.filterKeys('visibility','d');
  };  
  SchemaClass.prototype.editVisibleKeys = function()
  {
    return this.filterKeys('visibility','e');
  };  
  SchemaClass.prototype.allEditVisibleKeys = function()
  {
    var keys =  _.union(
      this.editVisibleKeys(),
      this.createKeys(),
      this.updateKeys());
    
    return Iccbl.sortOnOrdinal(keys,this.fields);
  };  
  SchemaClass.prototype.createKeys = function()
  {
    return this.filterKeys('editability','c');
  };  
  SchemaClass.prototype.updateKeys = function()
  {
    return this.filterKeys('editability','u');
  };  
  SchemaClass.prototype.filterKeys = function(select_field, visibility_term)
  {
    var self = this;
    var keys = Iccbl.sortOnOrdinal(
      _.keys(self.fields), self.fields)
    var detailKeys = _(keys).filter(function(key){
      return _.has(self.fields, key) && 
          _.has(self.fields[key], select_field) && 
          _.contains(self.fields[key][select_field], visibility_term);
    });
    return detailKeys;
  };  
  
  SchemaClass.prototype.groupedKeys = function(keys){
    var self = this;
    var groupedKeys = [];
    var groups = {};
    _.each(keys,function(key){
      var fi = self.fields[key];
      if (fi.display_options && fi.display_options.group){
        if (!_.has(groups,fi.display_options.group)){
          var newGroup = {
            title: fi.display_options.group,
            fields: []
          }
          groups[fi.display_options.group] = newGroup;
          groupedKeys.push(newGroup); 
        }
        group = groups[fi.display_options.group];
        group.fields.push(key);
      }else{
        groupedKeys.push(key);
      }
    });
    return groupedKeys;
  };

  var AppState = Backbone.Model.extend({
    
    defaults: {

      // TODO: deprecate these variables
      // use the REPORTS_API_URI, and DB_API_URI defined below
      root_url: '/lims',  // used for the backbone history
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
      _.bindAll(this,'error',
        'setCurrentUser','getResources','getVocabularies',
        'getAdminUserOptions','getUserOptions','getUserGroupOptions',
        'getPrincipalInvestigatorOptions','getLibraries','getLibraryOptions');
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
    
    /**
     * perform lazy fetch for Admin data
     */
    initializeAdminMode: function(callback) {
      var self = this;
      console.log('start app_state.js');
      // TODO: use jquery.queue() as in libraryscreening.js
      self.getAdminUserOptions(function(options){
        self.getUserOptions(function(options){
          self.getPrincipalInvestigatorOptions(function(options){
            self.getUserGroupOptions(function(options){
              callback();
            });
          });
        });
      });
    },
    
    setCurrentUser: function(callBack) {
      console.log('set the current user', window.user);
      var self = this;
      this.getModel('user', window.user, 
        function(model){
          console.log('user model:', model);
          self.currentUser = model.toJSON();
          console.log('Current user: ' + JSON.stringify(self.currentUser));
          
          if(callBack) callBack(); 
        }
      );
    },
    
    getCurrentUser: function(){
      return this.currentUser;
    },
    
    setSearch: function(searchId, search_object){
      console.log('setSearch', searchId, search_object)
      localStorage.setItem(''+searchId, JSON.stringify(search_object));
      this.getSearch(searchId);
    },
    
    getSearch: function(searchId){
      var obj = localStorage.getItem(''+searchId);
      if(obj) obj = JSON.parse(obj);
      console.log('getSearch', searchId, obj);
      return obj;
    },

    getLibraryOptions: function(callBack){
      var self = this;
      var prop = 'libraryOptions';
      var options = this.get(prop);
      if(!options){
        this.getLibraries(function(libraries){
          options = [];
          libraries.each(function(library){
            var short_name = library.get('short_name');
            var library_name = library.get('library_name');
            options.push({ val: short_name, label: library_name });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
      }
      return options;
    },

    getScreeningLibraryOptions: function(screen_type, callBack){
      var self = this;
      var prop = 'screeningLibraryOptions-'+screen_type;
      var options = this.get(prop);
      if(!options){
        this.getLibraries(function(libraries){
          options = [];
          libraries.each(function(library){
            if (_.contains(
                ['not_allowed','retired','not_yet_plated'],
                library.get('screening_status'))){
              return;
            }
            if (library.get('screen_type') != screen_type){
              return;
            }
            var short_name = library.get('short_name');
            var library_name = library.get('library_name');
            options.push({ val: short_name, label: library_name });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
      }
      return options;
    },

    /**
     * Generate a list of "options" suitable for use in a user multiselect.
     * [ { val: username, label: name:username }, ... ]
     */
    getAdminUserOptions: function(callBack){
      var self = this;
      var prop = 'adminUserOptions';
      var options = this.get(prop);
      if(!options){
        this.getAdminUsers(function(users){
          var options = [];
          users.each(function(user){
            var username = user.get('username');
            var name = user.get('name');
            options.push({ val: username, label: name + ':' + username });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
      }
      return options;
    },
    
    getUserOptions: function(callBack){
      var self = this;
      var prop = 'userOptions';
      var options = this.get(prop);
      if(!options){
        self.getUsers(function(users){
          var options = [];
          users.each(function(user){
            var username = user.get('username');
            var name = user.get('name');
            options.unshift({ val: username, label: name + ':' + username });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
      }
      return options;
    },
    getPrincipalInvestigatorOptions: function(callBack){
      var self = this;
      var prop = 'piOptions';
      var options = this.get(prop);
      if(!options){
        this.getPrincipalInvestigators(function(users){
          var options = [{ val:'',label:'' }];
          users.each(function(user){
            var username = user.get('username');
            var lab_name = user.get('lab_name');
            options.push({ val: username, label: lab_name });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
      }
      return options;
    },

    getUserGroupOptions: function(callBack){
      var self = this;
      var prop = 'userGroupOptions';
      var options = this.get(prop);
      if(!options){
        this.getUserGroups(function(usergroups){
          var options = [];
          usergroups.each(function(usergroup){
            var name = usergroup.get('name');
            options.unshift({ val: name, label: name });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
      }
      return options;
    },

    getPermissionsOptions: function(){
      return this.get('permissionOptions');
    },
    
    getLibraries: function(callback){
      data_for_get = { 
        exact_fields: [
           'short_name','library_name','start_plate','end_plate','copies',
           'screening_copies',
           'library_type','screening_status','screen_type'], 
        order_by: ['short_name']
      };
      return this.getCachedResourceCollection(
        'libraries', this.dbApiUri + '/library', data_for_get, callback );
    },

    getCachedResourceCollection: function(
        prop,
        url,
        data_for_fetch,
        callback,
        options
        ){
      var self = this;
      options = options || {};
      data_for_fetch = _.extend({ limit: 0 }, data_for_fetch);
      var flush = options.flush || false;
      var failCallback = options.failCallback;
      //, timeout - TODO

      var collection = this.get(prop);
      if(flush || _.isEmpty(collection)){
        console.log('get all '+prop+' from the server...');
        var CollectionClass = Iccbl.CollectionOnClient.extend({
          url: url 
        });
        var instance = new CollectionClass();
        instance.fetch({
          data: data_for_fetch,
          success: function(collection, response) {
            self.set(prop, collection);
            if (callback) callback(collection);
          }
        }).fail(function(){ 
          if (failCallback){
            failCallback.apply(this,arguments);
          }else{
            Iccbl.appModel.jqXHRfail.apply(this,arguments); 
          }
        });
      }
      else{
        if (callback) callback(collection);
      }
      return this.get(prop);
    },
    
    getPrincipalInvestigators: function(callback) {
      return this.getCachedResourceCollection(
        'principal_investigators', 
        this.dbApiUri + '/screensaveruser', 
        { 
          lab_head_affiliation__is_blank: false,
          classification__eq: 'principal_investigator',
          exact_fields: ['username','lab_name','name','email'], 
          order_by: ['name']
        }, 
        callback );
    },
    
    getUsers: function(callback) {
      return this.getCachedResourceCollection(
        'users', 
        this.dbApiUri + '/screensaveruser', 
        { 
          exact_fields: ['username','name','email'], 
          order_by: ['name']
        }, 
        callback );
    },
    
    getAdminUsers: function(callback) {
      return this.getCachedResourceCollection(
        'adminUsers', 
        this.dbApiUri + '/screensaveruser', 
        { 
          is_staff__eq: true,
          exact_fields: ['username','name','email'], 
          order_by: ['name']
        }, 
        callback );
    },

    getUserGroups: function(callback) {
      return this.getCachedResourceCollection(
        'usergroups', 
        this.reportsApiUri + '/usergroup', 
        { 
          order_by: ['name']
        }, 
        callback );
    },
    
    /**
     * Create a vocabulary hash, from the server:
     * { 
     *    v.scope: { 
     *      v.key : { scope: , key: , title: , ordinal: }
     *    },
     *    v.scope1: {},
     *    etc.,...
     */
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
        console.log('getVocabularies: done');
        callBack(vocabularies);
      });
    },    
    
    /**
     * Cache an "options" array for all permissions, for the editor UI
     * "options" are:
     * { val: 'key', label: 'label' }
     */
    setPermissionsOptions: function(resources){
      var self = this;
      var permissionOptions = [];
      var optionGroups = { resources: []};
      _.each(_.keys(resources), function(rkey){
        var resource = resources[rkey];
        if(!_.has(resource,'fields')){
          // skip these, only consider server side resources that have fields
          return;
        }
        optionGroups['resources'].unshift({ 
          val: ['resource',rkey,'read'].join('/'), 
          label: ['resource',rkey,'read'].join(':'), 
        });
        optionGroups['resources'].unshift({ 
          val: ['resource',rkey,'write'].join('/'), 
          label: ['resource',rkey,'write'].join(':'), 
        });
        var fieldGroup = rkey;
        var fields = resource['fields'];
        if(!_.has(optionGroups,fieldGroup)){
          optionGroups[fieldGroup] = [];
        }
        _.each(_.keys(fields),function(fkey){
          optionGroups[fieldGroup].unshift({
            val: ['fields.' + fieldGroup,fkey,'read'].join('/'),
            label: [fieldGroup,fkey,'read'].join(':')
          });
        });
      });
      _.each(_.keys(optionGroups),function(key){
        permissionOptions.unshift({ group: key, options: optionGroups[key] });
      });
      self.set('permissionOptions',permissionOptions);
    },

    
    /** 
     * Return an array of options for a vocabulary select.
     */
    getVocabularySelectOptions: function(scope){
      choiceHash = []
      try{
        var vocabulary = this.getVocabulary(scope);
        _.each(_.keys(vocabulary),function(choice){
          if(vocabulary[choice].is_retired){
            console.log('skipping retired vocab: ',choice,vocabulary[choice].title );
          }else{
            choiceHash.push({ val: choice, label: vocabulary[choice].title });
          }
        });
      }catch(e){
        var msg = 'Vocabulary unavailable: vocabulary_scope_ref: ' + scope;
        console.log(msg,e);
        this.error(msg);
      }
      return choiceHash;
    },
    
    /**
     * return the set of vocabulary items for a scope:
     *    v.scope: { 
     *      v.key1 : { scope: , key: , title: , ordinal: },
     *      v.key2 : {},
     *      etc.,...
     *    },
     * @param scope - scope to search for
     * - "scope" can also be a regex, and will be matched to all scopes using
     * String.prototype.match(candidateScope,^scope$)
     */
    getVocabulary: function(scope){
      if(!this.has('vocabularies')){
        throw "Vocabularies aren't initialized";
      }
      var vocabularies = this.get('vocabularies');
      if(!_.has(vocabularies, scope)) {
        // test for regex match/matches
        var matchedVocabularies = {};
        _.each(_.keys(vocabularies), function(candidateScope){
          if(candidateScope.match('^' + scope + '$')){
            console.log('matching: ' + '^' + scope + '$' + ', to: ' + candidateScope );
            _.extend(matchedVocabularies,vocabularies[candidateScope]);
          }
        });
        if(!_.isEmpty(matchedVocabularies)){
          console.log('matchedVocabularies', scope, matchedVocabularies );
          return matchedVocabularies;
        }
        throw "Unknown vocabulary: " + scope;
      }
      return vocabularies[scope];
    },
    
    /**
     * return the title for the vocabulary entry for the specified scope and val.
     */
    getVocabularyTitle: function(scope,val){
      try{
        vocabulary = this.getVocabulary(scope);
        if(_.has(vocabulary,val)){
          val = vocabulary[val].title;
        }else{
          var msg = Iccbl.formatString(
            'vocabulary: {vocabulary} is misconfigured: rawData: {rawData}',
            { vocabulary: _.result(this, "vocabulary_scope_ref"),
              rawData: val 
            });
          console.log(msg);
          this.error(msg);
        }
      }catch(e){
        this.error('unknown vocabulary: ' + scope);
      }
      return val;
    },

    parseSchemaField: function(rawField) {
      if(_.has(rawField,'display_options') && ! _.isEmpty(rawField.display_options)){
        var options= rawField['display_options'];
        options = options.replace(/'/g,'"');
        try{
          rawField['display_options'] = JSON.parse(options);
        }catch(e){
          console.log('warn: display_options is not JSON parseable, column: ' +
              rawField['key'] + ', options: ', options,e);      
        }
      } 
    },
    
    getResources: function(callBack){
      console.log('getResources from the server...')
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
          resource = _.extend(resource, schemaClass);
        });

        // 1. Create a resource hash, keyed by the resource id key
        // 2. Augment uiResources with the (api)resources
        var ui_resources = self.get('ui_resources');
        var resources = {};
        _.each(resourcesCollection, function(resource){
          resources[resource.key] = resource;
          _.each(_.values(resource.fields), self.parseSchemaField );
          
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
        
        // set up permissions for all the resources
        self.setPermissionsOptions(self.get('ui_resources'));
        
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
      return this.getResource(resourceId);
    },
    
    createNewModel: function(resourceId, defaults) {
      console.log('create new model for: ',resourceId, defaults );
      var self = this;
      var resource = self.getResource(resourceId);
      var defaults = defaults || {};
      _.each(resource.allEditVisibleKeys(), function(key){
        var field = resource.fields[key];
        if (key == 'resource_uri') {
          // TODO: not used
          defaults[key] = resource.apiUri; 
        } else if (key == 'id'){
          // nop always exclude the id field to signal create case to the server
        } else {
          if(! _.has(defaults, key)){
            if(field.default && !_.isNull(field.default)){
              try {
                if (field.data_type == 'string'){
                  defaults[key] = field.default;
                }else if (field.data_type == 'boolean'){
                  defaults[key] = JSON.parse(field.default.toLowerCase());
                }else{
                  // TODO: test all types here
                  defaults[key] = JSON.parse(field.default);
                }
              }catch(e){
                console.log('json parse error', key, field.default);
                if (field.data_type == 'date' && field.default == 'now'){
                  defaults[key] = Iccbl.getISODateString(new Date());
                }else{
                  self.error('Warning, unparseable default for field: ' 
                    + field.key + ', value: ' + field.default );
                  defaults[key] = '';
                }
              }
            }else{
              defaults[key] = ''; // placeholders for the edit form
            }
          }
        }
      });

      var NewModel = Backbone.Model.extend({
        urlRoot: resource.apiUri , defaults: defaults 
      });
      var newModel = new NewModel();
      newModel.resource = resource;
      console.log('new model', newModel);
      
      return newModel;
    },
    
    /**
     * Get a model from the server
     */
    getModel: function(resourceId, key, callBack, options) {
      console.log('getModel', resourceId, key);
      var self = this;
      var options = options || {};
      var failCallback = options.failCallback;
      var data_for_get = _.extend({ includes: '*' }, options.data_for_get );
      
      var resource = this.getResource(resourceId);
      if(_.isArray(key)){
        key = key.join('/');
      }
      var url = resource.apiUri + '/' + key;
      console.log('fetch model', url);
      var ModelClass = Backbone.Model.extend({
        url : url,
        defaults : { },
        parse: function(resp, options){
          // workaround for if the server returns the object in an "objects" array
          if(_.has(resp,'objects') && _.isArray(resp.objects)
              && resp.objects.length == 1 ){
            resp = resp.objects[0];
          }
          return resp;
        }
      });
      var instance = new ModelClass();
      instance.fetch({
        data: data_for_get,
        success : function(model) {
          console.log('model retrieved:', model);
          model.resource = resource;
          model.key = key;
          try{
            callBack(model);
          } catch (e) {
            console.log('uncaught error: ' + e);
            self.error('error displaying model: ' + model + ', '+ e);
          }
        }
      }).fail(function(){ 
        console.log('fail...', arguments, failCallback);
        if (failCallback){
          failCallback.apply(this,arguments);
        }else{
          Iccbl.appModel.jqXHRfail.apply(this,arguments); 
        }
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
        });
        menu.submenus = new_submenus;
      }
      return menu;
    },
        
    /**
     * Test if the current user has the resource/permission - 
     * - if permission is unset, 
     * will check if the user has *any* permission on the resource.
     */
    hasPermission: function(resource, permission){
      if (DEBUG) console.log('hasPermission', resource, permission);
      var self = this;
      if(self.getCurrentUser().is_superuser) return true;
      
      var r_perm = 'resource/'+ resource;
      if( permission=='edit'){
        permission = 'write';
      }
      if(!_.isUndefined(permission) && permission == 'write' ){
        // only check for 'write' permission; write implies read
        r_perm += '/'+ permission;
      }// otherwise, will return true if user has either permission
      var match = _.find(
          self.getCurrentUser().all_permissions, 
          function(p){
            if(p.indexOf(r_perm) > -1 ) {
              return true;
            }
          });
      console.log('permission', r_perm, !_.isUndefined(match));
      return !_.isUndefined(match);
    },
    
    
    /**
     * @return the SmallMoleculeUserAgreementLevel (group) assigned to this user
     */
    getSMUALevel: function(user){
      return _.find(user.get('usergroups'), function(usergroup){
        if(usergroup.toLowerCase().indexOf('smdsl') > -1) return true;
      });
    },
    /**
     * @return the SmallMoleculeUserAgreementLevel (group) assigned to this user
     */
    getRNAUALevel: function(user){
      return _.find(user.get('usergroups'), function(usergroup){
        if(usergroup.toLowerCase().indexOf('rnaidsl') > -1) return true;
      });
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
            "l",
            "d"
          ]
        }
     *    
     */
    getResource: function(resourceId){
      var uiResources = this.get('ui_resources');
      if(!_.has(uiResources, resourceId)) {
        var msg = "Unknown resource: " + resourceId;
        Iccbl.appModel.error(msg);
      }
      return uiResources[resourceId];
    },
        
    getResourceFromUrl: function(resourceId, schemaUrl, callback){
      var self = this;
      var uiResources = this.get('ui_resources');
      var ui_resource = {};
      if(_.has(uiResources, resourceId)) {
        ui_resource = uiResources[resourceId];
      }
      
      if (_.isEmpty(ui_resource) 
          || _.contains(['screenresult','well'], resourceId)){
        // Re-fetch the specific resource schema from the Resource endpoint:
        // - if the Resource has customizations of the schema
        // - i.e. "extraSelectorOptions
        
        var ModelClass = Backbone.Model.extend({
          url : schemaUrl,
          defaults : {}
        });
        var instance = new ModelClass();
        instance.fetch({
            success : function(model) {
              console.log('resource schema model', model.toJSON());
              schema = model.toJSON();
              _.each(_.values(schema.fields), self.parseSchemaField );
  
              ui_resource = _.extend({}, ui_resource, schema);
              var schemaClass = new SchemaClass();
              ui_resource = _.extend(ui_resource, schemaClass);
              
              callback(ui_resource)
            }
        }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      } else {
        callback(ui_resource);
      }      
    },
    
    /**
     * Process a jQuery.jqXHR.fail callback for the ajax() call.
     */
    jqXHRfail: function(jqXHR, textStatus, errorThrown){
      console.log('jqXHRfail', textStatus, errorThrown);
      var msg = textStatus || 'Error';
      if (errorThrown){
        msg += ': ' + errorThrown; 
      }
      if (this.url){
        // * url will be set for most backbone objects doing the fetching
        msg += ':\n ' + this.url;
      }
      if (jqXHR){
        if (jqXHR.status){
          msg += ':\n status:' + jqXHR.status;
        }
      }
      $(document).trigger('ajaxComplete');
      if (jqXHR && _.has(jqXHR,'responseJSON') && !_.isEmpty(jqXHR.responseJSON) ) {
        Iccbl.appModel.showJsonMessages(jqXHR.responseJSON);
      } else {
        Iccbl.appModel.error(msg);
      }
    },
    
    /**
     * Show a JSON object in a modal dialog:
     * - transform the object into a table using a depth-first traversal:
     * - each row contains the nodes traversed to each leaf node.
     */
    showJsonMessages: function(jsonObj){
      
      function dict_to_rows(dict){
        var rows = [];
        if (_.isObject(dict) && !_.isArray(dict)){
          _.each(_.keys(dict), function(key){
            _.each(dict_to_rows(dict[key]),function(row){
              if (_.isEmpty(row)){
                rows.push(key);
              }else{
                var keyrow = [key];
                if (!_.isArray(row)){
                  keyrow.push(row);
                }else{
                  keyrow = keyrow.concat(row);
                }
                rows.push(keyrow);
              }
            });
          });
        }else{
          console.log('dict: ', dict);
          if (_.isArray(dict)){
            return dict;
          }else{
            return [dict];
          }
        }
        return rows;
      }
      
      var title = "Messages";
      if(_.keys(jsonObj).length == 1){
        title = _.keys(jsonObj)[0];
        jsonObj = jsonObj[title];
      }
      
      var msg_rows = dict_to_rows(jsonObj);
      var bodyMsg = msg_rows;
      if (_.isArray(msg_rows) && msg_rows.length > 1){
        bodyMsg = _.map(msg_rows, function(msg_row){
          return msg_row.join(' - ');
        }).join('<br>');
      }
      
      Iccbl.appModel.showModalMessage({
        body: bodyMsg,
        title: title  
      });
      
    },
   
    error: function(msg){
      console.log('error: ', msg);
      var msgs = this.get('messages');
      if (msgs && msgs.length > 0) {
        msgs.push(msg);
        
        if(msgs.length > 5){
          msgs = msgs.splice(4, msgs.length-5);
        }
        // FIXME: consider a model attribute on app_state for messages, as this
        // pattern is needed for additions to the array
        this.set({'messages': msgs},{silent: true} );
        this.trigger('change:messages');
      }else{
        this.set('messages', [msg]);
      }
    }, 
    
    setUriStack: function(value){
      if(this.get('uriStack') == value ){
        // signal a change event if this method was called with the same value
        // - for the menu signaling
        this.trigger('change:uriStack', this);
      }else{
        this.set({ uriStack: value });
      }
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
      var val = callback || true;
      this.set({'pagePending': val});
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
      var options = options || {};
      var self = this;
      var callbackOk = options.ok;
      options.ok = function(){
        if(callbackOk) callbackOk();
        self.clearPagePending();
      };
      if(! self.isPagePending()){
        options.ok();
      }else{
        options.title = 'Please confirm';
        options.body = "Pending changes in the page: continue anyway?";
        var pendingFunction = this.get('pagePending');
        if(_.isFunction(pendingFunction)){
          options.cancel = pendingFunction;
        }
        self.showModal(options);
      }
    },
    
    /** Add a vocabulary term to the editForm & to the server:
     * @param vocabulary_scope_ref
     * @param name or title of the vocabulary
     * @param callback(new_vocabulary_item) callback recieves the new vocab term.
     **/
    addVocabularyItemDialog: function(
        vocabulary_scope_ref, vocabulary_name, callback, options){
      
      var self = this;
      var options = options || {};
      var description = options.description || 'Enter a ' + vocabulary_name;
      var choiceHash = {}
      var currentVocab, vocabulary;
      var formSchema = {};
      var formKey = vocabulary_scope_ref;  // name of the form field being added

      try{
        currentVocab = self.getVocabulary(vocabulary_scope_ref);
      }catch(e){
        console.log('on get vocabulary', e);
        self.error('Error locating vocabulary: ' + vocabulary_scope_ref);
        return;
      }
      formSchema[formKey] = {
        title: vocabulary_name,
        key: formKey,
        editorAttrs: { placeholder: description },
        type: 'Text',
        validators: ['required'],
        template: self._field_template 
      };
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        type: 'TextArea',
        template: self._field_template
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          var errs = {};
          var newVal = attrs[formKey];
          if (newVal){
            newVal = newVal.toLowerCase().replace(/\W+/g, '_');
            if(_.has(currentVocab,newVal)){
              errs[formKey] = '"'+ attrs[formKey] + '" is already used';
            }
          }
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: self._form_template
      });
      var _form_el = form.render().el;

      var dialog = self.showModal({
        okText: 'Create',
        view: _form_el,
        title: 'Create a new ' + vocabulary_name,
        ok: function(e){
          e.preventDefault();
          var errors = form.commit({ validate: true }); // runs schema and model validation
          if(!_.isEmpty(errors) ){
            _.each(_.keys(errors), function(key){
              $('[name="'+key +'"]').parents('.form-group,.input-group').addClass('has-error');
              $('[name="'+key +'"]').parents('.form-group,.input-group').append(errors[key].message);
            });
            return false;
          }else{
            var values = form.getValue();
            var resource = self.getResource('vocabularies');
            var key = values[formKey].toLowerCase().replace(/\W+/g, '_');
            var ordinal = currentVocab.length + 1;
            var max_item = _.max(currentVocab, function(item){ return item.ordinal });
            if (max_item){
              ordinal = max_item.ordinal + 1;
            }
            
            var data = {
              'scope': vocabulary_scope_ref,
              'key': key,
              'title': values[formKey],
              'description': values[formKey],
              'ordinal': ordinal,
              'comment': values['comment']
            };
            
            $.ajax({
              url: resource.apiUri,    
              data: JSON.stringify(data),
              contentType: 'application/json',
              method: 'POST',
              success: function(data){
                self.getVocabularies(function(vocabularies){
                  self.set('vocabularies', vocabularies);
                  callback(data);
                });
                self.showModalMessage({
                  title: 'New ' + vocabulary_name + ' created',
                  okText: 'Ok',
                  body: '"' + values[formKey] + '"',
                  ok: function(e){
                    e.preventDefault();
                  }
                });
              },
              done: function(model, resp){
                // FIXME: success is deprecated for done, but it must be chained
              },
            }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });
          
            return true;
          }
        }
      });
    
    }, //addVocabularyItem  
    

    download: function(url, resource, post_data){
      var self = this;
      
      if(url.search(/\?/) < 1 ){
        url = url + '?';
      }
      var self = this;
      var altCheckboxTemplate =  _.template('\
          <div class="form-group" style="margin-bottom: 0px;" > \
            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
              <label title="<%= help %>" for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>\
            </div>\
          </div>\
        ');
      var formSchema = {};
      
      formSchema['use_vocabularies'] = {
        title: 'Use vocabulary labels',
        help: 'If selected, vocabulary key values will be replaced with labels',
        key: 'use_vocabularies',
        type: 'Checkbox',
        template: altCheckboxTemplate
      };
      formSchema['use_titles'] = {
        title: 'Use column titles',
        help: 'If selected, column key values will be replaced with column titles',
        key: 'use_titles',
        type: 'Checkbox',
        template: altCheckboxTemplate
      };
      formSchema['raw_lists'] = {
        title: 'Export nested lists without list brackets',
        help: [ 'If selected, a nested list in a cell will be quoted, ',
                'but not denoted with brackets[]. '].join(''),
        key: 'raw_lists',
        type: 'Checkbox',
        template: altCheckboxTemplate
      };
      formSchema['content_type'] = {
        title: 'Download type',
        help: 'Select the data format',
        key: 'content_type',
        options: _.without(resource.content_types, 'json'), // never json
        type: 'Select',
        validators: ['required'],
        template: _.template([
          '<div class="input-group col-xs-6">',
          '   <label class="input-group-addon" for="<%= editorId %>" ',
          '         title="<%= help %>" ><%= title %></label>',
          '   <span  data-editor></span>',
          '</div>',
        ].join('')),
        editorClass: 'form-control'
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema
      });
      var formFields = new FormFields();
      formFields.set('use_vocabularies', true);
      formFields.set('use_titles', true);
      formFields.set('raw_lists', true);
      formFields.set('content_type', formSchema['content_type'].options[0]);
      var form = new Backbone.Form({
        model: formFields,
        template: _.template([
          "<div>",
          "<form data-fieldsets class='form-horizontal container' >",
          "</form>",
          // tmpFrame is a target for the download
          '<iframe name="tmpFrame" id="tmpFrame" width="1" height="1" ',
          ' style="visibility:hidden;position:absolute;display:none"></iframe>',
          "</div>"
          ].join(''))
      });
      
      form.listenTo(form, "change", function(e){
        console.log('change');
        var content_type = form.getValue('content_type');
        if(content_type == 'sdf'){
          form.$el.find('[name="use_vocabularies"]').prop('disabled', false);
          form.$el.find('[name="use_titles"]').prop('disabled', false);
          form.$el.find('[name="raw_lists"]').prop('disabled', true);
        }else{
          form.$el.find('[name="use_vocabularies"]').prop('disabled', false);
          form.$el.find('[name="use_titles"]').prop('disabled', false);
          form.$el.find('[name="raw_lists"]').prop('disabled', false);
        }
      });
      
      var el = form.render().el;
      
      var default_content = form.getValue('content_type');
      console.log('default_content: ' + default_content);
      if(default_content != 'csv' && default_content != 'xls'){
        $(el).find('[name="use_vocabularies"]').prop('disabled', true);
        $(el).find('[name="use_titles"]').prop('disabled', true);
        $(el).find('[name="raw_lists"]').prop('disabled', true);
      }
      
      self.showModal({
        view: el,
        title: 'Download',  
        ok: function(event){
          var intervalCheckTime = 1000; // 1s
          var maxIntervals = 3600;      // 3600s
          var limitForDownload = 0;
          
          var errors = form.commit({ validate: true }); // runs schema and model validation
          if(!_.isEmpty(errors) ){
            console.log('errors', errors);
            _.each(_.keys(errors), function(key){
              form.$el.find('[name="'+key +'"]').parents('.form-group,.input-group').addClass('has-error');
            });
            return false;
          }
          var values = form.getValue();

          url += '&format=' + values['content_type']

          if(values['use_vocabularies']){
            url += '&use_vocabularies=true';
          }
          if(values['use_titles']){
            url += '&use_titles=true';
          }
          if(values['raw_lists']){
            url += '&raw_lists=true';
          }
          
          // When downloading via AJAX, the "Content-Disposition: attachment" 
          // does not trigger a "load" (completed) event; this workaraound
          // will set a cookie "downloadID" from the server when the download happens.
          // How to trigger a download and notify JavaScript when finished:
          // send a downloadID to the server and wait for a response cookie to appear.
          // The code here was helpful:
          // http://www.bennadel.com/blog/2533-tracking-file-download-events-using-javascript-and-coldfusion.htm
          
          // When tracking the download, we're going to have
          // the server echo back a cookie that will be set
          // when the download Response has been received.
          var downloadID = ( new Date() ).getTime();
          // Add the "downloadID" parameter for the server
          // Server will set a cookie on the response to signal download complete
          url += "&downloadID=" + downloadID;

          if(post_data){
            console.log('post_data found for download', post_data);
            // create a form for posting
            var el = form.$el.find('#tmpFrame');
            
            var postform = $("<form />", {
              method: 'POST',
              action: url
            });

            _.each(_.keys(post_data), function(key){
              var hiddenField = $('<input/>',{
                type: 'hidden',
                name: key,
                value: encodeURI(post_data[key])
              });
              postform.append(hiddenField);
            });
            el.append(postform);
            console.log('submitting post form...' + url);
            postform.submit();
            
          }else{
            // simple GET request
            form.$el.find('#tmpFrame').attr('src', url);
          }
          
          
          $('#loading').fadeIn({duration:100});

          // The local cookie cache is defined in the browser
          // as one large string; we need to search for the
          // name-value pattern with the above ID.
          var cookiePattern = new RegExp( ( "downloadID=" + downloadID ), "i" );

          // Now, we need to start watching the local Cookies to
          // see when the download ID has been updated by the
          // response headers.
          var cookieTimer = setInterval( checkCookies, intervalCheckTime );

          var i = 0;
          function checkCookies() {
            if ( document.cookie.search( cookiePattern ) >= 0 ) {
              clearInterval( cookieTimer );
              $('#loading').fadeOut({duration:100});
              return(
                console.log( "Download complete!!" )
              );
            }else if(i >= maxIntervals){
              clearInterval( cookieTimer );
              window.alert('download abort after tries: ' + i);
              return(
                console.log( "Download abort!!" )
              );
            }
            console.log(
              "File still downloading...",
              new Date().getTime()
            );
            i++;
          }
        }
        
      });
      
    },
    
    /**
     * show a modal pop up with only an "ok" button, no "cancel" button.
     */
    showModalMessage: function(options){
      var modalDialog = this.showModal(options);
      modalDialog.$el.find('#modal-cancel').hide();
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
      var okText = (options && options.okText)? options.okText : 'Continue';
      var cancelText = (options && options.cancelText)? 
        options.cancelText : 'Cancel and return to page';
      console.log('options', options);
      var template = _.template(modalOkCancelTemplate)({ 
            body: options.body,
            title: options.title } );
      var modalDialog = new Backbone.View({
          el: template,
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
                  self.clearPagePending();
                  if(callbackOk(event)===false){
                    return;
                  }
                  $('#modal').modal('hide');
              }
          },
      });
      modalDialog.render();
      if(!_.isUndefined(options.view)){
        modalDialog.$el.find('.modal-body').append(options.view);
      }
      modalDialog.$el.find('#modal-cancel').html(cancelText);
      modalDialog.$el.find('#modal-ok').html(okText);
      $modal = $('#modal');
      $modal.empty();
      $modal.html(modalDialog.$el);
      $modal.modal({show:true, backdrop: 'static'});
      return modalDialog;
    },
    
    createSearchString: function(searchHash){
      return _.map(_.pairs(searchHash), 
        function(keyval) {
          return keyval.join('=');
        }).join(SEARCH_DELIMITER)
    }
  });

  var appState = new AppState();
  
  appState._form_template = _.template([
     "<div class='form-horizontal container' id='add_value_field' >",
     "<form data-fieldsets class='form-horizontal container' >",
     "</form>",
     "</div>"].join(''));      
  appState._field_template = _.template([
    '<div class="form-group" >',
    '    <label class="control-label " for="<%= editorId %>"><%= title %></label>',
    '    <div class="" >',
    '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
    '      <div data-error class="text-danger" ></div>',
    '      <div><%= help %></div>',
    '    </div>',
    '  </div>',
  ].join(''));
  
  
  
  appState.schemaClass = new SchemaClass(); // make accessible to outside world
  
  appState.resources = {};   // TO be retrieved from the server 
  
  appState.apiVersion = API_VERSION;
  appState.reportsApiUri = REPORTS_API_URI;
  appState.dbApiUri = DB_API_URI;
  appState.DEBUG = DEBUG;
  appState.LIST_ARGS = ['rpp','page','includes','order','log', 'children','search'];      
  appState.SEARCH_DELIMITER = SEARCH_DELIMITER;
  appState.HEADER_APILOG_COMMENT = 'X-APILOG-COMMENT';

  Iccbl.appModel = appState;
  
  return appState;
});