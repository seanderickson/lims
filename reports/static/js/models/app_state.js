define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'templates/modal_ok_cancel.html'    
], function($, _, Backbone, Backgrid, Iccbl, modalOkCancelTemplate ){
  
  var API_VERSION = 'api/v1';
  var REPORTS_API_URI = '/reports/' + API_VERSION;
  var DB_API_URI = '/db/' + API_VERSION;
  
  var DEBUG = false;
  
  var API_RESULT_META = 'meta';
  var API_RESULT_DATA = 'objects';
  var API_MSG_RESULT = 'Result';
    
  var SchemaClass = Iccbl.SchemaClass = function() {};
  SchemaClass.prototype.detailKeys = function()
  {
    return this.filterKeys('visibility','d');
  };  
  SchemaClass.prototype.adminKeys = function()
  {
    var self = this;
    var keys = Iccbl.sortOnOrdinal(
      _.keys(self.fields), self.fields)
    var adminKeys = _(keys).filter(function(key){
      return _.result(self.fields[key], 'is_admin', false);
    });
    return adminKeys;
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
    
    var do_not_display_keys = this.filterKeys('visibility','none');
    
    keys = _.difference(keys, do_not_display_keys);
    
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
  SchemaClass.prototype.allKeys = function()
  {
    return this.filterKeys('visibility', '-none');
  };
  SchemaClass.prototype.filterKeys = function(select_field, visibility_term)
  {
    var self = this;
    var keys = Iccbl.sortOnOrdinal(
      _.keys(self.fields), self.fields)
    if (visibility_term.charAt(0) == '-'){
      visibility_term = visibility_term.slice(1);
      return _(keys).filter(function(key){
        return _.has(self.fields, key) && 
            _.has(self.fields[key], select_field) && 
            ! _.contains(self.fields[key][select_field], visibility_term);
      });
    }else{
      return _(keys).filter(function(key){
        return _.has(self.fields, key) && 
            _.has(self.fields[key], select_field) && 
            _.contains(self.fields[key][select_field], visibility_term);
      });
    }
  };  
  
  /**
   * Groups the keys by the display_options "group":
   * @ return an array of: 
   * [ {
   *     group: name,
   *     fields: [fieldkey1, fieldkey2, fieldkey3... ]
   *    }, ...]
   */
  SchemaClass.prototype.groupedKeys = function(keys){
    var self = this;
    var groupedKeys = [];
    var groups = {}; // note, hash is being used only as a lazy cache, not returned
    _.each(keys,function(key){
      var fi = self.fields[key];
      if (fi && fi.display_options && fi.display_options.group){
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

      path: '',
      actualStack: [],
      
      // NOTE: in backbone, change notifications are only detected for the 
      // object identity.  For this reason, each of the following state items is
      // kept at the top level.
      current_view: 'home',
      current_resource_id: 'home',
      current_options: {},
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
      var self = this;
      
      self.userProps = {};
      this.on('change:users', function(){
        _.each(_.keys(self.userProps), function(prop){
          self.unset(prop);
        });
        //self.unset('adminUsers');
        //self.unset('principal_investigators');
        //self.unset('userOptions');
        //self.unset('adminUserOptions');
      });
      this.on('change:screens', function(){
        self.unset('screenOptions');
      });
      this.on('change:libraries', function(){
        self.unset('libraryOptions');
        self.unset('screeningLibraryOptions-small_molecule');
        self.unset('screeningLibraryOptions-rnai');
        self.unset('platedLibraryOptions');
      });
      this.on('change:usergroups', function(){
        self.unset('userGroupOptions');
      });
      this.on('change:platelocations', function(){
        self.unset('locationHash')
      });
      _.bindAll(this,'error',
        'setCurrentUser','getResources','getVocabularies',
        'getAdminUserOptions','getUserOptions','getUserGroupOptions',
        'getPrincipalInvestigatorOptions','getLibraries','getLibraryOptions',
        'getScreens','getScreenOptions','getPlateLocationTree',
        'getPlateLocations');
    },
        
    start: function(callBack) {
      var self = this;
      console.log('start app_state.js');
      self.getAppDataFromServer(function(){
        self.getResources(function(){
          self.setCurrentUser(function(){
            self.getVocabularies(function(vocabularies){
              self.set('vocabularies', vocabularies);
              self.initialized();
              callBack();
            });
          });
        });
      });
    },
    
    initialized: function() {
      var self = this;
      var jobResource = self.getResource('job');
      var JobCollection = Iccbl.CollectionOnClient.extend({
        // explicitly define the id so that collection compare & equals work
        modelId: function(attrs) {
          return Iccbl.getIdFromIdAttribute( attrs, jobResource);
        },
        url: jobResource.apiUri
      })
      self.jobCollection = new JobCollection();
      
      var currentUser = self.getCurrentUser();
      if (currentUser.is_staff || currentUser.is_superuser ){
        // fetch pending jobs
        data_for_fetch = {
          limit:0,
          username: currentUser.username,
          'state__in': ['pending','submitted', 'processing'],
        };
        // background jobs
        self.jobCollection.fetch({
          data: data_for_fetch,
          // prevent triggering global event handler: ajaxStart loading gif
          global: false 
        }).done(function(data, textStatus, jqXHR){
          console.log('initial job fetch', arguments);
        }).fail(function(){ 
          Iccbl.appModel.jqXHRfail.apply(this,arguments); 
        });
      }
    },
    
    readCookie: function(name) {
      var nameEQ = escape(name) + "=";
      var ca = document.cookie.split(';');
      for (var i = 0; i < ca.length; i++) {
          var c = ca[i];
          while (c.charAt(0) === ' ') c = c.substring(1, c.length);
          if (c.indexOf(nameEQ) === 0) return unescape(c.substring(nameEQ.length, c.length));
      }
      return null;
    },
    
    /**
     * perform lazy fetch for Admin data
     */
    initializeAdminMode: function(callback) {
      var self = this;
      console.log('initializeAdminMode...');
      // Pre-fetch options for the search_box
      $(this).queue([
         self.getAdminUserOptions,
         self.getUserOptions,
         self.getPrincipalInvestigatorOptions,
         self.getUserGroupOptions,
         callback]);
    },
    
    setCurrentUser: function(callBack) {
      console.log('set the current user', window.user);
      var self = this;
      // NOTE: can be replaced with "user" from the reports app
      // 'screensaveruser' is needed for the database app
      this.getModel('screensaveruser', window.user, 
        function(model){
          // NOTE that Backbone.Model.toJSON returns a shallow copy as a object,
          // not as a JSON string.
          self.currentUser = model.toJSON();
          self.set('currentUser', model);
          console.log('Current user: ' + JSON.stringify(self.currentUser));
          
          if(callBack) callBack(); 
        }
      );
    },
    
    getCurrentUser: function(){
      return this.currentUser;
    },
    
    /**
     * Use the AppModel to hold search state:
     * - search_box writes
     * - list views read
     */
    setSearch: function(searchId, search_object){
      if (DEBUG) console.log('setSearch', searchId, search_object)
      localStorage.setItem(''+searchId, JSON.stringify(search_object));
      this.getSearch(searchId);
    },
    
    getSearch: function(searchId){
      var obj = localStorage.getItem(''+searchId);
      if(obj) obj = JSON.parse(obj);
      return obj;
    },

    /** 
     * Find the search data in the uriStack, if it exists: either in:
     * URI_PATH_COMPLEX_SEARCH or,
     * URI_PATH_ENCODED_SEARCH
     * 
     * @return {
     *    API_PARAM_SEARCH: search data - resource specific search data
     *    search_id: integer for local storage
     *    uriStack: uriStack with the search component removed
     *    errors: [] empty if none
     *  }
     */
    findComplexSearch: function(uriStack){
      var result = {};
      if (_.isEmpty(uriStack)) return null;
      var appModel = Iccbl.appModel;
      
      // Complex search: uses local storage for large sets of search data
      if (_.contains(uriStack,appModel.URI_PATH_COMPLEX_SEARCH) ){
        var index = _.indexOf(uriStack,appModel.URI_PATH_COMPLEX_SEARCH);
        var searchId = uriStack[index+1];
        
        if (_.isEmpty(searchId)){
          result['errors'] = 'no complex search id found:', uriStack;
          return result;
        }
        if (!_.isNaN(parseInt(searchId))) {
          var search_data = appModel.getSearch(searchId);
          if(_.isUndefined(search_data) || _.isEmpty(search_data)){
            result['errors'] = 'Browser state no longer contains searchId:'+searchId 
              +'", search must be performed again';
            return result;
          }else{
            result[appModel.API_PARAM_SEARCH_ID] = searchId;
            result[appModel.API_PARAM_SEARCH] = search_data;
          }  
        }else{
          result['errors'] = 'search ID is not a valid number: ' + searchId;
          return result;
        }
        var newStack = uriStack.slice(0,index);
        result['uriStack'] = newStack.concat(uriStack.slice(index+2));
      
      // Encoded search: search data are encoded as part of the URI
      } else if (_.contains(uriStack, appModel.URI_PATH_ENCODED_SEARCH)){
        var index = _.indexOf(uriStack,appModel.URI_PATH_ENCODED_SEARCH);
        result[appModel.API_PARAM_SEARCH] = uriStack[index+1];
        var newStack = uriStack.slice(0,index);
        result['uriStack'] = newStack.concat(uriStack.slice(index+2));
      } else {
        console.log('no search data found on the uriStack');
      }
      return result;
    },
    
    getLabAffiliationOptions: function(callBack){
      if (DEBUG) console.log('getLabAffiliationOptions...', callBack);
      var self = this;
      var prop = 'labAffiliationOptions';
      var options = this.get(prop);
      if(!options){
        this.getLabAffiliations(function(affiliations){
          options = [];
          affiliations.each(function(affiliation){
            var id = affiliation.get('lab_affiliation_id');
            var name = affiliation.get('name');
            var category = affiliation.get('category')
            options.push({ val: id, label: name + ' (' + category + ')' });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
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
            
            var label = short_name;
            if (short_name != library_name){
              label += ': ' + library_name;
            }
            label += ' (' + library.get('start_plate');
            label += '-' + library.get('end_plate') + ')';
            
            options.push({ val: short_name, label: label });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
    },

    getPlateLocationTree: function(callBack, options){
      var self = this;
      var options = options || {};
      var flush = options.flush || false;
      
      var prop = 'locationHash';
      var locationHash = this.get(prop);
      if(!locationHash || flush){
        this.getPlateLocations(function(locations){
          var locationHash = {};
          locations.each(function(location){
            var room = location.get('room');
            if (!room) return;
            var roomHash = _.result(locationHash, room,{});
            locationHash[room] = roomHash;
            var freezer = location.get('freezer');
            if (!freezer) return;
            var freezerHash = _.result(roomHash, freezer, {});
            roomHash[freezer] = freezerHash;
            var shelf = location.get('shelf');
            if (!shelf) return;
            var shelfHash = _.result(freezerHash,shelf,{});
            freezerHash[shelf]=shelfHash;
            var bin = location.get('bin');
            if (!bin) return;
            shelfHash[bin] = bin;
          });
          self.set(prop,locationHash);
          if (callBack) callBack(locationHash);
        }, options);
      }else{
        if (callBack) callBack(locationHash);
        else return locationHash;
      }
    },

    getScreenOptions: function(callBack){
      var self = this;
      var prop = 'screenOptions';
      var options = this.get(prop);
      if(!options){
        this.getScreens(function(screens){
          options = [];
          screens.each(function(screen){
            var facility_id = screen.get('facility_id');
            var title = screen.get('title');
            options.push({ val: facility_id, label: facility_id + ': ' + title });
          });
          self.set(prop,options);
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
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
        else return options;
      }
    },

    getPlatedLibraryOptions: function(callBack){
      /** return libraries with plate copies **/
      var self = this;
      var prop = 'platedLibraryOptions';
      var options = this.get(prop);
      if(!options){
        this.getLibraries(function(libraries){
          options = [];
          libraries.each(function(library){
            if (! library.has('copies')){
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
        else return options;
      }
    },

    /** 
		 * Get the libraries covered for search:
		 * TODO: Test 20180322 
		 **/
    getLibrariesForSearch: function(searchArray, errors){
      var libraries = this.getLibraries();
      
      if (!libraries){
        console.log('libraries are not initialized');
        return;
      }
      
      function libraryForPlate(plate_number){
        return libraries.find(function(library){
          return ( library.get('start_plate') <= plate_number 
              && library.get('end_plate') >= plate_number);
        });
      }; 
      
      function librariesForRange(start,end){
        return libraries.filter(function(library){
          return ( library.get('start_plate') <= end
              && library.get('end_plate') >= start);
        });
      }
      
      var librariesFound = [];
      
      _.each(searchArray, function(parsedSearch){
        if (_.has(parsedSearch, 'well_ids')){
          _.each(parsedSearch['well_ids'], function(well_id){
            var plate_number = parseInt(well_id.split(':')[0]);
            var library = libraryForPlate(plate_number);
            if (!library){
              errors.push('No library found for plate: ' + plate_number);
            } else {
              librariesFound.unshift(library);
            }
          });
        }
        if (_.has(parsedSearch, 'plate_ranges')){
          _.each(parsedSearch['plate_ranges'], function(plate_range){
            var start_end = plate_range.split('-');
            var libraries = librariesForRange(
              parseInt(start_end[0]), parseInt(start_end[1]));
            if (_.isEmpty(libraries)){
              errors.push('No libraries found for range: ' + plate_range);
            }
            // consider checking that all libraries are same type
            librariesFound = librariesFound.concat(libraries);
          });
        }
        if (_.has(parsedSearch, 'plates')){
          _.each(parsedSearch['plates'], function(plate){
            var library = libraryForPlate(plate);
            if (!library){
              errors.push('No libraries found for plate number: ' + plate);
            }else{
              librariesFound.unshift(library);
            }
          });
        }
      });
      
      var librariesFound = _.unique(librariesFound);
      
      // Check that all libraries are of the same type
      var screen_types = {};
      _.each(librariesFound, function(library){
        var screen_type = library.get('screen_type');
        var list = _.result(screen_types, screen_type, []);
        list.push(library.get('short_name'));
        screen_types[screen_type] = list;
      });
      
      if(_.size(screen_types)>1){
        errors.push('Well searches may not contain different library types: ' 
          + _.keys(screen_types).join(', '));
      }
      
      return librariesFound;
    },
    
    
    /**
     * Generate a list of "options" suitable for use in a user multiselect.
     * [ { val: username, label: name:username }, ... ]
     */
    getAdminUserOptions: function(callBack, resource){
      var self = this;
      var prop = 'adminUserOptions';
      // Note: check for the "empty" object passed when the caller is JQuery.queue function call
      if (!_.isString(resource)){
        resource = null;
      }
      
      if (!_.isEmpty(resource)){
        prop += '-' + resource;
      }
      var options = this.get(prop);
      if(!options){
        this.getAdminUsers(function(users){
          var options = [];
          users.each(function(user){
            var username = user.get('username');
            var name = user.get('name');
            options.push({ val: username, label: name + ': ' + username });
          }, resource);
          self.userProps[prop] = new Date();
          self.set(prop,options);
          if (callBack) callBack(options);
        }, resource);
      }else{
        if (callBack) callBack(options);
        else return options;
      }
    },
    
    getUserOptions: function(callBack){
      var self = this;
      var is_super = this.getCurrentUser().is_superuser;
      var prop = 'userOptions';
      var options = this.get(prop);
      if(!options){
        self.getUsers(function(users){
          var options = [];
          users.each(function(user){
            // Construct the option label:
            // NOTE: for search box searching, must match the patttern:
            // appState.USER_OPTION_PATTERN
            // /([^\[]+)(\[(\w+)\])?(\:\s+(\d+))?/
            // where (1=name) (3=ecommons) (4=ssID)
            var user_id = user.get('screensaver_user_id');
            var username = user.get('username');
            var ecommons = user.get('ecommons')
            var name = user.get('name');
            var _label = name;
            if(!_.isEmpty(username)){
              _label += ' [' + username + ']';
            }else if (!_.isEmpty(ecommons)){
              _label += ' [' + ecommons + ']';
            }
            if (is_super){
              _label += ': ' + user_id;
            }
            options.unshift({ val: user_id, label: _label});
          });
          self.set(prop,options);
          self.userProps[prop] = new Date();
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
    },
    
    getUsernameOptions: function(callBack){
      var self = this;
      var prop = 'usernameOptions';
      var options = this.get(prop);
      if(!options){
        self.getUsers(function(users){
          var options = [];
          users.each(function(user){
            var user_id = user.get('screensaver_user_id');
            var username = user.get('username');
            if (_.isEmpty(username)){
              return;
            }
            var name = user.get('name');
            var _label = name;
            if(!_.isEmpty(username)) _label += ': ' + username;
            _label += '(' + user_id + ')';
            options.unshift({ val: username, label: _label });
          });
          self.set(prop,options);
          self.userProps[prop] = new Date();
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
    },

    getPrincipalInvestigatorOptions: function(callBack){
      var self = this;
      var prop = 'piOptions';
      var options = this.get(prop);
      if(!options){
        this.getPrincipalInvestigators(function(users){
          var options = [{ val:'',label:'' }];
          users.each(function(user){
            var user_id = user.get('screensaver_user_id');
            var lab_name = user.get('lab_name');
            options.push({ val: user_id, label: lab_name });
          });
          self.set(prop,options);
          self.userProps[prop] = new Date();
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
    },
    
    getUserIdsInGroupOptions: function(usergroup, callBack){
      return this.getUserInGroupOptions(
        usergroup, 'screensaver_user_id', '{name} ({screensaver_user_id})', callBack );
    },
    
    getUsernamesInGroupOptions: function(usergroup, callBack){
      return this.getUserInGroupOptions(
        usergroup, 'username', '{name} ({username})', callBack );
    },
    
    getUserInGroupOptions: function(usergroup, val_prop, label_prop, callBack){
      var self = this;
      var prop = ['usergroup', usergroup, val_prop, label_prop, 'options'].join('_');
      var options = this.get(prop);
      if(!options){
        this.getUsersInGroup(usergroup, function(users){
          var options = [{ val:'',label:'' }];
          users.each(function(user){
            options.push({ 
              val: user.get(val_prop), 
              label: Iccbl.formatString(label_prop, user) });
          });
          self.set(prop,options);
          self.userProps[prop] = new Date();
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
    },
    
    getUsersInGroup: function(usergroup, callBack) {
      var prop = usergroup + '_users';
      this.userProps[prop] = new Date();
      return this.getCachedResourceCollection(
        prop, 
        this.dbApiUri + '/screensaveruser', 
        { 
          usergroups__eq: usergroup,
          exact_fields: ['username','name','email'], 
          order_by: ['name']
        }, 
        callBack );
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
          self.userProps[prop] = new Date();
          if (callBack) callBack(options);
        });
      }else{
        if (callBack) callBack(options);
        else return options;
      }
//      return options;
    },

    getPermissionsOptions: function(){
      return this.get('permissionOptions');
    },
    
    getLabAffiliations: function(callback){
      data_for_get = { 
        exact_fields: [
           'lab_affiliation_id','name', 'category'],
        use_vocabularies: true,
        order_by: ['category','name']
      };
      return this.getCachedResourceCollection(
        'labAffiliations', this.dbApiUri + '/labaffiliation', data_for_get, callback );
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
    
    getPlateLocations: function(callback, options){
      data_for_get = { 
        exact_fields: ['room','freezer','shelf','bin']
      };
      return this.getCachedResourceCollection(
        'platelocations', this.dbApiUri + '/platelocation', data_for_get, 
        callback, options );
    },
    
    getScreens: function(callback){
      data_for_get = { 
        order_by: ['facility_id'],
        exact_fields: ['title','facility_id','data_sharing_level','screen_type'], 
        study_type__is_null: true
      };
      return this.getCachedResourceCollection(
        'screens', this.dbApiUri + '/screen', data_for_get, callback );
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
      var prop = 'principal_investigators';
      this.userProps[prop] = new Date();
      return this.getCachedResourceCollection(
        prop, 
        this.dbApiUri + '/screensaveruser', 
        { 
          lab_affiliation_id__is_blank: false,
          classification__eq: this.VOCAB_USER_CLASSIFICATION_PI,
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
          exact_fields: ['username','name','email',
            'sm_data_sharing_level','rnai_data_sharing_level'], 
          order_by: ['name']
        }, 
        callback );
    },
    
    /**
     * @param resource (optional) to filter for users having write permission 
     * on the resource
     */
    getAdminUsers: function(callback, resource) {
      var self = this;
      function filterCallback(users){
        if (!_.isEmpty(resource)){
          var filteredUsers = new Backbone.Collection(
            users.filter(function(user){
              return self._permissionMatch(
                resource, 'write', user.get('all_permissions'));
            })
          );
          callback(filteredUsers);
        } else {
          callback(users);
        }
      };
      var prop = 'adminUsers';
      self.userProps[prop] = new Date();
      return this.getCachedResourceCollection(
        prop, 
        this.dbApiUri + '/screensaveruser', 
        { 
          is_staff__eq: true,
          exact_fields: ['username','name','email', 'all_permissions'], 
          order_by: ['name'], 
          
        }, 
        filterCallback );
    },

    getUserGroups: function(callback) {
      var prop = 'usergroups';
      this.userProps[prop] = new Date();
      return this.getCachedResourceCollection(
        prop, 
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
      
      var resourceUrl = self.reportsApiUri + '/vocabulary'
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
     * [ 
     *   group: "user",
     *   options: [
     *     { val: 'key', label: 'label' }, { val: ... }
     *   ]
     * ]
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
     * Retrun an array of options for a Iccbl.SelectCell
     */
    getVocabularySelectCellArray: function(scope){
      var options = [];
      try{
        var vocabulary = Iccbl.appModel.getVocabulary(scope);
        _.each(_.keys(vocabulary),function(choice){
          options.push([vocabulary[choice].title,choice]);
        });
      }catch(e){
        var msg = 'Vocabulary unavailable: vocabulary_scope_ref: ' + scope;
        console.log(msg,e);
        this.error(msg);
      }
      return options;
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
            if (DEBUG) {
              console.log(
                'skipping retired vocab: ',choice,vocabulary[choice].title );
            }
          }else{
            choiceHash.push({ 
              val: choice, 
              label: vocabulary[choice].title, 
              ordinal: vocabulary[choice].ordinal });
          }
        });
      }catch(e){
        var msg = 'Vocabulary unavailable: vocabulary_scope_ref: ' + scope;
        console.log(msg,e);
        this.error(msg);
      }
      choiceHash = _.sortBy(choiceHash, function(choice){
        return choice.ordinal;
      });
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
            if (DEBUG){
                console.log(
                'matching: ' + '^' + scope + '$' + ', to: ' + candidateScope );
            }
            _.extend(matchedVocabularies,vocabularies[candidateScope]);
          }
        });
        if(!_.isEmpty(matchedVocabularies)){
          if (DEBUG){
            console.log('matchedVocabularies', scope, matchedVocabularies );
          }
          return matchedVocabularies;
        }
        throw "Unknown vocabulary: " + scope;
      }
      if (DEBUG){
        console.log('vocabulary: ' + scope, vocabularies[scope]);
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
            { vocabulary: scope,
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
      var resourceUrl = self.dbApiUri + '/resource'
      
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
        // 2. Extend the (api)resources with the uiResource value overrides
        var ui_resources = self.get('ui_resources');
        var resources = {};
        _.each(resourcesCollection, function(resource){
          resources[resource.key] = resource;
          _.each(_.values(resource.fields), self.parseSchemaField );
          
          if (_.has(ui_resources, resource.key)) {
            ui_resources[resource.key] = _.extend(
                {}, resource, ui_resources[resource.key]);
          } else {
            ui_resources[resource.key] = _.extend({},resource);
          }
        });

        // 3. for "virtual" uiResources, find the underlying apiResource, 
        // and extend it using the uiResource overrides
        _.each(_.keys(ui_resources), function(key){
          var ui_resource = ui_resources[key];
          // then augment the ui resources with their api_resource, if different
          if ( key !== ui_resource.api_resource ) {
            ui_resources[key] = _.extend(
                {}, resources[ui_resource.api_resource], ui_resource);
          }
        });
        
        var app_data = self.getAppData();
        ui_resources['home']['title'] = app_data.get('app_name');
        
        // set up permissions for all the resources
        self.setPermissionsOptions(self.get('ui_resources'));
        
        
        
        if(callBack) callBack();                
      }, null, function(jqXHR, textStatus, errorThrown){
        // pass a special error handler, because getResources runs before 
        // layout is rendered
        var msg = self.parseJqXHRfail.apply(this, arguments);
        if (jqXHR && _.has(jqXHR,'responseJSON') && !_.isEmpty(jqXHR.responseJSON) ) {
          msg += '\n\n' + self.print_json(jqXHR.responseJSON);
        }
        window.alert(msg);
      });
      
      console.log('finished getResources')
    },

    createNewModel: function(resourceId, defaults) {

      if (DEBUG) console.log('create new model for: ',resourceId, defaults );
      var self = this;
      var resource = self.getResource(resourceId);
      return this.newModelFromResource(resource,defaults);
    },
    
    get_field_defaults: function(fields){

      var defaults = {};
      _.each(_.keys(fields), function(key){
        var field = fields[key];
        if( (!_.isUndefined(field.default) && !_.isNull(field.default))
            && ( _.isNumber(field.default) || _.isBoolean(field.default)
                || !_.isEmpty(field.default)))
        {
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
            if (field.data_type == 'date' && field.default == 'now'){
              defaults[key] = Iccbl.getISODateString(new Date());
            }else{
              self.error('Warning, unparseable default for field: ' 
                + field.key + ', value: ' + field.default );
              defaults[key] = null;
            }
          }
        }else{
          defaults[key] = null; // placeholders for the edit form
        }
      });
      
      return defaults;
    },

    newModelFromResource: function(resource, defaults) {
      
      var finalDefaults = _.extend({},
        this.get_field_defaults(resource.fields), defaults);
      delete finalDefaults['resource_uri']
      delete finalDefaults['id'];
      var NewModel = Backbone.Model.extend({
        __classname: 'iccbl-model-' + resource.key,
        urlRoot: resource.apiUri , 
        defaults: finalDefaults,
        resource: resource,
        
        /**
         * Override parse:
         * - unnest the post_obj_response from { objects: [ post_obj_response ] }
         * - set the model.id; so that isNew() will report false
         */
        parse: function(resp,options){
          var serverAttrs = _.result(resp, Iccbl.appModel.API_RESULT_DATA);
          if (serverAttrs&& _.isArray(serverAttrs)){
            if (serverAttrs.length == 1){
              serverAttrs = serverAttrs[0];
            }
          }
          serverAttrs[this.idAttribute] = Iccbl.getIdFromIdAttribute(serverAttrs,resource );
          return serverAttrs;
        }
      });
      var newModel = new NewModel();
      newModel.resource = resource;
      
      return newModel;
      
    },
    
    newModelFromResourceOld: function(resource,defaults){

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
            if( (!_.isUndefined(field.default) && !_.isNull(field.default))
                && ( _.isNumber(field.default) || _.isBoolean(field.default)
                    || !_.isEmpty(field.default)))
            {
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
                  defaults[key] = null;
                }
              }
            }else{
              defaults[key] = null; // placeholders for the edit form
            }
          }
        }
      });

      var NewModel = Backbone.Model.extend({
        __classname: 'iccbl-model-' + resource.key,
        urlRoot: resource.apiUri , 
        defaults: defaults,
        resource: resource,
        /**
         * Override parse:
         * - unnest the post_obj_response from { objects: [ post_obj_response ] }
         * - set the model.id; so that isNew() will report false
         */
        parse: function(resp,options){
          var serverAttrs = _.result(resp, Iccbl.appModel.API_RESULT_DATA);
          if (serverAttrs&& _.isArray(serverAttrs)){
            if (serverAttrs.length == 1){
              serverAttrs = serverAttrs[0];
            }
          }
          serverAttrs[this.idAttribute] = Iccbl.getIdFromIdAttribute(serverAttrs,resource );
          return serverAttrs;
        }
      });
      var newModel = new NewModel();
      newModel.resource = resource;
      
      return newModel;
    },

    _get_user_screen_choices: function(userModel) {
        var screenOptions = this.getScreenOptions();
        var userScreens = _.union(
          userModel.get('screens_collaborator'),
          userModel.get('screens_lab_head'),
          userModel.get('screens_lead'));
        var options = _.filter(screenOptions, function(option){
          return _.contains(userScreens, option['val']);
        });
        if (DEBUG) {
          console.log('found user screens: ' + 
            userModel.get('screensaver_user_id'), options);
        }
        return options;
    },
    
    _get_screen_member_choices: function(screenModel) {
      
      var userOptions = this.getUserOptions();
      var screenMembers = [
        screenModel.get('lead_screener_id'),
        screenModel.get('lab_head_id')
      ];
      var collaborator_ids = screenModel.get('collaborator_ids');
      if (!_.isEmpty(collaborator_ids)){
        screenMembers = _.union(screenMembers, collaborator_ids);
      }
      screenMembers = _.map(screenMembers, function(id){
        return '' + id;
      });
      var options = _.filter(userOptions, function(option){
        return _.contains(screenMembers, '' + option['val']);
      });
      // reverse, default sorting
      options.reverse();
      if (DEBUG) {
        console.log('found user screens: ' + 
          userModel.get('screensaver_user_id'), options);
      }
      return options;
    },
    
    /**
     * Get a model from the server
     */
    getModel: function(resourceId, key, callBack, options) {
      var resource = this.getResource(resourceId);
      return this.getModelFromResource(resource, key, callBack, options);
    },
    
    getModelFromResource: function(resource, key, callBack, options) {
      var self = this;
      var options = options || {};
      var failCallback = options.failCallback;
      var data_for_get = _.extend({ includes: '*' }, options.data_for_get );
      if(_.isArray(key)){
        key = key.join('/');
      }
      var url = options.url || resource.apiUri + '/' + key;
      if (DEBUG) console.log('fetch model', url);
      var ModelClass = Backbone.Model.extend({
        url : url,
        defaults : { },
        parse: function(resp, options){
          return _.result(resp, Iccbl.appModel.API_RESULT_DATA, resp);
        },
        fetch: function(options) {
          var options = options || {};
          if (!_.has(options,'data')){
            options['data'] = data_for_get;
          }
          return Backbone.Model.prototype.fetch.apply(this, [options]);
        }        
      });
      var instance = new ModelClass();
      instance.fetch({
        success : function(model) {
          model.resource = resource;
          model.options = options;
          model.key = key;
          model.set('id', key); // used by model.isNew(), if set use PATCH on save
          try{
            callBack(model);
          } catch (e) {
            console.log('uncaught error: ', e);
            self.error('error displaying model: ' + model + ', '+ e);
          }
        }
      }).fail(function(){ 
        console.log('instance fetch fail:', arguments);
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
      
      if(!currentUser.is_staff && !currentUser.is_superuser ){
        menu = this.get('screener_menu');
      }
      
      else if(!currentUser.is_superuser) { 
        
        // Use permissions to show only allowed menus for the user
        var new_submenus = {}
        _.each(_.keys(menu.submenus), function(key){
          // 'reports' and 'admin' are visible for the admins
          if ( (key == 'reports' || key == 'admin')
            && _.contains(currentUser.usergroups, 'readEverythingAdmin')){

            new_submenus[key] = menu.submenus[key];
          }else{
            if(self.hasPermission(key,'read')){
              new_submenus[key] = menu.submenus[key];
            }else{
              console.log('user: ' + currentUser.username 
                  + ', doesnt have permission to view the menu: ' + key );
            }
          }
        });
        menu.submenus = new_submenus;
      }
      return menu;
    },
        
    isEditable: function(resourceId){
      try {
        return _.contains(_.result(this.getResource(resourceId),'visibility'), 'e');
      } catch (e) {
        return false;
      }
    },
    
    /**
     * Test if the current user has the resource/permission - 
     * - if permission is unset, 
     * will check if the user has *any* permission on the resource.
     */
    hasPermission: function(resource, permission){
      if (DEBUG) console.log('hasPermission', resource, permission);
      var self = this;
      if (self.getCurrentUser().is_superuser) return true;
      
      return self._permissionMatch(
        resource, permission, self.getCurrentUser().all_permissions);
    
    },
    
    _permissionMatch: function(resource, permission, permissions_to_search) {
      var r_perm = 'resource/'+ resource;
      if( permission=='edit'){
        permission = 'write';
      }
      if(!_.isUndefined(permission) && permission == 'write' ){
        // only check for 'write' permission; write implies read
        r_perm += '/'+ permission;
      }// otherwise, will return true if user has either permission
      var match = _.find(
          permissions_to_search, 
          function(p){
            if(p.indexOf(r_perm) > -1 ) {
              return true;
            }
          });
      return !_.isUndefined(match);
    },
    
    hasGroup: function(groupName) {
      var self = this;
      if(self.getCurrentUser().is_superuser) return true;
      return _.contains(self.getCurrentUser().usergroups, groupName);
    },
    
    /**
     * Retrieve a resource (ResourceResource) from the local cache:
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
          ],
          fields: {
            field1: {}
          }
        }
     *    
     */
    getResource: function(resourceId){
      var self = this;
      var uiResources = this.get('ui_resources');
      if(!_.has(uiResources, resourceId)) {
        throw "Unknown resource: " + resourceId;
      }
      var resource = this.cloneResource(uiResources[resourceId]);
      
      if (resourceId == 'my_screens'){
        var suResource = uiResources['screensaveruser'];
        var currentUser = self.getCurrentUser();
        resource.apiUri = [suResource.apiUri, currentUser.screensaver_user_id, 'screens'].join('/');
      }
      
      return resource;
    },
    
    /**
     * Deep copy of a Resource.
     */
    cloneResource: function(resource){
      resource = _.clone(resource);
      resource.fields = _.clone(resource.fields);
      _.each(_.keys(resource.fields), function(key){
        var field = resource.fields[key];
        resource.fields[key] = _.clone(field);
        _.each(_.keys(field), function(fk){
          field[fk] = _.clone(field[fk]);
        });
      });
      return resource;
    },

    getAppData: function(){
      return this.get('app_data');
    },
    
    getAppDataFromServer: function(callback){
      var self = this;
      var ModelClass = Backbone.Model.extend({
        url : [REPORTS_API_URI,'resource','app_data'].join('/')
      });
      var instance = new ModelClass();
      instance.fetch({
        success: function(model){
          self.set('app_data', model);
          callback(model);
        }
      });
    },
    
    getResourceFromUrl: function(schemaUrl, callback, options){
      var self = this;
      var options = options || {};
      var ModelClass = Backbone.Model.extend({
        url : schemaUrl,
        defaults : {},
        parse: function(resp, options){
          return _.result(resp, Iccbl.appModel.API_RESULT_DATA, resp);
        }
      });
      var instance = new ModelClass();
      instance.fetch({
        data: options,
        success : function(model) {
          schema = model.toJSON();
          _.each(_.values(schema.fields), self.parseSchemaField );
          var schemaClass = new SchemaClass();
          var ui_resource = _.extend(schema, schemaClass);
          ui_resource.apiUri = '/' + [
            ui_resource.api_name,self.apiVersion,ui_resource.key].join('/');
          callback(ui_resource)
        }
      }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
    },

    parseJqXHRfail: function(jqXHR, textStatus, errorThrown){
      console.log('jqXHRfail', textStatus, errorThrown);
      var msg = textStatus || 'Error';
      msg = msg.charAt(0).toUpperCase() + msg.slice(1);
      if (errorThrown){
        msg += ': ' + errorThrown; 
      }
      if (this.url){
        // * url will be set for most backbone objects doing the fetching
        if(_.isFunction(this.url)) msg += ':\n ' + this.url();
        else msg += ':\n ' + this.url;
      }
      if (jqXHR){
        if (jqXHR.status){
          msg += ':\n status:' + jqXHR.status;
        }
      }
      return msg;
    },
    
    /**
     * Process a jQuery.jqXHR.fail callback for the ajax() call.
     */
    jqXHRfail: function(jqXHR, textStatus, errorThrown){
      var msg = Iccbl.appModel.parseJqXHRfail.apply(this, arguments);
      $(document).trigger('ajaxComplete');
      if (jqXHR && _.has(jqXHR,'responseJSON') && !_.isEmpty(jqXHR.responseJSON) ) {
        Iccbl.appModel.showJsonMessages(jqXHR.responseJSON);
      } else {
        // Iccbl.appModel.error(msg);
        Iccbl.appModel.showModalMessage({
          buttons_on_top: false,
          body: msg,
          title: 'Error'
        });
      }
    },
    /**
     * Creates a hierarchical printout of an object using 
     * JSON.stringify
     */
    print_json: function(obj){
      var str;
      function replacer(key, val){
        // Convert single list values into strings
        if (_.isArray(val) && val.length==1){
          return val[0];
        }else{
          return val;
        }
      };
      try
      {
        str = JSON.stringify(obj, replacer, 2);
        
        // remove escapes around quotes
        str = str.replace(/\\(["'])/g, '$1');
        // remove python unicode literal
        str = str.replace(/u(['"])/g,'$1');
        // remove single quotes
        str = str.replace(/"/g,'');
        // remove object brackets
        str = str.replace(/[{}]+/g,'');
        // remove list brackets
        str = str.replace(/[\[\]]+/g, '');
        // remove single comma on a line
        str = str.replace(/^\s*,\s*$/gm,'')
        // remove empty lines
        str = str.replace(/^\s*$\n+/gm,'');
        // convert escaped lines
        str = str.replace(/\\n+/gm,'<br/>');
      }catch(e){
        console.log('print_json', e, obj);
        str = '' + obj;
      }
        
      return str;
    },
    
    /**
     * Process an error dict into single string for display to the end user.
     */
    print_dict: function(dict, sep){
      var sep = sep || '<br/>';
      var self = this;
      var output = '';
      if (_.isObject(dict) && !_.isArray(dict)){
        output = _.map(_.keys(dict), function(key){
          var err = key + ': ';
          var parts = self.print_dict(dict[key], sep);
          if (parts.length > 40){
            err += sep + parts;
          } else {
            err += parts;
          }
          return err;
        }).join(sep);
        
      } else if (_.isArray(dict)){
        var cumulative = '';
        dict = _.map(dict, function(val){
          if (_.isObject(val) && !_.isArray(val)){
            var parts = self.print_dict(val, sep);
            if (parts.length > 40){
              cumulative += sep + parts;
              val = sep + parts;
            } else {
              val = parts;
              cumulative += parts;
            }
          } else {
            cumulative += val;
            if (cumulative.length > 40 && cumulative != val){
              val = sep + val;
              cumulative = val;
            }
          }
          return val;
        });
        output += dict.join(', ');
      } else {
        output += dict;
      }
      return output;
    },
    
    /**
     * Parse connection result (sent as a JSON object)
     * 
     * @param options - usual options for showModalMessage, including:
     * - title
     * - ok function
     * - body - text to display if the data is not recognized
     */
    showConnectionResult: function(data, options){
      if (_.isObject(data) && !_.isString(data)){
        data = _.result(data,API_MSG_RESULT,data);
        data = _.result(data,API_RESULT_META,data);
        
        this.showJsonMessages(data, options);
        //var bodyMsg = this.print_dict(data);
        //this.showModalMessage(
        //  _.extend({}, options, {
        //    body: bodyMsg,
        //    buttons_on_top: true,
        //    buttons_on_bottom: false
        //  }));
      }else{
        console.log('Warn: data should have been parsed as json]', data);
        this.showModalMessage(
          _.extend({}, options, {
            okText: 'ok',
            body: 'Action Complete'
          }));
      }
    },
    
    /**
     * Show a JSON object in a modal dialog:
     * - transform the object into a table using a depth-first traversal:
     * - each row contains the nodes traversed to each leaf node.
     */
    showJsonMessages: function(jsonObj, options){
      
      if (DEBUG) console.log('showJsonMessages: ', jsonObj);
      var options = _.extend({
          buttons_on_top: false,
        }, options);
      
      if (!_.has(options,'title')){
        var title = "Messages";
        if(_.keys(jsonObj).length == 1){
          title = _.keys(jsonObj)[0];
          jsonObj = jsonObj[title];
          title = title.charAt(0).toUpperCase() + title.slice(1);
        }
        options['title'] = title;
      }
      var bodyMsg = this.print_json(jsonObj);
      var rowCount = (bodyMsg.match(new RegExp('\n', "g")) || []).length;
      if (rowCount > Iccbl.appModel.MAX_ROWS_IN_DIALOG_MSG){
        body = $('<textarea class="input-full" rows=' 
          + Iccbl.appModel.MAX_ROWS_IN_DIALOG_MSG + ' ></textarea>');
        body.val(bodyMsg);
        options['buttons_on_top'] = true;
        options['view'] = body;
      } else {
        bodyMsg = bodyMsg.replace(/\n/g, '<br/>');
        bodyMsg = bodyMsg.replace(/[\t ]/g,'&nbsp;');
        options['body'] = bodyMsg;
      }
      Iccbl.appModel.showModalMessage(options);
    },
    
    
    /**
     * Uses stock JSON.stringify to generate display information in a modal dialog.
     */
    showJsonDirect: function(jsonObj, options){
      if (DEBUG) console.log('showJsonDirect: ', jsonObj);
      var options = _.extend({
          buttons_on_top: false,
        }, options);
      if (!_.has(options,'title')){
        var title = "Messages";
        if(_.keys(jsonObj).length == 1){
          title = _.keys(jsonObj)[0];
          jsonObj = jsonObj[title];
          title = title.charAt(0).toUpperCase() + title.slice(1);
        }
        options['title'] = title;
      }
      var str = JSON.stringify(jsonObj, null, 2); // spacing level = 2
      var body = $('<textarea class="input-full" rows=' 
        + Iccbl.appModel.MAX_ROWS_IN_DIALOG_MSG + ' ></textarea>');
      body.val(str);
      options['buttons_on_top'] = true;
      options['view'] = body;
      Iccbl.appModel.showModalMessage(options);
    },
   
    clearErrors: function(){
      this.unset('messages');
    },
    
    error: function(msg){
      console.log('error: ', msg);
      var msgs = this.get('messages');
      if (msgs && msgs.length > 0) {
        msgs.push(msg);
        
        if(msgs.length > 5){
          msgs = msgs.splice(1, msgs.length-1);
        }
        // FIXME: consider a model attribute on app_state for messages, as this
        // pattern is needed for additions to the array
        this.set({'messages': msgs},{silent: true} );
        this.trigger('change:messages');
      }else{
        this.set('messages', [msg]);
      }
    }, 
    
    addBackgroundJob: function(jobData){
      console.log('set job', jobData);
      var jobResource = this.getResource('job');
      var url = jobResource.apiUri + '/' + jobData['id'];
      var ModelClass = Backbone.Model.extend({
        url : url,
        defaults : { },
        parse: function(resp, options){
          return _.result(resp, Iccbl.appModel.API_RESULT_DATA, resp);
        }
      });
      var job = new ModelClass(jobData);
      this.jobCollection.add(job);
      return job;
    },
    
    getJobCollection: function(){
      return this.jobCollection;
    },
    

    /**
     * setUriStack:
     * Initiates a "change:uriStack" event signal to view/controllers in the 
     * application:
     * Called by: UI events (menu click, tab actions)
     * Listeners: content.js, menu.js, search_box.js
     * Coding convention; the final designated target view (which may modify 
     * the stack as needed) will call "reportUriStack", initiating a cascade of 
     * "reportUriStack" calls to set the final URL and history.
     */  
    setUriStack: function(value, options){

      var options = _.extend( {'deferred_navigation': true }, options);
      
      // 20180425 - fix: window location automatically encodes the hash url,
      // and backbone router generates a second route change with the 
      // encoded url; ensure that uriStack contains only strings so that all
      // iterations match.
      value = _.map(value, function(fragment){
        return '' + fragment;
      });
      if(_.isEqual(this.get('uriStack'),value )){
        // signal a change event if this method was called with the same value
        // - for the menu signaling
          this.trigger('change:uriStack', this);
      }else{
        this.set({ uriStack: value }, options);
      }
    },
    
    /**
     * reportUriStack:
     * Reports the final URI stack for history management:
     * - router listens for "change:reportUriStack events.
     * 
     * UriStack parsing is performed by the hierarchy of controllers beginning 
     * at "content.js" and ending with a designated final view/controller. For
     * any given UriStack; there will be one final view/controller that is 
     * responsible for calling "reportUriStack" which triggers a cascade of
     * "reportUriStack" calls ending here.
     */
    reportUriStack: function(reportedUriStack, options) {
      // 20180425 - fix: window location automatically encodes the hash url,
      // and backbone router generates a second route change with the 
      // encoded url; ensure that uriStack contains only strings so that all
      // iterations match.
      reportedUriStack = _.map(reportedUriStack, function(fragment){
        return '' + fragment;
      });
      console.log('reportUriStack:', reportedUriStack, options);
      this.set({ 'uriStack': reportedUriStack }, { silent: true } ); 
      this.set({ 'reportedUriStack': reportedUriStack }, options);     
    },

    /**
     * Set flag to signal that the current page has pending changes;
     * (see setUriStack; modal confirm dialog will be triggered).
     */
    setPagePending: function(callback, message){
      var val = callback || true;
      this.set({'pagePending': val});
      if (!_.isUndefined(message)){
        this.set({'pagePendingMessage': message });
      }
    },
    clearPagePending: function(){
      this.unset('pagePending');
      this.unset('pagePendingMessage');
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
        options.title = options.title || 'Cancel pending changes?';
        options.cancelText = 'Return to save changes';
        options.okText = 'Cancel changes'
        var pendingMessage = self.get('pagePendingMessage') 
          || "Your changes are unsaved. Continue anyway?";
        options.body = options.body || pendingMessage;
        var pendingFunction = this.get('pagePending');
        if(_.isFunction(pendingFunction)){
          options.cancel = pendingFunction;
        }
        self.showModal(options);
      }

      // Clear out error messages after navigating away from page
      //      self.unset('messages');
    },

    /**
     * Show "comment" type ApiLogs in a separate table; ApiLogs with comments only
     **/
    createCommentTable: function(model,search_data, $target_el) {
      var self = this;
      
      if (!self.hasGroup('readEverythingAdmin')){
        return;
      }

      var search_data = _.extend({
        comment__is_blank: false
      }, search_data);
      
      var apilogResource = Iccbl.appModel.getResource('apilog');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: apilogResource.apiUri 
      });
      $target_el.empty();
      
      function build_table(collection) {
        if (collection.isEmpty()) {
          return;
        }

        // _.map(_.pairs(search_data...
        var search_ui_url = 'search/' + _.map(
          _.zip(_.keys(search_data),_.values(search_data)), 
          function(kv){
            return kv.join('=');
          }).join(Iccbl.appModel.SEARCH_DELIMITER);        
        var search_link = $('<a>',{
            text: 'Comments (Admins only)',
            target: '_blank',
            title: '',
            href: '#apilog/order/-date_time/' + search_ui_url
        });
        var addCommentButton = $([
          '<a class="btn btn-default btn-xs pull-down pull-right" ',
          'role="button" id="addCommentButton" href="#">',
          'Add</a>'
        ].join(''));
        
        addCommentButton.click(function(e){
          e.preventDefault();
          function processComment(formValues) {
            var url = [model.resource.apiUri,model.key].join('/');
            var headers = {};
            headers[Iccbl.appModel.HEADER_APILOG_COMMENT] = formValues['comments'];
            $.ajax({
              url: url,     
              cache: false,
              contentType: 'application/json', 
              dataType: 'json', // what is expected back from the server
              type: 'PATCH',
              headers: headers
            }).done(function(data, textStatus, jqXHR){
              Iccbl.appModel.createCommentTable(model,search_data, $target_el);
            }).fail(function(jqXHR, textStatus, errorThrown){
              Iccbl.appModel.jqXHRfail.apply(this,arguments); 
            });
            
          };
          Iccbl.appModel.showOkCommentForm({
            ok: processComment,
            okText: 'Save Comment',
            title: 'Add a comment to ' + model.resource.title + ': '
              + Iccbl.getTitleFromTitleAttribute(model,model.resource) + '?' 
          });
        });
        
        $target_el.append($([
          '<div class="col-xs-12" id="comment_table_title"></div>',
          '<div class="col-xs-12" id="comment_table_grid"/>'].join('')));
        
        $target_el.find('#comment_table_title').html(search_link);
        
        if (Iccbl.appModel.hasPermission(model.resource.key, 'write')){
          $target_el.find('#comment_table_title').append(addCommentButton);
        }
        
        collection.each(function(model) {
          model.set('date', Iccbl.getDateString(model.get('date_time')));
        });
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell.extend({
            className: 'admin-field'
          })
        };
        var columns = [
          _.extend({},colTemplate,{
            'name' : 'username',
            'label' : 'User',
            'description' : 'User performing the action',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.LinkCell.extend({
              hrefTemplate: '#screensaveruser/{username}',
              target: '_blank'
            })
          }),
          _.extend({},colTemplate,{
            'name' : 'date',
            'label' : 'Date',
            'description' : 'Date',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.LinkCell.extend({
              hrefTemplate: '#apilog/{log_uri}',
              target: '_blank'
            })
          }),
          _.extend({},colTemplate,{
            'name' : 'comment',
            'label' : 'Comment',
            'description' : 'Comment',
            'order': 1,
            'sortable': true,
            'cell': Iccbl.StringCell
            })
        ];
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();
        var grid = new Backgrid.Grid({
          columns: colModel,
          collection: collection,
          className: 'backgrid table-striped table-condensed table-hover '
        });
        
        $target_el.find('#comment_table_grid').html(grid.render().$el);
        
      }
      
      if (model.has('comment_data')) {
        build_table(new CollectionClass(model.get('comment_data')));
      } else {
        var comment_collection = new CollectionClass();
        
        comment_collection.fetch({
          data: _.extend({ 
              limit: 0,
              order_by: ['-date_time'],
              includes: ['log_uri']
            }, search_data),
          success: build_table
        }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
      }

    },
    
    showOkCommentForm: function(options){
      var self = this;
      var defaultOptions = {
        ok: function(){},
        okText: 'Ok',
        cancelText: 'Cancel and return to page',
        title: 'Save changes?',
        requireComment: false,
        commentValue: null
      };
      var options = _.extend({}, defaultOptions, options);
      var schema = {
        comments: {
          title: 'Comment (optional)',
          key: 'comments',
          type: 'TextArea',
          editorClass: 'input-full',
          template: self._field_template //altFieldTemplate
        }
      };
      
      if (options.requireComment == true) {
        schema['comments']['validators'] =['required'];
        schema['comments']['title'] = 'Comment';
      }
      var FormFields = Backbone.Model.extend({
        schema: schema
      });
      
      var defaultValue = {};
      if (options.commentValue){
        defaultValue['comments'] = options.commentValue;
      }
      var formFields = new FormFields(defaultValue);
      var form = new Backbone.Form({
        model: formFields,
        template: self._form_template //_.template(form_template.join(''))
      });
      var _form_el = form.render().el;
      
      if (!_.isEmpty(options.body)){
        _form_el.prepend(options.body);
      }
      options.view = _form_el;
      var okCallback = options.ok;
      options.ok = function(e){
        e.preventDefault();
        self.clearPagePending();
        var errors = form.commit();
        if(!_.isEmpty(errors)){
          console.log('form errors, abort submit: ' + JSON.stringify(errors));
          return false;
        }else{
          okCallback(form.getValue());            
        }
      };
      self.showModal(options);
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
        editorClass: 'col-xs-12',
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
            var resource = self.getResource('vocabulary');
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
    
    
    /**
     * Download the URL:
     * - use GET for simple URL (with no post_data)
     * - if post_data is provided, use POST with the post_data in a form.
     * @param $el document element to be used for a temporary iframe target for
     * the download.
     */
    downloadUrl: function(url, post_data) {
      
      // When downloading via AJAX, the "Content-Disposition: attachment" 
      // does not trigger a "load" (completed) event; 
      // to monitor the state of the download, a "downloadID" is sent to the 
      // server; the server responds with the downloadID in a response 
      // cookie when the download is complete.
      // The code here was helpful:
      // http://www.bennadel.com/blog/2533-tracking-file-download-events-using-javascript-and-coldfusion.htm

      // Set a downloadID that will be echoed back by the server as a cookie 
      // when the download is complete.
      var downloadID = ( new Date() ).getTime();
      var intervalCheckTime = 1000; // 1s
      var maxIntervals = 3600;      // 3600s
      var limitForDownload = 0;
      
      // Add the "downloadID" parameter for the server
      // Server will set a cookie on the response to signal download complete
      
      if (url.indexOf('?')>0){
        url += "&downloadID=" + downloadID;
      }else{
        url += "?downloadID=" + downloadID;
      }
      
      // tmpFrame is a target for the download
      var $iframe = $('#tmpFrame');
      $('#loading').fadeIn({duration:100});

      // Watch local document Cookies for the the download ID 
      // to be updated by the response headers.
      var i = 0;
      var cookiePattern = new RegExp( ( "downloadID=" + downloadID ), "i" );
      var cookieTimer = window.setInterval( checkCookies, intervalCheckTime );
      function checkCookies() {
        if ( document.cookie.search( cookiePattern ) >= 0 ) {
          window.clearInterval( cookieTimer );
          $('#loading').fadeOut({duration:100});
          return(
            console.log( "Download complete!!" )
          );
        }else if(i >= maxIntervals){
          window.clearInterval( cookieTimer );
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
      
      if(post_data){
        if (DEBUG) console.log('post_data found for download', post_data);
        console.log('create a form for posting...');
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
        // Add the Django "csrfmiddlewaretoken" to the post form
        postform.append($('<input/>', {
          type: 'hidden',
          name: "csrfmiddlewaretoken",
          value: this.readCookie('csrftoken')
        }));
        $iframe.append(postform);
        console.log('submitting post form...' + url);
        postform.submit();
        
      }else{
        // simple GET request
        $iframe.attr('src', url);
      }
    },

    /**
     * Display a download dialog for the url and resource.
     */
    download: function(url, resource, post_data){
      var self = this;
      
      var formSchema = {};
      
      var downloadOptions = self.getVocabularySelectOptions('resource.content_type');
      downloadOptions = _.filter(downloadOptions, function(option){
        return (
            option.val != 'json' // never json
            && _.contains(resource.content_types, option.val));
      });
      
      formSchema['use_vocabularies'] = {
        title: 'Use vocabulary labels',
        help: 'If selected, vocabulary key values will be replaced with labels',
        key: 'use_vocabularies',
        type: 'Checkbox',
        template: self._alt_checkbox_template
      };
      formSchema['use_titles'] = {
        title: 'Use column titles',
        help: 'If selected, column key values will be replaced with column titles',
        key: 'use_titles',
        type: 'Checkbox',
        template: self._alt_checkbox_template
      };
      formSchema['raw_lists'] = {
        title: 'Export nested lists without list brackets',
        help: [ 'If selected, a nested list in a cell will be quoted, ',
                'but not denoted with brackets[]. '].join(''),
        key: 'raw_lists',
        type: 'Checkbox',
        template: self._alt_checkbox_template
      };
      formSchema['data_interchange'] = {
        title: 'Download for data interchange',
        help: [ 'Use a format suitable for data upload ',
                '(vocabulary and column key values are used)'].join(''),
        key: 'data_interchange',
        type: 'Checkbox',
        template: self._alt_checkbox_template
      };
      formSchema['content_type'] = {
        title: 'Download type',
        help: 'Select the data format',
        key: 'content_type',
        options: downloadOptions,
//        options: _.without(resource.content_types, 'json'), // never json
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
      if (_.contains(resource.content_types, 'xls')){
        formFields.set('content_type','xls');
      } else {
        formFields.set('content_type', formSchema['content_type'].options[0]);
      }
      var form = new Backbone.Form({
        model: formFields,
        template: _.template([
          "<form data-fieldsets class='form-horizontal container' >",
          "</form>",
          ].join(''))
      });
      
      form.listenTo(form, "content_type:change", function(e){
        // NOTE: 20180910 - default option for sdf is data_interchange
        var content_type = form.getValue('content_type');
        if(content_type == 'sdf'){
          //form.$el.find('[name="use_vocabularies"]').prop('disabled', false);
          //form.$el.find('[name="use_titles"]').prop('disabled', false);
          //form.$el.find('[name="raw_lists"]').prop('disabled', false);
          form.setValue('data_interchange', true);
        }
        
      });
      
      form.listenTo(form, "change", function(e){
        
        // NOTE: 20180910 - allow options for sdf as well
        //var content_type = form.getValue('content_type');
        //if(content_type == 'sdf'){
        //  form.$el.find('[name="use_vocabularies"]').prop('disabled', false);
        //  form.$el.find('[name="use_titles"]').prop('disabled', false);
        //  form.$el.find('[name="raw_lists"]').prop('disabled', false);
        //  form.setValue('data_interchange', true);
        //}
        //else{
        //  form.$el.find('[name="use_vocabularies"]').prop('disabled', false);
        //  form.$el.find('[name="use_titles"]').prop('disabled', false);
        //  form.$el.find('[name="raw_lists"]').prop('disabled', false);
        //  form.$el.find('[name="data_interchange"]').prop('disabled', false);
        //}
        if (form.getValue('data_interchange') === true ){
          form.$el.find('[name="use_vocabularies"]').prop('disabled', true);
          form.$el.find('[name="use_titles"]').prop('disabled', true);
          form.$el.find('[name="raw_lists"]').prop('disabled', true);
        }else{
          form.$el.find('[name="use_vocabularies"]').prop('disabled', false);
          form.$el.find('[name="use_titles"]').prop('disabled', false);
          form.$el.find('[name="raw_lists"]').prop('disabled', false);
        }
      });
      var el = form.render().el;
      
      self.showModal({
        view: el,
        title: 'Download',  
        ok: function(event){
          
          var errors = form.commit({ validate: true }); // runs schema and model validation
          if(!_.isEmpty(errors) ){
            _.each(_.keys(errors), function(key){
              form.$el.find('[name="'+key +'"]')
                .parents('.form-group,.input-group').addClass('has-error');
            });
            return false;
          }
          var values = form.getValue();

          if(url.search(/\?/) < 1 ){
            url += '?';
          } else {
            url += '&';
          }
          url += 'format=' + values['content_type']

          if(values['use_vocabularies']){
            url += '&use_vocabularies=true';
          }
          if(values['use_titles']){
            url += '&use_titles=true';
          }
          if(values['raw_lists']){
            url += '&raw_lists=true';
          }
          if(values['data_interchange']){
            url += '&data_interchange=true';
          }
          
          self.downloadUrl(url, post_data);
          
        }
        
      });
      
    },
    
    showModalError: function(message){
      return this.showModal({
        ok_only: true,
        title: 'Error',
        body: message
      });
    },
    
    /**
     * show a modal pop up with only an "ok" button, no "cancel" button.
     */
    showModalMessage: function(options){
      options.ok_only = true;
      return this.showModal(options);
    },
    
    /**
     * Display a modal dialog
     */
    showModal: function(options){
      
      var self = this;
      var defaultOptions = {
        ok: function(){},
        cancel: function(){},
        okText: 'Continue',
        cancelText: 'Cancel and return to page',
        body: 'Body text',
        title: 'Title Text',
        buttons_on_top: false,
        buttons_on_bottom: true
        // Optional:
        // view: HTML or Jquery object to append to the modal body
        // width: override the modal-dialog width
        // css: (dict) set the .modal-dialog css
        // css_modal_content: (dict) set the .modal-content css
        // css_modal_body: (dict) set the .modal-body css
        // ok_only: boolean show only the "ok" button, hide the "cancel" button
        // 
      };
      
      var options = _.extend({}, defaultOptions, options);
      if (DEBUG) console.log('showModal options', options);
      
      if (!_.isUndefined(options.view)){
        delete options.body;
      }
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
                  options.cancel(event);
              },
              'click #modal-ok':function(event) {
                  console.log('ok button click event, '); 
                  event.preventDefault();
                  self.clearPagePending();
                  if(options.ok(event)===false){
                    return;
                  }
                  $('#modal').modal('hide');
              }
          },
      });
      modalDialog.render(); //.$el.css('display','block');
      if(!_.isUndefined(options.view)){
        modalDialog.$el.find('.modal-body').append(options.view);
      }
      modalDialog.$el.find('.modal-body').css('padding-top','0px');
      if (options.ok_only == true){
        modalDialog.$el.find('#modal-cancel').hide();
      } else {
        modalDialog.$el.find('#modal-cancel').html(options.cancelText);
      }
      modalDialog.$el.find('#modal-ok').html(options.okText);
      $modal = $('#modal');
      $modal.empty();
      $modal.html(modalDialog.$el);
//      modalDialog.$el.show();
      $modal.one('show.bs.modal', function () {
        // $('.modal-content').css('height',$( window ).height()*0.95);
        $('.modal-content').css('height', 'auto');
        $('.modal-content').css('max-height','100%');
      });   
      if (options.width){
        $('.modal-dialog').css('width', options.width);
      }
      if (options.css){
        $('.modal-dialog').css(options.css);
      }
      if (options.css_modal_content){
        $('.modal-content').css(options.css_modal_content);
      }
      if (options.css_modal_body){
        $('.modal-body').css(options.css_modal_body);
      }
      
      if (options.buttons_on_top == true){
        $('#top_modal_buttons').show();
      } else {
        $('#top_modal_buttons').hide();
      }
      if (options.buttons_on_bottom == true){
        $('#modal_footer').show();
      } else {
        $('#modal_footer').hide();
      }
      
      $.fn.drags = function(opt) {
        // Makes the target "handle" element draggable
        // See: https://css-tricks.com/snippets/jquery/draggable-without-jquery-ui/
        opt = $.extend({handle:"",cursor:"move"}, opt);

        if(opt.handle === "") {
            var $el = this;
        } else {
            var $el = this.find(opt.handle);
        }

        return $el.css('cursor', opt.cursor).on("mousedown", function(e) {
            if(opt.handle === "") {
                var $drag = $(this).addClass('draggable');
            } else {
                var $drag = $(this).addClass('active-handle').parent().addClass('draggable');
            }
            var z_idx = $drag.css('z-index'),
                drg_h = $drag.outerHeight(),
                drg_w = $drag.outerWidth(),
                pos_y = $drag.offset().top + drg_h - e.pageY,
                pos_x = $drag.offset().left + drg_w - e.pageX;
            $drag.css('z-index', 1000).parents().on("mousemove", function(e) {
                $('.draggable').offset({
                    top:e.pageY + pos_y - drg_h,
                    left:e.pageX + pos_x - drg_w
                }).on("mouseup", function() {
                    $(this).removeClass('draggable').css('z-index', z_idx);
                });
            });
            e.preventDefault(); // disable selection
        }).on("mouseup", function() {
            if(opt.handle === "") {
                $(this).removeClass('draggable');
            } else {
                $(this).removeClass('active-handle').parent().removeClass('draggable');
            }
        });

      };
      
      $modal.one('shown.bs.modal', function () {
        $('#modal').find('.form').find('input').first().focus();
        $("#modal").drags({
            handle: ".modal-header"
        });
      });
      $modal.modal({
        backdrop: 'static',
        keyboard: false, // prevent the escape key from closing the modal
        show: true
      });
      
      //$('#modal').one('hidden.bs.modal', '.modal', function () {
      //  console.log('hidden.bs.modal', arguments);
      //});
      return modalDialog;
    },
    
    /**
     * Create a search string for the URIstack
     * NOTE: items on the URIStack avoid looking like URL params; therefore
     * the "&" separator is not used.
     * @see SEARCH_DELIMITER
     */
    createSearchString: function(searchHash){
      return _.map(_.pairs(searchHash), 
        function(keyval) {
          return keyval.join('=');
        }).join(Iccbl.appModel.SEARCH_DELIMITER)
    }
  });

  var appState = new AppState();
  
  appState._form_template = _.template([
     '<div class="form-horizontal container" id="_form_template" >',
     '<form data-fieldsets class="form-horizontal container" autocomplete="off">',
     "</form>",
     '<div id="data-error" class="has-error" ></div>',
     "</div>"].join(''));      
  appState._field_template = _.template([
    '<div class="form-group" key="form-group-<%=key%>" >',
    '    <label class="control-label " for="<%= editorId %>" title="<%= help %>" ><%= title %></label>',
    '    <div class="" >',
    '      <div data-editor key="<%=key%>" style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
    '      <div data-error class="text-danger" ></div>',
    '    </div>',
    '  </div>',
  ].join(''));
  appState._alt_checkbox_template =  _.template('\
    <div class="form-group" style="margin-bottom: 0px;" > \
      <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
        <label title="<%= help %>" for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>\
      </div>\
    </div>\
  ');
  
  appState._horizontal_form_field_template = _.template([
    '<div class="form-group" >',
    '  <label class="control-label col-sm-6" for="<%= editorId %> "title="<%= help %>" ><%= title %></label>',
    '  <div class="col-sm-6" >',
    '    <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
    '    <div data-error class="text-danger" ></div>',
    //   <div><%= help %></div>',
    '  </div>',
    '</div>',
  ].join(''));
  appState._horizontal_form_2col_field_template = _.template([
    '<div class="form-group" >',
    '  <label class="control-label col-sm-2" for="<%= editorId %> "title="<%= help %>" ><%= title %></label>',
    '  <div class="col-sm-10" >',
    '    <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
    '    <div data-error class="text-danger" ></div>',
    //   <div><%= help %></div>',
    '  </div>',
    '</div>',
  ].join(''));
  
  appState.schemaClass = new SchemaClass(); // make accessible to outside world
  appState.resources = {};   // TO be retrieved from the server 
  
  appState.apiVersion = API_VERSION;
  appState.reportsApiUri = REPORTS_API_URI;
  appState.dbApiUri = DB_API_URI;
  appState.DEBUG = DEBUG;

  appState.HEADER_APILOG_COMMENT = 'X-APILOG-COMMENT';

  // Constants shared with the API
  // FIXME: Get shared constants from shared API_CONSTANTS file
  appState.API_RESULT_META = API_RESULT_META;
  appState.API_RESULT_DATA = API_RESULT_DATA;
  appState.API_MSG_RESULT = API_MSG_RESULT;
  appState.API_META_MSG_WARNING = 'Warning';

  appState.API_PARAM_SHOW_OTHER_REAGENTS = 'show_other_reagents';
  appState.API_PARAM_SHOW_ALTERNATE_SELECTIONS = 'show_alternate_selections';

  appState.API_PARAM_DC_IDS = 'dc_ids';
  
  appState.API_PARAM_SHOW_RESTRICTED = 'show_restricted'
  appState.API_PARAM_SHOW_ARCHIVED = 'show_archived'

  appState.API_PARAM_VOLUME_OVERRIDE = 'volume_override';
  appState.API_PARAM_SET_DESELECTED_TO_ZERO = 'set_deselected_to_zero';
  appState.API_PARAM_OVERRIDE = 'override';
  appState.API_PARAM_PATCH_PREVIEW_MODE = 'patch_with_preview';
  
  appState.MSG_SEARCH_SIZE_EXCEEDED = 
    'Maximum allowed search terms: {size_limit}' + 
    ', number of terms entered: {actual_size}';
  appState.API_MSG_LCPS_INSUFFICIENT_VOLUME = 'Insufficient volume';
  appState.VOCAB_USER_CLASSIFICATION_PI = 'principal_investigator';
  
  appState.USER_OPTION_PATTERN = /([^:\[]+)(\[(\w+)\])?(\:\s+(\d+))?/;
  //appState.USER_OPTION_PATTERN = /([A-Za-z\,\-\(\) ]+)(\[(\w+)\])?(\:\s+(\d+))?/;
  
  /**
   * URIStack search element SEARCH_DELIMITER
   * NOTE: items on the URIStack avoid looking like URL params; therefore
   * the "&" separator is not used.
   */
  appState.SEARCH_DELIMITER = ';';
  /** 
   * TODO:raw search line encode serves the same purpose as the SEARCH_DELIMITER
   * - line encoder was created for compound name search; may include semicolon?
   */
  appState.UI_PARAM_RAW_SEARCH_LINE_ENCODE = '|';

  appState.MAX_SEARCH_ARRAY_SIZE = 100;
  appState.MAX_RAW_SEARCHES_IN_URL = 10;
  appState.MAX_ROWS_IN_DIALOG_MSG = 20;
  
  // Use MAX_PRECISION to remove floating point math errors
  appState.MAX_PRECISION = 12;
  
  // API Param elements are used by the API
  appState.API_PARAM_SEARCH_ID = 'search_id'
  appState.API_PARAM_NESTED_SEARCH = 'nested_search_data';
  appState.API_PARAM_SEARCH = 'search';
  appState.API_PARAM_ENCODED_SEARCH = 'esearch';
  
  // URI path elements are used to create the URI stack
  appState.URI_PATH_COMPLEX_SEARCH = 'csearch';
  appState.URI_PATH_ENCODED_SEARCH = 'esearch';
  appState.URI_PATH_SEARCH = 'search';
  appState.LIST_ARGS = [
    appState.URI_PATH_ENCODED_SEARCH, appState.URI_PATH_COMPLEX_SEARCH,
    'rpp','page','includes','order','log', 'children',
    appState.URI_PATH_SEARCH];      
  
  appState.LIST_DELIMITER_SUB_ARRAY = '$';

  Iccbl.appModel = appState;
  
  return appState;
});