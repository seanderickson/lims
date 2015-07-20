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
    return this.filterKeys('detail');
  };  
  SchemaClass.prototype.filterKeys = function(visibility_term)
  {
    console.log('filter keys for: ' + visibility_term );
    var self = this;
    var keys = Iccbl.sortOnOrdinal(
      _.keys(self.fields), self.fields)
    var detailKeys = _(keys).filter(function(key){
      return _.has(self.fields, key) && 
          _.has(self.fields[key], 'visibility') && 
          _.contains(self.fields[key]['visibility'], visibility_term);
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
    
    /**
     * Generate a list of "options" suitable for use in a user multiselect.
     * [ { val: username, label: name:username }, ... ]
     */
    getUserOptions: function(callBack){
      this.getUsers(function(users){
        var options = [];
        users.each(function(user){
          var username = user.get('username');
          var name = user.get('name');
          options.unshift({ val: username, label: name + ':' + username });
        });
        callBack(options);
      });
    },
    
    getUsers: function(callBack) {
      var self = this;
      users = this.get('users');
      if(_.isEmpty(users)){
        console.log('get all users from the server...');
        var resourceUrl = self.dbApiUri + '/screensaveruser'
        Iccbl.getCollectionOnClient(resourceUrl, function(collection){
          self.set('users', collection);
          callBack(collection);
        });
      }
      else{
        callBack(users);
      }
    },
    
    getUserGroupOptions: function(callBack){
      this.getUserGroups(function(usergroups){
        var options = [];
        usergroups.each(function(usergroup){
          var name = usergroup.get('name');
          options.unshift({ val: name, label: name });
        });
        callBack(options);
      });
    },
    
    getUserGroups: function(callBack) {
      var self = this;
      usergroups = this.get('usergroups');
      if(_.isEmpty(usergroups)){
        console.log('get all UserGroups from the server...');
        var resourceUrl = self.reportsApiUri + '/usergroup'
        Iccbl.getCollectionOnClient(resourceUrl, function(collection){
          self.set('usergroups', collection);
          callBack(collection);
        });
      }
      else{
        callBack(usergroups);
      }
      
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
        callBack(vocabularies);
      });
    },    
    
    /**
     * Cache an "options" array for all permissions, for the editor UI
     */
    setPermissionsOptions: function(resources){
      var self = this;
      var permissionOptions = [];
      var optionGroups = { resources: []};
      _.each(_.keys(resources), function(rkey){
        var resource = resources[rkey];
        if(!_.has(resource,'schema')){
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
        var fields = resource['schema']['fields'];
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
     * return the set of vocabulary items for a scope:
     *    v.scope: { 
     *      v.key1 : { scope: , key: , title: , ordinal: },
     *      v.key2 : {},
     *      etc.,...
     *    },
     * @param scope - scope to search for
     * - "scope" can also be a regex, and will be matched to all scopes using
     * String.prototype.match(candidateScope,scope)
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
            if(candidateScope.match(scope)){
              console.log('matching: ' + scope + ', to: ' + candidateScope );
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
        defaults : {},
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
          // to force inclusion of all columns: data: { includes: '*' },
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
    
    /**
     * FIXME: a jquery ajax error function should handle the args:
     * Type: Function( jqXHR jqXHR, String textStatus, String errorThrown )
     */
    jqXHRerror: function(jqXHR, textStatus, errorThrown ) {
      this.jqXHRFail(jqXHR,textStatus,errorThrown);
//      //console.log('error fetching the model: '+ model + ', response:
//      // ' + JSON.stringify(response));
//      var url = options.url || model.url || '';
//      var msg = 'Error locating resource: ' + url;
//      this.error(msg);
//      
//      var sep = '\n';
//      if (!_.isUndefined(response.status))
//          msg += sep + response.status;
//      if (!_.isUndefined(response.statusText))
//          msg += sep + response.statusText;
//      if (!_.isEmpty(response.responseText))
//          msg += sep + response.responseText;
//      
//      if(DEBUG) window.alert(msg);
//      else console.log(msg);
//      // TODO: use Bootstrap inscreen alert classed message div
    },
    
    /**
     * TODO: a jquery ajax error function should handle the args:
     * Type: Function( jqXHR jqXHR, String textStatus, String errorThrown )
     */
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
        // FIXME: consider a model attribute on app_state for messages, as this
        // pattern is needed for additions to the array
        this.set({'messages': msgs},{silent: true} );
        this.trigger('change:messages');
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

    download: function(url, resource, post_data){
      var self = this;
      
      if(url.search(/\?/) < 1 ){
        url = url + '?';
      }
      var self = this;
      var altCheckboxTemplate =  _.template('\
          <div class="form-group" style="margin-bottom: 0px;" > \
            <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > \
              <label title="<%= help %>" for="<%= editorId %>"><span data-editor\><%= title %></label>\
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

      var form = new Backbone.Form({
        model: formFields,
        template: _.template([
                              "<div>",
          "<form data-fieldsets class='form-horizontal container' >",
          "</form>",
          // tmpFrame is a target for the download
          '<iframe name="tmpFrame" id="tmpFrame" width="1" height="1" style="visibility:hidden;position:absolute;display:none"></iframe>',
          "</div>"
          ].join(''))
      });
      
      form.listenTo(form, "change", function(e){
        console.log('change');
        var content_type = form.getValue('content_type');
        console.log('content_type: ' + content_type );
        if(content_type != 'csv' && content_type != 'xls'){
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
          
          form.commit();
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
  appState.DEBUG = DEBUG;
  appState.LIST_ARGS = ['rpp','page','includes','order','log', 'children','search'];      
  appState.SEARCH_DELIMITER = ';';
  appState.HEADER_APILOG_COMMENT = 'X-APILOG-COMMENT';

  Iccbl.appModel = appState;
  
  return appState;
});