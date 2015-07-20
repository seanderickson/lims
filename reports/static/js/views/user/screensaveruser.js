define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_layout',
    'views/list2',
    'views/user/user2',
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, ReportsUserView, layout) {

  var UserView = ReportsUserView.extend({

    screensaver_tabbed_resources: {
      userchecklistitem: {
        description: "User Checklist Items",
        title: "User Checklist Items",
        invoke: "setUserChecklistItems"
      }
    },
    
    initialize: function(args) {
      UserView.__super__.initialize.apply(this, arguments);      
      var self = this;
      this.tabbed_resources = _.extend({},
        this.tabbed_resources, this.screensaver_tabbed_resources);
      
      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail' && !appModel.hasPermission(
            self.tabbed_resources[key].resource,'read')){
          delete self.tabbed_resources[key];
        }
      });
//      
//      _.bindAll(this, 'click_tab');
    },
    
//    template: _.template(layout),
    
//    tabbed_resources: {
//        detail: { 
//          description: 'User Details', 
//          title: 'User Details', 
//          invoke: 'setDetail' },
//        usergroup: { 
//          description: 'User Groups', 
//          title: 'User Groups', 
//          invoke: 'setGroups', 
//          resource: 'usergroup' },
//        permission: { 
//          description: 'User Permissions', 
//          title: 'User Permissions', 
//          invoke: 'setPermissions',
//          resource: 'permission' },
//    },
    
//    events: {
//      'click ul.nav-tabs >li': 'click_tab',
//    },

//    /**
//     * Layoutmanager hook
//     */
//    serialize: function() {
//      return {
//        'tab_resources': this.tabbed_resources
//      }      
//    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      UserView.__super__.afterRender.apply(this, arguments);      
//      var viewId = 'detail';
//      if (!_.isEmpty(this.uriStack)){
//        viewId = this.uriStack.shift();
//        if(viewId == '+add'){
//          this.uriStack.unshift(viewId);
//          viewId = 'detail';
//        }
//        if (!_.has(this.tabbed_resources, viewId)){
//          var msg = 'could not find the tabbed resource: ' + viewId;
//          window.alert(msg);
//          throw msg;
//        }
//      }
//      this.change_to_tab(viewId);
    },
    
    
    click_tab : function(event){
      UserView.__super__.click_tab.apply(this, arguments);      
    },

    change_to_tab: function(key){
      UserView.__super__.change_to_tab.apply(this, arguments);      
    },
    
    setUserChecklistItems: function(delegateStack) {
      var self = this;
      var key = 'userchecklistitem';
      var resource = appModel.getResource('userchecklistitem');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'checklistitems'].join('/');

      var show_save_button = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="save_button" href="#">',
            'save</a>'
          ].join(''));
      var form_template = [
         "<form  class='form-horizontal container' >",
         "<div data-fields='comments'/>",
         "</form>"];
      var altFieldTemplate =  _.template('\
        <div class="form-group" > \
            <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
            <div class="col-sm-10" >\
              <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />\
              <div data-error class="text-danger" ></div>\
              <div><%= help %></div>\
            </div> \
          </div>\
        ');
      // Build the form model

      var FormFields = Backbone.Model.extend({
        schema: {
          comments: {
            title: 'Comments',
            key: 'comments',
            type: 'TextArea',
            validators: ['required'], 
            template: altFieldTemplate
          }
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: _.template(form_template.join(''))
      });
      var _form_el = form.render().el;
      
      var receiverFunction = function(collection){
        var Collection = Backbone.Collection.extend({
          url: url,
          toJSON: function(){
            return {
              objects: Collection.__super__.toJSON.apply(this) 
            };
          }
          
        });
        var changedCollection = new Collection();
        var MyModel = Backbone.Model.extend({
          url: url,
          initialize : function() {
            this.on('change', function(model, options) {
              // Prevent save on update
              if (options.save === false)
                  return;
              model.url = url;
              if(_.isEmpty(model.get('status_date'))){
                model.set('status_date', (new Date()).toISOString());
              }
              if(_.isEmpty(model.get('admin_username'))){
                model.set('admin_username', appModel.getCurrentUser().username);
              }
              changedCollection.add(model);
            });
          },
        });
        collection.model = MyModel;
        show_save_button.click(function(e){
          e.preventDefault();
          console.log('changed collection', changedCollection,changedCollection.url);
          
          if(changedCollection.isEmpty()){
            appModel.error('nothing changed');
            return;
          }
          
          appModel.showModal({
            okText: 'ok',
            ok: function(e){
              e.preventDefault();
              var errors = form.commit();
              if(!_.isEmpty(errors)){
                console.log('form errors, abort submit: ' + JSON.stringify(errors));
                return false;
              }else{
                Backbone.sync("patch",changedCollection,
                  {
                    error: function(){
                      appModel.jqXHRError.apply(this,arguments);
                      console.log('error, refetch', arguments);
                      collection.fetch();
                    },
                  }
                );
              }
            },
            view: _form_el,
            title: 'Save changes?'  
          })
        });
        
        view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: resource.schema,
          resource: resource,
          url: url,
          collection: collection,
          extraControls: [show_save_button]
        }});
        Backbone.Layout.setupView(view);
        self.consumedStack = [key]; 
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
      };
      
      
      Iccbl.getCollectionOnClient(url, receiverFunction);
    },    

  });

  return UserView;
});