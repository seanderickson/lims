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
    'views/generic_edit',
    'text!templates/generic-tabbed.html',
    'bootstrap-datepicker'
], function($, _, Backbone, Iccbl, layoutmanager, 
            appModel, DetailLayout, 
            ListView, ReportsUserView, EditView, layout) {

  var UserView = ReportsUserView.extend({

    screensaver_tabbed_resources: {
      userchecklistitem: {
        description: "User Checklist Items",
        title: "User Checklist Items",
        invoke: "setUserChecklistItems",
        resource: 'userchecklistitem'
      },
      attachedfile: {
        description: "Attached Files",
        title: "Attached Files",
        invoke: "setAttachedFiles",
        resource: 'attachedfile'
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
    },
    
    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      UserView.__super__.afterRender.apply(this, arguments);      
    },
    
    click_tab : function(event){
      UserView.__super__.click_tab.apply(this, arguments);      
    },

    change_to_tab: function(key){
      UserView.__super__.change_to_tab.apply(this, arguments);      
    },
    
    setAttachedFiles: function(delegateStack) {
      var self = this;
      var key = 'attachedfile';
      var resource = appModel.getResource('attachedfile');
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'attachedfiles'].join('/');
      var upload_attached_file = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="save_button" href="#">',
          'Add</a>'
        ].join(''));
      var show_delete = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="show_delete" href="#">',
            'Delete</a>'
          ].join(''));
      
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource.schema,
        resource: resource,
        url: url,
        extraControls: [upload_attached_file, show_delete]
      }});
      upload_attached_file.click(function(e){
        e.preventDefault();
        self.upload(view.collection)
      });
      show_delete.click(function(e){
        e.preventDefault();
        if (! view.grid.columns.findWhere({name: 'deletor'})){
          view.grid.columns.unshift({ 
            name: 'deletor', label: 'Delete', text:'X', 
            description: 'delete record', 
            cell: Iccbl.DeleteCell, sortable: false });
        }
      });

      Backbone.Layout.setupView(view);
      self.consumedStack = [key]; 
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      
    },

    upload: function(attachedfileCollection){
      var self = this;
      var form_template = [
         "<div class='form-horizontal container' id='upload_attached_file_form' >",
         "<form data-fieldsets class='form-horizontal container' >",
         "<div class='form-group' ><input type='file' name='fileInput' /></div>",
         "</form>",
         "</div>"].join('');      
      var choiceHash = {}
      try{
        var vocabulary = Iccbl.appModel.getVocabulary('attachedfiletype.user');
          _.each(_.keys(vocabulary),function(choice){
            choiceHash[choice] = vocabulary[choice].title;
          });
      }catch(e){
        console.log('on get vocabulary', e);
        self.appModel.error('Error locating vocabulary: ' + 'attachedfiletype.user');
      }
      
      var fieldTemplate = _.template([
        '<div class="form-group" >',
        '    <label class="control-label " for="<%= editorId %>"><%= title %></label>',
        '    <div class="" >',
        '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
        '      <div data-error class="text-danger" ></div>',
        '      <div><%= help %></div>',
        '    </div>',
        '  </div>',
      ].join(''));
      var fileDateTemplate = _.template([
        '<div class="form-group" >',
        '    <label class="control-label" for="<%= editorId %>"><%= title %></label>',
        '    <div class="" >',
        '      <div data-editor style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
        '      <div data-error class="text-danger" ></div>',
        '      <div><%= help %></div>',
        '    </div>',
        '</div>',
      ].join(''));
      
      var formSchema = {};
      formSchema['file_type'] = {
        title: 'File Type',
        key: 'file_type',
        type: 'Select',
        options: choiceHash,
        template: fieldTemplate
      };
      formSchema['file_date'] = {
        title: 'File Date',
        key: 'file_date',
        type: EditView.DatePicker,
        template: fileDateTemplate
      };
      formSchema['file_name'] = {
        title: 'Option 2: Name',
        key: 'file_name',
        type: 'TextArea',
        template: fieldTemplate
      };
      formSchema['contents'] = {
        title: 'Option 2: Contents',
        key: 'contents',
        type: 'TextArea',
        template: fieldTemplate
      };
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        type: 'TextArea',
        template: fieldTemplate
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          console.log('form validate', attrs);
          var errs = {};
          var file = $('input[name="fileInput"]')[0].files[0]; 
          if (file) {
            if (!_.isEmpty(attrs.contents)){
              console.log('error, multiple file uploads specified');
              errs.contents = 'Specify either file or contents, not both';
            }
          } else {
            console.log('2',attrs.contents, _.isEmpty(attrs.contents));
            if (_.isEmpty(attrs.contents)){
              errs.contents = 'Specify either file or contents';
            }else{
              if (_.isEmpty(attrs.file_name)){
                errs.file_name = 'Must specify a filename with the file contents';
              }
            }
          }
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: _.template(form_template)
      });
      var _form_el = form.render().el;

      var dialog = appModel.showModal({
          okText: 'upload',
          ok: function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }); // runs schema and model validation
            if(!_.isEmpty(errors) ){
              console.log('form errors, abort submit: ',errors);
              _.each(_.keys(errors), function(key){
                $('[name="'+key +'"').parents('.form-group').addClass('has-error');
              });
              return false;
            }else{
              var values = form.getValue();
              console.log('form values', values);
              var comments = values['comments'];
              var headers = {};
              headers[appModel.HEADER_APILOG_COMMENT] = comments;
              
              var data = new FormData();
              _.each(_.keys(values), function(key){
                data.append(key,values[key])
              });

              var file = $('input[name="fileInput"]')[0].files[0];
              if(file){
                data.append('attached_file',file);
              }
              console.log('file: ' , file);
              
              var url = [self.model.resource.apiUri, 
                         self.model.key,
                         'attachedfiles'].join('/');
              $.ajax({
                url: url,    
                data: data,
                cache: false,
                contentType: false,
                processData: false,
                type: 'PUT',
                headers: headers, 
                success: function(data){
                  attachedfileCollection.fetch({ reset: true });
                  appModel.showModalMessage({
                    title: 'Attached File uploaded',
                    okText: 'ok',
                    body: '"' + file.name + '"'
                  });
                },
                done: function(model, resp){
                  // TODO: done replaces success as of jq 1.8
                  console.log('done');
                },
                error: appModel.jqXHRError
              });
            
              return true;
            }
          },
          view: _form_el,
          title: 'Upload an Attached File'  });
      
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
      
      var PostCollection = Backbone.Collection.extend({
        url: url,
        toJSON: function(){
          return {
            objects: Collection.__super__.toJSON.apply(this) 
          };
        }
      });
      var changedCollection = new PostCollection();
      var MyModel = Backbone.Model.extend({
        url: url,
        initialize : function() {
          this.on('change', function(model, options) {
            // Prevent save on update
            if (options.save === false)
                return;
            model.url = url;
            if(_.isEmpty(model.get('status_date'))){
              model.set('status_date', Iccbl.getISODateString(new Date()));
            }
            if(_.isEmpty(model.get('admin_username'))){
              model.set('admin_username', appModel.getCurrentUser().username);
            }
            changedCollection.add(model);
            appModel.setPagePending();
          });
        },
      });

      var Collection = Iccbl.MyCollection.extend({
        url: url
      });
      collection = new Collection({
        url: url,
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
            
            appModel.clearPagePending();
            
            var errors = form.commit();
            if(!_.isEmpty(errors)){
              console.log('form errors, abort submit: ' + JSON.stringify(errors));
              return false;
            }else{
              var values = form.getValue();
              console.log('form values', values);
              var comments = values['comments'];
              var headers = {};
              headers[appModel.HEADER_APILOG_COMMENT] = comments;
              
              Backbone.sync("patch",changedCollection,
                {
                  headers: headers,
                  error: function(){
                    appModel.jqXHRError.apply(this,arguments);
                    console.log('error, refetch', arguments);
                    changedCollection.reset();
                    collection.fetch({ reset: true });
                  },
                }
              );
            }
          },
          view: _form_el,
          title: 'Save changes?'  
        });
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
    },
    
  });

  return UserView;
});