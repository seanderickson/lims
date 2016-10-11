define([
  'jquery',
  'underscore',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_edit'  
], 
function($, _, Iccbl, appModel, EditView) {

  var UploadDataFormPrototype = {
    
    /**
     * @ returns a Promise wrapping the jQuery.ajax() Promise:
     * jqXHR.done(function( data, textStatus, jqXHR ) {});
     * jqXHR.fail(function( jqXHR, textStatus, errorThrown ) {});
     * jqXHR.always(function( data|jqXHR, textStatus, jqXHR|errorThrown ) { });
     */
    postUploadFileDialog: function(url, content_types, callbackSuccess){

      var promise = $.Deferred();
      

      var self = this;
      var form_template = [
         "<div class='form-horizontal container' id='upload_data_form' >",
         "<form data-fieldsets class='form-horizontal container' >",
         //"<div class='form-group' ><input type='file' name='fileInput' /></div>",
         "</form>",
         "</div>"].join('');      
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
      var choiceHash = {}
      _.each(content_types,function(choice){
        choiceHash[choice] = choice;
      });
      var DisabledField = EditView.DisabledField.extend({
        tagName: 'p',
        className: 'form-control-static'
      });
      
      var formSchema = {};

      formSchema['fileInputPlaceholder'] = {
        title: 'File',
        key: 'file_input',
        type: DisabledField, 
        template: fieldTemplate
      };
      formSchema['type'] = {
        title: 'File Type',
        key: 'type',
        validators: ['required'],
        type: Backbone.Form.editors.Select.extend({
            className: 'form-control'
          }),
        options: choiceHash,
        template: fieldTemplate
      };
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        type: 'TextArea',
        editorClass: 'input-full',
        template: fieldTemplate
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          console.log('form validate', attrs);
          var errs = {};
          var file = $('input[name="fileInput"]')[0].files[0]; 
          if (!file) {
            errs.file_input = 'Must specify a file';
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
      form.$el.find("[name='fileInputPlaceholder']").append(
        '<label class="btn btn-default btn-file">' + 
        'Browse<input type="file" name="fileInput" style="display: none;"></label>'+ 
        '<p id="filechosen" class="form-control-static" ></p>');
      form.$el.on('change', ':file', function() {
        var input = $(this),
            file = input.get(0).files[0],
            label = input.val().replace(/\\/g, '/').replace(/.*\//, '');

        form.$el.find('#filechosen').html(label);
        
        // Try to set the upload type based on the file extension
        // - the upload type is used by the API deserialize method
        var match = /.*\.(.*)$/.exec(file.name)
        if (match && match[1]){
          var extension = match[1].toLowerCase();
          if(_.contains(['xls','xlsx'],extension)){
            form.setValue('type','xls');
          }else if (extension == 'csv'){
            form.setValue('type','csv');
          }else if (extension == 'sdf'){
            form.setValue('type','sdf');
          }else if (extension == 'json'){
            form.setValue('type', 'json');
          }else{
            form.setValue('type',null);
          }
        }
      });

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
            var file = $('input[name="fileInput"]')[0].files[0];
            var values = form.getValue();
            var comments = values['comments'];
            var headers = {};
            headers[appModel.HEADER_APILOG_COMMENT] = comments;
            
            var data = new FormData();
            _.each(_.keys(values), function(key){
              if(values[key]){
                data.append(key,values[key]);
              }
            });
            
            // TODO - add "download for update" option to the download dialog
            
            data.append(values['type'],file);
            
            $.ajax({
              url: url,    
              data: data,
              cache: false,
              contentType: false, // defaults to multipart/form-data when using formdata
              dataType: 'json', // what is expected back from the server
              processData: false, // do not process data being sent
              type: 'POST',
              headers: headers
            }).done(function(data, textStatus, jqXHR){ 
              console.log('success', data);
              
              if (_.isObject(data) && !_.isString(data)){
                data = _.result(_.result(data,'meta',data),'Result',data);
                var msg_rows = appModel.dict_to_rows(data);
                var bodyMsg = msg_rows;
                if (_.isArray(msg_rows) && msg_rows.length > 1){
                  bodyMsg = _.map(msg_rows, function(msg_row){
                    return msg_row.join(': ');
                  }).join('<br>');
                }
                var title = 'Upload success, File: "' + file.name + '"';
                appModel.showModalMessage({
                  body: bodyMsg,
                  title: title  
                });
              }else{
                console.log('data should have been parsed as json', data);
                appModel.showModalMessage({
                  title: 'data file uploaded',
                  okText: 'ok',
                  body: '"' + file.name + '", ' + data
                });
              }
              promise.resolveWith(this,arguments); 
            }).fail(function(jqXHR, textStatus, errorThrown){
              promise.rejectWith(this,arguments);
            })
          
            return true;
          }
        },
        view: _form_el,
        title: 'Upload data file'  
      });
      // return a restricted Promise
      return promise.promise();
    },

    /**
     * Create a dialog to choose the upload file, then execute a jQuery ajax()
     * call to POST the file.
     * @ returns a Promise wrapping the jQuery.ajax() Promise:
     * jqXHR.done(function( data, textStatus, jqXHR ) {});
     * jqXHR.fail(function( jqXHR, textStatus, errorThrown ) {});
     * jqXHR.always(function( data|jqXHR, textStatus, jqXHR|errorThrown ) { });
     */
    uploadAttachedFileDialog: function(url, vocabulary_ref){
  
      var self = this;
      
      var promise = $.Deferred();
      
      var form_template = [
         "<div class='form-horizontal container' id='uploadAttachedFileButton_form' >",
         "<form data-fieldsets class='form-horizontal container' >",
         //"<div class='form-group' ><input type='file' name='fileInput' /></div>",
         "</form>",
         "</div>"].join('');      
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
      var choiceHash = {}
      try{
        var vocabulary = appModel.getVocabulary(vocabulary_ref);
          _.each(_.keys(vocabulary),function(choice){
            choiceHash[choice] = vocabulary[choice].title;
          });
      }catch(e){
        console.log('on get vocabulary', e);
        appModel.error('Error locating vocabulary: ' + vocabulary_ref);
      }
      var DisabledField = EditView.DisabledField.extend({
        tagName: 'p',
        className: 'form-control-static'
      });
      
      var formSchema = {};
  
      formSchema['fileInputPlaceholder'] = {
        title: 'File',
        key: 'file_input',
        type: DisabledField, 
        template: fieldTemplate
      };
      formSchema['type'] = {
        title: 'File Type',
        validators: ['required'],
        key: 'type',
        type: Backbone.Form.editors.Select.extend({
            className: 'form-control'
          }),
        options: choiceHash,
        template: fieldTemplate
      };
      formSchema['file_date'] = {
        title: 'File Date',
        key: 'file_date',
        type: EditView.DatePicker,
        template: fieldTemplate
      };
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        type: 'TextArea',
        editorClass: 'input-full',
        template: fieldTemplate
      };
  
      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs){
          console.log('form validate', attrs);
          var errs = {};
          var file = $('input[name="fileInput"]')[0].files[0]; 
          if (! file) {
            errs.fileInputPlaceholder = 'Must specify a file';
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
      form.$el.find("[name='fileInputPlaceholder']").append(
        '<label class="btn btn-default btn-file">' + 
        'Browse<input type="file" name="fileInput" style="display: none;"></label>'+ 
        '<p id="filechosen" class="form-control-static" ></p>');
      form.$el.on('change', ':file', function() {
        var input = $(this),
            numFiles = input.get(0).files ? input.get(0).files.length : 1,
            label = input.val().replace(/\\/g, '/').replace(/.*\//, '');
        form.$el.find('#filechosen').html(label);
        var file = $('input[name="fileInput"]')[0].files[0];
        if (file) {
          form.setValue(
            'file_date',
            _.result(file,'lastModifiedDate', null));
        }
      });
  
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
            var comments = values['comments'];
            var headers = {};
            headers[appModel.HEADER_APILOG_COMMENT] = comments;
            
            var data = new FormData();
            _.each(_.keys(values), function(key){
              if(values[key]){
                data.append(key,values[key]);
              }
            });
            
            var file = $('input[name="fileInput"]')[0].files[0];
            data.append('attached_file',file);
            filename = file.name;
            if (_.isEmpty(values['file_date'])){
              data.append(
                'file_date',
                Iccbl.getISODateString(_.result(file,'lastModifiedDate', null)));
            }
            
            $.ajax({
              url: url,    
              data: data,
              cache: false,
              contentType: false,
              processData: false,
              dataType: 'json', // what is expected back from the server
              type: 'POST',
              headers: headers
            }).done(function(data, textStatus, jqXHR){ 
              promise.resolveWith(this,arguments); 
            }).fail(function(jqXHR, textStatus, errorThrown){
              promise.rejectWith(this,arguments);
            })
            return true;
          }
        },
        view: _form_el,
        title: 'Upload an Attached File'  
      });
      
      // return a restricted Promise
      return promise.promise();
    }
    
  };

  return UploadDataFormPrototype;
});



