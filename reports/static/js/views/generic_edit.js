define([
    'jquery',
    'underscore',
    'backbone',
//    'backbone_stickit',
    'backbone_forms',
    'iccbl_backgrid',
    'models/app_state',
    'text!templates/generic-edit.html',
    'text!templates/modal_ok_cancel.html',

], function( $, _, Backbone, backbone_forms, Iccbl, appModel,
            editTemplate, modalOkCancel ) {
	
  // like 'Select' editor, but will always return a boolean (true or false)
  Backbone.Form.editors.BooleanSelect = Backbone.Form.editors.Select.extend({
    
    initialize: function(options) {
      options.schema.options = [{ val: 'true', label: 'Yes' },
                                { val: 'false', label: 'No' }
                                ];
      Backbone.Form.editors.Select.prototype.initialize.call(this, options);
    },
    getValue: function() {
      return !!Backbone.Form.editors.Select.prototype.getValue.call(this);
    },
    setValue: function(value) {
      value = value ? 'true' : 'false';
      Backbone.Form.editors.Select.prototype.setValue.call(this, value);
    }
  });            
  
  
  var EditView = Backbone.Form.extend({

    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    
    initialize: function(args) {
      console.log('---- initialize EditView');
      this.uriStack = args.uriStack;
      this.consumedStack = []; 
      
      Backbone.Form.prototype.initialize.apply(this,args);
      
//      this.model.on('change', function(evt){
//        console.log('chg evt: ' + evt);
//      });
      this.listenTo(this, 'screen_type:change', function(target){
        console.log('chg evt:' + target);
////        if(evt) evt.preventDefault();
//        $( "#c32_screen_type-1" ).prop( "checked", true );
      });
      this.on('change', function(t){
        console.log('all:' + t);
      });
    },

    events: {
      'click button#save': 'save'
    },
    
    // build the backbone-forms schema
    // using custom templates to hold the editors,
    // control the layout with the "controls, control-group" classes
    altFieldTemplate:  _.template('\
      <div class="form-group" > \
            <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
            <div class="col-sm-10" >\
              <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />\
              <div data-error class="text-danger" ></div>\
              <div><%= help %></div>\
            </div> \
          </div>\
        '),
    
    // using custom templates to hold the editors,
    // control the layout with the "controls, control-group" classes
    altRadioFieldTemplate: _.template('\
      <div class="form-group"  ><fieldset> \
          <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
          <div class="col-sm-10" >\
            <div data-editor  ></div>\
            <div data-error></div>\
            <div><%= help %></div>\
          </div></fieldset>\
        </div>\
      '),
      
    schema: function() {
      var self = this;
      var schema = this.model.resource.schema;
      var keys = Iccbl.sortOnOrdinal(
          _.keys(this.model.attributes), schema.fields)
      
      var editKeys = _(keys).filter(function(key){
          return _.has(schema.fields, key) &&
              _.has(schema.fields[key], 'visibility') &&
              _.contains(schema.fields[key]['visibility'], 'edit');
      });
                  
      var editSchema = {};
      var itemcount = 0;
      // process the ui_types - convert to backbone-forms schema editor type
      _.each(editKeys, function(key){
        if( _(schema.fields).has(key)){
          var option = schema.fields[key];
//          console.log('key: ' + key + ', option: ' + JSON.stringify(option));
          var fieldSchema = editSchema[key] = {};

          var typeMap = {
              'boolean': 'Checkbox',
              'string': 'Text',
              'uri': 'Text',
              'float': 'Number',
              'integer': 'Number'
          }
          
          var temp = option.ui_type || 'Text';
          temp = temp.toLowerCase();
          fieldSchema['type'] = temp.charAt(0).toUpperCase() + temp.slice(1);
          if(_.has(typeMap, option.ui_type)){
            fieldSchema['type'] = typeMap[option.ui_type];
          }
          var validators = [];
          
          if(option.ui_type == 'Select' 
              || option.ui_type == 'Radio'
              || option.ui_type == 'Checkboxes' ){
            fieldSchema['options'] = option.choices || [];
            if(_.isEmpty(option.choices)){
              appModel.error('no choices defined for: ' + key);
            }
          }
          
          if(option.ui_type == 'Radio'){
            editSchema[key]['template'] = self.altRadioFieldTemplate;
          }else{
            editSchema[key]['template'] = self.altFieldTemplate;
          }

          // validation stuff
          if(fieldSchema['type']  == 'Number'){
            // TODO: check for the "min" "max","range" validation properties and implement
            if( !_.isUndefined(option.min)){
              var validator = function checkMin(value, formValues) {
                var err = {
                    type: 'Min',
                    message: 'must be >= ' + option.min
                };

                if (value <= option.min ) return err;
              };
              validators.unshift(validator);
            }
            if( !_.isUndefined(option.range)){
              var validator = function checkRange(value, formValues) {
                
                var last = '';
                var rangeMsg = '';
                var value_ok = false;
                var schema_lower = 0, schema_upper = 0
                for(var i=0; i<option.range.length; i++){
                  var schema_val = option.range[i]
                  if(i>0) rangeMsg += ', ';
                  if(i%2 == 0){
                    rangeMsg += '> ' + schema_val
                  }else{
                    rangeMsg += '< ' + schema_val
                  }
                }
                // compare range in pairs
                for(var i=0; i<option.range.length; i+=2){
                  schema_lower = parseInt(option.range[i])
                  if(option.range.length > i+1){
                    schema_upper = parseInt(option.range[i+1])
                    if(value >schema_lower && value<schema_upper){
                      value_ok = true;
                      break;
                    }
                  }else{
                    if(value > schema_lower){
                      value_ok = true;
                      break; // not nec
                    }
                  }
                }
                var result = {
                    type: 'Range',
                    message: 'value not in range: ' + value + ' range: ' + rangeMsg
                };

                if (!value_ok) return result;
              };
              validators.unshift(validator);
            }
          }
          if(option.required){
            validators.unshift('required');
          }
          if(!_.isUndefined(option.regex) && !_.isEmpty(option.regex)){
            var validator = { type: 'regexp', regexp: new RegExp(option.regex) };
            if(!_.isUndefined(option.validation_message) && !_.isEmpty(option.validation_message)){
              validator.message = option.validation_message;
              // TODO: rework, if req'd, to use tokenized strings (will need 
              // to reimplement backbone-forms
              //  function(value, formValues){
              //    //                TODO: figure out how to get the pending model
              //    return 'value: ' + value + ' is incorrect: ' + Iccbl.replaceTokens(new Backbone.Model(formValues), option.validation_message);
              //  };
            }
            validators.unshift(validator);
          }
          if(!_.isEmpty(validators)){
            editSchema[key].validators = validators;
          }
          
          editSchema[key]['maxlength'] = 50;

          if(itemcount++ == 0){
            // Set autofocus (HTML5) on the first field
            // NOTE: see
            // http://stackoverflow.com/questions/20457902/how-to-automatically-focus-first-backbone-forms-input-field
            // - we may want to revisit this for a more robust solution
            editSchema[key]['editorAttrs'] = { autofocus: 'autofocus'}
          }
        }
        console.log('editSchema[' + key + '] = ' + JSON.stringify(editSchema[key]));
      });      
      
      // Note: Enforced comment
      editSchema['comment'] = {
          type: 'TextArea',
          validators: ['required'], 
          template: self.altFieldTemplate
      };
      
      return editSchema;
    },
  
    /** 
     * Override the Backbone Forms templateData: this will take the place of
     * the serialize function, since we're overriding the 
     * Backbone Layoutmanager renderTemplate as well.
     */    
    templateData: function() {
      var schema = this.model.resource.schema;
      var keys = Iccbl.sortOnOrdinal(
          _.keys(this.model.attributes), schema.fields)
      
      var editKeys = _(keys).filter(function(key){
          return _.has(schema.fields, key) &&
              _.has(schema.fields[key], 'visibility') &&
              _.contains(schema.fields[key]['visibility'], 'edit');
      });
      
      editKeys.push('comment');
      schema.fields['comment'] = { key: 'comment', title: 'Comment', ui_type:'string'};
                  
      return {
        'fieldDefinitions': schema.fields,
        'keys': _.chain(editKeys)
      };      
    },	

    /** 
     * Override the Backbone Layoutmanager template rendering to use 
     * Backbone Forms
     */
    renderTemplate: function() {
      return Backbone.Form.prototype.render.apply(this);
    },
    
//    afterRender: function(){
//      console.log('afterRender called');
//      
//      $( "#c32_screen_type-1" ).change(function(e){
//        console.log('on change called: ' + e);
//      });
//
//    },

    template: _.template(editTemplate),
    
    save: function( event ) {
      event.preventDefault();
      var self = this;
      
      var errors = this.commit();
      
      if(errors){
        console.log(JSON.stringify(errors));
        _.each(_.keys(errors), function(key){
          var error = errors[key];
          
          $('[name="'+key +'"').parents('.form-group').addClass('has-error');
        });
        return;
      }
      
      // Fixup the URL - if it points to the model instance, make it point to 
      // the API resource only: tastypie wants this
      // Note: this is happening if the model was fetched specifically for this
      // page, and has the url used to fetch it, rather than the collection url.
      var key = Iccbl.getIdFromIdAttribute( self.model,self.model.resource.schema );
      var url = _.result(this.model, 'url');
      ////    if ( url && url.indexOf(key) != -1 ) {
      ////    url = url.substring( 0,url.indexOf(key) );
      ////  }    
      
      // TODO: this should be standard to have url end with '/'
      if( url && url.charAt(url.length-1) != '/'){
        url += '/';
      }
      
      var _patch = true;
      if (_.contains(this.uriStack, '+add')){
        _patch = false;
        // NOTE: don't set the key, since this is a create/POST to the resource URL
        //        url += key;
      }else{
        // TODO: check if creating new or updating here
        // set the id specifically on the model: backbone requires this to 
        // determine whether a "POST" or "PATCH" will be used
        this.model.id = key;
      }
      
      this.model.save(null, {
        url: url, // set the url property explicitly
        patch: _patch,
        // Note:
        // You have to send { dataType: 'text' } to have the success function 
        // work with jQuery and empty responses ( otherwise, fails on JSON.parse 
        // of the empty response).        
        //        dataType: 'text', 
        // The other solution: use "always_return_data" in the tastypie
        // resource definitions - which we are doing.
        headers: {
          'APILOG_COMMENT': self.model.get('comment')
        }
      })
      .success(function(model, resp){
        // note, not a real backbone model, just JSON
        model = new Backbone.Model(model);
        var key = Iccbl.getIdFromIdAttribute( model,self.model.resource.schema );
        //        appModel.set('routing_options', {trigger: true});
        //        self.reportUriStack([key]);
        appModel.router.navigate(self.model.resource.key + '/' + key, {trigger:true});
      })
      .done(function(model, resp){
        // TODO: done replaces success as of jq 1.8
        console.log('model saved');
      })
      .error(appModel.jqXHRFail)
      .always(function() {
        // always replaces complete as of jquery 1.8
        self.trigger('remove');
      });
    }

	});

	return EditView;
});