define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'backbone_forms',
    'iccbl_backgrid',
    'models/app_state',
    'text!templates/generic-edit.html',
    'text!templates/modal_ok_cancel.html',

], function( $, _, Backbone, stickit, backbone_forms, Iccbl, appModel,
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
    },

    events: {
      'click button#save': 'save'
    },
    
    // build the backbone-forms schema
    // using custom templates to hold the editors,
    // control the layout with the "controls, control-group" classes
    altFieldTemplate:  _.template('\
      <div class="form-group"> \
            <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
            <div class="col-sm-10" >\
              <div data-editor></span>\
              <div data-error class="text-danger" ></div>\
              <div><%= help %></div>\
            </div>\
          </div>\
        '),
    
    // using custom templates to hold the editors,
    // control the layout with the "controls, control-group" classes
    altRadioFieldTemplate: _.template('\
      <div class="form-group"> \
          <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
          <div class="col-sm-10" >\
            <span data-editor></span>\
            <div data-error></div>\
            <div><%= help %></div>\
          </div>\
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
                  
      // TODO: memoize?
      var editSchema = {};
      var itemcount = 0;
      _.each(editKeys, function(key){
        if( _(schema.fields).has(key)){
          option = schema.fields[key];
          if(option.ui_type == 'Select' 
              || option.ui_type == 'Radio'
              || option.ui_type == 'Checkboxes' ){
            var _optionsCollection = [];
            if(_.has(option, 'choices')){
              _optionsCollection = option.choices.map(function(choice){
                return choice;
              });
            }else{
              // TODO: use bootstrap alerts
              window.alert('Warning, no choices defined for: ' + key);
              // set choices so the template doesn't complain
              option.choices = _(_optionsCollection); 
            }
            editSchema[key] = { 
              type: option.ui_type, 
              options: _optionsCollection
            };

            if(option.ui_type == 'Checkboxes' ){ 
              editSchema[key]['defaults'] = self.model.get(key);
            } 
          }else if( option.ui_type == 'boolean'){
            editSchema[key] = {
              type: 'Checkbox'
            };
          }else if( option.ui_type.toLowerCase() == 'date'){
            editSchema[key] = {
              type: 'Date'
            };
          }else{
            editSchema[key] = {
              type: 'Text'
            };
          }
          if(option.ui_type == 'Radio'){
            // editSchema[key]['template'] = self.altRadioFieldTemplate;
          }else{
            editSchema[key]['template'] = self.altFieldTemplate;
          }

          if(itemcount++ == 0){
            // Set autofocus (HTML5) on the first field
            // NOTE: see
            // http://stackoverflow.com/questions/20457902/how-to-automatically-focus-first-backbone-forms-input-field
            // - we may want to revisit this for a more robust solution
            editSchema[key]['editorAttrs'] = { autofocus: 'autofocus'}
          }
        }
      });      
      
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