define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'backbone_forms',
    'iccbl_backgrid',
    'text!templates/generic-edit.html',
    'text!templates/modal_ok_cancel.html',

], function( $, _, Backbone, stickit, backbone_forms, Iccbl, editTemplate, 
            modalOkCancel ) {
	
  // like 'Select' editor, but will always return a boolean (true or
  // false)
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

    // build the backbone-forms schema
    // using custom templates to hold the editors,
    // control the layout with the "controls, control-group" classes
    altFieldTemplate:  _.template('\
      <div class="control-group"> \
            <label class="control-label" for="<%= editorId %>"><%= title %></label>\
            <div class="controls" >\
              <span data-editor></span>\
              <div data-error></div>\
              <div><%= help %></div>\
            </div>\
          </div>\
        '),
    
    // using custom templates to hold the editors,
    // control the layout with the "controls, control-group" classes
    altRadioFieldTemplate: _.template('\
      <div class="control-group"> \
          <label class="control-label" for="<%= editorId %>"><%= title %></label>\
          <div class="controls" >\
            <span data-editor></span>\
            <div data-error></div>\
            <div><%= help %></div>\
          </div>\
        </div>\
      '),
      
    /** 
     * Override the Backbone Forms templateData: this will take the place of
     * the serialize function, since we're overriding the 
     * Backbone Layoutmanager renderTemplate as well.
     */  
    schema: function() {
      var self = this;
      console.log('-----schema called -----');
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
      return editSchema;
    },
  
    /** 
     * Override the Backbone Forms templateData: this will take the place of
     * the serialize function, since we're overriding the 
     * Backbone Layoutmanager renderTemplate as well.
     */    
    templateData: function() {
      console.log('-----templateData called -----');
      var schema = this.model.resource.schema;
      var keys = Iccbl.sortOnOrdinal(
          _.keys(this.model.attributes), schema.fields)
      
      var editKeys = _(keys).filter(function(key){
          return _.has(schema.fields, key) &&
              _.has(schema.fields[key], 'visibility') &&
              _.contains(schema.fields[key]['visibility'], 'edit');
      });
                  
      return {
        'fieldDefinitions': schema.fields,
        'title': Iccbl.getTitleFromTitleAttribute(this.model, schema),
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

    template: _.template(editTemplate)

	});

	return EditView;
});