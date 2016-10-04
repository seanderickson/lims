define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_forms',
    'multiselect',
    'quicksearch',
    'bootstrap-3-typeahead',
    'iccbl_backgrid',
    'models/app_state',
    'views/generic_detail_stickit',
    'templates/generic-edit.html',
    'bootstrap',
    'bootstrap-datepicker',
    'chosen'
], function( $, _, Backbone, backbone_forms, 
            multiselect, quicksearch, typeahead, Iccbl, appModel,
            DetailView, editTemplate ) {
  var SIunitEditor = Backbone.Form.editors.Base.extend({
    
    tagname: 'siuniteditor',
    
    fieldTemplate: _.template([
      '<div data-editor class="form-control col-sm-10" title="<%= help %>"  >'
      ].join('')),
    unitFieldTemplate: _.template([
      '<div data-editor title="<%= help %>"  >'
      ].join('')),

    formTemplate: _.template([
      '<form class="iccbl-nested-form" >',
      '<div class="input-group ">',
      '   <div data-fields="number" />',
      '   <div class="input-group-addon iccbl-nested-form-unit" data-fields="unit"/>',
      '</div>',
      '</form>'
    ].join('')),

    siunits: [
      ['T', 1e12],
      ['G', 1e9],
      ['M', 1e6],
      ['k', 1e3],
      ['', 1],
      ['m', 1e-3],
      ['μ', 1e-6],
      ['n', 1e-9 ],
      ['p', 1e-12 ]
      ],

    events: {
        'change': function() {
            // The 'change' event should be triggered whenever something happens
            // that affects the result of `this.getValue()`.
            this.trigger('change', this);
        },
        'focus': function() {
            // The 'focus' event should be triggered whenever an input within
            // this editor becomes the `document.activeElement`.
            this.trigger('focus', this);
            // This call automatically sets `this.hasFocus` to `true`.
        },
        'blur': function() {
            // The 'blur' event should be triggered whenever an input within
            // this editor stops being the `document.activeElement`.
            this.trigger('blur', this);
            // This call automatically sets `this.hasFocus` to `false`.
        }
    },

    initialize: function(options) {
      var self = this;
      Backbone.Form.editors.Base.prototype.initialize.call(this, options);

      var _options = typeof options !== 'undefined' ?  options : {};
      options = _options;
      this.options = options;
      
      var formSchema = this.formSchema = {};
      var multiplier = this.multiplier = options.schema.multiplier || 1;
      var symbol = this.symbol = options.schema.symbol;
      var decimals = this.decimals = options.schema.decimals || 0;
      var defaultUnit = this.defaultUnit = options.schema.defaultUnit || 1;
      var units = this.units = [];
      
      if(! symbol){
        throw 'Error: SIUnitFormFilter requires a "symbol" option'; 
      }
      _.each(this.siunits,function(pair){
        if(options.schema.maxunit){
          if(options.schema.maxunit < pair[1]) return;
        }
        if(options.schema.minunit){
          if(options.schema.minunit > pair[1]) return;
        }
        units.push({ val: pair[1], label: pair[0] + self.symbol });
      });
      formSchema['unit'] = {
        title: '', 
        key:  'unit', 
        type: 'Select',
        options: units,
        editorClass: 'form-control',
        template: this.unitFieldTemplate
      };
      formSchema['number'] = {
        title: '', 
        key:  'number',
        type: 'Number',
        editorClass: 'form-control',
        template: this.fieldTemplate
      };
    },

    focus: function() {
      if (this.hasFocus) return;
      this.$el.focus();
    },

    blur: function() {
      if (!this.hasFocus) return;
      this.$el.blur();
    },
    
    render: function() {
      var formModel;
      if(this.value){
        formModel= new Backbone.Model(this._findNumberAndUnit(this.value));
      }else{
        formModel= new Backbone.Model({ number: 0, unit: this.defaultUnit });
      }
      this.nestedForm = new Backbone.Form({
        schema: this.formSchema,
        model: formModel,
        idPrefix: this.id + '_',
        Field: Backbone.Form.NestedField,
        template: this.formTemplate
      });
      this._observeFormEvents();
      this.$el.html(this.nestedForm.render().el);
      if (this.hasFocus) this.trigger('blur', this);
      return this;
    },
    
    
    /**
     * Returns the current editor value
     * @return {String}
     */
    getValue: function() {
      var self = this;
      if (this.nestedForm){
        return this._calculate(
            self.multiplier,
            this.nestedForm.getValue()['unit'],
            this.nestedForm.getValue()['number']);
      }
  
      return this.value;
    },

    setValue: function(value) {
      this.value = value;
      this.render();
    },

    _calculate: function(multiplier, si_mult, val){
      // run strip after every calculation to round out floating point math errors
      function strip(number) {
        return (parseFloat(number.toPrecision(12)));
      }
      val = strip(val * multiplier);
      if(si_mult > 0){ // if si unit is undefined, assume to be 1
        val = strip(val * si_mult);
      }
      return val;
    },
    
    _findNumberAndUnit: function(number){
      var self = this;
      function strip(number) {
        return (parseFloat(number.toPrecision(12)));
      }
      number = strip(number/self.multiplier);
      pair = _.find(this.siunits, function(pair){
        return pair[1] <= Math.abs(number); 
      });
      
      if(_.isUndefined(pair)){
        console.log('could not find units for the input number: ' + number);
        return { number:number, unit: ''};
      }
      
      var val = (1/pair[1])*number;
      val = Math.round(val*Math.pow(10,this.decimals))/Math.pow(10,this.decimals);
      return {number:val, unit: pair[1]};
    },
  
    remove: function() {
      this.nestedForm.remove();
      Backbone.View.prototype.remove.call(this);
    },
  
    validate: function() {
      return this.nestedForm.validate();
    },
  
    _observeFormEvents: function() {
      // redirect all of the nested events up
      if (!this.nestedForm) return;
      this.nestedForm.on('all', function() {
        // args = ["key:change", form, fieldEditor]
        var args = _.toArray(arguments);
        args[1] = this;
        // args = ["key:change", this=objectEditor, fieldEditor]
        this.trigger.apply(this, args);
      }, this);
    },
    
    /**
     * Update the embedded model, checking for nested validation errors and pass them up
     * Then update the main model if all OK
     *
     * @return {Error|null} Validation error or null
     */
    commit: function() {
      console.log('commit');
      var error = this.nestedForm.commit();
      if (error) {
        this.$el.addClass('error');
        return error;
      }
    }
    
  });
  
  var DatePicker = Backbone.Form.editors.Base.extend({

    tagName: 'datepicker',

    events: {
        'change': function() {
            // The 'change' event should be triggered whenever something happens
            // that affects the result of `this.getValue()`.
            this.trigger('change', this);
        },
        'focus': function() {
            // The 'focus' event should be triggered whenever an input within
            // this editor becomes the `document.activeElement`.
            this.trigger('focus', this);
            // This call automatically sets `this.hasFocus` to `true`.
        },
        'blur': function() {
            // The 'blur' event should be triggered whenever an input within
            // this editor stops being the `document.activeElement`.
            this.trigger('blur', this);
            // This call automatically sets `this.hasFocus` to `false`.
        }
    },

    initialize: function(options) {
        Backbone.Form.editors.Base.prototype.initialize.call(this, options);
        if (this.value){
          this.value = new Date(this.value);
        }
    },

    getValue: function() {
      var input = $('input', this.el);
      var date = input.datepicker('getDate');
      var val = Iccbl.getISODateString(date);
      return val;
    },

    setValue: function(value) {
      $('input', this.el).datepicker('setUTCDate', value);
    },

    focus: function() {
        if (this.hasFocus) return;
        this.$el.focus();
    },

    blur: function() {
        if (!this.hasFocus) return;
        this.$el.blur();
    },
    
    render: function() {
      var el = $(this.el);
      el.html([
         '<div class="input-group date" >',
         '  <input type="text" class="form-control">',
         '  <span class="input-group-addon" id="datepicker-icon" >',
         '    <i class="glyphicon glyphicon-th"  ></i>',
         '  </span>',
         '</div>'                                             
      ].join(''));

      var input = $('input', el);
      input.datepicker({
          dateFormat: 'dd/mm/yyyy',
          autoclose: true,
          todayBtn: 'linked',
          todayHighlight: true,
          orientation: "bottom auto"
      });
      
      // manually get the input-group-addon click
      $('#datepicker-icon',el).click(function(e){
        input.datepicker().focus();
      });
      this.setValue(this.value);
        
      return this;
    }
    
  });  

  var DisabledField = Backbone.Form.editors.Text.extend({

    tagName: 'p',

    initialize: function(options) {
        Backbone.Form.editors.Text.prototype.initialize.call(this, options);
    },

    getValue: function() {
      return this.value;
    },
    setValue: function(value) {
      this.value = value;
    },

    _getValue: function() {
      if(this.binding){
        var val = this.binding.onGet(this.value);
        if(_.isBoolean(val)){
          if(! val) return 'false';
        }
        return val;
      }else{
        console.log('no view binding found for disabled field, value', this.value);
        if(_.isArray(this.value)){
          return this.value.join(', ');
        }
        return this.value;
      }
    },

    render: function() {
      this.$el.html(this._getValue());
      return this;
    }
  });  

  // like 'Select' editor, but will always return a boolean (true or false)
  var BooleanSelect = Backbone.Form.editors.Select.extend({
    
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
  
  // FIXME: 20150624 - overriding the checkbox editor, to use "div" instead of "<ul><li>
  // - something in CSS/JS is overriding the click function and not allowing 
  // selection of checkboxes in the <ul>,<li> elements
  var Checkboxes = Backbone.Form.editors.Checkboxes.extend({

    tagName: 'div',

    groupNumber: 0,

    events: {
      'click input[type=checkbox]': function() {
        this.trigger('change', this);
      },
      'focus input[type=checkbox]': function() {
        if (this.hasFocus) return;
        this.trigger('focus', this);
      },
      'blur input[type=checkbox]':  function() {
        if (!this.hasFocus) return;
        var self = this;
        setTimeout(function() {
          if (self.$('input[type=checkbox]:focus')[0]) return;
          self.trigger('blur', self);
        }, 0);
      }
    },

    getValue: function() {
      var values = [];
      this.$('input[type=checkbox]:checked').each(function() {
        values.push($(this).val());
      });
      return values;
    },

    setValue: function(values) {
      if (!_.isArray(values)) values = [values];
      this.$('input[type=checkbox]').val(values);
    },

    focus: function() {
      if (this.hasFocus) return;

      this.$('input[type=checkbox]').first().focus();
    },

    blur: function() {
      if (!this.hasFocus) return;

      this.$('input[type=checkbox]:focus').blur();
    },

    /**
     * Create the checkbox list HTML
     * @param {Array}   Options as a simple array e.g. ['option1', 'option2']
     *                      or as an array of objects e.g. [{val: 543, label: 'Title for object 543'}]
     * @return {String} HTML
     */
    _arrayToHtml: function (array) {
      var html = $();
      var self = this;

      _.each(array, function(option, index) {
        var itemHtml = $('<div>');
        if (_.isObject(option)) {
          if (option.group) {
            var originalId = self.id;
            self.id += "-" + self.groupNumber++;
            itemHtml = $('<fieldset class="group">').append( $('<legend>').text(option.group) );
            itemHtml = itemHtml.append( self._arrayToHtml(option.options) );
            self.id = originalId;
            close = false;
          }else{
            var val = (option.val || option.val === 0) ? option.val : '';
            itemHtml.append( $('<input type="checkbox" name="'+self.getName()+'" id="'+self.id+'-'+index+'" />').val(val) );
            if (option.labelHTML){
              itemHtml.append( $('<label for="'+self.id+'-'+index+'">').html(option.labelHTML) );
            }
            else {
              itemHtml.append( $('<label for="'+self.id+'-'+index+'">').text(option.label) );
            }
          }
        }
        else {
          itemHtml.append( $('<input type="checkbox" name="'+self.getName()+'" id="'+self.id+'-'+index+'" />').val(option) );
          itemHtml.append( $('<label for="'+self.id+'-'+index+'">').text(option) );
        }
        html = html.add(itemHtml);
      });

      return html;
    }

  });
  
  var MultiSelect2 = Backbone.Form.editors.Select.extend({
    render: function() {
      this.$el.attr('multiple', 'multiple');
      this.setOptions(this.schema.options);
      return this;
    }
  });
  
  var ChosenMultiSelect = Backbone.Form.editors.Select.extend({
    render: function() {
      Backbone.Form.editors.Select.prototype.render.apply(this);
      
      this.$el.attr('multiple', 'multiple');
      this.setOptions(this.schema.options);
      
      if(this.schema.placeholder){
        this.$el.attr('data-placeholder', this.schema.placeholder);
      }
      return this;
    }
  });
  
  var ChosenSelect = Backbone.Form.editors.Select.extend({
    render: function() {
      var self = this;
      Backbone.Form.editors.Select.prototype.render.apply(this);
      
      this.setOptions(this.schema.options);
      if(this.schema.placeholder){
        this.$el.attr('data-placeholder', this.schema.placeholder);
      }
      
      // if the current value is not in the options, then display it as the placeholder
      var currVal = this.model.get(this.key);
      if(!_.find(this.schema.options,function(option){
        return option.val == currVal;
      })){
        this.$el.attr('data-placeholder', currVal);
      }else if(this.schema.placeholder){
        this.$el.attr('data-placeholder', this.schema.placeholder);
      }
      return this;
    }
  });
  
  
  var EditView = Backbone.Form.extend({
    initialize: function(args) {

      console.log('Editview initialize: ', args);
      var self = this;
      
      this.uriStack = args.uriStack;
      this.consumedStack = []; 
      this.saveCallBack = args.saveCallBack;

      this.modelSchema = args.modelSchema || this.model.resource;
      this.modelFields = args.modelFields || this.modelSchema.fields;
      this.editKeys = args.editKeys || this.modelSchema.allEditVisibleKeys();
      this.groupedKeys = this.modelSchema.groupedKeys(this.editKeys);
      this.editableKeys = args.editableKeys || this.modelSchema.updateKeys();
      if(args.isCreate 
          || _.isEmpty(_.compact(_.values(this.model.pick(this.model.resource['id_attribute']))))){
        this.editableKeys = _.union(this.editableKeys,this.modelSchema.createKeys());
      }
      // Add comment field
      if( _.propertyOf(this.model.resource,'require_comment_on_save') && 
          !_.contains(this.editableKeys, 'apilog_comment')){
        this.editableKeys.push('apilog_comment');
        this.modelFields = _.extend({
          'apilog_comment': {
            edit_type: 'textarea',
            title: 'Changelog Comment',
            is_required: true
          }
        }, this.modelFields);
      }
      this.finalEditableKeys = [];
      
      // look in the uriStack for preset values
      console.log('Editview; uriStack', this.uriStack);
      if(this.uriStack){
        for(var i=0; i<this.uriStack.length; ){
          var key = this.uriStack[i++];
          if (i<this.uriStack.length){
            var val = this.uriStack[i++];
            console.log('key', key, 'val', val);
            this.model.set(key,val);
          }
        }
      }
      
      // The delegateModel/View is used to display visible, but not editable fields
      var delegateModel = new Backbone.Model(_.clone(this.model.attributes ));
      delegateModel.resource = this.model.resource;
      this.delegateDetailView = new DetailView({ model: delegateModel });
      
      if (args.editTemplate){
        this.template = _.template(args.editTemplate);
      }else{
        this.template = _.template(editTemplate);
      }
      
      // NOTE: due to a phantomjs js bug, must convert arguments to a real array
      Backbone.Form.prototype.initialize.apply(this,Array.prototype.slice.apply(arguments));
      
    },
    
    events: {
      'click button#save': 'save',
      'click button#cancel': 'cancel'
    },
    
    altFieldTemplate:  _.template([
      '<div class="form-group" key="form-group-<%=key%>">',
      '    <label class="control-label col-sm-2" for="<%= editorId %>" ><%= title %></label>',
      '    <div class="<%= editorAttrs.widthClass%>" >',
      '      <div data-editor  key="<%=key%>" style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
      '      <div data-error class="text-danger" ></div>',
      '      <div><%= help %></div>',
      '    </div>',
      '  </div>'
    ].join('')),

    altTextAreaFieldTemplate:  _.template([
      '<div class="form-group" >',
      '    <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>',
      '    <div class="<%= editorAttrs.widthClass%>" >',
      '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
      '      <div data-error class="text-danger" ></div>',
      '      <div><%= help %></div>',
      '    </div>',
      '  </div>'
    ].join('')),

    altMultiselectFieldTemplate:  _.template([
      '<div class="form-group" >',
      '    <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>',
      '    <div class="<%= editorAttrs.widthClass%>" >',
      '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
      '      <div data-error class="text-danger" ></div>',
      '      <div><%= help %></div>',
      '    </div>',
      '  </div>'
    ].join('')),

    altMultiselect2FieldTemplate:  _.template('\
      <div class="form-group" > \
          <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
          <div class="<%= editorAttrs.widthClass%>" >\
            <div data-editor  class="multiselect2" style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />\
            <div data-error class="text-danger" ></div>\
            <div><%= help %></div>\
          </div> \
        </div>\
    '),

    altRadioFieldTemplate: _.template('\
      <div class="form-group"  ><fieldset> \
          <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
          <div class="<%= editorAttrs.widthClass%>" >\
            <div data-editor  ></div>\
            <div data-error></div>\
            <div><%= help %></div>\
          </div></fieldset>\
        </div>\
    '),

    schema: function() {
      
      var self = this;
      var editSchema = {};
      
      var typeMap = {
        'custom': 
          {
            type: DisabledField, // presumably, the custom view will replace this 
            editorAttrs: { widthClass: 'col-sm-6' }
          },
        'boolean': 
          {
            type: Backbone.Form.editors.Checkbox,
            editorAttrs: { widthClass: 'col-sm-10'}
          },
        'string': 
          {
            type: Backbone.Form.editors.Text,
            editorClass: 'form-control',
            editorAttrs: { widthClass: 'col-sm-4'}
          },
        'uri': 
          {
            type: Backbone.Form.editors.Text,
            editorAttrs: { widthClass: 'col-sm-4'}
          },
        'float': 
          {
            type: Backbone.Form.editors.Number,
            editorClass: 'form-control',
            editorAttrs: { widthClass: 'col-sm-4'}
          },
        'integer': 
          {
            type: Backbone.Form.editors.Number,
            editorClass: 'form-control',
            editorAttrs: { widthClass: 'col-sm-4'}
          },
        'decimal': 
          {
            type: Backbone.Form.editors.Number,
            editorClass: 'form-control',
            editorAttrs: { widthClass: 'col-sm-4'}
          },
        'list': 
          {
            // NOTE: list will only be a multi-select if the edit type is set 
            // TODO:test 20160202, added as a test with the plate ranges, but useful for other list fields
            type: Backbone.Form.editors.List, 
            editorAttrs: { widthClass: 'col-sm-6' }
          },
        'date': 
          {
            type: DatePicker,
            editorAttrs: { widthClass: 'col-sm-3'}
          },
        'textarea':
          {
            type: Backbone.Form.editors.TextArea.extend({
              initialize: function(){
                Backbone.Form.editors.TextArea.prototype.initialize.apply(this, arguments);
                if(_.has(this.schema,'rows')){
                  this.$el.attr('rows',this.schema['rows']);
                }
              }
            }),
            editorClass: 'form-control',
            template: self.altTextAreaFieldTemplate,
            editorAttrs: { widthClass: 'col-sm-10'} 
          },
        'siunit': 
          {
            type: SIunitEditor,
            editorClass: '',
            editorAttrs: { widthClass: 'col-sm-4'}
          },
        'select':
          {
            type: ChosenSelect,
            editorClass: 'chosen-select',
            editorAttrs: { widthClass: 'col-sm-5'}

          },
        'multiselect':
          {
            type: Checkboxes,
            template: self.altMultiselectFieldTemplate,
            editorAttrs: { widthClass: 'col-sm-10'}
          },
        'multiselect2':
          {
            type: MultiSelect2,
            template: self.altMultiselect2FieldTemplate,
            editorAttrs: { widthClass: 'col-sm-10'}
          },
        'multiselect3':
          {
            type: ChosenMultiSelect,
            editorClass: 'chosen-select',
            editorAttrs: { widthClass: 'col-sm-5'}
          },
        'typeahead': 
          {
            type: Backbone.Form.editors.Text,
            editorClass: 'form-control',
            editorAttrs: { widthClass: 'col-sm-4'}
          },
        'radio':
          {
            template: self.altRadioFieldTemplate,
            editorAttrs: { widthClass: 'col-sm-10'}
          }
      };
      var defaultFieldSchema = 
        {
           template: self.altFieldTemplate,
           editorAttrs: { widthClass: 'col-sm-10', maxlength: 50},
           fieldAttrs: {}
         };
      
      _.each(this.editKeys, function(key){

        var fi = self.modelFields[key];
        var cell_options = fi.display_options;

        var fieldSchema = editSchema[key] = _.extend({}, defaultFieldSchema);
        
        fieldSchema['title'] = fi.title;
        var tooltip = fi.description;
        if (fi.required){
          tooltip += ' (required)';
        }
        fieldSchema['fieldAttrs'] = { title: tooltip };
        
        if(!_.contains(self.editableKeys, key)){
        
          console.log('create disabled entry', key, fi['editability'])
          fieldSchema['type'] = DisabledField.extend({
            initialize: function(){
              DisabledField.__super__.initialize.apply(this,arguments);
              this.binding = self.delegateDetailView.createBinding(key,fi)
            }
          });
          fieldSchema['editorClass'] = 'form-control-disabled';

        }else{
          
          console.log('build edit schema for key: ',key);
          self.finalEditableKeys.push(key);
          if(_.has(typeMap, fi.data_type)){
            _.extend(fieldSchema, typeMap[fi.data_type]);
          }
          if(_.has(typeMap, fi.display_type)){
            _.extend(fieldSchema, typeMap[fi.display_type]);
          }
          if(_.has(typeMap, fi.edit_type)){
            _.extend(fieldSchema, typeMap[fi.edit_type]);
          }
          
          if (cell_options){
            _.extend(fieldSchema,cell_options);
            // editorAttrs are available as data for the template compilation
            // TODO: create editorAttrs nested property on cell_options to clean up html
            fieldSchema.editorAttrs = _.extend({},fieldSchema.editorAttrs,cell_options)
          }
          if(_.contains(['select','multiselect','multiselect2','multiselect3'],fi.edit_type)){
            fieldSchema['options'] = self._createVocabularyChoices(fi);
            if (!fieldSchema.placeholder){
              fieldSchema.placeholder = 'choose ' + fieldSchema.title + '...';
            }
          }
          fieldSchema.validators = self._createValidators(fi);
          
          console.log('editSchema for key created: ', key, editSchema[key]);
          
        }

      });      
      
      console.log('editSchema created');
                 
      return editSchema;
    },

    _createValidators: function(fi) {
        
      var validators = [];
      var validator;
      
      if (_.contains(['integer','float','decimal'],fi.data_type))
      {
        if (fi.min && !_.isUndefined(fi.min)){
          validator = function checkMin(value, formValues) {
            var err = {
                type: 'Min',
                message: 'must be >= ' + fi.min
            };
            if (value < fi.min ) return err;
          };
          validators.unshift(validator);
        }
        if( !_.isUndefined(fi.range)){
          validator = function checkRange(value, formValues) {
            
            var last = '';
            var rangeMsg = '';
            var value_ok = false;
            var schema_lower = 0, schema_upper = 0, schema_val;
            for(var i=0; i<fi.range.length; i++){
              schema_val = fi.range[i]
              if(i>0) rangeMsg += ', ';
              if(i%2 === 0){
                rangeMsg += '> ' + schema_val
              }else{
                rangeMsg += '< ' + schema_val
              }
            }
            // compare range in pairs
            for(i=0; i<fi.range.length; i+=2){
              schema_lower = parseInt(fi.range[i],10)
              if(fi.range.length > i+1){
                schema_upper = parseInt(fi.range[i+1],10)
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
            if (!value_ok){
              return {
                  type: 'Range',
                  message: 'value not in range: ' + value + ' range: ' + rangeMsg
              };
            }
          };
          validators.unshift(validator);
        }
      }

      if(fi.required){
        validators.unshift('required');
      }
      
      if(!_.isUndefined(fi.regex) && !_.isEmpty(fi.regex)){
        validator = { type: 'regexp', regexp: new RegExp(fi.regex) };
        console.log('create RegExp: ', fi.regex, validator );
        if(!_.isUndefined(fi.validation_message) && !_.isEmpty(fi.validation_message)){
          validator.message = fi.validation_message;
          // TODO: rework, if req'd, to use tokenized strings (will need 
          // to reimplement backbone-forms
          //  function(value, formValues){
          //    //                TODO: figure out how to get the pending model
          //    return 'value: ' + value + ' is incorrect: ' + Iccbl.formatString(fi.validation_message, formValues);
          //  };
        }
        validators.unshift(validator);
      }
      return validators;
    },
    
    _createVocabularyChoices: function(fi){
      var choiceHash = fi.choices || [];
      if(!_.isEmpty(fi.vocabulary_scope_ref)){
        choiceHash = []
        // replace the fi.choices with the vocabulary, if available
        try{
          var vocabulary = appModel.getVocabulary(fi.vocabulary_scope_ref);
          _.each(_.keys(vocabulary),function(choice){
            if(vocabulary[choice].is_retired){
              console.log('skipping retired vocab: ',choice,vocabulary[choice].title );
            }else{
              choiceHash.push({ val: choice, label: vocabulary[choice].title });
            }
          });
          if(fi.edit_type == 'select' && !fi.required ){
            choiceHash.unshift({ val: '', label: ''});
          }
        }catch(e){
          var msg = 'Vocabulary unavailable: field: ' + fi.key +  
            ', vocabulary_scope_ref: ' + fi.vocabulary_scope_ref;
          console.log(msg,e);
          appModel.error(msg);
        }
      }
      return choiceHash;
    },
    
    /** 
     * Override the Backbone Forms templateData: this will take the place of
     * the serialize function, since we're overriding the 
     * Backbone Layoutmanager renderTemplate as well.
     */    
    templateData: function() {
      return {
        'fieldDefinitions': this.modelFields,
        'keys': _.chain(this.finalEditableKeys),
        'editKeys': _.chain(this.editKeys),
        'groupedKeys': _.chain(this.groupedKeys)
      };      
    },	

    /** 
     * Override the Backbone Layoutmanager template rendering to use 
     * Backbone Forms
     */
    renderTemplate: function() {
      console.log('renderTemplate');
      return Backbone.Form.prototype.render.apply(this);
    },

    /**
     * Override after render to patch in the multi select
     */
    afterRender: function(){
      var self = this;
      console.log('afterRender...');

      console.log('setup single selects using chosen...');
      // See http://harvesthq.github.io/chosen/
      this.$el.find('.chosen-select').chosen({
        disable_search_threshold: 3,
        width: '100%',
        allow_single_deselect: true,
        search_contains: true
        });
      
      _.each(this.editKeys, function(key){

        var fi = self.modelFields[key];
        if (fi.edit_type == 'typeahead' && fi.choices){
          $('[name="'+key +'"]').typeahead({
            autoSelect: false,
            delay: 1,
            minLength: 0,
            items: 'all',
            source: fi.choices
          });
        }
      });
      // TODO: move this to the multiselect2 render()
      console.log('setup multiselect2 elements using loudev multiselect');
      // see: http://loudev.com/
      // see: https://github.com/lou/multi-select
      // update the multiselect2 with the multiselect enhancements
      // multiselect with search using:
      // https://github.com/riklomas/quicksearch
      this.$el.find('.multiselect2').find('select').multiSelect({
        selectableOptgroup: true,
        selectableHeader: "<input type='text' class='form-control' autocomplete='off'>",
        selectionHeader: "<input type='text' class='form-control' autocomplete='off'>",
        afterInit: function(ms){
          var that = this,
              $selectableSearch = that.$selectableUl.prev(),
              $selectionSearch = that.$selectionUl.prev(),
              selectableSearchString = '#'+that.$container.attr('id')+' .ms-elem-selectable:not(.ms-selected)',
              selectionSearchString = '#'+that.$container.attr('id')+' .ms-elem-selection.ms-selected';

          that.qs1 = $selectableSearch.quicksearch(selectableSearchString).on('keydown', function(e){
            if (e.which === 40){
              that.$selectableUl.focus();
              return false;
            }
          });

          that.qs2 = $selectionSearch.quicksearch(selectionSearchString).on('keydown', function(e){
            if (e.which == 40){
              that.$selectionUl.focus();
              return false;
            }
          });
        },
        afterSelect: function(){
          this.qs1.cache();
          this.qs2.cache();
        },
        afterDeselect: function(){
          this.qs1.cache();
          this.qs2.cache();
        }
      });
      
      // Set autofocus (HTML5) on the first field
      // see: http://stackoverflow.com/questions/20457902/how-to-automatically-focus-first-backbone-forms-input-field
      this.$el.find('input').first().attr('autofocus','autofocus');
      
      // cache initial values to detect changes later:
      this.initialValues = this.getValue();
      console.log('Editview initialized, initial values:', this.initialValues);
      
      _.each(this.modelFields, function(field){
        if (field.required){
          var key = field.key;
          $('[name="'+key +'"').parents('.form-group').addClass('required');
        }
      });
      
      console.log('afterRender finished');
    },
    

    /**
     * Filter Backbone's changeAttributes:
     * - only include fields that were editable ('create' or 'update' edit_visibility)
     * - remove nulls/empty strings
     * - sort lists to test equality (TODO: support for ordered lists)
     */
    _getChangedAttributes: function(model){
      var self = this;
      var changedAttributes = model.changedAttributes();
      console.log('original changedAttributes: ', changedAttributes);
      if(changedAttributes){
        changedAttributes = _.pick(changedAttributes, this.finalEditableKeys);
      }
      changedAttributes = _.omit(changedAttributes, function(value,key,object){
        var prev;
        if( _.has(self.initialValues,key)){
          // test if the editor conversion is equivalent
          if (self.initialValues[key]== value){
            return true;
          }
        }
        
        // equate null to empty string, object, or array
        if(_.isNull(value)){
          prev = model.previous(key);
          if(_.isNull(prev)) return true;
          if(_.isObject(prev) || _.isString(prev) || _.isArray(prev)){
            return _.isEmpty(prev);
          }
        }
        if(_.isObject(value) || _.isString(value) || _.isArray(value)){
          if(_.isEmpty(value)){
            prev = model.previous(key);
            if(_.isNull(prev)) return true;
            if(_.isObject(prev) || _.isString(prev) || _.isArray(prev)){
              return _.isEmpty(prev);
            }
          }
        }
        
        // Cleanup strings containing newlines: 
        // carriage-return,line-feed (0x13,0x10) may be converted to line feed only (0x10)
        // NOTE: JSON does not officially support control-characters, so 
        // newlines should be escaped/unescaped on send/receive in the API (TODO)
        if(_.isString(value) && _.isString(model.previous(key))){
          if( value.replace(/(\r\n|\n|\r)/gm,"\n") == 
              model.previous(key).replace(/(\r\n|\n|\r)/gm,"\n") ){
            return true;
          }
        }
        return false;
      });
      
      return changedAttributes;
    },
    
    cancel: function(e) {
      e.preventDefault();
      this.remove();
      appModel.router.back();
    }, 
    
    save: function( event ) {
      event.preventDefault();
      var self = this;
      var errors, changedAttributes, url,
        options = {};
      var headers = options['headers'] = {};
      
      $('.has-error').removeClass('has-error');
      $('[data-error]').empty();
      errors = this.commit({ validate: true });
      if(errors){
        console.log('errors in form:', errors);
        _.each(_.keys(errors), function(key){
          var error = errors[key];
          if (_.has(self.fields, key)){
            $('[name="'+key +'"').parents('.form-group').addClass('has-error');
            console.log('added error for: "', key, '", val: "', self.fields[key].getValue(), '"');
          } else if (key=='_others') {
            var other_errors = errors[key];
            _.each(other_errors, function(error_obj){
              console.log('other error', error_obj);
              _.each(_.keys(error_obj), function(other_key){
                other_error = error_obj[other_key];
                if (_.has(self.fields, other_key)){
                  $('[key="form-group-'+other_key +'"')
                    .find('[data-error]').append('<br/>' + other_error);
                } else {
                  self.$el.append(
                  '<div data-error class="text-danger">' + other_key + ': ' + other_error + '</div>');
                }
              });
            });
          }
        });
        return;
      }
      
      changedAttributes = self._getChangedAttributes(this.model);
      if (! changedAttributes || _.isEmpty(changedAttributes)){
        appModel.error('no changes were detected');
        return;
      }

      // Set up options for Backbone sync / AJAX 
      
      // Wait for the server before setting the new attributes on the model
      options['wait'] = true;
      
      // Fixup the URL - if it points to the model instance, make it point to 
      // the API resource only: tastypie wants this
      // Note: this is happening if the model was fetched specifically for this
      // page, and has the url used to fetch it, rather than the collection url.
      url = options['url'] || _.result(this.model, 'url');
      // TODO: this should be optional (for most resources, to have url end with '/'
      if( url && url.charAt(url.length-1) != '/'){
        url += '/';
      }
      
      // Determine PATCH (update) or POST (create)
      options['key'] = Iccbl.getIdFromIdAttribute( self.model,self.model.resource );
      if (!_.contains(this.uriStack, '+add') && options['key'] ){
        self.model.idAttribute = self.model.resource['id_attribute'];
        options['patch'] = true;
        // TODO: check if creating new or updating here
        // set the id specifically on the model: backbone requires this to 
        // determine whether a "POST" or "PATCH" will be used
        this.model.id = options['key'];
      }
      
      headers[appModel.HEADER_APILOG_COMMENT] = self.model.get('apilog_comment');
      
      if(!_.isUndefined(this.saveCallBack) && _.isFunction(this.saveCallBack)){
        this.saveCallBack(this.model,headers,options, url);
      }else{
        // Note:
        // You have to send { dataType: 'text' } to have the success function 
        // work with jQuery and empty responses ( otherwise, fails on JSON.parse 
        // of the empty response).        
        //        dataType: 'text', 
        // The other solution: use "always_return_data" in the tastypie
        // resource definitions - which we are doing.
        console.log('save, changedAttributes: ', changedAttributes);
        this.model.save(changedAttributes, options)
          .success(function(model, resp) {
            console.log('success');
            if(!options['patch']){
              // this is an +add event
              model = new Backbone.Model(model);
              var key = Iccbl.getIdFromIdAttribute( model,self.model.resource );
              model.key = self.model.resource.key + '/' + key;
              appModel.router.navigate(self.model.resource.key + '/' + key, {trigger:true});
            }else{
              console.log('trigger remove');
              self.trigger('remove');
            }
          })
          .done(function(model, resp) {
            // TODO: done replaces success as of jq 1.8
            console.log('model saved');
          })
          .fail(function(jqXHR, textStatus, errorThrown){ 
            
            if (jqXHR && _.has(jqXHR,'responseJSON') && !_.isEmpty(jqXHR.responseJSON) ) {
              var errors = _.result(jqXHR.responseJSON,'errors',null);
              if(errors){
                console.log('errors in response:', errors);
                _.each(_.keys(errors), function(key){
                  var error = errors[key];
                  if (_.has(self.fields, key)){
                    self.fields[key].setError(error);
                  }
                  $('[name="'+key +'"').parents('.form-group').addClass('has-error');
                  console.log('added error for: "', key, '", val: "', self.fields[key].getValue(), '"');
                });
                return;
              }
            }
            if (options['patch']){
              self.model.fetch();
            }else{
              self.remove();
              appModel.router.back();
            }
            Iccbl.appModel.jqXHRfail.apply(this,arguments); 
            console.log('trigger remove');
            self.trigger('remove');
          })
          .always(function() {
            // always replaces complete as of jquery 1.8
          });
      }
      
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    }
    
	});
  EditView.ChosenSelect = ChosenSelect;
  EditView.DatePicker = DatePicker;
  EditView.DisabledField = DisabledField;
  EditView.SIunitEditor = SIunitEditor;
  /**
   * TODO: move this to utils module or parent view module
   */
  EditView.uploadAttachedFileDialog = function(url, attachedfileCollection, vocabulary_ref){
      var self = this;
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
      formSchema['filename'] = {
        title: 'Option 2: Name',
        key: 'filename',
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
            if (_.isEmpty(attrs.contents)){
              errs.contents = 'Specify either file or contents';
            }else{
              if (_.isEmpty(attrs.filename)){
                errs.filename = 'Must specify a filename with the file contents';
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
      form.$el.find("[name='fileInputPlaceholder']").append(
        '<label class="btn btn-default btn-file">' + 
        'Browse<input type="file" name="fileInput" style="display: none;"></label>'+ 
        '<p id="filechosen" class="form-control-static" ></p>');
      form.$el.on('change', ':file', function() {
        var input = $(this),
            numFiles = input.get(0).files ? input.get(0).files.length : 1,
            label = input.val().replace(/\\/g, '/').replace(/.*\//, '');
        form.$el.find('#filechosen').html(label);
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
              var filename;
              if(file){
                data.append('attached_file',file);
                filename = file.name;
                if (!_.isEmpty(values['filename'])){
                  filename = values['filename'];
                }
                if (_.isEmpty(values['file_date'])){
                  data.append(
                    'file_date',
                    Iccbl.getISODateString(_.result(file,'lastModifiedDate', null)));
                }
              }else{
                filename = values['filename'];
              }
              for(var pair of data.entries()) {
                 console.log(pair[0]+ ', '+ pair[1]); 
              }
              
              $.ajax({
                url: url,    
                data: data,
                cache: false,
                contentType: false,
                processData: false,
                type: 'POST',
                headers: headers, 
                success: function(data){
                  attachedfileCollection.fetch({ reset: true });
                  appModel.showModalMessage({
                    title: 'Attached File uploaded',
                    okText: 'ok',
                    body: '"' + filename + '"'
                  });
                },
                done: function(model, resp){
                  // TODO: done replaces success as of jq 1.8
                  console.log('done');
                }
              }).fail(function(){ appModel.jqXHRfail.apply(this,arguments); });
            
              return true;
            }
          },
          view: _form_el,
          title: 'Upload an Attached File'  });
      
    };
    
  
  
	return EditView;
});