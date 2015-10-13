define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_forms',
    'multiselect',
    'quicksearch',
    'iccbl_backgrid',
    'models/app_state',
    'text!templates/generic-edit.html',
    'text!templates/modal_ok_cancel.html',
    'bootstrap',
    'bootstrap-datepicker',
    'chosen'
], function( $, _, Backbone, backbone_forms, multiselect, quicksearch, Iccbl, appModel,
            editTemplate, modalOkCancel ) {
  
  var SIunitEditor = Backbone.Form.editors.Base.extend({
    
    tagname: 'siuniteditor',
    
    fieldTemplate: _.template([
      '<div data-editor title="<%= help %>"  >',
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
      ['m', 1e-3,],
      ['μ', 1e-6,],
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
      console.log('siunit initialize');
      Backbone.Form.editors.Base.prototype.initialize.call(this, options);

      var options = this.options = options || {};
      var formSchema = this.formSchema = {};
      var multiplier = this.multiplier = options.schema.multiplier || 1;
      var symbol = this.symbol = options.schema.symbol;
      var decimals = this.decimals = options.schema.decimals || 0;
      var defaultUnit = this.defaultUnit = options.schema.defaultUnit || 1;
      var units = this.units = [];
      
      if(! symbol){
        throw 'Error: SIUnitFormFilter requires a "symbol" option' 
      }
      _.each(this.siunits,function(pair){
        units.push({ val: pair[1], label: pair[0] + self.symbol });
      });
      formSchema['unit'] = {
        title: '', 
        key:  'unit', 
        type: 'Select',
        options: units,
        editorClass: 'form-control',
        template: this.fieldTemplate
      };
      formSchema['number'] = {
        title: '', 
        key:  'number',
        type: 'Number',
        editorClass: 'form-control',
        template: this.fieldTemplate
      };
      
      console.log('siunit initialize complete');
    
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
      console.log('render siunit', this.value);
  
      this._observeFormEvents();

      this.nestedForm = new Backbone.Form({
        schema: this.formSchema,
        model: new Backbone.Model(this._findNumberAndUnit(this.value)),
        idPrefix: this.id + '_',
        Field: Backbone.Form.NestedField,
        template: this.formTemplate
      });

      this.$el.html(this.nestedForm.render().el);
      
      if (this.hasFocus) this.trigger('blur', this);

      console.log('render siunit done', this.nestedForm.getValue());
      return this;
    },
    
    
    /**
     * Returns the current editor value
     * @return {String}
     */
    getValue: function() {
      console.log('getvalue', this.nestedForm.getValue());
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
      console.log('calculate',multiplier,si_mult,val);
      // run strip after every calculation to round out floating point math errors
      function strip(number) {
        return (parseFloat(number.toPrecision(12)));
        };
      val = strip(val * multiplier);
      if(si_mult > 0){ // if si unit is undefined, assume to be 1
        val = strip(val * si_mult);
      }
      console.log('calculate',multiplier,si_mult,val);
      return val;
    },
    
    _findNumberAndUnit: function(number){
      console.log('findNumberAndUnit',number);
      var self = this;
      function strip(number) {
        return (parseFloat(number.toPrecision(12)));
        };
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
      
      console.log('findNumberAndUnit',val,pair);
      return {number:val, unit: pair[1]};
    },
  
    remove: function() {
      console.log('remove');
      this.nestedForm.remove();
  
      Backbone.View.prototype.remove.call(this);
    },
  
    validate: function() {
      console.log('v');
      return this.nestedForm.validate();
    },
  
    _observeFormEvents: function() {
      console.log('obs');
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
        if (!this.value) {
          var date = new Date();
          date.setSeconds(0);
          date.setMilliseconds(0);
          this.value = date;
        }else{
          this.value = new Date(this.value)
        }
    },

    getValue: function() {
      var input = $('input', this.el),
      date = input.datepicker('getDate');
      return Iccbl.getISODateString(date);
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
         '<div class="input-group date input-group-sm" >',
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
          todayBtn: true,
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

  var DisabledField = Backbone.Form.editors.Base.extend({

    tagName: 'p',

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
    },

    getValue: function() {
    },

    setValue: function(value) {
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
      // FIXME: delegate disabled rendering to generic_detail render
      var el = $(this.el);
      var val = this.value;
      if(_.isArray(this.value)){
        val = this.value.join(', ');
      }
      el.html(val);
      return this;
    }
    
  });  

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
  
  // FIXME: 20150624 - overriding the checkbox editor, to use "div" instead of "<ul><li>
  // - something in CSS/JS is overriding the click function and not allowing 
  // selection of checkboxes in the <ul>,<li> elements
  Backbone.Form.editors.Checkboxes = Backbone.Form.editors.Checkboxes.extend({

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
  
  
  Backbone.Form.editors.MultiSelect2 = Backbone.Form.editors.Select.extend({
    render: function() {
      this.$el.attr('multiple', 'multiple');
      this.setOptions(this.schema.options);
 
      return this;
    }
  });
  
  
  
  var EditView = Backbone.Form.extend({
    
    initialize: function(args) {
      console.log('Editview initialize: ', args);
      this.uriStack = args.uriStack;
      this.consumedStack = []; 
      this.saveCallBack = args.saveCallBack;
      
      // NOTE: due to a phantomjs js bug, must convert arguments to a real array
      Backbone.Form.prototype.initialize.apply(this,Array.prototype.slice.apply(arguments));
    },

    events: {
      'click button#save': 'save'
    },
    
    altFieldTemplate:  _.template([
      '<div class="form-group" >',
      '    <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>',
      '    <div class="col-sm-2" >',
      '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
      '      <div data-error class="text-danger" ></div>',
      '      <div><%= help %></div>',
      '    </div>',
      '  </div>',
    ].join('')),

    datepickerComponentTemplate:  _.template([
      '<div class="form-group" >',
      '    <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>',
      '    <div class="col-sm-2" >',
      '      <div data-editor style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
      '      <div data-error class="text-danger" ></div>',
      '      <div><%= help %></div>',
      '    </div>',
      '</div>',
    ].join('')),
    
    altMultiselect2FieldTemplate:  _.template('\
      <div class="form-group" > \
          <label class="control-label col-sm-2" for="<%= editorId %>"><%= title %></label>\
          <div class="col-sm-10" >\
            <div data-editor  class="multiselect2" style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />\
            <div data-error class="text-danger" ></div>\
            <div><%= help %></div>\
          </div> \
        </div>\
    '),

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
      console.log('schema...', this.model);
      var self = this;
      var schema = this.model.resource.schema;
      var keys = Iccbl.sortOnOrdinal(
          _.keys(schema.fields), schema.fields)
      var editKeys = _(keys).filter(function(key){
          return _.has(schema.fields, key) &&
              _.has(schema.fields[key], 'visibility') &&
              _.contains(schema.fields[key]['visibility'], 'edit');
      });
      
      var editSchema = {};
      var itemcount = 0;
      var typeMap = {
        'boolean': 'Checkbox',
        'string': 'Text',
        'uri': 'Text',
        'float': 'Number',
        'integer': 'Number',
        'decimal': 'Number',
        'list': 'Checkboxes', 
        'date': DatePicker
      };
      
      // process the data_types - convert to backbone-forms schema editor type
      _.each(keys, function(key){
        console.log('build edit schema for key: ',key);
        
        var validators = [];
        var fieldSchema = editSchema[key] = {};
        var fi = schema.fields[key];
        var cell_options;
        
        if(!_.isEmpty(fi['display_options'])){
          cell_options = fi['display_options'];
          cell_options = cell_options.replace(/'/g,'"');
          try{
            cell_options = JSON.parse(cell_options);
            _.extend(fieldSchema,cell_options);
          }catch(e){
            console.log('warn: display_options is not JSON parseable, column: ',
                key,', options: ',cell_options);
          }
        }
        
        fieldSchema['title'] = fi.title
        fieldSchema['template'] = self.altFieldTemplate;
        
        var data_type = fi.data_type || 'Text';
        data_type = data_type.toLowerCase();
        fieldSchema['type'] = data_type.charAt(0).toUpperCase() + data_type.slice(1);
        
        if(!_.contains(editKeys, key)){
          console.log('create disabled entry', key, fi['visibility'])
          fieldSchema['disabled'] = true;
          fieldSchema['type'] = DisabledField;
          fieldSchema['editorClass'] = 'form-control-static';
          return;
        }
        
        if(_.has(typeMap, fi.data_type)){
          fieldSchema['type'] = typeMap[fi.data_type];
        }
        
        if(fi.data_type == 'string' || fi.data_type == 'float' 
          || fi.data_type == 'decimal' || fi.data_type == 'integer' ){
          fieldSchema['editorClass'] = 'form-control';
        }
        
        if(fi.data_type == 'date'){
          fieldSchema['template'] = self.datepickerComponentTemplate;
        } 
        
        if(fi.display_type == 'siunit'){
          console.log('siunit field: ' + key);
          fieldSchema['type'] = SIunitEditor;
          fieldSchema['editorClass'] = '';
        }
        
        if(fi.edit_type == 'select'){
          fieldSchema['type'] = 'Select';
          fieldSchema['editorClass'] = 'chosen-select';
        } 
        else if(fi.edit_type == 'multiselect'){
          fieldSchema['type'] = 'Checkboxes';
        } 
        else if(fi.edit_type == 'multiselect2'){
          fieldSchema['type'] = 'MultiSelect2';
          fieldSchema['template'] = self.altMultiselect2FieldTemplate;
        } 
        
        if(fi.edit_type == 'select' 
          || fi.edit_type == 'multiselect'
          || fi.edit_type == 'multiselect2'){

          fieldSchema['options'] = fi.choices || [];
          if(!_.isEmpty(fi.vocabulary_scope_ref)){
            // replace the fi.choices with the vocabulary, if available
            try{
              var vocabulary = appModel.getVocabulary(fi.vocabulary_scope_ref);
              var choiceHash = {}
              _.each(_.keys(vocabulary),function(choice){
                choiceHash[choice] = vocabulary[choice].title;
              });
              if(fi.edit_type == 'select'){
                choiceHash[''] = '';
              }
              fieldSchema['options'] = choiceHash;
            }catch(e){
              var msg = 'Vocabulary unavailable: field: ' + fi.key +  
                ', vocabulary_scope_ref: ' + fi.vocabulary_scope_ref;
              console.log(msg,e);
              appModel.error(msg);
            }
          }
        }
        
        if(fi.edit_type == 'radio'){
          fieldSchema['template'] = self.altRadioFieldTemplate;
        }
        
        //// validation stuff ////

        if(fieldSchema['type']  == 'Number')
        {
          // TODO: check for the "min" "max","range" validation properties and implement
          if( fi.min && !_.isUndefined(fi.min)){
            var validator = function checkMin(value, formValues) {
              var err = {
                  type: 'Min',
                  message: 'must be >= ' + fi.min
              };

              if (value <= fi.min ) return err;
            };
            validators.unshift(validator);
          }
          if( !_.isUndefined(fi.range)){
            var validator = function checkRange(value, formValues) {
              
              var last = '';
              var rangeMsg = '';
              var value_ok = false;
              var schema_lower = 0, schema_upper = 0
              for(var i=0; i<fi.range.length; i++){
                var schema_val = fi.range[i]
                if(i>0) rangeMsg += ', ';
                if(i%2 == 0){
                  rangeMsg += '> ' + schema_val
                }else{
                  rangeMsg += '< ' + schema_val
                }
              }
              // compare range in pairs
              for(var i=0; i<fi.range.length; i+=2){
                schema_lower = parseInt(fi.range[i])
                if(fi.range.length > i+1){
                  schema_upper = parseInt(fi.range[i+1])
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
        if(fi.required){
          validators.unshift('required');
        }
        if(!_.isUndefined(fi.regex) && !_.isEmpty(fi.regex)){
          var validator = { type: 'regexp', regexp: new RegExp(fi.regex) };
          if(!_.isUndefined(fi.validation_message) && !_.isEmpty(fi.validation_message)){
            validator.message = fi.validation_message;
            // TODO: rework, if req'd, to use tokenized strings (will need 
            // to reimplement backbone-forms
            //  function(value, formValues){
            //    //                TODO: figure out how to get the pending model
            //    return 'value: ' + value + ' is incorrect: ' + Iccbl.replaceTokens(new Backbone.Model(formValues), fi.validation_message);
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
        console.log('editSchema for key created: ', key, editSchema[key]);
      });      
      
      console.log('editSchema created');
      
      if( ! _.has(editSchema, 'apilog_comment')){
        console.log('enforced apilog_comment');
        // Note: Enforced comment
        editSchema['apilog_comment'] = {
            type: 'TextArea',
            title: 'Changelog Comment',
            validators: ['required'], 
            template: self.altFieldTemplate
        };
      }
           
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
      
      if( ! _.contains(editKeys, 'apilog_comment')){
        editKeys.push('apilog_comment');
      }
                  
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

    /**
     * Override after render to patch in the multi select
     */
    afterRender: function(){
      console.log('afterRender...');
      
      this.$el.find('.chosen-select').chosen({
        disable_search_threshold: 3,
        width: '100%',
        allow_single_deselect: true,
        search_contains: true
        });
      
      
      // update the multiselect2 with the multiselect enhancements
      // multiselect with search
      this.$el.find('.multiselect2').find('select').multiSelect({
        selectableOptgroup: true,
        selectableHeader: "<input type='text' class='search-input' autocomplete='off'>",
        selectionHeader: "<input type='text' class='search-input' autocomplete='off'>",
        afterInit: function(ms){
          var that = this,
              $selectableSearch = that.$selectableUl.prev(),
              $selectionSearch = that.$selectionUl.prev(),
              selectableSearchString = '#'+that.$container.attr('id')+' .ms-elem-selectable:not(.ms-selected)',
              selectionSearchString = '#'+that.$container.attr('id')+' .ms-elem-selection.ms-selected';

          that.qs1 = $selectableSearch.quicksearch(selectableSearchString)
          .on('keydown', function(e){
            if (e.which === 40){
              that.$selectableUl.focus();
              return false;
            }
          });

          that.qs2 = $selectionSearch.quicksearch(selectionSearchString)
          .on('keydown', function(e){
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
      
      var headers = {};
      headers[appModel.HEADER_APILOG_COMMENT] = self.model.get('apilog_comment');
      
      if(!_.isUndefined(this.saveCallBack) && _.isFunction(this.saveCallBack)){
        this.saveCallBack(this.model,headers,_patch,url);
      }else{
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
          headers: headers,
          wait: true // wait for the server before setting the new attributes on the model
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
        .error(function(model,response,options){
          // TODO: investigate: wait:true does not work if the model was already updated
          self.model.set(self.model.previousAttributes());
          appModel.backboneFetchError(model,response,options);
        })
        .always(function() {
          // always replaces complete as of jquery 1.8
          self.trigger('remove');
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
  EditView.DatePicker = DatePicker;
  
	return EditView;
});