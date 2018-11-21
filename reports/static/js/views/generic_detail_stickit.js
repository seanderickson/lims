define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'backgrid',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'templates/generic-detail-stickit.html'
], function( $, _, Backbone, stickit, BackGrid, Iccbl, layoutmanager, 
      appModel, detailTemplate) {

  // NOTE: Webpack 3 patch:
  // Manually add the Stickit.ViewMixin object to the Backbone.Layout prototype
  // (the Stickit object is bound using the webpack assetsPluginInstance).
  _.extend(Backbone.Layout.prototype, Stickit.ViewMixin);
	
	var DetailView = Backbone.Layout.extend({


	  attributes: { id: 'generic-detail-stickit-container' },
    
	  /**
     * @param resource - the API resource schema
	   */
	  initialize: function(args) {
	    console.log('initialize generic_detail_stickit view');
	    var self = this;
	    this.args = args;
	    var resource = this.resource = args.resource || this.model.resource;
      
      var nestedModels = this.nestedModels = {};
      var nestedLists = this.nestedLists = {};
      var buttons = this.buttons = args.buttons || ['download','history','back','edit'];
      if (! appModel.isEditable(self.model.resource.key)
          || !appModel.hasPermission(resource.key, 'write')){
        
          this.buttons = _.without(this.buttons,'edit');
          this.buttons = _.without(this.buttons,'delete');
      }
      if (! appModel.hasPermission(self.model.resource.key, 'write')){
        this.buttons = _.without(this.buttons,'edit');
      }
      if(! appModel.getCurrentUser().is_superuser){
        this.buttons = _.without(this.buttons,'delete');
      }
      this.buttons_left = _.filter(this.buttons, function(button){
        return button == 'edit';
      });
      this.buttons_right = _.filter(this.buttons, function(button){
        return button != 'edit';
      });

      self.show_all_fields_control = $([
        '<label class="checkbox-inline" id="show_all_fields_control" ',
        '   title="Show all fields available for the resource" >',
        '  <input type="checkbox">show all available fields',
        '</label>'
        ].join(''));
      
      this.createBindings();
	  },
	  
	  get_show_extra_fields: function() {
	    return this.show_all_fields_control.find('input[type="checkbox"]').prop('checked');
	  },
	  
	  createBindings: function() {
	    var self = this;
	    var resource = this.resource;
	    var bindings = this.bindings = {};
	    var schemaBindings = this.schemaBindings = {};

      var detailKeys = self.args.detailKeys || resource.detailKeys(); 

      if (self.get_show_extra_fields()){
        detailKeys = self.resource.allKeys();
      }
      
      var adminKeys = this.adminKeys = self.model.resource.adminKeys();
      if (! appModel.hasGroup('readEverythingAdmin')) {
        detailKeys = _.difference(detailKeys, adminKeys);
      }
	    
	    if (!this.get_show_extra_fields()){
        // If "hideIfEmpty" then remove null attributes
        _.each(self.model.keys(), function(key){
          if(! self.model.has(key) && _.has(resource.fields,key)){
            var fi = resource.fields[key];
            if (fi.display_options && fi.display_options.hideIfEmpty === true){
              detailKeys = _.without(detailKeys, key);
            }
          }
        });
	    }
	    self.detailKeys = detailKeys;
      var groupedKeys = self.groupedKeys = resource.groupedKeys(detailKeys);

      if (appModel.DEBUG){
        console.log('final detailKeys', self.detailKeys, self.groupedKeys);
      }
      
      function create_title(fi){
        if (appModel.getCurrentUser().is_superuser
            && fi.vocabulary_scope_ref){
          // Show link to the vocab term for superuser
          return Iccbl.formatString(
            '<a href="#vocabulary/search/scope__exact={vocabulary_scope_ref}" '
            + ' class="" '
            + ' target=_blank >{title}</a>',fi);
        } else {
          return fi.title + ':';
        }
      };

      _.each(self.detailKeys, function(key) {
        bindings['#'+key] = self.createBinding(key,resource.fields[key]);
        schemaBindings['#title-'+key] = {
          observe: key,
          updateMethod: 'html',
          onGet: function(value) {
            if (value){
              return create_title(value);
            }
            else return 'Title for ' + key;
          },
          attributes: [{
            name: 'title', observe: key,
            onGet: function(value) {
              if (value) return value.description;
            }
          }]
        };
      });      
	  },
	  
	  createBinding: function(key, fi){
	    var self = this;
      var binding;
      var data_type = _.isEmpty(fi.data_type) ? 'string' : fi.data_type.toLowerCase();
      var display_type = _.isEmpty(fi.display_type) ? data_type : fi.display_type.toLowerCase();
      var data_type_formatters,display_type_formatters;
      var cell_options = fi.display_options;
      
      function baseGetter(value) {
        if (_.isUndefined(value) || _.isNull(value)) return '-';
        else if (_.isString(value) && value === '' ) return '-';
        else if (_.isObject(value) && _.isEmpty(value)) return '-';
        else{
          return value;
        }
      }     
	    binding = {
          observe: key,
          updateMethod: 'html',
          onGet: baseGetter
        };
      
      if(display_type == 'link'){
        binding.updateMethod = 'html';
        if( data_type == 'list') display_type = 'linklist';
      }

      function getVocabulary(){
        if( fi.vocabulary ){
          return fi.vocabulary;
        }else if(!_.isEmpty(fi.vocabulary_scope_ref)){
          // replace the fi.choices with the vocabulary, if available
          try{
            return Iccbl.appModel.getVocabulary(fi.vocabulary_scope_ref);
          }catch(e){
            var msg = 'Vocabulary unavailable: field: ' + key +  
              ', vocabulary_scope_ref: ' + fi.vocabulary_scope_ref;
            console.log(msg,e);
            appModel.error(msg);
          }
        }
      }
      
      function getTitle(vocabulary,value){
        if(_.isUndefined(value) || _.isNull(value)) return value;
        if (!_.isEmpty(vocabulary[value])){
          if(vocabulary[value].title){
            return vocabulary[value].title;
          }else if(_.isString(vocabulary[value])){
            return _.escape(vocabulary[value]);
            //.replace(/</g,'&lt;');//.replace(/>/g,'&gt').replace(/&/g,'&amp');
          }else{
            console.log('error: ' + fi.vocabulary_scope_ref + ', key: ' + key, fi);
            appModel.error('vocabulary misconfigured for: ' + 
              fi.vocabulary_scope_ref + ', field: ' + fi.key + ': ' + value);
          }
        }
        return value;
      };      
      
      // define "data_type" getters
      function defaultGetter(value){
        var vocabulary = getVocabulary();
        if(!_.isUndefined(value) && !_.isNull(value) && vocabulary){
          value = getTitle(vocabulary,value);
        }
        if(value && _.isString(value)){
          value = value.replace(/(\r\n|\n|\r)/gm,"<br/>")
        }
        return value;
      }
      
      function dateGetter(value){
        if (value && !_.isEmpty(value)) {
          try {
            return Iccbl.getDateString(value);
          } catch(e) {
            var msg = Iccbl.formatString(
              'unable to parse date value: {value}, for field: {field}',
              { field:key, value: value } );
            console.log(msg,e);
            appModel.error(msg);
          }
        }
        return value;
      }

      function listGetter(value){
        var finalValue = value;
        if(_.isArray(value)){
          var vocabulary = getVocabulary();
          if(vocabulary){
            finalValue = _.map(value,function(v){ 
              return getTitle(vocabulary,v);
            });
          }
          finalValue = finalValue.join(', ');
        }
        return finalValue;
      }
      function decimalGetter(value){
        if (cell_options){
          console.log('decimalGetter:', value);
          var formatted;
          var formatter = new Iccbl.DecimalFormatter(cell_options);
          if(_.isString(value) || _.isNumber(value)){
            formatted = formatter.fromRaw(value);
            console.log('decimal getter: v:', value, ', formatted: ', 
              formatted, ', options: ', cell_options);
            return formatted;
          }else{
            return value;
          }
        } else {
          return value;
        }
      }
      function integerGetter(value){
        var vocabulary = getVocabulary();
        if(!_.isUndefined(value) && !_.isNull(value) && vocabulary){
          value = getTitle(vocabulary,value);
          return value;
        }
        else if(_.isNumber(value)){
          var formatter = new Iccbl.IntegerFormatter(cell_options);
          formatted = formatter.fromRaw(value);
          if (appModel.DEBUG){
            console.log('integer getter: v:', value, ', formatted: ', 
              formatted, ', options: ', cell_options);
          }
          return formatted;
        }else{
          return value;
        }
      }
      function booleanGetter(value){
        if(_.isBoolean(value)){
          if (value === true) return 'True';
          else return 'False';
        }else{
          return value;
        }
      };
          
      function imageGetter(value){
        if (!_.isEmpty(value) && value != '-' ){
          console.log('image:', value);
          return '<img src="' + value + '" alt="structure image" />';
        }
      }
      
      data_type_formatters = {
        'date': dateGetter,
        'list': listGetter,
        'float': decimalGetter,
        'decimal': decimalGetter,
        'integer': integerGetter,
        //'string' : defaultGetter,
        'boolean' : booleanGetter
        //'datetime': defaultGetter,
      };
      
      // Define getters for "display_type"

      function siUnitGetter(value){
        var formatter = new Iccbl.SIUnitsFormatter(cell_options)
        if(value){
          return formatter.fromRaw(value);
        }else{
          return value;
        }
      };

      function linkGetter(value){
        var _options = _.extend(
          { hrefTemplate: '#', target: '_self' }, cell_options );                
        // use the display options if needed for backward compatibility
        if( _.has(fi,'display_options') && _options.hrefTemplate == '#' ) {
          _options.hrefTemplate = window.location.pathname + '#' + fi['display_options'];
          _options.target = '_self';
        } 

        if(value && !_.isNull(value) && value != '-' ){
          var interpolatedVal = Iccbl.formatString(_options.hrefTemplate,self.model, value);
          var _html = '<a ' + 
            'id="' + key + '" ' + 
            'href="' + interpolatedVal + '" ' +
            'target="' + _options.target + '" ' +
            'tabIndex=-1 ' +
            '>' + value + '</a>';
          return _html;
        }else{
          return value;
        }
      };
      
      function linkListGetter(values){  
        if(values){
          var modelValues = values;
          var _options = _.extend(
            { hrefTemplate: '#', target: '_self' }, cell_options );                
          // use the display options if needed for backward compatibility
          if( _.has(fi,'display_options') && _options.hrefTemplate == '#' ) {
            _options.hrefTemplate = window.location.pathname + '#' + fi['display_options'];
            _options.target = '_self';
          } 
          var vocabulary = getVocabulary();
          var output = [];
          _.each(modelValues, function(value){
            var text = _.result(vocabulary, value, value);
            if(value && !_.isNull(value) && value != '-' ){
              var interpolatedVal = Iccbl.formatString(
                _options.hrefTemplate,self.model, value);
              var _html = '<a ' + 
                'id="' + key + '" ' + 
                'href="' + interpolatedVal + '" ' +
                'target="' + _options.target + '" ' +
                'tabIndex=-1 ' +
                '>' + text + '</a>';
              output.push(_html);
            }else{
              output.push(value);
            }
          });
          var sep = ', ';
          if (fi.display_options && fi.display_options.separator ){
            sep = fi.display_options.separator;
          }
          return output.join(sep);
        }
        return values;
      };

      display_type_formatters = {
        'link': linkGetter,
        'linklist': linkListGetter,
        'siunit': siUnitGetter,
        'image': imageGetter
      };
      
      // compose getter hierarchy; default<-data_type<-display_type<-vocabulary
      
      if(_.has(data_type_formatters, data_type)){
        binding.onGet = _.compose(binding.onGet,data_type_formatters[data_type]);
      }else{
        binding.onGet = _.compose(binding.onGet, defaultGetter);
      }
      if(_.has(display_type_formatters, display_type)){
        if (appModel.DEBUG) 
          console.log('add getter for field: ', key, ', display_type: ', display_type);
        if (display_type == 'linklist'){
          // for the linklist, have to hack in the defaultGetter
          binding.onGet = display_type_formatters[display_type];
        }else if (display_type == 'siunit'){
          // for siunit, skip other (i.e. decimalGetter)
          binding.onGet = _.compose(baseGetter, display_type_formatters[display_type]);
        }else if (display_type == 'link'){
          // for linkGetter, run the other (i.e. dateGetter) first
          binding.onGet = _.compose(display_type_formatters[display_type], binding.onGet);
        }else{
          console.log('unknown display type', display_type, key);
          binding.onGet = _.compose(display_type_formatters[display_type], binding.onGet);
        }
      }

      // FIXME: "nested" and "nested_list" are unused 20150629 after ui_type refactor to data_type
      
      if(display_type === 'nested'){
        nestedModels[key] = new Backbone.Model(self.model.get(key));
        nestedModels[key].resourceId = fi.nested_resource || key;
      }else if(display_type === 'nested_list'){
        nestedLists[key] = {};
        nestedLists[key].list = self.model.get(key);
        nestedLists[key].resourceId = fi.nested_resource || key;
        if(!nestedLists[key].list || nestedLists[key].list.length === 0){
          delete nestedLists[key];
          self.model.unset(key); // to signal empty
        }else if (nestedLists[key].length==1) {
          nestedModels[key] = nestedLists[key][0];
          nestedModels[key].resourceId = nestedLists[key].resourceId;
          delete nestedLists[key];
          self.model.unset(key); // to signal empty
        }
      }
      
      return binding;
	  },
	  
    afterRender : function() {
      var self = this;
      this.stickit(this.model, this.bindings);
      this.schemaFieldsModel = new Backbone.Model(this.model.resource.fields);
      this.stickit(this.schemaFieldsModel, this.schemaBindings);

      var btnbindings = {};
      // FIXME: 20171114 - using stickit to manage the buttons is unnecessary
      // Better would be using backbone forms
      var buttonModel = new Backbone.Model();
      _.each(this.buttons, function(button){
        btnbindings['[name="' + button + '"]'] = button;
        var buttonName = button.charAt(0).toUpperCase() + button.slice(1);
        buttonModel.set(button,buttonName);
      });
      this.stickit(buttonModel, btnbindings);
      
      if (appModel.getCurrentUser().is_superuser ) {
        $('#generic-detail-bottom-buttonpanel-left').prepend(self.show_all_fields_control);
        self.show_all_fields_control.find('input[type="checkbox"]').change(function(e) {
          console.log('click show all', e);
          e.preventDefault();
          self.createBindings();
          self.render();
        });
      }
      
      if (!_.isEmpty(self.nestedModels)) {
        _.each(_.keys(self.nestedModels), function(key) {
          var nestedModel = self.nestedModels[key];
          nestedModel.resource = appModel.getResource(nestedModel.resourceId);
          nestedModel.key = null; // TODO
         
          var view = new DetailView({ model: nestedModel, buttons: [] });
          view = self.setView("#"+key, view ).render();
          view.$el.addClass('well');
          view.$el.addClass('nested');
        });
      }
      
      // redo
      if (!_.isEmpty(self.nestedLists)) {
        _.each(_.keys(self.nestedLists), function(key) {
          var resource = appModel.getResource(self.nestedLists[key].resourceId)
          var collection;
          collection = new Backbone.Collection(
              _.map(self.nestedLists[key].list, function(item) {
                model = new Backbone.Model(item);
                model.resource = resource;
                model.collection = collection;
                // custom handle click on cell
                model.clickHandler = function(model){
                  var id = Iccbl.getIdFromIdAttribute( model, model.resource);
                  appModel.router.navigate(model.resource.key + '/' + id, {trigger:true});
                };
                return model;
              }
              ));
          var commentFields = _.pick(resource.fields, ['username','date_time','comment']);
          var columns = Iccbl.createBackgridColModel(commentFields);
          var view = new Backgrid.Grid({
            columns: columns,
            collection: collection,
            resource: resource
          });
          // FIXME: this should work
          //            Backbone.Layout.setupView(view);
          //            view = self.setView("#"+key, view ).render();
          // FIXME: memory leak if left this way
          self.$('#'+key).html(view.render().$el);
          view.$el.addClass('well');
          view.$el.addClass('nested');
        });
      }

      console.log('generic detail stickit, afterRender done');
      return this;
    },

    serialize: function() {
      return {
        'buttons_right': _.chain(this.buttons_right), // TODO: buttons from the resource
        'buttons_left': _.chain(this.buttons_left), // TODO: buttons from the resource
        'groupedKeys': _.chain(this.groupedKeys),
        'keys': _.chain(this.detailKeys), // TODO: groupedKeys replaces detailKeys
        'adminKeys': this.adminKeys,
        'table_class': _.result(this.args,'table_class', 'col-sm-8'),
        'label_col_class': _.result(this.args,'label_col_class', 'col-xs-2'),
        'value_col_class': _.result(this.args,'value_col_class', 'col-xs-8')    
      };      
    },    
    
    template: _.template(detailTemplate)
    
	});

	return DetailView;
});