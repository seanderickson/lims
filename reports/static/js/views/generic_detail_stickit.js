define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'backgrid',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_stickit',
    'views/simple-list',
    'text!templates/generic-detail-stickit.html'
], function( $, _, Backbone, stickit, BackGrid, Iccbl, layoutmanager, 
      appModel, DetailView, SimpleListView, detailTemplate) {
	
	var DetailView = Backbone.Layout.extend({
	
	  /**
	   * args:
	   * this.model - implicit as the first argument the constructor
	   * schema - resource schema hash/object
	   * schema.detailKeys() returns the metahash field keys 'detail' in 'visibility'
	   */
	  initialize: function(args) {
	    console.log('initialize generic_detail_stickit view');
	    var self = this;
	    var schema = this.schema = args.schema || this.model.resource.schema;
      this.detailKeys = args.detailKeys || schema.detailKeys(); 
      var nestedModels = this.nestedModels = {};
      var nestedLists = this.nestedLists = {};
      
      var buttons = this.buttons = args.buttons || ['download','history','back','edit','delete'];
      if(! appModel.hasPermission(self.model.resource.key, 'edit')){
        this.buttons = _.without(this.buttons,'edit');
      }
      if(! appModel.getCurrentUser().is_superuser){
        this.buttons = _.without(this.buttons,'delete');
      }
      this.createBindings();
	  },
	  
	  createBindings: function() {
	    var self = this;
	    var keys = this.detailKeys;
	    var schema = this.schema;
      var bindings = this.bindings = {};
      var schemaBindings = this.schemaBindings = {};
      
      _.each(keys, function(key) {
        bindings['#'+key] = self.createBinding(key,schema.fields[key]);
        schemaBindings['#title-'+key] = {
          observe: key,
          onGet: function(value) {
            if (value) return value.title + ":";
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
	    var vocabulary;
      var cell_options;
      var binding;
      var data_type = _.isEmpty(fi.data_type) ? 'string' : fi.data_type.toLowerCase();
      var display_type = _.isEmpty(fi.display_type) ? data_type : fi.display_type.toLowerCase();
      var data_type_formatters,display_type_formatters;
      
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

      if(!_.isEmpty(fi.vocabulary_scope_ref)){
        // replace the fi.choices with the vocabulary, if available
        try{
          vocabulary = Iccbl.appModel.getVocabulary(fi.vocabulary_scope_ref);
        }catch(e){
          var msg = 'Vocabulary unavailable: field: ' + fi.key +  
            ', vocabulary_scope_ref: ' + fi.vocabulary_scope_ref;
          console.log(msg,e);
          appModel.error(msg);
        }
      }
      if( fi.vocabulary ){
        vocabulary = fi.vocabulary;
      }
      
      function getTitle(vocabulary,value){
        if (!_.isEmpty(vocabulary[value])){
          if(vocabulary[value].title){
            return vocabulary[value].title;
          }else if(_.isString(vocabulary[value])){
            return _.escape(vocabulary[value]);//.replace(/</g,'&lt;');//.replace(/>/g,'&gt').replace(/&/g,'&amp');
          }else{
            console.log('error: ' + fi.vocabulary_scope_ref + ', key: ' + key, fi);
            appModel.error('vocabulary misconfigured for: ' + 
              fi.vocabulary_scope_ref + ', field: ' + fi.key + ': ' + value);
          }
        }else{
          console.log('error: ' + fi.vocabulary_scope_ref + ', key: ' + key, fi);
          appModel.error('vocabulary not found for: ' + 
            fi.vocabulary_scope_ref + ', field: ' + fi.key + ': ' + value);
        }
        return value;
      };      

      if(!_.isEmpty(fi['display_options'])){
        cell_options = fi['display_options'];
        cell_options = cell_options.replace(/'/g,'"');
        try{
          cell_options = JSON.parse(cell_options);
        }catch(e){
          console.log('warn: display_options is not JSON parseable, column: ',
              key,', options: ',cell_options);
        }
      }
  
      // define "data_type" getters
      function defaultGetter(value){
        if(value && vocabulary){
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
            return Iccbl.getUTCDateString(new Date(value));
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
          if(vocabulary){
            finalValue = Iccbl.sortOnOrdinal(value,vocabulary);
            finalValue = _.map(finalValue,function(value){ 
              return getTitle(vocabulary,value);
            });
          }
          finalValue = finalValue.join(', ');
        }
        return finalValue;
      }
      function decimalGetter(value){
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
      }
      updateMethod: 'html',
      function booleanGetter(value){
        if(_.isBoolean(value)){
          return value;
        }else{
          return false;
        }
      }
            
      data_type_formatters = {
        'date': dateGetter,
        'list': listGetter,
        'float': decimalGetter,
        'decimal': decimalGetter
        //'integer': defaultGetter,
        //'string' : defaultGetter,
        //'boolean' : defaultGetter
        //'datetime': defaultGetter,
      };
      
      // Define getters for "display_type"

      function siUnitGetter(value){
        var formatter = new Iccbl.SIUnitsFormatter(cell_options)
        console.log('get siunit value', value && !_.isEmpty(value), formatter.fromRaw(value));
        if(value){
          console.log('get siunit value', value, formatter.fromRaw(value));
          return formatter.fromRaw(value);
        }else{
          return value;
        }
      }

      function linkGetter(value){
        var finalValue = value;
        var _options = _.extend(
          { hrefTemplate: '#', target: '_self' }, cell_options );                
        // use the display options if needed for backward compatibility
        if( _.has(fi,'display_options') && _options.hrefTemplate == '#' ) {
          _options.hrefTemplate = window.location.pathname + '#' + fi['display_options'];
          _options.target = '_self';
        } 
        if(vocabulary){
          finalValue = getTitle(vocabulary,value);
        }
        var interpolatedVal = Iccbl.replaceTokens(self.model,_options.hrefTemplate,value);
 
        if(value && !_.isEmpty(value)){
          if(vocabulary){
            finalValue = getTitle(vocabulary,value);
          }
          var _html = '<a ' + 
            'id="link-' + key + '" ' + 
            'href="' + interpolatedVal + '" ' +
            'target="' + _options.target + '" ' +
            'tabIndex=-1 ' +
            '>' + finalValue + '</a>';
          return _html;
        }else{
          return value;
        }
      }
      
      function linkListGetter(values){  
        var modelValues = values;
        var finalValues = values;
        if(values){
          if(vocabulary){
            finalValues = Iccbl.sortOnOrdinal(modelValues,vocabulary);
          }
          return _.map(finalValues, linkGetter).join(', ');
        }
        return finalValues;
      }

      display_type_formatters = {
        'link': linkGetter,
        'linklist': linkListGetter,
        'siunit': siUnitGetter
      };
      
      // compose getter hierarchy; default<-data_type<-display_type<-vocabulary
      
      if(_.has(data_type_formatters, data_type)){
        binding.onGet = _.compose(binding.onGet, data_type_formatters[data_type]);
      }else{
        binding.onGet = _.compose(binding.onGet, defaultGetter);
      }
      if(_.has(display_type_formatters, display_type)){
        console.log('add getter for field: ', key, ', display_type: ', display_type);
        binding.onGet = _.compose(baseGetter, display_type_formatters[display_type]);
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
      console.log('generic detail stickit, afterRender');
      var self = this;
      this.stickit(this.model, this.bindings);
      this.schemaFieldsModel = new Backbone.Model(this.model.resource.schema.fields);
      this.stickit(this.schemaFieldsModel, this.schemaBindings);

      var btnbindings = {};
      var buttonModel = new Backbone.Model();
      _.each(this.buttons, function(button){
        btnbindings['[name="' + button + '"]'] = button;
        buttonModel.set(button,button);
      });
      this.stickit(buttonModel, btnbindings);
      
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
                  var id = Iccbl.getIdFromIdAttribute( model, model.resource.schema);
                  appModel.router.navigate(model.resource.key + '/' + id, {trigger:true});
                };
                return model;
              }
              ));
          var commentFields = _.pick(resource.schema.fields, ['username','date_time','comment']);
          var columns = Iccbl.createBackgridColModel(commentFields);
          var view = new Backgrid.Grid({
            columns: columns,
            collection: collection,
            schemaResult: resource.schema,
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
        'buttons': _.chain(this.buttons), // TODO: buttons from the schema
        'title': Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema),
        'keys': _.chain(this.detailKeys)
      };      
    },    
    
    template: _.template(detailTemplate)

	});

	return DetailView;
});