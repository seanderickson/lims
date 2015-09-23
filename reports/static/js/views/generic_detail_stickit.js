define([
    'jquery',
    'underscore',
    'backbone',
    'backgrid',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'views/generic_detail_stickit',
    'views/simple-list',
    'text!templates/generic-detail-stickit.html'

], function( $, _, Backbone, BackGrid, Iccbl, layoutmanager, 
      appModel, DetailView, SimpleListView, detailTemplate) {
	
	var DetailView = Backbone.Layout.extend({
	
	  /**
	   * args:
	   * this.model - implicit as the first argument the constructor
	   * schema - resource schema hash/object
	   * schema.detailKeys() returns the metahash field keys 'detail' in 'visibility'
	   */
	  initialize: function(args) {
	    var self = this;
      var schema = args.schema || this.model.resource.schema;
      
      //      var buttons = this.buttons = 
      //          args.buttons || ['edit','delete','download','history','back'];
      var buttons = this.buttons = args.buttons || ['download','history','back'];
      
      if(appModel.hasPermission(self.model.resource.key, 'edit')){
        this.buttons.unshift('edit');
      }

      if(appModel.getCurrentUser().is_superuser){
        this.buttons.unshift('delete');
      }
      
      var keys = this.detailKeys = args.detailKeys || schema.detailKeys(); //this.model);
      var events = this.events = {};
      var bindings = this.bindings = {};
      var schemaBindings = this.schemaBindings = {};
      var nestedModels = this.nestedModels = {};
      var nestedLists = this.nestedLists = {};
      _.each(keys, function(key) {
        var fi = schema.fields[key];
        var data_type = _.isEmpty(fi.data_type)?'string':fi.data_type.toLowerCase();
        var display_type = _.isEmpty(fi.display_type)?data_type:fi.display_type.toLowerCase();
        var edit_type = _.isEmpty(fi.edit_type)?display_type:fi.edit_type.toLowerCase();
        
        // Default binding
        bindings['#'+key] = {
          observe: key,
          onGet: function(value) {
            if (_.isUndefined(value) || _.isNull(value)) return '-';
            else if (_.isString(value) && value === '' ) return '-';
            else if (_.isObject(value) && _.isEmpty(value)) return '-';
            else return value;
          }
        };
        
        // cell type specific render:
        // display_type has precedence
        var cell_options;
        if(!_.isEmpty(fi['display_options'])){
          cell_options = fi['display_options'];
          cell_options = cell_options.replace(/'/g,'"');
          try{
            cell_options = JSON.parse(cell_options);
          }catch(e){
            console.log('warn: display_options is not JSON parseable, column: ' +
                key + ', options: ' + cell_options);
          }
        }
        
        if(data_type == 'list'){
          var temp = bindings['#'+key].onGet;
          bindings['#'+key].onGet = function(value){
            var temp1 = temp(value);
            if(_.isArray(temp1)){
              return temp1.join(', ');
            }
            return temp1;
          };
        }
        
        if(display_type == 'date'){
          bindings['#'+key].onGet = function(value){
            console.log('date value', value);
            if(!value || _.isUndefined(value)) return '-';
            try{
              var date = new Date(value);
              var month = (date.getMonth()+1);
              if(month < 10) month = '0' + month;
              var day = date.getDate();
              if(day < 10) day = '0' + day;
              return date.getFullYear() 
                + '-' + month
                + '-' + day;
            }catch(e){
              console.log('warn: unable to parse date value: ' + key + ', ' + value );
            }
            return value;
          };
          
        }else if(display_type == 'siunit'){
          var formatter = new Iccbl.SciUnitsFormatter(cell_options)
          bindings['#'+key].onGet = function(value){
            return formatter.fromRaw(value);
          }
        }else if(display_type == 'decimal'){
          var formatter = new Iccbl.DecimalFormatter(cell_options);
          bindings['#'+key].onGet = function(value){
            return formatter.fromRaw(value);
          }
        }else if(display_type == 'link' && data_type != 'list' ){
            var c_options = _.extend({
              hrefTemplate: '#',
              target: '_self'
            }, cell_options );                
            
            // Backward compatability: Check for non-json backgrid_cell_options.
            // interpret as simple string to be interpolated, and not a JSON object.
            if( c_options.hrefTemplate == '#' && 
                _.has(fi,'display_options')) {
              // NOTE: format for backgrid cell options is "/{attribute_key}/"
              // NOTE: interpreting *all* links as *hash* values only, for now - 
              // TODO: make this switchable, using a flag in c_options
              c_options.hrefTemplate = window.location.pathname + '#' + fi['display_options'];
              c_options.target = '_self';
            } 
            var interpolatedVal = Iccbl.replaceTokens(self.model,c_options.hrefTemplate);
 
            bindings['#'+key].updateMethod = 'html';
            var originalOnGet = bindings['#'+key].onGet;
            bindings['#'+key].onGet = function(rawValue){
              rawValue = originalOnGet(rawValue);
              if(false && Iccbl.appModel.DEBUG) 
                console.log('urilist, raw value: ' + rawValue  + ', ' + interpolatedVal 
                  + ', hrefTemplate: ' + c_options.hrefTemplate);
              
              if(rawValue && rawValue != '-'){
                var _html = ( '<a ' + 
                    'id="link-' + key + '" ' + 
                    'href="' + interpolatedVal + '" ' +
                    'target="' + c_options.target + '" ' +
                    'tabIndex=-1 ' +
                    '>' + rawValue + '</a>'
                    );
                return _html;
              }
              return rawValue;
            };
          
        }else if(display_type == 'link' && data_type == 'list'){

          var c_options = _.extend({
            hrefTemplate: 'http//',
            target: '_blank'
          }, cell_options );
          bindings['#'+key].updateMethod = 'html';
          bindings['#'+key].onGet = function(rawValue){            
            if(rawValue && !_.isEmpty(rawValue) && rawValue != '-'){
              var i = 0;
              var _html = [];
              _.each(rawValue, function(val){
                var interpolatedVal = Iccbl.replaceTokens(self.model,c_options.hrefTemplate,val);
                _html.push( '<a ' + 
                    'id="link-' + key + '-'+ i + '" ' + 
                    'href="' + interpolatedVal + '" ' +
                    'target="' + c_options.target + '" ' +
                    '>' + val + '</a>'
                    );
                i++;
              });
              return _html.join(', ');
            }
            return '-';
          };
        }else if(display_type == 'image'){
          console.log('todo: implement ImageCell');
        
        // FIXME: "nested" and "nested_list" are unused 20150629 after ui_type refactor to data_type
        }else if(display_type === 'nested'){
            nestedModels[key] = new Backbone.Model(self.model.get(key));
            nestedModels[key].resourceId = fi.nested_resource || key;
        }else if(display_type === 'nested_list'){
            nestedLists[key] = {};
            nestedLists[key].list = self.model.get(key);
            nestedLists[key].resourceId = fi.nested_resource || key;
            if(!nestedLists[key].list || nestedLists[key].list.length == 0){
              delete nestedLists[key];
              self.model.unset(key); // to signal empty
            }else if (nestedLists[key].length==1) {
              nestedModels[key] = nestedLists[key][0];
              nestedModels[key].resourceId = nestedLists[key].resourceId;
              delete nestedLists[key];
              self.model.unset(key); // to signal empty
            }
        }else if(display_type == 'string' ){
          // string uses default display type
        }else{
          console.log('field', key,'using default display for display_type:',display_type);
        }
        
        if(edit_type == 'select' || edit_type == 'multiselect' ){
          try{
            if(fi.vocabulary_scope_ref){
              var vocabulary = Iccbl.appModel.getVocabulary(fi.vocabulary_scope_ref);
              if(vocabulary && !_.isUndefined(vocabulary)){
                bindings['#'+key].onGet = function(value){
                  if(_.has(vocabulary, value)){
                      return vocabulary[value].title;
                  }
                  console.log('no vocab: ' + key + ', ' + vocabulary );
                  return value;
                };
              }else{
                console.log('Warning, no vocabulary found for: ' 
                    + fi.key + ', ' + fi.vocabulary_scope_ref);
              }
            }
          }catch(e){
            console.log('error getting vocab', e);
            appModel.error('Warning, no vocabulary found for: ' 
              + fi.key + ', ' + fi.vocabulary_scope_ref);
          }

        }
        
        // Also build a binding hash for the titles
        schemaBindings['#title-'+key] = {
          observe: key,
          onGet: function(value) {
            if (value) return value.title + ":";
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
    
    afterRender : function() {
      var self = this;
      this.stickit();
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

      return this;
    },

    serialize: function() {
      return {
        'buttons': _.chain(this.buttons), // TODO: buttons from the schema
        'title': Iccbl.getTitleFromTitleAttribute(this.model, this.model.resource.schema),
        'keys': _.chain(this.detailKeys),
      };      
    },    
    
    template: _.template(detailTemplate)

	});

	return DetailView;
});