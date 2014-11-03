define([
    'jquery',
    'underscore',
    'backbone',
//    'backbone_modelbinder',
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
	
	  initialize: function(args) {
	    var self = this;
      var schema = this.model.resource.schema;
//      var buttons = this.buttons = 
//          args.buttons || ['edit','delete','download','history','back'];
      var buttons = this.buttons = args.buttons || ['download','history','back'];
      
      if(appModel.hasPermission(self.model.resource.key, 'edit')){
        this.buttons.unshift('edit');
      }

      if(appModel.getCurrentUser().is_superuser){
        this.buttons.unshift('delete');
      }
      
      var keys = this.detailKeys = schema.detailKeys(this.model);
      var events = this.events = {};
      var bindings = this.bindings = {};
      var schemaBindings = this.schemaBindings = {};
      var nestedModels = this.nestedModels = {};
      var nestedLists = this.nestedLists = {};
      _.each(keys, function(key) {
        var fi = schema.fields[key];
        if (fi.ui_type === 'nested'){
          nestedModels[key] = new Backbone.Model(self.model.get(key));
          nestedModels[key].resourceId = fi.nested_resource || key;
        }
        if (fi.ui_type === 'nested_list'){
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
        }

        bindings['#'+key] = {
          observe: key,
          onGet: function(value) {
            if (_.isUndefined(value) || _.isNull(value)) return '-';
            else if (_.isString(value) && value === '' ) return '-';
            else if (_.isObject(value) && _.isEmpty(value)) return '-';
            else return value;
          }
        };
        
        if(fi.backgrid_cell_type == 'Iccbl.LinkCell'){
          var _route = '#'
          if( _.has(fi,'backgrid_cell_options')) {
            // NOTE: format for backgrid cell options is "/{attribute_key}/"
            var backgrid_cell_options = fi['backgrid_cell_options'];
            _route = Iccbl.replaceTokens(self.model,backgrid_cell_options);
          }else{
            console.log('no options defined for link cell');
          }
          bindings['#' + key].update = function($el, val, model, options) {
            if(val != '-'){
              $el.html('<a id="link-' + key + '" href="' + _route + '" >'+ val + "</a>");
            }
          };
        }
        
        schemaBindings['#title-'+key] = {
          observe: key,
          onGet: function(value) {
            if (value) return value.title;
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
//    events: {
//      'click a#child_logsa':'testfunction'
//    },
    
    testfunction: function(event){
      event.preventDefault();
      console.log('testfunction called, ' + event);
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
          
//          var date_col = _.find(columns, function(column){
//            return column['name'] == 'date_time';
//          });
//          if(!_.isUndefined(date_col)){
////            date_col.cell =new Iccbl.EditCell();
//            date_col.cell = new Iccbl.EditCell({
//              clickHandler: function(model){
//                console.log('clickhandler for model, ' + model.get('toString'));
//              }
//            });
//          }
          
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

      _.each(self.detailKeys, function(key) {
        var schema = self.model.resource.schema;
        var fi = schema.fields[key];
        if(fi.backgrid_cell_type == 'Iccbl.LinkCell'){
          if( _.has(fi,'backgrid_cell_options')) {
            // NOTE: format for backgrid cell options is "/{attribute_key}/"
            var backgrid_cell_options = fi['backgrid_cell_options'];
            var _route = Iccbl.replaceTokens(self.model,backgrid_cell_options);
            
            // FIXME: would rather add this to the backbone "events" hash, but that is not working
            $('#link-'+key).click(function(e){
              e.preventDefault();
              console.log('route: ' + _route);
              appModel.router.navigate(_route, {trigger:true});
            });
          }else{
            console.log('no options defined for link cell');
          }
        }

      });
//      self.$('#child_logsa').click(function(event){ console.log('clicked...');});
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