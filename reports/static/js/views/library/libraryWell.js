define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_detail_stickit', 
  'views/generic_detail_layout',
  'views/generic_edit',
  'views/list2',
  'views/library/library',
  'templates/generic-tabbed.html',
  'templates/genericResource.html'
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, 
         DetailView, DetailLayout, EditView, ListView, LibraryView, layout,
         genericLayout ) {
  
  var LibraryWellView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      var self = this;
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.library = args.library;
      
      this.tabbed_resources = _.extend(
        {}, _.mapObject(this.tabbed_resources_template, function(val,key){
          return _.clone(val);
        }));
      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail'){
          var permission = self.tabbed_resources[key].permission;
          if (_.isUndefined(permission)){
            permission = self.tabbed_resources[key].resource;
          }
          if (!appModel.hasPermission(permission)){
            delete self.tabbed_resources[key];
          }
        }
      });
      if (! self.model.has('duplex_wells')){
        delete self.tabbed_resources['duplex_wells'];
      }
      
      _.bindAll(this, 'click_tab');
    },
    
    tabbed_resources_template: {
      detail: { 
        description: 'Well Details', 
        title: 'Well', 
        invoke: 'setDetail'
      },
      duplex_wells: { 
        description: 'Duplex Wells', 
        title: 'Duplex Wells', 
        invoke: 'setDuplexWells',
        resource: 'well'
      },
//      other_wells: { 
//        description: 'Other Wells', 
//        title: 'Other Wells', 
//        invoke: 'setOtherWells',
//        resource: 'well'
//      },
      annotations: { 
        description: 'Annotations', 
        title: 'Annotations', 
        invoke: 'setAnnotations',
        resource: 'well',
        permission: 'screen'
      }
    },      
    
    events: {
        'click ul.nav-tabs >li': 'click_tab',
    },

    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },

    /**
     * Layoutmanager hook
     */
    serialize: function() {
      console.log('serialize called...');
      var self = this;
      return {
        'base_url': [
           self.library.resource.key,self.library.key,self.model.resource.key,
           self.model.key].join('/'),
        'tab_resources': self.tabbed_resources
      }      
    }, 

    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      var self = this;
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if (viewId == 'edit'){
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
      console.log('afterRender, done.');
    },
    
    click_tab : function(event){
      event.stopPropagation();
      event.preventDefault();
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      this.change_to_tab(key);
    },

    change_to_tab: function(key){
      if(_.has(this.tabbed_resources, key)){
        this.$('li').removeClass('active');
        this.$('#' + key).addClass('active');
        if(key !== 'detail'){
          this.consumedStack = [key];
        }else{
          this.consumedStack = [];
        }
        var delegateStack = _.clone(this.uriStack);
        this.uriStack = [];
        var method = this[this.tabbed_resources[key]['invoke']];
        if (_.isFunction(method)) {
          method.apply(this,[delegateStack]);
        } else {
          throw "Tabbed resources refers to a non-function: " + this.tabbed_resources[key]['invoke']
        }
      }else{
        var msg = 'Unknown tab: ' + key;
        appModel.error(msg);
        throw msg;
      }
    },
    
    
    setDetail: function(delegateStack) {
      var key = 'detail';
      
      var view = this.tabViews[key];
      if (view) {
        this.removeView(this.tabViews[key]);
      }
      
      var detailView = DetailView.extend({
        afterRender: function(){
          var self = this;
          DetailView.prototype.afterRender.apply(this,arguments);
          // TODO: support for generic images
          if(this.model.has('structure_image')){
            self.$('#content').append(
                '<img style="position: absolute; top: 8em; right: 3em; height: 16em;" src="' 
                + self.model.get('structure_image') + '" alt="image" />')
          }
        }
      });
      
      view = new DetailLayout({ 
        model: this.model,
        uriStack: delegateStack, 
        buttons: ['download', 'history'],
        DetailView: detailView
      });
      this.tabViews[key] = view;

      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();
    },
    
    /** 
     * Retrieve the Study annotations for a well:
     *  Format is:
     *  [ { study 1 information, dc_1, dc_2, etc. }, { study 2 information ... }...]
     *  - where each study conforms to the study schema
     *  - each dc has:
     *      { dc_schema, value }
    **/
    setAnnotations: function(delegateStack){
      console.log('set annotations...');
      var self = this;
      var url = [self.model.resource.apiUri,
                 self.model.get('well_id'), 'annotations'].join('/');
      var studyResource = appModel.getResource('study'); 
      studyResource['fields'] = _.pick(studyResource['fields'],
        ['lead_screener_name','lab_name','title','facility_id', 
         'date_created', 'summary','study_type'])

      function showAnnotations(data){
        if (!data|| _.isEmpty(data)){
          console.log('empty annotation information');
          return;
        }
        console.log('process annotations...');
        var AnnotationView = Backbone.Layout.extend({
          template: _.template(genericLayout),
          afterRender: function(){
            $('#resource_content_title').html(
              'Annotations for ' + self.model.get('well_id'));
            $content = $('<div class="container" id="studies_container"></div>');
            $('#resource_content').html($content);

            _.each(data, function(studyData){

              var facility_id = data['facility_id']
              $content = $([
                '<div class="row">',
                '<div class="col-xs-6" id="study_info-'+facility_id + '"></div>',
                '<div class="col-xs-6" id="annotation_info-'+facility_id + '"></div>',
                '</div>',
                ].join(''));
              $('#studies_container').append($content);
              var model = new Backbone.Model(studyData);
              model.resource = studyResource;
              view = new DetailView({
                model: model,
                resource: studyResource,
                buttons: []
              });
              $content.find('#study_info-'+facility_id + '').append(view.render().$el);
              
              // Create a resource schema on the fly for the annotations
              var schema = {
                title_attribute: 'Annotation',
                key: 'annotation',
                fields: studyData.fields
              }
              _.each(_.values(schema.fields), appModel.parseSchemaField );
              schema = _.extend(schema, appModel.schemaClass);
              var model = new Backbone.Model(studyData.values);
              model.resource = schema;
              view = new DetailView({
                model: model,
                resource: schema,
                buttons: []
              });
              $content.find('#annotation_info-'+facility_id + '').append(view.render().$el);
            });
            
          }
        });
        var $el = self.setView('#tab_container', new AnnotationView()).render().$el;
      };
      $.ajax({
        type : "GET",
        url : url,
        data : "",
        dataType : "json",
        success : showAnnotations,
        fail: function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); }      
      });
      
      this.reportUriStack();
    },

    setDuplexWells: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'duplex_wells'].join('/')
      
      function showDuplexWells(response){
        console.log('process duplex data...', response);
        if (!response|| _.isEmpty(response)){
          console.log('empty duplex information');
          return;
        }
        var duplex_wells = _.result(response,'duplex_wells');
        var confirmed_positive_values = _.result(response,'confirmed_positive_values');
        if (!duplex_wells){
          appModel.error('no duplex well information in response');
          return;
        }
        if (!confirmed_positive_values){
          appModel.error('No confirmed positive values in the response');
          return;
        }
        var ColoredConfirmmationCell = Backgrid.Cell.extend({
          className: 'text-wrap-cell',
          render: function(){
            this.$el.empty();
            var key = this.column.get('name');
            if (key == 'screen_facility_id'){
              this.$el.html(this.model.get(key));
            }else{
              var confirmationValue = parseInt(this.model.get(key));
              if (confirmationValue == 3){
                this.$el.attr('style','background-color: blue;');
              }else if(confirmationValue == 2){ 
                this.$el.attr('style','background-color: red;');
              }else{
                this.$el.attr('style','background-color: white;');
              }
            }
            return this;
          }
        });
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell.extend({
            render: function(){
              this.$el.empty();
              var well_data = this.column.get("label");
              this.$el.append([well_data['vendor_id'],well_data['sequence'], 
                               well_data['well_id']].join('<br>'));
              return this;
            }
          })
        };
        var columns = [
            _.extend({},colTemplate,{
              'name' : 'screen_facility_id',
              'label' : 'Screen',
              'description' : 'Screen',
              'order': 1,
              'sortable': true,
              'cell': ColoredConfirmmationCell
            })
        ];
        
        _.each(duplex_wells, function(well_data){
          var well_id = well_data['well_id'];
          var label = well_data
          
          columns.push(
            _.extend({},colTemplate,{
              'name' : well_id,
              'label' : label,
              'description' : 'Duplex Well Data',
              'order': 1,
              'sortable': true,
              'cell': ColoredConfirmmationCell
            })
          );
        });
        var newCollection = new Backbone.Collection(
          confirmed_positive_values,
          { comparator: 'screen_facility_id' });
        var colModel = new Backgrid.Columns(columns);
        var _grid = new Backgrid.Grid({
          columns: colModel,
          collection: newCollection,
          className: 'backgrid table-striped table-condensed table-hover'
        });
        
        var GridView = Backbone.Layout.extend({
          template: _.template(genericLayout),
          afterRender: function(){
            $('#resource_content_title').html('Duplex Wells for ' + self.model.get('well_id'));
            
            $content = $([
              '<div class="row"><div class="col-xs-12" id="legend"></div></div>',
              '<div class="row"><div class="col-xs-12" id="grid"></div></div>'
              ].join(''));
            $legend = $([
              '<div class="well col-xs-12">',
              '<label  class="col-xs-1">Legend:</label>',
              '<span class="col-xs-9"><ul class="list-inline">',
              '<li><label><button style="background-color: red;" ',
                'class="btn btn-default" type="button">',
                '</button> pool result not confirmed</label></li>',
              '<li><label><button style="background-color: blue;" ',
                'class="btn btn-defaul" type="button">',
                '</button> pool result confirmed</label></li>',
              '<li><label><button style="background-color: white;" ',
                'class="btn btn-default" type="button">',
                '</button> inconclusive or no data</label></li>',
              '</ul></span></div>'
              ].join(''));
            $('#resource_content').html($content);
            $('#legend').html($legend);
            $('#grid').html(_grid.render().$el);
          }
        });
        
        var $el = self.setView('#tab_container', new GridView()).render().$el;
      };
      $.ajax({
        type : "GET",
        url : url,
        data : "",
        dataType : "json",
        success : showDuplexWells,
        fail: function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); }      
      });
      
      this.reportUriStack();
      
    }
  });

  return LibraryWellView;
});