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
  'utils/tabbedController',
  'templates/generic-tabbed.html',
  'templates/genericResource.html'
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, 
         DetailView, DetailLayout, EditView, ListView, LibraryView, 
         TabbedController, layout, genericLayout ) {

  var LibraryWellView = TabbedController.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      var self = this;
      this.__klazz = 'libraryWell';
      
      this.tabbed_resources = _.extend({}, 
        _.mapObject(this.tabbed_resources_template, function(val,key){
          return _.clone(val);
        }));
      
      if (! self.model.has('other_wells')){
          delete self.tabbed_resources['other_wells_with_reagent'];
        
      }
      TabbedController.prototype.initialize.apply(this,arguments);
      
      _.bindAll(this, 'click_tab');
    },
    
    tabbed_resources_template: {
      detail: { 
        description: 'Well Details', 
        title: 'Well Details', 
        invoke: 'setDetail'
      },
      duplex_wells: { 
        description: 'Duplex Wells', 
        title: 'Duplex Wells', 
        invoke: 'setDuplexWells',
        group: 'readEverythingAdmin'
      },
      annotations: { 
        description: 'Annotations', 
        title: 'Annotations', 
        invoke: 'setAnnotations',
        permission: '' // All logged in users
      },
      other_wells: { 
        description: 'Other Wells With the same reagent identifier', 
        title: 'Other Wells With Reagent', 
        invoke: 'setOtherWells',
        permission: '' // All logged in users
      }
    },      
    
    events: {
      'click ul.nav-tabs >li': 'click_tab',
    },

    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      var base_url = [
        'library', self.model.get('library_short_name'),'well',
        self.model.key].join('/');
      return {
        'base_url': base_url,
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
        } else if (viewId == appModel.API_PARAM_SHOW_RESTRICTED){
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }
        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }

      // re-fetch the model with all necessary fields
      var data_for_get = {
        // Include "none" fields
        includes: ['duplex_wells','pool_well', 
          'screening_mg_ml_concentration', 'screening_molar_concentration', 
          '*']
      };
      var resource = self.model.resource;
      if (self.model.get('screen_type') == 'small_molecule'){
        resource = appModel.getResource('smallmoleculereagent');
      } else if (self.model.get('screen_type') == 'rnai'){
        resource = appModel.getResource('silencingreagent');
      } else {
        resource = appModel.getResource('naturalproductreagent');
      }
      if (_.contains(self.uriStack, appModel.API_PARAM_SHOW_RESTRICTED)){
        data_for_get[appModel.API_PARAM_SHOW_RESTRICTED] = true;
      }

      appModel.getModel(resource.key, self.model.get('well_id'), 
        function(model){
          self.model = model;
          if (self.model.get('is_pool') != true){
            console.log('remove duplex_wells tab for non pool well', self.model, self.model.get('is_pool'));
            delete self.tabbed_resources['duplex_wells'];
            $('#duplex_wells').remove();
          } else if (!self.model.has('duplex_wells')){
            delete self.tabbed_resources['duplex_wells'];
            $('#duplex_wells').remove();
          }
          if (_.isEmpty(self.model.get('other_wells_with_reagent'))){
            console.log('remove other_wells tab for non pool well', self.model, self.model.get('other_wells'));
            delete self.tabbed_resources['other_wells'];
            $('#other_wells').remove();
          }
          if (self.model.get('library_well_type') != 'experimental'){
            delete self.tabbed_resources['annotations'];
            $('#annotations').remove();
          }
          self.change_to_tab(viewId, self.uriStack );
        },
        {
          data_for_get: data_for_get
        }
      );
      console.log('afterRender, done.');
    },
    
    setDetail: function(delegateStack) {
      var self = this;
      var key = 'detail';
      
      var view = this.tabViews[key];
      if (view) {
        this.removeView(this.tabViews[key]);
      }
      
      if (_.contains(delegateStack, appModel.API_PARAM_SHOW_RESTRICTED)){
        self.consumedStack.unshift(appModel.API_PARAM_SHOW_RESTRICTED);
        delegateStack = _.without(delegateStack, appModel.API_PARAM_SHOW_RESTRICTED);
      }
      
      var DetailLayoutWell = DetailLayout.extend({
        history: function(event) {
          event.preventDefault();
          var newUriStack = ['apilog','order','-date_time', appModel.URI_PATH_SEARCH];
          var search = {};
          search['ref_resource_name'] = 'well';
          search['key'] = encodeURIComponent(this.model.key);
          newUriStack.push(appModel.createSearchString(search));
          var route = newUriStack.join('/');
          console.log('history route: ' + route);
          appModel.router.navigate(route, {trigger: true});
        },
        download: function(e){
          e.preventDefault();
          e.stopPropagation();
          var url = [self.model.resource.apiUri,self.model.get('well_id'),
                     'report'].join('/');
          var show_restricted = _.contains(self.consumedStack, appModel.API_PARAM_SHOW_RESTRICTED);
          if (show_restricted){
            url += '?' + appModel.API_PARAM_SHOW_RESTRICTED + '=true';
          }
          appModel.download(url, self.model.resource);
          //url += '?format=xls&use_vocabularies=true&use_titles=true&raw_lists=true';
          //appModel.downloadUrl(url);
        }
      
      });
      
      var WellView = DetailView.extend({
        afterRender: function(){
          var detail_instance = this;
          DetailView.prototype.afterRender.apply(this,arguments);
          
          if (self.model.get('is_deprecated') != true){
            detail_instance.$el.find('#is_deprecated').parent().hide();
          }
          
          var is_restricted = ( self.model.get('is_restricted_structure')
              || self.model.get('is_restricted_sequence') );
          if (!is_restricted){
            detail_instance.$el.find('#is_restricted_sequence').parent().hide();
            detail_instance.$el.find('#is_restricted_structure').parent().hide();
          }
          
          if (appModel.hasPermission(self.model.resource.key, 'read')
              && is_restricted ){
            
            var show_restricted_control = $([
              '<label class="checkbox-inline pull-left" ',
              '   title="Show restricted structure or sequence information, if applicable" >',
              '  <input type="checkbox">Show restricted structure information</input>&nbsp;',
              '</label>'
              ].join(''));
            if (self.model.get('screen_type') == 'rnai'){
              show_restricted_control = $([
                '<label class="checkbox-inline pull-left" ',
                '   title="Show restricted structure or sequence information, if applicable" >',
                '  <input type="checkbox">Show restricted sequences</input>',
                '</label>'
                ].join(''));
            }
            $('#generic-detail-buttonpanel-left').append(show_restricted_control);
            
            var initialVal = _.contains(self.consumedStack, appModel.API_PARAM_SHOW_RESTRICTED);
            show_restricted_control.find('input[type="checkbox"]').prop('checked',initialVal);
            
            show_restricted_control.find('input[type="checkbox"]').change(function(e) {
              var data_for_get = {
                includes: ['duplex_wells','pool_well', 
                  'screening_mg_ml_concentration', 'screening_molar_concentration', 
                  '*']
              }

              if (e.target.checked) {
                data_for_get[appModel.API_PARAM_SHOW_RESTRICTED] = true;
                self.consumedStack = [appModel.API_PARAM_SHOW_RESTRICTED];
                // FIXME: DetailView needs a reportUriStack method; so just use
                // the containing reportUriStack instead.
                self.reportUriStack([]);
              } else {
                data_for_get[appModel.API_PARAM_SHOW_RESTRICTED] = false;
                self.consumedStack = [];
                self.reportUriStack([]);
              }
              self.model.fetch({
                data: data_for_get
              });
            });
          }
        }
      });
      
      view = new DetailLayoutWell({ 
        model: this.model,
        uriStack: delegateStack,
        DetailView: WellView,
        buttons: ['download', 'history']
      });
      
      this.tabViews[key] = view;
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.setView("#tab_container", view ).render();

    },

    setOtherWells: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,
                 self.model.key,
                 'other_wells'].join('/')
      var view = new ListView({
        uriStack: delegateStack,
        url: url,
        resource: self.model.resource
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
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
      var wellResource = appModel.getResource('well');
      var url = [wellResource.apiUri,
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
            var $content = $('<div class="container" id="studies_container"></div>');
            _.each(data, function(studyData){
              var facility_id = studyData['facility_id']
              var $studyContainer = $([
                '<div class="row">',
                '<div class="col-xs-6" id="study_info-'+facility_id + '"></div>',
                '<div class="col-xs-6" id="annotation_info-'+facility_id + '"></div>',
                '</div>',
                ].join(''));
              $content.append($studyContainer);
              var model = new Backbone.Model(studyData);
              model.resource = studyResource;
              var StudyView = DetailView.extend({
                afterRender: function(){
                  DetailView.prototype.afterRender.apply(this,arguments);
                  $('#show_all_fields_control').remove();
                }
              });
              var studyView = new StudyView({
                model: model,
                resource: studyResource,
                buttons: []
              });
              $studyContainer.find('#study_info-'+facility_id + '')
                .append(studyView.render().$el);
              
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
              var StudyDetailView = DetailView.extend({
                afterRender: function(){
                  DetailView.prototype.afterRender.apply(this,arguments);
                  $('#show_all_fields_control').remove();
                }
              });
              var studySpecificDetail = new StudyDetailView({
                model: model,
                resource: schema,
                buttons: []
              });
              $studyContainer.find('#annotation_info-'+facility_id + '')
                .append(studySpecificDetail.render().$el);
              $('#show_all_fields_control').remove();
            });
            this.$el.find('#resource_content').html($content);
          }
        });
        self.setView('#tab_container', new AnnotationView()).render();
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

    
    /** Generate a duplex well report for the given pool well **/
    setDuplexWells: function(delegateStack){

      var self = this;
      var url = [appModel.dbApiUri,'well',self.model.key,'duplex_wells'].join('/');
      
      console.log('setDuplexWells');
      var siResource = appModel.getResource('silencingreagent');
      
      function showDuplexWells(response){

        if (!response|| _.isEmpty(response)){
          console.log('empty duplex information');
          return;
        }
        var duplex_wells = _.result(response,'duplex_wells');
        var confirmed_positive_values = _.result(response,'confirmed_positive_values');
        if (!duplex_wells){
          appModel.error('duplex_wells not found in the response');
          return;
        }
        if (!confirmed_positive_values){
          appModel.error('confirmed_positive_values not found in the response');
          return;
        }
        
        
        var grid = $('<table>');
        grid.addClass('backgrid duplex-table table-striped table-condensed table-hover');
        
        
        var thead = $('<thead>');
        var tr = $('<tr>');
        var td = $('<th>');
        td.append([
          siResource.fields['vendor_identifier']['title'] + ':',
          siResource.fields['sequence']['title'] + ':',
          siResource.fields['well_id']['title'] + ':'
        ].join('<br/>'));
        tr.append(td);
        
        _.each(duplex_wells, function(duplex_well){
          var td = $('<th>');
          td.append([
            duplex_well['vendor_id'],duplex_well['sequence'], 
          ].join('<br>'));
          td.append('<br>')
          var hrefTemplate = '#library/{library_short_name}/well/{well_id}'
          var href = Iccbl.formatString(hrefTemplate,duplex_well);
          td.append($('<a>', {
            tabIndex : -1,
            href : href,
            target : '_blank',
          }).text(duplex_well['well_id']));              
          tr.append(td);
        });
        tr.append($('<th>').text('Screen'));
        
        thead.append(tr);
        grid.append(thead);
        console.log('confirmed_positive_values', confirmed_positive_values);
        if (!_.isEmpty(confirmed_positive_values)){
          for(var i=0; i<confirmed_positive_values.length; i++){
            console.log('i', i);
            tr = $('<tr>');
            if (i==0){
              td = $('<td>').html('<strong>Duplex Activity<br/>Confirmation Data</strong>');
              td.attr('rowspan','0');
              tr.append(td);
            }
            
            cpvals = confirmed_positive_values[i];
            
            _.each(duplex_wells, function(duplex_well){
              td = $('<td>');
              var well_id = duplex_well['well_id'];
              var confirmationValue = parseInt(cpvals[well_id]);
              if (confirmationValue == 3){
                td.attr('style','background-color: blue;');
              }else if(confirmationValue == 2){ 
                td.attr('style','background-color: red;');
              }else{
                td.attr('style','background-color: white;');
              }
              tr.append(td);
            });
            td = $('<td>');
            var link = $('<a>',{
              tabIndex: -1,
              href: '#screen/' + cpvals['screen_facility_id'],
              target: '_blank',
              title: cpvals['screen_title']
            }).text(cpvals['screen_facility_id']);
            td.append(link);
            tr.append(td);
            grid.append(tr);
          }
        }
        grid.find('th').addClass('renderable');
        grid.find('td').addClass('renderable');
      
        var $legend = $([
          '<label>Legend: ',
          '<label>&nbsp;<button style="background-color: red;" ',
            'class="btn btn-default btn-lg" type="button">',
            '</button> Pool result not confirmed </label>',
          '<label>&nbsp;<button style="background-color: blue;" ',
            'class="btn btn-default btn-lg" type="button">',
            '</button> Pool result confirmed </label>',
          '<label>&nbsp;<button style="background-color: white;" ',
            'class="btn btn-default btn-lg" type="button">',
            '</button> Inconclusive or no data</label>',
          '</label>'
          ].join(''));
        var $content = $([
          '<div class="row"><div class="col-lg-12" id="grid"></div></div>',
          '<div class="row"><div class="col-lg-12 text-center" id="legend"></div></div>'
          ].join(''));
        
        if (!_.isEmpty(confirmed_positive_values)){
          $content.find('#legend').html($legend);
        } else {
          $content.find('#legend').html('No confirmation data recorded for duplex wells');
        }
        $content.find('#grid').html(grid);
        self.$el.find('#tab_container').html($content);
      };
      $.ajax({
        type : "GET",
        url : url,
        data : "",
        dataType : "json",
        success : showDuplexWells,
        fail: function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); }      
      });
      
      this.reportUriStack([]);
    }
    
  });

  return LibraryWellView;
});