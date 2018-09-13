define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/generic_detail_stickit', 
  'views/list2',
  'views/screen/libraryScreening',
  'views/screen/plateRangeSearch',
  'views/screen/rawDataTransformer',
  'views/generic_edit',
  'utils/tabbedController',
  'utils/plateRangeTable',
  'templates/genericResource.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, DetailView,
            ListView, LibraryScreeningView, PlateRangeSearchView, 
            RawDataTransformer, EditView, TabbedController, 
            PlateRangePrototype, genericLayout) {
  
  var ScreenSummaryView = TabbedController.extend({
    
    initialize: function(args) {
      var self = this;
      this._classname = 'ScreenSummary';
      var tabbed_resources = this.tabbed_resources;
      TabbedController.prototype.initialize.apply(this,arguments);
      if (this.model.get('user_access_level_granted') > 1){
        // NOTE: if user_access_level_granted == 1, should not be able to see
        this.tabbed_resources['plates'] = tabbed_resources['plates'];
        this.tabbed_resources['library'] = tabbed_resources['library'];
      }
      if (this.model.get('user_access_level_granted') == 3){
        this.tabbed_resources['libraryscreening'] = tabbed_resources['libraryscreening'];
      }
      _.bindAll(this,'addLibraryScreening');
    },

    tabbed_resources: {
      detail: { 
        description: 'Details', 
        title: 'Details', 
        invoke: 'setDetail'
      },
      libraryscreening: {
        description : 'Visits (Library Screenings)',
        title : 'Visits',
        invoke : 'setScreenings',
        permission: 'libraryscreening'
      },
      plates: {
        description : 'Plates Screened',
        title : 'Plates Screened',
        invoke : 'setPlates',
        permission:'librarycopyplate'
      },
      library: {
        description : 'Libraries Screened',
        title : 'Libraries Screened',
        invoke : 'setLibraries',
        permission: 'screen'
      },
      plateranges: {
        description : 'Plate Range Search View',
        title : 'Plate Screening Inquiry',
        invoke : 'setPlateRangeSearch',
        permission: 'screen'
      },
      transformer: {
        description : 'Transform Raw Data',
        title : 'Transform Raw Data',
        invoke : 'setRawDataTransformer',
        permission: 'rawdatatransform'
      },
    },      
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      return {
        'base_url': self.model.resource.key + '/' + self.model.key + '/summary',
        'tab_resources': this.tabbed_resources
      }      
    },
    
    setRawDataTransformer: function(delegateStack) {
      console.log('setRawDataTransformer', delegateStack);
      var self = this;
      
      var rawdataResource = appModel.getResource('rawdatatransform');
      var schemaUrl = [rawdataResource.apiUri,
                       'schema'].join('/');
      function showTransformer(schemaResult) {
        var options = {
          failCallback: function(){
            var newModel = appModel.newModelFromResource(schemaResult);
            newModel.set('screen_facility_id', self.model.get('facility_id'));
            newModel.resource = schemaResult;
            showView(newModel);
          }
        };
        appModel.getModel(
          rawdataResource.key, self.model.key, 
          function(model){
            model.resource = schemaResult;
            if (!model.has('library_plate_size')){
              model.set('library_plate_size', schemaResult.fields['library_plate_size'].default );
            }
            showView(model);
          }, 
          options);
      };
      function showView(model){
        var view = new RawDataTransformer({
          screen: self.model,
          model: model,
          uriStack: delegateStack
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        this.consumedStack = ['transformer'];
      }
      
      appModel.getResourceFromUrl(schemaUrl, showTransformer);
      
    },
    
    setPlateRangeSearch: function(delegateStack) {
      var self = this;
      var delegateStack = _.clone(delegateStack);

      // Get the current library screenings
      var url = [self.model.resource.apiUri, 
                 self.model.key,
                 'libraryscreening'].join('/');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: url
      });
      var currentLibraryScreenings = this.currentLibraryScreenings = new CollectionClass();
      currentLibraryScreenings.fetch({
        data: { 
          limit: 0,
          order_by: ['-date_of_activity']
        }
      })
      .done(function() {

        delegateStack.push('show_existing');
        var view = new PlateRangeSearchView({
          model: self.model,
          uriStack: delegateStack,
          summaryView: self,
          currentLibraryScreenings: currentLibraryScreenings
        });
        this.consumedStack = ['plateranges'];
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        Backbone.Layout.setupView(view);
        self.setView("#tab_container", view ).render();
      })
      .fail(function() { 
        Iccbl.appModel.jqXHRfail.apply(this,arguments); 
      });      
      
    },
    
    setDetail: function(delegateStack) {
      console.log('setSummary...', delegateStack);
      var self = this;
      var key = 'summary';
      
      function viewLoadHistory(e) {
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', appModel.URI_PATH_SEARCH];
        var search = {};
        search['ref_resource_name'] = 'screenresult';
        search['key'] = self.model.key;
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        appModel.router.navigate(route, {trigger: true});
        self.remove();
      }
      
      var summaryKeys = self.model.resource.filterKeys('visibility', 'summary');
      var summaryModel = appModel.getModel(
        self.model.resource.key, self.model.key, 
        function(model) {
          view = new DetailView({ 
            events: {
              'click button#history': viewLoadHistory
            },
            model: model, 
            uriStack: _.clone(delegateStack),
            detailKeys: summaryKeys,
            buttons: ['history'],
            afterRender: function() {
              DetailView.prototype.afterRender.apply(this);
              this.$el.find('#libraries_screened_count').click(function(e) {
                e.preventDefault();
                self.change_to_tab('library');
              });
              this.$el.find('#library_plates_data_loaded').click(function(e) {
                e.preventDefault();
                self.consumedStack = ['plates'];
                self.showCopyPlatesLoaded(delegateStack);
              });

              if (self.model.get('has_screen_result')) {
                self._createPositivesSummary(this.$el.find('#positives_summary'));
              }
            }
          });
          self.tabViews[key] = view;
          self.reportUriStack([]);
          self.setView("#tab_container", view ).render();
        },{ data_for_get: { visibilities: ['summary']} }
      );
    },

    /**
     * Library Copy Plates Loaded view is a sub-view of Summary
     */
    showCopyPlatesLoaded: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'copyplatesloaded'].join('/');
      var resource = appModel.getResource('librarycopyplate');

      var fields_to_show = [
          'library_short_name', 'library_screening_status', 
          'library_comment_array','plate_number','comment_array', 
          'screening_count','assay_plate_count','copies_screened', 
          'last_date_screened','first_date_screened']
      var copies_screened_field = _.clone(resource.fields['copy_name']);
      copies_screened_field['visibility'] = ['l'];
      copies_screened_field['data_type'] = 'list'
      copies_screened_field['key'] = 'copies_screened';
      copies_screened_field['title'] = 'Copies Screened';
      copies_screened_field['description'] = 'Copies screened for this plate';
      delete resource.fields['copy_name']
      resource.fields['copies_screened'] = copies_screened_field;
      
      _.each(resource['fields'], function(field){
        if (!_.contains(fields_to_show, field['key'])){
          field['visibility'] = [];
        }else{
          field['visibility'] = ['l','d'];
        }
      });
      
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        resource: resource,
        url: url,
        extraControls: []
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Library Copy Plates loaded for Screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');
      self.$("#tab_container-title").hide();
    },
    
    
    _createPositivesSummary: function($target_el) {
      var self = this;

      var experimental_wells_loaded = self.model.get('experimental_well_count');
      if (!_.isNumber(experimental_wells_loaded) ||
          experimental_wells_loaded < 1) {
        return;
      }
      function createPositiveStat(raw_value) {
        var formatter = new Iccbl.DecimalFormatter({ decimals: 2 });
        if (!raw_value) raw_value = 0;
        if ( raw_value === 0 || !_.isNumber(experimental_wells_loaded) ||
            experimental_wells_loaded < 1) {
          return ''
        }
        return Iccbl.formatString(
            '{val} ({percent}%)', {
              val: raw_value,
              percent: formatter.fromRaw(
                100.0*raw_value/experimental_wells_loaded )
            });
      }
      
      var dcResource = appModel.getResource('datacolumn');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: dcResource.apiUri 
      });
      
      var dcCollection = new CollectionClass();
      dcCollection.fetch({
        data: { 
          limit: 0,
          screen_facility_id__eq: self.model.get('facility_id'),
          data_type__in: [
            'partition_positive_indicator','boolean_positive_indicator',
            'confirmed_positive_indicator'],
          includes: ['positives_count', 'strong_positives_count',
                     'medium_positives_count','weak_positives_count'],
          order_by: ['ordinal']
        },
        success: function(collection, response) {
          if (!collection || collection.isEmpty()) {
            return;
          }
          collection.each(function(dc) {
            dc.set('total_positives', createPositiveStat(dc.get('positives_count')));
            dc.set('strong_positives', createPositiveStat(dc.get('strong_positives_count')));
            dc.set('medium_positives', createPositiveStat(dc.get('medium_positives_count')));
            dc.set('weak_positives', createPositiveStat(dc.get('weak_positives_count')));
          });
          
          var colTemplate = {
            'cell' : 'string',
            'order' : -1,
            'sortable': false,
            'searchable': false,
            'editable' : false,
            'visible': true,
            'headerCell': Backgrid.HeaderCell
          };
          var columns = [
              _.extend({},colTemplate,{
                'name' : 'name',
                'label' : 'Data Column',
                'description' : 'Data Column',
                'order': 1,
                'sortable': true,
                'cell': Iccbl.StringCell
              }),
              _.extend({},colTemplate,{
                'name' : 'total_positives',
                'label' : 'Total Positives',
                'description' : 'Total Positives',
                'order': 2,
                'sortable': true,
                'cell': Iccbl.StringCell
              }),
              _.extend({},colTemplate,{
                'name' : 'strong_positives',
                'label' : 'Strong Positives',
                'description' : 'Strong Positives',
                'order': 3,
                'sortable': true,
                'cell': Iccbl.StringCell
              }),
              _.extend({},colTemplate,{
                'name' : 'medium_positives',
                'label' : 'Medium Positives',
                'description' : 'Medium Positives',
                'order': 4,
                'sortable': true,
                'cell': Iccbl.StringCell
              }),
              _.extend({},colTemplate,{
                'name' : 'weak_positives',
                'label' : 'Weak Positives',
                'description' : 'Weak Positives',
                'order': 5,
                'sortable': true,
                'cell': Iccbl.StringCell
              }),
          ];
          var colModel = new Backgrid.Columns(columns);
          colModel.comparator = 'order';
          colModel.sort();
          var positives_grid = new Backgrid.Grid({
            columns: colModel,
            collection: collection,
            className: 'backgrid table-striped table-condensed table-hover'
          });
          
          $target_el.empty();
          var cell = $('<div>',{ class: 'col-xs-4' });
          cell.html(positives_grid.render().$el);
          $target_el.append(cell);
        }
      }).fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); });
      
    },

    addLibraryScreening: function(e) {
      if (e) e.preventDefault();
      
      console.log('add library screening visit');
      var self = this;
      var defaults = {
        screen_facility_id: self.model.get('facility_id'),
        screen_type: self.model.get('screen_type')
      };
      var newModel = appModel.createNewModel('libraryscreening', defaults);

      var view = new LibraryScreeningView({ 
        model: newModel, 
        screen: self.model,
        uriStack: ['+add']
      });
      
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();

      $title = self.$el.find('#tab_container-title');
      $title.html('<H5 id="title">Add Library Screening Visit</H5>');
      $title.show();
      
      this.consumedStack = ['libraryscreening'];
      view.reportUriStack(['+add']);
    },
    
    setScreenings: function(delegateStack) {
      var self = this;
      var lsResource = appModel.getResource('libraryscreening'); 
      
      if (!_.isEmpty(delegateStack) && delegateStack[0]=='+add') {
        self.addLibraryScreening();
      }
      else if (!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
          !_.contains(appModel.LIST_ARGS, delegateStack[0])) {
        // Detail view
        
        var activityId = delegateStack.shift();
        self.consumedStack.push(activityId);
        var _key = [self.model.key,'libraryscreening',activityId ].join('/');
        
        appModel.getModel(lsResource.key, activityId, function(model) {
          var view = new LibraryScreeningView({ 
            model: model, 
            screen: self.model,
            uriStack: _.clone(delegateStack),
          });
          
          Backbone.Layout.setupView(view);
          self.listenTo(view , 'uriStack:change', self.reportUriStack);
          self.setView("#tab_container", view ).render();
          
          $title = self.$el.find('#tab_container-title');
          $title.html(view.getTitle());
          $title.show();
        });        
        return;
      } else {
        // List view
        var url = [self.model.resource.apiUri, 
                   self.model.key,
                   'libraryscreening'].join('/');

        var extraControls = [];
        
        if (appModel.hasPermission('libraryscreening','write')) {
          var $addLibraryScreeningButton = $(
            '<a class="btn btn-default btn-sm" role="button" \
            id="addLibraryScreening" href="#">Add Library Screening Visit</a>');
          $addLibraryScreeningButton.click(self.addLibraryScreening);
          extraControls.push($addLibraryScreeningButton);
        }
        lsResource.fields['activity_class']['visibility'] = ['l','d'];
        lsResource.fields['type']['visibility'] = [];
        
        _.each(_.values(lsResource.fields), function(fi){
          if (_.result(fi.display_options, 'optgroup')=='screen'){
            fi.visibility = [];
          }
        });
        
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          resource: lsResource,
          url: url,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        this.setView("#tab_container", view ).render();
        $title = self.$el.find('#tab_container-title');
        $title.empty();
        $title.hide();
      }
    },
    
    setLibraries: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'libraries'].join('/');
      var resource = appModel.getResource('library');
      resource.options = {};
      var includes = resource.options.includes = [];
      var fields = resource['fields'];
      var visible_fields = ['short_name','library_name','experimental_well_count',
                            'provider', 'screen_type', 'library_type',
                            'is_pool', 'start_plate', 'end_plate', 
                            'screening_status', 'date_screenable'];
      _.each(_.keys(fields),function(fieldkey){
        if (!_.contains(visible_fields,fieldkey)){
          fields[fieldkey]['visibility'] = [];
        }else{
          if (!_.contains(fields[fieldkey]['visibility'],'l')){
            fields[fieldkey]['visibility'] = ['l'];
            includes.unshift(fieldkey);
          }
          fields[fieldkey].ordinal = -visible_fields.length + _.indexOf(visible_fields,fieldkey);
        }
      });
      
      console.log('includes', includes);
      resource.fields['short_name']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'comment_array',
          title_function: function(model){
            return 'Comments for library: ' + model.get('short_name');
          }
        });
      
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        resource: resource,
        url: url,
        extraControls: [],
      });
      Backbone.Layout.setupView(view);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      this.$el.find('#tab_container-title').hide();
      
    },

    /**
     * Plates screened
     */
    setPlates: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'plates_screened'].join('/');
      var resource = appModel.getResource('librarycopyplate');
      
      var fields_to_show = [
            'library_short_name', 'library_screening_status',
            'plate_number','screening_count', 'assay_plate_count', 'copies_screened',
            'first_date_screened','last_date_screened'];
      
      // NOTE: API LibraryCopyPlateResource.build_screened_plate_response
      // generates a custom resource with default field visibilities for this view
      
      var copies_screened_field = _.clone(resource.fields['copy_name']);
      copies_screened_field['visibility'] = ['l'];
      copies_screened_field['data_type'] = 'list'
      copies_screened_field['key'] = 'copies_screened';
      copies_screened_field['title'] = 'Copies Screened';
      copies_screened_field['description'] = 'Copies screened for this plate';
      delete resource.fields['copy_name']
      resource.fields['copies_screened'] = copies_screened_field;
      
      _.each(resource['fields'], function(field){
        if (!_.contains(fields_to_show, field['key'])){
          field['visibility'] = [];
        }else{
          field['visibility'] = ['l','d'];
        }
      });
      
      resource.fields['library_short_name']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'library_comment_array',
          title_function: function(model){
            return 'Comments for library: ' + model.get('library_short_name');
          }
        });
      
      resource.fields['copies_screened']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'copy_comments',
          get_href: function() {
            var self = this;
            var template = '#library/{library_short_name}/copy'
            var href = Iccbl.formatString(template,self.model);
            href += '/search/copy_name__in=' + self.model.get('copies_screened').join(',')
            return href;
          },
          title_function: function(model){
            return 'Comments for Copies: ' + model.get('copies_screened').join(', ');
          }
        });
      
      resource.fields['plate_number']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'comment_array',
          title_function: function(model){
            return 'Comments for all Plate: ' 
              + model.get('library_short_name') + '/' 
              + model.get('plate_number') 
              + ', Copies: ' + model.get('copies_screened').join(',')
          },
          get_href: function(){
            var self = this;
            var template = '#library/{library_short_name}/plate'
            var href = Iccbl.formatString(template,self.model);
            
            href += '/search/plate_number__eq=' + self.model.get('plate_number');
            href += ';copy_name__in=' + self.model.get('copies_screened').join(',')
            
            return href;
          }
        });
      
      if(_.has(resource.fields,'screening_count')){
        // Note: screening_count may be restricted
        resource.fields['screening_count'].backgridCellType = 
          Iccbl.LinkCell.extend(_.extend({},
            resource.fields['screening_count'].display_options,
            {
              get_search_entry: function(){
                var search_entries = [];
                var plate_number = this.model.get('plate_number');
                console.log('copies_screened', this.model.get('copies_screened'));
                _.each(this.model.get('copies_screened'), function(copy_name){
                  search_entries.push(copy_name + '/' + plate_number);
                });
                var search_entry = 'library_plates_screened__contains=' + search_entries.join(',');
                return search_entry;
              },
              get_href: function(){
                return ['#screen', self.model.get('facility_id'),
                  'summary/libraryscreening',appModel.URI_PATH_SEARCH, 
                  this.get_search_entry()].join('/');
              },
              // implement linkCallBack for better UI response
              linkCallback: function(e){
                e.preventDefault();
                self.uriStack = [appModel.URI_PATH_SEARCH, this.get_search_entry()];
                self.change_to_tab('libraryscreening');
              }
            }));
      }
      
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        resource: resource,
        url: url,
        extraControls: []
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView("#tab_container", view ).render();
      this.$el.find('#tab_container-title').hide();
    },

  });
  
  return ScreenSummaryView;
});