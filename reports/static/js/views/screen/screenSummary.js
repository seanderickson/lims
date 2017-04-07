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
  'views/screen/libraryScreening2',
  'utils/tabbedController'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, DetailView,
            ListView, LibraryScreeningView, TabbedController) {
  
  var ScreenSummaryView = TabbedController.extend({
    
    initialize: function(args) {
      var self = this;
      this._classname = 'ScreenSummary';
      TabbedController.prototype.initialize.apply(this,arguments);
      
    },

    tabbed_resources: {
      detail: { 
        description: 'Details', 
        title: 'Details', 
        invoke: 'setDetail'
      },
      libraryscreening: {
        description : 'Library Screenings',
        title : 'Screenings',
        invoke : 'setScreenings',
        permission: 'libraryscreening'
      },
      library: {
        description : 'Libraries Screened',
        title : 'Libraries Screened',
        invoke : 'setLibraries',
        permission: 'screen'
      },
      plates: {
        description : 'Plates Screened',
        title : 'Plates Screened',
        invoke : 'setPlates',
        permission:'librarycopyplate'
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
    
    setDetail: function(delegateStack) {
      console.log('setSummary...', delegateStack);
      var self = this;
      var key = 'summary';
      
      function viewLoadHistory(e) {
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', 'search'];
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
                self.consumedStack = ['libraries'];
                self.showLibraries(delegateStack);
              });
              this.$el.find('#library_plates_screened').click(function(e) {
                e.preventDefault();
                self.consumedStack = ['copyplates'];
                self.showCopyPlates(delegateStack);
              });
              this.$el.find('#library_plates_data_loaded').click(function(e) {
                e.preventDefault();
                self.consumedStack = ['copyplates'];
                self.showCopyPlatesLoaded(delegateStack);
              });

              if (self.model.get('has_screen_result')) {
                self._createPositivesSummary(this.$el.find('#positives_summary'));
              }
            }
          });
          self.tabViews[key] = view;
          
          self.setView("#tab_container", view ).render();
          
          self.reportUriStack([]);
          
        },{ data_for_get: { visibilities: ['summary']} }
      );
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
          
//          var TextWrapCell = Backgrid.Cell.extend({
//            className: 'text-wrap-cell'
//          });
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
                'cell': Iccbl.TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'total_positives',
                'label' : 'Total Positives',
                'description' : 'Total Positives',
                'order': 2,
                'sortable': true,
                'cell': Iccbl.TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'strong_positives',
                'label' : 'Strong Positives',
                'description' : 'Strong Positives',
                'order': 3,
                'sortable': true,
                'cell': Iccbl.TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'medium_positives',
                'label' : 'Medium Positives',
                'description' : 'Medium Positives',
                'order': 4,
                'sortable': true,
                'cell': Iccbl.TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'weak_positives',
                'label' : 'Weak Positives',
                'description' : 'Weak Positives',
                'order': 5,
                'sortable': true,
                'cell': Iccbl.TextWrapCell
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

    addLibraryScreening: function() {
      console.log('add library screening');
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

      
      this.consumedStack = ['libraryscreening','+add'];
      self.reportUriStack([]);
    },
    
    setScreenings: function(delegateStack) {
      var self = this;
      var lsResource = appModel.getResource('libraryscreening'); 
      
      if (!_.isEmpty(delegateStack) && delegateStack[0]=='+add') {
        // do not allow navigation directly to +add
        delegateStack.shift();
      }
      
      if (!_.isEmpty(delegateStack) && !_.isEmpty(delegateStack[0]) &&
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
            id="addLibraryScreening" href="#">Add Library Screening</a>');
          $addLibraryScreeningButton.click(function(e) {
            e.preventDefault();
            self.addLibraryScreening();
          });
          extraControls.push($addLibraryScreeningButton);
        }
        
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: lsResource,
          resource: lsResource,
          url: url,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
//        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        this.setView("#tab_container", view ).render();
      }
    },

    setLibraries: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'libraries'].join('/');
      var resource = appModel.getResource('screened_library');

      var library_link_cell = Iccbl.LinkCell.extend({
          render: function(){
            var self = this;
            Iccbl.LinkCell.prototype.render.apply(this, arguments);
            var comments = this.model.get('comment_array');
            if (!_.isEmpty(comments)){
              comments = Iccbl.parseComments(comments);
              this.$el.attr('title',comments);
              this.$el.append(Iccbl.createCommentIcon(
                comments,
                'Comments for library: ' + self.model.get('library_short_name')));
            }
            return this;
          }
        });
      resource.fields['library_name']['backgridCellType'] =
        library_link_cell.extend(
          resource.fields['library_name'].display_options);
      resource.fields['short_name']['backgridCellType'] =
        library_link_cell.extend(
          resource.fields['short_name'].display_options);
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: []
      });
      Backbone.Layout.setupView(view);
//      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Libraries for screen: ' + self.model.key + '</H4>');
      });

    },
    
    /**
     * Plates screened
     */
    setPlates: function(delegateStack) {
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'plates_screened'].join('/');
      var resource = appModel.getResource('librarycopyplate');
      
      var fields_to_show = ['library_short_name', 'copy_name', 'plate_number',
                            'screening_count','assay_plate_count',
                            'last_date_screened','first_date_screened']
      
      _.each(_.keys(resource.fields), function(key){
        var field = resource.fields[key];
        
        if (_.contains(fields_to_show, key)){
          field.ordinal = -fields_to_show.length + _.indexOf(fields_to_show,key);
          field.visibility = ['l'];
        } else {
          field.visibility = [];
        }
        
      });
      
      resource.fields['library_short_name']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'library_comment_array',
          title_function: function(model){
            return 'Comments for library: ' + model.get('library_short_name');
          }
        });
      
      resource.fields['copy_name']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'copy_comments',
          title_function: function(model){
            return 'Comments for Copy: ' + model.get('library_short_name')
              + '/' + model.get('copy_name');
          }
        });
      
      resource.fields['plate_number']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'comment_array',
          title_function: function(model){
            return 'Comments for Plate: ' 
              + model.get('library_short_name') + '/' 
              + model.get('copy_name')  + '/'
              + model.get('plate_number');
          }
        });
      
      var view = new ListView({ 
        uriStack: _.clone(delegateStack),
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: []
      });
      Backbone.Layout.setupView(view);
//      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event) {
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Library Copy Plates for Screen: ' + self.model.key + '</H4>');
      });

    }
    
  });
  
  return ScreenSummaryView;
});