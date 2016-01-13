/**
 * Screen form/view
 */
define([
  'jquery', 
  'underscore', 
  'backbone', 
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/generic_detail_layout', 
  'views/generic_detail_stickit', 
  'views/list2', 
  'views/collectionColumns',
  'text!templates/generic-tabbed.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, DetailLayout, 
            DetailView, ListView, CollectionColumnView, tabbedTemplate){

  var ScreenView = Backbone.Layout.extend({

    initialize: function(args) {
      var self = this;
      this.tabViews = {}; // view cache
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      
      _.each(_.keys(this.tabbed_resources), function(key){
        if(key !== 'detail' && !appModel.hasPermission(self.tabbed_resources[key].resource)){
          delete self.tabbed_resources[key];
        }
      });
      _.bindAll(this, 'click_tab');
    },

    template: _.template(tabbedTemplate),

    tabbed_resources: {
      detail : {
        description : 'Screen Details',
        title : 'Screen Details',
        invoke : 'setDetail'
      },
      summary : {
        description : 'Screening Summary',
        title : 'Screening Summary',
        invoke : 'setSummary'
      },
      billingItems: {
        description : 'Billing information',
        title : 'Billing',
        invoke : 'setBilling'
      },
      datacolumns : {
        description : 'Data Columns',
        title : 'Data Columns',
        invoke : 'setDatacolumns'
      },
      results : {
        description : 'Screen Results',
        title : 'Screen Results',
        invoke : 'setResults'
      },
      cherryPicks: {
        description : 'Cherry Pick Requests',
        title : 'Cherry Picks',
        invoke : 'setCherryPicks'
      },
      activities: {
        description : 'Activities',
        title : 'Activities',
        invoke : 'setActivities'
      },
    },

    events: {
      // TODO: how to make this specific to this view? 
      // (it is also catching clicks on the table paginator)
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
      var self = this;
      return {
        'base_url': self.model.resource.key + '/' + self.model.key,
        'tab_resources': this.tabbed_resources
      }      
    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function(){
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if (viewId == '+add') {
          this.$('ul.nav-tabs > li').addClass('disabled');
          this.uriStack.unshift(viewId);
          viewId = 'detail';
        }else if (viewId == 'edit'){
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }else if (_.contains(['libraries'],viewId)){
          this.consumedStack = [viewId];
          this.showLibraries(this.uriStack);
          return;
        }else if (_.contains(['copyplates'],viewId)){
          this.consumedStack = [viewId];
          this.showCopyPlates(this.uriStack);
          return;
        }else if (_.contains(['copyplatesloaded'],viewId)){
          this.consumedStack = [viewId];
          this.showCopyPlatesLoaded(this.uriStack);
          return;
        }
        

        if (!_.has(this.tabbed_resources, viewId)){
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },

    /**
     * Libraries view is a sub-view of Summary
     */
    showLibraries: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'libraries'].join('/');
      var resource = appModel.getResource('library');
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource.schema,
        resource: resource,
        url: url,
        extraControls: []
      }});
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event){
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Libraries for screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');

    },
    
    /**
     * Library Copy Plates view is a sub-view of Summary
     */
    showCopyPlates: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'copyplates'].join('/');
      var resource = appModel.getResource('librarycopyplate');
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource.schema,
        resource: resource,
        url: url,
        extraControls: []
      }});
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event){
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Library Copy Plates for Screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');

    },
    
    /**
     * Library Copy Plates Loaded view is a sub-view of Summary
     */
    showCopyPlatesLoaded: function(delegateStack){
      var self = this;
      var url = [self.model.resource.apiUri,self.model.key,'copyplatesloaded'].join('/');
      var resource = appModel.getResource('librarycopyplate');
      var view = new ListView({ options: {
        uriStack: _.clone(delegateStack),
        schemaResult: resource.schema,
        resource: resource,
        url: url,
        extraControls: []
      }});
      Backbone.Layout.setupView(view);
      self.reportUriStack([]);
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#tab_container", view ).render();
      self.listenTo(view, 'afterRender', function(event){
        view.$el.find('#list-title').show().append(
          '<H4 id="title">Library Copy Plates loaded for Screen: ' + self.model.key + '</H4>');
      });
      this.$('li').removeClass('active');
      this.$('#summary').addClass('active');

    },    
    
    click_tab : function(event){
      var self = this;
      event.preventDefault();
      var key = event.currentTarget.id;
      if(_.isEmpty(key)) return;
      if(this.$('#'+key).hasClass('disabled')){
        return;
      }
      if(this.key && this.key === key){
        return;
      }
      appModel.requestPageChange({
        ok: function(){
          self.change_to_tab(key);
        }
      });
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

    render : function(){
      window.scrollTo(0, 0);
      return this;
    },

    setDetail: function(delegateStack){
      var self = this;
      var key = 'detail';
      // set up a custom vocabulary that joins username to name; will be 
      // used as the text of the linklist
      this.model.resource.schema.fields['collaborator_usernames'].vocabulary = (
          _.object(this.model.get('collaborator_usernames'),
            this.model.get('collaborator_names')));
      
      // onEditCallBack: wraps the edit display function
      // - lazy fetch of the expensive principal investigators hash
      // - perform post-render enhancement of the display
      // TODO: improve this by extending the EditView.afterRender function,
      // as with the detailView, shown below
      var onEditCallBack = function(displayFunction){
        console.log('on edit callback...');
        appModel.getPrincipalInvestigatorOptions(function(piOptions){
          appModel.getUserOptions(function(userOptions){
            self.model.resource.schema.fields['collaborator_usernames']['choices'] = userOptions;
            self.model.resource.schema.fields['lead_screener_username']['choices'] = (
                [{ val: '', label: ''}].concat(userOptions));
            self.model.resource.schema.fields['lab_head_username']['choices'] = piOptions;
            var editForm = displayFunction();

            // Render the editForm; then add the add cell line
            var temp = editForm.afterRender;
            editForm.afterRender = function(){
              
              self._addVocabularyButton(
                editForm, 'cell_lines', 'cell_line', 'Cell Line', { description: 'ATCC Designation' });
             
              self._addVocabularyButton(
                editForm, 'transfection_agent', 'transfection_agent', 'Transfection Agent');
             
              self._addVocabularyButton(
                editForm, 'species', 'screen.species', 'Screened Species');
              temp.call(editForm,arguments);
            };
          
          });
        });
      };

      var view = this.tabViews[key];
      if (view) {
        this.removeView(this.tabViews[key]);
      }
      
      var detailView = DetailView.extend({
        afterRender: function(){
          DetailView.prototype.afterRender.call(this,arguments);
          self.createStatusHistoryTable(this.$el.find('#status'));
        }
      })
      view = new DetailLayout({ 
        model: this.model, 
        uriStack: delegateStack,
        onEditCallBack: onEditCallBack,
        DetailView: detailView
      });

      this.tabViews[key] = view;
      
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      this.consumedStack = []; 
      this.setView("#tab_container", view ).render();
      
      
      return view;

//      if(!self.screen.get('has_screen_result')){
//        self.$('#results').hide();
//        self.$('#datacolumns').hide();
//        console.log('no results');
//      }

    },
      
    /**
     * Update the screen status with a status history table: populate
     * using the apilog history of the status attribute
     **/
    createStatusHistoryTable: function($target_el){
      var self = this;
      var apilogResource = appModel.getResource('apilog');
      var CollectionClass = Iccbl.CollectionOnClient.extend({
        url: apilogResource.apiUri 
      });
      var status_collection = new CollectionClass();
      status_collection.fetch({
        data: { 
          limit: 0,
          key: self.model.get('facility_id'),
          ref_resource_name: self.model.resource.key,
          diff_keys__icontains: '"status"',
          order_by: ['date_time']
        },
        success: function(collection, response) {
          collection.each(function(model){
            var diffs = JSON.parse(model.get('diffs'));
            console.log('diffs', diffs);
            model.set('status', diffs.status[1]);
          });
          var TextWrapCell = Backgrid.Cell.extend({
            className: 'text-wrap-cell'
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
                'name' : 'status',
                'label' : 'Status',
                'description' : 'Screen status',
                'order': 1,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'date_time',
                'label' : 'Date',
                'description' : 'Date',
                'order': 1,
                'sortable': true,
                'cell': 'Date'
              })];
          var colModel = new Backgrid.Columns(columns);
          colModel.comparator = 'order';
          colModel.sort();

          $target_el.empty();
          var cell = $('<div>',{ class: 'col-sm-4' });
          
          var status_grid = new Backgrid.Grid({
            columns: colModel,
            collection: collection,
            className: 'backgrid table-striped table-condensed table-hover'
          });
          cell.html(status_grid.render().$el);
          $target_el.append(cell);
        },
        error: Iccbl.appModel.backboneFetchError, 
      });
    },
    
    _addVocabularyButton: function(
        editForm, fieldKey, vocabulary_scope_ref, vocabulary_name, options){
      var options = options || {};
      var addButton = $([
        '<a class="btn btn-default btn-sm" ',
          'role="button" id="add_' + fieldKey + '_button" href="#">',
          'Add</a>'
        ].join(''));
      addButton.click(function(event){
        event.preventDefault();
        appModel.addVocabularyItemDialog(
          vocabulary_scope_ref, vocabulary_name,
          function(new_vocab_item){
            // manually add the new vocabulary options to the multiselect
            editForm.$el.find('[key="' + fieldKey + '"]')
              .find('.chosen-select').append($('<option>',{
                value: new_vocab_item['key']
              }).text(new_vocab_item['title']));
            editForm.$el.find('[key="' + fieldKey + '"]')
              .find('.chosen-select').trigger("chosen:updated");
          }, options);
      });
      editForm.$el.find('div[key="form-group-' + fieldKey + '"]').append(addButton);
    },
    
    setCherryPicks: function(delegateStack){
      var self = this;
      var key = 'cherryPicks';
      var view = this.tabViews[key];
      
      if (!view){
        var self = this;
        var url = [self.model.resource.apiUri,self.model.key,'cherrypicks'].join('/');
        var resource = appModel.getResource('cherrypickrequest');
        var view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: resource.schema,
          resource: resource,
          url: url,
          extraControls: []
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        this.$('li').removeClass('active');
        this.$('#'+key).addClass('active');

      }else{
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
    },
      
    setActivities: function(delegateStack){
      var self = this;
      var key = 'activities';
      var view = this.tabViews[key];
      
      if (!view){
        var self = this;
        var url = [self.model.resource.apiUri,self.model.key,'activities'].join('/');
        var resource = appModel.getResource('activity');
        
        var sa_vocab = appModel.getVocabulary('activity.type');
        resource.schema.fields['type'].vocabulary = 
          _.map(sa_vocab, function(v){
            return [v.title,v.key];
          }); // TODO: app model method for this
        
        console.log('combined vocab', resource.schema.fields['type'].vocabulary );
        
        var view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: resource.schema,
          resource: resource,
          url: url,
          extraControls: []
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        this.$('li').removeClass('active');
        this.$('#'+key).addClass('active');

      }else{
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
    },
      
    setSummary : function(delegateStack){
      var self = this;
      var key = 'summary';
      var view = this.tabViews[key];
      
      if (!view){
      
          var summaryKeys = self.model.resource.schema.filterKeys('visibility', 'summary');
          var summaryModel = appModel.getModel(
            self.model.resource.key, self.model.key, 
            function(model){
              view = new DetailView({ 
                model: model, 
                uriStack: _.clone(delegateStack),
                detailKeys: summaryKeys,
                buttons: ['history'],
                afterRender: function(){
                  DetailView.prototype.afterRender.apply(this);
                  this.$el.find('#libraries_screened_count').click(function(e){
                    e.preventDefault();
                    console.log('libraries screened click', e);
                    self.consumedStack = ['libraries'];
                    self.showLibraries(delegateStack);
                  });
                  this.$el.find('#library_plates_screened').click(function(e){
                    e.preventDefault();
                    console.log('library_plates_screend screened click', e);
                    self.consumedStack = ['copyplates'];
                    self.showCopyPlates(delegateStack);
                  });
                  this.$el.find('#library_plates_data_loaded').click(function(e){
                    e.preventDefault();
                    console.log('library_plates_data_loaded click', e);
                    self.consumedStack = ['copyplates'];
                    self.showCopyPlatesLoaded(delegateStack);
                  });

                  self.createPositivesSummary(this.$el.find('#positives_summary'));
                  
                  this.$el.prepend('<button>Add Library Screening</button>');
                }
              });
              self.tabViews[key] = view;
              
              self.setView("#tab_container", view ).render();
              
              // TODO: fixup the detaillayout & detail view so that well is on detailview
              view.$el.addClass('well');
              self.reportUriStack([]);
              
            },{ visibilities: ['summary']}
          );
        
      }else{
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
    },
    
    createPositivesSummary: function($target_el){
      var self = this;

      var experimental_wells_loaded = self.model.get('experimental_well_count');
      function createPositiveStat(raw_value){
        var formatter = new Iccbl.DecimalFormatter({ decimals: 2 });
        if (!raw_value) raw_value = 0;
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
          screen_facility_id: self.model.get('facility_id'),
          data_type__in: [
            'partition_positive_indicator','boolean_positive_indicator',
            'confirmed_positive_indicator'],
          order_by: ['ordinal']
        },
        success: function(collection, response) {
          collection.each(function(dc){
            dc.set('total_positives', createPositiveStat(dc.get('positives_count')));
            dc.set('strong_positives', createPositiveStat(dc.get('strong_positives_count')));
            dc.set('medium_positives', createPositiveStat(dc.get('medium_positives_count')));
            dc.set('weak_positives', createPositiveStat(dc.get('weak_positives_count')));
          });
          
          var TextWrapCell = Backgrid.Cell.extend({
            className: 'text-wrap-cell'
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
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'total_positives',
                'label' : 'Total Positives',
                'description' : 'Total Positives',
                'order': 2,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'strong_positives',
                'label' : 'Strong Positives',
                'description' : 'Strong Positives',
                'order': 3,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'medium_positives',
                'label' : 'Medium Positives',
                'description' : 'Medium Positives',
                'order': 4,
                'sortable': true,
                'cell': TextWrapCell
              }),
              _.extend({},colTemplate,{
                'name' : 'weak_positives',
                'label' : 'Weak Positives',
                'description' : 'Weak Positives',
                'order': 5,
                'sortable': true,
                'cell': TextWrapCell
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
          var cell = $('<div>',{ class: 'col-sm-4' });
          cell.html(positives_grid.render().$el);
          $target_el.append(cell);
        },
        error: Iccbl.appModel.backboneFetchError, 
      });
      
    },

    setBilling: function(delegateStack){
      var self = this;
      var key = 'billing';
      var view = this.tabViews[key];
      if (!view){
        var billingKeys = self.model.resource.schema.filterKeys('visibility', 'billing');
        var summaryModel = appModel.getModel(
          self.model.resource.key, self.model.key, 
          function(model){
            view = new DetailLayout({ 
              model: model, 
              uriStack: delegateStack,
              detailKeys: billingKeys,
              editKeys: billingKeys,
              editableKeys: billingKeys
            });
            self.tabViews[key] = view;
            
            self.listenTo(view , 'uriStack:change', this.reportUriStack);
            self.setView("#tab_container", view ).render();
            self.reportUriStack([]);
          },{ visibilities: ['billing']});
      }else{
        self.listenTo(view , 'uriStack:change', this.reportUriStack);
        self.setView("#tab_container", view ).render();
        self.reportUriStack([]);
      }
    },

    setDatacolumns: function(delegateStack){
      // pivot the datacolumn list
      
      var self = this;
      var datacolumnResource = appModel.getResource('datacolumn'); 
      var url = [self.model.resource.apiUri,self.model.key,'datacolumns'].join('/');
      
      // construct a collection-in-columns
      Iccbl.getCollectionOnClient(url, function(collection){
        
        // create a colModel for the list
        var columns = [];
        var TextWrapCell = Backgrid.Cell.extend({
          className: 'text-wrap-cell'
        })
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell
        };
        // setup the first column - schema field titles (that are visible)
        var col0 = _.extend({},colTemplate,{
          'name' : 'label',
          'label' : 'labels',
          'description' : 'Datacolumn field',
        });
        columns.push(col0);
        
        collection.each(function(datacolumn){
          var col = _.extend({},colTemplate,{
            'name' : datacolumn.get('key'),
            'label' : datacolumn.get('title'),
            'description' : datacolumn.get('description'),
            'order': datacolumn.get('ordinal'),
            'sortable': true,
            'cell': TextWrapCell
          });
          columns.push(col);
        });
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();
        
        // create the collection by pivoting
        orderedFields = _.sortBy(datacolumnResource.schema.fields,'ordinal');
        var pivotCollection = new Backbone.Collection();
        _.each(orderedFields, function(field){
          if(_.contains(field.visibility, 'l') && field.key != 'name' ){
            var row = {'key': field.key, 'label': field.title };
            collection.each(function(datacolumn){
              row[datacolumn.get('key')] = datacolumn.get(field.key);
            });
            pivotCollection.push(row);
          }
        });
        
        // now show as a list view
        var view = new Backgrid.Grid({
          columns: colModel,
          collection: pivotCollection,
          className: 'backgrid table-striped table-condensed table-hover'
        });
        var view = view.render();
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);

        self.setView("#tab_container", view ).render();
     });
      
    },
    
    setResults : function(delegateStack){
      var self = this;
      var screenResultResource = appModel.getResource('screenresult'); 
      var schemaUrl = [appModel.dbApiUri,
                       'screenresult',
                       self.model.key,
                       'schema'].join('/');
      var url = [appModel.dbApiUri,
                 'screenresult',
                 self.model.key].join('/');

      var _id = self.model.key;
      
      // TODO: have to add the "extra_control" because the list rendering is delayed
      var show_positives_control = $([
        '<form>',
        '<div class="checkbox">',
        '<label>',
        '  <input type="checkbox">positives',
        '</label>',
        '</div>',
        '</form>'
        ].join(''));
      var show_mutual_positives_control = $([
         '<form>',
         '<div class="checkbox">',
         '<label>',
         '  <input type="checkbox">mutual positives',
         '</label>',
         '</div>',
         '</form>'
         ].join(''));      
      var createResults = function(schemaResult){
        var initialSearchHash;
        view = new ListView({ options: {
          uriStack: _.clone(delegateStack),
          schemaResult: schemaResult,
          resource: screenResultResource,
          url: url,
          extraControls: [show_positives_control, show_mutual_positives_control]
        }});
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        initialSearchHash = view.listModel.get('search');
        if(_.has(initialSearchHash, 'is_positive__eq')
            && initialSearchHash.is_positive__eq.toLowerCase()=='true'){
          show_positives_control.find('input[type="checkbox"]').prop('checked',true);
        }
        if(_.has(initialSearchHash, 'show_mutual_positives')
            && initialSearchHash.show_mutual_positives.toLowerCase()=='true'){
          show_mutual_positives_control.find('input[type="checkbox"]').prop('checked',true);
        }
        
        show_positives_control.click(function(e){
          if(e.target.checked){
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['is_positive__eq'] = 'true';
            view.listModel.set('search',searchHash);
          }else{
            // make sure unset
            var searchHash = _.clone(view.listModel.get('search'));
            if(_.has(searchHash,'is_positive__eq')){
              delete searchHash['is_positive__eq'];
              view.listModel.set('search',searchHash);
            }
          }
        });
        show_mutual_positives_control.click(function(e){
          if(e.target.checked){
            view.show_mutual_positives(true);
          }else{
            view.show_mutual_positives(false);
          }
        });
      };
      Iccbl.getSchema(schemaUrl, createResults);
    },
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });

  return ScreenView;
});