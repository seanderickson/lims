/**
 * Screen form/view
 * 
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
  'views/list2', 
  'views/collectionColumns',
  'text!templates/generic-tabbed.html'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, DetailView, 
            ListView, CollectionColumnView, tabbedTemplate){

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
      return {
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
          this.uriStack.unshift(viewId);
          this.showAdd();
          return;
        }else if (viewId == 'edit'){
          this.uriStack.unshift(viewId);
          this.showEdit();
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

    click_tab : function(event){
      event.preventDefault();
      // Block clicks from the wrong elements
      // TODO: how to make this specific to this view? (it is also catching
      // clicks on the table paginator)
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

    render : function(){
      window.scrollTo(0, 0);
      return this;
    },

    showAdd: function() {
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack
      });
      Backbone.Layout.setupView(view);

      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },
    
    showEdit: function() {
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack, 
        buttons: ['download', 'upload']
      });
      Backbone.Layout.setupView(view);

      // NOTE: have to re-listen after removing a view
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },

    setDetail: function(delegateStack){
      var key = 'detail';
      
      var view = this.tabViews[key];
      if ( !view ) {
        view = new DetailView({ 
          model: this.model,
          uriStack: delegateStack, 
          buttons: ['download'] });
        this.tabViews[key] = view;
      } 
      // NOTE: have to re-listen after removing a view
      this.listenTo(view , 'uriStack:change', this.reportUriStack);
      // NOTE: if subview doesn't report stack, report it here
      //      this.reportUriStack([]);
      this.setView("#tab_container", view ).render();

//      if(!self.screen.get('has_screen_result')){
//        self.$('#results').hide();
//        self.$('#datacolumns').hide();
//        console.log('no results');
//      }

    },

    setSummary : function(){
      var self = this;
      var summaryKeys = self.model.resource.schema.filterKeys('summary')
      var detailView = new DetailView({
        model : self.model,
        detailKeys: summaryKeys
      });
      $('#tab_container').html(detailView.render().el);
    },

    setDatacolumns: function(delegateStack){
      // pivot the datacolumn list
      
      var self = this;
      var datacolumnResource = appModel.getResource('datacolumn'); 
      var url = datacolumnResource.apiUri + '?limit=0&screen_facility_id=' + self.model.key;
      
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
          if(_.contains(field.visibility, 'list') && field.key != 'name' ){
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
    
    setDatacolumnsOld : function(){
      var self = this;
      
      if(_.isUndefined(this.datacolumns)){
        var url = [appModel.dbApiUri,
                   'screenresult',
                   self.model.key,
                   'schema'].join('/');

        var createDataColumns = function(schemaResult){
          var columns = Iccbl.createBackgridColModel(schemaResult.fields,
            Iccbl.MyHeaderCell);// , col_options );
          var ListModel = Backbone.Model.extend({
            defaults : {
              rpp : 25,
              page : 1,
              order : {},
              search : {}
            }
          });
          var listModel = new ListModel();

          var collection = new Iccbl.MyCollection({
            'url' : url,
            currentPage : parseInt(listModel.get('page')),
            pageSize : parseInt(listModel.get('rpp')),
            listModel : listModel
          });

          collection.fetch({
            success : function(){
              console.log('create datacolumns, schemaResult: '
                + schemaResult);
              var datacolumnView = new CollectionColumnView({
                model : self.model
              }, {
                schemaResult : schemaResult,
                router : self.options.router,
                isEditMode : false,
                collection : collection
              });
              self.datacolumns = datacolumnView;
              $('#tab_container').html(self.datacolumns.render().el);
            },
            error : function(model, response, options){
              window.alert('Could not get: ' + usr + '\n'
                + Iccbl.formatResponseError(response));
            }
          });
        };

        Iccbl.getSchema(url, createDataColumns);
      }else{
        self.datacolumns.setElement(self.$('#tab_container')).render();
      }
      self.$('#datacolumn').addClass('active'); // first time not clicked so
                                                // set manually
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