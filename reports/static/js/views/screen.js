/**
 * Screen form/view
 *
 */
define([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid',
  'views/detail_stickit_backbone_forms',
  'views/list',
  'views/collectionColumns',
  'text!templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, DetailView, ListView, CollectionColumnView, tabbedTemplate) {
    var ScreenView = Backbone.View.extend({
      events: {
        // TODO: how to make this specific to this view? 
        // (it is also catching clicks on the table paginator)
//            'click .tabbable li': 'click_tab', 
          'click li': 'click_tab',
      },

      initialize : function(attributes, options) {
          this.options = options;
          Iccbl.assert( !_.isUndefined(options.url_root), 'ScreenView: options.url_root must be defined');
          Iccbl.assert( !_.isUndefined(options.current_options), 'ScreenView: options.current_options must be defined');

          if(_.isArray(options.current_options) || _.isString(options.current_options) ){
              options.current_options= { key: options.current_options };
          }
          console.log('initialize ScreenView: current options: ' + JSON.stringify(options.current_options));

          _.bindAll(this, 'setDetail', 'setSummary', 'setResults', 'setDatacolumns');

          this._tabbed_resources = {
              detail: { description: 'Screen Details', title: 'Screen Details', invoke : this.setDetail },
              summary: { description: 'Screening Summary', title: 'Screening Summary', invoke : this.setSummary },
              datacolumns: { description: 'Data Columns', title: 'Data Columns', invoke : this.setDatacolumns },
              results: { description: 'Screen Results', title: 'Screen Results', invoke : this.setResults },
          };

          var compiledTemplate = _.template( tabbedTemplate,
              { tab_resources: this._tabbed_resources, });

          this.$el.html(compiledTemplate);

          var tab = options.current_options['tab'];
          if(_.isUndefined(tab) || _.isEmpty(tab) ) tab = 'detail';

          this._tabbed_resources[tab].invoke();
          this.$('#' + tab).addClass('active');
          this.listenTo(this.model, 'change:current_detail', this.change_current_detail);
      },

      change_current_detail: function(){
          var _current_options = this.model.get('current_detail');
          console.log('change_current_detail: ' + JSON.stringify(_current_options) );

          //TODO: will there be a case where we have to reset all the options?
          var tab = _current_options['tab'];
          if(this.options.current_options['tab'] != tab)
              if(!_.isUndefined(tab)) this.change_to_tab(tab);
      },

      setDetail: function(){
          var self = this;
          if(_.isUndefined(this.detail)){
              var createDetail = function(schemaResult, model){
                  self.screen = model;
                  if(! self.screen.get('has_screen_result' )){
                      self.$('#results').hide();
                      self.$('#datacolumns').hide();
                      console.log('no results');
                  }
                  var detailView =
                      new DetailView({ model: model},
                          {
                              schemaResult:schemaResult,
                              router:self.options.router,
                              isEditMode: false
                          });
                  self.detail = detailView;
                  self.$('#tab_container').html(self.detail.render().el);
              };
              if(_.isUndefined(this.options.screen_schema ) ||_.isUndefined(this.options.screen_model)){  // allow reloading
                  var resource_url = this.options.url_root + '/screen';
                  var schema_url =  resource_url + '/schema';
                  var url = resource_url  + '/' + Iccbl.getKey(this.options.current_options);

                  Iccbl.getSchemaAndModel(schema_url, url, createDetail);

              }else{
                  createDetail(this.options.screen_schema,this.options.screen_model);
              }
          }else{
              self.$('#tab_container').html(this.detail.render().el);
          }

          self.$('#detail').addClass('active'); // first time not clicked so set manually
      },

      
      setSummary: function(){
          var self = this;
          if(_.isUndefined(this.summary)){

              var createDetail = function(schemaResult, model){
                      var detailView =
                          new DetailView({ model: model},
                              {
                                  schemaResult:schemaResult,
                                  router:self.options.router,
                                  isEditMode: false
                              });
                      self.summary = detailView;
                      $('#tab_container').html(self.summary.render().el);
                  };
              var resource_url = this.options.url_root + '/screensummary'; // to test
              var schema_url =  resource_url + '/schema';
              var url = resource_url  + '/'+ Iccbl.getKey(this.options.current_options);

              Iccbl.getSchemaAndModel(schema_url, url, createDetail);
          }else{
              $('#tab_container').html(self.summary.render().el);
          }
      },


      setDatacolumns: function(){
          var self = this;
          if(_.isUndefined(this.datacolumns)){
              var _id = Iccbl.getKey(self.options.current_options);
              var resource_url = this.options.url_root + '/datacolumn'; // to test
              var schema_url =  resource_url + '/schema';
              var url = resource_url + '/?facility_id=' + _id;

              var createDataColumns = function(schemaResult){

                  var columns = Iccbl.createBackgridColModel(schemaResult.fields, Iccbl.MyHeaderCell);//, col_options );
                  var ListModel = Backbone.Model.extend({
                      defaults: {
                          rpp: 25,
                          page: 1,
                          order: {},
                          search: {}}
                      });
                  var listModel = new ListModel();

                  var collection  = new Iccbl.MyCollection({
                          'url': url,
                          currentPage: parseInt(listModel.get('page')),
                          pageSize: parseInt(listModel.get('rpp')),
                          listModel: listModel
                      });

                  collection.fetch({
                      success: function(){
                          console.log('create datacolumns, schemaResult: ' + schemaResult);
                          var datacolumnView =
                              new CollectionColumnView({ model: self.model},
                                  {
                                      schemaResult:schemaResult,
                                      router:self.options.router,
                                      isEditMode: false,
                                      collection: collection
                                  });
                          self.datacolumns = datacolumnView;
                          $('#tab_container').html(self.datacolumns.render().el);
                      },
                      error: function(model, response, options){
                          window.alert('Could not get: ' + usr + '\n' + Iccbl.formatResponseError(response));
                      }
                  });
              };

              Iccbl.getSchema(schema_url, createDataColumns);
          }else{
              self.datacolumns.setElement(self.$('#tab_container')).render();
          }
          self.$('#datacolumn').addClass('active'); // first time not clicked so set manually
      },

      setResults: function(){
          var self = this;
          if(_.isUndefined(this.results)){
              var _id = Iccbl.getKey(self.options.current_options);
              var createResults = function(schemaResult){
                  var _url = self.options.url_root + '/screenresult/'+ _id;
                  var options = _.extend({},self.options.current_options,
                      {
                          schemaResult:schemaResult,
                          router:self.options.router,
                          ui_resource_id : 'screenresults',
                          url: _url,
                          title: 'Screen Results for ' + _id,
                          header_message: 'Screen Results for ' + _id
                      });

                  var listView = new ListView({ model: self.model},options);
                  self.results = listView;
                  self.results.setElement(self.$('#tab_container')).render();
              };
              var resource_url = self.options.url_root + '/screenresult/' +_id;
              if(resource_url.charAt(resource_url.length-1)!= '/') resource_url += '/'
              var schema_url =  resource_url + 'schema';
              Iccbl.getSchema(schema_url, createResults);
          }else{
              self.results.setElement(self.$('#tab_container')).render();
          }

          self.$('#results').addClass('active'); // first time not clicked so set manually
      },

      click_tab : function(event){
          event.preventDefault();

          // Block clicks from the wrong elements
          // TODO: how to make this specific to this view? (it is also catching clicks on the table paginator)
          var key = event.currentTarget.id;
          if(_.isEmpty(key)) return;

          console.log('tab click: ' + key + ', options: ' + JSON.stringify(this.options.current_options));
          this.change_to_tab(key);

          // notify the model listeners
          var _current_options = {'tab': key, 'key':Iccbl.getKey(this.options.current_options) };
          this.model.set({ current_options: {}, current_detail: {} },{ silent: true });
          this.model.set({
              current_options: _current_options,
              routing_options: {trigger: false, replace: false}
          }); // signal to the app_model that the current view has changed // todo: separate out app_model from list_model

      },

      change_to_tab: function(key){
          console.log('change to tab called: ' + key);
          this.$('li').removeClass('active');
          this.$('#' + key).addClass('active');

          if( _.has(this._tabbed_resources, key)){
              // NOTE: calling empty like this removes all of the event handlers from the views,
              // so we have to call delegateEvents when adding it back.  figure out a better way.
              $('#tab_container').empty();
              this._tabbed_resources[key].invoke();
          }else{
              window.alert('Unknown tab: ' + key);
          }
      },

      render : function() {
          window.scrollTo(0, 0);
          return this;
      }
    });

    return ScreenView;
});