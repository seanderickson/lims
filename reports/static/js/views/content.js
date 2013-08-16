define([
    'jquery',
    'underscore',
    'backbone',
    'bootstrap',
    'iccbl_backgrid',
    'views/list',
    'views/home',
    'views/detail_stickit'
], function($, _, Backbone, Bootstrap, Iccbl, ListView, HomeView, DetailView) {

    var ContentView = Backbone.View.extend({
        el: '#container',

        initialize: function(attributes, options) {
            console.log('ContentView initialize');
            // _.bindAll(this,'change_menu_item', 'render'); // binds all of the objects function properties to this instance
            // this.model.bind('change:menu_item', this.change_menu_item );
//            this.listenTo(this.model, 'change:current_ui_resource_id', this.change_current_ui_resource);
            this.listenTo(this.model, 'change:current_view', this.current_view);
            this._options = options;
            this.router = options.router;

            this.currentView = new HomeView({ model: this.model });
        },

        current_view : function(type) {
            var self = this;
            var current_view = this.model.get('current_view');
            Iccbl.assert( !_.isUndefined(current_view), 'current_view is not defined');

            var current_ui_resource_id = this.model.get('current_ui_resource_id');
            console.log('---- change_current_ui_resource: ' + current_ui_resource_id + ', view :' + current_view);
            Iccbl.assert( !_.isUndefined(current_ui_resource_id), 'current_ui_resource_id is not defined');

            var content_options = this.model.get('content_options');
            var current_ui_resource = this.model.get('ui_resources')[current_ui_resource_id];
            var api_root_url = this.model.get('api_root_url');

            Iccbl.assert( !_.isUndefined(content_options), 'content_options is not defined');
            Iccbl.assert( !_.isUndefined(current_ui_resource), 'current_ui_resource is not defined');
            Iccbl.assert( !_.isUndefined(api_root_url), 'api_root_url is not defined');

            if( !_.isUndefined(this.currentView)) this.currentView.close();

            if (current_view === 'home'){  // TODO: make into a "menu view"
                this.currentView = new HomeView({ model: this.model });
            }else if (current_view === 'list'){
                var list_options = this.model.get('list_defaults');
                var options = _.extend( {}, list_options, current_ui_resource, /*current_ui_resource['options'],*/ content_options ); // TODO: move the nested options up into the model
                options.ui_resource_id = current_ui_resource_id;

                options.router = this.router;

                options.url = api_root_url + '/' + options.api_resource;
                options.url_schema = options.url + '/schema';

                //console.log('list view options: ' + JSON.stringify(options) );
                // this.listView.setOptions(options);
                this.listView = new ListView({ model: this.model }, options);
                this.currentView = this.listView;
                this.render();
            }else if (current_view === 'detail'){
                var detail_options = this.model.get('detail_defaults');
                var options = _.extend( {}, detail_options, current_ui_resource, content_options ); // TODO: move the nested options up into the model
                var resource_url = api_root_url + '/' + options.api_resource;

                var createDetail = function(schemaResult, model){
                    console.log('sr: ' + schemaResult);
                    var detailView =
                        new DetailView({ model: model}, {
                            schemaResult:schemaResult,
                            router:self.router,
                            isEditMode: false, title: "Detail for " + current_ui_resource.title } );
                    console.log('detail view: ' + detailView  );
                    self.currentView = detailView;
                    self.render();
                };

                if(_.isUndefined(content_options.schemaResult) ||_.isUndefined(content_options.model)){
                    var schema_url =  resource_url + '/schema';
                    var _key = options.key;
                    // handle composite keys
                    if(_.isArray(options.key)){
                        _key = _.reduce(options.key, function(memo, item){
                            if(!_.isNull(item)) memo += item + '/';
                            return memo;
                            }, '');
                    }
                    var url = resource_url  + '/' + _key;

                    this.getSchema(schema_url, function(schemaResult) {
                        console.log('schemaResult callback: ' + schemaResult + ', ' + url);
                        if(_.isUndefined(content_options.model)){
                            self.getModel(schemaResult, url, createDetail);
                        }else{
                            createDetail(schemaResult,content_options.model);
                        }
                    });
                }else{
                    createDetail(content_options.schemaResult,content_options.model);
                }

                // $.ajax({
                    // type: "GET",
                    // url: schema_url, //options.url_schema,
                    // data: "",
                    // dataType: "json",
                    // success: function(schemaResult) {
                        // console.log('got the schema for detail view');
                        // //var field_defs = self.wrapAsUnderscore(result.fields);
                        // var ModelClass = Backbone.Model.extend({urlRoot: url, defaults: {} });
                        // var instance = new ModelClass();
                        // instance.fetch({
                            // success: function(model){
                                // console.log('fetch called: ');
                                // options = _.extend(options, { title: "Details for " + schemaResult['resource_definition']['title'],
                                    // fields:schemaResult.fields,
                                    // resource_definition: schemaResult['resource_definition'],
                                    // router: self.router } );
//
                            // },
                            // error: function(model, response, options){
                                // console.log('error fetching the model: '+ model + ', response: ' + JSON.stringify(response));
                                // var msg = 'Error locating resource: ' + url;
                                // var sep = '\n';
                                // if(!_.isUndefined(response.status)) msg += sep + response.status;
                                // if(!_.isUndefined(response.statusText)) msg += sep+ response.statusText;
                                // if(!_.isEmpty(response.responseText)) msg += sep+ response.responseText;
                                // window.alert(msg); // TODO: use Bootstrap inscreen alert classed message div
                            // }
                        // });
                    // }, // end success outer ajax call
                    // error: function(x, e) {
                        // alert(x.readyState + " "+ x.status +" "+ e.msg);
                    // }
                    // });

            }else{
                window.alert('unknown view: ' + content_options['view']);
            }
        },

        render: function() {
            this.$el.append(this.currentView.render().el);
        },


        getSchema: function (schema_url, callback) {
            $.ajax({
                type: "GET",
                url: schema_url, //options.url_schema,
                data: "",
                dataType: "json",
                success: function(schemaResult) {
                    console.log('got the schema for detail view');
                    //var field_defs = self.wrapAsUnderscore(result.fields);
                    callback(schemaResult);
                }, // end success outer ajax call
                error: function(x, e) {
                    alert(x.readyState + " "+ x.status +" "+ e.msg);
                }
            });
        },

        getModel: function(schemaResult, url, callback) {
            var ModelClass = Backbone.Model.extend({urlRoot: url, defaults: {} });
            var instance = new ModelClass();
            instance.fetch({
                success: function(model){
                    console.log('model retrieved ' + schemaResult);

                    callback(schemaResult, model);
                },
                error: function(model, response, options){
                    console.log('error fetching the model: '+ model + ', response: ' + JSON.stringify(response));
                    var msg = 'Error locating resource: ' + url;
                    var sep = '\n';
                    if(!_.isUndefined(response.status)) msg += sep + response.status;
                    if(!_.isUndefined(response.statusText)) msg += sep+ response.statusText;
                    if(!_.isEmpty(response.responseText)) msg += sep+ response.responseText;
                    window.alert(msg); // TODO: use Bootstrap inscreen alert classed message div
                }
            });
        },


    });



    return ContentView;
});