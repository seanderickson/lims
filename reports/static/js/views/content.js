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

            this.listenTo(this.model, 'change:current_view', this.current_view);
            this._options = options;
            this.router = options.router;

            this.currentView = new HomeView({ model: this.model });
            this.render();
        },

        current_view : function(type) {
            var self = this;

            var current_resource_id = this.model.get('current_resource_id');
            Iccbl.assert( !_.isUndefined(current_resource_id), 'list: current_resource_id is not defined');

            var current_view = this.model.get('current_view');
            Iccbl.assert( !_.isUndefined(current_view), 'list: current_view is not defined');

            var current_options = this.model.get('current_options');
            Iccbl.assert( !_.isUndefined(current_options), 'list: current_options is not defined');

            var current_ui_resource = this.model.get('ui_resources')[current_resource_id];
            Iccbl.assert( !_.isUndefined(current_ui_resource), 'list: current_ui_resource is not defined for current_resource_id: ' + current_resource_id );

            if( !_.isUndefined(this.currentView)) this.currentView.close();

            if (current_view === 'home'){  // TODO: make into a "menu view"
                this.currentView = new HomeView({ model: this.model });
            }else if (current_view === 'list'){
                var options = _.extend( {}, this.model.get('list_defaults'), current_ui_resource, current_options ); // TODO: move the nested options up into the model
                options.ui_resource_id = current_resource_id;
                options.router = this.router;
                options.url = options.url_root + '/' + options.api_resource;
                options.url_schema = options.url + '/schema';

                this.listView = new ListView({ model: this.model }, options);
                this.currentView = this.listView;
                this.render();
            }else if (current_view === 'detail'){
                var options = _.extend( {}, this.model.get('detail_defaults'), current_ui_resource ); // TODO: move the nested options up into the model

                var current_scratch = this.model.get('current_scratch');
                this.model.set({ current_scratch: {} });

                var createDetail = function(schemaResult, model){
                    var detailView =
                        new DetailView({ model: model},
                            {
                                schemaResult:schemaResult,
                                router:self.router,
                                isEditMode: false, title: "Detail for " + current_ui_resource.title
                            });
                    self.currentView = detailView;
                    self.render();
                };

                if(_.isUndefined(current_scratch.schemaResult) ||_.isUndefined(current_scratch.model)){  // allow reloading
                    var resource_url = current_ui_resource.url_root + '/' + current_ui_resource.api_resource;
                    var schema_url =  resource_url + '/schema';
                    var _key = current_options;
                    Iccbl.assert( !_.isEmpty(_key), 'content:detail: options.key required if not schemaResult, model supplied');
                    // handle composite keys
                    if(_.isArray(_key)){
                        _key = _.reduce(_key, function(memo, item){
                            if(!_.isNull(item)) memo += item + '/';
                            return memo;
                        }, '');
                    }
                    var url = resource_url  + '/' + _key;

                    this.getSchema(schema_url, function(schemaResult) {
                        console.log('schemaResult callback: ' + schemaResult + ', ' + url);
                        if(_.isUndefined(current_scratch.model)){
                            self.getModel(schemaResult, url, createDetail);
                        }else{
                            createDetail(schemaResult,current_scratch.model);
                        }
                    });
                }else{
                    createDetail(current_scratch.schemaResult,current_scratch.model);
                }

            }else{
                window.alert('unknown view: ' + current_view);
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
                    callback(schemaResult);
                }, // end success outer ajax call
                error: function(x, e) {
                    alert(x.readyState + " "+ x.status +" "+ e.msg);
                }
            });
        },

        getModel: function(schemaResult, url, callback) {
            var ModelClass = Backbone.Model.extend({url: url, defaults: {} });
            var instance = new ModelClass();
            instance.fetch({
                success: function(model){
                    callback(schemaResult, model);
                },
                error: function(model, response, options){
                    //console.log('error fetching the model: '+ model + ', response: ' + JSON.stringify(response));
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