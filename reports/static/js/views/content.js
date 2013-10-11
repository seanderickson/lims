define([
    'jquery',
    'underscore',
    'backbone',
    'bootstrap',
    'iccbl_backgrid',
    'views/list',
    'views/home',
    'views/detail_stickit',
    'views/screen',
    'views/user'
], function($, _, Backbone, Bootstrap, Iccbl, ListView, HomeView, DetailView, ScreenView, UserView) {

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
            console.log('current_view called');
            var self = this;

            var current_resource_id = this.model.get('current_resource_id');
            Iccbl.assert( !_.isUndefined(current_resource_id), 'list: current_resource_id is not defined');

            var current_view = this.model.get('current_view');
            Iccbl.assert( !_.isUndefined(current_view), 'list: current_view is not defined');

            var current_options = _.clone(this.model.get('current_options'));
            Iccbl.assert( !_.isUndefined(current_options), 'list: current_options is not defined');

            var current_ui_resource = this.model.get('ui_resources')[current_resource_id];
            Iccbl.assert( !_.isUndefined(current_ui_resource), 'list: current_ui_resource is not defined for current_resource_id: ' + current_resource_id );

            var current_scratch = this.model.get('current_scratch');
            this.model.set({ current_scratch: {} });

            if( !_.isUndefined(this.currentView)) this.currentView.close();

            if (current_view === 'home'){  // TODO: make into a "menu view"
                this.currentView = new HomeView({ model: this.model });
            }else if (current_view == 'menu'){
                console.log('todo: "menu" view!');

            }else if (current_view === 'list'){
                var options = _.extend( {}, this.model.get('list_defaults'), current_ui_resource, current_options ); // TODO: move the nested options up into the model
                if(!_.isUndefined(current_ui_resource['options'])){
                    var resource_defined_options = current_ui_resource['options'];
                    _.each(_.keys(resource_defined_options), function(key){
                        current_options[key] = _.extend({}, current_options[key], resource_defined_options[key]);
                    });
                    console.log('========== current_options: ' + JSON.stringify(current_options));
                    self.model.set({'current_options': current_options });
                }

                options.ui_resource_id = current_resource_id;
                options.router = this.router;
                options.url = options.url_root + '/' + options.api_resource;
                options.url_schema = options.url + '/schema';

                var createList = function(schemaResult){
                    options.schemaResult = schemaResult;
                    self.listView = new ListView({ model: self.model }, options);
                    self.currentView = self.listView;
                    self.render();
                };
                Iccbl.getSchema(options.url_schema, createList);

            }else if (current_view === 'detail'){

                if(current_resource_id == 'screen' ||
                    current_resource_id == 'small_molecule_screens' ||
                    current_resource_id == 'rnai_screens'){
                    var options = {
                        url_root: current_ui_resource.url_root,
                        current_options: current_options,
                        screen_model: current_scratch.model,
                        screen_schema: current_scratch.schemaResult,
                        router: this.router,
                    };

                    this.currentView = new ScreenView({ model: this.model }, options);
                    this.render();
                }else if(current_resource_id == 'users'){
                    // TODO: deal with the hoisting mess here & elsewhere
                    (function() {
                        console.log('setup user detail view');
                        var options = {
                            url_root: current_ui_resource.url_root,
                            current_options: current_options,
                            user_model: current_scratch.model,
                            user_schema: current_scratch.schemaResult,
                            router: self.router,
                        };

                        var setView = function(){
                            self.currentView = new UserView({ model: self.model }, options);
                            self.render();
                        };

                        if(_.isUndefined(options.user_schema ) ||_.isUndefined(options.user_model)){
                            var resource_url = options.url_root + '/user';
                            var id = Iccbl.getKey(options.current_options);

                            Iccbl.getSchemaAndModel2(resource_url, id, function(schemaResult, model){
                                options.user_model = model;
                                options.user_schema = schemaResult;
                                setView();
                            });
                        }else{
                            setView();
                        }
                    })();
                }else{
                    var createDetail = function(schemaResult, model){
                        var detailView =
                            new DetailView({ model: model},
                                {
                                    schemaResult:schemaResult,
                                    router:self.router,
                                    isEditMode: false
                                });
                        self.currentView = detailView;
                        self.render();
                    };

                    if(_.isUndefined(current_scratch.schemaResult) ||_.isUndefined(current_scratch.model)){  // allow reloading
                        var resource_url = current_ui_resource.url_root + '/' + current_ui_resource.api_resource;
                        var schema_url =  resource_url + '/schema';
                        var _key = Iccbl.getKey(current_options);
                        var url = resource_url  + '/' + _key;

                        Iccbl.getSchema(schema_url, function(schemaResult) {
                            console.log('schemaResult callback: ' + schemaResult + ', ' + url);
                            if(_.isUndefined(current_scratch.model)){
                                Iccbl.getModel(schemaResult, url, createDetail);
                            }else{
                                createDetail(schemaResult,current_scratch.model);
                            }
                        });
                    }else{
                        createDetail(current_scratch.schemaResult,current_scratch.model);
                    }
                }
            }else{
                window.alert('unknown view: ' + current_view);
            }
        },

        render: function() {
            this.$el.append(this.currentView.render().el);
        },

    });



    return ContentView;
});