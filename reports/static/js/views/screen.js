define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'views/detail_stickit',
    'views/list',
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, DetailView, ListView, tabbedTemplate) {
    var ScreenView = Backbone.View.extend({
        events: {
//            'click .tabbable li': 'click_tab', // TODO: how to make this specific to this view? (it is also catching clicks on the table paginator)
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

            _.bindAll(this, 'setDetail', 'setSummary', 'setResults');

            this._tabbed_resources = {
                detail: { description: 'Screen Details', title: 'Screen Details', invoke : this.setDetail },
                summary: { description: 'Screen Summary', title: 'Screen Summary', invoke : this.setSummary },
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
            // if(_.isArray(_current_options) || _.isString(_current_options) ){
                // _current_options= { key: _current_options };
            // }

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
                // this.detail.delegateEvents();
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
                var resource_url = this.options.url_root + '/screensaveruser'; // to test
                var schema_url =  resource_url + '/schema';
                var url = resource_url  + '/3797';

                Iccbl.getSchemaAndModel(schema_url, url, createDetail);
            }else{
                $('#tab_container').html(self.summary.render().el);
                // this.summary.delegateEvents();
            }
        },

        setResults: function(){
            var self = this;
            if(_.isUndefined(this.results)){
                var _id = Iccbl.getKey(self.options.current_options);
                var createResults = function(schemaResult){
                    var _url = self.options.url_root + '/screenresult/'+ _id;
                    var options = _.extend({},self.options.current_options, {
                                schemaResult:schemaResult,
                                router:self.options.router,
                                ui_resource_id : 'screenresults',
                                url: _url,
                                title: 'Screen Results for ' + _id,
                                header_message: 'Screen Results for ' + _id
                    });

                    var listView =
                        new ListView({ model: self.model},options);
                    self.results = listView;
                    self.$('#tab_container').html(self.results.render().el);
                // self.results.delegateEvents();
                };
                var resource_url = self.options.url_root + '/screenresult/' +_id;
                if(resource_url.charAt(resource_url.length-1)!= '/') resource_url += '/'
                var schema_url =  resource_url + 'schema';
                Iccbl.getSchema(schema_url, createResults);
            }else{
                self.$('#tab_container').html(self.results.render().el);
                // self.results.delegateEvents();
            }

            self.$('#results').addClass('active'); // first time not clicked so set manually
        },

        // setResults2 : function(){
            // console.log('-setResults');
            // var self = this;
            // if(_.isUndefined(this.results)){
                // console.log('-setResults: create');
                // var options = _.extend({}, this.options);
                // options.ui_resource_id = 'screenresults';
                // //options.router = this.options.router;
                // var _id = Iccbl.getKey(this.options.current_options);
                // options.url = this.options.url_root + '/screenresult/'+ _id;
                // if(options.url.charAt(options.url.length-1) != '/' ) options.url_schema = options.url + '/'; // TODO: cleanup code so no trailing slashes
                // else options.url_schema = options.url;
                // options.url_schema = options.url_schema + 'schema';
//
                // options.title = 'Screen Results for ' + _id;
                // options.header_message = 'Screen Results for ' + _id;
                // this.results = new ListView({ model: this.model }, options);
// //                self.results.render();
                // $('#tab_container').html(self.results.render().el);
//
                // // console.log('-----queue render');
                // // // TODO: obviously, somewhere still rendering in the wrong order.  will have to figure this one out!!
                // // window.setTimeout(function(){
                    // // console.log('---------------set tab container');
                // // $('#tab_container').html(self.results.render().el);
                // // }, 1000);
            // }else{
                // console.log('-setResults: redisplay');
                // $('#tab_container').html(self.results.render().el);
                // this.delegateEvents();
                // this.results.delegateEvents();
            // }
//
            // self.$('#results').addClass('active'); // first time not clicked so set manually
        // },

        click_tab : function(event){
            event.preventDefault();
            var key = event.currentTarget.id;
            if(_.isEmpty(key)) return; // TODO: how to make this specific to this view? (it is also catching clicks on the table paginator)

            console.log('tab click: ' + key);
            this.change_to_tab(key);

            // notify the model listeners
            var _current_options = _.extend({}, this.options.current_options, {'tab': key});
            this.model.set({ current_options: {} },{ silent: true });
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
            console.log('1c');
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