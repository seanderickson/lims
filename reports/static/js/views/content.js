define([
    'jquery',
    'underscore',
    'backbone',
    'bootstrap',
    'views/list',
    'views/home',
    'views/detail_stickit'
], function($, _, Backbone, Bootstrap, ListView, HomeView, DetailView) {

    var ContentView = Backbone.View.extend({
        el: '#container',

        initialize: function(attributes, options) {
            console.log('ContentView initialize');
            // _.bindAll(this,'change_menu_item', 'render'); // binds all of the objects function properties to this instance
            // this.model.bind('change:menu_item', this.change_menu_item );
            this.listenTo(this.model, 'change:current_submodel_id', this.change_current_submodel);
            this._options = options;
            this.router = options.router;

            // this.listView = new ListView({ model: this.model });
            // this.homeView = new HomeView({ model: this.model });
            this.currentView = new HomeView({ model: this.model });
        },

        // wrapAsUnderscore: function(obj) {
            // var self = this;
            // var newObj = {};
            // _.each(_.pairs(obj), function(pair){
                // if(_.isArray(pair[1])){
                    // newObj[pair[0]] = _(pair[1]);
                // }else if(_.isObject(pair[1])){
                    // newObj[pair[0]] = self.wrapAsUnderscore(pair[1]);
                // }else{
                    // newObj[pair[0]] = pair[1];
                // }
            // });
            // return newObj;
        // },

        change_current_submodel : function() {
            console.log('change_current_submodel');
            var self = this;
            var current_submodel_id = this.model.get('current_submodel_id');
            var content_options = this.model.get('content_options');
            var current_submodel = this.model.get('submodels')[current_submodel_id];
            console.log('Content View change_current_submodel: ' + current_submodel_id ); //JSON.stringify(current_submodel) );

            // TODO: activate the appropriate content view here
            if( !_.isUndefined(this.currentView)) this.currentView.close();

            if (content_options.view === 'home_view'){  // TODO: make into a "menu view"
                this.currentView = new HomeView({ model: this.model });
            }else if (content_options.view === 'list_view'){
                var list_options = this.model.get('list_defaults');
                var options = _.extend( {}, list_options, current_submodel, current_submodel['options'], content_options ); // TODO: move the nested options up into the model
                options.type = current_submodel_id;
                options.router = this.router;

                console.log('list view options: ' + JSON.stringify(options) );
                // this.listView.setOptions(options);
                this.listView = new ListView({ model: this.model }, options);
                this.currentView = this.listView;
                this.render();
            }else if (content_options.view === 'detail_view'){
                var detail_options = this.model.get('detail_defaults');
                var options = _.extend( { isEditMode: false, title: "Detail for " + current_submodel.title },
                    detail_options, current_submodel, content_options ); // TODO: move the nested options up into the model

                options.router = this.router;

                var _id = options.id;
                // handle composite keys
                if(_.isArray(options.id)){
                    _id = _.reduce(options.id, function(memo, item){
                        return memo += item + '/';
                        }, '');
                }
                var url = options.url  + '/' + _id;
                console.log('get the detail view for: ' + url);
                //url = '/reports/api/v1/metahash/metahash:fields/key/'; // for example
                $.ajax({
                    type: "GET",
                    url: options.url_schema,
                    data: "",
                    dataType: "json",
                    success: function(result) {
                        console.log('got the schema for detail view');
                        //var field_defs = self.wrapAsUnderscore(result.fields);
                        var ModelClass = Backbone.Model.extend({urlRoot: url, defaults: {} });
                        var instance = new ModelClass();
                        instance.fetch({
                            success: function(model){
                                console.log('fetch called: ');
                                options = _.extend(options, { fields:result.fields });
                                var detailView = new DetailView({ model: model}, options); // TODO: get the model "edit title" from the metainformation_hash
                                console.log('detail view: ' + detailView  );
                                self.currentView = detailView;
                                self.render();
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
                    }, // end success outer ajax call
                    error: function(x, e) {
                        alert(x.readyState + " "+ x.status +" "+ e.msg);
                    }
                });

            }else{
                window.alert('unknown view: ' + content_options['view']);
            }
        },

        render: function() {
            console.log('render');
            this.$el.append(this.currentView.render().el);
        }
    });

    return ContentView;
});