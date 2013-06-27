define([
    'jquery',
    'underscore',
    'backbone',
    'bootstrap',
    'views/list',
    'views/home',
], function($, _, Backbone, Bootstrap, ListView, HomeView) {

    var ContentView = Backbone.View.extend({
        el: '#container',

        initialize: function() {
            console.log('ContentView initialize');
            // _.bindAll(this,'change_menu_item', 'render'); // binds all of the objects function properties to this instance
            // this.model.bind('change:menu_item', this.change_menu_item );
            this.listenTo(this.model, 'change:current_submodel', this.change_current_submodel);

            this.listView = new ListView({ model: this.model });
            this.homeView = new HomeView({ model: this.model });
            this.currentView = this.homeView;
        },

        change_current_submodel : function() {
            var current_submodel_id = this.model.get('current_submodel');
            var current_submodel = this.model.get('submodels')[current_submodel_id];
            console.log('Content View change_current_submodel: ' + current_submodel_id ); //JSON.stringify(current_submodel) );

            // TODO: activate the appropriate content view here
            this.currentView.close();
            if (current_submodel_id === 'home'){
                this.currentView = this.homeView;
            }else if (current_submodel.view === 'ListView'){
                this.currentView = this.listView;
                var list_options = this.model.get('list_defaults');
                var content_options = this.model.get('content_options');
                var options = _.extend( {}, list_options, current_submodel, current_submodel['options'], content_options ); // TODO: move the nested options up into the model
                console.log('list view options: ' + JSON.stringify(options) );
                this.listView.setOptions(options);
            }else if (current_submodel.view === 'SemanticListView'){
                this.currentView = this.listView;
                var list_options = this.model.get('list_defaults');
                var content_options = this.model.get('content_options');
                var options = _.extend( {}, list_options, current_submodel, current_submodel['options'], content_options ); // TODO: move the nested options up into the model
                console.log('list view options: ' + JSON.stringify(options) );
                this.listView.setOptions(options);
            }else{
                window.alert('unknown menu item: ' + menu_item);
            }
            this.render();
        },

        render: function() {
            this.$el.append(this.currentView.render().el);
        }
    });

    return ContentView;
});