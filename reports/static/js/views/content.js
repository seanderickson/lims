define([
    'jquery',
    'underscore',
    'backbone',
    'views/list',
    'views/home',
], function($, _, Backbone, ListView, HomeView) {

    var ContentView = Backbone.View.extend({
        el: '#container',

        initialize: function() {
            console.log('ContentView initialize');
            _.bindAll(this,'change_menu_item', 'render'); // binds all of the objects function properties to this instance
            this.model.bind('change:menu_item', this.change_menu_item );
            this.listView = new ListView({ model: this.model }); // Maybe this should be the actual list model
            this.homeView = new HomeView({ model: this.model });
            this.currentView = this.homeView;
        },

        change_menu_item : function() {
            var menu_item = this.model.get('menu_item');
            console.log('Content View change_menu_item: ' + menu_item );

            // TODO: activate the appropriate content view here
            this.currentView.close();
            if (menu_item === 'home'){
                this.currentView = this.homeView;
            }else if (menu_item === 'list'){
                this.currentView = this.listView;
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