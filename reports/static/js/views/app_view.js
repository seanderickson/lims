define([
    'jquery',
    'underscore',
    'backbone',
    'views/menu',
    'views/content'
], function($, _, Backbone, MenuView, ContentView) {

    var AppView = Backbone.View.extend({
        el : $('#container'),

        initialize : function() {
            console.log('AppView initialize' + this.model );
            this.menuView = new MenuView({ model: this.model });
            this.contentView = new ContentView({ model: this.model });
        },

        render: function() {
            this.menuView.render();
            this.contentView.render();
        },

    });

    return AppView;
});