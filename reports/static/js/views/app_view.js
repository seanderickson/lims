define([
    'jquery',
    'underscore',
    'backbone',
    'views/menu',
    'views/content'
], function($, _, Backbone, MenuView, ContentView) {

    var AppView = Backbone.View.extend({
        el : $('#container'),

        initialize : function(attributes, options) {
            console.log('AppView initialize' + this.model );
            this.router = options.router;
            this.menuView = new MenuView({ model: this.model }, options);
            this.contentView = new ContentView({ model: this.model }, options);
        },

        render: function() {
            this.menuView.render();
            this.contentView.render();
        },

    });

    return AppView;
});