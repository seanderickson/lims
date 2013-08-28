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

            Backbone.history.on('route', function(router, route, params) {
                // console.log('route changed, set title: ' + JSON.stringify(params));
                document.title = route + ':' + JSON.stringify(params);
             }, this);

        },

        render: function() {
            this.menuView.render();
            this.contentView.render();
        },

    });

    return AppView;
});