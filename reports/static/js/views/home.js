// provisional...
define([
    'jquery',
    'underscore',
    'backbone',
], function($, _, Backbone) {
    var HomeView = Backbone.View.extend({

        initialize : function() {
            console.log('initialize HomeView');
        },

        render : function() {
            this.el.innerHTML = "Hello World!";
            return this;
        }
    });

    return HomeView;
});