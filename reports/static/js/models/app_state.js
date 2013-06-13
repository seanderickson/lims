define([
  'jquery',
  'underscore',
  'backbone',
], function($, _, Backbone ){
    var AppState = Backbone.Model.extend({
        defaults: {
            menu_item: 'home',
            route: '',
        },

        initialize : function() {
            console.log('AppState initialized');
        },
    });

    return AppState;
});