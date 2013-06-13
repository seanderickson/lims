// provisional...
define([
    'jquery',
    'underscore',
    'backbone',
    'text!templates/menu.html'
], function($, _, Backbone, menuTemplate) {
    var MenuView = Backbone.View.extend({
        el : $('#menu'),

        initialize : function() {
            console.log('MenuView initialize');
            _.bindAll(this,'change_menu_item'); // binds all of the objects function properties to this instance
            this.model.bind('change:menu_item', this.change_menu_item );
        },

        change_menu_item : function() {
            var menu_item = this.model.get('menu_item');
            console.log('MenuView: change_menu_item: ' + menu_item);
            this.$('li').removeClass('active');
            this.$('#' + menu_item).addClass('active');
        },

        render : function() {
            console.log('MenuView render');
            var menuitems = { items: _([
                { id:'home', text:'Home' },
                { id:'list', text:'List' } ]) };// todo: make these into a collection of menu models?
            var compiledTemplate = _.template(menuTemplate, menuitems);
            this.$el.append(compiledTemplate);
            var self = this;

            $('li').click(function() {
                var menuitem = $(this).attr('id');
                console.log('clicked: ' + menuitem);
                self.model.set({
                    menu_item : menuitem
                });
            });
        }
    });

    return MenuView;
});