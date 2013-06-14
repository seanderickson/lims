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
            var data = { items: _( this.model.get('menu_items') ) }; // Note, wrap it in underscore array to facilitate 'each' call in the template
            var compiledTemplate = _.template(menuTemplate, data);
            this.$el.append(compiledTemplate);
            var self = this;

            $('li').click(function(event) {
                event.preventDefault();
                var menuitem = $(this).attr('id');
                var route = self.model.get('menu_actions')[menuitem]['route'];
                console.log('clicked: ' + menuitem + ', route: ' + route );
                self.model.set({ content_options: {} });
                console.log("...........");
                self.model.set({
                    menu_item : menuitem,
                    route: route,
                    content_options: {}, // unset previous content options;
                    // TODO: this is a bit like magic; if there were a mediator to talk to the app_model, it might be better?
                });
                console.log('model: ' + JSON.stringify(self.model));
            });
        }
    });

    return MenuView;
});