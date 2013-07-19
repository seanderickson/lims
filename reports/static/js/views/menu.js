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
            //_.bindAll(this,'change_menu_item'); // binds all of the objects function properties to this instance
            // this.model.bind('change:menu_item', this.change_menu_item );
            this.listenTo(this.model, 'change:current_submodel', this.change_current_submodel);
        },

        change_current_submodel : function() {
            var current_submodel = this.model.get('current_submodel');
            console.log('MenuView: change_current_submodel: ' + current_submodel);
            this.$('li').removeClass('active');
            this.$('#' + current_submodel).addClass('active');
        },

        render : function() {
            console.log('MenuView render');
            // _.pairs converts the object to array of [key, value]
            var data = { items: _(_.map(_.pairs(this.model.get('submodels')), function(obj){ return { id: obj[0], text: obj[1]['title'] } } )) }; // Note, wrap it in underscore array to facilitate 'each' call in the template
            //console.log('--- data: ' + JSON.stringify(data));
            var compiledTemplate = _.template(menuTemplate, data);
            this.$el.append(compiledTemplate);
            var self = this;

            $('li').click(function(event) {
                event.preventDefault();
                var submodel = $(this).attr('id');
                var route = self.model.get('submodels')[submodel]['route'];
                console.log('clicked: ' + submodel + ', route: ' + route );
                self.model.set({ content_options: {} });
                self.model.set({
                    current_submodel : submodel,
                    route: route,
                    // TODO: this should work to suppress initial history record, but it's also messing up the back action
                    // routing_options: { replace: true }, // suppresses the url rewriting for this step, the collection will set it when adding pagination/etc. options
                    content_options: {}, // unset previous content options;
                    // TODO: this is a bit like magic; if there were a mediator to talk to the app_model, it might be better?
                });
                //console.log('model: ' + JSON.stringify(self.model));
            });
        }
    });

    return MenuView;
});