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

            // Note: Menu listens to current_submodel; view pages listen to model.route
            // this separation makes sense because the menu directly controls the ordered list selection
            // through the current_submodel property; and routes convey [control state + submodel] to views.
            // routes also represent explicitly addressable state of the app;
            this.listenTo(this.model, 'change:current_submodel_id', this.change_current_submodel);
        },

        change_current_submodel : function() {
            var current_submodel_id = this.model.get('current_submodel_id');
            console.log('MenuView: change_current_submodel_id: ' + current_submodel_id);
            this.$('li').removeClass('active');
            this.$('#' + current_submodel_id).addClass('active');
        },

        render : function() {
            console.log('MenuView render');

            var data = {
                menu: this.model.get('menu'),
                subModels: this.model.get('submodels')
            };
            var compiledTemplate = _.template(menuTemplate, data);
            this.$el.append(compiledTemplate);
            var self = this;

            var find_submenu = function(menu, id){
                if( _.has(menu, id) ) return menu[id];
                else if(_.has(menu, 'submenus')){
                    var pairs = _.pairs(menu['submenus']);
                    for(var i=0; i < pairs.length; i++){
                        var pair = pairs[i];
                        if(pair[0] == id ) return pair[1];
                        else{
                            if(_.has(pair[1], 'submenus')){
                                var temp = find_submenu(pair[1]['submenus'],id);
                                if( _.isObject(temp)) return temp;
                            }
                        }
                    }
                }
            };

            $('li').click(function(event) {
                event.preventDefault();
                var submodel_id = $(this).attr('id');
                var route = data.subModels[submodel_id]['route'];
                // var menu = data.menu[submodel_id];
                var menu = find_submenu(data.menu, submodel_id);
                if(_.isUndefined(menu)){
                    window.alert('unknown submenu: ' + submodel_id);
                    return;
                }
                console.log('clicked: ' + submodel_id + ', route: ' + route + ', view: ' + JSON.stringify(menu.view) );
                self.model.set({ content_options: {} });
                self.model.set({
                    current_submodel_id : submodel_id,
                    route: route,
                    // TODO: this should work to suppress initial history record, but it's also messing up the back action
                    // routing_options: { replace: true }, // suppresses the url rewriting for this step, the collection will set it when adding pagination/etc. options
                    content_options: {
                        view: menu.view
                    }, // unset previous content options;
                    // TODO: this is a bit like magic; if there were a mediator to talk to the app_model, it might be better?
                });
            });
        }
    });

    return MenuView;
});