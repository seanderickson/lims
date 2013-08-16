define([
    'jquery',
    'underscore',
    'backbone',
    'text!templates/menu.html'
], function($, _, Backbone, menuTemplate) {
    var MenuView = Backbone.View.extend({
        el : $('#menu'),

        initialize : function(attributes, options) {
            console.log('MenuView initialize');
            //_.bindAll(this,'change_menu_item'); // binds all of the objects function properties to this instance
            // this.model.bind('change:menu_item', this.change_menu_item );

            this.router = options.router;
            this.listenTo(this.model, 'change:current_ui_resource_id', this.change_ui_resource);
        },

        change_ui_resource : function() {
            var current_ui_resource_id = this.model.get('current_ui_resource_id');
            console.log('--MenuView: detected change_current_ui_resource_id: ' + current_ui_resource_id);
            this.$('li').removeClass('active');
            this.$('#' + current_ui_resource_id).addClass('active');
        },

        render : function() {
            console.log('MenuView render');

            var data = {
                menu: this.model.get('menu'),
                ui_resources: this.model.get('ui_resources')
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
                var ui_resource_id = $(this).attr('id');
                console.log('menu click: ' + ui_resource_id);
                var route = data.ui_resources[ui_resource_id]['route'];
                var menu = find_submenu(data.menu, ui_resource_id);
                if(_.isUndefined(menu)){
                    window.alert('unknown submenu: ' + ui_resource_id);
                    return;
                }
                console.log('clicked: ' + ui_resource_id + ', route: ' + route + ', view: ' + JSON.stringify(menu.view) );
                self.model.set({ content_options: {} });
                self.model.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
                self.model.set({
                    current_ui_resource_id : ui_resource_id,
                    current_view: menu.view,
                    content_options: {
                        view: menu.view
                    }, // unset previous content options;
                    // TODO: this is a bit like magic; if there were a mediator to talk to the app_model, it might be better?
                });

                console.log('-- set route: ' + route);
//                self.router.navigate(route, {trigger: true} );
            });
        }
    });

    return MenuView;
});