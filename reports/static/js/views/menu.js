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
            this.listenTo(this.model, 'change:current_resource_id', this.change_ui_resource);
        },

        change_ui_resource : function() {
            var current_resource_id = this.model.get('current_resource_id');
            console.log('--MenuView: detected change_current_ui_resource_id: ' + current_resource_id);
            this.$('li').removeClass('active');
            this.$('#' + current_resource_id).addClass('active');
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
                var ui_resource_id = event.currentTarget.id;
                console.log('menu click: ' + ui_resource_id);

                var menu = find_submenu(data.menu, ui_resource_id);
                if(_.isUndefined(menu)){
                    window.alert('unknown submenu: ' + ui_resource_id);
                    return;
                }

                self.model.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
                self.model.set({
                    current_resource_id: ui_resource_id,
                    current_view: menu.view,
                    current_options: {}
                });

            });
        }
    });

    return MenuView;
});