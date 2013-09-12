define([
    'jquery',
    'underscore',
    'backbone',
    'text!templates/menu.html'
], function($, _, Backbone, menuTemplate) {
    var MenuView = Backbone.View.extend({

        events: {
            'click li': 'menuClick'
        },

        el : $('#menu'),

        initialize : function(attributes, options) {
            console.log('MenuView initialize');
            //_.bindAll(this,'change_menu_item'); // binds all of the objects function properties to this instance
            // this.model.bind('change:menu_item', this.change_menu_item );

            this.router = options.router;
            this.listenTo(this.model, 'change:current_resource_id', this.change_ui_resource);
        },

        change_ui_resource : function() {
            var ui_resource_id = this.model.get('current_resource_id');
            console.log('--MenuView: detected change_current_ui_resource_id: ' + ui_resource_id);
            this.$('li').removeClass('active');

            var menus = this.model.get('menu');
            var found_menus = this.find_submenu_path(menus, ui_resource_id);
            if(_.isUndefined(found_menus)){
                window.alert('unknown submenu: ' + ui_resource_id);
                return;
            }else{
                _.each(found_menus, function(menu){
                    if( menu['expanded'] == false ){
                        menu['expanded'] = true;
                    }
                });
                this.model.set({'menu':menus});
                this.render();
            }
            this.$('#' + ui_resource_id).addClass('active');
        },

        find_submenu_path : function(menu, id){
            if( _.has(menu, id) ) return menu[id];
            else if(_.has(menu, 'submenus')){
                var pairs = _.pairs(menu['submenus']);
                for(var i=0; i < pairs.length; i++){
                    var pair = pairs[i];
                    if(pair[0] == id ) return pair[1];
                    else{
                        if(_.has(pair[1], 'submenus')){
                            var temp = this.find_submenu(pair[1]['submenus'],id);
                            if( _.isObject(temp)) return _.flatten([pair[1],temp]);
                        }
                    }
                }
            }
        },

        find_submenu : function(menu, id){
            if( _.has(menu, id) ) return menu[id];
            else if(_.has(menu, 'submenus')){
                var pairs = _.pairs(menu['submenus']);
                for(var i=0; i < pairs.length; i++){
                    var pair = pairs[i];
                    if(pair[0] == id ) return pair[1];
                    else{
                        if(_.has(pair[1], 'submenus')){
                            var temp = this.find_submenu(pair[1]['submenus'],id);
                            if( _.isObject(temp)) return temp;
                        }
                    }
                }
            }
        },

        menuClick: function(event){
            var self = this;
            event.preventDefault();
            var ui_resource_id = event.currentTarget.id;
            console.log('menu click: ' + ui_resource_id);

            var menus = this.model.get('menu');

            var menu = this.find_submenu(menus, ui_resource_id);
            if(_.isUndefined(menu)){
                window.alert('unknown submenu: ' + ui_resource_id);
                return;
            }

            // if menu doesn't have an "expanded" flag, then just do it's action
            if( ! _.has(menu, 'expanded')){
                self.model.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
                self.model.set({
                    current_resource_id: ui_resource_id,
                    current_view: menu.view,
                    current_options: {}
                });
            }else{
                // first click on a menu item expands it
                // second click on a menu item will cause its action, if it is expanded
                if( menu['expanded'] == false ){
                    menu['expanded'] = true;
                    this.model.set({'menu':menus});
                    this.render();
                }else if( menu['expanded'] == true ){
                    menu['expanded'] = false;
                    this.model.set({'menu':menus});
                    this.render();

                    self.model.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
                    self.model.set({
                        current_resource_id: ui_resource_id,
                        current_view: menu.view,
                        current_options: {}
                    });
                }
            }
            this.$('li').removeClass('active');
            this.$('#' + ui_resource_id).addClass('active');

        },

        render : function() {
            var menus = this.model.get('menu');
            console.log('MenuView render: ' + JSON.stringify(menus));

            var data = {
                menu: menus,
                ui_resources: this.model.get('ui_resources')
            };
            var compiledTemplate = _.template(menuTemplate, data);
            this.$el.html(compiledTemplate);
            var self = this;

        }
    });

    return MenuView;
});