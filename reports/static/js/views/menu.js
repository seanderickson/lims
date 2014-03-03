define([
    'jquery',
    'underscore',
    'backbone',
    'layoutmanager',
    'iccbl_backgrid',
    'models/app_state',
    'text!templates/menu.html'
], function($, _, Backbone, layoutmanager, Iccbl, appModel, menuTemplate) {
  
    var MenuView = Backbone.Layout.extend({

        events: {
          'click li': 'menuClick'
        },

        initialize: function(attributes, options) {
          console.log('MenuView initialize');
//          this.listenTo(appModel, 'change:current_resource_id', this.change_ui_resource);
          this.listenTo(appModel, 'change:uriStack', this.uriStackChange);
        },
        
        template: _.template(menuTemplate),
        
        serialize: function(){
          return {
            menu: appModel.get('menu'),
            ui_resources: appModel.get('ui_resources')
          };          
        },
        
        uriStackChange: function(model, val, options) {
          if(options.source === this){
            console.log('self generated uristack change');
            return;
          }else{
            this.change_ui_resource();
          }
        },
   

        change_ui_resource : function() {
          
          var uriStack = appModel.get('uriStack');
          var ui_resource_id = uriStack[0];
          console.log('got ui_resource: ' + ui_resource_id);
          
//            var ui_resource_id = appModel.get('current_resource_id');

            this.$('li').removeClass('active');

            var menus = appModel.get('menu');
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
                appModel.set({'menu':menus});
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

          var menus = appModel.get('menu');

          var menu = this.find_submenu(menus, ui_resource_id);
          if(_.isUndefined(menu)){
              window.alert('unknown submenu: ' + ui_resource_id);
              return;
          }

          // if menu doesn't have an "expanded" flag, then just do it's action
          if( ! _.has(menu, 'expanded')){
            appModel.set({ actualStack: [ui_resource_id] });
            appModel.set({ path: ui_resource_id });
//              appModel.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
//              appModel.set({
//                  current_resource_id: ui_resource_id,
//                  current_view: menu.view,
//                  current_options: {}
//              });
          }else{
              // first click on a menu item expands it
              // second click on a menu item will cause its action, if it is expanded
              if( menu['expanded'] == false ){
                  menu['expanded'] = true;
                  appModel.set({'menu':menus});
                  this.render();
              }else if( menu['expanded'] == true ){
                  menu['expanded'] = false;
                  appModel.set({'menu':menus});
                  this.render();

                  appModel.set({ uriStack: [ui_resource_id] });
                  appModel.set({ path: ui_resource_id });

//                  appModel.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
//                  appModel.set({
//                      current_resource_id: ui_resource_id,
//                      current_view: menu.view,
//                      current_options: {}
//                  });
              }
          }
          this.$('li').removeClass('active');
          this.$('#' + ui_resource_id).addClass('active');

      },

      menuClick_old: function(event){
          var self = this;
          event.preventDefault();
          var ui_resource_id = event.currentTarget.id;
          console.log('menu click: ' + ui_resource_id);

          var menus = appModel.get('menu');

          var menu = this.find_submenu(menus, ui_resource_id);
          if(_.isUndefined(menu)){
              window.alert('unknown submenu: ' + ui_resource_id);
              return;
          }

          // if menu doesn't have an "expanded" flag, then just do it's action
          if( ! _.has(menu, 'expanded')){
              appModel.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
              appModel.set({
                  current_resource_id: ui_resource_id,
                  current_view: menu.view,
                  current_options: {}
              });
          }else{
              // first click on a menu item expands it
              // second click on a menu item will cause its action, if it is expanded
              if( menu['expanded'] == false ){
                  menu['expanded'] = true;
                  appModel.set({'menu':menus});
                  this.render();
              }else if( menu['expanded'] == true ){
                  menu['expanded'] = false;
                  appModel.set({'menu':menus});
                  this.render();

                  appModel.set({ current_view: {} }, {silent: true} ); // since current view must change to trigger content view
                  appModel.set({
                      current_resource_id: ui_resource_id,
                      current_view: menu.view,
                      current_options: {}
                  });
              }
          }
          this.$('li').removeClass('active');
          this.$('#' + ui_resource_id).addClass('active');
      }

    });

    return MenuView;
});