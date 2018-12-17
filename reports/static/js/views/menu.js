define([
    'jquery',
    'underscore',
    'backbone',
    'layoutmanager',
    'iccbl_backgrid',
    'models/app_state',
    'templates/menu.html'
], function($, _, Backbone, layoutmanager, Iccbl, appModel, menuTemplate) {
  
    var MenuView = Backbone.Layout.extend({

      events: {
        'click li': 'menuClick'
      },

      initialize: function(attributes, options) {
        console.log('initialize menu.js');
        this._classname = 'MenuView';
        this.listenTo(appModel, 'change:uriStack', this.uriStackChange);
      },
      
      template: _.template(menuTemplate),
      
      serialize: function(){
        var menuData = appModel.getMenu();
        var uiResources = appModel.get('ui_resources');
        return {
          menu: menuData,
          ui_resources: uiResources
        };          
      },
      
      /**
       * Backbone.Model change event handler
       * @param options.source = the event source triggering view
       */
      uriStackChange: function(model, val, options) {
        if(options && options.source === this){
          console.log('menu: self generated uristack change');
          return;
        }else{
          this.changeUri();
        }
      },
      
      updateTopMenu: function(ui_resource_id) {
        
        $('ul.topnav > li').each(function(index){
          var id = $(this).attr('id');
          if( id != ui_resource_id ){
            $(this).removeClass('active');
          }else{
            $(this).addClass('active');
          }
        });        
        if(_.isUndefined(ui_resource_id)){
          $('ul.topnav > #home').addClass('active');
        }
      },

      changeUri: function() {
        var uriStack = appModel.get('uriStack');
        var ui_resource_id = uriStack[0];
        console.log('got ui_resource: ' + ui_resource_id);
        
        this.updateTopMenu(ui_resource_id);

        this.$('li').removeClass('active');
        var menus = appModel.getMenu();
        this.collapseAll(menus);

        if(_.isEmpty(uriStack)){
          this.render();
          return; // Home, for now
        }

        var found_menus = this.find_submenu_path(menus, ui_resource_id);
        if(_.isUndefined(found_menus)){
          console.log('unknown submenu: ' + ui_resource_id);
          return;
        }else{
          _.each(found_menus, function(menu){
            if (menu['expanded'] == false ){
              menu['expanded'] = true;
            }
          });
          appModel.set({'menu':menus});
          this.render();
          var found_menu = found_menus[found_menus.length-1];
          if (_.has(found_menu,'key')){
            this.$('#' + found_menu['key']).addClass('active');
            this.ui_resource_id = found_menu['key'];
          } else {
            this.$('#' + ui_resource_id).addClass('active');
            this.ui_resource_id = ui_resource_id;
          }
        }
      },
      
      afterRender: function(){
        if(this.ui_resource_id){
          this.$('#' + this.ui_resource_id).addClass('active');
        }
        window.scrollTo(0, 0);
        
      },

      find_submenu_path : function(menu, id){
        if( _.has(menu, id) ) return [menu[id]];
        else if (_.has(menu,'alias') && menu['alias'] == id) return menu;
        else if(_.has(menu, 'submenus')){
          var pairs = _.pairs(menu['submenus']);
          for(var i=0; i < pairs.length; i++){
            var pair = pairs[i];
            if(pair[0] == id ) return [pair[1]];
            else if (_.has(pair[1],'alias') && pair[1]['alias'] == id){
              pair[1]['key'] = pair[0];
              return [pair[1]];
            }
            else{
              if(_.has(pair[1], 'submenus')){
                var temp = this.find_submenu_path(pair[1],id);
                if( _.isObject(temp)) return _.flatten([pair[1],temp]);
              }
            }
          }
        }
      },

//      find_submenu : function(menu, id){
//        if( _.has(menu, id) ) return menu[id];
//        else if(_.has(menu, 'submenus')){
//          var pairs = _.pairs(menu['submenus']);
//          for(var i=0; i < pairs.length; i++){
//            var pair = pairs[i];
//            if(pair[0] == id ) return pair[1];
//            else{
//              if(_.has(pair[1], 'submenus')){
//                var temp = this.find_submenu(pair[1],id);
//                if( _.isObject(temp)) return temp;
//              }
//            }
//          }
//        }
//      },

      menuClick: function(event){
        var self = this;
        event.preventDefault();

        var ui_resource_id = event.currentTarget.id;
        console.log('menu click: ' + ui_resource_id);
        
        if(appModel.isPagePending()){
          appModel.requestPageChange({
            immediate: function(){
              self._menuAction(ui_resource_id)
            },
            ok: function(){
              self._menuAction(ui_resource_id, true)
            }
          });
        }else{
          self._menuAction(ui_resource_id, false)
        }
        
      },
      
//      _menuActionOld: function(ui_resource_id, isImmediate){
//        if (ui_resource_id==='home') {
//          appModel.setUriStack([]);
//          return;
//        }
//
//        var menus = appModel.getMenu();
//        var menu = this.find_submenu(menus, ui_resource_id);
//        if(_.isUndefined(menu)){
//          window.alert('unknown submenu: ' + ui_resource_id);
//          return;
//        }
//
//        // if menu doesn't have an "expanded" flag, then just do it's action
//        if( ! _.has(menu, 'expanded') || isImmediate ){
//          appModel.setUriStack([ui_resource_id]);
//        }else{
//          // first click on a menu item expands it
//          // second click on a menu item will cause its action, if it is expanded
//          if( menu['expanded'] == false ){
//              menu['expanded'] = true;
//              appModel.set({'menu':menus});
//              this.render();
//          }else if( menu['expanded'] == true ){
//              menu['expanded'] = false;
//              appModel.set({'menu':menus});
//              this.render();
//
//              if( ui_resource_id == 'reports'){
//                appModel.setUriStack([]);
//                return;
//              } else if (_.isEmpty(menu['view'])){
//                return;
//              }else{
//                appModel.setUriStack([ui_resource_id]);
//              }
//          }
//        }
//        this.$('li').removeClass('active');
//        // also the top navbar
//        $('.nav').children('li').removeClass('active');
//        this.$('#' + ui_resource_id).addClass('active');
//        this.updateTopMenu(ui_resource_id);
//
//        // Clear out error messages after navigating away from page
//        appModel.unset('messages');
//
//      },
      
      _menuAction: function(ui_resource_id, isImmediate){
        var self = this;
        if (ui_resource_id==='home') {
          appModel.setUriStack([]);
          return;
        }
        var menus = appModel.getMenu();
        var found_menus = this.find_submenu_path(menus, ui_resource_id);
        if(_.isUndefined(found_menus)){
          window.alert('unknown submenu: ' + ui_resource_id);
          return;
        }else{
          // first click on a menu item expands it
          // second click on a menu item will cause its action, if it is expanded
          
          var found_menu = found_menus[found_menus.length-1];
          if( ! _.has(found_menu, 'expanded') || isImmediate ){
            self.collapseAll(menus);
            _.each(found_menus, function(menu){
              if (_.has(menu,'expanded') && menu['expanded'] == false ){
                menu['expanded'] = true;
              }
            });
            appModel.setUriStack([ui_resource_id]);
          }else if( found_menu['expanded'] == false ){
            _.each(found_menus, function(menu){
              if (_.has(menu,'expanded') && menu['expanded'] == false ){
                menu['expanded'] = true;
              }
            });
            appModel.set({'menu':menus});
            this.render();
          } else if (found_menu['expanded'] == true ){
            found_menu['expanded'] = false;
            appModel.set({'menu':menus});
            this.render();
            if( ui_resource_id == 'reports'){
              appModel.setUriStack([]);
              return;
            } else if (_.isEmpty(found_menu['view'])){
              return;
            }else{
              appModel.setUriStack([ui_resource_id]);
            }
          }
        }
        
        this.$('li').removeClass('active');
        // also the top navbar
        $('.nav').children('li').removeClass('active');
        this.$('#' + ui_resource_id).addClass('active');
        this.updateTopMenu(ui_resource_id);

        // Clear out error messages after navigating away from page
        appModel.unset('messages');

      },
      
      collapseAll: function(menu) {
        var self = this;
        if (_.result(menu, 'expanded', false) == true){
          menu['expanded'] = false;
        }
        var submenus = _.result(menu,'submenus');
        if (submenus){
          _.each(_.values(submenus), function(menu){
            self.collapseAll(menu);
          });
        }
      }

      
      
   
    });

    return MenuView;
});