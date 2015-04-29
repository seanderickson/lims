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
        console.log('initialize menu.js');
        this.listenTo(appModel, 'change:uriStack', this.uriStackChange);
      },
      
      template: _.template(menuTemplate),
      
      serialize: function(){
        return {
          menu: appModel.getMenu(),
          ui_resources: appModel.get('ui_resources')
        };          
      },
      
      /**
       * Backbone.Model change event handler
       * @param options.source = the event source triggering view
       */
      uriStackChange: function(model, val, options) {
        if(options.source === this){
          console.log('self generated uristack change');
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
        if(_.isEmpty(uriStack)) return; // Home, for now

        var menus = appModel.get('menu');
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
        }
        this.ui_resource_id = ui_resource_id;
        this.$('#' + ui_resource_id).addClass('active');
      },
      
      afterRender: function(){
        if(this.ui_resource_id){
          this.$('#' + this.ui_resource_id).addClass('active');
        }
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

        // TODO: 20140618: test remove messages
        appModel.unset('messages');
        
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
      
      _menuAction: function(ui_resource_id, isImmediate){
        if (ui_resource_id==='home') {
          appModel.setUriStack([]);
          return;
        }

        var menus = appModel.get('menu');

        var menu = this.find_submenu(menus, ui_resource_id);
        if(_.isUndefined(menu)){
          window.alert('unknown submenu: ' + ui_resource_id);
          return;
        }

        // if menu doesn't have an "expanded" flag, then just do it's action
        if( ! _.has(menu, 'expanded') || isImmediate ){
          appModel.setUriStack([ui_resource_id]);
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

              appModel.setUriStack([ui_resource_id]);
          }
        }
        this.$('li').removeClass('active');
        // also the top navbar
        $('.nav').children('li').removeClass('active');
        this.$('#' + ui_resource_id).addClass('active');
        this.updateTopMenu(ui_resource_id);
      }
   
    });

    return MenuView;
});