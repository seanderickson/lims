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
//        _.bindAll(this,'_menuClick');
      },
      
      template: _.template(menuTemplate),
      
      serialize: function(){
        return {
          menu: appModel.get('menu'),
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
        
        if (ui_resource_id==='home') {
          appModel.setUriStack([]);
//          appModel.set({uriStack: []});
          return;
        }

        var menus = appModel.get('menu');

        var menu = this.find_submenu(menus, ui_resource_id);
        if(_.isUndefined(menu)){
          window.alert('unknown submenu: ' + ui_resource_id);
          return;
        }

        // if menu doesn't have an "expanded" flag, then just do it's action
        if( ! _.has(menu, 'expanded')){
//          appModel.set({ uriStack: [ui_resource_id] });
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

//              appModel.set({ uriStack: [ui_resource_id] });
              appModel.setUriStack([ui_resource_id]);
          }
        }
        this.$('li').removeClass('active');
        this.$('#' + ui_resource_id).addClass('active');
        this.updateTopMenu(ui_resource_id);
      },
      
//      _setUriStack: function(value){
//        var self = this;
//        var _deferred = function(){
//          appModel.set({ uriStack: value });
//        };
//        if(appModel.get('modal')){
//          self.showModal(_deferred);
//        }else{
//          _deferred();
//        }
//      },
//      
//      showModal: function(callback, initialEvent){
//        var self = this;
//        console.log('showModal: ' + JSON.stringify(this.model));
//        var modalDialog = new Backbone.View({
//            el: _.template(modalOkCancel, { 
//              body: "Pending changes in the page, continue or cancel and finish with changes?",
//              title: "Please confirm" } ),
//            events: {
//                'click #modal-cancel':function(event) {
//                    console.log('cancel button click event, '); 
//                    event.preventDefault();
//                    $('#modal').modal('hide'); 
//                },
//                'click #modal-ok':function(event) {
//                    console.log('ok button click event, '); 
//                    event.preventDefault();
//                    $('#modal').modal('hide');
////                    self.$el.empty();
////                    self.trigger('remove');
//                    appModel.set({'modal': false});
//                    callback(initialEvent);
//                }
//            },
//        });
//        modalDialog.render();
//        modalDialog.$el.find('#modal-ok').html('Continue');
//        $('#modal').empty();
//        $('#modal').html(modalDialog.$el);
//        $('#modal').modal({show:true, backdrop: 'static'});
//      }      
      

    });

    return MenuView;
});