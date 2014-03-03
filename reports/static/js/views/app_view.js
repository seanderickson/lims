define([
    'jquery',
    'underscore',
    'backbone',
    'layoutmanager',
    'iccbl_backgrid',
    'models/app_state',
    'views/menu',
    'views/content',
    'text!templates/app_layout.html'
], function($, _, Backbone, layoutmanager, Iccbl, appModel, MenuView, 
            ContentView, layout) {

    var AppView = Backbone.Layout.extend({
      
      views: {
        "#menu": new MenuView(),
        "#container": new ContentView()
      },
      
      
      template: _.template(layout)


    });

    return AppView;
});