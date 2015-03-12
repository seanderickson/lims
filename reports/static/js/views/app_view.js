define([
    'jquery',
    'underscore',
    'backbone',
    'layoutmanager',
    'iccbl_backgrid',
    'models/app_state',
    'views/menu',
    'views/content',
    'views/message',
    'views/search_box',
    'text!templates/app_layout.html'
], function($, _, Backbone, layoutmanager, Iccbl, appModel, MenuView, 
            ContentView, MessageView, SearchView, layout) {

    var AppView = Backbone.Layout.extend({
      el: '#application_div',
      
      initialize: function(args) {
        this.listenTo(appModel, 'change:messages', this.setMessages);
        _.bindAll(this,'setMessages');
      },
      
      views: {
        "#menu": new MenuView(),
        "#container": new ContentView({model: appModel, property: 'uriStack'}),
      },
      
      setMessages: function() {
        var messages = appModel.get('messages');
        if(!_.isEmpty(messages)){
          this.messageView = new MessageView({model: appModel});
          this.setView("#messages", this.messageView).render();
        }else if(this.messageView){
          this.messageView.remove();
        }
      },
      
      afterRender: function() {
        this.searchView = new SearchView(),
        Backbone.Layout.setupView(this.searchView);
        this.setView("#search_box", this.searchView ).render();
      },
      
      template: _.template(layout)

    });

    return AppView;
});

