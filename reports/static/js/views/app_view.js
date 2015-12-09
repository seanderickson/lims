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
        
        var addScreenButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="add_screen_button" href="#">',
            'Add Screen</a>'
          ].join(''));
        addScreenButton.click(function(e){
          e.preventDefault();
          var route = 'screen/+add';
          appModel.router.navigate(route, {trigger: true});
        });
        var addUserButton = $([
          '<a class="btn btn-default btn-sm pull-down" ',
            'role="button" id="add_user_button" href="#">',
            'Add User</a>'
          ].join(''));
        addUserButton.click(function(e){
          e.preventDefault();
          var route = 'screensaveruser/+add';
          appModel.router.navigate(route, {trigger: true});
        });
        $('#additional_buttons_box').append(addScreenButton);
        $('#additional_buttons_box').append(addUserButton);
        
        
        
      },
      
      template: _.template(layout)

    });

    return AppView;
});

