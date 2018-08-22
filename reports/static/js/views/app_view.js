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
    'views/backgroundJobPanel',
    'views/search_box',
    'templates/app_layout.html'
], function($, _, Backbone, layoutmanager, Iccbl, appModel, MenuView, 
            ContentView, MessageView, BackgroundJobPanel, SearchView, layout) {

    var AppView = Backbone.Layout.extend({
      el: '#application_div',
      
      initialize: function(args) {
        this.listenTo(appModel, 'change:messages', this.setMessages);
        
        _.bindAll(this,'setMessages','setBackgroundJobs');
      },
      
      views: {
        "#menu": new MenuView(),
        "#container": new ContentView({model: appModel, property: 'uriStack'}),
        "#search_box": new SearchView()
      },
      
      setMessages: function() {
        var messages = appModel.get('messages');
        console.log('set messages', messages);
        if(!_.isEmpty(messages)){
          this.messageView = new MessageView({model: appModel});
          this.setView("#messages", this.messageView).render();
        }else if(this.messageView){
          this.messageView.remove();
        }
      },
      
      setBackgroundJobs: function() {
        var jobs = appModel.getJobCollection();
        console.log('set background jobs', jobs);
        
        $jobsEl = this.$('#background_jobs');
        $jobsEl.empty();
        jobs.each(function(job){
          console.log('process job', job);
          var jobView = new BackgroundJobPanel({model: job});
          $jobsEl.append(jobView.render().$el)
        });
        
      },
      
      afterRender: function() {
        var self = this;
        function postRender() {
          if (appModel.hasPermission('screen','write')){
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
            $('#additional_buttons_box').append(addScreenButton);
          }
          if (appModel.hasPermission('user','write')){
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
            $('#additional_buttons_box').append(addUserButton);
          }
          if (appModel.hasPermission('library','write')){
            var addLibraryButton = $([
              '<a class="btn btn-default btn-sm pull-down" ',
                'role="button" id="add_library_button" href="#">',
                'Add Library</a>'
              ].join(''));
            addLibraryButton.click(function(e){
              e.preventDefault();
              var route = 'library/+add';
              appModel.router.navigate(route, {trigger: true});
            });
            $('#additional_buttons_box').append(addLibraryButton);
          }
          self.listenTo(
            appModel.getJobCollection(), 'add remove', self.setBackgroundJobs);

          $('#about').click(function(e){
            e.preventDefault();
            console.log('about click...');
            appModel.setUriStack(['about']);
          });
          $('#home').click(function(e){
            e.preventDefault();
            console.log('home click...');
            appModel.setUriStack(['home']);
          });
          
        
        };
        
        postRender();
      },
      
      template: _.template(layout)

    });

    return AppView;
});

