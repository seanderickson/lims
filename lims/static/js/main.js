/**
 * Application loading script for the Iccbl-lims app.
 */

requirejs(['require-config'], 
    function() { // First, load require.js and require-config
  
  require([ // now load application code
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'models/app_state',
    'views/app_view',
    'router',
    // NOTE: browser security restrictions prevent loading this from file:///
    'text!models/ui_resources_fixture.js', 
    'text!models/menu_fixture.js', 
    // TODO: verify this: Bootstrap does not return an object; it modifies the 
    // Jquery object with new methods
    'bootstrap'
  ],
      function($, _, Backbone, Iccbl, appModel, AppView, AppRouter, 
               ui_resources_raw, menu_raw ) {
    
    console.log('init screensaver/reports...')
  
    // Augment the view prototype to prevent memory leaks - 
    // See: http://lostechies.com/derickbailey/2011/09/15/zombies-run-managing-page-transitions-in-backbone-apps/
    // Todo: is this still needed with the Backbone.Layout extension managing all the views?
    Backbone.View.prototype.close = function(){
      this.remove();
      this.unbind();
      if (this.onClose){
        this.onClose();
      }
    };
  
    if(_.isUndefined(window.logged_in) || window.logged_in != 'True' ){
      window.location='/accounts/login/?next=' + 
        window.location.pathname + window.location.hash;
      return;
    }
    
    appModel.set('ui_resources', JSON.parse(ui_resources_raw));
    appModel.set('menu', JSON.parse(menu_raw));

    var appRouter = appModel.router = new AppRouter({ model: appModel });
    var appView = new AppView({ model: appModel },{ router: appRouter});
  
    appModel.start(function(){
      console.log('Render application')
      appView.$el.appendTo("#application_div")
      appView.render();
          
      Backbone.history = Backbone.history || new Backbone.History({});
      Backbone.history.start({ pushState: false, root: appModel.get('root_url') });
    });
    
    // Set the document title
    Backbone.history.on('route', function(router, route, params) {
      var title = _.reduce(
          params,
          function(memo, item){
            if(item){
              if (memo !== ' ') memo += ', ';
              memo += item;
            }
            return memo ;
          }, ' ');              
      
        document.title = 'Screensaver LIMS' + ':' + title;
     }, this);    
   
  });

});