require.config({
  // baseUrl: "/js/",
  paths: {
    jquery: 'libs/jquery',
    underscore: 'libs/underscore',
    backbone: 'libs/backbone',
    backbone_pageable: 'libs/backbone-pageable',
    backgrid: 'libs/backgrid',
    backgrid_filter: 'libs/backgrid-filter',
    backgrid_paginator: 'libs/backgrid-paginator',
    backgrid_select_all: 'libs/backgrid-select-all',
    backbone_stickit: 'libs/backbone.stickit',
    backbone_modelbinder: 'libs/Backbone.ModelBinder',
    backbone_forms: 'libs/backbone-forms', // TODO: evaluating vs. backbone.stickit
    layoutmanager: 'libs/backbone.layoutmanager',
    bootstrap: 'libs/bootstrap',
    lunr: 'libs/lunr',
    text: 'libs/text',
    router: 'router',
    app_state: 'models/app_state',
    iccbl_backgrid: 'iccbl-backgrid'

  },
  shim: {
    underscore: {
      exports: "_"
    },
    backbone: {
      deps: ['underscore', 'jquery'],
      exports: 'Backbone'
    },
    backbone_pageable: {
      deps: ['backbone', 'underscore', 'jquery'],
      exports: 'BackbonePageableCollection'
      // Note, object naming rules won't allow the dot separator, 
      // as in "Backbone.PageableCollection"
      // so the PageableCollection must be tied back in to the Backbone object 
      // with each module it is imported into
    },
    backgrid: {
      deps: ['backbone', 'underscore', 'jquery'],
      exports: 'Backgrid'
    },
    backgrid_filter: {
      deps: ['backgrid', 'backbone', 'underscore', 'jquery'],
      exports: 'BackgridFilter'
    },
    backgrid_paginator: {
      deps: ['backgrid', 'backbone', 'underscore', 'jquery'],
      exports: 'BackgridPaginator'
    },
    backgrid_select_all: {
        deps: ['backgrid', 'backbone', 'underscore', 'jquery'],
      exports: 'BackgridSelectAll'
      },
    bootstrap: {
        deps: ['jquery'],
        // exports: '$.fn.modal'
    },
    backbone_stickit: { deps: ["backbone"] },
    
    backbone_forms: { deps: ["backbone"] }, // TODO: evaluating vs. backbone_stickit
    
    'backbone-associations': { deps: ["backbone"] },

    iccbl_backgrid: {
        exports: 'Iccbl'
    }
  }
});
require([
  'jquery',
  'underscore',
  'backbone',
  'iccbl_backgrid',
  'models/app_state',
  'views/app_view',
  'router',
  'bootstrap' // Note: Bootstrap does not return an object; it modifies the Jquery object with new methods
],
function($, _, Backbone, Iccbl, appModel, AppView, AppRouter ) {
  console.log('init screensaver/reports...')
  // Augment the view prototype with this close utility function to 
  // prevent memory leaks
  // See: http://lostechies.com/derickbailey/2011/09/15/zombies-run-managing-page-transitions-in-backbone-apps/
  // Also, consider using the Marionette extension
  Backbone.View.prototype.close = function(){
    this.remove();
    this.unbind();
    if (this.onClose){
      this.onClose();
    }
  };

  if(_.isUndefined(window.logged_in) || window.logged_in != 'True' ){
  	window.location='/accounts/login/?next=' + window.location.pathname + window.location.hash;
  }

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
    
      document.title = route + ':' + title;
   }, this);    
 
});