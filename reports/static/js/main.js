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
    lunr: 'libs/lunr',
    text: 'libs/text',
    router: 'router',
    app_state: 'models/app_state',
    iccbl_backgrid: 'iccbl-backgrid',

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
      // Note, object naming rules won't allow the dot separator, as in "Backbone.PageableCollection"
      // so the PageableCollection must be tied back in to the Backbone object with each module it is imported into
      // see http://localhost:8000/reports/#list/screensaveruser/search/last_name=Y/order_by/first_name/rpp/25/page/2
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
    }
  }
});
require([
    'jquery',
    'underscore',
    'backbone',
    'models/app_state',
    'views/app_view',
    'router'
  ],
    function($, _, Backbone, AppModel, AppView, AppRouter ) {
        Backbone.View.prototype.close = function(){
          console.log('close called for' + this ); // not so informative
          this.remove();
          this.unbind();
          if (this.onClose){
            this.onClose();
          }
        };

        var appModel = new AppModel();
        var appView = new AppView({ model: appModel });
        var appRouter = new AppRouter({ model: appModel });
        appView.render();
    });