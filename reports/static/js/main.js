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
//    backbone_modelbinder: 'libs/Backbone.ModelBinder',
    backbone_stickit: 'libs/backbone.stickit',
    bootstrap: 'libs/bootstrap',
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
    bootstrap: {
        deps: ['jquery'],
        // exports: '$.fn.modal'
    },
    backbone_stickit: { deps: ["backbone"] }
    , iccbl_backgrid: {
        exports: 'Iccbl'
    }
    //    deps: ["backbone", "underscore", "jquery", "backgrid", "backgrid_filter", "backgrid_paginator"]
    //}
  }
});
require([
    'jquery',
    'underscore',
    'backbone',
    'models/app_state',
    'views/app_view',
    'router',
    'bootstrap' // Note: Bootstrap does not return an object; it modifies the Jquery object with new methods
  ],
    function($, _, Backbone, AppModel, AppView, AppRouter ) {

        // Augment the view prototype with this close utility function to prevent memory leaks
        // See: http://lostechies.com/derickbailey/2011/09/15/zombies-run-managing-page-transitions-in-backbone-apps/
        // Also, consider using the Marionette extension
        Backbone.View.prototype.close = function(){
          this.remove();
          this.unbind();
          if (this.onClose){
            this.onClose();
          }
        };

        var appModel = new AppModel();
        var appRouter = new AppRouter({ model: appModel });
        var appView = new AppView({ model: appModel });
        appView.render();
        Backbone.history = Backbone.history || new Backbone.History({});
        Backbone.history.start({ pushState: false, root: appModel.get('root_url') });

        //Backbone.history.start({ pushState: false, root: appModel.get('root_url') });  // TODO: root should not be necessary if not using pushState
    });