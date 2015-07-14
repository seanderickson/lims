require.config({
  // TODO: how to use the newest version of underscore.js? current version is 1.8.x
  // but we are using 1.5.2
  
  // Note: baseUrl is the current working directory.
  
  paths: {
    jquery: 'libs/jquery',
    underscore: 'libs/underscore',
    backbone: 'libs/backbone',
    'backbone.paginator': 'libs/backbone.paginator',
    backgrid: 'libs/backgrid',
    backgrid_filter: 'libs/backgrid-filter',
    backgrid_paginator: 'libs/backgrid-paginator',
    backbone_stickit: 'libs/backbone.stickit',
    backbone_modelbinder: 'libs/Backbone.ModelBinder',
    backbone_forms: 'libs/backbone-forms', // TODO: evaluating vs. backbone.stickit
    layoutmanager: 'libs/backbone.layoutmanager',
    bootstrap: 'libs/bootstrap',
    lunr: 'libs/lunr',
    text: 'libs/text',
    iccbl_backgrid: 'iccbl-backgrid',
    app_state: 'models/app_state',
    router: 'router',
    multiselect: 'libs/multiselect',
    quicksearch: 'libs/jquery.quicksearch'
  },
  
  shim: {
    underscore: {
      exports: "_"
    },
    backbone: {
      deps: ['underscore', 'jquery'],
      exports: 'Backbone'
    },
    'backbone.paginator': {
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
    },
    app_state: {
      deps: ['iccbl_backgrid']
    }
    
  }
  // see https://github.com/amdjs/amdjs-api/wiki/Common-Config
});