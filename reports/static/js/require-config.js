require.config({
  
  // Note: baseUrl is the current working directory.
  
  paths: {
    jquery: 'libs/jquery',
    underscore: 'libs/underscore',
    backbone: 'libs/backbone',
    'backbone.paginator' : 'libs/backbone.paginator',
    backgrid: 'libs/backgrid',
    backgrid_filter: 'libs/backgrid-filter',
    backgrid_paginator: 'libs/backgrid-paginator',
    backgrid_select_all: 'libs/backgrid-select-all',
    backbone_stickit: 'libs/backbone.stickit',
    backbone_modelbinder: 'libs/Backbone.ModelBinder',
    backbone_forms: 'libs/backbone-forms', 
    backbone_forms_list: 'libs/backbone-forms-list',
    layoutmanager: 'libs/backbone.layoutmanager',
    bootstrap: 'libs/bootstrap',
    lunr: 'libs/lunr',
    text: 'libs/text',
    router: 'router',
    app_state: 'models/app_state',
    iccbl_backgrid: 'iccbl-backgrid',
    multiselect: 'libs/multiselect',
    quicksearch: 'libs/jquery.quicksearch',
    'bootstrap-datepicker': 'libs/bootstrap-datepicker',
    chosen: 'libs/chosen'
  },
  shim: {
    underscore: {
      exports: "_"
    },
    backbone: {
      deps: ['underscore', 'jquery'],
      exports: 'Backbone'
    },
    'backbone.paginator' : {
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
    },
    'bootstrap-datepicker': {
      deps: ['jquery','bootstrap'],
      exports: 'bootstrap-datepicker'
    },
    backbone_stickit: { deps: ["backbone"] },
    
    backbone_forms: { deps: ["backbone"] }, 
    backbone_forms_list: { deps: ["backbone","backbone_forms"] }, 

    'backbone-associations': { deps: ["backbone"] },

    iccbl_backgrid: {
        exports: 'Iccbl'
    },
    app_state: {
      deps: ['iccbl_backgrid']
    }
  }
});