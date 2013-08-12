define([
  'jquery',
  'underscore',
  'backbone',
], function($, _, Backbone ){
    // TODO: feed AppState from the server
    var AppState = Backbone.Model.extend({
        defaults: {
            root_url: '/reports/',
            current_submodel_id: 'home',
            content_options: {},
            route: '',
            // menu's can only be two levels deep (todo: recursive cool menu thing)
            menu: {
                view: 'home_view',
                submenus:{
                'screensaveruser': {
                    view: 'list_view',
                    submenus: {
                    'screeners':{},
                    'staff':{},
                    }
                },
                'screen': {
                    view: 'list_view',
                    submenus: {
                    'small_molecule_screens':{},
                    'rnai_screens':{},
                    }
                },
                'admin': {
                    view: 'menu_view',  // TODO: not implemented; would show children menu items as links in the page
                    submenus: {
                        'metahash':{
                            view: 'list_view',
                        },
                        'resource':{
                            view: 'list_view',
                        },
                        'vocabularies':{
                            view: 'list_view',
                        },
                        'apilog':{
                            view: 'list_view',
                        },
                    }
                },
                }
            },

            submodels: {
                admin: {
                    title: 'Admin',
                    route: '',
                    view: 'HomeView',
                    content_header: 'Welcome',
                    description: 'Perform admin activities'
                },
                home: {
                    title: 'Screensaver Reporting',
                    route: '',
                    view: 'HomeView',
                    content_header: 'Welcome',
                    description: 'Menu starting point'
                },

                metahash: {
                    header_message: 'Define fields for display on detail and list views',
                    title: 'Field Information',
                    route: 'list/metahash',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    url_schema : '/reports/api/v1/metahash/schema',
                    url : '/reports/api/v1/metahash',
                    searchBy : 'scope=fields:metahash',
                    description: 'Control field settings'
                },

                resource: {
                    header_message: 'Define resources for display in the reporting interface',
                    title: 'Resource Information',
                    route: 'list/resource',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    url_schema : '/reports/api/v1/resource/schema',
                    url : '/reports/api/v1/resource',
                    description: 'Control resource information'
                },

                vocabularies: {
                    header_message: 'Define controlled vocabularies',
                    title: 'Application Vocabularies',
                    route: 'list/vocabularies',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    url_schema : '/reports/api/v1/vocabularies/schema',
                    url : '/reports/api/v1/vocabularies',
                    description: 'Enter controlled vocabularies'
                },

                apilog: {
                    header_message: 'View change logs',
                    title: 'Change logs',
                    route: 'list/apilog',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    url_schema : '/reports/api/v1/apilog/schema',
                    url : '/reports/api/v1/apilog',
                    description: 'Change logs'
                },

                screensaveruser: {
                    header_message: 'All users (Screeners and Staff)',
                    title: 'Screensaver Users',
                    route: 'list/screensaveruser',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    url_schema : '/db/api/v1/screensaveruser/schema' ,
                    url : '/db/api/v1/screensaveruser',
                    description: 'View user information'
                },
                screen: {
                    header_message: 'All screens (Small Molecule and RNAi)',
                    title: 'Screens',
                    route: 'list/screen',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    url_schema : '/db/api/v1/screen/schema' ,
                    url : '/db/api/v1/screen',
                    description: 'View screen information'
                }
            },
            list_defaults: {
                type: null, // one of the keys for the list_endpoints
                page: 1,
                pageSize: 25,
                orderBy: null,
                searchBy: null,
            },
            detail_defaults: {


            },

        },

        initialize : function() {
            console.log('AppState initialized');
        },
    });

    return AppState;
});