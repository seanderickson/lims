define([
  'jquery',
  'underscore',
  'backbone',
], function($, _, Backbone ){
    // TODO: feed AppState from the server
    var AppState = Backbone.Model.extend({
        defaults: {

            root_url: '/reports',  // used for the backbone history
            api_root_url: '/reports/api/v1',

            // current_ui_resource_id: 'home',
            // current_view: 'home',
            // current_route_options: {},
            // current_route_update: {},
            // content_options: {},
            // route: '',

            // NOTE: apparently, change notifications are only detected for 'simple' objects (i.e. one level)
            current_view: 'home',
            current_resource_id: 'home',
            current_options: {},
            routing_options: {},
            current_scratch: {},


            // menu's can only be two levels deep (todo: recursive cool menu thing)
            // - the menu tree here is composed of ui_resources, defined below.
            // - the topographical structure here represents the visual structure in the app.
            menu: {
                view: 'home',
                submenus:{
                'home': {
                    view: 'home'
                },
                'screensaveruser': {
                    view: 'list',
                    submenus: {
                    'screeners':{},
                    'staff':{},
                    }
                },
                'screen': {
                    view: 'list',
                    submenus: {
                    'small_molecule_screens':{},
                    'rnai_screens':{},
                    }
                },
                'admin': {
                    view: 'menu',  // TODO: not implemented; would show children menu items as links in the page
                    submenus: {
                        'metahash':{
                            view: 'list',
                        },
                        'resource':{
                            view: 'list',
                        },
                        'vocabularies':{
                            view: 'list',
                        },
                        'apilog':{
                            view: 'list',
                        },
                    }
                },
                }
            },

            ui_resources: {
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
                    api_resource: 'metahash',
                    url_root: '/reports/api/v1',
                    // url_schema : '/reports/api/v1/metahash/schema',
                    // url : '/reports/api/v1/metahash',
                    searchBy : 'scope=fields:metahash',
                    description: 'Control field settings'
                },

                resource: {
                    header_message: 'Define resources for display in the reporting interface',
                    title: 'Resource Information',
                    route: 'list/resource',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'resource',
                    url_root: '/reports/api/v1',
                    // url_schema : '/reports/api/v1/resource/schema',
                    // url : '/reports/api/v1/resource',
                    description: 'Control resource information'
                },

                vocabularies: {
                    header_message: 'Define controlled vocabularies',
                    title: 'Application Vocabularies',
                    route: 'list/vocabularies',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'vocabularies',
                    url_root: '/reports/api/v1',
                    // url_schema : '/reports/api/v1/vocabularies/schema',
                    // url : '/reports/api/v1/vocabularies',
                    description: 'Enter controlled vocabularies'
                },

                apilog: {
                    header_message: 'View change logs',
                    title: 'Change logs',
                    route: 'list/apilog',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'apilog',
                    url_root: '/reports/api/v1',
                    // url_schema : '/reports/api/v1/apilog/schema',
                    // url : '/reports/api/v1/apilog',
                    description: 'Change logs'
                },

                screensaveruser: {
                    header_message: 'All users (Screeners and Staff)',
                    title: 'Screensaver Users',
                    route: 'list/screensaveruser',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'screensaveruser',
                    url_root: '/db/api/v1',
                    // url_schema : '/db/api/v1/screensaveruser/schema' ,
                    // url : '/db/api/v1/screensaveruser',
                    description: 'View user information'
                },
                screen: {
                    header_message: 'All screens (Small Molecule and RNAi)',
                    title: 'Screens',
                    route: 'list/screen',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'screen',
                    url_root: '/db/api/v1',
                    // url_schema : '/db/api/v1/screen/schema' ,
                    // url : '/db/api/v1/screen',
                    description: 'View screen information'
                }
            },
            list_defaults: {
                page: 1,
                rpp: 25,
                order: null,
                search: null,
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