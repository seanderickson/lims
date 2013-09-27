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

            // NOTE: apparently, change notifications are only detected for 'simple' objects (i.e. one level)
            current_view: 'home',
            current_resource_id: 'home',
            current_options: {},
            routing_options: {},
            current_scratch: {},
            current_details: {},

            // menu's can only be two levels deep (todo: recursive cool menu thing)
            // - the menu tree here is composed of ui_resources, defined below.
            // - the topographical structure here represents the visual structure in the app.
            menu: {
                view: 'home',
                submenus:{
                    'library': {
                        expanded: false,
                        view: 'list',
                        submenus: {
                            'smallmoleculelibrary':{
                                view: 'list'
                            },
                            'rnalibrary':{
                                view: 'list'
                            },
                        }
                    },
                    'screensaveruser': {
                        expanded: false,
                        view: 'list',
                        submenus: {
                        'screeners':{
                            view: 'list'
                        },
                        'staff':{
                            view: 'list'
                        },
                        }
                    },
                    'screen': {
                        view: 'list',
                        expanded: false,
                        submenus: {
                        'small_molecule_screens':{
                            view: 'list',
                        },
                        'rnai_screens':{
                            view: 'list',
                        },
                        }
                    },
                    'admin': {
                        expanded: false,
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

            ui_resources: {   /* TODO: *all* of this will come from the resourceResource */
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
                    options: { search: {'scope': 'fields:metahash'}},
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
                    description: 'View user information'
                },

                screeners: {
                    header_message: 'Screening Users',
                    title: 'Screeners',
                    route: 'list/screeners',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'screensaveruser',
                    url_root: '/db/api/v1',
                    description: 'View user information',
                    options: { search: {'screeningroomuser__isnull': 'False'} }
                },

                staff: {
                    header_message: 'Staff',
                    title: 'Staff Users',
                    route: 'list/staff',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'screensaveruser',
                    url_root: '/db/api/v1',
                    description: 'View user information',
                    options: { search: {'administratoruser__isnull': 'False'} }
                },
                screen: {
                    header_message: 'All screens (Small Molecule and RNAi)',
                    title: 'Screens',
                    route: 'list/screen',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'screen',
                    url_root: '/db/api/v1',
                    description: 'View screen information'
                },
                small_molecule_screens: {
                    header_message: 'Small Molecule Screens',
                    title: 'Small Molecule',
                    route: 'list/small_molecule_screens',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'screen',
                    url_root: '/db/api/v1',
                    description: 'View small molecule screen information',
                    options: { search: { screen_type: 'Small Molecule'} }
                },
                rnai_screens: {
                    header_message: 'All screens (Small Molecule and RNAi)',
                    title: 'RNAi',
                    route: 'list/rnai_screens',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'screen',
                    url_root: '/db/api/v1',
                    description: 'View rnai screen information',
                    options: { search: { screen_type: 'RNAi'} }
                },
                library: {
                    header_message: 'All libraries (Small Molecule and RNAi)',
                    title: 'Libraries',
                    route: 'list/library',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'library',
                    url_root: '/db/api/v1',
                    description: 'View library information'
                },
                smallmoleculelibrary: {
                    header_message: 'Small Molecule Libraries',
                    title: 'Small Molecule',
                    route: 'list/smallmoleculelibrary',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'library',
                    url_root: '/db/api/v1',
                    description: 'View Small Molecule Library information',
                    options: { search: { screen_type: 'Small Molecule'} }
                },
                rnalibrary: {
                    header_message: 'RNAi Libraries',
                    title: 'RNAi',
                    route: 'list/rnalibrary',
                    list_view: 'ListView',
                    detail_view: 'DetailView',
                    api_resource: 'library',
                    url_root: '/db/api/v1',
                    description: 'View RNAi library information',
                    options: { search: { screen_type: 'RNAi'} }
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