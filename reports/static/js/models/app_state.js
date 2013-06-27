define([
  'jquery',
  'underscore',
  'backbone',
], function($, _, Backbone ){
    // TODO: feed AppState from the server
    var AppState = Backbone.Model.extend({
        defaults: {
            root_url: '/reports/',
            current_submodel: 'home',
            content_options: {},
            route: '',
            // menu_items: [
                // { id:'home', text:'Home'},
                // { id:'list_fieldinformation', text:'List FieldInformation' },
                // { id:'list_metainformation', text:'List MetaInformation' },
                // { id:'list_screensaveruser', text:'Show Users' },
            // ],// todo: make these into a collection of menu models?
            submodels: {
                home: {
                    type: 'home',
                    title: 'Home',
                    route: '',
                    view: 'HomeView',
                    content_header: 'Welcome',
                },

                fieldmetainformation: {
                    type: 'fieldmetainformation',
                    header_message: 'Field Meta Information',
                    title: 'Field Meta Information',
                    route: 'list/fieldmetainformation',
                    view: 'ListView',
                    url_schema : '/reports/api/v1/metahash/schema/',
                    url : '/reports/api/v1/metahash/',
                    searchBy : 'scope=metahash:fields',
                },

                vocabularies: {
                    type: 'vocabularies',
                    header_message: 'Application Vocabularies',
                    title: 'Application Vocabularies',
                    route: 'list/vocabularies',
                    view: 'ListView',
                    url_schema : '/reports/api/v1/vocabularies/schema/',
                    url : '/reports/api/v1/vocabularies/',
                },

                screensaveruser: {
                    type: 'screensaveruser',
                    header_message: 'Screensaver Users',
                    title: 'Users',
                    route: 'list/screensaveruser',
                    view: 'ListView',
                    url_schema : '/db/api/v1/screensaveruser/schema/' ,
                    url : '/db/api/v1/screensaveruser/'
                },
                screen: {
                    type: 'screen',
                    header_message: 'Screens',
                    title: 'Screens',
                    route: 'list/screen',
                    view: 'ListView',
                    url_schema : '/db/api/v1/screen/schema/' ,
                    url : '/db/api/v1/screen/'
                }
            },
            list_defaults: {
                type: null, // one of the keys for the list_endpoints
                page: 1,
                pageSize: 25,
                orderBy: null,
                searchBy: null,
            },
            // menu_actions: {
                // home: {
                    // route: '',
                    // view: 'HomeView',
                    // content_header: 'Welcome',
                // },
                // list_fieldinformation: {
                    // route: 'list/fieldinformation',
                    // view: 'ListView',
                    // options: { type: 'fieldinformation', header_message: 'Field Meta Information', title: 'Field Information Table' }
                // },
                // list_metainformation: {
                    // route: 'list/metainformation',
                    // view: 'ListView',
                    // options: { type: 'metainformation', header_message: 'Meta Information', title: 'Meta Information Table' }
                // },
                // list_screensaveruser: {
                    // route: 'list/screensaveruser',
                    // view: 'ListView',
                    // options: { type: 'screensaveruser', header_message: 'Screensaver Users', title: 'Screensaver User Table' }
                // }
            // },
        },

        initialize : function() {
            console.log('AppState initialized');
        },
    });

    return AppState;
});