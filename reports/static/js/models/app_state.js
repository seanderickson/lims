define([
  'jquery',
  'underscore',
  'backbone',
], function($, _, Backbone ){
    // TODO: feed AppState from the server
    var AppState = Backbone.Model.extend({
        defaults: {
            root_url: '/reports/',
            menu_items: [
                { id:'home', text:'Home'},
                { id:'home2', text:'Home2' },
                { id:'list_fieldinformation', text:'List FieldInformation' },
                { id:'list_screensaveruser', text:'Show Users' },
            ],// todo: make these into a collection of menu models?
            menu_actions: {
                home: {
                    route: '',
                    view: 'HomeView',
                    content_header: 'Welcome',
                },
                home2: {
                    route: 'hometest',
                    view: 'HomeView',
                },
                list_fieldinformation: {
                    route: 'list/fieldinformation',
                    view: 'ListView',
                    options: { type: 'fieldinformation', header_message: 'Field Meta Information' }
                },
                list_screensaveruser: {
                    route: 'list/screensaveruser',
                    view: 'ListView',
                    options: { type: 'screensaveruser', header_message: 'Screensaver Users' }
                }
            },
            menu_item: 'home',
            content_options: {},
            route: '',
            list_defaults: {
                type: null, // one of the keys for the list_endpoints
                page: 1,
                pageSize: 25,
                orderBy: null,
                searchBy: null,
            },
            list_endpoints: {
                fieldinformation: {
                    url_schema : '/reports/api/v1/fieldinformation/schema/',
                    url : '/reports/api/v1/fieldinformation/',
                },
                screensaveruser: {
                    url_schema : '/db/api/v1/screensaveruser/schema/' ,
                    url : '/db/api/v1/screensaveruser/'
                },
                screen: {
                    url_schema : '/db/api/v1/screen/schema/' ,
                    url : '/db/api/v1/screen/'
                }
            },
        },

        initialize : function() {
            console.log('AppState initialized');
        },
    });

    return AppState;
});