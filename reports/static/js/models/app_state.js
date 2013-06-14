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
                { id:'home', text:'Home' },
                { id:'home2', text:'Home2' },
                { id:'list_fieldinformation', text:'List FieldInformation' },
                { id:'list_screensaveruser', text:'Show Users' },
            ],// todo: make these into a collection of menu models?
            menu_actions: {
                home: {
                    route: '',
                    view: 'HomeView',
                },
                home2: {
                    route: 'hometest',
                    view: 'HomeView',
                },
                list_fieldinformation: {
                    route: 'list/fieldinformation',
                    view: 'ListView',
                    options: { type: 'fieldinformation' }
                },
                list_screensaveruser: {
                    route: 'list/screensaveruser',
                    view: 'ListView',
                    options: { type: 'screensaveruser' }
                }
            },
            menu_item: 'home',
            route_content_options: {},
            route: '',
            list_defaults: {
                type: null, // one of the keys for the list_endpoints
                page: 1,
                pageSize: 25,
                sort: {
                    sortKey: null, // a valid column
                    order: null, // 1 === descending, -1 asc
                    direction: null, // ascending or descending
                },
                search: {
                    searchBy: null, // the full search term, i.e. "name=foo"
                    searchColumn: null, // some searchable field of the api endpoint
                    searchTerm: null // the term to search for
                }
            },
            list_endpoints: {
                fieldinformation: {
                    url_schema : '/reports/api/v1/fieldinformation/schema/',
                    url : '/reports/api/v1/fieldinformation/?format=json',
                },
                screensaveruser: {
                    url_schema : '/db/api/v1/screensaveruser/schema/' ,
                    url : '/db/api/v1/screensaveruser/?format=json'
                },
                screen: {
                    url_schema : '/db/api/v1/screen/schema/' ,
                    url : '/db/api/v1/screen/?format=json'
                }
            },
        },

        initialize : function() {
            console.log('AppState initialized');
        },
    });

    return AppState;
});