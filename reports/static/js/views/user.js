define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'views/detail_stickit',
    'views/list',
    'views/simple-list',
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, Iccbl, DetailView, ListView, SimpleListView, tabbedTemplate) {
    var UserView = Backbone.View.extend({
        events: {
//            'click .tabbable li': 'click_tab', // TODO: how to make this specific to this view? (it is also catching clicks on the table paginator)
            'click li': 'click_tab',
        },

        initialize : function(attributes, options) {
            this.options = options;
            var self = this;
            Iccbl.assert( !_.isUndefined(options.url_root), 'UserView: options.url_root must be defined');
            Iccbl.assert( !_.isUndefined(options.current_options), 'UserView: options.current_options must be defined');
            Iccbl.assert( !_.isUndefined(options.user_model), 'UserView: options.user_model must be defined');
            Iccbl.assert( !_.isUndefined(options.user_schema), 'UserView: options.user_schema must be defined');

            console.log('initialize UserView: current options: ' + JSON.stringify(options.current_options));

            _.bindAll(this, 'setDetail', 'setGroups', 'setPermissions');

            this._tabbed_resources = {
                detail: { description: 'User Details', title: 'User Details', invoke : this.setDetail },
                groups: { description: 'User Groups', title: 'User Groups', invoke : this.setGroups },
                permissions: { description: 'User Permissions', title: 'User Permissions', invoke : this.setPermissions },
            };

            var compiledTemplate = _.template( tabbedTemplate, { tab_resources: this._tabbed_resources, });

            this.$el.html(compiledTemplate);

            var tab = options.current_options['tab'];
            if(_.isUndefined(tab) || _.isEmpty(tab) ) tab = 'detail';

            this.listenTo(this.model, 'change:current_detail', this.change_current_detail);

            self.user = self.options.user_model;
            self.userSchema = self.options.user_schema;
            self._tabbed_resources[tab].invoke();
            self.$('#' + tab).addClass('active');

            console.log('initialized....');
        },

        change_current_detail: function(){
            var _current_options = this.model.get('current_detail');
            console.log('change_current_detail: ' + JSON.stringify(_current_options) );

            //TODO: will there be a case where we have to reset all the options?
            var tab = _current_options['tab'];
            if(this.options.current_options['tab'] != tab)
                if(!_.isUndefined(tab)) this.change_to_tab(tab);
        },

        setDetail: function(){
            var self = this;

            var detailView =
                new DetailView({ model: self.user },
                    {
                        schemaResult:self.userSchema,
                        router:self.options.router,
                        isEditMode: false
                    });
            self.detail = detailView;
            self.$('#tab_container').html(self.detail.render().el);
            self.$('#detail').addClass('active'); // first time not clicked so set manually
        },

        setGroups: function(){
            var self = this;

            if(_.isUndefined(this.groups)){
                var createGroups = function(schemaResult, collection){
                    var columns = Iccbl.createBackgridColModel(schemaResult.fields, Iccbl.MyHeaderCell);//, col_options );

                    // make a checkbox column that will be used to manage the group->user relationship
                    columns.unshift({
                        column: "isForUser",
                        name: "isForUser",
                        label: self.user.get('username') + " groups",
                        cell: "boolean"
                    });

                    // update the collection models with a new attribute tied to the checkbox column
                    collection.each(function(model){
                        var currentUsers = model.get('users');
                        var currentUserUri = self.user.get('resource_uri');
                        console.log('evaluate group with users: ' + model.get('name') + ': ' + JSON.stringify(currentUsers) + ', for user: ' + currentUserUri)
                        model.set('isForUser', _.contains(currentUsers, currentUserUri));
                        console.log('is for: ' + model.get('isForUser'));

                        // when the checkbox is selected, use a custom model event to update the group model relationship and refetch the user model
                        model.on('change:isForUser', function(model){
                            var currentUsers = model.get('users');
                            var currentUserUri = self.user.get('resource_uri');
                            if(model.get('isForUser')){
                                currentUsers.unshift(currentUserUri);
                            }else{
                                currentUsers = _.without(currentUsers, currentUserUri);
                            }
                            model.save({'users':currentUsers},{patch: true,
                                success: function(){
                                    self.user.fetch();
                                }
                            });
                        });
                    });
                    self.groups = new SimpleListView({ model: self.model},{ columns: columns,collection: collection });
                    $('#tab_container').html(self.groups.render().el);
                };

                var resource_url = this.options.url_root + '/usergroup';
                var schema_url =  resource_url + '/schema';
                var url = resource_url;
                Iccbl.getSchemaAndCollection(schema_url, url, createGroups);
            }else{
                self.groups.setElement(self.$('#tab_container')).render();
            }
            self.$('#groups').addClass('active'); // first time not clicked so set manually
        },

        setPermissions: function(){
            var self = this;

            if(_.isUndefined(this.permissions)){
                var createPermissions = function(schemaResult, collection){
                    var columns = Iccbl.createBackgridColModel(schemaResult.fields, Iccbl.MyHeaderCell);//, col_options );

                    // make a checkbox column that will be used to manage the user->permission relationship
                    columns.unshift({
                        column: "isForUser",
                        name: "isForUser",
                        label: self.user.get('username') + " permissions",
                        cell: "boolean"
                    });

                    // update the collection models with a new attribute tied to the checkbox column
                    collection.each(function(model){
                        var currentPerms = self.user.get('permissions');
                        var currentPermissionUri = model.get('resource_uri');
                        model.set('isForUser', _.contains(currentPerms, currentPermissionUri));

                        // when the checkbox is selected, use a custom model event to update the group model relationship and refetch the user model
                        model.on('change:isForUser', function(model){
                            var currentPerms = self.user.get('permissions');
                            var permissionUri = model.get('resource_uri');
                            if(model.get('isForUser')){
                                currentPerms.push(permissionUri);
                            }else{
                                currentPerms = _.without(currentPerms, permissionUri);
                            }
                            self.user.save({'permissions': currentPerms },{
                                patch: true,
                                success: function(){
                                    console.log('refetch the permission model');
                                    model.fetch();
                                }
                            });
                        });
                    });
                    self.permissions = new SimpleListView({ model: self.model},{ columns: columns,collection: collection });
                    $('#tab_container').html(self.permissions.render().el);
                };

                // get all of the permissions first, then checkoff the ones that are for this user
                var resource_url = this.options.url_root + '/permission';
                var schema_url =  resource_url + '/schema';
                var url = resource_url;
                Iccbl.getSchemaAndCollection(schema_url, url, createPermissions);
            }else{
                self.permissions.setElement(self.$('#tab_container')).render();
            }
            self.$('#permissions').addClass('active'); // first time not clicked so set manually
        },


        click_tab : function(event){
            event.preventDefault();

            // Block clicks from the wrong elements
            // TODO: how to make this specific to this view? (it is also catching clicks on the table paginator)
            var key = event.currentTarget.id;
            if(_.isEmpty(key)) return;

            console.log('tab click: ' + key + ', options: ' + JSON.stringify(this.options.current_options));
            this.change_to_tab(key);

            // notify the model listeners
            var _current_options = {'tab': key, 'key':Iccbl.getKey(this.options.current_options) };
            this.model.set({ current_options: {}, current_detail: {} },{ silent: true });
            this.model.set({
                current_options: _current_options,
                routing_options: {trigger: false, replace: false}
            }); // signal to the app_model that the current view has changed // todo: separate out app_model from list_model

        },

        change_to_tab: function(key){
            console.log('change to tab called: ' + key);
            this.$('li').removeClass('active');
            this.$('#' + key).addClass('active');

            if( _.has(this._tabbed_resources, key)){
                // NOTE: calling empty like this removes all of the event handlers from the views,
                // so we have to call delegateEvents when adding it back.  figure out a better way.
                $('#tab_container').empty();
                this._tabbed_resources[key].invoke();
            }else{
                window.alert('Unknown tab: ' + key);
            }
        },

        render : function() {
            window.scrollTo(0, 0);
            return this;
        }
    });

    return UserView;
});