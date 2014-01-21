define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_pageable',
    'iccbl_backgrid',
    'views/detail_stickit',
    'views/list',
    'views/simple-list',
    'text!templates/generic-tabbed.html',
], function($, _, Backbone, BackbonePageableCollection, Iccbl, DetailView, ListView, SimpleListView, tabbedTemplate) {

    // for compatibility with require.js, attach PageableCollection in the right place on the Backbone object
    // see https://github.com/wyuenho/backbone-pageable/issues/62
    Backbone.PageableCollection = BackbonePageableCollection;

    var UserGroupView = Backbone.View.extend({
        events: {
            'click li': 'click_tab',
        },

        initialize : function(attributes, options) {
            this.options = options;
            var self = this;

            Iccbl.requireOptions(options, ['url_root', 'current_options', 'model', 'schema', 'router' ]);

            console.log('initialize UserGroupView: current options: ' + JSON.stringify(options.current_options));

            _.bindAll(this, 'setDetail', 'setUsers', 'setPermissions');

            this._tabbed_resources = {
                detail: { description: 'UserGroup Details', title: 'UserGroup Details', invoke : this.setDetail },
                users: { description: 'Users', title: 'Users', invoke : this.setUsers },
                permissions: { description: 'UserGroup Permissions', title: 'UserGroup Permissions', invoke : this.setPermissions },
            };

            var compiledTemplate = _.template( tabbedTemplate, { tab_resources: this._tabbed_resources, });

            this.$el.html(compiledTemplate);

            var tab = options.current_options['tab'];
            if(_.isUndefined(tab) || _.isEmpty(tab) ) tab = 'detail';

            this.listenTo(this.model, 'change:current_detail', this.change_current_detail);

            self.objectModel = self.options.model;
            self.objectSchema = self.options.schema;

            // set the id specifically on the model: backbone requires this to determine whether a "POST" or "PATCH" will be used
            self.objectModel.id = Iccbl.getIdFromIdAttribute(self.objectModel, self.objectSchema);
            self.title = Iccbl.getTitleFromTitleAttribute(self.objectModel, self.objectSchema);

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
                new DetailView({ model: self.objectModel },
                    {
                        schemaResult:self.objectSchema,
                        router:self.options.router,
                        isEditMode: false
                    });
            self.detail = detailView;
            self.$('#tab_container').html(self.detail.render().el);
            self.$('#detail').addClass('active'); // first time not clicked so set manually
        },

        setUsers: function(){
            var self = this;
            var resource_url = this.options.url_root + '/usergroup/' + self.objectModel.id + '/users';
            var schema_url =  this.options.url_root + '/user/schema';
            var url = resource_url;

            if(_.isUndefined(self.users)){
                var createUsers = function(schemaResult){
                    var columns = Iccbl.createBackgridColModel(schemaResult.fields, Iccbl.MyHeaderCell);//, col_options );

                    // make a checkbox column that will be used to manage the group->user relationship
                    columns.unshift({
                        column: "is_for_group",
                        name: "is_for_group",
                        label: "Select users for " + self.title,
                        cell: "boolean"
                    });

                     var ListModel = Backbone.Model.extend({
                        defaults: {
                            rpp: 25,
                            page: 1,
                            order: {},
                            search: {}}
                        });
                    var listModel = new ListModel();

                    var collection  = new Iccbl.MyCollection({
                            'url': url,
                            currentPage: parseInt(listModel.get('page')),
                            pageSize: parseInt(listModel.get('rpp')),
                            listModel: listModel
                        });

                    // collection.bind("change reset add remove", function(){
                    var currentUsers = self.objectModel.get('users');
                    collection.bind("sync", function(){
                        collection.each(function(model){  // each model is a user
                            var currentUserUri = model.get('resource_uri');
                            // model.set('is_for_group', Iccbl.containsByMatch(currentUsers, currentUserUri) );
                            model.set('is_for_group', _.contains(currentUsers, currentUserUri) );

                            model.on('change:is_for_group', function(model){
                                var checked = model.get('is_for_group');
                                var userUri = model.get('resource_uri');
                                if(checked){
                                  currentUsers.unshift(userUri);
                                }else{
                                  // NOTE: this will fail if userUri is a _localized_ URI and currentUsers is not, or vice versa.
                                  // See get_local_resource_uri in api.py
                                  currentUsers = _.without(currentUsers, userUri);
                                }
                                self.objectModel.save({'users':currentUsers}, {
                                  patch: true,
                                  success: function(){
                                    collection.fetch();
                                    // TODO: should not need this
                                    //self.objectModel.fetch();
                                  }
                                });
                            });
                        });
                    });

                    self.users = new SimpleListView({ model: listModel },{ columns: columns,collection: collection, schemaResult: schemaResult });
                    $('#tab_container').html(self.users.render().el);
                };

                Iccbl.getSchema(schema_url, createUsers);
            }else{
                self.users.setElement(self.$('#tab_container')).render();
            }
            self.$('#users').addClass('active'); // first time not clicked so set manually
        },

        setPermissions: function(){
            var self = this;
            // get all of the permissions first, then checkoff the ones that are for this user
            var resource_url = this.options.url_root + '/usergroup/' + self.objectModel.id + '/permissions';
            var schema_url =  this.options.url_root + '/permission/schema';
            var url = resource_url;

            if(_.isUndefined(this.permissions)){
                var createPermissions = function(schemaResult){
                    var columns = Iccbl.createBackgridColModel(schemaResult.fields, Iccbl.MyHeaderCell);//, col_options );

                    // make a checkbox column that will be used to manage the user->permission relationship
                    columns.unshift({
                        column: "is_for_group",
                        name: "is_for_group",
                        label: "Select permissions for " + self.title,
                        cell: "boolean"
                    });

                    var ListModel = Backbone.Model.extend({
                        defaults: {
                            rpp: 25,
                            page: 1,
                            order: {},
                            search: {}}
                        });
                    var listModel = new ListModel();

                    var collection  = new Iccbl.MyCollection({
                            'url': url,
                            currentPage: parseInt(listModel.get('page')),
                            pageSize: parseInt(listModel.get('rpp')),
                            listModel: listModel
                        });

                    collection.bind("sync", function(){
                        collection.each(function(model){ // each model is a permission
                            var currentPerms = self.objectModel.get('permissions');
                            var currentPermissionUri = model.get('resource_uri');
                            model.set('is_for_group', _.contains(currentPerms, currentPermissionUri));
                            // model.set('is_for_group', Iccbl.containsByMatch(currentPerms, currentPermissionUri));

                            // when the checkbox is selected, use a custom model event to update the group model relationship and refetch the user model
                            model.on('change:is_for_group', function(model){
                                var currentPerms = self.objectModel.get('permissions');
                                var permissionUri = model.get('resource_uri');
                                // var containsByMatch = Iccbl.containsByMatch(currentPerms, permissionUri);
                                var containsByMatch = _.contains(currentPerms, permissionUri);
                                var changed = false;
                                if(model.get('is_for_group') && ! containsByMatch ){
                                    currentPerms.push(permissionUri);
                                    changed = true;
                                }else if(!model.get('is_for_group')  && containsByMatch ) {
                                    currentPerms = _.without(currentPerms, permissionUri);
                                    changed = true
                                }
                                if(changed){
                                    self.objectModel.save({'permissions': currentPerms },{
                                        patch: true,
                                        success: function(){
                                            console.log('refetch the permission model');
                                            // model.fetch();
                                            collection.fetch(); // TODO: _this_ is why actions should be batch processed
                                        }
                                    });
                                }
                            });
                        });
                    });

                    collection.fetch();

                    self.permissions = new SimpleListView({ model: listModel },{ columns: columns,collection: collection, schemaResult: schemaResult });
                    $('#tab_container').html(self.permissions.render().el);
                };

                Iccbl.getSchema(schema_url, createPermissions);
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

    return UserGroupView;
});