define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'iccbl_backgrid',
    'text!templates/generic-detail.html',
    'text!templates/generic-form-stickit.html',
    'text!templates/modal_ok_cancel.html',

], function( $, _, Backbone, stickit, Iccbl, genericDetailTemplate, genericFormTemplate, modalOkCancel ) {
    var DetailView = Backbone.View.extend({

        // bindings: function(){
            // console.log('evaluating bindings...');
            // return {
            // '#screensaver_user_id': 'screensaver_user_id'
            // };
        // },

        events: {
            'click button#save': 'save',
            'click button#delete': 'delete',
            'click button#edit': 'edit',
            'click button#cancel': 'cancel',
            'click button#history': 'history'
        },

        initialize : function(attributes, options) {
            var self = this;
            Iccbl.assert( !_.isUndefined(options.schemaResult), 'detail view requires a schemaResult struct');
            Iccbl.assert( !_.isUndefined(options.schemaResult['resource_definition']), 'detail view schemaResult requires a resource_definition');
            Iccbl.assert( !_.isUndefined(options.router), 'detail view requires a router');

            this._schemaResult = options['schemaResult'];
            this._resource_definition = this._schemaResult['resource_definition'];
            this._router = options.router;
            this._options = options;

            this._id = this.model.get('id');
            if(_.has(this._schemaResult['resource_definition'], 'title_attribute')){
                console.log('create id from ' + this._schemaResult['resource_definition']['title_attribute']);
                this._id = _.reduce(this._schemaResult['resource_definition']['title_attribute'],
                        function(memo, item){
                            if( self.model.has(item) ) memo += self.model.get(item)
                            else memo += item
                            return memo ;
                        }, '');
            }else{
                console.log('Warn: schema for this type has no resource_definition,id_attribute; type: ' + JSON.stringify(this._schemaResult));
            }
            console.log('id: ' + this._id);

            this._keys = _(this.model.attributes).keys().sort(function(a,b){
                order_a = self._schemaResult.fields[a]['ordinal'];  // TODO: need an edit order by
                order_b = self._schemaResult.fields[b]['ordinal'];
                if(_.isNumber(order_a) && _.isNumber(order_b)){
                    return order_a - order_b;
                }else if(_.isNumber(order_a)){
                    return -1;
                }else if(_.isNumber(order_b)){
                    return 1;
                }else{
                    return 0;
                }
            });

            if(options.isEditMode){
                // template = genericFormTemplate;
                this.edit(null);
            }else{
                this.detail(null);
            }
        },

        detail: function(event) {
            // filter keys for detail view
            this._options.isEditMode = false; // TODO: lets split edit to a different view?
            var self=this;
            var detailKeys = _(this._keys).filter(function(key){
                return _.has(self._schemaResult.fields[key], 'visibility') && _.contains(self._schemaResult.fields[key]['visibility'], 'detail');
            });
            console.log('detail keys: ' + JSON.stringify(detailKeys));

            var compiledTemplate = _.template( genericDetailTemplate,
                { 'fieldDefinitions': this._schemaResult.fields ,
                  title: this._resource_definition['title']+ ': ' + this._id,
                  keys: _(detailKeys),
                  object: this.model.attributes
                });

            this.$el.html(compiledTemplate);

            this.listenTo(this.model, 'error', this.error);
            console.log('initialize DetailView');
        },

        edit: function(event) {
            this._options.isEditMode = true;
            //console.log(' template: ', genericFormTemplate, 'fields: ', this._options.fields['screensaver_user_id']['title'], ', k: ' , this._keys, ', t: ', this._options.title );
            bindings = {};
            var self = this;
            var editKeys = _(this._keys).filter(function(key){
                return _.has(self._schemaResult.fields[key], 'visibility') && _.contains(self._schemaResult.fields[key]['visibility'], 'edit');
            });

            _.each(editKeys, function(key){
                if( _(self._schemaResult.fields).has(key)){
                    option = self._schemaResult.fields[key];
                    if(option.ui_type == 'choice' || option.ui_type == 'multiselect' ){
                        var _optionsCollection = [];
                        if(_.has(option, 'choices')){
                            _optionsCollection = option.choices.map(function(choice){
                                return { label: choice, value: choice };
                            });
                        }else{
                            window.alert('Warning, no choices defined for: ' + key); // TODO: use bootstrap alerts!
                            option.choices = _(_optionsCollection); // so the template doesn't complain
                        }

                        if(option.ui_type == 'choice' ){ // radio type choice
                            // Note: stickit uses the radio button element class, not the id
                            bindings['.radio_' + key] = {
                                observe: key,
                                selectOptions: { collection: _optionsCollection } };
                        }
                        else if(option.ui_type == 'choice2' ){  // traditional drop down type choice
                            bindings['#' + key] = {
                                observe: key,
                                selectOptions: { collection: _optionsCollection } };
                        }
                        else if(option.ui_type == 'multiselect' ){ // checkbox type multiselect
                            // Note: stickit uses the checkbox element class, not the id
                            bindings['.checkbox_' + key] = {
                                observe: key,
                                selectOptions: { collection: _optionsCollection } };
                        }
                        else if(option.ui_type == 'multiselect2' ){ // traditional multiselect box
                            bindings['#' + key] = {
                                observe: key,
                                selectOptions: { collection: _optionsCollection } };
                        }
                    }
                    else if(option.ui_type == 'boolean' ){ //
                        bindings['.checkbox_' + key] = {
                                observe: key,
                                selectOptions: { collection: [{label: 'select', value: "True" }] }
                            };
                    }else{
                        bindings['#' + key] = key;
                    }
                }else{
                    bindings['#' + key] = key;
                }
            });
            // console.log('bindings: ' + JSON.stringify(bindings));
            // console.log('model:' + JSON.stringify(this.model.attributes));
            var compiledTemplate = _.template( genericFormTemplate,
                { 'fieldDefinitions': this._schemaResult.fields ,
                  title: this._resource_definition['title']+ ': ' + this._id,
                  keys: _(editKeys)
                });
            this.$el.html(compiledTemplate);
            this.stickit(this.model,  bindings);

            // this.modelBinder = new modelbinder();
            // this.modelBinder.bind(this.model, this.el);
        },

        save: function(event){
            event.preventDefault();
            var self = this;
            console.log('save: ' + JSON.stringify(this.model));
            this.model.save(null, {
                success: function(model, resp){
                    console.log('success');
                    self.$el.empty();
                    //self.trigger('remove');
                    self.detail(null);

                },
                error: function(model,xhr){
                    var re = /([\s\S]*)Request Method/;
                    var match = re.exec(xhr.responseText);
                    if (match) window.alert('error saving: ' + match[1] + ': ' + xhr.status + ':' + xhr.statusText );
                    else window.alert('error on saving: ' + xhr.status + ':' + xhr.statusText );
                }

            });
        },

        delete: function(event){
            event.preventDefault();
            var self = this;
            console.log('delete: ' + JSON.stringify(this.model));
            var model = this.model;
            var modalDialog = new Backbone.View({
                el: _.template(modalOkCancel, { body: "Please confirm deletion of record: '" + model.get('toString') + "'", title: "Please confirm deletion" } ),
                events: {
                    'click #modal-cancel':function(event) {
                        console.log('cancel button click event, '); // + JSON.stringify(fieldDefinitions));
                        event.preventDefault();
                        $('#modal').modal('hide'); // TODO: read-up on modal!  this is not ideal with the reference to template elements!
                    },
                    'click #modal-ok':function(event) {
                        console.log('ok button click event, '); // + JSON.stringify(fieldDefinitions));
                        event.preventDefault();
                        model.destroy();
                        $('#modal').modal('hide');
                        self.$el.empty();
                        self.trigger('remove');
                        self._router.back();
                    }
                },
            });
            modalDialog.render();
            $('#modal').empty();
            $('#modal').html(modalDialog.$el);
            $('#modal').modal();
        },

        cancel: function(event){
            event.preventDefault();
            // TODO: do we have to do anything to abort form changes?
            // this.close(); // can we do this
            // this.$el.remove();
            if(this._options.isEditMode){
                this.detail(null);
            }else{
                this.$el.empty();
                this.trigger('remove');
                this._router.back();
            }
        },

        history: function(event){
            console.log('history click');
            event.preventDefault();
            this.$el.empty();
            this.trigger('remove');
            var self = this;

            // construct the route to the history
            // TODO: consider changing the app state to get to the records (future benefit of not reloading models if cached?)
            var _history_search = "ref_resource_name=" + this._resource_definition['key'];
            var id = this.model.get('id');
            if(_.has(this._resource_definition, 'id_attribute')){
                //console.log('create id from ' + this._resource_definition['id_attribute']);
                id = _.reduce(this._resource_definition['id_attribute'],
                        function(memo, item){
                            if(!_.isEmpty(memo)) memo += '/';
                            return memo += self.model.get(item);}, '');
            }else{
                console.log('Warn: schema for this type has no resource_definition,id_attribute; type: ' + this._options.type);
            }
            _history_search += ',key='+id;
            // TODO: discuss whether full encoding of the search fragment is necessary.
            // to-date, we know that the forward slash messes up the backbone router parsing, but other URL chars do not,
            // and full encoding reduces the usability of the URL for the end user
            //_history_search = encodeURIComponent(_history_search);
            _history_search = _history_search.replace('\/','%2F');
            //console.log('history_search:' + _history_search);

            var _route = 'list/apilog/search/' + _history_search; //ref_resource_name=' + _resource_name + ',key='+ id;
            console.log('-- set route: ' + _route);
            this._router.navigate(_route, {trigger: true});
        },

        error: function(options){
            console.log('A model save error reported: ' + JSON.stringify(options));
        },

        render : function() {
            console.log('render detail_stickit');
            return this;
        },

        onClose: function(){
            console.log('...detail view close method');
            //this.$el.remove();
            //this.modelBinder.unbind();
        }
    });

    return DetailView;
});