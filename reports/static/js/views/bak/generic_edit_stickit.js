define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'iccbl_backgrid',
    'text!templates/generic-form-stickit.html',
    'text!templates/modal_ok_cancel.html',

], function( $, _, Backbone, stickit, Iccbl, genericFormTemplate, modalOkCancel ) {
    var EditView = Backbone.View.extend({


        events: {
            'click button#save': 'save',
            'click button#cancel': 'cancel',
            'click button#history': 'history'
        },

        initialize : function(attributes, options) {
            var self = this;
            Iccbl.requireOptions(options, ['resource', 'router']);
            Iccbl.assert( !_.isUndefined(options.resource['resource_definition']), 'detail view resource requires a resource_definition');

            this.resource = options['resource'];
            this._resource_definition = this.resource['resource_definition'];
            this._router = options.router;
            this._options = options;

            this._id = Iccbl.getTitleFromTitleAttribute(self.model, this.resource);
            console.log('id: ' + this._id);

            this._keys = Iccbl.sortOnOrdinal(_.keys(this.model.attributes), self.resource.fields);
            this._keys = _(this._keys).filter(function(key){
                return _.has(self.resource.fields[key], 'visibility') && _.contains(self.resource.fields[key]['visibility'], 'edit');
            });

            bindings = {};
            _.each(this._keys, function(key){
                if( _(self.resource.fields).has(key)){
                    option = self.resource.fields[key];
                    console.log('option: ' + JSON.stringify(option));
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
                            // Note: stickit uses the radio button element _class_, not the id
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
                        bindings['#' + key] = {
                            observe: key,
                            events: ['blur'] };
                    }
                }else{
                    bindings['#' + key] = {
                            observe: key,
                            events: ['blur'] };
                }
            });
            _.bindAll(this, 'render');
        },

        save: function(event){
            event.preventDefault();
            var self = this;
            console.log('save: ' + JSON.stringify(this.model));
            this.model.save(null, {
                success: function(model, resp){
                    console.log('success');
                    self.trigger('remove');
                    self._router.back();
                },
                error: function(model,xhr){
                    var re = /([\s\S]*)Request Method/;
                    var match = re.exec(xhr.responseText);
                    if (match) window.alert('error saving: ' + match[1] + ': ' + xhr.status + ':' + xhr.statusText );
                    else window.alert('error on saving: ' + xhr.status + ':' + xhr.statusText );
                }
            });
        },

        cancel: function(event){
            event.preventDefault();

            this.trigger('remove');
            this._router.back();
        },

        error: function(options){
            console.log('A model save error reported: ' + JSON.stringify(options));
        },

        render : function() {
            var self = this;
            console.log('render generic_edit_stickit');

            var compiledTemplate = _.template( genericFormTemplate,
                { 'fieldDefinitions': this.resource.fields ,
                  title: this._resource_definition['title']+ ': ' + this._id,
                  keys: _(this._keys)
                });
            this.$el.html(compiledTemplate);
            this.stickit(this.model,  bindings);


            return this;
        },

        onClose: function(){
            console.log('...edit view close method');
            this.$el.empty();
            //this.$el.remove();
            //this.modelBinder.unbind();
        }
    });

    return EditView;
});