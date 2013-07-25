define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'text!templates/generic-detail.html',
    'text!templates/generic-form-stickit.html',

], function( $, _, Backbone, stickit, genericDetailTemplate, genericFormTemplate ) {
    var DetailView = Backbone.View.extend({


        bindings: function(){
            console.log('evaluating bindings...');
            return {
            '#screensaver_user_id': 'screensaver_user_id'
            };
        },

        events: {
            'click button#save': 'save',
            'click button#edit': 'edit',
            'click button#cancel': 'cancel'
        },

        initialize : function(attributes, options) {
            // console.log('initialize: ' + JSON.stringify(attributes) + ', ' + JSON.stringify(options));
            // console.log('DetailView initializer: options: ' + JSON.stringify(options));

            keys = _(this.model.attributes).keys().sort(function(a,b){
                console.log('sorting: a: ' + a + ', ' + JSON.stringify(options.fields[a]) + ', b: ' + b + ', ' + JSON.stringify(options.fields[b]));
                order_a = options.fields[a]['ordinal'];  // TODO: need an edit order by
                order_b = options.fields[b]['ordinal'];
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
            this._keys = keys;  // TODO: need to put the sorted keys on the options and remove from this class
            template = genericDetailTemplate;
            this._options = options;
            if(options.isEditMode){
                // template = genericFormTemplate;
                this.edit(null);
            }else{
                //console.log(' template: ', template, 'fields: ', options.fields['screensaver_user_id']['title'] );
                // filter keys for detail view
                var detailKeys = _(keys).filter(function(key){
                    if(_.isUndefined(options.fields[key]['visibility']) || !options.fields[key]['visibility'].contains('detail')){
                        console.log('skip non-detail field: ' + key);
                        return false;
                    }else{
                        return true;
                    }
                });
                console.log('detail keys: ' + JSON.stringify(detailKeys));

                var compiledTemplate = _.template( template,
                    { 'fieldDefinitions': options.fields ,
                      title: options.title,
                      keys: _(detailKeys),
                      object: this.model.attributes
                    });

                this.$el.html(compiledTemplate);

                this.listenTo(this.model, 'error', this.error);
                console.log('initialize DetailView');
            }
        },

        edit: function(event) {
            //console.log(' template: ', genericFormTemplate, 'fields: ', this._options.fields['screensaver_user_id']['title'], ', k: ' , this._keys, ', t: ', this._options.title );
            bindings = {};
            var self = this;
            var editKeys = _(this._keys).filter(function(key){
                var options = self._options;
                if(_.isUndefined(options.fields[key]['visibility']) || !options.fields[key]['visibility'].contains('edit')){
                    console.log('skip non-edit field: ' + key);
                    return false;
                }else{
                    return true;
                }
            });

            _.each(editKeys, function(key){
                if( _(self._options.fields).has(key)){
                    option = self._options.fields[key];
                    // if(!option.visibility.contains('edit')){
                        // console.log('skip non-edit field: ' + key);
                        // return;
                    // }
                    if(option.ui_type == 'choice' || option.ui_type == 'multiselect' ){
                        console.log('--choice key: ' + key + ', ' + JSON.stringify(option));
                        var _optionsCollection = [];
                        if(_.has(option, 'choices')){
                            _optionsCollection = option.choices.map(function(choice){
                                return { label: choice, value: choice };
                            });
                        }else{
                            window.alert('Warning, no choices defined for: ' + key); // TODO: use bootstrap alerts!
                            option.choices = _(_optionsCollection); // so the template doesn't complain
                        }
                        console.log('-- optionsCollection for choice: ' + key + " , options: " + JSON.stringify(_optionsCollection));

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
                { 'fieldDefinitions': this._options.fields ,
                  title: this._options.title,
                  keys: _(editKeys)
                });
            this.$el.html(compiledTemplate);
            //this.stickit();
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
                    self.trigger('remove');
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
            // TODO: do we have to do anything to abort form changes?
            // this.close(); // can we do this
            // this.$el.remove();
            this.$el.empty();
            this.trigger('remove');
        },

        error: function(options){
            console.log('A model save error reported: ' + JSON.stringify(options));
        },

        render : function() {
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