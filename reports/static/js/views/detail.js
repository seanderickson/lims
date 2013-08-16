define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_modelbinder',
    'text!templates/generic-detail.html',
    'text!templates/generic-form.html',

], function( $, _, Backbone, modelbinder, genericDetailTemplate, genericFormTemplate ) {
    var DetailView = Backbone.View.extend({

        events: {
            'click button#save': 'save',
            'click button#edit': 'edit',
            'click button#cancel': 'cancel'
        },

        initialize : function(attributes, options) {
            // console.log('initialize: ' + JSON.stringify(attributes) + ', ' + JSON.stringify(options));
            // console.log('DetailView initializer: options: ' + JSON.stringify(options));

            keys = _(this.model.attributes).keys().sort(function(a,b){
                order_a = options.fields[a]['order_by'];  // TODO: need an edit order by
                order_b = options.fields[b]['order_by'];
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
                template = genericFormTemplate;
            }
            console.log(' template: ', template, 'fields: ', options.fields['screensaver_user_id']['title'] );
            var compiledTemplate = _.template( template,
                { 'fieldDefinitions': options.fields ,
                  title: options.title,
                  keys: _(keys),
                  object: this.model.attributes
                });

            this.$el.html(compiledTemplate);
            // this.modelBinder = new modelbinder();
            // this.modelBinder.bind(this.model, this.el);

            this.listenTo(this.model, 'error', this.error);
            _.bindAll(this, 'history');
            console.log('initialize DetailView');
        },

        edit: function(event) {
            console.log(' template: ', genericFormTemplate, 'fields: ', this._options.fields['screensaver_user_id']['title'], ', k: ' , this._keys, ', t: ', this._options.title );
            var compiledTemplate = _.template( genericFormTemplate,
                { 'fieldDefinitions': this._options.fields ,
                  title: this._options.title,
                  keys: _(this._keys)
                });
            this.$el.html(compiledTemplate);
            this.modelBinder = new modelbinder();
            this.modelBinder.bind(this.model, this.el);

        },

        save: function(event){
            event.preventDefault();
            var self = this;
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
            this.modelBinder.unbind();
        }
    });

    return DetailView;
});