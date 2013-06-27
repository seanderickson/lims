define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_modelbinder',
    'text!templates/generic-form.html',

], function( $, _, Backbone, modelbinder, genericFormTemplate ) {
    var DetailView = Backbone.View.extend({

        events: {
            'click button#save': 'save',
            'click button#cancel': 'cancel'
        },

        initialize : function(attributes, options) {
            // console.log('initialize: ' + JSON.stringify(attributes) + ', ' + JSON.stringify(options));
            // console.log('DetailView initializer: options: ' + JSON.stringify(options));
            var compiledTemplate = _.template( genericFormTemplate,
                { 'fieldDefinitions': options.fieldDefinitions ,
                  title: options.title,
                  object: _(_(this.model.attributes).keys())
                });
            this.$el.html(compiledTemplate);
            this.modelBinder = new modelbinder();
            this.modelBinder.bind(this.model, this.el);

            this.listenTo(this.model, 'error', this.error);
            console.log('initialize DetailView');
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