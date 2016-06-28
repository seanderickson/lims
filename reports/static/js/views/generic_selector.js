define([
    'jquery',
    'underscore',
    'backbone',
    'templates/generic-selector.html'
], function($, _, Backbone, genericSelectorTemplate) {

  var GenericSelector = Backbone.View.extend({
    
    /** default to span4, adjust as needed **/
    className : 'pull-down pull-right',
    
    initialize : function(attributes, options) {

          this.listenTo(this.model, 'change', this.changeNotified);
          this._options = options;
          _.bindAll(this, 'changeNotified');
      },

      events : {
          "change #generic_selector" : "selectorChanged"
      },

      selectorChanged : function(e) {
          e.preventDefault();
          var option = e.currentTarget.value;
          this.model.set({
              'selection' : option
          });
      },

      changeNotified : function() {
          var selection = this.model.get('selection');
          this.$('#generic_selector').val(String(selection));
      },

      render : function() {
          this.$el.empty();
          
          // NOTE: for Underscore 1.7.0, templates will have to be initialize in 
          // two steps, like this.  So 1. compiled, 2. interpolated
          //
          var temp = _.template(genericSelectorTemplate);
          temp = temp(              {
                label : this._options.label,
                options : _(this._options.options)
              });
          this.$el.html(temp);
          if (!_.isUndefined(this._options.selectorClass)) {
              this.$('#generic_selector').removeClass().addClass(this._options.selectorClass);
          };
          this.changeNotified();
          this.delegateEvents();
          return this;
      }
  });

  return GenericSelector;
});