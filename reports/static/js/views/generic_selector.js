define([
    'jquery',
    'underscore',
    'backbone',
    'text!templates/generic-selector.html'
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
//          if (!_.contains(this._options.options, '')) {
//              this._options.options.unshift('');
//              // create a blank entry
//          }
          this.$el.html(_.template(genericSelectorTemplate, {
              label : this._options.label,
              'options' : _(this._options.options)
          }));
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