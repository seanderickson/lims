define([
    'jquery',
    'underscore',
    'backbone',
    'layoutmanager',
    'text!templates/messages.html'
], function($, _, Backbone, layoutmanager, layout ) {
    
  var MessageView = Backbone.Layout.extend({

      template: _.template(layout),
      
      events: {
        'click button#close': 'close'
      },

      initialize : function() {
          console.log('initialize MessageView');
          _.bindAll(this, 'close');
      },
      
      serialize: function() {
        return {
          messages: _.chain(this.model.get('messages'))
        }
      },
      
      close: function() {
        this.model.unset('messages', {silent: true});
        this.remove();
      }

    });

    return MessageView;
});