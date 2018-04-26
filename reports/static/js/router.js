define([
  'jquery',
  'underscore',
  'backbone',
  'models/app_state'
], 
function($, _, Backbone, appModel) { 

  var AppRouter = Backbone.Router.extend({

    initialize : function() {
      var self = this;
      // send all routes to URIstack processing function
      this.route(/(.*)/, "toPath", this.toPath);
      this.routesHit = 0;
      Backbone.history.on(
        'route', 
        function(router, route, params)
        {
          self.routesHit++;
          console.log('detected route: ' + route + ', params: ' 
              + JSON.stringify(params) + ', routesHit:' + self.routesHit);
        }, this);

      this.listenTo(appModel, 'change:uriStack', this.uriStackChange);
      console.log('router initialized...');
    },
    
    /** 
     * Pull out complex keys in search - to allow for slashes in the keys
     * recursive to grab multiple search terms 
     **/
    toPath: function(path){
      function popKeysAndDecode(stack){
        if(!_.isEmpty(stack)){
          // 20180425 - fix: window location automatically encodes the hash url,
          // and backbone router generates a second route change with the 
          // encoded url; so clean it so that the final uriStack is the same.
          stack = _.map(stack, function(fragment){
            fragment = decodeURIComponent(fragment);
            fragment = ''+fragment;
            return fragment;
          });
          var searchIndex = _.indexOf(stack,appModel.URI_PATH_SEARCH);
          if(searchIndex > -1){
            var newStack = stack.slice(0,searchIndex+1);
            var keys = [];
            stack = stack.slice(searchIndex+1);
            while(!_.isEmpty(stack)){
              var temp = stack.shift();
              if (!_.contains(appModel.LIST_ARGS, temp)){
                keys.push(temp);
              }else{
                stack.unshift(temp);
                break;
              }
            }
            temp = [keys.join('/')];
            newStack = newStack.concat(temp);
            newStack = newStack.concat(popKeys(stack));
            return newStack;
          }
        }
        return stack;
      }
      var uriStack = [];
      if (path){
        uriStack = popKeysAndDecode(path.split('/'));
      }
      appModel.set({ uriStack: uriStack }, { source: this });
    },
    
    back: function() {  
      if(this.routesHit >= 1) {
        console.log('back, routesHit: ' + this.routesHit);
        // More than one route hit -> user did not land to current page directly
        this.routesHit--;
        window.history.back();
      } else {
        console.log('first route in site, back to home...');
        // Otherwise go to the home page. Use replaceState if available so
        // the navigation doesn't create an extra history entry
        this.navigate('/', {trigger:true, replace:true});
      }
    },

    /**
     * Generate a route that can be used by navigate, from the current state.
     */
    get_route: function(stack){
      return stack.join('/');
    },
    
    uriStackChange: function(model, vals, options){
      if(options.source === this){
        console.log('self generated uristack change');
        return;
      }else{
        this.model_set_route();
      }
    },
    
    model_set_route: function() {
      var uriStack = appModel.get('uriStack');
      console.log('model_set_route: ' + JSON.stringify(uriStack));
      var route = this.get_route(uriStack);
      // TODO: this mirrors the handler for route match in main.js
      document.title = 'Screensaver LIMS' + ': ' + route;

      // Trigger false to suppress further parsing, 
      // Replace false (default) to create browser history
      var options = { trigger: false, replace:false }; // , replace:false
      var routing_options = appModel.get('routing_options');
      appModel.set({ routing_options: {} });
      if(!_.isUndefined(routing_options)){
          options = _.extend(options, routing_options);
      }
      // Clear out error messages after navigating away from page
      appModel.clearErrors();
      console.log('route: ', route, options);
      
      this.navigate( route, options );
    }

  });

  return AppRouter;
});