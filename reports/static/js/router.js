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

      this.listenTo(appModel, 'change:reportedUriStack', this.reportUriStack);
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
            newStack = newStack.concat(popKeysAndDecode(stack));
            return newStack;
          }
        }
        return stack;
      }
      var uriStack = [];
      if (path){
        uriStack = popKeysAndDecode(path.split('/'));
      }
      appModel.set({ uriStack: uriStack });
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
    
    reportUriStack: function(model, vals, options){
      var uriStack = appModel.get('reportedUriStack');
      console.log('reportUriStack: ' + JSON.stringify(uriStack));
      var route = this.get_route(uriStack);

      // TODO: this mirrors the handler for route match in main.js
      document.title = 'Screensaver LIMS' + ': ' + route;

      console.log('route: ', route, options);
      
      // Trigger false to suppress further parsing, 
      // Replace false (default) to create browser history
      var options = _.extend({ trigger: false, replace:false }, options); 
      this.navigate( route, options );
    }

  });

  return AppRouter;
});