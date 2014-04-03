define([
  'jquery',
  'underscore',
  'backbone',
  'models/app_state'
], 
function($, _, Backbone, appModel) { 

  var AppRouter = Backbone.Router.extend({

    LIST_ROUTE_ORDER: ['rpp', 'page','order','search'],
    
    initialize : function() {
      // send all routes to URIstack processing function
      this.route(/(.*)/, "toPath", this.toPath);
      this.routesHit = 0;
      Backbone.history.on('route', function(router, route, params) {
              this.routesHit++;
              console.log('detected route: ' + route + ', params: ' 
                  + JSON.stringify(params) + ', routesHit:' + this.routesHit);
           }, this);

      this.listenTo(appModel, 'change:uriStack', this.uriStackChange);
      console.log('router initialized...');
    },
    
    toPath: function(path){
      console.log('toPath: ' + path);
      var uriStack = [];
      if (path) uriStack = path.split('/');
      appModel.set({ uriStack: uriStack }, { source: this });
    },

    back: function() {  
      if(this.routesHit >= 1) {
        console.log('back, routesHit: ' + this.routesHit);
        //more than one route hit -> user did not land to current page directly
        this.routesHit--;
        window.history.back();
      } else {
        console.log('first route in site, back to home...');
        //otherwise go to the home page. Use replaceState if available so
        //the navigation doesn't create an extra history entry
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

      // trigger false to suppress further parsing, 
      // replace false (default) to create browser history
      var options = { trigger: false }; // , replace:false
      var routing_options = appModel.get('routing_options');
      appModel.set({ routing_options: {} });
      if(!_.isUndefined(routing_options)){
          options = _.extend(options, routing_options);
      }
      this.navigate( route, options );
    }

  });

  return AppRouter;
});