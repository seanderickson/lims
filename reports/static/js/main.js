/**
 * Application loading script for the Iccbl-lims app.
 */

// workaround for libraries that need jquery $ in webpack
//window.$ = window.jQuery = require('jquery');  

require('css/external/bootstrap_custom_build/css/bootstrap.css');
//require('bootstrap/dist/css/bootstrap.css');

require('backgrid/lib/backgrid.css');
require('backgrid-paginator/backgrid-paginator.css');
require('backgrid-filter/backgrid-filter.css');
require('multiselect/css/multi-select.css');
require('bootstrap-datepicker/dist/css/bootstrap-datepicker3.css');
require('jquery-chosen/chosen.css');
require('bootstrap-chosen/bootstrap-chosen.css');
require('css/hmsiccbl.css');

require([ // now load application code
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'models/app_state',
    'views/app_view',
    'router',
    'models/reports_ui_resources_fixture.json', 
    'models/reports_menu_fixture.json', 
    'models/ui_resources_fixture.json', 
    'models/menu_fixture.json', 
    'bootstrap'
  ],
function($, _, Backbone, Iccbl, appModel, AppView, AppRouter, 
  reports_ui_resources_raw, reports_menu_raw, ui_resources_raw, menu_raw ) {
  
  console.log('init screensaver/reports...')
  
  /** 
   * NOTE: The ajax "traditional" setting is needed to serialize the ordering
   * array sent to the server and used to implement multisort with the 
   * tastypie server.
   * With traditional serialization, the array values are serialized as
   * order=val1&order=val1&order=...
   * see http://api.jquery.com/jQuery.param/c
   **/
  $.ajaxSettings.traditional = true;

  /**
   * For ajax POST operations: send the session id as a request
   * header (CSRF protection in SessionAuthentication requires).
   * (the other option is to use Basic Authentication)
   * see:
   * tastypie/authentication.is_authenticated
   * and
   * http://stackoverflow.com/questions/15388694/does-sessionauthentication-work-in-tastypie-for-http-post
   */
  function readCookie(name) {
    var nameEQ = escape(name) + "=";
    var ca = document.cookie.split(';');
    for (var i = 0; i < ca.length; i++) {
        var c = ca[i];
        while (c.charAt(0) === ' ') c = c.substring(1, c.length);
        if (c.indexOf(nameEQ) === 0) return unescape(c.substring(nameEQ.length, c.length));
    }
    return null;
  };
  // sending a csrftoken with every ajax request
  function csrfSafeMethod(method) {
      // these HTTP methods do not require CSRF protection
      return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
  };
  $.ajaxSetup({
      crossDomain: false, // obviates need for sameOrigin test
      beforeSend: function(xhr, settings) {
          if (!csrfSafeMethod(settings.type)) {
              xhr.setRequestHeader("X-CSRFToken", readCookie('csrftoken'));
          }
      }
  });
  
  // Prevent keypresses in input elements from propogating to the form submit button
  $(document).on("keypress", "input:not(textarea):not([type=submit])", function(event) {
    if (event.keyCode == 13) {
        event.preventDefault();
    }
  });
  /**
   * Augment the view prototype to prevent memory leaks -
   * See: http://lostechies.com/derickbailey/2011/09/15/zombies-run-managing-page-transitions-in-backbone-apps/
   * Todo: is this still needed with the Backbone.Layout extension managing all the views?
   **/
  Backbone.View.prototype.close = function(){
    this.remove();
    this.unbind();
    if (this.onClose){
      this.onClose();
    }
  };

  // Check logged in status
  if(_.isUndefined(window.logged_in) || window.logged_in != 'True' ){
    console.log('window.logged_in: ' + window.logged_in );
    window.location='/accounts/login/?next=' + 
      window.location.pathname + window.location.hash;
    return;
  }
  
  // Bootstrap the resources hash
  var ui_resources = JSON.parse(reports_ui_resources_raw);
  _.extend(ui_resources, JSON.parse(ui_resources_raw));
  appModel.set('ui_resources', ui_resources);

  // Bootstrap the menu hash
  var menu_resource = JSON.parse(menu_raw);
  _.extend(menu_resource['submenus'], JSON.parse(reports_menu_raw)['submenus']);
  appModel.set('menu', menu_resource);

  var loadCount = 0
  $(document).bind("ajaxSend", function(){
    $('#loading').fadeIn({duration:100});
    loadCount++;
    console.log('add to loadCount: ' + loadCount );
  }).bind("ajaxComplete", function(){
    loadCount--;
    if (loadCount <= 0 ){
      loadCount = 0;
      console.log('loadCount: ' + loadCount );
      $('#loading').fadeOut({duration:100});
    }else{
      $('#loading').show();
      console.log('loadCount: ' + loadCount );
    }
  }).bind("ajaxStop"), function(){
    loadCount--;
    if (loadCount <= 0 ){
      loadCount = 0;
      console.log('loadCount: ' + loadCount );
      $('#loading').fadeOut({duration:100});
    }else{
      $('#loading').show();
      console.log('loadCount: ' + loadCount );
    }
  };
  
  var appRouter = appModel.router = new AppRouter({ model: appModel });
  var appView = new AppView({ model: appModel },{ router: appRouter});

  appModel.start(function(){
    console.log('Render application')
    appView.$el.appendTo("#application_div")
    appView.render();
        
    Backbone.history = Backbone.history || new Backbone.History({});
    Backbone.history.start({ pushState: false, root: appModel.get('root_url') });
  });
  
  // Set the document title
  Backbone.history.on('route', function(router, route, params) {
    var title = _.reduce(
        params,
        function(memo, item){
          if(item){
            if (memo !== ' ') memo += ', ';
            memo += item;
          }
          return memo ;
        }, ' ');              
    
      document.title = 'Screensaver LIMS' + ':' + title;
   }, this);    
 
});
