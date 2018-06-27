/**
 * Application loading script for the Iccbl-lims app.
 */

require('css/external/bootstrap_custom_build/css/bootstrap.css');
//require('bootstrap/dist/css/bootstrap.css');

require('backgrid/lib/backgrid.css');
require('backgrid-paginator/backgrid-paginator.css');
require('backgrid-filter/backgrid-filter.css');
require('multiselect/css/multi-select.css');
require('bootstrap-datepicker/dist/css/bootstrap-datepicker3.css');
require('jquery-chosen/chosen.css');
require('bootstrap-chosen/bootstrap-chosen.css');
require('jquery-bonsai/jquery.bonsai.css');
require('css/jquery.bonsai.overrides.css');
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
    'models/menu_fixture_screener.json', 
    'bootstrap'
  ],
function($, _, Backbone, Iccbl, appModel, AppView, AppRouter, 
  reports_ui_resources_raw, reports_menu_raw, ui_resources_raw, menu_raw,
  screener_menu_raw ) {
  
  console.log('init screensaver/reports...')
  
  
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

  var screener_menu_resource = JSON.parse(screener_menu_raw);
  appModel.set('screener_menu', screener_menu_resource);

  var loadCount = 0

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

  // sending a csrftoken with every ajax request
  function csrfSafeMethod(method) {
      // these HTTP methods do not require CSRF protection
      return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
  };
  $.ajaxSetup({
      crossDomain: false, // obviates need for sameOrigin test
      beforeSend: function(xhr, settings) {
          if (!csrfSafeMethod(settings.type)) {
              xhr.setRequestHeader("X-CSRFToken", appModel.readCookie('csrftoken'));
          }
      },
      statusCode: {
      401: function(err){
        console.log('Login Failed.', err.responseJSON);
        authErrorHandler(err);
      }
    }      
  });

  $(document).bind("ajaxSend", function(event, jqxhr, settings){
    console.log('ajaxSend', arguments);
// Remove: "global: false" flag should  prevent trigger event handler 
// see: app_state.initialized()
//    if ( settings.url.indexOf(appModel.reportsApiUri + '/job') >= 0 ){
//      // suppress the spinner for job panel updates
//      // FIXME: show spinner if explicitly loading the jobs view itself (use router)
//      return;
//    }
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

  function authErrorHandler(){
    window.logged_in = false;
    window.location='/accounts/login/?next=' + 
      window.location.pathname + window.location.hash;
    
  };
  
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
