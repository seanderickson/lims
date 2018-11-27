/**
 * Application loading script for the Iccbl-lims app.
 */

require('external/bootstrap_custom_build/css/bootstrap.css');
//require('bootstrap/dist/css/bootstrap.css');

require('backgrid/lib/backgrid.css');
require('backgrid-paginator/backgrid-paginator.css');
require('backgrid-filter/backgrid-filter.css');
require('multiselect/css/multi-select.css');
require('bootstrap-datepicker/dist/css/bootstrap-datepicker3.css');
require('select2/dist/css/select2.css');
require('bootstrap-chosen/bootstrap-chosen.css');
require('jquery-bonsai/jquery.bonsai.css');
require('jquery.bonsai.overrides.css');
require('hmsiccbl.css');

require([ // now load application code
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'models/app_state',
    'views/app_view',
    'router',
    'models/reports_ui_resources_fixture.json', 
    'models/ui_resources_fixture.json', 
    'models/reports_menu_fixture.json', 
    'models/menu_fixture.json', 
    'models/menu_fixture_screener.json', 
    'bootstrap'
  ],
function($, _, Backbone, Iccbl, appModel, AppView, AppRouter, 
  reports_ui_resources, app_ui_resources, reports_menu, app_menu,
  screener_menu_resource ) {

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
  if(_.isUndefined(window.logged_in) || window.logged_in != true ){
    console.log('window.logged_in: ' + window.logged_in );
    window.location='/accounts/login/?next=' + 
      window.location.pathname + window.location.hash;
    return;
  }
  
  _.extend(reports_ui_resources, app_ui_resources);
  appModel.set('ui_resources', reports_ui_resources);

  _.extend(app_menu['submenus'], reports_menu['submenus']);
  appModel.set('menu', app_menu);

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

  var serverAccessTimeKey = 'serverAccessTimeKey';
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
    
    localStorage.setItem(serverAccessTimeKey, Date.now() );
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
    Backbone.history.start({ pushState: false, root: '/' });
  });
  
//  // Set the document title
//  Backbone.history.on('route', function(router, route, params) {
//    console.log('params', params);
//    var title = _.reduce(
//        params,
//        function(memo, item){
//          if(item){
//            if (memo !== ' ') memo += ', ';
//            memo += item;
//          }
//          return memo ;
//        }, ' ');              
//      title = title.trim();
//      if (!_.isEmpty(title)) title = ': ' + title;
//      title = appModel.getAppData().get('page_title') + title;
//      document.title = title;
//   }, this);   
  
  // Logout timer
  
  var defaultIdleTimeoutSec = 2 * 60 * 60; // 2 hours
  var defaultServerTimeoutSec = 8 * 60 * 60; // 8 hours
  var warnTimeout = 2 * 60 * 1000; // 2 minutes in ms
  var logoutUrl = window.logout_url;
  var idleTimer;
  var warnTimer;
  var serverTimeoutSeconds, idleTimeout;
  var timerKey = 'lastActivityTime';
  
  // Check the APP_PUBLIC_DATA from the server for settings and start the timer
  appModel.getAppDataFromServer(function(appData){
    serverTimeoutSeconds = defaultServerTimeoutSec;
    if (appData.has('SESSION_COOKIE_AGE')){
      serverTimeoutSeconds = parseInt(appData.get('SESSION_COOKIE_AGE'));
    }
    idleTimeout = defaultIdleTimeoutSec*1000;
    if (appData.has('BROWSER_IDLE_TIMEOUT_SECONDS')){
      idleTimeout = parseInt(appData.get('BROWSER_IDLE_TIMEOUT_SECONDS'))*1000;
    }
    
    console.log('startIdleTimer, serverTimeoutSeconds: ', serverTimeoutSeconds, 
      'idleTimeout (s)', idleTimeout/1000);
    startIdleTimer();
  });
  
  function logout() {
    window.location = logoutUrl + '?next=' + window.location.pathname + window.location.hash;
  }  

  function resetTimeOutTimer() {
    clearTimeout(warnTimer);
    clearTimeout(idleTimer);
    startIdleTimer();
    
    // Set the real time of last activity to local storage
    localStorage.setItem(timerKey, Date.now());
    // If time elapsed > server timeout, ping the server
    var serverAccessTime = localStorage.getItem(serverAccessTimeKey);
    if (!_.isUndefined(serverAccessTime) ){
      serverAccessTime = parseInt(serverAccessTime);
      var elapsed = Date.now() - serverAccessTime;
      elapsed = elapsed / 1000;
      if (elapsed > serverTimeoutSeconds){
        console.log('access server to renew session', elapsed, serverTimeoutSeconds);
        appModel.setCurrentUser();
      }
    }
  }
  
  function warnIdle() {
    console.log('warnIdle...');
    clearTimeout(idleTimer);
    
    // Check localStorage for the real time of last activity;
    // calculate the real time to idle timeout and restart the timer.
    
    var globalLastTime = localStorage.getItem(timerKey);
    var realTimeNow = Date.now();
    if (globalLastTime){
      console.log('globalLastTime', globalLastTime, ', realTimeNow', realTimeNow);
      globalLastTime = parseInt(globalLastTime);
      var realInterval = realTimeNow-globalLastTime;
      if ( realInterval < idleTimeout ){
        var realTimeout = idleTimeout - realInterval;
        idleTimer = setTimeout(warnIdle, realTimeout);
        console.log('activity detected in another window: realInterval(ms)', 
          realInterval, 'realTimeout from now(ms)', realTimeout);
        return;
      }
    }

    // If idle time is past, set warn timer and show modal
    
    warnTimer = setTimeout(logout, warnTimeout);
    
    appModel.showModal({
      title: 'Idle timeout warning',
      body: 'Select continue to renew your session',
      okText: 'Continue',
      ok: function(){
        resetTimeOutTimer();
        appModel.setCurrentUser();
      },
      cancelText: 'Logout now',
      cancel: function(){
        window.location = logoutUrl;
      }
    });
  }
  function startIdleTimer() {
    idleTimer = setTimeout(warnIdle, idleTimeout);
    localStorage.setItem(timerKey, Date.now());
  }

  document.onload = resetTimeOutTimer;
  document.onmousemove = resetTimeOutTimer;
  document.onmousedown = resetTimeOutTimer; // touchscreen presses
  document.ontouchstart = resetTimeOutTimer;
  document.onclick = resetTimeOutTimer;     // touchpad clicks
  document.onscroll = resetTimeOutTimer;    // scrolling with arrow keys
  document.onkeypress = resetTimeOutTimer;  

  // End: logout timer
  
});
