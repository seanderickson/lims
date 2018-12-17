define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView) {
  
  var LibraryListView = ListView.extend({
    
    initialize: function(args) {
      var self = this;
      this._class = 'LibraryListView';
      this._args = args;
      var resource = args.resource;
      console.log('LibraryListView', args);

      var show_archived_control = $([
        '<label class="checkbox-inline pull-left" ',
        '   title="Show archived libraries" >',
        '  <input type="checkbox">Show archived libraries</input>&nbsp;',
        '</label>'
        ].join(''));
      show_archived_control.css('margin-left','10px');
      
      if (appModel.hasPermission(resource.key, 'read')){
        var extraControls = args.extraControls || [];
        args.extraControls = extraControls;
      
        extraControls.push(show_archived_control);
      }      
      LibraryListView.__super__.initialize.apply(this, arguments);      
    
      if (appModel.hasPermission(resource.key, 'read')){
        var initialSearchHash = self.listModel.get(appModel.URI_PATH_SEARCH);
        var initial_show_archived = 
          _.result(initialSearchHash, appModel.API_PARAM_SHOW_ARCHIVED, false);
        if (initial_show_archived && initial_show_archived.toLowerCase()=='true'){
          show_archived_control.find('input[type="checkbox"]').prop('checked',true);
        }
        show_archived_control.find('input[type="checkbox"]').change(function(e) {
          var searchHash = _.clone(self.listModel.get(appModel.URI_PATH_SEARCH));
          if (e.target.checked) {
            searchHash[appModel.API_PARAM_SHOW_ARCHIVED] = 'true';
          } else {
            searchHash[appModel.API_PARAM_SHOW_ARCHIVED] = 'false';
          }
          self.listModel.set(appModel.URI_PATH_SEARCH,searchHash);
        });
        
      }    
    },
    
    afterRender: function(){
      ListView.prototype.afterRender.apply(this,arguments);

    }
    
    
  });

  return LibraryListView;
});