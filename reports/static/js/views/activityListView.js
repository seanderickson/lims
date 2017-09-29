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
  
  var ActivityView = ListView.extend({
    
    initialize: function(args) {
      this._args = args;
      var resource = args.resource;
      var linkCell = Iccbl.LinkCell.extend({

        render: function(){
          
          this.$el.empty();
          var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
          
          var hrefTemplate = this.hrefTemplate;
          
          var activityType = this.model.get('activity_class');
          if (_.contains(['libraryscreening','externallibraryscreening'],activityType)){
            hrefTemplate = '#screen/{screen_facility_id}' + 
              '/summary/libraryscreening/{activity_id}';
          }else if (_.contains(['cplt','cherrypickscreening'],activityType)){
            hrefTemplate = '#screen/{screen_facility_id}' + 
              '/cherrypickrequest/{cherry_pick_request_id}/cherrypickplates';
          }
          else {
            if (!_.isEmpty(this.model.get('serviced_user_id'))){
              hrefTemplate = '#screensaveruser/{serviced_user_id}/serviceactivity/{activity_id}';
            } 
            else if (!_.isEmpty(this.model.get('screen_facility_id'))){
              hrefTemplate = '#screen/{screen_facility_id}/activities/{activity_id}';
            } else {
              console.log('Activity does not have a serviced user or screen!', this.model)
            }
          }
          var interpolatedVal = Iccbl.formatString(hrefTemplate,this.model);
          this.$el.append($('<a>', {
            tabIndex : -1,
            href : interpolatedVal,
            target : this.target,
            title: this.title
          }).text(formattedValue));
          return this;
        }
      },resource.fields['activity_id'].display_options);

      resource.fields['activity_id']['backgridCellType'] = linkCell;
      resource.fields['date_of_activity']['backgridCellType'] = linkCell;

      ActivityView.__super__.initialize.apply(this, arguments);      
      
    }
    
  });

  return ActivityView;
});