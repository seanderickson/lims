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
      var self = this;
      this._args = args;
      this.screen = args.screen;
      this.user = args.user;
      var resource = args.resource;

      console.log('ActivityListView', args);
      var linkCell = Iccbl.LinkCell.extend({

        render: function(){
          
          this.$el.empty();
          var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
          var hrefTemplate = this.hrefTemplate;
          var classification = this.model.get('classification');
          if (classification == 'screening'){
            var activity_type = this.model.get('type');
            if (_.contains(['library_screening','ext_library_screening'], activity_type)){
              hrefTemplate = '#screen/{screen_facility_id}' + 
                '/summary/libraryscreening/{activity_id}';
            }else if (_.contains(['cp_transfer','cp_screening'], activity_type)){
              hrefTemplate = '#screen/{screen_facility_id}' + 
                '/cherrypickrequest/{cherry_pick_request_id}/cherrypickplates';
            }else{
              console.log('unknown activity_type for classifcation', 
                activity_type, classification, this.model);
            }
          }else{
            if (self.screen){
              hrefTemplate = '#screen/{screen_facility_id}/activities/{activity_id}';
            }
            else if (self.user){
              hrefTemplate = '#screensaveruser/{serviced_user_id}/activity/{activity_id}';
            } else {
              if (this.model.has('serviced_user_id')){
                hrefTemplate = '#screensaveruser/{serviced_user_id}/activity/{activity_id}';
              } 
              else if (this.model.has('screen_facility_id')){
                hrefTemplate = '#screen/{screen_facility_id}/activities/{activity_id}';
              } else {
                console.log('Activity does not have a serviced user or screen!', this.model)
              }
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
      resource.fields['date_of_activity']['backgridCellType'] = linkCell.extend({
        'formatter': _.extend({}, Iccbl.StringFormatter.prototype, {
          fromRaw(rawValue){
            return Iccbl.getDateString(rawValue);
          }
        })
      });

      ActivityView.__super__.initialize.apply(this, arguments);      
    }
  });

  return ActivityView;
});