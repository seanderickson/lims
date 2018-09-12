define([
    'jquery',
    'underscore',
    'backbone',
    'iccbl_backgrid',
    'layoutmanager',
    'models/app_state',
    'templates/contact.html'
], function( $, _, Backbone, Iccbl, layoutmanager, 
      appModel, contactTemplate) {

	var ContactView = Backbone.Layout.extend({

	  initialize: function(args) {
	    console.log('initialize contact view');
	    
	    var self = this;
	    
	    var app_data = appModel.getAppData();
	    console.log('app_data', app_data.toJSON());
	    this.contactInfo = [
        ['Screening Facility', 'Screening Facility Information', 
          Iccbl.formatString('<a href="{facility_url}">{facility_name}</a>', 
            app_data)],
        ['Site URL', 'Screensaver URL', 
          Iccbl.formatString('<a href="{site_url}">{site_url}</a>', app_data)],
        ['Feedback', 'Submit feedback and questions', 
          Iccbl.formatString(
            '<a href="mailto:{contact_feedback_email}">{contact_feedback_email}</a>', 
            app_data)],
	    ];
	    console.log('contactInfo', this.contactInfo);
	  },
	  
    afterRender : function() {
      var self = this;
      return this;
    },

    serialize: function() {
      console.log('serialize:', this.contactInfo);
      return {
        'data': this.contactInfo, 
      };      
    },    
    
    template: _.template(contactTemplate)
    
	});

	return ContactView;
});