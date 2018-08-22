define([
    'jquery',
    'underscore',
    'backbone',
    'layoutmanager',
    'iccbl_backgrid',
    'models/app_state',
    'templates/background_job_panel.html'
], function($, _, Backbone, layoutmanager,Iccbl, appModel, layout ) {
 
  var BackgroundJobPanelView = Backbone.Layout.extend({

      template: _.template(layout),
      
      JOB_CHECK_INTERVAL_MS: 10000, // 10 seconds
      
      // NOTE: it is good to synchronize the ellipsis interval with the 
      // JOB_CHECK_INTERVAL; see hmsiccbl.css .loading-ellipsis
      animated_ellipsis: [
        '&nbsp;<span class="loading-ellipsis"></span>'
      ].join(''),
      
      events: {
        'click button#close': 'close'
      },

      initialize : function() {
        self = this;
        console.log('initialize BackgroundJopPanel', this.model);
        _.bindAll(this, 'close','afterRender');
        
        this.jobResource = appModel.getResource('job');
        
        function retrieveModel(){
          console.log('retreive job state...')
          if (_.contains(['completed','failed'] , self.finalJobState)){
            clearInterval(self.jobTimer);
            return;
          }
          self.model.fetch({global: false })
            .done(function(model){
              console.log('new model', model);
              self.afterRender();
            })
            .fail(function(){
              var arguments = arguments;
              var failEl = $('<a>', {
                href: '#', title: 'show failure response', class: 'alert-link'
              }).text("Server Error");
              failEl.click(function(e){
                e.preventDefault();
                Iccbl.appModel.jqXHRfail.apply(this,arguments); 
              });
              self.$el.removeClass('alert-info');
              self.$el.addClass('alert alert-danger alert-dismissible');
              self.$el.append(failEl);
              Iccbl.appModel.jqXHRfail.apply(this,arguments); 
            });
        }
        
        this.jobTimer = setInterval(retrieveModel, self.JOB_CHECK_INTERVAL_MS); 
        
      },
      
      serialize: function() {
        return {}
      },
      
      update: function(){
        
      },
      
      afterRender: function(){
        
        var self = this;
        var $el = self.$el;
        $el.empty();
        
        function jobLink(){
          var jobId = self.model.get('id');
          var hrefTemplate = 
            self.jobResource['fields']['id']['display_options']['hrefTemplate'];
          var interpolatedVal = Iccbl.formatString(hrefTemplate,self.model);
          var title = 'Job: ' + jobId;
          try{ 
            var jobJson = self.model.pick('username','id','process_id',
              'uri','params','method','comment','state',
              'date_time_requested','date_time_submitted','date_time_processing',
              'date_time_completed','response_status_code');
            
            var sep = '\n';
            var sepRegex = new RegExp(sep, "g");
            
            function truncateLongElement(obj){
              // Truncate response if too long
              var testString = appModel.print_dict(obj, sep);
              var rowCount = (testString.match(sepRegex) || []).length;
              if (rowCount > appModel.MAX_ROWS_IN_DIALOG_MSG){
                testString = testString.split(sepRegex).slice(0,10);
                testString.push('-- truncated --');
                return testString.join(sep);
              }else{
                return obj;
              }
            }
            var response_content = JSON.parse(self.model.get('response_content'));
            if (!_.isEmpty(response_content)){
               var temp = truncateLongElement(response_content);
               if (temp !== response_content){
                 jobJson['response_content'] = temp;
               }else{
                 jobJson['response_content'] = response_content;
               }
            }
            var context_data = JSON.parse(self.model.get('context_data'));
            if (!_.isEmpty(context_data)){
               var temp = truncateLongElement(context_data);
               if (temp !== context_data){
                 jobJson['context_data'] = temp;
               }else{
                 jobJson['context_data'] = context_data;
               }
            }
            var finalJobJson = {};
            _.each(_.keys(jobJson), function(key){
              var title = self.jobResource['fields'][key].title;
              finalJobJson[title] = jobJson[key];
            });
            title = appModel.print_dict(finalJobJson, sep);
          }catch(e){
            console.log('error parsing job data', e, self.model.toJSON());
          }
          var link = $('<a>', {
            tabIndex : -1,
            href : interpolatedVal,
            class: 'alert-link',
            target : self.target,
            title: title
          }).text('Job: ' + jobId);
          return link;
        };
        
        function getState() {
          function showResponse(e) {
            e.preventDefault();
            var response_content = self.model.get('response_content');
            try{
              responseJSON = JSON.parse(response_content);
              responseJSON['Status Code'] = self.model.get('response_status_code');
              appModel.showJsonMessages(responseJSON, { title: 'Job Response'});
            }catch(e){
              console.log(
                'Error, unable to parse response as JSON: ', response_content);
              appModel.showModalMessage({
                buttons_on_top: false,
                body: response_content,
                title: 'Error'
              });
            }
          };
          var state = self.model.get('state');
          if (_.contains(['pending', 'submitted', 'processing'], state)){
            state += self.animated_ellipsis;
          }
          if (_.contains(['failed','completed'], state)){
            var stateEl = $('<a>', {
              href: '#', title: 'show response', class: 'alert-link'
            }).text(state);
            stateEl.click(showResponse);
            state = stateEl;
          }
          return state;
        };          
        
        var dismissLink = 
          $('<button type="button" class="close" aria-label="Close">&times;</button>');
        dismissLink.click(self.close);
        
        $el.append(jobLink());
        $el.append(
          Iccbl.formatString(': {method}: {uri}', self.model)
        );
        $el.append(', State: ');
        $el.append(getState());
        var state = self.model.get('state');
        if (state == 'pending'){
          $el.removeClass('alert-success alert-danger');
          $el.addClass('alert alert-info')
        }else if (state=='submitted'){
          $el.removeClass('alert-success alert-danger');
          $el.addClass('alert alert-info')
        }else if (state=='processing'){
          $el.removeClass('alert-success alert-danger');
          $el.addClass('alert alert-info')
        }else if (state=='completed'){
          $el.removeClass('alert-info');
          $el.addClass('alert alert-success alert-dismissible')
          $el.append(dismissLink);
          console.log('clear internal:', self.jobTimer);
          clearInterval(self.jobTimer);
        }else if (state=='failed'){
          $el.removeClass('alert-info');
          $el.addClass('alert alert-danger alert-dismissible');
          $el.append(dismissLink);
          console.log('clear internal:', self.jobTimer);
          clearInterval(self.jobTimer);
        }
      },
      
      close: function(e) {
        e.preventDefault();
        console.log('close panel...');
        this.model.collection.remove(this.model);
        this.remove();
      }

    });

    return BackgroundJobPanelView;
});