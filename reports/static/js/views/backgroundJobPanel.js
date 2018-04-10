define([
    'jquery',
    'underscore',
    'backbone',
    'layoutmanager',
    'iccbl_backgrid',
    'models/app_state',
    'templates/background_job_panel.html'
], function($, _, Backbone, layoutmanager,Iccbl, appModel, layout ) {
 
  // alert classes alert-warning alert-dismissible
  //  <button type="button" class="close" id="close" 
  //    data-dismiss="alert" aria-hidden="true">&times;</button>

  
  var BackgroundJobPanelView = Backbone.Layout.extend({

      template: _.template(layout),
      
      JOB_CHECK_INTERVAL_MS: 10000, // 10 seconds
      
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
          self.model.fetch()
            .success(function(model){
              console.log('new model', model);
              self.afterRender();
            });
        }
        
        this.jobTimer = setInterval(retrieveModel, self.JOB_CHECK_INTERVAL_MS); 
        
      },
      
      serialize: function() {
        return {
//          messages: _.chain(this.model.get('messages'))
        }
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
            var jobJson = self.model.toJSON();
            jobJson['response_content'] = JSON.parse(jobJson['response_content']);
            jobJson['context_data'] = JSON.parse(jobJson['context_data']);
            title = appModel.print_dict(jobJson, '\n');
          }catch(e){
            console.log('error parsing job data', self.model.toJSON());
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
              appModel.showJsonMessages(responseJSON);
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
          if (_.contains(['failed','completed'], state)){
            state = $('<a>', {
              href: '#', title: 'show response', class: 'alert-link'
            }).text(state);
            state.click(showResponse);
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
          clearInterval(self.jobTimer);
        }else if (state=='failed'){
          $el.removeClass('alert-info');
          $el.addClass('alert alert-danger alert-dismissible');
          $el.append(dismissLink);
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