define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_forms',
    'iccbl_backgrid',
    'models/app_state',
    'views/generic_edit',
    'views/screen/plateRangeSearch',
    'templates/search-box.html'
], function($, _, Backbone, backbone_forms, Iccbl, appModel, EditView, 
    PlateRangeSearchView, searchBoxTemplate) {
    
  
  var SearchView = Backbone.Layout.extend({

    template: _.template(searchBoxTemplate),

    formTemplate: _.template([
        "<form class='iccbl-headerfield-form' ",
        ">",
        '<div data-fieldsets/>',
        '<div id="data-error" class="has-error" ></div>',
        "</form>",
        ].join('')),

    fieldTemplate: _.template([
        '<div class="form-group" key="form-group-<%=key %>" >',
        '<label title="<%= help %>" for="<%= editorId %>"><%= title %></label>',
        '<div>',
        '  <span key="<%=key%>" placeholder="testing.." data-editor></span>',
    '      <div data-error class="text-danger" ></div>',
        '</div>',
        '</div>'
        ].join('')),
    
    choiceFieldTemplate: _.template([
      '<div class="form-group" key="form-group-<%=key%>" >',
      '    <div class="" >',
      '      <div data-editor  key="<%=key%>" style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
      '      <div data-error class="text-danger" ></div>',
      '      <div><%= help %></div>',
      '    </div>',
      '  </div>',
    ].join('')),
   
    initialize: function(args) {
      
      var self = this;
      self.objects_to_destroy = [];
      console.log('---- initialize search_box');

      var args = args || {};

//      var TextArea2 = Backbone.Form.editors.TextArea.extend({
//        render: function() {
//          TextArea2.__super__.render.apply(this,arguments);
//          this.$el.attr('placeholder', this.schema.placeholder);
//          return this;
//        },        
//      });
      
      var schema2 = {};
      function validatePlateSearch(value, formValues){
        var errors = [];
        var parsedData = Iccbl.parseRawPlateSearch(value,errors);
        if (_.isEmpty(parsedData)){
          errors.push('no values found for input');
        } else {
          console.log('parsedData', parsedData);
        }
        if (!_.isEmpty(errors)){
          return {
            type: 'copyplate',
            message: errors.join('; ')
          };
        }
      };
      schema2['copyplate'] = {
        title: 'Search for Plate & Copy',
        key: 'copyplate',
        help: 'enter a comma separated list',
        placeholder: 'e.g. 1000-1005 B\n2000 C',
        validators: ['required',validatePlateSearch],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var FormFields2 = Backbone.Model.extend({
        schema: schema2,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var formFields2 = new FormFields2({
      });
      
      var form2 = self.form2 = new Backbone.Form({
        model: formFields2,
        template: self.formTemplate,
        el: '#search-box-2'
      });
      
      ///// Well search
      var schema3 = {};
      schema3['well'] = {
        title: 'Search for Wells',
        key: 'well',
        help: 'enter a comma separated list',
        placeholder: 'e.g. 1000 A08,A09,A10\n1740:D12',
        validators: ['required'],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var FormFields3 = Backbone.Model.extend({
        schema: schema3,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var formFields3 = new FormFields3({
      });
      
      var form3 = self.form3 = new Backbone.Form({
        model: formFields3,
        template: self.formTemplate,
        el: '#search-box-3'
      });
      
      ///// Screening Inquiry
      var schema5 = {};
      function validateScreeningInquiry(value, formValues){
        var errors = [];
        var parsedData = Iccbl.parseRawScreeningInquiry(value,errors);
        if (_.isEmpty(parsedData)){
          errors.push('no values found for input');
        } else {
          console.log('parsedData', parsedData);
        }
        if (!_.isEmpty(errors)){
          return {
            type: 'screening_inquiry',
            message: errors.join('; ')
          };
        }
      };
      schema5['screening_inquiry'] = {
        title: 'Screening Inquiry',
        key: 'screening_inquiry',
        help: 'Enter a screening inquiry',
        placeholder: 'e.g. (Screen #) 1000-2000,2015,3017 100 nL x 2',
        validators: ['required',validateScreeningInquiry],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var FormFields5 = Backbone.Model.extend({
        schema: schema5,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var formFields5 = new FormFields5({
      });
      
      var form5 = self.form5 = new Backbone.Form({
        model: formFields5,
        template: self.formTemplate,
        el: '#search-box-5'
      });      
      
      ///// Cherry Pick Request
      var schema4 = {};
      schema4['cpr'] = {
        title: 'CPR #',
        key: 'cpr',
        help: 'enter a Cherry Pick ID',
        placeholder: 'e.g. 44443',
        key:  'search_value', // TODO: "key" not needed>?
        type: 'Number',
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var FormFields4 = Backbone.Model.extend({
        schema: schema4,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var formFields4 = new FormFields4({
        cpr: null
      });
      
      var form4 = self.form4 = new Backbone.Form({
        model: formFields4,
        template: self.formTemplate,
        el: '#search-box-4'
      });
      
      this.listenTo(appModel, 'change:uriStack' , this.uriStackChange );
      
      SearchView.__super__.initialize.apply(this, args);
      _.bindAll(this,'buildDynamicForms');
      console.log('initialize SearchView');
    },

    buildDynamicForms: function() {
      var self = this;
      
      var schema1 = {};
      schema1['user_id'] = {
        title: '',
        key: 'user_id',
        type: EditView.ChosenSelect,
        editorClass: 'chosen-select',
        placeholder: 'Find a User...',
        options: appModel.getUserOptions(),
        template: self.choiceFieldTemplate 
      };
      schema1['screen_facility_id'] = {
        title: '',
        key: 'screen_facility_id',
        type: EditView.ChosenSelect,
        editorClass: 'chosen-select',
        placeholder: 'Find a Screen...',
        options: appModel.getScreenOptions(),
        template: self.choiceFieldTemplate
      };
      schema1['library_short_name'] = {
        title: '',
        key: 'library_short_name',
        type: EditView.ChosenSelect,
        editorClass: 'chosen-select',
        placeholder: 'Find a Library...',
        options: appModel.getLibraryOptions(),
        template: self.choiceFieldTemplate
      };
      var FormFields1 = Backbone.Model.extend({
        schema: schema1,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var formFields1 = new FormFields1({
      });

      if (self.form1){
          self.form1.remove();
      }
      var form1 = self.form1 = new Backbone.Form({
        model: formFields1,
        template: self.formTemplate
      });
      self.form1.on("change", function(){

        appModel.requestPageChange({
          ok: function(){
            var user_id = self.form1.getValue('user_id');
            var screen_facility_id = self.form1.getValue('screen_facility_id');
            var library_short_name = self.form1.getValue('library_short_name');
            if(user_id){
              var route = 'screensaveruser/' + user_id;
              self.form1.setValue('user_id',null);
              self.form1.$el.find('[key="user_id"]')
                  .find('.chosen-select').trigger("chosen:updated");
              appModel.router.navigate(route, {trigger: true});
            }
            else if(screen_facility_id){
              var route = 'screen/' + screen_facility_id;
              self.form1.setValue('screen_facility_id',null);
              self.form1.$el.find('[key="screen_facility_id"]')
                  .find('.chosen-select').trigger("chosen:updated");
              appModel.router.navigate(route, {trigger: true});
            }
            else if(library_short_name){
              var route = 'library/' + library_short_name;
              self.form1.setValue('library_short_name',null);
              self.form1.$el.find('[key="library_short_name"]')
                  .find('.chosen-select').trigger("chosen:updated");
              appModel.router.navigate(route, {trigger: true});
            }
          }
        });
      });
      $('#search-box-1').html(this.form1.render().$el);
      
      console.log('setup single selects using chosen...');
      // See http://harvesthq.github.io/chosen/
      this.$el.find('.chosen-select').chosen({
        width: '100%',
        allow_single_deselect: true,
        search_contains: true
        });
      self.stopListening(appModel,'change:users');
      self.listenTo(appModel,'change:users', function(){
        console.log('user options changed');
        appModel.getUserOptions(function(options){
          self.buildDynamicForms();
          self.render();
        });
      });
      
      self.stopListening(appModel,'change:screens');
      self.listenTo(appModel,'change:screens', function(){
        console.log('screen options changed');
        appModel.getScreenOptions(function(options){
          self.buildDynamicForms();
          self.render();
        });
      });
      
      self.stopListening(appModel,'change:libraries');
      self.listenTo(appModel,'change:libraries', function(){
        console.log('library options changed');
        appModel.getLibraryOptions(function(options){
          self.buildDynamicForms();
          self.render();
        });
      });
    },
    
    afterRender: function() {
      var self=this;
      // Pre-fetch options for the search_box
      $(this).queue([
        appModel.getScreenOptions,
        appModel.getUserOptions,
        appModel.getLibraryOptions,
        self.buildDynamicForms
       ]);
      
      // TODO: 20170809 - lazy load the forms
      // 20170809 - see Select2 for ajax loading (instead of Chosen)
      
      //var lchosen = this.form1.$el.find('[key="library_short_name"]').find('.chosen-container');
      //lchosen.click(function(e){
      //  var lchosenselect = self.form1.$el.find('[key="library_short_name"]').find('.chosen-select');
      //  if (lchosenselect.find('option').length==0){
      //    e.preventDefault();
      //    appModel.getLibraryOptions(function(options){
      //      _.each(options, function(new_vocab_item){
      //        lchosenselect.append($('<option>',{
      //          value: new_vocab_item['key']
      //        }).text(new_vocab_item['title']));
      //      });
      //      self.form1.$el.find('[key="library_short_name"]').find('.chosen-select').trigger("chosen:updated");
      //      
      //    });
      //    
      //  }
      //});
      
      ///// CopyPlate search 2 - perform search on server
      var form2 = this.form2;
      var $form2 = this.form2.render().$el;
      if (this.form2_data){
        this.form2.setValue('copyplate', this.form2_data);
      }
      $('#search-box-2').html($form2);
      $form2.append([
          '<button type="submit" ',
          'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
          ].join(''));

      $form2.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        var errors = form2.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          console.log('form2 errors, abort submit: ' + JSON.stringify(errors));
          $form2.find('[name="copyplate"]').addClass(self.errorClass);
          return;
        }else{
          $form2.find('[name="copyplate"]').removeClass(self.errorClass);
        }
        var text_to_search = self.form2.getValue('copyplate');
        errors = [];
        var parsedSearchArray = Iccbl.parseRawPlateSearch(text_to_search,errors);
        if (!_.isEmpty(errors)){
          throw Exception('Unexpected errors after submit:', errors);
        }
        var resource = appModel.getResource('librarycopyplate');
        if (parsedSearchArray.length <= 5){
          // encode simple searches as a URL param
          var encodedPlateSearches = [];
          _.each(parsedSearchArray, function(parsedPlateSearch){
            encodedPlateSearches.push(_.map(
              parsedPlateSearch.combined, encodeURIComponent).join(','));
          });
          var uriStack = [resource.key, 'search', 
            'raw_search_data='+ encodedPlateSearches.join(encodeURIComponent(';'))];
          console.log('route: ', uriStack);
          appModel.router.navigate(uriStack.join('/'), {trigger:true});
        }else{
          // Send complex search data as a POST
          // must change the route, and create a post
          // TODO: key the search data using the searchID: 
          // this allows for browser "back" in the session
          // will also need to listen to URIStack changes and grab the data
          // from the search ID
          var searchID = ( new Date() ).getTime();
          appModel.setSearch(searchID,text_to_search);
          this.searchID = searchID;
          appModel.set('routing_options', {replace: false});  
          var _route = resource.key + '/search/'+ searchID;
          appModel.router.navigate(_route, {trigger:true});
        }
      });      
      
      ///// Well/Reagent search 2 - perform search on server
      var form3 = this.form3;
      var $form3 = this.form3.render().$el;
      if (this.form3_data){
        this.form3.setValue('well', this.form3_data);
      }
      $('#search-box-3').html($form3);
      $form3.append([
          '<button type="submit" ',
          'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
          ].join(''));

      $form3.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        var errors = form3.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          console.log('form3 errors, abort submit: ' + JSON.stringify(errors));
          $form3.find('[name="well"]').addClass(self.errorClass);
          // FIXME: add errors to the form
          return;
        }else{
          $form3.find('[name="well"]').removeClass(self.errorClass);
        }
        var text_to_search = self.form3.getValue()['well'];
        // must change the route, and create a post
        
        var searchID = ( new Date() ).getTime();
        appModel.setSearch(searchID,text_to_search);
        this.searchID = searchID;
        appModel.set('routing_options', {replace: false});  
        // Use the router to process the search: see content.js for handler
        var resource = appModel.getResource('reagent');
        var _route = resource.key + '/search/'+ searchID;
        appModel.router.navigate(_route, {trigger:true});
        // var newStack = [resource.key,'search',searchID];
        // NOTE: easier to control the router history using navigate: 
        // when using uristack, there is the problem of who set appModel.routing_options last:
        // a race condition is set up between list2.js and search_box.js
        //        appModel.set({'uriStack': newStack});     
      });      
      
      ///// Screening Inquiry
      if (appModel.hasPermission('libraryscreening', 'read')){
        var form5 = this.form5;
        var $form5 = this.form5.render().$el;
        if (this.form5_data){
          this.form5.setValue('screening_inquiry', this.form5_data);
        }
        $('#search-box-5').html($form5);
        $form5.append([
            '<button type="submit" ',
            'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
            ].join(''));
  
        $form5.find('[ type="submit" ]').click(function(e){
          e.preventDefault();
          var errors = form5.commit({ validate: true }); 
          if(!_.isEmpty(errors)){
            console.log('form5 errors, abort submit: ' + JSON.stringify(errors));
            $form5.find('[name="screening_inquiry"]').addClass(self.errorClass);
            return;
          }else{
            $form5.find('[name="screening_inquiry"]').removeClass(self.errorClass);
          }
  
          var searchValue = form5.getValue('screening_inquiry');
          var errorArray = [];
          var parsedSearch = Iccbl.parseRawScreeningInquiry(searchValue,errorArray);
          if (!_.isEmpty(errorArray)){
            throw Exception('Unexpected errors after submit:', errorArray);
          }
          var urlSearchParts = 
            PlateRangeSearchView.prototype.encodeFormData.call(this,parsedSearch);
          var uriStack = ['screen', parsedSearch.screen_facility_id,
                          'summary','plateranges','search',
                          urlSearchParts.join(appModel.SEARCH_DELIMITER)]
          console.log('route: ', uriStack);
          appModel.router.navigate(uriStack.join('/'), {trigger:true});
        });      
      }      

      ///// CPR
      if (appModel.hasPermission('cherrypickrequest', 'read')){
        var form4 = this.form4;
        var $form4 = this.form4.render().$el;
        $('#search-box-4').html($form4);
        $form4.append([
            '<button type="submit" ',
            'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
            ].join(''));
  
        $form4.find('[ type="submit" ]').click(function(e){
          e.preventDefault();
          var errors = form4.commit({ validate: true }); 
          if(!_.isEmpty(errors)){
            console.log('form4 errors, abort submit: ' + JSON.stringify(errors));
            $form4.find('#cpr').addClass(self.errorClass);
            return;
          }else{
            $form4.find('#cpr').removeClass(self.errorClass);
          }
          var cpr_id = self.form4.getValue()['cpr'];
          var resource = appModel.getResource('cherrypickrequest');
          var _route = ['#', resource.key,cpr_id].join('/');
          appModel.set('routing_options', {replace: false});  
          appModel.router.navigate(_route, {trigger:true});
        });      
      }

      this.uriStackChange();
    },

    /**
     * Backbone.Model change event handler
     * If the uriStackChange event is for a search handled with the search box,
     * attempt to set the search box value appropriately.
     * 
     * @param options.source = the event source triggering view
     */
    uriStackChange: function(model, val, options) {
      var self=this;
      self.cleanup();
      var uriStack = _.clone(appModel.get('uriStack'));
      var uiResourceId = uriStack.shift();
      
      if (_.contains(uriStack,'search')){
        var index = _.indexOf(uriStack,'search');
        var searchData = uriStack[index+1];
        
        // Test the search data as a simple search ("raw_search_data" param is
        // set as a URI param).
        if (searchData.indexOf('raw_search_data') > -1){
          var search_data = decodeURIComponent(searchData.split('=')[1]);
          console.log('setup search box for search data', search_data);
          if (uiResourceId=='librarycopyplate'){
            var errors = [];
            var parsedData = Iccbl.parseRawPlateSearch(search_data,errors);
            if (!_.isEmpty(errors)){
              console.log('Search data not parsed properly', errors);
              return;
            }
            parsedData = _.map(parsedData, function(parsedLine){
              return parsedLine.combined.join(' ');
            }).join('\n');
            this.form2.setValue('copyplate', parsedData);
            // form2 may not be rendered yet, store the value
            this.form2_data = parsedData;
          }else if (uiResourceId == 'reagent') {
            // TODO: "raw_search_data" not implemented yet for well/reagent 20170516
            this.form3.setValue('well', search_data);
            // form3 may not be rendered yet, store the value
            this.form3_data = search_data;
          }
          return;
        }
        
        if (self.form5 && uiResourceId == 'screen') {
          var screenId = uriStack.shift();
          if (_.contains(uriStack, 'plateranges')){
            if (_.contains(uriStack,'search')){
              var index = _.indexOf(uriStack,'search');
              var searchData = uriStack[index+1];
              var search_data = decodeURIComponent(searchData);
              
              var errors = [];
              var plateRangeSearchData = 
                Iccbl.parseScreeningInquiryURLParam(search_data, errors);
              if (!_.isEmpty(errors)){
                console.log('Search data not parsed properly', errors);
                return;
              }
              var si_formatter = new Iccbl.SIUnitsFormatter({ symbol: 'L' });
              parsedData = '(' + screenId + ') ';
              parsedData += plateRangeSearchData.plate_search;
              parsedData += ' ' + si_formatter.fromRaw(plateRangeSearchData.volume_required);
              parsedData += ' x' + plateRangeSearchData.replicate_count;
              self.form5.setValue('screening_inquiry',parsedData);
              // form3 may not be rendered yet, store the value
              self.form5_data = parsedData;
            }
          }
          return;
        }
        
        // Otherwise, test the searchData as a searchID for a complex search.
        
        if(searchData == this.searchID){
          console.log('self generated uristack change');
          return;
        }
        if(searchData && _.isNaN(parseInt(searchData))){
          console.log('search ID is not a valid number: ' + searchData);
          return;
        }
        var searchID = searchData;
        var search_data = appModel.getSearch(searchID);
        if(_.isUndefined(search_data) || _.isEmpty(search_data)){
          var msg = 'Search box requires a "search_data:'+searchID +'" in current browser state';
          console.log(msg);
          return;
        }else{
          
          console.log('setup search box for search data', search_data);
          if (uiResourceId=='librarycopyplate'){
            this.form2.setValue('copyplate', search_data);
            // form2 may not be rendered yet, store the value
            this.form2_data = search_data;
            this.searchID = searchID;
          }else if (uiResourceId == 'reagent') {
            this.form3.setValue('well', search_data);
            // form3 may not be rendered yet, store the value
            this.form3_data = search_data;
            this.searchID = searchID;
          }
        }
      }else{
        console.log('No search found: Search box requires the "search/[searchID]" param');
        return;
      }
    },
    
    // Backbone layoutmanager callback
    cleanup: function(){
      this.form2.setValue('copyplate', null);
      this.form3.setValue('well', null);
      this.form4.setValue('cpr', null);
      this.form5.setValue('screening_inquiry', null)
    },
    
  });

  return SearchView;
});