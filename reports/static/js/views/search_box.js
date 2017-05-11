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
      console.log('---- initialize search_box');

      var args = args || {};

      var TextArea2 = Backbone.Form.editors.TextArea.extend({
        render: function() {
          TextArea2.__super__.render.apply(this,arguments);
          this.$el.attr('placeholder', this.schema.placeholder);
          return this;
        },        
      });
      
      var schema2 = {};
      schema2['copyplate'] = {
        title: 'Search for Plate & Copy',
        key: 'copyplate',
        help: 'enter a comma separated list',
        placeholder: 'e.g. 1000-1005 B\n2000 C',
        key:  'search_value', // TODO: "key" not needed>?
        type: TextArea2,
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
        key:  'search_value', // TODO: "key" not needed>?
        type: TextArea2,
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
        type: TextArea2,
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
      
//      // Pre-fetch options for the search_box
//      $(this).queue([
//         appModel.getScreenOptions,
//         appModel.getUserOptions,
//         appModel.getLibraryOptions,
//         self.buildDynamicForms]);
      
      this.listenTo(appModel, 'change:uriStack' , this.uriStackChange );
      
      SearchView.__super__.initialize.apply(this, args);
      _.bindAll(this,'buildDynamicForms');
      console.log('initialize SearchView');
    },

    buildDynamicForms: function() {
      var self = this;
      
      var schema1 = {};
      schema1['username'] = {
        title: '',
        key: 'username',
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
            var username = self.form1.getValue('username');
            var screen_facility_id = self.form1.getValue('screen_facility_id');
            var library_short_name = self.form1.getValue('library_short_name');
            if(username){
              var route = 'screensaveruser/' + username;
              self.form1.setValue('username',null);
              self.form1.$el.find('[key="username"]')
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

      appModel.on('change:userOptions', function(){
        console.log('user options changed');
        appModel.getUserOptions(function(options){
          self.buildDynamicForms();
          self.render();
        });
      });
      
      appModel.on('change:screenOptions', function(){
        console.log('screen options changed');
        appModel.getScreenOptions(function(options){
          self.buildDynamicForms();
          self.render();
        });
      });
      
      appModel.on('change:libraryOptions', function(){
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
         self.buildDynamicForms]);
      
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
        var text_to_search = self.form2.getValue()['copyplate'];
        Iccbl.parseRawPlateSearch(text_to_search,errors);
        if(!_.isEmpty(errors)){
          console.log('form2 errors, abort submit: ' + JSON.stringify(errors));
          $form2.find('[name="copyplate"]').addClass(self.errorClass);
          // FIXME: add errors to the form
          return;
        }else{
          $form2.find('[name="copyplate"]').removeClass(self.errorClass);
        }
        // must change the route, and create a post
        
        var resource = appModel.getResource('librarycopyplate');
        // TODO: key the search data using the searchID: 
        // this allows for bfowser "back" in the session
        // will also need to listen to URIStack changes and grab the data
        // from the search ID
        var searchID = ( new Date() ).getTime();
        appModel.setSearch(searchID,text_to_search);
        this.searchID = searchID;
        appModel.set('routing_options', {replace: false});  
        var _route = resource.key + '/search/'+ searchID;
        appModel.router.navigate(_route, {trigger:true});
        // var newStack = [resource.key,'search',searchID];
        // NOTE: easier to control the router history using navigate: 
        // when using uristack, there is the problem of who set appModel.routing_options last:
        // a race condition is set up between list2.js and search_box.js
        //        appModel.set({'uriStack': newStack});     
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
      var form5 = this.form5;
      var $form5 = this.form5.render().$el;
      if (this.form5_data){
        this.form5.setValue('well', this.form5_data);
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
        var urlSearchParts = PlateRangeSearchView.prototype.encodeFormData.call(this,parsedSearch);
        var uriStack = ['screen', parsedSearch.screen_facility_id,
                        'summary','plateranges','search',
                        urlSearchParts.join(appModel.SEARCH_DELIMITER)]
        console.log('route: ', uriStack);
        appModel.router.navigate(uriStack.join('/'), {trigger:true});
      });      
      
      ///// CPR
      var form4 = this.form4;
      var $form4 = this.form4.render().$el;
      $('#search-box-4').html($form4);
      $form4.append([
          '<button type="submit" class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
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

    },

    /**
     * Backbone.Model change event handler
     * @param options.source = the event source triggering view
     */
    uriStackChange: function(model, val, options) {
      var self=this;

      var uriStack = _.clone(appModel.get('uriStack'));
      var uiResourceId = uriStack.shift();

      if(uriStack.length > 1 && uriStack[0] == 'search'){
        var stack_change_type = uriStack.shift();
        console.log('create a search_box with search data');
        var searchID = uriStack.shift();
        
        if(searchID == this.searchID){
          console.log('self generated uristack change');
          return;
        }
        if(searchID && _.isNaN(parseInt(uriStack[1]))){
          console.log('search ID is not a valid number: ' + searchID);
          return;
        }
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
        console.log('Search box requires the "search/[searchID]" param');
        return;
      }
    } 
    
//    uriStackChange: function(model, val, options) {
//      var self=this;
//
//      var uriStack = _.clone(appModel.get('uriStack'));
//      console.log('uristack',uriStack);
//      var uiResourceId = uriStack.shift();
//      if(uriStack.length > 1 && uriStack[0] == 'search'){
//        var stack_change_type = uriStack.shift();
//        console.log('create a search_box with search data');
//        var searchID = uriStack.shift();
//        
//        if(searchID == this.searchID){
//          console.log('self generated uristack change');
//          return;
//        }
//        if(searchID && _.isNaN(parseInt(uriStack[1]))){
//          console.log('search ID is not a valid number: ' + searchID);
//          return;
//        }
//        var search_data = appModel.getSearch(searchID);
//        if(_.isUndefined(search_data) || _.isEmpty(search_data)){
//          var msg = 'Search box requires a "search_data:'+searchID +'" in current browser state';
//          console.log(msg);
//          return;
//        }else{
//          
//          console.log('setup search box for search data', search_data);
////          var search_string = '';
////          
////          search_string = _.reduce(search_data, function(memo, orclause){
////            if(memo != '') memo +='\n';
////            var fieldKeys = [
////                'library_short_name__icontains', 'plate_number__in',
////                'plate_number__range', 'copy_name__in', 'well_name__in', 'well_id__in'];
////            var temp = '';
////            _.each(fieldKeys, function(key){
////              if(_.has(orclause, key)){
////                var val = orclause[key];
////                if(key == 'plate_number__range'){
////                  val = val.split(',').join('-');
////                }
////                if(temp !==''){
////                  temp += ' '
////                }
////                temp +=  val;
////              }
////            });
////            
////            return memo + temp;
////          }, '');
////          console.log('setValue: ' + search_string)
//          
//          if (uiResourceId=='librarycopyplate'){
//            this.form2.setValue('copyplate', search_data);
//            // form2 may not be rendered yet, store the value
//            this.form2_data = search_data;
//            this.searchID = searchID;
//          }else if (uiResourceId == 'reagent') {
//            this.form3.setValue('well', search_data);
//            this.searchID = searchID;
//          }
//        }
//      }else{
//        console.log('Search box requires the "search/[searchID]" param');
//        return;
//      }
//    } 
    
//    afterRenderOld: function() {
//      var self=this;
//
//      $('#search-box-1').html(this.form1.render().$el);
//      
//// Removed - client side plate-copy-search parsing
////      var form2 = this.form2;
////      var $form2 = this.form2.render().$el;
////      $('#search-box-2').html($form2);
////      $form2.append([
////          '<button type="submit" class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
////          ].join(''));
////
////      $form2.find('[ type="submit" ]').click(function(e){
////        e.preventDefault();
////        var errors = form2.commit({ validate: true }); 
////        if(!_.isEmpty(errors)){
////          $form2.find('#copyplate').addClass(self.errorClass);
////          return;
////        }else{
////          $form2.find('#copyplate').removeClass(self.errorClass);
////        }
////        var copy_plate_search_entry = self.form2.getValue()['copyplate'];
////        if (copy_plate_search_entry) {
////          var search_array = self.parse_search_val(copy_plate_search_entry);
////          if (search_array.length > appModel.MAX_SEARCH_ARRAY_SIZE ){
////            appModel.showModalMessage({
////              title: 'Please try your search again',
////              body: Iccbl.formatString(
////                appModel.MSG_SEARCH_SIZE_EXCEEDED, {
////                  size_limit: appModel.MAX_SEARCH_ARRAY_SIZE,
////                  actual_size: search_array.length
////                })
////            });
////            return;
////          }
////          matching_hash = {
////            matching_order: ['plate_number','copy_name'],
////            copy_name: {
////              title: 'Copy name',
////              pattern: /\b([a-zA-Z0-9]*[a-zA-Z][a-zA-Z0-9]*)\b/,
////              // pattern: /(\w+)\s+((\d+)-(\d+))/
////              help: 'letters and numbers with no spaces'            
////            },
////            plate_number: {
////              title: 'Plate number',
////              pattern: /^\d{1,5}$/,
////              help: 'a 1-5 digit number'
////            },
////          };
////          var searches_parsed = self.parse_searches(search_array, matching_hash);
////          var data_el = $form2.find('#data-error');
////          if (searches_parsed.errors) {
////            data_el.html('<div class="help-block"><h5>Search term not recognized:</h5>"' + 
////              searches_parsed.errors.join('" ,"') + '"</br>' +
////                "<h5>Allowed patterns:</h5>" + 
////                _.map(matching_hash, function(val,key){
////                  if (val.title){
////                    return '<dt>' + val.title + '</dt><dd>' + val.help + '</dd>'; 
////                  }
////                }).join('') + "</dl></div>");
////            data_el.show();
////            return;
////          }else{
////            data_el.hide();
////          }
////          self.submit_searches('librarycopyplate', searches_parsed.or_list);
////        }
////      });
//      
//      ///// CopyPlate search 2 - perform search on server
//      var form2 = this.form2;
//      var $form2 = this.form2.render().$el;
//      $('#search-box-2').html($form2);
//      $form2.append([
//          '<button type="submit" class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
//          ].join(''));
//
//      $form2.find('[ type="submit" ]').click(function(e){
//        e.preventDefault();
//        var errors = form2.commit({ validate: true }); 
//        if(!_.isEmpty(errors)){
//          console.log('form2 errors, abort submit: ' + JSON.stringify(errors));
//          $form2.find('[name="copyplate"]').addClass(self.errorClass);
//          return;
//        }else{
//          $form2.find('[name="copyplate"]').removeClass(self.errorClass);
//        }
//        var text_to_search = self.form2.getValue()['copyplate'];
//        // must change the route, and create a post
//        
//        var resource = appModel.getResource('librarycopyplate');
//        // TODO: key the search data using the searchID: 
//        // this allows for bfowser "back" in the session
//        // will also need to listen to URIStack changes and grab the data
//        // from the search ID
//        var searchID = ( new Date() ).getTime();
//        appModel.setSearch(searchID,text_to_search);
//        this.searchID = searchID;
//        appModel.set('routing_options', {replace: false});  
//        var _route = resource.key + '/search/'+ searchID;
//        appModel.router.navigate(_route, {trigger:true});
//        // var newStack = [resource.key,'search',searchID];
//        // NOTE: easier to control the router history using navigate: 
//        // when using uristack, there is the problem of who set appModel.routing_options last:
//        // a race condition is set up between list2.js and search_box.js
//        //        appModel.set({'uriStack': newStack});     
//      });      
//
//// REMOVED: client side well-search parsing
//      ///// Well
////      var form3 = this.form3;
////      var $form3 = this.form3.render().$el;
////      $('#search-box-3').html($form3);
////      $form3.append([
////          '<button type="submit" class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
////          ].join(''));
////
////      $form3.find('[ type="submit" ]').click(function(e){
////        e.preventDefault();
////        var errors = form3.commit({ validate: true }); 
////        if(!_.isEmpty(errors)){
////          console.log('form3 errors, abort submit: ' + JSON.stringify(errors));
////          $form3.find('#well').addClass(self.errorClass);
////          return;
////        }else{
////          $form3.find('#well').removeClass(self.errorClass);
////        }
////        var well_search_entry = self.form3.getValue()['well'];
////        if (well_search_entry) {
////          var search_array = self.parse_search_val(well_search_entry);
////          if (search_array.length > appModel.MAX_SEARCH_ARRAY_SIZE ){
////            appModel.showModalMessage({
////              title: 'Please try your search again',
////              body: Iccbl.formatString(
////                appModel.MSG_SEARCH_SIZE_EXCEEDED, {
////                  size_limit: appModel.MAX_SEARCH_ARRAY_SIZE,
////                  actual_size: search_array.length
////                })
////            });
////            return;
////          }
////          matching_hash = {
////            matching_order: ['plate_number','well_id','well_name'],
////            plate_number: {
////              title: 'Plate number',
////              pattern: /^\d{1,5}$/,
////              help: 'a 1-5 digit number'
////            },
////            well_id: {
////              title: 'Well ID',
////              pattern: /(\d{1,5})\:([a-zA-Z]{1,2})(\d{1,2})/,
////              help: '[plate_number]:[well_name]',
////              parser: function(val_list){
////                return _.map(val_list, function(term){
////                  var match = matching_hash.well_id.pattern.exec(term);
////                  if (! match) return null;
////                  var plate_number = match[1];
////                  if (plate_number.length < 5){
////                    var pad = "00000";
////                    (pad+plate_number).substring(plate_number.length);
////                  }
////                  var well_row = match[2].toUpperCase();
////                  var well_col = match[3];
////                  // check the col, 2 digits
////                  if(well_col.length == 1) well_col = '0' + well_col;
////                  return plate_number + ':' + well_row + well_col;
////                });
////              }
////            },
////            well_name: {
////              title: 'Well name',
////              pattern: /([a-zA-Z]{1,2})(\d{1,2})/,
////              help: '[row][col] in the form A01,A12,B01 etc.',
////              parser: function(val_list){
////                return _.map(val_list, function(term){
////                  var rowcol = matching_hash.well_name.pattern.exec(term);
////                  if (! rowcol) return null;
////                  // check the col, 2 digits
////                  if(rowcol[2].length == 1) term = rowcol[1] + '0'+rowcol[2];
////                  return term.toUpperCase();
////                });
////              }
////            }
////          };
////          var searches_parsed = self.parse_searches(search_array, matching_hash);
////          var data_el = $form3.find('#data-error');
////          if (searches_parsed.errors) {
////            var data_el = $form3.find('#data-error');
////            data_el.html('<div class="help-block"><h5>Search term not recognized:</h5>"' + 
////              searches_parsed.errors.join('" ,"') + '"</br>' +
////                "<h5>Allowed patterns:</h5>" + 
////                _.map(matching_hash, function(val,key){
////                  if (val.title){
////                    return '<dt>' + val.title + '</dt><dd>' + val.help + '</dd>'; 
////                  }
////                }).join('') + "</dl></div>");
////            data_el.show();
////            return;
////          }else{
////            data_el.hide();
////          }
////          self.submit_searches('reagent', searches_parsed.or_list);
////        }
////      });
//
//      ///// Well/Reagent search 2 - perform search on server
//      var form3 = this.form3;
//      var $form3 = this.form3.render().$el;
//      $('#search-box-3').html($form3);
//      $form3.append([
//          '<button type="submit" class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
//          ].join(''));
//
//      $form3.find('[ type="submit" ]').click(function(e){
//        e.preventDefault();
//        var errors = form3.commit({ validate: true }); 
//        if(!_.isEmpty(errors)){
//          console.log('form3 errors, abort submit: ' + JSON.stringify(errors));
//          $form3.find('[name="well"]').addClass(self.errorClass);
//          return;
//        }else{
//          $form3.find('[name="well"]').removeClass(self.errorClass);
//        }
//        var text_to_search = self.form3.getValue()['well'];
//        // must change the route, and create a post
//        
//        var resource = appModel.getResource('reagent');
//        // TODO: key the search data using the searchID: 
//        // this allows for bfowser "back" in the session
//        // will also need to listen to URIStack changes and grab the data
//        // from the search ID
//        var searchID = ( new Date() ).getTime();
//        appModel.setSearch(searchID,text_to_search);
//        this.searchID = searchID;
//        appModel.set('routing_options', {replace: false});  
//        var _route = resource.key + '/search/'+ searchID;
//        appModel.router.navigate(_route, {trigger:true});
//        // var newStack = [resource.key,'search',searchID];
//        // NOTE: easier to control the router history using navigate: 
//        // when using uristack, there is the problem of who set appModel.routing_options last:
//        // a race condition is set up between list2.js and search_box.js
//        //        appModel.set({'uriStack': newStack});     
//      });      
//      
//      ///// CPR
//      var form4 = this.form4;
//      var $form4 = this.form4.render().$el;
//      $('#search-box-4').html($form4);
//      $form4.append([
//          '<button type="submit" class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
//          ].join(''));
//
//      $form4.find('[ type="submit" ]').click(function(e){
//        e.preventDefault();
//        var errors = form4.commit({ validate: true }); 
//        if(!_.isEmpty(errors)){
//          console.log('form4 errors, abort submit: ' + JSON.stringify(errors));
//          $form4.find('#cpr').addClass(self.errorClass);
//          return;
//        }else{
//          $form4.find('#cpr').removeClass(self.errorClass);
//        }
//        var cpr_id = self.form4.getValue()['cpr'];
//        var resource = appModel.getResource('cherrypickrequest');
//        var _route = ['#', resource.key,cpr_id].join('/');
//        appModel.set('routing_options', {replace: false});  
//        appModel.router.navigate(_route, {trigger:true});
//      });      
//      
//      console.log('setup single selects using chosen...');
//      // See http://harvesthq.github.io/chosen/
//      this.$el.find('.chosen-select').chosen({
//        width: '100%',
//        allow_single_deselect: true,
//        search_contains: true
//        });
//    },
    
// 20170410 - moved copyplate and well search to the server
//    parse_searches: function(search_array, matching_hash){
//      var errors = [];
//      var or_list = [];
//      _.each(search_array, function(and_array){
//        var and_hash = {};
//        or_list.push(and_hash);
//        var found_matches = [];
//        _.each(and_array, function(terms){
//          // find out what kind of terms they are
//          var clause = ''
//          var type = 'in';
//          var entity = '';
//          
//          var testing_val = terms[0];
//          if(testing_val.indexOf('-')>0){
//            type='range';
//            var testing_val = testing_val.split('-')[0];
//          }
//          entity = _.find(matching_hash.matching_order,function(key){
//            if(matching_hash[key].pattern.exec(testing_val)){
//              if (!_.contains(found_matches,key)){
//                found_matches.push(key);
//                return true;
//              }
//            }
//          });
//          
//          if(_.isEmpty(entity)){
//            errors.push(testing_val);
//          }else{
//            
//            if(matching_hash[entity].parser){
//              terms = matching_hash[entity].parser(terms);
//            }
//            
//            if(matching_hash[entity].type) type = matching_hash[entity].type;
//            
//            var search_entry = entity + '__' + type;
//            
//            if(type == 'range'){
//               number_array = terms[0].split('-');
//               _.each(number_array,function(num_val){
//                 if(_.isNaN(parseInt(num_val))){
//                   errors.push('ranges must be numbers: ' + terms);
//                 }
//               });
//               and_hash[search_entry] = number_array.join(',');
//            }else{
//              and_hash[search_entry] = terms.join(',');
//            }
//          }
//        });
//      });
//      var result = { 'or_list': or_list };
//      if (!_.isEmpty(errors)){
//        result['errors'] = errors;
//      }
//      return result;
//    },

// 20170410 - moved copyplate and well search to the server
//    /**
//     * Process search data as an array of AND clauses that will be OR'd together
//     * FIXME: if the search_data are large, implement server side searching,
//     * as in the reagent search. 201703/sde
//     */
//    submit_searches: function(search_target, or_list){
//      console.log('submit search: ', search_target, or_list);
//      if(or_list && or_list.length > 1){
//        // then this a composite/or search
//        var searchID = ( new Date() ).getTime();
//        
//        // must change the route, and create a post
//        
//        var resource = appModel.getResource(search_target);
//        // TODO: key the search data using the searchID: 
//        // this allows for bfowser "back" in the session
//        // will also need to listen to URIStack changes and grab the data
//        // from the search ID
//        var _serialized = JSON.stringify(or_list);
//        appModel.setSearch(searchID,or_list);
//        this.searchID = searchID;
//        appModel.set('routing_options', {replace: false});  
//        var _route = search_target + '/search/'+ searchID;
//        appModel.router.navigate(_route, {trigger:true});
//        // var newStack = [resource.key,'search',searchID];
//        // NOTE: easier to control the router history using navigate: 
//        // when using uristack, there is the problem of who set appModel.routing_options last:
//        // a race condition is set up between list2.js and search_box.js
//        //        appModel.set({'uriStack': newStack});     
//      }else{
//        // this is a simple "AND" search -which be translated to the list for handling
//        var _route = search_target + '/search/';
//        var routes = [];
//        _.each(_.pairs(or_list[0]),function(pair){
//          var val = pair[1];
//          if(_.isArray(val)){
//            val = val.join(',');
//          }
//          routes.push(pair[0] + '=' + val);
//        });
//        _route += routes.join(';')
//        console.log('route', _route);
//        appModel.set('routing_options', {replace: false});  
//        appModel.router.navigate(_route, {trigger:true});
//      }
//    },
    
// 20170410 - moved copyplate and well search to the server
//    /**
//     * parse:
//     * - newline for "or"
//     * - space for "and"
//     * - comma for "in"
//     * creates search_array (example):
//     * [ [ [1124,1123],[A01,A02,D12,G14] ],
//     *   [ [1200-1300],[A01,A03] ]
//     * ]
//     * will translated and sent to server as:
//     * [ 
//     *   { plate_number__range=[1222,1332], well_id__in=[b02,c02,d03,N12] },
//     *   { plate_number__range=[1222,1332], well_id__in=[b02,c02,d03,N12] }
//     * ]
//     */
//    parse_search_val: function(val){
//      var search_array = []
//      var or_list = val.split('\n');
//      _.each(or_list, function(clause){
//        clause = clause.trim();
//        if(clause=='') return;
//        // split on space, but not space followed by a comma:
//        // (reverse) the string to do (negative) look ahead to exclude ',\s' 
//        // do this because JavaScript doesn't support negative lookbehind
//        clause = clause.split('').reverse().join('');
//        var and_list = clause.split(/\s+(?!,)/);
//        var and_array = [];
//        search_array.push(and_array);
//        _.each(and_list, function(reversedterms){
//          // now reverse each term to put it back in proper non-reversed form
//          var terms = reversedterms.split('').reverse().join('');
//          and_array.unshift(terms.split(/[,\s]+/));
//        });
//      });
//      return search_array;
//    },
    
    
  });

  return SearchView;
});