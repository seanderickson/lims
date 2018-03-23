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
    
//  var PARAM_RAW_SEARCH = appModel.API_PARAM_RAW_SEARCH;
  
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
      function validateWellSearch(value,formValues){
        var errors = [];
        var parsedData = Iccbl.parseRawWellSearch(value, errors);
        if (_.isEmpty(parsedData)){
          errors.push('no values found for input');
        } else {
          console.log('parsedData', parsedData);
        }
        if (!_.isEmpty(errors)){
          return {
            type: 'searchVal',
            message: errors.join('; ')
          };
        }
      };
      
      var schema3 = {};
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
            type: 'searchVal',
            message: errors.join('; ')
          };
        }
      };
      schema3['searchVal'] = {
        title: 'Search for Plate & Copy',
        key: 'searchVal',
        help: 'enter a comma separated list',
        placeholder: 'e.g. 1000-1005 B\n2000 C',
        validators: ['required',validatePlateSearch],
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
      
      ///// Well search
      var schema2 = {};
      schema2['searchVal'] = {
        title: 'Search for 384-well Plate Wells',
        key: 'searchVal',
        help: 'enter a comma separated list',
        placeholder: 'e.g. 1000 A08,A09,A10\n1740:D12',
        validators: ['required',validateWellSearch],
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
      
      ///// Compound Name, Facility or Vendor ID
      var schema2a = {};
      schema2a['searchVal'] = {
        title: 'Compound Name or Vendor ID',
        key: 'searchVal',
        help: 'enter each compound name or vendor ID on a separate line',
        placeholder: 'enter each compound name or vendor ID on a separate line',
        validators: ['required'],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var FormFields2a = Backbone.Model.extend({
        schema: schema2a,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var formFields2a = new FormFields2a({
      });
      
      var form2a = self.form2a = new Backbone.Form({
        model: formFields2a,
        template: self.formTemplate,
        el: '#search-box-2a'
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
            type: 'searchVal',
            message: errors.join('; ')
          };
        }
      };
      schema5['searchVal'] = {
        title: 'Screening Inquiry',
        key: 'searchVal',
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
      schema4['searchVal'] = {
        title: 'CPR #',
        key: 'searchVal',
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
      if (this.form2_data){
        this.form2.setValue('searchVal', this.form2_data);
      }
      if (this.form2a_data){
        this.form2a.setValue('searchVal', this.form2a_data);
      }
      if (this.form3_data){
        this.form3.setValue('searchVal', this.form3_data);
      }
      if (this.form5_data){
        this.form5.setValue('searchVal', this.form5_data);
      }
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
      
      ///// CopyPlate search 3 - perform search on server
      var form3 = this.form3;
      var $form3 = this.form3.render().$el;
      $('#search-box-3').html($form3);
      $form3.append([
          '<button type="submit" ',
          'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
          ].join(''));

      $form3.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        self.$('[data-error]').empty();
        var errors = form3.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          $form3.find('[name="searchVal"]').addClass(self.errorClass);
          return;
        }else{
          $form3.find('[name="searchVal"]').removeClass(self.errorClass);
        }
        var errors = self.processCopyPlateSearch(self.form3.getValue('searchVal'));
        if (!_.isEmpty(errors)){
          throw 'Unexpected errors after submit: ' + errors.join(', ');
        }
          
      });      
      
      ///// Well/Reagent search 2 - perform search on server
      var form2 = this.form2;
      var $form2 = this.form2.render().$el;
      $('#search-box-2').html($form2);
      $form2.append([
          '<button type="submit" ',
          'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
          ].join(''));

      $form2.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        self.$('[data-error]').empty();
        var errors = form2.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          $form2.find('[name="searchVal"]').addClass(self.errorClass);
          return;
        }else{
          $form2.find('[name="searchVal"]').removeClass(self.errorClass);
        }
        
        var errors = self.processWellSearch(self.form2.getValue('searchVal'));
        if (!_.isEmpty(errors)){
          $form2.find('[name="searchVal"]').addClass(self.errorClass);
          $form2.find('[data-error]').html(errors.join('<br/>'));
          return;
        }
      });
      var $help = $('<a>', {
        title: 'Help'
      }).text('?');
      $help.click(function(e){
        e.preventDefault();
        appModel.showModalMessage({
          title: 'Searching for wells',
          okText: 'Ok',
          body: [
            '<b>Search for single wells:</b>',
            '50001:P24, 50002:A01, 50002:A23',
            '<b>Equivalent search:</b>',
            '50001P24, 50002A1, 50002A23',
            '<b>Search for multiple wells in a plate:</b>',
            '50001P24 A1 A23',
            '50001 P24 A1 A23',
            '50001 P24,A1,A23',
            '<b>Search for multiple wells in a plate range:</b>',
            '50001-50020 A1,A2,A3',
            '50001-50020 A1 A2 A3',
            '<b>Search for many wells:</b>',
            '50001:P21','50001:P22','50001:P23','50001:P24',
            '<b>Search for all wells in a plate or plate range:</b>',
            '50001 50002 50003',
            '50010-50020',
          ].join('<br/>'),
          ok: function(e){
            e.preventDefault();
          }
        });
      });
      $form2.find('[key="form-group-searchVal"]').find('label').append('&nbsp;').append($help);

      ///// Compound name, Vendor ID search- perform search on server
      var form2a = this.form2a;
      var $form2a = this.form2a.render().$el;
      $('#search-box-2a').html($form2a);
      $form2a.append([
          '<button type="submit" ',
          'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
          ].join(''));

      $form2a.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        self.$('[data-error]').empty();
        var errors = form2a.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          $form2a.find('[name="searchVal"]').addClass(self.errorClass);
          return;
        }else{
          $form2a.find('[name="searchVal"]').removeClass(self.errorClass);
        }
        var errors = self.processCompoundSearch(self.form2a.getValue()['searchVal']);
        if (!_.isEmpty(errors)){
          throw 'Unexpected errors after submit: ' + errors.join(', ');
        }
      });      
      
      ///// Screening Inquiry
      if (appModel.hasPermission('libraryscreening', 'read')){
        var form5 = this.form5;
        var $form5 = this.form5.render().$el;
        $('#search-box-5').html($form5);
        $form5.append([
            '<button type="submit" ',
            'class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
            ].join(''));
  
        $form5.find('[ type="submit" ]').click(function(e){
          e.preventDefault();
          self.$('[data-error]').empty();
          var errors = form5.commit({ validate: true }); 
          if(!_.isEmpty(errors)){
            $form5.find('[name="searchVal"]').addClass(self.errorClass);
            return;
          }else{
            $form5.find('[name="searchVal"]').removeClass(self.errorClass);
          }
  
          var searchValue = form5.getValue('searchVal');
          var errors = [];
          var parsedSearch = Iccbl.parseRawScreeningInquiry(searchValue,errors);
          if (!_.isEmpty(errors)){
            throw 'Unexpected errors after submit: ' + errors.join(', ');
          }
          var urlSearchParts = 
            PlateRangeSearchView.prototype.encodeFormData.call(this,parsedSearch);
          // TODO: 20180312 - screening inquiry should use URI_PATH_ENCODED_SEARCH:
          // rework to support list args
          var uriStack = ['screen', parsedSearch.screen_facility_id,
                          'summary','plateranges',appModel.URI_PATH_SEARCH,
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
          self.$('[data-error]').empty();
          var errors = form4.commit({ validate: true }); 
          if(!_.isEmpty(errors)){
            console.log('form4 errors, abort submit: ' + JSON.stringify(errors));
            $form4.find('#cpr').addClass(self.errorClass);
            return;
          }else{
            $form4.find('#cpr').removeClass(self.errorClass);
          }
          var cpr_id = self.form4.getValue()['searchVal'];
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
      
      console.log('uriStackChange', arguments, uriStack);
      var uiResourceId = uriStack.shift();

      var complex_search = appModel.findComplexSearch(uriStack);
      console.log('complex search', complex_search);
      if (!_.isEmpty(complex_search)){
        if (_.has(complex_search,'search_id') 
            && complex_search['search_id']== this.searchId){
          console.log('self generated uristack change');
          return;
        }
        if (_.has(complex_search,'errors')){
          console.log('complex search errors: ', complex_search['errors']);
          return;
        }
        
        console.log('complex search found:', complex_search);
        if (uiResourceId=='librarycopyplate'){
          var errors = [];
          var search_data = complex_search[appModel.API_PARAM_SEARCH];
          var parsedData = Iccbl.parseRawPlateSearch(search_data,errors);
          if (!_.isEmpty(errors)){
            console.log('Search data not parsed properly', errors);
            return;
          }
          parsedData = _.map(parsedData, function(parsedLine){
            return parsedLine.combined.join(' ');
          }).join('\n');
          this.form3.setValue('searchVal', parsedData);
          // form3 may not be rendered yet, store the value
          this.form3_data = parsedData;
        }else if (_.contains(
          ['silencingreagent','smallmoleculereagent','reagent'], uiResourceId)) {
          var errors = [];
          var search_data = complex_search[appModel.API_PARAM_SEARCH];
          var parsedData = Iccbl.parseRawWellSearch(search_data,errors);
          if (!_.isEmpty(errors)){
            console.log('Search data not parsed properly', errors);
            return;
          }
          parsedData = _.map(parsedData, function(parsedLine){
            return parsedLine.combined.join(' ');
          }).join('\n');
          this.form2.setValue('searchVal', parsedData);
          // form3 may not be rendered yet, store the value
          this.form2_data = parsedData;
        }else if (uiResourceId == 'compound_search') {
          var search_data = complex_search[appModel.API_PARAM_SEARCH];
          var parsedData = Iccbl.parseCompoundVendorIDSearch(search_data,errors);
          if (!_.isEmpty(errors)){
            console.log('Search data not parsed properly', errors);
            return;
          }
          parsedData = parsedData.join('\n');
          this.form2a.setValue('searchVal', parsedData);
          this.form2a_data = parsedData;
        } else {
          throw 'unknown resource for search: ' + uiResourceId;
        }
        return;
      }
      else if (self.form5 && uiResourceId == 'screen') {
        var screenId = uriStack.shift();
        if (_.contains(uriStack, 'plateranges')){
          // TODO: 20180312 - screening inquiry: using URI_PATH_SEARCH: 
          // convert to URI_PATH_ENCODED_SEARCH; rework to enable list args
          if (_.contains(uriStack,appModel.URI_PATH_SEARCH)){
            var index = _.indexOf(uriStack,appModel.URI_PATH_SEARCH);
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
            self.form5.setValue('searchVal',parsedData);
            // form2 may not be rendered yet, store the value
            self.form5_data = parsedData;
          }
        }
        return;
      } // not screen/plate ranges

    },
    
    processCopyPlateSearch: function(text_to_search){
      errors = [];
      var parsedSearchArray = Iccbl.parseRawPlateSearch(text_to_search,errors);
      if (!_.isEmpty(errors)){
        return errors;
      }
      var resource = appModel.getResource('librarycopyplate');
      if (parsedSearchArray.length <= appModel.MAX_RAW_SEARCHES_IN_URL){
        // encode simple searches as a URL param
        var encodedPlateSearches = [];
        _.each(parsedSearchArray, function(parsedPlateSearch){
          encodedPlateSearches.push(_.map(
            parsedPlateSearch.combined, encodeURIComponent).join(','));
        });
        var uriStack = [resource.key, appModel.URI_PATH_ENCODED_SEARCH, 
          encodedPlateSearches.join(appModel.UI_PARAM_RAW_SEARCH_LINE_ENCODE)];
        console.log('route: ', uriStack);
        appModel.router.navigate(uriStack.join('/'), {trigger:true});
      }else{
        // Send complex search data as a POST
        // must change the route, and create a post
        // TODO: key the search data using the searchId: 
        // this allows for browser "back" in the session
        // will also need to listen to URIStack changes and grab the data
        // from the search ID
        var searchId = ( new Date() ).getTime();
        appModel.setSearch(searchId,text_to_search);
        this.searchId = searchId;
        appModel.set('routing_options', {replace: false});  
        var _route = [
          resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
          searchId].join('/');
        appModel.router.navigate(_route, {trigger:true});
      }
    },
    
    processWellSearch: function(text_to_search){
      errors = [];
      var parsedSearchArray = Iccbl.parseRawWellSearch(text_to_search,errors);
      if (!_.isEmpty(errors)){
        return errors;
      }
      var librariesForSearch = appModel.getLibrariesForSearch(parsedSearchArray, errors);
      if (!_.isEmpty(errors)){
        return errors;
      }
      
      var screenType = librariesForSearch[0].get('screen_type');
      var resource = appModel.getResource('reagent');
      
      if (screenType == 'small_molecule'){
        resource = appModel.getResource('smallmoleculereagent');
      } else if (screenType == 'rnai'){
        resource = appModel.getResource('silencingreagent')
      } else {
        throw 'Unknown library type for search: ' + screenType;
      }
      
      if (parsedSearchArray.length <= appModel.MAX_RAW_SEARCHES_IN_URL){
        // encode simple searches as a URL param
        var encodedSearches = [];
        _.each(parsedSearchArray, function(parsedSearch){
          encodedSearches.push(_.map(
            parsedSearch.combined, encodeURIComponent).join(','));
        });
        var uriStack = [resource.key, appModel.URI_PATH_ENCODED_SEARCH, 
          encodedSearches.join(appModel.UI_PARAM_RAW_SEARCH_LINE_ENCODE)];
        console.log('route: ', uriStack);
        appModel.router.navigate(uriStack.join('/'), {trigger:true});
      }else{
        // must change the route, and create a post
        var searchId = ( new Date() ).getTime();
        appModel.setSearch(searchId,text_to_search);
        this.searchId = searchId;
        appModel.set('routing_options', {replace: false});  
        // Use the router to process the search: see content.js for handler
        var resource = appModel.getResource('reagent');
        var _route = [
          resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
          searchId].join('/');
        appModel.router.navigate(_route, {trigger:true});
        // var newStack = [resource.key,appModel.URI_PATH_SEARCH,searchId];
        // NOTE: easier to control the router history using navigate: 
        // when using uristack, there is the problem of who set appModel.routing_options last:
        // a race condition is set up between list2.js and search_box.js
        //        appModel.set({'uriStack': newStack});     
      }
      
    },
    
    processCompoundSearch: function(text_to_search){
      errors = [];
      var parsedSearchArray = Iccbl.parseCompoundVendorIDSearch(text_to_search,errors);
      if (!_.isEmpty(errors)){
        return errors;
      }
      var resource = appModel.getResource('compound_search');
      if (parsedSearchArray.length <= appModel.MAX_RAW_SEARCHES_IN_URL){
        var encodedSearches = _.map(parsedSearchArray,encodeURIComponent);
        var uriStack = [resource.key, appModel.URI_PATH_ENCODED_SEARCH, 
          encodedSearches.join(appModel.UI_PARAM_RAW_SEARCH_LINE_ENCODE)];
        appModel.router.navigate(uriStack.join('/'), {trigger:true});
      }else{
        // must change the route, and create a post
        var parsedSearches = text_to_search.split('')
        var searchId = ( new Date() ).getTime();
        appModel.setSearch(searchId,text_to_search);
        this.searchId = searchId;
        appModel.set('routing_options', {replace: false});  
        // Use the router to process the search: see content.js for handler
        var resource = appModel.getResource('compound_search');
        var _route = [
          resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
          searchId].join('/');
        appModel.router.navigate(_route, {trigger:true});
      }
    },
    
    // Backbone layoutmanager callback
    cleanup: function(){
      this.form3.setValue('searchVal', null);
      this.form2.setValue('searchVal', null);
      this.form2a.setValue('searchVal', null);
      this.form4.setValue('searchVal', null);
      this.form5.setValue('searchVal', null)
    },
    
  });

  return SearchView;
});