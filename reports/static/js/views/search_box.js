define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_forms',
    'iccbl_backgrid',
    'models/app_state',
    'views/generic_edit',
    'views/screen/plateRangeSearch',
    'select2',
    'templates/search-box.html'
], function($, _, Backbone, backbone_forms, Iccbl, appModel, EditView, 
    PlateRangeSearchView, select2, searchBoxTemplate) {
    
  var SearchView = Backbone.Layout.extend({

    template: _.template(searchBoxTemplate),
    
    formTemplate: _.template([
        "<form class='iccbl-headerfield-form' ",
        ">",
        '<div data-fieldsets/>',
        '<div id="data-error" class="has-error" ></div>',
        "</form>",
        ].join('')),

    searchFormTemplate: _.template([
        '<form class="iccbl-headerfield-form" >',
        '  <div data-fieldsets class="col-xs-11" />',
        '  <div class="col-xs-1">',
        '    <button type="submit" class="btn btn-default pull-right btn-block"  >&gt;</button>',
        '  </div>',
        '  <div id="data-error" class="has-error" ></div>',
        "</form>",
        ].join('')),
    
    fieldTemplate: _.template([
        '<div key="form-field-<%=key%>" >',
        '  <span key="<%=key%>" title="<%= help %>" data-editor></span>',
        '  <div data-error class="col-sm-12 text-danger" ></div>',
        '</div>',
    ].join('')),
    
    choiceFieldTemplate: _.template([
      '<div class="form-group" key="form-group-<%=key%>" title="<%=title%>" >',
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
      
      var plateCopySchema = {};
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
      plateCopySchema['searchVal'] = {
        help: 'Plate and Copy Search',
        key: 'searchVal',
        placeholder: 'Plate & Copy...',
        validators: ['required',validatePlateSearch],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var PlateCopyFormModel = Backbone.Model.extend({
        schema: plateCopySchema,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var plateCopyFormModel = new PlateCopyFormModel({
      });
      
      var plateCopyForm = self.plateCopyForm = new Backbone.Form({
        model: plateCopyFormModel,
        template: self.searchFormTemplate,
        el: '#search-box-3'
      });
      
      ///// Well search

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
      
      var wellsSchema = {};
      wellsSchema['searchVal'] = {
        title: 'Find',
        key: 'searchVal',
        help: 'Plate Well Search',
        placeholder: 'Plate Wells...',
        validators: ['required',validateWellSearch],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var WellsFormModel = Backbone.Model.extend({
        schema: wellsSchema,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var wellsFormModel = new WellsFormModel({
      });
      
      var wellsForm = self.wellsForm = new Backbone.Form({
        model: wellsFormModel,
        template: self.searchFormTemplate,
        el: '#search-box-2-content'
      });
      
      // 20181029 - show the copywell form only from the reports menu
      ///// Copy Well search

      //function validateCopyWellSearch(value,formValues){
      //  var errors = [];
      //  var parsedData = Iccbl.parseRawCopyWellSearch(value, errors);
      //  if (_.isEmpty(parsedData)){
      //    errors.push('no values found for input');
      //  } else {
      //    console.log('parsedData', parsedData);
      //  }
      //  if (!_.isEmpty(errors)){
      //    return {
      //      type: 'searchVal',
      //      message: errors.join('; ')
      //    };
      //  }
      //};
      
      //var copyWellsSchema = {};
      //copyWellsSchema['searchVal'] = {
      //  title: 'Find',
      //  key: 'searchVal',
      //  help: 'Copy Well Search',
      //  placeholder: 'Copy Wells...',
      //  validators: ['required',validateCopyWellSearch],
      //  type: EditView.TextArea2,
      //  template: self.fieldTemplate,
      //  editorClass: 'form-control'
      //};
      //var CopyWellsFormModel = Backbone.Model.extend({
      //  schema: copyWellsSchema,
      //  validate: function(attrs) {
      //    var errors = {};
      //    if (!_.isEmpty(errors)) return errors;
      //  }
      //});
      //var copyWellsFormModel = new CopyWellsFormModel({
      //});
      //
      //var copyWellsForm = self.copyWellsForm = new Backbone.Form({
      //  model: copyWellsFormModel,
      //  template: self.searchFormTemplate,
      //  el: '#search-box-2b'
      //});
      
      ///// Compound Name, Facility or Vendor ID
      
      var compoundSchema = {};
      compoundSchema['searchVal'] = {
        help: 'Compound Name, Vendor ID Search',
        key: 'searchVal',
        placeholder: 'Compound Name or Vendor ID...',
        validators: ['required'],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var CompoundFormModel = Backbone.Model.extend({
        schema: compoundSchema,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var compoundFormModel = new CompoundFormModel({
      });
      
      var compoundNameForm = self.compoundNameForm = new Backbone.Form({
        model: compoundFormModel,
        template: self.searchFormTemplate,
        el: '#search-box-2a'
      });
      
      ///// Screening Inquiry
      
      var screeningInqurySchema = {};
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
      screeningInqurySchema['searchVal'] = {
        help: 'Screening Inquiry Search (enter Screen, Plates, Volume and Replicates desired)',
        key: 'searchVal',
        placeholder: 'Screening Inquiry...',
        validators: ['required',validateScreeningInquiry],
        type: EditView.TextArea2,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var ScreeningInquryFormModel = Backbone.Model.extend({
        schema: screeningInqurySchema,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var screeningInquiryModel = new ScreeningInquryFormModel({
      });
      
      var screeningInquiryForm = self.screeningInquiryForm = new Backbone.Form({
        model: screeningInquiryModel,
        template: self.searchFormTemplate,
        el: '#search-box-5'
      });      
      
      ///// Cherry Pick Request
      
      var cprSchema = {};
      cprSchema['searchVal'] = {
        help: 'Cherry Pick Request Search',
        key: 'searchVal',
        placeholder: 'Cherry Pick Request...',
        key:  'search_value', // TODO: "key" not needed>?
        type: EditView.NumberEditor,
        template: self.fieldTemplate,
        editorClass: 'form-control'
      };
      var CprFormModel = Backbone.Model.extend({
        schema: cprSchema,
        validate: function(attrs) {
          var errors = {};
          if (!_.isEmpty(errors)) return errors;
        }
      });
      var cprFormModel = new CprFormModel({
        cpr: null
      });
      
      var cprForm = self.cprForm = new Backbone.Form({
        model: cprFormModel,
        template: self.searchFormTemplate,
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
        title: 'Enter an eCommons ID (if available) or name of the user',
        key: 'user_id',
        type: Backbone.Form.editors.Select,
        editorClass: 'form-control user-select',
        options: appModel.getUserOptions(),
        template: self.choiceFieldTemplate 
      };
      schema1['screen_facility_id'] = {
        title: 'Enter the screen ID or text to search in the title',
        key: 'screen_facility_id',
        type: Backbone.Form.editors.Select,
        editorClass: 'form-control screen-select',
        options: appModel.getScreenOptions(),
        template: self.choiceFieldTemplate
      };
      
      var libraryOptions = appModel.getLibraryOptions();
      schema1['library_short_name'] = {
        title: 'Enter the library name or plate number',
        key: 'library_short_name',
        type: Backbone.Form.editors.Select,
        editorClass: 'form-control library-select',
        options: libraryOptions,
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
              var uriStack = ['screensaveruser', user_id];
              self.form1.setValue('user_id',null);
              self.form1.$el.find('[key="user_id"]')
                  .find('.chosen-select').trigger("chosen:updated");
              appModel.setUriStack(uriStack);
            }
            else if(screen_facility_id){
              var uriStack = ['screen',screen_facility_id];
              self.form1.setValue('screen_facility_id',null);
              self.form1.$el.find('[key="screen_facility_id"]')
                  .find('.chosen-select').trigger("chosen:updated");
              appModel.setUriStack(uriStack);
            }
            else if(library_short_name){
              var uriStack = ['library',library_short_name];
              self.form1.setValue('library_short_name',null);
              self.form1.$el.find('[key="library_short_name"]')
                  .find('.chosen-select').trigger("chosen:updated");
              appModel.setUriStack(uriStack);
            }
          }
        });
      });
      $('#search-box-1').html(this.form1.render().$el);

      //console.log('setup single selects using chosen...');
      //// See http://harvesthq.github.io/chosen/
      //this.$el.find('.chosen-select').chosen({
      //  width: '100%',
      //  allow_single_deselect: true,
      //  search_contains: true
      //  });
      
      // Select2
      function matchLibrary(params, data) {
        // User term: params.term
        // Option title: data.text
        
        if (_.isUndefined(data.text)){
          return null;
        }
        if (_.isUndefined(params.term)){
          return data;
        }
        var label = data.text.trim().toLowerCase();
        var term = params.term.trim().toLowerCase();
        
        if (_.isEmpty(label)) {
          return null;
        }
        if (_.isEmpty(term)){
          return data;
        }

        if (label.indexOf(term) > -1) {
          return data;
        }else{
          var val = parseInt(term);
          if (!_.isNaN(val)){
            // Match the special plate range text created for the options
            var PLATE_RANGE_PATTERN = /\((\d+)-(\d+)\)/;
            var rangeParts = PLATE_RANGE_PATTERN.exec(data.text);
            if (!_.isEmpty(rangeParts)){
              var start = parseInt(rangeParts[1]);
              var end = parseInt(rangeParts[2]);
              var result = val >= start && val <= end;
              if (result){
                return data;
              }
            }
          }
        }      
        return null;
      };
      function matchUser(params, data) {
        // User term: params.term
        // Option title: data.text
        if (_.isUndefined(data.text)){
          return null;
        }
        if (_.isUndefined(params.term)){
          return data;
        }
        var label = data.text.trim().toLowerCase();
        var term = params.term.trim().toLowerCase();
        
        if (_.isEmpty(label)) {
          return null;
        }
        if (_.isEmpty(term)){
          return data;
        }

        var parts = appModel.USER_OPTION_PATTERN.exec(label);
        if (!_.isEmpty(parts)){
          
          var name = parts[1];
          var username = parts[3];
          var ss_id = parseInt(parts[5]);
          var val = parseInt(params.term);
          if (!_.isNaN(val)){
            if (ss_id==val){
              return data;
            }
          }
          else if (username && username.indexOf(term)==0){
            return data;
          }else{
            var name = name.toLowerCase();
            var found = _.find(name.split(','), function(first_last){
              if (first_last){
                return first_last.trim().indexOf(term)==0;
              }
            });
            if (found) {
              return data;
            } else if (name.trim().indexOf(term) == 0){
              return data;
            }
          }
        }else if (label.indexOf(term) > -1) {
          return data;
        }
        return null;
      };
      self.$el.find('.library-select').select2({
        placeholder : { id:'-1', text:'Find a Library...'},
        matcher: matchLibrary,
        width: '100%'
      });

      self.$el.find('.screen-select').select2({
        placeholder : { id:'-1', text:'Find a Screen...'},
        width: '100%'
      });

      self.$el.find('.user-select').select2({
        placeholder : { id:'-1', text:'Find a User...'},
        matcher: matchUser,
        width: '100%'
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
      if (this.wells_form_data){
        this.wellsForm.setValue('searchVal', this.wells_form_data);
      }
      // 20181029 - show the copywell form only from the reports menu
      //if (this.copy_wells_form_data){
      //  this.copyWellsForm.setValue('searchVal', this.copy_wells_form_data);
      //}
      if (this.compound_name_form_data){
        this.compoundNameForm.setValue('searchVal', this.compound_name_form_data);
      }
      if (this.plateCopyForm_data){
        this.plateCopyForm.setValue('searchVal', this.plateCopyForm_data);
      }
      if (this.screening_inq_form_data){
        this.screeningInquiryForm.setValue('searchVal', this.screening_inq_form_data);
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
      var plateCopyForm = this.plateCopyForm;
      var $plateCopyForm = this.plateCopyForm.render().$el;
      $('#search-box-3').html($plateCopyForm);

      $plateCopyForm.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        self.$('[data-error]').empty();
        var errors = plateCopyForm.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          $plateCopyForm.find('[name="searchVal"]').addClass(self.errorClass);
          $plateCopyForm.find('[data-error]').html(errors.join('<br/>'));
          return;
        }else{
          $plateCopyForm.find('[name="searchVal"]').removeClass(self.errorClass);
        }
        var errors = self.processCopyPlateSearch(self.plateCopyForm.getValue('searchVal'));
        if (!_.isEmpty(errors)){
          throw 'Unexpected errors after submit: ' + errors.join(', ');
        }
          
      });      
      
      ///// Well/Reagent search 2 - perform search on server
      var wellsForm = this.wellsForm;
      var $wellsForm = this.wellsForm.render().$el;
      $('#search-box-2-content').html($wellsForm);

      $wellsForm.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        self.$('[data-error]').empty();
        var errors = wellsForm.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          $wellsForm.find('[name="searchVal"]').addClass(self.errorClass);
          return;
        }else{
          $wellsForm.find('[name="searchVal"]').removeClass(self.errorClass);
        }
        
        var errors = self.processWellSearch(self.wellsForm.getValue('searchVal'));
        if (!_.isEmpty(errors)){
          $wellsForm.find('[name="searchVal"]').addClass(self.errorClass);
          $wellsForm.find('[data-error]').html(errors.join('<br/>'));
          return;
        }
      });
      
      // 20181029 - show the copywell form only from the reports menu
      ///// CopyWell search - perform search on server
      //var copyWellsForm = this.copyWellsForm;
      //var $copyWellsForm = this.copyWellsForm.render().$el;
      //$('#search-box-2b').html($copyWellsForm);
      //
      //$copyWellsForm.find('[ type="submit" ]').click(function(e){
      //  e.preventDefault();
      //  self.$('[data-error]').empty();
      //  var errors = copyWellsForm.commit({ validate: true }); 
      //  if(!_.isEmpty(errors)){
      //    $copyWellsForm.find('[name="searchVal"]').addClass(self.errorClass);
      //    return;
      //  }else{
      //    $copyWellsForm.find('[name="searchVal"]').removeClass(self.errorClass);
      //  }
      //  
      //  var errors = self.processCopyWellSearch(self.copyWellsForm.getValue('searchVal'));
      //  if (!_.isEmpty(errors)){
      //    $copyWellsForm.find('[name="searchVal"]').addClass(self.errorClass);
      //    $copyWellsForm.find('[data-error]').html(errors.join('<br/>'));
      //    return;
      //  }
      //});
      
      ///// Compound name, Vendor ID search- perform search on server
      var compoundNameForm = this.compoundNameForm;
      var $compoundNameForm = this.compoundNameForm.render().$el;
      $('#search-box-2a').html($compoundNameForm);
      $compoundNameForm.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        self.$('[data-error]').empty();
        var errors = compoundNameForm.commit({ validate: true }); 
        if(!_.isEmpty(errors)){
          $compoundNameForm.find('[name="searchVal"]').addClass(self.errorClass);
          return;
        }else{
          $compoundNameForm.find('[name="searchVal"]').removeClass(self.errorClass);
        }
        var errors = self.processCompoundSearch(self.compoundNameForm.getValue()['searchVal']);
        if (!_.isEmpty(errors)){
          throw 'Unexpected errors after submit: ' + errors.join(', ');
        }
      });      
      
      ///// Screening Inquiry
      if (appModel.hasPermission('libraryscreening', 'read')){
        var screeningInquiryForm = this.screeningInquiryForm;
        var $screeningInquiryForm = this.screeningInquiryForm.render().$el;
        $('#search-box-5').html($screeningInquiryForm);
        $screeningInquiryForm.find('[ type="submit" ]').click(function(e){
          e.preventDefault();
          self.$('[data-error]').empty();
          var errors = screeningInquiryForm.commit({ validate: true }); 
          if(!_.isEmpty(errors)){
            $screeningInquiryForm.find('[name="searchVal"]').addClass(self.errorClass);
            return;
          }else{
            $screeningInquiryForm.find('[name="searchVal"]').removeClass(self.errorClass);
          }
  
          var searchValue = screeningInquiryForm.getValue('searchVal');
          var errors = [];
          var parsedSearch = Iccbl.parseRawScreeningInquiry(searchValue,errors);
          if (!_.isEmpty(errors)){
            throw 'Unexpected errors after submit: ' + errors.join(', ');
          }
          var urlSearchParts = 
            PlateRangeSearchView.prototype.encodeFormData.call(this,parsedSearch);
          // rework to support list args
          var uriStack = ['screen', parsedSearch.screen_facility_id,
                          'summary','plateranges',appModel.URI_PATH_SEARCH,
                          urlSearchParts.join(appModel.SEARCH_DELIMITER)];
          console.log('route: ', uriStack);
          appModel.setUriStack(uriStack);
        });      
      }      

      ///// CPR
      if (appModel.hasPermission('cherrypickrequest', 'read')){
        var cprForm = this.cprForm;
        var $cprForm = this.cprForm.render().$el;
        $('#search-box-4').html($cprForm);
        $cprForm.find('[ type="submit" ]').click(function(e){
          e.preventDefault();
          self.$('[data-error]').empty();
          var errors = cprForm.commit({ validate: true }); 
          if(!_.isEmpty(errors)){
            console.log('cprForm errors, abort submit: ' + JSON.stringify(errors));
            $cprForm.find('#cpr').addClass(self.errorClass);
            return;
          }else{
            $cprForm.find('#cpr').removeClass(self.errorClass);
          }
          var cpr_id = self.cprForm.getValue()['searchVal'];
          var resource = appModel.getResource('cherrypickrequest');
          appModel.getModelFromResource(resource, cpr_id, function(model){
            var uriStack = ['screen', model.get('screen_facility_id'), 
              resource.key,cpr_id];
              appModel.setUriStack(uriStack);
          });
        });      
      }

      $('#search-help').click(function(e){
        appModel.showModalMessage({
          title: 'Searching',
          okText: 'Ok',
          body: [
            '<br/><strong>Plate Well searching:</strong><br/>',
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
            '<br/><strong>Plate and Copy searching:</strong><br/>',
            '<b>Search for copies and plate range:</b>',
            '50001-50005 A',
            '50001-50005 A B C D',
            '50001,50002,50004 A B C D',
            '<br/><strong>Screening Inquiry (find available Plate Copies)</strong><br/>',
            'Syntax: &quot;(Screen #) &lt;plate numbers and ranges&gt; &lt;vol required&gt; nL x &lt;replicate count&gt;&quot;',
            'e.g.',
            '(1409) 3411-3414, 3448-3467, 3580 100 nL x 2',
            '(1292) 3560-3567, 1795-1811 100 nL x 2',
            'Note: text is ignored at the beginning:',
            'Tester Name (1292) 3535-3559 100 nL x 2'
          ].join('<br/>'),
          ok: function(e){
            e.preventDefault();
          }
        });
      });

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
      
      // For small viewports, the collapsed navbar has been shown by user action:
      // re-hide and allow the css to add the "show" class for larger viewports
      $('.navbar-collapse').removeClass('show');
      
      var uriStack = appModel.get('uriStack');
      if (!uriStack) return;
      var uriStack = _.clone(uriStack);
      var uiResourceId = uriStack.shift();
      var complex_search = appModel.findComplexSearch(uriStack);

      // Set the appropriate search box if there is a complex search in effect
      if (!_.isEmpty(complex_search)){
        if (_.has(complex_search,'search_id') 
            && complex_search['search_id']== this.searchId){
          return;
        }
        if (_.has(complex_search,'errors')){
          console.log('complex search errors: ', complex_search['errors']);
          return;
        }
        
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
          this.plateCopyForm.setValue('searchVal', parsedData);
          // plateCopyForm may not be rendered yet, store the value
          this.plateCopyForm_data = parsedData;
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
          this.wellsForm.setValue('searchVal', parsedData);
          // plateCopyForm may not be rendered yet, store the value
          this.wells_form_data = parsedData;
        }else if (uiResourceId == 'compound_search') {
          var search_data = complex_search[appModel.API_PARAM_SEARCH];
          var parsedData = Iccbl.parseCompoundVendorIDSearch(search_data,errors);
          if (!_.isEmpty(errors)){
            console.log('Search data not parsed properly', errors);
            return;
          }
          parsedData = parsedData.join('\n');
          this.compoundNameForm.setValue('searchVal', parsedData);
          this.compound_name_form_data = parsedData;
        } else {
          console.log('unknown resource for search: ' + uiResourceId);
        }
        return;
      }
      else if (self.screeningInquiryForm && uiResourceId == 'screen') {
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
            self.screeningInquiryForm.setValue('searchVal',parsedData);
            // wellsForm may not be rendered yet, store the value
            self.screening_inq_form_data = parsedData;
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
        appModel.setUriStack(uriStack);
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
        var uriStack = [
          resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
          searchId];
        appModel.setUriStack(uriStack);
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
        appModel.setUriStack(uriStack);
      }else{
        // must change the route, and create a post
        var searchId = ( new Date() ).getTime();
        appModel.setSearch(searchId,text_to_search);
        this.searchId = searchId;
        var uriStack = [
          resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
          searchId];
        appModel.setUriStack(uriStack);
      }
      
    },
    
    // 20181029 - show the copywell form only from the reports menu
    //processCopyWellSearch: function(text_to_search){
    //  console.log('processCopyWellSearch: ', text_to_search);
    //  errors = [];
    //  var parsedSearchArray = Iccbl.parseRawCopyWellSearch(text_to_search,errors);
    //  if (!_.isEmpty(errors)){
    //    return errors;
    //  }
    //
    //  var resource = appModel.getResource('copywell');
    //  
    //  if (parsedSearchArray.length <= appModel.MAX_RAW_SEARCHES_IN_URL){
    //    // encode simple searches as a URL param
    //    var encodedSearches = [];
    //    _.each(parsedSearchArray, function(parsedSearch){
    //      encodedSearches.push(_.map(
    //        parsedSearch.combined, encodeURIComponent).join(','));
    //    });
    //    var uriStack = [resource.key, appModel.URI_PATH_ENCODED_SEARCH, 
    //      encodedSearches.join(appModel.UI_PARAM_RAW_SEARCH_LINE_ENCODE)];
    //    appModel.setUriStack(uriStack);
    //  }else{
    //    // must change the route, and create a post
    //    var searchId = ( new Date() ).getTime();
    //    appModel.setSearch(searchId,text_to_search);
    //    this.searchId = searchId;
    //    var uriStack = [
    //      resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
    //      searchId];
    //    appModel.setUriStack(uriStack);
    //  }
    //},
    
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
        appModel.setUriStack(uriStack);
      }else{
        // must change the route, and create a post
        var parsedSearches = text_to_search.split('')
        var searchId = ( new Date() ).getTime();
        appModel.setSearch(searchId,text_to_search);
        this.searchId = searchId;
        var resource = appModel.getResource('compound_search');
        var uriStack = [
          resource.key, appModel.URI_PATH_COMPLEX_SEARCH, 
          searchId];
        appModel.setUriStack(uriStack);
      }
    },
    
    // Backbone layoutmanager callback
    cleanup: function(){
      this.plateCopyForm.setValue('searchVal', null);
      this.wellsForm.setValue('searchVal', null);
      this.compoundNameForm.setValue('searchVal', null);
      this.cprForm.setValue('searchVal', null);
      this.screeningInquiryForm.setValue('searchVal', null);
      if (this.form1){
        this.form1.setValue('user_id',null);
        this.form1.setValue('screen_facility_id',null);
        this.form1.setValue('library_short_name',null);
      }
    },
    
  });

  return SearchView;
});