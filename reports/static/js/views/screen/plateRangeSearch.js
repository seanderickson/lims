define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'utils/plateRangeTable',
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'templates/genericResource.html'
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, 
         PlateRangePrototype, DetailView, EditView, genericLayout ) {
  
  /**
   * For a given Screen, show current plate ranges, or perform a screening 
   * inquiry to search for plate ranges.
   * 
   * Screening inquiry form data may be set as a "search" param on the URI 
   * for the page to preload the form.
   */
  var PlateRangeSearchView = Backbone.Layout.extend({
    
    template: _.template(genericLayout),

    si_formatter: new Iccbl.SIUnitsFormatter({ symbol: 'L' }),
    DEFAULT_VOLUME: 33e-9,
    MAX_DAYS_FROM_LAST_SCREENING: 30,
    NUMBER_OF_PREVIOUS_PROTOCOLS_TO_CONSIDER: 3,
    ERR_MSG_LAST_SCREENING: 
      'Last visit was > {max_days} days ago ({last_screening_date})',
    ERR_MSG_PROTOCOL_REPLICATES: 
      'Replicate count does not match other visits: ',
    ERR_MSG_PROTOCOL_VOL: 'Volume does not match other visits: ',
    ERR_MSG_PIN_TRANSFER_NOT_APPROVED: 'The pin transfer approval date has not been entered',
    ERR_MSG_DATA_MEETING_NOT_COMPLETED: 'The data meeting completed date has not been entered',
    
    initialize: function(args) {

      var self = this;
      this.uriStack = args.uriStack;
      this.screenSummaryView = args.summaryView;
      this.currentLibraryScreenings = args.currentLibraryScreenings;
      this.consumedStack = [];
      
      this.defaultShowPlatesUri = Iccbl.formatString(
        '#screen/{facility_id}/summary/plates', self.model);
      var urlStackData = this.urlStackData = this.parseUrlStack(this.uriStack);
      var url = this.url = [
        self.model.resource.apiUri, self.model.get('facility_id'),
        'plate_range_search'].join('/')

      // create the search form
      
      var formSchema = {};
      var libraryScreeningResource = appModel.getResource('libraryscreening');
      var volumeField = _.result(
        libraryScreeningResource['fields'], 
        'volume_transferred_per_well_from_library_plates', {});
      var TextArea2 = Backbone.Form.editors.TextArea.extend({
        render: function() {
          TextArea2.__super__.render.apply(this,arguments);
          this.$el.attr('placeholder', 
            'Enter Plate Number ranges, followed by the Copy Name');
          return this;
        },        
      });
      function validatePlateSearch(value, formValues){
        var errors = [];
        var final_search_array = Iccbl.parseRawPlateSearch(value,errors);
        if (_.isEmpty(final_search_array)){
          errors.push('no values found for input');
        } else {
          console.log('final_search_array', final_search_array);
        }
        _.each(final_search_array, function(search_line){
          if (_.isEmpty(search_line.plates)
              &&_.isEmpty(search_line.plate_ranges)){
            errors.push('must specify a plate or plate-range: ' 
              + search_line.combined.join(', ') );
          }
        });
        if (!_.isEmpty(errors)){
          return {
            type: 'plate_search',
            message: errors.join('; ')
          };
        }
      };
      formSchema['plate_search'] = {
        title: 'Search for Plate Ranges',
        key: 'plate_search',
        editorClass: 'input-full form-control',
        validators: ['required', validatePlateSearch],
        type: TextArea2,
        template: appModel._field_template
      };
      formSchema['volume_required'] = {
        title: 'Volume Required',
        key: 'volume_required',
        validators: ['required',EditView.CheckPositiveNonZeroValidator],
        type: EditView.SIunitEditor,
        template: appModel._field_template
      };
      _.extend(
        formSchema['volume_required'],
        volumeField['display_options']);
      formSchema['replicate_count'] = {
        title: 'Replicates',
        key: 'replicate_count',
        editorClass: 'form-control',
        validators: ['required',EditView.CheckPositiveNonZeroValidator],
        type: Backbone.Form.editors.Number,
        template: appModel._field_template
      };
      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs) {
          var errs = {};
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      if (!_.isEmpty(urlStackData.plate_search)){
        formFields.set('plate_search', urlStackData.plate_search.join('\n'));
      }
      if (_.isNumber(parseFloat(urlStackData.volume_required))){
        formFields.set('volume_required', urlStackData.volume_required);
      } else {
        var lastLibraryScreening = self.getLastLibraryScreening();
        var volume_required = self.DEFAULT_VOLUME;
        if (lastLibraryScreening){
          volume_required = _.result(
            lastLibraryScreening.attributes,
            'volume_transferred_per_well_to_assay_plates',
            volume_required);
        }
        formFields.set('volume_required', volume_required);
      }
      if (_.isNumber(urlStackData.replicate_count)){
        formFields.set('replicate_count', urlStackData.replicate_count);
      } else {
        var replicate_count = 1;
        var lastLibraryScreening = self.getLastLibraryScreening();
        if (lastLibraryScreening){
          replicate_count = _.result(
            lastLibraryScreening.attributes,
            'number_of_replicates',
            replicate_count);
        }
        formFields.set('replicate_count', replicate_count);
      }
      
      this.form = new Backbone.Form({
        model: formFields,
        template: _.template([
          '<div>',
          '<form class="form-horizontal container" >',
          '<div data-fields="plate_search"></div>',
          '<div class="row">', 
          '<div class="col-sm-4">', 
          '<div data-fields="volume_required"></div>',
          '</div>',
          '<div class="col-sm-4">', 
          '<div data-fields="replicate_count"></div>',
          '</div>',
          '<div class="col-sm-4">', 
          '<label class="control-label " for="total_volume_required">Total Volume</label>',
          '<div id="total_volume_required"></div>',
          '</div>',
          '</div>',
          '<div id="plate_search_error" class="text-danger"></div>',
          '</form>',
          '</div>'
          ].join(''))
      });
      var showRetiredPlatesControl = this.showRetiredPlatesControl = $([
          '<label class="checkbox-inline" ', 
          ' style="margin-left: 10px;" ',
          'title="Show all \'available\' and \'retired\' plates" >',
          '  <input id="show_retired_plates_control" ',
          '    type="checkbox">Show available and retired plates',
          '</label>'
        ].join(''));
      showRetiredPlatesControl.find('input[type="checkbox"]')
        .prop('checked', urlStackData.show_retired_plates);
      var showFirstCopyOnly = this.showFirstCopyOnly = $([
          '<label class="checkbox-inline" ', 
          ' style="margin-left: 10px;" ',
          'title="Show the first active copy ',
          '(status is &quot;available&quot, located in freezer 1, ',
          'alphabetical sort, by copy name)" >',
          '  <input id="show_first_copy_only" ',
          '    type="checkbox">Show first active copy only',
          '</label>'
        ].join(''));
      showFirstCopyOnly.find('input[type="checkbox"]')
        .prop('checked', !urlStackData.show_retired_plates);
      var showExistingControl = this.showExistingControl = $([
          '<label class="checkbox-inline" ', 
          ' style="margin-left: 10px;" ',
          'title="Show plate ranges from existing screenings for the screen" >',
          '  <input id="show_existing" ',
          '    type="checkbox">Show existing plates',
          '</label>'
        ].join(''));
      showExistingControl.find('input[type="checkbox"]')
        .prop('checked', urlStackData.show_existing);

      var downloadButton = this.downloadButton = $([
        '<button type="button" class="btn btn-default btn-xs pull-right" ',
        'id="download_button" >download</button>',
      ].join(''));
      var showPlatesLink = this.showPlatesLink = $('<a>', {
        tabIndex : -1,
        href : self.defaultShowPlatesUri,
        target : '_blank',
        title: 'Display plates for the current page'
      }).text('show plates');
      showPlatesLink.addClass('btn btn-default btn-xs pull-right');

      // Set up the collection
      
      var Collection = Iccbl.CollectionOnClient.extend({
        // explicitly define the id so that collection compare & equals work
        lsPLATE_RANGE_SPECIFIER: '{library_screening_id}:' + Iccbl.PLATE_RANGE_KEY_SPECIFIER,
        url: self.url,
        modelId: function(attrs) {
          return Iccbl.formatString(this.lsPLATE_RANGE_SPECIFIER, attrs);
        }
      });
      var plate_range_collection = this.plate_range_collection = new Collection({});
      
      // FIXME: use this as the default sort to allow column sorting
      plate_range_collection.comparator = function(one,two){
        var ls1 = one.get('library_screening_id');
        var ls2 = two.get('library_screening_id');
        if (ls1!=ls2){
          return ls1 < ls2 ? -1 : 1;
        }
        var start1 = one.get('start_plate');
        var start2 = two.get('start_plate');
        if (start1!=start2){
          return start1 < start2 ? -1 : 1;
        }
        var copy1 = one.get('copy_name');
        var copy2 = two.get('copy_name');
        
        return copy1 < copy2 ? -1: 1;
        
      };
      this.listenTo(plate_range_collection,'add', function(model){
        model.collection = plate_range_collection; // backbone bug: this should not be necessary
      });
      
    },

    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },

    /**
     * Set the form values to the URL using:
     * Plate Search URI mini-language:
     * [rest_of_url]/search/[plate_search_mini_language]
     * where plate_search_mini_language is:
     * [plate_range_specifier 1];[plate_range_specifier n];[volume_specifier];[show_retired_plates]
     * where:
     * plate_range_specifier is:
     * [copy_name 1],[copy_name n],[plate 1],[plate n],[plate_range 1],[plate_range n]
     * volume_specifier is:
     * [fixed point number]|[exponential notation number]
     * "show_retired_plates is a literal flag
     */
    reportFormUri: function(search_value,volume_required,replicate_count, 
        show_retired_plates,show_existing){
      var self = this;
      var errors = [];
      
      var formData = {};
      
      formData.plate_ranges = Iccbl.parseRawPlateSearch(search_value,errors);
      formData.volume_required = volume_required;
      formData.replicate_count = replicate_count;
      
      var url_search_data = self.encodeFormData(formData);
      
      if (show_retired_plates=='true'){
        url_search_data.push('show_retired_plates');
      }
      if (show_existing == 'true'){
        url_search_data.push('show_existing');
      }
      var newStack = ['search', url_search_data.join(appModel.SEARCH_DELIMITER)];
      console.log('newStack:', newStack);
      self.reportUriStack(newStack);
    },
    
    /**
     *  plateSearchData { 
     *    screen_facility_id, volume_required, 
     *    plate_ranges: [[combined_plate-range-copy array]], 
     *    replicate_count 
     *  }
     */
    encodeFormData: function(formData){
      var self = this;
      var url_search_data = []

      _.each(formData.plate_ranges, function(plate_search_data){
        // encodeURIComponent to handle spaces and quoting
        var plate_searches = _.map(plate_search_data.combined, encodeURIComponent);
        url_search_data.push(plate_searches.join(','));
      });
      if (!_.isUndefined(formData.volume_required)){
        
        var si_formatter = new Iccbl.SIUnitsFormatter({ symbol: 'L' });
        var volText = si_formatter.fromRaw(formData.volume_required);
        // remove space
        volText = volText.replace(/\s+/,'');
        if (!_.isUndefined(formData.replicate_count)){
          volText = formData.replicate_count + 'x' + volText;
        }
        url_search_data.push(volText);
      }
      return url_search_data;
    },
    
    /**
     * Decode the form values from the URL "search" parameter, encoded using 
     * the above described "plate_search_mini_language"
     * 
     * @return urlStackData {
     *  plate_search, volume_required, replicate_count, show_retired, show_existing
     * }
     */
    parseUrlStack: function(delegateStack) {
      var self = this;
      console.log('parseUrlStack', delegateStack);

      var defaultUrlStackData = {
        show_retired_plates: false,
        show_existing: false
      };
      if (!_.isEmpty(delegateStack) 
          && delegateStack[0] == 'search'){
        
        var search_data = delegateStack[1];
        var errors = [];
        urlStackData = Iccbl.parseScreeningInquiryURLParam(search_data, errors);
        if (!_.isEmpty(errors)){
          appModel.error(
            'Error parsing search data on the URL: ' + errors.join(','));
        }
        urlStackData = _.extend({}, defaultUrlStackData, urlStackData);
        console.log('urlStackData', urlStackData);
      } else {
        urlStackData = defaultUrlStackData;
      }
      return urlStackData;
    },
    
    setPlatesLink: function() {
      var self = this;
      if (!self.plate_range_collection.isEmpty()){

        var resource = appModel.getResource('librarycopyplate');
        // convert the plate range table into explicit plate searches
        var plate_searches = [];
        self.plate_range_collection.each(function(plate_range_model){
          plate_searches.push(Iccbl.formatString(
            '{start_plate}-{end_plate},{copy_name}', plate_range_model));
        });
        if (_.isEmpty(plate_searches)){
          console.log('nothing to search! ', self.plate_range_collection);
        }
        var encodedSearch = encodeURIComponent(plate_searches.join(';'));
        var route = ['#'+resource.key, 'search', 'raw_search_data=' + encodedSearch];
        self.showPlatesLink.attr("href", route.join('/'));
      } else {
        self.showPlatesLink.attr("href",self.defaultShowPlatesUri);
      }
    },
    
    formDownload: function(){
      var self = this;
      var form = self.form;
      var showFirstCopyOnly = self.showFirstCopyOnly;
      var showExistingControl = self.showExistingControl;
      var showRetiredPlatesControl = self.showRetiredPlatesControl;

      // Note: the form does not have to be filled in for the download
      // var errors = form.commit({ validate: true }); 
      var post_data = {};
      post_data['volume_required'] = form.getValue('volume_required');
      post_data['raw_search_data'] = form.getValue('plate_search');
      var show_retired_plates = null;
      if (showRetiredPlatesControl.find('input[type=checkbox]').prop('checked')){
        show_retired_plates = 'true';
        post_data['show_retired_plates'] =show_retired_plates;
      }
      if (showFirstCopyOnly.find('input[type=checkbox]').prop('checked')){
        post_data['show_first_copy_only'] = 'true';
      }
      if (!showExistingControl.find('input[type=checkbox]').prop('checked')){
        if (!_.isEmpty(form.getValue('plate_search'))){
          // Only send the hide_existing flag if the form is filled
          post_data['hide_existing'] = 'true';
        }
      }
      
      var downloadFormSchema = {};
      downloadFormSchema['content_type'] = {
        title: 'Download type',
        help: 'Select the data format',
        key: 'content_type',
        options: _.without(self.model.resource.content_types, 'json'), // never json
        type: 'Select',
        validators: ['required'],
        template: _.template([
          '<div class="input-group col-xs-6">',
          '   <label class="input-group-addon" for="<%= editorId %>" ',
          '         title="<%= help %>" ><%= title %></label>',
          '   <span  data-editor></span>',
          '</div>',
        ].join('')),
        editorClass: 'form-control'
      };
      var dFormFields = Backbone.Model.extend({
        schema: downloadFormSchema
      });
      var dformFields = new dFormFields();
      dformFields.set('content_type', downloadFormSchema['content_type'].options[0]);
      var downloadForm = new Backbone.Form({
        model: dformFields,
        template: _.template([
          "<div>",
          "<form data-fieldsets class='form-horizontal container' >",
          "</form>",
          ].join(''))
      });
      var el = downloadForm.render().el;
      appModel.showModal({
        view: el,
        title: 'Download',  
        ok: function(event){
          var url = self.url + '?format=' + downloadForm.getValue('content_type');
          appModel.downloadUrl(url, post_data);
        }
      });
    },
    
    calculateTotalVolume: function(){
      var self = this;
      var form = self.form;
      var totalVolume = final_volume = form.getValue('volume_required');
      var replicate_count = form.getValue('replicate_count');
      if (replicate_count > 0){
        totalVolume*= replicate_count;
      }
      var volText = self.si_formatter.fromRaw(totalVolume);
      form.$el.find('#total_volume_required').html(volText);
      return totalVolume;
    },

    getLastLibraryScreening: function(){
      if (!this.currentLibraryScreenings.isEmpty()){
        this.currentLibraryScreenings.sortBy('-date_of_activity');
        return this.currentLibraryScreenings.findWhere({
          'is_for_external_library_plates': false});
      }
      return null;
    },
    
    checkLastScreeningDate: function(){
      if (!this.currentLibraryScreenings.isEmpty()){
        this.currentLibraryScreenings.sortBy('-date_of_activity');
        var lastScreening = this.currentLibraryScreenings.at(0);
        var lastDate = new Date(lastScreening.get('date_of_activity'));
        var currentDate = new Date();
        console.log('checkLastScreeningDate: ', lastDate, currentDate, 
          currentDate-lastDate, (currentDate-lastDate)/(24*60*60*1000));
        var max_interval_ms = this.MAX_DAYS_FROM_LAST_SCREENING*24*60*60*1000;
        if(currentDate-lastDate > max_interval_ms ){
          appModel.error(Iccbl.formatString(
            this.ERR_MSG_LAST_SCREENING, { 
              max_days: this.MAX_DAYS_FROM_LAST_SCREENING,
              last_screening_date: Iccbl.getDateString(lastDate) 
            }));
        }
      }
    },
    
    /**
     * Check the entered protocol against existing screenings and display an 
     * error message if there are differences.
     */
    checkProtocol: function(volume_per_assay_plate, replicate_count) {
      var self = this;
      if (!this.currentLibraryScreenings.isEmpty()){
        var errVol = {}
        var errReplicates = {};
        var screeningsForProtocolCheck = _.reject(
            this.currentLibraryScreenings.sortBy('-date_of_activity'),
            function(model){
              return model.get('is_for_external_library_plates');
            });
        
        if (_.isEmpty(screeningsForProtocolCheck)){
          return;
        }
        screeningsForProtocolCheck = screeningsForProtocolCheck.slice(
          0,this.NUMBER_OF_PREVIOUS_PROTOCOLS_TO_CONSIDER);
        
        console.log('screeningsForProtocolCheck', screeningsForProtocolCheck);
        _.each(screeningsForProtocolCheck,function(model){
          
          if (model.get('is_for_external_library_plates')){
            return;
          }
          
          var otherVol = model.get('volume_transferred_per_well_to_assay_plates');
          otherVol = parseFloat(otherVol);
          var otherReplicates = model.get('number_of_replicates');
          
          if (otherVol !== volume_per_assay_plate){
            console.log('other protocol vol', otherVol, volume_per_assay_plate);
            var otherVolText = self.si_formatter.fromRaw(otherVol);
            var listing = _.result(errVol,otherVolText,[]);
            listing.push(model.get('activity_id'));
            errVol[otherVolText] = listing;
          }
          if(otherReplicates != replicate_count) {
            console.log('other protocol replicates', otherReplicates, replicate_count);
            var listing = _.result(errReplicates,otherReplicates,[]);
             listing.push(model.get('activity_id'));
             errReplicates[otherReplicates] = listing;
          }
        });
        
        if (!_.isEmpty(errVol)){
          appModel.error(self.ERR_MSG_PROTOCOL_VOL
            + _.map(_.pairs(errVol),function(pair){
              console.log('pair', pair);
              return '(' + pair[0] + '): ' + pair[1].join(', ');
            }).join('; '));
        }
        if (!_.isEmpty(errReplicates)){
          appModel.error(self.ERR_MSG_PROTOCOL_REPLICATES 
            + _.map(_.pairs(errReplicates),function(pair){
              return '( count: ' + pair[0] + '): ' + pair[1].join(', ');
            }).join('; '));
        }
      }
      
      if (_.isEmpty(self.model.get('pin_transfer_date_approved'))){
        appModel.error(self.ERR_MSG_PIN_TRANSFER_NOT_APPROVED);
      }
      if (_.isEmpty(self.model.get('data_meeting_complete'))){
        appModel.error(self.ERR_MSG_DATA_MEETING_NOT_COMPLETED);
      }
    }, 
    
    formSubmit: function() {
      var self = this;
      var form = self.form;
      var showFirstCopyOnly = self.showFirstCopyOnly;
      var showRetiredPlatesControl = self.showRetiredPlatesControl;
      var showExistingControl = self.showExistingControl;
      appModel.clearErrors();
      var errors = form.commit({ validate: true }); 
      if(!_.isEmpty(errors)){
        if (_.contains(errors,'plate_search')){
          form.$el.find('#plate_search').addClass("has-error");
          form.$el.find('#plate_search_error').html(errors['plate_search'].message);
        }
        if (_.contains(errors,'volume_required')){
          form.$el.find('#volume_required').addClass("has-error");
          form.$el.find('#volume_required_error').html(errors['volume_required'].message);
        }
        return;
      }else{
        form.$el.find('#plate_search_error').empty();
        form.$el.find('#plate_search').removeClass("has-error");
      }
      var data = new FormData();
      var search_value = form.getValue('plate_search');
      data.append('raw_search_data', search_value);
      
      // NOTE/TODO: 20170504; 
      // The server API will now support the plate_search mini-language, 
      // which may be sent as an encoded URI param: e.g. 
      // raw_search_data=1814-1821,J&raw_search_data=1709,%22B%201%3A125%22
      // is a valid "raw_search_data" array param for plate searches
      // returned by:
      // Iccbl.parseRawPlateSearch(search_value); using 
      // encodeURIComponent for each line;
      // User data may also be posted, as-is as form data 
      // (application/x-www-form-urlencoded), which is done here, 
      // for the server to parse.
      
      var volume_required = form.getValue('volume_required');
      var replicate_count = form.getValue('replicate_count');
      var show_existing = null;
      
      data.append('volume_required', self.calculateTotalVolume());
      var show_retired_plates = null;
      if (showRetiredPlatesControl.find('input[type=checkbox]').prop('checked')){
        show_retired_plates = 'true';
        data.append('show_retired_plates', 'true');
      }
      if (showExistingControl.find('input[type=checkbox]').prop('checked')){
        show_existing = 'true';
      } else {
        data.append('hide_existing', 'true');
      }
      if (showFirstCopyOnly.find('input[type=checkbox]').prop('checked')){
        data.append('show_first_copy_only', 'true');
      }
      self.reportFormUri(search_value, volume_required, replicate_count, 
        show_retired_plates, show_existing);
      
      self.checkProtocol(volume_required, replicate_count);
      self.checkLastScreeningDate();
      
      $.ajax({
        url:  self.url,     
        data: data,
        cache: false,
        contentType: false,
        processData: false,
        dataType: 'json', // what is expected back from the server
        type: 'POST',
      }).done(function(data, textStatus, jqXHR){
        console.log('success', data);
        var objects = _.result(data, 'objects');
        if (!_.isEmpty(objects)){
          self.plate_range_collection.set(objects);
        } else {
          self.plate_range_collection.reset(null);
        }
        var meta = _.result(data,appModel.API_RESULT_META, {});
        var warning = _.result(meta, appModel.API_META_MSG_WARNING);
        if (warning){
          appModel.error(warning);
        }
        self.setPlatesLink();
      }).fail(function(jqXHR, textStatus, errorThrown) { 
        console.log('errors', arguments);
        appModel.jqXHRfail.apply(this,arguments); 
      });
    },
    
    afterRender: function(){
      var self = this;
      var form = self.form;
      var showFirstCopyOnly = self.showFirstCopyOnly;
      var showRetiredPlatesControl = self.showRetiredPlatesControl;
      var showExistingControl = self.showExistingControl;
      var urlStackData = self.urlStackData;
           
      var $form_toggle = $('<a>search >></a>');
      var form_shown = false;
      var $form_div = $('<div id="plate_search_form_div" style="display: none;" ></div>');
      $('#resource_content').append($form_toggle);
      $('#resource_content').append(self.showPlatesLink);
      $('#resource_content').append(self.downloadButton);
      $('#resource_content').append($form_div);
      function toggleForm(){
        if (form_shown){
          form_shown = false;
          $form_div.hide();
          $form_toggle.html('search >>');
        } else {
          form_shown = true;
          $form_div.show();
          $form_toggle.html('<< search');
        }
      };
      $form_toggle.click(function(e){
        e.preventDefault();
        toggleForm();
      });
      
      $form_div.append(form.render().el);
      form.$el.append([
        '<button type="button" class="btn btn-default btn-xs" ',
        'id="submit_button" >search</button>',
      ].join(''));
      $('#resource_content').append('<div id="plate_range_table"></div>');

      var helpLink = $('<a >&nbsp;?</a>');
      form.$el.find('[key="form-group-plate_search"]').find('label').append(helpLink);
      helpLink.click(function(e){
        e.preventDefault();
        var bodyMessage = [
          'Plate Range search, e.g.:',
          '1000',
          '1000-2000 A',
          'B 3000-4000',
          '5000-6000 A,B,C     ',
          '9000,2000,3000 5000-4000 A,"b-c",D'
        ];
        appModel.showModalMessage({
          title: 'Search for plates using patterns',
          body: bodyMessage.join('<br/>')
        });
      });
      form.$el.find('form.form-horizontal').append(showRetiredPlatesControl);
      form.$el.find('form.form-horizontal').append(showFirstCopyOnly);
      form.$el.find('form.form-horizontal').append(showExistingControl);
      
      PlateRangePrototype._createPlateRangeTable.call(this,
        self.plate_range_collection, $('#plate_range_table'), false, 
        ['library_screening_id', 'plate_locations', 'warnings','errors'],
        self.model.get('facility_id'));
      
      showRetiredPlatesControl.find('input[type=checkbox]').change(function(e) {
        if (showFirstCopyOnly.find('input[type=checkbox]').prop('checked')){
          showFirstCopyOnly.find('input[type="checkbox"]')
            .prop('checked',false);
        }
        self.formSubmit();
      });
      showFirstCopyOnly.find('input[type=checkbox]').change(function(e) {
        if (showRetiredPlatesControl.find('input[type=checkbox]').prop('checked')){
          showRetiredPlatesControl.find('input[type="checkbox"]')
            .prop('checked',false);
        }
        self.formSubmit();
      });
      showExistingControl.find('input[type=checkbox]').change(function(e) {
        self.formSubmit();
      });
      form.$el.find('#submit_button').click(function(e){
        e.preventDefault();
        self.formSubmit();
      });
      if (!_.isEmpty(urlStackData.plate_search) 
          || !_.isUndefined(urlStackData.volume_required)){
        toggleForm();
        if (!_.isEmpty(urlStackData.plate_search) 
            && !_.isUndefined(urlStackData.volume_required)){
          self.formSubmit();
        }
      } else {
        self.plate_range_collection.fetch();
      }
      self.downloadButton.click(function(e){
        e.preventDefault();
        self.formDownload();
      });
      
      self.showPlatesLink.click(function(e){
        if (self.defaultShowPlatesUri == self.showPlatesLink.attr('href')){
          e.preventDefault();
          self.screenSummaryView.change_to_tab('plates');          
        }
      });
      this.listenTo(
        form, "replicate_count:change", self.calculateTotalVolume);
      this.listenTo(
        form, "volume_required:change", 
        self.calculateTotalVolume);
      self.calculateTotalVolume();
    }
  });

  return PlateRangeSearchView;
});