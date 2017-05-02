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
   * Show a table of the plate ranges for the plates in the library 
   * screenings for the Screen.
   * Display a plate range search form above the table:
   * - plates for the search are concatenated with the existing table
   * - form data set as a "search" param on the URI for the page.
   * 
   */
  var PlateRangeSearchView = Backbone.Layout.extend({
    
    template: _.template(genericLayout),
    
    initialize: function(args) {

      var self = this;
      this.uriStack = args.uriStack;
      this.screenSummaryView = args.summaryView;
      this.consumedStack = [];

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
        var final_search_array = Iccbl.parsePlateSearch(value,errors);
        if (_.isEmpty(final_search_array)){
          errors.push('no values found for input');
        } else {
          console.log('final_search_array', final_search_array);
        }
        _.each(final_search_array, function(search_line){
          if (_.isEmpty(search_line.plates)
              &&_.isEmpty(search_line.plate_ranges)){
            errors.push('must specify a plate or plate-range: ' 
              + search_line.original_text );
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
      if (urlStackData.volume_required){
        formFields.set('volume_required', urlStackData.volume_required);
      }
      this.form = new Backbone.Form({
        model: formFields,
        template: _.template([
          '<div>',
          '<form class="form-horizontal container" >',
          '<div data-fields="plate_search"></div>',
          '<div data-fields="volume_required"></div>',
          '<div id="plate_search_error" class="text-danger"></div>',
          '</form>',
          '</div>'
          ].join(''))
      });
      var showRetiredPlatesControl = this.showRetiredPlatesControl = $([
          '<label class="checkbox-inline" ', 
          ' style="margin-left: 10px;" ',
          'title="Show \'retired\' plates" >',
          '  <input id="show_retired_plates_control" ',
          '    type="checkbox">Show retired plates',
          '</label>'
        ].join(''));
      showRetiredPlatesControl.find('input[type="checkbox"]')
        .prop('checked', urlStackData.show_retired_plates);
      var showFirstCopyOnly = this.showFirstCopyOnly = $([
          '<label class="checkbox-inline" ', 
          ' style="margin-left: 10px;" ',
          'title="Show the first copy (alphabetical, by copy name)" >',
          '  <input id="show_first_copy_only" ',
          '    type="checkbox">Show first matching copy only',
          '</label>'
        ].join(''));
      showFirstCopyOnly.find('input[type="checkbox"]')
        .prop('checked', !urlStackData.show_retired_plates);

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
      plate_range_collection.fetch();
      
      // FIXME: use this as the default sort to allow column sorting
      plate_range_collection.comparator = function(one,two){
        var ls1 = one.get('library_screening_id');
        var ls2 = two.get('library_screening_id');
        if (ls1!=ls2){
          return ls1 < ls2 ? -1 : 1;
        }
        var start1 = one.get('start_plate');
        var start2 = two.get('start_plate');
        return start1 < start2 ? -1 : start1 == start2 ? 0 : 1;
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

    reportFormUri: function(search_value,volume_required,show_retired_plates){
      var self = this;
      var errors = [];
      var url_search_data = []
      var plate_search_data_array = Iccbl.parsePlateSearch(search_value,errors);
      console.log('plate_search: ', plate_search_data_array);
      _.each(plate_search_data_array, function(plate_search_data){
        var plate_searches = [];
        _.each(plate_search_data.plates, function(plate){
          plate_searches.push(plate);
        });
        _.each(plate_search_data.plate_ranges, function(plate_range){
          plate_searches.push(plate_range);
        })
        var copy_searches = [];
        _.each(plate_search_data.copies,function(copy){
          copy_searches.push(encodeURIComponent(copy));
        });
        search_line = plate_searches.join(',');
        if (!_.isEmpty(copy_searches)){
          search_line = copy_searches.join(',') + ',' + search_line;
        }
        console.log('search line: ', search_line);
        url_search_data.push(search_line);
      });
      if (volume_required){
        url_search_data.push(volume_required);
      }
      if (show_retired_plates=='true'){
        url_search_data.push('show_retired_plates');
      }
      var newStack = ['search', url_search_data.join(';')];
      console.log('newStack:', newStack);
      self.reportUriStack(newStack);
    },

    parseUrlStack: function(delegateStack) {
      
      console.log('parseUrlStack1', delegateStack);
      var errors = [];
      var urlStackData = {
        show_retired_plates: false
      };
      if (!_.isEmpty(delegateStack) 
          && delegateStack[0] == 'search'){
        
        var url_plate_searches = [];
        var search_data = delegateStack[1];
        var extra_volumes = [];
        console.log('parseUrlStack: ', search_data);
        _.each(search_data.split(';'), function(element){
          console.log('parse element', element);
          if (element == 'show_retired_plates'){
            urlStackData.show_retired_plates = true;
            return;
          }
          // Detect decimal (group 1) or number with exponent (group 2) as a volume
          else if (/^((\d+\.?\d?e\-\d+)?|(\d?\.\d+)?)$/.exec(element)){
            if (urlStackData.volume_required){
              console.log('more than one volume has been defined, ignoring extra');
              extra_volumes.push(urlStackData.volume_required);
              extra_volumes.push(element);
            }else{
              urlStackData.volume_required = element;
            }
            return;
          }
          
          var copies = [];
          var plate_or_plate_ranges = [];
          var parts = element.split(',');
          _.each(parts, function(part){
            part = part.trim();
            if (_.isEmpty(part)) return;
            
            if (part.match(Iccbl.PLATE_PATTERN)){
              plate_or_plate_ranges.push(part);
            } else if (part.match(Iccbl.PLATE_RANGE_PATTERN)){
              plate_or_plate_ranges.push(part);
            } else {
              part = decodeURIComponent(part);
              if (part.match(Iccbl.COPY_NAME_PATTERN)) {
                if (part.match(/[ \,\-]/)){
                  part = '"' + part + '"';
                }
                copies.push(part);
              } else {
                console.log('url part not recognized:', part);
                errors.push('Not recognized: "' + part + '" in "' + element + '"');
                return;
              }
            }
          });
          if (_.isEmpty(plate_or_plate_ranges)){
            console.log('error, could not detect plates in search line: ', element);
            errors.push('Plate or plate range not detected: "'+ element + '"');
            return;
          }
          url_search = plate_or_plate_ranges.join(',');
          if (!_.isEmpty(copies)){
            url_search = copies.join(',') + ' ' + url_search;
          }
          url_plate_searches.push(url_search);
        });
        urlStackData['plate_search'] = url_plate_searches;
        
        if (!_.isEmpty(extra_volumes)){
          appModel.error(
            'More than one volume required specified in the URL; ' +
            'Extra specifiers are ignored: ' + extra_volumes.join(', '));
        }
        if (!_.isEmpty(errors)){
          appModel.error(
            'Error parsing earch data on the URL: ' + errors.join(','));
        }
        console.log('urlStackData', urlStackData);
      }
      return urlStackData;
    },
    
    /**
     * Show the LibraryCopyPlates for the current plate ranges
     */
    showPlates: function() {
      var self = this;
      var form = self.form;
      var showFirstCopyOnly = self.showFirstCopyOnly;
      var showRetiredPlatesControl = self.showRetiredPlatesControl;
      // Note: the form does not have to be filled in for the plates dialog
      // var errors = form.commit({ validate: true }); 
      var plate_search = form.getValue('plate_search');
      
      if (!_.isEmpty(plate_search) && !self.plate_range_collection.isEmpty()){

        var resource = appModel.getResource('librarycopyplate');
        // TODO: key the search data using the searchID: 
        // this allows for browser "back" in the session
        // will also need to listen to URIStack changes and grab the data
        // from the search ID
        var searchID = ( new Date() ).getTime();
        
        // convert the plate range table into explicit plate searches
        var plate_searches = [];
        self.plate_range_collection.each(function(plate_range_model){
          plate_searches.push(Iccbl.formatString(
            '{start_plate}-{end_plate} {copy_name}', plate_range_model));
        });
        console.log('plate_searches: ', plate_searches);
        if (_.isEmpty(plate_searches)){
          console.log('nothing to search! ', self.plate_collection);
        }
        appModel.setSearch(searchID,plate_searches.join('\n'));
        this.searchID = searchID;
        appModel.set('routing_options', {replace: false});  
        var _route = resource.key + '/search/'+ searchID;
        appModel.router.navigate(_route, {trigger:true});
        
      } else {
        self.screenSummaryView.change_to_tab('plates');
      }
    },
    
    formDownload: function(){
      var self = this;
      var form = self.form;
      var showFirstCopyOnly = self.showFirstCopyOnly;
      var showRetiredPlatesControl = self.showRetiredPlatesControl;

      // Note: the form does not have to be filled in for the download
      // var errors = form.commit({ validate: true }); 
      var post_data = {};
      post_data['volume_required'] = form.getValue('volume_required');
      post_data['raw_search_data'] = form.getValue('plate_search');
      console.log('post_data', post_data);
      var show_retired_plates = null;
      if (showRetiredPlatesControl.find('input[type=checkbox]').prop('checked')){
        show_retired_plates = 'true';
        post_data['show_retired_plates'] =show_retired_plates;
      }
      if (showFirstCopyOnly.find('input[type=checkbox]').prop('checked')){
        post_data['show_first_copy_only'] = 'true';
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
          appModel.downloadUrl(downloadForm.$el, url, post_data);
        }
      });
    },

    formSubmit: function() {
      var self = this;
      var form = self.form;
      var showFirstCopyOnly = self.showFirstCopyOnly;
      var showRetiredPlatesControl = self.showRetiredPlatesControl;

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
      var volume_required = form.getValue('volume_required');
      data.append('raw_search_data', search_value);
      data.append('volume_required', volume_required);
      var show_retired_plates = null;
      if (showRetiredPlatesControl.find('input[type=checkbox]').prop('checked')){
        show_retired_plates = 'true';
        data.append('show_retired_plates', 'true');
      }
      if (showFirstCopyOnly.find('input[type=checkbox]').prop('checked')){
        data.append('show_first_copy_only', 'true');
      }
      self.reportFormUri(search_value, volume_required, show_retired_plates);
      
      var headers = {}; // can be used to send a comment
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
        }
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
      var urlStackData = self.urlStackData;
      
      var $form_toggle = $('<a>search >></a>');
      var $downloadButton = $([
        '<button type="button" class="btn btn-default btn-xs pull-right" ',
        'id="download_button" >download</button>',
      ].join(''));
      var $showPlatesButton = $([
        '<button type="button" class="btn btn-default btn-xs pull-right" ',
        'id="show_plates_button" >show plates</button>',
      ].join(''));
      
      var form_shown = false;
      var $form_div = $('<div id="plate_search_form_div" style="display: none;" ></div>');
      $('#resource_content').append($form_toggle);
      $('#resource_content').append($showPlatesButton);
      $('#resource_content').append($downloadButton);
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
      
      PlateRangePrototype._createPlateRangeTable.call(this,
        self.plate_range_collection, $('#plate_range_table'), false, 
        ['library_screening_id', 'library_screening_status', 'warnings','errors'],
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
      form.$el.find('#submit_button').click(function(e){
        e.preventDefault();
        self.formSubmit();
      });
      if (!_.isEmpty(urlStackData.plate_search) || !_.isUndefined(urlStackData.volume_required)){
        console.log('toggle form on open', urlStackData);
        toggleForm();
        if (!_.isEmpty(urlStackData.plate_search) && !_.isUndefined(urlStackData.volume_required)){
          self.formSubmit();
        }
      }
      $downloadButton.click(function(e){
        e.preventDefault();
        self.formDownload();
      });
      $showPlatesButton.click(function(e){
        e.preventDefault();
        self.showPlates();
      });
    }
  });

  return PlateRangeSearchView;
});