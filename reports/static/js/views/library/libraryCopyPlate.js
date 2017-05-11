define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_detail_layout',
  'views/generic_detail_stickit', 
  'views/generic_edit',
  'utils/uploadDataForm',
  'templates/genericResource.html'
], 
function($, _, Backbone, Backgrid, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         DetailView, EditView, UploadDataForm, layout) {
  
  var LibraryCopyPlateView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      this._args = args;
      this._classname = 'libraryCopyPlate';
      this.resource = args.resource;
      this.uriStack = args.uriStack;
      this.library = args.library;
      this.copy = args.copy;
      this.consumedStack = [];
      console.log('uriStack', this.uriStack);
      _.bindAll(this, 'showDetail');
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    },
    
    afterRender: function() {
      var self = this;
      var uriStack = this.uriStack;
      var library = this.library;
      var copy = this.copy;
      var resourceId = 'librarycopyplate';
      var resource = appModel.getResource(resourceId);

      if (self.model){
        self.showDetail(self.model);
      }
      else if (!_.isEmpty(uriStack) &&
              !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        // Detail view
        var plate_number = uriStack.shift();
        this.consumedStack = [plate_number];
        _key = library.key + '/' + copy.get('copy_name')+ '/' + plate_number;
        appModel.getModel(resourceId, _key, this.showDetail );
      } else {
        // List view
        var url = resource.apiUri;
        if (self.library && self.copy) {
          url = [ 
              library.resource.apiUri,library.key,'copy',copy.get('copy_name'),
              'plate'].join('/');
        }
        console.log('url: ', url);
        this.consumedStack = [];
        this.showList(resource, url);
      }
    },     
    
    showDetail: function(model) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      
      var detailView = DetailView.extend({
        afterRender: function() {
          DetailView.prototype.afterRender.apply(this,arguments);
          var search_data = {
            key: model.key,
            ref_resource_name: model.resource.key,
            comment__is_blank: false
          };
          appModel.createCommentTable(model,search_data, $('#comment_table'));
        }
      });
      
      // TODO: 20170407 - Plate.is_active flag management
      // - disable if status != 'available'
      // - set to false if status is changed to !available
      // extend EditView and add a change listener/custom validation
      
      var view = new DetailLayout({ 
        model: model, 
        uriStack: uriStack,
        title: 'Plate ' + model.get('plate_number'),
        DetailView: detailView
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#resource_content', view).render();
    },

    showList: function(resource, url) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      resource = appModel.cloneResource(resource);
      
      var extraIncludes = [];
      var commentColumns = ['library_comment_array','comment_array','copy_comments'];
      console.log('uriStack', uriStack);
      var showEditLocationButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="batch_edit_locations_button" href="#">',
          'Edit Locations</a>'
        ].join(''));
      var showEditPlatesButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="batch_edit_plates_button" href="#">',
          'Edit Plates</a>'
        ].join(''));
      var showHistoryButton = $([
      '<a class="btn btn-default btn-sm pull-down" ',
        'role="button" id="showHistoryButton" href="#">',
        'History</a>'
      ].join(''));
      var showUploadButton = $([
      '<a class="btn btn-default btn-sm pull-down" ',
        'role="button" id="patch_resource" href="#">',
        'Upload data</a>'
      ].join(''));
      
      var showCommentsControl = $([
          '<label class="checkbox-inline" ', 
          ' style="margin-left: 10px;" ',
          'title="Show comments" >',
          '  <input id="show_comments" ',
          '    type="checkbox">show comments',
          '</label>'
        ].join(''));

      var extraControls = [];
      if (appModel.hasPermission(resource.key, 'write')){
        extraControls = [
         showEditLocationButton, showEditPlatesButton, showHistoryButton,
         showUploadButton];
      }
      extraControls.push(showCommentsControl);
      
      showHistoryButton.click(function(e) {
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', 'search'];
        var search = {};
        search['ref_resource_name'] = 'librarycopyplate';
        if (self.library && self.copy) {
          search['key__icontains'] = [
            self.library.get('short_name'),
            self.copy.get('copy_name')].join('/');
        }
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        console.log('history route: ' + route);
        appModel.router.navigate(route, {trigger: true});
      });

      // Set the visibility of concentration and volume fields based on copywell status
      if (self.library) {
        var concentration_types = self.library.get('concentration_types');
        if (concentration_types) {
          _.each(_.pick(resource['fields'], 
            ['mg_ml_concentration','min_mg_ml_concentration',
             'max_mg_ml_concentration','molar_concentration',
             'min_molar_concentration','max_molar_concentration'
            ]),
            function(field){
              field['visibility'] = [];
            });
          
          if (_.contains(concentration_types, 'mg_ml')){
            resource.fields['mg_ml_concentration']['visibility'] = ['l','d'];
            if (self.copy && self.copy.get('has_copywell_concentrations')) {
              resource.fields['min_mg_ml_concentration']['visibility'] = ['l','d'];
              resource.fields['max_mg_ml_concentration']['visibility'] = ['l','d'];
            }
          }
          if (_.contains(concentration_types, 'molar')){
            resource.fields['molar_concentration']['visibility'] = ['l','d'];
            if (self.copy && self.copy.get('has_copywell_concentrations')) {
              resource.fields['min_molar_concentration']['visibility'] = ['l','d'];
              resource.fields['max_molar_concentration']['visibility'] = ['l','d'];
            }
          }
        }
        _.each(_.pick(resource['fields'], 
          ['avg_remaining_volume','min_remaining_volume',
           'max_remaining_volume'
          ]), 
          function(field){
            if(self.copy && self.copy.get('has_copywell_volumes')){
              field['visibility'] = ['l','d'];
            } else {
              field['visibility'] = [];
            }
        });
        
      }  
      if (!_.isEmpty(self.library) || !_.isEmpty(self.copy)){
        extraIncludes = _.union(
          extraIncludes,commentColumns);
        showCommentsControl.find('input[type="checkbox"]')
          .prop('checked', true);
      }
      resource.fields['library_short_name']['backgridCellType'] =
        Iccbl.CommentArrayLinkCell.extend({
          comment_attribute: 'library_comment_array',
          title_function: function(model){
            return 'Comments for library: ' + model.get('library_short_name');
          }
        });
      resource.fields['copy_name']['backgridCellType'] =
        Iccbl.LinkCell.extend({
          'hrefTemplate': '#library/{library_short_name}/copy/{copy_name}',
          render: function(){
            var self = this;
            Iccbl.LinkCell.prototype.render.apply(this, arguments);
            var comments = this.model.get('copy_comments');
            if (!_.isEmpty(comments)){
              this.$el.attr('title', comments);
              this.$el.append(Iccbl.createCommentIcon(
                comments,
                'Comments for Copy: ' 
                  + self.model.get('library_short_name') + '/'
                  + self.model.get('source_copy_name')
                ));
            }
            return this;
          }
        });
      
      resource.fields['plate_number'].backgridCellType = 
        Iccbl.CommentArrayLinkCell.extend(_.extend({},
          resource.fields['plate_number'].display_options,
          {
            comment_attribute: 'comment_array',
            linkCallback: function(e){
              e.preventDefault();
              // re-fetch the full model
              plate_number = this.model.get('plate_number');
              this.consumedStack = [plate_number];
              _key = [this.model.get('library_short_name'),
                      this.model.get('copy_name'),
                      plate_number].join('/');
              appModel.getModel(resource.key, _key, self.showDetail );
              return;
            }
          }));
      
      // TODO: extending the passed args to get the "search_data", or other 
      // passed args, grab options explicitly instead 201608
      options = _.extend({}, self._args ,{
        uriStack: uriStack,
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: extraControls,
        extraIncludes: extraIncludes
        }) ;
      var view = new ListView(options);
      showEditLocationButton.click(function(e) {
        e.preventDefault();
        self.batchEditLocationsDialog(view);
      });
      showEditPlatesButton.click(function(e) {
        e.preventDefault();
        self.batchEditPlatesDialog(view);
      });
      showUploadButton.click(function(e){
        
        var collection = view.collection;
        e.preventDefault();
        var url = resource.apiUri;
        UploadDataForm.postUploadFileDialog(
          url, resource.content_types)
          .done(function(){
            collection.fetch({ reset: true });
          })
          .fail(function(){
            appModel.jqXHRfail.apply(this,arguments); 
          });
      });
      showCommentsControl.find('input[type=checkbox]').change(function(e) {
        e.preventDefault();
        if (showCommentsControl.find('input[type=checkbox]').prop('checked')){
          view.collection.extraIncludes = _.union(extraIncludes,commentColumns);
          view.collection.fetch();
        } else {
          view.collection.extraIncludes = [];
          view.collection.fetch();
        }
      });
      
      self.listenTo(view, 'update_title', function(val) {
        if(val) {
          this.$('#content_title').html(resource.title + ': <small>' + val + '</small>');
          $('#content_title_row').show();
        } else {
          this.$('#content_title').html(resource.title);
          $('#content_title_row').show();
        }
      });
        
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#resource_content', view ).render();
    },
    
    batchEditPlatesDialog: function(listView) {
      var self = this;
      console.log('batch edit plates dialog...');
      // find any extra search data
      post_data = {};
      if(_.has(listView._options, 'search_data')) {
        console.log('search data found on the list._options: ', listView._options.search_data);
        // endcode for the post_data arg; post data elements are sent 
        // as urlencoded values via a POST form for simplicity
        post_data['search_data'] = JSON.stringify(listView._options.search_data);  
      }
      // Construct the form
      var initfun = function() {
        var formSchema = {};
        var fieldTemplate = appModel._field_template;
        var formTemplate = appModel._form_template;
        
        formSchema['status'] = {
          title: 'Status',
          key: 'status',
          type: EditView.ChosenSelect,
          editorClass: 'form-control chosen-select',
          options: appModel.getVocabularySelectOptions('plate.status'),
          template: fieldTemplate 
        };
        
        formSchema['plate_type'] = {
          title: 'Plate Type',
          key: 'plate_type',
          type: EditView.ChosenSelect,
          editorClass: 'form-control chosen-select',
          options: appModel.getVocabularySelectOptions('plate.type'),
          template: fieldTemplate 
        };
        
        formSchema['remaining_well_volume'] = {
          title: 'Remaining Well Volume',
          key: 'well_volume',
          validators: [EditView.CheckPositiveNonZeroValidator],
          type: EditView.SIunitEditor,
          template: fieldTemplate 
        };
        _.extend(
          formSchema['remaining_well_volume'],
          self.resource['fields']['remaining_well_volume']['display_options']);
        if (self.copy && self.copy.get('has_copywell_volumes')) {
          delete formSchema['remaining_well_volume'];
        }
        
        formSchema['molar_concentration'] = {
          title: 'Molar Concentration',
          key: 'molar_concentration',
          validators: [EditView.CheckPositiveNonZeroValidator],
          type: EditView.SIunitEditor,
          template: fieldTemplate 
        };
        _.extend(
          formSchema['molar_concentration'],
          self.resource['fields']['molar_concentration']['display_options']);
        
        formSchema['mg_ml_concentration'] = {
          title: 'mg/ml Concentration',
          key: 'mg_ml_concentration',
          validators: [EditView.CheckPositiveNonZeroValidator],
          type: Backbone.Form.editors.Number,
          editorClass: 'form-control',
          template: fieldTemplate 
        };
        _.extend(
          formSchema['mg_ml_concentration'],
          self.resource['fields']['mg_ml_concentration']['display_options']);
        if (self.copy && self.copy.get('has_copywell_concentrations')) {
          delete formSchema['molar_concentration'];
          delete formSchema['mg_ml_concentration'];
        }
        
        formSchema['screening_count'] = {
          title: 'Screening Count',
          key: 'screening_count',
          validators: [EditView.CheckPositiveNonZeroValidator],
          type: Backbone.Form.editors.Number,
          editorClass: 'form-control',
          template: fieldTemplate 
        };
        _.extend(
          formSchema['screening_count'],
          self.resource['fields']['screening_count']['display_options']);
        
        formSchema['comments'] = {
          title: 'Comments',
          key: 'comments',
          validators: ['required'],
          type: 'TextArea',
          editorClass: 'input-full',
          template: fieldTemplate
        };
        
        var FormFields = Backbone.Model.extend({
          schema: formSchema,
          validate: function(attrs) {
            var errs = {};
            if (_.isEmpty(attrs)) {
              errs['error'] = 'Must enter values';
            }
            var mgml = attrs['mg_ml_concentration'];
            var molar = attrs['molar_concentration'];
            if(_.isNumber(mgml) && mgml != 0 
                && _.isNumber(molar) && molar != 0){
              var msg = 'Must enter either (mg/ml) or (molar)';
              errs['mg_ml_concentration'] = msg;
              errs['molar_concentration'] = msg;
            }
            if (!_.isEmpty(errs)) return errs;
          }
        });
        var formFields = new FormFields({
          'screening_count': null,
          'remaining_well_volume': null,
          'mg_ml_concentration': null,
          'molar_concentration': null
        });
        
        var form = new Backbone.Form({
          model: formFields,
          template: formTemplate
        });
        var formview = form.render();
        var _form_el = formview.el;
        var $form = formview.$el;
        
        function submitForm(values, callback) {
          var headers = {};
          headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
          delete values['comments']
          var post_data = {
            data: { 'plate_info': values },
            search_data: listView._options.search_data 
          }
          
          if (self.library && self.copy) {
            search_data = _.result(post_data,'search_data', {});
            search_data['library_short_name'] = self.library.get('short_name');
            search_data['copy_name'] = self.copy.get('copy_name');
            post_data['search_data'] = search_data;
          }            
          
          var url = self.resource.apiUri + '/batch_edit';
          
          var listParamString = listView.getCollectionUrl();
          if(!_.isUndefined(listParamString)
              && listParamString.indexOf('?')>-1) {
            url += '?' + listParamString.split('?')[1];
          }
          $.ajax({
            url: url,    
            data: JSON.stringify(post_data),
            cache: false,
            contentType: 'application/json',
            dataType: 'json', // what is expected back from the server
            processData: false, // do not process data being sent
            type: 'PATCH',
            headers: headers
          })
          .fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); })
          .done(function(data) {
            appModel.showConnectionResult(data, {
              title: 'Batch Edit Plates Success'
            });

            listView.collection.fetch();
            
            if (callback) {
              callback();
            }
          });
        }
        
        var dialog = appModel.showModal({
          okText: 'Submit',
          view: _form_el,
          title: 'Batch edit plate information',
          ok: function(e) {
            e.preventDefault();
            var errors = form.commit({ validate: true }) || {}; 
            var values = form.getValue();
            
            if (_.result(values, 'status', null) == null){
              delete values['status']
            }
            if (_.result(values, 'mg_ml_concentration', null) == null){
              delete values['mg_ml_concentration']
            }
            if (_.result(values, 'molar_concentration', null) == null){
              delete values['molar_concentration']
            }
            if (_.result(values, 'screening_count', null) == null){
              delete values['screening_count']
            }
            if (_.result(values, 'remaining_well_volume', null) == null){
              delete values['remaining_well_volume']
            }
            if (_.result(values, 'plate_type', null) == null){
              delete values['plate_type']
            }
            var valueTest = _.clone(values);
            delete valueTest['comments']
            if (_.isEmpty(valueTest)){
              errors['_others'] = [{'error': 'Must fill in at least one value'}];
            }
            form.$el.find('.form-group').removeClass('has-error');
            if (!_.isEmpty(errors) ) {
              _.each(_.keys(errors), function(key) {
                form.$el.find('[name="'+key +'"]').parents('.form-group').addClass('has-error');
                if (key=='_others') {
                  form.$el.append('<div class="text-danger">' + errors[key][0]['error'] + '</div>');
                }
              });
              return false;
            }            
            
            submitForm(values);
          }
        });
      };    
 
      console.log('call init...');
      $(this).queue([initfun]);
    },
    
    batchEditLocationsDialog: function(listView) {
      var self = this;
      console.log('batch edit locations dialog...');
      // find any extra search data
      post_data = {};
      if(_.has(listView._options, 'search_data')) {
        console.log('search data found on the list._options: ', listView._options.search_data);
        // endcode for the post_data arg; post data elements are sent 
        // as urlencoded values via a POST form for simplicity
        post_data['search_data'] = JSON.stringify(listView._options.search_data);  
      }
      // Construct the form
      var initfun = function() {
        var plateLocationTree = appModel.getPlateLocationTree();
        console.log('construct the batch edit form, ', plateLocationTree );
        var formSchema = {};
        var fieldTemplate = appModel._field_template;
        var formTemplate = appModel._form_template;
        formSchema['room'] = {
          title: 'Room',
          key: 'room',
          editorClass: 'form-control',
          type: Backbone.Form.editors.Text,
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['freezer'] = {
          title: 'Freezer',
          key: 'freezer',
          editorClass: 'form-control',
          type: Backbone.Form.editors.Text,
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['shelf'] = {
          title: 'Shelf',
          key: 'shelf',
          editorClass: 'form-control',
          type: Backbone.Form.editors.Text,
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['bin'] = {
          title: 'Bin',
          key: 'bin',
          editorClass: 'form-control',
          type: Backbone.Form.editors.Text,
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['comments'] = {
          title: 'Comments',
          key: 'comments',
          validators: ['required'],
          type: 'TextArea',
          editorClass: 'input-full',
          template: fieldTemplate
        };
        var FormFields = Backbone.Model.extend({
          schema: formSchema,
          validate: function(attrs) {
            var errs = {};
            if (!_.isEmpty(errs)) return errs;
          }
        });
        var formFields = new FormFields();
        
        var form = new Backbone.Form({
          model: formFields,
          template: formTemplate
        });
        var formview = form.render();
        var _form_el = formview.el;
        var $form = formview.$el;
        
        var subKey = 'room';
        $form.find('[name="'+subKey +'"]').typeahead({
          autoSelect: false,
          delay: 1,
          minLength: 0,
          items: 'all',
          source: _.keys(plateLocationTree),
          afterSelect: function(roomVal) {
            var subKey = 'freezer';
            form.setValue(subKey,null);
            form.setValue('shelf',null);
            form.setValue('bin',null);
            $form.find('[name="'+subKey +'"]').typeahead('destroy')
            $form.find('[name="'+subKey +'"]').typeahead({
              autoSelect: false,
              delay: 1,
              minLength: 0,
              items: 'all',
              source: _.keys(plateLocationTree[roomVal]),
              afterSelect: function(freezerVal) {
                var subKey = 'shelf';
                form.setValue(subKey,null);
                form.setValue('bin',null);
                $form.find('[name="'+subKey +'"]').typeahead('destroy')
                $form.find('[name="'+subKey +'"]').typeahead({
                  autoSelect: false,
                  delay: 1,
                  minLength: 0,
                  items: 'all',
                  source: _.keys(plateLocationTree[roomVal][freezerVal]),
                  afterSelect: function(shelfVal) {
                    var subKey = 'bin';
                    form.setValue(subKey,null);
                    $form.find('[name="'+subKey +'"]').typeahead('destroy')
                    $form.find('[name="'+subKey +'"]').typeahead({
                      autoSelect: false,
                      delay: 1,
                      minLength: 0,
                      items: 'all',
                      source: _.keys(plateLocationTree[roomVal][freezerVal][shelfVal])
                    });// bin
                  }
                });// shelf
              }
            });//freezer
          }
        });//room
        
        function findPlateLocation(room, freezer, shelf, bin) {
          return _.result(_.result(_.result(_.result(
            plateLocationTree, room),freezer),shelf),bin);
        }
        
        function submitForm(values, callback) {
          var headers = {};
          headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
          
          var post_data = {
            data: { 'plate_location': values },
            search_data: listView._options.search_data 
          }
          
          if (self.library && self.copy) {
            search_data = _.result(post_data,'search_data', {});
            search_data['library_short_name'] = self.library.get('short_name');
            search_data['copy_name'] = self.copy.get('copy_name');
            post_data['search_data'] = search_data;
          }            
          
          var url = self.resource.apiUri + '/batch_edit';
          
          var listParamString = listView.getCollectionUrl();
          if(!_.isUndefined(listParamString)
              && listParamString.indexOf('?')>-1) {
            url += '?' + listParamString.split('?')[1];
          }
          $.ajax({
            url: url,    
            data: JSON.stringify(post_data),
            cache: false,
            contentType: 'application/json',
            dataType: 'json', // what is expected back from the server
            processData: false, // do not process data being sent
            type: 'PATCH',
            headers: headers
          })
          .fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); })
          .done(function(data) {
            appModel.showConnectionResult(data, {
              title: 'Batch Edit Locations Success'
            });
            listView.collection.fetch();
            
            if (callback) {
              callback();
            }
          });
        }
        
        var dialog = appModel.showModal({
          okText: 'Submit',
          view: _form_el,
          title: 'Batch edit plates',
          ok: function(e) {
            e.preventDefault();
            var errors = form.commit({ validate: true }) || {}; 
            if (!_.isEmpty(errors) ) {
              _.each(_.keys(errors), function(key) {
                form.$el.find('[name="'+key +'"]').parents('.form-group').addClass('has-error');
              });
              return false;
            }            
            // TODO: construct the batch url and submit
            var values = form.getValue()
            
            if (_.isUndefined(findPlateLocation(
                values['room'],values['freezer'],values['shelf'],values['bin']))) {
              $('#modal').modal('hide');
              appModel.showModal({
                title: 'Plate location will be created',
                body: Iccbl.formatString(
                  'Plate location: {room}-{freezer}-{shelf}-{bin}',values),
                ok: function(e) {
                  submitForm(values, function(){
                    appModel.getPlateLocationTree(null, { flush:true });
                  });
                }
              });
              return false; // signal early exit for modal
            } else {
              submitForm(values);
            }
          }
        });
      };
      console.log('call init...');
      $(this).queue([appModel.getPlateLocationTree,initfun]);

    },
    
  });

  return LibraryCopyPlateView;
});