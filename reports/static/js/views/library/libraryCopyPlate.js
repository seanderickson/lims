define([
  'jquery',
  'underscore',
  'backbone',
  'layoutmanager',
  'iccbl_backgrid',
  'models/app_state',
  'views/list2',
  'views/generic_detail_layout',
  'views/generic_edit',
  'templates/genericResource.html'
], 
function($, _, Backbone, layoutmanager, Iccbl, appModel, ListView, DetailLayout, 
         EditView, layout) {
  
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

      if (!_.isEmpty(uriStack) && !_.isEmpty(uriStack[0]) &&
              !_.contains(appModel.LIST_ARGS, uriStack[0]) ) {
        // Detail view
        var plate_number = uriStack.shift();
        this.consumedStack = [plate_number];
        _key = library.key + '/' + copy.get('name')+ '/' + plate_number;
        appModel.getModel(resourceId, _key, this.showDetail );
      } else {
        // List view
        var url = resource.apiUri;
        if (self.library && self.copy) {
          url = [ 
              library.resource.apiUri,library.key,'copy',copy.get('name'),
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
      var view = new DetailLayout({ 
        model: model, 
        uriStack: uriStack,
//        buttons: ['download','history'],
        title: 'Plate ' + model.get('plate_number')
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#resource_content', view).render();
    },

    showList: function(resource, url) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
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
      showHistoryButton.click(function(e) {
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', 'search'];
        var search = {};
        search['ref_resource_name'] = 'librarycopyplate';
        if (self.library && self.copy) {
          search['key__icontains'] = [
            self.library.get('short_name'),
            self.copy.get('name')].join('/');
        }
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        console.log('history route: ' + route);
        appModel.router.navigate(route, {trigger: true});
      });

      // TODO: extending the passed args to get the "search_data", or other 
      // passed args, grab options explicitly instead 201608
      options = _.extend({
        uriStack: uriStack,
        schemaResult: resource,
        resource: resource,
        url: url,
        extraControls: [
          showEditLocationButton, showEditPlatesButton, showHistoryButton]
        }, self._args ) ;
      var view = new ListView(options);
      showEditLocationButton.click(function(e) {
        e.preventDefault();
        self.batchEditLocationsDialog(view);
      });
      showEditPlatesButton.click(function(e) {
        e.preventDefault();
        self.batchEditPlatesDialog(view);
      });
      self.listenTo(view, 'update_title', function(val) {
        if(val) {
          this.$('#content_title').html(resource.title + ': <small>' + val + '</small>');
        } else {
          this.$('#content_title').html(resource.title);
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
        formSchema['well_volume'] = {
          title: 'Initial Well Volume',
          key: 'well_volume',
          type: EditView.SIunitEditor,
          template: fieldTemplate 
        };
        _.extend(
          formSchema['well_volume'],
          self.resource['fields']['well_volume']['display_options']);
        formSchema['comments'] = {
          title: 'Comments',
          key: 'comments',
          validators: ['required'],
          type: 'TextArea',
          template: fieldTemplate
        };
        
        var FormFields = Backbone.Model.extend({
          schema: formSchema,
          validate: function(attrs) {
            var errs = {};
            if (attrs['status']==null 
                && attrs['plate_type']==null
                && attrs['well_volume']==0 ) {
              errs['error'] = 'Must enter Status, Plate Type, or Initial Well Volume';
            }
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
            search_data['copy_name'] = self.copy.get('name');
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
            listView.collection.fetch();
            
            if (_.isObject(data) && !_.isString(data)) {
              data = _.result(_.result(data,'meta',data),'Result',data);
              var msg_rows = appModel.dict_to_rows(data);
              var bodyMsg = msg_rows;
              if (_.isArray(msg_rows) && msg_rows.length > 1) {
                bodyMsg = _.map(msg_rows, function(msg_row) {
                  return msg_row.join(': ');
                }).join('<br>');
              }
              var title = 'Upload success';
              appModel.showModalMessage({
                body: bodyMsg,
                title: title  
              });
            } else {
              console.log('ajax should have been parsed data as json', data);
              appModel.showModalMessage({
                title: 'Batch edit success',
                okText: 'ok',
                body: 'Records: ' + listView.collection.state.totalRecords + ', ' + data
              });
            }
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
            if (!_.isEmpty(errors) ) {
              _.each(_.keys(errors), function(key) {
                form.$el.find('[name="'+key +'"]').parents('.form-group').addClass('has-error');
                if (key=='_others') {
                  form.$el.append('<div class="text-danger">' + errors[key][0]['error'] + '</div>');
                }
              });
              return false;
            }            
            var values = form.getValue()
            submitForm(values);
          }
        });
      };    
 
      console.log('call init...');
      $(this).queue([initfun]);
    },
    
    batchEditLocationsDialog: function(listView) {
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
            search_data['copy_name'] = self.copy.get('name');
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
            listView.collection.fetch();
            
            if (_.isObject(data) && !_.isString(data)) {
              data = _.result(_.result(data,'meta',data),'Result',data);
              var msg_rows = appModel.dict_to_rows(data);
              var bodyMsg = msg_rows;
              if (_.isArray(msg_rows) && msg_rows.length > 1) {
                bodyMsg = _.map(msg_rows, function(msg_row) {
                  return msg_row.join(': ');
                }).join('<br>');
              }
              var title = 'Upload success';
              appModel.showModalMessage({
                body: bodyMsg,
                title: title  
              });
            } else {
              console.log('ajax should have been parsed data as json', data);
              appModel.showModalMessage({
                title: 'Batch edit success',
                okText: 'ok',
                body: 'Records: ' + listView.collection.state.totalRecords + ', ' + data
              });
            }
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