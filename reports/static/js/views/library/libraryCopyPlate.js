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
  
  var VIEWS = {
    'ListView': ListView, 
    'DetailView': DetailLayout
  };
    
  var LibraryCopyPlateView = Backbone.Layout.extend({
    
    template: _.template(layout),
    
    initialize: function(args) {
      this.resource = args.resource || args.options.resource;
      this.uriStack = args.uriStack || args.options.uriStack;
      this.library = args.library  || args.options.library;
      this.copy = args.copy || args.options.copy;
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
    
    afterRender: function(){
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
        if (self.library && self.copy){
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
        buttons: ['download','history'],
        title: 'Plate ' + model.get('plate_number')
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView('#resource_content', view).render();
    },

    showList: function(resource, url) {
      var self = this;
      var uriStack = _.clone(this.uriStack);
      console.log('uriStack', uriStack);
      var showBatchEditButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="batch_edit_button" href="#">',
          'Batch Edit</a>'
        ].join(''));
      var showHistoryButton = $([
      '<a class="btn btn-default btn-sm pull-down" ',
        'role="button" id="showHistoryButton" href="#">',
        'History</a>'
      ].join(''));
      showHistoryButton.click(function(e){
        e.preventDefault();
        var newUriStack = ['apilog','order','-date_time', 'search'];
        var search = {};
        search['ref_resource_name'] = 'librarycopyplate';
        if (self.library && self.copy){
          search['key__icontains'] = [
            self.library.get('short_name'),
            self.copy.get('name')].join('/');
        }
        newUriStack.push(appModel.createSearchString(search));
        var route = newUriStack.join('/');
        console.log('history route: ' + route);
        appModel.router.navigate(route, {trigger: true});
        self.remove();
      });

      var view = new ListView({ options: {
        uriStack: uriStack,
        schemaResult: resource,
        resource: resource,
        url: url,
//        collection: collection,
        extraControls: [showBatchEditButton, showHistoryButton]
      }});
      showBatchEditButton.click(function(e){
        e.preventDefault();
        self.createBatchEditDialog(view);
      });
        
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      Backbone.Layout.setupView(view);
      self.setView('#resource_content', view ).render();
    },
    
    createBatchEditDialog: function(listView){
      var self = this;
      console.log('batch edit dialog...');
      // find any extra search data
      post_data = {};
      if(_.has(listView._options, 'search_data')){
        console.log('search data found on the list._options: ', listView._options.search_data);
        // endcode for the post_data arg; post data elements are sent 
        // as urlencoded values via a POST form for simplicity
        post_data['search_data'] = JSON.stringify(listView._options.search_data);  
      }
      // Construct the form
      var initfun = function(){
        var plateLocationTree = appModel.getPlateLocationTree();
        console.log('construct the batch edit form, ', plateLocationTree );
        var formSchema = {};
        var fieldTemplate = appModel._field_template;
        var formTemplate = appModel._form_template;
        roomOptions = [{val: '', label: '' }];
        _.each(_.keys(plateLocationTree),function(room){
            roomOptions.push({ val: room, label: room });
        });
        formSchema['room'] = {
          title: 'Room',
          key: 'room',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          options: roomOptions,
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['freezer'] = {
          title: 'Freezer',
          key: 'freezer',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          options: [],
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['shelf'] = {
          title: 'Shelf',
          key: 'shelf',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          options: [],
          validators: ['required'],
          template: fieldTemplate 
        };
        formSchema['bin'] = {
          title: 'Bin',
          key: 'bin',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          options: [],
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
          validate: function(attrs){
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
        formview.$el.find('.chosen-select').chosen({
          disable_search_threshold: 3,
          width: '100%',
          allow_single_deselect: true,
          search_contains: true
        });

        form.on("room:change", function(e){
          var room = form.getValue('room');
          console.log('change:room ' + room);
          var options = [{val: '', label: '' }];
          _.each(_.keys(plateLocationTree[room]),function(val){
              options.push({ val: val, label: val });
          });
          var fieldKey = 'freezer';
          var $chosen = form.$el.find('[name="' + fieldKey + '"]').parent()
              .find('.chosen-select');
          $chosen.empty();
          _.each(options, function(option){
            $chosen.append($('<option>',{
                value: option.val
              }).text(option.label));
          });
          $chosen.trigger("chosen:updated");
        });
        
        form.on("freezer:change", function(e){
          var room = form.getValue('room');
          var freezer = form.getValue('freezer');
          console.log('change:freezer ', room, freezer);
          var options = [{val: '', label: '' }];
          _.each(_.keys(plateLocationTree[room][freezer]),function(val){
              options.push({ val: val, label: val });
          });
          var fieldKey = 'shelf';
          var $chosen = form.$el.find('[name="' + fieldKey + '"]').parent()
              .find('.chosen-select');
          $chosen.empty();
          _.each(options, function(option){
            $chosen.append($('<option>',{
                value: option.val
              }).text(option.label));
          });
          $chosen.trigger("chosen:updated");
        });
        
        form.on("shelf:change", function(e){
          var room = form.getValue('room');
          var freezer = form.getValue('freezer');
          var shelf = form.getValue('shelf');
          console.log('change:shelf', room, freezer, shelf);
          var options = [{val: '', label: '' }];
          _.each(_.keys(plateLocationTree[room][freezer][shelf]),function(val){
              options.push({ val: val, label: val });
          });
          var fieldKey = 'bin';
          var $chosen = form.$el.find('[name="' + fieldKey + '"]').parent()
              .find('.chosen-select');
          $chosen.empty();
          _.each(options, function(option){
            $chosen.append($('<option>',{
                value: option.val
              }).text(option.label));
          });
          $chosen.trigger("chosen:updated");
        });
        
        var dialog = appModel.showModal({
          okText: 'Submit',
          view: _form_el,
          title: 'Batch edit plates',
          ok: function(e){
            e.preventDefault();
            var errors = form.commit({ validate: true }) || {}; 
            if(!_.isEmpty(errors) ){
              _.each(_.keys(errors), function(key){
                form.$el.find('[name="'+key +'"]').parents('.form-group').addClass('has-error');
              });
              return false;
            }            
            // TODO: construct the batch url and submit
            var values = form.getValue()
            var headers = {};
            headers[appModel.HEADER_APILOG_COMMENT] = values['comment'];
            
            var post_data = {
              data: values,
              search_data: listView._options.search_data 
            }
            
            if (self.library && self.copy){
              search_data = _.result(post_data,'search_data', {});
              search_data['library_short_name'] = self.library.get('short_name');
              search_data['copy_name'] = self.copy.get('name');
              post_data['search_data'] = search_data;
            }            
            
            var url = self.resource.apiUri + '/batch_edit';
            
            var listParamString = listView.getCollectionUrl();
            if(!_.isUndefined(listParamString)
                && listParamString.indexOf('?')>-1){
              url += '?' + listParamString.split('?')[1];
            }
            $.ajax({
              url: url,    
              data: JSON.stringify(post_data),
              cache: false,
              contentType: 'application/json',
              processData: false,
              type: 'PATCH',
              headers: headers
            })
            .fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); })
            .done(function(model, resp){
                // FIXME: should be showing a regular message
                appModel.error('success');
            });
            
         }
        });
      };
      console.log('call init...');
      $(this).queue([appModel.getPlateLocationTree,initfun]);

    }
    
  });

  return LibraryCopyPlateView;
});