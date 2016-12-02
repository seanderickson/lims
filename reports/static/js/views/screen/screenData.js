define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'layoutmanager',
  'models/app_state',
  'views/list2',
  'views/screen/libraryScreening',
  'utils/tabbedController'
], function($, _, Backbone, Backgrid, Iccbl, layoutmanager, appModel, 
            ListView, LibraryScreeningView, TabbedController) {
  
  var ScreenDataView = TabbedController.extend({
    
    initialize: function(args) {
      var self = this;
      this._classname = 'ScreenDataView';
      TabbedController.prototype.initialize.apply(this,arguments);
    },

    tabbed_resources: {
      detail: { 
        description: 'Result Data', 
        title: 'Screen Results', 
        invoke: 'setDetail'
      },
      datacolumns: {
        description : 'Data Columns',
        title : 'Data Columns',
        invoke : 'setDatacolumns',
        permission: 'screenresult'
      }
    },
    study_tabbed_resources: {
      detail: { 
        description: 'Reagent Annotations', 
        title: 'Reagent Annotations', 
        invoke: 'setDetail'
      },
      datacolumns: {
        description : 'Data Columns',
        title : 'Data Columns',
        invoke : 'setDatacolumns',
        permission: 'screenresult'
      }
    },
    
    events: {
      'click ul.nav-tabs >li': 'click_tab',
      'click button#loadScreenResults': 'loadScreenResults',
      'click button#deleteScreenResults': 'deleteScreenResults',
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      var displayed_tabbed_resources = _.extend({},this.tabbed_resources);
      if (self.model.get('project_phase') == 'annotation') {
        var displayed_tabbed_resources = _.extend({},this.study_tabbed_resources);
      }
      if (! self.model.get('has_screen_result')) {
        delete displayed_tabbed_resources['datacolumns'];
      }
      return {
        'base_url': self.model.resource.key + '/' + self.model.key + '/data',
        'tab_resources': displayed_tabbed_resources
      }      
    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function() {
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)) {
        viewId = this.uriStack.shift();

        if (_.contains(appModel.LIST_ARGS, viewId) ) {
          this.uriStack.unshift(viewId);
          viewId = 'detail'
        }

        if (!_.has(this.tabbed_resources, viewId)) {
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },    
    
    setDetail: function(delegateStack) {

      console.log('set screen result detail...', delegateStack);
      
      var self = this;
      var key = 'detail';

      var $loadScreenResultsButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
        id="loadScreenResults" >Load Data</button>');
      var $deleteScreenResultsButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
        id="deleteScreenResults" >Delete Data</button>');
      
      
      var screenResultResource = appModel.getResource('screenresult');
      var schemaUrl = [appModel.dbApiUri,
                       'screenresult',
                       self.model.key,
                       'schema'].join('/');
      var url = [appModel.dbApiUri,
                 'screenresult',
                 self.model.key].join('/');

      var _id = self.model.key;
      var screen_facility_id = self.model.get('facility_id');
      
      var createResults = function(schemaResult) {
        var extraControls = [];
        var show_positives_control = $([
          '<label class="checkbox-inline">',
          '  <input type="checkbox">positives',
          '</label>'
          ].join(''));
        var show_mutual_positives_control = $([
          '<label class="checkbox-inline">',
           '  <input type="checkbox">mutual positives',
           '</label>'
           ].join(''));
        if (self.model.get('project_phase') != 'annotation') {
          extraControls = extraControls.concat(
            show_positives_control, show_mutual_positives_control);
          extraControls.push($('<span>&nbsp;</span>'));
        }
        if (appModel.hasPermission('screenresult','write')) {
          extraControls = extraControls.concat(
            $deleteScreenResultsButton,$loadScreenResultsButton);
        }
        
        // create an option vocab for the exclued col, if needed
        if (_.has(schemaResult['fields'], 'excluded')) {
          var options = [];
          _.each(schemaResult['fields'],function(field) {
            var dc_col_key = 'dc_' + self.model.key + '_';
            if (field['key'].indexOf(dc_col_key)>-1) {
              var option_name = field['key'].substring(dc_col_key.length);
              options.unshift([option_name, option_name]);
            }
          });
          schemaResult['fields']['excluded']['vocabulary'] = options;
          console.log('custom exclude options', options);
        }
        var initialSearchHash;
        view = new ListView({ 
          uriStack: _.clone(delegateStack),
          schemaResult: schemaResult,
          resource: screenResultResource,
          url: url,
          extraControls: extraControls
        });
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);
        self.listenTo(view , 'uriStack:change', self.reportUriStack);
        self.setView("#tab_container", view ).render();
        initialSearchHash = view.listModel.get('search');
        
        if (_.has(initialSearchHash, 'is_positive__eq')
            && initialSearchHash.is_positive__eq.toLowerCase()=='true') {
          show_positives_control.find('input[type="checkbox"]').prop('checked',true);
        }
        if (_.has(initialSearchHash, 'show_mutual_positives')
            && initialSearchHash.show_mutual_positives.toLowerCase()=='true') {
          show_mutual_positives_control.find('input[type="checkbox"]').prop('checked',true);
          view.show_mutual_positives(screen_facility_id, true);
        }
        if (_.has(initialSearchHash, 'other_screens')) {
          view.show_other_screens(initialSearchHash['other_screens']);
        }
        show_positives_control.click(function(e) {
          if (e.target.checked) {
            var searchHash = _.clone(view.listModel.get('search'));
            searchHash['is_positive__eq'] = 'true';
            view.listModel.set('search',searchHash);
          } else {
            // make sure unset
            var searchHash = _.clone(view.listModel.get('search'));
            if (_.has(searchHash,'is_positive__eq')) {
              delete searchHash['is_positive__eq'];
              view.listModel.set('search',searchHash);
            }
          }
        });
        show_mutual_positives_control.click(function(e) {
          if (e.target.checked) {
            window.setTimeout(function() {
              view.show_mutual_positives(screen_facility_id, true);
            });
          } else {
            window.setTimeout(function() {
              view.show_mutual_positives(screen_facility_id, false);
            });
          }
        });
        
      };
      
      if (self.model.get('has_screen_result')) {
        appModel.getResourceFromUrl(schemaUrl, createResults);
      } else {
        if (appModel.hasPermission('screenresult','write')) {
          this.$('#tab_container').append("No data loaded ");
          this.$('#tab_container').append($loadScreenResultsButton);
          self.reportUriStack([]);
        }
      }
    },
    
    /**
     * Loads the screen results using ajax to post the file.
     * - because ajax cannot handle post response attachments (easily), set
     * the response type to 'application/json' and display the errors in a 
     * modal dialog.
     * 
     * NOTE: a "POST" form cannot be used - there is no standard way to signal
     * the result of the post to JavaScript.
     */
    loadScreenResults: function() {
      var self = this;
      var isStudy = self.model.get('project_phase') == 'annotation';
      var form_template = [
         "<div class='form-horizontal container' id='screenresult_form' >",
         "<form data-fieldsets class='form-horizontal container' >",
         "<div class='form-group' ><input type='file' name='fileInput' /></div>",
         "</form>",
         "</div>"].join('');      
      
      var fieldTemplate = _.template([
        '<div class="form-group" >',
        '    <label class="control-label " for="<%= editorId %>"><%= title %></label>',
        '    <div class="" >',
        '      <div data-editor  style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
        '      <div data-error class="text-danger" ></div>',
        '      <div><%= help %></div>',
        '    </div>',
        '  </div>',
      ].join(''));
      
      var formSchema = {};
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        validators: ['required'],
        editorClass: 'input-full',
        type: 'TextArea',
        template: fieldTemplate
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema,
        validate: function(attrs) {
          var errs = {};
          var file = $('input[name="fileInput"]')[0].files[0]; 
          if (! file) {
            errs.fileInput = 'Must specify a file';
          }
          if (!_.isEmpty(errs)) return errs;
        }
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields,
        template: _.template(form_template)
      });
      var _form_el = form.render().el;

      function saveFile() {

        var errors = form.commit({ validate: true }); // runs schema and model validation
        if (!_.isEmpty(errors)) {
          console.log('form errors, abort submit: ',errors);
          _.each(_.keys(errors), function(key) {
            $('[name="'+key +'"').parents('.form-group').addClass('has-error');
            if ( key == '_others') {
              _.each(errors[key], function(otherError) {
                // TODO: display file req'd msg
              });
            }
          });
          return false;
        } else {
          var url = [self.model.resource.apiUri,self.model.key,'screenresult'].join('/');
          var values = form.getValue();
          var file = $('input[name="fileInput"]')[0].files[0]; 
          var data = new FormData();
          // NOTE: 'xls' is sent as the key to the FILE object in the upload.
          // Use this as a non-standard way to signal the upload type to the 
          // serializer. 
          data.append('xls', file);
          var headers = {
            'Accept': 'application/json'
          };
          headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
          $.ajax({
            url: url,    
            data: data,
            cache: false,
            contentType: false,
            processData: false,
            type: 'POST',
            headers: headers
          })
          .fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); })
          .done(function(model, resp) {
            // FIXME: should be showing a regular message
            appModel.error('success');
            // Must refetch and render, once for the containing tabbed layout,
            // and then (during render), will refetch again for the table view
            self.model.fetch({
              success: function() {
                if (!isStudy) {
                  self.uriStack = ['detail'];
                  // remove the child view before calling render, to prevent
                  // it from being rendered twice, and calling afterRender twice
                  self.removeView('#tab_container');
                }
                self.render();
              }
            }).fail(function() { 
              Iccbl.appModel.jqXHRfail.apply(this,arguments); 
            });
          });
        }        
      }
      
      var dialog = appModel.showModal({
          ok: saveFile,
          view: _form_el,
          title: 'Select a Screen Result (xlsx workbook) file to upload'
      });
    },
    
    deleteScreenResults: function() {
      var self = this;
      var isStudy = self.model.get('project_phase') == 'annotation';
      var formSchema = {};
      formSchema['comments'] = {
        title: 'Comments',
        key: 'comments',
        editorClass: 'input-full',
        validators: ['required'],
        type: 'TextArea'
      };

      var FormFields = Backbone.Model.extend({
        schema: formSchema
      });
      var formFields = new FormFields();
      var form = new Backbone.Form({
        model: formFields
      });
      var _form_el = form.render().el;
      var dialog = appModel.showModal({
          ok: function() {
            var errors = form.commit({ validate: true }); 
            if (!_.isEmpty(errors)) {
              _.each(_.keys(errors), function(key) {
                $('[name="'+key +'"').parents('.form-group').addClass('has-error');
              });
              return false;
            } else {
              var url = [self.model.resource.apiUri,self.model.key,'screenresult'].join('/');
              var values = form.getValue();
              var headers = {
                'Accept': 'application/json'
              };
              headers[appModel.HEADER_APILOG_COMMENT] = values['comments'];
              $.ajax({
                url: url,    
                cache: false,
                contentType: false,
                processData: false,
                type: 'DELETE',
                headers: headers
              })
              .fail(function() { Iccbl.appModel.jqXHRfail.apply(this,arguments); })
              .done(function(model, resp) {
                appModel.error('success');
                // Must refetch and render, once for the containing tabbed layout,
                // and then (during render), will refetch again for the detail view
                self.model.fetch({
                  success: function() {
                    if (!isStudy) {
                      self.uriStack = ['detail'];
                      // remove the child view before calling render, to prevent
                      // it from being rendered twice
                      self.removeView('#tab_container');
                    }
                    self.render();
                  }
                }).fail(function() { 
                  Iccbl.appModel.jqXHRfail.apply(this,arguments); 
                });
              });
            }        
          },
          view: _form_el,
          title: 'Delete screen results for ' + self.model.get('facility_id')
      });
      
    },
    
    setDatacolumns: function(delegateStack) {
      
      var self = this;
      var datacolumnResource = appModel.getResource('datacolumn'); 
      var url = [self.model.resource.apiUri,self.model.key,'datacolumns'].join('/');
      
      // construct a collection-in-columns
      // pivot the datacolumn list
      Iccbl.getCollectionOnClient(url, function(collection) {
        // create a colModel for the list
        var columns = [];
        var TextWrapCell = Backgrid.Cell.extend({
          className: 'text-wrap-cell'
        })
        var colTemplate = {
          'cell' : 'string',
          'order' : -1,
          'sortable': false,
          'searchable': false,
          'editable' : false,
          'visible': true,
          'headerCell': Backgrid.HeaderCell.extend({
            render: function() {
              this.$el.empty();
              var label = Iccbl.createLabel(
                this.column.get("label"), 10);
              this.$el.append(label);
              return this;
            }
          })
        };
        // setup the first column - schema field titles (that are visible)
        var col0 = _.extend({},colTemplate,{
          'name' : 'fieldname',
          'label' : '',
          'description' : 'Datacolumn field',
          'cell': TextWrapCell.extend({
            className: '150_width'
          })
        });
        columns.push(col0);
        
        collection.each(function(datacolumn) {
          var col = _.extend({},colTemplate,{
            'name' : datacolumn.get('key'),
            'label' : datacolumn.get('title'),
            'description' : datacolumn.get('description'),
            'order': datacolumn.get('ordinal'),
            'sortable': true,
            'cell': TextWrapCell
          });
          columns.push(col);
          if (datacolumn.has('derived_from_columns')) {
            datacolumn.set(
              'derived_from_columns', 
              datacolumn.get('derived_from_columns').join(', '));
          }
        });
        var colModel = new Backgrid.Columns(columns);
        colModel.comparator = 'order';
        colModel.sort();
        
        // create the collection by pivoting
        orderedFields = _.sortBy(datacolumnResource.fields,'ordinal');
        var pivotCollection = new Backbone.Collection();
        _.each(orderedFields, function(field) {
          if (_.contains(field.visibility, 'l') 
              && !_.contains(['name','title','key'],field.key )) {
            var row = {'key': field.key, 'fieldname': field.title };
            collection.each(function(datacolumn) {
              row[datacolumn.get('key')] = datacolumn.get(field.key);
            });
            pivotCollection.push(row);
          }
        });
        
        // now show as a list view
        var view = new Backgrid.Grid({
          columns: colModel,
          collection: pivotCollection,
          className: 'backgrid table-striped table-condensed table-hover'
        });
        var view = view.render();
        Backbone.Layout.setupView(view);
        self.reportUriStack([]);

        self.setView("#tab_container", view ).render();
     });
      
    }    

  });
  
  return ScreenDataView;
});