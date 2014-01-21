define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'backbone_forms',
    'iccbl_backgrid',
    'text!templates/generic-detail.html',
    'text!templates/generic-form-backbone-forms.html',
    'text!templates/modal_ok_cancel.html',

], function( $, _, Backbone, stickit, backbone_forms, Iccbl, genericDetailTemplate, genericFormTemplate, modalOkCancel ) {
    //Backbone.Form = backbone_forms;
	
	var DetailView = Backbone.View.extend({
//		onDomRefresh: function() {
//		  this.focusFirstInput();
//		},
//
//		focusFirstInput: function() {
//			console.log('focusFirstInput' + this.$(':input:visible:enabled:first') );
//		  this.$(':input:visible:enabled:first').focus();
//		},

        events: {
            'click button#save': 'save',
            'click button#delete': 'delete',
            'click button#edit': 'edit',
            'click button#cancel': 'cancel',
            'click button#history': 'history'
        },

        initialize : function(attributes, options) {
            var self = this;
            Iccbl.requireOptions(options, ['schemaResult', 'router']);
            Iccbl.assert( !_.isUndefined(options.schemaResult['resource_definition']), 'detail view schemaResult requires a resource_definition');

            this._schemaResult = options['schemaResult'];
            this._resource_definition = this._schemaResult['resource_definition'];
            this._router = options.router;
            this._options = options;

            this._id = Iccbl.getTitleFromTitleAttribute(self.model, this._schemaResult);

            console.log('id: ' + this._id);

            this._keys = Iccbl.sortOnOrdinal(_.keys(this.model.attributes), self._schemaResult.fields);
            _.bindAll(this, 'render');
            this.model.on('change', this.render);

        },

        detail: function(event) {
            // filter keys for detail view
            this._options.isEditMode = false; // TODO: lets split edit to a different view?
            var self=this;
            var detailKeys = _(this._keys).filter(function(key){
                return _.has(self._schemaResult.fields[key], 'visibility') && _.contains(self._schemaResult.fields[key]['visibility'], 'detail');
            });
            console.log('detail keys: ' + JSON.stringify(detailKeys));

            var compiledTemplate = _.template( genericDetailTemplate,
                { 'fieldDefinitions': this._schemaResult.fields ,
                  title: this._resource_definition['title']+ ': ' + this._id,
                  keys: _(detailKeys),
                  object: this.model.attributes
                });

            this.$el.html(compiledTemplate);

            this.listenTo(this.model, 'error', this.error);
            console.log('initialize DetailView');
        },

        edit: function(event) {
            this._options.isEditMode = true;
            //console.log(' template: ', genericFormTemplate, 'fields: ', this._options.fields['screensaver_user_id']['title'], ', k: ' , this._keys, ', t: ', this._options.title );
            bindings = {};
            var self = this;
            var editKeys = _(this._keys).filter(function(key){
                return _.has(self._schemaResult.fields[key], 'visibility') 
                	&& _.contains(self._schemaResult.fields[key]['visibility'], 'edit');
            });
            
         // like 'Select' editor, but will always return a boolean (true or false)
            Backbone.Form.editors.BooleanSelect = Backbone.Form.editors.Select.extend({
                initialize: function(options) {
                    options.schema.options = [
                        { val: 'true', label: 'Yes' },
                        { val: 'false', label: 'No' }
                    ];
                    Backbone.Form.editors.Select.prototype.initialize.call(this, options);
                },
                getValue: function() {
                    return !!Backbone.Form.editors.Select.prototype.getValue.call(this);
                },
                setValue: function(value) {
                    value = value ? 'true' : 'false';
                    Backbone.Form.editors.Select.prototype.setValue.call(this, value);
                }
            });            
            

            // build the backbone-forms schema
            
            // using custom templates to hold the editors, 
            // control the layout with the "controls, control-group" classes
            var altFieldTemplate = _.template('\
            	<div class="control-group"> \
            	      <label class="control-label" for="<%= editorId %>"><%= title %></label>\
            	      <div class="controls" >\
            	        <span data-editor></span>\
            	        <div data-error></div>\
            	        <div><%= help %></div>\
            	      </div>\
            	    </div>\
            	  ');       
            
            // using custom templates to hold the editors, 
            // control the layout with the "controls, control-group" classes
            var altRadioFieldTemplate = _.template('\
            	<div class="control-group"> \
          	      <label class="control-label" for="<%= editorId %>"><%= title %></label>\
          	      <div class="controls" >\
          	        <span data-editor></span>\
          	        <div data-error></div>\
          	        <div><%= help %></div>\
          	      </div>\
          	    </div>\
          	  ');       
            

            
            var schema = {};
            var itemcount = 0;
            _.each(editKeys, function(key){
                if( _(self._schemaResult.fields).has(key)){
                    option = self._schemaResult.fields[key];
                    console.log('option: ' + JSON.stringify(option));
                    if(option.ui_type == 'Select' 
                    	|| option.ui_type == 'Radio'
                    	|| option.ui_type == 'Checkboxes' ){
                        var _optionsCollection = [];
                        if(_.has(option, 'choices')){
                            _optionsCollection = option.choices.map(function(choice){
                                return choice;
                            });
                        }else{
                            window.alert('Warning, no choices defined for: ' + key); // TODO: use bootstrap alerts!
                            option.choices = _(_optionsCollection); // so the template doesn't complain
                        }
                    	schema[key] = { 
                    			type: option.ui_type, 
                    			options: _optionsCollection
                			};


                        if(option.ui_type == 'Checkboxes' ){ 
                        	schema[key]['defaults'] = self.model.get(key);
                        } 
                    }else if( option.ui_type == 'boolean'){
                    	schema[key] = {
                    			type: 'Checkbox'
                    	}
                    }else if( option.ui_type.toLowerCase() == 'date'){
                        	schema[key] = {
                        			type: 'Date'
                        	}
                    }else{
                        schema[key] = {
                        		type: 'Text'
                    		}
                    }
                	if(option.ui_type == 'Radio'){
//                    	schema[key]['template'] = altRadioFieldTemplate;
                    }else{
                    	schema[key]['template'] = altFieldTemplate;
                    }

                    if(itemcount++ == 0){
                		// Set autofocus (HTML5) on the first field
                		// NOTE: see http://stackoverflow.com/questions/20457902/how-to-automatically-focus-first-backbone-forms-input-field
                		// - we may want to revisit this for a more robust solution
                		schema[key]['editorAttrs'] = { autofocus: 'autofocus'}
                	}
                }
            });
            console.log('compile template...');
            var compiledTemplate = _.template( genericFormTemplate );  //,
            var GenericForm = Backbone.Form.extend({
            	template: compiledTemplate,
            	schema: schema
            	
            });
            console.log('setup the form and render...');
            self._form = new GenericForm({
            	model: self.model,
            	templateData: { 
            		title: this._resource_definition['title']+ ': ' + this._id,
            		keys: _(editKeys),
            		fieldDefinitions: this._schemaResult.fields 
        		}
            });
            console.log('render....');
            self._form.render();
            this.$el.html(self._form.el);
        },

        save: function(event){
            event.preventDefault();
            var self = this;
            self._form.commit();
            console.log('save: ' + JSON.stringify(this.model));
            this.model.save(null, {
                success: function(model, resp){
                    console.log('success');
                    self.$el.empty();
                    //self.trigger('remove');
                    self.detail(null);

                },
                error: function(model,xhr){
                    var re = /([\s\S]*)Request Method/;
                    var match = re.exec(xhr.responseText);
                    if (match) window.alert('error saving: ' + match[1] + ': ' + xhr.status + ':' + xhr.statusText );
                    else window.alert('error on saving: ' + xhr.status + ':' + xhr.statusText );
                }

            });
        },

        delete: function(event){
            event.preventDefault();
            var self = this;
            console.log('delete: ' + JSON.stringify(this.model));
            var model = this.model;
            var modalDialog = new Backbone.View({
                el: _.template(modalOkCancel, { body: "Please confirm deletion of record: '" + model.get('toString') + "'", title: "Please confirm deletion" } ),
                events: {
                    'click #modal-cancel':function(event) {
                        console.log('cancel button click event, '); // + JSON.stringify(fieldDefinitions));
                        event.preventDefault();
                        $('#modal').modal('hide'); // TODO: read-up on modal!  this is not ideal with the reference to template elements!
                    },
                    'click #modal-ok':function(event) {
                        console.log('ok button click event, '); // + JSON.stringify(fieldDefinitions));
                        event.preventDefault();
                        model.destroy();
                        $('#modal').modal('hide');
                        self.$el.empty();
                        self.trigger('remove');
                        self._router.back();
                    }
                },
            });
            modalDialog.render();
            $('#modal').empty();
            $('#modal').html(modalDialog.$el);
            $('#modal').modal();
        },

        cancel: function(event){
            event.preventDefault();
            // TODO: do we have to do anything to abort form changes?
            // this.close(); // can we do this
            // this.$el.remove();
            if(this._options.isEditMode){
                this.detail(null);
            }else{
                this.$el.empty();
                this.trigger('remove');
                this._router.back();
            }
        },

        history: function(event){
            console.log('history click');
            event.preventDefault();
            this.$el.empty();
            this.trigger('remove');
            var self = this;

            // construct the route to the history
            // TODO: consider changing the app state to get to the records (future benefit of not reloading models if cached?)
            var _history_search = "ref_resource_name=" + this._resource_definition['key'];
            var id = this.model.get('id');
            if(_.has(this._resource_definition, 'id_attribute')){
                //console.log('create id from ' + this._resource_definition['id_attribute']);
                id = _.reduce(this._resource_definition['id_attribute'],
                        function(memo, item){
                            if(!_.isEmpty(memo)) memo += '/';
                            return memo += self.model.get(item);}, '');
            }else{
                console.log('Warn: schema for this type has no resource_definition,id_attribute; type: ' + this._options.type);
            }
            _history_search += ',key='+id;
            // TODO: discuss whether full encoding of the search fragment is necessary.
            // to-date, we know that the forward slash messes up the backbone router parsing, but other URL chars do not,
            // and full encoding reduces the usability of the URL for the end user
            //_history_search = encodeURIComponent(_history_search);
            _history_search = _history_search.replace('\/','%2F');
            //console.log('history_search:' + _history_search);

            var _route = 'list/apilog/search/' + _history_search; //ref_resource_name=' + _resource_name + ',key='+ id;
            console.log('-- set route: ' + _route);
            this._router.navigate(_route, {trigger: true});
        },

        error: function(options){
            console.log('A model save error reported: ' + JSON.stringify(options));
        },

        render : function() {
            console.log('render detail_stickit');

            if(this._options.isEditMode){
                // template = genericFormTemplate;
                this.edit(null);
            }else{
                this.detail(null);
            }

            return this;
        },

        onClose: function(){
            console.log('...detail view close method');
            //this.$el.remove();
            //this.modelBinder.unbind();
        }
    });

    return DetailView;
});