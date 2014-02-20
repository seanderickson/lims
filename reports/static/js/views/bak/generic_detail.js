define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'iccbl_backgrid',
    'text!templates/generic-detail.html',
    'text!templates/modal_ok_cancel.html',

], function( $, _, Backbone, stickit, Iccbl, genericDetailTemplate, modalOkCancel ) {
    var DetailView = Backbone.View.extend({


        events: {
            'click button#delete': 'delete',
            'click button#edit': 'edit',
            'click button#back': 'back',
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
            // this._id = this.model.get('id');
            // if(_.has(this._schemaResult['resource_definition'], 'title_attribute')){
                // console.log('create title_attribute from ' + this._schemaResult['resource_definition']['title_attribute']);
                // this._id = _.reduce(this._schemaResult['resource_definition']['title_attribute'],
                    // function(memo, item){
                        // if( self.model.has(item) ) memo += self.model.get(item)
                        // else memo += item
                        // return memo ;
                    // }, '');
            // }
            console.log('id: ' + this._id);

            var self=this;

            this._keys = Iccbl.sortOnOrdinal(_.keys(this.model.attributes), self._schemaResult.fields);
            this._keys = _(this._keys).filter(function(key){
                return _.has(self._schemaResult.fields[key], 'visibility') && _.contains(self._schemaResult.fields[key]['visibility'], 'detail');
            });
            console.log('detail keys: ' + JSON.stringify(this._keys));

            _.bindAll(this, 'render');
            this.model.on('change', this.render);

        },


        edit: function(event) {
            console.log('!! edit not implemented !!');
            // TODO
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

            this.$el.empty();
            this.trigger('remove');
            this._router.back();
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

            var compiledTemplate = _.template( genericDetailTemplate,
                { 'fieldDefinitions': this._schemaResult.fields ,
                  title: this._resource_definition['title']+ ': ' + this._id,
                  keys: _(this._keys),
                  object: this.model.attributes
                });

            this.$el.html(compiledTemplate);

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