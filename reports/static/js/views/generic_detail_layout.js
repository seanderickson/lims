define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_stickit',
    'backbone_forms',
    'iccbl_backgrid',
    'models/app_state',
    'views/generic_detail_stickit',
    'views/generic_edit',
    'text!templates/generic-detail-layout.html',

], function($, _, Backbone, stickit, backbone_forms, Iccbl, appModel,
            DetailView, EditView, layoutTemplate ) {
	
	var DetailLayout = Backbone.Layout.extend({
	  
	  initialize: function() {
	    console.log('---- initialize detail layout');
	    this.subviews = {};
      _.bindAll(this, 'showDetail');
	  },
    events: {
      'click button#save': 'save',
      'click button#delete': 'delete',
      'click button#edit': 'showEdit',
      'click button#cancel': 'cancel',
      'click button#back': 'back',
      'click button#history': 'history'
    },
    
    serialize: function() {
      return {
        'title': Iccbl.getTitleFromTitleAttribute(
            this.model,
            this.model.resource.schema),
      }      
    },    
  
    template: _.template(layoutTemplate),
    
    showDetail: function() {
      var view = this.subviews['detail'];
      if (!view) {
        view = new DetailView({ model: this.model});
        this.subviews['detail'] = view;
      }
      // NOTE: no subviews, so can notify on render
      this.trigger('uriStack:change',[]);
      this.setView("#detail_content", view ).render();
    },
    
    showEdit: function() {
      var view = this.subviews['edit'];
      if (!view) {
        view = new EditView({ model: this.model});
        Backbone.Layout.setupView(view);
        this.subviews['edit'] = view;
      }
      // NOTE: no subviews, so can notify on render
      this.trigger('uriStack:change',[]);
      this.setView("#detail_content", view ).render();
    },
    
    afterRender: function(){
      this.showDetail();
    },
    
    save: function( event ) {
      event.preventDefault();
      var self = this;
      
      var view = this.subviews['edit'];
      if ( ! view ) { 
        window.alert("edit view not found");
      }else{
        view.commit();
        
        // Fixup the URL - if it points to the model instance, make it point to 
        // the API resource only: tastypie wants this
        // Note: this is happening if the model was fetched specifically for this
        // page, and has the url used to fetch it, rather than the collection url.
        var key = Iccbl.getIdFromIdAttribute( self.model,self.model.resource.schema );
        var url = _.result(this.model, 'url');
        ////    if ( url && url.indexOf(key) != -1 ) {
        ////    url = url.substring( 0,url.indexOf(key) );
        ////  }        
        
        // TODO: check if creating new or updating here
        // set the id specifically on the model: backbone requires this to 
        // determine whether a "POST" or "PATCH" will be used
        this.model.id = key;

        this.model.save(null, {
          url: url // set the url property explicitly
        })
        .done(function(model, resp){
          // done replaces success as of jq 1.8
          console.log('model saved');
        })
        .fail(function(model,xhr){
          // fail replaces error as of jquery 1.8
          if ( _.has(xhr,'responseJSON' ) ) {
            if ( _.has( xhr.responseJSON, 'error_message') ) {
              window.alert( xhr.responseJSON.error_message );
            }
            else if ( _.has( xhr.responseJSON, 'error') ) {
              window.alert( xhr.responseJSON.error_message );
            }
          } else {
            var re = /([\s\S]*)Request Method/;
            var match = re.exec(xhr.responseText);
            if (match) {
              window.alert(
                'error saving: ' + match[1] + ': ' + xhr.status + ':' +
                xhr.statusText );
            } else {
              window.alert('error on saving: ' + xhr.status + ':' + 
                           xhr.statusText );
            }              
          }
        })
        .always(function() {
          // always replaces complete as of jquery 1.8
          self.showDetail();
        });
      }
    },    

    cancel: function(event){
      event.preventDefault();
      if ( this.getView("#detail_content") == this.subviews['detail'] ) {
        this.remove();
        appModel.back();
      } else {
        this.showDetail();
      }
    },    
    
    back: function(event){
      event.preventDefault();
      this.remove();
      appModel.router.back();
    },

    history: function(event) {
      event.preventDefault();
      var self = this;
      
      var newUriStack = ['apilog','search'];
      var search = {};
      search['ref_resource_name'] = this.model.resource.key;
      search['key'] = this.model.key;
      newUriStack.push(
          _.map(
            _.pairs(search), 
            function(keyval) {
              return keyval.join('=');
            }).join(','));
      
      var route = newUriStack.join('/');
      console.log('history route: ' + route);
      appModel.router.navigate(route, {trigger: true});
      this.remove();
    },

    onClose: function() {
      this.subviews = {};
    },
    
	});
	
	return DetailLayout;
});




//    delete: function(event){
//      event.preventDefault();
//      var self = this;
//      console.log('delete: ' + JSON.stringify(this.model));
//      var model = this.model;
//      var modalDialog = new Backbone.View({
//        el: _.template(
//            modalOkCancel, { body: "Please confirm deletion of record: '" + 
//            model.get('toString') + "'", title: "Please confirm deletion" } ),
//        events: {
//          'click #modal-cancel':function(event) {
//              event.preventDefault();
//              $('#modal').modal('hide'); // TODO: read-up on modal!
//                                          // this is not ideal with
//                                          // the reference to template
//                                          // elements!
//          },
//          'click #modal-ok':function(event) {
//              event.preventDefault();
//              model.destroy();
//              $('#modal').modal('hide');
//              self.$el.empty();
//              self.trigger('remove');
//              self._router.back();
//          }
//        },
//      });
//      modalDialog.render();
//      $('#modal').empty();
//      $('#modal').html(modalDialog.$el);
//      $('#modal').modal();
//    },
