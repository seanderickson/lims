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
      'click button#download': 'download',
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
      this.listenTo(view,'remove',function(){
        this.showDetail();
      });
      // NOTE: no subviews, so can notify on render
      this.trigger('uriStack:change',[]);
      this.setView("#detail_content", view ).render();
    },
    
    afterRender: function(){
      this.showDetail();
    },
 
    
    download: function(event){
      $('#tmpFrame').attr('src', this.model.url + '?format=csv' );
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
      newUriStack.push('order/-date_time');
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
