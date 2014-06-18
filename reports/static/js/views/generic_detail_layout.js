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
	  
	  initialize: function(args) {
	    console.log('---- initialize detail layout');
      this.uriStack = args.uriStack;
      this.consumedStack = [];
	    this.subviews = {};
      _.bindAll(this, 'showDetail');
      _.bindAll(this, 'showEdit');
	  },

	  events: {
      'click button#download': 'download',
      'click button#delete': 'delete',
      'click button#edit': 'clickEdit',
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
      this.setView("#detail_content", view ).render();
    },
    
    clickEdit: function(event){
      event.preventDefault();
      this.reportUriStack(['edit']);
      this.showEdit();
    },
    
    showEdit: function() {
      var view = this.subviews['edit'];
      if (!view) {
        view = new EditView({ model: this.model, uriStack: ['edit'] });
        Backbone.Layout.setupView(view);
        this.subviews['edit'] = view;
      }
      this.listenTo(view,'remove',function(){
        this.showDetail();
      });
      
      this.setView("#detail_content", view ).render();
    },    
    
    showAdd: function() {
      var view = this.subviews['add'];
      if (!view) {
        view = new EditView({ model: this.model, uriStack: ['add'] });
        Backbone.Layout.setupView(view);
        this.subviews['add'] = view;
      }
      this.listenTo(view,'remove',function(){
        this.showDetail();
      });
      this.setView("#detail_content", view ).render();
    },
    
    afterRender: function(){
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if (viewId == 'edit') {
          this.showEdit();
          return;
        }else if (viewId == '+add') {
          this.showAdd();
          return;
        }
      }
      this.showDetail();
    },
     
    download: function(event){
      $('#tmpFrame').attr('src', this.model.url + '?format=csv' );
    },

    cancel: function(event){
      event.preventDefault();
      this.remove();
      appModel.router.back();
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
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    }
    
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
