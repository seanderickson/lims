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
    'templates/generic-detail-layout.html',

], function($, _, Backbone, stickit, backbone_forms, Iccbl, appModel,
            DetailView, EditView, layoutTemplate ) {

  var DetailLayout = Backbone.Layout.extend({

    /**
     * @param resource - the API resource schema
     */
    initialize: function(args) {
      console.log('---- initialize detail layout', args);
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      this.subviews = {};
      this.args = args;
      this.DetailView = args.DetailView || DetailView;
      this.EditView = args.EditView || EditView;
      this.resource = args.resource || this.model.resource;
      this.modelFields = args.modelFields || this.resource.fields;
      this.title = args.title;
      
      _.bindAll(this, 'showDetail', 'showEdit');
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
            this.model.resource),
      }      
    },    
  
    template: _.template(layoutTemplate),
    
    showDetail: function() {
      var self = this;
      var view = this.subviews['detail'];
      if (!view) {
        view = new this.DetailView(
          _.extend({}, self.args, { model: this.model }));
        this.subviews['detail'] = view;
      }
      this.setView("#detail_content", view ).render();
      this.reportUriStack([]);
      return view;
    },
    
    clickEdit: function(event){
      event.preventDefault();
      this.reportUriStack(['edit']);
      this.showEdit();
    },

    showEdit: function() {
      console.log('showEdit: editView: ',EditView);
      var self = this;
      var view = new this.EditView(_.extend({}, self.args, 
        { 
          model: self.model, 
          uriStack: self.uriStack 
        }));
      Backbone.Layout.setupView(view);
      self.listenTo(view,'remove',function(){
        self.removeView(view);
        self.showDetail();
      });
      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      self.setView("#detail_content", view ).render();
      appModel.setPagePending();
      return view;
    },
    
    afterRender: function(){
      if (this.title){
        this.$el.find('#content_title').html(this.title).show();
        $('#content_title_row').show();
      }
      if (!_.isEmpty(this.uriStack)){
        viewId = this.uriStack.shift();
        if (viewId == 'edit') {
          this.uriStack.push(viewId);
          this.showEdit();
          return;
        }else if (viewId == '+add') {
          this.uriStack.push(viewId);
          this.showEdit();
          return;
        }
      }
      this.showDetail();
    },
     
    download: function(e){
      e.preventDefault();
      appModel.download(this.model.url, this.model.resource);
    },

    cancel: function(event){
      event.preventDefault();
      appModel.clearPagePending();
      this.showDetail();
    },    
    
    back: function(event){
      event.preventDefault();
      this.remove();
      appModel.router.back();
    },

    history: function(event) {
      event.preventDefault();
      var self = this;
      
      var newUriStack = ['apilog','order','-date_time', appModel.URI_PATH_SEARCH];
      var search = {};
      search['ref_resource_name'] = this.model.resource.key;

      //search['key'] = this.model.key;
      search['key'] = encodeURIComponent(this.model.key);
      
      newUriStack.push(appModel.createSearchString(search));
      var route = newUriStack.join('/');
      console.log('history route: ' + route);
      appModel.router.navigate(route, {trigger: true});
//      this.remove();
    },

    onClose: function() {
      this.subviews = {};
    },
    
    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      console.log('reportUriStack --- ', reportedUriStack, this.consumedStack);
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      this.trigger('uriStack:change', actualStack );
    }
    
	});
	
	return DetailLayout;
});
