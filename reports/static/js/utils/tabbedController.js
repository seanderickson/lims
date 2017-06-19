define([
  'jquery',
  'underscore',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'views/generic_edit',
  'templates/generic-tabbed.html',
  'views/generic_detail_layout'
], 
function($, _, Backgrid, Iccbl, appModel, EditView, tabbedTemplate, DetailLayout) {

  var TabbedController = Backbone.Layout.extend({
    
    template: _.template(tabbedTemplate),

    initialize: function(args) {
      var self = this;
      this.tabViews = {}; 
      this.uriStack = args.uriStack;
      this.consumedStack = [];
      
      var displayed_tabbed_resources = _.extend({},this.tabbed_resources);
      _.each(_.keys(displayed_tabbed_resources), function(key) {
        if (key !== 'detail') {
          var permission = displayed_tabbed_resources[key].permission;
          if (_.isUndefined(permission)) {
            permission = displayed_tabbed_resources[key].resource;
          }
          if (!appModel.hasPermission(permission)) {
            delete displayed_tabbed_resources[key];
          }
        }
      });
      this.tabbed_resources = displayed_tabbed_resources;
      
      _.bindAll(this, 'click_tab','change_to_tab');
    },

    tabbed_resources: {
      // To be defined in the subclass
    },      
    
    // NOTE: if overriding, must register click_tab in the subclass
    events: {
      'click ul.nav-tabs >li': 'click_tab',
    },

    /**
     * Child view bubble up URI stack change event
     */
    reportUriStack: function(reportedUriStack) {
      var consumedStack = this.consumedStack || [];
      var actualStack = consumedStack.concat(reportedUriStack);
      console.log('reportUriStack: actualStack', actualStack, arguments);
      this.trigger('uriStack:change', actualStack );
    },
    
    /**
     * Layoutmanager hook
     */
    serialize: function() {
      var self = this;
      return {
        'base_url': self.model.resource.key + '/' + self.model.key ,
        'tab_resources': self.tabbed_resources
      }      
    }, 
    
    /**
     * Layoutmanager hook
     */
    afterRender: function() {
      var viewId = 'detail';
      if (!_.isEmpty(this.uriStack)) {
        viewId = this.uriStack.shift();

        if(viewId == '+add'){
          console.log('disable tabs during add operation...');
          this.$('ul.nav-tabs > li').addClass('disabled');
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }
        if(viewId == 'edit'){
          this.uriStack.unshift(viewId); 
          viewId = 'detail';
        }

        if (!_.has(this.tabbed_resources, viewId)) {
          var msg = 'could not find the tabbed resource: ' + viewId;
          appModel.error(msg);
          throw msg;
        }
      }
      this.change_to_tab(viewId);
    },
    
    click_tab : function(event) {
      var self = this;
      event.preventDefault();
      event.stopPropagation();
      var key = event.currentTarget.id;
      if (_.isEmpty(key)) return;
      if (this.$('#'+key).hasClass('disabled')) {
        return;
      }
      appModel.requestPageChange({
        ok: function() {
          self.change_to_tab(key);
        }
      });
      
    },

    change_to_tab: function(key) {
      var self = this;
      console.log('change_to_tab: ' + key);
      appModel.clearErrors();
      
      if (_.has(this.tabbed_resources, key)) {
        this.removeView("#tab_container");
        this.$('li').removeClass('active');
        this.$('#' + key).addClass('active');
        if (key !== 'detail') {
          this.consumedStack = [key];
        } else {
          this.consumedStack = [];
        }
        var delegateStack = _.clone(this.uriStack);
        this.uriStack = [];
        var method = this[this.tabbed_resources[key]['invoke']];
        if (_.isFunction(method)) {
          method.apply(this,[delegateStack]);
        } else {
          throw ( 
            "Tabbed resources refers to a non-function: "
            + this.tabbed_resources[key]['invoke'] );
        }
      } else {
        var msg = 'Unknown tab: ' + key;
        appModel.error(msg);
        throw msg;
      }
    },

    /**
     * TODO: not used?//
     */
    showAdd: function() {
      console.log('add view');
      
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack
      });
      Backbone.Layout.setupView(view);

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },
    
    /**
     * TODO: not used?//
     */
    showEdit: function() {
      var self = this;
      var delegateStack = _.clone(this.uriStack);
      var view = new DetailLayout({
        model: self.model,
        uriStack: delegateStack, 
        buttons: ['download']
      });
      Backbone.Layout.setupView(view);

      self.listenTo(view , 'uriStack:change', self.reportUriStack);
      this.setView("#tab_container", view ).render();
      this.$('li').removeClass('active');
      this.$('#detail').addClass('active');
    },
    
    setDetail: function(delegateStack) {
      // must be implemented

      // Example, set title
      var title = self.model.resource.title;
      if (!self.model.isNew()){
        title += ': ' + Iccbl.getTitleFromTitleAttribute(self.model, self.model.resource);
      }
      $title = this.$el.find('#tab_container-title');
      $title.html(title);
      $title.show();
      
    },    
    
    onClose: function() {
      // TODO: is this necessary when using Backbone LayoutManager
      this.tabViews = {};
      this.remove();
    }

  });
  
  return TabbedController;
});



