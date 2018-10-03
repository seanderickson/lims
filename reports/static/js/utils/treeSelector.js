define([
  'jquery',
  'underscore',
  'backbone',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'templates/tree_selector.html',
  'jquery_bonsai'
], 
function($, _, Backgrid, Backbone, Iccbl, appModel,
    selector_tree_layout) {

  // NOTE: Webpack 3 patch:
  // Bind the Backbone.Layout object (and Backbone.Form) using the webpack 
  // assetsPluginInstance
  Backbone.Layout = LayoutManager;

  /**
   * Constructs a tree component for selecting objects from a collection.
   * @see https://github.com/aexmachina/jquery-bonsai
   * 
   * @param treeAttributes: object attributes to use for constructing the tree
   *   hierarchy.
   * 
   * TODO: (performance enhancement) consider rebuilding tree on each action
   * to contain only visible (searched, selected) nodes.
   */
  var TreeSelector = Backbone.Layout.extend({
    
    template: _.template(selector_tree_layout),
    /**
     * Name of the model property where the checked flag is stored.
     */
    PROP_CHECK_KEY: 'checked',
    
    initialize: function(options) {
      var self = this;
      console.log('initialize TreeSelector...', options);
      
      var options = options || {};
      self.options = options;
      self.collection = options.collection;
      self.treeAttributes = options.treeAttributes;
      self.extraControls = options.extraControls;
      if(_.has(options,'checked_key')){
        self.PROP_CHECK_KEY = options.checked_key;
      }
      
      /**
       * Build a data tree and an HTML list from the collection and the given
       * attributes.
       */
      function buildTree(collection, treeAttributes, chKey) {
        var dataTree = {}
        
        function makeNodes(dataTree, attrs, i, model){
          var attr = model.get(attrs[i]);
          
          if (_.isUndefined(attr)){
            // TODO: some columns in the collection may not be represented
            // with all the needed field schema attributes specified by treeAttributes
            // (this is a bug).
            // e.g. cherryPickRequest.setCherryPickPlates, the "selected" virtual
            // column is not meant to be seen in the treeSelector and does not have
            // all of the field schema fields.
            return;
          }
          
          attr = attr.replace(/\"/g,"'").trim();
          var data_name = _.map(
            attrs.slice(0,i+1), function(key){
              return model.get(key).replace(/\"/g,"'").trim();
            }).join('/');
          var node = _.result(dataTree, attr, { 
            text: attr,
            data_name: data_name,
            nodes: {}
          });
          if (attrs.length == i+1 ){
            node.id = collection.modelId(model);
            if (model.get(chKey) === true){
              node.checked = true;
            }
            if (model.has('description')){
              node.description = model.get('description');
            }
          } else {
            makeNodes(node.nodes, attrs, i+1, model)
          }
          dataTree[attr] = node;
        }
        collection.each(function(model){
          makeNodes(dataTree, treeAttributes, 0, model);
        });
        
        function makeList(nodes){
          var html = ''
          _.each(_.values(nodes), function(node){
            html += '<li data-name="' + node.data_name + '" ';
            if (node.id){
              html += ' id="' + node.id + '" ';
            }
            if (node.description){
              html += ' title="' + node.description +'" ';
            }
            if (node.checked===true){
              html += ' data-checked=true ';
            }
            
            html += ' >&nbsp;' + node.text;
            if (!_.isEmpty(node.nodes)){
              html += '<ol>';
              html += makeList(node.nodes)
              html += '</ol>';
            }
            html += '</li>';
          });
          return html;
        }
        html = makeList(dataTree);
        return { dataTree: dataTree, html: html };
      }; // buildTree
      self.treeData = buildTree(
        self.collection, self.treeAttributes, self.PROP_CHECK_KEY);
      
      console.log('initialized TreeSelector');
    },
    
    afterRender: function(){
      var self = this;
      var chKey = this.PROP_CHECK_KEY;

      var $treeControl = this.$treeControl = $('<ol id="tree"></ol>');
      var $selectedTreeControl = $('<ol id="selectedTree"></ol>');
      var $search = $('#treeSearch');
      var $searchButton = $('#treeSearchButton');
      var $clearSearchButton = $('#clearSearchButton');
      var $selectSearchHitsButton = $('#selectSearchHits');
      var $unselectSearchHitsButton = $('#unselectSearchHits');
      var $toggleExpanded = $('#toggleExpanded');
      var treeAttributes = self.treeAttributes;
      
      if(self.extraControls && !_.isEmpty(self.extraControls)){
        _.each(self.extraControls, function(control){
          $('#left_panel_controls').prepend(control);
        });
      }
      
      var $unselectSelectedButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
         id="unselectSelected" title="unselect the selected items" >Unselect</button>');
      $('#right_title').append($unselectSelectedButton);
      
      $treeControl.html(self.treeData.html);
      $("#left_panel").append($treeControl);
      
      $selectedTreeControl.html(self.treeData.html);
      $("#right_panel").append($selectedTreeControl);

      // TODO: allow customization of typeahead keys: currently only the top
      // level keys are used
      var typeAheadKeys = {};
      self.collection.each(function(model){
        var attributes = self.treeAttributes;
        if (self.options.treeAttributesForTypeAhead){
          attributes = self.options.treeAttributesForTypeAhead;
        }
        _.each(attributes, function(attr){
          typeAheadKeys[model.get(attr)] = attr;
        });
      });
//      $search.typeahead({
//        source: _.keys(self.treeData.dataTree)
//      });
      $search.typeahead({
        source: _.keys(typeAheadKeys)
      });
      
      $treeControl.bonsai({ 
        expandAll: false,
        checkboxes: true, // depends on jquery.qubit plugin
        handleDuplicateCheckboxes: true,         
        createInputs: 'checkbox'
      });
      var bonsai = $treeControl.data('bonsai');
      $selectedTreeControl.bonsai({ 
        expandAll: true,
        checkboxes: true, // depends on jquery.qubit plugin
        handleDuplicateCheckboxes: true,         
        createInputs: 'checkbox'
      });
      var bonsaiSelected = $selectedTreeControl.data('bonsai');

      function toggleExpanded(){
        var state = $toggleExpanded.attr('aria-pressed');
        $toggleExpanded.attr('aria-pressed',!state);
        if (state == 'true'){
          bonsai.collapseAll();
          $toggleExpanded.empty();
          $toggleExpanded.append('Expand all');
        } else {
          bonsai.expandAll();
          $toggleExpanded.empty();
          $toggleExpanded.append('Collapse')
        }
      };
      function showChecked($bonsaiTree){
        var items = $bonsaiTree.find('li');
        for (i=0; i<items.length; i++){
          var obj = $(items[i]);
          var id = obj.attr('id');
          var model = self.collection.get(id);
          if (model && model.get(chKey)===true){
            obj.show();
            obj.parents().show();
          } else {
            obj.hide();
          }
        }
      };

      $toggleExpanded.click(toggleExpanded);
      if (self.options.startExpanded === true){
        toggleExpanded();
      }
      
      showChecked($selectedTreeControl);

      $search.on('keyup', function (e) {
        if (e.keyCode == 13) {
          var searchedModels = self.search();
          self.collection.trigger('searchChange', searchedModels);
        }
      });      
      $searchButton.click(function(){
        var searchedModels = self.search();
        self.collection.trigger('searchChange', searchedModels);
      });

      $selectSearchHitsButton.click(function(){
        var searchedModels = self.getSearched();
        _.each(searchedModels, function(model){
          model.set(chKey,true, {silent: true});
        });
        if (!_.isEmpty(searchedModels)){
          // "block" does not work in modal
          // set loading before firing event; single thread blocks dom updates
          $('#treeloading').css( "display", "table" ); 
          window.setTimeout(function(){
            self.collection.trigger('bulkCheck',searchedModels);
          }, 100);
        }
        return false;
      });

      $unselectSearchHitsButton.click(function(){
        var searchedModels = self.getSearched();
        _.each(searchedModels, function(model){
          model.set(chKey,false, {silent: true});
        });
        if (!_.isEmpty(searchedModels)){
          // "block" does not work in modal
          // set loading before firing event; single thread blocks dom updates
          $('#treeloading').css( "display", "table" ); 
          window.setTimeout(function(){
            self.collection.trigger('bulkCheck',searchedModels);
          }, 100);
        }
      });

      $unselectSelectedButton.click(function(){
        var changedModels = self.collection.where({checked: true});
        _.each(changedModels, function(model){
          model.set(chKey, false, { silent: true });
        });
        // "block" does not work in modal
        // set loading before firing event; single thread blocks dom updates
        $('#treeloading').css( "display", "table" ); 
        window.setTimeout(function(){
          self.collection.trigger('bulkCheck',changedModels);
        }, 100);
      });

      $clearSearchButton.click(function(){
        $search.val('');
        // check if searchedModels is empty (sub class may override search)
        var searchedModels = self.search();
        if (_.isEmpty(searchedModels)){
          self.clearSearch();
        } else {
          self.collection.trigger('searchChange', searchedModels);
        }
      });

      // Store checked keys on treeChanged for batch update
      var bulkCheck = [];
      function treeChanged(e){
        var id = $(e.target).parent().attr('id');
        var model = self.collection.get(id);
        if (model){
          var prop = {};
          var selectedLi = $selectedTreeControl.find('#'+ id);
          if ($(e.target).prop(chKey) === true){
            prop[chKey] = true;
            model.set(prop, {silent: true});
            // set UI elements instead of triggering a change on the collection or elements
            // TODO: does this help?
            selectedLi.show();
            selectedLi.parents().show();
            // defer select notifications
            bulkCheck.unshift(model);
            //selectedLi.find('input').prop(chKey, true); //.trigger('change');
          } else {
            prop[chKey] = false;
            model.set(prop, {silent: true});
            // defer select notifications
            bulkCheck.unshift(model);
            //selectedLi.find('input').prop(chKey, false); //.trigger('change');
          }
          // defer notification for performance
          window.setTimeout(function(){
            self.collection.trigger('bulkCheck', bulkCheck.slice());
            bulkCheck.length = 0;
          },100);
        }
      };
      $treeControl.on('change', treeChanged);
      $selectedTreeControl.on('change', treeChanged);
      
      self.collection.on('change:'+chKey, function(model){
        var id = self.collection.modelId(model)
        var li = $treeControl.find('#'+ id)
        var selectedLi = $selectedTreeControl.find('#'+ id)
        if (li){
          if (model.get(chKey)===true){
            li.find('input').prop(chKey, true).trigger('change');
            bonsai.expandTo(id);
          } else {
            li.find('input').prop(chKey, false).trigger('change');
          }
        }
        if (selectedLi){
          if (model.get(chKey)===true){
            selectedLi.show();
            selectedLi.parents().show();
            selectedLi.find('input').prop(chKey, true).trigger('change');
          } else {
            selectedLi.find('input').prop(chKey, false).trigger('change');
          }
        }
        //showChecked($selectedTreeControl);
      });
      
      self.collection.on('bulkCheck', function(changedModels){
        if (_.isEmpty(changedModels)) return;
        // NOTE: setting display to "block" does not work in modal
        // NOTE2: event processing is single threaded, so "loading" will not 
        // display if set here; must be set in the event handler parent
        //$('#treeloading').css( "display", "table" ); 
        var storedParents1 = {};
        var storedParents2 = {};
        _.each(changedModels, function(model){
          var id = self.collection.modelId(model)
          var li = $treeControl.find('#'+ id)
          var selectedLi = $selectedTreeControl.find('#'+ id)
          if (li){
            if (model.get(chKey)===true){
              li.find('input').prop('checked', true);
              // NOTE: defer notification for performance
              //li.find('input').prop(chKey, true).trigger('change');
            } else {
              li.find('input').prop('checked', false);
            }
            storedParents1[model.get(self.treeAttributes[0])] = li;
          }
          if (selectedLi){
            if (model.get(chKey)===true){
              selectedLi.find('input').prop('checked', true);
            } else {
              selectedLi.find('input').prop('checked', false);
            }
            storedParents2[model.get(self.treeAttributes[0])] = selectedLi;
          }
        });
        $treeControl.bonsai('update');
        showChecked($selectedTreeControl);
        $selectedTreeControl.bonsai('update');
        
        // NOTE: nofify just once per branch of leaf nodes for performance
        // Send one change event per child; so that bonsai can render parent checkboxes
        _.each(_.values(storedParents1), function(li){
          li.find('input').trigger('change');
        });
        _.each(_.values(storedParents2), function(li){
          li.find('input').trigger('change');
        });
        $('#treeloading').css( "display", "none" );
      });
      
      self.collection.on('searchChange', function(searchedModels){

        $toggleExpanded.attr('aria-pressed',true);
        $toggleExpanded.empty();
        $toggleExpanded.append('Collapse');

        // NOTE: 20180926 - limit "select all" availability for performance:
        // TODO: improve tree selector performance: note, however, that 
        // selecting too many columns may also make the query non-performant.
        // A better solution may be to pass in a "max selections" setting, that
        // is based on the capabilities for the query servicing the list view.
        var MAX_SELECTIONS = 40;
        if (!_.isEmpty(searchedModels) ){
          if (_.size(searchedModels) < MAX_SELECTIONS ){
            $selectSearchHitsButton.show();
          }
          $unselectSearchHitsButton.show();
        } else {
          $selectSearchHitsButton.hide();
          $unselectSearchHitsButton.hide();
        }
        self.collection.each(function(model){
          model.unset('searchHit', {silent: true})
        });

        var ids = _.map(searchedModels, function(model){
          model.set({'searchHit': true}, {silent: true});
          return self.collection.modelId(model);
        });
        var items = $treeControl.find('li');
        for (i=0; i<items.length; i++){
          var obj = $(items[i]);
          
          if (_.contains(ids, obj.attr('id'))){
            obj.show();
            obj.parents().show();
            bonsai.expandTo(obj);
          } else {
            obj.hide();
          }
        }        
      });

    },
    
    getSearchVal: function(){
      return $('#treeSearch').val().trim()
    },

    clearSearch: function(){
      var self = this;
      var $search = $('#treeSearch');
      var $selectSearchHitsButton = $('#selectSearchHits');
      var $unselectSearchHitsButton = $('#unselectSearchHits');
      var $toggleExpanded = $('#toggleExpanded');
      var $treeControl = this.$treeControl;
      var bonsai = $treeControl.data('bonsai');
      
      $search.val('');
      this.collection.each(function(model){
        model.unset('searchHit', {silent: true})
      });
      $treeControl.find('li').show();
      bonsai.collapseAll();
      _.each(this.collection.where({checked: true}), function(model){
        var id = self.collection.modelId(model);
        bonsai.expandTo(id);
      });
      
      $selectSearchHitsButton.hide();
      $unselectSearchHitsButton.hide();
      
      $toggleExpanded.attr('aria-pressed',false);
      $toggleExpanded.empty();
      $toggleExpanded.append('Expand all');
    },
    
    /**
     * Search the designated "treeAttributes" of the collection's models.
     */
    search: function(){
      var self = this;
      var pattern = self.getSearchVal();
      var collection = self.collection;
      var treeAttributes = self.treeAttributes;
      searchedModels = [];
      if (_.isEmpty(pattern)){
        return [];
      }
      pattern = pattern.replace(/[-[\]{}()*+?.,\\^$|#\s]/g, '\\$&');
      var reSearch = new RegExp('.*' + pattern + '.*', "i");
      
      searchedModels = collection.filter(
        function(model){
          return _.find(treeAttributes, function(attr){
            return reSearch.exec(model.get(attr));
          });
        });
      return searchedModels;
    },
    
    getSearched: function(){
      return this.collection.where({searchHit: true});
    }, 
    
    getSelected: function(){
      return this.collection.where({checked: true});
    },
    
    getTreeControl: function(){
      return this.$treeControl;
    }
    
  });
  return TreeSelector;
  
});