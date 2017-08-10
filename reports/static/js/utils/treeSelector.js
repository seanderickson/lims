define([
  'jquery',
  'underscore',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'templates/tree_selector.html',
  'jquery_bonsai'
], 
function($, _, Backgrid, Iccbl, appModel,
    selector_tree_layout) {

  /**
   * FIXME: 20170809 - THIS FILE IS IN DEVELOPMENT:
   * TODO: refactor form controls into the template
   * TODO: refactor search method so that external clients may override search
   * (Screen - other columns view: filter mutual positive screens, filter
   * positive columns, filter study columns)
   */
  
  /**
   * Constructs a tree component for selecting objects from a collection
   * @param treeAttributes: object attributes to use for constructing the tree
   *   hierarchy.
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
        this.PROP_CHECK_KEY = options.checked_key;
      }
      console.log('initialized TreeSelector');
    },
    
    afterRender: function(){
      var self = this;
      var chKey = this.PROP_CHECK_KEY;

      
//      var $searchBlock = $([
//        '<div class="input-group input-group-sm">',
//        ' <input id="treeSearch" class="form-control" placeholder="Search...">',
//        '<div class="input-group-btn">',
//        ' <button type="button" class="btn btn-default btn-sm" ',
//        '   id="treeSearchButton" title="Search">',
//        '   <span class="glyphicon glyphicon-search"></span></button>',
//        ' <button type="button" class="btn btn-default btn-sm" ',
//        '   id="clearSearchButton" label="Clear" title="Clear Search">',
//        '   Clear</button>',
//        '</div>',
//        '</div>'].join(''));
      var $search = $('#treeSearch');
      var $searchButton = $('#treeSearchButton');
      var $clearSearchButton = $('#clearSearchButton');
      var $selectSearchHitsButton = $('#selectSearchHits');
//      var $selectSearchHitsButton = $(
//        '<button class="btn btn-default btn-sm" role="button" \
//         style="display: none; " title="select the searched items" \
//         id="selectSearchHits" >Select</button>');
      var $unselectSearchHitsButton = $('#unselectSearchHits');
//      var $unselectSearchHitsButton = $(
//        '<button class="btn btn-default btn-sm" role="button" \
//         style="display: none;" title="unselect the searched items" \
//         id="unselectSearchHits" >Unselect</button>');
      var $toggleExpanded = $('#toggleExpanded');
//      var $toggleExpanded = $(
//        '<button class="btn btn-default btn-sm" data-toggle="button" \
//          aria-pressed="false" autocomplete="off" \
//          title="expand (collapse) the entire tree" >Expand</button>');
      
//      $("#left_title").append($searchBlock);
//      $('#left_panel').append($selectSearchHitsButton);
//      $('#left_panel').append($unselectSearchHitsButton);
//      $('#left_panel').append($toggleExpanded);
      
      if(self.extraControls && !_.isEmpty(self.extraControls)){
        _.each(self.extraControls, function(control){
          $('#left_panel_controls').append(control);
        });
      }
      
      var $unselectSelectedButton = $(
        '<button class="btn btn-default btn-sm" role="button" \
         id="unselectSelected" title="unselect the selected items" >Unselect</button>');
      $('#right_title').append($unselectSelectedButton);
      
      var $treeControl = this.$treeControl = $('<ol id="tree"></ol>');
      var $selectedTreeControl = $('<ol id="selectedTree"></ol>');
      
      var treeAttributes = self.treeAttributes;
      var startChecked = [];
      
      function buildTree(collection) {
        var dataTree = {}
        
        function makeNodes(dataTree, attrs, i, model){
          var attr = model.get(attrs[i]);
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
              startChecked.push(node.id);
            }
          } else {
            //attrs = attrs.slice(1);
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
        // console.log('html', html);
        return { dataTree: dataTree, html: html };
      }//
      
      var treeData = buildTree(self.collection);
      $treeControl.html(treeData.html);
      $("#left_panel").append($treeControl);
      
      $selectedTreeControl.html(treeData.html);
      $("#right_panel").append($selectedTreeControl);
      $search.typeahead({
        source: _.keys(treeData.dataTree)
      })
      
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

      _.each(startChecked, function(id){
        bonsai.expandTo(id);
      });
      $toggleExpanded.click(function(){
        var state = $toggleExpanded.attr('aria-pressed');
        $toggleExpanded.attr('aria-pressed',!state);
        console.log('state', state);
        if (state == 'true'){
          bonsai.collapseAll();
          $toggleExpanded.empty();
          $toggleExpanded.append('Expand');
        } else {
          bonsai.expandAll();
          $toggleExpanded.empty();
          $toggleExpanded.append('Collapse')
        }
      });

      function showChecked($bonsaiTree){
        var bonsaiCtrl = $bonsaiTree.data('bonsai');
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
      }
      showChecked($selectedTreeControl);
      
      var searchedModels = [];
      function searchNodes($bonsaiTree){
        searchedModels.length = 0;
        var pattern = $search.val().trim();
        
        if (_.isEmpty(pattern)){
          return;
        }
        
        console.log('search for "' + pattern + '"');
        var reSearch = new RegExp('.*' + pattern + '.*', "i");
        var items = $bonsaiTree.find('li');
        for (i=0; i<items.length; i++){
          var obj = $(items[i]);
          var data_name = obj.attr('data-name');
          if (data_name){
            if (reSearch.exec(data_name)){
              obj.show();
              obj.parents().show();
              bonsai.expandTo(obj);
              var id = obj.attr('id');
              var model = self.collection.get(id);
              if (model){
                searchedModels.push(model);
              }
            }else{
              obj.hide();
            }
          }
        }
        if (!_.isEmpty(searchedModels)){
          $selectSearchHitsButton.show();
          $unselectSearchHitsButton.show();
        } else {
          $selectSearchHitsButton.hide();
          $unselectSearchHitsButton.hide();
        }
        console.log('search:"' + pattern + '", models:', searchedModels.length);
      }
      $search.on('keyup', function (e) {
        if (e.keyCode == 13) {
          searchNodes($treeControl);
        }
      });      
      $searchButton.click(function(){
        searchNodes($treeControl);
      });
      
      $selectSearchHitsButton.click(function(){
        _.each(searchedModels, function(model){
          model.set(chKey,true, {silent: true});
        });
        if (!_.isEmpty(searchedModels)){
          self.collection.trigger('bulkchange',searchedModels);
        }
      });
      $unselectSearchHitsButton.click(function(){
        _.each(searchedModels, function(model){
          model.set(chKey,false, {silent: true});
        });
        if (!_.isEmpty(searchedModels)){
          self.collection.trigger('bulkchange', searchedModels);
        }
      });

      function treeChanged(e){
        var id = $(e.target).parent().attr('id');
        var model = self.collection.get(id);
        if (model){
          if ($(e.target).prop(chKey) === true){
            model.set(chKey,true);
          } else {
            model.set(chKey, false);
          }
        }
      }
      $treeControl.on('change', treeChanged);
      $selectedTreeControl.on('change', treeChanged);

      $unselectSelectedButton.click(function(){
        var changedModels = self.collection.where({checked: true});
        _.each(changedModels, function(model){
          model.set(chKey, false, { silent: true });
        });
        self.collection.trigger('bulkchange', changedModels);
      });
      
      self.collection.on('change:checked', function(model){
        var id = self.collection.modelId(model)
        var li = $treeControl.find('#'+ id)
        var selectedLi = $selectedTreeControl.find('#'+ id)
        console.log('collection change', id, li, model);
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
        showChecked($selectedTreeControl);
      });
      
      self.collection.on('bulkchange', function(changedModels){
        console.log('collection change', changedModels);
        var storedParents1 = {};
        var storedParents2 = {};
        _.each(changedModels, function(model){
          var id = self.collection.modelId(model)
          var li = $treeControl.find('#'+ id)
          var selectedLi = $selectedTreeControl.find('#'+ id)
          if (li){
            if (model.get(chKey)===true){
              li.find('input').prop(chKey, true);
              // NOTE: defer notification for performance
              //li.find('input').prop(chKey, true).trigger('change');
            } else {
              li.find('input').prop(chKey, false);
            }
            storedParents1[li.parent().parent().attr('data-name')] = li;
          }
          if (selectedLi){
            if (model.get(chKey)===true){
              selectedLi.find('input').prop(chKey, true);
            } else {
              selectedLi.find('input').prop(chKey, false);
            }
            storedParents2[selectedLi.parent().parent().attr('data-name')] = selectedLi;
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
      });
      
      $clearSearchButton.click(function(){
        $search.val('');
        $treeControl.find('li').show();
        bonsai.collapseAll();
        _.each(self.collection.where({checked: true}), function(model){
          var id = self.collection.modelId(model);
          bonsai.expandTo(id);
        });
        searchedModels.length = 0;
        $selectSearchHitsButton.hide();
        $unselectSearchHitsButton.hide();
      });

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