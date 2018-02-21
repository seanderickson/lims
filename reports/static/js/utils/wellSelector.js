define([
  'jquery',
  'underscore',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'templates/genericResource.html'
], 
function($, _, Backgrid, Iccbl, appModel, genericLayout) {
  /**
   * WellSelector is a UI Grid Viewer that uses plate mapping row and column
   * assignment.
   * - If editable, grid squares (wells) may be selected by the user.
   */
  /** Clicking this cell toggles the entire grid **/
  var ToggleColumnHeader = Backgrid.HeaderCell.extend({
    
    /** @property */
    events: {
      "click": "cellClicked",
    },
    
    /** Clicking this cell toggles the entire grid **/
    cellClicked: function(e) {
      if (this.collection.editable === false) return;
      this.collection.each(function(model){
        model.set(_.object(_.map(
          _.keys(model.omit('row')),
          function(key){
            return [key, !model.get(key)];
          }
        )));
      });
    },
    
    render: function(){
      ToggleColumnHeader.__super__.render.apply(this);
      this.$el.attr('title', 'Toggle the grid values');
      return this;
    }
  });
  
  /** Clicking this cell sets/unsets the entire column **/
  var WellColumnSelectorHeader = Backgrid.HeaderCell.extend({
    
    className: 'wellselector',
    
    /** @property */
    events: {
      "click": "cellClicked",
      "click input[type=checkbox]": "cellClicked"
    },

    initialize : function(options) {
    
      var self = this;
      WellColumnSelectorHeader.__super__.initialize.apply(this, arguments);
      var colNumber = parseInt(this.column.get('name'));
      
      this.listenTo(this.collection, "select-column", function (colNumber, selected) {
        var colNumber = parseInt(colNumber);
        self.collection.each(function(model){
          model.set(colNumber, selected);
        });
      });
      
      this.listenTo(this.collection, 'change:' + colNumber, function(){
        self.renderSelection();
      });
    },
    
    /**
     * @return true if the column cells are all selected.
     */
    isColumnSelected: function(){
      var self = this;
      var colNumber = parseInt(this.column.get('name'));
      var hasFalse = self.collection.find(function(model){
        return model.get(colNumber)===false;
      });
      return _.isUndefined(hasFalse);
    },
    
    cellClicked: function(e) {
      if (this.collection.editable === false) return;
      var checked = this.checkbox().prop("checked");
      checked = !checked;
      this.checkbox().prop("checked", checked);
      
      // Because the col does not have a state model, fire a select-column event
      this.collection.trigger("select-column", this.column.get('name'), checked);
    },

    checkbox: function () {
      return this.$el.find("input[type=checkbox]");
    },
      
    renderSelection: function(){
      var self = this;
      var selected = self.isColumnSelected();
      self.checkbox().prop('checked', selected);
      this.$el.toggleClass("selected", selected);
    },
    
    render : function() {
      var self = this;
      WellColumnSelectorHeader.__super__.render.apply(this);
      var isSelected = self.isColumnSelected();
      var columnSelector = $('<input type="checkbox" style="display: none;" >');
      this.$el.prepend(columnSelector);
      self.renderSelection();
      return this;
    }
  });
  
  
  /** Clicking this cell sets/unsets the entire row **/
  var RowSelectorCell = Backgrid.Cell.extend({
  
    /** @property */
    className: "wellselector select-row-cell",
  
    /** @property */
    tagName: "td",
  
    /** @property */
    events: {
      "keydown input[type=checkbox]": "onKeydown",
      "change input[type=checkbox]": "onChange",
      "click": "cellClicked",
      "click input[type=checkbox]": "cellClicked"
    },
    
    
    initialize: function(options){
      var self = this;
      this.column = options.column;
      if (!(this.column instanceof Backgrid.Column)) {
        this.column = new Backgrid.Column(this.column);
      }
      
      // TODO: review: code copied from Backgrid.Cell
      // Related to keyboard navigation and displaying editmode for cells
      var column = this.column, model = this.model, $el = this.$el;
      this.listenTo(column, "change:renderable", function (column, renderable) {
        $el.toggleClass("renderable", renderable);
      });
  
      if (Backgrid.callByNeed(column.renderable(), column, model)) $el.addClass("renderable");
      ///////
      
      this.listenTo(model, "select-row", function (model, selected) {
        if (selected) {
          console.log('select row');
          self.model.set(_.object(_.map(
            _.keys(self.model.omit('row')),
            function(key){
              return [key, true];
            }
          )));
        } else {
          self.model.set(_.object(_.map(
            _.keys(self.model.omit('row')),
            function(key){
              return [key, false];
            }
          )));
        }
      });
      
      this.listenTo(model, 'change', function(model, options){
        self.renderSelection();
      });
      _.bindAll(this, 'isRowSelected','renderSelection');
    },
    
    renderSelection: function(){
      var isSelected = this.isRowSelected();
      this.checkbox().prop('checked',isSelected);
      this.$el.toggleClass("selected", isSelected);
    },

    /**
     * @return true if all the cells in the row are selected
     */
    isRowSelected: function() {
      var self = this;
      var hasFalse = _.find(self.model.keys(), function(key){
        if (key=='row') return false;
        if (self.model.get(key)===false) return true;
      });
      return _.isUndefined(hasFalse);
    },
    
    /**
       Returns the checkbox.
     */
    checkbox: function () {
      return this.$el.find("input[type=checkbox]");
    },
  
    /**
       Focuses the checkbox.
    */
    enterEditMode: function () {
      this.checkbox().focus();
    },
  
    /**
       Unfocuses the checkbox.
    */
    exitEditMode: function () {
      this.checkbox().blur();
    },
  
    /**
       Process keyboard navigation.
       ** FIXME: copied code from parent class - remove.
    */
    onKeydown: function (e) {
      var command = new Backgrid.Command(e);
      if (command.passThru()) return true; // skip ahead to `change`
      if (command.cancel()) {
        e.stopPropagation();
        this.checkbox().blur();
      }
      else if (command.save() || command.moveLeft() || command.moveRight() ||
               command.moveUp() || command.moveDown()) {
        e.preventDefault();
        e.stopPropagation();
        this.model.trigger("backgrid:edited", this.model, this.column, command);
      }
    },
  
    cellClicked: function(e) {
      if (this.model.collection.editable === false) return;
      var checked = this.checkbox().prop("checked");
      checked = !checked;
      this.checkbox().prop("checked", checked);
      
      // Because the row does not have a state model, fire a select-row event
      this.model.trigger("select-row", this.model, checked);
    },
    
    /**
       Renders a checkbox in a table cell.
    */
    render: function () {

      this.$el.empty();
      var rowSelector = $('<input tabindex="-1" type="checkbox" style="display: none;" />');
      this.$el.append(rowSelector);
      this.$el.append(Iccbl.rowToLetter(
        parseInt(this.model.get('row'))));
      this.renderSelection();

      this.delegateEvents();
      return this;
    }
  });
  
  var WellSelectorCell = Backgrid.Cell.extend({

    /** @property */
    tagName: "td",
  
    /** @property */
    className: "wellselector",

    /** @property */
    events: {
      'click': "cellClicked"
    },
    
    initialize: function(options){

      var self = this;
      this.column = options.column;
      
      if (!(this.column instanceof Backgrid.Column)) {
        this.column = new Backgrid.Column(this.column);
      }
      
      // TODO: review: actions copied from Backgrid.Cell initializer
      var column = this.column, model = this.model, $el = this.$el;
      this.listenTo(column, "change:renderable", function (column, renderable) {
        $el.toggleClass("renderable", renderable);
      });
  
      if (Backgrid.callByNeed(column.renderable(), column, model)) $el.addClass("renderable");

      this.listenTo(model, "change:" + column.get("name"), function () {
        self.renderSelection();
      });
      
      _.bindAll(this,'cellClicked', 'getValue');
    },
    
    cellClicked: function(e) {
      if (this.model.collection.editable === false) return;
      var self = this;
      // Toggle
      var columnNumber = parseInt(self.column.get("name"));
      self.model.set(columnNumber, !self.getValue());
    },
    
    getValue: function(){
      var self = this;
      var columnNumber = parseInt(self.column.get("name"));
      return self.model.get(columnNumber);
    },
    
    renderSelection: function() {
      this.$el.toggleClass("selected", this.getValue());
    },
    
    render: function () {
      
      var model = this.model;
      var row = this.model.get('row');
      var columnNumber = parseInt(this.column.get("name"));
      var val = this.model.get(columnNumber);
      
      if (columnNumber < 10) columnNumber = '0' + columnNumber;

      this.$el.empty();
      this.$el.append(Iccbl.rowToLetter(row) + columnNumber);
      this.renderSelection();
      return this;
    }
    
  });
  
  var WellSelector = Backbone.Layout.extend({
    
    template: _.template(genericLayout),

    initialize: function(options) {
      var self = this;
      console.log('initialize wellselector');
      
      var options = options || {};
      self.options = options;
      
      var plate_size = options.plateSize || 384;
      var nCols = 24;
      var nRows = 16;
      if (plate_size === 96){
        nCols = 12;
        nRows = 8;
      }
      self.nCols = nCols;
      self.nRows = nRows;
      
      var wellSelections = options.wellSelections || '';

      var RowModel = Backbone.Model.extend({
        idAttribute: 'row'
      })
      var RowCollection = Backbone.Collection.extend({
        model: RowModel
      });
      
      var rowCollection = this.rowCollection = new RowCollection();
      // Hack, store the editable flag here
      rowCollection.editable = self.options.editable;
      // Initialize all cells to false
      for (var i=0; i<nRows; i++){
        var row = { row: i };
        for (var j=1; j<=nCols; j++){
          row[j] = false;
        }
        rowCollection.add(new RowModel(row));
      };
      
      console.log('wellSelections', wellSelections);
      if (!_.isEmpty(wellSelections)){
        _.each(wellSelections.split(/\s*,\s*/), function(wellSelection){
          var colMatch = wellSelection.match(/col:\s*(\d{1,2})/i);
          if (colMatch !== null){
            var col = parseInt(colMatch[1]);
            if (col > nCols){
              throw 'Selection is out of the column range: ' + wellSelection;
            }
            rowCollection.each(function(rowModel){
              rowModel.set(col,true);
            });
            return;
          }
          var rowMatch = wellSelection.match(/row:\s*([a-zA-Z]{1,2})/i);
          if (rowMatch !== null){
            rowModel = rowCollection.find(function(model){
              if (Iccbl.rowToLetter(model.get('row'))==rowMatch[1].toUpperCase()){
                return true;
              }
              return false;
            });
            if (!_.isUndefined(rowModel)){
              _.each(rowModel.keys(), function(key){
                if (key!=='row'){
                  rowModel.set(key,true);
                }
              });
            } else {
              throw 'Selection is out of the row range: ' + wellSelection;
            }
            return;
          }
          
          var wellMatch = wellSelection.match(/([a-zA-Z]{1,2})(\d{1,2})/);
          if (wellMatch != null){
            var row = wellMatch[1].toUpperCase();
            var col = parseInt(wellMatch[2]);
            if (col > nCols){
              throw 'Selection is out of the column range: ' + wellSelection;
            }
            rowModel = rowCollection.find(function(model){
              if (Iccbl.rowToLetter(model.get('row'))==row){
                return true;
              }
              return false;
            });
            if (!_.isUndefined(rowModel)){
              rowModel.set(col, true);
            }
            return;
          }
          
          throw new Exception('Unknown well selection pattern: "' + wellSelection + '"');
        });
      }
      
      console.log('initialized...', this.rowCollection);
    },
    
    afterRender: function(){
      var self = this;
      console.log('afterRender grid...', self.rowCollection);
      var colTemplate = {
        'cell' : WellSelectorCell,
        'order' : -1,
        'sortable': false,
        'searchable': false,
        'editable' : false,
        'visible': true,
      };
      var columns = [
        _.extend({},colTemplate,{
          'name' : 'row',
          'label' : '',
          'description' : 'Select Row',
          'order': 1,
          'sortable': false,
          'editable': false,
          'cell': RowSelectorCell,
          'headerCell': ToggleColumnHeader
          
        }),
      ];
      for (var i=1; i<=self.nCols; i++){
        columns.push(
          _.extend({},colTemplate,{
            'name' : i,
            'label' : i,
            'description' : 'Column',
            'order': 1,
            'sortable': false,
            'cell': WellSelectorCell,
            'headerCell': WellColumnSelectorHeader
          })
        );
      }
      var colModel = new Backgrid.Columns(columns);
      var _grid = new Backgrid.Grid({
        className: 'wellselector ', 
        columns: colModel,
        collection: self.rowCollection
      });
      _grid.render();          
      
      var clearSelectionsButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="save_button_selected" href="#">',
          'Clear</a>'
        ].join(''));
      clearSelectionsButton.click(function(e){
        e.preventDefault();
        self.rowCollection.each(function(model){
          model.set(_.object(_.map(
            _.keys(model.omit('row')),
            function(key){
              return [key, false];
            }
          )));
        });
      });
      var selectAllButton = $([
        '<a class="btn btn-default btn-sm pull-down" ',
          'role="button" id="save_button_selected" href="#">',
          'Select all</a>'
        ].join(''));
      selectAllButton.click(function(e){
        e.preventDefault();
        self.rowCollection.each(function(model){
          model.set(_.object(_.map(
            _.keys(model.omit('row')),
            function(key){
              return [key, true];
            }
          )));
        });
      });
      
      $('#resource_content').html(_grid.el);
      if (self.options.editable){
        $('#resource_content').prepend(selectAllButton);
        $('#resource_content').prepend(clearSelectionsButton);
      }
      if (self.options.title){
        $('#content_title').html(self.options.title);
        $('#content_title_row').show();
      } else {
        $('#content_title_row').hide();
      }
    },
    
    getSelectedWells: function(){
      var self = this;
      var fullCols = [];
      var fullRows = [];
      var individualWells = [];
      
      self.rowCollection.each(function(rowModel){
        var row = Iccbl.rowToLetter(rowModel.get('row'));
        var isFalse = _.find(rowModel.keys(), function(key){
          return key != 'row' && rowModel.get(key)==false;
        });
        if (_.isUndefined(isFalse)) fullRows.push(row);
      });
      for(var i=1; i<=self.nCols; i++){
        var isFalse = self.rowCollection.find(function(rowModel){
          return rowModel.get(i)==false;
        });
        if (_.isUndefined(isFalse)) fullCols.push(i);
      }
      self.rowCollection.each(function(rowModel){
        var row = Iccbl.rowToLetter(rowModel.get('row'));
        if (!_.contains(fullRows,row)){
          _.each(rowModel.keys(), function(key){
            if (key!=='row' && rowModel.get(key)){
              var columnNumber = parseInt(key);
              if (!_.contains(fullCols, columnNumber)){
                if (columnNumber < 10) columnNumber = '0' + columnNumber;
                var wellName = row + columnNumber;
                individualWells.push(wellName);
              }
            }
          });
        }
      });
      fullRows = _.map(fullRows, function(row){ return 'Row:'+row; });
      fullCols = _.map(fullCols, function(col){ return 'Col:'+col });
      return fullCols.concat(fullRows,individualWells);
    }
    
  });
  
  return WellSelector;
  
});