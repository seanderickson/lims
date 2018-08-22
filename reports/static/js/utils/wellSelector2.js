define([
  'jquery',
  'underscore',
  'backgrid',
  'iccbl_backgrid',
  'models/app_state',
  'templates/genericResource.html',
  'google_palette'
], 
function($, _, Backgrid, Iccbl, appModel, genericLayout, palette) {
  /**
   * WellSelector is a UI Grid Viewer that uses plate mapping row and column
   * assignment.
   * - If editable, grid squares (wells) may be selected by the user.
   * 
   * V2: NamedRanges color support:
   * - palette may be specified, otherwise it will be generated:
   * see: https://github.com/google/palette.js/tree/master
   * 
   * TODO: wellSelector.js should be replaced by this version.
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
      this.collection.trigger('selection_update');
      
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

      var columnIndex = this.columnIndex = parseInt(this.column.get('name'))-1;
      
      this.on("select-column", function (colNumber, selected) {
        self.collection.each(function(model){
          model.set(columnIndex, selected);
        });
        self.collection.trigger('selection_update');
      });
      
      this.listenTo(this.collection, 'change:' + columnIndex, function(){
        console.log('change:(column index):' + columnIndex);
        self.renderSelection();
      });
    },
    
    /**
     * @return true if the column cells are all selected.
     */
    isColumnSelected: function(){
      var self = this;
      var hasFalse = self.collection.find(function(model){
        return model.get(self.columnIndex)===false;
      });
      return _.isUndefined(hasFalse);
    },
    
    cellClicked: function(e) {
      if (this.collection.editable === false) return;
      var checked = this.checkbox().prop("checked");
      checked = !checked;
      this.checkbox().prop("checked", checked);
      
      // Because the col does not have a state model, fire a select-column event
      this.trigger("select-column", this.column.get('name'), checked);
    },

    checkbox: function () {
      return this.$el.find("input[type=checkbox]");
    },
      
    renderSelection: function(){
      var selected = this.isColumnSelected();
      this.checkbox().prop('checked', selected);
      this.$el.toggleClass("selected", selected);
    },
    
    render : function() {
      WellColumnSelectorHeader.__super__.render.apply(this);
      var columnSelector = $('<input type="checkbox" style="display: none;" >');
      this.$el.prepend(columnSelector);
      this.renderSelection();
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
        self.model.collection.trigger('selection_update');
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
      var columnIndex = this.columnIndex = parseInt(this.column.get('name'))-1;
      var rowIndex = this.rowIndex = this.model.get('row');
      var wellName = this.wellName = Iccbl.getWellName(this.rowIndex,this.columnIndex);

      var column = this.column, model = this.model, $el = this.$el;
      this.listenTo(column, "change:renderable", function (column, renderable) {
        $el.toggleClass("renderable", renderable);
      });
  
      if (Backgrid.callByNeed(column.renderable(), column, model)) $el.addClass("renderable");

      this.listenTo(model, "change:" + columnIndex, function () {
        self.renderSelection();
      });
      
      _.bindAll(this,'cellClicked', 'getValue');
    },
    
    cellClicked: function(e) {
      // Toggle
      if (this.model.collection.editable === false) return;
      this.model.set(this.columnIndex, !this.getValue());
      this.model.collection.trigger('selection_update');
    },
    
    getValue: function(){
      return this.model.get(this.columnIndex);
    },
    
    renderSelection: function() {
      this.$el.toggleClass("selected", this.getValue());
    },
    
    render: function () {
      this.$el.empty();
      this.$el.append(this.wellName);
      this.renderSelection();
      return this;
    }
    
  });
  
  var NamedRange = Backbone.Model.extend({
    idAttribute: "ordinal"
  });
  
  var NamedRangeView = Backbone.View.extend({
    _template: [
        '<span class="input-group-addon" id="basic-addon1">',
        '<input type="radio" name="rangeControl"  id="radio_{ordinal}" >',
        '</span>',
        '<textarea rows=1 class="form-control textarea-noscroll" id="text_{ordinal}" readonly ',
        '  placeholder="Enter a label" value="{label}" aria-describedby="basic-addon1">',
        '{label}</textarea>',
        '<a class="btn input-group-addon disabled" id="clear_{ordinal}"  >x</a>',
      ].join(''),
    
    events: {
      'click input[type=radio]' : 'selectRange',
      'click a': 'clearSelected',
      'change input[type=text]' : 'changeLabel',
      'click textarea' : 'selectRange',
      'change textarea' : 'changeLabel',
    },
    
    clearSelected: function(e){
      e.preventDefault();
      if ($(e.target).hasClass('disabled')){
        //pass
      }else{
        console.log('clearSelected', e.target);
        this.model.set('wells', []);
      }
    },
    
    selectRange: function(){
      console.log('selected-range', this.model);
      this.model.trigger('selected-range',this.model);
    },
    
    changeLabel: function(e){
      var text = this.$('textarea').val();
      console.log('changeLabel', text);
      
      var extantModel = this.model.collection.findWhere({
        label: text
      });
      
      if (extantModel){
        this.model.trigger('set_errors', ['Label already used: "' + text + '"']);
        this.$('textarea').val('');
        return;
      }else{
        this.model.trigger('set_errors', null);
        this.model.set('label', text);
        console.log('changeLabel', this.model);
      }
    },
    
    getColor: function(){
      console.log('TODO: getColor must be defined');
      return null;
    },
    
    render: function() {
      this.$el.html(Iccbl.formatString(this._template, this.model));
      this.$el.addClass('input-group  col-sm-12');
      var color = this.getColor();
      this.$el.css('background-color', color);
      this.$el.find('input').css('background-color', color);
      this.$el.children().css('background-color', color);
      return this;
    }    
    
  });

  var NamedRangesView = Backbone.View.extend({
    
    allRangesEntry: [
        '<div class="input-group col-sm-12" id="show_all_input_group" >',
        '<span class="input-group-addon" id="basic-addon1">',
        '<input type="radio" name="rangeControl"  id="radio_show_all" checked >',
        '</span>',
        //'<input type="text" class="form-control" id="text_0" disabled="disabled" ',
        //'  value="All" >',
        //'</input>',
        '<textarea rows=1 class="form-control" id="text_0" readonly >',
        'All</textarea>',
        '<a class="btn input-group-addon" id="add_range_button" >+</a>',
        '</div>'
      ].join(''),
    
    events: {
      'click #add_range_button': 'addRange',
      'change #radio_show_all': 'showAll',
      'click #text_0' : 'showAll'
    },

    initialize: function(options){
      var self = this;
      this.collection = options.collection || new Backbone.Collection();
      this.editable = options.editable;
      
      _.bindAll(this,'addRange','showAll','selectRange','generate_palette');
      this.collection.on('selected-range', self.selectRange);
      this.collection.on('add remove',self.generate_palette);
      
      if (options.colorPalette){
        this.palette = options.colorPalette;
      } else {
        this.generate_palette();
      }
    },
    
    generate_palette: function(){
      
      var maxRange = 1;
      if (!this.collection.isEmpty()){
        maxRange = this.collection.max(function(range){ return range.id }).id;
      }
      // see: http://google.github.io/palette.js/ for palettes chosen
      if (maxRange < 11){
        this.palette = palette('cb-Accent', maxRange );
      } else {
        this.palette = palette('tol-rainbow', maxRange);
      }
      console.log('generate_palette', maxRange, this.palette);
    },
    
    get_palette: function(){
      return this.palette;
    },
    
    getColor: function(range){
      var self = this;
      var palette = self.get_palette();
      if (range.id > palette.length){
        self.generate_palette();
      }
      
      var color = self.palette[range.id-1];
      if (color.indexOf('#')<0){
        color = '#' + color;
      }
      // console.log('getColor:', range.id, range.get('label'), color);
      return color;
    },
    
    selectRange: function(namedRangeModel){
      var self = this;
      console.log('model selected-range', namedRangeModel);

      self.collection.each(function(model){
        model.set('selected', false);
      });
      namedRangeModel.set('selected', true);
      self.$el.find('textarea').prop('readonly', true);
      self.$el.find('a').addClass('disabled');
      self.$el.find('#radio_'+namedRangeModel.id).prop('checked', true);
      
      if (self.editable){
        self.$el.find('#text_'+namedRangeModel.id).prop('readonly', false);
        self.$el.find('#clear_'+namedRangeModel.id).removeClass('disabled');
      }
      self.collection.trigger('selection_update');
    },
    
    addRange: function(){
      console.log('add range');
      var newRange = new NamedRange({
        ordinal: this.collection.size()+1,
        wells: [],
        label: ''
      });
      console.log('newRange.id', newRange.id)
      this.collection.add(newRange);
      this.render();
      this.selectRange(newRange);
    },
    
    showAll: function(){
      var self = this;
      console.log('showAll');
      var selectedRange = this.collection.findWhere({'selected': true});
      if (selectedRange){
        selectedRange.set('selected',false);
        self.$el.find('textarea').prop('readonly', true);
        self.$el.find('a').addClass('disabled');
        self.$el.find('#radio_show_all').prop('checked',true);
        self.collection.trigger('selection_update');
      }
      if (self.editable){
        self.$el.find('#show_all_input_group').find('a').removeClass('disabled');
      }
    },
    
    render: function(){
      var self = this;
      this.$el.empty();
      this.$el.addClass('input-group col-sm-12');
      
      this.$el.append(this.allRangesEntry);
      
      var PaletteNamedRangeView = NamedRangeView.extend({
        getColor: function(){
          return self.getColor(this.model);
        }
      });
      
      this.collection.each(function(namedRangeModel){
        var namedRangeView = new PaletteNamedRangeView({
          model: namedRangeModel
        });
        self.$el.append(namedRangeView.render().el);
      });
      return this;
    }
  });
  
  var WellSelector = Backbone.Layout.extend({

    template: $([
      '<div id="content_title_row" class="row" style="display: none;">',
      '  <div class="panel col-sm-12" style="margin-bottom: 0px;">',
      '    <h4 id="content_title"></h4>',
      '  </div>',
      '</div>',
      '<div class="row" id="errors_div" class="row" style="display: none;"></div>',
      '<div class="row">',
      '<div class="col-sm-3" id="buttons_div_row" >',
      '<div class="row" id="buttons_div" ></div>',
      '</div>',
      '<div class="col-sm-9" id="grid_div_row" >',
      '<div class="row" id="grid_div"></div>',
      '</div>',
      '</div>'
    ].join('')),
    
    events: {
      'click #clear_button_selected': 'clearSelected'
    },

    initialize: function(options) {
      var self = this;
      console.log('initialize wellselector');
      
      var options = options || {
        useNamedRanges: false
      };
      self.options = options;
      var plate_size = self.plate_size = options.plateSize || 384;
      var nCols = self.nCols = Iccbl.getCols(plate_size);
      var nRows = self.nRows = Iccbl.getRows(plate_size);
      var useNamedRanges = self.useNamedRanges = options.useNamedRanges;
      var namedWellRangeInput = options.namedWellRangeInput || {};
      var errors = options.errors || [];

      console.log('namedWellRangeInput', namedWellRangeInput);

      var RowModel = Backbone.Model.extend({
        idAttribute: 'row'
      })
      var RowCollection = Backbone.Collection.extend({
        model: RowModel
      });
      var rowCollection = this.rowCollection = new RowCollection();

      // Initialize the grid to false (not selected)
      for (var i=0; i<nRows; i++){
        var row = { row: i };
        for (var j=0; j<nCols; j++){
          row[j] = false; 
        }
        rowCollection.add(new RowModel(row));
      };
      
//      if (!options.useNamedRanges != true){
//        if (_.size(namedWellRange)>0){
//          var testRange = this.namedWellRanges.at(0);
//          if (!_.isEmpty(testRange)){
//            appModel.showError
//          }
//        }
//      }
      
//      var errors = [];
//      var parsedRanges = Iccbl.parseNamedWellSelections(wellSelections, plate_size, errors);
      var namedWellRanges = this.namedWellRanges = 
        new Backbone.Collection(_.map(_.values(namedWellRangeInput), function(nr){
          return new NamedRange(nr);
        }),
        {comparator: 'ordinal' });

      console.log('initialize: namedWellRanges', namedWellRanges);
      if (!_.isEmpty(errors)){
        this.setErrors(errors);
      }
      
      self.setSelectedWellsFromRanges();
      
      this._createGrid();
      
      var namedRangesEditable = self.options.editable;
      if (options.useNamedRanges == true){
        this.namedRangesView = new NamedRangesView({ 
          collection: this.namedWellRanges,
          editable: self.options.editable,
          colorPalette: self.options.colorPalette
        });
      }
      
      _.bindAll(this,'setSelectedWellsFromRanges','setSelectedWells','_clearGridOnly',
        'recordGridSelections','setErrors','getErrors');

      rowCollection.on('selection_update', self.recordGridSelections);
      namedWellRanges.on('selection_update', self.setSelectedWellsFromRanges);
      namedWellRanges.on('change:wells', self.setSelectedWellsFromRanges);
      namedWellRanges.on('set_errors', self.setErrors);

      console.log('initialized...', this.rowCollection, this.namedWellRanges);
    },

    /**
     * Record the selected wells from the displayed grid to the selected
     * namedWellRange model.
     */
    recordGridSelections: function(){
      var self = this;
      console.log('recordGridSelections' );
      var newWells = [];
      self.rowCollection.each(function(rowModel){
        var rowIndex = rowModel.get('row');
        _.each(rowModel.keys(), function(key){
          if (key!=='row' && rowModel.get(key)){
            var columnIndex = parseInt(key);
            newWells.push(Iccbl.getWellName(rowIndex,columnIndex));
          }
        });
      });
      var selectedRange = this.namedWellRanges.at(0);
      if (this.namedWellRanges.size()>1){
        var selectedRange = this.namedWellRanges.findWhere({'selected': true});
      }
      newWells.sort();
      
      allWells = _.map(self.namedWellRanges.without(selectedRange), function(namedRange){
        return namedRange.get('wells');
      });
      allWells.push(newWells);
      var duplicates = Iccbl.find_duplicates(allWells);
      if (!_.isEmpty(duplicates)){
        self.setErrors(['wells are shared between ranges: ' + duplicates.join(', ')]);
      }else{
        self.setErrors(null);
        selectedRange.set({'wells': newWells}, {silent: true });
      }
    },
    
    /**
     * Clear the grid, without triggering updates to the named ranges.
     */
    _clearGridOnly: function(){
      this.rowCollection.each(function(model){
        model.set(_.object(_.map(
          _.keys(model.omit('row')),
          function(key){
            return [key, false];
          }
        )));
      });
    },
    
    /**
     * Set the wells as selected in the grid:
     * - does not deselect wells already selected
     */
    setSelectedWells: function(wells){
      var self = this;
      var rowCollection = this.rowCollection;
      _.each(wells, function(well){
        var row_col = Iccbl.getWellRowCol(well);
        var row = row_col[0];
        var col = row_col[1];
        
        rowModel = rowCollection.find(function(model){
          if (model.get('row')==row){
            model.set(col, true);
            return true;
          }
        });
        if (_.isUndefined(rowModel)){
          console.log('error, no well found', well);
        }
      });
    },
    
    /**
     * Set the well selections in the displayed grid from the currently
     * selected NamedWellRange, or for all ranges if no range is selected.
     */
    setSelectedWellsFromRanges: function() {
      console.log('setSelectedWellsFromRanges...')
      var self = this;
      var namedWellRanges = this.namedWellRanges;
      self._clearGridOnly();
      var selectedRange = namedWellRanges.findWhere({'selected': true});
      if (selectedRange){
        console.log('setSelectedWells:', selectedRange.get('label'));
        self.setSelectedWells(selectedRange.get('wells'));
        // grid is editable only if a named range is selected
        if (self.options.editable){
           self.rowCollection.editable = true;
        }
      }else{
        // grid is editable only if a named range is selected
        self.rowCollection.editable = false;

        namedWellRanges.each(function(namedWellRange){
          self.setSelectedWells(namedWellRange.get('wells'));
        });
      }
    },
    
    /**
     * Set and display errors for a selection operation
     */
    setErrors: function(errors){
      var errorsDiv = this.$el.find('#errors_div');
      errorsDiv.empty();
      if(!errors || _.isEmpty(errors)){
        errorsDiv.hide();
      } else {
        console.log('setErrors', errors);
        errorsDiv.show();
        _.each(errors, function(error){
          errorsDiv.append('<div class="alert alert-danger">' + error + '</div>');
        });
      }
      this.errors = errors;
    },
    
    getErrors: function(){
      return this.errors;
    },
    
    _createGrid: function(){
      var self = this;
      var WellSelectorCellFinal = WellSelectorCell;
      
      if(self.useNamedRanges){
        function colorCalculator(wellName){
          var selectedRange = self.namedWellRanges.findWhere({'selected': true});
          if (selectedRange){
            return self.namedRangesView.getColor(selectedRange);
          }
          var range = self.namedWellRanges.find(function(model){
            return _.contains(model.get('wells'), wellName);
          });
          if (range){
            return self.namedRangesView.getColor(range);
          }
          return '';
        };
        WellSelectorCellFinal = WellSelectorCell.extend({
          
          renderSelection: function() {
            if (this.getValue()){
              var row = this.model.get('row');
              var columnIndex = parseInt(this.column.get("name"))-1;
              var wellName = Iccbl.getWellName(row,columnIndex);
              this.$el.css('background-color', colorCalculator(wellName));
            }else{
              this.$el.removeClass('selected');
              this.$el.css('background-color', '');
            }
          },
        });
      }
      var colTemplate = {
        'cell' : WellSelectorCellFinal,
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
            'headerCell': WellColumnSelectorHeader
          })
        );
      }
      var className = 'wellselector';
      if (self.plateSize == 1536){
        className += ' ' + 'extra-small-table';
      }
        
      var colModel = new Backgrid.Columns(columns);
      var _grid = this._grid = new Backgrid.Grid({
        className: className, 
        columns: colModel,
        collection: self.rowCollection
      });
    },
    
    afterRender: function(){
      var self = this;
      console.log('afterRender grid...', self.rowCollection);
      
      this.$el.empty();
      this.$el.addClass('container-fluid');
      this.$el.append(this.template);
      
      var buttonDiv = this.$el.find('#buttons_div').empty();

      this._grid.render();          
      this.$el.find('#grid_div').empty().append(this._grid.el);

      if (this.options.useNamedRanges == true){
        buttonDiv.append(this.namedRangesView.render().el);
      } else {
        // TODO: add 'clear' and 'select-all' buttons
      }
      
      // call set_errors, to show any parse errors recorded
      self.setErrors(self.getErrors());
      
      if (self.options.title){
        $('#content_title').html(self.options.title);
        $('#content_title_row').show();
      } else {
        $('#content_title_row').hide();
      }
      
      if (this.plateSize == 1536){
        $('#buttons_div_row').removeClass();
        $('#grid_div_row').removeClass();
        $('#buttons_div_row').addClass('col-sm-2');
        $('#grid_div_row').addClass('col-sm-10');
      }else if (this.plateSize == 96){
        $('#buttons_div_row').removeClass();
        $('#grid_div_row').removeClass();
        $('#buttons_div_row').addClass('col-sm-4');
        $('#grid_div_row').addClass('col-sm-8');
      }else{
        $('#buttons_div_row').removeClass();
        $('#grid_div_row').removeClass();
        $('#buttons_div_row').addClass('col-sm-3');
        $('#grid_div_row').addClass('col-sm-9');
      }
      
    },
    
    getSelectedWells: function(){
      return this.namedWellRanges.toJSON();
    },
    
  });
  
  return WellSelector;
  
});