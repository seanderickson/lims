define(['jquery', 'underscore', 'backbone', 'backgrid','backbone_forms', 
        'backgrid_filter', 'backgrid_paginator', 'backgrid_select_all',
        'layoutmanager'],
    function($, _, Backbone, Backgrid, BackboneForms,
             BackgridFilter, BackgridSelectAll, layoutmanager ) {

var root = window;
  
var Iccbl = root.Iccbl = {
    VERSION : "0.0.1",
    appModel : "This value will be initialized on app start"
};
var ICCBL_DATE_RE = Iccbl.ICCBL_DATE_RE =  /^(\d{1,2})\/(\d{1,2})\/(\d{2,4})$/;
var DATE_RE = Iccbl.DATE_RE = /^([+\-]?\d{4})-(\d{2})-(\d{2})$/;
var TIME_RE = Iccbl.TIME_RE = /^(\d{2}):(\d{2}):(\d{2})(\.(\d{3}))?$/;
var ISO_SPLITTER_RE = Iccbl.ISO_SPLITTER_RE = /T|Z| +/;


/**
 * Convert a plate row index to a letter
 */
var rowToLetter = Iccbl.rowToLetter = function(i){
  if (i<26){
    return 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'[i];
  } else {
    var rem = i%26;
    var part = parseInt(i/26)-1;
    return rowToLetter(part) + rowToLetter(rem);
  }
};


/**
 * Utility function from Backgrid
 */
function lpad(str, length, padstr) {
  var paddingLen = length - (str + '').length;
  paddingLen =  paddingLen < 0 ? 0 : paddingLen;
  var padding = '';
  for (var i = 0; i < paddingLen; i++) {
    padding = padding + padstr;
  }
  return padding + str;
}

/**
 * Format a string with embedded replacement fields.
 * 
 * "replacment fields" are surrounded by braces '{}'. 
 * Each replacement field is used as a key to lookup values in the object.
 * 
 * @param object - either a Backbone.Model, or a object
 * @param defaul_val value to use if the matched token is not found in the model
 * - this can be used to replace any token with a given default value
 * - if default_val is not provided, the replacement_field is left in the string.
 */
var formatString = Iccbl.formatString = function(
    stringWithTokens, 
    object, 
    default_val) 
  {
  var isBackboneModel = object instanceof Backbone.Model;
  var interpolatedString = stringWithTokens.replace(/{([^}]+)}/g, 
    function (match) 
    {
      match = match.replace(/[{}]/g,'');
      if(isBackboneModel && !_.isUndefined(object.get(match))){
        return object.get(match);
      }else if(_.has(object, match)){
        return object[match];
      }else{
        if(!_.isUndefined(default_val)){
          return default_val;
        }else{
          return match;
        }
      }
    });
  return interpolatedString;
};

/**
 * Convert a string to a function.
 * @return function instance referred to by the string
 */
var stringToFunction = Iccbl.stringToFunction = function(str) {
  if (!str) return;
  
  var arr = str.split(".");

  var fn = (window || this);
  for (var i = 0, len = arr.length; i < len; i++) {
    fn = fn[arr[i]];
  }
  if ( typeof fn !== "function") {
    throw new ReferenceError("function not found: " + str);
  }
  return fn;
};


/**
 * Parse a string date value, ignoring the time and timezone.
 * 
 * @return Date
 * TODO: see Backgrid.DateFormatter.convert() method and refactor
 */
var dateParse = Iccbl.dateParse = function dateParse(rawData){

  if (_.isNull(rawData) || _.isUndefined(rawData)) return '';
  rawData = rawData.trim();
  if ((rawData + '').trim() === '') return null;

  if(ICCBL_DATE_RE.test(rawData)){
    var DDMMYYYY = ICCBL_DATE_RE.exec(rawData) || [];
    var jsDate = new Date(
      DDMMYYYY[3] * 1 || 0,
      DDMMYYYY[2] * 1 - 1 || 0,
      DDMMYYYY[1] * 1 || 0);
    console.log('date: raw: ', rawData, 'converted', jsDate);
    return jsDate;
  }else{
    // ISO date format, ignore timezone / time
    var temp = rawData.split(ISO_SPLITTER_RE)[0];
    if(DATE_RE.test(temp)){
      var YYYYMMDD = DATE_RE.exec(temp);
      if (Iccbl.appModel.DEBUG) console.log('YYYYMMDD', YYYYMMDD);
      var jsDate = new Date(YYYYMMDD[1]*1, YYYYMMDD[2]*1-1, YYYYMMDD[3]*1 )
      if (Iccbl.appModel.DEBUG) console.log('date: raw: ', rawData, 'converted', jsDate);
      return jsDate;
    }else{
      throw new Error('unrecognized date: ' + rawData );
    }
  }
};


/**
 * Generate an ISO date string from a JavaScript Date.
 * 
 * @param jsDate a JavaScript Date object
 * @return ISO date string for the Date object, ignoring the timezone.
 * - for internal representation of dates and communicating dates to the server.
 */
var getISODateString = Iccbl.getISODateString = function(jsDate){
  return jsDate && _.isDate(jsDate) ? jsDate.toISOString().split('T')[0] : jsDate;
  // equivalent:
  //  date = lpad(jsDate.getUTCFullYear(), 4, 0) 
  //    + '-' + lpad(jsDate.getUTCMonth() + 1, 2, 0) 
  //    + '-' + lpad(jsDate.getUTCDate(), 2, 0);
  //  return date;
};


/**
 * Generate a date string from a Javascript Date.
 * 
 * @param jsDate a JavaScript Date object
 * @return string representation of the Date object, ignoring the timezone.
 * - "ICCBL" format is "MM/DD/YYYY"
 */
var getIccblDateString = Iccbl.getDateString = function(jsDate){

  if (!jsDate) return jsDate;
  if (!_.isDate(jsDate)){
    // attempt to parse the date
    jsDate = Iccbl.dateParse(jsDate);
  }
  return ( 
      lpad(jsDate.getMonth() + 1, 2, 0) 
      + '/' + lpad(jsDate.getDate(), 2, 0) 
      + '/' + lpad(jsDate.getFullYear(), 4, 0) );
};

/**
 * Get ID keys from the model 
 * 
 * @param schema a resource definition as defined by the API
 * @param model Backbone.Model or Object described by the schema
 * @return array
 */
var getIdKeys = Iccbl.getIdKeys = function(model,schema) {
  
  if (_.has(schema, 'id_attribute')) {
  
    var id_attribute = schema['id_attribute'];
    console.log('create id from ' + id_attribute);
    var idList = [];
    
    _.each(id_attribute, function(item){
      var keyval;
      if (model instanceof Backbone.Model){
        keyval = model.get(item);
      }else{
        keyval = _.result(model,item);
      }
      if (_.isUndefined(keyval)){
        throw new TypeError('ID key: ' + item + ', not found on: ' + model);
      }
      idList.push(keyval);
    });
    return idList;
  } else {
    throw new TypeError("'id_attribute' not found on the schema: " 
            + JSON.stringify(schema)
            + ', for the model: ' + JSON.stringify(model.attributes));
  }
  
};

/**
 * Get the complete ID string for the model.
 * 
 * The "complete ID" is formed by joining the ID keys with the forward slash.
 */
var getIdFromIdAttribute = Iccbl.getIdFromIdAttribute = 
  function(model, schema){
    
    return Iccbl.getIdKeys(model,schema).join('/');
};

/**
 * Get the model key off the URI stack.
 *  
 * - pop one key, in order, for each of the key fields specified in 
 * the resource id_attribute.
 * 
 * @param resource - a resource definition as defined by the API
 * @param urlStack -
 *          array representation of the current unprocessed URI elements.
 * @param consumedStack -
 *          holds the items popped off the stack
 */
var popKeyFromStack = Iccbl.popKeyFromStack = function(
    resource, urlStack, consumedStack){
  
  var id  = '';
  var self = this;
  var checkStack = function(stack){
    if (_.isEmpty(urlStack)){
      var msg = 'not enough items on the URL to create the key for resource: ' + 
          resource.title + JSON.stringify(resource.id_attribute);
      throw msg;
    }
  };
  if(resource.key == 'apilog'){
    ref_resource_name = urlStack.shift();
    consumedStack.push(ref_resource_name);
    ref_resource = Iccbl.appModel.getResource(ref_resource_name);
    if(ref_resource && !_.isEmpty(ref_resource)){
      checkStack(urlStack);
      if (urlStack[0] == ref_resource_name){
        // if the key == resource, then this is a "parent log" (patch/post list)
        key = urlStack.shift();
        consumedStack.push(key);
      }else{
        key = self.popKeyFromStack(ref_resource,urlStack,consumedStack);
      }
    }
    checkStack(urlStack);
    date_time = urlStack.shift();
    consumedStack.push(date_time);
    id += [ref_resource_name,key,date_time].join('/');
  }else{
    _.each(resource.id_attribute, function(attribute){
      checkStack(urlStack);
      // don't care what the id is, just pop one for each
      var item = urlStack.shift();
      consumedStack.push(item);
      if ( id !== '' ){
        id += '/' + item;
      }
      else {
        id += item;
      }
    });
  }
  return id;
};

/**
 * Sort keys based on the associated ordinal for the key 
 * in the Resource.fields
 */
var sortOnOrdinal = Iccbl.sortOnOrdinal = function(keys, fieldHash) {
  var sorted = _(keys).sort(function(a, b) {
    if (!_.has(fieldHash, a) || !_.has(fieldHash, b)) {
      if (_.has(fieldHash, b)) {
        return -1;
      } else if (_.has(fieldHash, a)) {
        return 1;
      }
      return 0;
    }
    order_a = fieldHash[a]['ordinal'];
    // TODO: need an edit order by
    order_b = fieldHash[b]['ordinal'];
    if (_.isNumber(order_a) && _.isNumber(order_b)) {
      return order_a - order_b;
    } else if (_.isNumber(order_a)) {
      return -1;
    } else if (_.isNumber(order_b)) {
      return 1;
    } else {
      return 0;
    }
  });
  return sorted;
};

/**
 * Create a title from the schema "title_attribute".
 * 
 * Note: the title_attribute is an array of field specifiers and strings. 
 * If an array item is a field, the field value will be used, 
 * If an array item is not a field, then it will be concatenated directly.
 */
var getTitleFromTitleAttribute = Iccbl.getTitleFromTitleAttribute = 
  function(model, schema){
    var re_isQuoted = /['"]+/g;
    var fields = schema['fields'];
    if(_.has(schema, 'title_attribute')){
      var title = _.reduce(
        schema['title_attribute'],
        function(memo, item){
          if(item && item.match(re_isQuoted)){
            memo += item.replace(re_isQuoted, '');
          }else{
            if( model.has(item) ){
              var val = model.get(item);
              if (!_.has(fields,item)){
                throw Exception(
                  'Title property: ' + item + ', not present on model: ' 
                  + schema.key );
              }
              if (!_.isEmpty(fields[item].vocabulary_scope_ref)){
                val = Iccbl.appModel.getVocabularyTitle(
                  fields[item].vocabulary_scope_ref,val);
              }
              memo += val
            }else{
              memo += item
            }
          }
          return memo ;
        }, '');
      return title;
    }else{
      throw new TypeError("'title_attribute' not found on the schema: " + 
          JSON.stringify(schema)
          + ', for the model: ' + JSON.stringify(model.attributes));
    }
};


/**
 * Match URI fragment (matchstrings) to an array of specifiers.
 * 
 * Matches from the right to left; 
 * allowing URI fragments to match their parent URIs. 
 * Similar the the contains function, but using item.indexOf(matchString) ||
 * matchString.indexOf(item) for the truth test.
 */
var containsByMatch = Iccbl.containsByMatch = function(array, matchstring){
  return _.find(array, function(item) {
    var result = false;
    var index = matchstring.indexOf(item);
    if (index > -1 && index+item.length == matchstring.length ){
      result = true;
    }
    var index = item.indexOf(matchstring);
    if (!result && index > -1 && index+matchstring.length == item.length){
      result = true;
    }
    return result;
  });
};

/**
 * Utility Function:
 * createLabel: break long labels into multiple lines for display.
 * 
 * Split label strings on non-word characters and re-join into lines, 
 * where each line is less than max_line_length 
 * (or, if a part is longer than max_line_length, 
 * then that part becomes an entire line).
 */
var createLabel = Iccbl.createLabel = 
  function(original_label, max_line_length, break_char){
  
    var lines = [];
    var labelParts = original_label.split(/([\s_\-\.,]+)/);
    var line = '';
    _.each(labelParts, function(part){
      if(line.length > 0){
        var temp = line + part;
        if(temp.trim().length <= max_line_length){
          line = temp;
        }else{
          line = line.trim();
          if (line.length > 0) {
            lines.push(line);
          }
          line = part;
        }
      }else{
        if( part != ' '){
          line = part;
        }
      }
    });
    line = line.trim();
    lines.push(line);
    
    if(_.isUndefined(break_char)){
      break_char = '<br>';
    }
    return lines.join(break_char);
};


var formatResponseError = Iccbl.formatResponseError = function(response){
  
  var msg = '';
  var sep = '\n';
  if (!_.isUndefined(response.status))
      msg += response.status;
  if (!_.isUndefined(response.statusText))
      msg += sep + response.statusText;
  if (!_.isEmpty(response.responseText))
      msg += sep + response.responseText;
  return msg;
}

var UrlStack = Iccbl.UrlStack = Backbone.Model.extend({

  defaults: {
    // current app uri
    path: '',
    // current app uri, as array
    actualStack: [],
    // unprocessed uri elements
    currentStack: [],
    
    level: 0,
    // resources pop'd, one per level
    resources: [],
    // keys pop'd, one per level
    keys: []  
  },
  
  initialize: function(options) {
    this.path = options.path;
    this.actualStack = options.path.split('/');
    
  },
  
  pop: function() {
  },
  
});

var CollectionOnClient = Iccbl.CollectionOnClient = Backbone.Collection.extend({

  /**
   * Override collection parse method: Parse server response data.
   */
  parse : function(response) {
    console.log('Collection on client, parse called');
    if(_.has(response,'objects') && _.isArray(response.objects) ){
      return response.objects;
    }
    return response;
  },
  setSearch: function(){
    //nop, for now
  }
});


var getCollectionOnClient = Iccbl.getCollectionOnClient = 
  function(url, callback, options){
  
    var options = options || {};
    var data_for_get = options.data_for_get || {};
    data_for_get = _.extend({ limit: 0 },data_for_get);
    var CollectionClass = Iccbl.CollectionOnClient.extend({
      url: url 
    });
    var instance = new CollectionClass();
    instance.fetch({
      data: data_for_get,
      success: function(collection, response) {
        callback(collection);
      },
      always: function(){
        console.log('done: ');
      }
    }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
};

var getCollection = Iccbl.getCollection = 
  function(schemaResult, url, callback) {
    var CollectionClass = Iccbl.MyCollection.extend({
        url : url,
        defaults : {}
    });
    var instance = new CollectionClass();
    instance.fetch({
        success : function(collection) {
            callback(schemaResult, collection);
        },
        error : function(model, response, options) {
            var msg = 'Error locating resource: ' + url;
            var sep = '\n';
            if (!_.isUndefined(response.status))
                msg += sep + response.status;
            if (!_.isUndefined(response.statusText))
                msg += sep + response.statusText;
            if (!_.isEmpty(response.responseText))
                msg += sep + response.responseText;
            window.alert(msg);
            // TODO: use Bootstrap inscreen alert classed message div
        }
    });
};

var MyModel = Iccbl.MyModel = Backbone.Model.extend({

  url : function() {
    // Add trailing slash for Tastypie
    var url = Backbone.Model.prototype.url.call(this);
    return url + (url.charAt(url.length - 1) === '/' ? '' : '/');
  },

  initialize : function() {
    Backbone.Model.prototype.initialize.apply(this, arguments);
    var self = this;
    this.url = MyModel.prototype.url;
  },
});


///// Backgrid.Cell customizations /////

var BaseCell = Iccbl.BaseCell = Backgrid.Cell.extend({
  
  initialize: function (options) {
    
    Backgrid.Cell.prototype.initialize.apply(this, arguments);
    
    var self = this;
    this.model.on('change:'+this.column.get("name") , function(){
      // Block updates caused by adding columns
      if (!_.isUndefined(self.model.previous(self.column.get("name")))){
        self.$el.addClass('edited');
      }
    });
  }
});

var BooleanCell = Iccbl.BooleanCell = Backgrid.BooleanCell.extend({
  
  initialize: function(){
    
    BooleanCell.__super__.initialize.apply(this, arguments);
    var self = this;
    this.model.on('change:'+this.column.get("name") , function(){
      // Block updates caused by adding columns
      if (!_.isUndefined(self.model.previous(self.column.get("name")))){
        self.$el.addClass('edited');
      }
    });
  },
  
  // Set up to toggle the checkbox whenever the TD is clicked
  events: {
    'click': 'cellClicked'
  },
  
  cellClicked: function(){
    console.log('clicked...');
    var model = this.model, column = this.column;
    var checked = model.get(column.get("name"));
    model.set(column.get("name"), !checked);
  }
});
  
var StringFormatter = Iccbl.StringFormatter = function () {};

StringFormatter.prototype = new Backgrid.StringFormatter();

_.extend(StringFormatter.prototype, {
  /**
   * Extend Backgrid.StringFormatter to add spaces between values in arrays.
   */
  fromRaw: function (rawValue, model) {
    if (_.isUndefined(rawValue) || _.isNull(rawValue)) return '';
    if (_.isArray(rawValue)) {
      return rawValue.join(', ');
    } else {
      return rawValue + '';
    }
  }
});

var StringCell = Iccbl.StringCell = Backgrid.StringCell.extend({
  formatter: Iccbl.StringFormatter,
  initialize: function(){
    StringCell.__super__.initialize.apply(this, arguments);
    var self = this;
    this.model.on('change:'+this.column.get("name") , function(){
      // Block updates caused by adding columns
      if (!_.isUndefined(self.model.previous(self.column.get("name")))){
        self.$el.addClass('edited');
      }
    });
  }  
});

var LinkCell = Iccbl.LinkCell = Iccbl.BaseCell.extend({
  
  // TODO: redo the link cell like the UriListCell
  
  className : "link-cell",

  events : {
    'click A': 'linkCallback',
  },

  /**
   * @property {string} ["string with {model_key} values to interpolate"]
   */
  hrefTemplate: 'Http://', 

  /**
   * @property {string} [title] The title attribute of the generated anchor.
   *           It uses the display value formatted by the `formatter.fromRaw`
   *           by default.
   */
  title: null,

  /**
   * @property {string} [target="_blank"] The target attribute of the
   *           generated anchor.
   */
  target: "_self",
  
  initialize : function(options) {
    LinkCell.__super__.initialize.apply(this, arguments);
  },

  linkCallback: function(e){
    console.log('link clicked, override to handle', e);
  },
  
  render : function() {
    var self = this;
    this.$el.empty();
    var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
    var interpolatedVal = Iccbl.formatString(self.hrefTemplate,self.model);
    self.$el.append($('<a>', {
      tabIndex : -1,
      href : interpolatedVal,
      target : self.target,
      title: self.title
    }).text(formattedValue));
    return this;
  },
});


var UriListCell = Iccbl.UriListCell = Iccbl.BaseCell.extend({
  
  className : "uri-list-cell",

  /**
   * @property {string} ["string with {model_key} values to interpolate"]
   */
  hrefTemplate: 'Http://', 

  /**
   * @property {string} [title] The title attribute of the generated anchor. It
   *           uses the display value formatted by the `formatter.fromRaw` by
   *           default.
   */
  title: null,

  /**
   * @property {string} [target="_blank"] The target attribute of the generated
   *           anchor.
   */
  target: "_blank",
  
  initialize : function(options) {
    UriListCell.__super__.initialize.apply(this, arguments);
  },

  render : function() {
    var self = this;
    this.$el.empty();
    var rawValue = this.model.get(this.column.get("name"));
    
    if(rawValue && !_.isEmpty(rawValue)){
      var i = 0;
      _.each(rawValue, function(val){
        var interpolatedVal = Iccbl.formatString(self.hrefTemplate, self.model, val);
        if(i>0) self.$el.append(', ');
        self.$el.append($('<a>', {
          tabIndex : -1,
          href : interpolatedVal,
          title : val,
          target : self.target
        }).text(val));
        i++;
      });
    }
    return this;
  },

}); 

var ImageCell = Iccbl.ImageCell = Iccbl.BaseCell.extend({
  
  className : "image-cell",
  
  events : {
      "click #link" : "toLink",
  },

  initialize : function(options) {
    ImageCell.__super__.initialize.apply(this, arguments);
  },

  render : function() {
      this.$el.empty();
      this.$el.html(this.render_image());
      this.delegateEvents();
      return this;
  },
  
  render_image: function(){
    var val = this.model.get(this.column.get('name'));
    if (!_.isEmpty(val)){
      return '<img src="'+val+'" width="200" alt="" />';
    }else{
      return '';
    }
  },
      
});

var EditCell = Iccbl.EditCell = Iccbl.BaseCell.extend({
  
  className : "detail-cell",
  
  events : {
    "click #edit" : "editDetail",
  },

  initialize : function(options) {
    this.options = options;
    EditCell.__super__.initialize.apply(this, arguments);
  },

  render : function() {
    this.$el.empty();
    var formattedValue = this.formatter.fromRaw(
      this.model.get(this.column.get("name")));
    this.$el.append($("<a id='edit' >", {
      tabIndex : -1,
      href : '',
      title : formattedValue,
      // target : "_blank"
    }).text(formattedValue));

    this.delegateEvents();
    return this;
  },

  editDetail : function(e) {
    e.preventDefault();
    if(_.has(this.model, 'clickHandler')){
      this.model.clickHandler(this.model);
    }else{
      this.model.collection.trigger("MyCollection:detail", this.model);
    }
  },
});

var NumberFormatter = Backgrid.NumberFormatter;

var NumberCell = Iccbl.NumberCell = Backgrid.NumberCell.extend({
  
  initialize: function (options) {
    
    NumberCell.__super__.initialize.apply(this, arguments);
    
    var self = this;
    var model = this.model;
    var column = this.column;
    var currVal = model.get(column.get("name"));
    this.model.on('change:'+this.column.get("name") , function(){
      if (!_.isUndefined(self.model.previous(self.column.get("name")))){
        if (parseFloat(currVal) !== parseFloat(model.get(column.get("name")))) {
          self.$el.addClass('edited');
        }
      }
    });
  }
});

var DecimalFormatter = Iccbl.DecimalFormatter = function () {
  Backgrid.NumberFormatter.apply(this, arguments);
};
 
DecimalFormatter.prototype = new Backgrid.NumberFormatter();
 
_.extend(DecimalFormatter.prototype, {
 
  /**
   * Modified to pass non-numbers through.
   */
  fromRaw: function (number, model) {
    var args = [].slice.call(arguments, 1);
    if(_.isUndefined(number)){
      return '';
    }
    if(_.isNull(number)){
      return '';
    }
    if(!_.isNumber(number)){
      try{
        number = parseFloat(number);
      }catch(e){
        console.log('not a number: ' + number + ', ex:' + e);
        return number;
      }
    }
     
    args.unshift(number);
    return (NumberFormatter.prototype.fromRaw.apply(this, args) || "0");
  }
});  

/**
 * A number formatter that converts a floating point number, optionally
 * multiplied by a multiplier, to a units string and vice versa.
 * 
 * @class Backgrid.UnitsFormatter
 * @extends Backgrid.NumberFormatter
 * @constructor
 * @throws {RangeError}
 *           If decimals < 0 or > 20.
 */
var SIUnitsFormatter = Iccbl.SIUnitsFormatter = function () {
 Backgrid.NumberFormatter.apply(this, arguments);
};

SIUnitsFormatter.prototype = new Backgrid.NumberFormatter();

_.extend(SIUnitsFormatter.prototype, {

  /**
   * @member Backgrid.UnitsFormatter
   * @cfg {Object} options
   * 
   * @cfg {number} [options.multiplier=1] The number used to multiply the model
   *      value for display.
   * 
   * @cfg {string} [options.symbol='%'] The symbol to append to the Unitsage
   *      string.
   */
  defaults: _.extend({}, NumberFormatter.prototype.defaults, {
    siunits: [
      ['T', 1e12],
      ['G', 1e9],
      ['M', 1e6],
      ['k', 1e3],
      ['', 1],
      ['m', 1e-3,],
      ['μ', 1e-6,],
      ['u', 1e-6,],
      ['n', 1e-9 ],
      ['p', 1e-12 ]
      ],
    multiplier: 1
  }),

  /**
   * Extends Backgrid.NumberFormatter to support SI Units.
   * Takes a raw value from a model and returns an optionally formatted string.
   * 
   * NOTE: precision is lost; input string values are converted to floating 
   * point numbers.
   * 
   * Convert input to a floating point number, where the number is first multiplied by
   * `multiplier`, then converted to a formatted string like
   * NumberFormatter#fromRaw, then finally append `symbol` to the end.
   */
  fromRaw: function (number, model) {

    if (_.isNull(number) || _.isUndefined(number)) return '';

    return this.getUnit(number, this.multiplier, this.symbol, this.decimals);
  },

  /**
   * Convert the number to a siunit value
   * 
   * @return Math.round(number*multiplier,decimals) + symbol
   */
  getUnit: function(rawNumber, multiplier, symbol, decimals) {
     
    if(!_.isNumber(rawNumber)){
      try{
        number = parseFloat(rawNumber);
      }catch(e){
        console.log('not a number: ' + number + ', ex:' + e);
        return number;
      }
    }else{
      number = rawNumber;
    }
    if(number == 0 ) return number;
    var zeroPadding = 0;
    if (number<1){
      for (i=1; i<=rawNumber.length; i++) {
        var c = rawNumber.charAt(rawNumber.length-i);
        if (c=='.' || c!= '0'){
          break;
        }
        zeroPadding += 1;
      }
    }
    if(!_.isNumber(multiplier)){
      try{
        multiplier = parseFloat(multiplier);
      }catch(e){
        console.log('not a number - multiplier: ' + multiplier+ ', ex:' + e);
      }
    }

//    if(multiplier >= 1){
//      number = number * multiplier;
//    }else{
//      console.log("Error, DecimalCell multiplier < 1: " + multiplier);
//    }
     
    if(multiplier > 0){
      number = number * multiplier;
    }else{
      console.log("Error, DecimalCell multiplier ! > 0: " + multiplier);
    }
     
    pairUnit = _.find(this.siunits, function(pair){
      return pair[1] <= Math.abs(number); 
    });
     
    if(_.isUndefined(pairUnit)){
      console.log('could not find units for the input number: ' + number);
      return number;
    }
     
    var val = (1/pairUnit[1])*number;
    val = val.toFixed(decimals)*1;
    // val = Math.round(val*Math.pow(10,decimals))/Math.pow(10,decimals);

    var finalValue = '' + val;
    if (zeroPadding!=0){
      var newZeroPadding = 0;
      for (i=1; i<=finalValue.length; i++) {
        var c = finalValue.charAt(finalValue.length-i);
        if (c=='.' || c!= '0'){
          break;
        }
        newZeroPadding += 1;
      }
      var newDecimals = 0;
      if (finalValue.indexOf('.')>-1){
        newDecimals = finalValue.split('.')[1].length;
      }
      if (newZeroPadding < zeroPadding){
        if (newDecimals < decimals){
          newDecimals += zeroPadding;
        }
        if (newDecimals > decimals){
          newDecimals = decimals;
        }
        finalValue = val.toFixed(newDecimals);
      }
    }
    return finalValue + ' ' + pairUnit[0] + symbol;
  },
 
  
  /**
   * Extends Backgrid.NumberFormatter to support SI Units.
   * 
   * Takes a formatted string, usually from user input, and returns a
   * number for persistence in the model.
   * 
   * NOTE: precision is lost with zero padding; 
   * NumberFormatter converts string values to number values.
   */
  toRaw: function (formattedValue, model) {
    var self = this;
    var unit = '';
    var tokens;
    var rawValue, scaledRawValue;
    var unitMultiplier;
    
    if (formattedValue === '') return null;

    tokens = formattedValue.split(' ');

    rawValue = NumberFormatter.prototype.toRaw.call(this, tokens[0]);

    if (_.isUndefined(rawValue)) return rawValue;

    if (tokens.length == 2 ) {
      unit = tokens[1].trim().charAt(0);
    }
      
    var pairUnit = _.find(self.siunits, function(pair){
      return pair[0] == unit;
    });
    if (_.isUndefined(pairUnit)){
      console.log('unknown unit: ' + unit);
      // FIXME: indicate error to the user
      return null;
    }
    unitMultiplier = pairUnit[1];
    
    var originalPrecision = 0;
    if (tokens[0].indexOf('.')>-1){
      originalPrecision = (tokens[0] + "").split(".")[1].length;
    }
    
    // use the multiplier to scale the raw value
    newRawValue = rawValue / this.multiplier;
    // use the SIUnit to scale the raw value
    newRawValue= rawValue*unitMultiplier;
    
    // Truncate or pad the precision to the decimal setting
    var allowedPrecision = this.decimals;
    var desiredPrecision = originalPrecision;
    if (unitMultiplier < 1){
      allowedPrecision += (unitMultiplier + "").split(".")[1].length;
      desiredPrecision += (unitMultiplier + "").split(".")[1].length;
    }
    newRawValue = newRawValue.toFixed(allowedPrecision)*1;
    
    var newPrecision = 0;
    if ( (newRawValue+"").indexOf(".") > -1){
      newPrecision = (newRawValue+"").split(".")[1].length;
    }
    if (newPrecision < desiredPrecision){
      if (desiredPrecision < allowedPrecision){
        newRawValue = newRawValue.toFixed(desiredPrecision);
      }else{
        newRawValue = newRawValue.toFixed(allowedPrecision);
      }
    }
    
    // Convert to string to finalize precision with subsequent serialization
    return "" + newRawValue;
  },
 
});  

/**
 * A DecimalCell is another Backgrid.NumberCell that takes a floating number,
 * and showing a decimals number of digits.
 * 
 * @class Backgrid.DecimalCell
 * @extends Backgrid.NumberCell
 */
var DecimalCell = Iccbl.DecimalCell = NumberCell.extend({
  
  /** @property */
  className: "decimal-cell",

  /** @property {Backgrid.CellFormatter} [formatter=Backgrid.NumberFormatter] */
  formatter: DecimalFormatter,
 
  /**
   * Initializes this cell and the Units formatter.
   * 
   * @param {Object}
   *          options
   * @param {Backbone.Model}
   *          options.model
   * @param {Backgrid.Column}
   *          options.column
   */
  initialize: function () {
    DecimalCell.__super__.initialize.apply(this, arguments);
    var formatter = this.formatter;
    formatter.decimals = this.decimals;
  }
 
 });  

/**
 * A SIUnitsCell is another Backgrid.NumberCell that takes a floating number,
 * optionally multiplied by a multiplier, showing a decimals number of digits,
 * and displayed with a units symbol.
 * 
 * @class Backgrid.SIUnitsCell
 * @extends Backgrid.NumberCell
 */
var SIUnitsCell = Iccbl.SIUnitsCell = NumberCell.extend({

  /** @property */
  className: "units-cell",

  /** @property {number} [multiplier=1] */
  multiplier: SIUnitsFormatter.prototype.defaults.multiplier,

  /** @property {string} [symbol='%'] */
  symbol: SIUnitsFormatter.prototype.defaults.symbol,

  /** @property {Backgrid.CellFormatter} [formatter=Backgrid.UnitsFormatter] */
  formatter: SIUnitsFormatter,

  /**
   * Initializes this cell and the Units formatter.
   * 
   * @param {Object}
   *          options
   * @param {Backbone.Model}
   *          options.model
   * @param {Backgrid.Column}
   *          options.column
   */
  initialize: function () {
    SIUnitsCell.__super__.initialize.apply(this, arguments);
    var formatter = this.formatter;
    formatter.multiplier = this.multiplier;
    formatter.decimals = this.decimals;
    formatter.symbol = this.symbol;
  }

});  

var SelectCell = Iccbl.SelectCell = Backgrid.SelectCell.extend({
  /** 
   * override Backgrid.SelectCell:
   * - render to return the cell value if optionValues is malformed or missing 
   * the value
   */
  render: function () {
    this.$el.empty();

    var optionValues = _.result(this, "optionValues");
    var model = this.model;
    var rawData = this.formatter.fromRaw(
      model.get(this.column.get("name")), model);
    var selectedText = [];

    if (_.isArray(optionValues) &&  !_.isEmpty(optionValues)){
      for (var k = 0; k < rawData.length; k++) {
        var rawDatum = rawData[k];
        for (var i = 0; i < optionValues.length; i++) {
          var optionValue = optionValues[i];
          if (_.isArray(optionValue)) {
            var optionText  = optionValue[0];
            var optionValue = optionValue[1];
            if (optionValue.toLowerCase() == rawDatum.toLowerCase()){
              selectedText.push(_.escape(optionText));
            }
          }
          else if (_.isObject(optionValue)) {
            var optionGroupValues = optionValue.values;
            for (var j = 0; j < optionGroupValues.length; j++) {
              var optionGroupValue = optionGroupValues[j];
              if (optionGroupValue[1].toLowerCase() == rawDatum.toLowerCase()) {
                selectedText.push(_.escape(optionGroupValue[0]));
              }
            }
          }
          else {
            throw new TypeError;
          }
        }
      }
      
    }
    var isEmpty = ( _.isEmpty(rawData) 
        || (_.isArray(rawData) && rawData.length == 1 && _.isEmpty(rawData[0])))
    if( !isEmpty && _.isEmpty(selectedText)){
      selectedText = rawData;

      Iccbl.appModel.error(Iccbl.formatString(
        'column: {column}, vocabulary: {vocabulary} is misconfigured: rawData: {rawData}',
        { column: this.column.get("name"),
          vocabulary: _.result(this, "vocabulary_scope_ref"),
          rawData: rawData 
        }));
      console.log(Iccbl.formatString(
        'column: {column}, vocabulary: {vocabulary} is misconfigured,' 
          + 'rawData: {rawData}, optionValues: {optionValues}',
        { column: this.column.get("name"),
          vocabulary: _.result(this, "vocabulary_scope_ref"),
          rawData: rawData,
          optionValues: optionValues
        }));
    }
    var finalText = selectedText.join(this.delimiter);
    
    if(!_.isUndefined(this.hrefTemplate)){
      // hack - if hrefTemplate is defined, treat this like a link cell - 20150828
      var target = this.target || "_self";
      var interpolatedVal = Iccbl.formatString(this.hrefTemplate,this.model);
      this.$el.append($('<a>', {
        tabIndex : -1,
        href : interpolatedVal,
        target : this.target
      }).text(finalText));
      
    }else{
      this.$el.append(finalText);
    }

    this.delegateEvents();

    return this;
  }

});

/**
 * Override Backgrid DateTimeFormatter
 * - recognize user input in the format MM/DD/YYYY
 * - but also recognize ISO 8601 format, so values from the server are parsed
 */
var DatetimeFormatter = Iccbl.DatetimeFormatter = function (options) {
  _.extend(this, this.defaults, options || {});

  if (!this.includeDate && !this.includeTime) {
    throw new Error("Either includeDate or includeTime must be true");
  }
};

DatetimeFormatter.prototype = new Backgrid.DatetimeFormatter();
_.extend(DatetimeFormatter.prototype, {

  fromRaw: function (rawData, model) {
    if (_.isNull(rawData) || _.isUndefined(rawData)) return '';
    rawData = rawData.trim();
    if ((rawData + '').trim() === '') return null;
    
    return getIccblDateString(Iccbl.dateParse(rawData));
  },

  toRaw: function(formattedData, model){
    if (_.isNull(formattedData) || _.isUndefined(formattedData)) return '';
    if ((formattedData + '').trim() === '') return null;
    if(ICCBL_DATE_RE.test(formattedData)){
      var DDMMYYYY = ICCBL_DATE_RE.exec(formattedData) || [];
      var jsDate = new Date(
        DDMMYYYY[3] * 1 || 0,
        DDMMYYYY[2] * 1 - 1 || 0,
        DDMMYYYY[1] * 1 || 0);
      var temp = getISODateString(jsDate);
      return temp;
    }else{
      return;
    }
  },
  
  /**
   * Modify Backgrid DatetimeFormatter convert:
   * - ignore timezone - remove UTC conversion
   * - add in ICCBL_DATE_RE
   * TODO: refactor, using Iccbl.dateParse
   */
  _convert: function (data, validate) {
    if ((data + '').trim() === '') return null;

    var date, time = null;
    if (_.isNumber(data)) {
      var jsDate = new Date(data);
      date = lpad(jsDate.getFullYear(), 4, 0) 
        + '-' + lpad(jsDate.getMonth() + 1, 2, 0) 
        + '-' + lpad(jsDate.getDate(), 2, 0);
      time = lpad(jsDate.getHours(), 2, 0) 
        + ':' + lpad(jsDate.getMinutes(), 2, 0) 
        + ':' + lpad(jsDate.getSeconds(), 2, 0);
      // modified 20150831 - use local date/time
      //date = lpad(jsDate.getFullYear(), 4, 0) + '-' + lpad(jsDate.getMonth() + 1, 2, 0) + '-' + lpad(jsDate.getDate(), 2, 0);
      //time = lpad(jsDate.getHours(), 2, 0) + ':' + lpad(jsDate.getMinutes(), 2, 0) + ':' + lpad(jsDate.getUTCSeconds(), 2, 0);
    }
    else {
      data = data.trim();
      var parts = data.split(ISO_SPLITTER_RE) || [];
      date = ICCBL_DATE_RE.test(parts[0]) ? parts[0] : '';
      time = date && parts[1] ? parts[1] : TIME_RE.test(parts[0]) ? parts[0] : '';
    }
    // FIXME: review this 
    var DDMMYYYY = ICCBL_DATE_RE.exec(date) || [];
    var HHmmssSSS = TIME_RE.exec(time) || [];

    if (validate) {
      if (this.includeDate && _.isUndefined(DDMMYYYY[0])) return;
      if (this.includeTime && _.isUndefined(HHmmssSSS[0])) return;
      if (!this.includeDate && date) return;
      if (!this.includeTime && time) return;
    }

    var jsDate = new Date(DDMMYYYY[3] * 1 || 0,
                                   DDMMYYYY[2] * 1 - 1 || 0,
                                   DDMMYYYY[1] * 1 || 0,
                                   HHmmssSSS[1] * 1 || null,
                                   HHmmssSSS[2] * 1 || null,
                                   HHmmssSSS[3] * 1 || null,
                                   HHmmssSSS[5] * 1 || null);

    var result = '';

    if (this.includeDate) {
      result = ( getIccblDateString(jsDate) );
    }

    if (this.includeTime) {
      result = ( result 
          + (this.includeDate ? 'T' : '') 
          + lpad(jsDate.getHours(), 2, 0) 
          + ':' + lpad(jsDate.getMinutes(), 2, 0) 
          + ':' + lpad(jsDate.getSeconds(), 2, 0)
          );

      if (this.includeMilli) {
        result = result + '.' + lpad(jsDate.getMilliseconds(), 3, 0);
      }
    }

    if (this.includeDate && this.includeTime) {
      result += "Z";
    }
    return result;
  },
  
});

/**
 * Override DateCell
 * - set the format to MM/DD/YYYY
 */
var DateCell = Iccbl.DateCell = Backgrid.DateCell.extend({

  /** @property {Backgrid.CellFormatter} [formatter=Backgrid.DatetimeFormatter] */
  formatter: DatetimeFormatter,

  /**
  Initializes this cell and the datetime formatter.

  @param {Object} options
  @param {Backbone.Model} options.model
  @param {Backgrid.Column} options.column
  */
  initialize: function (options) {
    // Note __super__ == Backgrid.DateCell.prototype
    DateCell.__super__.initialize.apply(this, arguments);
    var formatter = this.formatter;
    //   formatter.includeDate = this.includeDate;
    //   formatter.includeTime = this.includeTime;
    //   formatter.includeMilli = this.includeMilli;

    var placeholder = "MM/DD/YYYY";
  
    this.editor = this.editor.extend({
      attributes: _.extend(
          {}, this.editor.prototype.attributes, this.editor.attributes, {
            placeholder: placeholder
          })
    });
  }
  
});

/**
 * uses the options.attributes.label
 */
var DeleteCell = Iccbl.DeleteCell = Iccbl.BaseCell.extend({
  className: "delete-cell",
  events : {
  "click #delete" : "delete"
  },

  initialize: function(options){
    DeleteCell.__super__.initialize.apply(this, arguments);
  },

  render: function () {
    this.$el.empty();
  
    this.$el.append("&nbsp;");
    this.$el.append(
      $("<a id='delete' >", {
        tabIndex : -1,
        href : '',
      }).text(this.column.attributes['text']));
    this.delegateEvents();
    return this;
  },

  delete: function(e) {
    e.preventDefault();
    this.model.collection.trigger("MyCollection:delete", this.model);
  }
});

var CollectionInColumns = Iccbl.CollectionInColumns = Backbone.Collection.extend({
  /**
   * Override collection parse method: Parse server response data.
   */
  parse : function(response) {
    console.log('Collection on client, parse called');
    var pivoted = {};
    var i = 0;
    _.each(response.objects, function(obj) {
      _.pairs(obj, function(pair) {
        if (_.has(pivoted, pair[0])) {
          pivoted[pair[0]] = {};
        }
        pivoted[pair[0]][i] = pair[1];
      });
      i++;
    });
    return _.values(pivoted);
  },
});

var UriContainerView = Iccbl.UriContainerView = Backbone.Layout.extend({
  initialize: function(args) {
    console.log('initialize UriContainerView');
    var model = this.model = args.model;
    var targetProperty = args.property || 'uriStack';
    this.listenTo(model, 'change:'+targetProperty , this.uriStackChange );
    
    Backbone.View.prototype.initialize.apply(this,arguments);
  },
  
  /**
   * This method will report URI stack change events from child views.
   */
  reportUriStack: function(reportedUriStack, options) {
    var options = options || {source: this};
    var consumedStack = this.consumedStack || [];
    var actualStack = consumedStack.concat(reportedUriStack);
    this.model.set({'uriStack': actualStack}, options);     
  },
  
  /**
   * Backbone.Model change event handler
   * 
   * @param options.source =
   *          the event source triggering view
   */
  uriStackChange: function(model, val, options) {
    if(options.source === this){
      console.log('self generated uristack change');
      return;
    }else{
      var uriStack = _.clone(this.model.get('uriStack'));
      try {
        this.changeUri(uriStack);
      }catch (e){
        console.log('error thrown: ',e);
        Iccbl.appModel.error('error: ' + e);
      }
    }
  },
  
  changeUri: function(uriStack) {
    window.alert(
      'ContentView changeUri function must be implemented. uriStack: ' +
      JSON.stringify(uriStack));
  }
});

var MultiSortBody = Iccbl.MultiSortBody = Backgrid.Body.extend({
  
  /**
   * See Backgrid.Body.sort: - created to solve the multisort case for the
   * server side backbone-pageable collection only. triggered by "backgrid:sort" -
   * sent from the column header cells.
   */
  sort: function (column, direction) {
    if (_.isString(column)) column = this.columns.findWhere({name: column});

    console.log('MultiSortBody.sort( ' + column.get('name') + ', ' + direction);

    var collection = this.collection;
    var order;
    if (direction === "ascending") order = -1;
    else if (direction === "descending") order = 1;
    else order = null;
    
    collection.setSorting(column.get("name"), order,
        {sortValue: column.sortValue()});
    collection.fetch({
      reset: true, 
      success: function () {
        console.log('fetch success, direction: ' + direction);
        collection.trigger("backgrid:sorted", column, direction, collection);
      }
    }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
    
    column.set("direction", direction);

    return this;
  }
});

var MyCollection = Iccbl.MyCollection = Backbone.PageableCollection.extend({

  initialize : function(options) {
    var self = this;
    this.options = options;
    Backbone.PageableCollection.prototype.initialize.apply(this, options);
  },

  mode: 'server',

  url : function() {
      return this.url;
  },
  searchHash : {},
  model : MyModel,
  state : {
    pageSize : 25,  // TODO: probably not necessary
  },
  
  // PageableCollection.fetch() uses the queryParams attribute to interpret
  // the server response and to determine the data hash sent to the server.
  // Adjust the query params for tastypie.
  queryParams : {
    pageSize : 'limit',
    offset : function() {
      return (this.state.currentPage - 1) * this.state.pageSize;
    },
    use_vocabularies: null, // this signals to the api to replace out
                            // vocabularies - FIXME: make this a setting?
    totalRecords : null, // unset for tastypie
    totalPages : null, // unset for tastypie
    sortKey : "order_by", // modified for tastypie
    order : null, // unset for tastypie
    order_by : function() {// modified for tastypie: use only
      if ( typeof this.state !== 'undefined' 
            && this.state.orderStack
            && this.state.orderStack.length ) {
        // Note: convert the orderStack using "traditional" array serialization
        // see: http://api.jquery.com/jQuery.param/
        // var val = this.state.orderStack.join('&order_by=');
        return this.state.orderStack;
      }
    }, 
    includes: function(){
      return this.listModel.get('includes');
    }
  },

  /**
   * Override
   */
  parseState : function(response, queryParams, state, options) {
    // hack the response for tastypie:
    // note, this is because the pageable collection doesn't work with the
    // backbone-tastypie.js fix
    var state = _.clone(state);
    if(!response.meta && response.meta.total_count){
      msg = 'error in response: no meta information';
      Iccbl.appModel.error(msg);
      throw Exception(msg);
    }
    state.totalRecords = response.meta.total_count;
    
    // FIXME: having to set this for the pre-fetched collections
    if(!_.isNumber(state.firstPage)) state.firstPage = 1;
    
    if (Math.ceil(state.totalRecords / state.pageSize) < state.currentPage) {
      console.log('adjust currentPage');
      state.currentPage = 1;
    }
    return state;
  },
  
  /**
   * Override
   */
  parseRecords: function (resp, options) {
    return resp.objects;
  },
      
  /**
   * Method for external callers to set the search, with fetch
   * @param options - for Backbone.Collection.fetch:
   * { reset: false } (default) - uses set to (intelligently) merge the fetched 
   * models ("add" events are fired),
   * {reset: true}, in which case the collection will be (efficiently) reset 
   * (no "add" events will be fired)
   */
  setSearch: function(searchHash, options) {
    var self = this;
    var searchHash = _.clone(searchHash);
    self.listModel.set('search', searchHash);
    
    // Tell all the header cells
    this.trigger("MyServerSideFilter:search", searchHash, this);

    // TODO: debug: "_data" is not needed
    // backbone.paginator should translate all queryparams into "data" in the 
    // fetch method
    //
    // Allow searches that aren't for a visible column:
    // - if the search key is not in the queryParams, then it is not a column
    // - this will add it manually to the queryParams (which are serialized in
    // the fetch to the server).
    var _data = {};
    _.each(_.keys(searchHash), function(key) {
      var val = searchHash[key]

      if(_.isEmpty("" + val)){
          delete self.queryParams[key];
      }else{
        // Check if param dne, or if param exists and is a value to be set.
        // The reason for the "isFunction" check is that the Backgrid-filter
        // defined params are function calls to get the current value in the
        // searchbox - so skip those as state is stored there.
        if (!_.has(self.queryParams, key) || !_.isFunction(self.queryParams[key])) {
        	_data[key]=val;
        	// Make the params persistent (if not a backgrid-filter)
          self.queryParams[key] = function () {
            var search = self.listModel.get('search');
            if (!_.isEmpty(search)){
              return _.result(search, key, null);
            }
            return null;
          };
        }
      }
    });

    //self.fetch({reset: true}).fail(
    self.fetch(options).fail(
      function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments);}
    );
  },

  /**
   * Proxy for the search elements to add search terms to the listModel,
   * if options.reset == true, then fetch, otherwise no fetch.
   */ 
  addSearch: function(searchHash, options) {
    console.log('addSearch: ' + JSON.stringify(searchHash) 
      + ', options: ' + JSON.stringify(options) );
    var self = this;
    var newSearchHash = _.clone(self.listModel.get('search'));
    newSearchHash = _.extend(newSearchHash, searchHash);
    if(options && options.reset && self.state.currentPage != 1){
      self.state.currentPage = 1;
    }
    self.listModel.set({
      'search' : newSearchHash
    });
// Removed 20170321
//    _.each(_.keys(searchHash), function(key){
//      // make the params persistent (if not a backgrid-filter)
//      self.queryParams[key] = function () {
//        return self.listModel.get('search')[key] || null;
//      };
//    });
// Removed 20161213, see search listener in list2.js
//    if(options && options.reset && self.state.currentPage != 1){
//      console.log('collection.addSearch: reset');
//      self.getFirstPage({reset: true, fetch: true}).fail(
//        function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments);}
//      );
//    }
  },

  /**
   * Proxy for the search elements to clear search terms from the listModel on
   * the collection.
   * if options.reset == true, then fetch, otherwise no fetch.
   */ 
  clearSearch: function(searchKeys, options) {
//    console.log('clearsearch: ' + JSON.stringify(searchKeys) 
//        + ', options: ' + JSON.stringify(options) );
    var self = this;
    var searchHash = {};
    var found = false;
    if (!_.isUndefined(searchKeys)) {
      searchHash = _.clone(self.listModel.get('search'));
      _.each(searchKeys, function(searchKey) {
        if(_.has(searchHash, searchKey)){
          delete searchHash[searchKey];
          found = true;
        }
      });
    }
    if(found){
      
      self.state.currentPage = 1;
      self.listModel.set({
        'search' : searchHash
      }, options);
// Removed 20161213, see search listener in list2.js
//      if(options && options.reset){
//        console.log('collection.clearSearch: reset');
//        self.getFirstPage({reset: true, fetch: true}).fail(
//          function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments);}
//        );
//      }
    }
  },

  /**
   * Override - HeaderCell.onClick -> backgrid:sort -> body.sort -> ->
   * (BackbonePageable)Collection.setSorting -> Collection.fetch() -> grab data
   * from state based on "queryParams"
   */
  setSorting : function(sortKey, order, options) {
    var state = this.state;
    
    var orderStack = state.orderStack || [];
    
    var newdir = order == 1 ? '-' : order == -1 ? '': null;
    
    var newStack = [];
    var found = false;
    
    _.each(orderStack, function(order_entry){
      var dir = order_entry.substring(0,1);
      var fieldname = order_entry;
      if(dir == '-'){
        fieldname = order_entry.substring(1);
      }else{
        dir = '';
      }
      if(fieldname == sortKey){
        found = true;
        if(newdir === null){
          // pop this off
        }else if(newdir == dir){
          // no change; push back on the stack
          newStack.push(order_entry);
        }else if (newdir !== dir){
          newStack.push(newdir + fieldname);
        }
      }else{
        newStack.push(order_entry);
      }
    });
    
    if(!found && newdir !== null) newStack.push(newdir + sortKey);
    
    console.log('Ordering update: old: ' + JSON.stringify(orderStack) 
        + ', new: ' + JSON.stringify(newStack));
    state.orderStack = newStack;
    
    // Backbone.PageableCollection.prototype.setSorting.call(this, sortKey,
    // order);
	  // TODO: Investigate why PageableCollection.setSorting is not triggering
	  // a 'sort' event (needed to clear old sort indicators).
	  // Sequence of a sort:
	  // Backgrid.HeaderCell.onClick-> collection.trigger('backgrid:sort')
	  // -> 'backgrid:sort' -> Backgrid.Body.sort()
	  // -> PageableCollection.setSorting(): sets state.sortKey
	  // -> if(fullCollection) (client mode) collection.sort()
	  // -> else
	  // -> Collection.fetch(reset:true)
	  // *so in this case, no sort(), if reset:false, then a "set"
	  // would be called, and a 'sort' triggered
	  // * also calls column.set to put the new indicator
	  // column.set('direction')
	  // without a sort, there is no erasing of the old sort indicators:
	  // Collection.sort() -> trigger('sort') ->
	  // -> Backgrid.HeaderCell.removeCellDirection
	  // Last note: this may be caused by not getting the sortKey from the
	  // queryParams on parseState.
  	this.trigger('sort',this); 
  },

});

//// Header Cell Definitions /////

var MultiSortHeaderCell = Iccbl.MultiSortHeaderCell = Backgrid.HeaderCell.extend({

  ___klass: 'MultiSortHeaderCell',
  
  initialize : function(options) {
    this.options = options;
    MultiSortHeaderCell.__super__.initialize.apply(this, arguments);

    this.fieldinformation = _.clone(this.column.get('fieldinformation'));
    
    this.listenTo(this.collection,"sort",this.collectionSorted);
    this.listenTo(this.collection,"Iccbl:clearSorts", this.removeCellDirection);
    _.bindAll(this, '_submit', 'clearSearch');
  },
  
  // Original onClick, for reference, from Backgrid
  // onClick: function (e) {
  // e.preventDefault();
  //
  // var column = this.column;
  // var collection = this.collection;
  // var event = "backgrid:sort";
  //
  // function cycleSort(header, col) {
  // if (column.get("direction") === "ascending") collection.trigger(event, col,
  // "descending");
  // else if (column.get("direction") === "descending")
  // collection.trigger(event, col, null);
  // else collection.trigger(event, col, "ascending");
  // }
  //
  // function toggleSort(header, col) {
  // if (column.get("direction") === "ascending") collection.trigger(event, col,
  // "descending");
  // else collection.trigger(event, col, "ascending");
  // }
  //
  // var sortable = Backgrid.callByNeed(column.sortable(), column,
  // this.collection);
  // if (sortable) {
  // var sortType = column.get("sortType");
  // if (sortType === "toggle") toggleSort(this, column);
  // else cycleSort(this, column);
  // }
  // },
      
  // TODO: debounced clicking - sort of working, but...
  // - at this time this method is debouncing on the cell instance:
  // --- no coordination with other headercells
  // 
  /**
   * Event handler for the `click` event on the cell's anchor. If the column is
   * sortable, clicking on the anchor will cycle through 3 sorting orderings -
   * `ascending`, `descending`, and default.
   */
  onClick: function (e) {
    var self=this;
    e.preventDefault();
    e.stopPropagation();
    
    var collection = this.collection;
    var event = "backgrid:sort"
    var column = this.column;
    
    if(_.isUndefined(this.tempdirection)){
      this.tempdirection = this.column.get("direction");
    }
    console.log('tempdirection1: ' + this.tempdirection);
    if(this.tempdirection == "ascending"){
      this.tempdirection = "descending";
    }else if(this.tempdirection == "descending"){
      this.tempdirection = "none";
    }else{
      this.tempdirection = "ascending";
    }
    
    this.setCellDirection(
      column, self.tempdirection=="none"?null:self.tempdirection );

    var args = arguments;
    
    var delayedClick = function(){
      console.log('delayedclick: tempdirection: ' + self.tempdirection);
      if(self.tempdirection !== self.lastExecutedVal){
        console.log('delayed click: ' + self.tempdirection );
        collection.trigger(
          event, column, self.tempdirection=="none"?null:self.tempdirection);
        self.lastExecutedVal = self.tempdirection;
      }else{
        console.log(
          'this.tempdirection == self.lastExecutedVal: ' + self.lastExecutedVal);
      }
    };
    _.debounce(delayedClick, 750)();
    // FIXME: both throttle and debounce seem to work the same;
    // that is they both are working like setTimeout
    // _.throttle(delayedClick, 5000, {leading: false})();
    
    console.log('onclick exit');
  },   
 
  collectionSorted: function(collection, options){
    var self = this;
    var name = this.column.get('name');
    var state = this.collection.state;

    var i = 0;
    if (state) {
      
      _.each(state.orderStack, function(order_entry){
        i++;
        var dir = order_entry.substring(0,1);
        var direction = null;
        var fieldname = order_entry;
        if(dir == '-'){
          fieldname = order_entry.substring(1);
          direction = 'descending';
        }else{
          dir = '';
          direction = 'ascending';
        }
        if(fieldname == name){
          self.$el.removeClass("ascending").removeClass("descending");
          self.$el.addClass(direction); 
  
          var sorter = self.$el.find('#sorter');
          sorter.empty();
          sorter.append(
           "<span style='margin-bottom: 2px;' class='badge pull-right'>" 
           + i + "<b class='sort-caret'></b></span>");
        }
      });
    }
  },
 
  /**
   * Event handler for the column's `change:direction` event. If this
   * HeaderCell's column is being sorted on, it applies the direction given as a
   * CSS class to the header cell. Removes all the CSS direction classes
   * otherwise.
   */
  setCellDirection: function (column, direction) {
    var self = this;
    var name = column.get('name');
    
    if(_.isUndefined(direction) || _.isNull(direction)){
      // this.$el.removeClass("ascending").removeClass("descending");
      // this.$el.find("#sorter").empty();
      this.removeCellDirection();
    }else{
      this.$el.removeClass("ascending").removeClass("descending");
      this.$el.addClass(direction);
       
      var num = 1;
      if(self.collection.state.orderStack &&
          !_.isEmpty(self.collection.state.orderStack)){
        var i = 0;
        var found = _.find(self.collection.state.orderStack, function(fieldname){
          i++;
          if(fieldname == name || fieldname == '-' + name){ 
            num = i;
            return true;
          }
        });
        if(!found){
          // TODO: is this still necessary?
          num = self.collection.state.orderStack.length+1; 
        }
        
      }
      sorterText = $("<span style='margin-bottom: 2px;' class='badge pull-right'>" 
          + num + "<b class='sort-caret'></b></span>");
       
      self.sorter.empty();
      self.sorter.append(sorterText);
    }
  },

  /**
   * Backgrid event handler for the PageableCollection 'sort' event.
   */
  removeCellDirection: function () {
    var self = this;
    this.$el.removeClass("ascending").removeClass("descending");
    if(self.sorter) self.sorter.empty();
  },

   /**
   * Renders a header cell with a sorter and a label.
   */
  render : function() {
    var self = this;
    this.$el.empty();
    var column = this.column;
    var sortable = Backgrid.callByNeed(column.sortable(), column, this.collection);
    if(sortable){
      var label = Iccbl.createLabel(column.get("label"), 10);
      self.sorter = $("<div id='sorter'></div>");
      label = $("<a>" + label +"</a>").append(self.sorter);
    } else {
      // NOTE: using anchor node to set the text color/style the same as other
      // cells
      // label = document.createTextNode(column.get("label"));
      label = $("<a>" + column.get("label") +"</a>");
    }
    this.$el.append(label);
    this.$el.addClass(column.get("direction"));
    this.$el.addClass(column.get("name"));
    this.delegateEvents();
    
    var mouseover = this.options['column']['attributes']["description"];  
    if (this.options['column'].has('mouseover')){
      mouseover = this.options['column'].get('mouseover');
    }
    this.$el.prop('title', mouseover);
    
    return this;
  }
    
}); // end MultiSortHeaderCell

var FilterHeaderCell = Iccbl.FilterHeaderCell = Iccbl.MultiSortHeaderCell.extend({
  
  filtericon_text : '<span class="pull-left glyphicon glyphicon-search" ' +
    ' id="filter-icon" ></span>',
  expandicon_text : '<span ' +
    ' class="pull-left glyphicon glyphicon-chevron-down" ' +
    ' id="expand-filter-icon" > ' +
    ' <span id="expand-filter-icon-text" ></span></span>',
  collapseicon_text : '<span ' +
    ' class="pull-left glyphicon glyphicon-chevron-up"'+
    ' id="collapse-filter-icon" ></span>',
  
  events: {
    'filter_clear': 'clearSearch'
  },
  
  
  initialize : function(options) {
    var self = this;
    FilterHeaderCell.__super__.initialize.apply(this, arguments);

    this.filterIcon = $(this.filtericon_text);
    this.collapseIcon = $(this.collapseicon_text);
    this.expandIcon = $(this.expandicon_text);
    this.expandIconText = this.expandIcon.find('#expand-filter-icon-text');
    
    this.fieldinformation = options.fieldinformation || this.fieldinformation;
    if (_.isUndefined(this.fieldinformation)){
      throw Exception('must define a fieldinformation');
    }
    this.serverSideFilter = options.serverSideFilter || this.serverSideFilter;
    if (_.isUndefined(this.serverSideFilter)){
      throw Exception('must define a serverSideFilter');
    }
    
    this._serverSideFilter = new this.serverSideFilter(_.extend({
      columnName: this.column.get('name'),
      fieldinformation: this.fieldinformation
    }, options));

    this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
    this.listenTo(this.collection,"Iccbl:clearSearches",this.clearSearch);

    _.bindAll(this, '_submit', 'clearSearch');
  },
  
  /**
   * Listen for router generated search events
   */
  _search: function(hash, collection){
    var self = this;
    var name = this.column.get('name');
    var searchHash = _.clone(hash);

    // TODO: could use form.isSet() instead of found
    var found = this._serverSideFilter._search(searchHash);
    
    if(found){
      self.$el.addClass('filtered');
      self.expandIconText.html(
        self._serverSideFilter._printSearchHash(
          self._serverSideFilter._getSearchHash()));
    }else{
      self.$el.removeClass('filtered');
      self.expandIconText.empty();
    }
  },  
  
  clearSearch: function(options){
    var self=this;
    self._serverSideFilter.clear();
    self._serverSideFilter.$el.hide();
    self.filterIcon.show();
    self.collapseIcon.hide();
    self.expandIconText.empty();
    self.expandIcon.hide();
  },      
  
  _submit: function(e){
    var self  = this;
    console.log('_submit called');
    if (e) e.preventDefault();      
  
    var searchHash = self._serverSideFilter._submit();
    if(!_.isEmpty(searchHash)){
      var possibleSearches = self._serverSideFilter.getPossibleSearches();
      self.collection.clearSearch(possibleSearches, {silent: true});
      
      console.log('server side filter add search: ' + 
          JSON.stringify(searchHash));
      this.collection.addSearch(searchHash,{reset: true});
    }else{
      console.log('nothing submitted');
    }
  },
  
  render : function() {
    var self = this;
    FilterHeaderCell.__super__.render.apply(this);
    if (_.result(this.fieldinformation,'filtering') !== true){
      return this;
    }
  
    this._serverSideFilter.render();
    this.$el.append(this._serverSideFilter.el);
    this._serverSideFilter.$el.hide();
    this.$el.append(this.filterIcon);
    this.$el.append(this.collapseIcon);
    this.$el.append(this.expandIcon);

    this._serverSideFilter.clearButton().click(function(e){
      e.preventDefault();
      e.stopPropagation();
      self.clearSearch();
    
      var possibleSearches = self._serverSideFilter.getPossibleSearches();
      self.collection.clearSearch(possibleSearches);
    });
    
    this._serverSideFilter.submitButton().click(function(e){
      e.preventDefault();
      self._submit();
    });
  
    this.filterIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.$el.addClass('expanded');
    });
  
    this.collapseIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.hide();
      self.$el.removeClass('expanded');
    });
  
    this.expandIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.$el.addClass('expanded');
    });
    return this;
  },  

}); // end FilterHeaderCell

///// Header Cell Filters /////

var BackgridFormFilter = Backbone.Form.extend({
  
  template: _.template([
    "<form class='form-horizontal container' >",
    '<div class="row center-block" style="margin: 0 0 0 0;" >',
    "<div data-fieldsets   />",
    "</div>",
    "</form>"].join('')),
   checkboxTemplate: _.template([
      '<div class="form-group" style="margin-bottom: 0px;" >',
      '    <div class="checkbox" style="text-align: left; min-height: 0px; padding-top: 0px;" > ',
      '      <label for="<%= editorId %>"><div><span data-editor\></div><%= title %></label>',
      '    </div>',
      '  </div>'
      ].join('')),
   
  initialize : function(options) {
    BackgridFormFilter.__super__.initialize.apply(this, arguments);
    var self = this;

    this.fieldinformation = options.fieldinformation || this.fieldinformation;
    if (_.isUndefined(this.fieldinformation)){
      throw Exception('must define a fieldinformation');
    }
    
    this.columnName = options.columnName || this.columnName;
    if (_.isUndefined(this.columnName)){
      throw Exception('must define a columnName member');
    }
    
  },

  _printSearchHash: function(searchHash){
    function lookupOperator(operator){
      var lookup = {
        gt: '>',
        lt: '<',
        gte: '>=',
        lte: '<=',
        eq: '='
      };
      return _.result(lookup, operator, operator);
    }
    
    return '&nbsp;' + _.map(              
      _.pairs(searchHash), 
      function(pair){
        var key_operator = pair[0].split('__');
        var val = '' + pair[1]
        var operator = '='
        if (key_operator.length == 2){
          operator = key_operator[1];
          operator = lookupOperator(operator);
        }
        if (key_operator[0].charAt(0)=='-'){
          operator = 'not ' + operator;
        }
        val = val.split(",").join(', ');
        return operator + '&nbsp;' + val;
      }).join('&');
  },
  
  _getSearchHash: function(){
    throw '_getSearchHash must be implemented';
  },
  
  _search: function(){
    throw '_search must be implemented';
  },
  
  clearSearch: function(){
    throw 'clearSearch must be implemented';
  },
  
  _submit: function(){
    throw '_submit must be implemented';
  },
 
  /**
   * - add a submit button 
   * - add a clear button
   */
  render: function () {
    var self = this;
    BackgridFormFilter.__super__.render.apply(this, arguments);

    this.$el.append([
      '<div id="form-last-row" class="iccbl-headerfield-form-last-row" >',
      '<div class="col-xs-6">',
      '<button type="submit" class="btn btn-default btn-block" > ok </input>',
      '</div>',
      '<div class="col-xs-6">',
      '<a class="backgrid-filter clear close" data-backgrid-action="clear"',
      ' href="#">&times;</a>',
      '</div>',
      '</div>'
      ].join(''));

    return this;
  },
  
  clearButton: function(){
    return this.$el.find("a[data-backgrid-action=clear]");
  },
  
  submitButton: function(){
    return this.$el.find(':submit');
  },
  
  /**
   * Override the Backbone Layoutmanager template rendering to use Backbone
   * Forms
   */
  renderTemplate: function() {
    return Backbone.Form.prototype.render.apply(this);
  },
  
  templateData: function() {
    return { years: 0, months: 0, dates: 0 };
  },
  
  /**
   * Convenience - determine if the form has been set with any values
   */
  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      // signal isSet for any field value set
      return values[key]>0 || !_.isEmpty(values[key]);
    });
    return !_.isEmpty(found);
  },
  
});


var CriteriumFormFilter = Iccbl.CriteriumFormFilter = BackgridFormFilter.extend({
  criterium: {'=':'eq'},
  errorClass: 'has-error',
  criteriaTemplate: 
    [
      '<span  data-editor ></span>'
    ].join(''),
  fieldTemplate: [
      '<div data-editor title="<%= help %>" class="" >',
    ].join(''),

  getPossibleSearches: function(){
    var self = this;
    var name = self.columnName;
    var possibleSearches = [];
    _.each(_.values(self.criterium), function(criteria){
      possibleSearches.push(name + '__' + criteria);
      possibleSearches.push('-'+name + '__' + criteria);
    });
    return possibleSearches;
  },
  
  /**
   * Determine if the form has been set with any values.
   */
  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      if(key == 'lower_criteria' ){
        if(values[key] == 'blank' || values[key] == 'not blank') return true;
        return false;
      }
      // signal isSet for any field value set
      return values[key]>0 || !_.isEmpty(values[key]);
    });
    return !_.isEmpty(found);
  },
  
  clear: function(){
    var self = this;
    _.each(_.keys(self.getValue()), function(key){
      self.setValue(key, null);
    });
  }
  
});

var TextFormFilter = CriteriumFormFilter.extend({
  
  criterium: {'=':'eq','contains':'contains','icontains':'icontains','<>':'ne', 
    'in': 'in','blank':'is_blank','not blank':'not_blank'},
    
  // provide a custom form template; use Bootstrap layout/styling
  template: _.template([
      '<form class="iccbl-headerfield-form" >',
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      "<div class='col-xs-12' >",
      '   <div data-fields="lower_criteria" ',
      '     class="iccbl-headerfield-text" for="lower_value"   />',
      '<div class="form-group" data-fields="form_textarea" />',
      '<div class=""  data-fields="invert_field" />',
      '</div>',
      '</div>',
      '</form>'
    ].join('')),

  initialize : function(options) {
    var self = this;
    var options = this.options = options || {};
    
    var formSchema = this.schema = {};
    formSchema['lower_criteria'] = {
        title: '', 
        key:  'lower_criteria', // TODO: "key" not needed>?
        type: 'Select',
        options: _.keys(self.criterium),
        template: _.template(self.criteriaTemplate),
        editorClass: 'form-control'
    };
    formSchema['form_textarea'] = {
        title: '',
        help: 'enter a comma separated list',
        key: 'form_textarea',
        type: 'TextArea',
        template: _.template(self.fieldTemplate),
        editorClass: 'form-control'
    };
    formSchema['invert_field'] = {
        title: 'invert',
        help: 'select this to invert the criteria',
        type: 'Checkbox',
        template: self.checkboxTemplate,
        editorClass: ''
    };    

    var FormFields = Backbone.Model.extend({
      schema: formSchema,
      validate: function(attrs) {
        var errs = {};
        if(!_.isEmpty(errs)) return errs;
      }
    });
    this.model = new FormFields();
    this.model.set('lower_criteria','='); // default
    this.selectedFields = ['lower_criteria','form_textarea','invert_field']; 
    
    TextFormFilter.__super__.initialize.apply(this, arguments);

    this.listenTo(this, "change", function(e){
      var criteria = self.getValue('lower_criteria');
      if(criteria == 'blank' || criteria == 'not blank'){
        self.$el.find('[data-fields="form_textarea"]').hide();
      }else{
        self.$el.find('[data-fields="form_textarea"]').show();
      }
    });
  },
  
  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      if(key == 'lower_criteria' ){
        if(values[key] == 'blank' || values[key] == 'not blank') return true;
        return false;
      }
      // signal isSet for any field value set
      return values[key]>0 || !_.isEmpty(values[key]);
    });
    return !_.isEmpty(found);
  },
  
  _search: function(hash){
    var self = this;
    var searchHash = _.clone(hash);
    
    var found = false;
    _.each(_.keys(self.criterium), function(criteriaKey){
      
      var criteria = self.criterium[criteriaKey];
      var searchTerm = self.columnName + '__' + criteria;
      var nsearchTerm = '-' + self.columnName + '__' + criteria;
      
      var searchVal = null;
      var negated = false;
      if(_.has(searchHash, searchTerm)){
        searchVal = searchHash[searchTerm];
      }else if(_.has(searchHash, nsearchTerm)){
        searchVal = searchHash[nsearchTerm];
        negated=true;
      }else if (_.has(searchHash,self.columnName)){
        searchVal = searchHash[self.columnName];
        criteria = 'eq';
        searchTerm = self.columnName + '__' + criteria;
        delete searchHash[self.columnName]
        searchHash[searchTerm] = searchVal;
      }else if (criteria == 'is_blank'){
        searchTerm = self.columnName + '__is_null';
        nsearchTerm = '-' + self.columnName + '__is_null';
        searchVal = _.result(
          searchHash,searchTerm, _.result(searchHash,nsearchTerm, null));
      }
      if(searchVal !== null){
        found = true;
        self.setValue('lower_criteria', criteriaKey);
        if(criteria == 'is_blank'){
          self.$el.find('[data-fields="form_textarea"]').hide();
          if(searchVal == 'false'){
            self.setValue('lower_criteria', 'not blank');
          }
        }else{
          self.$el.find('[data-fields="form_textarea"]').show();
          self.setValue('form_textarea', searchVal);
        }
        if(negated){
          self.setValue('invert_field', true);
        }
      }
    });
    return found;
  },  

  _getSearchHash: function(){
    var self = this;
    
    var searchHash = {};
    var values = self.getValue();
    var criteria = self.criterium[values['lower_criteria']];
    var searchKey = self.columnName + '__' + criteria;
    var searchVal = values['form_textarea'];
    if(criteria == 'not_blank'){
      searchKey = self.columnName + '__' + 'is_blank';
      searchVal = 'false'
    }else if(criteria == 'is_blank'){
      searchVal = 'true';
    }
    var invert = values['invert_field'];
    if(invert) searchKey = '-'+searchKey;
    searchHash[searchKey] = searchVal;
  
    return searchHash;
    
  },

  _submit: function(){
    var self = this;
    if(!self.isSet()) return;
    var self  = this;
    var searchHash = {};
    var errors = self.commit({ validate: true }); 
    if(!_.isEmpty(errors)){
      console.log('form errors, abort submit: ' + JSON.stringify(errors));
      return;
    }else{
      // this.$el.find('#range_upper_block').removeClass(self.errorClass);
    }
    
    return self._getSearchHash();
  },

  getPossibleSearches: function(){
    var possibleSearches = 
      CriteriumFormFilter.prototype.getPossibleSearches.apply(this,arguments);
    // TODO: add in the "=" (without __eq)
    possibleSearches.push(this.columnName)
    return possibleSearches;
  }

});

var DateEditor = Backbone.Form.editors.Date.extend({

  /**
   * need to extend because an error in the bbf initializer makes it impossible
   * to override the monthnames otherwise. (Using var Self = Form.editors.Date;
   * this.options = _.extend({ monthNames: Self.monthNames, why?
   */
  monthNames: [
    'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 
    'Nov', 'Dec'],

  initialize : function(options) {
    var self = this;
    DateEditor.__super__.initialize.apply(this, arguments);
    // now override monthnames
    self.options.monthNames = self.monthNames;
  },
});

var DateFormFilter = CriteriumFormFilter.extend({
  
  criterium: {'': 'unset', '=':'eq','>':'gt','>=':'gte','<':'lt','<=':'lte','<>':'ne',
    'between':'range', 'in': 'in','blank':'is_null','not blank':'not_blank'},

  template: _.template([
      '<form class="iccbl-headerfield-form" >',
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      '   <div data-fields="lower_criteria" ',
      '     class="form-control" for="lower_value"   />',
      '   <div class="input-group pull-right"  data-fields="lower_value"/>',
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      '   <div class="form-group" data-fields="form_textarea" style="display: none;" />',
      '</div>',
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      '   <div class="input-group" id="range_upper_block" style="display: none;" >',
      '     <span class="input-group-addon" for="upper_value"  style="width: 4em; ">and</span>',
      '     <span data-fields="upper_value"/>',
      '   </div>',
      '</div>',
      '   <div class=""  data-fields="invert_field" />',
      '</div>',
      '</form>'
    ].join('')),
      
  initialize : function(options) {
    var self = this;
    
    var options = this.options = options || {};
    var formSchema = this.schema = {};
    formSchema['lower_criteria'] = {
        title: '', 
        key:  'lower_criteria', // TODO: "key" not needed>?
        type: 'Select',
        options: _.keys(self.criterium),
        template: _.template(self.criteriaTemplate),
        editorClass: 'form-control'
    };
    formSchema['lower_value'] = {
        title: '',
        key: 'lower_value',
        type: DateEditor,
        template: _.template(self.fieldTemplate),
        editorClass: 'form-control',
        monthNames: self.monthNames

    };
    formSchema['form_textarea'] = {
        title: '',
        help: 'enter a comma separated list',
        key: 'form_textarea',
        type: 'TextArea',
        template: _.template(self.fieldTemplate),
        editorClass: 'form-control'
    };
    formSchema['upper_value'] = {
        title: '',
        key: 'upper_value',
        type: DateEditor,
        template: _.template(self.fieldTemplate),
        editorClass: 'form-control',
        monthNames: self.monthNames
    };
    formSchema['invert_field'] = {
        title: 'invert',
        help: 'select this to invert the criteria',
        type: 'Checkbox',
        template: self.checkboxTemplate,
        editorClass: ''
    };    

    var FormFields = Backbone.Model.extend({
      schema: formSchema,
      validate: function(attrs) {
        var errs = {};
        if(attrs.lower_criteria == 'in' 
          && !_.isEmpty(attrs.form_textarea) ){
          var datevals = attrs.form_textarea.split(',');
          var errmsgs = [];
          _.each(datevals, function(dateval){
            try{
              var v = new Date(dateval)
            }catch(e){
              errmsgs.push('not a date: '+dateval + ', err: ' + e);
            }
          })
        }
        if(!_.isEmpty(errmsgs)){
          errs['form_textarea'] = errmsgs;
          return errs;
        }
      }
    });
    this.model = new FormFields();
    this.model.set('lower_criteria','='); // default
    this.selectedFields = ['lower_criteria','lower_value','form_textarea',
                           'upper_value','invert_field']; 
    
    this.listenTo(this, "change", function(e){
      var criteria = self.getValue('lower_criteria');
      console.log('change:' + criteria)
      if(criteria == 'between'){
        self.$el.find('[data-fields="lower_value"]')
          .find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').show();
      }else if(criteria == 'in'){
        self.$el.find('[data-fields="lower_value"]')
          .find('input').prop('disabled', true);
        self.setValue('lower_value', '');
        self.$el.find('#range_upper_block').hide();
        self.$el.find('[data-fields="form_textarea"]').show();
      }else{
        self.$el.find('[data-fields="lower_value"]')
          .find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').hide();
      }
    });

    DateFormFilter.__super__.initialize.apply(this, arguments);
    
  },

  clear: function(){
    this.model.set('lower_criteria',null);
  },
  
  _search: function(hash){
    var self = this;
    var searchHash = _.clone(hash);
    
    var found = false;
    _.each(_.keys(self.criterium), function(criteriaKey){
      var criteria = self.criterium[criteriaKey];
      var searchTerm = self.columnName + '__' + criteria;
      var nsearchTerm = '-' + self.columnName + '__' + criteria;
      var searchVal = null;
      var negated = false;
      if(_.has(searchHash, searchTerm)){
        var searchVal = searchHash[searchTerm];
      }else if(_.has(searchHash, nsearchTerm)){
        var searchVal = searchHash[nsearchTerm];
        negated=true;
      }
      if(searchVal !== null){
        found = true;
        self.setValue('lower_criteria', criteriaKey);
        try{
          if(criteria == 'range'){
            self.$el.find('#range_upper_block').show();
            var vals = searchVal.split(',');
            if(vals.length < 2){
              throw "the range filter requires 2 date arguments separated by " +
              "a comma, given: " + searchVal
            }
            self.setValue('lower_value', new Date(vals[0]));
            self.setValue('upper_value', new Date(vals[1]));
          }else if(criteria == 'in'){
            self.$el.find('[data-fields="lower_value"]')
              .find('input').prop('disabled', true);
            self.setValue('lower_value', '');
            self.$el.find('[data-fields="form_textarea"]').show();
            self.setValue('form_textarea', searchVal);
          }else if(criteria == 'is_null'){
            if(searchVal == 'false'){
              self.setValue('lower_criteria', 'not blank');
            }
          }else{
            self.setValue('lower_value', new Date(searchVal));
          }
        }catch(e){
          var msg = 'Unable to parse date portion of the url, column: ' + 
              self.columnName +', searchVal:'+ searchVal + ', error: ' + e;
          console.log(msg);
          Iccbl.appModel.error(msg);
          return false;
        }
        
        if(negated){
          self.setValue('invert_field', true);
        }
      }
    });
    return found;
  },  

  _getSearchHash: function(){
    var self = this;
    var searchHash = {};
    var values = self.getValue();
    var name = self.columnName;

    
    var invert = values['invert_field'];
    if(invert) name = '-'+name;
    var criteria = self.criterium[values['lower_criteria']];
    var searchKey = name + '__' + criteria;
    
    if(criteria == 'in'){
      searchHash[searchKey] = values['form_textarea'];
    }else if(criteria == 'is_null'){
      searchHash[searchKey] = 'true';
    }else if(criteria == 'not_blank'){
      var searchKey = name + '__is_null';
      searchHash[searchKey] = 'false';
    }else if(_.isDate(values['lower_value']) ){
      if(criteria == 'range'){
        if(_.isDate(values['upper_value'])){
          var searchKey = name + '__' + criteria;
          searchHash[searchKey] = values['lower_value'].toISOString()
              + ',' + values['upper_value'].toISOString();
        }else{
          console.log('upper value not set; validation should have caught this');
        }
      }else{
        searchHash[searchKey] = values['lower_value'].toISOString();
      }
    }
    return searchHash
    
  },

  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;
    
    // validate:true: tells bbf to run model.validate(), in addition to
    // field[].validate()
    var errors = self.commit({ validate: true }); 
    if(!_.isEmpty(errors)){
      console.log('form errors, abort submit: ' + JSON.stringify(errors));
      this.$el.find('#range_upper_block').addClass(self.errorClass);
      return;
    }else{
      this.$el.find('#range_upper_block').removeClass(self.errorClass);
    }
    
    return self._getSearchHash();
  }
});

var BooleanFormFilter = CriteriumFormFilter.extend({
  criterium: {'': 'unset', 'true': 'true', 'false': 'false', 
    'blank':'is_null','not blank':'not_blank'},
    
  template: _.template([
      '<form class="iccbl-headerfield-form" >',
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      "<div class='col-xs-12' >",
      '   <div data-fields="lower_criteria" ',
      '     class="iccbl-headerfield-text" for="lower_value"   />',
      "</div>",
      "</form>"].join('')),
  
  initialize : function(options) {
    var self = this;
    var options = this.options = options || {};
    var formSchema = this.schema = {};
    formSchema['lower_criteria'] = {
        title: '', 
        key:  'lower_criteria', // TODO: "key" not needed>?
        type: 'Select',
        options: _.keys(self.criterium),
        template: _.template(self.criteriaTemplate),
        editorClass: 'form-control'
    };
    
    var FormFields = Backbone.Model.extend({
      schema: formSchema,
      validate: function(attrs) {
        var errs = {};
        if(!_.isEmpty(errs)) return errs;
      }
    });
    this.model = new FormFields();
    this.selectedFields = ['lower_criteria'] 
    
    BooleanFormFilter.__super__.initialize.apply(this, arguments);
  },

  _getSearchHash: function(){
    var self = this;
    var searchHash = {};
    var values = self.getValue();
    var name = self.columnName;

    var criteria = self.criterium[values['lower_criteria']];
    if(criteria == 'not_blank'){
      searchKey = name + '__' + 'is_null';
      searchHash[searchKey]='false';
    }else if(criteria == 'is_null'){
      searchKey = name + '__' + 'is_null';
      searchHash[searchKey]='true';
    }else if(criteria == 'false'){
      searchKey = name + '__eq';
      searchHash[searchKey]='false';
    }else if(criteria == 'true'){
      searchKey = name + '__eq';
      searchHash[searchKey]='true';
    }
    return searchHash;
  },
  
  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;

    return self._getSearchHash();
  },

  _search: function(hash){
    var self = this;
    var searchHash = _.clone(hash);
    
    var searchTerm = null;
    _.each(self.getPossibleSearches(), function(term){
      if(_.has(searchHash,term)) searchTerm = term;
    });
    var searchVal = searchHash[searchTerm];
    var name = this.columnName;
    if(searchTerm){
      if(searchTerm ==  name + '__is_null'){
        if(searchVal == 'true'){
          self.setValue('lower_criteria','blank');
        }else{
          self.setValue('lower_criteria', 'not blank');
        }
      }else if(searchTerm ==  name + '__eq'){
        if(searchVal == 'true'){
          self.setValue('lower_criteria','true');
        }else{
          self.setValue('lower_criteria', 'false');
        }
      }
    }
    return searchTerm;
  },  
  
  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      return values[key];
    });
    return found;
  },
  
  getPossibleSearches: function(){
    return [this.columnName + '__eq',
            this.columnName + '__is_null'];
  }
  
});

var SelectorFormFilter = CriteriumFormFilter.extend({

  criterium: {'': 'unset', 'blank':'is_null','not blank':'not_blank'},
  
  template: _.template([
      "<form  class='form-horizontal container ' >",
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      "<div data-fieldsets   />",
      "</div>",
      "</form>"].join('')),

   altFieldTemplate: _.template([
      '<div class="form-group" style="margin-bottom: 0px;" >',
      '    <div class="checkbox" style="text-align: left; ',
      '       min-height: 0px; padding-top: 0px;" > ',
      '      <label for="<%= editorId %>">',
      '      <div><span data-editor\></div><%= title %></label>',
      '    </div>',
      '  </div>'
      ].join('')),
 
  initialize: function(options){
    
    var self = this;
    this.fieldinformation = options.fieldinformation || this.fieldinformation;
    if (_.isUndefined(this.fieldinformation)){
      throw Exception('must define a fieldinformation');
    }
    
    var choiceHash = {}
    var vocabulary;
    if(_.isUndefined(this.fieldinformation.choices)){
      if (Iccbl.appModel.DEBUG)
        console.log([
            'Warn: fieldinformation for a selection field type must define a ',
            '"choices" list: field key: ' + this.column.get('name')].join(''));
      this.fieldinformation.choices = [];
    }
    if(!_.isEmpty(this.fieldinformation.vocabulary)){
      // TODO: vocabulary is using the titles as the key, 
      // because of how Backgrid.SelectCell initializes
      _.each(this.fieldinformation.vocabulary,function(pair){
        choiceHash[pair[1]] = pair[0];
      });
    }else{
      try{
        vocabulary = Iccbl.appModel.getVocabulary(
          this.fieldinformation.vocabulary_scope_ref);
        _.each(_.keys(vocabulary),function(choice){
          choiceHash[choice] = vocabulary[choice].title;
        });
      }catch(e){
        console.log('vocabulary error', this.column.get('name'),e);
      }
    }

    var formSchema = this.schema = {};  

    if(_.isUndefined(this.fieldinformation.choices)){
      if (Iccbl.appModel.DEBUG)
        console.log([
            'Warn: fieldinformation for a selection field type must define a ',
            '"choices" list: field key: ' + this.column.get('name')].join(''));
      this.fieldinformation.choices = [];
    }
    var selectedFields = _.clone(this.fieldinformation.choices);
    
    _.each(_.keys(choiceHash), function(choice){
      formSchema[choice] = { 
          title: choiceHash[choice],
          key:  choice, 
          type: 'Checkbox',
          template: self.altFieldTemplate };
    });

    formSchema['lower_criteria'] = {
      title: '', 
      key:  'lower_criteria', // TODO: "key" not needed>?
      type: 'Select',
      options: _.keys(self.criterium),
      template: _.template(self.criteriaTemplate),
      editorClass: 'form-control'
    };
    selectedFields.push('lower_criteria');
    
    formSchema['invert_field'] = {
        title: 'invert',
        help: 'select this to invert the criteria',
        type: 'Checkbox',
        template: self.checkboxTemplate,
        editorClass: ''
    };    
    selectedFields.push('invert_field');
    
    var FormFields = Backbone.Model.extend({
      schema: formSchema,
      validate: function(attrs) {
        var errs = {};
        if(!_.isEmpty(errs)) return errs;
      }
    });
    this.model = new FormFields();
    this.selectedFields = selectedFields; 
    
    SelectorFormFilter.__super__.initialize.apply(this, arguments);

    this.listenTo(this, "change", function(e){
      var criteria = self.getValue('lower_criteria');
      console.log('criteria: ' + criteria);
      if(criteria == 'blank'){
        self.$el.find('[data-fields]').find('input').prop('disabled', true);
      }else if(criteria == 'not blank'){
        self.$el.find('[data-fields]').find('input').prop('disabled', true);
      }else {
        self.$el.find('[data-fields]').find('input').prop('disabled', false);
      }
    });
    
  },
      
  clear: function(){
    SelectorFormFilter.__super__.clear.apply(this, arguments);
    this.$el.find('[data-fields]').find('input').prop('disabled', false);
  },

  _getSearchHash: function(){
    var self = this;
    var searchHash = {};
    var values = self.getValue();
    var name = self.columnName;

    var invert = values['invert_field'];
    if(invert) name = '-'+name;
    
    var criteria = self.criterium[values['lower_criteria']];
    var searchKey = name + '__' + criteria;
    if(criteria == 'not_blank'){
      searchKey = name + '__' + 'is_null';
      searchVal = 'false';
      searchHash[searchKey]=searchVal;
    }else if(criteria == 'is_null'){
      searchVal = 'true';
      searchHash[searchKey]=searchVal;
    }else{
      var selected = _.filter(_.keys(values), function(key){ 
        if(key !== 'invert_field') return values[key]; 
        return false;
      });
      searchHash[name +'__in'] = selected.join(',');
    }
    
    return searchHash;
  },
    
  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;

    return self._getSearchHash();
  },

  _search: function(hash){
    var self = this;
    var searchHash = _.clone(hash);
    var searchTerm = null;
    var name = this.columnName;
    var searchVal = null;

    _.each(self.getPossibleSearches(), function(term){
      if(_.has(searchHash,term)){
        searchTerm = term;
        searchVal = searchHash[term];
      }
    });
    
    if(searchTerm){
      if(searchTerm.charAt(0) == '-'){
        self.setValue('invert_field', true);
        name = '-' + name;
      }
      if(searchTerm ==  name + '__is_null'){
        if(searchVal == 'true' || searchVal == true){
          self.setValue('lower_criteria','blank');
        }else{
          self.setValue('lower_criteria', 'not blank');
        }
      }else{
        var searchVal = searchHash[searchTerm];
        if (!_.isArray(searchVal)){
          searchVal = searchVal.split(',');
        }
        _.each(searchVal, function(choice){
          self.setValue(choice, true);
        });
      }
    }
    return searchTerm;
  },  
  
  /**
   * SelectorFormFilter Convenience - determine if the form has been set with
   * any values
   */
  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      if(key == 'invert_field' ) return false; // skip invert field 
      return values[key];
    });
    return found;
  },
  
  getPossibleSearches: function(){
    return [this.columnName + '__in', '-'+this.columnName + '__in',
            this.columnName + '__is_null'];
  },
});

var NumberFormFilter = CriteriumFormFilter.extend({
  
  criterium: {
    '=':'eq','\u2248':'about','>':'gt', '>=':'gte','<':'lt','<=':'lte',
    '<>':'ne', 'x..y':'range', 'in': 'in','blank':'is_null',
    'not blank':'not_blank'
  },
  
  // FIXME: this template is a mix of adaptive and fixed styles, for instance:
  // for the first input-group: style="width: 50px"
  template: _.template([
      '<form class="iccbl-headerfield-form" >',
      '<div class="center-block" style="margin: 0 0 0 0;" >',
      '<div class="input-group" style="width: 50px" >',
      '   <div data-fields="lower_criteria" ',
      '     class="input-group-addon iccbl-headerfield-number" for="lower_value"   />',
      '   <div data-fields="lower_value"/>',
      '</div>',
      '<div class="form-group" data-fields="form_textarea" style="display: none;" />',
      '<div class="input-group" id="range_upper_block" style="display: none;" >',
      '   <span class="input-group-addon" for="upper_value" style="width: 4em; ">to</span>',
      '   <span data-fields="upper_value"/>',
      '</div>',
      '<div class=""  data-fields="invert_field" />',
      '</div>',
      '</form>'
    ].join('')),
      
  initialize : function(options) {
    var self = this;
    var options = this.options = options || {};
    var formSchema = options.schema = options.schema || {};
    var fields = options.fields = options.fields || [];

    formSchema['lower_criteria'] = {
        title: '', 
        key:  'lower_criteria',
        type: 'Select',
        options: _.keys(self.criterium),
        template: _.template(self.criteriaTemplate),
        editorClass: 'form-control'
    };
    formSchema['lower_value'] = {
        title: '',
        key: 'lower_value',
        type: 'Number',
        template: _.template(self.fieldTemplate),
        editorClass: 'form-control'
    };
    formSchema['form_textarea'] = {
        title: '',
        help: 'enter a comma separated list',
        key: 'form_textarea',
        type: 'TextArea',
        template: _.template(self.fieldTemplate),
        editorClass: 'form-control'
    };
    formSchema['upper_value'] = {
        title: '',
        key: 'upper_value',
        type: 'Number',
        template: _.template(self.fieldTemplate),
        editorClass: 'form-control '
    };
    formSchema['invert_field'] = {
        title: 'invert',
        help: 'select this to invert the criteria',
        type: 'Checkbox',
        template: self.checkboxTemplate,
        editorClass: ''
    };    

    var FormFields = Backbone.Model.extend({
      validate: function(attrs) {
        var errs = {};
        //        if(attrs.lower_criteria == '...' 
        //          && ( attrs.upper_value < 1 ) ){
        //          errs['upper_value'] = '!'
        //        }
        if(!_.isEmpty(errs)) return errs;
      }
    });
    this.model = options['model'] = new FormFields();
    this.model.set('lower_criteria','>'); // default
    
    options.fields = fields.concat(
        ['lower_criteria','lower_value','form_textarea',
         'upper_value','invert_field']); 
    
    this.listenTo(this, "change", function(e){
      var criteria = self.getValue('lower_criteria');
      if(criteria == 'x..y'){
        self.$el.find('[data-fields="lower_value"]')
          .find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').show();
      }else if(criteria == 'in'){
        self.$el.find('[data-fields="lower_value"]')
          .find('input').prop('disabled', true);
        self.setValue('lower_value', '');
        self.$el.find('#range_upper_block').hide();
        self.$el.find('[data-fields="form_textarea"]').show();
      }else{
        self.$el.find('[data-fields="lower_value"]')
          .find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').hide();
      }
    });

    NumberFormFilter.__super__.initialize.call(this, options);
    
  },

  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      return values[key] !== '';
    });
    return !_.isEmpty(found);
  },
  
  _search: function(hash){
    var self = this;
    var searchHash = _.clone(hash);
    var found = false;

    _.each(_.keys(self.criterium), function(criteriaKey){
      var criteria = self.criterium[criteriaKey];
      var searchTerm = self.columnName + '__' + criteria;
      var nsearchTerm = '-' + self.columnName + '__' + criteria;
      var searchVal = null;
      var negated = false;
      if(_.has(searchHash, searchTerm)){
        var searchVal = searchHash[searchTerm];
      }else if(_.has(searchHash, nsearchTerm)){
        var searchVal = searchHash[nsearchTerm];
        negated=true;
      }
      if(searchVal !== null){
        found = true;
        self.setValue('lower_criteria', criteriaKey);
        if(criteria == 'range'){
          self.$el.find('#range_upper_block').show();
          var vals = searchVal.split(',');
          self.setValue('lower_value', vals[0]);
          self.setValue('upper_value', vals[1]);
        }else if(criteria == 'in'){
          self.$el.find('[data-fields="lower_value"]')
            .find('input').prop('disabled', true);
          self.setValue('lower_value', '');
          self.$el.find('[data-fields="form_textarea"]').show();
          self.setValue('form_textarea', searchVal);
        }else if(criteria == 'is_null'){
          if(searchVal == 'false'){
            self.setValue('lower_criteria', 'not blank');
          }
        }else{
          self.setValue('lower_value', searchVal);
        }
        
        if(negated){
          self.setValue('invert_field', true);
        }
      }
    });
    return found;
  },  

  _getSearchHash: function(){
    var self = this;
    var searchHash = {};
    var values = self.getValue();
    var name = self.columnName;

    var searchHash = {};
    var invert = values['invert_field'];
    if(invert) name = '-'+name;
    var criteria = self.criterium[values['lower_criteria']];
    var searchKey = name + '__' + criteria;
    
    if(criteria == 'in'){
      searchHash[searchKey] = values['form_textarea'];
    }else if(criteria == 'is_null'){
      searchHash[searchKey] = 'true';
    }else if(criteria == 'not_blank'){
      var searchKey = name + '__is_null';
      searchHash[searchKey] = 'false';
    }else if(''+values['lower_value'] !== ''){
      if(criteria == 'range'){
        if(''+values['upper_value'] !== ''){
          searchKey = name + '__' + criteria;
          searchHash[searchKey] = values['lower_value'] + ',' + values['upper_value'];
        }else{
          console.log('upper value not set; validation should have caught this');
        }
      }else{
        searchHash[searchKey] = ''+values['lower_value'];
      }
    }
    return searchHash;
  },

  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;
    var errors = self.commit({ validate: true }); 
    if(!_.isEmpty(errors)){
      console.log('form errors, abort submit: ' + JSON.stringify(errors));
      this.$el.find('#range_upper_block').addClass(self.errorClass);
      return;
    }else{
      this.$el.find('#range_upper_block').removeClass(self.errorClass);
    }
    
    return self._getSearchHash();
  }
});

var SIUnitFormFilter = NumberFormFilter.extend({
  
  symbol: "",
  
  // provide a custom form template; use Bootstrap layout/styling
  template: _.template([
      '<form class="iccbl-headerfield-form" >',
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      '<div class="input-group  col-sm-2">',
      '   <div data-fields="lower_criteria" ',
      '     class="input-group-addon iccbl-headerfield-number" for="lower_value"   />',
      '   <div data-fields="lower_value"/>',
      '   <div class="input-group-addon iccbl-headerfield-number" data-fields="lower_siunit"/>',
      '</div>',
      '<div class="form-group" data-fields="form_textarea" style="display: none;" />',
      '<div class="input-group" id="range_upper_block" style="display: none;" >',
      '   <span class="input-group-addon" for="upper_value"  style="width: 4em; ">to</span>',
      '   <span data-fields="upper_value"/>',
      '   <div class="input-group-addon iccbl-headerfield-number" data-fields="upper_siunit"/>',
      '</div>',
      '</div>',
      '<div class=""  data-fields="invert_field" />',
      '</form>'
    ].join('')),

  siunits: [
      ['T', 1e12],
      ['G', 1e9],
      ['M', 1e6],
      ['k', 1e3],
      ['', 1],
      ['m', 1e-3,],
      ['μ', 1e-6,],
      ['n', 1e-9 ],
      ['p', 1e-12 ]
  ],

  initialize : function(options) {
    var self = this;

    this.fieldinformation = options.fieldinformation || this.fieldinformation;
    if (_.isUndefined(this.fieldinformation)){
      throw Exception('must define a fieldinformation');
    }
    
    var options = _.extend({},this.fieldinformation['display_options'],options);

    if(!options.symbol)
    {
      throw new Error('SIUnitHeaderCell: field information requires the '+
          '"symbol" option');
    }
    var multiplier = this.multiplier = options.multiplier || 1;
    var symbol = this.symbol = options.symbol;
    var defaultUnit = this.defaultUnit = options.defaultUnit;
    var units = this.units = [];
    var formSchema = options.schema = options.schema || {};
    
    if(! options.symbol){
      throw 'Error: SIUnitFormFilter requires a "symbol" option' 
    }
    this.defaultSymbol = null;
    _.each(this.siunits,function(pair){
      if(options.maxunit){
        if(options.maxunit < pair[1]) return;
      }
      if(options.minunit){
        if(options.minunit > pair[1]) return;
      }
      units.push({ val: pair[1], label: pair[0] + self.symbol });
      if (pair[1]==defaultUnit){
        self.defaultSymbol = pair[0] + self.symbol;
      }
    });
    formSchema['lower_siunit'] = {
      title: '', 
      key:  'lower_siunit', // TODO: "key" not needed>?
      type: 'Select',
      options: units,
      template: _.template(self.criteriaTemplate),
      editorClass: 'form-control'
    };
    
    formSchema['upper_siunit'] = {
      title: '', 
      key:  'upper_siunit', // TODO: "key" not needed>?
      type: 'Select',
      options: units,
      template: _.template(self.criteriaTemplate),
      editorClass: 'form-control'
    };
    
    options['fields'] = ['lower_siunit','upper_siunit']
    
    SIUnitFormFilter.__super__.initialize.call(this, options);
    
  },
  
  render: function(){
    SIUnitFormFilter.__super__.render.call(this, arguments);
    
    // Fixme: these values must be set after render, because inheritance is not
    // proper for this class.
    if (_.isNumber(this.defaultUnit)){
      this.setValue('lower_siunit',this.defaultUnit);
      this.setValue('upper_siunit',this.defaultUnit);
    }
    return this;
    
  },
  
  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      if(key == 'lower_criteria' 
        || key == 'lower_siunit'
        || key == 'upper_siunit' ) return false;
      // signal isSet for any field value set
      return values[key] !== '';
    });
    return !_.isEmpty(found);
  },
  
  _printSearchHash: function(searchHash){
    var self = this;
    _.each(_.keys(searchHash), function(key){
      var nu = self._findNumberAndUnit(searchHash[key]);
      searchHash[key] = nu.number + '&nbsp;' 
        + _.result(_.invert(_.object(self.siunits)),nu.unit,nu.unit)
        + self.symbol;
    });
    return SIUnitFormFilter.__super__._printSearchHash.call(this, searchHash);
  },
  
  _getSearchHash: function(){
    var self = this;
    var searchHash = {};
    var values = self.getValue();
    var name = self.columnName;

    var searchHash = SIUnitFormFilter.__super__._getSearchHash.call(this);
    if(!_.isEmpty(searchHash)){
      var searchKey = _.keys(searchHash)[0];
      var searchValue = searchHash[searchKey];
      var values = self.getValue();
      var criteria = self.criterium[values['lower_criteria']];

      if(criteria == 'range'){
        searchHash[searchKey] = [
            self._calculate(
              self.multiplier,values['lower_siunit'],values['lower_value']),
            self._calculate(
              self.multiplier,values['upper_siunit'],values['upper_value'])
            ].join(',');
      }else if(criteria == 'in'){
        var newvalues = _.map(searchValue.split(','),function(val){
          return ''+self._calculate(self.multiplier,values['lower_siunit'],val);
        });
        searchHash[searchKey] = newvalues.join(',');
      }else if(criteria == 'is_null'){
        // do nothing it's good already
      }else if(criteria == 'not_blank'){
        // do nothing it's good already
      }else{
        searchHash[searchKey] = ''+self._calculate(
            self.multiplier,values['lower_siunit'],values['lower_value']);
      }
      console.log('SIunit new search value: ' + searchHash[searchKey]);
    }
    return searchHash;
  },

  
  /**
   * SIUnitFormFilter Form submit handler
   */
  _submit: function(){
    var self = this;
    SIUnitFormFilter.__super__._submit.call(this);
    return self._getSearchHash();
  },
  
  _calculate: function(multiplier, sci_mult, val){
    // Run strip after every calculation to round out floating point math errors
    function strip(number) {
      return (parseFloat(number.toPrecision(12)));
      };
    val = strip(val / multiplier);
    if(sci_mult > 0){ // if sci unit is undefined, assume to be 1
      val = strip(val * sci_mult);
    }
    return val;
  },
  
  _findNumberAndUnit: function(number){
    var decimals = 3; // TODO: users will not be expected to enter values beyond
                      // 3 decimals
    var self = this;
    function strip(number) {
      return (parseFloat(number.toPrecision(12)));
      };
    number = strip(number/self.multiplier);
    pair = _.find(this.siunits, function(pair){
      return pair[1] <= Math.abs(number); 
    });
    
    if(_.isUndefined(pair)){
      console.log('could not find units for the input number: ' + number);
      return { number:number, unit: ''};
    }
    
    var val = (1/pair[1])*number;
    val = Math.round(val*Math.pow(10,decimals))/Math.pow(10,decimals);
    
    return {number:val, unit: pair[1]};
  },
  
  _search: function(hash){
    var self = this;
    var searchHash = _.clone(hash);
    var found = SIUnitFormFilter.__super__._search.call(this, hash);
    if(found){
      var values = self.getValue();
      
      if(values['lower_value'] !== ''){
        var numberAndUnit = self._findNumberAndUnit(values['lower_value']);
        self.setValue('lower_value', numberAndUnit.number);
        self.setValue('lower_siunit', numberAndUnit.unit);
      }
      if(values['upper_value'] !== ''){
        var numberAndUnit = self._findNumberAndUnit(values['upper_value']);
        self.setValue('upper_value', numberAndUnit.number);
        self.setValue('upper_siunit', numberAndUnit.unit);
      }
    }
    return found;
  }  
  
});

/**
 * Return an array for backgrid column descriptors.
 * 
 * @param {Object}
 *          prop - hash of field properties from REST metadata: field properties {
 *          visibility: [array of strings], title: a label for the field, order:
 *          display order of the field , data_type: determines the type of header
 *          cell }
 * @param {Object}
 *          optionalHeaderCell - a Backgrid.HeaderCell to use for each column
 * @param {string}
 *          key - field information key in the metahash
 * @param {array}
 *          orderStack - for rendering ordered columns
 */
var createBackgridColumn = Iccbl.createBackgridColumn = 
  function(key, prop, _orderStack, optionalHeaderCell ){
  var orderStack = _orderStack || [];
  var column = {};
  var visible = _.has(prop, 'visibility') && 
                    _.contains(prop['visibility'], 'l');
  var data_type = _.isEmpty(prop.data_type)?'string':prop.data_type.toLowerCase();
  var display_type = _.isEmpty(prop.display_type)?data_type:prop.display_type.toLowerCase();
  var cell_options = prop.display_options || {};
  var edit_type = _.isEmpty(prop.edit_type)?display_type:prop.edit_type.toLowerCase();
  var backgridCellType = StringCell;
  var typeMap = {
    'date': Iccbl.DateCell,
    'link': Iccbl.LinkCell,
    'siunit': Iccbl.SIUnitsCell,
    'float': Iccbl.NumberCell, //'Number',
    'decimal': Iccbl.DecimalCell,
    'image': Iccbl.ImageCell,
    'boolean': Iccbl.BooleanCell
  }
  
  if (_.has(prop, 'backgridCellType')){
    console.log('using specified "backgridCellType": ',key, prop.backgridCellType );
    backgridCellType = prop.backgridCellType;
  } else {
    if(_.has(typeMap,display_type)){
      if (Iccbl.appModel.DEBUG)
        console.log('field', key, display_type, 'typemap',typeMap[display_type])
      var backgridCellType = typeMap[display_type];
      
      if (display_type=='link' && data_type=='list'){
        backgridCellType = Iccbl.UriListCell;
      }
    }else{
      if (Iccbl.appModel.DEBUG)
        console.log('no special cell type for', display_type, 'data_type',data_type);
    }
  }
  if(!_.isEmpty(cell_options)){
    backgridCellType = backgridCellType.extend(cell_options);
  }
  
  if (data_type == 'list'){
    backgridCellType = backgridCellType.extend({
      formatter: Iccbl.StringFormatter
    });
  }
  column = _.extend(column, {
    'name' : key,
    'label' : prop['title'],
    'description' : prop['description'],
    'mouseover' : _.result(prop, 'mouseover', prop['description']),
    'cell' : backgridCellType,
    'order' : prop['ordinal'],
    'sortable': prop['ordering'],
    'searchable': prop['filtering'],
    'editable' : false,
    'visible': visible,
    'fieldinformation': prop
  });
  if(_.has(prop,'editability') && _.contains(prop['editability'],'l')){
    column['editable'] = true;
  }
  if(!_.isEmpty(prop.vocabulary)){
    // Note: this is probably backwards, since it appears SelectCell
    // wants [title,val], not [val,title]
    cell_options.optionValues = prop.vocabulary; 
    column['cell'] = Iccbl.SelectCell.extend(cell_options);
  }else if(!_.isEmpty(prop.vocabulary_scope_ref)){
    var optionValues = [];
    var vocabulary_scope_ref = prop.vocabulary_scope_ref;
    try{
      var vocabulary = Iccbl.appModel.getVocabulary(vocabulary_scope_ref);
        _.each(_.keys(vocabulary),function(choice){
          optionValues.push([vocabulary[choice].title,choice]);
        });
    }catch(e){
      console.log('build column errorr: e',e);
    }
    cell_options = _.extend(cell_options,{
      optionValues: optionValues,
      vocabulary_scope_ref: vocabulary_scope_ref
    } );
    column['cell'] = Iccbl.SelectCell.extend(cell_options);
  }
  if(orderStack && _.contains(orderStack, key)){
    column['direction'] = 'ascending';
  }
  else if(orderStack && _.contains(orderStack, '-' + key)){
    column['direction'] = 'descending';
  }
  
  var headerCellDefaults = {
    'fieldinformation': prop,
    'serverSideFilter': TextFormFilter
  };
  if (optionalHeaderCell) {
    column['headerCell'] = optionalHeaderCell.extend(headerCellDefaults);
  }else if (_.has(prop, 'headerCell')){
    column['headerCell'] = prop.headerCell.extend(headerCellDefaults);
  } 
  else{
    // Set up a more specific header cell, with filter
    if(data_type == 'string'){
      headerCellDefaults['serverSideFilter'] = TextFormFilter;
    }
    else if(data_type == 'integer'
      || data_type == 'float'
      || data_type == 'decimal' ){
      
      if(display_type == 'siunit'){
        headerCellDefaults['serverSideFilter'] = SIUnitFormFilter;
      } else {
        headerCellDefaults['serverSideFilter'] = NumberFormFilter;
      }
    }
    else if( data_type == 'date'){
      headerCellDefaults['serverSideFilter'] = DateFormFilter;
    }
    else if( data_type == 'boolean'){
      headerCellDefaults['serverSideFilter'] = BooleanFormFilter;
    }

    if( edit_type == 'select' 
      || edit_type == 'multiselect' 
      || edit_type == 'multiselect2' 
      || edit_type == 'multiselect3' 
        
    ){
      headerCellDefaults['serverSideFilter'] = SelectorFormFilter;
    }
    column['headerCell'] = FilterHeaderCell.extend(headerCellDefaults);
  }
    
  return column;
};


var createBackgridColModel = Iccbl.createBackgridColModel = 
  function(restFields, _orderStack, _searchHash, _manualIncludes) {
    
  console.log('--createBackgridColModel');
  var manualIncludes = _manualIncludes || [];
  var orderStack = _orderStack || [];
  var searchHash = _searchHash || {};

  var colModel = [];
  var i = 0;
  var _total_count = 0;
  _.each(_.pairs(restFields), function(pair) {
    var key = pair[0];
    var prop = pair[1];
    var column = createBackgridColumn(key, prop, orderStack);
    column.key = key;
    var visible = _.has(prop, 'visibility') && 
      _.contains(prop['visibility'], 'l');
    if (visible || _.contains(manualIncludes, key) ) {
      if(_.contains(manualIncludes, '-'+key)){
        console.log('Column: ' + key + ' is manually excluded');
      }else{
        colModel[i] = column;
        i++;
      }
    } else {
      var hashSearch = RegExp('^(' + key + ')(_{2}\w+)?$');
      var orderSearch = RegExp('^-?' + key + '$');
      if( 
        _.findKey(searchHash, function(val,hashkey){
          return hashSearch.test(hashkey);
        })
        ||  _.find(orderStack, function(orderkey){
          return orderSearch.test(orderkey);
        }))
      {
        colModel[i] = column;
        i++;
      }
    }
  });
  
  colModel = new Backgrid.Columns(colModel);
  colModel.comparator = 'order';
  colModel.sort();
  return colModel;
};


return Iccbl;
});
