/**
 * @summary ICCBL-lims Backgrid Extension Functions
 * @description Utility Functions for the iccbl-lims
 * @version 0.1
 * @file iccbl-backgrid.js
 * @author Sean Erickson
 * @contact sean_erickson “AT” hms.harvard.edu
 * 
 * @copyright 2014 Harvard University, all rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * This source file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the license files for details.
 * 
 * For details please refer to: https://github.com/hmsiccbl/lims
 */

// FIXME: lunr should not be a requirement - for server side filtering,
// and backgrid_filter is requiring it.
// Note: backbone_forms included for testing imports to work properly
// - bad imports can result in the:
// "Warning: PhantomJS timed out, possibly due to a missing Mocha run() call"

define(['jquery', 'underscore', 'backbone', 'backgrid','backbone_forms', 
        'backgrid_filter', 'backgrid_paginator', 'lunr', // 'backgrid_select_all',
        'layoutmanager'],
    function($, _, Backbone, Backgrid, BackboneForms,
             BackgridFilter, BackgridPaginator, lunr, // BackgridSelectAll,
             layoutmanager ) {

var root = window;
  
var Iccbl = root.Iccbl = {
    VERSION : "0.0.1",
    appModel : "This value will be initialized on app start"
};

// TODO: remove this
var assertIccbl = Iccbl.assert = function(condition, message) {
    if (!condition) {
        throw message || "Assertion failed";
    }
};

// TODO: deprecated
requireOptions = Iccbl.requireOptions = function(options,requireOptionKeys){
    for (var i = 0; i < requireOptionKeys.length; i++) {
      var key = requireOptionKeys[i];
      if (_.isUndefined(options[key])) {
        throw new TypeError("'" + key  + "' is required");
      }
    }
};

/**
 * Replace all {tokens} in the string with model attributes. - fallback to
 * "token" if model attribute is not set.
 * @param - defaul_val - value to use if the matched token is not found in the model
 * - this can be used to replace any token with a given default value
 */
var replaceTokens = Iccbl.replaceTokens = function(model,stringWithTokens, default_val) {
  var interpolatedString = stringWithTokens.replace(/{([^}]+)}/g, 
    function (match) 
    {
      match = match.replace(/[{}]/g,'');
      if(!_.isUndefined(model.get(match))){
        return model.get(match);
      }else{
        if(!_.isUndefined(default_val)){
          return default_val;
        }else{
          return match;
        }
      }
    });
  return interpolatedString;
}

/**
 * Format a string of containing "replacement fields" surrounded by 
 * curly braces '{}'. Replacment fields are used as keys to lookup the 
 * replacement values in the valueHash.
 * If a replacement field key is not found in the valueHash, and a defaul_val
 * is provided, the default_val will replace the field, otherwise, the 
 * replacement field key itself is used.
 */
var formatString = Iccbl.formatString = function(stringWithTokens, valueHash, default_val) {
  var interpolatedString = stringWithTokens.replace(/{([^}]+)}/g, 
    function (match) 
    {
      match = match.replace(/[{}]/g,'');
      if(_.has(valueHash, match)){
        return valueHash[match];
      }else{
        if(!_.isUndefined(default_val)){
          return default_val;
        }else{
          return match;
        }
      }
    });
  return interpolatedString;
}

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
 * Returns the ISO string date part of the Date object, ignoring the timezone.
 * - for internal representation of dates and communicating dates to the server.
 * @param jsDate a JavaScript Date object
 */
var getISODateString = Iccbl.getISODateString = function(jsDate){
  return jsDate ? jsDate.toISOString().split('T')[0] : jsDate;
  // equivalent:
  //  date = lpad(jsDate.getUTCFullYear(), 4, 0) 
  //    + '-' + lpad(jsDate.getUTCMonth() + 1, 2, 0) 
  //    + '-' + lpad(jsDate.getUTCDate(), 2, 0);
  //  return date;
}

/**
 * Returns the "ICCBL" formatted date part of the Date object, ignoring the timezone.
 * "ICCBL" format is "MM/DD/YYYY"
 * - for display of dates in the UI
 * @param jsDate a JavaScript Date object
 */
var getIccblUTCDateString = Iccbl.getUTCDateString = function(jsDate){
  if (!jsDate) return jsDate;
  return ( 
      lpad(jsDate.getUTCMonth() + 1, 2, 0) 
      + '/' + lpad(jsDate.getUTCDate(), 2, 0) 
      + '/' + lpad(jsDate.getUTCFullYear(), 4, 0) );
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

/**
 * Sort the list of field keys using the associated ordinal for the key in the 
 * fieldHash
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
 * convert the app model "current_options" object into a route fragment "key"
 */
var getKey = Iccbl.getKey = function(options) {
    var route_fragment = '';
    
    if (_.isString(options)) {
        route_fragment += options;
        
    // an array is just a set of keys, to be separated by slashes
    } else if (_.isArray(options)) {
        route_fragment = _.reduce(options, function(route, option) {
                if (!_.isNull(option)) {
                    if (!_.isEmpty(route))
                        route += '/';
                    route += option;
                }
                return route;
            }, route_fragment);
    } else if (_.isObject(options)) {// generic, option order not defined
        if (_.has(options, 'key')) {
            return Iccbl.getKey(options.key);
        }
    }

    return route_fragment;
};

/**
 * Get the key off the URI stack: - we don't know if the key is composite or
 * not; so we don't know how many items to pop off the stack; So look at the
 * resource definition id_attribute; which lists the keys.
 * 
 * @param resource -
 *          a resource definition as defined by the API
 * @param urlStack -
 *          array representation of the current unprocessed URI elements.
 * @param consumedStack -
 *          holds the items popped off the stack
 */
var popKeyFromStack = Iccbl.popKeyFromStack = function(resource, urlStack, consumedStack){
  var id  = '';
  var self = this;
  var checkStack = function(stack){
    if (_.isEmpty(urlStack)){
      var msg = 'not enough items on the URL to create the key for resource: ' + 
          resource.title;
      throw msg;
    }
  };
  if(resource.key == 'apilog'){
    ref_resource_name = urlStack.shift();
    consumedStack.push(ref_resource_name);
    ref_resource = Iccbl.appModel.getResource(ref_resource_name);
    if(ref_resource && !_.isEmpty(ref_resource)){
      checkStack(urlStack);
      key = self.popKeyFromStack(ref_resource,urlStack,consumedStack);
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

var getIdKeys = Iccbl.getIdKeys = function(model,schema) {
  if (_.has(schema['resource_definition'], 'id_attribute')) {
    var id_attribute = schema['resource_definition']['id_attribute'];
    console.log('create id from ' + id_attribute);
    var idList = [];
    _.each(id_attribute, function(item){
      if (model.has(item)){
        idList.push(model.get(item));
      }
    });
    return idList;
  } else {
    throw new TypeError("'id_attribute' not found on the schema: " 
            + JSON.stringify(schema)
            + ', for the model: ' + JSON.stringify(model.attributes));
  }
  
};

/**
 * Create a string ID from the 'id_attribute' of a schema resource definition. -
 * the id_attribute is an array of field specifiers - the 'ID' will be each of
 * these fields, concatenated with a forward slash. Note: the composite key is
 * the public key; although there is often an 'ID' field, it should not be used
 * for URL's, as it is a internal, possibly transient implementation detail.
 */
var getIdFromIdAttribute = Iccbl.getIdFromIdAttribute = 
  function(model, schema){
  
  if (_.has(schema['resource_definition'], 'id_attribute')) {
    var id_attribute = schema['resource_definition']['id_attribute'];
    console.log('create id from ' + id_attribute);
    var id = _.reduce(id_attribute, function(memo, item){
      if(!_.isEmpty(memo)) memo += '/';
      return memo += model.get(item);
    }, '');
    return id;
  } else {
    throw new TypeError("'id_attribute' not found on the schema: " 
        + JSON.stringify(schema)
        + ', for the model: ' + JSON.stringify(model.attributes));
  }
};

/**
 * Create an string 'title' from the 'id_attribute' of a schema resource
 * definition. - the title_attribute is an array of field specifiers and
 * strings. - if an array item is a field, the field value will be used, - if an
 * array item is not a field, then it will be concatenated directly.
 */
var getTitleFromTitleAttribute = Iccbl.getTitleFromTitleAttribute = 
    function(model, schema){
  var re_isQuoted = /['"]+/g;
  var fields = schema['fields'];
  if(_.has(schema['resource_definition'], 'title_attribute')){
    var title = _.reduce(
      schema['resource_definition']['title_attribute'],
      function(memo, item){
        if(item && item.match(re_isQuoted)){
          memo += item.replace(re_isQuoted, '');
        }else{
          if( model.has(item) ){
            var val = model.get(item);
            if (!_.isEmpty(fields[item].vocabulary_scope_ref)){
              val = Iccbl.appModel.getVocabularyTitle(fields[item].vocabulary_scope_ref,val);
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
 * Matches array items against the matchstring. - Matches from the right to
 * left; allowing URI fragments to match their parent URIs. Similar the the
 * contains function, but using item.indexOf(matchString) ||
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
 * Split label strings on non-word characters, then re-join into lines, where
 * each line is less than max_line_length (or, if a part is longer than
 * max_line_length, then that part becomes an entire line).
 */
var createLabel = Iccbl.createLabel = function(original_label, max_line_length, break_char){
  var lines = [];
  var labelParts = original_label.split(/([\s_\-\.,])/);
  var line = '';
  _.each(labelParts, function(part){
    if(line.length > 0){
      if(line.length+part.length <= max_line_length){
        line += part;
      }else{
        lines.push(line);
        line = part;
      }
    }else{
      if( part != ' '){
        line = part;
      }
    }
  });
  lines.push(line);
  
  if(_.isUndefined(break_char)){
    break_char = '<br>';
  }
  return lines.join(break_char);
};

// // Network Utilities ////

var getSchema = Iccbl.getSchema = function(schema_url, callback) {
  $.ajax({
    type : "GET",
    url : schema_url, // options.url_schema,
    data : "",
    dataType : "json",
    success : function(schemaResult) {
        callback(schemaResult);
    }, // end success outer ajax call
    error : function(x, e) {
        alert(x.readyState + " " + x.status + " " + e.msg);
        // TODO: use error div in Bootstrap
    }
  });
};

/**
 * @deprecated see app_state.getModel
 */
var getModel = Iccbl.getModel = function(schemaResult, url, callback) {
  var ModelClass = Backbone.Model.extend({
      url : url,
      defaults : {}
  });
  var instance = new ModelClass();
  instance.fetch({
      success : function(model) {
        model.resourceSchema = schemaResult;
        callback(schemaResult, model);
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

/**
 * Note use this version if the model will be updated to the server. - Backbone
 * will use "put" if the model has an ID set - also, use
 * model.save([attribute-list], {patch: true} ) to do incremental updates, see
 * http://backbonejs.org/#Model-save
 */
var getModel2 = Iccbl.getModel2 = function(schemaResult, urlRoot, id, callback) {
    var ModelClass = Backbone.Model.extend({
        urlRoot : urlRoot,
        defaults : {}
    });
    var instance = new ModelClass({
        id : id
    });
    instance.fetch({
        success : function(model) {
            callback(schemaResult, model);
        },
        error : function(model, response, options) {
            var msg = 'Error locating resource: ' + urlRoot + ', ' + id;
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

var CollectionOnClient = Iccbl.CollectionOnClient = Backbone.Collection.extend({
  /**
   * Override collection parse method: Parse server response data.
   */
  parse : function(response) {
    console.log('Collection on client, parse called');
    // hack the response for tastypie:
    return response.objects;
  },
  setSearch: function(){
    //nop, for now
  }
});

var getCollectionOnClient = Iccbl.getCollectionOnClient = function(url, callback, options){
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

var getCollection = Iccbl.getCollection = function(schemaResult, url, callback) {
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

var getCollection2 = Iccbl.getCollection2 = function(schemaResult, url, callback) {
    var CollectionClass = Iccbl.CollectionInColumns.extend({
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


var getSchemaAndModel = Iccbl.getSchemaAndModel = function(schema_url, url, callback) {
    Iccbl.getSchema(schema_url, function(schemaResult) {
        console.log('schemaResult callback: ' + schemaResult + ', ' + url);
        Iccbl.getModel(schemaResult, url, callback);
    });
};

/**
 * Note use this version if the model will be updated. - Backbone will use "put"
 * if the model has an ID set - also, use model.save([attribute-list], {patch:
 * true} ) to do incremental updates, see http://backbonejs.org/#Model-save
 */
var getSchemaAndModel2 = Iccbl.getSchemaAndModel2 = function(urlRoot, id, callback) {
    Iccbl.getSchema(urlRoot + '/schema', function(schemaResult) {
        console.log('schemaResult callback: ' + schemaResult + ', ' + urlRoot);
        Iccbl.getModel2(schemaResult, urlRoot, id, callback);
    });
};

var getSchemaAndCollection = Iccbl.getSchemaAndCollection = function(schema_url, url, callback) {
    Iccbl.getSchema(schema_url, function(schemaResult) {
        console.log('schemaResult callback: ' + schemaResult + ', ' + url);
        Iccbl.getCollection(schemaResult, url, callback);
    });
};

var getSchemaAndCollection2 = Iccbl.getSchemaAndCollection2 = function(schema_url, url, callback) {
    Iccbl.getSchema(schema_url, function(schemaResult) {
        console.log('schemaResult callback: ' + schemaResult + ', ' + url);
        Iccbl.getCollection2(schemaResult, url, callback);
    });
};


var MyModel = Iccbl.MyModel = Backbone.Model.extend({
    // TODO: we want to make sure there is a trailing slash, or tastypie
    // doesn't work.
    url : function() {
        var url = Backbone.Model.prototype.url.call(this);
        return url + (url.charAt(url.length - 1) === '/' ? '' : '/');
    },

    initialize : function() {
        Backbone.Model.prototype.initialize.apply(this, arguments);
        var self = this;
        // we want to make sure there is a trailing slash, or tastypie doesnt
        // work.
        // TODO: not sure why we have to override url function like this
        // this.url = function(){
        // var url = Backbone.Model.prototype.url.call(self);
        // console.log('---- url1: ' + url);
        // return url + (url.charAt(url.length - 1) === '/' ? '' : '/') ;
        // };
        // definition above should work, but doesn't.
        // however, when overriding url like above only, the function has to
        // be attached to the prototype manually here. why?
        this.url = MyModel.prototype.url;
    },
});


// TODO: redo the link cell like the UriListCell
var LinkCell = Iccbl.LinkCell = Backgrid.Cell.extend({
    className : "link-cell",

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
        Backgrid.Cell.prototype.initialize.apply(this, arguments);
    },

    render : function() {
      var self = this;
      this.$el.empty();
      var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
      var interpolatedVal = Iccbl.replaceTokens(self.model, self.hrefTemplate);
      self.$el.append($('<a>', {
        tabIndex : -1,
        href : interpolatedVal,
        target : self.target
      }).text(formattedValue));
      return this;
    },
});

  
var UriListCell = Iccbl.UriListCell = Backgrid.Cell.extend({
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
      Backgrid.Cell.prototype.initialize.apply(this, arguments);
  },

  render : function() {
    var self = this;
    this.$el.empty();
    var rawValue = this.model.get(this.column.get("name"));
    
    if(rawValue && !_.isEmpty(rawValue)){
      var i = 0;
      _.each(rawValue, function(val){
        var interpolatedVal = Iccbl.replaceTokens(self.model, self.hrefTemplate, val);
        if(i>0) self.$el.append(',');
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

var ImageCell = Iccbl.ImageCell = Backgrid.Cell.extend({
  className : "image-cell",
  events : {
      "click #link" : "toLink",
  },

  initialize : function(options) {
      Backgrid.Cell.prototype.initialize.apply(this, arguments);
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

    // var image_src_attr = this.column.get('backgrid_cell_options')
    //      
    // if (!_.isEmpty(image_src_attr)){
    // // NOTE: format for backgrid cell options is "/{attribute_key}/"
    // var src = Iccbl.replaceTokens(this.model,image_src_attr);
    // console.log('image src: ' + src);
    // return '<img src="'+src+'" width="100" alt="" />';
    // }else{
    // return '';
    // }
  },
      
});

var EditCell = Iccbl.EditCell = Backgrid.Cell.extend({
    className : "detail-cell",
    events : {
        "click #edit" : "editDetail",
    },

    initialize : function(options) {
        this.options = options;
        Backgrid.Cell.prototype.initialize.apply(this, arguments);
    },

    render : function() {
        this.$el.empty();
        var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
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
var NumberCell = Backgrid.NumberCell;

var DecimalFormatter = Iccbl.DecimalFormatter = function () {
  Backgrid.NumberFormatter.apply(this, arguments);
};
 
DecimalFormatter.prototype = new Backgrid.NumberFormatter();
 
_.extend(DecimalFormatter.prototype, {
 
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
              ['n', 1e-9 ],
              ['p', 1e-12 ]
              ],
    multiplier: 1
  }),

  /**
   * Takes a floating point number, where the number is first multiplied by
   * `multiplier`, then converted to a formatted string like
   * NumberFormatter#fromRaw, then finally append `symbol` to the end.
   * 
   * @member Backgrid.UnitsFormatter
   * @param {number}
   *          rawValue
   * @param {Backbone.Model}
   *          model Used for more complicated formatting
   * @return {string}
   */
  fromRaw: function (number, model) {
    return this.getUnit(number, this.multiplier, this.symbol, this.decimals);
  },

  /**
   * Convert the number to a siunit value; number = .0025 multiplier = 1000
   * 
   * @return Math.round(number*multiplier,decimals) + symbol
   */
  getUnit: function(number, multiplier, symbol, decimals) {
     
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
    if(number == 0 ) return number;
    if(!_.isNumber(multiplier)){
      try{
        multiplier = parseFloat(multiplier);
      }catch(e){
        console.log('not a number - multiplier: ' + multiplier+ ', ex:' + e);
      }
    }

    if(multiplier >= 1){
      number = number * multiplier;
    }else{
      console.log("Error, DecimalCell multiplier < 1: " + multiplier);
    }
     
    pair = _.find(this.siunits, function(pair){
      return pair[1] <= Math.abs(number); 
    });
     
    if(_.isUndefined(pair)){
      console.log('could not find units for the input number: ' + number);
      return number;
    }
     
    var val = (1/pair[1])*number;
    val = Math.round(val*Math.pow(10,decimals))/Math.pow(10,decimals);
    return '' + val + ' ' + pair[0] + symbol;
  },
 
  /**
   * FIXME: refactor using SIUnitFormFilter._findNumberAndUnit
   * 
   * @member Backgrid.UnitsFormatter
   * @param {string}
   *          formattedData
   * @param {Backbone.Model}
   *          model Used for more complicated formatting
   * @return {number|undefined} Undefined if the string cannot be converted to a
   *         number.
   */
  toRaw: function (formattedValue, model) {
    var tokens = formattedValue.split(this.symbol);
    if (tokens && tokens[0] && tokens[1] ) {
      var rawValue = NumberFormatter.prototype.toRaw.call(this, tokens[0]);
      if (_.isUndefined(rawValue)) return rawValue;
      var temp = rawValue / this.multiplier;
      var power_mult = _.find(self.siunits, function(pair){
        if(pair[0] == tokens[1]) return true;
      });
      return temp*power_mult;
    }else{
      console.log('error, number of tokens wrong:' + tokens );
    }
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

/**
 * uses the options.attributes.label
 */
var DeleteCell = Iccbl.DeleteCell = Backgrid.Cell.extend({
  className: "delete-cell",
  events : {
  "click #delete" : "delete"
  },

  initialize: function(options){
    Backgrid.Cell.prototype.initialize.apply(this, arguments);
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
    window.alert('ContentView changeUri function must be implemented. uriStack: ' +
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
    
    // Replaces:
    // collection.setSorting(order && column.get("name"), order,
    // {sortValue: column.sortValue()});
    
    console.log('call setSorting: order: ' + order + ', ' + direction );
    collection.setSorting(column.get("name"), order,
        {sortValue: column.sortValue()});
    console.log('post-sort collection fetch');
    collection.fetch({
      reset: true, 
      success: function () {
        console.log('fetch success, direction: ' + direction);
        collection.trigger("backgrid:sorted", column, direction, collection);
      }
//      error: Iccbl.appModel.backboneFetchError
    }).fail(function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments); });      
    
    column.set("direction", direction);

    return this;
  }
});


var MyCollection = Iccbl.MyCollection = Backbone.PageableCollection.extend({

  initialize : function(options) {
    var self = this;
    this.options = options;
    this.url = options.url;
    this.listModel = options.listModel;

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
  queryParams : {
    // adjust the query params for tastypie
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
   * Method for external callers to set the search
   */
  setSearch: function(searchHash) {
    var self = this;
    var searchHash = _.clone(searchHash);
    self.listModel.set('search', searchHash);
    
    // Tell all the header cells
    // TODO: test further - REMOVED 20150113
    this.trigger("MyServerSideFilter:search", searchHash, this);

    // Allow searches that aren't for a visible column:
    // if the search key is not in the queryParams, then it is not a column
    // search (TODO: verify).
    // this will add it manually to the queryParams (which are serialized in
    // the fetch to the server)
    var _data = {};
    _.each(_.keys(searchHash), function(key) {
      var val = searchHash[key]

      if(_.isEmpty(val)){
          delete self.queryParams[key];
      }else{
        // check if param dne, or if param exists and is a value to be set
        // the reason for the "isFunction" check is that the Backgrid-filter
        // defined params are function calls to get the current value in the
        // searchbox - so skip those as state is stored there.
        if (!_.has(self.queryParams, key) || !_.isFunction(self.queryParams[key])) {
        	_data[key]=val;
        	
        	// make the params persistent (if not a backgrid-filter)
          self.queryParams[key] = function () {
            return self.listModel.get('search')[key] || null;
          };
        }
      }
    });

    if(!_.isEmpty(_data)){
      self.fetch().fail(
        function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments);}
      );      
      // TODO: 2014-04-11: not sure why removing this works:
      // removed to fix sort not working when custom searches
      // are used. Custom searches seem to still work: further testing of search
      // may reveal problems. (problem was that sort badge indicators would not
      // render
      // when there was a custom search).
      // (this was cargo culted from backgrid header search)
      // 
      // self.fetch({data:_.clone(_data), reset: true});
    }else{
      // TODO: test further, part of the double network hit on search - REMOVED
      // 20150113
      self.fetch().fail(
        function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments);}
      );      
    }
  },

  /**
   * Proxy for the search elements to add search terms to the listModel
   */ 
  addSearch: function(searchHash) {
    var self = this;
    var oldsearchHash = _.clone(self.listModel.get('search'));
    console.log('collection addSearch: current: ' + 
    		JSON.stringify(oldsearchHash) + 
    		', adding: ' + JSON.stringify(searchHash));
    oldsearchHash = _.extend(oldsearchHash, searchHash);
    self.listModel.set({
        'search' : oldsearchHash
    });
    _.each(_.keys(searchHash), function(key){
      // make the params persistent (if not a backgrid-filter)
      self.queryParams[key] = function () {
        return self.listModel.get('search')[key] || null;
      };
    });
  },

  /**
   * Proxy for the search elements to clear search terms from the listModel on
   * the collection.
   */ 
  clearSearch: function(searchKeys, options) {
    console.log('clearsearch: ' + JSON.stringify(searchKeys) 
        + ', options: ' + JSON.stringify(options) );
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
      self.listModel.set({
        'search' : searchHash
      });
      if(options && options.reset){
        console.log('collection.clearSearch: reset');
        self.getFirstPage({reset: true, fetch: true});
      }
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



var MultiSortHeaderCell = Iccbl.MultiSortHeaderCell = Backgrid.HeaderCell.extend({
  
  filtericon_text : '<span class="pull-left glyphicon glyphicon-search" id="filter-icon" ></span>',
 
  initialize : function(options) {
    this.options = options;
    MultiSortHeaderCell.__super__.initialize.apply(this, arguments);
    this.listenTo(this.collection,"sort",this.collectionSorted);
    this.listenTo(this.collection,"Iccbl:clearSorts", this.removeCellDirection);
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
    
    this.setCellDirection(column, self.tempdirection=="none"?null:self.tempdirection );

    var args = arguments;
    
    var delayedClick = function(){
      console.log('delayedclick: tempdirection: ' + self.tempdirection);
      if(self.tempdirection !== self.lastExecutedVal){
        console.log('delayed click: ' + self.tempdirection );
        collection.trigger(event, column, self.tempdirection=="none"?null:self.tempdirection);
        self.lastExecutedVal = self.tempdirection;
      }else{
        console.log('this.tempdirection == self.lastExecutedVal: ' + self.lastExecutedVal);
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
      return this;
    }
   
  });


var BackgridFormFilter = Backbone.Form.extend({
  
  template: _.template([
    "<form class='form-horizontal container' >",
    '<div class="row center-block" style="margin: 0 0 0 0;" >',
    "<div data-fieldsets   />",
    "</div>",
    "</form>"].join('')),
  checkboxTemplate: [
    '<label data-editor title="<%= help %>" class="checkbox-inline" >',
    '<%= title %></label>'      
    ].join(''),
   
  initialize : function(options) {
    var self = this;
    var options = this.options = options || {};
      
    if(!options.columnName){
      throw "must define column name option for the BackgridFormFilter, options: " + options;
    }else{
      this.columnName = options.columnName;
    }
    BackgridFormFilter.__super__.initialize.apply(this, arguments);
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
  }

});


var CriteriumFormFilter = Iccbl.CriteriumFormFilter = BackgridFormFilter.extend({
  criterium: {'=':'eq'},
  errorClass: 'has-error',
  criteriaTemplate: 
    [
      '<span  data-editor></span>'
    ].join(''),
  fieldTemplate: [
      '<div data-editor title="<%= help %>" class="form-control" >',
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
   * Convenience - determine if the form has been set with any values
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
  }
  
});

var TextFormFilter = CriteriumFormFilter.extend({
  
  criterium: {'=':'eq','contains':'contains','icontains':'icontains','<>':'ne', 
    'in': 'in','blank':'is_null','not blank':'not_blank'},
    
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
        template: _.template(self.checkboxTemplate),
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
  
  /**
   * TextFormFilter handles router generated searches
   */
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
      if(searchVal){
        found = true;
        self.setValue('lower_criteria', criteriaKey);
        if(criteria == 'is_null'){
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

  /**
   * TextFormFilter Form submit handler
   */
  _submit: function(){
    var self  = this;

    if(!self.isSet()) return;
    
    var searchHash = {};
    
    // validate:true: tells bbf to run model.validate(), in addition to
    // field[].validate()
    var errors = self.commit({ validate: true }); 
    if(!_.isEmpty(errors)){
      console.log('form errors, abort submit: ' + JSON.stringify(errors));
      return;
    }else{
      // this.$el.find('#range_upper_block').removeClass(self.errorClass);
    }
    
    var values = self.getValue();
    var criteria = self.criterium[values['lower_criteria']];
    var searchKey = self.columnName + '__' + criteria;
    var searchVal = values['form_textarea'];
    if(criteria == 'not_blank'){
      searchKey = self.columnName + '__' + 'is_null';
      searchVal = 'false'
    }else if(criteria == 'is_null'){
      searchVal = 'true';
    }
    var invert = values['invert_field'];
    if(invert) searchKey = '-'+searchKey;
    searchHash[searchKey] = searchVal;
  
    return searchHash
  }
  
});

var DateEditor = Backbone.Form.editors.Date.extend({

  /**
   * need to extend because an error in the bbf initializer makes it impossible
   * to override the monthnames otherwise. (Using var Self = Form.editors.Date;
   * this.options = _.extend({ monthNames: Self.monthNames, why?
   */
  monthNames: ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],

  initialize : function(options) {
    var self = this;
    DateEditor.__super__.initialize.apply(this, arguments);
    // now override monthnames
    self.options.monthNames = self.monthNames;
  },
});

var DateFormFilter = CriteriumFormFilter.extend({
  
  criterium: {'=':'eq','>':'gt','>=':'gte','<':'lt','<=':'lte','<>':'ne',
    'between':'range', 'in': 'in','blank':'is_null','not blank':'not_blank'},

  // provide a custom form template; use Bootstrap layout/styling
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
        template: _.template(self.checkboxTemplate),
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
    this.selectedFields = ['lower_criteria','lower_value','form_textarea','upper_value','invert_field']; 
    
    this.listenTo(this, "change", function(e){
      var criteria = self.getValue('lower_criteria');
      console.log('change:' + criteria)
      if(criteria == 'between'){
        self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').show();
      }else if(criteria == 'in'){
        self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', true);
        self.setValue('lower_value', '');
        self.$el.find('#range_upper_block').hide();
        self.$el.find('[data-fields="form_textarea"]').show();
      }else{
        self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').hide();
      }
    });

    DateFormFilter.__super__.initialize.apply(this, arguments);
    
  },

  /**
   * DateFormFilter
   */
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
      if(searchVal){
        found = true;
        self.setValue('lower_criteria', criteriaKey);
        try{
          if(criteria == 'range'){
            self.$el.find('#range_upper_block').show();
            var vals = searchVal.split(',');
            if(vals.length < 2){
              throw "the range filter requires 2 date arguments separated by a comma, given: " + searchVal
            }
            self.setValue('lower_value', new Date(vals[0]));
            self.setValue('upper_value', new Date(vals[1]));
          }else if(criteria == 'in'){
            self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', true);
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

  /**
   * DateFormFilter Form submit handler
   */
  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;
    var searchHash = {};
    
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
  }
});


var TextHeaderCell = MultiSortHeaderCell.extend({

  initialize : function(options) {
    var self = this;
    
    this.options = options;
    TextHeaderCell.__super__.initialize.apply(this, arguments);
    this.fieldinformation = _.clone(this.column.get('fieldinformation'));
  
    this._serverSideFilter = new TextFormFilter({
      columnName: this.column.get('name')
    });
    this.filterIcon = $(self.filtericon_text);
    this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
    this.listenTo(this.collection,"Iccbl:clearSearches",this.clearSearch);
    _.bindAll(this, '_submit', 'clearSearch');
  },

  
  /**
   * TextHeaderCell Form submit handler
   */
  _submit: function(e){
    var self  = this;
    console.log('_submit called');
    if (e) e.preventDefault();      
  
    var searchHash = self._serverSideFilter._submit();
    if(!_.isEmpty(searchHash)){
      var possibleSearches = self._serverSideFilter.getPossibleSearches();
      self.collection.clearSearch(possibleSearches);
      
      console.log('server side filter add search: ' + 
          JSON.stringify(searchHash));
      this.collection.addSearch(searchHash);
    }else{
      console.log('nothing submitted');
    }
  },

  /**
   * TextHeaderCell Listen for router generated search events
   */
  _search: function(hash, collection){
    var self = this;
    var name = this.column.get('name');
    var searchHash = _.clone(hash);

    // TODO: could use form.isSet() instead of found
    var found = this._serverSideFilter._search(searchHash);
    
    if(found){
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    }else{
      self._serverSideFilter.$el.hide();
      self.filterIcon.show();
    }
  },  
  
  /**
   * TextHeaderCell clears all possible searches if options.reset, then trigger
   * a collection fetch
   */
  clearSearch: function(options){
    var self=this;
    self._serverSideFilter.$el.hide();
    self.filterIcon.show();
    
    var possibleSearches = self._serverSideFilter.getPossibleSearches();
    self.collection.clearSearch(possibleSearches);
  },      
  
  /**
   * TextHeaderCell
   */
  render : function() {
    var self = this;
    TextHeaderCell.__super__.render.apply(this);
  
    this.$el.append(this.filterIcon);
  
    // Backbone.View.prototype.manage = false;
    this._serverSideFilter.render();
    this.$el.append(this._serverSideFilter.el);
    this._serverSideFilter.$el.hide();
    // Backbone.View.prototype.manage = false;
    
    this._serverSideFilter.clearButton().click(function(e){
      e.preventDefault();
      e.stopPropagation();
      self.clearSearch();
    });
    
    this.filterIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    });
  
    this._serverSideFilter.submitButton().click(function(e){
      e.preventDefault();
      self._submit();
    });
  
    this.$el.prop('title', this.options['column']['attributes']["description"]);
    
    return this;
  },  
});


var DateHeaderCell = MultiSortHeaderCell.extend({

  initialize : function(options) {
    var self = this;
    
    this.options = options;
    DateHeaderCell.__super__.initialize.apply(this, arguments);
    this.fieldinformation = _.clone(this.column.get('fieldinformation'));
    
    this._serverSideFilter = new DateFormFilter({
      columnName: this.column.get('name')
    });
    this.filterIcon = $(self.filtericon_text);
    
    this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
    this.listenTo(this.collection,"Iccbl:clearSearches",this.clearSearch);
    _.bindAll(this, '_submit', 'clearSearch');
  },
  
  
  /**
   * DateHeaderCell Form submit handler
   */
  _submit: function(e){
    var self  = this;
    console.log('_submit called');
    if (e) e.preventDefault();      
    
    var name = this.column.get('name');
    var searchHash = self._serverSideFilter._submit();
    if(!_.isEmpty(searchHash)){
      var possibleSearches = self._serverSideFilter.getPossibleSearches();
      self.collection.clearSearch(possibleSearches);
      console.log('server side filter add search: ' + 
          JSON.stringify(searchHash));
      this.collection.addSearch(searchHash);
    }else{
      console.log('Date Cell nothing submitted');
    }
  },
  
  /**
   * DateHeaderCell Listen for router generated search events
   */
  _search: function(hash, collection){
    var self = this;
    var name = this.column.get('name');
    var searchHash = _.clone(hash);

    // TODO: could use form.isSet() instead of found
    var found = this._serverSideFilter._search(searchHash);
    if(found){
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    }else{
      self._serverSideFilter.$el.hide();
      self.filterIcon.show();
    }
  },  

  /**
   * DateHeaderCell clears all possible searches if options.reset, then trigger
   * a collection fetch
   */
  clearSearch: function(options){
    var self=this;
    self._serverSideFilter.$el.hide();
    self.filterIcon.show();
    
    var possibleSearches = self._serverSideFilter.getPossibleSearches();
    self.collection.clearSearch(possibleSearches);
  },      

  /**
   * DateHeaderCell
   */
  render : function() {
    var self = this;
    DateHeaderCell.__super__.render.apply(this);
  
    this.$el.append(this.filterIcon);
    this._serverSideFilter.render();
    this.$el.append(this._serverSideFilter.el);
    this._serverSideFilter.$el.hide();

    this._serverSideFilter.clearButton().click(function(e){
      e.preventDefault();
      e.stopPropagation();
      self.clearSearch();
    });
    this.filterIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    });
    this._serverSideFilter.submitButton().click(function(e){
      e.preventDefault();
      self._submit();
    });
  
    this.$el.prop('title', this.options['column']['attributes']["description"]);
    
    return this;
  },  
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

  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;
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
    
    return searchHash
  },

  /**
   * BooleanFormFilter Set the form from router generated event
   */
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
  
  /**
   * BooleanFormFilter Convenience - determine if the form has been set with any
   * values
   */
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

var BooleanHeaderCell = MultiSortHeaderCell.extend({

  initialize : function(options) {
    var self = this;
    
    this.options = options;
    BooleanHeaderCell.__super__.initialize.apply(this, arguments);
    this.fieldinformation = _.clone(this.column.get('fieldinformation'));
  
    this._serverSideFilter = new BooleanFormFilter({
      columnName: this.column.get('name')
    });
    this.filterIcon = $(self.filtericon_text);
    this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
    this.listenTo(this.collection,"Iccbl:clearSearches",this.clearSearch);
    _.bindAll(this, '_submit', 'clearSearch');
  },

  
  /**
   * BooleanHeaderCell Form submit handler
   */
  _submit: function(e){
    var self  = this;
    console.log('_submit called');
    if (e) e.preventDefault();      
  
    var searchHash = self._serverSideFilter._submit();
    if(!_.isEmpty(searchHash)){
      var possibleSearches = self._serverSideFilter.getPossibleSearches();
      self.collection.clearSearch(possibleSearches);
      console.log('server side filter add search: ' + 
          JSON.stringify(searchHash));
      this.collection.addSearch(searchHash);
    }else{
      console.log('nothing submitted');
    }
  },

  /**
   * BooleanHeaderCell Listen for router generated search events
   */
  _search: function(hash, collection){
    var self = this;
    var name = this.column.get('name');
    var searchHash = _.clone(hash);

    // TODO: could use form.isSet() instead of found
    var found = this._serverSideFilter._search(searchHash);
    
    if(found){
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    }else{
      self._serverSideFilter.$el.hide();
      self.filterIcon.show();
    }
  },  
  
  /**
   * BooleanHeaderCell clears all possible searches if options.reset, then
   * trigger a collection fetch
   */
  clearSearch: function(options){
    var self=this;
    var possibleSearches = self._serverSideFilter.getPossibleSearches();
    self._serverSideFilter.$el.hide();
    self.filterIcon.show();
    self.collection.clearSearch(possibleSearches);
  },      
  
  /**
   * BooleanHeaderCell
   */
  render : function() {
    var self = this;
    BooleanHeaderCell.__super__.render.apply(this);
    this.$el.append(this.filterIcon);
    this._serverSideFilter.render();
    this.$el.append(this._serverSideFilter.el);
    this._serverSideFilter.$el.hide();
    this._serverSideFilter.clearButton().click(function(e){
      e.preventDefault();
      e.stopPropagation();
      self.clearSearch();
    });
    
    this.filterIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    });
  
    this._serverSideFilter.submitButton().click(function(e){
      e.preventDefault();
      self._submit();
    });
  
    this.$el.prop('title', this.options['column']['attributes']["description"]);
    
    return this;
  },  
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
      '    <div class="checkbox" style="min-height: 0px; padding-top: 0px;" > ',
      '      <label for="<%= editorId %>"><span data-editor\><%= title %></label>',
      '    </div>',
      '  </div>'
      ].join('')),
    
  initialize : function(options) {
    var self = this;
    var options = this.options = options || {};
    var formSchema = this.schema = {};
    
    if(!options.choiceHash){
      throw "must define choiceHash option for the SelectorFormFilter, options: " + options;
    }else{
      var choiceHash = this.choiceHash = options.choiceHash;
    }
    if(!options.selectedFields){
      throw "must define selectedFields option for the SelectorFormFilter, options: " + options;
    }else{
      var selectedFields = this.selectedFields = options.selectedFields || [];
    }
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
        template: _.template(self.checkboxTemplate),
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

  /**
   * SelectorFormFilter Form submit handler
   */
  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;
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
    
    return searchHash
  },

  /**
   * SelectorFormFilter Set the form from router generated event
   */
  _search: function(hash){
    var self = this;
    var searchHash = _.clone(hash);
    var searchTerm = null;
    var name = this.columnName;

    _.each(self.getPossibleSearches(), function(term){
      if(_.has(searchHash,term)){
        searchTerm = term;
      }
    });
    
    if(searchTerm){
      if(searchTerm.charAt(0) == '-'){
        self.setValue('invert_field', true);
        name = '-' + name;
      }
      if(searchTerm ==  name + '__is_null'){
        if(searchVal == 'true'){
          self.setValue('lower_criteria','blank');
        }else{
          self.setValue('lower_criteria', 'not blank');
        }
      }else{
        var searchVal = searchHash[searchTerm];
        _.each(searchVal.split(','), function(choice){
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
      if(key == 'invert_field' ) return false;
      return values[key];
    });
    return found;
  },
  
  getPossibleSearches: function(){
    return [this.columnName + '__in', '-'+this.columnName + '__in',
            this.columnName + '__is_null'];
  },
});

var SelectorHeaderCell = MultiSortHeaderCell.extend({

  initialize : function(options) {
    SelectorHeaderCell.__super__.initialize.apply(this, arguments);

    var self = this;
    var choiceHash = {}
    var vocabulary;
    
    this.fieldinformation = _.clone(this.column.get('fieldinformation'));
    this.options = options;

    if(_.isUndefined(this.fieldinformation.choices)){
      console.log([
         'error: fieldinformation for a selection field type must define a ',
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
        vocabulary = Iccbl.appModel.getVocabulary(this.fieldinformation.vocabulary_scope_ref);
          _.each(_.keys(vocabulary),function(choice){
            choiceHash[choice] = vocabulary[choice].title;
          });
      }catch(e){
        console.log('vocabulary error', this.column.get('name'),e);
      }
    }
    this._serverSideFilter = new SelectorFormFilter({
        choiceHash: choiceHash,
        columnName: this.column.get('name'),
        selectedFields: _.clone(this.fieldinformation.choices)
      });
    
    this.filterIcon = $(self.filtericon_text);
    
    this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
    this.listenTo(this.collection,"Iccbl:clearSearches",this.clearSearch);
    _.bindAll(this, '_submit', 'clearSearch');
  },
  
  /**
   * SelectorHeaderCell Form submit handler
   */
  _submit: function(e){
    console.log('_submit called');
    if (e) e.preventDefault();      

    var self  = this;
    var searchHash = self._serverSideFilter._submit();
    var possibleSearches = self._serverSideFilter.getPossibleSearches();
    self.collection.clearSearch(possibleSearches);
    if(!_.isEmpty(searchHash)){
      console.log('server side filter add search: ' + 
          JSON.stringify(searchHash));
      this.collection.addSearch(searchHash);
    }else{
      self.collection.fetch().fail(
        function(){ Iccbl.appModel.jqXHRfail.apply(this,arguments);}
      );      
      
    }
  },
  
  /**
   * SelectorHeaderCell Listen for router generated search events
   */
  _search: function(hash, collection){
    var self = this;
    var searchHash = _.clone(hash);
    // TODO: could use form.isSet() instead of found
    var found = this._serverSideFilter._search(searchHash);
    if(found){
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    }else{
      self._serverSideFilter.$el.hide();
      self.filterIcon.show();
    }
  },  

  /**
   * SelectorHeaderCell clears all possible searches if options.reset, then
   * trigger a collection fetch
   */
  clearSearch: function(options){
    var self=this;
    self._serverSideFilter.$el.hide();
    self.filterIcon.show();
    
    var possibleSearches = self._serverSideFilter.getPossibleSearches();
    self.collection.clearSearch(possibleSearches);
  },      

  /**
   * SelectorHeaderCell
   */
  render : function() {
    var self = this;
    SelectorHeaderCell.__super__.render.apply(this);
  
    this.$el.append(this.filterIcon);
    this._serverSideFilter.render();
    this.$el.append(this._serverSideFilter.el);
    this._serverSideFilter.$el.hide();
    this._serverSideFilter.clearButton().click(function(e){
      e.preventDefault();
      e.stopPropagation();
      self.clearSearch();
    });
    
    this.filterIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    });

    this._serverSideFilter.submitButton().click(function(e){
      e.preventDefault();
      self._submit();
    });

    this.$el.prop('title', 
            this.options['column']['attributes']["description"]);
    
    return this;
  }
  
});


var NumberFormFilter = CriteriumFormFilter.extend({
  
  criterium: {'=':'eq','\u2248':'about','>':'gt', '>=':'gte','<':'lt','<=':'lte',
    '<>':'ne', 'x..y':'range', 'in': 'in','blank':'is_null','not blank':'not_blank'},

  // use Bootstrap layout/styling
  template: _.template([
      '<form class="iccbl-headerfield-form" >',
      '<div class="row center-block" style="margin: 0 0 0 0;" >',
      '<div class="input-group ">',
      '   <div data-fields="lower_criteria" ',
      '     class="input-group-addon iccbl-headerfield-number" for="lower_value"   />',
      '   <div data-fields="lower_value"/>',
      '</div>',
      '<div class="form-group" data-fields="form_textarea" style="display: none;" />',
      '<div class="input-group" id="range_upper_block" style="display: none;" >',
      '   <span class="input-group-addon" for="upper_value"  style="width: 4em; ">to</span>',
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
        template: _.template(self.checkboxTemplate),
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
    
    options.fields = fields.concat(
        ['lower_criteria','lower_value','form_textarea','upper_value','invert_field']); 
    
    this.listenTo(this, "change", function(e){
      var criteria = self.getValue('lower_criteria');
      if(criteria == 'x..y'){
        self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').show();
      }else if(criteria == 'in'){
        self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', true);
        self.setValue('lower_value', '');
        self.$el.find('#range_upper_block').hide();
        self.$el.find('[data-fields="form_textarea"]').show();
      }else{
        self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', false);
        self.$el.find('[data-fields="form_textarea"]').hide();
        self.$el.find('#range_upper_block').hide();
      }
    });

    NumberFormFilter.__super__.initialize.call(this, options);
    
  },

  /**
   * NumberFormFilter Convenience - determine if the form has been set with any
   * values
   */
  isSet: function(){
    var values = this.getValue();
    var found = _.find(_.keys(values), function(key){
      return values[key] !== '';
    });
    return !_.isEmpty(found);
  },
  
  /**
   * NumberFormFilter Set the form from router generated event
   */
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
      if(searchVal){
        found = true;
        self.setValue('lower_criteria', criteriaKey);
        if(criteria == 'range'){
          self.$el.find('#range_upper_block').show();
          var vals = searchVal.split(',');
          self.setValue('lower_value', vals[0]);
          self.setValue('upper_value', vals[1]);
        }else if(criteria == 'in'){
          self.$el.find('[data-fields="lower_value"]').find('input').prop('disabled', true);
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

  /**
   * NumberFormFilter Form submit handler
   */
  _submit: function(){
    var self  = this;
    if(!self.isSet()) return;
    var searchHash = {};
    var errors,values,name,invert,criteria,searchKey;

    // validate:true: tells bbf to run model.validate(), in addition to
    // field[].validate()
    errors = self.commit({ validate: true }); 
    if(!_.isEmpty(errors)){
      console.log('form errors, abort submit: ' + JSON.stringify(errors));
      this.$el.find('#range_upper_block').addClass(self.errorClass);
      return;
    }else{
      this.$el.find('#range_upper_block').removeClass(self.errorClass);
    }
    values = self.getValue();
    name = self.columnName;
    invert = values['invert_field'];
    if(invert) name = '-'+name;
    criteria = self.criterium[values['lower_criteria']];
    searchKey = name + '__' + criteria;
    
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
    return searchHash
  }

});

var NumberHeaderCell = MultiSortHeaderCell.extend({

  initialize : function(options) {
    var self = this;
    NumberHeaderCell.__super__.initialize.apply(this, arguments);
    this.options = options;

    this.fieldinformation = _.clone(this.column.get('fieldinformation'));
    
    this._serverSideFilter = new NumberFormFilter({
      columnName: this.column.get('name')
    });
    this.filterIcon = $(self.filtericon_text);
    
    this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
    this.listenTo(this.collection,"Iccbl:clearSearches",this.clearSearch);
    _.bindAll(this, '_submit', 'clearSearch');
  },
  
  
  /**
   * NumberHeaderCell Form submit handler
   */
  _submit: function(e){
    console.log('_submit called');
    if (e) e.preventDefault();      
    var self  = this;
    var name = this.column.get('name');
    var searchHash = self._serverSideFilter._submit();
    
    if(!_.isEmpty(searchHash)){
      self.collection.clearSearch(self._serverSideFilter.getPossibleSearches());
      console.log('server side filter add search: ' + 
          JSON.stringify(searchHash));
      this.collection.addSearch(searchHash);
    }else{
      console.log('nothing submitted');
    }
  },
  
  /**
   * NumberHeaderCell Listen for router generated search events
   */
  _search: function(hash, collection){
    var self = this;
    var name = this.column.get('name');
    var searchHash = _.clone(hash);
    // TODO: could use form.isSet() instead of found
    var found = this._serverSideFilter._search(searchHash);

    if(found){
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    }else{
      self._serverSideFilter.$el.hide();
      self.filterIcon.show();
    }
  },  

  /**
   * NumberHeaderCell clears all possible searches if options.reset, then
   * trigger a collection fetch
   */
  clearSearch: function(options){
    var self=this;
    self._serverSideFilter.$el.hide();
    self.filterIcon.show();
    self.collection.clearSearch(self._serverSideFilter.getPossibleSearches());
  },      
  
  /**
   * NumberHeaderCell
   */
  render : function() {
    var self = this;
    NumberHeaderCell.__super__.render.apply(this);
  
    this.$el.append(this.filterIcon);
    this._serverSideFilter.render();
    this.$el.append(this._serverSideFilter.el);
    this._serverSideFilter.$el.hide();
    this._serverSideFilter.clearButton().click(function(e){
      e.preventDefault();
      e.stopPropagation();
      self.clearSearch();
    });
    
    this.filterIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    });
  
    this._serverSideFilter.submitButton().click(function(e){
      e.preventDefault();
      self._submit();
    });
  
    this.$el.prop('title', this.options['column']['attributes']["description"]);
    
    return this;
  },  
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
    var options = this.options = options || {};
    var multiplier = this.multiplier = options.multiplier || 1;
    var symbol = this.symbol = options.symbol;
    var units = this.units = [];
    var formSchema = options.schema = options.schema || {};
    
    if(! options.symbol){
      throw 'Error: SIUnitFormFilter requires a "symbol" option' 
    }
    _.each(this.siunits,function(pair){
      if(options.maxunit){
        if(options.maxunit < pair[1]) return;
      }
      if(options.minunit){
        if(options.minunit > pair[1]) return;
      }
      units.push({ val: pair[1], label: pair[0] + self.symbol });
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
  
  /**
   * SIUnitFormFilter Convenience - determine if the form has been set with any
   * values
   */
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
  
  /**
   * SIUnitFormFilter Form submit handler
   */
  _submit: function(){
    var self = this;
    var searchHash = SIUnitFormFilter.__super__._submit.call(this);
    
    if(!_.isEmpty(searchHash)){
      var searchKey = _.keys(searchHash)[0];
      var searchValue = searchHash[searchKey];
      var values = self.getValue();
      var criteria = self.criterium[values['lower_criteria']];

      if(criteria == 'range'){
        searchHash[searchKey] = [
            self._calculate(self.multiplier,values['lower_siunit'],values['lower_value']),
            self._calculate(self.multiplier,values['upper_siunit'],values['upper_value'])
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
    return searchHash
  },
  
  _calculate: function(multiplier, sci_mult, val){
    // run strip after every calculation to round out floating point math errors
    function strip(number) {
      return (parseFloat(number.toPrecision(12)));
      };
    val = strip(val * multiplier);
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
  
  /**
   * SIUnitFormFilter Set the form from router generated event
   */
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


var SIUnitHeaderCell = MultiSortHeaderCell.extend({

  /** @property */
  symbol: "",
  
  initialize: function(options){

    SIUnitHeaderCell.__super__.initialize.apply(this,arguments);

    var self = this;
    var name = this.column.get('name');
    // FIXME: cannot pass initialization params, so sending them on the column
    var fi = this.fieldinformation = _.clone(this.column.get('fieldinformation'));
    this.options = options;
    var cell_options = fi['display_options'];
    if(!cell_options || !cell_options.symbol)
    {
      throw new Error('SIUnitHeaderCell: field information requires the '+
          '"symbol" option');
    }
    
    this._serverSideFilter = new SIUnitFormFilter(_.extend({
      columnName: this.column.get('name')}, cell_options ));
    this.filterIcon = $(self.filtericon_text);
    this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
    this.listenTo(this.collection,"Iccbl:clearSearches",this.clearSearch);
  },

  
  /**
   * SIUnitHeaderCell Form submit handler
   */
  _submit: function(e){
    var self  = this;
    console.log('_submit called');
    if (e) e.preventDefault();      
    
    var name = this.column.get('name');
    var searchHash = self._serverSideFilter._submit();
    if(!_.isEmpty(searchHash)){
      var possibleSearches = self._serverSideFilter.getPossibleSearches();
      self.collection.clearSearch(possibleSearches);
      console.log('server side filter add search: ' + 
          JSON.stringify(searchHash));
      this.collection.addSearch(searchHash);
    }else{
      console.log('nothing submitted');
    }
  },
  
  /**
   * SIUnitHeaderCell Listen for router generated search events
   */
  _search: function(hash, collection){
    var self = this;
    var name = this.column.get('name');
    var searchHash = _.clone(hash);

    // TODO: could use form.isSet() instead of found
    var found = this._serverSideFilter._search(searchHash);
    if(found){
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    }else{
      self._serverSideFilter.$el.hide();
      self.filterIcon.show();
    }
  },  

  /**
   * SIUnitHeaderCell clears all possible searches if options.reset, then
   * trigger a collection fetch
   */
  clearSearch: function(options){
    var self=this;
    self._serverSideFilter.$el.hide();
    self.filterIcon.show();
    
    var possibleSearches = self._serverSideFilter.getPossibleSearches();
    self.collection.clearSearch(possibleSearches);
  }, 
  
  /** SIUnitHeaderCell * */
  render : function() {
    SelectorHeaderCell.__super__.render.apply(this);
    var self = this;
  
    this.$el.append(this.filterIcon);

    this._serverSideFilter.render();
    this.$el.append(this._serverSideFilter.el);
    this._serverSideFilter.$el.hide();
    
    this._serverSideFilter.clearButton().click(function(e){
      e.preventDefault();
      e.stopPropagation();
      self.clearSearch();
    });
    
    this.filterIcon.click(function(e){
      e.stopPropagation();
      e.preventDefault();
      self._serverSideFilter.$el.show();
      self.filterIcon.hide();
    });

    this._serverSideFilter.submitButton().click(function(e){
      e.preventDefault();
      self._submit();
    });

    this.$el.prop('title', 
            this.options['column']['attributes']["description"]);
    
    // TODO: customize the default units values
    this._serverSideFilter.setValue('lower_siunit',1);
    this._serverSideFilter.setValue('upper_siunit',1);
    return this;
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
    var rawData = this.formatter.fromRaw(model.get(this.column.get("name")), model);

    var selectedText = [];

    if (_.isArray(optionValues) &&  !_.isEmpty(optionValues)){
      for (var k = 0; k < rawData.length; k++) {
        var rawDatum = rawData[k];

        for (var i = 0; i < optionValues.length; i++) {
          var optionValue = optionValues[i];

          if (_.isArray(optionValue)) {
            var optionText  = optionValue[0];
            var optionValue = optionValue[1];

            if (optionValue == rawDatum){
              selectedText.push(optionText);
            }
          }
          else if (_.isObject(optionValue)) {
            var optionGroupValues = optionValue.values;

            for (var j = 0; j < optionGroupValues.length; j++) {
              var optionGroupValue = optionGroupValues[j];
              if (optionGroupValue[1] == rawDatum) {
                selectedText.push(optionGroupValue[0]);
              }
            }
          }
          else {
            throw new TypeError;
          }
        }
      }
      
    }
    if(!_.isEmpty(rawData) && _.isEmpty(selectedText)){
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
      var interpolatedVal = Iccbl.replaceTokens(this.model, this.hrefTemplate);
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
  ICCBL_DATE_RE:  /^(\d{1,2})\/(\d{1,2})\/(\d{2,4})$/,
  // ISO DATE/TIME regex, from Backgrid
  DATE_RE: /^([+\-]?\d{4})-(\d{2})-(\d{2})$/,
  TIME_RE: /^(\d{2}):(\d{2}):(\d{2})(\.(\d{3}))?$/,
  ISO_SPLITTER_RE: /T|Z| +/,

  fromRaw: function (rawData, model) {
    if (_.isNull(rawData) || _.isUndefined(rawData)) return '';
    rawData = rawData.trim();
    if ((rawData + '').trim() === '') return null;

    if(this.ICCBL_DATE_RE.test(rawData)){
      var DDMMYYYY = this.ICCBL_DATE_RE.exec(rawData) || [];
      var jsDate = new Date(Date.UTC(DDMMYYYY[3] * 1 || 0,
        DDMMYYYY[2] * 1 - 1 || 0,
        DDMMYYYY[1] * 1 || 0));
      return getIccblUTCDateString(jsDate);
    }else{
      var temp = rawData.split('T')[0];
      if(this.DATE_RE.test(temp)){
        var YYYYMMDD = this.DATE_RE.exec(temp);
        var jsDate = new Date(Date.UTC(YYYYMMDD[1]*1, YYYYMMDD[2]*1-1, YYYYMMDD[3]*1 ))
        return getIccblUTCDateString(jsDate);
      }else{
        Iccbl.appModel.error('unrecognized date: ' + rawData );
      }
    }
  },

  toRaw: function(formattedData, model){
    if (_.isNull(formattedData) || _.isUndefined(formattedData)) return '';
    if ((formattedData + '').trim() === '') return null;
    if(this.ICCBL_DATE_RE.test(formattedData)){
      var DDMMYYYY = this.ICCBL_DATE_RE.exec(formattedData) || [];
      var jsDate = new Date(Date.UTC(DDMMYYYY[3] * 1 || 0,
        DDMMYYYY[2] * 1 - 1 || 0,
        DDMMYYYY[1] * 1 || 0));
      var temp = getISODateString(jsDate);
      return temp;
    }else{
      return;
    }
  },
  
  
  /**
   * Use Backgrid DatetimeFormatter convert, add in ICCBL_DATE_RE
   */
  _convert: function (data, validate) {
    if ((data + '').trim() === '') return null;

    var date, time = null;
    if (_.isNumber(data)) {
      var jsDate = new Date(data);
      date = lpad(jsDate.getUTCFullYear(), 4, 0) + '-' + lpad(jsDate.getUTCMonth() + 1, 2, 0) + '-' + lpad(jsDate.getUTCDate(), 2, 0);
      time = lpad(jsDate.getUTCHours(), 2, 0) + ':' + lpad(jsDate.getUTCMinutes(), 2, 0) + ':' + lpad(jsDate.getUTCSeconds(), 2, 0);
      // modified 20150831 - use local date/time
      //date = lpad(jsDate.getFullYear(), 4, 0) + '-' + lpad(jsDate.getMonth() + 1, 2, 0) + '-' + lpad(jsDate.getDate(), 2, 0);
      //time = lpad(jsDate.getHours(), 2, 0) + ':' + lpad(jsDate.getMinutes(), 2, 0) + ':' + lpad(jsDate.getUTCSeconds(), 2, 0);
    }
    else {
      data = data.trim();
      var parts = data.split(this.ISO_SPLITTER_RE) || [];
      date = this.ICCBL_DATE_RE.test(parts[0]) ? parts[0] : '';
      time = date && parts[1] ? parts[1] : this.TIME_RE.test(parts[0]) ? parts[0] : '';
    }
    // FIXME: review this 
    var DDMMYYYY = this.ICCBL_DATE_RE.exec(date) || [];
    var HHmmssSSS = this.TIME_RE.exec(time) || [];

    if (validate) {
      if (this.includeDate && _.isUndefined(DDMMYYYY[0])) return;
      if (this.includeTime && _.isUndefined(HHmmssSSS[0])) return;
      if (!this.includeDate && date) return;
      if (!this.includeTime && time) return;
    }

    var jsDate = new Date(Date.UTC(DDMMYYYY[3] * 1 || 0,
                                   DDMMYYYY[2] * 1 - 1 || 0,
                                   DDMMYYYY[1] * 1 || 0,
                                   HHmmssSSS[1] * 1 || null,
                                   HHmmssSSS[2] * 1 || null,
                                   HHmmssSSS[3] * 1 || null,
                                   HHmmssSSS[5] * 1 || null));

    var result = '';

    if (this.includeDate) {
      // Override:
      // result = (
      //    lpad(jsDate.getUTCFullYear(), 4, 0) 
      //    + '-' + lpad(jsDate.getUTCMonth() + 1, 2, 0) 
      //    + '-' + lpad(jsDate.getUTCDate(), 2, 0);
      result = ( getIccblUTCDateString(jsDate) );
      //          lpad(jsDate.getUTCMonth() + 1, 2, 0) 
      //          + '/' + lpad(jsDate.getUTCDate(), 2, 0) 
      //          + '/' + lpad(jsDate.getUTCFullYear(), 4, 0)
      //          );
    }

    if (this.includeTime) {
      result = ( result 
          + (this.includeDate ? 'T' : '') 
          + lpad(jsDate.getUTCHours(), 2, 0) 
          + ':' + lpad(jsDate.getUTCMinutes(), 2, 0) 
          + ':' + lpad(jsDate.getUTCSeconds(), 2, 0)
          );

      if (this.includeMilli) {
        result = result + '.' + lpad(jsDate.getUTCMilliseconds(), 3, 0);
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
    DateCell.__super__.initialize.apply(this, arguments);
    var formatter = this.formatter;
    //   formatter.includeDate = this.includeDate;
    //   formatter.includeTime = this.includeTime;
    //   formatter.includeMilli = this.includeMilli;

    var placeholder = "MM/DD/YYYY";
  
    this.editor = this.editor.extend({
      attributes: _.extend({}, this.editor.prototype.attributes, this.editor.attributes, {
        placeholder: placeholder
      })
    });
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
  console.log(key, 'data_type', data_type, 'display_type', display_type, prop);
  
  var cell_options = prop.display_options || {};
//  var cell_options = {};
//  if(_.has(prop,'display_options') && ! _.isEmpty(prop.display_options)){
//    cell_options = prop['display_options'];
//    cell_options = cell_options.replace(/'/g,'"');
//    try{
//      cell_options = JSON.parse(cell_options);
//    }catch(e){
//      console.log('warn: display_options is not JSON parseable, column: ' +
//          key + ', options: ' + cell_options,e);      
//    }
//  } 
  var edit_type = _.isEmpty(prop.edit_type)?display_type:prop.edit_type.toLowerCase();

  var backgridCellType = 'string';
  var typeMap = {
    'date': Iccbl.DateCell,
    'link': Iccbl.LinkCell,
    'siunit': Iccbl.SIUnitsCell,
    'float': 'Number',
    'decimal': Iccbl.DecimalCell,
    'image': Iccbl.ImageCell
  }
  
  if(_.has(typeMap,display_type)){
    console.log('field', key, display_type, 'typemap',typeMap[display_type])
    var backgridCellType = typeMap[display_type];
    if(!_.isEmpty(cell_options)){
      backgridCellType = backgridCellType.extend(cell_options);
    }
  }else{
    console.log('no special cell type for', display_type, 'data_type',data_type);
  }
  
  column = _.extend(column, {
    'name' : key,
    'label' : prop['title'],
    'description' : prop['description'],
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
    // TODO: this is probably backwards, since it appears SelectCell
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
  
  if (optionalHeaderCell) {
    column['headerCell'] = optionalHeaderCell;
  }else{
    // More specific header cell, if available
    if(data_type == 'string'){
      column['headerCell'] = TextHeaderCell;
    }
    else if(data_type == 'integer'
      || data_type == 'float'
      || data_type == 'decimal' ){
      var headerCell = NumberHeaderCell;
      
      if(display_type == 'siunit'){
        headerCell = SIUnitHeaderCell;
        if(!cell_options.symbol){
          var msg = ('Error constructing SIUnit header cell: ' +
              'missing required "symbol" backgrid_cell_option');
          Iccbl.appModel.error('error: ' + e);
        }
        
      }
      column['headerCell'] = headerCell;
    }
    else if( data_type == 'date'){
      column['headerCell'] = DateHeaderCell;
    }
    else if( data_type == 'boolean'){
      column['headerCell'] = BooleanHeaderCell;
    }

    if( edit_type == 'select' 
      || edit_type == 'multiselect' 
      || edit_type == 'multiselect2' 
      || edit_type == 'multiselect3' 
        
    ){
      column['headerCell'] = SelectorHeaderCell;
    }
    
    // FIXME: should be able to pass initialization options properly
    //if(cell_options){
    //  column['headerCell'] = column['headerCell'].extend(cell_options);
    //}
    if(!_.has(column, 'headerCell')){
      console.log('no special header cell found for dt: ', data_type,'edit_type', edit_type);
    }
    
  }
    
  return column;
};

var createBackgridColModel = Iccbl.createBackgridColModel = 
  function(restFields, _orderStack, _manualIncludes) {
    
  console.log('--createBackgridColModel');
  var manualIncludes = _manualIncludes || [];
  var orderStack = _orderStack || [];

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
        // console.log('field not visible in list view: ' + key)
    }
  });
  
  colModel = new Backgrid.Columns(colModel);
  colModel.comparator = 'order';
  colModel.sort();
  return colModel;
};


return Iccbl;
});
