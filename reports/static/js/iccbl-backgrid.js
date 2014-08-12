/**
 * @summary     ICCBL-lims Backgrid Extension Functions
 * @description Utility Functions for the iccbl-lims
 * @version     0.1
 * @file        iccbl-backgrid.js
 * @author      Sean Erickson
 * @contact     sean_erickson “AT” hms.harvard.edu
 *
 * @copyright 2014 Harvard University, all rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This source file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license files for details.
 *
 * For details please refer to: https://github.com/hmsiccbl/lims
 **/

// FIXME: lunr should not be a requirement - for server side filtering,
// and backgrid_filter is requiring it.

define(['jquery', 'underscore', 'backbone', 'backbone_pageable', 'backgrid', 
        'backgrid_filter', 'backgrid_paginator', 'backgrid_select_all', 'lunr',
        'layoutmanager'],
    function($, _, Backbone, BackbonePageableCollection, Backgrid, 
             BackgridFilter, BackgridPaginator, BackgridSelectAll, lunr,
             layoutmanager ) {

  // Attach PageableCollection to the right place on the Backbone object
  // for compatibility with require.js
	// see https://github.com/wyuenho/backbone-pageable/issues/62
  Backbone.PageableCollection = BackbonePageableCollection;

  var root = window;
  
  var Iccbl = root.Iccbl = {
      VERSION : "0.0.1",
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
   * Replace all {tokens} in the string with model attributes.
   * - fallback to "token" if model attribute is not set.
   */
  var replaceTokens = Iccbl.replaceTokens = function(model,stringWithTokens) {
    var interpolatedString = stringWithTokens.replace(/{([^}]+)}/g, 
      function (match) 
      {
        match = match.replace(/[{}]/g,'');
        return !_.isUndefined(model.get(match)) ? model.get(match) : match;
      });
    return interpolatedString;
  }
  
  var stringToFunction = Iccbl.stringToFunction = function(str) {
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
   * Get the key off the URI stack:
   * - we don't know if the key is composite or not; so we don't know how many 
   * items to pop off the stack; So look at the resource definition 
   * id_attribute; which lists the keys.
   * @param resource - a resource definition as defined by the API
   * @param urlStack - array representation of the current unprocessed URI 
   * elements.
   * @param consumedStack - holds the items popped off the stack
   */
  var popKeyFromStack = Iccbl.popKeyFromStack = function(resource, urlStack, consumedStack){
    var id  = '';
    _.each(resource.id_attribute, function(attribute){
      if (_.isEmpty(urlStack)){
        var msg = 'not enough items on the URL to create the key for resource: ' + 
            resource.title;
        window.alert(msg);
        throw msg;
      }
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
    return id;
  };
  
  var getIdKeys = Iccbl.getIdKeys = function(model,schema) {
    if (_.has(schema['resource_definition'], 'id_attribute')) {
      var id_attribute = schema['resource_definition']['id_attribute'];
      console.log('create id from ' + id_attribute);
      var idList = [];
      _.each(id_attribute, function(item){
        idList.push(model.get(item));
      });
      return idList;
    } else {
      throw new TypeError("'id_attribute' not found on the schema: " 
              + JSON.stringify(schema)
              + ', for the model: ' + JSON.stringify(model.attributes));
    }
    
  };
  
  /**
   * Create a string ID from the 'id_attribute' of a schema resource definition.
   * - the id_attribute is an array of field specifiers
   * - the 'ID' will be each of these fields, concatenated with a forward slash.
   * Note: the composite key is the public key; although there is often an
   * 'ID' field, it should not be used for URL's, as it is a internal, possibly
   * transient implementation detail.
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
   * definition.
   * - the title_attribute is an array of field specifiers and strings.
   * - if an array item is a field, the field value will be used,
   * - if an array item is not a field, then it will be concatenated directly.
   */
  var getTitleFromTitleAttribute = Iccbl.getTitleFromTitleAttribute = 
      function(model, schema){
  
    if(_.has(schema['resource_definition'], 'title_attribute')){
      var title_attribute = schema['resource_definition']['title_attribute'];
      console.log('create title from ' + title_attribute);
      var title = _.reduce(
        schema['resource_definition']['title_attribute'],
        function(memo, item){
          if( model.has(item) ) memo += model.get(item)
          else memo += item
          return memo ;
        }, '');
      console.log('extracted title: ' + title + ' from ' + model);
      return title;
    }else{
      throw new TypeError("'title_attribute' not found on the schema: " + 
          JSON.stringify(schema)
          + ', for the model: ' + JSON.stringify(model.attributes));
    }
  };

  /**
   * Matches array items against the matchstring.  
   * - Matches from the right to left; allowing URI fragments to match their 
   * parent URIs.
   * Similar the the contains function, 
   * but using item.indexOf(matchString) || matchString.indexOf(item)
   * for the truth test.
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

  var getSchema = Iccbl.getSchema = function(schema_url, callback) {
    $.ajax({
      type : "GET",
      url : schema_url, //options.url_schema,
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
            //console.log('error fetching the model: '+ model + ', response:
            // ' + JSON.stringify(response));
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
   * Note use this version if the model will be updated to the server.
   * - Backbone will use "put" if the model has an ID set
   * - also, use model.save([attribute-list], {patch: true} ) to do incremental
   * updates,
   * see http://backbonejs.org/#Model-save
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
              //console.log('error fetching the model: '+ model + ', response:
              // ' + JSON.stringify(response));
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

  var getCollectionOnClient = Iccbl.getCollectionOnClient = function(url, callback){
    var CollectionClass = Iccbl.CollectionOnClient.extend({
      url: url 
    });
    var instance = new CollectionClass();
    instance.fetch({
      data: { limit: 0 },
      success: function(collection, response) {
        console.log('success callback...');
        callback(collection);
      },
      error: function(model, response, options) {
          //console.log('error fetching the model: '+ model + ', response:
          // ' + JSON.stringify(response));
          var msg = 'Error locating resource: ' + url;
          var sep = '\n';
          if (!_.isUndefined(response.status))
              msg += sep + response.status;
          if (!_.isUndefined(response.statusText))
              msg += sep + response.statusText;
          if (!_.isEmpty(response.responseText))
              msg += sep + response.responseText;
          window.alert(msg);
          // TODO: 1. use Bootstrap inscreen alert classed message div
          // TODO: 2. jquery seems to swallow json parsing exceptions, fyi
          throw msg;
      },
      always: function(){
        console.log('done: ');
      }
    });
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
              //console.log('error fetching the model: '+ model + ', response:
              // ' + JSON.stringify(response));
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
              //console.log('error fetching the model: '+ model + ', response:
              // ' + JSON.stringify(response));
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
   * Note use this version if the model will be updated.
   * - Backbone will use "put" if the model has an ID set
   * - also, use model.save([attribute-list], {patch: true} ) to do incremental
   * updates,
   * see http://backbonejs.org/#Model-save
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

  /**
   * Return an array for backgrid column descriptors.
   * @param {Object} fields_from_rest - hash of fields for the current dataset:
   *    field properties { visibility: [array of strings], title: a label for
   *    the field, order: display order of the field }
   * @param {Object} optionalHeaderCell - a Backgrid.HeaderCell to use for each
   *    column
   * @param {Object} options - a hash of { fieldKey: [custom cell: extend
   *    Backgrid.Cell] } to map custom cell implementations to fields
   */
  var createBackgridColModel = Iccbl.createBackgridColModel = 
  	    function(restFields, optionalHeaderCell, orderStack) {
      
  	  console.log('--createBackgridColModel');
      //: restFields: ' + JSON.stringify(restFields));
      var colModel = [];
      var i = 0;
      var _total_count = 0;
      _.each(_.pairs(restFields), function(pair) {
          var key = pair[0];
          var prop = pair[1];

          var visible = _.has(pair[1], 'visibility') && _.contains(pair[1]['visibility'], 'list');
          if (visible) {

              var backgridCellType = 'string';
              if (!_.isEmpty(prop['backgrid_cell_type'])) {
                  backgridCellType = prop['backgrid_cell_type'];
                  try {
                      //                            console.log('look for ' +
                      // key + ', ' + prop['backgrid_cell_type']);
                      var klass = Iccbl.stringToFunction(prop['backgrid_cell_type']);
                      if (!_.isUndefined(klass)) {
                          console.log('----- class found: ' + prop['backgrid_cell_type']);
                          backgridCellType = klass;
                      }
                  } catch(ex) {
                      var msg = ('----for: field: ' + key + 
                          ', no Iccbl class found for type: ' + 
                          prop['backgrid_cell_type'] + 
                          ', this may be a Backgrid cell type');
                      console.log(msg + ': ' + JSON.stringify(ex));
                  }
              }
              
              colModel[i] = {
                  'name' : key,
                  'label' : prop['title'],
                  'description' : prop['description'],
                  cell : backgridCellType,
                  order : prop['ordinal'],
                  editable : false,
              };
              if(orderStack && _.contains(orderStack, key)){
                colModel[i]['direction'] = 'ascending';
              }
              else if(orderStack && _.contains(orderStack, '-' + key)){
                colModel[i]['direction'] = 'descending';
              }
              if (optionalHeaderCell) {
                  colModel[i]['headerCell'] = optionalHeaderCell;
              }
              i++;
          } else {
              //console.log('field not visible in list view: ' + key)
          }
      });

      colModel.sort(function(a, b) {
          if (_.isNumber(a['order']) && _.isNumber(b['order'])) {
              return a['order'] - b['order'];
          } else if (_.isNumber(a['order'])) {
              return -1;
          } else if (_.isNumber(b['order'])) {
              return 1;
          } else {
              return 0;
          }
      });
      return colModel;
  };

  var MyModel = Iccbl.MyModel = Backbone.Model.extend({
      // TODO: we want to make sure there is a trailing slash, or tastypie
      // doesn't work.
      url : function() {
          var url = Backbone.Model.prototype.url.call(this);
          // console.log('---- url1: ' + url);
          return url + (url.charAt(url.length - 1) === '/' ? '' : '/');
      },

      initialize : function() {
          Backbone.Model.prototype.initialize.apply(this, arguments);
          // console.log('x--- urlRoot: ' + this.urlRoot + ", " + this.id + ',
          // ' + this.collection.url);
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
          // be attached to the prototype manually here.  why?
          this.url = MyModel.prototype.url;
      },
  });


  var LinkCell = Iccbl.LinkCell = Backgrid.Cell.extend({
      className : "link-cell",
      events : {
          "click #link" : "toLink",
      },

      initialize : function(options) {
          Backgrid.Cell.prototype.initialize.apply(this, arguments);
      },

      render : function() {
          this.$el.empty();
          var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
          this.$el.append($("<a id='link' >", {
              tabIndex : -1,
              href : '',
              title : formattedValue,
              //target : "_blank"
          }).text(formattedValue));

          this.delegateEvents();
          return this;
      },

      toLink : function(e) {
          e.preventDefault();
          this.model.collection.trigger("MyCollection:link", this.model, this.column.get("name"));
      },
  });

  var EditCell = Iccbl.EditCell = Backgrid.Cell.extend({
      className : "detail-cell",
      events : {
          "click #edit" : "editDetail",
      },

      initialize : function(options) {
          Backgrid.Cell.prototype.initialize.apply(this, arguments);
      },

      render : function() {
          this.$el.empty();
          var formattedValue = this.formatter.fromRaw(this.model.get(this.column.get("name")));
          this.$el.append($("<a id='edit' >", {
              tabIndex : -1,
              href : '',
              title : formattedValue,
              //target : "_blank"
          }).text(formattedValue));

          this.delegateEvents();
          return this;
      },

      editDetail : function(e) {
          e.preventDefault();
          this.model.collection.trigger("MyCollection:detail", this.model);
      },
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
    this.$el.append($("<a id='delete' >", {
    tabIndex : -1,
    href : '',
    }).text(this.options.column.attributes['text']));
    this.delegateEvents();
    return this;
    },
  
    delete: function(e) {
        e.preventDefault();
        this.model.collection.trigger("MyCollection:delete", this.model);
    }
  });

var CollectionOnClient = Iccbl.CollectionOnClient = Backbone.Collection.extend({
    /**
     *  Override collection parse method:
     *      Parse server response data.
     */
    parse : function(response) {
        console.log('Collection on client, parse called');
        // hack the response for tastypie:
        return response.objects;
    },
});


var CollectionInColumns = Iccbl.CollectionInColumns = Backbone.Collection.extend({
    /**
     *  Override collection parse method:
     *      Parse server response data.
     * untested
     */
    parse : function(response) {
        console.log('Collection on client, parse called');
        // hack the response for tastypie:
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
    console.log('UriContainerView initialize');
    var model = this.model = args.model;
    var targetProperty = args.property || 'uriStack';
    this.listenTo(model, 'change:'+targetProperty , this.uriStackChange );
    
    Backbone.View.prototype.initialize.apply(this,arguments);
  },
  
  /**
   * Child view bubble up URI stack change event
   */
  reportUriStack: function(reportedUriStack, options) {
    var options = options || {source: this};
    var consumedStack = this.consumedStack || [];
    var actualStack = consumedStack.concat(reportedUriStack);
    this.model.set({'uriStack': actualStack}, options);     
  },
  
  /**
   * Backbone.Model change event handler
   * @param options.source = the event source triggering view
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
        // FIXME: better global error handling
        console.log('error thrown' + e);
        // FIXME: why is the ajaxStop handler in main.js not being called anyway?
        $('#loading').fadeOut({duration:100});
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
   * See Backgrid.Body.sort:
   * - created to solve the multisort case for the server side backbone-pageable
   * collection only.
   */
  sort: function (column, direction) {
    var collection = this.collection;

    var order;
    if (direction === "ascending") order = -1;
    else if (direction === "descending") order = 1;
    else order = null;
    
    // Replaces:    
    //    collection.setSorting(order && column.get("name"), order,
    //        {sortValue: column.sortValue()});
    collection.setSorting(column.get("name"), order,
        {sortValue: column.sortValue()});

    collection.fetch({reset: true, success: function () {
      collection.trigger("backgrid:sorted", column, direction, collection);
    }});
    
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
          //          var val = this.state.orderStack.join('&order_by=');
          return this.state.orderStack;
        }
      }
      //      directions : {
      //          "-1" : "asc",
      //          "1" : "desc"
      //      }
  },

  /**
   *  Override 
   */
  parseState : function(response, queryParams, state, options) {
    // hack the response for tastypie:
    // note, this is because the pageable collection doesn't work with the
    // backbone-tastypie.js fix
    var state = _.clone(state);
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
   *  Override 
   */
  parseRecords: function (resp, options) {
      return resp.objects;
  },
      
  /**
   *  Method for external callers to set the search
   */
  setSearch: function(searchHash) {
    var self = this;
    var searchHash = _.clone(searchHash);
    
    // Tell all the header cells
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
      self.fetch({data:_.clone(_data), reset: true});
    }else{
      self.fetch();
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
  },

  /**
   * Proxy for the search elements to clear search terms from the listModel
   */ 
  clearSearch: function(searchKeys) {
    console.log('clearsearch: ' + JSON.stringify(searchKeys));
    var self = this;
    var searchHash = {};
    if (!_.isUndefined(searchKeys)) {
        searchHash = _.clone(self.listModel.get('search'));
        _.each(searchKeys, function(searchKey) {
            delete searchHash[searchKey];
        });
    }
    self.listModel.set({
        'search' : searchHash
    });
  },

  /**
   *  Override - 
   *  HeaderCell.onClick -> backgrid:sort -> body.sort -> 
   *    -> (BackbonePageable)Collection.setSorting
   *    -> Collection.fetch() -> grab data from state based on "queryParams"
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
        }else if (newdir != dir){
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
    
    //  	Backbone.PageableCollection.prototype.setSorting.call(this, sortKey, order);
	  // TODO: Investigate why PageableCollection.setSorting is not triggering 
	  // a 'sort' event (needed to clear old sort indicators).
	  // Sequence of a sort:
	  // Backgrid.HeaderCell.onClick-> collection.trigger('backgrid:sort')
	  // 	-> 'backgrid:sort' -> Backgrid.Body.sort()
	  //    	-> PageableCollection.setSorting(): sets state.sortKey
	  //		-> if(fullCollection) (client mode) collection.sort()
	  //		-> else
	  //		-> Collection.fetch(reset:true) 
	  //			*so in this case, no sort(), if reset:false, then a "set" 
	  //				would be called, and a 'sort' triggered
	  //		* also calls column.set to put the new indicator 
	  //		column.set('direction')
	  // without a sort, there is no erasing of the old sort indicators:
	  // Collection.sort() -> trigger('sort') ->
	  //		-> Backgrid.HeaderCell.removeCellDirection
	  // Last note: this may be caused by not getting the sortKey from the 
	  // queryParams on parseState.
    	this.trigger('sort',this); 
  },

});


/**
 * Override so that we can keep a handle to the containing column name.
 * TODO: can handle this with events instead (so that the filter notifies the
 * containing headercell?)
 * 
 **/
var MyServerSideFilter = 
		  Iccbl.MyServerSideFilter = 
		  Backgrid.Extension.ServerSideFilter.extend({
    columnName : null, // TODO: use "name"
    
    // override
    className: "form-search",    
    
    // override, provide our own class
    template: _.template('<input type="search" <% if (placeholder) { %> placeholder="<%- placeholder %>" <% } %> name="<%- name %>" /><a class="backgrid-filter clear" data-backgrid-action="clear" href="#">&times;</a>', null, {variable: null}),

    initialize : function(options) {
      this.columnName = options.columnName;
      Backgrid.Extension.ServerSideFilter.prototype.initialize.apply(
      		this, [options]);
    },
    
    /**
     * Override - to allow fetch without returning to the first page.
     * (note: this has negative side effect that if the page is too high, a 
     * blank list is shown).
     */
    search: function(e) {
      if (e) e.preventDefault();    	
    	var searchHash = {};
    	searchHash[this.name] = this.searchBox().val();
    	console.log('server side filter add search: ' + 
    			JSON.stringify(searchHash));
    	this.collection.addSearch(searchHash);
    	
    	
    	// Note: not calling super, because it will force first page.
    	//    	Backgrid.Extension.ServerSideFilter.prototype.search.apply(this,e);
      if (e) e.preventDefault();

      var data = {};
      var query = this.searchBox().val();
      if (query) data[this.name] = query;

      var collection = this.collection;
      // We're overriding this behaviour:
  		//        // go back to the first page on search
  		//        if (Backbone.PageableCollection &&
  		//            collection instanceof Backbone.PageableCollection) {
  		//          collection.getFirstPage({data: data, reset: true, fetch: true});
  		//        }
      collection.fetch({data: data, reset: true});
      
      this.showClearButtonMaybe();
    },
    
    /**
     * Override
     */
    clear: function(e){
      if (e) e.preventDefault();
      this.remove();
      this.collection.clearSearch([this.name]); 
      Backgrid.Extension.ServerSideFilter.prototype.clear.apply(this,e);
    },    
});

/**
 * Override of the Backgrid.HeaderCell to:
 * - add search
 * - add multisort
 * 
 * TODO: this should be named "contains" header cell, because it searches using
 * "__contains"
 **/
var MyHeaderCell = Iccbl.MyHeaderCell = Backgrid.HeaderCell.extend({

    _serverSideFilter : null,

    initialize : function(options) {
        this.options = options;
        Backgrid.HeaderCell.prototype.initialize.apply(this, [options]);

        this._serverSideFilter = new MyServerSideFilter({
            // TODO: collection should be passed as an option
            collection : this.collection, 
            // Name of the URL query parameter for tastypie/django 
            // TODO: support other filters
            name : this.column.get("name") + "__contains", 
            // HTML5 placeholder for the search box
            placeholder : "Search " + this.column.get("label"), 
            columnName : this.column.get("name"),
        });

        this.listenTo(this.collection,"MyServerSideFilter:search",this._search);
        this.listenTo(this.collection,"sort",this.collectionSorted);
    },

    remove: function() {
        console.log('headercell remove called');
        this._serverSideFilter.remove();
        this._serverSideFilter.unbind();
        this._serverSideFilter.collection = null;
        this.unbind();

        Backgrid.HeaderCell.prototype.remove.apply(this);
    },
    
    collectionSorted: function(collection, options){
      var self = this;
      var name = this.column.get('name');
      var state = this.collection.state;
      if(_.isUndefined(state.orderStack)){
        console.log('Warning: grid does not support the ordered sort functionality');
        return;
      }
      var i = 0;
      _.each(state.orderStack, function(order_entry){
        i++;
        var dir = order_entry.substring(0,1);
        var direction = null;
        var fieldname = order_entry;
        if(dir == '-'){
          fieldname = order_entry.substring(1);
          direction = 'ascending';
        }else{
          dir = '';
          direction = 'descending';
        }
        if(fieldname == name){
          self.$el.removeClass("ascending").removeClass("descending");
          self.$el.addClass(direction);

          var sorter = self.$el.find('#sorter');
          sorter.empty();
          sorter.append("<span class='badge pull-right'>" + i + "<b class='sort-caret'></b></span>");
        }
      });
    },
    
    /**
     Event handler for the collection's `sort` event. Removes all the CSS
     direction classes.
    */
    removeCellDirection: function () {
      //   this.$el.removeClass("ascending").removeClass("descending");
      //   this.column.set("direction", null);
    },

    /**
      Event handler for the column's `change:direction` event. If this
      HeaderCell's column is being sorted on, it applies the direction given as a
      CSS class to the header cell. Removes all the CSS direction classes
      otherwise.
    */
    setCellDirection: function (column, direction) {
      if(! direction){
        this.$el.removeClass("ascending").removeClass("descending");
        this.$el.find("#sorter").empty();
      }
    },
    
    /**
     * Function to listen for router generated custom search event
     * "MyServerSideFilter:search"
     * - add search term and show box
     * - clear old search terms
     **/
    _search: function(searchHash, collection){
        var self = this;
        if (collection == this.collection){
            var found = false;
            _.each(_(searchHash).pairs(), function(pair){
                var key = pair[0];
                var val = pair[1];
                if (self._serverSideFilter.name == key){
                    console.log('--found search: ' + key + '=' + val + 
                    		', on: ' + self.column.get('name'));
                    found = true;
                    // create the DOM element
                    self.$el.append(self._serverSideFilter.render().el);
                    // set the search term
                    self._serverSideFilter.searchBox().val(val);
                    // the filter search method will call collection.fetch()
                    self._serverSideFilter.search(); 
                }
            });

            if (!found) {
                if (!_.isEmpty(this.$el.find("input[type=text]").val())) {
                    this._serverSideFilter.remove();
                    this.collection.clearSearch([self._serverSideFilter.name]);
                }
            }
        }
    },

    /**
     * Renders a header cell with a sorter and a label.
     */
    render : function() {
        Backgrid.HeaderCell.prototype.render.apply(this);
        var _handle = this;
        var filterIcon = $('<span class="glyphicon glyphicon-search"></span>');
        filterIcon.click(function(e) {
            _handle.$el.append(_handle._serverSideFilter.render().el);
        });
        this.$el.append(filterIcon);

        this.$el.prop('title', 
        			  this.options['column']['attributes']["description"]);
        
        // modify the sorting to use our badge-sort multisort 
        var column = this.column;
        var sortable = Backgrid.callByNeed(column.sortable(), column, this.collection);
        if(sortable){
          var sort = this.$el.find(".sort-caret");
          var parent = sort.parent();
          sort.remove();
          parent.append("<div id='sorter'></div>");
        }
        
        return this;
    }
});

return Iccbl;
});
