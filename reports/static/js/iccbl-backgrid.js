/**
 * @summary     ICCBL-lims Backgrid Extension Functions
 * @description Utility Functions for the iccbl-lims
 * @version     0.1
 * @file        iccbl-backgrid.js
 * @author      Sean Erickson
 * @contact     sean_erickson “AT” hms.harvard.edu
 *
 * @copyright Copyright 2013 Harvard University, all rights reserved.
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
 *
 **/
define(['jquery', 'underscore', 'backbone', 'backbone_pageable', 'backgrid', 
        'backgrid_filter', 'backgrid_paginator', 'backgrid_select_all', 'lunr',
        // TODO: lunr should not be a requirement - it is for client side filtering,
        // and backgrid_filter is requiring it.
        'text!templates/generic-selector.html'], 
    function($, _, Backbone, BackbonePageableCollection, Backgrid, 
		BackgridFilter, BackgridPaginator, BackgridSelectAll, lunr, 
		genericSelectorTemplate) {

    // Attach PageableCollection to the right place on the Backbone object
    // for compatibility with require.js
	// see https://github.com/wyuenho/backbone-pageable/issues/62
    Backbone.PageableCollection = BackbonePageableCollection;

    var root = window;
    
    var Iccbl = root.Iccbl = {
        VERSION : "0.0.1",

        // Extension: {},
        //
        // requireOptions: function (options, requireOptionKeys) {
        // for (var i = 0; i < requireOptionKeys.length; i++) {
        // var key = requireOptionKeys[i];
        // if (_.isUndefined(options[key])) {
        // throw new TypeError("'" + key  + "' is required");
        // }
        // }
        // },
        //
        // resolveNameToClass: function (name) {
        // if (_.isString(name)) {
        // // var key = _.map(name.split('-'), function (e) { return
        // capitalize(e); }).join('') + suffix;
        // var klass = Iccbl[key]; // || Backgrid.Extension[key];
        // if (_.isUndefined(klass)) {
        // throw new ReferenceError("Class '" + key + "' not found");
        // }
        // return klass;
        // }
        //
        // return name;
        // }
    };

    var assertIccbl = Iccbl.assert = function(condition, message) {
        if (!condition) {
            throw message || "Assertion failed";
        }
    };

    requireOptions = Iccbl.requireOptions = function (options, requireOptionKeys) {
        for (var i = 0; i < requireOptionKeys.length; i++) {
          var key = requireOptionKeys[i];
          if (_.isUndefined(options[key])) {
            throw new TypeError("'" + key  + "' is required");
          }
        }
    };

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

    var getKey = Iccbl.getKey = function(options) {
        var route_fragment = '';
        if (_.isString(options)) {
            route_fragment += options;
        } else if (_.isArray(options)) {// an array is just a set of keys, to be
            // separated by slashes
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

    var getIdFromIdAttribute = Iccbl.getIdFromIdAttribute = function(model, schema){
      // Used to set the id specifically on the model:
      // backbone requires this to determine whether a "POST" or "PATCH" will be used

      if(_.has(schema['resource_definition'], 'id_attribute')){
          var id_attribute = schema['resource_definition']['id_attribute'];
          console.log('create id from ' + id_attribute);
          var id = _.reduce(id_attribute, function(memo, item){
            if(!_.isEmpty(memo)) memo += '/';
            return memo += model.get(item);
          }, '');
          return id;
      }else{
        throw new TypeError("'id_attribute' not found on the schema: " + JSON.stringify(schema)
          + ', for the model: ' + JSON.stringify(model.attributes));
      }
    };

    var getTitleFromTitleAttribute = Iccbl.getTitleFromTitleAttribute = function(model, schema){
      // Used to set the id specifically on the model:
      // backbone requires this to determine whether a "POST" or "PATCH" will be used

      if(_.has(schema['resource_definition'], 'title_attribute')){
          var title_attribute = schema['resource_definition']['title_attribute'];
          console.log('create title from ' + title_attribute);
          var title = _.reduce(schema['resource_definition']['title_attribute'],
              function(memo, item){
                  if( model.has(item) ) memo += model.get(item)
                  else memo += item
                  return memo ;
              }, '');
          console.log('extracted title: ' + title + ' from ' + model);
          return title;
      }else{
        throw new TypeError("'title_attribute' not found on the schema: " + JSON.stringify(schema)
          + ', for the model: ' + JSON.stringify(model.attributes));
      }
    };

    /**
     * Similar the the contains function, but using item.indexOf(matchString) || matchString.indexOf(item)
     * for the truth test.
     */
    var containsByMatch = Iccbl.containsByMatch = function(collection, matchstring){
      return _.find(collection, function(item) {
          var result =  ( matchstring.indexOf(item) != -1 ) || ( item.indexOf(matchstring) != -1 );

          console.log('containsByMatch: ' + result + ', ' + matchstring + ', ' + JSON.stringify(collection));
          return result;
      });
    };

    // var getKey = Iccbl.getKey = function( key_array ){
    // var _key = key_array;
    // Iccbl.assert( !_.isEmpty(_key), 'content:detail: current_options must be
    // defined (as the key), if not schemaResult, model supplied');
    // // handle composite keys
    // if(_.isArray(_key)){
    // _key = _.reduce(_key, function(memo, item){
    // if(!_.isNull(item)) memo += item + '/';
    // return memo;
    // }, '');
    // }
    // return _key;
    // };

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

    var getModel = Iccbl.getModel = function(schemaResult, url, callback) {
        var ModelClass = Backbone.Model.extend({
            url : url,
            defaults : {}
        });
        var instance = new ModelClass();
        instance.fetch({
            success : function(model) {
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
     * Note use this version if the model will be updated.
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
     *
     * @param {Object} fields_from_rest - hash of fields for the current dataset:
     *      field properties { visibility: [array of strings], title: a label for
     * the field, order: display order of the field }
     * @param {Object} optionalHeaderCell - a Backgrid.HeaderCell to use for each
     * column
     * @param {Object} options - a hash of { fieldKey: [custom cell: extend
     * Backgrid.Cell] } to map custom cell implementations to fields
     */
    var createBackgridColModel = Iccbl.createBackgridColModel = 
    	function(restFields, optionalHeaderCell, orderHash) {
        
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
                        //                            console.log('got  ' +
                        // klass);
                        if (!_.isUndefined(klass)) {
                            console.log('----- class found: ' + klass);
                            backgridCellType = klass;
                        }
                    } catch(ex) {
                        var msg = '----for: field: ' + key + ', no Iccbl class found for type: ' + prop['backgrid_cell_type'] + ', this may be a Backgrid cell type';
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
                if(orderHash && _.has(orderHash, key)){
                	colModel[i]['direction'] = 'ascending';
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

    var GenericSelector = Iccbl.GenericSelector = Backbone.View.extend({

        initialize : function(attributes, options) {

            this.listenTo(this.model, 'change', this.changeNotified);
            this._options = options;
            _.bindAll(this, 'changeNotified');
        },

        events : {
            "change #generic_selector" : "selectorChanged"
        },

        selectorChanged : function(e) {
            e.preventDefault();
            var option = e.currentTarget.value;
            this.model.set({
                'selection' : option
            });
        },

        changeNotified : function() {
            var selection = this.model.get('selection');
            console.log('change notified for ' + this._options.label + ', selection: ' + selection + ', options: ' + JSON.stringify(this._options.options));
            this.$('#generic_selector').val(String(selection));
        },

        render : function() {
            this.$el.empty();
            if (!_.contains(this._options.options, '')) {
                this._options.options.unshift('');
                // create a blank entry
            }
            this.$el.append(_.template(genericSelectorTemplate, {
                label : this._options.label,
                'options' : _(this._options.options)
            }));
            if (!_.isUndefined(this._options.selectorClass)) {
                this.$('#generic_selector').removeClass().addClass(this._options.selectorClass);
            };
            this.changeNotified();
            this.delegateEvents();
            return this;
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

var CollectionOnClient = Iccbl.CollectionOnClient = Backbone.PageableCollection.extend({
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

var MyCollection = Iccbl.MyCollection = Backbone.PageableCollection.extend({

    initialize : function(options) {
        var self = this;
        this.options = options;

        Iccbl.requireOptions(options, ['url','listModel']);

        this.url = options.url;
        this.listModel = options.listModel;

//        this.listenTo(this.listModel, 'change:search', function() {
//            console.log('===--- collection detect: listModel change:search');
//            var searchHash = self.listModel.get('search');
//            self.setSearch(searchHash);
//        });

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
    // the server response.
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
            // "order_by=(-)field_name"
            if ( typeof this.state !== 'undefined' && this.state.sortKey && this.state.order) {
                var dir = "-";
                if (this.state.order < 0) {// according to docs, -1 == ascending
                    dir = "";
                }
                return dir + this.state.sortKey;
            }
        },
        directions : {
            "-1" : "asc",
            "1" : "desc"
        }
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
    setSearch : function(searchHash) {
        var self = this;
        var searchHash = _.clone(searchHash);
        console.log('collection.setSearch: trigger search to headers:' + 
        		JSON.stringify(searchHash));
        this.trigger("MyServerSideFilter:search", searchHash, this);

        // if the search key is not in the queryParams, then it is not a column
        // search.
        // this will add it manually to the queryParams (which are serialized in
        // the fetch to the server)
        _.each(_.keys(searchHash), function(key) {
            var val = searchHash[key]

            if(_.isEmpty(val)){
                delete self.queryParams[key];
            }else{
                // check if param dne, or if param exists and is a value to be set
                // the reason for the "isFunction" check is that the Backgrid-filter
                // defined params are function calls to get the current value in the
                // searchbox
                if (!_.has(self.queryParams, key) || !_.isFunction(self.queryParams[key])) {
                	var _data = {};
                	_data[key]=val;
                	collection.getFirstPage({data:_data, reset: true, fetch: true});

//                    self.queryParams[key] = val;
//                    // NOTE: setting the key on the queryParams means that it will be
//                    // sent to the server on the next fetch
                }
            }
        });

    },

    // Proxy for the search elements to add search terms to the listModel
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

    // Proxy for the search elements to clear search terms from the listModel
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
     *  Override - called from HeaderCells
     */
    setSorting : function(sortKey, order, options) {
    	Backbone.PageableCollection.prototype.setSorting.call(
			  this, sortKey, order);
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
 * TODO: probably going to have to implement full, instead of extending.  Too 
 * many side effects of extending ServerSideFilter
 **/
var MyServerSideFilter = 
		Iccbl.MyServerSideFilter = 
		Backgrid.Extension.ServerSideFilter.extend({
    columnName : null, // TODO: use "name"
    className: "form-search",    
    template: _.template('<input type="search" <% if (placeholder) { %> placeholder="<%- placeholder %>" <% } %> name="<%- name %>" /><a class="backgrid-filter clear" data-backgrid-action="clear" href="#">&times;</a>', null, {variable: null}),

    initialize : function(options) {
        this.columnName = options.columnName;
        Backgrid.Extension.ServerSideFilter.prototype.initialize.apply(
        		this, [options]);
    },
    
    /**
     * Override
     */
    search: function(e){
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
 * Override to add Search capability.
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

//        this._serverSideFilter['events']['click .close'] = function(e) {
//            if(e) e.preventDefault();
//
//            this.remove();
//            this.clear();
//            this.collection.clearSearch([this.name]);
//        };

//        this._serverSideFilter['events']['submit'] = function(e) {
//            var searchHash = {};
//            searchHash[this.name] = this.searchBox().val();
//            console.log('server side filter add search: ' + 
//            		JSON.stringify(searchHash));
//            this.collection.addSearch(searchHash);
//            this.search(e);
//        };
    },


    
    remove : function() {
        console.log('headercell remove called');
        this._serverSideFilter.remove();
        this._serverSideFilter.unbind();
        this._serverSideFilter.collection = null;
        this.unbind();

        Backgrid.HeaderCell.prototype.remove.apply(this);
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
//                    self._serverSideFilter.$el.find("input[type=text]").val(val);
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
        var filterIcon = $('<i class="icon-search"></i>');
        filterIcon.click(function(e) {
            _handle.$el.append(_handle._serverSideFilter.render().el);
        });
        this.$el.append(filterIcon);

        this.$el.prop('title', 
        			  this.options['column']['attributes']["description"]);
        return this;
    }
});

return Iccbl;
});
