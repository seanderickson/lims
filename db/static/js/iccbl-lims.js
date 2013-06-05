/**
 * @summary     ICCBL-lims Utility Functions
 * @description Utility Functions for the iccbl-lims
 * @version     0.1
 * @file        iccbl-lims.js
 * @author      Sean Erickson
 * @contact     sean_erickson “AT”hms.harvard.edu
 *
 * @copyright Copyright 2013 Harvard University, all rights reserved.
 *
 * This source file is free software, under either the GPL v2 license
 *
 * This source file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the license files for details.
 *
 * For details please refer to: https://github.com/hmsiccbl/lims
 **/


var iccbl = {
    debug: true,
}; // create a global namespace
// /**
 // * For working with history.js.  Include this before the datatables javascript file for inclusion in the pagination actions.
 // * @param {Object} $
 // */
// (function ($) {
    // $.pushOffsetHistory = function (offset) {
		// History.pushState({state:offset,rand:Math.random()}, "?offset="+offset, "?offset="+offset);
//
		// // if(!window.dontpush){
			// // console.log('button: ' + parseInt($('a', this).text(),10) + ',' + oSettings._iDisplayStart);
			// // window.updatingHistory = true;
			// // console.log('---------fnDrawCallback history state before xx: ' + JSON.stringify(History.getState().data) +
				// // ',' + oSettings._iDisplayStart );
			// // History.pushState({state:oSettings._iDisplayStart,rand:Math.random()}, null, "?offset="+oSettings._iDisplayStart);
			// // console.log('pushed: ' + JSON.stringify(History.getState().data));
			// // window.updatingHistory = false;
			// // console.log('done pushing');
		// // }else{
			// // console.log('unset window.dontpush');
			// // window.dontpush = false;
		// // }
	// };
// }(jQuery));
//
//
// iccbl.createColModel = function(fields_from_rest) {
	// var colModel = [];
	// var i = 0;
	// var _total_count = 0;
	// for (var field in fields_from_rest){
		// if (fields_from_rest.hasOwnProperty(field)) { // filter
			// var prop = fields_from_rest[field];
			// colModel[i] = { 'field':field, 'sTitle':prop['name'], 'bSortable':true, 'arrayOrder':prop['order'] };
			// i++;
		// }
	// }
//
	// colModel.sort(function(a,b){ return a['arrayOrder']-b['arrayOrder'] });
	// //var _colWidth = 1/i * _width;
	// //console.log('colWidth:' + _colWidth);
	// return colModel;
// };

iccbl.createBackgridColModel = function(fields_from_rest, optionalHeaderCell) {
	var colModel = [];
	var i = 0;
	var _total_count = 0;
	for (var field in fields_from_rest){
		if (fields_from_rest.hasOwnProperty(field)) { // filter
			var prop = fields_from_rest[field];
			colModel[i] = { 'name':field, 'label':prop['name'], cell: 'string', order: prop['order']};
			if (optionalHeaderCell){
				colModel[i]['headerCell'] = optionalHeaderCell;
			}
			i++;
		}
	}

	colModel.sort(function(a,b){
	    //console.log('sort: ' + a['order'] +',' + b['order']);
        if(typeof a['order'] !== 'undefined' && typeof b['order'] !== 'undefined'){
            return a['order']-b['order'];
        }else if(typeof a['order'] !== 'undefined'){
            return -1;
        }else if(typeof b['order'] !== 'undefined'){
            return 1;
        }else{
            return 0;
        }
    });
	//console.log('colModel: ' + JSON.stringify(colModel));
	//var _colWidth = 1/i * _width;
	//console.log('colWidth:' + _colWidth);
	return colModel;
};


iccbl.MyModel = Backbone.Model.extend({});

iccbl.MyCollection = Backbone.PageableCollection.extend({
  url: function(){
    return this.instanceUrl;
  },
  initialize: function(props){
    this.instanceUrl = props.url;
  },

      // initialize: function(options){
        // console.log('initialize: ' +  JSON.stringify(options));
        // if (!options || !options.url) {
            // if(typeof iccbl.debug !== 'undefined') console.log( JSON.stringify(options));
            // throw "InvalidConstructArgs";
        // }
        // this.set({ url: options.url });
        // Backbone.PageableCollection.prototype.initialize.call(this, null, options);
    // },

    searchBy: null,
    model: iccbl.MyModel,
//    url: url,
    state: {
        pageSize: 25,
    },
    queryParams: {
        // adjust the query params for tastypie
        pageSize: 'limit',
        offset: function(){
            return (this.state.currentPage-1) * this.state.pageSize;
        },
        totalRecords: null, // unset for tastypie
        totalPages: null, // unset for tastypie
        sortKey: "order_by", // modified for tastypie
        order: null, // unset for tastypie
        order_by: function() { // modified for tastypie: use only "order_by=(-)field_name"
            if (typeof this.state !== 'undefined' && this.state.sortKey && this.state.order) {
                var dir = "-";
                if (this.state.order<0){  // according to docs, -1 == ascending
                    dir = "";
                }
                return dir + this.state.sortKey;
            }
        },

    },
    parse: function(response) {
        // hack the response for tastypie:
        // note, this is because the pageable collection doesn't work with the backbone-tastypie.js fix
        //this.state.totalRecords = response.meta.total_count;
        var state = _.clone(this.state);
        state.totalRecords = response.meta.total_count;
        if(Math.ceil(state.totalRecords/state.pageSize) < state.currentPage ){
            console.log('adjust currentPage');
            state.currentPage=1;
        }
        this.state = this._checkState(state); // recalculate the state and do sanity checks.
        console.log('new state: ' + JSON.stringify(this.state));
        this.setRoutes();
        return response.objects;
    },
    /**
       @property {-1|0|1} [state.order=-1] The order to use for sorting. Specify
       -1 for ascending order or 1 for descending order. If 0, no client side
       sorting will be done and the order query parameter will not be sent to
     */

    setRoutes: function() {
        console.log('setRoutes: ' + this.searchBy + ', '
            + this.state.sortKey + ', ' + this.state.order + ', '+ this.state.pageSize + ', ', + this.state.currentPage);
        var route = '';
        if(this.searchBy !== null){
            route += 'search/'+this.searchBy ;
        }
        if(this.state.sortKey !== null){
            var sortKey = this.state.sortKey;
            if(this.state.order > 0){
                sortKey = '-' + sortKey;
            }
            if(route.length > 0) route += '/';

            route += 'order_by/' + sortKey;
            //router.navigate('order_by/' + sortKey +'/page/'+ this.state.currentPage);
        }
        if(route.length > 0) route += '/';

        route += 'rpp/' + this.state.pageSize + '/page/' + this.state.currentPage;

        // router.navigate('/page/'+ this.state.currentPage);
        router.navigate(route);
    },

    // Override
    setSorting: function(sortKey,order,options) { // override and hack in sorting URL support
        console.log('setSorting called: ' + sortKey+ ', order: ' + order + typeof(order) + ', options: ' + options );
        console.log('searchBy: ' + JSON.stringify(this.data));
        var dir = '-';
        var direction = 'descending';
        if(typeof order !== 'undefined' && order < 0){
            dir = '';
            direction = 'ascending';
        }
        var obj = Backbone.PageableCollection.prototype.setSorting.call(this,sortKey, order);

        //this.setRoutes();

        return obj;
    },

    // Override
    getPage: function(page) { // override and hack in paging URL support
        var obj = Backbone.PageableCollection.prototype.getPage.call(this,page);
        //this.setRoutes();
        return obj;
    },

    fetch: function(options) {
        model = this;
        if (typeof options !== 'undefined') { // TODO: review; find a better way to _always_ stop the ajax spinner
            options.error = function(resp){
                if (error) error(model, resp, options);
                window.App.ajaxComplete();
                window.alert(error);
            };
        }
        return Backbone.PageableCollection.prototype.fetch.call(this,options);
    }


});// end definition of Collection extension

