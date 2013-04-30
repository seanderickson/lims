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


var iccbl = {}; // create a global namespace
/**
 * For working with history.js.  Include this before the datatables javascript file for inclusion in the pagination actions.
 * @param {Object} $
 */
(function ($) {
    $.pushOffsetHistory = function (offset) {
		History.pushState({state:offset,rand:Math.random()}, "?offset="+offset, "?offset="+offset);
    	
		// if(!window.dontpush){
			// console.log('button: ' + parseInt($('a', this).text(),10) + ',' + oSettings._iDisplayStart);
			// window.updatingHistory = true;
			// console.log('---------fnDrawCallback history state before xx: ' + JSON.stringify(History.getState().data) + 
				// ',' + oSettings._iDisplayStart );
			// History.pushState({state:oSettings._iDisplayStart,rand:Math.random()}, null, "?offset="+oSettings._iDisplayStart);
			// console.log('pushed: ' + JSON.stringify(History.getState().data));
			// window.updatingHistory = false;
			// console.log('done pushing');
		// }else{
			// console.log('unset window.dontpush');
			// window.dontpush = false;
		// }
	};
}(jQuery));


iccbl.createColModel = function(fields_from_rest) {
	var colModel = [];
	var i = 0;
	var _total_count = 0;
	for (var field in fields_from_rest){
		if (fields_from_rest.hasOwnProperty(field)) { // filter
			var prop = fields_from_rest[field];
			colModel[i] = { 'field':field, 'sTitle':prop['name'], 'bSortable':true, 'arrayOrder':prop['order'] };
			i++;
		}
	}
	
	colModel.sort(function(a,b){ return a['arrayOrder']-b['arrayOrder'] });
	//var _colWidth = 1/i * _width;
	//console.log('colWidth:' + _colWidth);
	return colModel;
};
	