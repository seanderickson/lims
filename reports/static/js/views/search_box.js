define([
    'jquery',
    'underscore',
    'backbone',
    'backbone_forms',
    'iccbl_backgrid',
    'views/list2',
    'models/app_state',
    'views/generic_edit'
], function($, _, Backbone, backbone_forms, Iccbl, ListView, appModel, EditView) {
    
  
  var SearchView = Backbone.Form.extend({
    
    template: _.template([
        "<form class='iccbl-headerfield-form' ",
        "title='Look up items by plate number(s), plate range(s) and copy name(s);",
        "separate lists by commas and separate lookup types by spaces,",
        " use a newline to separate into subgroups'",
        ">",
        '<div data-fieldsets/>',
        '<div id="data-error" class="has-error" ></div>',
        "</form>"].join('')),

    fieldTemplate: _.template([
        '<div>',
        '<label title="<%= help %>" for="<%= editorId %>"><%= title %></label>',
        '<div>',
        '  <span placeholder="testing.." data-editor></span>',
        '</div>',
        '</div>'
        ].join('')),
    
    choiceFieldTemplate: _.template([
      '<div class="form-group" key="form-group-<%=key%>" >',
      '    <div class="" >',
      '      <div data-editor  key="<%=key%>" style="min-height: 0px; padding-top: 0px; margin-bottom: 0px;" />',
      '      <div data-error class="text-danger" ></div>',
      '      <div><%= help %></div>',
      '    </div>',
      '  </div>',
    ].join('')),
   
    initialize: function(args) {
      
      var self = this;
      console.log('---- initialize search_box');

      var args = args || {};
      var formSchema = this.schema = args.schema || {};
      
      var Select2 = Backbone.Form.editors.Select.extend({
        render: function() {
          Select2.__super__.render.apply(this,arguments);
//          this.$el.css('position', 'relative'); // kludgy; brings it on top
//          this.$el.css('z-index', 10000); // brings it on top
          return this;
        },        
      });
      
      var TextArea2 = Backbone.Form.editors.TextArea.extend({
        render: function() {
          TextArea2.__super__.render.apply(this,arguments);
          //this.$el.css('position', 'relative'); // kludgy; brings it on top
          //this.$el.css('z-index', 10000); // kludgy; brings it on top
          return this;
        },        
      });
      
      formSchema['username'] = {
          title: '',
          key: 'username',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          placeholder: 'Find a User...',
          options: appModel.getUserOptions(),
          template: self.choiceFieldTemplate 
        };
      
      formSchema['screen_facility_id'] = {
          title: '',
          key: 'screen_facility_id',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          placeholder: 'Find a Screen...',
          options: appModel.getScreenOptions(),
          template: self.choiceFieldTemplate
        };
      
      formSchema['library_short_name'] = {
          title: '',
          key: 'library_short_name',
          type: EditView.ChosenSelect,
          editorClass: 'chosen-select',
          placeholder: 'Find a Library...',
          options: appModel.getLibraryOptions(),
          template: self.choiceFieldTemplate
        };
      
      formSchema['search_target'] = {
        title: 'Search', 
        key:  'search_target', // TODO: "key" not needed>?
        type: EditView.ChosenSelect,
        placeholder: 'Search by...',
        options: [{ val: 'reagent', label: 'Wells'},
                  { val: 'librarycopyplate', label: 'Plate copies'},
                  { val: 'librarycopy', label: 'Copies'},
                  { val: 'copywell', label: 'Copy Wells'},
                  ],
        template: self.choiceFieldTemplate,
        editorClass: 'chosen-select'
      };
      
      formSchema['search_value'] = {
        title: 'for', 
        help: 'enter a comma separated list',
        key:  'search_value', // TODO: "key" not needed>?
        type: TextArea2,
        template: self.fieldTemplate,
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

      this.listenTo(appModel, 'change:uriStack' , this.uriStackChange );
      this.on("change", function(){
        var username = self.getValue('username');
        var screen_facility_id = self.getValue('screen_facility_id');
        var library_short_name = self.getValue('library_short_name');
        if(username){
          var route = 'screensaveruser/' + username;
          self.setValue('username',null);
          self.$el.find('[key="username"]')
              .find('.chosen-select').trigger("chosen:updated");
          appModel.router.navigate(route, {trigger: true});
        }
        else if(screen_facility_id){
          var route = 'screen/' + screen_facility_id;
          self.setValue('screen_facility_id',null);
          self.$el.find('[key="screen_facility_id"]')
              .find('.chosen-select').trigger("chosen:updated");
          appModel.router.navigate(route, {trigger: true});
        }
        else if(library_short_name){
          var route = 'library/' + library_short_name;
          self.setValue('library_short_name',null);
          self.$el.find('[key="library_short_name"]')
              .find('.chosen-select').trigger("chosen:updated");
          appModel.router.navigate(route, {trigger: true});
        }
      });

      SearchView.__super__.initialize.apply(this, args);
      console.log('initialize SearchView');
    },

    /**
     * Backbone.Model change event handler
     * @param options.source = the event source triggering view
     */
    uriStackChange: function(model, val, options) {
      var self=this;

      var uriStack = _.clone(appModel.get('uriStack'));
      console.log('uristack',uriStack);
      var uiResourceId = uriStack.shift();
      self.setValue('search_target', uiResourceId);
      
      if(uriStack.length > 1 && uriStack[0] == 'search'){
        var stack_change_type = uriStack.shift();
        console.log('create a search_box with search data');
        var searchID = uriStack.shift();
        
        if(searchID == this.searchID){
          console.log('self generated uristack change');
          return;
        }
        if(searchID && _.isNaN(parseInt(uriStack[1]))){
          console.log('search ID is not a valid number: ' + searchID);
          return;
        }
        var search_data = appModel.getSearch(searchID);
        if(_.isUndefined(search_data) || _.isEmpty(search_data)){
          var msg = 'Search box requires a "search_data:'+searchID +'" in current browser state';
          console.log(msg);
          return;
        }else{
          
          console.log('setup search box for search data', search_data);
          var search_string = '';
          
          search_string = _.reduce(search_data, function(memo, orclause){
            if(memo != '') memo +='\n';
            var fieldKeys = [
                'library_short_name__icontains', 'plate_number__in',
                'plate_number__range', 'copy_name__in', 'well_name__in', 'well_id__in'];
            var temp = '';
            _.each(fieldKeys, function(key){
              if(_.has(orclause, key)){
                var val = orclause[key];
                if(key == 'plate_number__range'){
                  val = val.split(',').join('-');
                }
                if(temp !==''){
                  temp += ' '
                }
                temp +=  val;
              }
            });
            
            return memo + temp;
          }, '');
          console.log('setValue: ' + search_string)
          self.setValue('search_value', search_string);
          this.searchID = searchID;
        }
      }else{
        console.log('Search box requires the "search/[searchID]" param');
        self.setValue('search_value', '');
        self.searchID = null;
        return;
      }
    },    
    
    /** 
     * Override the Backbone Layoutmanager template rendering to use 
     * Backbone Forms
     */
    renderTemplate: function() {
      return Backbone.Form.prototype.render.apply(this);
    },
    
    afterRender: function() {
      var self=this;
      this.$el.append([
          '<div class="col-xs-3">',
          '<button type="submit" class="btn btn-default btn-xs" style="width: 3em; " >ok</input>',
          '</div>',
          '<div class="col-xs-3">',
          '<a class="backgrid-filter clear" data-backgrid-action="clear"',
          ' href="#">&times;</a></div>',
          '</div>',
          ].join(''));

      this.$el.find('[ type="submit" ]').click(function(e){
        e.preventDefault();
        self._submit();
      });
      
      console.log('setup single selects using chosen...');
      // See http://harvesthq.github.io/chosen/
      this.$el.find('.chosen-select').chosen({
        disable_search_threshold: 3,
        width: '100%',
        allow_single_deselect: true,
        search_contains: true
        });
    },
    
    /**
     * parse:
     * - newline for "or"
     * - space for "and"
     * - comma for "in"
     * creates search_array (example):
     * [ [ [1124,1123],[A01,A02,D12,G14] ],
     *   [ [1200-1300],[A01,A03] ]
     * ]
     * will translated and sent to server as:
     * [ 
     *   { plate_number__range=[1222,1332], well_id__in=[b02,c02,d03,N12] },
     *   { plate_number__range=[1222,1332], well_id__in=[b02,c02,d03,N12] }
     * ]
     */
    parse_search_val: function(val){
      var search_array = []
      var or_list = val.split('\n');
      _.each(or_list, function(clause){
        clause = clause.trim();
        if(clause=='') return;
        // split on space, but not space followed by a comma:
        // (reverse) the string to do (negative) look ahead to exclude ',\s' 
        // do this because JavaScript doesn't support negative lookbehind
        clause = clause.split('').reverse().join('');
        var and_list = clause.split(/\s+(?!,)/);
        var and_array = [];
        search_array.push(and_array);
        _.each(and_list, function(reversedterms){
          // now reverse each term to put it back in proper non-reversed form
          var terms = reversedterms.split('').reverse().join('');
          and_array.push(terms.split(/[,\s]+/));
        });
      });
      return search_array;
    },
    
    _submit: function(){
      var self  = this;
      // validate:true: tells bbf to run model.validate(), in addition to field[].validate()
      var errors = self.commit({ validate: true }); 
      if(!_.isEmpty(errors)){
        console.log('form errors, abort submit: ' + JSON.stringify(errors));
        this.$el.find('#search_value').addClass(self.errorClass);
        return;
      }else{
        this.$el.find('#search_value').removeClass(self.errorClass);
      }
      
      var search_target = self.getValue('search_target');
      if(!search_target ||
          _.isUndefined(search_target)||_.isEmpty(search_target.trim())){
        console.log('cancelling submit with nothing entered');
        return;
      }
      search_target = search_target.trim();
      var search_array = self.parse_search_val(self.getValue('search_value'));
      console.log('search_target', search_target, 'search_array', search_array);
      
      var matching_hash = {};
      
      if(search_target == 'reagent'){
      
        matching_hash = {
          matching_order: ['plate_number','well_id','well_name'],
          plate_number: {
            title: 'Plate number',
            pattern: /^\d{1,5}$/,
            help: 'a 1-5 digit number'
          },
          well_id: {
            title: 'Well ID',
            pattern: /(\d{1,5})\:([a-zA-Z]{1,2}\d{1,2})/,
            help: '[plate_number]:[well_name]'
          },
          well_name: {
            title: 'Well name',
            pattern: /([a-zA-Z]{1,2})(\d{1,2})/,
            help: '[row][col] in the form A01,A12,B01 etc.',
            parser: function(val_list){
              return _.map(val_list, function(term){
                var rowcol = matching_hash.well_name.pattern.exec(term);
                // check the col, 2 digits
                if(rowcol[2].length == 1) term = rowcol[1] + '0'+rowcol[2];
                return term.toUpperCase();
              });
            }
          }
        };
      }else if(search_target == 'librarycopyplate'){
        matching_hash = {
          matching_order: ['plate_number','copy_name'],
          copy_name: {
            title: 'Copy name',
            pattern: /\b([a-zA-Z0-9]*[a-zA-Z][a-zA-Z0-9]*)\b/,
            help: 'letters and numbers with no spaces'            
          },
          plate_number: {
            title: 'Plate number',
            pattern: /^\d{1,5}$/,
            help: 'a 1-5 digit number'
          },
        };
      }else if(search_target == 'librarycopy'){
        matching_hash = {
          matching_order: ['name','library_short_name'],
          library_short_name: {
            title: 'Library short name',
            pattern: /\w+/,
            help: 'single word pattern to match',
            type: 'icontains'
          },
          name: {
            title: 'Copy name',
            pattern: /\b([a-zA-Z0-9]*[a-zA-Z][a-zA-Z0-9]*)\b/,
            help: 'letters and numbers with no spaces'            
          },
        };
      }else if(search_target == 'copywell'){
        matching_hash = {
          matching_order: ['well_name','plate_number','copy_name'],
          copy_name: {
            title: 'Copy name',
            pattern: /\b([a-zA-Z0-9]*[a-zA-Z][a-zA-Z0-9]*)\b/,
            help: 'letters and numbers with no spaces'            
          },
          plate_number: {
            title: 'Plate number',
            pattern: /^\d{1,5}$/,
            help: 'a 1-5 digit number'
          },
          well_name: {
            title: 'Well name',
            pattern: /([a-zA-Z]{1,2})(\d{1,2})/,
            help: '[row][col] in the form A01,A12,B01 etc.',
            parser: function(val_list){
              return _.map(val_list, function(term){
                var rowcol = matching_hash.well_name.pattern.exec(term);
                // check the col, 2 digits
                if(rowcol[2].length == 1) term = rowcol[1] + '0'+rowcol[2];
                return term.toUpperCase();
              });
            }
          }          
          
        };
      }else{
        var msg = 'Unknown search target: ' + search_target;
        console.log(msg);
        appModel.error(msg);
        return;
      }
      
      var errors = [];
      var or_list = [];
      _.each(search_array, function(and_array){
        var and_hash = {};
        or_list.push(and_hash);
        _.each(and_array, function(terms){
          // find out what kind of terms they are
          var clause = ''
          var type = 'in';
          var entity = '';
          
          var testing_val = terms[0];
          if(testing_val.indexOf('-')>0){
            type='range';
            var testing_val = testing_val.split('-')[0];
          }
          entity = _.find(matching_hash.matching_order,function(key){
            if(matching_hash[key].pattern.exec(testing_val)) return true;
          });
          
          if(_.isEmpty(entity)){
            errors.push(testing_val);
          }else{
            
            if(matching_hash[entity].parser){
              terms = matching_hash[entity].parser(terms);
            }
            
            if(matching_hash[entity].type) type = matching_hash[entity].type;
            
            var search_entry = entity + '__' + type;
//            // Hack because pattern for library_short_name is equivalent to the copy_name pattern
//            if( search_entry == 'library_short_name__icontains'){
//              if(!_.has(and_hash, 'copy_name__in')){
//                search_entry = 'copy_name__in';
//              }
//            }
            
            if(type == 'range'){
              and_hash[search_entry] = terms[0].split('-').join(',');
            }else{
              and_hash[search_entry] = terms.join(',');
            }
          }
        });
      });
      
      var data_el = this.$el.find('#data-error');
      if(!_.isEmpty(errors)){
          data_el.html('<div class="help-block"><h5>Unparseable search values:</h5>"' + 
              errors.join('" ,"') + '"</br>' +
              "<h5>Allowed patterns:</h5>" + 
              _.map(matching_hash, function(val,key){
                return '<dt>' + val.title + '</dt><dd>' + val.help + '</dd>'; 
              }).join('') + "</dl></div>");
          data_el.show();
          return;
      }else{
        data_el.hide();
      }
      
      console.log('or_list', or_list);
      var has_val = !_.isEmpty(or_list);
      if(has_val){
        has_val = _.find(or_list,function(and_hash){
          return !_.isEmpty(and_hash) && 
              !_.isUndefined(_.find(_.values(and_hash), function(val){
                return !_.isUndefined(val) && !_.isEmpty(val);
              }));
        });
      }
      if(!has_val){
        console.log('submit with nothing entered that is recognized', 
            self.getValue('search_value'));
        return;
      }
      if(or_list && or_list.length > 1){
        // then this a composite/or search
        var searchID = ( new Date() ).getTime();
        
        // must change the route, and create a post
        
        var resource = appModel.getResource(search_target);
        var newStack = [resource.key,'search',searchID];
        // TODO: key the search data using the searchID: 
        // this allows for bfowser "back" in the session
        // will also need to listen to URIStack changes and grab the data
        // from the search ID
        var _serialized = JSON.stringify(or_list);
        appModel.setSearch(searchID,or_list);
        this.searchID = searchID;
        appModel.set('routing_options', {replace: false});  
        var _route = search_target + '/search/'+ searchID;
        appModel.router.navigate(_route, {trigger:true});
        // NOTE: easier to control the router history using navigate: 
        // when using uristack, there is the problem of who set appModel.routing_options last:
        // a race condition is set up between list2.js and search_box.js
        //        appModel.set({'uriStack': newStack});     
      }else{
        // this is a simple "AND" search -which be translated to the list for handling
        var _route = search_target + '/search/';
        var routes = [];
        _.each(_.pairs(or_list[0]),function(pair){
          var val = pair[1];
          if(_.isArray(val)){
            val = val.join(',');
          }
          routes.push(pair[0] + '=' + val);
        });
        _route += routes.join(';')
        appModel.set('routing_options', {replace: false});  
        appModel.router.navigate(_route, {trigger:true});
      }
    }
    
  });

  return SearchView;
});