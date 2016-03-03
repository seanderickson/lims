define([
        'jquery',
        'underscore',
        'backbone',
        'sinon', 'fakeServer', 'chai',
        'iccbl_backgrid',
        'models/app_state',
        'views/generic_detail_stickit',
        'views/generic_edit',
        'text!test/models/detail_test_model.json',
        'text!test/models/detail_test_resource.json', 
        'text!test/models/detail_test_vocabularies.json',    
        'text!test/models/detail_test_users.json'    
], function($, _, Backbone, 
            sinon, fakeServer, chai, 
            Iccbl, appModel, DetailView, EditView,
            test_model_raw,resource_raw,test_vocabularies_raw, test_users_raw) {
  var expect = chai.expect;
  var assert = chai.assert;
  console.log('start testDetail');
  
  window.user = 'testuser';  
  
  describe("load the Detail view test", function(){
    before(function () {
      this.$fixture = $("<div id='detail-fixture'></div>");
      this.resource = JSON.parse(resource_raw);
      this.resource.schema = _.extend(this.resource.schema, new Iccbl.SchemaClass());
      _.extend(appModel.get("vocabularies"), JSON.parse(test_vocabularies_raw));
      appModel.set("users", JSON.parse(test_users_raw));
    });
    beforeEach(function () {
      this.$fixture.empty().appendTo($("#fixtures"));
      this.model = new Backbone.Model(JSON.parse(test_model_raw));
      this.model.resource = this.resource;
      this.view = new DetailView({ 
        useRAF: false,
        model: this.model
      });
      this.view.render();
      this.$fixture.append(this.view.el);
    });

    afterEach(function () {
      // Destroying the model also destroys the view.
      this.view.model.destroy();
    });

    it( "should load the Detail view with the model", function() {
      var self = this;
      
      _.each(this.model.keys(),function(key){

        var modelVal = self.model.get(key);
        var finalModelVal = modelVal
        var fi = self.resource.schema.fields[key];
        var data_type = fi['data_type'];
        var display_type = fi['display_type'];
        var vocabulary_scope_ref = fi['vocabulary_scope_ref'];
        var displayedVal,vocabulary,cell_options;
        var input = self.view.$('#'+key);
        cell_options = fi.display_options;
        
        // TODO: use the app_state method to get the resource/parse the display_options
        if (!_.isEmpty(fi.display_options)){
          cell_options = fi['display_options'];
          cell_options = cell_options.replace(/'/g,'"');
          try{
            cell_options = JSON.parse(cell_options);
          }catch(e){
            assert(false,'warn: display_options is not JSON parseable, column: ',
                key,', options: ',cell_options);
          }
        }
        
        if (!_.isEmpty(vocabulary_scope_ref)){
          vocabulary = Iccbl.appModel.getVocabulary(vocabulary_scope_ref);
        }

        console.log('detail testing: ', key, ', types: ',
          [data_type,display_type],', ',input.html(), ', ', modelVal );

        if(_.has(fi,'visibility') && _.contains(fi['visibility'],'detail')){
        
          expect(input.length, 'input for selector not found: ' + key);
          
          displayedVal = input.html();
          
          if (! modelVal ){
            expect(displayedVal).to.be.equal('-');
          }else{
            if(display_type =='link'){
              if(data_type=='list'){
                input = input.find("a");
                expect(input.length).to.be.equal(modelVal.length);
                // TODO: link list testing
              }else{
                if(vocabulary){
                  finalModelVal = vocabulary.getTitle(modelVal);
                }
                input = input.find("a");
                expect(input.length).to.be.equal(1);
                expect(input.html()).to.be.equal(finalModelVal);

                expect(input.attr('href')).to.be.equal(
                  Iccbl.formatString(
                    cell_options.hrefTemplate,
                    self.model,
                    modelVal));
              }
            }else if (data_type == 'date') {
              modelVal = Iccbl.getUTCDateString(new Date(modelVal));
              expect(modelVal).to.be.equal(displayedVal);
            }else if(data_type=='list'){
              if(vocabulary){
                modelVal = Iccbl.sortOnOrdinal(modelVal,vocabulary);
                modelVal = _.map(modelVal,vocabulary.getTitle);
                modelVal = modelVal.join(', ');
                console.log('sorted modelVal: ' + modelVal);
                expect(modelVal).to.be.equal(displayedVal,
                  'vocabulary_scope_ref: ' + vocabulary_scope_ref + ', ' + modelVal + ', ' + displayedVal);
              }else{
               expect(modelVal.join(', ')).to.be.deep.equal(displayedVal);
              }
            }else if(data_type=='decimal' || data_type=='float'){
              console.log('testing decimal');
              if (display_type == 'siunit'){
                console.log('testing siunit');
                expect(displayedVal).to.be.equal(
                  (new Iccbl.SIUnitsFormatter(cell_options).fromRaw(modelVal)));
              }else{
                expect(displayedVal).to.be.equal(
                  (new Iccbl.DecimalFormatter(cell_options).fromRaw(modelVal)));
              }
            }else{
              console.log('default data type:', data_type);
              if(vocabulary){
                modelVal = vocabulary.getTitle(modelVal);
                expect(modelVal).to.be.equal(displayedVal,
                  'vocabulary_scope_ref: ' + vocabulary_scope_ref);
              }else{
                expect('' + modelVal).to.be.equal(displayedVal)
              }
            }
          }// if (modelVal)
        }// if(fi)
        
      });
    });
    
  });
  
  describe("load Edit view test:", function(){
    before(function () {
      var self = this;
      console.log('edit, before');
      this.$fixture = $("<div id='detail-fixture'></div>");
      this.resource = JSON.parse(resource_raw);
      this.resource.schema = _.extend(this.resource.schema, new Iccbl.SchemaClass());
      _.extend(appModel.get("vocabularies"), JSON.parse(test_vocabularies_raw));
      try{
        if(appModel.has('users')){
          appModel.set("users", new Backbone.Collection(
            _.union(appModel.get("users"), JSON.parse(test_users_raw))));
        }else{
          appModel.set("users", new Backbone.Collection(
            JSON.parse(test_users_raw)));
        }
        appModel.getUserOptions(function(options){
          self.resource.schema.fields['field13']['choices'] = options;
        });
      }catch(e){
        console.log('error reading users,',e);
        throw e;
      }
      console.log('edit, before done');
    });

    beforeEach(function () {
      this.$fixture.empty().appendTo($("#fixtures"));
      this.model = new Backbone.Model(JSON.parse(test_model_raw));
      this.model.resource = this.resource;
//      var onEditCallBack = function(displayFunction){
//        appModel.getUserOptions(function(options){
//          self.model.resource.schema.fields['users']['choices'] = options;
//          displayFunction();
//        });
//      };
    
      this.view = new EditView({ 

        model: this.model
      });
      this.view.render();
      this.$fixture.append(this.view.el);
    });

    afterEach(function () {
      // Destroying the model also destroys the view.
      this.view.model.destroy();
    });
  
    it( "should load the Edit view with the model", function() {
      var self = this;
      
      _.each(this.model.keys(),function(key){
      
        var searchEl = '[name="{field}"]'.replace(/{\w+}/,key);
        var input = self.view.$(searchEl);
        var fi = self.resource.schema.fields[key];
        var data_type = fi['data_type'];
        var edit_type = fi['edit_type'];
        var val,temp;
        var modelVal = self.model.get(key);

        if(_.has(fi,'visibility') && _.contains(fi['visibility'],'edit')){
          val = input.val();
          
          console.log('Edit testing key:', key,',data_type:', data_type,
            ', edit_type: ', edit_type,
            ', modelVal: ', modelVal, 'el val:', val); 
          expect(input[0]).to.exist;
  
          if(data_type=='date'){
            
            val = input.datepicker('getDate');
            val = Iccbl.getISODateString(val);
            expect(modelVal).to.be.equal(val);
            
          }else if(data_type=='boolean'){
            
            expect(modelVal).to.be.a(
              'boolean','model val is not an boolean: ' + typeof(modelVal))
            val = input.is(':checked');
            expect(modelVal).to.be.equal(val,'boolean checkbox is not true');

          }else if(edit_type=='multiselect'){
            
            expect(modelVal).to.be.a(
              'array','model val is not an array: ' + typeof(modelVal))
            val = [];
            input.find('input:checked').each(function(){
              val.push($(this).val());
            });
            console.log('val', '' + val, ', ', '' + modelVal );
            expect(modelVal.sort()).to.be.deep.equal(val.sort(),'multiselect error: ' + val);
            
          }else if(edit_type=='multiselect2'){
            
            expect(modelVal).to.be.a(
              'array','model val is not an array: ' + typeof(modelVal));
            console.log('multiselect2: val', '' + val, ', ', '' + modelVal );
            expect(modelVal.sort()).to.be.deep.equal(val.sort(),'multiselect2 error, val: ' + val);
            
          }else if (fi['display_type']=='siunit'){

            val = input.find('input').val();
            val = parseFloat(val);
            temp = input.find('select').val();
            temp = parseFloat(temp);
            val = val * temp;
            expect('' + modelVal).to.equal('' + val);   
          
          }else{
            
            expect('' + modelVal).to.equal(val);   
            
          }
          console.log('input: ', key, ', ', ', modelVal: ', modelVal, ', formVal: ' + val );
          
        }
      });
    });
  });  
});