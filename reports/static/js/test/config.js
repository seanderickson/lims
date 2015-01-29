mocha.setup('bdd');

// Mocha run helper, used for grunt-mocha with amd, see
// see: https://gist.github.com/kmiyashiro/2655876
var runMocha = function() {
  //  chai.use(chaiJquery);
  console.log('mocha run1');
  mocha.run();
  console.log('mocha run2');
};

