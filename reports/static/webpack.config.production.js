const merge = require('webpack-merge');
const common = require('./webpack.common.js');

// TODO: 20180821: not tested on production yet
module.exports = merge.smart(common, {
  mode: 'production',
  devtool: 'inline-source-map',
  output: {
    publicPath: '../_static/'
  },
  devtool: 'source-map'
});

