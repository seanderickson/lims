var path = require('path')
var webpack = require("webpack");
const AssetsPlugin = require('assets-webpack-plugin');
const assetsPluginInstance = new AssetsPlugin({
    update: true,
    path: path.join(__dirname, '..','..','lims'),
    filename: 'webpack_bundle_hash_setting.py',  
    processOutput: function (assets) {
        console.log('assets: ' , assets);
        var hash = assets.main.js.match(/.*bundle\.([^.]+)\.js/)[1]
        console.log('AssetsPlugin parsed hash: ' + hash);
        return 'LIMS_BUNDLE_HASH="' + hash + '"\n';
    }
});
module.exports = {
  context: path.resolve(__dirname, 'js'),
  entry: './main',
  output: {
    filename: 'bundle.[hash].js',
    // with css?sourceMap; the public path must be a full URL for fonts to load
    // see: https://github.com/webpack/css-loader/issues/29
    publicPath: 'http://localhost:8000/_static/'
  },
  devtool: 'inline-source-map',
  resolve: {
    modulesDirectories: ['.', 'node_modules'],
    root: [
      path.resolve('./node_modules')
    ],
    alias: {
      iccbl_backgrid: 'iccbl-backgrid',
      backbone_forms: 'backbone-forms',
      backgrid_paginator: 'backgrid-paginator',
      backbone_stickit: 'backbone.stickit',
      backgrid_filter: 'backgrid-filter',
      backgrid_select_all: 'backgrid-select-all',
      chosen: 'jquery-chosen/chosen.jquery.min',
      quicksearch: 'jquery.quicksearch/dist/jquery.quicksearch.min.js',
      layoutmanager: 'layoutmanager/backbone.layoutmanager'
    },
    extensions: ['', '.js', '.json'] 
  },
  module: {
    loaders: [
      { test: require.resolve('jquery'), loader: 'expose?jQuery!expose?$' },
      { test: /\.(html|json)$/, loaders: ['raw'], exclude: /node_modules/ },
      { test: /\.css$/, loaders: ['style', 'css?sourceMap'] },
      // image-webpack-loader has dependency issues on orchestra
      //{
      //  test: /\.(jpe?g|png|gif|svg)$/i,
      //  loaders: [
      //    'file?hash=sha512&digest=hex&name=[hash].[ext]',
      //    'image-webpack?bypassOnDebug&optimizationLevel=7&interlaced=false'
      //  ]
      //},
      { test: /\.(jpe?g|png|gif|svg)$/i, loader: 'url-loader?limit=8192' }, // inline base64 URLs for <=8k images, direct URLs for the rest
      {
        test: /\.(eot|svg|ttf|woff|woff2)$/,
        loader: 'file?name=css/fonts/[name].[ext]'
      }    
    ]
  },
  plugins: [
    assetsPluginInstance
  ]
};