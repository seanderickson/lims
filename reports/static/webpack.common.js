var path = require('path')
var webpack = require("webpack");

// Plugin Configuration: Bundle hash name generator
//
// Use the "assets-webpack-plugin" plugin to generate a webpack bundle with
// a hash value in the bundle name, as in "bundle.[hash_value].js"
// (for cache busting).
// See: https://github.com/kossnocorp/assets-webpack-plugin
//
// - The configured assetsPluginInstance will generate the settings file:
//
// [project_dir]/reports/static/bundle_name.js, containing, e.g.:
// {"main":{"js":"http://localhost:8000/_static/bundle.c026caedeb210bd3ec0a.js"}}
//
const AssetsPlugin = require('assets-webpack-plugin');
// Implement the default processOutput for debug information
const assetsPluginInstance = new AssetsPlugin({
    filename: 'bundle_name.json',
    processOutput: function (assets) {
        console.log('assets: ' , assets);
        return JSON.stringify(assets);
    }
});

module.exports = {
  context: path.resolve(__dirname),
  entry: './js/main.js',
  output: {
    path : path.join(__dirname, './'),
    filename: './bundle.[hash].js',
  },
  resolve: {
    modules: [path.resolve(__dirname, 'js'), 'node_modules','js','css'],
    alias: {
      iccbl_backgrid: 'iccbl-backgrid',
      backbone_forms: 'backbone-forms',
      backgrid_paginator: 'backgrid-paginator',
      backbone_stickit: 'backbone.stickit',
      backgrid_filter: 'backgrid-filter',
      backgrid_select_all: 'backgrid-select-all',
      chosen: 'bootstrap-chosen/dist/chosen.jquery-1.4.2/chosen.jquery.min',
      quicksearch: 'jquery.quicksearch/dist/jquery.quicksearch.min.js',
      layoutmanager: 'layoutmanager/backbone.layoutmanager',
      jquery_bonsai: 'jquery-bonsai/jquery.bonsai.js',
      google_palette: 'google-palette/palette.js'
    },
    extensions: ['.js', '.json'] 
  },
  module: {
    rules: [
      {
        test: /\.css$/,
        use: [
          { loader: "style-loader" },
          { loader: "css-loader" },      
        ]
      },
      {
        test: require.resolve('jquery'),
        use: [
          {
            loader: 'expose-loader',
            options: 'jQuery'
          },          
          {
            loader: 'expose-loader',
            options: '$'
          }]
      },
      { 
        test: /\.(html)$/, 
        use: [{ loader: 'raw-loader' }],
        exclude: /node_modules/ 
      },
      {
        test: /\.(jpe?g|png|gif|svg)$/i, 
        use: [
          // inline base64 URLs for <=8k images, direct URLs for the rest
          { loader: 'url-loader?limit=8192' }] 
      }, 
      {
        test: /\.(eot|svg|ttf|woff|woff2)$/,
        use: [{ loader: 'file-loader?name=css/fonts/[name].[ext]' }]
      }
    ]
  },
  plugins: [
    assetsPluginInstance,
    new webpack.ProvidePlugin({
      'LayoutManager': 'layoutmanager'
    }),
    new webpack.ProvidePlugin({
      'BackboneForm': 'backbone_forms'
    }),
    new webpack.ProvidePlugin({
      'Stickit': 'backbone_stickit'
    })
  ]
};