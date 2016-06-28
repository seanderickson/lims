var path = require('path')
var webpack = require("webpack");
module.exports = {
  context: path.resolve(__dirname, 'js'),
  entry: './main',
  output: {
    filename: 'bundle.js',
    publicPath: '../static/'
  },
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
      chosen: 'jquery-chosen/chosen.jquery.min',
      quicksearch: 'jquery.quicksearch/dist/jquery.quicksearch.min.js',
      layoutmanager: 'backbone.layoutmanager'
    },
    extensions: ['', '.js', '.json'] 
  },
//  plugins: [
//      new webpack.ResolverPlugin(
//          new webpack.ResolverPlugin.DirectoryDescriptionFilePlugin(".bower.json", ["main"])
//      )
//  ],
  module: {
    loaders: [
      { test: /\.(html|json)$/, loaders: ['raw'], exclude: /node_modules/ }
      // TODO: follows is experimental - to help load the jquery obj to the global namespace
      //{ test: /vendor\/.+\.(jsx|js)$/,
      //  loader: 'imports?jQuery=jquery,$=jquery,this=>window'
      //}
    ]
  }
};