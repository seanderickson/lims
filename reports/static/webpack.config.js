var path = require('path')
var webpack = require("webpack");
module.exports = {
  context: path.resolve(__dirname, 'js'),
  entry: './main',
  output: {
    filename: 'bundle.js',
    publicPath: '../_static/'
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
      layoutmanager: 'layoutmanager/backbone.layoutmanager'
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
      { test: require.resolve('jquery'), loader: 'expose?jQuery!expose?$' },
      { test: /\.(html|json)$/, loaders: ['raw'], exclude: /node_modules/ },
      { test: /\.css$/, loaders: ['style', 'css'] },
      {
        test: /\.(jpe?g|png|gif|svg)$/i,
        loaders: [
          'file?hash=sha512&digest=hex&name=[hash].[ext]',
          'image-webpack?bypassOnDebug&optimizationLevel=7&interlaced=false'
        ]
      },
      {
        test: /\.(eot|svg|ttf|woff|woff2)$/,
        loader: 'file?name=css/fonts/[name].[ext]'
      }    
    ]
  }
};