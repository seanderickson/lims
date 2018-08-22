const merge = require('webpack-merge');
const common = require('./webpack.common.js');

module.exports = merge.smart(common, {
  mode: 'development',
  devtool: 'inline-source-map',
  output: {
    // with css?sourceMap; the public path must be a full URL for fonts to load
    // see: https://github.com/webpack/css-loader/issues/29
    publicPath: 'http://localhost:8000/_static/'
  },
  devtool: 'source-map',
  module: {
    rules: [
      {
        test: /\.css$/,
        use: [
          { loader: "style-loader",
            options: {
              sourceMap: true
            }
          },
          { loader: "css-loader",
            options: {
              sourceMap: true
            }
          }
        ]
      }
    ]
  }
});

