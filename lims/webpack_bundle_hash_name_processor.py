from django.conf import settings

#
# This Context Processor exposes the hash value used to make the webpack bundle
# name as the setting:
# "WEBPACK_BUNDLE_NAME=bundle.[hash_value].js" 
# 
# - webpack bundles will have the naming "bundle.[hash_value].js"
# - this setting is used in the index_require.html to create the src tag:
# <script src="{{ STATIC_URL }}{{ WEBPACK_BUNDLE_NAME }}" ></script>
# 
# Background: 
# see: reports/static/webpack.config:
# - Uses the "assets-webpack-plugin" plugin to generate a webpack bundle with
# a hash value in the bundle name, as in "bundle.[hash_value].js"
# (for cache busting).
# See: https://github.com/kossnocorp/assets-webpack-plugin
# 
# - The configured assetsPluginInstance will generate the settings file:
#
# [project_dir]/reports/static/bundle_name.js, e.g.:
# {"main":{"js":"http://localhost:8000/_static/bundle.c026caedeb210bd3ec0a.js"}}
# 
# 
import os
import json
import logging
import re

logger = logging.getLogger(__name__)

BUNDLE_NAME_PATTERN = re.compile(r'.*\/(bundle\..*js)')

def bundle_context_processor(request):
    
    bundle_filename = 'bundle_name.json'
    filename = os.path.join(
        settings.PROJECT_ROOT, '..','reports','static',bundle_filename)
    try:
        with open(filename) as f:
        
            assets = json.loads(f.read())
            
            logger.info('bundle assets: %r', assets )
            if 'main' not in assets:
                raise Exception(
                    'bundle filename: %r, does not contain a "main" entry: %r', 
                    bundle_filename, assets)
            if 'js' not in assets['main']:
                raise Exception(
                    'bundle filename: %r, does not contain a "[main][js]" entry: %r', 
                    bundle_filename, assets)
            
            match = BUNDLE_NAME_PATTERN.match(assets['main']['js'])
            
            if not match:
                raise Exception('bundle name does not match the pattern: %r', assets)
            
            return { 'WEBPACK_BUNDLE_NAME': match.group(1) }    
    except:
        logger.exception('filename not found: %s' % bundle_filename )
        return {}
