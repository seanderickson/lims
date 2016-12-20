from django.conf import settings

#
# This Context Processor exposes the hash value used to make the webpack bundle
# name as the setting:
# "LIMS_BUNDLE_HASH=[hash_value]" 
# 
# - webpack bundles will have the naming "bundle.[hash_value].js"
# - this setting is used in the index_require.html to generate:
# <script src="{{ STATIC_URL }}bundle.{{ LIMS_BUNDLE_HASH }}.js?" ></script>
# 
# Background: 
# see: reports/static/webpack.config:
# - Uses the "assets-webpack-plugin" plugin to generate a webpack bundle with
# a hash value in the bundle name, as in "bundle.[hash_value].js"
# (for cache busting).
# See: https://github.com/kossnocorp/assets-webpack-plugin
# 
# - The configured assetsPluginInstance will generate the settings file:
# lims/webpack_bundle_hash_setting.py:
# LIMS_BUNDLE_HASH=[some hash value]
# 
# - This settings file is loaded in the base_settings.py file
# 
# 

def bundle_context_processor(request):
    return {
        'LIMS_BUNDLE_HASH': settings.LIMS_BUNDLE_HASH,
}    

