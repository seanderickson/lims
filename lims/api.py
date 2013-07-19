
from tastypie.resources import ModelResource
from tastypie.serializers import Serializer

import json
import logging
        
logger = logging.getLogger(__name__)


class BackboneSerializer(Serializer):
    
    def from_json(self, content):
        """
        Given some JSON data, returns a Python dictionary of the decoded data.
        Override to quote attributes - the backbone client doesn't want to do this.
        """
        logger.info(str(("loading content:", content)))
        content = content.replace(r'(\w+):', r'"\1" :')
        logger.info(str(("loading content:", content)))
        return json.loads(content)

class PostgresSortingResource(ModelResource):

    def apply_sorting(self, obj_list, options):
        """
        Create a non-too-pretty workaround for the postgresql null sorting issue - nulls sort higher than values, 
        which is not desired.  We want nulls to sort lower than values.
        """ 
        
        obj_list = super(PostgresSortingResource, self).apply_sorting(obj_list, options)
        
        extra_select = {}
        extra_ordering = []
        for field in obj_list.query.order_by:
            is_null_dir = '-'  # default nulls first for ascending
            if field.startswith('-'):
                is_null_dir = ''
                field = field[1:]
            extra_select[field+"_null"]=field + ' is null'
            extra_ordering.append(is_null_dir + field+"_null")
        logger.info(str(('extra_select', extra_select)))
        obj_list = obj_list.extra(extra_select)

        # Note that the following doesn't work, something in the framework deletes the extra 
        # order_by clause when apply_sorting, or, if this is run last, it deletes the sorting applied in apply_sorting...
        #        obj_list = obj_list.extra(order_by=['-comments_null'])

        # note: this doesn't work because the "is null" field order by clauses
        # must be prepended so that they occur before their intended fields
        #        obj_list.query.add_ordering('comments_null')
        
        temp = obj_list.query.order_by;
        obj_list.query.clear_ordering()
        for xfield in extra_ordering:
            temp.insert(0,xfield)
        obj_list.query.add_ordering(*temp)
        
        return obj_list

