from django import template

register = template.Library()
# following: http://gnuvince.wordpress.com/2007/09/14/a-django-template-tag-for-the-current-active-page/
@register.simple_tag
def activeTag(request, pattern):
    """ 
    Return the string "active" if the request.path matches the pattern.  For use with the <select> HTML element.
    """
    import re
    if re.search(pattern, request.path):
        return 'active'
    return ''
