{   
    "admin": {
        "title": "Admin",
        "route": "",
        "view": "HomeView",
        "content_header": "Welcome",
        "description": "Perform admin activities"
    },
    "home": {
        "title": "Screensaver LIMS",
        "route": "/",
        "view": "HomeView",
        "content_header": "Welcome",
        "description": "Menu starting point"
    },
    
    "about": {
      "title": "ICCB-L Screensaver LIMS",
      "route": "about",
      "view": "AboutView",
      "content_header": "ICCB-L Screensaver LIMS",
      "description": "About page"
    },

    "metahash": {
        "header_message": "Define fields for display on detail and list views",
        "title": "Field Information",
        "route": "list/metahash",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "metahash",
        "url_root": "/reports/api/v1",
        "options": { 
          "search": {"scope__exact": "fields.metahash"}},
        "description": "Control field settings"
    },

    "resource": {
        "header_message": "Define resources for display in the reporting interface",
        "title": "Resource Information",
        "route": "list/resource",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "resource",
        "url_root": "/reports/api/v1",
        "description": "Control resource information"
    },

    "vocabularies": {
        "header_message": "Define controlled vocabularies",
        "title": "Application Vocabularies",
        "route": "list/vocabularies",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "vocabularies",
        "url_root": "/reports/api/v1",
        "description": "Enter controlled vocabularies"
    },

    "apilog": {
        "header_message": "View change logs",
        "title": "Change logs",
        "route": "list/apilog",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "apilog",
        "url_root": "/reports/api/v1",
        "description": "Change logs"
    },

    "user": {
      "header_message": "Reports users",
      "title": "Reports Users",
      "list_view": "ListView",
      "detailView": "UserAdminView",
      "api_resource": "user",
      "url_root": "/reports/api/v1",
      "description": "Reports user"
    },
    "usergroup": {
        "header_message": "User Groups",
        "title": "User Groups",
        "route": "list/groups",
        "list_view": "ListView",
        "detailView": "UserGroupAdminView",
        "api_resource": "usergroup",
        "url_root": "/reports/api/v1",
        "description": "User Groups"
    },
    "permission": {
        "header_message": "Django permissions",
        "title": "Permissions",
        "route": "list/permissions",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "permission",
        "url_root": "/reports/api/v1",
        "description": "Permissions"
    }
}