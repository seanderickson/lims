{   
    "admin": {
        "title": "Admin",
        "route": "",
        "view": "HomeView",
        "content_header": "Welcome",
        "description": "Perform admin activities"
    },
    "home": {
        "title": "Screensaver Reporting",
        "route": "/",
        "view": "HomeView",
        "content_header": "Welcome",
        "description": "Menu starting point"
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
    "users": {
        "header_message": "Django users",
        "title": "Users",
        "route": "list/users",
        "list_view": "ListView",
        "detailView": "UserAdminView",
        "api_resource": "user",
        "url_root": "/reports/api/v1",
        "description": "Django user"
    },
    "groups": {
        "header_message": "User Groups",
        "title": "User Groups",
        "route": "list/groups",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "usergroup",
        "url_root": "/reports/api/v1",
        "description": "User Groups"
    },
    "permissions": {
        "header_message": "Django permissions",
        "title": "Permissions",
        "route": "list/permissions",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "permission",
        "url_root": "/reports/api/v1",
        "description": "Permissions"
    },

    "screensaveruser": {
        "header_message": "All users (Screeners and Staff)",
        "title": "Screensaver Users",
        "route": "list/screensaveruser",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "screensaveruser",
        "url_root": "/db/api/v1",
        "description": "View user information"
    },

    "screeners": {
        "header_message": "Screening Users",
        "title": "Screeners",
        "route": "list/screeners",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "screensaveruser",
        "url_root": "/db/api/v1",
        "description": "View user information",
        "options": { "search": {"screeningroomuser__isnull": "False"} }
    },

    "staff": {
        "header_message": "Staff",
        "title": "Staff Users",
        "route": "list/staff",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "screensaveruser",
        "url_root": "/db/api/v1",
        "description": "View user information",
        "options": { "search": {"administratoruser__isnull": "False"} }
    },
    "screen": {
        "header_message": "All screens (Small Molecule and RNAi)",
        "title": "Screens",
        "route": "list/screen",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "screen",
        "url_root": "/db/api/v1",
        "description": "View screen information",
        "options": { "order_by": { "facility_id":"-"} }
    },
    "small_molecule_screens": {
        "header_message": "Small Molecule Screens",
        "title": "Small Molecule",
        "route": "list/small_molecule_screens",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "screen",
        "url_root": "/db/api/v1",
        "description": "View small molecule screen information",
        "options": { "search": { "screen_type": "small_molecule"} }
    },
    "rnai_screens": {
        "header_message": "All screens (Small Molecule and RNAi)",
        "title": "RNAi",
        "route": "list/rnai_screens",
        "list_view": "ListView",
        "detailView": "DetailView",
        "api_resource": "screen",
        "url_root": "/db/api/v1",
        "description": "View rnai screen information",
        "options": { "search": { "screen_type": "rnai"} }
    },
    "library": {
        "header_message": "All libraries (Small Molecule and RNAi)",
        "title": "Libraries",
        "route": "list/library",
        "list_view": "ListView",
        "detailView": "LibraryView",
        "api_resource": "library",
        "url_root": "/db/api/v1",
        "description": "View library information"
    },
    "smallmoleculelibrary": {
        "header_message": "Small Molecule Libraries",
        "title": "Small Molecule",
        "route": "list/smallmoleculelibrary",
        "list_view": "ListView",
        "detailView": "LibraryView",
        "api_resource": "library",
        "url_root": "/db/api/v1",
        "description": "View Small Molecule Library information",
        "options": { "search": { "screen_type": "small_molecule"} }
    },
    "rnalibrary": {
        "header_message": "RNAi Libraries",
        "title": "RNAi",
        "route": "list/rnalibrary",
        "list_view": "ListView",
        "detailView": "LibraryView",
        "api_resource": "library",
        "url_root": "/db/api/v1",
        "description": "View RNAi library information",
        "options": { "search": { "screen_type": "rnai"} }
    }
}