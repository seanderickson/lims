* These tests currently require a local server to run 
(to avoid "Cross origin requests are only supported for HTTP." message).
* a simple server can be run with:
/reports/static $ python -m SimpleHTTPServer
* or the grunt task may be used:
* grunt test 
* (will spin up the phantomjs server.)