#!/bin/sh

# run uglify-js, a dependency that may be installed with
#  > npm install uglify-js
# or
#  > npm install -g uglify-js

uglifyjs jmat.js -m -c --lint -v --comments "/(Copyright)|(Licensing)/" --source-map jmat.map >jmat.min.js