#!/bin/sh

# run uglify-js, a dependency that may be installed with
#  > npm install uglify-js
# or
#  > npm install -g uglify-js

uglifyjs jmat_real.js jmat_complex.js jmat_matrix.js jmat_quaternion.js jmat_special.js jmat_bigint.js jmat.js -m -c --lint -v --comments --source-map jmat.map >jmat.min.js
