#!/bin/sh

# compiler.jar is the closure compiler, available from https://github.com/google/closure-compiler

# The warnings that are turned off are the following:
# global this: closure and jsdoc do not support our dual constructor/factory methods
# usage of new: closure and jsdoc do not support our dual constructor/factory methods
# amount of arguments: because optional arguments are not annotated due to not using jsdoc
# duplicate var: vars work per function instead of scope, but for safe copypastable code we reuse var sometimes anyway as it it harmless
# Unfortunately, closure compiler does not warn about unused local varaibles, a warning we would want...

java -jar compiler.jar jmat_real.js jmat_complex.js jmat_matrix.js jmat_quaternion.js jmat_special.js jmat_bigint.js jmat.js --create_source_map jmat.map -W VERBOSE --jscomp_off={globalThis,checkTypes,duplicate} >jmat.min.js

# Run uglifyjs as well, because it gives warnings about unused local variables that the closure compiler doesn't give.
uglifyjs jmat_real.js jmat_complex.js jmat_matrix.js jmat_quaternion.js jmat_special.js jmat_bigint.js jmat.js -m -c --lint -v --comments --source-map /dev/null > /dev/null
