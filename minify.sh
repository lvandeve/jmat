#!/bin/sh

# This not only minifies but also lints the code.

# To run this on Archlinux. (But this info may still be helpful other OSes too)
# *) For Closure Compiler:
# compiler.jar is the closure compiler, available from https://github.com/google/closure-compiler, unzip its latest download and rename to compiler.jar
# Just have compiler.jar in this directory, no need to install closure-compiler from AUR in archlinux. But install java-runtime if needed.
# *) For uglifyjs: [no longer used here]
# pacman -S uglify-js
# *) For jslint: to install it: [not used here]
# run in this directory: npm install jslint. Then you can execute it with node node_modules/jslint/bin/jslint.js
# *) For jshint: to install it:
# run in this directory: npm install jshint. Then you can execute it with node node_modules/jshint/bin/jshint

# Run jshint
echo 'running jshint'
# This is horrible! jshint has no command line args to disable warnings! You must create a local file for that... Done here by overwriting .jshintrc:
# The disabled warnings (with reasoning why) are:
# W004: "... is already defined":
# --> I do use var in different scopes and don't want to not do it.
# W069: "... is better written in dot notation"
# --> sometimes string notation used for readability, don't want to have warning for that
# W083: "Functions declared within loops referencing an outer scoped variable may lead to confusing semantics.":
# --> I do use that pattern (and am aware that for callbacks stored for later you don't want this), and works and looks fine where used.
# W018: "Confusing use of '!'"
# --> warns about legit use of !! to convert something that can be undefined or boolean into boolean. So disabled.
# echo '{ "-W004":false, "-W069":false }' > .jshintrc
echo '{ "-W004":false, "-W069":false, "-W083":false, "-W018":false }' > .jshintrc

node node_modules/jshint/bin/jshint --verbose jmat_real.js jmat_complex.js jmat_matrix.js jmat_quaternion.js jmat_special.js jmat_bigint.js jmat.js jmat_test.js jmat_plot.js



# Unfortunately, since I'm not using JSDoc in most places, we cannot use the closure compiler linter and have to disable a lot of warnings.
# The warnings that are turned off are the following:
# global this: closure and jsdoc do not support our dual constructor/factory methods
# usage of new: closure and jsdoc do not support our dual constructor/factory methods
# amount of arguments: because optional arguments are not annotated due to not using jsdoc
# duplicate var: vars work per function instead of scope, but for safe copypastable code we reuse var sometimes anyway as it it harmless

echo 'running closure compiler'
java -jar compiler.jar jmat_real.js jmat_complex.js jmat_matrix.js jmat_quaternion.js jmat_special.js jmat_bigint.js jmat.js --create_source_map jmat.map -W VERBOSE --jscomp_off={globalThis,checkTypes,duplicate} >jmat.min.js

# uglifyjs disabled! We now have jshint instead

# Run uglifyjs as well, because it gives warnings about unused local variables that the closure compiler doesn't give.
#echo 'running uglifyjs'
# Filtering out "Dropping duplicated definition" because there are a lot of those, and personally I prefer to rewrite "var" in different scopes, and that's ok.
# Filtering out "unreachable" because uglifyjs prints warnings about it even though the code is reachable
# Filtering out "Condition always true" because sometimes variables are inited to true and tested with conditions to enable/disable behavior
# What we still find are unused variables...
#uglifyjs jmat_real.js jmat_complex.js jmat_matrix.js jmat_quaternion.js jmat_special.js jmat_bigint.js jmat.js -m -c --comments --warn -o /dev/null 2>&1 | grep -v "Dropping duplicated definition" | grep -v "unreachable" | grep -v "Condition always true"

#previous version of uglifyjs had different syntax:
#uglifyjs jmat_real.js jmat_complex.js jmat_matrix.js jmat_quaternion.js jmat_special.js jmat_bigint.js jmat.js -m -c --lint -v --comments --source-map /dev/null > /dev/null

