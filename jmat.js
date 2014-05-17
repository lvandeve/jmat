/*
Jmat

Copyright (c) 2011-2014, Lode Vandevenne
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

--------------------------------------------------------------------------------

Jmat is a numerical library in JavaScript for complex, matrix and statistical
arithmetic.

--------------------------------------------------------------------------------

Most of the code is written based on formulas from Wikipedia, Mathworld,
NIST DLMF, various papers, and books "Handbook of Mathematical Functions" and
"Numerical Linear Algebra".

Two algorithms use code from third party open source libraries, converted to
JavaScript. Their licenses and attribution are included at the implementation.
It are KISS FFT (see JmatM.kiss_...) and linpack zsvdc (see JmatM.zsvdc_). They
are included in this source file to keep it self contained.

--------------------------------------------------------------------------------

Features
--------

-Complex numbers: add, multiply, power, trigonometrics, etc...
-Special functions in complex domain: gamma, beta, bessel, airy, theta, zeta, erf, agm, hypergeometric, ...
-Numerical linear algebra on complex matrices: matrix multiplication, determinant, SVD, eigenvalues, pseudo inverse, ...
-Statistical distributions pdf, cdf and qf: normal, student t, chi square, gamma, beta, weibull, fisher, ...
-Numerical algorithms: root finding, quadrature, differentiation, ...
-Prime number functions: primality test, euler's totient, prime factors, prime counting, gcd, ...
-Plotting (with jmat_plot.js)
-...


Names
-----

The following names and namespaces are introduced by this library:

Main API:
-Jmat: This is a wrapper around JmatR, JmatC and JmatM, and is the most public and easiest to use API to Jmat
Types:
-JmatC: Complex type. Also, complex functions and numerical algorithms, such as gamma, erfc, besselj, ...
-JmatM: Matrix type. Also, matrix functions and numerical linear algebra, e.g. determinant, svd, eigenvalues, ...
Sub API:
-JmatR: Real functions. Only a few, mostly related to primes. JmatC has the main bunch of math functions.
-JmatG: internal graphics functions for plotting (in other .js file: jmat_plot.js) 
Convenience factory functions:
-jmatc: creates new JmatC complex number: can be used instead of "new JmatC" and supports more argument types
-jmatm: creates new JmatM matrix: can be used instead of "new JmatM" and supports more argument types

The main data types are JmatC (complex numbers), JmatM (matrices), and plain JS numbers.
The functions are in JmatC (complex functions), JmatM (matrix functions), JmatR (real functions on JS numbers), and Jmat (wrapper around JmatC, JmatM and JmatR)

It has such typenames to avoid name collisions. However, to make usage simpler when possible, the following aliases are suggested:

Complex = JmatC; // complex type
Matrix = JmatM; // matrix type
newc = jmatc; // construct new complex number
newm = jmatm; // construct new matrix
// Assign some JmatR functions to Math, to get e.g. Math.gamma(5)
['gamma', 'factorial', 'lambertw', 'erf', 'erfc'].map(function(fun) { Math[fun] = JmatR[fun]; });


Usage
-----

See the Jmat.#### function definitions at the beginning of the source code to see a list of most functions.
These are however wrappers around JmatR, JmatC and JmatM functions. See further below in the code to find
their implementation. Some other less exposed functions can also be found further down in JmatR, JmatC and JmatM,
but these are less stable API.

The rest of this manual shows usage examples:

*) Create a complex number, e.g. 1+2i:
var complex = new JmatC(1, 2)
shorthand:
var complex = jmatc(1, 2)
its fields are re and im:
complex.re --> 1
complex.im --> 2

Most math functions in JmatC always require complex number objects (not JS numbers), even if they are all real. Here are the ways and shortcuts of making a real number "3":
var real = new JmatC(3, 0)
var real = jmatc(3, 0)
var real = jmatc(3)

Unfortunately, JavaScript does not have operator overloading, so simply using "+" to add two complex numbers does not work and the "add" function is needed instead.

*) Add two complex numbers:
jmatc(1, 2).add(jmatc(3, 4))
alternatives:
JmatC.add(jmatc(1, 2), jmatc(3, 4))
Jmat.add(jmatc(1, 2), jmatc(3, 4))
other basic operators: sub, mul, div, pow

*) Adding a real number to a complex number:
var complex = jmatc(1, 2)
complex.addr(3) --> 4+2i. Prototype member function of JmatC that takes real JS number
complex.add(jmatc(3)) --> 4+2i. Prototype member function of JmatC that takes JmatC object
JmatC.addr(complex, 3) --> 4+2i. JmatC.addr requires first argument of type JmatC, second as primitive JS number
JmatC.add(complex, jmatc(3)) --> 4+2i. JmatC.add requires both arguments to be exactly of type JmatC
Jmat.add(complex, 3) --> 4+2i. Jmat.add supports various input types and casts internally

Similar for sub, mul, div, pow, ...

*) Output JmatC object as more readable string:
jmatc(1, 2).add(jmatc(3, 4)).toString()
--> 4+6i

*) Complex trigonometrics
var z = jmatc(1, 2);
var s = JmatC.sin(z);
var c = JmatC.cos(z);
s.mul(s).add(c.mul(c)).toString(8) --> 1. This is cos^2(z) + sin^2(z), the '8' argument rounds to 8 digits to hide some numerical imprecision
Jmat.acos(5).toString() --> 2.2924316695611733i. Arccos of large argument gives imaginary number.


*) Gamma function of 5:
JmatC.gamma(jmatc(5))
more convenient (but slightly less efficient):
Jmat.gamma(5) --> The Jmat wrapper automatically casts to JmatC if regular JS number is given

*) Example of a function with many arguments: hypergeometric
Jmat.hypergeometric(0.5, 2, 3, 4).toString()
--> 0.1667-0.866i

*) Student T distribution (others are e.g. normal, standardnormal, gamma, weibull, chi_square, fisher, cauchy, ...)
Jmat.pdf_studentt(5, 5).toString()
--> 0.0017574383788078426
Jmat.cdf_studentt(5, 5).toString()
--> 0.9979476420099738
Jmat.qf_studentt(0.5, 5).toString()
--> 0.002183661175527183

*) Create a new 2x2 matrix [[1, 2], [3, 4]]
var matrix = new JmatM(2, 2); // num rows, nuw columns
matrix.e[0][0] = jmatc(1);
matrix.e[0][1] = jmatc(2);
matrix.e[1][0] = jmatc(3);
matrix.e[1][1] = jmatc(4);
shorthand for all the above:
var matrix = jmatm([[1, 2], [3, 4]]) // 2D array, array of rows
var matrix = jmatm(2, 2, 1, 2, 3, 4) // height, width, 4 elements
print it:
matrix.toString()
--> [[1, 2], [3, 4]]
The fields of the matrix are:
matrix.h: 2 (number of rows or "height")
matrix.w: 2 (number of columns or "width")
matrix.e: 2D array of 2x2 elements, first index is row, second is column

*) Create a complex matrix [[1+2i], [3+4i], [5+6i], [7+8i]]
var matrix = new JmatM(2, 2);
matrix.e[0][0] = new JmatC(1, 2);
matrix.e[0][1] = new JmatC(3, 4);
matrix.e[1][0] = new JmatC(5, 6);
matrix.e[1][1] = new JmatC(7, 8);
shorthand for all the above:
var matrix = jmatm([[jmatc(1, 2), jmatc(3, 4)], [jmatc(5, 6), jmatc(7, 8)]])
var matrix = jmatm([[1, 3], [5, 7]], [[2, 4], [6, 8]])
var matrix = jmatm(2, 2, jmatc(1, 2), jmatc(3, 4), jmatc(5, 6), jmatc(7, 8))

*) Calculating with matrices
jmatm([[1, 2], [3, 4]]).mul(jmatm([[5, 6], [7, 8]])).toString()
--> [[19, 22], [43, 50]]
jmatm([[1, 2], [3, 4]]).add(jmatm([[5, 6], [7, 8]]).mulc(JmatC.I)).toString()
--> [[1+5i, 2+6i], [3+7i, 4+8i]]

*) Loop through the elements of a matrix
var a = jmatm([[1, 2], [3, 4]])
for(var y = 0; y < a.h; y++) { // row index
  for(var x = 0; x < a.w; x++) { // column index
    var complex = a.e[y][x]; // complex is of type JmatC
    console.log('x: ' + x + ' y: ' + y + ' re: ' + complex.re + ' im: ' + complex.im);
  }
}

*) Inverse of a matrix
JmatM.inv(jmatm([[1,2],[3,4]])).toString()
more convenient (but slightly less efficient):
Jmat.inv([[1,2],[3,4]]).toString()
--> [[-2, 1], [1.5, -0.5]]

*) SVD of a matrix
var svd = Jmat.svd([[1,2],[3,4]]);
'V:' + svd.v.toString() + ' S:' + svd.s.toString() + ' U:' + svd.u.toString();
--> V:[[0.576, -0.8174], [0.8174, 0.576]] S:[[5.465, 0], [0, 0.366]] U:[[0.4046, 0.9145], [0.9145, -0.4046]]

*) eigenvectors/values of a matrix:
var eig = Jmat.eig([[1,2],[3,4]]);
'values:' + eig.l.toString() + ' vectors:' + eig.v.toString();
--> values:[[5.3723], [-0.3723]] vectors:[[1, 1], [2.1861, -0.6861]]

*) submatrix
Jmat.submatrix([[1,2,3],[3,4,5],[6,7,8]], 0, 2, 0, 2).toString()
--> [[1, 2], [3, 4]]

*) Numerical integration (quadrature) of x^2 from 0 to 10 with 30 steps:
Jmat.integrate(jmatc(0), jmatc(10), function(x) { return x.mul(x); }, 30).toString();
--> 333.3333333333333

*) JmatR contains functions which operate on regular JS numbers, so could be seen as an extension to the standard JS Math library:
JmatR.gamma(Math.cos(2)) + JmatR.erfc(0.5) * JmatR.EM
--> -3.39388952436638
JmatR.nextPrime(17)
--> 19

*) Numerical integration (quadrature) of x^2 from 0 to 10 with 30 steps, with real numbers (no JmatC objects):
JmatR.integrate(0, 10, function(x) { return x * x; }, 30);
--> 333.3333333333333

*) Plotting with jmat_plot.js
Jmat.plot2D(Jmat.gamma_p, document.body);
Jmat.plotComplex(function(z) { return Jmat.polygamma(4, z); }, document.body, {p:1});
Jmat.plotReal(function(x) { return Jmat.cdf_studentt(x, 2); }, document.body, {p:1});
NOTE: document.body can be any parent element, e.g. a div, preferably further away from the left side to see more of left labels
NOTE: {p:1} sets to highest resolution
NOTE: if the arrows look like garbled text characters, specifying UTF-8 charset in the HTML page is needed.

*) Plot types:
-Real: function of 1 argument plotted with classic X/Y real plot. NOTE: other colors than black indicate NaN or complex values.
-Complex: function of 1 argument plotted in 2D for complex arguments with complex color wheel.
-2D: function of 2 arguments plotted for real arguments in 2D with result as complex color wheel.
Color wheel legend: red = positive, cyan = negative, black = zero, white = infinity, darker = lower abs,
    lighter = higher abs, acid green = positive imaginary, purple = negative imaginary,
    other hues = other complex arguments, grey = NaN


Precision
---------

The algorithms chosen are as precise as possible while being fast and simple.
However, this library focuses more on supporting the full complex domain for the
input of the functions. Often, approximations are used. Number of correct digits
is not specified. This means the precision may be low in some regions,
especially for extreme input values. In general, though, most functions are
reasonably precise for double precision floating point. If precision is
important, please check the results first in the applicable input domain.


Contact
-------

Feel free to contact me about bugs, improvements, comments, ...

My email address is first name DOT last name AT gmail DOT com, and my first and
last name are: Lode Vandevenne.
*/

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Jmat main API
////////////////////////////////////////////////////////////////////////////////

// Constructor
function Jmat() {
  // Emtpy, this is a namespace, no need to ever call this
}

/*
The Jmat functions below are the most public API (as well as any JmatR, JmatC or
JmatM functions with the same name).

Unlike their JmatC and JmatM counterparts, these versions are more tolerant to
various input types, e.g. Jmat.gamma can take both 5 and jmatc(5) as input,
while JmatC.gamma only accepts jmatc(5) and JmatR.gamma only accepts 5.

To work with complex numbers or matrices, you'll need to use JmatC and JmatM
objects anyway, create them with jmatc or jmatm.

Return values should always be treated as immutable objects. Do not modify the
re and im fields directly: doing so could result in changing the internal
constants (like JmatC.PI).

The JmatC and JmatM objects contain several functions in their prototype as
well, e.g. add, so it is possible to write jmatc(5).add(jmatc(6)) instead of
Jmat.add(jmatc(5), jmatc(6)). However, those prototype functions do not accept
input of a wrong type (must be JmatC, not a plain JS number).

The comments below use closure-style types, e.g. {number|JmatC} means the type
is either a plain JS number of a JmatC object, and {Array.<JmatC>} means an array
of JmatC objects.
*/

// Elementary operators

/* Add. x,y:{number|JmatC|JmatM}. returns {JmatC|JmatM}. */
Jmat.add = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return JmatM.add(JmatM.cast(x), JmatM.cast(y));
  return JmatC.add(JmatC.cast(x), JmatC.cast(y));
};
/* Subtract. x,y:{number|JmatC|JmatM}. returns {JmatC|JmatM}. */
Jmat.sub = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return JmatM.sub(JmatM.cast(x), JmatM.cast(y));
  return JmatC.sub(JmatC.cast(x), JmatC.cast(y));
};
/* Multiply. x,y:{number|JmatC|JmatM}. returns {JmatC|JmatM}. */
Jmat.mul = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return JmatM.mul(JmatM.cast(x), JmatM.cast(y));
  if(Jmat.matrixIn_(x)) return JmatM.mulc(JmatM.cast(x), JmatC.cast(y));
  if(Jmat.matrixIn_(y)) return JmatM.mulc(JmatM.cast(y), JmatC.cast(x));
  return JmatC.mul(JmatC.cast(x), JmatC.cast(y));
};
/* Division. x,y:{number|JmatC|JmatM}. returns {JmatC|JmatM}. */
Jmat.div = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return JmatM.div(JmatM.cast(x), JmatM.cast(y));
  if(Jmat.matrixIn_(x)) return JmatM.divc(JmatM.cast(x), JmatC.cast(y));
  if(Jmat.matrixIn_(y)) return JmatM.divc(JmatM.cast(y), JmatC.cast(x));
  return JmatC.div(JmatC.cast(x), JmatC.cast(y));
};

// Compare

/* Equal? x,y:{number|JmatC|JmatM}. returns {boolean}. */
Jmat.eq = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return JmatM.eq(JmatM.cast(x), JmatM.cast(y));
  return JmatC.eq(JmatC.cast(x), JmatC.cast(y));
};
/* Nearly equal? x,y:{number|JmatC|JmatM}. epsilon:{number}. returns {boolean}. */
Jmat.near = function(x, y, epsilon) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return JmatM.near(JmatM.cast(x), JmatM.cast(y), JmatR.caststrict(epsilon));
  return JmatC.near(JmatC.cast(x), JmatC.cast(y), JmatR.caststrict(epsilon));
};

// Power & Logarithms

/* Power, or matrix power. x:{number|JmatC|JmatM}. y:{number|JmatC}. returns {JmatC}. */
Jmat.pow = function(x, y) {
  if(Jmat.matrixIn_(x)) return JmatM.powc(JmatM.cast(x), JmatC.cast(y)); // matrix raised to any complex number power
  return JmatC.pow(JmatC.cast(x), JmatC.cast(y));
};
/* Square root, or matrix square root. z:{number|JmatC|JmatM}. returns {JmatC|JmatM}. */
Jmat.sqrt = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.sqrt(JmatM.cast(z));
  return JmatC.sqrt(JmatC.cast(z));
};
/* Exponential, or matrix exponential. z:{number|JmatC|JmatM}. returns {JmatC|JmatM}. */
Jmat.exp = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.exp(JmatM.cast(z)); // matrix exponential
  return JmatC.exp(JmatC.cast(z));
};
/* Logarithm, or matrix logarithm. z:{number|JmatC|JmatM}. returns {JmatC|JmatM}. */
Jmat.log = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.log(JmatM.cast(z));
  return JmatC.log(JmatC.cast(z));
};
/* log(z) + 1. z:{number|JmatC}. returns {JmatC}. */
Jmat.log1p = function(z) { return JmatC.log1p(JmatC.cast(z)); };
/* exp(z) - 1. z:{number|JmatC}. returns {JmatC}. */
Jmat.expm1 = function(z) { return JmatC.expm1(JmatC.cast(z)); };
/* Base-y logarithm of x. x,y:{number|JmatC}. returns {JmatC}. */
Jmat.logy = function(x, y) { return JmatC.logy(JmatC.cast(x), JmatC.cast(y)); };
/* Principal branch of LambertW. z:{number|JmatC}. returns {JmatC}. */
Jmat.lambertw = function(z) { return JmatC.lambertw(JmatC.cast(z)); };
/* Negative branch of LambertW. z:{number|JmatC}. returns {JmatC}. */
Jmat.lambertwm = function(z) { return JmatC.lambertwm(JmatC.cast(z)); };
/* Specific branch of LambertW. branch:{number} must be integer, z:{number|JmatC}. returns {JmatC}. */
Jmat.lambertwb = function(branch, z) { return JmatC.lambertwb(JmatR.caststrict(branch), JmatC.cast(z)); };
/* Tetration (power tower). x,y:{number|JmatC}. returns {JmatC} */
Jmat.tetration = function(x, y) { return JmatC.tetration(JmatC.cast(x), JmatC.cast(y)); };

// Elementary functions

/* Identity function. x:{number|JmatC}. returns {JmatC} */
Jmat.x = function(x) { return JmatC.cast(x); };
/* Negate. z:{number|JmatC|JmatM}. returns {JmatC|JmatM} */
Jmat.neg = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.neg(JmatM.cast(z));
  return JmatC.neg(JmatC.cast(z));
};
/* Reciproke. z:{number|JmatC|JmatM}. returns {JmatC|JmatM} */
Jmat.inv = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.inv(JmatM.cast(z));
  return JmatC.inv(JmatC.cast(z));
};
/* Complex conjugate. z:{number|JmatC|JmatM}. returns {JmatC|JmatM} */
Jmat.conj = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.conj(JmatM.cast(z));
  return JmatC.conj(JmatC.cast(z));
};

// Trigonometric functions

/* Sine, or matrix-sine. z:{number|JmatC|JmatM}. returns {JmatC|JmatM} */
Jmat.sin = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.sin(JmatM.cast(z));
  return JmatC.sin(JmatC.cast(z));
};
/* Cosine, or matrix-cosine. z:{number|JmatC|JmatM}. returns {JmatC|JmatM} */
Jmat.cos = function(z) {
  if(Jmat.matrixIn_(z)) return JmatM.cos(JmatM.cast(z));
  return JmatC.cos(JmatC.cast(z));
};
/* Tangent. z:{number|JmatC}. returns {JmatC} */
Jmat.tan = function(z) { return JmatC.tan(JmatC.cast(z)); };
/* Arcsine. z:{number|JmatC}. returns {JmatC} */
Jmat.asin = function(z) { return JmatC.asin(JmatC.cast(z)); };
/* Arccosine. z:{number|JmatC}. returns {JmatC} */
Jmat.acos = function(z) { return JmatC.acos(JmatC.cast(z)); };
/* Arctangent. z:{number|JmatC}. returns {JmatC} */
Jmat.atan = function(z) { return JmatC.atan(JmatC.cast(z)); };
/* Atan2 function. z:{number|JmatC}. returns {JmatC} */
Jmat.atan2 = function(x, y) { return JmatC.atan2(JmatC.cast(x), JmatC.cast(y)); };
/* Unnormalized sinc: sin(z)/z. z:{number|JmatC}. returns {JmatC} */
Jmat.sinc = function(z) { return JmatC.sinc(JmatC.cast(z)); };

// Hyperbolic functions

/* Hyperbolic sine. z:{number|JmatC}. returns {JmatC} */
Jmat.sinh = function(z) { return JmatC.sinh(JmatC.cast(z)); };
/* Hyperbolic cosine. z:{number|JmatC}. returns {JmatC} */
Jmat.cosh = function(z) { return JmatC.cosh(JmatC.cast(z)); };
/* Hyperbolic tangent. z:{number|JmatC}. returns {JmatC} */
Jmat.tanh = function(z) { return JmatC.tanh(JmatC.cast(z)); };
/* Hyperbolic arcsine. z:{number|JmatC}. returns {JmatC} */
Jmat.asinh = function(z) { return JmatC.asinh(JmatC.cast(z)); };
/* Hyperbolic arccosine. z:{number|JmatC}. returns {JmatC} */
Jmat.acosh = function(z) { return JmatC.acosh(JmatC.cast(z)); };
/* Hyperbolic arctangent. z:{number|JmatC}. returns {JmatC} */
Jmat.atanh = function(z) { return JmatC.atanh(JmatC.cast(z)); };

// Gamma and related functions

/* Gamma function. z:{number|JmatC}. returns {JmatC} */
Jmat.gamma = function(z) { return JmatC.gamma(JmatC.cast(z)); };
/* Factorial. z:{number|JmatC}. returns {JmatC} */
Jmat.factorial = function(z) { return JmatC.factorial(JmatC.cast(z)); };
/* Digamma function (psi). z:{number|JmatC}. returns {JmatC} */
Jmat.digamma = function(z) { return JmatC.digamma(JmatC.cast(z)); };
/* Trigamma function. z:{number|JmatC}. returns {JmatC} */
Jmat.trigamma = function(z) { return JmatC.trigamma(JmatC.cast(z)); };
/* Polygamma function. n,z:{number|JmatC}. returns {JmatC} */
Jmat.polygamma = function(n, z) { return JmatC.polygamma(JmatC.cast(n), JmatC.cast(z)); };
/* Logarithm of gamma function. z:{number|JmatC}. returns {JmatC} */
Jmat.loggamma = function(z) { return JmatC.loggamma(JmatC.cast(z)); };
/* Inverse of gamma function (not reciproke). z:{number|JmatC}. returns {JmatC} */
Jmat.gamma_inv = function(z) { return JmatC.gamma_inv(JmatC.cast(z)); };
/* Lower incomplete gamma function. s,z:{number|JmatC}. returns {JmatC} */
Jmat.incgamma_lower = function(s, z) { return JmatC.incgamma_lower(JmatC.cast(s), JmatC.cast(z)); };
/* Upper incomplete gamma function. s,z:{number|JmatC}. returns {JmatC} */
Jmat.incgamma_upper = function(s, z) { return JmatC.incgamma_upper(JmatC.cast(s), JmatC.cast(z)); };
/* Lower regularized incomplete gamma function "P". s,z:{number|JmatC}. returns {JmatC} */
Jmat.gamma_p = function(s, z) { return JmatC.gamma_p(JmatC.cast(s), JmatC.cast(z)); };
/* Upper regularized incomplete gamma function. s,z:{number|JmatC}. returns {JmatC} */
Jmat.gamma_q = function(s, z) { return JmatC.gamma_q(JmatC.cast(s), JmatC.cast(z)); };
/* Inverse of "P" in z. s,o:{number|JmatC}. returns {JmatC} */
Jmat.gamma_p_inv = function(s, p) { return JmatC.gamma_p_inv(JmatC.cast(s), JmatC.cast(p)); };
/* Beta function. x,y:{number|JmatC}. returns {JmatC} */
Jmat.beta = function(x, y) { return JmatC.beta(JmatC.cast(x), JmatC.cast(y)); };
/* Incomplete beta function. x,a,b:{number|JmatC}. returns {JmatC} */
Jmat.incbeta = function(x, a, b) { return JmatC.incbeta(JmatC.cast(x), JmatC.cast(a), JmatC.cast(b)); };
/* Regularized incomplete beta function "I". x,a,b:{number|JmatC}. returns {JmatC} */
Jmat.beta_i = function(x, a, b) { return JmatC.beta_i(JmatC.cast(x), JmatC.cast(a), JmatC.cast(b)); };
/* Inverse of "I" in x. x,a,b:{number|JmatC}. returns {JmatC} */
Jmat.beta_i_inv = function(x, a, b) { return JmatC.beta_i_inv(JmatC.cast(x), JmatC.cast(a), JmatC.cast(b)); };

// Number and complex number utility functions

/* Real part. z:{number|JmatC}. returns {JmatC} */
Jmat.re = function(z) { return jmatc(JmatC.cast(z).re); };
/* Imaginary part. z:{number|JmatC}. returns {JmatC} */
Jmat.im = function(z) { return jmatc(JmatC.cast(z).im); };
/* Absolute value or complex odulus. z:{number|JmatC}. returns {JmatC} */
Jmat.abs = function(z) { return JmatC.abs(JmatC.cast(z)); }; // absolute value, modulus
/* Complex argument or phase. z:{number|JmatC}. returns {JmatC} */
Jmat.arg = function(z) { return JmatC.arg(JmatC.cast(z)); };
/* Real sign. z:{number|JmatC}. returns {JmatC} */
Jmat.sign = function(z) { return JmatC.sign(JmatC.cast(z)); };
/* Complex sign. z:{number|JmatC}. returns {JmatC} */
Jmat.csgn = function(z) { return JmatC.csgn(JmatC.cast(z)); };
/* Floor. x:{number|JmatC}. returns {JmatC} */
Jmat.floor = function(x) { return JmatC.floor(JmatC.cast(x)); };
/* Ceiling. x:{number|JmatC}. returns {JmatC} */
Jmat.ceil = function(x) { return JmatC.ceil(JmatC.cast(x)); };
/* Round to integer. x:{number|JmatC}. returns {JmatC} */
Jmat.round = function(x) { return JmatC.round(JmatC.cast(x)); };
/* Truncate towards zero. x:{number|JmatC}. returns {JmatC} */
Jmat.trunc = function(x) { return JmatC.trunc(JmatC.cast(x)); };
/* Fractional part, always positive. x:{number|JmatC}. returns {JmatC} */
Jmat.frac = function(x) { return JmatC.frac(JmatC.cast(x)); };
/* Fractional part, negative for negative x. x:{number|JmatC}. returns {JmatC} */
Jmat.fracn = function(x) { return JmatC.fracn(JmatC.cast(x)); };
/* Rotate complex number by a radians. z:{number|JmatC}, a:{number}. returns {JmatC} */
Jmat.rotate = function(z, a) { return JmatC.rotate(JmatC.cast(z), JmatR.caststrict(a)); };
/* Approximate with integer numerator and denominator. x:{number|JmatC}, max:{number} maximum denominator. returns {Array.<JmatC>} numerator, denominator */
Jmat.decompose = function(x, max) { return JmatC.decompose(JmatC.cast(x), JmatR.cast(max)); };

// Cylindrical functions

/* Bessel function of the first kind.  n,z:{number|JmatC}. returns {JmatC} */
Jmat.besselj = function(n, z) { return JmatC.besselj(JmatC.cast(n), JmatC.cast(z)); };
/* Bessel function of the second kind.  n,z:{number|JmatC}. returns {JmatC} */
Jmat.bessely = function(n, z) { return JmatC.bessely(JmatC.cast(n), JmatC.cast(z)); };
/* Modified bessel function of the first kind.  n,z:{number|JmatC}. returns {JmatC} */
Jmat.besseli = function(n, z) { return JmatC.besseli(JmatC.cast(n), JmatC.cast(z)); };
/* Modified bessel function of the second kind.  n,z:{number|JmatC}. returns {JmatC} */
Jmat.besselk = function(n, z) { return JmatC.besselk(JmatC.cast(n), JmatC.cast(z)); };
/* Hankel function of the first kind.  n,z:{number|JmatC}. returns {JmatC} */
Jmat.hankel1 = function(n, z) { return JmatC.hankel1(JmatC.cast(n), JmatC.cast(z)); };
/* Hankel function of the second kind.  n,z:{number|JmatC}. returns {JmatC} */
Jmat.hankel2 = function(n, z) { return JmatC.hankel2(JmatC.cast(n), JmatC.cast(z)); };
/* Airy function Ai. z:{number|JmatC}. returns {JmatC} */
Jmat.airy = function(z) { return JmatC.airy(JmatC.cast(z)); };
/* Bairy function Bi. z:{number|JmatC}. returns {JmatC} */
Jmat.bairy = function(z) { return JmatC.bairy(JmatC.cast(z)); };
/* Derivative of Airy function Ai. z:{number|JmatC}. returns {JmatC} */
Jmat.airy_deriv = function(z) { return JmatC.airy_deriv(JmatC.cast(z)); };
/* Derivative of Bairy function Bi. z:{number|JmatC}. returns {JmatC} */
Jmat.bairy_deriv = function(z) { return JmatC.bairy_deriv(JmatC.cast(z)); };

// Theta functions

/* Jacobi theta function θ1. z,q:{number|JmatC}. returns {JmatC} */
Jmat.theta1 = function(z, q) { return JmatC.theta1(JmatC.cast(z), JmatC.cast(q)); };
/* Jacobi theta function θ2. z,q:{number|JmatC}. returns {JmatC} */
Jmat.theta2 = function(z, q) { return JmatC.theta2(JmatC.cast(z), JmatC.cast(q)); };
/* Jacobi theta function θ3. z,q:{number|JmatC}. returns {JmatC} */
Jmat.theta3 = function(z, q) { return JmatC.theta3(JmatC.cast(z), JmatC.cast(q)); };
/* Jacobi theta function θ4. z,q:{number|JmatC}. returns {JmatC} */
Jmat.theta4 = function(z, q) { return JmatC.theta4(JmatC.cast(z), JmatC.cast(q)); };

// Zeta and related functions

/* Riemann zeta function. z:{number|JmatC}. returns {JmatC} */
Jmat.zeta = function(z) { return JmatC.zeta(JmatC.cast(z)); };
/* Dirichlet eta function. z:{number|JmatC}. returns {JmatC} */
Jmat.eta = function(z) { return JmatC.eta(JmatC.cast(z)); };
/* Dirichlet lambda function. z:{number|JmatC}. returns {JmatC} */
Jmat.lambda = function(z) { return JmatC.lambda(JmatC.cast(z)); };
/* Dilogarithm. z:{number|JmatC}. returns {JmatC} */
Jmat.dilog = function(z) { return JmatC.dilog(JmatC.cast(z)); };
/* Trilogarithm. z:{number|JmatC}. returns {JmatC} */
Jmat.trilog = function(z) { return JmatC.trilog(JmatC.cast(z)); };
/* Polylogarithm. s,z:{number|JmatC}. returns {JmatC} */
Jmat.polylog = function(s, z) { return JmatC.polylog(JmatC.cast(s), JmatC.cast(z)); };
/* Hurwitz Zeta function. s,q:{number|JmatC}. returns {JmatC} */
Jmat.hurwitzzeta = function(s, q) { return JmatC.hurwitzzeta(JmatC.cast(s), JmatC.cast(q)); };

// Hypergeometric functions

/* Hypergeometric function 0F1. a,z:{number|JmatC}. returns {JmatC} */
Jmat.hypergeometric0F1 = function(a, z) { return JmatC.hypergeometric0F1(JmatC.cast(a), JmatC.cast(z)); };
/* Confluent Hypergeometric function (Kummer). a,b,z:{number|JmatC}. returns {JmatC} */
Jmat.hypergeometric1F1 = function(a, b, z) { return JmatC.hypergeometric1F1(JmatC.cast(a), JmatC.cast(b), JmatC.cast(z)); };
/* Hypergeometric function 2F1. a,b,c,z:{number|JmatC}. returns {JmatC} */
Jmat.hypergeometric = function(a, b, c, z) { return JmatC.hypergeometric(JmatC.cast(a), JmatC.cast(b), JmatC.cast(c), JmatC.cast(z)); };

// Error and related functions

/* Error function. z:{number|JmatC}. returns {JmatC} */
Jmat.erf = function(z) { return JmatC.erf(JmatC.cast(z)); };
/* Complementary error function. z:{number|JmatC}. returns {JmatC} */
Jmat.erfc = function(z) { return JmatC.erfc(JmatC.cast(z)); };
/* Scaled complementary error function. z:{number|JmatC}. returns {JmatC} */
Jmat.erfcx = function(z) { return JmatC.erfcx(JmatC.cast(z)); };
/* Inverse of error function (not reciproke). z:{number|JmatC}. returns {JmatC} */
Jmat.erf_inv = function(z) { return JmatC.erf_inv(JmatC.cast(z)); };
/* Inverse of complementary error function (not reciproke). z:{number|JmatC}. returns {JmatC} */
Jmat.erfc_inv = function(z) { return JmatC.erfc_inv(JmatC.cast(z)); };
/* Imaginary error function. x:{number|JmatC}. returns {JmatC} */
Jmat.erfi = function(x) { return x.re == undefined ? jmatc(JmatR.erfi(x)) : JmatC.erfi(JmatC.cast(x)); };
/* Dawson function D+(x). x:{number|JmatC}. returns {JmatC} */
Jmat.dawson = function(x) { return x.re == undefined ? jmatc(JmatR.dawson(x)) : JmatC.dawson(JmatC.cast(x)); };
/* Faddeeva function w(z). z:{number|JmatC}. returns {JmatC} */
Jmat.faddeeva = function(z) { return JmatC.faddeeva(JmatC.cast(z)); };

// Statistical distributions. pdf = probability density function, cdf = cumulative distribution function, qf = quantile function

/* Uniform distribution in range a-b. x,a,b:{number|JmatC}. returns {JmatC} */
Jmat.pdf_uniform = function(x, a, b) { return JmatC.pdf_uniform(JmatC.cast(x), JmatC.cast(a), JmatC.cast(b)); };
Jmat.cdf_uniform = function(x, a, b) { return JmatC.cdf_uniform(JmatC.cast(x), JmatC.cast(a), JmatC.cast(b)); };
Jmat.qf_uniform = function(x, a, b) { return JmatC.qf_uniform(JmatC.cast(x), JmatC.cast(a), JmatC.cast(b)); };
/* Standardnormal distribution. x:{number|JmatC}. returns {JmatC} */
Jmat.pdf_standardnormal = function(x) { return JmatC.pdf_standardnormal(JmatC.cast(x)); };
Jmat.cdf_standardnormal = function(x) { return JmatC.cdf_standardnormal(JmatC.cast(x)); };
Jmat.qf_standardnormal = function(x) { return JmatC.qf_standardnormal(JmatC.cast(x)); };
/* Normal distribution with mean mu and variance signa. x,mu,sigma:{number|JmatC}. returns {JmatC} */
Jmat.pdf_normal = function(x, mu, sigma) { return JmatC.pdf_normal(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(sigma)); };
Jmat.cdf_normal = function(x, mu, sigma) { return JmatC.cdf_normal(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(sigma)); };
Jmat.qf_normal = function(x, mu, sigma) { return JmatC.qf_normal(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(sigma)); };
/* Log-normal distribution with mean mu and variance signa. x,mu,sigma:{number|JmatC}. returns {JmatC} */
Jmat.pdf_lognormal = function(x, mu, sigma) { return JmatC.pdf_lognormal(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(sigma)); };
Jmat.cdf_lognormal = function(x, mu, sigma) { return JmatC.cdf_lognormal(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(sigma)); };
Jmat.qf_lognormal = function(x, mu, sigma) { return JmatC.qf_lognormal(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(sigma)); };
/* Cauchy distribution with location x0 and scale gamma. x,x0,gamma:{number|JmatC}. returns {JmatC} */
Jmat.pdf_cauchy = function(x, x0, gamma) { return JmatC.pdf_cauchy(JmatC.cast(x), JmatC.cast(x0), JmatC.cast(gamma)); };
Jmat.cdf_cauchy = function(x, x0, gamma) { return JmatC.cdf_cauchy(JmatC.cast(x), JmatC.cast(x0), JmatC.cast(gamma)); };
Jmat.qf_cauchy = function(x, x0, gamma) { return JmatC.qf_cauchy(JmatC.cast(x), JmatC.cast(x0), JmatC.cast(gamma)); };
/* Student's t distribution with degrees of freedom nu. x,nu:{number|JmatC}. returns {JmatC} */
Jmat.pdf_studentt = function(x, nu) { return JmatC.pdf_studentt(JmatC.cast(x), JmatC.cast(nu)); };
Jmat.cdf_studentt = function(x, nu) { return JmatC.cdf_studentt(JmatC.cast(x), JmatC.cast(nu)); };
Jmat.qf_studentt = function(x, nu) { return JmatC.qf_studentt(JmatC.cast(x), JmatC.cast(nu)); };
/* Chi square distribution with degrees of freedom k. x,k:{number|JmatC}. returns {JmatC} */
Jmat.pdf_chi_square = function(x, k) { return JmatC.pdf_chi_square(JmatC.cast(x), JmatC.cast(k)); };
Jmat.cdf_chi_square = function(x, k) { return JmatC.cdf_chi_square(JmatC.cast(x), JmatC.cast(k)); };
Jmat.qf_chi_square = function(x, k) { return JmatC.qf_chi_square(JmatC.cast(x), JmatC.cast(k)); };
/* Logistic distribution with location mu and scale s. x,mu,s:{number|JmatC}. returns {JmatC} */
Jmat.pdf_logistic = function(x, mu, s) { return JmatC.pdf_logistic(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(s)); };
Jmat.cdf_logistic = function(x, mu, s) { return JmatC.cdf_logistic(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(s)); };
Jmat.qf_logistic = function(x, mu, s) { return JmatC.qf_logistic(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(s)); };
/* Gamma distribution with shape k and scale theta. x,k,theta:{number|JmatC}. returns {JmatC} */
Jmat.pdf_gamma = function(x, k, theta) { return JmatC.pdf_gamma(JmatC.cast(x), JmatC.cast(k), JmatC.cast(theta)); };
Jmat.cdf_gamma = function(x, k, theta) { return JmatC.cdf_gamma(JmatC.cast(x), JmatC.cast(k), JmatC.cast(theta)); };
Jmat.qf_gamma = function(x, k, theta) { return JmatC.qf_gamma(JmatC.cast(x), JmatC.cast(k), JmatC.cast(theta)); };
/* Beta distribution with shape alpha and beta. x,alpha,beta:{number|JmatC}. returns {JmatC} */
Jmat.pdf_beta = function(x, alpha, beta) { return JmatC.pdf_beta(JmatC.cast(x), JmatC.cast(alpha), JmatC.cast(beta)); };
Jmat.cdf_beta = function(x, alpha, beta) { return JmatC.cdf_beta(JmatC.cast(x), JmatC.cast(alpha), JmatC.cast(beta)); };
Jmat.qf_beta = function(x, alpha, beta) { return JmatC.qf_beta(JmatC.cast(x), JmatC.cast(alpha), JmatC.cast(beta)); };
/* F-distribution with degrees of freedom d1 and d2. x,d1,d2:{number|JmatC}. returns {JmatC} */
Jmat.pdf_fisher = function(x, d1, d2) { return JmatC.pdf_fisher(JmatC.cast(x), JmatC.cast(d1), JmatC.cast(d2)); };
Jmat.cdf_fisher = function(x, d1, d2) { return JmatC.cdf_fisher(JmatC.cast(x), JmatC.cast(d1), JmatC.cast(d2)); };
Jmat.qf_fisher = function(x, d1, d2) { return JmatC.qf_fisher(JmatC.cast(x), JmatC.cast(d1), JmatC.cast(d2)); };
/* Weibull distribution with scale lambda, shape k. x,lambda,k:{number|JmatC}. returns {JmatC} */
Jmat.pdf_weibull = function(x, lambda, k) { return JmatC.pdf_weibull(JmatC.cast(x), JmatC.cast(lambda), JmatC.cast(k)); };
Jmat.cdf_weibull = function(x, lambda, k) { return JmatC.cdf_weibull(JmatC.cast(x), JmatC.cast(lambda), JmatC.cast(k)); };
Jmat.qf_weibull = function(x, lambda, k) { return JmatC.qf_weibull(JmatC.cast(x), JmatC.cast(lambda), JmatC.cast(k)); };
/* Exponential distribution with rate lambda. x,lambda:{number|JmatC}. returns {JmatC} */
Jmat.pdf_exponential = function(x, lambda) { return JmatC.pdf_exponential(JmatC.cast(x), JmatC.cast(lambda)); };
Jmat.cdf_exponential = function(x, lambda) { return JmatC.cdf_exponential(JmatC.cast(x), JmatC.cast(lambda)); };
Jmat.qf_exponential = function(x, lambda) { return JmatC.qf_exponential(JmatC.cast(x), JmatC.cast(lambda)); };
/* Laplace distribution with location mu, scale b. x,mu,b:{number|JmatC}. returns {JmatC} */
Jmat.pdf_laplace = function(x, mu, b) { return JmatC.pdf_laplace(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(b)); };
Jmat.cdf_laplace = function(x, mu, b) { return JmatC.cdf_laplace(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(b)); };
Jmat.qf_laplace = function(x, mu, b) { return JmatC.qf_laplace(JmatC.cast(x), JmatC.cast(mu), JmatC.cast(b)); };

// Combinatorics

/* Number of permutations. n,p:{number|JmatC}. returns {JmatC} */
Jmat.permutation = function(n, p) { return JmatC.permutation(JmatC.cast(n), JmatC.cast(p)); };
/* Binomial, number of combinations. n,p:{number|JmatC}. returns {JmatC} */
Jmat.binomial = function(n, p) { return JmatC.binomial(JmatC.cast(n), JmatC.cast(p)); };
/* Stirling number of hte second kind. n{number|JmatC} integer, p:{number|JmatC}. returns {JmatC} */
Jmat.stirling2 = function(n, k) { return JmatC.stirling2(JmatC.cast(n), JmatC.cast(k)); };

// Mean

/* Arithmetic-Geometric mean. a,b:{number|JmatC}. returns {JmatC} */
Jmat.agm = function(a, b) { return JmatC.agm(JmatC.cast(a), JmatC.cast(b)); };
/* Geometric-Harmonic mean. a,b:{number|JmatC}. returns {JmatC} */
Jmat.ghm = function(a, b) { return JmatC.ghm(JmatC.cast(a), JmatC.cast(b)); };

// Prime numbers. NOTE: Use the JmatR versions to get results as simple JS numbers rather than JmatC objects.

/* Primality test. x:{number|JmatC}. returns {boolean} */
Jmat.isPrime = function(x) { return JmatC.isPrime(JmatC.cast(x)); }
/* Smallest prime factor of integer. x:{number|JmatC} real. returns {JmatC} */
Jmat.smallestPrimeFactor = function(x) { return jmatc(JmatR.smallestPrimeFactor(JmatR.cast(x))); };
/* Factorize into prime factors, returned as array. x:{number|JmatC} integer. returns {Array.<JmatC>} */
Jmat.factorize = function(n) {
  var res = JmatR.factorize(JmatR.caststrict(n));
  for(var i = 0; i < res.length; i++) res[i] = jmatc(res[i]); // Convert to JmatC objects to be consistent with other Jmat functions. Just use JmatR.factorize instead to have simple JS numbers.
  return res;
};
/* Euler's Totient. x:{number|JmatC} real integer. returns {JmatC} */
Jmat.eulerTotient = function(x) { return jmatc(JmatR.eulerTotient(JmatR.caststrict(x))); };
/* Prime counting function. x:{number|JmatC} real. returns {JmatC} */
Jmat.primeCount = function(x) { return jmatc(JmatR.primeCount(JmatR.caststrict(x))); };
/* Nearest prime number. x:{number|JmatC}. returns {JmatC} */
Jmat.nearestPrime = function(x) { return jmatc(JmatR.nearestPrime(JmatR.cast(x))); };
/* Next larger prime number function. x:{number|JmatC}. returns {JmatC} */
Jmat.nextPrime = function(x) { return jmatc(JmatR.nextPrime(JmatR.cast(x))); };
/* Next smaller prime number function. x:{number|JmatC}. returns {JmatC} */
Jmat.previousPrime = function(x) { return jmatc(JmatR.previousPrime(JmatR.cast(x))); };
/* Greatest common divisor. x,y:{number|JmatC} real integer. returns {JmatC} */
Jmat.gcd = function(x, y) { return jmatc(JmatR.gcd(JmatR.caststrict(x), JmatR.caststrict(y))); };
/* Least common multiple. x,y:{number|JmatC} real integer. returns {JmatC} */
Jmat.lcm = function(x, y) { return jmatc(JmatR.lcm(JmatR.caststrict(x), JmatR.caststrict(y))); };

// Sexagesimal

/* Decimal degrees to degrees, minutes, seconds. x:{number|JmatC} real. returns {JmatC} */
Jmat.dms = function(a) { return jmatc(JmatR.dms(JmatR.caststrict(a))); }; // E.g. 1.6 becomes 1.36 (1 degree, 36 minutes)
/* Degrees, minutes, seconds to decimal degrees. x:{number|JmatC} real. returns {JmatC} */
Jmat.dd = function(a) { return jmatc(JmatR.dd(JmatR.caststrict(a))); }; // E.g. 1.36 becomes 1.6 (1 point 6 degrees)

// Bitwise

/* Bitwise not. x:{number|JmatC}. returns {JmatC} */
Jmat.bitnot = function(x) { return JmatC.bitnot(JmatC.cast(x)); };
/* Bitwise and. x,y:{number|JmatC}. returns {JmatC} */
Jmat.bitand = function(x, y) { return JmatC.bitand(JmatC.cast(x), JmatC.cast(y)); };
/* Bitwise or. x,y:{number|JmatC}. returns {JmatC} */
Jmat.bitor = function(x, y) { return JmatC.bitor(JmatC.cast(x), JmatC.cast(y)); };
/* Bitwise xor. x,y:{number|JmatC}. returns {JmatC} */
Jmat.bitxor = function(x, y) { return JmatC.bitxor(JmatC.cast(x), JmatC.cast(y)); };
/* Left bitshift. x:{number|JmatC}, y:{number|JmatC} real. returns {JmatC} */
Jmat.lshift = function(x, y) { return JmatC.lshift(JmatC.cast(x), JmatR.caststrict(y)); };
/* Right bitshift. x:{number|JmatC}, y:{number|JmatC} real. returns {JmatC} */
Jmat.rshift = function(x, y) { return JmatC.rshift(JmatC.cast(x), JmatR.caststrict(y)); };

// Modulo division and related functions

/* Modulo division. a,b:{number|JmatC}. returns {JmatC} */
Jmat.mod = function(a, b) { return JmatC.mod(JmatC.cast(a), JmatC.cast(b)); }; // result has sign of divisor (unlike JS '%' operator)
/* Remainder. a,b:{number|JmatC}. returns {JmatC} */
Jmat.rem = function(a, b) { return JmatC.rem(JmatC.cast(a), JmatC.cast(b)); }; // result has sign of divident (same result as JS '%' operator on real numbers)
/* Wrap x between from and to (to excluded). x,to,from:{number|JmatC}. returns {JmatC} */
Jmat.wrap = function(x, from, to) { return JmatC.wrap(JmatC.cast(x), JmatC.cast(from), JmatC.cast(to)); };
/* Clamp x between from and to (to included). x,to,from:{number|JmatC}. returns {JmatC} */
Jmat.clamp = function(x, from, to) { return JmatC.clamp(JmatC.cast(x), JmatC.cast(from), JmatC.cast(to)); };

// Other special functions

/* Minkowski's question mark function. x:{number|JmatC}. returns {JmatC} */
Jmat.minkowski = function(x) { return JmatC.minkowski(JmatC.cast(x)); };

// Matrix (NOTE: more are above: add, sub, mul, div, inv, neg, conj, exp, log, sqrt, cos, sin)

/* Eigenvalues and eigenvectors. m:{Array|JmatM}. returns {Object.<string, JmatM>} object with l:eigenvalues, v:eigenvectors */
Jmat.eig = function(m) { return JmatM.eig(JmatM.cast(m)); };
/* Moore-Penrose pseudo-inverse. m:{Array|JmatM}. returns {JmatM} */
Jmat.pseudoinverse = function(m) { return JmatM.pseudoinverse(JmatM.cast(m)); };
/* Determinant. m:{Array|JmatM}. returns {JmatC} */
Jmat.determinant = function(m) { return JmatM.determinant(JmatM.cast(m)); };
/* Transpose. m:{Array|JmatM}. returns {JmatM} */
Jmat.transpose = function(m) { return JmatM.transpose(JmatM.cast(m)); };
/* Transjugate. m:{Array|JmatM}. returns {JmatM} */
Jmat.transjugate = function(m) { return JmatM.transjugate(JmatM.cast(m)); };
/* Trace. m:{Array|JmatM}. returns {JmatC} */
Jmat.trace = function(m) { return JmatM.trace(JmatM.cast(m)); };
/* Rank. m:{Array|JmatM}. returns {JmatC} */
Jmat.rank = function(m) { return JmatM.rank(JmatM.cast(m)); };
/* Adjoint aka Adjugate. m:{Array|JmatM}. returns {JmatM} */
Jmat.adj = function(m) { return JmatM.adj(JmatM.cast(m)); };
/* Fourier transform. m:{Array|JmatM}. returns {JmatM} */
Jmat.fft = function(m) { return JmatM.fft(JmatM.cast(m)); };
/* Inverse Fourier transform. m:{Array|JmatM}. returns {JmatM} */
Jmat.ifft = function(m) { return JmatM.ifft(JmatM.cast(m)); };
/* Solve AX=B. a,b:{Array|JmatM}. returns {JmatM} */
Jmat.solve = function(a, b) { return JmatM.solve(JmatM.cast(a), JmatM.cast(b)); };
/* Dot product of two vectors. a,b:{Array|JmatM}. returns {JmatC} */
Jmat.dot = function(a, b) { return JmatM.dot(JmatM.cast(a), JmatM.cast(b)); };
/* Cross product of two size-3 vectors. m:{Array|JmatM}. returns {JmatM} */
Jmat.cross = function(a, b) { return JmatM.cross(JmatM.cast(a), JmatM.cast(b)); };
/* Minor. a:{Array|JmatM}, row,col:{number} integer. returns {JmatC} */
Jmat.minor = function(a, row, col) { return JmatM.minor(JmatM.cast(a), JmatR.caststrict(row), JmatR.caststrict(col)); };
/* Cofactor. a:{Array|JmatM}, row,col:{number} integer. returns {JmatC} */
Jmat.cofactor = function(a, row, col) { return JmatM.cofactor(JmatM.cast(a), JmatR.caststrict(row), JmatR.caststrict(col)); };
/* Submatrix, x1 and y1 excluded, 0-based coordinates. a:{Array|JmatM}, y0,y1,x0,x1:{number} integer. returns {JmatM} */
Jmat.submatrix = function(a, y0, y1, x0, x1) { return JmatM.submatrix(JmatM.cast(a), JmatR.caststrict(y0), JmatR.caststrict(y1), JmatR.caststrict(x0), JmatR.caststrict(x1)); };

// Matrix decompositions

/* Singular value decomposition. m:{Array|JmatM}. returns {Object.<string, JmatM>} object with u:left vectors, s:singular values, v:right vectors */
Jmat.svd = function(m) { return JmatM.svd(JmatM.cast(m)); };
/* Spectral decomposition. m:{Array|JmatM}. returns {Object.<string, JmatM>} object with v:eigenvectors, d:eigenvalues on diagonal */
Jmat.eigd = function(m) { return JmatM.eigd(JmatM.cast(m)); };
/* QR decomposition. m:{Array|JmatM}. returns {Object.<string, JmatM>} object with q:Q matrix, r:R matrix */
Jmat.qr = function(m) { return JmatM.qr(JmatM.cast(m)); };

// Matrix norms

/* Frobenius norm. m:{Array|JmatM}. returns {JmatC} */
Jmat.norm = function(m) { return JmatM.norm(JmatM.cast(m)); };
/* Spectral norm. m:{Array|JmatM}. returns {JmatC} */
Jmat.norm2 = function(m) { return JmatM.norm2(JmatM.cast(m)); };
/* Maximum column norm. m:{Array|JmatM}. returns {JmatC} */
Jmat.maxcolnorm = function(m) { return JmatM.maxcolnorm(JmatM.cast(m)); };
/* Maximum row norm. m:{Array|JmatM}. returns {JmatC} */
Jmat.maxrownorm = function(m) { return JmatM.maxrownorm(JmatM.cast(m)); };

// Numerical algorithms. The function parameters f and df must always work with JmatC objects, one input, one output.

/* Derivative of f at x. x:{number|JmatC}, f:{function(JmatC):JmatC}. returns {JmatC} */
Jmat.differentiate = function(x, f) { return JmatC.differentiate(JmatC.cast(x), f); }; // derivative
/* Quadrature: definite integral of f from x to y. x,y:{number|JmatC}, f:{function(JmatC):JmatC}, steps:{number} integer. returns {JmatC} */
Jmat.integrate = function(x, y, f, steps) { return JmatC.integrate(JmatC.cast(x), JmatC.cast(y), f, JmatR.caststrict(steps)); };
/* Secant method. f:{function(JmatC):JmatC}, z0:{number|JmatC}, maxiter:{number} integer. returns {JmatC} */
Jmat.rootfind_secant = function(f, z0, maxiter) { return JmatC.rootfind_secant(f, JmatC.cast(z0), JmatR.caststrict(maxiter)); };
/* Newton's method. f,df:{function(JmatC):JmatC} df is derivative of f, z0:{number|JmatC}, maxiter:{number} integer. returns {JmatC} */
Jmat.rootfind_newton = function(f, df, z0, maxiter) { return JmatC.rootfind_secant(f, df, JmatC.cast(z0), JmatR.caststrict(maxiter)); };

////////////////////////////////////////////////////////////////////////////////

// Test if input is a matrix, this makes it decide e.g. whether to call the matrix-specific or complex number specific function if there is functionname collision
Jmat.matrixIn_ = function(v) {
  return v && (v instanceof JmatM || v.length != undefined);
};

// debug string of anything, like the arrays of matrices returned by svd
Jmat.dbg_ = function(a) {
  if(!a) return '' + a;
  var result = '';
  // For arrays of known types
  if(a.length) {
    result += '[';
    for(var i = 0; i < a.length; i++) result += (Jmat.dbg_(a[i]) + (i + 1 == a.length ? '' : ', '));
    result += ']';
    return result;
  }
  // For those objects like the result of SVD or eig
  if(a.constructor == Object) {
    result += '{';
    var comma = false;
    for(var it in a) {
      if(!it || !a[it]) continue;
      if(comma) result += ', ';
      result += it + ': ';
      result += Jmat.dbg_(a[it]);
      comma = true;
    }
    result += '}';
    return result;
  }
  // Prefer toString implementation if available (that is where this function is better than JSON.stringify! :))
  result = a.toString ? a.toString() : ('' + a);
  return result;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Real math: works on JS numbers
////////////////////////////////////////////////////////////////////////////////

/*
JmatR: real math. This contains only a small amount of functions, most is
implemented in JmatC (and JmatM for the matrices) further down.

These work on actual JS numbers rather than JmatC objects though, so JmatR is
as easy to use as the standard JS Math library. Functions existing in Math, such
as cos and exp, are not implemented here for that reason.
*/

//Constructor, but not to be actually used, just a namespace for real functions
function JmatR() {
}

// cast all known numeric types to JS number
JmatR.cast = function(v) {
  if(v && v.re != undefined) return v.re;
  if(v == undefined) return 0;
  return v;
};

// cast all known numeric types to JS number, but only if real (so complex/imag gives NaN)
JmatR.caststrict = function(v) {
  if(v && v.re != undefined) return v.im == 0 ? v.re : NaN;
  if(v == undefined) return 0;
  return v;
};

JmatR.SQRT2 = Math.sqrt(2);
JmatR.SQRTPI = Math.sqrt(Math.PI); // gamma(0.5)
JmatR.EM = 0.57721566490153286060; // Euler-Mascheroni constant
JmatR.APERY = 1.2020569; // Apery's constant, zeta(3)

////////////////////////////////////////////////////////////////////////////////
// Categories
////////////////////////////////////////////////////////////////////////////////

JmatR.isInt = function(x) {
  return x == Math.floor(x);
};

JmatR.isPositiveInt = function(x) {
  return x == Math.floor(x) && x > 0;
};

JmatR.isNegativeInt = function(x) {
  return x == Math.floor(x) && x < 0;
};

JmatR.isPositiveIntOrZero = function(x) {
  return x == Math.floor(x) && x >= 0;
};

JmatR.isNegativeIntOrZero = function(x) {
  return x == Math.floor(x) && x <= 0;
};

// x is odd integer
JmatR.isOdd = function(x) {
  return Math.abs(x % 2) == 1; //works for negative x too
};

// x is even integer
JmatR.isEven = function(x) {
  return x % 2 == 0; //works for negative x too
};

//isnanorinf isinfornan
JmatR.isInfOrNaN = function(x) {
  return x == Infinity || x == -Infinity || isNaN(x);
};

////////////////////////////////////////////////////////////////////////////////

JmatR.dist = function(a, b) {
  return Math.abs(a - b);
};

// works for non-integers too
JmatR.mod = function(a, b) {
  /*
  mod in terms of rem (%). The table below compares the two operators.
     x    :   -4 -3 -2 -1  0  1  2  3  4
  x mod  3:    2  0  1  2  0  1  2  0  1
  x mod -3:   -1  0 -2 -1  0 -2 -1  0 -2
  x rem  3:   -1  0 -2 -1  0  1  2  0  1
  x rem -3:   -1  0 -2 -1  0  1  2  0  1
  The sign of mod is that of b, while that of rem is that of a.
  */

  //without branching: return (b + (a % b)) % b OR alternatively and more generic: return a-floor(a/b)*b.

  // with branching (advantage: only one % division is performed)
  // the -b ==/!= a conditions are important, the else part does not do that one correct
  if(a <= 0) {
    if(b < 0 || a == 0 || -b == a) return a % b;
    else return b - ((-a) % (-b));
  } else {
    if(b < 0 && -b != a) return b - ((-a) % (-b));
    else return a % b;
  }
};

// Remainder. This is just the % operator. It is here only for reference. Compare with JmatR.mod, which is different.
JmatR.rem = function(a, b) {
  return a % b;
};

// to is not included in the range
JmatR.wrap = function(x, from, to) {
  if(from == to) return from;
  var m0 = Math.min(from, to);
  var m1 = Math.max(from, to);
  return m0 + JmatR.mod(x - m0, m1 - m0);
};

// to is included in the range
JmatR.clamp = function(x, from, to) {
  var m0 = Math.min(from, to);
  var m1 = Math.max(from, to);
  return Math.max(m0, Math.min(m1, x));
};

//Inspired by Wikipedia, Lanczos approximation, precision is around 15 decimal places
JmatR.gamma = function(z) {
  // Return immediately for some common values, to avoid filling the cache with those
  if(z == Infinity) return Infinity;
  if(JmatR.useFactorialLoop_(z - 1)) {
    return JmatR.factorial(z - 1); //that one uses memoization
  }
  if(z == 0.5) return JmatR.SQRTPI;

  // The internal function that doesn't do internal checks
  var gamma_ = function(z) {
    if(z <= 0 && z == Math.round(z)) return /*NaN*/ Infinity; //gamma not defined for negative integers. TODO: this should be "undirected" infinity

    // reflection formula
    if(z < 0.5) {
      return Math.PI / (Math.sin(Math.PI * z) * gamma_(1 - z));
    }

    var g = 7;
    var p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    z -= 1;
    var x = p[0];
    for(var i = 1; i < g + 2; i++) {
      x += p[i] / (z + i);
    }
    var t = z + g + 0.5;
    return Math.sqrt(Math.PI * 2) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
  };

  return gamma_(z);
};

JmatR.factorialmem_ = [1]; //memoization for factorial of small integers

JmatR.useFactorialLoop_ = function(x) {
  return JmatR.isPositiveIntOrZero(x) && x < 200;
}

JmatR.factorial = function(a) {
  if(!JmatR.useFactorialLoop_(a)) {
    return JmatR.gamma(a + 1);
  }

  if(JmatR.factorialmem_[a]) return JmatR.factorialmem_[a];

  var result = JmatR.factorialmem_[JmatR.factorialmem_.length - 1];
  for(var i = JmatR.factorialmem_.length; i <= a; i++) {
    result *= i;
    JmatR.factorialmem_[i] = result;
  }
  return result;
};

JmatR.firstPrimes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

// Initial set up shared by several of the prime test functions.
// Returns 0 if not prime, 1 if prime, NaN if problem, -1 if unknown by this function
JmatR.isPrimeInit_ = function(n) {
  if(n == Infinity || n != n) return NaN;
  if(n != Math.round(n)) return 0;
  if(n < 2) return 0;
  if(n > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense (decimal: 9007199254740992)
  for(var i = 0; i < JmatR.firstPrimes_.length; i++) {
    if(n == JmatR.firstPrimes_[i]) return 1;
    if(n % JmatR.firstPrimes_[i] == 0) return 0;
  }
  return -1;
}

//Returns 1 if prime, 0 if not prime, NaN if error. Naive slow algorithm. However, faster than miller rabin for n < 1500000
JmatR.isPrimeSlow_ = function(n) {
  // Ensures the number is odd and integer in supported range, tested against first known primes
  var init = JmatR.isPrimeInit_(n);
  if(init != -1) return init;

  var p = JmatR.firstPrimes_[JmatR.firstPrimes_.length - 1];
  var s = Math.ceil(Math.sqrt(n)) + 6;
  p = Math.floor(p / 6) * 6;
  while(p < s) {
    if(n % (p - 1) == 0 || n % (p + 1) == 0) return 0;
    p += 6;
  }
  return 1;
};

//Deterministic Miller-Rabin primality test
//Returns 1 if prime, 0 if not prime, NaN if error.
//Supposedly fast, but only faster than the naive method for n > 1500000
JmatR.isPrimeMillerRabin_ = function(n) {
  // Ensures the number is odd and integer in supported range, tested against first known primes
  var init = JmatR.isPrimeInit_(n);
  if(init != -1) return init;

  // Miller-Rabin test
  var base = undefined;
  if(n < 1373653) base = [2, 3];
  else if(n < 9080191) base = [31, 73];
  else if(n < 4759123141) base = [2, 7, 61];
  else if(n < 1122004669633) base = [2, 13, 23, 1662803];
  else if(n < 2152302898747) base = [2, 3, 5, 7, 11];
  else if(n < 3474749660383) base = [2, 3, 5, 7, 11, 13];
  else if(n < 341550071728321) base = [2, 3, 5, 7, 11, 13, 17];
  else if(n < 3770579582154547) base = [2, 2570940, 880937, 610386380, 4130785767];
  else base = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]; //valid up to >2^64

  var d = Math.floor(n / 2);
  var s = 1;
  while(!(d & 1)) {
    d = Math.floor(d / 2);
    ++s;
  }

  // returns (a * b) % c, taking overflow into account
  var modmul = function(a,b,c) {
    var x = 0;
    var y = a % c;
    while(b > 0){
      if(b & 1) x = (x + y) % c;
      y = (y * 2) % c;
      b= Math.floor(b / 2);
    }
    return x % c;
  };

  // returns (a to the n) % c, taking overflow into account
  var modpow = function(a, n, mod) {
    var result = 1;
    while(n > 0) {
      if(n & 1) result = modmul(result, a, mod);//(result * a) % mod;
      a = modmul(a, a, mod);//(a * a) % mod;
      n = Math.floor(n / 2);
    }
    return result;
  };

  var witness = function(n, s, d, a) {
    var x = modpow(a, d, n);
    while(s) {
      var y = modmul(x, x, n); //(x * x) % n;
      if(y == 1 && x != 1 && x != n - 1) return false;
      x = y;
      s--;
    }
    if(y != 1) return false;
    return true;
  };

  for(var i = 0; i < base.length; i++) {
    if(!witness(n, s, d, base[i])) return 0;
  }
  return 1;
};
/*
test in console for the above function:
function benchfun(n) {
  var ra = 0; var ta0 = new Date().getTime(); for(var i = 0; i < n; i++) ra += JmatR.isPrimeMillerRabin_(i); var ta1 = new Date().getTime();
  var rb = 0; var tb0 = new Date().getTime(); for(var i = 0; i < n; i++) rb += JmatR.isPrimeSlow_(i); var tb1 = new Date().getTime();
  var rc = 0; var tc0 = new Date().getTime(); for(var i = 0; i < n; i++) rc += JmatR.isPrime(i); var tc1 = new Date().getTime();
  console.log('fast: ' + (ta1 - ta0) + ' slow: ' + (tb1 - tb0) + ' both: ' + (tc1 - tc0) + ' test: ' + ra + ' = ' + rb + ' = ' + rc);
};
benchfun(100000);
--> it will report that slow if faster than miller rabin. That's because miller rabin is only faster for very large numbers. E.g. here you can see that miller rabin is fater:
JmatR.isPrimeSlow_(4444280714420857)
JmatR.isPrimeMillerRabin_(4444280714420857)


function testfun(n) {
  for(var i = 0; i < n; i++) {
    var a = JmatR.isPrimeMillerRabin_(i);
    var b = JmatR.isPrimeSlow_(i);
    var c = JmatR.isPrime(i);
    if(a != b || a != c) console.log('error: ' + i + ' ' + a + ' ' + b + ' ' + c);
  }
};
testfun(100000);

Nice primes to test:
3770579582154547 --> NOT prime, but above this boundary, last "base" for miller rabin test is used
9007199254740992 --> NOT prime, but highest integer number that JavaScript supports
9007199254740881: just small enough for JS! ==> overflow with sum, does not work
4444280714420857: one of the largest that does not overflow in isPrimeMillerRabin_
311111111111113: for third last base
344555666677777: for second last base
*/

//Returns 1 if prime, 0 if not prime, NaN if error.
JmatR.isPrime = function(n) {
  // below that, the "slow" method is faster. For higher values, Miller Rabin becomes more and more significantly faster.
  // However, for values above 0x010000000000000, a sum in the miller rabin overflows, so does not work either
  // ==> JmatR.isPrime(9007199254740881) is noticably slower than JmatR.isPrime(4444280714420857)
  return (n < 1500000 || n > 0x010000000000000) ? JmatR.isPrimeSlow_(n) : JmatR.isPrimeMillerRabin_(n);
};

//for factorize
JmatR.smallestPrimeFactor = function(x) {
  if(x == Infinity || x != x) return NaN;
  if(x != Math.round(x)) return NaN;
  if(x < 1) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense
  if(x == 1) return 1;
  for(var i = 0; i < JmatR.firstPrimes_.length; i++) {
    if(x == JmatR.firstPrimes_[i]) return x;
    if(x % JmatR.firstPrimes_[i] == 0) return JmatR.firstPrimes_[i];
  }
  var p = JmatR.firstPrimes_[JmatR.firstPrimes_.length - 1];
  var s = Math.ceil(Math.sqrt(x));
  p = Math.floor(p / 6) * 6;
  while(p < s) {
    if(x % (p - 1) == 0) return p - 1;
    if(x % (p + 1) == 0) return p + 1;
    p += 6;
  }
  return x;
};

//factorize: returns prime factors as array of real integers. x must be real positive integer.
JmatR.factorize = function(x) {
  var x = Math.round(x);
  var result = [];
  for(;;) {
    if(x < 1) break;
    var y = JmatR.smallestPrimeFactor(x);
    result.push(y);
    if(x == y) break;
    x = Math.round(x / y);
  }
  return result;
};


JmatR.primeCount = function(value) {
  var primesN = [ 0, 2, 3, 5, 7, 11, 13, 17 ];
  // Nth prime (1-indexed: n=1 gives 2)
  var p = function(n) {
    if(n < primesN.length) return primesN[n];
    var i = primesN[primesN.length - 1] + 2;
    var count = primesN.length - 1;
    for(;;) {
      if(JmatR.isPrime(i)) {
        primesN.push(i);
        count++;
        if(count == n) return i;
      }
      i += 2;
    }
  };

  var phiCache = {};
  // number of natural numbers smaller than m which are not divisible by
  // the first n primes
  var phi = function(m, n) {
    if(n == 0) return Math.floor(m);
    else if(n == 1) return Math.floor((m + 1) / 2);
    else {
      if(phiCache[m] && phiCache[m][n] != undefined) return phiCache[m][n];
      var result = phi(m, n - 1) - phi(Math.floor(m / p(n)), n - 1);
      if(!phiCache[m]) phiCache[m] = {};
      phiCache[m][n] = result;
      return result;
    }
  };

  var piCache = {};
  var pi = function(v) {
    if(v > 1000000000) return NaN; //it starts giving rounding errors or so somewhere before 1050000000
    if(v < 2) return 0;
    if(v < 3) return 1;
    if(v < 5) return 2;
    var n = Math.floor(v);
    if(piCache[n]) return piCache[n];
    var a = Math.floor(pi(Math.pow(v, 1.0 / 4.0)));
    var b = Math.floor(pi(Math.sqrt(v)));
    var c = Math.floor(pi(Math.pow(v, 1.0 / 3.0)));
    var sum = phi(n, a) + Math.floor((b + a - 2) * (b - a + 1) / 2);
    for(var i = a + 1; i <= b; i++) {
      var w = n / p(i); //NOT integer division!
      sum -= pi(w);
      if(i <= c) {
        var bi = pi(Math.sqrt(w));
        for(var j = i; j <= bi; j++) {
          sum -= pi(w / p(j)) - j + 1;
        }
      }
    }
    piCache[n] = sum;
    return sum;
  };

  return pi(value);
};

JmatR.nearestPrime = function(value) {
  var x = Math.round(value);
  if(x < 2) return 2;
  if(x == Infinity || x != x) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense

  if(JmatR.isPrime(x)) return x;
  var i = x % 2 == 0 ? 1 : 2;
  while(true) {
    if(JmatR.isPrime(x + i)) return x + i;
    if(JmatR.isPrime(x - i)) return x - i;
    i += 2; //TODO: make faster (e.g. using += 6 at the least)
  }
};

JmatR.nextPrime = function(value) {
  var x = Math.floor(value);
  if(x < 2) return 2; //the calculations below would give 3 instead of 2
  if(x == Infinity || x != x) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense

  var i = x % 2 == 0 ? 1 : 2;
  while(true) {
    if(JmatR.isPrime(x + i)) return x + i;
    i += 2; //TODO: make faster (e.g. using += 6 at the least)
  }
};

JmatR.previousPrime = function(value) {
  var x = Math.ceil(value);
  if(x <= 3) return 2; //avoid infinite loop
  if(x == Infinity || x != x) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense

  var i = x % 2 == 0 ? 1 : 2;
  while(true) {
    if(JmatR.isPrime(x - i)) return x - i;
    i += 2; //TODO: make faster (e.g. using += 6 at the least)
  }
};

JmatR.eulerTotient = function(value) {
  if(value <= 0) return NaN;
  var n = Math.floor(value);
  var f = JmatR.factorize(n);
  var prev = -1;
  var result = n;
  for(var i = 0; i < f.length; i++) {
    if(prev == f[i]) continue; //must be unique factors
    if(f[i] == 1) break;
    prev = f[i];
    result *= (1 - (1 / f[i]));
  }
  return jmatresult;
};

// The first integer binomials, allows fast calculation of those by just looking up in the array, e.g. binomial(5, 8) = JmatC.pascal_triangle_cache_[5][8]
// some rows area prefilled to start it off (just prefilling the first one would be sufficient normally, the rest is just for the shows)
JmatR.pascal_triangle_cache_ = [
    [1],
    [1, 1],
    [1, 2, 1],
    [1, 3, 3, 1],
    [1, 4, 6, 4, 1],
    [1, 5, 10, 10, 5, 1],
    [1, 6, 15, 20, 15, 6, 1],
    [1, 7, 21, 35, 35, 21, 7, 1],
];

// A helper function for integer binomial. Uses cached pascal triangle, so is guaranteed to be O(1) once the cache is filled up.
// Only works for integers, and only works for n < 180. After that, the double precision numbers no longer recognise every integer number.
JmatR.pascal_triangle = function(n, p) {
  if(n < 0 || p < 0 || n < p) return NaN;
  if(n > 180) return NaN; //triangle values get too big for integers in double precision floating point
  //fill up cache if needed
  var t = JmatR.pascal_triangle_cache_;
  while(t.length <= n) {
    var l = t.length; //the 'n' of the new row
    var l2 = l + 1; // number of elements of this new row
    t[l] = [];
    for(var i = 0; i < l2; i++) {
      t[l][i] = (i == 0 || i == l2 - 1) ? 1 : (t[l-1][i-1] + t[l-1][i]);
    }
  }
  return t[n][p];
}

//geatest common divisor
JmatR.gcd = function(x, y) {
 //Euclid's algorithm
 while(true) {
   if(y == 0) return x;
   var z = JmatR.mod(x, y);
   x = y;
   y = z;
 }
};

//least common multiple
JmatR.lcm = function(x, y) {
 return Math.abs(x * y) / JmatR.gcd(x, y);
};

// Decomposes fraction (aka rational approximation): returns two integers [numerator, denominator] such that n/d = a.
// Very slow, too slow for inner loop of running programs (integrate or complex plot)...
// max = max value for denominator
// E.g. JmatR.decompose(Math.PI, 100) gives [22, 7], because 22/7 approximates pi.
JmatR.decompose = function(x, max) {
  if(!max) max = 100000;
  var neg = (x < 0);
  if(neg) x = -x;
  var f = Math.floor(x);
  var y = x - f;

  if(y == 0) return [x, 1]; //otherwise the loop will run max times for nothing, very inefficient

  var result;

  var a = 0
  var b = 1;
  var c = 1;
  var d = 1;

  //mediant of two fractions a/c and b/d is defined as (a+b)/(c+d)
  while (b <= max && d <= max) {
    var mediant = (a + c) / (b + d);
    if(y == mediant) {
      if(b + d <= max) result = [a + c, b + d];
      else if(d > b) result = [c, d];
      else result = [a, b];
      break;
    } else if(y > mediant) {
      a = a + c;
      b = b + d;
    } else {
      c = a + c;
      d = b + d;
    }
  }
  if (!result) {
    if (b > max) result = [c, d]
    else result = [a, b];
  }

  result[0] += f * result[1];
  if(neg) result[0] = -result[0];

  return result;
};

// Hybrid between decompose and decomposeFast
JmatR.decomposeSemiFast = function(x, max) {
  var maxslow = 1000;
  if(max < maxslow) {
    return JmatR.decompose(x, max);
  } else {
    var a = JmatR.decompose(x, maxslow);
    var ax = a[0] / a[1];
    if(ax == x) return a;
    var b = JmatR.decomposeFast(x, maxslow);
    var bx = b[0] / b[1];
    return (Math.abs(x - ax) < Math.abs(x - bx)) ? a : b;
  }
};

// Decomposes fraction (aka rational approximation): returns two integers [numerator, denominator] such that n/d = a.
// max = max value for denominator
// A lot faster, but less nice than JmatR.decompose (e.g. returns 83333/10000 instead of 1/12), and with high preference for decimal denominators. TODO: other bases than base 10
JmatR.decomposeFast = function(x, max) {
  if(!max) max = 100000;
  var max1 = max - 1;

  if(x <= max1 && x >= -max1 && (x < -1.0 / max1 || x > 1.0 / max1)) {
    var neg = (x < 0);
    if(neg) x = -x;
    var z = Math.floor(x);
    var n = x - z;
    var d = max;
    n = Math.floor(n * d);
    var g = JmatR.gcd(n, d);
    d /= g;
    n /= g;
    n += z * d;
    if(neg) n = -n;
    // n = numerator, d = denominator
    return [n, d];
  }
  return [x, 1];
};

JmatR.near = function(x, y, precision) {
  // works also for infinities
  return x >= y - precision && x <= y + precision;
};

// Fractional part of x, x - floor(x). NOTE: this variant gives positive results for negative x
JmatR.frac = function(x) {
  return x - Math.floor(x);
};

// Fractional part of x, x - int(x). NOTE: this variant gives negative results for negative x
JmatR.fracn = function(x) {
  return x > 0 ? (x - Math.floor(x)) : -(-x - Math.floor(-x));
};

// Only the principal branch for real values above -1/e
JmatR.lambertw = function(x) {
  if(isNaN(x)) return NaN;
  if(x == Infinity || x == -Infinity) return Infinity;

  if(x >= -1.0 / Math.E && x <= 703) {
    //Newton's method. Only works up to 703
    var wj = x < 10 ? 0 : Math.log(x) - Math.log(Math.log(x)); // Without good starting value, it requires hundreds of iterations rather than just 30.
    var num = Math.max(30, x > 0 ? 10 + Math.floor(x) : 30);
    for(var i = 0; i < num; i++) {
      var ew = Math.exp(wj);
      wj = wj - ((wj * ew - x) / (ew + wj * ew));
    }
    return wj;
  } else if (x > 0) {
    //Since the above method works only up to 703, use some kind of binary search instead (it's a monotously increasing function at this point)
    // TODO: probably just use Halley's method here instead
    var step = 1;
    var lastDir = 0;
    var result = Math.log(x) - Math.log(Math.log(x)); // good starting value speeds up iterations. E.g. only 76 instead of 292 for 7e100.
    for(;;) {
      if(step == 0 || step * 0.5 == step || result + step == result) return result; //avoid infinite loop
      var v = result * Math.exp(result);
      if(JmatR.near(v, x, 1e-15)) return result;
      if(v > x) {
        result -= step;
        if(lastDir == -1) step *= 0.5;
        lastDir = 1;
      } else {
        result += step;
        if(lastDir == 1) step *= 0.5;
        lastDir = -1;
      }
    }
    return result;
  }
  return NaN;
};

// Tetration
// Returns experimental (not mathematically correct) results unless x is an integer or Infinity
JmatR.tetration = function(a, x) {
  // if(a == 1) return 1; // Not true, e.g. 1^Infinity is not 1. Probably things are not like this for tetration either.
  if(x == 0) return 1; //by definition
  if(x == 1) return a;
  if(x == 2) return Math.pow(a, a);
  if(a >= 2 && x > 5) return Infinity; // too big for double
  if(a == 0 && JmatR.isPositiveInt(x)) {
    // Chosen definition that considers 0^0 = 1, so 0^^n is 1 if n is even, 0 if n is odd.
    return JmatC.isEven(x) ? 1 : 0;
  }

  // Power tower (infinitely iterated exponentiation)
  if(x == Infinity && a > 0) {
    // converges if a >= 0.066 && a <= 1.44
    var l = Math.log(a);
    return JmatR.lambertw(-l) / (-l);
  }

  var runloop = function(a, b, num, l) {
    var result = b;
    var last;
    for(var i = 0; i < num; i++) {
      if(l) result = JmatR.logy(result, a);
      else result = Math.pow(a, result);
      if(isNaN(result)) return result;
      if(result == Infinity) return result; // Actually redundant, result == last already checks that too
      if(result == last) return result; // E.g. a=1.01 already converges to 1.0101015237405409 after just 7 iterations. 1.44 converges to 2.393811748202943 after 244 iterations, so numbers near that take a while to loop. 1.45 and higher converge to Infinity.
      last = result;
      if(i > 1000) return NaN; //avoid infinite loop
    }
    return result;
  }

  // This supports power tower as well
  // It supports it better than the lambertw formula above because the domain here is real numbers only, the lambertw formula fails on, for example, -1 and returns nan instead of -1
  if(JmatR.isPositiveInt(x)) {
    return runloop(a, a, x - 1, false);
  }

  // Everything above is true tetration for those cases where possible. What follows below is intermediate tetration research, to return "something"

  // Linear approximation for the extension to real heights
  // a^^x = x+1 for x > -1 && x <= 0
  // a^^x = log_a(a^^(x+1)) for x <= -1  ==> a^^-1.5 = log_a(x+2), a^^-2.5 = log_a(log_a(x+3)), etc...
  // a^^x = a^(a^^(x-1)) for x > 0  ==> a^^0.5 = a^x, a^^1.5 = a^(a^(x-1)), a^^2.5 = a^(a^(a^(x-2))), etc...
  // examples: e^^(0.5*pi) ~= 5.868
  //           0.5^^(-4.3) ~= 4.033
  if(x > -1 && x <= 0) return 1 + x;
  if(x > 0) {
    var b = x - Math.floor(x); //always in range 0-1
    return runloop(a, b, Math.ceil(x), false);
  }
  if(x <= -1) {
    var b = x - Math.floor(x); //always in range 0-1
    return runloop(a, b, -Math.ceil(x), true);
  }

  return NaN;
};


//arbitrary log: log_y(x)
//warning: base y is second argument
JmatR.logy = function(x, y) {
  return Math.log(x) / Math.log(y);
};

// dawson function D+(x) aka F(x). erfi(x) = 2/sqrt(pi) * exp(x*x) * D(x)
// so while erfi overflows for large args, this one doesn't due to no exp(x*x)
JmatR.dawson = function(x) {
  var x2 = x * x;
  var x4 = x2 * x2;
  var x6 = x4 * x2;
  var x8 = x6 * x2;
  var x10 = x8 * x2;
  var x12 = x10 * x2;

  var p1 = 0.1049934947, p2 = 0.0424060604, p3 = 0.0072644182, p4 = 0.0005064034, p5 = 0.0001789971;
  var q1 = 0.7715471019, q2 = 0.2909738639, q3 = 0.0694555761, q4 = 0.0140005442, q5 = 0.0008327945;

  var p = 1 + p1 * x2 + p2 * x4 + p3 * x6 + p4 * x8 + p5 * x10;
  var q = 1 + q1 * x2 + q2 * x4 + q3 * x6 + q4 * x8 + q5 * x10 + 2 * p5 * x12;

  return p / q * x;
};

//erfi, the imaginary error function, on real argument: erfi(z) = -i erf(iz)
JmatR.erfi = function(x) {
  var neg = false;
  if(x < 0) {
    x = -x;
    neg = true;
  }
  var result = 0;
  var ps = 1.0 / JmatR.SQRTPI;

  if(x <= 0.5) {
    // only gives good approximation for x < 0.5 or so
    var x3 = x * x * x;
    var x5 = x3 * x * x;
    var x7 = x5 * x * x;
    var t = 2*x + 2/3*x3 + 1/5*x5 + 1/21*x7;
    result = ps * t;
  } else if (x >= 5) {
    // only gives good approximation for x > 5 or so
    var xi = 1 / x;
    var xi3 = xi * xi * xi;
    var xi5 = xi3 * xi * xi;
    var xi7 = xi5 * xi * xi;
    var e = Math.exp(x * x);
    var t = xi + 1/2*xi3 + 3/4*xi5 + 15/8*xi7;
    result = ps * e * t;
  } else {
    result = 2 / JmatR.SQRTPI * Math.exp(x * x) * JmatR.dawson(x);
  }

  if(neg) result = -result;
  return result;
};

JmatR.erf = function(x) {
  var neg = x < 0;
  if(neg) x = -x;

  if (x == 0) return 0;
  var t = 1 / (1 + 0.3275911 * x);
  var p = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  var result = 1.0 - p * Math.exp(-(x*x));

  if(neg) result = -result;
  return result;
};

JmatR.erfc = function(x) {
  var neg = x < 0;
  if(neg) x = -x;
  var result;

  if(x <= 0.5) {
    var x2 = x * x;
    var x3 = x * x2;
    var x5 = x3 * x2;
    var x7 = x5 * x2;
    result = 1 - 2 / JmatR.SQRTPI * (x - x3 / 3 + x5 / 10 + x7 / 42);
    //result = Math.exp(-x*x) / 6 + Math.exp(-0.75 * x * x) / 2;
  } else if (x >= 5) {
    // asymptotic expansion for large real x
    var x2 = x * x;
    var x4 = x2 * x2;
    var x6 = x4 * x2;
    var x8 = x6 * x2;
    result = Math.exp(-(x*x)) / (x * JmatR.SQRTPI) * (1 - 1/2.0/x2 + 3/4.0/x4 - 15/8.0/x6 + 105/16.0/x8);
  } else {
    var t = 1 / (1 + 0.3275911 * x);
    var p = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    result = p * Math.exp(-(x*x));
  }

  if(neg) result = 2 - result;
  return result;
};

//Minkowski's question mark function, from Wikipedia
JmatR.minkowski = function(x) {
  if(x != x) return NaN;
  var p = Math.floor(x);
  var q = 1, r = p + 1, s = 1, m, n;
  var d = 1.0, y = p;
  if(x < p || (p < 0) != (r <= 0)) return x; //out of range ?(x) =~ x
  while(true) //invariants: q*r-p*s==1 && p/q <= x && x < r/s
  {
    d /= 2; if(y + d == y) break; //reached max possible precision
    m = p + r; if((m < 0) != (p < 0)) break; //sum overflowed
    n = q + s; if(n < 0) break; //sum overflowed

    if(x < m / n) r = m, s = n;
    else y += d, p = m, q = n;
  }
  return y + d; //final round-off
};

// decimal to degrees/minutes/second (e.g. 1.5 degrees becomes 1.30 (1 degree and 30 minutes))
JmatR.dms = function(a) {
  var neg = a < 0;
  if(neg) a = -a;

  var deg = Math.floor(a);
  var mins = Math.floor((a * 60 - deg * 60));
  var sec = Math.floor(a * 3600 - deg * 3600 - mins * 60);

  var result = deg + mins / 100.0 + sec / 10000.0;

  if(neg) result = -result;
  return result;
};

// degrees/minutes/second to decimal degrees (e.g. 1.30 (1 deg 30 minutes) becomes 1.5 degrees)
JmatR.dd = function(a) {
  var neg = a < 0;
  if(neg) a = -a;

  var deg = Math.floor(a);
  var mins = Math.floor((a * 100 - deg * 100));
  var sec = Math.floor(a * 10000 - deg * 10000 - mins * 100);

  var result = deg + mins / 60.0 + sec / 3600.0

  if(neg) result = -result;
  return result;
};

////////////////////////////////////////////////////////////////////////////////

JmatR.isLeapYear = function(y) {
  return y % 400 == 0 || (y % 4 == 0 && y % 100 != 0);
};

JmatR.montharray_ = [-1 /*there is no month 0*/, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; //february must be adjusted for leap year

JmatR.monthLength = function(month, leap) {
  return (leap && month == 2) ? 29 : JmatR.montharray_[month];
};


//number of days since the first day of year 0. 1 january of the year 0 is 0. 2 january is 1, etc...
//only for Gregorian calendar, does not take Julian calendar into account
JmatR.numDaysSince0 = function(year, month, day) {
  //calculate number of leap years before this year (year 0 is considered leap)
  var num400 = Math.floor((year - 1) / 400) + 1;
  var num100 = Math.floor((year - 1) / 100) + 1 - num400;
  var num4 = Math.floor((year - 1) / 4) + 1 - num100 - num400;
  var numleap = num4 + num400;
  var yeardays = numleap * 366 + (year - numleap) * 365; //days of years before this year

  var monthdays = 0; //days of years before this month
  var leap = JmatR.isLeapYear(year);
  for (var i = 1; i < month; i++) monthdays += JmatR.monthLength(i, leap);

  return yeardays + monthdays + day - 1; //-1 because day 1 is in fact zero
};

//converts number of days since the year 0, to year/month/day
//returns array [year, month, day]
//only for Gregorian calendar, does not take Julian calendar into account
JmatR.daysSince0ToDate = function(days) {
  //every 400 years there are 97 leap years. So a year is 365.2425 days in average.
  var year = Math.floor(days / 365.2425);
  var leap = JmatR.isLeapYear(year);

  days -= JmatR.numDaysSince0(year, 1, 1);

  var month = 0;
  for (var i = 1; i <= 12; i++) {
    month++;
    var num = JmatR.monthLength(i, leap);
    if (days >= num /*>= because of the +1 done at the end to turn zero based into one based indexing*/) {
      days -= num;
    } else {
      break;
    }
  }

  return [year, month, days + 1 /*because month day starts at 1 instead of 0*/];
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Complex
////////////////////////////////////////////////////////////////////////////////

/*
Constructor.
Class representing a complex value with real and imaginary part.

This is the actual object used as complex number. In addition, most of the
functions are implemented as static functions in here.

The only sad thing is that Javascript doesn't support operator overloading
and nice expressions like a + b have to become a.add(b) instead.
*/
JmatC = function(re, im) {
  // Does not do any checks, to be "fast"
  this.re = re;
  this.im = im;
};

// TODO: define a bit better which combinations of Infinity/Nan/... in re and im mean what (E.g. re and im both Infinity means "undirected infinity", already used by gamma function but by nothing else)

// Create a new JmatC complex value. Copies JmatC if a JmatC is given as first argument
// with 0 arguments, creates zero value. With a and b numbers, creates complex number from it. With a JmatC object, copies it.
// the first parameter must be given and be number or JmatC. The second parameter is optional.
JmatC.make = function(a, b) {
  if(a.re == undefined) return new JmatC(a, b == undefined ? 0 : b);
  return new JmatC(a.re, a.im); // Copy value object
};
// Shortcut function because this is so common
jmatc = JmatC.make;

// Create a new JmatC value, real
JmatC.newr = function(re) {
  return new JmatC(re, 0);
};

// Create a new JmatC value, imaginary
JmatC.newi = function(im) {
  return new JmatC(0, im);
};

// Create a new JmatC value, polar
JmatC.polar = function(r, a) {
  return new JmatC(r * Math.cos(a), r * Math.sin(a));
};

// Casts the given number type to JmatC. If the given type is already of type JmatC, does not copy it but returns the input.
// TODO: also support strings of the form '5+6i', and be able to parse them
JmatC.cast = function(v) {
  if(v && v.re != undefined) return v;
  if(v == undefined) return jmatc(0);
  return jmatc(v);
};

// Only use these as constants, never modify these, never return them!
JmatC.ZERO = jmatc(0);
JmatC.ONE = jmatc(1);
JmatC.TWO = jmatc(2);
JmatC.I = JmatC.newi(1);
JmatC.PI = jmatc(Math.PI);
JmatC.E = jmatc(Math.E);
JmatC.SQRT2 = jmatc(Math.sqrt(2));
JmatC.SQRTPI = jmatc(Math.sqrt(Math.PI));
JmatC.INVSQRT2PI = jmatc(1 / Math.sqrt(2 * Math.PI)); //0.3989422804014327
JmatC.EM = jmatc(JmatR.EM); // Euler-Mascheroni constant
JmatC.APERY = jmatc(JmatR.APERY); // Apery's constant, zeta(3)

JmatC.real = function(z) {
  return jmatc(z.re);
};
JmatC.prototype.real = function() {
  return jmatc(this.re);
};

JmatC.imag = function(z) {
  return jmatc(z.im);
};
JmatC.prototype.imag = function() {
  return jmatc(this.im);
};

//Basic operators

JmatC.add = function(x, y) {
  return new JmatC(x.re + y.re, x.im + y.im);
};
JmatC.prototype.add = function(y) {
  return new JmatC(this.re + y.re, this.im + y.im);
};

JmatC.sub = function(x, y) {
  return new JmatC(x.re - y.re, x.im - y.im);
};
JmatC.prototype.sub = function(y) {
  return new JmatC(this.re - y.re, this.im - y.im);
};

JmatC.mul = function(x, y) {
  if(x.im == 0 && y.im == 0) {
    return new JmatC(x.re * y.re, 0);
  } else {
    var re = x.re * y.re - x.im * y.im;
    var im = x.im * y.re + x.re * y.im;
    return new JmatC(re, im);
  }
};
JmatC.prototype.mul = function(y) {
  return JmatC.mul(this, y);
};

JmatC.div = function(x, y) {
  if(x.im == 0 && y.im == 0) {
    return new JmatC(x.re / y.re, 0);
  } else {
    if(JmatC.isInf(x) && !JmatC.isInfOrNaN(y)) {
      // Result should be some infinity (because it's infinity divided through finite value), but the formula below would give a NaN somewhere.
      // 4 possible rotations of the infinity, based on quadrant of y (TODO: THIS IS IGNORED NOW!!)
      return x;
    }
    var d = y.re * y.re + y.im * y.im;
    if(d == Infinity || d == -Infinity) {
      // the calculations below would give Infinity/Infinity = NaN even though result should be 0.
      if(!JmatC.isInfOrNaN(x)) return jmatc(0);
    }
    if(d == 0 && !JmatC.isInfOrNaN(x) && (x.re != 0 || x.im != 0)) {
      // the calculations below would give 0/0 = NaN even though result should be some infinity.
      return new JmatC(x.re == 0 ? 0 : (x.re < 0 ? -Infinity : Infinity), x.im == 0 ? 0 : (x.im < 0 ? -Infinity : Infinity));
    }
    var re = (x.re * y.re + x.im * y.im) / d;
    var im = (x.im * y.re - x.re * y.im) / d;
    return new JmatC(re, im);
  }
};
JmatC.prototype.div = function(y) {
  return JmatC.div(this, y);
};

JmatC.addr = function(z, a) {
  return jmatc(z.re + a, z.im);
};
JmatC.prototype.addr = function(a) {
  return jmatc(this.re + a, this.im);
};

JmatC.subr = function(z, a) {
  return jmatc(z.re - a, z.im);
};
JmatC.prototype.subr = function(a) {
  return jmatc(this.re - a, this.im);
};
JmatC.rsub = function(a, z) {
  return jmatc(a - z.re, -z.im);
};
JmatC.prototype.rsub = function(a) {
  return jmatc(a - this.re, -this.im);
};

JmatC.mulr = function(z, a) {
  return jmatc(z.re * a, z.im * a);
};
JmatC.prototype.mulr = function(a) {
  return jmatc(this.re * a, this.im * a);
};

JmatC.divr = function(z, a) {
  return jmatc(z.re / a, z.im / a);
};
JmatC.prototype.divr = function(a) {
  return jmatc(this.re / a, this.im / a);
};
JmatC.rdiv = function(a, z) {
  return JmatC.div(jmatc(a), z);
};

//rotate complex number z by a radians. That is, change its argument. a is real (JS number).
JmatC.rotate = function(z, a) {
  if(a == 0) return z;
  return JmatC.polar(z.absr(), z.argr() + a);
};

//rotate complex number z by 2pi/n radians. This results in giving the next solution of the nth root.
JmatC.nextroot = function(z, n) {
  var result = JmatC.rotate(z, Math.PI * 2 / n);
  if(JmatR.near(result.im, 0, 1e-14)) result.im = 0;
  return result;
};

// mod operation, result has the sign of the divisor (unlike % operator in JS, Java and C99), so it's like wrapping x in range 0..y.
// works on real or complex numbers too, e.g. (6+4i) mod (3+5i) gives (-2+2i)
JmatC.mod = function(x, y) {
  if(x.im != 0 || y.im != 0) return x.sub(JmatC.floor(x.div(y)).mul(y));
  return jmatc(JmatR.mod(x.re, y.re));
};

// remainder operation, like the % operator in JS, Java and C99.
JmatC.rem = function(x, y) {
  if(x.im != 0 || y.im != 0) return x.sub(JmatC.trunc(x.div(y)).mul(y));
  return jmatc(x.re % y.re);
};

JmatC.wrap = function(x, from, to) {
  return new JmatC(JmatR.wrap(x.re, from.re, to.re), JmatR.wrap(x.im, from.im, to.im));
};

JmatC.clamp = function(x, from, to) {
  return new JmatC(JmatR.clamp(x.re, from.re, to.re), JmatR.clamp(x.im, from.im, to.im));
};

JmatC.bitnot = function(x) {
  var result = jmatc(0);
  result.re = ~x.re;
  //imaginary part not affected on purpose: otherwise it appears when bit-inverting real number, which is in 99.9% of the cases not wanted
  //result.im = ~x.im;
  return result;
};

JmatC.bitand = function(x, y) {
  var result = jmatc(0);
  result.re = x.re & y.re;
  result.im = x.im & y.im;
  return result;
};

JmatC.bitor = function(x, y) {
  var result = jmatc(0);
  result.re = x.re | y.re;
  result.im = x.im | y.im;
  return result;
};

JmatC.bitxor = function(x, y) {
  var result = jmatc(0);
  result.re = x.re ^ y.re;
  result.im = x.im ^ y.im;
  return result;
};

JmatC.lshift = function(x, y) {
  var result = jmatc(0);
  result.re = x.re << y.re;
  result.im = x.im << y.im;
  return result;
};

JmatC.rshift = function(x, y) {
  var result = jmatc(0);
  result.re = x.re >> y.re;
  result.im = x.im >> y.im;
  return result;
};

JmatC.neg = function(x) {
  return jmatc(-x.re, -x.im);
};
JmatC.prototype.neg = function() {
  return jmatc(-this.re, -this.im);
};

// Returns 0 if z is 0, 1 if z is positive, -1 if z is negative. For complex z, returns z / abs(z)
JmatC.sign = function(z) {
  if (z.im == 0) {
    if(z.re == 0) return jmatc(0);
    else if(z.re < 0) return jmatc(-1);
    return jmatc(1);
  }

  return z.div(JmatC.abs(z));
};

// Returns 0 if z is 0, 1 if z is positive, -1 if z is negative. For complex z, returns sign of z.im if z.re == 0, sign of z.re otherwise (that is, the function returns sqrt(z*z) / z, except for z=0)
JmatC.csgn = function(z) {
  if (JmatR.near(z.re, 0, 1e-15)) { //avoid numeric imprecisions for e.g. the values of e.g. acosh
    if(z.im == 0) return jmatc(0);
    else if(z.im < 0) return jmatc(-1);
    return jmatc(1);
  } else {
    if(z.re == 0) return jmatc(0);
    else if(z.re < 0) return jmatc(-1);
    return jmatc(1);
  }
};

JmatC.conj = function(x) {
  return jmatc(x.re, -x.im);
};
JmatC.prototype.conj = function() {
  return jmatc(this.re, -this.im);
};

JmatC.eq = function(x, y) {
  if(!x || !y) return x == y;
  return (x.re == y.re && x.im == y.im);
};
JmatC.prototype.eq = function(y) {
  return y && this.re == y.re && this.im == y.im;
};

JmatC.eqr = function(x, y) {
  return (x.re == y && x.im == 0);
};
JmatC.prototype.eqr = function(y) {
  return (this.re == y && this.im == 0);
};

JmatC.powr = function(z, a) {
  return JmatC.pow(z, jmatc(a));
};
JmatC.prototype.powr = function(a) {
  return JmatC.pow(this, jmatc(a));
};

JmatC.inv = function(z) {
  return JmatC.ONE.div(z);
};
JmatC.prototype.inv = function() {
  return JmatC.ONE.div(this);
};

//increment
JmatC.inc = function(z) {
  return new JmatC(z.re + 1, z.im);
};
JmatC.prototype.inc = function() {
  return new JmatC(this.re + 1, this.im);
};

//decrement
JmatC.dec = function(z) {
  return new JmatC(z.re - 1, z.im);
};
JmatC.prototype.dec = function(z) {
  return new JmatC(this.re - 1, this.im);
};

// absolute value squared, returned as real
JmatC.abssqr = function(x) {
  return x.re * x.re + x.im * x.im;
};
JmatC.prototype.abssqr = function() {
  return this.re * this.re + this.im * this.im;
};

// absolute value, aka modulus of complex number
JmatC.absr = function(x) {
  if(x.im == 0) return Math.abs(x.re);
  if(x.re == 0) return Math.abs(x.im);

  if(x.re == Infinity || x.re == -Infinity || x.im == Infinity || x.im == -Infinity) {
    return Infinity;
  }

  // Numerically more stable version of "Math.sqrt(x.re * x.re + x.im * x.im);"
  var sqr = function(x) {
    return x * x;
  };
  var absre = Math.abs(x.re);
  var absim = Math.abs(x.im);
  if(absre > absim) return absre * Math.sqrt(1 + sqr(x.im / x.re));
  else if(x.im == 0) return 0;
  else return absim * Math.sqrt(1 + sqr(x.re / x.im));
};
JmatC.prototype.absr = function() {
  return JmatC.absr(this);
};

JmatC.abs = function(x) {
  var result = jmatc(0);
  if(x.im == 0) result.re = Math.abs(x.re);
  else result.re = JmatC.absr(x);
  return result;
};
JmatC.prototype.abs = function() {
  return JmatC.abs(this);
};

// returns the complex argument in range -PI to +PI
JmatC.argr = function(x) {
  if(x.im == 0) return x.re < 0 ? Math.PI : 0;
  return Math.atan2(x.im, x.re);
};
JmatC.prototype.argr = function() {
  return JmatC.argr(this);
};

JmatC.arg = function(z) {
  return jmatc(JmatC.argr(z));
};
JmatC.prototype.arg = function() {
  return JmatC.arg(this);
};

//returns result in range 0-1 rather than -PI to PI. Useful for graphical representations, not for math. 0 matches 0 degrees, 0.5 matches 180 degrees, 0.999 matches around 359 degrees.
JmatC.argr1 = function(z) {
  var result = JmatC.argr(z);
  if(result < 0) result += 2 * Math.PI;
  result /= (2 * Math.PI);
  if(result < 0) result = 0;
  if(result > 1) result = 1;
  return result;
};


////////////////////////////////////////////////////////////////////////////////
// Categories
////////////////////////////////////////////////////////////////////////////////

JmatC.isReal = function(z) {
  return z.im == 0;
};

JmatC.isImaginary = function(z) {
  return z.re == 0;
};

JmatC.isInt = function(z) {
  return z.im == 0 && JmatR.isInt(z.re);
};

// Gaussian integer
JmatC.isGaussian = function(z) {
  return JmatR.isInt(z.re) && JmatR.isInt(z.im);
};

JmatC.isNaN = function(z) {
  return !z || isNaN(z.re) || isNaN(z.im);
};

//is infinite
JmatC.isInf = function(z) {
  return (z.re * z.re + z.im * z.im) == Infinity;
};

//isnanorinf isinfornan
JmatC.isInfOrNaN = function(z) {
  return !z || JmatR.isInfOrNaN(z.re) || JmatR.isInfOrNaN(z.im);
};

//real and strictly positive
JmatC.isPositive = function(z) {
  return z.re > 0 && z.im == 0;
};

//real and strictly negative
JmatC.isNegative = function(z) {
  return z.re < 0 && z.im == 0;
};

JmatC.isPositiveOrZero = function(z) {
  return z.re >= 0 && z.im == 0;
};

JmatC.isNegativeOrZero = function(z) {
  return z.re <= 0 && z.im == 0;
};

//strictly positive
JmatC.isPositiveInt = function(z) {
  return JmatC.isInt(z) && z.re > 0;
};

//strictly negative
JmatC.isNegativeInt = function(z) {
  return JmatC.isInt(z) && z.re < 0;
};

JmatC.isPositiveIntOrZero = function(z) {
  return JmatC.isInt(z) && z.re >= 0;
};

JmatC.isNegativeIntOrZero = function(z) {
  return JmatC.isInt(z) && z.re <= 0;
};

// z is odd integer
JmatC.isOdd = function(z) {
  return JmatC.isInt(z) && Math.abs(z.re % 2) == 1;
};

// z is even integer
JmatC.isEven = function(z) {
  return JmatC.isInt(z) && z.re % 2 == 0;
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pow = function(x, y) {
  if(JmatC.isReal(x) && JmatC.isReal(y) && (x.re >= 0 || y.re == Infinity || y.re == -Infinity || JmatR.isInt(y.re))) {
    //if(x.re == 0 && y.re == 0) return jmatc(NaN); // JS's pow returns 1 for 0^0
    // It is chosen to return 1 for 0^0, not NaN. NaN is mathematically more correct, however 0^0 is correct in many practical applications.
    return jmatc(Math.pow(x.re, y.re));
  } else {
    // This is just one branch. In fact it returns a complex result for -3 ^ (1/3),
    // the cube root of -3. To get the real result, use absolute value (and then negate) on it.
    // This is correct: the principal result of the cube root for this is a complex number.
    // Note: This returns incorrect values for a negative real to the power of Infinity: the result should be -Infinity for < -1, 0 for > -1, NaN for -1, but it always gives NaN. However, the "if" part above already handles that.
    var r = JmatC.absr(x);
    var t = JmatC.argr(x);
    var u = Math.pow(r, y.re) * Math.exp(-y.im * t);
    if(isNaN(u)) {
      u = Math.pow(1, y.re / r) * Math.exp(-y.im * t / r);
      if(u < 0) u = -Infinity;
      else if(u > 0) u = Infinity;
      else u = NaN;
    }
    var v = y.im * Math.log(r) + y.re * t;
    return jmatc(u * Math.cos(v), u * Math.sin(v));
  }
};
JmatC.prototype.pow = function(y) {
  return JmatC.pow(this, y);
};

JmatC.sin = function(z) {
  if(z.im == 0) return jmatc(Math.sin(z.re));

  var iz = jmatc(-z.im, z.re);
  var eiz = JmatC.exp(iz);
  var ieiz = JmatC.inv(eiz);
  return eiz.sub(ieiz).div(jmatc(0, 2));
};

//unnormalized sinc: sin(x) / x, but also defined for x = 0
JmatC.sinc = function(z) {
  if(z.eqr(0)) return jmatc(1);
  return JmatC.sin(z).div(z);
};

JmatC.cos = function(z) {
  if(z.im == 0) return jmatc(Math.cos(z.re));

  var iz = jmatc(-z.im, z.re);
  var eiz = JmatC.exp(iz);
  var ieiz = JmatC.inv(eiz);
  return eiz.add(ieiz).mulr(0.5);
};

JmatC.tan = function(z) {
  if(z.im == 0) return jmatc(Math.tan(z.re));

  var iz = jmatc(-z.im, z.re);
  var eiz = JmatC.exp(iz);
  var ieiz = JmatC.inv(eiz);
  return (eiz.sub(ieiz).div(jmatc(0, 2))).div(eiz.add(ieiz).mulr(0.5)); // JmatC.sin(z).div(JmatC.cos(z));
};

JmatC.asin = function(z) {
  if(z.im == 0 && z.re >= -1 && z.re <= 1) return jmatc(Math.asin(z.re));

  var s = JmatC.sqrt(JmatC.ONE.sub(z.mul(z)));
  var l = JmatC.log(jmatc(-z.im, z.re).add(s));
  return jmatc(l.im, -l.re);
};

JmatC.acos = function(z) {
  if(z.im == 0 && z.re >= -1 && z.re <= 1) return jmatc(Math.acos(z.re));

  //i * ln(x - i * sqrt(1-x^2))
  var s = JmatC.sqrt(JmatC.ONE.sub(z.mul(z))).mul(JmatC.I);
  var l = JmatC.log(z.add(s));
  return jmatc(l.im, -l.re);
};

JmatC.atan = function(z) {
  if(z.im == 0) return jmatc(Math.atan(z.re));

  var iz = jmatc(-z.im, z.re);
  var b = JmatC.ONE.sub(iz).div(iz.inc());
  var l = JmatC.log(b);
  return jmatc(-0.5 * l.im, 0.5 * l.re);
};

JmatC.atan2 = function(x, y) {
  if(!JmatC.isReal(x) || !JmatC.isReal(y)) {
    if(y.eqr(0)) return jmatc(Math.PI / 2);

    // For complex values, an alternate form of the defintion can be used:
    // 2 * atan(y / (sqrt(x^2+y^2)+x))
    return JmatC.atan(JmatC.sqrt(x.mul(x).add(y.mul(y))).sub(x).div(y)).mulr(2);
  } else {
    var result = jmatc(0);
    result.re = Math.atan2(x.re, y.re);
    return result;
  }
};

JmatC.sinh = function(z) {
  var e = JmatC.exp(z);
  var ei = JmatC.inv(e);
  return e.sub(ei).divr(2);
};

JmatC.cosh = function(z) {
  var e = JmatC.exp(z);
  var ei = JmatC.inv(e);
  return e.add(ei).divr(2);
};

JmatC.tanh = function(z) {
  var e = JmatC.exp(z);
  var ei = JmatC.inv(e);
  return e.sub(ei).div(e.add(ei));
};

JmatC.asinh = function(z) {
  return JmatC.log(z.add(JmatC.sqrt(z.mul(z).addr(1))));
};

JmatC.acosh = function(z) {
  // ln(x + sqrt(z-1)*sqrt(z+1))
  return JmatC.log(z.add(JmatC.sqrt(z.subr(1)).mul(JmatC.sqrt(z.addr(1)))));
};

JmatC.atanh = function(z) {
  // 0.5 * (ln(1+z) - ln(1-z))
  return JmatC.log(z.inc().div(z.dec())).mulr(0.5);
};

// This is NOT the logsine function (the intergral). It's simply ln(sin(z))
//ln(sin(z)), with good approximation for large |Im(z)|. The thing is, for large imaginary values, sin(z) becomes huge, because it involves an exponential of the imaginary parts
// For large imaginary part (or very small below 0), log(sin(x)) fails while this function is then very accurate
JmatC.logsin = function(z) {
  if(z.im > -10 && z.im < 10) return JmatC.log(JmatC.sin(z));

  var ln2i = jmatc(0.69314718056, 1.570796326795); // ln(2i)
  // This approximation is using a formula e^ix/2i or -e^(-ix)/2i, instead of the full (e^ix - e^(-ix) / 2i) = sin(x). This requires the real part to be exactly in range -pi/2, 3pi/2. So wrap, since it's periodic.
  var p = jmatc(JmatR.wrap(z.re, -Math.PI / 2, 3 * Math.PI / 2), z.im);
  if(z.im > 0) return JmatC.newi(JmatC.PI).sub(JmatC.I.mul(p)).sub(ln2i);
  else return JmatC.I.mul(p).sub(ln2i);
};

// See description of JmatC.logsin
JmatC.logcos = function(z) {
  return JmatC.logsin(z.rsub(Math.PI / 2));
};

JmatC.floor = function(x) {
  var result = jmatc(0);
  result.re = Math.floor(x.re);
  result.im = Math.floor(x.im);
  return result;
};

JmatC.ceil = function(x) {
  var result = jmatc(0);
  result.re = Math.ceil(x.re);
  result.im = Math.ceil(x.im);
  return result;
};

JmatC.round = function(x) {
  var result = jmatc(0);
  result.re = Math.round(x.re);
  result.im = Math.round(x.im);
  return result;
};

// truncate towards 0
JmatC.trunc = function(x) {
  var result = jmatc(0);
  result.re = x.re < 0 ? Math.ceil(x.re) : Math.floor(x.re);
  result.im = x.im < 0 ? Math.ceil(x.im) : Math.floor(x.im);
  return result;
};

// Fractional part of x, x - floor(x). NOTE: this variant gives positive results for negative x
JmatC.frac = function(x) {
  return jmatc(JmatR.frac(x.re), JmatR.frac(x.im));
};

// Fractional part of x, x - int(x). NOTE: this variant gives negative results for negative x
JmatC.fracn = function(x) {
  return jmatc(JmatR.fracn(x.re), JmatR.fracn(x.im));
};

JmatC.exp = function(x) {
  if(x.im == 0) {
    return jmatc(Math.exp(x.re));
  } else {
    var ea = Math.exp(x.re);
    return new JmatC(ea * Math.cos(x.im), ea * Math.sin(x.im));
  }
};

//exp(x) - 1, with better precision for x around 0
JmatC.expm1 = function(x) {
  if(JmatC.abssqr(x) < 1e-5) return x.add(x.mul(x).divr(2)).add(x.mul(x).mul(x).divr(6));
  else return JmatC.exp(x).subr(1);
};

//natural log (base e, ln)
JmatC.log = function(x) {
  if(x.eqr(-Infinity)) return jmatc(Infinity);

  if(JmatC.isReal(x) && x.re >= 0) {
    return jmatc(Math.log(x.re));
  }

  return jmatc(Math.log(JmatC.absr(x)), JmatC.argr(x));
};

//ln(x + 1), with better precision for x around 0
JmatC.log1p = function(x) {
  if(JmatC.abssqr(x) < 1e-8) return x.mulr(-0.5).addr(1).mul(x);
  else return JmatC.log(x.addr(1));
};

//arbitrary log: log_y(x)
//warning: base y is second argument
JmatC.logy = function(x, y) {
  return JmatC.log(x).div(JmatC.log(y));
};

JmatC.sqrt = function(x) {
  if(JmatC.isReal(x)) {
    var result = jmatc(0);
    if(x.re >= 0 || x.re != x.re) result.re = Math.sqrt(x.re);
    else result.im = Math.sqrt(-x.re);
    return result;
  } else return x.pow(jmatc(0.5));
};

// Because JS number toFixed appends zeros
JmatC.formatFloat_ = function(value, precision) {
  var power = Math.pow(10, precision || 0);
  return String(Math.round(value * power) / power);
};

//debugstring
JmatC.toString = function(value, opt_precision) {
  if(!value) return value == 0 ? 'invalid0' : ('' + value);
  var re = (opt_precision ? JmatC.formatFloat_(value.re, opt_precision) : ('' + value.re));
  var im = (opt_precision ? JmatC.formatFloat_(value.im, opt_precision) : ('' + value.im));

  if(value.im == 0 || im == '0') return '' + re;
  if(value.re == 0) return '' + im + 'i';
  if(value.im < 0) return '' + re + im + 'i';
  return '' + re + '+' + im + 'i';
};
JmatC.prototype.toString = function(opt_precision) {
  return JmatC.toString(this, opt_precision);
};

JmatC.toInt = function(value) {
  return Math.round(value.re);
};

// normalizes even if re or im are infinite, e.g. (Infinity, -Infinity) becomes (1, -1), (0, Infinity) becomes (0, 1). Without infinities, remains as-is. Des not normalize to length 1.
JmatC.infNormalize = function(value) {
  if (JmatC.isNaN(value)) return jmatc(NaN);

  if (value.re == Infinity) {
    if (value.im == Infinity) return jmatc(1, 1);
    if (value.im == -Infinity) return jmatc(1, -1);
    return jmatc(1, 0);
  }
  if (value.re == -Infinity) {
    if (value.im == Infinity) return jmatc(-1, 1);
    if (value.im == -Infinity) return jmatc(-1, -1);
    return jmatc(-1, 0);
  }
  if (value.im == Infinity) {
    if (value.re == Infinity) return jmatc(1, 1);
    if (value.re == -Infinity) return jmatc(-1, 1);
    return jmatc(0, 1);
  }
  if (value.im == -Infinity) {
    if (value.re == Infinity) return jmatc(1, -1);
    if (value.re == -Infinity) return jmatc(-1, -1);
    return jmatc(0, -1);
  }

  return value.divr(value.absr());
};

// Automatically cache last value. Useful for parameters of statistical distributions that are often the same in repeated calls.
// Cache must be an array (initially []), so that this function can modify it to set the necessary values.
// Function fun is called with z.
// n is cache size
// if n is given, cache contains alternating: index, input0, result0, input1, result1, input2, result2, ... where index is circular pointer to fill in new cache values
// if n is not given, cache contains: input, result
JmatC.calcCache_ = function(z, fun, cache, n) {
  if(n) {
    for(var i = 0; i < n; i++) if(z.eq(cache[i * 2 + 1])) return cache[i * 2 + 2];
    var index = cache[0] || 0;
    index++;
    if(index >= n) index = 0;
    var result = fun(z);
    cache[index * 2 + 1] = z;
    cache[index * 2 + 2] = result;
    cache[0] = index;
    return result;
  } else {
    if(z.eq(cache[0])) return cache[1];
    var result = fun(z);
    cache[0] = z;
    cache[1] = result;
    return result;
  }
};

//Inspired by Wikipedia, Lanczos approximation, precision is around 15 decimal places
JmatC.gamma = function(z) {
  if(z.re == Infinity) return jmatc(Infinity);
  if(JmatC.isNegativeIntOrZero(z)) return jmatc(Infinity, Infinity); // Undirected infinity
  if(z.im == 0) return jmatc(JmatR.gamma(z.re));

  // The internal function that doesn't do internal checks
  var gamma_ = function(z) {
    if(z.re < 0.5) {
      // Use the reflection formula, because, the approximation below is not accurate
      // for values around -6.5+0.1i
      // gamma(1-z)*gamma(z) = pi/sin(pi*z)
      var result = JmatC.PI.div(JmatC.sin(JmatC.PI.mul(z))).div(gamma_(JmatC.ONE.sub(z)));
      if(JmatC.isNaN(result)) result = jmatc(0); // For those values that it can't calculate, it's 0 on the negative side of the complex plane.
      return result;
    }

    var g = 7;
    var p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];

    z = z.subr(1);
    var x = jmatc(p[0]);
    for(var i = 1; i < g + 2; i++) {
      x = x.add(jmatc(p[i]).div(z.addr(i)));
    }
    var t = z.addr(g + 0.5);
    var pisq = Math.sqrt(Math.PI * 2);

    var w = t.pow(z.addr(0.5));
    var e = JmatC.exp(t.neg());
    
    var result = w.mul(e).mul(x).mulr(pisq);
    return result;
  };

  return gamma_(z);
};

JmatC.factorial = function(a) {
  return JmatC.gamma(JmatC.inc(a));
};

//using Stirling series
//logarithm of the gamma function, and more specific branch of the log
JmatC.loggamma = function(z) {
  //the result is way too imprecise if the real part of z is < 0, use the log of the reflection formula
  // loggamma(z) = log(pi/sin(pi*z)) - loggamma(1 - z)
  if(z.re < 0) {
    if(z.im == 0 && z.re == Math.floor(z.re)) return jmatc(NaN); // gamma does not exist here, so it shouldn't return bogus values
    var l = JmatC.log(JmatC.PI.div(JmatC.sin(JmatC.PI.mul(z))));
    // the complex sine goes out of bounds for example for (-4+120i), log(pi / sin(pi*z)) is very roughly approximated by log(2*pi*i) - i*z*pi*sign(z.im) (TODO: the approximation is not fully correct, e.g. for -160-116i)
    if(JmatC.isInfOrNaN(l)) l = JmatC.log(JmatC.newi(2 * Math.PI)).sub(JmatC.newi(-Math.PI).mul(z.im > 0 ? z : z.neg()));
    return l.sub(JmatC.loggamma(JmatC.ONE.sub(z)));
  }

  // The series below has a weird artefact for values near 0 with re > 0. Use actual log(gamma) for that
  if(z.im < 1 && z.im > -1 && z.re < 1 && z.re >= 0) return JmatC.log(JmatC.gamma(z));

  var result = jmatc(0.918938533205); //0.5 * ln(2pi)
  result = result.add(JmatC.subr(z, 0.5).mul(JmatC.log(z)));
  result = result.sub(z);
  result = result.add(z.mulr(12).inv());
  result = result.sub(JmatC.powr(z, 3).mulr(360).inv());
  result = result.add(JmatC.powr(z, 5).mulr(1260).inv());
  result = result.sub(JmatC.powr(z, 7).mulr(1680).inv());
  result = result.add(JmatC.powr(z, 9).mulr(1188).inv());
  return result;
};

//Inverse of the gamma function (not the reciproke, the inverse or "arc" function)
//source: http://mathforum.org/kb/message.jspa?messageID=342551&tstart=0
//lx = log(x + c) / sqrt(2pi); result = lx / lambertw(lx / e) + 0.5
//not very precise
JmatC.gamma_inv = function(value) {
  if(JmatC.isPositive(value) && x > 0.85) { //doesn't work for negative values, nor values smaller than somewhere around 0.85
    // Approximation for positive real x.
    var x = value.re;
    //c = sqrt(2 * pi) / e - gamma(k), where k = the positive zero of the digamma function (1.4616321449683623412626...)
    var c = 0.03653381448490041660;
    var lx = Math.log((x + c) / Math.sqrt(2 * Math.PI));
    return jmatc(lx / JmatR.lambertw(lx / Math.E) + 0.5);
  }

  // TODO: this has problems for |z| > 20. Use more stable root finding
  // TODO: since the current digamma implementation is nothing more than a derivative approximation, I might as well use finvert_secant. Try using digamma anyway if it's ever made more precise.
  var result = JmatC.finvert_newton(value, JmatC.gamma, function(z) {
    // Derivative of gamma function is: gamma function multiplied by digamma function
    return JmatC.gamma(z).mul(JmatC.digamma(z));
  });
  if(!JmatC.near(JmatC.gamma(result), value, 0.01)) return jmatc(NaN);
  return result;
};

// digamma function: psi_0(z)
JmatC.digamma = function(z) {
  // digamma(z) = gamma'(z) / gamma(z)
  // the derivtive gamma'(z) is approximated here as (gamma(z+epsilon) - gamma(z-epsilon)) / (2*epsilon). TODO: use better approximation
  // TODO: for real z, use an approximation such as the Euler Maclaurin formula which gives ln(x) -1/2 x - 1/12 xx + 1/120 xxxx + ....

  // METHOD A: simple approximate derivative of gamma function - works in fact quite well
  // There are two ways: the derivative of loggamma, or, the derivative of gamma divided through gamma. The first is faster [TODO: verify that], but does not work near any z.re that has fractional part 0 or 0.5 due to branch cuts of loggamma.
  var f = JmatR.frac(z.re);
  if(f < 0.001 || f > 0.999 || (f > 0.499 && f < 0.501)) return JmatC.gamma(z.addr(0.0001)).sub(JmatC.gamma(z.subr(0.0001))).divr(0.0002).div(JmatC.gamma(z));
  return JmatC.loggamma(z.addr(0.0001)).sub(JmatC.loggamma(z.subr(0.0001))).divr(0.0002);

  // METHOD B: series - does not work well, requires thousands of iterations and even then the approximate derivative works better
  /*var sum = JmatC.ZERO;
  for(var n = 0; n < 300; n++) {
    var d = z.addr(n).mulr(n + 1);
    sum = sum.add(z.dec().div(d));
  }
  // subtract Euler-Mascheroni constant
  return sum.subr(JmatR.EM);*/
};

JmatC.trigamma = function(z) {
  // The following is implemented in here, but then in such way to avoid redundant identical gamma calls
  //return JmatC.digamma(z.addr(0.0001)).sub(JmatC.digamma(z.subr(0.0001))).divr(0.0002);

  // There are two ways: the derivative of loggamma, or, the derivative of gamma divided through gamma. The first is faster [TODO: verify that], but does not work near any z.re that has fractional part 0 or 0.5 due to branch cuts of loggamma.
  var f = JmatR.frac(z.re);
  if(f < 0.001 || f > 0.999 || (f > 0.499 && f < 0.501)) {
    var a = JmatC.gamma(z.addr(-0.0002));
    var b = JmatC.gamma(z.addr(-0.0001));
    var c = JmatC.gamma(z);
    var d = JmatC.gamma(z.addr(0.0001));
    var e = JmatC.gamma(z.addr(0.0002));

    var d0 = a.sub(c).divr(0.0002).div(b);
    var d1 = c.sub(e).divr(0.0002).div(d);

    return d0.sub(d1).divr(0.0002);
  } else {
    var a = JmatC.loggamma(z.addr(-0.0002));
    var b = JmatC.loggamma(z);
    var c = JmatC.loggamma(z.addr(0.0002));

    var d0 = a.sub(b).divr(0.0002);
    var d1 = b.sub(c).divr(0.0002);

    return d0.sub(d1).divr(0.0002);
  }
};

JmatC.polygamma = function(n, z) {
  if(n.eqr(0)) return JmatC.digamma(z);
  if(n.eqr(1)) return JmatC.trigamma(z);

  // METHOD A: series: Only for positive integer n. Does not work well, requires too many iterations
  /*if(JmatC.isPositiveInt(n))
    var sum = JmatC.ZERO;
    for(var k = 0; k < 300; k++) {
      var d = z.addr(k).pow(n.inc());
      sum = sum.add(d.inv());
    }
    return jmatc(-1).pow(n.inc()).mul(JmatC.factorial(n)).mul(sum);
  }*/

  // METHOD B: with hurwitz zeta
  var h = JmatC.hurwitzzeta(n.inc(), z);
  return jmatc(-1).pow(n.inc()).mul(JmatC.factorial(n)).mul(h);
};

// lower incomplete gamma function
// lowercase gamma(s, x) 
JmatC.incgamma_lower = function(s, z) {
  // METHOD A: in terms of hypergeometric1F1 function
  //return z.pow(s).div(s).mul(JmatC.hypergeometric1F1(s, s.inc(), z.neg()));

  // METHOD B: series expansion - has some problems with division through zero
  // sum_k=0.oo ((-1)^k / k!) * (z^(s+k) / (s + k))
  var result = jmatc(0);
  var kk = jmatc(1);
  var zz = z.pow(s);
  var sign = 1;
  for(var k = 0; k < 30; k++) {
    if(k > 0) {
      kk = kk / k;
      sign = -sign;
      zz = zz.mul(z);
    }
    result = result.add(jmatc(sign * kk).mul(zz).div(s.addr(k)));
  }
  return result;
};

// upper incomplete gamma function
// uppercase GAMMA(s, x) 
JmatC.incgamma_upper = function(s, z) {
  // For negative integer s, gamma, and lower incomplete gamma, are not defined. But the upper is.
  if(JmatC.isNegativeIntOrZero(s)) {
    s = s.addr(1e-7); //twiddle it a bit for negative integers, so that formula in terms of lower gamma sort of works... TODO: use better approximation
  }

  return JmatC.gamma(s).sub(JmatC.incgamma_lower(s, z));
};

JmatC.gamma_p_cache_ = []; // cache used because s is often constant between gamma_p calls

// regularized lower incomplete gamma function (lower gamma_regularized, regularized_gamma)
// P(s, x) = gamma(s, z) / GAMMA(s)
// Note: the derivative of this function is: (e^(-z) * z^(s-1))/GAMMA(s)
JmatC.gamma_p = function(s, z) {
  if(JmatC.isNegativeIntOrZero(s)) return jmatc(1);
  var g = JmatC.calcCache_(s, JmatC.gamma, JmatC.gamma_p_cache_); // gamma(s)
  return JmatC.incgamma_lower(s, z).div(g);
};

// regularized upper incomplete gamma function (upper gamma_regularized, regularized_gamma)
// Q(s, x) = 1 - P(s, x) = GAMMA(s, z) / GAMMA(s)
// Note: the derivative of this function is: -(e^(-z) * z^(s-1))/GAMMA(s)
JmatC.gamma_q = function(s, z) {
  if(JmatC.isNegativeIntOrZero(s)) return jmatc(0);
  return JmatC.ONE.sub(JmatC.gamma_p(s, z));
};

// One possible approximation series for inverse gamma P. Valid for real values with p in range 0-1 and positive s. Possibly more.
JmatC.gamma_p_inv_series_1_ = function(s, p) {
  var s1 = s.inc();
  var s1s = s1.mul(s1);
  var s2 = s.addr(2);
  var s2s = s2.mul(s2);
  var ss = s.mul(s);
  var sss = ss.mul(s);
  var ssss = sss.mul(s);
  var s3 = s.addr(3);
  var s4 = s.addr(4);

  var c = [0, 1];
  c[2] = s.inc().inv();
  c[3] = s.mulr(3).addr(5).div(s1s.mul(s2).mulr(2));
  c[4] = ss.mulr(8).add(s.mulr(33)).addr(31).div(s1s.mul(s1).mul(s2).mul(s3).mulr(3));
  c[5] = ssss.mulr(125).add(sss.mulr(1179)).add(ss.mulr(3971)).add(s.mulr(5661)).addr(2888).div(s1s.mul(s1s).mul(s2s).mul(s3).mul(s4).mulr(24));

  var r = p.mul(JmatC.gamma(s1)).pow(s.inv());

  return JmatC.powerSeries(c, c.length, JmatC.ZERO, r);
};

// Inverse regularized gamma P
// Finds z for p == gamma_p(s, z)
// This is an *approximation* of inverse of incomplete regularized gamma (P) - it works for real values with p in range 0-1 and positive s. Good enough for qf of chi square distribution
JmatC.gamma_p_inv = function(s, p) {
  // Return NaN for unsupported values to prevent bogus results
  if(!JmatC.isReal(p) || !JmatC.isReal(s) || p.re < 0 || p.re > 1 || s.re < 0) return jmatc(NaN);

  // TODO: more complete support of the entire complex domain
  return JmatC.gamma_p_inv_series_1_(s, p);
};

// Calculates gamma(x) / gamma(y), and can cancel out negative integer arguments if both are negative integers in some cases. That is useful for functions like beta, binomial, ...
// Uses the formula gamma(a)/gamma(b) = (-1)^(a-b) * gamma(-b+1) / gamma(-a+1) if necessary
// It is also optimized to use for loops if x and y are nearby integers, ...
JmatC.gammaDiv_ = function(x, y) {
  if(x.eq(y)) return JmatC.ONE; // For the "combined" function, this is considered correct even for negative integers...

  if(JmatC.isInfOrNaN(x) || JmatC.isInfOrNaN(y)) {
    return jmatc(NaN);
  }

  if(JmatC.isNegativeIntOrZero(y) && (x.re > 0 || !JmatC.isInt(x))) return JmatC.ZERO; // division of non-infinity through infinity

  if(JmatC.isInt(x) && JmatC.isInt(y)) {
    if(x.re <= 0 && y.re <= 0) {
      var sign = JmatR.isOdd(x.re - y.re) ? -1 : 1;
      return JmatC.gammaDiv_(y.neg().addr(1), x.neg().addr(1)).mulr(sign);
    }

    if(x.re > 0 && y.re > 0 && JmatR.dist(x.re, y.re) < 16) {
      if(x.re > y.re) {
        var result = y.re
        for(var z = y.re + 1; z < x.re; z++) result *= z;
        return jmatc(result);
      } else {
        var result = 1 / x.re
        for(var z = x.re + 1; z < y.re; z++) result /= z;
        return jmatc(result);
      }
    }
  }

  return JmatC.gamma(x).div(JmatC.gamma(y));
};

// Similar to JmatC.gammaDiv_, but takes 3 arguments and cancels out more if possible. Returns gamma(a) / (gamma(b) * gamma(c))
JmatC.gammaDiv12_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(JmatC.isNegativeIntOrZero(b)) {
    if(JmatC.isNegativeIntOrZero(c) && (a.re > 0 || !JmatC.isInt(a))) return JmatC.ZERO; // division of non-infinity through infinity
    return JmatC.gammaDiv_(a, b).div(JmatC.gamma(c));
  } else {
    return JmatC.gammaDiv_(a, c).div(JmatC.gamma(b));
  }
};

// Similar to JmatC.gammaDiv_, but takes 3 arguments and cancels out more if possible. Returns (gamma(a) * gamma(b)) / gamma(c)
JmatC.gammaDiv21_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(JmatC.isNegativeIntOrZero(a)) {
    return JmatC.gammaDiv_(a, c).mul(JmatC.gamma(b));
  } else {
    return JmatC.gammaDiv_(b, c).mul(JmatC.gamma(a));
  }
};

// Similar to JmatC.gammaDiv_, but takes 4 arguments and cancels out more if possible. Returns (gamma(a) * gamma(b)) / (gamma(c) * gamma(d))
JmatC.gammaDiv22_ = function(a, b, c, d) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(JmatC.isNegativeIntOrZero(a) == JmatC.isNegativeIntOrZero(c)) {
    return JmatC.gammaDiv_(a, c).mul(JmatC.gammaDiv_(b, d));
  } else {
    return JmatC.gammaDiv_(a, d).mul(JmatC.gammaDiv_(b, c));
  }
};

// Calculates log(gamma(x) / gamma(y)) = loggamma(x) - loggamma(y), and can cancel out negative integer arguments if both are negative integers in some cases. That is useful for functions like beta, binomial, ...
// Uses the formula gamma(a)/gamma(b) = (-1)^(a-b) * gamma(-b+1) / gamma(-a+1) if necessary
JmatC.loggammaDiv_ = function(x, y) {
  if(x.eq(y)) return JmatC.ZERO; // For the "combined" function, this is correct even for negative integers...

  if(JmatC.isInfOrNaN(x) || JmatC.isInfOrNaN(y)) {
    return jmatc(NaN);
  }

  if(JmatC.isNegativeIntOrZero(y) && (x.re > 0 || !JmatC.isInt(x))) return jmatc(-Infinity); // division of non-infinity through infinity ==> log(0)

  if(JmatC.isInt(x) && JmatC.isInt(y)) {
    if(x.re <= 0 && y.re <= 0) {
      var sign = JmatR.isOdd(x.re - y.re) ? -1 : 1;
      var result = JmatC.loggammaDiv_(y.neg().addr(1), x.neg().addr(1));
      if(sign == -1) result = result.add(JmatC.newi(Math.PI)); // log(-x) = log(x) + i*pi
      return result;
    }
  }

  return JmatC.loggamma(x).sub(JmatC.loggamma(y));
};

// Similar to JmatC.loggammaDiv_, but takes 3 arguments and cancels out more if possible. Returns loggamma(a) - (loggamma(b) + loggamma(c))
JmatC.loggammaDiv12_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(JmatC.isNegativeIntOrZero(b)) {
    if(JmatC.isNegativeIntOrZero(c) && (a.re > 0 || !JmatC.isInt(a))) return jmatc(-Infinity); // division of non-infinity through infinity ==> log(0)
    return JmatC.loggammaDiv_(a, b).sub(JmatC.loggamma(c));
  } else {
    return JmatC.loggammaDiv_(a, c).sub(JmatC.loggamma(b));
  }
};

// Similar to JmatC.loggammaDiv_, but takes 3 arguments and cancels out more if possible. Returns (loggamma(a) + loggamma(b)) - loggamma(c)
JmatC.loggammaDiv21_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(JmatC.isNegativeIntOrZero(a)) {
    return JmatC.loggammaDiv_(a, c).add(JmatC.loggamma(b));
  } else {
    return JmatC.loggammaDiv_(b, c).add(JmatC.loggamma(a));
  }
};

// Similar to JmatC.loggammaDiv_, but takes 4 arguments and cancels out more if possible. Returns (loggamma(a) + loggamma(b)) - (loggamma(c) + loggamma(d))
// To have only 3 values, set a, b, c or d to 1 (it will be fast for that)
JmatC.loggammaDiv2_ = function(a, b, c, d) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(JmatC.isNegativeIntOrZero(a) == JmatC.isNegativeIntOrZero(c)) {
    return JmatC.loggammaDiv_(a, c).add(JmatC.loggammaDiv_(b, d));
  } else {
    return JmatC.loggammaDiv_(a, d).add(JmatC.loggammaDiv_(b, c));
  }
};

// Returns a particular non-principal complex branch of sqrt(a * b)
// For AGM and GHM, a particular branch of the complex sqrt must be returned, otherwise it hops to other branches and gives wrong result
JmatC.agmMulSqrt_ = function(a, b) {
  // Note: if a and b are real but negative, with a*b being positive, still a different branch is used.
  if(JmatC.isPositive(a) && JmatC.isPositive(b)) {
    return JmatC.sqrt(a.mul(b));
  } else {
    return JmatC.sqrt(JmatC.abs(a.mul(b))).mul(JmatC.exp(JmatC.I.mul(JmatC.arg(a).add(JmatC.arg(b)).divr(2))));
  }
};

//Arithmetic-Geometric mean (iteratively calculated)
JmatC.agm = function(a, b) {
  // Wikipedia: We have the following inequality involving the Pythagorean means {H, G, A} and iterated Pythagorean means {HG, HA, GA}:
  // min <= H <= HG <= (G == HA) <= GA <= A <= max
  // For completeness, also (with L = logarithmic mean, R = RMS, C = contraharmonic mean):
  // min <= H <= G <= L <= A <= R <= C <= max

  // The AGM is normally only defined for positive real numbers, but it can be extended to whole complex plane, as is done here.
  
  //avoid imprecisions for special cases
  if(a.eq(b.neg()) || a.eq(JmatC.ZERO) || b.eq(JmatC.ZERO)) return jmatc(0);
  var real = JmatC.isReal(a) && JmatC.isReal(b) && (a.re < 0) == (b.re < 0);

  var a2, b2;
  for(var i = 0; i < 60; i++) {
    if(a.eq(b)) break;
    a2 = a.add(b).divr(2);
    b2 = JmatC.agmMulSqrt_(a, b);
    a = a2;
    b = b2;
  }

  if(real) a.im = 0; // it may be 1e-16 or so due to numerical error (only if both inputs are real and have same sign)
  
  return a;
};

//Geometric-Harmonic mean (iteratively calculated)
JmatC.ghm = function(a, b) {
  // Not really defined for negative and complex numbers, but when using JmatC.agmMulSqrt_ it looks smooth in the 2D plot
  // NOTE: An alternative, that returns different values for neg/complex (but same for positive reals) is: return JmatC.agm(a.inv(), b.inv()).inv();
  var a2, b2;
  for(var i = 0; i < 60; i++) {
    if(a.eq(b)) break;
    a2 = JmatC.agmMulSqrt_(a, b);
    b2 = JmatC.TWO.div(a.inv().add(b.inv()));
    a = a2;
    b = b2;
  }
  return a;
};

// Bessel function of the first kind
// Mostly an approximation, there are problems with high z
// TODO: make faster and fix problems
JmatC.besselj = function(n, z) {
  var pi = Math.PI;

  if(z.im == 0 && n.im == 0 && z.re < 0 && z.re * z.re < (n.re + 1) / 10) {
    return JmatC.inv(JmatC.gamma(n.inc())).mul(z.divr(2).pow(n));
  } else if(JmatC.abs(z).re < 20) {
    // The gamma functions give NaN if n is negative integer < -1. TODO: also for n = 0 and neg won't help there, fix that
    var negintn = JmatC.isNegativeInt(n);
    if(negintn) n = n.neg();
    var result = jmatc(0);
    var m = 1;
    for(var i = 0; i < 50; i++) {
      var i1 = jmatc(i + 1);
      var d = JmatC.gamma(i1).mul(JmatC.gamma(n.add(i1)));
      var term = jmatc(m).div(d).mul(z.divr(2).pow(jmatc(2 * i).add(n)));
      m = -m;
      result = result.add(term);
    }
    if(negintn && JmatR.isOdd(n.re)) result = result.neg(); //besselj(-n, x) = (-1)^n * besselj(n, x)
    return result;
  } else if(z.im == 0) {
    // Something is wrong with this formula, see the 2D plot
    var a = JmatC.sqrt(jmatc(2/pi).div(z));
    var b = z.sub(n.mulr(pi/2)).subr(pi/4);
    return a.mul(JmatC.cos(b));
  } else {
    // Something is wrong with this formula, see the 2D plot
    var s = z.im > 0 ? -1 : 1
    var a = z.sub(n.mulr(pi/2)).subr(pi/4);
    var b = JmatC.sqrt(z.mulr(pi*2));
    return JmatC.exp(JmatC.newi(s).mul(a)).div(b);
  }
};

// Bessel function of the second kind
// Mostly an approximation, there are problems with high z
// TODO: make faster and fix problems
JmatC.bessely = function(n, z) {
  var pi = Math.PI;

  if(JmatC.abs(z).re < 15) {
    // Y_a(x) = (J_a(x)*cos(a*pi) - J_-a(x)) / sin(a*pi)
    // For integer n, Y_n(x) = lim_a_to_n(Y_a(x)
    if(n.re == Math.floor(n.re)) n = n.addr(0.000000001);

    var a = JmatC.besselj(n, z);
    var b = JmatC.cos(n.mulr(Math.PI));
    var c = JmatC.besselj(n.neg(), z);
    var d = JmatC.sin(n.mulr(Math.PI));
    return a.mul(b).sub(c).div(d);
  } else if(z.im == 0) {
    // Something is wrong with this formula, see the 2D plot
    var a = JmatC.sqrt(jmatc(2/pi).div(z));
    var b = z.sub(n.mulr(pi/2)).subr(pi/4);
    return a.mul(JmatC.sin(b));
  } else {
    // Something is wrong with this formula, see the 2D plot
    var s = z.im > 0 ? -1 : 1
    var a = z.sub(n.mulr(pi/2)).subr(pi/4);
    var b = JmatC.sqrt(z.mulr(pi*2));
    return JmatC.newi(-s).mul(JmatC.exp(JmatC.newi(s).mul(a)).div(b));
  }
};

// Hankel function of the first kind (approximation)
// TODO: this one is very slow. Make faster and more precise
JmatC.hankel1 = function(n, z) {
  // There are numerical imprecisions for e.g. hankel2(-8-8i, i).
  // So when needed, apply the formula: H1_-a(x) = exp(a*pi*i)*H1_a(x)
  if(n.im < 0) return JmatC.exp(n.mul(JmatC.newi(-Math.PI))).mul(JmatC.hankel1(n.neg(), z));

  return JmatC.besselj(n, z).add(JmatC.bessely(n, z).mul(JmatC.I));
};

// Hankel function of the second kind (approximation)
// TODO: this one is very slow. Make faster and more precise
JmatC.hankel2 = function(n, z) {
  // There are numerical imprecisions for e.g. hankel2(-8-8i, i).
  // So when needed, apply the formula: H2_-a(x) = exp(-a*pi*i)*H2_a(x)
  if(n.im > 0) return JmatC.exp(n.mul(JmatC.newi(Math.PI))).mul(JmatC.hankel2(n.neg(), z));

  return JmatC.besselj(n, z).sub(JmatC.bessely(n, z).mul(JmatC.I));
};

// Modified bessel function of the first kind (approximation)
// TODO: this one is very slow. Make faster and more precise
JmatC.besseli = function(n, z) {
  var result = JmatC.I.pow(n.neg()).mul(JmatC.besselj(n, z.mul(JmatC.I)));
  if(z.im == 0 && n.im == 0 && JmatR.near(result.im, 0, 1e-10)) result.im = 0;
  return result;
};

// Modified bessel function of the second kind (approximation)
// TODO: this one is very slow. Make faster and more precise
JmatC.besselk = function(n, z) {
  var result;
  if(z.im >= 0) {
    result = JmatC.I.pow(n.inc()).mulr(Math.PI / 2).mul(JmatC.hankel1(n, JmatC.I.mul(z)));
  } else {
    result = JmatC.newi(-1).pow(n.inc()).mulr(Math.PI / 2).mul(JmatC.hankel2(n, JmatC.newi(-1).mul(z)));
  }
  if(z.im == 0 && n.im == 0 && JmatR.near(result.im, 0, 1e-3 /*this one is super imprecise...*/)) result.im = 0;
  return result;
};

//pl = 3^(-2/3) for airy, 3^(-4/3) for bairy
//pr = 3^(-4/3) for airy, 3^(-5/6) for bairy
//s = -1 for airy, +1 for bairy
JmatC.airyloop_ = function(z, pl, pr, s) {
  var gl = 1.3541179394264004169; // GAMMA(2/3)
  var gr = 0.8929795115692492112; // GAMMA(4/3)
  var zzz = z.mul(z).mul(z);
  var zr = JmatC.ONE;
  var zl = z;
  var r = JmatC.ONE;
  var kk = 1;
  var result = JmatC.ZERO;
  for(var k = 0; k < 30; k++) {
    if (k > 0) {
      kk = kk * k;
      gl = gl * (k + 2/3 - 1); 
      gr = gr * (k + 4/3 - 1);
      pl /= 9; // add another -2k to the power
      pr /= 9; // add another -2k to the power
      zr = zr.mul(zzz);
      zl = zl.mul(zzz);
    }

    var rl = zr.mulr(pl / kk / gl);
    var rr = zl.mulr(pr / kk / gr);
    if(JmatC.isNaN(rl) || JmatC.isNaN(rr)) break;
    if(rl.eqr(0) && rr.eqr(0)) break;
    result = result.add(rl).add(rr.mulr(s));
  }

  return result;
};

// Airy function Ai(x)
JmatC.airy = function(z) {
  if(z.absr() > 8) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < 2 * Math.PI / 3) {
      var d = z.powr(1/4).mul(JmatC.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return JmatC.exp(zeta.neg()).div(d.mulr(2));
    } else {
      var zm = z.neg();
      var d = zm.powr(1/4).mul(JmatC.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return JmatC.sin(zeta.addr(Math.PI / 4)).div(d);
    }
  }

  // METHOD A: definition with hypergeometric0F1
  /*var zzz = z.mul(z).mul(z);
  var a = JmatC.hypergeometric0F1(jmatc(2/3), zzz.mulr(1/9)).mulr(Math.pow(3,-2/3) / JmatR.gamma(2/3));
  var b = JmatC.hypergeometric0F1(jmatc(4/3), zzz.mulr(1/9)).mul(z).mulr(Math.pow(3,-1/3) / JmatR.gamma(1/3));
  return a.sub(b);*/

  // METHOD B: summation:
  // SUM_k=0..oo 3^(-2k-2/3)/(k!*GAMMA(k+2/3))*x^(3k) - SUM_k=0..oo 3^(-2k-4/3)/(k!*GAMMA(k+4/3))*x^(3k+1)
  // This is basically the same as the hypergeometric definitions, but faster
  var pl = Math.pow(3, -2/3);
  var pr = Math.pow(3, -4/3);
  return JmatC.airyloop_(z, pl, pr, -1);
};

// Airy Bi function Bi(x)
JmatC.bairy = function(z) {
  if(z.absr() > 10) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < Math.PI / 3) {
      var d = z.powr(1/4).mul(JmatC.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return JmatC.exp(zeta).div(d);
    } else {
      var zm = z.neg();
      var d = zm.powr(1/4).mul(JmatC.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return JmatC.cos(zeta.addr(Math.PI / 4)).div(d);
    }
  }

  // SUM_k=0..oo 3^(-2k-1/6)/(k!*GAMMA(k+2/3))*x^(3k) + SUM_k=0..oo 3^(-2k-5/6)/(k!*GAMMA(k+4/3))*x^(3k+1)
  var pl = Math.pow(3, -1/6);
  var pr = Math.pow(3, -5/6);
  return JmatC.airyloop_(z, pl, pr, +1);
};

//pl = 3^(-1/3) for airy, 3^(+1/6) for bairy
//pr = 3^(-5/3) for airy, 3^(-7/6) for bairy
//s = -1 for airy, +1 for bairy
JmatC.airy_deriv_loop_ = function(z, pl, pr, s) {
  var gl = 2.6789385347077476336; // GAMMA(1/3)
  var gr = 0.9027452929509336112; // GAMMA(5/3)
  var zzz = z.mul(z).mul(z);
  var zr = JmatC.ONE;
  var zl = z.mul(z);
  var r = JmatC.ONE;
  var kk = 1;
  var result = JmatC.ZERO;
  for(var k = 0; k < 30; k++) {
    if (k > 0) {
      kk = kk * k;
      gl = gl * (k + 1/3 - 1); 
      gr = gr * (k + 5/3 - 1);
      pl /= 9; // add another -2k to the power
      pr /= 9; // add another -2k to the power
      zr = zr.mul(zzz);
      zl = zl.mul(zzz);
    }

    var rl = zr.mulr(pl / kk / gl);
    var rr = zl.mulr(pr / kk / gr);
    if(JmatC.isNaN(rl) || JmatC.isNaN(rr)) break;
    if(rl.eqr(0) && rr.eqr(0)) break;
    result = result.add(rl.mulr(s)).add(rr);
  }

  return result;
};

// Derivative Airy function Ai'(x)
JmatC.airy_deriv = function(z) {
  if(z.absr() > 8) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < 2 * Math.PI / 3) {
      var d = z.powr(-1/4).mul(JmatC.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return JmatC.exp(zeta.neg()).div(d.mulr(2)).neg();
    } else {
      var zm = z.neg();
      var d = zm.powr(-1/4).mul(JmatC.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return JmatC.cos(zeta.addr(Math.PI / 4)).div(d).neg();
    }
  }

  // -SUM_k=0..oo 3^(-2k-1/3)/(k!*GAMMA(k+1/3))*x^(3k) + SUM_k=0..oo 3^(-2k-5/3)/(k!*GAMMA(k+5/3))*x^(3k+2)
  var pl = Math.pow(3, -1/3);
  var pr = Math.pow(3, -5/3);
  return JmatC.airy_deriv_loop_(z, pl, pr, -1);
};

// Derivative Airy Bi function Bi'(x)
JmatC.bairy_deriv = function(z) {
  if(z.absr() > 10) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < Math.PI / 3) {
      var d = z.powr(-1/4).mul(JmatC.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return JmatC.exp(zeta).div(d);
    } else {
      var zm = z.neg();
      var d = zm.powr(-1/4).mul(JmatC.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return JmatC.sin(zeta.addr(Math.PI / 4)).div(d);
    }
  }

  // SUM_k=0..oo 3^(-2k+1/6)/(k!*GAMMA(k+1/3))*x^(3k) + SUM_k=0..oo 3^(-2k-7/6)/(k!*GAMMA(k+5/3))*x^(3k+2)
  var pl = Math.pow(3, +1/6);
  var pr = Math.pow(3, -7/6);
  return JmatC.airy_deriv_loop_(z, pl, pr, +1);
};


// This integral representation of riemann zeta works quite well and is valid for all s. This is quite slow though. And in practice it does not work for z.im < -42 or z.im > 42.
// That because the integrand fluctuates heavily between negative and positive values I think.
// Anyway, it is better than the eta series formula for the very particular values of x.re = 1 and x.im = a multiple of 9.0625, so used for that...
// The formula: zeta(s) = 2 * integrate_0..oo sin(s.atan(t)) / ((1+t^2)^(s/2) * (e^(2*pi*t) - 1) dt + 1/2 + 1/(s-1)
JmatC.zetaint_ = function(s) {
  var it = function(s, t) {
    var n = JmatC.sin(s.mul(JmatC.atan(t)));
    var ta = (t.mul(t).inc()).pow(s.divr(2));
    var tb = JmatC.exp(t.mulr(2*Math.PI)).dec();
    var t = ta.mul(tb);
    if(n.eq(t)) return JmatC.ONE; //avoid 0/0
    return n.div(t);
  };

  // Binary search to the place where the function starts becoming zero forever
  var q = JmatC.ONE;
  for(;;) {
    var q1 = it(s, q);
    var q2 = it(s, q.inc());
    if(q1.abssqr() < 1e-4 && q1.abssqr() < 1e-4) break;
    if(q.re > 1000) break; //don't do too much doublings or we're doing as much calculations in here as in the integral...
  
    q = q.mulr(2);
  }

  var maxt = q.re;
  var step = 0.01;

  var r = JmatC.ZERO;
  for(var t = 0; t < maxt; t+=step) {
    var a = it(s, jmatc(t)).mulr(step);
    r = r.add(a);
    step *= 1.01;
  }

  return r.mulr(2).addr(0.5).add(s.dec().inv());
};

// Riemann zeta function (infinite sum of reciprocal powers)
JmatC.zeta = function(s) {
  if(s.eq(JmatC.ZERO)) return jmatc(-0.5);
  if(s.eqr(-1)) return jmatc(-1/12);
  if(s.eqr(1)) return jmatc(+Infinity);
  if(s.eqr(3)) return jmatc(1.202056903159594); // Apery's constant
  if(s.re == +Infinity) return jmatc(1);
  if(JmatC.isEven(s) && s.re < 0) return jmatc(0); //it's 0 at all negative even integers

  // TODO: Use better algorithms:
  // Borwein algorithm for s close to real line
  // Riemann-Siegel formula for s with large imaginary party
  // Euler-Maclaurin summation in all other cases
  // The reflection formula where appliccable

  var a1 = s.dec().absr();
  if(a1 < 10) {
    // Laurent series around 1: 1/(s-1) + SUM_n=0..oo (-1)^n / n! * stieltjes_n * (s-1)^n
    // For 30 terms, seems correct in a radius of around 10-15 around the point 1
    var s1 = s.dec();
    var result = s1.inv();
    var ss = JmatC.ONE;
    var num = a1 < 5 ? 10 : 30;
    for(var i = 0; i < num; i++) {
      result = result.add(jmatc(JmatC.stieltjes_zeta[i]).mul(ss));
      ss = ss.mul(s1);
    }
    return result;
  } else if(s.re >= 5) {
    // The riemann zeta infinite series: SUM_n=1..oo n^(-s) = 1/1^s + 1/2^s + 1/3^s + ...
    // Converges for Re(s) > 1, but in practice only converges properly for bigger Re(s)
    var result = jmatc(0);
    for(var i = 1; i < 32; i++) {
      result = result.add(jmatc(i).pow(s).inv());
    }
    return result;
  } else if(s.re >= 0.5 /*the articles make it seem as if this loop is for s.re > 0, but it only converges well for s.re 0.5*/) {
    // This is a series that is convergent (but not absolutely) for s.re > 0 instead of s.re > 1. In practice, it does not converge well, and only somewhat for s.re > 0.5, and even there not good...
    // eta(s) = SUM_n=1..oo (-1)^(n+1) / n^s = (1-2^(1-s)) * zeta(s) ==> zeta(s) = 1/(1-2^(1-s)) * (SUM_n=1..oo (-1)^(n+1) / n^s)
    if(JmatC.near(s, jmatc(1, 18.125), 0.1)) return JmatC.zetaint_(s); // The code below has a problem for s of the form 1 + n*9.0625i, with n any negative or positive integer
    var result = jmatc(0);
    var p = JmatC.TWO.pow(JmatC.ONE.sub(s));
    var a = JmatC.ONE.div(JmatC.ONE.sub(p));
    var n = 1;
    for(var i = 1; i < 128; i++) {
      var r = jmatc(n).div(jmatc(i).pow(s));
      result = result.add(r);
      n = -n;
    }
    return result.mul(a);
  } else {
    // Reflection: The functional equation: zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
    var s2 = JmatC.ONE.sub(s);
    var a = JmatC.TWO.pow(s);
    var b = JmatC.PI.pow(s.subr(1));
    var c = JmatC.sin(s.mulr(Math.PI / 2));
    var d = JmatC.gamma(s2);
    var e = JmatC.zeta(s2);
    return a.mul(b).mul(c).mul(d).mul(e);
  }
};

// Dirichlet eta function
JmatC.eta = function(s) {
  //The calculation only works for s.re > 0.5, so use reflection formula if needed
  if(s.re < 0.5) {
    s = s.neg();

    var a = (JmatC.ONE.sub(JmatC.TWO.pow(s.neg().subr(1))));
    var b = (JmatC.ONE.sub(JmatC.TWO.pow(s.neg())));
    var c = a.div(b).mulr(2).mul(JmatC.PI.pow(s.neg().subr(1)));
    var d = JmatC.sin(JmatC.PI.mul(s).divr(2));
    var e = JmatC.gamma(s);
    var f = JmatC.eta(s.inc());
    return c.mul(s).mul(d).mul(e).mul(f);
  }

  // Borwein's method
  var n = 50;

  // function works on real values
  var d = function(k, n) {
    var result = 0;
    for(var i = 0; i <= k; i++) {
      var a = JmatR.factorial(n + i - 1) / JmatR.factorial(n - i);
      var b = Math.pow(4, i) / JmatR.factorial(2 * i);
      result += a * b;
    }
    return result * n;
  };

  var dn = d(n, n);

  var result = jmatc(0);
  var sign = jmatc(1); //alternates
  for(var k = 0; k < n; k++) {
    var a = sign.mulr(d(k, n) - dn);
    var b = jmatc(k + 1).pow(s);
    result = result.add(a.div(b));
    sign = sign.neg();
  }

  return result.mulr(-1 / dn);

  /*
  NOTE:
  Riemman zeta, Direchlet lambda and Direchlet eta are related as follows:
  zeta(s)/(2^s) = lambda(s)/(2^s - 1) = eta(s)/(2^s - 2)
  also:
  zeta(s) + eta(s) = 2*lambda(s)
  So one of these could be expressed as the other, so that a more precise one can be used.
  E.g. the Borwein's method here looks more precise than the current formulas used in JmatC.zeta...
  */
};

// Dirichlet lambda function
// TODO: use numerically more precise formula
JmatC.lambda = function(s) {
  // definition: lambda(s) = (1 - s^(-s))*zeta(s)
  return JmatC.ONE.sub(JmatC.TWO.pow(s.neg())).mul(JmatC.zeta(s));
};

/*
zeta(1 / (s - 1)) * SUM_n=0..oo 1/(n+1) SUM_k=0..n ((-1)^k binomial(n, k) (q+k)^(1-s))
            k=0  k=1  k=2  k=3
             +    -    +    -
n=0 1    *   1
n=1 1/2  *   1    1
n=2 1/3  *   1    2    1
n=3 1/4  *   1    3    3    1
           +----------------
     [2.0833333, -1.916666, 1.083333, -0.25] ---> the result for steps = 4
print out with JSON.stringify(JmatC.hurwitzzeta_generate_hasse_table_(30)), and use in JmatC.hurwitzzeta_hasse_series_.
*/
JmatC.hurwitzzeta_generate_hasse_table_ = function(steps) {
  var result = [];
  for(var n = 0; n < steps; n++) {
    result[n] = 0;
    var sign = 1;
    for(var k = 0; k <= n; k++) {
      result[k] += sign * JmatR.pascal_triangle(n, k) / (n + 1);
      sign = -sign;
    }
  }
  return result;
};

JmatC.hurwitzzeta_hasse_tables_ = [];

// Series by Helmut Hasse: zeta(1 / (s - 1)) * SUM_n=0..oo 1/(n+1) SUM_k=0..n ((-1)^k binomial(n, k) (q+k)^(1-s))
// Requires: in theory: complex s not equal to 1, real q > 0 (or > -1 according to other sources)
// In practice: any complex s with s.im > 0, complex q with q.re > 2 (for smaller q.re, it becomes very imprecise)
JmatC.hurwitzzeta_hasse_series_ = function(s, q) {
  // Precomputed table implementation, exact same result as the slower version but faster.
  // Generated with JSON.stringify(JmatC.hurwitzzeta_generate_binomial_table_(30))
  // Table size must exactly match the N in "k < N" in the loop below.
  // TODO: there is a weird noisy area in the center of the following plot. fix it? Jmat.plotComplex(function(z){return Jmat.hurwitzzeta(-5, z); })
  var N = 30;
  if(!JmatC.hurwitzzeta_hasse_tables_[N]) {
    JmatC.hurwitzzeta_hasse_tables_[N] = JmatC.hurwitzzeta_generate_hasse_table_(N);
  }
  var table = JmatC.hurwitzzeta_hasse_tables_[N];
  var result = JmatC.ZERO;
  for(var k = 0; k < N; k++) {
    result = result.add(q.addr(k).pow(JmatC.ONE.sub(s)).mulr(table[k]));
  }
  result = result.mul(s.dec().inv());
  return result;
};

// A series that works for negative s, but q only in range 0-1
// requirements: s.re < 0 && q.re > 0 && q.re <= 1 && q.im == 0
// TODO: this one actually has more potential than it seemed, it works best of all formulas for JmatC.hurwitzzeta_cos_series_(jmatc(-10,-10), jmatc(-2.1)). Investigate area where it works and make more use of it.
JmatC.hurwitzzeta_cos_series_ = function(s, q) {
  // Series representation 25.11.9 from NIST handbook
  var s1 = s.rsub(1);
  var g = JmatC.gamma(s1).mulr(2).div(jmatc(2 * Math.PI).pow(s1));
  var sum = JmatC.ZERO;
  for(var n = 1; n < 30; n++) {
    var c = JmatC.cos(s1.mulr(Math.PI / 2).sub(q.mulr(2 * Math.PI * n)));
    sum = sum.add(c.div(jmatc(n).pow(s1)));
  }
  return sum.mul(g);
};

// The series of the definition of hurwitz zeta.
// This algorithm works reliable for, according to my practical testing: s with s.re > 2 and complex q with q.re > 0
// Theoretically it works for s.re > 1 but it converges too slow there, and would work for any complex s with q.re > 0, but it seems to work for some (but not reliably) negative s as well, and work NOT well for |q| too big.
JmatC.hurwitzzeta_simple_series_ = function(s, q) {
  if(JmatC.isNegativeIntOrZero(q)) return jmatc(Infinity);
  // would work for s.re > 1, but too slow convergence

  var result = JmatC.ZERO;
  for(var n = 0; n < 30; n++) {
    var r = q.addr(n).pow(s).inv();
    result = result.add(r);
  }

  return result;
};

// Euler-Maclaurin algorithm for hurwitz zeta, basically a sped up version of the simple series
// See (7.1) in paper "efficient algorithm for accelerating the convergence of ..." by linas vepstas.
// This algorithm works for: complex s with s.re > 0, and complex q with q.re > -20 (or even > -1). For q with very large negative real part, it goes numerically very wrong.
JmatC.hurwitzzeta_euler_ = function(s, q) {
  var result = JmatC.ZERO;
  var N = 25; // D/2 + 10 where D is desired decimal digits of precision
  var p = 15; // TODO: find better value for it? It should not go above 15 though, the bernoulli number table ends at 30 (and k*2 is needed)
  var fn = q.addr(N).pow(s).inv(); // f(N)
  var ig = s.dec().inv().mul(q.addr(N).pow(s.dec()).inv());
  var sum1 = JmatC.ZERO;
  for(var k = 0; k < N; k++) {
    var r = q.addr(k).pow(s).inv(); // f(k)
    sum1 = sum1.add(r);
  }
  var ss = s; //for the derivative
  var kk = 1; //for the (2k)!
  var sum2 = JmatC.ZERO;
  for(var k = 1; k <= p; k++) {
    kk *= (2 * k - 1) * (2 * k);
    ss = ss.mul(s.addr(k * 2 - 1)).mul(s.addr(k * 2));
    var d = ss.div(q.addr(N).pow(s.addr(2 * k + 1))).neg();
    sum2 = sum2.add(d.mulr(JmatC.bernoulli[2 * k] / kk));
  }
  return sum1.sub(sum2).add(fn.divr(2)).add(ig);
};

// Hurwitz zeta function
// Currently works ok for complex s and q with real parts > 0. For negative real parts in s or q, there are some zones of very high imprecision, but works quite well for reasonable values.
// TODO: fix several problems, e.g. see the plot Jmat.plotComplex(function(z){return Jmat.hurwitzzeta(-10, z); })
// TODO: support the cases outside of the currently supported domain, e.g. large negative values, ...
JmatC.hurwitzzeta = function(s, q) {
  // Bernoulli polynomials for a few small negative integer s values
  // It's not that necessary to have this and the riemann zeta conversion below, but it helps as a check that the algorithms below are correct, if not these values would show up as different lines in 2D plots
  if(s.eqr(0)) return jmatc(0.5).sub(q);
  if(s.eqr(-1)) return jmatc(-1/12.0).add(q.mulr(0.5)).add(q.mul(q).mulr(-0.5));

  if(q.eqr(1)) return JmatC.zeta(s); // riemann zeta
  
  if(s.re > 2 && q.re > 0) {
    return JmatC.hurwitzzeta_simple_series_(s, q);
  }

  if(s.re >= 0 && q.re > -10) {
    return JmatC.hurwitzzeta_euler_(s, q);
  }

  // Only a very limited zone in which this works.....
  if(s.re < 0 && q.im == 0 && q.re > 0 && q.re <= 1) {
    // TODO: try the formula with both cos and sin (formula 8 at http://mathworld.wolfram.com/HurwitzZetaFunction.html)
    //       does support complex q. Then use the formulas for Distante neighbors at http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta2/17/01/ShowAll.html
    //       to bring real part of q in range 0-1, and see if that would make it support all complex q for all complex negative s...
    //       or just make q.re > 2 to make hurwitzzeta_hasse_series_ work ...
    return JmatC.hurwitzzeta_cos_series_(s, q);
  }

  if(q.re > 2) {
    return JmatC.hurwitzzeta_hasse_series_(s, q);
  }

  // The smaller, the more iterations needed, hence limited to such number...
  if(q.re > -40 && !q.eqr(0) /* 0 causes infinite recursion*/) {
    //zeta(s, q) = zeta(s, q + m) + SUM_n=0..(m-1) (1/(n+q)^s)
    //if m < 0:
    //zeta(s, a) = SUM_n=0..(m-1) (1/(n+a-m)^s) - zeta(s, a - m)
    // For s.re < 0 and real q, hurwitzzeta_cos_series_ works best, so bring q in range 0-1. Else, make q.re > 2 for hurwitzzeta_hasse_series_.
    var m = (s.re < 0 && q.im == 0) ? Math.ceil(-q.re) : Math.ceil(2 - q.re + 0.5);
    var sum = JmatC.ZERO;
    if(m > 0) {
      for(var n = 0; n < m; n++) {
        var r = q.addr(n).pow(s).inv();
        if(!JmatC.isNaN(r)) sum = sum.add(r);
      }
    } else {
      for(var n = 0; n < -m; n++) {
        var r = q.addr(n).addr(m).pow(s).inv();
        if(!JmatC.isNaN(r)) sum = sum.add(r);
      }
    }
    var h = JmatC.hurwitzzeta(s, q.addr(m));
    return m > 0 ? h.add(sum) : h.sub(sum);
  }

  // (TODO: use the functional equation, that will allow to support negative s for integer and rational q (with small denominator, otherwise it gets too many calculations))

  // TODO: support s.re < 0 for all q
  return jmatc(NaN); // not supported in this region :( e.g.
};

JmatC.beta = function(x, y) {
  // definition: beta(x, y) = gamma(x)*gamma(y) / gamma(x+y)
  //return JmatC.gamma(x).mul(JmatC.gamma(y)).div(JmatC.gamma(x.add(y)));
  //return JmatC.exp(JmatC.loggamma(x).add(JmatC.loggamma(y)).sub(JmatC.loggamma(x.add(y))));

  // For the negative integers (which cause the gamma function to return indeterminate results) for which beta is still defined:
  // beta(x, y), with x < 0 and y > 0, is, by the formula: gamma(y)*gamma(x) / gamma(x + y), rewritten to: gamma(y) / (gamma(x + y) / gamma(x)), for example:
  // beta(x, 2) = 1 / x(x+1),  beta(x, 3) = 2 / x(x+1)(x+2),  beta(x, 4) = 6 / x(x+1)(x+2)(x+3),  beta(x, 5) = 24 / x(x+1)(x+2)(x+3)(x+4), etc...
  // When x is negative integer, x(x+1)(x+2)...(x+y-1) can be rewritten as (-1^y) * |x|(|x|-1)(|x|-2)...(|x|-y-1), which is (-1^y) * gamma(|x|+1) / gamma(|x|-y+1)
  // And those last gamma functions get positive argument, so it works without NaN again, so we get the right answers
  // The JmatC.gammaDiv_ function is used to get similar effect. Use a negative input value as its numerator.

  if(x.re < 50 && x.re > -50 && y.re < 50 && y.re > -50) {
    // gamma rather than loggamma more precise here
    return JmatC.gammaDiv21_(x, y, x.add(y));
  } else {
    return JmatC.exp(JmatC.loggammaDiv21_(x, y, x.add(y)));
  }
};

// Incomplete beta function
// B_x(a, b)
JmatC.incbeta = function(x, a, b) {
  // for x = 1, incbeta(1, a, b) = beta(a, b)
  // otherwise: integrate(0..x, t^(a-1) * (1-t)^(b-1) * dt)

  // TODO: gives incorrect result for incbeta(5, 100, 100). Should be -1.4e127. Gives +1.1e141. Such values are needed by fisher F distribution though.

  if(x.eqr(1)) return JmatC.beta(a, b);

  // METHOD A: series representation (simpler case of hypergeometric series). Probably precise, but does not work for neg a/b, x out of range 0-1 (0.75 for better convergence), ...
  // The summation is z^a * SUM_(n=0..inf) ((1-b)_n / (a + n) * x^n / n!
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  if(JmatC.isPositive(x) && x.re < 0.75 && JmatC.isPositive(a) && JmatC.isPositive(b)) {
    var b1 = JmatC.ONE.sub(b);
    var r;
    var result;
    for(var n = 0; n < 30; n++) {
      if(n == 0) {
        r = jmatc(1);
        result = jmatc(0);
      } else {
        r = r.mul(b1.addr(n - 1));
        r = r.mul(x);
        r = r.divr(n);
      }
      if(r.eqr(0)) break;
      result = result.add(r.div(a.addr(n)));
    }
    return x.pow(a).mul(result);
  }

  // METHOD B: in terms of hypergeometric (computationally more expensive, but precise)
  return x.pow(a).div(a).mul(JmatC.hypergeometric(a, JmatC.ONE.sub(b), a.inc(), x));

  // METHOD C: with the integral definition. However this one is numerically very imprecise
  // The integral is numerically too imprecise, due to the degenerate value at t=0
  // the start should be 0, but it's degenerate there
  /*if(JmatC.isPositive(x) && x.re < 1 && JmatC.isPositive(a) && JmatC.isPositive(b)) {
    return JmatC.integrate(x.divr(3000), x, function(t) {
      var r = t.pow(a.subr(1)).mul(JmatC.ONE.sub(t).pow(b.subr(1)));
      return r; 
    }, 30);
  }*/
};

// Regularized incomplete beta: (incomplete beta / beta)
// I_x(a, b)
JmatC.beta_i = function(x, a, b) {
  if(JmatC.isNegativeIntOrZero(a) && !JmatC.isNegativeIntOrZero(b)) return jmatc(1);
  if(JmatC.isNegativeIntOrZero(b) && !JmatC.isNegativeIntOrZero(a)) return jmatc(0);
  if(JmatC.isNegativeIntOrZero(b.add(a)) && x.eqr(1)) return jmatc(1);
  return JmatC.incbeta(x, a, b).div(JmatC.beta(a, b));
};

// Inverse of regularized incomplete beta: I_x^-1(a, b)
JmatC.beta_i_inv = function(x, a, b) {
  if(!(JmatC.isPositiveOrZero(x) && JmatC.isPositiveOrZero(a) && JmatC.isPositiveOrZero(b))) return jmatc(NaN);

  var bab = JmatC.beta(a, b);
  var x2 = JmatC.ZERO;
  var a2 = JmatC.ZERO;
  var b2 = JmatC.ONE;

  while(!JmatC.near(a2, b2, 1e-6)) {
    x2 = a2.add(b2).divr(2);
    if(JmatC.incbeta(x2, a, b).div(bab).re > x.re) b2 = x2;
    else a2 = x2;
  }
  return x2;
};

// Confluent hypergeometric limit function 0F1(; a; z) --> yep, that's the notation of it
JmatC.hypergeometric0F1 = function(a, z) {
  // The summation is SUM_(n=0..inf) z^n / ((a)_n * n!)
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  // This loop converges very nicely and result is quite accurate for all complex input values (unlike for 2F1)
  var r = JmatC.ONE;
  var result = JmatC.ZERO;
  for(var n = 0; n < 30; n++) {
    if (n > 0) {
      if(!r.eqr(0)) r = r.div(a.addr(n - 1)); // the pochammer. The if avoids 0/0
      r = r.mul(z); // the z^n
      r = r.divr(n); // the n!
    }
    if(r.eqr(0)) break;
    result = result.add(r);
  }

  return result;
};

// Confluent hypergeometric function of the first kind 1F1(a; b; z), a.k.a. Kummer's function M (approximation)
// Try similar online: http://keisan.casio.com/exec/system/1349143651
JmatC.hypergeometric1F1 = function(a, b, z) {
  // The summation is SUM_(n=0..inf) ((a)_n / (b)_n) * z^n / n!
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  // This loop converges very nicely (unlike for 2F1)
  var r = a.div(b).mul(z).divr(1);
  var result = jmatc(1);
  for(var n = 1; n < 30; n++) {
    if (n > 1) {
      r = r.mul(a.addr(n - 1));
      if(!r.eqr(0)) r = r.div(b.addr(n - 1));
      r = r.mul(z);
      r = r.divr(n);
    }
    if(r.eqr(0)) break;
    result = result.add(r);
  }

  return result;
};

// Hypergeometric series 2F1(a, b; c; z), aka Gauss hypergeometric function (approximation)
// Simply called "hypergeometric" instead of "hypergeometric2F1"
// Try similar online: http://keisan.casio.com/exec/system/1349143084, or wolframalpha, HyperGeometric2F1[1.5, -0.5, 2.5, 0.25+i]
// I debugged this function to find all the areas in the complex or real input plane where it goes wrong, by 2D or complexdomain plotting this, and beta_i which uses this
JmatC.hypergeometric = function(a, b, c, z) {
  if(JmatC.absr(z) > 1.0001 /*not 1 to avoid infinite loop if abs 1/z is also > 1 due to numeric problems, e.g. for 0.9726962457337884+i0.23208191126279865*/) {
    // The series converges only for |z| < 1. But there are some linear transformations
    // that can convert it to a different z. There are conditions though, and some are
    // more complex than others (e.g. requiring gamma functions).
    // TODO: with only those two supported transformations below, there is probably a lot wrong,
    //       for example, the points at exp(pi*i/3) and exp(-pi*i/3) are known to be not covered by this
    //       To fix, do as in the paper "Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Poschl-Teller-Ginocchio potential wave functions"
    
    // Linear transformations to bring |z| in value < 1
    var z2 = z.div(z.dec());
    if(JmatC.absr(z2) < 0.75) { // the if can in theory do "< 1", but then a value like z=0.3+4i makes z2 very close to 1, and it has bad numeric results for such values
      // z / (z - 1)
      return JmatC.ONE.sub(z).pow(a.neg()).mul(JmatC.hypergeometric(a, c.sub(b), c, z2));
    } else {
      // 1 / z
      var quirckyGammaDiv22_ = function(a, b, c, d) {
        // For some cases where a, b, c - a or c - b are negative integers, the formula doesn't work and requires a rather complicated other formula for the solution.
        // For now, I temporarily instead twiddle the parameters a bit. TODO: that is evil, do it properly (but it is surprisingly somewhat accurate though... well, not for everything)
        var result = JmatC.gammaDiv22_(a, b, c, d);
        if(JmatC.isNaN(result)) {
          if(JmatC.isNegativeIntOrZero(a)) a = a.addr(1e-5);
          if(JmatC.isNegativeIntOrZero(b)) b = b.addr(1e-5);
          if(JmatC.isNegativeIntOrZero(c)) c = c.addr(1e-5);
          if(JmatC.isNegativeIntOrZero(d)) d = d.addr(1e-5);
          result = JmatC.gammaDiv22_(a, b, c, d);
        }
        return result;
      };

      var zi = z.inv();
      var za = z.neg().pow(a.neg());
      var zb = z.neg().pow(b.neg());
      var ga = quirckyGammaDiv22_(c, b.sub(a), b, c.sub(a));
      var gb = quirckyGammaDiv22_(c, a.sub(b), a, c.sub(b));
      var fa = JmatC.hypergeometric(a, JmatC.ONE.sub(c).add(a), JmatC.ONE.sub(b).add(a), zi);
      var fb = JmatC.hypergeometric(b, JmatC.ONE.sub(c).add(b), JmatC.ONE.sub(a).add(b), zi);
      var va = ga.mul(za).mul(fa);
      var vb = gb.mul(zb).mul(fb);
      return va.add(vb);
    }
  }

  var z2 = z.div(z.dec());
  if(JmatC.absr(z2) < JmatC.absr(z)) {
    // Same z / (z - 1) transform as above. Reason for doing this: the summation below converges faster for smaller absolute values of z. Without this, for e.g. z = -0.75 and c < -3, it converges only after hundreds of steps.
    // TODO: in fact it's almost always possible to make |z| < 0.5 with the linear transformations. Use them better.
    return JmatC.ONE.sub(z).pow(a.neg()).mul(JmatC.hypergeometric(a, c.sub(b), c, z2));
  }

  // TODO: avoid needing more than 30 iterations by always getting |z| < 0.5 with more clever use of transformations
  var num_iterations = 30;
  if(JmatC.absr(z) > 0.5) num_iterations = 60;
  if(JmatC.absr(z) > 0.75) num_iterations = 100;

  // The summation definition of the series, converges because |z| < 1

  // The summation is SUM_(n=0..inf) ((a)_n * (b)_n / (c)_n) * z^n / n!
  // Where (q)_n is the rising Pochhammer symbol x(x+1)*...*(x+n-1)
  // The variable r below gets updated every step with next factor of Pochhammer symbol, power and factorial values as the loop goes.
  var r = a.mul(b).div(c).mul(z).divr(1);
  var result = jmatc(1);
  for(var n = 1; n < num_iterations; n++) {
    if (n > 1) {
      r = r.mul(a.addr(n - 1));
      r = r.mul(b.addr(n - 1));
      if(!r.eqr(0)) r = r.div(c.addr(n - 1)); // The if is there so that if a==c or b==c and c is neg integer, we don't get 0/0 here but just 0
      r = r.mul(z);
      r = r.divr(n);
    }
    if(r.eqr(0)) break; // If a or b are negative integer, the loop will terminate early
    if(JmatC.near(r, JmatC.ZERO, 1e-15)) break; // In fact why not even break out if it's already converged very near zero
    result = result.add(r);
  }
  return result;

  // Alternative: integration formula (doesn't work well)
  /*var g = JmatC.gammaDiv12_(c, b, c.sub(b));
  var i = JmatC.integrate(jmatc(0), jmatc(1), function(t) {
    var n = t.pow(b.dec()).mul(JmatC.ONE.sub(t).pow(c.sub(b).dec()));
    var d = JmatC.ONE.sub(z.mul(t)).pow(a);
    return n.mul(d);
  }, 30);
  return g.mul(i);*/
};

// TODO: generalized hypergeometric

JmatC.PIPI6_ = jmatc(Math.PI * Math.PI / 6);

//Dilogarithm: Li_2(z)
JmatC.dilog = function(z) {
  if(z.eqr(0)) return JmatC.ZERO;
  if(z.eqr(1)) return JmatC.PIPI6_;
  if(z.eqr(+Infinity)) return jmatc(-Infinity);
  if(z.eqr(-Infinity)) return jmatc(-Infinity);

  // Do the series expansion for |z| < 1
  var summation = function(z) {
    var result = JmatC.ZERO;
    var zz = z;
    var N = a <= 0.5 ? 30 : 50;
    for(var i = 1; i < N; i++) {
      var r = zz.divr(i * i);
      result = result.add(r);
      if(JmatC.near(r, JmatC.ZERO, 1e-15)) break;
      zz = zz.mul(z);
    }
    return result;
  };

  var a = z.absr();
  if(a < 0.5) {
    // The series only converges for |z| < 1 (and only fast enough if < 0.75).
    // For other cases, identities in the various other ifs below are used to transform it to this size.
    return summation(z);
  }

  // The identity for 1/z
  if(1 / a < 0.5) {
    var d = summation(z.inv()).neg();
    var l = JmatC.log(z.neg());
    return d.sub(l.mul(l).divr(2)).sub(JmatC.PIPI6_);
  }

  // The above was done to avoid all the stuff below in the obvious good above cases. Now consider each option, including the above two and other identities, again for worse |z|

  // Three more identities
  var z1 = JmatC.ONE.sub(z);
  var a1 = z1.absr();
  var z2 = JmatC.ONE.sub(z).inv();
  var a2 = z2.absr();
  var z3 = z.div(z.dec());
  var a3 = z3.absr();

  var best = a;
  var bestindex = 0; //0 for z, 1 for 1/z, 2 for 1-z, 3 for 1/(1-z), 4 for z/(z-1)
  if(1/a < best) {
    best = 1/a;
    bestindex = 1;
  }
  if(a1 < best) {
    best = a1;
    bestindex = 2;
  }
  if(a2 < best) {
    best = a2;
    bestindex = 3;
  }
  if(a3 < best) {
    best = a3;
    bestindex = 4;
  }

  if(best < 0.8) {
    if(bestindex == 0) {
      return summation(z);
    }
    if(bestindex == 1) {
      var d = summation(z.inv()).neg();
      var l = JmatC.log(z.neg());
      return d.sub(l.mul(l).divr(2)).sub(JmatC.PIPI6_);
    }
    if(bestindex == 2) {
      var d = summation(z1).neg();
      var l1 = JmatC.log(z1);
      var l2 = JmatC.log(z);
      return d.sub(l1.mul(l2)).add(JmatC.PIPI6_);
    }
    if(bestindex == 3) {
      var d = summation(z2);
      var l1 = JmatC.log(JmatC.ONE.sub(z));
      var l2 = JmatC.log(z.neg());
      return d.add(l1.mul(l1).divr(2)).sub(l2.mul(l1)).sub(JmatC.PIPI6_);
    }
    if(bestindex == 4) {
      var d = summation(z3).neg();
      var l = JmatC.log(JmatC.ONE.sub(z));
      return d.sub(l.mul(l).divr(2));
    }
  }
  
  // Near values 0.5+i and 0.5-i the above is not working.
  // The following identity is: Li_2(z) = Li_2(z^2)/2 - Li_2(-z)
  // That is, the square relationship
  // This gives z's which are eventually not in the problem region, and thus we will not end up with infinite recursive call of this function.

  return JmatC.dilog(z.mul(z)).divr(2).sub(JmatC.dilog(z.neg()));
};


//Trilogarithm: Li_3(z)
JmatC.trilog = function(z) {
  if(z.eqr(0)) return JmatC.ZERO;
  if(z.eqr(1)) return JmatC.APERY;
  if(z.eqr(+Infinity)) return jmatc(-Infinity);
  if(z.eqr(-Infinity)) return jmatc(-Infinity);

  // Do the series expansion for |z| < 1
  var summation = function(z) {
    var result = JmatC.ZERO;
    var zz = z;
    var N = a < 0.5 ? 30 : 50;
    for(var i = 1; i < N; i++) {
      var r = zz.divr(i*i*i);
      result = result.add(r);
      if(JmatC.near(r, JmatC.ZERO, 1e-15)) break;
      zz = zz.mul(z);
    }
    return result;
  };

  var a = z.absr();
  if(a < 0.75) {
    // The series only converges for |z| < 1 (and only fast enough if < 0.75).
    // For other cases, identities in the various other ifs below are used to transform it to this size.
    return summation(z);
  }

  // The identity for 1/z
  if(1 / a < 0.75) {
    var d = summation(z.inv());
    var l = JmatC.log(z.neg());
    return d.sub(l.mul(l).mul(l).divr(6)).sub(JmatC.PIPI6_.mul(l));
  }

  // TODO: better implementation for this case
  return JmatC.polylog_integral_(jmatc(3), z);
};

// Bernoulli numbers B_m. Only supports a hardcoded first few dozens.
JmatC.bernoulli = [
    1, -1/2, 1/6, 0, //0-3
    -1/30, 0, 1/42, 0, //4-7
    -1/30, 0, 5/66, 0, //8-11
    -691/2730, 0, 7/6, 0, //12-15
    -3617/510, 0, 43867/798, 0, //16-19
    -174611/330, 0, 854513/138, 0, //20-23
    -23749461029/2730, 0, 8615841276005/6, 0, //24-27
    -7709321041217/870, 0, 2577687858367/14322, 0, //28-31
];

// Stieltjes constants. Only supports a hardcoded first few dozens.
JmatC.stieltjes = [
  0.577215664901532861, -0.0728158454836767249, -0.00969036319287231848, 0.00205383442030334587,
  0.00232537006546730006, 0.000793323817301062702, -0.000238769345430199610, -0.000527289567057751046,
  -0.000352123353803039509, -0.0000343947744180880482, 0.000205332814909064795, 0.000270184439543903527,
  0.000167272912105140193, -0.0000274638066037601589, -0.000209209262059299946, -0.000283468655320241447,
  -0.000199696858308969775, 0.0000262770371099183367, 0.000307368408149252827, 0.000503605453047355629,
  0.000466343561511559449, 0.000104437769756000116, -0.000541599582203997702, -0.00124396209040824578,
  -0.00158851127890356156, -0.00107459195273848882, 0.000656803518637154432, 0.00347783691361853821,
  0.00640006853170062946, 0.00737115177047223913, 0.00355772885557316095, -0.00751332599781522893,
];

// Stieltjes constants times (-1)^n / n!, this is for calculating riemann zeta. Only supports a hardcoded first few dozens (0-31).
JmatC.stieltjes_zeta = [
  0.577215664901532861, 0.0728158454836767249, -0.00484518159643615924, -0.000342305736717224311,
  0.0000968904193944708357, -6.61103181084218918e-6, -3.31624090875277236e-7, 1.04620945844791874e-7,
  -8.73321810027379736e-9, 9.47827778276235895e-11, 5.65842192760870797e-11, -6.76868986351369666e-12,
  3.49211593667203185e-13, 4.41042474175775338e-15, -2.39978622177099918e-15, 2.16773122007268285e-16,
  -9.54446607636696517e-18, -7.38767666053863650e-20, 4.80085078248806523e-20, -4.13995673771330564e-21,
  1.91682015939912339e-22, -2.04415431222621661e-24, -4.81849850110735344e-25, 4.81185705151256648e-26,
  -2.56026331031881494e-27, 6.92784089530466712e-29, 1.62860755048558674e-30, -3.19393756115325558e-31,
  2.09915158936342553e-32, -8.33674529544144048e-34, 1.34125937721921867e-35, 9.13714389129817200e-37,
];

// This is really a last resort, very inaccurate and slow
// Welcome to the land of randomness and bad results.
// Some cases are really good, others horrible.
// Unfortunately almost every integral has problems on the positive real axis of z.
// Since there is an awesome solution for all s with s.re < 0, only the formulas for s.re > 0 below are actually u sed
JmatC.polylog_integral_ = function(s, z) {
  // To test these, try e.g.:
  // complexDomainPlot(function(z){return JmatC.polylog(jmatc(15, 0.5), z);}, 2, 1);
  // complexDomainPlot(function(z){return JmatC.polylog(jmatc(0.5, 15), z);}, 2, 1);

  if(s.re > 1 && Math.abs(s.im) < s.re && Math.abs(z.argr()) > 0.1) {
    // Approximate with an integral representation (only works for real s.re, and has some serious problems with periodic things appearing near the positive real axis of z)
    // In practice, only works for s.re > 1 (!!) and s.im = 0, or very small |s.im| and s.re > 0
    var g = JmatC.gamma(s);
    var r = JmatC.integrate(jmatc(0), jmatc(20), function(t) {
      var result = t.pow(s.dec()).div(JmatC.exp(t).sub(z));
      if(JmatC.isNaN(result)) result = JmatC.ZERO;
      return result;
    }, 100);
    return z.div(g).mul(r);
  } else if(JmatC.isNegative(s) && Math.abs(z.argr()) > 0.1) {
    // Approximate with an integral representation (only works for s.re < 0
    var lzm = JmatC.log(z.neg());
    var r = JmatC.integrate(jmatc(0), jmatc(20), function(t) {
      var ta = t.pow(s.neg());
      var tb = JmatC.sin(s.mulr(Math.PI/2).sub(t.mul(lzm)));
      var na = JmatC.sinh(t.mulr(Math.PI));
      var result = ta.mul(tb).div(na);
      if(JmatC.isNaN(result)) result = JmatC.ZERO;
      return result;
    }, 100);
    return r;
  } else if(z.im <= 0 || (s.re > 0 && Math.abs(s.im) < s.re)) {
    // Integral formula for polylog that works for all complex s and z, except for positive real z if s.re < 0
    // Is from 0 to infinity, but seems to work well from 0-10 with only 100 steps (at least in the areas where this is used)
    // While in theory it works for all z and all s (except z near positive axis if s.re < 0), in practice it seems to work only for z.im <= 0 and s.re > 0 and |s.im| << s.re
    // --> I tried with z-plot for s=-15+0.5i, s=0.5+15i, and other variations like that
    var lzm = JmatC.log(z.neg());
    var f = function(t) {
      var ta = JmatC.sin(s.mul(JmatC.atan(t)).sub(t.mul(lzm)));
      var na = (JmatC.ONE.add(t.mul(t))).pow(s.divr(2));
      var nb = JmatC.sinh(t.mulr(Math.PI));
      var result = ta.div(na).div(nb);
      if(JmatC.isNaN(result)) result = JmatC.ZERO;
      return result;
    };
    var r = JmatC.ZERO;
    r = r.add(JmatC.integrate(jmatc(0), jmatc(5), f, 50));
    r = r.add(JmatC.integrate(jmatc(5), jmatc(20), f, 20));
    r = r.add(JmatC.integrate(jmatc(20), jmatc(100), f, 10));
    return z.mulr(0.5).add(z.mul(r));
  } else if(z.im > 0) { // Because something is broken for z.im <= 0. TODO: find out what
    // Similar integral, but with upper incomplete gamma.
    // It is theoretically better than the above because it supports even its last edge case.
    // But it seems broken. It does not work for z.im <= 0... At least it gets 
    var lz = JmatC.log(z);
    var f = function(t) {
      var ta = JmatC.sin(s.mul(JmatC.atan(t)).sub(t.mul(lz)));
      var na = (JmatC.ONE.add(t.mul(t))).pow(s.divrdivr(2));
      var nb = JmatC.exp(t.mulr(2*Math.PI)).sub(1);
      var result = ta.div(na).div(nb);
      if(JmatC.isNaN(result)) result = JmatC.ZERO;
      return result;
    };
    var r = JmatC.ZERO;
    r = r.add(JmatC.integrate(jmatc(0), jmatc(5), f, 50));
    r = r.add(JmatC.integrate(jmatc(5), jmatc(20), f, 20));
    r = r.add(JmatC.integrate(jmatc(20), jmatc(100), f, 10));
    var g = JmatC.incgamma_upper(JmatC.ONE.sub(s), lz.neg());
    var l = lz.neg().pow(JmatC.ONE.sub(s));
    return z.mulr(0.5).add(g.div(l)).add(z.mulr(2).mul(r));
  } else {
    return jmatc(NaN); //there is nothing we can do... :(
  }
};

// Borwein algorithm. Converges theoretically only for |z^2/(z-1)| < 3.7.
// In practice this zone of convergence is much smaller due to numerical problems. For extreme values of s (e.g. -20), it is as small as 0.5
// This algorithm is actually intended for arbitrary precision libraries, not for floating point.
// I'm keeping it here, because for s.re > -3, and with the adaptive n, it works better than most other polylog code here and is quite accurate... The function "polylog_borwein_ok_" returns whether it's usable for the given values.
// Uses the following formula:
// Li_s(z) = SUM_k=1..n(z^k / k^s) + 1/(1-z)^n * SUM_k=n+1..2n(z^k / k^s * SUM_j=0..2n-k((-z)^j * binomial(n, j))) + error(s, z)
// where the error term is negligable in the zone of convergence
// NOTE: changing z to -1 makes this the borwein algorithm for riemann zeta.
JmatC.polylog_borwein_ = function(s, z) {
  var kidney_radius = z.mul(z).div(z.dec()).absr(); //radius of the kidney shaped region of convergence. Theoretically works for < 4, in practice with double precision only for < 3 and then only if not too much n!
  if(kidney_radius >= 3.7) return jmatc(NaN); //yeah right... it sucks way before this

  // number of loops. NOTE: higher is better for arbitrary precision library (with 31 being perfect), but results in random garbage with floating point. So it is limited here instead.
  var n = Math.floor(Math.min(31, 16 / kidney_radius));


  var binomial_cache = [];
  binomial_cache[0] = JmatC.ONE;
  var bin_sum_term = JmatC.ONE; // Binomial summation term

  var zz = JmatC.ONE;
  var result0 = JmatC.ZERO;

  var barray = [];

  // First sum, and preparation for binomial sum
  for(var k = 1; k <= n; k++) {
    zz = zz.mul(z); // numerically critical...
    var r = zz.div(jmatc(k).pow(s));
    result0 = result0.add(r);

    // Binomial sum
    var bin = JmatC.binomial(jmatc(n), jmatc(k));
    var b = bin.mul(zz); // numerically critical...
    if(k % 2) b = b.neg(); // numerically critical...

    bin_sum_term = bin_sum_term.add(b); // numerically critical...

    // Store binomial sums in array for later reference
    binomial_cache[k] = bin_sum_term;
  }

  // Second sum, using the binomial sum preparation from above
  var result1 = JmatC.ZERO;
  for(k = n + 1; k <= 2 * n; k++) {
    zz = zz.mul(z);
    var r = zz.div(jmatc(k).pow(s));

    bin_sum_term = binomial_cache[2*n-k];
    r = r.mul(bin_sum_term);

    result1 = result1.add(r);
  }

  // 1/(z-1)^n
  var invpow = z.dec().inv().powr(n);

  result1 = result1.mul(invpow);
  var result;
  if(n % 2) result = result0.sub(result1);
  else result = result0.add(result1);

  return result;
};

// Is it ok to use borwein method?
JmatC.polylog_borwein_ok_ = function(s, z) {
  var kidney_radius = z.mul(z).div(z.dec()).absr();
  return (s.re >= -5 && z.im == 0 && kidney_radius <= 1.5) || (s.re >= 0 && Math.abs(s.im) < s.re && kidney_radius <=2) || z.absr() <= 0.5;
};

JmatC.polylog_residue_ = function(s, z) {
  // Sum of residues: LI_s(e^mu) = gamma(1-s) * SUM_k=-oo..oo (2*k*pi*i - mu)^(s-1)
  // Theoretically holds for Re(s) < 0 and e^mu != 1
  // This seems to work pretty well in practice! Even for large z. It works for all z except 1, for all s with negative real part (even things like JmatC.polylog(jmatc(-15,0.5), jmatc(-100,1)) match wolfram alpha polylog(-15+0.5i, -100+i)!!)

  //if(z.re < 0) return square(s, z); //numerically not stable

  var result_is_real = false;
  var origz = z;

  // The only thing this doesn't work well for is for real z with |z| > 2. Twiddle with slightly negative imaginary part to compensate. This because wolfram alpha matches the values for negative imaginary part there as well. NOTE: this is a branch cut line so of course it's tricky here
  if(s.im == 0 && z.im == 0) {
    z = z.add(JmatC.newi(-0.0001));
    if(JmatC.isNegativeInt(s) || z.re < 0) result_is_real = true; //if s is not a negative integer and z > 0, then despite s and z being real, the result can be complex. However in other cases, it is real, and the twiddling ruins that, so compensate it at the end.
  }

  var mu = JmatC.log(z);
  var g = JmatC.gamma(JmatC.ONE.sub(s));

  var r = JmatC.ZERO;
  for(var k = 0; k < 30; k++) {
    // The Wikipedia formula has sum from -Infinity to +Infinity. So there are two terms, a positive and a negative.
    // However, when I try it, for z.re > 0, that makes if off by a factor of two, and for z.re < 0 wrong in some different way. The first term works for z.im > 0, the second for z.im <= 0, I think.
    // So: negative-k term DISABLED!!!
    // TODO: figure out why that is. Because some things do still go wrong with this...
    if(origz.im > 0) r = r.add((JmatC.newi(2 * k * Math.PI).sub(mu)).pow(s.dec()));
    else r = r.add((JmatC.newi(-2 * k * Math.PI).sub(mu)).pow(s.dec()));
  }

  var result = g.mul(r);
  if(result_is_real && result.im != 0) result = jmatc(result.re);
  return result;
};

//Polylogarithm: Li_s(z)
JmatC.polylog = function(s, z) {
  if(JmatC.isInt(s)) {
    if(s.eqr(0)) {
      return z.div(JmatC.ONE.sub(z));
    } else if(s.eqr(1)) {
      return JmatC.log(JmatC.ONE.sub(z)).neg();
    } else if(s.eqr(2)) {
      return JmatC.dilog(z);
    } else if(s.eqr(3)) {
      return JmatC.trilog(z);
    } else if(s.eqr(-1)) {
      var z1 = JmatC.ONE.sub(z);
      return z.div(z1).div(z1);
    } else if(s.eqr(-2)) {
      var z1 = JmatC.ONE.sub(z);
      return z.mul(z.inc()).div(z1).div(z1).div(z1);
    } else if(s.eqr(-3)) {
      var z1 = JmatC.ONE.sub(z);
      var zz1 = z1.mul(z1);
      return z.mul(JmatC.ONE.add(z.mulr(4)).add(z.mul(z))).div(zz1).div(zz1);
    } else if(s.eqr(-4)) {
      var z1 = JmatC.ONE.sub(z);
      var zz1 = z1.mul(z1);
      return z.mul(z.inc()).mul(JmatC.ONE.add(z.mulr(10)).add(z.mul(z))).div(zz1).div(zz1).div(z1);
    }
  }

  // TODO: there are still problems with:
  // 1) pure imaginary s
  // 2) z near 1 (no matter what s)
  // 3) despite its goodness, there's something wrong with the Sum of residues formula. Normally both sum terms should be used, but some zones seem to require only one. Make it work both for complexDomainPlot(function(z){return JmatC.polylog(jmatc(-5, 15), z);}, 2, 1); (where using the other one, makes it black) and plot2D(JmatC.polylog, 10, 1) (where the complex values in lower left quadrant should be real), and yet some other zones would make it double its value if both are enabled...
  // 4) for s.re > 0, it is a lot less accurate than for s.re < 0, because for s.re < 0 a good series can be used, while for s.re > 0, the inaccurate integrals are used and no inversion formula for is available here

  if(s.re > 1 && z.eqr(0)) return jmatc(0);
  if(s.re > 1 && z.eqr(1)) return JmatC.zeta(s); // It's the riemann zeta function for z=1, but only if s.re > 1
  /*if(z.eqr(-1)) return JmatC.eta(s).neg();*/

  // The duplication formula (square relationship):
  // Li_s(-z) + Li_s(z) = 2^(1-s)*Li_s(z^2)
  // or:
  // Li_s(z) = 2^(1-s)*Li_s(z^2) - Li_s(-z)
  // Li_s(-z) = 2^(1-s)*Li_s(z^2) - Li_s(z)
  var square = function(s, z) {
    var a = JmatC.polylog(s, z.mul(z));
    var b = JmatC.polylog(s, z.neg());
    var c = JmatC.TWO.pow(JmatC.ONE.sub(s));
    return c.mul(a).sub(b);
  };

  var a = z.absr();

  if((a <= 0.5 && s.re >= -10) || (a < 0.75 && s.re >= -2) || (a < 0.9 && s.re > -1)) {
    // METHOD A: the series definition SUM_k=1..oo z^k / k^s. But only valid for |z| < 1 and converges slowly even then
    // In practice, only useful if |z| < 0.5 and s.re not too small, or higher |z| for larger s.re. If it works, this is one of the most accurate representations...
    var result = JmatC.ZERO;
    var zz = z;
    var N = a <= 0.5 ? 30 : 50;
    for(var i = 1; i < N; i++) {
      var r = zz.div(jmatc(i).pow(s));
      result = result.add(r);
      if(JmatC.near(r, JmatC.ZERO, 1e-15)) break;
      zz = zz.mul(z);
    }
    return result;
  } else if(JmatC.polylog_borwein_ok_(s, z)) {
    return JmatC.polylog_borwein_(s, z);
  } else if(s.re < 0) {
    return JmatC.polylog_residue_(s, z);
  } else if(JmatC.isNegativeInt(s) && a > 1) {
    // Inversion formula. Since the "sum of residue" above works so fantastic and is already for all negative s, this code probably never gets called anymore.
    var sign = JmatC.isOdd(s) ? 1 : -1; // (-1)^(s-1)
    return JmatC.polylog(s, z.inv()).mulr(sign);
  } else {
    // Last resort: the not very well working integrals :(
    // Polylog is supposed to be very interesting at s=0.5+t*i. But no code here is precise enough for that.
    return JmatC.polylog_integral_(s, z);
  }
  return jmatc(NaN);
};

//Jacobi theta1 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
JmatC.theta1 = function(z, q) {
  // theta1 is an odd function, and the series below has numerical problems when imaginary part is large, due to values of complex sine becoming huge there --> actually I don't think that helps, the same problem with sine is there for large negative imaginary value...
  // therefore, mirror it
  if(z.im > 1) return JmatC.theta1(z.neg(), q).neg();

  var result = jmatc(0);
  var s = jmatc(1);
  for(var n = 0; n < 20; n++) {
    var a = q.powr((n + 0.5) * (n + 0.5));
    var b = JmatC.sin(z.mulr(2 * n + 1));
    result = result.add(s.mul(a).mul(b));
    /*var a = JmatC.log(q.powr((n + 0.5) * (n + 0.5)));
    var b = JmatC.logsin(z.mulr(2 * n + 1));
    result = result.add(s.mul(JmatC.exp((a).add(b))));*/
    s = s.neg();
  }
  return result.mulr(2);
};

//Jacobi theta2 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
JmatC.theta2 = function(z, q) {
  // theta2 is an even function, and the series below has numerical problems when imaginary part is large, so mirror it --> actually I don't think that helps, the same problem with sine is there for large negative imaginary value...
  if(z.im > 1) return JmatC.theta2(z.neg(), q);

  var result = jmatc(0);
  for(var n = 0; n < 20; n++) {
    var a = q.powr((n + 0.5) * (n + 0.5));
    var b = JmatC.cos(z.mulr(2 * n + 1));
    result = result.add(a.mul(b));
  }
  return result.mulr(2);
};

//Jacobi theta3 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
JmatC.theta3 = function(z, q) {
  // theta3 is an even function, and the series below has numerical problems when imaginary part is large, so mirror it
  if(z.im > 1) return JmatC.theta3(z.neg(), q);

  var result = jmatc(0);
  for(var n = 0; n < 20; n++) {
    var a = q.powr(n * n);
    var b = JmatC.cos(z.mulr(2 * n));
    result = result.add(a.mul(b));
  }
  return result.mulr(2).addr(1);
};

//Jacobi theta4 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
JmatC.theta4 = function(z, q) {
  // theta4 is an even function, and the series below has numerical problems when imaginary part is large, so mirror it
  if(z.im > 1) return JmatC.theta4(z.neg(), q);

  var result = jmatc(0);
  var s = jmatc(1);
  for(var n = 0; n < 20; n++) {
    var a = q.powr(n * n);
    var b = JmatC.cos(z.mulr(2 * n));
    result = result.add(s.mul(a).mul(b));
    s = s.neg();
  }
  return result.mulr(2).addr(1);
};

// TODO: exponential integral, sine integral, cosine integral, logarithmic integral, elliptic integrals

// Returns 0 or 1
// The reason this function is here as well as in JmatR, while the other prime functions are only in JmatR, is that this one needs to look at the imaginary part and say it's not prime if it's not zero
JmatC.isPrime = function(value) {
  if(!JmatC.isReal(value)) return 0; //complex numbers are not prime
  return JmatR.isPrime(value.re);
};

// returns numerator and denominator of fraction
// max = max value for denominator (a real JS number)
JmatC.decompose = function(x, max) {
  if (Math.abs(x.re) >= Math.abs(x.im)) {
    var nd = JmatR.decompose(x.re, max);
    var im = Math.round(x.im * nd[1]);
    return [jmatc(nd[0], im), jmatc(nd[1])];
  } else {
    var nd = JmatR.decompose(x.im, max);
    var re = Math.round(x.re * nd[1]);
    return [jmatc(re, nd[0]), jmatc(nd[1])];
  }
};

// n! / (n-p)!
JmatC.permutation = function(n, p) {
  // gammaDiv_ is already optimized for integers near each other etc...
  return JmatC.gammaDiv_(n.inc(), n.sub(p).inc());
};

//Binomial coefficient, aka combination(s). Number of rows of p elements that can be made out of n elements, where order doesn't matter.
//    ( n )
//    ( p )
// n! / (p! * (n-p)!)
JmatC.binomial = function(n, p) {
  if(JmatC.isPositiveIntOrZero(n) && JmatC.isPositiveIntOrZero(p) && p.re <= n.re && n.re < 30) return jmatc(JmatR.pascal_triangle(n.re, p.re));
  
  // gammaDiv_ is already optimized for integers near each other etc...
  var result = JmatC.gammaDiv12_(n.inc(), p.inc(), n.sub(p).inc());
  // Round to integer if large result, it sometimes gets numerically a bit off.
  if(result.re > 100 && JmatC.isPositiveInt(n) && JmatC.isPositiveInt(p) && n.re > p.re) result = JmatC.round(result);
  return result;
};

//Stirling number of the second kind
//    { n }
//    { k }
// 1/k! * SUM_j=0..k((-1)^(k-j) * combination(k, j) * j^n)
JmatC.stirling2 = function(n, k) {
  if(!JmatC.isInt(k)) return jmatc(NaN); // only defined for integer k

  var result = JmatC.ZERO;
  var sign = JmatR.isOdd(k.re) ? -1 : 1;
  var j = jmatc(0);
  for(j.re = 0; j.re <= k.re; j.re++) {
    result = result.add(JmatC.binomial(k, j).mul(j.pow(n)).mulr(sign));
    sign *= -1;
  }
  return result.div(JmatC.factorial(k));
};

// manhattan distance of complex numbers, returned as a real number (JS float)
JmatC.manhattan = function(a, b) {
  return Math.max(Math.abs(a.re - b.re), Math.abs(a.im - b.im));
};

JmatC.near = function(x, y, precision) {
  //return JmatC.manhattan(x, y) <= precision;
  // Manhatton NOT used, because then this function returns false for equal infinities
  return x.re - precision <= y.re && x.re + precision >= y.re && x.im - precision <= y.im && x.im + precision >= y.im;
};

// Lambertw for branch (0 = principal branch Wp, -1 is also common (Wm))
// Branch is real integer, z is JmatC object (complex)
JmatC.lambertwb = function(branch, z) {
  if(JmatC.isReal(z) && z.re > -0.36 /*~ -1/e*/ && branch == 0) return jmatc(JmatR.lambertw(z));

  if(!JmatR.isInt(branch)) return jmatc(NaN);


  // Known special values
  if(JmatC.isNaN(z)) return NaN;
  if(JmatC.isInf(z)) return jmatc(Infinity); // any complex infinity gives positive infinity
  if(branch == 0 && z.re == 0 && z.im == 0) return jmatc(0);
  if(branch != 0 && z.re == 0 && z.im == 0) return jmatc(-Infinity); //at all other branch than the principal one, it's -infinity at 0

  /*
  Choosing a good starting value is important. jmatc(0) as starting value works
  most of the time, but does not work at some regions in the negative complex domain,
  e.g. around 5.4+0.1i, 5.5+0.1i, ... and that can be seen as mandelbrot-fractal-like
  circles around those regions in the complex domain plot.
  */
  var w = JmatC.log(z).add(jmatc(0, branch * Math.PI * 2));
  if(branch == 0 && z.absr() < 1.2 /*supposed to be 1/Math.E, but I still see problems beyond that in the complex domain plot...*/) {
    w = JmatC.sqrt(z.mulr(5.43656365691809047).addr(2)).add(jmatc(-1, branch * Math.PI * 2));
  }
  if(branch != 0 && z.im == 0) z.im += 1e-14; // Give it small imaginary part, otherwise it never gets there

  var num = 36;
  for(var i = 0; i < num; i++) {
    var ew = JmatC.exp(w);
    var wew = w.mul(ew);
    var t = wew.sub(z);
    var a = ew.mul(w.addr(1));
    var b = w.addr(2).mul(t).div(w.mulr(2).addr(2));
    w = w.sub(t.div(a.sub(b)));

    var ltest = JmatC.log(z.div(w)); //for testing if near (z = w*exp(w) OR ln(z/w) = w)
    if(JmatC.near(ltest, w, 1e-16) || JmatC.near(wew, z, 1e-16)) break;
    if(i + 1 == num && !(JmatC.near(ltest, w, 1) || JmatC.near(wew, z, 1))) return jmatc(NaN); // iteration could not finish and too far from result
  }

  // Remove numeric tiny imaginary part that appeared in error
  if(z.im == 0 && z.re >= 0) w.im = 0;

  return w;
};

// Principal branch of Lambert's W function: Wp, inverse (not reciprocal) of exp(x) * x
JmatC.lambertw = function(z) {
  return JmatC.lambertwb(0, z);
};

// Negative branch of Lambert's W function: Wm, inverse (not reciprocal) of exp(x) * x
JmatC.lambertwm = function(z) {
  // TODO: wrong. Look at the real plot. Fix this! Jmat.plotReal(JmatC.lambertwm)
  return JmatC.lambertwb(-1, z);
};

// Tetration
// Returns experimental (not mathematically correct) results unless z is a positive integer or Infinity
JmatC.tetration = function(a, z) {
  if(JmatC.isPositive(a) && JmatC.isPositiveInt(z) && z.re != Infinity) {
    return jmatc(JmatR.tetration(a.re, z.re));
  }

  //if(a.eqr(1)) return jmatc(1);  // Not true, e.g. 1^Infinity is not 1. Probably things are not like this for tetration either. Indeed it looks like e.g. 1^^0.5 is not 1.
  if(z.eqr(0)) return jmatc(1); //by definition
  if(z.eqr(1)) return a;
  if(z.eqr(2)) return a.pow(a);
  if(JmatC.isReal(a) && a.re >= 2 && z > 5) return jmatc(Infinity); // too big for double
  if(a.eqr(0) && JmatC.isPositiveInt(z)) {
    // Chosen definition that considers 0^0 = 1, so 0^^n is 1 if n is even, 0 if n is odd.
    return JmatC.isEven(z) ? jmatc(1) : jmatc(0);
  }

  // Power tower (infinitely iterated exponentiation)
  if(z.eqr(Infinity) /*&& JmatC.isPositive(a)*/) {
    if(a.eqr(1)) return jmatc(1); //0/0 ==> 1
    // converges if a >= 0.066 && a <= 1.44
    // when using "runloop" with high iterations here, it it indeed only converges there. The lambertw formula below has values everywhere though.
    var l = JmatC.log(a);
    return JmatC.lambertw(l.neg()).div(l.neg());
  }

  var runloop = function(a, b, num, l) {
    var result = b;
    var last;
    for(var i = 0; i < num; i++) {
      if(l) result = JmatC.logy(result, a);
      else result = a.pow(result);
      if(JmatC.isNaN(result)) return result;
      if(result.eq(last)) return result; // E.g. a=1.01 already converges to 1.0101015237405409 after just 7 iterations. 1.44 converges to 2.393811748202943 after 244 iterations, so numbers near that take a while to loop. 1.45 and higher converge to Infinity.
      last = result;
      if(i > 1000) return jmatc(NaN); //avoid infinite loop
    }
    return result;
  }

  // This supports power tower as well
  // It supports it better than the lambertw formula above because the domain here is real numbers only, the lambertw formula fails on, for example, -1 and returns nan instead of -1
  if(JmatC.isPositiveInt(z)) {
    return runloop(a, a, z.re - 1, false);
  }

  // Everything above is true tetration for those cases where possible. What follows below is intermediate tetration research, to return "something", because there is no more than that available.
  if (JmatC.isReal(z)) {
    // Linear approximation for the extension to real heights
    // a^^z = z+1 for z > -1 && z <= 0
    // a^^z = log_a(a^^(z+1)) for z <= -1  ==> a^^-1.5 = log_a(z+2), a^^-2.5 = log_a(log_a(z+3)), etc...
    // a^^z = a^(a^^(z-1)) for z > 0  ==> a^^0.5 = a^z, a^^1.5 = a^(a^(z-1)), a^^2.5 = a^(a^(a^(z-2))), etc...
    // examples: e^^(0.5*pi) ~= 5.868, 0.5^^(-4.3) ~= 4.033
    if(z.eqr(-1)) {
      return jmatc(0); //the formulas below would give -Infinity
    }
    if(z.re > -1 && z.re <= 0) {
      return z.inc();
    }
    if(z.re > 0) {
      var b = z.sub(JmatC.floor(z)); //always in range 0-1
      return runloop(a, b, Math.ceil(z.re), false);
    }
    if(z.re <= -1) {
      var b = z.sub(JmatC.floor(z)); //always in range 0-1
      return runloop(a, b, -Math.ceil(z.re), true);
    }
  }

  if (JmatC.near(a, JmatC.E, 1e-15)) {
    // This implementation, which only works for base e, is based on the paper: Tetration as special function, E. Kouznetsov

    // TODO: more coefficients
    var s = [0.30685281944005, 0.59176735125832, 0.39648321290170, 0.17078658150959,
             0.08516537613999, 0.03804195209047, 0.01734090876306, 0.00755271038865,
             0.00328476064839, 0.00139361740170, 0.00058758348148, 0.00024379186661,
             0.00024379186661, 0.00010043966462, 0.00001654344436, 0.00000663102846,
             0.00000264145664, 0.00000104446533, 0.00000041068839, 0.00000016048059,
             0.00000006239367, 0.00000002412797, 0.00000000928797, 0.00000000355850,
             0.00000000135774, 0.00000000051587];
    // TODO: does javascript always reexecute all these constructors? Make object contant outside of this function instead
    var t = [jmatc(0.37090658903229, 1.33682167078891), jmatc(0.01830048268799, 0.06961107694975),
             jmatc(-0.04222107960160, 0.02429633404907), jmatc(-0.01585164381085, -0.01478953595879),
             jmatc(0.00264738081895, -0.00657558130520), jmatc(0.00182759574799, -0.00025319516391),
             jmatc(0.00036562994770, 0.00028246515810), jmatc(0.00002689538943, 0.00014180498091),
             jmatc(-0.00003139436775, 0.00003583704949), jmatc(-0.00001376358453, -0.00000183512708),
             jmatc(-0.00000180290980, -0.00000314787679), jmatc(0.00000026398870, -0.00000092613311),
             jmatc(0.00000024961828, -0.00000013664223), jmatc(0.00000000637479, 0.00000002270476),
             jmatc(-0.00000000341142, 0.00000000512289), jmatc(-0.00000000162203, 0.00000000031619),
             jmatc(-0.00000000038743, -0.00000000027282), jmatc(-0.00000000001201, -0.00000000013440),
             jmatc(0.00000000002570, -0.00000000002543), jmatc(0.00000000000935, 0.00000000000045),
             jmatc(0.00000000000170, 0.00000000000186), jmatc(-0.00000000000005, 0.00000000000071),
             jmatc(-0.00000000000016, 0.00000000000012), jmatc(-0.00000000000005, -0.00000000000001),
             jmatc(-0.00000000000001, -0.00000000000001)];
    var an = [jmatc(0.3181315052047641353, 1.3372357014306894089), jmatc(1), jmatc(-0.1513148971556517359, -0.2967488367322413067),
              jmatc(-0.03697630940906762, 0.09873054431149697), jmatc(0.0258115979731401398, -0.017386962126530755), jmatc(-0.0079444196, 0.00057925018)];
    var fima = function(z) {
      var r = jmatc(1.0779614375280, -0.94654096394782);
      var beta = jmatc(0.12233176, -0.02366108);
      var l = an[0]; //fixed point of the logarithm
      var e = JmatC.exp(l.mul(z).add(r));
      var c = beta.mul(e).mul(JmatC.exp(z.mul(JmatC.newi(Math.PI * 2))));
      return JmatC.powerSeries(an, an.length, JmatC.ZERO, e).add(c);
    };

    var b;
    var z2 = jmatc(JmatR.fracn(z.re), z.im);

    if(z.im < -4.5) b = fima(z2.conj()).conj();
    else if(z.im < -1.5) b = JmatC.powerSeries(t, t.length, JmatC.newi(3), z2.conj()).conj();
    else if(z.im < 1.5) b = JmatC.log(z2.addr(2)).add(JmatC.powerSeries(s, s.length, JmatC.ZERO, z2));
    else if(z.im < 4.5) b = JmatC.powerSeries(t, t.length, JmatC.newi(3), z2);
    else b = fima(z2);

    if(z.re > 0) return runloop(a, b, Math.floor(z.re), false);
    else return runloop(a, b, -Math.ceil(z.re), true);
  }

  // TODO: for complex z and arbitrary base: implement something like "kneser" function

  return jmatc(NaN);

  // TODO: implement super logarithm (slog), and super root (sroot) as well. Though, so far the only formulas I have is slog of base e in the Kouznetsov paper, and the 2-super root, so no implementation using both parameters can be done so far.
};


// Faddeeva function, used as helper functions to calculate erf and related functions for certain locations in the complex plane
// Faddeeva(z) = exp(-z^2)*erfc(-iz).
// Also known as Faddeyeva, or as w(x), but that may be confusing with LambertW...
JmatC.faddeeva = function(z) {
  // METOD A: series 7.1.8 from Handbook of Mathematical Functions
  // smaller area of convergence than METHOD A, so not used
  /*var result = jmatc(0);
  var zi = JmatC.I.mul(z);
  var zz = JmatC.ONE;
  for(var n = 0; n < 30; n++) {
    result = result.add(zz.divr(JmatR.gamma(n/2 + 1)));
    zz = zz.mul(zi);
  }
  return result;*/

  var invsqrtpi2   = 2 / JmatR.SQRTPI;
  var eye = (z.re * z.re) + (z.im * z.im * 2);

  // METHOD B: series
  // A small eye-shaped region in which the series works
  if(eye < 3.5) {
    // Based on sum 7.1.5 from Handbook of Mathematical Functions
    // erf(z) = 2/sqrt(pi) * SUM_n=0..oo (-1)^n * z^(2n+1) / (n! * (2n+1))
    // and then w(z) = e^(-z^2) * (1 - erf(iz))
    var sum = JmatC.ZERO;
    var sign = 1.0;
    var nn = 1;
    var iz = JmatC.I.mul(z).neg();
    var izz = iz;
    for(var n = 0; n < 20; n++) {
      if(n > 0) {
        nn = nn * n; // n!
        sign = -sign; // (-1)^n
        izz = izz.mul(iz).mul(iz); // iz^(2n+1)
      }
      sum = sum.add(izz.mulr(sign / (nn * (2*n + 1))));
    }
    var e = JmatC.exp(z.mul(z).neg());
    return e.sub(e.mul(sum).mulr(invsqrtpi2));
  }

  // METHOD C: Laplace Continued Fraction
  var za = jmatc(Math.abs(z.re), Math.abs(z.im)); // Operate on positive re, positive im quadrant
  // requires quite a lot of iterations unfortunately
  var num = eye < 40 ? 40 : eye < 80 ? 20 : 10;
  var result = jmatc(0);
  for(var n = 0; n < num; n++) {
    var r = jmatc(za.im + result.re, za.re - result.im);
    result = r.mulr(0.5 / r.abssqr());
  }
  result = result.mulr(invsqrtpi2);
  // Fix for pure imaginary values with large negative imaginary part
  if(za.im == 0.0) result.re = Math.exp(-za.re * za.re);
  // Put the solution back in the original quadrant, using the transformations w(-z) = 2 * exp(-z*z) - w(z) and w(conj(z)) = cons(w(-z))
  if(z.im < 0.0) {
    var e = JmatC.exp(za.mul(za).neg()).mulr(2);
    result = e.sub(result);
    if(z.re > 0.0) result.im = -result.im;
  } else if(z.re < 0.0) {
    result.im = -result.im;
  }
  return result;
};


// erfcx(z) = exp(z^2) * erfc(z): the scaled complementary error function
JmatC.erfcx = function(z) {
  return JmatC.faddeeva(jmatc(-z.im, z.re)); //erfcx(z) = faddeeva(iz)
};

JmatC.erf = function(z) {
  if(z.im == 0) {
    return jmatc(JmatR.erf(z.re));
  } else if(z.re == 0) {
    return JmatC.I.mulr(JmatR.erfi(z.im));
  } else {
    var a = JmatC.exp(z.mul(z).neg()); // If abs of z is very large, and |im| > |re|, then this becomes some NaN or Infinity. That is ok, erf is also some unrepresentably huge value there.
    if (z.re >= 0) return JmatC.ONE.sub(a.mul(JmatC.faddeeva(z.mul(JmatC.I))));
    else return a.mul(JmatC.faddeeva(z.mul(JmatC.I.neg()))).sub(JmatC.ONE);

    // With integration, don't use.
    /*var ps2 = 2.0 / JmatR.SQRTPI;
    var result;
    result = JmatC.integrate(jmatc(0), z, function(z){ return JmatC.exp(z.mul(z).neg()); }, 100);
    result = result.mulr(ps2);
    return result;*/
  }
};

// erfc(x) = 1 - erf(x). This function gives numerically a better result if erf(x) is near 1.
JmatC.erfc = function(z) {
  if(z.im == 0) {
    return jmatc(JmatR.erfc(z.re));
  } else {
    var a = JmatC.exp(z.mul(z).neg());
    if (z.re >= 0) return a.mul(JmatC.faddeeva(z.mul(JmatC.I)));
    else return JmatC.TWO.sub(a.mul(JmatC.faddeeva(z.mul(JmatC.I.neg()))));
  }
};


// TODO: rewrite some of the rational approximations to not use so many multiplications
//a + bx + cxx + dxxx + exxxx = a + x * (b + x * (c + x * (d + x * e)))   etc...
//and that equals: x.mulr(e).addr(d).mul(x).addr(c).mul(x).addr(b).mul(x).addr(a) ...


JmatC.erf_inv = function(z) {
  if (z.im != 0 && Math.abs(z.re) > 1) {
    //this branch is taken for large complex numbers because the implementation below doesn't work well on those. This one isn't much better btw, but slightly less bad for those cases.
    //TODO: is incorrect on many complex numbers!! e.g. erf_inv(erf(5 + 5i)) gives way wrong result
    //var zzpi = z.mul(z).mulr(Math.PI);
    //return zzpi.mulr(34807/182476800.0).addr(4369/5806080.0).mul(zzpi).addr(127/40320.0).mul(zzpi).addr(7/480.0).mul(zzpi).addr(1/12.0).mul(zzpi).addr(1).mul(z).mul(JmatC.SQRTPI).mulr(0.5);

    // With newton method

    // derivative of erf is: 2/sqrt(pi) * e^(-x^2)
    var derf = function(x) {
      return JmatC.TWO.divr(JmatR.SQRTPI).mul(JmatC.exp(x.mul(x).neg()));
    }

    var neg = z.re < 0;
    if(neg) z = z.neg();

    // For abs(z) > 1 and z.re > 0, the following starting value works well: sqrt(-log(x * sqrt(pi) * (1 - x)))
    // NOTE: for abs(z) < 1, instead z*sqrtpi/2 would work, but other code already handles such case
    var start = JmatC.sqrt(JmatC.log(z.mulr(JmatR.SQRTPI).mul(JmatC.ONE.sub(z))).neg());
    // NOTE: erf_inv has multiple solutions, with this starting value only one particular one is returned.
    // e.g. with the chosen starting value, erf_inv(2+2i) gives 0.386507600275 + 1.320769860731i. But with starting value 0, it gives the also correct 2.947736167125 + 3.401249486995i.
    var result = JmatC.finvert_newton(z, JmatC.erf, derf, start);
    if(neg) result = result.neg();
    return result;
  } else {
    //if (a > 1) return jmatc(NaN); //only relevant for real numbers
    if (z.im == 0) {
      if (z.re == 0) return jmatc(0);
      if (z.re == 1) return jmatc(Infinity);
      if (z.re == -1) return jmatc(-Infinity);
    }

    var erf_inv_a_ = [0.886226899, -1.645349621, 0.914624893, -0.140543331];
    var erf_inv_b_ = [1, -2.118377725, 1.442710462, -0.329097515, 0.012229801];
    var erf_inv_c_ = [-1.970840454, -1.62490649, 3.429567803, 1.641345311];
    var erf_inv_d_ = [1, 3.543889200, 1.637067800];

    var a = JmatC.abs(z).re;
    if (a <= 0.7) {
      var z2 = z.mul(z);
      var r = z.mul(z2.mulr(erf_inv_a_[3]).addr(erf_inv_a_[2]).mul(z2).addr(erf_inv_a_[1]).mul(z2).addr(erf_inv_a_[0]));
      r = r.div(z2.mulr(erf_inv_b_[4]).addr(erf_inv_b_[3]).mul(z2).addr(erf_inv_b_[2]).mul(z2).addr(erf_inv_b_[1]).mul(z2).addr(erf_inv_b_[0]));
    }
    else {
      var y = JmatC.sqrt(JmatC.log(JmatC.ONE.sub(z).divr(2)).neg());
      var r = y.mulr(erf_inv_c_[3]).addr(erf_inv_c_[2]).mul(y).addr(erf_inv_c_[1]).mul(y).addr(erf_inv_c_[0]);
      r = r.div(y.mulr(erf_inv_d_[2]).addr(erf_inv_d_[1]).mul(y).addr(erf_inv_d_[0]));
    }

    return r;
  }
};

// inverse complementary error function.
JmatC.erfc_inv = function(z) {
  return JmatC.erf_inv(JmatC.ONE.sub(z));
};

//erfi(z) = -i erf(iz)
JmatC.erfi = function(z) {
  if(JmatC.isReal(z)) return jmatc(JmatR.erfi(z.re));
  return JmatC.erf(z.mul(JmatC.I)).mul(JmatC.I).neg();
};

// D+(x) aka F(x)
JmatC.dawson = function(z) {
  if(JmatC.isReal(z)) {
    return jmatc(JmatR.dawson(z.re));
  } else {
    var w = JmatC.faddeeva(z);
    var a = JmatC.exp(z.mul(z).neg());
    return a.sub(w).mul(JmatC.I.mulr(JmatR.SQRTPI / 2));
  }
};

//Minkowski's question mark function, from Wikipedia
JmatC.minkowski = function(z) {
  return jmatc(JmatR.minkowski(z.re), JmatR.minkowski(z.im));
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Numerical function algorithms: root finding, integration, differentiation, sum, product
////////////////////////////////////////////////////////////////////////////////

//find one zero between valuex and valuey in function f (real function, complex isn't supported)
//this is totally not an efficient algorithm like Newton's method let alone Brent's method. Also valuex and valuey are allowed to be anything and to have function value with the same sign.
//if no zero is found, it returns NaN
JmatC.rootfind_bisection = function(valuex, valuey, f, maxit, prec) {
  var steps = (maxit == undefined ? 256 : maxit);
  var precision = (prec == undefined ? 0.000000001 : prec); //for 'near'
  var x = valuex.re;
  var y = valuey.re;
  if(y < x) {
    var temp = x;
    x = y;
    y = temp;
  }

  if(y - x == Infinity) return jmatc(NaN);

  //find a positive and a negative value
  var p = NaN;
  var n = NaN;
  var lastp = 0;
  var lastn = 0;
  for(var i = 0; i < steps; i++) {
    var z = x + i * (y - x) / steps; //exclude the leftmost point by adding precision. Otherwise annoyingness can occur.
    var fz = f(jmatc(z)).re;
    if(JmatR.near(fz, 0, precision)) return jmatc(z);

    if(fz > lastp) {
      p = z;
      lastp = fz;
    }
    if(fz < lastn) {
      n = z;
      lastn = fz;
    }

    //line below commented out, instead keep searching for better values thanks to lastp and lastn
    //if(n == n && p == p) break; // not NaN
  }

  if(n != n || p != p) {
    return jmatc(NaN);
  }

  for(;;) {
    var z = (n + p) / 2;
    var fz = f(jmatc(z)).re;
    if(JmatR.near(fz, 0, precision)) {
      return jmatc(z);
    }

    if(fz > 0) p = z;
    if(fz < 0) n = z;

    if(JmatR.near(n, p, precision)) {
      return jmatc(NaN); //not found
    }
  }
};

//root finding, aka find zeroes (findZero)
//most parameters are optional. Based on which are given, a certain algorithm is chosen (newton, bisection, ...)
//f: the function to find root of
//o: object with the following optional values:
// o.z0: starting value, or, if z1 is given, lower value of range (and starting value is assumed in the center of both). Default: 0
// o.z1: end value, if range is given. Default: undefined
// o.real: Whether to only find real zeroes. If true, f is assumed to be a real function and complex roots are ignored. Default: false
// o.df: derivative of f. If given, something like newton's method can be used. Default: undefined
// o.maxit: max iterations for e.g. newton. Default: 30
// o.prec: precision for e.g. newton. Default: 0.000000001
JmatC.rootfind = function(f, o) {
  if(!o) o = {};
  var z0 = o.z0 || JmatC.ZERO;
  var maxit = o.maxit || 30;
  var prec = (o.prec == undefined ?  0.000000001 : o.prec);

  // TODO: work in progress. Find better start values. Try other root finding algorithms. Etc...
  if(o.real && o.z1 != undefined) return JmatC.rootfind_bisection(z0, z1, f, maxit, prec);
  if(o.df) return JmatC.rootfind_newton(f, o.df, z0, maxit);
  return JmatC.rootfind_secant(f, z0, maxit);
};

JmatC.newtonStartValues_ = [
    jmatc(0),
    jmatc(0.1), jmatc(-0.1), JmatC.newi(0.1), JmatC.newi(-0.1), 
    jmatc(0.1, 0.1), jmatc(0.1, -0.1), jmatc(-0.1, 0.1), jmatc(-0.1, -0.1), 
    jmatc(1), jmatc(-1), JmatC.newi(1), JmatC.newi(-1)
  ];

// This is really not that good. TODO: better root finding
JmatC.newtonStartValue_ = function(f) {
  var s = JmatC.newtonStartValues_;
  var bestdist = Infinity;
  var best = jmatc(NaN);

  for(var i = 0; i < s.length; i++) {
    var z = s[i];
    var fz = f(z);
    var m = JmatC.manhattan(fz, JmatC.ZERO);
    if(m < bestdist) {
      bestdist = m;
      best = z;
    }
  }

  return best;
};

JmatC.newtonStartValuesAround_ = [ jmatc(1), jmatc(-1), JmatC.newi(1), JmatC.newi(-1), ];

//Find a new value to start at after having accidently encountered a bad point like NaN or Inf
//dist: e.g. 0.1 or 0.01
// This is really not that good. TODO: better root finding
JmatC.newtonStartValueAround_ = function(f, z0, dist) {
  var s = JmatC.newtonStartValuesAround_;
  var bestdist = Infinity;
  var best = jmatc(NaN);

  for(var i = 0; i < s.length; i++) {
    var z = z0.add(s[i].mulr(dist));
    var fz = f(z);
    var m = JmatC.manhattan(fz, JmatC.ZERO);
    if(m < bestdist) {
      bestdist = m;
      best = z;
    }
  }

  return best;
};

//Newton-Raphson. Finds a complex root (zero) given function f, its derivative df, and an initial value z0
JmatC.rootfind_newton = function(f, df, z0, maxiter) {
  if (!z0) z0 = JmatC.ZERO;//JmatC.newtonStartValue_(f);//JmatC.ZERO;
  if (!maxiter) maxiter = 30;
  var z = z0;
  var prevz = z;
  var bestdist = Infinity;
  var best = jmatc(NaN);
  for (var i = 0; i < maxiter; i++) {
    var fz = f(z);
    var m = JmatC.manhattan(fz, JmatC.ZERO);
    if(JmatR.near(m, 0, 1e-15)) return z; // Near enough, stop iterations
    if(m < bestdist) {
      bestdist = m;
      best = z;
    }
    var zn = z.sub(fz.div(df(z)));
    if(JmatC.isInfOrNaN(zn)) {
      var m = JmatC.manhattan(z, prevz);
      zn = JmatC.newtonStartValueAround_(f, z, m ? m : 0.1);
    }
    z = zn;
    prevz = z;
  }

  return best;
};

//finds a complex root (zero) given function f, and an initial value z0 (no need to give the derivative)
JmatC.rootfind_secant = function(f, z0, maxiter) {
  return JmatC.rootfind_newton(f, function(x) {
    return JmatC.differentiate_stencil5(x, f);
  }, z0, maxiter);
};

//find result of inverse function using the newton method
JmatC.finvert_newton = function(z, f, df, z0, maxiter) {
  return JmatC.rootfind_newton(function(x) { return f(x).sub(z); }, df, z0, maxiter);
};

//find result of inverse function using the newton method (no need to give the derivative)
JmatC.finvert_secant = function(z, f, z0,  maxiter) {
  return JmatC.rootfind_secant(function(x) { return f(x).sub(z); }, z0, maxiter);
};


//numerical integration, aka quadrature
//integrate function f using simpsons rule, from x to y, with amount of steps given by steps parameter (higher = more precision, number of evaluations is throughly steps * 2)
//NOTE: this is the real version, real JS numbers only. Complex version is JmatC.integrate_simpson
JmatR.integrate_simpson = function(x, y, steps, f, stopLoop) {
  var step = (y - x) / steps;
  var result = 0;
  var fa = 0;

  for(var i = 0; i <= steps; i++) {
    var a = x + (i - 1) * step;
    var b = x + i * step;

    var fab = f((a + b) / 2);
    var fb = f(b);

    //Simpson's rule
    if(i > 0) result += ((b - a) / 6) * (fa + 4 * fab + fb);

    fa = fb;

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return NaN;
  }
  return result;
};

//numerical integration, aka quadrature
//NOTE: this is the real version, real JS numbers only. Complex version is JmatC.integrate
JmatR.integrate = function(x, y, f, steps) {
  if(!steps) steps = 30;
  return JmatR.integrate_simpson(x, y, steps, f);
};

//numerical integration, aka quadrature
//integrate function f using simpsons rule, from x to y, with amount of steps given by steps parameter (higher = more precision, number of evaluations is throughly steps * 2)
JmatC.integrate_simpson = function(x, y, steps, f, stopLoop) {
  var step = y.sub(x).divr(steps);
  var result = jmatc(0);
  var fa = null;

  var a, b;

  for(var i = 0; i <= steps; i++) {
    if(y.eqr(Infinity)) {
      a = (b ? b : x);
      b = (b ? b.mulr(2) : x.addr(0.1));
    } else {
      a = x.add(step.mulr(i - 1));
      b = x.add(step.mulr(i));
    }

    var fab = f(a.add(b).divr(2));
    var fb = f(b);

    if(i > 0) {
      //Simpson's rule
      var s = b.sub(a).divr(6).mul(fa.add(fab.mulr(4)).add(fb));
      result = result.add(s);
    }

    fa = fb;

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return jmatc(NaN);
  }
  return result;
};

//numerical integration, aka quadrature
JmatC.integrate = function(x, y, f, steps) {
  if(!steps) steps = 30;
  return JmatC.integrate_simpson(x, y, steps, f);
};

// differentiation with just two points (finite difference, or secant)
JmatC.differentiate_secant = function(x, f) {
  //var h = jmatc(0.0001);
  var h = Math.max(0.01, Math.abs(x.re)) / 1000;

  var f1 = f(x.addr(h / 2));
  var f2 = f(x.subr(h / 2));
  return f1.sub(f2).divr(h);
};

// differentiation with 5-point stencil
// Ignores imaginary direction (e.g. derivative of the function Im(x) is always 0 with this. Im(x) is not holomorphic though), but that seems to be alright, Wolfram Alpha does that too...
JmatC.differentiate_stencil5 = function(x, f) {
  //var h = jmatc(0.0001);
  var h = Math.max(0.01, Math.abs(x.re)) / 1000;

  var f1 = f(x.addr(h * 2)).neg();
  var f2 = f(x.addr(h)).mulr(8);
  var f3 = f(x.subr(h)).mulr(-8);
  var f4 = f(x.subr(h * 2));
  return f1.add(f2).add(f3).add(f4).divr(h * 12);
};

// second derivative with 5-point stencil
JmatC.differentiate2nd_stencil5 = function(x, f) {
  //var h = jmatc(0.0001);
  var h = jmatc(Math.max(0.01, Math.abs(x.re)) / 1000);

  var f1 = f(x.add(h.mulr(2))).neg();
  var f2 = f(x.add(h)).mulr(16);
  var f3 = f(x).mulr(-30);
  var f4 = f(x.sub(h)).mulr(16);
  var f5 = f(x.sub(h.mulr(2))).neg();

  return f1.add(f2).add(f3).add(f4).add(f5).div(h.mul(h).mulr(12)); // (f1+f2+f3+f4+f5) / (h*h*12)
};

// differentiation (derivative)
JmatC.differentiate = function(x, f) {
  return JmatC.differentiate_stencil5(x, f);
};

//do summation of discrete values of f using start value x, end value y (inclusive) and given step.
JmatC.doSummation = function(valuex, valuey, step, f, stopLoop) {
  if(step == 0) return jmatc(NaN);
  if(step < 0) return jmatc(NaN);
  if(!JmatC.isReal(valuex) || !JmatC.isReal(valuey)) return jmatc(NaN);
  var x = valuex.re;
  var y = valuey.re;
  var result = jmatc(0);
  // the step / 4 thing is to avoid numerical problems that may let it miss the last value
  for(var z = x; z <= y + step / 4; z += step) {
    var fz = f(jmatc(z));
    result = result.add(fz);

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return jmatc(NaN);
  }
  return result;
};

//do product of discrete values of f using start value x, end value y (inclusive) and given step.
JmatC.doProduct = function(valuex, valuey, step, f, stopLoop) {
  if(step == 0) return jmatc(NaN);
  if(step < 0) return jmatc(NaN);
  if(!JmatC.isReal(valuex) || !JmatC.isReal(valuey)) return jmatc(NaN);
  var x = valuex.re;
  var y = valuey.re;
  var result = JmatC.ONE;
  // the step / 4 thing is to avoid numerical problems that may let it miss the last value
  for(var z = x; z <= y + step / 4; z += step) {
    var fz = f(jmatc(z));
    result = result.mul(fz);

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return jmatc(NaN);
  }
  return result;
};

// n is max coeff.length. Sum is from 0..n-1
JmatC.powerSeries = function(coeff, n, z0, z) {
  var realcoeff = coeff[0].re == undefined;
  var result = JmatC.ZERO;
  var zz = JmatC.ONE;
  for(var i = 0; i < n; i++) {
    if(realcoeff) result = result.add(zz.mulr(coeff[i]));
    else result = result.add(zz.mul(coeff[i]));
    zz = zz.mul(z.sub(z0));
  }
  return result;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Statistical Distributions
////////////////////////////////////////////////////////////////////////////////

/*
For each distribution, the following function prefixes are available:

pdf = probability density function for continuous distributions, pmf = probability mass function for discrete distributions
cdf = cumulative distribution function, integral of the pdf
qf = quantile function, the inverse function of cdf
*/


// TODO: add discrete distributions like binomial, bernoulli, poisson, ...
// TODO: add characteristic functions cf_
// TODO: Pareto distribution

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_uniform = function(x, a, b) {
  if(x.re >= a.re && x.re <= b.re) return b.sub(a).inv();
  return jmatc(0);
};

JmatC.cdf_uniform = function(x, a, b) {
  if(x.re < a.re) return jmatc(0);
  if(x.re < b.re) return x.sub(a).div(b.sub(a));
  return jmatc(1);
};

JmatC.qf_uniform = function(x, a, b) {
  var r = a.add(x.mul(b.sub(a)));
  if(r.re >= a.re && r.re <= b.re) return r;
  return jmatc(NaN);
};

////////////////////////////////////////////////////////////////////////////////

//aka small phi
JmatC.pdf_standardnormal = function(x) {
  return JmatC.exp(x.mul(x).mulr(-0.5)).mul(JmatC.INVSQRT2PI);
};

//aka capital PHI
JmatC.cdf_standardnormal = function(x) {
  return JmatC.erf(x.div(JmatC.SQRT2)).addr(1).mulr(0.5);
};

//aka the probit function
JmatC.qf_standardnormal = function(x) {
  return JmatC.erf_inv(x.mulr(2).subr(1)).mul(JmatC.SQRT2);
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_normal = function(x, mu, sigma) {
  var a = JmatC.INVSQRT2PI.div(sigma);
  var b = x.sub(mu).mul(x.sub(mu)).div(sigma.mul(sigma).mulr(2)); // (x-mu)^2 / 2*sigma^2
  return a.mul(JmatC.exp(b.neg()));
};

JmatC.cdf_normal = function(x, mu, sigma) {
  var a = x.sub(mu).div(JmatC.abs(sigma).mul(JmatC.SQRT2)); // (x-mu) / sqrt(2*sigma^2)
  return JmatC.erf(a).addr(1).mulr(0.5);
};

JmatC.qf_normal = function(x, mu, sigma) {
  return JmatC.exp(mu.add(sigma.mul(JmatC.qf_standardnormal(x))));
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_lognormal = function(x, mu, sigma) {
  // (1 / (x * sqrt(2*pi) * sigma)) * exp(- (ln(x) - mu)^2 / (2*sigma^2))
  var a = JmatC.INVSQRT2PI.div(sigma).div(x);
  var b = JmatC.log(x).sub(mu);
  return a.mul(JmatC.exp(b.mul(b).div(sigma.mul(sigma).mulr(2)).neg()));
};

JmatC.cdf_lognormal = function(x, mu, sigma) {
  var a = JmatC.log(x).sub(mu).div(JmatC.SQRT2.mul(sigma));
  return JmatC.erf(a).addr(1).mulr(0.5);
};

JmatC.qf_lognormal = function(x, mu, sigma) {
  var a = JmatC.log(x).sub(mu).div(JmatC.SQRT2.mul(sigma));
  return JmatC.erf(a).addr(1).mulr(0.5);
};

////////////////////////////////////////////////////////////////////////////////

//gamma = scale parameter (HWHM)
JmatC.pdf_cauchy = function(x, x0, gamma) {
  var x2 = x.sub(x0);
  var d = x2.mul(x2).add(gamma.mul(gamma)).mulr(Math.PI);
  return gamma.div(d);
};

JmatC.cdf_cauchy = function(x, x0, gamma) {
  return JmatC.atan(x.sub(x0).div(gamma)).divr(Math.PI).addr(0.5);
};

JmatC.qf_cauchy = function(x, x0, gamma) {
  return x0.add(gamma.mul(JmatC.tan(x.subr(0.5).divr(Math.PI))));
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_studentt_cache_ = []; // cache used because nu is often constant between calls
JmatC.pdf_studentt_cachefun_ = function(nu) { var nu2 = nu.inc().divr(2); return JmatC.gammaDiv_(nu2, nu.divr(2)); };

// nu = degrees of freedom
JmatC.pdf_studentt = function(x, nu) {
  if(nu.eqr(1)) {
    return x.mul(x).addr(1).mulr(Math.PI).inv(); // cauchy distribution
  }
  if(nu.eqr(2)) {
    return x.mul(x).addr(2).powr(3/2).inv();
  }
  if(nu.eqr(3)) {
    var sqrt3x6 = 10.392304845413264; //6 * sqrt(3)
    var x2 = x.mul(x).addr(3);
    return jmatc(sqrt3x6).div(x2.mul(x2).mulr(Math.PI));
  }
  if(nu.eqr(Infinity)) {
    return JmatC.pdf_standardnormal(x);
  }

  var nu2 = nu.inc().divr(2);
  var g = JmatC.calcCache_(nu, JmatC.pdf_studentt_cachefun_, JmatC.pdf_studentt_cache_);
  var s = JmatC.sqrt(JmatC.PI.mul(nu)).inv();
  var gs = g.mul(s);
  if(JmatC.isNaN(gs) && nu.re > 100) gs = jmatc(JmatC.INVSQRT2PI); //this is what it is for nu = +Infinity
  
  var a = x.mul(x).div(nu).inc();
  return gs.mul(a.pow(nu2.neg()));
};

JmatC.cdf_studentt = function(x, nu) {
  if(nu.eqr(1)) {
    return JmatC.atan(x).divr(Math.PI).addr(0.5); // cauchy distribution
  }
  if(nu.eqr(2)) {
    return x.div(JmatC.sqrt(x.mul(x).addr(2)).mulr(2)).addr(0.5);
  }
  if(nu.eqr(Infinity)) {
    return JmatC.cdf_standardnormal(x);
  }

  // METHOD A: in terms of incomplete beta (probably slightly faster than with hypergeometric)
  if(x.eqr(0)) return jmatc(0.5);
  var nu2 = nu.inc().divr(2);
  var g = JmatC.calcCache_(nu, JmatC.pdf_studentt_cachefun_, JmatC.pdf_studentt_cache_);
  var b = JmatC.incbeta(x.mul(x).div(nu).neg(), jmatc(0.5), JmatC.ONE.sub(nu).mulr(0.5));
  var n = JmatC.I.mul(x).mul(b);
  var d = x.abs().mulr(2).mul(JmatC.SQRTPI);
  return jmatc(0.5).sub(g.mul(n).div(d));

  // METHOD B: in terms of hypergeometric
  /*var nu2 = nu.inc().divr(2);
  var g = JmatC.calcCache_(nu, JmatC.pdf_studentt_cachefun_, JmatC.pdf_studentt_cache_);
  var s = JmatC.sqrt(JmatC.PI.mul(nu)).inv();
  var gs = g.mul(s);
  if(JmatC.isNaN(gs) && nu.re > 100) gs = jmatc(JmatC.INVSQRT2PI); //this is what it is for nu = +Infinity
  var f = JmatC.hypergeometric(jmatc(0.5), nu.inc().mulr(0.5), jmatc(1.5), x.mul(x).div(nu).neg());
  return jmatc(0.5).add(x.mul(gs).mul(f));*/
};

//aka "tinv", t inverse cumulative distribution function
JmatC.qf_studentt= function(x, nu) {
  // test in console:
  // function testInvStudentt(x, nu) { return JmatC.qf_studentt(JmatC.cdf_studentt(jmatc(x), jmatc(nu)), jmatc(nu)).re; }
  // testInvStudentt(0.5,1.5)
  if(nu.eqr(1)) {
    return JmatC.tan(JmatC.PI.mul(x.subr(0.5))); // cauchy distribution
  }
  if(nu.eqr(2)) {
    var a = x.mulr(4).mul(JmatC.ONE.sub(x));
    return x.subr(0.5).mulr(2).mul(JmatC.sqrt(jmatc(2).div(a)));
  }
  if(nu.eqr(4) && JmatC.isReal(x)) {
    var a = x.mulr(4).mul(JmatC.ONE.sub(x));
    var as = JmatC.sqrt(a);
    var q = JmatC.cos(jmatc(1.0 / 3).mul(JmatC.acos(as))).div(as);
    return JmatC.sign(x.subr(0.5)).mulr(2).mul(JmatC.sqrt(q.dec()));
  }
  if(nu.eqr(Infinity)) {
    return JmatC.qf_standardnormal(x);
  }

  // Formula for student t inverse CDF in terms of inverse regularized beta function.
  // Works for real x in range 0-1. Extension to complex plane is fully of numerical problems here, so not supported.
  if(JmatR.near(x.im, 0, 1e-15)) x = jmatc(x.re);
  if(JmatR.near(nu.im, 0, 1e-15)) nu = jmatc(nu.re);
  if(JmatC.isPositive(nu) && JmatC.isPositive(x) && x.re < 1) {
    if(x.re < 0.5) {
      var i = JmatC.beta_i_inv(x.mulr(2), nu.divr(2), jmatc(0.5));
      return JmatC.sqrt(nu.mul(i.inv().subr(1))).neg();
    } else {
      var i = JmatC.beta_i_inv(JmatC.ONE.sub(x).mulr(2), nu.divr(2), jmatc(0.5));
      return JmatC.sqrt(nu.mul(i.inv().subr(1)));
    }
  }

  return jmatc(NaN); //TODO: support approximate inverse CDF for student t for arbitrary nu
};

////////////////////////////////////////////////////////////////////////////////

//k = degrees of freedom
JmatC.pdf_chi_square = function(x, k) {
  var k2 = k.divr(2);
  var a = k2.eqr(1) ? jmatc(2) : jmatc(2).pow(k2);
  var g = JmatC.isNegativeInt(k2) ? jmatc(Infinity) : JmatC.gamma(k2);
  var b = x.pow(k2.dec()).mul(JmatC.exp(x.divr(-2)));
  return a.mul(g).inv().mul(b);
};

JmatC.cdf_chi_square = function(x, k) {
  return JmatC.gamma_p(k.divr(2), x.divr(2));
};

// aka "cinv". x is in range 0-1
JmatC.qf_chi_square = function(x, k) {
  // Does NOT work with '0' as starting value - and not with '0.5' either. '1' seems to be the best for the range of operation (x = 0..1)
  // TODO: not very precise for higher k such as k=5. Fix this.
  // E.g. should be: k=1,x=0.4: 0.274997, k=0.5&x=0.5: 0.087347000000, k=5&x=0.8:7.28928

  return JmatC.gamma_p_inv(k.divr(2), x).mulr(2);
};

////////////////////////////////////////////////////////////////////////////////

// mu = location, s = scale
JmatC.pdf_logistic = function(x, mu, s) {
  var e = JmatC.exp(x.sub(mu).div(s).neg());
  var ee = e.inc();
  return e.div(s.mul(ee).mul(ee));
};

JmatC.cdf_logistic = function(x, mu, s) {
  return JmatC.tanh(x.sub(mu).div(s).divr(2)).mulr(1/2).addr(1/2);
};

// Note: the *derivative* of this is the logit function s/(x*(1-x))
JmatC.qf_logistic = function(x, mu, s) {
  var xx = x.div(JmatC.ONE.sub(x));
  return mu.add(s.mul(JmatC.log(xx)));
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_gamma_cache_ = [];

// k = shape, theta = scale
JmatC.pdf_gamma = function(x, k, theta) {
  var xk = x.pow(k.dec());
  var e = JmatC.exp(x.div(theta).neg());
  var t = theta.pow(k);
  var g = JmatC.calcCache_(k, JmatC.gamma, JmatC.pdf_gamma_cache_); // gamma(k)
  return xk.mul(e).div(t).div(g);
};

JmatC.cdf_gamma = function(x, k, theta) {
  return JmatC.gamma_p(k, x.div(theta));
};

JmatC.qf_gamma = function(x, k, theta) {
  return JmatC.gamma_p_inv(k, x).mul(theta);
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_beta = function(x, alpha, beta) {
  var xa = x.pow(alpha.dec());
  var xb = x.rsub(1).pow(beta.dec());
  var b = JmatC.beta(alpha, beta);
  return xa.mul(xb).div(b);
};

JmatC.cdf_beta = function(x, alpha, beta) {
  return JmatC.beta_i(x, alpha, beta);
};

JmatC.qf_beta = function(x, alpha, beta) {
  return JmatC.beta_i_inv(x, alpha, beta);
};

////////////////////////////////////////////////////////////////////////////////

// Fisher–Snedecor F distribution 
// d1 and d2 are the degrees of freedom parameters
// TODO: gives NaN for x = 0, it probably should give 0 instead?
JmatC.pdf_fisher = function(x, d1, d2) {
  var a = d1.mul(x).pow(d1).mul(d2.pow(d2));
  var b = d1.mul(x).add(d2).pow(d1.add(d2));
  var c = x.mul(JmatC.beta(d1.divr(2), d2.divr(2)));
  return JmatC.sqrt(a.div(b)).div(c);
};

JmatC.cdf_fisher = function(x, d1, d2) {
  var a = d1.mul(x).div(d1.mul(x).add(d2));
  return JmatC.beta_i(a, d1.divr(2), d2.divr(2));
};

// aka "finv"
JmatC.qf_fisher = function(x, d1, d2) {
  var a = JmatC.beta_i_inv(x, d1.divr(2), d2.divr(2));
  var b = JmatC.ONE.sub(a);
  return d2.mul(a).div(d1.mul(b));
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_weibull = function(x, lambda, k) {
  var a = k.div(lambda);
  var b = (x.div(lambda)).pow(k.dec());
  var c = x.div(lambda).pow(k).neg();
  return a.mul(b).mul(JmatC.exp(c));
};

JmatC.cdf_weibull = function(x, lambda, k) {
  var c = x.div(lambda).pow(k).neg();
  return JmatC.ONE.sub(JmatC.exp(c));
};

JmatC.qf_weibull = function(x, lambda, k) {
  var a = JmatC.ONE.sub(x);
  var b = JmatC.log(a).neg().pow(k.inv());
  return lambda.mul(b);
};

////////////////////////////////////////////////////////////////////////////////

JmatC.pdf_exponential = function(x, lambda) {
  if(x.re < 0) return jmatc(0);
  return lambda.mul(JmatC.exp(x.mul(lambda).neg()));
};

JmatC.cdf_exponential = function(x, lambda) {
  if(x.re < 0) return jmatc(0);
  return JmatC.ONE.sub(JmatC.exp(x.mul(lambda).neg()));
};

JmatC.qf_exponential = function(x, lambda) {
  return JmatC.log(x.rsub(1)).neg().div(lambda);
};

////////////////////////////////////////////////////////////////////////////////

// aka "double exponential"
JmatC.pdf_laplace = function(x, mu, b) {
  // 1/2b * exp(-|x-mu\/b)
  var e = JmatC.exp(JmatC.abs(x.sub(mu)).neg().div(b));
  return b.mulr(2).inv().mul(e);
};

JmatC.cdf_laplace = function(x, mu, b) {
  var e = JmatC.exp(JmatC.abs(x.sub(mu)).neg().div(b));
  var s = JmatC.sign(x.sub(mu));
  return e.rsub(1).mul(s).mulr(0.5).addr(0.5);
};

JmatC.qf_laplace = function(x, mu, b) {
  var l = JmatC.log(JmatC.abs(x.subr(0.5)).mulr(2).rsub(1));
  var s = JmatC.sign(x.subr(0.5));
  return mu.sub(b.mul(s).mul(l));
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// JmatM: Matrix math
////////////////////////////////////////////////////////////////////////////////

/*
Constructor
height first because that's the math convention: a 2x3 matrix is 2 rows high, 3 columns wide, and made as new JmatM(2, 3)

This is the actual object used as matrix. In addition, most of the matrix
functions are implemented as static functions in here.
*/
function JmatM(height, width) {
  this.h = height; //number of rows
  this.w = width; //number of columns
  this.e = []; //array of arrays. first index is row (y), second index is column (x)
  for(var y = 0; y < height; y++) {
    this.e[y] = [];
    for(var x = 0; x < width; x++) {
      this.e[y][x] = jmatc(0);
    }
  }
}

// Makes a matrix from many types of combinations of arguments
// numerical values can either be JS numbers, or JmatC complex numbers
// a = integer, b = integer: a*b matrix if 0's
// a = 1D/2D array of numerical values: column vector or 2D matrix
// a = 1D/2D array or undefined, b = 1D/2D array: column vector or 2D matrix with complex values made from real parts of a, imaginary parts of b
// a and b positive integers, var_arg = 1D array or implicit 1D array in arguments: 2D matrix with elements from var_arg, height a, width b. If size of array <= min(a,b), uses it as diagonal instead.
// a = numerical value: 1x1 matrix with that element
// a, b, and var_arg all contain JmatC complex numbers: square matrix with those elements, the square size is Math.ceil(Math.sqrt(arguments.length))
// a = a JmatM matrix object: copy the matrix to a new JmatM
JmatM.make = function(a, b, var_arg) {
  if(a instanceof JmatM) return JmatM.copy(a);

  // Tolerant to all kinds of unexisting array
  // Also supports a 1D array representing an Nx1 2D array
  var softget = function(a, y, x) {
    return (a && a[y]) ? JmatC.cast(a[y][x] == undefined ? a[y] : a[y][x]) : jmatc();
  };
  var softgetr = function(a, y, x) {
    return (a && a[y]) ? JmatR.cast(a[y][x] == undefined ? a[y] : a[y][x]) : 0;
  };
  var softget2 = function(a, b, y, x) {
    return new JmatC(softgetr(a, y, x), softgetr(b, y, x)); // real from a, imag from b
  };
  var arrayw = function(a) {
    if(!a || !a[0]) return 0; // empty array
    if(a[0].length == undefined) return 1; // this means it's a 1D array, such array has witdh 1, not width 0
    return a[0].length;
  };
  // a is array, b is optional second array (only supported for shift == -1)
  // shift -1 ==> a is 2D array
  // shift >= 0 ==> a is 1D array of which the 2D elements are read, first element at a[shift]
  // shift >= 0 and a.length - shift is min(h, w) ==> makes diagonal matrix with those elements on the diagonal instead
  var loop = function(h, w, a, shift, opt_b) {
    var result = new JmatM(h, w);
    if(shift >= 0 && a.length - shift <= Math.min(h, w)) {
      for(var x = 0; x < w && x < h; x++) {
        result.e[x][x] = JmatC.cast(a[x + shift]);
      }
      return result;
    }
    for(var y = 0; y < result.h; y++) {
      for(var x = 0; x < result.w; x++) {
        if(shift < 0) result.e[y][x] = (opt_b ? softget2(a, opt_b, y, x) : softget(a, y, x));
        else result.e[y][x] = JmatC.cast(a[y * w + x + shift]);
      }
    }
    return result;
  };

  // one or two arrays, make all elements from them, size defined by the arrays
  if((a && a.length) || (b && b.length)) {
    var h = Math.max((a && a.length) || 0, (b && b.length) || 0);
    var w = Math.max(arrayw(a), arrayw(b));
    return loop(h, w, a, -1, b);
  }

  // single number or JmatC
  if(a != undefined && b == undefined) {
    var result = new JmatM(1, 1);
    result.e[0][0] = JmatC.cast(a);
    return result;
  }

  // a and b contain dimensions, then elements in array or variable arguments
  if(a != undefined && b != undefined) {
    var h = a;
    var w = b;
    if(var_arg && var_arg.length) return loop(h, w, var_arg, 0);
    return loop(h, w, arguments, 2); // use JS function arguments, shifted by 2 because the first two are a and b
  }

  return new JmatM(0, 0);
};
// Shorcut function for JmatM.make because it's so common
// See comment of Jmat.make for arguments.
/*
here are some examples of how to make matrix, with jmatm JS notation on the left, and a sort of mathematical notation on the right
jmatm(2, 2).toString()                        --> [[0, 0], [0, 0]]
jmatm([[1, 2], [3, 4]]).toString()            --> [[1, 2], [3, 4]]
jmatm(undefined, [[1, 2], [3, 4]]).toString() --> [[1i, 2i], [3i, 4i]]
jmatm([[1, 2], [3, 4]]).mulc(JmatC.I).toString()            --> [[1i, 2i], [3i, 4i]]
jmatm(2, 2, 1, 2, 3, 4).toString()            --> [[1, 2], [3, 4]]
jmatm([1, 2, 3, 4]).toString()                --> [[1],[2],[3],[4]]: a column matrix
jmatm([[1, 2], [3, 4]], [[5, 6], [7, 8]]).toString()                       --> [[1+5i, 2+6i], [3+7i, 4+8i]]
jmatm([[jmatc(1, 5), jmatc(2, 6)], [jmatc(3, 7), jmatc(4, 8)]]).toString() --> [[1+5i, 2+6i], [3+7i, 4+8i]]
jmatm(4, 4, 1, 2, 3, 4).toString()            --> [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0], [0, 0, 0, 4]]: diagonal matrix
jmatm(0).toString()                           --> [0]
*/
jmatm = JmatM.make;

//debugstring
JmatM.toString = function(m) {
  // e.g in console: JmatM.toString(getMatrixFromMem(jmatc(100))
  if(!m) return '' + m;
  var s = '[';
  for(var y = 0; y < m.h; y++) {
    s += '[';
    for(var x = 0; x < m.w; x++) {
      var e = m.e[y][x];
      s += JmatC.toString(e);
      if(x + 1 < m.w) s += ', ';
    }
    s += ']';
    if(y + 1 < m.h) s += ', ';
  }
  return s + ']';
};
JmatM.prototype.toString = function() {
  return JmatM.toString(this);
};


// Does not copy if a is of type JmatM.
JmatM.cast = function(a) {
  return a instanceof JmatM ? a : JmatM.make(a);
};

JmatM.copy = function(a) {
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = jmatc(a.e[y][x].re, a.e[y][x].im);
    }
  }
  return result;
};

// Returns new nxn identity matrix
JmatM.identity = function(n) {
  var r = new JmatM(n, n);
  for(var i = 0; i < n; i++) r.e[i][i] = jmatc(1);
  return r;
};

////////////////////////////////////////////////////////////////////////////////

JmatM.add = function(a, b) {
  if(a.w != b.w || a.h != b.h) return null;
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].add(b.e[y][x]);
    }
  }
  return result;
};
JmatM.prototype.add = function(b) {
  return JmatM.add(this, b);
};

JmatM.sub = function(a, b) {
  if(a.w != b.w || a.h != b.h) return null;
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].sub(b.e[y][x]);
    }
  }
  return result;
};
JmatM.prototype.sub = function(b) {
  return JmatM.sub(this, b);
};

JmatM.mul = function(a, b) {
  if(a.w != b.h) return null;
  var result = new JmatM(a.h, b.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < b.w; x++) {
      var e = jmatc(0);
      for(var z = 0; z < a.w; z++) e = e.add(a.e[y][z].mul(b.e[z][x]));
      result.e[y][x] = e;
    }
  }
  return result;
};
JmatM.prototype.mul = function(b) {
  return JmatM.mul(this, b);
};

// mulScalar (c from complex number)
JmatM.mulc = function(a, s) {
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].mul(s);
    }
  }
  return result;
};
JmatM.prototype.mulc = function(s) {
  return JmatM.mulc(this, s);
};

JmatM.mulr = function(a, s) {
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].mulr(s);
    }
  }
  return result;
};
JmatM.prototype.mulr = function(s) {
  return JmatM.mulr(this, s);
};

// returns a/b = a * b^-1
// In other words, "divides" matrix through matrix
JmatM.div = function(a, b) {
  if(a.w != b.h) return null;
  var result = new JmatM(a.h, b.w);

  b = JmatM.inv(b); //TODO: use pseudo inverse?

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < b.w; x++) {
      var e = jmatc(0);
      for(var z = 0; z < a.w; z++) e = e.add(a.e[y][z].mul(b.e[z][x]));
      result.e[y][x] = e;
    }
  }
  return result;
};
JmatM.prototype.div = function(b) {
  return JmatM.div(this, b);
};

// returns a/b = b^-1 * a
JmatM.leftdiv = function(a, b) {
  if(a.w != b.h) return null;
  var result = new JmatM(a.h, b.w);

  b = JmatM.inv(b); //TODO: use pseudo inverse?

  for(var y = 0; y < b.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = jmatc(0);
      for(var z = 0; z < b.w; z++) e = e.add(b.e[y][z].mul(a.e[z][x]));
      result.e[y][x] = e;
    }
  }
  return result;
};
JmatM.prototype.leftdiv = function(b) {
  return JmatM.leftdiv(this, b);
};

// Divide through complex scalar
JmatM.divc = function(a, s) {
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].div(s);
    }
  }
  return result;
};
JmatM.prototype.divc = function(s) {
  return JmatM.divc(this, s);
};

// divide through real scalar (JS number)
JmatM.divr = function(a, s) {
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].divr(s);
    }
  }
  return result;
};
JmatM.prototype.divr = function(s) {
  return JmatM.divr(this, s);
};

////////////////////////////////////////////////////////////////////////////////
// Categories
////////////////////////////////////////////////////////////////////////////////

JmatM.isReal = function(a) {
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = a.e[y][x];
      if(e.im != 0) return false;
    }
  }
  return true;
};

//valid object, no infinities, no NaN elements
JmatM.isValid = function(a) {
  if(!a || !a.w || !a.h || !a.e || !a.e.length) return false;
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = a.e[y][x];
      if(e) {
        if(e.re == Infinity || e.re == -Infinity || e.re != e.re) return false;
        if(e.im == Infinity || e.im == -Infinity || e.im != e.im) return false;
      } else {
        return false;
      }
    }
  }
  return true;
};

// TODO: functions like isSymmetrical, isHermitian, isDiagonal, ...

////////////////////////////////////////////////////////////////////////////////

JmatM.transpose = function(a) {
  var result = new JmatM(a.w, a.h); //arguments inverted (ctor takes height first)

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[x][y] = a.e[y][x];
    }
  }
  return result;
};
JmatM.prototype.transpose = function() {
  return JmatM.transpose(this);
};

JmatM.neg = function(a) {
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].neg();
    }
  }
  return result;
};
JmatM.prototype.neg = function() {
  return JmatM.neg(this);
};

JmatM.conj = function(a) {
  var result = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = JmatC.conj(a.e[y][x]);
    }
  }
  return result;
};
JmatM.prototype.conj = function() {
  return JmatM.conj(this);
};

//transjugate = transposed conjugate, denoted A* (also called hermitian transpose)
JmatM.transjugate = function(a) {
  return JmatM.conj(JmatM.transpose(a));
};
JmatM.prototype.transjugate = function() {
  return JmatM.transjugate(this);
};

//TODO: LU Decomposition. Implement LUP or similar, regular LU would not work for [0 1][2 3]

// Submatrix with 1 row removed
JmatM.subrow = function(a, row) {
  if(a.h < 2) return null;
  var m = new JmatM(a.h - 1, a.w);
  for(var y = 0; y < a.h - 1; y++) {
    for(var x = 0; x < a.w; x++) {
      m.e[y][x] = a.e[y < row ? y : y + 1][x];
    }
  }
  return m;
};

// Submatrix with 1 column removed
JmatM.subcol = function(a, col) {
  if(a.w < 2) return null;
  var m = new JmatM(a.h, a.w - 1);
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w - 1; x++) {
      m.e[y][x] = a.e[y][x < col ? x : x + 1];
    }
  }
  return m;
};

// The submatrix for a minor, that is, with onw row and one column removed
JmatM.minorsub = function(a, row, col) {
  if(a.h < 2 || a.w < 2) return null;
  var m = new JmatM(a.h - 1, a.w - 1);
  for(var y = 0; y < a.h - 1; y++) {
    for(var x = 0; x < a.w - 1; x++) {
      m.e[y][x] = a.e[y < row ? y : y + 1][x < col ? x : x + 1];
    }
  }
  return m;
};

//submatrix defined by the rectangle y0:y1 x0:x1 (excluding y1 and x1)
//e.g. to get the bottom right 2x2 submatrix of a 5x5 matrix a, use:
//JmatM.submatrix(a, 3, 5, 3, 5)
//note that x and y are 0-indexed, so 5 is outside the matrix
JmatM.submatrix = function(a, y0, y1, x0, x1) {
  if(x0 < 0 || y0 < 0 || x0 > a.w || y0 > a.h) return null;
  if(x1 < 0 || y1 < 0 || x1 > a.w || y1 > a.h) return null;
  var w2 = x1 - x0;
  var h2 = y1 - y0;
  if(w2 <= 0 || h2 <= 0 || w2 > a.w || h2 > a.h) return null;

  var result = new JmatM(h2, w2);

  for(var y = 0; y < h2; y++) {
    for(var x = 0; x < w2; x++) {
      result.e[y][x] = a.e[y0 + y][x0 + x];
    }
  }

  return result;
};

// Requires square matrix
JmatM.minor = function(a, row, col) {
  if(a.h < 2 || a.w < 2 || a.w != a.h) return jmatc(NaN);
  var m = JmatM.minorsub(a, row, col);
  return JmatM.determinant(m);
};

// cofactor: minor with sign depending on alternating position
JmatM.cofactor = function(a, row, col) {
  var m = JmatM.minor(a, row, col);
  var sign = (row + col) % 2 == 0 ? 1 : -1;
  return m.mulr(sign);
};

JmatM.determinant = function(a) {
  if(a.w != a.h) return NaN; //square matrices only

  if(a.w == 1) return a.e[0][0];
  if(a.w == 2) return a.e[0][0].mul(a.e[1][1]).sub(a.e[0][1].mul(a.e[1][0]));

  var result = jmatc(0);

  for(var x = 0; x < a.w; x++) {
    result = result.add(a.e[0][x].mul(JmatM.cofactor(a, 0, x)));
  }

  return result;
};

//Adjugate aka adjoint matrix
JmatM.adj = function(a) {
  if(a.w != a.h) return NaN; //square matrices only
  if(a.w == 1) return JmatM.identity(1);

  //result matrix
  var r = new JmatM(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      //row and column are switched (transpose)
      r.e[y][x] = JmatM.cofactor(a, x, y);
    }
  }

  return r;
};

// Inverse of a matrix
JmatM.inv = function(a) {
  if(a.w != a.h) return null; //square matrices only

  //Cramer's rule
  return JmatM.mulc(JmatM.adj(a), JmatC.inv(JmatM.determinant(a)));
};

//forced pseudoinverse (does not try regular inverse first)
JmatM.pseudoinverse_ = function(a) {
  //TODO: instead of formula below, easier ways are available for scalars, vectors, lin. indep. rows/columns, ...

  var svd = JmatM.svd(a);
  var n = Math.min(svd.s.w, svd.s.h);
  var tolerance = 1e-15; // choose this such that e.g. the pseudoinverse of [[1,2],[1,2]] returns [[0.1,0.2],[0.1,0.2]] - if the tolerance is too small, some near zero gets inverted and the result very wrong (depends on the svd algorithm used of course)
  // Invert all the elements of s, except those that are zero (with some tolerance due to numerical problems)
  // Each element of s should be real, and it only has elements on the diagonal
  for(var i = 0; i < n; i++) svd.s.e[i][i] = (Math.abs(svd.s.e[i][i].re) < tolerance) ? svd.s.e[i][i] : JmatC.inv(svd.s.e[i][i]);
  svd.s = JmatM.transpose(svd.s);

  return JmatM.mul(JmatM.mul(svd.v, svd.s), JmatM.transjugate(svd.u));
};

//Moore-Penrose pseudoinverse: one unique solution for any matrix
JmatM.pseudoinverse = function(a) {
  /*
  Test in console:
  var result = JmatM.pseudoinverse(JmatM.make(2,2,1,2,1,2));
  var orig = JmatM.pseudoinverse(result);
  result.toString() + ' | ' + orig.toString();
  */

  //first try if regular inverse works, if so that is more accurate and faster
  var result = JmatM.inv(a);
  if(JmatM.isValid(result)) return result;

  return JmatM.pseudoinverse_(a);
};

JmatM.getFirstNonZeroDigit_ = function(v) {
  if(v == 0) return 0;
  if(v < 0) v = -v;
  for(var i = 0; i < 100; i++) {
    if(v < 1) v *= 10;
    else if(v >= 10) v /= 10;
    else return Math.floor(v);
  }
  return 0;
};

// Get the matrix representation as a single decimal number. E.g. a [[1,2],[3,4]] would become 22.1234. The dimensions before the comma, the first digit of each value after the comma.
JmatM.getDebugNumber = function(a) {
  var result = a.w + 10*a.h;
  var pos = 0.1;
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var v = JmatM.getFirstNonZeroDigit_(a.e[y][x].re);
      result += v * pos;
      pos /= 10;
    }
  }

  // fix some numerical problem
  pos *= 10;
  result = Math.round(result / pos) * pos;

  return jmatc(result);
};

/*
Matrix norms:
This list is shown here to ensure to not confuse the Frobenius norm with the 2-norm
oo = infinity

p-norms (a.k.a. induced norms or operator norms)
-------
1-norm: maximum absolute column sum of the matrix --> Jmat.maxcolnorm
2-norm: largest singular value, aka spectral norm --> JmatM.norm2
oo-norm: maximum absolute row sum of the matrix --> JmatM.maxrownorm

entrywise norms
---------------
entrywise 1-norm: sum of abs of all the elements --> (not implemented yet)
entrywise 2-norm: Frobenius norm, sqrt of sum of squares of the elements --> JmatM.norm
entrywise oo-norm: maximum of abs of all the elements --> (not implemented yet)

schatten norms
--------------
schatten 1-norm: sum of singular values, aka nuclear norm, trace norm or Ky Fan norm --> (not implemented yet)
schatten 2-norm: sqrt of sum of squares of singular values, results in same value as Frobenius norm --> JmatM.norm
schatten oo-norm: max of the singular values, results in same value as the spectral norm (2-norm) --> JmatM.norm2

*/

//Frobenius norm of the matrix (sqrt of sum of squares of modulus of all elements)
//For a vector, this is the Euclidean norm.
//TODO: since usually the more expensive to calculate 2-norm is meant by "the" norm of the matrix, maybe rename this function to "frobeniusnorm" or "frob"?
JmatM.norm = function(m) {
  var result = 0;
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      result += JmatC.abssqr(m.e[y][x]);
    }
  }
  result = Math.sqrt(result);
  return jmatc(result);
};

//Maximum absolute column sum norm
JmatM.maxcolnorm = function(m) {
  var result = 0;
  for(var x = 0; x < m.w; x++) {
    var current = 0;
    for(var y = 0; y < m.h; y++) {
      current += JmatC.abssqr(m.e[y][x]);
    }
    if (current > result) result = current;
  }
  return jmatc(result);
};

//Maximum absolute row sum norm
JmatM.maxrownorm = function(m) {
  var result = 0;
  for(var y = 0; y < m.h; y++) {
    var current = 0;
    for(var x = 0; x < m.w; x++) {
      current += JmatC.abssqr(m.e[y][x]);
    }
    if (current > result) result = current;
  }
  return jmatc(result);
};

//2-norm: largest singular value (= sqrt of largest eigenvalue of m^H * m). AKA spectral norm
JmatM.norm2 = function(m) {
  var svd = JmatM.svd(m);
  return svd.s.e[0][0];
};

//condition number: largest singular value divided through smallest singular value. Higher ==> more imprecise numerical calculations with this matrix
JmatM.conditionNumber = function(m) {
  var svd = JmatM.svd(m);
  var d = Math.min(m.w, m.h);
  return svd.s.e[0][0].div(svd.s.e[d - 1][d - 1]);
};

//Rank of matrix
JmatM.rank = function(m) {
  // TODO: use the faster RRQR? Or at least svd that returns only s?
  var s = JmatM.svd(m).s;
  var rank = 0;
  for(var i = 0; i < s.w; i++) {
    if(!JmatR.near(s.e[i][i].re, 0, 1e-14)) rank++;
  }
  return jmatc(rank);
};

// Only defined for square matrices
JmatM.trace = function(m) {
  if(m.w != m.h) return jmatc(NaN);
  var result = JmatC.ZERO;
  for(var x = 0; x < m.w; x++) result = result.add(m.e[x][x]);
  return result;
};

// Returns column as h*1 matrix
JmatM.col = function(m, col) {
  var r = new JmatM(m.h, 1);
  for(var y = 0; y < m.h; y++) r.e[y][0] = m.e[y][col];
  return r;
};

// Returns row as 1*w matrix
JmatM.row = function(m, row) {
  var r = new JmatM(1, m.w);
  for(var x = 0; x < m.w; x++) r.e[0][x] = m.e[row][x];
  return r;
};

// Add two non-equal sized matrices.
// b's top left element is at position (row, col) in a (0-indexed)
// so the size of the result matrix is:
// max(row + b.h, a.h) - min(0, row) by max(col + b.w, a.w) - min(0, col)
// any element not overlapped by a or b, will be zero.
JmatM.overlap = function(a, b, row, col) {
  var h = Math.max(row + b.h, a.h) - Math.min(0, row);
  var w = Math.max(col + b.w, a.w) - Math.min(0, col);

  var result = new JmatM(h, w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var rx = col < 0 ? x - col : x;
      var ry = row < 0 ? y - row : y;
      result.e[ry][rx] = a.e[y][x];
    }
  }

  for(var y = 0; y < b.h; y++) {
    for(var x = 0; x < b.w; x++) {
      var rx = col < 0 ? x : x + col;
      var ry = row < 0 ? y : y + row;
      result.e[ry][rx] = result.e[ry][rx].add(b.e[y][x]);
    }
  }

  return result;
};

// QR factorization of complex matrix (with householder transformations)
// m is h*w matrix with h >= w (however, it seems to work correct with h < w too ...)
// returns {q: Q, r: R}
// q is unitary matrix of h*h
// r is h*w upper triangular matrix
JmatM.qr = function(m) {
  /*
  Checks in console:
  var result = JmatM.qr(jmatm(2,2,1,2,3,4));
  JmatM.toString(result.q) + ' \n ' + JmatM.toString(result.r) + ' \n ' + JmatM.toString(JmatM.mul(result.q, result.r));

  var result = JmatM.qr(jmatm(3,3,1,2,3,4,5,6,7,8,9));
  JmatM.toString(result.q) + ' \n ' + JmatM.toString(result.r) + ' \n ' + JmatM.toString(JmatM.mul(result.q, result.r));

  var result = JmatM.qr(jmatm(3,3,12,-51,4,6,167,-68,-4,24,-41));
  JmatM.toString(result.q) + ' \n ' + JmatM.toString(result.r) + ' \n ' + JmatM.toString(JmatM.mul(result.q, result.r));
  */

  //if(m.h < m.w) return null; //seems to work anyway, so don't do this check. TODO: verify this

  var t = Math.min(m.h - 1, m.w);
  var real = JmatM.isReal(m);
  var a = JmatM.copy(m);
  var r = a;
  var q;

  for(var k = 0; k < t; k++) {
    var x = JmatM.col(a, 0);

    var xk = a.e[0][0];
    var normx = JmatM.norm(x);

    var alpha;
    if(xk.im != 0) alpha = JmatC.exp(JmatC.newi(JmatC.argr(xk))).neg().mul(normx);
    else if(xk.re < 0) alpha = normx;
    else alpha = normx.neg();

    var u = JmatM.col(a, 0);
    u.e[0][0] = u.e[0][0].sub(alpha);
    var normu = JmatM.norm(u);

    var v = JmatM.divc(u, normu);

    var vv;
    if(real) {
      vv = JmatM.mulr(JmatM.mul(v, JmatM.transpose(v)), 2);
    } else {
      var xhv = JmatM.mul(JmatM.transjugate(x), v);
      var vhx = JmatM.mul(JmatM.transjugate(v), x);
      var w = xhv.e[0][0].div(vhx.e[0][0]);
      vv = JmatM.mulc(JmatM.mul(v, JmatM.transjugate(v)), w.inc());
    }
    var id = JmatM.identity(a.h);
    var qk = JmatM.sub(id, vv); // here, qk*x = [alpha, 0,...,0]^T

    if (k + 1 < t) {
      a = JmatM.mul(qk, a);
      a = JmatM.minorsub(a, 0, 0);
    }

    if(k == 0) {
      q = JmatM.transjugate(qk);
    } else {
      qk = JmatM.overlap(JmatM.identity(k), qk, k, k);
      q = JmatM.mul(q, JmatM.transjugate(qk)); //Q1^h * Q2^h * ... * Qt^h
    }
    //r = JmatM.mul(qk, r); // Qt * ... * Q2 * Q1 * A //not needed to calculate here, done below instead, is q^t * m
  }
  r = JmatM.mul(JmatM.transjugate(q), m);

  // Solve numerical problem: 0-values of the upper triangular matrix sometimes become a tiny e-16 value
  for(var y = 0; y < r.h; y++) {
    for(var x = 0; x < r.w && x < y; x++) {
      // every value below diagonal should be 0, but leave big ones so that miscalculation bugs can be seen.
      if(JmatC.absr(r.e[y][x]) < 1e-15) r.e[y][x] = jmatc(0);
    }
  }

  return { q: q, r: r };
};

// eigenvalues and vectors of 1x1 matrix
JmatM.eig11 = function(m) {
  if(m.w != 1 || m.h != 1) return null;
  var result = {};
  result.l = new JmatM(1, 1);
  result.v = new JmatM(1, 1);
  result.v.e[0][0] = jmatc(1);
  return result;
};

// explicit algebraic formula for eigenvalues and vectors of 2x2 matrix
JmatM.eig22 = function(m) {
  if(m.w != 2 || m.h != 2) return null;
  var a = jmatc(1);
  var b = m.e[0][0].neg().sub(m.e[1][1]);
  var c = m.e[0][0].mul(m.e[1][1]).sub(m.e[0][1].mul(m.e[1][0]));
  var d = JmatC.sqrt(b.mul(b).sub(a.mul(c).mulr(4)));
  var l1 = b.neg().add(d).div(a.mulr(2));
  var l2 = b.neg().sub(d).div(a.mulr(2));

  var v12 = l1.sub(m.e[0][0]).div(m.e[0][1]);
  var v11 = jmatc(1);
  /*//normalize
  var n = JmatC.sqrt(v12.mul(v12).addr(1));
  v12 = v12.div(n);
  v11 = v11.div(n);*/
  var v22 = l2.sub(m.e[0][0]).div(m.e[0][1]);
  var v21 = jmatc(1);
  /*//normalize
  n = JmatC.sqrt(v12.mul(v12).addr(1));
  v12 = v12.div(n);
  v11 = v11.div(n);*/

  var result = {};
  result.l = new JmatM(2, 1);
  result.v = new JmatM(2, 2);
  result.l.e[0][0] = l1;
  result.l.e[1][0] = l2;
  result.v.e[0][0] = v11;
  result.v.e[1][0] = v12;
  result.v.e[0][1] = v21;
  result.v.e[1][1] = v22;
  return result;
};

// explicit algebraic formula for eigenvalues of 2x2 matrix
JmatM.eigval22 = function(m) {
  if(m.w != 2 || m.h != 2) return null;
  var a = jmatc(1);
  var b = m.e[0][0].add(m.e[1][1]);
  var c = m.e[0][0].mul(m.e[1][1]).sub(m.e[0][1].mul(m.e[1][0]));
  var d = JmatC.sqrt(b.mul(b).sub(a.mul(c).mulr(4)));
  var l1 = b.add(d).div(a.add(a));
  var l2 = b.sub(d).div(a.add(a));
  return [l1, l2];
};

// Returns the eigenvectors and eigenvalues of m as { l: eigenvalues, v: eigenvectors }
// eigenvalues as n*1 column vector, eigenvectors as n*n matrix
// for each column of v and corresponding eigenvalue: A*v = l*v (l represents lambda, A is m)
JmatM.eig = function(m) {
  /*
  Checks in console:
  var result = JmatM.eig(jmatm(2,2,1,2,3,4))
  JmatM.toString(result.l) + ' \n ' + JmatM.toString(result.v);
  */
  if(m.w != m.h || m.w < 1) return null;
  var n = m.w;
  if(n == 1) return JmatM.eig11(m);
  if(n == 2) return JmatM.eig22(m);

  // using the QR algorithm
  var a = JmatM.copy(m); //will contain the eigenvalues on the diagonal

  // Naive implicit QR without shifts, does not work in many cases, left commented out for demo purpose only
  //for(var i = 0; i < 30; i++) {
  //  var qr = JmatM.qr(a);
  //  a = JmatM.mul(qr.r, qr.q); // RQ instead of QR: A_k -> QR, A_(k+1) = RQ
  //}

  // QR with double shifting. This because with single shift or no shift, it does not support complex eigenvalues of real matrix, e.g. [[1,-1],[5,-1]]
  // TODO: this is slow, optimize like with Hessenberg form
  var id = JmatM.identity(n);
  for(var i = 0; i < 15; i++) {
    // var s = a.e[a.h - 1][a.w - 1]; //value that would be chosen for single shift
    // double shift: choose two sigma's, the eigenvalues of the bottom right 2x2 matrix (for which we have the explicit solution)
    var l = JmatM.eigval22(JmatM.submatrix(a, a.h - 2, a.h, a.w - 2, a.w));

    var si0 = id.mulc(l[0]);
    a = a.sub(si0);
    var qr = JmatM.qr(a);
    a = JmatM.mul(qr.r, qr.q).add(si0);

    var si1 = id.mulc(l[1]);
    a = a.sub(si1);
    qr = JmatM.qr(a);
    a = JmatM.mul(qr.r, qr.q).add(si1);
  }

  var v = new JmatM(n, n);

  // Find eigenvectors by solving system of linear equations.
  // TODO: this is not very efficient...
  // Normally, the product of all the qr.q's of the loop above should give the eigenvectors, but that applies only for symmetric matrices while this is supposed to support all
  // So, instead, solve system eqation (A - lamba * I) * x = 0, but with last element of 0 set to 1, and bottom row of (A - lamba * I) set to 0,0,...,0,1.
  // That makes the system solvable, and makes each vector have its last element be 1.
  for(var j = 0; j < n; j++) {
    var value = a.e[j][j];
    var e = JmatM.copy(m); //TODO: this makes it even slower, copy only the needed columns
    for(var i = 0; i < n; i++) e.e[i][i] = e.e[i][i].sub(value);
    for(var i = 0; i < n; i++) e.e[e.h - 1][i] = jmatc(i == n - 1 ? 1 : 0);
    var f = new JmatM(n, 1);
    f.e[f.h - 1][0] = jmatc(1);
    var g = JmatM.solve(e, f);
    for(var i = 0; i < n; i++) v.e[i][j] = g.e[i][0];
  }

  var l = new JmatM(m.w, 1);
  for(var i = 0; i < m.w; i++) l.e[i][0] = a.e[i][i];

  return { l: l, v: v };
};

// Returns the eigen decomposition of m as { v: V, d: D }
// If M is diagonizable, M = V * D * V^(-1)
// In other words: m == result.v.mul(result.d).mul(JmatM.inv(result.v))
// This function is very similar to JmatM.eig. v is the same, d is the same as l but put on the diagonal of a matrix
JmatM.eigd = function(m) {
  var eig = JmatM.eig(m);

  return { v: eig.v, d: JmatM.diag(eig.l) };
};

// Puts all the elements of d in a single diagonal matrix
JmatM.diag = function(d) {
  var n = d.w * d.h;
  var result = new JmatM(n, n);
  var i = 0;
  for(var y = 0; y < d.h; y++) {
    for(var x = 0; x < d.w; x++) {
      result.e[i][i] = d.e[y][x];
      i++;
    }
  }
  return result;
};

// Get element using a one-dimensional index. E.g. a 3x5 matrix has 15 elements, with index 0-14. The index is row by row.
JmatM.get1 = function(m, i) {
  var x = i % m.w;
  var y = Math.floor(i / m.w);
  return m.e[y][x];
};
JmatM.prototype.get1 = function(i) {
  return JmatM.get1(this, i);
}

// Set element using a one-dimensional index. See JmatM.get1.
JmatM.set1 = function(m, i, v) {
  var x = i % m.w;
  var y = Math.floor(i / m.w);
  m.e[y][x] = v;
};

// Cross product between two vectors of length 3, that is, 3x1 and/or 1x3 matrices. Other input dimensions are not accepted.
// Return value has dimensions of the first input (that is, if first input is row vector, it's row vector, otherwise column vector)
// JmatM.cross(jmatm([1,2,3]), jmatm([4,5,6])).toString()
JmatM.cross = function(a, b) {
  if(a.w * a.h != 3 || b.w * b.h != 3) return jmatm(NaN);
  var c = new JmatM(a.h, a.w);
  JmatM.set1(c, 0, a.get1(1).mul(b.get1(2)).sub(a.get1(2).mul(b.get1(1))));
  JmatM.set1(c, 1, a.get1(2).mul(b.get1(0)).sub(a.get1(0).mul(b.get1(2))));
  JmatM.set1(c, 2, a.get1(0).mul(b.get1(1)).sub(a.get1(1).mul(b.get1(0))));
  return c;
};

// Dot product of two vectors.
// Also supports it for matrices of same dimensions (it then is the Frobenius inner product)
// If vectors, row and column vectors may be mixed.
JmatM.dot = function(a, b) {
  if(a.w != b.w || a.h != b.h) {
    if(!(a.w == b.h && a.h == b.w && (a.w == 1 || a.h == 1))) return jmatm(NaN); // Do allow it for differently orientated vectors (h or w is 1)
  }
  var n = a.w * a.h;
  var result = jmatc(0);
  // TODO: use conjugate on b, for complex matrices/vectors?
  for(var i = 0; i < n; i++) result = result.add(a.get1(i).mul(b.get1(i)));
  return result;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Based on LINPACK zsvdc.f, converted to JavaScript.
// Licensing: This is from linpack, from http://www.netlib.org/linpack/
// The license is not mentioned directly in the source code or the website, but
// has been said to now be a variant of the BSD license: see
// https://bugzilla.redhat.com/show_bug.cgi?id=1000829.
// Here is the original author comment from zsvdc.f:
//     linpack. this version dated 03/19/79 .
//              correction to shift calculation made 2/85.
//     g.w. stewart, university of maryland, argonne national lab.
//
// Modifies argument arrays in place.
// Parameters:
//  x: input matrix, n*p, as a 1D array, row-wise
//  ldx: leading dimension of x, ldx >= n
//  n: matrix height
//  p: matrix width
//  s: output, vector of singular values (sorted)
//  e: output, matrix of possible error values
//  u: output, left singular vectors, n*n, as a 1D array, row-wise
//  ldu: leading dimension of u, ldu >= n
//  v: output, right singular vectors, p*p, as a 1D array, row-wise
//  ldv: leading dimension of v, ldv >= p
//  work: scratch array []
//  job: What to calculate. Decimal expansion AB. A=0:no u,1:u,>=2:partial u. B=0:no v,1:v. So to get everything, set job to 11.
// NOTE: "leading dimension" means first dimension of a 2D array, but the arrays are 1D.
//        Normally, leading dimension is the height, or higher if superfluous values were allocated in between.
JmatM.zsvdc_ = function(x, ldx, n, p, s, e, u, ldu, v, ldv, work, job) {
  var i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,
      mm,mm1,mp1,nct,ncu,nrt,info; // integers
  var maxit = 30;
  var t,r; //complex
  var b,c,cs,el,emm1,f,g,dznrm2,scale,shift,sl,sm,sn,
      smm1,t1,test,ztest; // double

  var dr; //drotg result

  var dreal = function(z) { return z.re; };
  var cabs1 = function(z) { return Math.abs(z.re) + Math.abs(z.im); };
  // returns value with absolute value of x, argument of y (transfers sign)
  var csign = function(x, y) { return y.eqr(0) ? jmatc(0) : y.mulr(x.absr() / y.absr()); };
  var sign = function(x, y) { return y == 0 ? 0 : (y < 0 ? -Math.abs(x) : Math.abs(x)); };
  // Euclidean norm of complex vector, n elements starting at index start
  var dznrm2 = function(n, arr, start) {
    var result = jmatc(0);
    for(var i = 0; i < n; i++) {
      var e = arr[start + i];
      result = result.add(e.mul(e));
    }
    return JmatC.sqrt(result);
  };
  // sets vector arry to alpha * arrx + arry
  var zaxpy = function(n, alpha, arrx, startx, arry, starty) {
    for(var i = 0; i < n; i++) {
      arry[starty + i] = arry[starty + i].add(arrx[startx + i].mul(alpha));
    }
  };
  // dot product
  var zdotc = function(n, arrx, startx, arry, starty) {
    var result = jmatc(0);
    for(var i = 0; i < n; i++) {
      result = result.add(arry[starty + i].mul(arrx[startx + i]));
    }
    return result;
  };
  // scales the vector with complex value alpha
  var zscal = function(n, alpha, arr, start) {
    for(var i = 0; i < n; i++) {
      arr[start + i] = arr[start + i].mul(alpha);
    }
  };
  // returns [r, z, c, s]
  var drotg = function(a, b) {
    var a2 = Math.abs(a);
    var b2 = Math.abs(b);
    var r = sign(Math.sqrt(a * a + b * b), (a2 > b2) ? a : b);
    var c = r == 0 ? 1 : a / r;
    var s = r == 0 ? 0 : b / r;
    var z = a2 > b2 ? s : r == 0 ? 0 : c == 0 ? 1 : (1 / c);
    return [r, z, c, s];
  };
  // plane rotation on vectors
  var zdrot = function(n, arrx, startx, arry, starty, c, s) {
    for(var i = 0; i < n; i++) {
      var ax = arrx[startx + i];
      var ay = arrx[starty + i];
      arrx[startx + i] = ax.mulr(c).add(ay.mulr(s));
      arry[starty + i] = ay.mulr(c).sub(ax.mulr(s));
    }
  };
  // swap two vectors
  var zswap = function(n, arrx, startx, arry, starty) {
    for(var i = 0; i < n; i++) {
      var temp = arrx[i + startx];
      arrx[i + startx] = arry[i + starty];
      arry[i + starty] = temp;
    }
  };
  var wantu = false;
  var wantv = false;
  jobu = Math.floor((job % 100) / 10);
  ncu = 1 < jobu ? Math.min(n, p) : n;
  if(jobu != 0) wantu = true;
  if((job % 10) != 0) wantv = true;
  // reduce x to bidiagonal form, storing the diagonal elements
  // in s and the super-diagonal elements in e.
  info = 0;
  nct = Math.min(n - 1, p);
  nrt = Math.max(0, Math.min(p - 2, n));
  lu = Math.max(nct, nrt);
  for(l = 0; l < lu; l++) {
    lp1 = l + 1;
    if(l < nct) {
      // compute the transformation for the l-th column and
      // place the l-th diagonal in s(l).
      s[l] = jmatc(dznrm2(n - l, x, l + l * ldx));
      if(cabs1(s[l]) != 0.0) {
        if(cabs1(x[l + l * ldx]) != 0.0) s[l] = csign(s[l], x[l + l * ldx]);
        t = jmatc(1.0).div(s[l]);
        zscal(n - l, t, x, l + l * ldx);
        x[l + l * ldx] = x[l + l * ldx].addr(1);
      }
      s[l] = s[l].neg();
    }
    for(j = lp1; j < p; j++) {
      if(l < nct) {
        if(cabs1(s[l]) != 0.0) {
          t = zdotc(n - l, x, l + l * ldx, x, l + j * ldx).neg().div(x[l + l * ldx]);
          zaxpy(n - l, t, x, l + l * ldx, x, l + j * ldx);
        }
      }
      // place the l-th row of x into e for the
      // subsequent calculation of the row transformation.
      e[j] = x[l + j * ldx].conj();
    }
    // place the transformation in u for subsequent back multiplication.
    if(wantu && l < nct) {
      for(i = l; i < n; i++) {
        u[i + l * ldu] = x[i + l * ldx];
      }
    }
    if(l < nrt) {
      // compute the l-th row transformation and place the
      // l-th super-diagonal in e(l).
      e[l] = jmatc(dznrm2(p - l - 1, e, lp1));

      if(cabs1(e[l]) != 0.0) {
        if(cabs1(e[lp1]) != 0.0) {
          e[l] = csign(e[l], e[lp1]);
        }
        t = jmatc(1.0) / e[l];
        zscal(p - l - 1, t, e, lp1);
        e[lp1] = jmatc(1.0) + e[lp1];
      }
      e[l] = e[l].conj().neg();
      // apply the transformation.
      if(lp1 < n && cabs1(e[l]) != 0.0) {
        for(j = lp1; j < n; j++) {
          work[j] = jmatc(0.0);
        }
        for(j = lp1; j < p; j++) {
          zaxpy(n - l - 1, e[j], x, lp1 + j * ldx, work, lp1);
        }
        for(j = lp1; j < p; j++) {
          zaxpy(n - l - 1, (-e[j] / e[lp1]).conj(), work, lp1, x, lp1 + j * ldx);
        }
      }
      // place the transformation in v for subsequent back multiplication.
      if(wantv) {
        for(i = lp1; i < p; i++) {
          v[i + l * ldv] = e[i];
        }
      }
    }
  }
  // set up the final bidiagonal matrix of order m.
  m = Math.min(p, n + 1);
  if(nct < p) s[nct] = x[nct + nct * ldx];
  if(n < m) s[m - 1] = jmatc(0.0);
  if(nrt + 1 < m) e[nrt] = x[nrt + (m - 1) * ldx];
  e[m - 1] = jmatc(0.0);
  // if required, generate u.
  if(wantu) {
    for(j = nct; j < ncu; j++) {
      for(i = 0; i < n; i++) {
        u[i + j * ldu] = jmatc(0.0);
      }
      u[j + j * ldu] = jmatc(1.0);
    }
    for(ll = 0; ll < nct; ll++) {
      l = nct - ll - 1;
      if(cabs1(s[l]) != 0.0) {
        lp1 = l + 1;
        for(j = lp1; j < ncu; j++) {
          t = zdotc(n - l, u, l + l * ldu, u, l + j * ldu).neg().div(u[l + l * ldu]);
          zaxpy(n - l, t, u, l + l * ldu, u, l + j * ldu);
        }
        zscal(n - l, jmatc(-1.0), u, l + l * ldu);
        u[l + l * ldu] = u[l + l * ldu].inc();
        for(i = 0; i < l; i++) {
          u[i + l * ldu] = jmatc(0.0);
        }
      } else {
        for(i = 0; i < n; i++) {
          u[i + l * ldu] = jmatc(0.0);
        }
        u[l + l * ldu] = jmatc(1.0);
      }
    }
  }
  // if it is required, generate v.
  if(wantv) {
    for(ll = 0; ll < p; ll++) {
      l = p - ll - 1;
      lp1 = l + 1;
      if(l < nrt) {
        if(cabs1(e[l]) != 0.0) {
          for(j = lp1; j < p; j++) {
            t = zdotc(p - l + 1, v, lp1 - 1 + l * ldv, v, lp1 + j * ldv).neg().div(v[lp1 + l * ldv]);
            zaxpy(p - l + 1, t, v, lp1 + l * ldv, v, lp1 + j * ldv);
          }
        }
      }
      for(i = 0; i < p; i++) {
        v[i + l * ldv] = jmatc(0.0);
      }
      v[l + l * ldv] = jmatc(1.0);
    }
  }
  // transform s and e so that they are real.
  for(i = 0; i < m; i++) {
    if(cabs1(s[i]) != 0.0) {
      t = jmatc(JmatC.absr(s[i]));
      r = s[i].div(t);
      s[i] = t;
      if(i + 1 < m) e[i] = e[i].div(r);
      if(wantu) zscal(n, r, u, i * ldu);
    }
    if(i + 1 == m) break;
    if(cabs1(e[i]) != 0.0) {
      t = jmatc(JmatC.absr(e[i]));
      r = t.div(e[i]);
      e[i] = t;
      s[i + 1] = s[i + 1].mul(r);
      if(wantv) zscal(p, r, v, (i + 1) * ldv);
    }
  }
  // main iteration loop for the singular values.
  mm = m;
  iter = 0;
  for(;;) {
    // quit if all the singular values have been found.
    if(m == 0) break;
    // if too many iterations have been performed, set flag and return.
    if(maxit <= iter) {
      info = m;
      break;
    }
    // this section of the program inspects for
    // negligible elements in the s and e arrays.  on
    // completion the variables kase and l are set as follows.
    //
    // kase = 1     if s(m) and e(l - 1) are negligible and l.lt.m
    // kase = 2     if s(l) is negligible and l.lt.m
    // kase = 3     if e(l - 1) is negligible, l.lt.m, and
    //              s(l), ..., s(m) are not negligible (qr step).
    // kase = 4     if e(m - 1) is negligible (convergence).
    for(ll = 1; ll <= m; ll++) {
      l = m - ll;
      if(l == 0) break;
      test = JmatC.absr(s[l - 1]) + JmatC.absr(s[l]);
      ztest = test + JmatC.absr(e[l - 1]);
      if(ztest == test) {
        e[l - 1] = jmatc(0.0);
        break;
      }
    }
    if(l == m - 1) {
      kase = 4;
    } else {
      lp1 = l + 1;
      mp1 = m + 1;
      for(lls = lp1; lls <= mp1; lls++) {
        ls = m - lls + lp1;
        if(ls == l) break;
        test = 0.0;
        if(ls != m) test = test + JmatC.absr(e[ls - 1]);
        if(ls != l + 1) test = test + JmatC.absr(e[ls - 2]);
        ztest = test + JmatC.absr(s[ls - 1]);
        if(ztest == test) {
          s[ls - 1] = jmatc(0.0);
          break;
        }
      }
      if(ls == l) {
        kase = 3;
      } else if(ls == m) {
        kase = 1;
      } else {
        kase = 2;
        l = ls;
      }
    }
    l++;
    // deflate negligible s(m).
    if(kase == 1) {
      mm1 = m - 1;
      f = dreal(e[m - 2]);
      e[m - 2] = jmatc(0.0);
      for(kk = 1; kk <= mm1; kk++) {
        k = mm1 - kk + l;
        t1 = dreal(s[k - 1]);
        dr = drotg(t1, f); t1 = dr[0]; f = dr[1]; cs = dr[2]; sn = dr[3];
        s[k - 1] = jmatc(t1);
        if(k != l) {
          f = -sn * dreal(e[k - 2]);
          e[k - 2] = cs * e[k - 2];
        }
        if(wantv) zdrot(p, v, (k - 1) * ldv, v, (m - 1) * ldv, cs, sn);
      }
    }
    // split at negligible s(l).
    else if(kase == 2) {
      f = dreal(e[l - 2]);
      e[l - 2] = jmatc(0.0);
      for(k = l; k <= m; k++) {
        t1 = dreal(s[k - 1]);
        dr = drotg(t1, f); t1 = dr[0]; f = dr[1]; cs = dr[2]; sn = dr[3];
        s[k - 1] = jmatc(t1);
        f = -sn * dreal(e[k - 1]);
        e[k - 1] = e[k - 1].mulr(cs);
        if(wantu) zdrot(n, u, (k - 1) * ldu, u, (l - 2) * ldu, cs, sn);
      }
    }
    // perform one qr step.
    else if(kase == 3) {
      // calculate the shift.
      scale = Math.max(Math.max(Math.max(Math.max(JmatC.absr(s[m - 1]),
              JmatC.absr(s[m - 2])), JmatC.absr(e[m - 2])), JmatC.absr(s[l - 1])),
              JmatC.absr(e[l - 1]));
      sm = dreal(s[m - 1]) / scale;
      smm1 = dreal(s[m - 2]) / scale;
      emm1 = dreal(e[m - 2]) / scale;
      sl = dreal(s[l - 1]) / scale;
      el = dreal(e[l - 1]) / scale;
      b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
      c = (sm * emm1) * (sm * emm1);
      shift = 0.0;
      if(b != 0.0 || c != 0.0) {
        shift = Math.sqrt(b * b + c);
        if(b < 0.0) shift = -shift;
        shift = c / (b + shift);
      }
      f =(sl + sm) * (sl - sm) + shift;
      g = sl * el;
      //  chase zeros.
      mm1 = m - 1;
      for(k = l; k <= mm1; k++) {
        dr = drotg(f, g); f = dr[0]; g = dr[1]; cs = dr[2]; sn = dr[3];
        if(k != l) e[k - 2] = jmatc(f);
        f = cs * dreal(s[k - 1]) + sn * dreal(e[k - 1]);
        e[k - 1] = e[k - 1].mulr(cs).sub(s[k - 1].mulr(sn));
        g = sn * dreal(s[k]);
        s[k] = s[k].mulr(cs);
        if(wantv) zdrot(p, v, (k - 1) * ldv, v, k * ldv, cs, sn);
        dr = drotg(f, g); f = dr[0]; g = dr[1]; cs = dr[2]; sn = dr[3];
        s[k - 1] = jmatc(f);
        f = cs * dreal(e[k - 1]) + sn * dreal(s[k]);
        s[k] = e[k - 1].mulr(-sn).add(s[k].mulr(cs));
        g = sn * dreal(e[k]);
        e[k] = e[k].mulr(cs);
        if(wantu && k < n) zdrot(n, u, (k - 1) * ldu, u, k * ldu, cs, sn);
      }
      e[m - 2] = jmatc(f);
      iter++;
    } else if(kase == 4) {
      // convergence.
      // make the singular value positive.
      if(dreal(s[l - 1]) < 0.0) {
        s[l - 1] = s[l - 1].neg();
        if(wantv) zscal(p, jmatc(-1.0), v, (l - 1) * ldv);
      }
      // order the singular values.
      while(l != mm) {
        if(dreal(s[l]) <= dreal(s[l - 1])) break;
        t = s[l - 1];
        s[l - 1] = s[l];
        s[l] = t;
        if(wantv && l < p) zswap(p, v, (l - 1) * ldv, v, l * ldv);
        if(wantu && l < n) zswap(n, u, (l - 1) * ldu, u, l * ldu);
        l = l++;
      }
      iter = 0;
      m--;
    }
  }
  return info;
};
// End of Linpack zsvdc
////////////////////////////////////////////////////////////////////////////////

// Singular value decomposition with Matrix objects
// input M, returns {u: U, s: S, v: V } such that U * W * V^T = M and S diagonal with singular values (^T means conjugate transpose here)
// Input allowed to be non-square. The size of "S" is same as the input matrix.
JmatM.svd = function(m) {
  /*
  Checks in console:
  function testSvd(m) {
    var result = JmatM.svd(jmatm(m));
    console.log(Jmat.dbg_(result) + ' | ' + JmatM.mul(JmatM.mul(result.u, result.s), JmatM.transpose(result.v)).toString());
  }
  testSvd(jmatm(2,2,1,2,3,4))
  testSvd(jmatm(2,2,1,2,3,4).mulc(JmatC.I))
  testSvd(jmatm(2,2,1,2,1,2))
  testSvd(jmatm([[1,2]]))
  testSvd(jmatm([[1],[2]]))
  */

  // 1D array representing the matrix
  var a = [];
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      a[y + x * m.h] = m.e[y][x];
    }
  }

  var s = []; //h*w
  var u = []; //h*h
  var v = []; //w*w
  var e = []; //h*w

  // function(x, ldx, n, p, s, e, u, ldu, v, ldv, job)
  JmatM.zsvdc_(a, m.h, m.h, m.w, s, e, u, m.h, v, m.w, [], 11);
  
  // Solve numerical problems: singular values < eta should be 0
  var eta = 1e-15; //TODO: check if this tolerance can be improved 
  for(i = 0; i < s.length; i++) if(Math.abs(s[i]) < eta) s[i] = 0;

  var result = { u: new JmatM(m.h, m.h), s: new JmatM(m.h, m.w), v: new JmatM(m.w, m.w) };


  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.h; x++) {
      result.u.e[y][x] = u[y + x * m.h];
    }
  }
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      result.s.e[y][x] = (x == y) ? jmatc(s[x]) : jmatc(0);
    }
  }
  for(var y = 0; y < m.w; y++) {
    for(var x = 0; x < m.w; x++) {
      result.v.e[y][x] = v[y + x * m.w];
    }
  }

  return result;
};

// equals
JmatM.eq = function(a, b) {
  if(a.w != b.w || a.h != b.h) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      return a.e[y][x].eq(b.e[y][x]);
    }
  }

  return true;
};

// nearly equal
JmatM.near = function(a, b, epsilon) {
  if(a.w != b.w || a.h != b.h) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var ea = a.e[y][x];
      var eb = b.e[y][x];
      if(Math.abs(ea.re - eb.re) > epsilon) return false;
      if(Math.abs(ea.im - eb.im) > epsilon) return false;
    }
  }

  return true;
};


//solves system of linear equations ax=b.
//Returns null if the system is inconsistent and has no solution, x otherwise.
//If multiple solutions are possible, returns the solution where the vector of free variables is 0.
//Uses the pseudoinverse if A is not invertible
//a: input matrix, h*w size
//b: input vector, h*1 size
JmatM.solve = function(a, b) {
  if(a.h != b.h) return null;
  var ag = JmatM.pseudoinverse(a);

  aag = JmatM.mul(a, ag);
  var aagb = JmatM.mul(aag, b);
  if(!JmatM.near(b, aagb)) return null; //inconsistent system with no solution

  return JmatM.mul(ag, b);
};


////////////////////////////////////////////////////////////////////////////////
/*
License of the kiss_ and kf_ functions below (converted from C to JavaScript):
Kiss FFT
Copyright (c) 2003-2010, Mark Borgerding
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
// Constructor
JmatM.kiss_fft_state_ = function(nfft, inverse) {
  this.nfft = nfft;
  this.inverse = inverse;
  this.factors = []; //int array, size 32
  this.twiddles = []; //complex JmatC array, size nfft-1
};
JmatM.kf_bfly2_ = function(Fout /*array of complex JmatC*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
  var j = 0;
  for(var i = 0; i < m; i++) {
    var t = Fout[Fout_index + i + m].mul(st.twiddles[j]);
    j += fstride;
    Fout[Fout_index + i + m] = Fout[Fout_index + i].sub(t);
    Fout[Fout_index + i] = Fout[Fout_index + i].add(t);
  }
};
JmatM.kf_bfly4_ = function(Fout /*array of complex JmatC*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
  var scratch = []; //size 6
  var m2=2*m;
  var m3=3*m;
  var j1 = 0;
  var j2 = 0;
  var j3 = 0;
  for(var i = 0; i < m; i++) {
    scratch[0] = Fout[Fout_index + i + m].mul(st.twiddles[j1]);
    scratch[1] = Fout[Fout_index + i + m2].mul(st.twiddles[j2]);
    scratch[2] = Fout[Fout_index + i + m3].mul(st.twiddles[j3]);
    scratch[5] = Fout[Fout_index + i].sub(scratch[1]);
    Fout[Fout_index + i] = Fout[Fout_index + i].add(scratch[1]);
    scratch[3] = scratch[0].add(scratch[2]);
    scratch[4] = scratch[0].sub(scratch[2]);
    Fout[Fout_index + i + m2] = Fout[Fout_index + i].sub(scratch[3]);
    j1 += fstride;
    j2 += fstride*2;
    j3 += fstride*3;
    Fout[Fout_index + i] = Fout[Fout_index + i].add(scratch[3]);
    if(st.inverse) {
      Fout[Fout_index + i + m].re = scratch[5].re - scratch[4].im;
      Fout[Fout_index + i + m].im = scratch[5].im + scratch[4].re;
      Fout[Fout_index + i + m3].re = scratch[5].re + scratch[4].im;
      Fout[Fout_index + i + m3].im = scratch[5].im - scratch[4].re;
    } else {
      Fout[Fout_index + i + m].re = scratch[5].re + scratch[4].im;
      Fout[Fout_index + i + m].im = scratch[5].im - scratch[4].re;
      Fout[Fout_index + i + m3].re = scratch[5].re - scratch[4].im;
      Fout[Fout_index + i + m3].im = scratch[5].im + scratch[4].re;
    }
  }
};
JmatM.kf_bfly3_ = function(Fout /*array of complex JmatC*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
  var k=m;
  var m2 = 2*m;
  var j1 = 0;
  var j2 = 0;
  var scratch = [];
  var epi3 = st.twiddles[fstride*m];
  for(var i = 0; i < k; i++) {
    scratch[1]=Fout[Fout_index + i+m].mul(st.twiddles[j1]);
    scratch[2]=Fout[Fout_index + i+m2].mul(st.twiddles[j2]);
    scratch[3]=scratch[1].add(scratch[2]);
    scratch[0]=scratch[1].sub(scratch[2]);
    j1 += fstride;
    j2 += fstride*2;
    Fout[Fout_index + i+m].re = Fout[Fout_index + i].re - scratch[3].re/2;
    Fout[Fout_index + i+m].im = Fout[Fout_index + i].im - scratch[3].im/2;
    scratch[0] = scratch[0].mulr(epi3.im);
    Fout[Fout_index + i] = Fout[Fout_index + i].add(scratch[3]);
    Fout[Fout_index + i+m2].re = Fout[Fout_index + i+m].re + scratch[0].im;
    Fout[Fout_index + i+m2].im = Fout[Fout_index + i+m].im - scratch[0].re;
    Fout[Fout_index + i+m].re -= scratch[0].im;
    Fout[Fout_index + i+m].im += scratch[0].re;
  }
};
JmatM.kf_bfly5_ = function(Fout /*array of complex JmatC*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
  var scratch = []; //size-13 complex array
  var ya = st.twiddles[fstride*m];
  var yb = st.twiddles[fstride*2*m];
  var m2 = 2 * m;
  var m3 = 3 * m;
  var m4 = 4 * m;
  for (var u=0; u<m; ++u ) {
    scratch[0] = jmatc(Fout[Fout_index + u]);
    scratch[1] = Fout[Fout_index + m+u].mul(st.twiddles[u*fstride]);
    scratch[2] = Fout[Fout_index + m2+u].mul(st.twiddles[2*u*fstride]);
    scratch[3] = Fout[Fout_index + m3+u].mul(st.twiddles[3*u*fstride]);
    scratch[4] = Fout[Fout_index + m4+u].mul(st.twiddles[4*u*fstride]);
    scratch[7] = scratch[1].add(scratch[4]);
    scratch[10]= scratch[1].sub(scratch[4]);
    scratch[8] = scratch[2].add(scratch[3]);
    scratch[9] = scratch[2].sub(scratch[3]);
    Fout[Fout_index + u].re += scratch[7].re + scratch[8].re;
    Fout[Fout_index + u].im += scratch[7].im + scratch[8].im;
    scratch[5] = jmatc(0);
    scratch[5].re = scratch[0].re + scratch[7].re*ya.re + scratch[8].re*yb.re;
    scratch[5].im = scratch[0].im + scratch[7].im*ya.re + scratch[8].im*yb.re;
    scratch[6] = jmatc(0);
    scratch[6].re = scratch[10].im*ya.im + scratch[9].im*yb.im;
    scratch[6].im = -scratch[10].re*ya.im - scratch[9].re*yb.im;
    Fout[Fout_index + m+u]=scratch[5].sub(scratch[6]);
    Fout[Fout_index + m4+u]=scratch[5].add(scratch[6]);
    scratch[11] = jmatc(0);
    scratch[11].re = scratch[0].re + scratch[7].re*yb.re + scratch[8].re*ya.re;
    scratch[11].im = scratch[0].im + scratch[7].im*yb.re + scratch[8].im*ya.re;
    scratch[12] = jmatc(0);
    scratch[12].re = -scratch[10].im*yb.im + scratch[9].im*ya.im;
    scratch[12].im = scratch[10].re*yb.im - scratch[9].re*ya.im;
    Fout[Fout_index + m2+u]=scratch[11].add(scratch[12]);
    Fout[Fout_index + m3+u]=scratch[11].sub(scratch[12]);
  }
};
// perform the butterfly for one stage of a mixed radix FFT
JmatM.kf_bfly_generic_ = function(Fout /*array of complex JmatC*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/, p /*int*/) {
  var u,k,q1,q; /*int*/
  var t; // complex JmatC
  var Norig = st.nfft;
  var scratch = [];
  for ( u=0; u<m; ++u ) {
    k=u;
    for ( q1=0 ; q1<p ; ++q1 ) {
      scratch[q1] = jmatc(Fout[Fout_index + k]);
      k += m;
    }
    k=u;
    for ( q1=0 ; q1<p ; ++q1 ) {
      var twidx=0;
      Fout[Fout_index + k] = scratch[0];
      for (q=1;q<p;++q ) {
        twidx += fstride * k;
        if (twidx>=Norig) twidx-=Norig;
        t = scratch[q].mul(st.twiddles[twidx] );
        Fout[Fout_index + k] = Fout[Fout_index + k].add(t);
      }
      k += m;
    }
  }
};
JmatM.kf_work_ = function(Fout /*array of complex JmatC*/, Fout_index, f /*array of complex JmatC*/,f_index,
    fstride /*int*/, in_stride /*int*/, factors /*int array*/, factors_index,st /*kiss_fft_state*/) {
  var p = factors[factors_index + 0]; /* the radix */
  var m = factors[factors_index + 1]; /* stage's fft length/p */
  var j = 0;

  if (m==1) {
    for(var i = 0; i < p*m; i++) {
      Fout[i + Fout_index] = jmatc(f[f_index + j]);
      j += fstride*in_stride;
    }
  }else{
    for(var i = 0; i < p*m; i += m) {
      // recursive call:
      // DFT of size m*p performed by doing p instances of smaller DFTs of size m, each one takes a decimated version of the input
      JmatM.kf_work_(Fout, Fout_index + i, f, f_index + j, fstride*p, in_stride, factors, factors_index + 2, st);
      j += fstride*in_stride;
    }
  }
  // recombine the p smaller DFTs
  switch (p) {
    case 2: JmatM.kf_bfly2_(Fout,Fout_index,fstride,st,m); break;
    case 3: JmatM.kf_bfly3_(Fout,Fout_index,fstride,st,m); break;
    case 4: JmatM.kf_bfly4_(Fout,Fout_index,fstride,st,m); break;
    case 5: JmatM.kf_bfly5_(Fout,Fout_index,fstride,st,m); break;
    default: JmatM.kf_bfly_generic_(Fout,Fout_index,fstride,st,m,p); break;
  }
};
// facbuf is populated by p1,m1,p2,m2, ... where p[i] * m[i] = m[i-1], m0 = n 
JmatM.kf_factor_ = function(n, facbuf) {
  var i = 0;
  var p=4;
  var floor_sqrt = Math.floor( Math.sqrt(n) );
  // factor out powers of 4, powers of 2, then any remaining primes
  do {
    while (n % p != 0) {
      switch (p) {
        case 4: p = 2; break;
        case 2: p = 3; break;
        default: p += 2; break;
      }
      if (p > floor_sqrt) p = n; // no more factors, skip to end
    }
    n = Math.floor(n / p);
    facbuf[i + 0] = p;
    facbuf[i + 1] = n;
    i += 2;
  } while (n > 1);
};
// returns kiss_fft_state object initialized for given size and inversion
JmatM.kiss_fft_alloc_ = function(nfft,inverse_fft) {
  var st = new JmatM.kiss_fft_state_(nfft, inverse_fft);
  for (var i=0;i<nfft;++i) {
    var phase = -Math.PI*2*i / nfft;
    if (st.inverse) phase *= -1;
    st.twiddles[i] = new JmatC(Math.cos(phase), Math.sin(phase));
  }
  JmatM.kf_factor_(nfft,st.factors);
  return st;
};
JmatM.kiss_fft_ = function(st/*kiss_fft_state object*/,fin/*complex JmatC array of size nfft*/,fout/*complex JmatC array of size nfft*/) {
  JmatM.kf_work_(fout,0, fin,0, 1, 1/*in_stride*/, st.factors,0,st);
};
// End of Kiss FFT
////////////////////////////////////////////////////////////////////////////////


JmatM.matrixfft_ = function(m, inverse) {
  var rowresult = new JmatM(m.h, m.w);

  // apply to each row
  if(m.w > 1) {
    for(var j = 0; j < m.h; j++) {
      var out = [];
      for(var i = 0; i < m.w; i++) out[i] = jmatc(0);
      var st = JmatM.kiss_fft_alloc_(m.w, inverse);
      JmatM.kiss_fft_(st, m.e[j], out);
      for(var i = 0; i < m.w; i++) rowresult.e[j][i] = out[i];
    }
  } else {
    rowresult = m;
  }

  var result = new JmatM(m.h, m.w);

  // apply to each column
  if (m.h > 1) {
    for(var j = 0; j < m.w; j++) {
      var col = JmatM.transpose(JmatM.col(rowresult, j));
      var out = [];
      for(var i = 0; i < m.h; i++) out[i] = jmatc(0);
      var st = JmatM.kiss_fft_alloc_(m.h, inverse);
      JmatM.kiss_fft_(st, col.e[0], out);
      for(var i = 0; i < m.h; i++) result.e[i][j] = out[i];
    }
  } else {
    result = rowresult;
  }

  var factor = 1.0 / Math.sqrt(m.w * m.h);
  result = JmatM.mulr(result, factor);

  return result;
};

// Discrete fourier transform
JmatM.fft = function(m) {
  // Debug in console:
  // JmatM.toString(JmatM.fft(jmatm(2,2,1,2,3,4)))
  // should give [5 -1][-2 0]
  return JmatM.matrixfft_(m, 0);
};

// Inverse discrete fourier transform
JmatM.ifft = function(m) {
  return JmatM.matrixfft_(m, 1);
};


//Matrix exponential
JmatM.exp = function(m) {
  if(m.h != m.w) return null; //must be square

  var result = m.add(JmatM.identity(m.w));
  var mm = m;
  var k = 1;

  // TODO: quit early if convergence reached
  for(var i = 2; i <= 20; i++) {
    k *= i;
    mm = mm.mul(m);
    result = result.add(mm.mulr(1 / k));
  }

  return result;
};

/*
Matrix cosine

debug in console:
var a = jmatm(2,2,1,2,3,4);
var c = JmatM.cos(a); var s = JmatM.sin(a);
JmatM.toString(c) + ' ' + JmatM.toString(s) + ' ' + JmatM.toString(c.mul(c).add(s.mul(s)));
--> c*c+s*s should give identity matrix
*/
JmatM.cos = function(m) {
  if(m.h != m.w) return null; //must be square

  var result = JmatM.identity(m.w);
  var mm = m.mul(m);
  var mmm = null;
  var k = 1;
  var sign = 1;

  // TODO: quit early if convergence reached
  for(var i = 0; i < 20; i++) {
    if(i == 0) {
      k = 2;
    } else {
      k *= (i * 2 + 1) * (i * 2 + 2); //e.g. when i = 1, k is (4!)
    }

    sign = -sign;
    mmm = (mmm == null) ? mm : mmm.mul(mm);
    result = result.add(mmm.mulr(sign / k));
  }

  return result;
};

// Matrix sine
JmatM.sin = function(m) {
  if(m.h != m.w) return null; //must be square

  var result = m;
  var mm = m.mul(m);
  var mmm = m;
  var k = 1;
  var sign = 1;

  // TODO: quit early if convergence reached
  for(var i = 0; i < 20; i++) {
    k *= (i * 2 + 2) * (i * 2 + 3); //e.g. when i = 1, k is (5!)

    sign = -sign;
    mmm = mmm.mul(mm);
    result = result.add(mmm.mulr(sign / k));
  }

  return result;
};

// Square root of matrix (returns a such that a * a = m)
// debug in console: var a = JmatM.sqrt(jmatm(3,3,1,1,0,0,0,1,1,0,1)); JmatM.toString(a) + ' ' + JmatM.toString(a.mul(a))
JmatM.sqrt = function(m) {
  if(m.h != m.w) return null; //must be square

  // Babylonian method. Does not work for [[1,2][3,4]], because that one has complex result. Left commented out for demonstration purpose only.
  /*var result = m.add(JmatM.identity(m.w)).mulr(0.5);
  for(var i = 0; i <= 30; i++) {
    result = result.add(m.div(result)).mulr(0.5);
  }
  return result;*/

  // With eigen decomposition: only the sqrt of the diagonals of D needs to be taken
  // TODO: this will only work for diagonizable matrix. Implement a way for non-diagonizable one as well (with Jordan normal form)
  var e = JmatM.eigd(m);
  var v = e.v;
  var d = e.d;
  for(var i = 0; i < d.w; i++) d.e[i][i] = JmatC.sqrt(d.e[i][i]);
  return v.mul(d).mul(JmatM.inv(v));
};

// Matrix logarithm (base e, ln)
// debug in console: var a = JmatM.log(jmatm(2,2,1,2,3,4)); JmatM.toString(a) + ' ' + JmatM.toString(JmatM.exp(a))
JmatM.log = function(m) {
  if(m.h != m.w) return null; //must be square

  // With eigen decomposition: only the log of the diagonals of D needs to be taken
  // TODO: this will only work for diagonizable matrix. Implement a way for non-diagonizable one as well (with Jordan normal form)
  var e = JmatM.eigd(m);
  var v = e.v;
  var d = e.d;
  for(var i = 0; i < d.w; i++) d.e[i][i] = JmatC.log(d.e[i][i]);
  return v.mul(d).mul(JmatM.inv(v));
};

// Matrix to any complex scalar power s
JmatM.powc = function(m, s) {
  if(m.h != m.w) return null; //must be square

  // With eigen decomposition: only the log of the diagonals of D needs to be taken
  // TODO: this will only work for diagonizable matrix. Implement a way for non-diagonizable one as well (with Jordan normal form)
  var e = JmatM.eigd(m);
  var v = e.v;
  var d = e.d;
  for(var i = 0; i < d.w; i++) d.e[i][i] = d.e[i][i].pow(s);
  return v.mul(d).mul(JmatM.inv(v));
};
