/*
Jmat.js

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

Jmat.js is a numerical library in JavaScript for complex, matrix and statistical
arithmetic.

--------------------------------------------------------------------------------

Most of the code is written based on formulas from Wikipedia, Mathworld,
NIST DLMF, various papers, and books "Handbook of Mathematical Functions" and
"Numerical Linear Algebra".

Two algorithms use code from third party open source libraries, converted to
JavaScript. Their licenses and attribution are included at the implementation.
It are KISS FFT (see Jmat.Matrix.kiss_...) and linpack zsvdc (see
Jmat.Matrix.zsvdc_). They are included in this source file to keep it self
contained.

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

Usage
-----

The mathematical types introduced are:
 - Complex or Jmat.Complex: complex number, contains re and im fields. In addition has lots of static functions.
 - Matrix or Jmat.Matrix: complex matrix, its fields are: h: number of rows, w: number of columns, e: 2D array of elements (row per row). In addition has lots of static functions.

The namespaces introduced are:
 - Jmat: main API, wrapper around functions defined in Real, Complex and Matrix
 - Real or Jmat.Real: contains some functions that operate on primitive JS numbers. Most math however, even real math, is instead in Complex.

See the Jmat.#### function definitions at the beginning of the source code to see a list of most functions.
These are however wrappers around Jmat.Real, Jmat.Complex and Jmat.Matrix functions. See further below in the code to find
their implementation. Some other less exposed functions can also be found further down in Jmat.Real, Jmat.Complex and Jmat.Matrix,
but these are less stable API.

The rest of this manual shows usage examples:

*) Create a complex number, e.g. 1+2i:
var complex = new Complex(1, 2)
shorthand (no new and more argument types supported):
var complex = Complex(1, 2)
var complex = Complex('1+2i')
its fields are re and im:
complex.re --> 1
complex.im --> 2

Most math functions in Complex always require complex number objects (not JS numbers), even if they are all real. Here are the ways and shortcuts of making a real number "3":
var real = new Complex(3, 0)
var real = Complex(3)

Unfortunately, JavaScript does not have operator overloading, so simply using "+" to add two complex numbers does not work and the "add" function is needed instead.

*) Add two complex numbers:
Complex(1, 2).add(Complex(3, 4))
alternatives:
Complex.add(Complex(1, 2), Complex(3, 4))
Jmat.add(Complex(1, 2), Complex(3, 4))
other basic operators: sub, mul, div, pow

*) Adding a real number to a complex number:
var complex = Complex(1, 2)
complex.addr(3) --> 4+2i. Prototype member function of Complex that takes real JS number
complex.add(Complex(3)) --> 4+2i. Prototype member function of Complex that takes Complex object
Complex.addr(complex, 3) --> 4+2i. Complex.addr requires first argument of type Complex, second as primitive JS number
Complex.add(complex, Complex(3)) --> 4+2i. Complex.add requires both arguments to be exactly of type Complex
Jmat.add(complex, 3) --> 4+2i. Jmat.add supports various input types and casts internally

Similar for sub, mul, div, pow, ...

*) Output Complex object as more readable string:
Complex(1, 2).add(Complex(3, 4)).toString()
--> 4+6i

*) Complex trigonometrics
var z = Complex(1, 2);
var s = Complex.sin(z);
var c = Complex.cos(z);
s.mul(s).add(c.mul(c)).toString(8) --> 1. This is cos^2(z) + sin^2(z), the '8' argument rounds to 8 digits to hide some numerical imprecision
Jmat.acos(5).toString() --> 2.2924316695611733i. Arccos of large argument gives imaginary number.


*) Gamma function of 5:
Complex.gamma(Complex(5))
more convenient (but slightly less efficient):
Jmat.gamma(5) --> The Jmat wrapper automatically casts to Complex if regular JS number is given

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
var matrix = new Matrix(2, 2); // num rows, nuw columns
matrix.e[0][0] = Complex(1);
matrix.e[0][1] = Complex(2);
matrix.e[1][0] = Complex(3);
matrix.e[1][1] = Complex(4);
shorthand for all the above:
var matrix = Matrix([[1, 2], [3, 4]]) // 2D array, array of rows
var matrix = Matrix(2, 2, 1, 2, 3, 4) // height, width, 4 elements
print it:
matrix.toString()
--> [[1, 2], [3, 4]]
The fields of the matrix are:
matrix.h: 2 (number of rows or "height")
matrix.w: 2 (number of columns or "width")
matrix.e: 2D array of 2x2 elements, first index is row, second is column

*) Create a complex matrix [[1+2i, 3+4i], [5+6i, 7+8i]]
var matrix = new Matrix(2, 2);
matrix.e[0][0] = new Complex(1, 2);
matrix.e[0][1] = new Complex(3, 4);
matrix.e[1][0] = new Complex(5, 6);
matrix.e[1][1] = new Complex(7, 8);
shorthands for all the above:
var matrix = Matrix([[Complex(1, 2), Complex(3, 4)], [Complex(5, 6), Complex(7, 8)]])
var matrix = Matrix([[1, 3], [5, 7]], [[2, 4], [6, 8]])
var matrix = Matrix(2, 2, Complex(1, 2), Complex(3, 4), Complex(5, 6), Complex(7, 8))
var matrix = Matrix('[[1+2i, 3+4i], [5+6i, 7+8i]]')

*) Calculating with matrices
Matrix([[1, 2], [3, 4]]).mul(Matrix([[5, 6], [7, 8]])).toString()
--> [[19, 22], [43, 50]]
Matrix([[1, 2], [3, 4]]).add(Matrix([[5, 6], [7, 8]]).mulc(Complex.I)).toString()
--> [[1+5i, 2+6i], [3+7i, 4+8i]]

*) Loop through the elements of a matrix
var a = Matrix([[1, 2], [3, 4]])
for(var y = 0; y < a.h; y++) { // row index
  for(var x = 0; x < a.w; x++) { // column index
    var complex = a.e[y][x]; // complex is of type Complex
    console.log('x: ' + x + ' y: ' + y + ' re: ' + complex.re + ' im: ' + complex.im);
  }
}

*) Inverse of a matrix
Matrix.inv(Matrix([[1,2],[3,4]])).toString()
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
Jmat.integrate(Complex(0), Complex(10), function(x) { return x.mul(x); }, 30).toString();
--> 333.3333333333333

*) Real contains functions which operate on regular JS numbers, so could be seen as an extension to the standard JS Math library:
Real.gamma(Math.cos(2)) + Real.erfc(0.5) * Real.EM
--> -3.39388952436638
Real.nextPrime(17)
--> 19

*) Numerical integration (quadrature) of x^2 from 0 to 10 with 30 steps, with real numbers (no Complex objects):
Real.integrate(0, 10, function(x) { return x * x; }, 30);
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

*) Plot a Mandelbrot set:
Jmat.plotComplex(function(c) {
  var i = 0; var z = Complex(0); for(;;) { if(z.absr() > 2) break; z = z.mul(z).add(c); i++; if(i > 60) return Complex(0); } return Complex.polar(1, i * Math.PI / 60);
}, plotContainerEl, {p:1, s:4});


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

Feel free to comment about bugs, improvements, experiences, ... as well as to
contribute.

GitHub page: https://github.com/lvandeve/jmat

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
  // Empty, this is a namespace, no need to ever call this
}

/*
The Jmat functions below are the most public API (as well as any Jmat.Real, Jmat.Complex or
Jmat.Matrix functions with the same name).

Unlike their Jmat.Complex and Jmat.Matrix counterparts, these versions are more tolerant to
various input types, e.g. Jmat.gamma can take both 5 and Jmat.Complex(5) as input,
while Jmat.Complex.gamma only accepts Jmat.Complex(5) and Jmat.Real.gamma only accepts 5.

To work with complex numbers or matrices, you'll need to use Jmat.Complex and Jmat.Matrix
objects anyway, create them with Jmat.Complex or Jmat.Matrix.

Return values should always be treated as immutable objects. Do not modify the
re and im fields directly: doing so could result in changing the internal
constants (like Jmat.Complex.PI).

The Jmat.Complex and Jmat.Matrix objects contain several functions in their prototype as
well, e.g. add, so it is possible to write Jmat.Complex(5).add(Jmat.Complex(6)) instead of
Jmat.add(Jmat.Complex(5), Jmat.Complex(6)). However, those prototype functions do not accept
input of a wrong type (must be Jmat.Complex, not a plain JS number).

The comments below use closure-style types, e.g. {number|Complex} means the type
is either a plain JS number of a Jmat.Complex object, and {Array.<Complex>} means an array
of Jmat.Complex objects. Jmat.Complex is aliased as Complex, Jmat.Matrix is aliased as Matrix.
The shorter aliases are used in these comments.
*/

// Elementary operators

/* Add. x,y:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.add = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.add(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  return Jmat.Complex.add(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Subtract. x,y:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.sub = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.sub(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  return Jmat.Complex.sub(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Multiply. x,y:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.mul = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.mul(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  if(Jmat.matrixIn_(x)) return Jmat.Matrix.mulc(Jmat.Matrix.cast(x), Jmat.Complex.cast(y));
  if(Jmat.matrixIn_(y)) return Jmat.Matrix.mulc(Jmat.Matrix.cast(y), Jmat.Complex.cast(x));
  return Jmat.Complex.mul(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Division. x,y:{number|Jmat.Complex|Jmat.Matrix}. returns {Jmat.Complex|Jmat.Matrix}. */
Jmat.div = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.div(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  if(Jmat.matrixIn_(x)) return Jmat.Matrix.divc(Jmat.Matrix.cast(x), Jmat.Complex.cast(y));
  if(Jmat.matrixIn_(y)) return Jmat.Matrix.divc(Jmat.Matrix.cast(y), Jmat.Complex.cast(x));
  return Jmat.Complex.div(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};

// Compare

/* Equal? x,y:{number|Complex|Matrix}. returns {boolean}. */
Jmat.eq = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.eq(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  return Jmat.Complex.eq(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Nearly equal? x,y:{number|Complex|Matrix}. epsilon:{number}. returns {boolean}. */
Jmat.near = function(x, y, epsilon) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.near(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y), Jmat.Real.caststrict(epsilon));
  return Jmat.Complex.near(Jmat.Complex.cast(x), Jmat.Complex.cast(y), Jmat.Real.caststrict(epsilon));
};

// Categories

Jmat.isNaN = function(x) { return Jmat.matrixIn_(x) ? Jmat.Matrix.isNaN(Jmat.Matrix.cast(x)) : Jmat.Complex.isNaN(Jmat.Complex.cast(x)); };

// Power & Logarithms

/* Power, or matrix power. x:{number|Complex|Matrix}. y:{number|Complex}. returns {Complex}. */
Jmat.pow = function(x, y) {
  if(Jmat.matrixIn_(x)) return Jmat.Matrix.powc(Jmat.Matrix.cast(x), Jmat.Complex.cast(y)); // matrix raised to any complex number power
  return Jmat.Complex.pow(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Square root, or matrix square root. z:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.sqrt = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.sqrt(Jmat.Matrix.cast(z));
  return Jmat.Complex.sqrt(Jmat.Complex.cast(z));
};
/* Exponential, or matrix exponential. z:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.exp = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.exp(Jmat.Matrix.cast(z)); // matrix exponential
  return Jmat.Complex.exp(Jmat.Complex.cast(z));
};
/* Logarithm, or matrix logarithm. z:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.log = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.log(Jmat.Matrix.cast(z));
  return Jmat.Complex.log(Jmat.Complex.cast(z));
};
/* log(z) + 1. z:{number|Complex}. returns {Complex}. */
Jmat.log1p = function(z) { return Jmat.Complex.log1p(Jmat.Complex.cast(z)); };
/* exp(z) - 1. z:{number|Complex}. returns {Complex}. */
Jmat.expm1 = function(z) { return Jmat.Complex.expm1(Jmat.Complex.cast(z)); };
/* Base-y logarithm of x. x,y:{number|Complex}. returns {Complex}. */
Jmat.logy = function(x, y) { return Jmat.Complex.logy(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };
/* Principal branch of LambertW. z:{number|Complex}. returns {Complex}. */
Jmat.lambertw = function(z) { return Jmat.Complex.lambertw(Jmat.Complex.cast(z)); };
/* Negative branch of LambertW. z:{number|Complex}. returns {Complex}. */
Jmat.lambertwm = function(z) { return Jmat.Complex.lambertwm(Jmat.Complex.cast(z)); };
/* Specific branch of LambertW. branch:{number} must be integer, z:{number|Complex}. returns {Complex}. */
Jmat.lambertwb = function(branch, z) { return Jmat.Complex.lambertwb(Jmat.Real.caststrict(branch), Jmat.Complex.cast(z)); };
/* Tetration (power tower). x,y:{number|Complex}. returns {Complex} */
Jmat.tetration = function(x, y) { return Jmat.Complex.tetration(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };

// Elementary functions

/* Identity function. x:{number|Complex}. returns {Complex} */
Jmat.x = function(x) { return Jmat.Complex.cast(x); };
/* Negate. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.neg = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.neg(Jmat.Matrix.cast(z));
  return Jmat.Complex.neg(Jmat.Complex.cast(z));
};
/* Reciproke. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.inv = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.inv(Jmat.Matrix.cast(z));
  return Jmat.Complex.inv(Jmat.Complex.cast(z));
};
/* Complex conjugate. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.conj = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.conj(Jmat.Matrix.cast(z));
  return Jmat.Complex.conj(Jmat.Complex.cast(z));
};

// Trigonometric functions

/* Sine, or matrix-sine. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.sin = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.sin(Jmat.Matrix.cast(z));
  return Jmat.Complex.sin(Jmat.Complex.cast(z));
};
/* Cosine, or matrix-cosine. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.cos = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.cos(Jmat.Matrix.cast(z));
  return Jmat.Complex.cos(Jmat.Complex.cast(z));
};
/* Tangent. z:{number|Complex}. returns {Complex} */
Jmat.tan = function(z) { return Jmat.Complex.tan(Jmat.Complex.cast(z)); };
/* Arcsine. z:{number|Complex}. returns {Complex} */
Jmat.asin = function(z) { return Jmat.Complex.asin(Jmat.Complex.cast(z)); };
/* Arccosine. z:{number|Complex}. returns {Complex} */
Jmat.acos = function(z) { return Jmat.Complex.acos(Jmat.Complex.cast(z)); };
/* Arctangent. z:{number|Complex}. returns {Complex} */
Jmat.atan = function(z) { return Jmat.Complex.atan(Jmat.Complex.cast(z)); };
/* Atan2 function. z:{number|Complex}. returns {Complex} */
Jmat.atan2 = function(x, y) { return Jmat.Complex.atan2(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };
/* Unnormalized sinc: sin(z)/z. z:{number|Complex}. returns {Complex} */
Jmat.sinc = function(z) { return Jmat.Complex.sinc(Jmat.Complex.cast(z)); };

// Hyperbolic functions

/* Hyperbolic sine. z:{number|Complex}. returns {Complex} */
Jmat.sinh = function(z) { return Jmat.Complex.sinh(Jmat.Complex.cast(z)); };
/* Hyperbolic cosine. z:{number|Complex}. returns {Complex} */
Jmat.cosh = function(z) { return Jmat.Complex.cosh(Jmat.Complex.cast(z)); };
/* Hyperbolic tangent. z:{number|Complex}. returns {Complex} */
Jmat.tanh = function(z) { return Jmat.Complex.tanh(Jmat.Complex.cast(z)); };
/* Hyperbolic arcsine. z:{number|Complex}. returns {Complex} */
Jmat.asinh = function(z) { return Jmat.Complex.asinh(Jmat.Complex.cast(z)); };
/* Hyperbolic arccosine. z:{number|Complex}. returns {Complex} */
Jmat.acosh = function(z) { return Jmat.Complex.acosh(Jmat.Complex.cast(z)); };
/* Hyperbolic arctangent. z:{number|Complex}. returns {Complex} */
Jmat.atanh = function(z) { return Jmat.Complex.atanh(Jmat.Complex.cast(z)); };

// Gamma and related functions

/* Gamma function. z:{number|Complex}. returns {Complex} */
Jmat.gamma = function(z) { return Jmat.Complex.gamma(Jmat.Complex.cast(z)); };
/* Factorial. z:{number|Complex}. returns {Complex} */
Jmat.factorial = function(z) { return Jmat.Complex.factorial(Jmat.Complex.cast(z)); };
/* Digamma function (psi). z:{number|Complex}. returns {Complex} */
Jmat.digamma = function(z) { return Jmat.Complex.digamma(Jmat.Complex.cast(z)); };
/* Trigamma function. z:{number|Complex}. returns {Complex} */
Jmat.trigamma = function(z) { return Jmat.Complex.trigamma(Jmat.Complex.cast(z)); };
/* Polygamma function. n,z:{number|Complex}. returns {Complex} */
Jmat.polygamma = function(n, z) { return Jmat.Complex.polygamma(Jmat.Complex.cast(n), Jmat.Complex.cast(z)); };
/* Logarithm of gamma function. z:{number|Complex}. returns {Complex} */
Jmat.loggamma = function(z) { return Jmat.Complex.loggamma(Jmat.Complex.cast(z)); };
/* Inverse of gamma function (not reciproke). z:{number|Complex}. returns {Complex} */
Jmat.gamma_inv = function(z) { return Jmat.Complex.gamma_inv(Jmat.Complex.cast(z)); };
/* Lower incomplete gamma function. s,z:{number|Complex}. returns {Complex} */
Jmat.incgamma_lower = function(s, z) { return Jmat.Complex.incgamma_lower(Jmat.Complex.cast(s), Jmat.Complex.cast(z)); };
/* Upper incomplete gamma function. s,z:{number|Complex}. returns {Complex} */
Jmat.incgamma_upper = function(s, z) { return Jmat.Complex.incgamma_upper(Jmat.Complex.cast(s), Jmat.Complex.cast(z)); };
/* Lower regularized incomplete gamma function "P". s,z:{number|Complex}. returns {Complex} */
Jmat.gamma_p = function(s, z) { return Jmat.Complex.gamma_p(Jmat.Complex.cast(s), Jmat.Complex.cast(z)); };
/* Upper regularized incomplete gamma function. s,z:{number|Complex}. returns {Complex} */
Jmat.gamma_q = function(s, z) { return Jmat.Complex.gamma_q(Jmat.Complex.cast(s), Jmat.Complex.cast(z)); };
/* Inverse of "P" in z. s,o:{number|Complex}. returns {Complex} */
Jmat.gamma_p_inv = function(s, p) { return Jmat.Complex.gamma_p_inv(Jmat.Complex.cast(s), Jmat.Complex.cast(p)); };
/* Beta function. x,y:{number|Complex}. returns {Complex} */
Jmat.beta = function(x, y) { return Jmat.Complex.beta(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };
/* Incomplete beta function. x,a,b:{number|Complex}. returns {Complex} */
Jmat.incbeta = function(x, a, b) { return Jmat.Complex.incbeta(Jmat.Complex.cast(x), Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };
/* Regularized incomplete beta function "I". x,a,b:{number|Complex}. returns {Complex} */
Jmat.beta_i = function(x, a, b) { return Jmat.Complex.beta_i(Jmat.Complex.cast(x), Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };
/* Inverse of "I" in x. x,a,b:{number|Complex}. returns {Complex} */
Jmat.beta_i_inv = function(x, a, b) { return Jmat.Complex.beta_i_inv(Jmat.Complex.cast(x), Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };

// Number and complex number utility functions

/* Real part. z:{number|Complex}. returns {Complex} */
Jmat.re = function(z) { return Jmat.Complex(Jmat.Complex.cast(z).re); };
/* Imaginary part. z:{number|Complex}. returns {Complex} */
Jmat.im = function(z) { return Jmat.Complex(Jmat.Complex.cast(z).im); };
/* Absolute value or complex modulus. z:{number|Complex}. returns {Complex} */
Jmat.abs = function(z) { return Jmat.Complex.abs(Jmat.Complex.cast(z)); }; // absolute value, modulus
/* Complex argument or phase. z:{number|Complex}. returns {Complex} */
Jmat.arg = function(z) { return Jmat.Complex.arg(Jmat.Complex.cast(z)); };
/* Real sign. z:{number|Complex}. returns {Complex} */
Jmat.sign = function(z) { return Jmat.Complex.sign(Jmat.Complex.cast(z)); };
/* Complex sign. z:{number|Complex}. returns {Complex} */
Jmat.csgn = function(z) { return Jmat.Complex.csgn(Jmat.Complex.cast(z)); };
/* Floor. x:{number|Complex}. returns {Complex} */
Jmat.floor = function(x) { return Jmat.Complex.floor(Jmat.Complex.cast(x)); };
/* Ceiling. x:{number|Complex}. returns {Complex} */
Jmat.ceil = function(x) { return Jmat.Complex.ceil(Jmat.Complex.cast(x)); };
/* Round to integer. x:{number|Complex}. returns {Complex} */
Jmat.round = function(x) { return Jmat.Complex.round(Jmat.Complex.cast(x)); };
/* Truncate towards zero. x:{number|Complex}. returns {Complex} */
Jmat.trunc = function(x) { return Jmat.Complex.trunc(Jmat.Complex.cast(x)); };
/* Fractional part, always positive. x:{number|Complex}. returns {Complex} */
Jmat.frac = function(x) { return Jmat.Complex.frac(Jmat.Complex.cast(x)); };
/* Fractional part, negative for negative x. x:{number|Complex}. returns {Complex} */
Jmat.fracn = function(x) { return Jmat.Complex.fracn(Jmat.Complex.cast(x)); };
/* Rotate complex number by a radians. z:{number|Complex}, a:{number}. returns {Complex} */
Jmat.rotate = function(z, a) { return Jmat.Complex.rotate(Jmat.Complex.cast(z), Jmat.Real.caststrict(a)); };
/* Approximate with integer numerator and denominator. x:{number|Complex}, max:{number} maximum denominator. returns {Array.<Complex>} numerator, denominator */
Jmat.decompose = function(x, max) { return Jmat.Complex.decompose(Jmat.Complex.cast(x), Jmat.Real.cast(max)); };

// Cylindrical functions

/* Bessel function of the first kind.  n,z:{number|Complex}. returns {Complex} */
Jmat.besselj = function(n, z) { return Jmat.Complex.besselj(Jmat.Complex.cast(n), Jmat.Complex.cast(z)); };
/* Bessel function of the second kind.  n,z:{number|Complex}. returns {Complex} */
Jmat.bessely = function(n, z) { return Jmat.Complex.bessely(Jmat.Complex.cast(n), Jmat.Complex.cast(z)); };
/* Modified bessel function of the first kind.  n,z:{number|Complex}. returns {Complex} */
Jmat.besseli = function(n, z) { return Jmat.Complex.besseli(Jmat.Complex.cast(n), Jmat.Complex.cast(z)); };
/* Modified bessel function of the second kind.  n,z:{number|Complex}. returns {Complex} */
Jmat.besselk = function(n, z) { return Jmat.Complex.besselk(Jmat.Complex.cast(n), Jmat.Complex.cast(z)); };
/* Hankel function of the first kind.  n,z:{number|Complex}. returns {Complex} */
Jmat.hankel1 = function(n, z) { return Jmat.Complex.hankel1(Jmat.Complex.cast(n), Jmat.Complex.cast(z)); };
/* Hankel function of the second kind.  n,z:{number|Complex}. returns {Complex} */
Jmat.hankel2 = function(n, z) { return Jmat.Complex.hankel2(Jmat.Complex.cast(n), Jmat.Complex.cast(z)); };
/* Airy function Ai. z:{number|Complex}. returns {Complex} */
Jmat.airy = function(z) { return Jmat.Complex.airy(Jmat.Complex.cast(z)); };
/* Bairy function Bi. z:{number|Complex}. returns {Complex} */
Jmat.bairy = function(z) { return Jmat.Complex.bairy(Jmat.Complex.cast(z)); };
/* Derivative of Airy function Ai. z:{number|Complex}. returns {Complex} */
Jmat.airy_deriv = function(z) { return Jmat.Complex.airy_deriv(Jmat.Complex.cast(z)); };
/* Derivative of Bairy function Bi. z:{number|Complex}. returns {Complex} */
Jmat.bairy_deriv = function(z) { return Jmat.Complex.bairy_deriv(Jmat.Complex.cast(z)); };

// Theta functions

/* Jacobi theta function θ1. z,q:{number|Complex}. returns {Complex} */
Jmat.theta1 = function(z, q) { return Jmat.Complex.theta1(Jmat.Complex.cast(z), Jmat.Complex.cast(q)); };
/* Jacobi theta function θ2. z,q:{number|Complex}. returns {Complex} */
Jmat.theta2 = function(z, q) { return Jmat.Complex.theta2(Jmat.Complex.cast(z), Jmat.Complex.cast(q)); };
/* Jacobi theta function θ3. z,q:{number|Complex}. returns {Complex} */
Jmat.theta3 = function(z, q) { return Jmat.Complex.theta3(Jmat.Complex.cast(z), Jmat.Complex.cast(q)); };
/* Jacobi theta function θ4. z,q:{number|Complex}. returns {Complex} */
Jmat.theta4 = function(z, q) { return Jmat.Complex.theta4(Jmat.Complex.cast(z), Jmat.Complex.cast(q)); };

// Zeta and related functions

/* Riemann zeta function. z:{number|Complex}. returns {Complex} */
Jmat.zeta = function(z) { return Jmat.Complex.zeta(Jmat.Complex.cast(z)); };
/* Dirichlet eta function. z:{number|Complex}. returns {Complex} */
Jmat.eta = function(z) { return Jmat.Complex.eta(Jmat.Complex.cast(z)); };
/* Dirichlet lambda function. z:{number|Complex}. returns {Complex} */
Jmat.lambda = function(z) { return Jmat.Complex.lambda(Jmat.Complex.cast(z)); };
/* Dilogarithm. z:{number|Complex}. returns {Complex} */
Jmat.dilog = function(z) { return Jmat.Complex.dilog(Jmat.Complex.cast(z)); };
/* Trilogarithm. z:{number|Complex}. returns {Complex} */
Jmat.trilog = function(z) { return Jmat.Complex.trilog(Jmat.Complex.cast(z)); };
/* Polylogarithm. s,z:{number|Complex}. returns {Complex} */
Jmat.polylog = function(s, z) { return Jmat.Complex.polylog(Jmat.Complex.cast(s), Jmat.Complex.cast(z)); };
/* Hurwitz Zeta function. s,q:{number|Complex}. returns {Complex} */
Jmat.hurwitzzeta = function(s, q) { return Jmat.Complex.hurwitzzeta(Jmat.Complex.cast(s), Jmat.Complex.cast(q)); };

// Hypergeometric functions

/* Hypergeometric function 0F1. a,z:{number|Complex}. returns {Complex} */
Jmat.hypergeometric0F1 = function(a, z) { return Jmat.Complex.hypergeometric0F1(Jmat.Complex.cast(a), Jmat.Complex.cast(z)); };
/* Confluent Hypergeometric function (Kummer). a,b,z:{number|Complex}. returns {Complex} */
Jmat.hypergeometric1F1 = function(a, b, z) { return Jmat.Complex.hypergeometric1F1(Jmat.Complex.cast(a), Jmat.Complex.cast(b), Jmat.Complex.cast(z)); };
/* Hypergeometric function 2F1. a,b,c,z:{number|Complex}. returns {Complex} */
Jmat.hypergeometric = function(a, b, c, z) { return Jmat.Complex.hypergeometric(Jmat.Complex.cast(a), Jmat.Complex.cast(b), Jmat.Complex.cast(c), Jmat.Complex.cast(z)); };

// Error and related functions

/* Error function. z:{number|Complex}. returns {Complex} */
Jmat.erf = function(z) { return Jmat.Complex.erf(Jmat.Complex.cast(z)); };
/* Complementary error function. z:{number|Complex}. returns {Complex} */
Jmat.erfc = function(z) { return Jmat.Complex.erfc(Jmat.Complex.cast(z)); };
/* Scaled complementary error function. z:{number|Complex}. returns {Complex} */
Jmat.erfcx = function(z) { return Jmat.Complex.erfcx(Jmat.Complex.cast(z)); };
/* Inverse of error function (not reciproke). z:{number|Complex}. returns {Complex} */
Jmat.erf_inv = function(z) { return Jmat.Complex.erf_inv(Jmat.Complex.cast(z)); };
/* Inverse of complementary error function (not reciproke). z:{number|Complex}. returns {Complex} */
Jmat.erfc_inv = function(z) { return Jmat.Complex.erfc_inv(Jmat.Complex.cast(z)); };
/* Imaginary error function. x:{number|Complex}. returns {Complex} */
Jmat.erfi = function(x) { return x.re == undefined ? Jmat.Complex(Jmat.Real.erfi(x)) : Jmat.Complex.erfi(Jmat.Complex.cast(x)); };
/* Dawson function D+(x). x:{number|Complex}. returns {Complex} */
Jmat.dawson = function(x) { return x.re == undefined ? Jmat.Complex(Jmat.Real.dawson(x)) : Jmat.Complex.dawson(Jmat.Complex.cast(x)); };
/* Faddeeva function w(z). z:{number|Complex}. returns {Complex} */
Jmat.faddeeva = function(z) { return Jmat.Complex.faddeeva(Jmat.Complex.cast(z)); };

// Statistical distributions. pdf = probability density function, cdf = cumulative distribution function, qf = quantile function

/* Uniform distribution in range a-b. x,a,b:{number|Complex}. returns {Complex} */
Jmat.pdf_uniform = function(x, a, b) { return Jmat.Complex.pdf_uniform(Jmat.Complex.cast(x), Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };
Jmat.cdf_uniform = function(x, a, b) { return Jmat.Complex.cdf_uniform(Jmat.Complex.cast(x), Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };
Jmat.qf_uniform = function(x, a, b) { return Jmat.Complex.qf_uniform(Jmat.Complex.cast(x), Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };
/* Standardnormal distribution. x:{number|Complex}. returns {Complex} */
Jmat.pdf_standardnormal = function(x) { return Jmat.Complex.pdf_standardnormal(Jmat.Complex.cast(x)); };
Jmat.cdf_standardnormal = function(x) { return Jmat.Complex.cdf_standardnormal(Jmat.Complex.cast(x)); };
Jmat.qf_standardnormal = function(x) { return Jmat.Complex.qf_standardnormal(Jmat.Complex.cast(x)); };
/* Normal distribution with mean mu and variance sigma. x,mu,sigma:{number|Complex}. returns {Complex} */
Jmat.pdf_normal = function(x, mu, sigma) { return Jmat.Complex.pdf_normal(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(sigma)); };
Jmat.cdf_normal = function(x, mu, sigma) { return Jmat.Complex.cdf_normal(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(sigma)); };
Jmat.qf_normal = function(x, mu, sigma) { return Jmat.Complex.qf_normal(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(sigma)); };
/* Log-normal distribution with mean mu and variance sigma. x,mu,sigma:{number|Complex}. returns {Complex} */
Jmat.pdf_lognormal = function(x, mu, sigma) { return Jmat.Complex.pdf_lognormal(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(sigma)); };
Jmat.cdf_lognormal = function(x, mu, sigma) { return Jmat.Complex.cdf_lognormal(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(sigma)); };
Jmat.qf_lognormal = function(x, mu, sigma) { return Jmat.Complex.qf_lognormal(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(sigma)); };
/* Cauchy distribution with location x0 and scale gamma. x,x0,gamma:{number|Complex}. returns {Complex} */
Jmat.pdf_cauchy = function(x, x0, gamma) { return Jmat.Complex.pdf_cauchy(Jmat.Complex.cast(x), Jmat.Complex.cast(x0), Jmat.Complex.cast(gamma)); };
Jmat.cdf_cauchy = function(x, x0, gamma) { return Jmat.Complex.cdf_cauchy(Jmat.Complex.cast(x), Jmat.Complex.cast(x0), Jmat.Complex.cast(gamma)); };
Jmat.qf_cauchy = function(x, x0, gamma) { return Jmat.Complex.qf_cauchy(Jmat.Complex.cast(x), Jmat.Complex.cast(x0), Jmat.Complex.cast(gamma)); };
/* Student's t distribution with degrees of freedom nu. x,nu:{number|Complex}. returns {Complex} */
Jmat.pdf_studentt = function(x, nu) { return Jmat.Complex.pdf_studentt(Jmat.Complex.cast(x), Jmat.Complex.cast(nu)); };
Jmat.cdf_studentt = function(x, nu) { return Jmat.Complex.cdf_studentt(Jmat.Complex.cast(x), Jmat.Complex.cast(nu)); };
Jmat.qf_studentt = function(x, nu) { return Jmat.Complex.qf_studentt(Jmat.Complex.cast(x), Jmat.Complex.cast(nu)); };
/* Chi square distribution with degrees of freedom k. x,k:{number|Complex}. returns {Complex} */
Jmat.pdf_chi_square = function(x, k) { return Jmat.Complex.pdf_chi_square(Jmat.Complex.cast(x), Jmat.Complex.cast(k)); };
Jmat.cdf_chi_square = function(x, k) { return Jmat.Complex.cdf_chi_square(Jmat.Complex.cast(x), Jmat.Complex.cast(k)); };
Jmat.qf_chi_square = function(x, k) { return Jmat.Complex.qf_chi_square(Jmat.Complex.cast(x), Jmat.Complex.cast(k)); };
/* Logistic distribution with location mu and scale s. x,mu,s:{number|Complex}. returns {Complex} */
Jmat.pdf_logistic = function(x, mu, s) { return Jmat.Complex.pdf_logistic(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(s)); };
Jmat.cdf_logistic = function(x, mu, s) { return Jmat.Complex.cdf_logistic(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(s)); };
Jmat.qf_logistic = function(x, mu, s) { return Jmat.Complex.qf_logistic(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(s)); };
/* Gamma distribution with shape k and scale theta. x,k,theta:{number|Complex}. returns {Complex} */
Jmat.pdf_gamma = function(x, k, theta) { return Jmat.Complex.pdf_gamma(Jmat.Complex.cast(x), Jmat.Complex.cast(k), Jmat.Complex.cast(theta)); };
Jmat.cdf_gamma = function(x, k, theta) { return Jmat.Complex.cdf_gamma(Jmat.Complex.cast(x), Jmat.Complex.cast(k), Jmat.Complex.cast(theta)); };
Jmat.qf_gamma = function(x, k, theta) { return Jmat.Complex.qf_gamma(Jmat.Complex.cast(x), Jmat.Complex.cast(k), Jmat.Complex.cast(theta)); };
/* Beta distribution with shape alpha and beta. x,alpha,beta:{number|Complex}. returns {Complex} */
Jmat.pdf_beta = function(x, alpha, beta) { return Jmat.Complex.pdf_beta(Jmat.Complex.cast(x), Jmat.Complex.cast(alpha), Jmat.Complex.cast(beta)); };
Jmat.cdf_beta = function(x, alpha, beta) { return Jmat.Complex.cdf_beta(Jmat.Complex.cast(x), Jmat.Complex.cast(alpha), Jmat.Complex.cast(beta)); };
Jmat.qf_beta = function(x, alpha, beta) { return Jmat.Complex.qf_beta(Jmat.Complex.cast(x), Jmat.Complex.cast(alpha), Jmat.Complex.cast(beta)); };
/* F-distribution with degrees of freedom d1 and d2. x,d1,d2:{number|Complex}. returns {Complex} */
Jmat.pdf_fisher = function(x, d1, d2) { return Jmat.Complex.pdf_fisher(Jmat.Complex.cast(x), Jmat.Complex.cast(d1), Jmat.Complex.cast(d2)); };
Jmat.cdf_fisher = function(x, d1, d2) { return Jmat.Complex.cdf_fisher(Jmat.Complex.cast(x), Jmat.Complex.cast(d1), Jmat.Complex.cast(d2)); };
Jmat.qf_fisher = function(x, d1, d2) { return Jmat.Complex.qf_fisher(Jmat.Complex.cast(x), Jmat.Complex.cast(d1), Jmat.Complex.cast(d2)); };
/* Weibull distribution with scale lambda, shape k. x,lambda,k:{number|Complex}. returns {Complex} */
Jmat.pdf_weibull = function(x, lambda, k) { return Jmat.Complex.pdf_weibull(Jmat.Complex.cast(x), Jmat.Complex.cast(lambda), Jmat.Complex.cast(k)); };
Jmat.cdf_weibull = function(x, lambda, k) { return Jmat.Complex.cdf_weibull(Jmat.Complex.cast(x), Jmat.Complex.cast(lambda), Jmat.Complex.cast(k)); };
Jmat.qf_weibull = function(x, lambda, k) { return Jmat.Complex.qf_weibull(Jmat.Complex.cast(x), Jmat.Complex.cast(lambda), Jmat.Complex.cast(k)); };
/* Exponential distribution with rate lambda. x,lambda:{number|Complex}. returns {Complex} */
Jmat.pdf_exponential = function(x, lambda) { return Jmat.Complex.pdf_exponential(Jmat.Complex.cast(x), Jmat.Complex.cast(lambda)); };
Jmat.cdf_exponential = function(x, lambda) { return Jmat.Complex.cdf_exponential(Jmat.Complex.cast(x), Jmat.Complex.cast(lambda)); };
Jmat.qf_exponential = function(x, lambda) { return Jmat.Complex.qf_exponential(Jmat.Complex.cast(x), Jmat.Complex.cast(lambda)); };
/* Laplace distribution with location mu, scale b. x,mu,b:{number|Complex}. returns {Complex} */
Jmat.pdf_laplace = function(x, mu, b) { return Jmat.Complex.pdf_laplace(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(b)); };
Jmat.cdf_laplace = function(x, mu, b) { return Jmat.Complex.cdf_laplace(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(b)); };
Jmat.qf_laplace = function(x, mu, b) { return Jmat.Complex.qf_laplace(Jmat.Complex.cast(x), Jmat.Complex.cast(mu), Jmat.Complex.cast(b)); };

// Combinatorics

/* Number of permutations. n,p:{number|Complex}. returns {Complex} */
Jmat.permutation = function(n, p) { return Jmat.Complex.permutation(Jmat.Complex.cast(n), Jmat.Complex.cast(p)); };
/* Binomial, number of combinations. n,p:{number|Complex}. returns {Complex} */
Jmat.binomial = function(n, p) { return Jmat.Complex.binomial(Jmat.Complex.cast(n), Jmat.Complex.cast(p)); };
/* Stirling number of hte second kind. n{number|Complex} integer, p:{number|Complex}. returns {Complex} */
Jmat.stirling2 = function(n, k) { return Jmat.Complex.stirling2(Jmat.Complex.cast(n), Jmat.Complex.cast(k)); };

// Mean

/* Arithmetic-Geometric mean. a,b:{number|Complex}. returns {Complex} */
Jmat.agm = function(a, b) { return Jmat.Complex.agm(Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };
/* Geometric-Harmonic mean. a,b:{number|Complex}. returns {Complex} */
Jmat.ghm = function(a, b) { return Jmat.Complex.ghm(Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };

// Prime numbers. NOTE: Use the Jmat.Real versions to get results as simple JS numbers rather than Jmat.Complex objects.

/* Primality test. x:{number|Complex}. returns {boolean} */
Jmat.isPrime = function(x) { return Jmat.Complex.isPrime(Jmat.Complex.cast(x)); };
/* Smallest prime factor of integer. x:{number|Complex} real. returns {Complex} */
Jmat.smallestPrimeFactor = function(x) { return Jmat.Complex(Jmat.Real.smallestPrimeFactor(Jmat.Real.cast(x))); };
/* Factorize into prime factors, returned as array. x:{number|Complex} integer. returns {Array.<Complex>} */
Jmat.factorize = function(n) {
  var res = Jmat.Real.factorize(Jmat.Real.caststrict(n));
  for(var i = 0; i < res.length; i++) res[i] = Jmat.Complex(res[i]); // Convert to Jmat.Complex objects to be consistent with other Jmat functions. Just use Jmat.Real.factorize instead to have simple JS numbers.
  return res;
};
/* Euler's Totient. x:{number|Complex} real integer. returns {Complex} */
Jmat.eulerTotient = function(x) { return Jmat.Complex(Jmat.Real.eulerTotient(Jmat.Real.caststrict(x))); };
/* Prime counting function. x:{number|Complex} real. returns {Complex} */
Jmat.primeCount = function(x) { return Jmat.Complex(Jmat.Real.primeCount(Jmat.Real.caststrict(x))); };
/* Nearest prime number. x:{number|Complex}. returns {Complex} */
Jmat.nearestPrime = function(x) { return Jmat.Complex(Jmat.Real.nearestPrime(Jmat.Real.cast(x))); };
/* Next larger prime number function. x:{number|Complex}. returns {Complex} */
Jmat.nextPrime = function(x) { return Jmat.Complex(Jmat.Real.nextPrime(Jmat.Real.cast(x))); };
/* Next smaller prime number function. x:{number|Complex}. returns {Complex} */
Jmat.previousPrime = function(x) { return Jmat.Complex(Jmat.Real.previousPrime(Jmat.Real.cast(x))); };
/* Greatest common divisor. x,y:{number|Complex} real integer. returns {Complex} */
Jmat.gcd = function(x, y) { return Jmat.Complex(Jmat.Real.gcd(Jmat.Real.caststrict(x), Jmat.Real.caststrict(y))); };
/* Least common multiple. x,y:{number|Complex} real integer. returns {Complex} */
Jmat.lcm = function(x, y) { return Jmat.Complex(Jmat.Real.lcm(Jmat.Real.caststrict(x), Jmat.Real.caststrict(y))); };

// Sexagesimal

/* Decimal degrees to degrees, minutes, seconds. x:{number|Complex} real. returns {Complex} */
Jmat.dms = function(a) { return Jmat.Complex(Jmat.Real.dms(Jmat.Real.caststrict(a))); }; // E.g. 1.6 becomes 1.36 (1 degree, 36 minutes)
/* Degrees, minutes, seconds to decimal degrees. x:{number|Complex} real. returns {Complex} */
Jmat.dd = function(a) { return Jmat.Complex(Jmat.Real.dd(Jmat.Real.caststrict(a))); }; // E.g. 1.36 becomes 1.6 (1 point 6 degrees)

// Bitwise

/* Bitwise not. x:{number|Complex}. returns {Complex} */
Jmat.bitnot = function(x) { return Jmat.Complex.bitnot(Jmat.Complex.cast(x)); };
/* Bitwise and. x,y:{number|Complex}. returns {Complex} */
Jmat.bitand = function(x, y) { return Jmat.Complex.bitand(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };
/* Bitwise or. x,y:{number|Complex}. returns {Complex} */
Jmat.bitor = function(x, y) { return Jmat.Complex.bitor(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };
/* Bitwise xor. x,y:{number|Complex}. returns {Complex} */
Jmat.bitxor = function(x, y) { return Jmat.Complex.bitxor(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };
/* Left bitshift. x:{number|Complex}, y:{number|Complex} real. returns {Complex} */
Jmat.lshift = function(x, y) { return Jmat.Complex.lshift(Jmat.Complex.cast(x), Jmat.Real.caststrict(y)); };
/* Right bitshift. x:{number|Complex}, y:{number|Complex} real. returns {Complex} */
Jmat.rshift = function(x, y) { return Jmat.Complex.rshift(Jmat.Complex.cast(x), Jmat.Real.caststrict(y)); };

// Modulo division and related functions

/* Modulo division. a,b:{number|Complex}. returns {Complex} */
Jmat.mod = function(a, b) { return Jmat.Complex.mod(Jmat.Complex.cast(a), Jmat.Complex.cast(b)); }; // result has sign of divisor (unlike JS '%' operator)
/* Remainder. a,b:{number|Complex}. returns {Complex} */
Jmat.rem = function(a, b) { return Jmat.Complex.rem(Jmat.Complex.cast(a), Jmat.Complex.cast(b)); }; // result has sign of dividend (same result as JS '%' operator on real numbers)
/* Wrap x between from and to (to excluded). x,to,from:{number|Complex}. returns {Complex} */
Jmat.wrap = function(x, from, to) { return Jmat.Complex.wrap(Jmat.Complex.cast(x), Jmat.Complex.cast(from), Jmat.Complex.cast(to)); };
/* Clamp x between from and to (to included). x,to,from:{number|Complex}. returns {Complex} */
Jmat.clamp = function(x, from, to) { return Jmat.Complex.clamp(Jmat.Complex.cast(x), Jmat.Complex.cast(from), Jmat.Complex.cast(to)); };

// Other special functions

/* Minkowski's question mark function. x:{number|Complex}. returns {Complex} */
Jmat.minkowski = function(x) { return Jmat.Complex.minkowski(Jmat.Complex.cast(x)); };

// Matrix (NOTE: more are above: add, sub, mul, div, inv, neg, conj, exp, log, sqrt, cos, sin)

/* Eigenvalues and eigenvectors. m:{Array|Matrix}. returns {Object.<string, Matrix>} object with l:eigenvalues, v:eigenvectors */
Jmat.eig = function(m) { return Jmat.Matrix.eig(Jmat.Matrix.cast(m)); };
/* Moore-Penrose pseudo-inverse. m:{Array|Matrix}. returns {Matrix} */
Jmat.pseudoinverse = function(m) { return Jmat.Matrix.pseudoinverse(Jmat.Matrix.cast(m)); };
/* Determinant. m:{Array|Matrix}. returns {Complex} */
Jmat.determinant = function(m) { return Jmat.Matrix.determinant(Jmat.Matrix.cast(m)); };
/* Transpose. m:{Array|Matrix}. returns {Matrix} */
Jmat.transpose = function(m) { return Jmat.Matrix.transpose(Jmat.Matrix.cast(m)); };
/* Transjugate. m:{Array|Matrix}. returns {Matrix} */
Jmat.transjugate = function(m) { return Jmat.Matrix.transjugate(Jmat.Matrix.cast(m)); };
/* Trace. m:{Array|Matrix}. returns {Complex} */
Jmat.trace = function(m) { return Jmat.Matrix.trace(Jmat.Matrix.cast(m)); };
/* Rank. m:{Array|Matrix}. returns {Complex} */
Jmat.rank = function(m) { return Jmat.Matrix.rank(Jmat.Matrix.cast(m)); };
/* Adjoint aka Adjugate. m:{Array|Matrix}. returns {Matrix} */
Jmat.adj = function(m) { return Jmat.Matrix.adj(Jmat.Matrix.cast(m)); };
/* Fourier transform. m:{Array|Matrix}. returns {Matrix} */
Jmat.fft = function(m) { return Jmat.Matrix.fft(Jmat.Matrix.cast(m)); };
/* Inverse Fourier transform. m:{Array|Matrix}. returns {Matrix} */
Jmat.ifft = function(m) { return Jmat.Matrix.ifft(Jmat.Matrix.cast(m)); };
/* Solve AX=B. a,b:{Array|Matrix}. returns {Matrix} */
Jmat.solve = function(a, b) { return Jmat.Matrix.solve(Jmat.Matrix.cast(a), Jmat.Matrix.cast(b)); };
/* Dot product of two vectors. a,b:{Array|Matrix}. returns {Complex} */
Jmat.dot = function(a, b) { return Jmat.Matrix.dot(Jmat.Matrix.cast(a), Jmat.Matrix.cast(b)); };
/* Cross product of two size-3 vectors. m:{Array|Matrix}. returns {Matrix} */
Jmat.cross = function(a, b) { return Jmat.Matrix.cross(Jmat.Matrix.cast(a), Jmat.Matrix.cast(b)); };
/* Minor. a:{Array|Matrix}, row,col:{number} integer. returns {Complex} */
Jmat.minor = function(a, row, col) { return Jmat.Matrix.minor(Jmat.Matrix.cast(a), Jmat.Real.caststrict(row), Jmat.Real.caststrict(col)); };
/* Cofactor. a:{Array|Matrix}, row,col:{number} integer. returns {Complex} */
Jmat.cofactor = function(a, row, col) { return Jmat.Matrix.cofactor(Jmat.Matrix.cast(a), Jmat.Real.caststrict(row), Jmat.Real.caststrict(col)); };
/* Submatrix, x1 and y1 excluded, 0-based coordinates. a:{Array|Matrix}, y0,y1,x0,x1:{number} integer. returns {Matrix} */
Jmat.submatrix = function(a, y0, y1, x0, x1) { return Jmat.Matrix.submatrix(Jmat.Matrix.cast(a), Jmat.Real.caststrict(y0), Jmat.Real.caststrict(y1), Jmat.Real.caststrict(x0), Jmat.Real.caststrict(x1)); };

// Matrix decompositions

/* Singular value decomposition. m:{Array|Matrix}. returns {Object.<string, Matrix>} object with u:left vectors, s:singular values, v:right vectors */
Jmat.svd = function(m) { return Jmat.Matrix.svd(Jmat.Matrix.cast(m)); };
/* Spectral decomposition. m:{Array|Matrix}. returns {Object.<string, Matrix>} object with v:eigenvectors, d:eigenvalues on diagonal */
Jmat.eigd = function(m) { return Jmat.Matrix.eigd(Jmat.Matrix.cast(m)); };
/* QR decomposition. m:{Array|Matrix}. returns {Object.<string, Matrix>} object with q:Q matrix, r:R matrix */
Jmat.qr = function(m) { return Jmat.Matrix.qr(Jmat.Matrix.cast(m)); };

// Matrix norms

/* Frobenius norm. m:{Array|Matrix}. returns {Complex} */
Jmat.norm = function(m) { return Jmat.Matrix.norm(Jmat.Matrix.cast(m)); };
/* Spectral norm. m:{Array|Matrix}. returns {Complex} */
Jmat.norm2 = function(m) { return Jmat.Matrix.norm2(Jmat.Matrix.cast(m)); };
/* Maximum column norm. m:{Array|Matrix}. returns {Complex} */
Jmat.maxcolnorm = function(m) { return Jmat.Matrix.maxcolnorm(Jmat.Matrix.cast(m)); };
/* Maximum row norm. m:{Array|Matrix}. returns {Complex} */
Jmat.maxrownorm = function(m) { return Jmat.Matrix.maxrownorm(Jmat.Matrix.cast(m)); };

// Numerical algorithms. The function parameters f and df must always work with Jmat.Complex objects, one input, one output.

/* Derivative of f at x. x:{number|Complex}, f:{function(Complex):Complex}. returns {Complex} */
Jmat.differentiate = function(x, f) { return Jmat.Complex.differentiate(Jmat.Complex.cast(x), f); }; // derivative
/* Quadrature: definite integral of f from x to y. x,y:{number|Complex}, f:{function(Complex):Complex}, steps:{number} integer. returns {Complex} */
Jmat.integrate = function(x, y, f, steps) { return Jmat.Complex.integrate(Jmat.Complex.cast(x), Jmat.Complex.cast(y), f, Jmat.Real.caststrict(steps)); };
/* Secant method. f:{function(Complex):Complex}, z0:{number|Complex}, maxiter:{number} integer. returns {Complex} */
Jmat.rootfind_secant = function(f, z0, maxiter) { return Jmat.Complex.rootfind_secant(f, Jmat.Complex.cast(z0), Jmat.Real.caststrict(maxiter)); };
/* Newton's method. f,df:{function(Complex):Complex} df is derivative of f, z0:{number|Complex}, maxiter:{number} integer. returns {Complex} */
Jmat.rootfind_newton = function(f, df, z0, maxiter) { return Jmat.Complex.rootfind_secant(f, df, Jmat.Complex.cast(z0), Jmat.Real.caststrict(maxiter)); };

////////////////////////////////////////////////////////////////////////////////

// Test if input is a matrix, this makes it decide e.g. whether to call the matrix-specific or complex number specific function if there is functionname collision
Jmat.matrixIn_ = function(v) {
  if(!v) return false;
  if(typeof v == 'string') {
    if(v.length < 2) return false; // smallest possible valid matrix string is []
    if(v[0] == '[') return true;
    if(v[0] == '-') return false;
    if(v.charCodeAt(0) >= 48 && v.charCodeAt(0) <= 57) return false;
    return v.indexOf('[') != -1;
  }
  return v && (v instanceof Jmat.Matrix || v.length != undefined);
};

// Nice string of any known object, like Complex and Matrix and the objects of named matrices returned by svd or array returned by factorize
Jmat.toString = function(a) {
  if(!a) return '' + a;
  var result = '';
  // For arrays of known types
  if(typeof a == 'object' && Array.isArray(a)) {
    result += '[';
    for(var i = 0; i < a.length; i++) result += (Jmat.toString(a[i]) + (i + 1 == a.length ? '' : ', '));
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
      result += Jmat.toString(a[it]);
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
Jmat.Real: real math. This contains only a small amount of functions, most is
implemented in Jmat.Complex (and Jmat.Matrix for the matrices) further down.

Aliased as simply "Real" at the end of the file - disable that if it causes name clashes

These work on actual JS numbers rather than Jmat.Complex objects though, so Jmat.Real is
as easy to use as the standard JS Math library. Functions existing in Math, such
as cos and exp, are not implemented here for that reason.
*/

//Constructor, but not to be actually used, just a namespace for real functions
Jmat.Real = function() {
};

// cast all known numeric types to JS number
Jmat.Real.cast = function(v) {
  if(v && v.re != undefined) return v.re;
  if(v == undefined) return 0;
  return v;
};

// cast all known numeric types to JS number, but only if real (so complex/imag gives NaN)
Jmat.Real.caststrict = function(v) {
  if(v && v.re != undefined) return v.im == 0 ? v.re : NaN;
  if(v == undefined) return 0;
  return v;
};

Jmat.Real.SQRT2 = Math.sqrt(2);
Jmat.Real.SQRTPI = Math.sqrt(Math.PI); // gamma(0.5)
Jmat.Real.EM = 0.57721566490153286060; // Euler-Mascheroni constant
Jmat.Real.APERY = 1.2020569; // Apery's constant, zeta(3)

////////////////////////////////////////////////////////////////////////////////
// Categories
////////////////////////////////////////////////////////////////////////////////

Jmat.Real.isInt = function(x) {
  return x == Math.floor(x);
};

Jmat.Real.isPositiveInt = function(x) {
  return x == Math.floor(x) && x > 0;
};

Jmat.Real.isNegativeInt = function(x) {
  return x == Math.floor(x) && x < 0;
};

Jmat.Real.isPositiveIntOrZero = function(x) {
  return x == Math.floor(x) && x >= 0;
};

Jmat.Real.isNegativeIntOrZero = function(x) {
  return x == Math.floor(x) && x <= 0;
};

// x is odd integer
Jmat.Real.isOdd = function(x) {
  return Math.abs(x % 2) == 1; //works for negative x too
};

// x is even integer
Jmat.Real.isEven = function(x) {
  return x % 2 == 0; //works for negative x too
};

//isnanorinf isinfornan
Jmat.Real.isInfOrNaN = function(x) {
  return x == Infinity || x == -Infinity || isNaN(x);
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Real.dist = function(a, b) {
  return Math.abs(a - b);
};

// works for non-integers too
Jmat.Real.mod = function(a, b) {
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

// Remainder. This is just the % operator. It is here only for reference. Compare with Jmat.Real.mod, which is different.
Jmat.Real.rem = function(a, b) {
  return a % b;
};

// to is not included in the range
Jmat.Real.wrap = function(x, from, to) {
  if(from == to) return from;
  var m0 = Math.min(from, to);
  var m1 = Math.max(from, to);
  return m0 + Jmat.Real.mod(x - m0, m1 - m0);
};

// to is included in the range
Jmat.Real.clamp = function(x, from, to) {
  var m0 = Math.min(from, to);
  var m1 = Math.max(from, to);
  return Math.max(m0, Math.min(m1, x));
};

//Inspired by Wikipedia, Lanczos approximation, precision is around 15 decimal places
Jmat.Real.gamma = function(z) {
  // Return immediately for some common values, to avoid filling the cache with those
  if(z == Infinity) return Infinity;
  if(Jmat.Real.useFactorialLoop_(z - 1)) {
    return Jmat.Real.factorial(z - 1); //that one uses memoization
  }
  if(z == 0.5) return Jmat.Real.SQRTPI;

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

Jmat.Real.factorialmem_ = [1]; //memoization for factorial of small integers

Jmat.Real.useFactorialLoop_ = function(x) {
  return Jmat.Real.isPositiveIntOrZero(x) && x < 200;
};

Jmat.Real.factorial = function(a) {
  if(!Jmat.Real.useFactorialLoop_(a)) {
    return Jmat.Real.gamma(a + 1);
  }

  if(Jmat.Real.factorialmem_[a]) return Jmat.Real.factorialmem_[a];

  var result = Jmat.Real.factorialmem_[Jmat.Real.factorialmem_.length - 1];
  for(var i = Jmat.Real.factorialmem_.length; i <= a; i++) {
    result *= i;
    Jmat.Real.factorialmem_[i] = result;
  }
  return result;
};

Jmat.Real.firstPrimes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

// Initial set up shared by several of the prime test functions.
// Returns 0 if not prime, 1 if prime, NaN if problem, -1 if unknown by this function
Jmat.Real.isPrimeInit_ = function(n) {
  if(n == Infinity || n != n) return NaN;
  if(n != Math.round(n)) return 0;
  if(n < 2) return 0;
  if(n > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense (decimal: 9007199254740992)
  for(var i = 0; i < Jmat.Real.firstPrimes_.length; i++) {
    if(n == Jmat.Real.firstPrimes_[i]) return 1;
    if(n % Jmat.Real.firstPrimes_[i] == 0) return 0;
  }
  return -1;
};

//Returns 1 if prime, 0 if not prime, NaN if error. Naive slow algorithm. However, faster than miller rabin for n < 1500000
Jmat.Real.isPrimeSlow_ = function(n) {
  // Ensures the number is odd and integer in supported range, tested against first known primes
  var init = Jmat.Real.isPrimeInit_(n);
  if(init != -1) return init;

  var p = Jmat.Real.firstPrimes_[Jmat.Real.firstPrimes_.length - 1];
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
Jmat.Real.isPrimeMillerRabin_ = function(n) {
  // Ensures the number is odd and integer in supported range, tested against first known primes
  var init = Jmat.Real.isPrimeInit_(n);
  if(init != -1) return init;

  // Miller-Rabin test
  var base;
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
    //if(y != 1) return false;
    return y != 1;
  };

  for(var i = 0; i < base.length; i++) {
    if(!witness(n, s, d, base[i])) return 0;
  }
  return 1;
};
/*
test in console for the above function:
function benchfun(n) {
  var ra = 0; var ta0 = new Date().getTime(); for(var i = 0; i < n; i++) ra += Jmat.Real.isPrimeMillerRabin_(i); var ta1 = new Date().getTime();
  var rb = 0; var tb0 = new Date().getTime(); for(var i = 0; i < n; i++) rb += Jmat.Real.isPrimeSlow_(i); var tb1 = new Date().getTime();
  var rc = 0; var tc0 = new Date().getTime(); for(var i = 0; i < n; i++) rc += Jmat.Real.isPrime(i); var tc1 = new Date().getTime();
  console.log('fast: ' + (ta1 - ta0) + ' slow: ' + (tb1 - tb0) + ' both: ' + (tc1 - tc0) + ' test: ' + ra + ' = ' + rb + ' = ' + rc);
};
benchfun(100000);
--> it will report that slow if faster than miller rabin. That's because miller rabin is only faster for very large numbers. E.g. here you can see that miller rabin is faster:
Jmat.Real.isPrimeSlow_(4444280714420857)
Jmat.Real.isPrimeMillerRabin_(4444280714420857)


function testfun(n) {
  for(var i = 0; i < n; i++) {
    var a = Jmat.Real.isPrimeMillerRabin_(i);
    var b = Jmat.Real.isPrimeSlow_(i);
    var c = Jmat.Real.isPrime(i);
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
Jmat.Real.isPrime = function(n) {
  // below that, the "slow" method is faster. For higher values, Miller Rabin becomes more and more significantly faster.
  // However, for values above 0x010000000000000, a sum in the miller rabin overflows, so does not work either
  // ==> Jmat.Real.isPrime(9007199254740881) is noticeably slower than Jmat.Real.isPrime(4444280714420857)
  return (n < 1500000 || n > 0x010000000000000) ? Jmat.Real.isPrimeSlow_(n) : Jmat.Real.isPrimeMillerRabin_(n);
};

//for factorize
Jmat.Real.smallestPrimeFactor = function(x) {
  if(x == Infinity || x != x) return NaN;
  if(x != Math.round(x)) return NaN;
  if(x < 1) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense
  if(x == 1) return 1;
  for(var i = 0; i < Jmat.Real.firstPrimes_.length; i++) {
    if(x == Jmat.Real.firstPrimes_[i]) return x;
    if(x % Jmat.Real.firstPrimes_[i] == 0) return Jmat.Real.firstPrimes_[i];
  }
  var p = Jmat.Real.firstPrimes_[Jmat.Real.firstPrimes_.length - 1];
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
Jmat.Real.factorize = function(x) {
  var x = Math.round(x);
  var result = [];
  for(;;) {
    if(x < 1) break;
    var y = Jmat.Real.smallestPrimeFactor(x);
    result.push(y);
    if(x == y) break;
    x = Math.round(x / y);
  }
  return result;
};


Jmat.Real.primeCount = function(value) {
  var primesN = [ 0, 2, 3, 5, 7, 11, 13, 17 ];
  // Nth prime (1-indexed: n=1 gives 2)
  var p = function(n) {
    if(n < primesN.length) return primesN[n];
    var i = primesN[primesN.length - 1] + 2;
    var count = primesN.length - 1;
    for(;;) {
      if(Jmat.Real.isPrime(i)) {
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

Jmat.Real.nearestPrime = function(value) {
  var x = Math.round(value);
  if(x < 2) return 2;
  if(x == Infinity || x != x) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense

  if(Jmat.Real.isPrime(x)) return x;
  var i = x % 2 == 0 ? 1 : 2;
  while(true) {
    if(Jmat.Real.isPrime(x + i)) return x + i;
    if(Jmat.Real.isPrime(x - i)) return x - i;
    i += 2; //TODO: make faster (e.g. using += 6 at the least)
  }
};

Jmat.Real.nextPrime = function(value) {
  var x = Math.floor(value);
  if(x < 2) return 2; //the calculations below would give 3 instead of 2
  if(x == Infinity || x != x) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense

  var i = x % 2 == 0 ? 1 : 2;
  while(true) {
    if(Jmat.Real.isPrime(x + i)) return x + i;
    i += 2; //TODO: make faster (e.g. using += 6 at the least)
  }
};

Jmat.Real.previousPrime = function(value) {
  var x = Math.ceil(value);
  if(x <= 3) return 2; //avoid infinite loop
  if(x == Infinity || x != x) return NaN;
  if(x > 0x020000000000000) return NaN; //too large for the floating point's integer precision, result will not make sense

  var i = x % 2 == 0 ? 1 : 2;
  while(true) {
    if(Jmat.Real.isPrime(x - i)) return x - i;
    i += 2; //TODO: make faster (e.g. using += 6 at the least)
  }
};

Jmat.Real.eulerTotient = function(value) {
  if(value <= 0) return NaN;
  var n = Math.floor(value);
  var f = Jmat.Real.factorize(n);
  var prev = -1;
  var result = n;
  for(var i = 0; i < f.length; i++) {
    if(prev == f[i]) continue; //must be unique factors
    if(f[i] == 1) break;
    prev = f[i];
    result *= (1 - (1 / f[i]));
  }
  return result;
};

// The first integer binomials, allows fast calculation of those by just looking up in the array, e.g. binomial(5, 8) = Jmat.Complex.pascal_triangle_cache_[5][8]
// some rows area pre-filled to start it off (just pre-filling the first one would be sufficient normally, the rest is just for the shows)
Jmat.Real.pascal_triangle_cache_ = [
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
Jmat.Real.pascal_triangle = function(n, p) {
  if(n < 0 || p < 0 || n < p) return NaN;
  if(n > 180) return NaN; //triangle values get too big for integers in double precision floating point
  //fill up cache if needed
  var t = Jmat.Real.pascal_triangle_cache_;
  while(t.length <= n) {
    var l = t.length; //the 'n' of the new row
    var l2 = l + 1; // number of elements of this new row
    t[l] = [];
    for(var i = 0; i < l2; i++) {
      t[l][i] = (i == 0 || i == l2 - 1) ? 1 : (t[l-1][i-1] + t[l-1][i]);
    }
  }
  return t[n][p];
};

//greatest common divisor
Jmat.Real.gcd = function(x, y) {
 //Euclid's algorithm
 while(true) {
   if(y == 0) return x;
   var z = Jmat.Real.mod(x, y);
   x = y;
   y = z;
 }
};

//least common multiple
Jmat.Real.lcm = function(x, y) {
 return Math.abs(x * y) / Jmat.Real.gcd(x, y);
};

// Decomposes fraction (aka rational approximation): returns two integers [numerator, denominator] such that n/d = a.
// Very slow, too slow for inner loop of running programs (integrate or complex plot)...
// max = max value for denominator
// E.g. Jmat.Real.decompose(Math.PI, 100) gives [22, 7], because 22/7 approximates pi.
Jmat.Real.decompose = function(x, max) {
  if(!max) max = 100000;
  var neg = (x < 0);
  if(neg) x = -x;
  var f = Math.floor(x);
  var y = x - f;

  if(y == 0) return [x, 1]; //otherwise the loop will run max times for nothing, very inefficient

  var result;

  var a = 0;
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
    if (b > max) result = [c, d];
    else result = [a, b];
  }

  result[0] += f * result[1];
  if(neg) result[0] = -result[0];

  return result;
};

// Hybrid between decompose and decomposeFast
Jmat.Real.decomposeSemiFast = function(x, max) {
  var maxslow = 1000;
  if(max < maxslow) {
    return Jmat.Real.decompose(x, max);
  } else {
    var a = Jmat.Real.decompose(x, maxslow);
    var ax = a[0] / a[1];
    if(ax == x) return a;
    var b = Jmat.Real.decomposeFast(x, maxslow);
    var bx = b[0] / b[1];
    return (Math.abs(x - ax) < Math.abs(x - bx)) ? a : b;
  }
};

// Decomposes fraction (aka rational approximation): returns two integers [numerator, denominator] such that n/d = a.
// max = max value for denominator
// A lot faster, but less nice than Jmat.Real.decompose (e.g. returns 83333/10000 instead of 1/12), and with high preference for decimal denominators. TODO: other bases than base 10
Jmat.Real.decomposeFast = function(x, max) {
  if(!max) max = 100000;
  var max1 = max - 1;

  if(x <= max1 && x >= -max1 && (x < -1.0 / max1 || x > 1.0 / max1)) {
    var neg = (x < 0);
    if(neg) x = -x;
    var z = Math.floor(x);
    var n = x - z;
    var d = max;
    n = Math.floor(n * d);
    var g = Jmat.Real.gcd(n, d);
    d /= g;
    n /= g;
    n += z * d;
    if(neg) n = -n;
    // n = numerator, d = denominator
    return [n, d];
  }
  return [x, 1];
};

Jmat.Real.near = function(x, y, precision) {
  // works also for infinities
  return x >= y - precision && x <= y + precision;
};

// Fractional part of x, x - floor(x). NOTE: this variant gives positive results for negative x
Jmat.Real.frac = function(x) {
  return x - Math.floor(x);
};

// Fractional part of x, x - int(x). NOTE: this variant gives negative results for negative x
Jmat.Real.fracn = function(x) {
  return x > 0 ? (x - Math.floor(x)) : -(-x - Math.floor(-x));
};

// Only the principal branch for real values above -1/e
Jmat.Real.lambertw = function(x) {
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
    //Since the above method works only up to 703, use some kind of binary search instead (it's a monotonously increasing function at this point)
    // TODO: probably just use Halley's method here instead
    var step = 1;
    var lastDir = 0;
    var result = Math.log(x) - Math.log(Math.log(x)); // good starting value speeds up iterations. E.g. only 76 instead of 292 for 7e100.
    for(;;) {
      if(step == 0 || step * 0.5 == step || result + step == result) return result; //avoid infinite loop
      var v = result * Math.exp(result);
      if(Jmat.Real.near(v, x, 1e-15)) return result;
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
  }
  return NaN;
};

// Tetration
// Returns experimental (not mathematically correct) results unless x is an integer or Infinity
Jmat.Real.tetration = function(a, x) {
  // if(a == 1) return 1; // Not true, e.g. 1^Infinity is not 1. Probably things are not like this for tetration either.
  if(x == 0) return 1; //by definition
  if(x == 1) return a;
  if(x == 2) return Math.pow(a, a);
  if(a >= 2 && x > 5) return Infinity; // too big for double
  if(a == 0 && Jmat.Real.isPositiveInt(x)) {
    // Chosen definition that considers 0^0 = 1, so 0^^n is 1 if n is even, 0 if n is odd.
    return Jmat.Complex.isEven(x) ? 1 : 0;
  }

  // Power tower (infinitely iterated exponentiation)
  if(x == Infinity && a > 0) {
    // converges if a >= 0.066 && a <= 1.44
    var l = Math.log(a);
    return Jmat.Real.lambertw(-l) / (-l);
  }

  var runloop = function(a, b, num, l) {
    var result = b;
    var last;
    for(var i = 0; i < num; i++) {
      if(l) result = Jmat.Real.logy(result, a);
      else result = Math.pow(a, result);
      if(isNaN(result)) return result;
      if(result == Infinity) return result; // Actually redundant, result == last already checks that too
      if(result == last) return result; // E.g. a=1.01 already converges to 1.0101015237405409 after just 7 iterations. 1.44 converges to 2.393811748202943 after 244 iterations, so numbers near that take a while to loop. 1.45 and higher converge to Infinity.
      last = result;
      if(i > 1000) return NaN; //avoid infinite loop
    }
    return result;
  };

  // This supports power tower as well
  // It supports it better than the lambertw formula above because the domain here is real numbers only, the lambertw formula fails on, for example, -1 and returns nan instead of -1
  if(Jmat.Real.isPositiveInt(x)) {
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
Jmat.Real.logy = function(x, y) {
  return Math.log(x) / Math.log(y);
};

// dawson function D+(x) aka F(x). erfi(x) = 2/sqrt(pi) * exp(x*x) * D(x)
// so while erfi overflows for large args, this one doesn't due to no exp(x*x)
Jmat.Real.dawson = function(x) {
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
Jmat.Real.erfi = function(x) {
  var neg = false;
  if(x < 0) {
    x = -x;
    neg = true;
  }
  var result = 0;
  var ps = 1.0 / Jmat.Real.SQRTPI;

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
    result = 2 / Jmat.Real.SQRTPI * Math.exp(x * x) * Jmat.Real.dawson(x);
  }

  if(neg) result = -result;
  return result;
};

Jmat.Real.erf = function(x) {
  var neg = x < 0;
  if(neg) x = -x;

  if (x == 0) return 0;
  var t = 1 / (1 + 0.3275911 * x);
  var p = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  var result = 1.0 - p * Math.exp(-(x*x));

  if(neg) result = -result;
  return result;
};

Jmat.Real.erfc = function(x) {
  var neg = x < 0;
  if(neg) x = -x;
  var result;

  if(x <= 0.5) {
    var x2 = x * x;
    var x3 = x * x2;
    var x5 = x3 * x2;
    var x7 = x5 * x2;
    result = 1 - 2 / Jmat.Real.SQRTPI * (x - x3 / 3 + x5 / 10 + x7 / 42);
    //result = Math.exp(-x*x) / 6 + Math.exp(-0.75 * x * x) / 2;
  } else if (x >= 5) {
    // asymptotic expansion for large real x
    var x2 = x * x;
    var x4 = x2 * x2;
    var x6 = x4 * x2;
    var x8 = x6 * x2;
    result = Math.exp(-(x*x)) / (x * Jmat.Real.SQRTPI) * (1 - 1/2.0/x2 + 3/4.0/x4 - 15/8.0/x6 + 105/16.0/x8);
  } else {
    var t = 1 / (1 + 0.3275911 * x);
    var p = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    result = p * Math.exp(-(x*x));
  }

  if(neg) result = 2 - result;
  return result;
};

//Minkowski's question mark function, from Wikipedia
Jmat.Real.minkowski = function(x) {
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
Jmat.Real.dms = function(a) {
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
Jmat.Real.dd = function(a) {
  var neg = a < 0;
  if(neg) a = -a;

  var deg = Math.floor(a);
  var mins = Math.floor((a * 100 - deg * 100));
  var sec = Math.floor(a * 10000 - deg * 10000 - mins * 100);

  var result = deg + mins / 60.0 + sec / 3600.0;

  if(neg) result = -result;
  return result;
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Real.isLeapYear = function(y) {
  return y % 400 == 0 || (y % 4 == 0 && y % 100 != 0);
};

Jmat.Real.montharray_ = [-1 /*there is no month 0*/, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; //february must be adjusted for leap year

Jmat.Real.monthLength = function(month, leap) {
  return (leap && month == 2) ? 29 : Jmat.Real.montharray_[month];
};


//number of days since the first day of year 0. 1 january of the year 0 is 0. 2 january is 1, etc...
//only for Gregorian calendar, does not take Julian calendar into account
Jmat.Real.numDaysSince0 = function(year, month, day) {
  //calculate number of leap years before this year (year 0 is considered leap)
  var num400 = Math.floor((year - 1) / 400) + 1;
  var num100 = Math.floor((year - 1) / 100) + 1 - num400;
  var num4 = Math.floor((year - 1) / 4) + 1 - num100 - num400;
  var numleap = num4 + num400;
  var yeardays = numleap * 366 + (year - numleap) * 365; //days of years before this year

  var monthdays = 0; //days of years before this month
  var leap = Jmat.Real.isLeapYear(year);
  for (var i = 1; i < month; i++) monthdays += Jmat.Real.monthLength(i, leap);

  return yeardays + monthdays + day - 1; //-1 because day 1 is in fact zero
};

//converts number of days since the year 0, to year/month/day
//returns array [year, month, day]
//only for Gregorian calendar, does not take Julian calendar into account
Jmat.Real.daysSince0ToDate = function(days) {
  //every 400 years there are 97 leap years. So a year is 365.2425 days in average.
  var year = Math.floor(days / 365.2425);
  var leap = Jmat.Real.isLeapYear(year);

  days -= Jmat.Real.numDaysSince0(year, 1, 1);

  var month = 0;
  for (var i = 1; i <= 12; i++) {
    month++;
    var num = Jmat.Real.monthLength(i, leap);
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

Aliased as simply "Complex" at the end of the file - disable that if it causes name clashes

This is the actual object used as complex number. In addition, most of the
functions are implemented as static functions in here.

The only sad thing is that Javascript doesn't support operator overloading
and nice expressions like a + b have to become a.add(b) instead.
*/
Jmat.Complex = function(re, im) {
  if(this instanceof Jmat.Complex) {
    // Keyword "new" in front. Does not do any checks, to be "fast"
    this.re = re;
    this.im = im;
  } else {
    // No keyword "new" in front, use the convenience factory function instead
    return Jmat.Complex.make(re, im); // This supports several argument types
  }
};

// TODO: define a bit better which combinations of Infinity/Nan/... in re and im mean what (E.g. re and im both Infinity means "undirected infinity", already used by gamma function but by nothing else)

// Create a new Jmat.Complex value. Copies Jmat.Complex if a Jmat.Complex is given as first argument
// with 0 arguments, creates zero value. With a and b numbers, creates complex number from it. With a Jmat.Complex object, copies it.
// the first parameter must be given and be number or Jmat.Complex. The second parameter is optional.
Jmat.Complex.make = function(a, b) {
  if(typeof a == 'number') return new Jmat.Complex(a, b == undefined ? 0 : b);
  if(typeof a == 'string') return Jmat.Complex.parse(a);
  return new Jmat.Complex(a.re, a.im); // Copy value object
};

// Create a new Jmat.Complex value, real
Jmat.Complex.newr = function(re) {
  return new Jmat.Complex(re, 0);
};

// Create a new Jmat.Complex value, imaginary
Jmat.Complex.newi = function(im) {
  return new Jmat.Complex(0, im);
};

// Create a new Jmat.Complex value, polar
Jmat.Complex.polar = function(r, a) {
  return new Jmat.Complex(r * Math.cos(a), r * Math.sin(a));
};

// Casts the given number type to Jmat.Complex. If the given type is already of type Jmat.Complex, does not copy it but returns the input.
// TODO: also support strings of the form '5+6i', and be able to parse them
Jmat.Complex.cast = function(v) {
  if(v && v.re != undefined) return v;
  if(v == undefined) return Jmat.Complex(0);
  return Jmat.Complex(v);
};

// Because JS number toFixed appends zeros
Jmat.Complex.formatFloat_ = function(value, precision) {
  var power = Math.pow(10, precision || 0);
  return String(Math.round(value * power) / power);
};

//debugstring
Jmat.Complex.toString = function(value, opt_precision) {
  if(!value) return value == 0 ? 'invalid0' : ('' + value);
  var re = (opt_precision ? Jmat.Complex.formatFloat_(value.re, opt_precision) : ('' + value.re));
  var im = (opt_precision ? Jmat.Complex.formatFloat_(value.im, opt_precision) : ('' + value.im));

  if(value.im == 0 || im == '0') return '' + re;
  if(value.re == 0) return '' + im + 'i';
  if(value.im < 0) return '' + re + im + 'i';
  return '' + re + '+' + im + 'i';
};
Jmat.Complex.prototype.toString = function(opt_precision) {
  return Jmat.Complex.toString(this, opt_precision);
};

// Parses strings of the form '5', '5+i', '5-2.3i', '1.25e-25+17.37e5i'
Jmat.Complex.parse = function(text) {
  var i = text.indexOf('i');
  if(i == -1) {
    return Jmat.Complex(parseFloat(text));
  } else {
    if(text == 'i') return Jmat.Complex(0, 1);
    text = text.substr(0, i); // remove the i and anything after it
    text = text.replace(/ /g, ''); // support forms with spaces like '5 + 2i' too
    // Make it handle it correctly if just 'i' without number in front is used
    if(text[i - 1] == '+' || text[i - 1] == '-') text += '1';
    
    // Find the + or - which is not after an 'e' or 'E'
    for(var j = 1; j < text.length; j++) {
      if((text[j] == '+' || text[j] == '-') && !(text[j - 1] == 'e' || text[j - 1] == 'E')) {
        return Jmat.Complex(parseFloat(text.substr(0, j)), parseFloat(text.substr(j)));
      }
    }
    return Jmat.Complex(0, parseFloat(text)); // pure imaginary
  }
};

// Only use these as constants, never modify these, never return them!
Jmat.Complex.ZERO = Jmat.Complex(0);
Jmat.Complex.ONE = Jmat.Complex(1);
Jmat.Complex.TWO = Jmat.Complex(2);
Jmat.Complex.I = Jmat.Complex.newi(1);
Jmat.Complex.PI = Jmat.Complex(Math.PI);
Jmat.Complex.E = Jmat.Complex(Math.E);
Jmat.Complex.SQRT2 = Jmat.Complex(Math.sqrt(2));
Jmat.Complex.SQRTPI = Jmat.Complex(Math.sqrt(Math.PI));
Jmat.Complex.INVSQRT2PI = Jmat.Complex(1 / Math.sqrt(2 * Math.PI)); //0.3989422804014327
Jmat.Complex.EM = Jmat.Complex(Jmat.Real.EM); // Euler-Mascheroni constant
Jmat.Complex.APERY = Jmat.Complex(Jmat.Real.APERY); // Apery's constant, zeta(3)

Jmat.Complex.real = function(z) {
  return Jmat.Complex(z.re);
};
Jmat.Complex.prototype.real = function() {
  return Jmat.Complex(this.re);
};

Jmat.Complex.imag = function(z) {
  return Jmat.Complex(z.im);
};
Jmat.Complex.prototype.imag = function() {
  return Jmat.Complex(this.im);
};

//Basic operators

Jmat.Complex.add = function(x, y) {
  return new Jmat.Complex(x.re + y.re, x.im + y.im);
};
Jmat.Complex.prototype.add = function(y) {
  return new Jmat.Complex(this.re + y.re, this.im + y.im);
};

Jmat.Complex.sub = function(x, y) {
  return new Jmat.Complex(x.re - y.re, x.im - y.im);
};
Jmat.Complex.prototype.sub = function(y) {
  return new Jmat.Complex(this.re - y.re, this.im - y.im);
};

Jmat.Complex.mul = function(x, y) {
  if(x.im == 0 && y.im == 0) {
    return new Jmat.Complex(x.re * y.re, 0);
  } else {
    var re = x.re * y.re - x.im * y.im;
    var im = x.im * y.re + x.re * y.im;
    return new Jmat.Complex(re, im);
  }
};
Jmat.Complex.prototype.mul = function(y) {
  return Jmat.Complex.mul(this, y);
};

Jmat.Complex.div = function(x, y) {
  if(x.im == 0 && y.im == 0) {
    return new Jmat.Complex(x.re / y.re, 0);
  } else {
    if(Jmat.Complex.isInf(x) && !Jmat.Complex.isInfOrNaN(y)) {
      // Result should be some infinity (because it's infinity divided through finite value), but the formula below would give a NaN somewhere.
      // 4 possible rotations of the infinity, based on quadrant of y (TODO: THIS IS IGNORED NOW!!)
      return x;
    }
    var d = y.re * y.re + y.im * y.im;
    if(d == Infinity || d == -Infinity) {
      // the calculations below would give Infinity/Infinity = NaN even though result should be 0.
      if(!Jmat.Complex.isInfOrNaN(x)) return Jmat.Complex(0);
    }
    if(d == 0 && !Jmat.Complex.isInfOrNaN(x) && (x.re != 0 || x.im != 0)) {
      // the calculations below would give 0/0 = NaN even though result should be some infinity.
      return new Jmat.Complex(x.re == 0 ? 0 : (x.re < 0 ? -Infinity : Infinity), x.im == 0 ? 0 : (x.im < 0 ? -Infinity : Infinity));
    }
    var re = (x.re * y.re + x.im * y.im) / d;
    var im = (x.im * y.re - x.re * y.im) / d;
    return new Jmat.Complex(re, im);
  }
};
Jmat.Complex.prototype.div = function(y) {
  return Jmat.Complex.div(this, y);
};

Jmat.Complex.addr = function(z, a) {
  return Jmat.Complex(z.re + a, z.im);
};
Jmat.Complex.prototype.addr = function(a) {
  return Jmat.Complex(this.re + a, this.im);
};

Jmat.Complex.subr = function(z, a) {
  return Jmat.Complex(z.re - a, z.im);
};
Jmat.Complex.prototype.subr = function(a) {
  return Jmat.Complex(this.re - a, this.im);
};
Jmat.Complex.rsub = function(a, z) {
  return Jmat.Complex(a - z.re, -z.im);
};
Jmat.Complex.prototype.rsub = function(a) {
  return Jmat.Complex(a - this.re, -this.im);
};

Jmat.Complex.mulr = function(z, a) {
  return Jmat.Complex(z.re * a, z.im * a);
};
Jmat.Complex.prototype.mulr = function(a) {
  return Jmat.Complex(this.re * a, this.im * a);
};

Jmat.Complex.divr = function(z, a) {
  return Jmat.Complex(z.re / a, z.im / a);
};
Jmat.Complex.prototype.divr = function(a) {
  return Jmat.Complex(this.re / a, this.im / a);
};
Jmat.Complex.rdiv = function(a, z) {
  return Jmat.Complex.div(Jmat.Complex(a), z);
};

//rotate complex number z by a radians. That is, change its argument. a is real (JS number).
Jmat.Complex.rotate = function(z, a) {
  if(a == 0) return z;
  return Jmat.Complex.polar(z.absr(), z.argr() + a);
};

//rotate complex number z by 2pi/n radians. This results in giving the next solution of the nth root.
Jmat.Complex.nextroot = function(z, n) {
  var result = Jmat.Complex.rotate(z, Math.PI * 2 / n);
  if(Jmat.Real.near(result.im, 0, 1e-14)) result.im = 0;
  return result;
};

// mod operation, result has the sign of the divisor (unlike % operator in JS, Java and C99), so it's like wrapping x in range 0..y.
// works on real or complex numbers too, e.g. (6+4i) mod (3+5i) gives (-2+2i)
Jmat.Complex.mod = function(x, y) {
  if(x.im != 0 || y.im != 0) return x.sub(Jmat.Complex.floor(x.div(y)).mul(y));
  return Jmat.Complex(Jmat.Real.mod(x.re, y.re));
};

// remainder operation, like the % operator in JS, Java and C99.
Jmat.Complex.rem = function(x, y) {
  if(x.im != 0 || y.im != 0) return x.sub(Jmat.Complex.trunc(x.div(y)).mul(y));
  return Jmat.Complex(x.re % y.re);
};

Jmat.Complex.wrap = function(x, from, to) {
  return new Jmat.Complex(Jmat.Real.wrap(x.re, from.re, to.re), Jmat.Real.wrap(x.im, from.im, to.im));
};

Jmat.Complex.clamp = function(x, from, to) {
  return new Jmat.Complex(Jmat.Real.clamp(x.re, from.re, to.re), Jmat.Real.clamp(x.im, from.im, to.im));
};

Jmat.Complex.bitnot = function(x) {
  var result = Jmat.Complex(0);
  result.re = ~x.re;
  //imaginary part not affected on purpose: otherwise it appears when bit-inverting real number, which is in 99.9% of the cases not wanted
  //result.im = ~x.im;
  return result;
};

Jmat.Complex.bitand = function(x, y) {
  var result = Jmat.Complex(0);
  result.re = x.re & y.re;
  result.im = x.im & y.im;
  return result;
};

Jmat.Complex.bitor = function(x, y) {
  var result = Jmat.Complex(0);
  result.re = x.re | y.re;
  result.im = x.im | y.im;
  return result;
};

Jmat.Complex.bitxor = function(x, y) {
  var result = Jmat.Complex(0);
  result.re = x.re ^ y.re;
  result.im = x.im ^ y.im;
  return result;
};

Jmat.Complex.lshift = function(x, y) {
  var result = Jmat.Complex(0);
  result.re = x.re << y.re;
  result.im = x.im << y.im;
  return result;
};

Jmat.Complex.rshift = function(x, y) {
  var result = Jmat.Complex(0);
  result.re = x.re >> y.re;
  result.im = x.im >> y.im;
  return result;
};

Jmat.Complex.neg = function(x) {
  return Jmat.Complex(-x.re, -x.im);
};
Jmat.Complex.prototype.neg = function() {
  return Jmat.Complex(-this.re, -this.im);
};

// Returns 0 if z is 0, 1 if z is positive, -1 if z is negative. For complex z, returns z / abs(z)
Jmat.Complex.sign = function(z) {
  if (z.im == 0) {
    if(z.re == 0) return Jmat.Complex(0);
    else if(z.re < 0) return Jmat.Complex(-1);
    return Jmat.Complex(1);
  }

  return z.div(Jmat.Complex.abs(z));
};

// Returns 0 if z is 0, 1 if z is positive, -1 if z is negative. For complex z, returns sign of z.im if z.re == 0, sign of z.re otherwise (that is, the function returns sqrt(z*z) / z, except for z=0)
Jmat.Complex.csgn = function(z) {
  if (Jmat.Real.near(z.re, 0, 1e-15)) { //avoid numeric imprecisions for e.g. the values of e.g. acosh
    if(z.im == 0) return Jmat.Complex(0);
    else if(z.im < 0) return Jmat.Complex(-1);
    return Jmat.Complex(1);
  } else {
    if(z.re == 0) return Jmat.Complex(0);
    else if(z.re < 0) return Jmat.Complex(-1);
    return Jmat.Complex(1);
  }
};

Jmat.Complex.conj = function(x) {
  return Jmat.Complex(x.re, -x.im);
};
Jmat.Complex.prototype.conj = function() {
  return Jmat.Complex(this.re, -this.im);
};

Jmat.Complex.eq = function(x, y) {
  if(!x || !y) return x == y;
  return (x.re == y.re && x.im == y.im);
};
Jmat.Complex.prototype.eq = function(y) {
  return y && this.re == y.re && this.im == y.im;
};

Jmat.Complex.eqr = function(x, y) {
  return (x.re == y && x.im == 0);
};
Jmat.Complex.prototype.eqr = function(y) {
  return (this.re == y && this.im == 0);
};

Jmat.Complex.powr = function(z, a) {
  return Jmat.Complex.pow(z, Jmat.Complex(a));
};
Jmat.Complex.prototype.powr = function(a) {
  return Jmat.Complex.pow(this, Jmat.Complex(a));
};

Jmat.Complex.inv = function(z) {
  return Jmat.Complex.ONE.div(z);
};
Jmat.Complex.prototype.inv = function() {
  return Jmat.Complex.ONE.div(this);
};

//increment
Jmat.Complex.inc = function(z) {
  return new Jmat.Complex(z.re + 1, z.im);
};
Jmat.Complex.prototype.inc = function() {
  return new Jmat.Complex(this.re + 1, this.im);
};

//decrement
Jmat.Complex.dec = function(z) {
  return new Jmat.Complex(z.re - 1, z.im);
};
Jmat.Complex.prototype.dec = function(z) {
  return new Jmat.Complex(this.re - 1, this.im);
};

// absolute value squared, returned as real
Jmat.Complex.abssqr = function(x) {
  return x.re * x.re + x.im * x.im;
};
Jmat.Complex.prototype.abssqr = function() {
  return this.re * this.re + this.im * this.im;
};

// absolute value, aka modulus of complex number
Jmat.Complex.absr = function(x) {
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
Jmat.Complex.prototype.absr = function() {
  return Jmat.Complex.absr(this);
};

Jmat.Complex.abs = function(x) {
  var result = Jmat.Complex(0);
  if(x.im == 0) result.re = Math.abs(x.re);
  else result.re = Jmat.Complex.absr(x);
  return result;
};
Jmat.Complex.prototype.abs = function() {
  return Jmat.Complex.abs(this);
};

// returns the complex argument in range -PI to +PI
Jmat.Complex.argr = function(x) {
  if(x.im == 0) return x.re < 0 ? Math.PI : 0;
  return Math.atan2(x.im, x.re);
};
Jmat.Complex.prototype.argr = function() {
  return Jmat.Complex.argr(this);
};

Jmat.Complex.arg = function(z) {
  return Jmat.Complex(Jmat.Complex.argr(z));
};
Jmat.Complex.prototype.arg = function() {
  return Jmat.Complex.arg(this);
};

//returns result in range 0-1 rather than -PI to PI. Useful for graphical representations, not for math. 0 matches 0 degrees, 0.5 matches 180 degrees, 0.999 matches around 359 degrees.
Jmat.Complex.argr1 = function(z) {
  var result = Jmat.Complex.argr(z);
  if(result < 0) result += 2 * Math.PI;
  result /= (2 * Math.PI);
  if(result < 0) result = 0;
  if(result > 1) result = 1;
  return result;
};


////////////////////////////////////////////////////////////////////////////////
// Categories
////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.isReal = function(z) {
  return z.im == 0;
};

Jmat.Complex.isImaginary = function(z) {
  return z.re == 0;
};

Jmat.Complex.isInt = function(z) {
  return z.im == 0 && Jmat.Real.isInt(z.re);
};

// Gaussian integer
Jmat.Complex.isGaussian = function(z) {
  return Jmat.Real.isInt(z.re) && Jmat.Real.isInt(z.im);
};

Jmat.Complex.isNaN = function(z) {
  return !z || isNaN(z.re) || isNaN(z.im);
};

//is infinite
Jmat.Complex.isInf = function(z) {
  return (z.re * z.re + z.im * z.im) == Infinity;
};

//isnanorinf isinfornan
Jmat.Complex.isInfOrNaN = function(z) {
  return !z || Jmat.Real.isInfOrNaN(z.re) || Jmat.Real.isInfOrNaN(z.im);
};

//real and strictly positive
Jmat.Complex.isPositive = function(z) {
  return z.re > 0 && z.im == 0;
};

//real and strictly negative
Jmat.Complex.isNegative = function(z) {
  return z.re < 0 && z.im == 0;
};

Jmat.Complex.isPositiveOrZero = function(z) {
  return z.re >= 0 && z.im == 0;
};

Jmat.Complex.isNegativeOrZero = function(z) {
  return z.re <= 0 && z.im == 0;
};

//strictly positive
Jmat.Complex.isPositiveInt = function(z) {
  return Jmat.Complex.isInt(z) && z.re > 0;
};

//strictly negative
Jmat.Complex.isNegativeInt = function(z) {
  return Jmat.Complex.isInt(z) && z.re < 0;
};

Jmat.Complex.isPositiveIntOrZero = function(z) {
  return Jmat.Complex.isInt(z) && z.re >= 0;
};

Jmat.Complex.isNegativeIntOrZero = function(z) {
  return Jmat.Complex.isInt(z) && z.re <= 0;
};

// z is odd integer
Jmat.Complex.isOdd = function(z) {
  return Jmat.Complex.isInt(z) && Math.abs(z.re % 2) == 1;
};

// z is even integer
Jmat.Complex.isEven = function(z) {
  return Jmat.Complex.isInt(z) && z.re % 2 == 0;
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pow = function(x, y) {
  if(Jmat.Complex.isReal(x) && Jmat.Complex.isReal(y) && (x.re >= 0 || y.re == Infinity || y.re == -Infinity || Jmat.Real.isInt(y.re))) {
    //if(x.re == 0 && y.re == 0) return Jmat.Complex(NaN); // JS's pow returns 1 for 0^0
    // It is chosen to return 1 for 0^0, not NaN. NaN is mathematically more correct, however 0^0 is correct in many practical applications.
    return Jmat.Complex(Math.pow(x.re, y.re));
  } else {
    // This is just one branch. In fact it returns a complex result for -3 ^ (1/3),
    // the cube root of -3. To get the real result, use absolute value (and then negate) on it.
    // This is correct: the principal result of the cube root for this is a complex number.
    // Note: This returns incorrect values for a negative real to the power of Infinity: the result should be -Infinity for < -1, 0 for > -1, NaN for -1, but it always gives NaN. However, the "if" part above already handles that.
    var r = Jmat.Complex.absr(x);
    var t = Jmat.Complex.argr(x);
    var u = Math.pow(r, y.re) * Math.exp(-y.im * t);
    if(isNaN(u)) {
      u = Math.pow(1, y.re / r) * Math.exp(-y.im * t / r);
      if(u < 0) u = -Infinity;
      else if(u > 0) u = Infinity;
      else u = NaN;
    }
    var v = y.im * Math.log(r) + y.re * t;
    return Jmat.Complex(u * Math.cos(v), u * Math.sin(v));
  }
};
Jmat.Complex.prototype.pow = function(y) {
  return Jmat.Complex.pow(this, y);
};

Jmat.Complex.sin = function(z) {
  if(z.im == 0) return Jmat.Complex(Math.sin(z.re));

  var iz = Jmat.Complex(-z.im, z.re);
  var eiz = Jmat.Complex.exp(iz);
  var ieiz = Jmat.Complex.inv(eiz);
  return eiz.sub(ieiz).div(Jmat.Complex(0, 2));
};

//unnormalized sinc: sin(x) / x, but also defined for x = 0
Jmat.Complex.sinc = function(z) {
  if(z.eqr(0)) return Jmat.Complex(1);
  return Jmat.Complex.sin(z).div(z);
};

Jmat.Complex.cos = function(z) {
  if(z.im == 0) return Jmat.Complex(Math.cos(z.re));

  var iz = Jmat.Complex(-z.im, z.re);
  var eiz = Jmat.Complex.exp(iz);
  var ieiz = Jmat.Complex.inv(eiz);
  return eiz.add(ieiz).mulr(0.5);
};

Jmat.Complex.tan = function(z) {
  if(z.im == 0) return Jmat.Complex(Math.tan(z.re));

  var iz = Jmat.Complex(-z.im, z.re);
  var eiz = Jmat.Complex.exp(iz);
  var ieiz = Jmat.Complex.inv(eiz);
  return (eiz.sub(ieiz).div(Jmat.Complex(0, 2))).div(eiz.add(ieiz).mulr(0.5)); // Jmat.Complex.sin(z).div(Jmat.Complex.cos(z));
};

Jmat.Complex.asin = function(z) {
  if(z.im == 0 && z.re >= -1 && z.re <= 1) return Jmat.Complex(Math.asin(z.re));

  var s = Jmat.Complex.sqrt(Jmat.Complex.ONE.sub(z.mul(z)));
  var l = Jmat.Complex.log(Jmat.Complex(-z.im, z.re).add(s));
  return Jmat.Complex(l.im, -l.re);
};

Jmat.Complex.acos = function(z) {
  if(z.im == 0 && z.re >= -1 && z.re <= 1) return Jmat.Complex(Math.acos(z.re));

  //i * ln(x - i * sqrt(1-x^2))
  var s = Jmat.Complex.sqrt(Jmat.Complex.ONE.sub(z.mul(z))).mul(Jmat.Complex.I);
  var l = Jmat.Complex.log(z.add(s));
  return Jmat.Complex(l.im, -l.re);
};

Jmat.Complex.atan = function(z) {
  if(z.im == 0) return Jmat.Complex(Math.atan(z.re));

  var iz = Jmat.Complex(-z.im, z.re);
  var b = Jmat.Complex.ONE.sub(iz).div(iz.inc());
  var l = Jmat.Complex.log(b);
  return Jmat.Complex(-0.5 * l.im, 0.5 * l.re);
};

Jmat.Complex.atan2 = function(x, y) {
  if(!Jmat.Complex.isReal(x) || !Jmat.Complex.isReal(y)) {
    if(y.eqr(0)) return Jmat.Complex(Math.PI / 2);

    // For complex values, an alternate form of the definition can be used:
    // 2 * atan(y / (sqrt(x^2+y^2)+x))
    return Jmat.Complex.atan(Jmat.Complex.sqrt(x.mul(x).add(y.mul(y))).sub(x).div(y)).mulr(2);
  } else {
    var result = Jmat.Complex(0);
    result.re = Math.atan2(x.re, y.re);
    return result;
  }
};

Jmat.Complex.sinh = function(z) {
  var e = Jmat.Complex.exp(z);
  var ei = Jmat.Complex.inv(e);
  return e.sub(ei).divr(2);
};

Jmat.Complex.cosh = function(z) {
  var e = Jmat.Complex.exp(z);
  var ei = Jmat.Complex.inv(e);
  return e.add(ei).divr(2);
};

Jmat.Complex.tanh = function(z) {
  var e = Jmat.Complex.exp(z);
  var ei = Jmat.Complex.inv(e);
  return e.sub(ei).div(e.add(ei));
};

Jmat.Complex.asinh = function(z) {
  return Jmat.Complex.log(z.add(Jmat.Complex.sqrt(z.mul(z).addr(1))));
};

Jmat.Complex.acosh = function(z) {
  // ln(x + sqrt(z-1)*sqrt(z+1))
  return Jmat.Complex.log(z.add(Jmat.Complex.sqrt(z.subr(1)).mul(Jmat.Complex.sqrt(z.addr(1)))));
};

Jmat.Complex.atanh = function(z) {
  // 0.5 * (ln(1+z) - ln(1-z))
  return Jmat.Complex.log(z.inc().div(z.dec())).mulr(0.5);
};

// This is NOT the logsine function (the integral). It's simply ln(sin(z))
//ln(sin(z)), with good approximation for large |Im(z)|. The thing is, for large imaginary values, sin(z) becomes huge, because it involves an exponential of the imaginary parts
// For large imaginary part (or very small below 0), log(sin(x)) fails while this function is then very accurate
Jmat.Complex.logsin = function(z) {
  if(z.im > -10 && z.im < 10) return Jmat.Complex.log(Jmat.Complex.sin(z));

  var ln2i = Jmat.Complex(0.69314718056, 1.570796326795); // ln(2i)
  // This approximation is using a formula e^ix/2i or -e^(-ix)/2i, instead of the full (e^ix - e^(-ix) / 2i) = sin(x). This requires the real part to be exactly in range -pi/2, 3pi/2. So wrap, since it's periodic.
  var p = Jmat.Complex(Jmat.Real.wrap(z.re, -Math.PI / 2, 3 * Math.PI / 2), z.im);
  if(z.im > 0) return Jmat.Complex.newi(Jmat.Complex.PI).sub(Jmat.Complex.I.mul(p)).sub(ln2i);
  else return Jmat.Complex.I.mul(p).sub(ln2i);
};

// See description of Jmat.Complex.logsin
Jmat.Complex.logcos = function(z) {
  return Jmat.Complex.logsin(z.rsub(Math.PI / 2));
};

Jmat.Complex.floor = function(x) {
  var result = Jmat.Complex(0);
  result.re = Math.floor(x.re);
  result.im = Math.floor(x.im);
  return result;
};

Jmat.Complex.ceil = function(x) {
  var result = Jmat.Complex(0);
  result.re = Math.ceil(x.re);
  result.im = Math.ceil(x.im);
  return result;
};

Jmat.Complex.round = function(x) {
  var result = Jmat.Complex(0);
  result.re = Math.round(x.re);
  result.im = Math.round(x.im);
  return result;
};

// truncate towards 0
Jmat.Complex.trunc = function(x) {
  var result = Jmat.Complex(0);
  result.re = x.re < 0 ? Math.ceil(x.re) : Math.floor(x.re);
  result.im = x.im < 0 ? Math.ceil(x.im) : Math.floor(x.im);
  return result;
};

// Fractional part of x, x - floor(x). NOTE: this variant gives positive results for negative x
Jmat.Complex.frac = function(x) {
  return Jmat.Complex(Jmat.Real.frac(x.re), Jmat.Real.frac(x.im));
};

// Fractional part of x, x - int(x). NOTE: this variant gives negative results for negative x
Jmat.Complex.fracn = function(x) {
  return Jmat.Complex(Jmat.Real.fracn(x.re), Jmat.Real.fracn(x.im));
};

Jmat.Complex.exp = function(x) {
  if(x.im == 0) {
    return Jmat.Complex(Math.exp(x.re));
  } else {
    var ea = Math.exp(x.re);
    return new Jmat.Complex(ea * Math.cos(x.im), ea * Math.sin(x.im));
  }
};

//exp(x) - 1, with better precision for x around 0
Jmat.Complex.expm1 = function(x) {
  if(Jmat.Complex.abssqr(x) < 1e-5) return x.add(x.mul(x).divr(2)).add(x.mul(x).mul(x).divr(6));
  else return Jmat.Complex.exp(x).subr(1);
};

//natural log (base e, ln)
Jmat.Complex.log = function(x) {
  if(x.eqr(-Infinity)) return Jmat.Complex(Infinity);

  if(Jmat.Complex.isReal(x) && x.re >= 0) {
    return Jmat.Complex(Math.log(x.re));
  }

  return Jmat.Complex(Math.log(Jmat.Complex.absr(x)), Jmat.Complex.argr(x));
};

//ln(x + 1), with better precision for x around 0
Jmat.Complex.log1p = function(x) {
  if(Jmat.Complex.abssqr(x) < 1e-8) return x.mulr(-0.5).addr(1).mul(x);
  else return Jmat.Complex.log(x.addr(1));
};

//arbitrary log: log_y(x)
//warning: base y is second argument
Jmat.Complex.logy = function(x, y) {
  return Jmat.Complex.log(x).div(Jmat.Complex.log(y));
};

Jmat.Complex.sqrt = function(x) {
  if(Jmat.Complex.isReal(x)) {
    var result = Jmat.Complex(0);
    if(x.re >= 0 || x.re != x.re) result.re = Math.sqrt(x.re);
    else result.im = Math.sqrt(-x.re);
    return result;
  } else return x.pow(Jmat.Complex(0.5));
};

Jmat.Complex.toInt = function(value) {
  return Math.round(value.re);
};

// normalizes even if re or im are infinite, e.g. (Infinity, -Infinity) becomes (1, -1), (0, Infinity) becomes (0, 1). Without infinities, remains as-is. Des not normalize to length 1.
Jmat.Complex.infNormalize = function(value) {
  if (Jmat.Complex.isNaN(value)) return Jmat.Complex(NaN);

  if (value.re == Infinity) {
    if (value.im == Infinity) return Jmat.Complex(1, 1);
    if (value.im == -Infinity) return Jmat.Complex(1, -1);
    return Jmat.Complex(1, 0);
  }
  if (value.re == -Infinity) {
    if (value.im == Infinity) return Jmat.Complex(-1, 1);
    if (value.im == -Infinity) return Jmat.Complex(-1, -1);
    return Jmat.Complex(-1, 0);
  }
  if (value.im == Infinity) {
    if (value.re == Infinity) return Jmat.Complex(1, 1);
    if (value.re == -Infinity) return Jmat.Complex(-1, 1);
    return Jmat.Complex(0, 1);
  }
  if (value.im == -Infinity) {
    if (value.re == Infinity) return Jmat.Complex(1, -1);
    if (value.re == -Infinity) return Jmat.Complex(-1, -1);
    return Jmat.Complex(0, -1);
  }

  return value.divr(value.absr());
};

// Automatically cache last value. Useful for parameters of statistical distributions that are often the same in repeated calls.
// Cache must be an array (initially []), so that this function can modify it to set the necessary values.
// Function fun is called with z.
// n is cache size
// if n is given, cache contains alternating: index, input0, result0, input1, result1, input2, result2, ... where index is circular pointer to fill in new cache values
// if n is not given, cache contains: input, result
Jmat.Complex.calcCache_ = function(z, fun, cache, n) {
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
Jmat.Complex.gamma = function(z) {
  if(z.re == Infinity) return Jmat.Complex(Infinity);
  if(Jmat.Complex.isNegativeIntOrZero(z)) return Jmat.Complex(Infinity, Infinity); // Undirected infinity
  if(z.im == 0) return Jmat.Complex(Jmat.Real.gamma(z.re));

  // The internal function that doesn't do internal checks
  var gamma_ = function(z) {
    if(z.re < 0.5) {
      // Use the reflection formula, because, the approximation below is not accurate
      // for values around -6.5+0.1i
      // gamma(1-z)*gamma(z) = pi/sin(pi*z)
      var result = Jmat.Complex.PI.div(Jmat.Complex.sin(Jmat.Complex.PI.mul(z))).div(gamma_(Jmat.Complex.ONE.sub(z)));
      if(Jmat.Complex.isNaN(result)) result = Jmat.Complex(0); // For those values that it can't calculate, it's 0 on the negative side of the complex plane.
      return result;
    }

    var g = 7;
    var p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];

    z = z.subr(1);
    var x = Jmat.Complex(p[0]);
    for(var i = 1; i < g + 2; i++) {
      x = x.add(Jmat.Complex(p[i]).div(z.addr(i)));
    }
    var t = z.addr(g + 0.5);
    var pisq = Math.sqrt(Math.PI * 2);

    var w = t.pow(z.addr(0.5));
    var e = Jmat.Complex.exp(t.neg());
    
    var result = w.mul(e).mul(x).mulr(pisq);
    return result;
  };

  return gamma_(z);
};

Jmat.Complex.factorial = function(a) {
  return Jmat.Complex.gamma(Jmat.Complex.inc(a));
};

//using Stirling series
//logarithm of the gamma function, and more specific branch of the log
Jmat.Complex.loggamma = function(z) {
  //the result is way too imprecise if the real part of z is < 0, use the log of the reflection formula
  // loggamma(z) = log(pi/sin(pi*z)) - loggamma(1 - z)
  if(z.re < 0) {
    if(z.im == 0 && z.re == Math.floor(z.re)) return Jmat.Complex(NaN); // gamma does not exist here, so it shouldn't return bogus values
    var l = Jmat.Complex.log(Jmat.Complex.PI.div(Jmat.Complex.sin(Jmat.Complex.PI.mul(z))));
    // the complex sine goes out of bounds for example for (-4+120i), log(pi / sin(pi*z)) is very roughly approximated by log(2*pi*i) - i*z*pi*sign(z.im) (TODO: the approximation is not fully correct, e.g. for -160-116i)
    if(Jmat.Complex.isInfOrNaN(l)) l = Jmat.Complex.log(Jmat.Complex.newi(2 * Math.PI)).sub(Jmat.Complex.newi(-Math.PI).mul(z.im > 0 ? z : z.neg()));
    return l.sub(Jmat.Complex.loggamma(Jmat.Complex.ONE.sub(z)));
  }

  // The series below has a weird artefact for values near 0 with re > 0. Use actual log(gamma) for that
  if(z.im < 1 && z.im > -1 && z.re < 1 && z.re >= 0) return Jmat.Complex.log(Jmat.Complex.gamma(z));

  var result = Jmat.Complex(0.918938533205); //0.5 * ln(2pi)
  result = result.add(Jmat.Complex.subr(z, 0.5).mul(Jmat.Complex.log(z)));
  result = result.sub(z);
  result = result.add(z.mulr(12).inv());
  result = result.sub(Jmat.Complex.powr(z, 3).mulr(360).inv());
  result = result.add(Jmat.Complex.powr(z, 5).mulr(1260).inv());
  result = result.sub(Jmat.Complex.powr(z, 7).mulr(1680).inv());
  result = result.add(Jmat.Complex.powr(z, 9).mulr(1188).inv());
  return result;
};

//Inverse of the gamma function (not the reciproke, the inverse or "arc" function)
//source: http://mathforum.org/kb/message.jspa?messageID=342551&tstart=0
//lx = log(x + c) / sqrt(2pi); result = lx / lambertw(lx / e) + 0.5
//not very precise
Jmat.Complex.gamma_inv = function(value) {
  if(Jmat.Complex.isPositive(value) && x > 0.85) { //doesn't work for negative values, nor values smaller than somewhere around 0.85
    // Approximation for positive real x.
    var x = value.re;
    //c = sqrt(2 * pi) / e - gamma(k), where k = the positive zero of the digamma function (1.4616321449683623412626...)
    var c = 0.03653381448490041660;
    var lx = Math.log((x + c) / Math.sqrt(2 * Math.PI));
    return Jmat.Complex(lx / Jmat.Real.lambertw(lx / Math.E) + 0.5);
  }

  // TODO: this has problems for |z| > 20. Use more stable root finding
  // TODO: since the current digamma implementation is nothing more than a derivative approximation, I might as well use finvert_secant. Try using digamma anyway if it's ever made more precise.
  var result = Jmat.Complex.finvert_newton(value, Jmat.Complex.gamma, function(z) {
    // Derivative of gamma function is: gamma function multiplied by digamma function
    return Jmat.Complex.gamma(z).mul(Jmat.Complex.digamma(z));
  });
  if(!Jmat.Complex.near(Jmat.Complex.gamma(result), value, 0.01)) return Jmat.Complex(NaN);
  return result;
};

// digamma function: psi_0(z)
Jmat.Complex.digamma = function(z) {
  // digamma(z) = gamma'(z) / gamma(z)
  // the derivative gamma'(z) is approximated here as (gamma(z+epsilon) - gamma(z-epsilon)) / (2*epsilon). TODO: use better approximation
  // TODO: for real z, use an approximation such as the Euler Maclaurin formula which gives ln(x) -1/2 x - 1/12 xx + 1/120 xxxx + ....

  // METHOD A: simple approximate derivative of gamma function - works in fact quite well
  // There are two ways: the derivative of loggamma, or, the derivative of gamma divided through gamma. The first is faster [TODO: verify that], but does not work near any z.re that has fractional part 0 or 0.5 due to branch cuts of loggamma.
  var f = Jmat.Real.frac(z.re);
  if(f < 0.001 || f > 0.999 || (f > 0.499 && f < 0.501)) {
    // Too small for even this formula: everything becomes numerically 0, Approximate with nearby value instead.
    if(z.re < -98) return Jmat.Complex.loggamma(z.addr(0.0003)).sub(Jmat.Complex.loggamma(z.addr(0.0001))).divr(0.0002);
    // Near integer or half integer real part, 3 gammas required
    var result = Jmat.Complex.gamma(z.addr(0.0001)).sub(Jmat.Complex.gamma(z.subr(0.0001))).divr(0.0002).div(Jmat.Complex.gamma(z));
    // Again, sometimes 0/0 happens simply because the numbers are underflowing. Approximate with nearby value.
    if(Jmat.Complex.isInfOrNaN(result) && !Jmat.Complex.isInf(z) && z.im != 0) return Jmat.Complex.loggamma(z.addr(0.0003)).sub(Jmat.Complex.loggamma(z.addr(0.0001))).divr(0.0002);
    return result;
  }
  return Jmat.Complex.loggamma(z.addr(0.0001)).sub(Jmat.Complex.loggamma(z.subr(0.0001))).divr(0.0002);

  // METHOD B: series - does not work well, requires thousands of iterations and even then the approximate derivative works better
  /*var sum = Jmat.Complex.ZERO;
  for(var n = 0; n < 300; n++) {
    var d = z.addr(n).mulr(n + 1);
    sum = sum.add(z.dec().div(d));
  }
  // subtract Euler-Mascheroni constant
  return sum.subr(Jmat.Real.EM);*/
};

Jmat.Complex.trigamma = function(z) {
  // The following is implemented in here, but then in such way to avoid redundant identical gamma calls
  //return Jmat.Complex.digamma(z.addr(0.0001)).sub(Jmat.Complex.digamma(z.subr(0.0001))).divr(0.0002);

  // There are two ways: the derivative of loggamma, or, the derivative of gamma divided through gamma. The first is faster [TODO: verify that], but does not work near any z.re that has fractional part 0 or 0.5 due to branch cuts of loggamma.
  var f = Jmat.Real.frac(z.re);
  if(f < 0.001 || f > 0.999 || (f > 0.499 && f < 0.501)) {
    var a = Jmat.Complex.gamma(z.addr(-0.0002));
    var b = Jmat.Complex.gamma(z.addr(-0.0001));
    var c = Jmat.Complex.gamma(z);
    var d = Jmat.Complex.gamma(z.addr(0.0001));
    var e = Jmat.Complex.gamma(z.addr(0.0002));

    var d0 = a.sub(c).divr(0.0002).div(b);
    var d1 = c.sub(e).divr(0.0002).div(d);

    return d0.sub(d1).divr(0.0002);
  } else {
    var a = Jmat.Complex.loggamma(z.addr(-0.0002));
    var b = Jmat.Complex.loggamma(z);
    var c = Jmat.Complex.loggamma(z.addr(0.0002));

    var d0 = a.sub(b).divr(0.0002);
    var d1 = b.sub(c).divr(0.0002);

    return d0.sub(d1).divr(0.0002);
  }
};

Jmat.Complex.polygamma = function(n, z) {
  if(n.eqr(0)) return Jmat.Complex.digamma(z);
  if(n.eqr(1)) return Jmat.Complex.trigamma(z);

  // METHOD A: series: Only for positive integer n. Does not work well, requires too many iterations
  /*if(Jmat.Complex.isPositiveInt(n))
    var sum = Jmat.Complex.ZERO;
    for(var k = 0; k < 300; k++) {
      var d = z.addr(k).pow(n.inc());
      sum = sum.add(d.inv());
    }
    return Jmat.Complex(-1).pow(n.inc()).mul(Jmat.Complex.factorial(n)).mul(sum);
  }*/

  // METHOD B: with hurwitz zeta
  var h = Jmat.Complex.hurwitzzeta(n.inc(), z);
  return Jmat.Complex(-1).pow(n.inc()).mul(Jmat.Complex.factorial(n)).mul(h);
};

// lower incomplete gamma function
// lowercase gamma(s, x) 
Jmat.Complex.incgamma_lower = function(s, z) {
  // METHOD A: in terms of hypergeometric1F1 function
  //return z.pow(s).div(s).mul(Jmat.Complex.hypergeometric1F1(s, s.inc(), z.neg()));

  // METHOD B: series expansion - has some problems with division through zero
  // sum_k=0.oo ((-1)^k / k!) * (z^(s+k) / (s + k))
  var result = Jmat.Complex(0);
  var kk = Jmat.Complex(1);
  var zz = z.pow(s);
  var sign = 1;
  for(var k = 0; k < 30; k++) {
    if(k > 0) {
      kk = kk / k;
      sign = -sign;
      zz = zz.mul(z);
    }
    result = result.add(Jmat.Complex(sign * kk).mul(zz).div(s.addr(k)));
  }
  return result;
};

// upper incomplete gamma function
// uppercase GAMMA(s, x) 
Jmat.Complex.incgamma_upper = function(s, z) {
  // For negative integer s, gamma, and lower incomplete gamma, are not defined. But the upper is.
  if(Jmat.Complex.isNegativeIntOrZero(s)) {
    s = s.addr(1e-7); //twiddle it a bit for negative integers, so that formula in terms of lower gamma sort of works... TODO: use better approximation
  }

  return Jmat.Complex.gamma(s).sub(Jmat.Complex.incgamma_lower(s, z));
};

Jmat.Complex.gamma_p_cache_ = []; // cache used because s is often constant between gamma_p calls

// regularized lower incomplete gamma function (lower gamma_regularized, regularized_gamma)
// P(s, x) = gamma(s, z) / GAMMA(s)
// Note: the derivative of this function is: (e^(-z) * z^(s-1))/GAMMA(s)
Jmat.Complex.gamma_p = function(s, z) {
  if(Jmat.Complex.isNegativeIntOrZero(s)) return Jmat.Complex(1);
  var g = Jmat.Complex.calcCache_(s, Jmat.Complex.gamma, Jmat.Complex.gamma_p_cache_); // gamma(s)
  return Jmat.Complex.incgamma_lower(s, z).div(g);
};

// regularized upper incomplete gamma function (upper gamma_regularized, regularized_gamma)
// Q(s, x) = 1 - P(s, x) = GAMMA(s, z) / GAMMA(s)
// Note: the derivative of this function is: -(e^(-z) * z^(s-1))/GAMMA(s)
Jmat.Complex.gamma_q = function(s, z) {
  if(Jmat.Complex.isNegativeIntOrZero(s)) return Jmat.Complex(0);
  return Jmat.Complex.ONE.sub(Jmat.Complex.gamma_p(s, z));
};

// One possible approximation series for inverse gamma P. Valid for real values with p in range 0-1 and positive s. Possibly more.
Jmat.Complex.gamma_p_inv_series_1_ = function(s, p) {
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

  var r = p.mul(Jmat.Complex.gamma(s1)).pow(s.inv());

  return Jmat.Complex.powerSeries(c, c.length, Jmat.Complex.ZERO, r);
};

// Inverse regularized gamma P
// Finds z for p == gamma_p(s, z)
// This is an *approximation* of inverse of incomplete regularized gamma (P) - it works for real values with p in range 0-1 and positive s. Good enough for qf of chi square distribution
Jmat.Complex.gamma_p_inv = function(s, p) {
  // Return NaN for unsupported values to prevent bogus results
  if(!Jmat.Complex.isReal(p) || !Jmat.Complex.isReal(s) || p.re < 0 || p.re > 1 || s.re < 0) return Jmat.Complex(NaN);

  // TODO: more complete support of the entire complex domain
  return Jmat.Complex.gamma_p_inv_series_1_(s, p);
};

// Calculates gamma(x) / gamma(y), and can cancel out negative integer arguments if both are negative integers in some cases. That is useful for functions like beta, binomial, ...
// Uses the formula gamma(a)/gamma(b) = (-1)^(a-b) * gamma(-b+1) / gamma(-a+1) if necessary
// It is also optimized to use for loops if x and y are nearby integers, ...
Jmat.Complex.gammaDiv_ = function(x, y) {
  if(x.eq(y)) return Jmat.Complex.ONE; // For the "combined" function, this is considered correct even for negative integers...

  if(Jmat.Complex.isInfOrNaN(x) || Jmat.Complex.isInfOrNaN(y)) {
    return Jmat.Complex(NaN);
  }

  if(Jmat.Complex.isNegativeIntOrZero(y) && (x.re > 0 || !Jmat.Complex.isInt(x))) return Jmat.Complex.ZERO; // division of non-infinity through infinity

  if(Jmat.Complex.isInt(x) && Jmat.Complex.isInt(y)) {
    if(x.re <= 0 && y.re <= 0) {
      var sign = Jmat.Real.isOdd(x.re - y.re) ? -1 : 1;
      return Jmat.Complex.gammaDiv_(y.neg().addr(1), x.neg().addr(1)).mulr(sign);
    }

    if(x.re > 0 && y.re > 0 && Jmat.Real.dist(x.re, y.re) < 16) {
      if(x.re > y.re) {
        var result = y.re;
        for(var z = y.re + 1; z < x.re; z++) result *= z;
        return Jmat.Complex(result);
      } else {
        var result = 1 / x.re;
        for(var z = x.re + 1; z < y.re; z++) result /= z;
        return Jmat.Complex(result);
      }
    }
  }

  return Jmat.Complex.gamma(x).div(Jmat.Complex.gamma(y));
};

// Similar to Jmat.Complex.gammaDiv_, but takes 3 arguments and cancels out more if possible. Returns gamma(a) / (gamma(b) * gamma(c))
Jmat.Complex.gammaDiv12_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(Jmat.Complex.isNegativeIntOrZero(b)) {
    if(Jmat.Complex.isNegativeIntOrZero(c) && (a.re > 0 || !Jmat.Complex.isInt(a))) return Jmat.Complex.ZERO; // division of non-infinity through infinity
    return Jmat.Complex.gammaDiv_(a, b).div(Jmat.Complex.gamma(c));
  } else {
    return Jmat.Complex.gammaDiv_(a, c).div(Jmat.Complex.gamma(b));
  }
};

// Similar to Jmat.Complex.gammaDiv_, but takes 3 arguments and cancels out more if possible. Returns (gamma(a) * gamma(b)) / gamma(c)
Jmat.Complex.gammaDiv21_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(Jmat.Complex.isNegativeIntOrZero(a)) {
    return Jmat.Complex.gammaDiv_(a, c).mul(Jmat.Complex.gamma(b));
  } else {
    return Jmat.Complex.gammaDiv_(b, c).mul(Jmat.Complex.gamma(a));
  }
};

// Similar to Jmat.Complex.gammaDiv_, but takes 4 arguments and cancels out more if possible. Returns (gamma(a) * gamma(b)) / (gamma(c) * gamma(d))
Jmat.Complex.gammaDiv22_ = function(a, b, c, d) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(Jmat.Complex.isNegativeIntOrZero(a) == Jmat.Complex.isNegativeIntOrZero(c)) {
    return Jmat.Complex.gammaDiv_(a, c).mul(Jmat.Complex.gammaDiv_(b, d));
  } else {
    return Jmat.Complex.gammaDiv_(a, d).mul(Jmat.Complex.gammaDiv_(b, c));
  }
};

// Calculates log(gamma(x) / gamma(y)) = loggamma(x) - loggamma(y), and can cancel out negative integer arguments if both are negative integers in some cases. That is useful for functions like beta, binomial, ...
// Uses the formula gamma(a)/gamma(b) = (-1)^(a-b) * gamma(-b+1) / gamma(-a+1) if necessary
Jmat.Complex.loggammaDiv_ = function(x, y) {
  if(x.eq(y)) return Jmat.Complex.ZERO; // For the "combined" function, this is correct even for negative integers...

  if(Jmat.Complex.isInfOrNaN(x) || Jmat.Complex.isInfOrNaN(y)) {
    return Jmat.Complex(NaN);
  }

  if(Jmat.Complex.isNegativeIntOrZero(y) && (x.re > 0 || !Jmat.Complex.isInt(x))) return Jmat.Complex(-Infinity); // division of non-infinity through infinity ==> log(0)

  if(Jmat.Complex.isInt(x) && Jmat.Complex.isInt(y)) {
    if(x.re <= 0 && y.re <= 0) {
      var sign = Jmat.Real.isOdd(x.re - y.re) ? -1 : 1;
      var result = Jmat.Complex.loggammaDiv_(y.neg().addr(1), x.neg().addr(1));
      if(sign == -1) result = result.add(Jmat.Complex.newi(Math.PI)); // log(-x) = log(x) + i*pi
      return result;
    }
  }

  return Jmat.Complex.loggamma(x).sub(Jmat.Complex.loggamma(y));
};

// Similar to Jmat.Complex.loggammaDiv_, but takes 3 arguments and cancels out more if possible. Returns loggamma(a) - (loggamma(b) + loggamma(c))
Jmat.Complex.loggammaDiv12_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(Jmat.Complex.isNegativeIntOrZero(b)) {
    if(Jmat.Complex.isNegativeIntOrZero(c) && (a.re > 0 || !Jmat.Complex.isInt(a))) return Jmat.Complex(-Infinity); // division of non-infinity through infinity ==> log(0)
    return Jmat.Complex.loggammaDiv_(a, b).sub(Jmat.Complex.loggamma(c));
  } else {
    return Jmat.Complex.loggammaDiv_(a, c).sub(Jmat.Complex.loggamma(b));
  }
};

// Similar to Jmat.Complex.loggammaDiv_, but takes 3 arguments and cancels out more if possible. Returns (loggamma(a) + loggamma(b)) - loggamma(c)
Jmat.Complex.loggammaDiv21_ = function(a, b, c) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(Jmat.Complex.isNegativeIntOrZero(a)) {
    return Jmat.Complex.loggammaDiv_(a, c).add(Jmat.Complex.loggamma(b));
  } else {
    return Jmat.Complex.loggammaDiv_(b, c).add(Jmat.Complex.loggamma(a));
  }
};

// Similar to Jmat.Complex.loggammaDiv_, but takes 4 arguments and cancels out more if possible. Returns (loggamma(a) + loggamma(b)) - (loggamma(c) + loggamma(d))
// To have only 3 values, set a, b, c or d to 1 (it will be fast for that)
Jmat.Complex.loggammaDiv2_ = function(a, b, c, d) {
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(Jmat.Complex.isNegativeIntOrZero(a) == Jmat.Complex.isNegativeIntOrZero(c)) {
    return Jmat.Complex.loggammaDiv_(a, c).add(Jmat.Complex.loggammaDiv_(b, d));
  } else {
    return Jmat.Complex.loggammaDiv_(a, d).add(Jmat.Complex.loggammaDiv_(b, c));
  }
};

// Returns a particular non-principal complex branch of sqrt(a * b)
// For AGM and GHM, a particular branch of the complex sqrt must be returned, otherwise it hops to other branches and gives wrong result
Jmat.Complex.agmMulSqrt_ = function(a, b) {
  // Note: if a and b are real but negative, with a*b being positive, still a different branch is used.
  if(Jmat.Complex.isPositive(a) && Jmat.Complex.isPositive(b)) {
    return Jmat.Complex.sqrt(a.mul(b));
  } else {
    return Jmat.Complex.sqrt(Jmat.Complex.abs(a.mul(b))).mul(Jmat.Complex.exp(Jmat.Complex.I.mul(Jmat.Complex.arg(a).add(Jmat.Complex.arg(b)).divr(2))));
  }
};

//Arithmetic-Geometric mean (iteratively calculated)
Jmat.Complex.agm = function(a, b) {
  // Wikipedia: We have the following inequality involving the Pythagorean means {H, G, A} and iterated Pythagorean means {HG, HA, GA}:
  // min <= H <= HG <= (G == HA) <= GA <= A <= max
  // For completeness, also (with L = logarithmic mean, R = RMS, C = contraharmonic mean):
  // min <= H <= G <= L <= A <= R <= C <= max

  // The AGM is normally only defined for positive real numbers, but it can be extended to whole complex plane, as is done here.
  
  //avoid imprecisions for special cases
  if(a.eq(b.neg()) || a.eq(Jmat.Complex.ZERO) || b.eq(Jmat.Complex.ZERO)) return Jmat.Complex(0);
  var real = Jmat.Complex.isReal(a) && Jmat.Complex.isReal(b) && (a.re < 0) == (b.re < 0);

  var a2, b2;
  for(var i = 0; i < 60; i++) {
    if(a.eq(b)) break;
    a2 = a.add(b).divr(2);
    b2 = Jmat.Complex.agmMulSqrt_(a, b);
    a = a2;
    b = b2;
  }

  if(real) a.im = 0; // it may be 1e-16 or so due to numerical error (only if both inputs are real and have same sign)
  
  return a;
};

//Geometric-Harmonic mean (iteratively calculated)
Jmat.Complex.ghm = function(a, b) {
  // Not really defined for negative and complex numbers, but when using Jmat.Complex.agmMulSqrt_ it looks smooth in the 2D plot
  // NOTE: An alternative, that returns different values for neg/complex (but same for positive reals) is: return Jmat.Complex.agm(a.inv(), b.inv()).inv();
  var a2, b2;
  for(var i = 0; i < 60; i++) {
    if(a.eq(b)) break;
    a2 = Jmat.Complex.agmMulSqrt_(a, b);
    b2 = Jmat.Complex.TWO.div(a.inv().add(b.inv()));
    a = a2;
    b = b2;
  }
  return a;
};

// Bessel function of the first kind
// Mostly an approximation, there are problems with high z
// TODO: make faster and fix problems
Jmat.Complex.besselj = function(n, z) {
  var pi = Math.PI;

  if(z.im == 0 && n.im == 0 && z.re < 0 && z.re * z.re < (n.re + 1) / 10) {
    return Jmat.Complex.inv(Jmat.Complex.gamma(n.inc())).mul(z.divr(2).pow(n));
  } else if(z.absr() < 20) {
    // The gamma functions give NaN if n is negative integer < -1. TODO: also for n = 0 and neg won't help there, fix that
    var negintn = Jmat.Complex.isNegativeInt(n);
    if(negintn) n = n.neg();
    var result = Jmat.Complex(0);
    var m = 1;
    var gn = Jmat.Complex.gamma(n.inc()); // during the loop: gamma(n + i + 1)
    for(var i = 0; i < 50; i++) {
      var i1 = Jmat.Complex(i + 1);
      var d = Jmat.Complex.gamma(i1).mul(gn); //calling gamma(i1) is no problem, all integer solutions in this range are cached
      var term = Jmat.Complex(m).div(d).mul(z.divr(2).pow(Jmat.Complex(2 * i).add(n)));
      m = -m;
      result = result.add(term);
      gn = gn.mul(n.add(i1));
    }
    if(negintn && Jmat.Real.isOdd(n.re)) result = result.neg(); //besselj(-n, x) = (-1)^n * besselj(n, x)
    return result;
  } else if(z.im == 0) {
    // Something is wrong with this formula, see the 2D plot
    var a = Jmat.Complex.sqrt(Jmat.Complex(2/pi).div(z));
    var b = z.sub(n.mulr(pi/2)).subr(pi/4);
    return a.mul(Jmat.Complex.cos(b));
  } else {
    // Something is wrong with this formula, see the 2D plot
    var s = z.im > 0 ? -1 : 1
    var a = z.sub(n.mulr(pi/2)).subr(pi/4);
    var b = Jmat.Complex.sqrt(z.mulr(pi*2));
    return Jmat.Complex.exp(Jmat.Complex.newi(s).mul(a)).div(b);
  }

  // METHOD B: in terms of confluent hypergeometric
  /*var neg = 1;
  if(Jmat.Complex.isNegativeInt(n)) {
    if(Jmat.Complex.isOdd(n)) neg = -1;
    n = n.neg(); // Formula does not work for negative n, but J_-n(z) = (-1)^n*J_n(z) for integer n
  }
  var a = z.mulr(0.5).pow(n).div(Jmat.Complex.gamma(n.inc()));
  var b = Jmat.Complex.hypergeometric0F1(n.inc(), z.mul(z).mulr(-0.25));
  return a.mul(b).mulr(neg);*/
};

// Bessel function of the second kind
// Mostly an approximation, there are problems with high z
// TODO: make faster and fix problems
Jmat.Complex.bessely = function(n, z) {
  var pi = Math.PI;

  if(Jmat.Complex.abs(z).re < 15) {
    // Y_a(x) = (J_a(x)*cos(a*pi) - J_-a(x)) / sin(a*pi)
    // For integer n, Y_n(x) = lim_a_to_n(Y_a(x)
    if(n.re == Math.floor(n.re)) n = n.addr(0.000001); // TODO: really fix this, find better formula for bessel Y of integer order

    var a = Jmat.Complex.besselj(n, z);
    var b = Jmat.Complex.cos(n.mulr(Math.PI));
    var c = Jmat.Complex.besselj(n.neg(), z);
    var d = Jmat.Complex.sin(n.mulr(Math.PI));
    return a.mul(b).sub(c).div(d);
  } else if(z.im == 0) {
    // Something is wrong with this formula, see the 2D plot
    var a = Jmat.Complex.sqrt(Jmat.Complex(2/pi).div(z));
    var b = z.sub(n.mulr(pi/2)).subr(pi/4);
    return a.mul(Jmat.Complex.sin(b));
  } else {
    // Something is wrong with this formula, see the 2D plot
    var s = z.im > 0 ? -1 : 1;
    var a = z.sub(n.mulr(pi/2)).subr(pi/4);
    var b = Jmat.Complex.sqrt(z.mulr(pi*2));
    return Jmat.Complex.newi(-s).mul(Jmat.Complex.exp(Jmat.Complex.newi(s).mul(a)).div(b));
  }
};

// Hankel function of the first kind (approximation)
// TODO: this one is very slow. Make faster and more precise
Jmat.Complex.hankel1 = function(n, z) {
  // There are numerical imprecisions for e.g. hankel2(-8-8i, i).
  // So when needed, apply the formula: H1_-a(x) = exp(a*pi*i)*H1_a(x)
  if(n.im < 0) return Jmat.Complex.exp(n.mul(Jmat.Complex.newi(-Math.PI))).mul(Jmat.Complex.hankel1(n.neg(), z));

  return Jmat.Complex.besselj(n, z).add(Jmat.Complex.bessely(n, z).mul(Jmat.Complex.I));
};

// Hankel function of the second kind (approximation)
// TODO: this one is very slow. Make faster and more precise
Jmat.Complex.hankel2 = function(n, z) {
  // There are numerical imprecisions for e.g. hankel2(-8-8i, i).
  // So when needed, apply the formula: H2_-a(x) = exp(-a*pi*i)*H2_a(x)
  if(n.im > 0) return Jmat.Complex.exp(n.mul(Jmat.Complex.newi(Math.PI))).mul(Jmat.Complex.hankel2(n.neg(), z));

  return Jmat.Complex.besselj(n, z).sub(Jmat.Complex.bessely(n, z).mul(Jmat.Complex.I));
};

// Modified bessel function of the first kind (approximation)
// TODO: this one is very slow. Make faster and more precise
Jmat.Complex.besseli = function(n, z) {
  var result = Jmat.Complex.I.pow(n.neg()).mul(Jmat.Complex.besselj(n, z.mul(Jmat.Complex.I)));
  if(z.im == 0 && n.im == 0 && Jmat.Real.near(result.im, 0, 1e-10)) result.im = 0;
  return result;
};

// Modified bessel function of the second kind (approximation)
// TODO: this one is very slow. Make faster and more precise
Jmat.Complex.besselk = function(n, z) {
  var result;
  if(z.im >= 0) {
    result = Jmat.Complex.I.pow(n.inc()).mulr(Math.PI / 2).mul(Jmat.Complex.hankel1(n, Jmat.Complex.I.mul(z)));
  } else {
    result = Jmat.Complex.newi(-1).pow(n.inc()).mulr(Math.PI / 2).mul(Jmat.Complex.hankel2(n, Jmat.Complex.newi(-1).mul(z)));
  }
  if(z.im == 0 && n.im == 0 && Jmat.Real.near(result.im, 0, 1e-3 /*this one is super imprecise...*/)) result.im = 0;
  return result;
};

//pl = 3^(-2/3) for airy, 3^(-4/3) for bairy
//pr = 3^(-4/3) for airy, 3^(-5/6) for bairy
//s = -1 for airy, +1 for bairy
Jmat.Complex.airyloop_ = function(z, pl, pr, s) {
  var gl = 1.3541179394264004169; // GAMMA(2/3)
  var gr = 0.8929795115692492112; // GAMMA(4/3)
  var zzz = z.mul(z).mul(z);
  var zr = Jmat.Complex.ONE;
  var zl = z;
  var r = Jmat.Complex.ONE;
  var kk = 1;
  var result = Jmat.Complex.ZERO;
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
    if(Jmat.Complex.isNaN(rl) || Jmat.Complex.isNaN(rr)) break;
    if(rl.eqr(0) && rr.eqr(0)) break;
    result = result.add(rl).add(rr.mulr(s));
  }

  return result;
};

// Airy function Ai(x)
Jmat.Complex.airy = function(z) {
  if(z.absr() > 8) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < 2 * Math.PI / 3) {
      var d = z.powr(1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return Jmat.Complex.exp(zeta.neg()).div(d.mulr(2));
    } else {
      var zm = z.neg();
      var d = zm.powr(1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return Jmat.Complex.sin(zeta.addr(Math.PI / 4)).div(d);
    }
  }

  // METHOD A: definition with hypergeometric0F1
  /*var zzz = z.mul(z).mul(z);
  var a = Jmat.Complex.hypergeometric0F1(Jmat.Complex(2/3), zzz.mulr(1/9)).mulr(Math.pow(3,-2/3) / Jmat.Real.gamma(2/3));
  var b = Jmat.Complex.hypergeometric0F1(Jmat.Complex(4/3), zzz.mulr(1/9)).mul(z).mulr(Math.pow(3,-1/3) / Jmat.Real.gamma(1/3));
  return a.sub(b);*/

  // METHOD B: summation:
  // SUM_k=0..oo 3^(-2k-2/3)/(k!*GAMMA(k+2/3))*x^(3k) - SUM_k=0..oo 3^(-2k-4/3)/(k!*GAMMA(k+4/3))*x^(3k+1)
  // This is basically the same as the hypergeometric definitions, but faster
  var pl = Math.pow(3, -2/3);
  var pr = Math.pow(3, -4/3);
  return Jmat.Complex.airyloop_(z, pl, pr, -1);
};

// Airy Bi function Bi(x)
Jmat.Complex.bairy = function(z) {
  if(z.absr() > 10) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < Math.PI / 3) {
      var d = z.powr(1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return Jmat.Complex.exp(zeta).div(d);
    } else {
      var zm = z.neg();
      var d = zm.powr(1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return Jmat.Complex.cos(zeta.addr(Math.PI / 4)).div(d);
    }
  }

  // SUM_k=0..oo 3^(-2k-1/6)/(k!*GAMMA(k+2/3))*x^(3k) + SUM_k=0..oo 3^(-2k-5/6)/(k!*GAMMA(k+4/3))*x^(3k+1)
  var pl = Math.pow(3, -1/6);
  var pr = Math.pow(3, -5/6);
  return Jmat.Complex.airyloop_(z, pl, pr, +1);
};

//pl = 3^(-1/3) for airy, 3^(+1/6) for bairy
//pr = 3^(-5/3) for airy, 3^(-7/6) for bairy
//s = -1 for airy, +1 for bairy
Jmat.Complex.airy_deriv_loop_ = function(z, pl, pr, s) {
  var gl = 2.6789385347077476336; // GAMMA(1/3)
  var gr = 0.9027452929509336112; // GAMMA(5/3)
  var zzz = z.mul(z).mul(z);
  var zr = Jmat.Complex.ONE;
  var zl = z.mul(z);
  var r = Jmat.Complex.ONE;
  var kk = 1;
  var result = Jmat.Complex.ZERO;
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
    if(Jmat.Complex.isNaN(rl) || Jmat.Complex.isNaN(rr)) break;
    if(rl.eqr(0) && rr.eqr(0)) break;
    result = result.add(rl.mulr(s)).add(rr);
  }

  return result;
};

// Derivative Airy function Ai'(x)
Jmat.Complex.airy_deriv = function(z) {
  if(z.absr() > 8) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < 2 * Math.PI / 3) {
      var d = z.powr(-1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return Jmat.Complex.exp(zeta.neg()).div(d.mulr(2)).neg();
    } else {
      var zm = z.neg();
      var d = zm.powr(-1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return Jmat.Complex.cos(zeta.addr(Math.PI / 4)).div(d).neg();
    }
  }

  // -SUM_k=0..oo 3^(-2k-1/3)/(k!*GAMMA(k+1/3))*x^(3k) + SUM_k=0..oo 3^(-2k-5/3)/(k!*GAMMA(k+5/3))*x^(3k+2)
  var pl = Math.pow(3, -1/3);
  var pr = Math.pow(3, -5/3);
  return Jmat.Complex.airy_deriv_loop_(z, pl, pr, -1);
};

// Derivative Airy Bi function Bi'(x)
Jmat.Complex.bairy_deriv = function(z) {
  if(z.absr() > 10) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.argr()) < Math.PI / 3) {
      var d = z.powr(-1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = z.powr(3/2).mulr(2/3);
      return Jmat.Complex.exp(zeta).div(d);
    } else {
      var zm = z.neg();
      var d = zm.powr(-1/4).mul(Jmat.Complex.SQRTPI);
      var zeta = zm.powr(3/2).mulr(2/3);
      return Jmat.Complex.sin(zeta.addr(Math.PI / 4)).div(d);
    }
  }

  // SUM_k=0..oo 3^(-2k+1/6)/(k!*GAMMA(k+1/3))*x^(3k) + SUM_k=0..oo 3^(-2k-7/6)/(k!*GAMMA(k+5/3))*x^(3k+2)
  var pl = Math.pow(3, +1/6);
  var pr = Math.pow(3, -7/6);
  return Jmat.Complex.airy_deriv_loop_(z, pl, pr, +1);
};


// This integral representation of riemann zeta works quite well and is valid for all s. This is quite slow though. And in practice it does not work for z.im < -42 or z.im > 42.
// That because the integrand fluctuates heavily between negative and positive values I think.
// Anyway, it is better than the eta series formula for the very particular values of x.re = 1 and x.im = a multiple of 9.0625, so used for that...
// The formula: zeta(s) = 2 * integrate_0..oo sin(s.atan(t)) / ((1+t^2)^(s/2) * (e^(2*pi*t) - 1) dt + 1/2 + 1/(s-1)
Jmat.Complex.zetaint_ = function(s) {
  var it = function(s, t) {
    var n = Jmat.Complex.sin(s.mul(Jmat.Complex.atan(t)));
    var ta = (t.mul(t).inc()).pow(s.divr(2));
    var tb = Jmat.Complex.exp(t.mulr(2*Math.PI)).dec();
    var t = ta.mul(tb);
    if(n.eq(t)) return Jmat.Complex.ONE; //avoid 0/0
    return n.div(t);
  };

  // Binary search to the place where the function starts becoming zero forever
  var q = Jmat.Complex.ONE;
  for(;;) {
    var q1 = it(s, q);
    var q2 = it(s, q.inc());
    if(q1.abssqr() < 1e-4 && q1.abssqr() < 1e-4) break;
    if(q.re > 1000) break; // do not do too many doublings or we are doing as many calculations in here as in the integral...
  
    q = q.mulr(2);
  }

  var maxt = q.re;
  var step = 0.01;

  var r = Jmat.Complex.ZERO;
  for(var t = 0; t < maxt; t+=step) {
    var a = it(s, Jmat.Complex(t)).mulr(step);
    r = r.add(a);
    step *= 1.01;
  }

  return r.mulr(2).addr(0.5).add(s.dec().inv());
};

// Riemann zeta function (infinite sum of reciprocal powers)
Jmat.Complex.zeta = function(s) {
  if(s.eq(Jmat.Complex.ZERO)) return Jmat.Complex(-0.5);
  if(s.eqr(-1)) return Jmat.Complex(-1/12);
  if(s.eqr(1)) return Jmat.Complex(+Infinity);
  if(s.eqr(3)) return Jmat.Complex(1.202056903159594); // Apery's constant
  if(s.re == +Infinity) return Jmat.Complex(1);
  if(Jmat.Complex.isEven(s) && s.re < 0) return Jmat.Complex(0); //it's 0 at all negative even integers

  // TODO: Use better algorithms:
  // Borwein algorithm for s close to real line
  // Riemann-Siegel formula for s with large imaginary party
  // Euler-Maclaurin summation in all other cases
  // The reflection formula where applicable

  var a1 = s.dec().absr();
  if(a1 < 10) {
    // Laurent series around 1: 1/(s-1) + SUM_n=0..oo (-1)^n / n! * stieltjes_n * (s-1)^n
    // For 30 terms, seems correct in a radius of around 10-15 around the point 1
    var s1 = s.dec();
    var result = s1.inv();
    var ss = Jmat.Complex.ONE;
    var num = a1 < 5 ? 10 : 30;
    for(var i = 0; i < num; i++) {
      result = result.add(Jmat.Complex(Jmat.Complex.stieltjes_zeta[i]).mul(ss));
      ss = ss.mul(s1);
    }
    return result;
  } else if(s.re >= 5) {
    // The riemann zeta infinite series: SUM_n=1..oo n^(-s) = 1/1^s + 1/2^s + 1/3^s + ...
    // Converges for Re(s) > 1, but in practice only converges properly for bigger Re(s)
    var result = Jmat.Complex(0);
    for(var i = 1; i < 32; i++) {
      result = result.add(Jmat.Complex(i).pow(s).inv());
    }
    return result;
  } else if(s.re >= 0.5 /*the articles make it seem as if this loop is for s.re > 0, but it only converges well for s.re 0.5*/) {
    // This is a series that is convergent (but not absolutely) for s.re > 0 instead of s.re > 1. In practice, it does not converge well, and only somewhat for s.re > 0.5, and even there not good...
    // eta(s) = SUM_n=1..oo (-1)^(n+1) / n^s = (1-2^(1-s)) * zeta(s) ==> zeta(s) = 1/(1-2^(1-s)) * (SUM_n=1..oo (-1)^(n+1) / n^s)
    if(Jmat.Complex.near(s, Jmat.Complex(1, 18.125), 0.1)) return Jmat.Complex.zetaint_(s); // The code below has a problem for s of the form 1 + n*9.0625i, with n any negative or positive integer
    var result = Jmat.Complex(0);
    var p = Jmat.Complex.TWO.pow(Jmat.Complex.ONE.sub(s));
    var a = Jmat.Complex.ONE.div(Jmat.Complex.ONE.sub(p));
    var n = 1;
    for(var i = 1; i < 128; i++) {
      var r = Jmat.Complex(n).div(Jmat.Complex(i).pow(s));
      result = result.add(r);
      n = -n;
    }
    return result.mul(a);
  } else {
    // Reflection: The functional equation: zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
    var s2 = Jmat.Complex.ONE.sub(s);
    var a = Jmat.Complex.TWO.pow(s);
    var b = Jmat.Complex.PI.pow(s.subr(1));
    var c = Jmat.Complex.sin(s.mulr(Math.PI / 2));
    var d = Jmat.Complex.gamma(s2);
    var e = Jmat.Complex.zeta(s2);
    return a.mul(b).mul(c).mul(d).mul(e);
  }
};

// Dirichlet eta function
Jmat.Complex.eta = function(s) {
  //The calculation only works for s.re > 0.5, so use reflection formula if needed
  if(s.re < 0.5) {
    s = s.neg();

    var a = (Jmat.Complex.ONE.sub(Jmat.Complex.TWO.pow(s.neg().subr(1))));
    var b = (Jmat.Complex.ONE.sub(Jmat.Complex.TWO.pow(s.neg())));
    var c = a.div(b).mulr(2).mul(Jmat.Complex.PI.pow(s.neg().subr(1)));
    var d = Jmat.Complex.sin(Jmat.Complex.PI.mul(s).divr(2));
    var e = Jmat.Complex.gamma(s);
    var f = Jmat.Complex.eta(s.inc());
    return c.mul(s).mul(d).mul(e).mul(f);
  }

  // Borwein's method
  var n = 50;

  // function works on real values
  var d = function(k, n) {
    var result = 0;
    for(var i = 0; i <= k; i++) {
      var a = Jmat.Real.factorial(n + i - 1) / Jmat.Real.factorial(n - i);
      var b = Math.pow(4, i) / Jmat.Real.factorial(2 * i);
      result += a * b;
    }
    return result * n;
  };

  var dn = d(n, n);

  var result = Jmat.Complex(0);
  var sign = Jmat.Complex(1); //alternates
  for(var k = 0; k < n; k++) {
    var a = sign.mulr(d(k, n) - dn);
    var b = Jmat.Complex(k + 1).pow(s);
    result = result.add(a.div(b));
    sign = sign.neg();
  }

  return result.mulr(-1 / dn);

  /*
  NOTE:
  Riemann zeta, Direchlet lambda and Direchlet eta are related as follows:
  zeta(s)/(2^s) = lambda(s)/(2^s - 1) = eta(s)/(2^s - 2)
  also:
  zeta(s) + eta(s) = 2*lambda(s)
  So one of these could be expressed as the other, so that a more precise one can be used.
  E.g. the Borwein's method here looks more precise than the current formulas used in Jmat.Complex.zeta...
  */
};

// Dirichlet lambda function
// TODO: use numerically more precise formula
Jmat.Complex.lambda = function(s) {
  // definition: lambda(s) = (1 - s^(-s))*zeta(s)
  return Jmat.Complex.ONE.sub(Jmat.Complex.TWO.pow(s.neg())).mul(Jmat.Complex.zeta(s));
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
print out with JSON.stringify(Jmat.Complex.hurwitzzeta_generate_hasse_table_(30)), and use in Jmat.Complex.hurwitzzeta_hasse_series_.
*/
Jmat.Complex.hurwitzzeta_generate_hasse_table_ = function(steps) {
  var result = [];
  for(var n = 0; n < steps; n++) {
    result[n] = 0;
    var sign = 1;
    for(var k = 0; k <= n; k++) {
      result[k] += sign * Jmat.Real.pascal_triangle(n, k) / (n + 1);
      sign = -sign;
    }
  }
  return result;
};

Jmat.Complex.hurwitzzeta_hasse_tables_ = [];

// Series by Helmut Hasse: zeta(1 / (s - 1)) * SUM_n=0..oo 1/(n+1) SUM_k=0..n ((-1)^k binomial(n, k) (q+k)^(1-s))
// Requires: in theory: complex s not equal to 1, real q > 0 (or > -1 according to other sources)
// In practice: any complex s with s.im > 0, complex q with q.re > 2 (for smaller q.re, it becomes very imprecise)
Jmat.Complex.hurwitzzeta_hasse_series_ = function(s, q) {
  // Precomputed table implementation, exact same result as the slower version but faster.
  // Generated with JSON.stringify(Jmat.Complex.hurwitzzeta_generate_binomial_table_(30))
  // Table size must exactly match the N in "k < N" in the loop below.
  // TODO: there is a weird noisy area in the center of the following plot. fix it? Jmat.plotComplex(function(z){return Jmat.hurwitzzeta(-5, z); })
  var N = 30;
  if(!Jmat.Complex.hurwitzzeta_hasse_tables_[N]) {
    Jmat.Complex.hurwitzzeta_hasse_tables_[N] = Jmat.Complex.hurwitzzeta_generate_hasse_table_(N);
  }
  var table = Jmat.Complex.hurwitzzeta_hasse_tables_[N];
  var result = Jmat.Complex.ZERO;
  for(var k = 0; k < N; k++) {
    result = result.add(q.addr(k).pow(Jmat.Complex.ONE.sub(s)).mulr(table[k]));
  }
  result = result.mul(s.dec().inv());
  return result;
};

// A series that works for negative s, but q only in range 0-1
// requirements: s.re < 0 && q.re > 0 && q.re <= 1 && q.im == 0
// TODO: this one actually has more potential than it seemed, it works best of all formulas for Jmat.Complex.hurwitzzeta_cos_series_(Jmat.Complex(-10,-10), Jmat.Complex(-2.1)). Investigate area where it works and make more use of it.
Jmat.Complex.hurwitzzeta_cos_series_ = function(s, q) {
  // Series representation 25.11.9 from NIST handbook
  var s1 = s.rsub(1);
  var g = Jmat.Complex.gamma(s1).mulr(2).div(Jmat.Complex(2 * Math.PI).pow(s1));
  var sum = Jmat.Complex.ZERO;
  for(var n = 1; n < 30; n++) {
    var c = Jmat.Complex.cos(s1.mulr(Math.PI / 2).sub(q.mulr(2 * Math.PI * n)));
    sum = sum.add(c.div(Jmat.Complex(n).pow(s1)));
  }
  return sum.mul(g);
};

// The series of the definition of hurwitz zeta.
// This algorithm works reliable for, according to my practical testing: s with s.re > 2 and complex q with q.re > 0
// Theoretically it works for s.re > 1 but it converges too slow there, and would work for any complex s with q.re > 0, but it seems to work for some (but not reliably) negative s as well, and work NOT well for |q| too big.
Jmat.Complex.hurwitzzeta_simple_series_ = function(s, q) {
  if(Jmat.Complex.isNegativeIntOrZero(q)) return Jmat.Complex(Infinity);
  // would work for s.re > 1, but too slow convergence

  var result = Jmat.Complex.ZERO;
  for(var n = 0; n < 30; n++) {
    var r = q.addr(n).pow(s).inv();
    result = result.add(r);
  }

  return result;
};

// Euler-Maclaurin algorithm for hurwitz zeta, basically a sped up version of the simple series
// See (7.1) in paper "efficient algorithm for accelerating the convergence of ..." by Linas Vepstas.
// This algorithm works for: complex s with s.re > 0, and complex q with q.re > -20 (or even > -1). For q with very large negative real part, it goes numerically very wrong.
Jmat.Complex.hurwitzzeta_euler_ = function(s, q) {
  var result = Jmat.Complex.ZERO;
  var N = 25; // D/2 + 10 where D is desired decimal digits of precision
  var p = 15; // TODO: find better value for it? It should not go above 15 though, the bernoulli number table ends at 30 (and k*2 is needed)
  var fn = q.addr(N).pow(s).inv(); // f(N)
  var ig = s.dec().inv().mul(q.addr(N).pow(s.dec()).inv());
  var sum1 = Jmat.Complex.ZERO;
  for(var k = 0; k < N; k++) {
    var r = q.addr(k).pow(s).inv(); // f(k)
    sum1 = sum1.add(r);
  }
  var ss = s; //for the derivative
  var kk = 1; //for the (2k)!
  var sum2 = Jmat.Complex.ZERO;
  for(var k = 1; k <= p; k++) {
    kk *= (2 * k - 1) * (2 * k);
    ss = ss.mul(s.addr(k * 2 - 1)).mul(s.addr(k * 2));
    var d = ss.div(q.addr(N).pow(s.addr(2 * k + 1))).neg();
    sum2 = sum2.add(d.mulr(Jmat.Complex.bernoulli[2 * k] / kk));
  }
  return sum1.sub(sum2).add(fn.divr(2)).add(ig);
};

// Hurwitz zeta function
// Currently works ok for complex s and q with real parts > 0. For negative real parts in s or q, there are some zones of very high imprecision, but works quite well for reasonable values.
// TODO: fix several problems, e.g. see the plot Jmat.plotComplex(function(z){return Jmat.hurwitzzeta(-10, z); })
// TODO: support the cases outside of the currently supported domain, e.g. large negative values, ...
Jmat.Complex.hurwitzzeta = function(s, q) {
  // Bernoulli polynomials for a few small negative integer s values
  // It's not that necessary to have this and the riemann zeta conversion below, but it helps as a check that the algorithms below are correct, if not these values would show up as different lines in 2D plots
  if(s.eqr(0)) return Jmat.Complex(0.5).sub(q);
  if(s.eqr(-1)) return Jmat.Complex(-1/12.0).add(q.mulr(0.5)).add(q.mul(q).mulr(-0.5));

  if(q.eqr(1)) return Jmat.Complex.zeta(s); // riemann zeta
  
  if(s.re > 2 && q.re > 0) {
    return Jmat.Complex.hurwitzzeta_simple_series_(s, q);
  }

  if(s.re >= 0 && q.re > -10) {
    return Jmat.Complex.hurwitzzeta_euler_(s, q);
  }

  // Only a very limited zone in which this works.....
  if(s.re < 0 && q.im == 0 && q.re > 0 && q.re <= 1) {
    // TODO: try the formula with both cos and sin (formula 8 at http://mathworld.wolfram.com/HurwitzZetaFunction.html)
    //       does support complex q. Then use the formulas for Distant neighbors at http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/Zeta2/17/01/ShowAll.html
    //       to bring real part of q in range 0-1, and see if that would make it support all complex q for all complex negative s...
    //       or just make q.re > 2 to make hurwitzzeta_hasse_series_ work ...
    return Jmat.Complex.hurwitzzeta_cos_series_(s, q);
  }

  if(q.re > 2) {
    return Jmat.Complex.hurwitzzeta_hasse_series_(s, q);
  }

  // The smaller, the more iterations needed, hence limited to such number...
  if(q.re > -40 && !q.eqr(0) /* 0 causes infinite recursion*/) {
    //zeta(s, q) = zeta(s, q + m) + SUM_n=0..(m-1) (1/(n+q)^s)
    //if m < 0:
    //zeta(s, a) = SUM_n=0..(m-1) (1/(n+a-m)^s) - zeta(s, a - m)
    // For s.re < 0 and real q, hurwitzzeta_cos_series_ works best, so bring q in range 0-1. Else, make q.re > 2 for hurwitzzeta_hasse_series_.
    var m = (s.re < 0 && q.im == 0) ? Math.ceil(-q.re) : Math.ceil(2 - q.re + 0.5);
    var sum = Jmat.Complex.ZERO;
    if(m > 0) {
      for(var n = 0; n < m; n++) {
        var r = q.addr(n).pow(s).inv();
        if(!Jmat.Complex.isNaN(r)) sum = sum.add(r);
      }
    } else {
      for(var n = 0; n < -m; n++) {
        var r = q.addr(n).addr(m).pow(s).inv();
        if(!Jmat.Complex.isNaN(r)) sum = sum.add(r);
      }
    }
    var h = Jmat.Complex.hurwitzzeta(s, q.addr(m));
    return m > 0 ? h.add(sum) : h.sub(sum);
  }

  // (TODO: use the functional equation, that will allow to support negative s for integer and rational q (with small denominator, otherwise it gets too many calculations))

  // TODO: support s.re < 0 for all q
  return Jmat.Complex(NaN); // not supported in this region :( e.g.
};

Jmat.Complex.beta = function(x, y) {
  // definition: beta(x, y) = gamma(x)*gamma(y) / gamma(x+y)
  //return Jmat.Complex.gamma(x).mul(Jmat.Complex.gamma(y)).div(Jmat.Complex.gamma(x.add(y)));
  //return Jmat.Complex.exp(Jmat.Complex.loggamma(x).add(Jmat.Complex.loggamma(y)).sub(Jmat.Complex.loggamma(x.add(y))));

  // For the negative integers (which cause the gamma function to return indeterminate results) for which beta is still defined:
  // beta(x, y), with x < 0 and y > 0, is, by the formula: gamma(y)*gamma(x) / gamma(x + y), rewritten to: gamma(y) / (gamma(x + y) / gamma(x)), for example:
  // beta(x, 2) = 1 / x(x+1),  beta(x, 3) = 2 / x(x+1)(x+2),  beta(x, 4) = 6 / x(x+1)(x+2)(x+3),  beta(x, 5) = 24 / x(x+1)(x+2)(x+3)(x+4), etc...
  // When x is negative integer, x(x+1)(x+2)...(x+y-1) can be rewritten as (-1^y) * |x|(|x|-1)(|x|-2)...(|x|-y-1), which is (-1^y) * gamma(|x|+1) / gamma(|x|-y+1)
  // And those last gamma functions get positive argument, so it works without NaN again, so we get the right answers
  // The Jmat.Complex.gammaDiv_ function is used to get similar effect. Use a negative input value as its numerator.

  if(x.re < 50 && x.re > -50 && y.re < 50 && y.re > -50) {
    // gamma rather than loggamma more precise here
    return Jmat.Complex.gammaDiv21_(x, y, x.add(y));
  } else {
    return Jmat.Complex.exp(Jmat.Complex.loggammaDiv21_(x, y, x.add(y)));
  }
};

// Incomplete beta function
// B_x(a, b)
Jmat.Complex.incbeta = function(x, a, b) {
  // for x = 1, incbeta(1, a, b) = beta(a, b)
  // otherwise: integrate(0..x, t^(a-1) * (1-t)^(b-1) * dt)

  // TODO: gives incorrect result for incbeta(5, 100, 100). Should be -1.4e127. Gives +1.1e141. Such values are needed by fisher F distribution though.

  if(x.eqr(1)) return Jmat.Complex.beta(a, b);

  // METHOD A: series representation (simpler case of hypergeometric series). Probably precise, but does not work for neg a/b, x out of range 0-1 (0.75 for better convergence), ...
  // The summation is z^a * SUM_(n=0..inf) ((1-b)_n / (a + n) * x^n / n!
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  if(Jmat.Complex.isPositive(x) && x.re < 0.75 && Jmat.Complex.isPositive(a) && Jmat.Complex.isPositive(b)) {
    var b1 = Jmat.Complex.ONE.sub(b);
    var r;
    var result;
    for(var n = 0; n < 30; n++) {
      if(n == 0) {
        r = Jmat.Complex(1);
        result = Jmat.Complex(0);
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
  return x.pow(a).div(a).mul(Jmat.Complex.hypergeometric(a, Jmat.Complex.ONE.sub(b), a.inc(), x));

  // METHOD C: with the integral definition. However this one is numerically very imprecise
  // The integral is numerically too imprecise, due to the degenerate value at t=0
  // the start should be 0, but it's degenerate there
  /*if(Jmat.Complex.isPositive(x) && x.re < 1 && Jmat.Complex.isPositive(a) && Jmat.Complex.isPositive(b)) {
    return Jmat.Complex.integrate(x.divr(3000), x, function(t) {
      var r = t.pow(a.subr(1)).mul(Jmat.Complex.ONE.sub(t).pow(b.subr(1)));
      return r; 
    }, 30);
  }*/
};

// Regularized incomplete beta: (incomplete beta / beta)
// I_x(a, b)
Jmat.Complex.beta_i = function(x, a, b) {
  if(Jmat.Complex.isNegativeIntOrZero(a) && !Jmat.Complex.isNegativeIntOrZero(b)) return Jmat.Complex(1);
  if(Jmat.Complex.isNegativeIntOrZero(b) && !Jmat.Complex.isNegativeIntOrZero(a)) return Jmat.Complex(0);
  if(Jmat.Complex.isNegativeIntOrZero(b.add(a)) && x.eqr(1)) return Jmat.Complex(1);
  return Jmat.Complex.incbeta(x, a, b).div(Jmat.Complex.beta(a, b));
};

// Inverse of regularized incomplete beta: I_x^-1(a, b)
Jmat.Complex.beta_i_inv = function(x, a, b) {
  if(!(Jmat.Complex.isPositiveOrZero(x) && Jmat.Complex.isPositiveOrZero(a) && Jmat.Complex.isPositiveOrZero(b))) return Jmat.Complex(NaN);

  var bab = Jmat.Complex.beta(a, b);
  var x2 = Jmat.Complex.ZERO;
  var a2 = Jmat.Complex.ZERO;
  var b2 = Jmat.Complex.ONE;

  while(!Jmat.Complex.near(a2, b2, 1e-6)) {
    x2 = a2.add(b2).divr(2);
    if(Jmat.Complex.incbeta(x2, a, b).div(bab).re > x.re) b2 = x2;
    else a2 = x2;
  }
  return x2;
};

// Confluent hypergeometric limit function 0F1(; a; z) --> yep, that's the notation of it
Jmat.Complex.hypergeometric0F1 = function(a, z) {
  // The summation is SUM_(n=0..inf) z^n / ((a)_n * n!)
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  // This loop converges very nicely and result is quite accurate for all complex input values (unlike for 2F1)
  var r = Jmat.Complex.ONE;
  var result = Jmat.Complex.ZERO;
  for(var n = 0; n < 30; n++) {
    if (n > 0) {
      if(!r.eqr(0)) r = r.div(a.addr(n - 1)); // the Pochhammer. The if avoids 0/0
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
Jmat.Complex.hypergeometric1F1 = function(a, b, z) {
  // The summation is SUM_(n=0..inf) ((a)_n / (b)_n) * z^n / n!
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  // This loop converges very nicely (unlike for 2F1)
  var r = a.div(b).mul(z).divr(1);
  var result = Jmat.Complex(1);
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
Jmat.Complex.hypergeometric = function(a, b, c, z) {
  if(Jmat.Complex.absr(z) > 1.0001 /*not 1 to avoid infinite loop if abs 1/z is also > 1 due to numeric problems, e.g. for 0.9726962457337884+i0.23208191126279865*/) {
    // The series converges only for |z| < 1. But there are some linear transformations
    // that can convert it to a different z. There are conditions though, and some are
    // more complex than others (e.g. requiring gamma functions).
    // TODO: with only those two supported transformations below, there is probably a lot wrong,
    //       for example, the points at exp(pi*i/3) and exp(-pi*i/3) are known to be not covered by this
    //       To fix, do as in the paper "Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Poschl-Teller-Ginocchio potential wave functions"
    
    // Linear transformations to bring |z| in value < 1
    var z2 = z.div(z.dec());
    if(Jmat.Complex.absr(z2) < 0.75) { // the if can in theory do "< 1", but then a value like z=0.3+4i makes z2 very close to 1, and it has bad numeric results for such values
      // z / (z - 1)
      return Jmat.Complex.ONE.sub(z).pow(a.neg()).mul(Jmat.Complex.hypergeometric(a, c.sub(b), c, z2));
    } else {
      // 1 / z
      var quirkyGammaDiv22_ = function(a, b, c, d) {
        // For some cases where a, b, c - a or c - b are negative integers, the formula doesn't work and requires a rather complicated other formula for the solution.
        // For now, I temporarily instead twiddle the parameters a bit. TODO: that is evil, do it properly (but it is surprisingly somewhat accurate though... well, not for everything)
        var result = Jmat.Complex.gammaDiv22_(a, b, c, d);
        if(Jmat.Complex.isNaN(result)) {
          if(Jmat.Complex.isNegativeIntOrZero(a)) a = a.addr(1e-5);
          if(Jmat.Complex.isNegativeIntOrZero(b)) b = b.addr(1e-5);
          if(Jmat.Complex.isNegativeIntOrZero(c)) c = c.addr(1e-5);
          if(Jmat.Complex.isNegativeIntOrZero(d)) d = d.addr(1e-5);
          result = Jmat.Complex.gammaDiv22_(a, b, c, d);
        }
        return result;
      };

      var zi = z.inv();
      var za = z.neg().pow(a.neg());
      var zb = z.neg().pow(b.neg());
      var ga = quirkyGammaDiv22_(c, b.sub(a), b, c.sub(a));
      var gb = quirkyGammaDiv22_(c, a.sub(b), a, c.sub(b));
      var fa = Jmat.Complex.hypergeometric(a, Jmat.Complex.ONE.sub(c).add(a), Jmat.Complex.ONE.sub(b).add(a), zi);
      var fb = Jmat.Complex.hypergeometric(b, Jmat.Complex.ONE.sub(c).add(b), Jmat.Complex.ONE.sub(a).add(b), zi);
      var va = ga.mul(za).mul(fa);
      var vb = gb.mul(zb).mul(fb);
      return va.add(vb);
    }
  }

  var z2 = z.div(z.dec());
  if(Jmat.Complex.absr(z2) < Jmat.Complex.absr(z)) {
    // Same z / (z - 1) transform as above. Reason for doing this: the summation below converges faster for smaller absolute values of z. Without this, for e.g. z = -0.75 and c < -3, it converges only after hundreds of steps.
    // TODO: in fact it's almost always possible to make |z| < 0.5 with the linear transformations. Use them better.
    return Jmat.Complex.ONE.sub(z).pow(a.neg()).mul(Jmat.Complex.hypergeometric(a, c.sub(b), c, z2));
  }

  // TODO: avoid needing more than 30 iterations by always getting |z| < 0.5 with more clever use of transformations
  var num_iterations = 30;
  if(Jmat.Complex.absr(z) > 0.5) num_iterations = 60;
  if(Jmat.Complex.absr(z) > 0.75) num_iterations = 100;

  // The summation definition of the series, converges because |z| < 1

  // The summation is SUM_(n=0..inf) ((a)_n * (b)_n / (c)_n) * z^n / n!
  // Where (q)_n is the rising Pochhammer symbol x(x+1)*...*(x+n-1)
  // The variable r below gets updated every step with next factor of Pochhammer symbol, power and factorial values as the loop goes.
  var r = a.mul(b).div(c).mul(z).divr(1);
  var result = Jmat.Complex(1);
  for(var n = 1; n < num_iterations; n++) {
    if (n > 1) {
      r = r.mul(a.addr(n - 1));
      r = r.mul(b.addr(n - 1));
      if(!r.eqr(0)) r = r.div(c.addr(n - 1)); // The if is there so that if a==c or b==c and c is neg integer, we don't get 0/0 here but just 0
      r = r.mul(z);
      r = r.divr(n);
    }
    if(r.eqr(0)) break; // If a or b are negative integer, the loop will terminate early
    if(Jmat.Complex.near(r, Jmat.Complex.ZERO, 1e-15)) break; // In fact why not even break out if it's already converged very near zero
    result = result.add(r);
  }
  return result;

  // Alternative: integration formula (doesn't work well)
  /*var g = Jmat.Complex.gammaDiv12_(c, b, c.sub(b));
  var i = Jmat.Complex.integrate(Jmat.Complex(0), Jmat.Complex(1), function(t) {
    var n = t.pow(b.dec()).mul(Jmat.Complex.ONE.sub(t).pow(c.sub(b).dec()));
    var d = Jmat.Complex.ONE.sub(z.mul(t)).pow(a);
    return n.mul(d);
  }, 30);
  return g.mul(i);*/
};

// TODO: generalized hypergeometric

Jmat.Complex.PIPI6_ = Jmat.Complex(Math.PI * Math.PI / 6);

//Dilogarithm: Li_2(z)
Jmat.Complex.dilog = function(z) {
  if(z.eqr(0)) return Jmat.Complex.ZERO;
  if(z.eqr(1)) return Jmat.Complex.PIPI6_;
  if(z.eqr(+Infinity)) return Jmat.Complex(-Infinity);
  if(z.eqr(-Infinity)) return Jmat.Complex(-Infinity);

  // Do the series expansion for |z| < 1
  var summation = function(z) {
    var result = Jmat.Complex.ZERO;
    var zz = z;
    var N = a <= 0.5 ? 30 : 50;
    for(var i = 1; i < N; i++) {
      var r = zz.divr(i * i);
      result = result.add(r);
      if(Jmat.Complex.near(r, Jmat.Complex.ZERO, 1e-15)) break;
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
    var l = Jmat.Complex.log(z.neg());
    return d.sub(l.mul(l).divr(2)).sub(Jmat.Complex.PIPI6_);
  }

  // The above was done to avoid all the stuff below in the obvious good above cases. Now consider each option, including the above two and other identities, again for worse |z|

  // Three more identities
  var z1 = Jmat.Complex.ONE.sub(z);
  var a1 = z1.absr();
  var z2 = Jmat.Complex.ONE.sub(z).inv();
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
      var l = Jmat.Complex.log(z.neg());
      return d.sub(l.mul(l).divr(2)).sub(Jmat.Complex.PIPI6_);
    }
    if(bestindex == 2) {
      var d = summation(z1).neg();
      var l1 = Jmat.Complex.log(z1);
      var l2 = Jmat.Complex.log(z);
      return d.sub(l1.mul(l2)).add(Jmat.Complex.PIPI6_);
    }
    if(bestindex == 3) {
      var d = summation(z2);
      var l1 = Jmat.Complex.log(Jmat.Complex.ONE.sub(z));
      var l2 = Jmat.Complex.log(z.neg());
      return d.add(l1.mul(l1).divr(2)).sub(l2.mul(l1)).sub(Jmat.Complex.PIPI6_);
    }
    if(bestindex == 4) {
      var d = summation(z3).neg();
      var l = Jmat.Complex.log(Jmat.Complex.ONE.sub(z));
      return d.sub(l.mul(l).divr(2));
    }
  }
  
  // Near values 0.5+i and 0.5-i the above is not working.
  // The following identity is: Li_2(z) = Li_2(z^2)/2 - Li_2(-z)
  // That is, the square relationship
  // This gives z's which are eventually not in the problem region, and thus we will not end up with infinite recursive call of this function.

  return Jmat.Complex.dilog(z.mul(z)).divr(2).sub(Jmat.Complex.dilog(z.neg()));
};


//Trilogarithm: Li_3(z)
Jmat.Complex.trilog = function(z) {
  if(z.eqr(0)) return Jmat.Complex.ZERO;
  if(z.eqr(1)) return Jmat.Complex.APERY;
  if(z.eqr(+Infinity)) return Jmat.Complex(-Infinity);
  if(z.eqr(-Infinity)) return Jmat.Complex(-Infinity);

  // Do the series expansion for |z| < 1
  var summation = function(z) {
    var result = Jmat.Complex.ZERO;
    var zz = z;
    var N = a < 0.5 ? 30 : 50;
    for(var i = 1; i < N; i++) {
      var r = zz.divr(i*i*i);
      result = result.add(r);
      if(Jmat.Complex.near(r, Jmat.Complex.ZERO, 1e-15)) break;
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
    var l = Jmat.Complex.log(z.neg());
    return d.sub(l.mul(l).mul(l).divr(6)).sub(Jmat.Complex.PIPI6_.mul(l));
  }

  // TODO: better implementation for this case
  return Jmat.Complex.polylog_integral_(Jmat.Complex(3), z);
};

// Bernoulli numbers B_m. Only supports a hardcoded first few dozens.
Jmat.Complex.bernoulli = [
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
Jmat.Complex.stieltjes = [
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
Jmat.Complex.stieltjes_zeta = [
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
Jmat.Complex.polylog_integral_ = function(s, z) {
  // To test these, try e.g.:
  // complexDomainPlot(function(z){return Jmat.Complex.polylog(Jmat.Complex(15, 0.5), z);}, 2, 1);
  // complexDomainPlot(function(z){return Jmat.Complex.polylog(Jmat.Complex(0.5, 15), z);}, 2, 1);

  if(s.re > 1 && Math.abs(s.im) < s.re && Math.abs(z.argr()) > 0.1) {
    // Approximate with an integral representation (only works for real s.re, and has some serious problems with periodic things appearing near the positive real axis of z)
    // In practice, only works for s.re > 1 (!!) and s.im = 0, or very small |s.im| and s.re > 0
    var g = Jmat.Complex.gamma(s);
    var r = Jmat.Complex.integrate(Jmat.Complex(0), Jmat.Complex(20), function(t) {
      var result = t.pow(s.dec()).div(Jmat.Complex.exp(t).sub(z));
      if(Jmat.Complex.isNaN(result)) result = Jmat.Complex.ZERO;
      return result;
    }, 100);
    return z.div(g).mul(r);
  } else if(Jmat.Complex.isNegative(s) && Math.abs(z.argr()) > 0.1) {
    // Approximate with an integral representation (only works for s.re < 0
    var lzm = Jmat.Complex.log(z.neg());
    var r = Jmat.Complex.integrate(Jmat.Complex(0), Jmat.Complex(20), function(t) {
      var ta = t.pow(s.neg());
      var tb = Jmat.Complex.sin(s.mulr(Math.PI/2).sub(t.mul(lzm)));
      var na = Jmat.Complex.sinh(t.mulr(Math.PI));
      var result = ta.mul(tb).div(na);
      if(Jmat.Complex.isNaN(result)) result = Jmat.Complex.ZERO;
      return result;
    }, 100);
    return r;
  } else if(z.im <= 0 || (s.re > 0 && Math.abs(s.im) < s.re)) {
    // Integral formula for polylog that works for all complex s and z, except for positive real z if s.re < 0
    // Is from 0 to infinity, but seems to work well from 0-10 with only 100 steps (at least in the areas where this is used)
    // While in theory it works for all z and all s (except z near positive axis if s.re < 0), in practice it seems to work only for z.im <= 0 and s.re > 0 and |s.im| << s.re
    // --> I tried with z-plot for s=-15+0.5i, s=0.5+15i, and other variations like that
    var lzm = Jmat.Complex.log(z.neg());
    var f = function(t) {
      var ta = Jmat.Complex.sin(s.mul(Jmat.Complex.atan(t)).sub(t.mul(lzm)));
      var na = (Jmat.Complex.ONE.add(t.mul(t))).pow(s.divr(2));
      var nb = Jmat.Complex.sinh(t.mulr(Math.PI));
      var result = ta.div(na).div(nb);
      if(Jmat.Complex.isNaN(result)) result = Jmat.Complex.ZERO;
      return result;
    };
    var r = Jmat.Complex.ZERO;
    r = r.add(Jmat.Complex.integrate(Jmat.Complex(0), Jmat.Complex(5), f, 50));
    r = r.add(Jmat.Complex.integrate(Jmat.Complex(5), Jmat.Complex(20), f, 20));
    r = r.add(Jmat.Complex.integrate(Jmat.Complex(20), Jmat.Complex(100), f, 10));
    return z.mulr(0.5).add(z.mul(r));
  } else if(z.im > 0) { // Because something is broken for z.im <= 0. TODO: find out what
    // Similar integral, but with upper incomplete gamma.
    // It is theoretically better than the above because it supports even its last edge case.
    // But it seems broken. It does not work for z.im <= 0... At least it gets 
    var lz = Jmat.Complex.log(z);
    var f = function(t) {
      var ta = Jmat.Complex.sin(s.mul(Jmat.Complex.atan(t)).sub(t.mul(lz)));
      var na = (Jmat.Complex.ONE.add(t.mul(t))).pow(s.divr(2));
      var nb = Jmat.Complex.exp(t.mulr(2*Math.PI)).sub(1);
      var result = ta.div(na).div(nb);
      if(Jmat.Complex.isNaN(result)) result = Jmat.Complex.ZERO;
      return result;
    };
    var r = Jmat.Complex.ZERO;
    r = r.add(Jmat.Complex.integrate(Jmat.Complex(0), Jmat.Complex(5), f, 50));
    r = r.add(Jmat.Complex.integrate(Jmat.Complex(5), Jmat.Complex(20), f, 20));
    r = r.add(Jmat.Complex.integrate(Jmat.Complex(20), Jmat.Complex(100), f, 10));
    var g = Jmat.Complex.incgamma_upper(Jmat.Complex.ONE.sub(s), lz.neg());
    var l = lz.neg().pow(Jmat.Complex.ONE.sub(s));
    return z.mulr(0.5).add(g.div(l)).add(z.mulr(2).mul(r));
  } else {
    return Jmat.Complex(NaN); //there is nothing we can do... :(
  }
};

// Borwein algorithm. Converges theoretically only for |z^2/(z-1)| < 3.7.
// In practice this zone of convergence is much smaller due to numerical problems. For extreme values of s (e.g. -20), it is as small as 0.5
// This algorithm is actually intended for arbitrary precision libraries, not for floating point.
// I'm keeping it here, because for s.re > -3, and with the adaptive n, it works better than most other polylog code here and is quite accurate... The function "polylog_borwein_ok_" returns whether it's usable for the given values.
// Uses the following formula:
// Li_s(z) = SUM_k=1..n(z^k / k^s) + 1/(1-z)^n * SUM_k=n+1..2n(z^k / k^s * SUM_j=0..2n-k((-z)^j * binomial(n, j))) + error(s, z)
// where the error term is negligible in the zone of convergence
// NOTE: changing z to -1 makes this the Borwein algorithm for riemann zeta.
Jmat.Complex.polylog_borwein_ = function(s, z) {
  var kidney_radius = z.mul(z).div(z.dec()).absr(); //radius of the kidney shaped region of convergence. Theoretically works for < 4, in practice with double precision only for < 3 and then only if not too much n!
  if(kidney_radius >= 3.7) return Jmat.Complex(NaN); //yeah right... it sucks way before this

  // number of loops. NOTE: higher is better for arbitrary precision library (with 31 being perfect), but results in random garbage with floating point. So it is limited here instead.
  var n = Math.floor(Math.min(31, 16 / kidney_radius));


  var binomial_cache = [];
  binomial_cache[0] = Jmat.Complex.ONE;
  var bin_sum_term = Jmat.Complex.ONE; // Binomial summation term

  var zz = Jmat.Complex.ONE;
  var result0 = Jmat.Complex.ZERO;

  var barray = [];

  // First sum, and preparation for binomial sum
  for(var k = 1; k <= n; k++) {
    zz = zz.mul(z); // numerically critical...
    var r = zz.div(Jmat.Complex(k).pow(s));
    result0 = result0.add(r);

    // Binomial sum
    var bin = Jmat.Complex.binomial(Jmat.Complex(n), Jmat.Complex(k));
    var b = bin.mul(zz); // numerically critical...
    if(k % 2) b = b.neg(); // numerically critical...

    bin_sum_term = bin_sum_term.add(b); // numerically critical...

    // Store binomial sums in array for later reference
    binomial_cache[k] = bin_sum_term;
  }

  // Second sum, using the binomial sum preparation from above
  var result1 = Jmat.Complex.ZERO;
  for(k = n + 1; k <= 2 * n; k++) {
    zz = zz.mul(z);
    var r = zz.div(Jmat.Complex(k).pow(s));

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

// Is it ok to use Borwein method?
Jmat.Complex.polylog_borwein_ok_ = function(s, z) {
  var kidney_radius = z.mul(z).div(z.dec()).absr();
  return (s.re >= -5 && z.im == 0 && kidney_radius <= 1.5) || (s.re >= 0 && Math.abs(s.im) < s.re && kidney_radius <=2) || z.absr() <= 0.5;
};

Jmat.Complex.polylog_residue_ = function(s, z) {
  // Sum of residues: LI_s(e^mu) = gamma(1-s) * SUM_k=-oo..oo (2*k*pi*i - mu)^(s-1)
  // Theoretically holds for Re(s) < 0 and e^mu != 1
  // This seems to work pretty well in practice! Even for large z. It works for all z except 1, for all s with negative real part (even things like Jmat.Complex.polylog(Jmat.Complex(-15,0.5), Jmat.Complex(-100,1)) match wolfram alpha polylog(-15+0.5i, -100+i)!!)

  //if(z.re < 0) return square(s, z); //numerically not stable

  var result_is_real = false;
  var origz = z;

  // The only thing this doesn't work well for is for real z with |z| > 2. Twiddle with slightly negative imaginary part to compensate. This because wolfram alpha matches the values for negative imaginary part there as well. NOTE: this is a branch cut line so of course it's tricky here
  if(s.im == 0 && z.im == 0) {
    z = z.add(Jmat.Complex.newi(-0.0001));
    if(Jmat.Complex.isNegativeInt(s) || z.re < 0) result_is_real = true; //if s is not a negative integer and z > 0, then despite s and z being real, the result can be complex. However in other cases, it is real, and the twiddling ruins that, so compensate it at the end.
  }

  var mu = Jmat.Complex.log(z);
  var g = Jmat.Complex.gamma(Jmat.Complex.ONE.sub(s));

  var r = Jmat.Complex.ZERO;
  for(var k = 0; k < 30; k++) {
    // The Wikipedia formula has sum from -Infinity to +Infinity. So there are two terms, a positive and a negative.
    // However, when I try it, for z.re > 0, that makes if off by a factor of two, and for z.re < 0 wrong in some different way. The first term works for z.im > 0, the second for z.im <= 0, I think.
    // So: negative-k term DISABLED!!!
    // TODO: figure out why that is. Because some things do still go wrong with this...
    if(origz.im > 0) r = r.add((Jmat.Complex.newi(2 * k * Math.PI).sub(mu)).pow(s.dec()));
    else r = r.add((Jmat.Complex.newi(-2 * k * Math.PI).sub(mu)).pow(s.dec()));
  }

  var result = g.mul(r);
  if(result_is_real && result.im != 0) result = Jmat.Complex(result.re);
  return result;
};

//Polylogarithm: Li_s(z)
Jmat.Complex.polylog = function(s, z) {
  if(Jmat.Complex.isInt(s)) {
    if(s.eqr(0)) {
      return z.div(Jmat.Complex.ONE.sub(z));
    } else if(s.eqr(1)) {
      return Jmat.Complex.log(Jmat.Complex.ONE.sub(z)).neg();
    } else if(s.eqr(2)) {
      return Jmat.Complex.dilog(z);
    } else if(s.eqr(3)) {
      return Jmat.Complex.trilog(z);
    } else if(s.eqr(-1)) {
      var z1 = Jmat.Complex.ONE.sub(z);
      return z.div(z1).div(z1);
    } else if(s.eqr(-2)) {
      var z1 = Jmat.Complex.ONE.sub(z);
      return z.mul(z.inc()).div(z1).div(z1).div(z1);
    } else if(s.eqr(-3)) {
      var z1 = Jmat.Complex.ONE.sub(z);
      var zz1 = z1.mul(z1);
      return z.mul(Jmat.Complex.ONE.add(z.mulr(4)).add(z.mul(z))).div(zz1).div(zz1);
    } else if(s.eqr(-4)) {
      var z1 = Jmat.Complex.ONE.sub(z);
      var zz1 = z1.mul(z1);
      return z.mul(z.inc()).mul(Jmat.Complex.ONE.add(z.mulr(10)).add(z.mul(z))).div(zz1).div(zz1).div(z1);
    }
  }

  // TODO: there are still problems with:
  // 1) pure imaginary s
  // 2) z near 1 (no matter what s)
  // 3) despite its goodness, there's something wrong with the Sum of residues formula. Normally both sum terms should be used, but some zones seem to require only one. Make it work both for complexDomainPlot(function(z){return Jmat.Complex.polylog(Jmat.Complex(-5, 15), z);}, 2, 1); (where using the other one, makes it black) and plot2D(Jmat.Complex.polylog, 10, 1) (where the complex values in lower left quadrant should be real), and yet some other zones would make it double its value if both are enabled...
  // 4) for s.re > 0, it is a lot less accurate than for s.re < 0, because for s.re < 0 a good series can be used, while for s.re > 0, the inaccurate integrals are used and no inversion formula for is available here

  if(s.re > 1 && z.eqr(0)) return Jmat.Complex(0);
  if(s.re > 1 && z.eqr(1)) return Jmat.Complex.zeta(s); // It's the riemann zeta function for z=1, but only if s.re > 1
  /*if(z.eqr(-1)) return Jmat.Complex.eta(s).neg();*/

  // The duplication formula (square relationship):
  // Li_s(-z) + Li_s(z) = 2^(1-s)*Li_s(z^2)
  // or:
  // Li_s(z) = 2^(1-s)*Li_s(z^2) - Li_s(-z)
  // Li_s(-z) = 2^(1-s)*Li_s(z^2) - Li_s(z)
  var square = function(s, z) {
    var a = Jmat.Complex.polylog(s, z.mul(z));
    var b = Jmat.Complex.polylog(s, z.neg());
    var c = Jmat.Complex.TWO.pow(Jmat.Complex.ONE.sub(s));
    return c.mul(a).sub(b);
  };

  var a = z.absr();

  if((a <= 0.5 && s.re >= -10) || (a < 0.75 && s.re >= -2) || (a < 0.9 && s.re > -1)) {
    // METHOD A: the series definition SUM_k=1..oo z^k / k^s. But only valid for |z| < 1 and converges slowly even then
    // In practice, only useful if |z| < 0.5 and s.re not too small, or higher |z| for larger s.re. If it works, this is one of the most accurate representations...
    var result = Jmat.Complex.ZERO;
    var zz = z;
    var N = a <= 0.5 ? 30 : 50;
    for(var i = 1; i < N; i++) {
      var r = zz.div(Jmat.Complex(i).pow(s));
      result = result.add(r);
      if(Jmat.Complex.near(r, Jmat.Complex.ZERO, 1e-15)) break;
      zz = zz.mul(z);
    }
    return result;
  } else if(Jmat.Complex.polylog_borwein_ok_(s, z)) {
    return Jmat.Complex.polylog_borwein_(s, z);
  } else if(s.re < 0) {
    return Jmat.Complex.polylog_residue_(s, z);
  } else if(Jmat.Complex.isNegativeInt(s) && a > 1) {
    // Inversion formula. Since the "sum of residue" above works so fantastic and is already for all negative s, this code probably never gets called anymore.
    var sign = Jmat.Complex.isOdd(s) ? 1 : -1; // (-1)^(s-1)
    return Jmat.Complex.polylog(s, z.inv()).mulr(sign);
  } else {
    // Last resort: the not very well working integrals :(
    // Polylog is supposed to be very interesting at s=0.5+t*i. But no code here is precise enough for that.
    return Jmat.Complex.polylog_integral_(s, z);
  }
  return Jmat.Complex(NaN);
};

//Jacobi theta1 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
Jmat.Complex.theta1 = function(z, q) {
  // theta1 is an odd function, and the series below has numerical problems when imaginary part is large, due to values of complex sine becoming huge there --> actually I don't think that helps, the same problem with sine is there for large negative imaginary value...
  // therefore, mirror it
  if(z.im > 1) return Jmat.Complex.theta1(z.neg(), q).neg();

  var result = Jmat.Complex(0);
  var s = Jmat.Complex(1);
  for(var n = 0; n < 20; n++) {
    var a = q.powr((n + 0.5) * (n + 0.5));
    var b = Jmat.Complex.sin(z.mulr(2 * n + 1));
    result = result.add(s.mul(a).mul(b));
    /*var a = Jmat.Complex.log(q.powr((n + 0.5) * (n + 0.5)));
    var b = Jmat.Complex.logsin(z.mulr(2 * n + 1));
    result = result.add(s.mul(Jmat.Complex.exp((a).add(b))));*/
    s = s.neg();
  }
  return result.mulr(2);
};

//Jacobi theta2 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
Jmat.Complex.theta2 = function(z, q) {
  // theta2 is an even function, and the series below has numerical problems when imaginary part is large, so mirror it --> actually I don't think that helps, the same problem with sine is there for large negative imaginary value...
  if(z.im > 1) return Jmat.Complex.theta2(z.neg(), q);

  var result = Jmat.Complex(0);
  for(var n = 0; n < 20; n++) {
    var a = q.powr((n + 0.5) * (n + 0.5));
    var b = Jmat.Complex.cos(z.mulr(2 * n + 1));
    result = result.add(a.mul(b));
  }
  return result.mulr(2);
};

//Jacobi theta3 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
Jmat.Complex.theta3 = function(z, q) {
  // theta3 is an even function, and the series below has numerical problems when imaginary part is large, so mirror it
  if(z.im > 1) return Jmat.Complex.theta3(z.neg(), q);

  var result = Jmat.Complex(0);
  for(var n = 0; n < 20; n++) {
    var a = q.powr(n * n);
    var b = Jmat.Complex.cos(z.mulr(2 * n));
    result = result.add(a.mul(b));
  }
  return result.mulr(2).addr(1);
};

//Jacobi theta4 function
//q = exp(pi * i * tau) ==> tau = ln(q) / (pi * i)
//|q| must be < 1
// TODO: look broken for imaginary values > 8
Jmat.Complex.theta4 = function(z, q) {
  // theta4 is an even function, and the series below has numerical problems when imaginary part is large, so mirror it
  if(z.im > 1) return Jmat.Complex.theta4(z.neg(), q);

  var result = Jmat.Complex(0);
  var s = Jmat.Complex(1);
  for(var n = 0; n < 20; n++) {
    var a = q.powr(n * n);
    var b = Jmat.Complex.cos(z.mulr(2 * n));
    result = result.add(s.mul(a).mul(b));
    s = s.neg();
  }
  return result.mulr(2).addr(1);
};

// TODO: exponential integral, sine integral, cosine integral, logarithmic integral, elliptic integrals, fresnel integrals

// Returns 0 or 1
// The reason this function is here as well as in Jmat.Real, while the other prime functions are only in Jmat.Real, is that this one needs to look at the imaginary part and say it's not prime if it's not zero
Jmat.Complex.isPrime = function(value) {
  if(!Jmat.Complex.isReal(value)) return 0; //complex numbers are not prime
  return Jmat.Real.isPrime(value.re);
};

// returns numerator and denominator of fraction
// max = max value for denominator (a real JS number)
Jmat.Complex.decompose = function(x, max) {
  if (Math.abs(x.re) >= Math.abs(x.im)) {
    var nd = Jmat.Real.decompose(x.re, max);
    var im = Math.round(x.im * nd[1]);
    return [Jmat.Complex(nd[0], im), Jmat.Complex(nd[1])];
  } else {
    var nd = Jmat.Real.decompose(x.im, max);
    var re = Math.round(x.re * nd[1]);
    return [Jmat.Complex(re, nd[0]), Jmat.Complex(nd[1])];
  }
};

// n! / (n-p)!
Jmat.Complex.permutation = function(n, p) {
  // gammaDiv_ is already optimized for integers near each other etc...
  return Jmat.Complex.gammaDiv_(n.inc(), n.sub(p).inc());
};

//Binomial coefficient, aka combination(s). Number of rows of p elements that can be made out of n elements, where order doesn't matter.
//    ( n )
//    ( p )
// n! / (p! * (n-p)!)
Jmat.Complex.binomial = function(n, p) {
  if(Jmat.Complex.isPositiveIntOrZero(n) && Jmat.Complex.isPositiveIntOrZero(p) && p.re <= n.re && n.re < 30) return Jmat.Complex(Jmat.Real.pascal_triangle(n.re, p.re));
  
  // gammaDiv_ is already optimized for integers near each other etc...
  var result = Jmat.Complex.gammaDiv12_(n.inc(), p.inc(), n.sub(p).inc());
  // Round to integer if large result, it sometimes gets numerically a bit off.
  if(result.re > 100 && Jmat.Complex.isPositiveInt(n) && Jmat.Complex.isPositiveInt(p) && n.re > p.re) result = Jmat.Complex.round(result);
  return result;
};

//Stirling number of the second kind
//    { n }
//    { k }
// 1/k! * SUM_j=0..k((-1)^(k-j) * combination(k, j) * j^n)
Jmat.Complex.stirling2 = function(n, k) {
  if(!Jmat.Complex.isInt(k)) return Jmat.Complex(NaN); // only defined for integer k

  var result = Jmat.Complex.ZERO;
  var sign = Jmat.Real.isOdd(k.re) ? -1 : 1;
  var j = Jmat.Complex(0);
  for(j.re = 0; j.re <= k.re; j.re++) {
    result = result.add(Jmat.Complex.binomial(k, j).mul(j.pow(n)).mulr(sign));
    sign *= -1;
  }
  return result.div(Jmat.Complex.factorial(k));
};

// Manhattan distance of complex numbers, returned as a real number (JS float)
Jmat.Complex.manhattan = function(a, b) {
  return Math.max(Math.abs(a.re - b.re), Math.abs(a.im - b.im));
};

Jmat.Complex.near = function(x, y, precision) {
  //return Jmat.Complex.manhattan(x, y) <= precision;
  // Manhattan NOT used, because then this function returns false for equal infinities
  return x.re - precision <= y.re && x.re + precision >= y.re && x.im - precision <= y.im && x.im + precision >= y.im;
};

// Lambertw for branch (0 = principal branch Wp, -1 is also common (Wm))
// Branch is real integer, z is Jmat.Complex object (complex)
Jmat.Complex.lambertwb = function(branch, z) {
  if(Jmat.Complex.isReal(z) && z.re > -0.36 /*~ -1/e*/ && branch == 0) return Jmat.Complex(Jmat.Real.lambertw(z));

  if(!Jmat.Real.isInt(branch)) return Jmat.Complex(NaN);


  // Known special values
  if(Jmat.Complex.isNaN(z)) return NaN;
  if(Jmat.Complex.isInf(z)) return Jmat.Complex(Infinity); // any complex infinity gives positive infinity
  if(branch == 0 && z.re == 0 && z.im == 0) return Jmat.Complex(0);
  if(branch != 0 && z.re == 0 && z.im == 0) return Jmat.Complex(-Infinity); //at all other branch than the principal one, it's -infinity at 0

  /*
  Choosing a good starting value is important. Jmat.Complex(0) as starting value works
  most of the time, but does not work at some regions in the negative complex domain,
  e.g. around 5.4+0.1i, 5.5+0.1i, ... and that can be seen as mandelbrot-fractal-like
  circles around those regions in the complex domain plot.
  */
  var w = Jmat.Complex.log(z).add(Jmat.Complex(0, branch * Math.PI * 2));
  if(branch == 0 && z.absr() < 1.2 /*supposed to be 1/Math.E, but I still see problems beyond that in the complex domain plot...*/) {
    w = Jmat.Complex.sqrt(z.mulr(5.43656365691809047).addr(2)).add(Jmat.Complex(-1, branch * Math.PI * 2));
  }
  if(branch != 0 && z.im == 0) z.im += 1e-14; // Give it small imaginary part, otherwise it never gets there

  var num = 36;
  for(var i = 0; i < num; i++) {
    var ew = Jmat.Complex.exp(w);
    var wew = w.mul(ew);
    var t = wew.sub(z);
    var a = ew.mul(w.addr(1));
    var b = w.addr(2).mul(t).div(w.mulr(2).addr(2));
    w = w.sub(t.div(a.sub(b)));

    var ltest = Jmat.Complex.log(z.div(w)); //for testing if near (z = w*exp(w) OR ln(z/w) = w)
    if(Jmat.Complex.near(ltest, w, 1e-16) || Jmat.Complex.near(wew, z, 1e-16)) break;
    if(i + 1 == num && !(Jmat.Complex.near(ltest, w, 1) || Jmat.Complex.near(wew, z, 1))) return Jmat.Complex(NaN); // iteration could not finish and too far from result
  }

  // Remove numeric tiny imaginary part that appeared in error
  if(z.im == 0 && z.re >= 0) w.im = 0;

  return w;
};

// Principal branch of Lambert's W function: Wp, inverse (not reciprocal) of exp(x) * x
Jmat.Complex.lambertw = function(z) {
  return Jmat.Complex.lambertwb(0, z);
};

// Negative branch of Lambert's W function: Wm, inverse (not reciprocal) of exp(x) * x
Jmat.Complex.lambertwm = function(z) {
  // TODO: wrong. Look at the real plot. Fix this! Jmat.plotReal(Jmat.Complex.lambertwm)
  return Jmat.Complex.lambertwb(-1, z);
};

// Tetration
// Returns experimental (not mathematically correct) results unless z is a positive integer or Infinity
Jmat.Complex.tetration = function(a, z) {
  if(Jmat.Complex.isPositive(a) && Jmat.Complex.isPositiveInt(z) && z.re != Infinity) {
    return Jmat.Complex(Jmat.Real.tetration(a.re, z.re));
  }

  //if(a.eqr(1)) return Jmat.Complex(1);  // Not true, e.g. 1^Infinity is not 1. Probably things are not like this for tetration either. Indeed it looks like e.g. 1^^0.5 is not 1.
  if(z.eqr(0)) return Jmat.Complex(1); //by definition
  if(z.eqr(1)) return a;
  if(z.eqr(2)) return a.pow(a);
  if(Jmat.Complex.isReal(a) && a.re >= 2 && z > 5) return Jmat.Complex(Infinity); // too big for double
  if(a.eqr(0) && Jmat.Complex.isPositiveInt(z)) {
    // Chosen definition that considers 0^0 = 1, so 0^^n is 1 if n is even, 0 if n is odd.
    return Jmat.Complex.isEven(z) ? Jmat.Complex(1) : Jmat.Complex(0);
  }

  // Power tower (infinitely iterated exponentiation)
  if(z.eqr(Infinity) /*&& Jmat.Complex.isPositive(a)*/) {
    if(a.eqr(1)) return Jmat.Complex(1); //0/0 ==> 1
    // converges if a >= 0.066 && a <= 1.44
    // when using "runloop" with high iterations here, it it indeed only converges there. The lambertw formula below has values everywhere though.
    var l = Jmat.Complex.log(a);
    return Jmat.Complex.lambertw(l.neg()).div(l.neg());
  }

  var runloop = function(a, b, num, l) {
    var result = b;
    var last;
    for(var i = 0; i < num; i++) {
      if(l) result = Jmat.Complex.logy(result, a);
      else result = a.pow(result);
      if(Jmat.Complex.isNaN(result)) return result;
      if(result.eq(last)) return result; // E.g. a=1.01 already converges to 1.0101015237405409 after just 7 iterations. 1.44 converges to 2.393811748202943 after 244 iterations, so numbers near that take a while to loop. 1.45 and higher converge to Infinity.
      last = result;
      if(i > 1000) return Jmat.Complex(NaN); //avoid infinite loop
    }
    return result;
  };

  // This supports power tower as well
  // It supports it better than the lambertw formula above because the domain here is real numbers only, the lambertw formula fails on, for example, -1 and returns nan instead of -1
  if(Jmat.Complex.isPositiveInt(z)) {
    return runloop(a, a, z.re - 1, false);
  }

  // Everything above is true tetration for those cases where possible. What follows below is intermediate tetration research, to return "something", because there is no more than that available.
  if (Jmat.Complex.isReal(z)) {
    // Linear approximation for the extension to real heights
    // a^^z = z+1 for z > -1 && z <= 0
    // a^^z = log_a(a^^(z+1)) for z <= -1  ==> a^^-1.5 = log_a(z+2), a^^-2.5 = log_a(log_a(z+3)), etc...
    // a^^z = a^(a^^(z-1)) for z > 0  ==> a^^0.5 = a^z, a^^1.5 = a^(a^(z-1)), a^^2.5 = a^(a^(a^(z-2))), etc...
    // examples: e^^(0.5*pi) ~= 5.868, 0.5^^(-4.3) ~= 4.033
    if(z.eqr(-1)) {
      return Jmat.Complex(0); //the formulas below would give -Infinity
    }
    if(z.re > -1 && z.re <= 0) {
      return z.inc();
    }
    if(z.re > 0) {
      var b = z.sub(Jmat.Complex.floor(z)); //always in range 0-1
      return runloop(a, b, Math.ceil(z.re), false);
    }
    if(z.re <= -1) {
      var b = z.sub(Jmat.Complex.floor(z)); //always in range 0-1
      return runloop(a, b, -Math.ceil(z.re), true);
    }
  }

  if (Jmat.Complex.near(a, Jmat.Complex.E, 1e-15)) {
    // This implementation, which only works for base e, is based on the paper: Tetration as special function, E. Kouznetsov

    // TODO: more coefficients
    var s = [0.30685281944005, 0.59176735125832, 0.39648321290170, 0.17078658150959,
             0.08516537613999, 0.03804195209047, 0.01734090876306, 0.00755271038865,
             0.00328476064839, 0.00139361740170, 0.00058758348148, 0.00024379186661,
             0.00024379186661, 0.00010043966462, 0.00001654344436, 0.00000663102846,
             0.00000264145664, 0.00000104446533, 0.00000041068839, 0.00000016048059,
             0.00000006239367, 0.00000002412797, 0.00000000928797, 0.00000000355850,
             0.00000000135774, 0.00000000051587];
    // TODO: does JavaScript always re-execute all these constructors? Make object constant outside of this function instead
    var t = [Jmat.Complex(0.37090658903229, 1.33682167078891), Jmat.Complex(0.01830048268799, 0.06961107694975),
             Jmat.Complex(-0.04222107960160, 0.02429633404907), Jmat.Complex(-0.01585164381085, -0.01478953595879),
             Jmat.Complex(0.00264738081895, -0.00657558130520), Jmat.Complex(0.00182759574799, -0.00025319516391),
             Jmat.Complex(0.00036562994770, 0.00028246515810), Jmat.Complex(0.00002689538943, 0.00014180498091),
             Jmat.Complex(-0.00003139436775, 0.00003583704949), Jmat.Complex(-0.00001376358453, -0.00000183512708),
             Jmat.Complex(-0.00000180290980, -0.00000314787679), Jmat.Complex(0.00000026398870, -0.00000092613311),
             Jmat.Complex(0.00000024961828, -0.00000013664223), Jmat.Complex(0.00000000637479, 0.00000002270476),
             Jmat.Complex(-0.00000000341142, 0.00000000512289), Jmat.Complex(-0.00000000162203, 0.00000000031619),
             Jmat.Complex(-0.00000000038743, -0.00000000027282), Jmat.Complex(-0.00000000001201, -0.00000000013440),
             Jmat.Complex(0.00000000002570, -0.00000000002543), Jmat.Complex(0.00000000000935, 0.00000000000045),
             Jmat.Complex(0.00000000000170, 0.00000000000186), Jmat.Complex(-0.00000000000005, 0.00000000000071),
             Jmat.Complex(-0.00000000000016, 0.00000000000012), Jmat.Complex(-0.00000000000005, -0.00000000000001),
             Jmat.Complex(-0.00000000000001, -0.00000000000001)];
    var an = [Jmat.Complex(0.3181315052047641353, 1.3372357014306894089), Jmat.Complex(1), Jmat.Complex(-0.1513148971556517359, -0.2967488367322413067),
              Jmat.Complex(-0.03697630940906762, 0.09873054431149697), Jmat.Complex(0.0258115979731401398, -0.017386962126530755), Jmat.Complex(-0.0079444196, 0.00057925018)];
    var fima = function(z) {
      var r = Jmat.Complex(1.0779614375280, -0.94654096394782);
      var beta = Jmat.Complex(0.12233176, -0.02366108);
      var l = an[0]; //fixed point of the logarithm
      var e = Jmat.Complex.exp(l.mul(z).add(r));
      var c = beta.mul(e).mul(Jmat.Complex.exp(z.mul(Jmat.Complex.newi(Math.PI * 2))));
      return Jmat.Complex.powerSeries(an, an.length, Jmat.Complex.ZERO, e).add(c);
    };

    var b;
    var z2 = Jmat.Complex(Jmat.Real.fracn(z.re), z.im);

    if(z.im < -4.5) b = fima(z2.conj()).conj();
    else if(z.im < -1.5) b = Jmat.Complex.powerSeries(t, t.length, Jmat.Complex.newi(3), z2.conj()).conj();
    else if(z.im < 1.5) b = Jmat.Complex.log(z2.addr(2)).add(Jmat.Complex.powerSeries(s, s.length, Jmat.Complex.ZERO, z2));
    else if(z.im < 4.5) b = Jmat.Complex.powerSeries(t, t.length, Jmat.Complex.newi(3), z2);
    else b = fima(z2);

    if(z.re > 0) return runloop(a, b, Math.floor(z.re), false);
    else return runloop(a, b, -Math.ceil(z.re), true);
  }

  // TODO: for complex z and arbitrary base: implement something like "kneser" function

  return Jmat.Complex(NaN);

  // TODO: implement super logarithm (slog), and super root (sroot) as well. Though, so far the only formulas I have is slog of base e in the Kouznetsov paper, and the 2-super root, so no implementation using both parameters can be done so far.
};


// Faddeeva function, used as helper functions to calculate erf and related functions for certain locations in the complex plane
// Faddeeva(z) = exp(-z^2)*erfc(-iz).
// Also known as Faddeyeva, or as w(x), but that may be confusing with LambertW...
Jmat.Complex.faddeeva = function(z) {
  // METOD A: series 7.1.8 from Handbook of Mathematical Functions
  // smaller area of convergence than METHOD A, so not used
  /*var result = Jmat.Complex(0);
  var zi = Jmat.Complex.I.mul(z);
  var zz = Jmat.Complex.ONE;
  for(var n = 0; n < 30; n++) {
    result = result.add(zz.divr(Jmat.Real.gamma(n/2 + 1)));
    zz = zz.mul(zi);
  }
  return result;*/

  var invsqrtpi2   = 2 / Jmat.Real.SQRTPI;
  var eye = (z.re * z.re) + (z.im * z.im * 2);

  // METHOD B: series
  // A small eye-shaped region in which the series works
  if(eye < 3.5) {
    // Based on sum 7.1.5 from Handbook of Mathematical Functions
    // erf(z) = 2/sqrt(pi) * SUM_n=0..oo (-1)^n * z^(2n+1) / (n! * (2n+1))
    // and then w(z) = e^(-z^2) * (1 - erf(iz))
    var sum = Jmat.Complex.ZERO;
    var sign = 1.0;
    var nn = 1;
    var iz = Jmat.Complex.I.mul(z).neg();
    var izz = iz;
    for(var n = 0; n < 20; n++) {
      if(n > 0) {
        nn = nn * n; // n!
        sign = -sign; // (-1)^n
        izz = izz.mul(iz).mul(iz); // iz^(2n+1)
      }
      sum = sum.add(izz.mulr(sign / (nn * (2*n + 1))));
    }
    var e = Jmat.Complex.exp(z.mul(z).neg());
    return e.sub(e.mul(sum).mulr(invsqrtpi2));
  }

  // METHOD C: Laplace Continued Fraction
  var za = Jmat.Complex(Math.abs(z.re), Math.abs(z.im)); // Operate on positive re, positive im quadrant
  // requires quite a lot of iterations unfortunately
  var num = eye < 40 ? 40 : eye < 80 ? 20 : 10;
  var result = Jmat.Complex(0);
  for(var n = 0; n < num; n++) {
    var r = Jmat.Complex(za.im + result.re, za.re - result.im);
    result = r.mulr(0.5 / r.abssqr());
  }
  result = result.mulr(invsqrtpi2);
  // Fix for pure imaginary values with large negative imaginary part
  if(za.im == 0.0) result.re = Math.exp(-za.re * za.re);
  // Put the solution back in the original quadrant, using the transformations w(-z) = 2 * exp(-z*z) - w(z) and w(conj(z)) = cons(w(-z))
  if(z.im < 0.0) {
    var e = Jmat.Complex.exp(za.mul(za).neg()).mulr(2);
    result = e.sub(result);
    if(z.re > 0.0) result.im = -result.im;
  } else if(z.re < 0.0) {
    result.im = -result.im;
  }
  return result;
};


// erfcx(z) = exp(z^2) * erfc(z): the scaled complementary error function
Jmat.Complex.erfcx = function(z) {
  return Jmat.Complex.faddeeva(Jmat.Complex(-z.im, z.re)); //erfcx(z) = faddeeva(iz)
};

Jmat.Complex.erf = function(z) {
  if(z.im == 0) {
    return Jmat.Complex(Jmat.Real.erf(z.re));
  } else if(z.re == 0) {
    return Jmat.Complex.I.mulr(Jmat.Real.erfi(z.im));
  } else {
    var a = Jmat.Complex.exp(z.mul(z).neg()); // If abs of z is very large, and |im| > |re|, then this becomes some NaN or Infinity. That is ok, erf is also some unrepresentable huge value there.
    if (z.re >= 0) return Jmat.Complex.ONE.sub(a.mul(Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I))));
    else return a.mul(Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I.neg()))).sub(Jmat.Complex.ONE);

    // With integration, don't use.
    /*var ps2 = 2.0 / Jmat.Real.SQRTPI;
    var result;
    result = Jmat.Complex.integrate(Jmat.Complex(0), z, function(z){ return Jmat.Complex.exp(z.mul(z).neg()); }, 100);
    result = result.mulr(ps2);
    return result;*/
  }
};

// erfc(x) = 1 - erf(x). This function gives numerically a better result if erf(x) is near 1.
Jmat.Complex.erfc = function(z) {
  if(z.im == 0) {
    return Jmat.Complex(Jmat.Real.erfc(z.re));
  } else {
    var a = Jmat.Complex.exp(z.mul(z).neg());
    if (z.re >= 0) return a.mul(Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I)));
    else return Jmat.Complex.TWO.sub(a.mul(Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I.neg()))));
  }
};


// TODO: rewrite some of the rational approximations to not use so many multiplications
//a + bx + cxx + dxxx + exxxx = a + x * (b + x * (c + x * (d + x * e)))   etc...
//and that equals: x.mulr(e).addr(d).mul(x).addr(c).mul(x).addr(b).mul(x).addr(a) ...


Jmat.Complex.erf_inv = function(z) {
  if (z.im != 0 && Math.abs(z.re) > 1) {
    //this branch is taken for large complex numbers because the implementation below doesn't work well on those. This one isn't much better btw, but slightly less bad for those cases.
    //TODO: is incorrect on many complex numbers!! e.g. erf_inv(erf(5 + 5i)) gives way wrong result
    //var zzpi = z.mul(z).mulr(Math.PI);
    //return zzpi.mulr(34807/182476800.0).addr(4369/5806080.0).mul(zzpi).addr(127/40320.0).mul(zzpi).addr(7/480.0).mul(zzpi).addr(1/12.0).mul(zzpi).addr(1).mul(z).mul(Jmat.Complex.SQRTPI).mulr(0.5);

    // With newton method

    // derivative of erf is: 2/sqrt(pi) * e^(-x^2)
    var derf = function(x) {
      return Jmat.Complex.TWO.divr(Jmat.Real.SQRTPI).mul(Jmat.Complex.exp(x.mul(x).neg()));
    };

    var neg = z.re < 0;
    if(neg) z = z.neg();

    // For abs(z) > 1 and z.re > 0, the following starting value works well: sqrt(-log(x * sqrt(pi) * (1 - x)))
    // NOTE: for abs(z) < 1, instead z*sqrtpi/2 would work, but other code already handles such case
    var start = Jmat.Complex.sqrt(Jmat.Complex.log(z.mulr(Jmat.Real.SQRTPI).mul(Jmat.Complex.ONE.sub(z))).neg());
    // NOTE: erf_inv has multiple solutions, with this starting value only one particular one is returned.
    // e.g. with the chosen starting value, erf_inv(2+2i) gives 0.386507600275 + 1.320769860731i. But with starting value 0, it gives the also correct 2.947736167125 + 3.401249486995i.
    var result = Jmat.Complex.finvert_newton(z, Jmat.Complex.erf, derf, start);
    if(neg) result = result.neg();
    return result;
  } else {
    //if (a > 1) return Jmat.Complex(NaN); //only relevant for real numbers
    if (z.im == 0) {
      if (z.re == 0) return Jmat.Complex(0);
      if (z.re == 1) return Jmat.Complex(Infinity);
      if (z.re == -1) return Jmat.Complex(-Infinity);
    }

    var erf_inv_a_ = [0.886226899, -1.645349621, 0.914624893, -0.140543331];
    var erf_inv_b_ = [1, -2.118377725, 1.442710462, -0.329097515, 0.012229801];
    var erf_inv_c_ = [-1.970840454, -1.62490649, 3.429567803, 1.641345311];
    var erf_inv_d_ = [1, 3.543889200, 1.637067800];

    var a = Jmat.Complex.abs(z).re;
    if (a <= 0.7) {
      var z2 = z.mul(z);
      var r = z.mul(z2.mulr(erf_inv_a_[3]).addr(erf_inv_a_[2]).mul(z2).addr(erf_inv_a_[1]).mul(z2).addr(erf_inv_a_[0]));
      r = r.div(z2.mulr(erf_inv_b_[4]).addr(erf_inv_b_[3]).mul(z2).addr(erf_inv_b_[2]).mul(z2).addr(erf_inv_b_[1]).mul(z2).addr(erf_inv_b_[0]));
    }
    else {
      var y = Jmat.Complex.sqrt(Jmat.Complex.log(Jmat.Complex.ONE.sub(z).divr(2)).neg());
      var r = y.mulr(erf_inv_c_[3]).addr(erf_inv_c_[2]).mul(y).addr(erf_inv_c_[1]).mul(y).addr(erf_inv_c_[0]);
      r = r.div(y.mulr(erf_inv_d_[2]).addr(erf_inv_d_[1]).mul(y).addr(erf_inv_d_[0]));
    }

    return r;
  }
};

// inverse complementary error function.
Jmat.Complex.erfc_inv = function(z) {
  return Jmat.Complex.erf_inv(Jmat.Complex.ONE.sub(z));
};

//erfi(z) = -i erf(iz)
Jmat.Complex.erfi = function(z) {
  if(Jmat.Complex.isReal(z)) return Jmat.Complex(Jmat.Real.erfi(z.re));
  return Jmat.Complex.erf(z.mul(Jmat.Complex.I)).mul(Jmat.Complex.I).neg();
};

// D+(x) aka F(x)
Jmat.Complex.dawson = function(z) {
  if(Jmat.Complex.isReal(z)) {
    return Jmat.Complex(Jmat.Real.dawson(z.re));
  } else {
    var w = Jmat.Complex.faddeeva(z);
    var a = Jmat.Complex.exp(z.mul(z).neg());
    return a.sub(w).mul(Jmat.Complex.I.mulr(Jmat.Real.SQRTPI / 2));
  }
};

//Minkowski's question mark function, from Wikipedia
Jmat.Complex.minkowski = function(z) {
  return Jmat.Complex(Jmat.Real.minkowski(z.re), Jmat.Real.minkowski(z.im));
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
Jmat.Complex.rootfind_bisection = function(valuex, valuey, f, maxit, prec) {
  var steps = (maxit == undefined ? 256 : maxit);
  var precision = (prec == undefined ? 0.000000001 : prec); //for 'near'
  var x = valuex.re;
  var y = valuey.re;
  if(y < x) {
    var temp = x;
    x = y;
    y = temp;
  }

  if(y - x == Infinity) return Jmat.Complex(NaN);

  //find a positive and a negative value
  var p = NaN;
  var n = NaN;
  var lastp = 0;
  var lastn = 0;
  for(var i = 0; i < steps; i++) {
    var z = x + i * (y - x) / steps; //exclude the leftmost point by adding precision. Otherwise annoying things can occur.
    var fz = f(Jmat.Complex(z)).re;
    if(Jmat.Real.near(fz, 0, precision)) return Jmat.Complex(z);

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
    return Jmat.Complex(NaN);
  }

  for(;;) {
    var z = (n + p) / 2;
    var fz = f(Jmat.Complex(z)).re;
    if(Jmat.Real.near(fz, 0, precision)) {
      return Jmat.Complex(z);
    }

    if(fz > 0) p = z;
    if(fz < 0) n = z;

    if(Jmat.Real.near(n, p, precision)) {
      return Jmat.Complex(NaN); //not found
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
Jmat.Complex.rootfind = function(f, o) {
  if(!o) o = {};
  var z0 = o.z0 || Jmat.Complex.ZERO;
  var maxit = o.maxit || 30;
  var prec = (o.prec == undefined ?  0.000000001 : o.prec);

  // TODO: work in progress. Find better start values. Try other root finding algorithms. Etc...
  if(o.real && o.z1 != undefined) return Jmat.Complex.rootfind_bisection(z0, o.z1, f, maxit, prec);
  if(o.df) return Jmat.Complex.rootfind_newton(f, o.df, z0, maxit);
  return Jmat.Complex.rootfind_secant(f, z0, maxit);
};

Jmat.Complex.newtonStartValues_ = [
    Jmat.Complex(0),
    Jmat.Complex(0.1), Jmat.Complex(-0.1), Jmat.Complex.newi(0.1), Jmat.Complex.newi(-0.1), 
    Jmat.Complex(0.1, 0.1), Jmat.Complex(0.1, -0.1), Jmat.Complex(-0.1, 0.1), Jmat.Complex(-0.1, -0.1), 
    Jmat.Complex(1), Jmat.Complex(-1), Jmat.Complex.newi(1), Jmat.Complex.newi(-1)
  ];

// This is really not that good. TODO: better root finding
Jmat.Complex.newtonStartValue_ = function(f) {
  var s = Jmat.Complex.newtonStartValues_;
  var bestdist = Infinity;
  var best = Jmat.Complex(NaN);

  for(var i = 0; i < s.length; i++) {
    var z = s[i];
    var fz = f(z);
    var m = Jmat.Complex.manhattan(fz, Jmat.Complex.ZERO);
    if(m < bestdist) {
      bestdist = m;
      best = z;
    }
  }

  return best;
};

Jmat.Complex.newtonStartValuesAround_ = [ Jmat.Complex(1), Jmat.Complex(-1), Jmat.Complex.newi(1), Jmat.Complex.newi(-1), ];

//Find a new value to start at after having accidentally encountered a bad point like NaN or Inf
//dist: e.g. 0.1 or 0.01
// This is really not that good. TODO: better root finding
Jmat.Complex.newtonStartValueAround_ = function(f, z0, dist) {
  var s = Jmat.Complex.newtonStartValuesAround_;
  var bestdist = Infinity;
  var best = Jmat.Complex(NaN);

  for(var i = 0; i < s.length; i++) {
    var z = z0.add(s[i].mulr(dist));
    var fz = f(z);
    var m = Jmat.Complex.manhattan(fz, Jmat.Complex.ZERO);
    if(m < bestdist) {
      bestdist = m;
      best = z;
    }
  }

  return best;
};

//Newton-Raphson. Finds a complex root (zero) given function f, its derivative df, and an initial value z0
Jmat.Complex.rootfind_newton = function(f, df, z0, maxiter) {
  if (!z0) z0 = Jmat.Complex.ZERO;//Jmat.Complex.newtonStartValue_(f);//Jmat.Complex.ZERO;
  if (!maxiter) maxiter = 30;
  var z = z0;
  var prevz = z;
  var bestdist = Infinity;
  var best = Jmat.Complex(NaN);
  for (var i = 0; i < maxiter; i++) {
    var fz = f(z);
    var m = Jmat.Complex.manhattan(fz, Jmat.Complex.ZERO);
    if(Jmat.Real.near(m, 0, 1e-15)) return z; // Near enough, stop iterations
    if(m < bestdist) {
      bestdist = m;
      best = z;
    }
    var zn = z.sub(fz.div(df(z)));
    if(Jmat.Complex.isInfOrNaN(zn)) {
      var m = Jmat.Complex.manhattan(z, prevz);
      zn = Jmat.Complex.newtonStartValueAround_(f, z, m ? m : 0.1);
    }
    z = zn;
    prevz = z;
  }

  return best;
};

//finds a complex root (zero) given function f, and an initial value z0 (no need to give the derivative)
Jmat.Complex.rootfind_secant = function(f, z0, maxiter) {
  return Jmat.Complex.rootfind_newton(f, function(x) {
    return Jmat.Complex.differentiate_stencil5(x, f);
  }, z0, maxiter);
};

//find result of inverse function using the newton method
Jmat.Complex.finvert_newton = function(z, f, df, z0, maxiter) {
  return Jmat.Complex.rootfind_newton(function(x) { return f(x).sub(z); }, df, z0, maxiter);
};

//find result of inverse function using the newton method (no need to give the derivative)
Jmat.Complex.finvert_secant = function(z, f, z0,  maxiter) {
  return Jmat.Complex.rootfind_secant(function(x) { return f(x).sub(z); }, z0, maxiter);
};


//numerical integration, aka quadrature
//integrate function f using simpsons rule, from x to y, with amount of steps given by steps parameter (higher = more precision, number of evaluations is thoroughly steps * 2)
//NOTE: this is the real version, real JS numbers only. Complex version is Jmat.Complex.integrate_simpson
Jmat.Real.integrate_simpson = function(x, y, steps, f, stopLoop) {
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
//NOTE: this is the real version, real JS numbers only. Complex version is Jmat.Complex.integrate
Jmat.Real.integrate = function(x, y, f, steps) {
  if(!steps) steps = 30;
  return Jmat.Real.integrate_simpson(x, y, steps, f);
};

//numerical integration, aka quadrature
//integrate function f using simpsons rule, from x to y, with amount of steps given by steps parameter (higher = more precision, number of evaluations is thoroughly steps * 2)
Jmat.Complex.integrate_simpson = function(x, y, steps, f, stopLoop) {
  var step = y.sub(x).divr(steps);
  var result = Jmat.Complex(0);
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

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return Jmat.Complex(NaN);
  }
  return result;
};

//numerical integration, aka quadrature
Jmat.Complex.integrate = function(x, y, f, steps) {
  if(!steps) steps = 30;
  return Jmat.Complex.integrate_simpson(x, y, steps, f);
};

// differentiation with just two points (finite difference, or secant)
Jmat.Complex.differentiate_secant = function(x, f) {
  //var h = Jmat.Complex(0.0001);
  var h = Math.max(0.01, Math.abs(x.re)) / 1000;

  var f1 = f(x.addr(h / 2));
  var f2 = f(x.subr(h / 2));
  return f1.sub(f2).divr(h);
};

// differentiation with 5-point stencil
// Ignores imaginary direction (e.g. derivative of the function Im(x) is always 0 with this. Im(x) is not holomorphic though), but that seems to be alright, Wolfram Alpha does that too...
Jmat.Complex.differentiate_stencil5 = function(x, f) {
  //var h = Jmat.Complex(0.0001);
  var h = Math.max(0.01, Math.abs(x.re)) / 1000;

  var f1 = f(x.addr(h * 2)).neg();
  var f2 = f(x.addr(h)).mulr(8);
  var f3 = f(x.subr(h)).mulr(-8);
  var f4 = f(x.subr(h * 2));
  return f1.add(f2).add(f3).add(f4).divr(h * 12);
};

// second derivative with 5-point stencil
Jmat.Complex.differentiate2nd_stencil5 = function(x, f) {
  //var h = Jmat.Complex(0.0001);
  var h = Jmat.Complex(Math.max(0.01, Math.abs(x.re)) / 1000);

  var f1 = f(x.add(h.mulr(2))).neg();
  var f2 = f(x.add(h)).mulr(16);
  var f3 = f(x).mulr(-30);
  var f4 = f(x.sub(h)).mulr(16);
  var f5 = f(x.sub(h.mulr(2))).neg();

  return f1.add(f2).add(f3).add(f4).add(f5).div(h.mul(h).mulr(12)); // (f1+f2+f3+f4+f5) / (h*h*12)
};

// differentiation (derivative)
Jmat.Complex.differentiate = function(x, f) {
  return Jmat.Complex.differentiate_stencil5(x, f);
};

//do summation of discrete values of f using start value x, end value y (inclusive) and given step.
Jmat.Complex.doSummation = function(valuex, valuey, step, f, stopLoop) {
  if(step == 0) return Jmat.Complex(NaN);
  if(step < 0) return Jmat.Complex(NaN);
  if(!Jmat.Complex.isReal(valuex) || !Jmat.Complex.isReal(valuey)) return Jmat.Complex(NaN);
  var x = valuex.re;
  var y = valuey.re;
  var result = Jmat.Complex(0);
  // the step / 4 thing is to avoid numerical problems that may let it miss the last value
  for(var z = x; z <= y + step / 4; z += step) {
    var fz = f(Jmat.Complex(z));
    result = result.add(fz);

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return Jmat.Complex(NaN);
  }
  return result;
};

//do product of discrete values of f using start value x, end value y (inclusive) and given step.
Jmat.Complex.doProduct = function(valuex, valuey, step, f, stopLoop) {
  if(step == 0) return Jmat.Complex(NaN);
  if(step < 0) return Jmat.Complex(NaN);
  if(!Jmat.Complex.isReal(valuex) || !Jmat.Complex.isReal(valuey)) return Jmat.Complex(NaN);
  var x = valuex.re;
  var y = valuey.re;
  var result = Jmat.Complex.ONE;
  // the step / 4 thing is to avoid numerical problems that may let it miss the last value
  for(var z = x; z <= y + step / 4; z += step) {
    var fz = f(Jmat.Complex(z));
    result = result.mul(fz);

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return Jmat.Complex(NaN);
  }
  return result;
};

// n is max coeff.length. Sum is from 0..n-1
Jmat.Complex.powerSeries = function(coeff, n, z0, z) {
  var realcoeff = coeff[0].re == undefined;
  var result = Jmat.Complex.ZERO;
  var zz = Jmat.Complex.ONE;
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

Jmat.Complex.pdf_uniform = function(x, a, b) {
  if(x.re >= a.re && x.re <= b.re) return b.sub(a).inv();
  return Jmat.Complex(0);
};

Jmat.Complex.cdf_uniform = function(x, a, b) {
  if(x.re < a.re) return Jmat.Complex(0);
  if(x.re < b.re) return x.sub(a).div(b.sub(a));
  return Jmat.Complex(1);
};

Jmat.Complex.qf_uniform = function(x, a, b) {
  var r = a.add(x.mul(b.sub(a)));
  if(r.re >= a.re && r.re <= b.re) return r;
  return Jmat.Complex(NaN);
};

////////////////////////////////////////////////////////////////////////////////

//aka small phi
Jmat.Complex.pdf_standardnormal = function(x) {
  return Jmat.Complex.exp(x.mul(x).mulr(-0.5)).mul(Jmat.Complex.INVSQRT2PI);
};

//aka capital PHI
Jmat.Complex.cdf_standardnormal = function(x) {
  return Jmat.Complex.erf(x.div(Jmat.Complex.SQRT2)).addr(1).mulr(0.5);
};

//aka the probit function
Jmat.Complex.qf_standardnormal = function(x) {
  return Jmat.Complex.erf_inv(x.mulr(2).subr(1)).mul(Jmat.Complex.SQRT2);
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pdf_normal = function(x, mu, sigma) {
  var a = Jmat.Complex.INVSQRT2PI.div(sigma);
  var b = x.sub(mu).mul(x.sub(mu)).div(sigma.mul(sigma).mulr(2)); // (x-mu)^2 / 2*sigma^2
  return a.mul(Jmat.Complex.exp(b.neg()));
};

Jmat.Complex.cdf_normal = function(x, mu, sigma) {
  var a = x.sub(mu).div(Jmat.Complex.abs(sigma).mul(Jmat.Complex.SQRT2)); // (x-mu) / sqrt(2*sigma^2)
  return Jmat.Complex.erf(a).addr(1).mulr(0.5);
};

Jmat.Complex.qf_normal = function(x, mu, sigma) {
  return Jmat.Complex.exp(mu.add(sigma.mul(Jmat.Complex.qf_standardnormal(x))));
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pdf_lognormal = function(x, mu, sigma) {
  // (1 / (x * sqrt(2*pi) * sigma)) * exp(- (ln(x) - mu)^2 / (2*sigma^2))
  var a = Jmat.Complex.INVSQRT2PI.div(sigma).div(x);
  var b = Jmat.Complex.log(x).sub(mu);
  return a.mul(Jmat.Complex.exp(b.mul(b).div(sigma.mul(sigma).mulr(2)).neg()));
};

Jmat.Complex.cdf_lognormal = function(x, mu, sigma) {
  var a = Jmat.Complex.log(x).sub(mu).div(Jmat.Complex.SQRT2.mul(sigma));
  return Jmat.Complex.erf(a).addr(1).mulr(0.5);
};

Jmat.Complex.qf_lognormal = function(x, mu, sigma) {
  var a = Jmat.Complex.log(x).sub(mu).div(Jmat.Complex.SQRT2.mul(sigma));
  return Jmat.Complex.erf(a).addr(1).mulr(0.5);
};

////////////////////////////////////////////////////////////////////////////////

//gamma = scale parameter (HWHM)
Jmat.Complex.pdf_cauchy = function(x, x0, gamma) {
  var x2 = x.sub(x0);
  var d = x2.mul(x2).add(gamma.mul(gamma)).mulr(Math.PI);
  return gamma.div(d);
};

Jmat.Complex.cdf_cauchy = function(x, x0, gamma) {
  return Jmat.Complex.atan(x.sub(x0).div(gamma)).divr(Math.PI).addr(0.5);
};

Jmat.Complex.qf_cauchy = function(x, x0, gamma) {
  return x0.add(gamma.mul(Jmat.Complex.tan(x.subr(0.5).divr(Math.PI))));
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pdf_studentt_cache_ = []; // cache used because nu is often constant between calls
Jmat.Complex.pdf_studentt_cachefun_ = function(nu) { var nu2 = nu.inc().divr(2); return Jmat.Complex.gammaDiv_(nu2, nu.divr(2)); };

// nu = degrees of freedom
Jmat.Complex.pdf_studentt = function(x, nu) {
  if(nu.eqr(1)) {
    return x.mul(x).addr(1).mulr(Math.PI).inv(); // cauchy distribution
  }
  if(nu.eqr(2)) {
    return x.mul(x).addr(2).powr(3/2).inv();
  }
  if(nu.eqr(3)) {
    var sqrt3x6 = 10.392304845413264; //6 * sqrt(3)
    var x2 = x.mul(x).addr(3);
    return Jmat.Complex(sqrt3x6).div(x2.mul(x2).mulr(Math.PI));
  }
  if(nu.eqr(Infinity)) {
    return Jmat.Complex.pdf_standardnormal(x);
  }

  var nu2 = nu.inc().divr(2);
  var g = Jmat.Complex.calcCache_(nu, Jmat.Complex.pdf_studentt_cachefun_, Jmat.Complex.pdf_studentt_cache_);
  var s = Jmat.Complex.sqrt(Jmat.Complex.PI.mul(nu)).inv();
  var gs = g.mul(s);
  if(Jmat.Complex.isNaN(gs) && nu.re > 100) gs = Jmat.Complex(Jmat.Complex.INVSQRT2PI); //this is what it is for nu = +Infinity
  
  var a = x.mul(x).div(nu).inc();
  return gs.mul(a.pow(nu2.neg()));
};

Jmat.Complex.cdf_studentt = function(x, nu) {
  if(nu.eqr(1)) {
    return Jmat.Complex.atan(x).divr(Math.PI).addr(0.5); // cauchy distribution
  }
  if(nu.eqr(2)) {
    return x.div(Jmat.Complex.sqrt(x.mul(x).addr(2)).mulr(2)).addr(0.5);
  }
  if(nu.eqr(Infinity)) {
    return Jmat.Complex.cdf_standardnormal(x);
  }

  // METHOD A: in terms of incomplete beta (probably slightly faster than with hypergeometric)
  if(x.eqr(0)) return Jmat.Complex(0.5);
  var nu2 = nu.inc().divr(2);
  var g = Jmat.Complex.calcCache_(nu, Jmat.Complex.pdf_studentt_cachefun_, Jmat.Complex.pdf_studentt_cache_);
  var b = Jmat.Complex.incbeta(x.mul(x).div(nu).neg(), Jmat.Complex(0.5), Jmat.Complex.ONE.sub(nu).mulr(0.5));
  var n = Jmat.Complex.I.mul(x).mul(b);
  var d = x.abs().mulr(2).mul(Jmat.Complex.SQRTPI);
  return Jmat.Complex(0.5).sub(g.mul(n).div(d));

  // METHOD B: in terms of hypergeometric
  /*var nu2 = nu.inc().divr(2);
  var g = Jmat.Complex.calcCache_(nu, Jmat.Complex.pdf_studentt_cachefun_, Jmat.Complex.pdf_studentt_cache_);
  var s = Jmat.Complex.sqrt(Jmat.Complex.PI.mul(nu)).inv();
  var gs = g.mul(s);
  if(Jmat.Complex.isNaN(gs) && nu.re > 100) gs = Jmat.Complex(Jmat.Complex.INVSQRT2PI); //this is what it is for nu = +Infinity
  var f = Jmat.Complex.hypergeometric(Jmat.Complex(0.5), nu.inc().mulr(0.5), Jmat.Complex(1.5), x.mul(x).div(nu).neg());
  return Jmat.Complex(0.5).add(x.mul(gs).mul(f));*/
};

//aka "tinv", t inverse cumulative distribution function
Jmat.Complex.qf_studentt= function(x, nu) {
  // test in console:
  // function testInvStudentt(x, nu) { return Jmat.Complex.qf_studentt(Jmat.Complex.cdf_studentt(Jmat.Complex(x), Jmat.Complex(nu)), Jmat.Complex(nu)).re; }
  // testInvStudentt(0.5,1.5)
  if(nu.eqr(1)) {
    return Jmat.Complex.tan(Jmat.Complex.PI.mul(x.subr(0.5))); // cauchy distribution
  }
  if(nu.eqr(2)) {
    var a = x.mulr(4).mul(Jmat.Complex.ONE.sub(x));
    return x.subr(0.5).mulr(2).mul(Jmat.Complex.sqrt(Jmat.Complex(2).div(a)));
  }
  if(nu.eqr(4) && Jmat.Complex.isReal(x)) {
    var a = x.mulr(4).mul(Jmat.Complex.ONE.sub(x));
    var as = Jmat.Complex.sqrt(a);
    var q = Jmat.Complex.cos(Jmat.Complex(1.0 / 3).mul(Jmat.Complex.acos(as))).div(as);
    return Jmat.Complex.sign(x.subr(0.5)).mulr(2).mul(Jmat.Complex.sqrt(q.dec()));
  }
  if(nu.eqr(Infinity)) {
    return Jmat.Complex.qf_standardnormal(x);
  }

  // Formula for student t inverse CDF in terms of inverse regularized beta function.
  // Works for real x in range 0-1. Extension to complex plane is fully of numerical problems here, so not supported.
  if(Jmat.Real.near(x.im, 0, 1e-15)) x = Jmat.Complex(x.re);
  if(Jmat.Real.near(nu.im, 0, 1e-15)) nu = Jmat.Complex(nu.re);
  if(Jmat.Complex.isPositive(nu) && Jmat.Complex.isPositive(x) && x.re < 1) {
    if(x.re < 0.5) {
      var i = Jmat.Complex.beta_i_inv(x.mulr(2), nu.divr(2), Jmat.Complex(0.5));
      return Jmat.Complex.sqrt(nu.mul(i.inv().subr(1))).neg();
    } else {
      var i = Jmat.Complex.beta_i_inv(Jmat.Complex.ONE.sub(x).mulr(2), nu.divr(2), Jmat.Complex(0.5));
      return Jmat.Complex.sqrt(nu.mul(i.inv().subr(1)));
    }
  }

  return Jmat.Complex(NaN); //TODO: support approximate inverse CDF for student t for arbitrary nu
};

////////////////////////////////////////////////////////////////////////////////

//k = degrees of freedom
Jmat.Complex.pdf_chi_square = function(x, k) {
  var k2 = k.divr(2);
  var a = k2.eqr(1) ? Jmat.Complex(2) : Jmat.Complex(2).pow(k2);
  var g = Jmat.Complex.isNegativeInt(k2) ? Jmat.Complex(Infinity) : Jmat.Complex.gamma(k2);
  var b = x.pow(k2.dec()).mul(Jmat.Complex.exp(x.divr(-2)));
  return a.mul(g).inv().mul(b);
};

Jmat.Complex.cdf_chi_square = function(x, k) {
  return Jmat.Complex.gamma_p(k.divr(2), x.divr(2));
};

// aka "cinv". x is in range 0-1
Jmat.Complex.qf_chi_square = function(x, k) {
  // Does NOT work with '0' as starting value - and not with '0.5' either. '1' seems to be the best for the range of operation (x = 0..1)
  // TODO: not very precise for higher k such as k=5. Fix this.
  // E.g. should be: k=1,x=0.4: 0.274997, k=0.5&x=0.5: 0.087347000000, k=5&x=0.8:7.28928

  return Jmat.Complex.gamma_p_inv(k.divr(2), x).mulr(2);
};

////////////////////////////////////////////////////////////////////////////////

// mu = location, s = scale
Jmat.Complex.pdf_logistic = function(x, mu, s) {
  var e = Jmat.Complex.exp(x.sub(mu).div(s).neg());
  var ee = e.inc();
  return e.div(s.mul(ee).mul(ee));
};

Jmat.Complex.cdf_logistic = function(x, mu, s) {
  return Jmat.Complex.tanh(x.sub(mu).div(s).divr(2)).mulr(1/2).addr(1/2);
};

// Note: the *derivative* of this is the logit function s/(x*(1-x))
Jmat.Complex.qf_logistic = function(x, mu, s) {
  var xx = x.div(Jmat.Complex.ONE.sub(x));
  return mu.add(s.mul(Jmat.Complex.log(xx)));
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pdf_gamma_cache_ = [];

// k = shape, theta = scale
Jmat.Complex.pdf_gamma = function(x, k, theta) {
  var xk = x.pow(k.dec());
  var e = Jmat.Complex.exp(x.div(theta).neg());
  var t = theta.pow(k);
  var g = Jmat.Complex.calcCache_(k, Jmat.Complex.gamma, Jmat.Complex.pdf_gamma_cache_); // gamma(k)
  return xk.mul(e).div(t).div(g);
};

Jmat.Complex.cdf_gamma = function(x, k, theta) {
  return Jmat.Complex.gamma_p(k, x.div(theta));
};

Jmat.Complex.qf_gamma = function(x, k, theta) {
  return Jmat.Complex.gamma_p_inv(k, x).mul(theta);
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pdf_beta = function(x, alpha, beta) {
  var xa = x.pow(alpha.dec());
  var xb = x.rsub(1).pow(beta.dec());
  var b = Jmat.Complex.beta(alpha, beta);
  return xa.mul(xb).div(b);
};

Jmat.Complex.cdf_beta = function(x, alpha, beta) {
  return Jmat.Complex.beta_i(x, alpha, beta);
};

Jmat.Complex.qf_beta = function(x, alpha, beta) {
  return Jmat.Complex.beta_i_inv(x, alpha, beta);
};

////////////////////////////////////////////////////////////////////////////////

// Fisher–Snedecor F distribution 
// d1 and d2 are the degrees of freedom parameters
// TODO: gives NaN for x = 0, it probably should give 0 instead?
Jmat.Complex.pdf_fisher = function(x, d1, d2) {
  var a = d1.mul(x).pow(d1).mul(d2.pow(d2));
  var b = d1.mul(x).add(d2).pow(d1.add(d2));
  var c = x.mul(Jmat.Complex.beta(d1.divr(2), d2.divr(2)));
  return Jmat.Complex.sqrt(a.div(b)).div(c);
};

Jmat.Complex.cdf_fisher = function(x, d1, d2) {
  var a = d1.mul(x).div(d1.mul(x).add(d2));
  return Jmat.Complex.beta_i(a, d1.divr(2), d2.divr(2));
};

// aka "finv"
Jmat.Complex.qf_fisher = function(x, d1, d2) {
  var a = Jmat.Complex.beta_i_inv(x, d1.divr(2), d2.divr(2));
  var b = Jmat.Complex.ONE.sub(a);
  return d2.mul(a).div(d1.mul(b));
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pdf_weibull = function(x, lambda, k) {
  var a = k.div(lambda);
  var b = (x.div(lambda)).pow(k.dec());
  var c = x.div(lambda).pow(k).neg();
  return a.mul(b).mul(Jmat.Complex.exp(c));
};

Jmat.Complex.cdf_weibull = function(x, lambda, k) {
  var c = x.div(lambda).pow(k).neg();
  return Jmat.Complex.ONE.sub(Jmat.Complex.exp(c));
};

Jmat.Complex.qf_weibull = function(x, lambda, k) {
  var a = Jmat.Complex.ONE.sub(x);
  var b = Jmat.Complex.log(a).neg().pow(k.inv());
  return lambda.mul(b);
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.pdf_exponential = function(x, lambda) {
  if(x.re < 0) return Jmat.Complex(0);
  return lambda.mul(Jmat.Complex.exp(x.mul(lambda).neg()));
};

Jmat.Complex.cdf_exponential = function(x, lambda) {
  if(x.re < 0) return Jmat.Complex(0);
  return Jmat.Complex.ONE.sub(Jmat.Complex.exp(x.mul(lambda).neg()));
};

Jmat.Complex.qf_exponential = function(x, lambda) {
  return Jmat.Complex.log(x.rsub(1)).neg().div(lambda);
};

////////////////////////////////////////////////////////////////////////////////

// aka "double exponential"
Jmat.Complex.pdf_laplace = function(x, mu, b) {
  // 1/2b * exp(-|x-mu\/b)
  var e = Jmat.Complex.exp(Jmat.Complex.abs(x.sub(mu)).neg().div(b));
  return b.mulr(2).inv().mul(e);
};

Jmat.Complex.cdf_laplace = function(x, mu, b) {
  var e = Jmat.Complex.exp(Jmat.Complex.abs(x.sub(mu)).neg().div(b));
  var s = Jmat.Complex.sign(x.sub(mu));
  return e.rsub(1).mul(s).mulr(0.5).addr(0.5);
};

Jmat.Complex.qf_laplace = function(x, mu, b) {
  var l = Jmat.Complex.log(Jmat.Complex.abs(x.subr(0.5)).mulr(2).rsub(1));
  var s = Jmat.Complex.sign(x.subr(0.5));
  return mu.sub(b.mul(s).mul(l));
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Jmat.Matrix: Matrix math
////////////////////////////////////////////////////////////////////////////////

/*
Constructor
height first because that's the math convention: a 2x3 matrix is 2 rows high, 3 columns wide, and made as new Jmat.Matrix(2, 3).
Does NOT initialize elements if keyword "new" is in front. If keyword "new" is _not_ in front, then uses Jmat.Matrix.make and all its options to initialize elements instead.

Aliased as simply "Matrix" at the end of the file - disable that if it causes name clashes

This is the actual object used as matrix. In addition, most of the matrix
functions are implemented as static functions in here.
*/
Jmat.Matrix = function(height, width, var_arg) {
  if(this instanceof Jmat.Matrix) {
    // Keyword "new" in front, use basic constructor
    this.h = height; //number of rows
    this.w = width; //number of columns
    this.e = []; //array of arrays. first index is row (y), second index is column (x)
    for(var y = 0; y < height; y++) {
      this.e[y] = []; // Make the rows available, but elements are NOT initialized. They must be initialized after using the "new" constructor.
    }
  } else {
    // No keyword "new" in front, use the convenience factory function instead
    return Jmat.Matrix.make.apply(this, arguments); // pass on the variable length arguments with apply. Note that "this" is not the Matrix here, but a browser Window or so.
  }
};

/*
Makes a matrix from many types of combinations of arguments
numerical values can either be JS numbers, or Jmat.Complex numbers
Here are all the combinations:

*) a = integer, b = integer: a*b matrix of _uninitialized_ elements (undefined, similar to the new Matrix constructor)
Jmat.Matrix(2, 2).toString()               --> [[undefined, undefined], [undefined, undefined]]

*) a = string: matrix parsed from string (if valid) - the string uses square brackets around the whole matrix and around each row, with rows and elements separated by commas.
Jmat.Matrix('[[1, 1e2], [3, 4i]]').toString()         --> [[1, 100], [3, 4i]]

*) a = 1D/2D array of numerical values: column vector or 2D matrix
Jmat.Matrix([[1, 2], [3, 4]]).toString()   --> [[1, 2], [3, 4]]
Jmat.Matrix([1, 2, 3, 4]).toString()       --> [[1],[2],[3],[4]]: a column matrix
Jmat.Matrix([[1, 2, 3, 4]]).toString()     --> [[1, 2, 3, 4]]: a row matrix
Jmat.Matrix([[Jmat.Complex(1, 5), Jmat.Complex(2, 6)], [Jmat.Complex(3, 7), Jmat.Complex(4, 8)]]).toString() --> [[1+5i, 2+6i], [3+7i, 4+8i]]

*) a = 1D/2D array or undefined, b = 1D/2D array: column vector or 2D matrix with complex values made from real parts of a, imaginary parts of b
Jmat.Matrix([[1, 2], [3, 4]], [[5, 6], [7, 8]]).toString() --> [[1+5i, 2+6i], [3+7i, 4+8i]]
Jmat.Matrix(undefined, [[1, 2], [3, 4]]).toString()        --> [[1i, 2i], [3i, 4i]]

*) a and b positive integers, var_arg = 1D array or implicit 1D array in arguments ...
**) ... with a*b elements: 2D matrix with elements from var_arg, height a, width b
Jmat.Matrix(2, 2, 1, 2, 3, 4).toString()   --> [[1, 2], [3, 4]]
**) ... with size of array == min(a,b), creates diagonal matrix from those elements
Jmat.Matrix(4, 4, 1, 2, 3, 4).toString()   --> [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0], [0, 0, 0, 4]]: diagonal matrix
**) ... with size of array == 1, creates diagonal matrix with that one element repeated
Jmat.Matrix(3, 3, 0).toString()            --> [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
Jmat.Matrix(3, 3, 1).toString()            --> [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
**) ... with any other size: not allowed, returns invalid
Jmat.Matrix(3, 3, 1, 2)                    --> null

*) a = numerical value: 1x1 matrix with that element
Jmat.Matrix(0).toString()                  --> [0]

*) a = a Jmat.Matrix object: copy the matrix to a new Jmat.Matrix
Jmat.Matrix(Jmat.Matrix(0)).toString()     --> [0]

TODO: parse function for matrices
*/
Jmat.Matrix.make = function(a, b, var_arg) {
  if(a instanceof Jmat.Matrix) return Jmat.Matrix.copy(a);

  if(typeof a == 'string') return Jmat.Matrix.parse(a);

  // Tolerant to all kinds of nonexisting array
  // Also supports a 1D array representing an Nx1 2D array
  var softget = function(a, y, x) {
    return (a && a[y]) ? Jmat.Complex.cast(a[y][x] == undefined ? a[y] : a[y][x]) : Jmat.Complex();
  };
  var softgetr = function(a, y, x) {
    return (a && a[y]) ? Jmat.Real.cast(a[y][x] == undefined ? a[y] : a[y][x]) : 0;
  };
  var softget2 = function(a, b, y, x) {
    return new Jmat.Complex(softgetr(a, y, x), softgetr(b, y, x)); // real from a, imag from b
  };
  var arrayw = function(a) {
    if(!a || !a[0]) return 0; // empty array
    if(a[0].length == undefined) return 1; // this means it's a 1D array, such array has width 1, not width 0
    return a[0].length;
  };
  // a is array, b is optional second array (only supported for shift == -1)
  // shift -1 ==> a is 2D array
  // shift >= 0 ==> a is 1D array of which the 2D elements are read, first element at a[shift]
  // shift >= 0 and a.length - shift is min(h, w) ==> makes diagonal matrix with those elements on the diagonal instead
  var loop = function(h, w, a, shift, opt_b) {
    var result;
    if(shift >= 0 && a.length < h * w) {
      result = Jmat.Matrix.zero(h, w);
      if(shift >= 0 && a.length - shift == 1) {
        // repeat same element on diagonal
        for(var x = 0; x < w && x < h; x++) {
          result.e[x][x] = Jmat.Complex.cast(a[shift]);
        }
        return result;
      }
      if(shift >= 0 && a.length - shift == Math.min(h, w)) {
        // fill diagonal
        for(var x = 0; x < w && x < h; x++) {
          result.e[x][x] = Jmat.Complex.cast(a[x + shift]);
        }
        return result;
      }
      return null; //invalid size
    }
    // full 2D matrix
    result = new Jmat.Matrix(h, w);
    for(var y = 0; y < result.h; y++) {
      for(var x = 0; x < result.w; x++) {
        if(shift < 0) result.e[y][x] = (opt_b ? softget2(a, opt_b, y, x) : softget(a, y, x));
        else result.e[y][x] = Jmat.Complex.cast(a[y * w + x + shift]);
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

  // single number or Jmat.Complex
  if(a != undefined && b == undefined) {
    var result = new Jmat.Matrix(1, 1);
    result.e[0][0] = Jmat.Complex.cast(a);
    return result;
  }

  // a and b contain dimensions, then elements in array or variable arguments
  if(a != undefined && b != undefined) {
    if(var_arg == undefined) return new Jmat.Matrix(a, b); // simple constructor with uninitialized elements
    var h = a;
    var w = b;
    if(var_arg && var_arg.length) return loop(h, w, var_arg, 0);
    return loop(h, w, arguments, 2); // use JS function arguments, shifted by 2 because the first two are a and b
  }

  return new Jmat.Matrix(0, 0);
};

//debugstring
Jmat.Matrix.toString = function(m) {
  // e.g in console: Jmat.Matrix.toString(getMatrixFromMem(Jmat.Complex(100))
  if(!m) return '' + m;
  var s = '[';
  for(var y = 0; y < m.h; y++) {
    s += '[';
    for(var x = 0; x < m.w; x++) {
      var e = m.e[y][x];
      s += Jmat.Complex.toString(e);
      if(x + 1 < m.w) s += ', ';
    }
    s += ']';
    if(y + 1 < m.h) s += ', ';
  }
  return s + ']';
};
Jmat.Matrix.prototype.toString = function() {
  return Jmat.Matrix.toString(this);
};

Jmat.Matrix.parse = function(text) {
  var e = [];
  var stack = [e];
  var text2 = '';
  for(var i = 1; i < text.length - 1 && stack.length > 0; i++) {
    var c = text.charAt(i);
    if(c == '[') {
      var a = [];
      stack[stack.length - 1].push(a);
      stack.push(a);
    } else if(c == ']') {
      if(text2 != '') {
        stack[stack.length - 1].push(Jmat.Complex.parse(text2));
        text2 = '';
      }
      stack.pop();
    } else if(c == ',') {
      if(text2 != '') {
        stack[stack.length - 1].push(Jmat.Complex.parse(text2));
        text2 = '';
      }
    } else {
      text2 += c;
    }
  }
  if(text2 != '') stack[stack.length - 1].push(Jmat.Complex.parse(text2));
  return Jmat.Matrix.make(e);
};


// Does not copy if a is of type Jmat.Matrix.
Jmat.Matrix.cast = function(a) {
  return a instanceof Jmat.Matrix ? a : Jmat.Matrix.make(a);
};

Jmat.Matrix.copy = function(a) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = Jmat.Complex(a.e[y][x].re, a.e[y][x].im);
    }
  }
  return result;
};

// Returns new h*w identity matrix
Jmat.Matrix.identity = function(h, w) {
  var r = new Jmat.Matrix(h, w);
  for(var y = 0; y < h; y++) {
    for(var x = 0; x < w; x++) {
      r.e[y][x] = Jmat.Complex(x == y ? 1 : 0);
    }
  }
  return r;
};

// Returns new h*w zero matrix
Jmat.Matrix.zero = function(h, w) {
  var r = new Jmat.Matrix(h, w);
  for(var y = 0; y < h; y++) {
    for(var x = 0; x < w; x++) {
      r.e[y][x] = Jmat.Complex(0);
    }
  }
  return r;
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Matrix.add = function(a, b) {
  if(a.w != b.w || a.h != b.h) return null;
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].add(b.e[y][x]);
    }
  }
  return result;
};
Jmat.Matrix.prototype.add = function(b) {
  return Jmat.Matrix.add(this, b);
};

Jmat.Matrix.sub = function(a, b) {
  if(a.w != b.w || a.h != b.h) return null;
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].sub(b.e[y][x]);
    }
  }
  return result;
};
Jmat.Matrix.prototype.sub = function(b) {
  return Jmat.Matrix.sub(this, b);
};

Jmat.Matrix.mul = function(a, b) {
  if(a.w != b.h) return null;
  var result = new Jmat.Matrix(a.h, b.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < b.w; x++) {
      var e = Jmat.Complex(0);
      for(var z = 0; z < a.w; z++) e = e.add(a.e[y][z].mul(b.e[z][x]));
      result.e[y][x] = e;
    }
  }
  return result;
};
Jmat.Matrix.prototype.mul = function(b) {
  return Jmat.Matrix.mul(this, b);
};

// mulScalar (c from complex number)
Jmat.Matrix.mulc = function(a, s) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].mul(s);
    }
  }
  return result;
};
Jmat.Matrix.prototype.mulc = function(s) {
  return Jmat.Matrix.mulc(this, s);
};

Jmat.Matrix.mulr = function(a, s) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].mulr(s);
    }
  }
  return result;
};
Jmat.Matrix.prototype.mulr = function(s) {
  return Jmat.Matrix.mulr(this, s);
};

// returns a/b = a * b^-1
// In other words, "divides" matrix through matrix
Jmat.Matrix.div = function(a, b) {
  if(a.w != b.h) return null;
  var result = new Jmat.Matrix(a.h, b.w);

  b = Jmat.Matrix.inv(b); //TODO: use pseudo inverse?

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < b.w; x++) {
      var e = Jmat.Complex(0);
      for(var z = 0; z < a.w; z++) e = e.add(a.e[y][z].mul(b.e[z][x]));
      result.e[y][x] = e;
    }
  }
  return result;
};
Jmat.Matrix.prototype.div = function(b) {
  return Jmat.Matrix.div(this, b);
};

// returns a/b = b^-1 * a
Jmat.Matrix.leftdiv = function(a, b) {
  if(a.w != b.h) return null;
  var result = new Jmat.Matrix(a.h, b.w);

  b = Jmat.Matrix.inv(b); //TODO: use pseudo inverse?

  for(var y = 0; y < b.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = Jmat.Complex(0);
      for(var z = 0; z < b.w; z++) e = e.add(b.e[y][z].mul(a.e[z][x]));
      result.e[y][x] = e;
    }
  }
  return result;
};
Jmat.Matrix.prototype.leftdiv = function(b) {
  return Jmat.Matrix.leftdiv(this, b);
};

// Divide through complex scalar
Jmat.Matrix.divc = function(a, s) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].div(s);
    }
  }
  return result;
};
Jmat.Matrix.prototype.divc = function(s) {
  return Jmat.Matrix.divc(this, s);
};

// divide through real scalar (JS number)
Jmat.Matrix.divr = function(a, s) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].divr(s);
    }
  }
  return result;
};
Jmat.Matrix.prototype.divr = function(s) {
  return Jmat.Matrix.divr(this, s);
};

////////////////////////////////////////////////////////////////////////////////
// Categories
////////////////////////////////////////////////////////////////////////////////

Jmat.Matrix.isReal = function(a) {
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = a.e[y][x];
      if(e.im != 0) return false;
    }
  }
  return true;
};

//valid object, no infinities, no NaN elements
Jmat.Matrix.isValid = function(a) {
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

//returns true if any NaN in matrix. For the rest, must be valid object.
Jmat.Matrix.isNaN = function(a) {
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if(Jmat.Complex.isNaN(a.e[y][x])) return true;
    }
  }
  return false;
};

// TODO: functions like isSymmetrical, isHermitian, isDiagonal, ...

////////////////////////////////////////////////////////////////////////////////

Jmat.Matrix.transpose = function(a) {
  var result = new Jmat.Matrix(a.w, a.h); //arguments inverted (ctor takes height first)

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[x][y] = a.e[y][x];
    }
  }
  return result;
};
Jmat.Matrix.prototype.transpose = function() {
  return Jmat.Matrix.transpose(this);
};

Jmat.Matrix.neg = function(a) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].neg();
    }
  }
  return result;
};
Jmat.Matrix.prototype.neg = function() {
  return Jmat.Matrix.neg(this);
};

Jmat.Matrix.conj = function(a) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = Jmat.Complex.conj(a.e[y][x]);
    }
  }
  return result;
};
Jmat.Matrix.prototype.conj = function() {
  return Jmat.Matrix.conj(this);
};

//transjugate = transposed conjugate, denoted A* (also called hermitian transpose)
Jmat.Matrix.transjugate = function(a) {
  return Jmat.Matrix.conj(Jmat.Matrix.transpose(a));
};
Jmat.Matrix.prototype.transjugate = function() {
  return Jmat.Matrix.transjugate(this);
};

//TODO: LU Decomposition. Implement LUP or similar, regular LU would not work for [0 1][2 3]

// Submatrix with 1 row removed
Jmat.Matrix.subrow = function(a, row) {
  if(a.h < 2) return null;
  var m = new Jmat.Matrix(a.h - 1, a.w);
  for(var y = 0; y < a.h - 1; y++) {
    for(var x = 0; x < a.w; x++) {
      m.e[y][x] = a.e[y < row ? y : y + 1][x];
    }
  }
  return m;
};

// Submatrix with 1 column removed
Jmat.Matrix.subcol = function(a, col) {
  if(a.w < 2) return null;
  var m = new Jmat.Matrix(a.h, a.w - 1);
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w - 1; x++) {
      m.e[y][x] = a.e[y][x < col ? x : x + 1];
    }
  }
  return m;
};

// The submatrix for a minor, that is, with one row and one column removed
Jmat.Matrix.minorsub = function(a, row, col) {
  if(a.h < 2 || a.w < 2) return null;
  var m = new Jmat.Matrix(a.h - 1, a.w - 1);
  for(var y = 0; y < a.h - 1; y++) {
    for(var x = 0; x < a.w - 1; x++) {
      m.e[y][x] = a.e[y < row ? y : y + 1][x < col ? x : x + 1];
    }
  }
  return m;
};

//submatrix defined by the rectangle y0:y1 x0:x1 (excluding y1 and x1)
//e.g. to get the bottom right 2x2 submatrix of a 5x5 matrix a, use:
//Jmat.Matrix.submatrix(a, 3, 5, 3, 5)
//note that x and y are 0-indexed, so 5 is outside the matrix
Jmat.Matrix.submatrix = function(a, y0, y1, x0, x1) {
  if(x0 < 0 || y0 < 0 || x0 > a.w || y0 > a.h) return null;
  if(x1 < 0 || y1 < 0 || x1 > a.w || y1 > a.h) return null;
  var w2 = x1 - x0;
  var h2 = y1 - y0;
  if(w2 <= 0 || h2 <= 0 || w2 > a.w || h2 > a.h) return null;

  var result = new Jmat.Matrix(h2, w2);

  for(var y = 0; y < h2; y++) {
    for(var x = 0; x < w2; x++) {
      result.e[y][x] = a.e[y0 + y][x0 + x];
    }
  }

  return result;
};

// Requires square matrix
Jmat.Matrix.minor = function(a, row, col) {
  if(a.h < 2 || a.w < 2 || a.w != a.h) return Jmat.Complex(NaN);
  var m = Jmat.Matrix.minorsub(a, row, col);
  return Jmat.Matrix.determinant(m);
};

// cofactor: minor with sign depending on alternating position
Jmat.Matrix.cofactor = function(a, row, col) {
  var m = Jmat.Matrix.minor(a, row, col);
  var sign = (row + col) % 2 == 0 ? 1 : -1;
  return m.mulr(sign);
};

Jmat.Matrix.determinant = function(a) {
  if(a.w != a.h) return NaN; //square matrices only

  if(a.w == 1) return a.e[0][0];
  if(a.w == 2) return a.e[0][0].mul(a.e[1][1]).sub(a.e[0][1].mul(a.e[1][0]));

  var result = Jmat.Complex(0);

  for(var x = 0; x < a.w; x++) {
    result = result.add(a.e[0][x].mul(Jmat.Matrix.cofactor(a, 0, x)));
  }

  return result;
};

//Adjugate aka Adjoint matrix
Jmat.Matrix.adj = function(a) {
  if(a.w != a.h) return NaN; //square matrices only
  if(a.w == 1) return Jmat.Matrix.identity(1, 1);

  //result matrix
  var r = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      //row and column are switched (transpose)
      r.e[y][x] = Jmat.Matrix.cofactor(a, x, y);
    }
  }

  return r;
};

// Inverse of a matrix
Jmat.Matrix.inv = function(a) {
  if(a.w != a.h) return null; //square matrices only

  //Cramer's rule
  return Jmat.Matrix.mulc(Jmat.Matrix.adj(a), Jmat.Complex.inv(Jmat.Matrix.determinant(a)));
};

//forced pseudoinverse (does not try regular inverse first)
Jmat.Matrix.pseudoinverse_ = function(a) {
  //TODO: instead of formula below, easier ways are available for scalars, vectors, lin. indep. rows/columns, ...

  var svd = Jmat.Matrix.svd(a);
  var n = Math.min(svd.s.w, svd.s.h);
  var tolerance = 1e-15; // choose this such that e.g. the pseudoinverse of [[1,2],[1,2]] returns [[0.1,0.2],[0.1,0.2]] - if the tolerance is too small, some near zero gets inverted and the result very wrong (depends on the svd algorithm used of course)
  // Invert all the elements of s, except those that are zero (with some tolerance due to numerical problems)
  // Each element of s should be real, and it only has elements on the diagonal
  for(var i = 0; i < n; i++) svd.s.e[i][i] = (Math.abs(svd.s.e[i][i].re) < tolerance) ? svd.s.e[i][i] : Jmat.Complex.inv(svd.s.e[i][i]);
  svd.s = Jmat.Matrix.transpose(svd.s);

  return Jmat.Matrix.mul(Jmat.Matrix.mul(svd.v, svd.s), Jmat.Matrix.transjugate(svd.u));
};

//Moore-Penrose pseudoinverse: one unique solution for any matrix
Jmat.Matrix.pseudoinverse = function(a) {
  /*
  Test in console:
  var result = Jmat.Matrix.pseudoinverse(Jmat.Matrix.make(2,2,1,2,1,2));
  var orig = Jmat.Matrix.pseudoinverse(result);
  result.toString() + ' | ' + orig.toString();
  */

  //first try if regular inverse works, if so that is more accurate and faster
  var result = Jmat.Matrix.inv(a);
  if(Jmat.Matrix.isValid(result)) return result;

  return Jmat.Matrix.pseudoinverse_(a);
};

Jmat.Matrix.getFirstNonZeroDigit_ = function(v) {
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
Jmat.Matrix.getDebugNumber = function(a) {
  var result = a.w + 10*a.h;
  var pos = 0.1;
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var v = Jmat.Matrix.getFirstNonZeroDigit_(a.e[y][x].re);
      result += v * pos;
      pos /= 10;
    }
  }

  // fix some numerical problem
  pos *= 10;
  result = Math.round(result / pos) * pos;

  return Jmat.Complex(result);
};

/*
Matrix norms:
This list is shown here to ensure to not confuse the Frobenius norm with the 2-norm
oo = infinity

p-norms (a.k.a. induced norms or operator norms)
-------
1-norm: maximum absolute column sum of the matrix --> Jmat.maxcolnorm
2-norm: largest singular value, aka spectral norm --> Jmat.Matrix.norm2
oo-norm: maximum absolute row sum of the matrix --> Jmat.Matrix.maxrownorm

entrywise norms
---------------
entrywise 1-norm: sum of abs of all the elements --> (not implemented yet)
entrywise 2-norm: Frobenius norm, sqrt of sum of squares of the elements --> Jmat.Matrix.norm
entrywise oo-norm: maximum of abs of all the elements --> (not implemented yet)

schatten norms
--------------
schatten 1-norm: sum of singular values, aka nuclear norm, trace norm or Ky Fan norm --> (not implemented yet)
schatten 2-norm: sqrt of sum of squares of singular values, results in same value as Frobenius norm --> Jmat.Matrix.norm
schatten oo-norm: max of the singular values, results in same value as the spectral norm (2-norm) --> Jmat.Matrix.norm2

*/

//Frobenius norm of the matrix (sqrt of sum of squares of modulus of all elements)
//For a vector, this is the Euclidean norm.
//TODO: since usually the more expensive to calculate 2-norm is meant by "the" norm of the matrix, maybe rename this function to "frobeniusnorm" or "frob"?
Jmat.Matrix.norm = function(m) {
  var result = 0;
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      result += Jmat.Complex.abssqr(m.e[y][x]);
    }
  }
  result = Math.sqrt(result);
  return Jmat.Complex(result);
};

//Maximum absolute column sum norm
Jmat.Matrix.maxcolnorm = function(m) {
  var result = 0;
  for(var x = 0; x < m.w; x++) {
    var current = 0;
    for(var y = 0; y < m.h; y++) {
      current += Jmat.Complex.abssqr(m.e[y][x]);
    }
    if (current > result) result = current;
  }
  return Jmat.Complex(result);
};

//Maximum absolute row sum norm
Jmat.Matrix.maxrownorm = function(m) {
  var result = 0;
  for(var y = 0; y < m.h; y++) {
    var current = 0;
    for(var x = 0; x < m.w; x++) {
      current += Jmat.Complex.abssqr(m.e[y][x]);
    }
    if (current > result) result = current;
  }
  return Jmat.Complex(result);
};

//2-norm: largest singular value (= sqrt of largest eigenvalue of m^H * m). AKA spectral norm
Jmat.Matrix.norm2 = function(m) {
  var svd = Jmat.Matrix.svd(m);
  return svd.s.e[0][0];
};

//condition number: largest singular value divided through smallest singular value. Higher ==> more imprecise numerical calculations with this matrix
Jmat.Matrix.conditionNumber = function(m) {
  var svd = Jmat.Matrix.svd(m);
  var d = Math.min(m.w, m.h);
  return svd.s.e[0][0].div(svd.s.e[d - 1][d - 1]);
};

//Rank of matrix
Jmat.Matrix.rank = function(m) {
  // TODO: use the faster RRQR? Or at least svd that returns only s?
  var s = Jmat.Matrix.svd(m).s;
  var rank = 0;
  for(var i = 0; i < s.w; i++) {
    if(!Jmat.Real.near(s.e[i][i].re, 0, 1e-14)) rank++;
  }
  return Jmat.Complex(rank);
};

// Only defined for square matrices
Jmat.Matrix.trace = function(m) {
  if(m.w != m.h) return Jmat.Complex(NaN);
  var result = Jmat.Complex.ZERO;
  for(var x = 0; x < m.w; x++) result = result.add(m.e[x][x]);
  return result;
};

// Returns column as h*1 matrix
Jmat.Matrix.col = function(m, col) {
  var r = new Jmat.Matrix(m.h, 1);
  for(var y = 0; y < m.h; y++) r.e[y][0] = m.e[y][col];
  return r;
};

// Returns row as 1*w matrix
Jmat.Matrix.row = function(m, row) {
  var r = new Jmat.Matrix(1, m.w);
  for(var x = 0; x < m.w; x++) r.e[0][x] = m.e[row][x];
  return r;
};

// Add two non-equal sized matrices.
// b's top left element is at position (row, col) in a (0-indexed)
// so the size of the result matrix is:
// max(row + b.h, a.h) - min(0, row) by max(col + b.w, a.w) - min(0, col)
// any element not overlapped by a or b, will be zero.
Jmat.Matrix.overlap = function(a, b, row, col) {
  var h = Math.max(row + b.h, a.h) - Math.min(0, row);
  var w = Math.max(col + b.w, a.w) - Math.min(0, col);

  var result = new Jmat.Matrix(h, w);

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
Jmat.Matrix.qr = function(m) {
  /*
  Checks in console:
  var result = Jmat.Matrix.qr(Jmat.Matrix(2,2,1,2,3,4));
  Jmat.Matrix.toString(result.q) + ' \n ' + Jmat.Matrix.toString(result.r) + ' \n ' + Jmat.Matrix.toString(Jmat.Matrix.mul(result.q, result.r));

  var result = Jmat.Matrix.qr(Jmat.Matrix(3,3,1,2,3,4,5,6,7,8,9));
  Jmat.Matrix.toString(result.q) + ' \n ' + Jmat.Matrix.toString(result.r) + ' \n ' + Jmat.Matrix.toString(Jmat.Matrix.mul(result.q, result.r));

  var result = Jmat.Matrix.qr(Jmat.Matrix(3,3,12,-51,4,6,167,-68,-4,24,-41));
  Jmat.Matrix.toString(result.q) + ' \n ' + Jmat.Matrix.toString(result.r) + ' \n ' + Jmat.Matrix.toString(Jmat.Matrix.mul(result.q, result.r));
  */

  //if(m.h < m.w) return null; //seems to work anyway, so don't do this check. TODO: verify this

  var t = Math.min(m.h - 1, m.w);
  var real = Jmat.Matrix.isReal(m);
  var a = Jmat.Matrix.copy(m);
  var r;
  var q;

  for(var k = 0; k < t; k++) {
    var x = Jmat.Matrix.col(a, 0);

    var xk = a.e[0][0];
    var normx = Jmat.Matrix.norm(x);

    var alpha;
    if(xk.im != 0) alpha = Jmat.Complex.exp(Jmat.Complex.newi(Jmat.Complex.argr(xk))).neg().mul(normx);
    else if(xk.re < 0) alpha = normx;
    else alpha = normx.neg();

    var u = Jmat.Matrix.col(a, 0);
    u.e[0][0] = u.e[0][0].sub(alpha);
    var normu = Jmat.Matrix.norm(u);

    var v = Jmat.Matrix.divc(u, normu);

    var vv;
    if(real) {
      vv = Jmat.Matrix.mulr(Jmat.Matrix.mul(v, Jmat.Matrix.transpose(v)), 2);
    } else {
      var xhv = Jmat.Matrix.mul(Jmat.Matrix.transjugate(x), v);
      var vhx = Jmat.Matrix.mul(Jmat.Matrix.transjugate(v), x);
      var w = xhv.e[0][0].div(vhx.e[0][0]);
      vv = Jmat.Matrix.mulc(Jmat.Matrix.mul(v, Jmat.Matrix.transjugate(v)), w.inc());
    }
    var id = Jmat.Matrix.identity(a.h, a.h);
    var qk = Jmat.Matrix.sub(id, vv); // here, qk*x = [alpha, 0,...,0]^T

    if (k + 1 < t) {
      a = Jmat.Matrix.mul(qk, a);
      a = Jmat.Matrix.minorsub(a, 0, 0);
    }

    if(k == 0) {
      q = Jmat.Matrix.transjugate(qk);
    } else {
      qk = Jmat.Matrix.overlap(Jmat.Matrix.identity(k, k), qk, k, k);
      q = Jmat.Matrix.mul(q, Jmat.Matrix.transjugate(qk)); //Q1^h * Q2^h * ... * Qt^h
    }
    //r = Jmat.Matrix.mul(qk, r); // Qt * ... * Q2 * Q1 * A //not needed to calculate here, done below instead, is q^t * m
  }
  r = Jmat.Matrix.mul(Jmat.Matrix.transjugate(q), m);

  // Solve numerical problem: 0-values of the upper triangular matrix sometimes become a tiny e-16 value
  for(var y = 0; y < r.h; y++) {
    for(var x = 0; x < r.w && x < y; x++) {
      // every value below diagonal should be 0, but leave big ones so that miscalculation bugs can be seen.
      if(Jmat.Complex.absr(r.e[y][x]) < 1e-15) r.e[y][x] = Jmat.Complex(0);
    }
  }

  return { q: q, r: r };
};

// eigenvalues and vectors of 1x1 matrix
Jmat.Matrix.eig11 = function(m) {
  if(m.w != 1 || m.h != 1) return null;
  var result = {};
  result.l = new Jmat.Matrix(1, 1);
  result.l.e[0][0] = Jmat.Complex(0);
  result.v = new Jmat.Matrix(1, 1);
  result.v.e[0][0] = Jmat.Complex(1);
  return result;
};

// explicit algebraic formula for eigenvalues and vectors of 2x2 matrix
Jmat.Matrix.eig22 = function(m) {
  if(m.w != 2 || m.h != 2) return null;
  var a = Jmat.Complex(1);
  var b = m.e[0][0].neg().sub(m.e[1][1]);
  var c = m.e[0][0].mul(m.e[1][1]).sub(m.e[0][1].mul(m.e[1][0]));
  var d = Jmat.Complex.sqrt(b.mul(b).sub(a.mul(c).mulr(4)));
  var l1 = b.neg().add(d).div(a.mulr(2));
  var l2 = b.neg().sub(d).div(a.mulr(2));

  var v12 = l1.sub(m.e[0][0]).div(m.e[0][1]);
  var v11 = Jmat.Complex(1);
  /*//normalize
  var n = Jmat.Complex.sqrt(v12.mul(v12).addr(1));
  v12 = v12.div(n);
  v11 = v11.div(n);*/
  var v22 = l2.sub(m.e[0][0]).div(m.e[0][1]);
  var v21 = Jmat.Complex(1);
  /*//normalize
  n = Jmat.Complex.sqrt(v12.mul(v12).addr(1));
  v12 = v12.div(n);
  v11 = v11.div(n);*/

  var result = {};
  result.l = new Jmat.Matrix(2, 1);
  result.l.e[0][0] = l1;
  result.l.e[1][0] = l2;
  result.v = new Jmat.Matrix(2, 2);
  result.v.e[0][0] = v11;
  result.v.e[1][0] = v12;
  result.v.e[0][1] = v21;
  result.v.e[1][1] = v22;
  return result;
};

// explicit algebraic formula for eigenvalues of 2x2 matrix
Jmat.Matrix.eigval22 = function(m) {
  if(m.w != 2 || m.h != 2) return null;
  var a = Jmat.Complex(1);
  var b = m.e[0][0].add(m.e[1][1]);
  var c = m.e[0][0].mul(m.e[1][1]).sub(m.e[0][1].mul(m.e[1][0]));
  var d = Jmat.Complex.sqrt(b.mul(b).sub(a.mul(c).mulr(4)));
  var l1 = b.add(d).div(a.add(a));
  var l2 = b.sub(d).div(a.add(a));
  return [l1, l2];
};

// Returns the eigenvectors and eigenvalues of m as { l: eigenvalues, v: eigenvectors }
// eigenvalues as n*1 column vector, eigenvectors as n*n matrix
// for each column of v and corresponding eigenvalue: A*v = l*v (l represents lambda, A is m)
Jmat.Matrix.eig = function(m) {
  /*
  Checks in console:
  var result = Jmat.Matrix.eig(Jmat.Matrix(2,2,1,2,3,4))
  Jmat.Matrix.toString(result.l) + ' \n ' + Jmat.Matrix.toString(result.v);
  */
  if(m.w != m.h || m.w < 1) return null;
  var n = m.w;
  if(n == 1) return Jmat.Matrix.eig11(m);
  if(n == 2) return Jmat.Matrix.eig22(m);

  // using the QR algorithm
  var a = Jmat.Matrix.copy(m); //will contain the eigenvalues on the diagonal

  // Naive implicit QR without shifts, does not work in many cases, left commented out for demo purpose only
  // for(var i = 0; i < 30; i++) {
  //  var qr = Jmat.Matrix.qr(a);
  //  a = Jmat.Matrix.mul(qr.r, qr.q); // RQ instead of QR: A_k -> QR, A_(k+1) = RQ
  // }

  // QR with double shifting. This because with single shift or no shift, it does not support complex eigenvalues of real matrix, e.g. [[1,-1],[5,-1]]
  // TODO: this is slow, optimize like with Hessenberg form
  var id = Jmat.Matrix.identity(n, n);
  for(var i = 0; i < 15; i++) {
    // var s = a.e[a.h - 1][a.w - 1]; //value that would be chosen for single shift
    // double shift: choose two sigma's, the eigenvalues of the bottom right 2x2 matrix (for which we have the explicit solution)
    var l = Jmat.Matrix.eigval22(Jmat.Matrix.submatrix(a, a.h - 2, a.h, a.w - 2, a.w));

    var si0 = id.mulc(l[0]);
    a = a.sub(si0);
    var qr = Jmat.Matrix.qr(a);
    a = Jmat.Matrix.mul(qr.r, qr.q).add(si0);

    var si1 = id.mulc(l[1]);
    a = a.sub(si1);
    qr = Jmat.Matrix.qr(a);
    a = Jmat.Matrix.mul(qr.r, qr.q).add(si1);
  }

  var v = new Jmat.Matrix(n, n);

  // Find eigenvectors by solving system of linear equations.
  // TODO: this is not very efficient...
  // Normally, the product of all the qr.q's of the loop above should give the eigenvectors, but that applies only for symmetric matrices while this is supposed to support all
  // So, instead, solve system equation (A - lambda * I) * x = 0, but with last element of 0 set to 1, and bottom row of (A - lambda * I) set to 0,0,...,0,1.
  // That makes the system solvable, and makes each vector have its last element be 1.
  for(var j = 0; j < n; j++) {
    var value = a.e[j][j];
    var e = Jmat.Matrix.copy(m); //TODO: this makes it even slower, copy only the needed columns
    for(var i = 0; i < n; i++) e.e[i][i] = e.e[i][i].sub(value);
    for(var i = 0; i < n; i++) e.e[e.h - 1][i] = Jmat.Complex(i == n - 1 ? 1 : 0);
    var f = Jmat.Matrix.zero(n, 1);
    f.e[f.h - 1][0] = Jmat.Complex(1);
    var g = Jmat.Matrix.solve(e, f);
    for(var i = 0; i < n; i++) v.e[i][j] = g.e[i][0];
  }

  var l = new Jmat.Matrix(m.w, 1);
  for(var i = 0; i < m.w; i++) l.e[i][0] = a.e[i][i];

  return { l: l, v: v };
};

// Returns the eigen decomposition of m as { v: V, d: D }
// If M is diagonizable, M = V * D * V^(-1)
// In other words: m == result.v.mul(result.d).mul(Jmat.Matrix.inv(result.v))
// This function is very similar to Jmat.Matrix.eig. v is the same, d is the same as l but put on the diagonal of a matrix
Jmat.Matrix.eigd = function(m) {
  var eig = Jmat.Matrix.eig(m);

  return { v: eig.v, d: Jmat.Matrix.diag(eig.l) };
};

// Puts all the elements of d in a single diagonal matrix
Jmat.Matrix.diag = function(d) {
  var n = d.w * d.h;
  var result = Jmat.Matrix.zero(n, n);
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
Jmat.Matrix.get1 = function(m, i) {
  var x = i % m.w;
  var y = Math.floor(i / m.w);
  return m.e[y][x];
};
Jmat.Matrix.prototype.get1 = function(i) {
  return Jmat.Matrix.get1(this, i);
};

// Set element using a one-dimensional index. See Jmat.Matrix.get1.
Jmat.Matrix.set1 = function(m, i, v) {
  var x = i % m.w;
  var y = Math.floor(i / m.w);
  m.e[y][x] = v;
};

// Cross product between two vectors of length 3, that is, 3x1 and/or 1x3 matrices. Other input dimensions are not accepted.
// Return value has dimensions of the first input (that is, if first input is row vector, it's row vector, otherwise column vector)
// Jmat.Matrix.cross(Jmat.Matrix([1,2,3]), Jmat.Matrix([4,5,6])).toString()
Jmat.Matrix.cross = function(a, b) {
  if(a.w * a.h != 3 || b.w * b.h != 3) return Jmat.Matrix(NaN);
  var c = new Jmat.Matrix(a.h, a.w);
  Jmat.Matrix.set1(c, 0, a.get1(1).mul(b.get1(2)).sub(a.get1(2).mul(b.get1(1))));
  Jmat.Matrix.set1(c, 1, a.get1(2).mul(b.get1(0)).sub(a.get1(0).mul(b.get1(2))));
  Jmat.Matrix.set1(c, 2, a.get1(0).mul(b.get1(1)).sub(a.get1(1).mul(b.get1(0))));
  return c;
};

// Dot product of two vectors.
// Also supports it for matrices of same dimensions (it then is the Frobenius inner product)
// If vectors, row and column vectors may be mixed.
Jmat.Matrix.dot = function(a, b) {
  if(a.w != b.w || a.h != b.h) {
    if(!(a.w == b.h && a.h == b.w && (a.w == 1 || a.h == 1))) return Jmat.Matrix(NaN); // Do allow it for differently orientated vectors (h or w is 1)
  }
  var n = a.w * a.h;
  var result = Jmat.Complex(0);
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
Jmat.Matrix.zsvdc_ = function(x, ldx, n, p, s, e, u, ldu, v, ldv, work, job) {
  var i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,
      mm,mm1,mp1,nct,ncu,nrt,info; // integers
  var maxit = 30;
  var t,r; //complex
  var b,c,cs,el,emm1,f,g,dznrm2,scale,shift,sl,sm,sn,
      smm1,t1,test,ztest; // double

  var dr; //drotg result

  var dreal = function(z) { return z.re; };
  var cabs1 = function(z) { return Math.abs(z.re) + Math.abs(z.im); };
  // returns value with absolute value of x, argument of y (transfers sign)
  var csign = function(x, y) { return y.eqr(0) ? Jmat.Complex(0) : y.mulr(x.absr() / y.absr()); };
  var sign = function(x, y) { return y == 0 ? 0 : (y < 0 ? -Math.abs(x) : Math.abs(x)); };
  // Euclidean norm of complex vector, n elements starting at index start
  var dznrm2 = function(n, arr, start) {
    var result = Jmat.Complex(0);
    for(var i = 0; i < n; i++) {
      var e = arr[start + i];
      result = result.add(e.mul(e));
    }
    return Jmat.Complex.sqrt(result);
  };
  // sets vector arry to alpha * arrx + arry
  var zaxpy = function(n, alpha, arrx, startx, arry, starty) {
    for(var i = 0; i < n; i++) {
      arry[starty + i] = arry[starty + i].add(arrx[startx + i].mul(alpha));
    }
  };
  // dot product
  var zdotc = function(n, arrx, startx, arry, starty) {
    var result = Jmat.Complex(0);
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
  ncu = jobu > 1 ? Math.min(n, p) : n;
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
      s[l] = Jmat.Complex(dznrm2(n - l, x, l + l * ldx));
      if(cabs1(s[l]) != 0.0) {
        if(cabs1(x[l + l * ldx]) != 0.0) s[l] = csign(s[l], x[l + l * ldx]);
        t = Jmat.Complex(1.0).div(s[l]);
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
      e[l] = Jmat.Complex(dznrm2(p - l - 1, e, lp1));

      if(cabs1(e[l]) != 0.0) {
        if(cabs1(e[lp1]) != 0.0) {
          e[l] = csign(e[l], e[lp1]);
        }
        t = Jmat.Complex(1.0) / e[l];
        zscal(p - l - 1, t, e, lp1);
        e[lp1] = Jmat.Complex(1.0) + e[lp1];
      }
      e[l] = e[l].conj().neg();
      // apply the transformation.
      if(lp1 < n && cabs1(e[l]) != 0.0) {
        for(j = lp1; j < n; j++) {
          work[j] = Jmat.Complex(0.0);
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
  if(n < m) s[m - 1] = Jmat.Complex(0.0);
  if(nrt + 1 < m) e[nrt] = x[nrt + (m - 1) * ldx];
  e[m - 1] = Jmat.Complex(0.0);
  // if required, generate u.
  if(wantu) {
    for(j = nct; j < ncu; j++) {
      for(i = 0; i < n; i++) {
        u[i + j * ldu] = Jmat.Complex(0.0);
      }
      u[j + j * ldu] = Jmat.Complex(1.0);
    }
    for(ll = 0; ll < nct; ll++) {
      l = nct - ll - 1;
      if(cabs1(s[l]) == 0.0) {
        for(i = 0; i < n; i++) {
          u[i + l * ldu] = Jmat.Complex(0.0);
        }
        u[l + l * ldu] = Jmat.Complex(1.0);
      } else {
        lp1 = l + 1;
        for(j = lp1; j < ncu; j++) {
          t = zdotc(n - l, u, l + l * ldu, u, l + j * ldu).neg().div(u[l + l * ldu]);
          zaxpy(n - l, t, u, l + l * ldu, u, l + j * ldu);
        }
        zscal(n - l, Jmat.Complex(-1.0), u, l + l * ldu);
        u[l + l * ldu] = u[l + l * ldu].inc();
        for(i = 0; i < l; i++) {
          u[i + l * ldu] = Jmat.Complex(0.0);
        }
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
        v[i + l * ldv] = Jmat.Complex(0.0);
      }
      v[l + l * ldv] = Jmat.Complex(1.0);
    }
  }
  // transform s and e so that they are real.
  for(i = 0; i < m; i++) {
    if(cabs1(s[i]) != 0.0) {
      t = Jmat.Complex(Jmat.Complex.absr(s[i]));
      r = s[i].div(t);
      s[i] = t;
      if(i + 1 < m) e[i] = e[i].div(r);
      if(wantu) zscal(n, r, u, i * ldu);
    }
    if(i + 1 == m) break;
    if(cabs1(e[i]) != 0.0) {
      t = Jmat.Complex(Jmat.Complex.absr(e[i]));
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
      test = Jmat.Complex.absr(s[l - 1]) + Jmat.Complex.absr(s[l]);
      ztest = test + Jmat.Complex.absr(e[l - 1]);
      if(ztest == test) {
        e[l - 1] = Jmat.Complex(0.0);
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
        if(ls != m) test = test + Jmat.Complex.absr(e[ls - 1]);
        if(ls != l + 1) test = test + Jmat.Complex.absr(e[ls - 2]);
        ztest = test + Jmat.Complex.absr(s[ls - 1]);
        if(ztest == test) {
          s[ls - 1] = Jmat.Complex(0.0);
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
      e[m - 2] = Jmat.Complex(0.0);
      for(kk = 1; kk <= mm1; kk++) {
        k = mm1 - kk + l;
        t1 = dreal(s[k - 1]);
        dr = drotg(t1, f); t1 = dr[0]; f = dr[1]; cs = dr[2]; sn = dr[3];
        s[k - 1] = Jmat.Complex(t1);
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
      e[l - 2] = Jmat.Complex(0.0);
      for(k = l; k <= m; k++) {
        t1 = dreal(s[k - 1]);
        dr = drotg(t1, f); t1 = dr[0]; f = dr[1]; cs = dr[2]; sn = dr[3];
        s[k - 1] = Jmat.Complex(t1);
        f = -sn * dreal(e[k - 1]);
        e[k - 1] = e[k - 1].mulr(cs);
        if(wantu) zdrot(n, u, (k - 1) * ldu, u, (l - 2) * ldu, cs, sn);
      }
    }
    // perform one qr step.
    else if(kase == 3) {
      // calculate the shift.
      scale = Math.max(Math.max(Math.max(Math.max(Jmat.Complex.absr(s[m - 1]),
              Jmat.Complex.absr(s[m - 2])), Jmat.Complex.absr(e[m - 2])), Jmat.Complex.absr(s[l - 1])),
              Jmat.Complex.absr(e[l - 1]));
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
        if(k != l) e[k - 2] = Jmat.Complex(f);
        f = cs * dreal(s[k - 1]) + sn * dreal(e[k - 1]);
        e[k - 1] = e[k - 1].mulr(cs).sub(s[k - 1].mulr(sn));
        g = sn * dreal(s[k]);
        s[k] = s[k].mulr(cs);
        if(wantv) zdrot(p, v, (k - 1) * ldv, v, k * ldv, cs, sn);
        dr = drotg(f, g); f = dr[0]; g = dr[1]; cs = dr[2]; sn = dr[3];
        s[k - 1] = Jmat.Complex(f);
        f = cs * dreal(e[k - 1]) + sn * dreal(s[k]);
        s[k] = e[k - 1].mulr(-sn).add(s[k].mulr(cs));
        g = sn * dreal(e[k]);
        e[k] = e[k].mulr(cs);
        if(wantu && k < n) zdrot(n, u, (k - 1) * ldu, u, k * ldu, cs, sn);
      }
      e[m - 2] = Jmat.Complex(f);
      iter++;
    } else if(kase == 4) {
      // convergence.
      // make the singular value positive.
      if(dreal(s[l - 1]) < 0.0) {
        s[l - 1] = s[l - 1].neg();
        if(wantv) zscal(p, Jmat.Complex(-1.0), v, (l - 1) * ldv);
      }
      // order the singular values.
      while(l != mm) {
        if(dreal(s[l]) <= dreal(s[l - 1])) break;
        t = s[l - 1];
        s[l - 1] = s[l];
        s[l] = t;
        if(wantv && l < p) zswap(p, v, (l - 1) * ldv, v, l * ldv);
        if(wantu && l < n) zswap(n, u, (l - 1) * ldu, u, l * ldu);
        l++;
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
Jmat.Matrix.svd = function(m) {
  /*
  Checks in console:
  function testSvd(m) {
    var result = Jmat.Matrix.svd(Jmat.Matrix(m));
    console.log(Jmat.toString(result) + ' | ' + Jmat.Matrix.mul(Jmat.Matrix.mul(result.u, result.s), Jmat.Matrix.transpose(result.v)).toString());
  }
  testSvd(Jmat.Matrix(2,2,1,2,3,4))
  testSvd(Jmat.Matrix(2,2,1,2,3,4).mulc(Jmat.Complex.I))
  testSvd(Jmat.Matrix(2,2,1,2,1,2))
  testSvd(Jmat.Matrix([[1,2]]))
  testSvd(Jmat.Matrix([[1],[2]]))
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
  Jmat.Matrix.zsvdc_(a, m.h, m.h, m.w, s, e, u, m.h, v, m.w, [], 11);
  
  // Solve numerical problems: singular values < eta should be 0
  var eta = 1e-15; //TODO: check if this tolerance can be improved 
  for(var i = 0; i < s.length; i++) if(Math.abs(s[i]) < eta) s[i] = 0;

  var result = { u: new Jmat.Matrix(m.h, m.h), s: new Jmat.Matrix(m.h, m.w), v: new Jmat.Matrix(m.w, m.w) };


  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.h; x++) {
      result.u.e[y][x] = u[y + x * m.h];
    }
  }
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      result.s.e[y][x] = (x == y) ? Jmat.Complex(s[x]) : Jmat.Complex(0);
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
Jmat.Matrix.eq = function(a, b) {
  if(a.w != b.w || a.h != b.h) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      return a.e[y][x].eq(b.e[y][x]);
    }
  }

  return true;
};

// nearly equal
Jmat.Matrix.near = function(a, b, epsilon) {
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
Jmat.Matrix.solve = function(a, b) {
  if(a.h != b.h) return null;
  var ag = Jmat.Matrix.pseudoinverse(a);

  var aag = Jmat.Matrix.mul(a, ag);
  var aagb = Jmat.Matrix.mul(aag, b);
  if(!Jmat.Matrix.near(b, aagb)) return null; //inconsistent system with no solution

  return Jmat.Matrix.mul(ag, b);
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
Jmat.Matrix.kiss_fft_state_ = function(nfft, inverse) {
  this.nfft = nfft;
  this.inverse = inverse;
  this.factors = []; //int array, size 32
  this.twiddles = []; //complex Jmat.Complex array, size nfft-1
};
Jmat.Matrix.kf_bfly2_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
  var j = 0;
  for(var i = 0; i < m; i++) {
    var t = Fout[Fout_index + i + m].mul(st.twiddles[j]);
    j += fstride;
    Fout[Fout_index + i + m] = Fout[Fout_index + i].sub(t);
    Fout[Fout_index + i] = Fout[Fout_index + i].add(t);
  }
};
Jmat.Matrix.kf_bfly4_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
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
Jmat.Matrix.kf_bfly3_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
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
Jmat.Matrix.kf_bfly5_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
  var scratch = []; //size-13 complex array
  var ya = st.twiddles[fstride*m];
  var yb = st.twiddles[fstride*2*m];
  var m2 = 2 * m;
  var m3 = 3 * m;
  var m4 = 4 * m;
  for (var u=0; u<m; ++u ) {
    scratch[0] = Jmat.Complex(Fout[Fout_index + u]);
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
    scratch[5] = Jmat.Complex(0);
    scratch[5].re = scratch[0].re + scratch[7].re*ya.re + scratch[8].re*yb.re;
    scratch[5].im = scratch[0].im + scratch[7].im*ya.re + scratch[8].im*yb.re;
    scratch[6] = Jmat.Complex(0);
    scratch[6].re = scratch[10].im*ya.im + scratch[9].im*yb.im;
    scratch[6].im = -scratch[10].re*ya.im - scratch[9].re*yb.im;
    Fout[Fout_index + m+u]=scratch[5].sub(scratch[6]);
    Fout[Fout_index + m4+u]=scratch[5].add(scratch[6]);
    scratch[11] = Jmat.Complex(0);
    scratch[11].re = scratch[0].re + scratch[7].re*yb.re + scratch[8].re*ya.re;
    scratch[11].im = scratch[0].im + scratch[7].im*yb.re + scratch[8].im*ya.re;
    scratch[12] = Jmat.Complex(0);
    scratch[12].re = -scratch[10].im*yb.im + scratch[9].im*ya.im;
    scratch[12].im = scratch[10].re*yb.im - scratch[9].re*ya.im;
    Fout[Fout_index + m2+u]=scratch[11].add(scratch[12]);
    Fout[Fout_index + m3+u]=scratch[11].sub(scratch[12]);
  }
};
// perform the butterfly for one stage of a mixed radix FFT
Jmat.Matrix.kf_bfly_generic_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/, p /*int*/) {
  var u,k,q1,q; /*int*/
  var t; // complex Jmat.Complex
  var Norig = st.nfft;
  var scratch = [];
  for ( u=0; u<m; ++u ) {
    k=u;
    for ( q1=0 ; q1<p ; ++q1 ) {
      scratch[q1] = Jmat.Complex(Fout[Fout_index + k]);
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
Jmat.Matrix.kf_work_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, f /*array of complex Jmat.Complex*/,f_index,
    fstride /*int*/, in_stride /*int*/, factors /*int array*/, factors_index,st /*kiss_fft_state*/) {
  var p = factors[factors_index + 0]; /* the radix */
  var m = factors[factors_index + 1]; /* stage's fft length/p */
  var j = 0;

  if (m==1) {
    for(var i = 0; i < p*m; i++) {
      Fout[i + Fout_index] = Jmat.Complex(f[f_index + j]);
      j += fstride*in_stride;
    }
  }else{
    for(var i = 0; i < p*m; i += m) {
      // recursive call:
      // DFT of size m*p performed by doing p instances of smaller DFTs of size m, each one takes a decimated version of the input
      Jmat.Matrix.kf_work_(Fout, Fout_index + i, f, f_index + j, fstride*p, in_stride, factors, factors_index + 2, st);
      j += fstride*in_stride;
    }
  }
  // recombine the p smaller DFTs
  switch (p) {
    case 2: Jmat.Matrix.kf_bfly2_(Fout,Fout_index,fstride,st,m); break;
    case 3: Jmat.Matrix.kf_bfly3_(Fout,Fout_index,fstride,st,m); break;
    case 4: Jmat.Matrix.kf_bfly4_(Fout,Fout_index,fstride,st,m); break;
    case 5: Jmat.Matrix.kf_bfly5_(Fout,Fout_index,fstride,st,m); break;
    default: Jmat.Matrix.kf_bfly_generic_(Fout,Fout_index,fstride,st,m,p); break;
  }
};
// facbuf is populated by p1,m1,p2,m2, ... where p[i] * m[i] = m[i-1], m0 = n 
Jmat.Matrix.kf_factor_ = function(n, facbuf) {
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
Jmat.Matrix.kiss_fft_alloc_ = function(nfft,inverse_fft) {
  var st = new Jmat.Matrix.kiss_fft_state_(nfft, inverse_fft);
  for (var i=0;i<nfft;++i) {
    var phase = -Math.PI*2*i / nfft;
    if (st.inverse) phase *= -1;
    st.twiddles[i] = new Jmat.Complex(Math.cos(phase), Math.sin(phase));
  }
  Jmat.Matrix.kf_factor_(nfft,st.factors);
  return st;
};
Jmat.Matrix.kiss_fft_ = function(st/*kiss_fft_state object*/,fin/*complex Jmat.Complex array of size nfft*/,fout/*complex Jmat.Complex array of size nfft*/) {
  Jmat.Matrix.kf_work_(fout,0, fin,0, 1, 1/*in_stride*/, st.factors,0,st);
};
// End of Kiss FFT
////////////////////////////////////////////////////////////////////////////////


Jmat.Matrix.matrixfft_ = function(m, inverse) {
  var rowresult = new Jmat.Matrix(m.h, m.w);

  // apply to each row
  if(m.w > 1) {
    for(var j = 0; j < m.h; j++) {
      var out = [];
      for(var i = 0; i < m.w; i++) out[i] = Jmat.Complex(0);
      var st = Jmat.Matrix.kiss_fft_alloc_(m.w, inverse);
      Jmat.Matrix.kiss_fft_(st, m.e[j], out);
      for(var i = 0; i < m.w; i++) rowresult.e[j][i] = out[i];
    }
  } else {
    rowresult = m;
  }

  var result = new Jmat.Matrix(m.h, m.w);

  // apply to each column
  if (m.h > 1) {
    for(var j = 0; j < m.w; j++) {
      var col = Jmat.Matrix.transpose(Jmat.Matrix.col(rowresult, j));
      var out = [];
      for(var i = 0; i < m.h; i++) out[i] = Jmat.Complex(0);
      var st = Jmat.Matrix.kiss_fft_alloc_(m.h, inverse);
      Jmat.Matrix.kiss_fft_(st, col.e[0], out);
      for(var i = 0; i < m.h; i++) result.e[i][j] = out[i];
    }
  } else {
    result = rowresult;
  }

  var factor = 1.0 / Math.sqrt(m.w * m.h);
  result = Jmat.Matrix.mulr(result, factor);

  return result;
};

// Discrete fourier transform
Jmat.Matrix.fft = function(m) {
  // Debug in console:
  // Jmat.Matrix.toString(Jmat.Matrix.fft(Jmat.Matrix(2,2,1,2,3,4)))
  // should give [5 -1][-2 0]
  return Jmat.Matrix.matrixfft_(m, 0);
};

// Inverse discrete fourier transform
Jmat.Matrix.ifft = function(m) {
  return Jmat.Matrix.matrixfft_(m, 1);
};


//Matrix exponential
Jmat.Matrix.exp = function(m) {
  if(m.h != m.w) return null; //must be square

  var result = m.add(Jmat.Matrix.identity(m.w, m.w));
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
var a = Jmat.Matrix(2,2,1,2,3,4);
var c = Jmat.Matrix.cos(a); var s = Jmat.Matrix.sin(a);
Jmat.Matrix.toString(c) + ' ' + Jmat.Matrix.toString(s) + ' ' + Jmat.Matrix.toString(c.mul(c).add(s.mul(s)));
--> c*c+s*s should give identity matrix
*/
Jmat.Matrix.cos = function(m) {
  if(m.h != m.w) return null; //must be square

  var result = Jmat.Matrix.identity(m.w, m.w);
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
Jmat.Matrix.sin = function(m) {
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
// debug in console: var a = Jmat.Matrix.sqrt(Jmat.Matrix(3,3,1,1,0,0,0,1,1,0,1)); Jmat.Matrix.toString(a) + ' ' + Jmat.Matrix.toString(a.mul(a))
Jmat.Matrix.sqrt = function(m) {
  if(m.h != m.w) return null; //must be square

  // Babylonian method. Does not work for [[1,2][3,4]], because that one has complex result. Left commented out for demonstration purpose only.
  /*var result = m.add(Jmat.Matrix.identity(m.w, m.w)).mulr(0.5);
  for(var i = 0; i <= 30; i++) {
    result = result.add(m.div(result)).mulr(0.5);
  }
  return result;*/

  // With eigen decomposition: only the sqrt of the diagonals of D needs to be taken
  // TODO: this will only work for diagonizable matrix. Implement a way for non-diagonizable one as well (with Jordan normal form)
  var e = Jmat.Matrix.eigd(m);
  var v = e.v;
  var d = e.d;
  for(var i = 0; i < d.w; i++) d.e[i][i] = Jmat.Complex.sqrt(d.e[i][i]);
  return v.mul(d).mul(Jmat.Matrix.inv(v));
};

// Matrix logarithm (base e, ln)
// debug in console: var a = Jmat.Matrix.log(Jmat.Matrix(2,2,1,2,3,4)); Jmat.Matrix.toString(a) + ' ' + Jmat.Matrix.toString(Jmat.Matrix.exp(a))
Jmat.Matrix.log = function(m) {
  if(m.h != m.w) return null; //must be square

  // With eigen decomposition: only the log of the diagonals of D needs to be taken
  // TODO: this will only work for diagonizable matrix. Implement a way for non-diagonizable one as well (with Jordan normal form)
  var e = Jmat.Matrix.eigd(m);
  var v = e.v;
  var d = e.d;
  for(var i = 0; i < d.w; i++) d.e[i][i] = Jmat.Complex.log(d.e[i][i]);
  return v.mul(d).mul(Jmat.Matrix.inv(v));
};

// Matrix to any complex scalar power s
Jmat.Matrix.powc = function(m, s) {
  if(m.h != m.w) return null; //must be square

  // With eigen decomposition: only the log of the diagonals of D needs to be taken
  // TODO: this will only work for diagonizable matrix. Implement a way for non-diagonizable one as well (with Jordan normal form)
  var e = Jmat.Matrix.eigd(m);
  var v = e.v;
  var d = e.d;
  for(var i = 0; i < d.w; i++) d.e[i][i] = d.e[i][i].pow(s);
  return v.mul(d).mul(Jmat.Matrix.inv(v));
};

////////////////////////////////////////////////////////////////////////////////

// Expose everything with easier names (disable this if it would cause name clashes with other libraries)

var Real = Jmat.Real;
var Complex = Jmat.Complex;
var Matrix = Jmat.Matrix;

// Expose some of the Real functions into JS Math
['gamma', 'factorial', 'lambertw', 'erf', 'erfc'].map(function(fun) { if(!Math[fun]) Math[fun] = Jmat.Real[fun]; });
