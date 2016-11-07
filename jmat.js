/** @license
Jmat.js

Copyright (c) 2011-2016, Lode Vandevenne
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
*/

// REQUIRES: jmat_real.js, jmat_complex.js, jmat_matrix.js, jmat_quaternion.js, jmat_special.js, jmat_bigint.js

/*
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

*) Example of a function with many arguments: hypergeometric2F1
Jmat.hypergeometric2F1(0.5, 2, 3, 4).toString()
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
  var i = 0; var z = Complex(0); for(;;) { if(z.abs() > 2) break; z = z.mul(z).add(c); i++; if(i > 60) return Complex(0); } return Complex.polar(1, i * Math.PI / 60);
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

NOTE: The input types {number|Complex} and {Array|Matrix} almost always also accept a string
formatted in such way that it can be parsed to a number or a matrix, e.g. '5+2i' or '[[1,2],[3,4+5i]]'
*/

// Elementary operators

// TODO: add all BigInt functions here, many are currently missing.

/* Add. x,y:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.add = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.add(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.add(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y));
  if(Jmat.bigIntIn_(x) || Jmat.bigIntIn_(y)) return Jmat.BigInt.add(Jmat.BigInt.cast(x), Jmat.BigInt.cast(y));
  return Jmat.Complex.add(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Subtract. x,y:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.sub = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.sub(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.sub(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y));
  if(Jmat.bigIntIn_(x) || Jmat.bigIntIn_(y)) return Jmat.BigInt.sub(Jmat.BigInt.cast(x), Jmat.BigInt.cast(y));
  return Jmat.Complex.sub(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Multiply. x,y:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.mul = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.mul(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  if(Jmat.matrixIn_(x)) return Jmat.Matrix.mulc(Jmat.Matrix.cast(x), Jmat.Complex.cast(y));
  if(Jmat.matrixIn_(y)) return Jmat.Matrix.mulc(Jmat.Matrix.cast(y), Jmat.Complex.cast(x));
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.mul(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y));
  if(Jmat.bigIntIn_(x) || Jmat.bigIntIn_(y)) return Jmat.BigInt.mul(Jmat.BigInt.cast(x), Jmat.BigInt.cast(y));
  return Jmat.Complex.mul(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Division. x,y:{number|Jmat.Complex|Jmat.Matrix}. returns {Jmat.Complex|Jmat.Matrix}. */
Jmat.div = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.div(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  if(Jmat.matrixIn_(x)) return Jmat.Matrix.divc(Jmat.Matrix.cast(x), Jmat.Complex.cast(y));
  if(Jmat.matrixIn_(y)) return Jmat.Matrix.divc(Jmat.Matrix.cast(y), Jmat.Complex.cast(x));
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.div(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y));
  if(Jmat.bigIntIn_(x) || Jmat.bigIntIn_(y)) return Jmat.BigInt.div(Jmat.BigInt.cast(x), Jmat.BigInt.cast(y));
  return Jmat.Complex.div(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};

// Compare

/* Equal? x,y:{number|Complex|Matrix}. returns {boolean}. */
Jmat.eq = function(x, y) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.eq(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y));
  if(Jmat.matrixIn_(x) || Jmat.matrixIn_(y)) return false; //one is matrix, the other is not
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.eq(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y));
  if(Jmat.bigIntIn_(x) || Jmat.bigIntIn_(y)) return Jmat.BigInt.eq(Jmat.BigInt.cast(x), Jmat.BigInt.cast(y));
  return Jmat.Complex.eq(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Nearly equal? x,y:{number|Complex|Matrix}. epsilon:{number}. returns {boolean}. */
Jmat.near = function(x, y, epsilon) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.near(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y), Jmat.Real.caststrict(epsilon));
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.near(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y), Jmat.Real.caststrict(epsilon));
  if(Jmat.bigIntIn_(x) || Jmat.bigIntIn_(y)) return Jmat.BigInt.eq(Jmat.BigInt.cast(x), Jmat.BigInt.cast(y));
  return Jmat.Complex.near(Jmat.Complex.cast(x), Jmat.Complex.cast(y), Jmat.Real.caststrict(epsilon));
};
/* Nearly equal? With relative precision. x,y:{number|Complex|Matrix}. precision:{number}. returns {boolean}.
   Precision must be near 0 but slightly larger, e.g. 0.001 for 3 digits of precision, 1e-5 for 5 digits, ...*/
Jmat.relnear = function(x, y, precision) {
  if(Jmat.matrixIn_(x) && Jmat.matrixIn_(y)) return Jmat.Matrix.relnear(Jmat.Matrix.cast(x), Jmat.Matrix.cast(y), Jmat.Real.caststrict(precision));
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.relnear(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y), Jmat.Real.caststrict(precision));
  if(Jmat.bigIntIn_(x) || Jmat.bigIntIn_(y)) return Jmat.BigInt.eq(Jmat.BigInt.cast(x), Jmat.BigInt.cast(y));
  return Jmat.Complex.relnear(Jmat.Complex.cast(x), Jmat.Complex.cast(y), Jmat.Real.caststrict(precision));
};

// Categories

Jmat.isNaN = function(x) { return Jmat.matrixIn_(x) ? Jmat.Matrix.isNaN(Jmat.Matrix.cast(x)) : Jmat.Complex.isNaN(Jmat.Complex.cast(x)); };

// Power & Logarithms

/* Power, or matrix power. x:{number|Complex|Matrix}. y:{number|Complex}. returns {Complex}. */
Jmat.pow = function(x, y) {
  if(Jmat.matrixIn_(x)) return Jmat.Matrix.powc(Jmat.Matrix.cast(x), Jmat.Complex.cast(y)); // matrix raised to any complex number power
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.pow(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y));
  return Jmat.Complex.pow(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Square root, or matrix square root. z:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.sqrt = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.sqrt(Jmat.Matrix.cast(z));
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.sqrt(Jmat.Quaternion.cast(z));
  if(Jmat.bigIntIn_(z)) return Jmat.BigInt.sqrt(Jmat.BigInt.cast(z));
  return Jmat.Complex.sqrt(Jmat.Complex.cast(z));
};
/* Exponential, or matrix exponential. z:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.exp = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.exp(Jmat.Matrix.cast(z)); // matrix exponential
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.exp(Jmat.Quaternion.cast(z));
  return Jmat.Complex.exp(Jmat.Complex.cast(z));
};
/* Logarithm, or matrix logarithm. z:{number|Complex|Matrix}. returns {Complex|Matrix}. */
Jmat.log = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.log(Jmat.Matrix.cast(z));
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.log(Jmat.Quaternion.cast(z));
  return Jmat.Complex.log(Jmat.Complex.cast(z));
};
/* log(z) + 1. z:{number|Complex}. returns {Complex}. */
Jmat.log1p = function(z) { return Jmat.Complex.log1p(Jmat.Complex.cast(z)); };
/* exp(z) - 1. z:{number|Complex}. returns {Complex}. */
Jmat.expm1 = function(z) { return Jmat.Complex.expm1(Jmat.Complex.cast(z)); };
/* Base-y logarithm of x. x,y:{number|Complex}. returns {Complex}. */
Jmat.logy = function(x, y) {
  if(Jmat.quaternionIn_(x) || Jmat.quaternionIn_(y)) return Jmat.Quaternion.logy(Jmat.Quaternion.cast(x), Jmat.Quaternion.cast(y));
  return Jmat.Complex.logy(Jmat.Complex.cast(x), Jmat.Complex.cast(y));
};
/* Base-2 logarithm of z:{number|Complex}. returns {Complex}. */
Jmat.log2 = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.log2(Jmat.Quaternion.cast(z));
  return Jmat.Complex.log2(Jmat.Complex.cast(z));
};
/* Base-10 logarithm of z. z:{number|Complex}. returns {Complex}. */
Jmat.log10 = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.log10(Jmat.Quaternion.cast(z));
  return Jmat.Complex.log10(Jmat.Complex.cast(z));
};
/* Principal branch of LambertW. z:{number|Complex}. returns {Complex}. */
Jmat.lambertw = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.lambertw(Jmat.Quaternion.cast(z));
  return Jmat.Complex.lambertw(Jmat.Complex.cast(z));
};
/* Negative branch of LambertW. z:{number|Complex}. returns {Complex}. */
Jmat.lambertwm = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.lambertwm(Jmat.Quaternion.cast(z));
  return Jmat.Complex.lambertwm(Jmat.Complex.cast(z));
};
/* Specific branch of LambertW. branch:{number} must be integer, z:{number|Complex}. returns {Complex}. */
Jmat.lambertwb = function(branch, z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.lambertwb(branch, Jmat.Quaternion.cast(z));
  return Jmat.Complex.lambertwb(Jmat.Real.caststrict(branch), Jmat.Complex.cast(z));
};
/* Tetration (power tower). x,y:{number|Complex}. returns {Complex} */
Jmat.tetration = function(x, y) { return Jmat.Complex.tetration(Jmat.Complex.cast(x), Jmat.Complex.cast(y)); };
Jmat.tetrational = function(z) { return Jmat.Complex.tetrational(Jmat.Complex.cast(z)); };

// Elementary functions

/* Identity function. x:{number|Complex}. returns {Complex} */
Jmat.x = function(x) { return Jmat.Complex.cast(x); };
/* Negate. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.neg = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.neg(Jmat.Matrix.cast(z));
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.neg(Jmat.Quaternion.cast(z));
  return Jmat.Complex.neg(Jmat.Complex.cast(z));
};
/* Reciproke. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.inv = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.inv(Jmat.Matrix.cast(z));
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.inv(Jmat.Quaternion.cast(z));
  return Jmat.Complex.inv(Jmat.Complex.cast(z));
};
/* Complex conjugate. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.conj = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.conj(Jmat.Matrix.cast(z));
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.conj(Jmat.Quaternion.cast(z));
  return Jmat.Complex.conj(Jmat.Complex.cast(z));
};

// Trigonometric functions

/* Sine, or matrix-sine. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.sin = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.sin(Jmat.Matrix.cast(z));
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.zin(Jmat.Quaternion.cast(z));
  return Jmat.Complex.sin(Jmat.Complex.cast(z));
};
/* Cosine, or matrix-cosine. z:{number|Complex|Matrix}. returns {Complex|Matrix} */
Jmat.cos = function(z) {
  if(Jmat.matrixIn_(z)) return Jmat.Matrix.cos(Jmat.Matrix.cast(z));
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.cos(Jmat.Quaternion.cast(z));
  return Jmat.Complex.cos(Jmat.Complex.cast(z));
};
/* Tangent. z:{number|Complex}. returns {Complex} */
Jmat.tan = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.tan(Jmat.Quaternion.cast(z));
  return Jmat.Complex.tan(Jmat.Complex.cast(z));
};
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
Jmat.factorial = function(z) { return Jmat.bigIntIn_(z) ? Jmat.BigInt.factorial(Jmat.BigInt.cast(z)) : Jmat.Complex.factorial(Jmat.Complex.cast(z)); };
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
Jmat.abs = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.abs(Jmat.Quaternion.cast(z));
  return Jmat.Complex.abs(Jmat.Complex.cast(z));
};
/* Complex argument or phase. z:{number|Complex}. returns {Complex} */
Jmat.arg = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.qrg(Jmat.Quaternion.cast(z));
  return Jmat.Complex.arg(Jmat.Complex.cast(z));
};
/* Real sign. z:{number|Complex}. returns {Complex} */
Jmat.sign = function(z) {
  if(Jmat.quaternionIn_(z)) return Jmat.Quaternion.sign(Jmat.Quaternion.cast(z));
  return Jmat.Complex.sign(Jmat.Complex.cast(z));
};
/* Euclidean distance. a,b:{number|Complex|Quaternion|Matrix|BigInt}. returns {number} */
Jmat.dist = function(a, b) {
  if(Jmat.quaternionIn_(a) || Jmat.quaternionIn_(b)) return Jmat.Quaternion.dist(Jmat.Quaternion.cast(a), Jmat.Quaternion.cast(b));
  else if(Jmat.bigIntIn_(a) || Jmat.bigIntIn_(b)) return Jmat.BigInt.dist(Jmat.BigInt.cast(a), Jmat.BigInt.cast(b));
  else if(Jmat.matrixIn_(a) || Jmat.matrixIn_(b)) return Jmat.Matrix.dist(Jmat.Matrix.cast(a), Jmat.Matrix.cast(b));
  return Jmat.Complex.dist(Jmat.Complex.cast(a), Jmat.Complex.cast(b));
};
/* Chebyshev distance. a,b:{number|Complex|Quaternion|Matrix|BigInt}. returns {number} */
Jmat.cheb = function(a, b) {
  if(Jmat.quaternionIn_(a) || Jmat.quaternionIn_(b)) return Jmat.Quaternion.cheb(Jmat.Quaternion.cast(a), Jmat.Quaternion.cast(b));
  else if(Jmat.bigIntIn_(a) || Jmat.bigIntIn_(b)) return Jmat.BigInt.cheb(Jmat.BigInt.cast(a), Jmat.BigInt.cast(b));
  else if(Jmat.matrixIn_(a) || Jmat.matrixIn_(b)) return Jmat.Matrix.cheb(Jmat.Matrix.cast(a), Jmat.Matrix.cast(b));
  return Jmat.Complex.cheb(Jmat.Complex.cast(a), Jmat.Complex.cast(b));
};
/* Manhattan distance. a,b:{number|Complex|Quaternion|Matrix|BigInt}. returns {number} */
Jmat.manhattan = function(a, b) {
  if(Jmat.quaternionIn_(a) || Jmat.quaternionIn_(b)) return Jmat.Quaternion.manhattan(Jmat.Quaternion.cast(a), Jmat.Quaternion.cast(b));
  else if(Jmat.bigIntIn_(a) || Jmat.bigIntIn_(b)) return Jmat.BigInt.manhattan(Jmat.BigInt.cast(a), Jmat.BigInt.cast(b));
  else if(Jmat.matrixIn_(a) || Jmat.matrixIn_(b)) return Jmat.Matrix.manhattan(Jmat.Matrix.cast(a), Jmat.Matrix.cast(b));
  return Jmat.Complex.manhattan(Jmat.Complex.cast(a), Jmat.Complex.cast(b));
};
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

// Interpolation

/* Linear interpolation.  a,b,x:{number|Complex}. returns {Complex} */
Jmat.lerp = function(a, b, x) { return Jmat.Complex.lerp(Jmat.Complex.cast(a), Jmat.Complex.cast(b), Jmat.Complex.cast(x)); };


// Cylindrical functions

/* Bessel function of the first kind.  nu,z:{number|Complex}. returns {Complex} */
Jmat.besselj = function(nu, z) { return Jmat.Complex.besselj(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
/* Bessel function of the second kind.  nu,z:{number|Complex}. returns {Complex} */
Jmat.bessely = function(nu, z) { return Jmat.Complex.bessely(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
/* Modified bessel function of the first kind.  nu,z:{number|Complex}. returns {Complex} */
Jmat.besseli = function(nu, z) { return Jmat.Complex.besseli(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
/* Modified bessel function of the second kind.  nu,z:{number|Complex}. returns {Complex} */
Jmat.besselk = function(nu, z) { return Jmat.Complex.besselk(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
/* Hankel function of the first kind.  nu,z:{number|Complex}. returns {Complex} */
Jmat.hankelh1 = function(nu, z) { return Jmat.Complex.hankelh1(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
/* Hankel function of the second kind.  nu,z:{number|Complex}. returns {Complex} */
Jmat.hankelh2 = function(nu, z) { return Jmat.Complex.hankelh2(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
/* Airy function Ai. z:{number|Complex}. returns {Complex} */
Jmat.airy = function(z) { return Jmat.Complex.airy(Jmat.Complex.cast(z)); };
/* Bairy function Bi. z:{number|Complex}. returns {Complex} */
Jmat.bairy = function(z) { return Jmat.Complex.bairy(Jmat.Complex.cast(z)); };
/* Derivative of Airy function Ai. z:{number|Complex}. returns {Complex} */
Jmat.airy_deriv = function(z) { return Jmat.Complex.airy_deriv(Jmat.Complex.cast(z)); };
/* Derivative of Bairy function Bi. z:{number|Complex}. returns {Complex} */
Jmat.bairy_deriv = function(z) { return Jmat.Complex.bairy_deriv(Jmat.Complex.cast(z)); };
/* Struve function H_nu. nu,z:{number|Complex}. returns {Complex} */
Jmat.struveh = function(nu, z) { return Jmat.Complex.struveh(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
Jmat.struvek = function(nu, z) { return Jmat.Complex.struvek(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
Jmat.struvel = function(nu, z) { return Jmat.Complex.struvel(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
Jmat.struvem = function(nu, z) { return Jmat.Complex.struvem(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
Jmat.angerj = function(nu, z) { return Jmat.Complex.angerj(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };
Jmat.webere = function(nu, z) { return Jmat.Complex.webere(Jmat.Complex.cast(nu), Jmat.Complex.cast(z)); };

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
/* Hypergeometric function 2F1. Also works outside of convergence zone. a,b,c,z:{number|Complex}. returns {Complex} */
Jmat.hypergeometric2F1 = function(a, b, c, z) { return Jmat.Complex.hypergeometric2F1(Jmat.Complex.cast(a), Jmat.Complex.cast(b), Jmat.Complex.cast(c), Jmat.Complex.cast(z)); };
/* Generalized hypergeometric function. a,b are arrays, e.g. a.length == 3, b.length == 2 for 3F2. Only valid if the basic series converges. a,b:{Array.<number|Complex>}. c:{number|Complex}. returns {Complex} */
Jmat.hypergeometric = function(a, b, z) {
  var a2 = []; for(var i = 0; i < a.length; i++) a2[i] = Jmat.Complex.cast(a[i]);
  var b2 = []; for(var i = 0; i < b.length; i++) b2[i] = Jmat.Complex.cast(b[i]);
  return Jmat.Complex.hypergeometric(a2, b2, Jmat.Complex.cast(z));
};

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

/* Bitwise not, negating. x:{number|Complex}. returns {Complex} */
Jmat.bitneg = function(x) { return Jmat.Complex.bitneg(Jmat.Complex.cast(x)); };
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

/* Modulo division. a,b:{number|Complex}. Result has sign of b. returns {Complex}. */
Jmat.mod = function(a, b) { return Jmat.Complex.mod(Jmat.Complex.cast(a), Jmat.Complex.cast(b)); }; // result has sign of divisor (unlike JS '%' operator)
/* Remainder. a,b:{number|Complex}. Result has sign of a. returns {Complex} */
Jmat.rem = function(a, b) { return Jmat.Complex.rem(Jmat.Complex.cast(a), Jmat.Complex.cast(b)); }; // result has sign of dividend (same result as JS '%' operator on real numbers)
/* Wrap x between from and to (to excluded). x,to,from:{number|Complex}. returns {Complex} */
Jmat.wrap = function(x, from, to) { return Jmat.Complex.wrap(Jmat.Complex.cast(x), Jmat.Complex.cast(from), Jmat.Complex.cast(to)); };
/* Clamp x between from and to (to included). x,to,from:{number|Complex}. returns {Complex} */
Jmat.clamp = function(x, from, to) { return Jmat.Complex.clamp(Jmat.Complex.cast(x), Jmat.Complex.cast(from), Jmat.Complex.cast(to)); };

// Other special functions

/* Minkowski's question mark function. x:{number|Complex}. returns {Complex} */
Jmat.minkowski = function(x) { return Jmat.Complex.minkowski(Jmat.Complex.cast(x)); };

// Matrix (NOTE: more are above: add, sub, mul, div, inv, neg, conj, exp, log, sqrt, cos, sin)

/* LUP decomposition. m:{Array|Matrix}. returns {Object.<string, Matrix>} object with p:P, l:L, u:U */
Jmat.lu = function(m) { return Jmat.Matrix.lu(Jmat.Matrix.cast(m)); };
/* QR decomposition. m:{Array|Matrix}. returns {Object.<string, Matrix>} object with q:Q, r:R */
Jmat.qr = function(m) { return Jmat.Matrix.qr(Jmat.Matrix.cast(m)); };
/* singular value decomposition. m:{Array|Matrix}. returns {Object.<string, Matrix>} object with u:U, s:S, v:V */
Jmat.svd = function(m) { return Jmat.Matrix.svd(Jmat.Matrix.cast(m)); };
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
/* Eigenvalue decomposition (spectral decomposition). m:{Array|Matrix}. returns {Object.<string, Matrix>} object with v:eigenvectors, d:eigenvalues on diagonal */
Jmat.evd = function(m) { return Jmat.Matrix.evd(Jmat.Matrix.cast(m)); };
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
/* Newton's method. f,df:{function(Complex):Complex} df is derivative of f, z0:{number|Complex}, maxiter:{number} integer. returns {Complex} */
Jmat.rootfind_newton = function(f, df, z0, maxiter) { return Jmat.Complex.rootfind_newton(f, df, Jmat.Complex.cast(z0), Jmat.Real.caststrict(maxiter)); };
/* Newton's method, but without giving derivative (it is calculated like the secant method). f:{function(Complex):Complex}, z0:{number|Complex}, maxiter:{number} integer. returns {Complex} */
Jmat.rootfind_newton_noderiv = function(f, z0, maxiter) { return Jmat.Complex.rootfind_newton_noderiv(f, Jmat.Complex.cast(z0), Jmat.Real.caststrict(maxiter)); };

Jmat.polyroots = function(coeffs) { return Jmat.Complex.polyroots(Jmat.Complex.castArray(coeffs)); }

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

// Test if input is a quaternion
Jmat.quaternionIn_ = function(v) {
  if(!v) return false;
  if(typeof v == 'string' && (v.indexOf('j') >= 0 || v.indexOf('k') >= 0)) {
    return true; //not for 'i', then it is just complex. This function should only return true if it's unambiguously a quaternion
  }
  return v && (v instanceof Jmat.Quaternion);
};

// Test if input is a BigInt
Jmat.bigIntIn_ = function(v) {
  if(!v) return false;
  if(typeof v == 'string' && v.length > 16) {
    var ok = true;
    for(var i = 1; i < v.length && i < 17; i++) { // test for only digits (except first may be '-'). biggest JS integer (2^53) has around 16 digits, no need to check for more to exclude matrix/complex/...
      var c = v.charCodeAt(i);
      if(c < 48 || v > 57) { ok = false ; break; };
    }
    if(ok) return true;
  }
  return v && (v instanceof Jmat.BigInt);
  // No test for array, that can already mean matrix
};

Jmat.toString_ = function(a, opt_render, opt_precision) {
  if(!a) return '' + a;
  var result = '';
  if(opt_render && a.render) {
    return ('\n' + a.render(opt_precision));
  }
  // For arrays of known types
  if(typeof a == 'object' && Array.isArray(a)) {
    result += '[';
    for(var i = 0; i < a.length; i++) result += (Jmat.toString_(a[i], opt_render, opt_precision) + (i + 1 == a.length ? '' : ', '));
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
      result += Jmat.toString_(a[it], opt_render, opt_precision);
      comma = true;
    }
    result += '}';
    return result;
  }
  // Prefer toString implementation if available (that is where this function is better than JSON.stringify! :))
  result = (opt_render && a.render) ? ('\n' + a.render(opt_precision)) : (a.toString ? a.toString(opt_precision) : ('' + a));
  return result;
};

// Nice string of any known object, like Complex and Matrix and the objects of named matrices returned by svd or array returned by factorize
Jmat.toString = function(a) {
  return Jmat.toString_(a, false);
};

Jmat.print = function(a) {
  console.log(Jmat.toString(a));
}

// Nice string of any known object, with possibly some multiline rendered matrices
Jmat.render = function(a, opt_precision) {
  return Jmat.toString_(a, true, opt_precision);
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Export API
////////////////////////////////////////////////////////////////////////////////

// Expose everything with easier names (disable this if it would cause name clashes with other libraries)

var Real = Jmat.Real; // R
var Complex = Jmat.Complex; // C
var Matrix = Jmat.Matrix; // M
var BigInt = Jmat.BigIntC || Jmat.BigInt; // B
var Quaternion = Jmat.Quaternion; // Q

// Expose some of the Real functions into JS Math
['gamma', 'factorial', 'lambertw', 'erf', 'erfc'].map(function(fun) { if(!Math[fun]) Math[fun] = Jmat.Real[fun]; });
// And a few more, but these are likely to become part of JS Math already in a future ECMAScript revision
['sinh', 'cosh', 'tanh', 'asinh', 'acosh', 'atanh', 'clz32', 'log2', 'log10', 'log1p', 'expm1'].map(function(fun) { if(!Math[fun]) Math[fun] = Jmat.Real[fun]; });
