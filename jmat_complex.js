/*
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

// REQUIRES: jmat_real.js

/*
Jmat.Complex: arithmetic on complex numbers

Overview of some functionality:
-elementary arithmetic: Complex.add, Complex.sub, Complex.mul, Complex.div
-mathematical functions: Complex.pow, Complex.exp, Complex.sqrt, Complex.log, Complex.cos, Complex.cosh, Complex.acos, ...
-special functions: Complex.erf, Complex.lambertw, Complex.gamma (more are in jmat_special.js)
-fft
*/

/*
Constructor, but also usable without new as factory function.

Class representing a complex value with real and imaginary part.

Serves as both the actual object used as complex number, and namespace for most
complex mathematical functions.

The only sad thing is that Javascript doesn't support operator overloading
and nice expressions like a + b have to become a.add(b) instead.

When used as factory function, makes it easy to construct values, e.g. if you
set C = Jmat.Complex, you can create real value 1 with C(1), or a complex value
with C('1+2i') or C(1, 2).
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
  if(a == undefined) return new Jmat.Complex(0, 0);
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

// Create a new Jmat.Complex value, polar, angle in degrees
Jmat.Complex.polarDeg = function(r, a) {
  a = Jmat.Real.degToRad(a);
  return new Jmat.Complex(r * Math.cos(a), r * Math.sin(a));
};

// Casts the given number type to Jmat.Complex. If the given type is already of type Jmat.Complex, does not copy it but returns the input.
Jmat.Complex.cast = function(v) {
  if(v && v.re != undefined) return v;
  if(v == undefined) return Jmat.Complex(0);
  return Jmat.Complex(v);
};

Jmat.Complex.castArray = function(a) {
  var result = [];
  for (var i = 0; i < a.length; i++) result[i] = Jmat.Complex.cast(a[i]);
  return result;
};

//aka clone
Jmat.Complex.copy = function(v) {
  return new Jmat.Complex(v.re, v.im);
};

// Because JS number toFixed appends zeros
Jmat.Complex.formatFloat_ = function(value, precision) {
  var power = Math.pow(10, precision || 0);
  return String(Math.round(value * power) / power);
};

//debugstring
Jmat.Complex.toString = function(value, opt_precision) {
  if(!value) return value == 0 ? 'invalid0' : ('' + value);
  var vre = value.re;
  var vim = value.im;
  if(!opt_precision && Math.abs(vre) < 1e-15 && Math.abs(vim) > 0.1) vre = 0; // avoid something hard to read like "2.220446049250313e-16+1i" due to numerical imprecision
  if(!opt_precision && Math.abs(vim) < 1e-15 && Math.abs(vre) > 0.1) vim = 0; // avoid something hard to read like "2.220446049250313e-16+1i" due to numerical imprecision

  var re = (opt_precision ? Jmat.Complex.formatFloat_(vre, opt_precision) : ('' + vre));
  var im = (opt_precision ? Jmat.Complex.formatFloat_(vim, opt_precision) : ('' + vim));
  if(im == '1') im = '';
  if(im == '-1') im = '-';

  if(vim == 0 || im == '0') return '' + re;
  if(vre == 0) return '' + im + 'i';
  if(vim < 0) return '' + re + im + 'i';
  return '' + re + '+' + im + 'i';
};
Jmat.Complex.prototype.toString = function(opt_precision) {
  return Jmat.Complex.toString(this, opt_precision);
};

Jmat.Complex.toStringPolar = function(value, opt_precision) {
  if(!value) return value == 0 ? 'invalid0' : ('' + value);
  var r = (opt_precision ? Jmat.Complex.formatFloat_(value.abs(), opt_precision) : ('' + value.abs()));
  var a = (opt_precision ? Jmat.Complex.formatFloat_(value.arg(), opt_precision) : ('' + value.arg()));

  return '' + r + '∠' + a;
};
Jmat.Complex.prototype.toStringPolar = function(opt_precision) {
  return Jmat.Complex.toStringPolar(this, opt_precision);
};

// Test out this rendering: var c = Complex.polarDeg(1, 45); var d = Complex(1); for(var i = 0; i < 10; i++) { console.log(d.toStringPolarDeg() + ' ' + d.toString()); d = d.mul(c); }
Jmat.Complex.toStringPolarDeg = function(value, opt_precision) {
  if(!value) return value == 0 ? 'invalid0' : ('' + value);
  var arg = Jmat.Real.radToDeg(value.arg());
  if(!opt_precision && Math.abs(arg - Math.round(arg)) < 1e-10) arg = Math.round(arg); // avoid something like "1∠89.99999999999999°" due to numerical imprecision


  var r = (opt_precision ? Jmat.Complex.formatFloat_(value.abs(), opt_precision) : ('' + value.abs()));
  var a = (opt_precision ? Jmat.Complex.formatFloat_(arg, opt_precision) : ('' + arg));

  return '' + r + '∠' + a + '°';
};
Jmat.Complex.prototype.toStringPolarDeg = function(opt_precision) {
  return Jmat.Complex.toStringPolarDeg(this, opt_precision);
};

// Parses strings of the form '5', '5+i', '5-2.3i', '1.25e-25+17.37e5i'
// Cannot (yet) parse the polar forms.
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
    d = 1.0 / d; // optimization: avoid multiple times the same division
    var re, im;
    if(d > 1) {
      re = (x.re * y.re + x.im * y.im) * d;
      im = (x.im * y.re - x.re * y.im) * d;
    } else {
      // the multiplications with d are in the center, to avoid overflow in case e.g. x.re*y.re overflows floating point
      re = x.re * d * y.re + x.im * d * y.im;
      im = x.im * d * y.re - x.re * d * y.im;
    }
    return new Jmat.Complex(re, im);
  }
};
Jmat.Complex.prototype.div = function(y) {
  return Jmat.Complex.div(this, y);
};

Jmat.Complex.addr = function(z, a) {
  return new Jmat.Complex(z.re + a, z.im);
};
Jmat.Complex.prototype.addr = function(a) {
  return new Jmat.Complex(this.re + a, this.im);
};
Jmat.Complex.addi = function(z, a) {
  return new Jmat.Complex(z.re, z.im + a);
};
Jmat.Complex.prototype.addi = function(a) {
  return new Jmat.Complex(this.re, this.im + a);
};

Jmat.Complex.subr = function(z, a) {
  return new Jmat.Complex(z.re - a, z.im);
};
Jmat.Complex.prototype.subr = function(a) {
  return new Jmat.Complex(this.re - a, this.im);
};
// Subtract z from real.
Jmat.Complex.rsub = function(a, z) {
  return new Jmat.Complex(a - z.re, -z.im);
};
// Subtract self from real. This operator exists because it's less awkward to write z.rsub(3) than Complex(3).sub(z) in long formulas
Jmat.Complex.prototype.rsub = function(a) {
  return new Jmat.Complex(a - this.re, -this.im);
};
Jmat.Complex.subi = function(z, a) {
  return new Jmat.Complex(z.re, z.im - a);
};
Jmat.Complex.prototype.subi = function(a) {
  return new Jmat.Complex(this.re, this.im - a);
};

Jmat.Complex.mulr = function(z, a) {
  return new Jmat.Complex(z.re * a, z.im * a);
};
Jmat.Complex.prototype.mulr = function(a) {
  return new Jmat.Complex(this.re * a, this.im * a);
};
// multiply with imaginary number given as real
Jmat.Complex.muli = function(z, a) {
  return new Jmat.Complex(-z.im * a, z.re * a);
};
Jmat.Complex.prototype.muli = function(a) {
  return new Jmat.Complex(-this.im * a, this.re * a);
};

Jmat.Complex.divr = function(z, a) {
  return new Jmat.Complex(z.re / a, z.im / a);
};
Jmat.Complex.prototype.divr = function(a) {
  return new Jmat.Complex(this.re / a, this.im / a);
};
// Divide real a through z.
Jmat.Complex.rdiv = function(a, z) {
  return new Jmat.Complex.div(Jmat.Complex(a), z);
};
// Divide real a through self. This operator exists because it's less awkward to write z.rsub(3) than Complex(3).div(z) in long formulas
Jmat.Complex.prototype.rdiv = function(a) {
  return Jmat.Complex.div(Jmat.Complex(a), this);
};
// divide through imaginary number given as real
Jmat.Complex.divi = function(z, a) {
  return new Jmat.Complex(z.im / a, -z.re / a);
};
Jmat.Complex.prototype.divi = function(a) {
  return new Jmat.Complex(this.im / a, -this.re / a);
};

//rotate complex number z by a radians. That is, change its argument. a is real (JS number).
Jmat.Complex.rotate = function(z, a) {
  if(a == 0) return z;
  return Jmat.Complex.polar(z.abs(), z.arg() + a);
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

//Like JS ~: returns -(x + 1), limited to 32-bit int
Jmat.Complex.bitneg = function(x) {
  var result = Jmat.Complex(0);
  result.re = ~x.re;
  //imaginary part not bit-negated on purpose: otherwise it appears when bit-inverting real number, which is in 99.9% of the cases not wanted
  //instead negated, to follow the formula -(x + 1)
  //result.im = ~x.im;
  result.im = -x.im
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
// Another name for this could be "normalize" as it makes the length of the "vector" 1
Jmat.Complex.sign = function(z) {
  if (z.im == 0) {
    if(z.re == 0) return Jmat.Complex(0);
    else if(z.re < 0) return Jmat.Complex(-1);
    return Jmat.Complex(1);
  }

  return z.divr(z.abs());
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

// Similar to sign, but returns 1 if the input is 0
Jmat.Complex.sign1 = function(z) {
  if (z.im == 0) {
    if(z.re < 0) return Jmat.Complex(-1);
    return Jmat.Complex(1);
  }

  return z.divr(z.abs());
};

// Similar to csgn, but returns 1 if the input is 0
Jmat.Complex.csgn1 = function(z) {
  if (Jmat.Real.near(z.re, 0, 1e-15)) { //avoid numeric imprecisions
    if(z.im < 0) return Jmat.Complex(-1);
    return Jmat.Complex(1);
  } else {
    if(z.re < 0) return Jmat.Complex(-1);
    Jmat.Complex(1);
  }
};

// applies sign of y to x, so the result has the magnitude of x but the argument of y
Jmat.Complex.copysign = function(x, y) {
  return Jmat.Complex.abs(x).mul(Jmat.Complex.sign(y));
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
Jmat.Complex.prototype.dec = function() {
  return new Jmat.Complex(this.re - 1, this.im);
};

// TODO: consider no longer have prototype.abs return real and Complex.abs return Complex. Use absr for real instead.
// absolute value, aka modulus of complex number, as a Jmat.Complex object (its imaginary part is 0)
Jmat.Complex.abs = function(x) {
  return Jmat.Complex(x.abs());
};
// absolute value, aka modulus of complex number, returned as real (regular JS number, to be similar to .re and .im)
Jmat.Complex.prototype.abs = function() {
  if(this.im == 0) return Math.abs(this.re);
  if(this.re == 0) return Math.abs(this.im);

  if(this.re == Infinity || this.re == -Infinity || this.im == Infinity || this.im == -Infinity) {
    return Infinity;
  }

  // Numerically more stable version of "Math.sqrt(x.re * x.re + x.im * x.im);"
  return Jmat.Real.hypot(this.re, this.im);
};

// absolute value squared, returned as Jmat.Complex object. This is faster than abs due to not taking sqrt.
Jmat.Complex.abssq = function(x) {
  return Jmat.Complex(x.re * x.re + x.im * x.im);
};
// absolute value squared, returned as real (regular JS number). This is faster than abs due to not taking sqrt.
Jmat.Complex.prototype.abssq = function() {
  return this.re * this.re + this.im * this.im;
};

// returns the complex argument in range -PI to +PI, as a Jmat.Complex object (its imaginary part is 0)
Jmat.Complex.arg = function(x) {
  return Jmat.Complex(x.arg());
};
// returns the complex argument in range -PI to +PI, as a real (regular JS number, to be similar to .re and .im)
Jmat.Complex.prototype.arg = function() {
  if(this.im == 0) return this.re < 0 ? Math.PI : 0;
  return Math.atan2(this.im, this.re);
};

//returns result in range 0-1 rather than -PI to PI, as a regular JS number. Useful for graphical representations, not for math. 0 matches 0 degrees, 0.5 matches 180 degrees, 0.999 matches around 359 degrees.
Jmat.Complex.arg1 = function(z) {
  var result = z.arg();
  if(result < 0) result += 2 * Math.PI;
  result /= (2 * Math.PI);
  if(result < 0) result = 0;
  if(result > 1) result = 1;
  return result;
};

//manhattan norm, returned as real
Jmat.Complex.abs1r = function(x) {
  return Math.abs(x.re) + Math.abs(x.im);
};

// returns sqrt(|x|^2 + |y|^2) with numerically more precise formula ; a companion to atan2
// if more than 2 arguments are given, calculates norm of all arguments
Jmat.Complex.hypot = function(x, y) {
  var C = Jmat.Complex;
  if(C.abs1r(y) > C.abs1r(x)) {
    var temp = x;
    x = y;
    y = temp;
  }
  if(C.isInf(x)) return Infinity;
  var t = y.div(x);
  return x.mul(C.sqrt(C.abssq(t).addr(1)));
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
  return Math.abs(z.re) == Infinity || Math.abs(z.im) == Infinity;
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
  var C = Jmat.Complex;
  if(C.isReal(x) && C.isReal(y) && (x.re >= 0 || y.re == Infinity || y.re == -Infinity || Jmat.Real.isInt(y.re))) {
    //if(x.re == 0 && y.re == 0) return C(NaN); // JS's pow returns 1 for 0^0
    // It is chosen to return 1 for 0^0, not NaN. NaN is mathematically more correct, however 0^0 is correct in many practical applications.
    return C(Math.pow(x.re, y.re));
  } else if(x.eqr(0)) {
    return y.re == 0 ? C(NaN) : C(y.re < 0 ? Infinity : 0);
  } else {
    // This is just one branch. In fact it returns a complex result for -3 ^ (1/3),
    // the cube root of -3. To get the real result, use absolute value (and then negate) on it.
    // This is correct: the principal result of the cube root for this is a complex number.
    // Note: This returns incorrect values for a negative real to the power of Infinity: the result should be -Infinity for < -1, 0 for > -1, NaN for -1, but it always gives NaN. However, the "if" part above already handles that.
    var r = x.abs();
    var t = x.arg();
    var u = Math.pow(r, y.re) * Math.exp(-y.im * t);
    if(isNaN(u)) {
      u = Math.pow(1, y.re / r) * Math.exp(-y.im * t / r);
      if(u < 0) u = -Infinity;
      else if(u > 0) u = Infinity;
      else u = NaN;
    }
    var v = y.im * Math.log(r) + y.re * t;
    return C(u * Math.cos(v), u * Math.sin(v));
  }
};
Jmat.Complex.prototype.pow = function(y) {
  return Jmat.Complex.pow(this, y);
};

Jmat.Complex.powr = function(z, a) {
  return Jmat.Complex.pow(z, Jmat.Complex(a));
};
Jmat.Complex.prototype.powr = function(a) {
  return Jmat.Complex.pow(this, Jmat.Complex(a));
};

// raise regular js number x, to complex power a
Jmat.Complex.rpow = function(x, a) {
  return Jmat.Complex.pow(Jmat.Complex(x), a);
};
// raise regular js number x, to this complex number
Jmat.Complex.prototype.rpow = function(x) {
  return Jmat.Complex.pow(Jmat.Complex(x), this);
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
  var C = Jmat.Complex;
  if(!C.isReal(x) || !C.isReal(y)) {
    if(y.eqr(0)) return C(Math.PI / 2);
    // As per the definition, return atan of (x / y)
    return C.atan(x.div(y));
  } else {
    var result = C(0);
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
  // at z.re > 709, exp gives infinity, at z.re > 304, it gets wrong value if z.im != 0
  if(z.re > 350) return Jmat.Complex(1);
  if(z.re < -350) return Jmat.Complex(-1);
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
  return Jmat.Complex.log(z.addr(1).div(z.rsub(1))).mulr(0.5);
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
  result.re = Jmat.Real.round(x.re);
  result.im = Jmat.Real.round(x.im);
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

// Linear interpolation
Jmat.Complex.lerp = function(a, b, x) {
  return x.rsub(1).mul(a).add(x.mul(b));
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
  if(x.abssq() < 1e-5) return x.add(x.mul(x).divr(2)).add(x.mul(x).mul(x).divr(6));
  else return Jmat.Complex.exp(x).subr(1);
};

//natural log (base e, ln)
Jmat.Complex.log = function(x) {
  if(x.eqr(-Infinity)) return Jmat.Complex(Infinity);

  if(Jmat.Complex.isReal(x) && x.re >= 0) {
    return Jmat.Complex(Math.log(x.re));
  }

  return Jmat.Complex(Math.log(x.abs()), x.arg());
};

//ln(x + 1), with better precision for x around 0
Jmat.Complex.log1p = function(x) {
  if(x.abssq() < 1e-8) return x.mulr(-0.5).addr(1).mul(x);
  else return Jmat.Complex.log(x.addr(1));
};

//arbitrary log: log_y(x), y is also complex number
//warning: base y is second argument
Jmat.Complex.logy = function(x, y) {
  return Jmat.Complex.log(x).div(Jmat.Complex.log(y));
};

//arbitrary log: log_y(x), where y is a regular JS number
//warning: base y is second argument
Jmat.Complex.logr = function(x, y) {
  return Jmat.Complex.log(x).divr(Math.log(y));
};

Jmat.Complex.log2 = function(x) {
  return Jmat.Complex.log(x).divr(Math.LN2);
};

Jmat.Complex.log10 = function(x) {
  return Jmat.Complex.log(x).divr(Math.LN10);
};

// Complex square root. Sqrt has two solutions, this function always returns the solution with nonnegative real part.
Jmat.Complex.sqrt = function(x) {
  if(Jmat.Complex.isReal(x)) {
    var result = Jmat.Complex(0);
    if(x.re >= 0 || x.re != x.re) result.re = Math.sqrt(x.re);
    else result.im = Math.sqrt(-x.re);
    return result;
  } else return x.pow(Jmat.Complex(0.5));
};

Jmat.Complex.root = function(x, y) {
  return x.pow(Jmat.Complex(Jmat.Complex.inv(y)));
};

Jmat.Complex.rootr = function(x, y) {
  return x.pow(Jmat.Complex(1 / y));
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

  return value.divr(value.abs());
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


Jmat.Complex.near = function(x, y, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  //return Jmat.Complex.manhattan(x, y) <= epsilon;
  // Manhattan NOT used, because then this function returns false for equal infinities
  return x.re - epsilon <= y.re && x.re + epsilon >= y.re && x.im - epsilon <= y.im && x.im + epsilon >= y.im;
};

// Near regular JS number y
Jmat.Complex.nearr = function(x, y, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  // Manhattan NOT used, because then this function returns false for equal infinities
  return x.re - epsilon <= y && x.re + epsilon >= y && x.im - epsilon <= 0 && x.im + epsilon >= 0;
};

// dist, cheb and manhattan all return regular real JS numbers for all types. In some types they are all the same, but not for e.g. Complex or Matrix.
// Euclidean distance
Jmat.Complex.dist = function(a, b) {
  return a.sub(b).abs();
};
//Chebyshev distance
Jmat.Complex.cheb = function(a, b) {
  return Math.max(Jmat.Real.dist(a.re, b.re), Jmat.Real.dist(a.im, b.im));
};
// Manhattan distance
Jmat.Complex.manhattan = function(a, b) {
  return Math.max(Math.abs(a.re - b.re), Math.abs(a.im - b.im));
};

/*
Precision must be near 0 but slightly larger, e.g. 0.001 for 3 digits of precision, 1e-5 for 5 digits, ...
That many digits must match, starting from the first non-zero digit.
That means, if one value is zero and the other is not, no matter how close to zero the other is, this function will always return false.
For complex numbers, allows different argument, as long as the distance between the two numbers is relatively tiny compared to their magnitude.
*/
Jmat.Complex.relnear = function(x, y, precision) {
  //return Jmat.Real.relnear(x.re, y.re, precision) && Jmat.Real.relnear(x.im, y.im, precision);
  if(x.eq(y)) return true;
  return x.sub(y).abs() < (Math.max(x.abs(), y.abs()) * precision);
};

// This works for quaternions as well.. set M to Jmat.Complex for complex number, or Jmat.Quaternion for quaternions.
Jmat.Complex.lambertwb_generic_ = function(M, branch, z) {
  if(M.isReal(z) && z.re > -0.36 /*~ -1/e*/ && branch == 0) return M(Jmat.Real.lambertw(z));

  if(!Jmat.Real.isInt(branch)) return M(NaN);


  // Known special values
  if(M.isNaN(z)) return NaN;
  if(M.isInf(z)) return M(Infinity); // any complex infinity gives positive infinity
  if(branch == 0 && z.eqr(0)) return M(0);
  if(branch != 0 && z.eqr(0)) return M(-Infinity); //at all other branch than the principal one, it's -infinity at 0

  /*
  Choosing a good starting value is important. M(0) as starting value works
  most of the time, but does not work at some regions in the negative complex domain,
  e.g. around 5.4+0.1i, 5.5+0.1i, ... and that can be seen as mandelbrot-fractal-like
  circles around those regions in the complex domain plot.
  */
  var w = M.log(z).add(M(0, branch * Math.PI * 2));
  if(branch == 0 && z.abs() < 1.2 /*supposed to be 1/Math.E, but I still see problems beyond that in the complex domain plot...*/) {
    w = M.sqrt(z.mulr(5.43656365691809047).addr(2)).add(M(-1, branch * Math.PI * 2)); //TODO: verify if this where ctor gets two arguments works correctly for quaternions?
  }
  if(branch != 0 && z.im == 0) z.im += 1e-14; // Give it small imaginary part, otherwise it never gets there // TODO: this does not work correctly for quaternions

  var num = 36;
  for(var i = 0; i < num; i++) {
    var ew = M.exp(w);
    var wew = w.mul(ew);
    var t = wew.sub(z);
    var a = ew.mul(w.addr(1));
    var b = w.addr(2).mul(t).div(w.mulr(2).addr(2));
    w = w.sub(t.div(a.sub(b)));

    var ltest = M.log(z.div(w)); //for testing if near (z = w*exp(w) OR ln(z/w) = w)
    if(M.near(ltest, w, 1e-16) || M.near(wew, z, 1e-16)) break;
    if(i + 1 == num && !(M.near(ltest, w, 1) || M.near(wew, z, 1))) return M(NaN); // iteration could not finish and too far from result
  }

  // Remove numeric tiny imaginary part that appeared in error
  if(z.im == 0 && z.re >= 0) w.im = 0;

  return w;
};

// Lambertw for branch (0 = principal branch Wp, -1 is also common (Wm))
// Branch is real integer, z is Jmat.Complex object (complex)
Jmat.Complex.lambertwb = function(branch, z) {
  return Jmat.Complex.lambertwb_generic_(Jmat.Complex, branch, z);
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

////////////////////////////////////////////////////////////////////////////////

// Faddeeva function, used as helper functions to calculate erf and related functions for certain locations in the complex plane
// Faddeeva(z) = exp(-z^2)*erfc(-iz).
// Also known as Faddeyeva, or as w(x), but that may be confusing with LambertW...
Jmat.Complex.faddeeva = function(z) {
  var f = Jmat.Real.faddeeva(z.re, z.im);
  return Jmat.Complex(f[0], f[1]);
};

// erfcx(z) = exp(z^2) * erfc(z): the scaled complementary error function
Jmat.Complex.erfcx = function(z) {
  return Jmat.Complex.faddeeva(Jmat.Complex(-z.im, z.re)); //erfcx(z) = faddeeva(iz)
};

Jmat.Complex.erf = function(z) {
  var a = Jmat.Complex.exp(z.mul(z).neg()); // If abs of z is very large, and |im| > |re|, then this becomes some NaN or Infinity. That is ok, erf is also some unrepresentable huge value there.
  var result;
  if (z.re >= 0) result = Jmat.Complex.ONE.sub(a.mul(Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I))));
  else return result = a.mul(Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I.neg()))).sub(Jmat.Complex.ONE);
  // fix numerical imprecisions in case something is known to be exactly zero
  if(z.re == 0) result.re = 0; // pure imaginary input must give pure imaginary output, but due to subtracting from 1 it may be near-zero
  return result;

  // With integration, don't use.
  /*var ps2 = 2.0 / Jmat.Real.SQRTPI;
  var result;
  result = Jmat.Complex.integrate(Jmat.Complex(0), z, function(z){ return Jmat.Complex.exp(z.mul(z).neg()); }, 100);
  result = result.mulr(ps2);
  return result;*/
};

// erfc(x) = 1 - erf(x). This function gives numerically a better result if erf(x) is near 1.
Jmat.Complex.erfc = function(z) {
  var a = Jmat.Complex.exp(z.mul(z)); // we divide through this, rather than taking exp(-x*x) amd multiplying, to match what faddeeva does better, to ensure getting '1' for the real part in case of pure imaginary input
  if (z.re >= 0) return Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I)).div(a);
  else return Jmat.Complex.TWO.sub(Jmat.Complex.faddeeva(z.mul(Jmat.Complex.I.neg())).div(a));
};


// TODO: rewrite some of the rational approximations to not use so many multiplications
//a + bx + cxx + dxxx + exxxx = a + x * (b + x * (c + x * (d + x * e)))   etc...
//and that equals: x.mulr(e).addr(d).mul(x).addr(c).mul(x).addr(b).mul(x).addr(a) ...


//erfi(z) = -i erf(iz)
Jmat.Complex.erfi = function(z) {
  return Jmat.Complex.erf(z.mul(Jmat.Complex.I)).mul(Jmat.Complex.I).neg();
};

// D+(x) aka F(x)
Jmat.Complex.dawson = function(z) {
  var w = Jmat.Complex.faddeeva(z);
  var a = Jmat.Complex.exp(z.mul(z).neg());
  return a.sub(w).mul(Jmat.Complex.I.mulr(Jmat.Real.SQRTPI / 2));
};

////////////////////////////////////////////////////////////////////////////////

// gives a random complex number by default inside the unit circle, or else with absolute value between r0 and r1
Jmat.Complex.random = function(r0, r1) {
  r0 = (r0 == undefined) ?  0 : r0;
  r1 = (r1 == undefined) ?  1 : r1;
  return Jmat.Complex.polar(Math.random() * (r1 - r0) + r0, Math.random() * Math.PI * 2);
};



////////////////////////////////////////////////////////////////////////////////
/** @license
License of Jmat.Complex.kiss_fft_ (converted from C to JavaScript):
Kiss FFT
Copyright (c) 2003-2010, Mark Borgerding
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
Jmat.Complex.kiss_fft_ = function(fin/*Jmat.Complex array of size nfft*/, fout/*Jmat.Complex array of size nfft*/, inverse) {
  // Constructor
  var KissFFTState = function(nfft, inverse) {
    this.nfft = nfft;
    this.inverse = inverse;
    this.factors = []; //int array, size 32
    this.twiddles = []; //complex Jmat.Complex array, size nfft-1
  };
  var kf_bfly2_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
    var j = 0;
    for(var i = 0; i < m; i++) {
      var t = Fout[Fout_index + i + m].mul(st.twiddles[j]);
      j += fstride;
      Fout[Fout_index + i + m] = Fout[Fout_index + i].sub(t);
      Fout[Fout_index + i] = Fout[Fout_index + i].add(t);
    }
  };
  var kf_bfly4_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
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
  var kf_bfly3_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
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
  var kf_bfly5_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/) {
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
  var kf_bfly_generic_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, fstride /*int*/, st /*kiss_fft_state*/, m /*int*/, p /*int*/) {
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
  var kf_work_ = function(Fout /*array of complex Jmat.Complex*/, Fout_index, f /*array of complex Jmat.Complex*/,f_index,
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
        kf_work_(Fout, Fout_index + i, f, f_index + j, fstride*p, in_stride, factors, factors_index + 2, st);
        j += fstride*in_stride;
      }
    }
    // recombine the p smaller DFTs
    switch (p) {
      case 2: kf_bfly2_(Fout,Fout_index,fstride,st,m); break;
      case 3: kf_bfly3_(Fout,Fout_index,fstride,st,m); break;
      case 4: kf_bfly4_(Fout,Fout_index,fstride,st,m); break;
      case 5: kf_bfly5_(Fout,Fout_index,fstride,st,m); break;
      default: kf_bfly_generic_(Fout,Fout_index,fstride,st,m,p); break;
    }
  };
  // facbuf is populated by p1,m1,p2,m2, ... where p[i] * m[i] = m[i-1], m0 = n
  var kf_factor_ = function(n, facbuf) {
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
  var kiss_fft_alloc_ = function(nfft,inverse_fft) {
    var st = new KissFFTState(nfft, inverse_fft);
    for (var i=0;i<nfft;++i) {
      var phase = -Math.PI*2*i / nfft;
      if (st.inverse) phase *= -1;
      st.twiddles[i] = new Jmat.Complex(Math.cos(phase), Math.sin(phase));
    }
    kf_factor_(nfft,st.factors);
    return st;
  };
  var kiss_fft_ = function(st/*kiss_fft_state object*/,fin/*complex Jmat.Complex array of size nfft*/,fout/*complex Jmat.Complex array of size nfft*/) {
    kf_work_(fout, 0, fin, 0, 1, 1/*in_stride*/, st.factors, 0, st);
  };
  var st = kiss_fft_alloc_(fin.length, inverse);
  kiss_fft_(st, fin, fout);
};
// End of Kiss FFT
////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.fft = function(a, inverse) {
  var out = [];
  for(var i = 0; i < a.length; i++) out[i] = Jmat.Complex(0);
  Jmat.Complex.kiss_fft_(a, out, inverse);
  var factor = 1.0 / Math.sqrt(a.length);
  for(var i = 0; i < a.length; i++) out[i] = out[i].mulr(factor);
  return out;
};

