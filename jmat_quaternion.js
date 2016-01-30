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

// REQUIRES: jmat_real.js, jmat_complex.js

/*
Jmat.Quaternion: quaternion mathematics

Overview of some functionality:
-elementary operators: Quaternion.add, Quaternion.sub, Quaternion.mul, Quaternion.div, Quaternion.rightdiv, Quaternion.inv
-mathematical functions: Quaternion.abs, Quaternion.pow, Quaternion.exp, Quaternion.log, Quaternion.log2, Quaternion.log10, Quaternion.cos, Quaternion.acos, Quaternion.cosh, ...
-special functions: Quaternion.lambertw
*/

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Jmat.Quaternion
////////////////////////////////////////////////////////////////////////////////

/*
Constructor, but also usable without new as factory function.

Class representing a quaternion value.

Aliased as simply "Quaternion" at the end of the file - disable that if it causes name clashes

The quaternion value is w + x*i + y*j + z*k with i*i = j*j = k*k = i*j*k = -1
*/
Jmat.Quaternion = function(w, x, y, z) {
  if(this instanceof Jmat.Quaternion) {
    // Keyword "new" in front. Does not do any checks, to be "fast"
    this.w = w; // real part, aka scalar part, aka s or t
    this.x = x; // vector x
    this.y = y; // vector y
    this.z = z; // vector z
  } else {
    // No keyword "new" in front, use the convenience factory function instead
    return Jmat.Quaternion.make(w, x, y, z); // This supports several argument types
  }
};

/*
Can make a quaternion from:
-nothing
-a single real number
-up to 4 real numbers for all components
-string representation
-one or two complex numbers
-2x2 matrix (as 2D array or object with array in e field (Jmat.Matrix))
-4x4 matrix (as 2D array or object with array in e field (Jmat.Matrix))
-another quaternion
*/
Jmat.Quaternion.make = function(w, x, y, z) {
  if(w == undefined) return new Jmat.Quaternion(0, 0, 0, 0);
  if(typeof w == 'number') return new Jmat.Quaternion(w, x || 0, y || 0, z || 0);
  if(typeof w == 'string') return Jmat.Quaternion.parse(w);
  if(w.re != undefined) return new Jmat.Quaternion(w.re, w.im || 0, (x && x.re) || 0, (x && x.im) || 0);
  if(w.e && w.w == 2 && w.h == 2) return Jmat.Quaternion.from2x2(w);
  if(w.e && w.w == 4 && w.h == 4) return Jmat.Quaternion.from4x4(w);
  return Jmat.Quaternion.copy(w);
};


Jmat.Quaternion.cast = function(v) {
  if(v && v.w != undefined) return v;
  if(v == undefined) return Jmat.Quaternion(0);
  return Jmat.Quaternion(v);
};

//aka clone
Jmat.Quaternion.copy = function(v) {
  return new Jmat.Quaternion(v.w, v.x, v.y, v.z);
};

Jmat.Quaternion.toString = function(value, opt_precision) {
  if(!value) return value == 0 ? 'invalid0' : ('' + value);
  var w = (opt_precision ? Jmat.Complex.formatFloat_(value.w, opt_precision) : ('' + value.w));
  var x = (opt_precision ? Jmat.Complex.formatFloat_(value.x, opt_precision) : ('' + value.x));
  var y = (opt_precision ? Jmat.Complex.formatFloat_(value.x, opt_precision) : ('' + value.y));
  var z = (opt_precision ? Jmat.Complex.formatFloat_(value.x, opt_precision) : ('' + value.z));

  var result = '';

  if(value.w != 0) result += w;
  if(value.x != 0) result += ((result.length == 0 || value.x < 0) ? (x) : ('+' + x)) + 'i';
  if(value.y != 0) result += ((result.length == 0 || value.y < 0) ? (y) : ('+' + y)) + 'j';
  if(value.z != 0) result += ((result.length == 0 || value.z < 0) ? (z) : ('+' + z)) + 'k';

  return result;
};
Jmat.Quaternion.prototype.toString = function(opt_precision) {
  return Jmat.Quaternion.toString(this, opt_precision);
};

// Parses strings of the form '5', '5+i+2j+3k', '5-2.3i', '1.25e-25+17.37e5i-2.31k', 'i + j - k'
Jmat.Quaternion.parse = function(text) {
  text = text.replace(/ /g, ''); //no whitespace
  var s = text.split(/\b/);
  var t = [];
  var e = false;
  for(var i = 0; i < s.length; i++) {
    var p = s[i];
    if(i == 0 || (!e && (p[0] == '-' || p[0] == '+'))) t.push('');
    t[t.length - 1] += p;
    e = (p[p.length - 1] == 'e' || p[p.length - 1] == 'E');
  }

  var w = 0;
  var x = 0;
  var y = 0;
  var z = 0;

  for(var i = 0; i < t.length; i++) {
    var p = t[i];
    var h = p[p.length - 1];
    if(h == 'i' || h == 'j' || h == 'k') {
      p = p.substr(0, p.length - 1);
      if(p == '' || p == '+') p = '1';
      if(p == '-') p = '-1';
      if(h == 'i') x = parseFloat(p);
      if(h == 'j') y = parseFloat(p);
      if(h == 'k') z = parseFloat(p);
    } else {
      w = parseFloat(p);
    }

  }

  return new Jmat.Quaternion(w, x, y, z);
};

//Basic operators

Jmat.Quaternion.add = function(x, y) {
  return new Jmat.Quaternion(x.w + y.w, x.x + y.x, x.y + y.y, x.z + y.z);
};
Jmat.Quaternion.prototype.add = function(y) {
  return new Jmat.Quaternion(this.w + y.w, this.x + y.x, this.y + y.y, this.z + y.z);
};

Jmat.Quaternion.sub = function(x, y) {
  return new Jmat.Quaternion(x.w - y.w, x.x - y.x, x.y - y.y, x.z - y.z);
};
Jmat.Quaternion.prototype.sub = function(y) {
  return new Jmat.Quaternion(this.w - y.w, this.x - y.x, this.y - y.y, this.z - y.z);
};

// Hamilton product
Jmat.Quaternion.mul = function(x, y) {
  var rw = x.w * y.w - x.x * y.x - x.y * y.y - x.z * y.z;
  var rx = x.w * y.x + x.x * y.w + x.y * y.z - x.z * y.y;
  var ry = x.w * y.y - x.x * y.z + x.y * y.w + x.z * y.x;
  var rz = x.w * y.z + x.x * y.y - x.y * y.x + x.z * y.w;
  return new Jmat.Quaternion(rw, rx, ry, rz);
};
Jmat.Quaternion.prototype.mul = function(y) {
  return Jmat.Quaternion.mul(this, y);
};

Jmat.Quaternion.mulr = function(x, y) {
  return new Jmat.Quaternion(x.w * y, x.x * y, x.y * y, x.z * y);
};
Jmat.Quaternion.prototype.mulr = function(y) {
  return new Jmat.Quaternion(this.w * y, this.x * y, this.y * y, this.z * y);
};

Jmat.Quaternion.divr = function(x, y) {
  return new Jmat.Quaternion(x.w / y, x.x / y, x.y / y, x.z / y);
};
Jmat.Quaternion.prototype.divr = function(y) {
  return new Jmat.Quaternion(this.w / y, this.x / y, this.y / y, this.z / y);
};

Jmat.Quaternion.addr = function(x, y) {
  return new Jmat.Quaternion(x.w + y, x.x, x.y, x.z);
};
Jmat.Quaternion.prototype.addr = function(y) {
  return new Jmat.Quaternion(this.w + y, this.x, this.y, this.z);
};

Jmat.Quaternion.subr = function(x, y) {
  return new Jmat.Quaternion(x.w - y, x.x, x.y, x.z);
};
Jmat.Quaternion.prototype.subr = function(y) {
  return new Jmat.Quaternion(this.w - y, this.x, this.y, this.z);
};

Jmat.Quaternion.rsub = function(x, y) {
  return new Jmat.Quaternion(y - x.w, -x.x, -x.y, -x.z);
};
Jmat.Quaternion.prototype.rsub = function(y) {
  return new Jmat.Quaternion(y - this.w, -this.x, -this.y, -this.z);
};

Jmat.Quaternion.neg = function(q) {
  return Jmat.Quaternion(-q.w, -q.x, -q.y, -q.z);
};
Jmat.Quaternion.prototype.neg = function() {
  return Jmat.Quaternion(-this.w, -this.x, -this.y, -this.z);
};

Jmat.Quaternion.conj = function(q) {
  return Jmat.Quaternion(q.w, -q.x, -q.y, -q.z);
};
Jmat.Quaternion.prototype.conj = function() {
  return Jmat.Quaternion(this.w, -this.x, -this.y, -this.z);
};

// absolute value, aka norm or modulus, as a Jmat.Quaternion object (with zero imaginary parts)
Jmat.Quaternion.abs = function(q) {
  return Jmat.Quaternion(q.abs());
};
// returned as real (regular JS number)
Jmat.Quaternion.prototype.abs = function() {
  return Math.sqrt(this.abssq());
};

// absolute value squared
Jmat.Quaternion.abssq = function(q) {
  return Jmat.Quaternion(q.abssq());
};
Jmat.Quaternion.prototype.abssq = function() {
  if(this.w == -Infinity || this.x == -Infinity || this.y == -Infinity || this.z == -Infinity) {
    return Infinity;
  }

  return this.w * this.w + this.x * this.x + this.y * this.y + this.z * this.z;
};

// dist, cheb and manhattan all return regular real JS numbers for all types. In some types they are all the same, but not for e.g. Complex or Matrix.
// Euclidean distance
Jmat.Quaternion.dist = function(a, b) {
  return a.sub(b).abs();
};
//Chebyshev distance
Jmat.Quaternion.cheb = function(a, b) {
  return Math.max(Math.max(Jmat.Real.dist(a.w, b.w), Jmat.Real.dist(a.x, b.x)),
                  Math.max(Jmat.Real.dist(a.y, b.y), Jmat.Real.dist(a.z, b.z)));
};
// Manhattan distance
Jmat.Quaternion.manhattan = function(a, b) {
  return Math.abs(a.w - b.w) + Math.abs(a.x - b.x) + Math.abs(a.y - b.y) + Math.abs(a.z - b.z);
};

// absolute value of the vector part (returned as Quaternion object)
Jmat.Quaternion.absv = function(q) {
  return Jmat.Quaternion(q.absv());
};
// returned as real (regular JS number)
Jmat.Quaternion.prototype.absv = function() {
  return Math.sqrt(this.absvsq());
};

// absolute value of vector part squared
Jmat.Quaternion.absvsq = function(q) {
  return Jmat.Quaternion(q.abssq());
};
Jmat.Quaternion.prototype.absvsq = function() {
  if(this.x == -Infinity || this.y == -Infinity || this.z == -Infinity) {
    return Infinity;
  }

  return this.x * this.x + this.y * this.y + this.z * this.z;
};

// argument: atan2(|v|, w) where v is the vector part
Jmat.Quaternion.arg = function(q) {
  return Jmat.Quaternion(q.arg());
};
Jmat.Quaternion.prototype.arg = function() {
  return Math.atan2(this.absv(), this.w);
};

// inverse aka reciproke
Jmat.Quaternion.inv = function(q) {
  return q.conj().divr(q.abssq());
};
Jmat.Quaternion.prototype.inv = function() {
  return this.conj().divr(this.abssq());
};


// returns a/b = a * b^-1
Jmat.Quaternion.div = function(a, b) {
  return a.mul(b.inv());
};
Jmat.Quaternion.prototype.div = function(b) {
  return this.mul(b.inv());
};


// returns a/b = b^-1 * a
Jmat.Quaternion.leftdiv = function(a, b) {
  return a.inv().mul(b);
};
Jmat.Quaternion.prototype.leftdiv = function(b) {
  return this.inv().mul(b);
};

//aka versor aka normalize
Jmat.Quaternion.sign = function(q) {
  var a = q.abs();
  if(a == 0) return q;
  return q.divr(a);
};

Jmat.Quaternion.normalize = Jmat.Quaternion.sign;
Jmat.Quaternion.prototype.normalize = function() {
  return Jmat.Quaternion.sign(this);
};

// Convert to 2x2 matrix representation (as 2D array, can be given to Jmat.Matrix ctor)
Jmat.Quaternion.to2x2 = function(q) {
  var C = Jmat.Complex;
  return [[C(q.w, q.x), C(q.y, q.z)],
          [C(-q.y, q.z), C(q.w, -q.x)]];
};

// Convert from 2x2 matrix to quaternion (matrix as 2D array or object with array in e field (Jmat.Matrix))
Jmat.Quaternion.from2x2 = function(m) {
  var e = m.e ? m.e : m
  var e0 = e[0][0];
  var e1 = e[0][1];
  return new Jmat.Quaternion(e0.re, e0.im, e1.re, e1.im)
};

// Convert to 4x4 matrix representation (as 2D array, can be given to Jmat.Matrix ctor)
Jmat.Quaternion.to4x4 = function(q) {
  return [[q.w, q.x, q.y, q.z],
          [-q.x, q.w, -q.z, q.y],
          [-q.y, q.z, q.w, -q.x],
          [-q.z, -q.y, q.x, q.w]];
};

// Convert from 4x4 matrix to quaternion (matrix as 2D array or object with array in e field (Jmat.Matrix))
Jmat.Quaternion.from4x4 = function(m) {
  var e = m.e ? m.e : m
  return new Jmat.Quaternion(e[0][0].re, e[0][1].re, e[0][2].re, e[0][3].re);
};

// Convert to 3x3 rotation matrix (with column vectors) (as 2D array, can be given to Jmat.Matrix ctor). q must be unary quaternion.
Jmat.Quaternion.to3x3rot = function(q) {
  var aa = q.w * q.w;
  var bb = q.x * q.x;
  var cc = q.y * q.y;
  var dd = q.z * q.z;
  return [[aa + bb - cc - dd, 2 * q.x * q.y - 2 * q.w * q.z, 2 * q.x * q.z + 2 * q.w * q.y],
          [2 * q.x * q.y + 2 * q.w * q.z, aa - bb + cc - dd, 2 * q.y * q.z - 2 * q.w * q.x],
          [2 * q.x * q.z - 2 * q.w * q.y, 2 * q.y * q.z + 2 * q.w * q.x, aa - bb - cc + dd]];
};

// Convert from 3x3 rotation matrix (with column vectors) (matrix as 2D array or object with array in e field (Jmat.Matrix))
Jmat.Quaternion.from3x3rot = function(m) {
  var e = m.e ? m.e : m;
  var t = 1 + e[0][0].re + e[1][1].re + e[2][2].re;
  if(t > 1e-7) {
    var s = Math.sqrt(t) * 2;
    return new Jmat.Quaternion(s * 0.25, (e[2][1].re - e[1][2].re) / s,
        (e[0][2].re - e[2][0].re) / s, (e[1][0].re - e[0][1].re) / s);
  }
  if(e[0][0].re > e[1][1].re && e[0][0].re > e[2][2].re) {
    var s = Math.sqrt(1 + e[0][0].re - e[1][1].re - e[2][2].re) * 2;
    return new Jmat.Quaternion((e[2][1].re - e[1][2].re) / s, s * 0.25,
        (e[1][0].re + e[0][1].re) / s, (e[0][2].re - e[2][0].re) / s);
  }
  if(e[1][1].re > e[2][2].re) {
    var s = Math.sqrt(1 + e[1][1].re - e[0][0].re - e[2][2].re) * 2;
    return new Jmat.Quaternion((e[0][2].re - e[2][0].re) / s, (e[1][0].re + e[0][1].re) / s,
        0.25 * s, (e[2][1].re + e[1][2].re) / s);
  }
  var s = Math.sqrt(1 + e[2][2].re - e[0][0].re - e[1][1].re) * 2;
  return new Jmat.Quaternion((e[1][0].re - e[0][1].re) / s, (e[0][2].re + e[2][0].re) / s,
      (e[2][1].re + e[1][2].re) / s, 0.5 * s);
};

// Gets the 3D vector as column vector (matrix as 2D array, can be given to Jmat.Matrix ctor)
Jmat.Quaternion.prototype.getVector = function() {
  return [[this.x], [this.y], [this.z]];
};

// Gets the 3D vector as quaternion with scalar part 0
Jmat.Quaternion.prototype.qvector = function() {
  return new Jmat.Quaternion(0, this.x, this.y, this.z);
};

Jmat.Quaternion.exp = function(q) {
  var w = Math.exp(q.w);
  var v = q.absv();
  var cv = Math.cos(v);
  var sv = Math.sin(v);
  var sva = w * sv / v;
  return new Jmat.Quaternion(w * cv, sva * q.x, sva * q.y, sva * q.z);
};

//natural logarithm (ln)
Jmat.Quaternion.log = function(q) {
  var n = q.abs();
  var v = q.absv();
  if (v == 0) return new Jmat.Quaternion(Math.log(q.w), 0, 0, 0);
  var a = Math.acos(q.w / n);
  var av = a / v;
  return new Jmat.Quaternion(Math.log(n), av * q.x, av * q.y, av * q.z);
};

//arbitrary log: log_y(x), y is also quaternion
//warning: base y is second argument
Jmat.Quaternion.logy = function(x, y) {
  return Jmat.Quaternion.log(x).div(Jmat.Quaternion.log(y));
};

//arbitrary log: log_y(x), where y is a regular JS number
//warning: base y is second argument
Jmat.Quaternion.logr = function(x, y) {
  return Jmat.Quaternion.log(x).divr(Math.log(y));
};

Jmat.Quaternion.log2 = function(q) {
  return Jmat.Quaternion.log(q).divr(Math.LN2);
};

Jmat.Quaternion.log10 = function(q) {
  return Jmat.Quaternion.log(q).divr(Math.LN10);
};

Jmat.Quaternion.pow = function(x, y) {
  var Q = Jmat.Quaternion;
  return Q.exp(Q.log(x).mul(y));
};
Jmat.Quaternion.prototype.pow = function(y) {
  return Jmat.Quaternion.pow(this, y);
};

Jmat.Quaternion.powr = function(x, y) {
  var Q = Jmat.Quaternion;
  return Q.exp(Q.log(x).mulr(y));
};
Jmat.Quaternion.prototype.powr = function(y) {
  return Jmat.Quaternion.powr(this, y);
};

Jmat.Quaternion.sqrt = function(q) {
  return q.powr(0.5);
};

Jmat.Quaternion.eq = function(x, y) {
  if(!x || !y) return x == y;
  return (x.w == y.w && x.x == y.x && x.y == y.y && x.z == y.z);
};
Jmat.Quaternion.prototype.eq = function(y) {
  return y && this.w == y.w && this.x == y.x && this.y == y.y && this.z == y.z;
};

// equal to real number? x is quaternion, y is JS number
Jmat.Quaternion.eqr = function(x, y) {
  if(!x || !y) return x == y;
  return (x.w == y && x.x == 0 && x.y == 0 && x.z == 0);
};
Jmat.Quaternion.prototype.eqr = function(y) {
  return this.w == y && this.x == 0 && this.y == 0 && this.z == 0;
};

Jmat.Quaternion.near = function(x, y, epsilon) {
  return x.w - epsilon <= y.w && x.w + epsilon >= y.w &&
         x.x - epsilon <= y.x && x.x + epsilon >= y.x &&
         x.y - epsilon <= y.y && x.x + epsilon >= y.y &&
         x.z - epsilon <= y.z && x.x + epsilon >= y.z;
};

// See Jmat.Complex.relnear
Jmat.Quaternion.relnear = function(x, y, precision) {
  if(x.eq(y)) return true;
  return x.sub(y).abs() < (Math.max(x.abs(), y.abs()) * precision);
};

Jmat.Quaternion.ONE = new Jmat.Quaternion(1, 0, 0, 0);
Jmat.Quaternion.I = new Jmat.Quaternion(0, 1, 0, 0);
Jmat.Quaternion.J = new Jmat.Quaternion(0, 0, 1, 0);
Jmat.Quaternion.K = new Jmat.Quaternion(0, 0, 0, 1);

Jmat.Quaternion.sin = function(q) {
  var v = q.absv();
  var sc = Math.sin(q.w) * Jmat.Real.cosh(v);
  var cs = Math.cos(q.w) * Jmat.Real.sinh(v);
  return new Jmat.Quaternion(sc, cs * q.x / v, cs * q.y / v, cs * q.z / v);
};

Jmat.Quaternion.cos = function(q) {
  var v = q.absv();
  if(v == 0) return Jmat.Quaternion(Math.cos(q.w), 0, 0, 0);
  var cc = Math.cos(q.w) * Jmat.Real.cosh(v);
  var ss = Math.sin(q.w) * Jmat.Real.sinh(v);
  return new Jmat.Quaternion(cc, -ss * q.x / v, -ss * q.y / v, -ss * q.z / v);
};

Jmat.Quaternion.tan = function(q) {
  var Q = Jmat.Quaternion;
  return Q.sin(q).div(Q.cos(q));
};

Jmat.Quaternion.asin = function(q) {
  var v = (q.absvsq() == 0) ? Jmat.Quaternion.I : q.qvector().divr(q.absv());
  return v.mul(Jmat.Quaternion.asinh(q.mul(v))).neg();
};

Jmat.Quaternion.acos = function(q) {
  var v = (q.absvsq() == 0) ? Jmat.Quaternion.I : q.qvector().divr(q.absv());
  return v.mul(Jmat.Quaternion.acosh(q)).neg();
};

Jmat.Quaternion.atan = function(q) {
  var v = (q.absvsq() == 0) ? Jmat.Quaternion.I : q.qvector().divr(q.absv());
  return v.conj().mul(Jmat.Quaternion.atanh(q.mul(v)));
};

Jmat.Quaternion.sinh = function(z) {
  var e = Jmat.Quaternion.exp(z);
  var ei = Jmat.Quaternion.inv(e);
  return e.sub(ei).divr(2);
};

Jmat.Quaternion.cosh = function(z) {
  var e = Jmat.Quaternion.exp(z);
  var ei = Jmat.Quaternion.inv(e);
  return e.add(ei).divr(2);
};

Jmat.Quaternion.tanh = function(z) {
  var e = Jmat.Quaternion.exp(z);
  var ei = Jmat.Quaternion.inv(e);
  return e.sub(ei).div(e.add(ei));
};

Jmat.Quaternion.asinh = function(z) {
  return Jmat.Quaternion.log(z.add(Jmat.Quaternion.sqrt(z.mul(z).addr(1))));
};

Jmat.Quaternion.acosh = function(z) {
  // ln(x + sqrt(z-1)*sqrt(z+1))
  return Jmat.Quaternion.log(z.add(Jmat.Quaternion.sqrt(z.subr(1)).mul(Jmat.Quaternion.sqrt(z.addr(1)))));
};

Jmat.Quaternion.atanh = function(z) {
  // 0.5 * (ln(1+z) - ln(1-z))
  return Jmat.Quaternion.log(z.addr(1).div(z.rsub(1))).mulr(0.5);
};

Jmat.Quaternion.isReal = function(z) {
  return z.x == 0 && z.y == 0 && z.z == 0;
};

Jmat.Quaternion.isNaN = function(z) {
  return isNaN(z.w) || isNaN(z.x) || isNaN(z.y) || isNaN(z.z);
};

//is infinite
Jmat.Quaternion.isInf = function(z) {
  return Math.abs(z.w) == Infinity || Math.abs(z.x) == Infinity || Math.abs(z.y) == Infinity || Math.abs(z.z) == Infinity;
};

//isnanorinf isinfornan
Jmat.Quaternion.isInfOrNaN = function(z) {
  return !z || Jmat.Quaternion.isNaN(z) || Jmat.Quaternion.isInf(z);
};

// Lambertw for branch (0 = principal branch Wp, -1 is also common (Wm))
// Branch is real integer, z is Jmat.Quaternion object
Jmat.Quaternion.lambertwb = function(branch, z) {
  // We can use the generic algorithm from Jmat.Complex, because Jmat.Quaternion also has all the operators it needs (mul, isReal, ...)
  return Jmat.Complex.lambertwb_generic_(Jmat.Quaternion, branch, z);
};

// Principal branch of Lambert's W function: Wp, inverse (not reciprocal) of exp(x) * x
Jmat.Quaternion.lambertw = function(z) {
  return Jmat.Quaternion.lambertwb(0, z);
};

// Negative branch of Lambert's W function: Wm, inverse (not reciprocal) of exp(x) * x
Jmat.Quaternion.lambertwm = function(z) {
  // TODO: wrong. Look at the real plot. Fix this! Jmat.plotReal(Jmat.Quaternion.lambertwm)
  return Jmat.Quaternion.lambertwb(-1, z);
};

// TODO: quaternion gamma function
