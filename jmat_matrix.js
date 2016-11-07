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
Jmat.Matrix: matrix and vector math

NOTE: Many routines assume an epsilon of 1e-15, considering values smaller than this to be zero
NOTE: There are also a few matrix algorithms in jmat_real.js. Those work on 2D arrays of real numbers,
      while here we work on a custom object with complex numbers.

Overview of some functionality:
-decompositions: Matrix.lu, Matrix.qr, Matrix.svd, Matrix.evd, Matrix.cholesky, Matrix.ldl
-inverse: Matrix.inv, Matrix.pseudoinverse
-solve: Matrix.solve, Matrix.rref
-eigen: Matrix.eig, Matrix.eig11, Matrix.eig22
-fourier transform: Matrix.fft, Matrix.ifft
-vectors: Matrix.cross, Matrix.dot
-elementary operators: Matrix.add, Matrix.sub, Matrix.mul, Matrix.div, Matrix.leftdiv, Matrix.mulc, Matrix.mulr, ...
-special functions: Matrix.exp, Matrix.sin, Matrix.cos, Matrix.log, Matrix.sqrt, Matrix.powc
-matrix operators: Matrix.transpose, Matrix.conj, Matrix.tranjugate
-determinants and minors: Matrix.determinant, Matrix.minor, Matrix.cofactor, Matrix.adj
-norms and ranks: Matrix.norm, Matrix.maxrownorm, Matrix.maxcolnorm, Matrix.norm2, Matrix.conditionNumber, Matrix.rank, Matrix.trace
-tests: Matrix.isReal, Matrix.isNaN, Matrix.isInfOrNaN, Matrix.eq, Matrix.near
-constructors: Jmat.Matrix, Matrix.make, Matrix.parse, Matrix.cast, Matrix.copy, Matrix.identity, Matrix.zero
-pretty print: Matrix.render, Matrix.toString, Matrix.summary
*/

/*
Constructor, but also usable without new as factory function.

height first because that's the math convention: a 2x3 matrix is 2 rows high, 3 columns wide, and made as new Jmat.Matrix(2, 3).
Does not initialize elements if keyword "new" is in front. If keyword "new" is not in front, then uses Jmat.Matrix.make and all its options to initialize elements instead.

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
Jmat.Matrix([[1, 2], [3, 4]]).toString()   --> [[1, 2], [3, 4]]: row [1, 2] and row [3, 4]
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
*/
Jmat.Matrix.make = function(a, b, var_arg) {
  if(!a) return new Jmat.Matrix(0, 0);
  if(a instanceof Jmat.Matrix) return Jmat.Matrix.copy(a);

  if(typeof a == 'string') return Jmat.Matrix.parse(a);

  if(a.constructor === Array && b != undefined) throw 'no further arguments needed if first is array. Use 2D array if first is array, or otherwise w and h first';

  // Tolerant to all kinds of nonexisting array
  // Also supports a 1D array representing an Nx1 2D array
  var softget = function(a, y, x) {
    return (a && a[y] != undefined) ? Jmat.Complex.cast(a[y][x] == undefined ? a[y] : a[y][x]) : Jmat.Complex();
  };
  var softgetr = function(a, y, x) {
    return (a && a[y] != undefined) ? Jmat.Real.cast(a[y][x] == undefined ? a[y] : a[y][x]) : 0;
  };
  var softget2 = function(a, b, y, x) {
    return new Jmat.Complex(softgetr(a, y, x), softgetr(b, y, x)); // real from a, imag from b
  };
  var arrayw = function(a) {
    if(!a || a[0] == undefined) return 0; // empty array
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
  if((("undefined" !== typeof a) && a.length) || (("undefined" !== typeof b) && b.length)) {
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
Jmat.Matrix.toString = function(m, opt_precision) {
  // e.g in console: Jmat.Matrix.toString(getMatrixFromMem(Jmat.Complex(100))
  if(!m) return '' + m;
  var s = '[';
  for(var y = 0; y < m.h; y++) {
    s += '[';
    for(var x = 0; x < m.w; x++) {
      var e = m.e[y][x];
      s += Jmat.Complex.toString(e, opt_precision);
      if(x + 1 < m.w) s += ', ';
    }
    s += ']';
    if(y + 1 < m.h) s += ', ';
  }
  return s + ']';
};
Jmat.Matrix.prototype.toString = function(opt_precision) {
  return Jmat.Matrix.toString(this, opt_precision);
};

//similar to toString, but using curly braces instead of square brackets
Jmat.Matrix.toCurly = function(m, opt_precision) {
  return Jmat.Matrix.toString(m, opt_precision).replace(/\[/g, '{').replace(/\]/g, '}');
};
Jmat.Matrix.prototype.toCurly = function(opt_precision) {
  return Jmat.Matrix.toCurly(this, opt_precision);
};

//similar to toString, but using matlab/octave notation with semicolons
Jmat.Matrix.toSemi = function(m, opt_precision) {
  if(!m) return '' + m;
  var s = '[';
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      var e = m.e[y][x];
      s += Jmat.Complex.toString(e, opt_precision);
      if(x + 1 < m.w) s += ' ';
    }
    if(y + 1 < m.h) s += '; ';
  }
  return s + ']';
};
Jmat.Matrix.prototype.toSemi = function(opt_precision) {
  return Jmat.Matrix.toSemi(this, opt_precision);
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

// Makes ascii art rendering of the matrix (requires fixed width font)
Jmat.Matrix.render = function(a, opt_precision) {
  if(!a) return '' + a; // e.g. 'null'
  //turn undefined into nan
  if (!Jmat.Matrix.isValid(a)) {
    a = Jmat.Matrix.copy(a);
    for (var y = 0; y < a.h; y++) {
      for (var x = 0; x < a.w; x++) {
        if(!a.e[y][x]) a.e[y][x] = Jmat.Complex(NaN);
      }
    }
  }
  opt_precision = opt_precision == undefined ? (Jmat.Matrix.isInteger(a) ? 0 : 3) : opt_precision;
  var result = '';
  var real = Jmat.Matrix.isReal(a);
  var strings = [];
  var longest = [];
  for (var y = 0; y < a.h; y++) {
    strings.push([]);
    for (var x = 0; x < a.w; x++) {
      if(y == 0) longest.push([0, 0]);
      var s = [a.e[y][x].re.toFixed(opt_precision), Math.abs(a.e[y][x].im).toFixed(opt_precision)];
      longest[x] = [Math.max(s[0].length, longest[x][0]), Math.max(s[1].length, longest[x][1])];
      strings[y].push(s);
    }
  }
  for (var y = 0; y < a.h; y++) {
    var line = '';
    line += '|' + (y + 1 == a.h ? '_' : ' ');
    for (var x = 0; x < a.w; x++) {
      var s = strings[y][x][0];
      while (s.length < longest[x][0]) s = ' ' + s;
      if (!real) {
        var neg = a.e[y][x].im < 0;
        var t = strings[y][x][1];
        while (t.length < longest[x][1]) t = '0' + t;
        t = (neg ? '-' : '+') + t;
        s += t + 'i';
      }

      line += s + ((y + 1 == a.h && x + 1 == a.w) ? '_' : ' ');
    }
    line += '|';
    if(y == 0) {
      var top = ' _';
      while(top.length + 2 < line.length) top += ' ';
      top += '_';
      result += top + '\n';
    }
    result += line + '\n';
  }
  return result;
};
Jmat.Matrix.prototype.render = function(opt_precision) {
  return Jmat.Matrix.render(this, opt_precision);
};


// Does not copy if a is of type Jmat.Matrix.
Jmat.Matrix.cast = function(a) {
  return a instanceof Jmat.Matrix ? a : Jmat.Matrix.make(a);
};

//aka clone
Jmat.Matrix.copy = function(a) {
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if (!a.e[y][x]) result.e[y][x] = a.e[y][x]; //null or undefined...
      else result.e[y][x] = Jmat.Complex(a.e[y][x].re, a.e[y][x].im);
    }
  }
  return result;
};

// Returns new h*w identity matrix. AKA "eye".
Jmat.Matrix.identity = function(h, opt_w) {
  var w = opt_w || h;
  var r = new Jmat.Matrix(h, w);
  for(var y = 0; y < h; y++) {
    for(var x = 0; x < w; x++) {
      r.e[y][x] = Jmat.Complex(x == y ? 1 : 0);
    }
  }
  return r;
};

// Returns new h*w zero matrix
Jmat.Matrix.zero = function(h, opt_w) {
  var w = opt_w || h;
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

// the iterative O(n^3) multiplication algorithm
Jmat.Matrix.n3mul_ = function(a, b) {
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

// the iterative O(n^3) multiplication algorithm
// slow version without cache optimization, left for reference and comparison
Jmat.Matrix.n3mul_nocache_ = function(a, b) {
  if(a.w != b.h) return null; // mathematically invalid
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

// the iterative O(n^3) multiplication algorithm
Jmat.Matrix.n3mul_ = function(a, b) {
  if(a.w != b.h) return null; // mathematically invalid
  var result = new Jmat.Matrix(a.h, b.w);
  var temp = [];
  for (var x = 0; x < b.w; x++) {
    for (var z = 0; z < a.w; z++) temp[z] = b.e[z][x]; // copy for better caching (faster)
    for (var y = 0; y < a.h; y++) {
      var e = Jmat.Complex(0);
      for (var z = 0; z < a.w; z++) e = e.add(a.e[y][z].mul(temp[z]));
      result.e[y][x] = e;
    }
  }
  return result;
};

// Strassen matrix multiplication algorithm
// Measurably faster in JS for 400x400 matrices and higher
Jmat.Matrix.strassen_ = function(a, b) {
  var M = Jmat.Matrix;
  if(a.w != b.h) return null; // mathematically invalid
  if(a.w < 2 || a.h < 2 || b.w < 2) return M.n3mul_(a, b); // doesn't support smaller size than that

  var n = Math.min(a.h, Math.min(a.w, b.w));
  if(n & 1) n--; // we need even size

  var n2 = Math.floor(n / 2);

  var a00 = M.submatrix(a, 0, n2, 0, n2);
  var a01 = M.submatrix(a, 0, n2, n2, n2 * 2);
  var a10 = M.submatrix(a, n2, n2 * 2, 0, n2);
  var a11 = M.submatrix(a, n2, n2 * 2, n2, n2 * 2);
  var b00 = M.submatrix(b, 0, n2, 0, n2);
  var b01 = M.submatrix(b, 0, n2, n2, n2 * 2);
  var b10 = M.submatrix(b, n2, n2 * 2, 0, n2);
  var b11 = M.submatrix(b, n2, n2 * 2, n2, n2 * 2);

  var m0 = (a00.add(a11)).mul(b00.add(b11));
  var m1 = (a10.add(a11)).mul(b00);
  var m2 = a00.mul(b01.sub(b11));
  var m3 = a11.mul(b10.sub(b00));
  var m4 = (a00.add(a01)).mul(b11);
  var m5 = (a10.sub(a00)).mul(b00.add(b01));
  var m6 = (a01.sub(a11)).mul(b10.add(b11));

  var c00 = m0.add(m3).sub(m4).add(m6); // instead of: a00.mul(b00).add(a01.mul(b10));
  var c01 = m2.add(m4);                 // instead of: a00.mul(b01).add(a01.mul(b11));
  var c10 = m1.add(m3);                 // instead of: a10.mul(b00).add(a11.mul(b10));
  var c11 = m0.sub(m1).add(m2).add(m5); // instead of: a10.mul(b01).add(a11.mul(b11));

  var c = M(a.h, b.w);
  M.insertInPlace(c, c00, 0, 0);
  M.insertInPlace(c, c01, 0, c00.w);
  M.insertInPlace(c, c10, c00.h, 0);
  M.insertInPlace(c, c11, c00.h, c00.w);

  // fix dynamic peeling. TODO: this means it's as slow as the n^3 algorithm for parts of very non-square matrices. Implement smarter solution.
  if(n != a.w || n != a.h || n != b.w) {
    var temp = [];
    for (var x = 0; x < b.w; x++) {
      for (var z = 0; z < a.w; z++) temp[z] = b.e[z][x]; // copy for better caching (faster)
      for (var y = 0; y < a.h; y++) {
        var e = Jmat.Complex(0);
        var z0 = 0;
        if(x < n && y < n) z0 = n;
        else c.e[y][x] = Jmat.Complex(0);
        for (var z = z0; z < a.w; z++) e = e.add(a.e[y][z].mul(temp[z]));
        c.e[y][x] = c.e[y][x].add(e);
      }
    }
  }

  return c;
};

Jmat.Matrix.mul = function(a, b) {
  var m = Math.min(a.w, Math.min(a.h, b.w));
  if(m < 350) return Jmat.Matrix.n3mul_(a, b);
  return Jmat.Matrix.strassen_(a, b);
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

//hadamard product: element-wise product
Jmat.Matrix.elmul = function(a, b) {
  if(a.w != b.w || a.h != b.h) return null;
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].mul(b.e[y][x]);
    }
  }
  return result;
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

//element-wise division
Jmat.Matrix.eldiv = function(a, b) {
  if(a.w != b.w || a.h != b.h) return null;
  var result = new Jmat.Matrix(a.h, a.w);

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result.e[y][x] = a.e[y][x].div(b.e[y][x]);
    }
  }
  return result;
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
// Categories and Properties
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

//returns true if any infinity or NaN in matrix. For the rest, must be valid object.
Jmat.Matrix.isInfOrNaN = function(a) {
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if(Jmat.Complex.isInfOrNaN(a.e[y][x])) return true;
    }
  }
  return false;
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Matrix.isSquare = function(a) {
  return a.w == a.h;
};

// Singular matrix: square matrix with determinant zero.
// Other properties: non-invertible, rows or colums linearly dependent, ...
Jmat.Matrix.isSingular = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  var c = Jmat.Matrix.conditionNumber(a);
  // >= is used, not >, so that if epsilon is 0 and c is Infinity, it correctly returns 'true' since c is Infinity for a numerically exact singular matrix.
  return Jmat.Matrix.isSquare(a) && c.abs() >= 1 / epsilon;
};

// Invertible matrix: square matrix with determinant non-zero. AKA Non-singular.
Jmat.Matrix.isInvertible = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  return Jmat.Matrix.isSquare(a) && !Jmat.Matrix.isSingular(a, epsilon);
};

Jmat.Matrix.isIdentity = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = (x == y) ? 1 : 0;
      if (!Jmat.Complex.nearr(a.e[y][x], e, epsilon)) return false;
    }
  }
  return true;
};

Jmat.Matrix.isDiagonal = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if (x == y) continue;
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

Jmat.Matrix.isZero = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

// Equal to its transpose
Jmat.Matrix.isSymmetric = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = y + 1; x < a.w; x++) {
      if (!Jmat.Complex.near(a.e[y][x], a.e[x][y], epsilon)) return false;
    }
  }
  return true;
};

// Equal to its transjugate, complex equivalent of symmetric
// Diagonal elements must be real as they have to be their own complex conjugate
Jmat.Matrix.isHermitian = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = y; x < a.w; x++) {
      if (!Jmat.Complex.near(a.e[y][x], Jmat.Complex.conj(a.e[x][y]), epsilon)) return false;
    }
  }
  return true;
};

// Equal to negative of its transpose
// Diagonal elements must be zero as they have to be their own negative
Jmat.Matrix.isSkewSymmetric = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = y; x < a.w; x++) {
      if (!Jmat.Complex.near(a.e[y][x], Jmat.Complex.neg(a.e[x][y]), epsilon)) return false;
    }
  }
  return true;
};

// AKA anti-hermitian
Jmat.Matrix.isSkewHermitian = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = y; x < a.w; x++) {
      if (!Jmat.Complex.near(a.e[y][x], Jmat.Complex.conj(a.e[x][y]).neg(), epsilon)) return false;
    }
  }
  return true;
};

Jmat.Matrix.isUpperTriangular = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < y; x++) {
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

Jmat.Matrix.isLowerTriangular = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = y + 1; x < a.w; x++) {
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

Jmat.Matrix.isStrictlyUpperTriangular = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x <= y; x++) {
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

Jmat.Matrix.isStrictlyLowerTriangular = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = y; x < a.w; x++) {
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

// upper triangular with ones on the diagonal
Jmat.Matrix.isUpperUnitriangular = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isUpperTriangular(a)) return false;

  for(var i = 0; i < a.h; i++) {
     if (!Jmat.Complex.nearr(a.e[i][i], 1, epsilon)) return false;
  }
  return true;
};

// lower triangular with ones on the diagonal
Jmat.Matrix.isLowerUnitriangular = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isLowerTriangular(a)) return false;

  for(var i = 0; i < a.h; i++) {
     if (!Jmat.Complex.nearr(a.e[i][i], 1, epsilon)) return false;
  }
  return true;
};

// almost triangular: elements right below the diagonal are also allowed to be non-zero
Jmat.Matrix.isUpperHessenberg = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x + 1 < y; x++) {
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

// almost triangular: elements right above the diagonal are also allowed to be non-zero
Jmat.Matrix.isLowerHessenberg = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = y + 2; x < a.w; x++) {
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

// both upper and lower hessenberg
Jmat.Matrix.isTridiagonal = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if (x == y || x + 1 == y || x == y + 1) continue;
      if (!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
    }
  }
  return true;
};

// A*A^T == I, for real A (isUnitary is the complex equivelent). Its transpose is its inverse, and the vectors form an orthonormal basis.
Jmat.Matrix.isOrthogonal = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;
  //if (!Jmat.Matrix.isReal(a)) return false; // Not sure if to be strict with this check or not...

  var aa = a.mul(a.transpose());
  return Jmat.Matrix.isIdentity(aa, epsilon);
};

// A*A^(*) == I
Jmat.Matrix.isUnitary = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  var aa = a.mul(a.transjugate());
  return Jmat.Matrix.isIdentity(aa, epsilon);
};

// Normal matrix: A * A^(*) == A^(*) * A
Jmat.Matrix.isNormal = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  var at = a.transjugate();
  return Jmat.Matrix.near(a.mul(at), at.mul(a), epsilon);
};

// Permutation matrix: binary, exactly one 1 on each row and column
Jmat.Matrix.isPermutation = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y < a.h; y++) {
    var num0 = 0;
    var num1 = 0;
    var numx = 0;
    for(var x = 0; x < a.w; x++) {
      if(!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) num0++;
      else if(!Jmat.Complex.nearr(a.e[y][x], 1, epsilon)) num1++;
      else numx++;
    }
    if (num1 != 1 || numx > 0) return false;
  }

  for(var x = 0; x < a.w; x++) {
    var num0 = 0;
    var num1 = 0;
    var numx = 0;
    for(var y = 0; y < a.h; y++) {
      if(!Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) num0++;
      else if(!Jmat.Complex.nearr(a.e[y][x], 1, epsilon)) num1++;
      else numx++;
    }
    if (num1 != 1 || numx > 0) return false;
  }

  return true;
};

// constant elements along diagonals
Jmat.Matrix.isToeplitz = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 0; y + 1 < a.h; y++) {
    for(var x = 0; x + 1 < a.w; x++) {
      if (!Jmat.Complex.near(a.e[y][x], a.e[y + 1][x + 1], epsilon)) return false;
    }
  }
  return true;
};

// constant elements along anti-diagonals
Jmat.Matrix.isHankel = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  for(var y = 1; y < a.h; y++) {
    for(var x = 0; x + 1 < a.w; x++) {
      if (!Jmat.Complex.near(a.e[y][x], a.e[y - 1][x + 1], epsilon)) return false;
    }
  }
  return true;
};

// like identity except one column may have arbitrary elements below the diagonal
Jmat.Matrix.isFrobenius = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!Jmat.Matrix.isSquare(a)) return false;

  var col = -1; //the one column that non-zero elements below the diagonal

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if (x == y && !Jmat.Complex.nearr(a.e[y][x], 1, epsilon)) return false;
      if (x > y && !Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) return false;
      if (x < y && !Jmat.Complex.nearr(a.e[y][x], 0, epsilon)) {
        if(col >= 0 && x != col) return false;
        col = x;
      }
    }
  }
  return true;
};

Jmat.Matrix.isInteger = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = a.e[y][x];
      if(!Jmat.Real.near(e.im, 0, epsilon)) return false;
      if(Math.abs(Math.round(e.re) - e.re) > epsilon) return false;
    }
  }
  return true;
};

// only 0's and 1's
Jmat.Matrix.isBinary = function(a, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var e = a.e[y][x];
      if(!Jmat.Real.near(e, 0, epsilon) && !Jmat.Real.near(e, 1, epsilon)) return false;
    }
  }
  return true;
};

// Involutory matrix: its own inverse, A*A = I
Jmat.Matrix.isInvolutory = function(a, opt_epsilon) {
  if (!Jmat.Matrix.isSquare(a)) return false;
  return Jmat.Matrix.isIdentity(a.mul(a), opt_epsilon);
};

// Idempotent matrix: A*A = A
Jmat.Matrix.isIdempotent = function(a, opt_epsilon) {
  if (!Jmat.Matrix.isSquare(a)) return false;
  return Jmat.Matrix.near(a, a.mul(a), opt_epsilon);
};

// Nilpotent matrix: A^k = 0 for some positive integer k. May require up to log(n) multiplications, so N^3*log(N) worse complexity but usually much faster due to fast early checks.
Jmat.Matrix.isNilpotent = function(a, opt_epsilon) {
  var M = Jmat.Matrix;
  var epsilon = (opt_epsilon == undefined) ? 1e-15 : opt_epsilon;
  if (!M.isSquare(a)) return false;
  var n = a.w;
  var k = 1;
  for(;;) {
    if(!M.trace(a).eqr(0, epsilon)) return false; // the trace must always be zero, so this is a fast check to quit early
    if(M.isZero(a, epsilon)) return true;
    if(k >= n) break; // it is known that k is <= n, so we can stop checking once over that
    a = a.mul(a);
    k *= 2;
  }
  return false;
};

// Returns an object with various named boolean and scalar properties of the given matrix
// TODO: add parameter to only return fast to calculate properties, max N^2 complexity (so no determinant, svd based condition number, rank, definiteness, ...)
Jmat.Matrix.getProperties = function(a) {
  var M = Jmat.Matrix;
  var result = {};

  result['dimensions'] = '' + a.h + 'x' + a.w;
  result['height'] = a.h;
  result['width'] = a.w;
  result['square'] = M.isSquare(a);
  result['zero'] = M.isZero(a);
  result['real'] = M.isReal(a);
  result['rank'] = M.rank(a);
  result['frobeniusNorm'] = M.norm(a);
  result['spectralNorm'] = M.norm2(a);
  result['conditionNumber'] = M.conditionNumber(a);
  result['nan'] = M.isNaN(a);

  // The following properties only make sense for square matrices
  result['identity'] = M.isIdentity(a);
  result['diagonal'] = M.isDiagonal(a);
  result['tridiagonal'] = M.isTridiagonal(a);
  result['symmetric'] = M.isSymmetric(a);
  result['hermitian'] = M.isHermitian(a);
  result['skewHermitian'] = M.isSkewHermitian(a);
  result['skewSymmetric'] = M.isSkewSymmetric(a);
  result['upperTriangular'] = M.isUpperTriangular(a);
  result['lowerTriangular'] = M.isLowerTriangular(a);
  result['strictlyUpperTriangular'] = M.isStrictlyUpperTriangular(a);
  result['strictlyLowerTriangular'] = M.isStrictlyLowerTriangular(a);
  result['upperUnitriangular'] = M.isUpperUnitriangular(a);
  result['lowerUnitriangular'] = M.isLowerUnitriangular(a);
  result['upperHessenberg'] = M.isUpperHessenberg(a);
  result['lowerHessenberg'] = M.isLowerHessenberg(a);
  result['singular'] = M.isSingular(a);
  result['invertible'] = M.isInvertible(a);
  result['determinant'] = M.determinant(a);
  result['trace'] = M.trace(a);
  result['orthogonal'] = M.isOrthogonal(a);
  result['unitary'] = M.isUnitary(a);
  result['normal'] = M.isNormal(a);
  result['permutation'] = M.isPermutation(a);
  result['toeplitz'] = M.isToeplitz(a);
  result['hankel'] = M.isHankel(a);
  result['frobenius'] = M.isFrobenius(a);
  result['integer'] = M.isInteger(a);
  result['binary'] = M.isBinary(a);
  result['involutory'] = M.isInvolutory(a);
  result['idempotent'] = M.isIdempotent(a);
  if(result['hermitian']) {
    var d = M.definiteness(a);
    if(d == M.INDEFINITE) result['indefinite'] = true;
    else if(d == M.POSITIVE_DEFINITE) result['positiveDefinite'] = result['positiveSemidefinite'] = true;
    else if(d == M.NEGATIVE_DEFINITE) result['negativeDefinite'] = result['negativeSemidefinite'] = true;
    else if(d == M.POSITIVE_SEMI_DEFINITE || result['zero']) result['positiveSemidefinite'] = true;
    else if(d == M.NEGATIVE_SEMI_DEFINITE || result['zero']) result['negativeSemidefinite'] = true;
  }

  return result;
};

// Gives a one-sentence summary of some interesting properites of the matrix. The more properties the matrix has, the longer the sentence (e.g. if it's square more properties appear, ...)
// Does not show redundant properties. E.g. if the matrix is 'identity', will not show 'symmetric', if it's 'normal', will not show 'orthogonal', etc...
// To see every single property instead, do "Jmat.toString(Jmat.Matrix.getProperties(a))"
Jmat.Matrix.summary = function(a) {
  var p = Jmat.Matrix.getProperties(a);

  var toName = function(name) {
    // convert camelCase to lower case with spaces
    name = name.replace(/([A-Z])/g, ' $1').toLowerCase();
    // But keep own names
    name = name.replace('hessenberg', 'Hessenberg');
    name = name.replace('frobenius', 'Frobenius');
    name = name.replace('toeplitz', 'Toeplitz');
    name = name.replace('hankel', 'Hankel');
    name = name.replace('frobenius', 'Frobenius');
    return name;
  };

  //order of non-square related properties
  var nonsquare = ['height', 'width', 'zero', 'real', 'nan',
                   'rank', 'frobeniusNorm', 'spectralNorm', 'conditionNumber', 'integer', 'binary'];
  //order of properties only applicable for square matrices
  var square = ['identity', 'symmetric', 'hermitian', 'skewSymmetric', 'skewHermitian', 'diagonal', 'tridiagonal',
                'upperTriangular', 'lowerTriangular', 'strictlyUpperTriangular', 'strictlyLowerTriangular', 'upperUnitriangular', 'lowerUnitriangular',
                'upperHessenberg', 'lowerHessenberg',
                'singular', 'invertible', 'determinant', 'trace', 'orthogonal', 'unitary', 'normal', 'permutation', 'toeplitz', 'hankel',
                'indefinite', 'positiveDefinite', 'negativeDefinite', 'positiveSemidefinite', 'negativeSemidefinite', 'frobenius', 'involutory', 'idempotent'];

  var opposite = { 'square' : 'non-square', 'real' : 'complex' };
  // these properties are added only to avoid some redundancy in summary output with the "sub" sytem
  p['small2x2'] = (a.w <= 2 && a.h <= 2);
  p['small1x1'] = (a.w <= 1 && a.h <= 1);
  p['realsym'] = p['real'] && p['symmetric'];
  p['realskewsym'] = p['real'] && p['skewSymmetric'];

  // pairs of child:parents, where child is always true if any of the parents is true, with the intention to not display child in a list if parent is already true as it's redundant
  var sub = {
    'strictlyUpperTriangular': ['zero'], 'strictlyLowerTriangular' : ['zero'],
    'upperUnitriangular': ['identity'], 'lowerUnitriangular' : ['identity'],
    'upperTriangular' : ['diagonal', 'strictlyUpperTriangular', 'upperUnitriangular'],
    'lowerTriangular' : ['diagonal', 'frobenius', 'strictlyLowerTriangular', 'lowerUnitriangular'],
    'upperHessenberg' : ['upperTriangular', 'tridiagonal'], 'lowerHessenberg' : ['lowerTriangular', 'tridiagonal'],
    'diagonal' : ['small1x1', 'identity', 'zero'], 'tridiagonal' : ['small2x2', 'diagonal'],
    'orthogonal' : ['normal', 'identity'], 'unitary' : ['normal'], 'normal' : ['identity', 'zero'],
    'hermitian' : ['normal', 'realsym'],  'skewHermitian' : ['realskewsym'],
    'symmetric' : ['diagonal'], 'skewSymmetric' : ['zero'],
    'permutation' : ['identity'], 'invertible' : ['identity'], 'singular' : ['zero'],
    'real' : ['integer'], 'toeplitz' : ['identity', 'zero'], 'hankel' : ['zero'], 'frobenius' : ['identity'],
    'positiveDefinite' : ['identity'], 'negativeSemidefinite' : ['zero', 'negativeDefinite'], 'positiveSemidefinite' : ['zero', 'positiveDefinite'],
    'integer': ['binary'], 'binary': ['identity', 'zero'], 'involutory': ['identity'], 'idempotent': ['identity']
  };

  var summary = p['dimensions'] + ', ' + (p['square'] ? 'square' : opposite['square']);
  for(var i = 0; i < nonsquare.length + square.length; i++) {
    var e = i < nonsquare.length ? nonsquare[i] : square[i - nonsquare.length];
    if (p[e] === true) {
      if (sub[e]) {
        var redundant = false;
        for(var j = 0; j < sub[e].length; j++) if(p[sub[e][j]] === true) redundant = true; //e.g. don't say "upper triangular" if it's already "strictly upper triangular"
        if(redundant) continue;
      }
      summary += ', ' + toName(e);
    }
    if (p[e] === false && opposite[e]) {
      summary += ', ' + toName(opposite[e]);
    }
  }
  var det = p['square'] ? (', determinant ' + p['determinant']) : ', no determinant';
  summary = '' + summary + ' matrix with rank ' + p['rank'] + det + ' and condition number ' + p['conditionNumber'] + '.';

  return summary;
};
Jmat.Matrix.prototype.summary = function() {
  return Jmat.Matrix.summary(this);
};

// render() followed by summary()
Jmat.Matrix.render_summary = function(a) {
  return a.render() + '\n' + a.summary();
};
Jmat.Matrix.prototype.render_summary = function() {
  return Jmat.Matrix.render_summary(this);
};

// TODO: functions like isSymmetric, isHermitian, isDiagonal, ...

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

//transjugate = transposed conjugate, denoted A^* or A^H (also called hermitian transpose or conjugate transpose)
Jmat.Matrix.transjugate = function(a) {
  return Jmat.Matrix.conj(Jmat.Matrix.transpose(a));
};
Jmat.Matrix.prototype.transjugate = function() {
  return Jmat.Matrix.transjugate(this);
};

// Internal algorithm for lu.
// Returns L and U merged into one matrix (without L's diagonal 1's), a pivot array (permutation) with element i the pivot row interchanged with row i, and the parity (0 or 1) of the permutation.
// TODO: usually the permutation format of this function is what you want, not the matrix that lu() returns, so make this public in some way
// TODO: support rectangular matrices
Jmat.Matrix.doolittle_lup_ = function(a) {
  if(a.h != a.w) return null; //must be square
  var M = Jmat.Matrix;
  var C = Jmat.Complex;
  a = M.copy(a); // we'll modify it to interchange rows and insert elements from L and U

  var pivot = [];
  var parity = 0;

  for(var k = 0; k < a.h; k++)
  {
    var y = k;
    var max = a.e[y][k].abs();
    for(var i = k + 1; i < a.h; i++) {
      if (a.e[i][k].abs() > max) {
        max = a.e[i][k].abs();
        y = i;
      }
    }
    pivot[k] = y;

    if(y != k) {
      // pivot, interchange two rows. After this, the row we need is at k, not at y, so continue with k after this
      for(var i = 0; i < a.h; i++) {
        var temp = a.e[y][i];
        a.e[y][i] = a.e[k][i];
        a.e[k][i] = temp;
      }
      parity ^= 1;
    }

    //Returning for singular commented out: still works, resulting U will be singular.
    //if(C.nearr(a.e[k][k], 0, 1e-15)) return null; // singular

    for(var i = k + 1; i < a.h; i++) {
      a.e[i][k] = a.e[i][k].div(a.e[k][k]);
      if(C.isNaN(a.e[i][k])) a.e[i][k] = C(0); // Set 0/0 to 0 for singular input matrix.
    }
    for(var i = k + 1; i < a.h; i++) {
      for(var j = k + 1; j < a.w; j++) {
        a.e[i][j] = a.e[i][j].sub(a.e[i][k].mul(a.e[k][j]));
      }
    }
  }

  return [a, pivot, parity];
};

// LUP decomposition. Returns object {p: P, l: L, u: U} such that A = P*L*U, with P a permutation matrix, L lower triangular with ones on diagonal, U upper triangular
// Uses Doolittle algorithm with row pivoting
// Returns if a is singular.
Jmat.Matrix.lu = function(a) {
  var M = Jmat.Matrix;
  var C = Jmat.Complex;
  var r = M.doolittle_lup_(a);
  if(!r) return r; // error
  a = r[0];
  var pivot = r[1];

  // the above algorithm stored L and U in A, and P in the pivot array
  // now instead turn these into the 3 output matrices

  var l = M.identity(a.w); // the implicit ones on its diagonal were not stored in A above.
  var u = M.zero(a.w);
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      if(x >= y) {
        u.e[y][x] = a.e[y][x];
      } else {
        l.e[y][x] = a.e[y][x];
      }
    }
  }

  var p = M.zero(a.w);
  var pivot2 = []; // list with index of each row in the permutation matrix
  for(var i = 0; i < a.h; i++) pivot2[i] = i;
  for(var i = 0; i < a.h; i++) {
    var temp = pivot2[i]; pivot2[i] = pivot2[pivot[i]]; pivot2[pivot[i]] = temp;
  }
  for(var i = 0; i < a.h; i++) {
    p.e[pivot2[i]][i] = C(1);
  }

  return {p: p, l: l, u: u};
};

// Cholesky decomposition of a into lower triangular matrix and its conjugate transpose
// a must be hermitian and positive definite
// returns l, the lower triangular matrix (a = ll*)
Jmat.Matrix.cholesky = function(a) {
  var M = Jmat.Matrix;
  var C = Jmat.Complex;
  if(!M.isHermitian(a)) return null;

  var l = M.zero(a.w);


  for(var i = 0; i < a.h; i++) {
    for(var j = 0; j <= i; j++) {
      var s = C(0);
      for(var k = 0; k < j; k++) s = s.add(l.e[i][k].mul(l.e[j][k]));
      s = a.e[i][j].sub(s);
      l.e[i][j] = (i == j) ? C.sqrt(s) : l.e[j][j].inv().mul(s);
    }
  }
  return l;
};

// LDL decomposition: similar to cholesky, but A = LDL* with D diagonal matrix and L unitriangular
Jmat.Matrix.ldl = function(a) {
  var M = Jmat.Matrix;
  var C = Jmat.Complex;
  if(!M.isHermitian(a)) return null;

  var l = M.identity(a.w);
  var d = M.zero(a.w);

  for(var i = 0; i < a.h; i++) {
    var v = [];
    for(var j = 0; j < i; j++) {
      v[j] = l.e[i][j].mul(d.e[j][j]);
    }
    v[i] = a.e[i][i];
    for(var j = 0; j < i; j++) {
      v[i] = v[i].sub(l.e[i][j].mul(v[j]));
    }
    d.e[i][i] = v[i];
    for(var j = i + 1; j < a.h; j++) {
      l.e[j][i] = a.e[j][i];
      for(var k = 0; k < i; k++) {
        l.e[j][i] = l.e[j][i].sub(l.e[j][k].mul(v[k]));
      }
      l.e[j][i] = l.e[j][i].div(v[i]);
    }
  }
  return {l: l, d: d};
};



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
  var sign = ((row + col) & 1) == 0 ? 1 : -1;
  return m.mulr(sign);
};

Jmat.Matrix.determinant = function(a) {
  if(a.w != a.h) return NaN; //square matrices only

  if(a.w == 1) return a.e[0][0];
  if(a.w == 2) return a.e[0][0].mul(a.e[1][1]).sub(a.e[0][1].mul(a.e[1][0]));

  // Laplace expansion, this is an O(n!) algorithm, so not used
  /*var result = Jmat.Complex(0);
  for(var x = 0; x < a.w; x++) {
    result = result.add(a.e[0][x].mul(Jmat.Matrix.cofactor(a, 0, x)));
  }*/

  // Calculate with LU decomposition
  var lu = Jmat.Matrix.doolittle_lup_(a);
  var result = Jmat.Complex(1);
  for(var i = 0; i < a.w; i++) {
    result = result.mul(lu[0].e[i][i]);
  }
  if(lu[2] & 1) result = result.neg();


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

induced norm (aka operator norm, matrix p-norm)
------------
1-norm: maximum absolute column sum of the matrix --> Jmat.maxcolnorm
2-norm: largest singular value, aka spectral norm or 2-norm --> Jmat.Matrix.norm2
oo-norm: maximum absolute row sum of the matrix --> Jmat.Matrix.maxrownorm

entrywise norms (vector p-norm)
---------------
entrywise 1-norm: sum of abs of all the elements --> Jmat.Matrix.vectorNorm with p=0
entrywise 2-norm: Frobenius norm, sqrt of sum of squares of the elements --> Jmat.Matrix.norm
entrywise oo-norm: maximum of abs of all the elements, Chebyshev norm --> Jmat.Matrix.chebNorm
arbitrary entrywise norm: --> Jmat.Matrix.vectorNorm
L2,1-norm: sum of Euclidean norms of columns --> Jmat.Matrix.lpqNorm with p=2, q=1

schatten norms
--------------
schatten 1-norm: sum of singular values, aka nuclear norm or trace norm --> Jmat.Matrix.schattenNorm with p = 1
schatten 2-norm: sqrt of sum of squares of singular values, results in same value as Frobenius norm --> Jmat.Matrix.norm
schatten oo-norm: max of the singular values, aka spectral norm or 2-norm --> Jmat.Matrix.norm2

Ky Fan norms
------------
first Ky Fan norm: max of the singular values, aka spectral norm or 2-norm --> Jmat.Matrix.norm2
last Ky Fan norm: sum of singular values, aka nuclear norm or trace norm --> Jmat.Matrix.schattenNorm with p = 1
*/

Jmat.Matrix.sumsq = function(m) {
  var result = 0;
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      var e = m.e[y][x];
      result += e.abssq();
    }
  }
  return Jmat.Complex(result);
};

//Frobenius norm of the matrix (sqrt of sum of squares of modulus of all elements)
//For a vector, this is the Euclidean norm.
//TODO: since usually the more expensive to calculate 2-norm is meant by "the" norm of the matrix, maybe rename this function to "frobeniusnorm" or "frob"?
Jmat.Matrix.norm = function(m) {
  var result = Math.sqrt(Jmat.Matrix.sumsq(m).re);
  return Jmat.Complex(result);
};

// divides through the Frobenius norm
Jmat.Matrix.normalize = function(m) {
  var norm = Jmat.Matrix.norm(m);
  return m.divr(norm.re);
};

//Maximum absolute column sum norm
Jmat.Matrix.maxcolnorm = function(m) {
  var result = 0;
  for(var x = 0; x < m.w; x++) {
    var current = 0;
    for(var y = 0; y < m.h; y++) {
      current += m.e[y][x].abssq();
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
      current += m.e[y][x].abssq();
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

// entrywise norm with arbitrary p (vector p-norm)
// works on all elements of the matrix (it does not have to be a vector)
// NOTE: p must be given as complex number, not regular JS number
// with p = 0 (handwavy): calculates hamming distance of elements to zero
// with p = 1: entriwise manhattan norm
// with p = 2: frobenius norm
// with p = Infinity: maximum absolute value of elements
// with other p: arbitrary p-norms with complex p
Jmat.Matrix.vectorNorm = function(m, p) {
  var C = Jmat.Complex;
  if(C.isReal(p)) {
    var result = 0;
    for(var y = 0; y < m.h; y++) {
      for(var x = 0; x < m.w; x++) {
        var e = m.e[y][x];
        if(p.eqr(0)) {
          if(!C.nearr(e, 0, 1e-15)) result++;
        } else if(p.eqr(1)) {
          result += e.abs();
        } else if(p.eqr(2)) {
          result += e.abssq();
        } else if(p.eqr(Infinity)) {
          result = Math.max(e.abs(), result);
        } else {
          result += Math.pow(e.abssq(), p.re / 2);
        }
      }
    }
    if(result == Infinity && p.re > 0) return Jmat.Matrix.vectorNorm(m, C(Infinity)); // overflow, approximate with max norm instead
    if(p.eqr(2)) {
      result = Math.sqrt(result);
    } else if(!p.eqr(0) && !p.eqr(1) && !p.eqr(Infinity)) {
      result = Math.pow(result, 1 / p.re);
    }
    return C(result);
  } else {
    var result = C(0);
    for(var y = 0; y < m.h; y++) {
      for(var x = 0; x < m.w; x++) {
        var e = C.abssq(m.e[y][x]);
        result = result.add(e.pow(p.divr(2)));
      }
    }
    if(result.eqr(0)) return result;
    return C.pow(result, p.inv());
  }
};

// Lp,q norm, e.g. L2,1 norm for p=2, q=1
Jmat.Matrix.lpqNorm = function(m, p, q) {
  var M = Jmat.Matrix;
  var a = M(1, m.w);
  for(var x = 0; x < m.w; x++) {
    a.e[0][x] = M.vectorNorm(M.col(m, x), p);
  }
  return M.vectorNorm(a, q);
};

// Schatten norm with arbitrary p
// NOTE: p must be given as complex number, not regular JS number
// with p = 1: sum of singular values: nuclear norm or trace norm
// with p = 2: sqrt of sum of squares of singular values, results in same value as Frobenius norm
// with p = Infinity: value of largest singular value
// with other p: arbitrary p-norm of the singular values
Jmat.Matrix.schattenNorm = function(m, p) {
  if(p.eqr(2)) return Jmat.Matrix.norm(m); // not needed to calculate singular values if it's two, as it's the same as frobenius norm
  if(p.eqr(Infinity)) return Jmat.Matrix.norm2(m); // spectral norm
  var M = Jmat.Matrix;
  var svd = M.svd(m);
  var d = M.diagToRow(svd.s);
  return M.vectorNorm(d, p);
};

//Maximum absolute element value
Jmat.Matrix.chebNorm = function(m) {
  var result = 0;
  for(var x = 0; x < m.w; x++) {
    for(var y = 0; y < m.h; y++) {
      result = Math.max(m.e[y][x].abs());
    }
  }
  return Jmat.Complex(result);
};

// dist, cheb and manhattan all return regular real JS numbers for all types. In some types they are all the same, but not for e.g. Complex or Matrix.
// Euclidean distance
Jmat.Matrix.dist = function(a, b) {
  return Jmat.Matrix.norm(a.sub(b)).re;
};
//Chebyshev distance
Jmat.Matrix.cheb = function(a, b) {
  var result = 0; // chebyshev norm, sup norm, max norm or infinity norm of a-b
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result = Math.max(result, Jmat.Complex.cheb(a.e[y][x], b.e[y][x]));
    }
  }
  return result;
};
// Manhattan distance
Jmat.Matrix.manhattan = function(a, b) {
  var result = 0; // chebyshev norm, sup norm, max norm or infinity norm of a-b
  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      result += Jmat.Complex.manhattan(a.e[y][x], b.e[y][x]);
    }
  }
  return result;
};

//condition number: largest singular value divided through smallest singular value. Higher ==> more imprecise numerical calculations with this matrix.
//Infinity means the matrix is singular. For numerical stability, consider singular if above e.g. 1/1e-15
Jmat.Matrix.conditionNumber = function(m) {
  var svd = Jmat.Matrix.svd(m);
  var d = Math.min(m.w, m.h);
  var result = svd.s.e[0][0].div(svd.s.e[d - 1][d - 1]);
  if (Jmat.Complex.isNaN(result)) result = Jmat.Complex(Infinity);
  return result;
};

//Rank of matrix
Jmat.Matrix.rank = function(m, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-14 : opt_epsilon;
  // TODO: use the faster RRQR? Or at least svd that returns only s?
  var s = Jmat.Matrix.svd(m).s;
  var rank = 0;
  var n = Math.min(s.w, s.h);
  for(var i = 0; i < n; i++) {
    if(!Jmat.Real.near(s.e[i][i].re, 0, epsilon)) rank++;
  }
  return Jmat.Complex(rank);
};

// Mathematically only defined for square matrices, but will also return the
// sum of diagonal elements of non-square matrix in this implementation
Jmat.Matrix.trace = function(m) {
  var result = Jmat.Complex.ZERO;
  for(var x = 0; x < m.w && x < m.h; x++) result = result.add(m.e[x][x]);
  return result;
};

// Returns column as h*1 matrix. x is 0-indexed
Jmat.Matrix.col = function(m, x) {
  var r = new Jmat.Matrix(m.h, 1);
  for(var y = 0; y < m.h; y++) r.e[y][0] = m.e[y][x];
  return r;
};

// Returns row as 1*w matrix. y is 0-indexed
Jmat.Matrix.row = function(m, y) {
  var r = new Jmat.Matrix(1, m.w);
  for(var x = 0; x < m.w; x++) r.e[0][x] = m.e[y][x];
  return r;
};

// sets column x of matrix m to the values of c, in-place
// x is 0-indexed, and c must be a h*1 matrix
Jmat.Matrix.setCol = function(m, c, x) {
  for(var y = 0; y < m.h; y++) m.e[y][x] = c.e[y][0];
};

// sets row y of matrix m to the values of r, in-place
// y is 0-indexed, and r must be a 1*w matrix
Jmat.Matrix.setRow = function(m, r, y) {
  for(var x = 0; x < m.w; x++) m.e[y][x] = r.e[0][x];
};

// Add two non-equal sized matrices.
// b's top left element is at position (row, col) in a (0-indexed)
// so the size of the result matrix is:
// max(row + b.h, a.h) - min(0, row) by max(col + b.w, a.w) - min(0, col)
// any element not overlapped by a or b, will be zero.
Jmat.Matrix.overlap = function(a, b, row, col) {
  var h = Math.max(row + b.h, a.h) - Math.min(0, row);
  var w = Math.max(col + b.w, a.w) - Math.min(0, col);

  var result = Jmat.Matrix.zero(h, w);

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

// given b shifted by row,col, insert it into a, overwriting the matching elements of a, leaving other elements of a untouched. Parts of b outside of a, are discarded. The result has the same size as a.
Jmat.Matrix.insertInPlace = function(a, b, row, col) {
  for(var y = 0; y < b.h; y++) {
    for(var x = 0; x < b.w; x++) {
      var rx = x + col;
      var ry = y + row;
      if(rx >= 0 && rx < a.w && ry >= 0 && ry < a.h) {
        a.e[ry][rx] = b.e[y][x];
      }
    }
  }
};

// given b shifted by row,col, insert it into a, overwriting the matching elements of a, leaving other elements of a untouched. Parts of b outside of a, are discarded. The result has the same size as a.
Jmat.Matrix.insert = function(a, b, row, col) {
  var result = Jmat.Matrix.copy(a);
  Jmat.Matrix.insertInPlace(result, b, row, col);
  return result;
};

// similar to insert, but will write outside the matrix if needed, increasing its size
// by default, appends to the right
Jmat.Matrix.augment = function(a, b, opt_row, opt_col) {
  var row = (opt_row == undefined ? 0 : opt_row);
  var col = (opt_col == undefined ? a.w : opt_col);
  var h = Math.max(row + b.h, a.h) - Math.min(0, row);
  var w = Math.max(col + b.w, a.w) - Math.min(0, col);

  var result = Jmat.Matrix.zero(h, w);

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
      result.e[ry][rx] = b.e[y][x];
    }
  }

  return result;
};

// projects v onto u with inner product. Must be vectors.
Jmat.Matrix.proj = function(v, u) {
  var M = Jmat.Matrix;
  var uu = M.sumsq(u);
  if(uu.re == 0) return M.zero(v.h, v.w);
  var vu = M.dot(v, u);
  var f = vu.div(uu);
  return u.mulc(f);
};

// does 1 step of gram-schmidt orthogonalization, for the column x of the matrix. us must have the x previously calculated vectors. The function returns the vector that you can fill in us[x].
Jmat.Matrix.gramSchmidtStep_ = function(m, x, us) {
  var M = Jmat.Matrix;
  var v = M.col(m, x);
  var u = v;

  //for(var j = 0; j < i; j++) u = u.sub(M.proj(v, us[j])); // "classical" way, commented out
  for(var j = 0; j < x; j++) u = u.sub(M.proj(u, us[j])); // more stable way, MGS (only difference is v is changed into previous u iteration)

  return u;
};


// Performs gram-schmidt orthogonalization of the column vectors given in matrix M
// NOTE: This is not for real usage but a demonstration of the Gram-Schmidt process.
// The q returned by Matrix.qr is similar, but calculated with more stable and efficient algorithms (it givens or householder).
Jmat.Matrix.gramSchmidt = function(m, opt_normalize) {
  var M = Jmat.Matrix;
  var vs = [];
  for(var i = 0; i < m.w; i++) vs[i] = M.col(m, i);
  var us = [];
  for(var i = 0; i < m.w; i++) {
    us[i] = M.gramSchmidtStep_(m, i, us);
  }

  if(opt_normalize) {
    // Gram-Schmidt orthonormalization instead of orthogonalization.
    for(var i = 0; i < us.length; i++) us[i] = M.normalize(us[i]);
  }
  var result = new M(m.h, m.w);
  for(var y = 0; y < m.h; y++) {
    for(var x = 0; x < m.w; x++) {
      result.e[y][x] = us[x].e[y][0];
    }
  }
  return result;
};

// given a partial column vector x of a matrix, returns [v, tau] with:
// v a normalized householder reflection vector
// tau the multiplication factor (0 if degenerate, 2 if real)
// opt_real: set to true if you know the matrix is real
Jmat.Matrix.getHouseholderVector_ = function(x, opt_real) {
  var M = Jmat.Matrix;
  var C = Jmat.Complex;

  // Calculate the householder reflection vector
  var s = x.e[0][0].eqr(0) ? C(-1) : C.sign(x.e[0][0]);
  var v = M.identity(x.h, 1).mulc(s.mul(M.norm(x))).add(x);
  var degenerate = M.isZero(v, 1e-30);
  if(!degenerate) v = v.divc(M.norm(v)); // normalize

  // Calculate the multiplication factor, taking complex matrices and degenerateness into account
  var tau;
  if(degenerate) {
    tau = C(0); // In case of a degenerate column, do no reflection by setting tau to zero
  } else if(opt_real) {
    tau = C(2);
  } else {
    // complex
    var xhv = M.mul(M.row(M.transjugate(x), 0), v);
    var vhx = M.mul(M.transjugate(v), M.col(x, 0));
    tau = xhv.e[0][0].div(vhx.e[0][0]).addr(1);
  }

  return [v, tau];
};

// Returns the c and s parameters for givens rotation as array [c, s].
Jmat.Matrix.getGivensParams_ = function(a, b) {
  if(b.eqr(0)) return [Jmat.Complex(1), Jmat.Complex(0), a];
  // It is very important that the hypot function is implemented as "t=(|x|>|y|)?|y/x|:|x/y|;return sqrt(t*t+1)", and not as the numerically less stable "return sqrt(|x|^2+|y|^2)",
  // or else numerical imprecisions will show up in some eigenvalues of some large matrices.
  var r = Jmat.Complex.hypot(a, b);
  var c = a.div(r);
  var s = b.div(r).neg();
  return [c, s, r];
};

// Pre-multiply for givens transformation G, in-place
// returns G^H * M, calculating only the affected elements
Jmat.Matrix.givensPre_ = function(m, i, j, c, s) {
  for(var x = 0; x < m.w; x++) {
    var a = m.e[i][x];
    var b = m.e[j][x];
    m.e[i][x] = a.mul(c.conj()).add(b.mul(s.neg().conj()));
    m.e[j][x] = a.mul(s).add(b.mul(c));
  }
};

// Post-multiply for givens transformation G, in-place
// returns M * G, calculating only the affected elements
Jmat.Matrix.givensPost_ = function(m, i, j, c, s) {
  for(var y = 0; y < m.h; y++) {
    var a = m.e[y][i];
    var b = m.e[y][j];
    m.e[y][i] = a.mul(c).add(b.mul(s.neg()));
    m.e[y][j] = a.mul(s.conj()).add(b.mul(c.conj()));
  }
};

// do householder reflections to bring a to similar a in upper hessenberg form
Jmat.Matrix.toHessenberg = function(a) {
  var M = Jmat.Matrix;
  var T = M.transjugate;
  if(a.h != a.w) return null;
  var real = Jmat.Matrix.isReal(a);
  var r = M.copy(a);
  for(var k = 0; k + 2 < a.w; k++) {
    var x = M.submatrix(r, k + 1, r.h, k, k + 1); // partial column vector
    var vt = M.getHouseholderVector_(x, real);
    var v = vt[0];
    var tau = vt[1];

    var rs = M.submatrix(r, k + 1, r.h, k, r.w);
    rs = rs.sub(v.mul(T(v)).mul(rs).mulc(tau));
    r = M.insert(r, rs, k + 1, k);

    rs = M.submatrix(r, 0, r.h, k + 1, r.w);
    rs = rs.sub(rs.mul(v).mul(T(v)).mulc(tau));
    r = M.insert(r, rs, 0, k + 1);
  }
  //ensure the elements are really zero
  for(var y = 0; y < r.h; y++) {
    for(var x = 0; x < Math.max(y - 1, 0); x++) {
      r.e[y][x] = Jmat.Complex(0);
    }
  }

  return r;
};

// faster qr algorithm, a must be hessenberg
Jmat.Matrix.qr_hessenberg_ = function(a) {
  var M = Jmat.Matrix;
  var C = Jmat.Complex;
  var n = Math.min(a.h - 1, a.w);

  var r = M.copy(a);
  var q = M.identity(a.h, a.h);
  for(var k = 0; k < n; k++) {
    var g = M.getGivensParams_(r.e[k][k], r.e[k + 1][k]);
    M.givensPre_(r, k, k + 1, g[0], g[1]);
    r.e[k + 1][k] = C(0); // make extra sure it's zero
    M.givensPost_(q, k, k + 1, g[0], g[1]);
  }

  return { q: q, r: r };
};

Jmat.Matrix.qr_general_ = function(a) {
  var M = Jmat.Matrix;
  var T = M.transjugate;
  var real = Jmat.Matrix.isReal(a);
  var h = a.h;
  var w = a.w;
  var v = []; // the householder reflection vectors
  var taus = []; // the multiplication values

  if(a.h < a.w) return null;

  var r = M.copy(a);
  for(var k = 0; k < w; k++) {
    var x = M.submatrix(r, k, r.h, k, k + 1); // partial column vector
    var vt = M.getHouseholderVector_(x, real);
    v[k] = vt[0];
    taus[k] = vt[1];

    // Calculate R: it is A left-multiplied by all the householder matrices
    var rs = M.submatrix(r, k, h, k, w);
    rs = rs.sub(v[k].mul(T(v[k])).mul(rs).mulc(taus[k]));
    r = M.insert(r, rs, k, k);
  }

  // Calculate Q: it is the product of all transjugated householder matrices
  var q = M.identity(h, h);
  for(var k = w - 1; k >= 0; k--) {
    var qs = M.submatrix(q, k, h, 0, h);
    qs = qs.sub(v[k].mul(T(v[k])).mul(qs).mulc(taus[k]));
    q = M.insert(q, qs, k, 0);
  }

  return { q: q, r: r };
};

// QR factorization of complex matrix (with householder transformations)
// requirement: m.h >= m.w
// returns {q: Q, r: R}
// q is h*h unitary matrix
// r is h*w upper triangular matrix
Jmat.Matrix.qr = function(a) {
  /*
  Tests in console:
  var qr = Jmat.qr([[1,2],[3,'4i']]); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  var qr = Jmat.qr([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  var qr = Jmat.qr(Jmat.Matrix(3,3,12,-51,4,6,167,-68,-4,24,-41)); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  degenerate matrix:
  var qr = Jmat.qr([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  */

  // The algorithm for hessenberg form is faster when applicable.
  return Jmat.Matrix.isUpperHessenberg(a, 1e-18) ?
      Jmat.Matrix.qr_hessenberg_(a) :
      Jmat.Matrix.qr_general_(a);
};

// eigenvalues and vectors of 1x1 matrix
Jmat.Matrix.eig11 = function(m) {
  if(m.w != 1 || m.h != 1) return null;
  var result = {};
  result.l = new Jmat.Matrix(1, 1);
  result.l.e[0][0] = Jmat.Complex(m.e[0][0]);
  result.v = new Jmat.Matrix(1, 1);
  result.v.e[0][0] = Jmat.Complex(1);
  return result;
};

// explicit algebraic formula for eigenvalues of 2x2 matrix
Jmat.Matrix.eigval22 = function(m) {
  if(m.w != 2 || m.h != 2) return null;
  var a = Jmat.Complex(1);
  var b = m.e[0][0].neg().sub(m.e[1][1]);
  var c = m.e[0][0].mul(m.e[1][1]).sub(m.e[0][1].mul(m.e[1][0]));
  var d = Jmat.Complex.sqrt(b.mul(b).sub(a.mul(c).mulr(4)));
  var l1 = b.neg().add(d).div(a.mulr(2));
  var l2 = b.neg().sub(d).div(a.mulr(2));
  if(l2.abssq() > l1.abssq()) return [l2, l1];
  return [l1, l2];
};

// explicit algebraic formula for eigenvalues and vectors of 2x2 matrix
// NOTE: the eigenvectors are to be read as columns of v, not rows.
Jmat.Matrix.eig22 = function(m) {
  if(m.w != 2 || m.h != 2) return null;

  var l = Jmat.Matrix.eigval22(m);
  var l1 = l[0];
  var l2 = l[1];

  var v11 = m.e[0][1].div(l1.sub(m.e[0][0]));
  var v12 = Jmat.Complex(1);
  var v21 = m.e[0][1].div(l2.sub(m.e[0][0]));
  var v22 = Jmat.Complex(1);

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

// calculates eigenvalues of complex upper hessenberg matrix, using the shifted QR algorithm with givens rotations
// h is nxn upper hessenberg matrix, and it gets destroyed during the process
Jmat.Matrix.eigval_ = function(h) {
  var C = Jmat.Complex;
  var M = Jmat.Matrix;
  var epsilon = 1.2e-16; // relative machine precision
  var n = h.w; // initially size of the matrix, gets reduced with every eigenvalue found
  var s = C(0); // shift
  var t = C(0); // undo shifts
  var x = C(0), y = C(0), z = C(0);

  var result = [];
  for(var i = 0; i < n; i++) result[i] = C(0);

  while(n > 0) {
    for(var num_it = 0; num_it <= 30; num_it++) {
      if(num_it == 30) return null; // fail after 30 iterations
      // find near-zero sub-diagonal element, at l
      var l;
      for (l = n - 1; l >= 1; l--) {
        if(C.abs1r(h.e[l][l-1]) <= epsilon * (C.abs1r(h.e[l-1][l-1]) + C.abs1r(h.e[l][l]))) {
          h.e[l][l-1] = C(0);  // fix possible numerical imprecisions
          break;
        }
      }
      if(l == n - 1) {
        // root found
        result[n-1] = h.e[n-1][n-1].add(t);
        n--;
        break;
      }
      // calculate shift for shifted QR step
      if(num_it == 10 || num_it == 20) {
        // use a special shift after 10 or 20 iterations
        s.re = Math.abs(h.e[n-1][n-2].re) + Math.abs(h.e[n-2][n-3].re);
        s.im = Math.abs(h.e[n-1][n-2].im) + Math.abs(h.e[n-2][n-3].im);
      } else {
        s = h.e[n-1][n-1];
        x = h.e[n-2][n-1].mul(h.e[n-1][n-2]);
        if(!x.eqr(0)) {
          y = (h.e[n-2][n-2].sub(s)).divr(2);
          z = C.sqrt(y.mul(y).add(x));
          if(y.re * z.re + y.im * z.im < 0) z = z.neg();
          x = x.div(y.add(z));
          s = s.sub(x);
        }
      }
      // apply shift
      for(var i = 0; i < n; i++) {
        h.e[i][i] = h.e[i][i].sub(s);
      }
      t = t.add(s);

      // fast QR step with givens rotations. Implicitely decomposes h into Q*R, then sets h to R*Q (also, only uses nxn sub-matrix of original h).
      var g = []; // remember each of the givens parameters
      for (var k = 0; k + 1 < n; k++) {
        g[k] = M.getGivensParams_(h.e[k][k], h.e[k + 1][k]);
        M.givensPre_(h, k, k + 1, g[k][0], g[k][1]);
        h.e[k + 1][k] = C(0); // make extra sure it's zero
        h.e[k][k] = g[k][2]; // also, the r returned by getGivensParams_ is more precise here
      }
      for (var k = 0; k + 1 < n; k++) {
        M.givensPost_(h, k, k + 1, g[k][0], g[k][1]);
      }
    }
  }
  return result;
};

// Returns the eigenvalues (aka spectrum) of m in an array, from largest to smallest
// The eigenvalues also give the characteristic polynomial: (x-l[0])*(x-l[1])*...*(x-l[n-1])
Jmat.Matrix.eigval = function(m) {
  var M = Jmat.Matrix;
  if(m.w != m.h || m.w < 1) return null;
  var n = m.w;
  if(n == 1) return [m.e[0][0]];
  if(n == 2) return M.eigval22(m);

  // TODO: for hermitian or symmetric matrix, use faster algorithm for eigenvalues

  var a = M.toHessenberg(m);
  var l = M.eigval_(a);

  // Fullfill our promise of eigenvalues sorted from largest to smallest magnitude, the eigenvalue algorithm usually has them somewhat but not fully correctly sorted
  l.sort(function(a, b) { return b.abssq() - a.abssq(); });

  return l;
};

// Returns the eigenvector matching the given eigenvalue of m as a column vector
// m must be a square matrix, and lambda a correct eigenvalue of it
// opt_normalize is how to normalize the eigenvector: 0: don't (length unspecified), 1: last element equals 1, 2: length 1. The default is "1".
Jmat.Matrix.eigenVectorFor = function(m, lambda, opt_normalize) {
  var M = Jmat.Matrix;
  if(m.w != m.h || m.w < 1) return null;
  var n = m.w;
  var normalize_mode = (opt_normalize == undefined) ? 1 : opt_normalize;

  // Find eigenvectors by solving system of linear equations.
  m = M.copy(m); //avoid changing user input
  for(var i = 0; i < n; i++) m.e[i][i] = m.e[i][i].sub(lambda);
  var f = M.zero(n, 1);
  var g = M.solve(m, f, 0.01); // a very large epsilon is used... because the eigenvalues are numerically not precise enough to give the singular matrix they should
  if(g) {
    if(normalize_mode == 2) g = g.divc(M.norm(g));
    if(normalize_mode == 1) if(!g.e[n - 1][0].eqr(0)) g = g.divc(g.e[n - 1][0]);
  } else {
    // failed to find the corresponding eigenvectors, avoid crash, set values to NaN
    g = M.zero(n, 1);
    for(var i = 0; i < n; i++) g.e[i][0] = Jmat.Complex(NaN, NaN); // The eigenvectors are stored as column vectors
  }
  return g;
};

// Returns eigenvalues and eigenvectors of real symmetric matrix with the Jacobi eigenvalue algorithm, as { l: eigenvalues, v: eigenvectors }
// For correct result, requires that m is square, real and symmetric.
Jmat.Matrix.jacobi_ = function(m, opt_epsilon, opt_normalize) {
  var a = [];
  var v = [];
  var n = m.w;
  for(var y = 0; y < n; y++) {
    a[y] = [];
    for(var x = 0; x < n; x++) {
      a[y][x] = m.e[y][x].re;
    }
  }
  Jmat.Real.matrix_jacobi(a, v, n, opt_epsilon);

  if(opt_normalize == 1 || opt_normalize == undefined) {
    for(var y = 0; y < n; y++) {
      if(v[y][n - 1] == 0) continue;
      for(var x = 0; x < n; x++) {
        v[y][x] /= v[y][n - 1];
      }
    }
  }
  v = Jmat.Matrix(v).transpose();
  var l = Jmat.Matrix(n, 1);
  for(var i = 0; i < n; i++) {
    l.e[i][0] = Jmat.Complex(a[i][i]);
  }

  return { l: l, v: v };
};

// Returns the eigenvectors and eigenvalues of m as { l: eigenvalues, v: eigenvectors }
// eigenvalues as n*1 column vector, eigenvectors as n*n matrix
// for each column of v and corresponding eigenvalue: A*v = l*v (l represents lambda, A is m)
// opt_normalize is how to normalize the eigenvectors: 0: don't (length unspecified), 1: last element equals 1, 2: length 1. The default is "1".
// NOTE: the eigenvectors are to be read as columns of v, not rows.
Jmat.Matrix.eig = function(m, opt_normalize) {
  var M = Jmat.Matrix;
  if(m.w != m.h || m.w < 1) return null;
  var n = m.w;
  if(n == 1) return M.eig11(m);
  if(n == 2) return M.eig22(m);

  if(M.isReal(m) && M.isSymmetric(m)) return Jmat.Matrix.jacobi_(m, opt_normalize);

  var l = M.eigval(m);
  for(var i = 0; i < l.length; i++) if(Jmat.Complex.nearr(l[i], 0, 1e-15)) l[i] = Jmat.Complex(0); // this avoids numerical instability problems with calculation of eigenvectors in case of eigenvalues that should be zero, but are close to it instead

  // Fullfill our promise of eigenvalues sorted from largest to smallest magnitude, the eigenvalue algorithm usually has them somewhat but not fully correctly sorted
  l.sort(function(a, b) { return b.abssq() - a.abssq(); });

  var v = null;
  // TODO: use more efficient algorithm for eigenvectors, e.g. something that produces them as side-effect of the eigenvalue calculation
  // TODO: for hermitian matrix, use faster algorithm for eigenvectors (like is already done with jacobi for real symmetric)
  // TODO: solve numerical imprecisions, e.g. see how high epsilon eigenVectorFor is using to circumvent various problems in unstable ways
  var v = new M(n, n);
  for(var j = 0; j < n; j++) {
    var g = M.eigenVectorFor(m, l[j], opt_normalize);
    M.setCol(v, g, j);
  }

  return { l: l, v: v };
};

// Returns the eigen decomposition (aka spectral decomposition) of m as { v: V, d: D }
// If M is diagonizable, M = V * D * V^(-1)
// In other words: m == result.v.mul(result.d).mul(Jmat.Matrix.inv(result.v))
// This function is very similar to Jmat.Matrix.eig. v is the same, d is the same as l but put on the diagonal of a matrix
Jmat.Matrix.evd = function(m) {
  var eig = Jmat.Matrix.eig(m);

  return { v: eig.v, d: Jmat.Matrix.diag(eig.l) };
};

// TODO: Matrix.jordan, calculating jordan decomposition, jordan normal form, jordan canonical form, giving {v:V, j:J} with the generalized eigenvectors in V. (similar to Matrix.evd, but supports any matrix rather than only diagonizable ones)

// returns the definiteness as an array of 3 booleans
// the first boolean means negative eigenvalues are present
// the second boolean means eigenvalues of zero are present
// the third boolean means positive eigenvalues are present
// ignores whether it's hermitian or not, but requires square matrix
Jmat.Matrix.definiteness_ = function(m, opt_epsilon) {
  var epsilon = (opt_epsilon == undefined) ? 1e-12 : opt_epsilon; // default is 1e-12 instead of the usual 1e-15 because eigenvalues of [[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]] is otherwise not precise enough, making it think that one is indefinite
  // TODO: faster method than eigenvalues
  var e = Jmat.Matrix.eigval(m);
  if(!e) return null;
  var result = [false, false, false];
  for (var i = 0; i < e.length; i++) {
    if(Jmat.Real.near(e[i].re, 0, epsilon)) result[1] = true;
    else if(e[i].re > 0) result[2] = true;
    else result[0] = true;
  }
  return result;
};

Jmat.Matrix.INDEFINITE = 0;
Jmat.Matrix.POSITIVE_DEFINITE = 1;
Jmat.Matrix.POSITIVE_SEMI_DEFINITE = 2;
Jmat.Matrix.NEGATIVE_DEFINITE = 3;
Jmat.Matrix.NEGATIVE_SEMI_DEFINITE = 4;

// For a hermitian matrix (symmetric if real), returns its definiteness using the constants above.
// For non-hermitian matrix, returns null to indicate invalid. Use (M + transjugate(M)) / 2 to get it anyway
// Returns only one value that best describes the matrix. In practice, multiple values may apply:
// -a POSITIVE_DEFINITE matrix is always also POSITIVE_SEMI_DEFINITE
// -a NEGATIVE_DEFINITE matrix is always also NEGATIVE_SEMI_DEFINITE
// -a null matrix is both POSITIVE_SEMI_DEFINITE and NEGATIVE_SEMI_DEFINITE but this function only returns one of those
Jmat.Matrix.definiteness = function(m, opt_epsilon) {
  if(!Jmat.Matrix.isHermitian(m, opt_epsilon)) return null;
  var bools = Jmat.Matrix.definiteness_(m, opt_epsilon);
  if(bools[2] && !bools[0] && !bools[1]) return Jmat.Matrix.POSITIVE_DEFINITE;
  if(bools[0] && !bools[2] && !bools[1]) return Jmat.Matrix.NEGATIVE_DEFINITE;
  if(!bools[0]) return Jmat.Matrix.POSITIVE_SEMI_DEFINITE;
  if(!bools[2]) return Jmat.Matrix.NEGATIVE_SEMI_DEFINITE;
  if(bools[0] && bools[2]) return Jmat.Matrix.INDEFINITE;
  return null; // unreachable
};

// converts an array of complex values or regular JS numbers, to a column vector matrix
Jmat.Matrix.arrayToCol = function(a) {
  var r = new Jmat.Matrix(a.length, 1);
  for(var i = 0; i < a.length; i++) r.e[i][0] = Jmat.Complex.cast(a[i]);
  return r;
};

// converts an array of complex values or regular JS numbers, to a row vector matrix
Jmat.Matrix.arrayToRow = function(a) {
  var r = new Jmat.Matrix(1, a.length);
  for(var i = 0; i < a.length; i++) r.e[0][i] = Jmat.Complex.cast(a[i]);
  return r;
};

// converts an array of complex values or regular JS numbers, to a diagonal matrix
Jmat.Matrix.arrayToDiag = function(a) {
  //var r = Jmat.Matrix.zero(a.length, a.length);
  var r = new Jmat.Matrix(a.length, a.length);
  for(var i = 0; i < a.length; i++) r.e[i][i] = Jmat.Complex.cast(a[i]);
  return r;
};

// Puts all the elements of matrix d in a single diagonal matrix
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

// Puts all diagonal elements from a into a single column vector
Jmat.Matrix.diagToCol = function(a) {
  var n = Math.min(a.h, a.w);
  var result = Jmat.Matrix(n, 1);
  for(var i = 0; i < n; i++) result.e[i][0] = a.e[i][i];
  return result;
};

// Puts all diagonal elements from a into a single row vector
Jmat.Matrix.diagToRow = function(a) {
  var n = Math.min(a.h, a.w);
  var result = Jmat.Matrix(1, n);
  for(var i = 0; i < n; i++) result.e[0][i] = a.e[i][i];
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

// Dot product of two vectors (aka inner product, scalar product, projection product).
// Also supports it for matrices of same dimensions (it then is the Frobenius inner product)
// If vectors, row and column vectors may be mixed.
Jmat.Matrix.dot = function(a, b) {
  if(a.w != b.w || a.h != b.h) {
    if(!(a.w == b.h && a.h == b.w && (a.w == 1 || a.h == 1))) return Jmat.Matrix(NaN); // Do allow it for differently orientated vectors (h or w is 1)
  }
  var n = a.w * a.h;
  var result = Jmat.Complex(0);
  for(var i = 0; i < n; i++) result = result.add(a.get1(i).mul(b.get1(i).conj()));
  return result;
};


/** @license
License of Jmat.Matrix.zsvdc_: this function is from linpack, from http://www.netlib.org/linpack/
The license is not mentioned directly in the source code or the website, but
has been said to now be a variant of the BSD license: see
https://bugzilla.redhat.com/show_bug.cgi?id=1000829.
Here is the original author comment from zsvdc.f:
    linpack. this version dated 03/19/79 .
             correction to shift calculation made 2/85.
    g.w. stewart, university of maryland, argonne national lab.
*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Based on LINPACK zsvdc.f, converted to JavaScript.
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
  var i,iter,j,jobu,k,kase,kk,l,ll,lls,lp1,ls,lu,m,
      mm,mm1,mp1,nct,ncu,nrt,info; // integers
  var maxit = 30;
  var t,r; //complex
  var b,c,cs,el,emm1,f,g,dznrm2,scale,shift,sl,sm,sn,
      smm1,t1,test,ztest; // double

  var dr; //drotg result

  var dreal = function(z) { return z.re; };
  var cabs1 = function(z) { return Math.abs(z.re) + Math.abs(z.im); };
  var nearzero = function(z) {
    // At each nearzero call below, the original zsvdc instead checked "cabs1(z) != 0.0". Usually a division then follows. That
    // complex division can still result in inf and NaN for tiny numbers (e.g. 1e-307). Hence replaced with this check here.
    return cabs1(z) < 1e-150;
  };
  // returns value with absolute value of x, argument of y (transfers sign)
  var csign = function(x, y) { return y.eqr(0) ? Jmat.Complex(0) : y.mulr(x.abs() / y.abs()); };
  var sign = function(x, y) { return y == 0 ? 0 : (y < 0 ? -Math.abs(x) : Math.abs(x)); };
  // Euclidean norm of complex vector, n elements starting at index start
  var dznrm2 = function(n, arr, start) {
    var result = Jmat.Complex(0);
    for(var i = 0; i < n; i++) {
      var e = arr[start + i];
      result = result.add(e.mul(e.conj()));
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
      result = result.add(arrx[startx + i].conj().mul(arry[starty + i]));
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
      var ay = arry[starty + i];
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
      if(!nearzero(s[l])) {
        if(!nearzero(x[l + l * ldx])) s[l] = csign(s[l], x[l + l * ldx]);
        t = Jmat.Complex(1.0).div(s[l]);
        zscal(n - l, t, x, l + l * ldx);
        x[l + l * ldx] = x[l + l * ldx].addr(1);
      }
      s[l] = s[l].neg();
    }
    for(j = lp1; j < p; j++) {
      if(l < nct) {
        if(!nearzero(s[l])) {
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

      if(!nearzero(e[l])) {
        if(!nearzero(e[lp1])) {
          e[l] = csign(e[l], e[lp1]);
        }
        t = Jmat.Complex(1.0).div(e[l]);
        zscal(p - l - 1, t, e, lp1);
        e[lp1] = Jmat.Complex(1.0).add(e[lp1]);
      }
      e[l] = e[l].conj().neg();
      // apply the transformation.
      if(lp1 < n && !nearzero(e[l])) {
        for(j = lp1; j < n; j++) {
          work[j] = Jmat.Complex(0.0);
        }
        for(j = lp1; j < p; j++) {
          zaxpy(n - l - 1, e[j], x, lp1 + j * ldx, work, lp1);
        }
        for(j = lp1; j < p; j++) {
          zaxpy(n - l - 1, (e[j].neg().div(e[lp1])).conj(), work, lp1, x, lp1 + j * ldx);
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
        if(!nearzero(e[l])) {
          for(j = lp1; j < p; j++) {
            t = zdotc(p - lp1, v, lp1 + l * ldv, v, lp1 + j * ldv).neg().div(v[lp1 + l * ldv]);
            zaxpy(p - lp1, t, v, lp1 + l * ldv, v, lp1 + j * ldv);
          }
        }
      }
      for(i = 0; i < p; i++) {
        v[i + l * ldv] = Jmat.Complex(0);
      }
      v[l + l * ldv] = Jmat.Complex(1);
    }
  }
  // transform s and e so that they are real.
  for(i = 0; i < m; i++) {
    if(!nearzero(s[i])) {
      t = Jmat.Complex.abs(s[i]);
      r = s[i].div(t);
      s[i] = t;
      if(i + 1 < m) e[i] = e[i].div(r);
      if(wantu) zscal(n, r, u, i * ldu);
    }
    if(i + 1 == m) break;
    if(!nearzero(e[i])) {
      t = Jmat.Complex.abs(e[i]);
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
      test = s[l - 1].abs() + s[l].abs();
      ztest = test + e[l - 1].abs();
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
        if(ls != m) test = test + e[ls - 1].abs();
        if(ls != l + 1) test = test + e[ls - 2].abs();
        ztest = test + s[ls - 1].abs();
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
      for(kk = l; kk <= mm1; kk++) {
        k = mm1 - kk + l;
        t1 = dreal(s[k - 1]);
        dr = drotg(t1, f); t1 = dr[0]; f = dr[1]; cs = dr[2]; sn = dr[3];
        s[k - 1] = Jmat.Complex(t1);
        if(k != l) {
          f = -sn * dreal(e[k - 2]);
          e[k - 2] = e[k - 2].mulr(cs);
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
      scale = Math.max(Math.max(Math.max(Math.max(s[m - 1].abs(),
              s[m - 2].abs()), e[m - 2].abs()), s[l - 1].abs()),
              e[l - 1].abs());
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
// input M, returns {u: U, s: S, v: V } such that U * S * V^T = M and S diagonal with singular values (^T means conjugate transpose here)
// Input allowed to be non-square. The size of "S" is same as the input matrix.
Jmat.Matrix.svd = function(m) {
  /*
  Checks in console:
  function testSvd(m) {
    var result = Jmat.Matrix.svd(Jmat.Matrix(m));
    console.log(Jmat.toString(result) + ' | ' + Jmat.Matrix.mul(Jmat.Matrix.mul(result.u, result.s), Jmat.Matrix.transjugate(result.v)).toString());
  }
  testSvd(Jmat.Matrix(2,2,1,2,3,4))
  testSvd(Jmat.Matrix(2,2,1,2,3,4).mulc(Jmat.Complex.I))
  testSvd(Jmat.Matrix(2,2,1,2,1,2))
  testSvd(Jmat.Matrix([[1,2]]))
  testSvd(Jmat.Matrix([[1],[2]]))

  var x = Jmat.Matrix([[1,2,3,4],[5,6,7,8]]); var svd = Jmat.svd(x); 'X:\n' + Jmat.Matrix.render(x) + '\nU:\n' + Jmat.Matrix.render(svd.u) + '\nS:\n' + Jmat.Matrix.render(svd.s) + '\nV:\n' + Jmat.Matrix.render(svd.v) + '\n reconstruct: \n' + Jmat.Matrix.render(svd.u.mul(svd.s).mul(svd.v.transjugate()))
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
      if(!a.e[y][x].eq(b.e[y][x])) return false;
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
      if(!Jmat.Complex.near(ea, eb, epsilon)) return false;
    }
  }

  return true;
};

// nearly equal (relative, precision determines number of decimal digits to match per number, e.g. 1e-3 for 3 digits)
Jmat.Matrix.relnear = function(a, b, precision) {
  if(a.w != b.w || a.h != b.h) return false;

  for(var y = 0; y < a.h; y++) {
    for(var x = 0; x < a.w; x++) {
      var ea = a.e[y][x];
      var eb = b.e[y][x];
      if(!Jmat.Complex.relnear(ea, eb, precision)) return false;
    }
  }

  return true;
};


//solves system of linear equations ax=b.
//Returns null if the system is inconsistent and has no solution, x otherwise.
//If multiple solutions are possible, returns the solution where the vector of free variables is 0.
//a: input matrix, h*w size (h = number of equations, w = number of unknowns)
//b: input vector, h*1 size
//result: w*1 size
//opt_epsilon: optional parameter, in case of homogenous system, how near the smallest singular value of a should be to return a solution
Jmat.Matrix.solve = function(a, b, opt_epsilon) {
  var M = Jmat.Matrix;
  if(a.h != b.h) return undefined; // input error
  var epsilon = (opt_epsilon == undefined) ? 1e-12 : opt_epsilon;

  if (M.isZero(b)) {
    // Homogenous system.
    var svd = M.svd(a);
    if (!Jmat.Complex.nearr(svd.s.e[svd.s.h - 1][svd.s.h - 1], 0, epsilon)) return null; //allow large imprecision, so that it works if eigenvector calculation is a bit imprecise.
    // A is assumed singular. Last column should correspond to singular value which is zero, corresponding to a solution.
    return M.col(svd.v, a.w - 1);
  }

  // TODO: more numerically stable algorithm?
  var ab = M.augment(a, b);
  var r = M.rref(ab);
  var result = M(a.w, 1, 0);
  var x0 = 0;
  for(var y = 0; y < r.h; y++) {
    var x = x0;
    for(; x < r.w; x++) {
      if(!r.e[y][x].eqr(0)) {
        if(x == r.w - 1) return null; // inconsistent system: row with all zeroes except in the augmented part
        result.e[x][0] = r.e[y][ab.w - 1];
        x0 = x + 1;
        break;
      }
    }
    if(x == r.w - 1) break; // done
  }
  return result;
};

// Gives the least squares solution to the (overdetermined) linear system ax=b, that is, minimizes ||ax-b||^2
Jmat.Matrix.leastSquares = function(a, b) {
  return Jmat.Matrix.pseudoinverse(a).mul(b);
};

// Returns the matrix in reduced row echolon form. Supports non-square and singular matrices. Is always the identity matrix if the input is invertible.
Jmat.Matrix.rref = function(a) {
  /*
  Check in console:
  Matrix.rref(Matrix([[1,3,1,9],[1,1,-1,1],[3,11,5,35]])).render()
  Matrix.solve(Matrix([[1,3,1],[1,1,-1],[3,11,5]]), Matrix([[9],[1],[35]])).render()
  */
  var C = Jmat.Complex;
  var h = a.h;
  var w = a.w;
  a = Jmat.Matrix.copy(a); // The rest all works in-place, so copy to not modify user input.

  //swaps rows y0 and y1 in place
  var swaprow = function(a, y0, y1) {
    for (var i = 0; i < w; i++) {
      var temp = a.e[y0][i]; a.e[y0][i] = a.e[y1][i]; a.e[y1][i] = temp;
    }
  };

  // only starts at x rather than from the beginning
  var mulrow = function(a, x, y, v) {
    for (var i = x; i < w; i++) {
      a.e[y][i] = a.e[y][i].mul(v);
    }
  };

  // subtracts a multiple of row y0 from row y1, modifying y1. Only starts at x rather than from the beginning
  var submul = function(a, x, y0, v, y1) {
    var w = a.e[0].length;
    for (var i = x; i < w; i++) {
      a.e[y1][i] = a.e[y1][i].sub(a.e[y0][i].mul(v));
    }
  };

  var pivots = []; // x coordinate of pivot in each row (except the zero rows at the end, so may have smaller length than h)

  // gaussian elimination
  var k2 = 0; //next row to fill in, equal to k unless there are zero-rows
  for(var k = 0; k < w; k++) {
    var n = Jmat.Real.argmax(k2, h, function(i) { return a.e[i][k].abssq(); });
    if (a.e[n][k].eqr(0)) continue; // singular, no pivot for this column
    if(n != k2) swaprow(a, k2, n);
    mulrow(a, k, k2, a.e[k2][k].inv()); // pivot is now 1
    // make corresponding elements of row below zero using row operations
    for (var i = k2 + 1; i < h; i++) {
      if(!(a.e[i][k].eqr(0))) {
        submul(a, k + 1, k2, a.e[i][k], i);
        a.e[i][k] = C(0); // make extra-sure it's 0, avoid numerical imprecision
      }
    }
    pivots.push(k);
    k2++;
    if(k2 >= h) break;
  }

  //now bring from row echolon form to reduced row echolon form
  for(var k = 0; k < pivots.length; k++) {
    var p = pivots[k];
    // make corresponding elements of row above zero using row operations
    for(var y = k - 1; y >= 0; y--) {
      if(!(a.e[y][p].eqr(0))) {
        submul(a, p + 1, k, a.e[y][p], y);
        a.e[y][p] = C(0); // make extra-sure it's 0, avoid numerical imprecision
      }
    }
  }

  return a;
};

// generates a random matrix with some properties.
// all parameters are optional.
// properties: properties object similar to what "getProperties" returns. Not all are implemented though.
//             Currently only 'real', 'integer', 'binary' and 'hermitian' are supported.
//             Default is object {'real':true}. Pass "undefined" to get the default, or '{}' to get complex matrices.
// h: number of rows. Default: random 2..12
// w: number of columns. Default: random 2..12, or h if h is defined
// r0: lowest random value. Default: 0
// r1: highest random value. Default: 1
// s: sparseness: give value in range 0-1. Default: 1
// e.g. Matrix.random(4, 4, 0, 10, 1, {'integer':true}).render_summary()
//      Matrix.random(3, 3, -10, 10, 1, {'real':false,'hermitian':true}).render_summary()
Jmat.Matrix.random = function(properties, h, w, r0, r1, s) {
  var C = Jmat.Complex;
  w = w || h || Math.floor(Math.random() * 10 + 2);
  h = h || Math.floor(Math.random() * 10 + 2);
  s = (s == undefined) ? 1 : s;
  r0 = (r0 == undefined) ?  0 : r0;
  r1 = (r1 == undefined) ?  1 : r1;
  properties = properties || {'real':true};
  var real = properties['real'] || properties['integer'];
  var integer = properties['integer'];
  var binary = properties['binary'];
  var hermitian = properties['hermitian'];
  var f = function(real) {
    if(s >= 1 || Math.random() < s) {
      if (binary) return Math.random() > 0.5 ? C(1) : C(0);
      var result = real ? C(Math.random() * (r1 - r0) + r0) : C.random(r0, r1);
      return integer ? C(Math.floor(result.re)) : result;
    }
    return C(0);
  };

  var result = Jmat.Matrix(h, w);
  for(var y = 0; y < h; y++) {
    for(var x = 0; x < w; x++) {
      if(hermitian && y > x) result.e[y][x] = result.e[x][y].conj();
      else if(hermitian && x == y) result.e[y][x] = f(true);
      else result.e[y][x] = f(real);
    }
  }

  return result;
};

Jmat.Matrix.matrixfft_ = function(m, inverse) {
  var rowresult = new Jmat.Matrix(m.h, m.w);

  // apply to each row
  if(m.w > 1) {
    for(var j = 0; j < m.h; j++) {
      var fft = Jmat.Complex.fft(m.e[j], inverse);
      for(var i = 0; i < m.w; i++) rowresult.e[j][i] = fft[i];
    }
  } else {
    rowresult = m;
  }

  var result = new Jmat.Matrix(m.h, m.w);

  // apply to each column
  if (m.h > 1) {
    for(var j = 0; j < m.w; j++) {
      var col = Jmat.Matrix.transpose(Jmat.Matrix.col(rowresult, j));
      var fft = Jmat.Complex.fft(col.e[0], inverse);
      for(var i = 0; i < m.h; i++) result.e[i][j] = fft[i];
    }
  } else {
    result = rowresult;
  }

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
  var e = Jmat.Matrix.evd(m);
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
  var e = Jmat.Matrix.evd(m);
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
  var e = Jmat.Matrix.evd(m);
  var v = e.v;
  var d = e.d;
  for(var i = 0; i < d.w; i++) d.e[i][i] = d.e[i][i].pow(s);
  return v.mul(d).mul(Jmat.Matrix.inv(v));
};

Jmat.Matrix.convolve = function(a, b, opt_grow) {
  var C = Jmat.Complex;
  // shift of b. TODO: optional configurable
  var sy = Math.floor(b.h / 2);
  var sx = Math.floor(b.w / 2);
  if (opt_grow) {
    var h2 = a.h + b.h - 1;
    var w2 = a.w + b.w - 1;
    var result = Jmat.Matrix(h2, w2);
    for (var y = 0; y < h2; y++) {
      for (var x = 0; x < w2; x++) {
        var r = C(0);
        for (var y2 = 0; y2 < b.h; y2++) {
          for (var x2 = 0; x2 < b.w; x2++) {
            var y3 = y - sy * 2 + y2;
            var x3 = x - sx * 2 + x2;
            if(x3 < 0 || x3 >= a.w || y3 < 0 || y3 >= a.h) continue;
            r = r.add(a.e[y3][x3].mul(b.e[y2][x2]));
          }
        }
        result.e[y][x] = r;
      }
    }
    return result;
  } else {
    var result = Jmat.Matrix(a.h, a.w);
    for (var y = 0; y < a.h; y++) {
      for (var x = 0; x < a.w; x++) {
        var r = C(0);
        for (var y2 = 0; y2 < b.h; y2++) {
          for (var x2 = 0; x2 < b.w; x2++) {
            var y3 = y - sy + y2;
            var x3 = x - sx + x2;
            if(x3 < 0 || x3 >= a.w || y3 < 0 || y3 >= a.h) continue;
            r = r.add(a.e[y3][x3].mul(b.e[y2][x2]));
          }
        }
        result.e[y][x] = r;
      }
    }
    return result;
  }
};

// Finds roots of a polynomial given its coefficients sorted from the one belonging to x^0 to the one belonging to x^(n-1)
// Returns n - 1 complex roots
// This function is in jmat_matrix.js because it needs to calculate eigenvalues to do this. However, we add the function to Jmat.Complex because it really belongs there from user perspective.
Jmat.Complex.polyroots = function(coeffs) {
  var C = Jmat.Complex;
  var M = Jmat.Matrix;
  if(coeffs.length <= 1) return []; // no roots for unexisting or horizontal function. We do not take into account the infinite roots for the function f(x)=0.
  var v = coeffs[coeffs.length - 1].inv();
  var m = coeffs.length - 1;
  var matrix = M(m, m, 0); // companion matrix
  for (var i = 1; i < m; i++) matrix.e[i][i - 1] = C(1);
  for (var i = 0; i < m; i++) matrix.e[i][m - 1] = v.mul(coeffs[i]).neg();

  return M.eigval(matrix);
};


