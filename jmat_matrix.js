/*
Jmat.js

Copyright (c) 2011-2015, Lode Vandevenne
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

Overview of some functionality:
-decompositions: Matrix.lu, Matrix.qr, Matrix.svd, Matrix.evd
-inverse: Matrix.inv, Matrix.pseudoinverse
-solve: Matrix.solve
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
-pretty print: Matrix.render
*/

/*
Constructor
height first because that's the math convention: a 2x3 matrix is 2 rows high, 3 columns wide, and made as new Jmat.Matrix(2, 3).
Does not initialize elements if keyword "new" is in front. If keyword "new" is not in front, then uses Jmat.Matrix.make and all its options to initialize elements instead.

Aliased as simply "Matrix" by jmat.js - disable that if it causes name clashes

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

//similar to toString, but using curly braces instead of square brackets
Jmat.Matrix.toCurly = function(m) {
  return Jmat.Matrix.toString(m).replace(/\[/g, '{').replace(/\]/g, '}');
};
Jmat.Matrix.prototype.toCurly = function() {
  return Jmat.Matrix.toCurly(this);
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
  opt_precision = opt_precision == undefined ? 3 : opt_precision;
  //turn undefined into nan
  if (!Jmat.Matrix.isValid(a)) {
    a = Jmat.Matrix.copy(a);
    for (var y = 0; y < a.h; y++) {
      for (var x = 0; x < a.w; x++) {
        if(!a.e[y][x]) a.e[y][x] = Jmat.Complex(NaN);
      }
    }
  }
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
Jmat.Matrix.isSymmetrical = function(a, opt_epsilon) {
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
Jmat.Matrix.isSkewSymmetrical = function(a, opt_epsilon) {
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

// Returns an object with various named boolean and scalar properties of the given matrix
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
  result['NaN'] = M.isNaN(a);

  // The following properties only make sense for square matrices
  result['identity'] = M.isIdentity(a);
  result['diagonal'] = M.isDiagonal(a);
  result['tridiagonal'] = M.isTridiagonal(a);
  result['symmetrical'] = M.isSymmetrical(a);
  result['hermitian'] = M.isHermitian(a);
  result['skewHermitian'] = M.isSkewHermitian(a);
  result['skewSymmetrical'] = M.isSkewSymmetrical(a);
  result['upperTriangular'] = M.isUpperTriangular(a);
  result['lowerTriangular'] = M.isLowerTriangular(a);
  result['strictlyUpperTriangular'] = M.isStrictlyUpperTriangular(a);
  result['strictlyLowerTriangular'] = M.isStrictlyLowerTriangular(a);
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
// Does not show redundant properties. E.g. if the matrix is 'identity', will not show 'symmetrical', if it's 'normal', will not show 'orthogonal', etc...
// To see every single property instead, do "Jmat.toString(Jmat.Matrix.getProperties(a))"
Jmat.Matrix.summary = function(a) {
  var p = Jmat.Matrix.getProperties(a);

  var toName = function(name) {
    // convert camelCase to lower case with spaces
    if(name != 'NaN') name = name.replace(/([A-Z])/g, ' $1').toLowerCase();
    // But keep own names
    name = name.replace('hessenberg', 'Hessenberg');
    name = name.replace('frobenius', 'Frobenius');
    name = name.replace('toeplitz', 'Toeplitz');
    name = name.replace('hankel', 'Hankel');
    name = name.replace('frobenius', 'Frobenius');
    return name;
  };

  //order of non-square related properties
  var nonsquare = ['height', 'width', 'zero', 'real', 'NaN',
                   'rank', 'frobeniusNorm', 'spectralNorm', 'conditionNumber'];
  //order of properties only applicable for square matrices
  var square = ['identity', 'symmetrical', 'hermitian', 'skewSymmetrical', 'skewHermitian', 'diagonal', 'tridiagonal',
                'upperTriangular', 'lowerTriangular', 'strictlyUpperTriangular', 'strictlyLowerTriangular', 'upperHessenberg', 'lowerHessenberg',
                'singular', 'invertible', 'determinant', 'trace', 'orthogonal', 'unitary', 'normal', 'permutation', 'toeplitz', 'hankel',
                'indefinite', 'positiveDefinite', 'negativeDefinite', 'positiveSemidefinite', 'negativeSemidefinite', 'frobenius'];

  var opposite = { 'square' : 'non-square', 'real' : 'complex' };
  // these properties are added only to avoid some redundancy in summary output with the "sub" sytem
  p['small2x2'] = (a.w <= 2 && a.h <= 2);
  p['small1x1'] = (a.w <= 1 && a.h <= 1);
  p['realsym'] = p['real'] && p['symmetrical'];
  p['realskewsym'] = p['real'] && p['skewSymmetrical'];
  // pairs of child:parents, where child is always true if any of the parents is true, with the intention to not display child in a list if parent is already true as it's redundant
  var sub = {
    'strictlyUpperTriangular': ['zero'], 'strictlyLowerTriangular' : ['zero'],
    'upperTriangular' : ['diagonal', 'strictlyUpperTriangular'], 'lowerTriangular' : ['diagonal', 'frobenius', 'strictlyLowerTriangular'],
    'upperHessenberg' : ['upperTriangular', 'tridiagonal'], 'lowerHessenberg' : ['lowerTriangular', 'tridiagonal'],
    'diagonal' : ['small1x1', 'identity', 'zero'], 'tridiagonal' : ['small2x2', 'diagonal'],
    'orthogonal' : ['normal', 'identity'], 'unitary' : ['normal'], 'normal' : ['identity', 'zero'],
    'hermitian' : ['normal'], 'hermitian' : ['realsym'],  'skewHermitian' : ['realskewsym'],
    'symmetrical' : ['diagonal'], 'skewSymmetrical' : ['zero'],
    'permutation' : ['identity'], 'invertible' : ['identity'], 'singular' : ['zero'],
    'real' : ['identity', 'zero'], 'toeplitz' : ['identity', 'zero'], 'hankel' : ['zero'], 'frobenius' : ['identity'],
    'positiveDefinite' : ['identity'], 'negativeSemidefinite' : ['zero', 'negativeDefinite'], 'positiveSemidefinite' : ['zero', 'positiveDefinite']
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
  var det = p['square'] ? (', determinant ' + p['determinant']) : '';
  summary = '' + summary + ' matrix with rank ' + p['rank'] + det + ' and condition number ' + p['conditionNumber'] + '.\n';

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

// Internal algorithm for lu.
// Returns L and U merged into one matrix (without L's diagonal 1's), and a pivot array with element i the pivot row interchanged with row i.
Jmat.Matrix.doolittle_lup_ = function(a) {
  if(a.h != a.w) return null; //must be square
  var M = Jmat.Matrix;
  var C = Jmat.Complex;
  a = M.copy(a); // we'll modify it to interchange rows and insert elements from L and U

  var pivot = [];

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
    }

    //Returning for singular commented out: still works, resulting U will be singular.
    //if(C.nearr(a.e[k][k], 0, 1e-15)) return null; // singular

    for(var i = k + 1; i < a.h; i++) {
      a.e[i][k] = a.e[i][k].div(a.e[k][k]);
      if(C.isNaN(a.e[i][k])) a.e[i][k] = C(0); // Set 0/0 to 0 for singular input matrix.
    }
    for(var i = k + 1; i < a.h; i++) {
      for(var j = k + 1; j < a.h; j++) {
        a.e[i][j] = a.e[i][j].sub(a.e[i][k].mul(a.e[k][j]));
      }
    }
  }

  return [a, pivot];
};

// LUP decomposition. Returns object {p: P, l: L, u: U} such that A = P*L*U, with P a permutation matrix, L lower triangular, U upper triangular
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
      var e = m.e[y][x];
      result += e.abssq();
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
Jmat.Matrix.insert = function(a, b, row, col) {
  var result = Jmat.Matrix.copy(a);

  for(var y = 0; y < b.h; y++) {
    for(var x = 0; x < b.w; x++) {
      var rx = x + col;
      var ry = y + row;
      if(rx >= 0 && rx < a.w && ry >= 0 && ry < a.h) {
        result.e[ry][rx] = b.e[y][x];
      }
    }
  }

  return result;
};

// QR factorization of complex matrix (with householder transformations)
// requirement: m.h >= m.w
// returns {q: Q, r: R}
// q is h*h unitary matrix
// r is h*w upper triangular matrix
Jmat.Matrix.qr = function(m) {
  /*
  Tests in console:
  var qr = Jmat.qr([[1,2],[3,'4i']]); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  var qr = Jmat.qr([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  var qr = Jmat.qr(Jmat.Matrix(3,3,12,-51,4,6,167,-68,-4,24,-41)); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  degenerate matrix:
  var qr = Jmat.qr([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]); console.log(Matrix.render(qr.q)); console.log(Matrix.render(qr.r)); console.log(Matrix.render(qr.q.mul(qr.r)));
  */

  if(m.h < m.w) return null;

  var M = Jmat.Matrix;
  var C = Jmat.Complex;
  var T = M.transjugate;
  var real = Jmat.Matrix.isReal(m);
  var a = M.copy(m);
  var h = a.h;
  var w = a.w;
  var v = []; // the reflection vectors
  var taus = []; // the multiplication values
  for(var k = 0; k < w; k++) {
    var x = M.submatrix(a, k, h, k, k + 1);
    var s = x.e[0][0].eqr(0) ? C(-1) : C.sign(x.e[0][0]);
    v[k] = M.identity(h - k, 1).mulc(s.mul(M.norm(x))).add(x);
    var normv = M.norm(v[k]);
    var degenerate = normv.eqr(0);
    var tau;
    if(degenerate) {
      // In case of a degenerate column, do no reflection by setting tau to zero
      tau = C(0);
    } else {
      v[k] = v[k].divc(normv);
      tau = C(2);
      if(!real) {
        var xhv = M.mul(M.transjugate(x), v[k]);
        var vhx = M.mul(M.transjugate(v[k]), x);
        tau = xhv.e[0][0].div(vhx.e[0][0]).addr(1);
      }
    }

    taus[k] = tau;
    var as = M.submatrix(a, k, h, k, w);
    as = as.sub(v[k].mul(T(v[k])).mul(as).mulc(tau));
    a = M.insert(a, as, k, k);
  }
  var r = a;

  var q = M.identity(h, h);
  for(var k = w - 1; k >= 0; k--) {
    var z = M.submatrix(q, k, h, 0, h);
    var z = M.submatrix(q, k, h, 0, h);
    z = z.sub(v[k].mul(T(v[k])).mul(z).mulc(taus[k]));
    q = M.insert(q, z, k, 0);
  }

  return { q: q, r: r };
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
  return [l1, l2];
};

// explicit algebraic formula for eigenvalues and vectors of 2x2 matrix
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

// Returns the eigenvalues of m in an array, from largest to smallest
Jmat.Matrix.eigval = function(m) {
  var M = Jmat.Matrix;
  if(m.w != m.h || m.w < 1) return null;
  var n = m.w;
  if(n == 1) return [m.e[0][0]];
  if(n == 2) return M.eigval22(m);

  // using the QR algorithm
  var a = M.copy(m); //will contain the eigenvalues on the diagonal

  // Naive implicit QR without shifts, does not work in many cases, left commented out for demo purpose only
  // for(var i = 0; i < 30; i++) {
  //  var qr = M.qr(a);
  //  a = M.mul(qr.r, qr.q); // RQ instead of QR: A_k -> QR, A_(k+1) = RQ
  // }

  // QR with double shifting. This because with single shift or no shift, it does not support complex eigenvalues of real matrix, e.g. [[1,-1],[5,-1]]
  // TODO: this is slow, optimize like with Hessenberg form
  var id = M.identity(n, n);
  var good = undefined; // good a (no NaNs)
  for(var i = 0; i < 15; i++) {
    // var s = a.e[a.h - 1][a.w - 1]; //value that would be chosen for single shift
    // double shift: choose two sigma's, the eigenvalues of the bottom right 2x2 matrix (for which we have the explicit solution)
    var l = M.eigval22(M.submatrix(a, a.h - 2, a.h, a.w - 2, a.w));

    var si0 = id.mulc(l[0]);
    a = a.sub(si0);
    var qr = M.qr(a);
    a = M.mul(qr.r, qr.q).add(si0);

    var si1 = id.mulc(l[1]);
    a = a.sub(si1);
    qr = M.qr(a);
    a = M.mul(qr.r, qr.q).add(si1);
    if(M.isValid(a)) {
      good = a;
    } else {
      // Sometimes QR starts returning NaNs, e.g. if some values are too small. E.g. happens with [[1,2,3],[4,5,6],[7,8,9]] at input. Breaking out if i is large enough is still ok.
      if(i < 3) return null; // Too bad.
      if(!good) good = a;
      break;
    }
  }

  a = good;
  var result = [];
  for(var i = 0; i < m.w; i++) result.push(a.e[i][i]);

  return result;
};

// Returns the eigenvectors and eigenvalues of m as { l: eigenvalues, v: eigenvectors }
// eigenvalues as n*1 column vector, eigenvectors as n*n matrix
// for each column of v and corresponding eigenvalue: A*v = l*v (l represents lambda, A is m)
// opt_normalize is how to normalize the eigenvectors: 0: don't (length unspecified), 1: last element equals 1, 2: length 1. The default is "1".
Jmat.Matrix.eig = function(m, opt_normalize) {
  var M = Jmat.Matrix;
  var normalize_mode = (opt_normalize == undefined) ? 1 : opt_normalize;
  if(m.w != m.h || m.w < 1) return null;
  var n = m.w;
  if(n == 1) return M.eig11(m);
  if(n == 2) return M.eig22(m);
  /*
  Checks in console:
  var result = Jmat.Matrix.eig(Jmat.Matrix(2,2,1,2,3,4))
  Jmat.Matrix.toString(result.l) + ' \n ' + Jmat.Matrix.toString(result.v);
  */

  var val = M.eigval(m);
  var l = new M(m.w, 1);
  for(var i = 0; i < m.w; i++) l.e[i][0] = val[i];

  // Find eigenvectors by solving system of linear equations.
  // TODO: this is not very efficient...
  // Normally, the product of all the qr.q's of the loop above should give the eigenvectors, but that applies only for symmetric matrices while this is supposed to support all
  // So, instead, solve system equation (A - lambda * I) * x = 0, but with last element of 0 set to 1, and bottom row of (A - lambda * I) set to 0,0,...,0,1.
  // That makes the system solvable (TODO: not really, what if a bottom row 0 makes a whole column 0?), and makes each vector have its last element be 1.
  var v = new M(n, n);
  for(var j = 0; j < n; j++) {
    var lambda = val[j];
    var e = M.copy(m); //TODO: this makes it even slower, copy only the needed columns
    for(var i = 0; i < n; i++) e.e[i][i] = e.e[i][i].sub(lambda);
    var f = M.zero(n, 1);
    var g = M.solve(e, f);
    if(g) {
      if(normalize_mode == 2) g = g.divc(M.norm(g));
      if(normalize_mode == 1) if(!g.e[n - 1][0].eqr(0)) g = g.divc(g.e[n - 1][0]);
      for(var i = 0; i < n; i++) v.e[i][j] = g.e[i][0]; // The eigenvectors are stored as column vectors
    } else {
      // failed to find the corresponding eigenvectors, avoid crash, set values to NaN
      for(var i = 0; i < n; i++) v.e[i][j] = Jmat.Complex(NaN, NaN); // The eigenvectors are stored as column vectors
    }
  }

  return { l: l, v: v };
};


// Returns the eigen decomposition of m as { v: V, d: D }
// If M is diagonizable, M = V * D * V^(-1)
// In other words: m == result.v.mul(result.d).mul(Jmat.Matrix.inv(result.v))
// This function is very similar to Jmat.Matrix.eig. v is the same, d is the same as l but put on the diagonal of a matrix
Jmat.Matrix.evd = function(m) {
  var eig = Jmat.Matrix.eig(m);

  return { v: eig.v, d: Jmat.Matrix.diag(eig.l) };
};

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
  for(var i = 0; i < n; i++) result = result.add(a.get1(i).mul(b.get1(i).conj()));
  return result;
};


/*
@license
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
// input M, returns {u: U, s: S, v: V } such that U * W * V^T = M and S diagonal with singular values (^T means conjugate transpose here)
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
//Uses the pseudoinverse if A is not invertible
//a: input matrix, h*w size
//b: input vector, h*1 size
Jmat.Matrix.solve = function(a, b) {
  var M = Jmat.Matrix;
  if(a.h != b.h) return undefined; // input error

  if (M.isZero(b)) {
    // Homogenous system.
    var svd = M.svd(a);
    if (!Jmat.Complex.nearr(svd.s.e[svd.s.h - 1][svd.s.h - 1], 0, 1e-5)) return null; //allow large imprecision, so that it works if eigenvector calculation is a bit imprecise.
    // A is assumed singular. Last column should correspond to singular value which is zero, corresponding to a solution.
    return M.col(svd.v, a.w - 1);
  }

  var r = M.doolittle_lup_(a);
  if(!r) return null; // error
  var lu = r[0];
  var pivot = r[1];

  var n = lu.h;
  var x = M(n, 1);

  for (var k = 0; k < n; k++) {
    if (pivot[k] != k) { var temp = b.e[k][0]; b.e[k][0] = b.e[pivot[k]][0]; b.e[pivot[k]][0] = temp; }
    x.e[k][0] = b.e[k][0];
    for(var i = 0; i < k; i++) x.e[k][0] = x.e[k][0].sub(x.e[i][0].mul(lu.e[k][i]));
  }
  for (var k = n-1; k >= 0; k--) {
    if (pivot[k] != k) { var temp = b.e[k][0]; b.e[k][0] = b.e[pivot[k]][0]; b.e[pivot[k]][0] = temp; }
    for(var i = k + 1; i < n; i++) x.e[k][0] = x.e[k][0].sub(x.e[i][0].mul(lu.e[k][i]));
    if (!lu.e[k][k].eqr(0)) x.e[k][0] = x.e[k][0].div(lu.e[k][k]);
  }

  return x;
};



////////////////////////////////////////////////////////////////////////////////
/*
@license
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
