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
Jmat.BigInt: Big Integer Math

IMPORTANT NOTE: This is not intended for crypto! E.g. JavaScript's Math.random is used, there is no protection against side-channel attacks, there may be bugs, etc...

For convenience:
*) Jmat.BigInt is defined here. Its functions operate on Jmat.BigInt objects only (unless where other types are applicable).
*) Jmat.BigIntC is the same, but its functions are enriched to also take other data types as arguments (internally converting them to Jmat.BigInt ojbects)
*) BigInt in the global namespace, exported by jmat.js, is actually Jmat.BigIntC.

Jmat.BigIntC functions can work with four types of data: Jmat.BigInt (BigInt) objects, strings, arrays, and plain JS numbers.
Any radix is supported, but the native radix is Jmat.BigInt.ARRAYBASE_.

Overview of some functionality:
-elementary operators: BigInt.add, BigInt.sub, BigInt.mul, BigInt.div, BigInt.half
-base conversions: BigInt.convertBase, BigInt.convertArrayBase, BigInt.convertStringBase
-logic operators: BigInt.bitand, BigInt.bitor, BigInt.bitxor, BigInt.bitnot, BigInt.lshift, BigInt.rshift, ...
-tests: BigInt.isOdd, BigInt.isEven, BigInt.perfectsquare, BigInt.{eq, neq, lt, gt, lte, gte}
-compare: BigInt.compare, BigInt.min, BigInt.max
-power and logarithms: BigInt.sqrt, BigInt.root, BigInt.pow, BigInt.modpow, BigInt.log2, BigInt.log10, ...
-primes: BigInt.isPrime, BigInt.factorize, BigInt.nextPrime, BigInt.previousPrime, BigInt.nearestPrime
-Euclidean Algorithm: BigInt.gcd, BigInt.egcd
-factorial: BigInt.factorial, BigInt.primorial
*/

/*
Constructor, but also usable without new as factory function.
All parameters are optional.
*/
Jmat.BigInt = function(a, b, minus) {
  if(this instanceof Jmat.BigInt) {
    this.a = a || []; // least significant in rightmost element
    this.radix = b || Jmat.BigInt.ARRAYBASE_;
    this.minus = minus || false; //if true, represents negative integer
  } else {
    // No keyword "new" in front, use the convenience factory function instead
    return Jmat.BigInt.make(a, b, minus); // This supports several argument types
  }
};

// List all functions without underscore: for(var i in Jmat.BigInt) if(i.indexOf('_') == -1) console.log(i)

// default array base, must be power of two (otherwise things like 'and' and 'rshift' break). Larger is faster.
// 32768 is the maximum supported, above that, JS cannot do bit operations on product of two numbers (breaking leemondiv_).
Jmat.BigInt.ARRAYBASE_ = 32768;
Jmat.BigInt.ARRAYBASE_BITS_ = 15; // log2(Jmat.BigInt.ARRAYBASE_)
// the default base for numbers given as string. TODO: support hex strings if they start with 0x
Jmat.BigInt.STRINGBASE_ = 10;

//Typically, a is array, string, number or BigInt, and b is base of output (that of input may be different, e.g. always 10 for string).
//minus is only used for array or no input
Jmat.BigInt.make = function(a, b, minus) {
  if(a == undefined) return new Jmat.BigInt(a, b, minus);
  if(typeof a == 'number') return Jmat.BigInt.fromInt(a, b);
  if(typeof a == 'string') return Jmat.BigInt.parse(a, Jmat.BigInt.STRINGBASE_, b);
  if(a.length != undefined) return new Jmat.BigInt(a, b, minus);
  if(a.re != undefined) return Jmat.BigInt.fromInt(a.re, b);
  if(a.w != undefined) return Jmat.BigInt.fromInt(a.w, b);
  if(!b || a.radix == b) return Jmat.BigInt.copy(a);
  return Jmat.BigInt.convertBase(a, b);
};


// casts to a BigInt if not already one.
// if opt_base is undefined, does not care about radix (won't base-convert if v already BigInt). If given, casts or converts to one with that base
//minus is only used for array or no input
Jmat.BigInt.cast = function(v, opt_base, opt_minus) {
  if(v == undefined) return new Jmat.BigInt(undefined, opt_base, opt_minus);
  if(v.a != undefined && (!opt_base || v.radix == opt_base)) return v;
  if(v.radix && opt_base) return Jmat.BigInt.convertBase(v, opt_base);
  return Jmat.BigInt(v, opt_base, opt_minus);
};

//aka clone
Jmat.BigInt.copy = function(v) {
  return new Jmat.BigInt(v.a.slice(0), v.radix, v.minus);
};

Jmat.BigInt.parse = function(s, opt_stringbase, opt_outbase) {
  var stringbase = opt_stringbase || Jmat.BigInt.STRINGBASE_;
  var outbase = opt_outbase || Jmat.BigInt.ARRAYBASE_;
  var minus = false;
  if(s[0] == '-') {
    minus = true;
    s = s.substr(1);
  }
  var a = Jmat.BigInt.stringToArray(s, outbase, stringbase);
  return new Jmat.BigInt(a, outbase, minus);
};

Jmat.BigInt.toString = function(value, opt_base) {
  if(!value || !value.a) return 'invalid';
  var s = (value.minus ? '-' : '');
  return s + Jmat.BigInt.arrayToString(value.a, value.radix, opt_base);
};
Jmat.BigInt.prototype.toString = function(opt_base) {
  return Jmat.BigInt.toString(this, opt_base);
};

//By default, this converts ARRAYBASE_ array to base-10 string. The optional abase and sbase parameters change this.
//must be positive (no '-' in the string)
Jmat.BigInt.stringToArray = function(s, opt_abase, opt_sbase) {
  var abase = opt_abase || Jmat.BigInt.ARRAYBASE_;
  var sbase = opt_sbase || Jmat.BigInt.STRINGBASE_;

  var result = [];
  for(var i = 0; i < s.length; i++) {
    var v = Jmat.BigInt.di_(s[i]);
    result[i] = v;
  }

  return Jmat.BigInt.convertArrayBase(result, sbase, abase);
};

Jmat.BigInt.arrayToString = function(v, opt_abase, opt_sbase) {
  var abase = opt_abase || Jmat.BigInt.ARRAYBASE_;
  var sbase = opt_sbase || Jmat.BigInt.STRINGBASE_;
  v = Jmat.BigInt.convertArrayBase(v, abase, sbase);

  var result = '';
  for(var i = 0; i < v.length; i++) {
    result += Jmat.BigInt.d_[v[i]];
  }
  if(v.length == 0) result = '0';

  return result;
};

// This can be used by functions from within convertArrayBase to avoid recursive call if a binary BigInt is needed...
Jmat.BigInt.getInBaseTwoWithoutConvert_ = function(v) {
  var a = [];
  while(v > 0) {
    a.push(v & 1);
    v = v >> 1;
  }
  Jmat.BigInt.mirror_(a);
  return new Jmat.BigInt(a, 2);
};

//may return input object if bases are equal
Jmat.BigInt.convertArrayBase = function(s, from, to, opt_powrcache_) {
  if(s.length > 8 && to == 10 && Jmat.Real.isPOT(from)) {
    // Using larger to-bases for the non-linear algos gives huge speedup.
    // TODO: do this for all such to-bases rather than hardcoded for 10.
    var a = Jmat.BigInt.convertArrayBase(s, from, 1000000);
    return Jmat.BigInt.convertArrayBase(a, 1000000, 10);
  }
  var R = Jmat.Real;
  var B = Jmat.BigInt;
  if(from == to) return s;
  var r = [];
  s = B.maybecopystrip_(s);

  // O(n) algorithm if both bases are powers of two
  if(from <= 2147483648 && to <= 2147483648 && R.isPOT(from) && R.isPOT(to)) {
    var fbits = R.ilog2(from);
    var tbits = R.ilog2(to);
    var rlen = Math.ceil(s.length * fbits / tbits);
    for(var i = 0; i < rlen; i++) r[i] = 0;

    if(from < to) {
      var pos = rlen - 1; // element pos in result array
      var bp = 0; // bit pos into result element (from right)
      for(var i = s.length - 1; i >= 0; i--) {
        if(bp + fbits > tbits) {
          var a = tbits - bp;
          if(a != 0) r[pos] |= ((s[i] & ((1 << a) - 1)) << bp);
          pos--;
          r[pos] = (s[i] >> a);
          bp = fbits - a;
        } else {
          r[pos] |= (s[i] << bp);
          bp += fbits;
        }
      }
    } else {
      var mask = to - 1;
      var pos = s.length - 1; // element pos in input array
      var bp = 0; // bit pos into input element (from right)
      for(var i = rlen - 1; i >= 0; i--) {
        if(bp + tbits > fbits) {
          var a = fbits - bp;
          r[i] = s[pos] >> bp;
          pos--;
          r[i] |= ((s[pos] & mask) << a);
          bp = tbits - a;
        } else {
          r[i] |= ((s[pos] >> bp) & mask);
          bp += tbits;
        }
      }
    }

    B.stripInPlace_(r);
    return r;
  }

  // If one base is a power of the other base, an O(n) algorithm can be used
  var p = R.isPowerOf(from, to);
  if(p) {
    for(var i = 0; i < s.length; i++) {
      var el = s[i];
      for(var j = 0; j < p; j++) {
        r.push(Math.floor(el / (p - j - 1)) % to);
      }
      for(var j = p - 1; j >= 0; j--) {
        r[i * p + j] = el % to;
        el = Math.floor(el / to);
      }
    }
    B.stripInPlace_(r);
    return r;
  }
  p = R.isPowerOf(to, from);
  if(p) {
    var num = Math.ceil(s.length / p);
    for(var i = 0; i < num; i++) {
      var m = 1;
      r[num - 1 - i] = 0;
      for(var j = 0; j < p; j++) {
        var index = s.length - 1 - (i * p) - j;
        if(index < 0) break;
        r[num - 1 - i] += s[index] * m;
        m *= from;
      }
    }
    B.stripInPlace_(r);
    return r;
  }

  if(s.length > 8) {
    // Divide and conquer radix conversion
    var h = R.idiv(s.length, 2);
    var high = [];
    var low = [];
    for(var i = 0; i < h; i++) high[i] = s[i];
    for(var i = h; i < s.length; i++) low[i - h] = s[i];

    var cache = opt_powrcache_ || [];

    var high2 = B.convertArrayBase(high, from, to, cache);
    var low2 = B.convertArrayBase(low, from, to, cache);

    var e = s.length - h;
    var p = cache[e];
    if (!p) {
      var eh = R.idiv(e, 2);
      var el = e - eh;
      var ph = cache[eh];
      var pl = cache[el];
      // getInBaseTwoWithoutConvert_ is used to make sure the pow function does not call us, use more primitive base convert
      if(ph && pl) p = ph.mul(pl);
      else if(ph) p = ph.mul(B.pow(B(from, to), B.getInBaseTwoWithoutConvert_(el)));
      else if(pl) p = pl.mul(B.pow(B(from, to), B.getInBaseTwoWithoutConvert_(eh)));
      else p = B.pow(B(from, to), B.getInBaseTwoWithoutConvert_(e));
    }
    cache[e] = p;

    return B(low2, to).add(B(high2, to).mul(p)).a;
  }


  for(var i = 0; i < s.length; i++) {
    r = B.baseloop_(r, 0, from, [], 0, 1, s[i], to, false);
  }
  if(r.length == 0) r = [0];
  return r;
};
Jmat.BigInt.convertStringBase = function(s, from, to) {
  if(from == to) return s;
  var a = Jmat.BigInt.stringToArray(s, to, from);
  return Jmat.BigInt.arrayToString(a, to, to);
};

// Returns a + b in the simplest case
// The basic loop for add, with potential shift, mul, ... a and b are arrays, not BigInt objects. No parameters are optional, use [], 0 or 1.
// Supports also e.g. multiplying and adding a plain number to a single array: Jmat.BigInt.baseloop_(array, 0, mul, [], 0, 1, add, base, false)
// Also ensures the result is stripped (no leading zeroes), to avoid some calculations accidently becoming slower and slower when using many subtracts
// Calculates (a << ashift) * amul + (b << bshift) * bmul + overflow, with the arrays a and b in the given base, and the shifting meaning shifting whole digits of that base.
// If keep_asize is false, strips the result of leading zeroes. True should be the default value to use.
// If keep_asize is true, will ensure the result array size is as big as input a's array size (unless impossible due to leading non-zero value), which may cause it to add leading zeroes, or remove some.
// The result must be positive, that is, when subtracting, ensure b <= a.
// For example: to multiply a with 5, do: baseloop_(a, 0, 5, [], 0, 0, 0, radix, false)
// For example: to calculate a - b, do: baseloop_(a, 0, 1, b, 0, -1, 0, radix, false)
// For example: to calculate a << 5 + b, do: baseloop_(a, 5, 1, b, 0, 1, 0, radix, false)
Jmat.BigInt.baseloop_ = function(a, ashift, amul, b, bshift, bmul, overflow, base, keep_asize) {
  var result = [];
  var l = Math.max(a.length + ashift, b.length + bshift);
  var asize = a.length;

  if(!b.length && !ashift) {
    for(var i = 0; i < l; i++) {
      overflow += amul * (a[a.length - i - 1]);
      result[i] = (overflow + base) % base; // to also support negative values
      overflow = Math.floor(overflow / base);
    }
  } else if(ashift == 0 && bshift == 0 && amul == 1 && bmul == 1) {
    for(var i = 0; i < l; i++) {
      if(i < a.length) overflow += a[a.length - i - 1];
      if(i < b.length) overflow += b[b.length - i - 1];
      result[i] = (overflow + base) % base; // to also support negative values
      overflow = Math.floor(overflow / base);
    }
  } else {
    for(var i = 0; i < l; i++) {
      if(i >= ashift && i < a.length + ashift) overflow += amul * (a[a.length - i - 1 + ashift]);
      if(i >= bshift && i < b.length + bshift) overflow += bmul * (b[b.length - i - 1 + bshift]);
      result[i] = (overflow + base) % base; // to also support negative values
      overflow = Math.floor(overflow / base);
    }
  }
  if(overflow == Infinity) throw 'infinite overflow';
  while(overflow > 0) { //"if" would suffice but with while it supports larger than base values in the array, so that e.g. [0] + [100] with base 2 converts 100 to binary [1, 1, 0, 0, 1, 0, 0]
    result.push((overflow + base) % base);
    overflow = Math.floor(overflow / base);
  }
  if(keep_asize) {
    while(result.length > asize && result[result.length - 1] == 0) result.length--;
    while(result.length < asize) result.push(0);
  } else {
    // strip leading zeroes
    while(result.length > 1 && result[result.length - 1] == 0) result.length--;
  }
  Jmat.BigInt.mirror_(result);
  return result;
};

Jmat.BigInt.intToArray = function(i, opt_base) {
  // Above JS float's integer precision. A user may think the literal is a correct big number, but JS silently changes it to JS float. Loudly show that, do not return actual results on this.
  if(i > Jmat.Real.BIGGESTJSINT) throw 'too large integer literal for JS'; //return undefined;
  var base = opt_base || Jmat.BigInt.ARRAYBASE_;
  var result = [];
  while(i > 0) {
    result.push(i % base);
    i = Math.floor(i / base);
  }
  if(result.length == 0) result = [0];
  if(result.length > 1) Jmat.BigInt.mirror_(result); //the if is there so that the Jmat.BigInt.ONE etc... assignments can be there without already having mirror_ function.
  return result;
};

// If the number is too big to represent in JS precision, returns float approximation or if even too high for that, Infinity.
Jmat.BigInt.arrayToInt = function(v, opt_base) {
  var base = opt_base || Jmat.BigInt.ARRAYBASE_;
  var result = 0;
  var m = 1;
  for(var i = 0; i < v.length; i++) {
    result += m * v[v.length - 1 - i];
    m *= base;
    if(result == Infinity || m == Infinity) return Infinity; //avoid it becoming NaN
  }
  return result;
};

Jmat.BigInt.fromInt = function(i, opt_base) {
  var base = opt_base || Jmat.BigInt.ARRAYBASE_;
  var minus = false;
  if(i < 0) {
    minus = true;
    i = -i;
  }
  return new Jmat.BigInt(Jmat.BigInt.intToArray(i, base), base, minus);
};

// If the number is too big to represent in JS precision, returns float approximation or if even too high for that, Infinity.
Jmat.BigInt.toInt = function(v) {
  return Jmat.BigInt.arrayToInt(v.a, v.radix) * (v.minus ? -1 : 1);
};
Jmat.BigInt.prototype.toInt = function() {
  return Jmat.BigInt.arrayToInt(this.a, this.radix) * (this.minus ? -1 : 1);
};

Jmat.BigInt.ZERO = Jmat.BigInt(0);
Jmat.BigInt.ONE = Jmat.BigInt(1);
Jmat.BigInt.TWO = Jmat.BigInt(2);

//may return input if bases are equal
//only supports BigInt object at input. For arrays or strings use convertArrayBase or convertStringBase.
Jmat.BigInt.convertBase = function(s, to) {
  if(s.radix == to) return s;
  return new Jmat.BigInt(Jmat.BigInt.convertArrayBase(s.a, s.radix, to), to, s.minus);
};

Jmat.BigInt.cloneArray = function(v) {
  return v.slice(0);
};

Jmat.BigInt.cloneArrayTo = function(v, target) {
  target.length = v.length;
  for(var i = 0; i < v.length; i++) target[i] = v[i];
};

// left shift by s digits (based on a's radix)
// s is regular JS number
Jmat.BigInt.lshift_radix = function(a, s) {
  if(s == 0) return a;
  if(s < 0) return Jmat.BigInt.rshift_radix(a, -s);
  var result = Jmat.BigInt.copy(a);
  for(var i = 0; i < s; i++) result.a.push(0);
  return result;
};
Jmat.BigInt.prototype.lshift_radix = function(s) {
  return Jmat.BigInt.lshift_radix(this, s);
};

// right shift by s digits (based on a's radix)
// s is regular JS number
Jmat.BigInt.rshift_radix = function(a, s) {
  if(s == 0) return a;
  if(s < 0) return Jmat.BigInt.lshift_radix(a, -s);
  var result = Jmat.BigInt.copy(a);
  result.a = result.a.slice(0, -s);
  return result;
};
Jmat.BigInt.prototype.rshift_radix = function(s) {
  return Jmat.BigInt.rshift_radix(this, s);
};

//shift left by b bits
//a is bigint, b is regular js number (integer)
//does not treat negative as two's complement
Jmat.BigInt.lshift = function(a, b) {
  var B = Jmat.BigInt;
  if(b == 0) return a;
  if(b < 0) return B.rshift(a, -b);
  a = B.cast(a, Jmat.BigInt.ARRAYBASE_); //ensure power of two base
  var result = new B([], a.radix);

  var byteshift = Math.floor(b / B.ARRAYBASE_BITS_);
  var bitshift = b % B.ARRAYBASE_BITS_;
  if(bitshift == 0) {
    for(var i = 0; i < a.a.length; i++) result.a[i] = a.a[i];
    for(var i = 0; i < byteshift; i++) result.a.push(0);
  } else {
    result.a = [];
    // if byte = xxxxxxxx and bitshift is e.g. 3, then lmask and rmask are indicated by: lllllrrr (with lsb on the right),
    // so lmask would be rmask 7, rmask 248 (0b11111000)
    var lmask = ((1 << (B.ARRAYBASE_BITS_ - bitshift)) - 1) << bitshift;
    var rmask = (1 << bitshift) - 1;
    for(var i = 0; i <= a.a.length; i++) {
      //xxxxxxxx xxxxxxxxx aaaaaaaa xxxxxxxx
      //xxxxxxxx xxxxxxaaa aaaaaxxx xxxxxxxx
      var al = (i > 0) ? a.a[i - 1] : 0;
      var ar = (i < a.a.length) ? a.a[i] : 0;
      var r = ((ar << bitshift) >> B.ARRAYBASE_BITS_) & rmask;
      var l = (al << bitshift) & lmask;
      var v = l | r;
      if(v || result.a.length) result.a.push(v);
    }
    for(var i = 0; i < byteshift; i++) result.a.push(0);
    if(result.a.length == 0) result.a = [0];
  }

  if(a.minus) result = result.neg();
  return result;
};
Jmat.BigInt.prototype.lshift = function(b) {
  return Jmat.BigInt.lshift(this, b);
};

//shift right by b bits
//a is bigint, b is regular js number (integer)
//does not treat negative as two's complement
Jmat.BigInt.rshift = function(a, b) {
  var B = Jmat.BigInt;
  if(b == 0) return a;
  if(b < 0) return B.lshift(a, -b);
  a = B.cast(a, Jmat.BigInt.ARRAYBASE_); //ensure power of two base
  var result = new B([], a.radix);

  var byteshift = Math.floor(b / B.ARRAYBASE_BITS_);
  var bitshift = b % B.ARRAYBASE_BITS_;
  if(bitshift == 0) {
    result.a = [];
    for(var i = 0; i < a.a.length - byteshift; i++) result.a.push(a.a[i]);
    if(result.a.length == 0) result.a = [0];
  } else {
    result.a = [];
    // if byte = xxxxxxxx and bitshift is e.g. 3, then lmask and rmask are indicated by: lllrrrrr (with lsb on the right),
    // so lmask would be rmask 31, rmask 224
    var lmask = ((1 << bitshift) - 1) << (B.ARRAYBASE_BITS_ - bitshift);
    var rmask = (1 << (B.ARRAYBASE_BITS_ - bitshift)) - 1;
    for(var i = 0; i < a.a.length - byteshift; i++) {
      //xxxxxxxx xxxxxxxxx aaaaaaaa xxxxxxxx
      //xxxxxxxx xxxxxxxxx xxxaaaaa aaaxxxxx
      var al = (i > 0) ? a.a[i - 1] : 0;
      var ar = (i < a.a.length) ? a.a[i] : 0;
      var l = ((al << B.ARRAYBASE_BITS_) >> bitshift) & lmask;
      var r = (ar >> bitshift) & rmask;
      result.a.push(l | r);
    }
    if(result.a.length == 0) result.a = [0];
  }

  if(a.minus) result = result.neg();
  B.stripInPlace_(result.a);
  return result;
};
Jmat.BigInt.prototype.rshift = function(b) {
  return Jmat.BigInt.rshift(this, b);
};

// result is always positive. If a or b is negative, the two's complement representation is taken.
// This allows using and for modulo and have correct positive result for negative input.
Jmat.BigInt.bitand = function(a, b) {
  var B = Jmat.BigInt;
  //ensure power of two base
  a = B.cast(a, B.ARRAYBASE_);
  b = B.cast(b, B.ARRAYBASE_);

  if(a.minus || b.minus) {
    // even though internal representation isn't so, use twos complement
    var n = Math.max(B.getNumBits(a), B.getNumBits(b));
    if(a.minus) a = a.addr(1).bitnot(n);
    if(b.minus) b = b.addr(1).bitnot(n);
  }

  var result = new B([], B.ARRAYBASE_);
  var num = Math.min(a.a.length, b.a.length);
  for(var i = 0; i < num; i++) {
    var ax = a.a[a.a.length - 1 - i] || 0;
    var bx = b.a[b.a.length - 1 - i] || 0;
    result.a[num - 1 - i] = ax & bx;
  }

  return result;
};
Jmat.BigInt.prototype.bitand = function(b) {
  return Jmat.BigInt.bitand(this, b);
};

// same as bitand, but with b a regular JS number (max 31 bits), and returning a regular JS number
Jmat.BigInt.bitandr = function(a, b) {
  var B = Jmat.BigInt;
  if(a.minus) {
    // TODO: make also this case faster
    return B.bitand(a, B(b)).toInt();
  }
  a = B.cast(a, (b < B.ARRAYBASE_) ? B.ARRAYBASE_ : 2147483648);
  return a.a[a.a.length - 1] & b;
};

// result is always positive. If a or b is negative, their two's complement representation is taken.
Jmat.BigInt.bitor = function(a, b) {
  var B = Jmat.BigInt;
  //ensure power of two base
  a = B.cast(a, B.ARRAYBASE_);
  b = B.cast(b, B.ARRAYBASE_);

  if(a.minus || b.minus) {
    // even though internal representation isn't so, use twos complement
    var n = Math.max(B.getNumBits(a), B.getNumBits(b));
    if(a.minus) a = a.addr(1).bitnot(n);
    if(b.minus) b = b.addr(1).bitnot(n);
  }

  var result = new B([], B.ARRAYBASE_);
  var num = Math.max(a.a.length, b.a.length);
  for(var i = 0; i < num; i++) {
    var ax = a.a[a.a.length - 1 - i] || 0;
    var bx = b.a[b.a.length - 1 - i] || 0;
    result.a[num - 1 - i] = ax | bx;
  }

  return result;
};
Jmat.BigInt.prototype.bitor = function(b) {
  return Jmat.BigInt.bitor(this, b);
};

// result is always positive. If a or b is negative, their two's complement representation is taken.
Jmat.BigInt.bitxor = function(a, b) {
  var B = Jmat.BigInt;
  //ensure power of two base
  a = B.cast(a, B.ARRAYBASE_);
  b = B.cast(b, B.ARRAYBASE_);

  if(a.minus || b.minus) {
    // even though internal representation isn't so, use twos complement
    var n = Math.max(B.getNumBits(a), B.getNumBits(b));
    if(a.minus) a = a.addr(1).bitnot(n);
    if(b.minus) b = b.addr(1).bitnot(n);
  }

  var result = new B([], B.ARRAYBASE_);
  var num = Math.max(a.a.length, b.a.length);
  for(var i = 0; i < num; i++) {
    var ax = a.a[a.a.length - 1 - i] || 0;
    var bx = b.a[b.a.length - 1 - i] || 0;
    result.a[num - 1 - i] = ax ^ bx;
  }

  return result;
};
Jmat.BigInt.prototype.bitxor = function(b) {
  return Jmat.BigInt.bitxor(this, b);
};

// Inverts all the bits of a, and treats as if it's in two's complement
// So gives -(a+1) as result. Result always has opposite sign
// If you want to only invert bits and keep positive, do one of the following:
// -use bitnot
// -use bitxor with as many 1's as you want to invert.
// -add N to the result, where N is 1<<n, where n is amount of bits
Jmat.BigInt.bitneg = function(a) {
  var B = Jmat.BigInt;
  a = B.cast(a, Jmat.BigInt.ARRAYBASE_);

  a = a.addr(1);
  a = a.neg();

  return a;
};
Jmat.BigInt.prototype.bitneg = function() {
  return Jmat.BigInt.bitneg(this);
};

//negates all bits up to the one more than the msb. Ignores (but keeps) sign.
//if opt_bits isn't given, abs of result is always larger: the 0 in front of the msb becomes 1.
//if opt_bits is given, instead negates that many bits
//see bitneg and its comment for alternatives.
Jmat.BigInt.bitnot = function(a, opt_bits) {
  if(opt_bits <= 0) return a;
  var B = Jmat.BigInt;
  a = B.cast(a, B.ARRAYBASE_);
  var ar = B.maybecopystrip_(a.a);
  var result;
  var xinv = B.ARRAYBASE_ - 1; //xor invertor

  if(opt_bits != undefined) {
    result = B(0);
    var b = opt_bits % B.ARRAYBASE_BITS_;
    var n = Math.ceil(opt_bits / B.ARRAYBASE_BITS_);
    var mask = (1 << b) - 1;
    for(var i = 0; i < n; i++) {
      var x = ar[ar.length - n + i] || 0;
      result.a[i] = x ^ xinv;
      if(i == 0 && b != 0) result.a[i] &= mask;
    }
  } else {
    if(!ar.length) return a;
    var b = Jmat.Real.ilog2(ar[0]) + 1;
    if(b == 0) {

      result = B(1); // the new MSB
      for(var i = 0; i < ar.length; i++) result.a.push(ar[i] ^ xinv);
    } else {
      var mask = (2 << b) - 1; // 2<< instead of 1<< because a new MSB gets added
      result = B(0);
      result.a[0] = (ar[0] ^ xinv) & mask;
      for(var i = 1; i < ar.length; i++) result.a[i] = ar[i] ^ xinv;
    }
  }
  result.minus = a.minus;

  return result;
};
Jmat.BigInt.prototype.bitnot = function(opt_bits) {
  return Jmat.BigInt.bitnot(this, opt_bits);
};

// does NOT copy the array, on purpose (cheap operation). May return the input object.
Jmat.BigInt.abs = function(a) {
  return a.abs();
};
Jmat.BigInt.prototype.abs = function() {
  if(this.minus) {
    return new Jmat.BigInt(this.a, this.radix, false);
  }
  return this;
};

// dist, cheb and manhattan all return regular real JS numbers for all types. In some types they are all the same, but not for e.g. Complex or Matrix.
// Euclidean distance
Jmat.BigInt.dist = function(a, b) {
  return a.sub(b).abs();
};
//Chebyshev distance
Jmat.BigInt.cheb = function(a, b) {
  return Jmat.BigInt.dist(a, b);
};
//Manhattan distance
Jmat.BigInt.manhattan = function(a, b) {
  return Jmat.BigInt.dist(a, b);
};

// does NOT copy the array, on purpose (cheap operation).
Jmat.BigInt.neg = function(a) {
  return a.neg();
};
Jmat.BigInt.prototype.neg = function() {
  return new Jmat.BigInt(this.a, this.radix, !this.minus);
};

// result is regular JS number, -1 or 1 (does NOT return 0 for zero, use "sign" for that)
Jmat.BigInt.getSign = function(a) {
  return a.getSign();
};
Jmat.BigInt.prototype.getSign = function() {
  return this.minus ? -1 : 1;
};

Jmat.BigInt.nonZero = function(a) {
  return a.nonZero();
};
Jmat.BigInt.prototype.nonZero = function() {
  for(var i = 0; i < this.a.length; i++) {
    if(this.a[i]) return true;
  }
  return false;
};

// result is regular JS number, -1, 0 or 1 (slower than "getSign")
Jmat.BigInt.sign = function(a) {
  return a.sign();
};
Jmat.BigInt.prototype.sign = function() {
  return this.nonZero() ? this.getSign() : 0;
};

Jmat.BigInt.add = function(a, b) {
  return a.add(b);
};
Jmat.BigInt.prototype.add = function(b) {
  if(this.minus != b.minus)  {
    return this.sub(b.neg());
  }
  b = Jmat.BigInt.cast(b, this.radix); //must have same radix

  return new Jmat.BigInt(Jmat.BigInt.baseloop_(this.a, 0, 1, b.a, 0, 1, 0, this.radix, false), this.radix, this.minus);
};
Jmat.BigInt.addr = function(a, b) {
  return a.addr(b);
};
Jmat.BigInt.prototype.addr = function(b) {
  // TODO: make more efficient by usign baseloop directly
  return this.add(Jmat.BigInt.fromInt(b));
};

Jmat.BigInt.sub = function(a, b) {
  return a.sub(b);
};
Jmat.BigInt.prototype.sub = function(b) {
  if(this.minus != b.minus)  {
    return this.add(b.neg());
  }
  b = Jmat.BigInt.cast(b, this.radix); //must have same radix

  if(this.abs().gte(b.abs())) {
    return new Jmat.BigInt(Jmat.BigInt.baseloop_(this.a, 0, 1, b.a, 0, -1, 0, this.radix, false), this.radix, this.minus);
  } else {
    return new Jmat.BigInt(Jmat.BigInt.baseloop_(b.a, 0, 1, this.a, 0, -1, 0, this.radix, false), this.radix, !this.minus);
  }
};
Jmat.BigInt.subr = function(a, b) {
  return a.subr(b);
};
Jmat.BigInt.prototype.subr = function(b) {
  // TODO: make more efficient by using baseloop directly
  return this.sub(Jmat.BigInt.fromInt(b));
};
Jmat.BigInt.rsub = function(a, b) {
  return a.rsub(b);
};
Jmat.BigInt.prototype.rsub = function(b) {
  return this.subr(b).neg();
};

Jmat.BigInt.mul = function(a, b) {
  return a.mul(b);
};
Jmat.BigInt.prototype.mul = function(b) {
  b = Jmat.BigInt.cast(b, this.radix); //must have same radix
  var result = Jmat.BigInt.karatsuba_(this.strip().a, b.strip().a, this.radix);
  return new Jmat.BigInt(result, this.radix, this.minus != b.minus);
};

Jmat.BigInt.mulr = function(a, b) {
  var result = new Jmat.BigInt([], a.radix, a.minus != (b < 0));
  result.a = Jmat.BigInt.baseloop_(a.a, 0, b, [], 0, 1, 0, result.radix, false);
  return result;
};
Jmat.BigInt.prototype.mulr = function(b) {
  var result = new Jmat.BigInt([], this.radix, this.minus != (b < 0));
  result.a = Jmat.BigInt.baseloop_(this.a, 0, b, [], 0, 1, 0, result.radix, false);
  return result;
};

// Karatsuba multiplication. a and b are arrays with leading zeroes stripped.
Jmat.BigInt.karatsuba_ = function(a, b, radix) {
  if(a.length <= 1 || b.length <= 1) return Jmat.BigInt.schoolmul_(a, b, radix);

  // ensure a is the one with longest length
  if(a.length < b.length) return Jmat.BigInt.karatsuba_(b, a, radix);

  // Karatsuba is only faster for really large numbers.
  if(a.length <= 20) return Jmat.BigInt.schoolmul_(a, b, radix);
  var m = Math.floor(a.length / 2);

  var mb = b.length - (a.length - m);

  var a0 = a.slice(0, m);
  var a1 = a.slice(m, a.length);
  var b0 = (mb <= 0 ? [0] : b.slice(0, mb));
  var b1 = (mb <= 0 ? b : b.slice(mb, b.length));

  if(mb <= 0) {
    // Size difference between the numbers very big. Only one is split in half.
    var x = Jmat.BigInt.karatsuba_(a0, b1, radix);
    var y = Jmat.BigInt.karatsuba_(a1, b1, radix);
    return Jmat.BigInt.baseloop_(x, a.length - m, 1, y, 0, 1, 0, radix, false);
  } else {
    var a01 = Jmat.BigInt.baseloop_(a0, 0, 1, a1, 0, 1, 0, radix, false);
    var b01 = Jmat.BigInt.baseloop_(b0, 0, 1, b1, 0, 1, 0, radix, false);

    var x = Jmat.BigInt.karatsuba_(a0, b0, radix);
    var y = Jmat.BigInt.karatsuba_(a1, b1, radix);
    var z = Jmat.BigInt.karatsuba_(a01, b01, radix);
    var xy = Jmat.BigInt.baseloop_(x, 0, 1, y, 0, 1, 0, radix, false);
    var k = Jmat.BigInt.baseloop_(z, 0, 1, xy, 0, -1, 0, radix, false);

    var s1 = a.length - m;
    var s2 = s1 * 2;
    var ky = Jmat.BigInt.baseloop_(k, s1, 1, y, 0, 1, 0, radix, false);
    return Jmat.BigInt.baseloop_(x, s2, 1, ky, 0, 1, 0, radix, false);
  }
};

// Schoolbook multiplication. a and b are arrays with leading zeroes stripped.
Jmat.BigInt.schoolmul_ = function(a, b, radix) {
  if(a.length == 1 && a[0] == 0) return [0];
  if(a.length == 1 && a[0] == 1) return b;
  if(b.length == 1 && b[0] == 0) return [0];
  if(b.length == 1 && b[0] == 1) return a;

  // ensure a is the one with longest length, the loop below is faster that way
  if(a.length < b.length) return Jmat.BigInt.schoolmul_(b, a, radix);

  var result = [0];
  var ashift = 0;
  for(var j = 0; j < b.length; j++) {
    var d = b[b.length - j - 1];
    result = Jmat.BigInt.baseloop_(a, ashift, d, result, 0, 1, 0, radix, false);
    ashift++; // left shift a
  }

  return result;
};

//returns 1 if a > b, -1 if b > a, 0 if equal.
Jmat.BigInt.compare = function(a, b) {
  return a.compare(b);
};
Jmat.BigInt.prototype.compare = function(b) {
  if(b.minus != this.minus) {
    var as = this.sign();
    var bs = b.sign();
    if(as == bs) return 0; //both signs are 0
    if(as < bs) return -1;
    return 1;
  }
  if(b.radix != this.radix) b = Jmat.BigInt.convertBase(b, this.radix);

  var l = Math.max(this.a.length, b.a.length);
  for(var i = 0; i < l; i++) {
    var ai = i - l + this.a.length;
    var bi = i - l + b.a.length;
    var av = this.a[ai] || 0; // 0 if out of bounds
    var bv = b.a[bi] || 0; // 0 if out of bounds
    if(av < bv) return -1 * this.getSign();
    if(av > bv) return 1 * this.getSign();
  }
  return 0;
};

//Same as compare, but b is regular JS number
Jmat.BigInt.comparer = function(a, b) {
  return a.comparer(b);
};
Jmat.BigInt.prototype.comparer = function(b) {
  if((b < 0) != this.minus) {
    var as = this.sign();
    var bs = Math.sign(b);
    if(as == bs) return 0; //both signs are 0
    if(as < bs) return -1;
    return 1;
  }

  var sign = this.getSign();
  b = Math.abs(b);

  var r = 0;
  for(var i = 0; i < this.a.length && r <= b; i++) {
    r *= this.radix;
    r += this.a[i];
  }
  return ((r < b) ? -1 : (r == b ? 0 : 1)) * sign;
};

Jmat.BigInt.eq = function(a, b) { return Jmat.BigInt.compare(a, b) == 0; };
Jmat.BigInt.prototype.eq = function(b) { return Jmat.BigInt.compare(this, b) == 0; };
Jmat.BigInt.eqr = function(a, b) { return Jmat.BigInt.comparer(a, b) == 0; };
Jmat.BigInt.prototype.eqr = function(b) { return Jmat.BigInt.comparer(this, b) == 0; };

Jmat.BigInt.neq = function(a, b) { return Jmat.BigInt.compare(a, b) != 0; };
Jmat.BigInt.prototype.neq = function(b) { return Jmat.BigInt.compare(this, b) != 0; };
Jmat.BigInt.neqr = function(a, b) { return Jmat.BigInt.comparer(a, b) != 0; };
Jmat.BigInt.prototype.neqr = function(b) { return Jmat.BigInt.comparer(this, b) != 0; };

Jmat.BigInt.gt = function(a, b) { return Jmat.BigInt.compare(a, b) > 0; };
Jmat.BigInt.prototype.gt = function(b) { return Jmat.BigInt.compare(this, b) > 0; };
Jmat.BigInt.gtr = function(a, b) { return Jmat.BigInt.comparer(a, b) > 0; };
Jmat.BigInt.prototype.gtr = function(b) { return Jmat.BigInt.comparer(this, b) > 0; };

Jmat.BigInt.lt = function(a, b) { return Jmat.BigInt.compare(a, b) < 0; };
Jmat.BigInt.prototype.lt = function(b) { return Jmat.BigInt.compare(this, b) < 0; };
Jmat.BigInt.ltr = function(a, b) { return Jmat.BigInt.comparer(a, b) < 0; };
Jmat.BigInt.prototype.ltr = function(b) { return Jmat.BigInt.comparer(this, b) < 0; };

Jmat.BigInt.gte = function(a, b) { return Jmat.BigInt.compare(a, b) >= 0; };
Jmat.BigInt.prototype.gte = function(b) { return Jmat.BigInt.compare(this, b) >= 0; };
Jmat.BigInt.gter = function(a, b) { return Jmat.BigInt.comparer(a, b) >= 0; };
Jmat.BigInt.prototype.gter = function(b) { return Jmat.BigInt.comparer(this, b) >= 0; };

Jmat.BigInt.lte = function(a, b) { return Jmat.BigInt.compare(a, b) <= 0; };
Jmat.BigInt.prototype.lte = function(b) { return Jmat.BigInt.compare(this, b) <= 0; };
Jmat.BigInt.lter = function(a, b) { return Jmat.BigInt.comparer(a, b) <= 0; };
Jmat.BigInt.prototype.lter = function(b) { return Jmat.BigInt.comparer(this, b) <= 0; };

// Returns integer square root (floor), or undefined if negative
Jmat.BigInt.sqrt = function(a) {
  var B = Jmat.BigInt;

  if(a.eqr(0)) return B(0);
  if(a.minus) return undefined;

  var low = B([0], a.radix);
  var high = Jmat.BigInt.copystrip_(a);
  if(high.a.length == 1 && high.a[0] == 1) low = B([1], a.radix); //the algorithm below fails on [1]
  high.a = high.a.slice(0, Math.ceil(high.a.length / 2) + 1); // initial estimate for max of sqrt: half amount of digits
  if(a.radix > 16) high.a[0] = Math.min(high.a[0], Math.ceil(Math.sqrt(high.a[0] + 1))); // further improve estimate by setting most significant digit to sqrt of itself

  var one = B([1], a.radix);

  var result;
  for (;;) {
    var mid = low.add(high).divr(2);
    var rr = mid.mul(mid);
    var c = B.compare(rr, a);
    if(c == 0) {
      result = mid;
      break;
    }
    else if(c < 0) low = mid;
    else high = mid;
    if(B.compare(high.sub(low), one) <= 0) {
      result = low;
      break;
    }
  }

  return result;
};

// x is odd integer
Jmat.BigInt.isOdd = function(x) {
  return x.bitand(Jmat.BigInt.ONE).eqr(1);
};

// x is even integer
Jmat.BigInt.isEven = function(x) {
  return x.bitand(Jmat.BigInt.ONE).eqr(0);
};


//Finds the nth root of a (e.g. sqrt if n is 2)
Jmat.BigInt.root = function(a, n) {
  var B = Jmat.BigInt;
  var r = n.toInt();
  if(r <= 0) return undefined;
  if(r > Jmat.Real.BIGGESTJSINT) {
    if(a.eqr(0)) return B(0);
    if(B.isEven(n) && a.minus) return undefined;
    // We can return 1 because:
    // -the above checks already checked for 0 and invalid
    // -nth root of a with n >= log2(a) is smaller than 2
    // -number of bits of a is assumed to fit in a regular JS number (2^53 bits), so n is definitely bigger than that if it doesn't fit in JS int
    // -so as a result, the floor of the nth root is 1
    return B(a.minus ? -1 : 1);
  }
  return Jmat.BigInt.rootr(a, r);
};


//Finds the nth root of a (e.g. sqrt if n is 2), with n a regular JS int
Jmat.BigInt.rootr = function(a, n) {
  var B = Jmat.BigInt;
  if(n <= 0) return undefined;
  if(a.eqr(0)) return B(0);
  if(n == 1) return a;
  if(n == 2) return B.sqrt(a);
  if(Jmat.Real.isEven(n) && a.minus) return undefined;
  //if n is bigger than log2(a), the result is smaller than 2 (and given checks above, larger than 1).
  if(n > B.log2(a.abs()).toInt()) return B(a.minus ? -1 : 1);

  var low = B([0], a.radix);
  var high = Jmat.BigInt.copystrip_(a);
  if(high.a.length == 1 && high.a[0] == 1) low = B([1], a.radix); //the algorithm below fails on [1]
  high.a = high.a.slice(0, Math.ceil(high.a.length / n) + 1); // initial estimate for max of sqrt: half amount of digits
  if(a.radix > 16) high.a[0] = Math.min(high.a[0], Math.ceil(Jmat.Real.root(high.a[0] + 1, n))); // further improve estimate by setting most significant digit to root of itself

  var one = B([1], a.radix);

  var result;
  for (;;) {
    var mid = low.add(high).divr(2);
    var rr = B.powr(mid, n);
    var c = B.compare(rr, a);
    if(c == 0) {
      result = mid;
      break;
    }
    else if(c < 0) low = mid;
    else high = mid;
    if(B.compare(high.sub(low), one) <= 0) {
      result = low;
      break;
    }
  }

  return result;
};

// tests whether the number is a perfect square. Returns null if not, its square root if it is
Jmat.BigInt.perfectsquare = function(a) {
  var B = Jmat.BigInt;
  if(a.minus) return null;
  var nibble = B.bitandr(a, 15);
  // only if the last nibble is 0, 1, 4 or 9, it's a square.
  if(!(nibble == 0 || nibble == 1 || nibble == 4 || nibble == 9)) {
    return null;
  }
  var s = B.sqrt(a);
  if(s.mul(s).eq(a)) return s;
  return null;
};

//Finds the nearest perfect power <= a, returns array with [base^exponent, base, exponent], with base >= 1 and exponent >= 2
//To test if number is perfect power, simply check whether result is same as input.
//Gives smallest possible base (with largest exponent), e.g. for 64, gives 2^6, not 8^2
//if opt_next is true, then instead finds the first perfect power greater than a
//if opt_base is given, then only uses exactly opt_base as base (much simpler)
//if opt_k is given, then only uses exactly opt_k as exponent (also much simpler)
// E.g. to list first 100 perfect powers: var p = BigInt(0); for(var i = 0; i < 100; i++) { var pp = BigInt.perfectpow(p, true); console.log('' + pp); p = pp[0]; }
Jmat.BigInt.perfectpow = function(a, opt_next, opt_base, opt_k) {
  //The exponent of the result will usually be 2, most perfect powers are squares, higher exponents are much rarer.
  var B = Jmat.BigInt;
  var k = opt_k ? B.cast(opt_k) : undefined;
  var b = opt_base ? B.cast(opt_base) : undefined;
  if(a.minus) return undefined;
  if(opt_base && opt_k) return [b.pow(k), b, k]; //not useful, but consistent
  if(opt_base) {
    if(a.lter(1) && opt_next) return [b.mul(b), b, B(2)];
    var r = (typeof opt_base != 'number') ? opt_base : b.toInt();
    var l = B.logr(a, r);
    if(opt_next) l = l.addr(1);
    var p = b.pow(l);
    return [p, b, l];
  }
  if(opt_k) {
    var r = (typeof opt_k != 'number') ? opt_k : k.toInt();
    var s = B.rootr(a, r);
    if(opt_next) s = s.addr(1);
    var p = s.powr(r);
    return [p, s, k];
  }
  // Algo below doesn't work for a < 5
  if(a.ltr(5)) {
    var r = a.toInt();
    if(opt_next) {
      if(r == 0) return [B(1), B(1), B(2)];
      if(r == 4) return [B(8), B(2), B(3)];
      return [B(4), B(2), B(2)];
    } else {
      if(r == 0) return undefined;
      if(a == 4) return [B(4), B(2), B(2)];
      return [B(1), B(1), B(2)];
    }
  }
  var l = B.log2(a).toInt();
  if(opt_next) l++;

  var best = a;
  var besti = 0;
  var bestbase = B(0);
  var bestpow = B(0);

  for(var i = 2; i <= l; i++) {
    if(i > 5 && i % 5 == 0) { i++; continue; }
    if(i > 3 && i % 3 == 0) { i++; continue; }

    var s = B.rootr(a, i); // the potential base

    // For low bases (such as 2 and 3, but especially 2), we will end up having
    // this base many times, but only the highest possible power gives the nearest answer.
    // Speed things up by setting the exponent immediately to the highest possible:
    // log_s(a) (already calculated as l in case of log2)
    if(s.eqr(2)) i = l;

    if(opt_next) s = s.addr(1);
    var p = s.powr(i);
    var diff = a.sub(p).abs();
    if(diff.lt(best)) {
      best = diff;
      besti = i;
      bestbase = s;
      bestpow = p;
      if(diff.eqr(0)) break;
    }

    if(i > 2) i++; //add one more, skip next even numbers
  }

  // We want biggest, not smallest, exponent, and smallest, nog biggest, base
  if(besti == 2 || besti == 3 || besti == 5) {
    var r2 = Jmat.BigInt.perfectpow(bestbase);
    if(r2[0].eq(bestbase)) {
      bestbase = r2[1];
      besti = besti * r2[2].toInt();
    }
  }

  return [bestpow, bestbase, B(besti)];
};

// Multiplies a list of numbers recursively, to end up with factors of more equal size, ending with large factors that take advantage of sub-quadratic multiplication.
// That is faster than just multiplying the list serially. Optimized for roughly sorted list (pairs value from begin with end etc...).
Jmat.BigInt.mulmany_ = function(a) {
  if(a.length == 0) return Jmat.BigInt.ONE;
  while(a.length > 1) {
    var h = Math.floor(a.length / 2);
    var r = [];
    for(var i = 0; i < h; i++) {
      r.push(a[i].mul(a[a.length - i - 1]));
    }
    if(a.length % 2 == 1) r.push(a[h]);
    a = r;
  }

  return a[0];
};


// Can calculate factorial of 100000 in 30 seconds. The result has 1.5 million bits, and printing it in decimal then takes another minute or so.
Jmat.BigInt.factorial = function(a) {
  var B = Jmat.BigInt;
  var b = a.toInt(); //if a does not fit in regular JS number, the result is unreasonably huge anyway.
  if(b == Infinity) return undefined;

  if(b < 19) return B(Jmat.Real.factorial(b));

  var primes = Jmat.Real.eratosthenes(b);
  var f = [];

  // For each prime, it's easy to calculate the prime factors of the factorial. They are all powers of a prime, stored in f.
  // E.g. for prime 2, this matches the list "0,1,3,4,7,8,10,11,15,16,18,19,..." in OEIS, and for others it's analogous.
  for(var i = 0; i < primes.length; i++) {
    var p = primes[i];
    var c = Math.floor(b / p);
    var pw = 0;
    while(c > 0) {
      pw += c;
      c = Math.floor(c / p);
    }

    f[i] = B.powr(B(p), pw);
  }


  return B.mulmany_(f);
};

// Primorial, with a being an upper bound on the prime (not the nth prime)
Jmat.BigInt.primorial = function(a) {
  var B = Jmat.BigInt;
  var b = a.toInt(); //if a does not fit in regular JS number, the result is unreasonably huge anyway.
  if(b == Infinity) return undefined;

  var primes = Jmat.Real.eratosthenes(b);
  return B.mulmany_(primes);
};

Jmat.BigInt.nearestPrime = function(n) {
  var B = Jmat.BigInt;
  if(n.lter(2)) return B(2);
  if(n.lter(4)) return B(3);

  var p = B.previousPrime(n);
  var diff = n.sub(p);
  return B.previousPrime(n.add(diff).subr(1));
};

Jmat.BigInt.nextPrime = function(n) {
  var B = Jmat.BigInt;
  if(n.ltr(2)) return B(2);
  if(n.ltr(3)) return B(3);

  var m = n.modr(6).toInt();
  var step = 2;
  if(m == 0 || m == 5) {
    n = n.addr(m == 0 ? 1 : 2);
    step = 4;
  } else {
    n  = n.addr(5 - m);
  }
  for(;;) {
    if(B.isPrime(n)) return n;
    n = n.addr(step);
    step ^= 6; //swap step between 2 and 4
  }
};

Jmat.BigInt.previousPrime = function(n) {
  var B = Jmat.BigInt;
  // Anything below 6 does not work with the calculations below.
  if(n.lter(7)) {
    if(n.lter(2)) return undefined; // there is no lower prime
    if(n.lter(3)) return B(2);
    if(n.lter(5)) return B(3);
    return B(5);
  }

  var m = n.modr(6).toInt();
  var step = 2;
  if(m == 0 || m == 1) {
    n = n.subr(m + 1);
    step = 4;
  } else {
    n  = n.subr(m - 1);
  }
  for(;;) {
    if(B.isPrime(n)) return n;
    n = n.subr(step);
    step ^= 6; //swap step between 2 and 4
  }
};

Jmat.BigInt.min = function(a, b) {
  return a.lt(b) ? a : b;
};

Jmat.BigInt.max = function(a, b) {
  return a.gt(b) ? a : b;
};

Jmat.BigInt.primeCache_ = [];

// Sieves in a given range to collect smooth numbers. res is the quadratic residues.
// start and end are BigInt, primes and res are regular JS numbers, primes has "-1" as first element to handle negative numbers.
// opt_maxnum is a limit for the amount of results you want, in case it would find too many results.
// Returns 3 things: array of b's, array of a's, array of factors of b's, where each b = a*a-n
Jmat.BigInt.qs_sieve_ = function(n, start, end, primes, res, opt_maxnum) {
  var B = Jmat.BigInt;
  var num = end.sub(start).toInt();
  var array = []; // array of approximate logarithms
  for (var i = 0; i < num; i++) array[i] = 0;

  // returns null if invalid, the factorization array if valid, which is used to avoid having to calculate a large sqrt later
  var trialdiv = function(b) {
    var result = [];
    if(b.minus) {
      result.push(-1);
      b = b.abs();
    }
    for(var i = 1; i < primes.length; i++) {
      if(b.eqr(0)) break;
      while(b.modr(primes[i]).eqr(0)) {
        b = b.divr(primes[i]);
        result.push(primes[i]);
      }
    }
    return b.eqr(1) ? result : null;
  };

  for (var i = 1; i < primes.length; i++) {
    var p = primes[i];
    if(p < 4) continue; // skip small primes (and also -1), compensate for it in the approximate log later...
    var lp = Math.log(p);
    var pos = -start.modr(p).toInt();
    while (pos < num) {
      var index0 = pos + res[i];
      var index1 = pos + (p - res[i]);
      if(index0 >= 0 && index0 < num) array[index0] += lp;
      if(res[i] != 0 && index1 >= 0 && index1 < num) array[index1] += lp;
      pos += p;
    }
  }

  var result = [[],[],[]];
  for(var i = 0; i < array.length; i++) {
    var a = start.addr(i);
    var b = a.mul(a).sub(n);
    var l = B.rlog(b.abs());
    if(array[i] >= 0.95 * l - 10) {  // Tweak approximate log comparison. TODO: finetune this
      var factors = trialdiv(b);
      if(factors) {
        result[0].push(b);
        result[1].push(a);
        result[2].push(factors);
        if(result[0].length >= opt_maxnum) break;
      }
    }
  }
  return result;
};


// Binary matrix rref for quadratic sieve
Jmat.BigInt.qs_rref_ = function(a) {
  var h = a.length;
  var w = a[0].length;

  var swaprow = function(matrix, y0, y1) {
    var temp = matrix[y0];
    matrix[y0] = matrix[y1];
    matrix[y1] = temp;
  };

  // xors row y1 with row y0 (modifying row y1), starting from x.
  var xorrow = function(matrix, x, y0, y1) {
    var w = matrix[0].length;
    for (var i = x; i < w; i++) {
      matrix[y1][i] ^= matrix[y0][i];
    }
  };

  var pivots = []; // x coordinate of pivot in each row (except the zero rows at the end, so may have smaller length than h)

  // gaussian elimination
  var k2 = 0; //next row, equal to k unless there are zero-rows
  for(var k = 0; k < w; k++) {
    var n = Jmat.Real.argmax(k2, h, function(i) { return a[i][k]; });
    if (a[n][k] == 0) continue; // singular, leave row as is
    if(k2 != n) swaprow(a, k2, n);
    for (var i = k2 + 1; i < h; i++) {
      if(a[i][k] != 0) xorrow(a, k, k2, i);
    }
    pivots.push(k);
    k2++;
    if(k2 >= h) break;
  }

  //now bring from row echolon form to reduced row echolon form
  for(var k = 0; k < pivots.length; k++) {
    var p = pivots[k];
    for(var y = k - 1; y >= 0; y--) {
      if(a[y][p] != 0) xorrow(a, p, k, y);
    }
  }
};

//Solves homogenous equation, with matrix given in rref form.
//Start with seed 1. For more solutions, increment seed. Once it returns null, no more new unseen solutions are available.
Jmat.BigInt.qs_solve_ = function(matrix, seed) {
  var h = matrix.length;
  var w = matrix[0].length;
  seed = seed || 1;

  var sol = [];
  for(var x = 0; x < w; x++) sol[x] = -1; //-1 means uninited

  var k = w; //how we're progressing to the left in the echelon matrix
  for(var y = h - 1; y >= 0; y--) {
    var f = 0; //first non-zero element in horizontal direction
    while(f < w && matrix[y][f] == 0) f++;
    if(f == w) continue; //row with all zeroes
    k = f;
    var odd = 0;
    for(var x = f + 1; x < w; x++) {
      if(matrix[y][x]) {
        if(sol[x] == -1) {
          sol[x] = (seed & 1); // can choose between 0 and 1, but note that always choosing 0 makes it end up with the trivial solution
          seed >>= 1;
        }
        if(sol[x] == 1) odd ^= 1;
      }
    }
    sol[f] = odd;
  }
  if(seed != 0) return null; // indicates there are no more unseen solutions: not all seed bits were used up

  var result = [];
  for(var i = 0; i < sol.length; i++) if(sol[i] > 0) result.push(i); //values that are -1 also count as 0
  return result;
};

// Quadratic sieve. Many TODO's, not best implementation yet, supports 30 digit numbers (two 15-digit prime factors).
Jmat.BigInt.qsieve_ = function(n, opt_verbose) {
  // Test in console: ''+Jmat.BigInt.qsieve_(BigInt.randomPrime(50).mul(BigInt.randomPrime(50)), true)
  var B = Jmat.BigInt;

  var pad = function(n, len) { len = len || 2; n = '' + n; while(n.length < len) n = '0' + n; return n; };
  var formattime = function(d) { return pad(d.getHours()) + ':' + pad(d.getMinutes()) + ':' + pad(d.getSeconds()) + '.' + pad(d.getMilliseconds(), 3); };
  var logwithtime = function(s) { if(opt_verbose) console.log('[' + formattime(new Date()) + ']: ' + s); };

  var l = B.rlog10(n);
  var sb = Math.floor(50 * l * l); // smoothness bound. TODO: tweak this better
  var s = B.sqrt(n);

  var primes = [-1].concat(Real.eratosthenes(sb));
  if (primes.length > 6000) return null; // Limit the input size this way, to avoid JS hanging too long. Without this return, it can factor 38!+1 in 12 minutes.

  var res = [0]; // residues
  var filtered = [-1];
  for(var i = 1; i < primes.length; i++) {
    var l = B.legendre(n, B(primes[i]));
    if(l == 0) return B(primes[i]); // found a small factor early
    if(l == 1) {
      filtered.push(primes[i]);
      res.push(B.ressol(n, B(primes[i])).toInt());
      if(filtered.length > 3000) break; // JS may crash when making 10k x 10k array. Limit to 3k x 3k.
    }
  }
  var num = filtered.length + 1; // the matrix must have at least one more smooth value than the amount of primes to algebraically guarantee a solution
  var step = num * 999; // TODO: tweak this amount
  logwithtime('input: ' + n + ', smoothness bound: ' + sb + ', num primes: ' + primes.length + ', num filtered primes: ' + filtered.length + ', sieving for ' + num + ' smooth numbers in batches of ' + step);
  primes = filtered;

  var pos = s.subr(step); // begin with a few a's under the sqrt (gives negative b's, handled with "-1" as factor): slightly more chance of smooth numbers there.
  var all_b = [], all_a = [], all_factored = [];
  // TODO: support multiple polynomials
  while (all_b.length < num) {
    var r = B.qs_sieve_(n, pos, pos.addr(step), primes, res, num - all_b.length);
    pos = pos.addr(step);
    all_b = all_b.concat(r[0]);
    all_a = all_a.concat(r[1]);
    all_factored = all_factored.concat(r[2]);
    logwithtime('found ' + r[0].length + ' for this batch, ' + all_b.length + ' total');
  }
  logwithtime('sieving done, found ' + all_b.length + ' smooth numbers');

  var reverseprimes = {};
  for(var i = 0; i < primes.length; i++) reverseprimes[primes[i]] = i;

  var matrix = []; // each column represents a factorization
  for(var j = 0; j < primes.length; j++) matrix.push([]);
  for(var i = 0; i < all_b.length; i++) {
    for(var j = 0; j < primes.length; j++) matrix[j][i] = 0;
    var f = all_factored[i];
    for(j = 0; j < f.length; j++) {
      matrix[reverseprimes[f[j]]][i] ^= 1; // xor instead of increment: binary matrix to find squares (even amount of factors)
    }
  }
  logwithtime('performing gaussian elemination...');
  B.qs_rref_(matrix);
  logwithtime('gaussian elemination done. Matrix size: ' + matrix.length + 'x' + matrix[0].length);

  for(var seed = 0; seed < 20; seed++) {
    logwithtime('solving...');
    var sol = B.qs_solve_(matrix, seed);
    if(!sol || !sol.length) break;

    var as = [], bs = []; // The a's and the b's, with b = a*a-n
    var factors = [];
    for(var i = 0; i < sol.length; i++) {
      bs.push(all_b[sol[i]]);
      as.push(all_a[sol[i]]);
      var f = all_factored[sol[i]];
      var fcount = [];
      for(var j = 0; j < primes.length; j++) fcount[j] = 0;
      for(var j = 0; j < f.length; j++) fcount[reverseprimes[f[j]]]++;
      for(var j = 0; j < fcount.length; j++) { factors[j] = (factors[j] || 0) + fcount[j]; } // sum fcount to factors
    }

    var mula = B.mulmany_(as);
    for(var i = 0; i < factors.length; i++) factors[i] >>= 1;
    var pfactors = [];
    for(var i = 0; i < factors.length; i++) {
      if(factors[i] != 0) pfactors.push(B.pow(B(primes[i]), B(factors[i])));
    }
    var s = B.mulmany_(pfactors); // s is now B.sqrt(mulb) with mulb = B.mulmany_(bs), but, that would be very slow, can be 100k bit number, so calculated from half of factors instead

    var gcd = B.gcd(n, mula.sub(s));
    logwithtime('gcd: ' + gcd);
    if(!gcd.eqr(1) && !gcd.eq(n)) return gcd; // found one!
  }
  return null;
};

// Returns prime factors in array.
// Will not factorize if the problem is too difficult, last element of result is 0 then to indicate the error.
// Result may be probabilistic since a probabilistic prime test is used (but such error is very unlikely).
// For inputs smaller than 2, includes the non-composite factors -1, 0 and 1 in the result
Jmat.BigInt.factorize = function(a) {
  var B = Jmat.BigInt;

  var result = [];
  if(a.minus) {
    a = B.neg(a);
    result.push(B(-1));
  }
  if(a.lter(2)) {
    if(result.length == 0 || !a.eqr(1)) result.push(a); // return [0] if x is 0, [1] if x is 1. Also avoids infinite loops.
    return result;
  }

  if (a.ltr(Jmat.Real.BIGGESTJSINT)) {
    var r = Jmat.Real.factorize(a);
    for(var i = 0; i < r.length; i++) result.push(B(r[i]));
    return result;
  }

  // returns a factor, or a itself if end reached (a is prime), or 0 if undetermined because the problem is too hard
  var f = function(a) {
    // Check 0: weed out low factors asap without expensive tests
    if (!(a.a[a.a.length - 1] & 1)) return B(2);
    if(a.modr(3).eqr(0)) return B(3);
    if(a.modr(5).eqr(0)) return B(5);
    if(a.modr(7).eqr(0)) return B(7);
    if(a.modr(11).eqr(0)) return B(11);
    if(a.modr(13).eqr(0)) return B(13);

    // Check 1: if it's prime, return self and stop.
    // Note that this may be probabilistic, so to be super sure, if the result has a very large prime factor, factorize a few more times.
    if(B.isPrime(a)) {
      return a;
    }

    // Check 2: simple trial division with primes.
    var num = Math.min(50000, B.sqrt(a).toInt());  // not too high, this is slow, the quadratic sieve will handle the rest
    var p = 1; //prime throughout the for loop
    var i = 0;
    for(;;) {
      // TODO: use a prime sieve
      if(i >= B.primeCache_.length) B.primeCache_[i] = Jmat.Real.nextPrime(p);
      p = B.primeCache_[i];
      i++;
      if(p <= 13) continue; // already tested above with the first few low factors
      if(p != p || p > num) break;
      if(a.modr(p).eqr(0)) {
        return B(p);
      }
    }

    // Check 3: if perfect power, its base is a factor
    var pw = B.perfectpow(a);
    if(pw[0].eq(a) && pw[1].lt(a)) {
      return f(pw[1]); // the perfect power base itself is not necessarily prime, hence recursive call
    }

    // Check 4: quadratic sieve
    var qs = B.qsieve_(a);
    if (qs) {
      return f(qs); // it's a factor, but not necessarily a prime one, so recurse
    }
    return B(0); // Not found
  }

  for(;;) {
    var b = f(a);
    if(b.eqr(0)) {
      // the quadratic sieve doesn't guarantee getting the smallest factor so must sort
      result.sort(function(a, b) { return B.compare(a, b); });
      result.push(B(0));
      return result;
    }
    result.push(b);
    if(b.eqr(a)) {
      // the quadratic sieve doesn't guarantee getting the smallest factor so must sort
      result.sort(function(a, b) { return B.compare(a, b); });
      return result;
    }
    a = a.div(b);
  }
};



// Returns legendre symbol for odd prime p. The legendre symbol is 0 if p
// divides a, 1 if a is quadratic residue mod p, p - 1 if it's not a quadratic
// residue (if p divides a it's technically also a quadratic residue).
Jmat.BigInt.legendre = function(a, p) {
  return Jmat.BigInt.modpow(a, p.subr(1).divr(2), p);
};


// Tonelli-Shanks algorithm aka "RESSOL": finds square root of n mod p, but p
// must be odd prime and n should be quadratic residue mod p (if it isn't, then
// it returns 0). That means it finds the x such that x*x = n (mod p).
// Returns one result r, but there are two solutions, and the other solution is
// always p - r (this all assumes legendre symbol is 1, when it's 0 then instead
// there is exactly one solution which is 0).
Jmat.BigInt.ressol = function(n, p) {
  var B = Jmat.BigInt;
  if(p.eqr(2)) return B(0);
  // TODO: have something precomputed that can do all these mod p's faster, and
  // similarly have opt_monred parameter for modpow.
  if (p.lter(2)) return null; // avoid infinite loops
  n = n.mod(p);
  var q = p.subr(1);
  var s = 0;
  while (B.bitandr(q, 1) == 0) {
    q = q.divr(2);
    s++;
  }
  if (s == 1) {
    var r = B.modpow(n, p.addr(1).divr(4), p);
    if (r.mul(r).mod(p).neq(n)) return B(0);
    return r;
  }

  var nr = B(1); // find a non-residue
  for(;;) {
    nr = nr.addr(1);
    if(B.legendre(nr, p).gtr(1)) break; // if legendre symbol is -1
  }

  var c = B.modpow(nr, q, p);
  var r = B.modpow(n, q.addr(1).divr(2), p);
  var t = B.modpow(n, q, p);
  var m = s;
  while (!t.eqr(1)) {
    var tmp = t;
    var i = 0;
    while(!tmp.eqr(1)) {
      tmp = tmp.mul(tmp).mod(p);
      i++;
      if(i == m) return 0;
    }
    var b = B.modpow(c, B.modpow(B(2), B(m - i - 1), p.subr(1)), p);
    tmp = b.mul(b).mod(p);
    r = r.mul(b).mod(p);
    t = t.mul(tmp).mod(p);
    c = tmp;
    m = i;
  }

  if (r.mul(r).mod(p).neq(n)) return B(0);
  return r;
};

// takes floor of base-y logarithm of x (y is also BigInt, but typically something like 2 or 10. logr is faster with immediately giving regular JS number.)
Jmat.BigInt.logy = function(x, y) {
  return Jmat.BigInt.logr(x, y.toInt());
};


// takes floor of base-y logarithm of x (with y a regular JS number, returns BigInt)
Jmat.BigInt.logr = function(x, y) {
  if(x.eqr(y)) return Jmat.BigInt(1);
  if(x.minus || y < 0) return undefined;
  if(x.eqr(0)) return undefined; //-Infinity
  if(y == 1) return x;
  if(y == 2) return Jmat.BigInt.log2(x);

  var c = Jmat.BigInt.convertBase(x, y);
  var result = Jmat.BigInt.fromInt(Jmat.BigInt.getNumDigits(c) - 1, x.radix);

  return result;
};

// integer log2 (floor)
Jmat.BigInt.log2 = function(x) {
  if(x.minus) return undefined;
  var bits = Jmat.BigInt.getNumBits(x);
  if(bits == 0) return undefined; //-Infinity
  return Jmat.BigInt.fromInt(bits - 1);
};

Jmat.BigInt.log10 = function(x) {
  return Jmat.BigInt.logr(x, 10);
};

// TODO: rename log2 to ilog2, and rlog2 to log2. Idem for log10, log, logy, ...
// returns log2 of BigInt x as a real value, including fractional part
Jmat.BigInt.rlog2 = function(x) {
  var bits = Jmat.BigInt.getNumBits(x);
  // Could actually work with up to 1022 bits: x.toInt() can return double precision float approximating the integer
  if(bits > 128) {
    x = Jmat.BigInt.rshift(x, bits - 128);
    return Jmat.Real.log2(x.toInt()) + bits - 128;
  }
  return Jmat.Real.log2(x.toInt());
};

// returns ln of BigInt x as a real value, including fractional part
Jmat.BigInt.rlog = function(x) {
  return Jmat.BigInt.rlog2(x) * Math.LN2;
};

// returns ln of BigInt x as a real value, including fractional part
Jmat.BigInt.rlog10 = function(x) {
  return Jmat.BigInt.rlog2(x) / 3.321928094887362 /*log2(10)*/;
};

/*
E.g. try:
Jmat.BigInt.div('1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139', '37975227936943673922808872755445627854565536638199');
Should give: '40094690950920881030683735292761468389214899724061'
*/
Jmat.BigInt.div = function(a, b) {
  return a.div(b);
};
Jmat.BigInt.prototype.div = function(b) {
  if(b.radix != this.radix) b = Jmat.BigInt.convertBase(b, this.radix);
  return Jmat.BigInt.divmod_(this, b)[0];
};

Jmat.BigInt.mod = function(a, b) {
  return a.mod(b);
};
Jmat.BigInt.prototype.mod = function(b) {
  if(b.radix != this.radix) b = Jmat.BigInt.convertBase(b, this.radix);
  return Jmat.BigInt.divmod_(this, b)[1];
};

//returns array with [quotient, mod]
Jmat.BigInt.divmod = function(a, b) {
  return a.divmod(b);
};
Jmat.BigInt.prototype.divmod = function(b) {
  if(b.radix != this.radix) b = Jmat.BigInt.convertBase(b, this.radix);
  if(b.abs().ltr(65536)) return Jmat.BigInt.divsmall_(this, b.toInt());
  return Jmat.BigInt.divmod_(this, b);
};

// Divide through regular (must be integer) js number
Jmat.BigInt.divr = function(a, b) {
  if(b == 0) return undefined;
  if(b == 1) return a;
  if(b == -1) return a.neg();

  if(a.abs().ltr(Math.abs(b))) return Jmat.BigInt(0);
  if(a.eqr(b)) return Jmat.BigInt(1);

  if(Math.abs(b) < 65536) return Jmat.BigInt.divsmall_(a, b)[0];
  return Jmat.BigInt.div(a, Jmat.BigInt(b));
};
Jmat.BigInt.prototype.divr = function(b) {
  return Jmat.BigInt.divr(this, b);
};

// Divides through two in linear time, in any radix (without doing any conversions).
// Returns array of a/2, a%2 (modulus)
Jmat.BigInt.half_ = function(a) {
  if(Jmat.Real.isPOT(a.radix)) {
    var odd = a.a[a.a.length - 1] % 2;
    return [Jmat.BigInt.rshift(a, 1), Jmat.BigInt(odd)];
  }

  var r = [];
  var radix = a.radix;
  var borrow = 0;
  for(var i = 0; i < a.a.length; i++) {
    var d = a.a[i] + borrow;
    borrow = Jmat.Real.isOdd(d) ? radix : 0;
    var d2 = Math.floor(d / 2);
    if(d2 > 0 || r.length > 0) r.push(d2);
  }
  return [new Jmat.BigInt(r, radix, a.minus), (borrow == 0) ? Jmat.BigInt(0, radix) : Jmat.BigInt(1, radix)];
};

// Divides through b in linear time, with b a regular js number (may be negative).
Jmat.BigInt.divsmall_ = function(a, b) {
  if(b == a.radix) return [a.rshift_radix(a), Jmat.BigInt(0)];
  if(b == 2) return Jmat.BigInt.half_(a);
  if(b == -2) {
    var result = Jmat.BigInt.half_(a);
    return [result[0].neg(), result[1].neg()];
  }

  var aminus = a.minus;
  var bminus = (b < 0);
  var minus = (a.minus != b < 0);
  b = Math.abs(b);

  var r = [];
  var radix = a.radix;
  var borrow = 0;
  for(var i = 0; i < a.a.length; i++) {
    var d = a.a[i] + borrow * radix;
    borrow = (d % b);
    var d2 = Math.floor(d / b);
    if(d2 > 0 || r.length > 0) r.push(d2);
  }

  var m = Jmat.BigInt(borrow, radix);

  //if(borrow != 0 && minus) m = m.addr(b); //but not if result was 0, then it should stay 0

  if(borrow) {
    if(aminus && bminus) m = m.neg();
    else if(bminus) m = m.subr(b);
    else if(aminus) m = m.subr(b).neg();
  }

  return [new Jmat.BigInt(r, radix, minus), m];
};

// Modulo with regular (must be integer) js number
// TODO: consider returning regular js number instead of BigInt as well, after all it fits in the range
Jmat.BigInt.modr = function(a, b) {
  if(b == 0) return undefined;
  if(b == 1 || b == -1) return Jmat.BigInt(0);

  if(Math.abs(b) < a.radix && a.radix % b == 0) {
    return Jmat.BigInt(Jmat.Real.mod(a.a[a.a.length - 1], b));
  }

  if(Math.abs(b) < 65536) return Jmat.BigInt.divsmall_(a, b)[1];
  return Jmat.BigInt.mod(a, Jmat.BigInt(b));
};
Jmat.BigInt.prototype.modr = function(b) {
  return Jmat.BigInt.modr(this, b);
};

/*
Division based on public domain Big Integer Library v. 5.0 by Leemon Baird, and
adapted to use the Jmat functions. The original used arrays mirrored compared to
Jmat and required leading zeroes, hence the elaborate array indices and spliced
zeroes. Both x and y must have the same radix, which must be power of two, with
bpe bits (max supported bpe = 15).
*/
Jmat.BigInt.leemondiv_ = function(x, y, bpe) {
  var B = Jmat.BigInt;

  var radix = x.radix;
  var mask = x.radix - 1;

  var kx, ky;
  var i, y1, y2, c, a, b;
  var r = B.copy(x);

  x = x.strip();
  y = B.copystrip_(y);

  ky = y.a.length;

  //normalize: ensure the most significant element of y has its highest bit set
  b = y.a[y.a.length - 1 - ky + 1];
  for(a = 0; b; a++) b >>= 1;
  a = bpe - a; //a is how many bits to shift so that the high order bit of y is leftmost in its array element

  y = y.lshift(a);
  r = r.lshift(a);
  y.a.splice(0, 0, 0);

  kx = Math.max(ky, r.a.length);

  var q = B(0, radix);
  for(var i = 0; i < x.a.length; i++) q.a[i] = 0;
  q.a.splice(0, 0, 0);
  for(;;) {
    var s = y.lshift_radix(kx - ky);
    if(s.gt(r)) break;
    r = r.sub(s);
    q.a[q.a.length - kx + ky - 1]++;
  }
  // the algorithm requires that r has the same length as x, prepend zeroes if needed
  while(r.a.length < x.a.length + 1) r.a.splice(0, 0, 0);

  for (var i = kx - 1; i >= ky; i--) {
    var rl = r.a.length - 1;
    var ql = q.a.length - 1;
    var yl = y.a.length - 1;

    if (r.a[rl - i] == y.a[yl - (ky - 1)]) {
      q.a[ql - (i - ky)] = mask;
    } else {
      q.a[ql - (i - ky)] = Math.floor((r.a[rl - i] * radix + r.a[rl - (i - 1)]) / y.a[yl - (ky - 1)]);
    }

    //The following for(;;) loop is equivalent to the commented while loop,
    //except that the uncommented version avoids overflow.
    //The commented loop comes from HAC, which assumes r[-1]==y[-1]==0
    //  while(q[ql-(i-ky)]*(y[yl-(ky-1)]*radix+y[yl-(ky-2)]) > r[rl-i]*radix*radix+r[rl-(i-1)]*radix+r[rl-(i-2)])
    //    q[ql-(i-ky)]--;
    for(;;) {
      y2 = (ky > 1 ? y.a[yl - (ky - 2)] : 0) * q.a[ql - (i - ky)];
      c = y2 >> bpe;
      y2 = y2 & mask;
      y1 = c + q.a[ql - (i - ky)] * y.a[yl - (ky - 1)];
      c = y1 >> bpe;
      y1 = y1 & mask;

      if(c == r.a[rl - i] ? y1 == r.a[rl - (i - 1)] ? y2 > (i > 1 ? r.a[rl - (i - 2)] : 0) : y1 > r.a[rl - (i - 1)] : c > r.a[rl - i]) {
        q.a[ql - (i - ky)]--;
      } else {
        break;
      }
    }

    var ys = y.lshift_radix(i - ky);
    var rs = ys.mulr(q.a[ql - (i - ky)]);
    if(r.lt(rs)) {
      r.a = Jmat.BigInt.baseloop_(r.a, 0, 1, ys.a, 0, 1, 0, radix, true);
      r.a = Jmat.BigInt.baseloop_(r.a, 0, 1, rs.a, 0, -1, 0, radix, true);
      q.a[ql-(i-ky)]--;
    } else {
      r.a = Jmat.BigInt.baseloop_(r.a, 0, 1, rs.a, 0, -1, 0, radix, true);
    }
  }

  r = r.rshift(a);
  return [q, r];
};

// Divides a and b and also returns remainder, but there are some requirements:
// They must already of type BigInt.
// If the optional bits per element, then radix of a and b are assumed to already be 2^(bpe), and no base conversions will be done.
Jmat.BigInt.divmod_ = function(a, b, opt_bpe) {
  var B = Jmat.BigInt;
  if(b.eqr(0)) return undefined;
  if(b.eqr(1)) return [a, B(0)];
  if(b.eqr(-1)) return [a.neg(), B(0)];
  if(b.abs().ltr(65536)) return Jmat.BigInt.divsmall_(a, b.toInt());

  if(b.abs().gt(a.abs())) return [B(0), a];
  if(b.eq(a)) return [B(1), B(0)];

  var minus = (a.minus != b.minus);
  a = a.abs();
  b = b.abs();

  if(!opt_bpe) {
    a = B.cast(a, B.ARRAYBASE_);
    b = B.cast(b, B.ARRAYBASE_);
  }

  var ar = B.leemondiv_(a, b, opt_bpe || B.ARRAYBASE_BITS_);

  var result = ar[0];
  var m = ar[1];
  B.stripInPlace_(result.a);
  B.stripInPlace_(m.a);

  // To test correctness: m must be equal to a.sub(result.mul(b)), and in range [0, b)
  // if(m.neq(a.sub(result.mul(b))) || m.ltr(0) || m.gte(b)) console.log('alert1: ' + a + ' ' + b + ' ' + result + ' ' + m);

  if(minus) {
    result = result.neg();
    m = m.add(b);
    if(m.minus && m.eqr(0)) m.minus = false;
  }
  return [result, m];
};

Jmat.BigInt.pow = function(a, b) {
  var origb = b;
  b = Jmat.BigInt.cast(b, 2);
  var radix = a.radix;

  if(b.minus) return Jmat.BigInt(0, radix);
  var minus = false;
  if(a.minus && b.modr(2).eqr(1)) minus = true;
  a = a.abs();

  // Montgomery's ladder
  var ba = Jmat.BigInt.maybecopystrip_(b.a, b != origb);
  var a1 = Jmat.BigInt(1, radix);
  var a2 = a;
  var l = ba.length;
  for(var i = 0; i < l; i++) {
    if(ba[i] == 0) {
      a2 = a1.mul(a2);
      a1 = a1.mul(a1);
    } else {
      a1 = a1.mul(a2);
      a2 = a2.mul(a2);
    }
  }
  if(minus) a1 = a1.neg();
  return a1;
};
Jmat.BigInt.prototype.pow = function(b) {
  return Jmat.BigInt.pow(this, b);
};

Jmat.BigInt.powr = function(a, b) {
  return Jmat.BigInt.pow(a, Jmat.BigInt.fromInt(b));
};
Jmat.BigInt.prototype.powr = function(b) {
  return Jmat.BigInt.pow(this, Jmat.BigInt.fromInt(b));
};

//greatest common divisor
Jmat.BigInt.gcd = function(x, y) {
  var B = Jmat.BigInt;
  x = x.abs();
  y = y.abs();
 //Euclid's algorithm
 for(;;) {
   if(y.eqr(0)) return x;
   var z = B.mod(x, y);
   x = y;
   y = z;
 }
};


// Extended Euclidean algorithm. Returns array of 5 values: [gcd, bezout coeff. x, bezout coeff. y, x / gcd, y / gcd]
Jmat.BigInt.egcd = function(x, y) {
  var B = Jmat.BigInt;
  var s = B(0);
  var olds = B(1);
  var t = B(1);
  var oldt = B(0);
  var r = x;
  var oldr = y;
  var temp;
  while(!r.eqr(0)) {
    var q = oldr.div(r);
    temp = r;
    r = oldr.sub(q.mul(r));
    oldr = temp;
    temp = s;
    s = olds.sub(q.mul(s));
    olds = temp;
    temp = t;
    t = oldt.sub(q.mul(t));
    oldt = temp;
  }
  if(t.sign() != x.sign() || s.sign() != y.sign()) {
    temp = t;
    t = s;
    s = temp;
    temp = oldt;
    oldt = olds;
    olds = temp;
    if(t.sign() != x.sign()) t = t.neg();
    if(s.sign() != y.sign()) s = s.neg();
  }
  // [gcd, bezout coeff. x, bezout coeff. y, x / gcd, y / gcd]
  return [oldr, olds, oldt, t, s];
};

// calculates a^-1 mod m, modular inverse, so that (a*result) mod m == 1
// result only exists if gcd(a, m) == 1 (they're coprime), and m > 2
Jmat.BigInt.invmod = function(a, m) {
  var B = Jmat.BigInt;
  var origm = m;

  var x = B(1);
  var y = B(0);
  var result;

  for(;;) {
    if(a.eqr(1)) { result = x; break; }
    if(a.eqr(0)) { result = B(0); break; }
    var d = m.divmod(a);
    y = y.sub(x.mul(d[0]));
    m = d[1];

    if(m.eqr(1)) { result = y; break; }
    if(m.eqr(0)) { result = B(0); break; }
    d = a.divmod(m);
    x = x.sub(y.mul(d[0]));
    a = d[1];
  }

  if(result.minus) result = origm.add(result);
  return result;
};

//montgomery reduction (aka REDC): calculates a/r mod m (the division of course in modulo, so doesn't necessarily make a smaller)
//a must be smaller than m*r-1
//r must be power of 2 and bits is its log2, mask is r-1 (all ones)
//mi must be -(m^(-1)) mod r (there's a minus sign in front, but since it's mod r it's positive again)
Jmat.BigInt.monred_ = function(a, bits, mask, m, mi) {
  var s = a.mul(mi).bitand(mask); // Could also apply mask to a and mi, but usually mask is bigger so not a speedup.
  var t = a.add(s.mul(m)).rshift(bits);
  if(t.gte(m)) {
    t = t.sub(m);
    if(t.gte(m)) return null; // error, probably a was larger than m*r-1
  }

  return t;
};

//generates a function that can do the montgomery reduction for some value for modulo m, that has all precomputed values bound in it (including r). m must be odd.
Jmat.BigInt.genmonred_ = function(m) {
  var bits = Jmat.BigInt.getNumBits(m);
  var r = Jmat.BigInt.lshift(Jmat.BigInt(1), bits);
  var mask = r.subr(1);
  var mi = Jmat.BigInt.invmod(m, r).neg().bitand(mask);
  var rrm = r.lshift(bits).mod(m);

  // without init, does montgomery reduction. With init, does conversion towards montgomery domain.
  // E.g. to multiply a and b: var am = monred(a, true); var bm = monred(b, true); var cm = monred(am.mul(bm)); var c = monred(cm);
  return function(a, init) {
    if(init) {
      //The below is equivalent to: a.lshift(bits).mod(m), but faster (no mod m)
      return Jmat.BigInt.monred_(a.mul(rrm), bits, mask, m, mi); //convert to montgomery domain
    } else {
      return Jmat.BigInt.monred_(a, bits, mask, m, mi); //montgomery reduction
    }
  };
};

//modular exponentiation: (a^b) mod m
//faster if m is odd, e.g. try Jmat.BigInt.modpow(3, '500000000000000000000000', 2001)
//opt_monred is an optional parameter to have a precalculated montgomery reduction object for m with genmonred_. This may be faster if reused multiple times for the same m.
Jmat.BigInt.modpow = function(a, b, m, opt_monred) {
  var B = Jmat.BigInt;
  var origb = b;
  b = B.cast(b, 2);
  var ba = B.maybecopystrip_(b.a, b != origb);
  var a1 = B.ONE;
  var a2 = a;
  if(a2.gt(m)) a2 = a2.mod(m); // otherwise the monred code doesn't work due to too large base
  var l = ba.length;

  if(m.modr(2).eqr(1) && l > 1) {
    // Odd m, so faster Montgomery reduction with power of two r possible

    var monredm = opt_monred || B.genmonred_(m);

    a1 = monredm(a1, true);
    a2 = monredm(a2, true);

    // Montgomery's ladder
    for(var i = 0; i < l; i++) {
      if(ba[i] == 0) {
        a2 = a1.mul(a2);
        a1 = a1.mul(a1);
      } else {
        a1 = a1.mul(a2);
        a2 = a2.mul(a2);
      }
      a1 = monredm(a1);
      a2 = monredm(a2);
    }

    a1 = monredm(a1);
    return a1;
  } else {
    // Even m, slower
    // Montgomery's ladder
    for(var i = 0; i < l; i++) {
      if(ba[i] == 0) {
        a2 = a1.mul(a2).mod(m);
        a1 = a1.mul(a1).mod(m);
      } else {
        a1 = a1.mul(a2).mod(m);
        a2 = a2.mul(a2).mod(m);
      }
    }
    return a1;
  }
};

// Get random number with that amount of bits
Jmat.BigInt.randomBits = function(bits) {
  var result = Jmat.BigInt([], Jmat.BigInt.ARRAYBASE_);
  var numbytes = Math.ceil(bits / Jmat.BigInt.ARRAYBASE_BITS_);
  for(var i = 0; i < numbytes; i++) {
    result.a[i] = Math.floor(Jmat.BigInt.ARRAYBASE_ * Math.random());
  }
  //bits that should be zeroed because the last byte has 8 bits but bits may be non multiple of 8
  var nonbits = (bits * (Jmat.BigInt.ARRAYBASE_BITS_ - 1)) % Jmat.BigInt.ARRAYBASE_BITS_;
  var mask = Jmat.BigInt.ARRAYBASE_ - 1;
  for(var i = 0; i < nonbits; i++) {
    mask = Math.floor(mask / 2);
  }
  if(mask + 1 < Jmat.BigInt.ARRAYBASE_) result.a[0] &= mask;
  return result;
};

// Returns a random probable prime with roughly the given amount of bits.
// E.g. BigInt.randomPrime(20) could return something like 156007 or 787477
Jmat.BigInt.randomPrime = function(bits) {
  return Jmat.BigInt.nextPrime(Jmat.BigInt.randomBits(bits));
};

Jmat.BigInt.firstPrimes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

// Tests for being one of the first primes or being divisible through one of them.
// Returns 0 if not prime, 1 if prime, -1 if unknown by this function
Jmat.BigInt.isPrimeSimple = function(n) {
  n = Jmat.BigInt.cast(n);
  if(n.ltr(2)) return 0;
  for(var i = 0; i < Jmat.BigInt.firstPrimes_.length; i++) {
    if(n.eqr(Jmat.BigInt.firstPrimes_[i])) return 1;
    if(n.modr(Jmat.BigInt.firstPrimes_[i]).toInt() == 0) return 0;
  }
  return -1;
};

//do rounds of Miller-Rabin primality test
//opt_base = the potential witnesses to try as an array, e.g. [2, 3]. The members of the array may be either regular JS number or BigInt. Optional parameter, if not given, chosen automatically (random)
//requires that n is big enough, at least 3. First use Jmat.BigInt.isPrimeSimple, and only use this function if that returns -1.
Jmat.BigInt.isPrimeMillerRabin = function(n, opt_base) {
  var B = Jmat.BigInt;
  var base = opt_base || B.chooseMillerRabinBase_(n);

  // choose s and odd d such that n = 2^s * d
  var d = n.divr(2);
  var s = B.ONE;
  while(B.bitand(d, B.ONE).eqr(0)) {
    d = d.divr(2);
    s = s.addr(1);
  }

  var monred = B.genmonred_(n);

  var witness = function(n, s, d, a) {
    var x = B.modpow(a, d, n, monred);
    var y;
    while(!s.eqr(0)) {
      //this mod could also use monred (use vars ym and xm in mon domain). However, it doesn't speed up much, so code kept simpler here.
      y = x.mul(x).mod(n);
      if(y.eqr(1) && !x.eqr(1) && !x.eq(n.subr(1))) return false;
      x = y;
      s = s.subr(1);
    }
    return y.eqr(1);
  };

  for(var i = 0; i < base.length; i++) {
    if(!witness(n, s, d, B.cast(base[i]))) return false; //proven to be composite by this witness
  }
  return true; //probably prime, at least no compositeness was proven with the given witnesses
};

//Choose witnesses for Miller-Rabin primality test, with reasonable defaults (deterministic if n small enough, random otherwise so not reproducible).
Jmat.BigInt.chooseMillerRabinBase_ = function(n) {
  var bits = Jmat.BigInt.getNumBits(n);

  var base;
  if(bits <= 64)  {
    if(n.ltr(1373653)) base = [2, 3];
    else if(n.ltr(9080191)) base = [31, 73];
    else if(n.ltr(4759123141)) base = [2, 7, 61];
    else if(n.ltr(1122004669633)) base = [2, 13, 23, 1662803];
    else if(n.ltr(2152302898747)) base = [2, 3, 5, 7, 11];
    else if(n.ltr(3474749660383)) base = [2, 3, 5, 7, 11, 13];
    else if(n.ltr(341550071728321)) base = [2, 3, 5, 7, 11, 13, 17];
    else if(n.ltr(3770579582154547)) base = [2, 2570940, 880937, 610386380, 4130785767];
    else base = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]; //valid up to >2^64
  } else {
    base = [2, 3];
    for(var i = 0; i < 18; i++) {
      base.push(Jmat.BigInt.randomBits(bits - 1));
    }
  }

  return base;
};

//Is *probably* prime at least.
//If you wish to control amount of miller rabin rounds and bases, use isPrimeMillerRabin directly with own bases (first use isPrimeSimple).
Jmat.BigInt.isPrime = function(n) {
  if(n.ltr(Jmat.Real.BIGGESTJSPRIME)) return Jmat.Real.isPrime(n.toInt());
  var init = Jmat.BigInt.isPrimeSimple(n);
  if(init != -1) return !!init;

  return Jmat.BigInt.isPrimeMillerRabin(n);
};

/*
little benchmark on a prime:

var ta0 = new Date().getTime(); var p = BigInt.isPrime('40094690950920881030683735292761468389214899724061'); var ta1 = new Date().getTime();
console.log(((ta1 - ta0) / 1000.0) + ' ' + p);

Numbers to test that are not prime (each one has 2 large prime factors):
124620366781718784065835044608106590434820374651678805754818788883289666801188210855036039570272508747509864768438458621054865537970253930571891217684318286362846948405301614416430468066875699415246993185704183030512549594371372159029236099
1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139
88357 (149*593), strong pseudoprime
*/

Jmat.BigInt.d_ = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
                  'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
                  'W', 'X', 'Y', 'Z'];

Jmat.BigInt.di_ = function(c) {
  var i = c.charCodeAt(0);
  if(i <= 57) return i - 48;
  if(i > 90) i -= 32;
  return i - 65 + 10;
};

//reverses the array in-place (unlike JS array.reverse())
Jmat.BigInt.mirror_ = function(array) {
  for(var i = 0; i < array.length / 2; i++) {
    var temp = array[i];
    array[i] = array[array.length - i - 1];
    array[array.length - i - 1] = temp;
  }
};

// Similar to log2, but result is one higher, ignores sign, supports 0, and result is regular JS number instead of BigInt
Jmat.BigInt.getNumBits = function(x) {
  var pot = x.radix && ((x.radix & (x.radix - 1)) == 0);
  if(!pot) x = Jmat.BigInt.cast(x, Jmat.BigInt.ARRAYBASE_);

  var ar = Jmat.BigInt.maybecopystrip_(x.a);
  if(ar.length == 0) return 0;
  var result = Jmat.Real.ilog2(ar[0]) + 1;
  if(ar.length > 0) {
    var n = Jmat.Real.ilog2(x.radix); //x.radix is power of two as enforced above
    result += n * (ar.length - 1);
  }
  return result;
};

//gets number of significant digits for the radix x is in (does not include leading zeroes, and returns zero for zero)
Jmat.BigInt.getNumDigits = function(x) {
  var i = 0;
  while(i < x.a.length && x.a[i] == 0) i++;
  return x.a.length - i;
};

//strips array a (removing leading zeroes unless it is zero), modifying the object
//NOTE: do not give a BigInt, a must be an array
Jmat.BigInt.stripInPlace_ = function(a) {
  while(a[0] == 0 && a.length > 1) a.shift(); //JS array shift function (a is array, not BigInt)
};

//makes a copy of BigInt a, and stripped of leading zeroes if any (unless it is zero)
Jmat.BigInt.copystrip_ = function(a) {
  var numzeroes = 0;
  for(var i = 0; i < a.a.length; i++) {
    if(a.a[i] != 0) break;
    numzeroes++;
  }
  var result = new Array(a.a.length - numzeroes);
  if(result.length == 0) return new Jmat.BigInt([0], a.radix);
  for(var i = 0; i < result.length; i++) result[i] = a.a[i + numzeroes];
  return Jmat.BigInt(result, a.radix);
};

//makes copy if needed object needs to be altered, or original if it doesn't have to be stripped
//NOTE: do not give a BigInt, a must be an array
Jmat.BigInt.maybecopystrip_ = function(a, allowmodify) {
  if(a.length <= 1 || a[0] != 0) return a;
  if(!allowmodify) a = Jmat.BigInt.copyarray_(a);
  Jmat.BigInt.stripInPlace_(a);
  return a;
};

//strips (trims) BigInt a of leading zeroes. May return input. Does not make copy if not needed. Does not modify input.
Jmat.BigInt.strip = function(a) {
  return a.strip();
};
Jmat.BigInt.prototype.strip = function() {
  if(this.a.length <= 1 || this.a[0] != 0) return this;
  return Jmat.BigInt.copystrip_(this);
};

Jmat.BigInt.copyarray_ = function(a) {
  return a.slice(0);
};

////////////////////////////////////////////////////////////////////////////////

//BigInt formats: turn BigInt functions into functions that also take string and regular JS numbers as input

// "int" (plain JS number) is not a format on purpose: may not be able to contain the result.
Jmat.BigInt.FORMAT_UNKNOWN_ = 0;
Jmat.BigInt.FORMAT_BIGINT_ = 1;
Jmat.BigInt.FORMAT_ARRAY_ = 2;
Jmat.BigInt.FORMAT_STRING_ = 3;

Jmat.BigInt.getFormat = function(v) {
  if(typeof v == 'string') return Jmat.BigInt.FORMAT_STRING_;
  if(Object.prototype.toString.call(v) === "[object Array]") return Jmat.BigInt.FORMAT_ARRAY_;
  if(v instanceof Jmat.BigInt) return Jmat.BigInt.FORMAT_BIGINT_;
  return Jmat.BigInt.FORMAT_UNKNOWN_;
};

//v may be any also supported format
//if opt_minus is given, overrides sign if the output format supports it
Jmat.BigInt.toFormat = function(v, format, opt_minus) {
  if(v == undefined) return undefined; //propagate error condition
  if(format == Jmat.BigInt.getFormat(v)) {
    if(opt_minus == undefined) return v;
    if(format == Jmat.BigInt.FORMAT_BIGINT_) return (v.minus == opt_minus) ? v : v.neg();
    if(format == Jmat.BigInt.FORMAT_STRING_) return ((v[0] == '-') == opt_minus) ? v : (opt_minus ? ('-' + v) : v.substr(1));
  }
  v = Jmat.BigInt.cast(v);
  if(opt_minus != undefined && v.minus != opt_minus) v = v.neg();
  if(format == Jmat.BigInt.FORMAT_BIGINT_) return v;
  if(format == Jmat.BigInt.FORMAT_ARRAY_) return v.a;
  if(format == Jmat.BigInt.FORMAT_STRING_) return v.toString();
  if(format == Jmat.BigInt.FORMAT_UNKNOWN_) return v.toString(); //unknown also as string, is typically regular JS number, string displays better in JS console
  return v;
};

//Add the ability to a member function to auto-convert parameters from int, string or array to BigInt, and back at the result, by replacing it
//fname: name of the function as member of Jmat.BigInt or Jmat.BigInt.prototype
//numb: number of BigInt arguments. Must be the first ones.
//type:
// 0: non-prototype, output not converted.
// 1: non-prototype, output is BigInt but converted back to first argument format
// 2: prototype, output not converted (kept as BigInt)
// 3: add both a non-prototype version, with non-converted output, and, a prototype version, with one numb less, and non-converted output
// 4: add both a non-prototype version, with converted output, and, a prototype version, with one numb less
// Prototype means that Jmat.BigInt.prototype is used instead of Jmat.BigInt
Jmat.BigInt.enrichFunction_ = function(object, fname, numb, type) {
  if(type == 3) {
    Jmat.BigInt.enrichFunction_(fname, numb, 0);
    if(numb > 1) Jmat.BigInt.enrichFunction_(fname, numb - 1, 2);
    return;
  }
  if(type == 4) {
    Jmat.BigInt.enrichFunction_(fname, numb, 1);
    if(numb > 1) Jmat.BigInt.enrichFunction_(fname, numb - 1, 2);
    return;
  }
  var f = ((type == 2) ? Jmat.BigInt.prototype[fname] : Jmat.BigInt[fname]);
  var f2 = function() {
    var format = Jmat.BigInt.getFormat(arguments[0]);
    for(var i = 0; i < numb; i++) {
      arguments[i] = Jmat.BigInt.cast(arguments[i]);
    }
    var res = f.apply(this, arguments);
    if(type == 1) return Jmat.BigInt.toFormat(res, format);
    return res;
  };
  if(f.length == 2 && numb == 2 && type == 1) {
    f2 = function(a, b) {
      var result = f(Jmat.BigInt.cast(a), Jmat.BigInt.cast(b));
      return Jmat.BigInt.toFormat(result, Jmat.BigInt.getFormat(a));
    };
  }
  if(f.length == 1 && numb == 1 && type == 2) {
    f2 = function(b) {
      return f.call(this, Jmat.BigInt.cast(b));
    };
  }
  if(type == 2) object.prototype[fname] = f2;
  else object[fname] = f2;
};

Jmat.BigInt.enrichFunctions_ = function() {
  // set up
  var prototoo = 0; //set to 3 to also replace prototype functions (slower), 0 otherwise
  var otherobject = true; //apply the enrichment to another object (Jmat.BigIntC, C from "converting") instead, to avoid penalizing performance and make debugging easier

  if(otherobject) {
    Jmat.BigIntC = function() {
      return Jmat.BigInt.apply(this, arguments);
    };
    for(var v in Jmat.BigInt) Jmat.BigIntC[v] = Jmat.BigInt[v];
    for(var v in Jmat.BigInt.prototype) Jmat.BigIntC.prototype[v] = Jmat.BigInt.prototype[v];
  }

  var object = otherobject ? Jmat.BigIntC : Jmat.BigInt;

  Jmat.BigInt.enrichFunction_(object, 'add', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'addr', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'sub', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'subr', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'mul', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'mulr', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'div', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'divr', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'mod', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'modr', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'divmod', 2, 0 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'lshift', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'rshift', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'bitand', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'bitor', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'bitxor', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'bitneg', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'bitnot', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'abs', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'neg', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'sign', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'getSign', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'nonZero', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'compare', 2, 0 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'comparer', 1, 0 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'sqrt', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'root', 2, 1);
  Jmat.BigInt.enrichFunction_(object, 'rootr', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'logy', 2, 1);
  Jmat.BigInt.enrichFunction_(object, 'logr', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'log2', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'log10', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'rlog2', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'rlog', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'rlog10', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'pow', 2, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'powr', 1, 1 + prototoo);
  Jmat.BigInt.enrichFunction_(object, 'gcd', 2, 1);
  Jmat.BigInt.enrichFunction_(object, 'egcd', 2, 0);
  Jmat.BigInt.enrichFunction_(object, 'invmod', 2, 1);
  Jmat.BigInt.enrichFunction_(object, 'modpow', 3, 1);
  Jmat.BigInt.enrichFunction_(object, 'isPrimeSimple', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'isPrimeMillerRabin', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'isPrime', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'perfectpow', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'nextPrime', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'factorize', 1, 0);
  Jmat.BigInt.enrichFunction_(object, 'factorial', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'primorial', 1, 1);
  Jmat.BigInt.enrichFunction_(object, 'legendre', 2, 0);
  Jmat.BigInt.enrichFunction_(object, 'ressol', 2, 1);
};

Jmat.BigInt.enrichFunctions_();

