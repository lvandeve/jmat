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

// REQUIRES: no dependencies, most basic file of jmat.js.

/*
Jmat.Real: real math operating on plain JS numbers. Similar to JS's Math library except with more functions and algorithms.

Aliased as simply "Real" by jmat.js - disable that if it causes name clashes

Overview of some functionality:
-most standard Math functions are also copied in here
-polyfill for functions of higher versions of JavaScript: Real.log2, Real.log10, Real.clz32, ...
-mod and remainder: Real.mod, Real.rem, Real.wrap, Real.clamp
-special functions. Real.gamma, Real.erf, Real.erfc, Real.lambertw, ...
-primes and factors: Real.isPrime, Real.eratosthenes, Real.factorize, Real.nextPrime, Real.previousPrime, Real.nearestPrime, Real.eulerTotient, Real.gcd, Real.lcm
-date and time: Real.isLeapYear, Real.dayOfWeek,
*/

/** @constructor
Namespace for all of Jmat. Defined in jmat_real.js as this is the first js file that everything else depends on.
*/
function Jmat() {
  // Empty, this is a namespace, no need to ever call this
}

/** @constructor
namespace for real functions
*/
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
Jmat.Real.BIGGESTJSINT = 9007199254740992; // largest number that JS (float64) can represent as integer: 2^53, 0x20000000000000, 9007199254740992
Jmat.Real.BIGGESTJSPRIME = 9007199254740881; // largest prime number that JS (float64) can represent as integer, that is, the biggest prime smaller than Jmat.Real.BIGGESTJSINT.

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

// x is power of two
Jmat.Real.isPOT = function(x) {
  return x != 0 && (x & (x - 1)) == 0;
};

//isnanorinf isinfornan
Jmat.Real.isInfOrNaN = function(x) {
  return x == Infinity || x == -Infinity || isNaN(x);
};

////////////////////////////////////////////////////////////////////////////////

// dist, cheb and manhattan all return regular real JS numbers for all types. In some types they are all the same, but not for e.g. Complex or Matrix.
// Euclidean distance
Jmat.Real.dist = function(a, b) {
  if(a == b) return 0; // this is to avoid subtracting Infinity - Infinity
  return Math.abs(a - b);
};
//Chebyshev distance
Jmat.Real.cheb = function(a, b) {
  return Jmat.Real.dist(a, b);
};
//Manhattan distance
Jmat.Real.manhattan = function(a, b) {
  return Jmat.Real.dist(a, b);
};

// Modulo operation. Different than JS's % operator in case of negative operands.
// Result has the sign of the divisor b.
// Works for non-integers too, similar to "fmod" in C (in case of positive arguments).
// Compare with rem: Different programs and programming languages use different
// names for this, there is no convention which one has which sign, though in
// languages with both a "mod" and "rem", the convention adopted here is most
// popular. See the table at https://en.wikipedia.org/wiki/Modulo_operation.
//
// mod in terms of rem (%). The table below compares the two operators.
// x    :   -4 -3 -2 -1  0  1  2  3  4
// x mod  3:    2  0  1  2  0  1  2  0  1
// x mod -3:   -1  0 -2 -1  0 -2 -1  0 -2
// x rem  3:   -1  0 -2 -1  0  1  2  0  1
// x rem -3:   -1  0 -2 -1  0  1  2  0  1
// The sign of mod is that of b, while that of rem is that of a.
//
// "mod" is the one that is mathematically more useful, while "rem" is the one
// matching the "%" operator in most programming languages.
// mod corresponds to floored division, while rem corresponds to truncated division.
Jmat.Real.mod = function(a, b) {
  return a - Math.floor(a / b) * b; // alternative: (b + (a % b)) % b
};

// Remainder. This is the same as the % operator.
// Result has the sign of the dividend a.
// Compare with Jmat.Real.mod, which is different and contains more description about the difference between rem and mod.
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

// floored integer division. Note that this is distinct from the truncated integer division used on many platforms.
Jmat.Real.idiv = function(a, b) {
  return Math.floor(a / b);
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

// checks whether a is a perfect power of b, e.g. 27 is a power of 3, but 10 is not a power of 5. Negative values are not supported.
// returns 0 if not power of (or the power is 0 - that is trivial if a is 1)
// returns the positive power otherwise
Jmat.Real.isPowerOf = function(a, b) {
  var R = Jmat.Real;
  if(a == b) return 1;
  if(b <= 0) return 0; // false
  if(a <= 0) return 0; // false
  if(a == 1) return 0; // true but it's 0
  if(b > a) return 0; // false
  if(R.isPOT(a) && R.isPOT(b)) {
    var la = R.ilog2(a);
    var lb = R.ilog2(b);
    if(la % lb == 0) return la / lb;
    return 0;
  }
  if(R.isPOT(a) != R.isPOT(b)) return 0; // false
  if(R.isEven(a) != R.isEven(b)) return 0; // false (or a is 1)
  var c = b;
  // Binary search with powers.
  var bs = [];
  var bb = b;
  var result = 1;
  while(c < a) {
    bs.push(bb);
    c *= bb;
    result *= 2;
    if(c == a) return result;
    bb = bb * bb;
  }
  if(c == Infinity) return 0;
  while(bs.length > 0) {
    var p = bs.pop();
    if(c > a) {
      c /= p;
      result -= (1 << bs.length);
    }
    else {
      c *= p;
      result += (1 << bs.length);
    }
    if(c == a) return result;
  }
  return 0;
};

Jmat.Real.firstPrimes_ = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

// Initial set up shared by several of the prime test functions.
// Returns 0 if not prime, 1 if prime, NaN if problem, -1 if unknown by this function
Jmat.Real.isPrimeInit_ = function(n) {
  if(n == Infinity || n != n) return NaN;
  if(n != Math.round(n)) return 0;
  if(n < 2) return 0;
  if(n > Jmat.Real.BIGGESTJSINT) return NaN; //too large for the floating point's integer precision, result will not make sense
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
//Not probabilistic, but relies on the generalized Riemann hypothesis
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

  // returns (a + b) % c, taking overflow into account (in JS, overflow means reaching a part in the floating point representation where it can no longer distinguish 1)
  var modadd = function(a, b, c) {
    if (a + b < Jmat.Real.BIGGESTJSINT) return (a + b) % c;
    if(a + b > c) {
      return (a - c + b) % c;
    }
    // This assumes that c < 4503599627370496 or a + b doesn't overflow
    return ((a % c) + (b % c)) % c;
  };

  // returns (a * b) % c, taking overflow into account
  var modmul = function(a, b, c) {
    if(a * b < Jmat.Real.BIGGESTJSINT) return (a * b) % c;
    var x = 0;
    var y = a % c;
    while(b > 0) {
      if(b & 1) x = modadd(x, y, c);
      y = modadd(y, y, c);
      b = Math.floor(b / 2);
    }
    return x % c;
  };

  // returns (a to the n) % mod, taking overflow into account
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
    return y == 1;
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
    if(a != b) console.log('error: ' + i + ' ' + a + ' ' + b);
  }
  console.log('ok: ' + n);
};
testfun(100000);

Nice primes to test:
3770579582154547 --> NOT prime, but above this boundary, last "base" for miller rabin test is used
9007199254740992 --> NOT prime, but highest integer number that JavaScript supports
9007199254740881: just small enough for JS! ==> overflow with sum, does not work
4444280714420857: largest for half JS bits
311111111111113: for third last base
344555666677777: for second last base
*/

//Returns 1 if prime, 0 if not prime, NaN if error.
Jmat.Real.isPrime = function(n) {
  // below that, the "slow" method is faster. For higher values, Miller Rabin becomes more and more significantly faster.
  return (n < 1500000) ? Jmat.Real.isPrimeSlow_(n) : Jmat.Real.isPrimeMillerRabin_(n);
};

// Sieve of Eratosthenes: returns array of all the primes up to n.
Jmat.Real.eratosthenes = function(n) {
  if(n < 2) return [];
  var result = [2];
  var a = [];
  var s = Math.floor(Math.sqrt(n));

  var num = Math.ceil(n / 2);
  // a[i] represents odd numbers: a[0] represents 1, a[1] represents 3, a[n] represents n*2 + 1, m is represented by a[floor(m / 2)]
  for(var i = 0; i < num; i++) a[i] = true;
  for(var m = 3; m <= s; m += 2) {
    var i = Math.floor(m / 2);
    if(!a[i]) continue;
    for(var j = i + m; j < num; j += m) a[j] = false;
  }

  for(var i = 1; i <= n; i++) {
    if(a[i]) result.push((i * 2) + 1);
  }
  return result;
};

//for factorize
Jmat.Real.smallestPrimeFactor = function(x) {
  if(x == Infinity || x != x) return NaN;
  if(x != Math.round(x)) return NaN;
  if(x < 1) return NaN;
  if(x > Jmat.Real.BIGGESTJSINT) return NaN; //too large for the floating point's integer precision, result will not make sense
  if(x == 1) return 1;
  for(var i = 0; i < Jmat.Real.firstPrimes_.length; i++) {
    if(x == Jmat.Real.firstPrimes_[i]) return x;
    if(x % Jmat.Real.firstPrimes_[i] == 0) return Jmat.Real.firstPrimes_[i];
  }
  var p = Jmat.Real.firstPrimes_[Jmat.Real.firstPrimes_.length - 1];
  var s = Math.ceil(Math.sqrt(x));
  p = Math.floor(p / 6) * 6;
  while(p < s + 5) {
    if(x % (p - 1) == 0) return p - 1;
    if(x % (p + 1) == 0) return p + 1;
    p += 6;
  }
  return x;
};

//factorize: returns prime factors as array of real integers, sorted from smallest to largest. x must be integer.
Jmat.Real.factorize = function(x) {
  if(x > Jmat.Real.BIGGESTJSINT) return undefined; //too large for the floating point's integer precision, will cause crash
  var x = Math.round(x);
  var result = [];
  if(x < 0) {
    x = -x;
    result.push(-1);
  }
  if(x <= 2) {
    if(result.length == 0 || x != 1) result.push(x); // return [0] if x is 0, [1] if x is 1
    return result;
  }
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
  // Anything below 6 does not work with the calculations below.
  if(x < 7) {
    if(x <= 2) return 2;
    if(x <= 4) return 3;
    return 5;
  }
  if(x == Infinity || x != x) return NaN;
  if(x >= 9007199254740881) return NaN; //largest supported prime in floating point precision, after this result is not correct because after Jmat.Real.BIGGESTJSINT isPrime gives NaN

  if(Jmat.Real.isPrime(x)) return x;
  var d = x % 6;
  var e = 6 - d;
  var i = 0;
  var result = 0;
  for(;;) {
    if(Jmat.Real.isPrime(x - i - d + 1)) result = x - i - d + 1;
    else if(Jmat.Real.isPrime(x - i - d - 1)) result = x - i - d - 1;

    if((!result || (x - result) > (i + e - 1)) && Jmat.Real.isPrime(x + i + e - 1)) result = x + i + e - 1;
    else if((!result || (x - result) > (i + e + 1)) && Jmat.Real.isPrime(x + i + e + 1)) result = x + i + e + 1;

    if(result) return result;

    i += 6;
  }
};

Jmat.Real.nextPrime = function(value) {
  var x = Math.floor(value);
  if(x < 2) return 2;
  if(x < 3) return 3;
  if(x == Infinity || x != x) return NaN;
  if(x >= 9007199254740881) return NaN; //largest supported prime in floating point precision, after this will infinite loop because after Jmat.Real.BIGGESTJSINT isPrime gives NaN

  var m = x % 6;
  var step = 2;
  if(m == 0 || m == 5) {
    x += (m == 0 ? 1 : 2);
    step = 4;
  } else {
    x += (5 - m);
  }
  for(;;) {
    if(Jmat.Real.isPrime(x)) return x;
    x += step;
    step ^= 6; //swap step between 2 and 4
  }
};

Jmat.Real.previousPrime = function(value) {
  var x = Math.ceil(value);
  if(x <= 2) return NaN; // there is no lower prime
  if(x <= 3) return 2;
  if(x <= 5) return 3;
  if(x <= 7) return 5;
  if(x == Infinity || x != x) return NaN;
  if(x > Jmat.Real.BIGGESTJSINT) return NaN; //too large for the floating point's integer precision, result will not make sense

  var m = x % 6;
  var step = 2;
  if(m == 0 || m == 1) {
    x -= (m + 1);
    step = 4;
  } else {
    x -= (m - 1);
  }
  for(;;) {
    if(Jmat.Real.isPrime(x)) return x;
    x -= step;
    step ^= 6; //swap step between 2 and 4
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
    [1, 7, 21, 35, 35, 21, 7, 1]
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
  if(!Jmat.Real.isInt(x) || !Jmat.Real.isInt(y)) return NaN; //prevents infinite loop if both x and y are NaN. Also, reals are not supported here.
  if(Math.abs(x) > Jmat.Real.BIGGESTJSINT || Math.abs(y) > Jmat.Real.BIGGESTJSINT) return NaN; // does not work above JS integer precision
 //Euclid's algorithm
 for(;;) {
   if(y == 0) return Math.abs(x); //if x or y are negative, the result is still positive by the definition
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

Jmat.Real.near = function(x, y, epsilon) {
  // works also for infinities
  return x >= y - epsilon && x <= y + epsilon;
};

/*
Precision must be near 0 but slightly larger, e.g. 0.001 for 3 digits of precision, 1e-5 for 5 digits, ...
That many digits must match, starting from the first non-zero digit.
That means, if one value is zero and the other is not, no matter how close to zero the other is, this function will always return false.
It also always returns false if the signs differ.
Examples:
Jmat.Real.relnear(1.25e-300, 1.26e-300, 1e-2) --> true
Jmat.Real.relnear(1.25e-300, 1.26e-300, 1e-3) --> false
*/
Jmat.Real.relnear = function(x, y, precision) {
  if(x == y) return true;
  if(x == 0 || y == 0) return false; // case were both are 0 already handled with previous comparison
  if((x < 0) != (y < 0)) return false;
  x = Math.abs(x);
  y = Math.abs(y);
  var d = (x > y) ? (x / y) : (y / x);
  return d < 1 + precision;
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


//arbitrary log: log_y(x)
//warning: base y is second argument
Jmat.Real.logy = function(x, y) {
  return Math.log(x) / Math.log(y);
};

// Returns the number of leading zero bits in the 32-bit binary representation of x
// Gives floor of log2 of x by doing 31 - clz32(x)
// Gives num bits of x by doing 32 - clz32(x)
// Only guaranteed to work for numbers less than 32 bits
Jmat.Real.clz32 = Math['clz32'] || function(x) {
  var result = 0;
  while(x > 0) {
    x = Math.floor(x / 2);
    result++;
  }
  return 32 - result;
}

//NOTE: floating point version. For integer log2 use ilog2,
//because e.g. on 8 this gives 2.9999999999999996 (official Math.log2 too)
Jmat.Real.log2 = Math.log2 || function(x) {
  return Math.log(x) / Math.LN2;
};

// Integer log2: the floor of log2
Jmat.Real.ilog2 = function(x) {
  if(x <= 0) return NaN;
  if(x < 2147483648) return 31 - Jmat.Real.clz32(x);
  return Math.floor(Jmat.Real.log2(Math.floor(x) + 0.5));
};

Jmat.Real.getNumBits = function(x) {
  return Jmat.Real.ilog2(Math.abs(x)) + 1;
};

Jmat.Real.log10 = Math.log10 || function(x) {
  return Math.log(x) / Math.LN10;
};

Jmat.Real.root = function(x, y) {
  return Math.pow(x, 1 / y);
};

////////////////////////////////////////////////////////////////////////////////

// Faddeeva function: w(z) = exp(-z^2)*erfc(-iz). Also known as Faddeyeva or w(z) (not to be confused with LambertW)
// Returns complex result as an array [re, im], for complex input z = i*x + y. The result is real for pure imaginary input (arbitrary complex otherwise)
// Note that Jmat.Real does not depend on Jmat.Complex so does not use that datatype to represent the complex number.
// This complex function is used in Jmat.Real because some real erf related functions use imaginary numbers, and, this function gives very good accuracy for everything erf related.
// Based on paper "More Efficient Computation of the Complex Error Function, G. P. M. POPPE and C. M. J. WIJERS "
Jmat.Real.faddeeva = function(x, y) {
  var R = Jmat.Real;
  var invsqrtpi2   = 2 / R.SQRTPI;

  // complex exponentiation exp(x + yi)
  var cexp = function(x, y) {
    var e = Math.exp(x);
    return [e * Math.cos(y), e * Math.sin(y)];
  };

  // complex multiplication (a + bi) * (c + di)
  var cmul = function(a, b, c, d) {
    return [a * c - b * d, a * d + b * c];
  };

  // reciproke of complex number (x + yi)
  var cinv = function(x, y) {
    var d = x * x + y * y;
    return [x / d, -y / d];
  };

  // square of rho, used to determine which algorithm to use and what tweaking
  // parameters inside the algorithm. All magic numbers related to rho are as
  // suggested in the paper.
  var rho2 = (x / 6.3 * x / 6.3) + (y / 4.4 * y / 4.4);

  // Methods 2 and 3 below require positive imaginary part y, so transform
  // using the transformation w(-x) = 2 * exp(-x*x) - w(x).
  if(y < 0 && rho2 >= 0.292 * 0.292) {
    // For large negative pure imaginary values starting at -26.64, the code
    // starts returning NaN. Return Infinity instead (it is large positive overflow).
    if(x == 0 && y < -26.64) return [Infinity, 0]

    var e = cexp(y * y - x * x, -2 * x * y); // exp(-z*z)
    var f = R.faddeeva(-x, -y);
    return [2 * e[0] - f[0], 2 * e[1] - f[1]];
  }

  var result = [0, 0];

  if(rho2 < 0.292 * 0.292) {
    // Method 1: Power series
    // Based on sum 7.1.5 from Handbook of Mathematical Functions
    // erf(z) = 2/sqrt(pi) * SUM_n=0..oo (-1)^n * z^(2n+1) / (n! * (2n+1))
    // and then w(z) = e^(-z^2) * (1 - erf(iz))
    var s = (1 - 0.85 * y / 4.4) * Math.sqrt(rho2);
    var n = Math.ceil(6 + 72 * s); // ideal number of iterations
    var kk = 1;
    var zz = [y * y - x * x, -2 * x * y]; // -z*z
    var t = [y, -x]; // holds iz^(2k+1)
    for(var k = 0; k < n; k++) {
      if(k > 0) {
        kk *= -k; // (-1)^k * k!
        t = cmul(t[0], t[1], zz[0], zz[1]);
      }
      result[0] += t[0] / (kk * (2 * k + 1));
      result[1] += t[1] / (kk * (2 * k + 1));
    }
    var e = cexp(zz[0], zz[1]); // exp(-z*z)
    result = cmul(e[0], e[1], result[0], result[1]);
    result[0] = e[0] - result[0] * invsqrtpi2;
    result[1] = e[1] - result[1] * invsqrtpi2;
  } else if(rho2 < 1.0) {
    // Method 2: Taylor series
    // The continued fraction is used for derivatives of faddeeva function. More
    // info about the continued fraction is in Method 3.
    // h is the heuristically chosen point at which the taylor series is considered.
    // The derivative of w(z) is w'(z) = -2z*w(z)+2i/sqrt(pi)), but also
    // with w_n(z) defined as exp(-z*z) * i^n * erfc(-iz), the nth derivative
    // of w(z) is equal to (2i)^n * n! * w_n(z), and in the tayler expansion the
    // factorial gets canceled out: w(z) = SUM_0..oo (2h)^n * w_n(z + ih).
    var s = (1 - y / 4.4) * Math.sqrt(1 - rho2);
    var nu = Math.ceil(16 + 26 * s) + 1; // ideal number of iterations for continued fraction
    var n = Math.ceil(7  + 34 * s) + 1; // ideal number of iterations for taylor series
    var h = 1.88 * s;

    // The first iterations only warm up the w_n's with continued fraction
    var w = [0, 0]; // w_n's
    for (var k = nu; k > n; k--) {
      w = cinv(2 * (y + k * w[0] + h), 2 * (k * w[1] - x)); // 0.5/(h - i*z + k*w)
    }
    // The next iterations run the taylor series, while keeping updating the continued fraction
    var hh = Math.pow(h * 2, n - 1);
    for (var k = n; k > 0; k--) {
      w = cinv(2 * (y + k * w[0] + h), 2 * (k * w[1] - x)); // 0.5/(h - i*z + k*w)
      result = cmul(result[0] + hh, result[1], w[0], w[1]); // (result + hh) * w
      hh /= (h * 2);
    }
    result[0] *= invsqrtpi2;
    result[1] *= invsqrtpi2;
  } else {
    // Method 3: Continued fraction
    // The continued fraction is evaluated as r_nu = 0, r_(n-1) = 0.5 / (-iz + (n + 1)r_n), r_0 is the final approximate result
    var nu = Math.ceil(3 + (1442 / (26 * Math.sqrt(rho2) + 77))) + 1; // ideal number of iterations
    for (var k = nu; k > 0; k--) {
      result = cinv(2 * (y + k * result[0]), 2 * (k * result[1] - x)); //  0.5/(-i*z + k*result)
    }
    result[0] *= invsqrtpi2;
    result[1] *= invsqrtpi2;
  }

  if(x == 0) result[1] = 0; // for pure imaginary input, result is pure real. Fix potential numerical problems, and cases of "-0".
  if(y == 0) result[0] = Math.exp(-x * x); // for pure real input, the real part of the output is exactly exp(-x * x). Fix numerical imprecisions when near zero.
  return result;
};

// erfcx(x) = exp(x^2) * erfc(x): the scaled complementary error function
Jmat.Real.erfcx = function(x) {
  return Jmat.Real.faddeeva(0, x)[0]; //erfcx(x) = faddeeva(ix)
};

Jmat.Real.erf = function(x) {
  /*
  For verification, comparing with following 24 digit precision (our function gets about 14 digits correct):
  erf(0.00001) = 0.000011283791670578999349
  erf(0.1) = 0.112462916018284892203275
  erf(0.5) = 0.520499877813046537682746
  erf(1.5) = 0.966105146475310727066976
  erf(2.5) = 0.999593047982555041060435
  erf(3.5) = 0.999999256901627658587254
  erf(4.5) = 0.999999999803383955845711
  erf(5.5) = 0.999999999999992642152082
  erf(6.5) = 0.999999999999999999961578
  erf(3.5i) = 35282.287715171685310157997216i
  erf(3.5+3.5i) = 0.887129271239584272207414 + 0.015026380322129921373706i
  */
  var a = Math.exp(-x * x);
  if (x >= 0) return 1 - a * Jmat.Real.faddeeva(0, x)[0];
  else return a * Jmat.Real.faddeeva(0, -x)[0] - 1;
};

// erfc(x) = 1 - erf(x). This function gives numerically a better result if erf(x) is near 1.
Jmat.Real.erfc = function(x) {
  /*
  For verification, comparing with following 24 digit precision (our function gets about 14 digits correct):
  erfc(0.00001) = 0.999988716208329421000650
  erfc(0.1) = 0.887537083981715107796724
  erfc(0.5) = 0.479500122186953462317253
  erfc(1.5) = 0.033894853524689272933023
  erfc(2.5) = 0.000406952017444958939564
  erfc(3.5) = 7.430983723414127455236837 * 10^-7
  erfc(4.5) = 1.966160441542887476279160 * 10^-10
  erfc(5.5) = 7.357847917974398063068362 * 10^-15
  erfc(6.5) = 3.842148327120647469875804 * 10^-20
  erfc(3.5i) = 1 - 35282.287715171685310157997216i
  erfc(3.5+3.5i) = 0.112870728760415727792585 - 0.015026380322129921373706i
  */
  var a = Math.exp(-x * x);
  if (x >= 0) return a * Jmat.Real.faddeeva(0, x)[0];
  else return 2 - a * Jmat.Real.faddeeva(0, -x)[0];
};

//erfi(x) = -i erf(iz)
Jmat.Real.erfi = function(x) {
  var a = Math.exp(x * x);
  return a * Jmat.Real.faddeeva(x, 0)[1];
};

// D+(x) aka F(x)
Jmat.Real.dawson = function(x) {
  var a = Math.exp(-x * x);
  var w = Jmat.Real.faddeeva(x, 0)[1];
  return -(a - w) * (Jmat.Real.SQRTPI / 2);
};

// fast but inaccurate
Jmat.Real.erf_fast_ = function(x) {
  var neg = x < 0;
  if(neg) x = -x;

  if (x == 0) return 0;
  var t = 1 / (1 + 0.3275911 * x);
  var p = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
          t * (-1.453152027 + t * 1.061405429))));
  var result = 1.0 - p * Math.exp(-(x*x));

  if(neg) result = -result;
  return result;
};

// fast but inaccurate
Jmat.Real.erfc_fast_ = function(x) {
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

// fast but inaccurate
Jmat.Real.dawson_fast_ = function(x) {
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

// fast but inaccurate
Jmat.Real.erfi_fast_ = function(x) {
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

////////////////////////////////////////////////////////////////////////////////

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

Jmat.Real.degToRad = function(a) {
  return Math.PI * 2 * a / 360;
};

Jmat.Real.radToDeg = function(a) {
  return 360 * a / (Math.PI * 2);
};

// Like Math.round, but for fractional parts of 0.5, it is rounded to the nearest even value.
Jmat.Real.round = function(x) {
  // return Math.round(x);
  var l = Math.floor(x);
  var f = x - l;
  if(f == 0.5) return (l % 2 == 0) ? l : (l + 1)
  return (f < 0.5) ? l : (l + 1);
};

// Truncates towards zero
Jmat.Real.trunc = Math.trunc || function(x) {
  return (x < 0) ? Math.ceil(x) : Math.floor(x);
};

// Linear interpolation from a to b
Jmat.Real.lerp = function(a, b, x) {
  return (1 - x) * a + x * b;
};

// ECMAScript 5 doesn't have it
Jmat.Real.sinh = Math.sinh || function(x) {
  return (Math.exp(x) - Math.exp(-x)) / 2;
};

// ECMAScript 5 doesn't have it
Jmat.Real.cosh = Math.cosh || function(x) {
  return (Math.exp(x) + Math.exp(-x)) / 2;
};

// ECMAScript 5 doesn't have it
Jmat.Real.tanh = Math.tanh || function(x) {
  if(x > 354) return 1; // exp overflow
  return (Math.exp(2 * x) - 1) / (Math.exp(2 * x) + 1);
};

// ECMAScript 5 doesn't have it
Jmat.Real.asinh = Math.asinh || function(x) {
  if(x == -Infinity) {
    return x;
  } else {
    return Math.log(x + Math.sqrt(x * x + 1));
  }
};

// ECMAScript 5 doesn't have it
Jmat.Real.acosh = Math.acosh || function(x) {
  return Math.log(x + Math.sqrt(x * x - 1));
};

// ECMAScript 5 doesn't have it
Jmat.Real.atanh = Math.atanh || function(x) {
  return Math.log((1 + x) / (1 - x)) / 2;
};

// returns sqrt(x^2 + y^2), avoiding numerical underflow or overflow ; a companion to atan2
// Unlike Math.hypot from the JavaScript ES6 standard, this function does not support multiple arguments, only exactly two.
Jmat.Real.hypot = function(x, y) {
  x = Math.abs(x);
  y = Math.abs(y);
  var t = Math.min(x, y);
  x = Math.max(x, y);
  if(x == Infinity) return Infinity;
  t /= x;
  return x * Math.sqrt(1 + t * t);
};

//exp(x) - 1, with better precision for x around 0
Jmat.Real.expm1 = function(x) {
  if(Math.abs(x) < 1e-5) return x + x * x / 2 + x * x * x / 6;
  else return Math.exp(x) - 1;
};

////////////////////////////////////////////////////////////////////////////////

// Replicate the rest of JS Math library.

Jmat.Real.abs = Math.abs;
Jmat.Real.floor = Math.floor;
Jmat.Real.ceil = Math.ceil;
Jmat.Real.min = Math.min;
Jmat.Real.max = Math.max;
Jmat.Real.exp = Math.exp;
Jmat.Real.log = Math.log;
Jmat.Real.sqrt = Math.sqrt;
Jmat.Real.pow = Math.pow;
Jmat.Real.sin = Math.sin;
Jmat.Real.cos = Math.cos;
Jmat.Real.tan = Math.tan;
Jmat.Real.asin = Math.asin;
Jmat.Real.acos = Math.acos;
Jmat.Real.atan = Math.atan;
Jmat.Real.atan2 = Math.atan2;

////////////////////////////////////////////////////////////////////////////////

Jmat.Real.isLeapYear = function(y) {
  return (y % 400 == 0) || (y % 4 == 0 && y % 100 != 0);
};

Jmat.Real.montharray_ = [-1 /*there is no month 0*/, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; //february must be adjusted for leap year

Jmat.Real.monthLength = function(month, leap) {
  return (leap && month == 2) ? 29 : Jmat.Real.montharray_[month];
};


//number of days since the first day of year 0. 1 january of the year 0 is 0. 2 january is 1, etc...
//only for Gregorian calendar, does not take Julian calendar into account
Jmat.Real.numDaysSince0 = function(year, month, day) {
  var R = Jmat.Real;

  //number of leap years before this year (year 0 is considered leap)
  var numleap = year == 0 ? 0 : (R.idiv(year - 1, 4) - R.idiv(year - 1, 100) + R.idiv(year - 1, 400) + 1);
  var yeardays = year * 365 + numleap; //days of years before this year

  var feb = month > 2 ? (Jmat.Real.isLeapYear(year) ? 1 : 2) : 0; //days to subtract from formula below due to february
  var aug = (month > 8 && month % 2 == 1) ? 1 : 0; //correction because august shifts who has 31 days
  var monthdays = (month - 1) * 30 + R.idiv(month, 2) - feb + aug;

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

  // TODO: replace the for loop with shorter expressions
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

//determines day of week (0=sun, 1=mon, 2=tue, 3=wed, 4-thu, 5-fri, 6=sat), given year (e.g. 2014), month (1-12), day (1-31).
Jmat.Real.dayOfWeek = function(y, m, d) {
  var R = Jmat.Real;
  d += (m < 3) ? (y--) : (y - 2);
  return (R.idiv(23 * m, 9) + d + 4 + R.idiv(y, 4) - R.idiv(y, 100) + R.idiv(y, 400)) % 7;
};

////////////////////////////////////////////////////////////////////////////////

/*
Here are a few matrix algorithms in Jmat.Real. Much more algorithms are in
Jmat.Matrix. However, the ones here work only on real numbers, and do not use
any special class for matrices, just 2D arrays. The length of the array is the
matrix height, sub-arrays are rows. Sometimes a column vector is given as a 1D
array.
*/

Jmat.Real.matrix_add = function(a, b) {
  if(a.length != b.length || a[0].length != b[0].length) return undefined;
  var c = [];
  for(var y = 0; y < a.length; y++) {
    c[y] = [];
    for(var x = 0; x < a[y].length; x++) {
      c[y][x] = a[y][x] + b[y][x];
    }
  }
  return c;
};

Jmat.Real.matrix_sub = function(a, b) {
  if(a.length != b.length || a[0].length != b[0].length) return undefined;
  var c = [];
  for(var y = 0; y < a.length; y++) {
    c[y] = [];
    for(var x = 0; x < a[y].length; x++) {
      c[y][x] = a[y][x] - b[y][x];
    }
  }
  return c;
};

Jmat.Real.matrix_mul = function(a, b) {
  // TODO: add strassen algorithm
  var m = a.length;
  var n = a[0].length;
  var p = b[0].length;
  if(n != b.length) return undefined;
  var result = [];
  for (var y = 0; y < m; y++) result[y] = [];
  var temp = [];
  for (var x = 0; x < p; x++) {
    for (var z = 0; z < n; z++) temp[z] = b[z][x]; // copy for better caching (faster)
    for (var y = 0; y < m; y++) {
      var e = 0;
      for (var z = 0; z < n; z++) e += a[y][z] * temp[z];
      result[y][x] = e;
    }
  }
  return result;
};

Jmat.Real.matrix_mulr = function(a, v) {
  var r = [];
  for(var y = 0; y < a.length; y++) {
    r[y] = [];
    for(var x = 0; x < a[y].length; x++) {
      r[y][x] = a[y][x] * v;
    }
  }
  return r;
};

Jmat.Real.matrix_divr = function(a, v) {
  var r = [];
  for(var y = 0; y < a.length; y++) {
    r[y] = [];
    for(var x = 0; x < a[y].length; x++) {
      r[y][x] = a[y][x] / v;
    }
  }
  return r;
};

Jmat.Real.matrix_transpose = function(m) {
  var result = [];
  for(var y = 0; y < m[0].length; y++) {
    result[y] = [];
    for(var x = 0; x < m.length; x++) {
      result[y][x] = m[x][y];
    }
  }
  return result;
};

// Bring a to reduced row echelon form, in-place (modifies the input object)
Jmat.Real.matrix_rref = function(a) {
  var h = a.length;
  var w = a[0].length;

  var swaprow = function(matrix, y0, y1) {
    var temp = matrix[y0];
    matrix[y0] = matrix[y1];
    matrix[y1] = temp;
  };

  // subtracts f * y0 from y1 (modifying row y1), starting from x.
  var subrow = function(matrix, x, y0, y1, f) {
    var w = matrix[0].length;
    for (var i = x; i < w; i++) {
      matrix[y1][i] -= f * matrix[y0][i];
    }
  };

  // only starts at x rather than from the beginning
  var mulrow = function(matrix, x, y, v) {
    for (var i = x; i < w; i++) {
      matrix[y][i] = matrix[y][i] * v;
    }
  };

  var pivots = []; // x coordinate of pivot in each row (except the zero rows at the end, so may have smaller length than h)

  // gaussian elimination
  var k2 = 0; //next row, equal to k unless there are zero-rows
  for(var k = 0; k < w; k++) {
    var n = Jmat.Real.argmax(k2, h, function(i) { return Math.abs(a[i][k]); });
    if (a[n][k] == 0) continue; // singular, leave row as is
    if(k2 != n) swaprow(a, k2, n);
    mulrow(a, k, k2, 1 / a[k2][k]); // pivot is now 1
    for (var i = k2 + 1; i < h; i++) {
      if(a[i][k] != 0) subrow(a, k, k2, i, a[i][k]);
      a[i][k] = 0; // make extra-sure it's 0, avoid numerical imprecision
    }
    pivots.push(k);
    k2++;
    if(k2 >= h) break;
  }

  //now bring from row echolon form to reduced row echolon form
  for(var k = 0; k < pivots.length; k++) {
    var p = pivots[k];
    for(var y = k - 1; y >= 0; y--) {
      if(a[y][p] != 0) subrow(a, p, k, y, a[y][p]);
      a[y][p] = 0; // make extra-sure it's 0, avoid numerical imprecision
    }
  }

  return a;
};

// solves A*X = B. B and result are column vector given as 1D array.
Jmat.Real.matrix_solve = function(a, b) {
  var aug = [];
  for(var y = 0; y < a.length; y++) {
    aug[y] = [];
    for(var x = 0; x < a[y].length; x++) aug[y][x] = a[y][x];
    aug[y].push(b[y] || 0);
  }

  var r = Jmat.Real.matrix_rref(aug);

  // If we got a non-square a as input, our output size must be the width, not height, of a.
  var result = [];
  for(var i = 0; i < r[0].length - 1; i++) result[i] = r[i][r[i].length - 1];
  return result;
};

// This is a matrix algorithm, but is in Jmat.Real because it operates on real elements, and you can use the algorithm without Matrix class.
// Jacobi eigenvalue algorithm for real symmetric matrix
// a is n*n 2D array with input and output matrix (real symmetric), contains eigenvalues on diagonal after the algorithm (sorted)
// v is n*n 2D output array (may be initialized as []), matrix which will contain eigenvectors as rows after the algorithm (normalized)
// n is matrix size
// opt_epsilon is precision for when to stop the iterations (default 1e-15)
Jmat.Real.matrix_jacobi = function(a, v, n, opt_epsilon) {
  var epsilon = opt_epsilon == undefined ? 1e-15 : opt_epsilon;

  // Make identity
  for(var y = 0; y < n; y++) {
    v[y] = [];
    for(var x = 0; x < n; x++) {
      v[y][x] = (x == y) ? 1 : 0;
    }
  }

  // Sum of squares of all off-diagonal elements
  var off2 = 0;
  for(var y = 0; y < n; y++) {
    for(var x = y + 1; x < n; x++) {
      if(x != y) {
        off2 += 2 * a[y][x] * a[y][x];
      }
    }
  }

  while(off2 > epsilon) {
    for(var y = 0; y < n; y++) {
      for(var x = y + 1; x < n; x++) {
        if(a[y][x] * a[y][x] <= off2 / (2 * n * n)) continue; // Too small
        off2 -= 2 * a[y][x] * a[y][x];
        // Jacobi rotation coefficients
        var beta = (a[x][x] - a[y][y]) / (2 * a[y][x]);
        var t = Math.sign(beta) / (Math.abs(beta) + Math.sqrt(beta * beta + 1));
        var s = 1 / (Math.sqrt(t * t + 1));
        var c = s * t;
        // Rotate rows of A
        for(var k = 0; k < n; k++) {
          var tmp = a[k][y];
          a[k][y] = s * a[k][x] + c * tmp;
          a[k][x] = c * a[k][x] - s * tmp;
        }
        // Rotate columns of A and V
        for(var k = 0; k < n; k++) {
          var tmp = a[y][k];
          a[y][k] = s * a[x][k] + c * tmp;
          a[x][k] = c * a[x][k] - s * tmp ;
          tmp = v[y][k];
          v[y][k] = s * v[x][k] + c * tmp;
          v[x][k] = c * v[x][k] - s * tmp;
        }
      }
    }
  }

  // Sort eigenvalues if needed
  for(var k = 0; k < n; k++) {
    var m = k;
    for(var l = k + 1; l < n; l++) {
      if(a[l][l] > a[m][m]) m = l;
    }
    if(k != m) {
      var tmp = a[m][m];
      a[m][m] = a[k][k];
      a[k][k] = tmp;
      for(var l = 0; l < n; l++) {
        var tmp = v[m][l];
        v[m][l] = v[k][l];
        v[k][l] = tmp;
      }
    }
  }
};

////////////////////////////////////////////////////////////////////////////////

//f is function taking integer index as parameter and returning a real
//returns index belonging to max return value of f in index range [s, e)
Jmat.Real.argmax = function(s, e, f) {
  var m = f(s);
  var b = s;
  for(var i = s + 1; i < e; i++) {
    var mi = f(i);
    if(mi > m) {
      m = mi;
      b = i;
    }
  }
  return b;
};
