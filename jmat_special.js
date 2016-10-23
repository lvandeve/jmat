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
jmat_special.js extends Jmat.Real and Jmat.Complex with more functions. Has three categories of functions:
1. special mathematical functions such as hypergeometric and polylog (very special ones, more basic "special" mathematical functions are in jmat_real.js and jmat_complex.js)
2. numerical algorithms for complex functions such as root finding, integration, differentiation
3. statistical distributions (pdf/pmf, cdf and qf of each)

IMPORTANT NOTE: Admittedly, this file contains the least numerically stable functions of Jmat.js.
                Some only work well for a certain input range. So verify for your usecase before trusting these.

Some more basic special functions are in jmat_complex.js and jmat_real.js, that includes erf, gamma and lambertw.

Overview of some functionality, this lists function names in "Jmat.Complex.":
-Tetration: tetration, tetrational
-Minkowski question mark: minkowski
-gamma function related: loggamma, gamma_inv, digamma, trigamma, polygamma
-incomplete gamma: incgamma_lower, incgamma_upper, gamma_p, gamma_q
-means: agm, ghm
-bessel: besselj, bessely, besseli, besselk, hankelh1, hankelh2, struveh, struvek, struvel, struvem
-airy: airy, bairy, airy_deriv, airy_deriv
-zeta: zeta, eta, lambda, hurwitzzeta
-beta: beta, incbeta, betai
-hypergeometric: hypergeometric0F1, hypergeometric1F1, hypergeometric2F1, hypergeometric, regularized_hypergeometric
-polylog: dilog, trilog, polylog
-theta function: theta1, theta2, theta3, theta4
-error function related: erf_inv, erfc_inv
-algorithms: rootfind, integrate, differentiate, doSummation, doProduct, powerSeries
-statistical distributions: {pdf,pmf,cdf,qf}_{uniform,standardnormal,normal,lognormal,cauchy,studentt,chi_square,logistic,gamma,beta,fisher,weibull,exponential,laplace,bernoulli,binomial,poisson}
*/

////////////////////////////////////////////////////////////////////////////////

// Confluent hypergeometric limit function 0F1(; a; z) --> no "upper" term
Jmat.Complex.hypergeometric0F1 = function(a, z) {
  var C = Jmat.Complex;

  // It's infinity for negative integer a. However, the loop fails to encounter that term if rather large negative integer due to not enough steps.
  // Even though it seems as if it should be a finite value because all non integers around it converge to the same apparently...
  if(C.isNegativeInt(a)) return C(Infinity, Infinity); //undirected infinity

  // The series does not converge well for large negative z with smallish a, because then terms with opposite signs are supposed to cancel out.
  // Use this expansion into two 2F0 calls instead.
  // Some difficult cases that this handles: hypergeometric0F1(-23.5, -126.5625) = 692.923, hypergeometric0F1(6, -400) = 4.59650497414 * 10^-6
  // The values -80 and z.abssq() / 3 have been chosen to make it work for BesselJ in the zone where it uses 0F1. This may be need tweaking for other applications.
  if(z.re < -80 && a.abssq() < z.abssq() / 3) {
    var g = C.gamma(a).divr(2 * C.SQRTPI.re);
    var s = C.sqrt(z.neg());
    var c0 = C.exp(s.muli(-2)).mul(s.muli(-1).pow(a.rsub(0.5)));
    var h0 = C.hypergeometric([a.subr(0.5), a.rsub(3/2)], [], s.muli(4).inv().neg());
    var c1 = C.exp(s.muli(2)).mul(s.muli(1).pow(a.rsub(0.5)));
    var h1 = C.hypergeometric([a.subr(0.5), a.rsub(3/2)], [], s.muli(4).inv());
    return g.mul(c0.mul(h0).add(c1.mul(h1)));
  }

  // The summation is SUM_(n=0..inf) z^n / ((a)_n * n!)
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  // This loop converges very nicely and result is quite accurate for all complex input values (unlike for 2F1)
  var r = C.ONE;
  var result = C.ZERO;
  for(var n = 0; n < 30; n++) {
    if (n > 0) {
      if(!r.eqr(0)) r = r.div(a.addr(n - 1)); // the Pochhammer. The if avoids 0/0
      r = r.mul(z).divr(n); // the z^n / n!
    }
    if(r.eqr(0)) break;
    result = result.add(r);
  }

  return result;
};

// The basic series for 1F1
Jmat.Complex.hypergeometric1F1_series_ = function(a, b, z) {
  var C = Jmat.Complex;

  // The summation is SUM_(n=0..inf) ((a)_n / (b)_n) * z^n / n!
  // See for loop in hypergeometric function for explanation on above notation and how implementation works
  // This loop converges very nicely (unlike for 2F1)
  var r = a.div(b).mul(z).divr(1);
  var result = C(1);
  // It should break out much sooner than that many loops, unless a is very big (1000) or b is tiny (0.1)
  for(var n = 1; n < 200; n++) {
    if (n > 1) {
      r = r.mul(a.addr(n - 1));
      if(!r.eqr(0)) r = r.div(b.addr(n - 1));
      r = r.mul(z);
      r = r.divr(n);
    }
    if(r.eqr(0) || r.abssq() / result.abssq() < 1e-28) break;
    result = result.add(r);
  }

  return result;
};

// The asymptotic expansion for 1F1
// Does NOT work for negative integer b. Even the gammaDiv_ cannot fix it if a and b are both negative integer.
Jmat.Complex.hypergeometric1F1_asymp_ = function(a, b, z) {
  var C = Jmat.Complex;

  if(C.isNegativeInt(b)) return C(NaN); // not supported, even if a is compensating negative int

  // asymptotic expansion, abramowitz&stegun 13.5.1
  var sum0 = C.ONE;
  var sum1 = C.ONE;
  var s0 = C.ONE;
  var s1 = C.ONE;
  // The loop should end sooner than that many loops except for high values, in which case this diverged and gives bad result
  for(var i = 1; i < 30; i++) {
    s0 = s0.mul(a.addr(i - 1)).mul(a.sub(b).addr(i)).divr(i).div(z.neg());
    s1 = s1.mul(b.sub(a).addr(i - 1)).mul(a.rsub(i)).divr(i).div(z);
    sum0 = sum0.add(s0);
    sum1 = sum1.add(s1);
    if(s0.abssq() < 1e-28 && s1.abssq() < 1e-28) break;
  }

  if(s0.abssq() > 1e-10 || s1.abssq() > 1e-10) return C(NaN); // no convergence

  //var sign = (z.arg() <= -Math.PI / 2) ? -1 : 1; // This is what it should be according to abramowitz&stegun ...
  var sign = (z.arg() <= 0) ? -1 : 1; // ... but this is what gives the correct effect here instead for some reason

  if(z.abs() < 100) {
    var g0 = C.gammaDiv_(b, b.sub(a)).mul(C.exp(a.muli(sign * Math.PI))).div(z.pow(a));
    var g1 = C.gammaDiv_(b, a).mul(C.exp(z)).mul(z.pow(a.sub(b)));
    return g0.mul(sum0).add(g1.mul(sum1));
  } else {
    // Written in logarithmic form so that it works for very large z too
    var g0 = C.loggammaDiv_(b, b.sub(a)).add(a.muli(sign * Math.PI)).sub(C.log(z).mul(a));
    var g1 = C.loggammaDiv_(b, a).add(z).add(C.log(z).mul(a.sub(b)));
    return C.exp(g0.add(C.log(sum0))).add(C.exp(g1.add(C.log(sum1))));
  }
};

// Rational approximation of 1F1 (Luke Y.L. 1977 chapter XV)
// For extreme values, requires too much iterations (over 300), but slightly better in some places
Jmat.Complex.hypergeometric1F1_rational_ = function(a, b, z) {
  var C = Jmat.Complex;

  var e = function(a, c, z) {
    //var z2 = C(2).pow(z.mul(z)); // 2^(z^2)
    //var z2 = C(2).pow(z).powr(2); // (2^z)^2
    var z2 = z.mul(z);
    var z3 = z2.mul(z);
    var c2 = c.mul(c.inc()); // (c)_2, rising pochhammer of c

    var b0 = C(1);
    var b1 = C(1).add(a.inc().mul(z).div(c).divr(2));
    var b2 = C(1).add(a.addr(2).mul(z).div(c.inc()).divr(2)).
                  add(a.inc().mul(a.addr(2)).mul(z2).div(c2).divr(12));

    var azc = a.mul(z).div(c);
    var a0 = C(1);
    var a1 = b1.sub(azc);
    var a2_a = a.addr(2).mul(z).div(c.inc()).divr(2).inc();
    var a2_b = a.mul(a.inc()).mul(z2).div(c2).divr(2);
    var a2 = b2.sub(azc.mul(a2_a)).add(a2_b);

    //var b3 = C(1).add(a.addr(3).mul(z).div(c.addr(2)).divr(2))
    //             .add(a.addr(2).mul(a.addr(3)).mul(z2).div(c.addr(1)).div(c.addr(2)).divr(10))
    //             .add(a.addr(1).mul(a.addr(2)).mul(a.addr(3)).mul(z3).div(c).div(c.addr(1)).div(c.addr(2)).divr(120));

    var an, bn;

    for(var n = 3; n < 50; n++) {
      var n2 = n * 2;
      var f1 = a.rsub(n-2).div(
               c.addr(n-1).mulr(2 * (n2-3)));
      var f2 = a.addr(n).mul(a.addr(n-1)).div(
               c.addr(n-1).mul(c.addr(n-2)).mulr(4 * (n2-3) * (n2-1)));
      var f3 = a.addr(n-2).mul(a.addr(n-1)).mul(a.rsub(n-2)).div(
               c.addr(n-3).mul(c.addr(n-2)).mul(c.addr(n-1)).mulr(8*(n2-3)*(n2-3)*(n2-5))).neg();
      var f4 = a.addr(n-1).mul(c.rsub(n-1)).div(
               c.addr(n-2).mul(c.addr(n-1)).mulr(2 * (n2-3))).neg();

      var e1 = f1.mul(z).inc();
      var e2 = f2.mul(z).add(f4).mul(z);
      var e3 = f3.mul(z3);
      an = a2.mul(e1).add(a1.mul(e2)).add(a0.mul(e3));
      bn = b2.mul(e1).add(b1.mul(e2)).add(b0.mul(e3));

      a0 = a1;
      a1 = a2;
      a2 = an;
      b0 = b1;
      b1 = b2;
      b2 = bn;
    }

    var test = a1.div(b1);
    var result = an.div(bn);

    if(test.sub(result).abs() > 1e-10) return C(NaN); // no good convergence

    return result;
  };

  return e(a, b, z.neg());
};


// Confluent hypergeometric function of the first kind 1F1(a; b; z), a.k.a. Kummer's function M (approximation)
// Try similar online: http://keisan.casio.com/exec/system/1349143651
Jmat.Complex.hypergeometric1F1 = function(a, b, z) {
  var C = Jmat.Complex;

  // Kummer's transformation to make |a| a bit smaller if possible
  if(b.sub(a).abssq() < a.abssq()) return C.hypergeometric1F1(b.sub(a), b, z.neg()).mul(C.exp(z));

  // Special values
  if(a.eq(b)) return C.exp(z); // if a and b are equal, they cancel out and this becomes 0F0, which is the same as exp(z)
  if(C.isNegativeInt(b) && !(C.isNegativeInt(a) && a.re >= b.re)) return C(Infinity);
  if(a.eqr(0)) return C(1);
  if(z.eqr(0)) return C(1);
  if(z.eqr(Infinity)) return C(Infinity, Infinity); // not for neg or complex infinite z

  var r1 = z.mul(a).div(b).abs(); // |z*a/b| dictates how difficult the calculation will be. The lower the easier. If low, the regular series works.

  if(r1 < 8 /*|| C.isNegativeInt(b)*/) {
    return C.hypergeometric1F1_series_(a, b, z);
  }

  var result = C.hypergeometric1F1_asymp_(a, b, z); // returns NaN if it didn't converge
  if(C.isNaN(result)) result = C.hypergeometric1F1_rational_(a, b, z);

  return result;
};

// Hypergeometric series 2F1(a, b; c; z), aka Gauss hypergeometric function or "The hypergeometric". Approximation.
// Try similar online: http://keisan.casio.com/exec/system/1349143084, or wolframalpha, HyperGeometric2F1[1.5, -0.5, 2.5, 0.25+i]
// I debugged this function to find all the areas in the complex or real input plane where it goes wrong, by 2D or complexdomain plotting this, and beta_i which uses this
// If c is a non-positive integer, the function is undefined and usually returns NaN except for some cases where a is integer.
// NOTE: For some inputs this function is very imprecise. For example for |z| near 1 due to slow convergence, and,
//   in case of combinations of a,b,c giving integers causing singularities that should cancel out where our
//   non-arbitrary precision makes it hard to calculate precise limits. Many functions defined in terms of 2F1 happen to
//   hit exactly those cases.
// TODO: make more precise, e.g. hypergeometric2F1(1.1, -0.9, 2.1, 3) and hypergeometric2f1(6,6.5,7,0.4444444444444444)
Jmat.Complex.hypergeometric2F1 = function(a, b, c, z) {
  var C = Jmat.Complex;

  if(z.abs() > 1.0001 /*not 1 to avoid infinite loop if abs 1/z is also > 1 due to numeric problems, e.g. for 0.9726962457337884+i0.23208191126279865*/) {
    // The series converges only for |z| < 1. But there are some linear transformations
    // that can convert it to a different z. There are conditions though, and some are
    // more complex than others (e.g. requiring gamma functions).
    // TODO: with only those two supported transformations below, there is probably a lot wrong,
    //       for example, the points at exp(pi*i/3) and exp(-pi*i/3) are known to be not covered by this
    //       To fix, do as in the paper "Fast computation of the Gauss hypergeometric function with all its parameters complex with application to the Poschl-Teller-Ginocchio potential wave functions"

    // Linear transformations to bring |z| in value < 1
    var z2 = z.div(z.dec());
    if(z2.abs() < 0.75) { // the if can in theory do "< 1", but then a value like z=0.3+4i makes z2 very close to 1, and it has bad numeric results for such values
      // z / (z - 1)
      return C.ONE.sub(z).pow(a.neg()).mul(C.hypergeometric2F1(a, c.sub(b), c, z2));
    } else {
      // 1 / z
      var t = C.isNegativeIntOrZero;
      // This if probably covers more cases than needed, but there are a lot of cases. Even gammaDiv22_ doesn't help out. E.g. if a,b,c are all 0.5, the limit is different than gammaDiv22_.
      if(t(a) || t(b) || t(b) || C.isInt(a.sub(b)) || t(c.sub(a)) || t(c.sub(b))) {
        // Twiddle to take the limit to avoid singularities (they cancel it out, but, finding exact formula here is too hard)
        // TODO: improve this since this is very imprecise
        a = a.addr(0.1e-7);
        b = b.subr(0.2e-7);
        c = c.addr(0.3e-7);
      }
      var zi = z.inv();
      var za = z.neg().pow(a.neg());
      var zb = z.neg().pow(b.neg());
      var ga = C.gammaDiv22_(c, b.sub(a), b, c.sub(a));
      var gb = C.gammaDiv22_(c, a.sub(b), a, c.sub(b));
      var fa = C.hypergeometric2F1(a, C.ONE.sub(c).add(a), C.ONE.sub(b).add(a), zi);
      var fb = C.hypergeometric2F1(b, C.ONE.sub(c).add(b), C.ONE.sub(a).add(b), zi);
      var va = ga.mul(za).mul(fa);
      var vb = gb.mul(zb).mul(fb);
      return va.add(vb);
    }
  }

  var z2 = z.div(z.dec());
  if(z2.abs() < z.abs()) {
    // Same z / (z - 1) transform as above. Reason for doing this: the summation below converges faster for smaller absolute values of z.
    // Without this, for e.g. z = -0.75 and c < -3, it converges only after hundreds of steps.
    // TODO: in fact it's almost always possible to make |z| < 0.5 with the linear transformations. Use them better.
    return C.ONE.sub(z).pow(a.neg()).mul(C.hypergeometric2F1(a, c.sub(b), c, z2));
  }

  // TODO: avoid needing more than 30 iterations by always getting |z| < 0.5 with more clever use of transformations
  var num_iterations = 30;
  if(z.abs() > 0.5) num_iterations = 60;
  if(z.abs() > 0.75) num_iterations = 100;

  // The summation definition of the series, converges because |z| < 1

  // The summation is SUM_(n=0..inf) ((a)_n * (b)_n / (c)_n) * z^n / n!
  // Where (q)_n is the rising Pochhammer symbol x(x+1)*...*(x+n-1)
  // The variable r below gets updated every step with next factor of Pochhammer symbol, power and factorial values as the loop goes.
  var r = a.mul(b).div(c).mul(z);
  var result = C(1);
  for(var n = 1; n < num_iterations; n++) {
    if (n > 1) {
      r = r.mul(a.addr(n - 1));
      r = r.mul(b.addr(n - 1));
      // The if(!r.eqr(0)) is there so that if a==c or b==c and c is neg integer, we don't get 0/0 here but just 0
      if(!r.eqr(0)) r = r.div(c.addr(n - 1));
      r = r.mul(z);
      r = r.divr(n);
    }
    if(r.eqr(0)) break; // If a or b are negative integer, the loop will terminate early
    if(C.near(r, C.ZERO, 1e-15)) break; // In fact why not even break out if it's already converged very near zero
    result = result.add(r);
  }
  return result;

  // Alternative: integration formula (doesn't work well)
  /*var g = C.gammaDiv12_(c, b, c.sub(b));
  var i = C.integrate(C(0), C(1), function(t) {
    var n = t.pow(b.dec()).mul(C.ONE.sub(t).pow(c.sub(b).dec()));
    var d = C.ONE.sub(z.mul(t)).pow(a);
    return n.mul(d);
  });
  return g.mul(i);*/
};

// The core loop for generalized hypergeometric.
Jmat.Complex.generalized_hypergeometric_ = function(a, b, z) {
  // TODO: equal elements in a and b cancel out. Sort the arrays and check that. If the arrays are equal, it's the same as 0F0
  // TODO: if a contains more zeros than b, the result is always 1
  // TODO: return NaN or so ("unsupported") if diverging. Note that for a.length > b.length, even if in theory it only converges for 0, it returns practical values for tiny values like 1e-4, so give no NaN there.

  var r = z;
  for(var i = 0; i < a.length; i++) r = r.mul(a[i]);
  for(var i = 0; i < b.length; i++) r = r.div(b[i]);
  var result = Jmat.Complex(1);
  for(var n = 1; n < 40; n++) {
    if (n > 1) {
      for(var i = 0; i < a.length; i++) r = r.mul(a[i].addr(n - 1));
      for(var i = 0; i < b.length; i++) r = (r.eqr(0) ? r : r.div(b[i].addr(n - 1)));
      r = r.mul(z);
      r = r.divr(n);
    }
    if(r.eqr(0)) break;
    result = result.add(r);
  }

  return result;
};

// Ordinary generalized hypergeometric series
// a and b are arrays of complex values. E.g. Jmat.Complex.hypergeometric2F1([a, b, c], [d, e], z) would mean: hypergeometric3F2(a, b, c; d, e; z)
// Unless specific implementation exists, this uses basic series of the definition of the hypergeometric
// So if the series diverges (a.length > b.length), it will only be correct for very small |z| (near zero), if a.length == b.length, it will only converge for |z| < 1
// NOTE: This is the generalized hypergeometric! For the namesake 2F1, use hypergeometric2F1 instead.
Jmat.Complex.hypergeometric = function(a, b, z) {
  if(z.eqr(0)) return Jmat.Complex(1);
  if(a.length == 0 && b.length == 0) return Jmat.Complex.exp(z);
  if(a.length == 1 && b.length == 0) return Jmat.Complex.pow(z.rsub(1), a[0].neg());
  if(a.length == 0 && b.length == 1) return Jmat.Complex.hypergeometric0F1(b[0], z);
  if(a.length == 1 && b.length == 1) return Jmat.Complex.hypergeometric1F1(a[0], b[0], z);
  if(a.length == 2 && b.length == 1) return Jmat.Complex.hypergeometric2F1(a[0], a[1], b[0], z);

  return Jmat.Complex.generalized_hypergeometric_(a, b, z);
};

// regularized hypergeometric is the hypergeometric divided through product of gammas of b's
Jmat.Complex.regularized_hypergeometric = function(a, b, z) {
  var C = Jmat.Complex;
  var g = C(1);
  var twiddle = false; // twiddle negative integers a bit. TODO: handle this here and in other twiddle parts of the code more precise
  for (var i = 0; i < b.length; i++) {
    if(!twiddle && C.isInt(b[i]) && b[i].re <= 0) {
      var a2 = [];
      var b2 = [];
      for(var j = 0; j < a.length; j++) a2[j] = C.isInt(a[j]) ? a[j].addr(1e-7) : a[j];
      for(var j = 0; j < b.length; j++) b2[j] = C.isInt(b[j]) ? b[j].addr(1e-7) : b[j];
      a = a2;
      b = b2;
      twiddle = true;
    }
    g = g.mul(C.gamma(b[i]));
  }
  var h = C.hypergeometric(a, b, z);
  return h.div(g);
};


////////////////////////////////////////////////////////////////////////////////

// Associated legendre function P^mu_nu(z). Set mu to zero for Legendre function. Have nu integer for Legendre polynomial.
// opt_type is 1, 2 or 3, and gives a slightly different formula, see comment in function (1 and 2 are identical).
// NOTE: precision of this function is not high, could be less than 4 digits in some cases. TODO: improve
// NOTE: see inside the gaussian integration code for a recursive implementation of legendre polynomials
Jmat.Complex.legendrep = function(nu, mu, z, opt_type) {
  var C = Jmat.Complex;
  if(mu.eqr(0) && C.isInt(nu) && nu.re >= 0) {
    if(nu.re == 0) return C(1);
    if(nu.re == 1) return z;
    if(nu.re == 2) return z.mul(z).mulr(3).subr(1).divr(2);
  }
  // Type 1: same as type 2, but only defined in unit circle (just returns same as type 2 here)
  // Type 2: (1+z)^(mu/2) / (1-z)^(mu/2) * 2F1reg(-nu, nu+1, 1-mu, (1-z)/2)
  // Type 3: (1+z)^(mu/2) / (z-1)^(mu/2) * 2F1reg(-nu, nu+1, 1-mu, (1-z)/2)
  opt_type = opt_type || 2;
  if(opt_type < 1 || opt_type > 3) return undefined;
  var a = z.addr(1).pow(mu.divr(2));
  var b = opt_type == 3 ? z.subr(1).pow(mu.divr(2)) : z.rsub(1).pow(mu.divr(2));
  var h = C.regularized_hypergeometric([nu.neg(), nu.addr(1)], [mu.rsub(1)], z.rsub(1).divr(2));
  return a.div(b).mul(h);
};

// Associated legendre function Q^mu_nu(z). Set mu to zero for Legendre function. Have nu integer for Legendre polynomial.
// opt_type is 1, 2 or 3, and gives a slightly different formula (1 and 2 are identical).
// NOTE: precision of this function is not high, could be less than 4 digits in some cases. TODO: improve
Jmat.Complex.legendreq = function(nu, mu, z, opt_type) {
  var C = Jmat.Complex;

  opt_type = opt_type || 2;
  if(opt_type < 1 || opt_type > 3) return undefined;

  // returns exp(z*pi*i), adjusted for numerical precision problems (if z is a multiple of 0.5)
  var exppii = function(z) {
    if(z.im == 0 && Jmat.Real.isInt(z.re * 2)) {
      var n = Jmat.Real.mod(Math.round(z.re * 2), 4);
      if(n == 0) return C(1, 0);
      if(n == 1) return C(0, 1);
      if(n == 2) return C(-1, 0);
      if(n == 3) return C(0, -1);
    }
    return C.exp(z.muli(Math.PI));
  };

  if(opt_type == 3) {
    var s = nu.add(mu).addr(1);
    var a = C(Jmat.Real.SQRTPI).mul(C.gamma(s));
    var b = C(2).pow(nu.addr(1));
    var c = exppii(mu);
    var d = z.pow(s);
    var e = z.addr(1).pow(mu.divr(2)).mul(z.subr(1).pow(mu.divr(2)));
    var f = C.regularized_hypergeometric([s.divr(2), s.addr(1).divr(2)], [nu.addr(1.5)], z.mul(z).inv());
    return a.div(b).mul(c).div(d).mul(e).mul(f);
  } else {
    // In terms of legendrep. More expensive as two 2F1 calls, and, problems with mu equal to integer...
    if(C.isInt(mu)) mu = mu.addr(1e-8);
    var a = C(Math.PI).div(C.sin(mu.mulr(Math.PI)).mulr(2));
    var b = C.cos(mu.mulr(Math.PI));
    var p1 = C.legendrep(nu, mu, z, opt_type);
    var g = C.gammaDiv_(mu.add(nu).addr(1), C(1).sub(mu).add(nu));
    var p2 = C.legendrep(nu, mu.neg(), z, opt_type);
    return a.mul(b.mul(p1).sub(g.mul(p2)));
  }

  // The definition from Abramowitz and Stegun, and Wikipedia. Does not match implementations of various software, so not used.
  /*var s = nu.add(mu).addr(1);
  var a = C(Jmat.Real.SQRTPI).mul(C.gamma(s));
  var b = C(2).pow(nu.addr(1));
  var c = exppii(mu); // this exp is left out in some formulas. For integer mu, this is always -1 or +1, even that is sometimes left out.
  var d = z.pow(s);
  var e = z.mul(z).subr(1).pow(mu.divr(2));
  var f = C.regularized_hypergeometric([s.divr(2), s.addr(1).divr(2)], [nu.addr(1.5)], z.mul(z).inv());
  return a.div(b).mul(c).div(d).mul(e).mul(f);*/
};


////////////////////////////////////////////////////////////////////////////////

//using Stirling series
//natural logarithm of the gamma function, and more specific branch of the log
//to use this for log2 of the factorial (common in computing theory), use e.g.: function log2fac(n) { return Jmat.loggamma(n + 1).re / Math.log(2); }
// TODO: check if this small series is precise enough everywhere
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
  if(Jmat.Complex.isPositive(value) && value.re > 0.85) { //doesn't work for negative values, nor values smaller than somewhere around 0.85
    // Approximation for positive real x.
    var x = value.re;
    //c = sqrt(2 * pi) / e - gamma(k), where k = the positive zero of the digamma function (1.4616321449683623412626...)
    var c = 0.03653381448490041660;
    var lx = Math.log((x + c) / Math.sqrt(2 * Math.PI));
    return Jmat.Complex(lx / Jmat.Real.lambertw(lx / Math.E) + 0.5);
  }

  // TODO: this has problems for |z| > 20. Use more stable root finding
  // TODO: since the current digamma implementation is nothing more than a derivative approximation, I might as well use finvert_newton_noderiv. Try using digamma anyway if it's ever made more precise.
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
  // Has some noise here and there, but works better than series representation for large values
  if(z.re > 0) {
    var result = z.pow(s).div(s).mul(Jmat.Complex.hypergeometric1F1(s, s.inc(), z.neg()));
    if (!Jmat.Complex.isNaN(result)) return result;
  }


  // METHOD B: series expansion - has some problems with division through zero
  // sum_k=0.oo ((-1)^k / k!) * (z^(s+k) / (s + k))
  var result = Jmat.Complex(0);
  var kk = Jmat.Complex(1);
  var zz = z.pow(s);
  var sign = 1;
  var numit = z.abs() > 5 ? 50 : 30;
  for(var k = 0; k < numit; k++) {
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

Jmat.Complex.gamma_p_cache_ = []; // cache used because a is often constant between gamma_p calls

// regularized lower incomplete gamma function (lower gamma_regularized, regularized_gamma)
// P(a, x) = gamma(a, z) / GAMMA(a)
// Note: the derivative of this function is: (e^(-z) * z^(a-1))/GAMMA(a)
Jmat.Complex.gamma_p = function(a, z) {
  if(Jmat.Complex.isNegativeIntOrZero(a)) return Jmat.Complex(1);
  var g = Jmat.Complex.calcCache_(a, Jmat.Complex.gamma, Jmat.Complex.gamma_p_cache_); // gamma(a)
  return Jmat.Complex.incgamma_lower(a, z).div(g);
};

// regularized upper incomplete gamma function (upper gamma_regularized, regularized_gamma)
// Q(a, x) = 1 - P(a, x) = GAMMA(a, z) / GAMMA(a)
// Note: the derivative of this function is: -(e^(-z) * z^(a-1))/GAMMA(a)
Jmat.Complex.gamma_q = function(a, z) {
  if(Jmat.Complex.isNegativeIntOrZero(a)) return Jmat.Complex(0);
  return Jmat.Complex.ONE.sub(Jmat.Complex.gamma_p(a, z));
};

// One possible approximation series for inverse gamma P. Valid for real values with p in range 0-1 and positive a. Possibly more.
Jmat.Complex.gamma_p_inv_series_1_ = function(a, p) {
  var a1 = a.inc();
  var a1a = a1.mul(a1);
  var a2 = a.addr(2);
  var a2a = a2.mul(a2);
  var aa = a.mul(a);
  var aaa = aa.mul(a);
  var aaaa = aaa.mul(a);
  var a3 = a.addr(3);
  var a4 = a.addr(4);

  var c = [0, 1];
  c[2] = a.inc().inv();
  c[3] = a.mulr(3).addr(5).div(a1a.mul(a2).mulr(2));
  c[4] = aa.mulr(8).add(a.mulr(33)).addr(31).div(a1a.mul(a1).mul(a2).mul(a3).mulr(3));
  c[5] = aaaa.mulr(125).add(aaa.mulr(1179)).add(aa.mulr(3971)).add(a.mulr(5661)).addr(2888).div(a1a.mul(a1a).mul(a2a).mul(a3).mul(a4).mulr(24));

  var r = p.mul(Jmat.Complex.gamma(a1)).pow(a.inv());

  return Jmat.Complex.powerSeries(c, c.length, Jmat.Complex.ZERO, r);
};

// Inverse regularized gamma P
// Finds z for p == gamma_p(a, z)
// This is an *approximation* of inverse of incomplete regularized gamma (P) - it works for real values with p in range 0-1 and positive a. Good enough for qf of chi square distribution
Jmat.Complex.gamma_p_inv = function(a, p) {
  // Return NaN for unsupported values to prevent bogus results
  if(!Jmat.Complex.isReal(p) || !Jmat.Complex.isReal(a) || p.re < 0 || p.re > 1 || a.re < 0) return Jmat.Complex(NaN);

  // TODO: more complete support of the entire complex domain
  return Jmat.Complex.gamma_p_inv_series_1_(a, p);
};

//inverse of gamma_q, but in a instead of in z
//unfortunately very inaccurate. E.g. gamma_q(n, 0.4) returns 0.999... for nearly all n > 4. So anything near that range could give a wide variety of responses.
Jmat.Complex.gamma_q_inva = function(p, z) {
  // This is really not stable and reliable. I need better rootfinding methods and start values.
  return Jmat.Complex.rootfind_bisection(Jmat.Complex(0), Jmat.Complex(100), function(x) {
    return Jmat.Complex.gamma_q(x, z).sub(p);
  });
};

// Calculates gamma(x) / gamma(y), and can cancel out negative integer arguments if both are nonpositive integers in some cases.
// That is useful for functions like beta, binomial, hypergeometric, ...
// Uses the formula gamma(a)/gamma(b) = (-1)^(a-b) * gamma(-b+1) / gamma(-a+1) if necessary
// It is also optimized to use for loops if x and y are nearby integers, ...
Jmat.Complex.gammaDiv_ = function(x, y) {
  var C = Jmat.Complex;
  if(x.eq(y)) return C.ONE; // For the "combined" function, this is considered correct even for negative integers...

  if(C.isInfOrNaN(x) || C.isInfOrNaN(y)) {
    return C(NaN);
  }

  if(C.isNegativeIntOrZero(y) && (x.re > 0 || !C.isInt(x))) return C.ZERO; // division of non-infinity through infinity

  if(C.isInt(x) && C.isInt(y)) {
    if(x.re <= 0 && y.re <= 0) {
      var sign = Jmat.Real.isOdd(x.re - y.re) ? -1 : 1;
      return C.gammaDiv_(y.rsub(1), x.rsub(1)).mulr(sign);
    }

    if(x.re > 0 && y.re > 0 && Jmat.Real.dist(x.re, y.re) < 16) {
      if(x.re > y.re) {
        var result = y.re;
        for(var z = y.re + 1; z < x.re; z++) result *= z;
        return C(result);
      } else {
        var result = 1 / x.re;
        for(var z = x.re + 1; z < y.re; z++) result /= z;
        return C(result);
      }
    }
  }

  return C.gamma(x).div(C.gamma(y));
};

// Similar to Jmat.Complex.gammaDiv_, but takes 3 arguments and cancels out more if possible. Returns gamma(a) / (gamma(b) * gamma(c))
Jmat.Complex.gammaDiv12_ = function(a, b, c) {
  var C = Jmat.Complex;
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(C.isNegativeIntOrZero(b)) {
    if(C.isNegativeIntOrZero(c) && (a.re > 0 || !C.isInt(a))) return C.ZERO; // division of non-infinity through infinity
    return C.gammaDiv_(a, b).div(C.gamma(c));
  } else {
    return C.gammaDiv_(a, c).div(C.gamma(b));
  }
};

// Similar to Jmat.Complex.gammaDiv_, but takes 3 arguments and cancels out more if possible. Returns (gamma(a) * gamma(b)) / gamma(c)
Jmat.Complex.gammaDiv21_ = function(a, b, c) {
  var C = Jmat.Complex;
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(C.isNegativeIntOrZero(a)) {
    return C.gammaDiv_(a, c).mul(C.gamma(b));
  } else {
    return C.gammaDiv_(b, c).mul(C.gamma(a));
  }
};

// Similar to Jmat.Complex.gammaDiv_, but takes 4 arguments and cancels out more if possible. Returns (gamma(a) * gamma(b)) / (gamma(c) * gamma(d))
Jmat.Complex.gammaDiv22_ = function(a, b, c, d) {
  var C = Jmat.Complex;
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(C.isNegativeIntOrZero(a) == C.isNegativeIntOrZero(c)) {
    return C.gammaDiv_(a, c).mul(C.gammaDiv_(b, d));
  } else {
    return C.gammaDiv_(a, d).mul(C.gammaDiv_(b, c));
  }
};

// Calculates log(gamma(x) / gamma(y)) = loggamma(x) - loggamma(y), and can cancel out negative integer arguments if both are negative integers in some cases. That is useful for functions like beta, binomial, ...
// Uses the formula gamma(a)/gamma(b) = (-1)^(a-b) * gamma(-b+1) / gamma(-a+1) if necessary
Jmat.Complex.loggammaDiv_ = function(x, y) {
  var C = Jmat.Complex;
  if(x.eq(y)) return C.ZERO; // For the "combined" function, this is correct even for negative integers...

  if(C.isInfOrNaN(x) || C.isInfOrNaN(y)) {
    return C(NaN);
  }

  if(C.isNegativeIntOrZero(y) && (x.re > 0 || !C.isInt(x))) return C(-Infinity); // division of non-infinity through infinity ==> log(0)

  if(C.isInt(x) && C.isInt(y)) {
    if(x.re <= 0 && y.re <= 0) {
      var sign = Jmat.Real.isOdd(x.re - y.re) ? -1 : 1;
      var result = C.loggammaDiv_(y.rsub(1), x.rsub(1));
      if(sign == -1) result = result.add(C.newi(Math.PI)); // log(-x) = log(x) + i*pi
      return result;
    }
  }

  return C.loggamma(x).sub(C.loggamma(y));
};

// Similar to Jmat.Complex.loggammaDiv_, but takes 3 arguments and cancels out more if possible. Returns loggamma(a) - (loggamma(b) + loggamma(c))
Jmat.Complex.loggammaDiv12_ = function(a, b, c) {
  var C = Jmat.Complex;
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(C.isNegativeIntOrZero(b)) {
    if(C.isNegativeIntOrZero(c) && (a.re > 0 || !C.isInt(a))) return C(-Infinity); // division of non-infinity through infinity ==> log(0)
    return C.loggammaDiv_(a, b).sub(C.loggamma(c));
  } else {
    return C.loggammaDiv_(a, c).sub(C.loggamma(b));
  }
};

// Similar to Jmat.Complex.loggammaDiv_, but takes 3 arguments and cancels out more if possible. Returns (loggamma(a) + loggamma(b)) - loggamma(c)
Jmat.Complex.loggammaDiv21_ = function(a, b, c) {
  var C = Jmat.Complex;
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(C.isNegativeIntOrZero(a)) {
    return C.loggammaDiv_(a, c).add(C.loggamma(b));
  } else {
    return C.loggammaDiv_(b, c).add(C.loggamma(a));
  }
};

// Similar to Jmat.Complex.loggammaDiv_, but takes 4 arguments and cancels out more if possible. Returns (loggamma(a) + loggamma(b)) - (loggamma(c) + loggamma(d))
// To have only 3 values, set a, b, c or d to 1 (it will be fast for that)
Jmat.Complex.loggammaDiv2_ = function(a, b, c, d) {
  var C = Jmat.Complex;
  // Try to combine two negative-integer-ones to have the errors cancel out
  if(C.isNegativeIntOrZero(a) == C.isNegativeIntOrZero(c)) {
    return C.loggammaDiv_(a, c).add(C.loggammaDiv_(b, d));
  } else {
    return C.loggammaDiv_(a, d).add(C.loggammaDiv_(b, c));
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

////////////////////////////////////////////////////////////////////////////////

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

  // METHOD B: in terms of hypergeometric2F1 (computationally more expensive, but precise)
  return x.pow(a).div(a).mul(Jmat.Complex.hypergeometric2F1(a, Jmat.Complex.ONE.sub(b), a.inc(), x));

  // METHOD C: with the integral definition. However this one is numerically very imprecise
  // The integral is numerically too imprecise, due to the degenerate value at t=0
  // the start should be 0, but it's degenerate there
  /*if(Jmat.Complex.isPositive(x) && x.re < 1 && Jmat.Complex.isPositive(a) && Jmat.Complex.isPositive(b)) {
    return Jmat.Complex.integrate(x.divr(3000), x, function(t) {
      var r = t.pow(a.subr(1)).mul(Jmat.Complex.ONE.sub(t).pow(b.subr(1)));
      return r;
    });
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

////////////////////////////////////////////////////////////////////////////////

// The bessel0 and bessel1 functions below are inspired by "Computation of Special Functions" by Shanjie Zhang and Jianming Jin.

//besselj_0 and bessely_0 for large values. Only supports positive z.re.
Jmat.Complex.bessel0big_ = function(z) {
  var C = Jmat.Complex;
  var a = [-0.703125e-1, 0.112152099609375, -0.5725014209747314, 0.6074042001273483,
           -0.1100171402692467e3, 0.03038090510922384e4, -0.1188384262567832e6, 0.06252951493434797e7,
           -0.4259392165047669e9, 0.03646840080706556e11, -0.3833534661393944e13, 0.04854014686852901e15];
  var b = [0.732421875e-1, -0.2271080017089844, 0.1727727502584457e1, -0.2438052969955606e2,
           0.5513358961220206e3, -0.1825775547429318e5, 0.8328593040162893e6, -0.5006958953198893e8,
           0.3836255180230433e10, -0.3649010818849833e12, 0.4218971570284096e14, -0.5827244631566907e16];
  var za = z.abs();
  var num = za >= 50 ? 8 : za >= 35 ? 10 : 12; // number of iterations

  var ca = C.ONE;
  var cb = C(-0.125).div(z);
  var zi = z.inv();
  var zi2 = zi.mul(zi);
  var zz = C.ONE; // zz = z ^ (-2k-2) throughout the loop
  for(var k = 0; k < num; k++) {
    zz = zz.mul(zi2);
    ca = ca.add(zz.mulr(a[k]));
    cb = cb.add(zz.mul(zi).mulr(b[k]));
  }

  var sz = C.sqrt(z.mulr(Math.PI / 2).inv());
  var p = z.subr(Math.PI * 0.25);
  var j0 = sz.mul(ca.mul(C.cos(p)).sub(cb.mul(C.sin(p))));
  var y0 = sz.mul(ca.mul(C.sin(p)).add(cb.mul(C.cos(p))));
  return [j0, y0];
};

// besselj_0 for any complex z
Jmat.Complex.besselj0_ = function(z) {
  var C = Jmat.Complex;
  if(z.eqr(0)) return C(1);
  if(z.re < 0) z = z.neg();

  if(z.abs() < 12) {
    var j0 = C.ONE;
    var r = C.ONE;
    var z2 = z.mul(z);
    for(var k = 1; k <= 40; k++) {
      r = r.mul(z2).mulr(-0.25 / (k * k));
      j0 = j0.add(r);
      if(r.abssq() < j0.abssq() * 1e-30) break;
    }
    return j0;
  } else {
    return Jmat.Complex.bessel0big_(z)[0];
  }
};

// bessely_0 for any complex z
Jmat.Complex.bessely0_ = function(z) {
  var C = Jmat.Complex;
  if(z.eqr(0)) return C(-Infinity);
  var j0, y0;

  var neg = z.re < 0;
  if(neg) z = z.neg();

  if(z.abs() < 12) {
    j0 = C.besselj0_(z);
    var w = 0;
    var r = C.ONE;
    var s = C.ZERO;
    var z2 = z.mul(z);
    for(var k = 1; k <= 40; k++) {
      w = w + 1.0 / k;
      r = r.mulr(-0.25 / (k * k)).mul(z2);
      var cp = r.mulr(w);
      s = s.add(cp);
      if(cp.abssq() < s.abssq() * 1e-30) break;
     }
     var p = C(2 / Math.PI);
     y0 = p.mul(C.log(z.divr(2)).add(C.EM)).mul(j0).sub(p.mul(s));
  } else {
    var jy = Jmat.Complex.bessel0big_(z);
    j0 = jy[0];
    y0 = jy[1];
  }
  if(neg) {
    if(z.im <= 0) y0 = y0.add(C(0, 2).mul(j0));
    else y0 = y0.sub(C(0, 2).mul(j0));
  }

  return y0;
};

//besselj_1 and bessely_1 for large values. Only supports positive z.re.
Jmat.Complex.bessel1big_ = function(z) {
  var C = Jmat.Complex;
  var a = [0.1171875,-0.144195556640625, 0.6765925884246826,-0.6883914268109947e1,
           0.1215978918765359e3, -0.3302272294480852e4, 0.1276412726461746e6, -0.6656367718817688e7,
           0.4502786003050393e9, -0.3833857520742790e11, 0.4011838599133198e13, -0.5060568503314727e15];
  var b = [-0.1025390625,.2775764465332031, -0.1993531733751297e1,0.2724882731126854e2,
           -0.6038440767050702e3, 0.1971837591223663e5, -0.8902978767070678e6, 0.5310411010968522e8,
           -0.4043620325107754e10, 0.3827011346598605e12, -0.4406481417852278e14, 0.6065091351222699e16];
  var za = z.abs();
  var num = za >= 50 ? 8 : za >= 35 ? 10 : 12; // number of iterations

  var ca = C.ONE;
  var cb = C(0.375).div(z);
  var zi = z.inv();
  var zi2 = zi.mul(zi);
  var zz = C.ONE; // zz = z ^ (-2k-2) throughout the loop
  for(var k = 0; k < num; k++) {
    zz = zz.mul(zi2);
    ca = ca.add(zz.mulr(a[k]));
    cb = cb.add(zz.mul(zi).mulr(b[k]));
  }

  var sz = C.sqrt(z.mulr(Math.PI / 2).inv());
  var p = z.subr(Math.PI * 0.75);
  var j1 = sz.mul(ca.mul(C.cos(p)).sub(cb.mul(C.sin(p))));
  var y1 = sz.mul(ca.mul(C.sin(p)).add(cb.mul(C.cos(p))));
  return [j1, y1];
};

// besselj_1 for any complex z
Jmat.Complex.besselj1_ = function(z) {
  var C = Jmat.Complex;
  if(z.eqr(0)) return C(0);
  var j1;

  var neg = z.re < 0;
  if(neg) z = z.neg();

  if(z.abs() < 12) {
    j1 = C.ONE;
    var r = C.ONE;
    var z2 = z.mul(z);
    for(var k = 1; k <= 40; k++) {
      r = r.mul(z2).mulr(-0.25 / (k * (k + 1)));
      j1 = j1.add(r);
      if(r.abssq() < j1.abssq() * 1e-30) break;
    }
    j1 = j1.mul(z).mulr(0.5);
  } else {
    j1 = Jmat.Complex.bessel1big_(z)[0];
  }

  if(neg) j1 = j1.neg();
  return j1;
};

// bessely_1 for any complex z
Jmat.Complex.bessely1_ = function(z) {
  var C = Jmat.Complex;
  if(z.eqr(0)) return C(-Infinity);
  var j1, y1;

  var neg = z.re < 0;
  if(neg) z = z.neg();

  if(z.abs() < 12) {
    j1 = C.besselj1_(z);
    var w = 0;
    var r = C.ONE;
    var s = C.ONE;
    var z2 = z.mul(z);
    for(var k = 1; k <= 40; k++) {
      w = w + 1.0 / k;
      r = r.mulr(-0.25 / (k * (k + 1))).mul(z2);
      var cp = r.mulr(2 * w + 1.0 / (k + 1));
      s = s.add(cp);
      if(cp.abssq() < s.abssq() * 1e-30) break;
     }
     var p = C(2 / Math.PI);
     y1 = p.mul((C.log(z.divr(2)).add(C.EM)).mul(j1).sub(z.inv()).sub(z.mul(s).mulr(0.25)));
  } else {
    var jy = Jmat.Complex.bessel1big_(z);
    j1 = jy[0];
    y1 = jy[1];
  }

  if(neg) {
    if(z.im <= 0) y1 = y1.add(C(0, 2).mul(j1)).neg();
    else y1 = y1.sub(C(0, 2).mul(j1)).neg();
  }

  return y1;
};

// End of bessel0 and bessel1 functions

// This returns sqrt(z), with a different branch for negative real z.
Jmat.Complex.bessel_sqrt_ = function(z) {
  if(Jmat.Complex.isNegative(z)) {
    return Jmat.Complex.sqrt(z).neg();
  } else {
    return Jmat.Complex.sqrt(z);
  }
};

// This returns sqrt(2/(pi*z)), with a different branch for negative real z.
// This is an often occurring value in Bessel formulas.
Jmat.Complex.bessel_sqrt2piz_ = function(z) {
  return Jmat.Complex.bessel_sqrt_(z.rdiv(2/Math.PI));
};

// besselj for integer n >= 2
// To avoid too many loops, do not call for large n. There are probably faster approximations in such zones.
Jmat.Complex.besselj_miller_ = function(n, z) {
  if(z.eqr(0)) return Jmat.Complex(0);
  if(n.eqr(1)) return Jmat.Complex.besselj1_(z);

  // Miller's algorithm

  if(n.re < z.abs() / 4) {
    // Ascending (forwards recurrence)
    var j0 = Jmat.Complex.besselj0_(z);
    var j1 = Jmat.Complex.besselj1_(z);

    for(var i = 1; i < n.re; i++) {
      var j = j1.mulr(2 * i).div(z).sub(j0);
      j0 = j1;
      j1 = j;
    }
    return j1;
  }

  // Descending (backwards recurrence)
  var j0 = Jmat.Complex.besselj0_(z);
  //var j1 = Jmat.Complex.besselj1_(z);

  var jn0 = Jmat.Complex(1);
  var jn1 = Jmat.Complex.ZERO;
  var jn;
  var num = Math.ceil(Math.max(z.abs(), n.abs()));

  for(var i = n.re + num; i > 0; i--) {
    var j = jn0.mulr(2 * i).div(z).sub(jn1); // j = J_(i-1)(z)
    jn1 = jn0;
    jn0 = j;
    if(i - 1 == n.re) jn = j;
  }

  var p = j0.div(jn0); // After backwards reccurrence all the way to j0, calculate ratio of backwards j0 and real j0. Then apply same ratio to backwards jn: that is the real result. NOTE: this also gives the result of all other j_n's seen in the series.
  return jn.mul(p);
};

Jmat.Complex.besselj_series_ = function(nu, z) {
  var C = Jmat.Complex;

  // Works for all complex nu except 0 (0 is already handled above)
  // The gamma functions give NaN if nu is negative integer < -1.
  var negintn = C.isNegativeInt(nu);
  if(negintn) nu = nu.neg();
  var result = C(0);
  var m = 1;
  var gn = C.gamma(nu.inc()); // during the loop: gamma(nu + i + 1)
  for(var i = 0; i < 50; i++) {
    var i1 = C(i + 1);
    var d = C.gamma(i1).mul(gn); //calling gamma(i1) is no problem, all integer solutions in this range are cached
    var term = C(m).div(d).mul(z.divr(2).pow(C(2 * i).add(nu)));
    m = -m;
    result = result.add(term);
    gn = gn.mul(nu.add(i1));
  }
  if(negintn && Jmat.Real.isOdd(nu.re)) result = result.neg(); //besselj(-nu, x) = (-1)^nu * besselj(nu, x)
  return result;
};

Jmat.Complex.besselj_hankelexpansion_ = function(nu, z) {
  var C = Jmat.Complex;
  var pi = Math.PI;
  // The hankel expansion. Works for large z, but does not work well for large nu.
  // J_nu(z) ~ sqrt(2/(pi*z)) * [ cos(w) * SUM_k=0..oo (-1)^k*a_2k(nu)/z^2k - sin(w) * SUM_k=0..oo (-1)^k*a_(2k+1)(nu)/z^(2k+1) ]
  // a_k(nu) = ((4*nu^2 - 1^2) * ... * (4*nu^2 - (2k-1)^2)) / (k! * 8^k)
  // w = z - pi/2 * nu - pi/4

  var result;
  var neg = Math.abs(z.arg()) > 3; //pi - some delta
  if(neg) z = z.neg();

  var ak = Jmat.Complex.ONE; // a_k(nu)
  var zz = C.ONE;
  var w = z.sub(nu.mulr(pi / 2)).subr(pi / 4); // omega
  var s = C.ZERO, c = C.ZERO; // the sum for the sine and the sum for the cosine
  var nu4 = nu.mul(nu).mulr(4); // 4 * nu^2
  var prevt;
  var num = Math.max(4, Math.min(60, Math.abs(nu.re)));
  var t;
  for(var k = 0; k < num; k++) {
    // there are two sums, one uses 2k, the other uses 2k+1. This loop instead alternates between a term of each sum. So k here has the meaning of 2k and 2k+1 in the formula.
    if(k > 0) {
      var ak2 = nu4.subr((2 * k - 1) * (2 * k - 1));
      ak = ak.mul(ak2).divr(8 * k); // divided through the 8^k * k!
      zz = zz.mul(z);
      prevt = t;
    }
    var sign = ((k % 4) < 2) ? 1 : -1;
    t = ak.div(zz).mulr(sign);
    if(t.abssq() < 1e-28) break;

    if(k % 2 == 0)  c = c.add(t);
    else s = s.add(t);
  }

  if(t.abssq() > 1e-2) return C(NaN); // no convergence

  c = C.cos(w).mul(c);
  s = C.sin(w).mul(s);
  result = C.bessel_sqrt2piz_(z).mul(c.sub(s));
  // J_nu(z) = z^nu * (-z)^(-nu) * J_nu(-z) = e^(-nu*pi*i) * J_nu(-z). This alters phase. Note that z was negated above.
  if(neg) result = result.mul(C.exp(nu.muli(pi)));
  return result;
};

// Bessel function in terms of hypergeometrics. The Bessel function happens to exactly use the most difficult parameters for hypergeometric function, so this doesn't extend the range too much.
Jmat.Complex.besselj_hypergeom_ = function(nu, z) {
  var C = Jmat.Complex;

  /*// A: in terms of confluent hypergeometric limit function 0F1.
  // In terms of confluent hypergeometric. This is a last resort, it handles large nu better than the hankel expansion. TODO: use a better high-nu formula. E.g. J(110, 100) is wrong.
  var negnu = 1;
  if(C.isNegativeInt(nu)) {
    if(C.isOdd(nu)) negnu = -1;
    nu = nu.neg(); // Formula does not work for negative nu, but J_-nu(z) = (-1)^nu*J_n(z) for integer nu
  }
  var a = z.mulr(0.5).pow(nu).div(C.gamma(nu.inc()));
  var b = C.hypergeometric0F1(nu.inc(), z.mul(z).mulr(-0.25));
  result = a.mul(b).mulr(negnu);*/

  // B: in terms of confluent hypergeometric 1F1
  var negnu = 1;
  if(C.isNegativeInt(nu)) {
    if(C.isOdd(nu)) negnu = -1;
    nu = nu.neg(); // Formula does not work for negative nu, but J_-nu(z) = (-1)^nu*J_n(z) for integer nu
  }

  if(C.isNegativeInt(nu.addr(0.5))) {
    // The formula does not work if nu is a negative half-integer like -2.5. Not even on wolfram alpha: try besselj(-2.5,21) <--> "z^nu/(2^nu * gamma(nu + 1) * exp(i*z)) * hypergeometric1f1(nu+0.5, 2*nu+1, 2i * z) with z=21 and nu=-2.5"
    nu = (nu.abssq() < 5) ? nu.addr(1e-12) : nu.addr(1e-5); // For now this is the best fix I know. This emulates the limit towards the point with negative integer b of 1F1.
  }

  var h = C.hypergeometric1F1(nu.addr(0.5), nu.mulr(2).inc(), z.muli(2));

  /*var a = z.mulr(0.5).pow(nu).div(C.gamma(nu.inc())).mul(C.exp(z.muli(-1)));
  return a.mul(h).mulr(negnu);*/
  // Written in logarithmic form, otherwise it fails for larger values
  var a = C.log(z.mulr(0.5)).mul(nu).sub(C.loggamma(nu.inc())).add(z.muli(-1));
  return C.exp(a.add(C.log(h))).mulr(negnu);

};

Jmat.Complex.besselj_large_nu_ = function(nu, z) {
  var C = Jmat.Complex;
  // Simple asymptotic expansion which works only if |nu| much bigger than |z|, for nu going to positive infinity
  var a = C.bessel_sqrt2piz_(nu).divr(2);
  var b = C.E.mul(z).div(nu).divr(2);
  return a.mul(b.pow(nu));
}

// Bessel function of the first kind
// Works fine unless both |z| and |nu| are high. Returns NaN if it could not calculate a result with the current algorithms
// TODO: Find fast algorithm that reliably supports high |z| and |nu|, that does not require arbitrary precision arithmetic and does not require too much iterations.
Jmat.Complex.besselj = function(nu, z) {
  var C = Jmat.Complex;

  // Infinities. If any of nu or z is any infinity, the answer is 0
  if(C.isInf(nu)) return C(0);
  if(C.isInf(z)) return C(0);
  // If z is 0, result is 1 for nu=0; NaN for nu strict imaginary; 0 for nu.re > 0; undirected infinity for nu.re < 0.
  if(z.eqr(0)) return (nu.re == 0) ? (nu.im == 0 ? C(1) : C(NaN)) : (nu.re < 0 ? C(Infinity, Infinity) : C(0));

  if(nu.eqr(0)) return C.besselj0_(z);
  if(nu.eqr(1)) return C.besselj1_(z);
  if(nu.eqr(0.5)) return C.bessel_sqrt2piz_(z).mul(C.sin(z));
  if(nu.eqr(-0.5)) return C.bessel_sqrt2piz_(z).mul(C.cos(z));

  if(nu.re > 300 && Math.abs(nu.im) < nu.re && z.abssq() < 50*50) return C(0); // underflow, 0 is approximate. Also, for integer n would cause too much loops in miller algorithm
  if(nu.re < -300 && C.isInt(nu) && z.abssq() < 50*50) return C(0); //also 0, however for negative nu if there is the slightest deviation of integer, it becomes huge result instead

  // TODO: check if this function can work for complex z and/or nu. It performs bad for Jmat.besselj('112.5', '-5.5i'), besselj_series_ handles that one much better.
  if(nu.re > 50 && C.isReal(z) && nu.re > z.abs() * 16) {
    return Jmat.Complex.besselj_large_nu_(nu, z);
  }

  if(C.isInt(nu) && Math.abs(nu.re) < 50) {
    var result;
    // Supports integer nu not equal to 0 or 1
    //J_nu(z) = (-1)^nu * J_-nu(z), for integer nu, and besselj_miller_ only supports nu >= 2
    if(nu.re < 0) result = C.besselj_miller_(nu.neg(), z).mulr(C.isOdd(nu) ? -1 : 1);
    else result = C.besselj_miller_(nu, z);
    return result;
  } else if(z.abs() < 25) {
    return C.besselj_series_(nu, z);
  } else {
    if(z.abssq() > nu.abssq()) {
      return C.besselj_hankelexpansion_(nu, z);
    } else if(nu.re > z.abs() * 2) {
      return Jmat.Complex.besselj_large_nu_(nu, z);
    } else {
      // last resort. Will return NaN if that too could not approximate the result
      return C.besselj_hypergeom_(nu, z);
    }
  }
};

// bessely for integer n >= 2
// To avoid too many loops, do not call for large n. There are probably faster approximations in such zones.
Jmat.Complex.bessely_miller_ = function(n, z) {
  var y0 = Jmat.Complex.bessely0_(z);
  var y1 = Jmat.Complex.bessely1_(z);

  // Miller's algorithm
  for(var i = 1; i < n.re; i++) {
    var y = y1.mulr(2 * i).div(z).sub(y0);
    y0 = y1;
    y1 = y;
  }
  return y1;
};

Jmat.Complex.bessely_with_besselj_ = function(nu, z) {
  var C = Jmat.Complex;

  if(C.isInt(nu.addr(0.5))) {
    // For negative half-integers, the cos is zero (cos(pi/2)), but due to numerical imprecision it isn't, and when multiplied by the huge 'a' that gives a completely incorrect result.
    // It's not even necessary to calculate a, and the sin is 1 or -1 so that one neither.
    return C.besselj(nu.neg(), z).mulr(C.isOdd(nu.addr(0.5)) ? -1 : 1);
  }

  // Supports all complex nu except integers (and negative half-integers due to numerical problems), for small z
  // Y_a(x) = (J_a(x)*cos(a*pi) - J_-a(x)) / sin(a*pi)
  // (For integer nu, Y_n(x) = lim_a_to_n(Y_a(x), but that is not supported here)
  var a = C.besselj(nu, z);
  var b = C.cos(nu.mulr(Math.PI));
  if(C.isInt(b.addr(0.5))) b = C(0);
  var c = C.besselj(nu.neg(), z);
  var d = C.sin(nu.mulr(Math.PI));
  return a.mul(b).sub(c).div(d);
};

Jmat.Complex.bessely_hankelexpansion_ = function(nu, z) {
  var C = Jmat.Complex;
  var pi = Math.PI;
  // The hankel expansion
  // TODO: remove some of this code duplication with besselJ hankel expansion
  // TODO: for large nu, still another formula is needed, this one gives wrong values for that
  // J_nu(z) ~ sqrt(2/(pi*z)) * [ sin(w) * SUM_k=0..oo (-1)^k*a_2k(nu)/z^2k + cos(w) * SUM_k=0..oo (-1)^k*a_(2k+1)(nu)/z^(2k+1) ]
  // a_k(nu) = ((4*nu^2 - 1^2) * ... * (4*nu^2 - (2k-1)^2)) / (k! * 8^k)
  // w = z - pi/2 * nu - pi/4
  var neg = Math.abs(z.arg()) > 3; //pi - some delta
  if(neg) z = z.neg();

  var ak = Jmat.Complex.ONE; // a_k(nu)
  var zz = C.ONE;
  var w = z.sub(nu.mulr(pi / 2)).subr(pi / 4); // omega
  var s = C.ZERO, c = C.ZERO; // the sum for the sine and the sum for the cosine
  var nu4 = nu.mul(nu).mulr(4); // 4 * nu^2
  var prevt;
  var num = Math.max(4, Math.min(60, Math.abs(nu.re)));
  var t;
  for(var k = 0; k < num; k++) {
    // there are two sums, one uses 2k, the other uses 2k+1. This loop instead alternates between a term of each sum. So k here has the meaning of 2k and 2k+1 in the formula.
    if(k > 0) {
      var ak2 = nu4.subr((2 * k - 1) * (2 * k - 1));
      ak = ak.mul(ak2).divr(8 * k); // divided through the 8^k * k!
      zz = zz.mul(z);
      prevt = t;
    }
    var sign = ((k % 4) < 2) ? 1 : -1;
    t = ak.div(zz).mulr(sign);
    if(t.abssq() < 1e-28) break;

    if(k % 2 == 0)  c = c.add(t);
    else s = s.add(t);
  }

  if(t.abssq() > 1e-2) return C(NaN); // no convergence

  c = C.sin(w).mul(c);
  s = C.cos(w).mul(s);
  var result = C.bessel_sqrt2piz_(z).mul(c.add(s));
  if(neg) {
    // J_nu(z) = e^(-nu*pi*i) * Y_nu(-z) + i*z*cos(pi*nu) * J_nu(z)
    var j = C.besselj(nu, z);
    result = result.mul(C.exp(nu.muli(-pi))).add(C.cos(nu.mulr(pi)).mul(j).muli(2));
  }
  return result;
};

// Bessel function of the second kind
// Mostly an approximation, there are problems with high z
Jmat.Complex.bessely = function(nu, z) {
  var C = Jmat.Complex;

  // Special values
  if(C.isInf(nu)) return C.isInf(z) ? C(NaN) : (C.isPositive(z) ? C(-Infinity) : C(Infinity, Infinity)); // strictly positive z gives -Infinity, anything else, including 0 and complex values, gives complex infinity
  if(C.isInf(z)) return C(0);
  if(z.eqr(0)) return (nu.re == 0) ? (nu.im == 0 ? C(-Infinity) : C(NaN)) : C(Infinity, Infinity); // If z is 0, result is -inf for nu=0; NaN for nu strict imaginary (in some software); undirected infinity for other real or complex nu.

  if(nu.eqr(0)) return C.bessely0_(z);
  if(nu.eqr(1)) return C.bessely1_(z);

  if(C.isInt(nu)) {
    // Supports integer nu not equal to 0 or 1
    var neg = z.re < 0;
    if(neg) z = z.neg();
    var result;

    //Y_n(z) = (-1)^nu * Y_-nu(z), for integer nu, and bessely_miller_ only supports nu >= 2
    if(nu.re < 0) result = C.bessely_miller_(nu.neg(), z).mulr(C.isOdd(nu) ? -1 : 1);
    else result = C.bessely_miller_(nu, z);

    if(neg) {
      // For integer nu: Y_n(z) = (-1)^nu * (Y_n(-z) +- i*2*J_n(-z)) --> the +- is sign of z.im, or + if real
      var j = C.besselj(nu, z);
      if(z.im <= 0) result = result.add(C(0, 2).mul(j));
      else result = result.sub(C(0, 2).mul(j));
      if(C.isOdd(nu)) result = result.neg();
    }

    return result;
  } else if(z.abs() < 20) {
    return C.bessely_with_besselj_(nu, z);
  } else if(z.abssq() > nu.abssq()) {
    return C.bessely_hankelexpansion_(nu, z);
  } else {
    return C(NaN); // unfortunately not supported. TODO: find stable fast algorithm that does not require arbitrary precision that can solve bessely for large |z| and large |nu|
  }
};

// Hankel function of the first kind (approximation)
Jmat.Complex.hankelh1 = function(nu, z) {
  var C = Jmat.Complex;

  // Special values
  if(C.isInf(nu)) return C(Infinity, Infinity);
  if(C.isInf(z)) return C(0);
  if(z.eqr(0)) return (nu.re == 0 && nu.im != 0) ? C(NaN) : C(Infinity, Infinity); // undirected infinity, or nan for strict imaginary nu

  // There are numerical imprecisions for e.g. hankelh2(-8-8i, i).
  // So when needed, apply the formula: H1_-a(x) = exp(a*pi*i)*H1_a(x)
  if(nu.im < 0) return C.exp(nu.mul(C.newi(-Math.PI))).mul(C.hankelh1(nu.neg(), z));

  return C.besselj(nu, z).add(C.bessely(nu, z).mul(C.I));
};

// Hankel function of the second kind (approximation)
Jmat.Complex.hankelh2 = function(nu, z) {
  var C = Jmat.Complex;

  // Special values
  if(C.isInf(nu)) return C(Infinity, Infinity);
  if(C.isInf(z)) return C(0);
  if(z.eqr(0)) return (nu.re == 0 && nu.im != 0) ? C(NaN) : C(Infinity, Infinity); // undirected infinity, or nan for strict imaginary nu

  // There are numerical imprecisions for e.g. hankelh2(-8-8i, i).
  // So when needed, apply the formula: H2_-a(x) = exp(-a*pi*i)*H2_a(x)
  if(nu.im > 0) return C.exp(nu.mul(C.newi(Math.PI))).mul(C.hankelh2(nu.neg(), z));

  return C.besselj(nu, z).sub(C.bessely(nu, z).mul(C.I));
};

// Modified bessel function of the first kind (approximation)
Jmat.Complex.besseli = function(nu, z) {
  var C = Jmat.Complex;

  // Special values
  if(C.isInf(nu)) return C.isInf(z) ? C(NaN) : C(0);
  if(z.eqr(Infinity)) return C(Infinity);
  if(C.isInf(z)) return C(Infinity, Infinity); // other infinities for z
  // If z is 0, result is 1 for nu=0; NaN for nu strict imaginary; 0 for nu.re > 0; undirected infinity for nu.re < 0.
  if(z.eqr(0)) return (nu.re == 0) ? (nu.im == 0 ? C(1) : C(NaN)) : (nu.re < 0 ? C(Infinity, Infinity) : C(0));

  var result = C.I.pow(nu.neg()).mul(C.besselj(nu, z.mul(C.I)));
  if(z.im == 0 && nu.im == 0 && Jmat.Real.near(result.im, 0, 1e-10)) result.im = 0;
  return result;
};

// Modified bessel function of the second kind (approximation)
Jmat.Complex.besselk = function(nu, z) {
  var C = Jmat.Complex;

  // Special values
  if(C.isInf(nu)) return C.isInf(z) ? C(NaN) : (C.isPositive(z) ? C(Infinity) : C(Infinity, Infinity)); // strictly positive z gives Infinity, anything else, including 0 and complex values, gives complex infinity
  if(z.eqr(Infinity)) return C(0);
  if(C.isInf(nu)) return C(Infinity, Infinity); // other infinities for z
  if(z.eqr(0)) return (nu.re == 0) ? (nu.im == 0 ? C(Infinity) : C(NaN)) : C(Infinity, Infinity); // If z is 0, result is inf for nu=0; NaN for nu strict imaginary (in some software); undirected infinity for other real or complex nu.

  var result;
  if(z.im >= 0) {
    result = C.I.pow(nu.inc()).mulr(Math.PI / 2).mul(C.hankelh1(nu, C.I.mul(z)));
  } else {
    result = C.newi(-1).pow(nu.inc()).mulr(Math.PI / 2).mul(C.hankelh2(nu, C.newi(-1).mul(z)));
  }
  if(z.im == 0 && nu.im == 0 && Jmat.Real.near(result.im, 0, 1e-3 /*this one is super imprecise...*/)) result.im = 0;
  return result;
};

Jmat.Complex.spherical_besselj = function(nu, z) {
  var C = Jmat.Complex;
  var factor = C.sqrt(C(Math.PI / 2).div(z));
  var j = C.besselj(nu.addr(0.5), z);
  return factor.mul(j);
};

Jmat.Complex.spherical_bessely = function(nu, z) {
  var C = Jmat.Complex;
  var factor = C.sqrt(C(Math.PI / 2).div(z));
  var j = C.bessely(nu.addr(0.5), z);
  return factor.mul(j);
};

Jmat.Complex.spherical_hankelh1 = function(nu, z) {
  var C = Jmat.Complex;
  var factor = C.sqrt(C(Math.PI / 2).div(z));
  var j = C.hankelh1(nu.addr(0.5), z);
  return factor.mul(j);
};

Jmat.Complex.spherical_hankelh2 = function(nu, z) {
  var C = Jmat.Complex;
  var factor = C.sqrt(C(Math.PI / 2).div(z));
  var j = C.hankelh2(nu.addr(0.5), z);
  return factor.mul(j);
};

//Struve function H_nu(z)
//NOTE: has errors for large z (|z| > 20), at least there are regions where it's really bad
Jmat.Complex.struveh = function(nu, z) {
  var C = Jmat.Complex;

  if(C.abs(z) > 30 && C.abs(z) < 50) {
    // Test if this works with struveh(10, 31) = 26245.845183638363381068944066845095005243622259212620
    var factor = z.divr(2).pow(nu.addr(1));
    var f = C.regularized_hypergeometric([C(1)], [C(1.5), nu.addr(1.5)], z.mul(z).divr(-4));
    return factor.mul(f);
  }

  if(C.abs(z) > 30) {
    var g1mul = C(0.5);
    var g1 = C.gamma(g1mul); // gamma(k + 1/2)
    var g2mul = nu.addr(0.5);
    var g2 = C.gamma(g2mul); // gamma(nu + 1/2 - k)
    var zn = C(1); // (z/2)^(2k)
    var z2 = z.mul(z).divr(4);

    var result = C(0);
    var prev = result;
    for(var n = 0; n < 30; n++) {
      result = result.add(g1.div(g2).div(zn));
      if(result.eq(prev)) break;
      prev = result;
      g1 = g1.mul(g1mul);
      g1mul = g1mul.addr(1);
      g2mul = g2mul.subr(1);
      g2 = g2.div(g2mul);
      zn = zn.mul(z2);
    }

    var factor = z.mulr(0.5).pow(nu.subr(1)).divr(Math.PI);

    return factor.mul(result);
  }

  var g1mul = C(1.5);
  var g1 = C.gamma(g1mul); // gamma(n + 3/2)
  var g2mul = C(nu.addr(1.5));
  var g2 = C.gamma(g2mul); // gamma(n + nu + 3/2)
  var m = C(1); // (-1)^n
  var zn = C(1); // (z/2)^(2*n)
  var z2 = z.mul(z).divr(4);

  var result = C(0);
  var prev = result;
  for(var n = 0; n < 100; n++) {
    result = result.add(m.mul(zn).div(g1).div(g2));
    if(result.eq(prev)) break;
    prev = result;
    g1 = g1.mul(g1mul);
    g2 = g2.mul(g2mul);
    g1mul = g1mul.addr(1);
    g2mul = g2mul.addr(1);
    m = m.neg();
    zn = zn.mul(z2);
  }

  var factor = z.mulr(0.5).pow(nu.addr(1));

  return factor.mul(result);
};

//Struve function K_nu(z)
Jmat.Complex.struvek = function(nu, z) {
  var C = Jmat.Complex;
  return C.struveh(nu, z).sub(C.bessely(nu, z));
};

//Modified Struve function L_nu(z)
Jmat.Complex.struvel = function(nu, z) {
  var C = Jmat.Complex;
  var factor = C.exp(nu.muli(-Math.PI * 0.5)).muli(-1);
  return factor.mul(C.struveh(nu, z.muli(1)));
};

//Modified Struve function M_nu(z)
Jmat.Complex.struvem = function(nu, z) {
  var C = Jmat.Complex;
  return C.struvel(nu, z).sub(C.besseli(nu, z));
};

// Anger's J function
Jmat.Complex.angerj = function(nu, z) {
  var C = Jmat.Complex;
  return C.integrate(C(0), C(Math.PI), function(theta) {
    return C.cos(nu.mul(theta).sub(z.mul(C.sin(theta))));
  }).divr(Math.PI);
};

// Weber's E function
Jmat.Complex.webere = function(nu, z) {
  var C = Jmat.Complex;
  return C.integrate(C(0), C(Math.PI), function(theta) {
    return C.sin(nu.mul(theta).sub(z.mul(C.sin(theta))));
  }).divr(Math.PI);
};

////////////////////////////////////////////////////////////////////////////////

//pl = 3^(-2/3) for airy, 3^(-4/3) for bairy
//pr = 3^(-4/3) for airy, 3^(-5/6) for bairy
//s = -1 for airy, +1 for bairy
Jmat.Complex.airyloop_ = function(z, pl, pr, s) {
  var gl = 1.3541179394264004169; // GAMMA(2/3)
  var gr = 0.8929795115692492112; // GAMMA(4/3)
  var zzz = z.mul(z).mul(z);
  var zr = Jmat.Complex.ONE;
  var zl = z;
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
  if(z.abs() > 8) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.arg()) < 2 * Math.PI / 3) {
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
  if(z.abs() > 10) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.arg()) < Math.PI / 3) {
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
  if(z.abs() > 8) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.arg()) < 2 * Math.PI / 3) {
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
  if(z.abs() > 10) {
    // asymptotic expansion for large |z|
    // NOTE: only uses the first term. TODO use more?
    if(Math.abs(z.arg()) < Math.PI / 3) {
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


////////////////////////////////////////////////////////////////////////////////

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
    if(q1.abssq() < 1e-4 && q1.abssq() < 1e-4) break;
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

  var a1 = s.dec().abs();
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
for n = 0 to steps: 1/(n+1) SUM_k=0..n ((-1)^k binomial(n, k)
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
  return Jmat.Complex(NaN); // not supported in this region
};

////////////////////////////////////////////////////////////////////////////////


Jmat.Complex.lerchphi_pow_ = function(z, s, opt_var) {
  if(opt_var) return z.mul(z).pow(s.divr(2));
  else return z.pow(s);
};

// This is slower than Jmat.Complex.hurwitzzeta_generate_hasse_table_, because this depends on a and s, so it is only
// cached for constant a, due to the reason that in the lerchphi series with binomial, some value depends on z, while for
// hurwitz zeta that is not the case.
Jmat.Complex.lerchphi_generate_binomial_table_ = function(s, a, steps, opt_var) {
  var C = Jmat.Complex;
  var result = [];
  for(var n = 0; n < steps; n++) {
    result[n] = C(0);
    var sign = 1;
    for(var k = 0; k <= n; k++) {
      var p = C.lerchphi_pow_(a.addr(k), s, opt_var);
      result[n] = result[n].add(p.inv().mulr(sign * Jmat.Real.pascal_triangle(n, k)));
      sign = -sign;
    }
  }
  return result;
};

Jmat.Complex.lerchphi_binomial_table_a_ = undefined;
Jmat.Complex.lerchphi_binomial_table_s_ = undefined;
Jmat.Complex.lerchphi_binomial_table_var_ = undefined;
Jmat.Complex.lerchphi_binomial_table_ = undefined;

Jmat.Complex.lerchphi_binomial_series_ = function(z, s, a, opt_var) {
  var C = Jmat.Complex;
  var N = 30;
  if(!C.lerchphi_binomial_table_ || !C.lerchphi_binomial_table_a_.eq(a) || !C.lerchphi_binomial_table_s_.eq(s) || C.lerchphi_binomial_table_var_ != !!opt_var) {
    C.lerchphi_binomial_table_ = C.lerchphi_generate_binomial_table_(s, a, N);
    C.lerchphi_binomial_table_a_ = a;
    C.lerchphi_binomial_table_s_ = s;
    C.lerchphi_binomial_table_var_ = !!opt_var;
  }
  var z2 = z.neg().div(z.rsub(1));
  var table = C.lerchphi_binomial_table_;
  var result = C.ZERO;
  var zz = C(1);
  for(var k = 0; k < N; k++) {
    var r = result.add(zz.mul(table[k]));
    if(C.near(r, result, 1e-15)) break;
    result = r;
    zz = zz.mul(z2);
  }
  result = result.mul(z.rsub(1).inv());
  return result;
};

Jmat.Complex.lerchphi_series_ = function(z, s, a, opt_var) {
  var C = Jmat.Complex;
  var result = C(0);
  var zz = C(1);
  for(var k = 0; k < 30; k++) {
    var p = C.lerchphi_pow_(a.addr(k), s, opt_var);
    var r = zz.div(p);
    if(C.near(result.add(r), result, 1e-15)) break; // converged
    result = result.add(r);
    zz = zz.mul(z); // z^k
  }
  return result;
};


//NOTE: the precision of this is low
Jmat.Complex.lerchphi_integral_ = function(z, s, a, opt_var) {
  var C = Jmat.Complex;

  if(a.re < 1) {
    if(C.isNegativeIntOrZero(a)) return C(Infinity, Infinity);
    var m = Math.ceil(1 - a.re);
    if(m > 1000) return C(NaN); // a bit too much
    var v = C(0);
    var zz = C(1);
    for(var k = 0; k < m; k++) {
      v = v.add(zz.div(C.lerchphi_pow_(a.addr(k), s)));
      zz = zz.mul(z);
    }
    return zz.mul(C.lerchphi(z, s, a.addr(m), opt_var)).add(v);
  }

  // The code below requires a.re >= 1.
  // No lerchphi_pow_ necessary: only makes difference for negative a.re

  if(z.eqr(1)) z = C(1.0000000000001); // avoid singularity. TODO: do this better (this introduces small imaginary part): find out how incgamma multiplied by the log cancel out.
  var l = C.log(z);
  var result = a.pow(s).mulr(2).inv();
  var g = C.incgamma_upper(s.rsub(1), a.mul(l).neg());
  result = result.add(g.mul(l.neg().pow(s.subr(1))).div(z.pow(a)));
  var f = function(t) {
    var u = C.sin(C.atan(t.div(a)).mul(s).sub(t.mul(l)));
    var v = a.mul(a).add(t.mul(t)).pow(s.divr(2));
    var w = C.expm1(t.mulr(Math.PI * 2));
    var r = u.div(v).div(w);
    return r;
  };
  var q = C.integrate(C(0), C(Infinity), f, 30, 1).mulr(2);
  result = result.add(q);
  return result;
};

// Lerch transcendent. Precise if |z| < 0.75, numerically imprecise in many other cases
// opt_var enables a variation (different branch cut), differing only for a.re < 0:
// Without opt_var, it is SUM_k=0..oo z^k / (a+k)^s. This is the default, but is sometimes called hurwitzhlerchphi
// With opt_var, it is SUM_k=0..oo z^k / ((a+k)^2)^(s/2)
// TODO: with opt_var true, then with a=0 it should returns some result, not infinity (and it's not just a limit)
Jmat.Complex.lerchphi = function(z, s, a, opt_var) {
  var C = Jmat.Complex;
  if(z.abs() < 0.75) return Jmat.Complex.lerchphi_series_(z, s, a, opt_var);
  //if(z.re < 0.5) return Jmat.Complex.lerchphi_binomial_series_(z, s, a, opt_var); //the integral is more precise than this one, disabled
  return C.lerchphi_integral_(z, s, a, opt_var);
};


////////////////////////////////////////////////////////////////////////////////

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

  var a = z.abs();
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
  var a1 = z1.abs();
  var z2 = Jmat.Complex.ONE.sub(z).inv();
  var a2 = z2.abs();
  var z3 = z.div(z.dec());
  var a3 = z3.abs();

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

  var a = z.abs();
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
    -7709321041217/870, 0, 2577687858367/14322, 0 //28-31
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
  0.00640006853170062946, 0.00737115177047223913, 0.00355772885557316095, -0.00751332599781522893
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
  2.09915158936342553e-32, -8.33674529544144048e-34, 1.34125937721921867e-35, 9.13714389129817200e-37
];

// This is really a last resort, very inaccurate and slow
// Welcome to the land of randomness and bad results.
// Some cases are really good, others horrible.
// Unfortunately almost every integral has problems on the positive real axis of z.
// Since there is an awesome solution for all s with s.re < 0, only the formulas for s.re > 0 below are actually u sed
Jmat.Complex.polylog_integral_ = function(s, z) {
  var C = Jmat.Complex;
  // To test these, try e.g.:
  // complexDomainPlot(function(z){return C.polylog(C(15, 0.5), z);}, 2, 1);
  // complexDomainPlot(function(z){return C.polylog(C(0.5, 15), z);}, 2, 1);

  if(s.re > 1 && Math.abs(s.im) < s.re && Math.abs(z.arg()) > 0.1) {
    // Approximate with an integral representation (only works for real s.re, and has some serious problems with periodic things appearing near the positive real axis of z)
    // In practice, only works for s.re > 1 (!!) and s.im = 0, or very small |s.im| and s.re > 0
    var g = C.gamma(s);
    var f = function(t) {
      var result = t.pow(s.dec()).div(C.exp(t).sub(z));
      if(C.isNaN(result)) result = C.ZERO;
      return result;
    };
    var r = C.integrate(C(0), C(Infinity), f, 50);
    return z.div(g).mul(r);
  } else if(C.isNegative(s) && Math.abs(z.arg()) > 0.1) {
    // Approximate with an integral representation (only works for s.re < 0
    var lzm = C.log(z.neg());
    var f = function(t) {
      var ta = t.pow(s.neg());
      var tb = C.sin(s.mulr(Math.PI/2).sub(t.mul(lzm)));
      var na = C.sinh(t.mulr(Math.PI));
      var result = ta.mul(tb).div(na);
      if(C.isNaN(result)) result = C.ZERO;
      return result;
    };
    var r = C.integrate(C(0), C(Infinity), f, 30);
    return r;
  } else if(z.im <= 0 || (s.re > 0 && Math.abs(s.im) < s.re)) {
    // Integral formula for polylog that works for all complex s and z, except for positive real z if s.re < 0
    // Is from 0 to infinity, but seems to work well from 0-10 with only 100 steps (at least in the areas where this is used)
    // While in theory it works for all z and all s (except z near positive axis if s.re < 0), in practice it seems to work only for z.im <= 0 and s.re > 0 and |s.im| << s.re
    // --> I tried with z-plot for s=-15+0.5i, s=0.5+15i, and other variations like that
    var lzm = C.log(z.neg());
    var f = function(t) {
      var ta = C.sin(s.mul(C.atan(t)).sub(t.mul(lzm)));
      var na = (C.ONE.add(t.mul(t))).pow(s.divr(2));
      var nb = C.sinh(t.mulr(Math.PI));
      var result = ta.div(na).div(nb);
      if(C.isNaN(result)) result = C.ZERO;
      return result;
    };
    // TODO: better integral. The interval is in reality from 0 to infinity, but the gaus laguerre implementation (or f) here returns a better result if you end it at 10 than at infinity...
    var r = C.integrate(C(0), C(10), f, 30);
    return z.mulr(0.5).add(z.mul(r));
  } else if(z.im > 0) { // Because something is broken for z.im <= 0. TODO: find out what
    // Similar integral, but with upper incomplete gamma.
    // It is theoretically better than the above because it supports even its last edge case.
    // But it seems broken. It does not work for z.im <= 0... At least it gets
    var lz = C.log(z);
    var f = function(t) {
      var ta = C.sin(s.mul(C.atan(t)).sub(t.mul(lz)));
      var na = (C.ONE.add(t.mul(t))).pow(s.divr(2));
      var nb = C.exp(t.mulr(2*Math.PI)).subr(1);
      var result = ta.div(na).div(nb);
      if(C.isNaN(result)) result = C.ZERO;
      return result;
    };
    var r = C.integrate(C(0), C(Infinity), f, 30);
    var g = C.incgamma_upper(C.ONE.sub(s), lz.neg());
    var l = lz.neg().pow(C.ONE.sub(s));
    return z.mulr(0.5).add(g.div(l)).add(z.mulr(2).mul(r));
  } else {
    return C(NaN); //there is nothing we can do... :(
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
  var kidney_radius = z.mul(z).div(z.dec()).abs(); //radius of the kidney shaped region of convergence. Theoretically works for < 4, in practice with double precision only for < 3 and then only if not too much n!
  if(kidney_radius >= 3.7) return Jmat.Complex(NaN); //yeah right... it sucks way before this

  // number of loops. NOTE: higher is better for arbitrary precision library (with 31 being perfect), but results in random garbage with floating point. So it is limited here instead.
  var n = Math.floor(Math.min(31, 16 / kidney_radius));


  var binomial_cache = [];
  binomial_cache[0] = Jmat.Complex.ONE;
  var bin_sum_term = Jmat.Complex.ONE; // Binomial summation term

  var zz = Jmat.Complex.ONE;
  var result0 = Jmat.Complex.ZERO;

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
  var kidney_radius = z.mul(z).div(z.dec()).abs();
  return (s.re >= -5 && z.im == 0 && kidney_radius <= 1.5) || (s.re >= 0 && Math.abs(s.im) < s.re && kidney_radius <=2) || z.abs() <= 0.5;
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
  /*var square = function(s, z) {
    var a = Jmat.Complex.polylog(s, z.mul(z));
    var b = Jmat.Complex.polylog(s, z.neg());
    var c = Jmat.Complex.TWO.pow(Jmat.Complex.ONE.sub(s));
    return c.mul(a).sub(b);
  };*/

  var a = z.abs();

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
    // TODO: the integral with gaussian laguerre handles this one better than polylog_residue_, either fix polylog_residue_ or use the integral
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
};

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

// Returns a particular non-principal complex branch of sqrt(a * b)
// For elliptic integrals, and AGM and GHM, a particular branch of the complex sqrt must be returned, otherwise it hops to other branches and gives wrong result
Jmat.Complex.mulSqrtBranch_ = function(a, b) {
  // Note: if a and b are real but negative, with a*b being positive, still a different branch is used.
  if(Jmat.Complex.isPositive(a) && Jmat.Complex.isPositive(b)) {
    return Jmat.Complex.sqrt(a.mul(b));
  } else {
    // This gives a different result than Jmat.Complex.sqrt(a.mul(b))
    return Jmat.Complex.sqrt(a).mul(Jmat.Complex.sqrt(b));
    //the below gives the same result, but is slower
    //return Jmat.Complex(Math.sqrt(a.mul(b).abs())).mul(Jmat.Complex.exp(Jmat.Complex.I.mulr((a.arg() + b.arg()) / 2)));
  }
};

//Elliptic integral: Carlson's form R_F
Jmat.Complex.ellipticrf = function(x, y, z) {
  // Homogeneity: R_F(kx, ky, kz) = k^(-0.5) * R_F(x, y, z)
  // Duplication theorem: R_F(x, y, z) = R_F((x+l)/4, (y+l)/4, (z+l)/4) with l the lambda in the algorithm below.
  // The algorithm is based on the duplication theorem.
  var C = Jmat.Complex;
  var count = 0;
  var ex, ey, ez;
  var a;
  for(;;) {
    a = x.add(y).add(z).divr(3);
    if(C.nearr(a, 0, 1e-14)) return C(0); // avoid NaNs, result is 0 in this case
    ex = x.div(a).rsub(1);
    ey = y.div(a).rsub(1);
    ez = z.div(a).rsub(1);
    var error = Math.max(Math.max(ex.abs(), ey.abs()), ez.abs());
    if(error < 0.0025) break; // converged
    var sx = C.sqrt(x);
    var sy = C.sqrt(y);
    var sz = C.sqrt(z);
    // The square roots must be calculated as sqrt(x)*sqrt(y), not as sqrt(x*y), to ensure correct branch (given that sqrt returns the solution with positive real part) (see mulSqrtBranch_)
    var lambda = sx.mul(sy).add(sy.mul(sz)).add(sz.mul(sx));
    x = x.add(lambda).divr(4);
    y = y.add(lambda).divr(4);
    z = z.add(lambda).divr(4);
    count++;
    if(count > 1000) break; // avoid infinite loop
  }
  var e2 = ex.mul(ey).sub(ez.mul(ez));
  var e3 = ex.mul(ey).mul(ez);
  // extra error term (should be near 1)
  var extra = C(1).sub(e2.divr(10)).add(e3.divr(14)).add(e2.mul(e2).divr(24)).sub(e2.mul(e3).mulr(3.0 / 44));
  return a.powr(-0.5).mul(extra);
};

//Elliptic integral: Carlson's form R_C
Jmat.Complex.ellipticrc = function(x, y) {
  // R_C(x, y) == R_F(x, y, y), but the implementation below is more efficient
  var C = Jmat.Complex;
  var w;
  if(y.re > 0) {
    w = C(1);
  } else {
    // Cauchy principal value
    var xt = x.sub(y);
    w = C.sqrt(x).div(C.sqrt(xt));
    x = xt;
    y = y.neg();
  }
  var count = 0;
  var ey;
  var a;
  for(;;) {
    a = x.add(y).add(y).divr(3);
    ey = y.div(a).subr(1);
    if(ey.abs() < 0.0012) break;
    // The square roots must be calculated as sqrt(x)*sqrt(y), not as sqrt(x*y), to ensure correct branch (given that sqrt returns the solution with positive real part) (see mulSqrtBranch_)
    var lambda = C.sqrt(x).mul(C.sqrt(y)).mulr(2).add(y);
    x = x.add(lambda).divr(4);
    y = y.add(lambda).divr(4);
    count++;
    if(count > 1000) break; // avoid infinite loop
  }
  // extra error term (should be near 1)
  var ey2 = ey.mul(ey);
  var ey4 = ey2.mul(ey2);
  var ey6 = ey4.mul(ey2);
  var extra = C(1).add(ey2.mulr(3.0/10)).add(ey2.mul(ey).divr(7)).add(ey4.mulr(3.0/8)).add(ey4.mul(ey).mulr(9.0/22)).add(ey6.mulr(159.0/208)).add(ey6.mul(ey).mulr(9.0/8));
  return a.powr(-0.5).mul(w).mul(extra);
};

//Elliptic integral: Carlson's form R_J
Jmat.Complex.ellipticrj = function(x, y, z, p) {
  var C = Jmat.Complex;
  // Homogeneity: R_F(kx, ky, kz) = k^(-0.5) * R_F(x, y, z)
  // Duplication theorem: R_F(x, y, z) = R_F((x+l)/4, (y+l)/4, (z+l)/4) with l the lambda in the algorithm below.
  // The algorithm is based on the duplication theorem.

  var origp = p;
  var ca, cb, rcx;
  if(origp.re <= 0) {
    // Cauchy principal value
    if(x.re > y.re) { var t = x; x = y; y = t; }
    if(y.re > z.re) { var t = y; y = z; z = t; }
    if(x.re > y.re) { var t = x; x = y; y = t; }
    ca = y.sub(origp).rdiv(1);
    cb = ca.mul(z.sub(y)).mul(y.sub(x));
    p = y.add(cb);
    rcx = C.ellipticrc(x.mul(z).div(y), origp.mul(p).div(y));
  }

  var count = 0;
  var ex, ey, ez, ep;

  var sqr = function(x) { return x.mul(x); };

  var sum = C(0);
  var fac = C(1);
  for(;;) {
    var a = x.add(y).add(z).add(p).add(p).divr(5);
    ex = x.div(a).rsub(1);
    ey = y.div(a).rsub(1);
    ez = z.div(a).rsub(1);
    ep = p.div(a).rsub(1);
    var error = Math.max(Math.max(ex.abs(), ey.abs()), Math.max(ez.abs(), ep.abs()));
    if(error < 0.0015) break; // converged
    var sx = C.sqrt(x);
    var sy = C.sqrt(y);
    var sz = C.sqrt(z);
    // The square roots must be calculated as sqrt(x)*sqrt(y), not as sqrt(x*y), to ensure correct branch (given that sqrt returns the solution with positive real part) (see mulSqrtBranch_)
    var lambda = sx.mul(sy).add(sy.mul(sz)).add(sz.mul(sx));
    var alpha = sqr(p.mul(sx.add(sy).add(sz)).add(sx.mul(sy).mul(sz)));
    var beta = p.mul(sqr(p.add(lambda)));
    sum = sum.add(fac.mul(C.ellipticrc(alpha, beta)));
    fac = fac.divr(4);
    x = x.add(lambda).divr(4);
    y = y.add(lambda).divr(4);
    z = z.add(lambda).divr(4);
    p = p.add(lambda).divr(4);
    count++;
    if(count > 1000) break; // avoid infinite loop
  }

  var xyz = ex.mul(ey).mul(ez);
  var ppp = ep.mul(ep).mul(ep);
  var e2 = ex.mul(ey).add(ex.mul(ez)).add(ey.mul(ez)).sub(ep.mul(ep).mulr(3));
  var e3 = xyz.add(e2.mul(ep).mulr(2)).add(ppp.mulr(4));
  var e4 = xyz.mulr(2).add(e2.mul(ep)).add(ppp.mulr(3)).mul(ep);
  var e5 = xyz.mul(ep).mul(ep);
  // extra error term (should be near 1)
  var extra = C(1).sub(e2.mulr(3.0/14)).add(e3.divr(6.0)).add(e2.mul(e2).mulr(9.0/88)).sub(e4.mulr(3.0/22)).sub(e2.mul(e3).mulr(9.0/52)).add(e5.mulr(3.0/26));
  var result = sum.mulr(3).add(fac.mul(extra).mul(a.powr(-1.5)));
  if(origp.re <= 0) {
    var rfx = C.ellipticrf(x, y, z);
    result = ca.mul(cb.mul(result).add(rcx.sub(rfx).mulr(3)));
  }
  return result;
};

//Elliptic integral: Carlson's form R_D
Jmat.Complex.ellipticrd = function(x, y, z) {
  // R_D(x, y) == R_J(x, y, z, z), but the implementation below is more efficient
  // Homogeneity: R_F(kx, ky, kz) = k^(-0.5) * R_F(x, y, z)
  // Duplication theorem: R_F(x, y, z) = R_F((x+l)/4, (y+l)/4, (z+l)/4) with l the lambda in the algorithm below.
  // The algorithm is based on the duplication theorem.
  var C = Jmat.Complex;
  var count = 0;
  var ex, ey, ez;
  var a;
  var sum = C(0);
  var fac = C(1);
  for(;;) {
    a = x.add(y).add(z.mulr(3)).divr(5);
    if(C.nearr(a, 0, 1e-14)) return C(0); // avoid NaNs, result is 0 in this case
    ex = x.div(a).rsub(1);
    ey = y.div(a).rsub(1);
    ez = z.div(a).rsub(1);
    var error = Math.max(Math.max(ex.abs(), ey.abs()), ez.abs());
    if(error < 0.0015) break; // converged
    var sx = C.sqrt(x);
    var sy = C.sqrt(y);
    var sz = C.sqrt(z);
    // The square roots must be calculated as sqrt(x)*sqrt(y), not as sqrt(x*y), to ensure correct branch (given that sqrt returns the solution with positive real part) (see mulSqrtBranch_)
    var lambda = sx.mul(sy).add(sy.mul(sz)).add(sz.mul(sx));
    sum = sum.add(fac.div(sz.mul(z.add(lambda))));
    fac = fac.divr(4);
    x = x.add(lambda).divr(4);
    y = y.add(lambda).divr(4);
    z = z.add(lambda).divr(4);
    count++;
    if(count > 1000) break; // avoid infinite loop
  }
  var e2 = ex.mul(ey).sub(ez.mul(ez).mulr(6));
  var e3 = ex.mul(ey).mulr(3).sub(ez.mul(ez).mulr(8)).mul(ez);
  var e4 = ex.mul(ey).sub(ez.mul(ez)).mulr(3).mul(ez).mul(ez);
  var e5 = ez.mul(ey).mul(ez).mul(ez).mul(ez);
  // extra error term (should be near 1)
  var extra = C(1).sub(e2.mulr(3.0/14)).add(e3.divr(6)).add(e2.mul(e2).mulr(9.0/88)).sub(e4.mulr(3.0/22)).sub(e2.mul(e3).mulr(9.0 / 52)).add(e5.mulr(3.0/26));
  return fac.mul(a.powr(-1.5)).mul(extra).add(sum.mulr(3));
};

// Elliptic integral K(k)
// BE CAREFUL! There exists another convention where the argument is m = k*k. However, here we use the sqrt of m.
Jmat.Complex.elliptick = function(k) {
  var C = Jmat.Complex;
  return C.ellipticrf(C(0), k.mul(k).rsub(1), C(1));
};

// Incomplete elliptic integral F(phi, k)
// BE CAREFUL! There exists another convention where the argument is m = k*k. However, here we use the sqrt of m.
Jmat.Complex.incellipticf = function(phi, k) {
  var C = Jmat.Complex;
  var s = C.sin(phi);
  var c = C.cos(phi);
  return s.mul(C.ellipticrf(c.mul(c), k.mul(k).mul(s).mul(s).rsub(1), C(1)));
};

// Elliptic integral E(k)
// BE CAREFUL! There exists another convention where the argument is m = k*k. However, here we use the sqrt of m.
Jmat.Complex.elliptice = function(k) {
  var C = Jmat.Complex;
  var m = k.mul(k);
  var f = C.ellipticrf(C(0), m.rsub(1), C(1));
  var d = C.ellipticrd(C(0), m.rsub(1), C(1));
  return f.sub(d.mul(m).divr(3));
};

// Incomplete elliptic integral E(phi, k)
// BE CAREFUL! There exists another convention where the argument is m = k*k. However, here we use the sqrt of m.
Jmat.Complex.incelliptice = function(phi, n, k) {
  var C = Jmat.Complex;
  var s = C.sin(phi);
  var c = C.cos(phi);
  var ss = s.mul(s);
  var sss = ss.mul(s);
  var cc = c.mul(c);
  var m = k.mul(k);
  var f = C.ellipticrf(cc, m.mul(ss).rsub(1), C(1));
  var d = C.ellipticrd(cc, m.mul(ss).rsub(1), C(1));
  return s.mul(f).add(n.mul(sss).mul(d).mulr(3));
};

// Elliptic integral PI(n, k)
// BE CAREFUL! There exists another convention where the argument is m = k*k. However, here we use the sqrt of m.
Jmat.Complex.ellipticpi = function(n, k) {
  var C = Jmat.Complex;
  var m = k.mul(k);
  var f = C.ellipticrf(C(0), m.rsub(1), C(1));
  var j = C.ellipticrj(C(0), m.rsub(1), C(1), n.rsub(1));
  return f.add(j.mul(n).divr(3));
};

// Incomplete elliptic integral PI(phi, n, k)
// BE CAREFUL! There exists another convention where the argument is m = k*k. However, here we use the sqrt of m.
Jmat.Complex.incellipticpi = function(phi, n, k) {
  var C = Jmat.Complex;
  var s = C.sin(phi);
  var c = C.cos(phi);
  var ss = s.mul(s);
  var sss = ss.mul(s);
  var cc = c.mul(c);
  var m = k.mul(k);
  var f = C.ellipticrf(cc, m.mul(ss).rsub(1), C(1));
  var j = C.ellipticrj(cc, m.mul(ss).rsub(1), C(1), n.mul(ss).rsub(1));
  return s.mul(f).add(n.mul(sss).mul(j).mulr(3));
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
    b2 = Jmat.Complex.mulSqrtBranch_(a, b);
    a = a2;
    b = b2;
  }

  if(real) a.im = 0; // it may be 1e-16 or so due to numerical error (only if both inputs are real and have same sign)

  return a;
};

//Geometric-Harmonic mean (iteratively calculated)
Jmat.Complex.ghm = function(a, b) {
  // Not really defined for negative and complex numbers, but when using Jmat.Complex.mulSqrtBranch_ it looks smooth in the 2D plot
  // NOTE: An alternative, that returns different values for neg/complex (but same for positive reals) is: return Jmat.Complex.agm(a.inv(), b.inv()).inv();
  var a2, b2;
  for(var i = 0; i < 60; i++) {
    if(a.eq(b)) break;
    a2 = Jmat.Complex.mulSqrtBranch_(a, b);
    b2 = Jmat.Complex.TWO.div(a.inv().add(b.inv()));
    a = a2;
    b = b2;
  }
  return a;
};


////////////////////////////////////////////////////////////////////////////////

// Loop for tetration: raises b to power a num times (typically a == b). if l is true, takes logarithm instead of power, to loop in the opposite direction.
Jmat.Real.tetration_loop_ = function(a, b, num, l) {
  var result = b;
  var last;
  for(var i = 0; i < num; i++) {
    if(l) result = Jmat.Real.logy(result, a);
    else result = Math.pow(a, result);
    if(Jmat.Real.isInfOrNaN(result)) return result;
    if(result == Infinity) return result; // Actually redundant, result == last already checks that too
    if(result == last) return result; // E.g. a=1.01 already converges to 1.0101015237405409 after just 7 iterations. 1.44 converges to 2.393811748202943 after 244 iterations, so numbers near that take a while to loop. 1.45 and higher converge to Infinity.
    last = result;
    if(i > 1000) return NaN; //avoid infinite loop
  }
  return result;
}

// Tetration
// Returns experimental (not mathematically correct) results unless x is an integer or Infinity
Jmat.Real.tetration = function(a, x) {
  var R = Jmat.Real;
  // if(a == 1) return 1; // Not true, e.g. 1^Infinity is not 1. Probably things are not like this for tetration either.
  if(x == 0) return 1; //by definition
  if(x == 1) return a;
  if(x == 2) return Math.pow(a, a);
  if(a >= 2 && x > 5) return Infinity; // too big for double
  if(a == 0 && R.isPositiveInt(x)) {
    // Chosen definition that considers 0^0 = 1, so 0^^n is 1 if n is even, 0 if n is odd.
    return Jmat.Complex.isEven(x) ? 1 : 0;
  }

  // Power tower (infinitely iterated exponentiation)
  if(x == Infinity && a > 0) {
    // converges if a >= 0.066 && a <= 1.44
    var l = Math.log(a);
    return R.lambertw(-l) / (-l);
  }

  // This supports power tower as well
  // It supports it better than the lambertw formula above because the domain here is real numbers only, the lambertw formula fails on, for example, -1 and returns nan instead of -1
  if(R.isPositiveInt(x)) {
    return R.tetration_loop_(a, a, x - 1, false);
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
    return R.tetration_loop_(a, b, Math.ceil(x), false);
  }
  if(x <= -1) {
    var b = x - Math.floor(x); //always in range 0-1
    return R.tetration_loop_(a, b, -Math.ceil(x), true);
  }

  return NaN;
};

// Loop for tetration: raises b to power a num times (typically a == b). if l is true, takes logarithm instead of power, to loop in the opposite direction.
Jmat.Complex.tetration_loop_ = function(a, b, num, l) {
  var C = Jmat.Complex;
  var result = b;
  var last;
  for(var i = 0; i < num; i++) {
    if(l) result = C.logy(result, a);
    else result = a.pow(result);
    if(C.isInfOrNaN(result)) return result;
    // E.g. a=1.01 already converges to 1.0101015237405409 after just 7 iterations. 1.44 converges to 2.393811748202943 after 244 iterations, so numbers near that take a while to loop. 1.45 and higher converge to Infinity.
    if(result.eq(last)) return result;
    last = result;
    if(i > 1000) return C(NaN); //avoid infinite loop (check is here, not before the loop with num, because the loop may converge)
  }
  return result;
};

// complex tetration e^^z, aka "tet" or "fse" (fast superexponential)
// This implementation, which only works for base e, is based on the paper: Tetrational as special function, E. Kouznetsov
Jmat.Complex.tetrational = function(z) {
  var C = Jmat.Complex;

  // TODO: more coefficients
  var s = [0.30685281944005, 0.59176735125832, 0.39648321290170, 0.17078658150959,
           0.08516537613999, 0.03804195209047, 0.01734090876306, 0.00755271038865,
           0.00328476064839, 0.00139361740170, 0.00058758348148, 0.00024379186661,
           0.00010043966462, 0.00004090111776, 0.00001654344436, 0.00000663102846,
           0.00000264145664, 0.00000104446533, 0.00000041068839, 0.00000016048059,
           0.00000006239367, 0.00000002412797, 0.00000000928797, 0.00000000355850,
           0.00000000135774, 0.00000000051587];
  // TODO: does JavaScript always re-execute all these constructors? Make object constant outside of this function instead
  var t = [C(0.37090658903229, 1.33682167078891), C(0.01830048268799, 0.06961107694975),
           C(-0.04222107960160, 0.02429633404907), C(-0.01585164381085, -0.01478953595879),
           C(0.00264738081895, -0.00657558130520), C(0.00182759574799, -0.00025319516391),
           C(0.00036562994770, 0.00028246515810), C(0.00002689538943, 0.00014180498091),
           C(-0.00003139436775, 0.00003583704949), C(-0.00001376358453, -0.00000183512708),
           C(-0.00000180290980, -0.00000314787679), C(0.00000026398870, -0.00000092613311),
           C(0.00000024961828, -0.00000013664223), C(0.00000000637479, 0.00000002270476),
           C(-0.00000000341142, 0.00000000512289), C(-0.00000000162203, 0.00000000031619),
           C(-0.00000000038743, -0.00000000027282), C(-0.00000000001201, -0.00000000013440),
           C(0.00000000002570, -0.00000000002543), C(0.00000000000935, 0.00000000000045),
           C(0.00000000000170, 0.00000000000186), C(-0.00000000000005, 0.00000000000071),
           C(-0.00000000000016, 0.00000000000012), C(-0.00000000000005, -0.00000000000001),
           C(-0.00000000000001, -0.00000000000001)];
  var an = [C(0.3181315052047641353, 1.3372357014306894089), C(1), C(-0.1513148971556517359, -0.2967488367322413067),
            C(-0.03697630940906762, 0.09873054431149697), C(0.0258115979731401398, -0.017386962126530755), C(-0.0079444196, 0.00057925018)];

  // fima = Fast approximation at large IMaginary part of the Argument
  var fima = function(z) {
    var r = C(1.0779614375280, -0.94654096394782);
    var beta = C(0.12233176, -0.02366108);
    var l = an[0]; //fixed point of the logarithm
    var e = C.exp(l.mul(z).add(r));
    var c = beta.mul(e).mul(C.exp(z.mul(C.newi(Math.PI * 2))));
    return C.powerSeries(an, an.length, C.ZERO, e).add(c);
  };

  // McLaurin expansion
  var maclo = function(z) {
    return C.log(z.addr(2)).add(C.powerSeries(s, s.length, C.ZERO, z));
  };

  // Taylor expansion
  var tai = function(z) {
    return C.powerSeries(t, t.length, C.newi(3), z);
  };

  var b;
  var z2 = C(Jmat.Real.fracn(z.re), z.im);

  if(z.im < -4.5) b = fima(z2.conj()).conj();
  else if(z.im < -1.5) b = tai(z2.conj()).conj();
  else if(z.im < 1.5) b = maclo(z2);
  else if(z.im < 4.5) b = tai(z2);
  else b = fima(z2);

  if(z.re > 0) return C.tetration_loop_(C.E, b, Math.floor(z.re), false);
  else return C.tetration_loop_(C.E, b, -Math.ceil(z.re), true);
};

// Tetration
// Returns experimental (not mathematically correct) results unless z is a positive integer or Infinity
Jmat.Complex.tetration = function(a, z) {
  var C = Jmat.Complex;

  if(C.isPositive(a) && C.isPositiveInt(z) && z.re != Infinity) {
    return C(Jmat.Real.tetration(a.re, z.re));
  }

  //if(a.eqr(1)) return C(1);  // Not true, e.g. 1^Infinity is not 1. Probably things are not like this for tetration either. Indeed it looks like e.g. 1^^0.5 is not 1.
  if(z.eqr(0)) return C(1); //by definition
  if(z.eqr(1)) return a;
  if(z.eqr(2)) return a.pow(a);
  if(C.isReal(a) && a.re >= 2 && z > 5) return C(Infinity); // too big for double
  if(a.eqr(0) && C.isPositiveInt(z)) {
    // Chosen definition that considers 0^0 = 1, so 0^^n is 1 if n is even, 0 if n is odd.
    return C.isEven(z) ? C(1) : C(0);
  }

  // Power tower (infinitely iterated exponentiation)
  if(z.eqr(Infinity) /*&& C.isPositive(a)*/) {
    if(a.eqr(1)) return C(1); //0/0 ==> 1
    // converges if a >= 0.066 && a <= 1.44
    // when using "C.tetration_loop_" with high iterations here, it it indeed only converges there. The lambertw formula below has values everywhere though.
    var l = C.log(a);
    return C.lambertw(l.neg()).div(l.neg());
  }

  // This supports power tower as well
  // It supports it better than the lambertw formula above because the domain here is real numbers only, the lambertw formula fails on, for example, -1 and returns nan instead of -1
  if(C.isPositiveInt(z)) {
    return C.tetration_loop_(a, a, z.re - 1, false);
  }

  // Everything above is true tetration for those cases where possible. What follows below is intermediate tetration research, to return "something", because there is no more than that available.
  if (C.isReal(z)) {
    // Linear approximation for the extension to real heights
    // a^^z = z+1 for z > -1 && z <= 0
    // a^^z = log_a(a^^(z+1)) for z <= -1  ==> a^^-1.5 = log_a(z+2), a^^-2.5 = log_a(log_a(z+3)), etc...
    // a^^z = a^(a^^(z-1)) for z > 0  ==> a^^0.5 = a^z, a^^1.5 = a^(a^(z-1)), a^^2.5 = a^(a^(a^(z-2))), etc...
    // examples: e^^(0.5*pi) ~= 5.868, 0.5^^(-4.3) ~= 4.033
    if(z.eqr(-1)) {
      return C(0); //the formulas below would give -Infinity
    }
    if(z.re > -1 && z.re <= 0) {
      return z.inc();
    }
    if(z.re > 0) {
      var b = z.sub(C.floor(z)); //always in range 0-1
      return C.tetration_loop_(a, b, Math.ceil(z.re), false);
    }
    if(z.re <= -1) {
      var b = z.sub(C.floor(z)); //always in range 0-1
      return C.tetration_loop_(a, b, -Math.ceil(z.re), true);
    }
  }

  if (C.near(a, C.E, 1e-15)) {
    return C.tetrational(z);
  }

  // TODO: for complex z and arbitrary base: implement something like "kneser" function

  return C(NaN);

  // TODO: implement super logarithm (slog), and super root (sroot) as well. Though, so far the only formulas I have is slog of base e in the Kouznetsov paper, and the 2-super root, so no implementation using both parameters can be done so far.
};

////////////////////////////////////////////////////////////////////////////////

// Exponential integral Ei
Jmat.Complex.ei = function(z) {
  var C = Jmat.Complex;
  var za = z.abs();

  if(za == 0) return C(-Infinity);
  if(za == Infinity) {
    if(z.re > 0) return C(Infinity);
    return C(0, Math.sign(z.im) * Math.PI);
  }

  if(za > 5) {
    // asymptotic series
    var result = C(0);
    var term = C(1);
    var stop = Math.min(30, Math.floor(z.abs()) + 1);
    for(var i = 1; i <= stop; i++) {
      var prev = term;
      term = term.mul(z.rdiv(i));
      if(C.nearr(term, 0, 1e-15)) {
        //return C.exp(z).div(z).mul(result.addr(1)).add(C(0, Math.sign(z.im) * Math.PI));
        var t = C(0, Math.sign(z.im) * Math.PI);
        result = result.addr(1).div(z);
        if(z.abs() > 100) {
          // this avoids overflow in the exp call
          return C.exp(z.add(C.log(result))).add(t);
        }
        return C.exp(z).mul(result).add(t);
      }
      if(term.abs() > prev.abs()) {
        // diverging
        result = result.sub(prev);
        //return C.exp(z).div(z).mul(result.addr(1)).add(C(0, Math.sign(z.im) * Math.PI));
        break;
      }
      result = result.add(term);
    }
    // no return

  }

  if(za > 1 && (z.re < 0 || Math.abs(z.im) > 1)) {
    // continued fraction
    var result = C(0, Math.sign(z.im) * Math.PI);
    var c = C(0), d = C(0);
    var e = C.exp(z);
    if(z.im == 0) {
      d = z.rsub(1).inv();
      result = d.mul(e).neg();
    } else {
      c = z.rsub(1).sub(e.div(result)).inv();
      d = z.rsub(1).inv();
      result = result.mul(d.div(c));
    }
    for(var i = 1; i < 30; i++) {
      c = z.rsub(i * 2 + 1).sub(c.mulr(i * i)).inv();
      d = z.rsub(i * 2 + 1).sub(d.mulr(i * i)).inv();
      var prev = result;
      result = result.mul(d.div(c));
      if(C.near(result, prev, 1e-15)) return result;
    }
    return result;
  }

  // power series
  var a = 0;
  if(z.im != 0) a = Math.sign(z.im) * Math.abs(C.arg(z).re);
  var result = C.EM.add(C.log(C.abs(z))).add(C(0, a));
  var zz = z;
  var ii = 1;
  for(var i = 1; i < 30; i++) {
    var next = result.add(zz.divr(i * ii));
    if(result.eq(next)) break; // converged
    result = next;
    zz = zz.mul(z);
    ii = ii * (i + 1);
  }
  return result;

};

// Exponential integral E1
Jmat.Complex.e1 = function(z) {
  var C = Jmat.Complex;
  var result = C.ei(z.neg()).neg();
  if(z.im < 0) result = result.add(C(0, Math.PI));
  else if(z.im > 0 || z.re < 0) result = result.sub(C(0, Math.PI));
  return result;
};

// Logarithmic integral: li
Jmat.Complex.li = function(z) {
  var C = Jmat.Complex;
  return C.ei(C.log(z));
};

// Offset logarithmic integral: Li(z) = li(z) - li(2)
Jmat.Complex.li2 = function(z) {
  var li_2 = 1.045163780117492784844588889194613136522615578151201575832909;
  return Jmat.Complex.li(z).subr(li_2);
};

// only a good approximation for |z| < 4
Jmat.Complex.si_pade_ = function(z) {
  var z2 = z.mul(z);
  var a = z2.mulr(-6.05338212010422477e-16).addr(7.08240282274875911e-13).
      mul(z2).addr(-3.53201978997168357e-10).mul(z2).addr(9.43280809438713025e-8).
      mul(z2).addr(-1.41018536821330254e-5).mul(z2).addr(1.15457225751016682e-3).
      mul(z2).addr(-4.54393409816329991e-2).mul(z2).addr(1);
  var b = z2.mulr(3.21107051193712168e-16).addr(4.5049097575386581e-13).
      mul(z2).addr(3.28067571055789734e-10).mul(z2).addr(1.55654986308745614e-7).
      mul(z2).addr(4.99175116169755106e-5).mul(z2).addr(1.01162145739225565e-2).
      mul(z2).addr(1);
  return z.mul(a).div(b);
};

// only a good approximation for |z| < 4
Jmat.Complex.ci_pade_ = function(z) {
  var C = Jmat.Complex;
  var z2 = z.mul(z);
  var a = z2.mulr(-9.93728488857585407e-15).addr(1.06480802891189243e-11).
      mul(z2).addr(-4.68889508144848019e-9).mul(z2).addr(1.05297363846239184e-6).
      mul(z2).addr(-1.27528342240267686e-4).mul(z2).addr(7.51851524438898291e-3).
      mul(z2).addr(-0.25).mul(z2);
  var b = z2.mulr(1.39759616731376855e-18).addr(1.89106054713059759e-15).
      mul(z2).addr(1.38536352772778619e-12).mul(z2).addr(6.97071295760958946e-10).
      mul(z2).addr(2.55533277086129636e-7).mul(z2).addr(6.72126800814254432e-5).
      mul(z2).addr(1.1592605689110735e-2).mul(z2).addr(1);
  return C.log(z).addr(0.577215664901532861).add(a.div(b));
};

// only a good approximation for large |z| > 4, and |z.re| > 1 as there are some "warts" on the imaginary axis
Jmat.Complex.f_pade_ = function(z) {
  var y = z.mul(z).inv();
  var a = y.mulr(-4.94701168645415959931e11).addr(4.94816688199951963482e12).
      mul(y).addr(1.00795182980368574617e13).mul(y).addr(4.20968180571076940208e12).
      mul(y).addr(6.40533830574022022911e11).mul(y).addr(4.33736238870432522765e10).
      mul(y).addr(1.43073403821274636888e9).mul(y).addr(2.37750310125431834034e7).
      mul(y).addr(1.96396372895146869801e5).mul(y).addr(7.44437068161936700618e2).
      mul(y).addr(1);
  var b = y.mulr(1.11535493509914254097e13).addr(1.43468549171581016479e13).
      mul(y).addr(5.06084464593475076774e12).mul(y).addr(7.08501308149515401563e11).
      mul(y).addr(4.58595115847765779830e10).mul(y).addr(1.47478952192985464958e9).
      mul(y).addr(2.41535670165126845144e7).mul(y).addr(1.97865247031583951450e5).
      mul(y).addr(7.46437068161927678031e2).mul(y).addr(1);
  return a.div(b).div(z);
};

// only a good approximation for large |z| > 4, and |z.re| > 1 as there are some "warts" on the imaginary axis
Jmat.Complex.g_pade_ = function(z) {
  var y = z.mul(z).inv();
  var a = y.mulr(-1.36517137670871689e12).addr(6.43291613143049485e12).
      mul(y).addr(1.81004487464664575e13).mul(y).addr(7.57664583257834349e12).
      mul(y).addr(1.09049528450362786e12).mul(y).addr(6.83052205423625007e10).
      mul(y).addr(2.06297595146763354e9).mul(y).addr(3.12557570795778731e7).
      mul(y).addr(2.35239181626478200e5).mul(y).addr(8.1359520115168615e2).mul(y).addr(1);
  var b = y.mulr(3.99653257887490811e13).addr(4.01839087307656620e13).
      mul(y).addr(1.17164723371736605e13).mul(y).addr(1.39866710696414565e12).
      mul(y).addr(7.87465017341829930e10).mul(y).addr(2.23355543278099360e9).
      mul(y).addr(3.26026661647090822e7).mul(y).addr(2.40036752835578777e5).
      mul(y).addr(8.19595201151451564e2).mul(y).addr(1);
  return y.mul(a).div(b);
};

// Sine integral: Si
Jmat.Complex.si = function(z) {
  var C = Jmat.Complex;

  if(z.abs() < 4) {
    return C.si_pade_(z);
  }

  if(Math.abs(z.re) > 1) {
    var f = C.f_pade_(z);
    var c = C.cos(z);
    var g = C.g_pade_(z);
    var s = C.sin(z);
    // l is the same as log(-ix) - log(ix), which is equal to pi or -pi depending on complex quadrant
    var l = C.I.mulr(z.re < 0 || (z.re == 0 && z.im < 0) ? Math.PI : -Math.PI);
    return l.muli(0.5).sub(f.mul(c)).sub(g.mul(s));
  }

  // In terms of e1
  var a = C.e1(z.mul(C.I).neg());
  var b = C.e1(z.mul(C.I));
  // l is the same as log(-ix) - log(ix), which is equal to pi or -pi depending on complex quadrant
  var l = C.I.mulr(z.re < 0 || (z.re == 0 && z.im < 0) ? Math.PI : -Math.PI);
  return a.sub(b).add(l).mul(C.I).mulr(0.5);
};

// Cosine integral: Ci
Jmat.Complex.ci = function(z) {
  var C = Jmat.Complex;

  if(z.abs() < 4) return C.ci_pade_(z);

  if(Math.abs(z.re) > 1) {
    var f = C.f_pade_(z);
    var c = C.cos(z);
    var g = C.g_pade_(z);
    var s = C.sin(z);
    var lz = C.log(z);
    var l = (z.re > 0 || (z.re == 0 && z.im > 0)) ? C(0) : C(0, (z.im < 0 ? -1 : 1) * Math.PI);
    return f.mul(s).sub(g.mul(c)).add(l);
  }

  // In terms of e1
  var a = C.e1(z.mul(C.I).neg());
  var b = C.e1(z.mul(C.I));
  var lz = C.log(z);
  // l is the same as log(-ix) + log(ix), which is equal to log(z) or log(-z) depending on complex quadrant
  var l = ((z.re > 0 || (z.re == 0 && z.im > 0)) ? lz : lz.subi((z.im < 0 ? -1 : 1) * Math.PI)).mulr(2);
  return a.add(b).add(l).mulr(-0.5).add(lz);
};

// Hyperbolic sine integral: Shi
Jmat.Complex.shi = function(z) {
  return Jmat.Complex.si(z.muli(1)).muli(1).neg();
};

// Hyperbolic cosine integral: Chi
Jmat.Complex.chi = function(z) {
  var C = Jmat.Complex;
  // l is same as log(z) - log(i*z)
  var l = new C(0, (z.re < 0 && z.im >= 0 ? (3 * Math.PI / 2) : (-Math.PI / 2)));
  return C.ci(z.muli(1)).add(l);
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Complex.erf_inv = function(z) {
  var C = Jmat.Complex;
  if (z.im != 0 && Math.abs(z.re) > 1) {
    //this branch is taken for large complex numbers because the implementation below doesn't work well on those. This one isn't much better btw, but slightly less bad for those cases.
    //TODO: is incorrect on many complex numbers!! e.g. erf_inv(erf(5 + 5i)) gives way wrong result
    //var zzpi = z.mul(z).mulr(Math.PI);
    //return zzpi.mulr(34807/182476800.0).addr(4369/5806080.0).mul(zzpi).addr(127/40320.0).mul(zzpi).addr(7/480.0).mul(zzpi).addr(1/12.0).mul(zzpi).addr(1).mul(z).mul(C.SQRTPI).mulr(0.5);

    // With newton method

    // derivative of erf is: 2/sqrt(pi) * e^(-x^2)
    var derf = function(x) {
      return C.TWO.divr(Jmat.Real.SQRTPI).mul(C.exp(x.mul(x).neg()));
    };

    var neg = z.re < 0;
    if(neg) z = z.neg();

    // For abs(z) > 1 and z.re > 0, the following starting value works well: sqrt(-log(x * sqrt(pi) * (1 - x)))
    // NOTE: for abs(z) < 1, instead z*sqrtpi/2 would work, but other code already handles such case
    var start = C.sqrt(C.log(z.mulr(Jmat.Real.SQRTPI).mul(C.ONE.sub(z))).neg());
    // NOTE: erf_inv has multiple solutions, with this starting value only one particular one is returned.
    // e.g. with the chosen starting value, erf_inv(2+2i) gives 0.386507600275 + 1.320769860731i. But with starting value 0, it gives the also correct 2.947736167125 + 3.401249486995i.
    var result = C.finvert_newton(z, C.erf, derf, start);
    if(neg) result = result.neg();
    return result;
  } else {
    //if (a > 1) return C(NaN); //only relevant for real numbers
    if (z.im == 0) {
      if (z.re == 0) return C(0);
      if (z.re == 1) return C(Infinity);
      if (z.re == -1) return C(-Infinity);
    }

    var erf_inv_a_ = [0.886226899, -1.645349621, 0.914624893, -0.140543331];
    var erf_inv_b_ = [1, -2.118377725, 1.442710462, -0.329097515, 0.012229801];
    var erf_inv_c_ = [-1.970840454, -1.62490649, 3.429567803, 1.641345311];
    var erf_inv_d_ = [1, 3.543889200, 1.637067800];

    var a = z.abs();
    if (a <= 0.7) {
      var z2 = z.mul(z);
      var r = z.mul(z2.mulr(erf_inv_a_[3]).addr(erf_inv_a_[2]).mul(z2).addr(erf_inv_a_[1]).mul(z2).addr(erf_inv_a_[0]));
      r = r.div(z2.mulr(erf_inv_b_[4]).addr(erf_inv_b_[3]).mul(z2).addr(erf_inv_b_[2]).mul(z2).addr(erf_inv_b_[1]).mul(z2).addr(erf_inv_b_[0]));
    }
    else {
      var y = C.sqrt(C.log(C.ONE.sub(z).divr(2)).neg());
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

// Fresnel integral S(z)
Jmat.Complex.fresnels = function(z) {
  var C = Jmat.Complex;
  var R = Jmat.Real;
  var a = C.erf(C(1, 1).mulr(R.SQRTPI / 2).mul(z));
  var b = C.erf(C(1, -1).mulr(R.SQRTPI / 2).mul(z));
  return a.sub(b.muli(1)).mul(C(0.25, 0.25));
};

// Fresnel integral C(z)
Jmat.Complex.fresnelc = function(z) {
  var C = Jmat.Complex;
  var R = Jmat.Real;
  var a = C.erf(C(1, 1).mulr(R.SQRTPI / 2).mul(z));
  var b = C.erf(C(1, -1).mulr(R.SQRTPI / 2).mul(z));
  return a.add(b.muli(1)).mul(C(0.25, -0.25));
};

////////////////////////////////////////////////////////////////////////////////

//Minkowski's question mark function, from Wikipedia
Jmat.Real.minkowski = function(x) {
  if(x != x) return NaN;
  var p = Math.floor(x);
  var q = 1, r = p + 1, s = 1, m, n;
  var d = 1.0, y = p;
  if(x < p || (p < 0) != (r <= 0)) return x; //out of range ?(x) =~ x
  for(;;) { //invariants: q*r-p*s==1 && p/q <= x && x < r/s
    d /= 2; if(y + d == y) break; //reached max possible precision
    m = p + r; if((m < 0) != (p < 0)) break; //sum overflowed
    n = q + s; if(n < 0) break; //sum overflowed

    if(x < m / n) r = m, s = n;
    else y += d, p = m, q = n;
  }
  return y + d; //final round-off
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
//this is totally not an efficient algorithm like Newton's method let alone Brent's method.
//valuex and valuey are allowed to be anything (real) and to have function value with the same sign.
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

//TODO: does not work well
Jmat.Complex.rootfind_secant = function(z0, z1, f, maxit, prec) {
  maxit = maxit || 30;
  prec = prec || 1e-15;

  var f0 = f(z0);
  var f1 = f(z1);

  for(var i = 0; i < maxit; i++) {
    var zn = z0.mul(f1).sub(z1.mul(f0)).div(f1.sub(f0));
    if(Jmat.Complex.isInfOrNaN(zn)) {
      //TODO: fix z0 and z1 in some way to avoid this problem...
    }
    var fn = f(zn);
    if(Jmat.Real.near(fn, 0, prec)) return zn;
    z0 = z1;
    f0 = f1;
    z1 = zn;
    f1 = fn;
  }

  return z1; //not found, give best guess
};

//root finding, aka find zeroes (findZero)
//most parameters are optional. Based on which are given, a certain algorithm is chosen (newton, bisection, ...)
//f: the function to find root of
//o: object with the following optional values:
// o.z0: starting value, or, if z1 is given, lower value of range (and starting value is assumed in the center of both). Default: Jmat.Complex.ZERO
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
  return Jmat.Complex.rootfind_newton_noderiv(f, z0, maxit);
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

Jmat.Complex.newtonStartValuesAround_ = [ Jmat.Complex(1), Jmat.Complex(-1), Jmat.Complex.newi(1), Jmat.Complex.newi(-1) ];

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

//Like Jmat.Real.rootfind_newton, but fdf returns both f(x) and f'(x) in an array
Jmat.Real.rootfind_newton2 = function(fdf, z0, maxiter) {
  if (!z0) z0 = 0;
  if (!maxiter) maxiter = 30;
  var z = z0;
  var prevz = z;
  var bestdist = Infinity;
  var best = NaN;
  for (var i = 0; i < maxiter; i++) {
    var prevz = z;
    var v = fdf(z);
    z -= v[0] / v[1];
    var dist = Infinity;
    if(Jmat.Real.isInfOrNaN(z)) z = prevz - 0.1; // get out of singularity. TODO: improve this to choose correct direction
    else dist = Math.abs(z - prevz);
    if(dist < bestdist) {
      i = 0;
      best = z;
      bestdist = dist;
      if(dist < 1e-15) break; // Near enough, stop iterations
    }
  }

  return best;
};

//Newton-Raphson. Finds a root (zero) given function f, its derivative df, and an initial value z0
Jmat.Real.rootfind_newton = function(f, df, z0, maxiter) {
  return Jmat.Real.rootfind_newton2(function(t) { return [f(t), df(t)]}, z0, maxiter);
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
Jmat.Complex.rootfind_newton_noderiv = function(f, z0, maxiter) {
  return Jmat.Complex.rootfind_newton(f, function(x) {
    return Jmat.Complex.differentiate_stencil5(x, f);
  }, z0, maxiter);
};

//find result of inverse function using the newton method
Jmat.Complex.finvert_newton = function(z, f, df, z0, maxiter) {
  return Jmat.Complex.rootfind_newton(function(x) { return f(x).sub(z); }, df, z0, maxiter);
};

//find result of inverse function using the newton method (no need to give the derivative)
Jmat.Complex.finvert_newton_noderiv = function(z, f, z0,  maxiter) {
  return Jmat.Complex.rootfind_newton_noderiv(function(x) { return f(x).sub(z); }, z0, maxiter);
};

/*
Orthogonal polynomials for Gaussian quadrature:

Uses the n roots (with n the order) of a polynomial in the relevant interval,
with a weight for each root, and evaluates the function at each root.

The polynomial is evaluated with a recursive definition in each case, and the
resulting roots and weights are cached for reuse of later integrations with the
same order.

Legendre: -1..1
P0(x) = 1
P1(x) = x
Pn(x) = ((2n-1)*x*P(n-1)(x) - (n-1)*P(n-2)(x)) / n
Pn'(x) = (n*P(n-1)(x) - n*x*Pn(x)) / (1-x*x)
xi = roots of the polynomial in the interval
wi = 2 / ((1 - xi^2) * Pn'(xi)^2)

Laguerre: 0..oo
L0(x) = 1
L1(x) = 1-x
Ln(x) = ((2n-1-x)*L(n-1)(x) - (n-1)*L(n-2)(x))/n
Ln'(x) = (n*Ln(x) - n*L(n-1)(x)) / x
xi = roots of the polynomial in the interval
wi = -1 / (order * Ln'(xi) * L(n-1)(xi))

TODO: TanhSinh quadrature
*/

Jmat.Real.gausleg_cache_ = [];

// returns roots and weights for gauss-legendre quadrature
// caches the calculation result for the given order
Jmat.Real.gauss_legendre_roots_ = function(order) {
  if (Jmat.Real.gausleg_cache_[order]) return Jmat.Real.gausleg_cache_[order];

  // returns Pn(x) and Pn'(x), with n the order and P legendre polynomial
  var legendre_eval = function(x, order) {
     var p0 = 1;
     var p1 = x;
     for(var i = 2; i <= order; i++) {
       var pi = ((2 * i - 1) * x * p1 - (i - 1) * p0) / i;
       p0 = p1;
       p1 = pi;
     }
     var dp1 = (order * p0 - order * x * p1) / (1 - x * x);
     return [p1, dp1];
  };

  var roots = [];
  var weights = [];
  for (var i = 0; i < order; i++) {
    // initial root guess heuristics
    var x = Math.cos(Math.PI * (i + 0.75) / (order + 0.5));
    x = Jmat.Real.rootfind_newton2(function(t) { return legendre_eval(t, order); }, x, 100);
    var p = legendre_eval(x, order);
    roots[i] = x;
    weights[i] = 2 / ((1 - x * x) * p[1] * p[1]);
  }

  Jmat.Real.gausleg_cache_[order] = [roots, weights];
  return Jmat.Real.gausleg_cache_[order];
};

Jmat.Real.gauslag_cache_ = [];

// returns roots and weights for gauss-laguerre quadrature
// caches the calculation result for the given order
Jmat.Real.gauss_laguerre_roots_ = function(order) {
  if (Jmat.Real.gauslag_cache_[order]) return Jmat.Real.gauslag_cache_[order];

  // returns value of Ln(x), Ln'(x) and L(n-1)(x), with n the order and L laguerre polynomial
  var laguerre_eval = function(x, order) {
    if(order == 0) return [1, 0, 0];
    var l0 = 1;
    var l1 = 1 - x;
    for(var i = 2; i <= order; i++) {
      var li = ((2 * i - 1 - x) * l1 - (i - 1) * l0) / i;
      l0 = l1;
      l1 = li;
    }
    var dl1 = order * (l1 - l0) / x;
    return [l1, dl1, l0];
  };

  var roots = [];
  var weights = [];
  for(var i = 0; i < order; i++) {
    // initial root guess heuristics
    var x = 0;
    if (i == 0) {
      x = 3 / (1 + 2.4 * order);
    } else if(i == 1) {
      x = roots[i - 1] + 15 / (1 + 2.5 * order);
    } else {
      x = roots[i - 1] + (1 + 2.55 * (i - 1)) / (1.9 * (i - 1)) * (roots[i - 1] - roots[i - 2]);
    }
    x = Jmat.Real.rootfind_newton2(function(t) { return laguerre_eval(t, order); }, x, 100);
    var p = laguerre_eval(x, order);
    roots[i] = x;
    weights[i] = -1 / (order * p[1] * p[2]);
  }

  Jmat.Real.gauslag_cache_[order] = [roots, weights];
  return Jmat.Real.gauslag_cache_[order];
};

//numerical integration, aka quadrature
//integrate function f using simpsons rule, from x to y, with amount of steps given by steps parameter (higher = more precision, number of evaluations is thoroughly steps * 2)
//NOTE: this is the real version, real JS numbers only. Complex version is Jmat.Complex.integrate_simpson
Jmat.Real.integrate_simpson = function(x, y, f, steps, stopLoop) {
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

Jmat.Real.integrate_gaussian_legendre = function(x, y, f, order) {
  var R = Jmat.Real;
  if(order < 1 || order != Math.round(order)) return NaN;
  if(order > 100) return NaN; // too slow

  var l = R.gauss_legendre_roots_(order); // roots, weights
  var roots = l[0];
  var weights = l[1];

  var x2 = (y - x) / 2;
  var y2 = (y + x) / 2;
  var sum = 0;
  for (var i = 0; i < order; i++) {
    sum += f(x2 * roots[i] + y2) * weights[i];
  }
  return x2 * sum;
};


// Integrates from a to infinity
Jmat.Real.integrate_gaussian_laguerre = function(a, f, order) {
  var R = Jmat.Real;
  if(order < 1 || order != Math.round(order)) return NaN;
  if(order > 100) return NaN; // too slow

  var l = R.gauss_laguerre_roots_(order); // roots, weights
  var roots = l[0];
  var weights = l[1];

  var sum = 0;
  for (var i = 0; i < roots.length; i++) {
    var w = Math.exp(-roots[i]); // weighing function for laguerre-gaussian
    var v = f(roots[i] + a);
    sum += v / w * weights[i];
  }
  return sum;
};

// Gaussian quadrature (legendre for -1..1, laguerre for 0..oo)
// E.g. try Jmat.Real.integrate_gaussian(0, Infinity, function(t) { return 1 / Math.exp(t); }, 5);
Jmat.Real.integrate_gaussian = function(x, y, f, order) {
  var R = Jmat.Real;

  if(y == Infinity) {
    if(x == -Infinity) {
      // TODO: hermite
      var a = R.integrate_gaussian_laguerre(0, function(t) { return f(-t); }, order);
      var b = R.integrate_gaussian_laguerre(0, f, order);
      return a + b;
    }
    return R.integrate_gaussian_laguerre(x, f, order);
  }
  if(x == -Infinity) {
    return R.integrate_gaussian_laguerre(-x, function(t) { return f(-t); }, order);
  }

  return R.integrate_gaussian_legendre(x, y, f, order);
};

//numerical integration, aka quadrature
//opt_type: 0: simpson, 1: gaussian. Default: 1
//NOTE: this is the real version, real JS numbers only. Complex version is Jmat.Complex.integrate
//NOTE: using this function requires experimentation for your usecase, it depends a lot on the function what parameters work best. With gaussian, steps > 30 usually makes it worse, not better.
Jmat.Real.integrate = function(x, y, f, steps, opt_type) {
  if(opt_type == undefined) opt_type = 1;
  if(!steps) steps = opt_type == 1 ? 20 : 30;

  if(opt_type == 0) return Jmat.Real.integrate_simpson(x, y, f, steps);
  return Jmat.Real.integrate_gaussian(x, y, f, steps);
};

Jmat.Complex.integrate_gaussian_legendre = function(x, y, f, order) {
  var C = Jmat.Complex;
  if(order < 1 || order != Math.round(order)) return C(NaN);
  if(order > 100) return C(NaN); // too slow

  var l = Jmat.Real.gauss_legendre_roots_(order); // roots, weights
  var roots = l[0];
  var weights = l[1];

  var x2 = (y.sub(x)).divr(2);
  var y2 = (y.add(x)).divr(2);
  var sum = C(0);
  for (var i = 0; i < order; i++) {
    sum = sum.add(f(x2.mulr(roots[i]).add(y2)).mulr(weights[i]));
  }
  return x2.mul(sum);
};


// Integrates from a to infinity
Jmat.Complex.integrate_gaussian_laguerre = function(a, f, order) {
  var C = Jmat.Complex;
  if(order < 1 || order != Math.round(order)) return C(NaN);
  if(order > 100) return C(NaN); // too slow

  var l = Jmat.Real.gauss_laguerre_roots_(order); // roots, weights
  var roots = l[0];
  var weights = l[1];

  var sum = C(0);
  for (var i = 0; i < roots.length; i++) {
    var w = Math.exp(-roots[i]); // weighing function for laguerre-gaussian
    var v = f(C(roots[i]).add(a));
    sum = sum.add(v.divr(w).mulr(weights[i]));
  }
  return sum;
};

// Gaussian quadrature (legendre for -1..1, laguerre for 0..oo)
// NOTE: order over 20 or 30 is not useful. Sometimes less order is more
// precise, e.g. if f is quadratic or cubic, then order=2 is most precise
// Tricks for more precision: the smaller the derivatives of the function, the more precise, so try to transform to other function with lower derivatives
Jmat.Complex.integrate_gaussian = function(x, y, f, order) {
  var C = Jmat.Complex;

  if(y.eqr(Infinity)) {
    if(x.eqr(-Infinity)) {
      // TODO: hermite
      var a = C.integrate_gaussian_laguerre(C(0), function(t) { return f(t.neg()); }, order);
      var b = C.integrate_gaussian_laguerre(C(0), f, order);
      return a.add(b);
    }
    return C.integrate_gaussian_laguerre(x, f, order);
  }
  if(x.eqr(-Infinity)) {
    return C.integrate_gaussian_laguerre(x.neg(), function(t) { return f(t.neg()); }, order);
  }

  return C.integrate_gaussian_legendre(x, y, f, order);
};

//numerical integration, aka quadrature
//integrate function f using simpsons rule, from x to y, with amount of steps given by steps parameter (higher = more precision, number of evaluations is thoroughly steps * 2)
Jmat.Complex.integrate_simpson = function(x, y, f, steps, stopLoop) {
  var step = y.sub(x).divr(steps);
  var result = Jmat.Complex(0);
  var fa = null;

  var a, b;
  var infmul = 2;
  var inflow = 1e-12;
  var infhigh = 1e12;
  if(y.eqr(Infinity)) infmul = Math.pow(infhigh / inflow, 1.0 / steps);

  for(var i = 0; i <= steps; i++) {
    if(y.eqr(Infinity)) {
      a = (b ? b : x);
      b = (b ? b.mulr(infmul) : x.addr(inflow));
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
//opt_type: 0: simpson, 1: gaussian. Default: 1
//NOTE: using this function requires experimentation for your usecase, it depends a lot on the function what parameters work best. With gaussian, steps > 30 usually makes it worse, not better.
Jmat.Complex.integrate = function(x, y, f, steps, opt_type) {
  if(opt_type == undefined) opt_type = 1;
  if(!steps) steps = opt_type == 1 ? 20 : 30;

  if(opt_type == 0) return Jmat.Complex.integrate_simpson(x, y, f, steps);
  return Jmat.Complex.integrate_gaussian(x, y, f, steps);
};

// differentiation with just two points (finite difference, or secant)
Jmat.Complex.differentiate_newton_noderiv = function(x, f) {
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
  var i = 0;
  for(var z = x; z <= y + step / 4; z += step) {
    var fz = f(Jmat.Complex(z));
    result = result.add(fz);

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return Jmat.Complex(NaN);
    i++;
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
  var i = 0;
  // the step / 4 thing is to avoid numerical problems that may let it miss the last value
  for(var z = x; z <= y + step / 4; z += step) {
    var fz = f(Jmat.Complex(z));
    result = result.mul(fz);

    if(!!stopLoop && i % 50 == 49 && stopLoop()) return Jmat.Complex(NaN);
    i++;
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

pdf = probability density function for continuous distributions
pmf = probability mass function for discrete distributions
cdf = cumulative distribution function, integral of the pdf
qf = quantile function, the inverse function of cdf
*/


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
  var a = x.sub(mu).divr(sigma.abs()).div(Jmat.Complex.SQRT2); // (x-mu) / sqrt(2*sigma^2)
  return Jmat.Complex.erf(a).addr(1).mulr(0.5);
};

Jmat.Complex.qf_normal = function(x, mu, sigma) {
  return mu.add(sigma.mul(Jmat.Complex.qf_standardnormal(x)));
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
  var g = Jmat.Complex.calcCache_(nu, Jmat.Complex.pdf_studentt_cachefun_, Jmat.Complex.pdf_studentt_cache_);
  var b = Jmat.Complex.incbeta(x.mul(x).div(nu).neg(), Jmat.Complex(0.5), Jmat.Complex.ONE.sub(nu).mulr(0.5));
  var n = Jmat.Complex.I.mul(x).mul(b);
  var d = Jmat.Complex(x.abs()).mulr(2).mul(Jmat.Complex.SQRTPI);
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

// Fisher-Snedecor F distribution
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
  // 1/2b * exp(-|x-mu|/b)
  var e = Jmat.Complex.exp(Jmat.Complex(x.sub(mu).abs()).neg().div(b));
  return b.mulr(2).inv().mul(e);
};

Jmat.Complex.cdf_laplace = function(x, mu, b) {
  // 1/2b * exp(-|x-mu|/b)
  var e = Jmat.Complex.exp(Jmat.Complex(x.sub(mu).abs()).neg().div(b));
  var s = Jmat.Complex.sign(x.sub(mu));
  return e.rsub(1).mul(s).mulr(0.5).addr(0.5);
};

Jmat.Complex.qf_laplace = function(x, mu, b) {
  var l = Jmat.Complex.log(Jmat.Complex(x.subr(0.5).abs()).mulr(2).rsub(1));
  var s = Jmat.Complex.sign(x.subr(0.5));
  return mu.sub(b.mul(s).mul(l));
};

////////////////////////////////////////////////////////////////////////////////

//p in range 0-1, k integer 0 or 1
Jmat.Complex.pmf_bernoulli = function(k, p) {
  if(k.eqr(0)) return p.rsub(1);
  if(k.eqr(1)) return p;
  return p.pow(k).mul(p.rsub(1).pow(k.rsub(1))); //outside of definition but let's return something continuous...
};

Jmat.Complex.cdf_bernoulli = function(k, p) {
  if(k.re < 0) return Jmat.Complex.ZERO;
  if(k.re >= 1) return Jmat.Complex.ONE;
  return p.rsub(1); //1 - p
};

Jmat.Complex.qf_bernoulli = function(k, p) {
  if(k.re < 1-p.re) return Jmat.Complex.ZERO;
  else return Jmat.Complex.ONE;
};

////////////////////////////////////////////////////////////////////////////////

// returns probability of getting exactly k successes out of n trials with probability p each (p in range 0.0-1.0, k and n integers)
Jmat.Complex.pmf_binomial = function(k, n, p) {
  var b = Jmat.Complex.binomial(n, k);
  return b.mul(p.pow(k)).mul(p.rsub(1).pow(n.sub(k)));
};

Jmat.Complex.cdf_binomial = function(k, n, p) {
  return Jmat.Complex.beta_i(p, n.sub(k), k.addr(1));
};

// inverse of cdf_binomial in k. Works only for integer result in range 0-n
Jmat.Complex.qf_binomial = function(k, n, p) {
  if(p.eqr(0)) return Jmat.Complex.ZERO;
  if(p.eqr(1)) return n;
  var C = Jmat.Complex;
  var result = Jmat.Complex.rootfind_bisection(C(0), n, function(z) {
    return Jmat.Complex.cdf_binomial(z, n, p).sub(k);
  }, 100, 1e-14);
  return Jmat.Complex.round(result);  // Binomial pmf/cdf takes integer k, so qf, which is inverse of cdf in k, is expected to return integer.
};

////////////////////////////////////////////////////////////////////////////////

//Poisson distribution, with k integer (but still given as Jmat.Complex object)

Jmat.Complex.pmf_poisson = function(k, lambda) {
  var a = lambda.pow(k);
  var b = Jmat.Complex.factorial(k);
  var c = Jmat.Complex.exp(lambda.neg());
  return a.mul(c).div(b);
};

Jmat.Complex.cdf_poisson = function(k, lambda) {
  return Jmat.Complex.gamma_q(k.addr(1), lambda);
};

// TODO: make way more accurate
// Does not work well.
Jmat.Complex.qf_poisson = function(k, lambda) {
  return Jmat.Complex.gamma_q_inva(lambda, k).subr(1);
};
