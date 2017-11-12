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

// REQUIRES: jmat_real.js, jmat_complex.js, jmat_matrix.js, jmat_quaternion.js, jmat_special.js, jmat_bigint.js, jmat.js

// Unit tests for Jmat.js

// constructor
Jmat.Test = function() {
  // empty, it's a namespace
};

Jmat.Test.expectTrue = function(value, opt_message) {
  var message = opt_message ? ('fail: ' + opt_message) : 'fail';
  if (!value) throw message;
}

Jmat.Test.expectFalse = function(value, opt_message) {
  var message = opt_message ? ('fail: ' + opt_message) : 'fail';
  if (value) throw message;
}

// Works both for Complex or Matrix objects.
// Precision is number of decimal digits that should match
Jmat.Test.expectNear = function(e, a, precision) {
  var errorText = 'fail: expected ' + Jmat.toString(e) + ' got ' + Jmat.toString(a) + '. Expected precision: ' + precision;

  if(Jmat.eq(e, 0) || Jmat.eq(a, 0) || Jmat.matrixIn_(e) || Jmat.matrixIn_(a)) {
    // for 0, allow the other to be absolute rather than relative near it. For matrices, always use absolute epsilon as well.
    if(!Jmat.near(e, a, precision)) throw errorText;
  } else if (Jmat.bigIntIn_(e) || Jmat.bigIntIn_(a)) {
    if(!Jmat.eq(e, a)) throw errorText;
  } else if(typeof e == 'number') {
    if(Jmat.isNaN(e) && Jmat.isNaN(a)) return; //both NaN is ok for test
    if(!Jmat.relnear(e, a, precision)) throw errorText;
  } else {
    if(!Jmat.relnear(e, a, precision)) throw errorText;
  }
};

Jmat.Test.accrep = ''; // accuracy report string
Jmat.Test.accworst = 0; // worst accuracy seen

// Expect that the result of the mathematical function f with the arguments of var_arg, is near the expected result. Some numerical intolerance is allowed.
// The name must be given to report about accuracy.
// Set epsilon to NaN to not test anything, just report.
Jmat.Test.testFunction = function(expected, epsilon, f, var_arg) {
  var a = Array.prototype.slice.call(arguments).slice(3); /*var_arg*/
  var result = f.apply(this, a);
  if(!Jmat.bigIntIn_(result)) {
    var e = Jmat.cheb(expected, result);
    if(epsilon != Infinity && e > 0) {
      if(e > Jmat.Test.accworst) Jmat.Test.accworst = e;
      Jmat.Test.accrep += (f.testName ? f.testName : 'unk') + '(' + var_arg + ') = ' + expected + ', got: ' + result + ', eps: ' + e + '\n';
    }
  }
  if(!isNaN(epsilon)) Jmat.Test.expectNear(expected, result, epsilon);
};

//u,s,v = expected values
//m = input matrix
Jmat.Test.testSVD = function(u, s, v, epsilon, m) {
  // TODO: This test should tolerate some differences in signs of vectors in
  // u and v, because multiple solutions are possible and different software
  // returns different variants.
  var svd = Matrix.svd(Matrix.cast(m));
  Jmat.Test.expectNear(u, svd.u, epsilon);
  Jmat.Test.expectNear(s, svd.s, epsilon);
  Jmat.Test.expectNear(v, svd.v, epsilon);
};

//l,v = expected values. Note that eigenvectors in v are columns, not rows.
//m = input matrix
Jmat.Test.testEIG = function(l, v, epsilon, m) {
  var eig = Matrix.eig(Matrix.cast(m));
  Jmat.Test.expectNear(l, eig.l, epsilon);
  Jmat.Test.expectNear(v, eig.v, epsilon);
};

Jmat.Test.annotateFunctionNames = function() {
  var objects = [Jmat, Jmat.Real, Jmat.Complex, Jmat.Quaternion, Jmat.BigInt, Jmat.BigIntC, Jmat.Matrix];
  for(var i = 0; i < objects.length; i++) {
    for(var p in objects[i]) {
      if(!objects[i].hasOwnProperty(p)) continue;
      var o = objects[i][p];
      if(o) o.testName = p;
    }
  }
};

Jmat.Test.initAccuracy = function() {
  Jmat.Test.annotateFunctionNames();
  Jmat.Test.accrep = '';
  Jmat.Test.accworst = 0;
};

Jmat.Test.finalizeAccuracy = function(opt_verbose) {
  if(opt_verbose) console.log(Jmat.Test.accrep);
  console.log('precision: ' + Jmat.Test.accworst);
};

Jmat.Test.testPolylog = function() {
  Jmat.Test.initAccuracy();
  var polylog = function(a, b) { return Jmat.Complex.polylog(Jmat.Complex.cast(a), Jmat.Complex.cast(b)); };
  Jmat.Test.testFunction(Infinity, NaN, polylog, 1, 1);
  Jmat.Test.testFunction('-0.69314718055994530941', NaN, polylog, 1, -1);
  Jmat.Test.testFunction('-3.14159265358979323846i', NaN, polylog, 1, 2);
  Jmat.Test.testFunction('0.80612672304285226132', NaN, polylog, 0.5, 0.5);
  Jmat.Test.testFunction('-0.37375223798097305644', NaN, polylog, 0.5, -0.5);
  Jmat.Test.testFunction('1.34725375273575069219', NaN, polylog, -0.5, 0.5);
  Jmat.Test.testFunction('-0.28301281074650602317', NaN, polylog, -0.5, -0.5);
  Jmat.Test.testFunction('10.12002396823740075704-0.01575172198981096218i', NaN, polylog, 10, 10);
  Jmat.Test.testFunction('-0.09236458657786628703+9.98692954551758154675i', NaN, polylog, 10, '10i');
  Jmat.Test.testFunction('-43.94797353817581944710+13.70036401689160455461i', NaN, polylog, -10, '10i');
  //Jmat.Test.testFunction('-771.09626291551731193169-29255.534755954886987832i', NaN, polylog, '0.5+10i', '0.5+10i'); // result very large numbers, eps 0.15 but that clouds the other tests too much
  Jmat.Test.testFunction('-4.17886432849095879367+1.50745493019209405537i', NaN, polylog, '0.1+2i', '0.1+2i');
  Jmat.Test.finalizeAccuracy(true);
};

Jmat.doEigenPrecisionBenchmark = function() {
  var restore_seed = Jmat.Real.seed;
  Jmat.Real.seed = 5;

  var biggestError = 0;
  var biggestErrorData = null;

  var sumError = 0;
  var numError = 0;

  for (var i = 0; i < 100; i++) {
    var m = M.random({'square':true, 'real':true});
    var e = M.eig(m);

    //console.log(Jmat.render(m));
    //console.log(Jmat.render(e));

    for (var j = 0; j < m.w; j++) {

      var v = M.col(e.v, j);
      var v2 = m.mul(v); // multiplied with matrix
      var v3 = v.mulc(e.l[j]); // multiplied with eigenvalue. we expect this to be the same vector

      var error = Jmat.dist(v2, v3);
      //console.log('error: ' + error);

      //console.log(Jmat.render(v2));
      //console.log(Jmat.render(v3));

      if(error > biggestError) {
        biggestError = error;
        biggestErrorData = {col: j, m: m, l: e.l, v: e.v};
      }
      sumError += error;
      numError++;;
    }
  }

  console.log('average error: ' + (sumError / numError) +
              ', biggest error: ' + biggestError +
              ', data: ' + Jmat.render(biggestErrorData));

  Jmat.Real.seed = restore_seed;
};

Jmat.Test.doUnitTestLogGamma = function() {
  var R = Jmat.Real;
  var C = Jmat.Complex;
  var eps = 1e-10;
  Jmat.Test.testFunction(C('1.8989912736759002200831-0.8274647077730757440276i'), eps, C.loggamma, C('0.1+0.1i'));
  Jmat.Test.testFunction(C('-7857.23612011+37585.1962512i'), eps, C.loggamma, C('0.01+5000i'));
  Jmat.Test.testFunction(C('0.8023383351114638215-8.3060215023381284915i'), eps, C.loggamma, C('5-5i'));

  Jmat.Test.testFunction(0, eps, R.loggamma, 1);
  Jmat.Test.testFunction(2.25271265173, eps, R.loggamma, 0.1);
  Jmat.Test.testFunction(0.572364942925, eps, R.loggamma, 0.5);
  Jmat.Test.testFunction(6.63762397347e-02, eps, R.loggamma, 0.9);
  Jmat.Test.testFunction(-0.0982718364218, eps, R.loggamma, 1.25);
  Jmat.Test.testFunction(-0.0844011210205, eps, R.loggamma, 1.75);
  Jmat.Test.testFunction(0, eps, R.loggamma, 2);
  Jmat.Test.testFunction(0.693147180560, eps, R.loggamma, 3);
  Jmat.Test.testFunction(3.17805383035, eps, R.loggamma, 5);
  Jmat.Test.testFunction(1.28018274801e+01, eps, R.loggamma, 10);
  Jmat.Test.testFunction(3.93398841872e+01, eps, R.loggamma, 20);
  Jmat.Test.testFunction(3.59134205370e+02, eps, R.loggamma, 100);
  Jmat.Test.testFunction(5.90522042321e+03, eps, R.loggamma, 1000);
  Jmat.Test.testFunction(Infinity, eps, R.loggamma, Infinity);
  Jmat.Test.testFunction(1.20247578639, eps, R.loggamma, -1.3);
  Jmat.Test.testFunction(Infinity, eps, R.loggamma, -3);
  Jmat.Test.testFunction(-1.30900668499, eps, R.loggamma, -3.5);
  Jmat.Test.testFunction(-4.16245710210e+02, eps, R.loggamma, -111.5);
};

Jmat.Test.doUnitTestRealDistributions = function() {
  var R = Jmat.Real;
  var C = Jmat.Complex;
  var T = Jmat.Test;
  var that = this;
  var testDistribution = function(expected_pdf, expected_cdf, precision, precision_qf, name, var_arg) {
    var a = Array.prototype.slice.call(arguments).slice(5); /*var_arg*/
    var pdf_f = R['pdf_' + name];
    if(!pdf_f) pdf_f = R['pmf_' + name];
    var cdf_f = R['cdf_' + name];
    var qf_f = R['qf_' + name];

    T.testFunction.apply(that, [expected_pdf, precision, pdf_f].concat(a));
    T.testFunction.apply(that, [expected_cdf, precision, cdf_f].concat(a));
    var expected_qf = a[0];
    var qf_in = cdf_f.apply(that, a);
    var amod = a.slice(0);
    amod[0] = qf_in;
    T.testFunction.apply(that, [expected_qf, precision_qf, qf_f].concat(amod));

    // and now complex
    for(var i = 0; i < a.length; i++) a[i] = C(a[i]);
    expected_pdf = Complex(expected_pdf);
    expected_cdf = Complex(expected_cdf);
    pdf_f = C['pdf_' + name];
    if(!pdf_f) pdf_f = C['pmf_' + name];
    cdf_f = C['cdf_' + name];
    qf_f = C['qf_' + name];

    T.testFunction.apply(that, [expected_pdf, precision, pdf_f].concat(a));
    T.testFunction.apply(that, [expected_cdf, precision, cdf_f].concat(a));
    var expected_qf = a[0];
    var qf_in = cdf_f.apply(that, a);
    var amod = a.slice(0);
    amod[0] = qf_in;
    T.testFunction.apply(that, [expected_qf, precision_qf, qf_f].concat(amod));
  }

  var eps = 1e-10;
  var eps_qf = 1e-5;

  testDistribution(0.33333333333333, 0.33333333333333, eps, eps_qf, 'uniform', 2, 1, 4);
  testDistribution(0, 0, eps, Infinity, 'uniform', -1, 0, 1); // qf can't know what x was here, so its eps disabled
  testDistribution(0, 1, eps, Infinity, 'uniform', 5, 0, 2); // qf can't know what x was here, so its eps disabled

  testDistribution(1.486719514734297707908e-6, 2.86651571879193911674e-7, eps, eps_qf, 'standardnormal', -5);
  testDistribution(0.3520653267642994777747, 0.3085375387259868963623, eps, eps_qf, 'standardnormal', -0.5);
  testDistribution(0.39894228040143267794, 0.5, eps, eps_qf, 'standardnormal', 0);
  testDistribution(0.3520653267642994777747, 0.6914624612740131036377, eps, eps_qf, 'standardnormal', 0.5);
  testDistribution(1.486719514734297707908E-6, 0.9999997133484281208061, eps, eps_qf, 'standardnormal', 5);

  testDistribution(0.2515888184619954426786, 0.6305586598182363617272, eps, eps_qf, 'normal', 1, 0.5, 1.5);
  testDistribution(0.0381387815460524085608, 0.3820885778110473626935, eps, eps_qf, 'normal', 5, 8, 10);
  testDistribution(0.01713685920478073569636, 0.09680048458561033315201, eps, eps_qf, 'normal', -5, 8, 10);

  testDistribution(0.2515888184619954426786, 0.3694413401817636382728, eps, eps_qf, 'lognormal', 1, 0.5, 1.5);
  testDistribution(0.006505170494444215303193, 0.2613931832392658068243, eps, eps_qf, 'lognormal', 5, 8, 10);

  testDistribution(0.1591549430919, 0.75, eps, eps_qf, 'cauchy', 1, 0, 1);
  testDistribution(0.21220659078919, 0.5, eps, eps_qf, 'cauchy', 1.5, 1.5, 1.5);
  testDistribution(0.07639437268411, 0.79516723530087, eps, eps_qf, 'cauchy', 2.5, 0.5, 1.5);

  testDistribution(3.4633157382171e-4, 0.98949853491578, eps, eps_qf, 'studentt', 30.3, 1);
  testDistribution(3.5830637408089E-5, 0.99945627926805, eps, eps_qf, 'studentt', 30.3, 2);
  testDistribution(3.899041614152E-6, 0.99996051678271, eps, eps_qf, 'studentt', 30.3, 3);
  testDistribution(4.647805715373E-7, 0.99999646651138, eps, eps_qf, 'studentt', 30.3, 4);
  testDistribution(0.19245008972988, 0.78867513459481, eps, eps_qf, 'studentt', 1, 2);
  testDistribution(1.10015649E-7, 0.99999998054293, eps, eps_qf, 'studentt', 5.5, 10000);
  testDistribution(0.1176685042172, 0.87608177345685, eps, eps_qf, 'studentt', 1.5, 2.5);
  testDistribution(0.1176685042172, 0.12391822654315, eps, eps_qf, 'studentt', -1.5, 2.5);

  testDistribution(2.350110795843e-16, 5.198100275214e-17, eps, eps_qf, 'chi_square', 5, 50);
  testDistribution(0.207553748710, 0.427593295529, eps, eps_qf, 'chi_square', 2, 3);
  testDistribution(0.122041521349, 0.5841198130044, eps, eps_qf, 'chi_square', 5, 5);

  testDistribution(0.196611933241, 0.731058578630, eps, eps_qf, 'logistic', 1, 0, 1);
  testDistribution(3.105521163621e-6, 0.9999962733607, eps, eps_qf, 'logistic', 20, 5, 1.2);

  testDistribution(0.180447044315, 0.142876539501, eps, eps_qf, 'gamma', 2, 4, 1);
  testDistribution(0.096438800200, 0.554768248693, eps, eps_qf, 'gamma', 5.1, 1.7, 3.3);

  testDistribution(1.536, 0.1808, eps, eps_qf, 'beta', 0.2, 2, 3);
  testDistribution(0.16141880840106, 0.017820363993139, eps, eps_qf, 'beta', 0.77, 5, 0.1);

  testDistribution(0.641336686303, 0.132305485692, eps, eps_qf, 'fisher', 0.2, 2.5, 1.5);
  testDistribution(0.0044460064995, 0.97496898418, eps, eps_qf, 'fisher', 10, 1, 5);
  testDistribution(4.0989816415e-6, 0.99982905242425, eps, eps_qf, 'fisher', 100, 1, 5);

  testDistribution(1, 0, eps, eps_qf, 'weibull', 0, 1, 1);
  testDistribution(0.367879441171, 0.6321205588, eps, eps_qf, 'weibull', 1, 1, 1);
  testDistribution(Infinity, 0, eps, eps_qf, 'weibull', 0, 1, 0.5);
  testDistribution(0.183939720585, 0.6321205588, eps, eps_qf, 'weibull', 1, 1, 0.5);

  testDistribution(0.13533528323, 0.8646647167, eps, eps_qf, 'exponential', 2, 1);
  testDistribution(0.066887171122, 0.632120558828, eps, eps_qf, 'exponential', 5.5, 1 / 5.5);

  testDistribution(0.194700195767851, 0.610599608464297, eps, eps_qf, 'laplace', 0.5, 0, 2);
  testDistribution(0.134064009207127, 0.335160023017819, eps, eps_qf, 'laplace', 0.5, 1.5, 2.5);

  testDistribution(0.9, 0.9, eps, eps_qf, 'bernoulli', 0, 0.1);
  testDistribution(0.1, 1, eps, eps_qf, 'bernoulli', 1, 0.1);
  testDistribution(0.5, 0.5, eps, eps_qf, 'bernoulli', 0, 0.5);
  testDistribution(0.5, 1, eps, eps_qf, 'bernoulli', 1, 0.5);
  testDistribution(0.1, 0.1, eps, eps_qf, 'bernoulli', 0, 0.9);
  testDistribution(0.9, 1, eps, eps_qf, 'bernoulli', 1, 0.9);

  testDistribution(0.0014880348, 0.9998530974, eps, eps_qf, 'binomial', 5, 10, 0.1);
  testDistribution(0.0014880348, 0.0016349374, eps, eps_qf, 'binomial', 5, 10, 0.9);
  testDistribution(0.03125, 1, eps, eps_qf, 'binomial', 5, 5, 0.5);

  testDistribution(0.121714452424092, 0.158597619825332, eps, eps_qf, 'poisson', 1, 3.3);
  // precision for qf not tested below (by setting its eps to Infinity), it is low, partially because all high k makes the value very close to 1, so any tiny change in that value affects the inverse result enormously
  testDistribution(0.175467369767851, 0.615960654833063, eps, Infinity, 'poisson', 5, 5);
  testDistribution(1.0137771196303e-7, 0.999999989952234, eps, Infinity, 'poisson', 10, 1);
  testDistribution(2.80684572494333e-108, 1, eps, Infinity, 'poisson', 100, 3.3);
}

// throws on fail, prints 'success' on success
Jmat.doUnitTest = function(opt_verbose) {
  // check that the test framework itself can actually fail --> disabled because annoying when using "pause on exception"
  /*var thrown = false;
  try {
    Jmat.Test.testFunction(3, 0, Jmat.add, 1, 1);
  } catch(error) {
    thrown = true; // this is expected
  }
  if(!thrown) throw 'that should have thrown error!';*/

  Jmat.Test.initAccuracy();

  var eps = 1e-10;
  // basic operators
  Jmat.Test.testFunction(5, eps, Jmat.add, 2, 3);
  Jmat.Test.testFunction(-1, eps, Jmat.sub, 2, 3);
  Jmat.Test.testFunction(6, eps, Jmat.mul, 2, 3);
  Jmat.Test.testFunction(0.666666666666666666, eps, Jmat.div, 2, 3);
  Jmat.Test.testFunction(8, eps, Jmat.pow, 2, 3);
  Jmat.Test.testFunction(0.20787957635076190854695561, eps, Jmat.pow, 'i', 'i');

  // advanced functions
  Jmat.Test.testFunction(0, eps, Jmat.sin, Math.PI);
  Jmat.Test.testFunction(-1, eps, Jmat.cos, Math.PI);

  // special functions
  Jmat.Test.testFunction(24, eps, Jmat.gamma, 5);
  Jmat.Test.testFunction('-0.15494982830181-0.498015668118356i', eps, Jmat.gamma, 'i');
  Jmat.Test.testFunction('0.9303796037430951+0.0389361908951213i', 1e-5, Jmat.erf, '5+5i');
  Jmat.Test.testFunction('0.2074861066333588576972787235', 1e-9, Jmat.besselj, '10', '10');
  Jmat.Test.testFunction('0.2068008998147143416959879887', 1e-9, Jmat.besselj, '10.1', '10.1');
  Jmat.Test.testFunction('9.59012e-135-9.59012e-135i', 1e-6, Jmat.besselj, '112.5', '-5.5i');
  Jmat.Test.testFunction('-0.359814152183402722051986577', 1e-12, Jmat.bessely, '10', '10');
  Jmat.Test.testFunction('6.3618456410625559136428432181', 1e-12, Jmat.hypergeometric1F1, 1, 2, 3);
  Jmat.Test.testFunction('0.506370', 1e-5, Jmat.hypergeometric2F1, 0.5, 0.5, 0.5, -2.9);
  Jmat.Test.testFunction('0.493865', 1e-5, Jmat.hypergeometric2F1, 0.5, 0.5, 0.5, -3.1);
  Jmat.Test.doUnitTestLogGamma();

  //other
  Jmat.Test.testFunction(3581, 0, Real.smallestPrimeFactor, 12830723);

  // mod
  Jmat.Test.testFunction(0, 0, Jmat.mod, -6, 3);
  Jmat.Test.testFunction(1, 0, Jmat.mod, -5, 3);
  Jmat.Test.testFunction(2, 0, Jmat.mod, -4, 3);
  Jmat.Test.testFunction(0, 0, Jmat.mod, -3, 3);
  Jmat.Test.testFunction(1, 0, Jmat.mod, -2, 3);
  Jmat.Test.testFunction(2, 0, Jmat.mod, -1, 3);
  Jmat.Test.testFunction(0, 0, Jmat.mod, 0, 3);
  Jmat.Test.testFunction(1, 0, Jmat.mod, 1, 3);
  Jmat.Test.testFunction(2, 0, Jmat.mod, 2, 3);
  Jmat.Test.testFunction(0, 0, Jmat.mod, 3, 3);
  Jmat.Test.testFunction(1, 0, Jmat.mod, 4, 3);
  Jmat.Test.testFunction(2, 0, Jmat.mod, 5, 3);
  Jmat.Test.testFunction(0, 0, Jmat.mod, 6, 3);

  Jmat.Test.testFunction(0, 0, Jmat.mod, -6, -3);
  Jmat.Test.testFunction(-2, 0, Jmat.mod, -5, -3);
  Jmat.Test.testFunction(-1, 0, Jmat.mod, -4, -3);
  Jmat.Test.testFunction(-0, 0, Jmat.mod, -3, -3);
  Jmat.Test.testFunction(-2, 0, Jmat.mod, -2, -3);
  Jmat.Test.testFunction(-1, 0, Jmat.mod, -1, -3);
  Jmat.Test.testFunction(-0, 0, Jmat.mod, 0, -3);
  Jmat.Test.testFunction(-2, 0, Jmat.mod, 1, -3);
  Jmat.Test.testFunction(-1, 0, Jmat.mod, 2, -3);
  Jmat.Test.testFunction(-0, 0, Jmat.mod, 3, -3);
  Jmat.Test.testFunction(-2, 0, Jmat.mod, 4, -3);
  Jmat.Test.testFunction(-1, 0, Jmat.mod, 5, -3);
  Jmat.Test.testFunction(-0, 0, Jmat.mod, 6, -3);

  // distributions
  Jmat.Test.testFunction(0.274997, 1e-2, Jmat.qf_chi_square, 0.4, 1); // This one is very imprecise currently :(
  Jmat.Test.testFunction(0.198964, 1e-6, Jmat.pdf_studentt, 0.5, 0.5); // This one is very imprecise currently :(
  Jmat.Test.doUnitTestRealDistributions();

  // matrix basic operators
  Jmat.Test.testFunction([[6,8],[10,12]], eps, Jmat.add, [[1,2],[3,4]], [[5,6],[7,8]]);
  Jmat.Test.testFunction([[19,22],[43,50]], eps, Jmat.mul, [[1,2],[3,4]], [[5,6],[7,8]]);

  // matrix advanced operators
  Jmat.Test.testFunction([[-2,1],[1.5,-0.5]], eps, Jmat.inv, [[1,2],[3,4]]);
  Jmat.Test.testFunction([[5,-1],[-2,0]], eps, Jmat.fft, [[1,2],[3,4]]);
  Jmat.Test.testSVD([[0.404554, 0.914514], [0.914514, -0.404554]],
                    [[5.46499, 0],[0,0.365966]],
                    [[0.576048, -0.817416],[0.817416, 0.576048]],
                    1e-5, [[1,2],[3,4]]);
  Jmat.Test.testSVD([[1]],
                    [[2.23607, 0]],
                    [[0.447214, -0.894427], [0.894427, 0.447214]],
                    1e-5, [[1,2]]);
  Jmat.Test.testSVD([[0.447214, -0.894427], [0.894427, 0.447214]],
                    [[2.23607], [0]],
                    [[1]],
                    1e-5, [[1],[2]]);
  Jmat.Test.testSVD([[0.214837, 0.887231, -0.408248], [0.520587, 0.249644, 0.816497], [0.826338, -0.387943, -0.408248]],
                    [[16.8481, 0, 0],[0, 1.06837, 0],[0,0,0]],
                    [[0.479671, -0.776691, 0.408248],[0.572368, -0.0756865, -0.816497], [0.665064, 0.625318, 0.408248]],
                    1e-5, [[1,2,3],[4,5,6],[7,8,9]]);
  Jmat.Test.testEIG([1], [[1]], 1e-5, [[1]]);
  Jmat.Test.testEIG([5.37228, -0.372281], [[0.457427, -1.45743],[1, 1]], 1e-5, [[1,2],[3,4]]);
  Jmat.Test.testEIG([16.1168, -1.11684, 0], [[0.283349, -1.28335, 1],[0.641675, -0.141675, -2], [1, 1, 1]], 1e-4, [[1,2,3],[4,5,6],[7,8,9]]); //wolfram|alpha only gave 4 digits
  Jmat.Test.testEIG([80.1273, -5.15639, -1.18974, -0.781182],
                    [[0.0562274,-0.277836,-0.258731,-0.767583],[0.274212,-0.449235,1.15276,1.77101],[0.615571,-0.559067,-1.92247,-2.08569],[1,1,1,1]],
                    1e-4,
                    [[0,1,2,3],[5,7,11,13],[17,19,23,29],[31,37,41,43]]);
  Jmat.Test.testEIG([-23.6091, '11.3789+5.15137i', '11.3789-5.15137i', 0.851391],
                    [[0.100247,'-0.243-0.04689i','-0.243+0.04689i',82.8674],[-0.148667,'1.3376-0.39635i','1.33764+0.39635i',76.2564],[0.390964,'-0.4334-0.6958i','-0.4334+0.6958i', -1.35194],[1,1,1,1]],
                    1e-3,
                    [[0,1,2,-3],[-4,5,-6,7],[-8,9,10,-11],[-12,13,-14,-15]]);
  Jmat.Test.testEIG([3.38923, 1.85542, 0.755356],
                    [[0.268262/0.871621, -0.223801/0.470604, 0.936989/-0.137143], [0.410258/0.871621, -0.85349/0.470604, -0.321315/-0.137143], [1, 1, 1]],
                    1e-5,
                    [[1,0.5,0.5],[0.5,2,0.5],[0.5,0.5,3]]);
  Jmat.Test.testEIG([3.38923, 1.85542, 0.755356],
                    [[-0.871621/-0.268262,0.470604/-0.223801,-0.137143/0.936989],[-0.410258/-0.268262,-0.85349/-0.223801,-0.321315/0.936989],[1,1,1]],
                    1e-5,
                    [[3,0.5,0.5],[0.5,2,0.5],[0.5,0.5,1]]);
  Jmat.Test.testEIG([-2, 0, 0],
                    [[-1,0,0],[3,1,1],[1,0,0]],
                    1e-8,
                    [[-1,0,1],[3,0,-3],[1,0,-1]]);
  Jmat.Test.expectNear(Jmat.Matrix([[0.5],[-1],[0.5]]), Jmat.Matrix.solve(Jmat.Matrix([[0,1,2],[3,5,7],[11,13,17]]), Jmat.Matrix([[0],[0],[1]])), eps);

  // matrix parsing
  Jmat.Test.testFunction([[1,2],[3,4]], eps, Jmat.Matrix.parse, '[[1,2],[3,4]]');
  Jmat.Test.testFunction([[1,2,3,4]], eps, Jmat.Matrix.parse, '[[1,2,3,4]]');
  Jmat.Test.testFunction([[1],[2],[3],[4]], eps, Jmat.Matrix.parse, '[1,2,3,4]');
  Jmat.Test.testFunction([[1],[2],[3],[Complex(0, 4)]], eps, Jmat.Matrix.parse, '[1,2,3,4i]');

  // numerical algorithms
  Jmat.Test.testFunction(333.33333333333, eps, Jmat.integrate, 0, 10, function(z) { return z.mul(z); });

  // quaternions
  Jmat.Test.testFunction('-28+4i+6j+8k', eps, Jmat.mul, '1+2i+3j+4k', '1+2i+3j+4k');
  Jmat.Test.testFunction('-28+4i+6j+8k', eps, Jmat.pow, '1+2i+3j+4k', '2+0i+0j+0k');
  Jmat.Test.testFunction('0.778283161474+0.281490511227i+0.281490511227j+0.281490511227k', eps, Jmat.lambertw, '1+i+j+k');

  // bigints
  Jmat.Test.testFunction('40094690950920881030683735292761468389214899724061', 0, Jmat.div,
                         '1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139',
                         '37975227936943673922808872755445627854565536638199');
  Jmat.Test.testFunction(BigInt('32769'), 0, Jmat.div, BigInt('4653329978'), BigInt('142001'));
  Jmat.Test.testFunction('1250000000', 0, Jmat.div, '50000000000', '40');
  Jmat.Test.testFunction('20000000000', 0, BigInt.sqrt, '400000000000000000000');
  Jmat.Test.testFunction('5', 0, BigInt.log2, '63');
  Jmat.Test.testFunction('6', 0, BigInt.log2, '64');
  Jmat.Test.testFunction('30414093201713378043612608166064768844377641568960512000000000000', 0, BigInt.factorial, 50);
  Jmat.Test.expectTrue(BigInt.log2(BigInt.fromInt(63, 2)).toString() == '5');
  Jmat.Test.expectTrue(BigInt.log2(BigInt.fromInt(64, 2)).toString() == '6');
  Jmat.Test.expectTrue(BigInt.log2(BigInt.fromInt(63, 4)).toString() == '5');
  Jmat.Test.expectTrue(BigInt.log2(BigInt.fromInt(64, 4)).toString() == '6');
  Jmat.Test.expectTrue(BigInt.isPrime('671998030559713968361666935769'));
  Jmat.Test.expectFalse(BigInt.isPrime('19923108241787117701'));

  // 0d matrix generation
  var vec1 = Jmat.Matrix.parse('[1]');
  Jmat.Test.expectNear(vec1, Jmat.Matrix.make([1]), eps);
  Jmat.Test.expectNear(Jmat.sub(vec1,vec1), Jmat.Matrix.make([0]), eps);

  // 1d matrix generation
  var vec01 = Jmat.Matrix.parse('[0,1]');
  Jmat.Test.expectNear(vec01, Jmat.Matrix.make([0,1]), eps);
  Jmat.Test.expectNear(vec01, Jmat.Matrix.subcol(Jmat.Matrix.make([[2,0],[0,1]]), 0), eps);

  Jmat.Test.expectNear(vec01.transpose(), Jmat.Matrix.make([[0,1]]), eps);
  Jmat.Test.expectNear(vec01.transpose(), Jmat.Matrix.subrow(Jmat.Matrix.make([[2,0],[0,1]]), 0), eps);

  Jmat.Test.finalizeAccuracy(opt_verbose);
  console.log('success');
  return 'success';
};
