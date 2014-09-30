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
*/

// Unit tests for Jmat.js

// constructor
Jmat.Test = function() {
  // empty, it's a namespace
};

Jmat.Test.expectTrue = function(value, opt_message) {
  var message = opt_message ? ('fail: ' + opt_message) : 'fail';
  if (!value) throw message;
}

// Works both for Complex or Matrix objects.
// Precision is number of decimal digits that should match
Jmat.Test.expectNear = function(e, a, precision) {
  if(Jmat.isNaN(e) && Jmat.isNaN(a)) return; //both NaN is ok for test
  if(Jmat.eq(e, 0) || Jmat.eq(a, 0) || Jmat.matrixIn_(e) || Jmat.matrixIn_(a)) {
    // for 0, allow the other to be absolute rather than relative near it. For matrices, always use absolute epsilon as well.
    if(!Jmat.near(e, a, precision)) throw 'fail: expected ' + Jmat.toString(e) + ' got ' + Jmat.toString(a) + '. Expected precision: ' + precision;
  } else {
    if(!Jmat.relnear(e, a, precision)) throw 'fail: expected ' + Jmat.toString(e) + ' got ' + Jmat.toString(a) + '. Expected precision: ' + precision;
  }
};

// Expect that the result of the mathematical function f with the arguments of var_arg, is near the expected result. Some numerical intolerance is allowed.
Jmat.Test.testFunction = function(expected, epsilon, f, var_arg) {
  var result = f.apply(this, Array.prototype.slice.call(arguments).slice(3) /*var_arg*/);
  Jmat.Test.expectNear(expected, result, epsilon);
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

//l,v = expected values
//m = input matrix
Jmat.Test.testEIG = function(l, v, epsilon, m) {
  var eig = Matrix.eig(Matrix.cast(m));
  Jmat.Test.expectNear(l, eig.l, epsilon);
  Jmat.Test.expectNear(v, eig.v, epsilon);
};

// throws on fail, prints 'success' on success
Jmat.doUnitTest = function() {
  // check that the test framework itself can actually fail
  var thrown = false;
  try {
    Jmat.Test.testFunction(3, eps, Jmat.add, 1, 1);
  } catch(error) {
    thrown = true; // this is expected
  }
  if(!thrown) throw 'that should have thrown error!';

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
  Jmat.Test.testEIG([[1]], [[1]], 1e-5, [[1]]);
  Jmat.Test.testEIG([[5.37228], [-0.372281]], [[0.457427, -1.45743],[1, 1]], 1e-5, [[1,2],[3,4]]);
  Jmat.Test.testEIG([[16.1168], [-1.11684], [0]], [[0.283349, -1.28335, 1],[0.641675, -0.141675, -2], [1, 1, 1]], 1e-4, [[1,2,3],[4,5,6],[7,8,9]]); //wolfram|alpha only gave 4 digits

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

  // bignums
  Jmat.Test.testFunction('40094690950920881030683735292761468389214899724061', 0, Jmat.div, '1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139', '37975227936943673922808872755445627854565536638199');
  Jmat.Test.testFunction('1250000000', 0, Jmat.div, '50000000000', '40');
  Jmat.Test.testFunction('20000000000', 0, Jmat.BigNum.sqrt, '400000000000000000000');
  Jmat.Test.testFunction('5', 0, Jmat.BigNum.log2, '63');
  Jmat.Test.testFunction('6', 0, Jmat.BigNum.log2, '64');
  Jmat.Test.expectTrue(Jmat.BigNum.log2(Jmat.BigNum.fromInt(63, 2)).toString() == '5');
  Jmat.Test.expectTrue(Jmat.BigNum.log2(Jmat.BigNum.fromInt(64, 2)).toString() == '6');
  Jmat.Test.expectTrue(Jmat.BigNum.log2(Jmat.BigNum.fromInt(63, 4)).toString() == '5');
  Jmat.Test.expectTrue(Jmat.BigNum.log2(Jmat.BigNum.fromInt(64, 4)).toString() == '6');

  console.log('success');
  return 'success';
};
