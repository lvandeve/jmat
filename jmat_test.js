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

// Works both for Complex or Matrix objects.
Jmat.Test.expectNear = function(e, a, epsilon) {
  if(Jmat.isNaN(e) && Jmat.isNaN(a)) return; //both NaN is ok for test
  if(!Jmat.near(e, a, epsilon)) throw 'fail: expected ' + Jmat.toString(e) + ' got ' + Jmat.toString(a);
};

// Expect that the result of the mathematical function f with the arguments of var_arg, is near the expected result. Some numerical intolerance is allowed.
Jmat.Test.testFunction = function(expected, epsilon, f, var_arg) {
  var result = f.apply(this, Array.prototype.slice.call(arguments).slice(3) /*var_arg*/);
  Jmat.Test.expectNear(expected, result, epsilon);
};

//u,s,v = expected values
//m = input matrix
Jmat.Test.testSVD = function(u, s, v, epsilon, m) {
  var svd = Matrix.svd(Matrix.cast(m));
  Jmat.Test.expectNear(u, svd.u, epsilon);
  Jmat.Test.expectNear(s, svd.s, epsilon);
  Jmat.Test.expectNear(v, svd.v, epsilon);
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
  Jmat.Test.testFunction(6, eps, Jmat.mul, 2, 3);

  // advanced functions
  Jmat.Test.testFunction(0, eps, Jmat.sin, Math.PI);

  // special functions
  Jmat.Test.testFunction(24, eps, Jmat.gamma, 5);
  Jmat.Test.testFunction('-0.15494982830181-0.498015668118356i', eps, Jmat.gamma, 'i');
  Jmat.Test.testFunction('0.9303796037430951+0.0389361908951213i', 1e-5, Jmat.erf, '5+5i');

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

  // matrix parsing
  Jmat.Test.testFunction([[1,2],[3,4]], eps, Jmat.Matrix.parse, '[[1,2],[3,4]]');
  Jmat.Test.testFunction([[1,2,3,4]], eps, Jmat.Matrix.parse, '[[1,2,3,4]]');
  Jmat.Test.testFunction([[1],[2],[3],[4]], eps, Jmat.Matrix.parse, '[1,2,3,4]');

  // numerical algorithms
  Jmat.Test.testFunction(333.33333333333, eps, Jmat.integrate, 0, 10, function(z) { return z.mul(z); });
  

  console.log('success');
  return 'success';
};
