/*
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
function makeElement(parent, tag) {
  var el =  document.createElement(tag);
  parent.appendChild(el);
  return el;
}

function makeAbsElement(px, py, parent, tag) {
  var el =  document.createElement(tag);
  el.style.position = 'absolute';
  el.style.left = '' + Math.floor(px) + 'px';
  el.style.top = '' + Math.floor(py) + 'px';
  parent.appendChild(el);
  return el;
}

function makeRelElement(px, py, parent, tag) {
  var el =  document.createElement(tag);
  el.style.position = 'relative';
  el.style.left = '' + Math.floor(px) + 'px';
  el.style.top = '' + Math.floor(py) + 'px';
  parent.appendChild(el);
  return el;
}

function makeRelDivAt(px, py, parent) {
  return makeRelElement(px, py, parent, 'div');
}

function makeBlockDiv(parent) {
  var el =  document.createElement('div');
  parent.appendChild(el);
  return el;
}

function makeBlockDivAt(px, py, parent) {
  var el = makeBlockDiv(parent);
  el.style.paddingLeft = '' + Math.floor(px) + 'px';
  el.style.paddingTop = '' + Math.floor(py) + 'px';
  parent.appendChild(el);
  return el;
}

function makeRelDiv(parent) {
  var el =  document.createElement('div');
  el.style.position = 'relative';
  parent.appendChild(el);
  return el;
}

function makeDropDown(x, y, options, parent) {
  var sel = makeAbsElement(x, y, parent, 'select');
  for(var i = 0; i < options.length; i++) {
    makeElement(sel, 'option').innerHTML = options[i];
  }
  return sel;
}

function makeLabeledDropDown(x, y, label, options, parent) {
  var div = makeAbsElement(x, y, parent, 'div');
  div.innerHTML = label;
  return makeDropDown(x, y + 16, options, parent);
}

var fun1d = [
    '--',
    'ei', 'e1', 'li', 'si', 'ci', 'shi', 'chi',
    'gamma', 'digamma', 'trigamma', 'loggamma', 'gamma_inv',
    'airy', 'bairy', 'airy_deriv', 'bairy_deriv',
    'erf', 'erfc', 'erf_inv', 'erfc_inv', 'erfi', 'dawson', 'faddeeva',
    'fresnels', 'fresnelc',
    'zeta', 'eta', 'dilog', 'trilog',
    'exp', 'log', 'sqrt',
    'lambertw', 'lambertwm', 'tetrational',
    'sin', 'cos', 'tan',
    'asin', 'acos', 'atan',
    'sinh', 'cosh', 'tanh',
    'asinh', 'acosh', 'atanh',
    'sinc',
    'abs', 'arg', 'frac', 'fracn',
    'bitneg',
    'minkowski',
    'x', 'inv',
];
var fun2d = [
    '--',
    'struveh', 'struvek', 'struvel', 'struvem', 'angerj', 'webere',
    'gamma_p', 'gamma_q', 'incgamma_lower', 'incgamma_upper', 'polygamma',
    'beta',
    'besselj', 'bessely', 'besseli', 'besselk', 'hankelh1', 'hankelh2',
    'lambertwb',
    'polylog', 'hurwitzzeta',
    'binomial', 'permutation', 'stirling2',
    'theta1', 'theta2', 'theta3', 'theta4',
    'hypergeometric0F1',
    'atan2',
    'agm', 'ghm',
    'tetration',
    'mod', 'rem',
    'bitand', 'bitor', 'bitxor', 'lshift', 'rshift',
    'add', 'sub', 'mul', 'div', 'pow', 'logy',
];

var headerEl = makeBlockDiv(document.body);
headerEl.innerHTML = '<h2>Jmat.js</h2>';
headerEl.style.backgroundColor = '#2b6';
headerEl.style.padding = '2px';
headerEl.style.textAlign = 'center';
headerEl.style.color = 'white';

var bodyEl = makeRelDivAt(0, 0, document.body);
bodyEl.style.margin = '50px';

var el = makeBlockDivAt(0, 0, bodyEl);
el.innerHTML = 'Jmat.js is a mathematics library in JavaScript, supporting complex special functions, matrices and statistical distributions. This demo allows plotting several of the Jmat functions, and evaluating arbitrary Jmat code.';
el.innerHTML += '<p/><a href="https://github.com/lvandeve/jmat">GitHub page</a>';
var plotTitle = makeBlockDiv(bodyEl);
plotTitle.innerHTML = '<h2>Plot demo</h2>';
plotTitle.style.margin = '5px';

var plotContainerEl = makeRelDivAt(15, 0, bodyEl);
plotContainerEl.style.width = '400px';
plotContainerEl.style.height = '400px';

var plotDropdownsBlock = makeBlockDiv(bodyEl);

var plotDropdowns = makeRelDiv(plotDropdownsBlock);
plotDropdowns.style.height = '70px';

var el0 = makeLabeledDropDown(0, 0, 'real', fun1d, plotDropdowns);
el0.onchange = function() {
  Jmat.stopPlotting();
  var f = fun1d[el0.selectedIndex];
  if(Jmat[f]) Jmat.plotReal(Jmat[f], plotContainerEl, {p:1}, f);
  else if(Jmat.Complex[f]) Jmat.plotReal(function(z) { return Jmat.Complex[f](Jmat.Complex.cast(z)); }, plotContainerEl, {p:1}, f);
  else console.log('function ' + f + ' not found');
};

var el1 = makeLabeledDropDown(130, 0, 'complex', fun1d, plotDropdowns);
el1.onchange = function() {
  Jmat.stopPlotting();
  var f = fun1d[el1.selectedIndex];
  if(Jmat[f]) Jmat.plotComplex(Jmat[f], plotContainerEl, {p:1}, f);
  else if(Jmat.Complex[f]) Jmat.plotComplex(function(z) { return Jmat.Complex[f](Jmat.Complex.cast(z)); }, plotContainerEl, {p:1}, f);
  else console.log('function ' + f + ' not found');
};

var el2 = makeLabeledDropDown(260, 0, '2d', fun2d, plotDropdowns);
el2.onchange = function() {
  Jmat.stopPlotting();
  var f = fun2d[el2.selectedIndex];
  if(Jmat[f]) Jmat.plot2D(Jmat[f], plotContainerEl, {p:1}, f);
  else if(Jmat.Complex[f]) Jmat.plot2D(function(x, y) { return Jmat.Complex[f](Jmat.Complex.cast(x), Jmat.Complex.cast(y)); }, plotContainerEl, {p:1}, f);
  else console.log('function ' + f + ' not found');
};

Jmat.plotComplex(Jmat['gamma'], plotContainerEl, {p:1}, 'gamma');

var plotInfo = makeBlockDiv(bodyEl);


var info = makeRelDivAt(0, 0, plotInfo);
info.innerHTML =
    'Plot types: "Real": function of 1 argument plotted with classic X/Y real plot. "Complex": function of 1 argument plotted in 2D for complex arguments with complex color wheel. "2D": function of 2 arguments plotted for real arguments in 2D with result as complex color wheel.<p/>' +
    'Complex color wheel legend: red = positive, cyan = negative, black = zero, white = infinity, darker = lower abs, lighter = higher abs, acid green = positive imaginary, purple = negative imaginary, other hues = other complex arguments, grey = NaN<p/>';


var evalTitle = makeBlockDiv(bodyEl);
evalTitle.innerHTML = '<h2>Eval demo</h2>';
evalTitle.style.margin = '5px';

var evalContainer = makeBlockDiv(bodyEl);

var area = makeElement(evalContainer, 'textarea');
area.style.width = '700px';
area.style.height = '50px';
area.value = 'Jmat.fft([[1,2],[3,4]]).toString()';
var button = makeRelDivAt(0, 5, evalContainer);
button.innerHTML = 'eval';
button.style.textAlign = 'center';
button.style.backgroundColor = '#eee';
button.style.border = '1px solid black';
button.style.width = '50px';
button.style.height = '20px';
var result = makeBlockDivAt(0, 10, evalContainer);
result.innerHTML = '_';
result.style.width = '500px';
button.onclick = function() {
  result.innerHTML = 'answer: ' + eval(area.value);
};
var examples = makeBlockDivAt(0, 10, evalContainer);
examples.style.color = '#aaa';
examples.style.width = '1000px';
examples.style.fontSize = 'small';
examples.innerHTML = 'Examples:<br/>' +
    '<ul>' +
    '<li>Complex(\'1+2i\').add(Complex(\'3+4i\')) </li>' +
    '<li>Jmat.besselj(5, Complex(5.5, 1)) </li>' +
    '<li>Jmat.fft([[1,2],[3,4]]).toString()</li>' +
    '<li>Jmat.eig([[1,2],[3,4]]).l.toString()</li>' +
    '<li>Jmat.gamma(5.5).add(Jmat.trigamma(5.5)) </li>' +
    '<li>Jmat.factorize(30030) </li>' +
    '<li>Jmat.plot2D(Jmat.gamma_p, plotContainerEl); </li>' +
    '<li>Jmat.plotComplex(function(z) { return Jmat.polygamma(4, z); }, plotContainerEl, {p:1}); </li>' +
    '<li>Jmat.plotReal(function(x) { return Jmat.cdf_studentt(x, 2); }, plotContainerEl, {p:1}); </li>' +
    '<li>Jmat.plotComplex(function(z) { return Jmat.theta2(z, 0.2); }, plotContainerEl, {p:1}); </li>' +
    '<li>Jmat.plotReal(function(x) { return Jmat.pdf_laplace(x, 0, 1); }, plotContainerEl, {xsize:10, ysize:2});  </li>' +
    '<li>Jmat.plotComplex(function(c) { var i = 0; var z = Complex(0); for(;;) { if(z.abs() > 2) break; z = z.mul(z).add(c); i++; if(i > 60) return Complex(0); } return Complex.polar(1, i * Math.PI / 60); }, plotContainerEl, {p:1, s:4}); </li>' +
    '</ul>' +
    'ops: add, sub, mul, div, pow, sqrt, cos, sign, ceil, zeta, ... <br/>' +
    'matrix ops: determinant, eig, evd, rank, trace, norm, norm2, svd, qr, pseudoinverse, ... <br/>' +
    'stat ops: pdf_, cdf_ or qf_ + , normal, cauchy, studentt, chi_square, gamma, ... <br/>';


B = Jmat.BigInt;
C = Jmat.Complex;
Q = Jmat.Quaternion;
M = Jmat.Matrix;
R = Jmat.Real;

info = makeBlockDivAt(0, 10, bodyEl);
info.innerHTML = 'This software is provided \'as-is\', without any express or implied warranty. In no event will the authors be held liable for any damages arising from the use of this software.<p/>' +
    'Copyright (c) 2011-2016 by Lode Vandevenne.';
