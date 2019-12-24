/*
Jmat.js

Copyright (c) 2011-2019, Lode Vandevenne
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

// NOTE: requires charset="utf-8" to render navigation arrow characters correctly

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Plotting - public API
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// TODO: support zooming in/out with mousewheel (around mouse cursor as center)
// TODO: support dragging to move
// TODO: support showing of max value of function or auto scale to it
// TODO: for the 'real' plot (where real refers to the input value in fact), have choice between showing complex values as abs/arg with color, or separate re/im

// Jmat graphics, graphing and plotting library
Jmat.Plot = function() {
};

// e.g.: new Jmat.Plot.Params({p: 1}), see the the field assignments for each of the possible parameter names
// constructor
Jmat.Plot.Params = function(o) {
  if(!o) o = {};

  this.p = 2; //pixel size (==> resolution). 1 = highest resolution, 2 = half resolution, etc...

  // xsize and ysize define the zoom area
  this.xsize = 20; // e.g. an xsize of 10 makes it go from x=-5 to +5 if shift is 0
  this.ysize = 20;
  this.xshift = 0; // point which you want in the center of the plot
  this.yshift = 0;

  // these are only used for 2D plots
  this.xshift_im = 0;
  this.yshift_im = 0;
  this.twiddle = false; // when zoomed far out, use half integers instead of integers for pixels

  // The value at which the complex color wheel has highest saturation (pure red for positive real). Higher value gives more white color, lower gives darker color.
  // Not used by real plot.
  this.v = 1;

  // size of the plot area (excluding controls and labels)
  this.w = 640;
  this.h = 640;

  if (o) Jmat.Plot.copyParams_(o, this);
};

//plots a real plot with 1 input and 1 output with x/y axes (still works with Jmat.Complex objects, and uses color to indicate complex output)
//fun = mathematical function taking 1 Jmat.Complex argument, e.g. Jmat.sin
//params = parameter object with plot size, resolution, etc.... See Jmat.Plot.Params.
//parent = HTML parent element, best of type div. All necessary elements (e.g. canvas) will be created inside of it.
Jmat.Plot.plotReal = function(fun, label, params, parent) {
  params = new Jmat.Plot.Params(params || Jmat.Plot.defaultParams); // always make copy, since zooming etc... may change the object
  if(!parent) parent = Jmat.Plot.defaultParent || document.body;
  Jmat.Plot.plotReal_(fun, params, parent, label);
};

//plots a complex function with 1 input and 1 output using domain coloring
//fun = mathematical function taking 1 Jmat.Complex argument, e.g. Jmat.gamma
//params = parameter object with plot size, resolution, etc.... See Jmat.Plot.Params.
//parent = HTML parent element, best of type div. All necessary elements (e.g. canvas) will be created inside of it.
Jmat.Plot.plotComplex = function(fun, label, params, parent) {
  params = new Jmat.Plot.Params(params || Jmat.Plot.defaultParams); // always make copy, since zooming etc... may change the object
  if(!parent) parent = Jmat.Plot.defaultParent || document.body;
  Jmat.Plot.plot2D_(function(x, y) {
    return fun(Jmat.Complex(x.re, y.re));
  }, params, parent, label, 're', 'im', 1);
};

//plots a complex function with 2 inputs and 2 outputs using domain coloring
//uses 2 reals as input by default, which can be shifted to get imaginary parts, so makes a 2D slice through what is a function with 4D input (2 complex inputs)
//fun = mathematical function taking 2 Jmat.Complex arguments, e.g. Jmat.besselj
//params = parameter object with plot size, resolution, etc.... See Jmat.Plot.Params.
//parent = HTML parent element, best of type div. All necessary elements (e.g. canvas) will be created inside of it.
Jmat.Plot.plot2D = function(fun, label, params, parent) {
  params = new Jmat.Plot.Params(params || Jmat.Plot.defaultParams); // always make copy, since zooming etc... may change the object
  if(!parent) parent = Jmat.Plot.defaultParent || document.body;
  Jmat.Plot.plot2D_(fun, params, parent, label, undefined, undefined, 2);
};

// Stop plot2D or plotComplex, if they are taking a very long time. May still
// calculate several pixels before actually stopping: they do several lines at the time.
Jmat.Plot.stopPlotting = function() {
  Jmat.Plot.stopIndex_ = (Jmat.Plot.stopIndex_ ? Jmat.Plot.stopIndex_ + 1 : 1);
};


// If set, and no params given, this can override the default
Jmat.Plot.defaultParams = undefined;

// If set, and no DOM element given, this can override the default
Jmat.Plot.defaultParent = undefined;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Graphics (color & DOM) helper functions
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Jmat.Plot.copyParams_ = function(o, dest) {
  if(!o) return;

  for(var k in o) {
    if(!o.hasOwnProperty(k)) continue;
    dest[k] = o[k];
  }
};

// pixel screen coordinates to mathematical coordinates
Jmat.Plot.fromPixel = function(params, px, py) {
  py = params.h - py; // doesn't need -1 due to dividing through params.h instead of (params.h + 1) further on. TODO: check which we really want.
  var dx = px - params.w / 2;
  var dy = py - params.h / 2;
  var x = params.xshift + dx * params.xsize / params.w;
  var y = params.yshift + dy * params.ysize / params.h;
  if(params.xsize >= params.w) x = Math.floor(x) + (params.twiddle ? 0.5 : 0);
  if(params.ysize >= params.h) y = Math.floor(y) + (params.twiddle ? 0.5 : 0);

  return [Jmat.Complex(x, params.xshift_im), Jmat.Complex(y, params.yshift_im)];
};

//input and output in range [0-255]
Jmat.Plot.hslToRgb = function(h, s, l) {
  var r, g, b;
  var temp1, temp2, tempr, tempg, tempb;
  h /= 256.0;
  s /= 256.0;
  l /= 256.0;
  if(s == 0) {
    r = g = b = l;
  } else {
    if(l < 0.5) temp2 = l * (1 + s);
    else temp2 = (l + s) - (l * s);
    temp1 = 2 * l - temp2;
    tempr = h + 1.0 / 3.0;
    if(tempr > 1) tempr--;
    tempg = h;
    tempb = h - 1.0 / 3.0;
    if(tempb < 0) tempb++;

    //Red
    if(tempr < 1.0 / 6.0) r = temp1 + (temp2 - temp1) * 6.0 * tempr;
    else if(tempr < 0.5) r = temp2;
    else if(tempr < 2.0 / 3.0) r = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempr) * 6.0;
    else r = temp1;
    //Green
    if(tempg < 1.0 / 6.0) g = temp1 + (temp2 - temp1) * 6.0 * tempg;
    else if(tempg < 0.5) g = temp2;
    else if(tempg < 2.0 / 3.0) g = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempg) * 6.0;
    else g = temp1;
    //Blue
    if(tempb < 1.0 / 6.0) b = temp1 + (temp2 - temp1) * 6.0 * tempb;
    else if(tempb < 0.5) b = temp2;
    else if(tempb < 2.0 / 3.0) b = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempb) * 6.0;
    else b = temp1;
  }
  return [Math.floor(r * 255), Math.floor(g * 255), Math.floor(b * 255)];
};

//input and output in range [0-255]
Jmat.Plot.hsvToRgb = function(h, s, v) {
  var r, g, b;
  h /= 256.0;
  s /= 256.0;
  v /= 256.0;
  if(s == 0) {
    r = g = b = v;
  } else {
    h *= 6;
    var i = Math.floor(h);
    var f = h - i;
    var p = v * (1 - s);
    var q = v * (1 - (s * f));
    var t = v * (1 - (s * (1 - f)));
    if(i == 0) { r = v; g = t; b = p; }
    if(i == 1) { r = q; g = v; b = p; }
    if(i == 2) { r = p; g = v; b = t; }
    if(i == 3) { r = p; g = q; b = v; }
    if(i == 4) { r = t; g = p; b = v; }
    if(i == 5) { r = v; g = p; b = q; }
  }
  return [Math.floor(r * 255), Math.floor(g * 255), Math.floor(b * 255)];
};

// rgb: [r, g, b], each in range 0-255
Jmat.Plot.rgbToCss = function(rgb) {
  var r = rgb[0].toString(16);
  var g = rgb[1].toString(16);
  var b = rgb[2].toString(16);

  if(r.length < 2) r = '0' + r;
  if(g.length < 2) g = '0' + g;
  if(b.length < 2) b = '0' + b;
  if(r.length > 2) r = 'ff';
  if(g.length > 2) g = 'ff';
  if(b.length > 2) b = 'ff';

  return '#' + r + g + b;
};

Jmat.Plot.complexColorFormula_ = 5; //0 = "tweaked wikipedia", 1 = "wikipedia", 2 = "original", 3 = "re/im", 4 = "re/im 2", 5=red/green/yellow/blue

Jmat.Plot.complexColorLS_ = function(y, maxval) {
  var s = 255;
  var l;
  var a = y.abs() / maxval;

  var m = 254; //max lightness for non-infinity (e.g. 240 or 250)

  // Original formula. Use with Jmat.Plot.hslToRgb
  if(Jmat.Plot.complexColorFormula_ == 2) {
    l = (m / 255) * a / (a + 1);
    s = 255;
  }

  // The Wikipedia formula. Use with hsvToRgb
  if(Jmat.Plot.complexColorFormula_ == 1) {
    l = 1 - 1 / (1.1 + 5*Math.log(a + 1));
    s = 255 / (1 + 0.3*Math.log(a + 1));
  }

  // Tweaked version of the Wikipedia formula. Use with hsvToRgb
  // Advantage over "original" formula: prettier, less visible "transition lines"
  // Disadvantage over "original" formula: less clear, less difference between magnitudes, value 1 is not perfectly #ff0000
  if(Jmat.Plot.complexColorFormula_ == 0) {
    var mm = 1 - m/255;
    l = (1 - 4*mm) - (1 - 8*mm) / (1 + 15*Math.log(a + 1));
    s = 255 / (1 + 0.3*Math.log(a + 1));
  }

  l *= 255;
  if(l < (255-m) && a > 0) l = (255-m);
  if(l > m) l = m;

  return [l, s];
};

// For "Color wheel graphs of complex functions" (not for "Domain coloring" which has repeating pattern of lightness rather than black for 0, white for infinity)
//maxvalue matches the halfway brightness of the HSL color model. Higher values go towards white, lower values go towards black, but are capped.
//hue is the argument. Positive real values are red, negative real values are cyan, positive imag values are grassgreen, negative imag values are purple, other colors are complex.
Jmat.Plot.getComplexColor1 = function(y, maxval) {
  var rgb;

  if(Jmat.Complex.isNaN(y)) {
    rgb = [128, 128, 128];
  } else {
    // Intention of the color (if maxval = 1):
    // +1 = pure red (255,0,0)
    // -1 = pure cyan (0,255,255)
    // +-Infinity = white (small mod applied to that below: set to lightness 254 instead of 255, to show complex argument of it)
    // 0 = black
    // hue = complex argument: red for 0 deg, yellow/green for 90 deg, cyan for 180 deg, purple for 270 deg.
    // |z| = 1: best possible saturation, lightness of hsv model == 128
    // |z| > 1: brigher than that
    // |z| < 1: darker than that

    var h = Jmat.Complex.arg1(y) * 255;

    var ls = Jmat.Plot.complexColorLS_(y, maxval);
    var l = ls[0];
    var s = ls[1];

    if(Jmat.Plot.complexColorFormula_ == 2) rgb = Jmat.Plot.hslToRgb(h, s, l);
    else rgb = Jmat.Plot.hsvToRgb(h, s, l);
  }

  return rgb;
};

// this one is based on re/im instead of abs/arg
// here, gray is 0, red indicates real, green indicates imaginary, with some b added that is average of r and g
Jmat.Plot.getComplexColor2 = function(y, maxval) {
  var rgb;

  if(Jmat.Complex.isNaN(y)) {
    rgb = [0, 0, 255];
  } else {
    var r = Math.abs(y.re / maxval / 2);
    var g = Math.abs(y.im / maxval / 2);
    var b = 0;

    var m = 254; //max lightness for non-infinity (e.g. 240 or 250)

    var mm = 1 - m/255;
    //r = (1 - 4*mm) - (1 - 8*mm) / (1 + 15*Math.log(r + 1));
    //r = 1 - 1 / (1.1 + 5*Math.log(r + 1));
    r = Math.log(r + 1);
    r = Math.min(Math.max(0, r), 1);
    //g = (1 - 4*mm) - (1 - 8*mm) / (1 + 15*Math.log(g + 1));
    //g = 1 - 1 / (1.1 + 5*Math.log(g + 1));
    g = Math.log(g + 1);
    g = Math.min(Math.max(0, g), 1);

    r = 128 + r * (y.re > 0 ? 1 : -1) * 127;
    g = 128 + g * (y.im > 0 ? 1 : -1) * 127;

    b = (r + g) / 2;

    rgb = [r, g, b];
  }

  return rgb;
};

// similar to getComplexColor2 but g and b swapped, this gives orange, white, blue and black as the hues.
Jmat.Plot.getComplexColor3 = function(y, maxval) {
  var rgb = Jmat.Plot.getComplexColor2(y, maxval);
  return [rgb[0], rgb[2], rgb[1]];
};


// This one has similarities to getComplexColor1, but different hues:
// positive = red
// negative = green (darker than "full bright" green to make it more distinguishable from yellow)
// pos imaginary = yellow (made a little bit darker, but not too much for this one since darker yellow looks too drab)
// neg imaginary = blue (made a bit brighter because otherwise it looks too dark compared to other colors of similar strength)
// the advantage: the hues are more memorable: two sets of "classical" opposing colors: red-green for the real axis, yellow-blue for the imaginary axis. And each time the "warm" color for positive (red and yellow) and the "cold" color for negative (green and blue)
// NOTE: sometimes, like in traffic lights "green" means positive and "red" means negative, but here we use the electrical convention where red wire is positive (and it also matches the standard hue pattern of getComplexColor1 better, less confusion)
Jmat.Plot.getComplexColor4 = function(y, maxval) {
  if(Jmat.Complex.isNaN(y)) {
    return [128, 128, 128];
  }
  if(Jmat.Complex.isInf(y)) {
    return [255, 255, 255];
  }

  var l = 0;

  var r = y.abs();


  // When it starts becoming more white from being very high value, depends on the exponent of the power here.
  // TODO: take maxval into account for this instead
  var l = r;
  l = 1 - Math.pow(2, -Math.pow(r, 0.33333333));

  // the eye is very sensitive to any discontinuity between different piecewise linear functions, for both hue and lightness.
  // Match the edges of a smoothstep function with the discontinuities to make it look better.
  // However, the cubicspline below can do even better.
  var smoothstep = function(x) { return x * x * (3 - 2 * x); };
  var smootherstep = function(x) { return x * x * x * (x * (x * 6 - 15) + 10); };
  var softsmooth = function(x) { return smoothstep(x) * 0.5 + x * 0.5; };
  if(l > 0.5) l = 0.5 + ((l - 0.5) * 0.98); // distinguish it from full white infinity, allow to slightly still see the hue
  l *= 255;

  var cubicspline = function(x, x0, y0, x1, y1, x2, y2) {
    var a12 = 1 / (x1 - x0);
    var a21 = a12;
    var a11 = 2 * a12;
    var a23 = 1 / (x2 - x1);
    var a32 = a23;
    var a33 = 2 * a23;
    var a22 = a11 + a33;
    var b1 = 3 * (y1 - y0) / ((x1 - x0) * (x1 - x0));
    var b3 = 3 * (y2 - y1) / ((x2 - x1) * (x2 - x1));
    var b2 = b1 + b3;

    var a = [[a11, a12, 0], [a21, a22, a23], [0, a32, a33]];
    var b = [b1, b2, b3];
    // TODO: perhaps this all can be significantly simplified
    var k = Jmat.Real.matrix_solve(a, b);
    var k0 = k[0];
    var k1 = k[1];
    var k2 = k[2];

    var a1 = k0 * (x1 - x0) - (y1 - y0);
    var b1 = -k1 * (x1 - x0) + (y1 - y0);
    var a2 = k1 * (x2 - x1) - (y2 - y1);
    var b2 = -k2 * (x2 - x1) + (y2 - y1);

    var t = x;
    if(t < x1) {
      t = (x - x0) / (x1 - x0);
      return (1 - t) * y0 + t * y1 + t * (1 - t) * ((1 - t) * a1 + t * b1);
    } else {
      t = (x - x1) / (x2 - x1);
      return (1 - t) * y1 + t * y2 + t * (1 - t) * ((1 - t) * a2 + t * b2);
    }
  };


  var a = Jmat.Complex.arg1(y);
  if(a >= 1) a -= 1; // avoid some potential artefacts in rendering

  var data = [
    [255, 0, 0], // pos (red)
    [224, 224, 0], // pos imag (yellow)
    [0, 192, 0], // neg (green)
    [64, 64, 255], // neg imag (blue)
  ];
  /*var data = [
    [255, 0, 0], // pos (red)
    [255, 160, 0],
    [255, 255, 0], // pos imag (yellow)
    [160, 224, 0],
    [0, 192, 0], // neg (green)
    [0, 128, 160],
    [64, 64, 255], // neg imag (blue)
    [128, 0, 192],
  ];*/
  a *= data.length;
  data.push(data[0]);
  var i = Math.floor(a);
  a -= i;
  if(data.length > 4) a = softsmooth(a);
  else a = smoothstep(a);
  var r0 = data[i][0];
  var g0 = data[i][1];
  var b0 = data[i][2];
  var r1 = data[i + 1][0];
  var g1 = data[i + 1][1];
  var b1 = data[i + 1][2];

  var r = r0 * (1 - a) + r1 * a;
  var g = g0 * (1 - a) + g1 * a;
  var b = b0 * (1 - a) + b1 * a;

  var applyLightness = function(rgb, l) {
    var result = [0, 0, 0];
    for(var c = 0; c < 3; c++) {
      var s = cubicspline(l, 0, 0, 128, rgb[c], 255, 255);
      result[c] = s;
    }
    return result;
  };
  var rgb = applyLightness([r, g, b], l);

  return rgb;
};

Jmat.Plot.getComplexColor_ = function(y, maxval) {
  // TODO: use nonnumeric names, the numeric codes for Jmat.Plot.complexColorFormula_ don't match the names
  if(Jmat.Plot.complexColorFormula_ <= 2) return Jmat.Plot.getComplexColor1(y, maxval);
  else if(Jmat.Plot.complexColorFormula_ == 3) return Jmat.Plot.getComplexColor2(y, maxval);
  else if(Jmat.Plot.complexColorFormula_ == 4) return Jmat.Plot.getComplexColor3(y, maxval);
  else return Jmat.Plot.getComplexColor4(y, maxval);
};

Jmat.Plot.getComplexColor = function(y, maxval) {
  var rgb = Jmat.Plot.getComplexColor_(y, maxval);
  return [Math.floor(rgb[0]), Math.floor(rgb[1]), Math.floor(rgb[2])];
};

Jmat.Plot.makeElement = function(parent, tag) {
  var el =  document.createElement(tag);
  parent.appendChild(el);
  return el;
};

Jmat.Plot.makeSizedElement = function(parent, tag, x, y, w, h) {
  var el =  document.createElement(tag);
  el.style.position = 'absolute';
  el.style.left = '' + Math.floor(x) + 'px';
  el.style.top = '' + Math.floor(y) + 'px';
  el.style.width = Math.floor(w) + 'px';
  el.style.height = Math.floor(h) + 'px';
  parent.appendChild(el);
  return el;
};

Jmat.Plot.makeSizedDiv = function(parent, x, y, w, h) {
  return Jmat.Plot.makeSizedElement(parent, 'div', x, y, w, h);
};

//'align' meaning: 0 = left/top, 1 = center, 2 = right/bottom
// if width is given (> 0), it is multiline text. Else it is single line text
// fontSize is a CSS value like 'small'
// rot: if true, the text is rotated 90 degrees counterclockwise. alignx/aligny then operate in respectively y/x direction on screen
Jmat.Plot.makeAlignedText = function(parent, text, width, x, y, alignx, aligny, fontSize, rot) {
  var div = document.createElement('div');
  if(fontSize) div.style.fontSize = fontSize;
  div.innerHTML = text;
  div.style.position = 'absolute';
  div.style.textAlign = alignx == 0 ? 'left' : alignx == 1 ? 'center' : 'right';
  if(!width) div.style.overflow = 'hidden';
  if(!width) div.style.whiteSpace = 'nowrap';

  if(width) div.style.width = width + 'px';
  parent.appendChild(div);
  var finalwidth = width;
  var h = div.clientHeight;
  var w = div.clientWidth;
  if(!width) finalwidth = w;

  if(aligny == 0) div.style.top = '' + Math.floor(y) + 'px';
  if(aligny == 1) div.style.top = '' + Math.floor(y - h / 2) + 'px';
  if(aligny == 2) div.style.top = '' + Math.floor(y - h) + 'px';

  if(alignx == 0) div.style.left = '' + Math.floor(x) + 'px';
  if(alignx == 1) div.style.left = '' + Math.floor(x - finalwidth / 2) + 'px';
  if(alignx == 2) div.style.left = '' + Math.floor(x - finalwidth) + 'px';

  if(rot) {
    // TODO: check if this is the correct rotation center in all cases
    var origx = (aligny == 0) ? 'top' : (aligny == 1 ? 'center' : 'bottom');
    var origy = (alignx == 0) ? 'left' : (alignx == 1 ? 'center' : 'right');
    div.style.transformOrigin = origx + ' ' + origy;
    div.style.transform = 'rotate(-90deg)';
  }

  return div;
};

//adds text vertically and horizontally centered, multiline depending on width.
//returns the text element
// if width is given (> 0), it is multiline text. Else it is single line text
Jmat.Plot.makeCenteredText = function(parent, text, width, x, y, fontSize, rot) {
  return Jmat.Plot.makeAlignedText(parent, text, width, x, y, 1, 1, fontSize, rot);
};

// w and h are the area for the plot itself, the canvas will be larger to add the labels
Jmat.Plot.initCanvas_ = function(parent, w, h, clickfun) {
  var canvas = document.createElement('canvas');
  canvas.style.position = 'absolute';
  canvas.style.left = 0;
  canvas.style.top = 0;
  canvas.width = parseInt(parent.style.width);
  canvas.height = parseInt(parent.style.height);
  var ctx = canvas.getContext("2d");
  var id = ctx.createImageData(w, h);
  var data  = id.data;
  parent.idd = id;
  parent.data = data;
  parent.ctx = ctx;
  parent.appendChild(canvas);
  if(clickfun) canvas.onclick = clickfun;
};

// Plot a pixel or rectangle to the given element. E.g. if w and h are 2, it's a 2x2 pixel.
// Position: x, y
// Size: w, h
// Color: [r, g, b], each in range 0-255
Jmat.Plot.rect = function(parent, x, y, w, h, rgb) {

  var id = parent.idd;
  var data = parent.data;
  var ctx = parent.ctx;

  // using fillRect goes faster than using 1-pixel data and "ctx.putImageData(id, x, y);"
  // TODO: however, using full size (or line by line) putImageData(id, x, y) may be even faster, consider that
  ctx.fillStyle = "rgba(" + rgb[0] + "," + rgb[1] + "," + rgb[2] + ",1)";
  ctx.fillRect(x, y, w, h);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Plotting - Internal implementation
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// TODO: render these on the canvas instead so they're part of the image when saving the canvas as image
Jmat.Plot.addPlotLabels_ = function(xlabel, ylabel, params, parent) {
  var width = params.w;
  var height = params.h;
  var plotx0 = 32;
  var plotx1 = plotx0 + width + 2;
  var ploty0 = 32;
  var ploty1 = ploty0 + height + 2;

  var xy0 = Jmat.Plot.fromPixel(params, 0, params.h);
  var xy1 = Jmat.Plot.fromPixel(params, params.w, 0);
  var xyc = Jmat.Plot.fromPixel(params, params.w >> 1, params.h >> 1);
  var x0 = xy0[0].toString(6);
  var y0 = xy0[1].toString(6);
  var x1 = xy1[0].toString(6);
  var y1 = xy1[1].toString(6);
  var xc = xyc[0].toString(6);
  var yc = xyc[1].toString(6);

  Jmat.Plot.makeSizedDiv(parent, plotx0 - 8, ploty0, 8, 1).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, plotx0 - 8, (ploty0 + ploty1) / 2, 8, 1).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, plotx0 - 8, ploty1 - 1, 8, 1).style.backgroundColor = 'black';
  Jmat.Plot.makeAlignedText(parent, '' + y1, 0, plotx0 - 8, ploty0,                    1, 2, 'small', true);
  Jmat.Plot.makeAlignedText(parent, '' + yc, 0, plotx0 - 8, (ploty0 + ploty1) / 2,     1, 2, 'small', true);
  Jmat.Plot.makeAlignedText(parent, '' + y0, 0, plotx0 - 8, ploty1,                    0, 2, 'small', true);
  Jmat.Plot.makeAlignedText(parent,  ylabel, 0, plotx0 - 8, (ploty0 + ploty1 * 3) / 4, 1, 2, 'small', true);


  Jmat.Plot.makeSizedDiv(parent, plotx0, ploty1, 1, 8).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, (plotx0 + plotx1) / 2, ploty1, 1, 8).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, plotx1 - 1, ploty1, 1, 8).style.backgroundColor = 'black';
  Jmat.Plot.makeAlignedText(parent, '' + x0, 0, plotx0 - 4,                ploty1 + 8, 0, 0, 'small');
  Jmat.Plot.makeAlignedText(parent, '' + xc, 0, (plotx0 + plotx1) / 2,     ploty1 + 8, 1, 0, 'small');
  Jmat.Plot.makeAlignedText(parent, '' + x1, 0, plotx1,                    ploty1 + 8, 1, 0, 'small');
  Jmat.Plot.makeAlignedText(parent,  xlabel, 0, (plotx0 * 3 + plotx1) / 4, ploty1 + 8, 1, 0, 'small');
};

//type: 0=real, 1=complex, 2=2D
//returns the element in which result of clicking can be displayed
Jmat.Plot.addControls_ = function(params, parent, plotfun, type) {
  var L = 32;
  var d;

  /*var width = params.w;
  var height = params.h;
  var size = params.xsize / 2;*/

  var color = '#111';

  var s = 24; // button size
  var s2 = s + 6; // button size + separation

  var makeButton = function(parent, label, tooltip, x, y, fun) {
    var div = Jmat.Plot.makeSizedDiv(parent, x, y, s, s);
    div.style.border = '1px solid #aaa';
    var text = Jmat.Plot.makeCenteredText(div, label, 0, s >> 1, s >> 1);
    text.style.color = '#111';
    div.onclick = fun;
    div.title = tooltip;
  };

  var makeField = function(parent, label, tooltip, x, y, value, fun) {
    //var div = Jmat.Plot.makeSizedDiv(parent, x, y, s * 6, s);
    var el = Jmat.Plot.makeSizedElement(parent, 'input', x, y, s * 6, s - 2);
    el.style.border = '1px solid #aaa';
    el.style.padding = '1px';
    el.value = value;
    el.onchange = function() {
      fun(el.value);
    };
    el.title = tooltip;
    return el;
  };

  var xstart = L;
  var ystart = params.h + L + 32;

  var updatefun = function() {
    Jmat.Plot.stopPlotting();
    xfield.value = Jmat.Complex(params.xshift, params.xshift_im).toString();
    yfield.value = Jmat.Complex(params.yshift, params.yshift_im).toString();
    xsizefield.value = params.xsize;
    ysizefield.value = params.ysize;
    plotfun();
  };

  makeButton(parent, '←', 'shift left', xstart + s2 * 2, ystart + s2 * 0, function() {
    params.xshift -= params.xsize / 5; updatefun();
  });
  makeButton(parent, '→', 'shift right', xstart + s2 * 3, ystart + s2 * 0, function() {
    params.xshift += params.xsize / 5; updatefun();
  });
  var xfield = makeField(parent, 'x', 'center horizontal coordinate. Supports complex values for 2D plots.', xstart + s2 * 4, ystart + s2 * 0, Jmat.Complex(params.xshift, params.xshift_im).toString(), function(value) {
    var v = Jmat.Complex(value);
    params.xshift = v.re;
    params.xshift_im = v.im;
    updatefun();
  });

  makeButton(parent, '↑', 'shift up', xstart + s2 * 2, ystart + s2 * 1, function() {
    params.yshift += params.ysize / 5; updatefun();
  });
  makeButton(parent, '↓', 'shift down', xstart + s2 * 3, ystart + s2 * 1, function() {
    params.yshift -= params.ysize / 5; updatefun();
  });
  var yfield = makeField(parent, 'y', 'center vertical coordinate. Supports complex values for 2D plots.', xstart + s2 * 4, ystart + s2 * 1, Jmat.Complex(params.yshift, params.yshift_im).toString(), function(value) {
    var v = Jmat.Complex(value);
    params.yshift = v.re;
    params.yshift_im = v.im;
    updatefun();
  });

  makeButton(parent, '+', 'zoom in', xstart + s2 * 10, ystart + s2 * 0, function() {
    params.xsize /= 2; params.ysize /= 2; updatefun();
  });
  var xsizefield = makeField(parent, 'x', 'area size in horizontal direction', xstart + s2 * 11, ystart + s2 * 0, params.xsize, function(value) {
    var v = Jmat.Complex(value);
    params.xsize = v.re;
    updatefun();
  });
  makeButton(parent, '-', 'zoom out', xstart + s2 * 10, ystart + s2 * 1, function() {
    params.xsize *= 2; params.ysize *= 2; updatefun();
  });
  var ysizefield = makeField(parent, 'y', 'area size in vertical direction', xstart + s2 * 11, ystart + s2 * 1, params.ysize, function(value) {
    var v = Jmat.Complex(value);
    params.ysize = v.re;
    updatefun();
  });

  var origparams = new Jmat.Plot.Params(params);
  makeButton(parent, 'R', 'reset', xstart + s2 * 0, ystart + s2 * 0, function() {
    Jmat.Plot.copyParams_(origparams, params);
    updatefun();
  });

  makeButton(parent, 't', 'twiddle pixel integer position when zoomed out', xstart + s2 * 0, ystart + s2 * 1, function() {
    params.twiddle = !params.twiddle;
    updatefun();
  });

  d = Jmat.Plot.makeCenteredText(parent, '        ', 0, L + 16, L - 16);
  d.style.fontSize = 'small';
  return d;
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Plot.makeRealPixel_ = function(div, params, px, y, prevy, rgb) {
  var p = params.p;
  var ysize = params.ysize;
  var width = params.w;
  var height = params.h;

  var py = Math.floor(height / 2 - ((y+params.yshift) / ysize * height) - 1);
  var prevpy = Math.floor(height / 2 - ((prevy+params.yshift) / ysize * height) - 1);

  if(py < 0 && prevpy < 0) return;
  if(py > height && prevpy > height) return;

  if(isNaN(y)) {
    var d = Jmat.Plot.rect(div, px, 0, p, height, [160,160,160]);
    return;
  }

  if(py >= 0 && py < height) {
    var d = Jmat.Plot.rect(div, px, py, p, p, rgb);
  }

  if(py < 0) py = 0;
  if(py > height) py = height;
  if(prevpy < 0) prevpy = 0;
  if(prevpy > height) prevpy = height;

  if(prevpy + p < py) {
    var d = Jmat.Plot.rect(div, px, prevpy + p, p, (py - prevpy), rgb);
  } else if(prevpy > py + p) {
    var d = Jmat.Plot.rect(div, px, py, p, (prevpy - py), rgb);
  }
};

Jmat.Plot.plotReal_ = function(fun, params, parent, label) {
  var width = params.w;
  var height = params.h;
  var p = params.p;

  var L = 32;
  var steps = Math.floor(width / p);



  parent.innerHTML = '';
  var div = Jmat.Plot.makeSizedDiv(parent, L, L, width, height);
  var labeldiv = Jmat.Plot.makeSizedDiv(parent, 0, 0, 1, 1);
  //parent.style.backgroundColor = 'white'; // TODO: make element inside it instead to alter style
  var plotfun = function() {
    var prevy;

    div.innerHTML = '';

    labeldiv.innerHTML = '';
    Jmat.Plot.addPlotLabels_('x', 'y', params, labeldiv);
    if(label) Jmat.Plot.makeAlignedText(labeldiv, label, 0, L + width, L, 2, 2);

    // axis lines
    var d = Jmat.Plot.makeSizedDiv(div, 0, width / 2, height, 2);
    d.style.backgroundColor = '#ccc';
    d = Jmat.Plot.makeSizedDiv(div, width / 2, 0, 2, height);
    d.style.backgroundColor = '#ccc';

    Jmat.Plot.initCanvas_(div, width, height, function(e) {
      var x = e.offsetX;
      var y = e.offsetY;
      var xy = Jmat.Plot.fromPixel(params, x, y);
      var value = fun(xy[0]);
      var label = '' + Jmat.toString(xy[0], 4) + ', ' + Jmat.toString(xy[1], 4) + ': ' + value.toString(6);
      el.innerText = label;
    });

    var xsize = params.xsize;
    var ysize = params.ysize;


    for(var i = 0; i < steps; i++) {
      var px = i * p;
      var re = -xsize / 2 + (i / steps * xsize) + params.xshift;

      var x = Jmat.Complex(re, 0);
      var y = fun(x);
      if(!prevy) prevy = y;
      if(y.re == Infinity && y.im == Infinity) y = Jmat.Complex(NaN); // plot undirected infinity as NaN (should show up as vertical line)

      if(y.im == 0 || (Math.abs(y.im) < Math.abs(y.re) * 1e-10 /*imag likely due to numerical imprecisions*/)) {
        Jmat.Plot.makeRealPixel_(div, params, px, y.re, prevy.re, [0,0,0], Jmat.Complex.toString(x) + ': ' + Jmat.Complex.toString(y));
      } else {
        // Abs and arg-color-wheel, always in positive zone
        //var h = Jmat.Complex.arg1(y) * 255;
        //var rgb = Jmat.Plot.hslToRgb(h, 255, 128);
        var temp = Jmat.Complex.polar(1, y.arg());
        var rgb = Jmat.Plot.getComplexColor(temp, 1);
        var a = y.abs();
        var pa = prevy.abs();
        Jmat.Plot.makeRealPixel_(div, params, px, a, pa, rgb, Jmat.Complex.toString(x) + ': ' + Jmat.Complex.toString(y));
      }

      prevy = y;
    }
  };

  div.style.backgroundColor = '#eee';
  div.style.border = '1px solid black';

  var el = Jmat.Plot.addControls_(params, parent, plotfun, 0);


  plotfun();
};

//p = pixel cell size
// For "Color wheel graphs of complex functions" (not for "Domain coloring" which has repeating pattern of lightness rather than black for 0, white for infinity)
Jmat.Plot.plotColorPixel = function(y, maxval, p, px, py, div) {
  var rgb = Jmat.Plot.getComplexColor(y, maxval);
  Jmat.Plot.rect(div, px, py, p, p, rgb);
};

////////////////////////////////////////////////////////////////////////////////

//p = pixel cell size
Jmat.Plot.plot2DPixel_ = function(fun, params, px, py, div) {
  /*var x = -params.xsize * 0.5 + (px / steps * params.xsize);
  var y = params.ysize * 0.5 - (py / steps * params.ysize);

  var sx = Jmat.Complex(x + params.xshift, params.xshift_im);
  var sy = Jmat.Complex(y + params.yshift, params.yshift_im);*/

  var xy = Jmat.Plot.fromPixel(params, px, py);
  var sx = xy[0];
  var sy = xy[1];

  var z = fun(sx, sy);

  Jmat.Plot.plotColorPixel(z, params.v, params.p, px, py, div);
};


Jmat.Plot.plot2DLineTimeout_ = function(fun, params, py, div) {
  var stopindex = Jmat.Plot.stopIndex_;

  // This is for first rendering fast, and only then at full resolution
  var stage = 1;
  var params1 = params;

  if(params.p <= 2) {
    stage = 2;
    params = new Jmat.Plot.Params(params);
    params.p *= 4;
  }

  var linefun = function(py) {
    window.setTimeout(function() {
      for(var i = 0; i < 4; i++) {
        if(py >= params.h) {
          stage--;
          if(stage <= 0) return;
          params = params1;
          py = 0;
        }
        for(var px = 0; px < params.w; px += params.p) {
          Jmat.Plot.plot2DPixel_(fun, params, px, py, div);
        }
        py += params.p;
      }
      if(stopindex != Jmat.Plot.stopIndex_) return;
      linefun(py);
    }, 0);
  };
  linefun(py);
};

Jmat.Plot.plot2DNonBlocking_ = function(fun, params, div) {
  Jmat.Plot.plot2DLineTimeout_(fun, params, 0, div);
};

//type: 1=complex, 2=2d
Jmat.Plot.plot2D_ = function(fun, params, parent, label, xlabel, ylabel, type) {
  //parent.style.backgroundColor = 'white'; // TODO: make element inside it instead to alter style
  if(!xlabel) xlabel = 'x';
  if(!ylabel) ylabel = 'y';


  var width = params.w;
  var height = params.h;
  var size = params.xsize / 2;
  var p = params.p;
  var maxval = params.v;

  var L = 32;
  var steps = Math.floor(width / p);
  parent.innerHTML = '';
  var div = Jmat.Plot.makeSizedDiv(parent, L, L, steps * p, steps * p);

  var labeldiv = Jmat.Plot.makeSizedDiv(parent, 0, 0, 1, 1);

  var plotfun = function() {
    labeldiv.innerHTML = '';
    Jmat.Plot.addPlotLabels_(xlabel, ylabel, params, labeldiv);
    if(label) Jmat.Plot.makeAlignedText(labeldiv, label, 0, L + width, L, 2, 2);
    Jmat.Plot.plot2DNonBlocking_(fun, params, div);
  };


  div.style.backgroundColor = '#888888';
  div.style.border = '1px solid black';

  var el = Jmat.Plot.addControls_(params, parent, plotfun, type);

  Jmat.Plot.initCanvas_(div, width, height, function(e) {
    var x = e.offsetX;
    var y = e.offsetY;
    var xy = Jmat.Plot.fromPixel(params, x, y);
    var value = fun(xy[0], xy[1]);
    var label = '' + Jmat.toString(xy[0], 4) + ', ' + Jmat.toString(xy[1], 4) + ': ' + value.toString(6);
    el.innerText = label;
  });

  plotfun();
};

