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

// NOTE: requires charset="utf-8" to render navigation arrow characters correctly

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Plotting - public API
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// TODO: support any size rather than always 320x320 pixels

// e.g.: new Jmat.PlotParams({p: 1})
// constructor
Jmat.PlotParams = function(o) {
  if(!o) o = {};

  this.p = o.p != undefined ? o.p : 2; //pixel size (==> resolution). 1 = highest resolution, 2 = half resolution, etc...

  this.xsize = o.xsize != undefined ? o.xsize : 20; // e.g. an xsize of 10 makes it go from x=-5 to +5 if shift is 0
  this.ysize = o.ysize != undefined ? o.ysize : 20;
  if(o.s && !o.xsize && !o.ysize) this.xsize = this.ysize = o.s; // 's' is shortcut to make both that size
  this.xshift = o.xshift != undefined ? o.xshift : 0; // point which you want in the center of the plot
  this.yshift = o.yshift != undefined ? o.yshift : 0;

  // these are only used for 2D plot
  this.xshift_im = 0;
  this.yshift_im = 0;
  this.transpose = false;

  // The value at which the complex color wheel has highest saturation (pure red for positive real). Higher value gives more white color, lower gives darker color.
  // Not used by real plot.
  this.v = o.v != undefined ? o.v : 1;
};

//fun = mathematical function taking 1 Jmat.Complex arguments, e.g. Jmat.sin
//params = parameter object with plot size, resolution, etc.... See Jmat.PlotParams.
//parent = HTML parent element, best of type div. All necessary elements (e.g. canvas) will be created inside of it.
Jmat.plotReal = function(fun, parent, params, label) {
  if(!params) params = new Jmat.PlotParams();
  if(!(params instanceof Jmat.PlotParams)) {
    params = new Jmat.PlotParams(params);
  }
  if(!parent) parent = document.body;
  Jmat.Plot.plotReal_(fun, params, parent, label);
};

//fun = mathematical function taking 1 Jmat.Complex arguments, e.g. Jmat.gamma
//params = parameter object with plot size, resolution, etc.... See Jmat.PlotParams.
//parent = HTML parent element, best of type div. All necessary elements (e.g. canvas) will be created inside of it.
Jmat.plotComplex = function(fun, parent, params, label) {
  if(!params) params = new Jmat.PlotParams({p: 4});
  if(!(params instanceof Jmat.PlotParams)) {
    params = new Jmat.PlotParams(params);
  }
  if(!parent) parent = document.body;
  Jmat.Plot.plot2D_(function(x, y) {
    return fun(Jmat.Complex(x.re, y.re));
  }, params, parent, label, 're', 'im');
};

//fun = mathematical function taking 2 Jmat.Complex arguments, e.g. Jmat.besselj
//params = parameter object with plot size, resolution, etc.... See Jmat.PlotParams.
//parent = HTML parent element, best of type div. All necessary elements (e.g. canvas) will be created inside of it.
Jmat.plot2D = function(fun, parent, params, label) {
  if(!params) params = new Jmat.PlotParams();
  if(!(params instanceof Jmat.PlotParams)) {
    params = new Jmat.PlotParams(params);
  }
  if(!parent) parent = document.body;
  Jmat.Plot.plot2D_(fun, params, parent, label);
};

// Stop plot2D or plotComplex, if they are taking a very long time. May still
// calculate several pixels before actually stopping: they do several lines at the time.
Jmat.stopPlotting = function() {
  Jmat.Plot.stopIndex_ = (Jmat.Plot.stopIndex_ ? Jmat.Plot.stopIndex_ + 1 : 1);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Graphics (color & DOM) helper functions
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// Jmat graphics, graphing and plotting library
// internal functions are grouped in here
Jmat.Plot = function() {
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

Jmat.Plot.complexColorFormula_ = 0; //0 = "tweaked wikipedia", 1 = "wikipedia", 2 = "original"

// For "Color wheel graphs of complex functions" (not for "Domain coloring" which has repeating pattern of lightness rather than black for 0, white for infinity)
//maxvalue matches the halfway brightness of the HSL color model. Higher values go towards white, lower values go towards black, but are capped.
//hue is the argument. Positive real values are red, negative real values are cyan, positive imag values are grassgreen, negative imag values are purple, other colors are complex.
Jmat.Plot.getComplexColor = function(y, maxval) {
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
    if(l < (255-m) & a > 0) l = (255-m);
    if(l > m) l = m;

    if(Jmat.Plot.complexColorFormula_ == 2) rgb = Jmat.Plot.hslToRgb(h, s, l);
    else rgb = Jmat.Plot.hsvToRgb(h, s, l);
  }

  return rgb;
};

Jmat.Plot.makeSizedDiv = function(parent, x, y, w, h) {
  var el =  document.createElement('div');
  el.style.position = 'absolute';
  el.style.left = '' + Math.floor(x) + 'px';
  el.style.top = '' + Math.floor(y) + 'px';
  el.style.width = Math.floor(w) + 'px';
  el.style.height = Math.floor(h) + 'px';
  parent.appendChild(el);
  return el;
};

//'align' meaning: 0 = left/top, 1 = center, 2 = right/bottom
// if width is given (> 0), it is multiline text. Else it is single line text
// fontSize is a CSS value like 'small'
Jmat.Plot.makeAlignedText = function(parent, text, width, x, y, alignx, aligny, fontSize) {
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

  return div;
};

//adds text vertically and horizontally centered, multiline depending on width.
//returns the text element
// if width is given (> 0), it is multiline text. Else it is single line text
Jmat.Plot.makeCenteredText = function(parent, text, width, x, y, fontSize) {
  return Jmat.Plot.makeAlignedText(parent, text, width, x, y, 1, 1, fontSize)
};

Jmat.Plot.useHTML5canvas_ = true;

// Plot a pixel or rectangle to the given element. E.g. if w and h are 2, it's a 2x2 pixel.
// Position: x, y
// Size: w, h
// Color: [r, g, b], each in range 0-255
Jmat.Plot.rect = function(parent, x, y, w, h, rgb, label) {
  // NON-canvas version (slower)
  if(!Jmat.Plot.useHTML5canvas_) {
    var el =  Jmat.Plot.makeSizedDiv(parent, x, y, w, h);
    el.style.backgroundColor = Jmat.Plot.rgbToCss(rgb);
    el.title = label;
    return el;
  }

  // canvas version (faster)
  var data = parent.data;
  var id = parent.idd;
  var ctx = parent.ctx;

  if(!parent.labeldata) parent.labeldata = [];
  for(var y2 = y; y2 < y + h; y2++) {
    if(!parent.labeldata[y2]) parent.labeldata[y2] = [];
    for(var x2 = x; x2 < x + w; x2++) {
      parent.labeldata[y2][x2] = label;
    }
  }

  if(!data) {
    var canvas =  document.createElement('canvas');
    canvas.style.position = 'absolute';
    canvas.style.left = 0;
    canvas.style.top = 0;
    canvas.width = parseInt(parent.style.width);
    canvas.height = parseInt(parent.style.height);
    ctx = canvas.getContext("2d");
    id = ctx.createImageData(1,1);
    data  = id.data;
    parent.idd = id;
    parent.data = data;
    parent.ctx = ctx;
    parent.appendChild(canvas);

    canvas.onclick = function(e) {
      var x = e.offsetX;
      var y = e.offsetY;
      // TODO: somehow display this in a tooltip on the canvas instead
      console.log(parent.labeldata[y][x]);
    };
  }

  if(w == 1 && h == 1) {
    data[0] = rgb[0];
    data[1] = rgb[1];
    data[2] = rgb[2];
    data[3] = 255;
    ctx.putImageData(id, x, y);
  } else {
    ctx.fillStyle = "rgba(" + rgb[0] + "," + rgb[1] + "," + rgb[2] + ",1)";
    ctx.fillRect(x, y, w, h);
  }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Plotting - Internal implementation
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Jmat.Plot.addPlotLabels_ = function(xlabel, ylabel, x0, x1, y0, y1, parent) {
  var plotx0 = 32;
  var plotx1 = plotx0 + 322;
  var ploty0 = 32;
  var ploty1 = ploty0 + 322;

  var yc = (y0 + y1) / 2;
  Jmat.Plot.makeSizedDiv(parent, plotx0 - 8, ploty0, 8, 1).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, plotx0 - 8, (ploty0 + ploty1) / 2, 8, 1).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, plotx0 - 8, ploty1 - 1, 8, 1).style.backgroundColor = 'black';
  Jmat.Plot.makeAlignedText(parent, '' + y1, 0, plotx0 - 8 - 2, ploty0, 2, 1, 'small');
  Jmat.Plot.makeAlignedText(parent, '' + yc, 0, plotx0 - 8 - 2, (ploty0 + ploty1) / 2, 2, 1, 'small');
  Jmat.Plot.makeAlignedText(parent, '' + y0, 0, plotx0 - 8 - 2, ploty1, 2, 1, 'small');
  Jmat.Plot.makeAlignedText(parent, ylabel, 0, plotx0 - 2, (ploty0 + ploty1 * 3) / 4, 2, 1, 'small');

  var xc = (x0 + x1) / 2;
  Jmat.Plot.makeSizedDiv(parent, plotx0, ploty1, 1, 8).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, (plotx0 + plotx1) / 2, ploty1, 1, 8).style.backgroundColor = 'black';
  Jmat.Plot.makeSizedDiv(parent, plotx1 - 1, ploty1, 1, 8).style.backgroundColor = 'black';
  Jmat.Plot.makeAlignedText(parent, '' + x0, 0, plotx0, ploty1 + 8, 1, 0, 'small');
  Jmat.Plot.makeAlignedText(parent, '' + xc, 0, (plotx0 + plotx1) / 2, ploty1 + 8, 1, 0, 'small');
  Jmat.Plot.makeAlignedText(parent, '' + x1, 0, plotx1, ploty1 + 8, 1, 0, 'small');
  Jmat.Plot.makeAlignedText(parent, xlabel, 0, (plotx0 * 3 + plotx1) / 4, ploty1, 1, 0, 'small');
};

////////////////////////////////////////////////////////////////////////////////

Jmat.Plot.makeRealPixel_ = function(div, width, height, params, px, y, prevy, rgb, label) {
  var p = params.p;
  var ysize = params.ysize;

  var py = Math.floor(height / 2 - ((y+params.yshift) / ysize * height) - 1);
  var prevpy = Math.floor(height / 2 - ((prevy+params.yshift) / ysize * height) - 1);

  if(py < 0 && prevpy < 0) return;
  if(py > height && prevpy > height) return;

  if(isNaN(y)) {
    var d = Jmat.Plot.rect(div, px, 0, p, height, [160,160,160], label);
    return;
  }

  if(py >= 0 && py < height) {
    var d = Jmat.Plot.rect(div, px, py, p, p, rgb, label);
  }

  if(py < 0) py = 0;
  if(py > height) py = height;
  if(prevpy < 0) prevpy = 0;
  if(prevpy > height) prevpy = height;

  if(prevpy + p < py) {
    var d = Jmat.Plot.rect(div, px, prevpy + p, p, (py - prevpy), rgb, label);
  } else if(prevpy > py + p) {
    var d = Jmat.Plot.rect(div, px, py, p, (prevpy - py), rgb, label);
  }
};

Jmat.Plot.plotReal_ = function(fun, params, parent, label) {
  //parent.style.backgroundColor = 'white'; // TODO: make element inside it instead to alter style
  var plotfun = function() {
    var width = 320;
    var height = 320;
    var xsize = params.xsize;
    var ysize = params.ysize;
    var p = params.p;

    var L = 32;
    var steps = Math.floor(width / p);
    parent.innerHTML = '';
    var div = Jmat.Plot.makeSizedDiv(parent, L, L, width, height);
    div.style.backgroundColor = '#eee';
    div.style.border = '1px solid black';

    // axis lines
    var d = Jmat.Plot.makeSizedDiv(div, 0, width / 2, height, 2);
    d.style.backgroundColor = '#ccc';
    d = Jmat.Plot.makeSizedDiv(div, width / 2, 0, 2, height);
    d.style.backgroundColor = '#ccc';

    Jmat.Plot.addPlotLabels_('x', 'y', -xsize/2+params.xshift, xsize/2+params.xshift, -ysize/2-params.yshift, ysize/2-params.yshift, parent);
    if(label) Jmat.Plot.makeAlignedText(parent, label, 0, L + width, L, 2, 2);

    d = Jmat.Plot.makeCenteredText(parent, '←', 0, L + width / 2 - 35, L - 10);
    d.onclick = function() { params.xshift -= params.xsize / 8; plotfun(); };
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, '→', 0, L + width / 2 - 15, L - 10);
    d.onclick = function() { params.xshift += params.xsize / 8; plotfun(); };
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, '[+]', 0, L + width / 2 + 15, L - 10);
    d.onclick = function() { params.xsize /= 2; params.ysize /= 2; plotfun(); };
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, '[-]', 0, L + width / 2 + 35, L - 10);
    d.onclick = function() { params.xsize *= 2; params.ysize *= 2; plotfun(); };
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, '+', 0, L + width + 8, L + 70);
    d.onclick = function() { params.ysize /= 2; plotfun(); };
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, '-', 0, L + width + 8, L + 90);
    d.onclick = function() { params.ysize *= 2; plotfun(); };
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, '↑', 0, L + width + 8, L + 50);
    d.onclick = function() { params.yshift -= params.ysize / 8; plotfun(); };
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, '↓', 0, L + width + 8, L + 110);
    d.onclick = function() { params.yshift += params.ysize / 8; plotfun(); };
    d.style.color = '#ddd';

    var prevy;

    for(var i = 0; i < steps; i++) {
      var px = i * p;
      var re = -xsize / 2 + (i / steps * xsize) + params.xshift;

      var x = Jmat.Complex(re, 0);
      var y = fun(x);
      if(!prevy) prevy = y;
      if(y.re == Infinity && y.im == Infinity) y = Jmat.Complex(NaN); // plot undirected infinity as NaN (should show up as vertical line)

      if(y.im == 0 || (Math.abs(y.im) < Math.abs(y.re) * 1e-10 /*imag likely due to numerical imprecisions*/)) {
        Jmat.Plot.makeRealPixel_(div, 320, 320, params, px, y.re, prevy.re, [0,0,0], Jmat.Complex.toString(x) + ': ' + Jmat.Complex.toString(y));
      } else {
        // Abs and arg-color-wheel, always in positive zone
        var h = Jmat.Complex.arg1(y) * 255;
        var rgb = Jmat.Plot.hslToRgb(h, 255, 128);
        var a = y.abs();
        var pa = prevy.abs();
        Jmat.Plot.makeRealPixel_(div, 320, 320, params, px, a, pa, rgb, Jmat.Complex.toString(x) + ': ' + Jmat.Complex.toString(y));
      }

      prevy = y;
    }
  };

  plotfun();
};

//p = pixel cell size
// For "Color wheel graphs of complex functions" (not for "Domain coloring" which has repeating pattern of lightness rather than black for 0, white for infinity)
Jmat.Plot.plotColorPixel = function(y, maxval, p, px, py, div, label) {
  var rgb = Jmat.Plot.getComplexColor(y, maxval);
  var d = Jmat.Plot.rect(div, px * p, py * p, p, p, rgb, label);
  return d;
};

////////////////////////////////////////////////////////////////////////////////

//p = pixel cell size
Jmat.Plot.plot2DPixel_ = function(fun, size, steps, params, px, py, div, label) {
  var x = -size + (px / steps * size * 2);
  var y = size - (py / steps * size * 2);

  var sx = Jmat.Complex(x + params.xshift, params.xshift_im);
  var sy = Jmat.Complex(y + params.yshift, params.yshift_im);

  var z;
  if(params.transpose) {
    z = fun(sy, sx);
  } else {
    z = fun(sx, sy);
  }

  var label = sx + ', ' + sy + ': ' + Jmat.Complex.toString(z);
  var d = Jmat.Plot.plotColorPixel(z, params.v, params.p, px, py, div, label);
};


Jmat.Plot.plot2DLineTimeout_ = function(fun, size, steps, params, py, div) {
  var stopindex = Jmat.Plot.stopIndex_;

  // This is for first rendering fast, and only then at full resolution
  var stage = 1;
  var params1 = params;
  var steps1 = steps;

  if(params.p <= 2) {
    stage = 0;
    params = JSON.parse(JSON.stringify(params));
    params.p *= 4;
    steps /= 4;
  }

  var linefun = function(py) {
    window.setTimeout(function() {
      for(var i = 0; i < 4; i++) {
        if(py == steps) { stage++; params = params1; steps = steps1; py = 0; }
        if(stage == 2) return;
        if(py == steps) return;
        if(stage == 2) return;
        for(var px = 0; px < steps; px++) {
          Jmat.Plot.plot2DPixel_(fun, size, steps, params, px, py, div);
        }
        py++;
      }
      if(stopindex != Jmat.Plot.stopIndex_) return;
      linefun(py);
    }, 0);
  };
  linefun(py);
};

Jmat.Plot.plot2DNonBlocking_ = function(fun, size, steps, params, div) {
  Jmat.Plot.plot2DLineTimeout_(fun, size, steps, params, 0, div);
};

Jmat.Plot.plot2D_ = function(fun, params, parent, label, xlabel, ylabel) {
  //parent.style.backgroundColor = 'white'; // TODO: make element inside it instead to alter style
  if(!xlabel) xlabel = 'x';
  if(!ylabel) ylabel = 'y';
  var plotfun = function() {
    var width = 320;
    var height = 320;
    var size = params.xsize / 2;
    // TODO: support x and y size
    //var xsize = params.xsize;
    //var ysize = params.ysize;
    var p = params.p;
    var maxval = params.v;

    var L = 32;
    var steps = Math.floor(width / p);
    parent.innerHTML = '';
    var div = Jmat.Plot.makeSizedDiv(parent, L, L, steps * p, steps * p);
    div.style.backgroundColor = '#888888';
    div.style.border = '1px solid black';

    Jmat.Plot.addPlotLabels_(xlabel, ylabel, -size+params.xshift, size+params.xshift, -size+params.yshift, size+params.yshift, parent);
    if(label) Jmat.Plot.makeAlignedText(parent, label, 0, L + width, L, 2, 2);

    d = Jmat.Plot.makeCenteredText(parent, '←', 0, L + width / 2 - 35, L - 10);
    d.onclick = function() { params.xshift -= params.xsize / 5; Jmat.stopPlotting(); plotfun(); };
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, '→', 0, L + width / 2 - 15, L - 10);
    d.onclick = function() { params.xshift += params.xsize / 5; Jmat.stopPlotting(); plotfun(); };
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, '[+]', 0, L + width / 2 + 15, L - 10);
    d.onclick = function() { params.xsize /= 2; params.ysize /= 2; Jmat.stopPlotting(); plotfun(); };
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, '[-]', 0, L + width / 2 + 35, L - 10);
    d.onclick = function() { params.xsize *= 2; params.ysize *= 2; Jmat.stopPlotting(); plotfun(); };
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, '↑', 0, L + width + 8, L + 70);
    d.onclick = function() { params.yshift += params.ysize / 5; Jmat.stopPlotting(); plotfun(); };
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, '↓', 0, L + width + 8, L + 90);
    d.onclick = function() { params.yshift -= params.ysize / 5; Jmat.stopPlotting(); plotfun(); };
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, 'x_im+', 0, L + width + 30, L + 110);
    d.onclick = function() { params.xshift_im += 0.1; Jmat.stopPlotting(); plotfun(); };
    d.title = params.xshift_im;
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, 'x_im-', 0, L + width + 30, L + 125);
    d.onclick = function() { params.xshift_im -= 0.1; Jmat.stopPlotting(); plotfun(); };
    d.title = params.xshift_im;
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, 'y_im+', 0, L + width + 30, L + 140);
    d.onclick = function() { params.yshift_im += 0.1; Jmat.stopPlotting(); plotfun(); };
    d.title = params.yshift_im;
    d.style.color = '#ddd';
    d = Jmat.Plot.makeCenteredText(parent, 'y_im-', 0, L + width + 30, L + 155);
    d.onclick = function() { params.yshift_im -= 0.1; Jmat.stopPlotting(); plotfun(); };
    d.title = params.yshift_im;
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, 't', 0, L + width + 30, L + 170);
    d.onclick = function() { params.transpose = !params.transpose; Jmat.stopPlotting(); plotfun(); };
    d.style.color = '#ddd';

    d = Jmat.Plot.makeCenteredText(parent, 'r', 0, L + width / 2 - 70, L - 10);
    d.onclick = function() {
      // This button shifts by 1 half. Some functions have different value, or algorithm, for integers (e.g. negative bessel J). This reveals it when zoomed out.
      var s = params.xshift;
      params.xshift = Math.floor(params.xshift);
      params.yshift = Math.floor(params.yshift);
      if(s == Math.floor(s)) {
        params.xshift += 0.5;
        params.yshift += 0.5;
      };
      Jmat.stopPlotting();
      plotfun();
    }
    d.style.color = '#ddd';

    Jmat.Plot.plot2DNonBlocking_(fun, size, steps, params, div);
  }
  plotfun();
};

