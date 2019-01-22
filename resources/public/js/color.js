// translate a single component of an RGB triple from 0..1 to #00...ff
function componentToHex(n) {
  var c = Math.floor(n * 255);
  var hex = c.toString(16);
  return hex.length == 1 ? "0" + hex : hex;
}

// translate an RGB triple of the form [0..1, 0..1, 0..1] into a hex value #000000...ffffff
function rgbToHex(rgb) {
  return "#" + componentToHex(rgb[0]) + componentToHex(rgb[1]) + componentToHex(rgb[2]);
}

// build a function that given a value between 0 and 1, returns an RGB triple scaled between
// the supplied bottom and top RGB triples (0 is full bottom, 1 is full top).
function colorScale(bottom, top) {
  return function(value) {
    return _.map(_.zip(bottom, top), function(spread) {
      return spread[1] - value * (spread[1] - spread[0]);
    });
  }
}

// generate a periodic function that goes from lower to upper at the given rate
function periodic(lower, upper, steps) {
  return function(index) {
    var step = index % steps;
    var spread = (upper - lower);
    return step * spread / steps + lower;
  }
}

// create a periodic function for each color, both bottom and top bounds, given rates for
// each of low, mid, and high
function buildBounds(low, mid, high) {
  var redBottom = periodic(0, 0.4, low);
  var greenBottom = periodic(0, 0.4, mid);
  var blueBottom = periodic(0, 0.4, high);

  var redTop = periodic(0.6, 1.0, mid);
  var greenTop = periodic(0.6, 1.0, high);
  var blueTop = periodic(0.6, 1.0, low);

  return {
    bottom: [redBottom, greenBottom, blueBottom],
    top: [redTop, greenTop, blueTop]
  }
}

// given a set of scales and an index, create the periodic functions for both bottom and top
// bounds based on the given index
function colorBounds(scales, index) {
  return {
    bottom: _.map(scales.bottom, function(scale) {return scale(index)}),
    top: _.map(scales.top, function(scale) {return scale(index)}),
  }
}

// given a set of keys and three rates (low, mid and high), generate a function that creates
// the bottom and top bounds for each key
function buildColorScales(keys, low, mid, high) {
  var generators = buildBounds(low, mid, high);
  var scales = {}
  _.each(keys, function(key, index) {
    var bounds = colorBounds(generators, index);
    scales[key] = colorScale(bounds.bottom, bounds.top);
  });

  return scales;
}
