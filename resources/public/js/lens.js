// establish the various constants used in the visualization
TAU = 2 * Math.PI;
INVERSE_TAU = 1.0 / TAU;
VISUALIZATION_WIDTH = 1000;
PATCHES_PER_EDGE = 10;
DISPLAY_PRECISION = 4;
GLUCOSE_KEY = 'GLC[p]'

// generate a uuid
function uuid() {
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    var r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
    return v.toString(16);
  });
}

// given the volume of a cell and its radius, calculate the length of the cell.
function volumeToLength(radius, volume) {
  var ratio = (4/3) * Math.PI * Math.pow(radius, 3);
  var area = Math.PI * Math.pow(radius, 2);
  var cylinder = (volume - ratio) / area;
  return cylinder + 2 * radius;
}

// find the closest rotation between two points on a circle (0..TAU)
function relativeRotation(before, after) {
  var above = after - before;
  var below = after < before ? (after + TAU) - before : after - (before + TAU);
  return Math.abs(above) < Math.abs(below) ? above : below;
}

// translate a given theta into degrees
function rotationTransform(theta) {
  return (TAU - theta) * INVERSE_TAU * 360;
}

// given an object containing the various svg elements of a cell visualization, apply the new
// data to each of its components
function updateCell(cell, data, born) {
  // calculate length from volume
  var length = volumeToLength(data.width, data.volume);

  // calculate the offset from the corner of the rectangle to the center, given
  // the cell's orientation.
  var originX = data.width * 0.5
  var originY = length * 0.5;
  var sin = Math.sin(-data.orientation);
  var cos = Math.cos(-data.orientation);
  var offsetX = originX * cos - originY * sin;
  var offsetY = originX * sin + originY * cos;

  // apply the offset
  var cx = data.location[1] - originX // - offsetX;
  var cy = data.location[0] - originY // - offsetY;

  var previousRotation = /r([^,]+)/.exec(cell.whole.transform()['string']);
  var orientation = data.orientation;
  if (previousRotation) {
    var previous = parseFloat(previousRotation[1]) / 360.0;
    var delta = relativeRotation(previous, orientation);
    orientation = previous + delta;
  }

  console.log(data.location)
  console.log(length)

  var translation = [cx * data.scale, cy * data.scale];
  var transform = 't' + translation[0]
      + ',' + translation[1]
      + 'r' + rotationTransform(orientation)
      + ',' + originX * data.scale
      + ',' + originY * data.scale;
  console.log(transform)

  // transform the svg group as a whole
  var sizeUpdate = {
    width: data.width * data.scale,
    height: length * data.scale,
  }

  born = true;
  if (born) {
    cell.whole.attr({transform: transform})
    cell.membrane.attr(sizeUpdate);
  } else {
    cell.whole.animate({transform: transform}, 1000)
    cell.membrane.animate(sizeUpdate, 1000);
  }

  // // translate the center point to the center of the membrane
  // var nucleoid = born ? cell.nucleoid : cell.nucleoid.animate()
  // nucleoid
  //   .cx(originX * PATCH_WIDTH)
  //   .cy(originY * PATCH_WIDTH)
}

// build a cell from data, initializing the svg group and passing the result to updateCell
// to apply the rest of the transformations
function buildCell(lens, draw, id, data) {
  // make the container
  var whole = draw
      .group()
      .click(function() {
        console.log('membrane click! ' + id);
        event = {agent_id: id};

        if (lens.keyboard.Backspace) {
          event.event = 'REMOVE_AGENT';
        } else if (lens.keyboard.Space) {
          event.event = 'DIVIDE_CELL';
        }

        if (event.event) {
          lens.send(event)
        }
      })

  // create the rectangle representing the outer bounds of the cell
  var membrane = whole.rect(0, 0, data.width * data.scale, data.scale)
  membrane.attr({
    rx: 0.3 * data.scale,
    ry: 0.3 * data.scale,
    fill: rgbToHex([0.1, 0.1, 0.1])})
      
  // // create the center point
  // var nucleoid = whole
  //     .circle(30)
  //     .fill(rgbToHex([1, 1, 1]));

  // create an object containing a reference to each component of the svg group
  var cell = {
    whole: whole,
    membrane: membrane
    // nucleoid: nucleoid
  }

  // apply the transformations implied by the supplied data.
  updateCell(cell, data, true);
  return cell;
}

// build a bunch of random cells (for demonstration only)
function buildCells(draw, n, scale, radius) {
  var cells = {};

  for (c = 0; c < n; c++) {
    var x = Math.random() * scale;
    var y = Math.random() * scale;
    var volume = Math.random() + radius;
    var orientation = Math.random() * TAU;
    var data = {
      location: [x, y],
      orientation: orientation,
      scale: scale,
      width: radius,
      volume: volume};

    var cell = buildCell({}, draw, c, data);
    cells[c] = cell;
  }

  return cells;
}

// build the lattice representing the patches of molecular concentrations
function buildLattice(lens, draw, m, n, width) {
  var rows = [];
  for (i = 0; i < m; i++) {
    var column = [];
    for (j = 0; j < n; j++) {
      var rect = draw.rect(0, 0, Math.ceil(width), Math.ceil(width))
      rect.attr({
        x: i*width,
        y: j*width,
        fill: rgbToHex([i / m, 1 - (j / m), 1])
      })

      rect.click(function(patch) {
        return function(event) {
          var x = event.offsetY * lens.edgeLength / VISUALIZATION_WIDTH;
          var y = event.offsetX * lens.edgeLength / VISUALIZATION_WIDTH;
          // var x = event.offsetY * VISUALIZATION_WIDTH / lens.edgeLength;
          // var y = event.offsetX * VISUALIZATION_WIDTH / lens.edgeLength;

          console.log('lattice click! ' + [x, y]);

          lens.send({
            event: 'ADD_AGENT',
            agent_id: uuid(),
            agent_type: lens.data.agent_type,
            agent_config: {
              'outer_id': lens.data.outer_id,
              'location': [x, y]
            }
          });
        }
      } ([i, j]))

      column.push(rect);
    }

    rows.push(column);
  }

  return rows;
}

function buildConcentrationBar(draw, lens, molecule, bounds) {
  var concentration = draw.group();
  var colors = lens.colors[molecule];
  var gradient = draw.gradient('l(0, 1, 0, 0)'+rgbToHex(colors(0))+'-'+rgbToHex(colors(1)));

  var height = VISUALIZATION_WIDTH * 0.56;
  var width = 60;
  var top = 220;
  var left = VISUALIZATION_WIDTH + 20;
  var bottom = top + height;
  var right = left + width;

  var bar = concentration.rect(0, 0, width, height)
  bar.attr({
    x: left,
    y: top,
    fill: gradient});

  var min = concentration.text(left, bottom+20, bounds[0].toPrecision(DISPLAY_PRECISION))
  var max = concentration.text(left, top-10, bounds[1].toPrecision(DISPLAY_PRECISION))

  return {
    bar: bar,
    min: min,
    max: max,
    gradient: gradient
  };
}

function setStops(gradient, colors) {
  var stops = gradient.stops();
  _.each(_.zip(colors, stops), function(pair) {
    var color = pair[0];
    var stop = pair[1];
    stop.node.setAttribute('stop-color', color);
  })
}

function updateConcentrationBar(draw, lens, molecule, bounds) {
  var colors = lens.colors[molecule];
  var hexes = _.map(_.range(2), function(index) {
    return rgbToHex(colors(index));
  });

  setStops(lens.bar.gradient, hexes)
  lens.bar.min.attr({text: bounds[0].toPrecision(DISPLAY_PRECISION)});
  lens.bar.max.attr({text: bounds[1].toPrecision(DISPLAY_PRECISION)});
}

function removeLattice(lattice) {
  var width = lattice.length;
  var height = lattice[0].length;

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      lattice[i][j].remove();
    }
  }
}

function checkLattice(draw, lens, data) {
  var width = data.lattice[GLUCOSE_KEY].length;
  var height = data.lattice[GLUCOSE_KEY][0].length;
  if (lens.plate.length !== width ||
      lens.plate[0].length !== height) {
    removeLattice(lens.plate);
    var patchWidth = VISUALIZATION_WIDTH / Math.max(width, height);
    lens.plate = buildLattice(lens, draw, width, height, patchWidth);
    _.each(_.values(lens.cells), function(cell) {
      lens.draw.add(cell.whole);
    });
  }
}

// update state of lattice with chosen molecule concentrations
function updateLattice(draw, lens, data) {
  if (data && data.lattice) {
    checkLattice(draw, lens, data);

    // get concentrations and colors for current molecule
    var rows = data.lattice[lens.molecule];
    var colors = lens.colors[lens.molecule];

    // calculate the minimum and maximum levels of the concentrations in every patch
    var max = -Infinity;
    var min = Infinity;
    for (var i = 0; i < rows.length; i++) {
      max = Math.max(max, Math.max.apply(null, rows[i]));
      min = Math.min(min, Math.min.apply(null, rows[i]));
    }

    var spread = max - min;
    if (isNaN(spread)) {
      spread = 0
    }

    // apply a color to each patch given its level relative to the minimum and maximum detected
    for (var i = 0; i < rows.length; i++) {
      var column = rows[i];
      for (var j = 0; j < column.length; j++) {
        var value = spread === 0 ? 1 : (rows[i][j] - min) / spread;
        lens.plate[j][i].animate({
          fill: rgbToHex(colors(value))}, 1000);
      }
    }

    var bounds = [min, max];
    updateConcentrationBar(draw, lens, lens.molecule, bounds);
  }
}

// given the new state of the lens from the simulation, apply all of the
// updates to the various elements of the visualization
function updateLens(draw, lens, data) {
  lens.edgeLength = data.edge_length;

  // remove any cells that are no longer referenced in the data
  _.each(_.keys(lens.cells), function(key) {
    if (!_.has(data.simulations, key)) {
      lens.cells[key].whole.remove();
      lens.cells = _.omit(lens.cells, [key]);
    }
  })

  // add any new cells that have appeared since the last set of data and update
  // existing cells with their new state
  var scale = VISUALIZATION_WIDTH / data.edge_length;
  _.each(_.keys(data.simulations), function(key) {
    var simulation = data.simulations[key];
    simulation.width = data.cell_radius;
    simulation.scale = scale

    if (!_.has(lens.cells, key)) {
      var cell = buildCell(lens, draw, key, simulation);
      lens.cells[key] = cell;
    } else {
      updateCell(lens.cells[key], simulation);
    }
  })

  lens.trigger.innerHTML = data.running ? 'pause' : 'trigger';

  // keep a reference to the new simulation states
  lens.data = data;

  updateLattice(draw, lens, data);
}

function buildMoleculeSelect(select, molecules) {
  _.each(molecules, function(molecule) {
    option = new Option(molecule, molecule);
    option.selected = molecule === GLUCOSE_KEY;
    select.append(option);
  });
}

// initialize the visualization
function bootLens(lens) {
  var lens = {edgeLength: 10.0};

  // set up SVG
  var draw = Snap('#lens');
  var plate = buildLattice(
    lens,
    draw,
    PATCHES_PER_EDGE,
    PATCHES_PER_EDGE,
    VISUALIZATION_WIDTH / PATCHES_PER_EDGE);

  // initialize websockets
  var socket = new WebSocket('ws://localhost:33332/ws');

  // initialize keyboard tracking
  var keyboard = initializeKeyboard();

  // initialize_buttons
  var trigger = document.getElementById('trigger');
  trigger.addEventListener('click', function(event) {
    if (lens.data && lens.data.running) {
      socket.send(JSON.stringify({'event': 'PAUSE_ALL'}));
      lens.trigger.innerHTML = 'trigger';
    } else {
      socket.send(JSON.stringify({'event': 'TRIGGER_ALL'}));
      lens.trigger.innerHTML = 'pause';
    }
  });

  // assign glucose color scale directly
  glucoseGradient = colorScale(
    [0, 0.3, 0],
    [0.9, 1, 0.5]);

  // set up molecule select
  var moleculeSelect = document.getElementById('molecule');
  var select = buildMoleculeSelect(moleculeSelect, moleculeNames);
  var colors = buildColorScales(moleculeNames, 7, 5, 4);
  colors[GLUCOSE_KEY] = glucoseGradient;

  // set up lens state
  lens.draw = draw;
  lens.plate = plate;
  lens.cells = {};
  lens.centers = {};
  lens.socket = socket;
  lens.send = function(message) {
    console.log(message);
    socket.send(JSON.stringify(message))
  }
  lens.initialized = false;
  lens.colors = colors;
  lens.trigger = trigger;
  lens.molecule = moleculeSelect.value;
  lens.keyboard = keyboard;
  lens.bar = buildConcentrationBar(draw, lens, GLUCOSE_KEY, [0, 0]);

  // add event listeners
  socket.addEventListener('open', function(event) {
    socket.send(JSON.stringify({'event': 'VISUALIZATION_INITIALIZE'}));
  });

  socket.addEventListener('message', function(event) {
    var data = JSON.parse(event.data)
    console.log(data);

    if (data['environment-state']) {
      updateLens(draw, lens, data['environment-state']);
    }
  });

  moleculeSelect.addEventListener('change', function() {
    lens.molecule = moleculeSelect.value;
    updateLattice(draw, lens, lens.data);
  });

  return lens;
}
