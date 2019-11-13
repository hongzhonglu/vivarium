// establish the various constants used in the visualization
TAU = 2 * Math.PI;
INVERSE_TAU = 1.0 / TAU;
VISUALIZATION_WIDTH = 1000;
PATCHES_PER_EDGE = 10;
DISPLAY_PRECISION = 8;
DEFAULT_COLOR = [0.6, 0.4, 0.3]
RGB_SIZE = 255
MEMBRANE_COLOR = [102/RGB_SIZE, 102/RGB_SIZE, 255/RGB_SIZE]
COMPARTMENT_OFFSET = {
  periplasm : 0,
  cytoplasm: 0.09
}

// generate a uuid
function uuid() {
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
    var r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
    return v.toString(16);
  });
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

  var length = data.length;
  var width = data.width;
  var cx_corner = data.corner_location[1];
  var cy_corner = data.corner_location[0];

  var previousRotation = /r([^,]+)/.exec(cell.whole.transform()['string']);
  var orientation = data.orientation;
  if (previousRotation) {
    var previous = parseFloat(previousRotation[1]) / 360.0;
    var delta = relativeRotation(previous, orientation);
    orientation = previous + delta;
  }

  console.log(data.location)
  console.log(length)

  // create the matrix representing the successive application of rotation and translation
  var transform = new SVG.Matrix()
      .translate(cx_corner * data.scale, cy_corner * data.scale)
      .rotate(rotationTransform(data.orientation))

  function animateExisting(obj) {
    return born ? obj : obj.animate();
  }

  function animateCapsule(group, data, offset, color) {
	    animateExisting(group)
	    .attr({
	      width: data.scale * width * (1 - offset),
	      height: data.scale * length * (1 - 0.5*offset),  // why does 0.5*offset work here?
	      x: data.scale * 0.5*offset,
	      y: data.scale * 0.5*offset
	    })
	    .fill(rgbToHex(color));
	  }

  // transform the svg group as a whole
  animateExisting(cell.whole)
    .transform(transform)
    .fill(rgbToHex(data.color || DEFAULT_COLOR));

  // increase the length of the outerMembrane
  animateCapsule(cell.periplasm, data, COMPARTMENT_OFFSET.periplasm, data.color)
  animateCapsule(cell.cytoplasm, data, COMPARTMENT_OFFSET.cytoplasm, data.color)

  // heads up display
  var hudX = data.location[1];  // - (0.4 * data.width);
  var hudY = data.location[0];  // - (0.4 * data.length);
  var translate = new SVG.Matrix()
      .translate(hudX * data.scale, hudY * data.scale);

  cell.hud.text(data.volume.toPrecision(DISPLAY_PRECISION));
  animateExisting(cell.hud)
    .transform(translate);
}

// build a cell from data, initializing the svg group and passing the result to updateCell
// to apply the rest of the transformations
function buildCell(lens, draw, id, data) {
  // make the container
  var container = draw.group();

  var whole = container
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

  // create the rectangle representing the compartments of the cell
  var periplasm = capsuleShape(whole, COMPARTMENT_OFFSET.periplasm, data, data.color)
  var cytoplasm = capsuleShape(whole, COMPARTMENT_OFFSET.cytoplasm, data, data.color)

  var hud = container
      .text("hello")
      .fill(rgbToHex([0.9, 0.9, 0.9]))
      .attr({'fill-opacity': 0.0});

  // // create the center point
  // var nucleoid = whole
  //     .circle(30)
  //     .fill(rgbToHex([1, 1, 1]));

  // create an object containing a reference to each component of the svg group
  var cell = {
    container: container,
    whole: whole,
    periplasm: periplasm,
    cytoplasm: cytoplasm,
    hud: hud,
    hovering: false
    // nucleoid: nucleoid
  }

  container
    .mouseover(function() {
      console.log('mouseover');
      if (!cell.hovering) {
        cell.hovering = true;
        cell.hud.attr({'fill-opacity': 1.0});
        // cell.hud.animate(200).attr({'fill-opacity': 1.0});
      }
    })
    .mouseout(function() {
      console.log('mouseout');
      if (cell.hovering) {
        cell.hovering = false;
        cell.hud.attr({'fill-opacity': 0.0});
        // cell.hud.animate(200).attr({'fill-opacity': 0.0});
      }
    })

  // apply the transformations implied by the supplied data.
  updateCell(cell, data, true);
  return cell;
}

function capsuleShape(group, offset, data, color) {
	  // create the rectangle representing the outer bounds of the capsule
	  return group
	      .rect(data.width * data.scale, data.scale)  // width, length
	      // .x(data.scale * offset)
	      // .y(data.scale * offset)
	      .rx(0.6 * data.scale)
	      .ry(0.6 * data.scale)
	      .attr({
	          fill: rgbToHex(color || DEFAULT_COLOR),
	          stroke: rgbToHex(MEMBRANE_COLOR),
	          'stroke-width': 2,
	          'stroke-opacity': 0.9,
	        })
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
      var rect = draw
          .rect(width, width)
          .x(i*width)
          .y(j*width)
          .fill(rgbToHex([i / m, 1 - (j / m), 1]))
          .width(width*1.05)  // Helps reduce grid-like artifact for certain lattice size
          .height(width*1.05)
          .click(function(patch) {
            return function(event) {
              var x = event.offsetY * lens.edgeLength / VISUALIZATION_WIDTH;
              var y = event.offsetX * lens.edgeLength / VISUALIZATION_WIDTH;
              console.log(event);
              console.log('lattice click! ' + [x, y]);
              lens.socket.send(JSON.stringify({
                event: 'ADD_AGENT',
                agent_id: uuid(),
                agent_type: lens.data.agent_type,
                agent_config: {
                  'outer_id': lens.data.outer_id,
                  'seed': Math.round(37 * x + 3 * y),
                  'location': [x, y]
                }
              }));
            }
          } ([i, j]))
          .back();

      column.push(rect);
    }

    rows.push(column);
  }

  return rows;
}

function buildConcentrationBar(draw, lens, molecule, bounds) {
  var concentration = draw.group();
  var colors = lens.colors[molecule];
  var gradient = draw.gradient('linear', function(stop) {
    stop.at(0, rgbToHex(colors(0)));
    stop.at(1, rgbToHex(colors(1)));
  }).attr({
    'x1': 0,
    'x2': 0,
    'y1': 1,
    'y2': 0
  });

  var height = VISUALIZATION_WIDTH * 0.56;
  var width = 60;
  var top = 220;
  var left = VISUALIZATION_WIDTH + 20;
  var bottom = top + height;
  var right = left + width;

  var bar = concentration
      .rect(width, height)
      .x(left)
      .y(top)
      .fill(gradient);

  var min = concentration
      .text(bounds[0].toPrecision(DISPLAY_PRECISION))
      .x(left)
      .y(bottom + 5);

  var max = concentration
      .text(bounds[1].toPrecision(DISPLAY_PRECISION))
      .x(left)
      .y(top-20);

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

  lens.bar.gradient.update(function(stop) {
    stop.at(0, rgbToHex(colors(0)));
    stop.at(1, rgbToHex(colors(1)));
  });

  lens.bar.min.text(bounds[0].toPrecision(DISPLAY_PRECISION));
  lens.bar.max.text(bounds[1].toPrecision(DISPLAY_PRECISION));
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
  var width = data.lattice[lens.moleculeNames[0]].length;
  var height = data.lattice[lens.moleculeNames[0]][0].length;
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
        lens.plate[j][i].animate()
          .attr({fill: rgbToHex(colors(value))});
        // lens.plate[j][i].animate({
        //   fill: rgbToHex(colors(value))}, 1000);
      }
    }

    var bounds = [min, max];
    updateConcentrationBar(draw, lens, lens.molecule, bounds);
  }
}

// given the new state of the lens from the simulation, apply all of the
// updates to the various elements of the visualization
function updateLens(draw, lens, data) {
  lens.edgeLength = data.edge_length_x;  // TODO use data.edge_length_x, data.edge_length_y

  // remove any cells that are no longer referenced in the data
  _.each(_.keys(lens.cells), function(key) {
    if (!_.has(data.simulations, key)) {
      lens.cells[key].container.remove();
      lens.cells = _.omit(lens.cells, [key]);
    }
  })

  // get moleculeNames from environment-state, set in lens
  var moleculeNames = _.keys(data.lattice);
  if (!_.isEqual(moleculeNames, lens.moleculeNames)) {
    lens.moleculeNames = moleculeNames;
    select = setMoleculeSelect(moleculeNames);
    lens.select = select.select;
    lens.colors = select.colors;
    lens.molecule = moleculeNames[0];
  }

  // add any new cells that have appeared since the last set of data and update
  // existing cells with their new state
  var scale = VISUALIZATION_WIDTH / data.edge_length_x;  // TODO use data.edge_length_x, data.edge_length_y
  _.each(_.keys(data.simulations), function(key) {
    var simulation = data.simulations[key];
    simulation.scale = scale

    if (!_.has(lens.cells, key)) {
      var cell = buildCell(lens, draw, key, simulation);
      lens.cells[key] = cell;
    } else {
      updateCell(lens.cells[key], simulation);
    }
  })

  lens.trigger.innerHTML = data.running ? 'Pause' : 'Run';

  // keep a reference to the new simulation states
  lens.data = data;

  updateLattice(draw, lens, data);
}

// set up molecule select
function setMoleculeSelect(moleculeNames) {
  // TODO (Eran) -- get keys from environment_state to use in place of 'molecules
  // assign glucose color scale directly
  glucoseGradient = colorScale(
    [0, 0.3, 0],
    [0.9, 1, 0.5]);

  var moleculeSelect = document.getElementById('molecule');
  var select = buildMoleculeSelect(moleculeSelect, moleculeNames);
  var colors = buildColorScales(moleculeNames, 7, 5, 4);
  colors[moleculeNames[0]] = glucoseGradient;
  
  return {
    colors: colors,
    select: moleculeSelect,
  }
}

function buildMoleculeSelect(select, molecules) {
  var firstMolecule = molecules[0];
  _.each(select.options, function(option) {select.remove(0)});
  _.each(molecules, function(molecule) {
    option = new Option(molecule, molecule);
    option.selected = molecule === firstMolecule;
    select.append(option);
  });
}

// initialize the visualization
function bootLens(lens) {
  var lens = {edgeLength: 10.0};

  // set up SVG
  var draw = SVG('lens');
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
      lens.trigger.innerHTML = 'Run';
    } else {
      socket.send(JSON.stringify({'event': 'TRIGGER_ALL'}));
      lens.trigger.innerHTML = 'Pause';
    }
  });

  // make the molecule selection menu
  lens.moleculeNames = ['molecule'];
  var select = setMoleculeSelect(lens.moleculeNames);
  lens.select = select.select;
  lens.colors = select.colors;
  lens.molecule = lens.select.value;

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
  lens.trigger = trigger;
  lens.keyboard = keyboard;
  lens.bar = buildConcentrationBar(draw, lens, lens.moleculeNames[0], [0, 0]);

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

  lens.select.addEventListener('change', function() {
    lens.molecule = lens.select.value;
    updateLattice(draw, lens, lens.data);
  });

  return lens;
}
