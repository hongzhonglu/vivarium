function testRotation() {
  var draw = Snap('#test');
  var box = draw.rect(0, 0, 100, 161.5);
  var rotation = 268;
  var increment = 1;
  var transform = 't200,200r'+rotation+',50,80.75';
  console.log(transform)
  box.attr({
    fill: '#000000',
    transform: transform});

  function rotate() {
    if (rotation > 271) increment = -1;
    if (rotation < 269) increment = 1;
    rotation += increment;
    var transform = 't200,200r'+rotation+',50,80.75';
    console.log(rotation);
    box.animate({transform: transform}, 1000);
    setTimeout(rotate, 1000);
  }

  setTimeout(rotate, 1000)
}
