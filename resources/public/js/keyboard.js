function initializeKeyboard() {
  var keyboard = {};

  document.addEventListener('keydown', function(event) {
    if (!keyboard[event.code]) {
      console.log('key down! ' + event.code);
      keyboard[event.code] = true;
    }
  });

  document.addEventListener('keyup', function(event) {
    console.log('key up..... ' + event.code);
    keyboard[event.code] = null;
  });

  return keyboard;
}
