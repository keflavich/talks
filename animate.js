// Animate items
Reveal.addEventListener('slidechanged', function( event ) {
    if(event.currentSlide.classList.contains('to-animate')) {
        var vids = [event.currentSlide];
    
        for (var i = 0; i < vids.length; i++) {
            var vid = vids[i];
            vid.classList.add('animate');
            vid.classList.remove('to-animate');
        }
    }

    // reset slides when leaving
    if(event.previousSlide.classList.contains('animate')) {
        var vids = [event.previousSlide];
    
        for (var i = 0; i < vids.length; i++) {
            var vid = vids[i];
            vid.classList.add('to-animate');
            vid.classList.remove('animate');
        }
    }

});
