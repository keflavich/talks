
    // Function to check if CSS loaded and add fallback if not
    function addFallbackCSS(localId, fallbackUrl, fallbackId) {
    var localLink = document.getElementById(localId);
    if (localLink) {
        // Create a test element to check if styles are applied
        var testEl = document.createElement('div');
        testEl.style.cssText = 'position:absolute;top:-9999px;left:-9999px;';
        document.body.appendChild(testEl);

        // Use setTimeout to allow CSS to load
        setTimeout(function() {
        var isLoaded = false;
        try {
            // Check if the stylesheet loaded by testing for rules
            var sheets = document.styleSheets;
            for (var i = 0; i < sheets.length; i++) {
            if (sheets[i].href && sheets[i].href.indexOf(localLink.href.split('/').pop()) !== -1) {
                try {
                // Try to access cssRules to confirm the sheet loaded
                var rules = sheets[i].cssRules || sheets[i].rules;
                if (rules && rules.length > 0) {
                    isLoaded = true;
                    break;
                }
                } catch(e) {
                // Cross-origin or other error - assume it didn't load
                }
            }
            }
        } catch(e) {}

        document.body.removeChild(testEl);

        if (!isLoaded) {
            // Add fallback CSS
            var fallbackLink = document.createElement('link');
            fallbackLink.rel = 'stylesheet';
            fallbackLink.href = fallbackUrl;
            fallbackLink.id = fallbackId;
            document.head.appendChild(fallbackLink);
        }
        }, 100);
    }
    }

    // Add fallbacks for RevealJS CSS files
    document.addEventListener('DOMContentLoaded', function() {
    addFallbackCSS('reveal-css', 'https://cdn.jsdelivr.net/npm/reveal.js@5.2.1/dist/reveal.css', 'reveal-css-fallback');
    addFallbackCSS('zenburn-css', 'https://cdn.jsdelivr.net/npm/reveal.js@5.2.1/plugin/highlight/zenburn.css', 'zenburn-css-fallback');
    });

    function loadFallbackScript(fallbackUrl, fallbackId) {
    if (!document.getElementById(fallbackId)) {
        var fallbackScript = document.createElement('script');
        fallbackScript.src = fallbackUrl;
        fallbackScript.id = fallbackId;
        document.head.appendChild(fallbackScript);
    }
    }