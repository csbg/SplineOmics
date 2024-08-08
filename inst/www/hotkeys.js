// Array to store positions where 't' key is pressed
let savedPosition = null;

document.addEventListener('DOMContentLoaded', (event) => {
    console.log('DOMContentLoaded event triggered');

    document.addEventListener('keydown', function(event) {
        console.log('Key pressed:', event.key);
	
        if (event.key === 't') {
            // Store the current scroll position
		    savedPosition = window.scrollY;
        	console.log('Position saved:', savedPosition);

            // Ensuring correct element scrolls into view
            const tocElement = document.getElementById('toc');
            if (tocElement) {
                tocElement.scrollIntoView();
            } else {
                console.log('TOC element not found');
            }
        } else if (event.key === 's') {
            // Store the current scroll position
		    savedPosition = window.scrollY;
        	console.log('Position saved:', savedPosition);
        } else if (event.key === 'b') {
            // Jump back to the saved scroll position
            if (savedPosition !== null) {
                window.scrollTo(0, savedPosition);
                console.log('Jumped back to position:', savedPosition);
            } else {
                console.log('No position saved.');
            }
        } else if (event.key === 'e') {
            var email = '{{email}}';
            var fullName = '{{name}}';
            var subject = `SplineOmics HTML report`;
            var body = `Dear ${fullName},`;

            var mailto_url = `mailto:${email}?subject=${encodeURIComponent(subject)}&body=${encodeURIComponent(body)}`;
            var new_tab = window.open(mailto_url, '_blank');
            if (new_tab) {
                new_tab.opener = null;  // Prevents a possible security vulnerability
            }
        } else if (event.key === 'd') {
            console.log('Download embedded files');
            downloadEmbeddedFiles();
        }
    });
});

function downloadEmbeddedFiles() {
    var zip = new JSZip();
    var files = document.querySelectorAll('.embedded-file');
    console.log('Found embedded files:', files.length);

    if (files.length === 0) {
        console.error('No embedded files found.');
        return;
    }

    var promises = [];

    files.forEach(function(file, index) {
        var link = file.getAttribute('href');
        var filename = file.getAttribute('download');
        console.log('Processing file:', filename, link);

        var promise = fetch(link)
            .then(res => {
                if (!res.ok) {
                    throw new Error('Network response was not ok for ' + filename);
                }
                return res.blob();
            })
            .then(blob => {
                zip.file(filename, blob);
            })
            .catch(error => {
                console.error('Error fetching file:', error);
            });

        promises.push(promise);
    });

    Promise.all(promises)
        .then(() => {
            console.log('Generating zip file');
            return zip.generateAsync({type: 'blob'});
        })
        .then(content => {
            console.log('Saving zip file');
            saveAs(content, 'embedded_files.zip');
        })
        .catch(error => {
            console.error('Error generating zip file:', error);
        });
}

