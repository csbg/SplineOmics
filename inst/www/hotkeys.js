document.addEventListener('DOMContentLoaded', (event) => {
    document.addEventListener('keydown', function(event) {
        console.log('Key pressed:', event.key);
        if (event.key === 't') {
            document.getElementById('toc').scrollIntoView();
        }
        if (event.key === 'e') {
            var email = '{{email}}';
            var fullName = '{{name}}';
            var datetime = '{{datetime}}';
            var subject = `SplineOmics HTML report Date-Time: ${datetime}`;
            var body = `Dear ${fullName},\n\n`;

            var mailto_url = `mailto:${email}?subject=${encodeURIComponent(subject)}&body=${encodeURIComponent(body)}`;
            var new_tab = window.open(mailto_url, '_blank');
            if (new_tab) {
                new_tab.opener = null;  // Prevents a possible security vulnerability
            }
        }
        if (event.key === 'd') {
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

