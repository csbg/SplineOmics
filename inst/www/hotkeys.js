// hotkeys.js
document.addEventListener('DOMContentLoaded', (event) => {
    document.addEventListener('keydown', function(event) {
        if (event.key === 't') {
            document.getElementById('toc').scrollIntoView();
        }
        if (event.key === 'e') {
            var email = 'your.email@example.com';
            var mailto_url = 'mailto:' + email;
            var new_tab = window.open(mailto_url, '_blank');
            if (new_tab) {
                new_tab.opener = null;  // Prevents a possible security vulnerability
            }
        }
        if (event.key === 'd') {
            downloadEmbeddedFiles();
        }
    });
});

function downloadEmbeddedFiles() {
    var zip = new JSZip();
    var files = [
        // Add your embedded files here
        {name: 'file1.txt', content: 'Content of file 1'},
        {name: 'file2.txt', content: 'Content of file 2'}
    ];

    files.forEach(function(file) {
        zip.file(file.name, file.content);
    });

    zip.generateAsync({type: 'blob'}).then(function(content) {
        saveAs(content, 'embedded_files.zip');
    });
}

