process orthofinder {
    cpus = 8
    memory = '32 GB'
    time = '4h'

    label 'python_orth'
    //conda 'anaconda::scipy anaconda::numpy'

    //publishDir 'orthofinder', mode: 'copy', overwrite: true, pattern: 'Orthofinder_Results'

    input:
    file(proteins)
    output:
    tuple file(proteins), file('Orthofinder_Results')    

    script:
    """
    #!/bin/bash
    echo 'running orthofinder'
    mkdir proteins
    mv *faa proteins
    python ${baseDir}/software/OrthoFinder_source/orthofinder.py -f proteins
    mv proteins/OrthoFinder/Results_* Orthofinder_Results
    mv proteins/*faa .
    """
}
