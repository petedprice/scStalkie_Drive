process contig_names {
    cpus = 1
    memory = '4 GB'
    time = '4h'


    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("*.txt")

    script:
    """
    #!/bin/bash

    mkdir contigs
    cat ${params.gtf_dir}/${ref}.gtf | cut -f 1 | uniq > contigs.txt
    for contig in \$(cat contigs.txt)
    do
    echo \$contig > \${contig}_${species}.txt
    done
    rm contigs.txt
    """
}
