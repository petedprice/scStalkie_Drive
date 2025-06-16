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
    cat ${params.gff_dir}/${ref}.gff | grep gene | cut -f 1 | uniq > contigs.txt
    #cat ${params.fasta_dir}/${ref}.fna | grep '>' | cut -d ' ' -f1 | cut -c2- > contigs.txt
    for contig in \$(cat contigs.txt)
    do
    echo \$contig > \${contig}_${species}.txt
    done
    rm contigs.txt
    """
}
