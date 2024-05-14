process split_bam {

    label 'samtoolsetc'
    errorStrategy 'retry'

    cpus { 1 * task.attempt }
    errorStrategy 'retry'
    maxRetries 16
    memory { 8.GB * task.attempt }

    input:
    tuple val(species), val(sample),  file("${sample}_crdata"), file("contig.txt")

    output:
    tuple val(species), val(sample), env(contig), file("subset.bam"), file("subset.bam.bai")

    script:
    """
    #!/bin/bash

    echo dog
    contig=\$(cat contig.txt)
    echo \$contig
    samtools view -bh ${sample}_crdata/outs/possorted_genome_bam.bam \$contig > subset.bam
    samtools index subset.bam
    """
}

