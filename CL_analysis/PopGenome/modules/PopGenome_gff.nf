process PopGenome_gff {
    errorStrategy 'retry'

    maxRetries 8
    cpus = { 4 * task.attempt }
    memory = { 32.GB * task.attempt }
    time = {4.hour } //* task.attempt }
    label 'R'

    tag {'PopGenome_' + contig}

    publishDir 'PopGenome_Stats', mode: 'copy', overwrite: true, pattern: '*_PopGenome_Stats'

    input:
    tuple val(contig), file("${contig}_filt_tbx.vcf.gz"), file("${contig}_filt_tbx.vcf.gz.tbi"), val(contig_length)

    output:
    tuple val(contig), file("${contig}_PopGenome_Stats")

    script:
    """
    #!/bin/bash


    Rscript ${projectDir}/scripts/PopGenome_gff.R ${contig}_filt_tbx.vcf.gz $contig $contig_length ${params.sample_info} default ${projectDir}/R_packages ${params.gff} ${params.fasta}
    mkdir ${contig}_PopGenome_Stats
    mv *csv ${contig}_PopGenome_Stats

    """
}
