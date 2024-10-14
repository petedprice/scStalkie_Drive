process PopGenome {
    errorStrategy 'retry'

    maxRetries 8
    cpus = { 4 * task.attempt }
    memory = { 32.GB * task.attempt }
    time = {4.hour } //* task.attempt }
    label 'R'

    tag {'PopGenome_' + contig + '_' + bedtype}

    publishDir 'PopGenome_Stats', mode: 'copy', overwrite: true, pattern: '*_PopGenome_Stats.csv'

    input:
    tuple val(contig), val(ctg_len), val(bedtype), file(vcf), file(index)

    output:
    tuple val(contig), val(bedtype), file("${bedtype}_${contig}_PopGenome_Stats.csv")

    script:
    """
    #!/bin/bash

    echo 'rerunning again but i think this is fine now'

    Rscript ${projectDir}/scripts/PopGenome.R $vcf $contig $ctg_len ${params.sample_info} $bedtype ${projectDir}/R_packages ${params.gff}
    
    """
}
