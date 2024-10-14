process PH_PopGenome {
//    errorStrategy { 'retry' : 'terminate' }
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    maxRetries 8
    cpus = { 1 }
    memory = { 128.GB * task.attempt }
    time = {4.hour } //* task.attempt }
    label 'R'

    tag {'PopGenome_' + contig + '_' + bedtype}

    publishDir 'PopGenome_Stats', mode: 'copy', overwrite: true, pattern: '*_PopGenome_Stats.csv'

    input:
    tuple val(contig), val(ctg_len), val(bedtype), file(vcf), file(index)

    output:
    tuple val(contig), val(bedtype), file("subsettedAll_PH_${bedtype}_${contig}_PopGenome_Stats.csv"), file(vcf)

    script:
    """
    #!/bin/bash

    echo 'after shuffling the genotypes but keeping all samples for autosomes and Z but correct numbers but only 8As and 8Bs for Z'

    Rscript ${projectDir}/scripts/PH_PopGenome.R $vcf $contig $ctg_len ${params.sample_info} $bedtype ${projectDir}/R_packages ${params.gff}

    
    """
}
