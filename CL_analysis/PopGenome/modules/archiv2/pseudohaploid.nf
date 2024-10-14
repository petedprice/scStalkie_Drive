process pseudohaploid {

    errorStrategy 'retry'

    maxRetries 8
    cpus = { 4 * task.attempt }
    memory = { 32.GB * task.attempt }
    time = {4.hour } //* task.attempt }
    label 'R'

    tag {'pseudohaploidisation_' + contig }

    input:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_filt.vcf.gz")

    output:
    tuple val(vcf), val(contig), val(ctg_len), file("${contig}_filt.vcf.gz")    
    
    script:
    """
    #!/bin/bash

    echo 'dog again4'

    Rscript ${projectDir}/scripts/shuffle_pseudohaploid.R ${contig}_filt.vcf.gz contig ${projectDir}/sample_info.csv ${contig}_filt_PH.vcf.gz
    mv ${contig}_filt_PH.vcf.gz ${contig}_filt.vcf.gz


    """
}
